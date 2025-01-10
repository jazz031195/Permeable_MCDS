import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns

import dipy.reconst.dki as dki

from dipy.core.gradients import gradient_table
import nibabel as nib
import glob
from Watson import WMTI_Watson
from useful_functions import read_and_extract_parameters, read_binary_file, array_to_nifti, get_files_from_folder
from DKI import read_and_extract_parameters, read_and_extract_dwi, array_to_nifti_replicated, dki_from_file, get_total_icvf
from scipy.stats import linregress

# Constants
GIRO_CONSTANT = 267.51525e3  # rad/msT

def get_c2(file_path):
    """
    Extract the "C2" value from a file.
    """
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("C2"):
                return float(line.split()[-1])
    return None  # Return None if not found

def assemble_intra_extra(directory, scheme_file):
    files = get_files_from_folder(directory)

    files_txt = get_files_from_folder(directory, False)
    for file_extra in files:

        if "img" not in file_extra and "info" not in file_extra and "extra" in file_extra and "swc" not in file_extra and "new" not in file_extra:
            print(file_extra)
            # replace extra by intra in file
            file_intra = file_extra.replace("extra", "intra")

            if file_intra in files:

                file_info_voxel = file_extra.split("extra")[0] + "info.txt"

                if file_info_voxel in files_txt:

                    print(f"Assembling intra and extra files for : {file_extra} and {file_intra}") 

                    icvf = get_total_icvf(file_info_voxel)
                    c2 = get_c2(file_info_voxel)

                    if icvf > 1 or icvf < 0:
                        print(f"icvf value not valid for voxel : {file_info_voxel}")
                        continue

                    # read the extra file
                    extra = read_binary_file(file_extra)

                    min_extra = np.min(extra)
                    if min_extra < 0:
                 
                        #extra = np.array([e - min_extra for e in extra])
                        extra = np.array([e if e > 0 else 0 for e in extra], dtype=extra.dtype)
                        # rewrite the extra file with the new values
                        new_file_extra = file_extra.replace(".bfloat", "_new.bfloat")

                        with open(new_file_extra, "wb") as f:
                            extra.tofile(f)
                        
                    else:
                        new_file_extra = file_extra

                    # read the intra file
                    intra = read_binary_file(file_intra)
                    min_intra = np.min(intra)
                    if min_intra < 0:
            
                        #intra = np.array([i - min_intra for i in intra])
                        intra = np.array([e if e > 0 else 0 for e in intra], dtype=intra.dtype)
                        # rewrite the intra file with the new values
                        new_file_intra = file_intra.replace(".bfloat", "_new.bfloat")
                        with open(new_file_intra, "wb") as f:
                            intra.tofile(f)
                    else:
                        new_file_intra = file_intra

                    sum =[] 
                    for e in range(len(extra)):
                        if (extra[e]*(1-icvf)+intra[e]*icvf) <0 :
                            print(f"negative value ! {extra[e]*(1-icvf)} + {intra[e]*icvf}")
                            assert False
                        sum.append(extra[e]*(1-icvf)+intra[e]*icvf)
                    sum = np.array(sum)

                    filename = file_extra.replace("_extra", "_extra_intra") 
                    filename = filename.replace(".bfloat", ".txt")
                    filename = filename.replace("_new", "")
                    with open(filename, "w") as f:
                        for item in sum:
                            f.write(f"{item}\n")

                    print(filename)
                    
                    _, MD_extra, AD_extra, RD_extra,_,_, _= dki_from_file(new_file_extra, scheme_file)
                    

                    _,MD_intra, AD_intra, RD_intra,_,_,_ = dki_from_file(new_file_intra, scheme_file)

                    filename_ = file_extra.split("_DWI")[0]
                    filename_ = filename_.replace("_extra", "")
                    # create SMI_ouput file 
                    with open(f"{filename_}_SMI_output.txt", "w") as file:
                        file.write(f" real f : {icvf}\n")
                        file.write(f" real c2 : {c2}\n")
                        file.write(f" real Da : {AD_intra}\n")
                        file.write(f" real Deperp : {RD_extra}\n")
                        file.write(f" real Depar : {AD_extra}\n")        
                else:
                    print(f"no info file for voxel : {file_info_voxel}")

            else:
                print(f"no intra file : {file_intra}")
                

def WMTI_from_file(file_path, scheme_file, binary):
    # Existing code to extract parameters and DWI data
    bs, Gs, deltas, Deltas, vectors, TEs = read_and_extract_parameters(scheme_file)
    folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff"

    b_values = bs

    # delete all files in /home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff/parameters_WMTI
    for file in glob.glob(f"{folder}/parameters_WMTI/*"):
        os.remove(file)

    with open(f"{folder}/bvals.txt", "w") as f:
        f.write("\n".join(map(str, b_values)) + "\n")

    np.savetxt(f"{folder}/vectors.txt", vectors)

    dwi = read_and_extract_dwi(file_path, binary)

    if (dwi.shape[-1] != len(b_values)):
        print(dwi)
        print("The number of gradients in the scheme file does not match the number of volumes in the DWI data.")
        assert False

    # Save as NIfTI
    img = array_to_nifti_replicated(dwi)
    nib.save(img, f"{folder}/dwi.nii.gz")

    # Convert NIfTI to MIF
    os.system(f"mrconvert {folder}/dwi.nii.gz {folder}/dwi.mif -force")
    os.system(f"mrconvert {folder}/dwi.mif {folder}/dwi.mif -fslgrad {folder}/vectors.txt {folder}/bvals.txt -force")

    # Create a fake sigma map (optional)
    sigma_map = create_fake_sigma_map(img.shape, sigma_value=0.001)

    save_sigma_map(sigma_map, folder +"/sigmas/")

    # Run the SMI calculation using Docker and the sigma map
    os.system(f"docker run -v /home/localadmin/Documents/MCDS/Permeable_MCDS/output:/data nyudiffusionmri/designer2:main tmi -DKI -echo_time {TEs[0]} -sigma /data/DKI_stuff/sigmas/sigma.nii.gz /data/DKI_stuff/dwi.mif /data/DKI_stuff/parameters_WMTI")

    for file in glob.glob(f"{folder}/parameters_WMTI/*"):
        name = file.split("/")[-1]
        if "fa" in name:
            FA = nib.load(file).get_fdata()
        elif "md" in name:
            MD = nib.load(file).get_fdata()
        elif "ad" in name:
            AD = nib.load(file).get_fdata()
        elif "rd" in name:
            RD = nib.load(file).get_fdata()
        elif "mk" in name:
            MK = nib.load(file).get_fdata()
        elif "ak" in name:
            AK = nib.load(file).get_fdata()
        elif "rk" in name:
            RK = nib.load(file).get_fdata()

    wmti = WMTI_Watson(f"{folder}/parameters_WMTI/", params='invivo', nodes=4)
    # Fit
    wmti.fit()
    f, Da, Depar, Deperp, c2 = wmti.maps()

    return f[0][0][0], Da[0][0][0], Depar[0][0][0], Deperp[0][0][0], c2[0][0][0]


def create_fake_sigma_map(dwi_shape, sigma_value=0.001):
    """
    Create a fake sigma map with a constant sigma value for all voxels.
    
    Parameters:
    -----------
    dwi_shape : tuple
        The shape of the DWI data (x, y, z, num_gradients).
    sigma_value : float, optional
        The constant sigma value to assign to each voxel. Default is 0.001.
    
    Returns:
    --------
    sigma_map : numpy array
        A 3D sigma map with the same shape as the DWI spatial dimensions.
    """
    print(dwi_shape)
    x_size, y_size, z_size, _ = dwi_shape
    sigma_map = np.full((x_size, y_size, z_size), sigma_value)
    return sigma_map

def save_sigma_map(sigma_map, folder, filename="sigma.nii.gz"):
    """
    Save the sigma map as a NIfTI file.
    
    Parameters:
    -----------
    sigma_map : numpy array
        The 3D sigma map.
    folder : str
        The folder where the sigma map should be saved.
    filename : str, optional
        The name of the sigma map file. Default is 'sigma.nii.gz'.
    """
    sigma_img = nib.Nifti1Image(sigma_map, affine=np.eye(4))
    nib.save(sigma_img, os.path.join(folder, filename))


def calculate_SMI(path_to_DWI, path_scheme, binary):
    # Existing code to extract parameters and DWI data
    bs, Gs, deltas, Deltas, vectors, TEs = read_and_extract_parameters(path_scheme)
    folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff"

    b_values = bs

    # delete all files in /home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff/parameters_DKI
    for file in glob.glob(f"{folder}/parameters_SMI/*"):
        os.remove(file)

    with open(f"{folder}/bvals.txt", "w") as f:
        f.write("\n".join(map(str, b_values)) + "\n")

    np.savetxt(f"{folder}/vectors.txt", vectors)

    dwi = read_and_extract_dwi(path_to_DWI, binary)


    if (dwi.shape[-1] != len(b_values)):
        print(dwi)
        print("The number of gradients in the scheme file does not match the number of volumes in the DWI data.")
        assert False

    # Save as NIfTI
    img = array_to_nifti_replicated(dwi)
    nib.save(img, f"{folder}/dwi.nii.gz")

    # Convert NIfTI to MIF
    os.system(f"mrconvert {folder}/dwi.nii.gz {folder}/dwi.mif -force")
    os.system(f"mrconvert {folder}/dwi.mif {folder}/dwi.mif -fslgrad {folder}/vectors.txt {folder}/bvals.txt -force")

    # Create a fake sigma map (optional)
    sigma_map = create_fake_sigma_map(img.shape, sigma_value=0.001)

    save_sigma_map(sigma_map, folder +"/sigmas/")

    # Run the SMI calculation using Docker and the sigma map
    os.system(f"docker run -v /home/localadmin/Documents/MCDS/Permeable_MCDS/output:/data nyudiffusionmri/designer2:main tmi -SMI -echo_time {TEs[0]} -sigma /data/DKI_stuff/sigmas/sigma.nii.gz /data/DKI_stuff/dwi.mif /data/DKI_stuff/parameters_SMI -compartments EAS,IAS")

    # Parse the results
    for file in glob.glob(f"{folder}/parameters_SMI/*"):
        name = file.split("/")[-1]
        if "f_smi" in name:
            f = nib.load(file).get_fdata()
        elif "Da" in name:
            Da = nib.load(file).get_fdata()
        elif "DePar" in name:
            Depar = nib.load(file).get_fdata()
        elif "DePerp" in name:
            Deperp = nib.load(file).get_fdata()
        elif "p2" in name:
            p2 = nib.load(file).get_fdata()


    return f[0][0][0], Da[0][0][0], Depar[0][0][0], Deperp[0][0][0], p2[0][0][0]


def process_SMI_files(directory):
    f_values = []
    real_f_values = []
    Deperps = []
    Depars = []
    Das = []
    p2s = []

    for filename in os.listdir(directory):
        if filename.endswith(".txt") and "SMI_output" in filename:
            print(filename)
            file_path = os.path.join(directory, filename)
            f_value, real_f_value, Da, Depar, Deperp, p2 = extract_SMI_output(file_path)
            if f_value is not None and real_f_value is not None:
                print(f"f: {f_value}, real f: {real_f_value}")
                f_values.append(f_value)
                real_f_values.append(real_f_value)
                Deperps.append(Deperp)
                Depars.append(Depar)
                Das.append(Da)
                p2s.append(p2)
    
    return f_values, real_f_values, Das, Depars, Deperps, p2s


def extract_SMI_output(file_path):
    SMI_f_value = None
    real_f_value = None
    SMI_Da = None
    SMI_Depar = None
    SMI_Deperp = None
    SMI_p2 = None
    WMTI_f_value = None
    WMTI_Da = None
    WMTI_Depar = None
    WMTI_Deperp = None
    WMTI_c2 = None
    real_Da = None
    real_Depar = None
    real_Deperp = None
    real_c2 = None
    

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("SMI f :"):
                SMI_f_value = float(line.split(":")[1].strip())
            elif line.startswith(" real f :"):
                real_f_value = float(line.split(":")[1].strip())
            elif line.startswith("SMI Da :"):
                SMI_Da = float(line.split(":")[1].strip())
            elif line.startswith(" real Da :"):
                real_Da = float(line.split(":")[1].strip())
            elif line.startswith("SMI Depar :"):
                SMI_Depar = float(line.split(":")[1].strip())
            elif line.startswith(" real Depar :"):
                real_Depar = float(line.split(":")[1].strip())
            elif line.startswith("SMI Deperp :"):
                SMI_Deperp = float(line.split(":")[1].strip())
            elif line.startswith(" real Deperp :"):
                real_Deperp = float(line.split(":")[1].strip())
            elif line.startswith("SMI p2 :"):
                SMI_p2 = float(line.split(":")[1].strip())
            elif line.startswith("WMTI f :"):
                WMTI_f_value = float(line.split(":")[1].strip())
            elif line.startswith("WMTI Da :"):
                WMTI_Da = float(line.split(":")[1].strip())
            elif line.startswith("WMTI Depar :"):
                WMTI_Depar = float(line.split(":")[1].strip())
            elif line.startswith("WMTI Deperp :"):
                WMTI_Deperp = float(line.split(":")[1].strip())
            elif line.startswith("WMTI c2 :"):
                WMTI_c2 = float(line.split(":")[1].strip())
            elif line.startswith(" real c2 :"):
                real_c2 = float(line.split(":")[1].strip())
    return SMI_f_value, real_f_value, SMI_Da, SMI_Depar, SMI_Deperp, SMI_p2, WMTI_f_value, WMTI_Da, WMTI_Depar, WMTI_Deperp, WMTI_c2, real_Da, real_Depar, real_Deperp, real_c2


def read_SMI_files(directory):
    SMI_f_values = []
    real_f_values = []
    SMI_Deperps = []
    SMI_Depars = []
    SMI_Das = []
    SMI_p2s = []
    real_Das = []
    real_Depars = []
    real_Deperps = []
    WMTI_f_values = []
    WMTI_Das = []
    WMTI_Depars = []
    WMTI_Deperps = []
    WMTI_c2s = []
    filenames = []
    real_c2s = []

    for filename in os.listdir(directory):
        if filename.endswith(".txt") and "SMI_output" in filename:
            print(filename)
            file_path = os.path.join(directory, filename)
            SMI_f_value, real_f_value, SMI_Da, SMI_Depar, SMI_Deperp, SMI_p2, WMTI_f_value, WMTI_Da, WMTI_Depar, WMTI_Deperp, WMTI_c2, real_Da, real_Depar, real_Deperp, real_c2 = extract_SMI_output(file_path)
            SMI_f_values.append(SMI_f_value)
            real_f_values.append(real_f_value)
            SMI_Deperps.append(SMI_Deperp)
            SMI_Depars.append(SMI_Depar)
            SMI_Das.append(SMI_Da)
            SMI_p2s.append(SMI_p2)
            real_Das.append(real_Da)
            real_Depars.append(real_Depar)
            real_Deperps.append(real_Deperp)
            WMTI_f_values.append(WMTI_f_value)
            WMTI_Das.append(WMTI_Da)
            WMTI_Depars.append(WMTI_Depar)
            WMTI_Deperps.append(WMTI_Deperp)
            WMTI_c2s.append(WMTI_c2)
            filenames.append(filename)
            real_c2s.append(real_c2)

    
    return filenames, SMI_f_values, real_f_values, SMI_Das, SMI_Depars, SMI_Deperps, SMI_p2s, WMTI_f_values, WMTI_Das, WMTI_Depars, WMTI_Deperps, WMTI_c2s, real_Das, real_Depars, real_Deperps, real_c2s


def plot_SMI_f_vs_real_f(directory, cyl = False):
    
    filenames, SMI_f_values, real_f_values, SMI_Das, SMI_Depars, SMI_Deperps, SMI_p2s, WMTI_f_values, WMTI_Das, WMTI_Depars, WMTI_Deperps, WMTI_c2s, real_Das, real_Depars, real_Deperps, real_c2s = read_SMI_files(directory)

    dataframe = pd.DataFrame()
    dataframe['SMI f'] = SMI_f_values
    dataframe['real f'] = real_f_values
    dataframe['SMI Da'] = SMI_Das
    dataframe['SMI Depar'] = SMI_Depars
    dataframe['SMI Deperp'] = SMI_Deperps
    dataframe['SMI p2'] = SMI_p2s
    dataframe['DKI Da'] = real_Das
    dataframe['DKI Depar'] = real_Depars
    dataframe['DKI Deperp'] = real_Deperps
    dataframe['WMTI f'] = WMTI_f_values
    dataframe['WMTI Da'] = WMTI_Das
    dataframe['WMTI Depar'] = WMTI_Depars
    dataframe['WMTI Deperp'] = WMTI_Deperps
    dataframe['WMTI c2'] = WMTI_c2s
    dataframe['filename'] = filenames
    dataframe['real c2'] = real_c2s
    if not cyl:
        dataframe['real p2'] = [(3*c-1)/2 for c in real_c2s]
    else:
        dataframe["real p2"] = np.ones(len(filenames))
        dataframe['real Da'] = np.ones(len(filenames))*2.5
        dataframe['real Depar'] = np.ones(len(filenames))*1.5



    # Create a figure with subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    axes = axes.flatten()  # Flatten for easy indexing

    # Handle the first plot separately
    sns.scatterplot(data=dataframe, x='real f', y="SMI f", label="SMI f", color='black', s=40, ax=axes[0])
    sns.lineplot(data=dataframe, x='real f', y="real f", label="GT f", color='grey', ax=axes[0])
    axes[0].set_title("Estimated f with SMI vs Real f")
    axes[0].set_xlabel("Real f")
    axes[0].set_ylabel("Estimated f")
    axes[0].legend()

    # 3. Plot "Estimated p2"
    sns.scatterplot(data=dataframe, x='real f', y="SMI p2", label="SMI p2", color='black', s=40, ax=axes[1])
    sns.lineplot(data=dataframe, x='real f', y="real p2", label="GT p2", color='grey', ax=axes[1])
    axes[1].set_title("Estimated p2 with SMI vs Real f")
    axes[1].set_xlabel("Real f")
    axes[1].set_ylabel("SMI p2")
    axes[1].legend()


    if not cyl:
        # 2. Combine "SMI Da", "SMI Depar", and "SMI Deperp" in one subplot
        sns.scatterplot(data=dataframe, x='real f', y="SMI Da", label="SMI Da (·10⁻⁹ m²/s)", color='skyblue', s=40, ax=axes[2])
        sns.lineplot(data=dataframe, x='real f', y="DKI Da", label="DKI Da (·10⁻⁹ m²/s)", color='blue', ax=axes[2])

        sns.scatterplot(data=dataframe, x='real f', y="SMI Depar", label="SMI De∥ (·10⁻⁹ m²/s)", color='lightgreen', s=40, ax=axes[2])
        sns.lineplot(data=dataframe, x='real f', y="DKI Depar", label="DKI De∥ (·10⁻⁹ m²/s)", color='green', ax=axes[2])

        sns.scatterplot(data=dataframe, x='real f', y="SMI Deperp", label="SMI De⊥ (·10⁻⁹ m²/s)", color='#FFB347', s=40, ax=axes[2])
        sns.lineplot(data=dataframe, x='real f', y="DKI Deperp", label="DKI De⊥ (·10⁻⁹ m²/s)", color='orange', ax=axes[2])

        axes[2].set_title("Estimated Da, Depar, Deperp with SMI vs Real f")
        axes[2].set_xlabel("Real f")
        axes[2].set_ylabel("Values")
        axes[2].legend()
    else:
        # 2. Combine "SMI Da", "SMI Depar", and "SMI Deperp" in one subplot
        sns.scatterplot(data=dataframe, x='real f', y="SMI Da", label="SMI Da (·10⁻⁹ m²/s)", color='skyblue', s=40, ax=axes[2])
        sns.lineplot(data=dataframe, x='real f', y="real Da", label="GT Da (·10⁻⁹ m²/s)", color='blue', ax=axes[2])

        sns.scatterplot(data=dataframe, x='real f', y="SMI Depar", label="SMI De∥ (·10⁻⁹ m²/s)", color='lightgreen', s=40, ax=axes[2])
        sns.lineplot(data=dataframe, x='real f', y="real Depar", label="GT De∥ (·10⁻⁹ m²/s)", color='green', ax=axes[2])

        sns.scatterplot(data=dataframe, x='real f', y="SMI Deperp", label="SMI De⊥ (·10⁻⁹ m²/s)", color='#FFB347', s=40, ax=axes[2])
        sns.lineplot(data=dataframe, x='real f', y="DKI Deperp", label="De⊥ (·10⁻⁹ m²/s)", color='orange', ax=axes[2])

        axes[2].set_title("Estimated Da, Depar, Deperp with SMI vs Real f")
        axes[2].set_xlabel("Real f")
        axes[2].set_ylabel("Values")
        axes[2].legend()

    # Adjust layout
    plt.tight_layout()
    plt.show()



def plot_perm(directory):
    # Load the data
    filenames, SMI_f_values, real_f_values, SMI_Das, SMI_Depars, SMI_Deperps, SMI_p2s, WMTI_f_values, WMTI_Das, WMTI_Depars, WMTI_Deperps, WMTI_c2s, real_Das, real_Depars, real_Deperps, real_c2s = read_SMI_files(directory)

    dataframe = pd.DataFrame({
        'f': SMI_f_values,
        'GT f': real_f_values,
        'Da (·10⁻⁹ m²/s)': SMI_Das,
        'De∥ (·10⁻⁹ m²/s)': SMI_Depars,
        'De⊥ (·10⁻⁹ m²/s)': SMI_Deperps,
        'p2': SMI_p2s,
        'DKI Da (·10⁻⁹ m²/s)': real_Das,
        'DKI De∥ (·10⁻⁹ m²/s)': real_Depars,
        'DKI De⊥ (·10⁻⁹ m²/s)': real_Deperps,
        'GT c2': real_c2s,
        "GT p2": [(3 * c - 1) / 2 for c in real_c2s],
        'filename': filenames,
    })

    # Add permeability as a column
    dataframe['permeability'] = [int(filename.split("_")[1]) for filename in filenames]

    # Melt dataframe for combined metrics
    melted_combined = pd.melt(
        dataframe, 
        id_vars="permeability", 
        value_vars=["Da (·10⁻⁹ m²/s)", "De∥ (·10⁻⁹ m²/s)", "De⊥ (·10⁻⁹ m²/s)"], 
        var_name="Metric", 
        value_name="Value"
    )

    # Set up the figure with subplots
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))

    # Subplot 1: SMI f vs Real f
    sns.stripplot(data=dataframe, x='permeability', y='f', color='black', size=6, label="f", ax=axes[0])
    sns.lineplot(data=dataframe, x='permeability', y='GT f', color='grey', label="gt f", ax=axes[0])
    axes[0].set_title("SMI f and Real f vs Permeability")
    axes[0].set_xlabel("Permeability")
    axes[0].set_ylabel("Value")
    axes[0].legend()

    # Subplot 2: SMI p2 vs Real p2
    sns.stripplot(data=dataframe, x='permeability', y='p2', color='black', size=6,  label="p2", ax=axes[1])
    sns.lineplot(data=dataframe, x='permeability', y='GT p2', color='grey',  label="GT p2", ax=axes[1])
    axes[1].set_title("SMI p2 and Real p2 vs Permeability")
    axes[1].set_xlabel("Permeability")
    axes[1].set_ylabel("Value")
    axes[1].legend()

    # Subplot 3: Combined SMI Da, Depar, and Deperp
    colours = ["skyblue", "lightgreen", "#FFB347"]
    sns.scatterplot(
        data=melted_combined, 
        x="permeability", 
        y="Value", 
        hue="Metric", 
        s=40, 
        ax=axes[2],
        palette=colours
    )
    axes[2].set_title("Comparison of SMI Da, Depar, and Deperp")
    axes[2].set_xlabel("Permeability")
    axes[2].set_ylabel("Value")
    axes[2].legend(title="Metric")

    # Adjust layout for better spacing
    plt.tight_layout()
    plt.show()


def apply_SMI(directory, path_scheme, separate_intra_extra):
    
    icvf =None
    c2 = None
    if separate_intra_extra:
        files = get_files_from_folder(directory, False)
    else:
        files = get_files_from_folder(directory, True)

    for file in files:
        if separate_intra_extra and "_extra_intra" in file :
            print(file)
            # read file
            SMI_f, SMI_Da, SMI_Dpar, SMI_Dperp, SMI_p2= calculate_SMI(file, path_scheme, binary = False)
            WMTI_f, WMTI_Da, WMTI_Depar, WMTI_Deperp, WMTI_c2 = WMTI_from_file(file, path_scheme, binary = False)
            file_ = file.replace("_extra_intra", "")
            file_ = file_.split("_DWI")[0]
            with open(f"{file_}_SMI_output.txt", "a") as file:
                file.write(f"SMI f : {SMI_f}\n")
                file.write(f"SMI Da : {SMI_Da}\n")
                file.write(f"SMI Deperp : {SMI_Dperp}\n")
                file.write(f"SMI Depar : {SMI_Dpar}\n")
                file.write(f"SMI p2 : {SMI_p2}\n")
                file.write(f"WMTI f : {WMTI_f}\n")
                file.write(f"WMTI Da : {WMTI_Da}\n")
                file.write(f"WMTI Deperp : {WMTI_Deperp}\n")
                file.write(f"WMTI Depar : {WMTI_Depar}\n")
                file.write(f"WMTI c2 : {WMTI_c2}\n")
        elif not separate_intra_extra and "DWI.bfloat" in file:
            print(file)
            for file_ in get_files_from_folder(directory, False):
                if "info.txt" in file_ and "simulation" not in file_:
                    icvf = get_total_icvf(file_)
                    c2 = get_c2(file_)
                else :
                    continue
            SMI_f, SMI_Da, SMI_Dpar, SMI_Dperp, SMI_p2= calculate_SMI(file, path_scheme, binary = True)
            WMTI_f, WMTI_Da, WMTI_Depar, WMTI_Deperp, WMTI_c2 = WMTI_from_file(file, path_scheme, binary = True)
            file_ = file.replace("_extra_intra", "")
            file_ = file_.split("_DWI")[0]
            with open(f"{file_}_SMI_output.txt", "a") as file:
                file.write(f"SMI f : {SMI_f}\n")
                file.write(f"SMI Da : {SMI_Da}\n")
                file.write(f"SMI Deperp : {SMI_Dperp}\n")
                file.write(f"SMI Depar : {SMI_Dpar}\n")
                file.write(f"SMI p2 : {SMI_p2}\n")
                file.write(f"WMTI f : {WMTI_f}\n")
                file.write(f"WMTI Da : {WMTI_Da}\n")
                file.write(f"WMTI Deperp : {WMTI_Deperp}\n")
                file.write(f"WMTI Depar : {WMTI_Depar}\n")
                file.write(f"WMTI c2 : {WMTI_c2}\n")
                file.write(f" real f : {icvf}\n")
                file.write(f" real c2 : {c2}\n")


def delete_all_SMI(directory):
    for file in get_files_from_folder(directory, False):
        if "SMI_output" in file:
            os.remove(file)

def delete_all_extra_intra_files(directory):
    for file in get_files_from_folder(directory, True):
        if "extra_intra" in file:
            os.remove(file)
    for file in get_files_from_folder(directory, False):
        if "extra_intra" in file:
            os.remove(file)

def compute_slope_scipy(x, y):
    """
    Compute the slope of the best-fit line between two time series using scipy.

    Parameters:
        x (array-like): First time series (independent variable).
        y (array-like): Second time series (dependent variable).

    Returns:
        float: The slope of the best-fit line.
    """
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return slope

def plot_SMI_corrs(directory):
    
    SMI_f_values, real_f_values, SMI_Das, SMI_Depars, SMI_Deperps, SMI_p2s, WMTI_f_values, WMTI_Das, WMTI_Depars, WMTI_Deperps, WMTI_c2s, real_Das, real_Depars, real_Deperps, real_c2s = read_SMI_files(directory)

    dataframe = pd.DataFrame()
    dataframe['SMI f'] = SMI_f_values
    dataframe['real f'] = real_f_values
    dataframe['SMI Da'] = SMI_Das
    dataframe['SMI Depar'] = SMI_Depars
    dataframe['SMI Deperp'] = SMI_Deperps
    dataframe['SMI p2'] = SMI_p2s
    dataframe['real Da'] = real_Das
    dataframe['real Depar'] = real_Depars
    dataframe['real Deperp'] = real_Deperps
    dataframe['WMTI f'] = WMTI_f_values
    dataframe['WMTI Da'] = WMTI_Das
    dataframe['WMTI Depar'] = WMTI_Depars
    dataframe['WMTI Deperp'] = WMTI_Deperps
    dataframe['WMTI c2'] = WMTI_c2s
    dataframe['real c2'] = real_c2s
    dataframe['real p2'] = [(3*c-1)/2 for c in real_c2s]

    # sort dataframe by real f
    dataframe = dataframe.sort_values(by='real f')
    columns_x = [col for col in dataframe.columns if "real" in col]
    columns_y = [col for col in dataframe.columns if "SMI" in col]

    x = dataframe[columns_x]
    y = dataframe[columns_y]
    correlation_matrix = np.zeros((len(columns_x), len(columns_y)))
    for e_x, col in enumerate(columns_x):
        for e_y, col_ in enumerate(columns_y):
            mean_timeseries_x = np.array(x[col])
            #normalise
            mean_timeseries_x = (mean_timeseries_x - np.min(mean_timeseries_x)) / (np.max(mean_timeseries_x) - np.min(mean_timeseries_x))
            mean_timeseries_y = np.array(y[col_])
            #normalise
            mean_timeseries_y = (mean_timeseries_y - np.min(mean_timeseries_y)) / (np.max(mean_timeseries_y) - np.min(mean_timeseries_y))
            # linear slope between the two timeseries
            correlation_matrix[e_x, e_y] = compute_slope_scipy(mean_timeseries_x, mean_timeseries_y)

    sns.heatmap(correlation_matrix, annot=True, xticklabels=columns_y, yticklabels=columns_x, cmap='coolwarm', center =1)
    plt.show()




if __name__ == "__main__":
    
    directory= f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/permeable_axons/"
    path_scheme = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/SMI.scheme"
    permeability = True

    if not permeability:   
        delete_all_SMI(directory)
        delete_all_extra_intra_files(directory)
        assemble_intra_extra(directory, path_scheme)
        apply_SMI(directory, path_scheme, separate_intra_extra = True)
        plot_SMI_f_vs_real_f(directory, False)
        plot_SMI_corrs(directory)
    else:
        delete_all_SMI(directory)
        apply_SMI(directory, path_scheme, separate_intra_extra = False)
        plot_perm(directory)