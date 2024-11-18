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
from useful_functions import read_and_extract_parameters, read_binary_file, array_to_nifti


def get_files_from_folder(folder_path, binaray = True):
    if binaray:
        txt_files = glob.glob(os.path.join(folder_path, f"*.bfloat"))
    else:
        txt_files = glob.glob(os.path.join(folder_path, f"*.txt"))
    txt_files.sort()

    return txt_files

def read_and_extract_bvals(file_path):
    column_4 = []  # Initialize an empty list to store the values from the 4th column
    try:
            with open(file_path, 'r') as file:
                for line in file:
                    # Split each line into space-separated values and extract the 4th column
                    columns = line.strip().split()
                    if len(columns) >= 5:
                        G = float(columns[3]) # Gradient strength given in T/m
                        giro = 2.6751525e5 #Gyromagnetic radio given in rad/(ms*T)
                        delta = float(columns[5])
                        Delta = float(columns[4])
                        b = pow(G * giro * delta, 2) * (Delta - delta/ 3) 
                        column_4.append(b)  # Assuming columns are 0-based

    except FileNotFoundError:
        print(f"File not found: {file_path}")
        
    return np.array(column_4)



def read_and_extract_dwi(file_path, binary = True):
    column = []  # Initialize an empty list to store the values 
    if binary:
        lines = np.fromfile(file_path, dtype="float32")
        for line in lines:
  
            column.append(float(line))  # Assuming columns are 0-based
    else:
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    # Split each line into space-separated values and extract 
                    columns = line.strip().split()
                    column.append(float(columns[0]))  # Assuming columns are 0-based

        except FileNotFoundError:
            print(f"File not found: {file_path}")
    
    return np.array(column)

def extract_TE(file_path):
    with open(file_path, 'r') as file:
        for e, line in enumerate(file):
            if e != 0:
                # Extract the value after the label
                TE = float(line.split()[-1])
                return TE
    return None  # Return None if "TE" not found

def read_and_extract_bvecs(file_path):
    columns_1_to_3 = []  # Initialize an empty list to store the values from the first three columns

    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Split each line into space-separated values and extract the first three columns
                columns = line.strip().split()[:3]  # Assuming columns are 0-based
                if len(columns) == 3:
                    b_value = float(line.strip().split()[3])
                    if b_value != 0:
                        columns = [float(val) for val in columns]  # Convert to float if needed
                        columns_1_to_3.append(columns)
                    else :

                        columns_1_to_3.append([0,0,0])

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    
    return np.array(columns_1_to_3)

def txt_to_nifti(path_to_DWI):
    dwi_array = read_and_extract_dwi(path_to_DWI)

    # Create an empty 4x4 affine matrix with ones on the diagonal
    affine = np.eye(4)

    img = nib.Nifti1Image(dwi_array, affine)
    return img

def dki_from_file(file_name, scheme_file, binary= True):

    if binary:
        dwi = read_binary_file(file_name)
    else:
        dwi = read_and_extract_dwi(file_name, binary = False)

    bvecs = read_and_extract_bvecs(scheme_file)
    bvals = read_and_extract_bvals(scheme_file)
    TE = extract_TE(scheme_file)

    FA, MD, AD, RD, MK, AK, RK = calculate_DKI(dwi, bvals, bvecs, TE)

    return FA, MD, AD, RD, MK, AK, RK

def dti_from_file(file_name, scheme_file, binary= True):

    if binary:
        dwi = read_binary_file(file_name)
    else:
        dwi = read_and_extract_dwi(file_name, binary = False)

    bvecs = read_and_extract_bvecs(scheme_file)
    bvals = read_and_extract_bvals(scheme_file)
    TE = extract_TE(scheme_file)

    MD, AD, RD= calculate_DTI(dwi, bvals, bvecs, TE)

    return MD, AD, RD
def array_to_nifti_replicated(dwi_array):

    #replicate the array so that it has a shape of (2,2, 2,shape of dwi_array)
    dwi_array = np.tile(dwi_array, (2,2,2,1))

    print(dwi_array.shape)

    # Create affine
    affine = np.eye(4)

    # Ensure the array is of a numerical data type
    dwi_array = dwi_array.astype(np.float32)

    img = nib.Nifti1Image(dwi_array, affine)
    return img


def extract_and_save_dirs_bvals_DWI(paths_to_DWI, path_scheme, path_to_folder):
    
        bvals = read_and_extract_bvals(path_scheme)

        #save as text file
        np.savetxt(f"{path_to_folder}/pgse_21_dir.bval", bvals)
    
        bvecs = read_and_extract_bvecs(path_scheme)

        #save as text file
        np.savetxt(f"{path_to_folder}/pgse_21_dir.bvec", bvecs)

        nifti_array = np.zeros((2,2,2,len(bvecs)))

        print(nifti_array.shape)

        for i, path_to_DWI in enumerate(paths_to_DWI):
            #read .bfloat file
            dwi = read_and_extract_dwi(path_to_DWI, binary = True)
            if "1" in path_to_DWI:
                dwi = dwi[:252]
            if i <2:
                nifti_array[i,0,0] = dwi/50000
            else:
                nifti_array[i-2,1,0] = dwi/50000

        
        nifti_array = array_to_nifti(nifti_array)

        #save as nifti file
        nib.save(nifti_array, f"{path_to_folder}/dwi_voxel_0.nii")
    


def calculate_DTI(dwi, bvals, bvecs, TE):

    print("TE : ", TE)

    folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff"
    b_values = bvals[bvals<2999]
    vectors = bvecs[bvals<2999]
    dwi = dwi[bvals<2999]
    dwi = np.array([float(format(val, "f")) for val in dwi])

    # delete all files in /home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff/parameters_DKI
    for file in glob.glob(f"{folder}/parameters_DTI/*"):
        os.remove(file)

    # Write the integers back to bvals.txt (optional, if needed)
    with open(f"{folder}/bvals.txt", "w") as f:
        f.write("\n".join(map(str, b_values)) + "\n")

    # vectors to txt
    np.savetxt(f"{folder}/vectors.txt", vectors)

    # save as nifti
    img = array_to_nifti_replicated(dwi)
    # save img
    nib.save(img, f"{folder}/dwi.nii.gz")
    # Convert NIfTI to MIF
    os.system(f"mrconvert {folder}/dwi.nii.gz {folder}/dwi.mif -force")
    os.system(f"sudo chown -R localadmin:localadmin {folder}/parameters_DTI")
    # Add bval information to the MIF file
    os.system(f"mrconvert {folder}/dwi.mif {folder}/dwi_with_bvals.mif -fslgrad {folder}/vectors.txt {folder}/bvals.txt -force")
    #os.system(f"mrconvert {folder}/dwi.mif {folder}/dwi_with_bvals.mif -bvec {folder}/vectors.txt")
    os.system(f"docker run -v /home/localadmin/Documents/MCDS/Permeable_MCDS/output:/data nyudiffusionmri/designer2:main tmi -echo_time {TE} -sigma /data/DKI_stuff/sigmas/sigma.nii.gz -DTI /data/DKI_stuff/dwi_with_bvals.mif /data/DKI_stuff/parameters_DTI")

    for file in glob.glob(f"{folder}/parameters_DKI/*"):
        name = file.split("/")[-1]
        if "md" in name:
            MD = nib.load(file).get_fdata()
        elif "ad" in name:
            AD = nib.load(file).get_fdata()
        elif "rd" in name:
            RD = nib.load(file).get_fdata()

    return MD[0][0][0], AD[0][0][0], RD[0][0][0]



def calculate_DKI(dwi, bvals, bvecs, TE):

    print("TE : ", TE)

    folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff"
    b_values = bvals[bvals<=3000]
    vectors = bvecs[bvals<=3000]
    dwi = dwi[bvals<=3000]
    dwi = np.array([float(format(val, "f")) for val in dwi])

    # delete all files in /home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff/parameters_DKI
    for file in glob.glob(f"{folder}/parameters_DKI/*"):
        os.remove(file)

    # Write the integers back to bvals.txt (optional, if needed)
    with open(f"{folder}/bvals.txt", "w") as f:
        f.write("\n".join(map(str, b_values)) + "\n")

    # vectors to txt
    np.savetxt(f"{folder}/vectors.txt", vectors)

    # save as nifti
    img = array_to_nifti_replicated(dwi)
    # save img
    nib.save(img, f"{folder}/dwi.nii.gz")
    # Convert NIfTI to MIF
    os.system(f"mrconvert {folder}/dwi.nii.gz {folder}/dwi.mif -force")
    os.system("sudo chown -R localadmin:localadmin /home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff/parameters_DKI")
    # Add bval information to the MIF file
    os.system(f"mrconvert {folder}/dwi.mif {folder}/dwi_with_bvals.mif -fslgrad {folder}/vectors.txt {folder}/bvals.txt -force")
    #os.system(f"mrconvert {folder}/dwi.mif {folder}/dwi_with_bvals.mif -bvec {folder}/vectors.txt")
    os.system(f"docker run -v /home/localadmin/Documents/MCDS/Permeable_MCDS/output:/data nyudiffusionmri/designer2:main tmi -echo_time {TE} -sigma /data/DKI_stuff/sigmas/sigma.nii.gz -DKI /data/DKI_stuff/dwi_with_bvals.mif /data/DKI_stuff/parameters_DKI")

    for file in glob.glob(f"{folder}/parameters_DKI/*"):
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


    return FA[0][0][0], MD[0][0][0], AD[0][0][0], RD[0][0][0], MK[0][0][0], AK[0][0][0], RK[0][0][0]


def calculate_DKI_WMTI_extra_intra(path_to_DWI_extra, path_to_DWI_intra, path_scheme):

    bs, Gs, deltas, Deltas, vectors, TEs = read_and_extract_parameters(path_scheme)
    folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff"

    
    b_values = bs[bs<3000]
    vectors = vectors[bs<3000]
    TEs = TEs[bs<3000]

    # Write the integers back to bvals.txt (optional, if needed)
    with open(f"{folder}/bvals.txt", "w") as f:
        f.write("\n".join(map(str, b_values)) + "\n")

    # vectors to txt
    np.savetxt(f"{folder}/vectors.txt", vectors)

    dwi_intra = read_and_extract_dwi(path_to_DWI_intra)
    dwi_extra = read_and_extract_dwi(path_to_DWI_extra)
    # addition of extra and intra
    dwi = dwi_intra*0.1 + dwi_extra
    dwi = dwi[bs<3000]

    # save as nifti
    img = array_to_nifti_replicated(dwi)
    # save img
    nib.save(img, f"{folder}/dwi.nii.gz")
    # Convert NIfTI to MIF
    os.system(f"mrconvert {folder}/dwi.nii.gz {folder}/dwi.mif -force")
    os.system("sudo chown -R localadmin:localadmin /home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff/parameters")
    # Add bval information to the MIF file
    os.system(f"mrconvert {folder}/dwi.mif {folder}/dwi_with_bvals.mif -fslgrad {folder}/vectors.txt {folder}/bvals.txt -force")
    #os.system(f"mrconvert {folder}/dwi.mif {folder}/dwi_with_bvals.mif -bvec {folder}/vectors.txt")
    os.system(f"docker run -v /home/localadmin/Documents/MCDS/Permeable_MCDS/output:/data nyudiffusionmri/designer2:main tmi -DKI -fit_constraints 1,1,1 /data/DKI_stuff/dwi_with_bvals.mif /data/DKI_stuff/parameters")
    # Initialize the wrapper class
    wmti = WMTI_Watson(f"{folder}/parameters/", params='invivo', nodes=4)
    # Fit
    wmti.fit()
    f, Da, Depar, Deperp, c2 = wmti.maps()

    return f, Da, Depar, Deperp, c2


def calculate_DKI_WMTI(path_to_DWI, path_scheme):
    bs, Gs, deltas, Deltas, vectors, TEs = read_and_extract_parameters(path_scheme)
    folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff"

    # Filter out b-values greater than 3000
    filter_mask = bs <= 3000
    b_values = bs[filter_mask] 
    vectors = vectors[filter_mask]
    TEs = TEs[filter_mask]
    dwi = read_and_extract_dwi(path_to_DWI)[filter_mask]
    mask = np.ones((1))
    # save mask
    img_mask = array_to_nifti_replicated(mask)

    nib.save(img_mask, f"{folder}/mask.nii.gz")

    # Write b-values to bvals.txt
    with open(f"{folder}/bvals.txt", "w") as f:
        f.write("\n".join(map(str, b_values)) + "\n")

    # Save vectors to vectors.txt
    np.savetxt(f"{folder}/vectors.txt", vectors)

    # Print shapes for debugging
    print("DWI shape:", dwi.shape)
    print("b-values shape:", b_values.shape)
    print("Vectors shape:", vectors.shape)
    print("Mask shape:", mask.shape)

    # Save DWI data as NIfTI
    img = array_to_nifti_replicated(dwi)
    print(img.get_fdata())
    nib.save(img, f"{folder}/dwi.nii.gz")

    print(b_values)

    # Convert NIfTI to MIF
    os.system(f"mrconvert {folder}/dwi.nii.gz {folder}/dwi.mif -force")
    os.system("sudo chown -R localadmin:localadmin /home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff/parameters_DKI")
    os.system(f"mrconvert {folder}/dwi.mif {folder}/dwi.mif -fslgrad {folder}/vectors.txt {folder}/bvals.txt -force")

    # Run the DKI fitting using the Docker container
    os.system(f"docker run -v /home/localadmin/Documents/MCDS/Permeable_MCDS/output:/data nyudiffusionmri/designer2:main tmi -DKI -fit_constraints 0,1,0 -echo_time 0.0715 -bshape 1 -mask /data/DKI_stuff/mask.nii.gz /data/DKI_stuff/dwi.mif /data/DKI_stuff/parameters_DKI")

    # Initialize the WMTI wrapper class
    wmti = WMTI_Watson(f"{folder}/parameters_DKI/", params='invivo', nodes=4)
    wmti.fit()
    f, Da, Depar, Deperp, c2 = wmti.maps()

    return f[0][0][0] , Da[0][0][0], Depar[0][0][0], Deperp[0][0][0], c2[0][0][0]


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


def calculate_SMI(path_to_DWI, path_scheme):
    # Existing code to extract parameters and DWI data
    bs, Gs, deltas, Deltas, vectors, TEs = read_and_extract_parameters(path_scheme)
    folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff"

    b_values = bs

    with open(f"{folder}/bvals.txt", "w") as f:
        f.write("\n".join(map(str, b_values)) + "\n")

    np.savetxt(f"{folder}/vectors.txt", vectors)

    dwi = read_and_extract_dwi(path_to_DWI)

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
    os.system(f"docker run -v /home/localadmin/Documents/MCDS/Permeable_MCDS/output:/data nyudiffusionmri/designer2:main tmi -SMI -echo_time {TEs[0]} -sigma /data/DKI_stuff/sigmas/sigma.nii.gz /data/DKI_stuff/dwi.mif /data/DKI_stuff/parameters_SMI")

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
        elif "fw_smi" in name:
            fw = nib.load(file).get_fdata()

    return f[0][0][0], Da[0][0][0], Depar[0][0][0], Deperp[0][0][0], p2[0][0][0], fw[0][0][0]

def add_to_data(dwi,  path_scheme, factor, repetition):


    img  = array_to_nifti(dwi)

    FA, MD, AD, RD, MK, AK, RK = calculate_DKI(path_scheme, img)

    if factor is not None:
        col_name = f"R/{factor}"
    else :
        col_name = f"cylinder"
    df = pd.DataFrame()
    df["overlapping_distance"] = [col_name]
    df["repetition"] = [repetition]
    df["MD"] = [MD]
    df["FA"] = [FA]
    df["AD"] = [AD]
    df["RD"] = [RD]
    df["MK"] = [float(MK)]
    df["AK"] = [float(AK)]
    df["RK"] = [float(RK)]

    return df


def create_df_overlapping(path_scheme):

    folder = f"/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/instructions/axons_vs_cylinders/data2/"

    files = get_files_from_folder(folder)


    files = [file for file in files if "info" not in file]
    repetitions = ["rep_00", "rep_01", "rep_02", "rep_03"]
    datas = []
    factors = [2,4,8,16,32]

    dwi =None

    for factor in factors:
        print("factor : ", factor)

        for repetition in repetitions:

            for file in files :

                if "img" not in file and "info" not in file and f"factor_{factor}_" in file and str(repetition) in file:
                    print(file)

                    dwi = read_and_extract_dwi(file)
                    
                if dwi is not None:
                    df= add_to_data(dwi,  path_scheme, factor, repetition)
                    if "intra" in file:
                        df["Location"] = ["intra"]*len(df)
                    elif "extra" in file:
                        df["Location"] = ["extra"]*len(df)
                    datas.append(df)
                dwi = None

        for file in files:
            if "rep" not in file and "img" not in file and "info" not in file and f"factor_{factor}_" in file :
                print(file)
                dwi = read_and_extract_dwi(file)
                    
            if dwi is not None:
                df= add_to_data(dwi,  path_scheme, factor, "rep_")
                if "intra" in file:
                        df["Location"] = ["intra"]*len(df)
                elif "extra" in file:
                    df["Location"] = ["extra"]*len(df)
                datas.append(df)
            dwi = None


    print("cylinder")
    for repetition in repetitions:
            for file in files :
                if "img" not in file and "info" not in file and "_factor_" not in file and str(repetition) in file:
                    print(file)
                    dwi = read_and_extract_dwi(file)
                    
                if dwi is not None:
                    df= add_to_data(dwi,  path_scheme, None, repetition)
                    if "intra" in file:
                        df["Location"] = ["intra"]*len(df)
                    elif "extra" in file:
                        df["Location"] = ["extra"]*len(df)
                    datas.append(df)
                dwi = None


    for file in files:
        if "rep" not in file and "img" not in file and "info" not in file and "factor" not in file :
            print(file)

            dwi = read_and_extract_dwi(file)

        if dwi is not None:
            df= add_to_data(dwi,  path_scheme, None, "rep_")
            if "intra" in file:
                df["Location"] = ["intra"]*len(df)
            elif "extra" in file:
                df["Location"] = ["extra"]*len(df)
            datas.append(df)
        dwi = None



    data = pd.concat(datas)
    for col in ["MD","FA","AD","RD"]:
        data = normalise_column(col, data)

    return data


def normalise_column(col, data):
    # Calculate the average "adc [um²/ms]" for type = "cylinders" for each combination
    data_cylinder = data[data['overlapping_distance'] == "cylinder"][[col, "Location"]]
    mean_intra = data_cylinder.loc[data_cylinder["Location"] == "intra"][col].mean()
    mean_extra = data_cylinder.loc[data_cylinder["Location"] == "extra"][col].mean()

    list(map(lambda x,location : x/mean_intra if location == "intra" else x/mean_extra, list(data[col]), list(data["Location"])))
    # Divide "adc [um²/ms]" by the average value
    title = 'normalized_'+col
    data[title] = list(map(lambda x,location : x/mean_intra if location == "intra" else x/mean_extra, list(data[col]), list(data["Location"])))
    return data



def plot():

    path ="/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/docs/scheme_files/PGSE_sample_scheme_21_dir.scheme"
    df = create_df_overlapping(path)

    # List of numerical column names to plot
    numerical_columns = [
        "normalized_RD",
        "normalized_AD",
        "normalized_MD",
        "normalized_FA"
        #"normalized_MK",
        #"normalized_AK",
        #"normalized_RK"
    ]



    # Create subplots for each numerical column
    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(16, 8))
    #fig.suptitle("Percentage of Volume Increase (Swelling) vs. Diffusion and Kurtosis Parameters", fontsize=16)

    # Flatten the axes for easy iteration
    # axes = axes.flatten()

    # Set a palette for the plots
    palette0 = sns.color_palette("Set2", len(numerical_columns))

    palette1 = sns.color_palette("crest", len(numerical_columns))

    palette2 = sns.color_palette("magma", len(numerical_columns))

    # Loop through numerical columns and create scatter plots with regression lines
    for i, column in enumerate(numerical_columns):
        location = "intra"

        if (i < 3):
            palette = palette1
        elif (i == 3):
            palette = palette0
        else:
            palette = palette2
        ax = axes[0,i]
        
        # Create a scatter plot
        sns.boxplot(data= df.loc[df["Location"]==location], x="overlapping_distance", y=column, ax=ax, color=palette[i])
        
        # Add a regression line
        #sns.regplot(data=df, x="overlapping_distance", y=column, ax=ax, color=palette[i], scatter=False)
        title_name = f"{location} :{column}"
        ax.set_title(title_name)
        ax.set_xlabel("overlapping_distance")
        ax.set_ylabel(column)

    for i, column in enumerate(numerical_columns):
        location = "extra"
        if (i < 3):
            palette = palette1
        elif (i == 3):
            palette = palette0
        else:
            palette = palette2
        ax = axes[1,i]
        
        # Create a scatter plot
        sns.boxplot(data= df.loc[df["Location"]=="extra"], x="overlapping_distance", y=column, ax=ax, color=palette[i])
        
        # Add a regression line
        #sns.regplot(data=df, x="overlapping_distance", y=column, ax=ax, color=palette[i], scatter=False)
        
        title_name = f"{location} :{column}"
        ax.set_title(title_name)
        ax.set_xlabel("overlapping_distance")
        ax.set_ylabel(column)

    # Remove any remaining empty subplots
    for j in range(len(numerical_columns), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()

def get_total_icvf(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("Total icvf"):
                # Extract the value after the label
                total_icvf = float(line.split()[-1])
                return total_icvf
    return None  # Return None if "Total icvf" not found


def SMI_on_data(path_to_folder):
    path_scheme = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/SMI.scheme"
    paths_to_DWI = get_files_from_folder(path_to_folder)
    for path in paths_to_DWI:
        if "straight_cyl" in path and "test" not in path:
            if "info" not in path and "img" not in path:
                name = path.split("/")[-1].split("_DWI")[0]
                if "rep" in name:
                    info_name = name.split("_rep")[0]+"_info.txt"
                else:
                    info_name = name+"_info.txt"
                packing = get_total_icvf(f"{path_to_folder}/{info_name}")
                #f, Da, Depar, Deperp, c2 = calculate_DKI_WMTI(path, path_scheme)
                #print("WMTI")
                #print(f" f : {f} \n Da : {Da} \n Depar : {Depar} \n Deperp : {Deperp} \n c2 : {c2}")
                #with open(f"{path_to_folder}/f_{packing}/{name}_WMTI_output.txt", "w") as file: 
                #    file.write(f" f : {f} \n Da : {Da} \n Depar : {Depar} \n Deperp : {Deperp} \n c2 : {c2} \n  real f : {packing}")
                print("SMI")
                f, Da, Depar, Deperp, p2, fw = calculate_SMI(path, path_scheme)
                print(f" f : {f} \n Da : {Da} \n Depar : {Depar} \n Deperp : {Deperp} \n p2 : {p2} \n fw : {fw}")
                # save to text file
                with open(f"{path_to_folder}{name}_SMI_output.txt", "w") as file: 
                    file.write(f" f : {f} \n Da : {Da} \n Depar : {Depar} \n Deperp : {Deperp} \n p2 : {p2} \n fw : {fw} \n real f : {packing}")

def extract_SMI_output(file_path):
    f_value = None
    real_f_value = None

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(" f :"):
                f_value = float(line.split(":")[1].strip())
            elif line.startswith(" real f :"):
                real_f_value = float(line.split(":")[1].strip())
            elif line.startswith(" Da :"):
                Da = float(line.split(":")[1].strip())
            elif line.startswith(" Depar :"):
                Depar = float(line.split(":")[1].strip())
            elif line.startswith(" Deperp :"):
                Deperp = float(line.split(":")[1].strip())
            elif line.startswith(" p2 :"):
                p2 = float(line.split(":")[1].strip())
            
    
    return f_value, real_f_value, Da, Depar, Deperp, p2

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

def plot_SMI_f_vs_real_f(f_values, real_f_values, Da, Depar, Deperp, p2s):
    plt.figure(figsize=(8, 6))
    plt.scatter(real_f_values, f_values, color='blue', label='f vs. real f')
    plt.xlabel('real f')
    plt.ylabel('f')
    plt.title('Plot of f vs. real f')
    plt.grid(True)
    # plot line y = x
    plt.plot([0, 1], [0, 1], color='red', label='y = x')
    plt.legend()
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.scatter(real_f_values, Da, color='blue', label='Da vs. real f')
    plt.xlabel('real f')
    plt.ylabel('Da')
    plt.title('Plot of Da vs. real f')
    plt.grid(True)
    plt.legend()
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.scatter(real_f_values, Depar, color='blue', label='Depar vs. real f')
    plt.xlabel('real f')
    plt.ylabel('Depar')
    plt.title('Plot of Depar vs. real f')
    plt.grid(True)
    plt.legend()
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.scatter(real_f_values, p2s, color='blue', label='P2 vs. real f')
    plt.xlabel('real f')
    plt.ylabel('P2')
    plt.title('Plot of P2 vs. real f')
    plt.grid(True)
    plt.legend()
    plt.show()

def comparison_smi_preds():
    path_to_folder_= f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/"
    folders = os.listdir(path_to_folder_)

    types = []
    f_values_ = []
    real_f_values_ = []
    Das = []
    Depars = []
    Deperps = []
    p2s = []

    for folder in folders:
        if "straight_cyl" in folder:
            path_to_folder = f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/{folder}/"
            f_values, real_f_values, Da, Depar, Deperp, p2 = process_SMI_files(path_to_folder)
            if (len(f_values) == 0):
                continue

            types.extend([folder]*len(f_values))
            f_values_.extend(f_values)
            real_f_values_.extend(real_f_values)
            Das.extend(Da)
            Depars.extend(Depar)
            Deperps.extend(Deperp)
            p2s.extend(p2)
            
    data = pd.DataFrame()
    data["type"] = types
    data["f_values"] = f_values_
    data["real_f_values"] = real_f_values_
    data["Da"] = Das
    data["Depar"] = Depars
    data["Deperp"] = Deperps
    data["p2"] = p2s
    # plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=data, hue="type", y="f_values", x = "real_f_values")
    # plot x= y line
    plt.plot([0, 1], [0, 1], color='red', label='y = x')

    plt.xlabel('real f')
    plt.ylabel('f')
    plt.title('Plot of f vs. real f')
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=data, hue="type", y="Da", x = "real_f_values")
    plt.xlabel('real f')
    plt.ylabel('Da')
    plt.title('Plot of Da vs. real f')
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=data, hue="type", y="Deperp", x = "real_f_values")
    plt.xlabel('real f')
    plt.ylabel('Deperp')
    plt.title('Plot of Deperp vs. real f')
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=data, hue="type", y="Depar", x = "real_f_values")
    plt.xlabel('real f')
    plt.ylabel('Depar')
    plt.title('Plot of Depar vs. real f')
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=data, hue="type", y="p2", x = "real_f_values")
    plt.xlabel('real f')
    plt.ylabel('P2')
    plt.title('Plot of P2 vs. real f')
    plt.grid(True)
    plt.show()


if __name__ == "__main__":

    path_to_folder_= f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/"
    folders = os.listdir(path_to_folder_)
    for folder in folders:
        path_to_folder = f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/{folder}/"
        SMI_on_data(path_to_folder)

    comparison_smi_preds()
    