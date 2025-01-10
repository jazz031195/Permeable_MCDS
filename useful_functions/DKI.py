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
from useful_functions import read_and_extract_parameters, read_binary_file, array_to_nifti, get_files_from_folder, read_and_extract_dwi

def txt_to_nifti(path_to_DWI):
    dwi_array = read_and_extract_dwi(path_to_DWI)

    # Create an empty 4x4 affine matrix with ones on the diagonal
    affine = np.eye(4)

    img = nib.Nifti1Image(dwi_array, affine)
    return img


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



def calculate_DTI(dwi, bvals, bvecs, TE):

    print("TE : ", TE)

    folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/DKI_stuff"
    b_values = bvals[bvals<1001]
    vectors = bvecs[bvals<1001]
    dwi = dwi[bvals<1001]
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


def dki_from_file(file_name, scheme_file, binary= True):

    if binary:
        dwi = read_binary_file(file_name)
    else:
        dwi = read_and_extract_dwi(file_name, binary = False)
    bvals, _, _, _, bvecs, TE = read_and_extract_parameters(scheme_file)

    FA, MD, AD, RD, MK, AK, RK = calculate_DKI(dwi, bvals, bvecs, TE)

    return FA, MD, AD, RD, MK, AK, RK

def dti_from_file(file_name, scheme_file, binary= True):

    if binary:
        dwi = read_binary_file(file_name)
    else:
        dwi = read_and_extract_dwi(file_name, binary = False)

    bvals, _, _, _, bvecs, TE = read_and_extract_parameters(scheme_file)

    MD, AD, RD= calculate_DTI(dwi, bvals, bvecs, TE)

    return MD, AD, RD