import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")
from DKI import calculate_DKI, array_to_nifti
import glob

def read_binary_file(file_name):
    """
    Reads a binary file and returns the data as a numpy array
    """
    with open(file_name, "rb") as f:
        data = np.fromfile(f, dtype=np.float32)
    return data

def get_files_from_folder(folder_path, binary=True):
    if binary:
        binary_files = glob.glob(os.path.join(folder_path, "*.bfloat"))
    else:
        binary_files = glob.glob(os.path.join(folder_path, "*.txt"))
    binary_files.sort()

    return binary_files

def dki_from_file(file_name, scheme_file):

    dwi = read_binary_file(file_name)

    img  = array_to_nifti(dwi)

    FA, MD, AD, RD, MK, AK, RK = calculate_DKI(scheme_file, img)

    return FA, MD, AD, RD, MK, AK, RK

def get_info_files_from_path(folder_path):

    files = glob.glob(os.path.join(folder_path, "*.txt"))
    return files

def plot_sphere_from_directions_in_scheme():
    """
    Plots a sphere with the directions from the scheme file.
    """
    scheme_file = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/docs/scheme_files/PGSE_sample_scheme_21_dir.scheme"
    
    scheme = pd.read_csv(scheme_file, sep=" ", header=None, skiprows=1)
    scheme = scheme.dropna(axis=1)
    scheme = scheme.to_numpy()
    gradient_strength = scheme[:, 3]
    scheme = scheme[:, :3]
    #normalise the scheme
    scheme = scheme / np.linalg.norm(scheme, axis=1)[:, None]

    scheme = scheme * gradient_strength[:, None]
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(scheme[:, 0], scheme[:, 1], scheme[:, 2])
    # label x, y, z axis
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

    scheme_file = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/docs/scheme_files/PGSE_sample_scheme_21_dir.scheme"

