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
from useful_functions import read_binary_file, get_files_from_folder, get_info_files_from_path
import re

cur_path = os.getcwd()
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
scheme_file = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/docs/scheme_files/PGSE_sample_scheme_21_dir.scheme"


def organise_data_in_lists(file, locations_list, types_list, factors_list):
    """
    Organizes data from a file into separate lists based on the file name.

    Parameters:
    file (str): The name of the file.
    locations_list (list): The list to store the locations.
    types_list (list): The list to store the types.
    factors_list (list): The list to store the factors.

    Returns:
    tuple: A tuple containing the updated locations_list, types_list, and factors_list.
    """
    print(file)
    if "cylinders_" in file:
        types_list.append("cylinder")
    elif "axons_" in file:
        types_list.append("axons")
    else:
        print("Error, no type found")
    if "intra" in file :
        locations_list.append("intra")
    elif "extra" in file:   
        locations_list.append("extra") 
    else:
        print("Error, no location found")

    if "factor_2" in file:
        factors_list.append(2)
    elif "factor_4" in file:
        factors_list.append(4)
    elif "factor_8" in file:
        factors_list.append(8)
    elif "factor_16" in file:
        factors_list.append(16)
    elif "factor_32" in file:
        factors_list.append(32)
    else:
        factors_list.append(0)

    return locations_list, types_list, factors_list

def dki_from_file(file_name, scheme_file):

    dwi = read_binary_file(file_name)

    img  = array_to_nifti(dwi)

    FA, MD, AD, RD, MK, AK, RK, axial_diffusion_orientation = calculate_DKI(scheme_file, img)

    print("axial direction : ", axial_diffusion_orientation)
    return FA, MD, AD, RD, MK, AK, RK
    
def add_kurtosis_to_lists(FA, MD, AD, RD, MK, AK, RK, MDs, ADs, RDs, MKs, AKs, RKs, FAs):
    MDs.append(MD)
    ADs.append(AD)
    RDs.append(RD)
    MKs.append(MK)
    AKs.append(AK)
    RKs.append(RK)
    FAs.append(FA)

    return MDs, ADs, RDs, MKs, AKs, RKs, FAs

def create_dataframe_from_files(files, scheme_file):
    """
    Create a pandas DataFrame from a list of files.

    Args:
        path_to_files (str): The path to the directory containing the files.
        scheme_file (str): The file containing the scheme information.

    Returns:
        pandas.DataFrame: The DataFrame containing the data extracted from the files.
    """


    locations_list = []
    types_list = []
    factors_list = []
    MDs = []
    ADs = []
    RDs = []
    MKs = []
    AKs = []
    RKs = []
    FAs = []

    for file in files:
        if "img" not in file and "traj" not in file and "bhdr" not in file:
            locations_list, types_list, factors_list = organise_data_in_lists(file, locations_list, types_list, factors_list)
            FA, MD, AD, RD, MK, AK, RK = dki_from_file(file, scheme_file)
            MDs, ADs, RDs, MKs, AKs, RKs, FAs = add_kurtosis_to_lists(FA, MD, AD, RD, MK, AK, RK, MDs, ADs, RDs, MKs, AKs, RKs, FAs)

    df = pd.DataFrame(columns = ["location", "type", "factor", "FA", "MD", "AD", "RD", "MK", "AK", "RK"])
    df["location"] = locations_list
    df["type"] = types_list
    df["factor"] = factors_list
    df["FA"] = FAs
    df["MD"] = MDs
    df["AD"] = ADs
    df["RD"] = RDs
    df["MK"] = MKs
    df["AK"] = AKs
    df["RK"] = RKs

    return df

def plot_data(df, location):
    """
    Plot the data in the DataFrame.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data.
    """
    df = df.loc[df["location"] == location]
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.5)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    sns.boxplot(x="type", y="AD", hue="factor", data=df, ax=axs[0])
    sns.boxplot(x="type", y="RD", hue="factor", data=df, ax=axs[1])
    fig.suptitle(location)
    axs[0].set_ylabel("AD")
    axs[1].set_ylabel("RD")
    plt.tight_layout()
    plt.show()

def plt_simulation_time_from_files(files):
    """
    Plot the simulation time from a list of files.

    Args:
        files (list): The list of files to plot the simulation time from.
    """
    factors = []
    simulation_times = []
    for file in files:
        if "info" in file:
            if "factor_2" in file:
                factors.append(2)
            elif "factor_4" in file:
                factors.append(4)
            elif "factor_8" in file:
                factors.append(8)
            elif "factor_16" in file:
                factors.append(16)
            elif "factor_32" in file:   
                factors.append(32)
            elif "cylinders_" in file:
                factors.append(0)
            else:
                print("Error, no factor found")

            simulation_times.append(extract_simulation_time(file))

    df = pd.DataFrame(columns = ["factor", "simulation_time"])
    df["factor"] = factors
    df["simulation_time"] = simulation_times
    print(df.loc[df["factor"] == 2])
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.5)
    fig, ax = plt.subplots()
    sns.barplot(x="factor", y="simulation_time", data=df, ax=ax)
    ax.set_ylabel("Simulation time (s)")
    ax.set_xlabel("Overlapping factor (0 for cylinders)")
    plt.tight_layout()

    plt.show()



def extract_simulation_time(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if "simulations" in line:
                time = line.split(":")[1]
                time = time.split("in ")[0]
                time_minutes = int(time.split(" minutes")[0])
                time_seconds = int(time.split("and ")[1].split("seconds")[0])/60
                time = time_minutes + time_seconds
                return time
    return None




def plot_overlapping_spheres():
    path = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/instructions/axons_vs_cylinders/data/"
    path_to_files = get_files_from_folder(path)
    df = create_dataframe_from_files(path_to_files, scheme_file)
    plot_data(df, "extra")
    plot_data(df, "intra")


def plot_time_for_run():
    path = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/instructions/axons_vs_cylinders/data/"
    path_to_files = get_info_files_from_path(path)
    plt_simulation_time_from_files(path_to_files)

plot_time_for_run()