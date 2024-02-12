import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")
import glob
from useful_functions import read_binary_file, get_files_from_folder, dki_from_file
from overlapping_spheres_analysis import add_kurtosis_to_lists

def organise_data_in_lists(concentrations, timesteps, locations, file):

    if "intra" in file :
        locations.append("intra")
    elif "extra" in file:
        locations.append("extra")
    else:   
        print("Error, no location found")
    if "conc_1000000000" in file:
        concentrations.append(1000000000)
    elif "conc_100000000" in file:
        concentrations.append(100000000)
    elif "conc_10000000" in file:
        concentrations.append(10000000)
    else:
        print("Error, no concentration found")
 
    if "timesteps_3000" in file:
        timesteps.append(3000)
    elif "timesteps_6000" in file:
        timesteps.append(6000)
    elif "timesteps_12000" in file:
        timesteps.append(12000)
    else:
        print("Error, no timestep found")
    
    return concentrations, timesteps, locations

def plot_results(df, location):  
    df = df[df["Locations"] == location]
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)
    print(df)
    g = sns.boxplot( data=df,   
                        x="Timesteps", y="RD", hue="Timesteps")

    plt.show()

def main():
    directory = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/instructions/sanity_check/"
    scheme_file = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/docs/scheme_files/PGSE_sample_scheme_21_dir.scheme"

    binary_files = get_files_from_folder(directory, binary=True)

    concentrations =[]
    timesteps = []
    locations = []
    FAs = []
    MDs = []
    ADs = []
    RDs = []
    MKs = []
    AKs = []
    RKs = []

    for file in binary_files:
        if "img" not in file and "info" not in file:
            concentrations, timesteps, locations = organise_data_in_lists(concentrations, timesteps, locations, file)
            FA, MD, AD, RD, MK, AK, RK = dki_from_file(file, scheme_file)
            MDs, ADs, RDs, MKs, AKs, RKs, FAs = add_kurtosis_to_lists(FA, MD, AD, RD, MK, AK, RK, MDs, ADs, RDs, MKs, AKs, RKs, FAs)

    df = pd.DataFrame({"Concentration": concentrations, "Timesteps": timesteps, "Locations": locations})
    df["FA"] = FAs
    df["MD"] = MDs
    df["AD"] = ADs
    df["RD"] = RDs
    df["MK"] = MKs
    df["AK"] = AKs
    df["RK"] = RKs

    plot_results(df, "intra")
    plot_results(df, "extra")

main()
