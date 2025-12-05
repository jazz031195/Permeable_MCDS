"""
    File that plots the signal
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
if sys.platform == "linux":
    sys.path.insert(1, '/home/localadmin/Documents/analytical_formula/')
import pandas as pd
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
from utils import get_bvals, get_bvectors, calculate_DKI, get_dwi, get_psge, create_data, analytical_solutions

cur_path    = os.getcwd()
scheme_file = "/Users/ideriedm/Documents/MCDS_perm/Permeable_MCDS/docs/128_dir_12_b_9_td.scheme"
giro        = 2.6751525e8 # Gyromagnetic radio [rad/(s*T)]

def get_info(simu_info):
    with open(simu_info, 'r') as file:
        # Read the file line by line
        for line in file:
            if 'Number of particles:' in line:
                N = int(line.split('-')[-1])
            if 'Number of steps:' in line:
                T = int(line.split('-')[-1])
            
        d = {'N': [N], 'T': [T]}
    return pd.DataFrame(d)

def create_df_all(DWI_folder, scheme_file_path):
    """
    Creates a dataframe with all the simulations together (e.g. neuron 1 + neuron 2 + ...)
    Args:
        experience_folder (pathlib.PosixPath) : folder where all the experiences (all substrate, all repetitions) are stored
        scheme_file_path          (str) : path of the scheme file
    Returns:
        df_all_data  (pd.DataFrame) : Dataframe containing all the data needed for the plots, statistics, etc
        df_crossings (pd.DataFrame) : Dataframe containing the crossings information
        
    """
    # timestep 
    T = 15500
    
    df_all_data  = pd.DataFrame()
    df_crossings = pd.DataFrame()
    # for neuron in os.listdir(DWI_folder):
    #     if os.path.isdir(DWI_folder / neuron):
    #         # Iterate through the files in the folder
    for subdir in os.listdir(DWI_folder):
        if os.path.isdir(DWI_folder / subdir) and subdir != "output" and subdir != "old":
            for filename in os.listdir(DWI_folder / subdir):

                # Check if the filename contains "_rep_" and "DWI"
                if "DWI_img" in filename:
                    simu_info = str(DWI_folder / subdir / filename).replace("DWI_img", "simulation_info")
                    info = get_info(simu_info)

                    # Name of the experience
                    name         = ('_').join(filename.split('_')[:-1])
                    extension    = filename.split('_')[-1].split('.')[-1]
                    SNR          = np.inf
                    for SNR in [np.inf, 30]:
                        print(DWI_folder / subdir / filename)
                        data_one_exp = create_data(DWI_folder / subdir, SNR, name, extension, scheme_file)

                        data = data_one_exp.copy()
                        data["b_shell"] = data["b [ms/um²]"].round(1)
                        grouped = (
                            data
                            .groupby(["Delta [ms]", "b_shell"], as_index=False)
                            .agg({
                                "Sb/So": "mean",
                                "adc [ms/um²]": "mean",   # si tu veux garder l'ADC moyen
                                "MD": "mean",   
                                "AD": "mean",   
                                "FA": "mean",   
                                "RD": "mean",   
                                "MK": "mean",   
                                "RK": "mean",   
                                "AK": "mean",   
                            }))

                        grouped["loc"] = "intra"
                        grouped["N"]   = info["N"].values[0]
                        grouped["T"]   = info["T"].values[0]

                        if subdir.split("_")[-2] == "0":
                            conf = "_".join(subdir.split("_")[:-2]) + "_ref"
                        else:
                            conf = "_".join(subdir.split("_")[:-2])

                        grouped["configuration"] = conf
                        grouped["configuration_all"] = "_".join(subdir.split("_")[:-1])
                        grouped["SNR"] = SNR
                        # On renomme pour coller à ce que tu veux dans d
                        grouped = grouped.rename(columns={
                            "b_shell": "b [ms/um²]",
                        })
                        df_avg_data = grouped[["loc", "N", "T", "Sb/So", "b [ms/um²]", "Delta [ms]", "configuration", "configuration_all", "MD", "RD", "AD", "FA", "MK", "RK", "AK", "SNR"]]
                        df_all_data = pd.concat([df_all_data, df_avg_data], ignore_index=True)
    return df_all_data, df_crossings

MEDIUM_SIZE = 25
BIGGER_SIZE = 25
plt.rc('font',   size=MEDIUM_SIZE)       # controls default text sizes
plt.rc('axes',   titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes',   labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick',  labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick',  labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

DWI_folder = Path("/Users/ideriedm/Downloads/results_geo")
df_all_data, df_crossings = create_df_all(DWI_folder, scheme_file)

df_all_data.to_csv("/Users/ideriedm/Documents/MCDS_perm/Permeable_MCDS/results/OHBM26/data.csv")
df_all_data = pd.read_csv("/Users/ideriedm/Documents/MCDS_perm/Permeable_MCDS/results/OHBM26/data.csv")

df_all_data["Delta"] = df_all_data["Delta [ms]"]
b_labels    = df_all_data["b [ms/um²]"].unique()
Deltas      = df_all_data["Delta"].unique()
df_all_data = df_all_data[(df_all_data['b [ms/um²]'] > 0)]
volume_neuron = 1357.93 # µm3
for D in Deltas:
    fig, ax = plt.subplots(1, 1, figsize=(12,10))
    fig.subplots_adjust(top=0.7)
    df = df_all_data[(df_all_data['Delta'] == D) & ~(df_all_data["configuration"].str.contains("ref"))].copy()
    df_ref = df_all_data[(df_all_data['Delta'] == D) & (df_all_data["configuration"].str.contains("ref"))].copy()

    # Compute mean and std for each (b, configuration)
    stats = (
        df.groupby(['b [ms/um²]', 'configuration'])
        .agg(mean_S=('Sb/So', 'mean'),
            std_S=('Sb/So', 'std'))
        .reset_index()
    )

    # Compute systematic jitter *within each b*
    configs = ['neuron', 'neuron_tortuous', 'neuron_beaded', 'neuron_tortuous_beaded']
    n_configs = len(configs)

    # symmetric offsets, e.g. for 3 groups → [-0.1, 0, +0.1]
    offsets = np.linspace(-0.15, 0.15, n_configs)   # adjust 0.15 to tune spread
    offset_map = dict(zip(configs, offsets))

    # Apply jitter: offset added TO EACH b value
    # stats['b_jittered'] = stats['b [ms/um²]'] + stats['configuration'].map(offset_map)
    stats['b_jittered'] = stats['b [ms/um²]'] 
    hue_order = ["neuron", "neuron_tortuous", "neuron_beaded", "neuron_tortuous_beaded"]
    # Error bars
    ax.errorbar(
        stats['b_jittered'],
        stats['mean_S'],
        yerr=stats['std_S'],
        fmt='none',
        ecolor='gray',
        elinewidth=2,
        capsize=3
    )

    # Plot
    g = sns.scatterplot(
        data=stats,
        x='b_jittered',
        y='mean_S',
        hue='configuration',
        hue_order=hue_order,
        s=60,
        ax=ax
    )

    sns.lineplot(data=df_ref,
                 x="b [ms/um²]",
                 y="Sb/So",
                 hue='configuration',
                 hue_order=[h + "_ref" for h in hue_order],
                 ax=ax)

    ax.set_ylabel("Sb/S0")
    ax.set_ylim([0, 1])
    ax.set_xlabel("b [ms/um²]")
    handles, labels = ax.get_legend_handles_labels()

    labels = [label.replace("neuron_", "").replace("_", " ") for label in labels]
    labels[0] = "straight"
    ax.legend(handles, labels, ncol=2, bbox_to_anchor=(1, 1.5)).set_title('')
    if not os.path.exists("results/OHBM26/signal_decay/"):
        os.makedirs("results/OHBM26/signal_decay/")
    fig.savefig(f"results/OHBM26/signal_decay/signal_decay_{D}.png")
    fig.savefig(f"results/OHBM26/signal_decay/signal_decay_{D}.pdf")