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
scheme_file = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/test/test_juliette/128_dir_12_b_9_td.scheme"
giro        = 2.6751525e8 # Gyromagnetic radio [rad/(s*T)]

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
            
                # Read the simulation_info.txt to have crossings information
                if "simu" in filename:
                    with open(DWI_folder / subdir / filename, 'r') as file:
                        # Read the file line by line
                        for line in file:
                            if 'Number of particles:' in line:
                                N = float(line.split('-')[-1])
                            # Check if the line contains the relevant information
                            if 'Number of particles eliminated due crossings' in line:
                                # Split the line to get the number of particles as the last element
                                num_particles_crossings = int(line.split()[-1])
                                # Break the loop, as we have found the information we need
                                break
                        d = {'nb_crossings': [num_particles_crossings], 'N': [N], 'T': [T]}
                        df_avg_crossings = pd.DataFrame(d)
                        df_crossings     = pd.concat([df_crossings, df_avg_crossings])
                
                # Check if the filename contains "_rep_" and "DWI"
                if "DWI_img" in filename:
                    # Name of the experience
                    name         = ('_').join(filename.split('_')[:-1])
                    N = 1000
                    # Number of walkers
                    # N            = int(subdir.split('_')[0])
                    # Number of timesteps
                    #T            = int(subdir.split('_')[3])
                    extension    = filename.split('_')[-1].split('.')[-1]
                    SNR          = np.inf
                    data_one_exp = create_data(DWI_folder / subdir, SNR, name, extension, scheme_file)
                    data = data_one_exp.copy()
                    data["b_shell"] = data["b [ms/um²]"].round(1)
                    grouped = (
                        data
                        .groupby(["Delta [ms]", "b_shell"], as_index=False)
                        .agg({
                            "Sb/So": "mean",
                            "adc [ms/um²]": "mean"   # si tu veux garder l'ADC moyen
                        }))
                    grouped["loc"] = "intra"
                    grouped["N"]   = N
                    grouped["T"]   = T    # tu le définis quelque part avant
                    # On renomme pour coller à ce que tu veux dans d
                    grouped = grouped.rename(columns={
                        "b_shell": "b [ms/um²]",
                    })
                    df_avg_data = grouped[["loc", "N", "T", "Sb/So", "b [ms/um²]", "Delta [ms]"]]
                    df_all_data = pd.concat([df_all_data, df_avg_data], ignore_index=True)
    return df_all_data, df_crossings

MEDIUM_SIZE = 19
BIGGER_SIZE = 19
plt.rc('font',   size=MEDIUM_SIZE)       # controls default text sizes
plt.rc('axes',   titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes',   labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick',  labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick',  labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

DWI_folder = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/test")
df_all_data, df_crossings = create_df_all(DWI_folder, scheme_file)

df_all_data.to_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/test/test_juliette/data2.csv")
# df_all_data = pd.read_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/test/test_juliette/data.csv")
print(df_all_data)

df_all_data["Delta"] = df_all_data["Delta [ms]"]
b_labels    = df_all_data["b [ms/um²]"].unique()
Deltas      = df_all_data["Delta"].unique()
df_all_data = df_all_data[(df_all_data['b [ms/um²]'] > 0)]
volume_neuron = 1357.93 # µm3
for D in Deltas:
    fig, ax = plt.subplots(2, 5, figsize=(45,15))
    ax = ax.ravel()
    for i, b in enumerate(b_labels[1:-1]):
        sns.boxplot(data=df_all_data[(df_all_data['b [ms/um²]'] == b) & (df_all_data['Delta'] == D)], x='N', y='Sb/So', ax=ax[i], hue="T", dodge=True)
        # Get current x-axis tick labels and values
        xticks_labels = ax[i].get_xticklabels()
        # Modify the tick labels by dividing their values by 1000
        # new_xticks_labels = [int(float(label.get_text()) / 1000) for label in xticks_labels]
        new_xticks_labels = [float(float(label.get_text()) / volume_neuron) for label in xticks_labels]
        new_xticks_labels = [f"{lab:.1f}" for lab in new_xticks_labels]
        ax[i].set_title(f'b = {b:.1f} ms/um²')
        ax[i].legend().set_visible(False)
        if i not in [0, 5]:
            ax[i].set_ylabel("")
        if i < 5:
            ax[i].set_xticklabels("")
            ax[i].set_xlabel("")
        else:
            ax[i].set_xticklabels(new_xticks_labels)
            ax[i].set_xlabel("Density [Walkers / um³]")
    handles, labels = ax[0].get_legend_handles_labels()
    ph = [plt.plot([],marker="", ls="")[0]] # Canvas
    handles = ph + handles
    import re
    with open('results/ISMRM24/Benchmark/overlap4/n1/N_5000_T_5000/_simulation_info.txt', 'r') as file:
        # Read the file line by line
        for line in file:
            if 'Particle dynamics duration' in line:
                TE = float(line.split()[-2]) / 1000 # s
            if 'Intra Diffusivity' in line:
                pattern = r'-?\d*\.?\d+(?:e[-+]?\d+)?'
                # Search for the pattern in the line
                match = re.search(pattern, line)
                if match:
                    D0 = abs(float(match.group())) # m²/s
    # step_length = np.sqrt(6 * D0 * TE / T)
    labels = [np.sqrt(6 * D0 * TE / float(lab))*1e6 for lab in labels]
    labels = [f"{lab:.2f}" for lab in labels]
    plt.legend(handles, ["Step length [um] : "] + labels, bbox_to_anchor=(-2.2, 2.3), loc="lower center", borderaxespad=0., ncol=4, frameon=False)
    fig.savefig(f"benchmark {D}.png")