"""
    To compare the "de-mean" signal (signal - analytical formula) 
    with and without diffusion between soma & dendrites, at different SNR 
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
if sys.platform == "linux":
    sys.path.insert(1, '/home/localadmin/Documents/analytical_formula/')
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")
from my_murdaycotts import my_murdaycotts
import math
import statannot
import nibabel as nib
from dipy.core.gradients import gradient_table
import dipy.reconst.dki as dki
import dipy.reconst.dti as dti
import math
import scipy
from utils import get_bvals, get_bvectors, calculate_DKI, get_dwi, get_psge, create_data, analytical_solutions

cur_path    = os.getcwd()
scheme_file = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/docs/scheme_files/PGSE_21_dir_12_b.scheme"

def create_df_all(DWI_folder, scheme_file_path):
    """
    Creates a dataframe with all the simulations together (e.g. neuron 1 + neuron 2 + ...)

    Args:
        DWI_folder (pathlib.PosixPath) : folder where all the experiences (all substrate, all repetitions) are stored
        scheme_file_path                (str) : path of the scheme file

    Returns:
        df_all_data  (pd.DataFrame) : Dataframe containing all the data needed for the plots, statistics, etc
        df_crossings (pd.DataFrame) : Dataframe containing the crossings information
        
    """

    df_all_data  = pd.DataFrame()
    df_crossings = pd.DataFrame()
    for overlap in os.listdir(DWI_folder):
        if os.path.isdir(DWI_folder / overlap):
            for subcase in os.listdir(DWI_folder / overlap):
                if os.path.isdir(DWI_folder / overlap / subcase):
                    # Iterate through the files in the folder
                    for neuron in os.listdir(DWI_folder / overlap / subcase):
                        for subdir in os.listdir(DWI_folder / overlap / subcase / neuron):
                            if os.path.isdir(DWI_folder / overlap / subcase / neuron / subdir):
                                for filename in os.listdir(DWI_folder / overlap / subcase / neuron / subdir):
                                    # Read the simulation_info.txt to have crossings information
                                    if "simu" in filename:
                                        with open(DWI_folder / overlap / subcase / neuron / subdir / filename, 'r') as file:
                                            N = int(subdir.split('_')[1])
                                            T = int(subdir.split('_')[3])
                                            # Read the file line by line
                                            for line in file:
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
                                        # Number of walkers
                                        N            = int(subdir.split('_')[1])
                                        # Number of timesteps
                                        T            = int(subdir.split('_')[3])
                                        extension    = filename.split('_')[-1].split('.')[-1]

                                        for SNR in [np.inf, 100, 60, 20]:
                                            data_one_exp = create_data(DWI_folder / overlap / subcase / neuron / subdir, SNR, name, extension, scheme_file_path)

                                            # For each b, iterate over all directions, store the data, and average them (powder-average)
                                            nb_b   = len(data_one_exp["b [ms/um²]"].unique())
                                            nb_dir = int(len(data_one_exp["x"].values) / nb_b)
                                            for i in range(nb_b):
                                                sb_so     = []
                                                for j in range(nb_dir):
                                                    sb_so.append(data_one_exp.iloc[nb_b*j + i, :]["Sb/So"])
                                                    bval = data_one_exp.iloc[nb_b*j + i, :]["b [ms/um²]"]
                                                
                                                mean_Sb_so  = np.mean(sb_so)
                                                d = {'loc': "intra", 
                                                    'N': [N], 
                                                    'T': [T], 
                                                    'Sb/So': [mean_Sb_so], 
                                                    'SNR': [str(SNR)], 
                                                    'b [ms/um²]': [bval], 
                                                    'neuron': [neuron], 
                                                    'case': [subcase]}
                                                df_avg_data = pd.DataFrame(d)
                                                df_all_data = pd.concat([df_all_data, df_avg_data])

    return df_all_data, df_crossings


MEDIUM_SIZE = 20
BIGGER_SIZE = 22

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

experience_folder = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/ISMRM24/exch")
# df_all_data, df_crossings = create_df_all(experience_folder, scheme_file)

# df_all_data.to_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/ISMRM24/exch/SNR.csv")
df_all_data = pd.read_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/ISMRM24/exch/SNR.csv")
df_all_data = df_all_data[(df_all_data["case"] == "soma_dendrites_ex") | (df_all_data["case"] == "soma_dendrites")]

df_all_data = df_all_data[df_all_data['b [ms/um²]'] > 0]

# Analytical solutions
Delta     = np.array([0.05])   # in [s]
delta     = np.array([0.0165]) # in [s]
D0        = 2.5e-9             # [m²/s]
# If we keep b=0, it becomes significant 
# (but S in b=0 equals 1 all the time so it's cheating...)
bvals     = df_all_data['b [ms/um²]'].values * 1e9 # in [s/m²]

# Calculate soma and neurite fraction
r_soma           = 10e-6 # [m]
volume_neurites  = 8784.68 # 11368.4 # 0.57um dendrite # 8784.68 # in [um³] (3 branching)
volume_soma      = 4/3 * np.pi * r_soma**3 # in [m³]
volume_soma      = volume_soma * 1e18 # in [um³]
volume_neuron    = volume_neurites + volume_soma
neurite_fraction = volume_neurites / volume_neuron # without unit
soma_fraction    = volume_soma / volume_neuron     # without unit

log = False
soma_signal, neurites_signal, both_signal = analytical_solutions(bvals, Delta, delta, r_soma, D0, log, soma_fraction, neurite_fraction)


# Subtract the analytical value of the corresponding b-value to the signal ("de-meaning")
df_all_data['Sb/So_no_mean'] = df_all_data['Sb/So'].values -  both_signal
df_all_data['SNR'] = df_all_data['SNR'].astype(str)
print(df_all_data.case.unique())
# Plotting
fig, ax = plt.subplots(1, 1, figsize=(11, 9))
sns.violinplot(data=df_all_data, 
               x='SNR', 
               y='Sb/So_no_mean', 
               hue='case', 
               ax=ax, 
               order=['inf', '100.0', '60.0', '20.0'])

sns.stripplot(data=df_all_data, 
              x='SNR', 
              y='Sb/So_no_mean', 
              hue='case', 
              order=['inf', '100.0', '60.0', '20.0'],
              dodge=True) 

ax.legend().set_visible(False)

for collection in ax.collections:
    if isinstance(collection, matplotlib.collections.PolyCollection):
        collection.set_edgecolor(collection.get_facecolor())
        collection.set_facecolor(collection.get_facecolor())
        collection.set_alpha(0.5)

plt.ylabel("Sb/So - analytical signal (soma + dendrites)")
plt.axhline(y=0, linestyle="dashed", color='gray')
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[:int(len(labels)/2)], ['Soma-dendrites (disconnected)', 'Soma-dendrites (connected)'], loc='upper left', frameon=False)

couples = []
couples_end = []
for b in df_all_data['SNR'].unique():
    for i, branch in enumerate(df_all_data['case'].unique()):
        couples.append((b, branch))    
        print(b, branch, "H0 : population average is 0 (the simulated signal = the analytical one)")
        # scipy.stats.wilcoxon
        print(scipy.stats.ttest_1samp(df_all_data[(df_all_data["SNR"] == b) & (df_all_data["case"] == branch)]['Sb/So_no_mean'].values, popmean=0., nan_policy='omit'))
        print("\n")

for i in range(1, len(couples) + 1):
    if i % 2 == 0:
        couples_end.append((couples[i-2], couples[i-1]))

statannot.add_stat_annotation(
    ax,
    data=df_all_data,
    y='Sb/So_no_mean', 
    x='SNR',
    hue='case',
    hue_order=['soma_dendrites', 'soma_dendrites_ex'],
    box_pairs=couples_end,
    test="t-test_ind",
    text_format="star",
    loc="inside"
    )

plt.show()
