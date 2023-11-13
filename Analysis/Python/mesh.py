"""
    File to compare the mesh with the analytical solutions for different 
    decimation. Plots the mean diffusivity (MD) for different decimations
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
scheme_file = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/docs/scheme_files/PGSE_21_dir_12_b.scheme"
giro        = 2.6751525e8 # Gyromagnetic radio [rad/(s*T)]

def create_df_all(experience_folder, scheme_file_path):
    """
    Creates a dataframe with all the simulations together (e.g. neuron 1 + neuron 2 + ...)

    Args:
        experience_folder (pathlib.PosixPath) : folder where all the experiences (all substrate, all repetitions) are stored
        scheme_file_path          (str) : path of the scheme file

    Returns:
        df_all_data  (pd.DataFrame) : Dataframe containing all the data needed for the plots, statistics, etc
        df_crossings (pd.DataFrame) : Dataframe containing the crossings information
        
    """

    df_all_data  = pd.DataFrame()
    df_crossings = pd.DataFrame()
    for subcase in os.listdir(experience_folder):
        # Iterate through the files in the folder
        for neuron in os.listdir(experience_folder / subcase):
            if neuron == "n1":
                for filename in os.listdir(experience_folder / subcase / neuron):
                    # Read the simulation_info.txt to have crossings information
                    if "simu" in filename:
                        N = int(filename.split('_')[1])
                        T = int(filename.split('_')[3])
                        with open(experience_folder / subcase / neuron / filename, 'r') as file:
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
                        N            = int(name.split('_')[1])
                        # Number of timesteps
                        T            = int(name.split('_')[3])
                        extension    = filename.split('_')[-1].split('.')[-1]
                        SNR          = np.inf
                        data_one_exp = create_data(experience_folder / subcase / neuron, SNR, name, extension, scheme_file_path)
                        FA = data_one_exp["FA"][0]
                        MD = data_one_exp["MD"][0]
                        AD = data_one_exp["AD"][0]
                        RD = data_one_exp["RD"][0]
                        MK = data_one_exp["MK"][0]
                        AK = data_one_exp["AK"][0]
                        RK = data_one_exp["RK"][0]
                        
                        # For each b, iterate over all directions, store the data, and average them (powder-average)
                        nb_b   = len(data_one_exp["b [ms/um²]"].unique())
                        nb_dir = int(len(data_one_exp["x"].values) / nb_b)
                        for i in range(nb_b):
                            sb_so = []
                            adc   = []
                            for j in range(nb_dir):
                                sb_so.append(data_one_exp.iloc[nb_b * j + i, :]["Sb/So"])
                                adc.append(data_one_exp.iloc[nb_b * j + i, :]["adc [ms/um²]"])
                                bval = data_one_exp.iloc[nb_b * j + i, :]["b [ms/um²]"]
             
                            # Powder-average signal
                            mean     = np.mean(sb_so)
                            # Powder-average ADC
                            mean_adc = np.mean(adc)
                            d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, 'log_Sb/So': np.log(mean), "adc [ms/um²]": mean_adc, 
                                'b [ms/um²]': bval, 'neuron': neuron, 'case': subcase,
                                'FA': FA, 'MD': MD, 'RD': RD, 'AD': AD, 'MK': MK, 'RK': RK, 'AK': AK}
                            df_avg_data = pd.DataFrame(d, index=[i])
                            df_all_data = pd.concat([df_all_data, df_avg_data])

    return df_all_data, df_crossings


log = True

if log:
    y_lim_min = -3
    y_lim_max = 1
else:
    y_lim_min = 0.
    y_lim_max = 1.1

MEDIUM_SIZE = 19
BIGGER_SIZE = 19

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

experience_folder = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/exch")
df_all_data, df_crossings = create_df_all(experience_folder, scheme_file)

b_labels  = df_all_data["b [ms/um²]"].unique()
means     = df_all_data[((df_all_data['b [ms/um²]'] > 0) & (df_all_data["T"] == 15000))].groupby(['b [ms/um²]', 'case'])['Sb/So'].mean().reset_index()

fig, ax = plt.subplots(1, 1, figsize=(17, 10))
if not log:
    g = sns.scatterplot(data=means, 
                        x='b [ms/um²]', 
                        y='Sb/So', 
                        hue='case', 
                        hue_order=['soma', 'dendrites', 'soma_dendrites', 'soma_dendrites_ex', 'mesh_005', 'mesh_010', 'mesh_025', 'mesh_050', 'mesh_075', 'mesh_100'], 
                        ax=ax, 
                        style='case',
                        s=200, 
                        palette=['b', 'orange', 'g', 'g', 'k', 'k', 'k', 'k', 'k', 'k'])
else:
    means_log = df_all_data[((df_all_data['b [ms/um²]'] > 0) & (df_all_data["T"] == 15000))].groupby(['b [ms/um²]', 'case'])['log_Sb/So'].mean().reset_index()
    g = sns.scatterplot(data=means_log, 
                    x='b [ms/um²]', 
                    y='log_Sb/So', 
                    hue='case', 
                    hue_order=['soma', 'dendrites', 'soma_dendrites', 'soma_dendrites_ex', 'mesh_005', 'mesh_010', 'mesh_025', 'mesh_050', 'mesh_075', 'mesh_100'], 
                    ax=ax, 
                    style='case',
                    s=200, 
                    palette=['b', 'orange', 'g', 'g', 'k', 'k', 'k', 'k', 'k', 'k'])
    g.set_ylabel('log(Sb/So)')
    
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, ['Soma', 'Dendrites', 'Soma-Dendrites (disconnected)', 'Soma-Dendrites (connected)', 
                    'Mesh 5%', 'Mesh 10%', 'Mesh 25%', 'Mesh 50%', 'Mesh 75%', 'Mesh 100%'], loc='upper right', title='Intra signal')


# Analytical solutions
Delta     = np.array([0.05])  # in [s]
delta     = np.array([0.0165])# in [s]
D0        = 2.5e-9 # [m²/s]
bvals     = np.linspace(0.2, 10, 100) * 1e9 # in [s/m²]

r_soma           = 10e-6 # [m]
volume_neurites  = 8784.68 # in [um³] (3 branching)
volume_soma      = 4/3 * np.pi * r_soma**3 # in [m³]
volume_soma      = volume_soma * 1e18 # in [um³]
volume_neuron    = volume_neurites + volume_soma
neurite_fraction = volume_neurites / volume_neuron
soma_fraction    = volume_soma / volume_neuron
print("soma volume {:e}".format((volume_soma*1e18)))
print("neurites volume {:e}".format((volume_neurites*1e18)))
print("neuron {:e}".format((volume_neuron*1e18)))
print("soma fraction {:e}".format(soma_fraction))

soma_signal, neurites_signal, both_signal = analytical_solutions(bvals, Delta, delta, r_soma, D0, log, soma_fraction, neurite_fraction)

ax2 = ax.twinx()
ax2.plot(bvals*1e-9, soma_signal, label=f"Soma", color='b', linestyle="dotted")
ax2.plot(bvals*1e-9, neurites_signal, label=f"Dendrites", color='orange', linestyle="dotted")
ax2.plot(bvals*1e-9, both_signal, label=f"Soma & dendrites", color='g', linestyle="dotted")
if log:
    ax2.plot(bvals*1e-9, -bvals*D0, label="Water free diffusion, D = 2.5 [ms/um²]")
ax2.legend(title='Analytical solution', loc=3)
ax2.set_yticklabels([])
ax2.set_yticks([])
ax2.set_ylim([y_lim_min, y_lim_max])
ax.set_ylim([y_lim_min, y_lim_max])
plt.show()


# Plot Mean diffusivity for several decimation
df = df_all_data[(df_all_data['b [ms/um²]'] > 0) & (df_all_data["T"] == 15000) & ((df_all_data["case"] != "soma") & (df_all_data["case"] != "dendrites") & (df_all_data["case"] != "soma_dendrites") & (df_all_data["case"] != "soma_dendrites_ex"))]
fig, ax = plt.subplots(1,1, figsize=(10, 7))
g = sns.boxplot(data=df, 
                x='case', 
                y='MD', 
                palette="Blues", 
                order=['mesh_005', 'mesh_010', "mesh_025", 'mesh_050', 'mesh_075', 'mesh_100'])
# 5% decimation means that we keep only 5% of the faces
# 100% decimation means no decimation (original mesh)
labels = ['5%', '10%', '25%', '50%', '75%', '100%']
ax.set_xticklabels(labels)
ax.set_ylabel('Mean diffusivity')
ax.set_ylim([0.69, 0.74])
ax.set_xlabel('Mesh decimation')
plt.show()