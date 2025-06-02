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
from utils import get_bvals, get_bvectors, calculate_DKI, get_dwi, get_psge, analytical_solutions


cur_path    = os.getcwd()
scheme_file = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/PGSE_21_dir_12_b_6_td.scheme"
giro        = 2.6751525e8 # Gyromagnetic radio [rad/(s*T)]

def create_data(data_folder, SNR, name, extension, scheme_file_path):
    """
    Function that creates the dataframe of the data from one experience (e.g. neuron 1, 5 repetitions). 
    
    Args:
        data_folder (pathlib.PoxisPath) : path of the folder where the data are stored
        SNR                (np.float64) : Signal to noise ratio, to be added as Gaussian noise (sigma=1/SNR) on the real & imaginary part of the signal
        name                      (str) : Name of the experiment
        extension                 (str) : extension of the file name
        scheme_file_path          (str) : path of the scheme file

    Returns:
        data_dwi (pd.DataFrame) : Dataframe with columns = "x", "y", "z", "G [T/um]", "Delta [ms]", "delta [ms]", "TE [ms], "Sb/So","log(Sb/So)", 
                                                           "adc [ms/um²]", "FA", "MD", "AD", "RD", "MK", "AK", "RK"

    """

    dwi_real      = get_dwi(data_folder / f"{name}.{extension}")
    # There is an imaginary part to the signal
    if os.path.exists(data_folder / f"{name}_img.{extension}"):
        dwi_imaginary = get_dwi(data_folder / f"{name}_img.{extension}")
        dwi_no_noise  = np.sqrt(dwi_real**2 + dwi_imaginary**2)

        # Add Gaussian noise on the real and imaginary part => same as rician noise. Should converge to sqrt(pi/2)*sigma
        sigma = 1/SNR
        dwi_noise = np.sqrt((dwi_real/dwi_real[0] + np.random.randn(1, dwi_real.shape[0])*sigma)**2 
                            + (dwi_imaginary/dwi_real[0] + np.random.randn(1, dwi_real.shape[0])*sigma)**2)
    else:
        dwi_no_noise = dwi_real
        # Add Gaussian noise on the real and imaginary part => same as rician noise. Should converge to sqrt(pi/2)*sigma
        sigma = 1/SNR
        dwi_noise = (dwi_real/dwi_real[0] + np.random.randn(1, dwi_real.shape[0])*sigma)
        warnings.warn("Warning...........The signal is purely real")

    data_psge          = get_psge(scheme_file_path)
    Sb_So              = list(np.squeeze(dwi_noise.reshape((-1, 1))))

    data_psge["Sb/So"] = Sb_So
    data_psge["SNR"]   = [SNR] * len(Sb_So)
        
    # Number of b-values
    nb_b     = len(data_psge["b [ms/um²]"].round().unique())
    # Number of directions
    nb_td    = 6
    nb_dir   = int(len(dwi_noise.reshape((-1, 1))) / (nb_b * nb_td))
    # Data with all the directions
    data_dwi = pd.DataFrame()
    for i in range(nb_dir * nb_td):
        # Data for one direction
        data_dir   = data_psge.iloc[i * nb_b : (i + 1) * nb_b]
        b0         = list(data_dir["b [ms/um²]"])[0]
        data_dir["log(Sb/So)"]   = list(map(lambda Sb : np.log(Sb), list(data_dir["Sb/So"])))
        adc                      = list(map(lambda b,Sb : -np.log(Sb)/(b-b0) if b != b0 else np.nan, list(data_dir["b [ms/um²]"]), list(data_dir["Sb/So"])))
        data_dir["adc [ms/um²]"] = adc

        data_dwi = pd.concat([data_dwi, data_dir])

    return data_dwi

def create_df_all(DWI_folder, scheme_file_path, deltas):
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
    for funnel in os.listdir(DWI_folder):
        if os.path.isdir(DWI_folder / funnel / "overlap4"):
            for subcase in os.listdir(DWI_folder / funnel / "overlap4"):
                # Iterate through the files in the folder
                if os.path.isdir(DWI_folder / funnel / "overlap4" / subcase):
                    for neuron in os.listdir(DWI_folder / funnel / "overlap4" / subcase):
                        for filename in os.listdir(DWI_folder / funnel / "overlap4" / subcase / neuron):
                            # Read the simulation_info.txt to have crossings information
                            if "simu" in filename:
                                with open(DWI_folder / funnel / "overlap4" / subcase /neuron / filename, 'r') as file:
                                    # Read the file line by line
                                    for line in file:
                                        if ' Number of particles:' in line:
                                            N = int(line.split('-')[-1][:-1])
                                        if 'Number of steps:' in line:
                                            T = int(line.split('-')[-1])
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
                                extension    = filename.split('_')[-1].split('.')[-1]
                                SNR          = np.inf
                                data_one_exp = create_data(DWI_folder / funnel / "overlap4" / subcase / neuron, SNR, name, extension, scheme_file_path)
                                # For each b, iterate over all directions, store the data, and average them (powder-average)
                                nb_b   = len(data_one_exp["b [ms/um²]"].round().unique())
                                nb_td  = len(deltas) 
                                nb_dir = int(len(data_one_exp["x"].values) / (nb_b * nb_td))
                
                                for i in range(nb_b):
                                    for k in range(nb_td):
                                        sb_so = []
                                        adc   = []
                                        for j in range(nb_dir):
                                            sb_so.append(data_one_exp.iloc[nb_b * j + i + k*nb_dir*nb_b, :]["Sb/So"])
                                            adc.append(data_one_exp.iloc[nb_b * j + i + k*nb_dir*nb_b, :]["adc [ms/um²]"])
                                            bval = data_one_exp.iloc[nb_b * j + i + k*nb_dir*nb_b, :]["b [ms/um²]"].round()
                                        
                                        # Powder-average signal
                                        mean     = np.mean(sb_so)
                                        # Powder-average ADC
                                        mean_adc = np.mean(adc)
                                        if "funnel" in funnel:
                                            funnel_bool = True
                                        else:
                                            funnel_bool = False
                                        d = {'loc': "intra", 'N': 50000, 'T': 15000, 'Sb/So': mean, 
                                            'b [ms/um²]': bval, 'neuron': neuron, 'case': subcase,
                                            'Delta': deltas[k], 'funnel': funnel_bool}
                                        df_avg_data = pd.DataFrame(d, index=[i])
                                        df_all_data = pd.concat([df_all_data, df_avg_data])

    return df_all_data, df_crossings


# branching = "branching"

log  = False

if log:
    y_lim_min = -5
    y_lim_max = 0.1
else:
    y_lim_min = -0.05
    y_lim_max = 0.9

MEDIUM_SIZE = 19
BIGGER_SIZE = 19

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

DWI_folder = Path("results/diff_times_single_neuron/")
deltas = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07]
# df_all_data, df_crossings = create_df_all(DWI_folder, scheme_file, deltas)

# df_all_data.to_csv(DWI_folder / "data.csv")
df_all_data = pd.read_csv(DWI_folder / "data.csv")

b_labels    = df_all_data["b [ms/um²]"].unique()
cases       = df_all_data["case"].unique()

analytical_df = pd.DataFrame()
r_soma           = 10e-6 # [m]
volume_neurites  = 8784.68 # 11368.4 # 0.57um dendrite # 8784.68 # in [um³] (3 branching)
volume_soma      = 4/3 * np.pi * r_soma**3 # in [m³]
volume_soma      = volume_soma * 1e18 # in [um³]
volume_neuron    = volume_neurites + volume_soma
neurite_fraction = volume_neurites / volume_neuron
soma_fraction    = volume_soma / volume_neuron
print("soma volume {:e}".format((volume_soma*1e18)))
print("neurites volume {:e}".format((volume_neurites*1e18)))
print("neuron {:e}".format((volume_neuron*1e18)))
print("soma fraction {:e}".format(soma_fraction))
delta     = np.array([0.0165])# in [s]
D0        = 2.5e-9 # [m²/s]
bvals     = np.linspace(1, 10, 100) * 1e9 # in [s/m²]
np.array([0.6446, 0.515386, 0.41371])
np.array([0.579155, 0.514542, 0.459714])
for td in deltas:
    # Analytical solutions
    Delta = np.array([td])  # in [s]
    soma_signal, neurites_signal, both_signal = analytical_solutions(bvals, Delta, delta, r_soma, D0, log, soma_fraction, neurite_fraction)
    d = {'case': 'soma',
         'b [ms/um²]': bvals*1e-9,
         'Delta': td,
         'Sb/So': soma_signal}
    df = pd.DataFrame(d)
    analytical_df = pd.concat([analytical_df, df])
    d = {'case': 'dendrites',
         'b [ms/um²]': bvals*1e-9,
         'Delta': td,
         'Sb/So': neurites_signal}
    df = pd.DataFrame(d)
    analytical_df = pd.concat([analytical_df, df])
    d = {'case': 'soma_dendrites',
         'b [ms/um²]': bvals*1e-9,
         'Delta': td,
         'Sb/So': both_signal}
    df = pd.DataFrame(d)
    analytical_df = pd.concat([analytical_df, df])

for case in cases:
    funnel = True

    df_tmp = df_all_data[df_all_data.funnel == funnel]
    means = df_tmp[(df_tmp['b [ms/um²]'] > 0) & (df_tmp['case'] == case)].groupby(['b [ms/um²]', 'case', 'Delta'])['Sb/So'].mean().reset_index()

    fig, ax = plt.subplots(1, 1, figsize=(15,15))
    g = sns.scatterplot(data=means, x='b [ms/um²]', y='Sb/So', hue='Delta', ax=ax, s=200)
    # ax.set_xticklabels([f'{float(blab):.1f}' for blab in b_labels[1:]])
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles, ['Soma', 'Dendrites', 'Soma-Dendrites (disconnected)', 'Soma-Dendrites (connected)'], loc='upper right', title='Intra signal', markerscale=2.5)

    # replace labels
    # new_labels = [['Soma', 'Dendrites', 'Soma-Dendrites (disconnected)', 'Soma-Dendrites (connected)']]
    # for t, l in zip(g._legend.texts, new_labels):
    #     t.set_text(l)

    ax2 = ax.twinx()
    b_plot = (bvals*1e-9).round(2)
    analytical = analytical_df[analytical_df["case"] == case]
    if case == "soma_dendrites_ex":
        analytical = analytical_df[analytical_df["case"] == "soma_dendrites"]
    if funnel:
        linestyle = "dashed"
    else:
        linestyle = "solid"
    sns.lineplot(data=analytical, x='b [ms/um²]', y='Sb/So', hue='Delta', ax=ax2, linestyle=linestyle)
    
    ax.legend(title='Delta [s]', loc="upper right")
    ax2.legend(title='Analytical solution', loc=3)
    ax2.set_yticklabels([])
    ax2.set_yticks([])
    ax2.set_ylabel("")
    ax2.set_ylim([y_lim_min, y_lim_max])
    ax.set_ylim([y_lim_min, y_lim_max])

    # plt.show()
    plt.savefig(DWI_folder / f"diff_{case}.png")
