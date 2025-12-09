"""
    File that plots the signal
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
if sys.platform == "linux":
    sys.path.insert(1, '/home/localadmin/Documents/analytical_formula/')
    sys.path.insert(1, '/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS')
import pandas as pd
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
from utils import get_bvals, get_bvectors, calculate_DKI, get_dwi, get_psge, analytical_solutions

from fig_style import fig_size

cur_path    = os.getcwd()
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
    
    df_all_data  = pd.DataFrame()
    for filename in os.listdir(DWI_folder):
        # Check if the filename contains "_rep_" and "DWI"
        if "DWI_img" in filename:
            simu_info = DWI_folder / filename.replace("DWI_img.bfloat", "simulation_info.txt").replace("DWI_img.txt", "simulation_info.txt")
            with open(simu_info, 'r') as file:
                # Read the file line by line
                for line in file:
                    if ' Number of particles:' in line:
                        N = int(line.split('-')[-1][:-1])
                    if 'Number of steps:' in line:
                        T = int(line.split('-')[-1])

            # Name of the experience
            name         = ('_').join(filename.split('_')[:-1])
            extension    = filename.split('_')[-1].split('.')[-1]
            extension    = "txt"
            SNR          = np.inf
            data_one_exp = create_data(DWI_folder , SNR, name, extension, scheme_file_path)
            # For each b, iterate over all directions, store the data, and average them (powder-average)
            nb_b   = len(data_one_exp["b [ms/um²]"].round(1).unique())
            data_psge = get_psge(scheme_file_path)
            Deltas = data_psge['Delta [ms]'].unique()
            nb_td  = len(Deltas) 
            nb_dir = int(len(data_one_exp["x"].values) / (nb_b * nb_td))
            
            for i in range(nb_b):
                for k in range(nb_td):
                    sb_so = []
                    adc   = []
                    for j in range(nb_dir):
                        sb_so.append(data_one_exp.iloc[nb_b * j + i + k*nb_dir*nb_b, :]["Sb/So"])
                        adc.append(data_one_exp.iloc[nb_b * j + i + k*nb_dir*nb_b, :]["adc [ms/um²]"])
                        bval = data_one_exp.iloc[nb_b * j + i + k*nb_dir*nb_b, :]["b [ms/um²]"].round(1)
                    
                    # Powder-average signal
                    mean     = np.mean(sb_so)
                    # Powder-average ADC
                    mean_adc = np.mean(adc)

                    d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, 
                        'b [ms/um²]': bval, 'Delta [ms]': Deltas[k]}
                    df_avg_data = pd.DataFrame(d, index=[i])
                    df_all_data = pd.concat([df_all_data, df_avg_data])

    return df_all_data


# branching = "branching"

log  = False

if log:
    y_lim_min = -5
    y_lim_max = 0.1
else:
    y_lim_min = -0.05
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

DWI_folder  = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/test_juliette")
scheme_file = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/instructions/scheme/128_dir_12_b_5_td.scheme"
# df_all_data = create_df_all(DWI_folder, scheme_file)

# df_all_data.to_csv(DWI_folder / "data.csv")
df_all_data = pd.read_csv(DWI_folder / "data.csv")

data_psge = get_psge(scheme_file)
b_labels  = data_psge["b [ms/um²]"].unique()
Deltas    = data_psge['Delta [ms]'].unique()
delta     = data_psge['delta [ms]'].unique()
# 2.0e-9 [m²/s] -> 2.0 [um²/ms]
D0        = 2.0 # [um²/ms]
df_tmp = df_all_data

# means = df_tmp.groupby(['b [ms/um²]', 'Delta [ms]'])['Sb/So'].mean().reset_index()
# means['Delta [ms]'] = means['Delta [ms]'].astype('category')
for Delta in df_tmp['Delta [ms]'].unique():
    df_tmp = df_tmp[df_tmp['Delta [ms]'] == Delta]
    print(df_tmp)
    fig, ax = plt.subplots(1, 1, figsize=fig_size(fraction=0.42, height_ratio=0.8))
    g = sns.boxplot(data=df_tmp, x='b [ms/um²]', y='Sb/So', ax=ax)
    ax.legend(title='Delta [ms]', loc="upper right")
    ax.set_ylim([y_lim_min, y_lim_max])

    plt.show()
    plt.savefig(DWI_folder / f"diff_{delta}.png")
    plt.savefig(DWI_folder / f"diff_{delta}.pdf")