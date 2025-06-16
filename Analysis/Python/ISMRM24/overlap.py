""" 
    File that compares MD and normalized signal for different 
    overlap between adjacent spheres
"""

import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

if sys.platform == "linux":
    sys.path.insert(1, '/home/localadmin/Documents/analytical_formula/')
else:
    sys.path.insert(1, '/Users/ideriedm/Documents/analytical_formula/')

from my_murdaycotts import my_murdaycotts
import statannot
from utils import get_bvals, get_bvectors, calculate_DKI, get_dwi, get_psge, create_data

cur_path    = os.getcwd()
scheme_file = cur_path + "/results/funnel/overlap_4/n1/PGSE_21_dir_12_b.scheme"
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
    
    df_all_data  = pd.DataFrame()
    df_crossings = pd.DataFrame()
    for overlap in os.listdir(DWI_folder):
        if os.path.isdir(DWI_folder / overlap) and "overlap" in overlap:
            for neuron in os.listdir(DWI_folder / overlap):
                if os.path.isdir(DWI_folder / overlap / neuron):
                    # Iterate through the files in the folder
                    for subdir in os.listdir(DWI_folder / overlap / neuron):
                        if os.path.isdir(DWI_folder / overlap / neuron / subdir):
                            for filename in os.listdir(DWI_folder / overlap / neuron / subdir):
                            
                                # Read the simulation_info.txt to have crossings information
                                if "simu" in filename:
                                    N = int(subdir.split('_')[1])
                                    T = int(subdir.split('_')[3])
                                    with open(DWI_folder / overlap / neuron / subdir / filename, 'r') as file:
                                        # Read the file line by line
                                        for line in file:
                                            if 'Particle dynamics duration' in line:
                                                D0 = float(line.split()[-2])
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
                                    SNR          = np.inf
                                    data_one_exp = create_data(DWI_folder / overlap / neuron /subdir, SNR, name, extension, scheme_file_path)
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
                                        if i == 1:
                                            d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc, 
                                                'b [ms/um²]': bval, 'neuron': neuron, 'overlap': int(overlap[7:]),
                                                'FA': FA, 'MD': MD, 'RD': RD, 'AD': AD, 'MK': MK, 'RK': RK, 'AK': AK}
                                        else:
                                            d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc, 
                                                'b [ms/um²]': bval, 'neuron': neuron, 'overlap': int(overlap[7:])}
                                        df_avg_data = pd.DataFrame(d, index=[i])
                                        df_all_data = pd.concat([df_all_data, df_avg_data])

    return df_all_data, df_crossings



y_lim_min = 0.
y_lim_max = 1.1

MEDIUM_SIZE = 19
BIGGER_SIZE = 19

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)   # fontsize of the figure title


experience_folder = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/ISMRM24/overlaps/")
# df_data_all, df_crossings = create_df_all(experience_folder, scheme_file)

# df_data_all.to_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/ISMRM24/overlaps/data.csv")
df_data_all = pd.read_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/ISMRM24/overlaps/data.csv")

b_labels = df_data_all["b [ms/um²]"].unique()

fig, _ = plt.subplots(1, 1, figsize=(10, 7))
ax2 = sns.violinplot(data=df_data_all[(df_data_all['b [ms/um²]'] > 4)], 
               x='b [ms/um²]', 
               y='Sb/So',
               hue='overlap', 
               dodge=True,
               palette="Blues")

sns.stripplot(data=df_data_all[(df_data_all['b [ms/um²]'] > 4)], 
              x='b [ms/um²]', 
              y='Sb/So',
              hue='overlap',
              palette=[(0, 0, 0)],
              dodge=True,
              size=3) 

ax2.legend().set_visible(False)

for collection in ax2.collections:
    if isinstance(collection, matplotlib.collections.PolyCollection):
        collection.set_edgecolor(collection.get_facecolor())
        collection.set_facecolor(collection.get_facecolor())
        collection.set_alpha(0.5)

# Change b-values so that they have only 1 decimals
ax2.set_xticklabels([f'{float(blab):.1f}' for blab in b_labels[6:]])
handles, labels = ax2.get_legend_handles_labels()
labels          = ['R/' + item for item in labels]
plt.legend(handles[:int(len(labels)/2)], labels[:int(len(labels)/2)], loc='upper right', frameon=False, title="Overlap")

# couples = []
# couples_end = []
# for b in df_dwi[df_dwi['b [ms/um²]'] > 0]['b [ms/um²]'].unique():
#     for i, branch in enumerate(df_dwi['overlap'].unique()):
#         couples.append((b, branch))    

# for i in range(1, len(couples) + 1):
#     if i % 6 == 0:
#         couples_end.append((couples[i-6], couples[i-5]))
#         couples_end.append((couples[i-6], couples[i-4]))
#         couples_end.append((couples[i-6], couples[i-3]))
#         couples_end.append((couples[i-6], couples[i-2]))
#         couples_end.append((couples[i-6], couples[i-1]))
#         couples_end.append((couples[i-5], couples[i-4]))
#         couples_end.append((couples[i-5], couples[i-3]))
#         couples_end.append((couples[i-5], couples[i-2]))
#         couples_end.append((couples[i-5], couples[i-1]))
#         couples_end.append((couples[i-4], couples[i-3]))
#         couples_end.append((couples[i-4], couples[i-2]))
#         couples_end.append((couples[i-4], couples[i-1]))
#         couples_end.append((couples[i-3], couples[i-2]))
#         couples_end.append((couples[i-3], couples[i-1]))
#         couples_end.append((couples[i-2], couples[i-1]))

# statannot.add_stat_annotation(
#     ax,
#     data=df_dwi[(df_dwi['b [ms/um²]'] > 0) & (df_dwi["T"] == 15000)],
#     y='Sb/So', x='b [ms/um²]',
#     hue='overlap',
#     box_pairs=couples_end,
#     test="Mann-Whitney",
#     text_format="star",
#     loc="inside"
#     )
# ax.set_title(f'N = {N_labels[0]}, T = {T_labels[0]}')


fig, _ = plt.subplots(1,1, figsize=(10, 7))
ax = sns.boxplot(data=df_data_all[(df_data_all['b [ms/um²]'] > 0)], 
                 x='overlap', 
                 y='MD', 
                 palette="Blues")
sns.stripplot(data=df_data_all[(df_data_all['b [ms/um²]'] > 0)], 
              x='overlap', 
              y='MD', 
              color=(0.1350634371395617, 0.35630911188004605, 0.5703575547866206, 1))

from matplotlib.patches import PathPatch

box_patches = [p for p in ax.patches if isinstance(p, PathPatch)]
for patch in box_patches:
    patch.set_edgecolor(patch.get_facecolor())
    patch.set_facecolor(patch.get_facecolor())
    patch.set_alpha(0.5)
labels = ['R/' + item.get_text() for item in ax.get_xticklabels()]
ax.set_xticklabels(labels)
ax.set_xlabel('Distance between overlapping spheres')
ax.set_ylabel('Mean diffusivity')
plt.show()