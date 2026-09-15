"""
    File to compare the funnel vs no funnel data
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
import scipy
import matplotlib
import warnings
warnings.filterwarnings("ignore")
import statannot
from utils import get_bvals, get_bvectors, calculate_DKI, get_dwi, get_psge, create_data, analytical_solutions
import json

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
                                        N = int(subdir.split('_')[1])
                                        T = int(subdir.split('_')[3])
                                        with open(DWI_folder / overlap / subcase / neuron / subdir / filename, 'r') as file:
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
                                                d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, 'SNR': [str(SNR)], 
                                                    'b [ms/um²]': bval, 'neuron': neuron, 'funnel': subcase}
                                                df_avg_data = pd.DataFrame(d, index=[i])
                                                df_all_data = pd.concat([df_all_data, df_avg_data])

    return df_all_data, df_crossings

log  = False

if log:
    y_lim_min = -5
    y_lim_max = 0.1
else:
    y_lim_min = 0.
    y_lim_max = 1

MEDIUM_SIZE = 19
BIGGER_SIZE = 19

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

DWI_folder = Path("results/ISMRM24/funnel")

# Adjacent spheres are distant of R/overlap from each other
overlap = 4
# df_dwi, df_crossings_no_funnel = create_df_all(DWI_folder, scheme_file)
# df_dwi.to_csv(DWI_folder / "data_SNR.csv")
df_dwi = pd.read_csv(DWI_folder / "data_SNR.csv")


df_all = df_dwi
df_all = df_all[df_all['b [ms/um²]'] > 0]

fig, ax = plt.subplots(1, 1, figsize=(8.5, 6))
sns.violinplot(data=df_all[df_all['SNR'] == np.inf], 
               x='b [ms/um²]', 
               y='Sb/So', 
               hue='funnel', 
               ax=ax)

for collection in ax.collections:
    collection.set_edgecolor(collection.get_facecolor())
ax.set_xticks([0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, ['No tapering', 'Tapering'], loc='upper right', title='')
if log:
    ax.set_ylabel("ln(S/S0)")
else:
    ax.set_ylabel(r"$S_b/S_0$")

# Analytical solutions & Mesh
Delta     = np.array([0.05])  # in [s]
delta     = np.array([0.0165])# in [s]
D0        = 2.5e-9 # [m²/s]
bvals     = np.linspace(0.2, 10, 100) * 1e9 # in [s/m²]

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

soma_signal, neurites_signal, both_signal = analytical_solutions(bvals, Delta, delta, r_soma, D0, log, soma_fraction, neurite_fraction)

ax2 = ax.twinx()
ax2.plot(bvals*1e-9, both_signal, label=f"Soma & dendrites", color='k', linestyle="dotted")
if log:
    sns.lineplot(bvals*1e-9, -bvals*D0, label="D = 2.5 [ms/um²]")

ax2.set_yticklabels([])
ax2.set_yticks([])
ax2.set_ylim([y_lim_min, y_lim_max])
ax.set_ylim([y_lim_min, y_lim_max])

plt.show()











bvals = df_all['b [ms/um²]'].values * 1e9
soma_signal, neurites_signal, both_signal = analytical_solutions(bvals, Delta, delta, r_soma, D0, log, soma_fraction, neurite_fraction)
# Subtract the analytical value of the corresponding b-value to the signal ("de-meaning")
df_all['Sb/So_no_mean'] = df_all['Sb/So'].values -  both_signal
df_all['SNR'] = df_all['SNR'].astype(str)
# Plotting
fig, ax = plt.subplots(1, 1, figsize=(12, 9))
sns.violinplot(data=df_all, 
               x='SNR', 
               y='Sb/So_no_mean', 
               hue='funnel', 
               ax=ax, 
               order=['inf', '100.0', '60.0', '20.0'])

sns.stripplot(data=df_all, 
              x='SNR', 
              y='Sb/So_no_mean', 
              hue='funnel', 
              order=['inf', '100.0', '60.0', '20.0'],
              dodge=True) 

ax.legend().set_visible(False)

for collection in ax.collections:
    if isinstance(collection, matplotlib.collections.PolyCollection):
        collection.set_edgecolor(collection.get_facecolor())
        collection.set_facecolor(collection.get_facecolor())
        collection.set_alpha(0.5)

plt.ylabel(r"$S_b/S_0$ - analytical signal (soma + dendrites)")
plt.axhline(y=0, linestyle="dashed", color='gray')
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[:int(len(labels)/2)], ['No tapering', 'Tapering'], loc='upper left', frameon=False)

couples = []
couples_end = []
for b in df_all['SNR'].unique():
    for i, branch in enumerate(df_all['funnel'].unique()):
        couples.append((b, branch))    
        print(b, branch, "H0 : population average is 0 (the simulated signal = the analytical one)")
        # scipy.stats.wilcoxon
        print(scipy.stats.ttest_1samp(df_all[(df_all["SNR"] == b) & (df_all['funnel'] == branch)]['Sb/So_no_mean'].values, popmean=0., nan_policy='omit'))
        print("\n")

for i in range(1, len(couples) + 1):
    if i % 2 == 0:
        couples_end.append((couples[i-2], couples[i-1]))

statannot.add_stat_annotation(
    ax,
    data=df_all,
    y='Sb/So_no_mean', 
    x='SNR',
    hue='funnel',
    hue_order=['nofunnel', 'funnel'],
    box_pairs=couples_end,
    test="t-test_ind",
    text_format="star",
    loc="inside"
    )

plt.show()