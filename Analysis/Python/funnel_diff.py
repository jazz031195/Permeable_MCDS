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
import warnings
warnings.filterwarnings("ignore")
import statannot
from utils import get_bvals, get_bvectors, calculate_DKI, get_dwi, get_psge, create_data, analytical_solutions
import json

cur_path    = os.getcwd()
scheme_file = cur_path + "/results/funnel/overlap_4/n1/PGSE_21_dir_12_b.scheme"
giro        = 2.6751525e8 # Gyromagnetic radio [rad/(s*T)]

log  = False

# if log:
#     y_lim_min = -5
#     y_lim_max = 0.1
# else:
#     y_lim_min = 0.
#     y_lim_max = 1

MEDIUM_SIZE = 19
BIGGER_SIZE = 19

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

DWI_folder = Path("results/diff_times_single_neuron")

# Adjacent spheres are distant of R/overlap from each other
overlap = 4
deltas = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07]
# df_all_data, df_crossings = create_df_all(DWI_folder, scheme_file, deltas)

# df_all_data.to_csv(DWI_folder / "data.csv")
df_all = pd.read_csv(DWI_folder / "data.csv")

# Analytical solutions & Mesh
delta     = np.array([0.0165])# in [s]
D0        = 2.5e-9 # [m²/s]

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

analytical_df = pd.DataFrame()
bvals     = np.linspace(1, 10, 10) * 1e9 # in [s/m²]
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

means = df_all[(df_all['b [ms/um²]'] > 0)].groupby(['b [ms/um²]', 'case', 'Delta', 'neuron', 'funnel'])['Sb/So'].mean().reset_index()
means['Sb/So_no_mean'] = means['Sb/So']
for case in means.case.unique():
    if "soma_dendrites" in case:
        case_ = "soma_dendrites"
    else:
        case_ = case
    for b in bvals:
        for td in deltas:
            analytical_solution = analytical_df[(analytical_df['b [ms/um²]'] == b*1e-9) & (analytical_df['Delta'] == td) & (analytical_df['case'] == case_)]['Sb/So'].values
            mask = (means['b [ms/um²]'] == b*1e-9) & (means['Delta'] == td) & (means['case'] == case)
            means.loc[mask, 'Sb/So_no_mean'] = means.loc[mask, 'Sb/So_no_mean'] - analytical_solution



for funnel in [True, False]:
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    sns.violinplot(data=means[((means.case=="soma_dendrites") | (means.case=="soma_dendrites_ex")) & (means['funnel'] == funnel)], 
                x='Delta', 
                y='Sb/So_no_mean',
                hue='case', 
                ax=ax)

    sns.stripplot(data=means[((means.case=="soma_dendrites") | (means.case=="soma_dendrites_ex")) & (means['funnel'] == funnel)], 
                x='Delta', 
                y='Sb/So_no_mean',
                hue='case',
                dodge=True,
                edgecolor='darkgray',
                linewidth=1,
                zorder=1,
                alpha=0.8) 

    ax.legend().set_visible(False)

    for collection in ax.collections:
        if isinstance(collection, matplotlib.collections.PolyCollection):
            collection.set_edgecolor(collection.get_facecolor())
            collection.set_facecolor(collection.get_facecolor())
            collection.set_alpha(0.5)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, ['Soma - Dendrites (no exch)', 'Soma - Dendrites (exch)'], loc='upper left', title='')
    ax.set_ylabel("S/S0 - analytical solution")
    ax.set_ylim([-0.03, 0.04])
    ax.set_xlabel("Delta [s]")
    ax.axhline(y=0, xmin=0, xmax=len(deltas)-1, color='k', linestyle="dotted")
    if funnel:
        ax.set_title("Funnel")
    else:
        ax.set_title("No funnel")

    # couples = []
    # couples_end = []
    # for b in means['funnel'].unique():
    #     for i, branch in enumerate(means['case'].unique()):
    #         couples.append((b, branch))    

    # for i in range(1, len(couples) + 1):
    #     if i % 2 == 0:
    #         couples_end.append((couples[i-2], couples[i-1]))
    # print(couples, couples_end)
    # statannot.add_stat_annotation(
    #     ax,
    #     data=means,
    #     x='Delta', 
    #     y='Sb/So_no_mean', 
    #     hue='funnel', 
    #     box_pairs=couples_end,
    #     test="Mann-Whitney",
    #     text_format="star",
    #     loc="inside"
    #     )


    plt.show()

for funnel in [True, False]:
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    sns.violinplot(data=means[((means.case=="soma") | (means.case=="dendrites")) & (means['funnel'] == funnel)], 
                x='Delta', 
                y='Sb/So_no_mean',
                hue='case', 
                ax=ax)

    sns.stripplot(data=means[((means.case=="soma") | (means.case=="dendrites")) & (means['funnel'] == funnel)], 
                x='Delta', 
                y='Sb/So_no_mean',
                hue='case',
                dodge=True,
                edgecolor='darkgray',
                linewidth=1,
                zorder=1,
                alpha=0.8) 

    ax.legend().set_visible(False)

    for collection in ax.collections:
        if isinstance(collection, matplotlib.collections.PolyCollection):
            collection.set_edgecolor(collection.get_facecolor())
            collection.set_facecolor(collection.get_facecolor())
            collection.set_alpha(0.5)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, ['Dendrites', 'Soma'], loc='upper left', title='')
    ax.set_ylabel("S/S0 - analytical solution")
    # ax.set_ylim([-0.03, 0.04])
    ax.set_xlabel("Delta [s]")
    ax.axhline(y=0, xmin=0, xmax=len(deltas)-1, color='k', linestyle="dotted")
    if funnel:
        ax.set_title("Funnel")
    else:
        ax.set_title("No funnel")

    # couples = []
    # couples_end = []
    # for b in means['funnel'].unique():
    #     for i, branch in enumerate(means['case'].unique()):
    #         couples.append((b, branch))    

    # for i in range(1, len(couples) + 1):
    #     if i % 2 == 0:
    #         couples_end.append((couples[i-2], couples[i-1]))
    # print(couples, couples_end)
    # statannot.add_stat_annotation(
    #     ax,
    #     data=means,
    #     x='Delta', 
    #     y='Sb/So_no_mean', 
    #     hue='funnel', 
    #     box_pairs=couples_end,
    #     test="Mann-Whitney",
    #     text_format="star",
    #     loc="inside"
    #     )


    plt.show()