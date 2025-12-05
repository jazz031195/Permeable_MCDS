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

df_all_data = pd.read_csv("/Users/ideriedm/Documents/MCDS_perm/Permeable_MCDS/results/OHBM26/data.csv")
df_all_data = df_all_data[df_all_data["SNR"] == np.inf]
df_all_data["Delta"] = df_all_data["Delta [ms]"]
b_labels    = df_all_data["b [ms/um²]"].unique()
Deltas      = df_all_data["Delta"].unique()
df_all_data = df_all_data[(df_all_data['b [ms/um²]'] > 0)]
volume_neuron = 1357.93 # µm3
params = ["MD", "RD", "FA", "AD", "MK", "RK", "AK"]
hue_order = ["neuron", "neuron_tortuous", "neuron_beaded", "neuron_tortuous_beaded"]

for D in Deltas:
    for p in params:
        fig, ax = plt.subplots(1, 1, figsize=(6,5))

        df     = df_all_data[(df_all_data['Delta'] == D) & ~(df_all_data["configuration"].str.contains("ref"))].copy()
        df_ref = df_all_data[(df_all_data['Delta'] == D) & (df_all_data["configuration"].str.contains("ref"))].copy()
        
        means_ref = df_ref.groupby(["configuration_all"])[p].mean()
        d = {p: means_ref.values,
             "conf": [v[:-2] for v in means_ref.index.values]}
        to_plot_ref = pd.DataFrame(d)

        means     = df.groupby(["configuration_all"])[p].mean()
        d = {p: means.values,
             "conf": [v[:-2] for v in means.index.values]}
        to_plot = pd.DataFrame(d)


        # Plot
        g = sns.boxplot(
            data=to_plot,
            x='conf',
            y=p,
            hue='conf',
            hue_order=hue_order,
            order=hue_order,
            ax=ax,
            showfliers=False
        )

        sns.stripplot(
            data=to_plot,
            x='conf',
            y=p,
            order=hue_order,
            ax=ax,
            color="white",
            edgecolor='darkgray',
            linewidth=2,
            size=5
        )

        sns.stripplot(
            data=to_plot_ref,
            x='conf',
            y=p,
            hue='conf',
            ax=ax,
            marker='D',
            edgecolor='darkgray',
            linewidth=2,
            color="k",
            size=10
        )

        print(df_ref.configuration, df_ref.columns, hue_order)
        sns.boxplot(
            data=df_ref,
            x='configuration',
            y=p,
            hue='configuration',
            order=[h + "_ref" for h in hue_order],
            ax=ax
        )

        ax.set_xticks([])
        ax.set_xlabel("")
        ax.set_xticklabels([])
    
        if not os.path.exists(f"results/OHBM26/{p}/"):
            os.makedirs(f"results/OHBM26/{p}/")
        plt.tight_layout()
        fig.savefig(f"results/OHBM26/{p}/{p}_{D}.png")
        fig.savefig(f"results/OHBM26/{p}/{p}_{D}.pdf")