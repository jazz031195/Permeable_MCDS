"""
    File that plots the signal
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

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
from matplotlib.ticker import FuncFormatter
import pingouin as pg

def two_sig_fig(x, pos):
    return f"{x:.2f}"

def get_complexity(cfg):
    tort = "tortuous" in cfg
    bead = "beaded" in cfg

    if tort and bead:
        return "undulated-beaded"
    elif tort:
        return "undulated"
    elif bead:
        return "beaded"
    else:
        return "straight"
    
cur_path    = os.getcwd()
scheme_file = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/instructions/scheme/128_dir_12_b_9_td.scheme"
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

df_all_data = pd.read_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/OHBM26/results_geo/data.csv")

volume_neuron = 1357.93 # µm3
params = ["MD", "RD", "FA", "AD", "MK", "RK", "AK"]
params = ["MD", "MK"]
hue_order = [
            "neuron", 
             "neuron_ref", 
             "neuron_tortuous", 
             "neuron_tortuous_ref", 
             "neuron_beaded", 
             "neuron_beaded_ref", 
             "neuron_tortuous_beaded",
             "neuron_tortuous_beaded_ref"
             ]

print(df_all_data.columns)
for SNR in [np.inf, 30]:
    b_labels    = df_all_data["b [ms/um²]"].unique()
    Deltas      = [25.5, 105.5] # df_all_data["Delta"].unique()
    df_all_data = df_all_data[(df_all_data['b [ms/um²]'] > 0)]
    for D in Deltas:
        for p in params:

            fig, ax = plt.subplots(1, 1, figsize=(14, 10))

            df = df_all_data[(df_all_data["Delta [ms]"] == D) & (df_all_data["SNR"] == SNR) & (df_all_data["cells"].str.contains("neuron_0_MCDS"))].copy()

            sns.boxplot(
                data=df,
                x="cells",
                y=p,
                hue="cells",
                ax=ax,
                showfliers=False,
                dodge=True,
                palette=["lightgray", "lightgray"],
            )

            # Points
            sns.stripplot(
                x="cells",
                data=df,
                y=p,
                hue="cells",
                ax=ax,
                dodge=True,
                color="white",
                edgecolor="darkgray",
                linewidth=2,
                zorder=2,
                size=5,
                legend=False,
            )
            if ax.get_legend() is not None:
                ax.get_legend().remove()
                
            plt.show()

            