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


MEDIUM_SIZE = 40
BIGGER_SIZE = 35
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
             "neuron", 
             "neuron_tortuous", 
             "neuron_tortuous", 
             "neuron_beaded", 
             "neuron_beaded", 
             "neuron_tortuous_beaded",
             "neuron_tortuous_beaded"
             ]

colors = [
            "tab:blue", 
            "tab:orange", 
            "tab:green", 
            "tab:red", 
             ]



print(df_all_data.columns)
for SNR in [np.inf, 30]:
    b_labels    = df_all_data["b [ms/um²]"].unique()
    Deltas      = df_all_data["Delta [ms]"].unique()
    df_all_data = df_all_data[(df_all_data['b [ms/um²]'] > 0)]
    for b in ["branched", "non-branched"]:
        for p in params:
            order_complexity = [
                "straight",
                "beaded",
                "undulated",
                "undulated-beaded",
            ]
            xtick_labels     = ["Straight", "Beaded", "Undulated", "Und.-Bead."]

            colors_comp = {
                "straight": "tab:blue",
                "beaded": "tab:orange",
                "undulated": "tab:green",
                "undulated-beaded": "tab:red",
            }


            fig, ax = plt.subplots(1, 1, figsize=(8, 10))

            df = df_all_data[(df_all_data["Delta [ms]"].isin(Deltas)) & (df_all_data["SNR"] == SNR)].copy()

            if b == "non-branched":
                df = df[df.configuration.str.contains("_ref")]
            else:
                df = df[~df.configuration.str.contains("_ref")]

            means = df.groupby(["cells", "configuration", "Delta [ms]"])[p].mean()

            d = {
                p: means.values,
                "configuration": means.index.get_level_values("configuration"),
                "Delta [ms]": means.index.get_level_values("Delta [ms]")
            }
            to_plot = pd.DataFrame(d)

            to_plot["complexity"] = to_plot["configuration"].apply(get_complexity)

            sns.lineplot(
                data=to_plot,
                x="Delta [ms]",
                y=p,
                hue="complexity",
                hue_order=order_complexity,
                ax=ax,
                markers=True,
                style="complexity"
            )


            if ax.get_legend() is not None:
                ax.get_legend().remove()

            
            if p == "MD":
                ax.set_ylim(0.25, 0.45)
            elif p== "MK":
                ax.set_ylim(1.4, 2)

            y_min, y_max = ax.get_ylim()
            y_text = y_max + 0.01 * (y_max - y_min)

            ax.set_xlabel("Δ [ms]")
            # ax.set_xticklabels(xtick_labels, rotation=60)

            ax.set_ylabel(p)

            plt.tight_layout()
            if not os.path.exists(f"results/OHBM26/SNR_{SNR}/{p}/"):
                os.makedirs(f"results/OHBM26/SNR_{SNR}/{p}/")
            plt.grid()
            plt.tight_layout()
            print(f"results/OHBM26/SNR_{SNR}/{p}/{p}_{b}.png")
            
            fig.savefig(f"{p}_{b}_{SNR}.png")
            fig.savefig(f"{p}_{b}_{SNR}.pdf")
            fig.savefig(f"results/OHBM26/SNR_{SNR}/{p}/{p}_{b}.png")
            fig.savefig(f"results/OHBM26/SNR_{SNR}/{p}/{p}_{b}.pdf")