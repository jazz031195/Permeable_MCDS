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
    Deltas      = [25.5, 105.5] # df_all_data["Delta"].unique()
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

            sns.boxplot(
                data=to_plot,
                x="complexity",
                y=p,
                hue="Delta [ms]",
                order=order_complexity,
                hue_order=Deltas,
                ax=ax,
                showfliers=False,
                dodge=True,
                palette=["lightgray", "lightgray"],
            )

            # Points
            sns.stripplot(
                data=to_plot,
                x="complexity",
                y=p,
                hue="Delta [ms]",
                order=order_complexity,
                hue_order=Deltas,
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

            boxes = [p for p in ax.patches if isinstance(p, mpl.patches.PathPatch)]

            boxes_sorted = sorted(
                boxes,
                key=lambda b: np.mean(b.get_path().vertices[:, 0])
            )

            n_comp   = len(order_complexity)
            n_struct = len(Deltas)

            for g_idx in range(n_comp):
                comp = order_complexity[g_idx]
                base_color = mpl.colors.to_rgba(colors_comp[comp])

                for j in range(n_struct): 
                    patch = boxes_sorted[g_idx * n_struct + j]
                    alpha = 1.0 if Deltas[j] == 25.5 else 0.25
                    patch.set_facecolor((*base_color[:3], alpha))
                    patch.set_edgecolor(base_color[:3])


            
            if p == "MD":
                ax.set_ylim(0.25, 0.45)
            elif p== "MK":
                ax.set_ylim(1.4, 2)

            y_min, y_max = ax.get_ylim()
            y_text = y_max + 0.01 * (y_max - y_min)

            for c, x_center in enumerate(ax.get_xticks()):
                ax.text(x_center - 0.24, y_text, "Δ$_1$",
                        ha="center", va="bottom", fontsize=MEDIUM_SIZE - 2, color=colors[c])
                ax.text(x_center + 0.24, y_text, "Δ$_2$",
                        ha="center", va="bottom", fontsize=MEDIUM_SIZE - 2, color=colors[c])
            ax.set_xlabel("")
            ax.set_xticklabels(xtick_labels, rotation=60)

            ax.set_ylabel(p)

            plt.tight_layout()
            if not os.path.exists(f"results/OHBM26/SNR_{SNR}/{p}/"):
                os.makedirs(f"results/OHBM26/SNR_{SNR}/{p}/")
            plt.grid()
            plt.tight_layout()
            print(f"results/OHBM26/SNR_{SNR}/{p}/{p}_{b}.png")
            fig.savefig(f"results/OHBM26/SNR_{SNR}/{p}/{p}_{b}.png")
            fig.savefig(f"results/OHBM26/SNR_{SNR}/{p}/{p}_{b}.pdf")