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
from sklearn.linear_model import LinearRegression
from scipy.stats import linregress

cur_path    = os.getcwd()
scheme_file = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/instructions/scheme/128_dir_12_b_9_td.scheme"
giro        = 2.6751525e8 # Gyromagnetic radio [rad/(s*T)]

def get_info(simu_info):
    with open(simu_info, 'r') as file:
        # Read the file line by line
        for line in file:
            if 'Number of particles:' in line:
                N = int(line.split('-')[-1])
            if 'Number of steps:' in line:
                T = int(line.split('-')[-1])
            
        d = {'N': [N], 'T': [T]}
    return pd.DataFrame(d)

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
    # timestep 
    T = 15500
    
    df_all_data  = pd.DataFrame()
    df_crossings = pd.DataFrame()
    # for neuron in os.listdir(DWI_folder):
    #     if os.path.isdir(DWI_folder / neuron):
    #         # Iterate through the files in the folder
    for subdir in os.listdir(DWI_folder):
        if os.path.isdir(DWI_folder / subdir) and subdir != "output" and subdir != "old":
            for filename in os.listdir(DWI_folder / subdir):

                # Check if the filename contains "_rep_" and "DWI"
                if "DWI_img" in filename:
                    simu_info = str(DWI_folder / subdir / filename).replace("DWI_img", "simulation_info")
                    info = get_info(simu_info)

                    # Name of the experience
                    name         = ('_').join(filename.split('_')[:-1])
                    extension    = filename.split('_')[-1].split('.')[-1]
                    SNR          = np.inf
                    for SNR in [np.inf]:
                        print(filename)
                        data_one_exp = create_data(DWI_folder / subdir, SNR, name, extension, scheme_file)

                        data = data_one_exp.copy()
                        data["b_shell"] = data["b [ms/um²]"].round(1)
                        grouped = (
                            data
                            .groupby(["Delta [ms]", "b_shell"], as_index=False)
                            .agg({
                                "Sb/So": ["mean", "std"],
                                "adc [ms/um²]": "mean",
                                "MD": "mean",   
                                "AD": "mean",   
                                "FA": "mean",   
                                "RD": "mean",   
                                "MK": "mean",   
                                "RK": "mean",   
                                "AK": "mean",   
                            }))
                        
                        # print(grouped)
                        # print(grouped["Sb/So"])
                        # print(grouped["b_shell"])
                        # plt.figure(figsize=(12,8))

                        # x = grouped[grouped["Delta [ms]"] == 25.5]["b_shell"].values
                        # y = grouped[grouped["Delta [ms]"] == 25.5][("Sb/So", "mean")].values
                        # y_std = grouped[grouped["Delta [ms]"] == 25.5][("Sb/So", "std")].values
                        
                        # # line plot for the mean
                        # plt.plot(x, y)

                        # # shaded area for std
                        # plt.fill_between(x, y - y_std, y + y_std, alpha=0.3)

                        # plt.xlabel("x")
                        # plt.ylabel("Sb/So")
                        # plt.legend()
                        # plt.tight_layout()
                        # plt.savefig("signal.png")
                        # plt.savefig("signal.pdf")
                        # assert(0)

                        conf = "neuron"
                        if "tortuous_beaded" in subdir:
                            conf += "_tortuous_beaded"
                        elif "tortuous" in subdir:
                            conf += "_tortuous"
                        elif "beaded" in subdir:
                            conf += "_beaded"

                        if "_0_" in subdir:
                            conf += "_ref"
                            

                        grouped["configuration"] = conf
                        grouped["cells"] = subdir
                        grouped["SNR"] = SNR
                        grouped["N"] = info.N
                        grouped["T"] = info.T

                        grouped = grouped.rename(columns={
                            "b_shell": "b [ms/um²]",
                        })
                        df_avg_data = grouped[["N", "T", "Sb/So", "b [ms/um²]", "Delta [ms]", "configuration", "cells", "MD", "RD", "AD", "FA", "MK", "RK", "AK", "SNR"]]
                        df_all_data = pd.concat([df_all_data, df_avg_data], ignore_index=True)
    return df_all_data, df_crossings

MEDIUM_SIZE = 35
BIGGER_SIZE = 35
plt.rc('font',   size=MEDIUM_SIZE)       # controls default text sizes
plt.rc('axes',   titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes',   labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick',  labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick',  labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

DWI_folder = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/OHBM26/results_geo/")
# df_all_data, df_crossings = create_df_all(DWI_folder, scheme_file)

# df_all_data.to_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/OHBM26/results_geo/data.csv")
df_all_data = pd.read_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/OHBM26/results_geo/data.csv")
df_all_data = df_all_data[df_all_data.SNR == np.inf]

def PlotSignal(df_all_data) :

    df_all_data["Delta"] = df_all_data["Delta [ms]"]
    b_labels    = df_all_data["b [ms/um²]"].unique()
    Deltas      = df_all_data["Delta"].unique()
    df_all_data = df_all_data[(df_all_data['b [ms/um²]'] > 0)]
    volume_neuron = 1357.93 # µm3

    # Analytical solutions
    delta     = np.array([0.0165])# in [s]
    D0        = 2.0e-9 # [m²/s]
    bvals     = np.linspace(0.2, 10, 100) * 1e9 # in [s/m²]

    r_soma           = 4e-6 # [m]
    volume_soma      = 4/3 * np.pi * r_soma**3 # in [m³]
    volume_soma      = volume_soma * 1e18 # in [um³]
    volume_neuron    = 1357.93 
    volume_neurites  = volume_neuron - volume_soma
    neurite_fraction = volume_neurites / volume_neuron
    soma_fraction    = volume_soma / volume_neuron
    for D in Deltas:

        soma_signal, neurites_signal, both_signal = analytical_solutions(bvals, D/1000, delta, r_soma, D0, False, soma_fraction, neurite_fraction)

        fig, ax = plt.subplots(1, 1, figsize=(14,10))
        fig.subplots_adjust(top=0.7)
        df = df_all_data[(df_all_data['Delta'] == D) & ~(df_all_data["configuration"].str.contains("ref")) & (df_all_data["SNR"] == np.inf)].copy()
        
        # Compute mean and std for each (b, configuration)
        stats = (
            df.groupby(['b [ms/um²]', 'configuration'])
            .agg(mean_S=('Sb/So', 'mean'),
                std_S=('Sb/So', 'std'))
            .reset_index()
        )

        hue_order = [
                     "neuron", 
                     "neuron_beaded", 
                     "neuron_tortuous", 
                     "neuron_tortuous_beaded",
                     ]
        colors = [
            "tab:blue", 
            "tab:orange", 
            "tab:green", 
            "tab:red", 
             ]

        sns.scatterplot(data=stats,
                    x="b [ms/um²]",
                    y="mean_S",
                    hue='configuration',
                    hue_order=hue_order,
                    ax=ax,
                    marker = 'x',
                    linewidths =105,
                    s=500,
                    )
        

        ax.set_ylabel("S$_b$/S$_0$")
        ax.set_ylim([0, 1])
        ax.set_xlabel("b [ms µm⁻²]")
        # ax2 = ax.twinx()
        ax.plot(bvals*1e-9, both_signal, label=f"Sphere & sticks", color='k', linestyle="solid")
        # ax.plot(bvals*1e-9, neurites_signal, label=f"Sticks", color='k', linestyle="dashed")
        # ax.plot(bvals*1e-9, soma_signal, label=f"Sphere", color='k', linestyle="dotted")
        handles, labels = ax.get_legend_handles_labels()

        labels = [label.replace("neuron_", "").replace("_", " ") for label in labels]
        labels[0] = "straight"
        # ax.legend(handles, labels, ncol=2, bbox_to_anchor=(1, 1.5)).set_title('')
        plt.grid()

        fig.tight_layout()

        if not os.path.exists("results/OHBM26/signal_decay/"):
            os.makedirs("results/OHBM26/signal_decay/")
        fig.savefig(f"results/OHBM26/signal_decay/signal_decay_{D}.png")
        print(f"results/OHBM26/signal_decay/signal_decay_{D}.pdf")
        fig.savefig(f"results/OHBM26/signal_decay/signal_decay_{D}.pdf")

def PlotSignalDiff(df_all_data):

    df_all_data["Delta"] = df_all_data["Delta [ms]"]
    b_labels    = df_all_data["b [ms/um²]"].unique()
    Deltas      = df_all_data["Delta"].unique()
    df_all_data = df_all_data[(df_all_data['b [ms/um²]'] > 0)]
    volume_neuron = 1357.93 # µm3
    for D in Deltas:
        fig, ax = plt.subplots(1, 1, figsize=(14,10))
        df = df_all_data[(df_all_data['Delta'] == D) & ~(df_all_data["configuration"].str.contains("ref"))].copy()
        df_ref = df_all_data[(df_all_data['Delta'] == D) & (df_all_data["configuration"].str.contains("ref"))].copy()

        # Compute mean and std for each (b, configuration)
        stats = (
            df.groupby(['b [ms/um²]', 'configuration'])
            .agg(mean_S=('Sb/So', 'mean'),
                std_S=('Sb/So', 'std'))
            .reset_index()
        )

        stats_ref = (
            df_ref.groupby(['b [ms/um²]', 'configuration'])
            .agg(mean_S=('Sb/So', 'mean'),
                std_S=('Sb/So', 'std'))
            .reset_index()
        )

        configs = ['neuron', 'neuron_beaded', 'neuron_tortuous', 'neuron_tortuous_beaded']

        palette = sns.color_palette(n_colors=len(configs))
        palette_map = dict(zip(configs, palette))

        stats_ref['config_base'] = stats_ref['configuration'].str.replace('_ref', '', regex=False)
        stats['config_base'] = stats['configuration']

        merged = stats.merge(
            stats_ref[['b [ms/um²]', 'config_base', 'mean_S', 'std_S']],
            on=['b [ms/um²]', 'config_base'],
            suffixes=('', '_ref')
        )

        merged['diff'] = merged['mean_S'] - merged['mean_S_ref']
        merged['diff_std'] = np.sqrt(merged['std_S']**2 + merged['std_S_ref']**2)

        # x = merged[merged["b [ms/um²]"] > 0.2]["b [ms/um²]"].values
        # y = merged[merged["b [ms/um²]"] > 0.2]['diff'].values
 
        # result = linregress(x, y)
        # x = merged["b [ms/um²]"].values
        # y_fit = result.slope * x + result.intercept
        # print("Slope:", result.slope)
        # print("Intercept:", result.intercept)
        # print("R-squared:", result.rvalue**2)
        # print("p-value:", result.pvalue)
        # print("Std. error:", result.stderr)


        sns.scatterplot(
            data=merged,
            x='b [ms/um²]',
            y='diff',
            hue='config_base',
            marker = 'x',
            linewidths =105,
            s=500,
            hue_order=configs,
            palette = palette_map,
            ax=ax
        )

        # ax.plot(x, y_fit, color='k', label=f"{result.slope:.2f}x + {result.intercept:.2f}\nR2: {result.rvalue**2:.2f}, p: {result.pvalue:.2e}")
        ax.axhline(0, color='black', linestyle='--', linewidth=1)

        ax.set_xlabel("b [ms µm⁻²]")
        ax.set_ylabel("Δ (Branch. - Non-branch.)")
        # ax.set_title(f"Delta = {D} ms")
        ax.set_ylim([-0.011, 0.075])
        ax.grid(True)
        handles, labels = ax.get_legend_handles_labels()
        labels = [label.replace("neuron_", "").replace("_", " ").replace("tortuous", "ondulated") for label in labels]
        labels[0] = "straight"

        leg = ax.legend(
            handles,
            labels,
            ncol=2,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.05),
            frameon=True,
        )
        leg.set_title('')

        fig.tight_layout()

        if not os.path.exists("results/OHBM26/signal_decay_diff/"):
            os.makedirs("results/OHBM26/signal_decay_diff/")
        fig.savefig(f"results/OHBM26/signal_decay_diff/signal_decay_{D}.png")
        fig.savefig(f"results/OHBM26/signal_decay_diff/signal_decay_{D}.pdf")

PlotSignal(df_all_data)
PlotSignalDiff(df_all_data)