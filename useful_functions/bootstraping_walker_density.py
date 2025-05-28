import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")
import glob
import re
import nibabel as nib
from useful_functions import read_and_extract_dwi, get_scheme_info
import random
from scipy.stats import sem

def create_group_averages(df, group_size, num_groups, total_reps, walkers_count):
    # Shuffle the indices and ensure they don't repeat across groups
    indices = list(range(total_reps))
    random.shuffle(indices)

    # Ensure we have enough reps to create non-overlapping groups
    if group_size * num_groups > total_reps:
        raise ValueError("Not enough unique reps to create the requested number of non-overlapping groups")

    all_group_averages = []

    for i in range(num_groups):
        group = indices[i * group_size : (i + 1) * group_size]
        df_group = df[df["rep"].isin(group)]
        avg = df_group.groupby("bvals", as_index=False).mean()
        avg["rep"] = i
        avg["nbr_walkers"] = walkers_count
        all_group_averages.append(avg)

    return pd.concat(all_group_averages, ignore_index=True)

def repeated_group_averages(df, group_size, num_groups, total_reps, walkers_count, n_repeats=100):
    """Runs the create_group_averages function multiple times and aggregates results."""
    repeated_results = []

    for _ in range(n_repeats):
        avg = create_group_averages(df.copy(), group_size, num_groups, total_reps, walkers_count)
        avg["norm_DWI"] = avg.groupby(["rep", "nbr_walkers"])["DWI"].transform(lambda x: x / x.max())
        avg = avg.groupby(["nbr_walkers", "bvals"], as_index=False).agg(
            mean_DWI=("norm_DWI", "mean"),
            std_DWI=("norm_DWI", "std")
        )

        repeated_results.append(avg)

    all_runs_df = pd.concat(repeated_results, ignore_index=True)

    return all_runs_df

# === MAIN CODE ===
if __name__ == "__main__":
    path_to_scheme = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/SMI.scheme"
    path_to_folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/bootstrapping_walker_density/"

    bvals, directions = get_scheme_info(path_to_scheme)

    all_files = [f for f in glob.glob(os.path.join(path_to_folder, "*.bfloat")) if "img" not in f]
    all_data = []

    for e, file in enumerate(all_files):
        dwi_array = read_and_extract_dwi(file, True)
        for bval, vec, dwi in zip(bvals, directions, dwi_array):
            all_data.append({
                "bvals": bval,
                "vecs": vec,
                "DWI": dwi,
                "rep": e
            })

    df = pd.DataFrame(all_data)

    nbr_group = 3

    # Group-wise powder averages
    summary_10_4 = repeated_group_averages(df, group_size=1, num_groups=nbr_group, total_reps=len(all_files), walkers_count=1e4)
    summary_2_10_4 = repeated_group_averages(df, group_size=2, num_groups=nbr_group, total_reps=len(all_files), walkers_count=2e4)
    summary_5_10_4 = repeated_group_averages(df, group_size=5, num_groups=nbr_group, total_reps=len(all_files), walkers_count=5e4)
    summary_7_10_4 = repeated_group_averages(df, group_size=7, num_groups=nbr_group, total_reps=len(all_files), walkers_count=7e4)
    summary_10_5 = repeated_group_averages(df, group_size=10, num_groups=nbr_group, total_reps=len(all_files), walkers_count=1e5)
    summary_2_10_5 = repeated_group_averages(df, group_size=20, num_groups=nbr_group, total_reps=len(all_files), walkers_count=2e5)
    summary_5_10_5 = repeated_group_averages(df, group_size=50, num_groups=nbr_group, total_reps=len(all_files), walkers_count=5e5)

    final_df = pd.concat([summary_10_4, summary_2_10_4, summary_5_10_4, summary_7_10_4, summary_10_5, summary_2_10_5, summary_5_10_5], ignore_index=True)
    
    palette = sns.color_palette("husl", len(final_df["bvals"].unique()))
    # Plot Mean
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=final_df, x="nbr_walkers", y="mean_DWI", hue="bvals", marker="o", ci=None, palette=palette)
    plt.xscale("log")
    plt.xlabel("Number of Walkers (log scale)", fontsize=14)
    plt.ylabel(f"Mean powder-averaged DWI", fontsize=14)
    plt.title("Mean Signal per Number of Walkers")
    plt.grid(True)
    # no labels visible
    plt.legend().set_visible(False)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    ax = sns.lineplot(data=final_df, x="nbr_walkers", y="std_DWI", hue="bvals", marker="o", ci=None, palette=palette)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Number of Walkers (log scale)", fontsize=14)
    plt.ylabel(f"Powder-averaged DWI standard deviation (log scale)", fontsize=14)
    plt.title("Signal Standard deviation vs. Number of Walkers")
    plt.grid(True)

    # Customise legend labels
    handles, labels = ax.get_legend_handles_labels()
    new_labels = [f"b-value = {label}" for label in labels]
    ax.legend(handles=handles, labels=new_labels, fontsize=14, title="")

    plt.tight_layout()
    plt.show()
