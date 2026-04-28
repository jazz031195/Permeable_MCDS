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
import numpy as np
import pandas as pd
import random

def make_groups(indices, group_size, num_groups, with_replacement=False, rng=None):
    rng = rng or random
    groups = []
    pool = list(indices)
    if with_replacement:
        for _ in range(num_groups):
            groups.append([rng.choice(pool) for _ in range(group_size)])
    else:
        if group_size * num_groups > len(pool):
            raise ValueError("Not enough reps for non-overlapping groups")
        rng.shuffle(pool)
        for g in range(num_groups):
            groups.append(pool[g*group_size:(g+1)*group_size])
    return groups

def create_group_averages(df, group_size, num_groups, total_reps, walkers_count,
                          with_replacement=False, rng=None):
    
    # make 10 random groups of given size
    groups = make_groups(list(range(total_reps)), group_size, num_groups,
                         with_replacement=with_replacement, rng=rng)


    all_group_averages = []
    for i, group in enumerate(groups):
        df_group = df[df["rep"].isin(group)].copy()

        
        avg = df_group.groupby("bvals", as_index=False).mean(numeric_only=True)
        # normalize within group across bvals (powder-normalization)
        max_per = avg["DWI"].max()
        avg["norm_DWI"] = avg["DWI"] / max_per if max_per > 0 else avg["DWI"]
        avg["group_id"] = i
        avg["nbr_walkers"] = walkers_count
        all_group_averages.append(avg)

    return pd.concat(all_group_averages, ignore_index=True)

def repeated_group_averages(df, group_size, num_groups, total_reps, walkers_count,
                            n_repeats=100, with_replacement=False, seed=0):
    """
    For each repeat:
      - build num_groups groups,
      - compute per-group normalized DWI by bval,
      - compute mean across groups and std across groups for each bval.
    Return both per-repeat stats and an aggregated summary.
    """
    per_repeat = []
    rng = random.Random(seed)

    for r in range(n_repeats):
        # advance RNG reproducibly
        local_rng = random.Random(rng.random())
        grp_avg = create_group_averages(df, group_size, num_groups, total_reps, walkers_count,
                                        with_replacement=with_replacement, rng=local_rng)
        # stats across groups within this repeat
        stats = grp_avg.groupby(["bvals", "nbr_walkers"], as_index=False).agg(
            mean_DWI=("norm_DWI", "mean"),
            std_groups=("norm_DWI", "std")
        )
        stats["cov"] = stats["std_groups"] / stats["mean_DWI"]
        stats["repeat"] = r
        per_repeat.append(stats)

    per_repeat_df = pd.concat(per_repeat, ignore_index=True)

    # aggregate across repeats
    summary = per_repeat_df.groupby(["nbr_walkers", "bvals"], as_index=False).agg(
        mean_of_means=("mean_DWI", "mean"),
        se_of_means=("mean_DWI", lambda x: x.std(ddof=1)/np.mean(x)),
        mean_of_stds=("std_groups", "mean"),
        se_of_stds=("std_groups",  lambda x: x.std(ddof=1)/np.mean(x)),
        cov=("cov", "mean"),
    )

    return per_repeat_df, summary


# === MAIN CODE ===
if __name__ == "__main__":
    path_to_scheme = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/SMI.scheme"
    path_to_folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/complex_axons/simulations/f_0.5_extra/"

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

    f = 0.5

    nbr_walker_per_rep = 133000*f

    nbr_group = 10
    per_rep_1, summary_1 = repeated_group_averages(
        df, group_size=1, num_groups=nbr_group,
        total_reps=len(all_files),
        walkers_count=nbr_walker_per_rep,
        n_repeats=100, with_replacement=True, seed=42
    )
    per_rep_2, summary_2 = repeated_group_averages(
        df, group_size=2, num_groups=nbr_group,
        total_reps=len(all_files),
        walkers_count=2*nbr_walker_per_rep,
        n_repeats=100, with_replacement=True, seed=42
    )
    per_rep_5, summary_5 = repeated_group_averages(
        df, group_size=5, num_groups=nbr_group,
        total_reps=len(all_files),
        walkers_count=5*nbr_walker_per_rep,
        n_repeats=100, with_replacement=True, seed=42
    )
    per_rep_7, summary_7 = repeated_group_averages(
        df, group_size=7, num_groups=nbr_group,
        total_reps=len(all_files),
        walkers_count=7*nbr_walker_per_rep,
        n_repeats=100, with_replacement=True, seed=42
    )
    per_rep_10, summary_10 = repeated_group_averages(
        df, group_size=10, num_groups=nbr_group,
        total_reps=len(all_files),
        walkers_count=10*nbr_walker_per_rep,
        n_repeats=100, with_replacement=True, seed=42
    )
    per_rep_15, summary_15 = repeated_group_averages(
        df, group_size=15, num_groups=nbr_group,
        total_reps=len(all_files),
        walkers_count=15*nbr_walker_per_rep,
        n_repeats=100, with_replacement=True, seed=42
    )
    per_rep_20, summary_20 = repeated_group_averages(
        df, group_size=20, num_groups=nbr_group,
        total_reps=len(all_files),
        walkers_count=20*nbr_walker_per_rep,
        n_repeats=100, with_replacement=True, seed=42
    )
    per_rep_30, summary_30 = repeated_group_averages(
        df, group_size=30, num_groups=nbr_group,
        total_reps=len(all_files),
        walkers_count=30*nbr_walker_per_rep,
        n_repeats=100, with_replacement=True, seed=42
    )



    final_df = pd.concat([summary_1, summary_2, summary_5, summary_7, summary_10, summary_15, summary_20, summary_30], ignore_index=True)

    voxel_length = 110  # μm
    f = 0.5
    final_df["walkers/μm³"] = final_df["nbr_walkers"] / (voxel_length**3 * f)
    final_df["bvals_label"] = final_df["bvals"].apply(lambda b: f"{b} ms/µm²")


    # find elbow in final_df for cov vs walkers/μm³ for b=5000
    from kneed import KneeLocator
    df_b1000 = final_df[final_df["bvals"] == 5]
    print(df_b1000)
    kn = KneeLocator(df_b1000["walkers/μm³"], df_b1000["cov"], curve='convex', direction='decreasing')
    elbow_point = kn.knee

    # Plot mean
    plt.figure(figsize=(10,6))
    sns.lineplot(data=final_df, x="walkers/μm³", y="mean_of_means", hue="bvals_label", marker="o", ci=None, palette="tab10", linewidth=2.5)
    plt.xlabel("Density of Walkers")
    plt.ylabel("Mean powder-averaged DWI")
    plt.title("Mean Signal vs Number of Walkers")
    # x axis font size
    fontsize=14
    plt.xticks(fontsize=fontsize); plt.yticks(fontsize=fontsize)
    # legend font size
    plt.grid(True); plt.legend(fontsize= fontsize).set_title(""); plt.tight_layout(); plt.show()

    # Plot std across groups (averaged over repeats)
    plt.figure(figsize=(10,6))
    sns.lineplot(data=final_df, x="walkers/μm³", y="cov", hue="bvals_label", marker="o", ci=None, palette="tab10", linewidth=2.5)
    # highlight elbow point
    plt.axvline(x=elbow_point, color='grey', linestyle='--', label=f'Elbow at {elbow_point:.2f} walkers/μm³')
    plt.xlabel("Density of Walkers (walkers/μm³)")
    plt.ylabel(f"CV across {nbr_group} groups (powder-avg DWI)")
    plt.title("Monte-Carlo Variability vs Number of Walkers")
    plt.xticks(fontsize=fontsize); plt.yticks(fontsize=fontsize)
    plt.grid(True); plt.legend(fontsize= fontsize).set_title(""); plt.tight_layout(); plt.show()