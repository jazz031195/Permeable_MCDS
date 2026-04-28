import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import dipy.reconst.dki as dki
import nibabel as nib
import glob
from useful_functions import read_binary_file, get_scheme_info, get_files_from_folder
from diffusion_from_trajectories import diffusion
from scipy.stats import pearsonr
from scipy.stats import pearsonr, ttest_ind
from scipy.stats import ttest_1samp
nbr_subjects = 13

def create_noisy_data():
    # ------------------------
    # Paths & configuration
    # ------------------------
    directory = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/ADC_blood/"
    folders = ["Active_Voxel", "reduced_blood_vessel_radius", "reduced_axons_radius", "reduced_all_radius"]
    compartments = ["intra", "extra"]

    scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/ADC_only_x.scheme"
    b_values, directions = get_scheme_info(scheme_path)

    SNR = 50
    n_synthetic_rois = 1000  # Number of final averaged data points per repetition
    voxels_per_roi = 10      # Number of voxels averaged to boost SNR
    

    avg_samples = True

    # ------------------------
    # Load all DWI data
    # ------------------------
    rows = []

    for folder in folders:
        print("Processing folder:", folder)
        for compartment in compartments:
            print("  Compartment:", compartment)

            subfolder = f"{folder}_{compartment}/"
            path = os.path.join(directory, subfolder)
            print("    Path:", path)

            DWI_files = get_files_from_folder(path, binary=True)
            print("    Found", len(DWI_files), "DWI files.")

            traj_files = glob.glob(os.path.join(path, "*.traj"))

            header_files = [traj_file.split(".traj")[0] + ".bhdr" for traj_file in traj_files]
            print(header_files)


            nbr_paralell_runs = [int(traj_file.split(".traj")[0].split("_")[-1]) for traj_file in traj_files]
            nbr_paralell_runs = np.max(nbr_paralell_runs) + 1

            for i, DWI_file in enumerate(DWI_files):
                # skip .img sidecars etc., if present
                if "img" in DWI_file:
                    continue

                dwi = read_binary_file(DWI_file)

                prefix = DWI_file.split("DWI.bfloat")[0]
                print("    Processing DWI file:", DWI_file)

                diffusion_time = 0.07

                axial_diffusions = []
                radial_diffusions = []

                for header_file, traj_file in zip(header_files, traj_files):
                    if "rep" not in prefix and "rep" in traj_file:
                        continue
                    elif prefix not in traj_file:
                        continue
                    print("************************")
                    bhdr = np.fromfile(header_file, dtype="float32")
                    bhdr = bhdr[np.isfinite(bhdr)]

                    nbr_walkers = int(bhdr[1])
                    prefix_2 = traj_file.split(".traj")[0]
                    diff_file = prefix_2 + f"diffusion_AD_RD.txt"
                    if not os.path.exists(diff_file):
                        AD, RD = diffusion(traj_file, nbr_walkers, diffusion_time, compartment)
                        np.savetxt(diff_file, np.array([AD, RD]))
                    else:
                        AD, RD = np.loadtxt(diff_file)

                    axial_diffusions.append(AD)
                    radial_diffusions.append(RD)
                mean_AD = np.mean(axial_diffusions)
                mean_RD = np.mean(radial_diffusions)
                # sanity check: number of values must match number of b-values
                if len(dwi) != len(b_values):
                    raise ValueError(
                        f"File {DWI_file}: len(dwi)={len(dwi)} != len(b_values)={len(b_values)}"
                    )

                for b, s in zip(b_values, dwi):
                    rows.append({
                        "Group": folder,
                        "Compartment": compartment,
                        "Repetition": i,
                        "b_value": b,
                        "DWI": s,
                        "Axial_Diffusion_traj": mean_AD,
                        "Radial_Diffusion_traj": mean_RD
                    })

    df = pd.DataFrame(rows)
    print("Total rows loaded:", len(df))
    

    # sum intra and extra compartments for each group and repetition
    df = (
        df.groupby(["Group", "b_value"], as_index=False)
        .agg({"DWI": "sum", "Compartment": lambda x: "intra+extra", "Axial_Diffusion_traj" : "mean", "Radial_Diffusion_traj": "mean", "Repetition": "first"})
    )

    del df["Compartment"]  # we only have intra+extra now, so this column is redundant
    del df["Repetition"]  # will be created later after merging b0 means
    # ------------------------
    # Normalise by b=0 signal
    # ------------------------
    # Extract b0 and compute mean per (Group, Compartment, Repetition)
    df_b0 = (
        df[df["b_value"] == 0.0]
        .groupby(["Group", "Axial_Diffusion_traj", "Radial_Diffusion_traj"], as_index=False)["DWI"]
        .mean()
        .rename(columns={"DWI": "DWI_b0"})
    )

    # Merge back
    df = df.merge(
        df_b0,
        on=["Group",  "Axial_Diffusion_traj", "Radial_Diffusion_traj"],
        how="left"
    )

    # Normalised signal
    df["DWI_normalized"] = df["DWI"] / df["DWI_b0"]

    # Track original row index for later pairing
    df["OrigRow"] = np.arange(len(df))

    # ------------------------
    # Monte Carlo noise simulation (Rician)
    # ------------------------
    # Define how many final data points you want, and how many voxels make up an ROI
    
    total_samples_per_row = n_synthetic_rois * voxels_per_roi # e.g., 10,000

    # Replicate each row 
    df_rep = df.loc[df.index.repeat(total_samples_per_row)].copy()
    df_rep.reset_index(drop=True, inplace=True)

    # Assign a unique sample ID for noise generation
    df_rep["Sample"] = df_rep.groupby("OrigRow").cumcount()
    
    # Assign an ROI ID to group them into chunks of 10
    # Samples 0-9 become ROI 0, 10-19 become ROI 1, etc.
    df_rep["ROI_ID"] = df_rep["Sample"] // voxels_per_roi

    # Add complex Gaussian noise (σ = 1/SNR)
    sigma = 1.0 / SNR
    n = len(df_rep)
    df_rep["Real_Noise"] = np.random.normal(0, sigma, n)
    df_rep["Imag_Noise"] = np.random.normal(0, sigma, n)

    df_rep["Noisy_DWI"] = np.sqrt(
        (df_rep["DWI_normalized"] + df_rep["Real_Noise"])**2 +
        df_rep["Imag_Noise"]**2
    )

    # ------------------------
    # Compute ADC between two b-values (e.g. 0.2 and 1.0)
    # ------------------------
    b1 = 0.2
    b2 = 1.0

    df_b1 = df_rep[df_rep["b_value"] == b1].copy()
    df_b2 = df_rep[df_rep["b_value"] == b2].copy()

    # Merge signals at b1 and b2 for the same measurement + sample
    merge_keys = ["Group", "Sample", "ROI_ID", "Axial_Diffusion_traj", "Radial_Diffusion_traj"]
    df_adc = df_b1.merge(
        df_b2,
        on=merge_keys,
        suffixes=(f"_b{b1}", f"_b{b2}")
    )

    S1 = df_adc[f"Noisy_DWI_b{b1}"]
    S2 = df_adc[f"Noisy_DWI_b{b2}"]

    # ADC = -ln(S(b1)/S(b2)) / (b1 - b2)
    denom = (b1 - b2)
    df_adc["ADC"] = np.where(
        S2 > 0,
        -np.log(S1 / S2) / denom,
        np.nan
    )

    # Keep relevant columns
    df_adc = df_adc[["Group",  "ROI_ID", "ADC", "Axial_Diffusion_traj", "Radial_Diffusion_traj"]]

    if avg_samples:
        # ------------------------
        # Collapse the 10 voxels within each synthetic ROI
        # ------------------------
        # This groups by ROI_ID, averaging the 10 ADC values to boost SNR
        df_adc = (
            df_adc
            .groupby(["Group",  "ROI_ID", "Axial_Diffusion_traj", "Radial_Diffusion_traj"], as_index=False)["ADC"]
            .mean()
        )

    # ------------------------
    # Relabel groups for plotting
    # ------------------------
    df_adc["Group"] = df_adc["Group"].map({
        "Active_Voxel": "Active",
        "reduced_blood_vessel_radius": "Swollen axons",
        "reduced_axons_radius": "Swollen vessels",
        "reduced_all_radius": "Rest"
    })

    # ------------------------
    # Save & plot
    # ------------------------
    out_csv = os.path.join(directory, "ADC_blood_all_data.csv")
    df_rep.to_csv(out_csv, index=False)
    print("Saved full noisy dataset to:", out_csv)

    # reorder the categories for plotting
    df_adc["Group"] = pd.Categorical(df_adc["Group"], categories=["Active", "Swollen axons", "Swollen vessels", "Rest"], ordered=True)

    # ------------------------
    # Save & plot
    # ------------------------
    out_csv = os.path.join(directory, "ADC_blood_all_data.csv")
    df_rep.to_csv(out_csv, index=False)
    print("Saved full noisy dataset to:", out_csv)

    return df_adc



def simulate_adc_connectivity(df, value_col, active_group="Swollen axons", rest_group="Rest", 
                              n_timepoints=400, iterations=1000):
    """
    Simulates two independent time series drawing from empirical ADC distributions.
    Compares the Pearson r of time series sharing a temporal pattern vs independent patterns.
    """
    active_dist = df[df["Group"] == active_group][value_col].dropna().values
    rest_dist = df[df["Group"] == rest_group][value_col].dropna().values
    
    # ==========================================
    # SCENARIO 1: SHARED TEMPORAL PATTERN
    # ==========================================
    corr_shared = []
    
    for _ in range(iterations):
        sequence_labels = np.random.choice([active_group, rest_group], size=n_timepoints)
        active_idx = np.where(sequence_labels == active_group)[0]
        rest_idx = np.where(sequence_labels == rest_group)[0]
        
        ts1 = np.zeros(n_timepoints)
        ts2 = np.zeros(n_timepoints)
        
        ts1[active_idx] = np.random.choice(active_dist, size=len(active_idx), replace=True)
        ts1[rest_idx] = np.random.choice(rest_dist, size=len(rest_idx), replace=True)
        
        ts2[active_idx] = np.random.choice(active_dist, size=len(active_idx), replace=True)
        ts2[rest_idx] = np.random.choice(rest_dist, size=len(rest_idx), replace=True)
        
        r_val, _ = pearsonr(ts1, ts2)
        corr_shared.append(r_val)
        
    corr_shared = np.array(corr_shared)

    print(f"--- Shared Pattern ({iterations} iterations) ---")
    print(f"Mean r: {np.mean(corr_shared):.4f} ± {np.std(corr_shared):.4f}")


    # ==========================================
    # SCENARIO 2: INDEPENDENT TEMPORAL PATTERNS
    # ==========================================
    corr_independent = []
    
    for _ in range(iterations):
        sequence_labels1 = np.random.choice([active_group, rest_group], size=n_timepoints)
        sequence_labels2 = np.random.choice([active_group, rest_group], size=n_timepoints)
        
        # BUG FIXED HERE (added '1' to sequence_labels)
        active_idx1 = np.where(sequence_labels1 == active_group)[0]
        rest_idx1 = np.where(sequence_labels1 == rest_group)[0]

        active_idx2 = np.where(sequence_labels2 == active_group)[0]
        rest_idx2 = np.where(sequence_labels2 == rest_group)[0]

        ts1 = np.zeros(n_timepoints)
        ts2 = np.zeros(n_timepoints)

        ts1[active_idx1] = np.random.choice(active_dist, size=len(active_idx1), replace=True)
        ts1[rest_idx1] = np.random.choice(rest_dist, size=len(rest_idx1), replace=True)

        ts2[active_idx2] = np.random.choice(active_dist, size=len(active_idx2), replace=True)
        ts2[rest_idx2] = np.random.choice(rest_dist, size=len(rest_idx2), replace=True)

        r_val, _ = pearsonr(ts1, ts2)
        corr_independent.append(r_val)
    
    corr_independent = np.array(corr_independent)

    print(f"--- Independent Patterns ({iterations} iterations) ---")
    print(f"Mean r: {np.mean(corr_independent):.4f} ± {np.std(corr_independent):.4f}\n")


    # ==========================================
    # STATISTICAL TESTING (Fisher Z & T-Test)
    # ==========================================
    # 1. Fisher Z-transform the correlation coefficients
    z_shared = np.arctanh(corr_shared)
    z_independent = np.arctanh(corr_independent)
    
    # 2. Independent t-test (one-sided: is shared > independent?)
    t_stat, p_val = ttest_ind(z_shared, z_independent, alternative='greater')
    
    # 3. Calculate Cohen's d for effect size
    pooled_std = np.sqrt((np.std(z_shared, ddof=1)**2 + np.std(z_independent, ddof=1)**2) / 2)
    cohens_d = (np.mean(z_shared) - np.mean(z_independent)) / pooled_std

    print(f"--- Statistical Comparison ---")
    print(f"T-statistic: {t_stat:.4f}")
    print(f"P-value:     {p_val:.4e}")
    print(f"Cohen's d:   {cohens_d:.4f}")
    
    if p_val < 0.05:
        print(">> RESULT: The shared temporal pattern yields significantly higher functional connectivity.")
    else:
        print(">> RESULT: No significant difference detected between shared and independent patterns.")


    return corr_shared, corr_independent


def plot_correlation_distributions(corr_shared, corr_independent, cohens_d=None):
    """
    Plots overlaid kernel density estimates (KDE) and histograms for the 
    shared vs independent temporal pattern simulations.
    """
    # 1. Package into a DataFrame
    df_plot = pd.DataFrame({
        "Pearson_r": np.concatenate([corr_shared, corr_independent]),
        "Condition": ["Shared Pattern"] * len(corr_shared) + ["Independent Pattern"] * len(corr_independent)
    })

    # 2. Set up a professional academic plot style
    sns.set_theme(style="ticks", context="paper", font_scale=1.2)
    fig, ax = plt.subplots(figsize=(8, 5))

    # 3. Plot the distributions (Histogram + KDE curve)
    sns.histplot(
        data=df_plot, 
        x="Pearson_r", 
        hue="Condition", 
        palette={"Shared Pattern": "#d62728", "Independent Pattern": "#1f77b4"},
        fill=True, 
        kde=True,          # Adds the smooth curve
        line_kws={'linewidth': 2},
        alpha=0.4, 
        edgecolor=None,
        ax=ax
    )

    # 4. Calculate and plot the means as vertical dashed lines
    mean_shared = np.mean(corr_shared)
    mean_indep = np.mean(corr_independent)
    
    ax.axvline(mean_indep, color="#1f77b4", linestyle="--", linewidth=2)
    ax.axvline(mean_shared, color="#d62728", linestyle="--", linewidth=2)
    ax.axvline(0, color="black", linestyle="-", linewidth=1, alpha=0.5) # The zero-line

    # 5. Formatting for thesis-grade aesthetics
    ax.set_title("Simulated Functional Connectivity: Shared vs. Independent Patterns", fontweight='bold', pad=15)
    ax.set_xlabel("Pearson Correlation ($r$)", fontweight='bold')
    ax.set_ylabel("Density", fontweight='bold')
    
    # 6. Add text box with statistics if provided
    if cohens_d is not None:
        stats_text = f"Mean Shift: {mean_shared - mean_indep:.3f}\nCohen's d: {cohens_d:.2f}"
        ax.text(0.95, 0.85, stats_text, 
                transform=ax.transAxes, 
                fontsize=11, 
                verticalalignment='top', 
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

    # Clean up the legend
    sns.move_legend(ax, "upper left", title=None, frameon=True)
    sns.despine()

    plt.tight_layout()
    plt.show()

def simulate_adc_connectivity_multiple_r(df, value_col, active_group="Swollen axons", rest_group="Rest", 
                                         n_timepoints=400, iterations=10000):
    """
    Simulates time series with specific underlying pattern correlations (1.0, 0.8, 0.6, 0.4, 0.0).
    Outputs and compares the Fisher Z-transformed distributions.
    """
    active_dist = df[df["Group"] == active_group][value_col].dropna().values
    rest_dist = df[df["Group"] == rest_group][value_col].dropna().values
    
    target_correlations = [1.0, 0.8, 0.6, 0.4, 0.0]
    results = {}
    
    for r_target in target_correlations:
        z_list = []
        flip_prob = (1.0 - r_target) / 2.0 
        
        for _ in range(iterations):

            zvals = []

            for i in range(nbr_subjects):

                # 1. Generate Base Sequence
                seq1 = np.random.choice([active_group, rest_group], size=n_timepoints)
                
                # 2. Generate Sequence 2 by flipping bits
                flips = np.random.rand(n_timepoints) < flip_prob
                seq2 = np.copy(seq1)
                for i in range(n_timepoints):
                    if flips[i]:
                        seq2[i] = rest_group if seq1[i] == active_group else active_group
                
                # 3. Find indices
                a_idx1, r_idx1 = np.where(seq1 == active_group)[0], np.where(seq1 == rest_group)[0]
                a_idx2, r_idx2 = np.where(seq2 == active_group)[0], np.where(seq2 == rest_group)[0]
                
                # 4. Sample the empirical noise
                ts1, ts2 = np.zeros(n_timepoints), np.zeros(n_timepoints)
                
                ts1[a_idx1] = np.random.choice(active_dist, size=len(a_idx1), replace=True)
                ts1[r_idx1] = np.random.choice(rest_dist, size=len(r_idx1), replace=True)
                
                ts2[a_idx2] = np.random.choice(active_dist, size=len(a_idx2), replace=True)
                ts2[r_idx2] = np.random.choice(rest_dist, size=len(r_idx2), replace=True)
                
                # 5. Measure Pearson r, safely clip it, and immediately transform to Fisher Z
                r_val, _ = pearsonr(ts1, ts2)
                r_val_clipped = np.clip(r_val, -0.9999, 0.9999)
                z_val = np.arctanh(r_val_clipped)

                zvals.append(z_val)
            
            mean_zval = np.mean(zvals)
            
            z_list.append(mean_zval)
            
        # Store the Fisher Z values directly
        results[f"Underlying r={r_target}"] = np.array(z_list)
        
        print(f"--- Underlying Pattern r = {r_target} ---")
        print(f"Measured mean Fisher z: {np.mean(z_list):.4f} ± {np.std(z_list):.4f}\n")

    # ==========================================
    #   STATISTICAL SIGNIFICANCE & EFFECT SIZE
    # ==========================================
    print("\n==========================================")
    print("   STATISTICAL SIGNIFICANCE & EFFECT SIZE")
    print("==========================================")
    
    # The distributions are ALREADY Fisher z-transformed!
    z_zero_dist = results["Underlying r=0.0"]
    
    for r_target in target_correlations:
        if r_target == 0.0:
            continue
            
        z_target_dist = results[f"Underlying r={r_target}"]
        
        # 1. Perform the one-sided T-test
        t_stat, p_val = ttest_ind(z_target_dist, z_zero_dist, alternative='greater')
        
        # 2. Calculate Cohen's d
        var_target = np.std(z_target_dist, ddof=1)**2
        var_zero = np.std(z_zero_dist, ddof=1)**2
        pooled_std = np.sqrt((var_target + var_zero) / 2)
        
        cohens_d = (np.mean(z_target_dist) - np.mean(z_zero_dist)) / pooled_std
        
        # 3. Interpret effect size
        if cohens_d >= 0.8:
            effect_label = "Strong / Robust"
        elif cohens_d >= 0.5:
            effect_label = "Moderate"
        elif cohens_d >= 0.2:
            effect_label = "Small"
        else:
            effect_label = "Negligible"

        print(f"Condition: True Neural r={r_target} vs. Noise (r=0.0)")
        print(f"  -> T-statistic: {t_stat:.4f}")
        print(f"  -> P-value:     {p_val:.4e}")
        print(f"  -> Cohen's d:   {cohens_d:.4f} ({effect_label})\n")

    # ==========================================
    #   TYPE I & TYPE II ERROR ANALYSIS
    # ==========================================
    print("\n==========================================")
    print("   TYPE I & TYPE II ERROR ANALYSIS")
    print("==========================================")
    
    alpha = 0.05 
    
    # Calculate critical threshold directly on the Z null distribution
    z_critical = np.percentile(z_zero_dist, 100 * (1 - alpha))
 
    print(f"Fixed Alpha (Type I Error Rate): {alpha:.2f} (5%)")
    print(f"Critical Threshold required for detection: Fisher z = {z_critical:.4f}\n")
    
    # Calculate Beta and Power for each True correlation
    for r_target in target_correlations:
        z_target_dist = results[f"Underlying r={r_target}"]
        
        # Beta is the proportion that falls BELOW the critical threshold
        beta = np.mean(z_target_dist < z_critical)
        power = 1.0 - beta
        
        print(f"Condition: True Neural r = {r_target}")
        print(f"  -> Beta (Type II Error): {beta:.4f}  ({beta*100:.1f}% false negative rate)")
        print(f"  -> Statistical Power:    {power:.4f}  ({power*100:.1f}% detection rate)\n")

    return results

def plot_multiple_correlation_distributions(results_dict):
    """
    Plots the density distributions for all simulated correlation levels.
    """
    # 1. Unpack dictionary into a format Seaborn loves (long-form DataFrame)
    data = []
    for condition_name, r_values in results_dict.items():
        for r in r_values:
            data.append({"Measured Pearson r": r, "Condition": condition_name})
    df_plot = pd.DataFrame(data)

    # 2. Set up the plot
    sns.set_theme(style="ticks", context="paper", font_scale=1.2)
    fig, ax = plt.subplots(figsize=(9, 6))

    # A nice color palette progressing from Blue (noise) to Red (perfectly shared)
    palette = {
        "Underlying r=1.0": "#d62728", # Red
        "Underlying r=0.8": "#ff7f0e", # Orange
        "Underlying r=0.6": "#2ca02c", # Green
        "Underlying r=0.4": "#9467bd", # Purple
        "Underlying r=0.0": "#1f77b4"  # Blue
    }

    # 3. Plot the KDE density curves
    sns.kdeplot(
        data=df_plot, 
        x="Measured Pearson r", 
        hue="Condition", 
        palette=palette,
        fill=True, 
        alpha=0.3, 
        linewidth=2.5,
        ax=ax
    )

    # 4. Add vertical mean lines
    for condition, color in palette.items():
        mean_val = df_plot[df_plot["Condition"] == condition]["Measured Pearson r"].mean()
        ax.axvline(mean_val, color=color, linestyle="--", linewidth=1.5)

    ax.axvline(0, color="black", linestyle="-", linewidth=1, alpha=0.5) # The zero-line

    # 5. Formatting
    ax.set_title("Measured ADC Connectivity by Underlying Neural Correlation", fontweight='bold', pad=15)
    ax.set_xlabel("Measured Fisher Transformed Pearson Correlation ($z$)", fontweight='bold')
    ax.set_ylabel("Probability Density", fontweight='bold')
    
    sns.move_legend(ax, "upper right", title="True Signal Pattern", frameon=True)
    sns.despine()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":


    df_adc = create_noisy_data()
    results_dict = simulate_adc_connectivity_multiple_r(df_adc, value_col="ADC")
    plot_multiple_correlation_distributions(results_dict)

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df_adc, x="Group", y="ADC")
    plt.title("ADC values by Group and Compartment")
    plt.ylabel("ADC (μm²/ms)")
    plt.xlabel("Group")
    plt.legend(title="Compartment")
    plt.tight_layout()
    plt.show()

