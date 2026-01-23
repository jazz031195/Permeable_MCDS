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
if __name__ == "__main__":

    # ------------------------
    # Paths & configuration
    # ------------------------
    directory = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/ADC_blood/"
    folders = ["Active_Voxel", "reduced_blood_vessel_radius", "reduced_axons_radius", "reduced_all_radius"]
    compartments = ["intra", "extra"]

    scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/ADC_only_x.scheme"
    b_values, directions = get_scheme_info(scheme_path)

    SNR = 1000000000
    nbr_samples = 1000  # Monte Carlo repeats per measurement

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
                        "Axial_Diffusion": mean_AD,
                        "Radial_Diffusion": mean_RD
                    })

    df = pd.DataFrame(rows)
    print("Total rows loaded:", len(df))

    # ------------------------
    # Normalise by b=0 signal
    # ------------------------
    # Extract b0 and compute mean per (Group, Compartment, Repetition)
    df_b0 = (
        df[df["b_value"] == 0.0]
        .groupby(["Group", "Compartment", "Repetition", "Axial_Diffusion", "Radial_Diffusion"], as_index=False)["DWI"]
        .mean()
        .rename(columns={"DWI": "DWI_b0"})
    )

    # Merge back
    df = df.merge(
        df_b0,
        on=["Group", "Compartment", "Repetition", "Axial_Diffusion", "Radial_Diffusion"],
        how="left"
    )

    # Normalised signal
    df["DWI_normalized"] = df["DWI"] / df["DWI_b0"]

    # Track original row index for later pairing
    df["OrigRow"] = np.arange(len(df))

    # ------------------------
    # Monte Carlo noise simulation (Rician)
    # ------------------------
    # Replicate each row nbr_samples times
    df_rep = df.loc[df.index.repeat(nbr_samples)].copy()
    df_rep.reset_index(drop=True, inplace=True)

    # Sample index within each original measurement
    df_rep["Sample"] = df_rep.groupby("OrigRow").cumcount()

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
    merge_keys = ["Group", "Compartment", "Repetition", "Sample", "Axial_Diffusion", "Radial_Diffusion"]
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
    df_adc = df_adc[["Group", "Compartment", "Repetition", "Sample", "ADC", "Axial_Diffusion", "Radial_Diffusion"]]

    if avg_samples:
        # ------------------------
        # Optional: collapse Monte Carlo samples (e.g. mean ADC per repetition)
        # ------------------------
        # If you want to keep the full distribution, comment this out.
        df_adc = (
            df_adc
            .groupby(["Group", "Compartment", "Repetition", "Axial_Diffusion", "Radial_Diffusion"], as_index=False)["ADC"]
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

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df_adc, x="Group", y="ADC", hue="Compartment" )
    plt.title("ADC values by Group and Compartment")
    plt.ylabel("ADC (μm²/ms)")
    plt.xlabel("Group")
    plt.legend(title="Compartment")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df_adc, x="Group", y="Axial_Diffusion", hue="Compartment" )
    plt.title("Axial Diffusion values by Group and Compartment")
    plt.ylabel("Axial Diffusion (μm²/ms)")
    plt.xlabel("Group")
    plt.legend(title="Compartment")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df_adc, x="Group", y="Radial_Diffusion", hue="Compartment" )
    plt.title("Radial Diffusion values by Group and Compartment")
    plt.ylabel("Radial Diffusion (μm²/ms)")
    plt.xlabel("Group")
    plt.legend(title="Compartment")
    plt.tight_layout()
    plt.show()