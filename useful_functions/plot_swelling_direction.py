import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings

warnings.filterwarnings("ignore")

from useful_functions import read_binary_file, get_files_from_folder, get_scheme_info, get_scheme_info_iso


def calculate_ADC(b0 ,b1, scheme_path, path_to_data):
    
    all_b_values = [] 
    all_directions = []
    all_DWIs = []
    all_swellings = []
    isos = []
    files = get_files_from_folder(path_to_data)
    for file in files:
        if "img" not in file and "info" not in file:
            print(file)
            DWI = read_binary_file(file) 
            DWI =[float(i) for i in DWI] 
            if "_0_" in file:
                all_swellings.extend(np.zeros(len(DWI)))
            elif "_0.25_" in file:
                all_swellings.extend(np.ones(len(DWI))*0.25)
            elif "_0.5_" in file:
                all_swellings.extend(np.ones(len(DWI))*0.5)
            elif "_0.75_" in file:
                all_swellings.extend(np.ones(len(DWI))*0.75)
            elif "_1_" in file:
                all_swellings.extend(np.ones(len(DWI))*1)
            else:
                print("Error, no swelling found")
                assert(0)
            if "iso" not in file:
                if "_1_" in file:
                    scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/PGSE_22_dir_12_b.scheme"
                else:
                    scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/PGSE_21_dir_12_b.scheme"
                b_values, directions = get_scheme_info(scheme_path)
                isos.extend([False]*len(directions))
            else:
                scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/iso_waveform_vec_b.txt"
                b_values, directions = get_scheme_info_iso(scheme_path)
                isos.extend([True]*len(directions))
            all_DWIs.extend(DWI)
            all_b_values.extend(b_values)
            all_directions.extend(directions)
    df = pd.DataFrame({"DWI": all_DWIs, "b_value": all_b_values, "direction": all_directions, "swelling": all_swellings, "waveform": isos})

    # calculate angle between direction and vector (0,0,1)
    df["angle (rad)"] = df["direction"].apply(lambda x: np.arccos(np.dot(x, [0,0,1])/np.linalg.norm(x)))
    
    # angle must be from 0 to pi
    df["angle (rad)"] = df["angle (rad)"].apply(lambda x: x if x <= np.pi/2 else np.pi-x)
    
    df_waveform = df.loc[df["waveform"] == True]
    df_pgse = df.loc[df["waveform"] == False]


    df_b0 = df_pgse.loc[df_pgse["b_value"] == b0]
    df_b0 = df_b0.groupby(["angle (rad)", "swelling"]).sum()
    df_b0 = df_b0.groupby(["angle (rad)", "swelling"]).sum()
    df_b1 = df_pgse.loc[df_pgse["b_value"] == b1]
    df_b1 = df_b1.groupby(["angle (rad)", "swelling"]).sum()
    df_b1 = df_b1.groupby(["angle (rad)", "swelling"]).sum()
    # calculate adc list form the two b values
    adc = np.log(df_b0["DWI"].values/df_b1["DWI"].values)/(b1-b0)
    df_b0["ADC"] = adc
    df_b0["ADC"] = df_b0["ADC"].replace([np.inf, -np.inf], np.nan)
    df_b0 = df_b0.dropna().reset_index()
    df_b0 = df_b0.loc[df_b0["angle (rad)"] > 0]

    df_b0_waveform = df_waveform.loc[df_waveform["b_value"] == b0]
    df_b0_waveform = df_b0_waveform.groupby(["angle (rad)", "swelling"]).sum()
    df_b0_waveform = df_b0_waveform.groupby(["angle (rad)", "swelling"]).sum()
    df_b1 = df_waveform.loc[df_waveform["b_value"] == b1]
    df_b1 = df_b1.groupby(["angle (rad)", "swelling"]).sum()
    df_b1 = df_b1.groupby(["angle (rad)", "swelling"]).sum()
    # calculate adc list form the two b values
    adc = np.log(df_b0_waveform["DWI"].values/df_b1["DWI"].values)/(b1-b0)
    df_b0_waveform["ADC"] = adc
    df_b0_waveform["ADC"] = df_b0_waveform["ADC"].replace([np.inf, -np.inf], np.nan)
    df_b0_waveform = df_b0_waveform.dropna().reset_index()
    df_b0_waveform = df_b0_waveform.loc[df_b0_waveform["angle (rad)"] > 0]

    return df_b0, df_b0_waveform

def plot_with_respect_to_direction(df_pgse, df_waveform):

    # multiplie swelling by 100 to get percentage 
    df_pgse["swelling (%)"] = df_pgse["swelling"]
    df_waveform["swelling (%)"] = df_waveform["swelling"]

    # Calculate ADC_relative
    def calculate_relative(group):
        base_adc = group.loc[group['swelling'] == 0, 'ADC'].values[0]
        group['ADC_relative'] = group['ADC'] / base_adc
        return group

    df_pgse = df_pgse.groupby('angle (rad)').apply(calculate_relative)
    df_waveform = df_waveform.groupby('angle (rad)').apply(calculate_relative)
    print(df_pgse)
    print(df_waveform)
    # save data
    folder = "/home/localadmin/Documents/CATERPillar/arthurs_analysis"
    df_pgse.to_csv(f"{folder}/df_pgse.csv")
    df_waveform.to_csv(f"{folder}/df_waveform.csv")
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.5)
    sns.lmplot(hue="swelling (%)", y="ADC_relative", x="angle (rad)", data=df_pgse)
    plt.ylabel('Relative ADC')
    plt.xlabel('Angle (rad)')
    plt.title('Relative ADC decrease with swelling')
    plt.show()

    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.5)
    sns.lmplot(hue="swelling (%)", y="ADC_relative", x="angle (rad)", data=df_waveform)
    plt.ylabel('Relative ADC')
    plt.xlabel('Angle (rad)')
    plt.title('Relative ADC decrease with swelling')
    plt.show()


def main():
    b0 = 0.2
    b1 = 1
    scheme_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/PGSE_21_dir_12_b.scheme"
    path_to_data = "/home/localadmin/Documents/CATERPillar/arthurs_analysis"
    df_pgse, df_waveform = calculate_ADC(b0 ,b1, scheme_path, path_to_data)
    plot_with_respect_to_direction(df_pgse, df_waveform)

main()