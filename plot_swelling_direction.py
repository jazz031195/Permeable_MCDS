import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

from useful_functions import read_binary_file, get_files_from_folder
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)

def get_scheme_info(scheme_path):
    scheme = pd.read_csv(scheme_path, sep=" ", header=None, skiprows=1)
    scheme = scheme.dropna(axis=1)
    scheme = scheme.to_numpy()
    gradient_strength = scheme[:, 3]
    Delta, delta, TE = scheme[:, 4], scheme[:, 5], scheme[:, 6]
    b_values = ((giro*delta*gradient_strength)**2)*(Delta-delta/3)/1e9
    directions = scheme[:, :3]
    b_values = [round(float(i),2) for i in b_values]
    return b_values, directions

def calculate_ADC(b0 ,b1, scheme_path, path_to_data):
    b_values, directions = get_scheme_info(scheme_path)
    all_b_values = [] 
    all_directions = []
    all_DWIs = []
    all_swellings = []
    all_locations = []
    files = get_files_from_folder(path_to_data)
    for file in files:
        if "img" not in file and "simulation" not in file:
            DWI = read_binary_file(file) 
            DWI =[float(i) for i in DWI] 
            all_DWIs.extend(DWI)
            all_b_values.extend(b_values)
            all_directions.extend(directions)
            if "intra" in file:
                all_locations.extend(["intra"]*len(DWI))
            elif "extra" in file:
                all_locations.extend(["extra"]*len(DWI))
            else:
                print("Error, no location found")
                assert(0)
            if "swell_0_" in file:
                all_swellings.extend(np.zeros(len(DWI)))
            elif "swell_0.01_" in file:
                all_swellings.extend(np.ones(len(DWI))*0.01)
            elif "swell_0.005_" in file:
                all_swellings.extend(np.ones(len(DWI))*0.005)
            elif "swell_0.0025_" in file:
                all_swellings.extend(np.ones(len(DWI))*0.0025)
            elif "swell_0.0075_" in file:
                all_swellings.extend(np.ones(len(DWI))*0.0075)
            else:
                print("Error, no swelling found")
                assert(0)
    df = pd.DataFrame({"DWI": all_DWIs, "b_value": all_b_values, "direction": all_directions, "swelling": all_swellings, "location": all_locations})
    print(df)
    # calculate angle between direction and vector (0,0,1)
    df["angle (rad)"] = df["direction"].apply(lambda x: np.arccos(np.dot(x, [0,0,1]))) 
    # angle must be from 0 to pi
    df["angle (rad)"] = df["angle (rad)"].apply(lambda x: x if x <= np.pi/2 else np.pi-x)
    

    df_b0 = df.loc[df["b_value"] == b0]
    df_b0 = df_b0.groupby(["angle (rad)", "swelling", "location"]).sum()
    df_b0 = df_b0.groupby(["angle (rad)", "swelling"]).sum()
    print(df_b0)
    df_b1 = df.loc[df["b_value"] == b1]
    df_b1 = df_b1.groupby(["angle (rad)", "swelling", "location"]).sum()
    df_b1 = df_b1.groupby(["angle (rad)", "swelling"]).sum()
    print(df_b1)
    
    # calculate adc list form the two b values
    adc = np.log(df_b0["DWI"].values/df_b1["DWI"].values)/(b1-b0)
    df_b0["ADC"] = adc
    df_b0["ADC"] = df_b0["ADC"].replace([np.inf, -np.inf], np.nan)
    df_b0 = df_b0.dropna()

    
    print(df_b0)

    return df_b0.reset_index()

def plot_with_respect_to_direction(df):
    # df = df.loc[df.location == "intra"] 
     #df = df.loc[df.location == "extra"] 
    # df = df.loc[df.angle > 0.4*np.pi] 

    # 3d plot, x = angle between direction and (0,0,1), y = swelling, z = ADC
    # surface plotting with the mean ADC for each angle and swelling
        
    # Extract axes data from DataFrame
    x = df['angle (rad)']
    y = df['swelling']
    z = df['ADC']

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot surface
    ax.plot_trisurf(x, y, z, cmap='viridis')

    # Set labels and title
    ax.set_xlabel('Angle difference (rad)')
    ax.set_ylabel('Percentage of swelling')
    ax.set_zlabel('ADC (um2/ms)')
    ax.set_title('ADC decrease with swelling')

    # Show plot
    plt.show()
    # multiplie swelling by 100 to get percentage 
    df["swelling"] = df["swelling"]*100
    print(df)
    # normalise ADC with respect to adc at swelling = 0 for each direction
    df["ADC"] = df.groupby("angle (rad)")["ADC"].apply(lambda x: x/x.iloc[0]) 
 
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.5)
    sns.lineplot(hue="angle (rad)", y="ADC", x="swelling", data=df)
    plt.ylabel('Relative ADC (um2/ms)')
    plt.xlabel('Swelling (%)')
    plt.title('Relative ADC decrease with swelling')
    plt.show()


def main():
    b0 = 0.2
    b1 = 1
    scheme_path = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/docs/scheme_files/PGSE_sample_scheme_21_dir.scheme"
    path_to_data = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/instructions/axons/icvf_70_vox_50/"
    df = calculate_ADC(b0 ,b1, scheme_path, path_to_data)
    plot_with_respect_to_direction(df)

main()