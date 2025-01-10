import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
import copy
from scipy.optimize import curve_fit
from useful_functions import read_DWI, read_scheme, get_files_from_folder, extract_simulation_time, diffusion_kurtosis_model

giro = 2.6751525e5 #rad/msXT

def calculate_difference(group):
    # Get the DWI value for the 'cylinder' factor
    cylinder_dwi = np.mean(group.loc[group['factor'] == 'cylinder', 'DWI/DWI0'].values)
    
    # Calculate the difference for each row in the group
    group['DWI_diff'] = group['DWI/DWI0'] - cylinder_dwi
    
    return group


def create_df(DWIs, data, calculate_diffusion= True):

    data["DWI"] = DWIs

    # data as numerical values
    for column in data.columns:
        data[column] = pd.to_numeric(data[column])

    directions = data[["x", "y", "z"]].drop_duplicates()
    all_data = []
    for direction in directions.iterrows():
        x = direction[1]["x"]
        y = direction[1]["y"]
        z = direction[1]["z"]
        
        data_direction = data.loc[(data["x"] == x) & (data["y"] == y) & (data["z"] == z)]
        data_direction["direction"] = [f"{x}_{y}_{z}"]*len(data_direction)
        data_direction["DWI/DWI_0"] = data_direction["DWI"]/data_direction.loc[data_direction["G"] == 0]["DWI"].values[0]
        data_direction["diffusion_time (s)"] = data_direction["Delta"] - data_direction["delta"]/3
        data_direction["b_value"] = (data_direction["G"]*data_direction["delta"]*giro)*(data_direction["G"]*data_direction["delta"]*giro) * data_direction["diffusion_time (s)"] 
        data_direction["b_value"] = np.round(data_direction["b_value"], 2)

        if (calculate_diffusion):

            data_direction = data_direction.loc[data_direction["b_value"]<= 3001]
    
            # Initial guesses for D and K
            initial_guesses = [1e-3, 1e-3]

            # Fit the DKI model to the data
            params, _ = curve_fit(diffusion_kurtosis_model, np.array(list(data_direction["b_value"])), np.array(list(data_direction["DWI/DWI_0"])), p0=initial_guesses, bounds=(0, np.inf))

            # Extract the optimal values of D and K
            D_opt, K_opt = params

            # Add the optimal values to the data frame
            data_direction["D [mm²/s]"] = [D_opt]*len(data_direction)
            data_direction["K"] = [K_opt]*len(data_direction)
            #data_direction["D/D0"] = data_direction["D [mm²/s]"]/2e-3
            data_direction["D [um²/ms]"] = data_direction["D [mm²/s]"]*1e3


        all_data.append(data_direction)

    data_final = pd.concat(all_data)

    return data_final


def converge_overlapping_factor(path_to_folder, file_path_scheme):

    files = get_files_from_folder(path_to_folder)
    data_scheme = read_scheme(file_path_scheme)
    
    ############################## duffion parameters ########################################
    all_dfs = []
    for file in files:
        if "img" not in file and "info" not in file and "swc" not in file:
            print(file)

            if "factor_1_" in file:
                factor = "1"
            elif "factor_2_" in file:
                factor = "2"
            elif "factor_4_" in file:
                factor = "4"
            elif "factor_8_" in file:
                factor = "8"
            elif "factor_16_" in file:
                factor = "16"
            elif "cylinder_" in file:
                factor = "cylinder"
            else:
                print("Error, no factor found")
                continue
            
            if "intra" in file:
                location= True
            elif "extra" in file:
                location= False
            else:
                print("Error, no location found")
                continue
            
            DWIs = read_DWI(file)
            df = create_df(DWIs, data_scheme, True)
            df["factor"] = [factor]*len(df)
            df["isintra"] = [location]*len(df)

            df_copy = copy.deepcopy(df)
            all_dfs.append(df_copy)
        
        final_df = pd.concat(all_dfs)
    print(final_df)
    final_df = final_df.loc[final_df["direction"] == "0_0_1"]
    final_df["factor"] = pd.Categorical(final_df["factor"], categories=["1","2", "4", "8",  "cylinder"])
    df_intra = final_df.loc[final_df["isintra"] == False]
    del df_intra["direction"]
    sns.boxplot(x="factor", y="D [um²/ms]",hue = "factor", data=df_intra, order=["1","2","4","8","cylinder"] )
    plt.show()
    sns.boxplot(x="factor", y="K",hue = "factor", data=df_intra, order=["1","2","4","8","cylinder"] )
    plt.show()
    ############################## DWI comparison ########################################
    all_dfs = []
    for file in files:
        if "img" not in file and "info" not in file and "swc" not in file:
            print(file)

            if "factor_1_" in file:
                factor = "1"
            elif "factor_2_" in file:
                factor = "2"
            elif "factor_4_" in file:
                factor = "4"
            elif "factor_8_" in file:
                factor = "8"
            elif "factor_16_" in file:
                factor = "16"
            elif "cylinder_" in file:
                factor = "cylinder"
            else:
                print("Error, no factor found")
                continue
            
            if "intra" in file:
                location= True
            elif "extra" in file:
                location= False
            else:
                print("Error, no location found")
                continue
            
            DWIs = read_DWI(file)
            df = create_df(DWIs, data_scheme, False)
            df["factor"] = [factor]*len(df)
            df["isintra"] = [location]*len(df)

            df_copy = copy.deepcopy(df)
            all_dfs.append(df_copy)
        
        final_df = pd.concat(all_dfs)

    final_df = final_df.loc[final_df["direction"] == "0_0_1"]
    final_df["factor"] = pd.Categorical(final_df["factor"], categories=["1","2", "4", "8",  "cylinder"])
    df_intra = final_df.loc[final_df["isintra"] == False]
    del df_intra["direction"]
    df_intra["DWI/DWI0"] = df_intra["DWI"]/1e5
    df_intra =  df_intra.groupby('b_value').apply(calculate_difference).reset_index(drop=True)
    df_intra["MSE"] = df_intra["DWI_diff"]**2
    print(df_intra)
    sns.boxplot(x="factor", y="MSE",hue = "factor", data=df_intra, order=["1","2","4","8", "cylinder"] )
    plt.show()
    sns.lineplot(x="b_value", y="DWI/DWI0",hue = "factor", data=df_intra )
    #b_values = df_intra["b_value"].unique()
    #signal = []
    #for b in b_values:
    #    D = 2e-3
    #    S =  np.exp(-b * D)
    #    signal.append(S)
    #sns.lineplot(x=b_values, y=signal, color="black")
    plt.show()
    ############################## time for run ########################################


    df_times = []
    files = get_files_from_folder(path_to_folder, binary = False)
  
    for file in files:
        if "simulation_info" in file:
            nbr_simulations, time = extract_simulation_time(file)
            if "factor_1_" in file:
                factor = "1"
            elif "factor_2_" in file:
                factor = "2"
            elif "factor_4_" in file:
                factor = "4"
            elif "factor_8_" in file:
                factor = "8"
            elif "factor_16_" in file:
                factor = "16"
            elif "cylinder_" in file:
                factor = "cylinder"
            else:
                print("Error, no factor found")
                continue
            
            if "intra" in file:
                location= True
            elif "extra" in file:
                location= False
            else:
                print("Error, no location found")
                continue
            df_time = pd.DataFrame()
            df_time["time"] = [float(time)]
            df_time["nbr_simulations"] = [float(nbr_simulations)]
            df_time["factor"] = [factor]
            df_time["isintra"] = [location]
            df_copy = copy.deepcopy(df_time)
            df_times.append(df_copy)
    df_time = pd.concat(df_times)
    df_time = df_time.dropna()
    df_time["time (48 cores)"] = df_time["time"]*df_time["nbr_simulations"]/48
    df_time["log_time"] = [math.log(t) for t in df_time["time (48 cores)"]]
    df_time = df_time.loc[df_time["isintra"] == True]
    print(df_time.groupby("factor").mean())
    sns.scatterplot(x="factor", y="log_time", hue = "factor", data=df_time, hue_order=["1","2", "4", "8","cylinder"])
    plt.ylabel("Log(time (s))")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    file_path_scheme ="/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/PGSE_overlapping_spheres_analysis.scheme"
    path_to_folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/overlapping_factor"
    converge_overlapping_factor(path_to_folder, file_path_scheme)
    