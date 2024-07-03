import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
from scipy.optimize import least_squares
from DKI import calculate_DKI
from useful_functions import extract_simulation_info
import copy

def read_DWI(file_path):

    DWI_string = np.fromfile(file_path, dtype="float32")
    DWI_array = []
    for t in DWI_string:
        DWI_array.append(t)
    return DWI_array
    

def read_scheme(file_path):

    with open(file_path, 'r') as file:
        content = file.read()

    content = content.split("\n")
    content = content[1:]
    content = [x.split() for x in content]

    xs = []
    ys = []
    zs = []
    Gs = []
    Deltas = []
    deltas = []
    TEs = []
    for line in content:
        if len(line) == 0:
            continue

        xs.append(line[0])
        ys.append(line[1])
        zs.append(line[2])
        Gs.append(line[3])
        Deltas.append(line[4])
        deltas.append(line[5])
        TEs.append(line[6])

    # dataset
    data = {'x': xs,
            'y': ys,
            'z': zs,
            'G': Gs,
            'Delta': Deltas,
            'delta': deltas,
            'TE': TEs}
    data = pd.DataFrame(data)

    return data

def residuals(params, b, ln_S_S0):
    D, K = params
    return ln_S_S0 - (-b * D + (b ** 2 * D ** 2 * K) / 6)



def create_df_parameters_1d(DWIs, data):

    data["DWI"] = DWIs

    # data as numerical values
    for column in data.columns:
        data[column] = pd.to_numeric(data[column])

    giro = 2.6751525e5 #rad/msXT

    # create DWI/DWI_0 column which is DWI divided by DWI value when G = 0 for each diffusion time
    data = data.loc[data["DWI"] > 0]
    
    data["DWI/DWI_0"] = data["DWI"]/data.loc[data["G"] == 0]["DWI"].values[0]
    data["ln(DWI/DWI_0)"] = [math.log(d) for d in data["DWI/DWI_0"]]

    data["diffusion_time (s)"] = data["Delta"] - data["delta"]/3
    data["diffusion_time (ms)"] = data["diffusion_time (s)"]*1000                             

    data["1/sqrt(t)"] = [1.0/math.sqrt(float(t)) for t in data["diffusion_time (ms)"]]
    data["b-value"] = (data["G"]*data["delta"]*giro)*(data["G"]*data["delta"]*giro) * data["diffusion_time (s)"] 
    data["b-value"] = np.round(data["b-value"], 2)
    data = data.loc[data["b-value"] < 5000]
    # diffusion times
    diffusion_times = data["diffusion_time (ms)"].unique()

    all_data = []
    for t in diffusion_times:
        data_t = data[data["diffusion_time (ms)"] == t]
        
        # Initial guesses for D and K
        initial_guesses = [1e-3, 1e-3]

        # Set bounds for the parameters: D and K
        # D can be any value, but K must be non-negative
        lower_bounds = [0, 0]
        upper_bounds = [np.inf, np.inf]

        # Perform the least squares optimization with bounds
        result = least_squares(residuals, initial_guesses, args=(np.array(list(data_t["b-value"])), np.array(list(data_t["ln(DWI/DWI_0)"]))), bounds=(lower_bounds, upper_bounds))

        # Extract the optimal values of D and K
        D_opt, K_opt = result.x

        # Add the optimal values to the data frame
        data_t["D [mm²/s]"] = [D_opt]*len(data_t)
        data_t["K"] = [K_opt]*len(data_t)

        all_data.append(data_t)

    data_final = pd.concat(all_data)

    data_final["D/D0"] = data_final["D [mm²/s]"]/2e-3
    data_final["D [um²/ms]"] = data_final["D [mm²/s]"]*1e3

    return data_final



def create_df_parameters_21d(DWIs, data, path_to_DWI, path_scheme):

    data["DWI"] = DWIs

    # data as numerical values
    for column in data.columns:
        data[column] = pd.to_numeric(data[column])

    giro = 2.6751525e5 #rad/msXT

    # create DWI/DWI_0 column which is DWI divided by DWI value when G = 0 for each diffusion time
    #data = data.loc[data["DWI"] > 0]
    
    #data["DWI/DWI_0"] = data["DWI"]/data.loc[data["G"] == 0]["DWI"].values[0]
    #data["ln(DWI/DWI_0)"] = [math.log(d) for d in data["DWI/DWI_0"]]

    data["diffusion_time (s)"] = data["Delta"] - data["delta"]/3                         

    data["b-value"] = (data["G"]*data["delta"]*giro)*(data["G"]*data["delta"]*giro) * (data["Delta"] - data["delta"]/3)
    data["b-value"] = np.round(data["b-value"], 2)

    bvals = np.array(data["b-value"])/1000

    bvecs = np.array(data[["x", "y", "z"]])

    FA, MD, AD, RD, MK, AK, RK = calculate_DKI(path_scheme, path_to_DWI, bvals, bvecs)

    data["FA"] = FA
    data["MD"] = MD
    data["AD"] = AD
    data["RD"] = RD
    data["MK"] = MK
    data["AK"] = AK
    data["RK"] = RK
    #data["axial_diffusion_orientation"] = [axial_diffusion_orientation]*len(data)


    return data

def main():
    file_path_scheme ="/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/PGSE_21_dir_12_b.scheme"
    data = read_scheme(file_path_scheme)
    #################################### BASIC TEST ####################################
    #folders = ["/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/"]

    #nbr_repetitions = 1

    #all_dfs = []

    #for folder in folders:
    #    for rep in range(nbr_repetitions):
    #        if rep == 0:
    #            info_path = f"{folder}/test_simulation_info.txt"
    #            DWI_path = f"{folder}/test_DWI.bfloat"
    #        else:
    #            info_path = f"{folder}/test_rep_0{rep-1}_simulation_info.txt"
    #            DWI_path = f"{folder}/test_rep_0{rep-1}_DWI.bfloat"
    #        N, steps = extract_simulation_info(info_path)
    #        DWIs = read_DWI(DWI_path)
    #        df = create_df_parameters_21d(DWIs, data, DWI_path, file_path_scheme)
        
    #        df["N"] = [N]*len(df)
    #        df["steps"] = [steps]*len(df)
    #        df["rep"] = [rep]*len(df)

    #        all_dfs.append(df)

    #df_final = pd.concat(all_dfs)
    #print(df_final)



    #################################### CONVERGENCE ANALYSIS NUMBER OF WALKERS ####################################

    #Ns = [10000, 50000, 100000]

    #all_dfs = []
    #for N in Ns:
    #    file_path = f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_walkers/N_{N}/_DWI.bfloat"
    #    file_path_info = f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_walkers/N_{N}/_simulation_info.txt"
    #    DWIs = read_DWI(file_path)
    #    df = create_df_parameters_21d(DWIs, data, file_path, file_path_scheme)
    #    all_dfs.append(df)

    #df_final = pd.concat(all_dfs)
    #print(df_final)


    #################################### CONVERGENCE ANALYSIS NUMBER OF STEPS ####################################

    #folders = ["/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.3", "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.25", "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.2", "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.15", "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.1", "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.05"]
    folders = ["/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.3", "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.25", "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.2", "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_steps/step_length_0.15"]

    nbr_repetitions = 5

    all_dfs = []

    for folder in folders:
        for rep in range(nbr_repetitions):
            if rep == 0:
                info_path = f"{folder}/_simulation_info.txt"
                DWI_path = f"{folder}/_DWI.bfloat"
            else:
                info_path = f"{folder}/_rep_0{rep-1}_simulation_info.txt"
                DWI_path = f"{folder}/_rep_0{rep-1}_DWI.bfloat"

            N, steps = extract_simulation_info(info_path)
            DWIs = read_DWI(DWI_path)
            df = create_df_parameters_21d(DWIs, data, DWI_path, file_path_scheme)
            
            df["N"] = N
            df["steps"] = steps
            df["rep"] = rep

            # Make a deep copy of the DataFrame before appending
            df_copy = copy.deepcopy(df)
            all_dfs.append(df_copy)

    df_final = pd.concat(all_dfs)
    print(df_final)
    sns.scatterplot(x="steps", y="MD", data=df_final, hue = "rep")
    plt.show()

    sns.scatterplot(x="steps", y="FA", data=df_final, hue = "rep")
    plt.show()

    sns.scatterplot(x="steps", y="MK", data=df_final, hue = "rep")
    plt.show()
    #################################### TIME DEPENDENCE ANALYSIS ####################################

    #i = Ns[-1]
    #file_path = f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence/N_{i}/_DWI.bfloat"
    #DWIs = read_DWI(file_path)
    #df = create_df_parameters(DWIs, data)
    #sns.lmplot(x="1/sqrt(t)", y="D/D0", data=df)
    #plt.show()
    #sns.lineplot(x="1/sqrt(t)", y="K", data=df)
    #plt.show()
