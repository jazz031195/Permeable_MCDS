import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
from scipy.optimize import least_squares
from DKI import calculate_DKI, calculate_DKI_WMTI, calculate_DTI
from useful_functions import extract_simulation_info, get_files_from_folder, extract_simulation_time, get_csv_files_from_folder, read_and_extract_parameters
import copy
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import os

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
    return ln_S_S0 - (-b * D + ((b ** 2) * (D ** 2) * K) / 6)

def residuals_(params, b, D, ln_S_S0):
    K = params
    return ln_S_S0 - (-b * D + ((b ** 2) * (D ** 2) * K) / 6)

def model_function(diffusion_time, A, C):
    return A + C * (diffusion_time ** (-0.5))


def diffusion_kurtosis_model(b, D, K):
    """
    Diffusion Kurtosis Imaging (DKI) model.
    
    Parameters:
        b (float or np.array): b-value(s).
        D (float): Diffusion coefficient (ADC).
        K (float): Kurtosis coefficient.
        
    Returns:
        S/S0 (np.array): The normalized signal for the given b-values.
    """
    return np.exp(-b * D + (b ** 2) * D ** 2 * K / 6)

def diffusion_kurtosis_model_(b, D, K):
    """
    Diffusion Kurtosis Imaging (DKI) model.
    
    Parameters:
        b (float or np.array): b-value(s).
        D (float): Diffusion coefficient (ADC).
        K (float): Kurtosis coefficient.
        
    Returns:
        S/S0 (np.array): The normalized signal for the given b-values.
    """
    return -b * D + (b ** 2) * D ** 2 * K / 6


def create_df_parameters_1d(DWIs, data, calculate_diffusion= True):


    data["DWI"] = DWIs


    # data as numerical values
    for column in data.columns:
        data[column] = pd.to_numeric(data[column])

    giro = 2.6751525e5 #rad/msXT

    # create DWI/DWI_0 column which is DWI divided by DWI value when G = 0 for each diffusion time
     #data = data.loc[data["DWI"] > 0]
    # find all different directions that are combinations of x, y and z
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


def create_2d_df_time_dependence(DWIs, data):

    data["DWI"] = DWIs


    # data as numerical values
    for column in data.columns:
        data[column] = pd.to_numeric(data[column])

    giro = 2.6751525e5 #rad/msXT

    # create DWI/DWI_0 column which is DWI divided by DWI value when G = 0 for each diffusion time
    data = data.loc[data["DWI"] > 0]
    data["diffusion_time (ms)"] = (data["Delta"] - data["delta"]/3)*1000

    # find all different directions that are combinations of x, y and z
    diffusion_times = data["diffusion_time (ms)"].drop_duplicates()
    all_data = []
    for time in diffusion_times:

        data_direction = data.loc[(data["diffusion_time (ms)"] == time)]
        data_direction["b_value"] = (data_direction["G"]*data_direction["delta"]*giro)*(data_direction["G"]*data_direction["delta"]*giro) * ((data["Delta"] - data["delta"]/3))
        data_direction["b_value"] = np.round(data_direction["b_value"], 2)
        data_direction["DWI/DWI_0"] = data_direction["DWI"]/data_direction.loc[data_direction["G"] == 0]["DWI"].values[0]

        #S0_row = data_direction[data_direction['b_value'] == 0]
        # There might be multiple measurements with b-value = 0; take the mean if so
        #S0 = S0_row['DWI'].mean()

        #diffusion_data = data_direction[data_direction['b_value'] == 1000]

        # Extract signals and gradient directions
        #S_i = diffusion_data['DWI'].values  # Shape (65,)
        #gradients = diffusion_data[['x', 'y', 'z']].values  # Shape (65, 3)

        # Normalize gradient directions to ensure they are unit vectors
        #gradients /= np.linalg.norm(gradients, axis=1)[:, np.newaxis]

        #b_value = 1000  # Since all b-values here are 1000 s/mm²

        #n = gradients.shape[0]  # Number of measurements (65)
        #X = np.zeros((n, 6))

        #g = gradients  # For brevity in code

        #X[:, 0] = -b_value * g[:, 0] ** 2          # Dxx component
        #X[:, 1] = -b_value * g[:, 1] ** 2          # Dyy component
        #X[:, 2] = -b_value * g[:, 2] ** 2          # Dzz component
        #X[:, 3] = -2 * b_value * g[:, 0] * g[:, 1]  # Dxy component
        #X[:, 4] = -2 * b_value * g[:, 0] * g[:, 2]  # Dxz component
        #X[:, 5] = -2 * b_value * g[:, 1] * g[:, 2]  # Dyz component

        # Ensure S_i and S0 are properly formatted
        #Y = np.log(S_i / S0)  # Shape (65,)

        # Solve the linear system X * beta = Y
        #beta, residuals, rank, s = np.linalg.lstsq(X, Y, rcond=None)
        #Dxx, Dyy, Dzz, Dxy, Dxz, Dyz = beta

        #D = np.array([[Dxx, Dxy, Dxz],
        #      [Dxy, Dyy, Dyz],
        #      [Dxz, Dyz, Dzz]])
        
        #eigvals, eigvecs = np.linalg.eigh(D)

        # Sort eigenvalues in descending order
        #idx = eigvals.argsort()[::-1]
        #eigvals = eigvals[idx]
        #eigvecs = eigvecs[:, idx]

        #lambda1, lambda2, lambda3 = eigvals

        #RD = (lambda2 + lambda3) / 2
        #AD = lambda1

        # fit DTI model
        FA, MD, AD, RD, MK, AK, RK =  calculate_DKI(np.array(data_direction["DWI/DWI_0"]), np.array(data_direction["b_value"]), np.array(data_direction[["x", "y", "z"]]))

        # Add the optimal values to the data frame
        data_direction["RD"] = [RD]*len(data_direction)
        data_direction["AD"] = [AD]*len(data_direction)
        data_direction["MD"] = [MD]*len(data_direction)
        data_direction["FA"] = [FA]*len(data_direction)
        data_direction["MK"] = [MK]*len(data_direction)
        data_direction["AK"] = [AK]*len(data_direction)
        data_direction["RK"] = [RK]*len(data_direction)



        all_data.append(data_direction)

    data_final = pd.concat(all_data)



    return data_final


def create_df_time_dependence(DWIs, data):

    data["DWI"] = DWIs


    # data as numerical values
    for column in data.columns:
        data[column] = pd.to_numeric(data[column])

    giro = 2.6751525e5 #rad/msXT

    # create DWI/DWI_0 column which is DWI divided by DWI value when G = 0 for each diffusion time
    data = data.loc[data["DWI"] > 0]
    data["diffusion_time (ms)"] = (data["Delta"] - data["delta"]/3)*1000
    # data = data.loc[data["diffusion_time (ms)"]< 100]

    data["1/sqrt(t)"] = 1/np.sqrt(data["diffusion_time (ms)"])

    # find all different directions that are combinations of x, y and z
    diffusion_times = data["diffusion_time (ms)"].drop_duplicates()
    all_data = []
    for time in diffusion_times:

        data_direction = data.loc[(data["diffusion_time (ms)"] == time)]
        data_direction["DWI/DWI_0"] = data_direction["DWI"]/data_direction.loc[data_direction["G"] == 0]["DWI"].values[0]
        data_direction["ln(DWI/DWI_0)"] = [math.log(d) if d > 0 else np.nan for d in data_direction["DWI/DWI_0"]]

        data_direction["b_value"] = (data_direction["G"]*data_direction["delta"]*giro)*(data_direction["G"]*data_direction["delta"]*giro) * ((data["Delta"] - data["delta"]/3))
        data_direction["b_value"] = np.round(data_direction["b_value"], 2)

        data_direction_ = data_direction.loc[data_direction["b_value"]< 2500].copy()

        initial_guesses = [0.0, 0.0]

         # Fit the DKI model to the data
        params, _ = curve_fit(diffusion_kurtosis_model, np.array(list(data_direction_["b_value"])), np.array(list(data_direction_["DWI/DWI_0"])), p0=initial_guesses, bounds=(-np.inf, np.inf))
        D_opt, K_opt = params

        # Add the optimal values to the data frame
        data_direction["D [mm²/s]"] = [D_opt]*len(data_direction)
        data_direction["K"] =[K_opt]*len(data_direction)  
        data_direction["fit"] = diffusion_kurtosis_model_(data_direction["b_value"], data_direction["D [mm²/s]"], data_direction["K"])

        all_data.append(data_direction)

    data_final = pd.concat(all_data)

    data_final["D/D0"] = data_final["D [mm²/s]"]/2e-3
    data_final["D [um²/ms]"] = data_final["D [mm²/s]"]*1e3


    return data_final


def create_df_parameters_21d(DWIs, data, path_to_DWI, calculate_DKI_ = True):

    data["DWI"] = DWIs

    # data as numerical values
    for column in data.columns:
        if column != "factor":
            data[column] = pd.to_numeric(data[column])

    giro = 2.6751525e5 #rad/msXT

    data["diffusion_time (s)"] = data["Delta"] - data["delta"]/3                         

    data["b_value"] = (data["G"]*data["delta"]*giro)*(data["G"]*data["delta"]*giro) * (data["Delta"] - data["delta"]/3)
    data["b_value"] = np.round(data["b_value"], 2)

    if calculate_DKI_:
        bvecs = np.array(data[["x", "y", "z"]])
        bvals = np.array(data["b_value"])

        FA, MD, AD, RD, MK, AK, RK = calculate_DKI(path_to_DWI, bvals, bvecs)

        data["b_value"] = data["b_value"]/1000

        data["FA"] = FA
        data["MD"] = MD
        data["AD"] = AD
        data["RD"] = RD
        data["MK"] = MK
        data["AK"] = AK
        data["RK"] = RK

    return data

def convert_steps_to_time(file_path_scheme, folder):

    files = get_files_from_folder(folder)
    all_dfs = []

    for file in files:
        if "img" in file or "info" in file or "swc" in file:
            continue
        print(file)
        data = read_scheme(file_path_scheme)
        DWIs = read_DWI(file)
        data = create_df_parameters_1d(DWIs, data)
        data["DWI"] = DWIs
        #bs, Gs, deltas, Deltas, vectors, TEs = read_and_extract_parameters(file_path_scheme)
        #FA, MD, AD, RD, MK, AK, RK = calculate_DKI(file, bs, vectors)
        #data["FA"] = FA
        #data["MD"] = MD
        #data["AD"] = AD
        #data["RD"] = RD
        #data["MK"] = MK
        #data["AK"] = AK
        #data["RK"] = RK



        title = file.split("_DWI")[0]
        info_path = f"{title}_simulation_info.txt"
        N, steps = extract_simulation_info(info_path)
        print(steps)

        data["N"] = N
        data["steps"] = steps

        # Make a deep copy of the DataFrame before appending
        df_copy = copy.deepcopy(data)
        all_dfs.append(df_copy)

    df_final = pd.concat(all_dfs)
    print(df_final)
    sns.boxplot(x="steps", y="D [um²/ms]", data=df_final, hue= "direction")
    plt.show()

    sns.boxplot(x="steps", y="K", data=df_final, hue= "direction")
    plt.show()


    

def calculate_difference(group):
    # Get the DWI value for the 'cylinder' factor
    cylinder_dwi = np.mean(group.loc[group['factor'] == 'cylinder', 'DWI/DWI0'].values)
    
    # Calculate the difference for each row in the group
    group['DWI_diff'] = group['DWI/DWI0'] - cylinder_dwi
    
    return group

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
            df = create_df_parameters_1d(DWIs, data_scheme, True)
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
            df = create_df_parameters_1d(DWIs, data_scheme, False)
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

            
def time_analysis():
    # folder2 = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/time_dependence_tortuosity"
    #folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/time_dependence_beading"
    folder2 = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/time_dependence_intra_f_0.1"
    #file_path_scheme = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/PGSE_time_dependance_high_b.scheme"
    file_path_scheme = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/time_dependence_.scheme"
    files = get_files_from_folder(folder2)
    #files.extend(get_files_from_folder(folder2))
    df_all = []

    for file in files:
        if "img" in file or "info" in file or "swc" in file or "cylinders" in file:
            continue

        data = read_scheme(file_path_scheme)
        DWIs = read_DWI(file)
        print(file)
        df = create_df_time_dependence(DWIs, data)
        title = file.split("_DWI")[0]
        title = title.split("/")[-1]
        if "rep" in file :
            df["type"] = title.split("_rep_")[0]
            if (title.split("_rep_")[0]) == "beading":
                df["type"] = "beading_0.3"
            df["rep"] = title.split("_rep_")[1][:2]
        else:
            df["type"] = title.split("_DWI")[0]
            if (title.split("_DWI")[0]) == "beading":
                df["type"] = "beading_0.3"
            df["rep"] = 0
        # Make a deep copy of the DataFrame before appending
        df_copy = copy.deepcopy(df)
 
        df_all.append(df_copy)

    df = pd.concat(df_all)
    df_ = df.copy()
    df_["diffusion_time (ms)"] = np.round(df_["diffusion_time (ms)"].astype(float), 3)

    df_["1/t"] = 1/df_["diffusion_time (ms)"]
    # df_ = df_.loc[(df_["diffusion_time (ms)"] <= 80) & (df_["diffusion_time (ms)"] >= 15)]

    # Set up the FacetGrid
    g = sns.FacetGrid(df_, col="type", hue="diffusion_time (ms)", sharey=True, palette="Set1")

    # Plot the original data points as circles
    g.map_dataframe(sns.scatterplot, x="b_value", y="ln(DWI/DWI_0)")

    # Overlay the fitted data as crosses, with the same hue
    g.map_dataframe(sns.lineplot, x="b_value", y="fit")

    # Add grid lines to each subplot
    for ax in g.axes.flat:
        ax.grid(True)

    # Add the legend for the hue
    g.add_legend(title="Type")

    # Display the plot
    plt.show()



    sns.lmplot(x="1/sqrt(t)", y="D [um²/ms]",hue = "type", data=df_)
    plt.show()


    # Define the model D = a / sqrt(t)
    def model_sqrt(t, c, d):
        return d + c / np.sqrt(t)

    # Define the model D = b / t
    def model_inverse(t, c, d):
        return d+ c / t
    
    def model_independent(t,d):
        return d + 0*t
    
    for type in df_["type"].unique():

        _df_ = df_.loc[df_["type"] == type]
        # Assuming your dataframe is named df
        t = _df_["diffusion_time (ms)"].values
        D = _df_["D [um²/ms]"].values
        # Fit model 1: D = a / sqrt(t)
        popt_sqrt, _ = curve_fit(model_sqrt, t, D)

        # Fit model 2: D = b / t
        popt_inverse, _ = curve_fit(model_inverse, t, D)

        popt_independent, _ = curve_fit(model_independent, t, D)

        # Predicted D values from both models
        D_pred_sqrt = model_sqrt(t, *popt_sqrt)
        D_pred_inverse = model_inverse(t, *popt_inverse)
        D_pred_independent = model_independent(t, *popt_independent)

        # Calculate R² for both models
        r2_sqrt = r2_score(D, D_pred_sqrt)
        r2_inverse = r2_score(D, D_pred_inverse)
        r2_independent = r2_score(D, D_pred_independent)

        plt.figure(figsize=(10, 6))
        print(type)
        if type == "beading_0.3":
            color = sns.color_palette()[0]
        elif type == "periodic_beading":
            color = sns.color_palette()[2]
        elif type == "straight":
            color = sns.color_palette()[3]
        elif type == "beading_tortuous":
            color = sns.color_palette()[1]
  
        # Plot the original data
        sns.scatterplot(data = _df_ , x="diffusion_time (ms)", y="D [um²/ms]", color = color, marker='x')


        # reorder D_pred values according to t
        D_pred_sqrt = D_pred_sqrt[np.argsort(t)]
        D_pred_inverse = D_pred_inverse[np.argsort(t)]
        # sort t
        t = np.sort(t)

        # Plot the fitted models
        plt.plot(t, D_pred_sqrt, label=f"1/sqrt(t) (p=0), R²={r2_sqrt:.4f}", color='grey')
        print(t)
        plt.plot(t, D_pred_inverse, label=f"1/t (p=1), R²={r2_inverse:.4f}", color='black')
        #plt.scatter(t, D_pred_independent, label=f"Const (p=inf), R²={r2_independent:.4f}", color='green', marker='x')
        plt.xlabel("Diffusion Time (ms)")
        plt.ylabel("D [um²/ms]")
        plt.legend()
        plt.title(f" {type}")
        plt.show()






    # sort type column by alphabetical order
    df_ = df_.sort_values(by='type')
    
    #  df_ = df_.loc[(df_["diffusion_time (ms)"] <= 80)] 
    

    for type in df_["type"].unique():

        _df_ = df_.loc[df_["type"] == type]
        # Extract the relevant columns
        diffusion_time = _df_['diffusion_time (ms)'].values
        D_values = _df_['D [um²/ms]'].values

        # Initial guess for the parameters [A, C, theta]
        initial_guess = [2.0, 1.0]

        # Perform the curve fitting
        popt, pcov = curve_fit(model_function, diffusion_time, D_values, p0=initial_guess)

        # popt contains the best-fitting parameters: [A, C, theta]
        A, C = popt
        print("Type: ", type)
        print(f"Fitted parameters:")
        print(f"D_inf = {A}")
        print(f"c = {C}")



  
def time_analysis_extra():

    folder2 = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/time_dependence_extra"
    file_path_scheme = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/DTI_multi_shell.scheme"
    files = get_files_from_folder(folder2)

    
    csv_path = f"{folder2}/time_dependence_extra.csv"

    if os.path.exists(csv_path):
        df_ = pd.read_csv(csv_path)
    else:
        df_all = []

        for file in files:
            if "img" in file or "info" in file or "swc" in file or "cylinders" in file:
                continue
 
            data = read_scheme(file_path_scheme)
            DWIs = read_DWI(file)
            print(file)
            df = create_2d_df_time_dependence(DWIs, data)
            title = file.split("_DWI")[0]
            title = title.split("/")[-1]
            if "rep" in file :
                df["type"] = title.split("_rep_")[0]
                if (title.split("_rep_")[0]) == "beading":
                    df["type"] = "beading_0.3"
                df["rep"] = title.split("_rep_")[1][:2]
            else:
                df["type"] = title.split("_DWI")[0]
                if (title.split("_DWI")[0]) == "beading":
                    df["type"] = "beading_0.3"
                df["rep"] = 0
            # Make a deep copy of the DataFrame before appending
            df_copy = copy.deepcopy(df)
    
            df_all.append(df_copy)

        df = pd.concat(df_all)
        df_ = df.copy()
        # save the data
        df_.to_csv(csv_path, index=False)

    def fit_ln(t, A, B):
        """
        The model function: D(t) = A + (B * ln(t / C)) / t
        """
        return A + (B * np.log(t)) / t
    
    def fit_1_t(t, A, B):
        """
        The model function: D(t) = A + (B * ln(t / C)) / t
        """
        return A + (B  / t)

    
    for type in df_["type"].unique():
        _df_ = df_.loc[df_["type"] == type]
        # Assuming your dataframe is named df
        t = _df_["diffusion_time (ms)"].values
        D = _df_["RD"].values
        # Fit model 1/ln(t)
        popt_1, _ = curve_fit(fit_ln, t, D)

        # Fit model 1/t
        popt_2, _ = curve_fit(fit_1_t, t, D)

        # Predicted D values from both models
        D_pred_1 = fit_ln(t, *popt_1)
        D_pred_2 = fit_1_t(t, *popt_2)

        # Calculate R² for both models
        r2_1 = r2_score(D, D_pred_1)
        r2_2 = r2_score(D, D_pred_2)

        plt.figure(figsize=(10, 6))

        if type == "beading_0.3":
            color = sns.color_palette()[0]
        elif type == "straight":
            color = sns.color_palette()[3]
        elif type == "beading_tortuous":
            color = sns.color_palette()[1]
  
        # Plot the original data
        sns.scatterplot(data = _df_ , x="diffusion_time (ms)", y="RD", color = color, marker='x')

        # reorder D_pred values according to t
        D_pred_1 = D_pred_1[np.argsort(t)]
        D_pred_2 = D_pred_2[np.argsort(t)]
        # sort t
        t = np.sort(t)

        # Plot the fitted models
        plt.plot(t, D_pred_1, label=f"1/sqrt(t) (p=0), R²={r2_1:.4f}", color='grey')
        plt.plot(t, D_pred_2, label=f"1/t (p=1), R²={r2_2:.4f}", color='black')
        #plt.scatter(t, D_pred_independent, label=f"Const (p=inf), R²={r2_independent:.4f}", color='green', marker='x')
        plt.xlabel("Diffusion Time (ms)")
        plt.ylabel("D [um²/ms]")
        plt.legend()
        plt.title(f" {type}")
        plt.show()

    df_["diffusion_time (ms)"] = np.round(df_["diffusion_time (ms)"].astype(float), 3)
    
    palette = [sns.color_palette()[3], sns.color_palette()[0], sns.color_palette()[1]]
    
    sns.lineplot(x="diffusion_time (ms)", y="RD",hue = "type", data=df_, hue_order=["straight", "beading_0.3", "beading_tortuous"], palette=palette)
    plt.show()


if __name__ == "__main__":
    file_path_scheme ="/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/DTI.scheme"
    path_to_folder = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/convergence_timesteps"
    #convert_steps_to_time(file_path_scheme, path_to_folder)
    time_analysis_extra()
    #time_analysis()
    #converge_overlapping_factor(path_to_folder, file_path_scheme)
