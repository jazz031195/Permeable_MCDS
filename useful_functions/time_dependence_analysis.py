import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
from DKI import calculate_DKI
from useful_functions import get_files_from_folder, diffusion_kurtosis_model, read_DWI, read_scheme
import copy
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import os

def residuals(params, b, ln_S_S0):
    D, K = params
    return ln_S_S0 - (-b * D + ((b ** 2) * (D ** 2) * K) / 6)

def residuals_(params, b, D, ln_S_S0):
    K = params
    return ln_S_S0 - (-b * D + ((b ** 2) * (D ** 2) * K) / 6)

def model_function(diffusion_time, A, C):
    return A + C * (diffusion_time ** (-0.5))


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



def create_df_extra(DWIs, data):

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


def create_df_intra(DWIs, data):

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
        df = create_df_intra(DWIs, data)
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
            df = create_df_extra(DWIs, data)
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
    #time_analysis_extra()
    time_analysis()