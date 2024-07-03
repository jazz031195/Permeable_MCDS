import numpy as np
import pandas as pd
from useful_functions import extract_simulation_info
from convergence import create_df_parameters_21d, read_scheme, read_DWI
import matplotlib.pyplot as plt


def signal_attenuation_in_cylinder(b_values, S0, D, R):
    """
    Calculate signal attenuation in a cylinder for given b-values.

    Parameters:
    - b_values: array-like, list or numpy array of b-values (s/mm^2)
    - S0: float, signal intensity without diffusion weighting (b = 0)
    - D: float, diffusion coefficient of water molecules (in mm^2/s)
    - R: float, radius of the cylinder (in mm)

    Returns:
    - S: numpy array, signal intensities for each b-value
    """
    lambda_n = 2.4048  # dimensionless constant for a cylinder
    ADC = D * lambda_n / (R**2)
    S = S0 * np.exp(-b_values * ADC)
    return S



def read_swc_file(file_path):
    columns = ["ax_id","sph_id", "branch_id", "type", "x", "y", "z", "Rin","Rout", "P"]
    df = pd.read_csv(file_path, sep=' ', names=columns)

    df = df.iloc[1:]
    df["x"] = [float(i) for i in list(df["x"])]
    df["y"] = [float(i) for i in list(df["y"])]
    df["z"] = [float(i) for i in list(df["z"])]
    df["Rout"] = [float(i) for i in list(df["Rout"])]
    df["Rin"] = [float(i) for i in list(df["Rin"])]
    df["ax_id"] = [float(i) for i in list(df["ax_id"])]
    df["sph_id"] = [float(i) for i in list(df["sph_id"])]
    df["P"] = [float(i) for i in list(df["P"])]
    df["branch_id"] = [float(i) for i in list(df["branch_id"])]
    return df

def get_cylinders_R(swc_file_path):
    df = read_swc_file(swc_file_path)
    df = df.loc[df["P"] == -1]
    R = df["Rout"].values
    return R

def signal_attenuation_all_cylinders(swc_file_path, b_values, S0):
    Rs = get_cylinders_R(swc_file_path)
    Rs = [R*1e-3 for R in Rs] # mm
    D = 2.0e-3   # Example diffusion coefficient in mm^2/s
    signals = []
    for R in Rs:
        s = signal_attenuation_in_cylinder(b_values, S0, D, R)
        signals.append(s)
    signal = np.mean(signals, axis=0)

    return signal


# Example usage:
if __name__ == "__main__":

    file_path_scheme ="/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/test.scheme"
    data = read_scheme(file_path_scheme)
    info_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/analytical_solution/cylinder__simulation_info.txt"
    N, steps = extract_simulation_info(info_path)
    DWI_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/analytical_solution/cylinder__DWI.bfloat"
    DWIs = read_DWI(DWI_path)
    df_cylinder = create_df_parameters_21d(DWIs, data, DWI_path, file_path_scheme) 

    info_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/analytical_solution/_simulation_info.txt"
    N, steps = extract_simulation_info(info_path)
    DWI_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/analytical_solution/_DWI.bfloat"
    DWIs = read_DWI(DWI_path)
    df_spheres = create_df_parameters_21d(DWIs, data, DWI_path, file_path_scheme) 


    b_values = np.array(df_cylinder["b-value"].unique())
    S0 = N
    signal = signal_attenuation_all_cylinders("/home/localadmin/Documents/MCDS/Permeable_MCDS/output/analytical_solution/growth_vox_200_factor_2_0.swc", b_values, S0)
    plt.plot(b_values, signal, label='analytical')
    plt.plot(df_cylinder["b-value"], df_cylinder["DWI"], 'ro', label='simulation cylinder')
    plt.plot(df_spheres["b-value"], df_spheres["DWI"], 'bo', label='simulation spheres')
    plt.xlabel('b-value (s/mm^2)')
    plt.ylabel('Signal Intensity')
    plt.title('DWI Signal Attenuation in a Cylinder')
    plt.grid(True)
    plt.legend()
    plt.show()



