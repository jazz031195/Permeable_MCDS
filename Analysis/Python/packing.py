"""
    File to compare the mesh with the analytical solutions for different
    decimation. Plots the mean diffusivity (MD) for different decimations
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
if sys.platform == "linux":
    sys.path.insert(1, '/home/localadmin/Documents/analytical_formula/')
else:
    sys.path.insert(1, '/Users/ideriedm/Documents/analytical_formula/')
import pandas as pd
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
from utils import get_bvals, get_bvectors, calculate_DKI, get_dwi, get_psge, create_data, analytical_solutions
scheme_file_path = "docs/scheme_files/PGSE_21_dir_12_b.scheme"
folder_intra = Path("results/ISMRM24/packing/overlap2/ICVF_0_63/intra")
df_all_intra = pd.DataFrame()
for filename in os.listdir(folder_intra):
    # Check if the filename contains "_rep_" and "DWI"
    if "DWI_img" in filename:
        # Name of the experience
        name         = ('_').join(filename.split('_')[:-1])
        # Number of walkers
        N            = 630000
        # Number of timesteps
        T            = 50000
        extension    = filename.split('_')[-1].split('.')[-1]
        SNR          = np.inf
        data_one_exp = create_data(folder_intra, SNR, name, extension, scheme_file_path)
        FA = data_one_exp["FA"][0]
        MD = data_one_exp["MD"][0]
        AD = data_one_exp["AD"][0]
        RD = data_one_exp["RD"][0]
        MK = data_one_exp["MK"][0]
        AK = data_one_exp["AK"][0]
        RK = data_one_exp["RK"][0]
        # For each b, iterate over all directions, store the data, and average them (powder-average)
        nb_b   = len(data_one_exp["b [ms/um²]"].unique())
        nb_dir = int(len(data_one_exp["x"].values) / nb_b)
        for i in range(nb_b):
            sb_so = []
            adc   = []
            for j in range(nb_dir):
                sb_so.append(data_one_exp.iloc[nb_b * j + i, :]["Sb/So"])
                adc.append(data_one_exp.iloc[nb_b * j + i, :]["adc [ms/um²]"])
                bval = data_one_exp.iloc[nb_b * j + i, :]["b [ms/um²]"]
            # Powder-average signal
            mean     = np.mean(sb_so)
            # Powder-average ADC
            mean_adc = np.mean(adc)
            if i == 1:
                d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc,
                    'b [ms/um²]': bval, 'overlap': 2,
                    'FA': FA, 'MD': MD, 'RD': RD, 'AD': AD, 'MK': MK, 'RK': RK, 'AK': AK}
            else:
                d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc,
                    'b [ms/um²]': bval, 'overlap': 2}
            df_avg_data = pd.DataFrame(d, index=[i])
            df_all_intra = pd.concat([df_all_intra, df_avg_data])

folder_extra = Path("results/ISMRM24/packing/overlap2/ICVF_0_63/extra")
df_all_extra = pd.DataFrame()
for filename in os.listdir(folder_extra):
    # Check if the filename contains "_rep_" and "DWI"
    if "DWI_img" in filename:
        # Name of the experience
        name         = ('_').join(filename.split('_')[:-1])
        # Number of walkers
        N            = 370000
        # Number of timesteps
        T            = 50000
        extension    = filename.split('_')[-1].split('.')[-1]
        SNR          = np.inf
        data_one_exp = create_data(folder_extra, SNR, name, extension, scheme_file_path)
        FA = data_one_exp["FA"][0]
        MD = data_one_exp["MD"][0]
        AD = data_one_exp["AD"][0]
        RD = data_one_exp["RD"][0]
        MK = data_one_exp["MK"][0]
        AK = data_one_exp["AK"][0]
        RK = data_one_exp["RK"][0]
        # For each b, iterate over all directions, store the data, and average them (powder-average)
        nb_b   = len(data_one_exp["b [ms/um²]"].unique())
        nb_dir = int(len(data_one_exp["x"].values) / nb_b)
        for i in range(nb_b):
            sb_so = []
            adc   = []
            for j in range(nb_dir):
                sb_so.append(data_one_exp.iloc[nb_b * j + i, :]["Sb/So"])
                adc.append(data_one_exp.iloc[nb_b * j + i, :]["adc [ms/um²]"])
                bval = data_one_exp.iloc[nb_b * j + i, :]["b [ms/um²]"]
            # Powder-average signal
            mean     = np.mean(sb_so)
            # Powder-average ADC
            mean_adc = np.mean(adc)
            if i == 1:
                d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc,
                    'b [ms/um²]': bval, 'overlap': 2,
                    'FA': FA, 'MD': MD, 'RD': RD, 'AD': AD, 'MK': MK, 'RK': RK, 'AK': AK}
            else:
                d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc,
                    'b [ms/um²]': bval, 'overlap': 2}
            df_avg_data = pd.DataFrame(d, index=[i])
            df_all_extra = pd.concat([df_all_extra, df_avg_data])

# print(df_all_extra[['adc [ms/um²]', 'FA', 'MD', 'AD', 'RD', 'b [ms/um²]']])

print("MD : ", df_all_intra.MD.mean(skipna=True)*0.63 + 0.37*df_all_extra.MD.mean(skipna=True))
print("FA : ", df_all_intra.FA.mean(skipna=True)*0.63 + 0.37*df_all_extra.FA.mean(skipna=True))
print("AD : ", df_all_intra.AD.mean(skipna=True)*0.63 + 0.37*df_all_extra.AD.mean(skipna=True))
print("RD : ", df_all_intra.RD.mean(skipna=True)*0.63 + 0.37*df_all_extra.RD.mean(skipna=True))


folder_intra = Path("results/ISMRM24/packing/overlap2/ICVF_0_63/intra_shrink")
df_all_intra = pd.DataFrame()
for filename in os.listdir(folder_intra):
    # Check if the filename contains "_rep_" and "DWI"
    if "DWI_img" in filename:
        # Name of the experience
        name         = ('_').join(filename.split('_')[:-1])
        # Number of walkers
        N            = 630000
        # Number of timesteps
        T            = 50000
        extension    = filename.split('_')[-1].split('.')[-1]
        SNR          = np.inf
        data_one_exp = create_data(folder_intra, SNR, name, extension, scheme_file_path)
        FA = data_one_exp["FA"][0]
        MD = data_one_exp["MD"][0]
        AD = data_one_exp["AD"][0]
        RD = data_one_exp["RD"][0]
        MK = data_one_exp["MK"][0]
        AK = data_one_exp["AK"][0]
        RK = data_one_exp["RK"][0]
        # For each b, iterate over all directions, store the data, and average them (powder-average)
        nb_b   = len(data_one_exp["b [ms/um²]"].unique())
        nb_dir = int(len(data_one_exp["x"].values) / nb_b)
        for i in range(nb_b):
            sb_so = []
            adc   = []
            for j in range(nb_dir):
                sb_so.append(data_one_exp.iloc[nb_b * j + i, :]["Sb/So"])
                adc.append(data_one_exp.iloc[nb_b * j + i, :]["adc [ms/um²]"])
                bval = data_one_exp.iloc[nb_b * j + i, :]["b [ms/um²]"]
            # Powder-average signal
            mean     = np.mean(sb_so)
            # Powder-average ADC
            mean_adc = np.mean(adc)
            if i == 1:
                d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc,
                    'b [ms/um²]': bval, 'overlap': 2,
                    'FA': FA, 'MD': MD, 'RD': RD, 'AD': AD, 'MK': MK, 'RK': RK, 'AK': AK}
            else:
                d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc,
                    'b [ms/um²]': bval, 'overlap': 2}
            df_avg_data = pd.DataFrame(d, index=[i])
            df_all_intra = pd.concat([df_all_intra, df_avg_data])

folder_extra = Path("results/ISMRM24/packing/overlap2/ICVF_0_63/extra_shrink")
df_all_extra = pd.DataFrame()
for filename in os.listdir(folder_extra):
    # Check if the filename contains "_rep_" and "DWI"
    if "DWI_img" in filename:
        # Name of the experience
        name         = ('_').join(filename.split('_')[:-1])
        # Number of walkers
        N            = 370000
        # Number of timesteps
        T            = 50000
        extension    = filename.split('_')[-1].split('.')[-1]
        SNR          = np.inf
        data_one_exp = create_data(folder_extra, SNR, name, extension, scheme_file_path)
        FA = data_one_exp["FA"][0]
        MD = data_one_exp["MD"][0]
        AD = data_one_exp["AD"][0]
        RD = data_one_exp["RD"][0]
        MK = data_one_exp["MK"][0]
        AK = data_one_exp["AK"][0]
        RK = data_one_exp["RK"][0]
        # For each b, iterate over all directions, store the data, and average them (powder-average)
        nb_b   = len(data_one_exp["b [ms/um²]"].unique())
        nb_dir = int(len(data_one_exp["x"].values) / nb_b)
        for i in range(nb_b):
            sb_so = []
            adc   = []
            for j in range(nb_dir):
                sb_so.append(data_one_exp.iloc[nb_b * j + i, :]["Sb/So"])
                adc.append(data_one_exp.iloc[nb_b * j + i, :]["adc [ms/um²]"])
                bval = data_one_exp.iloc[nb_b * j + i, :]["b [ms/um²]"]
            # Powder-average signal
            mean     = np.mean(sb_so)
            # Powder-average ADC
            mean_adc = np.mean(adc)
            if i == 1:
                d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc,
                    'b [ms/um²]': bval, 'overlap': 2,
                    'FA': FA, 'MD': MD, 'RD': RD, 'AD': AD, 'MK': MK, 'RK': RK, 'AK': AK}
            else:
                d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc,
                    'b [ms/um²]': bval, 'overlap': 2}
            df_avg_data = pd.DataFrame(d, index=[i])
            df_all_extra = pd.concat([df_all_extra, df_avg_data])

# print(df_all_extra[['adc [ms/um²]', 'FA', 'MD', 'AD', 'RD', 'b [ms/um²]']])

print("MD shrink : ", df_all_intra.MD.mean(skipna=True)*0.63 + 0.37*df_all_extra.MD.mean(skipna=True))
print("FA shrink : ", df_all_intra.FA.mean(skipna=True)*0.63 + 0.37*df_all_extra.FA.mean(skipna=True))
print("AD shrink : ", df_all_intra.AD.mean(skipna=True)*0.63 + 0.37*df_all_extra.AD.mean(skipna=True))
print("RD shrink : ", df_all_intra.RD.mean(skipna=True)*0.63 + 0.37*df_all_extra.RD.mean(skipna=True))