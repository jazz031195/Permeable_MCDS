import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
import dipy.reconst.dki as dki
import dipy.reconst.dti as dti
from dipy.core.gradients import gradient_table
import nibabel as nib
import glob


def get_files_from_folder(folder_path, binaray = True):
    if binaray:
        txt_files = glob.glob(os.path.join(folder_path, f"*.bfloat"))
    else:
        txt_files = glob.glob(os.path.join(folder_path, f"*.txt"))
    txt_files.sort()

    return txt_files

def read_and_extract_bvals(file_path):
    column_4 = []  # Initialize an empty list to store the values from the 4th column
    try:
            with open(file_path, 'r') as file:
                for line in file:
                    # Split each line into space-separated values and extract the 4th column
                    columns = line.strip().split()
                    if len(columns) >= 5:
                        G = float(columns[3])*1e-3
                        giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
                        delta = float(columns[5])
                        Delta = float(columns[4])
                        b = pow(G * giro * delta, 2) * (Delta - delta/ 3) / 1000
                        column_4.append(b)  # Assuming columns are 0-based

    except FileNotFoundError:
        print(f"File not found: {file_path}")
        
    return np.array(column_4)


def read_and_extract_dwi(file_path, binary = True):
    column = []  # Initialize an empty list to store the values 
    if binary:
        lines = np.fromfile(file_path, dtype="float32")
        for line in lines:
  
            column.append(float(line))  # Assuming columns are 0-based
    else:
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    # Split each line into space-separated values and extract 
                    columns = line.strip().split()
                    column.append(float(columns[0]))  # Assuming columns are 0-based

        except FileNotFoundError:
            print(f"File not found: {file_path}")
    
    return np.array(column)

def read_and_extract_bvecs(file_path):
    columns_1_to_3 = []  # Initialize an empty list to store the values from the first three columns

    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Split each line into space-separated values and extract the first three columns
                columns = line.strip().split()[:3]  # Assuming columns are 0-based
                if len(columns) == 3:
                    b_value = float(line.strip().split()[3])
                    if b_value != 0:
                        columns = [float(val) for val in columns]  # Convert to float if needed
                        columns_1_to_3.append(columns)
                    else :

                        columns_1_to_3.append([0,0,0])

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    
    return np.array(columns_1_to_3)

def txt_to_nifti(path_to_DWI):
    dwi_array = read_and_extract_dwi(path_to_DWI)

    # Create an empty 4x4 affine matrix with ones on the diagonal
    affine = np.eye(4)

    img = nib.Nifti1Image(dwi_array, affine)
    return img

def array_to_nifti(dwi_array):

    # Create an empty 4x4 affine matrix with ones on the diagonal
    affine = np.eye(4)

    img = nib.Nifti1Image(dwi_array, affine)
    return img


def calculate_DKI(path_scheme, dwi):

    bvals = read_and_extract_bvals(path_scheme)

    bvecs = read_and_extract_bvecs(path_scheme)

    gtab = gradient_table(bvals, bvecs)

    #dwi = txt_to_nifti(path_to_DWI)

    # build model
    #dkimodel = dki.DiffusionKurtosisModel(gtab)
    dtimodel = dti.TensorModel(gtab)
    tenfit = dtimodel.fit(dwi.get_fdata())
    #dkifit = dkimodel.fit(dwi.get_fdata())
    # save maps

    #FA = dkifit.fa
    #MD = dkifit.md
    #AD = dkifit.ad
    #RD = dkifit.rd
    #MK = dkifit.mk(0, 10)
    #AK = dkifit.ak(0, 10)
    #RK = dkifit.rk(0, 10)

    FA = tenfit.fa
    MD = tenfit.md
    AD = tenfit.ad
    RD = tenfit.rd

    eigenvectors = tenfit.evecs  # Eigenvectors of the diffusion tensor

    # Each column of eigenvectors corresponds to an eigenvector
    # The orientation of the axial diffusion is given by the eigenvector corresponding to the largest eigenvalue
    axial_diffusion_orientation = eigenvectors[..., 0]  # Assuming the axial diffusion direction corresponds to the eigenvector with the largest eigenvalue

    RK= 0
    AK = 0
    MK = 0


    return FA, MD, AD, RD, MK, AK, RK, axial_diffusion_orientation

def add_to_data(dwi,  path_scheme, factor, repetition):


    img  = array_to_nifti(dwi)

    FA, MD, AD, RD, MK, AK, RK = calculate_DKI(path_scheme, img)

    if factor is not None:
        col_name = f"R/{factor}"
    else :
        col_name = f"cylinder"
    df = pd.DataFrame()
    df["overlapping_distance"] = [col_name]
    df["repetition"] = [repetition]
    df["MD"] = [MD]
    df["FA"] = [FA]
    df["AD"] = [AD]
    df["RD"] = [RD]
    df["MK"] = [float(MK)]
    df["AK"] = [float(AK)]
    df["RK"] = [float(RK)]


    return df


def create_df_overlapping(path_scheme):

    folder = f"/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/instructions/axons_vs_cylinders/data2/"

    files = get_files_from_folder(folder)


    files = [file for file in files if "info" not in file]
    repetitions = ["rep_00", "rep_01", "rep_02", "rep_03"]
    datas = []
    factors = [2,4,8,16,32]

    dwi =None

    for factor in factors:
        print("factor : ", factor)

        for repetition in repetitions:

            for file in files :

                if "img" not in file and "info" not in file and f"factor_{factor}_" in file and str(repetition) in file:
                    print(file)

                    dwi = read_and_extract_dwi(file)
                    
                if dwi is not None:
                    df= add_to_data(dwi,  path_scheme, factor, repetition)
                    if "intra" in file:
                        df["Location"] = ["intra"]*len(df)
                    elif "extra" in file:
                        df["Location"] = ["extra"]*len(df)
                    datas.append(df)
                dwi = None

        for file in files:
            if "rep" not in file and "img" not in file and "info" not in file and f"factor_{factor}_" in file :
                print(file)
                dwi = read_and_extract_dwi(file)
                    
            if dwi is not None:
                df= add_to_data(dwi,  path_scheme, factor, "rep_")
                if "intra" in file:
                        df["Location"] = ["intra"]*len(df)
                elif "extra" in file:
                    df["Location"] = ["extra"]*len(df)
                datas.append(df)
            dwi = None


    print("cylinder")
    for repetition in repetitions:
            for file in files :
                if "img" not in file and "info" not in file and "_factor_" not in file and str(repetition) in file:
                    print(file)
                    dwi = read_and_extract_dwi(file)
                    
                if dwi is not None:
                    df= add_to_data(dwi,  path_scheme, None, repetition)
                    if "intra" in file:
                        df["Location"] = ["intra"]*len(df)
                    elif "extra" in file:
                        df["Location"] = ["extra"]*len(df)
                    datas.append(df)
                dwi = None


    for file in files:
        if "rep" not in file and "img" not in file and "info" not in file and "factor" not in file :
            print(file)

            dwi = read_and_extract_dwi(file)

        if dwi is not None:
            df= add_to_data(dwi,  path_scheme, None, "rep_")
            if "intra" in file:
                df["Location"] = ["intra"]*len(df)
            elif "extra" in file:
                df["Location"] = ["extra"]*len(df)
            datas.append(df)
        dwi = None



    data = pd.concat(datas)
    for col in ["MD","FA","AD","RD"]:
        data = normalise_column(col, data)

    return data


def normalise_column(col, data):
    # Calculate the average "adc [um²/ms]" for type = "cylinders" for each combination
    data_cylinder = data[data['overlapping_distance'] == "cylinder"][[col, "Location"]]
    mean_intra = data_cylinder.loc[data_cylinder["Location"] == "intra"][col].mean()
    mean_extra = data_cylinder.loc[data_cylinder["Location"] == "extra"][col].mean()

    list(map(lambda x,location : x/mean_intra if location == "intra" else x/mean_extra, list(data[col]), list(data["Location"])))
    # Divide "adc [um²/ms]" by the average value
    title = 'normalized_'+col
    data[title] = list(map(lambda x,location : x/mean_intra if location == "intra" else x/mean_extra, list(data[col]), list(data["Location"])))
    return data



def plot():

    path ="/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/docs/scheme_files/PGSE_sample_scheme_21_dir.scheme"
    df = create_df_overlapping(path)

    # List of numerical column names to plot
    numerical_columns = [
        "normalized_RD",
        "normalized_AD",
        "normalized_MD",
        "normalized_FA"
        #"normalized_MK",
        #"normalized_AK",
        #"normalized_RK"
    ]



    # Create subplots for each numerical column
    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(16, 8))
    #fig.suptitle("Percentage of Volume Increase (Swelling) vs. Diffusion and Kurtosis Parameters", fontsize=16)

    # Flatten the axes for easy iteration
    # axes = axes.flatten()

    # Set a palette for the plots
    palette0 = sns.color_palette("Set2", len(numerical_columns))

    palette1 = sns.color_palette("crest", len(numerical_columns))

    palette2 = sns.color_palette("magma", len(numerical_columns))

    # Loop through numerical columns and create scatter plots with regression lines
    for i, column in enumerate(numerical_columns):
        location = "intra"

        if (i < 3):
            palette = palette1
        elif (i == 3):
            palette = palette0
        else:
            palette = palette2
        ax = axes[0,i]
        
        # Create a scatter plot
        sns.boxplot(data= df.loc[df["Location"]==location], x="overlapping_distance", y=column, ax=ax, color=palette[i])
        
        # Add a regression line
        #sns.regplot(data=df, x="overlapping_distance", y=column, ax=ax, color=palette[i], scatter=False)
        title_name = f"{location} :{column}"
        ax.set_title(title_name)
        ax.set_xlabel("overlapping_distance")
        ax.set_ylabel(column)

    for i, column in enumerate(numerical_columns):
        location = "extra"
        if (i < 3):
            palette = palette1
        elif (i == 3):
            palette = palette0
        else:
            palette = palette2
        ax = axes[1,i]
        
        # Create a scatter plot
        sns.boxplot(data= df.loc[df["Location"]=="extra"], x="overlapping_distance", y=column, ax=ax, color=palette[i])
        
        # Add a regression line
        #sns.regplot(data=df, x="overlapping_distance", y=column, ax=ax, color=palette[i], scatter=False)
        
        title_name = f"{location} :{column}"
        ax.set_title(title_name)
        ax.set_xlabel("overlapping_distance")
        ax.set_ylabel(column)

    # Remove any remaining empty subplots
    for j in range(len(numerical_columns), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()



