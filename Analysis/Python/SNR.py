# To compare ADC between multiple DWI files with varying N values

import numpy as np
import json
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")
sys.path.insert(1, '/home/localadmin/Documents/analytical_formula/')
from my_murdaycotts import my_murdaycotts
import math
import statannot
import array
import binascii
import nibabel as nib
from dipy.core.gradients import gradient_table
import dipy.reconst.dki as dki
import dipy.reconst.dti as dti
import math
import scipy

cur_path = os.getcwd()
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
scheme_file = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/docs/scheme_files/PGSE_21_dir_12_b.scheme"
icvf = 0.38

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

def calculate_DKI(path_scheme, dwi):

    bvals = read_and_extract_bvals(path_scheme)
    idx   = bvals <= 1
    bvals = bvals[idx]    
    bvecs = read_and_extract_bvecs(path_scheme)
    bvecs = bvecs[idx]    

    gtab = gradient_table(bvals, bvecs)

    #dwi = txt_to_nifti(path_to_DWI)

    # build model
    dkimodel = dki.DiffusionKurtosisModel(gtab)
    #dtimodel = dti.TensorModel(gtab)
    #tenfit = dtimodel.fit(dwi.get_fdata())
    dkifit = dkimodel.fit(dwi.get_fdata())
    # save maps

    FA = dkifit.fa
    MD = dkifit.md
    AD = dkifit.ad
    RD = dkifit.rd
    MK = dkifit.mk(0, 10)
    AK = dkifit.ak(0, 10)
    RK = dkifit.rk(0, 10)

    #FA = tenfit.fa
    #MD = tenfit.md
    #AD = tenfit.ad
    #RD = tenfit.rd

    #RK= 0
    #AK = 0
    #MK = 0


    return FA, MD, AD, RD, MK, AK, RK

def get_dwi_array(dwi_path):
    # create array with dwi values
    return np.fromfile(dwi_path, dtype="float32")

def get_psge_data():
    data_dwi = pd.DataFrame(columns = ["x", "y","z","G","Delta","delta","TE"])
    x, y, z, G, Delta, delta, TE = [], [], [], [], [], [], []
    with open(scheme_file) as f:
        for line in f.readlines():
            if len(line.split(' ')) > 3:
                for e, element in enumerate(line.split(' ')):
                    if e == 0:
                        x.append(float(element))
                    elif e == 1:
                        y.append(float(element))
                    elif e == 2:
                        z.append(float(element))
                    elif e == 3:
                        G.append(float(element)*1e-3)
                    elif e == 4:
                        Delta.append(float(element))
                    elif e == 5:
                        delta.append(float(element))
                    elif e == 6:
                        TE.append(float(element[:-1]))
    data_dwi["x"] = x
    data_dwi["y"] = y
    data_dwi["z"] = z
    data_dwi["G"] = G
    data_dwi["Delta"] = Delta
    data_dwi["delta"] = delta
    data_dwi["TE"] = TE
    data_dwi["b [ms/um²]"] = pow(data_dwi["G"]*giro*data_dwi["delta"],2) * (data_dwi["Delta"]-data_dwi["delta"]/3)/1000

    return data_dwi

def create_data(dwi_path, name, extension):
    dwi_signal_re = get_dwi_array(dwi_path / f"{name}.{extension}")
    dwi_signal_im = get_dwi_array(dwi_path / f"{name}_img.{extension}")
    dwi_signal_no_noise = np.sqrt(dwi_signal_re**2 + dwi_signal_im**2)/dwi_signal_re[0]
    data_psge  = get_psge_data()

    for SNR in [20, 60, 100]:
        sigma = 1/SNR
        dwi_signal_noise = np.sqrt((dwi_signal_re/dwi_signal_re[0]+ np.random.randn(1, dwi_signal_re.shape[0])*sigma)**2 
                                + (dwi_signal_im/dwi_signal_re[0]+ np.random.randn(1, dwi_signal_re.shape[0])*sigma)**2)
        data_psge[f"Sb/So_{SNR}"] = list(dwi_signal_noise.tolist()[0])


    # Create an empty 4x4 affine matrix with ones on the diagonal
    affine = np.eye(4)

    bvals = read_and_extract_bvals(scheme_file)
    idx   = bvals <= 1
    img   = nib.Nifti1Image(dwi_signal_no_noise[idx], affine)
    FA, MD, AD, RD, MK, AK, RK = calculate_DKI(scheme_file, img)


    # data_psge["Sb/So"] = list(dwi_signal_noise.reshape((-1, 1)))
    data_psge["Sb/So"] = list(dwi_signal_no_noise)
    nb_G = len(data_psge["G"].unique())
    nb_dir = int(len(dwi_signal_no_noise.reshape((-1, 1))) / nb_G)
    data_dwi = pd.DataFrame()
    for i in range(nb_dir):
        data_dir = data_psge.iloc[i*nb_G:(i+1)*nb_G]
        b0 = list(data_dir["b [ms/um²]"])[0]
        # Sb0 = list(data_dir.loc[data_dir["b [ms/um²]"]== b0]["DWI"])[0]
        # signal = list(map(lambda Sb : Sb/Sb0, list(data_dir["DWI"])))
        signal_log = list(map(lambda Sb : np.log(Sb), list(data_dir["Sb/So"])))
        adc = list(map(lambda b,Sb : -np.log(Sb)/(b-b0) if b!= b0 else np.nan, list(data_dir["b [ms/um²]"]),list(data_dir["Sb/So"])))
        # data_dir["Sb/So"] = signal
        data_dir["log(Sb/So)"] = signal_log
        data_dir["adc [ms/um²]"] = adc
        data_dir["FA"] = [FA] * len(adc)
        data_dir["MD"] = [MD] * len(adc)
        data_dir["AD"] = [AD] * len(adc)
        data_dir["RD"] = [RD] * len(adc)
        data_dir["MK"] = [MK] * len(adc)
        data_dir["AK"] = [AK] * len(adc)
        data_dir["RK"] = [RK] * len(adc)
        data_dwi = pd.concat([data_dwi, data_dir])
    return data_dwi


def create_df(DWI_folder):

    df_dwi = pd.DataFrame()
    df_crossings = pd.DataFrame()
    for subcase in os.listdir(DWI_folder):
        if subcase == "soma_dendrites" or subcase == "soma_dendrites_ex":
            # Iterate through the files in the folder
            for neuron in os.listdir(DWI_folder / subcase):
                for filename in os.listdir(DWI_folder / subcase / neuron):
                
                    if "simu" in filename:
                        N = int(filename.split('_')[1])
                        T = int(filename.split('_')[3])
                        with open(DWI_folder / subcase / neuron / filename, 'r') as file:
                            # Read the file line by line
                            for line in file:
                                # Check if the line contains the relevant information
                                if 'Number of particles eliminated due crossings' in line:
                                    # Split the line to get the number of particles as the last element
                                    num_particles_crossings = int(line.split()[-1])
                                    # Break the loop, as we have found the information we need
                                    break
                            d = {'nb_crossings': [num_particles_crossings], 'N': [N], 'T': [T]}
                            df_avg_crossings = pd.DataFrame(d)
                            df_crossings = pd.concat([df_crossings, df_avg_crossings])
                    # Check if the filename contains "_rep_" and "DWI"
                    if "DWI_img" in filename:
                        name = ('_').join(filename.split('_')[:-1])
                        N = int(name.split('_')[1])
                        T = int(name.split('_')[3])
                        extension = filename.split('_')[-1].split('.')[-1]
                        dwi_intra = create_data(DWI_folder / subcase / neuron , name, extension)
                        nb_G   = len(dwi_intra["G"].unique())
                        nb_dir = int(len(dwi_intra["x"].values) / nb_G)
                        for i in range(nb_G):
                            sb_so     = []
                            sb_so_20  = []
                            sb_so_60  = []
                            sb_so_100 = []
                            for j in range(nb_dir):
                                sb_so.append(dwi_intra.iloc[nb_G*j + i, :]["Sb/So"])
                                sb_so_20.append(dwi_intra.iloc[nb_G*j + i, :]["Sb/So_20"])
                                sb_so_60.append(dwi_intra.iloc[nb_G*j + i, :]["Sb/So_60"])
                                sb_so_100.append(dwi_intra.iloc[nb_G*j + i, :]["Sb/So_100"])
                                b_lab = dwi_intra.iloc[nb_G*j + i, :]["b [ms/um²]"]
                            mean     = np.mean(sb_so)
                            mean_20  = np.mean(sb_so_20)
                            mean_60  = np.mean(sb_so_60)
                            mean_100 = np.mean(sb_so_100)
                            means    = np.array([mean, mean_20, mean_60, mean_100])
                            SNR      = np.array(['inf', '20', '60', '100'])
                            b_labels = dwi_intra["b [ms/um²]"].unique()
                            d = {'loc': "intra", 'N': [N]*SNR.shape[0], 'T': [T]*SNR.shape[0], 'Sb/So': means, 'SNR': SNR, 'b [ms/um²]': [b_lab]*SNR.shape[0], 'neuron': [neuron]*SNR.shape[0], 'case': [subcase]*SNR.shape[0]}
                            df_avg_dwi = pd.DataFrame(d)
                            df_dwi = pd.concat([df_dwi, df_avg_dwi])

    return df_dwi, df_crossings


branching = "branching"

plot = True
log  = False

if log:
    y_lim_min = -5
    y_lim_max = 0.1
else:
    y_lim_min = 0.
    y_lim_max = 1.1

if plot:
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 22

    plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

DWI_folder = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/exch")
df_dwi, df_crossings = create_df(DWI_folder)
T_labels = df_dwi['T'].unique()
N_labels = df_dwi['N'].unique()
b_labels = df_dwi["b [ms/um²]"].unique()

# Analytical solutions
G         = np.array([0.015, 0.034, 0.048, 0.059, 0.068, 0.076, 0.083, 0.090, 0.096, 0.102, 0.107]) # in T/m
Delta     = np.array([0.05] * G.size)  # in s
delta     = np.array([0.0165] * G.size)# in s
TE        = np.array([0.067] * G.size) # in s
D0        = 2.5e-9 #m²/s
gamma     = 2.6751525e8 #rad/(s*T)
# bb        = gamma**2 * G**2 * delta**2 * (Delta - delta/3) # rad² * s / m²
# print((D0/(gamma*G))**(1/3))
# print("b val ", bb)

nb_neurites     = 20
r_soma    = 10e-6 #m
r_neurite = 0.5e-6 #m
if branching == 'branching':
    nb_branching = 3
    l_neurite       = 80e-6 #m
    volume_neurites = nb_neurites * (2**nb_branching + 1) * np.pi*r_neurite**2*l_neurite # in m³
else:
    l_neurite       = 240e-6 # m
    volume_neurites = nb_neurites * np.pi*r_neurite**2*l_neurite # in m³
volume_neurites = 8784.68 #5.64898e-06 (2 branching)
volume_soma     = 4/3 * np.pi * r_soma**3 # in m³
volume_soma     = volume_soma * 1e18
volume_neuron   = volume_neurites + volume_soma
neurite_fraction= volume_neurites / volume_neuron
soma_fraction   = volume_soma / volume_neuron

both_signal = []
for b_ in df_dwi['b [ms/um²]'].values:
    b_ = b_*1e9
    mlnS, mlnSneuman, mlnSnp, bardelta, b = my_murdaycotts(np.array([0.05]), np.array([0.0165]), r_soma, D0, np.array([b_]))
    Ain = np.sqrt(np.pi/(4 * b_ * D0)) * math.erf(np.sqrt(b_ * D0))
    both_signal.append(neurite_fraction * Ain + soma_fraction * math.exp(-mlnS))

df_dwi['Sb/So_norm'] = df_dwi['Sb/So'].values -  both_signal

fig, ax = plt.subplots(1, 1, figsize=(10, 7))
sns.violinplot(data=df_dwi, x='SNR', y='Sb/So_norm', hue='case', ax=ax, order=['inf', '100', '60', '20'])
plt.ylabel("Sb/So - analytical signal (soma + dendrites)")
plt.axhline(y=0, linestyle="dashed", color='gray')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, ['Soma-dendrites (disconnected)', 'Soma-dendrites (connected)'], loc='upper right')

couples = []
couples_end = []
for b in df_dwi['SNR'].unique():
    for i, branch in enumerate(df_dwi['case'].unique()):
        couples.append((b, branch))    
        print(b, branch)
        # print(scipy.stats.wilcoxon(df_dwi[(df_dwi["SNR"] == b) & (df_dwi["case"] == branch)]['Sb/So_norm'].values, nan_policy='omit'))
        print(scipy.stats.ttest_1samp(df_dwi[(df_dwi["SNR"] == b) & (df_dwi["case"] == branch)]['Sb/So_norm'].values, popmean=0., nan_policy='omit'))
      
for i in range(1, len(couples) + 1):
    if i % 2 == 0:
        couples_end.append((couples[i-2], couples[i-1]))

statannot.add_stat_annotation(
    ax,
    data=df_dwi,
    y='Sb/So_norm', x='SNR',
    hue='case',
    hue_order=['soma_dendrites', 'soma_dendrites_ex'],
    box_pairs=couples_end,
    test="t-test_ind",
    text_format="star",
    loc="inside"
    )

plt.show()
