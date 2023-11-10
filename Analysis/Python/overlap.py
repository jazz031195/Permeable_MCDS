# To compare ADC between multiple DWI files with varying N values

import numpy as np
import json
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
import statannot
import dipy.reconst.dki as dki
from dipy.core.gradients import gradient_table
import nibabel as nib

cur_path    = os.getcwd()
scheme_file = cur_path + "/results/funnel/overlap_4/n1/PGSE_21_dir_12_b.scheme"
giro        = 2.6751525e8 # Gyromagnetic radio [rad/(s*T)]

def get_bvals(scheme_file_path):
    """
    Function that gets b-values from scheme file

    Args:
        scheme_file_path (str) : path of the scheme file

    Returns:
        (np.ndarray) : b-values [ms/um²]
    """

    b_values = []  # Initialize an empty list to store the values from the 4th column

    try:
        with open(scheme_file_path, 'r') as file:
            for line in file:
                # Split each line into space-separated values and extract the 4th column
                columns = line.strip().split()
                if len(columns) >= 5:
                    G     = float(columns[3]) * 1e-6 # [T/um]
                    giro  = 2.6751525e8 * 1e-3 # Gyromagnetic radio [rad/(ms*T)]
                    delta = float(columns[5]) * 1e3 # [ms]
                    Delta = float(columns[4]) * 1e3 # [ms]
                    b     = pow(G * giro * delta, 2) * (Delta - delta/3) # [ms/um²]
                    b_values.append(b)  # Assuming columns are 0-based

    except FileNotFoundError:
        print(f"File not found: {scheme_file_path}")
    
    return np.array(b_values)

def get_bvectors(scheme_file_path):
    """
    Function that gets b-vectors from scheme file

    Args:
        scheme_file_path (str) : path of the scheme file

    Returns:
        (np.ndarray) : Unitary b_vectors, 3D
    """

    b_vectors = []  # Initialize an empty list to store the values from the first three columns

    try:
        with open(scheme_file_path, 'r') as file:
            for line in file:
                # Split each line into space-separated values and extract the first three columns
                b_vec = line.strip().split()[:3]  # Assuming columns are 0-based
                # Skip the header
                if len(b_vec) == 3:
                    b_value = float(line.strip().split()[3])
                    if b_value != 0:
                        b_vec = [float(val) for val in b_vec]  # Convert to float if needed
                        b_vectors.append(b_vec)
                    else :
                        b_vectors.append([0 ,0, 0])

    except FileNotFoundError:
        print(f"File not found: {scheme_file_path}")
    
    return np.array(b_vectors)

def calculate_DKI(scheme_file_path, dwi):
    """
    Function that calculates DTI / DKI metrics. 
    For DTI, 2 b-vals (<= 1) and 6 directions are needed. 
    For DKI, 3 b-vals (<= 2-3) and 21 directions are needed
    
    Args:
        scheme_file_path (str) : path of the scheme file
        dwi (np.ndarray)       : 4D DWI image

    Returns:
        FA (np.float64) : Fractional Anisotropy
        MD (np.float64) : Mean diffusivity
        AD (np.float64) : Axial diffusivity
        RD (np.float64) : Radial diffusivity
        MK (np.float64) : Mean Kurtosis
        AK (np.float64) : Axial Kurtosis
        RK (np.float64) : Radial Kurtosis

    """

    # Get b-values <= 1 [ms/um^2] (DTI has Gaussian assumption => small b needed)
    bvals = get_bvals(scheme_file_path)      
    bvecs = get_bvectors(scheme_file_path)
   
    
    # DTI fit
    idx       = bvals <= 1
    bvals_dti = bvals[idx]  
    bvecs_dti = bvecs[idx]    
    gtab      = gradient_table(bvals_dti, bvecs_dti)
    # build model
    dkimodel  = dki.DiffusionKurtosisModel(gtab)
    # Create an empty 4x4 affine matrix with ones on the diagonal
    affine = np.eye(4)
    dwi_nii   = nib.Nifti1Image(dwi[idx], affine)
    dkifit    = dkimodel.fit(dwi_nii.get_fdata())
    # save maps
    FA = dkifit.fa
    MD = dkifit.md
    AD = dkifit.ad
    RD = dkifit.rd

    # DKI fit
    idx       = bvals <= 3
    bvals_dki = bvals[idx]  
    bvecs_dki = bvecs[idx]    
    gtab      = gradient_table(bvals_dki, bvecs_dki)
    # build model
    dkimodel  = dki.DiffusionKurtosisModel(gtab)
    dki_nii   = nib.Nifti1Image(dwi[idx], affine)
    dkifit    = dkimodel.fit(dki_nii.get_fdata())
    MK = dkifit.mk(0, 10)
    AK = dkifit.ak(0, 10)
    RK = dkifit.rk(0, 10)

    return FA, MD, AD, RD, MK, AK, RK

def get_dwi(dwi_path):
    """
    Function that gets the DWI values from the simulation output file
    
    Args:
        dwi_path (pathlib.PoxisPath) : path of the output DWI file
        
    Returns:
        (np.ndarray) : DWI values
    """

    if ".bfloat" in str(dwi_path):
        return np.fromfile(dwi_path, dtype="float32")
    elif ".txt" in str(dwi_path):
        signal = []
        with open(dwi_path) as f:
            [signal.append(float(line)) for line in f.readlines()]
        return np.array(signal)

def get_psge(scheme_file_path):
    """
    Function that gets the Pulse-Gradient Spin-Echo (PGSE) parameters from the scheme file. 
    
    Args:
        scheme_file_path (str) : path of the scheme file

    Returns:
        PGSE_params (pd.DataFrame) : Dataframe with columns = "x", "y", "z", "G [T/um]", "Delta [ms]", "delta [ms]", "TE [ms]"

    """
     
    PGSE_params = pd.DataFrame(columns = ["x", "y", "z", "G [T/um]", "Delta [ms]", "delta [ms]", "TE [ms]"])
    x, y, z, G, Delta, delta, TE = [], [], [], [], [], [], []
    with open(scheme_file_path) as f:
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
                        G.append(float(element) * 1e-6)      # [T/um]
                    elif e == 4:
                        Delta.append(float(element) * 1e3)   # [ms]
                    elif e == 5:
                        delta.append(float(element) * 1e3)   # [ms]
                    elif e == 6:
                        TE.append(float(element[:-1]) * 1e3) # [ms]
    PGSE_params["x"]          = x
    PGSE_params["y"]          = y
    PGSE_params["z"]          = z
    PGSE_params["G [T/um]"]   = G
    PGSE_params["Delta [ms]"] = Delta
    PGSE_params["delta [ms]"] = delta
    PGSE_params["TE [ms]"]    = TE
    PGSE_params["b [ms/um²]"] = pow(PGSE_params["G [T/um]"] * giro*1e-3 * PGSE_params["delta [ms]"], 2) * (PGSE_params["Delta [ms]"] - PGSE_params["delta [ms]"]/3)

    return PGSE_params

def create_data(data_folder, SNR, name, extension):
    """
    Function that creates the dataframe of the data from one experience (e.g. neuron 1, 5 repetitions). 
    
    Args:
        data_folder (pathlib.PoxisPath) : path of the folder where the data are stored
        SNR                (np.float64) : Signal to noise ratio, to be added as Gaussian noise (sigma=1/SNR) on the real & imaginary part of the signal
        name                      (str) : Name of the experiment
        extension                 (str) : extension of the file name

    Returns:
        data_dwi (pd.DataFrame) : Dataframe with columns = "x", "y", "z", "G [T/um]", "Delta [ms]", "delta [ms]", "TE [ms], "Sb/So","log(Sb/So)", 
                                                           "adc [ms/um²]", "FA", "MD", "AD", "RD", "MK", "AK", "RK"

    """

    dwi_real      = get_dwi(data_folder / f"{name}.{extension}")
    dwi_imaginary = get_dwi(data_folder / f"{name}_img.{extension}")
    dwi_no_noise  = np.sqrt(dwi_real**2 + dwi_imaginary**2)

    # Add Gaussian noise on the real and imaginary part => same as rician noise. Should converge to sqrt(pi/2)*sigma
    sigma = 1/SNR
    dwi_noise = np.sqrt((dwi_real/dwi_real[0]+ np.random.randn(1, dwi_real.shape[0])*sigma)**2 
                             + (dwi_imaginary/dwi_real[0]+ np.random.randn(1, dwi_real.shape[0])*sigma)**2)

    FA, MD, AD, RD, MK, AK, RK = calculate_DKI(scheme_file, dwi_no_noise)


    data_psge          = get_psge(scheme_file)
    Sb_So              = list(np.squeeze(dwi_noise.reshape((-1, 1))))

    data_psge["Sb/So"] = Sb_So
    # Number of b-values
    nb_b     = len(data_psge["b [ms/um²]"].unique())
    # Number of directions
    nb_dir   = int(len(dwi_noise.reshape((-1, 1))) / nb_b)
    # Data with all the directions
    data_dwi = pd.DataFrame()
    for i in range(nb_dir):
        # Data for one direction
        data_dir   = data_psge.iloc[i * nb_b : (i + 1) * nb_b]
        b0         = list(data_dir["b [ms/um²]"])[0]
        data_dir["log(Sb/So)"]   = list(map(lambda Sb : np.log(Sb), list(data_dir["Sb/So"])))
        adc                      = list(map(lambda b,Sb : -np.log(Sb)/(b-b0) if b != b0 else np.nan, list(data_dir["b [ms/um²]"]), list(data_dir["Sb/So"])))
        data_dir["adc [ms/um²]"] = adc
        data_dir["FA"]           = [FA] * nb_b
        data_dir["MD"]           = [MD] * nb_b
        data_dir["AD"]           = [AD] * nb_b
        data_dir["RD"]           = [RD] * nb_b
        data_dir["MK"]           = [MK] * nb_b
        data_dir["AK"]           = [AK] * nb_b
        data_dir["RK"]           = [RK] * nb_b

        data_dwi = pd.concat([data_dwi, data_dir])

    return data_dwi


def create_df_all(experience_folder):
    """
    Creates a dataframe with all the simulations together (e.g. neuron 1 + neuron 2 + ...)

    Args:
        experience_folder (pathlib.PosixPath) : folder where all the experiences (all substrate, all repetitions) are stored

    Returns:
        df_all_data  (pd.DataFrame) : Dataframe containing all the data needed for the plots, statistics, etc
        df_crossings (pd.DataFrame) : Dataframe containing the crossings information
        
    """

    df_all_data  = pd.DataFrame()
    df_crossings = pd.DataFrame()
    for overlap in os.listdir(experience_folder):
        print(overlap)
        for neuron in os.listdir(experience_folder / overlap):
            if neuron != "n0" and neuron != "n6":
                print(neuron)

                # Read useful parameters from file
                f    = open(experience_folder / overlap / neuron / "params.json")
                data = json.load(f)
                # N                 = data.get("N")                 # Number of Walkers / water particles
                # T                 = data.get("T")                 # Number of timesteps
                duration          = data.get("duration")          # Simulation duration in s
                diff_intra        = data.get("diffusivity_intra") # Simulation diffusivity in m²/s
                diff_extra        = data.get("diffusivity_extra") # Simulation diffusivity in m²/s
                sphere_overlap    = data.get("sphere_overlap")    # the spheres are radius/sphere_overlap appart
                funnel            = data.get("funnel")            # boolean, do a funnel between soma & dendrites
                
                # Iterate through the files in the folder
                for filename in os.listdir(experience_folder / overlap / neuron):
                    # Read the simulation_info.txt to have crossings information
                    if "simu" in filename:
                        with open(experience_folder / overlap / neuron / filename, 'r') as file:
                            N = int(filename.split('_')[1])
                            T = int(filename.split('_')[3])
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
                            df_crossings     = pd.concat([df_crossings, df_avg_crossings])
                    
                    # Check if the filename contains "_rep_" and "DWI"
                    if "DWI_img" in filename:
                        # Name of the experience
                        name         = ('_').join(filename.split('_')[:-1])
                        # Number of walkers
                        N            = int(name.split('_')[1])
                        # Number of timesteps
                        T            = int(name.split('_')[3])
                        extension    = filename.split('_')[-1].split('.')[-1]
                        SNR          = np.inf
                        data_one_exp = create_data(experience_folder / overlap / neuron, SNR, name, extension)
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
                            d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, "adc [ms/um²]": mean_adc, 
                                'b [ms/um²]': bval, 'neuron': neuron, 'overlap': int(overlap.split('_')[-1]),
                                'FA': FA, 'MD': MD, 'RD': RD, 'AD': AD, 'MK': MK, 'RK': RK, 'AK': AK}
                            df_avg_data = pd.DataFrame(d, index=[i])
                            df_all_data = pd.concat([df_all_data, df_avg_data])

    return df_all_data, df_crossings


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
    MEDIUM_SIZE = 19
    BIGGER_SIZE = 19

    plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)   # fontsize of the figure title


experience_folder = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/no_funnel/")
df_data_all, df_crossings = create_df_all(experience_folder)

T_labels = df_data_all['T'].unique()
N_labels = df_data_all['N'].unique()
b_labels = df_data_all["b [ms/um²]"].unique()

fig, ax = plt.subplots(1, 1)
sns.violinplot(data=df_data_all[(df_data_all['b [ms/um²]'] > 0) & (df_data_all["T"] == 15000) & (df_data_all["overlap"] != 1)], 
               x='b [ms/um²]', 
               y='Sb/So',
               hue='overlap', 
               ax=ax)
# Change b-values so that they have only 1 decimals
ax.set_xticklabels([f'{float(blab):.1f}' for blab in b_labels[1:]])
handles, labels = ax.get_legend_handles_labels()
labels          = ['R/' + item for item in labels]
ax.legend(handles, labels, loc='upper right', frameon=False, title="Overlap")

# couples = []
# couples_end = []
# for b in df_dwi[df_dwi['b [ms/um²]'] > 0]['b [ms/um²]'].unique():
#     for i, branch in enumerate(df_dwi['overlap'].unique()):
#         couples.append((b, branch))    

# for i in range(1, len(couples) + 1):
#     if i % 6 == 0:
#         couples_end.append((couples[i-6], couples[i-5]))
#         couples_end.append((couples[i-6], couples[i-4]))
#         couples_end.append((couples[i-6], couples[i-3]))
#         couples_end.append((couples[i-6], couples[i-2]))
#         couples_end.append((couples[i-6], couples[i-1]))
#         couples_end.append((couples[i-5], couples[i-4]))
#         couples_end.append((couples[i-5], couples[i-3]))
#         couples_end.append((couples[i-5], couples[i-2]))
#         couples_end.append((couples[i-5], couples[i-1]))
#         couples_end.append((couples[i-4], couples[i-3]))
#         couples_end.append((couples[i-4], couples[i-2]))
#         couples_end.append((couples[i-4], couples[i-1]))
#         couples_end.append((couples[i-3], couples[i-2]))
#         couples_end.append((couples[i-3], couples[i-1]))
#         couples_end.append((couples[i-2], couples[i-1]))

# statannot.add_stat_annotation(
#     ax,
#     data=df_dwi[(df_dwi['b [ms/um²]'] > 0) & (df_dwi["T"] == 15000)],
#     y='Sb/So', x='b [ms/um²]',
#     hue='overlap',
#     box_pairs=couples_end,
#     test="Mann-Whitney",
#     text_format="star",
#     loc="inside"
#     )
# ax.set_title(f'N = {N_labels[0]}, T = {T_labels[0]}')





# # Analytical solutions & Mesh
# if plot:
#     # Analytical solutions
#     G         = np.array([0.015, 0.034, 0.048, 0.059, 0.068, 0.076, 0.083, 0.090, 0.096, 0.102, 0.107]) # in T/m
#     Delta     = np.array([0.05] * G.size)  # in s
#     delta     = np.array([0.0165] * G.size)# in s
#     TE        = np.array([0.067] * G.size) # in s
#     D0        = 2.5e-9 #m²/s
#     gamma     = 2.6751525e8 #rad/(s*T)
#     bb        = gamma**2 * G**2 * delta**2 * (Delta - delta/3) # rad² * s / m²
#     # print((D0/(gamma*G))**(1/3))
#     # print("b val ", bb)

#     nb_neurites     = 20
#     r_soma    = 10e-6 #m
#     r_neurite = 0.5e-6 #m
#     if branching == 'branching':
#         nb_branching = 3
#         l_neurite       = 80e-6 #m
#         volume_neurites = nb_neurites * (2**nb_branching + 1) * np.pi*r_neurite**2*l_neurite # in m³
#     else:
#         l_neurite       = 240e-6 # m
#         volume_neurites = nb_neurites * np.pi*r_neurite**2*l_neurite # in m³
#     volume_neurites = 8784.68 #5.64898e-06 (2 branching)
#     volume_soma     = 4/3 * np.pi * r_soma**3 # in m³
#     volume_soma     = volume_soma * 1e18
#     volume_neuron   = volume_neurites + volume_soma
#     neurite_fraction= volume_neurites / volume_neuron
#     soma_fraction   = volume_soma / volume_neuron
#     print("soma volume {:e}".format((volume_soma*1e18)))
#     print("neurites volume {:e}".format((volume_neurites*1e18)))
#     print("neuron {:e}".format((volume_neuron*1e18)))
#     print("soma fraction {:e}".format(soma_fraction))

#     soma_signal   = []
#     soma_signal_neuman   = []
#     neurites_signal = []
#     both_signal     = []
#     for i in range(bb.size):
#         mlnS, mlnSneuman, mlnSnp, bardelta, b = my_murdaycotts(Delta[i], delta[i], r_soma, D0, bb[i])
#         if log:
#             soma_signal.append(-mlnS)
#             soma_signal_neuman.append(-mlnSneuman)
#         else:
#             soma_signal.append(math.exp(-mlnS))
#             soma_signal_neuman.append(math.exp(-mlnSneuman))

  
#         # Signal intra sticks
#         if log:
#             Ain = np.log(np.sqrt(np.pi/(4 * bb[i] * D0)) * math.erf(np.sqrt(bb[i] * D0)))
#         else:
#             Ain = np.sqrt(np.pi/(4 * bb[i] * D0)) * math.erf(np.sqrt(bb[i] * D0))
#         neurites_signal.append(Ain)
#         if log:
#             both_signal.append(neurite_fraction * Ain + soma_fraction * -mlnS)
#         else:
#             both_signal.append(neurite_fraction * Ain + soma_fraction * math.exp(-mlnS))


#     i = 0
#     for t_i, t in enumerate(T_labels):

#         if plot:
#             ax2 = ax.twinx()
#             # Replace the NaN corresponding to b=0 to 1
#             ax2.plot(b_labels[1:], soma_signal, label=f"Soma (analytic)", color='b', linestyle="dotted")
#             ax2.plot(b_labels[1:], neurites_signal, label=f"Neurites (analytic)", color='orange', linestyle="dotted")
#             ax2.plot(b_labels[1:], both_signal, label=f"Neurites & soma (analytic)", color='g', linestyle="dotted")
#             # sns.lineplot(b_labels, soma_signal, label=f"Soma (analytic)", color='b', ax=ax)
#             # # ax2.errorbar([b_lab + 0.05 for b_lab in b_labels], soma_signal_neuman, 
#             # #                 yerr=[0], label=f"Soma (analytic, Neuman)", fmt='o', color='blue')
#             # sns.lineplot(b_labels, neurites_signal, label=f"Neurites (analytic)", color='orange', ax=ax)
#             # sns.lineplot(b_labels, both_signal, label=f"Neurites & soma (analytic)", color='g', ax=ax)
#             # if log:
#             #     sns.lineplot(b_labels, [-b_lab*D0*1e9 for b_lab in b_labels], label="D = 2.5 [ms/um²]")
#             ax2.legend(loc=3)
#             ax2.set_yticklabels([])
#             ax2.set_ylim([y_lim_min, y_lim_max])
#             ax.set_ylim([y_lim_min, y_lim_max])
#             # step_length = np.sqrt(6 * D0 * TE[0] / int(t))
#             # # ax2.set_title(f"T = {T_labels[i]}, step length = {step_length*1e6:.3f} um")
#             # i = i + 1
# if log:
#     fig.suptitle('ln(S/S0) average over 21 directions, ' + branching, y=0.95)
# else:
#     fig.suptitle('S/S0 average over 21 directions, ' + branching, y=0.95)

# plt.show()

# fig, axes = plt.subplots(2, 2)
# axes = axes.ravel()
# sns.violinplot(data=df_dwi[df_dwi['b [ms/um²]'] > 0], x='overlap', y='FA', hue='neuron', hue_order=['n1', 'n2', 'n3', 'n4', 'n5'], inner="points", ax=axes[0])
# sns.violinplot(data=df_dwi[df_dwi['b [ms/um²]'] > 0], x='overlap', y='MD', hue='neuron', hue_order=['n1', 'n2', 'n3', 'n4', 'n5'], inner="points", ax=axes[1])
# sns.violinplot(data=df_dwi[df_dwi['b [ms/um²]'] > 0], x='overlap', y='AD', hue='neuron', hue_order=['n1', 'n2', 'n3', 'n4', 'n5'], inner="points", ax=axes[2])
# sns.violinplot(data=df_dwi[df_dwi['b [ms/um²]'] > 0], x='overlap', y='RD', hue='neuron', hue_order=['n1', 'n2', 'n3', 'n4', 'n5'], inner="points", ax=axes[3])
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, labels, loc='center right', title='neuron')
# axes[0].legend([],[], frameon=False)
# axes[1].legend([],[], frameon=False)
# axes[2].legend([],[], frameon=False)
# axes[3].legend([],[], frameon=False)
# ax.set_xticklabels([f'{float(blab):.3f}' for blab in b_labels[1:]])

# couples = []
# couples_end = []
# for b in df_dwi[df_dwi['b [ms/um²]'] > 0]['b [ms/um²]'].unique():
#     for i, branch in enumerate(df_dwi['overlap'].unique()):
#         couples.append((b, branch))    

# for i in range(1, len(couples) + 1):
#     if i % 2 == 0:
#         couples_end.append((couples[i-2], couples[i-1]))

# statannot.add_stat_annotation(
#     ax,
#     data=df_dwi[df_dwi['b [ms/um²]'] > 0],
#     y='adc [ms/um²]', x='b [ms/um²]',
#     hue='overlap',
#     box_pairs=couples_end,
#     test="Mann-Whitney",
#     text_format="star",
#     loc="inside"
#     )
# ax.set_title(f'N = {N_labels[0]}, T = {T_labels[0]}')
# plt.show()


# fig, axes = plt.subplots(2, 2)
# axes = axes.ravel()
# sns.boxplot(data=df_dwi[df_dwi['b [ms/um²]'] > 0], x='overlap', y='FA', ax=axes[0])
# sns.boxplot(data=df_dwi[df_dwi['b [ms/um²]'] > 0], x='overlap', y='MD', ax=axes[1])
# sns.boxplot(data=df_dwi[df_dwi['b [ms/um²]'] > 0], x='overlap', y='AD', ax=axes[2])
# sns.boxplot(data=df_dwi[df_dwi['b [ms/um²]'] > 0], x='overlap', y='RD', ax=axes[3])
# plt.show()

fig, _ = plt.subplots(1,1)
ax = sns.boxplot(data=df_data_all[(df_data_all['b [ms/um²]'] > 0) & (df_data_all["T"] == 15000) & (df_data_all["overlap"] != 1)], 
                 x='overlap', 
                 y='MD', 
                 palette="Blues")
labels = ['R/' + item.get_text() for item in ax.get_xticklabels()]
ax.set_xticklabels(labels)
ax.set_xlabel('Distance between overlapping spheres')
ax.set_ylabel('Mean diffusivity')
plt.show()