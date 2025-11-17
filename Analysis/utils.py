import numpy as np
import dipy.reconst.dti as dti
import dipy.reconst.dki as dki
from dipy.core.gradients import gradient_table
import nibabel as nib
import pandas as pd
from my_murdaycotts import my_murdaycotts
import math
import os
import warnings

giro = 2.6751525e8 # Gyromagnetic radio [rad/(s*T)]

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
    bvalues = get_bvals(scheme_file_path)      
    bvecs = get_bvectors(scheme_file_path)
   
    
    # DTI fit
    idx       = bvalues <= 1
    bvals_dti = bvalues[idx]  
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
    idx       = bvalues <= 3
    bvals_dki = bvalues[idx]  
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

def create_data(data_folder, SNR, name, extension, scheme_file_path):
    """
    Function that creates the dataframe of the data from one experience (e.g. neuron 1, 5 repetitions). 
    
    Args:
        data_folder (pathlib.PoxisPath) : path of the folder where the data are stored
        SNR                (np.float64) : Signal to noise ratio, to be added as Gaussian noise (sigma=1/SNR) on the real & imaginary part of the signal
        name                      (str) : Name of the experiment
        extension                 (str) : extension of the file name
        scheme_file_path          (str) : path of the scheme file

    Returns:
        data_dwi (pd.DataFrame) : Dataframe with columns = "x", "y", "z", "G [T/um]", "Delta [ms]", "delta [ms]", "TE [ms], "Sb/So","log(Sb/So)", 
                                                           "adc [ms/um²]", "FA", "MD", "AD", "RD", "MK", "AK", "RK"

    """

    dwi_real      = get_dwi(data_folder / f"{name}.{extension}")
    # There is an imaginary part to the signal
    if os.path.exists(data_folder / f"{name}_img.{extension}"):
        dwi_imaginary = get_dwi(data_folder / f"{name}_img.{extension}")
        dwi_no_noise  = np.sqrt(dwi_real**2 + dwi_imaginary**2)

        # Add Gaussian noise on the real and imaginary part => same as rician noise. Should converge to sqrt(pi/2)*sigma
        sigma = 1/SNR
        dwi_noise = np.sqrt((dwi_real/dwi_real[0] + np.random.randn(1, dwi_real.shape[0])*sigma)**2 
                            + (dwi_imaginary/dwi_real[0] + np.random.randn(1, dwi_real.shape[0])*sigma)**2)
    else:
        dwi_no_noise = dwi_real
        # Add Gaussian noise on the real and imaginary part => same as rician noise. Should converge to sqrt(pi/2)*sigma
        sigma = 1/SNR
        dwi_noise = (dwi_real/dwi_real[0] + np.random.randn(1, dwi_real.shape[0])*sigma)
        warnings.warn("Warning...........The signal is purely real")

    FA, MD, AD, RD, MK, AK, RK = calculate_DKI(scheme_file_path, dwi_no_noise)


    data_psge          = get_psge(scheme_file_path)
    Sb_So              = list(np.squeeze(dwi_noise.reshape((-1, 1))))

    data_psge["Sb/So"] = Sb_So
    data_psge["SNR"]   = [SNR] * len(Sb_So)
        
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


def analytical_solutions(bvalues, Delta, delta, sphere_radius, D0, log, sphere_fraction, stick_fraction):
    """
    Function that calculates the analytical solutions of diffusion in :
        - A sphere of radius sphere_radius (Signal_sphere)
        - Infinitely long randomly-oriented sticks (ROS) (Signal_ROS)
        - Both of the two signals : 
            sphere_fraction * Signal_sphere + stick_fraction * Signal_ROS (Signal_both)
    
    Args:
        bvalues    (np.ndarray) : b-values [s/m²]
        Delta      (np.ndarray) : Delta [s]
        delta      (np.ndarray) : delta [s]
        sphere_radius   (float) : radius of the sphere [m]
        D0              (float) : water diffusivity [m²/s]
        log              (bool) : True => log(Signal), False => Signal
        sphere_fraction (float) : fraction composed of spheres
        stick_fraction  (float) : fraction composed of sticks

    Returns:
        signal_sphere (list) : analytical signal from diffusion in a sphere
        signal_sticks (list) : analytical signal from diffusion in ROS
        signal_both   (list) : analytical signal from diffusion in both sphere + ROS
    """
    
    signal_sphere = []
    signal_sticks = []
    signal_both   = []
    for i in range(bvalues.shape[0]):
        # Signal intra sphere
        mlnS, _, _, _, _ = my_murdaycotts(Delta, delta, sphere_radius, D0, bvalues[i])
        if log:
            signal_sphere.append(-mlnS)
        else:
            signal_sphere.append(math.exp(-mlnS))

        # Signal intra sticks
        if bvalues[i] != 0:
            Ain = np.sqrt(np.pi/(4 * bvalues[i] * D0)) * math.erf(np.sqrt(bvalues[i] * D0))
        else:
            Ain = 1

        if log:
            Ain = np.log(Ain)
        
        signal_sticks.append(Ain)
        
        # Signal intra both
        if log:
            signal_both.append(stick_fraction * Ain + sphere_fraction * -mlnS)
        else:
            signal_both.append(stick_fraction * Ain + sphere_fraction * math.exp(-mlnS))

    return signal_sphere, signal_sticks, signal_both