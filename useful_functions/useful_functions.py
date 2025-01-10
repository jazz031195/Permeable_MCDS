import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")
import glob
import re
import nibabel as nib

def array_to_nifti(dwi_array):

    # Create an empty 4x4 affine matrix with ones on the diagonal
    affine = np.eye(4)

    img = nib.Nifti1Image(dwi_array, affine)
    return img

def extract_simulation_time(file_path):
    # Regular expression to find the duration line
    duration_line_pattern = re.compile(r"All (\d+) simulations ended after: (\d+) minutes and (\d+) seconds")

    try:
        # Read the content of the text file
        with open(file_path, 'r') as file:
            text_content = file.read()

        # Search for the pattern in the text content
        match = duration_line_pattern.search(text_content)

        if match:
            nbr_simulations = int(match.group(1))
            minutes = int(match.group(2))
            seconds = int(match.group(3))
            total_seconds = minutes * 60 + seconds
            return nbr_simulations, total_seconds
        else:
            return np.nan, np.nan
    except FileNotFoundError:
        return np.nan, np.nan
    
def extract_simulation_info(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Regular expressions to extract the number of particles and steps
    particles_pattern = r'Number of particles:\s*-+\s*(\d+)'
    steps_pattern = r'Number of steps:\s*-+\s*(\d+)'

    # Find the matches
    particles_match = re.search(particles_pattern, content)
    steps_match = re.search(steps_pattern, content)

    # Extract the values
    number_of_particles = int(particles_match.group(1)) if particles_match else None
    number_of_steps = int(steps_match.group(1)) if steps_match else None

    return number_of_particles, number_of_steps

def get_scheme_info_iso(path_to_info):
    with open(path_to_info, "r") as f:
        lines = f.readlines()
    vecs = []
    bs = []
    for i, line in enumerate(lines):
        vec = np.array([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
        b = line.split()[3]
        vecs.append(vec)
        bs.append(int(b)/1000)
    return np.array(bs), np.array(vecs)

def get_scheme_info(scheme_path):
    giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
    scheme = pd.read_csv(scheme_path, sep=" ", header=None, skiprows=1)
    scheme = scheme.dropna(axis=1)
    scheme = scheme.to_numpy()
    gradient_strength = scheme[:, 3]
    Delta, delta, TE = scheme[:, 4], scheme[:, 5], scheme[:, 6]
    b_values = ((giro*delta*gradient_strength)**2)*(Delta-delta/3)/1e9
    directions = scheme[:, :3]
    b_values = [round(float(i),2) for i in b_values]
    return b_values, directions

def read_and_extract_parameters(file_path):
    bs = []  # Initialize an empty list to store the values from the 4th column
    Gs = []
    deltas = []
    Deltas = []
    vectors = []
    TEs = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Split each line into space-separated values and extract the 4th column
                columns = line.strip().split()
                if len(columns) >= 5:
                    G = float(columns[3]) # Gradient strength given in T/m
                    giro = 2.6751525e5 #Gyromagnetic radio given in rad/(ms*T)
                    delta = float(columns[5])
                    Delta = float(columns[4])
                    TE = float(columns[-1])
                    b = int(pow(G * giro * delta, 2) * (Delta - delta/ 3))
                    #round b
                    b = round(b, 2)
                    bs.append(b)  # Assuming columns are 0-based
                    Gs.append(G)
                    deltas.append(delta)
                    Deltas.append(Delta)
                    Vector = [float(val) for val in columns[:3]]
                    vectors.append(Vector)
                    TEs.append(TE)

    except FileNotFoundError:
        print(f"File not found: {file_path}")
        
    return np.array(bs), np.array(Gs), np.array(deltas), np.array(Deltas), np.array(vectors), np.array(TEs)

def read_binary_file(file_name):
    """
    Reads a binary file and returns the data as a numpy array
    """
    with open(file_name, "rb") as f:
        data = np.fromfile(f, dtype=np.float32)
    return data

def get_files_from_folder(folder_path, binary=True):
    if binary:
        binary_files = glob.glob(os.path.join(folder_path, "*.bfloat"))
    else:
        binary_files = glob.glob(os.path.join(folder_path, "*.txt"))
    binary_files.sort()

    return binary_files

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

def get_csv_files_from_folder(folder_path):

    binary_files = glob.glob(os.path.join(folder_path, "*.csv"))
    binary_files.sort()

    return binary_files

def get_traj_files_from_folder(folder_path):

    files = glob.glob(os.path.join(folder_path, "*.traj"))
    files.sort()

    return files

def get_info_files_from_path(folder_path):

    files = glob.glob(os.path.join(folder_path, "*.txt"))
    return files

def plot_sphere_from_directions_in_scheme():
    """
    Plots a sphere with the directions from the scheme file.
    """
    scheme_file = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/docs/scheme_files/PGSE_sample_scheme_21_dir.scheme"
    
    scheme = pd.read_csv(scheme_file, sep=" ", header=None, skiprows=1)
    scheme = scheme.dropna(axis=1)
    scheme = scheme.to_numpy()
    gradient_strength = scheme[:, 3]
    scheme = scheme[:, :3]
    #normalise the scheme
    scheme = scheme / np.linalg.norm(scheme, axis=1)[:, None]

    scheme = scheme * gradient_strength[:, None]
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(scheme[:, 0], scheme[:, 1], scheme[:, 2])
    # label x, y, z axis
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

    scheme_file = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/docs/scheme_files/PGSE_sample_scheme_21_dir.scheme"


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

def get_total_icvf(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("Total icvf"):
                # Extract the value after the label
                total_icvf = float(line.split()[-1])
                return total_icvf
    return None  # Return None if "Total icvf" not found
