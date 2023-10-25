# To compare ADC between multiple DWI files with varying N values

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
sys.path.insert(1, '/home/localadmin/Documents/analytical_formula/')
from my_murdaycotts import my_murdaycotts
import math
import statannot

cur_path = os.getcwd()
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
scheme_file = cur_path + "/results/funnel/overlap_4/n1/PGSE_21_dir_12_b.scheme"
icvf = 0.38
def get_dwi_array(dwi_path):

    if ".bfloat" in str(dwi_path):
        # create array with dwi values
        return np.fromfile(dwi_path, dtype="float32")
    else:
        signal = []
        with open(dwi_path) as f:
            [signal.append(float(line)) for line in f.readlines()]
        return np.array(signal)
        

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
    print(dwi_path / f"{name}.{extension}")
    dwi_signal_re = get_dwi_array(dwi_path / f"{name}.{extension}")
    dwi_signal    = dwi_signal_re
    data_psge     = get_psge_data()
    data_psge["DWI"] = list(dwi_signal)
    nb_G = len(data_psge["G"].unique())
    nb_dir = int(len(dwi_signal) / nb_G)
    data_dwi = pd.DataFrame()
    for i in range(nb_dir):
        data_dir = data_psge.iloc[i*nb_G:(i+1)*nb_G]
        b0 = list(data_dir["b [ms/um²]"])[0]
        Sb0 = list(data_dir.loc[data_dir["b [ms/um²]"]== b0]["DWI"])[0]
        signal = list(map(lambda Sb : Sb/Sb0, list(data_dir["DWI"])))
        signal_log = list(map(lambda Sb : np.log(Sb/Sb0), list(data_dir["DWI"])))
        adc = list(map(lambda b,Sb : -np.log(Sb/Sb0)/(b-b0) if b!= b0 else np.nan, list(data_dir["b [ms/um²]"]),list(data_dir["DWI"])))
        data_dir["Sb/So"] = signal
        data_dir["log(Sb/So)"] = signal_log
        data_dir["adc [ms/um²]"] = adc
        data_dwi = pd.concat([data_dwi, data_dir])
    return data_dwi


def create_df(DWI_folder):
    print(DWI_folder)
    df_dwi = pd.DataFrame()
    df_crossings = pd.DataFrame()
    for neuron in os.listdir(DWI_folder):
        print(neuron)

        # Open the file in read mode
        
        # Iterate through the files in the folder
        for filename in os.listdir(DWI_folder / neuron):
            if "simu" in filename:
                T = int(filename.split("_")[3])
                N = int(filename.split("_")[1])
                if T == 5000:
                    with open(DWI_folder / neuron / filename, 'r') as file:
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
            if ("DWI.txt" in filename and "soma_dendrites_ex" in filename) or "DWI.bfloat" in filename:
                T = int(filename.split("_")[3])
                N = int(filename.split("_")[1])
                print(filename)
                if T == 5000:
                    if "soma_dendrites_ex" in filename:
                        name = filename.split('.')[0]
                        old  = True
                    else:
                        name = filename.split('.')[0]
                        old  = False
                    extension = filename.split('.')[-1]
                    dwi_intra = create_data(DWI_folder / neuron, name, extension)
                    nb_G   = len(dwi_intra["G"].unique())
                    nb_dir = int(len(dwi_intra["x"].values) / nb_G)
                    for i in range(nb_G):
                        sb_so = []
                        for j in range(nb_dir):
                            sb_so.append(dwi_intra.iloc[nb_G*j + i, :]["Sb/So"])
                            b_lab = dwi_intra.iloc[nb_G*j + i, :]["b [ms/um²]"]
                        mean = np.mean(sb_so)
                        b_labels = dwi_intra["b [ms/um²]"].unique()
                        d = {'loc': "intra", 'N': N, 'T': T, 'Sb/So': mean, 'b [ms/um²]': b_lab, 'neuron': neuron, 'old': old}
                        df_avg_dwi = pd.DataFrame(d, index=[i])
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
    y_lim_max = 1

if plot:
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 16

    plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

DWI_folder_old = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/exchange/")
DWI_folder_new = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/no_funnel/overlap_4")

df_dwi_old, _ = create_df(DWI_folder_old)
df_dwi_new, _ = create_df(DWI_folder_new)
df_dwi_old["old"] = True
df_dwi_new["old"] = False

df_all = df_dwi_old.append(df_dwi_new)

T_labels = df_all['T'].unique()
N_labels = df_all['N'].unique()
b_labels = df_all["b [ms/um²]"].unique()

fig, ax = plt.subplots(1, 1)
sns.violinplot(data=df_all[df_all['b [ms/um²]'] > 0], x='b [ms/um²]', y='Sb/So', hue='old', hue_order=[True, False], ax=ax)
ax.set_xticklabels([f'{float(blab):.3f}' for blab in b_labels[1:]])

couples = []
couples_end = []
for b in df_all[df_all['b [ms/um²]'] > 0]['b [ms/um²]'].unique():
    for i, branch in enumerate(df_all['old'].unique()):
        couples.append((b, branch))    

for i in range(1, len(couples) + 1):
    if i % 2 == 0:
        couples_end.append((couples[i-2], couples[i-1]))

statannot.add_stat_annotation(
    ax,
    data=df_all[df_all['b [ms/um²]'] > 0],
    y='Sb/So', x='b [ms/um²]',
    hue='old',
    hue_order=[True, False],
    box_pairs=couples_end,
    test="Mann-Whitney",
    text_format="star",
    loc="inside"
    )

# Analytical solutions & Mesh
if plot:
    # Analytical solutions
    G         = np.array([0.015, 0.034, 0.048, 0.059, 0.068, 0.076, 0.083, 0.090, 0.096, 0.102, 0.107]) # in T/m
    Delta     = np.array([0.05] * G.size)  # in s
    delta     = np.array([0.0165] * G.size)# in s
    TE        = np.array([0.067] * G.size) # in s
    D0        = 2.5e-9 #m²/s
    gamma     = 2.6751525e8 #rad/(s*T)
    bb        = gamma**2 * G**2 * delta**2 * (Delta - delta/3) # rad² * s / m²
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
    volume_neurites = 3767.95 #8782.71 #5.64898e-06 (2 branching)
    volume_soma     = 4/3 * np.pi * r_soma**3 # in m³
    volume_soma     = volume_soma * 1e18
    volume_neuron   = volume_neurites + volume_soma
    neurite_fraction= volume_neurites / volume_neuron
    soma_fraction   = volume_soma / volume_neuron
    print("soma volume {:e}".format((volume_soma*1e18)))
    print("neurites volume {:e}".format((volume_neurites*1e18)))
    print("neuron {:e}".format((volume_neuron*1e18)))
    print("soma fraction {:e}".format(soma_fraction))

    soma_signal   = []
    soma_signal_neuman   = []
    neurites_signal = []
    both_signal     = []
    for i in range(bb.size):
        mlnS, mlnSneuman, mlnSnp, bardelta, b = my_murdaycotts(Delta[i], delta[i], r_soma, D0, bb[i])
        if log:
            soma_signal.append(-mlnS)
            soma_signal_neuman.append(-mlnSneuman)
        else:
            soma_signal.append(math.exp(-mlnS))
            soma_signal_neuman.append(math.exp(-mlnSneuman))

  
        # Signal intra sticks
        if log:
            Ain = np.log(np.sqrt(np.pi/(4 * bb[i] * D0)) * math.erf(np.sqrt(bb[i] * D0)))
        else:
            Ain = np.sqrt(np.pi/(4 * bb[i] * D0)) * math.erf(np.sqrt(bb[i] * D0))
        neurites_signal.append(Ain)
        if log:
            both_signal.append(neurite_fraction * Ain + soma_fraction * -mlnS)
        else:
            both_signal.append(neurite_fraction * Ain + soma_fraction * math.exp(-mlnS))


    i = 0
    for t_i, t in enumerate(T_labels):

        if plot:
            ax2 = ax.twinx()
            # Replace the NaN corresponding to b=0 to 1
            ax2.plot(b_labels[1:], soma_signal, label=f"Soma (analytic)", color='b')
            ax2.plot(b_labels[1:], neurites_signal, label=f"Neurites (analytic)", color='orange')
            ax2.plot(b_labels[1:], both_signal, label=f"Neurites & soma (analytic)", color='g')
            if log:
                sns.lineplot(b_labels, [-b_lab*D0*1e9 for b_lab in b_labels], label="D = 2.5 [ms/um²]")
            ax2.legend(loc=3)
            ax2.set_yticklabels([])
            ax2.set_ylim([y_lim_min, y_lim_max])
            ax.set_ylim([y_lim_min, y_lim_max])
            step_length = np.sqrt(6 * D0 * TE[0] / int(t))
            ax2.set_title(f"T = {T_labels[i]}, step length = {step_length*1e6:.3f} um")
            i = i + 1
if log:
    fig.suptitle('ln(S/S0) average over 21 directions, ' + branching, y=0.95)
else:
    fig.suptitle('S/S0 average over 21 directions, ' + branching, y=0.95)



plt.show()
