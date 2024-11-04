import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import pandas as pd
import plotly.graph_objects as go

def get_traj(traj_path):
    """
    Function that gets the traj values from the simulation output file
    
    Args:
        traj_path (pathlib.PoxisPath) : path of the output DWI file
        
    Returns:
        (np.ndarray) : traj values
    """

    if ".txt" in str(traj_path):
        signal = []
        with open(traj_path) as f:
            values = f.readlines()
            for value in values:
                if value != '\n':
                    signal.append(float(value))
        return np.array(signal)

plot_traj  = True
neuron_file = '/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/T_ex/neurons_list_22_5_perc.txt'
traj_file   = ['/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/T_ex/_rep_00_16.traj.txt']



neuron_id = 0
df_all = pd.DataFrame(columns=["neuron_id", "x", "y", "z", "r"])
df_tmp2 = pd.DataFrame(columns=["neuron_id", "x", "y", "z", "r"])
if plot_traj:
    lines = get_traj(traj_file[0])
    xp = []
    yp = []
    zp = []
    for i in range(int(len(lines))):
        if i%3 == 0:
            xp.append(float(lines[i]))
        elif i%3 == 1:
            yp.append(float(lines[i]))
        elif i%3 == 2:
            zp.append(float(lines[i]))

    min_ = [min(xp[:1000]), min(yp[:1000]), min(zp[:1000])]
    max_ = [max(xp[:1000]), max(yp[:1000]), max(zp[:1000])]
    d = {'x':xp[:5000], 'y': yp[:5000], 'z': zp[:5000], 'r':0.1/1000, 'traj': 1}
    df_tmp2 = pd.concat([df_tmp2, pd.DataFrame(d)])

with open(neuron_file) as f:
    lines = f.readlines()
    lines = lines[2:]


    for i in range(len(lines)):
        coords = lines[i].split(' ')
        # If it is the soma, plot it in any case
        if "end" in coords[0]:
            neuron_id += 1
        elif len(coords) > 3:
            coords = [float(coord)/1000 for coord in coords]

            count = 0
            for k in range(3):
                if((coords[k+3] > min_[k] - 10/1000) & (coords[k+3] < max_[k] + 10/1000)):
                    count += 1
            count = 3
            if count == 3:
                print(coords[k+3], min_[k], max_[k])
                d = {"neuron_id": neuron_id, 
                        "x": coords[3], 
                        "y": coords[4], 
                        "z": coords[5], 
                        "r": coords[6],
                        "traj": 0}

                df_avg_data = pd.DataFrame(d, index=[neuron_id])
                df_all = pd.concat([df_all, df_avg_data])

print(df_all)
df_all = pd.concat([df_all, df_tmp2])
fig = go.Figure()    
fig.add_trace(go.Scatter3d(
                            x=df_all["x"],
                            y=df_all["y"],
                            z=df_all["z"],
                            type="scatter3d",
                            mode="markers",
                            marker=dict(
                                sizemode="diameter",
                                size=df_all["r"]*15000,
                                color=df_all["traj"],
                                line=dict(
                                    color="rgba(0, 0, 0, 0)",
                                    width=0
                                )
                            )
                        )
                )

fig.show()

