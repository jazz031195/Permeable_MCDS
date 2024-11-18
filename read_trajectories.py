import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import plotly.graph_objects as go
import plotly.colors as colors

def read_bin_file(file_path):
    traj_part = np.fromfile(file_path, dtype="float32")
    chunks = 4
    nbr_steps=64300
    valid_traj_part = traj_part[np.isfinite(traj_part)]  # Ensure all entries are valid floats
    
    scatters = []
    for i in range(chunks):
        start = (nbr_steps+1) * i
        end = nbr_steps * (i + 1)-1
        
        chunk = valid_traj_part[3*start:3*end]
        xs, ys, zs = chunk[0::3], chunk[1::3], chunk[2::3]

        print("done")
        # Example scatter plot
        colours = colors.qualitative.Plotly[:10]
        c = colours[0]  # Assuming e is defined somewhere in your code

        scatter = go.Scatter3d(
            x=[i for i in xs],
            y=[i for i in ys],
            z=[i for i in zs],
            mode="markers+lines",  # Include both markers and lines
            name=f"Axon",
            marker=dict(
                sizemode="diameter",
                size=1,  # Set the size of scatter points for the axon
                color=c,  # Set the color of scatter points for the axon
                line=dict(
                    color="rgba(0, 0, 0, 0.6)",  # Set color to semi-transparent black
                    width=2  # Set the width of the lines
                )
            ),
            line=dict(
                color="rgba(0, 0, 0, 0.6)",  # Set color to semi-transparent black
                width=2  # Set the width of the lines
            )
        )

        layout = go.Layout(
            scene=dict(
                xaxis=dict(title='X [mm]'),
                yaxis=dict(title='Y [mm]'),
                zaxis=dict(title='Z [mm]')
            )
        )
        scatters.append(scatter)
    print(len(scatters))
    # Create the figure
    fig = go.Figure(data=scatters, layout=layout)

    # Show the figure
    fig.show()


def read_bin_file_2d(file_path):
    traj_part = np.fromfile(file_path, dtype="float32")
    chunks =1
    nbr_steps=64300
    valid_traj_part = traj_part[np.isfinite(traj_part)]  # Ensure all entries are valid floats
    
    alls_xs = []
    alls_ys = []
    for i in range(chunks):
        start = (nbr_steps+1) * i
        end = nbr_steps * (i + 1)-1
        
        chunk = valid_traj_part[3*start:3*end]
        xs, ys = chunk[0::3], chunk[1::3]
        alls_xs.extend(xs)
        alls_ys.extend(ys)
    
   # 2d plot of all_xs and all_ys
    fig, ax = plt.subplots()
    ax.plot(alls_xs, alls_ys, 'o', color='black', markersize=1)
    ax.set(xlabel='X [mm]', ylabel='Y [mm]',
           title='2D plot of all X and Y coordinates')
    ax.grid()
    plt.show()


file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/test_rep_00_0.traj"
read_bin_file(file)