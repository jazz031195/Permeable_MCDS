import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import plotly.graph_objects as go
import plotly.colors as colors

def read_bin_file(file_path):
    nbr_steps = 11171

    i = 200
    start = nbr_steps*i
    limit = nbr_steps*(i+1)-1

    traj_part = np.fromfile(file_path, dtype="float32")

    xs =[]
    ys =[]
    zs =[]  
    e = 0  

    for t in traj_part:
        if t != " ":
            if e%3 == 0:
                xs.append(t)
            elif e%3 == 1:
                ys.append(t)
            elif e%3 == 2:
                zs.append(t)

            e= e+ 1

    # Example scatter plot
    colours = colors.qualitative.Plotly[:10]
    c = colours[0]  # Assuming e is defined somewhere in your code

    scatter = go.Scatter3d(
        x=[i for i in xs[start:limit]],
        y=[i for i in ys[start:limit]],
        z=[i for i in zs[start:limit]],
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

    # Create the figure
    fig = go.Figure(data=[scatter], layout=layout)

    # Show the figure
    fig.show()


file = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/instructions/axons_vs_cylinders/data/cylinders_icvf_0.70_cap_24_vox_50_straight_intra_1.traj"
read_bin_file(file)