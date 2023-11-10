"""
File that calculates and plot the efficiency of the algorithm, in terms of memory load and running time
"""
import matplotlib.pyplot as plt
import os 
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import re


MEDIUM_SIZE = 17
BIGGER_SIZE = 19

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Sizes of the neuron file, in Mo
sizes = [6.6e-3, 0.74, 1.4, 2.9, 5.8, 11.5, 23.0, 3.3, 6.7, 17.6, 35.9, 54.2, 227.5]
plt.bar(list(range(len(sizes))), sizes)
plt.ylabel("Neuron file size [Mo]")
plt.xticks(list(range(len(sizes))), ["SWC", "R/1", "R/2", "R/4", "R/8", "R/16", "R/32", "Mesh 0.05", "Mesh 0.10", "Mesh 0.25", "Mesh 0.50", "Mesh 0.75", "Mesh 1"])
plt.show()

# Time overlapping spheres
folder = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/no_funnel/"
df = pd.DataFrame(columns=['overlap', 'time [min.]'])
for overlap in os.listdir(folder):
    print(overlap)
    for neuron in os.listdir(folder + overlap):
        if neuron == "n6":
            for rep in os.listdir(folder + overlap + "/" + neuron):
                if "sim" in rep and "T_15000" in rep:
                    path = folder + overlap + "/" + neuron + "/" + rep

                    # Find last modification time (time of run end)
                    last_modif_time = datetime.fromtimestamp(os.path.getmtime(path)).strftime('%H:%M:%S')
                    last_modif_time = datetime.strptime(last_modif_time, "%H:%M:%S")

                    # Find creation time in the simulation_info.txt
                    # Define a regular expression pattern to match the creation time line
                    pattern = r'Date and Time:.*?(\d{2}-\d{2}-\d{4} \(\d{2}:\d{2}:\d{2}\))'
                    # Read the text file
                    with open(path, 'r') as file:
                        content = file.read()
                    # Search for the pattern in the content
                    match = re.search(pattern, content)
                    # If a match is found, extract the creation time
                    creation_time = match.group(0)
                    creation_time = creation_time[-9:-1]
                    creation_time = datetime.strptime(creation_time, "%H:%M:%S")

                    run_time         = last_modif_time - creation_time
                    run_time_seconds = run_time.total_seconds()
                    if run_time_seconds/60 > 0 :
                        df.loc[len(df)] = {'overlap': overlap, 'time [min.]': run_time_seconds/60}

# Time mesh
folder = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/exch/"
for mesh in os.listdir(folder):
    if "mesh" in mesh:
        print(mesh)
        for neuron in os.listdir(folder + mesh):
            for rep in os.listdir(folder + mesh + "/" + neuron):
                if "sim" in rep and "T_15000" in rep:
                    path = folder + mesh + "/" + neuron + "/" + rep
                    print(path)
                    last_modif_time = datetime.fromtimestamp(os.path.getmtime(path)).strftime('%H:%M:%S')
                    last_modif_time = datetime.strptime(last_modif_time, "%H:%M:%S")

                    # Define a regular expression pattern to match the creation time line
                    pattern = r'Date and Time:.*?(\d{2}-\d{2}-\d{4} \(\d{2}:\d{2}:\d{2}\))'
                    # Read the text file
                    with open(path, 'r') as file:
                        content = file.read()
                    # Search for the pattern in the content
                    match = re.search(pattern, content)
                    # If a match is found, extract the creation time
                    creation_time = match.group(0)
                    creation_time = creation_time[-9:-1]
                    creation_time = datetime.strptime(creation_time, "%H:%M:%S")

                    run_time         = last_modif_time - creation_time
                    run_time_seconds = run_time.total_seconds()
                    if run_time_seconds/60 > 0 and run_time_seconds/60 < 80:
                        df.loc[len(df)] = {'overlap': mesh, 'time [min.]': run_time_seconds/60}

print(df.groupby(['overlap'])['time [min.]'].mean().reset_index())
print(df.groupby(['overlap'])['time [min.]'].std().reset_index())
sns.boxplot(data=df, x='overlap', y='time [min.]', order=['overlap_1', 'overlap_2', 'overlap_4', 'overlap_8', 'overlap_16', 'overlap_32', 'mesh_005', 'mesh_010', 'mesh_025', 'mesh_050', 'mesh_075', 'mesh_100'])
plt.xticks(list(range(len(df.overlap.unique()))), ["R/1", "R/2", "R/4", "R/8", "R/16", "R/32", "Mesh 0.05", "Mesh 0.10", "Mesh 0.25", "Mesh 0.50", "Mesh 0.75", "Mesh 1"])
plt.show()

