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
from pathlib import Path



MEDIUM_SIZE = 17
BIGGER_SIZE = 19

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# # Sizes of the neuron file, in Mo
# sizes = [6.6e-3, 1.6, 2.1, 4.1, 6.1, 8.2]
# plt.bar(list(range(len(sizes))), sizes)
# plt.ylabel("Neuron file size [Mo]")
# plt.xticks(list(range(len(sizes))), ["SWC", "Mesh 0.01", "Mesh 0.01325", "Mesh 0.0255", "Mesh 0.0375", "Mesh 0.05"])
# plt.show()

# Time overlapping spheres
folder = Path("results/ISMRM24/overlaps/")
df = pd.DataFrame(columns=['overlap', 'time [min.]'])
for overlap in os.listdir(folder):
        if (os.path.isdir(folder / overlap)) and ("overlap" in overlap):
            for neuron in os.listdir(folder / overlap):
                if os.path.isdir(folder / overlap / neuron):
                    # Iterate through the files in the folder
                    for subdir in os.listdir(folder / overlap / neuron):
                        if os.path.isdir(folder / overlap / neuron / subdir):
                            for filename in os.listdir(folder / overlap / neuron / subdir):
                            
                                # Read the simulation_info.txt to have crossings information
                                if "simu" in filename:
                                    path = folder / overlap / neuron / subdir / filename

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
for mesh in os.listdir(folder):
        if (os.path.isdir(folder / mesh)) and ("mesh" in mesh):
            for neuron in os.listdir(folder / mesh):
                if os.path.isdir(folder / mesh / neuron):
                    # Iterate through the files in the folder
                    for subdir in os.listdir(folder / mesh / neuron):
                        if os.path.isdir(folder / mesh / neuron / subdir):
                            for filename in os.listdir(folder / mesh / neuron / subdir):
                                if "sim" in filename:
                                    path = folder / mesh / neuron / subdir / filename
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
sns.boxplot(data=df, x='overlap', y='time [min.]', order=['overlap2', 'overlap4', 'overlap8', 'overlap16', 'overlap32', 'mesh_001', 'mesh_001325', 'mesh_002550', 'mesh_003775', 'mesh_005'])
plt.xticks(list(range(len(df.overlap.unique()))), ["R/2", "R/4", "R/8", "R/16", "R/32", "Mesh 1%", "Mesh 1.325%", "Mesh 2.55%", "Mesh 3.775%", "Mesh 5%"])
plt.show()

