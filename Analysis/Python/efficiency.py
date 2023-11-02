import matplotlib.pyplot as plt
import numpy as np
import os 
import time
from datetime import datetime
from datetime import timedelta
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# sizes = [6.6e-3, 3.3, 6.7, 17.6, 35.9, 54.2, 227.5]
# plt.bar(list(range(len(sizes))), sizes)
# plt.ylabel("Neuron file size [Mo]")
# plt.xticks(list(range(len(sizes))), ["Spheres", "Mesh 0.05", "Mesh 0.10", "Mesh 0.25", "Mesh 0.50", "Mesh 0.75", "Mesh 1"])
# plt.show()

folder = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/no_funnel/"

df = pd.DataFrame(columns=['overlap', 'time [min.]'])
for overlap in os.listdir(folder):
    print(overlap)
    for neuron in os.listdir(folder + overlap):
        if neuron == "n0":
            for rep in os.listdir(folder + overlap + "/" + neuron):
                if "sim" in rep and "T_15000" in rep:
                    path = folder + overlap + "/" + neuron + "/" + rep
                    print(path)
                    last = datetime.fromtimestamp(os.path.getmtime(path)).strftime('%H:%M:%S')

                    import re

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

                    time1 = datetime.strptime(last, "%H:%M:%S")
                    time2 = datetime.strptime(creation_time, "%H:%M:%S")
                    diff = time1 - time2
                    total_seconds = diff.total_seconds()
                    if total_seconds/60 > 0 :
                        df.loc[len(df)] = {'overlap': overlap, 'time [min.]': total_seconds/60}

folder = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/exch/"
for mesh in os.listdir(folder):
    if "mesh" in mesh:
        print(mesh)
        for neuron in os.listdir(folder + mesh):
            for rep in os.listdir(folder + mesh + "/" + neuron):
                if "sim" in rep and "T_15000" in rep:
                    path = folder + mesh + "/" + neuron + "/" + rep
                    print(path)
                    last = datetime.fromtimestamp(os.path.getmtime(path)).strftime('%H:%M:%S')

                    import re

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

                    time1 = datetime.strptime(last, "%H:%M:%S")
                    time2 = datetime.strptime(creation_time, "%H:%M:%S")
                    diff = time1 - time2
                    total_seconds = diff.total_seconds()
                    if total_seconds/60 > 0 and total_seconds/60 < 80:
                        df.loc[len(df)] = {'overlap': mesh, 'time [min.]': total_seconds/60}
print(df.groupby(['overlap'])['time [min.]'].mean().reset_index())
print(df.groupby(['overlap'])['time [min.]'].std().reset_index())
sns.boxplot(data=df, x='overlap', y='time [min.]', order=['overlap_1', 'overlap_2', 'overlap_4', 'overlap_8', 'overlap_16', 'overlap_32', 'mesh_005', 'mesh_010', 'mesh_025', 'mesh_050', 'mesh_075', 'mesh_100'])
plt.show()

