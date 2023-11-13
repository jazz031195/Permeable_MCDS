import re
from pathlib import Path
import os 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

DWI_folder = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/no_funnel/overlap_4/n1/")
df         = pd.DataFrame(columns = ["Soma begin", "Soma end", "Dendrites begin","Dendrites end", "rep"], dtype="object")

for file in os.listdir(DWI_folder):
    if "txt" in file and not "neurons_list" in file:
        N = float(file.split('_')[1])
        T = float(file.split('_')[3])
        if (N == 50000) & (T == 15000) & ("count_walker" in file):
            if "_rep_00" in file:
                rep = 1
            elif "_rep_01" in file:
                rep = 2
            elif "_rep_02" in file:
                rep = 3
            elif "_rep_03" in file:
                rep = 4
            else:
                rep = 0
            f = open(DWI_folder / file, "r")
            # Find numbers corresponding to the number of walkers starting in soma, ending in soma, ...
            numbers = re.findall(r'\d+', f.read())
            list_ = numbers[:4]
            list_.append(rep)
            df.loc[len(df)] = list_
        N = 50000

df["Soma begin"]      = df["Soma begin"].astype(int)
df["Soma end"]        = df["Soma end"].astype(int)
df["Dendrites begin"] = df["Dendrites begin"].astype(int)
df["Dendrites end"]   = df["Dendrites end"].astype(int)

sum = df.groupby(['rep']).sum() 

MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)   # fontsize of the figure title

# Calculate the percentage of walkers that start in soma or in dendrites
perc  = sum / N * 100
means = perc.mean().values
stds  = perc.std().values

fig, ax = plt.subplots()
ind   = np.arange(len(means) / 2)  # the x locations for the groups
width = 0.35  # the width of the bars
rects1  = ax.bar(ind - width/2, 
                 [means[0] / (means[0] + means[2]) * 100, means[2]/ (means[1] + means[3])* 100], 
                 width, 
                 yerr=[stds[0], stds[2]],
                 label='Begin')
rects2  = ax.bar(ind + width/2, 
                 [means[1] / (means[0] + means[2]) * 100, means[3]/ (means[1] + means[3])* 100], 
                 width, 
                 yerr=[stds[1], stds[3]],
                 label='End')

ax.set_ylabel('% water molecules')
ax.set_xticks(ind)
ax.set_xticklabels(('Soma', 'Dendrites'))
ax.legend()
plt.show()

print("T : ", T, ", step length :", np.sqrt(6 * 2.5 * 65 / T), ", decrease %", means[0] / (means[0] + means[2])* 100 - means[1] / (means[0] + means[2])* 100)

# Calculate the walkers' density in soma and dendrites 
r_soma          = 10e-6 # in [m]
volume_neurites = 8782.71 # in [um³] (3 branching)
volume_soma     = 4/3 * np.pi * r_soma**3 # in [m³]
volume_soma     = volume_soma * 1e18 # in [um³]
volume_neuron   = volume_neurites + volume_soma

means = sum.mean().values
stds  = sum.std().values

means[[0, 1]] = means[[0, 1]] / volume_soma
means[[2, 3]] = means[[2, 3]] / volume_neurites
stds[[0, 1]]  = stds[[0, 1]]  / volume_soma
stds[[2, 3]]  = stds[[2, 3]]  / volume_neurites
# 4.18879e-06 8.78271e-06
print(means, volume_soma, volume_neurites)

fig, ax = plt.subplots(figsize=(7, 7))
rects1  = ax.bar(ind - width/2, [means[0], means[2]], width, 
                 yerr=[stds[0], stds[2]],
                 label='Begin')
rects2  = ax.bar(ind + width/2, [means[1], means[3]], width, 
                 yerr=[stds[1], stds[3]],
                 label='End')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Water molecules density [part/um³]')
ax.set_xticks(ind)
ax.set_xticklabels(('Soma', 'Dendrites'))
# Place legend outside the plot at the top middle
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
plt.show()
