import re
from pathlib import Path
import os 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

# branching = 'branching'
# DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/exchange")
# df = pd.DataFrame(columns = ["Soma begin", "Soma end", "Dendrites begin","Dendrites end", "neuron"], dtype="object")
# for neuron in os.listdir(DWI_folder):
#     for file in os.listdir(DWI_folder / neuron / "perc_exchange"):
#         N = float(file.split('_')[1])
#         T = float(file.split('_')[3])
#         f = open(DWI_folder / neuron / f"perc_exchange/{file}", "r")
#         numbers = re.findall(r'\d+', f.read())
#         list_ = numbers[:4]
#         list_.append(neuron)
#         df.loc[len(df)] = list_

branching = "branching"
DWI_folder = Path("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/funnel/overlap_4/n1/")
df = pd.DataFrame(columns = ["Soma begin", "Soma end", "Dendrites begin","Dendrites end", "rep"], dtype="object")
for file in os.listdir(DWI_folder):
    number = 0
    if "txt" in file and not "neurons_list" in file:
        N = float(file.split('_')[1])
        T = float(file.split('_')[3])
        if (N == 50000) & ("count_walker" in file):
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
            numbers = re.findall(r'\d+', f.read())
            list_ = numbers[:4]
            list_.append(rep)
            df.loc[len(df)] = list_
        N = 50000

df["Soma begin"]      = df["Soma begin"].astype(int)
df["Soma end"]        = df["Soma end"].astype(int)
df["Dendrites begin"] = df["Dendrites begin"].astype(int)
df["Dendrites end"]   = df["Dendrites end"].astype(int)

# sum = df.groupby(['neuron']).sum() / 5
sum = df.groupby(['rep']).sum() 
# sum = df.sum() 

MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


perc  = sum / N * 100
# means = sum.values
means = perc.mean().values
# stds  = np.zeros_like(means)
stds  = perc.std().values
# plt.errorbar(["Soma begin", "Soma end", "Dendrites begin", "Dendrites end"], perc.mean().values, yerr=perc.std().values, fmt='.')
# plt.ylabel('% walker')
# plt.show()

ind = np.arange(len(means) / 2)  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, [means[0] / (means[0] + means[2]) * 100, means[2]/ (means[1] + means[3])* 100], width, yerr=[stds[0], stds[2]],
                label='Begin')
rects2 = ax.bar(ind + width/2, [means[1] / (means[0] + means[2]) * 100, means[3]/ (means[1] + means[3])* 100], width, yerr=[stds[1], stds[3]],
                label='End')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('% water molecules')
ax.set_xticks(ind)
ax.set_xticklabels(('Soma', 'Dendrites'))
ax.legend()
plt.show()
print("T : ", T, ", step length :", np.sqrt(6 * 2.5 * 65 / T), ", decrease %", means[0] / (means[0] + means[2])* 100 - means[1] / (means[0] + means[2])* 100)

nb_neurites = 20
r_soma      = 10e-6 #m
r_neurite   = 0.5e-6 #m
if branching == 'branching':
    nb_branching = 3
    l_neurite       = 80e-6 #m
    volume_neurites = nb_neurites * (2**nb_branching + 1) * np.pi*r_neurite**2*l_neurite # in m³
else:
    l_neurite       = 240e-6 # m
    volume_neurites = nb_neurites * np.pi*r_neurite**2*l_neurite # in m³
volume_soma     = 4/3 * np.pi * r_soma**3 # in m³
volume_soma     = volume_soma * 1e18
volume_neurites = volume_neurites * 1e18
volume_neurites = 8782.71
volume_neuron   = volume_neurites + volume_soma

# means = sum.values.astype('float') 
means = sum.mean().values
# stds  = np.zeros_like(means)
stds  = sum.std().values

means[[0, 1]] = means[[0, 1]] / volume_soma
means[[2, 3]] = means[[2, 3]] / volume_neurites
stds[[0, 1]]  = stds[[0, 1]]  / volume_soma
stds[[2, 3]]  = stds[[2, 3]]  / volume_neurites
# 4.18879e-06 8.78271e-06
# print(means, volume_soma, volume_neurites)

fig, ax = plt.subplots()
rects1  = ax.bar(ind - width/2, [means[0], means[2]], width, yerr=[stds[0], stds[2]],
                label='Begin')
rects2  = ax.bar(ind + width/2, [means[1], means[3]], width, yerr=[stds[1], stds[3]],
                label='End')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Water molecules density [part/um³]')
ax.set_xticks(ind)
ax.set_xticklabels(('Soma', 'Dendrites'))
ax.legend()
plt.show()

# decrease = np.array([[5000.0, 10.065652850524717],
#                      [20000.0, 7.7623999999999995],
#                      [40000.0, 5.972399999999997],
#                      [60000.0, 4.927199999999999],
#                      [80000.0, 4.30599999999999],
#                      [90500.0, 4.135500000000004]])

# model = LinearRegression()
# x = decrease[:, 0].reshape(-1, 1)
# y = decrease[:, 1].reshape(-1, 1)
# model.fit(x, y)
# r_sq = model.score(x, y)
# print(f"coefficient of determination: {r_sq}")
# x = np.linspace(5000, 100000, 1000)
# a = model.coef_[0][0]
# b = model.intercept_[0]
# y = np.squeeze(b + a * x)
# fig, ax = plt.subplots()
# ax.plot(decrease[:, 0], decrease[:, 1], '.')
# ax.plot(x, y)
# ax.set_xlabel("Number of timesteps")
# ax.set_ylabel("% walkers flow soma - dendrites")
# ax.annotate(f'y = {a:.1e}x + {b:.1e} \n R² = {r_sq:.2f}', xy =(50000, 9),
#                 xytext =(50000, 9))
# plt.show()

# decrease[:, 0] = np.sqrt(6 * 2.5 * 67 / decrease[:, 0])
# model = LinearRegression()
# x = decrease[:, 0].reshape(-1, 1)
# y = decrease[:, 1].reshape(-1, 1)
# model.fit(x, y)
# r_sq = model.score(x, y)
# x = np.linspace(0.1, 0.5, 1000)
# a = model.coef_[0][0]
# b = model.intercept_[0]
# y = np.squeeze(b + a * x)
# fig, ax = plt.subplots()
# ax.plot(decrease[:, 0], decrease[:, 1], '.')
# ax.plot(x, y)
# ax.set_xlabel("Step length [um]")
# ax.set_ylabel("% walkers flow soma - dendrites")
# ax.annotate(f'y = {a:.1e}x + {b:.1e} \n R² = {r_sq:.2f}', xy =(0.1, 9),
#                 xytext =(0.1, 9))
# plt.show()