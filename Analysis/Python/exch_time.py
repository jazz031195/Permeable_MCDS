import numpy as np
import os 
import matplotlib.pyplot as plt
import pandas as pd

ini_pos = "/home/localadmin/Bureau/simu_templates/simu_templates/diff_times_single_neuron/exch/overlap4/soma_dendrites_ex/n1/ini_pos_file.txt"
pos = pd.read_csv(ini_pos, header=None, sep=' ')

# fig, ax = plt.subplots(1,1)
# ax.scatter(pos.iloc[:, 0], pos.iloc[:, 1], pos.iloc[:, 2])
# plt.show()

folders = ["/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/soma_dendrites_time/exch/overlap4/n1",
           "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/soma_dendrites_time/exch_funnel/overlap4/n1"]
T = 15000
dt = 87/T # [ms]

for folder in folders:
    if "funnel" in folder:
        funnel = True
    else:
        funnel = False

    filenames = os.listdir(folder)
    filenames = [file for file in filenames if "idx" in file]
    
    d = {'soma_dendrites_t': [], 'dendrites_soma_t': []}
    for file in filenames:

        data = np.fromfile(f"{folder}/{file}", dtype="float32")
        data = data.reshape((-1, 3))
        for walker in range(int(data.shape[0] / T)):
            data_walker = data[walker*T:(walker+1)*T, :] 

            idx_soma = data_walker[:, 0]

            changes = np.diff(idx_soma) != 0
            change_indices = np.where(changes)[0] + 1

            segments = np.split(idx_soma, change_indices)

            if len(segments) > 1:
                value = segments[0][0]
                duration = len(segments[0])
                if value == -1:
                    d['dendrites_soma_t'].append(duration*dt)
                else:
                    d['soma_dendrites_t'].append(duration*dt)
            else:
                value = segments[0][0]
                if value == -1:
                    d['dendrites_soma_t'].append(1000)
                else:
                    d['soma_dendrites_t'].append(1000)


    fig, ax = plt.subplots(1, 2, sharey=True, tight_layout=True)
    ax[0].hist(d['dendrites_soma_t'], bins=50)
    ax[0].set_title("Dendrites to soma")
    ax[0].set_xlabel("Diffusion time [ms]")
    ax[0].set_ylabel("Count")
    ax[1].hist(d['soma_dendrites_t'], bins=50)
    ax[1].set_title("Soma to dendrites")
    ax[1].set_xlabel("Diffusion time [ms]")
    ax[1].set_ylabel("Count")
    if funnel:
        fig.suptitle("Funnel")
        fig.savefig(f"/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/soma_dendrites_time/funnel.jpg")
    else:
        fig.suptitle("No funnel")
        fig.savefig(f"/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/soma_dendrites_time/no_funnel.jpg")
    
    fig, ax = plt.subplots(1, 2, sharey=True, tight_layout=True)
    ax[0].hist(d['dendrites_soma_t'], bins=50, range=(0, 90))
    ax[0].set_title("Dendrites to soma")
    ax[0].set_xlabel("Diffusion time [ms]")
    ax[0].set_ylabel("Count")
    ax[1].hist(d['soma_dendrites_t'], bins=50, range=(0, 90))
    ax[1].set_title("Soma to dendrites")
    ax[1].set_xlabel("Diffusion time [ms]")
    ax[1].set_ylabel("Count")
    if funnel:
        fig.suptitle("Funnel")
        fig.savefig(f"/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/soma_dendrites_time/funnel_zoomed.jpg")
    else:
        fig.suptitle("No funnel")
        fig.savefig(f"/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/soma_dendrites_time/no_funnel_zoomed.jpg")

