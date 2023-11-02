import matplotlib.pyplot as plt
import random
import numpy as np
import os
import plotly.graph_objs as go
wd = os.getcwd()
plot_3d = True
plot_traj = False
projection = False
z_slice = [0.02, 0.04, 0.06, 0.08]
max_lim = 0.75
min_lim = 0.25
neuron = "n1"
neuron_file = f'instructions/neurons/_rep_03_0_neurons_list_5.txt'
min_ = np.ones(3)
max_ = np.zeros(3)
with open(neuron_file) as f:
    lines = f.readlines()
    if plot_3d:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('X [mm]')
        ax.set_ylabel('Y [mm]')
        ax.set_zlabel('Z [mm]')
    elif projection:
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    else:
        fig, axs = plt.subplots(1, len(z_slice), figsize=(15, 5))
        for k in range(len(z_slice)):
            axs[k].set_xlabel("X [mm]")
            axs[k].set_ylabel("Y [mm]")
    for i in range(len(lines)):
        for j, slice_ in enumerate(z_slice):
            coords = lines[i].split(' ')
            # If it is the soma, plot it in any case
            if "Soma" in coords[0]:
                coords = lines[i-1].split(' ')
                coords = [float(coord) for coord in coords]
                # draw sphere
                u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
                y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
                z = np.cos(v)*float(coords[3]) + float(coords[2])
                r = float(coords[3])
                if plot_3d:
                    ax.plot_wireframe(x, y, z, color="cornflowerblue")
                elif projection:
                    axs[0].plot(x, y, color="cornflowerblue")
                    axs[1].plot(x, z, color="cornflowerblue")
                    axs[2].plot(z, y, color="cornflowerblue")
                else:
                    if ((coords[2] - coords[3]) < slice_) and ((coords[2] + coords[3]) > slice_):
                        idx = ((x - coords[0])**2 + (y - coords[1])**2 < coords[3]**2 - (slice_ - coords[2])**2)
                        x = x[idx]
                        y = y[idx]
                        x = np.squeeze(x.reshape(1, -1))
                        y = np.squeeze(y.reshape(1, -1))
                        z = np.ones_like(x)*slice_
                        if z[0] == slice_:
                            axs[j].plot(x, y, '.',color="darkgray")
                    else:
                        z = np.squeeze(z.reshape(1, -1))
            elif len(coords) > 2 and len(coords) < 6:
                coords = [float(coord) for coord in coords]
                if coords[0] < min_[0]:
                    min_[0] = coords[0]
                if coords[1] < min_[1]:
                    min_[1] = coords[1]
                if coords[2] < min_[2]:
                    min_[2] = coords[2]
                if coords[0] > max_[0]:
                    max_[0] = coords[0]
                if coords[1] > max_[1]:
                    max_[1] = coords[1]
                if coords[2] > max_[2]:
                    max_[2] = coords[2]

                # draw sphere
                u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
                y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
                z = np.cos(v)*float(coords[3]) + float(coords[2])
                r = float(coords[3])

                if i % 16 == 0:
                    
                    if plot_3d:
                        ax.plot_wireframe(x, y, z, color="cornflowerblue")
                    elif projection:
                        axs[0].plot(x, y, color="cornflowerblue")
                        axs[1].plot(x, z, color="cornflowerblue")
                        axs[2].plot(z, y, color="cornflowerblue")
                    else:
                        if ((coords[2] - coords[3]) < slice_) and ((coords[2] + coords[3]) > slice_):
                            idx = ((x - coords[0])**2 + (y - coords[1])**2 < coords[3]**2 - (slice_ - coords[2])**2)
                            x = x[idx]
                            y = y[idx]
                            x = np.squeeze(x.reshape(1, -1))
                            y = np.squeeze(y.reshape(1, -1))
                            z = np.ones_like(x)*slice_
                            if z[0] == slice_:
                                axs[j].plot(x, y, color="darkgray")
                        else:
                            z = np.squeeze(z.reshape(1, -1))
            elif len(coords) == 6:
                coords = [float(coord) for coord in coords]
                xmin = coords[0]
                xmax = coords[1]
                ymin = coords[2]
                ymax = coords[3]
                zmin = coords[4]
                zmax = coords[5]

                if plot_3d:
                    ax.plot(xmin, ymin, zmin, color="cornflowerblue")
                    ax.plot(xmax, ymax, zmax, color="cornflowerblue")
                elif projection:
                    axs[0].plot(xmin, ymin, color="cornflowerblue")
                    axs[0].plot(xmax, ymax, color="cornflowerblue")
                    axs[1].plot(xmin, zmin, color="cornflowerblue")
                    axs[1].plot(xmax, zmax, color="cornflowerblue")
                    axs[2].plot(zmin, ymin, color="cornflowerblue")
                    axs[2].plot(zmax, ymax, color="cornflowerblue")
    print(min_ - max_)
    if not plot_3d:
        distance_from_borders = 0.0007
    if plot_traj:
        with open(traj_file) as f:
            lines = f.readlines()
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
            if not plot_3d:
                axs[0].plot(xp, yp, 'b.', markersize=1)
                axs[1].plot(xp, zp, 'b.', markersize=1)
                axs[2].plot(zp, yp, 'b.', markersize=1)
                axs[0].set_xlim([0, 0.1])
                axs[0].set_ylim([0, 0.1])
                axs[0].set_xlabel("X [mm]")
                axs[0].set_ylabel("Y [mm]")
                axs[0].axvline(distance_from_borders, 0, 1)
                axs[0].axvline(0.1-distance_from_borders, ymin=0, ymax=1)
                axs[0].axhline(distance_from_borders, xmin=0, xmax=1)
                axs[0].axhline(y=0.1-distance_from_borders, xmin=0, xmax=1)
                axs[1].set_xlim([0, 0.1])
                axs[1].set_ylim([0, 0.1])
                axs[1].set_xlabel("X [mm]")
                axs[1].set_ylabel("Z [mm]")
                axs[1].axvline(distance_from_borders, 0, 1)
                axs[1].axvline(0.1-distance_from_borders, ymin=0, ymax=1)
                axs[1].axhline(distance_from_borders, xmin=0, xmax=1)
                axs[1].axhline(y=0.1-distance_from_borders, xmin=0, xmax=1)
                axs[2].set_xlim([0, 0.1])
                axs[2].set_ylim([0, 0.1])
                axs[2].set_xlabel("Z [mm]")
                axs[2].set_ylabel("Y [mm]")
                axs[2].axvline(distance_from_borders, 0, 1)
                axs[2].axvline(0.1-distance_from_borders, ymin=0, ymax=1)
                axs[2].axhline(distance_from_borders, xmin=0, xmax=1)
                axs[2].axhline(y=0.1-distance_from_borders, xmin=0, xmax=1)
            else:
                ax.scatter(xp, yp, zp, color="cornflowerblue")
    if plot_3d:
        ax.set_xlim([min_lim, max_lim])
        ax.set_ylim([min_lim, max_lim])
        ax.set_zlim([min_lim, max_lim])
    elif projection:
        for k in range(3):
            axs[k].set_xlim([min_lim, max_lim])
            axs[k].set_ylim([min_lim, max_lim])
    else:
        for k in range(len(z_slice)):
            axs[k].set_xlim([min_lim, max_lim])
            axs[k].set_ylim([min_lim, max_lim])
    plt.show()