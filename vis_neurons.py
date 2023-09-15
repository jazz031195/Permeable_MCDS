import matplotlib.pyplot as plt
import random
import numpy as np
import os
import plotly.graph_objs as go

wd         = os.getcwd()
plot_3d    = False
plot_traj  = True
projection = True
z_slice = [0.02, 0.04, 0.06, 0.08]
# position = np.array([0.044359383652472821, 0.073838872992544963, 0.053730584153105686])
# position2 = np.array([0.0450215294826024, 0.073325599802571709, 0.053024812566757319])
max_lim = 0.75
min_lim = 0.25

neuron_file = wd + f'/instructions/neurons/_rep_10_neurons_list.txt'
traj_file = wd + '/instructions/neurons/_rep_16_0.traj.txt'
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
    # lines_plot = []
    for i in range(len(lines)):
        for j, slice_ in enumerate(z_slice):
            coords = lines[i].split(' ')
            # If it is the soma, plot it in any case
            if "Soma" in coords[0]:
                coords = lines[i-1].split(' ')
                coords = [float(coord) for coord in coords]
                # if max_lim < coords[0]:
                #     max_lim = coords[0]
                # draw sphere
                u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
                y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
                z = np.cos(v)*float(coords[3]) + float(coords[2])
                r = float(coords[3])
                # if np.linalg.norm(coords[:3] - position) <= r:
                #     print("soma", coords)
                if plot_3d:
                    ax.plot_wireframe(x, y, z, color="r")
                    # # Creating the plot
                    # line_marker = dict(color='red', width=2)
                    # for l, m, n in zip(x, y, z):
                    #     lines_plot.append(go.Scatter3d(x=l, y=m, z=n, mode='lines', line=line_marker))
                elif projection:
                    axs[0].plot(x, y, color="r")
                    axs[1].plot(x, z, color="r")
                    axs[2].plot(z, y, color="r")
                    # axs[0].plot(position[0], position[1], ".b", markersize=5)
                    # axs[1].plot(position[0], position[2], ".b", markersize=5)
                    # axs[2].plot(position[2], position[1], ".b", markersize=5)
                    # axs[0].plot(position2[0], position2[1], ".g", markersize=5)
                    # axs[1].plot(position2[0], position2[2], ".g", markersize=5)
                    # axs[2].plot(position2[2], position2[1], ".g", markersize=5)
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
                # if max_lim < coords[0]:
                #     max_lim = coords[0]
                # Plot only one sphere out of four for the dendrites (otherwise, too expensive)
                # a = random.randint(1, 4)
                # draw sphere
                u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
                y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
                z = np.cos(v)*float(coords[3]) + float(coords[2])
                r = float(coords[3])
                # if np.linalg.norm(coords[:3] - position) <= r:
                #     print("d", coords)
                if i % 8 == 0:
                    
                    if plot_3d:
                        ax.plot_wireframe(x, y, z, color="r")
                        # # Creating the plot
                        # line_marker = dict(color='blue', width=2)
                        # for l, m, n in zip(x, y, z):
                        #     lines_plot.append(go.Scatter3d(x=l, y=m, z=n, mode='lines', line=line_marker))
                    elif projection:
                        axs[0].plot(x, y, color="r")
                        axs[1].plot(x, z, color="r")
                        axs[2].plot(z, y, color="r")
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
                # print(coords)

                if plot_3d:
                    ax.plot(xmin, ymin, zmin, color="g")
                    ax.plot(xmax, ymax, zmax, color="g")
                    # fig.savefig("/home/localadmin/Images/ESMRMB_23/" + neuron, format="pdf")
                        # # Creating the plot
                        # line_marker = dict(color='blue', width=2)
                        # for l, m, n in zip(x, y, z):
                        #     lines_plot.append(go.Scatter3d(x=l, y=m, z=n, mode='lines', line=line_marker))
                elif projection:
                    axs[0].plot(xmin, ymin, color="g")
                    axs[0].plot(xmax, ymax, color="g")
                    axs[1].plot(xmin, zmin, color="g")
                    axs[1].plot(xmax, zmax, color="g")
                    axs[2].plot(zmin, ymin, color="g")
                    axs[2].plot(zmax, ymax, color="g")
    print(min_ - max_)
    if not plot_3d:
        distance_from_borders = 0.0007
    # if plot_3d:
    #     layout = go.Layout(
    #         title='Wireframe Plot',
    #         scene=dict(
    #             xaxis=dict(
    #                 gridcolor='rgb(255, 255, 255)',
    #                 zerolinecolor='rgb(255, 255, 255)',
    #                 showbackground=True,
    #                 backgroundcolor='rgb(230, 230,230)'
    #             ),
    #             yaxis=dict(
    #                 gridcolor='rgb(255, 255, 255)',
    #                 zerolinecolor='rgb(255, 255, 255)',
    #                 showbackground=True,
    #                 backgroundcolor='rgb(230, 230,230)'
    #             ),
    #             zaxis=dict(
    #                 gridcolor='rgb(255, 255, 255)',
    #                 zerolinecolor='rgb(255, 255, 255)',
    #                 showbackground=True,
    #                 backgroundcolor='rgb(230, 230,230)'
    #             )
    #         ),
    #         showlegend=False,
    #     )
    #     fig = go.Figure(data=lines_plot, layout=layout)
    #     fig.show()
    if plot_traj:
        with open(traj_file) as f:
            lines = f.readlines()
            xp = []
            yp = []
            zp = []
            for i in range(int(len(lines))):
                if i % 4 == 0:
                    xp.append(float(lines[i]))
                elif i % 4 == 1:
                    yp.append(float(lines[i]))
                elif i % 4 == 2:
                    zp.append(float(lines[i]))
            if not plot_3d:
                axs[0].plot(xp, yp, 'g.', markersize=1)
                axs[1].plot(xp, zp, 'g.', markersize=1)
                axs[2].plot(zp, yp, 'g.', markersize=1)
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
                ax.scatter(xp, yp, zp, color="g")
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