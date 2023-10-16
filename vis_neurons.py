import matplotlib.pyplot as plt
import random
import numpy as np
import os
import plotly.graph_objs as go

def read_swc(neuron_file_name):

    f = open ( neuron_file_name, 'r' )
    lines = f.readlines()
    point_dict = {}
    segments = []
    num_total_segments = 0
    num_segs_limit = 0
    first_count = 0
    soma = []
    for i, l in enumerate(lines):
        l = l.strip()
        if len(l) > 5:
            if l[0] != "#":
                fields = l.split()
                point_dict[fields[0]] = fields
                if first_count == 0 and fields[1] == '1':
                    soma = [float(u) for u in fields[2:-1]]
                first_count = first_count + 1
    point_keys = sorted ( [ k for k in point_dict.keys() ] )

    num_lines_in_file = len(point_keys)
    num_nodes_in_file = len(point_keys)

    # Next create the list of segments - one for each child that has a parent
    for k in point_keys:
        child_fields = point_dict[k]
        if child_fields[6] in point_keys:
            # This point has a parent, so make a segment from parent to child
            parent_fields = point_dict[child_fields[6]]
            px = float(parent_fields[2])
            py = float(parent_fields[3])
            pz = float(parent_fields[4])
            pr = float(parent_fields[5])
            cx = float(child_fields[2])
            cy = float(child_fields[3])
            cz = float(child_fields[4])
            cr = float(child_fields[5])
            segments = segments + [ [ [px, py, pz, pr], [cx, cy, cz, cr] ] ]
            num_total_segments += 1

    if num_segs_limit > 0:
        # Limit the number of segments
        segments = segments[0:self.num_segs_limit]
        num_total_segments = len(segments)

    num_segments_in_file = num_total_segments
    return soma, segments

def plot_skeleton_swc(soma, segment, projection):

    if not projection:
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    else:
        fig = plt.figure()
        axs = fig.add_subplot(111, projection='3d')

    lim_max = 1
    for segm in segments:
        #axis and radius
        p0 = np.array(segm[0][:-1])
        p1 = np.array(segm[1][:-1])
    
        #plot axis
        if not projection:
            axs[0].plot([p0[0], p1[0]], [p0[1], p1[1]], color = 'black', linewidth=1)
            axs[1].plot([p0[0], p1[0]], [p0[2], p1[2]], color = 'black', linewidth=1)
            axs[2].plot([p0[2], p1[2]], [p0[1], p1[1]], color = 'black', linewidth=1) 
        else:
            #plot axis
            axs.plot(*zip(p0, p1), color = 'black', linewidth=1)


    if len(soma) > 0:
        # Define the radius of the sphere
        radius = soma[-1]
        # Define the center of the sphere
        center = np.array(soma[:-1]) + np.array([0, 0, 0])
        # Create a sphere
        # draw sphere
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = np.cos(u)*np.sin(v)*radius + center[0]
        y = np.sin(u)*np.sin(v)*radius + center[1]
        z = np.cos(v)*radius + center[2]
        
        if not projection:
            # Add the sphere to the plot
            axs[0].plot(x, y, color="r")
            axs[1].plot(x, z, color="r")
            axs[2].plot(z, y, color="r")
        else:
            # Add the sphere to the plot
            axs.plot_surface(x, y, z)


    if not projection:
        axs[0].set_xlim([0, 1])
        axs[0].set_ylim([0, 1])
        axs[1].set_xlim([0, 1])
        axs[1].set_ylim([0, 1])
        axs[2].set_xlim([0, 1])
        axs[2].set_ylim([0, 1])

    return fig, axs

wd         = os.getcwd()
plot_3d    = False
plot_traj  = True

max_lim = 0.75
min_lim = 0.25

neuron_file = wd + f'/instructions/neurons/_neurons_list.swc'
traj_file   = wd + '/instructions/neurons/_rep_00_0.traj.txt'

# Load the SWC file
soma, segments = read_swc(neuron_file)
fig, axs = plot_skeleton_swc(soma, segments, plot_3d)

min_ = np.ones(3)
max_ = np.zeros(3)
distance_from_borders = 0.0007
with open(neuron_file) as f:
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
                axs[0].set_xlabel("X [mm]")
                axs[0].set_ylabel("Y [mm]")
                axs[0].axvline(distance_from_borders, 0, 1)
                axs[0].axvline(1-distance_from_borders, ymin=0, ymax=1)
                axs[0].axhline(distance_from_borders, xmin=0, xmax=1)
                axs[0].axhline(y=1-distance_from_borders, xmin=0, xmax=1)
                axs[1].set_xlabel("X [mm]")
                axs[1].set_ylabel("Z [mm]")
                axs[1].axvline(distance_from_borders, 0, 1)
                axs[1].axvline(1-distance_from_borders, ymin=0, ymax=1)
                axs[1].axhline(distance_from_borders, xmin=0, xmax=1)
                axs[1].axhline(y=1-distance_from_borders, xmin=0, xmax=1)
                axs[2].set_xlabel("Z [mm]")
                axs[2].set_ylabel("Y [mm]")
                axs[2].axvline(distance_from_borders, 0, 1)
                axs[2].axvline(1-distance_from_borders, ymin=0, ymax=1)
                axs[2].axhline(distance_from_borders, xmin=0, xmax=1)
                axs[2].axhline(y=1-distance_from_borders, xmin=0, xmax=1)
            else:
                axs.scatter(xp, yp, zp, color="g")
    
    if not plot_3d:
        for k in range(3):
            axs[k].set_xlim([min_lim, max_lim])
            axs[k].set_ylim([min_lim, max_lim])
    plt.show()