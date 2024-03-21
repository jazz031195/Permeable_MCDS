import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import pandas as pd
import plotly.graph_objects as go

def draw_neurons(neuron_ids):

    fig = go.Figure()

    for neuron_id in neuron_ids:
        df = pd.read_csv(f"/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/neuron_1.csv")
        # Filter out glial axons
        neurons = df # df[(df["neuron_id"] == 15) | (df["neuron_id"] == 18) | (df["neuron_id"] == 29)]

        fig.add_trace(go.Scatter3d(
                                    x=neurons["x"],
                                    y=neurons["y"],
                                    z=neurons["z"],
                                    type="scatter3d",
                                    mode="markers",
                                    marker=dict(
                                        sizemode="diameter",
                                        size=neurons["r"]*15000,
                                        color="blue",
                                        line=dict(
                                            color="rgba(0, 0, 0, 0)",
                                            width=0
                                        )
                                    )
                                )
                    )
    
        layout = go.Layout(
            scene=dict(
                xaxis=dict(title='X [mm]', range=[0, 0.1]),
                yaxis=dict(title='Y [mm]', range=[0, 0.1]),
                zaxis=dict(title='Z [mm]', range=[0, 0.1])
            )
        )

    # Show the figure
    fig.show()

def create_neurons_df(neuron_file):
    neuron_id = 1 
    df = pd.DataFrame(columns=["neuron_id", "x", "y", "z", "r"])
    
    with open(neuron_file) as f:
        lines = f.readlines()
        for i in range(len(lines)):
            coords = lines[i].split(' ')
            # If it is the soma, plot it in any case
            if ("Neuron" in lines[i].split(' ')[0]):
                if len(df.index) > 1:
                    df.to_csv(f"/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/neuron_{neuron_id}.csv")
                    df = pd.DataFrame(columns=["neuron_id", "x", "y", "z", "r"])
                neuron_id = neuron_id + 1
                coords = lines[i].split(' ')
                coords = [float(coord) for coord in coords]
                # if ((coords[0] < 0.075 and coords[0] > 0.025) and 
                #     (coords[1] < 0.075 and coords[1] > 0.025) and
                #     (coords[2] < 0.075 and coords[2] > 0.025)):
                d = {"neuron_id": neuron_id, 
                    "x": coords[0], 
                    "y": coords[1], 
                    "z": coords[2], 
                    "r": coords[3]}
                df_avg_data = pd.DataFrame(d, index=[neuron_id])
                df = pd.concat([df, df_avg_data])
            elif len(coords) > 2 and len(coords) < 6:
                coords = [float(coord) for coord in coords]
                # if ((coords[0] < 0.075 and coords[0] > 0.025) and 
                #     (coords[1] < 0.075 and coords[1] > 0.025) and
                #     (coords[2] < 0.075 and coords[2] > 0.025)):
                d = {"neuron_id": neuron_id, 
                    "x": coords[0], 
                    "y": coords[1], 
                    "z": coords[2], 
                    "r": coords[3]}
            
                df_avg_data = pd.DataFrame(d, index=[neuron_id])
                df = pd.concat([df, df_avg_data])

    return df

neuron_file = 'results/ISMRM24/branching/overlap8/nonbranching_560/n1/neurons_list.txt'
# df = create_neurons_df(neuron_file)
# df.to_csv("/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/packing/data.csv")
# print(df[df.r == 0.01])
draw_neurons([1])


# plot_3d    = True
# plot_traj  = True
# projection = False
# max_lim = 0.3
# min_lim = 0.5
# neuron = "n1"
# neuron_file = '/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/instructions/neurons/_rep_00_0_neurons_list.txt'
# traj_file   = ['/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/instructions/neurons/_rep_00_0.traj.txt',
#                '/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/instructions/neurons/_rep_03_0.traj.txt',
#                '/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/instructions/neurons/_rep_05_0.traj.txt']
# min_ = np.ones(3)
# max_ = np.zeros(3)
# i = 0
# with open(neuron_file) as f:
#     lines = f.readlines()
#     if plot_3d:
#         fig = plt.figure(figsize=(20, 20))
#         ax = fig.add_subplot(111, projection='3d')
#         # ax.set_xlabel('X [mm]')
#         # ax.set_ylabel('Y [mm]')
#         # ax.set_zlabel('Z [mm]')
#         ax.grid(False)
#         ax.xaxis.set_ticklabels([])
#         ax.yaxis.set_ticklabels([])
#         ax.zaxis.set_ticklabels([])
#         ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
#         ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
#         ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
#         ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#         ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#         ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#         ax.set_xticks([])
#         ax.set_yticks([])
#         ax.set_zticks([])
#         # for line in ax.xaxis.get_ticklines():
#         #     line.set_visible(False)
#         # for line in ax.yaxis.get_ticklines():
#         #     line.set_visible(False)
#         # for line in ax.zaxis.get_ticklines():
#         #     line.set_visible(False)
#     elif projection:
#         fig, axs = plt.subplots(1, 3, figsize=(15, 5))
#     for i in range(len(lines)):
#         coords = lines[i].split(' ')
#         # If it is the soma, plot it in any case
#         if "Soma" in coords[0]:
#             coords = lines[i-1].split(' ')
#             coords = [float(coord) for coord in coords]
#             # draw sphere
#             u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
#             x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
#             y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
#             z = np.cos(v)*float(coords[3]) + float(coords[2])
#             r = float(coords[3])
#             if plot_3d:
#                 ax.plot_surface(x, y, z, color="cornflowerblue", alpha=0.1)
#             elif projection:
#                 axs[0].plot(x, y, color="cornflowerblue")
#                 axs[1].plot(x, z, color="cornflowerblue")
#                 axs[2].plot(z, y, color="cornflowerblue")
#         elif len(coords) > 2 and len(coords) < 6:
#             coords = [float(coord) for coord in coords]
#             if coords[0] < min_[0]:
#                 min_[0] = coords[0]
#             if coords[1] < min_[1]:
#                 min_[1] = coords[1]
#             if coords[2] < min_[2]:
#                 min_[2] = coords[2]
#             if coords[0] > max_[0]:
#                 max_[0] = coords[0]
#             if coords[1] > max_[1]:
#                 max_[1] = coords[1]
#             if coords[2] > max_[2]:
#                 max_[2] = coords[2]

#             if i % 8 == 0:
#                 # draw sphere
#                 u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
#                 x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
#                 y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
#                 z = np.cos(v)*float(coords[3]) + float(coords[2])
#                 r = float(coords[3])

#                 if ((coords[0] < 0.54 and coords[0] > 0.48) and
#                     (coords[1] < 0.54 and coords[1] > 0.48) and
#                     (coords[2] < 0.54 and coords[2] > 0.48)):
#                     if plot_3d:
#                         ax.plot_surface(x, y, z, color="cornflowerblue", alpha=0.1)
#                     elif projection:
#                         axs[0].plot(x, y, color="cornflowerblue")
#                         axs[1].plot(x, z, color="cornflowerblue")
#                         axs[2].plot(z, y, color="cornflowerblue")
#             i += 1
            

#         elif len(coords) == 6:
#             coords = [float(coord) for coord in coords]
#             xmin = coords[0]
#             xmax = coords[1]
#             ymin = coords[2]
#             ymax = coords[3]
#             zmin = coords[4]
#             zmax = coords[5]

#             if plot_3d:
#                 ax.plot(xmin, ymin, zmin, color="cornflowerblue")
#                 ax.plot(xmax, ymax, zmax, color="cornflowerblue")
#             elif projection:
#                 axs[0].plot(xmin, ymin, color="cornflowerblue")
#                 axs[0].plot(xmax, ymax, color="cornflowerblue")
#                 axs[1].plot(xmin, zmin, color="cornflowerblue")
#                 axs[1].plot(xmax, zmax, color="cornflowerblue")
#                 axs[2].plot(zmin, ymin, color="cornflowerblue")
#                 axs[2].plot(zmax, ymax, color="cornflowerblue")
#     print(min_ - max_)
#     if not plot_3d:
#         distance_from_borders = 0.0007
#     if plot_traj:
#         with open(traj_file[2]) as f:
#             lines = f.readlines()
#             xp = []
#             yp = []
#             zp = []
#             for i in range(int(len(lines))):
#                 if i%4 == 0:
#                     xp.append(float(lines[i]))
#                 elif i%4 == 1:
#                     yp.append(float(lines[i]))
#                 elif i%4 == 2:
#                     zp.append(float(lines[i]))
#             xp = np.array(xp)
#             yp = np.array(yp)
#             zp = np.array(zp)
#             print(np.min(xp), np.max(xp), np.min(yp), np.max(yp), np.min(zp), np.max(zp))
#             if not plot_3d:
#                 axs[0].plot(xp, yp, 'b.', markersize=1)
#                 axs[1].plot(xp, zp, 'b.', markersize=1)
#                 axs[2].plot(zp, yp, 'b.', markersize=1)
#                 axs[0].set_xlim([0, 0.1])
#                 axs[0].set_ylim([0, 0.1])
#                 axs[0].set_xlabel("X [mm]")
#                 axs[0].set_ylabel("Y [mm]")
#                 axs[0].axvline(distance_from_borders, 0, 1)
#                 axs[0].axvline(0.1-distance_from_borders, ymin=0, ymax=1)
#                 axs[0].axhline(distance_from_borders, xmin=0, xmax=1)
#                 axs[0].axhline(y=0.1-distance_from_borders, xmin=0, xmax=1)
#                 axs[1].set_xlim([0, 0.1])
#                 axs[1].set_ylim([0, 0.1])
#                 axs[1].set_xlabel("X [mm]")
#                 axs[1].set_ylabel("Z [mm]")
#                 axs[1].axvline(distance_from_borders, 0, 1)
#                 axs[1].axvline(0.1-distance_from_borders, ymin=0, ymax=1)
#                 axs[1].axhline(distance_from_borders, xmin=0, xmax=1)
#                 axs[1].axhline(y=0.1-distance_from_borders, xmin=0, xmax=1)
#                 axs[2].set_xlim([0, 0.1])
#                 axs[2].set_ylim([0, 0.1])
#                 axs[2].set_xlabel("Z [mm]")
#                 axs[2].set_ylabel("Y [mm]")
#                 axs[2].axvline(distance_from_borders, 0, 1)
#                 axs[2].axvline(0.1-distance_from_borders, ymin=0, ymax=1)
#                 axs[2].axhline(distance_from_borders, xmin=0, xmax=1)
#                 axs[2].axhline(y=0.1-distance_from_borders, xmin=0, xmax=1)
#             else:
#                 # ax.scatter(xp, yp, zp, color="cornflowerblue", s=1)


#                 df = pd.DataFrame({"time": list(range(xp.shape[0])) ,"x" : xp, "y" : yp, "z" : zp})

#                 def update_graph(num):
#                     data=df[df['time']==num]
#                     graph.set_data(data.x, data.y)
#                     graph.set_3d_properties(data.z)
#                     title.set_text('3D Test, time={}'.format(num))
#                     return title, graph, 

#                 title = ax.set_title('3D Test')

#                 data=df[df['time']==0]
#                 graph, = ax.plot(data.x, data.y, data.z, linestyle="", marker="o")

#                 ani = FuncAnimation(fig, update_graph, 1000, 
#                                             interval=2, blit=True)

#                 plt.show()



#     # if plot_3d:
#     #     ax.set_xlim([min_lim, max_lim])
#     #     ax.set_ylim([min_lim, max_lim])
#     #     ax.set_zlim([min_lim, max_lim])
#     # elif projection:
#     #     for k in range(3):
#     #         axs[k].set_xlim([min_lim, max_lim])
#     #         axs[k].set_ylim([min_lim, max_lim])
#     # else:
#     #     for k in range(len(z_slice)):
#     #         axs[k].set_xlim([min_lim, max_lim])
#     #         axs[k].set_ylim([min_lim, max_lim])
#     plt.show()