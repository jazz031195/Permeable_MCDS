import numpy as np
import pyvista as pv
import random
import pandas as pd
import matplotlib.colors as mcolors

def plot_bfloat_points(file_paths, binary, location):
    """
    Reads X, Y, Z coordinates from a .bfloat file and plots them using PyVista.

    Parameters:
    - file_path (str): Path to the .bfloat file.

    Returns:
    - None (Displays a 3D plot)
    """
    # Plot the points
    plotter = pv.Plotter()
    for file_path in file_paths:

        if (not location and binary):
            print(f"Reading .bfloat file: {file_path}")
            # Read the binary float32 data
            data = np.fromfile(file_path, dtype="float32")

            # Ensure the data can be reshaped into (N, 3) format
            if len(data) % 3 != 0:
                raise ValueError("Invalid .bfloat file: Data length is not a multiple of 3 (expected XYZ triplets).")

            # Reshape into N x 3 (X, Y, Z)
            points = data.reshape(-1, 3)

            print(f"Number of points read: {points.shape[0]}")

            points = points*1000

            # Create a PyVista point cloud
            cloud = pv.PolyData(points)
            plotter.add_mesh(cloud, color="blue", point_size=1, render_points_as_spheres=True, ambient = 0.5, diffuse = 0.5, specular = 0.5)
        elif location:
            # Read the file (replace 'data.txt' with your filename)
            data = pd.read_csv(file_path, delim_whitespace=True)

            data_intra = data.loc[data['location'] == "intra"]
            x = data_intra['x'].values
            y = data_intra['y'].values
            z = data_intra['z'].values

            points = np.column_stack((x, y, z))

            points = points*1000

            cloud = pv.PolyData(points)
            plotter.add_mesh(cloud, color="red", point_size=1, render_points_as_spheres=True, ambient = 0.5, diffuse = 0.5, specular = 0.5)

            data_extra = data.loc[data['location'] == "extra"]
            x = data_extra['x'].values
            y = data_extra['y'].values
            z = data_extra['z'].values

            points = np.column_stack((x, y, z))

            points = points*1000

            cloud = pv.PolyData(points)
            plotter.add_mesh(cloud, color="blue", point_size=1, render_points_as_spheres=True, ambient = 0.5, diffuse = 0.5, specular = 0.5)

    return plotter


def get_random_element(dictionary, seed):
    random.seed(seed)
    key, value = random.choice(list(dictionary.items()))
    return key, value


def load_spheres_to_dataframe(file_path):
    """
    Reads the spheres data from a text file and filters rows where x < 50 and y < 50.

    Args:
        file_path (str): Path to the spheres file.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only relevant spheres.
    """
    # Define column names based on the provided data structure
    col_names = ["id_ax", "id_sph", "id_branch", "Type", "X", "Y", "Z", "Rin", "Rout", "P"]
    
    # Load data into a DataFrame
    df = pd.read_csv(file_path, delim_whitespace=True, names=col_names, skiprows=1)

    df = df.loc[df["Type"] != "axon"]

    df = df.loc[df["id_ax"] == 2]

    return df

def plot_cells(file_path, plotter):
    """
    Reads a file, extracts X, Y, Z, and Rout values, and plots them as a point cloud.
    
    Parameters:
    - file_path (str): Path to the input file.

    Returns:
    - None (Displays a 3D scatter plot)
    """

    # Load data into a Pandas DataFrame
    df = load_spheres_to_dataframe(file_path)
    colors = mcolors.CSS4_COLORS
    df["color"] = df["id_ax"].apply(lambda x: get_random_element(colors, seed = x)[1])
    df["color_rgb"] = df["color"].apply(lambda c: mcolors.to_rgb(c)) 

    # Extract X, Y, Z coordinates and radius (Rout)
    x_values = df["X"].values
    y_values = df["Y"].values
    z_values = df["Z"].values
    R_values = df["Rout"].values  # Sphere radii
    cols = df["color_rgb"].values

    # Create a MultiBlock container (to store multiple sphere meshes)
    spheres = pv.MultiBlock()

    # Generate a sphere for each point
    for (x, y, z, r,color) in zip(x_values, y_values, z_values, R_values, cols):
        sphere = pv.Sphere(radius=r, center=(x, y, z),theta_resolution = 1 ,phi_resolution = 1)  # Create sphere at (x, y, z) with radius r
        spheres.append(sphere)

    
    # Add all spheres to the plot
    plotter.add_mesh(
        spheres,
        opacity=0.5,
        show_edges=False,
        color="green"
    )

    # Set camera position for better visualization
    plotter.view_isometric()

    # Show the plot
    plotter.show()


# Example usage:
if __name__ == "__main__":

    nbr_trajectories = 1
    binary = True
    location = False
    if binary:
        trajectory_paths = [f"/home/localadmin/Documents/MCDS/Permeable_MCDS/test_{i}.traj" for i in range(nbr_trajectories)]
    else:
        trajectory_paths = [f"/home/localadmin/Documents/MCDS/Permeable_MCDS/test_{i}.traj" for i in range(nbr_trajectories)]
    plotter = plot_bfloat_points(trajectory_paths, binary, location)
    plotter.show()

    #swc_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/axons_astrocytes/astrocytes_0.06.swc"
    #plot_cells(swc_file, plotter)

