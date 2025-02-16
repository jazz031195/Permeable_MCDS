import numpy as np
import pyvista as pv


def plot_bfloat_points(file_paths):
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
        print(f"Reading .bfloat file: {file_path}")
        # Read the binary float32 data
        data = np.fromfile(file_path, dtype="float32")

        # Ensure the data can be reshaped into (N, 3) format
        if len(data) % 3 != 0:
            raise ValueError("Invalid .bfloat file: Data length is not a multiple of 3 (expected XYZ triplets).")

        # Reshape into N x 3 (X, Y, Z)
        points = data.reshape(-1, 3)

        # Create a PyVista point cloud
        cloud = pv.PolyData(points)
        plotter.add_mesh(cloud, color="blue", point_size=1, render_points_as_spheres=True, ambient = 0.5, diffuse = 0.5, specular = 0.5)
    plotter.show(title="3D Scatter Plot of .bfloat Points")

# Example usage:
if __name__ == "__main__":

    nbr_trajectories = 10
    trajectory_paths = [f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/verif_astrocytes/with_astrocytes_intra_{i}.traj" for i in range(nbr_trajectories)]
    plot_bfloat_points(trajectory_paths)

