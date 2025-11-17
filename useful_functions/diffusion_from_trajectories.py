import numpy as np

def diffusion(trajectory_path, nbr_walkers, diffusion_time, nbr_dimesions, layout = 'space-major'):

    traj = np.fromfile(trajectory_path, dtype="float32")
    traj = traj[np.isfinite(traj)]

    diffusion_time = diffusion_time*1000 # Convert to ms

    print(f"Trajectory length: {traj.size}")

    length_traj = len(traj)
    if length_traj % (nbr_walkers * 3) != 0:
        raise ValueError("Trajectory length is not compatible with the number of walkers.")
    
    n_points = length_traj // 3
    T = n_points // nbr_walkers 

    print(f"Number of points per walker: {T}")

    arr = traj.reshape(nbr_walkers, T, 3)
    displacements = []
    for walker in range(nbr_walkers):
        first_step = arr[walker, 0, :]
        last_step = arr[walker, -1, :]
        displacement = np.linalg.norm(last_step - first_step) # mm
        displacement = displacement* 1000  # Convert to um
        displacements.append(displacement)
    root_mean_square_displacement = np.mean(np.square(displacements))

    diffusion_coefficient = root_mean_square_displacement  / (2 * nbr_dimesions * diffusion_time)
    return diffusion_coefficient
        


if __name__ == "__main__":

    trajectory_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/axons_astrocytes/test_0.traj"
    nbr_walkers = 1000
    diffusion_time = 0.077  # seconds
    nbr_dimesions = 3
    diffusion_coefficient = diffusion(trajectory_path, nbr_walkers, diffusion_time, nbr_dimesions)
    print(f"Diffusion coefficient: {diffusion_coefficient:.4f} um^2/s")