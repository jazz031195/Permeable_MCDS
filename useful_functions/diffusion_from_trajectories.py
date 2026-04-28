import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def diffusion(trajectory_path, nbr_walkers, diffusion_time, compartment):

    traj = np.fromfile(trajectory_path, dtype="float32")
    traj = traj[np.isfinite(traj)]

    diffusion_time = diffusion_time*1000 # Convert to ms

    print(f"Trajectory length: {traj.size}")

    length_traj = len(traj)

    T = 0

    if length_traj % (3) != 0:
        raise ValueError("Trajectory length is not compatible with the number of coordinates.")
    
    else:
        print(f"Trajectory length is compatible with the number of coordinates: {length_traj // 3} points.")

    if length_traj % (nbr_walkers * 3) != 0:
        raise ValueError("Trajectory length is not compatible with the number of walkers.")
    
    n_points = length_traj // 3
    T = int(n_points // nbr_walkers) 

    print(f"Number of points per walker: {T}")

    arr = traj.reshape(nbr_walkers, T, 3)
    displacements_axial = []
    displacements_radial = []
    for walker in range(nbr_walkers):
        first_step = arr[walker, 0, :]
        last_step = arr[walker, -1, :]

        if compartment == "intra":
            displacement_axial =  np.linalg.norm(last_step - first_step)
            displacement_radial = 0.0 
            nbr_dimesions_axial = 1
            nbr_dimesions_radial = 2
        else:

            axis = "z"
            if axis == "z":
                displacement_axial =  np.linalg.norm(last_step[2] - first_step[2])
                nbr_dimesions_axial = 1
                displacement_radial = np.linalg.norm(last_step[:2] - first_step[:2])
                nbr_dimesions_radial = 2
            elif axis == "x":
                displacement_axial =  np.linalg.norm(last_step[0] - first_step[0])
                nbr_dimesions_axial = 1
                displacement_radial = np.linalg.norm(last_step[1:-1] - first_step[1:-1])
                nbr_dimesions_radial = 2
            elif axis == "y":
                displacement_axial =  np.linalg.norm(last_step[1] - first_step[1])
                nbr_dimesions_axial = 1
                displacement_radial = np.linalg.norm(last_step[[0,2]] - first_step[[0,2]])
                nbr_dimesions_radial = 2
            
        displacement_axial = displacement_axial* 1000  # Convert to um
        displacement_radial = displacement_radial* 1000  # Convert to um
        displacements_axial.append(displacement_axial)
        displacements_radial.append(displacement_radial)

    root_mean_square_displacement_axial = np.mean(np.square(displacements_axial))
    diffusion_coefficient_axial = root_mean_square_displacement_axial  / (2 * nbr_dimesions_axial * diffusion_time)

    root_mean_square_displacement_radial = np.mean(np.square(displacements_radial))
    diffusion_coefficient_radial = root_mean_square_displacement_radial  / (2 * nbr_dimesions_radial * diffusion_time)
    
    return diffusion_coefficient_axial, diffusion_coefficient_radial

def create_diffusion_files(packings, compartments, src_dir):

    no_compartments = False
    if compartments is None:
        no_compartments = True
        compartments = [""]
    elif len(compartments) == 2:
        no_compartments = False

    for packing in packings:

        for compartment in compartments:
            if not no_compartments:
                trajectory_path = f"{src_dir}/combined/{packing}/trajectories_{compartment}.traj"
                bhdr_file = f"{src_dir}/combined/{packing}/trajectories_{compartment}.bhdr"
            else:
                trajectory_path = f"{src_dir}/combined/{packing}/trajectories.traj"
                bhdr_file = f"{src_dir}/combined/{packing}/trajectories.bhdr"
            
            bhdr = np.fromfile(bhdr_file, dtype="float32")
            bhdr = bhdr[np.isfinite(bhdr)]

            nbr_walkers = int(bhdr[1])
            print(f"Number of walkers: {nbr_walkers}")

            if packing != "diff":
                diffusion_time = 0.077  # seconds
            else:
                diffusion_time = 0.15  # seconds

            diffusion_coefficient_axial, diffusion_coefficient_radial = diffusion(trajectory_path, nbr_walkers, diffusion_time, compartment)

            # write in text file
            if no_compartments:
                name_output = f"/{src_dir}/combined/{packing}/diffusion_coefficients.txt"
            else:
                name_output = f"/{src_dir}/combined/{packing}/diffusion_coefficients_{compartment}.txt"
            with open(name_output, "w") as f:
                f.write(f"Diffusion coefficient radial: {diffusion_coefficient_radial:.4f} um^2/s\n")
                f.write(f"Diffusion coefficient axial: {diffusion_coefficient_axial:.4f} um^2/s\n")

def check_diffusion_outputs():
    compartments = ["intra", "extra"]
    packings = ["f_0.1", "f_0.2", "f_0.3", "f_0.4"]

    dzs = []
    dxys = []
    comps = []
    fs = []

    for packing in packings:
        for compartment in compartments:
            diff = f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/complex_axons/simulations/{packing}_{compartment}/combined/diffusion_coefficients.txt"
            d_xy = None
            d_z = None
            with open(diff, "r") as f:
                lines = f.readlines()
                for line in lines:
                    if "xy" in line:
                        d_xy = float(line.split(":")[1].strip().split(" ")[0])
                    elif "z" in line:
                        d_z = float(line.split(":")[1].strip().split(" ")[0])
            if d_xy is not None and d_z is not None:
                dxys.append(d_xy)
                dzs.append(d_z)
                comps.append(compartment)
                fs.append(packing)
    df = pd.DataFrame({"D_xy": dxys, "D_z": dzs, "Compartment": comps, "Packing": fs})
    sns.lineplot(data=df, x="Packing", y="D_xy", hue="Compartment", marker="o")
    plt.title("Diffusion Coefficient D_xy")
    plt.ylabel("D_xy (um²/ms)")
    plt.show()

    sns.lineplot(data=df, x="Packing", y="D_z", hue="Compartment", marker="o")
    plt.title("Diffusion Coefficient D_z")
    plt.ylabel("D_z (um²/ms)")
    plt.show()



