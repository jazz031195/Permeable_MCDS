import itertools
import math
import numpy as np

giro = 267.51525e3 #rad/msXT

def get_directions(file_path):
    # Initialize dictionaries to store x_vals, y_vals, z_vals for each shell
    shells = {}
    all_x_vals = []
    all_y_vals = []
    all_z_vals = []

    # Read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Process each line (skip comment lines and header)
    for line in lines:
        # Skip comment lines (those that start with '#')
        if line.startswith("#"):
            continue

        # Split the line into components
        parts = line.split()
        
        # If the line has valid data (4 elements: shell, u_x, u_y, u_z)
        if len(parts) == 4:
            shell = int(parts[0])
            x_val = float(parts[1])
            y_val = float(parts[2])
            z_val = float(parts[3])

            # Add the values to the respective shell list
            if shell not in shells:
                shells[shell] = {'x_vals': [], 'y_vals': [], 'z_vals': []}
            shells[shell]['x_vals'].append(x_val)
            shells[shell]['y_vals'].append(y_val)
            shells[shell]['z_vals'].append(z_val)

    # Example of accessing and printing the lists for each shell
    for shell, values in shells.items():
        x_vals = values['x_vals']
        y_vals = values['y_vals']
        z_vals = values['z_vals']
        all_x_vals.append(x_vals)
        all_y_vals.append(y_vals)
        all_z_vals.append(z_vals)
    
    return all_x_vals, all_y_vals, all_z_vals
        






def fibonacci_sphere(samples=1, hemi=True):
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(samples * 2):  # Generate twice as many to ensure we get enough unique directions
        y = 1 - (i / float(samples * 2 - 1)) * 2  # y goes from 1 to -1
        if hemi and y < 0:
            continue
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points.append((x, y, z))

        if len(points) == samples:  # Stop once we have the desired number of points
            break

    return np.array(points)

def generate_combinations_b_vals(x_vals, y_vals, z_vals, b_vals, delta_vals, te_vals, delta_small_vals, output_file):
    with open(output_file, 'w') as f:
        f.write("VERSION: STEJSKALTANNER\n")
        for (x, y, z), b, (Delta, TE, delta) in itertools.product(
                zip(x_vals, y_vals, z_vals), b_vals, zip(delta_vals, te_vals, delta_small_vals)):
            G = (math.sqrt(b/(Delta- (delta/3))))/(delta*giro) #T/m
            # round G 
            G = round(G, 8)
            print("b : ", b)
            print("G : ",G)
            Delta = round(Delta, 8)
            TE = round(TE, 8)
            delta = round(delta, 8)

            f.write(f"{x} {y} {z} {G} {Delta} {delta} {TE}\n")

def generate_combinations_from_generated_directions(x_vals, y_vals, z_vals, b_vals, delta_vals, te_vals, delta_small_vals, output_file):
    with open(output_file, 'w') as f:
        f.write("VERSION: STEJSKALTANNER\n")
        for (x, y, z, b), (Delta, TE, delta) in itertools.product(
                zip(x_vals, y_vals, z_vals, b_vals), zip(delta_vals, te_vals, delta_small_vals)):
            G = (math.sqrt(b/(Delta- (delta/3))))/(delta*giro) #T/m
            # round G 
            G = round(G, 8)
            print("b : ", b)
            print("G : ",G)
            Delta = round(Delta, 8)
            TE = round(TE, 8)
            delta = round(delta, 8)

            f.write(f"{x} {y} {z} {G} {Delta} {delta} {TE}\n")

def generate_combinations_G(x_vals, y_vals, z_vals, g_vals, delta_vals, te_vals, delta_small_vals, output_file):
    with open(output_file, 'w') as f:
        f.write("VERSION: STEJSKALTANNER\n")
        for (x, y, z), G, (Delta, TE, delta) in itertools.product(
                zip(x_vals, y_vals, z_vals), g_vals, zip(delta_vals, te_vals, delta_small_vals)):

            # round G 
            G = round(G, 4)
            Delta = round(Delta, 4)
            TE = round(TE, 4)
            delta = round(delta, 4)

            f.write(f"{x} {y} {z} {G} {Delta} {delta} {TE}\n")

directions_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/directions/SMI_directions.txt"
# b0
x_vals = [0]
y_vals = [0]
z_vals = [1]
bs =[0]



all_x_vals, all_y_vals, all_z_vals = get_directions(directions_path)
b_values =[500, 1000, 1500, 2000, 2500] * len(all_x_vals)
for i in range(len(all_x_vals)):
    x_vals.extend(all_x_vals[i])
    y_vals.extend(all_y_vals[i])
    z_vals.extend(all_z_vals[i])
    bs.extend([b_values[i]]*len(all_x_vals[i]))


sqrt_diff_times = [0.08, 0.1, 0.15, 0.2, 0.25, 0.3]
diff_times = [1/(x*x) for x in sqrt_diff_times]
#diff_times = [50]
diff_times = [te/1000 for te in diff_times]
deltas = [0.004]*len(diff_times) # narrow pulse
#deltas = [0.0165]*len(diff_times)
Deltas = [diff_times[i]+deltas[i]/3 for i in range(len(diff_times))]
# tes = diff_times
tes = [Deltas[i]+deltas[i]+0.005 for i in range(len(Deltas))]

output_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/DTI_multi_shell.scheme"
#generate_combinations_b_vals(x_vals, y_vals, z_vals, bs, Deltas, tes, deltas, output_file)
generate_combinations_from_generated_directions(x_vals, y_vals, z_vals, bs, Deltas, tes, deltas, output_file)

