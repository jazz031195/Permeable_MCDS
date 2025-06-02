import itertools
import math
import numpy as np

# Gyromagnetic ratio for hydrogen nuclei in rad/ms*T
GYROMAGNETIC_RATIO = 267.51525e3  # rad/ms*T


def parse_direction_file(file_path, b_values):

    """
    Parse a file containing shell and directional data into x, y, z values, and associate each shell with its corresponding b-value.

    Args:
        file_path (str): Path to the input file containing directional data. The file should have four columns: 
                         shell index, x-coordinate, y-coordinate, and z-coordinate.
        b_values (list): List of b-values, where each b-value corresponds to a shell index in the input file.

    Returns:
        tuple: 
            - list of float: x-coordinates for all directions.
            - list of float: y-coordinates for all directions.
            - list of float: z-coordinates for all directions.
            - list of float: b-values repeated for each direction in their respective shells.

    Raises:
        ValueError: If the number of b-values does not match the number of shells in the input file.
    """

    all_x_vals, all_y_vals, all_z_vals, all_b_values = [], [], [], []

    # Read the file and process lines
    with open(file_path, 'r') as file:
        for line in file:
            # Skip comment lines
            if line.startswith("#"):
                continue

            # Split the line into components
            parts = line.split()
            if len(parts) == 3:
                x_val, y_val, z_val = map(float, parts)

                all_x_vals.append([x_val] * len(b_values))
                all_y_vals.append([y_val] * len(b_values))
                all_z_vals.append([z_val] * len(b_values))
                all_b_values.append(b_values)
                
    return all_x_vals, all_y_vals, all_z_vals, all_b_values


def write_b_value_combinations(
    x_vals, y_vals, z_vals, b_values, delta_values, te_values, small_delta_values, output_file
):
    """
    Generate combinations of b-values and corresponding gradient strengths.

    Args:
        x_vals, y_vals, z_vals (list): Directional components.
        b_values (list): b-values.
        delta_values, te_values, small_delta_values (list): Timing parameters.
        output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as file:
        file.write("VERSION: STEJSKALTANNER\n")
        for (x, y, z), b_value, (Delta, TE, delta) in itertools.product(
            zip(x_vals, y_vals, z_vals), b_values, zip(delta_values, te_values, small_delta_values)
        ):
            # Compute the gradient strength G
            G = (math.sqrt(b_value / (Delta - (delta / 3)))) / (delta * GYROMAGNETIC_RATIO)
            G = round(G, 8)
            Delta, TE, delta = map(round, [Delta, TE, delta], [8, 8, 8])

            # Write to file
            file.write(f"{x} {y} {z} {G} {Delta} {delta} {TE}\n")


def write_combined_directions_with_b_values(
    x_vals, y_vals, z_vals, b_values, delta_values, te_values, small_delta_values, output_file
):
    """
    Generate combinations of directions and b-values.

    Args:
        x_vals, y_vals, z_vals (list): Directional components.
        b_values (list): b-values.
        delta_values, te_values, small_delta_values (list): Timing parameters.
        output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as file:
        file.write("VERSION: STEJSKALTANNER\n")
        for (Delta, TE, delta) in zip(delta_values, te_values, small_delta_values):
            for (xs, ys, zs, bs) in zip(x_vals, y_vals, z_vals, b_values):
                for (x, y, z, b_value) in zip(xs, ys, zs, bs):
                    # Compute the gradient strength G
                    G = (math.sqrt(b_value / (Delta - (delta / 3)))) / (delta * GYROMAGNETIC_RATIO)
                    G = round(G, 8)
                    Delta, TE, delta = map(round, [Delta, TE, delta], [8, 8, 8])
                    
                    b = GYROMAGNETIC_RATIO **2 * G ** 2 * delta **2 * (Delta - delta / 3)
                    print(f"{x} {y} {z} G {G} D {Delta} d {delta} TE {TE} b {b}")
                    # Write to file
                    file.write(f"{x} {y} {z} {G} {Delta} {delta} {TE}\n")


def write_combinations_with_fixed_gradients(
    x_vals, y_vals, z_vals, gradient_values, delta_values, te_values, small_delta_values, output_file
):
    """
    Generate combinations of directions and fixed gradient strengths.

    Args:
        x_vals, y_vals, z_vals (list): Directional components.
        gradient_values (list): Fixed gradient strengths.
        delta_values, te_values, small_delta_values (list): Timing parameters.
        output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as file:
        file.write("VERSION: STEJSKALTANNER\n")
        for (x, y, z), G, (Delta, TE, delta) in itertools.product(
            zip(x_vals, y_vals, z_vals), gradient_values, zip(delta_values, te_values, small_delta_values)
        ):
            G, Delta, TE, delta = map(round, [G, Delta, TE, delta], [4, 4, 4, 4])

            # Write to file
            file.write(f"{x} {y} {z} {G} {Delta} {delta} {TE}\n")


def time_dependence_narrow_pulse_dki() :
    directions_path = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/PGSE_21_dir.txt"
    output_file = "/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/PGSE_21_dir_12_b_9_td.scheme"

    # B values and directions
    b_values =[0, 200, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    x_vals, y_vals, z_vals, b_values = parse_direction_file(directions_path, b_values)
    print("x_vals: ", len(x_vals))

    delta_values = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    # TE, big_delta, small_delta
    small_delta_values = [0.0165]*len(delta_values) # narrow pulse
    te_values = [0.117]*len(delta_values)

    write_combined_directions_with_b_values(x_vals, y_vals, z_vals, b_values, delta_values, te_values, small_delta_values, output_file)

if __name__ == "__main__":

    time_dependence_narrow_pulse_dki()