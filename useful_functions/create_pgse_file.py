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

    shell_data = {}
    all_x_vals, all_y_vals, all_z_vals, all_b_values = [], [], [], []

    # Read the file and process lines
    with open(file_path, 'r') as file:
        for line in file:
            # Skip comment lines
            if line.startswith("#"):
                continue

            # Split the line into components
            parts = line.split()
            if len(parts) == 4:
                shell = int(parts[0])
                x_val, y_val, z_val = map(float, parts[1:])

                # Initialize shell if not already present
                if shell not in shell_data:
                    shell_data[shell] = {'x_vals': [], 'y_vals': [], 'z_vals': []}
                shell_data[shell]['x_vals'].append(x_val)
                shell_data[shell]['y_vals'].append(y_val)
                shell_data[shell]['z_vals'].append(z_val)

    if (len(b_values) != len(shell_data)):
        raise ValueError("Number of b-values does not match number of shells in text file : ", file_path)

    # Consolidate all shells into separate x, y, and z lists
    for shell, values in shell_data.items():
        all_x_vals.extend(values['x_vals'])
        all_y_vals.extend(values['y_vals'])
        all_z_vals.extend(values['z_vals'])
        all_b_values.extend([b_values[shell-1]] * len(values['x_vals']))

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
    x_vals, y_vals, z_vals, b_values, delta_values, te_values, small_delta_values, output_file, intervals = None
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
        if intervals is not None:
            file.write("VERSION: PGSE_INTERVALS\n")
        else:
            file.write("VERSION: STEJSKALTANNER\n")
        for (x, y, z, b_value), (Delta, TE, delta) in itertools.product(
            zip(x_vals, y_vals, z_vals, b_values), zip(delta_values, te_values, small_delta_values)
        ):
            # Compute the gradient strength G
            G = (math.sqrt(b_value / (Delta - (delta / 3)))) / (delta * GYROMAGNETIC_RATIO)
            G = round(G, 8)
            Delta, TE, delta = map(round, [Delta, TE, delta], [8, 8, 8])

            # Write to file
            if intervals is not None:
                file.write(f"{x} {y} {z} {G} {Delta} {delta} {TE} {intervals}\n")
            else:
                file.write(f"{x} {y} {z} {G} {Delta} {delta} {TE}\n")



def write_scheme_file(
    x_vals, y_vals, z_vals, b_values, delta_values, te_values, small_delta_values, output_file, intervals = None
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
        if intervals is not None:
            file.write("VERSION: PGSE_INTERVALS\n")
        else:
            file.write("VERSION: STEJSKALTANNER\n")
        for (x, y, z, b_value, Delta, TE, delta) in zip(x_vals, y_vals, z_vals, b_values, delta_values, te_values, small_delta_values):
            # Compute the gradient strength G
            G = (math.sqrt(b_value / (Delta - (delta / 3)))) / (delta * GYROMAGNETIC_RATIO)
            G = round(G, 8)
            Delta, TE, delta = map(round, [Delta, TE, delta], [8, 8, 8])

            # Write to file
            if intervals is not None:
                file.write(f"{x} {y} {z} {G} {Delta} {delta} {TE} {intervals}\n")
            else:
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
    directions_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/directions/SMI_directions.txt"
    output_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/time_dependence_narrow_pulse_dki.scheme"

    # B values and directions
    b_values =[0, 500, 1000, 2000, 3000]
    x_vals, y_vals, z_vals, b_values = parse_direction_file(directions_path, b_values)
    print("x_vals: ", len(x_vals))

    delta_values = [0.15758333, 0.11033333, 0.04577778, 0.02633333, 0.01733333, 0.01244444]
    # TE, big_delta, small_delta
    small_delta_values = [0.004]*len(delta_values) # narrow pulse
    te_values = [delta_values[i]+small_delta_values[i]+0.005 for i in range(len(delta_values))]

    write_combined_directions_with_b_values(x_vals, y_vals, z_vals, b_values, delta_values, te_values, small_delta_values, output_file)


def fixed_direction() :
    output_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/quentins_exp_d_20.scheme"

    # B values and directions
    direction = [0.0, 0.0, 1.0]
    b_value =[0, 200, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
    delta = [0.08]
    small_delta = 0.0165

    nbr_combinations = len(b_value) * len(delta)
    x_vals = [direction[0]] * nbr_combinations
    y_vals = [direction[1]] * nbr_combinations
    z_vals = [direction[2]] * nbr_combinations
    b_values = list(itertools.chain.from_iterable([[b] * len(delta) for b in b_value]))
    print(b_values)
    delta_values = [delta]*len(b_values)
    #flatten delta_values
    delta_values = list(itertools.chain.from_iterable(delta_values))
    print(delta_values)
    small_delta_values = [small_delta] * len(delta_values)
    te_values = [delta_values[i]+small_delta_values[i]+0.005 for i in range(len(delta_values))]
    write_scheme_file(x_vals, y_vals, z_vals, b_values, delta_values, te_values, small_delta_values, output_file, intervals = 1000)

if __name__ == "__main__":

    fixed_direction()