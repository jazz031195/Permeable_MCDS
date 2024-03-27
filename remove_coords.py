def read_coordinates_from_file(filename):
    coordinates = []
    with open(filename, 'r') as file:
        for line in file:
            if len(line.split()) > 0:
                x, y, z = map(float, line.split())
                coordinates.append((x, y, z))
    return coordinates

def write_coordinates_to_file(filename, coordinates):
    with open(filename, 'w') as file:
        for coord in coordinates:
            file.write(f"{coord[0]} {coord[1]} {coord[2]}\n")

def remove_coordinates(coordinates, coordinates_to_remove):
    return [coord for coord in coordinates if coord not in coordinates_to_remove]

# File paths
input_file = "results/ISMRM24/exch/overlap8/mesh_0005/n5/ini_pos_file.txt"
coordinates_to_remove_file = "results/ISMRM24/exch/overlap8/mesh_0005/n5/N_50000_T_15000/_0_extra.txt"
output_file = "results/ISMRM24/exch/overlap8/mesh_0005/n5/ini_pos_file_updated.txt"

# Read coordinates from files
coordinates = read_coordinates_from_file(input_file)
coordinates_to_remove = read_coordinates_from_file(coordinates_to_remove_file)

# Remove coordinates
updated_coordinates = remove_coordinates(coordinates, coordinates_to_remove)

# Write updated coordinates to output file
write_coordinates_to_file(output_file, updated_coordinates)