# Monte Carlo Simulator modified for Overlapping Spheres

## Introduction

## Configuration Files - README

The simulation configuration is managed through `.conf` files located in the `instructions/conf/` directory. These files define the parameters for running simulations and include the following sections:

### General Parameters
The main simulation settings are defined as key-value pairs in the configuration file:

- **`N`**: Number of water molecules to simulate.
- **`T`**: Total number of steps in the simulation.
- **`duration`**: Duration of the simulation in seconds.
- **`diffusivity_intra`**: Intracellular diffusivity in \(m^2/s\).
- **`diffusivity_extra`**: Extracellular diffusivity in \(m^2/s\).
- **`scheme_file`**: Path to the scheme file used for simulation.
- **`exp_prefix`**: Path to the folder and desired prefix for output files, e.g., `path/to/folder/my_file`.
- **`scale_from_stu`**: Use standard units (`0` for no, `1` for yes).
- **`write_txt`**: Output DWI data as text files (`0` for no, `1` for yes).
- **`write_bin`**: Output DWI data as binary files (`.bfloat`) (`0` for no, `1` for yes).
- **`write_traj_file`**: Save water molecule trajectories (`0` for no, `1` for yes).
- **`num_process`**: Number of simulations to run simultaneously. It is recommended to set this to the number of CPU cores available.

### Cell Configuration
To include cells in the simulation, specify their properties in the following format:

```xml
<obstacle>
<axons_list>
path/to/swc/file
permeability global desired_permeability
</axons_list>
</obstacle>

## Voxel Size Adjustment

To configure the voxel size, specify the minimum and maximum values for each axis (`x`, `y`, `z`) in millimeters. The format is as follows:

```xml
<voxel>
xmin ymin zmin
xmax ymax zmax
</voxel>

## Sampling Area

The sampling area defines the region where the water molecules originate. Specify the minimum and maximum values for each axis as follows:

```xml
<sampling_area>
xmin ymin zmin
xmax ymax zmax
</sampling_area>

