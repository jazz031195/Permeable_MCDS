# Monte Carlo Simulator modified for Overlapping Spheres

## Introduction

Monte Carlo simulations are a powerful computational technique used to model and analyze the diffusion of water molecules in biological tissues, as observed in diffusion-weighted imaging (DWI). These simulations rely on random sampling to approximate the physical and statistical properties of diffusion, making them particularly suitable for exploring complex, heterogeneous environments such as brain tissue. This tool was modified from https://github.com/jonhrafe/MCDC_Simulator_public to enable diffusioin of water molecules through overlapping spheres that model cells.

### What is DWI?

Diffusion-weighted imaging is a magnetic resonance imaging (MRI) technique that measures the diffusion of water molecules within tissues. This diffusion is influenced by the microstructural properties of the tissue. DWI provides insights into tissue structure, making it essential for studying diseases like stroke, multiple sclerosis, and brain tumors.

### Why Use Monte Carlo Simulations?

Monte Carlo simulations are employed in DWI to:

- Model Microstructures: Simulate the behavior of water molecules in complex tissue environments, including intracellular, extracellular, and restricted diffusion.
- Validate Models: Test analytical diffusion models against ground-truth simulations.
- Predict Signal: Generate synthetic DWI signals for various tissue configurations and experimental parameters.

### How It Works

- Initialization: Define the tissue environment, such as the dimensions of axons, cell membranes, and extracellular spaces. Set physical properties like diffusivity and permeability.
- Random Walks: Simulate the movement of individual water molecules over time using random walks. Each step accounts for diffusion, collisions with structures, and interactions with barriers.
- Signal Generation: Apply gradient pulses, as defined by a DWI scheme, and calculate the resulting signal attenuation.
- Analysis: Aggregate the simulated signals to predict DWI metrics, such as fractional anisotropy (FA), mean diffusivity (MD), and kurtosis.

### Benefits of Monte Carlo Simulations

- Accuracy: Capture the nuances of diffusion in non-ideal, heterogeneous environments.
- Flexibility: Model various tissue geometries and experimental conditions.
- Insights: Explore relationships between microstructure and measured DWI signals.

## Configuration Files 

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
- **`ini_walkers_pos`**: Initial compartment in which the molecules start. Can be : intra or extra. If this is not given, the water molecules can be inside or outside the cells.

Other parameters can be found in /src/parameters.h

### Cell Configuration
To include axons in the simulation, specify their properties in the following format:

```xml
<obstacle>
<axons_list>
path/to/swc/file
permeability global desired_permeability
</axons_list>
</obstacle>
```

To include glial cells or neurons, specify their properties in the following format:

```xml
<obstacle>
<glials_list>
path/to/swc/file
permeability global desired_permeability
</glials_list>
</obstacle>
```

If you want both glial cells and axons from the same substrate, write : 

```xml
<obstacle>
<axons_list>
path/to/swc/file
permeability global desired_permeability
</axons_list>
<glials_list>
path/to/swc/file
permeability global desired_permeability
</glials_list>
</obstacle>
```
The path/to/swc/file should be the same twice. 

Set desired_permeability to 0 for no permeability.

### Voxel Size Adjustment

To configure the voxel size, specify the minimum and maximum values for each axis (`x`, `y`, `z`) in millimeters. The format is as follows:

```xml
<voxel>
xmin ymin zmin
xmax ymax zmax
</voxel>
```

### Sampling Area

The sampling area defines the region where the water molecules originate. Specify the minimum and maximum values for each axis as follows:

```xml
<sampling_area>
xmin ymin zmin
xmax ymax zmax
</sampling_area>
```

## Scheme Files

Scheme files are located in the `instructions/scheme/` directory. These files are essential for defining the parameters of your PGSE (Pulsed Gradient Spin Echo) sequences. Each line in a scheme file describes a single gradient direction and associated parameters required for the simulation or analysis. The parameters are specified in the following order: X Y Z G Δ δ TE

### Explanation of PGSE Parameters
- **X, Y, Z**: Components of the gradient direction in 3D space (mm).
- **G**: Gradient amplitude (in millitesla per meter, mT/m).
- **Δ (big delta)**: Time interval between the diffusion gradients (in seconds).
- **δ (small delta)**: Duration of the diffusion gradients (in seconds).
- **TE**: Echo time (in seconds).

If you want multiple DWI values during the PGSE sequence, add an additional parameter:
- **num_intervals**: Number of DWI values you wish to obtain during sequence.

### Example PGSE Scheme File

```xml
VERSION: STEJSKALTANNER 
0 0 1 0.0 0.05 0.0165 0.07
0 0 1 0.01518807 0.05 0.0165 0.07
0 0 1 0.02401444 0.05 0.0165 0.07
```
### Example if you want DWI outputs during PGSE sequence 
```xml
VERSION: PGSE_INTERVALS
0 0 1 0.0 0.05 0.0165 0.07 1000
0 0 1 0.01518807 0.05 0.0165 0.07 1000
0 0 1 0.02401444 0.05 0.0165 0.07 1000
```

### Example for waveforms
```xml
VERSION: WAVEFORM
0.06696
6697
11
0.0 0.0 0.0
0.0 0.0 0.0003316466017446737
0.0 0.0 0.0006632932034893485
0.0 0.0 0.00099493980523402
0.0 0.0 0.0013265864069786936
```

### Generating Isotropic Directions

If isotropic directions are required for your simulation, you can generate them using external tools. One recommended tool is Emmanuel Caruyer’s q-space sampling tool (http://www.emmanuelcaruyer.com/q-space-sampling.php). This tool allows you to produce text files containing isotropic gradient directions.

### Automating Scheme File Creation

Once you have the directional data, you can use the create_pgse_file Python script located in the /useful_functions/ directory to generate a .scheme file with the desired parameters. This script automates the process, ensuring consistency and reducing manual errors.

## SWC files

To run a Monte Carlo simulation on a specific substrate, first create it using the tool of your choice (ex: CATERPillar). The path to the output that has the .swc must be put in the configuration file. An example of such a substrate can be found in /example_swc_files/. Information on this substrate can be found in the .txt file. 

## Running the MCDS 

Once the configuration and scheme files are created, simply execute :

```xml
./run_mcds.sh
```