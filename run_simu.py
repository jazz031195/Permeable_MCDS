import json
import argparse
from pathlib import Path
import os 

def create_job(exp_path, number_of_rep, exec_type, conf_path):
    """
    Function that creates job.sh, the file to be launch on Urblauna

    Args:
        exp_path (pathlib.PosixPath) : path of the experience to be run (e.g. exp1/dendrites/n1)
        number_of_rep          (int) : number of times each simulation should be repeated
        exec_type              (str) : type of the executable (either "soma", "dendrites", "soma_dendrites", "release")
                                       soma : walkers start and stay in soma only
                                       dendrites : walkers start and stay in dendrites only
                                       soma_dendrites : walkers start in soma and dendrites, but don't exchange compartments
                                       release : walkers start in soma and dendrites, and can exchange compartments
        conf_path (pathlib.PosixPath): path of the conf file (e.g. exp1/dendrites/n1/neurons.conf)
    """

    with open (exp_path / 'job.sh', 'w') as file:  
        file.write('#!/bin/bash -l\n')
        file.write('#SBATCH --job-name=sim_job\n')
        file.write('#SBATCH --nodes=1\n')
        file.write('#SBATCH --mem-per-cpu=16G\n')
        file.write('#SBATCH --time=71:59:00\n')
        file.write('#SBATCH --ntasks=1\n')
        file.write('#SBATCH --cpus-per-task=48\n')
        file.write('#SBATCH -o OUTPUTS/out/%j.%a.%N.%x.out\n')
        file.write('#SBATCH -e OUTPUTS/err/%j.%a.%N.%x.err\n')
        file.write('#SBATCH --array=1-5\n')

        file.write('\n')

        file.write(f'for ((i=1;i<={number_of_rep};i++));\n')
        file.write('do\n')
        file.write('\t echo "Starting creation of substrate ..."\n')
        file.write('\t echo " N : $n"\n')
        file.write('\t echo " T : $t"\n')
        file.write(f'\t chmod u+x ./build/MC-DC_Simulator_{exec_type}\n')
        file.write(f'\t ./build/MC-DC_Simulator_{exec_type} {conf_path}\n')
        file.write(f'done \n')

    

def create_conf(exp_path, N, T):
    # Read useful parameters from file
    f    = open(exp_path / "params.json")
    data = json.load(f)
    duration          = data.get("duration")          # Simulation duration in s
    diff_intra        = data.get("diffusivity_intra") # Simulation diffusivity in m²/s (around 3e-9 m²/s = 3 um²/ms for water, at room temperature)
    diff_extra        = data.get("diffusivity_extra") # Simulation diffusivity in m²/s (around 3e-9 m²/s = 3 um²/ms for water, at room temperature)
    scale             = data.get("scale")             # Boolean, True if the scheme file is in standard unit m, s
    write_txt         = data.get("write_txt")         # Boolean, write outputs as txt files
    write_bin         = data.get("write_bin")         # Boolean, write outputs as binary files
    write_traj_file   = data.get("write_traj_file")   # Boolean, store walkers' trajectories
    write_hit_file    = data.get("write_hit_file")    # Boolean, store walkers' hit TODO [ines] : what is that ?
    write_full_c_file = data.get("write_full_c_file") # Boolean, store full c file TODO [ines] : what is that ?
    # ini_pos           = data.get("ini_pos")         # Str, either "intra" or "extra", for intracellular and extracellular, respectively
    ini_pos_file      = data.get("ini_pos_file")      # Str, either "intra" or "extra", for intracellular and extracellular, respectively
    Number_processes  = data.get("Number_processes")  # Int, number of parallel processes
    # Number_processes  = 24  # Int, number of parallel processes
    Path_to_scheme    = data.get("Path_to_scheme")    # Str, relative path to scheme file
    Path_to_neurons   = data.get("Path_to_neurons")   # Str, relative path to neurons_list.swc or neurons_list.txt
    Path_to_conf      = data.get("Path_to_conf")      # Str, relative path to neurons.conf
    permeability      = data.get("permeability")      # Float, permability value
    sphere_overlap    = data.get("sphere_overlap")    # the spheres are radius/sphere_overlap appart
    voxel_lb          = data.get("voxel_lb")          # [x, y, z], list of floats, lower bound of the voxel
    voxel_ub          = data.get("voxel_ub")          # [x, y, z], list of floats, upper bound of the voxel
    funnel            = data.get("funnel")            # Boolean, do a funnel between soma & dendrites
    seed              = data.get("seed")              # Int, seed number

    # TODO [ines] : add ini_walkers_file

    with open (exp_path / Path_to_conf, 'w') as file:  
        file.write(f'N {N}\n')
        file.write(f'T {T}\n')
        file.write(f'duration {duration}\n')
        file.write(f'diffusivity_intra {diff_intra}\n')
        file.write(f'diffusivity_extra {diff_extra}\n')
        
        file.write('\n')

        file.write(f'exp_prefix {exp_path}/\n')

        file.write('\n')

        file.write(f'scheme_file {exp_path}/{Path_to_scheme}\n')
        file.write(f'scale_from_stu {scale}\n')

        file.write('\n')

        file.write(f'write_txt {write_txt}\n')
        file.write(f'write_bin {write_bin}\n')
        file.write(f'write_traj_file {write_traj_file}\n')
        file.write(f'write_hit_file {write_hit_file}\n')
        file.write(f'write_full_c_file {write_full_c_file}\n')

        file.write('\n')

        if ("swc" in Path_to_neurons) or ("txt" in Path_to_neurons):
            file.write('<obstacle>\n')
            file.write('<neurons_list>\n')
            file.write(f'{exp_path.parent}/{Path_to_neurons}\n')
            file.write(f'permeability global {permeability}\n')
            file.write(f'sphere_overlap {sphere_overlap}\n')
            file.write(f'funnel {funnel}\n')
            file.write('</neurons_list>\n')
            file.write('</obstacle>\n')
        else:
            file.write('<obstacle>\n')
            file.write('<ply>\n')
            file.write(f'{exp_path.parent}/{Path_to_neurons}\n')
            file.write(f'permeability {permeability}\n')
            file.write('ply_scale 0.001\n')
            file.write('</ply>\n')
            file.write('</obstacle>\n')

        file.write('\n')

        file.write('<voxels>\n')
        file.write(f'{voxel_lb[0]} {voxel_lb[1]} {voxel_lb[2]}\n')
        file.write(f'{voxel_ub[0]} {voxel_ub[1]} {voxel_ub[2]}\n')
        file.write('</voxels>\n')

        file.write('\n')

        # file.write(f'ini_walkers_pos {ini_pos}\n')
        file.write(f'ini_walkers_file {exp_path.parent}/{ini_pos_file}.txt\n')

        file.write('\n')

        # file.write(f'seed {seed}\n')
        file.write(f'num_process {Number_processes}\n')

        file.write('\n')

        file.write('<END>')

    return 0

parser = argparse.ArgumentParser()
parser.add_argument('-p', required=True, help='path of the parameter file for the simulation')

args = parser.parse_args()

f    = open(args.p)
data = json.load(f)
simu_to_launch = data.get("simu_to_launch") # List of str, simulations to launch, e.g. "exp1/dendrites/overlap4/n1"
number_of_rep  = data.get("number_of_rep")  # int, number of repetitions per simulation
simu_type      = data.get("simu_type")      # str, type of simulation ("soma", "dendrites", "soma_dendrites", "release"):
                                            # soma : walkers start and stay in soma only
                                            # dendrites : walkers start and stay in dendrites only
                                            # soma_dendrites : walkers start in soma and dendrites, but don't exchange compartments
                                            # release : walkers start in soma and dendrites, and can exchange compartments
Ns             = data.get("N")
Ts             = data.get("T")

# For all simulations
for simu in simu_to_launch:
    for i in range(number_of_rep):
        for N in Ns:
            for T in Ts:
                os.system(f"mkdir {simu}/N_{N}_T_{T}")
                os.system(f"cp {simu}/params.json {simu}/N_{N}_T_{T}")
                os.system(f"cp {simu}/PGSE_21_dir_12_b.scheme {simu}/N_{N}_T_{T}")
                # Create the conf file
                create_conf(Path(f"{simu}/N_{N}_T_{T}"), N, T)
                # # Create the job.sh to be launched
                # create_job(Path(f"{simu}/N_{N}_T_{T}"), number_of_rep, simu_type, Path(f"{simu}/N_{N}_T_{T}") / "neurons.conf")

                # Launch the job.sh on the cluster
                # os.system(f"sbatch {simu}/N_{N}_T_{T}/job.sh")

                os.system(f"chmod u+x ./build/MC-DC_Simulator_{simu_type}")
                # print(f"{simu}/N_{N}_T_{T}/neurons.conf")
                os.system(f"./build/MC-DC_Simulator_{simu_type} {simu}/N_{N}_T_{T}/neurons.conf")

    # "N": [5000, 10000, 25000, 50000, 75000, 100000, 125000, 150000],
    # "T": [5000, 10000, 15000, 20000]