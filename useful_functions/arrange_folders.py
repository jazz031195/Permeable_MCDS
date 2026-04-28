import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob
import shutil
from useful_functions import read_and_extract_parameters, read_binary_file, array_to_nifti, get_files_from_folder, get_total_icvf, powder_avg 


def reorganise_files():
    path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/complex_axons/trajectories"
    signal_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/complex_axons/signal"

    os.makedirs(path, exist_ok=True)
    os.makedirs(signal_path, exist_ok=True)

    # Move *.bfloat from trajectories -> signal
    for name in sorted(os.listdir(path)):
        if not name.lower().endswith(".bfloat"):
            continue
        src = os.path.join(path, name)
        if not os.path.isfile(src):
            continue
        dst = os.path.join(signal_path, name)
        if os.path.exists(dst):
            print(f"{name} already exists in {signal_path}, skipping.")
            continue
        try:
            shutil.move(src, dst)
            print(f"Moved {name} -> {signal_path}")
        except Exception as e:
            print(f"Failed to move {name}: {e}")

    # Move *.txt from signal -> trajectories
    for name in sorted(os.listdir(path)):
        if not name.lower().endswith(".txt"):
            continue
        src = os.path.join(path, name)
        if not os.path.isfile(src):
            continue
        dst = os.path.join(signal_path, name)
        try:
            shutil.move(src, dst)
            print(f"Moved {name} -> {signal_path}")
        except Exception as e:
            print(f"Failed to move {name}: {e}")

def combine_trajectories():
    src_dir = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/complex_axons/trajectories"
    files = sorted(glob.glob(os.path.join(src_dir, "*.traj")))
    all_bhdr_files = sorted(glob.glob(os.path.join(src_dir, "*.bhdr")))

    compartments = ["intra", "extra"]
    packings = ["f_0.1", "f_0.2", "f_0.3", "f_0.4", "f_0.5", "f_0.6"]

    for packing in packings:
        for compartment in compartments:
            output_file = f"/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/complex_axons/combined_trajectories/{compartment}_{packing}.traj"
            combined_files = []
            bhdr_files = []
            for file in files:
                if compartment in file and packing in file:
                    combined_files.append(file)
            for file in all_bhdr_files:
                if compartment in file and packing in file:
                    bhdr_files.append(file)

            arrays = [np.fromfile(f, dtype="<f4") for f in combined_files]
            combined = np.concatenate(arrays, axis=0)
            combined.astype("<f4").tofile(output_file)
            print(f"Wrote {output_file} with {combined.size} float32 values")

            bhdr_files = [np.fromfile(f, dtype="<f4") for f in bhdr_files]
            # sum bhdr[1]
            total_walkers = sum(bhdr[1] for bhdr in bhdr_files)
            diffusion_time = sum(bhdr[0] for bhdr in bhdr_files) / len(bhdr_files) 
            nbr_steps = sum(bhdr[2] for bhdr in bhdr_files) / len(bhdr_files) 
            print(f"Total walkers for {output_file}: {total_walkers}")
            print(f"Average diffusion time for {output_file}: {diffusion_time:.4f} ms")
            # write new bhdr file
            bhdr_output_file = output_file.replace(".traj", ".bhdr")
            with open(bhdr_output_file, "wb") as f:
                f.write(np.array([diffusion_time, total_walkers, nbr_steps], dtype="<f4").tobytes())
            
def smi_output(packings, compartments, src_dir):

    if compartments is None :
        no_compartment = True
        compartments = [""]
    elif len(compartments) == 2:
        no_compartment = False
    else:
        print("compartments should be None or ['intra', 'extra']")
        return
    

    for packing in packings:
            
        for compartment in compartments:
            if not no_compartment:
                folder = src_dir + f"/simulations/{packing}_{compartment}/"
            else:
                folder = src_dir + f"/simulations/{packing}/"

            print(f"Processing folder: {folder}")
            traj_files = sorted(glob.glob(os.path.join(folder, "*.traj")))
            paralle_job_nbr = [traj_files.split("_")[-1].split(".")[0] for traj_files in traj_files]

            nbr_parallel_jobs = len(set(paralle_job_nbr))
            print(f"Number of parallel jobs detected: {nbr_parallel_jobs}")

            number_repetitions_required = 10

            def delete_replica_files(files): 
                rep_nbrs = [] 
                files_to_delete = [] 

                rep_in_all_files = True
                for i, file in enumerate(files): 
                    if "rep" in file :
                        continue
                    else:
                        rep_in_all_files = False
                        break
                
                if rep_in_all_files :
                    N = number_repetitions_required+1
                else:
                    N = number_repetitions_required
                for i, file in enumerate(files): 
                    if "rep" in file : 
                        rep_nbr = file.split("rep_")[-1].split("_")[0] 
                        rep_nbr = float(rep_nbr) 
                        if rep_nbr not in rep_nbrs: 
                            rep_nbrs.append(rep_nbr) 
                        if len(rep_nbrs) >= N: 
                            files_to_delete.append(i) 
                files_ = [] 
                for i, file in enumerate(files): 
                    if i not in files_to_delete: 
                        files_.append(file) 

                return files_

            bhdr_files = sorted(glob.glob(os.path.join(folder, "*.bhdr")))

            traj_files = delete_replica_files(traj_files)
            print(len(traj_files))
            bhdr_files = delete_replica_files(bhdr_files)
            print(len(bhdr_files))

            if int(len(traj_files)/int(nbr_parallel_jobs)) != number_repetitions_required:
                print(f"Warning: Expected {number_repetitions_required} repetitions files in {folder}, found {int(len(traj_files)/int(nbr_parallel_jobs))}.")
                assert(False)

            if (not no_compartment):
                output_traj_file = f"{src_dir}/combined/{packing}/trajectories_{compartment}.traj"
                output_bhdr_file = f"{src_dir}/combined/{packing}/trajectories_{compartment}.bhdr"
            else:
                output_traj_file = f"{src_dir}/combined/{packing}/trajectories.traj"
                output_bhdr_file = f"{src_dir}/combined/{packing}/trajectories.bhdr"
            
            #if folder doesnt exist create it 
            combined_folder = f"{src_dir}/combined/{packing}/"
                
            if not os.path.exists(combined_folder):
                os.makedirs(combined_folder)       

            arrays = [np.fromfile(f, dtype="<f4") for f in traj_files]
            traj = np.concatenate(arrays, axis=0)
            traj.astype("<f4").tofile(output_traj_file)
            print(f"Wrote {output_traj_file} with {traj.size} float32 values")

            bhdr_files_data = [np.fromfile(f, dtype="<f4") for f in bhdr_files]
            if (len(bhdr_files_data) == 0):
                print(f"No bhdr files found in {folder}, skipping.")
                continue
            for i, bhdr in enumerate(bhdr_files_data):
                if len(bhdr) < 3:
                    print(f"bhdr file in {folder} is too short, skipping.")
                    print(bhdr)
                    print("path : ", bhdr_files[i])
                    assert False
            # sum bhdr[1]
            total_walkers = sum(bhdr[1] for bhdr in bhdr_files_data)
            diffusion_time = sum(bhdr[0] for bhdr in bhdr_files_data) / len(bhdr_files_data) 
            nbr_steps = sum(bhdr[2] for bhdr in bhdr_files_data) / len(bhdr_files_data) 
            print(f"Total walkers for {output_bhdr_file}: {total_walkers}")
            print(f"Average diffusion time for {output_bhdr_file}: {diffusion_time:.4f} ms")
            # write new bhdr file
            with open(output_bhdr_file, "wb") as f:
                f.write(np.array([diffusion_time, total_walkers, nbr_steps], dtype="<f4").tobytes())

            DWI_files = sorted(glob.glob(os.path.join(folder, "*DWI.bfloat")))
            for DWI_file in DWI_files:
                if not no_compartment:
                    DWI_output_file = f"{src_dir}/combined/{packing}/DWI_{compartment}.bfloat"
                else:
                    DWI_output_file = f"{src_dir}/combined/{packing}/DWI.bfloat"

                arrays = [np.fromfile(f, dtype="<f4") for f in DWI_files]
                DWI = np.sum(arrays, axis=0)
                DWI.astype("<f4").tofile(DWI_output_file)
                print(f"Wrote {DWI_output_file} with {DWI.size} float32 values")
    if not no_compartment:
        create_DWI_total(packings, src_dir)


def create_DWI_total(packings, src_dir):

    for packing in packings:
        DWI_intra_path = f"{src_dir}/combined/{packing}/DWI_intra.bfloat"
        DWI_extra_path = f"{src_dir}/combined/{packing}/DWI_extra.bfloat"

        info_file = f"{src_dir}/{packing}_info.txt"

        icvf, ecvf = get_total_icvf(info_file)

        array_intra = np.fromfile(DWI_intra_path, dtype="<f4")
        array_extra = np.fromfile(DWI_extra_path, dtype="<f4")

        normalised_intra = (array_intra / np.max(array_intra)) * icvf
        normalised_extra = (array_extra / np.max(array_extra)) * ecvf
        DWI_total = np.sum([normalised_intra,normalised_extra], axis=0)
        min_value = np.min(DWI_total)

        if (min_value < 0):
            DWI_total = DWI_total - min_value

        DWI_output_path = f"{src_dir}/combined/{packing}/DWI.bfloat"
        DWI_total.astype("<f4").tofile(DWI_output_path)
        print(f"Wrote {DWI_output_path} with {DWI_total.size} float32 values")


def one_folder(folder_path):
    traj_files = sorted(glob.glob(os.path.join(folder_path, "*.traj")))
    print(traj_files)
    bhdr_files = sorted(glob.glob(os.path.join(folder_path, "*.bhdr")))
    output_traj_file = os.path.join(folder_path, "combined_trajectories.traj")
    output_bhdr_file = os.path.join(folder_path, "combined_trajectories.bhdr")
    arrays = [np.fromfile(f, dtype="<f4") for f in traj_files]
    traj = np.concatenate(arrays, axis=0)
    traj.astype("<f4").tofile(output_traj_file)
    print(f"Wrote {output_traj_file} with {traj.size} float32 values")
    bhdr_files = [np.fromfile(f, dtype="<f4") for f in bhdr_files]
    # sum bhdr[1]
    total_walkers = sum(bhdr[1] for bhdr in bhdr_files)
    diffusion_time = sum(bhdr[0] for bhdr in bhdr_files) / len(bhdr_files)
    nbr_steps = sum(bhdr[2] for bhdr in bhdr_files) / len(bhdr_files)
    print(f"Total walkers for {output_bhdr_file}: {total_walkers}")
    print(f"Average diffusion time for {output_bhdr_file}: {diffusion_time:.4f} ms")
    # write new bhdr file
    with open(output_bhdr_file, "wb") as f:
        f.write(np.array([diffusion_time, total_walkers, nbr_steps], dtype="<f4").tobytes())


            

