#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

# icvf
c="1000000000";
# number of steps
t="11170";
straight="false";
declare -a locs=("intra" "extra"); 

declare -a perms=("0" "1e-5" "2.5e-5" "5e-5"); 

#path="/scratch/jnguyend/Simulation/MCDS-dFMRI/MCDC_Simulator_public-master"
path="/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public"
./compile.sh

for ((i=0;i<1;i++)); 
do
    for p in "${perms[@]}";
    do
        for l in "${locs[@]}";
        do
            echo "Starting creation of substrate ... "
            echo " Permeability : $p"
            echo " Location : $l"

            # run code to create axons list file
	        chmod u+x make_conf_file.sh
            input="axons/growth_icvf_0.50_cap_1_vox_50_factor_8_swell_0.01.swc"
            name="axons/perm_${p}_${l}"
            ./make_conf_file.sh -l "$l" -p 1 -P "$path" -S "$straight" -c "$c" -t "$t" -n "$name" -I "$input" -k "$p"
            ./MC-DC_Simulator "${path}/instructions/axons/configuration_run.conf"
        done
        
    done
done
  

