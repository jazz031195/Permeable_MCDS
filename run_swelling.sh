#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

straight="false";
declare -a locs=("intra" "extra"); 


#declare -a time_steps=("1000" "5000" "10000" "30000"); 
declare -a time_steps=("11170"); 

declare -a conc=("1000000000"); 
#declare -a conc=("10000000000"); 

#path="/scratch/jnguyend/Simulation/MCDS-dFMRI/MCDC_Simulator_public-master"
path="/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public"
./compile.sh
for ((i=0;i<5;i++)); 
    do
    declare -a swelling=("0" "0.0025" "0.005" "0.0075" "0.01");

    for s in "${swelling[@]}";
    do
        for l in "${locs[@]}";
            do
            for t in "${time_steps[@]}";
            do
                for c in "${conc[@]}";
                do
                    echo "Starting creation of substrate ... "
                    echo " Location : $l"
                    echo " Concentration : $c"
                    echo " Time steps : $t"

                    # run code to create axons list file
                    chmod u+x make_conf_file.sh
                    name="instructions/axons/icvf_70_vox_50/axons_${l}_swell_${s}"
                    input="instructions/axons/icvf_70_vox_50/growth_icvf_0.70_cap_20_vox_50_factor_8_0_swell_${s}.swc"
                            
                    ./make_conf_file.sh -l "$l" -p 20 -P "$path" -S "$straight" -c "$c" -t "$t" -n "$name" -I "$input" -k "0"
                        
                    ./MC-DC_Simulator "${path}/instructions/configuration_run.conf"
                        
             
                done
            done
        done
    done
done
