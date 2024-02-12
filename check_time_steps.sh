#!/bin/bash -l
# chmod for access rights
# chmod u+x make_conf_file.sh
straight="false";
declare -a locs=("intra" "extra"); 
declare -a time_steps=("3000" "6000" "12000"); 
declare -a conc=("10000000"); 

path="/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public"
./compile.sh

for ((i=0;i<10;i++));
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
                name="instructions/sanity_check/loc_${l}_timesteps_${t}_conc_${c}"
                input="growth_icvf_0.70_cap_24_vox_50_straight_factor_2_0.swc"
                            
                ./make_conf_file_cylinders.sh -l "$l" -p 20 -P "$path" -S "$straight" -c "$c" -t "$t" -n "$name" -I "$input" -k "0"
                
                ./MC-DC_Simulator "${path}/instructions/configuration_run.conf"
                
            done
        done
    done
done