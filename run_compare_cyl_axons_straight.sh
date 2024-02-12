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
for ((i=0;i<4;i++)); 
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
                name="instructions/axons_vs_cylinders/data/cylinders_icvf_0.70_cap_24_vox_50_straight_${l}"
                input="growth_icvf_0.70_cap_24_vox_50_straight_factor_2_0.swc"
                           
                ./make_conf_file_cylinders.sh -l "$l" -p 20 -P "$path" -S "$straight" -c "$c" -t "$t" -n "$name" -I "$input" -k "0"
                
                ./MC-DC_Simulator "${path}/instructions/configuration_run.conf"
                
            done
        done
        
    done


    declare -a factors=("2" "4" "8" "16" "32"); 
    #declare -a factors=("4" "8" "16" "32"); 
    for f in "${factors[@]}";
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
                    name="instructions/axons_vs_cylinders/data/axons_icvf_0.70_cap_24_vox_50_straight_factor_${f}_${l}"
                    input="growth_icvf_0.70_cap_24_vox_50_straight_factor_${f}_0.swc"
                            
                    ./make_conf_file.sh -l "$l" -p 20 -P "$path" -S "$straight" -c "$c" -t "$t" -n "$name" -I "$input" -k "0"
                        
                    ./MC-DC_Simulator "${path}/instructions/configuration_run.conf"
                        
             
                done
                
            done
        done
    done
done
