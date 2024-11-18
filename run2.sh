#!/bin/bash -l
./compile.sh
for i in {1..2}
do
    echo "Running simulation $i"
    ./MC-DC_Simulator "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/conf/model3.conf"
done
