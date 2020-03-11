#!/bin/bash

fldr=$1
file_out=${fldr}/ens_err.txt 

echo "Processing results in: ${fldr}"

#for D in ${fldr}/* ; do ./check_sim.x $D ; done > ${file_out}

# Initialize file with table header
echo 'sim   time     rmse_H      rmse_uxy     rmse_uxy_log' > ${file_out} 

# Loop over subfolders in fldr 
for D in ${fldr}/*
do
    if [[ -d "$D" ]]
    then
        #./check_sim.x $D >> ${file_out}
        out=`./check_sim.x $D` 
        echo $out 
        echo "$out" >> ${file_out}
    fi
done 

echo "Results written to: ${file_out}"
