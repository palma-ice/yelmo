#!/bin/bash

date_and_time=$(date '+%Y-%m-%d-%H-%M') 

fldr="output/unit_tests/${date_and_time}"
fldr_ref='output/unit_tests/v1.12.5'

echo "Unit tests: ${fldr}"
echo "Comparing with: ${fldr_ref}"

# Compile the benchmarks code
#make benchmarks

# EISMINT1 moving margin, EXPA and EXPF
#./runylmo -r -e benchmarks -o ${fldr}/moving -n par-gmd/yelmo_EISMINT_moving.nml -p ctrl.dt2D_out=1000 ctrl.time_end=20000
#./runylmo -r -e benchmarks -o ${fldr}/expa   -n par-gmd/yelmo_EISMINT_expa.nml   -p ctrl.dt2D_out=1000 ctrl.time_end=20000
#./runylmo -r -e benchmarks -o ${fldr}/expf   -n par-gmd/yelmo_EISMINT_expf.nml   -p ctrl.dt2D_out=1000 ctrl.time_end=20000

