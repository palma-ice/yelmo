#!/bin/bash

fldr=$1

for D in ${fldr}/* ; do ./check_sim.x $D ; done > ${fldr}/ens_err.txt