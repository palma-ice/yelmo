#!/bin/bash

fldr = 'output/test'

make benchmarks
python run_yelmo.py -s -e benchmarks ${fldr}/expa par/gmd/yelmo_EISMINT_expa.nml
python run_yelmo.py -s -e benchmarks ${fldr}/expf par/gmd/yelmo_EISMINT_expf.nml
python run_yelmo.py -s -e benchmarks ${fldr}/expssa par/yelmo_EISMINT_ssa.nml
job run --shell -f -o ${fldr}/halfar -p eismint.dx=0.5,1.0,2.0,3.0,4.0,5.0,8.0,10.0 -- python run_yelmo.py -x -s -e benchmarks {} par/gmd/yelmo_HALFAR.nml

make initmip
python run_yelmo.py -r -e initmip ${fldr}/ant-pd par/gmd/yelmo_Antarctica_pd.nml
python run_yelmo.py -r -e initmip ${fldr}/ant-lgm par/gmd/yelmo_Antarctica_lgm.nml
python run_yelmo.py -s -e initmip ${fldr}/grl par/yelmo_Greenland_initmip.nml

make mismip
job run --shell -f -o ${fldr}/mismip/default -p ydyn.beta_gl_scale=0 ydyn.beta_gl_sep=0  ydyn.beta_gl_stag=0 mismip.dx=2.5,5.0,10.0,20.0 -- python run_yelmo.py -x -s -e mismip {} par/gmd/yelmo_MISMIP3D.nml
job run --shell -f -o ${fldr}/mismip/subgrid -p ydyn.beta_gl_scale=0 ydyn.beta_gl_sep=-1 ydyn.beta_gl_stag=3 mismip.dx=2.5,5.0,10.0,20.0 -- python run_yelmo.py -x -s -e mismip {} par/gmd/yelmo_MISMIP3D.nml
job run --shell -f -o ${fldr}/mismip/scaling -p ydyn.beta_gl_scale=2 ydyn.beta_gl_sep=-1 ydyn.beta_gl_stag=3 mismip.dx=2.5,5.0,10.0,20.0 -- python run_yelmo.py -x -s -e mismip {} par/gmd/yelmo_MISMIP3D.nml

make trough 
python run_yelmo.py -s -e trough ${fldr}/trough par/yelmo_TROUGH-F17.nml
