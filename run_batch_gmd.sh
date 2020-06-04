#!/bin/bash

fldr='tmp/yelmo1.1'


### BENCHMARK TESTS ###

make benchmarks

# EISMINT1 moving margin, EXPA and EXPF
python run_yelmo.py -s -e benchmarks ${fldr}/moving par-gmd/yelmo_EISMINT_moving.nml
python run_yelmo.py -s -e benchmarks ${fldr}/expa par-gmd/yelmo_EISMINT_expa.nml
python run_yelmo.py -s -e benchmarks ${fldr}/expf par-gmd/yelmo_EISMINT_expf.nml

# EISMINT1 EXPA with SSA velocity turned on for testing symmetry (not part of GMD suite of tests)
python run_yelmo.py -s -e benchmarks ${fldr}/expssa par/yelmo_EISMINT_ssa.nml

# Ensemble of HALFAR simulations with various values of 
# dx to test numerical convergence with analytical solution
job run --shell -f -o ${fldr}/halfar -p eismint.dx=0.5,1.0,2.0,3.0,4.0,5.0,8.0 -- python run_yelmo.py -x -s -e benchmarks {} par-gmd/yelmo_HALFAR.nml

# Ensemble of EISMINT1-moving simulations with various values of 
# dx and pc_eps to test adaptive timestepping
# Note: make sure to specify: eismint.time_end=25e3 yelmo.log_timestep=True ytherm.method='fixed'
job run --shell -f -o ${fldr}/moving_dts -p eismint.time_end=25e3 yelmo.log_timestep=True ytherm.method='fixed' eismint.dx=5.0,10.0,25.0,50.0,60.0 yelmo.pc_eps=1e-2,1e-1,1e0 -- python run_yelmo.py -x -s -e benchmarks {} par-gmd/yelmo_EISMINT_moving.nml

#yelmo.pc_tol=1.5

### INITMIP TESTS ### 

make initmip

# Antarctica present-day and LGM simulations
# In par-gmd/yelmo_Antarctica.nml, set control.clim_nm="clim_pd"
python run_yelmo.py -s -e initmip ${fldr}/ant-pd  par-gmd/yelmo_Antarctica.nml
# In par-gmd/yelmo_Antarctica.nml, set control.clim_nm="clim_lgm"
python run_yelmo.py -s -e initmip ${fldr}/ant-lgm par-gmd/yelmo_Antarctica.nml

# Or to run via batch call:
job run --shell -f -o ${fldr}/ant -a -p control.clim_nm="clim_pd","clim_lgm" -- python run_yelmo.py -x -s -e initmip {} par-gmd/yelmo_Antarctica.nml


# Greenland present-day simulation (not part of GMD suite of tests)
python run_yelmo.py -s -e initmip ${fldr}/grl par/yelmo_Greenland_initmip.nml

job run --shell -f -o tmp/perform2 -a -p control.clim_nm="clim_pd_grl" yelmo.domain="Greenland"  ydyn.solver="diva","hybrid" yelmo.grid_name="GRL-32KM","GRL-16KM","GRL-8KM" -- python run_yelmo.py -x -s -q priority -w 24 -e initmip {} par/yelmo_initmip.nml 
job run --shell -f -o tmp/perform2 -a -p control.clim_nm="clim_pd_ant" yelmo.domain="Antarctica" ydyn.solver="diva","hybrid" yelmo.grid_name="ANT-32KM","ANT-16KM" -- python run_yelmo.py -x -s -q priority -w 24 -e initmip {} par/yelmo_initmip.nml 


### MISMIP TESTS ###

make mismip

# For faster, less high-resolution simulations:
job run --shell -f -o ${fldr}/mismip/default -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=0 mismip.dx=5.0,10.0,20.0 -- python run_yelmo.py -x -s -q short -w 24 -e mismip {} par-gmd/yelmo_MISMIP3D.nml
job run --shell -f -o ${fldr}/mismip/subgrid -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=3 mismip.dx=5.0,10.0,20.0 -- python run_yelmo.py -x -s -q short -w 24 -e mismip {} par-gmd/yelmo_MISMIP3D.nml
job run --shell -f -o ${fldr}/mismip/scaling -p ydyn.beta_gl_scale=2 ydyn.beta_gl_stag=3 mismip.dx=5.0,10.0,20.0 -- python run_yelmo.py -x -s -q short -w 24 -e mismip {} par-gmd/yelmo_MISMIP3D.nml

# For high-resolution sims only (medium:
job run --shell -f -o ${fldr}/mismip-hi/default -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=0 mismip.dx=2.5 -- python run_yelmo.py -x -s -q medium -w 60 -e mismip {} par-gmd/yelmo_MISMIP3D.nml
job run --shell -f -o ${fldr}/mismip-hi/subgrid -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=3 mismip.dx=2.5 -- python run_yelmo.py -x -s -q medium -w 60 -e mismip {} par-gmd/yelmo_MISMIP3D.nml
job run --shell -f -o ${fldr}/mismip-hi/scaling -p ydyn.beta_gl_scale=2 ydyn.beta_gl_stag=3 mismip.dx=2.5 -- python run_yelmo.py -x -s -q medium -w 60 -e mismip {} par-gmd/yelmo_MISMIP3D.nml

### ISMIP-HOM TESTS ###

make ismiphom 

job run --shell -f -o ${fldr}/ismiphom/expa/hybrid -p control.experiment="EXPA" control.L=5,10,20,40,80,160 ydyn.solver="hybrid" -- python run_yelmo.py -x -r -e ismiphom {} par/yelmo_ISMIPHOM.nml
job run --shell -f -o ${fldr}/ismiphom/expa/diva   -p control.experiment="EXPA" control.L=5,10,20,40,80,160 ydyn.solver="diva"   -- python run_yelmo.py -x -r -e ismiphom {} par/yelmo_ISMIPHOM.nml
job run --shell -f -o ${fldr}/ismiphom/expc/hybrid -p control.experiment="EXPC" control.L=5,10,20,40,80,160 ydyn.solver="hybrid" -- python run_yelmo.py -x -r -e ismiphom {} par/yelmo_ISMIPHOM.nml
job run --shell -f -o ${fldr}/ismiphom/expc/diva   -p control.experiment="EXPC" control.L=5,10,20,40,80,160 ydyn.solver="diva"   -- python run_yelmo.py -x -r -e ismiphom {} par/yelmo_ISMIPHOM.nml

# Trough simulation following Feldmann and Levermann (2017)  (not part of GMD suite of tests)
make trough

python run_yelmo.py -s -e trough ${fldr}/trough par/yelmo_TROUGH-F17.nml

### AGE TESTS ###

# In test_icetemp.f90, set experiment="eismint", and check that the following values are set:
    # t_start = 0.0       ! [yr]
    # t_end   = 300e3     ! [yr]
    # dt      = 0.5_prec  ! [yr]
    # dt_out  = 1000.0    ! [yr] 

    # T_pmp_beta = 9.7e-8         ! [K Pa^-1] EISMINT2 value (beta1 = 8.66e-4 [K m^-1])

    # call init_eismint_summit(ice1,smb=0.1_prec)
# Next, set prec=dp in yelmo_defs.f90 
# Compile and run 3 experiments (nz=12,32,52), corresponding to 10,30 and 50 internal ice points
# All output go to nc files in the folder "output/"
make clean
make icetemp 
./libyelmo/bin/test_icetemp.x 12 
./libyelmo/bin/test_icetemp.x 32 
./libyelmo/bin/test_icetemp.x 52 
# Also set prec=sp in yelmo_defs.f90, and run the nz=32 case again.
make clean 
make icetemp 
./libyelmo/bin/test_icetemp.x 32


