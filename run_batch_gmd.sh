#!/bin/bash

fldr='tmp/53823a3c'


### BENCHMARK TESTS ###

make benchmarks

# EISMINT1 moving margin, EXPA and EXPF
./runylmo -s -q priority -w 1 -e benchmarks -o ${fldr}/moving -n par-gmd/yelmo_EISMINT_moving.nml
./runylmo -s -q priority -w 5 -e benchmarks -o ${fldr}/expa   -n par-gmd/yelmo_EISMINT_expa.nml
./runylmo -s -q priority -w 5 -e benchmarks -o ${fldr}/expf   -n par-gmd/yelmo_EISMINT_expf.nml

# EISMINT1 EXPA with SSA velocity turned on for testing symmetry (not part of GMD suite of tests)
./runylmo -s -e benchmarks -o ${fldr}/expssa -n par/yelmo_EISMINT_ssa.nml

# Ensemble of HALFAR simulations with various values of 
# dx to test numerical convergence with analytical solution
jobrun ./runylmo -s -e benchmarks -n par-gmd/yelmo_HALFAR.nml -- -o ${fldr}/halfar -p ctrl.dx=0.5,1.0,2.0,3.0,4.0,5.0,8.0

# Ensemble of EISMINT1-moving simulations with various values of 
# dx and pc_eps to test adaptive timestepping
# Note: make sure to specify: eismint.time_end=25e3 yelmo.log_timestep=True ytherm.method='fixed'
jobrun ./runylmo -s -e benchmarks -n par-gmd/yelmo_EISMINT_moving.nml -p eismint.time_end=25e3 yelmo.log_timestep=True ytherm.method='fixed' -- -o ${fldr}/moving_dts -p eismint.dx=5.0,10.0,25.0,50.0,60.0 yelmo.pc_eps=1e-2,1e-1,1e0


### INITMIP TESTS ### 

make initmip

# Antarctica present-day and LGM simulations (now with ydyn.solver='diva' by default)
./runylmo -s -q short -w 5 -e initmip -o ${fldr}/ant-pd  -n par-gmd/yelmo_Antarctica.nml -p ctrl.clim_nm="clim_pd"
./runylmo -s -q short -w 5 -e initmip -o ${fldr}/ant-lgm -n par-gmd/yelmo_Antarctica.nml -p ctrl.clim_nm="clim_lgm"

# Or to run via batch call:
jobrun ./runylmo -s -e initmip -n par-gmd/yelmo_Antarctica.nml -- -a -o ${fldr}/ant -p ctrl.clim_nm="clim_pd","clim_lgm"


# Greenland present-day simulation (not part of GMD suite of tests)
./runylmo -s -e initmip -o ${fldr}/grl -n par/yelmo_Greenland_initmip.nml


### MISMIP TESTS ###

make mismip

# For faster, less high-resolution simulations:
jobrun ./runylmo -s -q short -w 24 -e mismip -n par-gmd/yelmo_MISMIP3D.nml -- -o ${fldr}/mismip/default -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=0 ctrl.dx=2.5,5.0,10.0,20.0
jobrun ./runylmo -s -q short -w 24 -e mismip -n par-gmd/yelmo_MISMIP3D.nml -- -o ${fldr}/mismip/subgrid -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=3 ctrl.dx=2.5,5.0,10.0,20.0
jobrun ./runylmo -s -q short -w 24 -e mismip -n par-gmd/yelmo_MISMIP3D.nml -- -o ${fldr}/mismip/scaling -p ydyn.beta_gl_scale=2 ydyn.beta_gl_stag=3 ctrl.dx=2.5,5.0,10.0,20.0

# Trough simulation (not part of GMD suite of tests)
# MISMIP+ benchmark tests
# Feldmann and Levermann (2017)

make trough

# MISMIP+
./runylmo -s -e trough -o ${fldr}/mismip+ -n par/yelmo_MISMIP+.nml

# MISMIP+ ensemble (hybrid,diva), run on the cluster:
jobrun ./runylmo -s -q priority -w 5 -e trough -n par/yelmo_MISMIP+.nml -- -o ${fldr}/mismip+ -p ydyn.solver="hybrid","diva"


# F17
./runylmo -s -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough

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


