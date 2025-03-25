#!/bin/bash

fldr='tmp/yelmo-bench-v1.14'

runopt='-r'
#runopt='-rs'

### BENCHMARK TESTS ###

make benchmarks

# EISMINT1 moving margin, EXPA and EXPF
./runme ${runopt} -q 12h -w 1:00:00 -e benchmarks -o ${fldr}/moving -n par-gmd/yelmo_EISMINT_moving.nml
./runme ${runopt} -q 12h -w 5:00:00 -e benchmarks -o ${fldr}/expa   -n par-gmd/yelmo_EISMINT_expa.nml
./runme ${runopt} -q 12h -w 5:00:00 -e benchmarks -o ${fldr}/expf   -n par-gmd/yelmo_EISMINT_expf.nml

# EISMINT1-moving margin with DIVA solvers
# Note ydyn.solver="diva-noslip" is broken, as well as too high friction values, eg ydyn.beta_const=1e6
./runme ${runopt} -q 12h -w  1:00:00 -e benchmarks -o ${fldr}/moving-diva-noslip -n par-gmd/yelmo_EISMINT_moving.nml -p ydyn.solver="diva-noslip" ctrl.time_end=30e3 ctrl.dt2D_out=200
./runme ${runopt} -q 12h -w  1:00:00 -e benchmarks -o ${fldr}/moving-diva -n par-gmd/yelmo_EISMINT_moving.nml -p ydyn.solver="diva" ydyn.beta_method=0 ydyn.beta_const=1e4 ctrl.time_end=30e3 ctrl.dt2D_out=200

# EISMINT1 EXPA with SSA velocity turned on for testing symmetry (not part of GMD suite of tests)
./runme ${runopt} -e benchmarks -o ${fldr}/expssa -n par/yelmo_EISMINT_ssa.nml

# Ensemble of HALFAR simulations with various values of 
# dx to test numerical convergence with analytical solution
jobrun ./runme ${runopt} -e benchmarks -n par-gmd/yelmo_HALFAR.nml -o ${fldr}/halfar -p ctrl.dx=0.5,1.0,2.0,3.0,4.0,5.0,8.0

# Ensemble of EISMINT1-moving simulations with various values of 
# dx and pc_eps to test adaptive timestepping
# Note: make sure to specify: eismint.time_end=25e3 yelmo.log_timestep=True ytherm.method='fixed'
jobrun ./runme ${runopt} -e benchmarks -n par-gmd/yelmo_EISMINT_moving.nml -p eismint.time_end=25e3 yelmo.log_timestep=True ytherm.method='fixed' -o ${fldr}/moving_dts -p eismint.dx=5.0,10.0,25.0,50.0,60.0 yelmo.pc_eps=1e-2,1e-1,1e0




### INITMIP TESTS ### 

make initmip

# Antarctica present-day and LGM simulations (now with ydyn.solver='diva' by default)
./runme ${runopt} -q short -w 5:00:00 -e initmip -o ${fldr}/ant-pd  -n par-gmd/yelmo_Antarctica.nml -p ctrl.clim_nm="clim_pd"
./runme ${runopt} -q short -w 5:00:00 -e initmip -o ${fldr}/ant-lgm -n par-gmd/yelmo_Antarctica.nml -p ctrl.clim_nm="clim_lgm"

# Or to run via batch call:
jobrun ./runme -rs -e initmip -n par-gmd/yelmo_Antarctica.nml -a -o ${fldr}/ant -p ctrl.clim_nm="clim_pd","clim_lgm"


# Greenland present-day simulation (not part of GMD suite of tests)
#./runme -rs -q short -e initmip -o ${fldr}/grl -n par/yelmo_Greenland_initmip.nml
./runme -rs -q short -e initmip -o ${fldr}/grl -n par/yelmo_initmip.nml

### MISMIP TESTS ###

make mismip

# For faster, less high-resolution simulations:
jobrun ./runme -rs -q short -w 24:00:00 -e mismip -n par-gmd/yelmo_MISMIP3D.nml -o ${fldr}/mismip/default -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=0 ctrl.dx=2.5,5.0,10.0,20.0
jobrun ./runme -rs -q short -w 24:00:00 -e mismip -n par-gmd/yelmo_MISMIP3D.nml -o ${fldr}/mismip/subgrid -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=3 ctrl.dx=2.5,5.0,10.0,20.0
jobrun ./runme -rs -q short -w 24:00:00 -e mismip -n par-gmd/yelmo_MISMIP3D.nml -o ${fldr}/mismip/scaling -p ydyn.beta_gl_scale=2 ydyn.beta_gl_stag=3 ctrl.dx=2.5,5.0,10.0,20.0

# Trough simulation (not part of GMD suite of tests)
# MISMIP+ benchmark tests
# Feldmann and Levermann (2017)

make trough

# MISMIP+
./runme -rs -e trough -o ${fldr}/mismip+ -n par/yelmo_MISMIP+.nml

# MISMIP+ ensemble (hybrid,diva), run on the cluster:
jobrun ./runme -rs -q priority -w 5:00:00 -e trough -n par/yelmo_MISMIP+.nml -o ${fldr}/mismip+ -p ydyn.solver="hybrid","diva"


# F17
./runme -rs -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough

# ssa, dx=2km, [u0=100, cf_ref=5,10,20]:
./runme -rs -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-ssa -p ydyn.solver="ssa"
./runme -rs -q medium -w 48:00:00 -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-dx1 -p ctrl.dx=1.0
./runme -rs -q short  -w 24:00:00 -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-dx2 -p ctrl.dx=2.0
./runme -rs -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-u0.100-cf5.0  -p ydyn.beta_u0=100 ytill.cf_ref=5.0
./runme -rs -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-u0.100-cf10.0 -p ydyn.beta_u0=100 ytill.cf_ref=10.0
./runme -rs -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-u0.100-cf20.0 -p ydyn.beta_u0=100 ytill.cf_ref=20.0


### SLAB-S06 ###

# One run for testing
./runme -r -e trough -n par/yelmo_SLAB-S06.nml -o ${fldr}/slab06-test -p ctrl.dx=4 ydyn.ssa_iter_max=10


# Standard ensemble runs
jobrun ./runme -r -e trough -n par/yelmo_SLAB-S06.nml -o ${fldr}/slab06 -p ctrl.dx=0.5,1,2,4,8 

# Constant randomly-chosen fixed viscosity test
jobrun ./runme -r -e trough -n par/yelmo_SLAB-S06.nml -o ${fldr}/slab06-1 -p ctrl.dx=1,2,4,8 ydyn.visc_method=0 ydyn.visc_const=2e6

# Constant analytical velocity imposed
# (make sure to comment out call to ssa solver in velocity_ssa.f90)
jobrun ./runme -r -e trough -n par/yelmo_SLAB-S06.nml -o ${fldr}/slab06-2 -p ctrl.dx=0.5,1,2,4,8

# Standard ensemble, but different flow law exponent 
jobrun ./runme -r -e trough -n par/yelmo_SLAB-S06.nml -o ${fldr}/slab06-3 -p ctrl.dx=0.5,1,2,4,8 ydyn.beta_q=0.5

### SLAB Robinson 2022 ###

make slab 

# One simulation with solver of choice for each case
# No ice thickness evolution
jobrun ./runme ${runopt} -e slab -n par/yelmo_slab.nml -o ${fldr}/slab/weak   -p yelmo.pc_method="FE-SBE" ytopo.topo_fixed=True ctrl.dtt=0.1 ctrl.H0=1000 ctrl.H_stdev=0.0 ydyn.visc_const=1e5 ydyn.beta_const=1e3 ydyn.solver="ssa" ydyn.visc_method=3 ymat.rf_method=0
jobrun ./runme ${runopt} -e slab -n par/yelmo_slab.nml -o ${fldr}/slab/strong -p yelmo.pc_method="FE-SBE" ytopo.topo_fixed=True ctrl.dtt=0.1 ctrl.H0=500  ctrl.H_stdev=0.0 ydyn.visc_const=4e5 ydyn.beta_const=30  ydyn.solver="ssa" ydyn.visc_method=3 ymat.rf_method=0

# Full ensemble:
jobrun ./runme ${runopt} -e slab -n par/yelmo_slab.nml -o ${fldr}/slab-sd0.1/weak   -p yelmo.pc_method="FE-SBE" yelmo.pc_use_H_pred=True ctrl.dtt=0.0 ctrl.H0=1000 ctrl.H_stdev=0.1 ydyn.visc_const=1e5 ydyn.beta_const=1e3 ydyn.solver="diva","hybrid","l1l2","ssa","sia"
jobrun ./runme ${runopt} -e slab -n par/yelmo_slab.nml -o ${fldr}/slab-sd0.1/strong -p yelmo.pc_method="FE-SBE" yelmo.pc_use_H_pred=True ctrl.dtt=0.0 ctrl.H0=500  ctrl.H_stdev=0.1 ydyn.visc_const=4e5 ydyn.beta_const=30  ydyn.solver="diva","hybrid","l1l2","ssa","sia"

### CalvingMIP ###
#runopt='-rs  -q 12h -w 05:00:00'
#output/bench-2024-12-01

./runme ${runopt} -e calving -n par/yelmo_calving.nml -o ${fldr}/calvmip-exp01 -p ctl.exp="exp1" ctl.dt2D_out=200 ctl.time_end=10e3

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


