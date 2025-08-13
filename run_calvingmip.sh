#!/bin/bash

fldr='tmp/yelmo-calvingmip-2025-07-23'
ylmo='/Users/jablas001/Models/yelmo-ucm/yelmo_main_lsf'

### CALVINGMIP TESTS ###

make calving

# CalvingMIP steady state experiments.
#./runme -r -e calving -o ${fldr}/exp1 -n par/yelmo_calvingmip.nml -p ctl.exp="exp1" ycalv.calv_flt_method="exp1"
#./runme -r -e calving -o ${fldr}/exp3 -n par/yelmo_calvingmip.nml -p ctl.exp="exp3" ycalv.calv_flt_method="exp1"

# Running experiments. Be careful, we need first an equilibrated steady state.
restart_path=${ylmo}'/'${fldr}'/exp1/yelmo_restart.nc'

./runme -r -e calving -o ${fldr}/exp2 -n par/yelmo_calvingmip.nml -p ctl.exp="exp2" ctl.time_end=1000 ctl.dtt=1 ctl.dt2D_out=100 \
                            yelmo.restart=${restart_path} yelmo.restart_z_bed=True yelmo.restart_H_ice=True \
                            ycalv.calv_flt_method="exp2"

restart_path=${ylmo}'/'${fldr}'/exp3/yelmo_restart.nc'       

./runme -r -e calving -o ${fldr}/exp4 -n par/yelmo_calvingmip.nml -p ctl.exp="exp4" ctl.time_end=1000 ctl.dtt=1 ctl.dt2D_out=100 \
                            yelmo.restart=${restart_path} yelmo.restart_z_bed=True yelmo.restart_H_ice=True \
                            ycalv.calv_flt_method="exp2"

./runme -r -e calving -o ${fldr}/exp5 -n par/yelmo_calvingmip.nml -p ctl.exp="exp5" ctl.time_end=1000 ctl.dtt=1 ctl.dt2D_out=100 \
                            yelmo.restart=${restart_path} yelmo.restart_z_bed=True yelmo.restart_H_ice=True \
                            ycalv.calv_flt_method="exp5" ycalv.Hc_ref_flt=275.0
