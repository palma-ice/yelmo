#!/bin/bash

# Define user options

domain=Greenland
grid_name_src=GRL-64KM
grid_name_tgt=GRL-16KM

# Determine input grid file using options

nc_src=../ice_data/${domain}/${grid_name_src}/${grid_name_src}_REGIONS.nc 

if [ $grid_name_src = ERA5 ]
then
  nc_src=../ice_data/ERA5/era5_orography.nc 
fi

# Call cdo command to generate map weights between source grid and target grid:

cdo gencon,grid_${grid_name_tgt}.txt -setgrid,grid_${grid_name_src}.txt ${nc_src} scrip-con_${grid_name_src}_${grid_name_tgt}.nc

