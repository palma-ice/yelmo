#!/bin/bash

# Define user options

domain=Laurentide
grid_name_src=LIS-32KM
grid_name_tgt=LIS-16KM

# Determine input grid file using options

nc_src=../ice_data/${domain}/${grid_name_src}/${grid_name_src}_REGIONS.nc 


# Call cdo command to generate map weights between source grid and target grid:

cdo gencon,grid_${grid_name_tgt}.txt -setgrid,grid_${grid_name_src}.txt ${nc_src} scrip-con_${grid_name_src}_${grid_name_tgt}.nc

