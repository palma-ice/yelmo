#!/bin/bash

# Define user options

domain=Greenland
grid_name_src=GRL-64KM

# Determine input grid file using options

nc_src=../ice_data/${domain}/${grid_name_src}/${grid_name_src}_REGIONS.nc 

if [ $grid_name_src = ERA5 ]
then
  nc_src=../ice_data/ERA5/era5_orography.nc 
fi

# Call cdo command to generate grid description file for the source grid
cdo griddes $nc_src > grid_${grid_name_src}.txt

text_to_print="Done. Note: for projected grids, after calling \`cdo griddes\`, it will be necessary to delete the first grid definition that appears in the file that is 'curvilinear'. We only need the second grid definition which is the projected grid itself. Also, you should change the units of the grid from 'kilometers' to 'km'."

echo ""
echo $text_to_print
echo ""

