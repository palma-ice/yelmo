## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using Oceananigans
using NCDatasets
using CairoMakie

include("YelmoIO.jl")
using .YelmoIO

# Define grid from NetCDF file
path = "/Users/alrobi001/models/yelmox/output/16KM/test/restart-0.000-kyr/yelmo_restart.nc"

# Load Yelmo Grids
grid2d, grid3d = load_grids_from_restart(path)

# Load Yelmo Fields
fields = load_fields_from_restart(path, grid2d, grid3d)

# At this point I have a Yelmo object stored in `fields`, and the associated 3D `grid2d` and 3D `grid3d` grids.

# Calculate velocity magnitude
uxy_s = @at (Center, Center, Center) sqrt(fields["ux_s"]^2 + fields["uy_s"]^2)