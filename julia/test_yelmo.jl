## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using CairoMakie

include("Yelmo.jl")

# Initialize Yelmo
yelmo_init("../par/yelmo_initmip.nml", "file", 0.0)

# Initialize Yelmo state
yelmo_init_state(0.0, "robin-cold")

# Advance to new time!
yelmo_step(5.0)

# Get some data!
grd = yelmo_get_grid_info()
H = yelmo_get_var2D(grd.nx, grd.ny, "tpo_H_ice")
zs = yelmo_get_var2D(grd.nx, grd.ny, "tpo_z_srf")
T = yelmo_get_var3D(grd.nx, grd.ny, grd.nz_aa, "thrm_T_ice")
uxy = yelmo_get_var3D(grd.nx, grd.ny, grd.nz_aa, "dyn_uxy")

uxy_s = yelmo_get_var2D(grd.nx, grd.ny, "dyn_uxy_s")

# Plot some data
heatmap(log10.(uxy_s))


# Coupled time loop
# for t in 0.0:dt:t_end
#     # Pull fields from Yelmo
#     H = yelmo_get_H_ice(nx, ny)

#     # Run your Julia-side models using H, etc.
#     #bmb = run_basal_model(H, ...)

#     # Push fields back into Yelmo's boundary conditions
#     yelmo_set_bmb(bmb)

#     # Advance Yelmo one step
#     yelmo_step(t)
# end