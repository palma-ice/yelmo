## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using CairoMakie

include("Yelmo.jl")

# Initialize Yelmo
ylmo = yelmo_init("../par/yelmo_initmip.nml", "file", 0.0)

# Populate boundary fields
yelmo_set_var2D!("bnd_H_sed", fill(100.0,ylmo.g.nx,ylmo.g.ny) )

# [TO DO]

# Initialize Yelmo state
yelmo_init_state(0.0, "robin-cold")

# Advance to new time!
yelmo_step(5.0)

# Update data in julia
yelmo_update_variable_set!(ylmo.tpo, ylmo.v.tpo, "tpo")

# Update boundary fields
yelmo_set_var2D!("bnd_H_sed", fill(200.0,ylmo.g.nx,ylmo.g.ny) )

# Advance to new time!
yelmo_step(10.0)

# Get some data!
yelmo_get_var2D!(ylmo.tpo.H,    "tpo_H_ice")
yelmo_get_var2D!(ylmo.tpo.zs,   "tpo_z_srf")
yelmo_get_var3D!(ylmo.thrm.T,   "thrm_T_ice")
yelmo_get_var3D!(ylmo.dyn.uxy,  "dyn_uxy")
yelmo_get_var2D!(ylmo.dyn.uxy_s,"dyn_uxy_s")

# Plot some data
heatmap(log10.(ylmo.dyn.uxy_s))


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