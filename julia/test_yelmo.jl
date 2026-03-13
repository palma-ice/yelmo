## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using CairoMakie

include("Yelmo.jl")

# Initialize Yelmo
yelmo_init("../par/yelmo_initmip.nml", "file", 0.0)

nx = 106
ny = 181
H = yelmo_get_H_ice(nx, ny)

heatmap(H)


# Coupled time loop
for t in 0.0:dt:t_end
    # Pull fields from Yelmo
    H = yelmo_get_H_ice(nx, ny)

    # Run your Julia-side models using H, etc.
    #bmb = run_basal_model(H, ...)

    # Push fields back into Yelmo's boundary conditions
    yelmo_set_bmb(bmb)

    # Advance Yelmo one step
    yelmo_step(t)
end