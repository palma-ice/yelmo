## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using CairoMakie

include("Yelmo.jl")

# Initialize Yelmo
ylmo = Yelmo("../par/yelmo_initmip.nml", "file", 0.0);

# Populate boundary fields
yelmo_set_var2D!("bnd_H_sed", fill(100.0,ylmo.g.nx,ylmo.g.ny) )

# Initialize Yelmo state
init_state!(ylmo,0.0, "robin-cold");

time_init, time_end, dt = 0.0, 5.0, 1.0;

for t in time_init:dt:time_end

    # Advance by dt
    time_step!(ylmo,dt);

    # Update boundary fields
    yelmo_set_var2D!("bnd_H_sed", fill(200.0,ylmo.g.nx,ylmo.g.ny) )

end

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