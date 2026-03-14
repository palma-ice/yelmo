## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using CairoMakie

include("Yelmo.jl")

# Initialize Yelmo
ylmo = Yelmo("../par/yelmo_initmip.nml", "file", 0.0);
ylmo2 = Yelmo("../par/yelmo_initmip.nml", "file", 0.0; alias="ylmo2");

# Populate boundary fields
ylmo.bnd.H_sed .= 100.0
yelmo_sync(ylmo)

ylmo2.bnd.H_sed .= 200.0
yelmo_sync(ylmo2)

# Initialize Yelmo state
init_state!(ylmo,0.0, "robin-cold");
init_state!(ylmo2,0.0, "robin-cold");

time_init, time_end, dt = 0.0, 5.0, 1.0;

for t in time_init:dt:time_end

    # Advance by dt
    time_step!(ylmo,dt);

    # Update boundary fields
    ylmo.bnd.H_sed .= 100.0 .+ t
    yelmo_sync(ylmo)

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