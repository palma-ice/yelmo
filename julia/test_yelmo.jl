## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using CairoMakie

include("Yelmo.jl")

# Initialize Yelmo
ylmo = Yelmo("../par/yelmo_initmip.nml", "file", 0.0);

# Populate boundary fields
ylmo.bnd.H_sed .= 100.0
yelmo_sync(ylmo)

# Initialize Yelmo state
init_state!(ylmo, 0.0, "robin-cold");

time_init, time_end, dt = 0.0, 5.0, 1.0;

for t in time_init:dt:time_end

    # Advance by dt
    time_step!(ylmo,t-ylmo.time);

    # Update boundary fields
    ylmo.bnd.H_sed .= 100.0 .+ t
    yelmo_sync(ylmo)

end

# Plot some data
heatmap(log10.(ylmo.dyn.uxy_s))

