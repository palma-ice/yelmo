## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using CairoMakie

include("Yelmo.jl")

# Northern hemisphere #
begin
    ylmo_nh = Yelmo("../par/yelmo_initmip.nml", "file", 0.0);

    # Populate boundary fields
    ylmo_nh.bnd.H_sed .= 100.0
    yelmo_sync(ylmo_nh);

    init_state!(ylmo_nh, 0.0, "robin-cold");
end;

# Southern hemisphere #
begin
    ylmo_sh = Yelmo("../par/yelmo_initmip.nml", "file", 0.0; alias="ylmo2");

    ylmo_sh.bnd.H_sed .= 200.0
    yelmo_sync(ylmo_sh);

    init_state!(ylmo_sh, 0.0, "robin-cold");
end;

## TIME LOOP ##

time_init, time_end, dt = 0.0, 5.0, 1.0;

for t in time_init:dt:time_end

    # Advance by dt
    time_step!(ylmo_nh,t-ylmo_nh.time);
    
    # Update boundary fields
    ylmo_nh.bnd.z_bed .+= 100.0
    yelmo_sync(ylmo_nh)

end

# Compare surface velocity fields
heatmap(ylmo_nh.dyn.uxy_s .- ylmo_sh.dyn.uxy_s)

