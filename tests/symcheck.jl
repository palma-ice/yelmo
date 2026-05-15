#!/usr/bin/env julia
# Quantify L<->R, T<->B and 180-deg-rotation asymmetry of H_ice (lithk)
# in one or more yelmo2D.nc outputs.
#
# Intended use: cross-check that a symmetric setup (e.g. CalvingMIP exp1 on
# its circular domain) produces a reflection-symmetric ice thickness field.
# Sentinel fills (|H| > 9000) are masked as missing; cells where either side
# of the comparison is missing are excluded from the mean.
#
# Usage:
#   julia --project tests/symcheck.jl path/to/yelmo2D.nc [more.nc ...]
#
# Dependencies: NCDatasets, Statistics, Printf (stdlib).
#   julia> import Pkg; Pkg.add("NCDatasets")
#
# Output: per-file shape, Hmax, ice cell count, and Linf/L1 of the
# asymmetry under each of three reflections, normalised by Hmax.

using NCDatasets, Statistics, Printf

const FILL_THRESHOLD = 9000.0   # |H| above this counts as a sentinel fill

function diag(path::AbstractString; label::AbstractString = path)
    if !isfile(path)
        println(@sprintf("%32s: MISSING (%s)", label, path))
        return
    end
    H = NCDataset(path, "r") do nc
        var = haskey(nc, "lithk") ? "lithk" : "H_ice"
        Array{Float64}(nc[var][:, :, end])      # last timestep, (xc, yc)
    end
    H[abs.(H) .> FILL_THRESHOLD] .= NaN          # treat sentinel fills as missing
    Hmax = maximum(filter(!isnan, abs.(H)); init = 0.0)
    if Hmax == 0
        println(@sprintf("%s: all-zero H_ice", label)); return
    end
    nice = count(!isnan, H)

    function stats(d)
        a = abs.(d) ./ Hmax
        finite = .!isnan.(a)
        any(finite) || return (NaN, NaN)
        (maximum(a[finite]), mean(a[finite]))
    end

    H_lr  = reverse(H, dims = 1)                 # (xc, yc) so dim 1 is x
    H_tb  = reverse(H, dims = 2)
    H_rot = reverse(reverse(H, dims = 1), dims = 2)
    lr_inf, lr_l1   = stats(H .- H_lr)
    tb_inf, tb_l1   = stats(H .- H_tb)
    rot_inf, rot_l1 = stats(H .- H_rot)

    Hmean = mean(filter(!isnan, H))
    println(label)
    println(@sprintf("  shape=%s  Hmax=%8.2f  Hmean=%8.2f  ice_cells=%d",
                     string(size(H)), Hmax, Hmean, nice))
    println(@sprintf("  asym L-R   : Linf=%.3e  L1=%.3e", lr_inf, lr_l1))
    println(@sprintf("  asym T-B   : Linf=%.3e  L1=%.3e", tb_inf, tb_l1))
    println(@sprintf("  asym 180rot: Linf=%.3e  L1=%.3e", rot_inf, rot_l1))
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        println(stderr, "Usage: julia tests/symcheck.jl path/to/yelmo2D.nc [more.nc ...]")
        exit(2)
    end
    for arg in ARGS
        # Match the Python script's "..../tmp/<...>" label trimming when present.
        label = occursin("/tmp/", arg) ? split(arg, "/tmp/")[end] : arg
        diag(arg; label = label)
    end
end
