const yelmolib = "../libyelmo/include/libyelmo_c_api.so"

function yelmo_init(filename::String, grid_def::String, time::Float64)
    ccall((:yelmo_init, yelmolib), Cvoid,
        (Ptr{UInt8}, Ptr{UInt8}, Float64),
        filename * "\0", grid_def * "\0", time)
    
    # Return Yelmo object with current state
    g = yelmo_get_grid_info()

    bnd = (
        zb = yelmo_get_var2D(g.nx, g.ny, "tpo_z_bed"),
    )

    tpo = (
        H = yelmo_get_var2D(g.nx, g.ny, "tpo_H_ice"),
        zs = yelmo_get_var2D(g.nx, g.ny, "tpo_z_srf")
    )

    thrm = (
        T = yelmo_get_var3D(g.nx, g.ny, g.nz_aa, "thrm_T_ice"),
    )

    dyn = (
        uxy = yelmo_get_var3D(g.nx, g.ny, g.nz_aa, "dyn_uxy"),
        uxy_s = yelmo_get_var2D(g.nx, g.ny, "dyn_uxy_s")
    )

    return (g=g, tpo=tpo, thrm=thrm, dyn=dyn, bnd=bnd)
end

function yelmo_init_state(time::Float64, thrm_method::String)
    ccall((:yelmo_init_state, yelmolib), Cvoid,
        (Float64, Ptr{UInt8}),
        time, thrm_method * "\0")
end

function yelmo_step(time::Float64)
    ccall((:yelmo_step, yelmolib), Cvoid, (Float64,), time)
end

function yelmo_get_grid_info()

    # Step 1: get sizes
    nx    = Ref{Cint}(0)
    ny    = Ref{Cint}(0)
    nz_aa = Ref{Cint}(0)
    nz_ac = Ref{Cint}(0)

    ccall((:yelmo_get_grid_sizes, yelmolib), Cvoid,
          (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}),
          nx, ny, nz_aa, nz_ac)

    # Step 2: allocate buffers
    xc      = Vector{Cdouble}(undef, nx[])
    yc      = Vector{Cdouble}(undef, ny[])
    zeta_aa = Vector{Cdouble}(undef, nz_aa[])
    zeta_ac = Vector{Cdouble}(undef, nz_ac[])

    # Step 3: fill buffers
    ccall((:yelmo_get_grid_info, yelmolib), Cvoid,
          (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          xc, yc, zeta_aa, zeta_ac)

    return (
        nx    = Int(nx[]),
        ny    = Int(ny[]),
        nz_aa = Int(nz_aa[]),
        nz_ac = Int(nz_ac[]),
        xc      = xc,
        yc      = yc,
        zeta_aa = zeta_aa,
        zeta_ac = zeta_ac,
    )

end

function yelmo_get_var2D!(v2D::Array{Float64,2}, name::String)
    nx, ny = size(v2D)
    ccall((:yelmo_get_var2D, yelmolib), Cvoid,
        (Ptr{Float64}, Int32, Int32, Ptr{UInt8}),
        v2D, Int32(nx), Int32(ny), name * "\0")
    return v2D
end

function yelmo_get_var2D(nx::Int, ny::Int, name::String)
    v2D = Matrix{Float64}(undef, nx, ny)
    yelmo_get_var2D!(v2D, name)
    return v2D
end

function yelmo_get_var3D!(v3D::Array{Float64,3}, name::String)
    nx, ny, nz = size(v3D)
    ccall((:yelmo_get_var3D, yelmolib), Cvoid,
        (Ptr{Float64}, Int32, Int32, Int32, Ptr{UInt8}),
        v3D, Int32(nx), Int32(ny), Int32(nz), name * "\0")
    return v3D
end

function yelmo_get_var3D(nx::Int, ny::Int, nz::Int, name::String)
    v3D = Array{Float64}(undef, nx, ny, nz)
    yelmo_get_var3D!(v3D, name)
    return v3D
end

function yelmo_set_var2D!(name::String, v2D::Array{Float64,2})

    nx, ny = size(v2D)

    ccall((:yelmo_set_var2D, yelmolib), Cvoid,
          (Ptr{Cdouble}, Cint, Cint, Cstring),
          v2D, nx, ny, name)

    return nothing

end