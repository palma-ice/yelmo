const yelmolib = "../libyelmo/include/libyelmo_c_api.so"

function yelmo_init(filename::String, grid_def::String, time::Float64)

    # Call yelmo_init fortran function
    ccall((:yelmo_init, yelmolib), Cvoid,
        (Ptr{UInt8}, Ptr{UInt8}, Float64),
        filename * "\0", grid_def * "\0", time)
    
    # Populate Julia version of Yelmo object with info from fortran
    g = yelmo_get_grid_info()

    # Load variable meta information
    v = (
        bnd = parse_variable_table("input/yelmo-variables-ybound.md"),
        dta = parse_variable_table("input/yelmo-variables-ydata.md"),
        dyn = parse_variable_table("input/yelmo-variables-ydyn.md"),
        mat = parse_variable_table("input/yelmo-variables-ymat.md"),
        thrm = parse_variable_table("input/yelmo-variables-ytherm.md"),
        tpo = parse_variable_table("input/yelmo-variables-ytopo.md"),
    )

    bnd = yelmo_get_variable_set(v.bnd,"bnd",g.nx,g.ny,g.nz_aa,g.nz_ac,g.nzr_aa,g.nzr_ac)
    dta = yelmo_get_variable_set(v.dta,"dta",g.nx,g.ny,g.nz_aa,g.nz_ac,g.nzr_aa,g.nzr_ac)
    dyn = yelmo_get_variable_set(v.dyn,"dyn",g.nx,g.ny,g.nz_aa,g.nz_ac,g.nzr_aa,g.nzr_ac)
    mat = yelmo_get_variable_set(v.mat,"mat",g.nx,g.ny,g.nz_aa,g.nz_ac,g.nzr_aa,g.nzr_ac)
    thrm = yelmo_get_variable_set(v.thrm,"thrm",g.nx,g.ny,g.nz_aa,g.nz_ac,g.nzr_aa,g.nzr_ac)
    tpo = yelmo_get_variable_set(v.tpo,"tpo",g.nx,g.ny,g.nz_aa,g.nz_ac,g.nzr_aa,g.nzr_ac)
    
    # tpo = (
    #     H = yelmo_get_var2D(g.nx, g.ny, "tpo_H_ice"),
    #     zs = yelmo_get_var2D(g.nx, g.ny, "tpo_z_srf")
    # )

    # thrm = (
    #     T = yelmo_get_var3D(g.nx, g.ny, g.nz_aa, "thrm_T_ice"),
    # )

    # dyn = (
    #     uxy = yelmo_get_var3D(g.nx, g.ny, g.nz_aa, "dyn_uxy"),
    #     uxy_s = yelmo_get_var2D(g.nx, g.ny, "dyn_uxy_s")
    # )

    return (g=g, v=v, tpo=tpo, thrm=thrm, dyn=dyn, bnd=bnd)
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
    nzr_aa = Ref{Cint}(0)
    nzr_ac = Ref{Cint}(0)

    ccall((:yelmo_get_grid_sizes, yelmolib), Cvoid,
          (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}),
          nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac)

    # Step 2: allocate buffers
    xc      = Vector{Cdouble}(undef, nx[])
    yc      = Vector{Cdouble}(undef, ny[])
    zeta_aa = Vector{Cdouble}(undef, nz_aa[])
    zeta_ac = Vector{Cdouble}(undef, nz_ac[])
    zeta_r_aa = Vector{Cdouble}(undef, nzr_aa[])
    zeta_r_ac = Vector{Cdouble}(undef, nzr_ac[])

    # Step 3: fill buffers
    ccall((:yelmo_get_grid_info, yelmolib), Cvoid,
          (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          xc, yc, zeta_aa, zeta_ac, zeta_r_aa, zeta_r_ac)

    return (
        nx    = Int(nx[]),
        ny    = Int(ny[]),
        nz_aa = Int(nz_aa[]),
        nz_ac = Int(nz_ac[]),
        nzr_aa = Int(nzr_aa[]),
        nzr_ac = Int(nzr_ac[]),
        xc      = xc,
        yc      = yc,
        zeta_aa = zeta_aa,
        zeta_ac = zeta_ac,
        zeta_r_aa = zeta_r_aa,
        zeta_r_ac = zeta_r_ac,
    )

end

function yelmo_get_variable_set(vlist, prefix, nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac)
    NamedTuple{keys(vlist)}(
        _load_var(vlist[k], prefix, k, nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac) for k in keys(vlist)
    )
end

function _load_var(meta, prefix, k, nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac)
    varname = "$(prefix)_$(k)"
    dims = meta.dimensions
    if :zeta in dims
        return yelmo_get_var3D(nx, ny, nz_aa, varname)
    elseif :zeta_ac in dims
        return yelmo_get_var3D(nx, ny, nz_ac, varname)
    elseif :zeta_rock in dims
        return yelmo_get_var3D(nx, ny, nzr_aa, varname)
    elseif :zeta_rock_ac in dims
        return yelmo_get_var3D(nx, ny, nzr_ac, varname)
    else
        return yelmo_get_var2D(nx, ny, varname)
    end
end

function yelmo_update_variable_set!(dat, vlist, prefix)
    for k in keys(vlist)
        _update_var!(dat[k], vlist[k], prefix, k)
    end
    return dat
end

function _update_var!(arr, meta, prefix, k)
    varname = "$(prefix)_$(k)"
    dims = meta.dimensions
    if :zeta_aa in dims
        yelmo_get_var3D!(arr, varname)
    elseif :zeta_ac in dims
        yelmo_get_var3D!(arr, varname)
    else
        yelmo_get_var2D!(arr, varname)
    end
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

#### Yelmo variables ####

# --- Types ---

struct VariableMeta
    id         :: Int
    dimensions :: Tuple{Vararg{Symbol}}
    units      :: String
    long_name  :: String
end

# Pretty printing
Base.show(io::IO, v::VariableMeta) =
    print(io, "VariableMeta(\"$(v.long_name)\" [$(v.units)], dims=$(v.dimensions))")

# --- Parser ---

"""
    parse_variable_table(markdown::AbstractString) -> NamedTuple

Parse a Markdown variable table into a NamedTuple of `VariableMeta` entries,
keyed by variable name (as Symbols). The table must have columns:
`id | variable | dimensions | units | long_name`.

# Example
```julia
meta = parse_variable_table(raw_table_string)
meta.T_srf        # => VariableMeta("Surface temperature" [K], dims=(:xc, :yc))
meta.T_srf.units  # => "K"
```
"""
function parse_variable_table(filename::AbstractString)

    # Read the file
    markdown = read(filename, String)

    rows = VariableMeta[]
    names = Symbol[]

    for line in eachline(IOBuffer(markdown))
        # Skip header and separator lines
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "| id") && continue
        startswith(stripped, "|-")   && continue
        startswith(stripped, "# ")   && continue
        startswith(stripped, '|') || continue

        # Split on '|', strip whitespace, drop empty edge tokens
        cols = strip.(split(stripped, '|'))
        filter!(!isempty, cols)
        length(cols) >= 5 || continue

        id        = parse(Int, cols[1])
        varname   = Symbol(strip(cols[2]))
        dims      = Tuple(Symbol(strip(d)) for d in split(cols[3], ','))
        units     = strip(cols[4])
        long_name = strip(cols[5])

        push!(names, varname)
        push!(rows, VariableMeta(id, dims, units, long_name))
    end

    return NamedTuple{Tuple(names)}(Tuple(rows))
end

##################################################################################
