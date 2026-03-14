const yelmolib = "../libyelmo/include/libyelmo_c_api.so"

const VERTICAL_DIMS = (:zeta, :zeta_ac, :zeta_rock, :zeta_rock_ac)

mutable struct Yelmo
    alias::String
    time::Float64
    g::NamedTuple
    v::NamedTuple
    bnd::NamedTuple
    dta::NamedTuple
    dyn::NamedTuple
    mat::NamedTuple
    thrm::NamedTuple
    tpo::NamedTuple
end

function Yelmo(filename::String, grid_def::String, time::Float64; alias::String="ylmo1")

    # First call yelmo init to initialize model in fortran
    ccall((:yelmo_init, yelmolib), Cvoid,
        (Ptr{UInt8}, Ptr{UInt8}, Float64, Ptr{UInt8}),
        filename * "\0", grid_def * "\0", time, alias * "\0")

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
    
    return Yelmo(alias,time,g,v,bnd,dta,dyn,mat,thrm,tpo)
end

function init_state!(ylmo::Yelmo, time::Float64, thrm_method::String; alias::String="ylmo1")
    
    # call yelmo_init_state in fortran
    ccall((:yelmo_init_state, yelmolib), Cvoid,
        (Float64, Ptr{UInt8}, Ptr{UInt8}),
        time, thrm_method * "\0", alias * "\0")
    
    # Update yelmo in julia
    yelmo_get_variables!(ylmo)

    # Update time
    ylmo.time = time

    return ylmo
end

function time_step!(ylmo::Yelmo, dt::Float64; alias::String="ylmo1")

    # Update time
    ylmo.time += dt
    
    # Call yelmo_step in fortran
    ccall((:yelmo_step, yelmolib), Cvoid, (Float64, Ptr{UInt8}), ylmo.time, alias * "\0")

    # Update yelmo in julia
    yelmo_get_variables!(ylmo)

    return ylmo
end

function yelmo_get_grid_info(; alias::String="ylmo1")

    # Step 1: get sizes
    nx    = Ref{Cint}(0)
    ny    = Ref{Cint}(0)
    nz_aa = Ref{Cint}(0)
    nz_ac = Ref{Cint}(0)
    nzr_aa = Ref{Cint}(0)
    nzr_ac = Ref{Cint}(0)

    ccall((:yelmo_get_grid_sizes, yelmolib), Cvoid,
          (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ptr{UInt8}),
          nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac, alias * "\0")

    # Step 2: allocate buffers
    xc      = Vector{Cdouble}(undef, nx[])
    yc      = Vector{Cdouble}(undef, ny[])
    zeta_aa = Vector{Cdouble}(undef, nz_aa[])
    zeta_ac = Vector{Cdouble}(undef, nz_ac[])
    zeta_r_aa = Vector{Cdouble}(undef, nzr_aa[])
    zeta_r_ac = Vector{Cdouble}(undef, nzr_ac[])

    # Step 3: fill buffers
    ccall((:yelmo_get_grid_info, yelmolib), Cvoid,
          (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}),
          xc, yc, zeta_aa, zeta_ac, zeta_r_aa, zeta_r_ac, alias * "\0")

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

function yelmo_get_variables!(ylmo)

    yelmo_get_variable_set!(ylmo.bnd, ylmo.v.bnd, "bnd")
    yelmo_get_variable_set!(ylmo.dta, ylmo.v.dta, "dta")
    yelmo_get_variable_set!(ylmo.dyn, ylmo.v.dyn, "dyn")
    yelmo_get_variable_set!(ylmo.mat, ylmo.v.mat, "mat")
    yelmo_get_variable_set!(ylmo.thrm, ylmo.v.thrm, "thrm")
    yelmo_get_variable_set!(ylmo.tpo, ylmo.v.tpo, "tpo")
    
    return ylmo
end

function yelmo_get_variable_set!(dat, vlist, prefix)
    for k in keys(vlist)
        _get_var!(dat[k], vlist[k], prefix, k)
    end
    return dat
end

function _get_var!(arr, meta, prefix, k)
    varname = "$(prefix)_$(k)"
    dims = meta.dimensions
    is3D = any(d -> d in VERTICAL_DIMS, dims)
    if is3D
        yelmo_get_var3D!(arr, varname)
    else
        yelmo_get_var2D!(arr, varname)
    end
end

function yelmo_get_variable_set(vlist, prefix, nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac)
    dat = NamedTuple{keys(vlist)}(
        _alloc_var(vlist[k], nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac) for k in keys(vlist)
    )
    yelmo_get_variable_set!(dat, vlist, prefix)
    return dat
end

function _alloc_var(meta, nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac)
    dims = meta.dimensions
    nz = nothing
    for (d, n) in ((:zeta, nz_aa), (:zeta_ac, nz_ac), (:zeta_rock, nzr_aa), (:zeta_rock_ac, nzr_ac))
        d in dims && (nz = n; break)
    end
    return isnothing(nz) ? zeros(Float64, nx, ny) : zeros(Float64, nx, ny, nz)
end

function yelmo_get_var2D!(v2D::Array{Float64,2}, name::String; alias::String="ylmo1")
    nx, ny = size(v2D)
    ccall((:yelmo_get_var2D, yelmolib), Cvoid,
        (Ptr{Float64}, Int32, Int32, Ptr{UInt8}, Ptr{UInt8}),
        v2D, Int32(nx), Int32(ny), name * "\0", alias * "\0")
    return v2D
end

function yelmo_get_var2D(nx::Int, ny::Int, name::String; alias::String="ylmo1")
    v2D = Matrix{Float64}(undef, nx, ny)
    yelmo_get_var2D!(v2D, name; alias)
    return v2D
end

function yelmo_get_var3D!(v3D::Array{Float64,3}, name::String; alias::String="ylmo1")
    nx, ny, nz = size(v3D)
    ccall((:yelmo_get_var3D, yelmolib), Cvoid,
        (Ptr{Float64}, Int32, Int32, Int32, Ptr{UInt8}, Ptr{UInt8}),
        v3D, Int32(nx), Int32(ny), Int32(nz), name * "\0", alias * "\0")
    return v3D
end

function yelmo_get_var3D(nx::Int, ny::Int, nz::Int, name::String; alias::String="ylmo1")
    v3D = Array{Float64}(undef, nx, ny, nz)
    yelmo_get_var3D!(v3D, name, alias)
    return v3D
end

function yelmo_set_var2D!(name::String, v2D::Array{Float64,2}; alias::String="ylmo1")

    nx, ny = size(v2D)

    ccall((:yelmo_set_var2D, yelmolib), Cvoid,
          (Ptr{Cdouble}, Cint, Cint, Cstring, Ptr{UInt8}),
          v2D, nx, ny, name, alias)

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
