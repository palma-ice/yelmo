module YelmoIO

using Oceananigans
using NCDatasets

export load_grids_from_restart
export load_field_from_dataset_2D
export load_field_from_dataset_3D
export load_fields_from_restart

function matches_patterns(name, patterns)
    any(p -> occursin(p, name), patterns)
end

function load_grids_from_restart(filename)

    # Open the NetCDF file
    ds = NCDataset(filename)

    xc = ds["xc"][:]
    yc = ds["yc"][:]
    
    Nx = length(xc)
    Ny = length(yc)
    xlims = extrema(xc)
    ylims = extrema(yc)

    # Expand x- and ylims so that nodes represent centers of cells
    dx = xc[2]-xc[1]
    dy = yc[2]-yc[1]
    xlims = (xlims[1] - dx/2, xlims[2] + dx/2)
    ylims = (ylims[1] - dy/2, ylims[2] + dy/2)
    
    # First, create your grid
    grid2d = RectilinearGrid(
        size = (Nx, Ny),
        x = xlims,
        y = ylims,
        topology = (Bounded, Bounded,Flat)
    )

    # Create 3D grid if zeta exists
    grid3d = nothing
    if haskey(ds, "zeta")
        zeta = ds["zeta"][:]
        zeta_ac = ds["zeta_ac"][:]
        Nz = length(zeta)
        zlims = extrema(zeta)
        
        grid3d = RectilinearGrid(
            size = (Nx, Ny, Nz),
            x = xlims,
            y = ylims,
            #z = zlims,
            z = zeta_ac,
            topology = (Bounded, Bounded, Bounded)
        )
    end

    close(ds)

    return grid2d, grid3d
end

function load_field_from_dataset_2D(ds::NCDataset, varname::Union{AbstractString,Symbol}, grid2d)

    varname = String(varname)

    xface_variables = ["ux_s","ux_b", r".*_acx$"]
    yface_variables = ["uy_s","uy_b", r".*_acy$"]
    zface_variables = []

    # Create the field and store in interior of field
    if matches_patterns(varname, xface_variables)
        field = XFaceField(grid2d)
        interior(field)[2:end,:] .= ds[varname][:,:]
        interior(field)[1,:] .= ds[varname][1,:]
    elseif matches_patterns(varname, yface_variables)
        field = YFaceField(grid2d)
        interior(field)[:,2:end] .= ds[varname][:,:]
        interior(field)[:,1] .= ds[varname][:,1]
    elseif matches_patterns(varname, zface_variables)
        # Maybe this case doesn't exist in 2d?
        field = ZFaceField(grid2d)
    else
        field = CenterField(grid2d)
        interior(field) .= ds[varname][:,:]
    end

    #interior(field) .= ds[varname][:,:]

    return (field)
end

function load_field_from_dataset_3D(ds::NCDataset, varname::Union{AbstractString,Symbol}, grid3d)

    varname = String(varname)

    xface_variables = ["ux", r".*_acx$"]
    yface_variables = ["uy", r".*_acy$"]
    zface_variables = ["uz","uz_star","jvel_dzx","jvel_dzy","jvel_dzz"]

    # Create the field and store in interior of field
    if matches_patterns(varname, xface_variables)
        field = XFaceField(grid3d)
        interior(field)[2:end,:,:] .= ds[varname][:,:,:]
        interior(field)[1,:,:] .= ds[varname][1,:,:]
    elseif matches_patterns(varname, yface_variables)
        field = YFaceField(grid3d)
        interior(field)[:,2:end,:] .= ds[varname][:,:,:]
        interior(field)[:,1,:] .= ds[varname][:,1,:]
    elseif matches_patterns(varname, zface_variables)
        field = ZFaceField(grid3d)
        interior(field)[:,:,:] .= ds[varname][:,:,:]
    else
        field = CenterField(grid3d)
        interior(field) .= ds[varname][:,:,:]
    end

    return (field)
end

function load_field_from_dataset_2D(filename, varname::Union{AbstractString,Symbol}, grid2d)
    ds = NCDataset(filename)
    field = load_field_from_dataset_2D(ds,varname,grid2d)
    close(ds)
    return field
end

function load_field_from_dataset_3D(filename::AbstractString, varname::Union{AbstractString,Symbol}, grid3d)
    ds = NCDataset(filename)
    field = load_field_from_dataset_3D(ds,varname,grid3d)
    close(ds)
    return field
end

function load_fields_from_restart(filename,grid2d,grid3d)

    # Open the NetCDF file
    ds = NCDataset(filename)

    # Load variables into dictionary
    dat = Dict{String, Field}()
    dat = Dict{String, Field}()
    
    # Skip several variables for now
    variables_to_skip = [
        r"^jvel.*",
        r"^strs.*",
        r"^strn.*",
    ]

    # Also limit 3D variables to the following...
    variables_to_load_3d = [
        r"^ux.*",
        r"^uy.*",
        r"^uz.*",
        r"^T.*",
    ]
            
    for varname in keys(ds)

        # Skip variables as needed
        if matches_patterns(varname,variables_to_skip)
            continue
        end

        # Get current dimension names
        dimnames_now = dimnames(ds[varname])

        # Skip variables that do not at least contain 2D information
        if length(dimnames_now) < 2
            continue
        end

        if dimnames_now[1:2] != ("xc", "yc")
            continue
        end
        
        dims = size(ds[varname])
        
        println("$varname : $dims")

        if length(dims) == 2 || dimnames_now[3] == "time"
            # 2D variable
            dat[varname] = load_field_from_dataset_2D(ds,varname,grid2d)

        elseif (dimnames_now[3] == "zeta" || dimnames_now[3] == "zeta_ac") && grid3d !== nothing
            # Only load some variables to avoid segfault...

            if matches_patterns(varname,variables_to_load_3d)
                # 3D variable on zeta grid
                dat[varname] = load_field_from_dataset_3D(ds,varname,grid3d)
            end
        end

    end
    
    close(ds)

    return dat
end

end