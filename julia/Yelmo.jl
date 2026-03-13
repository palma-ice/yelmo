const yelmolib = "../libyelmo/include/libyelmo_c_api.so"

function yelmo_init(filename::String, grid_def::String, time::Float64)
    ccall((:yelmo_init, yelmolib), Cvoid,
        (Ptr{UInt8}, Ptr{UInt8}, Float64),
        filename * "\0", grid_def * "\0", time)
end

function yelmo_init_state(time::Float64, thrm_method::String)
    ccall((:yelmo_init_state, yelmolib), Cvoid,
        (Float64, Ptr{UInt8}),
        time, thrm_method * "\0")
end

function yelmo_step(time::Float64)
    ccall((:yelmo_step, yelmolib), Cvoid, (Float64,), time)
end

function yelmo_get_H_ice(nx::Int, ny::Int)
    H = Matrix{Float64}(undef, nx, ny)
    ccall((:yelmo_get_H_ice, yelmolib), Cvoid,
        (Ptr{Float64}, Int32, Int32),
        H, Int32(nx), Int32(ny))
    return H
end

function yelmo_set_bmb(bmb::Matrix{Float64})
    nx, ny = size(bmb)
    ccall((:yelmo_set_bmb, yelmolib), Cvoid,
        (Ptr{Float64}, Int32, Int32),
        bmb, Int32(nx), Int32(ny))
end