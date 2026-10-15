
using Base.Threads
println( "Number of threads: ", nthreads() )

include( "network.jl" )

# Global default interaction radius and time-step.
const Δt = defFloat( 1 )

# Update the full list of agents in the Vicsek model.
function vicsek(params::Dict{Char,Real}, A::SparseMatrixCSC{defInt,defInt},
    xlist::Vector{defFloat}, ylist::Vector{defFloat}, θlist::Vector{defFloat};
    δt::defFloat=Δt, type::Symbol=:forward)::Tuple{Vector{defFloat}, Vector{defFloat}, Vector{defFloat}}
    # Unpack parameters.
    N = Int( params['N'] )
    L = defFloat( params['L'] )
    v = defFloat( params['v'] )
    η = defFloat( params['η'] )

    # Parallel orientation update
    θnew = similar( θlist )
    @threads for i in 1:N
        # Initialize summand.
        sinθ = defFloat( 0.0 )
        cosθ = defFloat( 0.0 )

        # Loop through all neighbors (including particle i itself)
        @inbounds for idx in A.colptr[i]:(A.colptr[i+1] - 1)
            j = A.rowval[idx]
            sinθ = sinθ + sin( θlist[j] )
            cosθ = cosθ + cos( θlist[j] )
        end
        θbar = atan( sinθ, cosθ )

        # Add uniform noise in [-ηπ, ηπ] and wrap to [-π, π)
        noise = (η*π)*(2.0*rand( defFloat ) - 1.0)
        @inbounds θnew[i] = mod( θbar + noise + π, 2π ) - π
    end

    # Compute new positions using the specified update method.
    θn = (type === :forward) ? θnew : θlist
    xnew = @. mod( xlist + v*cos( θn )*δt, L )
    ynew = @. mod( ylist + v*sin( θn )*δt, L )

    # Return the new state (wrapped to the boundaries).
    return xnew, ynew, θnew
end

# Compute the alignment (order) of particles in the system.
function order(N, θlist::Vector{defFloat})::defFloat
    r = sum( exp.(im .* θlist) ) / N
    return abs( r )
end

# Compute the critical noise from order parameter data.
function criticalnoise(xlist::AbstractArray{defFloat}, ylist::AbstractArray{defFloat}; φc::defFloat=0.5)::defFloat
    # Append zero-noise order combination to lists.
    ηlist = [0.0; xlist; 1.0]
    φlist = [1.0; ylist; 0.0]

    # Extract upper and lower nearest points.
    ip = argmin( φlist[φc .< φlist] )
    id = argmax( φlist[φlist .< φc] )

    # Extract the upper and lower (noise,order) pair.
    ηp = ηlist[φc .< φlist][ip];  φp = φlist[φc .< φlist][ip]
    ηd = ηlist[φlist .< φc][id];  φd = φlist[φlist .< φc][id]

    # Compute the slope and intercept of the line between the points.
    m = (φd - φp)/(ηd - ηp)
    b = φp - m*ηp

    # Return critical point estimate.
    return (φc - b)/m
end

println( "Loaded model.jl class file." )
