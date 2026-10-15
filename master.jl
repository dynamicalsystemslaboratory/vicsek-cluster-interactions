
include( "args.jl" )

using DifferentialEquations

# Helper structure for simulating the master equations.
mutable struct Params
    N::defInt
    Ω::Matrix{defFloat}
    coll::Function
    frag::Function
end

# Compute the steady state CMD from the master equations.
function clustermassdistr(n::Vector{defFloat}, params::Params)::Vector{defFloat}
    Δn = sum( (1:params.N).*n )
    return (1:params.N).*n./Δn
end

# Compute the steady state CSD from the master equations.
function clustersizedistr(n::Vector{defFloat}, params::Params)::Vector{defFloat}
    Δn = sum( n )
    return n./Δn
end

# Master equation for cluster collision and fragmentation dynamics.
function model!(ṅ, n, params::Params, t::defFloat)
    # Reset the derivatives.
    fill!( ṅ, 0.0 )

    # For each potential cluster size.
    @threads for i ∈ 1:params.N
        # Cluster j > i fragments into clusters of size i and j-i.
        ṅ[i] += i < params.N ? sum( params.frag( j )*params.Ω[j,i]*n[j] for j ∈ i+1:params.N )          : 0.0

        # Cluster of size i fragments into clusters j and i-j.
        ṅ[i] -= i > 1        ? params.frag( i )/2*sum( params.Ω[i,j] for j ∈ 1:i-1 )*n[i]               : 0.0

        # Cluster of size j < i collides with cluster i-j to create i.
        ṅ[i] += i > 1        ? sum( params.coll( j,i-j )*j/i*params.Ω[j,i]*n[j]*n[i-j] for j ∈ 1:i-1 )  : 0.0

        # Cluster of size i collides with cluster j ≤ N-i to create i+j.
        ṅ[i] -= i < params.N ? sum( params.coll( i,j )*params.Ω[i,i+j]*n[i]*n[j] for j ∈ 1:params.N-i ) : 0.0
    end

    # Optional conservation of mass check. Good when simulating small systems.
    # Σṅ = sum( ṅ )
    # @assert Σṅ < 1e-6 "Conservation of mass is not preserved ($(round( Σṅ, digits=7 )))."
end

# Master equations reparameterized to a log-space.
function logmodel!(ṁ, m, params::Params, t)
    # Reset the derivatives.
    fill!( ṁ, 0.0 )

    # For each potential cluster size.
    @threads for i ∈ 1:params.N
        # Cluster j > i fragments into clusters of size i and j-i.
        ṁ[i] += i < params.N ?
            sum( params.frag( j )*params.Ω[j,i]*exp( m[j] - m[i] ) for j ∈ i+1:params.N ) : 0.0

        # Cluster of size i fragments into clusters j and i-j.
        ṁ[i] -= i > 1 ?
            params.frag( i )/2*sum( params.Ω[i,j] for j ∈ 1:i-1 ) : 0.0

        # Cluster of size j < i collides with cluster i-j to create i.
        ṁ[i] += i > 1 ?
            sum( (j/i)*params.coll( j,i-j )*params.Ω[j,i]*exp( m[j] + m[i-j] - m[i] ) for j ∈ 1:i-1 )  : 0.0

        # Cluster of size i collides with cluster j ≤ N-i to create i+j.
        ṁ[i] -= i < params.N ?
            sum( params.coll( i,j )*params.Ω[i,i+j]*exp( m[j] ) for j ∈ 1:params.N-i ) : 0.0
    end
end

# Simulate the model and return the cluster mass distribution.
function simulate(params::Params;
    f::Vector{defFloat}=fill( 1/params.N, params.N ),
    T::Real=10.0, verbose::Bool=false)

    # Initial conditions (from the given initial distribution).
    n0 = params.N/sum( i*f[i] for i ∈ 1:params.N )*f

    # Build the ODE problem and solve.
    problem = ODEProblem( model!, n0, T, params )
    results = solve( problem, Vern9(); isoutofdomain=(n,p,t) -> any( n .< 0 ),
        save_everystep=false, save_start=false, save_end=true )

    # Unpack the final state.
    n = results.u[end]

    # Return the final step of the simulation.
    Δn = clustermassdistr( n, params )
    return verbose ? (results, Δn) : Δn
end

# Simulate the model and return the cluster mass distribution.
function logsimulate(params::Params;
    f::Vector{defFloat}=fill( 1/params.N, params.N ),
    T::Real=10.0, verbose::Bool=false)
    # Initial conditions (from the given initial distribution).
    n0 = params.N/sum( i*f[i] for i ∈ 1:params.N )*f

    # Transform the initial conditions into log-space.
    m0 = log.( n0 )

    # Build the ODE problem and solve.
    problem = ODEProblem( logmodel!, m0, T, params )
    results = solve( problem, Vern9();  save_everystep=false, save_start=false, save_end=true )

    # Unpack the final state.
    m = results.u[end]
    n = exp.( m )

    # Return the final step of the simulation.
    Δn = clustersizedistr( n, params )
    return verbose ? (results, m, Δn) : Δn
end

println( "Loaded master.jl class file." )
