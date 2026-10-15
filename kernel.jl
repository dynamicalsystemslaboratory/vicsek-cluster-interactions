
include( "args.jl" )
include( "master.jl" )
include( "prob.jl" )

# Homogenous kernel functions.
function collision(ℓ::Vector{defFloat}, i::defInt, j::defInt; c0::defFloat=1.0, γ::defFloat=1.5)::defFloat
    return c0*(ℓ[i] + ℓ[j])^γ
end

function fragment(A::Vector{defFloat}, i::defInt; f0::defFloat=1.0, γ::defFloat=0.0)::defFloat
    return f0*A[i]^γ
end


# Create kernels from the calibration parameters.
function createkernels(N::defInt, A::Vector{defFloat}, ℓ::Vector{defFloat}, γ::Vector{defFloat}; f0::defFloat=100.0)::Tuple{Function,Function}
    # Construct coagulation kernel functions.
    coll(i::defInt, j::defInt)::defFloat = collision( ℓ, i, j; c0=γ[1], γ=1.0 )#γ[2] )

    # Construct fragmentation kernel functions.
    frag(i::defInt)::defFloat = fragment( A, i; f0=f0, γ=0.0 )#, γ=θ[4] )

    # Return function variables.
    return coll, frag
end

function kernelloss(F::Vector{defInt}, γ::Vector{defFloat}, data;
    T::defFloat=10.0, n::Union{Nothing,defInt}=nothing, verbose::Bool=false)
    # Unpack data variables.
    N, Ω, A, ℓ = data
    n = isnothing( n ) ? N : n

    # Define local kernels.
    coll, frag = createkernels( N, A, ℓ, γ )
    params = Params( N, Ω, coll, frag )

    # Compute the exponent related to F and initialize.
    λ = 10^optimize( createzipfloss( N, F ), [log10( 1.5 )] ).minimizer[1]

    # Simulate the master equations.
    f̂ = zipf.( 1:N; λ=λ, b=N )
    f = logsimulate( params; f=f̂, T=T )

    # Return the negative log-likelihood.
    F̂ = round.( defInt, sum( F )*f̂ )
    NLL = nll( f[1:n]./sum( f[1:n] ), F̂[1:n] )
    return verbose ? (NLL, f) : NLL
end

function createkernelloss(F::Vector{defInt}, N::defInt, Ω::Matrix{defFloat},
    A::Vector{defFloat}, ℓ::Vector{defFloat}; T::defFloat=10.0, n::Union{Nothing,defInt}=nothing)::Function
    # Initialize the loss function.
    loss(γ̂::Vector{defFloat}) = kernelloss( F, 10.0.^γ̂, (N, Ω, A, ℓ); T=T, n=n )

    # Return function.
    return loss
end

println( "Loaded kernel.jl class file." )
