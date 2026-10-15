
include( "args.jl" )

using Optim

# Negative log-likelihood (from frequency data).
function nll(P::Vector{defFloat}, F::AbstractArray{<:Real})::defFloat
    ι = (0.0 .< P) .& (0 .< F)
    return sum( -F[ι].*log.( P[ι] ) )/sum( F[ι] )
end

# ---------------------- Distribution Helpers ---------------------- #
# Zipf distribution function.
function zipf(j::defInt; a::defInt=1, b::defInt=128, λ::defFloat=0.15, ξ::defFloat=Inf)::defFloat
    H = sum( 1.0/(i^λ*exp( (i/ξ)^2 )) for i ∈ a:b )
    return a ≤ j ≤ b ? 1.0/(H*j^λ*exp( (j/ξ)^2 )) : 0.0
end

# Double-Zipf distribution with contiguous transition.
function doublezipf(b::defInt, j::defInt; a::defInt=1, c::defInt=128, λ::defFloat=0.15, λ̂::defFloat=0.15)::defFloat
    a == c && return j == a ? 1.0 : 0.0

    # Compute the normalization constant.
    H1 = sum( 1/i^λ for i ∈ a:b )
    H2 = b == c ? 1.0 : sum( 1/i^λ̂ for i ∈ b+1:c )
    H = H1 + (b == c ? 0.0 : b^(λ̂ - λ)*H2)

    # Return the PMF value.
    (a ≤ j ≤ b) && return 1/(H*j^λ)
    (b+1 ≤ j ≤ c) && return b^(λ̂ - λ)/(H*j^λ̂)
    return 0.0
end

# Weibull distribution.
function weibull(j::defInt; a::defInt=1, b::defInt=128, λ::defFloat=0.15, ξ::defFloat=64.0)::defFloat
    H = exp( -((a/ξ)^λ) ) - exp( -(((b + 1)/ξ)^λ) )
    return a ≤ j ≤ b ? (exp( -(j/ξ)^λ ) - exp( -((j + 1)/ξ)^λ ))/H : 0.0
end

# Distribution(s) to describe fragmentations.
function P1(j::defInt; N=128, i::defInt=1, λ::defFloat=0.15)::defFloat
    return 1 ≤ j ≤ i-1 ? zipf( j; b=(i - 1), λ=λ ) : 0.0
end

function P2(j::defInt; N::defInt=128, i::defInt=1, λ::defFloat=0.15)::defFloat
    return 1 ≤ j ≤ i-1 ? zipf( i - j; b=(i - 1), λ=λ ) : 0.0
end

# Distribution(s) to describe collisions.
function P3(j::defInt; N::defInt=128, N̂::defInt=N, i::defInt=1, λr::defFloat=0.15, λr̂::defFloat=λr)::defFloat
    n̂ = N̂ > i ? N̂-i : 1
    return i+1 ≤ j ≤ N ? doublezipf( n̂, j - i; c=(N - i), λ=λr, λ̂=λr̂ ) : 0.0
end

# Combined distribution function.
function P(j::defInt; N::defInt=128, i::defInt=1, λℓ::defFloat=0.15, λr::defFloat=0.15, κ::defFloat=0.66)::defFloat
    i == 1 && (return P3( j; N=N, i=i, λr=λr ))
    i == N && (return 1/2*(P1( j; N=N, i=i, λ=λℓ ) + P2( j; N=N, i=i, λ=λℓ )))
    return (κ/2*(P1( j; N=N, i=i, λ=λℓ ) + P2( j; N=N, i=i, λ=λℓ ))
        + (1 - κ)*P3( j; N=N, i=i, λr=λr ))
    # return 0.0
end

# ---------------------- Optimization Helpers ---------------------- #
# MLE cost function.
function cost(θ::Vector{defFloat}, F::Matrix{defInt}; N::defInt=1024, i::defInt=1)::defFloat
    # Compute the probability distribution using the parameters.
    κ = i == 1 ? 0.0 : (i == N ? 1.0 : θ[3])
    f = P.( 1:N; N=N, i=i, λℓ=θ[1], λr=θ[2], κ=κ )

    # Compute the negative log-likelihood.
    return nll( f, F[i,:] )
end

# Cost function generator.
function convert(θ̂::Vector{defFloat})::Vector{defFloat}
    return [10.0^θ̂[1], 10.0^θ̂[2], 10.0^(-10.0^θ̂[3])]
end

function createloss(F::Matrix{defInt}, N::defInt, i::defInt)
    loss(θ̂::Vector{defFloat}) = cost( convert( θ̂ ), F; N=N, i=i )
    return loss
end

# Run optimization.
function runopt(i::defInt, F::Matrix{defInt}, N::defInt, θ0::Vector{defFloat})
    θ̂0 = log10.( [θ0[1], θ0[2], -log10( θ0[3] )] )
    return optimize( createloss( F, N, i ), θ̂0, NelderMead() )
end

# --------------- Parameter Extrapolation Helpers ------------------ #
function adj(N::defInt, i::defInt, β::defFloat)::defFloat
    j = i - 1;  N̂ = N - 1
    x̂ = 2j/N̂ - 1
    return x̂/(1 + β*N̂*(1 - x̂^2))
end

function arcsine(N::defInt, i::defInt, β::defFloat)::defFloat
    return asin( adj( N, i, β ) )/π + 1/2
end

function adjarcsine(N::defInt, i::defInt, α̂::Vector{defFloat})::defFloat
    α = exp.( α̂ );  #β = α[2:3]
    return arcsine( N, i, α[2] ) + α[1]*(1 - adj( N, i, α[2] )^2)
end

function estimatetrend(N::defInt, x::Vector{defFloat}, α0::Vector{defFloat};
    f::Function=adjarcsine, verbose::Bool=false
    )::Union{Vector{defFloat},Tuple{Vector{defFloat},Any}}
    # Construct cost function using the optimized rate parameters.
    N̄ = length( x )
    cost(α̂::Vector{defFloat}) = sum( (x .- [f( N, i, α̂ ) for i ∈ 1:N̄]).^2 )/N̄

    # Perform the hyper-optimization.
    results = optimize( cost, log.( α0 ) )

    # Return the estimated points.
    x̂ = [f( N, i, Optim.minimizer( results ) ) for i ∈ 1:N]
    return verbose ? (x̂, results) : x̂
end

function fitline(N::defInt, c::Vector{defInt}, λ::Vector{defFloat})
    # Extract the non-empty indeces.
    N̂ = length( λ )
    c1 = c[1];  c2 = c[2]
    ĉ2 = N == N̂ ? c2 : (N̂ < N - c2 ? 0 : N̂ - (N - c2))

    # Compile X and y vectors for linear fit.
    X = hcat( (c1:N̂-ĉ2)./N, ones( N̂ - (c1 + ĉ2) + 1 ) )
    y = λ[c1:end-ĉ2]

    # Return the linear fit.
    return X \ y
end

# --------------- Cluster Mass Distribution Helpers ---------------- #
function kld(p::Vector{defFloat}, q::Vector{defFloat})::defFloat
    # Make vectors equal length.
    Np = length( p );  Nq = length( q )
    p̂ = Np < Nq ? [p;zeros( Nq-Np )] : p;  p̂ = p̂./sum( p̂ )
    q̂ = Nq < Np ? [q;zeros( Np-Nq )] : q;  q̂ = q̂./sum( q̂ )

    # Compute the Kullback-Leibler distance.
    i = (0.0 .< p̂) .&  (0.0 .< q̂)
    return sum( p̂[i].*log2.( p̂[i]./q̂[i] ) )
end

function jsd(p::Vector{defFloat}, q::Vector{defFloat})::defFloat
    # Check that distributions are normalized.
    @assert isapprox( sum( p ), 1.0; atol=1e-3 ) "P distribution is not normalized to 1: $(round( sum( p ), digits=3 ))."
    @assert isapprox( sum( q ), 1.0; atol=1e-3 ) "Q distribution is not normalized to 1: $(round( sum( q ), digits=3 ))."

    # Make vectors equal length.
    Np = length( p );  Nq = length( q )
    p̂ = Np < Nq ? [p;zeros( Nq-Np )] : p
    q̂ = Nq < Np ? [q;zeros( Np-Nq )] : q

    # Compute the Jenson-Shannon divergence.
    m = 1/2*(p̂ .+ q̂)
    return (kld( p̂, m ) + kld( q̂, m ))/2
end

# Fitting the CMD helpers.
function zipfloss(N::defInt, F::AbstractArray{<:Real}, λ::defFloat)::defFloat
    # Compute the Zipf law.
    f = zipf.( 1:N; b=N, λ=λ )

    # Return the negative log-likelihood.
    return nll( f, F )
end

function createzipfloss(N::defInt, F::AbstractArray{<:Real})::Function
    loss(θ̂::Vector{defFloat})::defFloat = zipfloss( N, F, 10.0.^θ̂[1]  )
    return loss
end

# Fitting the CMD helpers.
function zipfcutloss(N::defInt, Δi::Vector{defInt}, λ::defFloat, ξ::defFloat; xmin::defInt=1)::defFloat
    # Compute the Zipf law.
    Δn̂ = zipf.( xmin:N; b=N, λ=λ, ξ=ξ )

    # Return the negative log-likelihood.
    Σn̂ = sum( Δn̂ )
    return nll( Δn̂./Σn̂, Δi[xmin:N] )
end

function createzipfcutloss(N::defInt, Δi::Vector{defInt}; λ::Union{Nothing,defFloat}=nothing, xmin::defInt=1)::Function
    loss(θ̂::Vector{defFloat})::defFloat = zipfcutloss( N, Δi, isnothing( λ ) ? 10.0^θ̂[2] : λ, N/(1 + 10.0^θ̂[1]); xmin=xmin )
    return loss
end

# Helper for constructing master equation components.
function importfittedmap(N::defInt, folder::String; N̂::defInt=N)::Matrix{defFloat}
    # Construct map from individual exponent results.
    θlist = readdlm( folder*"optimal-rates-hat.txt" )

    # Compute the transition map using the calibrated coefficients.
    Ω = [P( j; N=N̂, i=i, λℓ=θlist[i,1], λr=θlist[i,2], κ=θlist[i,3] ) for i in 1:N̂, j in 1:N̂];
    return Ω
end

println( "Loaded prob.jl class file." )
