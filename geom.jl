``
include( "args.jl" )
include( "model.jl" )

using StatsBase
using LazySets: VPolygon, area
using LinearAlgebra
using Distributions

function perimeter(P::VPolygon{defFloat,Vector{defFloat}})::defFloat
    # Unpack vertices.
    vlist = P.vertices

    # Extract the number of vertices in the polygon.
    Nv = length( vlist )
    Nv < 2 && return 0.0

    # Compute perimeter length.
    ℓ = 0.0
    for i ∈ 1:Nv
        # Unpack focal vertices.
        v1, v2 = vlist[i], vlist[mod1( i+1, Nv )]

        # Add constribution to perimeter.
        ℓ += hypot( v2[1] - v1[1], v2[2] - v1[2] )
    end
    return ℓ
end

function countgeomlist!(L::defFloat, r::defFloat,
    Â::Dict{defInt,Vector{defFloat}}, ℓ::Dict{defInt,Vector{defFloat}},
    xlist::Vector{defFloat}, ylist::Vector{defFloat},
    A::SparseMatrixCSC{defInt, defInt})::Nothing
    # Extract clusters from adjacency.
    clist = connected_components( Graph( A ) )

    # For each cluster, compute the convex hull area.
    for c ∈ clist
        # Extract size of cluster and its shape.
        i = length( c )
        i > 1 || continue

        # Extract points from the cluster.
        xc = xlist[c];  yc = ylist[c]

        # Skip if the cluster is within radial proximity of the lower and left boundaries.
        #   to avoid clusters broken in Euclidian space by the periodic boundary.
        any( (xc.^2 .< r^2) .| (yc.^2 .< r^2) ) && continue

        # Extract agent polygon and compute area.
        P = VPolygon( permutedims( hcat( xc, yc ) ) )
        push!( Â[i], area( P ) )
        push!( ℓ[i], perimeter( P ) )
    end
end

function fitpower(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})::Vector{defFloat}
    x̂ = log10.( x );  ŷ = log10.( y )
    return [x̂ ones( defFloat, length( x̂ ) )]\ŷ
end

function powerlaw(i::defInt; λ::Vector{defFloat}=[1.0,1.0])::defFloat
    return 10.0^λ[2]*i^λ[1]
end

function fitlogpower(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; cut::Real=2)::Vector{defFloat}
    i = x .> cut
    x̂1 = log10.( x[i] .+ 1 )
    x̂2 = log10.( log10.( x[i] .+ 1 ) )
    ŷ  = log10.( y[i] )
    return [x̂1 x̂2 ones( defFloat, length( x̂1 ) )] \ ŷ
end

function logpowerlaw(i::defInt; γ::Vector{defFloat}=[1.0,1.0,1.0])::defFloat
    return 10.0^γ[3]*(i + 1)^γ[1]*log10( i + 1 )^γ[2]
end

function fitpolylog(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
    cut::Real=2, verbose::Bool=false, CI::defFloat=0.95)
    # Data transformation.
    i = x .> cut
    x̂ = log10.( log10.( x[i] .+ 1 ) )
    ŷ  = log10.( y[i] )

    # Build design matrix and return.
    X = hcat( x̂, ones( defFloat, length( x̂ ) ) )
    θ = X\ŷ
    !verbose && (return θ)

    # Compute the noise variance.
    residuals = ŷ - X*θ
    n, p = size( X )
    σ = sum( residuals.^2 )/(n - p)

    # Compute the covariance.
    Σ = σ*inv( X'*X )
    se = .√diag( Σ )

    # Return the confidence interval width.
    z = quantile( TDist( n - p ), (1 + CI)/2 )
    return θ, Σ, z*se
end

function polyloglaw(i::defInt; β::Vector{defFloat}=[1.0,1.0])::defFloat
    return 10.0^β[2]*log10( i + 1 )^β[1]
end

println( "Loaded geom.jl class files." )
