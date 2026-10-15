
include( "args.jl" )

using StatsBase
using SparseArrays
using Graphs

# Extract cluster information from an adjacency matrix.
function clusterinfo(N::defInt, A::SparseMatrixCSC{defInt, defInt}
    )::Tuple{Vector{defInt}, Vector{defInt}, Vector{defInt}}
    # Turn into graph variable and extract clusters.
    G = Graph( A )
    C = connected_components( G )

    # Extract cluster sizes.
    mlist = defInt.( [length( c ) for c ∈ C] )              # Size of each connected component of A.

    # Extract the size of the cluster each agent is a member of.
    ιlist = [findfirst( c -> i ∈ c, C ) for i ∈ 1:N]        # Cluster index of each node.
    ilist = [mlist[ι] for ι ∈ ιlist]                        # Cluster size of each node.

    return mlist, ιlist, ilist
end

# Extract steps where at least one cluster changes in size.
function extractsteps(N::defInt, T::defInt, ιdata::Matrix{defInt})::Vector{defInt}
    # Extract the steps from the simulation where there was a change in clusters.
    tlist = Vector{defInt}()
    for t ∈ 1:T-1
        sum( ιdata[t,:] .== ιdata[t+1,:] ) == N || push!( tlist, t )
    end

    # Return list of time steps.
    return tlist
end

# Extract the unique changes in cluster size and the associated magnitude.
function uniquechanges(ι1::Vector{defInt}, ι2::Vector{defInt}; zero::Bool=false)::Union{Tuple{Vector{defInt}, Vector{defInt}},Vector{defInt}}
    ιseen  = Set{defInt}()
    ilist  = Vector{defInt}()   # List of cluster sizes.
    δilist = Vector{defInt}()   # List of changes in the cluster size.
    ιlist  = Vector{defInt}()   # List of identities for each cluster.

    for ι in 1:maximum( ι1 )
        ι ∉ ιseen || continue

        # Identify members of cluster ι at time t.
        p1 = ι1 .== ι
        i1 = sum( p1 )

        # Get representative agent in cluster ι.
        p = findfirst( p1 )

        # Identify size of that agent's new cluster at t+1.
        p2 = ι2 .== ι2[p]
        i2 = sum( p2 )
        i2 == i1 && (push!( ιlist, ι ); continue)

        # Record initial size and change.
        !zero && (push!( ilist, i1 ); push!( δilist, i2 - i1 ))

        # If the cluster grew, mark the other member of the collision as visited.
        i1 < i2 && push!( ιseen, unique( ι1[(.!p1) .& p2] )... )
    end

    return zero ? ιlist : (ilist, δilist)
end

# Count the frequency of events under certain combination conditions.
#   1. i joins with δi
#       1.1. f(i|+δi) += 1
#       1.2. f(δi|+i) += 1
#   2. i breaks into -δi and i+δi
#       2.1. f(i|δi) += 1
#       2.2. f(i|-(i+δi)) += 1
# Count frequency list.
function countfreq!(N, f, i, δi; include0=false)
    # Check current size and change, return error if invalid.
    @assert (1 ≤ i+δi ≤ N) "ERROR: The change in cluster size is invalid, got $(i + δi)."

    # If change is zero, return original P.
    include0 || δi ≠ 0 || return f

    # Get secondary cluster.
    j = δi ≥ 0 ? δi : i

    # Push instances onto frequency dictionary.
    push!( f[i], δi ≥ 0 ? j : δi )
    push!( f[j], δi ≥ 0 ? i : -(i+δi) )

    # Return updated frequency lists.
    return f
end

# Append the list of (i,δi)-pairs to the frequency list f.
function countfreqlist!(N, f, ilist, δilist; include0=false)
    for (i, δi) ∈ zip( ilist, δilist )
        countfreq!( N, f, i, δi; include0=include0 )
    end
end

# Construct probability map from cluster transition data.
function probabilitymap(N, f)::Matrix{defFloat}
    # Create probability map.
    Ω = zeros( defFloat, (N, N) )
    for i ∈ 1:N
        # Compute change in mass probabilities.
        jfreq = countmap( f[i] )
        jkeys = collect( keys( jfreq ) )

        # Pack frequencies into map.
        for j ∈ jkeys[jkeys .!= 0]
            Ω[i,i+j] = jfreq[j]
        end

        # Divide by the number of occurrences.
        Ω[i,:] /= sum( Ω[i,:] )
    end

    # Return the probability map.
    return Ω
end

println( "Loaded cluster.jl class file." )
