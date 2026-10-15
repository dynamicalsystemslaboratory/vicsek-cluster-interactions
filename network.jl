
include( "args.jl" )

using SparseArrays

function periodicdistance(L::defFloat, r::defFloat, dx::defFloat, dy::defFloat)
    # Compute raw Euclidian distance.
    dx = abs( dx )
    dy = abs( dy )

    # Adjust for periodic boundary.
    dx = dx > L/2 ? L - dx : dx
    dy = dy > L/2 ? L - dy : dy

    # Return distance values.
    return (dx^2 + dy^2) < r^2
end

function createbuckets(N::defInt, L::defFloat, r::defFloat,
    xlist::Vector{defFloat}, ylist::Vector{defFloat})

    # Grid index for each agent.
    n = max( 1, floor( defInt, L/r ) )
    n0 = L/n

    # Create dictionary of bucket-index combinations.
    bdict = Dict{Tuple{defInt,defInt}, Vector{defInt}}()
    for k in 1:N
        # Extract bucket indeces.
        i = mod1( floor( defInt, xlist[k]/n0 ) + 1, n )
        j = mod1( floor( defInt, ylist[k]/n0 ) + 1, n )

        # For each agent, save its index to the appropriate grid bucket.
        push!( get!( () -> defInt[], bdict, (i, j) ), k )
    end

    # Return the bin information and the compiled dictionary.
    return n, bdict
end

function adjacency(N::defInt, L::defFloat, r::defFloat,
    xlist::Vector{defFloat}, ylist::Vector{defFloat}, θlist::Vector{defFloat})
    # Partition agent list into grid-based buckets.
    n, bdict = createbuckets( N, L, r, xlist, ylist )

    # Initialize list of column/row indeces for connection tracking.
    rlist, clist = Vector{defInt}(), Vector{defInt}()

    # For each unique grid ID contained in the buckets dictionary.
    for (key, ilist) in bdict
        # Unpack current grid bucket.
        i, j = key

        # For each grid bucket adjacent to the current bucket.
        for di in -1:1, dj in -1:1
            # Get current grid bucket to check.
            nkey = (mod1( i + di, n ), mod1( j + dj, n ))

            # If grid box is empty, skip loop iteration.
            haskey( bdict, nkey ) || continue

            # Unpack neighbors in adjacent grid bucket.
            jlist = bdict[nkey]
            for i in ilist, j in jlist
                # Form self-connections.
                i == j && (push!( rlist, i ); push!( clist, j ); continue)

                # Avoid double-computations.
                i < j || continue

                # Compute raw distance.
                dx = xlist[i] - xlist[j]
                dy = ylist[i] - ylist[j]

                # Skip pairs that are not within distance.
                periodicdistance( L, r, dx, dy ) || continue

                # Update index of connections list.
                push!( rlist, i );  push!( clist, j )
                push!( rlist, j );  push!( clist, i )
            end
        end
    end

    return sparse( rlist, clist, ones( defInt, length( rlist ) ), N, N )
end

println( "Loaded network.jl class file." )
