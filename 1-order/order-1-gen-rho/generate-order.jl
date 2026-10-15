
include( "../../model.jl" )
include( "../../visual.jl" )
include( "../../cluster.jl" )

using ProgressMeter
const showprogress = true

# base = "/scratch/mdn7721/vicsek-dynamic-data/system/"
_, base = expfolder( "system" )

# Hyper parameters.
N = 1024
r = 1.0

# Extract the speed argument.
include( "../../0-params/rho-list.jl" )
ρ = ρlist[parse( defInt, ARGS[1] )]
L = √(N/ρ)

# Speed from user input.
include( "../../0-params/spd-list.jl" )
s = slist[parse( defInt, ARGS[2] )]

# Create folderlist.
ηlist = [0.01;0.05:0.05:0.6]  # TODO: intermediate set (0.075:0.05:0.75)
Nη = length( ηlist )

# Print details.
println( "System size:         $(N)" )
println( "Particle speed:      $(round( s, digits=3 ))" )
println( "System-wide density: $(round( ρ, digits=3 ))" )

# Time till "steady state".
T = 10_000  # TODO: check behavior of order parameter

# Create parameter variable.
paramslist = [
    Dict{Char,Real}( 'N' => N, 'L' => L, 'r' => r, 'v' => s, 'η' => η )
    for η ∈ ηlist]

# Initialize progress meter.
showprogress && (meter = Progress( Nη*T, "Running to steady state."; dt=1.0 ))

# Initialize state and order parameter lists.
φdata = [Vector{defFloat}( undef, T + 1 ) for _ ∈ 1:Nη]
xlist = [Vector{defFloat}( undef, N ) for _ ∈ 1:Nη]
ylist = [Vector{defFloat}( undef, N ) for _ ∈ 1:Nη]
θlist = [Vector{defFloat}( undef, N ) for _ ∈ 1:Nη]

# Run simulation.
@threads for m ∈ 1:Nη
    # Unpack simulation parameters.
    params = paramslist[m]

    # Initial conditions and order parameter.
    x = L*rand( defFloat, N )
    y = L*rand( defFloat, N )
    θ = 2*π*rand( defFloat, N ) .- π
    φdata[m][1] = order( N, θ )

    for t ∈ 1:T
        A = adjacency( N, L, r, x, y, θ )
        x, y, θ = vicsek( params, A, x, y, θ )
        φdata[m][t+1] = order( N, θ )

        # Update progress meter.
        showprogress && next!( meter )
    end

    # Save state variables.
    xlist[m] = x
    ylist[m] = y
    θlist[m] = θ
end

# Initialize variables for saving.
Tsave = 10_000  # TODO: check behavior of order parameter
φsave = [Vector{defFloat}( undef, Tsave + 1 ) for _ ∈ 1:Nη]

# Initialize new progress meter.
showprogress && (meter = Progress( Nη*Tsave, "Running simulations."; dt=1.0 ))

# Set first position to final from initialization.
for m ∈ 1:Nη
    φsave[m][1] = φdata[m][end]
end

# Initialize cluster-related data variables.
clist = [zeros( defInt, N ) for _ ∈ 1:Nη];

# Run simulation from steady state.
@threads for m ∈ 1:Nη
    # Unpack parameters.
    params = paramslist[m]

    # Unpack state variables.
    x = xlist[m];  y = ylist[m];  θ = θlist[m]

    # Run simulation.
    for t ∈ 1:Tsave
        # Compute the current adjacency matrix.
        A = adjacency( N, L, r, x, y, θ )

        # Extract cluster information from adjacency.
        n, _, _ = clusterinfo( N, A )

        # Update cluster distribution count.
        for ni ∈ n;  clist[m][ni] += 1;  end

        # Upadate the position of each agent.
        x, y, θ = vicsek( params, A, x, y, θ )

        # Compute the order parameter.
        φsave[m][t+1] = order( N, θ )

        # Update progress meter.
        showprogress && next!( meter )
    end

    # Save final positions for posterity.
    xlist[m] = x;  ylist[m] = y;  θlist[m] = θ

    # Extract cluster information for the last time step.
    n, _, _ = clusterinfo( N, adjacency( N, L, r, x, y, θ ) )
    for ni ∈ n;  clist[m][ni] += 1;  end
end

# Compute probability map and cluster distribution.
Δclist = [c./sum( c ) for c ∈ clist];

# Save results of the simulation.
@threads for m ∈ 1:Nη
    # Data folder name.
    folder = createfolder( findfolder( N, ρ, s, ηlist[m]; base=base, warn=false ) )
    println( folder )

    # Save order parameter lists.
    writedlm( folder*"order_0-$(Tsave).txt", φdata[m] )
    writedlm( folder*"order_$(Tsave)-end.txt", φsave[m] )

    # Save cluster information.
    writedlm( folder*"cluster-frequency.txt", clist[m] )
    writedlm( folder*"cluster-distribution.txt", Δclist[m] )
end
