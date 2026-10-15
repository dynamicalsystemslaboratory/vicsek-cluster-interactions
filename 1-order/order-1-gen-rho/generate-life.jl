
include( "../../args.jl" )
include( "../../model.jl" )
include( "../../visual.jl" )
include( "../../cluster.jl" )

using ProgressMeter

# base = "/scratch/mdn7721/vicsek-dynamic-data/system/"
_, base = expfolder( "system" )

# Hyper parameters.
N = 128
s = 0.251
r = 1.0

# Extract the speed argument.
include( "../../0-params/rho-list.jl" )
ρ = ρlist[parse( defInt, ARGS[1] )]
L = √(N/ρ)

# Create folderlist.
ηlist = [0.01;0.025:0.025:0.4]
Nη = length( ηlist )

# Print details.
println( "Particle speed:      $(round( s, digits=3 ))" )
println( "System-wide density: $(round( ρ, digits=3 ))" )

# Time till "steady state".
T = 100_000

# Create parameter variable.
paramslist = [
    Dict{Char,Real}( 'N' => N, 'L' => L, 'r' => r, 'v' => s, 'η' => η )
    for η ∈ ηlist]

# Initialize progress meter.
P = Progress( Nη*T, "Running to steady state."; dt=1.0 )

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
    x = L*rand( N, )
    y = L*rand( N, )
    θ = 2*π*rand( N, ) .- π
    φdata[m][1] = order( N, θ )

    for t ∈ 1:T
        A = adjacency( N, L, r, x, y, θ )
        x, y, θ = vicsek( params, A, x, y, θ )
        φdata[m][t+1] = order( N, θ )

        # Update progress meter.
        next!( P )
    end

    # Save state variables.
    xlist[m] = x
    ylist[m] = y
    θlist[m] = θ
end

# Initialize variables for saving.
Tsave = 1_000_000
φsave = [Vector{defFloat}( undef, Tsave + 1 ) for _ ∈ 1:Nη]

# Initialize new progress meter.
P = Progress( Nη*Tsave, "Running simulations."; dt=1.0 )

# Set first position to final from initialization.
for m ∈ 1:Nη
    φsave[m][1] = φdata[m][end]
end

# Initialize cluster-related data variables.
τlist = [Dict{defInt,Vector{defInt}}( i => Vector{defInt}() for i ∈ 1:N ) for _ ∈ 1:Nη];

# Run simulation from steady state.
@threads for m ∈ 1:Nη
    # Unpack parameters.
    params = paramslist[m]

    # Unpack state variables.
    x = xlist[m];  y = ylist[m];  θ = θlist[m]

    # Initialize step-specific trackers.
    n, ιp, ip = clusterinfo( N, adjacency( N, L, r, x, y, θ ) )
    τc = zeros( defInt, N )

    # Run simulation.
    for t ∈ 1:Tsave
        # Compute the current adjacency matrix.
        A = adjacency( N, L, r, x, y, θ )

        # Extract cluster information from adjacency.
        n, ιc, ic = clusterinfo( N, A )

        # Update cluster distribution count.
        for ni ∈ n;  clist[m][ni] += 1;  end

        # Extact the unique changes between the steps and update trackers.
        ι = uniquechanges( ιp, ιc; zero=true )

        # Add cluster life length information to frequency list.
        countlifelist!( N, τlist[m], ip, ι, ιp, τc )
        ip = ic
        ιp = ιc

        # Update the position of each agent.
        x, y, θ = vicsek( params, A, x, y, θ )

        # Compute the order parameter.
        φsave[m][t+1] = order( N, θ )

        # Update progress meter.
        next!( P )
    end

    # Save final positions for posterity.
    xlist[m] = x;  ylist[m] = y;  θlist[m] = θ

    # Extract cluster information for the last time step.
    n, ιc, ic = clusterinfo( N, adjacency( N, L, r, x, y, θ ) )
    ι = uniquechanges( ιp, ιc; zero=true );

    # Update trackers.
    countlifelist!( N, τlist[m], ip, ι, ιp, τc )
end

# Save results of the simulation.
@threads for m ∈ 1:Nη
    # Data folder name.
    folder = createfolder( findfolder( N, ρ, s, ηlist[m]; base=base, warn=false ) )
    println( folder )

    # Save cluster information.
    writedlm( folder*"order_$(T)-end.txt", φsave[m] )
    savedict( folder*"cluster-life-times.txt", τlist[m] )
end
