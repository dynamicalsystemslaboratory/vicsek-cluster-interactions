
include( "../../model.jl" )
include( "../../visual.jl" )
include( "../../cluster.jl" )

using ProgressMeter
showmeter = true

# base = "/scratch/mdn7721/vicsek-dynamic-data/system/"
_, base = expfolder( "system" )

# Speed and density variables.
include( "../size-0-params.jl" )
ρ = ρlist[parse( defInt, ARGS[1] )]
s = slist[parse( defInt, ARGS[2] )]

# System size.
N0 = 2
N = round( defInt, 2.0^(parse( defInt, ARGS[3] ) + N0) )

# Environment parameters.
L = √(N/ρ)
r = 1.0

# Create folderlist.
# ηlist = [0.01;0.025:0.025:0.6]
ηlist = 0.425:0.025:0.6
Nη = length( ηlist )

# Print details.
println( "System size:         $(N)" )
println( "Particle speed:      $(round( s, digits=3 ))" )
println( "System-wide density: $(round( ρ, digits=3 ))" )

# Time till "steady state".
T = 100_000

# Create parameter variable.
paramslist = [
    Dict{Char,Real}( 'N' => N, 'L' => L, 'r' => r, 'v' => s, 'η' => η )
    for η ∈ ηlist]

# Initialize progress meter.
showmeter && (meter = Progress( Nη*T, "Running to steady state."; dt=0.1 ))

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
        showmeter && next!( meter )
    end

    # Save state variables.
    xlist[m] = x
    ylist[m] = y
    θlist[m] = θ
end

# Initialize variables for saving.
Tsave = 100_000
φsave = [Vector{defFloat}( undef, Tsave + 1 ) for _ ∈ 1:Nη]

# Set first position to final from initialization.
for m ∈ 1:Nη
    φsave[m][1] = φdata[m][end]
end

# Initialize cluster-related data variables.
clist = [zeros( defInt, N ) for _ ∈ 1:Nη]

# Initialize cluster-related data variables.
τlist = [Dict{defInt,Vector{defInt}}( i => Vector{defInt}() for i ∈ 1:N ) for _ ∈ 1:Nη];

# Progress meter.
showmeter && (meter = Progress( Nη*Tsave, "Running simulation."; dt=0.1 ))

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
        showmeter && next!( meter )
    end

    # Save final positions for posterity.
    xlist[m] = x;  ylist[m] = y;  θlist[m] = θ

    # Extract cluster information for the last time step.
    n, ιc, ic = clusterinfo( N, adjacency( N, L, r, x, y, θ ) )
    for ni ∈ n;  clist[m][ni] += 1;  end
    ι = uniquechanges( ιp, ιc; zero=true );

    # Update trackers.
    countlifelist!( N, τlist[m], ip, ι, ιp, τc )
end

# Compute probability map and cluster distribution.
Δclist = [c./sum( c ) for c ∈ clist];

@threads for m ∈ 1:Nη
    # Data folder name.
    folder = createfolder( findfolder( N, ρ, s, ηlist[m]; base=base, warn=false ) )
    println( folder )

    # Save cluster information.
    writedlm( folder*"order_0-$(T).txt", φdata[m] )
    writedlm( folder*"order_$(T)-end.txt", φsave[m] )

    # Save cluster information.
    writedlm( folder*"cluster-frequency.txt", clist[m] )
    writedlm( folder*"cluster-distribution.txt", Δclist[m] )

    # Save persistence time.
    savedict( folder*"cluster-life-times.txt", τlist[m] )
end
