
include( "../../geom.jl" )
include( "../../cluster.jl" )
include( "../../visual.jl" )

using ProgressMeter

# base = "/scratch/mdn7721/vicsek-dynamic-data/system/"
_, base = expfolder( "system" )

# Speed and density variables.
include( "../size-0-params.jl" )
ρ = ρlist[parse( defInt, ARGS[1] )]

# System size.
N0 = 2
N = round( defInt, 2.0^(parse( defInt, ARGS[2] ) + N0) )

# Environment parameters.
L = √(N/ρ)
r = 1.0

# Create folder list and extract critical noise.
folderlist = [findfolder( N, ρ, s; base=base, warn=false, crit=true ) for s ∈ slist]
ηlist = first.( readdlm.( folderlist.*"../order-0.5.txt" ) )

# Print details.
println( "System size:         $(N)" )
println( "System-wide density: $(round( ρ, digits=3 ))" )
println( "Critical noise: $(ηlist)" )

# Time till "steady state".
T = 100_000

# Create parameter variable.
paramslist = [
    Dict{Char,Real}( 'N' => N, 'L' => L, 'r' => r, 'v' => slist[m], 'η' => ηlist[m] )
    for m ∈ 1:Ns]

# Initialize progress meter.
meter = Progress( Ns*T, "Running to steady state."; dt=1.0 )

# Initialize state and order parameter lists.
φdata = [Vector{defFloat}( undef, T + 1 ) for _ ∈ 1:Ns]
xlist = [Vector{defFloat}( undef, N ) for _ ∈ 1:Ns]
ylist = [Vector{defFloat}( undef, N ) for _ ∈ 1:Ns]
θlist = [Vector{defFloat}( undef, N ) for _ ∈ 1:Ns]

# Run simulation.
@threads for m ∈ 1:Ns
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
        next!( meter )
    end

    # Save state variables.
    xlist[m] = x
    ylist[m] = y
    θlist[m] = θ
end

# Initialize variables for saving.
Tsave = 100_000
φsave = [Vector{defFloat}( undef, Tsave + 1 ) for _ ∈ 1:Ns]

# Set first position to final from initialization.
for m ∈ 1:Ns
    φsave[m][1] = φdata[m][end]
end

# Initialize cluster information.
Alist = [Dict{defInt,Vector{defFloat}}( [i => Vector{defFloat}() for i ∈ 2:N] ) for _ ∈ 1:Ns]
ℓlist = [Dict{defInt,Vector{defFloat}}( [i => Vector{defFloat}() for i ∈ 2:N] ) for _ ∈ 1:Ns];

# Progress meter.
meter = Progress( Ns*Tsave, "Running simulation." )

# Run simulation from steady state.
@threads for m ∈ 1:Ns
    # Unpack parameters.
    params = paramslist[m]

    # Unpack state variables.
    x = xlist[m];  y = ylist[m];  θ = θlist[m]

    # Run simulation.
    for t ∈ 1:Tsave
        # Compute the current adjacency matrix.
        A = adjacency( N, L, r, x, y, θ )

        # Extract cluster information from adjacency.
        countgeomlist!( L, r, Alist[m], ℓlist[m], x, y, A )

        # Update the position of each agent.
        x, y, θ = vicsek( params, A, x, y, θ )

        # Compute the order parameter.
        φsave[m][t+1] = order( N, θ )

        # Update progress meter.
        next!( meter )
    end

    # Save final positions for posterity.
    xlist[m] = x;  ylist[m] = y;  θlist[m] = θ
end

@threads for m ∈ 1:Ns
    # Data folder name.
    folder = folderlist[m]
    println( folder )

    # Save area distribution information.
    savedict( folder*"convex-hull-area.txt", Alist[m] )
    savedict( folder*"convex-hull-perimeter.txt", ℓlist[m] )
end
