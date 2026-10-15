
include( "../../geom.jl" )
include( "../../cluster.jl" )
include( "../../visual.jl" )

using ProgressMeter

# base = "/scratch/mdn7721/vicsek-dynamic-data/system/"
_, base = expfolder( "system" )

# Speed and density variables.
include( "../../0-params/rho-list.jl" )
include( "../../0-params/spd-list.jl" )

ρ = ρlist[parse( defInt, ARGS[1] )]
s = slist[parse( defInt, ARGS[2] )]

# Hyper parameters.
N = 1024
L = √(N/ρ)
r = 1.0
println( "System size: $(N)." )

# Create folder and extract critical noise.
folder = createfolder( findfolder( N, ρ, s; base=base, warn=false, crit=true ) )
η = readdlm( folder*"../order-0.5.txt" )[1]

# Create parameter variable.
params = Dict{Char,Real}( 'N' => N, 'L' => L, 'r' => r, 'v' => s, 'η' => η )

# Print details.
println( "Particle speed:      $(round( s, digits=3 ))." )
println( "System-wide density: $(round( ρ, digits=3 ))." )
println( "Critical noise:      $(round( η, digits=3 ))." )

# Time till "steady state".
T = 10_000

# Initialize progress meter.
meter = Progress( T, "Running to steady state."; dt=1.0 )

# Initial conditions.
global x = L*rand( N, )
global y = L*rand( N, )
global θ = 2*π*rand( N, ) .- π

# Run simulation.
φlist = Vector{defFloat}( undef, T + 1 )
φlist[1] = order( N, θ )
for t ∈ 1:T
    A = adjacency( N, L, r, x, y, θ )
    global x, y, θ = vicsek( params, A, x, y, θ )
    φlist[t+1] = order( N, θ )
    next!( meter )
end

# Initialize steady state variables for saving.
Tsave = 100_000
φsave = Vector{defFloat}( undef, Tsave + 1 )

# Initial conditions.
φsave[1] = φlist[end];

# Initialize cluster information.
Â = Dict{defInt,Vector{defFloat}}( [i => Vector{defFloat}() for i ∈ 2:N] )
ℓ = Dict{defInt,Vector{defFloat}}( [i => Vector{defFloat}() for i ∈ 2:N] );

# Cluster information for initial conditions.
global ιp = clusterinfo( N, adjacency( N, L, r, x, y, θ ) )[2];

# Progress meter.
meter = Progress( Tsave, "Running simulation."; dt=1.0 )

# Run simulation from steady state.
for t ∈ 1:Tsave
    # Compute the current adjacency matrix.
    A = adjacency( N, L, r, x, y, θ )

    # Extract cluster information from adjacency.
    if (t - 1)%4 == 0
        countgeomlist!( L, r, Â, ℓ, x, y, A )
    end

    # Upadate the position of each agent.
    global x, y, θ = vicsek( params, A, x, y, θ )

    # Update progress meter.
    next!( meter )
end

# Save results of the simulation.
println( createfolder( folder ) )

# Save area distribution information.
savedict( folder*"convex-hull-area.txt", Â )
savedict( folder*"convex-hull-perimeter.txt", ℓ )
