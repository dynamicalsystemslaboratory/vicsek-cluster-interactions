
include( "../../args.jl" )
include( "../../model.jl" )
include( "../../cluster.jl" )
include( "../../visual.jl" )

using ProgressMeter

# _, base = expfolder( "system" )
base = "/scratch/mdn7721/vicsek-dynamic-data/system/"

# Speed and density variables.
include( "../../0-params/rho-list.jl" )
include( "../../0-params/spd-list.jl" )

s = slist[parse( defInt, ARGS[1] )]
ρ = ρlist[parse( defInt, ARGS[2] )]

# Hyper parameters.
N = 2^10
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
T = 100_000

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
Tsave = 1_000_000
φsave = Vector{defFloat}( undef, Tsave + 1 )

# Initial conditions.
φsave[1] = φlist[end];

# Initialize cluster-related data variables.
f̂ = Vector{defInt}();
f = Dict{defInt,Vector{defInt}}( i => Vector{defInt}() for i ∈ 1:N )
c = zeros( defInt, N )

# Cluster information for initial conditions.
global ιp = clusterinfo( N, adjacency( N, L, r, x, y, θ ) )[2];

# Progress meter.
meter = Progress( Tsave, "Running simulation."; dt=1.0 )

# Run simulation from steady state.
for t ∈ 1:Tsave
    # Compute the current adjacency matrix.
    A = adjacency( N, L, r, x, y, θ )

    # Extract cluster information from adjacency.
    n, ιc, _ = clusterinfo( N, A )
    push!( f̂, n... )

    # Update cluster distribution count.
    for ni ∈ n;  c[ni] += 1;  end

    # Extact the unique changes between the steps and update trackers.
    i, δi = uniquechanges( ιp, ιc )
    global ιp = ιc

    # Add cluster change information to frequency list.
    countfreqlist!( N, f, i, δi )

    # Upadate the position of each agent.
    global x, y, θ = vicsek( params, A, x, y, θ )

    # Compute the order parameter.
    φsave[t+1] = order( N, θ )

    # Update progress meter.
    next!( meter )
end

# Extract cluster information for the last time step.
n, ιc, _ = clusterinfo( N, adjacency( N, L, r, x, y, θ ) )
i, δi = uniquechanges( ιp, ιc );

# Update trackers.
for ni ∈ n;  c[ni] += 1;  end
countfreqlist!( N, f, i, δi )

# Compute probability map and cluster distribution.
Δc = c./sum( c )
ω = probabilitymap( N, f );

# Extract frequency map.
bins = -1/2:N-1/2
F = hcat( [fit( Histogram, f[i] .+ (i - 1), bins ).weights for i ∈ 1:N]... )'

# Data folder name.
println( createfolder( folder ) )

# Save order parameter lists.
writedlm( folder*"order_0-$(T).txt", φlist )
writedlm( folder*"order_$(T)-end.txt", φsave )

# Save cluster information.
writedlm( folder*"cluster-frequency.txt", c )
writedlm( folder*"cluster-distribution.txt", Δc )
writedlm( folder*"frequency-map.txt", F )
