
include( "../../geom.jl" )
include( "../../kernel.jl" )

_, base = expfolder( "system" )

# Import parameters.
include( "../size-0-params.jl" )

# Unpack density.
ρ = ρlist[parse( defInt, ARGS[1] )]
println( "Density: ", ρ )

# Unpack speed.
s = slist[parse( defInt, ARGS[2] )]
println( "Speed:   ", s )

# System parameters.
N = 128;  N̂ = 1024

# Extract critical noise.
folder = findfolder( N, ρ, s; base=base, crit=true )
f̂older = findfolder( N̂, ρ, s; base=base, crit=true )

# Unpack the kernel parameters.
f0 = 100.0
γ = first( readdlm( folder*"kernel-parameters-1.txt" ) )

# Import area and perimeter scaling values.
β = first( readdlm( folder*"area-perimeter-scaling.txt" ), 2 )

# Import the cluster change frequency values.
θ = first( transpose( readdlm( folder*"optimal-rates-hat.txt" ) ), 3 )

# Extract the power law.
λ = first( readdlm( folder*"power-law.txt" ) )

# Unpack the probability map data.
ℓ = polyloglaw.( 1:N̂; β=β )
A = (ℓ.^2)./(4π)

# Construct the probability map.
Ω = [P( k2; N=N̂, i=k1, λℓ=θ[1], λr=θ[2], κ=θ[3] ) for k1 in 1:N̂, k2 in 1:N̂]

# Import cluster mass distribution extracted from data.
F = vec( readdlm( f̂older*"cluster-frequency.txt" ) )
f = F./sum( F )

# Import kernel coefficient parameters.
α = vec( readdlm( folder*"kernel-decay-strength.txt" ) )

# Adjust the coagulation-fragmentation ratio for the kernel strength.
γ̂ = f0*(2.0^first( α )*N̂.^last( α ))

# Initialize parameter variables.
params = Params( N̂, Ω, createkernels( N̂, A, ℓ, [γ̂] )... )

# Simulate the master equations in the larger system size.
T = 2.0
f̂ = logsimulate( params; T=T, f=zipf.( 1:N̂; b=N̂, λ=λ ) )

# Save the extrapolated distribution.
writedlm( f̂older*"cluster-distribution-hat.txt", f̂ )
println( "Simulation finished and cluster mass distribution saved." )
