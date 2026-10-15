
include( "../../geom.jl" )
include( "../../kernel.jl" )

include( "../../0-params/rho-list.jl" )
include( "../../0-params/spd-list.jl" )

_, base = expfolder( "system" );

# System size and cutoff.
N = 128
n = round( defInt, N/2 )

# Calibration speed.
s = slist[parse( defInt, ARGS[1] )]
println( "Calibrating for speed s = ", round( s; digits=3 ) )

# Extract critical noise.
folderlist = [findfolder( N, ρlist[i], s; base=base, crit=true ) for i ∈ 1:Nρ]

# Import area and perimeter scaling values.
βlist = [vec( readdlm( folder*"area-perimeter-scaling.txt" ) )[1:2]
    for folder ∈ folderlist]

# Unpack the probability map data.
ℓlist = [polyloglaw.( 1:N; β=β ) for β ∈ βlist]
Alist = [(ℓ.^2)./(4π) for ℓ ∈ ℓlist]

# Import the cluster change frequency values.
θdata = [readdlm( folder*"optimal-rates-hat.txt" ) for folder ∈ folderlist]

# Construct the probability map.
Ωlist = [[P( j; N=N, i=i, λℓ=θlist[i,1], λr=θlist[i,2], κ=θlist[i,3] ) for i in 1:N, j in 1:N] for θlist ∈ θdata]

# Import cluster mass distribution extracted from data.
Flist = [round.( defInt, vec( readdlm( folder*"cluster-frequency.txt" ) ) )
    for folder ∈ folderlist]
flist = [F./sum( F ) for F ∈ Flist]

# Initialize progress meter.
meter = Progress( Nρ, "Calibrating for each ρ."; dt=1.0 )

# Compute the cluster mass distribution.
T = 1.0;  γ0 = [5.0]
resultslist = Vector{Any}( undef, Nρ )

# Run optimization.
for i ∈ 1:Nρ
    resultslist[i] = optimize( createkernelloss( Flist[i], N, Ωlist[i], Alist[i], ℓlist[i]; T=T, n=n ), log10.( γ0 ) )
    next!( meter )
end

# Unpack solutions.
γlist = [10.0.^results.minimizer for results ∈ resultslist]

# Save the kernel exponents and coefficients.
for i ∈ 1:Nρ
    writedlm( folderlist[i]*"kernel-parameters-cut-$(length( γlist[i] )).txt", γlist[i] )
end
