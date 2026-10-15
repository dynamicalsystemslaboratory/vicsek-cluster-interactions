
δρ = 0.1;  ρmin = -2.0;  ρmax = 0.0
ρlist = 10.0.^(ρmin:δρ:ρmax)
Nρ = length( ρlist )
println( "Import $(Nρ) density: $(round.( 10.0.^[ρmin,ρmax], digits=3 )) (log spacing)." )
