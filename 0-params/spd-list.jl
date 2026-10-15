
δs = 0.1;  smin = -2.0;  smax = 0.0
slist = 10.0.^(smin:δs:smax)
Ns = length( slist )
println( "Import $(Ns) speeds: $(round.( 10.0.^[smin,smax], digits=3 )) (log spacing)." )
