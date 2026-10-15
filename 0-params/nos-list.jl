
# Variable effective temperature.
δη = 0.05;  ηlist = [0.01;δη:δη:0.6]
ηmin = minimum( ηlist );  ηmax = maximum( ηlist )
Nη = length( ηlist )
println( "Import $(length( ηlist )) noise: $(round.( [ηmin,ηmax], digits=3 )) (linear spacing)." )
