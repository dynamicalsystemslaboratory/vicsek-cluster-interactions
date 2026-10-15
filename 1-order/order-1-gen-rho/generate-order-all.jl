
using Base.Threads

include( "../../args.jl" )

include( "../../0-params/rho-list.jl" )
include( "../../0-params/spd-list.jl" )

# j = parse( defInt, ARGS[1] )
for i ∈ 1:Nρ, j ∈ 1:Ns
    run( `julia --threads=24 generate-order.jl $(i) $(j)` )
end
