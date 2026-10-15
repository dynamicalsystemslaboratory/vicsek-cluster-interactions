
include( "../size-0-params.jl" )

i = 7
for j ∈ 1:Ns
    run( `julia --threads=24 generate-extrap.jl $(i) $(j)` )
end
