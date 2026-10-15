
include( "../size-0-params.jl" )

Nk = 5
for k ∈ 1:Nk
    for i ∈ 1:Nρ, j ∈ 1:Ns
        run( `julia --threads=24 generate-cluster.jl $(i) $(j) $(k)` )
    end
end
