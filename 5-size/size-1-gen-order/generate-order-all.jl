
include( "../size-0-params.jl" )

NN = 5

for k ∈ NN:NN
    for i ∈ Nρ-1:Nρ, j ∈ 1:Ns
        run( `julia --threads=24 generate-order.jl $(i) $(j) $(k)` )
    end
end
