
include( "../size-0-params.jl" )

Nk = 5
for k ∈ 1:Nk, i ∈ 1:Nρ
    run( `julia --threads=24 generate-geometry.jl $(i) $(k)` )
end
