
const defInt = Int64
const defFloat = Float64

using DelimitedFiles

using Plots
using LaTeXStrings
using Measures

using StatsBase

using ProgressMeter

figurefolder = "~/Documents/papers/vicsek-dynamic/figures/";

# Define a custom color palette and figure defaults.
colorlist = [
    :cornflowerblue,  # 6699FF
    :indianred,       # CD5C5C
    :olivedrab,       # 6B8E23
    :darkgoldenrod,
    :mediumpurple,    # 9370DB
    :palevioletred,
    :lightseagreen,
    :cornflowerblue,  # 6699FF
    :indianred,       # CD5C5C
    :olivedrab,       # 6B8E23
    :darkgoldenrod,
    :mediumpurple,    # 9370DB
    :palevioletred,
    :lightseagreen,
    :cornflowerblue,  # 6699FF
    :indianred,       # CD5C5C
    :olivedrab,       # 6B8E23
    :darkgoldenrod,
    :mediumpurple,    # 9370DB
    :palevioletred,
    :lightseagreen,
    :cornflowerblue,  # 6699FF
    :indianred,       # CD5C5C
    :olivedrab,       # 6B8E23
    :darkgoldenrod,
    :mediumpurple,    # 9370DB
    :palevioletred,
    :lightseagreen
]
default(
    guidefont=font(10, "Computer Modern"),  # axis labels
    tickfont=font(10, "Computer Modern"),   # tick labels
    legendfont=font(10, "Computer Modern"), # legend
    titlefont=font(10, "Computer Modern"),  # title
    palette=colorlist,
    xformatter=:plain, dpi=100,
    left_margin=5pt, right_margin=5pt, bottom_margin=5pt, top_margin=5pt,
)

colormap1 = cgrad( :ice )
colormap2 = cgrad( [:black, :darkolivegreen, :yellowgreen, :white], [0.0, 0.20, 0.45, 0.90, 1.0] )

function expfolder(system::String; drive::String="napoli-fast/")::Tuple{String,String}
    user = ENV["USER"]
    home = "/media/"*user*"/"*drive*"vicsek-dynamic-data/"
    return home, home*system*"/"
end

function findfolder(N::defInt, ρ::defFloat, s::defFloat; base::String="", warn::Bool=true, crit::Bool=false)::String
    folder = "N-$(N)/rho-$(round( ρ, digits=3 ))/spd-$(round( s, digits=3 ))/"
    noise = crit ? "eta-$(round( readdlm( base*folder*"order-0.5.txt" )[1], digits=3 ))/" : ""
    warn && (@assert isdir( base*folder*noise ) "Folder ("*base*folder*noise*") has not been initialized.")
    return base*folder*noise
end

function findfolder(N::defInt, ρ::defFloat, s::defFloat, η::defFloat; base::String="", warn::Bool=true)::String
    folder = findfolder( N, ρ, s; base=base, warn=warn, crit=false )*"eta-$(round( η, digits=3 ))/"
    return folder
end

function findorder(folder::String, φ::defFloat)::String
    return folder*"../order-$(round( φ, digits=3 )).txt"
end

function createfolder(folder::String)::String
    if !isdir( folder )
        mkpath( folder )
    end
    return folder
end

function saveplot(plt, filename::String;
    background::Symbol=:transparent, legend_background_color::Symbol=:white, args...)
    plot!( plt; background=background, legend_background_color=legend_background_color, args... )
    savefig( plt,  filename )
    plot!( plt; background=:white )
    return plt
end

function savedict(filename::String, f::Any)
    open( filename, "w" ) do io
        for (k, v) in f
            line = join( [k; typeof( v ) == defInt ? round.( defInt, v ) : round.( v, digits=6 )], "," )
            println( io, line )
        end
    end
end

function loaddict(filename::String, dtype::DataType; M::Union{Nothing,defInt}=nothing)::Dict{defInt, Vector{dtype}}
    d = Dict{defInt, Vector{dtype}}()
    for line in eachline( filename )
        vals = parse.( Float64, split( line, "," ) )
        nums = dtype <: Integer ? round.( dtype, vals ) : vals
        k, v = nums[1], nums[2:end]
        ι = isnothing( M ) ? (1:length( v )) : sample( 1:length( v ), min( M, length( v ) ); replace=false )
        d[k] = v[ι]
    end
    return d
end

function meanbindata(bin::Tuple{defFloat,defFloat},
    xlist::Vector{defFloat}, ylist::Vector{defFloat}
    )::Tuple{defFloat,defFloat}
    Y = ylist[bin[1] .< xlist .< bin[2]]
    μ = mean( Y );  σ = std( Y )
    return (μ, σ)
end;
