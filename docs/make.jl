push!(LOAD_PATH,"../src/")
using Dicke
using Documenter
makedocs(
         sitename = "Dicke.jl",
         modules  = [Dicke],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/USERNAME/VegaGraphs.jl",
)
