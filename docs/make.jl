push!(LOAD_PATH,"../src/")
using Dicke
using Documenter

makedocs(
         sitename = "Dicke.jl",
         modules  = [Dicke],
         pages=[
                "The Dicke.jl package" => "index.md",
                "Documentation" => [
                    "ClassicalDicke" => "classicaldicke.md"
                ]

               ])
deploydocs(;
    repo="github.com/saulpila/Dicke.jl",
)
