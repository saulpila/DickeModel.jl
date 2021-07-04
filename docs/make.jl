push!(LOAD_PATH,"../src/")
using Dicke
using Documenter,Bibliography
#bib = CitationBibliography(Base.Filesystem.joinpath(@__DIR__,"refs.bib"))
include("citations.jl")

bib= CitationBibliography(import_bibtex(Base.Filesystem.joinpath(@__DIR__,"refs.bib")))

makedocs(bib,
         format = Documenter.HTML(
               prettyurls = get(ENV, "CI", nothing) == "true",
               assets = ["assets/customcss.css"]
           ),
         sitename = "Dicke.jl",
         modules  = [Dicke],
         pages=[
                "The Dicke.jl package" => "index.md",
                "Examples" => [
                    "Classical Dicke" => "ClassicalDickeExamples.md",
                    "Truncated Wigner Approximation" => "TruncatedWignerApproximationExamples.md",
                ],
                "Documentation" => [
                    "ClassicalDicke" => "ClassicalDicke.md",
                    "ClassicalSystems" => "ClassicalSystems.md",
                    "DickeBCE" => "DickeBCE.md",
                    "UPOS" => "UPOS.md",
                    "TruncatedWignerApproximation" => "TruncatedWignerApproximation.md",
                    "DickeHusimiProjections" => "DickeHusimiProjections.md",
                    "ClassicalLMG" => "ClassicalLMG.md",
                    "PhaseSpaces" => "PhaseSpaces.md"
                ],
                "References" => "references.md"
               ])
deploydocs(;
    repo="github.com/saulpila/Dicke.jl",
)
