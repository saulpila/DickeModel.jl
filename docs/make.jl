push!(LOAD_PATH,"../src/")
using Dicke
using Documenter,DocumenterCitations
bib = CitationBibliography(Base.Filesystem.joinpath(@__DIR__,"refs.bib"))
makedocs(bib,
         format = Documenter.HTML(
               prettyurls = get(ENV, "CI", nothing) == "true"
           ),
         sitename = "Dicke.jl",
         modules  = [Dicke],
         pages=[
                "The Dicke.jl package" => "index.md",
                "Examples" => "examples.md",
                "Documentation" => [
                    "ClassicalDicke" => "ClassicalDicke.md",
                    "ClassicalSystems" => "ClassicalSystems.md",
                    "DickeBCE" => "DickeBCE.md",
                    "UPOS" => "UPOS.md",
                    "TruncatedWignerApproximation" => "TruncatedWignerApproximation.md",
                    "FOTOCTWA" => "FOTOCTWA.md",
                    "DickeHusimiProjections" => "DickeHusimiProjections.md",
                    "ClassicalLMG" => "ClassicalLMG.md",
                    "PhaseSpaces" => "PhaseSpaces.md"
                ],
                "References" => "references.md"
               ])
deploydocs(;
    repo="github.com/saulpila/Dicke.jl",
)
