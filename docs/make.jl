push!(LOAD_PATH,"../src/")
using Dicke
using Documenter,DocumenterCitations
bib = CitationBibliography(Base.Filesystem.joinpath(@__DIR__,"refs.bib"))
makedocs(bib,
         sitename = "Dicke.jl",
         modules  = [Dicke],
         pages=[
                "The Dicke.jl package" => "index.md",
				"Examples" => "examples.md",
                "Documentation" => [
                    "ClassicalDicke" => "classicaldicke.md",
					"ClassicalSystems" => "classicalsystems.md"
                ],
				"References" => "references.md"
               ])
deploydocs(;
    repo="github.com/saulpila/Dicke.jl",
)
