push!(LOAD_PATH,"../src/")
using Dicke
using Documenter,Bibliography
#bib = CitationBibliography(Base.Filesystem.joinpath(@__DIR__,"refs.bib"))
include("citations.jl")

bib= CitationBibliography(import_bibtex(Base.Filesystem.joinpath(@__DIR__,"refs.bib")))

#Documenter.Selectors.runner(::Type{Documenter.Expanders.ExampleBlocks}, x, page, doc)= (page.mapping[x] = Documents.RawNode(:hey, x.code))
#Documenter.Selectors.runner(::Type{Documenter.Expanders.SetupBlocks}, x, page, doc)= (page.mapping[x] = Documents.RawNode(:hey, x.code))

makedocs(bib,
         format = Documenter.HTML(
               prettyurls = get(ENV, "CI", nothing) == "true",
               assets = ["assets/customcss.css"]
           ),
         sitename = "Dicke.jl",
         modules  = [Dicke.ClassicalDicke, Dicke.ClassicalSystems,  
                    Dicke.DickeBCE,  Dicke.UPOS,  Dicke.ClassicalLMG,  Dicke.PhaseSpaces,  
                    Dicke.TruncatedWignerApproximation.Weyl,
                    Dicke.TruncatedWignerApproximation, Dicke.DickeHusimiProjections],
         pages=[
                "The Dicke.jl package" => "index.md",

                "Documentation" => [
                    "ClassicalDicke" => "ClassicalDicke.md",
                    "DickeBCE" => "DickeBCE.md",
                    "UPOS" => "UPOS.md",
                    "TruncatedWignerApproximation" => "TruncatedWignerApproximation.md",
                    "DickeHusimiProjections" => "DickeHusimiProjections.md",
                    "ClassicalLMG" => "ClassicalLMG.md",
                    "ClassicalSystems" => "ClassicalSystems.md",
                    "PhaseSpaces" => "PhaseSpaces.md"
                ],
                "Examples" => [
                    "ClassicalDicke" => "ClassicalDickeExamples.md",
                    "DickeBCE (Quantum Dicke)" => "DickeBCEExamples.md",
                    "DickeHusimiProjections" => "DickeHusimiProjectionsExamples.md",
                    "TruncatedWignerApproximation" => "TruncatedWignerApproximationExamples.md",
                    "ClassicalLMG" => "ClassicalLMGExamples.md",
                     "UPOS" => "UPOSExamples.md",
                ],
                "References" => "references.md"
               ])




deploydocs(;
    repo="github.com/saulpila/Dicke.jl",
)
