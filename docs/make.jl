push!(LOAD_PATH,"../src/")
using DickeModel
using Documenter,Bibliography
using Distributed
addprocs(2,exeflags=`--project=docs/`)
@everywhere push!(LOAD_PATH,"../src/")
@everywhere using DickeModel
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
         sitename = "DickeModel.jl",
         modules  = [DickeModel.ClassicalDicke, DickeModel.ClassicalSystems,  
                    DickeModel.DickeBCE,  DickeModel.UPOs,  DickeModel.ClassicalLMG,  DickeModel.PhaseSpaces,  
                    DickeModel.TWA.Weyl,
                    DickeModel.TWA, DickeModel.EnergyShellProjections],
         pages=[
                "The DickeModel.jl package" => "index.md",

                "Documentation" => [
                    "ClassicalDicke" => "ClassicalDicke.md",
                    "DickeBCE" => "DickeBCE.md",
                    "UPOs" => "UPOs.md",
                    "TWA" => "TWA.md",
                    "EnergyShellProjections" => "EnergyShellProjections.md",
                    "ClassicalLMG" => "ClassicalLMG.md",
                    "ClassicalSystems" => "ClassicalSystems.md",
                    "PhaseSpaces" => "PhaseSpaces.md"
                ],
                "Examples" => [
                    "ClassicalDicke" => "ClassicalDickeExamples.md",
                    "DickeBCE (Quantum Dicke)" => "DickeBCEExamples.md",
                    "EnergyShellProjections" => "EnergyShellProjectionsExamples.md",
                    "TWA" => "TWAExamples.md",
                    "ClassicalLMG" => "ClassicalLMGExamples.md",
                     "UPOs" => "UPOsExamples.md",
                ],
                "References" => "references.md"
               ])




deploydocs(;
    repo="github.com/saulpila/DickeModel.jl",
)
