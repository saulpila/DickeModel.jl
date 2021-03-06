push!(LOAD_PATH,"../src/")
using Base.Filesystem
using Revise
using DickeModel
using Documenter,DocumenterCitations
while true
print("Press enter to compile docs")
readline()
try
	rm("build", recursive=true,force=true)
catch
	print("folder locked, close tabs")
	continue
end
revise(DickeModel)
sleep(3)
bib = CitationBibliography("refs.bib")

try
makedocs(bib,
         sitename = "DickeModel.jl",
         modules  = [DickeModel],
         pages=[
                "The DickeModel.jl package" => "index.md",
                "Documentation" => [
                    "ClassicalDicke" => "classicaldicke.md",
					"ClassicalSystems" => "classicalsystems.md"
                ],
				"References" => "references.md"

               ])
catch e
        showerror(stdout, e)
end
end
