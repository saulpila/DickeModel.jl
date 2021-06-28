push!(LOAD_PATH,"../src/")
using Base.Filesystem
using Revise
using Dicke
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
revise(Dicke)
sleep(3)
bib = CitationBibliography("refs.bib")

try
makedocs(bib,
         sitename = "Dicke.jl",
         modules  = [Dicke],
         pages=[
                "The Dicke.jl package" => "index.md",
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
