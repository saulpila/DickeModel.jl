push!(LOAD_PATH,"../src/")
using Base.Filesystem
using Revise
using Dicke
using Documenter
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
try
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
catch e
        showerror(stdout, e)
end
end
