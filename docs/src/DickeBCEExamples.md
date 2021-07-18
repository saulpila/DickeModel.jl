# Examples for DickeBCE 

```@setup examples
push!(LOAD_PATH,"../../src")
on_github=get(ENV, "CI", nothing) == "true"
cache_fold_name="./diags"
use_current_dir_for_diags=on_github
using DickeModel
```
The module [`DickeModel.DickeBCE`](@ref DickeModel.DickeBCE) works with the quantum Dicke model
using a very efficient basis known as the coherent efficient basis (BCE for its acronym in Spanish).
See Refs. [Bastarrachea2014PSa](@cite) and [Bastarrachea2014PSb](@cite) for a detailed explanation on how and why it works. 
Throughout this examples, we will work with a system size of `j = 30`, but  using this module you can easily go up
to `j = 100` (as done in Refs. [Pilatowsky2021](@cite), [Pilatowsky2021NatCommun](@cite), [Pilatowsky2021Identification](@cite), [Villasenor2021](@cite)) and beyond.

## Diagonalizing the Dicke Hamiltonian
Let us start by defining our parameters:
```@setup examples
@info "Starting example: Diagonalization"
```
```@example examples
using DickeModel.DickeBCE, DickeModel.ClassicalDicke
system = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j=30, Nmax=120)
nothing; #hide
```

To load the eigenbasis, simply use [`diagonalization`](@ref DickeModel.DickeBCE.diagonalization):

```@example examples
if !use_current_dir_for_diags #hide
@time eigenenergies,eigenstates = diagonalization(system)
else #hide
@time eigenenergies,eigenstates = diagonalization(system, cache_folder=cache_fold_name, load_cache=false)  #hide
end #hide
nothing; #hide
```

Diagonalizing the Hamiltonian is an expensive operation. For `j = 100` and `Nmax = 300` it
can take up to a day, but the function [`diagonalization`](@ref DickeModel.DickeBCE.diagonalization) saves the result to disk, so the second
time you call it with the same parameters it just loads it:
```@example examples
system = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j = 30) 
if on_github && use_current_dir_for_diags #hide
eigenenergies,eigenstates =  diagonalization(system, cache_folder=cache_fold_name,verbose=false) #hide
end #hide
if !use_current_dir_for_diags #hide
@time eigenenergies,eigenstates =  diagonalization(system)
else #hide
@time eigenenergies,eigenstates =  diagonalization(system,cache_folder=cache_fold_name) #hide
end #hide
nothing; #hide
```
Note that we did not have to pass `Nmax` this time, it loaded it from disk 
(see more details on the documentation of [`QuantumDickeSystem`](@ref DickeModel.DickeBCE.QuantumDickeSystem)). 

You can change the default folder, or disable caching altogether by passing extra arguments to [`diagonalization`](@ref DickeModel.DickeBCE.diagonalization).

The resulting `eigenstates` form a matrix. To get the ``k``th eigenstate,
simply call `state_k = eigenstates[:,k]` (or, even better, `state_k = @view eigenstates[:,k]`,
which [avoids unnecessary memory allocations.](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-views))

## [Local density of states of a coherent state](@id quantumldoscoherentstateex)

In this example, we obtain the eigenenergy components of a coherent state.
The function [`coherent_state`](@ref DickeModel.DickeBCE.coherent_state) will give us
a coherent state in the BCE, then we project it into the eigenbasis by left-multiplying
by `eigenstates'`, which is short for  [`adjoint`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#Base.adjoint)`(eigenstates)`.
Finally, we broadcast the [vectorized version of](https://docs.julialang.org/en/v1/manual/functions/#man-vectorized) the function [`abs2`](https://docs.julialang.org/en/v1/base/math/#BaseW.abs2) to extract all the coefficients.
```@setup examples
@info "Starting example: Coherent LDoS"
```
```@example examples
using Plots
using DickeModel.DickeBCE, DickeModel.ClassicalDicke
j = 30
system = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j = j, Nmax=120)
if false #hide
eigenenergies,eigenstates = diagonalization(system) 
end #hide

ϵₓ = -0.5
x = Point(system, Q=-1, P=0, p=0, ϵ=ϵₓ)
coh_state = coherent_state(system, x)
coherent_state_eigenbasis = eigenstates'*coh_state 
abscₖ²=abs2.(coherent_state_eigenbasis)
ϵₖs = eigenenergies/j

   
plot(histogram(ϵₖs, weights=abscₖ²,
        ylabel="Probability density", xticks=:none,normed=true,nbins=50),
    
    scatter(ϵₖs, abscₖ², ylabel="|cₖ|²", xlabel="ϵ"),
    
   size=(700,500), key=false, layout=(@layout [°;°]),
   xlim=(ϵₖs[1],ϵₖs[end]))
savefig("LDoS_quantum.svg");nothing #hide
```
![](LDoS_quantum.svg)

See [this example](@ref semiclassicalLDoS) for a semiclassical
computation of the envelope of this function.

## [Efficient Husimi functions](@id ExampleEfficientHusimiFunctions)

The functions [`DickeBCE.husimi`](@ref DickeModel.DickeBCE.husimi),  [`DickeBCE.coherent_overlap`](@ref DickeModel.DickeBCE.coherent_overlap),
and [`DickeBCE.coherent_state`](@ref DickeModel.DickeBCE.coherent_state) all accept a `chop` argument, which allows to significally speed
up computation time at the cost of slight numerical precision [Pilatowsky2020Notes](@cite). In this example we show how significant this speedup can be.
Let us construct a big system:
```@setup examples
@info "Starting example: Efficient Husimi"
```
```@example randstateHusimi
using DickeModel
using DickeModel.DickeBCE
using DickeModel.ClassicalDicke
using LinearAlgebra

j = 600
Nmax = 1200
system = QuantumDickeSystem(ω₀=1, ω=1, γ=1, j=j, Nmax=Nmax);
nothing; #hide
```
**Do not try to diagonalize such a big system! Your computer 💻 might explode 💥!**

For the sake of example, let us construct some random states in a simple manner (although
if you are interested in building random states in the eigenbasis, check the function
[`DickeBCE.random_state`](@ref DickeModel.DickeBCE.random_state) and see [this example](@ref ExampleRenyiOccupationsRandomStates)).

```@example randstateHusimi
n = 3 #how many random vectors
D = dimension(system)
random_vectors = rand(ComplexF64,(D,n))
for i in 1:n
    random_vectors[:,i] /= norm(@view random_vectors[:,i]) #normalize each one
end 
nothing; #hide
```
`random_vectors` is a matrix with `n` columns (states). Let us fix a point in the
phase space:
```@example randstateHusimi
x = Point(Q=0.6, P=-0.1, p=-0.2, q=-0.8)
husimi(system,x,random_vectors) #hide
nothing; #hide
```
We may call [`husimi`](@ref DickeModel.DickeBCE.husimi)`(system, x, random_vectors)`, which will
return an array with `n` elements. The `i`th element is the result of evaluating the Husimi function 
of the `i`th state (column) at the point `x`.
```@repl randstateHusimi
@time husimi(system, x, random_vectors, chop = 0)
```
By passing `chop = 0` we are allowing for no optimization. The code has to build all
the coefficients of the coherent state and then multiply them by each coefficient in `random_vectors`.
However, if we set `chop` to be slightly bigger, things will speed up significantly:
```@repl randstateHusimi
@time husimi(system, x, random_vectors, chop = 1e-14)
```
Note that the results barely changed, but this time it used a lot less memory and time.
The `chop` argument tells the code it can *chop* a portion of that size off the 
tails of the distribution of the coherent state (see Ref. [Pilatowsky2020Notes](@cite) for details). 
You loose almost no information, and you gain a lot of time. The default is `chop = 1e-6`, although you may increase it if you need more precision:
```@repl randstateHusimi
@time husimi(system, x, random_vectors) #default chop = 1e-6
```
That's fast!

## [Projected Wigner function of a cat state](@id wignerfuncexample)
    
Using [`DickeBCE.WignerProjqp`](@ref DickeModel.DickeBCE.WignerProjqp), we may compute
the Wigner function of a state, projected onto the atomic plane. We do this for a [cat state](https://en.wikipedia.org/wiki/Cat_state#Cat_states_in_single_modes) composed of two coherent states centered at `x` and `y`, which, taking advantage of [Julia's Unicode
capabilities](https://docs.julialang.org/en/v1/manual/unicode-input/), we name `🐱` (write `\:cat:` + Tab).
```@setup examples
@info "Starting example: Wigner Cat"
```
```@example examples
using DickeModel.ClassicalDicke, DickeModel.DickeBCE
system = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j=10, Nmax=50) 
res=0.025
if !on_github res=0.2 end #hide
Qs=Ps=-2:res:2
pts=[[Q,P] for Q in Qs, P in Ps if Q^2+P^2 <= 4]

x = Point(Q=-1.0, P=0, p=0, q=0)
y = Point(Q= 1.0, P=0, p=0, q=0)
🐱 = 1/sqrt(2) * (coherent_state(system, x) + coherent_state(system, y))

W=DickeBCE.WignerProjqp(system, 
                    [🐱], 
                    pts
                    ,show_progress = false, #hide
                    )[1]
d=Dict(zip(pts,W))
function mW(Q,P)
    if [Q,P] in pts
        return d[[Q,P]]
    else
        return NaN
    end
end
heatmap(Qs, Ps, mW, size=(600,600),
    xlabel = "Q", ylabel = "P",
    c=cgrad(:bwr, rev = true), clim=(-2.5,2.5))
savefig("catWigner.svg");nothing #hide
```
![](catWigner.svg)

!!! note
    The functions for computing Wigner functions are not thoroughly tested nor thoroughly optimized.
    They are based on these notes [Pilatowsky2019Notes](@cite), but they have room for improvement.
    
## Plotting the semiclassical density of states
Using [`DickeBCE.density_of_states`](@ref), we plot the semiclassical 
density of states originally calculated in Ref. [Bastarrachea2014](@cite).
Note that this function does not require diagonalization, so we can have
`j` as large as we want.
```@setup examples
@info "Starting example: DoS"
```
```@example examples
using DickeModel.ClassicalDicke
using DickeModel.DickeBCE

using Plots
system = QuantumDickeSystem(ω=1, γ=1, ω₀=1, j=100)

ν(ϵ) = density_of_states(system, ϵ)
ϵgs = minimum_energy(system)
plot(ν, ϵgs:0.01:2, xlabel="ϵ", ylabel="Density of States")
plot!(key=false) #hide
savefig("density_of_states.svg"); nothing #hide
```
![](density_of_states.svg)

This is precisely the red line in Fig. A1. of Ref. [Villasenor2020](@cite).