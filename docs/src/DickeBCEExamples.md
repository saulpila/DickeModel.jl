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
to `j = 100`, as done in Refs. [Pilatowsky2021](@cite), [Pilatowsky2021NatCommun](@cite),  [Villasenor2021](@cite).

## Diagonalizing the Dicke Hamiltonian
Let us start by defining our parameters:

```@example examples
using DickeModel.DickeBCE, DickeModel.ClassicalDicke
systemQ = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j=30, Nmax=120)
nothing; #hide
```

To load the eigenbasis, simply use [`diagonalization`](@ref DickeModel.DickeBCE.diagonalization):

```@example examples
if !use_current_dir_for_diags #hide
@time eigenenergies,eigenstates =  diagonalization(systemQ)
else #hide
@time eigenenergies,eigenstates =  diagonalization(systemQ, cache_folder=cache_fold_name, load_cache=false)  #hide
end #hide
nothing; #hide
```

This saves the diagonalization to disk, so next time you can do:
```@example examples
systemQ = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j = 30) 
eigenenergies,eigenstates =  0,0 #hide
if !use_current_dir_for_diags #hide
@time eigenenergies,eigenstates =  diagonalization(systemQ)
else #hide
@time eigenenergies,eigenstates =  diagonalization(systemQ,cache_folder=cache_fold_name) #hide
end #hide
nothing; #hide
```
(Note that we did not have to pass `Nmax` this time, it loaded it from disk.)

The resulting `eigenstates` form a matrix. To get the ``k``th eigenstate,
simply do `state_k = eigenstates[:,k]`.

## [Local density of states of a coherent state](@id quantumldoscoherentstateex)
In this example, we obtain the eigenenergy components of a coherent state.
```@example examples
using Plots
using DickeModel.DickeBCE, DickeModel.ClassicalDicke
j = 30
systemC = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)
systemQ = QuantumDickeSystem(systemC, j = j, Nmax=120)
if false #hide
eigenenergies,eigenstates = diagonalization(systemQ) 
end #hide

ϵₓ = -0.5
x = Point(systemC, Q=-1, P=0, p=0, ϵ=ϵₓ)
coh_state = coherent_state(systemQ, x)
coherent_state_eigenbasis = transpose(eigenstates)*coh_state 
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

## [Evolution of coherent state vs TWA](@id TWAvsQuantum)
Let us compare the evolution of a quantum state with that given by truncated Wigner
approximation (see Refs. [Villasenor2020](@cite), [Pilatowsky2020](@cite)).
```@example examples
using Plots
using DickeModel.TWA
using DickeModel.DickeBCE, DickeModel.ClassicalDicke
using LinearAlgebra
j = 30
systemC = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)
systemQ = QuantumDickeSystem(systemC, j = j, Nmax=120)
if false #hide
eigenenergies,eigenstates =  diagonalization(systemQ) 
end #hide

x = Point(systemC, Q=1.75, P=0, p=0, ϵ=-0.5)
coh_state=coherent_state(systemQ, x)
W = coherent_Wigner_HWxSU2(x,j=j)
nothing #hide 
```

First, we compare the expectation value of the observable ``\hat{J}_z^2``.
Note that we use [`Weyl.Jz²`](@ref DickeModel.TWA.Weyl.Jz²)`(j)`,
which is not the same as [`Weyl.Jz`](@ref DickeModel.TWA.Weyl.Jz)`(j)^2`.
```@example examples
ts= 0:0.05:40
evolution = evolve(ts, coh_state, 
                eigenstates=eigenstates,
                eigenenergies=eigenenergies);
Jz²=DickeBCE.Jz(systemQ)^2
exvals = [real(dot(v,Jz²,v)) for v in eachcol(evolution)]

N=20000
if !on_github N=1000 end #hide
TWAJz2 = TWA.average(systemC,
                 distribution = W,
                 show_progress=false, #hide
                 observable = Weyl.Jz²(j), 
                 ts = ts,
                 N = N)

plot(ts, [exvals TWAJz2], 
    size=(700,350), label=["Quantum" "TWA"],
    left_margin=2Plots.mm,
    xlabel = "time", ylabel="Jz²")
savefig("Jz2_QvsTWA.svg");nothing #hide
```
![](Jz2_QvsTWA.svg)

Now let us take a look at the survival probability (see Ref. [Villasenor2020](@cite)).
```@example examples
ts=exp10.(-2:0.01:3)

N=20000 # Higher values reduce numerical noise at the cost of speed
if !on_github N=1000 end #hide
classical_SP = TWA.survival_probability(
    systemC; 
    distribution = W,
    show_progress=false, #hide
    N=N, ts=ts
)
quantum_SP = DickeBCE.survival_probability(
    ts,
    state=coh_state, 
    eigenstates=eigenstates, 
    eigenenergies=eigenenergies
)

plot(ts, [quantum_SP classical_SP], 
    yscale=:log10, xscale=:log10, 
    ylim=(1e-4,1), label=["Quantum" "TWA"], 
    xlabel="time", ylabel="Survival probability")
savefig("SP_QvsTWA.svg");nothing #hide
```
![](SP_QvsTWA.svg)

## [Efficient Husimi functions](@id exampletolhusimis)

The functions [`DickeBCE.husimi`](@ref DickeModel.DickeBCE.husimi),  [`DickeBCE.coherent_overlap`](@ref DickeModel.DickeBCE.coherent_overlap),
and [`DickeBCE.coherent_state`](@ref DickeModel.DickeBCE.coherent_state) all accept a `tol` argument, which allows to significally speed
up computation time at the cost of slight numerical precision [Pilatowsky2020Notes](@cite). In this example we show how significant this speedup can be.
Let us construct a big system:
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
**Do not try to diagonalize such a big system! Your computer might explode!**

For the sake of example, let us construct some random states in a simple manner (although
if you are interested in building random states in the eigenbasis, check the function
[`DickeBCE.random_state`](@ref DickeModel.DickeBCE.random_state)).
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
@time husimi(system, x, random_vectors, tol=0)
```
By passing `tol=0` we are allowing for no optimization. The code has to build all
the coefficients of the coherent state and then multiply them by each coefficient in `random_vectors`.
However, if we allow `tol` to be slightly bigger, things will speed up significantly:
```@repl randstateHusimi
@time husimi(system, x, random_vectors, tol=1e-14)
```
Note that the results barely changed, but this time it used a lot less memory and time. The `tol` argument tells the code it can *chop* a portion of size `tol` off the tails of the distribution of the coherent state (see Ref. [Pilatowsky2020Notes](@cite) for details). You loose almost no information, 
and you gain a lot of time. The default is `tol = 1e-6`, which gives enough precision
for most purposes (although you may increase it if you need more precision):
```@repl randstateHusimi
@time husimi(system, x, random_vectors) #default tol = 1e-6
```
That's fast!

## [Projected Wigner function of a cat state](@id wignerfuncexample)
Using [`DickeBCE.WignerProjqp`](@ref DickeModel.DickeBCE.WignerProjqp), we may compute
the Wigner function of a state, projected onto the atomic plane.
Note: the functions for computing Wigner functions are not thoroughly tested.
They are based on these notes [Pilatowsky2019Notes](@cite)

```@example examples
using DickeModel.ClassicalDicke, DickeModel.DickeBCE
systemC = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)
systemQ = QuantumDickeSystem(systemC, j=10, Nmax=50) 
res=0.025
if !on_github res=0.2 end #hide
Qs=Ps=-2:res:2
pts=[[Q,P] for Q in Qs, P in Ps if Q^2+P^2 <= 4]

x = Point(Q=-1.0, P=0, p=0, q=0)
y = Point(Q= 1.0, P=0, p=0, q=0)
🐱 = 1/sqrt(2) * (coherent_state(systemQ, x) + coherent_state(systemQ, y))

W=DickeBCE.WignerProjqp(systemQ, 
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