# Examples for DickeBCE 

```@setup examples
push!(LOAD_PATH,"../../src")
on_github=get(ENV, "CI", nothing) == "true"
cache_fold_name="./diags"
use_current_dir_for_diags_and_remove=on_github
on_github=false
using Dicke
```
The module [`Dicke.DickeBCE`](@ref Dicke.DickeBCE) works with the quantum Dicke model
using a very efficient basis known as the coherent efficient basis (BCE for its acronym in Spanish).
See Refs. [Bastarrachea2014PSa](@cite) and [Bastarrachea2014PSb](@cite) for a detailed explanation on how and why it works. 
Throughout this examples, we will work with a system size of `j = 30`, but  using this module you can easily go up
to `j = 100` (and probably beyond, depending on the energy regime you are interested in studying).
## Diagonalizing the Dicke Hamiltonian
Let us start by defining our parameters:

```@example examples
using Dicke.DickeBCE, Dicke.ClassicalDicke
systemQ = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j=30, Nmax=120)
nothing; #hide
```

To load the eigenbasis, simply use [`diagonalization`](@ref Dicke.DickeBCE.diagonalization):

```@example examples
if !use_current_dir_for_diags_and_remove #hide
@time eigenenergies,eigenstates =  diagonalization(systemQ)
else #hide
@time eigenenergies,eigenstates =  diagonalization(systemQ, cache_folder=cache_fold_name)  #hide
end #hide
nothing; #hide
```

This saves the diagonalization to disk, so next time you can do:
```@example examples
systemQ = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j = 30) 
eigenenergies,eigenstates =  0,0 #hide
if !use_current_dir_for_diags_and_remove #hide
@time eigenenergies,eigenstates =  diagonalization(systemQ)
else #hide
@time eigenenergies,eigenstates =  diagonalization(systemQ,cache_folder=cache_fold_name) #hide
rm(cache_fold_name,recursive=true)#hide
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
using Dicke.DickeBCE, Dicke.ClassicalDicke
j = 30
systemC = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)
systemQ = QuantumDickeSystem(systemC, j = j, Nmax=120)
if false #hide
eigenenergies,eigenstates = diagonalization(systemQ) 
end #hide

ϵₓ = -0.5
x = Point(systemC, Q=-1, P=0, p=0, ϵ=ϵₓ)
coherent_state = coherentBCE(systemQ, x)
coherent_state_eigenbasis = transpose(eigenstates)*coherent_state 
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

## Survival probability of coherent state
Let us compare the survival probability of a coherent state, with
the truncated Wigner approximation. (See Ref. [Villasenor2020](@cite))
```@example examples
using Plots
import Dicke.TruncatedWignerApproximation
using Dicke.DickeBCE, Dicke.ClassicalDicke
j = 30
systemC = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)
systemQ = QuantumDickeSystem(systemC, j = j, Nmax=120)
if false #hide
eigenenergies,eigenstates =  diagonalization(systemQ) 
end #hide

x = Point(systemC, Q=1.75, P=0, p=0, ϵ=-0.5)
coherent_state=coherentBCE(systemQ, x)
ts=exp10.(-2:0.01:3)

N=20000
if !on_github N=1000 end
classical_SP = TruncatedWignerApproximation.survival_probability(
    systemC; 
    distribution=TruncatedWignerApproximation.coherent_Wigner_HWxSU2(x,j=j),
    N=N, ts=ts
)
quantum_SP = DickeBCE.survival_probability(
    ts,
    state=coherent_state, 
    eigenstates=eigenstates, 
    eigenenergies=eigenenergies
)

plot(ts, quantum_SP, yscale=:log10, xscale=:log10, 
    label="TWA", xlabel="time", ylabel="Survival probability")
plot!(ts, classical_SP, label="Quantum", ylim=(1e-4,1))
savefig("SP_QvsTWA.svg");nothing #hide
```
![](SP_QvsTWA.svg)

## [Efficient Husimi functions](@id exampletolhusimis)

The functions [`DickeBCE.Husimi`](@ref Dicke.DickeBCE.Husimi),  [`DickeBCE.coherentOverlap`](@ref Dicke.DickeBCE.coherentOverlap),
and [`DickeBCE.coherentBCE`](@ref Dicke.DickeBCE.coherentBCE) all accept a `tol` argument, which allows to significally speed
up computation time at the cost of slight numerical precision [Pilatowsky2020Notes](@cite). In this example we show how significant this speedup can be.
Let us construct a big system (Do not try to diagonalize such a big system! Your computer might explode!)
```@example randstateHusimi
using Dicke
using Dicke.DickeBCE
using Dicke.ClassicalDicke
using LinearAlgebra

j = 600
Nmax = 1200
system = QuantumDickeSystem(ω₀=1, ω=1, γ=1, j=j, Nmax=Nmax);
nothing; #hide
```
For the sake of example, let us construct some random in a simple manner (although
if you are interested in building random states in the eigenbasis, check the function
[`DickeBCE.randomState`](@ref Dicke.DickeBCE.randomState)).
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
nothing; #hide
```
We may call [`Husimi`](@ref Dicke.DickeBCE.Husimi)`(system, x, random_vectors)`, which will
return an array with `n` elements. The `i`th element is the result of evaluating the Husimi function 
of the `i`th state (column) at the point `x`.
```@setup randstateHusimi
Husimi(system,x,random_vectors)
```
```@repl randstateHusimi
@time Husimi(system, x, random_vectors, tol=0)
```
By passing `tol=0` we are allowing for no optimization. The code has to build all
the coefficients of the coherent state and then multiply them by each coefficient in `random_vectors`.
However, if we allow `tol` to be slightly bigger, things will speed up significantly:
```@repl randstateHusimi
@time Husimi(system, x, random_vectors, tol=1e-14)
```
Note that the results barely changed, but this time it used a lot less memory and time. The `tol` argument tells the code it can *chop* a portion of size `tol` off the tails of the distribution of the coherent state (see Ref. [Pilatowsky2020Notes](@cite) for details). You loose almost no information, 
and you gain a lot of time. The default is `tol = 1e-6`, which gives enough precision
for most purposes (although you may increase it if you need more precision):
```@repl randstateHusimi
@time Husimi(system, x, random_vectors) #default tol = 1e-6
```
That's fast!

## [Projected Wigner function of a cat state](@id wignerfuncexample)
Using [`DickeBCE.WignerProjqp`](@ref Dicke.DickeBCE.WignerProjqp), we may compute
the Wigner function of a state, projected onto the atomic plane.
Note: the functions for computing Wigner functions are not thoroughly tested.
They are based on these notes [Pilatowsky2019Notes](@cite)

```@example examples
using Dicke.ClassicalDicke, Dicke.DickeBCE
systemC = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)
systemQ = QuantumDickeSystem(systemC, j=10, Nmax=50) 
res=0.05
if !on_github res=0.2 end #end
Qs=Ps=-2:res:2
pts=[[Q,P] for Q in Qs, P in Ps if Q^2+P^2 <= 4]

x = Point(Q=-1.0, P=0, p=0, q=0)
y = Point(Q= 1.0, P=0, p=0, q=0)
🐱 = 1/sqrt(2) * (coherentBCE(systemQ, x) + coherentBCE(systemQ, y))

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