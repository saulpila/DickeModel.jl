# Examples for EnergyShellProjections
```@setup examples
push!(LOAD_PATH,"../../src")
on_github=get(ENV, "CI", nothing) == "true"
cache_fold_name="./diags"
use_current_dir_for_diags=on_github
using DickeModel
```
The module [`DickeModel.EnergyShellProjections`](@ref DickeModel.EnergyShellProjections) allows to integrate functions
in the classical energy shells of the Dicke model, and has specialized functions for projections of the
Husimi function and its moments.

## Projections of eigenstates

We generate a quantum system. The function [`EnergyShellProjections.proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.proj_husimi_QP_matrix) that 
we use below uses multiple workers if available, so let us load some as well:
```@example examples
using DickeModel
using DickeModel.DickeBCE
using DickeModel.ClassicalDicke
using DickeModel.EnergyShellProjections
using Distributed
using Plots
j = 30
Nmax = 120
system = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j=j, Nmax=Nmax)
if false #hide
eigenenergies,eigenstates = diagonalization(system) 
end #hide
if !use_current_dir_for_diags #hide
eigenenergies,eigenstates =  diagonalization(system,verbose=false) #hide
else #hide
eigenenergies,eigenstates =  diagonalization(system, cache_folder=cache_fold_name,verbose=false)  #hide
end #hide
ϵs=eigenenergies/j
if false #hide
addprocs(2) #we add 2 workers. Add as many as there are cores in your computer.
@everywhere using DickeModel
end #hide
nothing; #hide
```
The function [`EnergyShellProjections.proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.proj_husimi_QP_matrix) will make use
of all the available workers (you may disable this by passing `parallelize = true`)
!!! warning 
    The line 
    ```julia 
    @everywhere using DickeModel
    ```
    is necessary to load the module `DickeModel` in all workers. You will get errors if you omit it.
    
Let us first consider the projection
```math
    \iint \text{d} q\text{d} p \,\delta(h_\text{cl}(Q,q,P,p)-\epsilon_k)\, \mathcal{Q}_{k}(Q,q,P,p)
```
of the Husimi function ``\mathcal{Q}_{k}(\mathbf{x}) = \left | \left \langle \mathbf{x} \middle | E_k \right \rangle \right |^2``
of an eigenstate ``\left | E_k \right \rangle`` into the atomic plane intersected with the energy shell at ``\epsilon_k=E_k/j`` (see Ref. [Pilatowsky2021NatCommun](@cite)).

This can be easily plotted using  [`EnergyShellProjections.proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.proj_husimi_QP_matrix):

```@example examples
k = 600
state = eigenstates[:,k]  
ϵₖ = ϵs[k]
res = 0.02
if !on_github #hide
    res = 0.5 #hide 
end #hide
heatmap(
    EnergyShellProjections.proj_husimi_QP_matrix(system,
        state,
        ϵ = ϵₖ,
        symmetricQP = true,
        show_progress = false,#hide
        res = res)...,
    xlabel = "Q",
    ylabel = "P")
savefig("k600stateprojhu.svg");nothing #hide
```
![](k600stateprojhu.svg)

Moreover, using this function we can plot the projection of the ``\alpha \geq 0``-moments
```math
    \iint \text{d} q\text{d} p \,\delta(h_\text{cl}(\mathbf{x})-\epsilon_k)\, \mathcal{Q}_{k}(\mathbf{x})^\alpha.
```

```@example examples
k = 750
state = eigenstates[:,k]  
ϵₖ = ϵs[k]
powers = [0.5,1,2,3,4] 
res = 0.04
if !on_github #hide
    res = 0.5 #hide 
end #hide
Qs,Ps,matrices=EnergyShellProjections.proj_husimi_QP_matrix(system,
    state,
    ϵ = ϵₖ,
    show_progress = false, #hide
    matrix_powers = powers,
    symmetricQP = true,
    res = res)

plot(
    [heatmap(Qs,Ps,mat,
            key = false,
            xlabel = "Q",
            ylabel = if α==powers[1] "P" else "" end,
            yticks =  α==powers[1],
                left_margins = if α==powers[1] 20Plots.px else -5Plots.px end,
            right_margins= -5Plots.px,
            bottom_margins= 20Plots.px,
            clim = (0,NaN),
            title = "α = $α",
            
        ) 
        for (α,mat) in zip(powers,matrices)
    ]...,
    layout=(1,length(powers)),
    size=(1000,200),
    titlefontsize=10,
    tickfontsize=6,
)
savefig("k700momentsstateprojhu.svg");nothing #hide
```
![](k700momentsstateprojhu.svg)

The function [`proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.proj_husimi_QP_matrix) can also receive multiple
states as columns in a matrix. Even more, if these states are just vectors of length 4 `[Q,q,P,p]`, the function assumes you want the coherent 
state centered at `[Q,q,P,p]`. There is an analytical formula for the Husimi function of a coherent state (See [`husimi_of_coherent`](@ref DickeModel.DickeBCE.husimi_of_coherent)), so the coefficients of the coherent state are not even calculated, and the result is much faster:

```@example examples
ts = range(-π+0.3, -0.4, length = 30) ∪ range(0.2, π-0.6, length = 20)
💗(t) = (1.5sin(t)^3,(13cos(t) - 5cos(2t) -2cos(3t) - cos(4t))/10)
ϵ = 1
res = 0.02
if !on_github #hide
    res = 0.5 #hide 
end #hide
coherents = hcat([Point(system.classical_system,Q=Q,P=P,p=0,ϵ=ϵ) for (Q,P) in 💗.(ts)]...)
heatmap(EnergyShellProjections.proj_husimi_QP_matrix(system,coherents;
    mix_states = true,
    ϵ = ϵ,
    show_progress = false, #hide
    res = res),
    size = (600,600),
    color= :RdPu_9)
savefig("heartofcoherents.svg");nothing #hide
```
![](heartofcoherents.svg)

Nota that above, we passed `mix_states = true` to [`proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.proj_husimi_QP_matrix). This tells
the code to average together all of the
Husimis of the states (using [`Statistics.mean`](https://docs.julialang.org/en/v1/stdlib/Statistics/#Statistics.mean)).
You may even pass a  more complicated `mix_function` to add weights (see the documentation of [`proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.proj_husimi_QP_matrix) for details). If we had set `mix_states = false` (default), we would have obtained a matrix for each state. 

The fact that [`proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.proj_husimi_QP_matrix) may return the projection of multiple 
states at the same time allows to create really nice animations. We evolve the state using [`DickeBCE.evolve`](@ref DickeModel.DickeBCE.evolve).
```@example examples
ϵ = -0.5
x = Point(system.classical_system, Q=1, P=1, p=0, ϵ=ϵ)
coherent_x = coherent_state(system,x)
ts = 0:0.1:20
res = 0.1
if !on_github #hide
    ts = 0:2 #hide
    res= 0.5 #hide
end #hide
evolution = evolve(ts,coherent_x,eigenstates=eigenstates,eigenenergies=eigenenergies)
Qs,Ps,matrices=EnergyShellProjections.proj_husimi_QP_matrix(system,
    evolution,
    show_progress = false, #hide
    res = res,
    ϵ = ϵ)

animation=@animate for mat in matrices
   heatmap(Qs, Ps, mat,
        color = cgrad(:gist_heat, rev=true), size=(600,600),
        xlabel="Q", ylabel="P", key=false)
end

mp4(animation,
    "animation_of_evolution_Husimi.mp4",
    show_msg=false, #hide
    fps=15)
nothing; #hide
```
![](animation_of_evolution_Husimi.mp4)

!!! tip
    If you want better resolution, you may decrease `res` above. Computation time grows
    as the inverse cube of `res`. So twice the resolution will increase the computation time
    eightfold.


## Rényi occupation of random states
In this example, we construct a set of random states from the Gaussian Orthogonal
Ensemble (GOE) of Random Matrix Theory in the positive parity sector of the Dicke model 
using the function [`DickeBCE.random_state`](@ref DickeModel.DickeBCE.random_state). Then we study average Rényi
 Occupation (see Ref. [Villasenor2021](@cite)) using [`EnergyShellProjections.rényi_occupation`](@ref DickeModel.EnergyShellProjections.rényi_occupation).
 
```@example examples
using DickeModel
using DickeModel.DickeBCE
using DickeModel.ClassicalDicke
using DickeModel.EnergyShellProjections
using LinearAlgebra
using Statistics
using Plots
j = 30
Nmax = 120
system = QuantumDickeSystem(ω₀=1, ω=1, γ=1, j=j, Nmax=Nmax)
if false #hide
eigenenergies,eigenstates =  diagonalization(system)
end #hide 

ϵ = -0.5
n_of_states = 20
if !on_github #hide
    n_of_states = 2
end #hide
r_states = random_state(system, n_of_states,
    eigenenergies = eigenenergies,
    eigenstates = eigenstates,
    ϵ = ϵ, #center of gaussian envelope
    σ = 0.5,  #width of gaussian envelope
    ensemble = :GOE, #could be :GUE, try it!
    parity = +) #could be - or nothing, try it!

res = 0.05
if !on_github #hide
    res = 0.5
end #hide
αs=0:0.1:4
𝔏αs = rényi_occupation(system,
    r_states,
    ϵ = ϵ, 
    show_progress = false, #hide
    symmetricQP = true, # in Q because defined parity, in P because GOE
    res = res,
    α = αs)
av_𝔏αs = [mean(𝔏α) for 𝔏α in 𝔏αs]

plot(αs,av_𝔏αs,
    key = false,
    ylabel = "⟨𝔏α⟩",
    xlabel = "α",
    guidefont = "times")
savefig("randomStatesL.svg");nothing #hide
```
![](randomStatesL.svg)

!!! tip
    If you plan to compute both  [`EnergyShellProjections.rényi_occupation`](@ref DickeModel.EnergyShellProjections.rényi_occupation) and
    [`EnergyShellProjections.proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.proj_husimi_QP_matrix), you should use the 
    combined call  [`EnergyShellProjections.rényi_occupation_and_proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.rényi_occupation_and_proj_husimi_QP_matrix),
    which computes them both in the same routine and is faster than calling them separately.