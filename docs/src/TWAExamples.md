# Examples for TWA 

```@setup examples
push!(LOAD_PATH,"../../src")
on_github=get(ENV, "CI", nothing) == "true"
cache_fold_name="./diags"
use_current_dir_for_diags=on_github
using DickeModel
```
## Classical evolution of coherent states

The module [`DickeModel.TWA`](@ref DickeModel.TWA) is a very powerful tool to study the classical evolution
of distributions in the phase space.

Let us load the module `DickeModel.TWA`, together with [`Distributed`](https://docs.julialang.org/en/v1/stdlib/Distributed/) which allows parallelization.  
```@setup examples
@info "Starting example: TWA"
```
```@example examples
using Distributed
using Plots,Plots.PlotMeasures
using DickeModel.TWA
using DickeModel.ClassicalDicke
if false #hide
addprocs(2) #we add 2 workers. Add as many as there are cores in your computer.
@everywhere using DickeModel
end #hide
```
The functions from [`DickeModel.TWA`](@ref DickeModel.TWA) will make use
of all the available workers.
!!! warning 
    The line 
    ```julia 
    @everywhere using DickeModel
    ```
    is necessary to load the module `DickeModel` in all workers. You will get errors if you omit it.
    
For our first example, let us consider the Wigner function of a coherent state,
evolve it classically using the truncated Wigner approximation (TWA) (see Refs. [Villasenor2020](@cite), [Pilatowsky2020](@cite)), and then look at 
the expected value of the [Weyl symbol](https://en.wikipedia.org/wiki/Wigner%E2%80%93Weyl_transform#The_inverse_map) 
of the observable ``\hat{j}_z=\hat{J}_z/j`` in time with 
[`TWA.average`](@ref). Note the usage of [`Weyl.Jz`](@ref DickeModel.TWA.Weyl.Jz).
```@setup examples
@info "Starting example: TWA jz"
```
```@example examples
system = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)
x = Point(system, Q=1, P=1, p=0, ϵ=0.5)
j = 300
W = coherent_Wigner_HWxSU2(x,j=j)
jz = Weyl.Jz(j)/j 
times = 0:0.05:40
if !on_github  times=0:1:4 end #hide

N = 20000
if !on_github N=200 end #hide
jz_average = average(system; 
    observable = jz, 
    show_progress = false, #hide
    distribution = W, 
    N = N, 
    ts=times)
plot(times,jz_average, xlabel="time", 
    ylabel="jz",key=false,ylim=(-1,1), size=(700,350))
savefig("average_jz_TWA.svg");nothing #hide
```
![](average_jz_TWA.svg)

Okay, but we can do more. Let's see how the whole distribution of ``{j}_z``
evolves classically using [`TWA.calculate_distribution`](@ref).
```@setup examples
@info "Starting example: jz distribution"
```
```@example examples
y_axis_values = -1.1:0.01:1.1
matrix = calculate_distribution(system; 
    distribution = W, 
    N = N,
    x=:t,  
    show_progress=false, #hide
    ts = times,
    y = jz, 
    ys = y_axis_values)
heatmap(times, y_axis_values, matrix,
    size=(700,350), color=cgrad(:gist_heat, rev=true),
    xlabel="time", ylabel="jz")
savefig("distribution_jz_TWA.svg");nothing #hide
```
![](distribution_jz_TWA.svg)

See [this](@ref TWAvsQuantum) example for a comparison between exact quantum evolution
and TWA.

We can chain several computations using [`TWA.mcs_chain`](@ref). 
For example, let's see the evolution of ``q`` and ``p`` for the same coherent state evolving in time, along with the time-averaged
distribution in the plane ``q,p``.
```@setup examples
@info "Starting example: TWA qp"
```
```@example examples
qs = -4.3:0.02:4.3
ps = -2.4:0.02:2.4
if !on_github  qs=-4.2:0.5:4.2;ps=-2.4:0.5:2.4 end #hide
N = 50000
if !on_github N=200 end #hide

mcs=mcs_chain(
    mcs_for_distributions(
        system; N = N, 
        distribution=W,
        y=:t,  ts=times,
        x=:q, xs=qs),
    mcs_for_distributions(
            system; N = N,
            distribution=W,
            y=:p, ys=ps,
            x=:t,  ts=times),
    mcs_for_distributions(
        system; N = N,
        distribution=W,
        x=:q,  xs=qs,
        y=:p, ys=ps,ts=times)
)
matrix_q_vs_t, 
matrix_t_vs_p, 
matrix_q_vs_p = monte_carlo_integrate(system, 
                      mcs;
                      ts = times,
                      N = N,
                      show_progress = false, #hide
                      distribution = W)

    
plot(heatmap(qs,times, matrix_q_vs_t,
        color=cgrad(:gist_heat, rev=true),
        ylabel="time", xlabel="q", xmirror =true, 
        ymirror =true, bottom_margin = -15mm),
    heatmap(times,ps, matrix_t_vs_p,
        color=cgrad(:gist_heat, rev=true),
        xlabel="time", ylabel="p", right_margin = -15mm),
    heatmap(qs,ps, matrix_q_vs_p,
        color=cgrad(:gist_heat, rev=true),
        ticks=:none,size=(400,400),margin = -15mm),
    layout=(@layout [_ °; 
                     ° °]), 
    color=cgrad(:gist_heat, rev=true),
    size=(800,800),colorbar=:none,link=:both)
savefig("distribution_qptime_TWA.svg");nothing #hide
```
![](distribution_qptime_TWA.svg)

The function [`calculate_distribution`](@ref TWA.calculate_distribution) can even
animate the evolution (with a little help from the wonderful [`@animate` from Plots](https://docs.juliaplots.org/latest/animations/)).
```@setup examples
@info "Starting example: TWA qp animation"
```
```@example examples
N = 200000
qs = -4.3:0.04:4.3
ps = -2.4:0.04:2.4
if !on_github  qs=-4.2:0.5:4.2;ps=-2.4:0.5:2.4 end #hide

times = 0:0.1:40
if !on_github N=200 end #hide
if !on_github times=0:1 end #hide

matrices = calculate_distribution(system; 
    distribution=W, 
    N = N,
    x = :q, 
    y = :p,
    xs = qs, 
    ys = ps, 
    ts = times, 
    show_progress = false, #hide
    animate = true)
animation=@animate for mat in matrices
    heatmap(qs, ps, mat,
        color = cgrad(:gist_heat, rev=true), size=(600,600), 
        xlabel="q", ylabel="p", key=false)
end
mp4(animation,
    "animation_of_evolution.mp4",
    show_msg=false, #hide
    fps=30)
nothing; #hide
```
![](animation_of_evolution.mp4)
!!! note
    Computing animations with [`calculate_distribution(..., animate = true, ...)`](@ref TWA.calculate_distribution)
    may need a lot of RAM. You can estimate the maximum amount of RAM needed using the shorthand formula  
    
    ``(\text{\# of workers}) \times (\text{length of }`` `xs` ``)\times (\text{length of }`` `ys` ``)\times (\text{length of }`` `ts` ``) \times (64 \text{ bits})``.
    
    But this number would only be reached if trajectories filled all of the matrices in all of the workers at all the timesteps. You may stay much
    below this number by passing `maxNBatch` to [`calculate_distribution`](@ref TWA.calculate_distribution) (or
    to  [`monte_carlo_integrate`](@ref TWA.monte_carlo_integrate)). This parameter limits the number of trajectories
    that are calculated in batch in each  worker. Between batches, data is flushed to the main worker, which takes time, but liberates RAM. 
    If generating the animation is filling up your RAM, try to decrease `maxNBatch`.
    

## [TWA vs exact quantum evolution](@id TWAvsQuantum)
Let us compare the evolution of a quantum state with that given by the TWA. 

We select a point in the phase space and load
the Wigner function of a coherent state centered at that point using [`TWA.coherent_Wigner_HWxSU2`](@ref TWA.coherent_Wigner_HWxSU2)
```@setup examples
@info "Starting example: Quantum vs TWA"
```
```@example examples
using Plots
using DickeModel.TWA
using DickeModel.DickeBCE, DickeModel.ClassicalDicke
using LinearAlgebra
j = 30
system = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j = j, Nmax=120)
if false #hide
eigenenergies,eigenstates = diagonalization(system) 
end #hide
if !use_current_dir_for_diags #hide
eigenenergies,eigenstates =  diagonalization(system,verbose=false) #hide
else #hide
eigenenergies,eigenstates =  diagonalization(system, cache_folder=cache_fold_name,verbose=false)  #hide
end #hide

x = Point(system, Q=1.75, P=0, p=0, ϵ=-0.5)
coh_state=coherent_state(system, x)
W = coherent_Wigner_HWxSU2(x,j=j)
nothing #hide
```

First, we compare the expectation value of the observable ``\hat{J}_z^2``.
We load its Weyl symbol using [`Weyl.Jz²(j)`](@ref DickeModel.TWA.Weyl.Jz²), and its
matrix form in the BCE using [`DickeBCE.Jz(system)`](@ref DickeModel.TWA.Weyl.Jz²)`^2`.
```@example examples
ts= 0:0.05:40
evolution = evolve(ts, coh_state, 
                eigenstates=eigenstates,
                eigenenergies=eigenenergies);
Jz²=DickeBCE.Jz(system)^2

               #〈 v｜Jz²｜v 〉for each column v in evolution
exvals = [real(dot(v, Jz² ,v)) for v in eachcol(evolution)]

N=20000
if !on_github N=1000 end #hide
TWAJz2 = TWA.average(system,
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
We use the function [`DickeBCE.survival_probability`](@ref) to compute the exact
survival probability, and [`TWA.survival_probability`](@ref) gives us the result from
the TWA.
```@example examples
ts=exp10.(-2:0.01:3)

N=20000 # Higher values reduce numerical noise at the cost of speed
if !on_github N=1000 end #hide
classical_SP = TWA.survival_probability(
    system; 
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

## Fidelity out-of-time order correlator (FOTOC)

The FOTOC is a quantum equivalent of the classical Lyapunov exponent. It is just
the variance ``\text{var}(Q)+\text{var}(q)+\text{var}(P)+\text{var}(p)`` as a function
of time. It may be calculated using the TWA. 
(See Ref. [Pilatowsky2020](@cite) and references therein).
```@setup examples
@info "Starting example: TWA FOTOC Dicke"
```
```@example examples
using DickeModel.ClassicalDicke
using DickeModel.ClassicalSystems
using DickeModel.TWA
using Plots
system = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)

ts = 0:0.1:50
j = 1000
x = Point(system, Q=1, P=0, p=0, ϵ=-0.6)
W = coherent_Wigner_HWxSU2(x, j=j)
N = 10000
if !on_github N=200 end #hide
if !on_github times=0:1:10 end #hide

FOTOC=sum.(
    variance(system; 
        observable = [:Q,:q,:P,:p], 
        distribution = W, 
        N=N, 
        ts=ts, 
        show_progress = false, #hide
        tol=1e-8)
            )

plot(ts, FOTOC, 
     xlabel="time", ylabel="FOTOC",
     label="FOTOC", yscale=:log10)
if false #hide
lyapunov = lyapunov_exponent(system, x)
else #hide
lyapunov = lyapunov_exponent(system, x; verbose = false) #hide
end #hide
plot!(t->exp(2*lyapunov*t)*2/j, 0, 14, 
    label="exp(2λt)*2ħ", key=:bottomright)
savefig("FOTOC_TWA.svg");nothing #hide
```
![](FOTOC_TWA.svg)

## [Energy profiles of a coherent state](@id semiclassicalLDoS)
We have a semiclassical formula for the energy width of a coherent state, given in App. A of Ref. [Lerma2018](@cite),
and implemented in [`DickeBCE.energy_width_of_coherent_state`](@ref). Let's check
this formula against the semiclassical local density of states given by Eq. (E.3) of Ref. [Villasenor2020](@cite).
```@setup examples
@info "Starting example: Semiclassical LDoS"
```
```@example examples
using DickeModel.ClassicalDicke
using DickeModel.ClassicalSystems
using DickeModel.TWA
using DickeModel.DickeBCE
using Distributions
using Plots
j = 1000
system = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j = j)

ts = 0:0.1:50
ϵₓ = -0.5
x = Point(system, Q=-1, P=0, p=0, ϵ=ϵₓ)
W = coherent_Wigner_HWxSU2(x, j=j)
N = 500000

ϵ_binsize = 0.02
if !on_github N=1000 end #hide
if !on_github ϵ_binsize=0.1 end #hide

ϵs = -0.8:ϵ_binsize:0
ρ = calculate_distribution(system, 
    distribution = W, 
    N = N, 
    show_progress = false, #hide
    x = ClassicalDicke.hamiltonian(system), 
    xs = ϵs)'
ρ /= sum(ρ)*ϵ_binsize #normalization
σₓ = energy_width_of_coherent_state(system, x)
gaussian=Distributions.Normal(ϵₓ, σₓ)

plot(ϵs,ρ, label="∫ w(x) δ(ϵ - h(x)) dx")
plot!(ϵ->pdf(gaussian,ϵ), ϵs, 
    label="normal(σ)", linestyle=:dash,
    xlabel="ϵ", ylabel="Probability density")
savefig("LDoS_classical.svg");nothing #hide
```
![](LDoS_classical.svg)

See [this example](@ref quantumldoscoherentstateex) for the
full quantum computation of the energy spectrum of a coherent state.