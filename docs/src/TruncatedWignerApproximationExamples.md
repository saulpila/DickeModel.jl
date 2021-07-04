# Examples for TrucatedWingerApproximation 

```@setup examples
push!(LOAD_PATH,"../../src")
on_github=get(ENV, "CI", nothing) == "true"
using Dicke
```
## Classical evolution of coherent states

The module [`Dicke.TruncatedWignerApproximation`](@ref Dicke.TruncatedWignerApproximation) is a very powerful tool to study the classical evolution
of distributions in the phase space.

Let us load the module `Dicke.TruncatedWignerApproximation`, together with [`Distributed`](https://docs.julialang.org/en/v1/stdlib/Distributed/) which allows paralelization.  
```@example examples
using Distributed
using Plots,Plots.PlotMeasures
using Dicke,Dicke.TruncatedWignerApproximation, Dicke.ClassicalDicke
if false #hide
addprocs(2) #we add 2 workers. Add as many as there are cores in your computer.
@everywhere using Dicke,Dicke.TruncatedWignerApproximation
end #hide
```
The functions from [`Dicke.TruncatedWignerApproximation`](@ref Dicke.TruncatedWignerApproximation) will make use
of all the available workers.

For our first example, let us consider the Wigner function of a coherent state,
evolve it classically using the truncated Wigner approximation, and then look at 
the expected value of the Weyl symbol of the observable ``\hat{j}_z=\hat{J}_z/j`` in time with 
[`TruncatedWignerApproximation.average`](@ref). Note the usage of [`Weyl.Jz`](@ref Dicke.TruncatedWignerApproximation.Weyl.Jz).

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
jz_average=average(system;observable=jz, distribution=W, N = N, ts=times)
plot(times,jz_average, xlabel="time", 
    ylabel="jz",key=false,ylim=(-1,1), size=(700,350))
savefig("average_jz_TWA.svg");nothing #hide
```
![](average_jz_TWA.svg)

Okay, but we can do more. Let's see how the whole distribution of ``{j}_z``
evolves classically using [`TruncatedWignerApproximation.calculate_distribution`](@ref).

```@example examples
y_axis_values = -1.1:0.01:1.1
matrix = calculate_distribution(system; distribution=W, N = N,
    x=:t,  ts=times,
    y=jz, ys=y_axis_values)
heatmap(times, y_axis_values, matrix,
    size=(700,350), color=cgrad(:gist_heat, rev=true),
    xlabel="time", ylabel="jz")
savefig("distribution_jz_TWA.svg");nothing #hide
```
![](distribution_jz_TWA.svg)

We can chain several computations using [`TruncatedWignerApproximation.mcs_chain`](@ref). 
For example, let's see the evolution of ``q`` and ``p`` for the same coherent state evolving in time, along with the time-averaged
distribution in the plane ``q,p``.
```@example examples
qs=-4.2:0.05:4.2
ps=-2.4:0.05:2.4
if !on_github  qs=-4.2:0.5:4.2;ps=-2.4:0.5:2.4 end #hide
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
matrix_q_vs_t,matrix_t_vs_p,matrix_q_vs_p = monte_carlo_integrate(system,
    mcs;ts=times,N=N,distribution=W,tolerate_errors=false)
    
    
plot(heatmap(qs,times, matrix_q_vs_t,
        color=cgrad(:gist_heat, rev=true),
        ylabel="time", xlabel="q",xmirror =true,ymirror =true, bottom_margin = -15mm),
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
savefig("distribution_qptime_TWA.svg");nothing #hide.
```
![](distribution_qptime_TWA.svg)

The function [`calculate_distribution`](@ref TruncatedWignerApproximation.calculate_distribution) can even
animate the evolution (with a little help from the wonderful [`@animate` from Plots](https://docs.juliaplots.org/latest/animations/)).
```@example examples
N = 1000000
times = 0:0.1:40
if !on_github N=200 end #hide
if !on_github times=0:1 end #hide

matrices = calculate_distribution(system; distribution=W, N = N,
    x=:q,y=:p,xs=qs, ys=ps,ts=times,animate=true,maxNBatch=200000);
animation=@animate for mat in matrices
    heatmap(qs, ps, mat,
        color=cgrad(:gist_heat, rev=true),size=(600,600),xlabel="q",ylabel="p",key=false)
end
mp4(animation,
    "animation_of_evolution.mp4",
    show_msg=false, #hide
    fps=30)
nothing; #hide
```
![](animation_of_evolution.mp4)
## Fidelity out-of-time order correlator (FOTOC)

The FOTOC is a quantum-equivalent of the classcal Lyapunov exponent. It is just
the variance ``\text{var}(Q)+\text{var}(q)+\text{var}(P)+\text{var}(p)`` as a function
of time. It may be calculated using the TWA. 
(See Ref. [Pilatowsky2020](@cite) and references therein).

```@example examples
using Dicke.ClassicalDicke
using Dicke.ClassicalSystems
using Dicke.TruncatedWignerApproximation
using Plots
system = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)

ts = 0:0.1:50
j = 1000
x = Point(system, Q=-1, P=0, p=0, ϵ=-0.6)
W = coherent_Wigner_HWxSU2(x, j=j)
N = 10000
if !on_github N=200 end #hide
if !on_github times=0:1:10 end #hide

FOTOC=sum.(variance(system; observable=[:Q,:p,:P,:p], 
                    distribution=W, N=N, ts=ts, tol=1e-8))

plot(ts, FOTOC, 
     xlabel="time", ylabel="FOTOC",
     label="FOTOC", yscale=:log10)
lyapunov = lyapunov_exponent(system, u₀=x)
plot!(t->exp(2*lyapunov*t)*2/j, 0, 14, 
    label="exp(2λt)*2ħ", key=:bottomright)
savefig("FOTOC_TWA.svg");nothing #hide.
```
![](FOTOC_TWA.svg)

## Energy profiles of coherent states
We have a semiclassical formula for the energy width of a coherent state, given in App. A of Ref. [Lerma2018](@cite),
and implemented in [`ClassicalDicke.energy_width_of_coherent_state`](@ref). Let's check
this formula against the semiclassical local density of states given by Eq. (E.3) of Ref. [Villasenor2020](@cite).

```@example examples
using Dicke.ClassicalDicke
using Dicke.ClassicalSystems
using Dicke.TruncatedWignerApproximation
using Distributions
using Plots
system = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)

ts = 0:0.1:50
j = 1000
ϵₓ = -0.5
x = Point(system, Q=-1, P=0, p=0, ϵ=ϵₓ)
W = coherent_Wigner_HWxSU2(x, j=j)
N = 1000000

hamiltonian = ClassicalDicke.hamiltonian(system)
ϵ_binsize = 0.01
if !on_github N=1000 end #hide
if !on_github ϵ_binsize=0.1 end #hide

ϵs = -0.8:ϵ_binsize:0
ρ = calculate_distribution(system, distribution=W, 
    N=N,x=hamiltonian,xs=ϵs)'
ρ /= sum(ρ)*ϵ_binsize #normalization
σₓ = energy_width_of_coherent_state(system, x, j)
gaussian=Distributions.Normal(ϵₓ, σₓ)

plot(ϵs,ρ, label="∫ w(x) δ(ϵ - h(x)) dx")
plot!(ϵ->pdf(gaussian,ϵ), ϵs, 
    label="normal(σ)", linestyle=:dash,
    xlabel="ϵ", ylabel="Probability density")
savefig("LDoS_classical.svg");nothing #hide.
```
![](LDoS_classical.svg)