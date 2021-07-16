# Dicke.ClassicalDicke
```@meta
CurrentModule = Dicke.ClassicalDicke
```
```@docs 
ClassicalDickeSystem
hamiltonian
```

## Phase space
```@docs 
Point
Pointθϕ
energy_shell_volume
phase_space_dist_squared
classical_path_random_sampler
```

## Energy minimizing
```@docs 
minimum_energy
minimum_ϵ_for
minimum_energy_point
normal_frequency
```

## Semiclassical approximations
```@docs 
density_of_states
energy_width_of_coherent_state
```

## Roots in ``q``
```@docs 
discriminant_of_q_solution
q_of_ϵ
q_sign
Point(::ClassicalDickeSystem)
```
## Atomic boundary of energy shells

```@setup diagram
push!(LOAD_PATH,"../../src")
using Dicke.ClassicalDicke
using Plots
system=ClassicalDicke.ClassicalDickeSystem(ω=1, γ=1, ω₀=1)
Qs=-2:0.01:2
Ps=-1.5:0.01:1.5
ϵ=-1.1
contour(Qs,Ps,(Q,P)-> ClassicalDicke.minimum_ϵ_for(system,p=0,P=P,Q=Q),
    levels=[ϵ],xlabel="Q",ylabel="P",ylim=(-1.5,2.5),color=:black,key=false,linewidth=2,
    title="Available phase space at ϵ")
Qmin=minimum_nonnegative_Q_for_ϵ(system,ϵ)
Qmax=maximum_Q_for_ϵ(system,ϵ)
Pmax=maximum_P_for_ϵ(system,ϵ)
plot!([(Qmin,-2),(Qmin,2)],label="+ minimum_nonnegative_Q_for_ϵ",key=:top)
plot!([(-Qmin,-2),(-Qmin,2)],label="− minimum_nonnegative_Q_for_ϵ",linestyle=:dash)

plot!([(Qmax,-2),(Qmax,2)],label="+ maximum_Q_for_ϵ")
plot!([(-Qmax,-2),(-Qmax,2)],label="− maximum_Q_for_ϵ",linestyle=:dash)

plot!([(-2,Pmax),(2,Pmax)],label="+ maximum_P_for_ϵ")
diaplot=plot!([(-2,-Pmax),(2,-Pmax)],label="− maximum_P_for_ϵ",linestyle=:dash)
```
```@example diagram
diaplot #hide
```

```@docs 
minimum_nonnegative_Q_for_ϵ
maximum_Q_for_ϵ
maximum_P_for_ϵ
```
