# Examples
```@setup examples
push!(LOAD_PATH,"../../src")
using Dicke
```
## (Semi)classical Dicke model

### Drawing contours of the available phase space
We may use the function [`ClassicalDicke.minimum_ε_for`](@ref) to draw the contour of the available phase space on the variables
``(Q,P)``.
```@example examples
using Dicke
using Plots
system=  ClassicalDicke.ClassicalSystem(ω=1, γ=1, ω₀=1)
Qs = Ps = -2:0.01:2
εgs = ClassicalDicke.minimum_energy(system)
contour(Qs, Ps,
        (Q,P) -> ClassicalDicke.minimum_ε_for(system, p=0, P=P, Q=Q),
        levels=10, clim=(εgs,1), xlabel="Q", ylabel="P")
savefig("contourQP.svg"); nothing #hide
```
![](contourQP.svg)
### Plotting the density of states
Here is a plot of the semiclassical density of states

```@example examples
using Dicke
using Plots
system = ClassicalDicke.ClassicalSystem(ω=1, γ=1, ω₀=1)
ν(ε) = ClassicalDicke.density_of_states(system, j=100, ε=ε)
εgs = ClassicalDicke.minimum_energy(system)
plot(ν, εgs, 2, xlabel="ε", ylabel="Density of States")
plot!(key=false) #hide
savefig("density_of_states.svg"); nothing #hide
```
![](density_of_states.svg)

This is precisely the red line in Fig. A1. of Ref. [Villasenor2020](@cite).

### Drawing a Poincaré surface

Here is a way to draw a [Poincaré surface](https://en.wikipedia.org/wiki/Poincar%C3%A9_map) for the Dicke model. We use [`ClassicalSystems.integrate`](@ref) to integrate a bunch of initial conditions. Using the [callback system of DifferentialEquations](https://diffeq.sciml.ai/stable/features/callback_functions/#DiffEqBase.ContinuousCallback), we save the points where ``p=0``.
```@example examples
using Dicke
using Plots
using DiffEqBase

system=ClassicalDicke.ClassicalSystem(ω=0.8, γ=0.8, ω₀=1)
mplot=scatter(fmt=:png,key=false,markerstrokecolor = :transparent,
  markersize=1,legend=false,size=(500,500),
  color_palette=:darkrainbow,xlabel="Q",ylabel="P")

pts=Tuple{Float64, Float64}[] #a list of points (Q,P)
callback=ContinuousCallback((x,t,_)-> x[4], #when p=x[4] is 0,
    state->push!(pts,(state.u[1],state.u[3])),  #save Q=u[1], P=u[3]
    nothing;
    save_positions=(false,false),abstol=1e-6)
ε=-1.35
maxQ = ClassicalDicke.maximum_Q_for_ε(system,ε)
for Q in 0:0.02:maxQ #for a bunch of initial Qs,
        if ClassicalDicke.minimum_ε_for(system,P=0,p=0,Q=Q) > ε #we are outside the bounds
            continue
        end
        initial_condition = ClassicalDicke.Point(system, ε=ε, P=0, p=0, Q=Q)
        ClassicalSystems.integrate(system,u₀=initial_condition,
            t=10000,cb=callback,save_everystep=false,tol=1e-9)
        scatter!(pts)
        empty!(pts)
end
mplot
savefig("poincare_surface.png");nothing #hide
```
![](poincare_surface.png)
