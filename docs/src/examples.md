# Examples
```@setup examples
push!(LOAD_PATH,"../../src")
on_github=get(ENV, "CI", nothing) == "true"
on_github=false
using Dicke
```
## (Semi)classical Dicke model

### Drawing contours of the available phase space
We may use the function [`ClassicalDicke.minimum_ε_for`](@ref) to draw the contour of the available phase space on the variables
``(Q,P)``.
```@example examples
using Dicke.ClassicalDicke
using Plots
system =  ClassicalDicke.ClassicalSystem(ω=1, γ=1, ω₀=1)
Qs = Ps = -2:0.01:2
εgs = minimum_energy(system)
contour(Qs, Ps,
        (Q,P) -> minimum_ε_for(system, p=0, P=P, Q=Q),
        levels=10, clim=(εgs,1), xlabel="Q", ylabel="P")
savefig("contourQP.svg"); nothing #hide
```
![](contourQP.svg)
### Plotting the density of states
Here is a plot of the semiclassical density of states

```@example examples
using Dicke.ClassicalDicke
using Plots
system = ClassicalDicke.ClassicalSystem(ω=1, γ=1, ω₀=1)
ν(ε) = density_of_states(system, j=100, ε=ε)
εgs = minimum_energy(system)
plot(ν, εgs, 2, xlabel="ε", ylabel="Density of States")
plot!(key=false) #hide
savefig("density_of_states.svg"); nothing #hide
```
![](density_of_states.svg)

This is precisely the red line in Fig. A1. of Ref. [Villasenor2020](@cite).

### Drawing a Poincaré surface

Here is a way to draw a [Poincaré surface](https://en.wikipedia.org/wiki/Poincar%C3%A9_map) for the Dicke model. We use [`ClassicalSystems.integrate`](@ref) to integrate a bunch of initial conditions. Using the [callback system of DifferentialEquations](https://diffeq.sciml.ai/stable/features/callback_functions/#DiffEqBase.ContinuousCallback), we save the points where ``p=0``.
```@example examples
using Dicke, Dicke.ClassicalDicke, Dicke.ClassicalSystems
using Plots
using DifferentialEquations

system=ClassicalDicke.ClassicalSystem(ω=0.8,γ=0.8,ω₀=1)
mplot=scatter(fmt=:png, key=false, markersize=1, legend=false,
    size=(500,500), color_palette=:darkrainbow, xlabel="Q", ylabel="P") 

pts = Tuple{Float64, Float64}[] #a list of points (Q,P)
function save(state) #this function saves (Q,P) to pts if q = q₊ (and not q₋).
    Q,q,P,p = state.u  
    q₊ = q_of_ε(system; Q=Q, P=P, p=p, ε=ε, signo=+)
    q₋ = q_of_ε(system; Q=Q, P=P, p=p, ε=ε, signo=-)
    if abs(q - q₊) < abs(q - q₋) #we have to compare approximatelly because the equalities...
        push!(pts, (Q,P))    # ...may not be exact due to numerical error.
    end                     
end
    
callback=ContinuousCallback((x, t, _) -> x[4], #when p=x[4] is 0,
    save; #execute the function save
    save_positions=(false,false), abstol=1e-6)
ε = -1.35
for Q in 0:0.02:maximum_Q_for_ε(system, ε) #for a bunch of initial Qs,
        if minimum_ε_for(system, P=0, p=0, Q=Q) > ε
            continue
        end
        initial_condition = Point(system, ε=ε, P=0, p=0, Q=Q)
        integrate(system, u₀=initial_condition,
            t=10000, callback=callback, save_everystep=false)
        scatter!(pts)
        empty!(pts)

end
mplot
savefig("poincare_surface.png");nothing #hide
```
![](poincare_surface.png)

### Drawing a Lyapunov exponent map

Let us plot the Lyapunov exponents for the Poincaré map of the previous example.
```@example examples
using Dicke, Dicke.ClassicalDicke, Dicke.ClassicalSystems
using Plots
system = ClassicalDicke.ClassicalSystem(ω=0.8, γ=0.8, ω₀=1)
ε = -1.35
function λ(Q,P)
    if minimum_ε_for(system, P=P, p=0, Q=Q) > ε
        return NaN
    end
    point = Point(system, Q=Q, P=P, p=0, ε=ε)
    return lyapunov_exponent(system, t=5000, u₀ = point, tol = 1e-8)
end
resolution = 0.02 #making this smaller will make a smoother plot,
                 #but it may take time!
if !on_github resolution=0.5 end #hide
maxQ = maximum_Q_for_ε(system,ε) 
maxP = maximum_P_for_ε(system,ε) 

Qs = 0.6:resolution:maxQ
Ps=  0:resolution:maxP #we only calculate the top plane, and later
                       #we mirror it.

mat=[λ(Q,P) for P in Ps, Q in Qs]

mat=vcat(mat[end:-1:2,1:end], mat) #mirror in P
Ps=vcat(-Ps[end:-1:2], Ps)

heatmap(Qs, Ps, mat, xlabel="Q", ylabel="P")
savefig("lyapunov_map.svg");nothing #hide
```
![](lyapunov_map.svg)

It looks very similar to the Poincaré map of the previous example! Notice how there are regular (black) and chaotic (yellow) regions.