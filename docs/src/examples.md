# Examples
```@setup examples
push!(LOAD_PATH,"../../src")
on_github=get(ENV, "CI", nothing) == "true"
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
using DiffEqBase

system=ClassicalDicke.ClassicalSystem(ω=0.8,γ=0.8,ω₀=1)
mplot=scatter(fmt=:png, key=false, markersize=1, legend=false,
    size=(500,500), color_palette=:darkrainbow, xlabel="Q", ylabel="P") 

pts = Tuple{Float64, Float64}[] #a list of points (Q,P)
function save(state) #this function saves (Q,P) to pts if q = q₊ (and not q₋).
    if q_sign(system,state.u,ε) == + 
        Q,q,P,p = state.u 
        push!(pts, (Q,P))  
    end                     
end
callback=ContinuousCallback((x, t, _) -> x[4], #when p=x[4] is 0,
    save; #execute the function save
    save_positions=(false,false), abstol=1e-3)
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
        if !on_github break end #hide


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
using DiffEqBase

system = ClassicalDicke.ClassicalSystem(ω=0.8, γ=0.8, ω₀=1)
ε = -1.35

resolution = 0.01 #making this smaller will make a smoother plot,
                  #but it may take time!
maxQ = maximum_Q_for_ε(system,ε) 
maxP = maximum_P_for_ε(system,ε) 

Qs = 0.6:resolution:maxQ
Ps = 0.0:resolution:maxP  #we only compute the top half of 
                          #the plane, and later mirror it

#this matrix contains NaNs (outside of bounds) and tuples [λ,n] 
#(inside of bounds), where λ is an average of all the Lyapunov 
#exponents of the n trajectories that have passed through that square
matrix = [if minimum_ε_for(system, P=P, p=0, Q=Q) > ε 
            NaN else [0.0,0] end 
                for Q in Qs, P in Ps] 

pts = Tuple{Float64, Float64}[] #a list of points (Q,P) for temporary storage

function save(state) #this function saves (Q,P) to pts if q = q₊ (and not q₋).
    if q_sign(system,state.u,ε) == + 
        Q,q,P,p = state.u 
        push!(pts, (Q,P))  
    end                     
end
callback=ContinuousCallback((x, t, _) -> x[4], #when p=x[4] is 0,
    save; #execute the function save
    save_positions=(false,false), abstol=1e-3)
#an auxiliary function, which gives the index k so that r[k] ≥ v >r[k]
function index_of_subinterval(v,r) 
    if r[end] <v
        return length(r)
    end
    if v < r[1] 
        return 1
    end
    return Int(round(((v-r[1])/(r[end]-r[1]))*(length(r)-1) +1))
end
#we iterate over the matrix and the values of Q,P
for ((Q,P),element) in zip(Iterators.product(Qs,Ps),matrix)
    #if we are out of bounds or we have information on this Lyapunov exponent,
    if element===NaN || element[2] > 0 
        continue #we skip
    end
    
    #the following 4 lines will populate pts with all the points (mQ,mP) in the 
    #Poincaré map that the trajectory starting at (Q,P) passes through.
    empty!(pts)
    push!(pts,(Q,P))
    point  = Point(system, Q=Q, P=P, p=0, ε=ε)
    λ = lyapunov_exponent(system, t=5000, u₀=point, callback=callback)
    
    
    for (mQ,mP) in pts
        #we save the Lyapunov into all of the squares that the trajectory visited.
        el = matrix[index_of_subinterval(mQ,Qs),index_of_subinterval(abs(mP),Ps)]
        #if we are inside bounds
        if el !== NaN
           
            el[1] = (el[1]*el[2] + λ)/(el[2] + 1)  #we update the average
            el[2] += 1 #we update the count
        end
    end
    if !on_github break end #hide
end

mat=transpose([v[1] for v in matrix]) #we take the average Lyapunovs, and tranpose
#because heatmap takes transposed matrices.

mat=vcat(mat[end:-1:2,1:end], mat) #mirror in P
Ps=vcat(-Ps[end:-1:2], Ps) #update Ps with negative values

heatmap(Qs, Ps, mat, xlabel="Q", ylabel="P",size=(550,500))
savefig("lyapunov_map.svg");nothing #hide
```
![](lyapunov_map.svg)

It looks very similar to the Poincaré map of the previous example! Notice how there are regular (black) and chaotic (colored) regions.