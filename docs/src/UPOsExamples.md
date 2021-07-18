# Examples for UPOs
```@setup examples
push!(LOAD_PATH,"../../src")
on_github=get(ENV, "CI", nothing) == "true"
cache_fold_name="./diags"
use_current_dir_for_diags=on_github
using DickeModel
```
The module [`DickeModel.UPOs`](@ref DickeModel.UPOs) provides a toolbox for studying periodic orbits (POs) 
in the classical limit of the Dicke model [Pilatowsky2021](@cite), and their relation with the ubiquitous quantum
scarring they produce [Pilatowsky2021NatCommun](@cite).


## The fundamental families of periodic orbits.
In Ref. [Pilatowsky2020](@cite), it was shown that two families of periodic orbits
emanate from the normal modes of the ground state configuration of the classical Dicke
model. These families are called family ``\mathcal{A}`` and family ``\mathcal{B}``. They
are implemented in [`UPOs.family_A`](@ref DickeModel.UPOs.family_A) and [`UPOs.family_B`](@ref DickeModel.UPOs.family_B), respectively. Let
us look at family ``\mathcal{A}``.
```@setup examples
@info "Starting example: Family A"
```
```@setup examples
using Plots
plotly() #first call generates an info message that we want to hide
```
```@example examples
using DickeModel
using DickeModel.UPOs
using DickeModel.ClassicalDicke
using Plots
plotly() #interactive plots
Plots.isijulia() = true #hide

system = ClassicalDickeSystem(ω=1.0, γ=1.0, ω₀=1.0)
fam_A = UPOs.family_A(system)
fam_A = UPOs.family_A(system, verbose=false) #hide
nothing; #hide
```
You may retrieve the PO at a given energy ``\epsilon`` by calling `fam_A(ϵ)`,
```@repl examples
po = fam_A(-0.5)
```
Having an instance of [`PO`](@ref DickeModel.UPOs.PO), you may extract many properties 
of the periodic orbit:
```@repl examples
UPOs.lyapunov(po) #Lyapunov exponent 
po.T #period
UPOs.action(po) #action
UPOs.jz_PO_average(po) #average of jz along orbit
```
Let us plot the orbits of ``\mathcal{A}`` from the ground state energy up to ``\epsilon = 0``.
Using [`QPp`](@ref DickeModel.UPOs.QPp), we may extract the coordinates ``(Q,P,p)`` and then simply
pass them to `plot!` to generate a 3D plot. We also draw the mirrored
family ``\widetilde{\mathcal{A}}`` (see  [Pilatowsky2020](@cite)), which we obtain
by applying  [`mirror_Qq`](@ref DickeModel.UPOs.mirror_Qq).
```@example examples
ϵ₀ = minimum_energy(system)

orbits_plot = plot(xlabel="Q", ylabel="P", zlabel="p", size=(800,500))

energies = ϵ₀:0.1:0
for ϵ in energies
    plot!(orbits_plot, QPp(fam_A(ϵ)), 
        label = "ϵ = $(round(ϵ,digits=3))") 
    plot!(orbits_plot, QPp(mirror_Qq(fam_A(ϵ))), 
        label = "ϵ = $(round(ϵ,digits=3))",
        primary=false,linestyle=:dot) #mirrored orbit

end
gr() #hide
orbits_plot
```

```@setup examples
Plots.isijulia() = false 
```
Now try with family ``\mathcal{B}``!
## Finding periodic orbits with the Husimi function of the eigenstates.

In this example, we will find a family of periodic orbits by looking at the scars it
produces in a quantum eigenstates. This procedure was used in Ref. [Pilatowsky2021Identification](@cite) with
 ``j=100``. In this example we will work with a smaller system size of ``j=30``.

First, let us setup our quantum system.
```@setup examples
@info "Starting example: UPOS from Husimi"
```
```@example examples
using DickeModel
using DickeModel.DickeBCE
using DickeModel.ClassicalDicke
using DickeModel.EnergyShellProjections
using Plots

j = 30
Nmax = 120
system = QuantumDickeSystem(ω=1.0, γ=1.0, ω₀=1.0, j = j, Nmax=120)

if false #hide
eigenenergies,eigenstates = diagonalization(system) 
end #hide
if !use_current_dir_for_diags #hide
eigenenergies,eigenstates =  diagonalization(system,verbose=false) #hide
else #hide
eigenenergies,eigenstates =  diagonalization(system, cache_folder=cache_fold_name,verbose=false)  #hide
end #hide
ϵs = eigenenergies/j
nothing; #hide
```

If you want to speed up the computation of the Husimi projections, you may load more modules:

```julia
using Distributed
addprocs(4) #you may add as many as processors in your computer
@everywhere using DickeModel
```
In Ref. [Pilatowsky2021Identification](@cite) it was found
that the following procedure yields periodic orbits for eigenstates with a low value of the Rényi occupation. Having looked at
the occupations for many eigenstates beforehand, we selected `k = 559`. You can use [`EnergyShellProjections.rényi_occupation`](@ref DickeModel.EnergyShellProjections.rényi_occupation) to compute the  Rényi occupation of many eigenstates and then select
those with low values yourself (see Ref. [Pilatowsky2021Identification](@cite) for more details.)

Let's take a look at its Husimi projection
using [`EnergyShellProjections.proj_husimi_QP_matrix`](@ref DickeModel.EnergyShellProjections.proj_husimi_QP_matrix).
```@setup examples
@info "Starting example: UPOS from Husimi, plot"
```
```@example examples
k = 559
state = @view eigenstates[:,k]
res = 0.04
if !on_github #hide
    res = 0.1 #hide
end #hide
Qs,Ps,mat = proj_husimi_QP_matrix(system, 
    state, 
    ϵ = ϵs[k], 
    res = res, 
    show_progress = false, #hide
    symmetricQP = true, 
    chop = 1e-2)
heatmap(Qs, Ps, mat, xlabel="Q", ylabel="P")
```

It has a nice scar-like shape. Let us find the periodic orbit that is generating 
this scar. To do this, we will increase the moment of the Husimi function, allowing
us to better identify the scar. Also, the projection above is a sum of both the positive and negative q sides of the
energy shell, but it is easier to stay only on one side, so let us stay only the positive q solution, loosing the symmetry
in ``Q``. We add arguments `matrix_powers`, `onlyqroot`, and `symmetricP` (not `symmetricQP`):
```@setup examples
@info "Starting example: UPOS from Husimi, half plot"
```
```@example examples
Qs₊, Ps₊, mat₊ = proj_husimi_QP_matrix(system, 
    state, 
    ϵ = ϵs[k],
    show_progress = false, #hide
    matrix_powers=2, 
    res = res, 
    symmetricP = true, 
    chop = 1e-2,
    onlyqroot = +)

heatmap(Qs₊, Ps₊, mat₊, xlabel="Q", ylabel="P")
```

Okay, the procedure is to identify the peaks in this matrix. Let us define a function
that finds peaks in a matrix:
```@setup examples
@info "Starting example: UPOS from Husimi, peaks"
```
```@example examples
function find_2D_peaks(mat)
    pks=[] #a list to save the peaks
    ny,nx = size(mat) #dimentions of the matrix
    for x in 2:(nx-1), y in 2:(ny-1) #for each internal point
            if isnan(mat[y,x]) #if is outside  we skip
                continue
            end
            #If any of the neighboring coordinates has a
            #greater value, we skip
            if any(!isnan(mat[y+j,x+i]) && mat[y,x]<mat[y+j,x+i] 
                for (i,j) in [(1,0),(0,1),(-1,0),(0,-1)])
                continue
            end
            #if we got to here, habemus peak: 
            push!(pks,[x,y])
    end
    #We sort the peaks by their intensity
    pks = sort!(pks, by = a ->  mat[a[2],a[1]], rev=true)
    return pks
end
nothing; #hide
```
Lets apply this function to `mat₊`:

```@example examples
#we take only the first 19 peaks: there are many small ones!
peaks = find_2D_peaks(mat₊)[1:19] 
peaks_plot = heatmap(Qs₊, Ps₊, mat₊, 
    key=false, xlabel="Q", 
    ylabel="P", size=(700,500))
    
for peak in peaks
    Q,P = Qs₊[peak[1]], Ps₊[peak[2]]
    scatter!([(Q,P)], 
        series_annotations = text("$((round(Q,digits=5),round(P,digits=5)))  ", 
                                    :white, :right, 7)
        )
end
peaks_plot
```
We have many peaks. They idea is to try to use them as initial conditions for our
periodic orbit. We have to choose one peak, and identify the coordinates for  ``p``.
There are many choices, and you will have to do some trial and error.
Let us select the peak at ``Q=0.08, P=0.6``. This one will work. 

To identify the coordinate of the peak in  ``p``, we will find the local maxima of the
Husimi function with varying  ``p``. To do this, we will use  [`Optim.jl`](https://julianlsolvers.github.io/Optim.jl/stable/)
```@setup examples
@info "Starting example: UPOS from Husimi, optim"
```
```@example examples
using Optim
Q0=0.08
P0=0.6

#this function gives the Husimi function of the state evaluated at the
#point [Q,q,P,p], where q is on the + side of the energy shell. If
#the point is outside the energy shell, it returns 0.
function hus(;p,Q,P)
    let point
        try 
            point=Point(system, Q=Q, P=P, p=p, ϵ=ϵs[k], signo=+)
        catch #the point is outside the energy shell.
            return 0
        end

        return DickeBCE.husimi(system,point,state,chop=1e-3)
    end
end

#let's vary the p coordinate
h(p) = hus(p=p, Q=Q0, P=P0)
plot(h, -2, 2, label="Husimi(Q0,P0,p,q₊)", xlabel="p")
#we find a local maxima between -2 and 2 (you may move this to 
#    identify the peaks you want)
pmax = Optim.maximizer(maximize(h,-2,2))

plot!([pmax,pmax], [0,h(pmax)], label = "pmax=$pmax")
```

Okay, we now have our coordinates ``Q,P,p``. Before using them as initial conditions, 
recall that we got ``Q,P`` from our projection, which did not have the best of the
resolutions. We can refine them a little better by maximizing the Husimi function varying
both ``P`` and ``p``:
```@repl examples
pmax1, P1 = Optim.maximizer(maximize(x -> hus(p=x[1], P=x[2], Q=Q0), [pmax,P0]))
```
We now can define an initial condition
```@repl examples
u0 = Point(system, Q=Q0, P=P1, p=pmax1, ϵ=ϵs[k])
```

It is time to use the module [`UPOs`](@ref DickeModel.UPOs). Let us first see how the evolution of this
point looks like with [`UPOs.QP`](@ref DickeModel.UPOs.QP). Let's evolve it for, say, `T = 10`:
```@setup examples
@info "Starting example: UPOS from Husimi, evolution test"
```
```@example examples
using DickeModel.UPOs
matrixPlot = heatmap(Qs, Ps, mat, xlabel="Q", ylabel="P", key=false)
scatter(matrixPlot, [(u0[1],u0[3])])
plot!(QP(PO(system, u0, 10)))
```

It looks like it is following the scar. To estimate the period, we can use [`UPOs.approximate_period`](@ref DickeModel.UPOs.approximate_period).
We give it a `bound = 0.5`, and it will tell us the time it takes for the periodic condition to come back
to the same `p` plane inside a neighborhood of radius `bound` around `u0`:
```@repl examples
T = UPOs.approximate_period(system, u0; bound=0.5)
```
It found a period! The evolution up to that period looks like this:
```@example examples
plot(matrixPlot, QP(PO(system, u0, T)))
```
Well, it doesn't come back completely, and it looks a bit wonky, but the monodromy method [DeAguiar1988](@cite), [Baranger1988](@cite), [Pilatowsky2021](@cite), [Simonovi1999](@cite) will come to our aid. This algorithm converges to a true periodic orbit given a good enough initial guess, like ours.
It is implemented in [`UPOs.monodromy_method_constant_energy`](@ref DickeModel.UPOs.monodromy_method_constant_energy), which conserves the energy
of the initial condition.
```@setup examples
@info "Starting example: UPOS from Husimi, monodromy"
```
```@repl examples
po = monodromy_method_constant_energy(system, u0, T)
```
```@example examples
plot(matrixPlot, QP(po))
```
That is a nice periodic orbit!

We can measure that it is actually scarring the eigenstate using the measure ``\mathcal{P}`` from Ref. [Pilatowsky2021](@cite), which
is implemented in [`UPOs.scarring_measure`](@ref DickeModel.UPOs.scarring_measure):

```@repl examples
scarring_measure(po,state)
```
This tells us that the state is ``\sim 12`` times more likely to be found near the PO than
a totally delocalized state.

We can also compute its Lyapunov exponent using [`UPOs.lyapunov`](@ref DickeModel.UPOs.lyapunov)
```@repl examples
lyapunov(po)
```
Interestingly, this orbit is stable, with a zero Lyapunov exponent. However, the eigenenergy of this
eigenstate, ``\epsilon_k = -0.806`` is in the chaotic region, so we found a little island of stability in the sea of chaos.

Moreover, this orbit is symmetric with respect to the change of parity ``(Q,q)\mapsto (-Q,-p)``:
```@repl examples
mirror_Qq(po) == po
```
!!! note 
    The method [`po1 == po2`](@ref ==(::DickeModel.UPOs.PO,::DickeModel.UPOs.PO)) actually checks whether the orbits are equal by integration. The comparison
    `po1 == po2` is `true` if the orbits are the same, even if `po1.u != po2.u`.

We close this example using the function [`UPOs.follow_PO_family_from_energy`](@ref DickeModel.UPOs.follow_PO_family_from_energy) to perturb
this orbit into higher energies, obtaining a family of POs.
```@setup examples
@info "Starting example: UPOS from Husimi, family plot"
```
```@example examples
plotly() #interactive plots
Plots.isijulia() = true #hide
family = follow_PO_family_from_energy(po)
orbits_plot = plot(xlabel="Q", ylabel="P", zlabel="p", size=(800,600))
energies = -0.8:0.03:0
for ϵ in energies
    plot!(orbits_plot, QPp(family(ϵ)), label = "ϵ = $ϵ") 
end
gr() #hide
orbits_plot
```
```@setup examples
Plots.isijulia() = false 
```
This family of UPOs seems to contain the UPO ``\mathcal{O}_\text{B}`` from Ref. [Pilatowsky2021Identification](@cite). See Fig 3. (B2).

We can see how the stability of this family changes as a function of the energy (like Fig. 1 of Ref. [Pilatowsky2021](@cite)):
```@setup examples
@info "Starting example: UPOS from Husimi, family Lyapunov"
```
```@example examples
energies = -0.8:0.01:0
plot(energies, lyapunov.(family.(energies)),
    xlabel="energy", ylabel="Lyapunov exponent", key = false)
```
