module ClassicalSystems

export integrate,lyapunov_spectrum,lyapunov_exponent,ClassicalSystem

using DiffEqCallbacks,OrdinaryDiffEq,DiffEqBase
using ProgressMeter
using RecursiveArrayTools
using LinearAlgebra

"""
```julia
abstract struct ClassicalSystem
```
This abstract object represents a classical system that may be passed to multiple functions
in this module. To generate a system, use the subtypes [`ClassicalDicke.ClassicalDickeSystem`](@ref Main.ClassicalDicke.ClassicalDickeSystem)
or use [`ClassicalLMG.ClassicalLMGSystem`](@ref Main.ClassicalDicke.ClassicalLMGSystem).
"""
abstract type ClassicalSystem end

function parameters(system::ClassicalSystem)
    error("Please define a method ClassicalSystems.parameters(system::$(typeof(system)) that returns the list
    of parameters for this system")
end
function step(system::ClassicalSystem)
    error("Please define a method ClassicalSystems.step(system::$(typeof(system)) that returns ode system
    of equations")
end
function out_of_bounds(system::ClassicalSystem)
    error("Please define a method ClassicalSystems.out_of_bounds(system::$(typeof(system)) that returns a function 
    out_of_bounds(u,_,t) that determines if u is out of bounds")
end
function varnames(system::ClassicalSystem)
    error("Please define a method ClassicalSystems.varnames(system::$(typeof(system)) that returns a function 
    varnames(u,_,t) that reurns the coordinate names in order")
end
function nvars(system::ClassicalSystem,vars...)

    vs=[]
    d=Dict(zip(varnames(system),1:length(varnames(system))))
    for v in vars
        push!(vs,d[v])
    end
    return vs
end
eye(n)=Matrix(Diagonal(ones(n)))

"""
```julia
function integrate(system::ClassicalSystem;t::Real,
                   u₀::AbstractArray{<:Real, 1},
                   t₀=0.0::Real,
                   tol=1e-12::Real,
                   get_fundamental_matrix::Bool=false,
                   integrator_alg=TsitPap8(),
                   use_big_numbers::Bool=false,
                   integate_backwards::Bool=false,
                   kargs...)
```
This function integrates initial condition `u₀` from `t₀` to `t` under the Hamiltonian system determined by `system`, returning an instance of [`OrdinaryDiffEq.ODESolution`](https://diffeq.sciml.ai/dev/basics/solution/).

# Arguments:
- `system` is an instance of [`ClassicalSystems.ClassicalSystem`](@ref).
- `u₀` is an array which codifies the initial condition `[Q,q,P,p]` for `ClassicalDicke` and `[Q,P]` for `ClassicalLMG`.
- `t` is the start time of the integration.
- `t₀` is the start of the integration (defaults to `t₀ = 0.0`) [`ClassicalSystems.ClassicalSystem`](@ref).
- `tol` is the tolerance for the integration, which determines both `abstol` and `reltol` in [`OrdinaryDiffEq.solve`](https://diffeq.sciml.ai/stable/basics/common_solver_opts/)
- `get_fundamental_matrix` determines whether to also compute the fundametal matrix of the system. If `true`, the result at each time is an [`ArrayPartition(x,Φ)`](https://diffeq.sciml.ai/stable/features/diffeq_arrays/#ArrayPartitions), so that `x=result.x[2]` retrieves the coordinate and `Ψ=result.x[2]` retrieves the fundamental matrix. Default is `false`. Note that the integration is consideribly slowed down if this parameter is set to `true`.
- `integrator_alg` is the integration algorithm to use. Defaults to `TsitPap8` (Tsitouras-Papakostas 8/7 Runge-Kutta method). See [the `DifferentialEquations` documentation](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Full-List-of-Methods) for other options.
- `use_big_numbers` forces the integration to be performed with `BigFloat` instead of `Float`, allowing for infinite numerical precision, but hindering speed substantially. Defaults to `false`.
- `integate_backwards` tells the integrator to integrate back in time, from `-t₀` to `-t`. Defaults to  `false`.
- Additional `kargs` are passed to [`OrdinaryDiffEq.solve`](https://diffeq.sciml.ai/stable/basics/common_solver_opts/).
"""
 function integrate(system::ClassicalSystem;t::Real,
                            u₀::AbstractArray{<:Real, 1},
                            t₀=0.0::Real,
                            tol=1e-12::Real,
                            get_fundamental_matrix::Bool=false,
                            integrator_alg=TsitPap8(),
                            use_big_numbers::Bool=false,
                            integate_backwards::Bool=false,
                            kargs...)
        if t<t₀
           error("Use integate_backwards=true instead of negative time")
        end
        if integate_backwards==true
            integate_backwards=-1
        else
            integate_backwards=1
        end
        if length(u₀)!=length(varnames(system))
            error("u₀ should have length $(length(varnames(system))), which is the dimension of the phase space of system")
        end
        if out_of_bounds(system)(u₀,nothing,0)
            error("The initial condition u₀=$u₀ is out of bounds")
        end

        u₀=u₀
        if use_big_numbers
            u₀=big.(u₀)
        end
        sys=nothing

        if get_fundamental_matrix
            n=length(u₀)
            I=eye(n)
            cacheJ=eye(n)
            if use_big_numbers
                I=big.(I)
                cacheJ=big.(cacheJ)
            end

            u₀ = ArrayPartition(u₀,I)
            function sistemaconmatriz(du,u,p,t)
                step(system)(du.x[1],u.x[1],p,t)
                step(system).jac(cacheJ,u.x[1],p,t)
                mul!(du.x[2], cacheJ, u.x[2])
            end
            sys = ODEFunction(sistemaconmatriz)


        else
            sys=step(system)
        end

        tspan = (t₀,t)
        params=[parameters(system);integate_backwards]
        prob = ODEProblem(sys,u₀,tspan,params)
        s= solve(prob;alg=integrator_alg,isoutofdomain=out_of_bounds(system),progress=false,dtmin=tol,force_dtmin=true,abstol=tol,reltol=tol,kargs...)
        return s
    end

    function lyapunov_spectrum(point,t)
        #Ψ=big.(point.x[2])
        Ψ=point.x[2]
        λs=[convert(Float64,real(log(α)/(t))) for α in LinearAlgebra.svd(Ψ).S]
        return λs
    end

    """
    ```julia
    function lyapunov_exponent(system::ClassicalSystem;kargs...)
    ```
    Returns the maximal [Lyapunov exponent](https://en.wikipedia.org/wiki/Lyapunov_exponent) for `system`.
    # Arguments
    - `system` is an instance of [`ClassicalSystems.ClassicalSystem`](@ref).
    - `kargs...` are redirected to [`ClassicalSystems.integrate`](@ref). In particular, you should pass `u₀`. The value `t` is set to `10000`
      but you may change it: lower values increase speed, higher values increase precision. 
    """
    function lyapunov_exponent(system::ClassicalSystem;kargs...)
        return max(lyapunov_spectrum(system;kargs...)...)

     end
    """
    ```julia
    function lyapunov_spectrum(system::ClassicalSystem;kargs...)
    ```
    Same as [`ClassicalSystems.lyapunov_exponent`](@ref), but returns the whole [Lyapunov spectrum](https://en.wikipedia.org/wiki/Lyapunov_exponent#Definition_of_the_Lyapunov_spectrum).
    """
    function lyapunov_spectrum(system::ClassicalSystem;callback=nothing,t=10000,kargs...)
        v=false
        function save(uvar,t,integrator)
            v=(copy(uvar),t)
        end
        cb=FunctionCallingCallback(save;func_everystep=true)
        if callback===nothing
            callback=cb
        else
            callback=CallbackSet(callback,cb)
        end
        r=integrate(system;save_everystep=false,get_fundamental_matrix=true,callback=callback,t=t,verbose=false,kargs...) 
        return lyapunov_spectrum(v...)
     end
end
