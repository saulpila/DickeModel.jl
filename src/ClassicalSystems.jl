module ClassicalSystems

export ClassicalSystem,integrate,lyapunov_spectrum,lyapunov_exponent

using DiffEqCallbacks,OrdinaryDiffEq,DiffEqBase
using ProgressMeter
using RecursiveArrayTools
using LinearAlgebra

struct ClassicalSystem
       params
       step
       out_of_bounds
       varnames
end
function nvars(sistema::ClassicalSystem,vars...)

    vs=[]
    d=Dict(zip(sistema.varnames,1:length(sistema.varnames)))
    for v in vars
        push!(vs,d[v])
    end
    return vs
end
eye(n)=Matrix(Diagonal(ones(n)))

"""
```julia
function integrate(sistema::ClassicalSystem;t,
                           u₀,
                           t₀=0.0,
                           tol=1e-12,
                           show_progress=false,
                           save_intermediate_steps=nothing,
                           saveat=Array{Float64,1}(),
                           cb=nothing,
                           get_fundamental_matrix=false,
                           integrator_alg=TsitPap8(),
                           use_big_numbers=false,
                           integate_backwards=false,
                           kargs...)
```
Test ref Ref. [Pilatowsky2020](@cite).

# Arguments:
- `system` is an instance of [`ClassicalSystems.ClassicalSystem`](@ref).
"""
 function integrate(sistema::ClassicalSystem;t,
                            u₀,
                            t₀=0.0,
                            tol=1e-12,
                            show_progress=false,
                            save_intermediate_steps=nothing,
                            saveat=Array{Float64,1}(),
                            cb=nothing,
                            get_fundamental_matrix=false,
                            integrator_alg=TsitPap8(),
                            use_big_numbers=false,
                            integate_backwards=false,
                            kargs...)
        if t<t₀
           error("Use integate_backwards=true instead of negative time")
        end
        if integate_backwards==true
            integate_backwards=-1
        else
            integate_backwards=1
        end


        if save_intermediate_steps==nothing
            save_intermediate_steps = (length(saveat)==0)
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
                sistema.step(du.x[1],u.x[1],p,t)
                sistema.step.jac(cacheJ,u.x[1],p,t)
                mul!(du.x[2], cacheJ, u.x[2])
            end
            sys = ODEFunction(sistemaconmatriz)


        else
            sys=sistema.step
        end

        prog=nothing
        if show_progress
            prog = Progress(Int(ceil(t-t₀)), 1)
            function progfunc(u,nt,integrator)
                update!(prog,  Int(ceil(nt-t₀)))
            end
            cbprog=FunctionCallingCallback(progfunc;func_everystep=true)
            if cb!=nothing
                cb=CallbackSet(cb,cbprog)
            else
                cb=cbprog
            end
        end
        tspan = (t₀,t)
        params=[sistema.params;integate_backwards]
        prob = ODEProblem(sys,u₀,tspan,params)
        s= solve(prob;alg=integrator_alg,isoutofdomain=sistema.out_of_bounds,callback=cb,saveat=saveat,save_everystep=save_intermediate_steps,progress=false,dtmin=tol,force_dtmin=true,abstol=tol,reltol=tol,kargs...)
        if show_progress
            finish!(prog)
        end

        return s
    end


 function lyapunov_spectrum(punto,t)
        #Ψ=big.(punto.x[2])
        Ψ=punto.x[2]
        λs=[convert(Float64,real(log(α)/(t))) for α in LinearAlgebra.svd(Ψ).S]
        return λs
    end
    function lyapunov_exponent(sistema::ClassicalSystem;kargs...)
        return max(lyapunov_spectrum(sistema;kargs...)...)

     end
    function lyapunov_spectrum(sistema::ClassicalSystem;kargs...)
        v=false
        function save(uvar,t,integrator)

            v=(copy(uvar),t)
        end
        cb=FunctionCallingCallback(save;func_everystep=true)
        integrate(sistema;save_intermediate_steps=false,cb=cb,get_fundamental_matrix=true,verbose=false,kargs...)

        return lyapunov_spectrum(v...)
     end
end
