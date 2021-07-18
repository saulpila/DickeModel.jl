module UPOs

export PO,search_in_interval,approximate_period,
        get_period,follow_PO_family_from_period,
        find_p_zero,integrate,monodromy_method_constant_period,
        monodromy_method_constant_energy,po_coordinates,QP,qp,QPp,QPq,mirror_Qq,
        mirror_Pp,scarring_measure,lyapunov,follow_PO_family_from_energy,energy,
        family_A,family_B
    using ..ClassicalSystems
    using ..ClassicalDicke
    using LinearAlgebra
    using DiffEqCallbacks,OrdinaryDiffEq, DiffEqBase
    using ProgressMeter
    using ..PhaseSpaces
    using ..DickeBCE
    using ..EnergyShellProjections
    Id=Matrix{Float64}(I, 4, 4);
    """
    ```julia
    struct PO
    ```
    This struct represents a periodic orbit in the classical Dicke model. To generate,
    use 
    ```julia
    po = PO(system, u, T)
    ```
    where `system` is an instance of [`ClassicalDicke.DickeSystem`](@ref),
    `u` is a real vector of an initial condition in the form `[Q,q,P,p]`, and `T`
    is a real number representing the period. You can retrieve this values using
    `po.system`, `po.u`, and `po.T`.
    """
    struct PO{System <: DickeSystem}
        system::System
        u::Vector{Float64}
        T::Real
        PO{System}(system::System, u::AbstractVector{<:Real},
            T::Real) where {System<:DickeSystem} = new(system,Float64.(u),T)
    end
    PO(system::System, u::AbstractVector{<:Real},
        T::Real) where {System<:DickeSystem} = PO{System}(system,u,T)
        
    function Base.show(io::IO, m::PO)
         print(io,string("PO @ u=",m.u,", T=", m.T))
    end
    """
    ```julia
    Base.:(==)(po1::PO,po2::PO;Ttol::Real=1e-5,tol::Real=1e-6)
    ```
    The comparison `po1 == po2` returns `true` if `po1` and `po2` represent the same periodic orbit, regardless
    of which initial condition they have.
    This is done by evaluating `∈(po1.u,  po2, tol = tol) && abs(po1.T-po2.T) < Ttol`
    """
    function Base.:(==)(po1::PO,po2::PO;Ttol::Real=1e-5,tol::Real=1e-6)
        if abs(po1.T-po2.T)>= Ttol
            return false
        end
        return ∈(po1.u,  po2,tol = tol)
    end
    """
    ```julia
    Base.in(u::AbstractVector{<:Number}, po::PO;tol::Real=1e-6)
    Base.∈(u::AbstractVector{<:Number}, po::PO;tol::Real=1e-6)
    ```
    
    You may write `u ∈ po`.
    
    Returns `true` if the point `u` is part of the periodic orbit `po`. This
    is done by integrating `po.u` and seeing if it comes within `tol` of `u`.
    """
    function Base.in(u::AbstractArray{<:Number,1}, po::PO;tol::Real=1e-6)
        if sum(abs2,u-po.u)<tol
            return true
        end

        result=false

        function distance(uvar,t,integrator)
            m=sum(abs2,uvar-u)
            if m<tol
                result=true
                terminate!(integrator)
            end
            return m*sign(uvar[4]-u[4])
        end
        cb=ContinuousCallback(distance,integrator->distance(integrator.u,integrator.t,integrator),save_positions=(false,false),rootfind=true,interp_points=10,reltol=10^-8,abstol=10^-8)
        ClassicalSystems.integrate(po;tol=1e-8,save_everystep=false,callback=cb)
        return result
    end
    """
    ```julia
    function PO(system::DickeSystem, u::AbstractVector{<:Real})
    ```
    Generates a [`PO`](@ref), given a periodic condition `u = [Q,q,P,p]`, where 
    the period is calculated using [`approximate_period(system,u,bound=1e-2)`](@ref).
    """
    PO(system::DickeSystem,u::AbstractVector{<:Real})=PO(system,u,approximate_period(system,u,bound=1e-2))
    
    """
    ```julia
    function integrate(po::PO;tol::Real=1e-16, kargs...)
    ```
    Returns an instance of [`OrdinaryDiffEq.ODESolution`](https://diffeq.sciml.ai/dev/basics/solution/), resulting
    from integrating `po.u` from `t = 0` to `t = po.T`. 
    
    # Arguments
    - `po` should be an instance of [`PO`](@ref).
    # Keyword arguments
    - `tol` is a real number indicating the precision for the integrator. Defaults to `1e-16`.
    - `kargs` are redirected to [`ClassicalSystems.integrate`](@ref ClassicalSystems.integrate(::ClassicalSystems.ClassicalSystem,::AbstractVector{<:Real},::Real)).
    """
    ClassicalSystems.integrate(po::PO;tol::Real=1e-16, kargs...)=ClassicalSystems.integrate(po.system,po.u,po.T;tol=tol, kargs...);
    
    """
    ```julia
    function action(po::PO;kargs...)
    ```
    Numerically computes the [action](https://en.wikipedia.org/wiki/Action_(physics)#Abbreviated_action_(functional))
    ```math
        S = \\int p\\,\\text{d} q + P \\, \\text{d} Q
    ```
    along the periodic orbit `po`.
    # Arguments
    - `po` should be an instance of [`PO`](@ref)
    # Keyword arguments
    - `kargs` are redirected to [`integrate`](@ref integrate(::ClassicalSystem,::AbstractVector{<:Real},::Real))
    """
    function action(po::PO;kargs...)
        us=integrate(po;kargs...).u
        return sum(us[i][3]*(us[i][1]-us[i-1][1]) + us[i][4]*(us[i][2]-us[i-1][2]) for i in 2:length(us))
    end
    """
    ```julia
    function average_over_PO(po::Union{OrdinaryDiffEq.ODESolution,PO},
        f::Function;
        kargs...)
    ```
    Computes the average of a function ``f`` over the periodic orbit `po`
    ```math
        \\frac{1}{T} \\int_0^T f(u(t)) \\, \\text{d} t,
    ```
    where ``T=`` `po.T` and ``u=`` `po.u`.
    # Arguments
    - `po` should be an instance of [`PO`](@ref) or an instance of [`OrdinaryDiffEq.ODESolution`](https://diffeq.sciml.ai/dev/basics/solution/),
      in which case the integration is ommited.
    - `f` should be a function with a method `f([Q,q,P,p])`, which returns values that
      may be added together and multiplied by scalars (e.g. numbers or arrays).
    # Keyword arguments
    - `kargs` are redirected to  [`ClassicalSystems.integrate`](@ref ClassicalSystems.integrate(::ClassicalSystems.ClassicalSystem,::AbstractVector{<:Real},::Real)).
    """
    function average_over_PO(po::Union{OrdinaryDiffEq.ODESolution,PO},
        f::Function;
        tol::Real=1e-12,
        kargs...)
        
        tot=nothing

        if isa(po,PO)
            function guardar(u,t,integrator)
                v=f(u)*(t-integrator.tprev)/po.T
                if tot===nothing
                    tot=v
                else
                    tot+=v
                end
            end
            cb=FunctionCallingCallback(guardar;func_start = false,func_everystep =true)
            ClassicalSystems.integrate(po;kargs...,tol=tol,save_everystep=false,callback=cb)
        else
            T=po.t[end]-po.t[1]
            for i in 2:length(po.t)
                v=f(po.u[i])*(po.t[i]-po.t[i-1])/T
                if tot===nothing
                    tot=v
                else
                    tot+=v
                end
            end
        end
        return tot
    end
    """
    ```julia
    function average_over_PO_QP(po::Union{OrdinaryDiffEq.ODESolution,PO},
        f::Function;
        kargs...)
    ```
    Same as [`average_over_PO`](@ref), but `f` only depends on `Q,P`, that is, the signature
    of `f` is `f(Q,P)`.
    """
    average_over_PO_QP(PO::PO,f::Function;kargs...)=average_over_PO(PO,x->f(x[1],x[3]);kargs...)
    """
    ```julia
    function average_over_PO_qp(po::Union{OrdinaryDiffEq.ODESolution,PO},
        f::Function;
        kargs...)
    ```
    Same as [`average_over_PO`](@ref), but `f` only depends on `q,p`, that is, the signature
    of `f` is `f(q,p)`.
    """
    average_over_PO_qp(PO::PO,f::Function;kargs...)=average_over_PO(PO,x->f(x[2],x[4]);kargs...)

    """
    ```julia
    function jz_PO_average(po::Union{OrdinaryDiffEq.ODESolution,PO},
        f::Function;
        kargs...)
    ```
    Returns `average_over_PO_QP(PO,PhaseSpaces.jz;kargs...)`.
    """
    jz_PO_average(PO::PO;kargs...)=average_over_PO_QP(PO,PhaseSpaces.jz;kargs...)
    
    """
    ```julia
    function jy_PO_average(po::Union{OrdinaryDiffEq.ODESolution,PO},
        f::Function;
        kargs...)
    ```
    Returns `average_over_PO_QP(PO,PhaseSpaces.jy;kargs...)`.
    """
    jy_PO_average(PO::PO;kargs...)=average_over_PO_QP(PO,PhaseSpaces.jy;kargs...)
    
    """
    ```julia
    function jx_PO_average(po::Union{OrdinaryDiffEq.ODESolution,PO},
        f::Function;
        kargs...)
    ```
    Returns `average_over_PO_QP(PO,PhaseSpaces.jx;kargs...)`.
    """
    jx_PO_average(PO::PO;kargs...)=average_over_PO_QP(PO,PhaseSpaces.jx;kargs...)
    
    """
    ```julia
    function mirror_Qq(po::PO)
    ```
    Returns a periodic orbit `po1`, which results from changing the sign
    of the ``Q`` and ``q`` coordinates of the periodic orbit `po`.
    """
    mirror_Qq(po::PO)=PO(po.system,po.u.*[-1,-1,1,1],po.T)
    
    """
    ```julia
    function mirror_Pp(po::PO)
    ```
    Returns a periodic orbit `po1`, which results from changing the sign
    of the ``P`` and ``p`` coordinates of the periodic orbit `po`.
    """
    mirror_Pp(po::PO)=PO(po.system,po.u.*[1,1,-1,-1],po.T)
    
    """
    ```julia
    function approximate_period(system::DickeSystem,
        u₀::AbstractVector{<:Real};
        bound::Real = 0.1,
        tol::Real = 1e-8,
        max_p_crossings::Real = Inf,
        verbose::Bool = true)
    ```
    Computes the time necessary for the initial condition `u₀` to come back to the
    same `p`-plane within a neighborhood of radius `bound` of `u₀`. This is useful
    to approximate the period of an *almost* periodic condition `u₀`. The maximum time 
    of integration is `10000`.
    # Arguments
    - `system` should be an instance of [`ClassicalDicke.DickeSystem`](@ref).
    - `u₀` is a point `[Q,q,P,p]` in the phase space of the Dicke model.
    # Keyword arguments
    - `bound` is a positive real number indicating how close does the evolution has 
      to come back to `u₀`. Defaults to `0.1`.
    - `tol` is the integration tolerance. Not to be confused with `bound`. Defaults to  1e-8.
    - `max_p_crossings` is the maximum number of times the condition may cross the `p` plane
      before aborting. Defaults to `Inf` (no maximum).
    - `verbose` is a boolean indicating whether to log information messages. Defaults to `true`.
    """
    function approximate_period(system::DickeSystem,
        u₀::AbstractVector{<:Real};
        bound::Real = 0.1,
        tol::Real = 1e-8,
        max_p_crossings::Real = Inf,
        verbose::Bool = true)

        notreturnedcount=0
        lastsavet=0
        function guardar(integrator)
            c=integrator.u[:]

            returned=bound > norm(c-u₀)
            if !returned
                notreturnedcount+=1
                if notreturnedcount>=max_p_crossings
                    terminate!(integrator)
                    if verbose 
                        @info "Exceeded max_p_crossings = $max_p_crossings"
                    end
                    return nothing
                end
            else
                lastsavet = integrator.t
                terminate!(integrator)
                return nothing
                    
            end
            u_modified!(integrator,false)
            return nothing
        end
        cb=ContinuousCallback((uvar,t,integrator)-> (uvar[4]-u₀[4]),nothing,guardar,save_positions=(false,false),rootfind=true,interp_points=3,reltol=10^-8,abstol=10^-8)
        t=ClassicalSystems.integrate(system,u₀,10000;callback=cb,save_everystep=false,tol=tol).t[end]
        return lastsavet

    end

    function monodromy_method_step_constant_period(system::DickeSystem,
        u₀::AbstractVector{<:Real},
        T::Real;
        tol::Real=1e-12)
        
        u₁,M=ClassicalSystems.integrate(system,u₀,T;save_everystep=false,get_fundamental_matrix=true,tol=tol).u[end].x
        return u₀-(M-Id)^-1*(u₁-u₀),norm(u₁-u₀)
    end
    function hamiltonian_gradient(system::DickeSystem,u::AbstractVector{<:Real})
        Fu₁=Float64[0,0,0,0]
        ClassicalSystems.step(system)(Fu₁,u,[ClassicalSystems.parameters(system);1.0],1.0)
        gradH=[-Fu₁[3],-Fu₁[4],Fu₁[1],Fu₁[2]]
        return gradH
    end

    function monodromy_method_step_constant_energy(system::DickeSystem,u₀::AbstractVector{<:Real},T::Real;tol::Real=1e-12)
        u₁,M=ClassicalSystems.integrate(system,u₀,T;save_everystep=false,get_fundamental_matrix=true,tol=tol).u[end].x
        gradH=hamiltonian_gradient(system,u₁)
        ξ=Float64[0,0,0,1]
        Fu₁=-[-gradH[3],-gradH[4],gradH[1],gradH[2]]
        Mat=[M-Id Fu₁;
            transpose(gradH) 0;
            transpose(ξ) 0]
        res=pinv(Mat)*[u₀-u₁;0;0]
        du=res[1:4]
        dT=res[5]
        return u₀+du,T+dT,norm(u₁-u₀)
    end


    """
    ```julia
    function monodromy_method_constant_period(system::DickeSystem,
        u₀::AbstractVector{<:Real},
        T::Real;
        maxiters::Integer=100,
        tol::Real=1e-8,
        inttol::Real=tol/10)
    ```
    Returns a [`PO`](@ref) with a given period `T` that is found by iteratively perturbating the initial condition `u₀` in
    the direction that minimizes the distance between ``u_0`` and ``u_0(T)``. The perturbation
    is obtained from Eq. (31) of Ref. [Simonovi1999](@cite), taking ``\\Lambda`` as the fundamental
    matrix evaluated at time `T`. Note that the energy of `u₀` will change. To conserve the energy and
    vary `T`, use [`monodromy_method_constant_energy`](@ref) instead.
    
    # Arguments
    - `system` should be an instance of [`ClassicalDicke.DickeSystem`](@ref).
    - `u₀` is a point `[Q,q,P,p]` in the phase space of the Dicke model, which is used 
      as the starting point to find an orbit.
    - `T` is a positive real number, indicating the desired period.
    # Keyword arguments
    - `maxiters` is a integer which sets the maximum number of iterations. Defaults
      to `100`.
    - `tol` is a tolerance. If  ``\\|u-u(T)\\|<`` `tol`,
      ``u`` is considered a periodic condition with period `T`. The smaller the tolerance, the
      more iterations are needed to converge. Default is `1e-8`
    - `inttol` is the tolerance to be passed to [`ClassicalSystems.integrate`](@ref ClassicalSystems.integrate(::ClassicalSystems.ClassicalSystem,::AbstractVector{<:Real},::Real)). Default is
      `tol/10`.
    """
    function monodromy_method_constant_period(system::DickeSystem,
        u₀::AbstractVector{<:Real},
        T::Real;
        maxiters::Integer=100,
        tol::Real=1e-8,
        inttol::Real=tol/10)
    u₁=u₀
    for j in 1:maxiters
        u₁,Δ=monodromy_method_step_constant_period(system,u₁,T,tol=inttol)
        if Δ<tol
            break
        end
        if  j==maxiters
            @error("Maximum iterations reached")
        end
    end

    PO(system,u₁,T)
    end
    """
    ```julia
    function monodromy_method_constant_energy(system::DickeSystem,
        u₀::AbstractVector{<:Real},
        T::Real;
        maxiters::Integer=100,
        tol::Real=1e-8,
        inttol::Real=tol/10,
        correct_energy::Bool=true)
    ```
    Returns a [`PO`](@ref) with energy approximatelly (see below) the same energy as ``u₀`` obtained 
    by iteratively perturbating period `T` and the initial condition `u₀`  in
    the direction that minimizes the distance between ``u_0`` and ``u_0(T)`` and is perpendicular to the
    Hamiltonian gradient. The algorithm is described in App. A, section A.1. of Ref. [Pilatowsky2021](@cite).
    
    # Arguments
    - `system` should be an instance of [`ClassicalDicke.DickeSystem`](@ref).
    - `u₀` is a point `[Q,q,P,p]` in the phase space of the Dicke model, which is used 
      as the starting point to find an orbit.
    - `T` is a positive real number, which is used as the starting point to find an orbit.
    # Keyword arguments
    - `maxiters` is a integer which sets the maximum number of iterations. Defaults
      to `100`.
    - `tol` is a tolerance. If  ``\\|u-u(T)\\|<`` `tol`,
      ``u`` is considered a periodic condition with period `T`. The smaller the tolerance, the
      more iterations are needed to converge. Default is `1e-8`
    - `inttol` is the tolerance to be passed to [`ClassicalSystems.integrate`](@ref ClassicalSystems.integrate(::ClassicalSystems.ClassicalSystem,::AbstractVector{<:Real},::Real)). Default is
      `tol/10`.
    - `correct_energy` is a Boolean. As described in  Ref. [Pilatowsky2021](@cite), this algorithm
      only approximately conserves energy. If `correct_energy == true`, the initial condition 
      is projected back to the energy shell between each iteration. This allows energy to be truly
      conserved, however, it makes the algorithm more unstable. The default is `true`. 
    """
    function monodromy_method_constant_energy(system::DickeSystem,
        u₀::AbstractVector{<:Real},
        T::Real;
        maxiters::Integer=100,
        tol::Real=1e-8,
        inttol::Real=tol/10, 
        correct_energy::Bool=true)
    u₁=u₀
    T₁=T
    ϵ = ClassicalDicke.hamiltonian(system)(u₀)
    for j in 1:maxiters
        if correct_energy
            sgn = ClassicalDicke.q_sign(system,u₁)
            u₁=Point(system,Q=u₁[1],P=u₁[3],p=u₁[4],sgn=sgn,ϵ=ϵ)
        end
        u₁,T₁,Δ=monodromy_method_step_constant_energy(system,u₁,T₁;tol=inttol)
        if T₁<0
            error("Convergence error")
        end
        if Δ<tol
            break
        end
        if j==maxiters
            error("Maximum iterations")
        end
    end

    PO(system,u₁,T₁)
    end
    
    """
    ```julia
    function find_p_zero(system::DickeSystem,
        u₀::AbstractVector{<:Real};
        tol::Real=1e-6,
        negative::Bool=false)
    ```
    Evolves `u₀` until it crosses the ``p=0`` plane from negative to positive and returns the result.
    # Arguments
    - `system` should be an instance of [`ClassicalDicke.DickeSystem`](@ref).
    - `u₀` is a point `[Q,q,P,p]` in the phase space of the Dicke model.
    # Keyword arguments
    - `tol` is the numerical tolerance.
    - If `negative` is `true`, then the crossing is from positive to negative. Default is `false`.
    """
    function find_p_zero(system::DickeSystem,
        u₀::AbstractVector{<:Real};
        tol::Real=1e-6,
        negative::Bool=false)
        nu=Float64[0,0,0,0]

        function guardar(integrator)

            nu=integrator.u

            terminate!(integrator)
            return nothing
        end
        try

            negguardar=nothing
            if negative
                negguardar,guardar=guardar,negguardar
            end
            cb=ContinuousCallback((uvar,t,integrator)-> (uvar[4]-0),guardar,negguardar,save_positions=(false,false),rootfind=true,interp_points=3,abstol=tol)
            ClassicalSystems.integrate(system,u₀,10000;callback=cb,save_everystep=false,tol=tol)
            return nu
        catch
            if abs(u₀[4])< tol  #hay un error cuando el punto ya satisface la condicion por juliaDiff
                return u₀
            end
            return NaN
        end
    end
    """
    ```julia
    function lyapunov(po::PO)
    ```
    Returns the Lyapunov exponent of the periodic orbit `po`, as given by Eq. (B3) of
    Ref. [Pilatowsky2021](@cite).
    """
     function lyapunov(po::PO)
         system=po.system
         T=po.T
         maximum(log.(abs.(
            LinearAlgebra.eigvals(ClassicalSystems.integrate(po;
                    get_fundamental_matrix=true,save_everystep=false)[end].x[2])
            ))/T)
    end
    """
    ```julia
    energy(po::PO) = ClassicalDicke.hamiltonian(po.system)(po.u)
    ```
    Returns the energy of the periodic orbit `po`.
    
    """
    energy(po::PO) =ClassicalDicke.hamiltonian(po.system)(po.u)
    """
    ```julia
    function follow_PO_family_from_period(po::PO;
        step::Real=0.1,
        tol::Real=1e-5,
        initalperturbation::AbstractVector{<:Real}=Float64[0,1,0,0],
        verbose::Bool =true)
    ```
    Returns a function `T -> po1` that returns a PO from the same family as `po`
    but with period `T`. This algorithm applies the function [`monodromy_method_constant_period`](@ref)
    repeatedly, increasing or decreasing the period in small perturbations to reach the target.
    
    # Arguments
    - `po` should be an instance of [`PO`](@ref)
    # Keyword arguments
    - `step` is the initial size of the perturbations in period, this is decreased and increased
      dynamically (default is `0.1`).
    - `tol` is the tolerance to pass to [`monodromy_method_constant_period`](@ref) (default is `1e-5`).
    - `verbose` is a Boolean indicating whether to print the progress. Default is `true`.
    
    The function returned accepts the following keyword arguments:
    - `tol` overrides `tol` above.
    - `minstep` is the minimum step (which is varied dynamically). Default is `step/1000`.
    - `maxstep` is the maximum step. Default is `step`.
    - `step` overrides `step` above. Default is `maxstep`.
    - `maxiters` is the maximum number of iterations. Default is `1000`.
    """
    function follow_PO_family_from_period(po::PO;
        step::Real=0.1,
        tol::Real=1e-5,
        initalperturbation::AbstractVector{<:Real}=Float64[0,1,0,0],
        verbose::Bool =true)
        
        system=po.system
        u₀=po.u 
        period = po.T
        us=[u₀]
    

        periods=Float64[period]
    
    
        function get(finalperiod::Real;tol::Real=tol, minstep::Real=step/1000,maxstep::Real=step,
            step::Real=maxstep,maxiters::Integer=1000)
            index =NaN
            if verbose
                prog=ProgressUnknown("Remaining T")
            end
            _it=0
            while true
                _it+=1
                if _it>maxiters
                    error("Maximum iterations")
                end

                index=findmin(abs.(finalperiod.-periods))[2]
                last_p=periods[index]
                falta=abs(last_p-finalperiod)

                if falta==0
                    break
                end
                step=min(step,falta)
                direction=sign(finalperiod-last_p)
                newu=us[index]
                if verbose
                    prog.desc="Current period = $(round(last_p,digits=3)), step = 10^$(round(log10(step),digits=1)), iterations ="
                    next!(prog)

                end
               try
                    ΔT=direction*step
                    o= monodromy_method_constant_period(system,newu,last_p+ΔT,inttol=min(1e-12,tol/100),tol=tol)
                    
                    newu = o.u
                    period = o.T
                catch e
                    if step/1.2<minstep# (!converror && H(newu) in energies) ||

                            error("The algorithm is not converging. Got to ",PO(system,us[index],periods[index]))
                        
                    end
                    step/=1.2
                    continue
                end

                step=min(maxstep,step*1.1)
    
                nindex=searchsortedlast(periods,period) 
                insert!(us,nindex + 1,newu)
                insert!(periods,nindex + 1,period)
                index=nindex+1
    
            end
            if verbose
                finish!(prog)
            end
            return PO(system,us[index],periods[index])
        end
        return get
    end
    

    """
    ```julia
    function follow_PO_family_from_energy(po::PO;
        step::Real=0.1,
        tol::Real=1e-5,
        energytol::Real=1e-3,
        initalperturbation::AbstractVector{<:Real}=Float64[0,1,0,0],
        verbose::Bool =true,
        dontallowlessenergy::Bool=false)
    ```
    Returns a function `ϵ -> po1` that returns a PO from the same family as `po`
    but with energy `ϵ`. This algorithm is detailed in App. A.2. of Ref. [Pilatowsky2021](@cite)
    
    # Arguments
    - `po` should be an instance of [`PO`](@ref).
    # Keyword arguments
    - `step` is the initial size of the perturbations in energy, this is decreased and increased
      dynamically (default is `0.1`). If the orbits seem to suddenly jump, try to decrease this number.
    - `tol` is the tolerance to pass to [`monodromy_method_constant_energy`](@ref) (default is `1e-5`).
    - `energytol` is the precision in energy below which a PO is returned   (default is `1e-3`).
    - `initalperturbation` is de direction of the first perturbation in the form `[Q,q,P,p]` (default is `[0,1,0,0]`).
    - `verbose` is a Boolean indicating whether to print the progress. Default is `true`.
    - `dontallowlessenergy` forbids evaluation below the energy of `po`. Default is `false`.
    
    The function returned accepts the following keyword arguments:
    - `tol` overrides  `tol` above.
    - `energytol` overrides  `energytol` above.
    - `maxstep` is the maximum step. Default is `step`.
    - `step` overrides `step` above. Default is `maxstep`.
    - `minstep` is the minimum step (which is varied dynamically). Default is `min(step,energytol/10)`.
    - `correct_energy` is passed to [`monodromy_method_constant_energy`](@ref). Default is `true`.
    - `maxiters` is the maximum number of iterations. Default is `1000`.
    - `energywiggletolerance` is the tolerance for oscillations. Oscillating behaviour indicates
      instability, but sometimes it happens and there must be some tolerance. The default is `1e-2`.
    """
    function follow_PO_family_from_energy(po::PO;
        step::Real=0.1,
        tol::Real=1e-5,
        energytol::Real=1e-3,
        initalperturbation::AbstractVector{<:Real}=Float64[0,1,0,0],
        verbose::Bool =true,
        correct_energy = true,
        dontallowlessenergy::Bool=false)
        system=po.system
        u₀=po.u 
        period = po.T
        H = ClassicalDicke.hamiltonian(system)
        us=[u₀]


        E₀=H(u₀)

        energies=Float64[H(u₀)]
        periods=Float64[period]

        cacheJ=ClassicalSystems.eye(4)
        Λ=[0 0 -1 0
           0 0 0 -1
           1 0 0 0
           0 1 0 0]





        function get(finalEnergy::Real;tol::Real=tol,
            energytol::Real=energytol,
            maxstep=step::Real,step::Real=step, minstep::Real=min(step,energytol/10),
            maxiters::Integer=1000,energywiggletolerance::Real=1e-2)
            index =NaN
            _it=0
            if verbose
                prog=ProgressUnknown("Current Energy")
            end
            if dontallowlessenergy
                finalEnergy=max(finalEnergy,E₀)
            end
            while true #abs(energies[index]-finalEnergy)>energytol
                _it+=1
                if _it>maxiters
                    error("Maximum iterations")
                end
                #index=findmin(abs.(energies.-finalEnergy))[2]
                index=max(searchsortedlast(energies,finalEnergy),1)
                last_e=energies[index]
                falta=abs(last_e-finalEnergy)
                if falta<energytol
                    break
                end
                step=min(step,falta)
                direction=sign(finalEnergy-last_e)

                newu=us[index]
                if verbose
                    prog.desc="Current energy = $(round(last_e,digits=3)), step = 10^$(round(log10(step),digits=1)), iterations ="
                    next!(prog)
                end
                converror=false
               try
                    ClassicalSystems.step(system).jac(cacheJ,newu,[ClassicalSystems.parameters(system);1.0],1.0)
                    grad=hamiltonian_gradient(system,newu)
                    hessian=Λ*cacheJ
                    if index==1
                        Δu=initalperturbation
                    else
                
                       Δu=hamiltonian_gradient(system,newu)
                    end
        
                    Δu=Δu/norm(Δu)
                    ΔE=direction*step


                    a=0.5*transpose(Δu)*hessian*Δu
                    b=dot(grad,Δu)
                    c=- ΔE
                    #Taylor series: -c=Δϵ=b*s + a s^2
                    s=(-b+sqrt(b^2-4*a*c))/(2*a)
                    if abs(H(newu + s*Δu) - (last_e + ΔE)) >abs(H(newu - s*Δu) - (last_e + ΔE))
                        s=-s
                    end
                    lastu=newu
                    E=H(newu+s*Δu)
                    periodGuess=periods[index]
                    if index>1
                        periodGuess+= (E-energies[index-1])*(periods[index]-periods[index-1])/(energies[index]-energies[index-1])
                    end

                    o= monodromy_method_constant_energy(system,newu+s*Δu,periodGuess,inttol=min(1e-12,tol/100),tol=tol,correct_energy=correct_energy)
                    newu = o.u
                    period = o.T
                    if(E-H(newu))>energywiggletolerance
                        if index>1
                            deleteat!(energies,index)
                            deleteat!(us,index)
                            deleteat!(periods,index)
                        end
                        converror=true
                        index-=1
                    end

                catch e
                    converror=true
                end

                notbetter=abs(finalEnergy-H(newu))> falta 
                if  converror || H(newu) in energies || notbetter
        
                    if step/1.2<minstep# (!converror && H(newu) in energies) ||

                    
                        error("The algorithm is not converging. Got to ",PO(system,us[index],periods[index]), " with energy ",H(us[index]) )
                        
                    end
                    step/=1.2
                    continue


                end
                step=min(maxstep,step*1.1)

                nindex=searchsortedlast(energies,finalEnergy) 
                insert!(us,nindex+1,newu)
                insert!(energies,nindex+1,H(newu))
                insert!(periods,nindex+1,period)
                index=nindex+1

            end
            if verbose
                finish!(prog)
            end
            return PO(system,us[index],periods[index])
        end
        return get
    end
    """
    ```julia
    function family_A(system::DickeSystem;kargs...)
    ```
    Returns a function `ϵ -> po` that returns a PO from family ``\\mathcal{A}`` of Ref. [Pilatowsky2021](@cite). This is done
    by passing the ground state fixed point with the positive normal frequency as period to [`follow_PO_family_from_energy`](@ref).
    
    # Arguments
    - `system` should be an instance of [`ClassicalDicke.DickeSystem`](@ref ClassicalDicke.DickeSystem). 
      The system must be in the supperradiant regime. 
    # Keyword arguments
    - `kargs...` are redirected to [`follow_PO_family_from_energy`](@ref)
    """
    family_A(system::DickeSystem;kargs...)=follow_PO_family_from_energy(PO(system,
        ClassicalDicke.minimum_energy_point(system,+),2*pi/ClassicalDicke.normal_frequency(system,+));
        tol=1e-8,
        correct_energy =false, #for some reason, the algorithm with correct_energy=true crashes at energy -1.4 for family_A...
        energytol=1e-3,
        kargs...)
    """
    ```julia
    function family_B(system::DickeSystem;kargs...)
    ```
    Returns a function `ϵ -> po` that returns a PO from family ``\\mathcal{B}`` of Ref. [Pilatowsky2021](@cite). This is done
    by passing the ground state fixed point with the negative normal frequency as period to [`follow_PO_family_from_energy`](@ref).
    
    # Arguments
    - `system` should be an instance of [`ClassicalDicke.DickeSystem`](@ref ClassicalDicke.DickeSystem). 
      The system must be in the supperradiant regime. 
    # Keyword arguments
    - `kargs...` are redirected to [`follow_PO_family_from_energy`](@ref)
    """
    family_B(system::DickeSystem;kargs...)=follow_PO_family_from_energy(PO(system,ClassicalDicke.minimum_energy_point(system,+),2*pi/ClassicalDicke.normal_frequency(system,-));
        tol=1e-8,energytol=1e-3,kargs...)

    """
    ```julia
    function po_coordinates(coords::AbstractArray{<:Integer}, 
        po::PO;tol::Real=1e-12)
    ```
    Returns an array `[Tuple(x1[coords]), Tuple(x2[coords]), ...]`, where is an array containing up
    to 4 integers from 1 to 4, and `xi` are the points in `po` obtained with numerical 
    integration with tolerance `tol`.
    """
    function po_coordinates(coords::AbstractArray{<:Integer}, 
        po::PO;tol::Real=1e-12)
        u=integrate(po)
        us=u.u
        return [Tuple(u[i] for i in coords) for u in us];
    end
    """
    ```julia
    QP(po::PO; tol::Real=1e-12) = po_coordinates([1,3], po, tol=tol)
    ```
    Returns an array `[(Q1, P1), (Q2, P2), ...]` with the coordinates of `po`, integrated with tolerance `tol`.
    This can be directly passed to [`Plots.plot`](https://docs.juliaplots.org/latest/tutorial/)
    """
    QP(po::PO;tol::Real=1e-12)=po_coordinates([1,3],po,tol=tol)
    """
    ```julia
    qp(po::PO; tol::Real=1e-12) = po_coordinates([2,4], po, tol=tol)
    ```
    Returns an array `[(q1, p1), (q2, p2), ...]` with the coordinates of `po`, integrated with tolerance `tol`.
    This can be directly passed to [`Plots.plot`](https://docs.juliaplots.org/latest/tutorial/)
    """
    qp(po::PO;tol::Real=1e-12) = po_coordinates([2,4],po,tol=tol)
    """
    ```julia
    QPp(po::PO; tol::Real=1e-12) = po_coordinates([1,3,4], po, tol=tol)
    ```
    Returns an array `[(Q1, P1, p1), (Q2, P2, p2), ...]` with the coordinates of `po`, integrated with tolerance `tol`.
    This can be directly passed to [`Plots.plot`](https://docs.juliaplots.org/latest/tutorial/) to generate a 3D plot.
    """
    QPp(po::PO; tol::Real=1e-12)=po_coordinates([1,3,4],po,tol=tol)
    """
    ```julia
    QPq(po::PO; tol::Real=1e-12) = po_coordinates([1,3,2], po, tol=tol)
    ```
    Returns an array `[(Q1, P1, q1), (Q2, P2, q2), ...]` with the coordinates of `po`, integrated with tolerance `tol`.
    This can be directly passed to [`Plots.plot`](https://docs.juliaplots.org/latest/tutorial/) to generate a 3D plot.
    """
    QPq(po::PO; tol::Real=1e-12)=po_coordinates([1,3,2], po, tol=tol)

    """
    ```julia
    function overlap_of_tube_with_homogenous_state(po::PO{DickeBCE.QuantumDickeSystem};
                                            time_integral_tolerance::Real=1e-7,
                                            phase_space_integral_resolution::Real=0.1)
    ```
    Returns the overlap ``\\text{tr}(\\hat{\\rho}_\\epsilon \\hat{\\rho}_{\\mathcal{O}} )``
    of a tubular state ``\\hat{\\rho}_{\\mathcal{O}}`` around the periodic orbit ``\\mathcal{O}=`` `po`, (Eq. (15) of Ref. [Pilatowsky2021](@cite)) 
    with a totally delocalized state ``\\hat{\\rho}_\\epsilon`` (Eq. (16) of Ref. [Pilatowsky2021](@cite)).
    # Arguments
    - `po` should be an instance of [`PO`](@ref). The system passed to create `po` 
      should have been a [`QuantumDickeSystem`](@ref DickeBCE.QuantumDickeSystem).
    # Keyword arguments
    - `time_integral_tolerance` is the numerical tolerance for the integral in Eq. (15) of Ref. [Pilatowsky2021](@cite). Default is `1e-7`.
    - `phase_space_integral_resolution` is the phase space resolution for the integral in Eq. (16) of Ref. [Pilatowsky2021](@cite), that is, `res` in
      [`EnergyShellProjections.energy_shell_average`](@ref). Default is `0.1`.
    """
    function overlap_of_tube_with_homogenous_state(po::PO{DickeBCE.QuantumDickeSystem};
        time_integral_tolerance::Real=1e-7,
        phase_space_integral_resolution::Real=0.1)
        
        system=po.system
        orbit=integrate(po,tol=time_integral_tolerance)
        ∫dtHtx(x)=average_over_PO(orbit,u-> DickeBCE.husimi_of_coherent(system,u,x))
        res=phase_space_integral_resolution
        ϵ=energy(po)
        symmP = (po == mirror_Pp(po))
        symmQP = (po == mirror_Qq(po)) && symmP
        return EnergyShellProjections.energy_shell_average(po.system,
                ϵ=ϵ,
                f=∫dtHtx,
                res=res,
                symmetricP =symmP, 
                symmetricQP=symmQP)
    end
    """
    ```julia
    function scarring_measure(
        po::PO{DickeBCE.QuantumDickeSystem},
        quantum_state::AbstractVector{<:Number};
        chop::Real=1e-3,
        kargs...)
    ```
    Returns the scarring measure ``\\mathcal{P}(\\mathcal{O},\\hat{\\rho})`` as defined in Eq. (17) of Ref. [Pilatowsky2021](@cite).
    # Arguments
    - `po` should be an instance of [`PO`](@ref) representing ``\\mathcal{O}`` above.
      The system passed to create `po` should have been a [`QuantumDickeSystem`](@ref DickeBCE.QuantumDickeSystem).
    - `quantum_state` should be a vector representing the quantum state  ``\\hat{\\rho}`` in the coherent efficient basis. 
    # Keyword arguments
    - `chop` is the tolerance to be passed to [`DickeBCE.husimi`](@ref) (see the documentation of [`DickeBCE.coherent_overlap`](@ref))
    - `kargs` are redirected to [`overlap_of_tube_with_homogenous_state`](@ref), which gives the denominator in Eq. (17) of Ref. [Pilatowsky2021](@cite).
    """
    function scarring_measure(
        po::PO{DickeBCE.QuantumDickeSystem},
        quantum_state::AbstractVector{<:Number};
        chop::Real=1e-3,
        kargs...)
        system=po.system
        po_and_state = average_over_PO(po,x -> DickeBCE.husimi(system,x,quantum_state;chop=chop))
        po_and_hom_state = overlap_of_tube_with_homogenous_state(po;kargs...)
        return po_and_state/po_and_hom_state
    end
    
end
