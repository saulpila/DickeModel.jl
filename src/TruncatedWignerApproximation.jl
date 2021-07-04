module TruncatedWignerApproximation
    export calculate_distribution,average,variance,
    mcs_for_survival_probability,survival_probability,Weyl,mcs_chain,
    mcs_for_distributions,monte_carlo_integrate,MonteCarloSystem,coherent_Wigner_HWxSU2

    using Distributed
    using Distributions
    using ProgressMeter
    using DiffEqCallbacks
    using SymEngine
    using ..ClassicalSystems
    using DataStructures
    using ..PhaseSpaces
    """
    ```julia
    struct PhaseSpaceDistribution
    ```
    This object represents a probability distribution in the phase space.
    Currently, the only  implementation is through [`coherent_Wigner_HWxSU2`](@ref TruncatedWignerApproximation.coherent_Wigner_HWxSU2).
    """
    struct PhaseSpaceDistribution
        probability_density
        sample
        ħ::Real
    end
    """
    ```julia
    struct MonteCarloSystem
    ```
    This object may be passed to [`monte_carlo_integrate`](@ref TruncatedWignerApproximation.monte_carlo_integrate).
    Use [`mcs_for_averaging`](@ref TruncatedWignerApproximation.mcs_for_averaging), 
    [`mcs_for_variance`](@ref TruncatedWignerApproximation.mcs_for_variance), and
    [`mcs_for_survival_probability`](@ref TruncatedWignerApproximation.mcs_for_survival_probability) to generate them, 
    and [`mcs_chain`](@ref TruncatedWignerApproximation.mcs_chain) to join them together.
    """
    struct MonteCarloSystem
        f_initial
        f_loop
        reducer
        f_final
        N::Integer
        ts::AbstractVector{<:Real}
        distribution::PhaseSpaceDistribution
        function MonteCarloSystem(f_initial,f_loop,reducer,f_final,N::Integer,ts::Union{Nothing,AbstractVector{<:Real}},distribution::PhaseSpaceDistribution)
            if ts===nothing
                ts=[0.0]
            end
            return new(f_initial,f_loop,reducer,f_final,N,ts,distribution)
        end
    end
    
    function scale_to_index(v,r::LinRange)
            if v<r.start || r.stop <v
                return NaN
            end
            if r.len==1
                return 1
            end
            return Int(round(
                    ((v-r.start)/(r.stop-r.start))*(r.len-1)
                    +1))

        end
    """
    ```julia
    function monte_carlo_integrate(system::ClassicalSystems.ClassicalSystem,
        mc_system::MonteCarloSystem;
        tolerate_errors=true,
        maxNBatch=Inf,kargs...)
    ```
    This function is the backbone of this module. It performs a Monte Carlo integration-type procedure.
    The argument `mc_system` determines a `distribution::`[`PhaseSpaceDistribution`](@ref), an integer `N`, and a list of times `ts`.
    This function calls an initializing function also determined by `mc_system`. Then, it samples
    `N` random points `x` from `distribution`, which are integrated using the Hamiltonian given by `system`.
    For each trajectory and time in `ts`, an operation is performed, which is determined again by `mc_system`.
    A final operation is then performed and the result is returned. This may seem a bit abstract,
    but it is a very flexible system. The applications are more concrete; for example, if you
    generate `mc_system` using [`mcs_for_averaging(... observable=f ...)`](@ref TruncatedWignerApproximation.mcs_for_averaging), 
    then the initial operation is to generate an array of zeroes the same length as `ts`. Then, for each initial condition
    `x`, the result of `f(x(t))` is added to each element of the array. Finally, the array is overall divided by `N` and then
    returned. This is exactly a Monte Carlo integration of the function `f(x)` over the classical evolution 
    of the distribution.
    
    This function uses all the available workers, but make sure to import this module in all of them.
    
    # Arguments:
    - `system` should be an instance of [`ClassicalSystems.ClassicalSystem`](@ref).
    - `mc_system` should be an instance of  [`MonteCarloSystem`](@ref TruncatedWignerApproximation.MonteCarloSystem).
    - `tolerate_errors` indicates that some errors in the integration may be ignored.
      This is useful because sometimes a one-in-a-million numerical instability may arise, and
      you may want to ignore it. If more than `100` errors occur consecutively, then then the
      procedure is stopped. Defaults to `true`.
    - `maxNBatch` is the maximum number of batch-sizes sent to each worker. Defaults to `inf`.
    - `kargs` are redirected to [`ClassicalSystems.integrate`].
    """
    function monte_carlo_integrate(system::ClassicalSystems.ClassicalSystem,
        mc_system::MonteCarloSystem;
        tolerate_errors=true,
        maxNBatch=Inf,kargs...)
        
        f_inicial=mc_system.f_initial
        f_loop=mc_system.f_loop
        reductora= mc_system.reducer
        finalizadora= mc_system.f_final

        distribution_sampler = mc_system.distribution.sample
        ts=mc_system.ts
        N=mc_system.N
        show_progress= nothing
        trabajadores=length(workers())
        if trabajadores==0
            trabajadores=1
        end
        Neff=min(Int(ceil(N/trabajadores)),maxNBatch)
        workerload=[]
        rest=N
        while rest>0
            load=Neff
            rest=rest-load
            extra=-min(0,rest)
            rest=max(rest,0)
            push!(workerload,load-extra)
        end
        channel = RemoteChannel(()->Channel{Int8}(N), 1)
        res=0

        workersready=0
        missedprogressloop=trabajadores
        @sync begin
            # this task prints the progress bar
            @async while true
                v= take!(channel)
                if v==0
                    if show_progress!=nothing
                        finish!(show_progress)
                    end
                    break
                elseif v==1
                    if workersready>=trabajadores
                        if show_progress==nothing
                            show_progress=Progress(Int(N-missedprogressloop), 1)
                        end
                        next!(show_progress)
                    else
                        missedprogressloop+=1
                    end

                else
                    workersready+=1
                end
            end

            # this task does the computation
            @async begin
               try
                   res=@distributed reductora for load in workerload
                        val=f_inicial()
                        ti=1
                        ti_error=1
                        nerrores=0
                        call_loop= (u,t,_)-> begin
                            if ti>=ti_error
                                f_loop(val,u,ti)
                            end
                            ti+=1
                        end


                        j=1
                        notimeevolution= (ts==[0.0])
                        cb=nothing
                        if !notimeevolution
                             cb=FunctionCallingCallback(call_loop;funcat=ts,func_start = false)
                        end
                        while j <=load
                            u₀=distribution_sampler()

                            ti=1
                            try
                                if notimeevolution
                                    call_loop(u₀,0.0,0)
                                else
                                    ClassicalSystems.integrate(system,t=last(ts),save_everystep=false,u₀=u₀,callback=cb;kargs...)
                                end
                                ti_error=1

                                put!(channel,  j==1 ? 2 : 1 )
                                j+=1
                                nerrores=0

                            catch e
                                if tolerate_errors==false
                                    error(e)
                                end
                                if nerrores>100
                                    println("Too many errors")
                                    error(e)
                                    
                                end
                                nerrores+=1
                                ti_error=ti
                            end

                        end
                        val
                    end


                finally
                    put!(channel, 0) # this tells the printing task to finish

                end
            end
        end

        return finalizadora(res)
    end
    """
    ```julia
    function mcs_chain(mcs1::MonteCarloSystem, mcs2::MonteCarloSystem, ...)
    ```
    Functionally equivalent to [`mcs_chain`](@ref mcs_chain(::AbstractArray{MonteCarloSystem,1}))`([mcs1, mcs2, ...])`.
    """
    function mcs_chain(mcss::Vararg{MonteCarloSystem,N}) where {N}
        mcs_chain(mcss)
    end
    function mcs_chain(iterable)
        mcs_chain(collect(iterable))
    end
    """
    ```julia
    function mcs_chain(Fs::AbstractArray{MonteCarloSystem,1})
    ```
    Generates a [`MonteCarloSystem`](@ref) by chaining together those in the array `Fs = [monteCarloSystem1, monteCarloSystem2, ...]`.
    The output of the system generated is an array that has the outputs of each system, in the same order as `Fs`.
    """
    function mcs_chain(Fs::AbstractArray{MonteCarloSystem,1}) 
        ts= Fs[1].ts
        N= Fs[1].N
        distribution=Fs[1].distribution
        if any(N!=F.N for F in Fs)
            error("All systems must have the same N")
        end
        if any(ts!=F.ts for F in Fs)
            error("All systems must have the same ts")
        end
        if any(distribution!=F.distribution for F in Fs)
            error("All systems must have the same phase space distribution")
        end
        f_inicial=function()
            [F.f_initial() for F in Fs]
        end
        f_loop=function(valores,u,t)
            for i in 1:length(Fs)
                Fs[i].f_loop(valores[i],u,t)
            end
        end
        reductora=function(last,next)
             [Fs[i].reducer(last[i],next[i]) for i in 1:length(Fs)]

        end
        f_final=function(res)
            [Fs[i].f_final(res[i]) for i in 1:length(Fs)]
        end
        return MonteCarloSystem(f_inicial,f_loop,reductora,f_final,N,ts,distribution)
    end
    function function_from_expression(system::ClassicalSystems.ClassicalSystem,expression)
        varnames=ClassicalSystems.varnames(system)
        if typeof(expression)!=SymEngine.Basic && !isempty(methods(expression))#es una funcion
            if any(length(m.sig.parameters)-1==length(varnames) for m in methods(expression))
                return expression
            elseif any(length(m.sig.parameters)-1==1 for m in methods(expression))
                f(x...)=expression(x)
                return f
            else
                error("You should pass a function that takes a n-vector or n arguments, where n is the dimension of the phase space.")
            end

        end

        if !(Symbol.(SymEngine.free_symbols(SymEngine.Basic(expression))) ⊆ varnames)
            error("The expression $(expression) has unknown variables. Write it in terms only of $(ClassicalSystems.varnames(system))")
        end
        if  typeof(expression)==SymEngine.Basic
            expression=SymEngine.toString(expression) #esto no es lo mejor, seria mejor usar expresiones y no strings pero algo pasa con las raices cuadradas...

        end
        return Base.eval(Main,Meta.parse("($(join(varnames,',')))->($expression)"))
    end
    function functions_from_expressions(system::ClassicalSystems.ClassicalSystem,expressions)
        funciones_varianza=[function_from_expression(system,s) for s in expressions]
        return funciones_varianza
    end
    """
    ```julia
    function mcs_for_averaging(
        system::ClassicalSystems.ClassicalSystem;
        observable,
        ts::AbstractVector{<:Real}=[0.0],
        N::Integer,
        distribution::PhaseSpaceDistribution)
    ```
    Generates a [`MonteCarloSystem`](@ref) that computes
    ```math
        \\frac{1}{N} \\sum_{i=1}^N \\text{observable}(\\mathbf{x}_i(t))
    ```
    for each ``t`` in `ts`, where ``\\mathbf{x}_i`` will be sampled from `distribution`. The system produces an array
    the same size as `ts`, containing the result for each time.
    
    # Arguments
    - `system` should be an instance of [`ClassicalSystems.ClassicalSystem`](@ref),
    - `observable` can a function in the form `f(x::Vector)` or `f(x1,x2 ,..., x_n)` with 
      n the dimension of the phase space determined by `system`. `observable` can also be
      an expression determining the operations between the `varnames` determined by system. For example,
      if `system` were a instance of [`ClassicalDicke.ClassicalDickeSystem`](@ref Main.ClassicalDicke.ClassicalDickeSystem), then `observable`
      could be `:(q+p^2 +Q)`, `f(Q, q, P, p) = q + p^2 - Q` or `f(x) = x[2] + x[4]^2 - x[1]`, which
      are all equivalent.
      See also the submodule [`Weyl`](@ref Dicke.TruncatedWignerApproximation.Weyl), which produces expressions corresponding to 
      the [Weyl symbols](https://en.wikipedia.org/wiki/Wigner%E2%80%93Weyl_transform#The_inverse_map) of quantum observables.  
      Note: `observable` can also be an array of observables, in which case an array is returned for each time.
    - `ts` should be a sorted array of times.
    - `N` is the number of points to sample. The bigger the more accurate.
    - `distribution` should be an instance of [`PhaseSpaceDistribution`](@ref)
    """
    function mcs_for_averaging(
        system::ClassicalSystems.ClassicalSystem;
        observable,
        ts::AbstractVector{<:Real}=[0.0],
        N::Integer,
        distribution::PhaseSpaceDistribution)
        onevar=false

        if !isa(observable,Array)
            observable=[observable]
            onevar=true
        end

        lengthexpr=length(observable)
        f_inicial= function()
            [zeros(lengthexpr) for i in 1:length(ts)]
        end
        varnames=ClassicalSystems.varnames(system)
        funciones_a_promediar=functions_from_expressions(system,observable)
        f_loop=function(valores,u,i)

            for fi in 1:lengthexpr
                valores[i][fi]+=funciones_a_promediar[fi](u...)
            end
        end
        f_final=function(res)
            res=res/N
            if onevar
                res=[re[1] for re in res]
            end
            return res
        end

        return MonteCarloSystem(f_inicial,f_loop,+,f_final,N,ts,distribution)
    end
    """
    ```julia
    function mcs_for_distributions(system::ClassicalSystems.ClassicalSystem;
        x,
        xs::Union{AbstractRange{<:Real},Nothing}=nothing,
        y=nothing,
        ts::Union{AbstractRange{<:Real},Nothing}=nothing,
        animate::Bool=false,
        ys::Union{AbstractRange{<:Real},Nothing}=nothing,
        N::Integer,
        distribution::PhaseSpaceDistribution)
    ```
    Generates a [`MonteCarloSystem`](@ref) that produces a multidimensional histogram
    allowing to visualize the expected value of observables under [`PhaseSpaceDistribution`](@ref) `distribution`.
    
    Several behaviors can be produced.
    - If `x` is `:t`, a matrix of dimensions `length(ts)` ``\\times`` `length(ys)` 
      is produced. The coordinate `(tᵢ,yⱼ)` gives the probability of finding the observable
      `y` between `ys[j]` and `ys[j+1]` at the time `ts[i]`.
    - If `y` is `:t`, the same as above, changing `x` by `y`.
    - If neither `x` nor `y` are `:t`, and `animate = false`, then a matrix of dimensions `length(xs)` ``\\times`` `length(ys)` 
      is produced.  The coordinate `(xᵢ,yⱼ)` gives the probability of finding the observable
      `x` between `xs[i]` and `xs[i+1]` and the observable
      `y` between `ys[j]` and `ys[j+1]`, averaging over all the times in `ts`
    - If neither `x` nor `y` are `t`, and `animate = true` then an array of matrices of dimensions `length(xs)` ``\\times`` `length(ys)` 
      is produced, the same size as `ts`. For the `n`th matrix in this array, the coordinate `(xᵢ,yⱼ)` gives the probability of finding the observable
      `x` between `xs[i]` and `xs[i+1]` and the observable
      `y` between `ys[j]` and `ys[j+1]` at time `ts[n]`.
    
    # Arguments
    - `system` should be an instance of [`ClassicalSystems.ClassicalSystem`](@ref),
    - `x` and `y` can be `:t`, or an observable as described in the arguments of [`mcs_for_averaging`](@ref).
    - `xs` and `ys` should be given if  `x` and  `y` are not  `:t`, in wich case they should be
      range objects (e.g. `0:0.1:1`) containing the bins for the histogram.
    - `ts` should be a range object (e.g. `0:0.1:1`) of times.
    - `animate` should be a `Bool`, which determines the behavior if neither `x` nor `y` are `t`. Defaults to `true`
    - `N` is the number of points to sample. The bigger the more accurate.
    - `distribution` should be an instance of [`PhaseSpaceDistribution`](@ref)
    
    Note: If `ts` and `y` are both not passed, then `y` is set to `:t` and `ts` to `0:0`.
    """
    function mcs_for_distributions(system::ClassicalSystems.ClassicalSystem;x,xs::Union{AbstractRange{<:Real},Nothing}=nothing,
        y=nothing,ts::Union{AbstractRange{<:Real},Nothing}=nothing,animate::Bool=false,ys::Union{AbstractRange{<:Real},Nothing}=nothing,N::Integer,distribution::PhaseSpaceDistribution)
        onematrix=false
        #vars=[]
        if ts==nothing && y==nothing
            y=:t
            ts=0:0
        end
        if  ts==nothing
            error("You have to pass :ts")
        end

        if x==:t || y==:t || animate == false
            onematrix=true
        
            if x==:t
                xs=ts
            end
            if y==:t
                ys=ts
            end
        end
        xs=LinRange(xs)
        ys=LinRange(ys)
        f_inicial= function()
            if onematrix
                matrices=zeros(length(ys),length(xs))
            else
                matrices= [DefaultDict{Tuple{UInt16,UInt16},UInt32}(0) for i in ts]
            end
            return matrices
        end
        fx=nothing
        fy=nothing
        if x!=:t
            fx=functions_from_expressions(system,[x])[1]
        end
        if y!=:t
            fy=functions_from_expressions(system,[y])[1]
        end
        f_loop=nothing
        if onematrix==true
            if x==:t
                f_loop=function(matrices,u,i)
                    indu=scale_to_index(fy(u...),ys)
                    if !isnan(indu)
                        matrices[indu,i]+=1
                    end
                end
            elseif y==:t
                f_loop=function(matrices,u,i)
                    indu=scale_to_index(fx(u...),xs)
                    if !isnan(indu)
                        matrices[i,indu]+=1
                    end
                end
            else
                f_loop=function(matrices,u,i)
                    iy=scale_to_index(fy(u...),ys)
                    ix=scale_to_index(fx(u...),xs)
                    if !isnan(iy) && !isnan(ix)
                        matrices[iy,ix]+=1
                    end
                end
            end

        else
            f_loop=function(matrices,u,i)
                iy=scale_to_index(fy(u...),ys)
                ix=scale_to_index(fx(u...),xs)

            
                if !isnan(iy) && !isnan(ix)

                    matrices[i][iy,ix]+=1
                end
            end
            
        end
        if onematrix!=true

            datachannel = RemoteChannel(()->Channel{Any}(Inf), 1)

            resdat=[zeros(length(ys),length(xs)) for i in 1:length(ts)]
            flush_task=@async while true
                x= take!(datachannel)
                if x==nothing
                    close(datachannel)
                    break
                else
                    ra=rand(1:100)
                    for ri in length(resdat):-1:1
                        r=resdat[ri]
                        m=pop!(x)
                        for (k,v) in m
                            r[k...]+=v/N
                        end
                        yield()
                    end
                end
            end
       end

        function finalizacion(x)
            if onematrix==true
                return x/N
            end
            if length(workers())<=1 #if there is 1 worker, add does not get called
                add(x,nothing)
            end 
            put!(datachannel,nothing)
            wait(flush_task)
            return resdat
        end
        function add(a,b)
            if onematrix==true
                a.+=b
                return a
            end

            for x in [a,b]
                if x===nothing
                    continue
                end
                
                put!(datachannel,x)
            end
            return nothing
        end
        return MonteCarloSystem(f_inicial,f_loop,add,finalizacion,N,ts,distribution)
    end
    """
    ```julia
    function calculate_distribution(system::ClassicalSystems.ClassicalSystem, 
        [[[same kargs as mcs_for_distributions]]],
        kargs...)
    ```
    Calls [`mcs_for_distributions`](@ref) and then [`monte_carlo_integrate`](@ref) on the resulting `MonteCarloSystem`. Extra `kargs` are sent to the latter.
    """
    function calculate_distribution(system::ClassicalSystems.ClassicalSystem;x,N::Integer,
        animate::Bool=false,xs::Union{AbstractRange{<:Real},Nothing}=nothing,y=nothing,
        ts::Union{AbstractRange{<:Real},Nothing}=nothing,ys::Union{AbstractRange{<:Real},Nothing}=nothing,distribution::PhaseSpaceDistribution,kargs...)
        return monte_carlo_integrate(system,mcs_for_distributions(system;x=x,y=y,ts=ts,xs=xs,ys=ys,N=N,distribution=distribution,animate=animate);kargs...)
    end
    """
    ```julia
    function average(system::ClassicalSystems.ClassicalSystem, 
        [[[same kargs as mcs_for_averaging]]],
        kargs...)
    ```
    Calls [`mcs_for_averaging`](@ref) and then [`monte_carlo_integrate`](@ref) on the resulting `MonteCarloSystem`. Extra `kargs` are sent to the latter.
    """
    function average(system::ClassicalSystems.ClassicalSystem;observable,N::Integer,ts::AbstractArray{<:Real}=[0.0],distribution::PhaseSpaceDistribution,kargs...)
        return monte_carlo_integrate(system,
            mcs_for_averaging(system;observable=observable,ts=ts,N=N,distribution=distribution);kargs...)
    end
    """
    ```julia
    function mcs_for_variance(
        system::ClassicalSystems.ClassicalSystem;
        observable,
        ts::AbstractVector{<:Real}=[0.0],
        N::Integer,
        distribution::PhaseSpaceDistribution,
        return_average::Bool = false)
    ```
    Generates a [`MonteCarloSystem`](@ref) that computes the variance
    ```math
        \\frac{1}{N} \\sum_{i=1}^N \\text{observable}^2(\\mathbf{x}_i(t)) - \\left (\\frac{1}{N}\\sum_{i=1}^N \\text{observable}(\\mathbf{x}_i(t))\\right )^2
    ```
    for each ``t`` in `ts`, where ``\\mathbf{x}_i`` will be sampled from `distribution`. The system produces an array
    the same size as `ts`, containing the result for each time.
    
    # Arguments
      - See arguments for [`mcs_for_averaging`](@ref).
      - If `return_average=false` (default), only produces the variance, else it 
        will produce a tuple, where the first element is the variance and the second
        the average.
    """
    function mcs_for_variance(sistema::ClassicalSystems.ClassicalSystem;observable,N::Integer,ts::AbstractArray{<:Real}=[0.0],
        return_average::Bool=false,distribution::PhaseSpaceDistribution)
        onevar=false

        if !isa(observable,Array)
            observable=[observable]
            onevar=true
        end
        tam=length(observable)

        apromediar=[observable;[(x)->f(x...)^2 for f in TruncatedWignerApproximation.function_from_expression.((sistema,),observable)]]
        mcsa=TruncatedWignerApproximation.mcs_for_averaging(sistema,observable=apromediar,ts=ts,N=N,distribution=distribution)
        mf_final=function(res)
            r=mcsa.f_final(res)
            vars=[[re[i+tam] - re[i]^2 for i in 1:tam] for re in r]
            if onevar
                vars=[v[1] for v in vars]
            end
            if return_average
                proms=[[re[i] for i in 1:tam] for re in r]
                if onevar
                    proms=[p[1] for p in proms]
                end
                vars=vars,proms
            end

            return vars
        end
        return MonteCarloSystem(mcsa.f_initial,mcsa.f_loop,+,mf_final,N,ts,distribution)
    end
    """
    ```julia
    function variance(system::ClassicalSystems.ClassicalSystem, 
        [[[same kargs as mcs_for_variance]]],
        kargs...)
    ```
    Calls [`mcs_for_variance`](@ref) and then [`monte_carlo_integrate`](@ref) on the resulting `MonteCarloSystem`. Extra `kargs` are sent to the latter.
    """
    function variance(sistema::ClassicalSystems.ClassicalSystem;observable,N::Integer,ts::AbstractArray{<:Real}=[0.0],distribution::PhaseSpaceDistribution,return_average=false,kargs...)
        return TruncatedWignerApproximation.monte_carlo_integrate(sistema,
            mcs_for_variance(sistema;observable=observable,ts=ts,N=N,distribution=distribution,return_average = return_average);kargs...)
    end
    """
    ```julia
    function mcs_for_survival_probability(
        system::ClassicalSystems.ClassicalSystem;
        N::Integer,
        ts::AbstractArray{<:Real},
        distribution::PhaseSpaceDistribution)
    ```
    Generates a [`MonteCarloSystem`](@ref) that computes the survival probability
    through Eq. (C.7) of Ref. [Villasenor2020](@cite) (with ``M=`` `N`, ``w=`` `distribution`).
    
    # Arguments
    - `system` should be an instance of [`ClassicalSystems.ClassicalSystem`](@ref),
    - `ts` should be a sorted array of times.
    - `N` is the number of points to sample. The bigger the more accurate.
    - `distribution` should be an instance of [`PhaseSpaceDistribution`](@ref)
    """
    function mcs_for_survival_probability(system::ClassicalSystems.ClassicalSystem;
        N::Integer,
        ts::AbstractArray{<:Real},
        distribution::PhaseSpaceDistribution)
        f_inicial= function()
            zeros(length(ts))
        end
        f_loop=function(valores,u,i)
                valores[i]+=distribution.probability_density(u)
        end
        f_final=function(res)
            L=length(ClassicalSystems.varnames(system))/2
            return res*(2*π*distribution.ħ)^L/N
        end

        return MonteCarloSystem(f_inicial,f_loop,+,f_final,N,ts,distribution)
    end
    """
    ```julia
    function survival_probability(system::ClassicalSystems.ClassicalSystem, 
        [[[same kargs as mcs_for_survival_probability]]],
        kargs...)
    ```
    Calls [`mcs_for_survival_probability`](@ref) and then [`monte_carlo_integrate`](@ref) on the resulting `MonteCarloSystem`. Extra `kargs` are sent to the latter.
    """
    function survival_probability(system::ClassicalSystems.ClassicalSystem;distribution::PhaseSpaceDistribution,N::Integer,ts::AbstractArray{<:Real},kargs...)
        return TruncatedWignerApproximation.monte_carlo_integrate(system,
            mcs_for_survival_probability(system;ts=ts,N=N,distribution=distribution);integate_backwards=true,kargs...)
    end
    """
    ```julia
    function coherent_Wigner_HW(;q₀::Real,p₀::Real,j::Real=1,ħ::Real=1/j)
    ```
    Returns a [`PhaseSpaceDistribution`](@ref) corresponding to the two-dimensional Wigner
    function of a coherent state of the Heisenberg-Weyl algebra (i.e. a standard coherent state) 
    centered at `q₀,p₀`. A value of `ħ` may be passed, or pass `j` to set `ħ = 1/j`. (See Eq. (B.1) of Ref. [Villasenor2020](@cite))
    """
    function coherent_Wigner_HW(;q₀::Real,p₀::Real,j::Real=1,ħ::Real=1/j)
        σ=sqrt(ħ/2)
        gaussiana_qp=Distributions.MvNormal([q₀,p₀], σ)
        function W((q,p))
                return Distributions.pdf(gaussiana_qp,[q,p])
        end
        function sample()
            q,p=Distributions.rand(gaussiana_qp)
            return [q,p]
        end
        return PhaseSpaceDistribution(W,sample,ħ)
    end
    """
    ```julia
    function coherent_Wigner_HW(;Q₀::Real,P₀::Real,j::Real=1,ħ::Real=1/j)
    ```
    Same as [`coherent_Wigner_HW`](@ref), but for a coherent state for the SU(2) 
    algebra (a Bloch coherent state).  The approximation given by Eq. (B.4) of Ref. [Villasenor2020](@cite) is used.
    """
    function coherent_Wigner_SU2(;Q₀::Real,P₀::Real,j::Real=1,ħ::Real=1/j)
        σ=sqrt(ħ/2)
        gaussiana_Θ=Distributions.Normal(0, σ)
        Rayleigh_Θ=Distributions.Rayleigh(σ)
        R(θ,ϕ,θ0,ϕ0)= #rotacion que manda el polo norte a θ0,ϕ0 aplicada a θ,ϕ
            [atan(sqrt((cos(θ0)*cos(ϕ)*sin(θ)-cos(θ)*sin(θ0))^2+sin(θ)^2*sin(ϕ)^2),cos(θ)*cos(θ0)+cos(ϕ)*sin(θ)*sin(θ0)),
            mod(atan(-cos(ϕ0)*sin(θ)*sin(ϕ)-cos(θ0)*cos(ϕ)*sin(θ)*sin(ϕ0)+cos(θ)*sin(θ0)*sin(ϕ0),-cos(θ0)*cos(ϕ)*cos(ϕ0)*sin(θ)+cos(θ)*cos(ϕ0)*sin(θ0)+sin(θ)*sin(ϕ)*sin(ϕ0)),2*pi)]


        centroθ=PhaseSpaces.θ_of_QP(Q₀,P₀)
        centroϕ=PhaseSpaces.ϕ_of_QP(Q₀,P₀)
        function Wθ(θ,ϕ)
            try
                Θ=acos(cos(θ)*cos(centroθ)+ sin(θ)*sin(centroθ)*cos(ϕ - centroϕ))
                return Distributions.pdf(gaussiana_Θ,Θ)*sqrt(2*π)*σ *(1/(2*π*σ^2))
            catch
                return 0
            end
        end
        function W((Q,P))
            try

                 θ=PhaseSpaces.θ_of_QP(Q,P)
                 ϕ=PhaseSpaces.ϕ_of_QP(Q,P)
                 return Wθ(θ,ϕ)
            catch
                return 0
            end
        end


        function sample()
            Θ=Distributions.rand(Rayleigh_Θ)
            θ,ϕ = R(Θ,rand()*2*pi,centroθ,centroϕ)
            return [PhaseSpaces.Q_of_θϕ(θ,ϕ),PhaseSpaces.P_of_θϕ(θ,ϕ)]
        end
        return PhaseSpaceDistribution(W,sample,ħ)
    end
    """
    ```julia
    function coherent_Wigner_HWxSU2(;Q₀::Real,q₀::Real,P₀::Real,p₀::Real,j::Real=1,ħ::Real=1/j)
    ```
    Produces the Wigner function corresponding to the tensor product of [`coherent_Wigner_HW`](@ref) and [`coherent_Wigner_SU2`](@ref).
    """
    function coherent_Wigner_HWxSU2(;Q₀::Real,q₀::Real,P₀::Real,p₀::Real,j::Real=1,ħ::Real=1/j)
        WSU2=coherent_Wigner_SU2(Q₀=Q₀,P₀=P₀,j=j,ħ=ħ)
        WHW=coherent_Wigner_HW(q₀=q₀,p₀=p₀,j=j,ħ=ħ)
        W((Q,q,P,p))=WSU2.probability_density((Q,P))*WHW.probability_density((q,p))
        function sample()
            Q,P=WSU2.sample()
            q,p=WHW.sample()
            return [Q,q,P,p]
        end
        return PhaseSpaceDistribution(W,sample,ħ)
    end
    """
    ```julia
    function  coherent_Wigner_HWxSU2(u₀::AbstractVector{<:Real},j::Real=1,ħ::Real=1/j)
    ```
    Same as [`coherent_Wigner_HWxSU2`](@ref coherent_Wigner_HWxSU2()), taking `u₀ = [Q₀, q₀, P₀, p₀]`.
    """
    coherent_Wigner_HWxSU2(u₀::AbstractVector{<:Real};j::Real=1,ħ::Real=1/j)=coherent_Wigner_HWxSU2(Q₀=u₀[1],q₀=u₀[2],P₀=u₀[3],p₀=u₀[4],j=j,ħ=ħ)
    
    """
    ```julia
    function  coherent_Wigner_SU2(u₀::AbstractVector{<:Real},j::Real=1,ħ::Real=1/j)
    ```
    Same as [`coherent_Wigner_SU2`](@ref coherent_Wigner_SU2()), taking `u₀ = [Q₀, P₀]`.
    """
    coherent_Wigner_SU2(u₀;j::Real=1,ħ::Real=1/j)=coherent_Wigner_SU2(Q₀=u₀[1],P₀=u₀[2],j=j,ħ=ħ)
    
    """
    ```julia
    function  coherent_Wigner_HW(u₀::AbstractVector{<:Real},j::Real=1,ħ::Real=1/j)
    ```
    Same as [`coherent_Wigner_HW`](@ref coherent_Wigner_HW()), taking `u₀ = [q₀, p₀]`.
    """    
    coherent_Wigner_HW(u₀;j::Real=1,ħ::Real=1/j)=coherent_Wigner_HW(q₀=u₀[1],p₀=u₀[2],j=j,ħ=ħ)
    
    module Weyl
        using SymEngine
        SymEngine.@vars q p
        """
        ```julia
        function n(j::Real)
        ```
        Returns the Weyl symbol of the number operator ``W(n̂)``, `(p²+  q² - ħ)/2`, where `ħ = 1/j`.  
        """    
        function n(j::Real)
            SymEngine.@vars q p
            return (j*(p^2 + q^2) - 1)/2 #Weyl(n̂)=(p²+q²- ħ)/2 y nuestras q,p están escaladas por sqrt(j)
        end
        """
        ```julia
        function n²(j::Real)
        ```
        Returns the Weyl symbol of the number operator squared ``W(n̂^2) = W(n̂)^2 - \\hbar^2/4``, where ``\\hbar = 1/j``. (Note that the extra term ``\\hbar^2/4`` 
        appears due to the non-commutativity of ``q̂`` and ``p̂``. See Ref. [Polkovnikov2013](@cite) for details on how to compute these expressions.
        """    
        function n²(j::Real)
            return n(j)^2 - 1/4 
        end


        #"The Moyal Representation for Spin" JOSEPH C. VÁRILLY AND JOSÉ M. GRACIA-BONDÍA Ann. of Phys. 190, 107-148 (1989)
        #https://doi.org/10.1016/0003-4916(89)90262-5
        macro SU2angles(ex)
            q= quote
                SymEngine.@vars Q P
                cosθ=(P^2+Q^2)/2-1
                sinθ=sqrt(1-cosθ^2)
                cosϕ=Q/sqrt(P^2 + Q^2)
                sinϕ=-P/sqrt(P^2 + Q^2)
                $ex
            end
            return esc(q)
        end
        """
        ```julia
        function Jz(j::Real)
        ```
        Returns the Weyl symbol of the operator ``\\hat{J}_z``. (See p. 114 of Ref. [Varilly1989](@cite))
        """    
        function Jz(j::Real)
            return @SU2angles sqrt(j*(j+1))*cosθ
        end
        """
        ```julia
        function Jx(j::Real)
        ```
        Returns the Weyl symbol of the operator ``\\hat{J}_x``. (See p. 114 of Ref. [Varilly1989](@cite))
        """
        function Jx(j::Real)
            return @SU2angles sqrt(j*(j+1))*sinθ*cosϕ
        end
        """
        ```julia
        function Jy(j::Real)
        ```
        Returns the Weyl symbol of the operator ``\\hat{J}_y``. (See p. 114 of Ref. [Varilly1989](@cite))
        """
        function Jy(j::Real)
            return @SU2angles sqrt(j*(j+1))*sinθ*sinϕ
        end
        """
        ```julia
        function Jz²(j::Real)
        ```
        Returns the Weyl symbol of the operator ``\\hat{J}_z^2``. (See bottom of p. 128 of Ref. [Varilly1989](@cite))
        """
        function Jz²(j::Real)
            bj=sqrt(j*(j+1)*(2*j-1)*(2*j+3))
            return @SU2angles (j*(j+1))/3 + (bj/2)*(cosθ^2 - 1/3)
        end
        """
        ```julia
        function Jx²(j::Real)
        ```
        Returns the Weyl symbol of the operator ``\\hat{J}_x^2``. (See bottom of p. 128 of Ref. [Varilly1989](@cite))
        """
        function Jx²(j::Real)
            bj=sqrt(j*(j+1)*(2*j-1)*(2*j+3))
            return @SU2angles (j*(j+1))/3 + (bj/2)*((sinθ*cosϕ)^2 - 1/3)
        end
        """
        ```julia
        function Jy²(j::Real)
        ```
        Returns the Weyl symbol of the operator ``\\hat{J}_y^2``. (See bottom of p. 128 of Ref. [Varilly1989](@cite))
        """
        function Jy²(j::Real)
            bj=sqrt(j*(j+1)*(2*j-1)*(2*j+3))
            return @SU2angles (j*(j+1))/3 + (bj/2)*((sinθ*sinϕ)^2 - 1/3)
        end
    end
end




