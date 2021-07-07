module DickeHusimiProjections

export ∫∫dqdpδϵ,proj_husimi_QP_matrix,loc_measure
import ..DickeBCE
import ..ClassicalSystems
import Distributed
using ..ClassicalDicke
using Statistics
using Random
import ProgressMeter

    """
    ```julia
    function ∫∫dqdpδϵ(system::ClassicalDickeSystem;
        ϵ::Real,
        Q::Real,
        P::Real,
        f::Function,
        p_res::Real=0.01,
        nonvalue=nothing,
        onlyqroot::Union{typeof(-),typeof(+),Nothing}=nothing)
    ```
    Computes the integral  
    ```math 
        \\int \\int dq dp f(\\mathbf{x}) \\delta(h_\\text{cl}(\\mathbf{x}) - \\epsilon)
    ```
    using Eq. (8) of Ref. [Pilatowsky2021NatCommun](@cite) and a 
    [Chebyshev–Gauss quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature)
    method.  
    
    Note: to type `∫` in Julia, type `\\int` + Tab.
    # Arguments
    - `system` should be an instance of [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `ϵ` is a real number corresponding to ``\\epsilon`` above.
    - `Q` and `P` are the ``Q`` and ``P`` coordinates of ``\\mathbf{x}`` above.
    - `f` should be function with signature `f([Q,q,P,p])` that returns values that
      may be added together (Numbers, Arrays etc...).
    - `p_res` is the resolution of the integral. The amount of nodes for the
      quadrature will be equal to `2*floor(p₊/p_res) + 1`, where `p₊` is defined
      below Eq. (9) of Ref. [Pilatowsky2021NatCommun](@cite). Default is `0.01`.
    - `nonvalue` is the value to return if ``Q,P`` are outside of the energy shell
      at ``\\epsilon``. Default is `NaN`.
    - `onlyqroot` may be `+`, `-`, or `nothing` (default). If it is `+` (`-`),
      the integral is only computed over the positive (negative) roots in `q`, that is,
      only the ``q_+`` (``q_-``) term is taken in Eq. (8) of Ref. [Pilatowsky2021NatCommun](@cite).
      If it is `nothing`, both terms are added.
    """
    function ∫∫dqdpδϵ(system::ClassicalDickeSystem,
        ϵ::Real,
        Q::Real,
        P::Real,
        f::Function,
        p_res::Real=0.01,
        nonvalue = NaN,
        onlyqroot::Union{typeof(-),typeof(+),Nothing}=nothing)
        if Q^2+P^2>4
            return nonvalue
        end
        ω₀, ω, γ=ClassicalSystems.parameters(system)
        val=nonvalue #variable donde vamos a guardar la suma

        A= (4*γ^2*Q^2*(1-(Q^2+P^2)/4)  - ω*ω₀*(P^2 + Q^2 - 2) + 2*ω*ϵ)
        if A<0
            return nonvalue #estamos fuera del espacio fase permitido
        end
        p₊=sqrt(A)/ω-1e-10   #p₊ y p₋=-p₊ son los limites de integracion sobre p #el 1e-12 es para evitar errores numericos
        if p₊<0
            return nonvalue
        end
        n=Int(2*floor(p₊/p_res) + 1)

        for signo in [+,-]
            if onlyqroot!==nothing && onlyqroot!==signo
                continue
            end
            mf(p)=f(ClassicalDicke.Point(system,Q=Q,P=P,p=p,ϵ=ϵ,signo=signo))
            for k in 1:n
                mv= (π/(n*ω)) .* mf(p₊*cos(π*(2*k - 1)/(2n)))
                if val===nonvalue
                    val=deepcopy(mv) #copiamos por si mv es un cache (si no es un cache shame on you, usa un cache)
                else
                    if !isa(val, Array)
                        val+=mv
                    elseif isa(val[1], Array)
                        for (a,b) in zip(val,mv)
                            a.+=b
                        end
                    else
                        val.+=mv
                    end
                end
            end
        end
        return val
    end
    macro _nothingmacro(ex)
        ex
    end
    function _generate_QPs(symmetricQP,symmetricP,res,sistema,ϵ)
        minP=-ceil(ClassicalDicke.maximum_P_for_ϵ(sistema,ϵ)/res)*res
        minQ=-ceil(ClassicalDicke.maximum_Q_for_ϵ(sistema,ϵ)/res)*res

        lastQ=-minQ
        lastP=-minP
        if symmetricQP || symmetricP
            lastP=0
        end
        if symmetricQP
            lastQ=0
        end
        if (2/res)%1!=0.0
            error("res has to be an integer fraction of 2")
        end

        Ps=minP:res:lastP
        Qs=minQ:res:lastQ
        return Qs,Ps
    end
    """
    ```julia
    function matrix_QP∫∫dqdpδϵ(system::ClassicalDickeSystem;
        f::Function,
        ϵ::Real,
        res::Real=0.1,
        symmetricQP::Bool=false,
        symmetricP::Bool=symmetricQP,
        paralelize::Bool=(Distributed.myid()==1),
        showprogress::Bool = true,
        pbatch_size::Integer=Int(min(ceil((4/res)^2/Distributed.nprocs()/10),50)))
    ```
    Returns  a tuple `(Qs, Ps, mat)` where `Qs` and `Ps` are real vectors spanning
    all possible values of ``P`` and ``Q`` with some step `res`, and `mat` 
    is a matrix whose entries are given by
    `mat[i,j] = `[`∫∫dqdpδϵ(system; ϵ=ϵ, Q=Qs[i], P=P[j], p_res=res,kargs...)`](@ref ∫∫dqdpδϵ)

    # Arguments
    - `system` should be an instance of [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `ϵ` is a real number (see arguments for [`∫∫dqdpδϵ`](@ref)).
    - `f` should be function with signature `f([Q,q,P,p])` that returns values that
      may be added together (Numbers, Arrays etc...).
    - `res` determines the separation between the elements of `Qs` and `Ps`. It also
      determines `p_res` in the calls to [`∫∫dqdpδϵ`](@ref).
    - `symmetricQP` indicates that `f([Q,q,P,p]) = f([-Q,q,P,p]) = f([Q,q,-P,p])` for all `Q` and `P`. If `true`
      only the non-positive entries of `Q` and `P` will be computed, and the other entries will be mirrored. 
      Default is `false`.
    - `symmetricP` indicates that `f([Q,q,-P,p]) = f([Q,q,P,p])` for all `P`. If `true`, only the side with 
      non-positive `P` coordinate will be comptued, and the other entries will be mirrored. The default is to 
      be the same that `symmetricQP`.
    - `paralelize` indicates whether to use all available workers. Defaults to `true` if this function is called 
      from worker `1`, and `false` else.
    - `showprogress` turns the progress bar on/off. Default is `true`.
    - `pbatch_size` is the maximum number of batches to run in a single run in a single worker. The default value is automatically
      optimized depending on `res` and the number of workers.
    - `kargs...` are redirected to [`∫∫dqdpδϵ`](@ref).
    """
    function matrix_QP∫∫dqdpδϵ(system::ClassicalDickeSystem;
        f::Function,
        ϵ::Real,
        res::Real=0.1,
        symmetricQP::Bool=false,
        symmetricP::Bool=symmetricQP,
        paralelize::Bool=(Distributed.myid()==1),
        showprogress::Bool = true,
        pbatch_size::Integer=Int(min(ceil((4/res)^2/Distributed.nprocs()/10),50)),
        kargs...)
        
        Qs,Ps=_generate_QPs(symmetricQP,symmetricP,res,system,ϵ)
        QPs=[(Q,P) for P in Ps, Q in Qs]
        if showprogress #revolvemos la matriz para que el ETA del proceso sea más acertada
            is=1:length(Ps)*length(Qs)
            sh=shuffle(is)
            revsh=[findfirst(x->x==i,sh) for i in is]
            ind(i)=(div(i-1,length(Qs))+1,mod(i-1,length(Qs))+1)
            iind(i,j)=(i-1)*length(Qs)+j
            QPs=[QPs[ind(sh[iind(i,j)])...] for i in 1:length(Ps), j in 1:length(Qs)]
        end

        addbatchsize(f) =function(args...)
            f(batch_size=pbatch_size,args...)
        end
        mat= (showprogress ? (paralelize ? addbatchsize(ProgressMeter.progress_pmap) : ProgressMeter.map) : (paralelize ? addbatchsize(Distributed.pmap) : map))(QPs) do (Q,P)
            ∫∫dqdpδϵ(system;ϵ=ϵ,Q=Q,P=P,f=f,p_res=res,kargs...)
        end
        if showprogress #desrevolvemos
            mat=[mat[ind(revsh[iind(i,j)])...] for i in 1:length(Ps), j in 1:length(Qs)]
        end
        if symmetricQP || symmetricP
            mat=vcat(mat,mat[end-1:-1:1,1:end]) #hacemos espejo sobre P
            Ps=[Ps; -Ps[end-1:-1:1]]
        end
        if symmetricQP
            mat=hcat(mat,mat[1:end,end-1:-1:1]) #hacemos espejo sobre Q
            Qs=[Qs; -Qs[end-1:-1:1]]
        end
        return Qs,Ps,mat
    end
    function _Renyi_pow(m::Number,α::Real)
        if α<0
            error("α has to be a non negative number")
        end
        if m==0
            return m
        end
        if α==1
            return m*log(m)
        else
            return m^α
        end
    end
    #states puede ser una matriz de estados, o una matriz de altura 4, en cuyo caso se interpreta a los estados como coherentes centrados en Q,q,P,p (las filas)
    function _husimi_Renyi_powers(systemQ,pt,state,cache,averageallstates;tol=1e-4,averagingfunction=mean,statesAreCentersOfCoherentStates=size(state)[1]==4,αs::AbstractArray{<:Real,1})
        if statesAreCentersOfCoherentStates
            cache[1].=DickeBCE.HusimiOfCoherent.((systemQ,),eachcol(state),(pt,))
        else
            cache[1].=DickeBCE.Husimi(systemQ,pt,state;tol=tol,datacache=cache[1]) #asignamos porque puede pasar que cache[1] sea un numero,
                                                                                #en cuyo caso datacache no se usa
        end
        if averageallstates
            m=averagingfunction(cache[1])
        #    if ! isa(m, Array)
        #        m=[m]
        #    end
            return [[m];[_Renyi_pow.(m,α) for α in αs]]
        end
        for (i,α) in enumerate(αs)
            cache[i+1].=_Renyi_pow.(cache[1],α)
        end
        return cache
    end
    """
    ```julia
    function ∫∫dqdpδϵ_Husimi_Renyi_powers(
        systemQ::DickeBCE.QuantumDickeSystem;
        states::Union{Array{<:Number,1},Array{<:Number,2}},
        res::Real,
        Htol::Real=1e-4,
        nonvalue=NaN,
        averageallstates::Bool=false,
        averagingfunction::Function=mean,
        α::Union{AbstractArray{<:Real,1},Real}=2,
        kargs...)
    ```
    Calls matrix_QP∫∫dqdpδϵ PENDING

    """
    function ∫∫dqdpδϵ_Husimi_Renyi_powers(systemQ::DickeBCE.QuantumDickeSystem;
        states::Union{Array{<:Number,1},Array{<:Number,2}},
        res::Real,
        Htol::Real=1e-4,
        nonvalue=NaN,
        averageallstates::Bool=false,
        averagingfunction::Function=mean,
        α::Union{AbstractArray{<:Real,1},Real}=2,
        kargs...)
        
        if ! isa(α, AbstractArray)
            α=[α]
        end
        function gencache()
            if length(size(states))==1
                return [0.0im]
            else
                return zeros(Complex{Float64},size(states)[2])
            end
        end
        cache=[gencache() for i in 1:(length(α)+1)]
        Qs,Ps,mat= matrix_QP∫∫dqdpδϵ(;f=(pt -> _husimi_Renyi_powers(systemQ,pt,states,cache,averageallstates;tol=Htol,averagingfunction=averagingfunction,αs=α)),res=res,nonvalue=nonvalue,kargs...)
        matlist=[[m===nonvalue ? nonvalue : real.(m[i]) for m in mat] for i in  1:(length(α)+1)]
        return Qs,Ps,matlist
    end
    function loc_measure_and_proj_husimi_QP_matrix(systemQ::DickeBCE.QuantumDickeSystem;
        states::Union{Array{<:Number,1},Array{<:Number,2}},
        res::Real,
        Htol::Real=1e-4,
        nonvalue=NaN,
        averageallstates::Bool=false,
        averagingfunction::Function=mean,
        α::Union{AbstractArray{<:Real,1},Real}=2,
        matrix_powers::Union{AbstractArray{<:Real,1},Real}=1,
        _transposematrixoflistintolistofmatrices::Bool=true,
        kargs...)
        if ! isa(α, AbstractArray)
            α=[α]
        end
        if ! isa(matrix_powers, AbstractArray)
            matrix_powers=[matrix_powers]
        end
        if matrix_powers ⊈ α∪[1]
            error("matrix_powers should be a subset of αs ∪ {1}")
        end
        if length(size(states))==1
            resultlength=1
        else
            resultlength=size(states)[2]
            if averageallstates
                d=averagingfunction(zeros(resultlength))
                try
                    if length(size(d))>1
                        throw("")
                    end
                    if !isa(d,Number)
                        if any(x->!isa(x,Number),d)
                            throw("")
                        end
                    end
                catch
                    error("averagingfunction should return a scalar or a 1D array of scalars")
                end
                resultlength=length(d)
            end
        end
        Qs,Ps,matlist=  ∫∫dqdpδϵ_Husimi_Renyi_powers(systemQ;states=states,res=res,Htol=Htol,nonvalue=nonvalue,averageallstates=averageallstates,averagingfunction=averagingfunction,α=α,kargs...)
        L=loc_measure_from_matrices(matlist,nonvalue=nonvalue,α=α)
        if _transposematrixoflistintolistofmatrices
            matproj=[[Union{typeof(nonvalue),Float64}[m===nonvalue ? nonvalue : m[i] for m in matlist[findfirst(α->α==mp,[1;α])]] for mp in matrix_powers] for i in 1:resultlength]
            if length(matrix_powers)==1
                matproj=[m[1] for m in matproj]
            end
        else
            matproj=matlist[1]
        end
        if resultlength==1
            if length(α)==1
                L=L[1]
            else
                L=[l[1] for l in L]
            end
            matproj=matproj[1]
        end
        return L,(Qs,Ps,matproj)

    end
    function loc_measure_from_matrices(matlist;nonvalue=NaN,α::Union{AbstractArray{<:Real,1},Real}=2)
        if ! isa(α, AbstractArray)
            α=[α]
        end
        if length(matlist)!=length(α)+1
            error("Should pass n+1 matrices, where n is the length of α")
        end
        if length(matlist)==1
            return []
        end
        ∫1=0.0
        ∫Q=nothing
        ∫Qα=Any[nothing for a in α]
        for v in zip(matlist...)
            if all(x->x===nonvalue,v)
                continue
            elseif any(x->x===nonvalue,v)
                error("Inconsistent matrices")
            end
            if ∫Q===nothing
                    ∫Q=v[1]*0
                    for i in 1:length(α)
                        ∫Qα[i]=v[1]*0
                    end
           end
           ∫1+=2*pi
           ∫Q+=v[1]
           for i in 2:length(v)
               ∫Qα[i-1]+=v[i]
           end
        end
        res=[_loc_measure_from_integrals.(∫Qα[i],∫Q,∫1,α) for (i,α) in enumerate(α)]
        if length(α)==1
            res=res[1]
        end
        return res
    end
    function _loc_measure_from_integrals(∫Qα::Number,∫Q::Number,∫1::Number,α::Real)
        if α==1
            r= ∫Q*exp(-∫Qα/∫Q)
        elseif α==0
            r= ∫Qα
        elseif α>0
            r= (∫Q^(α/(α-1)))*(∫Qα^(1/(1-α)))
        else
            error("α must be non negative")
        end
        return real(r/∫1)
    end
    function loc_measure(;kargs...)
        return loc_measure_and_proj_husimi_QP_matrix(;_transposematrixoflistintolistofmatrices=false,kargs...)[1]
    end
    function proj_husimi_QP_matrix(;kargs...)
        return loc_measure_and_proj_husimi_QP_matrix(;α=Real[],kargs...)[2]
    end
    function energy_shell_average(system::ClassicalDickeSystem;
        ϵ::Real,
        f::Function,
        res::Real=0.01,
        symmetricQP::Bool=false,
        symmetricP::Bool=symmetricQP)
        
        t=0.0
        ft=nothing

        Qs,Ps=_generate_QPs(symmetricQP,symmetricP,res,system,ϵ)
        for Q in Qs, P in Ps
            v=∫∫dqdpδϵ(system,ϵ=ϵ,Q=Q,P=P,f=f,p_res=res,nonvalue=nothing)
            if !(v===nothing)
                add=2pi
                if Q!=0 && symmetricQP
                    add*=2
                    v*=2
                end
                if P!=0 && (symmetricP || symmetricQP)
                    add*=2
                    v*=2
                end
                t+=add
                if ft===nothing
                    ft=v
                else
                    ft+=v
                end
            end
        end
        return ft/t
    end

end

 

 
