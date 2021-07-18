module DickeBCE

export hamiltonian, husimi, husimi_of_coherent, QuantumDickeSystem, diagonalization, coherent_state,
    coherent_state!, coherent_overlap, eigenstate_parities, Wigner, evolve, WignerProjQP, WignerProjqp,
    factor_R_of_coherent_state, participation_ratio, random_state!, random_state, random_cₖ_generator,
    survival_probability,dimension
    import ..ClassicalSystems
    using LinearAlgebra
    using ProgressMeter
    using SparseArrays
    using Distributed
    import ..PhaseSpaces
    using GSL
    using WignerSymbols
    using Distributed
    using CSV
    using ..ClassicalDicke
    using Distributions
    import Random
    using Statistics
    using DataFrames,Tables
    
    """
    ```julia
    mutable struct QuantumDickeSystem
    ```
    This object represents the quantum Dicke model. It stores the parameters of
    the system, and it may be passed to multiple functions
    in this module. To generate it, use
    
    ```julia
    function QuantumDickeSystem(classical_system::ClassicalDicke.ClassicalDickeSystem;
        j::Real,
        Nmax::Union{Integer,Nothing}=nothing)
    ```
  
    # Arguments
    - `classical_system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `j` is the value of ``j``. It must be a positive half-integer.
    - `Nmax - 1` is the maximum excitation of the modified bosonic sector in the efficient coherent
      basis (see [Bastarrachea2014PSa](@cite), [Bastarrachea2014PSb](@cite)). `Nmax` must be a positive integer. 
      It may be omitted if there is a saved diagonalization. 
      In that case, a call to [`diagonalization`](@ref DickeBCE.diagonalization) will populate this value with the 
      greatest available in the saved cache.    
    """
    mutable struct QuantumDickeSystem
       j::Real
       Nmax::Union{Integer,Nothing}
       classical_system::ClassicalDicke.ClassicalDickeSystem
       function QuantumDickeSystem(classical_system::ClassicalDicke.ClassicalDickeSystem;
           j::Real,
           Nmax::Union{Integer,Nothing}=nothing)
           if !isinteger(2*j) || j < 0
               error("j must be a positive half-integer. Got j = $j.")
           end
           if Nmax!== nothing && Nmax <= 0
               error("Nmax must be a positive integer. Got Nmax = $Nmax.")
           end
           if isinteger(j)
               j= Int(j)
           end
           return new(j,Nmax,classical_system)
       end
    end
    """
    ```julia
    function QuantumDickeSystem(;ω₀,ω,γ,j,Nmax=nothing)
    ```
    Shorthand for 
    [`QuantumDickeSystem`](@ref)`(ClassicalDicke.ClassicalDickeSystem(ω₀=ω₀, ω=ω, γ=γ), j=j, Nmax=Nmax)`.
    """
    function QuantumDickeSystem(;ω₀::Real,ω::Real,γ::Real, j::Real,Nmax::Union{Integer,Nothing}=nothing) 
        return QuantumDickeSystem(ClassicalDicke.ClassicalDickeSystem(ω₀=ω₀,ω=ω,γ=γ), j= j, Nmax=Nmax)
    end
    function Base.show(io::IO, qds::QuantumDickeSystem)
        ω₀,ω,γ,j,Nmax=get_params(qds)
        print(io,"QuantumDickeSystem(ω₀ = $ω₀, ω = $ω, γ = $γ, j = $j, Nmax = $Nmax)")
    end
    """
    ```julia
    function dimension(system::QuantumDickeSystem)
    ```
    Returns the integer ``(2j + 1)\\times N_\\text{max}``, which gives the dimension
    of the Hilbert space.
    """
    dimension(system::QuantumDickeSystem) = Int(2*system.j +1)*system.Nmax
    ind(n,m,N,j)= Int(n)*Int(N+1)+Int(m+j) +1 #+1 because Julia counts from 1.
    function logfact(n::Real)::Float64
        if n==0
            return 0.0
        end
        return sum(log(i) for i in 1:n)
    end
    function _sqrt_factorial(n::Real)
        if n==0
            return 1
        end
        prod(sqrt(i) for i in 1:n)
    end
    function _sqrt_binomial(n::Real,k::Real)
        return exp(_sqrt_binomiallog(n,k))
    end
    function _sqrt_binomiallog(n::Real,k::Real)
        return (logfact(n)-logfact(k)-logfact(n-k))/2
    end
    function get_params(system::QuantumDickeSystem;warning_on_NmaxNothing=true)
        j=system.j
        Nmax=system.Nmax 
        ω₀,ω,γ =ClassicalSystems.parameters(system.classical_system)
        if Nmax!==nothing
            Nmax=Int(Nmax)
        else
            if warning_on_NmaxNothing
                @warn "You did not pass Nmax and did not call diagonalization to load it from disk. Nmax is undetermined and this probably will cause an error below."
            end
        end
        return ω₀,ω,γ,j,Nmax
    end
    """
    ```julia
    function parity_operator(system::QuantumDickeSystem)
    ```
    Returns a sparse matrix corresponding to the parity operator ``\\hat{\\Pi}=e^{i\\pi(\\hat{a}^\\dagger\\hat{a}+ \\hat{J}_z + j)}`` in the
    efficient coherent basis.
    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    """
    function parity_operator(system::QuantumDickeSystem)
        ω₀,ω,γ,j,Nmax=get_params(system)
        N=Int(2*j)
        I=Int64[]
        J=Int64[]
        V=Int8[]
        function pushv(i,j,v)
            push!(I,i)
            push!(J,j)
            push!(V,v)
        end
        
        for indexN in 0:(Nmax-1), m in -j:j
            pushv(ind(indexN,m,N,j),ind(indexN,-m,N,j),(-1)^indexN)
        end
        return sparse(I,J,V,Nmax*(N+1), Nmax*(N+1))
    end
    """
    ```julia
    function eigenstate_parities(system::QuantumDickeSystem,eigenstates)
    ```
    Returns a vector of `-1`s and `1`s contaning the parities of all of the eigenstates.
    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    - `eigenstates` is the matrix of eigenstates. (See  [`diagonalization`](@ref DickeBCE.diagonalization)).
    """
    function eigenstate_parities(system::QuantumDickeSystem,eigenstates)
        Π=parity_operator(system)
        ΠE=Π*eigenstates
        return Int8[Integer(round(dot(view(ΠE,:,k),view(eigenstates,:,k)))) for k in 1:size(eigenstates)[2]]
    end
    function population_of_level_n(system,state,n)
        ω₀,ω,γ,j,Nmax=get_params(system)
        N=Int(2*j)
        return sum(abs(state[ind(n,m,N,j)]) for m in -j:j)
    end
    """
    ```julia
    function eigenenergies(system::QuantumDickeSystem;kargs...)
    ```
    Returns a vector containing the eigenenergies of the system. 
    
    This function is just shorthand for [`diagonalization`](@ref DickeBCE.diagonalization)`(system;only_eigenenergies = true, kargs...)`. 
    """
    function eigenenergies(system::QuantumDickeSystem;kargs...)
        return diagonalization(system;only_eigenenergies=true,kargs...)
    end
    """
    ```julia
    function diagonalization(system::QuantumDickeSystem;
            load_cache = true,
            save_cache = true,
            cache_folder = joinpath(homedir(),"dicke_diagonalizations"),
            maxϵ::Real = 5.0,
            onlyload::Union{AbstractVector{<:Integer},Nothing} = nothing,
            only_eigenenergies = false,
            verbose::Bool = true,
            converged_tolerance=1e-3)
    ```
    Diagonalizes the Dicke Hamiltonian up to a maximum energy `maxϵ`. The resulting eigenstates 
    are guaranteed to be converged, with a tolerance determined by `converged_tolerance`.
    Numerical degeneracies are also corrected, to ensure that the eigenstates have parity +1 or -1.
    If `only_eigenenergies` is `false` (default), a tuple `(eigenenergies,eigenstates)` is returned, where
    `eigenenergies` is real vector containing the eigenenergies and `eigenstates` is a real matrix that
    contains the eigenstates as columns. 
    
    If `load_cache = true` it will try to load saved diagonalizations from `cache_folder`.


    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    - `load_cache` is a boolean indicating whether to try to load from cache folder. Defaults to `true`.
    - `save_cache` determines if the results of a computed diagonalization are saved to the cache folder. Defaults to `true`.
    - `cache_folder` is the cache folder where diagonalizations are saved. Default is `%HOME%/dicke_diagonalizations`
    - `maxϵ` is the maximum energy up to which we diagonalize. Keep this number 
      higher than the maximum converged regime you want.
    - `onlyload` may be a vector of integers, indicating the indices of the eigenstates
      to load, or the default value, `nothing`, indicates that all eigenstates should be loaded.
    - `only_eigenenergies` should be `true` if you only to load the eigenenergies. Defaults
      to `false` (return eigenstates and eigenenergies).
    - `verbose` is `true` (default) if you want to see info messages.
    - `converged_tolerance` determines how strict we are in saying that an eigenstate is converged (see [Bastarrachea2014PSb](@cite)).
      The default value `1e-3` is usually the best.
    """
    function diagonalization(system::QuantumDickeSystem;
            load_cache = true,
            save_cache = true,
            cache_folder = joinpath(homedir(),"dicke_diagonalizations"),
            maxϵ::Real = 5.0,
            onlyload::Union{AbstractVector{<:Integer},Nothing} = nothing,
            only_eigenenergies = false,
            verbose::Bool = true,
            converged_tolerance=1e-3)
        
        ω₀,ω,γ,j,Nmax=get_params(system,warning_on_NmaxNothing=false)
        
        filename(ω₀,ω,γ,j,Nmax)="dicke_bce_ω₀=$(ω₀),ω=$(ω),γ=$(γ),j=$(j),Nmax=$(Nmax)"
        filelocation(ω₀,ω,γ,j,Nmax)=joinpath(cache_folder,filename(ω₀,ω,γ,j,Nmax))
        if Nmax===nothing
            if save_cache==false
                error("Nmax is nothing,  so we can't diagonalize, but save_cache is false. Give Neff or set save_cache to true")
            end
            reg=Regex(string(filename(ω₀,ω,γ,j,"(\\d+)"),"_E.csv"))
            diags=readdir(cache_folder)
            for d in diags
                m=match(reg,d)
                if m===nothing
                    continue
                end
                poss_nmax=parse(Int64,m.captures[1])
                flname=filelocation(ω₀,ω,γ,j,poss_nmax)
                if isfile("$(flname)_E.csv") && isfile("$(flname)_V.csv")
                    if Nmax===nothing
                        Nmax=poss_nmax
                    else
                        Nmax=max(poss_nmax,Nmax)
                    end
                end
            end
            if Nmax===nothing
                error("No diagonalizations found, and Nmax is nothing. Pass a value for Nmax to diagonalize.")
            else
                system.Nmax=Nmax
            end
        end

        flname=filelocation(ω₀,ω,γ,j,Nmax)
        if load_cache

             if isfile("$(flname)_E.csv") && isfile("$(flname)_V.csv")
                if verbose
                @info "Loading diagonalization: Nmax=$Nmax"
                end
                eigenenergies=DataFrame(CSV.File("$(flname)_E.csv")).E[:]
                if only_eigenenergies
                    return eigenenergies
                end
                eigenstates=Tables.matrix(DataFrame(CSV.File("$(flname)_V.csv",select=onlyload,type=Float64,threaded =true,tasks=4)))
                return eigenenergies,eigenstates
             end
        end


        Hq=hamiltonian(system)

        if verbose
            @info "Diagonalizing..."
        end
        FullH=Symmetric(Matrix(Hq));
        Emin=-10*j
        Emax=maxϵ*j #maxenergy
        eigenenergies,eigenstates=eigen(FullH,Real(Emin),Real(Emax));

        convstates=(length(eigenenergies)-1)
        for k in 1:(length(eigenenergies)-1)
            if population_of_level_n(system,eigenstates[:,k],Nmax-1) > converged_tolerance
                convstates=k
                break
            end
        end
        eigenstates=eigenstates[:,1:convstates]
        eigenenergies=eigenenergies[1:convstates]
        if verbose
        @info "$convstates converged states were obtained up to ϵ=$(eigenenergies[end]/j)."

        @info "Fixing numerical degeneracies (correcting parity)."
        end
        #corregimos por degeneraciones numéricas
        Π=parity_operator(system)
        for k in 1:(length(eigenenergies)-1)
            if(eigenenergies[k+1]-eigenenergies[k]<1e-5)
                _correct_eigenstates(eigenstates,k,k+1,Π)
            end
        end
        if save_cache
            mkpath(cache_folder)
            CSV.write("$(flname)_V.csv",DataFrame(eigenstates,:auto))
            CSV.write("$(flname)_E.csv",DataFrame(E=eigenenergies))
        end
        return eigenenergies,eigenstates
    end
    function _correct_eigenstates(eigenstates,k1,k2,Π)
        v1=eigenstates[:,k1]
        v2=eigenstates[:,k2]
        eigenvals,mat =eigen(Transpose([v1 v2])*Π*[v1 v2])
        if !issetequal([1.0,-1.0],round.(eigenvals,digits=7))
            return
        end
        λ_1=mat[1,1]
        λ_2=mat[2,1]
        if λ_1<0
            λ_1,λ_2=-λ_1,-λ_2
        end
        μ_1=mat[1,2]
        μ_2=mat[2,2]
        if μ_2<0
            μ_1,μ_2=-μ_1,-μ_2
        end
        eigenstates[:,k1].= λ_1*v1+λ_2*v2;
        eigenstates[:,k2].=μ_1*v1+μ_2*v2;
        return true
    end
    """
    ```julia
    function hamiltonian(system::QuantumDickeSystem)
    ```
    Returns a sparse matrix corresponding to the Dicke Hamiltonian (See Eq. (1) of [Pilatowsky2021NatCommun](@cite)) in the
    efficient coherent basis. (See Ref. [Bastarrachea2014PSa](@cite), [Bastarrachea2014PSb](@cite))
    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    """
    function hamiltonian(system::QuantumDickeSystem;_ignorediagonalterms=false)
        ω₀,ω,γ,j,Nmax=get_params(system)
        N=Int(2*j)
        G=2*γ/(ω*sqrt(N))
        function DNNp(N,Np)

            f(k)=exp((logfact(N)+logfact(Np))/2 - (logfact(Np-k)+logfact(N-k)+logfact(k)))
            return exp(-G^2/2)*sum((f(k)*(-1)^k) * (G^(N+Np-2*k)) for k in 0:Np)
        end
        I=Int64[]
        J=Int64[]
        V=Float64[]
        function pushv(i,j,v)
            push!(I,i)
            push!(J,j)
            push!(V,v)
        end
        for indexN in 0:(Nmax-1)
            for indexNprime in 0:indexN
                D=DNNp(indexN,indexNprime)
                for m in -j:j
                    if !_ignorediagonalterms
                        if indexNprime==indexN
                            pushv(ind(indexN,m,N,j),ind(indexN,m,N,j),ω*(indexN-G^2*m^2))
                        end
                    end
                    if m<j
                        o=(-1)^indexN*D
                        v=-ω₀*sqrt(j*(j+1)-m*(m+1))*o/2 
                        pushv(ind(indexNprime,m+1,N,j),ind(indexN,m,N,j),v)

                    end
                    if m>-j
                        o=(-1)^indexNprime*D
                        v=-ω₀*sqrt(j*(j+1)-m*(m-1))*o/2 
                        pushv(ind(indexNprime,m-1,N,j),ind(indexN,m,N,j),v)

                    end
                end
            end
        end
        D=dimension(system)
        H = sparse(I,J,V, D,D)
        H=sparse(Symmetric(H));
        return H
    end
    """
    ```julia
    function Jz(system::QuantumDickeSystem)
    ```
    Returns a sparse matrix representing the operator ``\\hat{J}_z`` in the
    efficient coherent basis.
    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    """
    function Jz(system::QuantumDickeSystem)
        ω₀,ω,γ,j,Nmax=get_params(system)
        return hamiltonian(QuantumDickeSystem(ClassicalDicke.ClassicalDickeSystem(ω₀=1,ω=1,γ=1),j=j,Nmax=Nmax),_ignorediagonalterms=true)
    end
    """
    ```julia
    function Jx(system::QuantumDickeSystem)
    ```
    Returns a diagonal matrix representing the operator ``\\hat{J}_x`` in the
    efficient coherent basis.
    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    """
    function Jx(system::QuantumDickeSystem)
        ω₀,ω,γ,j,Nmax=get_params(system)
        N=Int(2j)
        Jx = zeros(Float64,dimension(system))
        for indexN in 0:(Nmax-1)
            for m in -j:j
               Jx[ind(indexN,m,N,j)]=m
            end
        end
        return Diagonal(Jx)
    end
    """
    ```julia
    function participation_ratio(state::AbstractVector{<:Number};
        eigenstates::AbstractMatrix{<:Number},
        eigenenergies::Union{AbstractVector{<:Real},Nothing}=nothing,
        count_degeneracies::Bool=eigenenergies!=nothing,
        degentol::Real=1e-5)
    ```
    Returns the Participation Ratio (PR) of `state` in the eigenbasis. If `count_degeneracies` is `true`,
    then Eq. (7) of [Villasenor2020](@cite) is used. If `count_degeneracies` is `false`, then the inverse
    of Eq. (19) of [Villasenor2020](@cite) is used.
    # Arguments
    - `state` should be a complex vector representing the state in the efficient coherent basis.
    - `eigenstates` should be a matrix containing the eigenstates.
    - `eigenenergies` should be passed if `count_degeneracies` is  `true`. It is a list
      containing the eigenenergies.
    - `count_degeneracies` -- whether to modify the PR definition to account for degeneracies (see above). 
      Default is `false` if `eigenenergies` is not passed, else it is `true`.
    - `degentol` minimum energy separation below which two eigenstates are considered degenerate. Default
      is `1e-5`.
    """
    function participation_ratio(state::AbstractVector{<:Number};
        eigenstates::AbstractMatrix{<:Number},
        eigenenergies::Union{AbstractVector{<:Real},Nothing}=nothing,
        count_degeneracies::Bool=eigenenergies!=nothing,
        degentol::Real=1e-5)
        eigenbasis_state=transpose(eigenstates)*state
        if count_degeneracies
            IPR=0.0
            degsbspace=0.0
            for i in 1:length(eigenenergies)

                if  i!=1 && eigenenergies[i]-eigenenergies[i-1]>degentol
                    IPR+=degsbspace^2
                    degsbspace=0.0
                end
                degsbspace+=abs2(eigenbasis_state[i])
                if i==length(eigenenergies)
                    IPR+=degsbspace^2
                end
            end
            return 1/IPR
        else
            return 1/sum(abs2(c)^2 for c in eigenbasis_state)
        end
    end
    """
    ```julia
    function factor_R_of_coherent_state(system::QuantumDickeSystem,
            x::AbstractVector{<:Real};
            eigenstates::AbstractMatrix{<:Number},
            eigenenergies::AbstractVector{<:Real},
            state::AbstractVector{<:Number} = coherent_state(system,x))
    ```
    Computes ``R`` for a coherent state, as defined in Eq. (30) of [Villasenor2020](@cite).
    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    - `x` is an array `[Q,q,P,p]` representing the center of the coherent state.
    - `eigenstates` should be a matrix containing the eigenstates as columns.
    - `eigenenergies`  is a list containing the eigenenergies.
    - `state` is a complex vector representing the coherent state ``\\left | \\mathbf{x} \\right \\rangle``.
      If it is not passed, it is computed using [`DickeBCE.coherent_state`](@ref)
    """
    function factor_R_of_coherent_state(system::QuantumDickeSystem,
            x::AbstractVector{<:Real};
            eigenstates::AbstractMatrix{<:Number},
            eigenenergies::AbstractVector{<:Real},
            state::AbstractVector{<:Number} = coherent_state(system,x))
        systemC = system.classical_system
        ω₀,ω,γ,j,Nmax=get_params(system)
        ϵ=ClassicalDicke.hamiltonian(systemC)(x)
        σ=ClassicalDicke.energy_width_of_coherent_state(systemC, x; j)
        ν(ϵ)=ClassicalDicke.density_of_states(systemC, ϵ; j)
        n=Distributions.Normal(ϵ,σ)
        f=1/(2*sqrt(pi)*σ*j*ν(ϵ))
        
        return 2*participation_ratio(state,eigenstates=eigenstates,eigenenergies=eigenenergies,count_degeneracies=true)*f 
    end
    """
    ```julia
    function survival_probability(t::Union{Real,AbstractArray{<:Real}};
        state::AbstractVector{<:Number},
        eigenstates::AbstractMatrix{<:Number},
        eigenenergies::AbstractVector{<:Real},
        normwarning::Real=0.99)
    ```
    Computes the survival probability ``| \\langle`` `state` ``| e^{-i \\hat{H} t} |`` `state` ``\\rangle | ^2``
    at the times given.
    # Arguments
    - `t` may be a single time, or a vector of times. In the first case, a number will be 
      returned, and in the second, a vector.
    # Keyword arguments
    - `state` is a vector representing the state in the BCE.
    - `eigenstates` and `eigenenergies` are the ones returned by [`diagonalization`](@ref).
    - `normwarning` is a tolerance (defaults to `0.99`). If the norm of the state
      in the eigenbasis is below this number, a warning will be thrown. This usually
      happens if `Nmax` is too small for the energy regime you are working on.
    """
    function survival_probability(t::Union{Real,AbstractArray{<:Real}};
        state::AbstractVector{<:Number},
        eigenstates::AbstractMatrix{<:Number},
        eigenenergies::AbstractVector{<:Real},
        normwarning::Real=0.99)
        if typeof(t)<:Real
            return survival_probability([t],state=state,eigenstates=eigenstates,eigenenergies=eigenenergies)[1]
        end
        eigenbasis_state=transpose(eigenstates)*state
        if norm(eigenbasis_state)<normwarning
            @warn "The state is not converged in the eigenbasis, this is probably because Nmax is too small"
        end
        return [abs2(sum(exp(-im*ti*Ek)*abs2(ck) for (ck,Ek) in zip(eigenbasis_state,eigenenergies))) for ti in t]
     end
     """
     ```julia
     function evolve(t::Union{Real,AbstractArray{<:Real}}
         state::AbstractVector{<:Number};
         eigenstates::AbstractMatrix{<:Number},
         eigenenergies::AbstractVector{<:Real},
         normwarning::Real=0.99)
     ```
     Computes the evolution ``e^{-i \\hat{H} t}|`` `state` ``\\rangle`` under the
     Dicke Hamiltonian ``\\hat{H}``. The result is returned in the BCE.
     
     # Arguments
     - `t` may be a single time, or a vector of times. In the first case, a vector will be 
       returned, and in the second, a matrix with each column corresponding to each time.
     - `state` is a vector representing the state in the BCE.
     # Keyword arguments
     - `eigenstates` and `eigenenergies` are the ones returned by [`diagonalization`](@ref).
     - `normwarning` is a tolerance (defaults to 0.99). If the norm of the state
       in the eigenbasis is below this number, a warning will be thrown. This usually
       happens if `Nmax` is too small for the energy regime you are working on.
     """
     function evolve(t::Union{Real,AbstractArray{<:Real}},
         state::AbstractVector{<:Number};
         eigenstates::AbstractMatrix{<:Number},
         eigenenergies::AbstractVector{<:Real},
         normwarning::Real=0.99)
        if typeof(t)<:Real
            return evolve([t],state,eigenstates=eigenstates,eigenenergies=eigenenergies)
        end
        eigenbasis_state=transpose(eigenstates)*state
        if norm(eigenbasis_state)<normwarning
            @warn "The state is not converged in the eigenbasis, this is probably because Nmax is too small"
        end
        eigenstates*hcat((Diagonal(exp.(-im*ti*eigenenergies))*eigenbasis_state for ti in t)...)
    end
    """
    ```julia
    function coherent_state(system::QuantumDickeSystem,
        x::AbstractVector{<:Real};
        normwarning::Real=0.99,
        extra_phase::Complex=1.0+0im,
        tol::Real=0.0)
    ```
    Calls [`coherent_state!`](@ref) passing `data` as a new vector and returns the result.
    """
    function coherent_state(system::QuantumDickeSystem,
        x::AbstractVector{<:Real};
        normwarning::Real=0.99,
        extra_phase::Complex=1.0+0im,
        tol::Real=0.0)
        
        ω₀,ω,γ,j,Nmax=get_params(system)

        data = zeros(Complex{Float64},Nmax*(Int(2*j)+1))
        coherent_state!(system,x,data;normwarning=normwarning,extra_phase=extra_phase,tol=tol)
        return data
    end
    """
    ```julia
    function coherent_state!(system::QuantumDickeSystem,
                x::AbstractArray{<:Real,1},
                data::AbstractVector{<:Number};
                normwarning=0.99,
                add_state_to_data::Bool=false,
                extra_phase::Complex=1.0+0im,
                tol::Real=1e-6)
    ```
    This function computes the coefficients of a coherent state centered at 
    `x` in the BCE, and stores the result in `data`. (See Ref. [Villasenor2020](@cite))
    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    - `x` is a vector with the coordinates `[Q, q, P, p]`.
    - `data` should be a complex vector of length [`dimension`](@ref)`(system)`.
    # Keyword arguments
    - `normwarning` is a tolerance (defaults to 0.99). If the norm of the resulting
      coherent state is below this number, a warning will be thrown. This usually
      happens if Nmax is too small for the energy regime you are working on.
    - `add_state_to_data` is a boolean. If it is `true`, the coefficients will be added
      to `data`. If it is `false` (default), this function will override any preexisting
      values in `data`. This is useful if you want to add many coherent states together
      without having to allocate that much memory.   
    - `extra_phase` is a complex number that multiplies the resulting state overall. 
      Defaults to `1`.
    - `tol` is a numerical tolerance. If `tol=0`, then all the coefficients all computed.
      However, if `ε = 1 -  tol > 0.0`, then some coefficients at the tail of the distribution (whose total squared norm does not exceed `ε`) will be treated as `0.0`. This allows to 
      significantly reduce computation time, but introduces a numerical error of order `ε`. The default is `1e-6`.
      See [this example](@ref exampletolhusimis) and Ref. [Pilatowsky2020Notes](@cite).
    """
    coherent_state!(system::QuantumDickeSystem,
                x::AbstractArray{<:Real,1},
                data::AbstractVector{<:Complex};
                normwarning=0.99,
                add_state_to_data::Bool=false,
                extra_phase::Complex=1.0+0im,
                tol::Real=1e-6)=
        _coherent_overlap_or_coherent_state(system,x;
                        datacache=data,tol=tol,normwarning=normwarning,add_state_to_data=add_state_to_data,extra_phase=extra_phase)

    function _Cohpart1(n,m,alpha,G)
        s=alpha+G*m
        if n==0 && s==0.0
            return 0.0im #para evitar 0*log(0)
        end
        return n*log(s) - logfact(n)/2
    end
    function _Cohpart2(m,N,j,w,alpha,G)
        a= _sqrt_binomiallog(N,Int(j + m)) + log(w)*m -abs2(G*m)/2 -G*m*(alpha)
        return a
    end

    function qtlmin(dist,tol,mn)
        vl=mn
        try
            dt=Distributions.quantile(dist,tol)
            vl=max(vl,dt)
        catch
        end
        return vl
    end
    function qtlmax(dist,tol,mx)
        vl=mx
        try
            dt=Distributions.quantile(dist,tol)
            vl=min(vl,dt)
        catch
        end
        return vl
    end
    """
    ```julia
    function husimi_of_coherent(system::QuantumDickeSystem,
        x::AbstractArray{<:Real,1},
        y::AbstractArray{<:Real,1})
    ```
    Returns  ``\\left | \\left \\langle \\mathbf{x}\\middle | \\mathbf{y}\\right \\rangle \\right | ^2``. 
    Equation (3.14b) of Ref. [Arecchi1972](@cite) is used for the overlap of the Bloch coherent states.
    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    - `x` and `y` are vectors in the form `[Q, q, P, p]`.
    """
    function husimi_of_coherent(system::QuantumDickeSystem,x::AbstractArray{<:Real,1},y::AbstractArray{<:Real,1})
        ω₀,ω,γ,j,Nmax=get_params(system)
        Q1,q1,P1,p1=x
        Q2,q2,P2,p2=y
        return exp(-(j/2)*((q1-q2)^2 + (p1-p2)^2))*cos(PhaseSpaces.arc_between_QP(Q1,P1,Q2,P2)/2)^(4j)
    end
    onelistcache=[0.0im]
    """
    ```julia
    function coherent_overlap(system::QuantumDickeSystem,
        x::AbstractVector{<:Real},
        state::AbstractVector{<:Number};
        tol::Real=1e-6,
        normwarning::Real=0.99,
        datacache::Union{AbstractArray{Complex{Float64},1},Nothing}=nothing)
    ```
    Returns  the overlap ``\\langle \\mathbf{x} | `` `state` ``\\rangle``,
    where ``\\left | \\mathbf{x} \\right \\rangle`` is a coherent state centered at `x`.
    # Arguments
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    - `x` is a vector in the form `[Q, q, P, p]`.
    - `state` is a vector in the BCE.
    # Keyword arguments
    - `tol` is a numerical tolerance. If `tol=0`, then all the products of coefficients are computed.
      However, if `ε = 1 -  tol > 0.0`, then some coefficients of the coherent state at the tail of the distribution  (whose total squared norm does not exceed `ε`)
      will be treated as `0.0`. This allows to 
      significantly reduce computation time, but introduces a numerical error of order `ε`. The default is `1e-6`.
      See [this example](@ref exampletolhusimis)  and Ref. [Pilatowsky2020Notes](@cite).
    - `normwarning` is a tolerance (defaults to 0.99). If the norm of the coherent  state is below this number, a warning will be thrown. 
      This usually happens if Nmax is too small for the energy regime you are working on.
    - `datacache` in an array where to store the result, or `nothing` for the result to be returned. 
    """
    function coherent_overlap(system::QuantumDickeSystem,
        x::AbstractVector{<:Real},
        state::AbstractVector{<:Number};
        tol::Real=1e-6,
        normwarning::Real=0.99,
        datacache::Union{AbstractVector{<:Complex},Nothing}=nothing)
        if datacache===nothing
            datacache=onelistcache
        end
        return _coherent_overlap(system,x,state;datalength=1,datacache=datacache,tol=tol,normwarning=normwarning)[1]
    end
    """
    ```julia
    function coherent_overlap(system::QuantumDickeSystem,
        x::AbstractVector{<:Real},
        states::AbstractMatrix{<:Number};
        tol::Real=1e-6,
        normwarning::Real=0.99,
        datacache::Union{AbstractVector{<:Complex},Nothing}=nothing)
    ```
    Same as [`coherent_overlap(..., state,...)`](@ref coherent_overlap(::QuantumDickeSystem,::AbstractVector{<:Real},::AbstractVector{<:Number})), but
    instead of one `state`, it allows `states` to be a matrix with multiple states as columns. It is functionally equivalent to
    ```julia
    [coherent_overlap(..., states[:,k], ...) for k in 1:size(states)[2]]
    ```
    but it is (much) faster.
    The result is stored in `datacache` if provided.
    """
    function coherent_overlap(system::QuantumDickeSystem,
        x::AbstractVector{<:Real},
        states::AbstractMatrix{<:Number};
        tol::Real=1e-6,
        normwarning::Real=0.99,
        datacache::Union{AbstractVector{<:Complex},Nothing}=nothing)
        d=size(states)[2]
        if datacache===nothing
            datacache=zeros(Complex{Float64},d)
        end
        return _coherent_overlap(system,x,states;datacache=datacache,datalength=d,tol=tol,normwarning=normwarning)
    end
    #if you pass estados it calculates overlap with coherent, else  coefficients of coherent.
    function _coherent_overlap_or_coherent_state(system::QuantumDickeSystem,punto::AbstractArray{<:Real,1};
                        datacache,datalength=nothing,tol=1e-6,normwarning=0.99,
                        estados::Union{AbstractArray{<:Number,1},AbstractArray{<:Number,2},Nothing}=nothing,add_state_to_data=false,extra_phase::Complex=1.0+0im)
        #@fastmath begin
            if !add_state_to_data
                datacache.=0.0im
            end

            ω₀,ω,γ,j,Nmax=get_params(system)
            N=Int(2*j)
            G=2*γ/(ω*sqrt(N))
            Q,q,P,p=punto
            theta=PhaseSpaces.θ_of_QP(Q,P)
            phi=PhaseSpaces.φ_of_QP(Q,P)

            z=exp(-phi*im)*tan(theta/2)
            alpha=(q+im*p)*sqrt(j/2)
            m=-j


            nrm=0.0
            w=(1+z)/(1-z)
            part0=-abs2(alpha)/2 +j*log(w/(1+abs2(w)))

            if tol==0
               mmin=-j
               mmax=j
            else
                atomicDist=Distributions.Binomial(N,abs2(w)/(1+abs2(w)))
                mmin=qtlmin(atomicDist,tol/4,0)-j
                mmax=qtlmax(atomicDist,1-tol/4,N)-j
            end
            for m in mmin:mmax
            #for m in -j:j
                part2m=part0+_Cohpart2(m,N,j,w,alpha,G)
                if tol==0
                    nmin=0
                    nmax=Nmax-1
                else
                    efftol=tol/((mmax-mmin+1) *Distributions.pdf(atomicDist,m+j))
                    d=Distributions.Poisson(abs(G*m+alpha)^2)
                    nmin= qtlmin(d,efftol/4,0)
                    nmax=qtlmax(d,1-efftol/4,Nmax-1)
                end
                for n in nmin:nmax
                    v=conj(exp(_Cohpart1(n,m,alpha,G)+part2m))
                    if estados!==nothing
                        for e in 1:datalength
                            datacache[e]+=estados[ind(n,m,N,j),e]*v
                        end

                    else
                       datacache[ind(n,m,N,j)]+=conj(v)*extra_phase
                    end
                    nrm+=abs2(v)
                end
            end
            if nrm<normwarning
                @warn "The coherent state is not converged (norm=$nrm). You are probably working in an energy regime that is too high. You should make Nmax bigger."
            end
        #end
        return datacache
    end
    _coherent_overlap(system::QuantumDickeSystem,
    punto::AbstractArray{<:Real,1},
    estados::Union{AbstractArray{<:Number,1},AbstractArray{<:Number,2}};
    datacache,
    datalength,
    tol=1e-6,
    normwarning=0.99) =_coherent_overlap_or_coherent_state(system,punto;
                                                        datacache=datacache,
                                                        datalength=datalength,
                                                        tol=tol,
                                                        estados=estados,
                                                        normwarning=normwarning)
    """
    ```julia
    function husimi(args...;kargs...)
    ```
    Computes `abs2.(`[`coherent_overlap`](@ref)`(args...;kargs...))`. The arguments
    and behavior are the same as [`coherent_overlap`](@ref).
    """
    function husimi(args...;kargs...)
        r=coherent_overlap(args...;kargs...)
        return abs2.(r)
    end

    #esto no funciona para j grande, me cambie a WignerSymbols
    function clebsch(j1, j2, j3, m1, m2, m3)

        if m3 != m1 + m2
            return 0
        end
        vmin = Int(max(-j1 + j2 + m3, -j1 + m1, 0))
        vmax = Int(min(j2 + j3 + m1, j3 - j1 + j2, j3 + m3))

        C = exp((log(2.0 * j3 + 1.0) + logfact(j3 + j1 - j2) +
                    logfact(j3 - j1 + j2) + logfact(j1 + j2 - j3) +
                    logfact(j3 + m3) + logfact(j3 - m3) -
                    (logfact(j1 + j2 + j3 + 1) +
                    logfact(j1 - m1) + logfact(j1 + m1) +
                    logfact(j2 - m2) + logfact(j2 + m2)))/2)
        S = 0
        for v in vmin:vmax
            S += (-1.0) ^ (v + j2 + m2)* exp( -logfact(v) +
                logfact(j2 + j3 + m1 - v) + logfact(j1 - m1 + v) -
                logfact(j3 - j1 + j2 - v) - logfact(j3 + m3 - v) -
                logfact(v + j1 - j2 - m3))
        end
        C = C * S
        return C
    end
    function sph_harm(l,m,theta,phi;legendre_sphPlm_cache=:nothing)
        fact=1
        if m<0
            fact=(-1)^(-m)
        end
        d=0
        if legendre_sphPlm_cache==:nothing
            d=GSL.sf_legendre_sphPlm(l,abs(m),cos(theta))
        else
            d=legendre_sphPlm_cache[l,abs(m),theta]
        end
        return fact*d*exp(im*m*phi)
    end
    function generate_hermite_cache(n_max,pts)
        cache=Dict{NTuple{2,Real},Complex}()
        points=unique(pts)
        for x in points

            data=GSL.sf_hermite_phys_array(n_max,x)
            for (n,d) in enumerate(data)
                cache[n-1,x]=d
            end
        end
        return cache
    end
    function generate_legendre_sphPlm_cache(l_max,pts)
        cache=Dict{NTuple{3,Real},Complex}()
        thetas=unique(x[3] for x in pts)
        for theta in thetas

            data=GSL.sf_legendre_array_e(GSL.GSL_SF_LEGENDRE_SPHARM,l_max,cos(theta),-1)
            for l in 0:l_max,m in 0:l
                cache[l,m,theta]=data[GSL.sf_legendre_array_index(l, m)+1]
            end
        end
        return cache
    end
    function laguerre(n,a,x)

        if a<0
            return x^(-a)/GSL.sf_poch(-n,-a)*laguerre(n+a,-a,x)
        end
        return GSL.sf_laguerre_n(n,a,x)
    end
    function RotZtoX(θ,φ)
        a=-π/2
        z=[cos(θ/2) -exp(-im*φ)*sin(θ/2)]*[cos(a/2) -sin(a/2);
                                        sin(a/2) cos(a/2)]
        return 2*atan(abs(z[2]),abs(z[1])),angle(z[2])-angle(z[1])
    end
    function classical_QqPp_to_quantum_rot_qpθφ(j,QqPp)
        Q,q,P,p=QqPp
        q=sqrt(j)*q
        p=sqrt(j)*p
        θ,φ=RotZtoX(PhaseSpaces.θ_of_QP(Q,P),PhaseSpaces.φ_of_QP(Q,P))
        return [q,p,θ,φ]
    end

    function sph_harm_pt(k,i,pt;legendre_sphPlm_cache=:nothing)
        return sph_harm(k, i, pt[3],pt[4];legendre_sphPlm_cache=legendre_sphPlm_cache)
    end

    function Wigner_Fock(n,np,q,p)
        two_r_squared=(q^2+p^2)*2
        if two_r_squared==0
            two_r_squared=1e-30
        end
        return (exp((logfact(np)-logfact(n))/2 +
            im*(np-n)*atan(p,q)- two_r_squared/2)*(-1)^(np)/π * two_r_squared^((n-np)/2)*laguerre(np,n-np,two_r_squared))
    end
    function Wigner_Fock_Recorrida(n,np,m,mp,G,pt)
        q=pt[1]
        p=pt[2]

        qrec=q+G*(m+mp)/2 #por el desplazamiento de alpha_m
        return exp(-im*p*G*(mp-m))*Wigner_Fock(n,np,qrec,p)
    end
    """
    ```julia
    function Wigner(system::QuantumDickeSystem,
        state::AbstractArray{<:Number},
        points::Vector{<:AbstractVector{<:Real}})
    ```
    Evaluates the Wigner function of `state` (which is a vector in the BCE) in all `points = [[Q1,q1,P1,p1], [Q2,q2,P2,p2], ...]`,
    returning an array with the results.
    
    Note: This function has not been thoroughly tested. It is based on Ref. [Pilatowsky2019Notes](@cite).
    """
    function Wigner(system::QuantumDickeSystem,
        state::AbstractArray{<:Complex},
        points::Vector{<:AbstractVector{<:Real}})
        
        ω₀,ω,γ,j,Nmax=get_params(system)
        puntos=points
        
        N=Int(2*j)
        G=2*γ/(ω*sqrt(N))

        puntos=classical_QqPp_to_quantum_rot_qpθφ.(j,puntos)
        legendre_sphPlm_cache=generate_legendre_sphPlm_cache(Int(2*j),puntos)
        #θ,φ=theta,phi


        function B(n,l,np,lp)
            altsm(n-l)*altsm(np-lp)*exp((logfact(n)+logfact(np))/2 -
                (logfact(l)+logfact(lp)))
        end





        W=fill!(similar(puntos,Complex{Float64}),0)
        _Wd=fill!(similar(puntos,Complex{Float64}),0)
        for m in -j:j, mp in -j:j

            mWd=Wigner_Dicke!(m,mp,j,puntos,_Wd,legendre_sphPlm_cache=legendre_sphPlm_cache)
            for n in 0:(Nmax-1),np in 0:n
                _factor=state[ind(n,m,N,j)]*conj(state[ind(np,mp,N,j)])
                if norm(_factor)<1e-6
                    continue
                end
                mW=_factor*mWd
                mW.*=Wigner_Fock_Recorrida.(n,np,m,mp,G,puntos)
                W+= real.(mW*(n==np ? 1 : 2))
             #   if n!=np
             #       W+=conj.(mW)
              #  end
            end


        end
        return real.(W)
    end
    function Wigner_Dicke!(m::Int64,mp::Int64,j::Int64,puntos::AbstractVector{<:AbstractVector{<:Real}},result::Array{<:Complex,1};legendre_sphPlm_cache=:nothing)

        result.=0.0+0.0im

        for k in abs(m-mp):Int(2 * j)
            i=m-mp
            _f=(-1)^(j - m - i) * clebschgordan(j,m, j,-mp, k,i)
            if _f !=0
                 result +=   _f*sph_harm_pt.(k, i, puntos,legendre_sphPlm_cache=legendre_sphPlm_cache)
            end
        end
        return result*sqrt((2*j+1)/(4*pi)) #normalizamos
    end
    function realmult(w2::Complex{Float64},factor)
        return real(factor*w2)
    end
    """
    ```julia
    function WignerProjqp(system::QuantumDickeSystem,
        states::AbstractArray{<:AbstractArray{<:Number}},
        ptsQP::Vector{<:AbstractVector{<:Real}};
        show_progress::Bool=true)
    ```
    Evaluates Wigner function, integrated over the bosonic variables ``q,p``, of each of the states in `states` 
    (which is a vector of complex vectors) 
    in all `ptsQP = [[Q1,P1], [Q2,P2], ...]`, returning an array with the results.
    See [this example](@ref wignerfuncexample). If `show_progress = true` it shows a progress bar.
        
    Note: This function has not been thoroughly tested. It is based on Ref. [Pilatowsky2019Notes](@cite).
    """
    function WignerProjqp(system::QuantumDickeSystem,
        states::AbstractArray{<:AbstractArray{<:Complex}},
        ptsQP::Vector{<:AbstractVector{<:Real}};
        show_progress::Bool=true)
        ω₀,ω,γ,j,Nmax=get_params(system)
        N=Int(2*j)
        G=2*γ/(ω*sqrt(N))
        puntos=[classical_QqPp_to_quantum_rot_qpθφ(j,[Q,0.0,P,0.0]) for (Q,P) in ptsQP]
        tmD=0.0
        tmF=0.0
        tmM=0.0
        _Wd=similar(puntos,Complex{Float64})
        AllWs=[similar(puntos,Float64) for s in states]
        fill!.(AllWs,0.0)

        legendre_sphPlm_cache=generate_legendre_sphPlm_cache(Int(2*j),puntos)
        if show_progress
            p = Progress(Int((2*j+1)^2*(Nmax-1)*Nmax/2),1)
        end
        for m in -j:j, mp in -j:j

            _Wd=Wigner_Dicke!(m,mp,j,puntos,_Wd,legendre_sphPlm_cache=legendre_sphPlm_cache)

            for n in 0:(Nmax-1),np in 0:n
                if np==n && mp==m
                    lag_fac=1.0
                else
                lag_fac=exp((logfact(np)-logfact(n))/2 - abs(G*(m-mp))^2/2 + (n-np)*log(G*abs(m-mp)))*sign(m-mp)^(n-np)*laguerre(np,n-np,abs(G*(m-mp))^2)

                end
                for i in 1:length(states)
                    _factor=states[i][ind(n,m,N,j)]*conj(states[i][ind(np,mp,N,j)])
                    finalfactor=_factor*lag_fac*(n==np ? 1 : 2)
                    if abs(finalfactor) ==0
                        continue
                    end
                 #   AllWs[i] = AllWs[i]+sumrealpart.(AllWs[i],_Wd,finalfactor) #forma eficiente de hacer W+=real.(finalfactor*_Wd)
                    AllWs[i].+=realmult.(_Wd,finalfactor)
                end
                if show_progress
                    next!(p)
                end

            end

        end

        if show_progress
            finish!(p)
        end
        return AllWs
    end

    function WignerProjp(system,states,puntos)
        ω₀,ω,γ,j,Nmax=get_params(system)
        N=Int(2*j)
        G=2*γ/(ω*sqrt(N))
        puntos=classical_QqPp_to_quantum_rot_qpθφ.(j,puntos)
        tmD=0.0
        tmF=0.0
        tmM=0.0
        _Wd=similar(puntos,Complex{Float64})
        _FactorHermite=similar(puntos,Float64)
        AllWs=[similar(puntos,Float64) for s in states]
        fill!.(AllWs,0.0)

        legendre_sphPlm_cache=generate_legendre_sphPlm_cache(Int(2*j),puntos)
        hermite_cache=generate_hermite_cache(Nmax-1,Iterators.flatten((x[1]+G*m for m in -j:j) for x in puntos))
        function Fock_n_np_wavefunction(n::Int64,m::Int64,np::Int64,mp::Int64,x::Array{Float64,1})
            return exp(-(logfact(np)+logfact(n)+(np+n)*log(2))/2 - (x[1]+G*m)^2/2 - (x[1]+G*mp)^2/2)*hermite_cache[n,x[1]+G*m]*hermite_cache[np,x[1]+G*mp]/sqrt(pi)
        end
        p = Progress(Int((2*j+1)^2*(Nmax-1)*Nmax/2),1)
        for m in -j:j, mp in -j:j

            _Wd=Wigner_Dicke!(m,mp,j,puntos,_Wd,legendre_sphPlm_cache=legendre_sphPlm_cache)

            for n in 0:(Nmax-1),np in 0:n

                _FactorHermite.=Fock_n_np_wavefunction.(n,m,np,mp,puntos)


                for i in 1:length(states)
                    _factor=states[i][ind(n,m,N,j)]*conj(states[i][ind(np,mp,N,j)])
                    finalfactor=_factor*(n==np ? 1 : 2)
                    if abs(finalfactor) ==0
                        continue
                    end
                    AllWs[i].+=realmult.(_Wd,finalfactor).*_FactorHermite
                end
                next!(p)

            end

        end

        finish!(p)
        return AllWs
    end
    """
    ```julia
    function WignerProjQP(system::QuantumDickeSystem,
        state::AbstractArray{<:Number},
        ptsqp::Vector{<:AbstractVector{<:Real}})
    ```
    Evaluates Wigner function, integrated over the atomic variables ``Q,P``, of `state` 
    in all `ptsqp = [[q1,p1], [q2,p2], ...]`, returning an array with the results.
        
    Note: This function has not been thoroughly tested. It is based on Ref. [Pilatowsky2019Notes](@cite).
    """
    function WignerProjQP(system::QuantumDickeSystem,
        state::AbstractArray{<:Number},
        ptsqp::Vector{<:AbstractVector{<:Real}})
        ω₀,ω,γ,j,Nmax=get_params(system)

        N=Int(2*j)
        G=2*γ/(ω*sqrt(N))
        puntos=[classical_QqPp_to_quantum_rot_qpθφ(j,[0.0,q,0.0,p]) for (q,p) in ptsqp]
        #θ,φ=theta,phi

        W=fill!(similar(puntos,Complex),0)
        _Wd=fill!(similar(puntos,Complex),0)
        for m in -j:j
            for n in 0:(Nmax-1),np in 0:n
                _factor=state[ind(n,m,N,j)]*conj(state[ind(np,m,N,j)])
                if norm(_factor)<1e-20
                    continue
                end
                mW=_factor*Wigner_Fock_Recorrida.(n,np,m,m,G,puntos)

                W+= real.(mW*(n==np ? 1 : 2))

            end


        end
        return real.(W)
    end
    """
    ```julia
    function random_state(system::QuantumDickeSystem,
                        n::Int=1;
                        kargs...)
    ```
    Generates `n` random states, returning them as columns in a matrix
    (or a vector if `n = 1`, which is the default). `kargs` are redirected to [`random_state!`](@ref): some
    additional keyword arguments are required: see [`random_state!`](@ref).
    """
    function random_state(system::QuantumDickeSystem,
        n::Int=1;kargs...)
        ω₀,ω,γ,j,Nmax=get_params(system)
        D=dimension(system)
        if n==1
            data = zeros(Complex{Float64},D)
        else
            data=zeros(Complex{Float64},D,n)
        end
        return random_state!(system,data;kargs...)
    end
    """
    ```julia
    function random_cₖ_generator(system::QuantumDickeSystem;
                                kargs...)
    ```
    Generates coefficients to build random states, which are given by
    ```math
        c_k(\\epsilon_k) \\sim e^{i \\theta_k} \\sqrt{\\frac{r_k\\rho(\\epsilon_k)}{\\nu(\\epsilon_k)}}.
    ```
    (See Refs. [Lerma2019](@cite), [Villasenor2020](@cite), and [Pilatowsky2021NatCommun](@cite)).
    It returns a function `cₖ(ϵ;number=1,cache=nothing)`, which produces n=`number` coefficients, storing them
    in `cache` (an n-vector), if given.
    
    The parameters of the equation above are determined as follows:

    # Required keyword arguments 
    - To determine  ``\\rho(\\epsilon)`` above, pass either: 
        - `σ` and `ϵ`, both real, in which case ``\\rho`` will be taken as a normal distribution
          centered at `ϵ` with standard deviation `σ`,
        - OR `envelope` as any of [`Distributions.UnivariateDistribution`](https://juliastats.org/Distributions.jl/v0.25/univariate/#Continuous-Distributions).
    - To determine ``r_k``, pass either:
        - `ensemble = :GUE`, which takes ``r_k`` from an exponential distribution `Exponential(0.91)` or `ensemble = :GOE`, which takes ``r_k`` from a ``\\chi^2``
          distribution with one degree of freedom,
        - OR `rₖ_distribution` as any of [`Distributions.UnivariateDistribution`](https://juliastats.org/Distributions.jl/v0.25/univariate/#Discrete-Distributions).
    - To determine ``\\theta_k``, pass either:
        - `phases = :real`, which takes ``\\theta_k`` to be random from ``\\{ -\\pi, \\pi \\}`` or `phases = :complex`, which takes ``\\theta_k`` to be uniform in ``[0, 2\\pi)``,
        - OR `ensemble = :GOE`, which automatically sets `phases = :real`, or `ensemble = :GUE`, which automatically sets `phases = :complex`.
    # Optional keyword arguments 
    - `divide_by_DoS` is a boolean. If `true` (default), leaves the density of states, ``\\nu(\\epsilon_k)`` in the equation above. A value of `false` removes it.
    """
    function random_cₖ_generator(system::QuantumDickeSystem;
                                σ::Real,
                                ϵ::Real,
                                envelope::UnivariateDistribution=Distributions.Normal(ϵ,σ),
                                ensemble::Symbol=:GUE,
                                rₖ_distribution::UnivariateDistribution=
                                if ensemble==:GUE
                                        Exponential(0.91)
                                elseif ensemble==:GOE
                                        Chisq(1)
                                else
                                    error("ensemble must be :GOE or :GUE, or pass rₖ_distribution and phases directly")
                                end,
                                phases::Symbol=if ensemble==:GUE
                                        :complex
                                elseif ensemble==:GOE
                                        :real
                                else
                                    error("ensemble must be :GOE or :GUE, or pass rₖ_distribution and phases directly")
                                end,
                                divide_by_DoS::Bool=true)

        ω₀,ω,γ,j,Nmax=get_params(system)
        systemC = system.classical_system
        ρ(ϵ)=Distributions.pdf(envelope,ϵ)
        ν(ϵ)=ClassicalDicke.density_of_states(systemC, ϵ; j)
        rₖ()=Distributions.rand(rₖ_distribution)
        rphase()=if phases==true || phases ==:complex
            exp(im*Random.rand()*2*pi)
        elseif phases==:real
            rand((1,-1))
        else
            1
        end
        function cₖ(ϵ;number=1,cache=nothing)
            if divide_by_DoS
                dos=ν(ϵ)
            else
                dos=1
            end
            if number>1 && cache==nothing
                cache=zeros(ComplexF64,number)
            end
            if dos==0
                if number==1
                    return 0
                else
                    cache*=0
                    return cache
                end
            else
                f()=sqrt(rₖ()*ρ(ϵ)/dos)*rphase()
                if number==1
                    return f()
                else
                    cache.=f.()
                    return cache
                end

            end
        end
        return cₖ
    end
    """
    ```julia
    function random_state!(system::QuantumDickeSystem,
        data::AbstractVecOrMat{<:Number};
        kargs...)
    ```
    Generates as many random states as columns in `data` saving the result to `data`.
    The general form of the random states in the eigenbasis is
    ``\\left | R \\right \\rangle = \\sum_k c_k \\left | E_k \\right \\rangle``, where
    the ``c_k`` are computed using [`random_cₖ_generator`](@ref). You may choose to select only
    positive or negative parity eigenstates (see below).
    # Arguments 
    - `system` should be an instance of [`DickeBCE.QuantumDickeSystem`](@ref).
    - `data` should be a vector or matrix of numbers, where the result is stored. It
      must have as many columns as random numbers you want and [`dimension`](@ref)`(system)` rows.
    # Required keyword arguments 
    - Please pass all the required keyword arguments of [`random_cₖ_generator`](@ref).
    - `eigenstates` and `eigenenergies` shoud be passed, being the ones returned by [`diagonalization`](@ref).
    # Optional keyword arguments 
    - All the optional keyword arguments of [`random_cₖ_generator`](@ref) may be passed.
    - `tol` may be a real number which determines the tolerance for the state convergence.
      A value `tol=0` builds the state with all of the eigenstates, and a positive number
      cuts the tails of the energy envelope up to that number. Higher values are faster, 
      but decrease precision (Defaults to `1e-6`)
    - `parity` may be `+`, `-`, or `nothing` (default). If it is not `nothing`, the random
      will be composed of only the eigenstates with this parity.
    - `parities` should be passed if `parity` is not `nothing`. It should be a vector of `-1`s 
      and `1`s containing the parities of all of the eigenstates. If it is not passed, it will be 
      computed with [`eigenstate_parities`](@ref), but this is slow. If you plan to use this
      function repeatedly, you should precompute `parities` by calling [`eigenstate_parities`](@ref) 
      yourself and then pass it to `parities` to avoid repeated calls to [`eigenstate_parities`](@ref).
    """
    function random_state!(system::QuantumDickeSystem,
        data::AbstractVecOrMat{<:Number};
        eigenstates::AbstractMatrix{<:Number},
        eigenenergies::AbstractVector{<:Real},
        σ::Real,
        ϵ::Real,
        envelope::UnivariateDistribution=Distributions.Normal(ϵ,σ),
        tol::Real=1e-6,
        parity::Union{Nothing,typeof(-),typeof(+)}=nothing,
        parities::AbstractArray{<:Integer,1}=if parity==nothing
            Int8[]
        else #@warn "If you pass parity=+-, it is more efficient for you to compute parities=eigenstate_parities(system,eigenstates), and pass it to random_state!(... parities=parities, ...). This prevents repeated calls to eigenstate_parities wich can be slow."
            eigenstate_parities(system,eigenstates) end,kargs...)
            
        if length(size(data))==2
           for i in 1:size(data)[2]
                random_state!(system,view(data,:,i);
                                eigenenergies=eigenenergies,eigenstates=eigenstates,σ=σ,ϵ=ϵ,envelope=envelope,tol=tol,parity=parity,parities=parities,kargs...)

           end
           return data
        end
        ω₀,ω,γ,j,Nmax=get_params(system)
        cₖ=random_cₖ_generator(system;ϵ=ϵ,σ=σ,envelope=envelope,kargs...)
        data.=0.0
        nrm=0.0
        ϵ0,ϵf=Distributions.quantile(envelope,tol/2),Distributions.quantile(envelope,1-tol/2)
        if parity!=nothing
            par=parity(1)
            is_correct_parity = k -> parities[k]==par #debes hacerlo con anonymus functions si no julia tiene un bug bien raro https://github.com/JuliaLang/julia/issues/15602

        else
            is_correct_parity = k -> true #idem
        end
        for (k,E) in enumerate(eigenenergies)

            ϵ=E/j
            if ϵ0<=ϵ<=ϵf && is_correct_parity(k)
                ck=cₖ(ϵ)
                nrm+= abs2(ck)
                data.+=ck.*view(eigenstates,:,k)
            end

        end
        data./=sqrt(nrm)
        return data
    end


    """
    ```julia
    function random_coherent_states_in_energy_shell(
        system::QuantumDickeSystem;
        ϵ::Real,
        N::Integer,
        dt::Real=3,
        tol::Real=1e-6,
        cachedata::Union{AbstractVector{<:Complex},Nothing}=nothing)
    ```
    Samples `N` points from the energy shell at `ϵ` using [`ClassicalDicke.classical_path_random_sampler`](@ref ClassicalDicke.classical_path_random_sampler),
    and then constructs the N-cat state of all the coherent states centered at those points. 
    """
    function random_coherent_states_in_energy_shell(
        system::QuantumDickeSystem;
        ϵ::Real,
        N::Integer,
        dt::Real=3,
        tol::Real=1e-6,
        cachedata::Union{AbstractVector{<:Complex},Nothing}=nothing)
        
        systemC = system.classical_system
        ω₀,ω,γ,j,Nmax=get_params(system)

        if cachedata===nothing
            cachedata=zeros(Complex,(Int(2j)+1)*(Nmax))
        end
        cachedata.=0.0im
        s=ClassicalDicke.classical_path_random_sampler(systemC,ϵ=ϵ,dt=dt)
        for i in 1:N
            θ=rand()*2*π
            p=s()
            coherent_state!(system,p,cachedata,extra_phase=exp(im*θ),tol=tol,add_state_to_data=true)
        end
        cachedata./=norm(cachedata)
        return cachedata
     end
end
     