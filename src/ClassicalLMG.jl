module ClassicalLMG

export hamiltonian,ClassicalLMGSystem,P_of_ϵ
    import ..ClassicalSystems
    using ..PhaseSpaces
    using ParameterizedFunctions

    LMG_step = @ode_def begin
        dQ = integate_backwards*P*(Ω - ξ*Q^2/2) # =dH/dP
        dP = integate_backwards*(Q*(ξ*P^2/2 - (2*ξ + Ω)) + ξ*Q^3) # =-dH/dQ
    end Ω ξ integate_backwards
    function out_of_bounds_LMG(u,_,t)
        return 4-u[2]^2-u[1]^2<=0
    end
    varnames=[:Q,:P]
    """
    ```julia
    struct ClassicalLMGSystem <: ClassicalSystems.ClassicalSystem
    ```
    Subtype of [`ClassicalSystems.ClassicalSystem`](@ref) which represents the
    classical LMG model with the given parameters ``Ω`` and ``ξ``. See Eq. (2) of Ref. [Pilatowsky2020](@cite).
    To generate this struct, use the constructor
    ```julia
        ClassicalLMGSystem(;Ω::Real,ξ::Real)
    ```
    For example, `system = ClassicalLMGSystem(Ω=1, ξ=1)`.
    
    This struct may be passed to all functions in this module that require an instance of `ClassicalLMG.ClassicalLMGSystem`,
    as well as functions in other modules that require the abstract [`ClassicalSystems.ClassicalSystem`](@ref), 
    such as  [`ClassicalSystems.integrate`](@ref ClassicalSystems.integrate(::ClassicalSystems.ClassicalSystem,::AbstractVector{<:Real},::Real)).
    """
    struct ClassicalLMGSystem <: ClassicalSystems.ClassicalSystem
        parameters::Vector{Float64}
        ClassicalLMGSystem(;Ω::Real,ξ::Real)= new([Float64(Ω),Float64(ξ)])

    end
    function ClassicalSystems.parameters(system::ClassicalLMGSystem)
        return system.parameters
    end
    function ClassicalSystems.step(system::ClassicalLMGSystem)
        return LMG_step
    end
    function ClassicalSystems.out_of_bounds(system::ClassicalLMGSystem)
        return out_of_bounds_LMG
    end
    function ClassicalSystems.varnames(system::ClassicalLMGSystem)
        return varnames
    end
    function Base.show(io::IO, system::ClassicalLMGSystem)
        Ω,ξ=ClassicalSystems.parameters(system)
        print(io,"ClassicalLMGSystem(Ω = $Ω, ξ = $ξ)")
    end
    """
    ```julia
    function hamiltonian(system::ClassicalLMGSystem)
    ```
    Returns a classical Hamiltonian function `h(x)` where `x=[Q,P]`, which is given by
    Eq. (2) of Ref. [Pilatowsky2020](@cite).

    # Arguments
    - `system` should be generated with [`ClassicalLMG.ClassicalLMGSystem`](@ref).
    """
    function hamiltonian(system::ClassicalLMGSystem)
        Ω,ξ=ClassicalSystems.parameters(system)
        function H(u)
            Q,P=u
            Ω*(Q^2 + P^2)/2 + ξ*(Q^2 - P^2*Q^2/4 - Q^4/4) - Ω
        end
        return H
    end
    """
    ```julia
    function discriminant_of_P_solution(system::ClassicalLMGSystem,Q::Real,ϵ::Real)
    ```
    Returns the discriminant of the second degree equation in ``P`` given by
    ```math
        h_\\text{cl}(Q,P)=\\epsilon,
    ```
    where ``h_\\text{cl}`` is given by Eq. (2) of Ref. [Pilatowsky2020](@cite).

    # Arguments
    - `system` should be generated with [`ClassicalLMGSystem`](@ref).
    - `Q` and `ϵ` are the values of ``Q`` and ``\\epsilon``, respectively.
    """
    function discriminant_of_P_solution(system::ClassicalLMGSystem,Q::Real,ϵ::Real)
        Ω,ξ=ClassicalSystems.parameters(system)
        return (-4*ϵ + 4*Q^2*ξ - Q^4*ξ - 4*Ω + 2*Q^2*Ω)/(Q^2*ξ - 2*Ω)
    end
    """
    ```julia
    function P_of_ϵ(system::ClassicalLMGSystem;
        Q::Real,
        ϵ::Real,
        sgn::Union{typeof(-),typeof(+)}=+,
        returnNaNonError::Bool=true)
    ```
    Returns the solutions ``P_\\pm`` of the second degree equation in ``P`` given by
    ```math
        h_\\text{cl}(Q,P)=\\epsilon,
    ```
    where ``h_\\text{cl}`` is given by Eq. (2) of Ref. [Pilatowsky2020](@cite).

    # Arguments
    - `system` should be generated with [`ClassicalLMGSystem`](@ref).
    # Keyword arguments
    - `Q` and `ϵ` are values of ``Q`` and ``\\epsilon``, respectively.
    - `sgn` is `+` for ``P_+`` and `-` for ``P_-``
    - If `returnNaNonError` is `true`, then `NaN` is returned if there are no solutions. If it is `false`, and error is raised.
    """
    function P_of_ϵ(system::ClassicalLMGSystem;
        Q::Real,
        ϵ::Real,
        sgn::Union{typeof(-),typeof(+)}=+,
        returnNaNonError::Bool=true)
        Δ=discriminant_of_P_solution(system,Q,ϵ)
        if Δ<0
            if returnNaNonError
                return NaN
            else
                ϵmin= minimum_ϵ_for(system;Q=Q)
                error("The minimum energy for the given parameters is $ϵmin. There is no P such that (Q=$Q,P) has ϵ=$ϵ.")
            end
        end
        return sgn(0.0,sqrt(Δ))
    end
    """
    ```julia
    function minimum_ϵ_for(system::ClassicalLMGSystem;
        Q::Union{Real,Nothing}=nothing,
        P::Union{Real,Nothing}=nothing)
    ```
    Returns the minimum energy ``\\epsilon`` when constraining the system to 
    one fixed value of the coordinates ``Q`` or ``P``.

    # Arguments
    - `system` should be generated with [`ClassicalLMGSystem`](@ref).
    # Keyword arguments
    - You may pass either ``P`` or ``Q``.
    """
    function minimum_ϵ_for(system::ClassicalLMGSystem;
        Q::Union{Real,Nothing}=nothing,
        P::Union{Real,Nothing}=nothing)
        if count(x->x==nothing,[Q,P])!=1
            error("You have to pass one of Q or P")
        end
        Ω,ξ=ClassicalSystems.parameters(system)
        if P==nothing
            return (2*Ω*Q^2 + 4*Q^2*ξ - Q^4*ξ)/4 - Ω
        else
            return Ω*P^2/2 - Ω
        end
    end



    """
    ```julia
    function Point(;Q::Real,P::Real)
    ```
    Returns the list `[Q,P]`
    """
    Point(;Q::Real,P::Real)=Float64[Q,P]
    """
    ```julia
    function Pointθϕ(;θ::Real,ϕ::Real)
    ```
    Returns a list `[Q,P]`, where `Q` and `P` are calculated from `θ` and `ϕ` using [`PhaseSpaces.Q_of_θϕ`](@ref) and  [`PhaseSpaces.P_of_θϕ`](@ref).
    """
    Pointθϕ(;θ::Real,ϕ::Real)=Point(Q=Q_of_θφ(θ,ϕ),P=P_of_θφ(θ,ϕ))
    Pointθφ(;θ::Real,φ::Real)=Pointθϕ(θ=θ,ϕ=φ)
    """
    ```julia
    function Point(system::ClassicalLMGSystem;
        Q::Real,
        ϵ::Real,
        sgn::Union{typeof(-),typeof(+)} = +)
    ```
    Returns a list `[Q,P]`, where `P` is calculated with [`P_of_ϵ`](@ref). 
    If there are no solutions for ``P``, an error is raised.
    """
    Point(system::ClassicalLMGSystem;
        Q::Real,
        ϵ::Real,
        sgn::Union{typeof(-),typeof(+)}=+)=Point(Q=Q,P=P_of_ϵ(system;Q=Q,ϵ=ϵ,sgn=sgn,returnNaNonError=false))
end
