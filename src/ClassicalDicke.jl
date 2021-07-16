module ClassicalDicke
export q_of_ϵ,q_sign,minimum_ϵ_for,Point,Pointθϕ,Pointθφ,discriminant_of_q_solution,
       energy_shell_volume,density_of_states,energy_width_of_coherent_state,
       maximum_P_for_ϵ,minimum_nonnegative_Q_for_ϵ,maximum_Q_for_ϵ,minimum_energy,minimum_energy_point,
       normal_frequency,phase_space_dist_squared,ClassicalDickeSystem,hamiltonian
    using ..ClassicalSystems
    using ..PhaseSpaces
    using ParameterizedFunctions
    using QuadGK
    
    dicke_step = @ode_def begin
        dQ = integate_backwards*(P*ω₀-γ*P*q*Q/sqrt(4-P^2-Q^2)) # =dH/dP
        dq = integate_backwards*p*ω # =dH/dp
        dP = integate_backwards*(γ*q*Q^2/sqrt(4-P^2-Q^2)-γ*q*sqrt(4-P^2-Q^2)-Q*ω₀) # =-dH/dQ
        dp = integate_backwards*(-γ*Q*sqrt(4-P^2-Q^2)-q*ω) #-dH/dq

    end ω₀ ω γ integate_backwards
    function dicke_out_of_bounds(u,_,t)
        return 4-u[3]^2-u[1]^2<=0
    end
    varnames=[:Q,:q,:P,:p]
    """
    ```julia
    struct ClassicalDickeSystem <: ClassicalSystems.ClassicalSystem
    ```
    Subtype of [`ClassicalSystems.ClassicalSystem`](@ref) which represents the
    classical Dicke model with the given parameters ``ω_0``, ``ω``, and ``γ``. See Eq. (5) of Ref. [Pilatowsky2021](@cite).
    To generate this struct, use the constructor
    ```julia
    ClassicalDickeSystem(;ω₀::Real,ω::Real,γ::Real)
    ```
    For example, `system = ClassicalDickeSystem(ω₀=1, ω=1, γ=1)`.
    
    This struct may be passed to all functions in this module that require an instance of `ClassicalDicke.ClassicalDickeSystem`,
    as well as functions in other modules that require the abstract [`ClassicalSystems.ClassicalSystem`](@ref), 
    such as  [`ClassicalSystems.integrate`](@ref).
    """
    struct ClassicalDickeSystem <: ClassicalSystems.ClassicalSystem
        parameters::Vector{Float64}
        ClassicalDickeSystem(;ω₀::Real,ω::Real,γ::Real)= new(Float64.([ω₀,ω,γ]))

    end
    function ClassicalSystems.parameters(system::ClassicalDickeSystem)
        return system.parameters
    end
    function ClassicalSystems.step(system::ClassicalDickeSystem)
        return dicke_step
    end
    function ClassicalSystems.out_of_bounds(system::ClassicalDickeSystem)
        return dicke_out_of_bounds
    end
    function ClassicalSystems.varnames(system::ClassicalDickeSystem)
        return varnames
    end
    function Base.show(io::IO, cds::ClassicalDickeSystem)
        ω₀,ω,γ=ClassicalSystems.parameters(cds)
        print(io,"ClassicalDickeSystem(ω₀ = $ω₀, ω = $ω, γ = $γ)")
    end
    """
    ```julia
    function hamiltonian(system::ClassicalDickeSystem)
    ```
    Returns a classical Hamiltonian function `h(x)` where `x=[Q,q,P,p]`, which is given by
    Eq. (5) of Ref. [Pilatowsky2021](@cite).

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    """
    function hamiltonian(system::ClassicalDickeSystem)
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        function H(u)
            Q,q,P,p=u
            r=1-(Q^2+P^2)/4
            (ω₀/2)*(Q^2+P^2)+(ω/2)*(q^2+p^2)+2*γ*q*Q*sqrt(r)-ω₀
        end
        return H
    end
    """
    ```julia
    function density_of_states(system::ClassicalDickeSystem;j::Real,ϵ::Real)
    ```
    Returns the semiclassical density of states (DoS) ``ν(ϵ)``, in units of ``1/ϵ``. This is computed
    with an expression similar to Eq. (19) of Ref. [Bastarrachea2014](@cite), where
    we add an additional factor of ``j`` to have units of ``1/ϵ`` instead of ``1/E``, and the integral is performed
    with a change of variable.

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `j` is the value of ``j``
    - `ϵ` is the scaled energy ``ϵ=E/j``

    See [Plotting the density of states](@ref) for an example.
    """
    function density_of_states(system::ClassicalDickeSystem;j::Real,ϵ::Real)
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        energy_shell_volume(system,ϵ)*(j/(2pi))^2
    end

    """
    ```julia
    function energy_shell_volume(system::ClassicalDickeSystem,ϵ::Real)
    ```
    Returns the volume of the classical energy shell in the phase space, that is, ``\\mathcal{V}(\\mathcal{M}_ϵ)`` in Eq. (27) of Ref. [Villasenor2021](@cite).

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `ϵ` is the scaled energy ``ϵ=E/j``
    """
    function energy_shell_volume(system::ClassicalDickeSystem,ϵ::Real)
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        if ϵ>=ω₀
            return 8*pi^2/ω
        end
        Pmax(Q)=sqrt(max((ω*(2*ϵ + 2*ω₀) - Q^2*((-4 + Q^2)*γ^2 + ω*ω₀))/(Q^2*γ^2 + ω*ω₀),0))
        return quadgk(Pmax,minimum_nonnegative_Q_for_ϵ(system,ϵ),maximum_Q_for_ϵ(system,ϵ))[1]*2*pi*4/ω
    end

    """
    ```julia
    function normal_frequency(system::ClassicalDickeSystem,
        signo::Union{typeof(-),typeof(+)}=+)
    ```
    Returns the ground-state normal frequency, that is, ``\\Omega_{\\epsilon_{\\text{GS}}}^{A,B}`` at the bottom of page 3 of Ref. [Pilatowsky2021](@cite).

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `signo` is `-` for ``\\Omega^A`` and `+` for ``\\Omega^B``.

    Note: This function currently only works for the supperadiant phase.
    """
    function normal_frequency(system::ClassicalDickeSystem,signo::Union{typeof(-),typeof(+)}=+)
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        γc=sqrt(ω₀*ω)/2
        if γc>γ
            error("No methods for the normal phase (please add it!)")
            return
        end
        return sqrt((1/(2*γc^4))*signo(ω^2*γc^4 + ω₀^2*γ^4, sqrt((ω₀^2*γ^4 -ω^2*γc^4)^2 + 4*γc^8*ω^2*ω₀^2)))
    end

    """
    ```julia
    function minimum_energy_point(system::ClassicalDickeSystem,
        Qsign::Union{typeof(-),typeof(+)}=+)
    ```
    Returns the ground-state coordinate, that is,  `(0, 0, 0, 0)` for the normal
    phase and ``\\mathbf{x}_\\text{GS}`` below Eq. (7) of Ref. [Pilatowsky2021](@cite) for the
    superradiant phase.

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `Qsign` toggles the sign of the ``Q`` coordinate in the superradiant phase, 
      that is, `+` for ``\\mathbf{x}_\\text{GS}`` and  `-` for ``\\widetilde{\\mathbf{x}}_\\text{GS}``.
    """
    function minimum_energy_point(system::ClassicalDickeSystem,signoQ::Union{typeof(-),typeof(+)}=+)
        local signoq
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        γc=sqrt(ω₀*ω)/2
        if γc >= γ
            return Point(Q=0,P=0,p=0,q=0)
        end
        if signoQ==+
            signoq=-
        else
            signoq=+
        end
        return Point(Q=signoQ(sqrt(2-ω*ω₀/(2*γ^2))),q=signoq(sqrt(4*γ^2/ω^2 - ω₀^2/(4*γ^2))),P=0,p=0)
    end

    """
    ```julia
    function minimum_energy(system::ClassicalDickeSystem)
    ```
    Returns the ground-state energy, as given by multiplying Eq. (15) of Ref. [Bastarrachea2014](@cite) by
    ``\\omega_0`` (because in that work they use ``\\epsilon =E/(\\omega_0 j)``, and in this module we
    use ``\\epsilon = E/j``).

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    """
    function minimum_energy(system::ClassicalDickeSystem)
        ω₀,ω,γ = ClassicalSystems.parameters(system)
        γc=sqrt(ω₀*ω)/2
        
        if γ ≤ γc
            return -ω₀
        else
            return - ω₀*( (γc/γ)^2 + (γ/γc)^2 )/2
        end
    end

    """
    ```julia
    function energy_width_of_coherent_state(system::ClassicalDickeSystem,
        x::AbstractVector{<:Real},
        j::Real)
    ```
    Returns the energy width ``\\sigma`` of the coherent state ``\\left | \\mathbf{x}\\right \\rangle``, in units of ``\\epsilon``. This
    quantity is given by ``\\sigma_D/j`` with ``\\sigma_D`` as in App. A of Ref. [Lerma2018](@cite).

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `x` is the coordinate ``\\mathbf{x}`` of the coherent state in the format `[Q,q,P,p]`.
    - `j` is the value of ``j``.
    """
    function energy_width_of_coherent_state(system::ClassicalDickeSystem,x::AbstractVector{<:Real},j::Real)
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        Q,q,P,p=x
        θ=PhaseSpaces.θ_of_QP(Q,P)
        ϕ=PhaseSpaces.ϕ_of_QP(Q,P)

        Ω₁=j*(ω^2*(q^2+p^2)/2+ω₀^2*sin(θ)^2/2 +2*γ^2*((sin( θ)^2*sin(ϕ)^2+cos(θ)^2)*q^2 + sin(θ)^2*cos(ϕ)^2) +2*γ*q*(ω*cos(ϕ) + ω₀*cos(θ)*cos(ϕ))*sin(θ))
        #note that the sign of the third term is flipped with respect to Lerma2018, because they use cos θ = jz and we use cos θ = - jz.

        Ω₂=γ^2*(sin(θ)^2*sin(ϕ)^2 + cos(θ)^2)
        return sqrt(Ω₁ + Ω₂)/j
    end

    """
    ```julia
    function discriminant_of_q_solution(system::ClassicalDickeSystem; 
        Q::Real,
        P::Real,
        p::Real,
        ϵ::Real)
    ```
    Returns the discriminant of the second degree equation in ``q`` given by
    ```math
        h_\\text{cl}(Q,q,P,p)=\\epsilon,
    ```
    where ``h_\\text{cl}`` is given by Eq. (5) of Ref. [Pilatowsky2021](@cite).

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `Q`, `P`, `p`, and `ϵ` are the values of ``Q``, ``P``, ``p``, and ``\\epsilon``, respectively.
    """
    function discriminant_of_q_solution(system::ClassicalDickeSystem; Q::Real,P::Real,p::Real,ϵ::Real)
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        r=1-(Q^2+P^2)/4
        return 4*γ^2*Q^2*r - p^2*ω^2 - ω*ω₀*(P^2 + Q^2 - 2) + 2*ω*ϵ
    end


    """
    ```julia
    function q_of_ϵ(system::ClassicalDickeSystem;
        Q::Real,
        P::Real,
        p::Real,
        ϵ::Real,
        signo::Union{typeof(-),typeof(+)}=+,
        returnNaNonError::Bool=true)
    ```
    Returns the solutions ``q_\\pm`` of the second degree equation in ``q`` given by
    ```math
        h_\\text{cl}(Q,q,P,p)=\\epsilon,
    ```
    where ``h_\\text{cl}`` is given by Eq. (5) of Ref. [Pilatowsky2021](@cite).

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `Q`, `P`, `p`, and `ϵ` are values of ``Q``, ``P``, ``p``, and ``\\epsilon``, respectively.
    - `signo` is `+` for ``q_+`` and `-` for ``q_-``
    - If `returnNaNonError` is `true`, then `NaN` is returned if there are no solutions. If it is `false`, and error is raised.
    """
    function q_of_ϵ(system::ClassicalDickeSystem;
        Q::Real,
        P::Real,
        p::Real,
        ϵ::Real,
        signo::Union{typeof(-),typeof(+)}=+,
        returnNaNonError::Bool=true)
        Δ=discriminant_of_q_solution(system;Q=Q,P=P,p=p,ϵ=ϵ)
        if Δ<0
            if returnNaNonError
                return NaN
            else
                ϵmin= minimum_ϵ_for(system;Q=Q,P=P,p=p)
                error("The minimum energy for the given parameters is $ϵmin. There is no q such that (Q=$Q,q,P=$P,p=$p) has ϵ=$ϵ.")
            end
        end
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        return signo(-2*γ*Q*sqrt(1-(Q^2+P^2)/4),sqrt(Δ))/ω
    end
    """
    ```julia
    function q_sign(system::ClassicalDickeSystem,
        x::AbstractVector{<:Real},
        ϵ::Real=hamiltonian(system)(x))
    ```
    Returns the sign of the root of the second degree equation in ``q`` given by
    ```math
        h_\\text{cl}(Q,q,P,p)=\\epsilon.
    ```
    That is, this function returns `+` if `q=x[2] ≈ q_of_ϵ(system;Q=x[1],P=x[3],p=x[4],ϵ=ϵ,signo=+)`
    and returns `-` if `q=x[2] ≈ q_of_ϵ(system;Q=x[1],P=x[3],p=x[4],ϵ=ϵ,signo=-)`.
    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `x` is a vector in the form `[Q,q,P,p]`.
    - `ϵ` should be the energy of `x`. If this is not passed, it is computed using [`hamiltonian`](@ref ClassicalDicke.hamiltonian)`(system)(x)`. 
    """
    function q_sign(system::ClassicalDickeSystem,
        x::AbstractVector{<:Real},
        ϵ::Real=hamiltonian(system)(x))
        Q,q,P,p = x
        q₊ = q_of_ϵ(system; Q=Q, P=P, p=p, ϵ=ϵ, signo=+)
        q₋ = q_of_ϵ(system; Q=Q, P=P, p=p, ϵ=ϵ, signo=-)
        if abs(q - q₊) < abs(q - q₋) 
          return +
        else
          return -
        end
    end
    """
    ```julia
    function minimum_ϵ_for(system::ClassicalDickeSystem;
        Q::Union{Real,Nothing}=nothing,
        q::Union{Real,Nothing}=nothing,
        P::Union{Real,Nothing}=nothing,
        p::Union{Real,Nothing}=nothing)
    ```
    Returns the minimum energy ``\\epsilon`` when constraining the system to three fixed values of the coordinates ``Q``, ``q``, ``P``, ``p``.

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - You may pass either ``(Q,q,P)`` or ``(q,P,p)``. The other combinanations are not implemented.

    This function can be especially useful to draw contours of the available phase space (see [Drawing contours of the available phase space](@ref))
    """
    function minimum_ϵ_for(system::ClassicalDickeSystem;
        Q::Union{Real,Nothing}=nothing,
        q::Union{Real,Nothing}=nothing,
        P::Union{Real,Nothing}=nothing,
        p::Union{Real,Nothing}=nothing)
        if count(x->x==nothing,[q,Q,p,P])!=1
            error("You have to give three of Q,q,P,p")
        end
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        if q==nothing
            r=1-(Q^2+P^2)/4
            return -(4*γ^2*Q^2*r - p^2*ω^2 - ω*ω₀*(P^2 + Q^2 - 2))/(2*ω)
        end
        if Q==nothing
            return (2*p^2*ω + 2*q^2*ω + P^2*ω₀- sqrt((-4 + P^2)^2* (4*q^2*γ^2 + ω₀^2)))/4
        end
        error("Sorry, that combinanation of variables is not implemented")
    end
    """
    ```julia
    function maximum_P_for_ϵ(system::ClassicalDickeSystem,ϵ::Real)
    ```
    Computes the maximum value of the parameter ``P`` accessible to the system at energy ``\\epsilon``.

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `ϵ` is the scaled energy ``ϵ=E/j``.
    """
    function maximum_P_for_ϵ(system::ClassicalDickeSystem,ϵ::Real)
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        if ϵ>=ω₀
            return 2.0
        elseif ω*ω₀>sqrt(2)*γ*sqrt(ω*(ω₀-ϵ))
            return sqrt(2*(ω₀+ϵ)/ω₀)
        else
            return sqrt(clamp((2*√(2*γ^6*ω*(ω₀-ϵ)))/(-γ^4) + ω*ω₀/(γ^2)+4,0,4))
        end
    end
    """
    ```julia
    function maximum_Q_for_ϵ(system::ClassicalDickeSystem,ϵ::Real)
    ```
    Computes the maximum value of the parameter ``Q`` accessible to the system at energy ``\\epsilon``.

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `ϵ` is the scaled energy ``ϵ=E/j``.
    """
    function maximum_Q_for_ϵ(system::ClassicalDickeSystem,ϵ::Real)
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        if ϵ>=ω₀
            return 2.0
        end
         return sqrt(clamp((4γ^2-ω*ω₀)/(2*γ^2) + (1/2)*sqrt(max(0,(16*γ^4+8γ^2*ϵ*ω+ω^2*ω₀^2)/γ^4)),0,4))
    end
    """
    ```julia
    function minimum_nonnegative_Q_for_ϵ(system::ClassicalDickeSystem,ϵ::Real)
    ```
    Computes the minimum nonnegative value of the parameter ``Q`` accessible to the system at energy ``\\epsilon``.

    # Arguments:
    - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
    - `ϵ` is the scaled energy ``ϵ=E/j``.
    """
    function minimum_nonnegative_Q_for_ϵ(system::ClassicalDickeSystem,ϵ::Real)
        ω₀,ω,γ=ClassicalSystems.parameters(system)
        if ϵ >=- ω₀
            return 0.0
        end
         return sqrt(clamp((4γ^2-ω*ω₀)/(2*γ^2) - (1/2)*sqrt(max(0,(16*γ^4+8γ^2*ϵ*ω+ω^2*ω₀^2)/γ^4)),0,4))
    end




    """
    ```julia
    function Point(;Q::Real,q::Real,P::Real,p::Real)
    ```
    Returns the array `[Q,q,P,p]`. This ordering is used throughout this package.
    """
    Point(;Q::Real,q::Real,P::Real,p::Real) = Float64[Q,q,P,p]

    """
    ```julia
    function Point(system::ClassicalDickeSystem;
        Q::Real,
        P::Real,
        p::Real,
        ϵ::Real,
        signo::Union{typeof(-),typeof(+)} = +)
    ```
    Returns a list `[Q,q,P,p]`, where `q` is calculated with [`ClassicalDicke.q_of_ϵ`](@ref). 
    If there are no solutions for ``q``, an error is raised.
    """
    Point(system::ClassicalDickeSystem;
        Q::Real,
        P::Real,
        p::Real,
        ϵ::Real,
        signo::Union{typeof(-),typeof(+)} = +) = Point(q=q_of_ϵ(system::ClassicalDickeSystem,Q=Q,P=P,p=p,ϵ=ϵ,signo=signo,returnNaNonError=false),
                                                        P=P, p=p, Q=Q)

    """
    ```julia
    function Pointθϕ(;θ::Real,ϕ::Real,q::Real,p::Real)
    ```
    Returns a list `[Q,q,P,p]`, where `Q` and `P` are calculated from `θ` and `ϕ` using [`PhaseSpaces.Q_of_θϕ`](@ref) and  [`PhaseSpaces.P_of_θϕ`](@ref).
    """
    Pointθϕ(;θ::Real,ϕ::Real,q::Real,p::Real)=Point(Q=Q_of_θϕ(θ,ϕ),q=q,P=P_of_θϕ(θ,ϕ),p=p)
    Pointθφ(;θ::Real,φ::Real,q::Real,p::Real)=Pointθϕ(θ=θ,ϕ=φ,q=q,p=p)

    """
    ```julia
    function phase_space_dist_squared(x::AbstractVector{<:Real},y::AbstractVector{<:Real})
    ```
    Returns the phase-space distance ``d_{\\mathcal{M}}(\\mathbf{x},\\mathbf{y})`` (See App. C of Ref [Pilatowsky2021](@cite)),
    where `x` and `y` are vectors in the form `[Q,q,P,p]`.
    """
    function phase_space_dist_squared(x::AbstractVector{<:Real},y::AbstractVector{<:Real})
        Q1,q1,P1,p1=x
        Q2,q2,P2,p2=y
        return (q1-q2)^2 + (p1-p2)^2 + arc_between_QP(Q1,P1,Q2,P2)^2
    end
    # """
    # ```julia
    # function WignerHWxSU2_fixed_ϵ(system::ClassicalDickeSystem;Q,P,p,ϵ,j,signoq=+)
    # ```
    # Returns an instance of [`TruncatedWignerApproximation.PhaseSpaceDistribution`](@ref) that samples points from the classical energy shell
    # using a [Random Walk Metropolis-Hastings algorithm](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm) implemented in [Mamba](https://mambajl.readthedocs.io/en/latest/samplers/rwm.html?highlight=RWMVariate).
    # """
    #  function WignerHWxSU2_fixed_ϵ(system::ClassicalDickeSystem;Q,P,p,ϵ,j,signoq=+)
    #         Δ(Q,P)=discriminant_of_q_solution(system;Q=Q,P=P,p=p,ϵ=ϵ)
    #         x₀=Point(system,Q=Q,P=P,p=p,ϵ=ϵ,signo=signoq);
    #         dst=TruncatedWignerApproximation.coherent_Wigner_HWxSU2(x₀,j);
    #         function probability_density(u)
    #             Q,q,P,p=u
    #             return dst.probability_density(Point(system,Q=Q,P=P,p=p,ϵ=ϵ,signo=signoq))/sqrt(Δ(Q,P))
    #         end
    #         function logf(u)
    #             θ,ϕ,p=u
    #             Q=PhaseSpaces.Q_of_θϕ(θ,ϕ)
    #             P=PhaseSpaces.P_of_θϕ(θ,ϕ)
    #             try
    #                 return log(probability_density([Q,NaN,P,p]))#*sin(θ)?
    #             catch
    #                 return -Inf
    #             end
    #           end
    #
    #         m = RWMVariate([PhaseSpaces.θ_of_QP(Q,P),PhaseSpaces.ϕ_of_QP(Q,P),p], [1/sqrt(j),1/sqrt(j),1/sqrt(j)],logf,proposal = Normal)
    #         function sample()
    #                 θ,ϕ,p=sample!(m)
    #                 return Point(system,Q=PhaseSpaces.Q_of_θϕ(θ,ϕ),P=PhaseSpaces.P_of_θϕ(θ,ϕ),p=p,ϵ=ϵ,signo=signoq)
    #         end
    #         for i in 1:10000
    #             sample!(m) #bake
    #         end
    #         TruncatedWignerApproximation.PhaseSpaceDistribution(probability_density,sample)
    #     end

        """
        ```julia
        function classical_path_random_sampler(system::ClassicalDickeSystem;ϵ::Real,dt::Real=3)
        ```
        This function returns a function `sample()`, which produces random points within the classical
        energy shell at energy `ϵ`. The function `sample()` returns points from a fixed chaotic 
        trajectory picked at random separated by a fixed time interval  `dt`. If the classical dynamics are ergodic, 
        **which only happens in chaotic regions,** this produces random points in the energy shell.
        # Arguments:
        - `system` should be generated with [`ClassicalDicke.ClassicalDickeSystem`](@ref).
        - `ϵ` is the scaled energy ``ϵ=E/j`` of the energy shell from where to sample.
        - `dt` is the fixed time interval that separates the points that are returned by `sample()`.
        """
        function classical_path_random_sampler(system::ClassicalDickeSystem;ϵ::Real,dt::Real=3)
             ω₀,ω,γ=ClassicalSystems.parameters(system)
             s=sqrt(16*γ^4 + 8*γ^2*ϵ*ω + ω^2*ω₀^2)/(2*γ^2)
             a=2- ω*ω₀/(2*γ^2)
             mx=sqrt(a+s)
             if a<s
                 mn=0
             else
                 mn=sqrt(a-s)
             end
             if mn==mx
                 Q=mn #estamos en la energia del estado base
             else
                 Q=Distributions.rand(Distributions.Uniform(mn,mx))
             end
             Q=Q*rand((-1,1))
             u0=ClassicalDicke.Point(system,Q=Q,P=0,p=0,ϵ=ϵ)
             function sample()
                 u0=ClassicalSystems.integrate(system; t=dt, u₀=u0, save_everystep=false,show_progress=false).u[end]
                 if rand((true,false))
                     u0[1]=-u0[1]
                     u0[2]=-u0[2]
                 end
                 return u0
             end
             for i in 1:10
                 sample()
             end
             return sample
         end

end
           