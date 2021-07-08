module PhaseSpaces

export Q_of_θφ,Q_of_θϕ,P_of_θφ,P_of_θϕ,θ_of_QP,φ_of_QP,ϕ_of_QP,jz,jx,jy,arc_between_θφ,arc_between_θϕ,arc_between_QP
    """
    ```julia
    function Q_of_θϕ(θ,ϕ)
    ```
    
    Returns
    ```math
        Q=\\sqrt{2 (1 - \\cos(\\theta))}\\cos(\\phi).
    ```
    
    Also [`Q_of_θφ(θ,φ)`](@ref Q_of_θφ).
    """
    Q_of_θϕ(θ::Real,ϕ::Real)=sqrt(2*(1-cos(θ)))*cos(ϕ)
    
    """
    See [`Q_of_θϕ`](@ref)
    """
    Q_of_θφ(θ::Real,φ::Real)=Q_of_θϕ(θ,φ)

    """
    ```julia
    function P_of_θϕ(θ,ϕ)
    ```
    
    Returns
    ```math
        P=- \\sqrt{2 (1 - \\cos(\\theta))}\\sin(\\phi).
    ```
    Also [`P_of_θφ(θ,φ)`](@ref P_of_θφ).
    """
    P_of_θϕ(θ::Real,ϕ::Real)=-sqrt(2*(1-cos(θ)))*sin(ϕ)
    
    """
    See [`P_of_θϕ`](@ref)
    """
    P_of_θφ(θ::Real,φ::Real)=P_of_θϕ(θ,φ)
    """
    ```julia
    function θ_of_QP(Q,P)
    ```
    Returns
    ```math
        \\theta =- \\arccos\\left(1 - \\frac{Q^2 + P^2}{2}\\right).
    ```
    """
    θ_of_QP(Q::Real,P::Real)=acos(1-(Q^2+P^2)/2)
    """
    ```julia
    function ϕ_of_QP(Q,P)
    ```
    
    Returns
    ```math
        \\phi =- \\text{arctan2}(-P,Q)\\in [0,2 \\pi ],
    ```
    where ``\\text{arctan2}`` is the [2-argument arctangent](https://en.wikipedia.org/wiki/Atan2).
    
    Also [`φ_of_QP(Q,P)`](@ref φ_of_QP).
    """
    ϕ_of_QP(Q::Real,P::Real)=mod(atan(-P,Q),2*pi)
    
    """
    See [`ϕ_of_QP`](@ref)
    """
    φ_of_QP(Q::Real,P::Real)=ϕ_of_QP(Q,P)

    """
    ```julia
    function jz(Q,P)
    ```
    Returns
    ```math
        j_z = -  \\cos(\\theta(Q,P)).
    ```
    where ``\\text{arctan2}`` is the [2-argument arctangent](https://en.wikipedia.org/wiki/Atan2).
    """
    jz(Q::Real,P::Real)=-cos(θ_of_QP(Q,P))
    """
    ```julia
    function jx(Q,P)
    ```
    Returns
    ```math
        j_x = \\sin(\\theta(Q,P))\\cos(\\phi(Q,P)).
    ```
    """
    jx(Q::Real,P::Real)=sin(θ_of_QP(Q,P))*cos(φ_of_QP(Q,P))
    """
    ```julia
    function jy(Q,P)
    ```
    Returns 
    ```math
        j_y = \\sin(\\theta(Q,P))\\sin(\\phi(Q,P)).
    ```
    """
    jy(Q::Real,P::Real)=sin(θ_of_QP(Q,P))*sin(φ_of_QP(Q,P))
    """
    ```julia
    function arc_between_θϕ(θ1,φ1,θ2,φ2)
    ```
    Returns 
    ```math
        \\Theta = \\arccos(\\cos(\\theta_1)\\cos(\\theta_2)+ \\sin(\\theta_1)\\sin(\\theta_2)\\cos(\\phi_1 - \\phi_2)).
    ```
    Also [`arc_between_θφ(θ1,φ1,θ2,φ2)`](@ref arc_between_θφ).

    """
    arc_between_θϕ(θ1::Real,ϕ1::Real,θ2::Real,ϕ2::Real)=acos(clamp(cos(θ1)*cos(θ2)+ sin(θ1)*sin(θ2)*cos(ϕ1 - ϕ2),-1,1))
    
    """
    See [`arc_between_θϕ`](@ref).
    """
    arc_between_θφ(θ1::Real,φ1::Real,θ2::Real,φ2::Real)=arc_between_θϕ(θ1,φ1,θ2,φ2)
    """
    ```julia
    function arc_between_QP(Q1,P1,Q2,P2)
    ```
    Returns 
    ```julia
        arc_between_θϕ(θ_of_QP(Q1,P1),ϕ_of_QP(Q1,P1),θ_of_QP(Q2,P2),ϕ_of_QP(Q2,P2)).
    ```
    """
    arc_between_QP(Q1::Real,P1::Real,Q2::Real,P2::Real)=arc_between_θϕ(θ_of_QP(Q1,P1),ϕ_of_QP(Q1,P1),θ_of_QP(Q2,P2),ϕ_of_QP(Q2,P2))
end
