module PhaseSpaces

export Q_of_θϕ,P_of_θϕ,θ_of_QP,ϕ_of_QP,jz,jx,jy,arc_between_θϕ,arc_between_QP
    Q_of_θϕ(θ,ϕ)=sqrt(2*(1-cos(θ)))*cos(ϕ)
    P_of_θϕ(θ,ϕ)=-sqrt(2*(1-cos(θ)))*sin(ϕ)

    θ_of_QP(Q,P)=acos(1-(Q^2+P^2)/2)
    ϕ_of_QP(Q,P)=mod(atan(-P,Q),2*pi)

    jz(Q,P)=-cos(θ_of_QP(Q,P))
    jx(Q,P)=sin(θ_of_QP(Q,P))*cos(ϕ_of_QP(Q,P))
    jy(Q,P)=sin(θ_of_QP(Q,P))*sin(ϕ_of_QP(Q,P))

    arc_between_θϕ(θ1,ϕ1,θ2,ϕ2)=acos(clamp(cos(θ1)*cos(θ2)+ sin(θ1)*sin(θ2)*cos(ϕ1 - ϕ2),-1,1))
    arc_between_QP(Q1,P1,Q2,P2)=arc_between_θϕ(θ_of_QP(Q1,P1),ϕ_of_QP(Q1,P1),θ_of_QP(Q2,P2),ϕ_of_QP(Q2,P2))
end
