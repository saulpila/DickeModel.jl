module ClassicalDicke

export hamiltonian, q_of_ε,minimum_ε_for,Point,discriminant_of_q_solution,energy_shell_volume,density_of_states
    using ..ClassicalSystems
    using ..PhaseSpaces
    using ParameterizedFunctions
    using ..TruncatedWignerApproximation
    using Mamba
    using QuadGK
    function hamiltonian(sistema::ClassicalSystems.ClassicalSystem)
        ω₀,ω,γ=sistema.params
        function H(u)
            Q,q,P,p=u
            r=1-(Q^2+P^2)/4
            (ω₀/2)*(Q^2+P^2)+(ω/2)*(q^2+p^2)+2*γ*q*Q*sqrt(r)-ω₀
        end
        return H
    end
    function density_of_states(sistemaC::ClassicalSystems.ClassicalSystem;j,ε)
        ω₀,ω,γ=sistemaC.params
        energy_shell_volume(sistemaC,ε)*(j/(2pi))^2
    end
    function energy_shell_volume(sistemaC::ClassicalSystems.ClassicalSystem,ε)
        ω₀,ω,γ=sistemaC.params
        if ε>=ω₀
            return 8*pi^2/ω
        end
        Pmax(Q)=sqrt(max((ω*(2*ε + 2*ω₀) - Q^2*((-4 + Q^2)*γ^2 + ω*ω₀))/(Q^2*γ^2 + ω*ω₀),0))
        return quadgk(Pmax,0,maximum_Q_for_ε(sistemaC,ε))[1]*2*pi*4/ω
    end
    function normal_frequency(sistema::ClassicalSystems.ClassicalSystem,signo::Union{typeof(-),typeof(+)}=+)
        ω₀,ω,γ=sistema.params
        γc=sqrt(ω₀*ω)/2
        if γc>γ
            error("Aún no se ponen métodos para la fase normal")
            return
        end
        return sqrt((1/(2*γc^4))*signo(ω^2*γc^4 + ω₀^2*γ^4, sqrt((ω₀^2*γ^4 -ω^2*γc^4)^2 + 4*γc^8*ω^2*ω₀^2)))
    end
    function minimum_energy_point(sistema::ClassicalSystems.ClassicalSystem,signoQ::Union{typeof(-),typeof(+)}=+)
        local signoq
        ω₀,ω,γ=sistema.params
        γc=sqrt(ω₀*ω)/2
        if γc>γ
            error("Aún no se ponen métodos para la fase normal")
            return
        end
        if signoQ==+
            signoq=-
        else
            signoq=+
        end
        return Point(Q=signoQ(sqrt(2-ω*ω₀/(2*γ^2))),q=signoq(sqrt(4*γ^2/ω^2 - ω₀^2/(4*γ^2))),P=0,p=0)
    end
    minimum_energy(sistema::ClassicalSystems.ClassicalSystem)= hamiltonian(sistema)(minimum_energy_point(sistema))
    function energy_width_of_coherent_state(sistema::ClassicalSystems.ClassicalSystem,punto,j::Real)
        ω₀,ω,γ=sistema.params
        Q,q,P,p=punto
        θ=PhaseSpaces.θ_of_QP(Q,P)
        ϕ=PhaseSpaces.ϕ_of_QP(Q,P)
        Ω₁=j*(ω^2*(q^2+p^2)/2+ω₀^2*sin(θ)^2/2 +2*γ^2*((sin( θ)^2*sin(ϕ)^2+cos(θ)^2)*q^2 + sin(θ)^2*cos(ϕ)^2) +2*γ*q*(ω*cos(ϕ) + ω₀*cos(θ)*cos(ϕ))*sin(θ))
        Ω₂=γ^2*(sin(θ)^2*sin(ϕ)^2 + cos(θ)^2)
        return sqrt(Ω₁ + Ω₂)/j
    end
    function discriminant_of_q_solution(sistema::ClassicalSystems.ClassicalSystem; Q,P,p,ε)
        ω₀,ω,γ=sistema.params
        r=1-(Q^2+P^2)/4
        return 4*γ^2*Q^2*r - p^2*ω^2 - ω*ω₀*(P^2 + Q^2 - 2) + 2*ω*ε
    end
    function q_of_ε(sistema::ClassicalSystems.ClassicalSystem;Q,P,p,ε,signo::Union{typeof(-),typeof(+)}=+,returnNaNonError=true)
        Δ=discriminant_of_q_solution(sistema;Q=Q,P=P,p=p,ε=ε)
        if Δ<0
            if returnNaNonError
                return NaN
            else
                εmin= minimum_ε_for(sistema;Q=Q,P=P,p=p)
                error("The minimum energy for the given parameters is $εmin. There is no q such that (Q=$Q,q,P=$P,p=$p) has ε=$ε.")
            end
        end
        ω₀,ω,γ=sistema.params
        return signo(-2*γ*Q*sqrt(1-(Q^2+P^2)/4),sqrt(Δ))/ω
    end
    function minimum_ε_for(sistema::ClassicalSystems.ClassicalSystem;Q=:nothing,q=:nothing,P=:nothing,p=:nothing)
        if count(x->x==:nothing,[q,Q,p,P])!=1
            error("You have to give two of Q,q,P,p")
        end
        ω₀,ω,γ=sistema.params
        if q==:nothing
            r=1-(Q^2+P^2)/4
            return -(4*γ^2*Q^2*r - p^2*ω^2 - ω*ω₀*(P^2 + Q^2 - 2))/(2*ω)
        end
        if Q==:nothing
            return (2*p^2*ω + 2*q^2*ω + P^2*ω₀- sqrt((-4 + P^2)^2* (4*q^2*γ^2 + ω₀^2)))/4
        end
        error("Sorry, that combinanation of variables is not implemented")
    end
    """maximum_P_for_ε(sistema::ClassicalSystems.ClassicalSystem,ε)
    Computes the maximum value of the parameter P accessible to the system at energy ε.
    """
    function maximum_P_for_ε(sistema::ClassicalSystems.ClassicalSystem,ε)
        ω₀,ω,γ=sistema.params
        if ε>=ω₀
            return 2.0
        elseif ω*ω₀>sqrt(2)*γ*sqrt(ω*(ω₀-ε))
            return sqrt(2*(ω₀+ε)/ω₀)
        else
            return sqrt(clamp((2*√(2*γ^6*ω*(ω₀-ε)))/(-γ^4) + ω*ω₀/(γ^2)+4,0,4))
        end
    end
    """maximum_Q_for_ε(sistema::ClassicalSystems.ClassicalSystem,ε)
    Computes the maximum value of the parameter Q accessible to the system at energy ε.
    """
    function maximum_Q_for_ε(sistema::ClassicalSystems.ClassicalSystem,ε)
        ω₀,ω,γ=sistema.params
        if ε>=ω₀
            return 2.0
        end
         return sqrt(clamp((4γ^2-ω*ω₀)/(2*γ^2) + (1/2)*sqrt(max(0,(16*γ^4+8γ^2*ε*ω+ω^2*ω₀^2)/γ^4)),0,4))
    end

    dicke_step = @ode_def begin
        dQ = integate_backwards*(P*ω₀-γ*P*q*Q/sqrt(4-P^2-Q^2)) # =dH/dP
        dq = integate_backwards*p*ω # =dH/dp
        dP = integate_backwards*(γ*q*Q^2/sqrt(4-P^2-Q^2)-γ*q*sqrt(4-P^2-Q^2)-Q*ω₀) # =-dH/dQ
        dp = integate_backwards*(-γ*Q*sqrt(4-P^2-Q^2)-q*ω) #-dH/dq

    end ω₀ ω γ integate_backwards
    function out_of_bounds(u,_,t)
        return 4-u[3]^2-u[1]^2<=0
    end
    varnames=[:Q,:q,:P,:p]
    ClassicalSystem(;ω₀,ω,γ)=ClassicalSystems.ClassicalSystem([ω₀,ω,γ],dicke_step,out_of_bounds,varnames)

    function Point(u::Array)
        u
    end
    Point(;Q,q,P,p)=Point([Q,q,P,p])
    Point(sistema::ClassicalSystems.ClassicalSystem;Q,P,p,ε,signo::Union{typeof(-),typeof(+)}=+) = Point([Q,q_of_ε(sistema::ClassicalSystems.ClassicalSystem,Q=Q,P=P,p=p,ε=ε,signo=signo,returnNaNonError=false),P,p])
    Pointθϕ(;θ,ϕ,q,p)=Point(Q=Q_of_θϕ(θ,ϕ),q=q,P=P_of_θϕ(θ,ϕ),p=p)


    function phase_space_dist_squared(x,y)
        Q1,q1,P1,p1=x
        Q2,q2,P2,p2=y
        return (q1-q2)^2 + (p1-p2)^2 + arc_between_QP(Q1,P1,Q2,P2)^2
    end

     function WignerHWxSU2_fixed_ε(sistemaC;Q,P,p,ε,j,signoq=+)
            Δ(Q,P)=discriminant_of_q_solution(sistemaC;Q=Q,P=P,p=p,ε=ε)
            x₀=Point(sistemaC,Q=Q,P=P,p=p,ε=ε,signo=signoq);
            dst=TruncatedWignerApproximation.WignerHWxSU2(x₀,j);
            function probability_density(u)
                Q,q,P,p=u
                return dst.probability_density(Point(sistemaC,Q=Q,P=P,p=p,ε=ε,signo=signoq))/sqrt(Δ(Q,P))
            end
            function logf(u)
                θ,ϕ,p=u
                Q=PhaseSpaces.Q_of_θϕ(θ,ϕ)
                P=PhaseSpaces.P_of_θϕ(θ,ϕ)
                try
                    return log(probability_density([Q,NaN,P,p]))#*sin(θ)?
                catch
                    return -Inf
                end
              end

            m = RWMVariate([PhaseSpaces.θ_of_QP(Q,P),PhaseSpaces.ϕ_of_QP(Q,P),p], [1/sqrt(j),1/sqrt(j),1/sqrt(j)],logf,proposal = Normal)
            function sample()
                    θ,ϕ,p=sample!(m)
                    return Point(sistemaC,Q=PhaseSpaces.Q_of_θϕ(θ,ϕ),P=PhaseSpaces.P_of_θϕ(θ,ϕ),p=p,ε=ε,signo=signoq)
            end
            for i in 1:10000
                sample!(m) #bake
            end
            TruncatedWignerApproximation.PhaseSpaceDistribution(probability_density,sample)
        end

        #esta funcion no funciona para la fase normal
         function classicalPathRandomSampler(sistema::ClassicalSystems.ClassicalSystem;ε,dt=3)
             ω₀,ω,γ=sistema.params
             s=sqrt(16*γ^4 + 8*γ^2*ε*ω + ω^2*ω₀^2)/(2*γ^2)
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
             u0=ClassicalDicke.Point(sistema,Q=Q,P=0,p=0,ε=ε)
             function sample()
                 u0=ClassicalSystems.integrate(sistema; t=dt, u₀=u0, save_intermediate_steps=false,show_progress=false).u[end]
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
