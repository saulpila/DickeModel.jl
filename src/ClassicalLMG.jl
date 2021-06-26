module ClassicalLMG

export hamiltonian, Point,ClassicalSystem
    import ..ClassicalSystems
    using ..PhaseSpaces
    using ParameterizedFunctions

    function hamiltonian(sistema::ClassicalSystems.ClassicalSystem)
        α,=sistema.params
        function H(u)
            Q,P=u
            (Q^2 + P^2)/2 + α*(Q^2 - P^2*Q^2/4 - Q^4/4)
        end
        return H
    end

    function discriminant_of_P_solution(sistema::ClassicalSystems.ClassicalSystem,Q,ε)
        α,=sistema.params
        return (-Q^4*α+Q^2*(2+4*α)-4*ε)/(-2+Q^2*α)
    end
    function P_of_ε(sistema::ClassicalSystems.ClassicalSystem;Q,ε,signo::Union{typeof(-),typeof(+)}=+,returnNaNonError=true)
        Δ=discriminant_of_P_solution(sistema,Q,ε)
        if Δ<0
            if returnNaNonError
                return NaN
            else
                εmin= minimum_ε_for(sistema;Q=Q)
                error("The minimum energy for the given parameters is $εmin. There is no P such that (Q=$Q,P) has ε=$ε.")
            end
        end
        α,=sistema.params
        return signo(0.0,sqrt(Δ))
    end
    function minimum_ε_for(sistema::ClassicalSystems.ClassicalSystem;Q=:nothing,P=:nothing)
        if count(x->x==:nothing,[Q,P])!=1
            error("Tienes que una una de Q,P")
        end
        α,=sistema.params
        if P==:nothing
            return (2*Q^2 + 4*Q^2*α - Q^4*α)/4
        end
        error("Formula no implementada")
    end
    lipkin_step = @ode_def begin
        dQ = integate_backwards*P*(1 - α*Q^2/2) # =dH/dP
        dP = integate_backwards*(Q*(α*P^2/2 - (2*α + 1)) + α*Q^3) # =-dH/dQ
    end α integate_backwards
    function out_of_bounds(u,_,t)
        return 4-u[2]^2-u[1]^2<=0
    end
    ClassicalSystem(;α)=ClassicalSystems.ClassicalSystem([α],lipkin_step,out_of_bounds,[:Q,:P])

    function Point(u::Array)
        u
    end
    Point(;Q,P)=Point([Q,P])
    Pointθϕ(;θ,ϕ)=Point(Q=Q_of_θϕ(θ,ϕ),P=P_of_θϕ(θ,ϕ))
    Point(sistema::ClassicalSystems.ClassicalSystem;Q,ε,signo::Union{typeof(-),typeof(+)}=+)=Point(Q=Q,P=P_of_ε(sistema;Q=Q,ε=ε,signo=signo,returnNaNonError=false))
end
