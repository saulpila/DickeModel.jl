module ClassicalLMG

export hamiltonian, Point,ClassicalLMGSystem
    import ..ClassicalSystems
    using ..PhaseSpaces
    using ParameterizedFunctions

    LMG_step = @ode_def begin
        dQ = integate_backwards*P*(1 - α*Q^2/2) # =dH/dP
        dP = integate_backwards*(Q*(α*P^2/2 - (2*α + 1)) + α*Q^3) # =-dH/dQ
    end α integate_backwards
    function out_of_bounds_LMG(u,_,t)
        return 4-u[2]^2-u[1]^2<=0
    end
    varnames=[:Q,:P]
    struct ClassicalLMGSystem <: ClassicalSystems.ClassicalSystem
        parameters::Vector{Float64}
        ClassicalLMGSystem(;α::Real)= new([Float64(α)])

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
    
    function hamiltonian(sistema::ClassicalLMGSystem)
        α,=ClassicalSystens.params(sistema)
        function H(u)
            Q,P=u
            (Q^2 + P^2)/2 + α*(Q^2 - P^2*Q^2/4 - Q^4/4)
        end
        return H
    end

    function discriminant_of_P_solution(sistema::ClassicalLMGSystem,Q,ϵ)
        α,=ClassicalSystens.params(sistema)
        return (-Q^4*α+Q^2*(2+4*α)-4*ϵ)/(-2+Q^2*α)
    end
    function P_of_ϵ(sistema::ClassicalLMGSystem;Q,ϵ,signo::Union{typeof(-),typeof(+)}=+,returnNaNonError=true)
        Δ=discriminant_of_P_solution(sistema,Q,ϵ)
        if Δ<0
            if returnNaNonError
                return NaN
            else
                ϵmin= minimum_ϵ_for(sistema;Q=Q)
                error("The minimum energy for the given parameters is $ϵmin. There is no P such that (Q=$Q,P) has ϵ=$ϵ.")
            end
        end
        α,=ClassicalSystens.params(sistema)
        return signo(0.0,sqrt(Δ))
    end
    function minimum_ϵ_for(sistema::ClassicalLMGSystem;Q=:nothing,P=:nothing)
        if count(x->x==:nothing,[Q,P])!=1
            error("Tienes que una una de Q,P")
        end
        α,=ClassicalSystens.params(sistema)
        if P==:nothing
            return (2*Q^2 + 4*Q^2*α - Q^4*α)/4
        end
        error("Formula no implementada")
    end


    function Point(u::Array)
        u
    end
    Point(;Q,P)=Point([Q,P])
    Pointθϕ(;θ,ϕ)=Point(Q=Q_of_θϕ(θ,ϕ),P=P_of_θϕ(θ,ϕ))
    Point(sistema::ClassicalLMGSystem;Q,ϵ,signo::Union{typeof(-),typeof(+)}=+)=Point(Q=Q,P=P_of_ϵ(sistema;Q=Q,ϵ=ϵ,signo=signo,returnNaNonError=false))
end


