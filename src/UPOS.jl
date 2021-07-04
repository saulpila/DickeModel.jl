module UPOS

export search_in_interval,find_orbit,get_period,follow_UPO_family,find_p_0,fix_UPO_keep_P,fix_UPO_keep_energy
    using ..ClassicalSystems
    using ..ClassicalDicke
    using LinearAlgebra
    using DiffEqCallbacks,OrdinaryDiffEq, DiffEqBase
    using ProgressMeter
    using ..PhaseSpaces
    using ..DickeBCE
    using ..DickeHusimiProjections
    Id=Matrix{Float64}(I, 4, 4);
    struct PO
        sistema::ClassicalDickeSystem
        u
        T
    end

    function Base.show(io::IO, m::PO)
         print(io,string("PO @ u=",m.u,", T=", m.T))
    end
    PO(sistema::ClassicalDickeSystem,u)=PO(sistema,u,get_period(sistema,u))
    integrate_PO(po::PO;tol=1e-16)=ClassicalSystems.integrate(po.sistema;t=po.T,u₀=po.u,tol=tol);
    function action(po::PO;tol=1e-16)
        us=integrate_PO(po,tol=tol).u
        return sum(us[i][3]*(us[i][1]-us[i-1][1]) + us[i][4]*(us[i][2]-us[i-1][2]) for i in 2:length(us))
    end
    function average_over_PO(po::Union{OrdinaryDiffEq.ODESolution,PO},f;tol=1e-12)
        tot=nothing

        if isa(po,PO)
            function guardar(u,t,integrator)
                v=f(u)*(t-integrator.tprev)/po.T
                if tot===nothing
                    tot=v
                else
                    tot+=v
                end
            end
            cb=FunctionCallingCallback(guardar;func_start = false,func_everystep =true)
            ClassicalSystems.integrate(po.sistema;t=po.T,u₀=po.u,tol=tol,save_everystep=false,callback=cb)
        else
            T=po.t[end]-po.t[1]
            for i in 2:length(po.t)
                v=f(po.u[i])*(po.t[i]-po.t[i-1])/T
                if tot===nothing
                    tot=v
                else
                    tot+=v
                end
            end
        end
        return tot
    end

    average_over_PO_QP(PO::PO,f;tol=1e-12)=average_over_PO(PO,x->f(x[1],x[3]);tol=tol)
    average_over_PO_qp(PO::PO,f;tol=1e-12)=average_over_PO(PO,x->f(x[2],x[4]);tol=tol)

    jz_PO_average(PO::PO;tol=1e-12)=average_over_PO_QP(PO,PhaseSpaces.jz;tol=tol)
    jy_PO_average(PO::PO;tol=1e-12)=average_over_PO_QP(PO,PhaseSpaces.jy;tol=tol)
    jx_PO_average(PO::PO;tol=1e-12)=average_over_PO_QP(PO,PhaseSpaces.jx;tol=tol)
    mirror_PO(po::PO)=PO(po.sistema,po.u.*[-1,-1,1,1],po.T)

    function find_orbit(sistema,u₀;bound=0.1,n=10,tol=1e-8,max_order=1,knownorder=0)
      #  @show Q
       # @show q
        counts=n
        recurrences=0
        boundexceeded=false
        count=0
        notreturnedcount=0
        lastsavet=0
        function guardar(integrator)
            count+=1
            c=integrator.u[:]

            counts-=1
            returned=bound>norm(c-u₀)
            if !returned
                notreturnedcount+=1
            else
                if knownorder == 0
                    knownorder=notreturnedcount+1
                    notreturnedcount=0
                    recurrences+=1
                    lastsavet=integrator.t
                else
                    if knownorder==notreturnedcount+1
                        recurrences+=1
                        lastsavet=integrator.t
                        notreturnedcount=0
                    else
                        terminate!(integrator)
                        return nothing
                    end
                end
            end
            if knownorder==0 && notreturnedcount>=max_order
                terminate!(integrator)
                return nothing
            end
            if knownorder*recurrences>n
                terminate!(integrator)
                return nothing
            end
            u_modified!(integrator,false)

            return nothing

        end
       # try
            cb=ContinuousCallback((uvar,t,integrator)-> (uvar[4]-u₀[4]),nothing,guardar,save_positions=(false,false),rootfind=true,interp_points=3,reltol=10^-8,abstol=10^-8)
            t=ClassicalSystems.integrate(sistema;t=10000,u₀=u₀,callback=cb,save_everystep=false,tol=tol).t[end]
            return knownorder,recurrences,lastsavet/recurrences

       # catch
       #     return NaN,NaN
      #  end
    end

    function get_period(sistema,u₀;tol=1e-8,othervartol=1e-2)
        t=0
        count=0
        function guardar(integrator)
            if maximum(abs.(integrator.u - u₀))<othervartol
                t=integrator.t

                terminate!(integrator)
                return nothing
            end
        end
        try
        cb=ContinuousCallback((uvar,t,integrator)-> uvar[end]-u₀[end],nothing,guardar,save_positions=(false,false),rootfind=true,interp_points=3,reltol=tol,abstol=tol)
            u=ClassicalSystems.integrate(sistema;t=10000,u₀=u₀,callback=cb,save_everystep=false,tol=tol).u[end]
            return t
        catch
            return NaN
        end
    end

    function findpeaks(list)
       inds = []
       prev=0
       for i in 2:length(list)
           if prev>0
               if list[i]>list[prev]
                    prev=i
                elseif  list[i] <list[prev]
                    length=i-prev
                    push!(inds,(Int(floor((prev+i)/2)),length))
                    prev=0
                end
            else
                if list[i]>list[i-1]
                    prev=i
                end
            end
         end
        return inds
    end
    closestodd(i)=Int(2*ceil((i+1)/2)-1)
    function search_in_interval(sistema,ϵ,Q=1,s=10000,ds=2;tol=1e-7,n=20,_start=true,refine=false)
        function f(Q,P,bound)
            try
                return find_orbit(sistema,ClassicalDicke.Point(sistema,P=P,Q=Q,p=0.0,ϵ=ϵ),n=n,tol=1e-6,bound=bound,max_order=3)
            catch
                return NaN,NaN
            end
        end
        found=[]
        Qs=range(Q-ds/2,stop=Q+ds/2,length=closestodd(s))
        pts=f.(Qs,0,0.1)
        peaks=findpeaks([pt[2] for pt in pts])
        if _start
            display(length(peaks))
            if !refine
                return [(Qs[i],pts[i][3]) for (i,size) in peaks]
            end
        end
        for (i,size) in peaks

            Qi=Qs[i]
            w=(ds*(size+1)/length(Qs))
            if w<tol
                if 4<pts[i][2]
                    display((Qi,pts[i][1],pts[i][2]))

                    push!(found,(Qi,pts[i][1],pts[i][2]))
                end
            else
                found=[found;searchInterval(Qi,100,w,tol=tol,_start=false)]
            end
        end
        for pt in found[:]
            for other in found[:]
                if other==pt
                    continue
                end

                if abs(pt[1]-other[1])<tol
                    to_remove=pt
                    if other[3]<pt[3]
                        to_remove=other
                    end
                    try
                        deleteat!(found,findfirst(x->x==to_remove,found))

                    catch
                    end

                end
            end
        end
        return found
    end
    function next(sistema,u₀,period;tol=1e-12)
        u₁,M=ClassicalSystems.integrate(sistema,t=period,u₀=u₀,save_everystep=false,get_fundamental_matrix=true,tol=tol).u[end].x
        return u₀-(M-Id)^-1*(u₁-u₀),norm(u₁-u₀)
    end
    function hamiltonian_gradient(sistema,u)
        Fu₁=[0.0,0,0,0]
        ClassicalSystems.step(sistema)(Fu₁,u,[ClassicalDickeSystem.parameters(sistema);1.0],1.0)
        gradH=[-Fu₁[3],-Fu₁[4],Fu₁[1],Fu₁[2]]
        return gradH
    end

    function next_constant_energy(sistema,u₀,T;tol=1e-12)
        u₁,M=ClassicalSystems.integrate(sistema,t=T,u₀=u₀,save_everystep=false,get_fundamental_matrix=true,tol=tol).u[end].x
        gradH=hamiltonian_gradient(sistema,u₁)
        ξ=[0.0,0,0,1]
        Fu₁=-[-gradH[3],-gradH[4],gradH[1],gradH[2]]
        Mat=[M-Id Fu₁;
            transpose(gradH) 0;
            transpose(ξ) 0]
        res=pinv(Mat)*[u₀-u₁;0;0]
        du=res[1:4]
        dT=res[5]
        return u₀+du,T+dT,norm(u₁-u₀)
    end

    function next_constant_P(sistema,u₀,T;tol=1e-12)
        u₁,M=ClassicalSystems.integrate(sistema,t=T,u₀=u₀,save_everystep=false,get_fundamental_matrix=true).u[end].x
        ξ=[0.0,0,0,1]
        gradH=hamiltonian_gradient(sistema,u₁)
        Fu₁=[-gradH[3],-gradH[4],gradH[1],gradH[2]]
        Mat=[M-Id Fu₁;
            transpose(ξ) 0]
        res=(Mat^-1)*[u₀-u₁;0]
        du=res[1:4]
        dT=res[5]
        return u₀+du,T+dT,norm(u₁-u₀)
    end
    function fix_UPO(sistema,u₀,period;maxiters=100,tol=1e-8,inttol=tol/10)
    u₁=u₀
    for j in 1:maxiters
        u₁,Δ=next(sistema,u₁,period,tol=inttol)
        if Δ<tol
            break
        end

    end
    u₁
    end

    function fix_UPO_keep_energy(sistema,u₀,period;maxiters=100,tol=1e-8,inttol=tol/10)
    u₁=u₀
    T=period

    for j in 1:maxiters
        u₁,T,Δ=next_constant_energy(sistema,u₁,T;tol=inttol)
        if T<0
            error("Error de convergencia")
        end
        if Δ<tol
            break
        end
        if(j==maxiters)
            @show "MAXITERS"
        end

    end
    u₁,T
    end

    function fix_UPO_keep_P(sistema,u₀,period;maxiters=100,tol=1e-8,inttol=tol/10)
    u₁=u₀
    T=period
    for j in 1:maxiters
        u₁,T,Δ=next_constant_P(sistema,u₁,T,tol=1e-8)
        if Δ<tol
            break
        end

    end
    u₁,T
    end

    function find_p_0(sistema,u₀;tol=1e-6,negative=false)
        nu=[0,0,0,0]

        function guardar(integrator)

            nu=integrator.u

            terminate!(integrator)
            return nothing
        end
        try

            negguardar=nothing
            if negative
                negguardar,guardar=guardar,negguardar
            end
            cb=ContinuousCallback((uvar,t,integrator)-> (uvar[4]-0),guardar,negguardar,save_positions=(false,false),rootfind=true,interp_points=3,abstol=tol)
            ClassicalSystems.integrate(sistema;t=10000,u₀=u₀,callback=cb,save_everystep=false,tol=1e-6)
            return nu
        catch
            if abs(u₀[4])< tol  #hay un error cuando el punto ya satisface la condicion por juliaDiff
                return u₀
            end
            return NaN
        end
    end

    function follow_UPO_family(sistema,u₀,period,finalperiod;step=0.1,savestep=:nothing,tol=1e-8)
        us=[]
        if savestep==:nothing
            savestep=step
        end
        if  savestep%step!=0
            error("savestep debe ser un multiplo de step")
        end
        periods=period:(period<finalperiod ? 1 : -1)*step:finalperiod
        saveperiods=period:(period<finalperiod ? 1 : -1)*savestep:finalperiod
        newu=u₀
        prog=Progress(length(periods), 1)
        for p in periods
            next!(prog)
            newu=fix_UPO(sistema,newu,p,tol=tol,inttol=min(1e-12,tol/100))
            newu=find_p_0(sistema,newu,tol=tol)
            if p in saveperiods
                push!(us,newu)
            end
            realp=get_period(sistema,newu,tol=tol)
            @show newu
            @show p,realp
        end
        finish!(prog)
        return [PO(sistema,us[i],saveperiods[i]) for i in 1:length(saveperiods)]
    end
     function lyapunov_exponent(PO)
         sistema=PO.sistema
         T=PO.T
         maximum(log.(abs.(
            LinearAlgebra.eigvals(ClassicalSystems.integrate(sistema,u₀=PO.u,t=T,
                    get_fundamental_matrix=true,save_everystep=false)[end].x[2])
            ))/T)
    end
    function _max0(x)
        if x>0
            return x
        else
            return Inf
        end
    end

    function follow_UPO_family_energies(sistema,u₀,period,H;step=0.1,tol=1e-5,energytol=1e-3,initalperturbation=[0,1,0,0],dontallowlessenergy=false)
        us=[u₀]


        E₀=H(u₀)

        energies=Float64[H(u₀)]
        periods=Float64[period]

        cacheJ=ClassicalSystems.eye(4)
        Λ=[0 0 -1 0
           0 0 0 -1
           1 0 0 0
           0 1 0 0]





        function get(finalEnergy;tol=tol,energytol=energytol,minstep=energytol/10,maxstep=0.1,p0positiveside=false,step=maxstep,_prog=ProgressThresh(0.0,"Remaining Energy"),maxiters=1000,_it=0,energywiggletolerance=1e-2)
            index =NaN
            prog=_prog
            if dontallowlessenergy
                finalEnergy=max(finalEnergy,E₀)
            end
            while true #abs(energies[index]-finalEnergy)>energytol
                _it+=1
                if _it>maxiters
                    error("Max iters")
                end
                #index=findmin(abs.(energies.-finalEnergy))[2]
                index=findmin(_max0.(finalEnergy.-energies))[2]
                last_e=energies[index]
                falta=abs(last_e-finalEnergy)
                if falta<energytol
                    break
                end
                step=min(step,falta)
                direction=sign(finalEnergy-last_e)

                newu=us[index]
                prog.desc="Remaining Energy (step=10^$(round(log10(step),digits=1))"

                converror=false
                try
                    ClassicalSystems.step(sistema).jac(cacheJ,newu,[ClassicalSystems.parameters(sistema);1.0],1.0)
                    grad=hamiltonian_gradient(sistema,newu)
                    hessian=Λ*cacheJ
                    if index==1
                        Δu=initalperturbation
                    else
                    #   Δu=us[index]-us[index-1]
                    #   newu=find_p_0(sistema,us[index]-,tol=tol,negative=!p0positiveside)
                      #  u1=find_p_0(sistema,us[index-1],tol=tol,negative=!p0positiveside)
                     #   periodMultiple=0.333
#
                      #  newu=ClassicalSystems.integrate(sistema,u₀=newu,t= periods[index]*0.1,tol=tol,save_everystep=false).u[end]
                     #   u1=ClassicalSystems.integrate(sistema,u₀=u1,t=periods[index-1]*periodMultiple,tol=tol,save_everystep=false).u[end]
                       # Δu= newu - u1
                       Δu=hamiltonian_gradient(sistema,newu)
                    end
                #    eigenvals,eigenvects=LinearAlgebra.eigen(log(ClassicalSystems.integrate(sistema,u₀=newu,t=periods[index],
                #    get_fundamental_matrix=true,save_everystep=false,tol=tol)[end].x[2]))
                 #   nullvects=nullspace(transpose(real.(eigenvects)), atol=0.01)
                 #   allnullvects=[nullvects[:,i]/norm(nullvects[:,i]) for i in 1:size(nullvects)[2]]
                   # minp=minimum(abs(v[4]) for v in allnullvects)
                  #  Δu=-allnullvects[findfirst(x->abs(x[4])== minp,allnullvects)]


                    Δu=Δu/norm(Δu)
                  #  @show Δu,allnullvects
                    ΔE=direction*step


                    a=0.5*transpose(Δu)*hessian*Δu
                    b=dot(grad,Δu)
                    c=- ΔE


                    #hacemos Taylor
                    #-c=Δϵ=b*s + a s^2
                    s=(-b+sqrt(b^2-4*a*c))/(2*a)
                    if abs(H(newu + s*Δu) - (last_e + ΔE)) >abs(H(newu - s*Δu) - (last_e + ΔE))
                        s=-s
                    end
                    lastu=newu
                    E=H(newu+s*Δu)
                    periodGuess=periods[index]
                    if index>1
                        periodGuess+= (E-energies[index-1])*(periods[index]-periods[index-1])/(energies[index]-energies[index-1])
                    end

                    newu,period= fix_UPO_keep_energy(sistema,newu+s*Δu,periodGuess,inttol=min(1e-12,tol/100),tol=tol)

                    if(E-H(newu))>energywiggletolerance
                        if index>1
                            deleteat!(energies,index)
                            deleteat!(us,index)
                            deleteat!(periods,index)
                        end
                        converror=true
                        index-=1
                    end
                  #  if sign(finalEnergy-H(newu))!=direction
                  #      if step/2<minstep
                  #          converror=true
                  #       else
                  #          step=step/4
                  #      end
                  #  end
                catch e
                    #error("ij")
                    converror=true
                end

                notbetter=abs(finalEnergy-H(newu))> falta # && sign(finalEnergy-H(newu))==direction
               #@show notbetter,H(newu),energies[index]
                if  converror || H(newu) in energies || notbetter
                #    if notbetter
              #          @show index, us[index]
             ##           us[index]=ClassicalSystems.integrate(sistema,u₀=newu,t= periods[index]*0.1,tol=tol,save_everystep=false).u[end]
              #          @show , us[index]
              #      end
                    if step/1.2<minstep# (!converror && H(newu) in energies) ||

                        if p0positiveside==true
                            p0positiveside=false
                            step=maxstep*1.1

                        else
                            cancel(prog)
                            error("El algoritmo no está convergiendo. Llegó hasta ",PO(sistema,us[index],periods[index]), " y energía ",H(us[index]) )
                        end
                    end
                    step/=1.2
                    continue
                    #return get(finalEnergy;step=step/1.1,tol=tol,energytol=energytol,p0positiveside=p0positiveside,minstep=minstep,maxstep=maxstep,_prog=prog,maxiters=maxiters,_it=_it)


                end
                update!(prog,abs(last_e-finalEnergy))

                step=min(maxstep,step*1.1)

                nindex=findmin(_max0.(H(newu).-energies))[2]
                insert!(us,nindex+1,newu)
                insert!(energies,nindex+1,H(newu))
                insert!(periods,nindex+1,period)
                index=nindex+1

            end
            finish!(prog)
            return PO(sistema,us[index],periods[index])#,energies[index]
        end
        return get
    end
    family_A(sistemaC::ClassicalDickeSystem)=follow_UPO_family_energies(sistemaC,ClassicalDicke.minimum_energy_point(sistemaC,+),2*pi/ClassicalDicke.normal_frequency(sistemaC,+),ClassicalDicke.hamiltonian(sistemaC),tol=1e-8,energytol=1e-3)
 family_B(sistemaC::ClassicalDickeSystem)=follow_UPO_family_energies(sistemaC,ClassicalDicke.minimum_energy_point(sistemaC,+),2*pi/ClassicalDicke.normal_frequency(sistemaC,-),ClassicalDicke.hamiltonian(sistemaC),tol=1e-8,energytol=1e-3)
    function plot_PO_QP(sis::ClassicalDickeSystem,F::Array{PO,1};p=plot(),H=:nothing,opts...)
        for po in F

            u=integrate_PO(po)
            us=u.u
            ts=u.t
            plot!(p,[(i[1],i[3]) for i in us];xlabel="Q",ylabel="P",fontfamily="Times",opts...)
        end
        if H!=:nothing
            contour!(range(-2,stop=2,length=100),range(-2,stop=2,length=100),(Q,P)->ClassicalDicke.minimum_ϵ_for(sis;Q=Q,P=P,p=0),levels=[maximum(H(f.u) for f in F)],linewidth=1,linecolor=:black)
        end
        return p
    end
    # function plot_PO_QPp(sis::ClassicalDickeSystem,F::Array{PO,1};p=plot(),H=:nothing,opts...)
    #     for po in F
    #
    #         u=integrate_PO(po)
    #         us=u.u
    #         ts=u.t
    #         plot!(p,[i[1] for i in us],[i[3] for i in us],[i[4] for i in us];xlabel="Q",ylabel="P",zlabel="p",fontfamily="Times",opts...)
    #     end
    #     return p
    # end
    # function plot_PO_qp(sis::ClassicalDickeSystem,F::Array{PO,1};H=:nothing,p=plot(),opts...)
    #     for po in F
    #
    #         u=integrate_PO(po)
    #         us=u.u
    #         ts=u.t
    #         plot!(p,[(i[2],i[4]) for i in us];xlabel="q",ylabel="p",fontfamily="Times",opts...)
    #     end
    #     if H!=:nothing
    #
    #         contour!(range(-4,stop=4,length=100),range(-2,stop=2,length=100),(q,p)->ClassicalDicke.minimum_ϵ_for(sis;q=q,p=p,P=0),levels=[maximum(H(f.u) for f in F)],linewidth=1,linecolor=:black)
    #     end
    #     return p
    # end
    # function plot_PO_Qq(sis::ClassicalDickeSystem,F::Array{PO,1},H;p=plot(),opts...)
    #     for po in F
    #
    #         u=integrate_PO(po)
    #         us=u.u
    #         ts=u.t
    #         plot!(p,[(i[1],i[2]) for i in us];xlabel="Q",ylabel="q",fontfamily="Times",opts...)
    #     end
    #     return p
    # end
    function overlap_of_tube_with_homogenous_state(sistemaQ::DickeBCE.QuantumSystem,po::PO;time_integral_tolerance=1e-7,phase_space_integral_resolution=0.1)
        orbit=integrate_PO(po,tol=time_integral_tolerance)
        ∫dtHtx(x)=average_over_PO(orbit,u-> DickeBCE.HusimiOfCoherent(sistemaQ,u,x))
        res=phase_space_integral_resolution
        ϵ=ClassicalDicke.hamiltonian(po.sistema)(po.u)
        insidepoints=0
        function f(Q,P)
            v=DickeHusimiProjections.∫∫dqdpδϵ(sistemaC=po.sistema,ϵ=ϵ,Q=Q,P=P,f=∫dtHtx,nonvalue=NaN,p_res=res)
            if isnan(v)
                return 0.0
            end
            insidepoints+=1
            return v
        end
        sm=sum( f(Q,P) for Q in -2:res:2, P in -2:res:2)
        return sm/(insidepoints*2*pi)
    end
    function scarring_measure(sistemaQ::DickeBCE.QuantumSystem,quantum_state::AbstractArray{<:Number,1},po::PO;Htol=1e-3,kargs...)
        POandState=average_over_PO(po,x->DickeBCE.Husimi(sistemaQ,x,quantum_state;tol=Htol))
        POandHomState=overlap_of_tube_with_homogenous_state(sistemaQ,po;phase_space_integral_resolution=0.1,kargs...)
        return POandState/POandHomState
    end
    function Base.in(u::AbstractArray{<:Number,1}, po::PO)
        tol=1e-6
        if sum(abs2,u-po.u)<tol
            return true
        end

        result=false

        function distance(uvar,t,integrator)
            m=sum(abs2,uvar-u)
            if m<tol
                result=true
                terminate!(integrator)
            end
            return m*sign(uvar[4]-u[4])
        end
        cb=ContinuousCallback(distance,integrator->distance(integrator.u,integrator.t,integrator),save_positions=(false,false),rootfind=true,interp_points=10,reltol=10^-8,abstol=10^-8)
        ClassicalSystems.integrate(po.sistema;t=po.T,u₀=po.u,tol=1e-8,save_everystep=false,callback=cb)
        return result
    end
    function Base.:(==)(po1::PO,po2::PO)
        if abs(po1.T-po2.T)>1e-5
            return false
        end
        return po1.u ∈ po2
    end
end
