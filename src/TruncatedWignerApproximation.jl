module TruncatedWignerApproximation
    export calculate_distribution, muestreoWignerSU2,muestreoWignerHW,muestreoWignerHWxSU2,functions_for_survival_probability,survival_probability

    using Distributed
    using Distributions
    using ProgressMeter
    using DiffEqCallbacks
    using SymEngine
    using ..ClassicalSystems
    using DataStructures
    using ..PhaseSpaces
    function scale_to_index(v,r::LinRange)
            if v<r.start || r.stop <v
                return NaN
            end
            if r.len==1
                return 1
            end
            return Int(round(
                    ((v-r.start)/(r.stop-r.start))*(r.len-1)
                    +1))

        end
    function monte_carlo_over_distribution(sistema::ClassicalSystems.ClassicalSystem,f_inicial,f_loop,reductora,finalizadora;ts=[0.0],funcion_muestreo,N,ignorar_errores=true,maxNBatch=Inf,kargs...)
        if ts==:nothing #algunas funciones pasan nothing aaaaah
           ts =[0.0]
        end
        show_progress= nothing
        trabajadores=length(workers())-1
        if trabajadores==0
            trabajadores=1
        end
        Neff=min(Int(ceil(N/trabajadores)),maxNBatch)
        workerload=[]
        rest=N
        while rest>0
            load=Neff
            rest=rest-load
            extra=-min(0,rest)
            rest=max(rest,0)
            push!(workerload,load-extra)
        end
        channel = RemoteChannel(()->Channel{Int8}(N), 1)
        res=0

        workersready=0
        missedprogressloop=trabajadores
        @sync begin
            # this task prints the progress bar
            @async while true
                v= take!(channel)
                if v==0
                    if show_progress!=nothing
                        finish!(show_progress)
                    end
                    break
                elseif v==1
                    if workersready>=trabajadores
                        if show_progress==nothing
                            show_progress=Progress(Int(N-missedprogressloop), 1)
                        end
                        next!(show_progress)
                    else
                        missedprogressloop+=1
                    end

                else
                    workersready+=1
                end
            end

            # this task does the computation
            @async begin
               try
                   res=@distributed reductora for load in workerload
                        val=f_inicial()
                        ti=1
                        ti_error=1
                        nerrores=0
                        call_loop= (u,t,_)-> begin
                            if ti>=ti_error
                                f_loop(val,u,ti)
                            end
                            ti+=1
                        end


                        j=1
                        notimeevolution= (ts==[0.0])
                        cb=:nothing
                        if !notimeevolution
                             cb=FunctionCallingCallback(call_loop;funcat=ts,func_start = false)
                        end
                        while j <=load
                            u₀=funcion_muestreo()

                            ti=1
                            try
                                if notimeevolution
                                    call_loop(u₀,0.0,0)
                                else
                                    ClassicalSystems.integrate(sistema,t=last(ts),save_everystep=false,u₀=u₀,callback=cb;kargs...)
                                end
                                ti_error=1

                                put!(channel,  j==1 ? 2 : 1 )
                                j+=1
                                nerrores=0

                            catch e
                                if ignorar_errores==false
                                    error(e)
                                end
                                if nerrores>100
                                    error("Demasiados errores")
                                end
                                nerrores+=1
                                ti_error=ti
                            end

                        end
                        val
                    end


                finally
                    put!(channel, 0) # this tells the printing task to finish

                end
            end
        end

        return finalizadora(res)
    end
    function encadenar(Fs)
        f_inicial=function()
            [F[1]() for F in Fs]
        end
        f_loop=function(valores,u,t)
            for i in 1:length(Fs)
                Fs[i][2](valores[i],u,t)
            end
        end
        reductora=function(last,next)
             [Fs[i][3](last[i],next[i]) for i in 1:length(Fs)]

        end
        f_final=function(res)
            [Fs[i][4](res[i]) for i in 1:length(Fs)]
        end
        return f_inicial,f_loop,reductora,f_final
    end
    function function_from_expression(sistema::ClassicalSystems.ClassicalSystem,expression)
        varnames=sistema.varnames
        if typeof(expression)!=SymEngine.Basic && !isempty(methods(expression))#es una funcion
            if any(length(m.sig.parameters)-1==length(varnames) for m in methods(expression))
                return expression
            elseif any(length(m.sig.parameters)-1==1 for m in methods(expression))
                f(x...)=expression(x)
                return f
            else
                error("Debes pasar una funcion que tome un vector, o n argumentos donde n es la dimension del espacio fase")
            end

        end

        if !(Symbol.(SymEngine.free_symbols(SymEngine.Basic(expression))) ⊆ varnames)
            error("La expresión $(expression) tiene variables desconocidas. Escribe todo en términos de $(sistema.varnames)")
        end
        if  typeof(expression)==SymEngine.Basic
            expression=SymEngine.toString(expression) #esto no es lo mejor, seria mejor usar expresiones y no strings pero algo pasa con las raices cuadradas...

        end
        return Base.eval(Main,Meta.parse("($(join(varnames,',')))->($expression)"))
    end
    function functions_from_expressions(sistema::ClassicalSystems.ClassicalSystem,expressions)
        funciones_varianza=[function_from_expression(sistema,s) for s in expressions]
        return funciones_varianza
    end
    function functions_for_averaging(sistema::ClassicalSystems.ClassicalSystem;expr,ts,N)
        onevar=false

        if !isa(expr,Array)
            expr=[expr]
            onevar=true
        end

        lengthexpr=length(expr)
        f_inicial= function()
            [zeros(lengthexpr) for i in 1:length(ts)]
        end
        varnames=sistema.varnames
        funciones_a_promediar=functions_from_expressions(sistema,expr)
        f_loop=function(valores,u,i)

            for fi in 1:lengthexpr
                valores[i][fi]+=funciones_a_promediar[fi](u...)
            end
        end
        f_final=function(res)
            res=res/N
            if onevar
                res=[re[1] for re in res]
            end
            return res
        end

        return f_inicial,f_loop,+,f_final
    end
    function functions_for_calculating_distributions(sistema::ClassicalSystems.ClassicalSystem;x,y=:nothing,z=:nothing,zs=:nothing,ts=:nothing,xs=-2:0.1:2,ys=:nothing,N)
        onematrix=false
        #vars=[]
        if ts==:nothing && y==:nothing
            y=:t
            ys=0:0
        end
        xs=LinRange(xs)
        ys=LinRange(ys)
        if zs!=:nothing
            zs=LinRange(zs)
        end
        if x==:t || y==:t
            onematrix=true
            if ts!=:nothing
                error("No usar el parametro ts al graficar el tiempo en un eje, usar el parametro del eje correspondiente")
            end

            if x==:t
                ts=xs
                #vars=ClassicalSystems.nvars(sistema,y)
            else
                ts=ys
               # vars=ClassicalSystems.nvars(sistema,x)
            end

      #  else
         #   vars=ClassicalSystems.nvars(sistema,x,y)
        end
        f_inicial= function()
            if onematrix
                matrices=zeros(length(ys),length(xs))
            else
                if zs!=:nothing
                    #matrices=[zeros(length(ys),length(xs),length(zs)) for i in ts]
                    matrices= [DefaultDict{Tuple{UInt16,UInt16,UInt16},UInt32}(0) for i in ts]
                else
                    #matrices=[zeros(length(ys),length(xs)) for i in ts]
                    matrices= [DefaultDict{Tuple{UInt16,UInt16},UInt32}(0) for i in ts]
                end

            end
            return matrices
        end
        fx=:nothing
        fy=:nothing
        fz=:nothing
        if x!=:t
            fx=functions_from_expressions(sistema,[x])[1]
        end
        if y!=:t
            fy=functions_from_expressions(sistema,[y])[1]
        end
        if z!=:nothing
            fz=functions_from_expressions(sistema,[z])[1]
        end
        f_loop=:nothing
        if onematrix==true
            if x==:t
                f_loop=function(matrices,u,i)
                    indu=scale_to_index(fy(u...),ys)
                    if !isnan(indu)
                        matrices[indu,i]+=1
                    end
                end
            else
                f_loop=function(matrices,u,i)
                    indu=scale_to_index(fx(u...),xs)
                    if !isnan(indu)
                        matrices[i,indu]+=1
                    end
                end
            end

        else
            f_loop=function(matrices,u,i)
            iy=scale_to_index(fy(u...),ys)
            ix=scale_to_index(fx(u...),xs)

            if fz!=:nothing
                iz=scale_to_index(fz(u...),zs)
                if !isnan(iy) && !isnan(ix) && !isnan(iz)

                    matrices[i][iy,ix,iz]+=1
                end
            else
                if !isnan(iy) && !isnan(ix)

                    matrices[i][iy,ix]+=1
                end
            end
            end
        end
        if onematrix!=true

            datachannel = RemoteChannel(()->Channel{Any}(Inf), 1)
            if zs!=:nothing
                resdat=[zeros(length(ys),length(xs),length(zs)) for i in 1:length(ts)]
            else
                resdat=[zeros(length(ys),length(xs)) for i in 1:length(ts)]
            end
            flush_task=@async while true
                x= take!(datachannel)
                if x==nothing
                    close(datachannel)
                    break
                else
                    ra=rand(1:100)
                    for ri in length(resdat):-1:1
                        r=resdat[ri]
                        m=pop!(x)
                        for (k,v) in m
                            r[k...]+=v/N
                        end
                        yield()
                    end
                end
            end
       end
       #     addtorsfname=Symbol(string("res_com",rand(1:10000000000)))
       #     resdat=Symbol(string("res_",rand(1:10000000000)))

       #     Base.eval(Main,
       #     quote
       #         if $(zs!=:nothing)
       #             $resdat=[zeros($(length(ys)),$(length(xs)),$(length(zs))) for i in 1:$(length(ts))]
       #         else
       #             $resdat=[zeros($(length(ys)),$(length(xs))) for i in 1:$(length(ts))]
       #         end
       #         $addtorsfname=function(x)
       #             @show "flushing"
       #             for (r,m) in zip($resdat,x)
       #                 for (k,v) in m
       #                     r[k...]+=v/$N
       #                 end
       #             end
       #         end
       #     end
       #     )
       # end
        function finalizacion(x)
            if onematrix==true
                return x/N
            end
            #if x!==nothing
            #    add(x,nothing)
            #end
            #res= Base.eval(Main,:($resdat))
            #Base.eval(Main,quote
            ##    $resdat=nothing
            #    $addtorsfname=nothing
            #end)
            put!(datachannel,nothing)
            wait(flush_task)
            return resdat
        end
        function add(a,b)
            if onematrix==true
                a.+=b
                return a
            end


            for x in [a,b]
                if x===nothing
                    continue
                end
                put!(datachannel,x)
            end
            return nothing
        end
        return f_inicial,f_loop,add,finalizacion
    end
    function calculate_distribution(sistema::ClassicalSystems.ClassicalSystem;distribucion,x,xs,N,y=:nothing,z=:nothing,ts=:nothing,ys=:nothing,zs=:nothing,kargs...)
        ots=ts
        if z==:t
            error("No puedes usar z para el tiempo")
        end
        if x==:t || y==:t
            if ts!=:nothing
                error("No usar el parametro ts al graficar el tiempo en un eje, usar el parametro del eje correspondiente")
            end
            if x==:t
                ts=xs
            else
                ts=ys
            end
        end
        return monte_carlo_over_distribution(sistema,functions_for_calculating_distributions(sistema;x=x,y=y,z=z,zs=zs,ts=ots,xs=xs,ys=ys,N=N)...;ts=ts,N=N,funcion_muestreo=distribucion.sample,kargs...)
    end
    function average(sistema::ClassicalSystems.ClassicalSystem;distribucion,expr,N=1000,ts=[0],kargs...)
        return monte_carlo_over_distribution(sistema,
            functions_for_averaging(sistema;expr=expr,ts=ts,N=N)...;ts=ts,N=N,funcion_muestreo=distribucion.sample,kargs...)
    end
    struct PhaseSpaceDistribution
        probability_density
        sample
    end

    function WignerHW(;q₀,p₀,ħ=:nothing,j=:nothing)
        if j!=:nothing
            ħ=1/j
        end
        σ=sqrt(ħ/2)
        gaussiana_qp=Distributions.MvNormal([q₀,p₀], σ)
        function W(q,p)
                return Distributions.pdf(gaussiana_qp,[q,p])
        end
        function sample()
            q,p=Distributions.rand(gaussiana_qp)
            return [q,p]
        end
        return PhaseSpaceDistribution(W,sample)
    end
    function WignerSU2(;Q₀,P₀,j)
        σ=1/sqrt(2*j)
        gaussiana_Θ=Distributions.Normal(0, σ)
        Rayleigh_Θ=Distributions.Rayleigh(σ)
        R(θ,ϕ,θ0,ϕ0)= #rotacion que manda el polo norte a θ0,ϕ0 aplicada a θ,ϕ
            [atan(sqrt((cos(θ0)*cos(ϕ)*sin(θ)-cos(θ)*sin(θ0))^2+sin(θ)^2*sin(ϕ)^2),cos(θ)*cos(θ0)+cos(ϕ)*sin(θ)*sin(θ0)),
            mod(atan(-cos(ϕ0)*sin(θ)*sin(ϕ)-cos(θ0)*cos(ϕ)*sin(θ)*sin(ϕ0)+cos(θ)*sin(θ0)*sin(ϕ0),-cos(θ0)*cos(ϕ)*cos(ϕ0)*sin(θ)+cos(θ)*cos(ϕ0)*sin(θ0)+sin(θ)*sin(ϕ)*sin(ϕ0)),2*pi)]


        centroθ=PhaseSpaces.θ_of_QP(Q₀,P₀)
        centroϕ=PhaseSpaces.ϕ_of_QP(Q₀,P₀)
        function Wθ(θ,ϕ)
            try
                Θ=acos(cos(θ)*cos(centroθ)+ sin(θ)*sin(centroθ)*cos(ϕ - centroϕ))
                return Distributions.pdf(gaussiana_Θ,Θ)*sqrt(2*π)*σ *(1/(2*π*σ^2))
            catch
                return 0
            end
        end
        function W(Q,P)
            try

                 θ=PhaseSpaces.θ_of_QP(Q,P)
                 ϕ=PhaseSpaces.ϕ_of_QP(Q,P)
                 return Wθ(θ,ϕ)
            catch
                return 0
            end
        end


        function sample()
            Θ=Distributions.rand(Rayleigh_Θ)
            θ,ϕ = R(Θ,rand()*2*pi,centroθ,centroϕ)
            return [PhaseSpaces.Q_of_θϕ(θ,ϕ),PhaseSpaces.P_of_θϕ(θ,ϕ)]
        end
        return PhaseSpaceDistribution(W,sample)
    end
    function WignerHWxSU2(;Q₀,q₀,P₀,p₀,j)
        WSU2=WignerSU2(Q₀=Q₀,P₀=P₀,j=j)
        WHW=WignerHW(q₀=q₀,p₀=p₀,j=j)
        W(Q,q,P,p)=WSU2.probability_density(Q,P)*WHW.probability_density(q,p)
        function sample()
            Q,P=WSU2.sample()
            q,p=WHW.sample()
            return [Q,q,P,p]
        end
        return PhaseSpaceDistribution(W,sample)
    end
    function WignerHWxSU2(u₀,j)

        W=WignerHWxSU2(Q₀=u₀[1],q₀=u₀[2],P₀=u₀[3],p₀=u₀[4],j=j)
        ev(u)=W.probability_density(u...)
        return PhaseSpaceDistribution(ev,W.sample)
    end
    function WignerSU2(u₀,j)

        W=WignerSU2(Q₀=u₀[1],P₀=u₀[2],j=j)
        ev(u)=W.probability_density(u...)
        return PhaseSpaceDistribution(ev,W.sample)
    end
     function WignerHW(u₀;kargs...)

        W=WignerHW(q₀=u₀[1],p₀=u₀[2];kargs...)
        ev(u)=W.probability_density(u...)
        return PhaseSpaceDistribution(ev,W.sample)
    end
    function functions_for_survival_probability(sistema::ClassicalSystems.ClassicalSystem;ts,N,funcion_distribucion,ħ)
        f_inicial= function()
            zeros(length(ts))
        end
        f_loop=function(valores,u,i)
                valores[i]+=funcion_distribucion(u)
        end
        f_final=function(res)
            L=length(sistema.varnames)/2
            return res*(2*π*ħ)^L/N
        end

        return f_inicial,f_loop,+,f_final
    end
    function survival_probability(sistema::ClassicalSystems.ClassicalSystem;ħ,distribucion,ts,N,kargs...)
        return TruncatedWignerApproximation.monte_carlo_over_distribution(sistema,
            functions_for_survival_probability(sistema;ts=ts,N=N,funcion_distribucion=distribucion.probability_density,ħ=ħ)...;ts=ts,N=N,funcion_muestreo=distribucion.sample,integate_backwards=true,kargs...)
    end
    module Weyl

        using SymEngine
        function n(j)
            SymEngine.@vars q p
            return (j*(p^2 + q^2) - 1)/2 #Weyl(n̂)=(p²+q²- ħ)/2 y nuestras q,p están escaladas por sqrt(j)
        end
        function n²(j)
            return n(j)^2 - 1/4 #Weyl(n̂²)=Weyl(n̂)² - ħ²/4 (ver notas de polkovnikov para la boulder summer school)
        end


        #"The Moyal Representation for Spin" JOSEPH C. VÁRILLY AND JOSÉ M. GRACIA-BONDÍA Ann. of Phys. 190, 107-148 (1989)
        #https://doi.org/10.1016/0003-4916(89)90262-5
        macro SU2angles(ex)
            q= quote
                SymEngine.@vars Q P
                cosθ=(P^2+Q^2)/2-1
                sinθ=sqrt(1-cosθ^2)
                cosϕ=Q/sqrt(P^2 + Q^2)
                sinϕ=-P/sqrt(P^2 + Q^2)
                $ex
            end
            return esc(q)
        end
        function Jz(j)
            return @SU2angles sqrt(j*(j+1))*cosθ
        end
        function Jx(j)
            return @SU2angles sqrt(j*(j+1))*sinθ*cosϕ
        end
        function Jy(j)
            return @SU2angles sqrt(j*(j+1))*sinθ*sinϕ
        end
        function Jz²(j)
            bj=sqrt(j*(j+1)*(2*j-1)*(2*j+3))
            return @SU2angles (j*(j+1))/3 + (bj/2)*(cosθ^2 - 1/3)
        end
        function Jx²(j)
            bj=sqrt(j*(j+1)*(2*j-1)*(2*j+3))
            return @SU2angles (j*(j+1))/3 + (bj/2)*((sinθ*cosϕ)^2 - 1/3)
        end
        function Jy²(j)
            bj=sqrt(j*(j+1)*(2*j-1)*(2*j+3))
            return @SU2angles (j*(j+1))/3 + (bj/2)*((sinθ*sinϕ)^2 - 1/3)
        end
    end

end
