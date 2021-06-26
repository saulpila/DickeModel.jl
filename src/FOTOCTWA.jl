module FOTOCTWA

export calculate_log_slope,functions_for_variance,variance

    import ..TruncatedWignerApproximation
    import ..ClassicalSystems
    using Statistics
    using LsqFit
    function calculate_log_slope(sistema;distribucion,vars,tf,sattime,deltat,N,returnall=false,kargs...)
        tsat=tf-sattime
        ts=0:deltat:tf
        varianzas= variance(sistema;distribucion=distribucion,ts=ts,vars=vars,N=N,kargs...)
        varprom=[mean(v) for v in varianzas]
        itsat=Int(round(tsat/deltat))

        valsat=mean(varprom[itsat:end])
        logsat=Float64(log(valsat))
        firstmorethansat=findfirst(varprom .> valsat)-1
        if minimum(varprom[firstmorethansat:end])<=varprom[1]
            if !returnall
                return 0
            else
                return 0,0,logsat,varprom,ts
            end
        end
        valsatoscilations=log(maximum(varprom[tf-sattime:end])/minimum(varprom[tf-sattime:end]))
        top=log(valsat) - valsatoscilations/2
        lastacceptable=findfirst(varprom .> exp(top))-1
        @. model(x,p)=p[1]*x+p[2]
        p0 = [0.5,varprom[1]]
        fit= curve_fit(model, ts[1:lastacceptable], Float64.(log.(varprom[1:lastacceptable])), p0)
        m0= fit.param[1]
        if !returnall
            return m0
        else
            return m0,fit.param[2],logsat,varprom,ts
        end
    end

    function functions_for_variance(sistema::ClassicalSystems.ClassicalSystem;expr,ts=[0.0],N,con_promedio=false)
        onevar=false

        if !isa(expr,Array)
            expr=[expr]
            onevar=true
        end
        tam=length(expr)

        apromediar=[expr;[(args...)->f(args...)^2 for f in TruncatedWignerApproximation.vars_from_expression(sistema,expr)]]

        f_inicial,f_loop,_,f_final=TruncatedWignerApproximation.functions_for_averaging(sistema,expr=apromediar,ts=ts,N=N)
        mf_final=function(res)
            r=f_final(res)
            vars=[[re[i+tam] - re[i]^2 for i in 1:tam] for re in r]
            if onevar
                vars=[v[1] for v in vars]
            end
            if con_promedio
                proms=[[re[i] for i in 1:tam] for re in r]
                if onevar
                    proms=[p[1] for p in proms]
                end
                vars=vars,proms
            end

            return vars
        end
        return f_inicial,f_loop,+,mf_final
    end
    function variance(sistema::ClassicalSystems.ClassicalSystem;distribucion,expr,ts=[0.0],N=1000,con_promedio=false,kargs...)
        return TruncatedWignerApproximation.monte_carlo_over_distribution(sistema,
            functions_for_variance(sistema;expr=expr,ts=ts,N=N,con_promedio=con_promedio)...;ts=ts,N=N,funcion_muestreo=distribucion.sample,kargs...)
    end

end
