module DickeBCE

export H_BCE,QuantumSystem,coherentBCE,coherentBCE!,Wigner,WignerProjQP,WignerProjqp
    import ClassicalSystems
    using LinearAlgebra
    using ProgressMeter
    using SparseArrays
    using Distributed
    import PhaseSpaces
    using GSL
    using WignerSymbols
    using Distributed
    using CSV
    using Dicke
    using Distributions
    import Random
    using Statistics
    using DataFrames,Tables
    #using Memoize
    struct QuantumSystem
       params
    end
    function QuantumSystem(sistema::ClassicalSystems.ClassicalSystem;j,Nmax=nothing)
    #si Nmax es nothing, se busca la diagonalizacion en disco con el Nmax mas grande
        QuantumSystem([sistema.params;[j,Nmax]])
    end
    ind(n,m,N,j)= n*(N+1)+Int(m+j) +1 #el +1 es por julia que cuenta desde 1
    function logfact(n::Int)::Float64
        if n==0
            return 0.0
        end
        return sum(log(i) for i in 1:n)
    end
    function _sqrt_factorial(n::Int)
        if n==0
            return 1
        end
        prod(sqrt(i) for i in 1:n)
    end
    function _sqrt_binomial(n::Int,k::Int)
        return exp(_sqrt_binomiallog(n,k))
    end
    function _sqrt_binomiallog(n::Int,k::Int)
        return (logfact(n)-logfact(k)-logfact(n-k))/2
    end
    function get_params(sistema::QuantumSystem)
        ω₀,ω,γ,j,Nmax=sistema.params
        if Nmax!==nothing
            Nmax=Int(Nmax)
        end
        return ω₀,ω,γ,j,Nmax
    end
    function Π_BCE(sistema::QuantumSystem)
        ω₀,ω,γ,j,Nmax=get_params(sistema)
        N=Int(2*j)
        Π=spzeros(Float64,Nmax*(N+1), Nmax*(N+1))
        for indexN in 0:(Nmax-1), m in -j:j
            Π[ind(indexN,m,N,j),ind(indexN,-m,N,j)]=(-1)^indexN
        end
        return Π
    end
    function eigenstate_parities(sistema::QuantumSystem,eigenstates)
        Π=Π_BCE(sistema)
        ΠE=Π*eigenstates
        return Int8[Integer(round(dot(view(ΠE,:,k),view(eigenstates,:,k)))) for k in 1:size(eigenstates)[2]]
    end
    function PopulationOfLeveln(sistemaQ,state,n)
        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
        N=Int(2*j)
        return sum(abs(state[ind(n,m,N,j)]) for m in -j:j)
    end
    function eigenenergies(sistemaQ::QuantumSystem;args...)
        return diagonalization(sistemaQ;only_eigenenergies=true,args...)
    end
    function diagonalization(sistemaQ::QuantumSystem;load_cache=true,save_cache=true,cache_folder=joinpath(homedir(),"dicke_diagonalizations"),maxε::Real=5.0,onlyload::Union{AbstractVector{<:Integer},Nothing}=nothing,only_eigenenergies=false)
        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
        filename(ω₀,ω,γ,j,Nmax)="dicke_bce_ω₀=$(ω₀),ω=$(ω),γ=$(γ),j=$(j),Nmax=$(Nmax)"
        filelocation(ω₀,ω,γ,j,Nmax)=joinpath(cache_folder,filename(ω₀,ω,γ,j,Nmax))
        if Nmax===nothing
            if save_cache==false
                error("You have to pass Nmax if save_cache==false")
            end
            reg=Regex(string(filename(ω₀,ω,γ,j,"(\\d+)"),"_E.csv"))
            diags=readdir(cache_folder)
            for d in diags
                m=match(reg,d)
                if m===nothing
                    continue
                end
                poss_nmax=parse(Int64,m.captures[1])
                flname=filelocation(ω₀,ω,γ,j,poss_nmax)
                if isfile("$(flname)_E.csv") && isfile("$(flname)_V.csv")
                    if Nmax===nothing
                        Nmax=poss_nmax
                    else
                        Nmax=max(poss_nmax,Nmax)
                    end
                end
            end
            if Nmax===nothing
                error("No diagonalizations found, please pass a value for Nmax to diagonalize")
            else
                sistemaQ.params[5]=Nmax
            end
        end

        flname=filelocation(ω₀,ω,γ,j,Nmax)
        if load_cache

             if isfile("$(flname)_E.csv") && isfile("$(flname)_V.csv")
                display("Loading diagonalization: Nmax=$Nmax")
                eigenenergies=DataFrame(CSV.File("$(flname)_E.csv")).E[:]
                if only_eigenenergies
                    return eigenenergies
                end
                eigenstates=Tables.matrix(DataFrame(CSV.File("$(flname)_V.csv",select=onlyload,type=Float64,threaded =true,tasks=4)))
                return eigenenergies,eigenstates
             end
        end


        Hq=H_BCE(sistemaQ)


        display("Diagonalizing")
        FullH=Symmetric(Matrix(Hq));
        Emin=-10*j
        Emax=maxε*j #maxenergy
        eigenenergies,eigenstates=eigen(FullH,Real(Emin),Real(Emax));

        convstates=(length(eigenenergies)-1)
        for k in 1:(length(eigenenergies)-1)
            if PopulationOfLeveln(sistemaQ,eigenstates[:,k],Nmax-1) >1e-3
                convstates=k
                break
            end
        end
        eigenstates=eigenstates[:,1:convstates]
        eigenenergies=eigenenergies[1:convstates]
        display(" $convstates converged states were obtained up to ε=$(eigenenergies[end]/j)")

        display("Fixing numerical degeneracies (correcting parity)")
        #corregimos por degeneraciones numéricas
        Π=Π_BCE(sistemaQ)
        for k in 1:(length(eigenenergies)-1)
            if(eigenenergies[k+1]-eigenenergies[k]<1e-5)
                _correct_eigenstates(eigenstates,k,k+1,Π)
            end
        end
        if save_cache
            mkpath(cache_folder)
            CSV.write("$(flname)_V.csv",DataFrame(eigenstates))
            CSV.write("$(flname)_E.csv",DataFrame(E=eigenenergies))
            #esto es para que dropbox no sincronice estos archivos que pueden pesar un monton
            #if Sys.iswindows()
            #    for fl in ["$(flname)_V.csv","$(flname)_E.csv"]
            #        absfile=Base.Filesystem.abspath(fl)
            #        cmd="& {Set-Content -Path \"$absfile\"  -Stream com.dropbox.ignored -Value 1}"
            #        run(`powershell -Command $cmd`)
            #    end
            #end
        end
        return eigenenergies,eigenstates
    end
    function _correct_eigenstates(eigenstates,k1,k2,Π)
        v1=eigenstates[:,k1]
        v2=eigenstates[:,k2]
        eigenvals,mat =eigen(Transpose([v1 v2])*Π*[v1 v2])
        if !issetequal([1.0,-1.0],round.(eigenvals,digits=7))
            return
        end
        λ_1=mat[1,1]
        λ_2=mat[2,1]
        if λ_1<0
            λ_1,λ_2=-λ_1,-λ_2
        end
        μ_1=mat[1,2]
        μ_2=mat[2,2]
        if μ_2<0
            μ_1,μ_2=-μ_1,-μ_2
        end
        eigenstates[:,k1].= λ_1*v1+λ_2*v2;
        eigenstates[:,k2].=μ_1*v1+μ_2*v2;
        return true
    end
    function H_BCE(sistema::QuantumSystem;_ignorediagonalterms=false)
        ω₀,ω,γ,j,Nmax=get_params(sistema)
        N=Int(2*j)
        G=2*γ/(ω*sqrt(N))
        function DNNp(N,Np)

            f(k)=exp((logfact(N)+logfact(Np))/2 - (logfact(Np-k)+logfact(N-k)+logfact(k)))
            return exp(-G^2/2)*sum((f(k)*(-1)^k) * (G^(N+Np-2*k)) for k in 0:Np)
        end
        H=spzeros(Float64,Nmax*(N+1), Nmax*(N+1))
        p = Progress(Int(Nmax*(Nmax-1)/2))
        for indexN in 0:(Nmax-1)
            for indexNprime in 0:indexN
                D=DNNp(indexN,indexNprime)
                next!(p)
                for m in -j:j
                    if !_ignorediagonalterms
                        if indexNprime==indexN
                            H[ind(indexN,m,N,j),ind(indexN,m,N,j)]=ω*(indexN-G^2*m^2)
                        end
                    end
                    if m<j
                        o=(-1)^indexN*D
                        v=-ω₀*sqrt(j*(j+1)-m*(m+1))*o/2 #queondaconeseentre2
                        H[ind(indexNprime,m+1,N,j),ind(indexN,m,N,j)]=v

                    end
                    if m>-j
                        o=(-1)^indexNprime*D
                        v=-ω₀*sqrt(j*(j+1)-m*(m-1))*o/2 #queondaconeseentre2
                       H[ind(indexNprime,m-1,N,j),ind(indexN,m,N,j)]=v
                    end
                end
            end
        end
        H=Symmetric(H);
        finish!(p)
        return H
    end
    function Jz(sistema::QuantumSystem)
        ω₀,ω,γ,j,Nmax=get_params(sistema)
        return H_BCE(QuantumSystem(Dicke.ClassicalSystem(ω₀=1,ω=1,γ=1),j=j,Nmax=Nmax),_ignorediagonalterms=true)
    end
    function Jx(sistema::QuantumSystem)
        ω₀,ω,γ,j,Nmax=get_params(sistema)
        N=Int(2*j)
        Jx=spzeros(Float64,Nmax*(N+1), Nmax*(N+1))
        for indexN in 0:(Nmax-1)
            for m in -j:j
               Jx[ind(indexN,m,N,j),ind(indexN,m,N,j)]=m
            end
        end
        return Diagonal(Jx)
    end

    function ParticipationRatio(state;eigenstates,eigenenergies=nothing,count_degeneracies=eigenenergies!=nothing,degentol=1e-5)
        eigenbasis_state=transpose(eigenstates)*state
        if count_degeneracies
            IPR=0.0
            degsbspace=0.0
            for i in 1:length(eigenenergies)

                if  i!=1 && eigenenergies[i]-eigenenergies[i-1]>degentol
                    IPR+=degsbspace^2
                    degsbspace=0.0
                end
                degsbspace+=abs2(eigenbasis_state[i])
                if i==length(eigenenergies)
                    IPR+=degsbspace^2
                end
            end
            return 1/IPR
        else
            return 1/sum(abs2(c)^2 for c in eigenbasis_state)
        end
    end
    function RFactor(sistemaC::ClassicalSystems.ClassicalSystem,sistemaQ::QuantumSystem,x,state;eigenstates,eigenenergies,exact=false)
        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
        ε=Dicke.hamiltonian(sistemaC)(x)
        σ=Dicke.energy_width_of_coherent_state(sistemaC,x,j)
        ν(ε)=Dicke.density_of_states(sistemaC,j=j,ε=ε)
        n=Distributions.Normal(ε,σ)
        ρ(ε)=Distributions.pdf(n,ε)

        if exact
            ρoverν(ε)=ρ(ε)/ν(ε)
            f=mean(ρoverν(Distributions.rand(n)) for i in 1:10000)/j #monte carlo integrate ρ^2/ν sampleando de ρ
        else
            f=1/(2*sqrt(pi)*σ*j*ν(ε))
        end
        return 2*ParticipationRatio(state,eigenstates=eigenstates,eigenenergies=eigenenergies,count_degeneracies=true)*f #cuidado esto solo funciona sin degeneraciones, hay que contar sobre subespacios degenerados...
    end
    function SurvivalProbability(t::Union{Real,AbstractArray{<:Real}};state,eigenstates,eigenenergies,normwarning=0.99)
        if typeof(t)<:Real
            return SurvivalProbability([t],state=state,eigenstates=eigenstates,eigenenergies=eigenenergies)[1]
        end
        eigenbasis_state=transpose(eigenstates)*state
        if norm(eigenbasis_state)<normwarning
            @warn "The coherent state is not converged, make Nmax bigger"
        end
        return [abs2(sum(exp(-im*ti*Ek)*abs2(ck) for (ck,Ek) in zip(eigenbasis_state,eigenenergies))) for ti in t]
     end
     function evolve(t::Union{Real,AbstractArray{<:Real}};state,eigenstates,eigenenergies,normwarning=0.99)
        if typeof(t)<:Real
            return evolve([t],state=state,eigenstates=eigenstates,eigenenergies=eigenenergies)
        end
        eigenbasis_state=transpose(eigenstates)*state
        if norm(eigenbasis_state)<normwarning
            @warn "The coherent state is not converged, make Nmax bigger"
        end
        eigenstates*hcat((Diagonal(exp.(-im*ti*eigenenergies))*eigenbasis_state for ti in t)...)
    end

    function coherentBCE(sistema::QuantumSystem,punto;normwarning=0.99,extra_phase::Complex=1.0+0im,tol::Real=0.0)
        ω₀,ω,γ,j,Nmax=get_params(sistema)

        data = zeros(Complex{Float64},Nmax*(Int(2*j)+1))
        coherentBCE!(sistema,punto,data;normwarning=normwarning,extra_phase=extra_phase,tol=tol)
        return data
    end
    coherentBCE!(sistema::QuantumSystem,punto,data;normwarning=0.99,add_state_to_data=false,extra_phase::Complex=1.0+0im,tol=1e-5)=
        _coherentOverlapOrCoherentState(sistema,punto;
                        datacache=data,tol=tol,normwarning=normwarning,add_state_to_data=add_state_to_data,extra_phase=extra_phase)


    function _Cohpart1(n,m,alpha,G)
        s=alpha+G*m
        if n==0 && s==0.0
            return 0.0im #para evitar 0*log(0)
        end
        return n*log(s) - logfact(n)/2
    end
    function _Cohpart2(m,N,j,w,alpha,G)
        a= _sqrt_binomiallog(N,Int(j + m)) + log(w)*m -abs2(G*m)/2 -G*m*(alpha)
        return a
    end

    function qtlmin(dist,tol,mn)
        vl=mn
        try
            dt=Distributions.quantile(dist,tol)
            vl=max(vl,dt)
        catch
        end
        return vl
    end
    function qtlmax(dist,tol,mx)
        vl=mx
        try
            dt=Distributions.quantile(dist,tol)
            vl=min(vl,dt)
        catch
        end
        return vl
    end

    function HusimiOfCoherent(sistemaQ::QuantumSystem,x::AbstractArray{<:Real,1},y::AbstractArray{<:Real,1})
        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
        Q1,q1,P1,p1=x
        Q2,q2,P2,p2=y
        return exp(-(j/2)*((q1-q2)^2 + (p1-p2)^2))*cos(PhaseSpaces.arc_between_QP(Q1,P1,Q2,P2)/2)^(4j)
    end
    onelistcache=[0.0im]
    function coherentOverlap(sistemaQ::QuantumSystem,punto::AbstractArray{<:Real,1},estado::AbstractArray{<:Number,1};tol=1e-5,normwarning=0.99,datacache=nothing)
        if datacache===nothing || !isa(datacache, AbstractArray)
            datacache=onelistcache
        end
        return _coherentOverlap(sistemaQ,punto,estado;datalength=1,datacache=datacache,tol=tol,normwarning=normwarning)[1]
    end
    function coherentOverlap(sistemaQ::QuantumSystem,punto::AbstractArray{<:Real,1},estados::AbstractArray{<:Number,2};tol=1e-5,normwarning=0.99,datacache::Union{AbstractArray{Complex{Float64},1},Nothing}=nothing)
        d=size(estados)[2]
        if datacache===nothing
            datacache=zeros(Complex{Float64},d)
        end
        return _coherentOverlap(sistemaQ,punto,estados;datacache=datacache,datalength=d,tol=tol,normwarning=normwarning)
    end
#si das estados hace el overlap, si no construye un coherente
#si add_state_to_data
function _coherentOverlapOrCoherentState(sistemaQ::QuantumSystem,punto::AbstractArray{<:Real,1};
                        datacache,datalength=nothing,tol=1e-5,normwarning=0.99,
                        estados::Union{AbstractArray{<:Number,1},AbstractArray{<:Number,2},Nothing}=nothing,add_state_to_data=false,extra_phase::Complex=1.0+0im)
        #@fastmath begin
            if !add_state_to_data
                datacache.=0.0im
            end

            ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
            N=Int(2*j)
            G=2*γ/(ω*sqrt(N))
            Q,q,P,p=punto
            theta=PhaseSpaces.θ_of_QP(Q,P)
            phi=PhaseSpaces.ϕ_of_QP(Q,P)

            z=exp(-phi*im)*tan(theta/2)
            alpha=(q+im*p)*sqrt(j/2)
            m=-j


            nrm=0.0
            w=(1+z)/(1-z)
            part0=-abs2(alpha)/2 +j*log(w/(1+abs2(w)))
            #  @show w,z


            if tol==0
               mmin=-j
               mmax=j
            else
                atomicDist=Distributions.Binomial(N,abs2(w)/(1+abs2(w)))
                mmin=qtlmin(atomicDist,tol/4,0)-j
                mmax=qtlmax(atomicDist,1-tol/4,N)-j
            end
            for m in mmin:mmax
            #for m in -j:j
                part2m=part0+_Cohpart2(m,N,j,w,alpha,G)
                if tol==0
                    nmin=0
                    nmax=Nmax-1
                else
                    efftol=tol/((mmax-mmin+1) *Distributions.pdf(atomicDist,m+j))
                    d=Distributions.Poisson(abs(G*m+alpha)^2)
                    nmin= qtlmin(d,efftol/4,0)
                    nmax=qtlmax(d,1-efftol/4,Nmax-1)
                end
                for n in nmin:nmax
                    v=conj(exp(_Cohpart1(n,m,alpha,G)+part2m))
                    if estados!==nothing
                        for e in 1:datalength
                            datacache[e]+=estados[ind(n,m,N,j),e]*v
                        end

                    else
                       datacache[ind(n,m,N,j)]+=conj(v)*extra_phase
                    end
                    nrm+=abs2(v)
                end
            end
            if nrm<normwarning
                @warn "The coherent state is not converged, make Nmax bigger or tol smaller"
            end
        #end
        return datacache
    end
    _coherentOverlap(sistemaQ::QuantumSystem,
    punto::AbstractArray{<:Real,1},
    estados::Union{AbstractArray{<:Number,1},AbstractArray{<:Number,2}};
    datacache,
    datalength,
    tol=1e-5,
    normwarning=0.99) =_coherentOverlapOrCoherentState(sistemaQ,punto;
                                                        datacache=datacache,
                                                        datalength=datalength,
                                                        tol=tol,
                                                        estados=estados,
                                                        normwarning=normwarning)
    function Husimi(args...;kargs...)
        r=coherentOverlap(args...;kargs...)
        return abs2.(r)
    end

    #esto no funciona para j grande, me cambie a WignerSymbols
    function clebsch(j1, j2, j3, m1, m2, m3)

        if m3 != m1 + m2
            return 0
        end
        vmin = Int(max(-j1 + j2 + m3, -j1 + m1, 0))
        vmax = Int(min(j2 + j3 + m1, j3 - j1 + j2, j3 + m3))

        C = exp((log(2.0 * j3 + 1.0) + logfact(j3 + j1 - j2) +
                    logfact(j3 - j1 + j2) + logfact(j1 + j2 - j3) +
                    logfact(j3 + m3) + logfact(j3 - m3) -
                    (logfact(j1 + j2 + j3 + 1) +
                    logfact(j1 - m1) + logfact(j1 + m1) +
                    logfact(j2 - m2) + logfact(j2 + m2)))/2)
        S = 0
        for v in vmin:vmax
            S += (-1.0) ^ (v + j2 + m2)* exp( -logfact(v) +
                logfact(j2 + j3 + m1 - v) + logfact(j1 - m1 + v) -
                logfact(j3 - j1 + j2 - v) - logfact(j3 + m3 - v) -
                logfact(v + j1 - j2 - m3))
        end
        C = C * S
        return C
    end
    function sph_harm(l,m,theta,phi;legendre_sphPlm_cache=:nothing)
        fact=1
        if m<0
            fact=(-1)^(-m)
        end
        d=0
        if legendre_sphPlm_cache==:nothing
            d=GSL.sf_legendre_sphPlm(l,abs(m),cos(theta))
        else
            d=legendre_sphPlm_cache[l,abs(m),theta]
        end
        return fact*d*exp(im*m*phi)
    end
    function generate_hermite_cache(n_max,pts)
        cache=Dict{NTuple{2,Real},Complex}()
        points=unique(pts)
        for x in points

            data=GSL.sf_hermite_phys_array(n_max,x)
            for (n,d) in enumerate(data)
                cache[n-1,x]=d
            end
        end
        return cache
    end
    function generate_legendre_sphPlm_cache(l_max,pts)
        cache=Dict{NTuple{3,Real},Complex}()
        thetas=unique(x[3] for x in pts)
        for theta in thetas

            data=GSL.sf_legendre_array_e(GSL.GSL_SF_LEGENDRE_SPHARM,l_max,cos(theta),-1)
            for l in 0:l_max,m in 0:l
                cache[l,m,theta]=data[GSL.sf_legendre_array_index(l, m)+1]
            end
        end
        return cache
    end
    function laguerre(n,a,x)

        if a<0
            return x^(-a)/GSL.sf_poch(-n,-a)*laguerre(n+a,-a,x)
        end
        return GSL.sf_laguerre_n(n,a,x)
    end
    function RotZtoX(θ,ϕ)
        a=-π/2
        z=[cos(θ/2) -exp(-im*ϕ)*sin(θ/2)]*[cos(a/2) -sin(a/2);
                                        sin(a/2) cos(a/2)]
        return 2*atan(abs(z[2]),abs(z[1])),angle(z[2])-angle(z[1])
    end
    function classical_QqPp_to_quantum_rot_qpθϕ(j,QqPp)
        Q,q,P,p=QqPp
        q=sqrt(j)*q
        p=sqrt(j)*p
        θ,ϕ=RotZtoX(PhaseSpaces.θ_of_QP(Q,P),PhaseSpaces.ϕ_of_QP(Q,P))
        return [q,p,θ,ϕ]
    end

    function sph_harm_pt(k,i,pt;legendre_sphPlm_cache=:nothing)
        return sph_harm(k, i, pt[3],pt[4];legendre_sphPlm_cache=legendre_sphPlm_cache)
    end

    function Wigner_Fock(n,np,q,p)
        two_r_squared=(q^2+p^2)*2
        if two_r_squared==0
            two_r_squared=1e-30
        end
        return (exp((logfact(np)-logfact(n))/2 +
            im*(np-n)*atan(p,q)- two_r_squared/2)*(-1)^(np)/π * two_r_squared^((n-np)/2)*laguerre(np,n-np,two_r_squared))
    end
    function Wigner_Fock_Recorrida(n,np,m,mp,G,pt)
        q=pt[1]
        p=pt[2]

        qrec=q+G*(m+mp)/2 #por el desplazamiento de alpha_m
        return exp(-im*p*G*(mp-m))*Wigner_Fock(n,np,qrec,p)
    end
    function Wigner(sistemaQ,state,puntos)
        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)

        N=Int(2*j)
        G=2*γ/(ω*sqrt(N))

        puntos=classical_QqPp_to_quantum_rot_qpθϕ.(j,puntos)
        legendre_sphPlm_cache=generate_legendre_sphPlm_cache(Int(2*j),puntos)
        #θ,ϕ=theta,phi


        function B(n,l,np,lp)
            altsm(n-l)*altsm(np-lp)*exp((logfact(n)+logfact(np))/2 -
                (logfact(l)+logfact(lp)))
        end





        W=fill!(similar(puntos,Complex{Float64}),0)
        _Wd=fill!(similar(puntos,Complex{Float64}),0)
        for m in -j:j, mp in -j:j

            mWd=Wigner_Dicke!(m,mp,j,puntos,_Wd,legendre_sphPlm_cache=legendre_sphPlm_cache)
            for n in 0:(Nmax-1),np in 0:n
                _factor=state[ind(n,m,N,j)]*conj(state[ind(np,mp,N,j)])
                if norm(_factor)<1e-6
                    continue
                end
                mW=_factor*mWd
                mW.*=Wigner_Fock_Recorrida.(n,np,m,mp,G,puntos)
                W+= real.(mW*(n==np ? 1 : 2))
             #   if n!=np
             #       W+=conj.(mW)
              #  end
            end


        end
        return real.(W)
    end
    function Wigner_Dicke!(m::Int64,mp::Int64,j::Int64,puntos::Array{Array{Float64,1},1},result::Array{Complex{Float64},1};legendre_sphPlm_cache=:nothing)

        result.=0.0+0.0im

        for k in abs(m-mp):Int(2 * j)
            i=m-mp
            _f=(-1)^(j - m - i) * clebschgordan(j,m, j,-mp, k,i)
            if _f !=0
                 result +=   _f*sph_harm_pt.(k, i, puntos,legendre_sphPlm_cache=legendre_sphPlm_cache)
            end
        end
        return result*sqrt((2*j+1)/(4*pi)) #normalizamos
    end
    function realmult(w2::Complex{Float64},factor)
        return real(factor*w2)
    end
    function WignerProjqp(sistemaQ,states,puntos)
        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
        N=Int(2*j)
        G=2*γ/(ω*sqrt(N))
        puntos=classical_QqPp_to_quantum_rot_qpθϕ.(j,puntos)
        tmD=0.0
        tmF=0.0
        tmM=0.0
        _Wd=similar(puntos,Complex{Float64})
        AllWs=[similar(puntos,Float64) for s in states]
        fill!.(AllWs,0.0)

        legendre_sphPlm_cache=generate_legendre_sphPlm_cache(Int(2*j),puntos)

        p = Progress(Int((2*j+1)^2*(Nmax-1)*Nmax/2),1)
        for m in -j:j, mp in -j:j

            _Wd=Wigner_Dicke!(m,mp,j,puntos,_Wd,legendre_sphPlm_cache=legendre_sphPlm_cache)

            for n in 0:(Nmax-1),np in 0:n
                if np==n && mp==m
                    lag_fac=1.0
                else
                lag_fac=exp((logfact(np)-logfact(n))/2 - abs(G*(m-mp))^2/2 + (n-np)*log(G*abs(m-mp)))*sign(m-mp)^(n-np)*laguerre(np,n-np,abs(G*(m-mp))^2)

                end
                for i in 1:length(states)
                    _factor=states[i][ind(n,m,N,j)]*conj(states[i][ind(np,mp,N,j)])
                    finalfactor=_factor*lag_fac*(n==np ? 1 : 2)
                    if abs(finalfactor) ==0
                        continue
                    end
                 #   AllWs[i] = AllWs[i]+sumrealpart.(AllWs[i],_Wd,finalfactor) #forma eficiente de hacer W+=real.(finalfactor*_Wd)
                    AllWs[i].+=realmult.(_Wd,finalfactor)
                end
                next!(p)
             #   if n!=np
             #       W+=conj.(mW)
              #  end
            end

        end

     #   @show tmF,tmD,tmM,tmD+tmF+tmM
        finish!(p)
        return AllWs
    end

    function WignerProjp(sistemaQ,states,puntos)
            ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
            N=Int(2*j)
            G=2*γ/(ω*sqrt(N))
            puntos=classical_QqPp_to_quantum_rot_qpθϕ.(j,puntos)
            tmD=0.0
            tmF=0.0
            tmM=0.0
            _Wd=similar(puntos,Complex{Float64})
            _FactorHermite=similar(puntos,Float64)
            AllWs=[similar(puntos,Float64) for s in states]
            fill!.(AllWs,0.0)

            legendre_sphPlm_cache=generate_legendre_sphPlm_cache(Int(2*j),puntos)
            hermite_cache=generate_hermite_cache(Nmax-1,Iterators.flatten((x[1]+G*m for m in -j:j) for x in puntos))
            function Fock_n_np_wavefunction(n::Int64,m::Int64,np::Int64,mp::Int64,x::Array{Float64,1})
                return exp(-(logfact(np)+logfact(n)+(np+n)*log(2))/2 - (x[1]+G*m)^2/2 - (x[1]+G*mp)^2/2)*hermite_cache[n,x[1]+G*m]*hermite_cache[np,x[1]+G*mp]/sqrt(pi)
            end
            p = Progress(Int((2*j+1)^2*(Nmax-1)*Nmax/2),1)
            for m in -j:j, mp in -j:j

                _Wd=Wigner_Dicke!(m,mp,j,puntos,_Wd,legendre_sphPlm_cache=legendre_sphPlm_cache)

                for n in 0:(Nmax-1),np in 0:n

                    _FactorHermite.=Fock_n_np_wavefunction.(n,m,np,mp,puntos)


                    for i in 1:length(states)
                        _factor=states[i][ind(n,m,N,j)]*conj(states[i][ind(np,mp,N,j)])
                        finalfactor=_factor*(n==np ? 1 : 2)
                        if abs(finalfactor) ==0
                            continue
                        end
                        AllWs[i].+=realmult.(_Wd,finalfactor).*_FactorHermite
                    end
                    next!(p)
                 #   if n!=np
                 #       W+=conj.(mW)
                  #  end
                end

            end

         #   @show tmF,tmD,tmM,tmD+tmF+tmM
            finish!(p)
            return AllWs
        end
    function WignerProjQP(sistemaQ,state,puntos)
        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)

        N=Int(2*j)
        G=2*γ/(ω*sqrt(N))
        puntos=classical_QqPp_to_quantum_rot_qpθϕ.(j,puntos)
        #θ,ϕ=theta,phi

        W=fill!(similar(puntos,Complex),0)
        _Wd=fill!(similar(puntos,Complex),0)
        for m in -j:j
          #  @show m
            for n in 0:(Nmax-1),np in 0:n
                _factor=state[ind(n,m,N,j)]*conj(state[ind(np,m,N,j)])
                if norm(_factor)<1e-20
                    continue
                end
                mW=_factor*Wigner_Fock_Recorrida.(n,np,m,m,G,puntos)

                W+= real.(mW*(n==np ? 1 : 2))
             #   if n!=np
             #       W+=conj.(mW)
              #  end
            end


        end
        return real.(W)
    end
    function randomState(sistemaC::ClassicalSystems.ClassicalSystem,sistemaQ::QuantumSystem;number_of_states=1,kargs...)
        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
        D=Nmax*(Int(2*j)+1)
        if number_of_states==1
            data = zeros(Complex{Float64},D)
        else
            data=zeros(Complex{Float64},D,number_of_states)
        end
        return randomState!(sistemaC,sistemaQ,data;kargs...)
    end
    function random_cₖ_generator(sistemaC::ClassicalSystems.ClassicalSystem,sistemaQ::QuantumSystem;
                                σ=0.3,ε=-0.5,envelope::UnivariateDistribution=Distributions.Normal(ε,σ),ensemble::Symbol=:GUE,
                                rₖ_distribution::UnivariateDistribution=
                                if ensemble==:GUE
                                        Exponential(0.91)
                                elseif ensemble==:GOE
                                        Chisq(1)
                                else
                                    error("ensemble must be :GOE or :GUE, or pass rₖ_distribution and phases directly")
                                end,
                                phases=if ensemble==:GUE
                                        :complex
                                elseif ensemble==:GOE
                                        :real
                                else
                                    error("ensemble must be :GOE or :GUE, or pass rₖ_distribution and phases directly")
                                end,divide_by_DoS=true)

        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
        ρ(ε)=Distributions.pdf(envelope,ε)
        ν(ε)=Dicke.density_of_states(sistemaC,j=j,ε=ε)
        rₖ()=Distributions.rand(rₖ_distribution)
        rphase()=if phases==true || phases ==:complex
            exp(im*Random.rand()*2*pi)
        elseif phases==:real
            rand((1,-1))
        else
            1
        end
        function cₖ(ε;number=1,cache=nothing)
            if divide_by_DoS
                dos=ν(ε)
            else
                dos=1
            end
            if number>1 && cache==nothing
                cache=zeros(ComplexF64,number)
            end
            if dos==0
                if number==1
                    return 0
                else
                    cache*=0
                    return cache
                end
            else
                f()=sqrt(rₖ()*ρ(ε)/dos)*rphase()
                if number==1
                    return f()
                else
                    cache.=f.()
                    return cache
                end

            end
        end
        return cₖ
    end

    function randomState!(sistemaC::ClassicalSystems.ClassicalSystem,sistemaQ::QuantumSystem,data;eigenenergies,eigenstates,σ=0.3,ε=-0.5,
        envelope::UnivariateDistribution=Distributions.Normal(ε,σ),tol=1e-5,
        parity::Union{Nothing,typeof(-),typeof(+)}=nothing,
        parities::AbstractArray{<:Integer,1}=if parity==nothing
            Int8[]
        else @warn "If you pass parity=+-, it is more efficient for you to compute parities=eigenstate_parities(sistemaQ,eigenstates), and pass it to randomState!(... parities=parities, ...). This prevents repeated calls to eigenstate_parities wich can be slow."
            eigenstate_parities(sistemaQ,eigenstates) end,kargs...)
        if length(size(data))==2
           for i in 1:size(data)[2]
                randomState!(sistemaC,sistemaQ,view(data,:,i);
                                eigenenergies=eigenenergies,eigenstates=eigenstates,σ=σ,ε=ε,envelope=envelope,tol=tol,parity=parity,parities=parities,kargs...)

           end
           return data
        end
        ω₀,ω,γ,j,Nmax=get_params(sistemaQ)
        cₖ=random_cₖ_generator(sistemaC,sistemaQ;ε=ε,σ=σ,envelope=envelope,kargs...)
        data.=0.0
        nrm=0.0
        ε0,εf=Distributions.quantile(envelope,tol/2),Distributions.quantile(envelope,1-tol/2)
        if parity!=nothing
            par=parity(1)
            is_correct_parity = k -> parities[k]==par #debes hacerlo con anonymus functions si no julia tiene un bug bien raro https://github.com/JuliaLang/julia/issues/15602

        else
            is_correct_parity = k -> true #idem
        end
        for (k,E) in enumerate(eigenenergies)

            ε=E/j
            if ε0<=ε<=εf && is_correct_parity(k)
                ck=cₖ(ε)
                nrm+= abs2(ck)
                data.+=ck.*view(eigenstates,:,k)
            end

        end
        data./=sqrt(nrm)
        return data
    end



    function randomCoherentStatesInEnergyShell(sistemaC::ClassicalSystems.ClassicalSystem,sistemaQ::QuantumSystem;ε::Real,N,dt=3,tol=1e-5,cachedata=nothing)
        ω₀,ω,γ,j,Nmax=sistemaQ.params

        if cachedata===nothing
            cachedata=zeros(Complex,(Int(2j)+1)*(Nmax))
        end
        cachedata.=0.0im
        s=Dicke.classicalPathRandomSampler(sistemaC,ε=ε,dt=dt)
        for i in 1:N
            θ=rand()*2*π
            p=s()
            coherentBCE!(sistemaQ,p,cachedata,extra_phase=exp(im*θ),tol=tol,add_state_to_data=true)
        end
        cachedata./=norm(cachedata)
        return cachedata
     end
end
