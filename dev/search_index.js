var documenterSearchIndex = {"docs":
[{"location":"classicaldicke/#Dicke.ClassicalDicke","page":"ClassicalDicke","title":"Dicke.ClassicalDicke","text":"","category":"section"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"Modules = [Dicke.ClassicalDicke]\nOrder   = [:function, :type]","category":"page"},{"location":"classicaldicke/#Dicke.ClassicalDicke.ClassicalSystem-Tuple{}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.ClassicalSystem","text":"function ClassicalSystem(;ω₀,ω,γ)\n\nGenerates an instance of ClassicalSystems.ClassicalSystem which represents the  classical Dicke model with the given parameters ω_0, ω, and γ. See Eq. (5) of Ref. Saúl Pilatowsky-Cameo, David Villaseñor, Miguel A. Bastarrachea-Magnani, Sergio Lerma-Hernández, Lea F. Santos, Jorge G. Hirsch (2021).\n\nThe returned value may be passed to all functions in this module that require an instance of ClassicalSystems.ClassicalSystem, as well as functions in other modules, such as  ClassicalSystems.integrate.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.Point-Tuple{Dicke.ClassicalSystems.ClassicalSystem}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.Point","text":"function Point(system::ClassicalSystems.ClassicalSystem;Q,P,p,ε,signo::Union{typeof(-),typeof(+)}=+)\n\nReturns a list [Q,q,P,p], where q is calculated with ClassicalDicke.q_of_ε. See that function for details on the arguments.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.Pointθϕ-Tuple{}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.Pointθϕ","text":"function Pointθϕ(;θ,ϕ,q,p)\n\nReturns a list [Q,q,P,p], where Q and P are calculated from θ and ϕ using PhaseSpaces.Q_of_θϕ and  PhaseSpaces.P_of_θϕ.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.WignerHWxSU2_fixed_ε-Tuple{Dicke.ClassicalSystems.ClassicalSystem}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.WignerHWxSU2_fixed_ε","text":"function WignerHWxSU2_fixed_ε(system::ClassicalSystems.ClassicalSystem;Q,P,p,ε,j,signoq=+)\n\nReturns an instance of TruncatedWignerApproximation.PhaseSpaceDistribution that samples points from the classical energy shell  using a Random Walk Metropolis-Hastings algorithm implemented in Mamba.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.density_of_states-Tuple{Dicke.ClassicalSystems.ClassicalSystem}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.density_of_states","text":"function density_of_states(system::ClassicalSystems.ClassicalSystem;j,ε)\n\nReturns the semiclassical density of states (DoS) ν(ϵ), in units of 1ϵ. This is computed with an expression similar to Eq. (19) of Ref. M. A. Bastarrachea-Magnani, S. Lerma-Hernández, J. G. Hirsch (2014), where we add an additional factor of j to have units of 1ϵ instead of 1E, and the integral is performed with a change of variable.\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\nj is the value of j\nε is the scaled energy ϵ=Ej\n\nSee Plotting the density of states for an example.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.discriminant_of_q_solution-Tuple{Dicke.ClassicalSystems.ClassicalSystem}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.discriminant_of_q_solution","text":"function discriminant_of_q_solution(system::ClassicalSystems.ClassicalSystem; Q,P,p,ε)\n\nReturns the discriminant of the second degree equation in q given by\n\n    h_textcl(QqPp)=epsilon\n\nwhere h_textcl is given by Eq. (5) of Ref. Saúl Pilatowsky-Cameo, David Villaseñor, Miguel A. Bastarrachea-Magnani, Sergio Lerma-Hernández, Lea F. Santos, Jorge G. Hirsch (2021).\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\nQ, P, p, and ε are the values of Q, P, p, and epsilon, respectively.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.energy_shell_volume-Tuple{Dicke.ClassicalSystems.ClassicalSystem, Any}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.energy_shell_volume","text":"function energy_shell_volume(system::ClassicalSystems.ClassicalSystem;ε)\n\nReturns the volume of the classical energy shell in the phase space, that is, mathcalV(mathcalM_ϵ) in Eq. (27) of Ref. Saúl Pilatowsky-Cameo, Jorge Chávez-Carlos, Miguel A. Bastarrachea-Magnani, Pavel Stránský, Sergio Lerma-Hernández, Lea F. Santos, Jorge G. Hirsch (2020).\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\nε is the scaled energy ϵ=Ej\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.energy_width_of_coherent_state-Tuple{Dicke.ClassicalSystems.ClassicalSystem, Any, Real}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.energy_width_of_coherent_state","text":"function energy_width_of_coherent_state(system::ClassicalSystems.ClassicalSystem,x,j::Real)\n\nReturns the energy width sigma of the coherent state left  mathbfxright rangle, in units of epsilon. This  quantity is given by sigma_Dj with sigma_D as in App. A of Ref. Sergio Lerma-Hernández, Jorge Chávez-Carlos, Miguel A. Bastarrachea-Magnani, Lea F. Santos, Jorge G. Hirsch (2018). \n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\nx is the coordinate mathbfx of the coherent state in the format [Q,q,P,p].\nj is the value of j.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.hamiltonian-Tuple{Dicke.ClassicalSystems.ClassicalSystem}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.hamiltonian","text":"function hamiltonian(system::ClassicalSystems.ClassicalSystem)\n\nReturns a classical Hamiltonian function h(x) where x=[Q,q,P,p], which is given by Eq. (5) of Ref. Saúl Pilatowsky-Cameo, David Villaseñor, Miguel A. Bastarrachea-Magnani, Sergio Lerma-Hernández, Lea F. Santos, Jorge G. Hirsch (2021).\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.maximum_P_for_ε-Tuple{Dicke.ClassicalSystems.ClassicalSystem, Any}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.maximum_P_for_ε","text":"function maximum_P_for_ε(system::ClassicalSystems.ClassicalSystem,ε)\n\nComputes the maximum value of the parameter P accessible to the system at energy epsilon.\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\nε is the scaled energy ϵ=Ej.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.maximum_Q_for_ε-Tuple{Dicke.ClassicalSystems.ClassicalSystem, Any}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.maximum_Q_for_ε","text":"function maximum_Q_for_ε(system::ClassicalSystems.ClassicalSystem,ε)\n\nSee ClassicalDicke.maximum_P_for_ε.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.minimum_energy-Tuple{Dicke.ClassicalSystems.ClassicalSystem}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.minimum_energy","text":"function minimum_energy(system::ClassicalSystems.ClassicalSystem)\n\nReturns the energy of the ground-state coordinate given by ClassicalDicke.minimum_energy_point.\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\n\nNote: This function currently only works for the supperadiant phase.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.minimum_energy_point","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.minimum_energy_point","text":"function minimum_energy_point(system::ClassicalSystems.ClassicalSystem,Qsign::Union{typeof(-),typeof(+)}=+)\n\nReturns the ground-state coordinate, that is, mathbfx_textGS below Eq. (7) of Ref. Saúl Pilatowsky-Cameo, David Villaseñor, Miguel A. Bastarrachea-Magnani, Sergio Lerma-Hernández, Lea F. Santos, Jorge G. Hirsch (2021).\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\nQsign toggles the sign of the Q coordinate, that is, + for mathbfx_textGS and  - for widetildemathbfx_textGS.\n\nNote: This function currently only works for the supperadiant phase.\n\n\n\n\n\n","category":"function"},{"location":"classicaldicke/#Dicke.ClassicalDicke.minimum_ε_for-Tuple{Dicke.ClassicalSystems.ClassicalSystem}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.minimum_ε_for","text":"function minimum_ε_for(system::ClassicalSystems.ClassicalSystem;Q=:nothing,q=:nothing,P=:nothing,p=:nothing)\n\nReturns the minimum energy epsilon when constraining the system to three fixed values of the coordinates Q, q, P, p.\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\nYou may pass either (QqP) or (qPp). The other combinanations are not implemented.\n\nThis function can be especially useful to draw contours of the available phase space (see Drawing contours of the available phase space)\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.normal_frequency","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.normal_frequency","text":"function normal_frequency(system::ClassicalSystems.ClassicalSystem,signo::Union{typeof(-),typeof(+)}=+)\n\nReturns the ground-state normal frequency, that is, Omega_epsilon_textGS^AB at the bottom of page 3 of Ref. Saúl Pilatowsky-Cameo, David Villaseñor, Miguel A. Bastarrachea-Magnani, Sergio Lerma-Hernández, Lea F. Santos, Jorge G. Hirsch (2021).\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\nsigno is - for Omega^A and + for Omega^B.\n\nNote: This function currently only works for the supperadiant phase.\n\n\n\n\n\n","category":"function"},{"location":"classicaldicke/#Dicke.ClassicalDicke.phase_space_dist_squared-Tuple{Any, Any}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.phase_space_dist_squared","text":"function phase_space_dist_squared(x,y)\n\nReturns the phase-space distance d_mathcalM(mathbfxmathbfy) (See App. C of Ref Saúl Pilatowsky-Cameo, David Villaseñor, Miguel A. Bastarrachea-Magnani, Sergio Lerma-Hernández, Lea F. Santos, Jorge G. Hirsch (2021)), where x and y are in the form [Q,q,P,p].\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/#Dicke.ClassicalDicke.q_of_ε-Tuple{Dicke.ClassicalSystems.ClassicalSystem}","page":"ClassicalDicke","title":"Dicke.ClassicalDicke.q_of_ε","text":"function q_of_ε(system::ClassicalSystems.ClassicalSystem;Q,P,p,ε,signo::Union{typeof(-),typeof(+)}=+,returnNaNonError=true)\n\nReturns the solutions q_pm of the second degree equation in q given by\n\n    h_textcl(QqPp)=epsilon\n\nwhere h_textcl is given by Eq. (5) of Ref. Saúl Pilatowsky-Cameo, David Villaseñor, Miguel A. Bastarrachea-Magnani, Sergio Lerma-Hernández, Lea F. Santos, Jorge G. Hirsch (2021).\n\nArguments:\n\nsystem should be generated with ClassicalDicke.ClassicalSystem.\nQ,P,p,ε are values of Q,P,p,epsilon, respectively.\nsigno is + for q_+ and - for q_-\nIf returnNaNonError is true, then NaN is returned if there is no solutions. If it is false, and error is raised.\n\n\n\n\n\n","category":"method"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"push!(LOAD_PATH,\"../../src\")\nusing Dicke","category":"page"},{"location":"classicaldicke/#Examples","page":"ClassicalDicke","title":"Examples","text":"","category":"section"},{"location":"classicaldicke/#Drawing-contours-of-the-available-phase-space","page":"ClassicalDicke","title":"Drawing contours of the available phase space","text":"","category":"section"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"We may use the function ClassicalDicke.minimum_ε_for to draw the contour of the available phase space on the variables (QP).","category":"page"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"using Dicke\nusing Plots\nsystem=  ClassicalDicke.ClassicalSystem(ω=1, γ=1, ω₀=1)\nQs = Ps = -2:0.01:2\nεgs = ClassicalDicke.minimum_energy(system)\ncontour(Qs, Ps, \n        (Q,P) -> ClassicalDicke.minimum_ε_for(system, p=0, P=P, Q=Q),\n        levels=10, clim=(εgs,1), xlabel=\"Q\", ylabel=\"P\")\nsavefig(\"contourQP.svg\"); nothing #hide","category":"page"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"(Image: )","category":"page"},{"location":"classicaldicke/#Plotting-the-density-of-states","page":"ClassicalDicke","title":"Plotting the density of states","text":"","category":"section"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"Here is a plot of the semiclassical density of states","category":"page"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"using Dicke\nusing Plots\nsystem = ClassicalDicke.ClassicalSystem(ω=1, γ=1, ω₀=1)\nν(ε) = ClassicalDicke.density_of_states(system, j=100, ε=ε)\nεgs = ClassicalDicke.minimum_energy(system)\nplot(ν, εgs, 2, xlabel=\"ε\", ylabel=\"Density of States\")\nplot!(key=false) #hide\nsavefig(\"density_of_states.svg\"); nothing #hide","category":"page"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"(Image: )","category":"page"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"This is precisely the red line in Fig. A1. of Ref. David Villaseñor, Saúl Pilatowsky-Cameo, Miguel A Bastarrachea-Magnani, Sergio Lerma, Lea F Santos, Jorge G Hirsch (2020).","category":"page"},{"location":"classicaldicke/#Drawing-a-Poincaré-surface","page":"ClassicalDicke","title":"Drawing a Poincaré surface","text":"","category":"section"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"Here is a way to draw a Poincaré surface for the Dicke model. We use ClassicalSystems.integrate to integrate a bunch of initial conditions. Using the callback system of DifferentialEquations, we save the points where p=0.","category":"page"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"using Dicke\nusing Plots\nusing DifferentialEquations\n\nsystem=ClassicalDicke.ClassicalSystem(ω=0.8, γ=0.8, ω₀=1)\nmplot=scatter(fmt=:png,key=false,markerstrokecolor = :transparent,\n  markersize=1,legend=false,size=(500,500),\n  color_palette=:darkrainbow,xlabel=\"Q\",ylabel=\"P\") \n\npts=Tuple{Float64, Float64}[] #a list of points (Q,P)\ncallback=ContinuousCallback((x,t,_)-> x[4], #when p=x[4] is 0,\n    state->push!(pts,(state.u[1],state.u[3])),  #save Q=u[1], P=u[3] \n    nothing;\n    save_positions=(false,false),abstol=1e-6)\nε=-1.35\nmaxQ = ClassicalDicke.maximum_Q_for_ε(system,ε)\nfor Q in 0:0.02:maxQ #for a bunch of initial Qs,\n        if ClassicalDicke.minimum_ε_for(system,P=0,p=0,Q=Q) > ε #we are outside the bounds\n            continue\n        end\n        initial_condition = ClassicalDicke.Point(system, ε=ε, P=0, p=0, Q=Q)\n        ClassicalSystems.integrate(system,u₀=initial_condition,\n            t=10000,cb=callback,save_everystep=false,tol=1e-9)\n        scatter!(pts)\n        empty!(pts)\nend\nmplot\nsavefig(\"poincare_surface.png\");nothing #hide","category":"page"},{"location":"classicaldicke/","page":"ClassicalDicke","title":"ClassicalDicke","text":"(Image: )","category":"page"},{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"","category":"page"},{"location":"classicalsystems/#Dicke.ClassicalSystems","page":"ClassicalSystems","title":"Dicke.ClassicalSystems","text":"","category":"section"},{"location":"classicalsystems/","page":"ClassicalSystems","title":"ClassicalSystems","text":"Modules = [Dicke.ClassicalSystems]\nOrder   = [:function, :type]","category":"page"},{"location":"classicalsystems/#Dicke.ClassicalSystems.integrate-Tuple{Dicke.ClassicalSystems.ClassicalSystem}","page":"ClassicalSystems","title":"Dicke.ClassicalSystems.integrate","text":"function integrate(sistema::ClassicalSystem;t,\n                           u₀,\n                           t₀=0.0,\n                           tol=1e-12,\n                           show_progress=false,\n                           save_intermediate_steps=nothing,\n                           saveat=Array{Float64,1}(),\n                           cb=nothing,\n                           get_fundamental_matrix=false,\n                           integrator_alg=TsitPap8(),\n                           use_big_numbers=false,\n                           integate_backwards=false,\n                           kargs...)\n\nTest ref Ref. Saúl Pilatowsky-Cameo, Jorge Chávez-Carlos, Miguel A. Bastarrachea-Magnani, Pavel Stránský, Sergio Lerma-Hernández, Lea F. Santos, Jorge G. Hirsch (2020).\n\nArguments:\n\nsystem is an instance of ClassicalSystems.ClassicalSystem.\n\n\n\n\n\n","category":"method"},{"location":"","page":"The Dicke.jl package","title":"The Dicke.jl package","text":"This package is a result of more than two years of investigation of the Dicke Model. It contains numerical methods that were used in the following publications:","category":"page"},{"location":"","page":"The Dicke.jl package","title":"The Dicke.jl package","text":"Quantum localization measures in phase space.  Physical Review E 103 052214 (2021) DOI: 10.1103/physreve.103.052214\nUbiquitous quantum scarring does not prevent ergodicity.  Nature Communications 12 852 (2021) DOI: 10.1038/s41467-021-21123-5\nQuantum scarring in a spin-boson system: fundamental families of periodic orbits.  New Journal of Physics 23 033045 (2021) DOI: 10.1088/1367-2630/abd2e6\nQuantum vs classical dynamics in a spin-boson system: manifestations of spectral correlations and scarring.  New Journal of Physics 22 063036 (2020) DOI: 10.1088/1367-2630/ab8ef8\nPositive quantum Lyapunov exponents in experimental systems with a regular classical limit.  Physical Review E 101 010202(R) (2020) DOI: 10.1103/PhysRevE.101.010202","category":"page"},{"location":"","page":"The Dicke.jl package","title":"The Dicke.jl package","text":"It is split into several submodules:","category":"page"},{"location":"","page":"The Dicke.jl package","title":"The Dicke.jl package","text":"ClassicalDicke allows to compute classical dynamics of the Dicke model, including a wide range of semiclassical aproximations to quantum properties.\nClassicalSystems provides a general framework for computing classical Hamiltonian dynamics, inlcuding Lyapunov exponents. It is mostly used for the Dicke model, but in principle it can be expanded to other Hamiltonians.\nDickeBCE provides multiple functions for analyzing the quantum Dicke model. It uses an efficient basis known as the Efficient Coherent Basis (BCE). See Refs. above.\nUPOS contains a set of functions to find unstable periodic orbits (UPOs) in the classical Dicke model. (See Ref. 3)\nTruncatedWignerApproximation allows to perform semiclassical calculations using the truncated Wigner approximation (TWA). (See Ref 2)\nFOTOCTWA is a small module providing functions to compute the fidelity out-of-order-time correlator using the TWA.  (See Ref 5)\nDickeHusimiProjections contains a set of functions to integrate functions over the classical energy shells of the Dicke Model. It also contains specialized functions to compute these integrals for the Husimi functions of quantum states, which define phase-space localization measures known as Rényi occupations (See Refs. 1 and 2)\nClassicalLMG provides very basic functions for the classical Lipkin-Meshkov-Glick model. (See Ref. 5)\nPhaseSpaces provides some canonical transformations of the Bloch-Sphere.","category":"page"}]
}
