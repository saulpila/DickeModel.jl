This package is a result of more than two years of investigation of the [Dicke Model](https://en.wikipedia.org/wiki/Dicke_model). It contains numerical methods that were used in the following publications:

1. *Quantum localization measures in phase space*.  Physical Review E **103** 052214 (2021) DOI: [10.1103/physreve.103.052214](https://doi.org/10.1103/physreve.103.052214)
2. *Ubiquitous quantum scarring does not prevent ergodicity*.  Nature Communications **12** 852 (2021) DOI: [10.1038/s41467-021-21123-5](https://doi.org/10.1038/s41467-021-21123-5)
3. *Quantum scarring in a spin-boson system: fundamental families of periodic orbits*.  New Journal of Physics **23** 033045 (2021) DOI: [10.1088/1367-2630/abd2e6](https://doi.org/10.1088/1367-2630/abd2e6)
4. *Quantum vs classical dynamics in a spin-boson system: manifestations of spectral correlations and scarring*.  New Journal of Physics **22** 063036 (2020) DOI: [10.1088/1367-2630/ab8ef8](https://doi.org/10.1088/1367-2630/ab8ef8)
5. *Positive quantum Lyapunov exponents in experimental systems with a regular classical limit*.  Physical Review E **101** 010202(R) (2020) DOI: [10.1103/PhysRevE.101.010202](https://doi.org/10.1103/PhysRevE.101.010202)

It is split into several submodules:
* [`ClassicalDicke`](@ref DickeModel.ClassicalDicke) allows to compute classical dynamics of the Dicke model, including a wide range of semiclassical aproximations to quantum properties.
* [`ClassicalSystems`](@ref DickeModel.ClassicalSystems) provides a general framework for computing classical Hamiltonian dynamics, inlcuding Lyapunov exponents. It is mostly used for the Dicke model, but in principle it can be expanded to other Hamiltonians.
* [`DickeBCE`](@ref DickeModel.DickeBCE) provides multiple functions for analyzing the quantum Dicke model. It uses an efficient basis known as the Efficient Coherent Basis (BCE). See Refs. above.
*  [`UPOs`](@ref DickeModel.UPOs) contains a set of functions to find unstable periodic orbits (UPOs) in the classical Dicke model. (See Ref. 3)
*  [`TWA`](@ref DickeModel.TWA) allows to perform semiclassical calculations using the truncated Wigner approximation (TWA). (See Refs 2,5)
* [`EnergyShellProjections`](@ref DickeModel.EnergyShellProjections) contains a set of functions to integrate functions over the classical energy shells of the Dicke Model. It also contains specialized functions to compute these integrals for the Husimi functions of quantum states, which define phase-space localization measures known as Rényi occupations (See Refs. 1 and 2)
* [`ClassicalLMG`](@ref DickeModel.ClassicalLMG) provides very basic functions for the classical Lipkin-Meshkov-Glick model. (See Ref. 5)
* [`PhaseSpaces`](@ref DickeModel.PhaseSpaces) provides some canonical transformations of the Bloch-Sphere.
