# [DickeModel.jl](https://github.com/saulpila/DickeModel.jl)

*A toolkit for the quantum and classical Dicke model in Julia.*

```@raw html
<img src="banner.png" id="bannerimage" width="100%" style="margin-bottom: 1em;"/>
<script> 
    document.getElementById("bannerimage").style.filter = "hue-rotate("+Math.random()+"turn)"
</script> 
``` 

This package contains numerical methods that allow to efficiently compute properties of the quantum and classical [Dicke model](https://en.wikipedia.org/wiki/Dicke_model), a fundamental model in quantum optics describing atoms interacting with light.

## Installation
 To install the package, use the following command inside the Julia REPL:
```julia
using Pkg
Pkg.add("DickeModel")
```
    
To load the package, run
```julia
using DickeModel
```

## Getting started

There are many things you can do using the different submodules:

* [`ClassicalDicke`](@ref DickeModel.ClassicalDicke) allows to compute classical dynamics of the Dicke model, including a wide range of properites of its phase space (used in [Pilatowsky2020](@cite), [Pilatowsky2021](@cite), [Pilatowsky2021NatCommun](@cite), [Pilatowsky2021Identification](@cite), [Villasenor2020](@cite),  [Villasenor2021](@cite)).\
  See [Examples for ClassicalDicke](@ref).
* [`DickeBCE`](@ref DickeModel.DickeBCE) provides multiple functions for analyzing the quantum Dicke model. It uses an efficient basis known as the Efficient Coherent Basis (BCE) [Bastarrachea2014PSa](@cite), [Bastarrachea2014PSb](@cite) for exact computations, and it also contains several semilcassical methods (used in [Pilatowsky2020](@cite), [Pilatowsky2021](@cite), [Pilatowsky2021NatCommun](@cite), [Pilatowsky2021Identification](@cite), [Villasenor2020](@cite),  [Villasenor2021](@cite)).\
  See [Examples for DickeBCE](@ref).
* [`UPOs`](@ref DickeModel.UPOs) contains a set of functions to find unstable periodic orbits (UPOs) in the classical Dicke model and to study quantum scars (used in  [Pilatowsky2021](@cite), [Pilatowsky2021NatCommun](@cite), [Pilatowsky2021Identification](@cite)).\
  See [Examples for UPOs](@ref).
* [`TWA`](@ref DickeModel.TWA) allows to perform semiclassical calculations using the truncated Wigner approximation (TWA) (used in [Pilatowsky2020](@cite), [Villasenor2020](@cite)).\
  See [Examples for TWA](@ref).
* [`EnergyShellProjections`](@ref DickeModel.EnergyShellProjections) contains a set of functions to integrate functions over the classical energy shells of the Dicke Model. It also contains specialized functions to compute these integrals for the Husimi functions of quantum states, which define phase-space localization measures known as Rényi occupations [Villasenor2021](@cite) (used in [Pilatowsky2021NatCommun](@cite), [Pilatowsky2021Identification](@cite), [Villasenor2021](@cite)).\
  See [Examples for EnergyShellProjections](@ref).
* [`ClassicalLMG`](@ref DickeModel.ClassicalLMG) provides very basic functions for the classical Lipkin-Meshkov-Glick model.\
  See [Examples for ClassicalLMG](@ref).

The following modules provide basic functionallity for the rest of the modules:
* [`ClassicalSystems`](@ref DickeModel.ClassicalSystems) provides a general framework for computing classical Hamiltonian dynamics, inlcuding Lyapunov exponents. It is mostly used for the Dicke model, but in principle it can be expanded to other Hamiltonians.
* [`PhaseSpaces`](@ref DickeModel.PhaseSpaces) provides some canonical transformations of the Bloch-Sphere.
