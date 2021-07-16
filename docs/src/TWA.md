# DickeModel.TWA
```@meta
CurrentModule = DickeModel.TWA
```
## Monte-Carlo systems
```@docs
monte_carlo_integrate
MonteCarloSystem
mcs_chain
```
## Wigner distribution of coherent states
```@docs
PhaseSpaceDistribution
coherent_Wigner_HW
coherent_Wigner_SU2
coherent_Wigner_HWxSU2
```
## Averages
```@docs
mcs_for_averaging
average
```
## Variances
```@docs
mcs_for_variance
variance
```
## Survival probability
```@docs
mcs_for_survival_probability
survival_probability
```
## Matrix distributions
```@docs
mcs_for_distributions
calculate_distribution
```


## [Weyl symbols](@id DickeModel.TWA.Weyl)
The submodule `TWA.Weyl` generates classical phase-space
observables that may be passed as the argument `observable` of [`average`](@ref average), [`calculate_distribution`](@ref calculate_distribution), etc. All the following functions return [SymEngine](https://juliahub.com/docs/SymEngine) objects,
so they may be operated as if they were numbers. 
```@autodocs
Modules = [TWA.Weyl]
Order   = [:type,:function]
```