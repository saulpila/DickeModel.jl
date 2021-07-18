# DickeModel.DickeBCE
For examples of usage, go to  [Examples for DickeBCE](@ref).
```@meta
CurrentModule = DickeModel.DickeBCE
```
```@docs
QuantumDickeSystem
QuantumDickeSystem()
dimension
```
## Diagonalization
```@docs
diagonalization
eigenenergies
eigenstate_parities
```
### Evolution
```@docs
evolve
survival_probability
```
## Operators
```@docs
DickeBCE.Jx
DickeBCE.Jz
hamiltonian
parity_operator
```
### Participation Ratio
```@docs
participation_ratio
factor_R_of_coherent_state
```

## Coherent states
```@docs
coherent_state!
coherent_state
```
### Overlaps 
```@docs
coherent_overlap
husimi
husimi_of_coherent
```

## Random states 
```@docs
random_cₖ_generator
random_state!
random_state
random_coherent_states_in_energy_shell
```
## Wigner functions
```@docs
Wigner 
WignerProjQP
WignerProjqp
```