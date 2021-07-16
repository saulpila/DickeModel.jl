# Dicke.UPOS
```@meta
CurrentModule = Dicke.UPOS
```
## Periodic Orbit (PO) type
```@docs 
PO
PO(::ClassicalDickeSystem, ::AbstractVector{<:Real})
Base.:(==)(::PO,::PO)
Base.in(::AbstractVector{<:Number},::PO)
```
## Integrals over POs
```@docs 
integrate_PO
average_over_PO
average_over_PO_QP
average_over_PO_qp
```

## Properties of POs
```@docs 
action
jx_PO_average
jy_PO_average
jz_PO_average
lyapunov
energy
```

## Operations on POs
### Mirror operations
```@docs 
mirror_Qq
mirror_Pp
```
### Coordinate extraction

```@docs 
po_coordinates
```
The following methods are shorthand implementations of [`po_coordinates`](@ref). They are useful to
plot a periodic orbit. For example,
`plot(QP(po))` plots the orbit in the 
atomic plane.

```@docs 
QP
qp
QPp
QPq
```

## Finding POs
```@docs 
approximate_period
find_p_zero
```

### Monodromy methods
```@docs 
monodromy_method_constant_period
monodromy_method_constant_energy
```

### PO familly followers
```@docs
follow_PO_family_from_period
follow_PO_family_from_energy

family_A
family_B
```
## Scarring measure
```@docs
scarring_measure
overlap_of_tube_with_homogenous_state
```
