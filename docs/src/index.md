# Chevie
```@docs
Chevie
```

Chevie uses its rich total infrastructure to provide extensions to
several of its infrastructure packages.

#### Extensions to Laurent and Puiseux polynomials
```@docs
FFfac.factor(::Pol{FFE{p}}, Any) where p
Fact.factor(::Pol{T}) where T<:Union{Integer, Rational{<:Integer}}
factor(::Mvp{T, N}) where {T, N}
```
#### Arithmetic and finite fields
```@docs
LaurentPolynomials.valuation(::Integer,::Integer)
FiniteFields.FFE{p}(Cyc)where p
```

#### Extensions to groups
```@docs
abelian_gens
abelian_invariants
Combinat.blocks(::Group,::Integer)
Base.parent(::Group)
```
#### Extensions to linear algebra
```@docs
eigmat
```
#### Useful macros
```@docs
@forward
```
