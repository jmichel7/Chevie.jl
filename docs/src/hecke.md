# Hecke algebras
```@index
Pages=["hecke.md"]
```

```@docs
HeckeAlgebras
hecke(::Group,::Vector{<:Vector{C}}) where C
CharTable(::HeckeAlgebra)
central_monomials
class_polynomials
char_values
alt
α(::HeckeTElt)
schur_elements
factorized_schur_element
factorized_schur_elements
HeckeAlgebras.FactSchur
representation(::HeckeAlgebra,::Integer)
representations(::Union{HeckeAlgebra,HeckeCoset})
isrepresentation
reflrep(::HeckeAlgebra)
HeckeCoset
hecke(::HeckeCoset)
hecke(::Spets,::HeckeAlgebra)
```
