# Coxeter groups
```@index
Pages=["coxgroups.md"]
```
```@docs
CoxGroups
isleftdescent(::CoxSym,::Any,::Int)
firstleftdescent
leftdescents(::CoxeterGroup,w)
word(::CoxeterGroup,::Any)
length(::CoxeterGroup,w)
longest
elements(::CoxeterGroup,::Int)
PermGroups.Groups.words(::CoxeterGroup{T},::T) where T
PermGroups.Groups.words(::CoxeterGroup)
bruhatless
coxeter_matrix
cartan(::AbstractMatrix)
coxeter_symmetric_group
reflection_subgroup(::CoxSym,::AbstractVector{Int})
cartan(::CoxSym)
reduced
standard_parabolic_class
coxeter_group(::AbstractMatrix)
Poset(::CoxeterGroup,w=longest(W))
```
