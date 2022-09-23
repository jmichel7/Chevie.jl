# Coxeter groups
```@index
Modules=[CoxGroups]
```
```@docs
CoxGroups
isleftdescent(::CoxSym,::Any,::Int)
firstleftdescent
leftdescents(::CoxeterGroup,w)
reduced
word(::CoxeterGroup,::Any)
length(::CoxeterGroup,w)
elements(::CoxeterGroup,::Int)
bruhatless
CoxSym
reflection_subgroup(::CoxSym,::AbstractVector{Int})
longest
coxmat
cartan(::AbstractMatrix)
cartan(::CoxSym)
standard_parabolic_class
coxgroup(::AbstractMatrix)
PermGroups.Groups.words(::CoxeterGroup{T},::T) where T
PermGroups.Groups.words(::CoxeterGroup)
```
