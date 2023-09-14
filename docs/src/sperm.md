# Signed permutations
```@index
Pages=["sperm.md"]
```

```@docs
SPerms
SPerm
Perm(::SPerm)
@sperm_str
permute(::AbstractVector,::SPerm)
orbit(::SPerm,::Integer)
order(::SPerm)
cycles(::SPerm)
cycletype(::SPerm)
Matrix(::SPerm)
coxeter_hyperoctaedral_group
reflection_subgroup(::CoxHyp,::AbstractVector{Int})
sstab_onmats
SPerms.SPerm_onmats
```
