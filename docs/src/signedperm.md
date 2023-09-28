# Signed permutations
```@index
Pages=["signedperm.md"]
```

```@docs
SignedPerms
SPerm
SPerm(::Integer...)
SignedPerms.@sperm_str(::String)
Perms.Perm(::SPerm)
signs
Perms.permute(::AbstractVector,::SPerm)
Perms.orbit(::SPerm,::Integer)
Perms.order(::SPerm)
Perms.cycles(::SPerm)
Perms.cycletype(::SPerm)
Matrix(::SPerm)
SPerm(::AbstractMatrix{<:Integer})
SPerm(::AbstractVector,::AbstractVector)
sstab_onmats
SPerm(::AbstractMatrix,::AbstractMatrix)
hyperoctaedral_group
SignedPerms.SPerm_onmats
```
