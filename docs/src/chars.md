# Classes/characters of reflection groups
```@index
Pages=["chars.md"]
```
```@docs
Chars
CharTable
decompose(::CharTable,::AbstractVector)
scalar_product(::CharTable,::AbstractVector,::AbstractVector)
on_chars
charinfo
charnames(io::IO, W::Union{Coset, Group})
classnames
classinfo
fakedegree(::Any,::Any,::Any)
fakedegrees
representation(::Union{Chars.Hastype,FiniteCoxeterGroup},::Integer)
representations(::Union{Spets, FiniteCoxeterGroup, PermRootGroup})
induction_table
j_induction_table
J_induction_table
schur_functor
detPerm
conjPerm
```
