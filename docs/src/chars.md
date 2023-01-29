# Classes/characters of reflection groups
```@index
Pages=["chars.md"]
```
```@docs
Chars
CharTable
on_chars
charinfo
charnames(io::IO, W::Union{Coset{T, TW}, Group}) where {T, TW}
classnames
classinfo
fakedegree
fakedegrees
representation(::Union{Chars.Hastype,FiniteCoxeterGroup},::Integer)
representations(::Union{Spets, FiniteCoxeterGroup, PermRootGroup})
InductionTable
jInductionTable
JInductionTable
detPerm
conjPerm
```