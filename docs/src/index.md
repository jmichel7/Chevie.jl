# Gapjm.jl Documentation
```@docs
Gapjm
```
```@contents
```
# Perms.jl Documentation
```@docs
Perms
Perm
largest_moved_point
smallest_moved_point
order
Groups.orbit(a::Perm,i::Integer,check=false)
Groups.orbits(a::Perm,domain=1:length(a.d);trivial=true,check=false)
cycles
cycletype
sign
permuted
restricted
```
# Groups.jl Documentation
```@docs
Groups
orbit
orbits
transversal
centralizer
word(G::Group,w)
elements
length
class_reps
minimal_words
```
# PermGroups.jl Documentation
```@docs
PermGroups
base
centralizers
transversals
symmetric_group
```
# Cycs.jl Documentation
```@docs
Cycs
galois
ER
Quadratic
```
# Pols.jl Documentation
```@docs
Pols
divrem
divrem1
gcd
```
# Mvps.jl Documentation
```@docs
Mvps
variables
coefficients
valuation
degree
Mvp(;arg...)
```
# CoxGroups.jl Documentation
```@docs
CoxGroups
reduced
word(W::CoxeterGroup,w)
bruhatless
coxsym
longest
nref
```
# Weyl.jl Documentation
```@docs
Weyl
cartan
two_tree
reflection_subgroup
coxgroup
rootdatum
```
# PermRoot.jl Documentation
```@docs
PermRoot
reflection
```
# Hecke.jl Documentation
```@docs
Hecke
hecke
```
# KL.jl Documentation
```@docs
KL
KLPol
Tbasis
```
# Garside.jl Documentation
```@docs
Garside
left_divisors
DualBraidMonoid
fraction
word(b::Garside.GarsideElm)
representative_operation
centralizer_generators
shrink
```
# Chars.jl Documentation
```@docs
Chars
charinfo
classinfo
fakedegree
fakedegrees
representation
representations
```
# Uch.jl Documentation
```@docs
Uch
UnipotentCharacters
```
# Ucl.jl Documentation
```@docs
Ucl
UnipotentClasses
ICCTable
```
# HasType.jl Documentation
```@docs
```
# Symbols.jl Documentation
```@docs
Symbols
shiftβ
βset
partβ
ranksymbol
fegsymbol
```
# Util.jl Documentation
```@docs
Util
groupby
constant
blocks
format
prime_residues
phi
primitiveroot
echelon!
```
# ModuleElts.jl Documentation
```@docs
ModuleElts
```
# Cycpols.jl Documentation
```@docs
CycPols
```
# Posets.jl Documentation
```@docs
Posets
lcm_partitions
gcd_partitions
transitive_closure
linear_extension
hasse
incidence
reverse
partition
restricted(p::Poset,ind::AbstractVector{<:Integer})
is_join_lattice
is_meet_lattice
Poset(W::CoxeterGroup,w=longest(W))
```

# Dictionary from Chevie
The dictionary from Chevie is as follows:
```
     CartanMat("A",5)                      →  cartan(:A,5) 
     CoxeterElements(W[,l])                →  elements(W[,l])
     CoxeterGroup("A",5)                   →  coxgroup(:A,5) 
     CoxeterGroupHyperoctaedralGroup(n)    →  CoxHyperoctaedral(n)
     CoxeterGroupSymmetricGroup(n)         →  CoxSym(n)
     CoxeterLength(W,w)                    →  length(W,w)
     CoxeterWord(W,w)                      →  word(W,w)
     ElementWithInversions(W,l)            →  with_inversions(W,l)
     FiniteCoxeterTypeFromCartanMat(m)     →  type_cartan(m) 
     FirstLeftDescending(W,w)              →  firstleftdescent(W,w)
     ForEachElement(W,f)                   →  for w in W f(w) end 
     HyperplaneOrbits                      →  hyperplane_orbits
     IndependentRoots                      →  independent_roots
     Inversions                            →  inversions 
     IsLeftDescending(W,w,i)               →  isleftdescent(W,w,i)
     LeftDescentSet(W,w)                   →  leftdescents(W,w)
     LongestCoxeterElement(W)              →  longest(W)
     MatXPerm(W,p)                         →  matX(W,p)
     .orbitRepresentativeElement           →  simple_conjugating_element
     .orbitRepresentative                  →  simple_representatives
     PositionClass                         →  position_class
     PrintDiagram(W)                       →  Diagram(W) 
     ReducedExpressions(W,w)               →  reduced_words(W,w)
     ReducedInRightCoset(W,w)              →  reduced(W,w)
     ReducedRightCosetRepresentatives(W,H) →  reduced(H,W)
     ReflectionCharacter                   →  reflchar
     ReflectionDegrees(W)                  →  degrees(W) 
     ReflectionEigenvalues                 →  refleigen
     ReflectionLength(W,w)                 →  reflength(W,w)
     Reflection                            →  reflection 
     Reflections                           →  reflections
     ReflectionSubgroup                    →  reflection_subgroup
     ReflectionType                        →  refltype
     RightDescentSet(W,w)                  →  rightdescents(W,w)
     RootsCartan(m)                        →  roots(m) 
     SemiSimpleRank(W)                     →  coxrank(W)
     Size(W)                               →  length(W) 
     StandardParabolic                     →  standard_parabolic
     TwoTree(m)                            →  twotree(m) 
     W.matgens[i]                          →  matX(W,i)
     W.N                                   →  nref(W)
     W.orbitRepresentative[i]              →  simple_representative(W,i) 
```
