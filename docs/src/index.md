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
CoxSym
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
Weyl.DescribeInvolution
Weyl.standard_parabolic
with_inversions
```
# PermRoot.jl Documentation
```@docs
PermRoot
reflection
hyperplane_orbits
```
# HeckeAlgebras.jl Documentation
```@docs
HeckeAlgebras
hecke
central_monomials
```
# KL.jl Documentation
```@docs
KL
KLPol
Tbasis
character
representation(c::LeftCell,H)
LeftCells
LeftCell
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
conjcat
endomorphisms
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
InductionTable
WGraphToRepresentation
```
# Cosets.jl Documentation
```@docs
Cosets
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
Ucl.InducedLinearForm
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
# SPerms.jl Documentation
```@docs
SPerms
SPerm
Perm(p::SPerm)
permuted(l::AbstractVector,a::SPerm)
Matrix
CoxHyperoctaedral
```
# Util.jl Documentation
```@docs
Util
groupby
constant
format
prime_residues
phi
primitiveroot
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
# GLinearAlgebra.jl Documentation
```@docs
GLinearAlgebra
GLinearAlgebra.echelon!
BigCellDecomposition
blocks
GLinearAlgebra.Transporter
```
# Eigenspaces.jl Documentation
```@docs
Eigenspaces
relative_degrees
regular_eigenvalues
```

# Dictionary from GAP3/Chevie
The dictionary from GAP3/Chevie is as follows:
```
CartanMat("A",5)                       → cartan(:A,5)
ChevieCharInfo                         → charinfo
ChevieClassInfo                        → classinfo
Coefficient(p,i)                       → p[i]
CoxeterElements(W[,l])                 → elements(W[,l])
CoxeterGroup("A",5)                    → coxgroup(:A,5)
CoxeterGroupHyperoctaedralGroup(n)     → CoxHyperoctaedral(n)
CoxeterGroupSymmetricGroup(n)          → CoxSym(n)
CoxeterLength(W,w)                     → length(W,w)
CoxeterWord(W,w)                       → word(W,w)
CyclotomicPolynomial(R,i)              → cyclotomic_polynomial(i)
Degree(p)                              → degree(p)
ElementWithInversions(W,l)             → with_inversions(W,l)
FiniteCoxeterTypeFromCartanMat(m)      → type_cartan(m)
FirstLeftDescending(W,w)               → firstleftdescent(W,w)
ForEachElement(W,f)                    → for w in W f(w) end
HeckeCentralMonomials                  → central_monomials
HyperplaneOrbits                       → hyperplane_orbits
IndependentRoots                       → independent_roots
Inversions                             → inversions
IsLeftDescending(W,w,i)                → isleftdescent(W,w,i)
LeadingCoefficient(p)                  → p[end]
LeftDescentSet(W,w)                    → leftdescents(W,w)
ListPerm(p)                            → vec(p)
LongestCoxeterElement(W)               → longest(W)
MatXPerm(W,p)                          → matX(W,p)
OnTuples(l,p)                          → l.^p
PermList(v)                            → Perm(v)
PermListList(l1,l2)                    → Perm(l1,l2)
Permuted(v,p)                          → permuted(v,p)
PositionClass                          → position_class
PrintDiagram(W)                        → Diagram(W)
ProportionalityCoefficient(v,w)        → ratio(v,w)
ReducedExpressions(W,w)                → reduced_words(W,w)
ReducedInRightCoset(W,w)               → reduced(W,w)
ReducedRightCosetRepresentatives(W,H)  → reduced(H,W)
Reflection                             → reflection
ReflectionCharacter                    → reflchar
ReflectionDegrees(W)                   → degrees(W)
ReflectionEigenvalues                  → refleigen
ReflectionLength(W,w)                  → reflength(W,w)
ReflectionSubgroup                     → reflection_subgroup
ReflectionType                         → refltype
Reflections                            → reflections
RestrictedPerm(p,d)                    → restricted(p,d)
RightDescentSet(W,w)                   → rightdescents(W,w)
RootsCartan(m)                         → roots(m)
SemiSimpleRank(W)                      → coxrank(W)
Size(W)                                → length(W)
StandardParabolic                      → standard_parabolic
TwoTree(m)                             → twotree(m)
Valuation(p)                           → valuation(p)
Value(p,x)                             → p(x)
W.N                                    → nref(W)
W.matgens[i]                           → matX(W,i)
W.orbitRepresentativeElement           → simple_conjugating_element(W,i)
W.orbitRepresentative[i]               → simple_representative(W,i)
```
