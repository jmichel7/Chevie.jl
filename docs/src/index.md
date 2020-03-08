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
Perms.order
Perms.orbit
Perms.orbits
cycles
cycletype
sign
Perms.restricted
reflength
```
# Groups.jl Documentation
```@docs
Groups
orbit
orbits
elements(Group)
transversal
centralizer
Groups.word
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
Root1
```
# Pols.jl Documentation
```@docs
Pols
divrem
divrem1
gcd
cyclotomic_polynomial
```
# Mvps.jl Documentation
```@docs
Mvps
Mvp
variables
Mvps.coefficients
Mvps.valuation
Mvps.degree
```
# CoxGroups.jl Documentation
```@docs
CoxGroups
firstleftdescent
reduced
word(W::CoxeterGroup,w)
bruhatless
CoxSym
longest
nref
braid_relations
coxmat
```
# Weyl.jl Documentation
```@docs
Weyl
cartan
two_tree
reflection_subgroup
coxgroup
rootdatum
describe_involution
Weyl.standard_parabolic
with_inversions
torus
SubTorus
fundamental_group
```
# PermRoot.jl Documentation
```@docs
PermRoot
reflection
hyperplane_orbits
bipartite_decomposition
catalan
matX
reflections
simple_conjugating_element
simple_representatives
invariant_form
```
# HeckeAlgebras.jl Documentation
```@docs
HeckeAlgebras
hecke
central_monomials
class_polynomials
char_values
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
leftgcd
α
DualBraidMonoid
fraction
word(b::Garside.GarsideElm)
image
representative_operation
centralizer_generators
conjcat
endomorphisms
shrink
```
# Chars.jl Documentation
```@docs
Chars
CharTable
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
UniChar
DLChar
AlmostChar
DLLefschetz
```
# Ucl.jl Documentation
```@docs
Ucl
UnipotentClasses
ICCTable
induced_linear_form
```
# HasType.jl Documentation
```@docs
ComplexReflectionGroup
Drinfeld_double
Family
family_imprimitive
schur_elements
traces_words_mats
```
# Symbols.jl Documentation
```@docs
Symbols
shiftβ
βset
partβ
ranksymbol
fegsymbol
degree_feg_symbol
degree_gendeg_symbol
valuation_feg_symbol
valuation_gendeg_symbol
tableaux
```
# SPerms.jl Documentation
```@docs
SPerms
SPerm
Perm(p::SPerm)
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
gcd_repr
```
# Combinat.jl Documentation
```@docs
arrangements
combinations
compositions
conjugate_partition
dominates
partitions
submultisets
```
# ModuleElts.jl Documentation
```@docs
ModuleElts
```
# Cycpols.jl Documentation
```@docs
CycPols
CycPol
```
# Posets.jl Documentation
```@docs
Posets
Poset
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
bigcell_decomposition
diagblocks
ratio
exterior_power
schur_functor
GLinearAlgebra.Transporter
```
# Eigenspaces.jl Documentation
```@docs
Eigenspaces
relative_degrees
regular_eigenvalues
eigenspace_projector
position_regular_class
split_levis
```

# Dictionary from GAP3/Chevie
The dictionary from GAP3/Chevie is as follows:
```
AlgebraicCentre                        → algebraic_centre
Arrangements                           → arrangements
AsFraction                             → fraction
AsReflection                           → reflection
AsRootOfUnity                          → Root1
AsWord                                 → word
AssociatedPartition                    → conjugate_partition
BetaSet                                → βset
BigCellDecomposition                   → bigcell_decomposition
BipartiteDecomposition                 → bipartite_decomposition
Braid                                  → BraidMonoid
BraidMonoid                            → BraidMonoid
BraidRelations                         → braid_relations
Bruhat                                 → bruhatless
BruhatPoset                            → Poset
BruhatSmaller                          → bruhatless
CartanMat("A",5)                       → cartan(:A,5)
CartanMatFromCoxeterMatrix             → cartan
Catalan                                → catalan
CentralizerGenerators                  → centralizer_generators
CharNames                              → charnames
CharRepresentationWords                → traces_words_mats
ChevieCharInfo                         → charinfo
ChevieClassInfo                        → classinfo
Coefficient(p,i)                       → p[i]
Combinations                           → combinations
ComplexConjugate                       → conj
ComplexReflectionGroup                 → ComplexReflectionGroup
Compositions                           → compositions
ConjugacySet(b[,F][,type])             → conjcat(b[,type[,F]]).obj
ConjugatePartition                     → conjugate_partition
CoxeterElements(W[,l])                 → elements(W[,l])
CoxeterGroup("A",5)                    → coxgroup(:A,5)
CoxeterGroupHyperoctaedralGroup(n)     → CoxHyperoctaedral(n)
CoxeterGroupSymmetricGroup(n)          → CoxSym(n)
CoxeterLength(W,w)                     → length(W,w)
CoxeterMatrix                          → coxmat
CoxeterMatrixFromCartanMat             → coxmat
CoxeterWord(W,w)                       → word(W,w)
CoxeterWords(W[,l])                    → word.(Ref(W),elements(W[,l])
CycPol                                 → CycPol
CycPolFakeDegreeSymbol                 → fegsymbol
CyclotomicPolynomial(R,i)              → cyclotomic_polynomial(i)
DecomposedMat                          → diagblocks
DefectSymbol                           → defectsymbol
Degree(p)                              → degree(p)
DescribeInvolution                     → describe_involution
Dominates                              → dominates
DualBraid                              → DualBraidMonoid
DualBraidMonoid                        → DualBraidMonoid
EigenspaceProjector                    → eigenspace_projector
ElementWithInversions(W,l)             → with_inversions(W,l)
EltBraid                               → image
EltWord(W,w)                           → W(w...)
ExteriorPower                          → exterior_power
FakeDegree                             → fakedegree
FakeDegrees                            → fakedegrees
FamiliesClassical                      → FamiliesClassical
Family                                 → Family
FamilyImprimitive                      → family_imprimitive
FiniteCoxeterTypeFromCartanMat(m)      → type_cartan(m)
FirstLeftDescending(W,w)               → firstleftdescent(W,w)
ForEachCoxeterWord(W,f)                → for w in W f(word(W,w)) end
ForEachElement(W,f)                    → for w in W f(w) end
GarsideAlpha                           → α
GarsideWords                           → elements
GcdPartitions                          → gcd_partitions
Hasse                                  → hasse
HeckeCentralMonomials                  → central_monomials
HeckeCharValues                        → char_values
HeckeClassPolynomials                  → class_polynomials
HighestPowerFakeDegreeSymbol           → degree_feg_symbol
HighestPowerGenericDegreeSymbol        → degree_gendeg_symbol
HyperplaneOrbits                       → hyperplane_orbits
ICCTable                               → ICCTable
Incidence                              → incidence
IndependentRoots                       → independent_roots
InducedLinearForm                      → induced_linear_form
InductionTable                         → InductionTable
IntListToString                        → joindigits
Inversions                             → inversions
IsCycPol(p)                            → p isa CycPol
IsFamily(f)                            → f isa Family
IsJoinLattice                          → is_join_lattice
IsLeftDescending(W,w,i)                → isleftdescent(W,w,i)
IsMeetLattice                          → is_meet_lattice
Join                                   → join
KazhdanLusztigPolynomial               → KLPol
LcmPartitions                          → lcm_partitions
LeadingCoefficient(p)                  → p[end]
LeftCell                               → LeftCell
LeftCells                              → LeftCells
LeftDescentSet(W,w)                    → leftdescents(W,w)
LeftDivisorsSimple                     → left_divisors
LeftGcd                                → leftgcd
LinearExtension                        → linear_extension
ListPerm(p)                            → vec(p)
LongestCoxeterElement(W)               → longest(W)
LongestCoxeterWord(W)                  → word(W,longest(W))
LowestPowerFakeDegreeSymbol            → valuation_feg_symbol
LowestPowerGenericDegreeSymbol         → valuation_gendeg_symbol
MatXPerm(W,p)                          → matX(W,p)
OnTuples(l,p)                          → l.^p
ParabolicRepresentatives               → parabolic_representatives
PartBeta                               → partβ
Partition                              → partition
PartitionTuples                        → partition_tuples
Partitions                             → partitions
PermList(v)                            → Perm(v)
PermListList(l1,l2)                    → Perm(l1,l2)
PermutationMat(p,dim)                  → Matrix(p,dim)
Permuted(v,p)                          → v^p
Poset                                  → Poset
PositionClass                          → position_class
PositionRegularClass                   → position_regular_class
PrintDiagram(W)                        → Diagram(W)
ProportionalityCoefficient(v,w)        → ratio(v,w)
Rank                                   → rank
RankSymbol                             → ranksymbol
ReducedCoxeterWord(W,w)                → word(W,W(w...))
ReducedExpressions(W,w)                → words(W,w)
ReducedInRightCoset(W,w)               → reduced(W,w)
ReducedRightCosetRepresentatives(W,H)  → reduced(H,W)
Reflection                             → reflection
ReflectionCharacter                    → reflchar
ReflectionCoDegrees(W)                 → codegrees(W)
ReflectionDegrees(W)                   → degrees(W)
ReflectionEigenvalues                  → refleigen
ReflectionLength(W,w)                  → reflength(W,w)
ReflectionSubgroup                     → reflection_subgroup
ReflectionType                         → refltype
Reflections                            → reflections
RegularEigenvalues                     → regular_eigenvalues
RelativeDegrees                        → relative_degrees
Representations                        → representations
RepresentativeConjugation              → representative_operation
Restricted                             → restricted
RestrictedPerm(p,d)                    → restricted(p,d)
Reversed                               → reverse
RightDescentSet(W,w)                   → rightdescents(W,w)
RightGcd                               → rightgcd
RightLcm                               → rightlcm
RootDatum                              → rootdatum
RootsCartan(m)                         → roots(m)
Rotation(v,i)                          → circshift(v,-i)
Rotations(v)                           → circshift.(Ref(a),length(a):-1:1)
SchurElements                          → schur_elements
SemisimpleRank                         → semisimplerank
SemisimpleRank(W)                      → coxrank(W)
ShiftBeta                              → shiftβ
ShrinkGarsideGeneratingSet             → shrink
SignedPerm                             → SPerm
SignedPermListList                     → SPerm
Size(W)                                → length(W)
StandardParabolic                      → standard_parabolic
StandardParabolicClass                 → standard_parabolic_class
SubTorus                               → SubTorus
Tableaux                               → tableaux
Torus                                  → torus
TransitiveClosure                      → transitive_closure
TwoTree(m)                             → twotree(m)
UnipotentCharacters                    → UnipotentCharacters
UnipotentClasses                       → UnipotentClasses
Valuation(p)                           → valuation(p)
Value(p,x)                             → p(x)
W.N                                    → nref(W)
W.matgens[i]                           → matX(W,i)
W.orbitRepresentativeElement           → simple_conjugating_element(W,i)
W.orbitRepresentative[i]               → simple_representative(W,i)
```
