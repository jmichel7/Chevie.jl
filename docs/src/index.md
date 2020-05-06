# Gapjm Documentation
```@docs
Gapjm
```
```@contents
```
# Permutations
```@docs
Perms
Perm
Perm(::Integer...)
Perm(::AbstractVector,::AbstractVector)
largest_moved_point
smallest_moved_point
Base.:^(::AbstractVector,::Perm) 
orbit(::Perm,::Integer,::Any)
orbits(::Perm,::Any)
order(::Perm)
cycles
cycletype
sign
Base.Matrix(::Perm,n)
Base.:^(::AbstractMatrix,::Perm)
restricted(::Perm,::AbstractVector{<:Integer})
reflength(::Perm)
mappingPerm
```
# Groups
```@docs
Groups
orbit(::Vector,::Any)
orbits(::Group,::AbstractVector)
elements(::Group)
transversal
centralizer
stabilizer
word(::Group,w)
length(::Group)
class_reps(::Group)
minimal_words
transporting_elt
```
# Permutation groups
```@docs
PermGroups
base
centralizers
transversals
symmetric_group
stab_onmat
perm_onmat
perm_rowcolmat
```
# Cyclotomic numbers
```@docs
Cycs
E
galois
ER
Quadratic
Root1
```
# Univariate (Laurent) polynomials
```@docs
Pols
divrem
divrem1
gcd
cyclotomic_polynomial
```
# Multivariate (Puiseux) polynomials
```@docs
Mvps
variables
Mvps.coefficients
Mvps.valuation
Mvps.value
Mvps.degree
Mvps.conj
Mvps.factor
derivative
laurent_denominator
scal
```
# Coxeter groups
```@docs
CoxGroups
isleftdescent
firstleftdescent
leftdescents
reduced
word(::CoxeterGroup,w)
length(::CoxeterGroup,w)
elements(::CoxeterGroup)
CoxGroups.words
bruhatless
CoxSym
reflection_subgroup(::CoxSym,::AbstractVector{Int})
longest
nref
braid_relations
coxmat
standard_parabolic_class
gencox
```
# Finite Coxeter groups and Weyl groups
```@docs
Weyl
cartan
two_tree
reflection_subgroup
coxgroup
rootdatum
describe_involution
Weyl.standard_parabolic
inversions
with_inversions
torus
SubTorus
fundamental_group
relative_group
```
# Finite reflection groups
```@docs
PermRoot
reflection
cartan(::PermRootGroup)
Diagram
hyperplane_orbits
rank
semisimplerank
degrees(::Group)
codegrees
bipartite_decomposition
catalan
refrep(::PRG,w)
PermX
reflections
refleigen
reflchar
simple_conjugating_element
simple_representatives
invariant_form
invariants
generic_order
torus_order
parabolic_representatives
ComplexReflectionGroup
```
# Hecke algebras
```@docs
HeckeAlgebras
hecke
central_monomials
class_polynomials
char_values
schur_elements
FactorizedSchurElement
FactorizedSchurElements
isrepresentation
refrep(::HeckeAlgebra)
```
# Kazhdan-Lusztig polynomials and bases
```@docs
KL
KLPol
Tbasis
KL.getCp
character
representation(::LeftCell,H)
LeftCells
LeftCell
```
# Garside monoids and groups, braids.
```@docs
Garside
left_divisors
leftgcd
α
DualBraidMonoid
fraction
word(::Garside.GarsideMonoid,w)
word(::Garside.GarsideElm)
elements(::Garside.LocallyGarsideMonoid,l)
image
conjugating_elt
centralizer_generators
conjcat
endomorphisms
shrink
```
# Classes/characters of reflection groups
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
# Reflection cosets
```@docs
Cosets
degrees(::Spets)
```
# Unipotent characters
```@docs
Uch
UnipotentCharacters
degrees(::UnipotentCharacters,q)
Uch.CycPolUnipotentDegrees
UniChar
DLChar
AlmostChar
DLLefschetz
LusztigInduce
LusztigRestrict
LusztigInductionTable
Families
Family
drinfeld_double
ndrinfeld_double
family_imprimitive
FamiliesClassical
```
# Unipotent classes of reductive groups
```@docs
Ucl
UnipotentClasses
ICCTable
induced_linear_form
```
# Symbols
```@docs
Symbols
shiftβ
βset
partβ
ranksymbol
defectsymbol
fegsymbol
gendeg_symbol
degree_feg_symbol
degree_gendeg_symbol
valuation_feg_symbol
valuation_gendeg_symbol
tableaux
symbols
```
# Signed permutations
```@docs
SPerms
SPerm
Perm(::SPerm)
orbit(::SPerm,::Integer)
order(::SPerm)
Matrix
CoxHyperoctaedral
reflection_subgroup(::CoxHyperoctaedral,::AbstractVector{Int})
stab_onsmat
perm_onsmat
```
# Utilities
```@docs
Util
groupby
tally
collectby
constant
format
prime_residues
phi
primitiveroot
gcd_repr
cut
```
# Combinatorics
```@docs
arrangements
combinations
compositions
conjugate_partition
dominates
partitions
submultisets
```
# Module Elements
```@docs
ModuleElts
```
# Cyclotomic polynomials
```@docs
CycPols
CycPol
```
# Posets
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
restricted(::Poset,::AbstractVector{<:Integer})
is_join_lattice
is_meet_lattice
Poset(::CoxeterGroup,w=longest(W))
```
# Linear algebra on any field/ring
```@docs
GLinearAlgebra
GLinearAlgebra.echelon!
bigcell_decomposition
diagblocks
blocks
ratio
exterior_power
permanent
symmetric_power
schur_functor
transporter
diagconj_elt
traces_words_mats
```
# Eigenspaces
```@docs
Eigenspaces
relative_degrees
regular_eigenvalues
eigenspace_projector
position_regular_class
split_levis
cuspidal_unipotent_characters
```

# Dictionary from GAP3/Chevie
The dictionary from GAP3/Chevie is as follows:
```
AlgebraicCentre                             → algebraic_centre
AlmostCharacter                             → AlmostChar
Arrangements                                → arrangements
AsFraction                                  → fraction
AsReflection                                → reflection
AsRootOfUnity                               → Root1
AsWord                                      → word
AssociatedPartition                         → conjugate_partition
BetaSet                                     → βset
BigCellDecomposition                        → bigcell_decomposition
Binomial                                    → binomial
BipartiteDecomposition                      → bipartite_decomposition
BlocksMat                                   → blocks
Braid                                       → BraidMonoid
BraidMonoid                                 → BraidMonoid
BraidRelations                              → braid_relations
Bruhat                                      → bruhatless
BruhatPoset                                 → Poset
BruhatSmaller                               → bruhatless
CartanMat("A",5)                            → cartan(:A,5)
CartanMatFromCoxeterMatrix                  → cartan
Catalan                                     → catalan
CentralizerGenerators                       → centralizer_generators
CharNames                                   → charnames
CharParams(W)                               → charinfo(W)[:charparams]
CharRepresentationWords                     → traces_words_mats
CheckHeckeDefiningRelations                 → isrepresentation
ChevieCharInfo                              → charinfo
ChevieClassInfo                             → classinfo
Coefficient(p,i)                            → p[i]
CollectBy                                   → collectby
Collected                                   → tally
Combinations                                → combinations
ComplexConjugate                            → conj
ComplexReflectionGroup                      → ComplexReflectionGroup
Compositions                                → compositions
ConcatenationString(s...)                   → prod(s)
ConjugacySet(b[,F][,type])                  → conjcat(b[,F],ss=type).obj
ConjugatePartition                          → conjugate_partition
CoxeterCoset                                → spets
CoxeterElements(W[,l])                      → elements(W[,l])
CoxeterGroup("A",5)                         → coxgroup(:A,5)
CoxeterGroupByCartanMatrix(C)               → gencox(C)
CoxeterGroupByCoxeterMatrix                 → gencox(cartan(C))
CoxeterGroupHyperoctaedralGroup(n)          → CoxHyperoctaedral(n)
CoxeterGroupSymmetricGroup(n)               → CoxSym(n)
CoxeterLength(W,w)                          → length(W,w)
CoxeterMatrix                               → coxmat
CoxeterMatrixFromCartanMat                  → coxmat
CoxeterSubCoset                             → subspets
CoxeterWord(W,w)                            → word(W,w)
CoxeterWords(W[,l])                         → word.(Ref(W),elements(W[,l]))
CuspidalUnipotentCharacters                 → cuspidal_unipotent_characters
CycPol                                      → CycPol
CycPolFakeDegreeSymbol                      → fegsymbol
CycPolGenericDegreeSymbol                   → gendeg_symbol
CycPolUnipotentDegrees                      → CycPolUnipotentDegrees
Cycle                                       → orbit
Cycles                                      → orbits
CyclotomicPolynomial(R,i)                   → cyclotomic_polynomial(i)
DecomposedMat                               → diagblocks
DefectSymbol                                → defectsymbol
Degree(p)                                   → degree(p)
DeligneLusztigCharacter                     → DLChar
DeligneLusztigLefschetz                     → DLLeftschetz
DescribeInvolution                          → describe_involution
Digits                                      → digits
Dominates                                   → dominates
DrinfeldDouble                              → drinfeld_double
Drop                                        → deleteat!
DualBraid                                   → DualBraidMonoid
DualBraidMonoid                             → DualBraidMonoid
EigenspaceProjector                         → eigenspace_projector
ElementWithInversions(W,l)                  → with_inversions(W,l)
Elements                                    → elements
EltBraid                                    → image
EltWord(W,w)                                → W(w...)
ExteriorPower                               → exterior_power
FakeDegree                                  → fakedegree
FakeDegrees                                 → fakedegrees
FamiliesClassical                           → FamiliesClassical
Family                                      → Family
FamilyImprimitive                           → family_imprimitive
FiniteCoxeterTypeFromCartanMat(m)           → type_cartan(m)
FirstLeftDescending(W,w)                    → firstleftdescent(W,w)
ForEachCoxeterWord(W,f)                     → for w in W f(word(W,w)) end
ForEachElement(W,f)                         → for w in W f(w) end
FullSymbol                                  → fullsymbol
FundamentalGroup                            → fundamental_group
GaloisCyc                                   → galois
GarsideAlpha                                → α
GarsideWords                                → elements
GcdPartitions                               → gcd_partitions
GcdRepresentation                           → gcd_repr
GenericOrder                                → generic_order
GenericSign                                 → generic_sign
GetRoot                                     → root
Hasse                                       → hasse
HeckeCentralMonomials                       → central_monomials
HeckeCharValues                             → char_values
HeckeClassPolynomials                       → class_polynomials
HeckeReflectionRepresentation               → refrep
HighestPowerFakeDegreeSymbol                → degree_feg_symbol
HighestPowerGenericDegreeSymbol             → degree_gendeg_symbol
HyperplaneOrbits                            → hyperplane_orbits
ICCTable                                    → ICCTable
Incidence                                   → incidence
IndependentLines(M)                         → echelon(M)[2]
IndependentRoots                            → independent_roots
InducedLinearForm                           → induced_linear_form
InductionTable                              → InductionTable
IntListToString                             → joindigits
Intersection                                → intersect
InvariantForm                               → invariant_form
Inversions                                  → inversions
IsAbelian                                   → isabelian
IsCycPol(p)                                 → p isa CycPol
IsFamily(f)                                 → f isa Family
IsJoinLattice                               → is_join_lattice
IsLeftDescending(W,w,i)                     → isleftdescent(W,w,i)
IsMeetLattice                               → is_meet_lattice
IsSubset(a,b)                               → issubset(b,a)
Join                                        → join
KazhdanLusztigPolynomial                    → KLPol
KroneckerProduct                            → kron
LcmPartitions                               → lcm_partitions
LeadingCoefficient(p)                       → p[end]
LeftCell                                    → LeftCell
LeftCells                                   → LeftCells
LeftDescentSet(W,w)                         → leftdescents(W,w)
LeftDivisorsSimple                          → left_divisors
LeftGcd                                     → leftgcd
LinearExtension                             → linear_extension
ListPerm(p)                                 → vec(p)
LongestCoxeterElement(W)                    → longest(W)
LongestCoxeterWord(W)                       → word(W,longest(W))
LowestPowerFakeDegreeSymbol                 → valuation_feg_symbol
LowestPowerGenericDegreeSymbol              → valuation_gendeg_symbol
LusztigInduction                            → LusztigInduce
LusztigInductionTable                       → LusztigInductionTable
LusztigRestriction                          → LusztigRestrict
MappingPermListList                         → mappingPerm
MatStab                                     → stab_onmat
MatXPerm(W,p)                               → refrep(W,p)
NrDrinfeldDouble                            → ndrinfeld_double
NrPartitionTuples                           → npartition_tuples
NrPartitions                                → npartitions
OnFamily(f,p::Int)                          → galois(f,p)
OnFamily(f,p::Perm)                         → f^p
OnMatrices(m,p)                             → ^(m,p;dims=(1,2))
OnPolynomials(m,p)                          → p^m
OnSets(s,g)                                 → unique!(sort(s.^g))
OnTuples(l,p)                               → l.^p
ParabolicRepresentatives                    → parabolic_representatives
PartBeta                                    → partβ
Partition                                   → partition
PartitionTuples                             → partition_tuples
Partitions                                  → partitions
PermList(v)                                 → Perm(v)
PermListList(l1,l2)                         → Perm(l1,l2)
PermMatMat                                  → perm_onmat
PermMatX                                    → PermX
PermutationMat(p,dim)                       → Matrix(p,dim)
Permuted(v,p)                               → v^p
PermutedByCols(m,p)                         → ^(m,p;dims=2)
Poset                                       → Poset
PositionClass                               → position_class
PositionRegularClass                        → position_regular_class
PrintDiagram(W)                             → Diagram(W)
ProportionalityCoefficient(v,w)             → ratio(v,w)
Rank                                        → rank
RankSymbol                                  → ranksymbol
ReducedCoxeterWord(W,w)                     → word(W,W(w...))
ReducedExpressions(W,w)                     → words(W,w)
ReducedInRightCoset(W,w)                    → reduced(W,w)
ReducedRightCosetRepresentatives(W,H)       → reduced(H,W)
Reflection                                  → reflection
ReflectionCharacter                         → reflchar
ReflectionCoDegrees(W)                      → codegrees(W)
ReflectionDegrees(W)                        → degrees(W)
ReflectionEigenvalues                       → refleigen
ReflectionLength(W,w)                       → reflength(W,w)
ReflectionSubgroup                          → reflection_subgroup
ReflectionType                              → refltype
Reflections                                 → reflections
RegularEigenvalues                          → regular_eigenvalues
RelativeDegrees                             → relative_degrees
Representations                             → representations
RepresentativeConjugation(b,b'[,F][,type])  → conjugating_elt(b,b'[,F],ss=type)
RepresentativeDiagonalConjugation           → diagconj_elt
RepresentativeOperation                     → transporting_elt
RepresentativeRowColPermutation             → perm_rowcolmat
Restricted                                  → restricted
RestrictedPerm(p,d)                         → restricted(p,d)
Reversed                                    → reverse
RightDescentSet(W,w)                        → rightdescents(W,w)
RightGcd                                    → rightgcd
RightLcm                                    → rightlcm
RootDatum                                   → rootdatum
RootsCartan(m)                              → roots(m)
Rotation(v,i)                               → circshift(v,-i)
Rotations(v)                                → circshift.(Ref(a),length(a):-1:1)
ScalMvp                                     → scal
SchurElements                               → schur_elements
SchurFunctor                                → schur_functor
SemisimpleRank                              → semisimplerank
SemisimpleRank(W)                           → coxrank(W)
ShiftBeta                                   → shiftβ
ShrinkGarsideGeneratingSet                  → shrink
SignedMatStab                               → stab_onsmat
SignedPerm                                  → SPerm
SignedPermListList                          → SPerm
SignedPermMatMat                            → perm_onsmat
Size(W)                                     → length(W)
SolutionMat                                 → solutionmat
Spets                                       → spets
SplitLevis                                  → split_levis
StandardParabolic                           → standard_parabolic
StandardParabolicClass                      → standard_parabolic_class
SubSpets                                    → subspets
SubTorus                                    → SubTorus
SymmetricDifference                         → symdiff
SymmetricPower                              → symmetric_power
Tableaux                                    → tableaux
Torus                                       → torus
TransitiveClosure                           → transitive_closure
Transporter                                 → transporter
TriangulizeMat                              → echelon!
Twistings                                   → twistings
TwoTree(m)                                  → twotree(m)
UnipotentCharacter                          → UniChar
UnipotentCharacters                         → UnipotentCharacters
UnipotentClasses                            → UnipotentClasses
UnipotentDegrees(W,q)                       → degrees(UnipotentCharacters(W),q)
UnorderedTuples                             → submultisets
Valuation(p)                                → valuation(p)
Value(p,x)                                  → p(x)
W.N                                         → nref(W)
W.matgens[i]                                → refrep(W,i)
W.orbitRepresentativeElement                → simple_conjugating_element(W,i)
W.orbitRepresentative[i]                    → simple_representative(W,i)
WeightInfo                                  → weightinfo
```
