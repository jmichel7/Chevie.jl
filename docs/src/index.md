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
@perm_str
largest_moved_point
smallest_moved_point
Base.:^(::AbstractVector,::Perm) 
sortPerm
Perms.orbit
Perms.orbits
Perms.order
cycles
cycletype
support
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
classreps(::Group)
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
stab_onmats
Perm_onmats
Perm_rowcolmat
```
# Cyclotomic numbers
```@docs
Cycs
conductor
coefficients(c::Cyc)
E
galois
ER
Quadratic
Root1
Cycs.root
```
# Univariate (Laurent) polynomials
```@docs
Pols
divrem
gcd
cyclotomic_polynomial
```
# Multivariate (Puiseux) polynomials
```@docs
Mvps
Mvp
variables
Mvps.coefficients(::Mvp,::Symbol)
Mvps.valuation
Mvps.value
Mvps.degree
Mvps.conj
factor(::Mvp)
derivative
laurent_denominator
scal
```
# Cyclotomic polynomials
```@docs
CycPols
CycPol
```
# Utilities
```@docs
Util
@forward
groupby
tally
collectby
constant
format
cut
prime_residues
phi
primitiveroot
```
# Combinatorics
```@docs
arrangements
combinations
compositions
conjugate_partition
dominates
partitions
npartitions
partition_tuples
npartition_tuples
submultisets
partitions_set
npartitions_set
bell
stirling2
```
# Module Elements
```@docs
ModuleElts
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
Posets.restricted(::Poset,::AbstractVector{<:Integer})
is_join_lattice
is_meet_lattice
Poset(::CoxeterGroup,w=longest(W))
```
# Signed permutations
```@docs
SPerms
SPerm
Perm(::SPerm)
@sperm_str
orbit(::SPerm,::Integer)
order(::SPerm)
Matrix
CoxHyperoctaedral
reflection_subgroup(::CoxHyperoctaedral,::AbstractVector{Int})
sstab_onmats
SPerm_onmats
```
# Linear algebra on any field/ring
```@docs
GLinearAlgebra
GLinearAlgebra.echelon!
GLinearAlgebra.echelon
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
solutionmat
```
# Integral matrices and lattices
```@docs
DiaconisGraham
```
# Finite fields
```@docs
FFields
mod(::Cyc,::Any)
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
GenCox
```
# Finite Coxeter groups and Weyl groups
```@docs
Weyl
cartan
two_tree
roots
reflection_subgroup
coxgroup
rootdatum
describe_involution
Weyl.standard_parabolic
inversions
with_inversions
torus
SubTorus
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
reflrep(::PermRootGroup,w)
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
alt
α(::HeckeTElt)
schur_elements
FactorizedSchurElement
FactorizedSchurElements
isrepresentation
reflrep(::HeckeAlgebra)
```
# Kazhdan-Lusztig polynomials and bases
```@docs
KL
KLPol
Tbasis
KL.getCp
Cbasis
Cpbasis
character
representation(::LeftCell,H)
LeftCells
LeftCell
Lusztigaw
LusztigAw
AsymptoticAlgebra
```
# Garside monoids and groups, braids.
```@docs
Garside
left_divisors
leftgcd
rightgcd
rightlcm
α(::Garside.GarsideElm)
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
charnames
classinfo
fakedegree
fakedegrees
representation
representations
InductionTable
jInductionTable
JInductionTable
detPerm
WGraphToRepresentation
pblocks
```
# Reductive groups, semisimple elements
```@docs
Semisimple
fundamental_group
QuasiIsolatedRepresentatives
is_isolated
torsion_subgroup
algebraic_centre
SScentralizer_representatives
StructureRationalPointsConnectedCentre
```
# Reflection cosets
```@docs
Cosets
degrees(::Spets)
spets
twistings
subspets
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
fourier
drinfeld_double
ndrinfeld_double
family_imprimitive
FamiliesClassical
*(f::Family, g::Family)
fusion_algebra
```
# Unipotent classes of reductive groups
```@docs
Ucl
UnipotentClasses
ICCTable
induced_linear_form
special_pieces
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
# Eigenspaces
```@docs
Eigenspaces
relative_degrees
regular_eigenvalues
eigenspace_projector
position_regular_class
split_levis
cuspidal
```
# Classtypes
```@docs
ClassTypes
```
# Unipotent Elements
```@docs
Urad
UnipotentGroup
Urad.norm
Urad.abelianpart
Urad.decompose
```
# Dictionary from GAP3/Chevie
The dictionary from GAP3/Chevie is as follows:
```
AlgebraicCentre                             algebraic_centre
AlmostCharacter                             AlmostChar
Arrangements                                arrangements
AsFraction                                  fraction
AsReflection                                reflection
AsRootOfUnity                               Root1
AsWord                                      word
AssociatedPartition                         conjugate_partition
AsymptoticAlgebra                           AsymptoticAlgebra
BetaSet                                     βset
BigCellDecomposition                        bigcell_decomposition
Binomial                                    binomial
BipartiteDecomposition                      bipartite_decomposition
BlocksMat                                   blocks
Braid                                       BraidMonoid
BraidMonoid                                 BraidMonoid
BraidRelations                              braid_relations
Bruhat                                      bruhatless
BruhatPoset                                 Poset
BruhatSmaller                               bruhatless
CartanMat("A",5)                            cartan(:A,5)
CartanMatFromCoxeterMatrix                  cartan
Catalan                                     catalan
CentralizerGenerators                       centralizer_generators
CharFFE(x)                                  field(x).p
CharNames                                   charnames
CharParams(W)                               charinfo(W)[:charparams]
CharRepresentationWords                     traces_words_mats
CharTable                                   CharTable
CheckHeckeDefiningRelations                 isrepresentation
ChevieCharInfo                              charinfo
ChevieClassInfo                             classinfo
ClassTypes                                  ClassTypes
Coefficient(p,i)                            p[i]
CollectBy(l,f)                              collectby(f,l)
Collected                                   tally
Combinations                                combinations
ComplexConjugate                            conj
ComplexReflectionGroup                      ComplexReflectionGroup
Compositions                                compositions
Concatenation(s::Vector...)                 vcat(s...)
ConcatenationString(s...)                   prod(s)
ConjugacySet(b[,F][,type])                  conjcat(b[,F],ss=type).obj
ConjugatePartition                          conjugate_partition
CoxeterCoset                                spets
CoxeterElements(W[,l])                      elements(W[,l])
CoxeterGroup("A",5)                         coxgroup(:A,5)
CoxeterGroupByCartanMatrix(C)               gencox(C)
CoxeterGroupByCoxeterMatrix(C)              gencox(cartan(C))
CoxeterGroupHyperoctaedralGroup(n)          CoxHyperoctaedral(n)
CoxeterGroupSymmetricGroup(n)               CoxSym(n)
CoxeterLength(W,w)                          length(W,w)
CoxeterMatrix                               coxmat
CoxeterMatrixFromCartanMat                  coxmat
CoxeterSubCoset                             subspets
CoxeterWord(W,w)                            word(W,w)
CoxeterWords(W[,l])                         word.(Ref(W),elements(W[,l]))
CuspidalUnipotentCharacters(W)              cuspidal(UnipotentCharacters(W))
CycPol                                      CycPol
CycPolFakeDegreeSymbol                      fegsymbol
CycPolGenericDegreeSymbol                   gendeg_symbol
CycPolUnipotentDegrees                      CycPolUnipotentDegrees
Cycle                                       orbit
Cycles                                      orbits
CyclotomicPolynomial(R,i)                   cyclotomic_polynomial(i)
DecomposedMat                               diagblocks
DefectSymbol                                defectsymbol
Degree(p)                                   degree(p)
DegreeFFE(x)                                field(x).n
DeligneLusztigCharacter                     DLChar
DeligneLusztigLefschetz                     DLLeftschetz
DescribeInvolution                          describe_involution
DetPerm                                     detPerm
Digits                                      digits
Dominates                                   dominates
DrinfeldDouble                              drinfeld_double
Drop                                        deleteat!
DualBraid                                   DualBraidMonoid
DualBraidMonoid                             DualBraidMonoid
EigenspaceProjector                         eigenspace_projector
ElementWithInversions(W,l)                  with_inversions(W,l)
Elements                                    elements
EltBraid                                    image
EltWord(W,w)                                W(w...)
ExteriorPower                               exterior_power
FactorizedSchurElement                      FactorizedSchurElement
FactorizedSchurElements                     FactorizedSchurElements
FakeDegree                                  fakedegree
FakeDegrees                                 fakedegrees
FamiliesClassical                           FamiliesClassical
Family                                      Family
FamilyImprimitive                           family_imprimitive
FiniteCoxeterTypeFromCartanMat(m)           type_cartan(m)
FirstLeftDescending(W,w)                    firstleftdescent(W,w)
ForEachCoxeterWord(W,f)                     for w in W f(word(W,w)) end
ForEachElement(W,f)                         for w in W f(w) end
FullSymbol                                  fullsymbol
FundamentalGroup                            fundamental_group
FusionAlgebra                               fusion_algebra
FusionConjugacyClasses                      fusion_conjugacy_classes
GaloisCyc                                   galois
GarsideAlpha                                α
GarsideWords                                elements
GcdPartitions                               gcd_partitions
GcdRepresentation(x,y)                      gcdx(x,y)[2:3]
GenericOrder                                generic_order
GenericSign                                 generic_sign
GetRoot                                     root
Hasse                                       hasse
Hecke                                       hecke
HeckeCentralMonomials                       central_monomials
HeckeCharValues                             char_values
HeckeClassPolynomials                       class_polynomials
HeckeReflectionRepresentation               reflrep
HighestPowerFakeDegreeSymbol                degree_feg_symbol
HighestPowerFakeDegrees(W)                  charinfo(W)[:B]
HighestPowerGenericDegreeSymbol             degree_gendeg_symbol
HighestPowerGenericDegrees(W)               charinfo(W)[:A]
HyperplaneOrbits                            hyperplane_orbits
ICCTable                                    ICCTable
Incidence                                   incidence
IndependentLines(M)                         echelon(M)[2]
IndependentRoots                            independent_roots
InducedLinearForm                           induced_linear_form
InductionTable                              InductionTable
IntFFE                                      Int
IntListToString                             joindigits
Intersection                                intersect
InvariantForm                               invariant_form
Invariants                                  invariants
Inversions                                  inversions
IsAbelian                                   isabelian
IsCycPol(p)                                 p isa CycPol
IsFFE(x)                                    x isa FFE
IsFamily(f)                                 f isa Family
IsIsolated                                  is_isolated
IsJoinLattice                               is_join_lattice
IsLeftDescending(W,w,i)                     isleftdescent(W,w,i)
IsMeetLattice                               is_meet_lattice
IsSubset(a,b)                               issubset(b,a)
IsUnipotentElement(x)                       x isa UnipotentElement
IsomorphismType                             IsomorphismType
JInductionTable                             JInductionTable
Join                                        join
KazhdanLusztigPolynomial                    KLPol
KroneckerProduct                            kron
LargestMovedPoint                           largest_moved_point
LcmPartitions                               lcm_partitions
LeadingCoefficient(p)                       p[end]
LeftCell                                    LeftCell
LeftCells                                   LeftCells
LeftDescentSet(W,w)                         leftdescents(W,w)
LeftDivisorsSimple                          left_divisors
LeftGcd                                     leftgcd
LinearExtension                             linear_extension
List(ConjugacyClasses(G),Representative)    classreps(G)
ListBlist(a,b)                              a[b]
ListPerm(p)                                 vec(p)
LogFFE                                      log
LongestCoxeterElement(W)                    longest(W)
LongestCoxeterWord(W)                       word(W,longest(W))
LowestPowerFakeDegreeSymbol                 valuation_feg_symbol
LowestPowerFakeDegrees(W)                   charinfo(W)[:b]
LowestPowerGenericDegreeSymbol              valuation_gendeg_symbol
LowestPowerGenericDegrees(W)                charinfo(W)[:a]
LusztigAw                                   LusztigAw
LusztigInduction                            LusztigInduce
LusztigInductionTable                       LusztigInductionTable
LusztigRestriction                          LusztigRestrict
Lusztigaw                                   Lusztigaw
MappingPermListList                         mappingPerm
MatStab                                     stab_onmats
MatXPerm(W,p)                               reflrep(W,p)
MatYPerm                                    matY
MovedPoints                                 support
Mvp("x")                                    Mvp(:x)
NrArrangements                              narrangements
NrDrinfeldDouble                            ndrinfeld_double
NrPartitionTuples                           npartition_tuples
NrPartitions                                npartitions
OnFamily(f,p::Int)                          galois(f,p)
OnFamily(f,p::Perm)                         f^p
OnMatrices(m,p)                             ^(m,p;dims=(1,2))
OnPolynomials(m,p)                          p^m
OnSets(s,g)                                 unique!(sort(s.^g))
OnTuples(l,p)                               l.^p
OrderFFE                                    order
OrderMod(n,m)                               order(Mod{m}(n))
ParabolicRepresentatives                    parabolic_representatives
PartBeta                                    partβ
Partition                                   partition
PartitionTuples                             partition_tuples
Partitions                                  partitions
PermList(v)                                 Perm(v)
PermListList(l1,l2)                         Perm(l1,l2)
PermMatMat                                  Perm_onmats
PermMatX                                    PermX
PermutationMat(p,dim)                       Matrix(p,dim)
Permuted(v,p)                               v^p
PermutedByCols(m,p)                         ^(m,p;dims=2)
Poset                                       Poset
Position(l,x)                               findfirst(==(x),l)
PositionCartesian(a,b)                      LinearIndices(reverse(Tuple(a)))[CartesianIndices(Tuple(b))]
PositionClass                               position_class
PositionDet                                 charinfo(W)[:PositionDet]
PositionId                                  charinfo(W)[:PositionId]
PositionRegularClass                        position_regular_class
PrintDiagram(W)                             Diagram(W)
ProportionalityCoefficient(v,w)             ratio(v,w)
QuasiIsolatedRepresentatives                QuasiIsolatedRepresentatives
Rank                                        rank
RankSymbol                                  ranksymbol
ReducedCoxeterWord(W,w)                     word(W,W(w...))
ReducedExpressions(W,w)                     words(W,w)
ReducedInRightCoset(W,w)                    reduced(W,w)
ReducedRightCosetRepresentatives(W,H)       reduced(H,W)
Reflection                                  reflection
ReflectionCharacter                         reflchar
ReflectionCoDegrees(W)                      codegrees(W)
ReflectionDegrees(W)                        degrees(W)
ReflectionEigenvalues                       refleigen
ReflectionGroup                             reflection_group
ReflectionLength(W,w)                       reflength(W,w)
ReflectionSubgroup                          reflection_subgroup
ReflectionType                              refltype
Reflections                                 reflections
RegularEigenvalues                          regular_eigenvalues
RelativeDegrees                             relative_degrees
Representations                             representations
RepresentativeConjugation(b,b'[,F][,type])  conjugating_elt(b,b'[,F],ss=type)
RepresentativeDiagonalConjugation           diagconj_elt
RepresentativeOperation                     transporting_elt
RepresentativeRowColPermutation             Perm_rowcolmat
Restricted                                  restricted
RestrictedPerm(p,d)                         restricted(p,d)
Reversed                                    reverse
RightDescentSet(W,w)                        rightdescents(W,w)
RightGcd                                    rightgcd
RightLcm                                    rightlcm
RootDatum                                   rootdatum
RootsCartan(m)                              roots(m)
Rotation(v,i)                               circshift(v,-i)
Rotations(v)                                circshift.(Ref(a),length(a):-1:1)
ScalMvp                                     scal
SchurElements                               schur_elements
SchurFunctor                                schur_functor
SemisimpleCentralizerRepresentatives        SScentralizer_representatives
SemisimpleElement                           SS
SemisimpleRank                              semisimplerank
SemisimpleRank(W::CoxeterGroup)             coxrank(W)
SemisimpleSubgroup                          torsion_subgroup
ShiftBeta                                   shiftβ
ShrinkGarsideGeneratingSet                  shrink
SignedMatStab                               sstab_onmats
SignedPerm                                  SPerm
SignedPermListList                          SPerm
SignedPermMatMat                            SPerm_onmats
Size(W)                                     length(W)
SmallestMovedPoint                          smallest_moved_point
SolutionMat                                 solutionmat
Spets                                       spets
SplitLevis                                  split_levis
StandardParabolic                           standard_parabolic
StandardParabolicClass                      standard_parabolic_class
StructureRationalPointsConnectedCentre      StructureRationalPointsConnectedCentre
SubSpets                                    subspets
SubTorus                                    SubTorus
SymmetricDifference                         symdiff
SymmetricPower                              symmetric_power
Tableaux                                    tableaux
Torus                                       torus
TransitiveClosure                           transitive_closure
Transporter                                 transporter
TriangulizeMat                              echelon!
Twistings                                   twistings
TwoTree(m)                                  twotree(m)
UnipotentAbelianPart                        abelianpart
UnipotentCharacter                          UniChar
UnipotentCharacters                         UnipotentCharacters
UnipotentClasses                            UnipotentClasses
UnipotentDecompose                          decompose
UnipotentDegrees(W,q)                       degrees(UnipotentCharacters(W),q)
UnipotentGroup                              UnipotentGroup
UnorderedTuples                             submultisets
Valuation(p)                                valuation(p)
Value(p,x)                                  p(x)
W.N                                         nref(W)
W.matgens[i]                                reflrep(W,i)
W.orbitRepresentative                       simple_representatives(W)
W.orbitRepresentativeElements[i]            simple_conjugating_element(W,i)
W.rootLengths                               rootlengths(W)
WeightInfo                                  weightinfo
jInductionTable                             jInductionTable
last                                        ans
```
