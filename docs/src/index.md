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
largest_moved_point(::Perm)
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
Perm_rowcol
```
# Groups
```@docs
Groups
Group
gens
ngens
orbit(::Vector,::Any)
orbits(::Group,::AbstractVector)
elements(::Group)
transversal
centralizer
centre
stabilizer
normalizer
word(::Group,w)
comm
length(::Group)
classreps(::Group)
conjugacy_classes
conjugacy_class
nconjugacy_classes
position_class
fusion_conjugacy_classes
minimal_words
words(::Group)
transporting_elt
isabelian
iscyclic
abelian_gens
Hom
kernel
blocks(G,p)
```
# Permutation groups
```@docs
PermGroups
largest_moved_point(::PermGroup)
base
centralizers
transversals
in(::Perm,::PermGroup)
on_classes
symmetric_group
onmats
stab_onmats
Perm_onmats
```
# Extensions to Laurent and Puiseux polynomials
```@docs
FFfac.factor(::Pol{FFE{p}}where p, F)
Fact.factor(f::Pol{var"#s162"} where var"#s162"<:Union{Integer, Rational})
factor(p::Mvp{T, N}) where {T, N}
```
# Cyclotomic polynomials
```@docs
cyclotomic_polynomial
CycPols
CycPol
eigmat
```
# Utilities
```@docs
Util
@forward
@GapObj
showtable
cut
```
# Combinatorics
```@docs
Combinat
combinations
Combinat.Combinations
arrangements
partitions
Combinat.Partitions
partition_tuples
compositions
multisets
lcm_partitions
gcd_partitions
conjugate_partition
dominates
tableaux
robinson_schensted
bell
stirling1
stirling2
catalan(::Integer)
bernoulli
groupby
tally
collectby
cartesian
unique_sorted!
diagblocks
blocks(m::AbstractMatrix)
prime_residues
primitiveroot
divisors
```
# Posets
```@docs
Posets
Poset
hasse
incidence
transitive_closure
linear_extension
reverse
partition
covering_chains
Posets.restricted(::Poset,::AbstractVector{<:Integer})
is_join_lattice
is_meet_lattice
moebius
moebiusmatrix
minimum(P::Poset)
maximum(P::Poset)
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
ratio
exterior_power
charpoly
comatrix
permanent
symmetric_power
schur_functor
transporter
diagconj_elt
traces_words_mats
solutionmat
sum_rowspace
intersect_rowspace
lnullspace
GLinearAlgebra.det
GLinearAlgebra.rank
```
# Integral matrices and lattices
```@docs
smith
smith_transforms
hermite
hermite_transforms
col_hermite
col_hermite_transforms
baseInt
lnullspaceInt
complementInt
solutionmatInt
diaconis_graham
```
# Finite fields
```@docs
FFields
FFE
FFE(i::Integer)
Z(::Any)
```
# Presentations
```@docs
Presentations
AbsWord
@AbsWord
FpGroup
Presentation(::FpGroup)
relators
simplify
conjugate
tryconjugate
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
coxmat
standard_parabolic_class
GenCox
```
# Finite Coxeter groups and Weyl groups
```@docs
Weyl
cartan(::Symbol,::Integer,::Integer)
cartan(::AbstractMatrix)
roots(::AbstractMatrix)
two_tree
reflection_subgroup(::Weyl.FCG,::AbstractVector{<:Integer})
coxgroup
rootdatum
describe_involution
badprimes
Weyl.standard_parabolic
inversions
with_inversions
torus
istorus
SubTorus
relative_group
```
# Finite reflection groups
```@docs
PermRoot
reflection
cartan(::PermRootGroup)
cartan(::PermRootGroup,i,j)
Diagram
hyperplane_orbits
rank
semisimplerank
degrees(::Group)
codegrees
nref
nhyp
roots
coroots
coroot
simpleroots
simplecoroots
braid_relations
bipartite_decomposition
catalan(W,m)
reflrep(::PermRootGroup,w)
reflrep(::PRG)
reflrep(::PRG,::Integer)
PermX
reflections
refleigen
reflchar
simple_conjugating
simple_reps
invariant_form
invariants
generic_order
torus_order
parabolic_reps
parabolic_closure
is_parabolic
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
HeckeCoset
```
# Kazhdan-Lusztig polynomials and bases
```@docs
KL
KLPol
Tbasis
Cbasis
Cpbasis
character
representation(::LeftCell,H)
Wgraph
LeftCells
LeftCell
Lusztigaw
LusztigAw
AsymptoticAlgebra
```
# Garside monoids and groups, braids.
```@docs
Garside
LocallyGarsideMonoid
GarsideMonoid
left_divisors
leftgcd
rightgcd
leftlcm
rightlcm
α(::Garside.LocallyGarsideElt)
α(::Garside.LocallyGarsideElt,::AbstractVector)
Brieskorn_normal_form
BraidMonoid
DualBraidMonoid
hurwitz
fraction
word(::Garside.GarsideMonoid,w)
word(::Garside.GarsideElt)
elements(::Garside.LocallyGarsideMonoid,l)
image
conjugating_elt
centralizer_gens
conjcat
endomorphisms
Presentation(::GarsideMonoid)
shrink
```
# Classes/characters of reflection groups
```@docs
Chars
CharTable
on_chars
charinfo
charnames
classnames
classinfo
fakedegree
fakedegrees
representation
representations
InductionTable
jInductionTable
JInductionTable
detPerm
conjPerm
WGraphToRepresentation
```
# Reductive groups, semisimple elements
```@docs
Semisimple
fundamental_group
intermediate_group
QuasiIsolatedRepresentatives(::FiniteCoxeterGroup)
is_isolated
torsion_subgroup
algebraic_centre
weightinfo
weights
coweights
centralizer(::FiniteCoxeterGroup,::SemisimpleElement)
SScentralizer_reps
StructureRationalPointsConnectedCentre
```
# Reflection cosets
```@docs
Cosets
degrees(::Spets)
spets
twistings
graph_automorphisms
subspets
```
# Non-connected reductive groups
```@docs
Sscoset
centralizer(::Spets,::SemisimpleElement{Root1})
QuasiIsolatedRepresentatives(::Spets)
is_isolated(::Spets,::SemisimpleElement{Root1})
```
# Unipotent characters
```@docs
Uch
UnipotentCharacters
degrees(::UnipotentCharacters,q)
Uch.CycPolUnipotentDegrees
UniChar
DLChar
almostChar
on_unipotents
DLLefschetz
LusztigInduce
LusztigRestrict
LusztigInductionTable
Families
Family
galois(f::Family,p::Int)
fourier
drinfeld_double
ndrinfeld_double
family_imprimitive
FamiliesClassical
*(f::Family, g::Family)
fusion_algebra
```
# d-Harish-Chandra series
```@docs
dSeries
cuspidal_data
Series
ennola
```
# Unipotent classes of reductive groups
```@docs
Ucl
UnipotentClasses
ICCTable
XTable
GreenTable
UnipotentValues
induced_linear_form
special_pieces
distinguished_parabolics
```
# Symbols
```@docs
Symbols
shiftβ
βset
partβ
symbols
ranksymbol
defectsymbol
symbol_partition_tuple
fegsymbol
gendeg_symbol
degree_fegsymbol
degree_gendeg_symbol
valuation_fegsymbol
valuation_gendeg_symbol
XSP
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
closed_subsystems
ClassTypes
```
# Unipotent Elements
```@docs
Urad
UnipotentGroup
Urad.reorder
Urad.abelianpart
Urad.decompose
```
# Decomposition Matrices
```@docs
decomposition_matrix
generic_decomposition_matrix
InducedDecompositionMatrix
```
# Dictionary from GAP3/Chevie
The dictionary from GAP3/Chevie is as follows:
```
AbelianGenerators                           abelian_gens
AlgebraicCentre                             algebraic_centre
AlmostCharacter                             AlmostChar
Arrangements                                arrangements
AsFraction                                  fraction
AsReflection                                reflection
AsRootOfUnity                               Root1
AssociatedPartition                         conjugate_partition
AsWord                                      word
AsymptoticAlgebra                           AsymptoticAlgebra
BadPrimes                                   badprimes
BaseIntMat                                  baseInt
BetaSet                                     βset
BigCellDecomposition                        bigcell_decomposition
Binomial                                    binomial
BipartiteDecomposition                      bipartite_decomposition
BlocksMat                                   blocks
Braid                                       BraidMonoid
BraidMonoid                                 BraidMonoid
BraidRelations                              braid_relations
BrieskornNormalForm                         Brieskorn_normal_form
Bruhat                                      bruhatless
BruhatPoset                                 Poset
BruhatSmaller                               bruhatless
CartanMat("A",5)                            cartan(:A,5)
CartanMatFromCoxeterMatrix                  cartan
Cartesian                                   cartesian
CartesianAt                                 lin2cart
Catalan                                     catalan
CentralizerGenerators                       centralizer_gens
CharFFE(x)                                  field(x).p
CharNames                                   charnames
CharParams(W)                               charinfo(W)[:charparams]
CharRepresentationWords                     traces_words_mats
CharTable                                   CharTable
CheckHeckeDefiningRelations                 isrepresentation
ChevieCharInfo                              charinfo
ChevieClassInfo                             classinfo
ClassName                                   see ClassNames
ClassTypes                                  ClassTypes
Coefficient(p,i)                            p[i]
CollectBy(l,f)                              collectby(f,l)
Collected                                   tally
Combinations                                combinations
Comm                                        comm
ComplementIntMat                            complementInt
ComplexConjugate                            conj
ComplexReflectionGroup                      ComplexReflectionGroup
Compositions                                compositions
Concatenation(s::Vector...)                 vcat(s...)
ConcatenationString(s...)                   prod(s)
ConjugacyClasses                            conjugacy_classes
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
CuspidalPairs                               cuspidal_data
CuspidalUnipotentCharacters(W[,d])          cuspidal(UnipotentCharacters(W)[,d])
Cycle                                       orbit
Cycles                                      orbits
CyclotomicModP(c,p)                         FFE{p}(c)
CyclotomicPolynomial(R,i)                   cyclotomic_polynomial(i)
CycPol                                      CycPol
CycPolFakeDegreeSymbol                      fegsymbol
CycPolGenericDegreeSymbol                   gendeg_symbol
CycPolUnipotentDegrees                      CycPolUnipotentDegrees
DecomposedMat                               diagblocks
DefectSymbol                                defectsymbol
Degree(p)                                   degree(p)
DegreeFFE(x)                                field(x).n
DeligneLusztigCharacter                     DLChar
DeligneLusztigLefschetz                     DLLeftschetz
DescribeInvolution                          describe_involution
DetPerm(W)                                  vec(detPerm(W))
DiaconisGraham                              diaconis_graham
Digits                                      digits
DistinguishedParabolicSubgroups             distinguished_parabolics
Dominates                                   dominates
DrinfeldDouble                              drinfeld_double
Drop                                        deleteat!
DualBraid                                   DualBraidMonoid
DualBraidMonoid                             DualBraidMonoid
EigenspaceProjector                         eigenspace_projector
EigenvaluesMat                              eigmat
Elements                                    elements
ElementWithInversions(W,l)                  with_inversions(W,l)
EltBraid                                    image
EltWord(W,w)                                W(w...)
ER                                          root
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
FormatTable                                 showtable
Frobenius                                   Frobenius
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
GraphAutomorphisms                          graph_automorphisms
Hasse                                       hasse
Hecke                                       hecke
HeckeCentralMonomials                       central_monomials
HeckeCharValues                             char_values
HeckeClassPolynomials                       class_polynomials
HeckeReflectionRepresentation               reflrep
HermiteNormalFormIntegerMat                 hermite
HermiteNormalFormIntegerMatTransforms(m)    hermite_transforms(m)
HighestPowerFakeDegrees(W)                  charinfo(W)[:B]
HighestPowerFakeDegreeSymbol                degree_fegsymbol
HighestPowerGenericDegrees(W)               charinfo(W)[:A]
HighestPowerGenericDegreeSymbol             degree_gendeg_symbol
HyperplaneOrbits                            hyperplane_orbits
ICCTable                                    ICCTable
Incidence                                   incidence
IndependentLines(M)                         echelon(M)[2]
IndependentRoots                            independent_roots
InducedLinearForm                           induced_linear_form
InductionTable                              InductionTable
Inherit                                     look at merge for hashes
IntermediateGroup                           intermediate_group
Intersection                                intersect
IntFFE                                      Int
IntListToString                             joindigits
InvariantForm                               invariant_form
Invariants                                  invariants
Inversions                                  inversions
IsAbelian                                   isabelian
IsCyclic                                    iscyclic
IsCycPol(p)                                 p isa CycPol
IsFamily(f)                                 f isa Family
IsFFE(x)                                    x isa FFE
IsIsolated                                  is_isolated
IsJoinLattice                               is_join_lattice
IsLeftDescending(W,w,i)                     isleftdescent(W,w,i)
IsMeetLattice                               is_meet_lattice
IsomorphismType                             isomorphism_type
IsParabolic                                 is_parabolic
IsSubset(a,b)                               issubset(b,a)
IsUnipotentElement(x)                       x isa UnipotentElement
JInductionTable                             JInductionTable
jInductionTable                             jInductionTable
Join                                        join
KazhdanLusztigPolynomial                    KLPol
KroneckerProduct                            kron
LargestMovedPoint                           largest_moved_point
last                                        ans
LcmPartitions                               lcm_partitions
LeadingCoefficient(p)                       p[end]
LeftCell                                    LeftCell
LeftCells                                   LeftCells
LeftDescentSet(W,w)                         leftdescents(W,w)
LeftDivisorsSimple                          left_divisors
LeftGcd                                     leftgcd
LeftLcm                                     leftlcm
LinearExtension                             linear_extension
List(ConjugacyClasses(G),Representative)    classreps(G)
ListBlist(a,b)                              a[b]
ListPerm(p)                                 vec(p)
LogFFE                                      log
LongestCoxeterElement(W)                    longest(W)
LongestCoxeterWord(W)                       word(W,longest(W))
LowestPowerFakeDegrees(W)                   charinfo(W)[:b]
LowestPowerFakeDegreeSymbol                 valuation_fegsymbol
LowestPowerGenericDegrees(W)                charinfo(W)[:a]
LowestPowerGenericDegreeSymbol              valuation_gendeg_symbol
Lusztigaw                                   Lusztigaw
LusztigAw                                   LusztigAw
LusztigInduction                            LusztigInduce
LusztigInductionTable                       LusztigInductionTable
LusztigRestriction                          LusztigRestrict
M.ToOrdinary(i)                             B(M,i)
MappingPermListList                         mappingPerm
MatStab                                     stab_onmats
MatXPerm(W,p)                               reflrep(W,p)
MatYPerm                                    matY
Mod1                                        modZ
MovedPoints                                 support
Mvp("x")                                    Mvp(:x)
NrArrangements                              narrangements
NrCombinations                              ncombinations
NrConjugacyClasses                          nconjugacy_classes
NrDrinfeldDouble                            ndrinfeld_double
NrPartitions                                npartitions
NrPartitionsSet                             npartitions
NrPartitionTuples                           npartition_tuples
NrRestrictedPartitions                      npartitions
NullMat(m[,n])                              zeros(Int,m,m) resp. zeros(Int,m,n)
NullspaceIntMat                             lnullspaceInt
OnFamily(f,p::Int)                          galois(f,p)
OnFamily(f,p::Perm)                         f^p
OnMatrices(m,p)                             ^(m,p;dims=(1,2))
OnPolynomials(m,p)                          p^m
OnSets(s,g)                                 unique!(sort(s.^g))
OnTuples(l,p)                               l.^p
OrderedPartitions                           compositions
OrderFFE                                    order
OrderMod(n,m)                               order(Mod{m}(n))
ParabolicRepresentatives                    parabolic_reps
PartBeta                                    partβ
Partition                                   partition
Partitions                                  partitions
PartitionsSet                               partitions
PartitionTuples                             partition_tuples
PermCosetsSubgroup(H,W)                     D=vcat(reduced(H,W)...);map(s->Perm(reduced.(Ref(H),D.*s),D),gens(W))
PermList(v)                                 Perm(v)
PermListList(l1,l2)                         Perm(l1,l2)
PermMatMat                                  Perm_onmats
PermMatX                                    PermX
PermutationMat(p,dim)                       Matrix(p,dim)
Permuted(v,p)                               v^p
PermutedByCols(m,p)                         ^(m,p;dims=2)
Poset                                       Poset
Position(l,x)                               findfirst(==(x),l)
PositionCartesian                           cart2lin
PositionCartesian(a,b)                      LinearIndices(reverse(Tuple(a)))[CartesianIndices(Tuple(b))]
PositionClass                               position_class
PositionDet                                 charinfo(W)[:positionDet]
PositionId                                  charinfo(W)[:positionId]
PositionProperty(l,f)                       findfirst(f,l)
PositionRegularClass                        position_regular_class
Positions(l,x)                              findall(==(x),l)
PositionsProperty(l,f)                      findall(f,l)
PowerRoot(x,y)                              (Root1(;r=x)^y).r
Presentation                                Presentation
PrintDiagram(W)                             Diagram(W)
ProportionalityCoefficient(v,w)             ratio(v,w)
QuasiIsolatedRepresentatives                QuasiIsolatedRepresentatives
QuoInt                                      div
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
Reflections                                 reflections
ReflectionSubgroup                          reflection_subgroup
ReflectionType                              refltype
RegularEigenvalues                          regular_eigenvalues
RelativeDegrees                             relative_degrees
RelativeGroup                               relative_group
Replace                                     replace
Representations                             representations
RepresentativeConjugation(b,b'[,F][,type])  conjugating_elt(b,b'[,F],ss=type)
RepresentativeDiagonalConjugation           diagconj_elt
RepresentativeOperation                     transporting_elt
RepresentativeRowColPermutation             Perm_rowcol
Restricted                                  restricted
RestrictedPartitions                        partitions
RestrictedPerm(p,d)                         restricted(p,d)
Reversed                                    reverse
RightDescentSet(W,w)                        rightdescents(W,w)
RightGcd                                    rightgcd
RightLcm                                    rightlcm
RootDatum                                   rootdatum
RootsCartan(m)                              roots(m)
Rotation(v,i)                               circshift(v,-i)
Rotations(v)                                circshift.(Ref(v),length(v):-1:1)
ScalarProduct                               scalarproduct
ScalMvp                                     scalar
SchurElements                               schur_elements
SchurFunctor                                schur_functor
SemisimpleCentralizerRepresentatives        SScentralizer_reps
SemisimpleElement                           SS
SemisimpleRank                              semisimplerank
SemisimpleSubgroup                          torsion_subgroup
ShiftBeta                                   shiftβ
ShrinkGarsideGeneratingSet                  shrink
SignedMatStab                               sstab_onmats
SignedPerm                                  SPerm
SignedPermListList                          SPerm
SignedPermMatMat                            SPerm_onmats
Size(W)                                     length(W)
SmallestMovedPoint                          smallest_moved_point
SmithNormalFormIntegerMat                   smith
SmithNormalFormIntegerMatTransforms(m)      smith_transforms(m)
SolutionIntMat                              solutionmatInt
SolutionMat                                 solutionmat
SpecialPieces                               special_pieces
Spets                                       spets
SplitLevis                                  split_levis
StandardParabolic                           standard_parabolic
StandardParabolicClass                      standard_parabolic_class
StructureRationalPointsConnectedCentre      StructureRationalPointsConnectedCentre
SubSpets                                    subspets
SubTorus                                    SubTorus
SumIntersectionMat(m,n)                     (sum_rowspace(m,n),intersect_rowspace(m,n))
Symbols                                     HasType.BDSymbols
SymbolsDefect(e,r,def,ct)                   symbols(e,r,ct,def)
SymmetricDifference                         symdiff
SymmetricPower                              symmetric_power
Tableaux                                    tableaux
Torus                                       torus
TorusOrder                                  torus_order
TransitiveClosure                           transitive_closure
Transporter                                 transporter
TransposedMat                               transpose or permutedims
Transversals                                related to transversals
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
UnorderedTuples                             multisets
Valuation(p)                                valuation(p)
Value(p,x)                                  p(x)
W.matgens                                   reflrep(W)
W.matgens[i]                                reflrep(W,i)
W.N                                         nref(W)
W.orbitRepresentative                       simple_reps(W)
W.orbitRepresentativeElements[i]            simple_conjugating(W,i)
W.rootInclusion                             inclusion(W)
W.rootLengths                               rootlengths(W)
W.rootRestriction                           restriction(W)
W.roots                                     W.rootdec
W.simpleCoroots                             simplecoroots(W)
W.simpleRoots                               simpleroots(W)
WeightInfo                                  weightinfo
WGraph                                      Wgraph
WGraphToRepresentation                      WGraphToRepresentation
```
