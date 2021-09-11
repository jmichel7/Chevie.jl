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
Group
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
abelian_gens
```
# Permutation groups
```@docs
PermGroups
base
centralizers
transversals
on_classes
symmetric_group
stab_onmats
Perm_onmats
Perm_rowcolmat
Base.in(::Perm,::PermGroup)
```
# Cyclotomic numbers
```@docs
Cycs
conductor
coefficients(c::Cyc)
denominator(c::Cyc{Rational})
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
Pol
divrem
gcd(::Pol,::Pol)
gcdx(::Pol,::Pol)
cyclotomic_polynomial
```
# Multivariate (Puiseux) polynomials
```@docs
Mvps
Mvp
variables
Mvps.coefficients(::Mvp,::Symbol)
Mvps.coefficient
Pol(::Mvp)
Pol(::Mvp,::Symbol)
valuation
degree
Mvps.value
Mvps.degree(::Mvp)
Mvps.conj
factor(::Mvp)
Mvps.derivative
laurent_denominator
gcd(::Mvp,::Mvp)
lcm(::Mvp,::Mvp)
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
showtable
cut
prime_residues
phi
primitiveroot
```
# Combinatorics
```@docs
Combinat
arrangements
combinations
partitions
partition_tuples
restrictedpartitions
compositions
submultisets
partitions_set
lcm_partitions
gcd_partitions
conjugate_partition
dominates
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
hasse
incidence
transitive_closure
linear_extension
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
charpoly
comatrix
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
smith_normal_form
baseInt
complementInt
DiaconisGraham
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
badprimes
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
left_divisors
leftgcd
rightgcd
leftlcm
rightlcm
α(::Garside.LocallyGarsideElm)
α(::GarsideElm,::AbstractVector)
Brieskorn_normal_form
BraidMonoid
DualBraidMonoid
hurwitz
fraction
word(::Garside.GarsideMonoid,w)
word(::Garside.GarsideElm)
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
intermediate_group
QuasiIsolatedRepresentatives(::FiniteCoxeterGroup)
is_isolated
torsion_subgroup
algebraic_centre
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
cuspidal_pairs
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
```
# Symbols
```@docs
Symbols
shiftβ
βset
partβ
ranksymbol
defectsymbol
symbol_partition_tuple
fegsymbol
gendeg_symbol
degree_fegsymbol
degree_gendeg_symbol
valuation_fegsymbol
valuation_gendeg_symbol
tableaux
symbols
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
AbelianGenerators                           abelian_gens
AlgebraicCentre                             algebraic_centre
AlmostCharacter                             AlmostChar
Arrangements                                arrangements
AsFraction                                  fraction
AsReflection                                reflection
AsRootOfUnity                               Root1
AsWord                                      word
AssociatedPartition                         conjugate_partition
AsymptoticAlgebra                           AsymptoticAlgebra
BadPrimes                                   badprimes
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
DetPerm(W)                                  vec(detPerm(W))
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
HighestPowerFakeDegreeSymbol                degree_fegsymbol
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
Inherit                                     look at merge for hashes
IntFFE                                      Int
IntListToString                             joindigits
IntermediateGroup                           intermediate_group
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
LeftLcm                                     leftlcm
LinearExtension                             linear_extension
List(ConjugacyClasses(G),Representative)    classreps(G)
ListBlist(a,b)                              a[b]
ListPerm(p)                                 vec(p)
LogFFE                                      log
LongestCoxeterElement(W)                    longest(W)
LongestCoxeterWord(W)                       word(W,longest(W))
LowestPowerFakeDegreeSymbol                 valuation_fegsymbol
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
NrCombinations                              ncombinations
NrConjugacyClasses                          nconjugacy_classes
NrDrinfeldDouble                            ndrinfeld_double
NrPartitionTuples                           npartition_tuples
NrPartitions                                npartitions
NrPartitionsSet                             npartitions_set
OnFamily(f,p::Int)                          galois(f,p)
OnFamily(f,p::Perm)                         f^p
OnMatrices(m,p)                             ^(m,p;dims=(1,2))
OnPolynomials(m,p)                          p^m
OnSets(s,g)                                 unique!(sort(s.^g))
OnTuples(l,p)                               l.^p
OrderFFE                                    order
OrderMod(n,m)                               order(Mod{m}(n))
OrderedPartitions                           compositions
ParabolicRepresentatives                    parabolic_reps
PartBeta                                    partβ
Partition                                   partition
PartitionTuples                             partition_tuples
Partitions                                  partitions
PartitionsSet                               partitions_set
PermCosetsSubgroup(H,W)                     D=reduced(H,W);map(s->Perm(reduced.(Ref(H),D.*s),D),gens(W))
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
PositionDet                                 charinfo(W)[:positionDet]
PositionId                                  charinfo(W)[:positionId]
PositionRegularClass                        position_regular_class
Presentation                                Presentation
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
Replace                                     replace
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
Rotations(v)                                circshift.(Ref(v),length(v):-1:1)
ScalMvp                                     scal
ScalarProduct                               scalarproduct
SchurElements                               schur_elements
SchurFunctor                                schur_functor
SemisimpleCentralizerRepresentatives        SScentralizer_reps
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
SpecialPieces                               special_pieces
Spets                                       spets
SplitLevis                                  split_levis
StandardParabolic                           standard_parabolic
StandardParabolicClass                      standard_parabolic_class
StructureRationalPointsConnectedCentre      StructureRationalPointsConnectedCentre
SubSpets                                    subspets
SubTorus                                    SubTorus
Symbols                                     HasType.BDSymbols
SymbolsDefect(e,r,def,ct)                   symbols(e,r,ct,def)
SymmetricDifference                         symdiff
SymmetricPower                              symmetric_power
Tableaux                                    tableaux
Torus                                       torus
TorusOrder                                  torus_order
TransitiveClosure                           transitive_closure
Transporter                                 transporter
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
UnorderedTuples                             submultisets
Valuation(p)                                valuation(p)
Value(p,x)                                  p(x)
W.N                                         nref(W)
W.matgens[i]                                reflrep(W,i)
W.orbitRepresentative                       simple_reps(W)
W.orbitRepresentativeElements[i]            simple_conjugating(W,i)
W.rootInclusion                             inclusion(W)
W.rootLengths                               rootlengths(W)
W.rootRestriction                           restriction(W)
WGraph                                      Wgraph
WGraphToRepresentation                      WGraphToRepresentation
WeightInfo                                  weightinfo
jInductionTable                             jInductionTable
last                                        ans
```
