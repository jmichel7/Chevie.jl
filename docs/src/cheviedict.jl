ChevieDict::Dict{String,String}=Dict(
"AbelianGenerators"=>"abelian_gens",
"AbelianInvariants"=>"abelian_invariants",
"Add"=>"push!",
"Affine"=>"affine",
#AffineRootAction
"AlgebraicCentre"=>"algebraic_center",
"AlmostCharacter"=>"almost_character or almostchar",
"Append"=>"append!",
"ApplyFunc(f,l)"=>"f(l...)",
"Arrangements"=>"arrangements",
"AsReflection"=>"reflection",
"AsFraction"=>"fraction",
"AsRootOfUnity"=>"Root1",
"AssociatedPartition"=>"conjugate_partition",
"AsymptoticAlgebra"=>"AsymptoticAlgebra",
"AsWord"=>"word",
"BadPrimes"=>"badprimes",
"BaseIntMat"=>"baseInt",
"Basis"=>"basis",
"BetaSet"=>"βset",
"BigCellDecomposition"=>"bigcell_decomposition",
"Binomial"=>"binomial",
"BipartiteDecomposition"=>"bipartite_decomposition",
"BlocksMat"=>"blocks",
"Braid"=>"BraidMonoid",
"BraidMonoid"=>"BraidMonoid",
"BraidRelations"=>"braid_relations",
"BrieskornNormalForm"=>"Brieskorn_normal_form",
"Bruhat"=>"bruhatless",
"BruhatPoset"=>"bruhatPoset",
"BruhatSmaller"=>"bruhatless",
"Catalan"=>"catalan",
"CartanMat(\"A\",5)"=>"cartan(:A,5)",
"CartanMatrix"=>"cartan",
"Cartesian"=>"cartesian",
"CartesianAt"=>"lin2cart",
"CartanMatFromCoxeterMatrix"=>"cartan",
"CentralIdempotents"=>"centralidempotents",
"Centralizer"=>"centralizer",
"CentralizerGenerators"=>"centralizer_gens",
"CharFFE"=>"char",
#CharName
"CharNames"=>"charnames",
"CharParams(W)"=>"charinfo(W).charparams",
"CharRepresentationWords"=>"traces_words_mats",
"CharTable"=>"CharTable",
"CheckHeckeDefiningRelations"=>"isrepresentation",
"ChevieClassInfo"=>"classinfo",
"ChevieCharInfo"=>"charinfo",
"ClassName"=>"see ClassNames",
"ClassTypes"=>"ClassTypes",
"Coefficient(p,i)"=>"p[i]",
"Collected"=>"tally",
"CollectBy(l,f)"=>"collectby(f,l)",
"Comm"=>"comm or commutator",
"Combinations"=>"combinations",
"ComplementIntMat"=>"complementInt",
"ComplexConjugate"=>"conj",
"ComplexReflectionGroup"=>"complex_reflection_group or crg",
"Compositions"=>"compositions",
"Concatenation(s::Vector)"=>"vcat(s...)",
"ConcatenationString(s...)"=>"prod([s...])",
"ConjugacyClasses"=>"conjugacy_classes",
"ConjugacySet(b[,F][,type])"=>"conjcat(b[,F],ss=type).obj",
"ConjugatePartition"=>"conjugate_partition",
"CoxeterCoset"=>"spets",
"CoxeterSubCoset"=>"subspets",
"CoxeterElements(W[,l])"=>"elements(W[,l])",
"CoxeterGroup(\"A\",5)"=>"coxeter_group(:A,5) or coxgroup",
"CoxeterGroupByCoxeterMatrix(C)"=>"coxeter_group(cartan(C)) or coxgroup",
"CoxeterGroupByCartanMatrix(C)"=>"coxeter_group(C) or coxgroup",
"CoxeterGroupHyperoctaedralGroup(n)"=>"coxeter_hyperoctaedral_group(n) or coxhyp",
"CoxeterGroupSymmetricGroup(n)"=>"coxeter_symmetric_group(n) or coxsym",
"CoxeterLength(W,w)"=>"length(W,w)",
"CoxeterMatrix"=>"coxmat or coxeter_matrix",
"CoxeterMatrixFromCartanMat"=>"coxmat or coxeter_matrix",
"CoxeterNumber"=>"coxnum or coxeter_number",
"CoxeterWord(W,w)"=>"word(W,w)",
"CoxeterWords(W[,l])"=>"word.(Ref(W),elements(W[,l]))",
"CuspidalPairs"=>"cuspidal_data",
"CuspidalUnipotentCharacters(W[,d])"=>"cuspidal(UnipotentCharacters(W)[,d])",
"CyclotomicPolynomial(R,i)"=>"cyclotomic_polynomial(i)",
"Cycle"=>"orbit",
"Cycles"=>"orbits",
"CyclotomicModP(c,p)"=>"FFE{p}(c)",
"CycPol"=>"CycPol",
"CycPolFakeDegreeSymbol"=>"fegsymbol",
"CycPolGenericDegreeSymbol"=>"gendeg_symbol",
"CycPolUnipotentDegrees(W)"=>"CycPoldegrees(UnipotentCharacters(W))",
"DecomposedMat"=>"diagblocks",
"DefectSymbol"=>"defectsymbol",
"Degree(p)"=>"degree(p)",
"DegreeFFE"=>"degree",
"M::GarsideMonoid.delta"=>"M.δ",
"DeligneLusztigCharacter"=>"deligne_lusztig_character or dlchar",
"DeligneLusztigLefschetz"=>"deligne_lusztig_leftschetz or dlleftschetz",
"DescribeInvolution"=>"describe_involution",
"DetPerm(W)"=>"perm(detPerm(W))",
"DiaconisGraham"=>"diaconis_graham",
"DiagonalMat"=>"Diagonal or cat",
"DiagonalOfMat"=>"diag",
#DifferenceMultiSet
"Digits"=>"digits",
"Dimension"=>"dim or dimension",
#Discriminant
"DistinguishedParabolicSubgroups"=>"distinguished_parabolics",
"Dominates"=>"dominates",
"DrinfeldDouble"=>"drinfeld_double",
"Drop"=>"deleteat!",
#Dual
"DualBraid"=>"DualBraidMonoid",
"DualBraidMonoid"=>"DualBraidMonoid",
"EigenspaceProjector"=>"eigenspace_projector",
"EigenvaluesMat"=>"eigmat",
"Elements"=>"elements",
"ElementWithInversions(W,l)"=>"with_inversions(W,l)",
"EltBraid"=>"image",
"EltWord(W,w)"=>"W(w...)",
"ER"=>"root",
"ExteriorPower"=>"exterior_power",
"FactorizedSchurElement"=>"FactorizedSchurElement",
"FactorizedSchurElements"=>"FactorizedSchurElements",
"FakeDegree"=>"fakedegree",
"FakeDegrees"=>"fakedegrees",
"FamiliesClassical"=>"FamiliesClassical",
"Family"=>"Family",
"FamilyImprimitive"=>"family_imprimitive",
"Filtered(l,f)"=>"filter(f,l)",
"FiniteCoxeterTypeFromCartanMat(m)"=>"type_cartan(m)",
"FirstLeftDescending(W,w)"=>"firstleftdescent(W,w)",
"ForAll(l::list,f::function)"=>"all(f,l)",
"ForAny(l::list,f::function)"=>"any(f,l)",
"ForEachCoxeterWord(W,f)"=>"for w in W f(word(W,w)) end",
"ForEachElement(W,f)"=>"for w in W f(w) end",
"FormatTable"=>"showtable",
"Frobenius"=>"Frobenius",
"FullSymbol"=>"fullsymbol",
"FundamentalGroup"=>"fundamental_group",
"FusionAlgebra"=>"Zbasedring",
"FusionConjugacyClasses"=>"fusion_conjugacy_classes",
"GaloisCyc"=>"galois",
"GarsideAlpha"=>"α",
"GarsideWords"=>"elements",
"GcdPartitions"=>"gcd_partitions",
"GcdRepresentation(x,y)"=>"gcdx(x,y)[2:3]",
"W.generators"=>"gens(W) or generators(W)",
"Length(W.generators)"=>"ngens(W) or number_of_generators(W)",
"List(l::list,f::function)"=>"map(f,l)",
"GenericOrder"=>"generic_order",
"GenericSign"=>"generic_sign",
"GetRoot"=>"root",
#GoodCoxeterWord
"GraphAutomorphisms"=>"graph_automorphisms",
"Hasse"=>"hasse",
"Hecke"=>"hecke",
"HeckeCharValues"=>"char_values",
#HeckeCharValuesGood
"HeckeCentralMonomials"=>"central_monomials",
"HeckeClassPolynomials"=>"class_polynomials",
"HeckeReflectionRepresentation"=>"reflection_representation or reflrep",
#HeckeSubAlgebra
"HermiteNormalFormIntegerMat"=>"hermite",
"HermiteNormalFormIntegerMatTransforms(m)"=>"hermite_transforms(m)",
"HighestPowerFakeDegrees(W)"=>"charinfo(W).B",
"HighestPowerFakeDegreeSymbol"=>"degree_fegsymbol",
"HighestPowerGenericDegrees(W)"=>"charinfo(W).A",
"HighestPowerGenericDegreeSymbol"=>"degree_gendeg_symbol",
"HighestShortRoot"=>"highest_short_root",
"KazhdanLusztigPolynomial"=>"KLPol",
"HyperplaneOrbits"=>"hyperplane_orbits",
"ICCTable"=>"ICCTable",
"Idempotents"=>"idempotents",
"IdentityMat(nn)"=>"Matrix(1I,n,n)",
"Incidence"=>"incidence",
"IndependentLines"=>"independent_rows",
"IndependentRoots"=>"independent_roots",
"InducedLinearForm"=>"induced_linear_form",
"InductionTable"=>"induction_table",
"Inherit"=>"look at merge for hashes",
"Intersection"=>"intersect",
"IntermediateGroup"=>"intermediate_group",
"IntFFE"=>"Int",
"IntListToString"=>"joindigits",
"InvariantForm"=>"invariant_form",
"Invariants"=>"invariants",
"Inversions"=>"inversions",
"IsAbelian"=>"isabelian",
"IsAssociative"=>"isassociative",
"IsCyclic"=>"iscyclic",
"IsCycPol(p)"=>"p isa CycPol",
"IsFamily(f)"=>"f isa Family",
"IsFFE(x)"=>"x isa FFE",
"IsIsolated"=>"isisolated",
"IsJoinLattice"=>"is_join_semilattice",
"IsMeetLattice"=>"is_meet_semilattice",
"IsLeftDescending(W,w,i)"=>"isleftdescent(W,w,i)",
"IsSubset(a,b)"=>"issubset(b,a)",
"IsParabolic"=>"isparabolic",
#IsNormalizing
#IsQuasiIsolated
"IsomorphismType"=>"isomorphism_type",
"IsUnipotentElement(x)"=>"x isa UnipotentElement",
"jInductionTable"=>"j_induction_table",
"JInductionTable"=>"J_induction_table",
"Join"=>"join",
"KroneckerProduct"=>"kron",
"LargestMovedPoint"=>"last_moved",
"last"=>"ans",
"LcmPartitions"=>"lcm_partitions",
"LeadingCoefficient(p)"=>"p[end]",
"LeftCell"=>"LeftCell",
"LeftCells"=>"left_cells",
"LeftDescentSet(W,w)"=>"leftdescents(W,w)",
"LeftDivisorsSimple"=>"left_divisors",
"LeftGcd"=>"leftgcd",
"LeftLcm"=>"leftlcm",
"M.LeftLcmSimples(x...)"=>"leftlcm(M,...)",
"LinearExtension"=>"linear_extension",
"List(l,f)"=>"map(f,l)",
"ListBlist(a,b)"=>"a[b]",
"ListPerm(p)"=>"perm(p)",
"List(ConjugacyClasses(G),Representative)"=>"classreps(G) or class_representatives(G)",
"LogFFE"=>"log",
"LongestCoxeterElement(W)"=>"longest(W)",
"LongestCoxeterWord(W)"=>"word(W,longest(W))",
"LoewyLength"=>"loewylength",
"LowestPowerFakeDegrees(W)"=>"charinfo(W).b",
"LowestPowerFakeDegreeSymbol"=>"valuation_fegsymbol",
"LowestPowerGenericDegrees(W)"=>"charinfo(W).a",
"LowestPowerGenericDegreeSymbol"=>"valuation_gendeg_symbol",
"Lusztigaw"=>"Lusztigaw",
"LusztigAw"=>"LusztigAw",
"LusztigInduction"=>"lusztig_induce",
"LusztigInductionTable"=>"lusztig_induction_table",
"LusztigRestriction"=>"lusztig_restrict",
"MappingPermListList"=>"mappingPerm",
"W.matgens"=>"reflection_representation(W) or reflrep",
"W.matgens[i]"=>"reflection_representation(W,i) or reflrep",
"MatStab"=>"stab_onmats",
"MatXPerm(W,p)"=>"reflection_representation(W,p) or reflrep",
"MatYPerm"=>"YMatrix",
"Maximum"=>"max or maximum",
"Minimum"=>"min or minimum",
"Mod1"=>"modZ",
"MovedPoints"=>"support",
"Mvp(\"x\")"=>"Mvp(:x)",
"W.N"=>"nref(W) or number_of_reflections(W)",
"W.Nhyp"=>"nhyp(W) or number_of_hyperplanes(W)",
"NrArrangements"=>"narrangements",
"NrCombinations"=>"ncombinations",
"NrConjugacyClasses"=>"nconjugacy_classes",
"NrDrinfeldDouble"=>"ndrinfeld_double",
"NrPartitions"=>"npartitions",
"NrPartitionsSet"=>"npartitions",
"NrPartitionTuples"=>"npartition_tuples",
"NrRestrictedPartitions"=>"npartitions",
"NullMat(m[,n])"=>"zeros(Int,m,m) resp. zeros(Int,m,n)",
"NullspaceIntMat"=>"lnullspaceInt",
"W.Nhyp"=>"number_of_hyperplanes(W) or nyp(W)",
"Number(l,f)"=>"count(f,l)",
"OnFamily(f,p::Perm)"=>"f^p",
"OnFamily(f,p::Int)"=>"galois(f,p)",
"OnMatrices"=>"onmats",
"OnSets"=>"onsets",
"OnTuples"=>"ontuples",
"OnPolynomials(m,p)"=>"p^m",
"W.orbitRepresentative"=>"simple_reps(W)",
"W.orbitRepresentativeElements[i]"=>"simple_conjugating(W,i)",
"OrderedPartitions"=>"compositions",
"OrderFFE"=>"order",
"OrderMod(n,m)"=>"order(Mod{m}(n))",
"W.OrdersGeneratingReflections"=>"ordergens(W) or orders_of_generators(W)",
"ParabolicClosure"=>"parabolic_closure",
"ParabolicRepresentatives"=>"parabolic_reps",
#ParabolicSubgroups
"PartBeta"=>"partβ",
"Partition"=>"partition",
"Partitions"=>"partitions",
"PartitionsSet"=>"partitions",
"PartitionTuples"=>"partition_tuples",
"PartitionTupleToString"=>"string_partition_tuple",
"PermCosetsSubgroup(H,W)"=>"D=vcat(reduced(H,W)...);map(s->Perm(reduced.(Ref(H),D.*s),D),gens(W))",
"PermListList(l1,l2)"=>"Perm(l1,l2)",
"PermList(v)"=>"Perm(v)",
"PermMatMat(m,n)"=>"Perm(m,n;dims=(1,2))",
"PermMatX"=>"PermX",
"PermMatY"=>"PermY",
"PermutationMat(p,dim)"=>"Matrix(p,dim)",
"PermutationOnClasses"=>"on_classes",
"PermutationOnCharacters"=>"on_chars",
"PermutationOnUnipotents"=>"on_unipotents",
"Permuted(v,p)"=>"invpermute(v,p)",
"PermutedByCols(m,p)"=>"invpermute(m,p;dims=2)",
#PoincarePolynomial
"Poset"=>"Poset",
"Position(l,x)"=>"findfirst(==(x),l)",
"PositionClass"=>"position_class",
"PositionCartesian(a,b)"=>"cart2lin(a,b)=LinearIndices(reverse(Tuple(a)))[reverse(b)...]",
"PositionDet"=>"charinfo(W).positionDet",
"PositionId"=>"charinfo(W).positionId",
"PositionRegularClass"=>"position_regular_class",
"Positions(l,x)"=>"findall(==(x),l)",
"PositionProperty(l,f)"=>"findfirst(f,l)",
"PositionsProperty(l,f)"=>"findall(f,l)",
"PowerRoot(x,y)"=>"(Root1(;r=x)^y).r",
"Presentation"=>"Presentation",
"PrintDiagram(W)"=>"diagram(W)",
"Product"=>"prod",
"ProportionalityCoefficient(v,w)"=>"ratio(v,w)",
"QuasiIsolatedRepresentatives"=>"quasi_isolated_reps",
"QuoInt"=>"div",
"Radical"=>"radical",
"RadicalPower"=>"radicalpower",
"Rank"=>"rank",
"RankSymbol"=>"ranksymbol",
"RecFields"=>"propertynames",
"ReducedCoxeterWord(W,w)"=>"word(W,W(w...))",
"ReducedExpressions(W,w)"=>"words(W,w)",
"ReducedInRightCoset(W,w)"=>"reduced(W,w)",
"ReducedRightCosetRepresentatives(W,H)"=>"reduced(H,W)",
"Reflection"=>"refls(W,i) or reflectionMatrix(root,coroot)",
"ReflectionCharacter"=>"reflection_character or reflchar",
"ReflectionCharValue"=>"tr(reflrep(W,w))",
#ReflectionCoset
"ReflectionDegrees(W)"=>"degrees(W)",
"ReflectionCoDegrees(W)"=>"codegrees(W)",
"ReflectionEigenvalues"=>"refleigen",
"ReflectionGroup"=>"reflection_group",
"ReflectionLength(W,w)"=>"reflength(W,w)",
#ReflectionWord
"ReflectionName(W)"=>"repr(W;context=:limit=>true)",
"Reflections"=>"Perm.(reflections(W)[1:nhyp(W)])",
#ReflectionSubCoset
"ReflectionSubgroup"=>"reflection_subgroup",
"ReflectionType"=>"refltype",
"RegularEigenvalues"=>"regular_eigenvalues",
"RelativeDegrees"=>"relative_degrees",
"RelativeGroup"=>"relative_group",
"Replace"=>"replace",
"Representations"=>"representations",
"RepresentativeConjugation(b,b'[,F][,type])"=>"conjugating_elt(b,b'[,F],ss=type)",
"RepresentativeDiagonalConjugation"=>"diagconj_elt",
"RepresentativeOperation"=>"transporting_elt or transporting_element",
"RepresentativeRowColPermutation"=>"Perm_rowcol",
"Restricted"=>"restricted",
"RestrictedPartitions"=>"partitions",
"RestrictedPerm(p,d)"=>"restricted(p,d)",
"Reversed"=>"reverse",
"ReversedWord"=>"reverse",
"RightDescentSet(W,w)"=>"rightdescents(W,w)",
"RightGcd"=>"rightgcd",
"RightLcm"=>"rightlcm",
"M.RightLcmSimples(x...)"=>"rightlcm(M,...)",
"RootDatum"=>"rootdatum",
"W.rootLengths"=>"rootlengths(W)",
"W.roots"=>"W.rootdec",
"RootsCartan(m)"=>"roots(m)",
"W.rootInclusion"=>"inclusion(W)",
"W.rootRestriction"=>"restriction(W)",
"Rotation(v,i)"=>"circshift(v,-i)",
"Rotations(v)"=>"circshift.(Ref(v),length(v):-1:1)",
"ScalarProduct"=>"Chars.scalarproduct",
"ScalMvp"=>"scalar",
#SchurElement
"SchurElements"=>"schur_elements",
"SchurFunctor"=>"schur_functor",
"SemisimpleCentralizerRepresentatives"=>"semisimple_centralizer_representatives or sscentralizer_reps",
"SemisimpleElement"=>"ss",
"SemisimpleRank"=>"semisimplerank",
"SemisimpleSubgroup"=>"torsion_subgroup",
"ShiftBeta"=>"shiftβ",
"ShrinkGarsideGeneratingSet"=>"shrink",
"SignedMatStab"=>"sstab_onmats",
"SignedPerm"=>"SPerm",
"SignedPermListList"=>"SPerm",
"SignedPermMatMat(M,N)"=>"SPerm(M,N;dims=(1,2))",
"W.simpleCoroots"=>"simplecoroots(W)",
"W.simpleRoots"=>"simpleroots(W)",
"Size(W)"=>"length(W)",
"SmallestMovedPoint"=>"first_moved",
"SmithNormalFormIntegerMat"=>"smith",
"SmithNormalFormIntegerMatTransforms(m)"=>"smith_transforms(m)",
"SolutionMat"=>"solutionmat",
"SolutionIntMat"=>"solutionmatInt",
"Sort"=>"sort!",
"SortBy(l,f)"=>"sort!(l,by=f)",
"SortingPerm(a)"=>"inv(sortPerm(a))",
"SortParallel(a,b)"=>"b=b[sortperm(a)];sort!(a)",
"SpecialPieces"=>"special_pieces",
"Spets"=>"spets",
"Split"=>"split",
"SplitLevis"=>"split_levis",
"Sprint"=>"string",
"Stabilizer"=>"stabilizer",
"StandardParabolic"=>"standard_parabolic",
"StandardParabolicClass"=>"standard_parabolic_class",
"StructureRationalPointsConnectedCentre"=>"structure_rational_points_connected_centre",
"SubSpets"=>"subspets",
"SubTorus"=>"SubTorus",
"Sum"=>"sum",
"SumIntersectionMat(m,n)"=>"(rowspace(vcat(m,n)),intersect_rowspace(m,n))",
"Symbols"=>"BDSymbols",
"SymbolsDefect(e,r,def,ct)"=>"symbols(e,r,ct,def)",
"SymmetricDifference"=>"symdiff",
"SymmetricPower"=>"symmetric_power",
"Tableaux"=>"tableaux",
"Tensored(c,d)"=>"vec([i.*j for i in c,j ind])",
"M.ToOrdinary(i)"=>"B(M,i)",
"Torus"=>"torus",
"TorusOrder"=>"torus_order",
"TraceMat"=>"tr",
"TransitiveClosure"=>"transitive_closure",
"Transporter"=>"transporter",
"TransposedMat"=>"transpose or permutedims",
"Transversals"=>"related to transversal and orbits",
"TriangulizeMat"=>"echelon!",
"Twistings"=>"twistings",
"TwoTree(m)"=>"twotree(m)",
"UnipotentAbelianPart"=>"abelianpart",
"UnipotentCharacter"=>"unipotent_character or unichar",
"UnipotentCharacters"=>"UnipotentCharacters",
"UnipotentClasses"=>"UnipotentClasses",
"UnipotentDecompose"=>"decompose",
"UnipotentDegrees(W,q)"=>"degrees(UnipotentCharacters(W),q)",
"UnipotentGroup"=>"UnipotentGroup",
"UnorderedTuples"=>"multisets",
"Valuation(p)"=>"valuation(p)",
"Value(p,x)"=>"p(x)",
"ValuePol"=>"evalpoly",
"WeightInfo"=>"weightinfo",
"WGraph"=>"Wgraph",
"WGraphToRepresentation"=>"WGraphToRepresentation",
"Zip(a1,..,an,f)"=>"map(f,a1,...,an)"
)

function gap(s)
  if !(s isa Regex) s=Regex(s,"i") end
  kk=filter(x->occursin(s,x),keys(ChevieDict))
  if isempty(kk) 
    println("no match")
    return
  end
  pad=maximum(length(k) for k in kk)+2
  print(join(sort([rpad(k,pad)*"=>  "*ChevieDict[k]*"\n" for k in kk]),""))
end

function fixdoc()
  pad=maximum(length(k) for k in keys(ChevieDict))+2
  u=[rpad(k,pad)*v*"\n" for (k,v) in ChevieDict]
  sort!(u,by=lowercase)
  open("dict.md","w")do f
    write(f,"""
# Dictionary from GAP3/Chevie
The dictionary from GAP3/Chevie is as follows:
```
""")
    write(f,join(u))
    write(f,"```")
  end
  length(u)
end
