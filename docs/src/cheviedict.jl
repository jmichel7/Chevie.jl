const ChevieDict=Dict(
"AbelianGenerators"=>"abelian_gens",
"AbelianInvariants"=>"abelian_invariants",
"Affine"=>"affine",
#AffineRootAction
"AlgebraicCentre"=>"algebraic_center",
"AlmostCharacter"=>"AlmostChar",
"Arrangements"=>"arrangements",
"AsReflection"=>"reflection",
"AsFraction"=>"fraction",
"AsRootOfUnity"=>"Root1",
"AssociatedPartition"=>"conjugate_partition",
"AsymptoticAlgebra"=>"AsymptoticAlgebra",
"AsWord"=>"word",
"BadPrimes"=>"badprimes",
"BaseIntMat"=>"baseInt",
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
"BruhatPoset"=>"Poset",
"BruhatSmaller"=>"bruhatless",
"Catalan"=>"catalan",
"CartanMat(\"A\",5)"=>"cartan(:A,5)",
"Cartesian"=>"cartesian",
"CartesianAt"=>"lin2cart",
"CartanMatFromCoxeterMatrix"=>"cartan",
"Centralizer"=>"centralizer",
"CentralizerGenerators"=>"centralizer_gens",
"CharFFE(x)"=>"field(x).p",
#CharName
"CharNames"=>"charnames",
"CharParams(W)"=>"charinfo(W)[:charparams]",
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
"Comm"=>"comm",
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
"CoxeterGroupHyperoctaedralGroup(n)"=>"CoxHyperoctaedral(n)",
"CoxeterGroupSymmetricGroup(n)"=>"CoxSym(n)",
"CoxeterLength(W,w)"=>"length(W,w)",
"CoxeterMatrix"=>"coxmat",
"CoxeterMatrixFromCartanMat"=>"coxmat",
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
"CycPolUnipotentDegrees"=>"CycPolUnipotentDegrees",
"DecomposedMat"=>"diagblocks",
"DefectSymbol"=>"defectsymbol",
"Degree(p)"=>"degree(p)",
"DegreeFFE(x)"=>"field(x).n",
"DeligneLusztigCharacter"=>"DLChar",
"DeligneLusztigLefschetz"=>"DLLeftschetz",
"DescribeInvolution"=>"describe_involution",
"DetPerm(W)"=>"vec(detPerm(W))",
"DiaconisGraham"=>"diaconis_graham",
"DiagonalMat"=>"Diagonal or cat",
"DiagonalOfMat"=>"diag",
#DifferenceMultiSet
"Digits"=>"digits",
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
"FiniteCoxeterTypeFromCartanMat(m)"=>"type_cartan(m)",
"FirstLeftDescending(W,w)"=>"firstleftdescent(W,w)",
"ForEachCoxeterWord(W,f)"=>"for w in W f(word(W,w)) end",
"ForEachElement(W,f)"=>"for w in W f(w) end",
"FormatTable"=>"showtable",
"Frobenius"=>"Frobenius",
"FullSymbol"=>"fullsymbol",
"FundamentalGroup"=>"fundamental_group",
"FusionAlgebra"=>"fusion_algebra",
"FusionConjugacyClasses"=>"fusion_conjugacy_classes",
"GaloisCyc"=>"galois",
"GarsideAlpha"=>"α",
"GarsideWords"=>"elements",
"GcdPartitions"=>"gcd_partitions",
"GcdRepresentation(x,y)"=>"gcdx(x,y)[2:3]",
"W.generators"=>"gens(W)",
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
"HighestPowerFakeDegrees(W)"=>"charinfo(W)[:B]",
"HighestPowerFakeDegreeSymbol"=>"degree_fegsymbol",
"HighestPowerGenericDegrees(W)"=>"charinfo(W)[:A]",
"HighestPowerGenericDegreeSymbol"=>"degree_gendeg_symbol",
"HighestShortRoot"=>"highest_short_root",
"KazhdanLusztigPolynomial"=>"KLPol",
"HyperplaneOrbits"=>"hyperplane_orbits",
"ICCTable"=>"ICCTable",
"Incidence"=>"incidence",
"IndependentLines(M)"=>"echelon(M)[2]",
"IndependentRoots"=>"independent_roots",
"InducedLinearForm"=>"induced_linear_form",
"InductionTable"=>"InductionTable",
"Inherit"=>"look at merge for hashes",
"Intersection"=>"intersect",
"IntermediateGroup"=>"intermediate_group",
"IntFFE"=>"Int",
"IntListToString"=>"joindigits",
"InvariantForm"=>"invariant_form",
"Invariants"=>"invariants",
"Inversions"=>"inversions",
"IsAbelian"=>"isabelian",
"IsCyclic"=>"iscyclic",
"IsCycPol(p)"=>"p isa CycPol",
"IsFamily(f)"=>"f isa Family",
"IsFFE(x)"=>"x isa FFE",
"IsIsolated"=>"isisolated",
"IsJoinLattice"=>"isjoin_lattice",
"IsMeetLattice"=>"ismeet_lattice",
"IsLeftDescending(W,w,i)"=>"isleftdescent(W,w,i)",
"IsSubset(a,b)"=>"issubset(b,a)",
"IsParabolic"=>"isparabolic",
#IsNormalizing
#IsQuasiIsolated
"IsomorphismType"=>"isomorphism_type",
"IsUnipotentElement(x)"=>"x isa UnipotentElement",
"jInductionTable"=>"jInductionTable",
"JInductionTable"=>"JInductionTable",
"Join"=>"join",
"KroneckerProduct"=>"kron",
"LcmPartitions"=>"lcm_partitions",
"LargestMovedPoint"=>"largest_moved_point",
"last"=>"ans",
"LeadingCoefficient(p)"=>"p[end]",
"LeftCell"=>"LeftCell",
"LeftCells"=>"left_cells",
"LeftDescentSet(W,w)"=>"leftdescents(W,w)",
"LeftDivisorsSimple"=>"left_divisors",
"LeftGcd"=>"leftgcd",
"LeftLcm"=>"leftlcm",
"M.LeftLcmSimples(x...)"=>"leftlcm(M,...)",
"LinearExtension"=>"linear_extension",
"ListBlist(a,b)"=>"a[b]",
"ListPerm(p)"=>"vec(p)",
"List(ConjugacyClasses(G),Representative)"=>"classreps(G)",
"LogFFE"=>"log",
"LongestCoxeterElement(W)"=>"longest(W)",
"LongestCoxeterWord(W)"=>"word(W,longest(W))",
"LowestPowerFakeDegrees(W)"=>"charinfo(W)[:b]",
"LowestPowerFakeDegreeSymbol"=>"valuation_fegsymbol",
"LowestPowerGenericDegrees(W)"=>"charinfo(W)[:a]",
"LowestPowerGenericDegreeSymbol"=>"valuation_gendeg_symbol",
"Lusztigaw"=>"Lusztigaw",
"LusztigAw"=>"LusztigAw",
"LusztigInduction"=>"LusztigInduce",
"LusztigInductionTable"=>"LusztigInductionTable",
"LusztigRestriction"=>"LusztigRestrict",
"MappingPermListList"=>"mappingPerm",
"W.matgens"=>"reflection_representation(W) or reflrep",
"W.matgens[i]"=>"reflection_representation(W,i) or reflrep",
"MatStab"=>"stab_onmats",
"MatXPerm(W,p)"=>"reflection_representation(W,p) or reflrep",
"MatYPerm"=>"matY",
"Mod1"=>"modZ",
"MovedPoints"=>"support",
"Mvp(\"x\")"=>"Mvp(:x)",
"W.N"=>"nref(W) or number_of_reflections(W)",
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
"OnFamily(f,p::Perm)"=>"f^p",
"OnFamily(f,p::Int)"=>"galois(f,p)",
"OnMatrices(m,p)"=>"^(m,p;dims=(1,2))",
"OnSets(s,g)"=>"unique!(sort(s.^g))",
"OnTuples(l,p)"=>"l.^p",
"OnPolynomials(m,p)"=>"p^m",
"W.orbitRepresentative"=>"simple_reps(W)",
"W.orbitRepresentativeElements[i]"=>"simple_conjugating(W,i)",
"OrderedPartitions"=>"compositions",
"OrderFFE"=>"order",
"OrderMod(n,m)"=>"order(Mod{m}(n))",
"ParabolicClosure"=>"parabolic_closure",
"ParabolicRepresentatives"=>"parabolic_reps",
#ParabolicSubgroups
"PartBeta"=>"partβ",
"Partition"=>"partition",
"Partitions"=>"partitions",
"PartitionsSet"=>"partitions",
"PartitionTuples"=>"partition_tuples",
#PartitionTupleToString
"PermCosetsSubgroup(H,W)"=>"D=vcat(reduced(H,W)...);map(s->Perm(reduced.(Ref(H),D.*s),D),gens(W))",
"PermListList(l1,l2)"=>"Perm(l1,l2)",
"PermList(v)"=>"Perm(v)",
"PermMatMat(m,n)"=>"Perm(m,n;dims=(1,2))",
"PermMatX"=>"PermX",
#PermMatY
"PermutationMat(p,dim)"=>"Matrix(p,dim)",
"PermutationOnClasses"=>"on_classes",
"PermutationOnCharacters"=>"on_chars",
"PermutationOnUnipotents"=>"on_unipotents",
"Permuted(v,p)"=>"permute(v,p)",
"PermutedByCols(m,p)"=>"permute(m,p;dims=2)",
#PoincarePolynomial
"Poset"=>"Poset",
"Position(l,x)"=>"findfirst(==(x),l)",
"PositionClass"=>"position_class",
"PositionCartesian(a,b)"=>"LinearIndices(reverse(Tuple(a)))[CartesianIndices(Tuple(b))]",
"PositionCartesian"=>"cart2lin",
"PositionDet"=>"charinfo(W)[:positionDet]",
"PositionId"=>"charinfo(W)[:positionId]",
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
"Rank"=>"rank",
"RankSymbol"=>"ranksymbol",
"ReducedCoxeterWord(W,w)"=>"word(W,W(w...))",
"ReducedExpressions(W,w)"=>"words(W,w)",
"ReducedInRightCoset(W,w)"=>"reduced(W,w)",
"ReducedRightCosetRepresentatives(W,H)"=>"reduced(H,W)",
"Reflection"=>"refls or reflectionmat",
"ReflectionCharacter"=>"reflchar",
#ReflectionCharValue
#ReflectionCoset
"ReflectionDegrees(W)"=>"degrees(W)",
"ReflectionCoDegrees(W)"=>"codegrees(W)",
"ReflectionEigenvalues"=>"refleigen",
"ReflectionGroup"=>"reflection_group",
"ReflectionLength(W,w)"=>"reflength(W,w)",
#ReflectionWord
"ReflectionName(W)"=>"repr(W;context=:limit=>true)",
"Reflections"=>"refls",
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
"RepresentativeOperation"=>"transporting_elt",
"RepresentativeRowColPermutation"=>"Perm_rowcol",
"Restricted"=>"restricted",
"RestrictedPartitions"=>"partitions",
"RestrictedPerm(p,d)"=>"restricted(p,d)",
"Reversed"=>"reverse",
#ReversedWord
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
"ScalarProduct"=>"scalarproduct",
"ScalMvp"=>"scalar",
#SchurElement
"SchurElements"=>"schur_elements",
"SchurFunctor"=>"schur_functor",
"SemisimpleCentralizerRepresentatives"=>"SScentralizer_reps",
"SemisimpleElement"=>"SS",
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
"SmallestMovedPoint"=>"smallest_moved_point",
"SmithNormalFormIntegerMat"=>"smith",
"SmithNormalFormIntegerMatTransforms(m)"=>"smith_transforms(m)",
"SolutionMat"=>"solutionmat",
"SolutionIntMat"=>"solutionmatInt",
"SpecialPieces"=>"special_pieces",
"Spets"=>"spets",
"SplitLevis"=>"split_levis",
"Stabilizer"=>"stabilizer",
"StandardParabolic"=>"standard_parabolic",
"StandardParabolicClass"=>"standard_parabolic_class",
"StructureRationalPointsConnectedCentre"=>"StructureRationalPointsConnectedCentre",
"SubSpets"=>"subspets",
"SubTorus"=>"SubTorus",
"Sum"=>"sum",
"SumIntersectionMat(m,n)"=>"(sum_rowspace(m,n),intersect_rowspace(m,n))",
"Symbols"=>"HasType.BDSymbols",
"SymbolsDefect(e,r,def,ct)"=>"symbols(e,r,ct,def)",
"SymmetricDifference"=>"symdiff",
"SymmetricPower"=>"symmetric_power",
"Tableaux"=>"tableaux",
"M.ToOrdinary(i)"=>"B(M,i)",
"Torus"=>"torus",
"TorusOrder"=>"torus_order",
"TraceMat"=>"tr",
"TransitiveClosure"=>"transitive_closure",
"Transporter"=>"transporter",
"TransposedMat"=>"transpose or permutedims",
"Transversals"=>"related to transversals",
"TriangulizeMat"=>"echelon!",
"Twistings"=>"twistings",
"TwoTree(m)"=>"twotree(m)",
"UnipotentAbelianPart"=>"abelianpart",
"UnipotentCharacter"=>"UniChar",
"UnipotentCharacters"=>"UnipotentCharacters",
"UnipotentClasses"=>"UnipotentClasses",
"UnipotentDecompose"=>"decompose",
"UnipotentDegrees(W,q)"=>"degrees(UnipotentCharacters(W),q)",
"UnipotentGroup"=>"UnipotentGroup",
"UnorderedTuples"=>"multisets",
"Valuation(p)"=>"valuation(p)",
"Value(p,x)"=>"p(x)",
"WeightInfo"=>"weightinfo",
"WGraph"=>"Wgraph",
"WGraphToRepresentation"=>"WGraphToRepresentation",
)

function gap(s)
  s=Regex(s,"i")
  kk=filter(x->occursin(s,x),keys(ChevieDict))
  if isempty(kk) 
    println("no match")
    return
  end
  pad=maximum(length(k) for k in kk)+2
  print(join(sort([rpad(k,pad)*"=>  "*ChevieDict[k]*"\n" for k in kk]),""))
end

function fixdoc()
  s=read("index.md",String)
  pad=maximum(length(k) for k in keys(ChevieDict))+2
  u=[rpad(k,pad)*v*"\n" for (k,v) in ChevieDict]
  u=join(sort(u,by=lowercase),"")
  s=replace(s,r"The dictionary from GAP3/Chevie is as follows:\n```(.*)```"s=>
           "The dictionary from GAP3/Chevie is as follows:\n```\n"*u*"```")
  open("index.md","w")do f
    write(f,s)
  end
end
