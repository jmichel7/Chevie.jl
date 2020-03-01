const ChevieDict=Dict(
#AbelianGenerators
#Affine
#AffineRootAction
"AlgebraicCentre"=>"algebraic_centre",
#AlmostCharacter
"Arrangements"=>"arrangements",
"AsReflection"=>"reflection",
"AsFraction"=>"fraction",
"AsRootOfUnity"=>"Root1",
"AssociatedPartition"=>"conjugate_partition",
#AsymptoticAlgebra
"AsWord"=>"word",
#BadPrimes
"BetaSet"=>"βset",
"BigCellDecomposition"=>"bigcell_decomposition",
"Binomial"=>"binomial",
"BipartiteDecomposition"=>"bipartite_decomposition",
#BlocksMat
"Braid"=>"BraidMonoid",
"BraidMonoid"=>"BraidMonoid",
"BraidRelations"=>"braid_relations",
#BrieskornNormalForm
"Bruhat"=>"bruhatless",
"BruhatPoset"=>"Poset",
"BruhatSmaller"=>"bruhatless",
"Catalan"=>"catalan",
"CartanMat(\"A\",5)"=>"cartan(:A,5)",
"CartanMatFromCoxeterMatrix"=>"cartan",
"CentralizerGenerators"=>"centralizer_generators",
#CharName
"CharNames"=>"charnames",
#CharParams
"CharRepresentationWords"=>"traces_words_mats",
"ChevieClassInfo"=>"classinfo",
"ChevieCharInfo"=>"charinfo",
#ClassTypes
"Coefficient(p,i)"=>"p[i]",
"Combinations"=>"combinations",
"ComplexConjugate"=>"conj",
"ComplexReflectionGroup"=>"ComplexReflectionGroup",
"Compositions"=>"compositions",
"ConjugacySet(b[,F][,type])"=>"conjcat(b[,type[,F]]).obj",
"ConjugatePartition"=>"conjugate_partition",
#CoxeterCoset
#CoxeterSubCoset
"CoxeterElements(W[,l])"=>"elements(W[,l])",
"CoxeterGroup(\"A\",5)"=>"coxgroup(:A,5)",
#CoxeterGroupByCoxeterMatrix
#CoxeterGroupByCartanMatrix
"CoxeterGroupHyperoctaedralGroup(n)"=>"CoxHyperoctaedral(n)",
"CoxeterGroupSymmetricGroup(n)"=>"CoxSym(n)",
"CoxeterLength(W,w)"=>"length(W,w)",
"CoxeterMatrix"=>"coxmat",
"CoxeterMatrixFromCartanMat"=>"coxmat",
"CoxeterWord(W,w)"=>"word(W,w)",
"CoxeterWords(W[,l])"=>"word.(Ref(W),elements(W[,l])",
#CuspidalUnipotentCharacters
"CyclotomicPolynomial(R,i)"=>"cyclotomic_polynomial(i)",
"CycPol"=>"CycPol",
"CycPolFakeDegreeSymbol"=>"fegsymbol",
"CycPolGenericDegreeSymbol"=>"gendeg_symbol",
#CycPolUnipotentDegrees
"DecomposedMat"=>"diagblocks",
"DefectSymbol"=>"defectsymbol",
"Degree(p)"=>"degree(p)",
#DeligneLusztigCharacter
#DeligneLusztigLefschetz
"DescribeInvolution"=>"describe_involution",
#DetPerm
#Dictionary
#DifferenceMultiSet
#Discriminant
"Dominates"=>"dominates",
#DrinfeldDouble
#Dual
"DualBraid"=>"DualBraidMonoid",
"DualBraidMonoid"=>"DualBraidMonoid",
"EigenspaceProjector"=>"eigenspace_projector",
#EigenvaluesMat
"ElementWithInversions(W,l)"=>"with_inversions(W,l)",
"EltBraid"=>"image",
"EltWord(W,w)"=>"W(w...)",
"ExteriorPower"=>"exterior_power",
#FactorizedSchurElement
#FactorizedSchurElements
"FakeDegree"=>"fakedegree",
"FakeDegrees"=>"fakedegrees",
"FamiliesClassical"=>"FamiliesClassical",
"Family"=>"Family",
"FamilyImprimitive"=>"family_imprimitive",
"FiniteCoxeterTypeFromCartanMat(m)"=>"type_cartan(m)",
"FirstLeftDescending(W,w)"=>"firstleftdescent(W,w)",
"ForEachCoxeterWord(W,f)"=>"for w in W f(word(W,w)) end",
"ForEachElement(W,f)"=>"for w in W f(w) end",
#FormatTable
#Frobenius
"FundamentalGroup"=>"fundamental_group",
#FusionAlgebra
"GarsideAlpha"=>"α",
"GarsideWords"=>"elements",
"GcdPartitions"=>"gcd_partitions",
#GenericDegrees
#GenericOrder
#GetRoot
#GoodCoxeterWord
#GraphAutomorphisms
"Hasse"=>"hasse",
"HeckeCharValues"=>"char_values",
#HeckeCharValuesGood
"HeckeCentralMonomials"=>"central_monomials",
"HeckeClassPolynomials"=>"class_polynomials",
#HeckeReflectionRepresentation
#HeckeSubAlgebra
#HighestPowerFakeDegrees
"HighestPowerFakeDegreeSymbol"=>"degree_feg_symbol",
#HighestPowerGenericDegrees
"HighestPowerGenericDegreeSymbol"=>"degree_gendeg_symbol",
#HighestShortRoot
"KazhdanLusztigPolynomial"=>"KLPol",
"HyperplaneOrbits"=>"hyperplane_orbits",
"ICCTable"=>"ICCTable",
"Incidence"=>"incidence",
#IndependentLines
"IndependentRoots"=>"independent_roots",
"InducedLinearForm"=>"induced_linear_form",
"InductionTable"=>"InductionTable",
#Inherit
#IntermediateGroup
"IntListToString"=>"joindigits",
#InvariantForm
#Invariants
"Inversions"=>"inversions",
"IsCycPol(p)"=>"p isa CycPol",
"IsFamily(f)"=>"f isa Family",
#IsIsolated
"IsJoinLattice"=>"is_join_lattice",
"IsMeetLattice"=>"is_meet_lattice",
"IsLeftDescending(W,w,i)"=>"isleftdescent(W,w,i)",
#IsNormalizing
#IsQuasiIsolated
#IsomorphismType
#IsUnipotentElement
#jInductionTable
#JInductionTable
"Join"=>"join",
"LcmPartitions"=>"lcm_partitions",
"LeadingCoefficient(p)"=>"p[end]",
"LeftCell"=>"LeftCell",
"LeftCells"=>"LeftCells",
"LeftDescentSet(W,w)"=>"leftdescents(W,w)",
"LeftDivisorsSimple"=>"left_divisors",
"LeftGcd"=>"leftgcd",
#LeftLcm
"LinearExtension"=>"linear_extension",
"ListPerm(p)"=>"vec(p)",
"LongestCoxeterElement(W)"=>"longest(W)",
"LongestCoxeterWord(W)"=>"word(W,longest(W))",
#LowestPowerFakeDegrees
"LowestPowerFakeDegreeSymbol"=>"valuation_feg_symbol",
#LowestPowerGenericDegrees
"LowestPowerGenericDegreeSymbol"=>"valuation_gendeg_symbol",
#Lusztigaw
#LusztigAw
#LusztigInduction
#LusztigInductionTable
#LusztigRestriction
#MatStab
"MatXPerm(W,p)"=>"matX(W,p)",
#MatYPerm
#NrDrinfeldDouble
#OnFamily
#OnMatrices
"OnTuples(l,p)"=>"l.^p",
"ParabolicRepresentatives"=>"parabolic_representatives",
#ParabolicSubgroups
"PartBeta"=>"partβ",
"Partition"=>"partition",
"Partitions"=>"partitions",
"PartitionTuples"=>"partition_tuples",
#PartitionTupleToString
#PermCosetsSubgroup
"PermListList(l1,l2)"=>"Perm(l1,l2)",
"PermList(v)"=>"Perm(v)",
#PermMatMat
#PermMatX
#PermMatY
"PermutationMat(p,dim)"=>"Matrix(p,dim)",
"Permuted(v,p)"=>"v^p",
#PermutedByCols
#PoincarePolynomial
"Poset"=>"Poset",
"PositionClass"=>"position_class",
#PositionDet
#PositionId
"PositionRegularClass"=>"position_regular_class",
#Presentation
"PrintDiagram(W)"=>"Diagram(W)",
"ProportionalityCoefficient(v,w)"=>"ratio(v,w)",
#QuasiIsolatedRepresentatives
"Rank"=>"rank",
"RankSymbol"=>"ranksymbol",
"ReducedCoxeterWord(W,w)"=>"word(W,W(w...))",
"ReducedExpressions(W,w)"=>"words(W,w)",
"ReducedInRightCoset(W,w)"=>"reduced(W,w)",
"ReducedRightCosetRepresentatives(W,H)"=>"reduced(H,W)",
"Reflection"=>"reflection",
"ReflectionCharacter"=>"reflchar",
#ReflectionCharValue
#ReflectionCoset
"ReflectionDegrees(W)"=>"degrees(W)",
"ReflectionCoDegrees(W)"=>"codegrees(W)",
"ReflectionEigenvalues"=>"refleigen",
"ReflectionLength(W,w)"=>"reflength(W,w)",
#ReflectionWord
#ReflectionName
"Reflections"=>"reflections",
#ReflectionSubCoset
"ReflectionSubgroup"=>"reflection_subgroup",
"ReflectionType"=>"refltype",
"RegularEigenvalues"=>"regular_eigenvalues",
"RelativeDegrees"=>"relative_degrees",
#Replace
"Representations"=>"representations",
"RepresentativeConjugation"=>"representative_operation",
#RepresentativeDiagonalConjugaction
#RepresentativeRowColPermutation
"Restricted"=>"restricted",
"RestrictedPerm(p,d)"=>"restricted(p,d)",
"Reversed"=>"reverse",
#ReversedWord
"RightDescentSet(W,w)"=>"rightdescents(W,w)",
"RightGcd"=>"rightgcd",
"RightLcm"=>"rightlcm",
"RootDatum"=>"rootdatum",
"RootsCartan(m)"=>"roots(m)",
"Rotation(v,i)"=>"circshift(v,-i)",
"Rotations(v)"=>"circshift.(Ref(a),length(a):-1:1)",
#SchurElement
"SchurElements"=>"schur_elements",
#SchurFunctor
#SemisimpleCentralizerRepresentatives
#SemisimpleElement
"SemisimpleRank(W)"=>"coxrank(W)",
"SemisimpleRank"=>"semisimplerank",
#SemisimpleSubgroup
"ShiftBeta"=>"shiftβ",
"ShrinkGarsideGeneratingSet"=>"shrink",
#SignedMatStab
"SignedPerm"=>"SPerm",
"SignedPermListList"=>"SPerm",
#SignedPermMatMat
"Size(W)"=>"length(W)",
#SpecialPieces
"Spets"=>"spets",
"SplitLevis"=>"split_levis",
"StandardParabolic"=>"standard_parabolic",
"StandardParabolicClass"=>"standard_parabolic_class",
#StructureRationalPointsConnectedCentre
#SubSpets
"SubTorus"=>"SubTorus",
#Symbols
#SymbolsDefect
#SymmetricPower
"Tableaux"=>"tableaux",
"Torus"=>"torus",
#TorusOrder
"TransitiveClosure"=>"transitive_closure",
#Transporter
#Transversals
#Twistings
"TwoTree(m)"=>"twotree(m)",
#UnipotentAbelianPart
#UnipotentCharacter
"UnipotentCharacters"=>"UnipotentCharacters",
"UnipotentClasses"=>"UnipotentClasses",
#UnipotentDecompose
"UnipotentDegrees(W,q)"=>"degrees(UnipotentCharacters(W),q)",
#UnipotentGroup
"UnorderedTuples"=>"submultiset",
"Valuation(p)"=>"valuation(p)",
"Value(p,x)"=>"p(x)",
"WeightInfo"=>"weightinfo",
#WGraph
#WGraphToRepresentation
"W.matgens[i]"=>"matX(W,i)",
"W.N"=>"nref(W)",
"W.orbitRepresentative[i]"=>"simple_representative(W,i)",
"W.orbitRepresentativeElement"=>"simple_conjugating_element(W,i)",
)

function gap(s)
  s=Regex(s,"i")
  kk=filter(x->occursin(s,x),keys(ChevieDict))
  pad=maximum(length(k) for k in kk)+2
  for k in kk
    println(rpad(k,pad),"=>  ",ChevieDict[k])
  end
end

function fixdoc()
  s=read("index.md",String)
  pad=maximum(length(k) for k in keys(ChevieDict))+2
  u=[rpad(k,pad)*"→ "*v*"\n" for (k,v) in ChevieDict]
  u=join(sort(u),"")
  s=replace(s,r"The dictionary from GAP3/Chevie is as follows:\n```(.*)```"s=>
           "The dictionary from GAP3/Chevie is as follows:\n```\n"*u*"```")
  open("index.md","w")do f
    write(f,s)
  end
end
