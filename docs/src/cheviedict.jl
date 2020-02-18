const ChevieDict=Dict(
"AlgebraicCentre"=>"algebraic_centre",
"AsReflection"=>"reflection",
"AsFraction"=>"fraction",
"AsWord"=>"word",
#BadPrimes
"BipartiteDecomposition"=>"bipartite_decomposition",
"Braid"=>"BraidMonoid",
"BraidMonoid"=>"BraidMonoid",
"BraidRelations"=>"braid_relations",
# BrieskonNormalForm
"Bruhat"=>"bruhatless",
"BruhatPoset"=>"Poset",
"BruhatSmaller"=>"bruhatless",
"Catalan"=>"catalan",
"CartanMat(\"A\",5)"=>"cartan(:A,5)",
"CentralizerGenerators"=>"centralizer_generators"
"CharNames"=>"charnames"
#CharParams
"ChevieClassInfo"=>"classinfo",
"ChevieCharInfo"=>"charinfo",
"Coefficient(p,i)"=>"p[i]",
"ComplexReflectionGroup"=>"ComplexReflectionGroup",
#ConjugacySet
"CoxeterElements(W[,l])"=>"elements(W[,l])",
"CoxeterGroup(\"A\",5)"=>"coxgroup(:A,5)",
#CoxeterGroupByCoxeterMatrix
#CoxeterGroupByCartanMatrix
"CoxeterGroupHyperoctaedralGroup(n)"=>"CoxHyperoctaedral(n)",
"CoxeterGroupSymmetricGroup(n)"=>"CoxSym(n)",
"CoxeterLength(W,w)"=>"length(W,w)",
"CoxeterMatrix"=>"coxmat",
"CoxeterWord(W,w)"=>"word(W,w)",
"CoxeterWords(W[,l])"=>"word.(Ref(W),elements(W[,l])",
"CyclotomicPolynomial(R,i)"=>"cyclotomic_polynomial(i)",
"Degree(p)"=>"degree(p)",
#DescribeInvolution
#DetPerm
#Discriminant
#Dual
"DualBraid"=>"DualBraidMonoid",
"DualBraidMonoid"=>"DualBraidMonoid",
"ElementWithInversions(W,l)"=>"with_inversions(W,l)",
#EltBraid
"EltWord(W,w)"=>"W(w...)",
"FakeDegree"=>"fakedegree",
"FakeDegrees"=>"fakedegrees",
"FiniteCoxeterTypeFromCartanMat(m)"=>"type_cartan(m)",
"FirstLeftDescending(W,w)"=>"firstleftdescent(W,w)",
"ForEachCoxeterWord(W,f)"=>"for w in W f(word(W,w)) end",
"ForEachElement(W,f)"=>"for w in W f(w) end",
#FundamentalGroup
"GarsideAlpha"=>"α",
"GarsideWords"=>"elements",
#GenericOrder
#GoodCoxeterWord
"HeckeCentralMonomials"=>"central_monomials",
#HighestPowerFakeDegrees
#HighestPowerGenericDegrees
#HighestShortRoot
"HyperplaneOrbits"=>"hyperplane_orbits",
"IndependentRoots"=>"independent_roots",
"InductionTable"=>"InductionTable",
#IntermediateGroup
#InvariantForm
#Invariants
"Inversions"=>"inversions",
#IsIsolated
"IsLeftDescending(W,w,i)"=>"isleftdescent(W,w,i)",
#IsQuasiIsolated
#IsomorphismType
#jInductionTable
#JInductionTable
"LeadingCoefficient(p)"=>"p[end]",
"LeftDescentSet(W,w)"=>"leftdescents(W,w)",
"LeftDivisorsSimple"=>"left_divisors",
"LeftGcd"=>"leftgcd",
#LeftLcm
"ListPerm(p)"=>"vec(p)",
"LongestCoxeterElement(W)"=>"longest(W)",
"LongestCoxeterWord(W)"=>"word(W,longest(W))",
#LowestPowerFakeDegrees
#LowestPowerGenericDegrees
"MatXPerm(W,p)"=>"matX(W,p)",
#MatYPerm
"OnTuples(l,p)"=>"l.^p",
"ParabolicRepresentatives"=>"parabolic_representatives",
#ParabolicSubgroups
#PermCosetsSubgroup
"PermListList(l1,l2)"=>"Perm(l1,l2)",
"PermList(v)"=>"Perm(v)",
#PermMatX
#PermMatY
"Permuted(v,p)"=>"permuted(v,p)",
"PositionClass"=>"position_class",
#PositionDet
#Presentation
"PrintDiagram(W)"=>"Diagram(W)",
"ProportionalityCoefficient(v,w)"=>"ratio(v,w)",
#QuasiIsolatedRepresentatives
"Rank"=>"rank",
#StandardParabolicClass
"ReducedCoxeterWord(W,w)"=>"word(W,W(w...))"
"ReducedExpressions(W,w)"=>"words(W,w)",
"ReducedInRightCoset(W,w)"=>"reduced(W,w)",
"ReducedRightCosetRepresentatives(W,H)"=>"reduced(H,W)",
"Reflection"=>"reflection",
"ReflectionCharacter"=>"reflchar",
#ReflectionCharValue
"ReflectionDegrees(W)"=>"degrees(W)",
"ReflectionCoDegrees(W)"=>"codegrees(W)",
"ReflectionEigenvalues"=>"refleigen",
"ReflectionLength(W,w)"=>"reflength(W,w)",
#ReflectionWord
#ReflectionName
"Reflections"=>"reflections",
"ReflectionSubgroup"=>"reflection_subgroup",
"ReflectionType"=>"refltype",
"Representations"=>"representations",
"RepresentativeConjugation"=>"representative_operation",
"RestrictedPerm(p,d)"=>"restricted(p,d)",
#ReversedWord
"RightDescentSet(W,w)"=>"rightdescents(W,w)",
"RightGcd"=>"rightgcd",
"RightLcm"=>"rightlcm"
"RootDatum"=>"rootdatum",
"RootsCartan(m)"=>"roots(m)",
#SemisimpleCentralizerRepresentatives
#SemisimpleElement
"SemisimpleRank(W)"=>"coxrank(W)",
"SemisimpleRank"=>"semisimplerank",
#SemisimpleSubgroup
"ShrinkGarsideGeneratingSet"=>"shrink",
"Size(W)"=>"length(W)",
"StandardParabolic"=>"standard_parabolic",
"StandardParabolicClass"=>"standard_parabolic_class",
#SubTorus
"TwoTree(m)"=>"twotree(m)",
#Torus
#TorusOrder
"Valuation(p)"=>"valuation(p)",
"Value(p,x)"=>"p(x)",
#WeightInfo
"W.matgens[i]"=>"matX(W,i)",
"W.N"=>"nref(W)",
"W.orbitRepresentative[i]"=>"simple_representative(W,i)",
"W.orbitRepresentativeElement"=>"simple_conjugating_element(W,i)",
)

function chevie(s)
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
