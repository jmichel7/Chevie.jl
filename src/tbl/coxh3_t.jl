#  tbl/coxh3.jl                 CHEVIE library      Meinolf Geck, Jean Michel
#  Copyright (C) 1992 - 2004  The CHEVIE Team
chevieset(:H3, :ReflectionDegrees, [2, 6, 10])

chevieset(:H3, :Size, 120)

# from Humphreys, "Reflection Groups and Coxeter Groups"
chevieset(:H3, :GeneratingRoots, function ()
  a=-E(5,2)-E(5,3)
  [[a,-1,a-1],[-a,1,a-1],[1,a-1,-a]]//2
end)

# The Galois-invariant model is PRG(roots,coroots) given below
chevieset(:H3, :InvariantModel, Dict{Symbol, Any}(
  :Aut=>[[1,2,1,2,3,2,1,2,1],[3],[2]],
  :roots=>[[(5-root(5))//2,-root(5),-1+root(5)],[-2root(5),1,0],[2root(5),1,0]],
  :coroots=>[[(5+root(5))//40,(-1-root(5))//4,(-3+root(5))//16],
     [-3root(5)//20,1//2,(2-root(5))//8],[3root(5)//20,1//2,(2+root(5))//8]],
  :conj=>[[(5-root(5))//2,-root(5),-1+root(5)],[-2root(5),1,0],[2root(5),1,0]]))
# usual matgens^conj=invariant matgens

chevieset(:H3, :NrConjugacyClasses, 10)

chevieset(:H3,:CartanMat,[2 -(1+root(5))//2 0;-(1+root(5))//2 2 -1;0 -1 2])

chevieset(:H3, :PowerMaps, [nothing, [1, 1, 7, 1, 5, 3, 3, 5, 7, 1], [1, 2, 7, 4, 1, 9, 3, 10, 6, 10], nothing, [1, 2, 1, 4, 5, 10, 1, 8, 10, 10], nothing, [1, 2, 7, 4, 5, 9, 3, 8, 6, 10]])

chevieset(:H3, :WordsClassRepresentatives, [Int[],[1],[1,2],[1,3],[2,3],[1,2,3],
  [1,2,1,2],[1,2,1,2,3],[1,2,1,2,3,2,1,2,3],[1,2,1,2,1,3,2,1,2,1,3,2,1,2,3]])

chevieset(:H3, :ParabolicRepresentatives, s->
  [[Int[]],[[1]],[[1,2],[1,3],[2,3]],[1:3]][s+1])

chevieset(:H3, :ClassInfo, function ()
  res=Dict{Symbol,Any}(:classtext=>chevieget(:H3,:WordsClassRepresentatives),
    :orders => [1, 2, 5, 2, 3, 10, 5, 6, 10, 2],
    :classes => [1, 15, 12, 15, 20, 12, 12, 20, 12, 1])
  res[:classnames]=joindigits.(res[:classtext])
  res[:classnames][1]="."
  res[:classparams]=res[:classnames]
  res
end)

chevieset(:H3, :CharInfo, function ()
  res=Dict(
    :charparams=>[[1,15],[1,0],[5,5],[5,2],[3,6],[3,8],[3,1],[3,3],[4,3],[4,4]],
    :gp=>["1_r'","1_r","5_r'","5_r","3_s","\\overline{3}_s","3_s'","\\overline{3}_s'","4_r'","4_r"],
    :hgal=> perm"(9,10)", :extRefl => [2, 7, 5, 1])
  res[:b]=map(x->x[2],res[:charparams])
  res[:charnames]=GAPENV.exceptioCharName.(res[:charparams])
  res
end)

chevieset(:H3, :vpolheckeirreducibles,
 [[Pol([1]), Pol([-1]), Pol([1]), Pol([1]), Pol([1]),
 Pol([-1]), Pol([1]), Pol([-1]), Pol([-1]), Pol([-1])], [Pol([1]),
 Pol([1],2), Pol([1],4), Pol([1],4), Pol([1],4), Pol([1],6), Pol([1],8),
 Pol([1],10), Pol([1],18), Pol([1],30)], [Pol([5]), Pol([-3, 2]),
 Pol([1, -1]), Pol([2, -2, 1]), Pol([1, -2]), Pol([0]), Pol([1, 0, -1]),
 Pol([1],4), Pol([0]), Pol([-5],12)], [Pol([5]), Pol([-2, 3]),
 Pol([-1, 1],2), Pol([1, -2, 2]), Pol([-2, 1],2), Pol([0]), Pol([-1, 0, 1],4),
 Pol([-1],6), Pol([0]), Pol([5],18)], Pol[Pol([3]), Pol([-2, 1]),
 Pol([1, E(5)+E(5,4)]), Pol([1, -2]), Pol([1, -1]),
 Pol([-E(5)-E(5,4)],2), Pol([1, 0, E(5,2)+E(5,3)]),
 Pol([0]), Pol([(1+root(5))//2],6), Pol([3],10)],
 Pol[Pol([3]), Pol([-2, 1]), Pol([1, E(5,2)+E(5,3)]), Pol([1, -2]),
 Pol([1, -1]), Pol([(1+root(5))//2],2), Pol([1,
 0, E(5)+E(5,4)]), Pol([0]), Pol([-E(5)-E(5,4)],6), Pol([3],10)],
 Pol[Pol([3]), Pol([-1, 2]), Pol([E(5)+E(5,4), 1],2), Pol([-2, 1],
2), Pol([-1, 1],2), Pol([E(5)+E(5,4)],4), Pol([E(5,2)+E(5,
3), 0, 1],4), Pol([0]), Pol([E(5,2)+E(5,3)],12), Pol([-3],20)],
 Pol[Pol([3]), Pol([-1, 2]), Pol([E(5,2)+E(5,3), 1],2), Pol([-2, 1],
2), Pol([-1, 1],2), Pol([E(5,2)+E(5,3)],4), Pol([E(5)+E(5,
4), 0, 1],4), Pol([0]), Pol([E(5)+E(5,4)],12), Pol([-3],20)],
 [Pol([4]), Pol([-2, 2]), Pol([-1],2), Pol([1, -2, 1]), Pol([1, -1,
 1]), Pol([1],3), Pol([-1],4), Pol([-1],5), Pol([1],9), Pol([-4],15)],
 [Pol([4]), Pol([-2, 2]), Pol([-1],2), Pol([1, -2, 1]), Pol([1, -1,
 1]), Pol([-1],3), Pol([-1],4), Pol([1],5), Pol([-1],9), Pol([4],15)]]
)

chevieset(:H3, :CycPolSchurElements, [
  [1,-15,2,2,2,3,5,6,10],
  [1,0,2,2,2,3,5,6,10],
  [1,-5,2,2,2,3,6],
  [1,-2,2,2,2,3,6],
  [(5+root(5))//2,-6,2,2,2,2//5,3//5,1//10,9//10],
  [(5-root(5))//2,-6,2,2,2,1//5,4//5,3//10,7//10],
  [(5+root(5))//2,-1,2,2,2,2//5,3//5,1//10,9//10],
  [(5-root(5))//2,-1,2,2,2,1//5,4//5,3//10,7//10],
  [2,-3,3,5],
  [2,-3,3,5]])

chevieset(:H3, :sparseFakeDegrees, 
  [[1,15],[1,0],[1,5,1,7,1,9,1,11,1,13],[1,2,1,4,1,6,1,8,1,10],[1,6,1,10,1,14],
   [1,8,1,10,1,12],[1,1,1,5,1,9],[1,3,1,5,1,7],[1,3,1,7,1,9,1,11],
   [1,4,1,6,1,8,1,12]])

chevieset(:H3, :HeckeCharTable, function (param, sqrtparam)
  a=(1+root(5))//2
  q=-param[1][1]//param[1][2]
  if sqrtparam[1]===nothing v=root(q)
  else v=-sqrtparam[1]//param[1][2]
  end
  ci=chevieget(:H3, :ClassInfo)()
  tbl=Dict{Symbol, Any}(:identifier => "H(H3)",
    :text=>"the representing matrices are those of Lusztig(1981)",
    :parameter => [q, q, q], :cartan => chevieget(:H3, :CartanMat),
    :size => 120, :order => 120, :powermap => chevieget(:H3, :PowerMaps),
    :irreducibles => map(i->map(function(j)
                res=evalpoly(q,j.c)
                if iseven(valuation(j)) res*=q^div(valuation(j),2)
                else res*=v^valuation(j)
                end
                return res
            end,i), chevieget(:H3, :vpolheckeirreducibles)),
    :irredinfo => chevieget(:H3, :IrredInfo))
  merge!(tbl, ci)
  tbl[:centralizers]=div.(tbl[:size],tbl[:classes])
  chevieget(:compat,:AdjustHeckeCharTable)(tbl, param)
  tbl
end)

chevieset(:H3,:Representation,i->
  chevieget(:H3,:HeckeRepresentation)([[1,-1],[1,-1],[1,-1]],[1,1,1],i))

# W-graphs given by Ivan Marin. see ?WGraphToRepresentation for the format.
chevieset(:H3, :WGraphs, 
  [[[[1,2,3]],[]],1, 
   [[[2],[1,2],[1,3],[1,3],[2,3]],[[-1,[[1,3],[2,4],[3,5],[4,5]]]]],3, 
   [[[1,2],[1,3],[2,3]],[[-1,[[1,2]]],[E(5,2)+E(5,3),[[2,3]]]]],
   [[[1,2],[1,3],[2,3]],[[-1,[[1,2]]],[E(5)+E(5,4),[[2,3]]]]],5,6,
   [[[1], [2], [1, 3], [2, 3]], [[1, [[1, 2, 3], [2, 3, 4], [3, 4]]]]], 9])

chevieset(:H3, :WGraph, function(i)
  gr=chevieget(:H3, :WGraphs)
  if gr[i] isa Integer DualWGraph(3, gr[gr[i]])
  else gr[i]
  end
end)

chevieset(:H3, :HeckeRepresentation, function (param, sqrtparam, i)
  if sqrtparam[1]===nothing v=root(-param[1][1]//param[1][2])
  else v=-sqrtparam[1]//param[1][2]
  end
  -param[1][2]*WGraphToRepresentation(3,chevieget(:H3, :WGraph)(i),v)
end)

chevieset(:H3, :UnipotentCharacters, function ()
  Dict{Symbol, Any}(:harishChandra=>[
    Dict(:relativeType=>TypeIrred(;series=:H,indices=1:3,rank=3),
      :levi=>[], :eigenvalue=>1, :parameterExponents=>[1, 1, 1],
      :cuspidalName=>"", :charNumbers=>1:10),
    Dict(:relativeType=>TypeIrred(;series=:A,indices=[3],rank=1), :levi=>1:2,
         :eigenvalue=>E(5, 2), :parameterExponents=>[5],
         :cuspidalName=>"I_2(5)[1,3]", :charNumbers=>[11, 13]), 
    Dict(:relativeType=>TypeIrred(;series=:A,indices=[3],rank=1), :levi=>1:2,
         :eigenvalue=>E(5, 3), :parameterExponents=>[5],
         :cuspidalName=>"I_2(5)[1,2]", :charNumbers=>[12, 14]),
    Dict(:relativeType=>TypeIrred(;series=:A,indices=[],rank=0), :levi=>1:3,
         :eigenvalue=>E(4), :qEigen=>1 // 2, :parameterExponents=>[],
         :cuspidalName=>"H_3[i]", :charNumbers=>[15]), 
    Dict(:relativeType=>TypeIrred(;series=:A,indices=[],rank=0), :levi=>1:3,
         :eigenvalue=>-(E(4)), :qEigen=>1 // 2, :parameterExponents=>[],
         :cuspidalName=>"H_3[-i]", :charNumbers=>[16])],
    :families=>[Family("C1", [2]),
       Family(CHEVIE[:families][:Dihedral](5), [7, 8, 14, 13],
              Dict{Symbol, Any}(:ennola=>-1)),
       Family("C1", [4]),
       Family(CHEVIE[:families][:TQZ](2,-1,[1,-1]), [9, 10, 16, 15],
              Dict{Symbol, Any}(:cospecial=>2, :ennola=>4)),
       Family("C1", [3], Dict{Symbol, Any}(:ennola=>-1)),
       Family(CHEVIE[:families][:Dihedral](5), [5, 6, 12, 11],
              Dict{Symbol, Any}(:ennola=>1)),
       Family("C1", [1], Dict{Symbol, Any}(:ennola=>-1))],
    :a=>[15, 0, 5, 2, 6, 6, 1, 1, 3, 3, 6, 6, 1, 1, 3, 3],
    :A=>[15, 0, 13, 10, 14, 14, 9, 9, 12, 12, 14, 14, 9, 9, 12, 12])
    end)

chevieset(:H3, :Discriminant, function ()
  function(a,b,c)131835937500*a*b^3*c-100195312500*a^2*b*c^2+395507812500*c^3-
   28369140625*a^3*b^4+1371093750*a^4*b^2*c+175781250000*b^5+
   1191796875*a^5*c^2+1162187500*a^6*b^3-74250000*a^7*b*c-22233750*a^9*b^2+
   438750*a^10*c+213700*a^12*b-829*a^15
  end
end)

chevieset(:H3, :KLeftCellRepresentatives, [
  Dict{Symbol,Any}(:character=>[2],:duflo=>[1,2,3],:reps=>""),
  Dict{Symbol,Any}(:character=>[1],:duflo=>[16,17,18],:reps=>""),
  Dict{Symbol,Any}(:character=>[3],:duflo=>[1,24,3],:reps=>""),
  Dict{Symbol,Any}(:character=>[4],:duflo=>[2,1,28],:reps=>""),
  Dict{Symbol,Any}(:character=>[6,5],:duflo=>[1,20,18],:reps=>[[7,19,24]]),
  Dict{Symbol,Any}(:character=>[8,7],:duflo=>[1,6,18],:reps=>[[9,2,27]]),
  Dict{Symbol,Any}(:character=>[10,9],:duflo=>[8,18,17],
                   :reps=>[[11,17,25],[11,27,10],[14,30,4]]),
  Dict{Symbol,Any}(:character=>[10,9],:duflo=>[13,30,8],
                   :reps=>[[10,29,5],[12,21,22],[13,22,23]])])
