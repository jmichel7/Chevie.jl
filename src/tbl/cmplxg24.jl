#  tbl/cmplxg24.jl            CHEVIE library       Gunter Malle, Jean Michel
#  Copyright (C) 1998-2018

chevieset(:G24, :simpleroots,
          [1 root(-7)  0;
           1 -root(-7)  0;
           (-1-root(-7))//2  (-7-3*root(-7))//6  -4//3])
chevieset(:G24, :simplecoroots,
          [1  (-3*root(-7))//7  0;
           1  (3*root(-7))//7  0;
           (-1+root(-7))//2  (-7+3*root(-7))//14  -1//2]//2)

chevieset(:G24,:ordergens,fill(2,3))

chevieset(:G24,:BraidRelations, [[[1, 2, 1], [2, 1, 2]],[[1, 3, 1], [3, 1, 3]],
  [[3, 2, 3, 2], [2, 3, 2, 3]],
  [[2, 3, 1, 2, 3, 1, 2, 3, 1], [3, 2, 3, 1, 2, 3, 1, 2, 3]]])
# 14/1/2009 JM: changed generators to get presentation no.2 of  Bessis-Michel

# of interest are the alternate presentations  (1,3 of Bessis-Michel)
chevieset(:G24, :AltPres, [Dict{Symbol, Any}(:gens => [[1], [2, 3, -2], [2]],
:rels => [[[1, 2, 1, 2], [2, 1, 2, 1]], [[2, 3, 2, 3], [3, 2, 3, 2]],
      [[1, 3, 1], [3, 1, 3]], [[2, 1, 2, 3, 1, 2, 3], [1, 2, 3, 1, 2, 3, 1]]]),
                  Dict{Symbol, Any}(:gens => [[2], [3], [-3, -2, 1, 2, 3]],
:rels => [[[1, 3, 1, 3], [3, 1, 3, 1]], [[3, 2, 3, 2], [2, 3, 2, 3]],
          [[1, 2, 1, 2], [2, 1, 2, 1]], [[2, 3, 1, 2, 3, 1, 2],
  [1, 2, 3, 1, 2, 3, 1]], [[2, 3, 1, 2, 3, 1, 2], [3, 1, 2, 3, 1, 2, 3]]])])
#  For the last pres., one of the 3rd or 5th relation may be omitted

chevieset(:G24, :ReflectionDegrees, [4, 6, 14])

chevieset(:G24, :NrConjugacyClasses, 12)

chevieset(:G24, :ParabolicRepresentatives, s->
       [[[]], [[1]], [[1, 2], [2, 3]], [1:3]][s+1])

chevieset(:G24, :ClassNames, [".", "1", "23", "13", "ccc", "c", "2323", "cc1",
                     "cccccc", "cc", "ccc12", "z"])
# z=c^7 where c=123 is coxeter --- in the braid group (cc1)^6=zz

chevieset(:G24, :WordsClassRepresentatives, map(x->
 collect(replace(x,"."=>"","z"=>"123"^7,"c"=>"123")).-'0',
 chevieget(:G24, :ClassNames)))

chevieset(:G24, :PowerMaps, [nothing, [1, 1, 7, 4, 9, 10, 1, 4, 9, 10, 7, 1], [1, 2, 3, 1, 6, 5, 7, 12, 10, 9, 11, 12], nothing, [1, 2, 3, 4, 6, 5, 7, 8, 10, 9, 11, 12], nothing, [1, 2, 3, 4, 12, 12, 7, 8, 1, 1, 11, 12], nothing, nothing, nothing, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], nothing, [1, 2, 3, 4, 6, 5, 7, 8, 10, 9, 11, 12]])

chevieset(:G24, :ClassInfo, Dict{Symbol, Any}(
  :classtext=>chevieget(:G24, :WordsClassRepresentatives),
  :classnames => chevieget(:G24, :ClassNames),
  :classparams => chevieget(:G24, :ClassNames),
  :orders => [1, 2, 4, 3, 14, 14, 2, 6, 7, 7, 4, 2],
  :classes => [1, 21, 42, 56, 24, 24, 21, 56, 24, 24, 42, 1],
  :powermaps=>chevieget(:G24,:PowerMaps)))

chevieset(:G24, :CharInfo, function ()
  res=Dict{Symbol,Any}(:charparams=>[[1,0],[1,21],[3,8],[3,1],[3,10],[3,3],
                                     [6,2],[6,9],[7,6],[7,3],[8,4],[8,5]],
   :hgal => perm"(11,12)", :extRefl => [1, 4, 5, 2])
  res[:b]=map(x->x[2],res[:charparams])
  res[:charnames]=exceptioCharName.(res[:charparams])
  res
end)

chevieset(:G24, :CycPolSchurElements, 
  [[1, 0, 2, 2, 2, 3, 4, 6, 7, 14], 
   [1, -21, 2, 2, 2, 3, 4, 6, 7, 14], 
   [2*root(-7), -8, 2, 2, 2, 1//7, 2//7, 4//7, 3//14, 5//14, 13//14], 
   [-2*root(-7), -1, 2, 2, 2, 3//7, 5//7, 6//7, 1//14, 9//14, 11//14], 
   [-2*root(-7), -8, 2, 2, 2, 3//7, 5//7, 6//7, 1//14, 9//14, 11//14], 
   [2*root(-7), -1, 2, 2, 2, 1//7, 2//7, 4//7, 3//14, 5//14, 13//14], 
   [2, -1, 2, 4, 7], 
   [2, -8, 2, 4, 7], 
   [1, -6, 2, 2, 2, 3, 4, 6], 
   [1, -3, 2, 2, 2, 3, 4, 6], 
   [2, -4, 3, 7], 
   [2, -4, 3, 7]])

chevieset(:G24, :sparseFakeDegrees, [[1, 0], [1, 21], [1, 8, 1, 16, 1, 18], 
  [1, 1, 1, 9, 1, 11], [1, 10, 1, 12, 1, 20], [1, 3, 1, 5, 1, 13], 
  [1, 2, 1, 4, 1, 6, 1, 8, 1, 10, 1, 12], 
  [1, 9, 1, 11, 1, 13, 1, 15, 1, 17, 1, 19], 
  [1, 6, 1, 8, 1, 10, 1, 12, 1, 14, 1, 16, 1, 18], 
  [1, 3, 1, 5, 1, 7, 1, 9, 1, 11, 1, 13, 1, 15], 
  [1, 4, 1, 6, 1, 8, 1, 10, 1, 12, 2, 14, 1, 16], 
  [1, 5, 2, 7, 1, 9, 1, 11, 1, 13, 1, 15, 1, 17]])

chevieset(:G24,:root,function(para,rootpara)
  if ismissing(rootpara[1]) root(-prod(para[1])) else rootpara[1] end
end)

# Computed JM oct. 2005
chevieset(:G24, :HeckeCharTable, function (para,rootpara)
  f1(r)=[1, r, r^2, r^2, r^9, r^3, r^4, r^7, r^18, r^6, r^11, r^21]
  f3(p,r,b)=[3,2p+r,p^2,p*r+p^2,b*p^6*r^3,(-b-1)*p^2*r,
             -2*p^2*r^2+p^4,0,b*p^12*r^6,(-b-1)*p^4*r^2,-p^7*r^4,3*p^14*r^7]
  f6(r,p)=[6,2p+4r,2*p*r+2*r^2,2*p*r+2*r^2,p^3*r^6,p*r^2,2*r^4,0,-p^6*r^12,
           -p^2*r^4,0,-6*p^7*r^14]
  f7(p,r)=[7,4p+3r,2*p*r+p^2,2*p*r+2*p^2+r^2,0,0,-2*p^2*r^2+p^4,p^4*r^3,0,0,
           -p^6*r^5,7*p^12*r^9]
  f8(p,r,u)=[8, 4p+4r, 2*p*r+p^2+r^2, 3*p*r+p^2+r^2, u*p^4*r^4, -p*r*u,
   -2*p^2*r^2+p^4+r^4, p^3*r^3*u, -p^9*r^9, -p^3*r^3, 0, 8*p^10*r^10*u]
  r,p=para[1]
  b=E(7,3)+E(7,5)+E(7,6)
  u=chevieget(:G24,:root)(para,rootpara)
  tbl=Dict(:identifier=>"H(G24)",:name=>"H(G24)",:size=>336,:order=>336,
    :powermap=>chevieget(:G24,:PowerMaps),
    :irreducibles=>toM([f1(r),f1(p),
    f3(p,r,b), f3(r,p,b), f3(p,r,conj(b)), f3(r,p,conj(b)),
    f6(r,p),f6(p,r),f7(p,r),f7(r,p),f8(p,r,u),f8(p,r,-u)]),
           :galomorphisms=>Group(perm"(5,6)(9,10)"))
  merge!(tbl, chevieget(:G24, :ClassInfo))
  merge!(tbl, chevieget(:G24, :CharInfo)())
  tbl[:centralizers]=div.(tbl[:size],tbl[:classes])
  tbl
end)

chevieset(:G24,:galomorphisms,perm"(5,6)(9,10)")

chevieset(:G24, :HeckeRepresentation, function (para, rootpara, i)
  f1(r)=map(x->[r;;],1:3)
  f3(p,r,b)=WGraph2Representation([[[2,3],[1,2],[1,3]],[[1,2,p,-r],
    [1,3,p,-r],[2,3,r*(b+1),p*b]]],[p,r]).*(p^0*r^0)
  f7(p,r)=WGraph2Representation([[[2,3],[2,3],[1,3],[1,3],[1,2],[1,2]],
    [[1,4,r,-p],[1,5,r,-p],[2,3,r,-p],[2,6,p,-r],[3,5,-p,0],[3,6,-2p,r],
     [4,5,r,0],[4,6,2r,0]]],[r,p]).*(p^0*r^0)
  f9(r,p)=WGraph2Representation([[[1],[1,2],[1,3],[2],[2],[3],[3]],
   [[1,2,0,-r],[1,3,0,p],[1,4,p,-r],[1,5,0,-r],[1,6,-p,r],[2,5,-p,0],[2,7,-p,r],
    [3,4,-p,0],[3,5,p,-r],[3,6,p,0],[3,7,p,0],[4,6,0,-p],[4,7,-r,p],[5,6,-r,p],
    [5,7,-r,0]]], [p,r]).*(p^0*r^0)
  function f11(x, y, v)
    expandrep(3,8,Tuple{typeof(v),Vector{Int64}}[(-v*y,[105]),(-v*y+x*y,
      [79]),(v*y,[129]),(v*y-x^2-x*y, [153]), (-v, [113, 150]), (-v-y, [156,
      160]),(v,[74,77, 185]), (-x^2, [151]), (-x*y, [9, 152]), (-x, [37, 169]),
      (x,[2,24,29,30, 55, 57, 84, 88, 89, 101, 102, 136, 137, 161, 168, 192]),
      (x+y,[28,56,190]),(-y, [132, 170]), (y, [3, 22, 82, 83, 100, 110, 111,
      138, 163, 165, 191]), (-1, [65]), (1, [68])])
  end
  r,p=para[1]
  b=E(7,3)+E(7,5)+E(7,6)
  if     i==1  f1(r)
  elseif i==2  f1(p)
  elseif i==3  f3(p,r,b)
  elseif i==4  f3(r,p,b)
  elseif i==5  f3(p,r,conj(b))
  elseif i==6  f3(r,p,conj(b))
  elseif i==7  f7(p,r)
  elseif i==8  f7(r,p)
  elseif i==9  f9(p,r)
  elseif i==10 f9(r,p)
  elseif i==11 f11(p,r,chevieget(:G24,:root)(para,rootpara))
  elseif i==12 f11(p,r,-chevieget(:G24,:root)(para,rootpara))
  end
end)

CHEVIE.families[:X7]=Family(Dict{Symbol, Any}(:name=>"X7",:fourierMat=>
   [-1 1 root(-7) root(-7) -2 -2 -2; 
    1 -1 root(-7) root(-7) 2 2 2; 
    root(-7) root(-7) root(-7) -root(-7) 0 0 0; 
    root(-7) root(-7) -root(-7) root(-7) 0 0 0; 
    -2 2 0 0 -2E(7)-2E(7,6) -2E(7,2)-2E(7,5) -2E(7,3)-2E(7,4);
    -2 2 0 0 -2E(7,2)-2E(7,5) -2E(7,3)-2E(7,4) -2E(7)-2E(7,6);
    -2 2 0 0 -2E(7,3)-2E(7,4) -2E(7)-2E(7,6) -2E(7,2)-2E(7,5)]//2root(-7),
   :eigenvalues => [1, 1, 1, -1, E(7,4), E(7,2), E(7)],
   :explanation => "mystery G24", :special => 1, :cospecial => 2))

chevieset(:G24, :UnipotentCharacters,
Dict{Symbol, Any}(
 :a => [0, 21, 8, 1, 8, 1, 1, 8, 6, 3, 4, 4, 1, 1, 1, 1, 4, 4, 8, 8, 8, 8],
 :harishChandra=>Dict[Dict(:parameterExponents=>[1, 1, 1], :cuspidalName => "",
 :charNumbers => 1:12, :relativeType =>
 TypeIrred(;series=:ST,indices=1:3,ST=24,rank=3),
 :eigenvalue => 1, :levi => Int[]),
 Dict(:parameterExponents => [7], :cuspidalName => "B_2", :charNumbers => [19,
 13], :relativeType =>
 TypeIrred(;series=:A,indices=[1],rank=1), :eigenvalue => -1, :levi => [2, 3]),
 mkcuspidal("G_{24}",17,E(4);qeig=1//2),
 mkcuspidal("G_{24}",18,E(4,3);qeig=1//2),
 mkcuspidal("G_{24}",20,E(7,3)),
 mkcuspidal("G_{24}",21,E(7,5)),
 mkcuspidal("G_{24}",22,E(7,6)),
 mkcuspidal("G_{24}",14,E(7,4)),
 mkcuspidal("G_{24}",15,E(7,2)),
 mkcuspidal("G_{24}",16,E(7))],
 :families => Family[Family(:C1,[1]),
   Family(:X7,[4, 6, 7, 13, 14, 15, 16],ennola=2,cospecial=2),
   Family(:C1,[10],ennola=-1),
   Family(Family(:TQZ)(2,-1,(1, -1)),[11, 12, 18, 17],ennola=-4,cospecial=2),
   Family(:C1,[9]), 
   conj(Family(:X7,[3, 5, 8, 19, 20, 21, 22],ennola=-2,cospecial=2)),
   Family(:C1,[2],ennola=-1)],
 :curtis=>[2,1,6,5,4,3,8,7,10,9,12,11,19,-20,-21,-22,-18,-17,13,-14,-15,-16],
 :A=>[0,21,20,13,20,13,13,20,18,15,17,17,13,13,13,13,17,17,20,20,20,20]))

# Invariants are f, g:=-Det(Hessian(f))/243 and bord(f,g)/145152
chevieset(:G24, :Invariants, [function (x, y, z)
  return -42*x^2*y*z-12*x^2*y^2+21//2*x^2*z^2-9//2*y^2*z^2-6*y^3*z+14*x^4+18//7*y^4-21//8*z^4 end,
  function(x,y,z) -1960*x^2*y*z^3+840*x^2*y^2*z^2-1120*x^2*y^3*z+1760*x^2*y^4-
    1225*x^2*z^4+525*y^2*z^4-280*y^3*z^3+3920*x^4*y*z+1120*x^4*y^2-980*x^4*z^2-     180*y^4*z^2-240*y^5*z+1568*x^6-416//7*y^6-49//2*z^6 end,
  function (x, y, z) return -857157//4*x^2*y*z^11-4847619//32*x^2*y^2*z^10+
  1596665//8*x^2*y^3*z^9+18321345//16*x^2*y^4*z^8+179046*x^2*y^5*z^7+
  576093//2*x^2*y^6*z^6-440118*x^2*y^7*z^5+1608075*x^2*y^8*z^4-
  633080*x^2*y^9*z^3+269760*x^2*y^10*z^2+48576*x^2*y^11*z+785728//7*x^2*y^12-
  1327753//128*x^2*z^12+569037//128*y^2*z^12-122451//4*y^3*z^11-
  11176655//16*x^4*y*z^9+432180*x^4*y^2*z^8-2088870*x^4*y^3*z^7-
  2922360*x^4*y^4*z^6-24696*x^4*y^5*z^5-5735940*x^4*y^6*z^4-
  4210080*x^4*y^7*z^3+2688840*x^4*y^8*z^2+148960*x^4*y^9*z-
  203456*x^4*y^10+11311111//64*x^4*z^10+2077551//64*y^4*z^10+
  684285//16*y^5*z^9+2924418*x^6*y*z^7+15047067//2*x^6*y^2*z^6-
  16696554*x^6*y^3*z^5+3755850*x^6*y^4*z^4+7721616*x^6*y^5*z^3+
  12098688*x^6*y^6*z^2-470400*x^6*y^7*z+2546880*x^6*y^8+17798613//16*x^6*z^8-
  396459//8*y^6*z^8+76734*y^7*z^7+8319465*x^8*y*z^5+432180*x^8*y^2*z^4+
  8643600*x^8*y^3*z^3+24572520*x^8*y^4*z^2-3457440*x^8*y^5*z-2511936*x^8*y^6-
  8621991//4*x^8*z^6-424809//4*y^8*z^6-114513*y^9*z^5+9008552*x^10*y*z^3-
  2304960*x^10*y^2*z^2+7222208*x^10*y^3*z-8978368*x^10*y^4+6537923*x^10*z^4-
  40392*y^10*z^4+14928*y^11*z^3-537824*x^12*y*z-153664*x^12*y^2+
  134456*x^12*z^2+92712//7*y^12*z^2+30816//7*y^13*z+1382976*x^14+
  210624//343*y^14+7203//256*z^14 end])

chevieset(:G24, :BasicDerivations,()->function (x, y, z)
  [[x, 3*y^2, 7z-9*x^2*y], [3y, 1792z, 64*x*y^2+3136*x^4],
   [7z,64*x*y^3+5376*x^2*z+3136*x^4*y,
    287//2*x*y*z-35//4*x^3*y^2+21//256*y^4-1568*x^6]]
end)

chevieset(:G24, :Discriminant, ()->
  function (x, y, z)
    18*x*y^4*z+5632*x^2*y*z^2-1024*z^3-67*x^3*y^5-4352*x^4*y^2*z-5504*x^6*y^3-
    27//3136*y^7-229376*x^7*z-114688*x^9*y
end)
