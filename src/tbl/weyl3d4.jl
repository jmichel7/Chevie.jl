# tbl/weyl3d4.jl        CHEVIE library     Frank Luebeck, Jean Michel
# Copyright (C) 1994 - 2001  The CHEVIE Team

chevieset(Symbol("3D4"),:cyclestructure, 
  [[2=>2,6=>3],[3=>6],[2=>3,6=>3],[2=>3,3=>2,6=>2],[12=>2],[3=>8],[6=>4]])

chevieset(Symbol("3D4"), :generators, [
  perm"(1,13)(3,5)(6,8)(7,9)(10,11)(15,17)(18,20)(19,21)(22,23)",
  perm"(2,14)(3,6)(5,8)(7,10)(9,11)(15,18)(17,20)(19,22)(21,23)",
  perm"(1,5)(2,6)(3,15)(4,7)(11,12)(13,17)(14,18)(16,19)(23,24)",
  perm"(3,7)(4,16)(5,9)(6,10)(8,11)(15,19)(17,21)(18,22)(20,23)"])

chevieset(Symbol("3D4"), :phi,
          perm"(1,2,4)(5,6,7)(8,10,9)(13,14,16)(17,18,19)(20,22,21)")

chevieset(Symbol("3D4"),:CartanMat,[2 0 -1 0;0 2 -1 0;-1 -1 2 -1;0 0 -1 2])

chevieset(Symbol("3D4"), :positionclass,function(x)
  findfirst(==(tally(cycletype(x*chevieget(Symbol("3D4"),:phi)))),
              chevieget(Symbol("3D4"), :cyclestructure))
end)

chevieset(Symbol("3D4"),:ClassInfo,function()
  res=Dict{Symbol, Any}(
    :classtext=>[[1],Int[],[1,2,3,1,2,3],[3],[1,3],[1,2,3,1,2,4,3,2],[1,2,3,2]],
    :classnames=>["C_3","\\tilde A_2","C_3+A_1","\\tilde A_2+A_1","F_4",
                  "\\tilde A_2+A_2","F_4(a_1)"],
    :orders=>[6,3,6,6,12,3,6],:classes=>[48,16,16,48,48,8,8])
  res[:classparams] = res[:classnames]
  res
end)

chevieset(Symbol("3D4"),:NrConjugacyClasses,7)

chevieset(Symbol("3D4"),:CharInfo,function()
  res=Dict{Symbol, Any}(:extRefl=>[1, 5, 4, 6, 2],
    :charparams=>[[Int[],[4]],[Int[],[1,1,1,1]],[Int[],[2,2]],[[1,1],[2]],
                  [[1],[3]],[[1],[1,1,1]],[[1],[2,1]]],
    :charRestrictions=>[13,4,10,5,11,3,6],:nrGroupClasses=>13,
    :b=>[0,12,4,4,1,7,3],:B=>[0,12,8,8,5,11,9])
  res[:charnames]=string_partition_tuple.(res[:charparams])
  res
end)

chevieset(Symbol("3D4"), :HeckeCharTable, function (param, sqrtparam)
  q=-param[1][1]//param[1][2]
  tbl=Dict{Symbol, Any}(:identifier => "H(3D4)", :parameter=>[q,q,q,q],
    :sqrtparameter=>[],:cartan=>chevieget(Symbol("3D4"),:CartanMat),:size=>192,
    :irreducibles=> 
    [q 1 q^6 q q^2 q^8 q^4; 
     -1 1 1 -1 1 1 1; 
     q-1 2 2*q^3 q-1 -q -q^4 -q^2; 
     0 0 q^4-2*q^3+q^2 0 -q 3*q^4 3*q^2; 
     q 1 q^5-2*q^4 -1 0 -2*q^6 2*q^3; 
     -1 1 -2*q^2+q q 0 -2*q^2 2q; 
     q-1 2 -q^4-q^2 q-1 0 2*q^4 -2*q^2]*q^0,
    :irredinfo => chevieget(Symbol("3D4"), :IrredInfo))
  merge!(tbl,chevieget(Symbol("3D4"),:ClassInfo)())
  tbl[:centralizers]=div.(tbl[:size],tbl[:classes])
  AdjustHeckeCharTable(tbl, param)
end)

chevieset(Symbol("3D4"), :PhiFactors, [1, E(3), 1, E(3, 2)])

chevieset(Symbol("3D4"), :Representation,i->chevieget(Symbol("3D4"),
  :HeckeRepresentation)([[1,-1],[1,-1],[1,-1],[1,-1]],[1,1,1,1],i))

chevieset(Symbol("3D4"), :HeckeRepresentation,function(param,sqrtparam,i)
  q=-param[1][1]//param[1][2]
  if sqrtparam[1]===nothing v=root(q)
  else v=-sqrtparam[1]//param[1][2]
  end
  if i==1     (gens=[[v^2;;],[v^2;;],[v^2;;],[v^2;;]],F=[1;;])
  elseif i==2 (gens=[[-1;;],[-1;;],[-1;;],[-1;;]],F=[1;;])
  elseif i==3 
   (gens=[[-1 v;0 v^2],[-1 v;0 v^2],[v^2 0;v -1],[-1 v;0 v^2]],F=[1 0;0 1])
  elseif i==4
  (gens=[[v^2 0 0 0 0 0; 
          0 -1 0 0 0 0; 
          0 0 v^2 0 0 0; 
          0 0 0 v^2 0 0; 
          0 0 v 0 -1 0; 
          0 0 0 v 0 -1],
         [-1 0 0 v 0 0; 
          0 v^2 0 0 0 0; 
          0 0 -1 0 0 0; 
          0 0 0 v^2 0 0; 
          0 v 0 0 -1 0; 
          0 0 0 0 0 v^2], 
         [v^2 0 0 0 0 0; 
          0 -1 0 0 v v; 
          -v 0 -1 0 v 0; 
          v 0 0 -1 0 v; 
          0 0 0 0 v^2 0; 
          0 0 0 0 0 v^2], 
         [-1 0 -v 0 0 0; 
          0 v^2 0 0 0 0; 
          0 0 v^2 0 0 0; 
          0 0 0 -1 0 0; 
          0 0 0 0 v^2 0; 
          0 v 0 0 0 -1]],
       F=[0 0 0 0 0 -1; 
          0 0 -1 0 0 0; 
          0 0 0 1 0 0; 
          0 -1 0 0 0 0; 
          1 0 0 0 0 0; 
          0 0 0 0 -1 0])
  elseif i==5
  (gens=[[-1 v 0 0;0 v^2 0 0;0 0 v^2 0;0 0 0 v^2],
         [v^2 0 0 0;0 v^2 0 0;0 v -1 0;0 0 0 v^2],
         [v^2 0 0 0;v -1 v v;0 0 v^2 0;0 0 0 v^2],
         [v^2 0 0 0;0 v^2 0 0;0 0 v^2 0;0 v 0 -1]],
       F=[0 0 1 0;0 1 0 0;0 0 0 1;1 0 0 0])
  elseif i==6
  (gens=[[-1 0 0 0;0 v^2 0 0;0 v -1 0;0 0 0 -1],
         [v^2 0 0 0;0 -1 0 0;v 0 -1 0;0 0 0 -1],
         [-1 0 v 0;0 -1 v 0;0 0 v^2 0;0 0 v -1],
         [-1 0 0 0;0 -1 0 0;0 0 -1 v;0 0 0 v^2]],
      F=[0 0 0 1;1 0 0 0;0 0 1 0;0 1 0 0])
  elseif i==7 (gens=[
         [-1 v 0 0 0 0 0 2v;
          0 v^2 0 0 0 0 0 0;
          0 v -1 0 0 0 0 0;
          0 0 0 v^2 0 0 0 0;
          0 0 0 0 v^2 0 0 0;
          0 0 0 v 0 -1 0 0;
          0 0 0 0 v 0 -1 0;
          0 0 0 0 0 0 0 v^2],
         [-1 v 0 0 0 0 2v 0;
          0 v^2 0 0 0 0 0 0;
          0 0 v^2 0 0 0 0 0;
          0 v 0 -1 0 0 0 0;
          0 0 0 0 v^2 0 0 0;
          0 0 v 0 0 -1 0 0;
          0 0 0 0 0 0 v^2 0;
          0 0 0 0 v 0 0 -1],
         [v^2 0 0 0 0 0 0 0;
          v -1 0 0 0 0 0 0;
          0 0 -1 0 0 v v 0;
          0 0 0 -1 0 v 0 v;
          0 0 0 0 -1 0 v v;
          0 0 0 0 0 v^2 0 0;
          0 0 0 0 0 0 v^2 0;
          0 0 0 0 0 0 0 v^2],
         [-1 v 0 0 0 2v 0 0;
          0 v^2 0 0 0 0 0 0;
          0 0 v^2 0 0 0 0 0;
          0 0 0 v^2 0 0 0 0;
          0 v 0 0 -1 0 0 0;
          0 0 0 0 0 v^2 0 0;
          0 0 v 0 0 0 -1 0;
          0 0 0 v 0 0 0 -1]],
       F=[1 0 0 0 0 0 0 0;
          0 1 0 0 0 0 0 0;
          0 0 0 1 0 0 0 0;
          0 0 0 0 1 0 0 0;
          0 0 1 0 0 0 0 0;
          0 0 0 0 0 0 0 1;
          0 0 0 0 0 1 0 0;
          0 0 0 0 0 0 1 0])
  end
end)

chevieset(Symbol("3D4"),:sparseFakeDegrees,[[1,0],[1,12],[1,4,1,8],
  [-1,4,2,6,-1,8],[1,1,-1,3,1,5],[1,7,-1,9,1,11],[1,3,1,9]])

chevieset(Symbol("3D4"),:UnipotentCharacters,
  Dict{Symbol,Any}(
    :harishChandra=>[Dict{Symbol,Any}(:relativeType=>
      TypeIrred(;series=:G,indices=[3,1],rank=2),:levi=>Int[],:eigenvalue=>1,
      :parameterExponents=>[1,3],:cuspidalName=>"",
      :charNumbers=>[1,2,5,6,7,3]),
     Dict{Symbol,Any}(:relativeType=>TypeIrred(;series=:A,indices=[],rank=0),
      :levi=>1:4,:eigenvalue=>1,:parameterExponents=>[],
      :cuspidalName=>"{}^3D_4[1]",:charNumbers=>[4]),
     Dict{Symbol,Any}(:relativeType=>TypeIrred(;series=:A,indices=[],rank=0),
      :levi=>1:4,:eigenvalue=>-1,:parameterExponents=>[],
      :cuspidalName=>"{}^3D_4[-1]",:charNumbers =>[8])],
   :almostHarishChandra=>[Dict{Symbol, Any}(:relativeType=>
     TypeIrred(;orbit=[TypeIrred(;series=:D,indices=1:4,rank=4)],
       twist=perm"(1,2,4)"),:levi=>Int[],:eigenvalue=>1,
     :cuspidalName=>"",:charNumbers=>1:7),
     Dict{Symbol,Any}(:relativeType=>
      TypeIrred(;orbit=[TypeIrred(;series=:A,indices=Int[],rank=0)],
                twist=Perm()),:levi=>1:4,:eigenvalue=>-1,
        :cuspidalName=>"D_4",:charNumbers=>[8])],
   :families=>[Family("C1",[1]),Family("C1",[2]),
               Family("C1",[5],Dict{Symbol,Any}(:ennola=>-1)),
               Family("C1", [6],Dict{Symbol,Any}(:ennola=>-1)),
               Family("C2",[7,4,3,8],Dict{Symbol,Any}(:ennola=>-4))],
   :a=>[0,12,3,3,1,7,3,3],
   :A=>[0,12,9,9,5,11,9,9])
)

chevieset(Symbol("3D4"), :UnipotentClasses, function(p)
  class(n)=uc[:classes][findfirst(x->x[:name]==n,uc[:classes])]
  uc=deepcopy(chevieget(:D, :UnipotentClasses)(4,p))
  for c in [["11111111", perm"(1,2,4)"],
            ["221111", perm"(1,2,3)"]]
    if p!=2 class(c[1])[:red]=spets(class(c[1])[:red],c[2]) end
  end
  c=class("3311")
  g=coxgroup(:G,2)
  c[:red]=torus(g,5)
  c[:F]=g(1,2,1,2)
  c[:AuAction]=ExtendedReflectionGroup(Group(c[:red]),longest(g))
  uc
end)
