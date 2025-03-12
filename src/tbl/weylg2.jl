#  weylg2.jl   CHEVIE library         Meinolf Geck, Jean Michel
chevieset(:G2, :CartanMat,(type_=1)->[2 -type_;-3//type_ 2])

chevieset(:G2, :ParabolicRepresentatives, s->
  chevieget(:imp,:ParabolicRepresentatives)(6,6,2,s))

# in dim 3, as in Bourbaki
chevieset(:G2, :simpleroots, [1 -1 0;-2 1 1])

chevieset(:G2, :HyperplaneRepresentatives, [1, 2])

chevieset(:G2, :ReflectionDegrees, [2, 6])

chevieset(:G2, :NrConjugacyClasses, 6)

chevieset(:G2, :CharInfo, function ()
  res=Dict{Symbol,Any}(:charparams=>[[1,0],[1,6],[1,3,1],[1,3,2],[2,1],[2,2]],
                       :extRefl=>[1,5,2],:a=>[0,6,1,1,1,1],:A=>[0,6,5,5,5,5])
  res[:b]=map(x->x[2],res[:charparams])
  res[:B]=[0,6,3,3,5,4]
  # charnames in Spaltenstein's "Sous-groupes de Borel et classes unipotentes"
  res[:spaltenstein]=["1","\\varepsilon","\\varepsilon_l","\\varepsilon_c","\\theta'","\\theta''"]
  res[:lusztig]=["1","\\varepsilon","\\varepsilon'","\\varepsilon''","\\theta'","\\theta''"]
  res[:charnames]=exceptioCharName.(res[:charparams])
  res[:charSymbols]=[[[0],[0],[0],[0],[0],[2]],
   [[0,1],[0,1],[0,1],[0,1],[0,1],[1,2]],[[0],[0],[1],2,0],[[0],[0],[1],2,1],
   [[0],[0],[0],[0],[1],[1]],[[0],[0],[0],[1],[0],[1]]]
  res
end)

chevieset(:G2,:ClassNames,["A_0","\\tilde A_1","A_1","G_2","A_2","A_1+\\tilde A_1"])

chevieset(:G2, :ClassInfo, Dict{Symbol, Any}(
  :classtext => [[], [2], [1], [1, 2], [1, 2, 1, 2], [1, 2, 1, 2, 1, 2]],
  :classnames => chevieget(:G2, :ClassNames),
  :classparams => chevieget(:G2, :ClassNames),
  :orders => [1, 2, 2, 6, 3, 2], :classes => [1, 3, 3, 2, 2, 1]))

chevieset(:G2, :PowerMaps, [nothing, [1, 1, 1, 5, 5, 1], [1, 2, 3, 6, 1, 6]])

chevieset(:G2,:sparseFakeDegrees,[[1,0],[1,6],[1,3],[1,3],[1,1,1,5],[1,2,1,4]])

chevieset(:G2, :ClassParameter, w->chevieget(:G2,:ClassNames)[
  findfirst(x->w in x,[[Int[]],[[2],[1,2,1],[2,1,2,1,2]],
    [[1],[2,1,2],[1,2,1,2,1]],[[2,1],[1,2]],[[2,1,2,1],[1,2,1,2]],
    [[1,2,1,2,1,2]]])])

#  schur_elements, CharTable and representations of Hecke algebras
#  need the square root of the product of the parameters for 2-dimensional
#  representations.
#  If  sqrtpara are not bound, the function squv knows how to compute this
#  root only if u and v are equal or one is the cube of the other one.
chevieset(:G2, :squv, function (para, sqrtpara)
  u=prod(para[1])
  v=prod(para[2])
  if u==v u
  elseif u==v^3 -v^2
  elseif v==u^3 -u^2
  elseif sqrtpara[1]!==nothing && sqrtpara[2]!==nothing sqrtpara[1]*sqrtpara[2]
  else root(u*v)
  end
end)

chevieset(:G2, :HeckeCharTable, function (para, sqrtpara)
  x,y=para[1]
  z,t=para[2]
  one=(x*y*z*t)^0
  f1(u,v)=[1,v,u,v*u,v^2*u^2,v^3*u^3]*one
  function f2(x,y,z,t,eps)
    squv=eps*chevieget(:G2,:squv)(para,sqrtpara)
    [2,z+t,x+y,-squv,-x*y*z*t,2*squv^3]*one
  end
  tbl=Dict{Symbol, Any}(:identifier => "H(G2)", :parameter=>[[x,y],[z,t]],
    :size => 12, :powermap => chevieget(:G2, :PowerMaps),
    :irreducibles=>toM([f1(x,z),f1(y,t),f1(y,z),f1(x,t),f2(x,y,z,t,1),f2(x,y,z,t,-1)]),
    :irredinfo => chevieget(:G2, :IrredInfo))
  merge!(tbl,chevieget(:G2,:ClassInfo))
  tbl[:centralizers]=div.(tbl[:size],tbl[:classes])
  tbl
end)

chevieset(:G2, :HeckeRepresentation, function(para, sqrtpara, i)
  one = prod(para[1])^0*prod(para[2])^0
  x,y=para[1]
  z,t=para[2]
  if i==1 return [[x;;],[z;;]]*one
  elseif i==2 return [[y;;],[t;;]]*one
  elseif i==3 return [[y;;],[z;;]]*one
  elseif i==4 return [[x;;],[t;;]]*one
  else squv=chevieget(:G2,:squv)(para, sqrtpara)
    if i==6 squv=-squv end
    [[y -1;0 x],[z 0;squv+y*z+x*t t]] * one
  end
end)

chevieset(:G2, :PoincarePolynomial, function(param)
  u=-param[1][1]//param[1][2]
  v=-param[2][1]//param[2][2]
  (1+u)*(v+1)*(1+u*v+u^2*v^2)
end)

chevieset(:G2, :SchurModels, Dict{Symbol, Any}(
  :f1=>Dict{Symbol,Any}(:vcyc=>[[[1,-1,0,0],1],[[0,0,1,-1],1],[[1,-1,1,-1],3]]),
  :f2=>Dict{Symbol,Any}(:coeff=>-2,:root=>[1,-1,1,-1]//2,:factor=>[-1,1,0,0],
              :vcyc=>[[[0,0,0,0,1],3],[[0,0,-1,1,1],3]])))

chevieset(:G2, :SchurData, [
  Dict{Symbol, Any}(:name => "f1", :order => [1, 2, 3, 4]), 
  Dict{Symbol, Any}(:name => "f1", :order => [2, 1, 4, 3]),
  Dict{Symbol, Any}(:name => "f1", :order => [2, 1, 3, 4]),
  Dict{Symbol, Any}(:name => "f1", :order => [1, 2, 4, 3]),
  Dict{Symbol, Any}(:name => "f2", :order => [1, 2, 3, 4], :rootPower => -1),
  Dict{Symbol, Any}(:name => "f2", :order => [1, 2, 3, 4], :rootPower => 1)])

chevieset(:G2, :SchurElement, function (phi, para, sqrtpara)
  u=-para[1][1]//para[1][2]
  v=-para[2][1]//para[2][2]
  p=findfirst(==(phi),chevieget(:G2, :CharInfo)()[:charparams])
  if p==1 return (1+u)*(v+1)*(u^2*v^2+u*v+1)
  elseif p==2 return (1+u)*(v+1)*(u^2*v^2+u*v+1)//u^3//v^3
  elseif p==3 return (u^2+v^2+u*v)*(1+u)*(v+1)//u^3
  elseif p==4 return (u^2+v^2+u*v)*(1+u)*(v+1)//v^3
  end
  squv=chevieget(:G2, :squv)(para, sqrtpara)//para[1][2]//para[2][2]
  if p==6 squv=-squv end
  2*(u*v)^-1*(u*v+1+squv)*(u+v-squv)
end)

chevieset(:G2, :UnipotentCharacters, function ()
  Dict{Symbol, Any}(:harishChandra => 
  [Dict{Symbol,Any}(:relativeType=>TypeIrred(;series=:G,indices=1:2,rank=2),
     :levi => [], :parameterExponents => [1, 1], :charNumbers => 1:6,
     :eigenvalue => 1, :cuspidalName => ""), 
   mkcuspidal("G_2",10,E(3,2)),
   mkcuspidal("G_2",7,-1),
   mkcuspidal("G_2",9,E(3)),
   mkcuspidal("G_2",8,1)], 
  :families=>[Family(:S3,[5,6,4,3,8,7,9,10],ennola=-5),
              Family("C1", [1]), Family("C1", [2])], 
  :a => [0, 6, 1, 1, 1, 1, 1, 1, 1, 1], 
  :A => [0, 6, 5, 5, 5, 5, 5, 5, 5, 5], 
  :charSymbols=> [[[0],[0],[0],[0],[0],[2]],
     [[0,1],[0,1],[0,1],[0,1],[0,1],[1,2]],[[0],[0],[1],2,0],[[0],[0],[1],2,1],
     [[0],[0],[0],[0],[1],[1]],[[0],[0],[0],[1],[0],[1]],
     [[0,1],[0],[0,1],Int[],[0],Int[]],[[0,1],[0,1],[0],Int[],Int[],[0]],
     [[0,1],[0],[0],[0,1],Int[],Int[]],[[0,1],[0,1],Int[],[0],[0],Int[]]])
end)

chevieset(:G2, :Invariants, [
  (x,y)->-3x*y+3x^2+y^2, 
  (x,y)->x^2*y^4-6x^3*y^3+13x^4*y^2-12x^5*y+4x^6])

chevieset(:G2, :Discriminant,()->function(x,y) return 4*x^3*y-27*y^2 end)

chevieset(:G2, :UnipotentClasses, function(p,cartantype) # cartantype not used
  if p==0 p=1 end
  Z(n)=crg(n,1,1)
  uc=Dict{Symbol,Any}(:classes =>[
    Dict(:name=>"1",:succ=>["A1"],:dynkin=>[0,0],:balacarter=>Int[],
         :red=>coxgroup(:G,2),:rep=>Int[]), 
    Dict(:name=>"A_1",:succ=>["~A1"],:dynkin=>[1,0],:balacarter=>[1],
         :red=> Z(2),:rep=>[6]), 
    Dict(:name=>"\\tilde A_1",:succ=>["G2(a1)"],:dynkin=>[0,1],:balacarter=>[2],
         :red=>Z(2-div(gcd(p,3)-1,2)),:rep=>[4]),
    Dict(:name=>"G_2(a_1)",:succ=>["G2"],:dynkin=>[2,0],:balacarter=>[1,-2],
         :Au=>coxgroup(:A,2-div(gcd(p,3)-1,2)),:rep=>[1,5]), 
    Dict(:name=>"G_2",:succ=>String[],:dynkin=>[2,2],:Au=>Z(gcd(p,6)),
         :balacarter=>[1,2],:rep=>[2,1])],
  :springerSeries=>[
     Dict{Symbol,Any}(:relgroup=>coxgroup(:G,2),:levi=>Int[],:Z=>Int[],
                      :locsys=>[[5,1],[1,1],[4,2],[2,1],[4,3],[3,1]]), 
     Dict{Symbol, Any}(:relgroup=>coxgroup(),:levi=>[1,2],:Z=>Int[],
       :locsys=>[[4,1]],:hc=>5)])# Fourier transform of 8th unip. character
  if p==2
    uc[:springerSeries][1][:locsys][1]=[5,2]
    push!(uc[:springerSeries], 
     Dict{Symbol, Any}(:relgroup=>coxgroup(),:levi=>[1,2],:Z=>Int[],
                       :locsys=>[[5,1]],:hc=>3))
  elseif p==3
    push!(uc[:classes],Dict{Symbol, Any}(:name=>"(\\tilde A_1)_3",
     :succ=>["~A1"],:dimBu=>3,:red=>Z(2),:Au=>coxgroup(),:rep=>[4,6]))
    push!(uc[:classes][1][:succ], "(~A1)3")
    uc[:classes][3][:dimBu]=2
    delete!(uc[:classes][3],:dynkin)
    uc[:springerSeries][1][:locsys][[3,5]]=[[6,1],[4,2]]
    for c in [2, 3]
      push!(uc[:springerSeries],Dict{Symbol, Any}(:relgroup=>coxgroup(),
        :levi=>[1,2],:Z=>Int[],:locsys=>[[5,c]],:hc=>2c-2))
    end
  end
  uc[:orderClasses]=map(c->map(n->findfirst(c->Ucl.nameclass(c)==n,
                     uc[:classes]),c[:succ]),uc[:classes])
  for c in uc[:classes]
    delete!(c,:succ)
    if !(haskey(c,:red)) c[:red]=Z(1) end
    if !(haskey(c,:Au)) c[:Au]=Z(1) end
    c[:AuAction]=ExtendedReflectionGroup(c[:red],map(x->
       Matrix(1I,rank(c[:red]),rank(c[:red])), 1:semisimplerank(c[:Au])))
  end
  uc
end)

chevieset(:G2, :KLeftCellRepresentatives, [
  Dict{Symbol,Any}(:character=>[1],:duflo=>[1,2],:reps=>Vector{Int}[]), 
  Dict{Symbol,Any}(:character=>[2],:duflo=>[7,8],:reps=>Vector{Int}[]),
  Dict{Symbol,Any}(:character=>[3,5,6],:duflo=>[5,8],:reps=>[[6,10],[12,3]]),
  Dict{Symbol,Any}(:character=>[4,5,6],:duflo=>[7,3],:reps=>[[5,10],[12,4]])])
