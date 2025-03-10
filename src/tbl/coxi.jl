# Hand-translated part of chevie/tbl/coxi.jl
# (C) 1992-2010  Meinolf Geck & Jean Michel
# Data for W(I_m)
#
chevieset(:I,:CartanMat,function(bond,cartantype=iseven(bond) ? 1 : E(2bond)+E(2bond,-1))
  m=[2*E(1) 0;0 2]
  if bond==2 return m end
  m[1,2]=-cartantype
  m[2,1]=(2+E(bond)+E(bond,-1))//m[1,2]
  m
end)

chevieset(:I, :ReflectionName, function(bond,opt,cartantype=iseven(bond) ? 1 : E(2bond)+E(2bond,-1))
  c=xrepr(cartantype^2//(2+E(bond)+E(bond,-1));opt...)
  if cartantype==1
    if haskey(opt, :TeX) string("I_2(", bond, ")")
    elseif haskey(opt, :arg) string("\"I\",2,", bond)
    else string("I2(", bond, ")")
    end
  elseif cartantype == E(2bond) + E(2bond, -1)
    if mod(bond, 2) == 1
      if haskey(opt, :TeX) string("I_2(", bond, ")")
      elseif haskey(opt, :arg) string("\"I\",2,", bond)
      else string("I2(", bond, ")")
      end
    else
      if haskey(opt, :TeX) string("I_{\\hbox{sym}2}(", bond, ")")
      elseif haskey(opt, :arg) string("\"Isym\",2,", bond)
      else string("Isym2(", bond, ")")
      end
    end
  elseif haskey(opt, :TeX) string("I_?(",c,")(",bond,")")
  elseif haskey(opt, :arg) string("\"Isym\",2,",bond,",",c)
  else string("I?(",c,")(",bond,")")
  end
end)

chevieset(:I, :SemisimpleRank, 2)

chevieset(:I, :simpleroots, function(m)
  a=E(2m,m-1)
  b=conj(a)
  if iseven(m) r=root(m//2)
  else r=1
  end
  [1 0;r*(a+b)//2 r*(a-b)//2//E(4)]
end)

chevieset(:I,:ordergens,m->fill(2,2))

chevieset(:I,:ReflectionDegrees,m->[2,m])

chevieset(:I,:NrConjugacyClasses,m->div(m+3,2)+mod(m+1,2)*2)

chevieset(:I,:ParabolicRepresentatives, function (m, s)
  chevieget(:imp, :ParabolicRepresentatives)(m, m, 2, s)
end)

chevieset(:I, :WordsClassRepresentatives, function(m)
  r=iseven(m) ? [Int[], [1], [2]] : [Int[], [1]]
  x = [1, 2]
  for i in 1:div(m,2) 
    push!(r,copy(x))
    append!(x,[1,2])
  end
  r
end)

chevieset(:I, :ClassInfo, function (m)
  r=chevieget(:I, :WordsClassRepresentatives)(m)
  clnp=joindigits.(r)
  g1=Perm()
  i=2
  while 2i<=m+1
    g1*=Perm(i,m-i+2)
    i+=1
  end
  g2=Perm()
  i=1
  while 2i<=m
    g2*=Perm(i,m-i+1)
    i+=1
  end
  gen=[g1,g2]
  perm(l)=isempty(l) ? Perm() : prod(gen[l])
  m1=div(m,2)
  if iseven(m)
    cl=[1,m1,m1]
    append!(cl,fill(2,m1-1))
    push!(cl,1)
  else
    cl=[1,m]
    append!(cl,fill(2,m1))
  end
  Dict{Symbol, Any}(:classtext=>r,:classnames=>clnp,:classparams=>clnp,
                   :orders=>map(i->order(perm(i)),r),:classes=>cl)
end)

chevieset(:I, :HeckeCharTable, function (m, para, rootpara)
  u=-para[1][1]//para[1][2]
  v=-para[2][1]//para[2][2]
  if isodd(m) squv=u
  elseif rootpara[1]!==nothing && rootpara[2]!==nothing
      squv=rootpara[1]*rootpara[2]
  else squv=root(u*v)
  end
  ct=[[u,v]]
  if iseven(m)
    append!(ct,[[u,-u^0],[-v^0,v]])
  end
  push!(ct,[-v^0,-v^0])
  cl=chevieget(:I,:ClassInfo)(m)
  r=cl[:classtext]
  ct=map(i->map(x->prod(i[x]),r),ct)*E(1)
  append!(ct, map(1:div(m-1,2))do j
    map(1:length(r))do i
      if r[i]==Int[] 2*v^0
      elseif r[i]==[1] u-1
      elseif r[i]==[2] v-1
      else k=div(length(r[i]),2);squv^k*(E(m,k*j)+E(m,-k*j))
      end
    end
  end)
  ci=chevieget(:I, :CharInfo)(m)
  tbl=Dict{Symbol, Any}(:identifier => string("H(I2(", m, "))"), 
    :cartan=>cartan(:I,2,m),:size=>2m,:irredinfo=>map((x,y)->
   Dict{Symbol,Any}(:charparam=>x,:charname=>y),ci[:charparams],ci[:charnames]),
   :parameter=>[u,v],:powermap=>[],:irreducibles=>ct*v^0)
  merge!(tbl, cl)
  tbl[:centralizers]= div.(tbl[:size],tbl[:classes])
  AdjustHeckeCharTable(tbl, para)
end)

#  The  characters of I_2(m) are uniquely  parametrized by [d,b] where d is
#  their  degree and b  is the valuation  of the fake  degrees, except that
#  when  m is  even, there  are two  characters with [d,b]=[1,m/2]. The one
#  which  maps the generators to [1,-1]  is denoted [1,m/2,1] and the one
#  which maps them to [-1,1] is denoted [1,m/2,2].
chevieset(:I, :CharInfo, function(m)
  res=Dict{Symbol, Any}()
  charparams=[[1,0]]
  if iseven(m)
    res[:extRefl]=[1,5,4]
    m1=div(m,2)
    push!(charparams,[1,m1,1],[1,m1,2])
  else
    res[:extRefl]=[1,3,2]
  end
  push!(charparams,[1,m])
  append!(charparams,map(i->[2,i],1:div(m-1,2)))
  res[:charnames]=exceptioCharName.(charparams)
  res[:charparams]=charparams
  res[:b]=map(x->x[2], charparams)
  res[:B]=map(phi->phi[1]==1 ? phi[2] : m-phi[2], charparams)
  charSymbols=[] # need type Any for m even
  push!(charSymbols,map(i->i==m ? [2] : [0],1:m))
  if iseven(m)
    v=vcat(map(i->i==m1 ? [1] : [0],1:m1),[2,0])
    push!(charSymbols,v)
    v=copy(v);v[m1+2]=1;push!(charSymbols,v)
  end
  push!(charSymbols,map(i->i==m ? [1,2] : [0,1],1:m))
  append!(charSymbols,map(1:div(m-1,2))do l
    map(i->i==1 ? [1] : i==l+1 ? [1] : [0],1:m)
  end)
  res[:charSymbols]=charSymbols
  res[:A]=degree_gendeg_symbol.(charSymbols)
  res[:a]=valuation_gendeg_symbol.(charSymbols)
  # malleParams are the partitiontuples for the symbols
  res[:malleParams]=map(x->map(partβ,fullsymbol(x)),charSymbols)
  if iseven(m)
    res[:malleParams]=convert.(Vector{Any},res[:malleParams])
    res[:malleParams][2]=append!(res[:malleParams][2][1:m1],[2,0])
    res[:malleParams][3]=append!(res[:malleParams][3][1:m1],[2,1])
  end
  return res
end)

chevieset(:I, :Representation, (m,i)->
  chevieget(:I, :HeckeRepresentation)(m,[[1,-1],[1,-1]],[1,1],i))

chevieset(:I, :HeckeRepresentation, function (m, param, rootparam, i)
  if i==1 return [[param[1][1];;], [param[2][1];;]] end
  if iseven(m) i-=2 end
  if i==0 
    [[param[1][1];;],[param[2][2];;]]
  elseif i==1 
    [[param[1][2];;],[param[2][1];;]]
  elseif i==2 
    [[param[1][2];;],[param[2][2];;]]
  else
    u=-param[1][1]//param[1][2]
    v=-param[2][1]//param[2][2]
    if isodd(m) squv=u
    elseif rootparam[1]!==nothing && rootparam[2]!==nothing
      squv=rootparam[1]*rootparam[2]
    else squv=root(u*v)
    end
    [-[-u^0 u^0;0u u]*param[1][2],
     -[v 0v;u+v+squv*(E(m,i-2)+E(m,2-i)) -v^0]*param[2][2]]
  end
end)

chevieset(:I, :Frobenius,(m, sqrtu, j)->
  [0 1//sqrtu//(E(2m,j)+E(2m,-j));sqrtu*(E(2m,j)+E(2m,-j)) 0]*sqrtu^0)

chevieset(:I, :PoincarePolynomial, function (m, param)
  u=-param[1][1]//param[1][2]
  v=-param[2][1]//param[2][2]
  if iseven(m) sum(i->(u*v)^(i-1),1:m//2)*(u+1)*(v+1)
  else sum(i->u^(i-1),1:m)*(u+1)
  end
end)

chevieset(:I,:SchurElement, function (m,phi,para,rootpara)
  if isodd(m)
    ci=chevieget(:I,:CharInfo)(m)
    ci=ci[:malleParams][findfirst(==(phi),ci[:charparams])]
    return chevieget(:imp,:SchurElement)(m,1,2,ci,[E.(m,0:m-1),para[2]],[])//m
  end
  u=-para[1][1]//para[1][2]
  v=-para[2][1]//para[2][2]
  m2=div(m,2)
  if phi[1]==1
    if phi[2]==m2
      e=sum(i->(u//v)^i,0:m2-1)*(u+1)*(v+1)//v
      if phi[3]==1 return e
      else return (v//u)^m2*e
      end
    else
      e=sum(i->(u*v)^i,0:m2-1)*(u+1)*(v+1)
      if phi[2]==0 return e
      else return (u*v)^-m2*e
      end
    end
  else
    e=E(m,phi[2])+E(m,-phi[2])
    if all(i->rootpara[i]!==nothing,[1,2]) ruv=prod(rootpara)
    else ruv=root(u*v)
    end
    -m*(u*v+1-ruv*e)*(u+v+e*ruv)//(u*v*(e^2-4))
  end
end)

chevieset(:I,:FakeDegree,(m,phi,q)->phi[1]==1 ? q^phi[2] : q^phi[2]+q^(m-phi[2]))

chevieset(:I,:CharTable,function(m)
  res=chevieget(:I,:HeckeCharTable)(m,[[1,-1],[1,-1]],[1,1])
  res[:identifier]=string("W(I2(",m,"))")
  res
end)

chevieset(:I, :FactorizedSchurElement, function(bond,phi,arg...)
  ci=chevieget(:I, :CharInfo)(bond)
  ci=ci[:malleParams][findfirst(==(phi),ci[:charparams])]
  chevieget(:imp, :FactorizedSchurElement)(bond,bond,2,ci,arg[1],1)
end)

chevieset(:I,:Invariants,function(e,typ=iseven(e) ? 1 : -E(e,(e+1)//2)-E(e,(e+3)//2))
  m=Diagonal([1+E(e,-1),-typ])*chevieget(:imp,:simpleroots)(e,e,2)
  map(f->(arg...)->f((transpose(collect(arg))*m)...),chevieget(:imp,:Invariants)(e,e,2))
end)

chevieset(:I, :DecompositionMatrix, function (n, p)
  T=chevieget(:I, :CharTable)(n)
  T[:name]=T[:identifier]
  m=DecompositionMatrix(mod(T, p))
  map(c->[c[1], m[c[1]][c[2]]], blocks(m))
end)

# parameters only for cuspidal symbols
chevieset(:I, :SymbolToParameter, function (S)
  if S[1]!=[0,1] || !any(isempty,S) return false end
  if isodd(length(S)) S=reverse(S) end
  a=findfirst(isempty,S)
  if isodd(length(S)) 
    [a,findfirst(==([0,1]),S)-a]
  else (a-1).+[-findfirst(==([0,1]),S[2:end]),0]
  end
end)

# The symbols returned are rotations of those given by Gunter.
# They are reduced in the sense of SymbolsDefect(e,2,0,0)
chevieset(:I,:ParameterToSymbol,function(e,p)
  if p==[0]
    S=map(x->[0],1:e)
    S[e]=[2]
  elseif p==[1]
    S=map(x->[0,1],1:e)
    S[e]=[1,2]
  elseif length(p)==3
    S=Vector{Any}(map(x->[0],1:e//2-1))
    append!(S,[[1],2,(p[3]+1)//2])
  elseif iseven(e)
    S=map(x->[0],1:e)
    if p[1]==0 S[[e,e-p[2]]]=[[1],[1]]
    else
      S[[0,mod(p[2]-p[1],e)].+1]=[[0,1],[0,1]]
      S[[mod(-p[1],e),p[2]].+1]=[Int[],Int[]]
    end
  else
    S=map(i->[0],1:e)
    if p[1]!=0
      S[[1,1+mod(-sum(p),e)]]=[[0,1],[0,1]]
      S[map(x->mod(x,e),-p).+1]=[Int[],Int[]]
    else
      S[[e-mod(p[2]-p[1],e),e]]=[[1],[1]]
    end
  end
  S
end)

chevieset(:I, :UnipotentCharacters, function(e)
  f=div(e,2)
  uc=Dict{Symbol, Any}()
  uc[:harishChandra] = [Dict{Symbol, Any}(:relativeType => 
    TypeIrred(;series=:I,indices=[1,2],rank=2,bond=e),
    :parameterExponents=>[1,1],:levi=>Int[],:eigenvalue=>1,:cuspidalName=>"")]
  uc[:harishChandra][1][:charNumbers]=
    isodd(e) ? (1:f+2) : vcat([1,3,4,2],4 .+(1:f-1))
# For I2(e) there are 3 families: Id, St and a big one.
# in the big one the cuspidal chars are S(k,l) where 0<k<l<e-k
  cusp=vcat(map(k->map(l->[k,l],k+1:e-k-1),1:f-1)...)
  f=f+1-mod(e,2) # number of principal series chars in big family
  append!(uc[:harishChandra],map(x->Dict{Symbol, Any}(:relativeType=>
    TypeIrred(;series=:A,indices=[],rank=0), 
    :parameterExponents=>[],:levi=>[1,2],:eigenvalue=>E(e,-prod(cusp[x])),
    :cuspidalName=>string("I_2(",e,")[",join(cusp[x],","),"]"),
    :charNumbers=>[x+f+2]),1:length(cusp)))
  uc[:families] = [Family(CHEVIE[:families][:Dihedral](e),
   (1:length(cusp)+f).+2),Family("C1", [1]), Family("C1", [2])]
  uc[:parameters]=vcat([[0],[1]],uc[:families][1][:parameters])
  uc[:charSymbols]=map(p->chevieget(:I,:ParameterToSymbol)(e,p),uc[:parameters])
  # for S(k,l) the b is min(k+l,e-k-l)
  uc[:a]=vcat([0,e],map(x->1,uc[:families][1][:parameters]))
  uc[:A]=vcat([0,e],map(x->e-1,uc[:families][1][:parameters]))
  if e==5 uc[:curtis]=[2, 1, 3, 4, 6, 5] end
  uc
end)

chevieset(:I,:Ennola,function(e)
  if isodd(e) return SPerm() end
  uc=chevieget(:I,:UnipotentCharacters)(e)
  l=uc[:charSymbols]
  SPerm(map(1:length(l))do i
    s=EnnolaSymbol(l[i])
    if s[end] isa AbstractVector u=circshift.(Ref(s),length(s):-1:1)
    else u=map(x->vcat(x,s[[end-1,end]]), 
               circshift.(Ref(s[1:end-2]),length(s)-2:-1:1))
    end
    for a in u
      p=findfirst(==(a),l)
      if !isnothing(p) return p*(-1)^uc[:A][i] end
    end
  end)
end)
