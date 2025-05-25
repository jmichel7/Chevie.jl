# Hand-translated part of chevie/tbl/weyld.g
# (C) Jean Michel, Frank Luebeck, Goetz Pfeiffer 1994-2001
# Data for type D

chevieset(:D,:CartanMat,n->cartan(:D,n))

chevieset(:D, :simpleroots, function (l)
  r=fill(0,l,l)
  for i in 1:l-1
    r[i,i:i+1]=[1,-1]
  end
  r[l,l-1:l]=[1,1]
  reverse(r,dims=1)
end)

chevieset(:D, :ReflectionDegrees,n->vcat(2:2:2n-2,n))

chevieset(:D, :WeightInfo,function(n)
  M=Matrix(1I,n,n)
  if isodd(n)
    for i in 3:n-1 M[i,n]=-mod(i-2,2) end
    M[1,2]=-1
    for i in 3:n M[i,2]=-2*mod(i-2,2) end
    Dict{Symbol, Any}(:minusculeWeights=>[1,2,n],:decompositions=>[[3],[1],[2]],
                      :moduli=>[4],:chosenAdaptedBasis=>M)
  else
    for i in 4:n-1 M[i,n]=-mod(i-3,2) end
    M[1,2]=M[1,n]=-1
    Dict{Symbol, Any}(:minusculeWeights=>[1,2,n],
   :decompositions=>[[1, 1],[1,0], [0,1]],:moduli=>[2,2],:chosenAdaptedBasis=>M)
  end
end)

chevieset(:D, :ParabolicRepresentatives,(l, s)->
  chevieget(:imp, :ParabolicRepresentatives)(2,2,l,s))

chevieset(:D, :WordsClassRepresentatives,function(n,param=partition_tuples(n,2))
  res=Vector{Int}[]
  for pi in param
    if ((pi[2] isa AbstractVector) && isodd(length(pi[2]))) || pi[2]=='-'
      continue
    end
    w=Int[]
    i=1
    if pi[2]!='+'
      for l in reverse(pi[2])
        if i==1 append!(w,2:l)
        else
          append!(w,i:-1:3)
          append!(w,1:i+l-1)
        end
        i+=l
      end
    end
    for l in pi[1]
      r=mod(l,2)
      append!(w,i.+vcat(1:2:l-1-r, 2:2:l+r-2))
      i+=l
    end
    if !isempty(w) && w[1]==2 w[1]=1 end # cosmetics for lexicographics.
    # classes are labelled with '+', if they have representatives
    # in parabolic subgroup of type A_{l-1}, given by {1,3,4,..}
    if (isempty(pi[2]) || pi[2]=='+') && all(iseven,pi[1])
      push!(res,w)
      w=copy(w)
      w[1]=2
    end
    push!(res, w)
  end
  res
end)

using Primes: primes
chevieset(:D, :ClassInfo, function (n)
  res=chevieget(:imp,:ClassInfo)(2,2,n)
  l=maximum(res[:orders])
  pmaps=Vector{Any}(fill(nothing,l))
  pp=res[:classparams]
  splits(S)=all(iseven,S[1]) && isempty(S[2])
  for pw in primes(l)
    pmaps[pw]=map(pp)do x
      px=chevieget(:imp,:pow)(x[1:2],pw)
      if splits(px) px=vcat(px,x[end]) end
      findfirst(y->y==px,pp)
    end
  end
  res[:powermaps]=pmaps
  res[:classparams]=map(x->x[end] isa Number ? [x[1],x[end]==0 ? '+' : '-'] : x,
                        res[:classparams])
  res[:classtext]=chevieget(:D,:WordsClassRepresentatives)(n,res[:classparams])
  res
end)

chevieset(:D, :NrConjugacyClasses,n->
  if isodd(n) div(npartition_tuples(n,2),2)
  else div(npartition_tuples(n,2)+3*npartitions(div(n,2)),2)
end)

chevieset(:D, :CharInfo, n->chevieget(:imp, :CharInfo)(2, 2, n))

#  If l is even then some of the classes  and restrictions split into two
#  classes or  characters, respectively. Their  values  are given by  the
#  character  values    for W(B_l) and   those  for  the  symmetric group
#  S_(l/2). This is described in [Pfeiffer, G., Character Tables of  Weyl
#  Groups in GAP]. 
chevieset(:D,:HeckeCharTable,function(n,para,root)
   ci=chevieget(:D,:ClassInfo)(n)
   cl=chevieget(:D,:CharInfo)(n)
   chp=cl[:charparams]
   q=-para[1][1]//para[1][2]
   tbl=Dict{Symbol,Any}(:name=>"H(D_$n)")
   tbl[:identifier]=tbl[:name]
   tbl[:parameter]=fill(q,n)
function chard(n,q)
  if n%2==0
    n1=div(n,2)-1
    AHk=chevieget(:A,:HeckeCharTable)(n1,fill([q^2,-1],n1),[])[:irreducibles]
    pA=partitions(n1+1)
    Airr(x,y)=AHk[findfirst(==(x),pA),findfirst(==(y),pA)]
  end
  BHk=isone(q) ? chevieget(:B,:CharTable)(n) :
      chevieget(:B,:HeckeCharTable)(n,vcat([[1,-1]],fill([q,-1],n)),[])
  pB=chevieget(:B,:CharInfo)(n)[:charparams]
  Birr(x,y)=BHk[:irreducibles][findfirst(==(x),pB),findfirst(==(y),pB)]
  function value(lambda,mu)
    if length(lambda)>2
      delta=[lambda[1], lambda[1]]
      if !(mu[2] isa Vector)
        vb=Birr(delta,[mu[1],Int[]])//2
        va=(q+1)^length(mu[1])//2*Airr(lambda[1],div.(mu[1],2))
        if "+-"[lambda[3]+1]==mu[2] val=vb+va
        else val=vb-va
        end
      else val=Birr(delta,mu)//2
      end
    else
      if !(mu[2] isa Vector) val=Birr(lambda, [mu[1], Int[]])
      else val=Birr(lambda, mu)
      end
    end
    return val
  end
  [[value(lambda,mu) for mu in ci[:classparams]] for lambda in chp]
end
  tbl[:irreducibles]=chard(n,q)
  tbl[:size]=prod(chevieget(:D,:ReflectionDegrees)(n))
#  tbl[:irredinfo]=List(chevieget(:D,:CharInfo)(n).charparams,p->
#     rec(charparam:=p,charname:=string_partition_tuple(p)));
  merge!(tbl,ci)
  merge!(tbl,cl)
  AdjustHeckeCharTable(tbl,para)
end)

chevieset(:D,:CharTable,n->chevieget(:D,:HeckeCharTable)(n,fill([1,-1],n),[]))

chevieset(:D,:CycPolPoincarePolynomial,n->CycPol(Pol()^n-1)*
          prod(i->CycPol(Pol()^2i-1),1:n-1)//CycPol(Pol()-1)^n)

chevieset(:D,:SchurElement, function (n, phi, para, sqrtparam)
  (chevieget(:D, :CycPolPoincarePolynomial)(n)//
   chevieget(:D, :CycPolGenericDegree)(phi))(-para[1][1]//para[1][2])
end)

chevieset(:D,:FactorizedSchurElement, function (n,p,arg...)
  if p[2] in "+-" p = [p[1], p[1]] end
  chevieget(:imp, :FactorizedSchurElement)(2, 2, n, p, arg[1], [])
end)

chevieset(:D,:HeckeRepresentation, function (n,para,rt,i,arg...)
  p=chevieget(:D, :CharInfo)(n)[:charparams][i]
  if p[end]==0 i+=1
  elseif p[end]==1 i-=1
  end
  chevieget(:imp, :HeckeRepresentation)(2, 2, n, para, [], i)
end)

chevieset(:D,:Representation, function (n, i)
  p=chevieget(:D, :CharInfo)(n)[:charparams][i]
  if p[end]==0 i+=1
  elseif p[end]==1 i-=1
  end
  return chevieget(:imp, :Representation)(2, 2, n, i)
end)

chevieset(:D,:PoincarePolynomial, function (n, para)
  q=-para[1][1]//para[1][2]
  sum(k->q^k,0:n-1)*prod(i->(q^i+1)*sum(k->q^k,0:i-1),1:n-1)
end)

chevieset(:D,:Invariants, function(n)
  m=chevieget(:imp,:simpleroots)(2,2,n)
  map(f->function(arg...) return f((transpose(collect(arg))*m)...)
    end,chevieget(:imp,:Invariants)(2,2,n))
end)

chevieset(:D,:CycPolGenericDegree,c->gendeg(Symbol_partition_tuple(c,0)))

chevieset(:D,:FakeDegree,(n,c,q)->fakedegree(Symbol_partition_tuple(c, 0))(q))

chevieset(:D,:UnipotentCharacters,function(rank)
  uc=Dict{Symbol, Any}(:harishChandra=>[],:charSymbols=>[])
  for d in 0:4:4*isqrt(div(rank,4))
    r=div(d^2,4)
    s=Dict{Symbol, Any}(:relativeType=> 
      TypeIrred(;series=:B,indices=1+r:rank,rank=rank-r),:levi=>1:r,
      :eigenvalue=>(-1)^div(d+1,4), 
      :parameterExponents=>vcat(d,fill(1,max(0, rank-1-r))))
    s[:cuspidalName]="D"*stringind(rio(TeX=true),r)
    if d==0
      s[:relativeType].series=:D
      s[:cuspidalName]=""
      s[:parameterExponents][1]=1
    end
    push!(uc[:harishChandra], s)
    if d==0
      # order differs from Symbols.Symbolsshape(n,[d,0]) for d=0
      symbols=map(x->Symbol_partition_tuple(x,d),
        chevieget(:imp,:CharInfo)(2,2,rank)[:charparams])
    else symbols=Symbols.Symbolsshape(rank,[d,0])
    end
    s[:charNumbers]=(1:length(symbols)).+length(uc[:charSymbols])
    FixRelativeType(s)
    append!(uc[:charSymbols],symbols)
  end
  uc[:a]=valuation_gendeg.(uc[:charSymbols])
  uc[:A]=degree_gendeg.(uc[:charSymbols])
  uc[:families]=FamiliesClassical(uc[:charSymbols])
  uc
end)

chevieset(:D,:Ennola,function(n)
  if isodd(n) return SPerm() end
  uc=chevieget(:D,:UnipotentCharacters)(n)
  l=uc[:charSymbols]
  SPerm(map(function(i)
    if !(l[i][2] isa AbstractVector) return i*(-1)^uc[:A][i] end
    s=ennola(l[i])
    p=findfirst(==(s),l)
    if isnothing(p) p=findfirst(==(reverse(s)),l) end
    p*(-1)^uc[:A][i]
  end, 1:length(l)))
end)

# References for unipotent classes:
# [Lu] G.Lusztig, Character sheaves on disconnected groups, II 
#   Representation Theory 8 (2004) 72--124
#
# [GM]  M.Geck and G.Malle, On the existence of a unipotent support for the
# irreducible  characters of  a finite  group of  Lie type,  Trans. AMS 352
# (1999) 429--456
# 
# [S]  N.Spaltenstein,  Classes  unipotentes  et  sous-groupes  de
# Borel, Springer LNM 946 (1982)
# 
chevieset(:D,:UnipotentClasses,function(n,char)
  function addSpringer(s, i, cc)
    ss=first(x for x in uc[:springerSeries] 
                     if x[:defect]==defectsymbol(s[:symbol]))
    if s.sp in [[Int[],[1]],[Int[],Int[]]] p=1
    elseif s.sp==[[1], Int[]] p=2
    else p=findfirst(==([s.sp]),charinfo(ss[:relgroup]).charparams)
    end
    ss[:locsys][p]=[i,findfirst(==(map(x->x ? [1,1] : [2],s.Au)),
                                           charinfo(cc[:Au]).charparams)]
  end
  function partition2DR(part)
    p=sort(reduce(vcat,map(x->1-x:2:x-1, part)))
    p=p[1+div(length(p),2):end]
    vcat([p[1]+p[2]], map(i->p[i+1]-p[i], 1:length(p)-1))
  end
  if char==2
    ss=XSP(4, 0, n, 1)
    symbol2partition=function(S) # see [GM] 2.17
      c=sort(reduce(vcat,S))
      i = 1
      part = Int[]
      ex = Int[]
      while i <= length(c)
        l=2(c[i]-2(i-1))
        if i==length(c) || c[i+1]-c[i]>1
          push!(part,l+2)
          i+=1
        elseif c[i+1]-c[i]>0
          append!(part,[l+1,l+1])
          i+=2
        else
          append!(part,[l,l])
          i+=2
          push!(ex,l)
        end
      end
      [reverse(filter(!iszero,sort(part))), ex]
    end
  else # see [GM] 2.10
    ss=XSP(2, 0, n, 1)
    symbol2partition=function(S)
      c=sort(reduce(vcat,S))
      i = 1
      part = Int[]
      while i <= length(c)
        l = 2 * (c[i]-(i-1))
        if i == length(c) || c[i+1]-c[i]>0
          push!(part, l+1)
          i+=1
        else
          part = append!(part, [l, l])
          i+=2
        end
      end
      reverse(filter(!iszero,sort(part)))
    end
  end
  l=union(map(c->map(x->[defectsymbol(x.symbol),
                         sum(sum,fullsymbol(x.sp))],c),ss)...)
  sort!(l, by=x->[abs(x[1]), -sign(x[1])])
  uc = Dict{Symbol, Any}(:classes => [], :springerSeries => map(function(d)
      res = Dict{Symbol, Any}(:defect=>d[1], :levi=>1:n-d[2])
      if mod(n-d[2],4)==0 || char==2 res[:Z]=mod(n,2)==0  ? [1, 1] : [1]
      else   res[:Z]=mod(n,2)==0  ? [-1,1] : [-1]
      end
      res[:relgroup]=coxgroup(d[1]==0 ? :D : :B, d[2])
      res[:locsys]=[[0,0] for i in 1:nconjugacy_classes(res[:relgroup])]
      res
  end, l))
  for cl in ss
    cc = Dict{Symbol, Any}(:parameter=>symbol2partition(cl[1].symbol))
    if char==2
      cc[:dimBu] = cl[1].dimBu
      cc[:name] = join(map(reverse(tally(cc[:parameter][1])))do x
        res=joindigits(fill(x[1], max(0, x[2])), "[]")
        if x[1] in cc[:parameter][2] return string("(", res, ")") end
        res
      end)
    else
      cc[:dynkin]=partition2DR(cc[:parameter])
      cc[:name] = joindigits(cc[:parameter])
    end
    cc[:Au] = isempty(cl[1].Au) ? coxgroup() : 
       prod(coxgroup(:A,1) for i in eachindex(cl[1].Au))
    if char != 2
      cc[:red] = coxgroup()
      j = cc[:parameter]
      for j = tally(j)
        if mod(j[1], 2) == 0
          cc[:red]*=coxgroup(:C, div(j[2],2))
        elseif mod(j[2], 2) != 0
          if j[2]>1 cc[:red]*=coxgroup(:B, div(j[2]-1,2)) end
        elseif j[2]>2 cc[:red]*=coxgroup(:D, div(j[2],2))
        else cc[:red]*=torus(1)
        end
      end
    end
    if !(cl[1].sp[2] isa Vector) cl[1].sp[3]=1-mod(div(n,2),2) end
    push!(uc[:classes], cc)
    for s in cl addSpringer(s, length(uc[:classes]), cc) end
    if !(cl[1].sp[2] isa Vector)
      cl[1].sp[3]=1-cl[1].sp[3]
      cc[:name]*="+"
      cc=deepcopy(cc)
      cc[:name]=replace(cc[:name],r".$"=>"-")
      if haskey(cc, :dynkin) cc[:dynkin][[1, 2]] = cc[:dynkin][[2, 1]] end
      push!(uc[:classes], cc)
      for s in cl addSpringer(s, length(uc[:classes]), cc) end
    end
  end
  if char==2 # see [Spaltenstein] 2.10 page 24
    uc[:orderClasses]=hasse(Poset(uc[:classes])do x,y
      xp=x[:parameter] 
      yp=y[:parameter] 
      m=max(xp[1][1], yp[1][1])
      f=x->map(i->sum(filter(<(i),x))+i*count(>=(i),x),1:m)
      fx=f(xp[1])
      fy=f(yp[1])
      for i in 1:m
        if fx[i]<fy[i] return false
        elseif fx[i]==fy[i] && i in yp[2]
          if i in setdiff(xp[1], xp[2]) return false end
          if i<m && mod(fx[i+1]-fy[i+1],2)==1 return false end
        end
      end
      xp!=yp || x==y
    end)
  else
    uc[:orderClasses]=hasse(Poset(uc[:classes])do x,y
      xp=x[:parameter] 
      yp=y[:parameter] 
      dominates(yp,xp) && (xp!=yp || x==y)
    end)
  end
  if char!=2
    d = 0
    while 4*d^2-d <= n
     i = 4*d^2-d
     if mod(n-d,2)==0
       l=vcat(1:i,i+2:2:n)
       s=Dict(:relgroup=>coxgroup(:B, div(n-i,2)),:levi=>l)
       if mod(n, 2) == 0 s[:Z]=[-1,-1] else s[:Z] = [E(4)] end
       s[:locsys]=[[0,0] for i in 1:nconjugacy_classes(s[:relgroup])]
       push!(uc[:springerSeries], s)
       if d==0 l=vcat([1], 4:2:n) end
       s=Dict(:relgroup => coxgroup(:B, div(n-i,2)),:levi=>l)
       if mod(n,2)==0 s[:Z] = [1,-1]
       else s[:Z] = [-(E(4))]
       end
       s[:locsys]=[[0,0] for i in 1:nconjugacy_classes(s[:relgroup])]
       push!(uc[:springerSeries], s)
       i=4d^2+d
       if d != 0 && i <= n
         l = vcat(1:i, i+2:2:n)
         s = Dict(:relgroup=>coxgroup(:B, div(n-i,2)),:levi=>l)
         if mod(n,2)==0 s[:Z]=[-1,-1] else s[:Z]=[E(4)] end
         s[:locsys]=[[0,0] for i in 1:nconjugacy_classes(s[:relgroup])]
         push!(uc[:springerSeries], s)
         s=Dict(:relgroup=>coxgroup(:B,div(n-i,2)),:levi=>l)
         if mod(n, 2) == 0 s[:Z] = [1, 1] else s[:Z] = [-(E(4))] end
         s[:locsys]=[[0,0] for i in 1:nconjugacy_classes(s[:relgroup])]
         push!(uc[:springerSeries], s)
       end
     end
     d+=1
  end
  function LuSpin(p) # see [Lusztig] 14.2
    sort!(p)
    a=Int[]
    b=Int[]
    d=[0,1,0,-1][map(x->1 + mod(x, 4), p)]
    i=1
    while i<=length(p)
      l=p[i]
      t=sum(d[1:i-1])
      if 1==mod(l,4)
        push!(a, div(l-1,4)-t)
        i+=1
      elseif 3==mod(l,4)
        push!(b, div(l-3,4)+t)
        i+=1
      else
        j=i
        while i<=length(p) && p[i]==l i+=1 end
        j=fill(0, max(0,div(i-j,2)))
        append!(a,j.+div(l+mod(l,4),4).-t)
        append!(b,j.+div(l-mod(l,4),4).+t)
      end
    end
    a=reverse(filter(!iszero,a))
    b=reverse(filter(!iszero,b))
    sum(d)>=1 ? [a, b] : [b, a]
  end
  function addSpringer1(f, i, s, k)
    ss=first(x for x in uc[:springerSeries] if f(x))
    if s in [[Int[], [1]], [Int[], Int[]]] p = 1
    elseif s == [[1], Int[]] p = 2
    else p = findfirst(==([s]),charinfo(ss[:relgroup]).charparams)
    end
    ss[:locsys][p] = [i, k]
  end
  function trspringer(i, new)
    for ss in uc[:springerSeries]
      for c in ss[:locsys] if c[1]==i c[2]=new[c[2]] end end
    end
  end
  l=filter(i->all(c->mod(c[1],2)==0 || c[2]==1,
                  tally(uc[:classes][i][:parameter])),eachindex(uc[:classes]))
  for i in l
     cl=uc[:classes][i]
     s=LuSpin(cl[:parameter])
     if length(cl[:Au]) == 1
       cl[:Au] = coxgroup(:A, 1)
       trspringer(i, [2])
       k = [1, 1]
     elseif length(cl[:Au]) == 2
       cl[:Au] = coxgroup(:A, 1)*coxgroup(:A,1)
       trspringer(i, [2, 4])
       k = [1, 3]
     elseif length(cl[:Au]) == 8
       cl[:Au] = coxgroup(:A,1)*coxgroup(:B,2)
       trspringer(i, [1, 6, 8, 5, 10, 3, 4, 9])
       k = [2, 7]
     else
       error("Au non-commutative of order ",length(cl[:Au])*2,"  ! implemented")
     end
     if !('-' in cl[:name])
      addSpringer1(ss->ss[:Z] in [[(-1)^(div(n,2)+1),-1], [E(4)]] && 
                     rank(ss[:relgroup]) == sum(sum,s) , i, s, k[1])
     end
     if !('+' in cl[:name])
      addSpringer1(ss->ss[:Z] in [[(-1)^div(n,2),-1], [-(E(4))]] && 
                     rank(ss[:relgroup]) == sum(sum,s), i, s, k[2])
     end
  end
end
  for ss in uc[:springerSeries] if !all(isone,ss[:Z]) ss[:hc]=0 end end
  uc
end)

function decode(l)
  res=Int[]
  for (n,m) in enumerate(l)
   if !isnothing(m) append!(res,fill(n+1,m)) end
  end
  tally(res)
end

# Used by 'ClassParamD' for distinguishing classes with '+' or '-' in label;
# precomputed for D_n with n=4,6,8
chevieset(:D, :gensMODA,Dict(4=>[[perm"(1,2)(7,8)", perm"(3,4)(5,6)", perm"(2,3)(6,7)", perm"(3,5)(4,6)"],[[2 => 4], [4 => 2]], [[2 => 2], [2 => 1, 4 => 1]]],
6=>[[perm"(1,2)(8,11)(12,14)(15,17)(16,18)(19,21)(22,25)(31,32)", 
     perm"(3,4)(5,6)(7,9)(10,13)(20,23)(24,26)(27,28)(29,30)", 
     perm"(2,3)(6,8)(9,12)(13,16)(17,20)(21,24)(25,27)(30,31)", 
     perm"(3,5)(4,6)(12,15)(14,17)(16,19)(18,21)(27,29)(28,30)", 
     perm"(5,7)(6,9)(8,12)(11,14)(19,22)(21,25)(24,27)(26,28)", 
     perm"(7,10)(9,13)(12,16)(14,18)(15,19)(17,21)(20,24)(23,26)"], 
    [[2 => 16], [2 => 4, 4 => 6], [2 => 1, 6 => 5]],
    [[2 => 12], [2 => 2, 4 => 6], [3 => 2, 6 => 4]]],
8=>[[perm"(1,2)(8,11)(12,15)(16,20)(17,21)(22,26)(23,27)(28,33)(29,34)(30,35)(36,41)(37,42)(43,50)(44,51)(52,59)(60,68)(61,69)(70,77)(78,85)(79,86)(87,92)(88,93)(94,99)(95,100)(96,101)(102,106)(103,107)(108,112)(109,113)(114,117)(118,121)(127,128)", 
     perm"(3,4)(5,6)(7,9)(10,13)(14,18)(19,24)(25,31)(32,38)(39,45)(40,46)(47,53)(48,54)(49,55)(56,62)(57,63)(58,64)(65,71)(66,72)(67,73)(74,80)(75,81)(76,82)(83,89)(84,90)(91,97)(98,104)(105,110)(111,115)(116,119)(120,122)(123,124)(125,126)", 
     perm"(2,3)(6,8)(9,12)(13,17)(18,23)(20,25)(24,30)(26,32)(33,39)(34,40)(41,48)(42,49)(50,57)(51,58)(53,61)(59,67)(62,70)(68,76)(71,78)(72,79)(80,87)(81,88)(89,95)(90,96)(97,103)(99,105)(104,109)(106,111)(112,116)(117,120)(121,123)(126,127)", 
     perm"(3,5)(4,6)(12,16)(15,20)(17,22)(21,26)(23,29)(27,34)(30,37)(35,42)(39,47)(45,53)(48,56)(54,62)(57,65)(58,66)(63,71)(64,72)(67,75)(73,81)(76,84)(82,90)(87,94)(92,99)(95,102)(100,106)(103,108)(107,112)(109,114)(113,117)(123,125)(124,126)", 
     perm"(5,7)(6,9)(8,12)(11,15)(22,28)(26,33)(29,36)(32,39)(34,41)(37,44)(38,45)(40,48)(42,51)(46,54)(49,58)(55,64)(65,74)(71,80)(75,83)(78,87)(81,89)(84,91)(85,92)(88,95)(90,97)(93,100)(96,103)(101,107)(114,118)(117,121)(120,123)(122,124)", 
     perm"(7,10)(9,13)(12,17)(15,21)(16,22)(20,26)(25,32)(31,38)(36,43)(41,50)(44,52)(48,57)(51,59)(54,63)(56,65)(58,67)(62,71)(64,73)(66,75)(70,78)(72,81)(77,85)(79,88)(86,93)(91,98)(97,104)(103,109)(107,113)(108,114)(112,117)(116,120)(119,122)", 
     perm"(10,14)(13,18)(17,23)(21,27)(22,29)(26,34)(28,36)(32,40)(33,41)(38,46)(39,48)(45,54)(47,56)(52,60)(53,62)(59,68)(61,70)(67,76)(69,77)(73,82)(75,84)(81,90)(83,91)(88,96)(89,97)(93,101)(95,103)(100,107)(102,108)(106,112)(111,116)(115,119)", 
     perm"(14,19)(18,24)(23,30)(27,35)(29,37)(34,42)(36,44)(40,49)(41,51)(43,52)(46,55)(48,58)(50,59)(54,64)(56,66)(57,67)(62,72)(63,73)(65,75)(70,79)(71,81)(74,83)(77,86)(78,88)(80,89)(85,93)(87,95)(92,100)(94,102)(99,106)(105,111)(110,115)"], 
  [[2=>64],[2=>16,4=>24],[4=>32],[2=>4,6=>20],[8=>16]],
  [[2=>56],[2=>12,4=>24],[2=>6,4=>28],[2=>2,3=>4,6=>18],[2=>1,4=>3,8=>14]]]))

chevieset(:D, :ClassParameter, function (n, w)
  res=cycletype(prod(i->i==1 ? SPerm(1,-2) : SPerm(i-1,i),w;init=SPerm()),n)
  if isempty(res[2]) && all(iseven, res[1])
  # for classes with '+' or '-' we use the cycle type for the
  # permutation representation on the cosets of the parabolic
  # subgroup 2:n (gens(Dₙ) in this representation are 
  # stored for n<=8 in the Dict 'chevieget(:D,:gensMODA)')
    u=get!(chevieget(:D,:gensMODA),n)do
      Dn=coxgroup(:D,n)
      R=reflection_subgroup(Dn,2:n)
      D=vcat(reduced(R,Dn)...)
      Dgens=map(s->Perm(reduced.(Ref(R),D.*s),D),gens(Dn))
      ci=conjugacy_classes(Dn)
      H=map(('+','-'))do c
        l=word.(ci)[filter(i->ci[i].name[end]==c,eachindex(ci))]
        map(a->tally(cycletype(prod(Dgens[a]))),l)
      end
      [Dgens,H[1],H[2]]
    end
    tmp=tally(cycletype(prod(u[1][w])))
    if tmp in u[2] && !(tmp in u[3]) res=[res[1],'+']
    elseif !(tmp in u[2]) && tmp in u[3] res=[res[1],'-']
    end
  end
  sort!(res[1],rev=true)
  if res[2] isa Vector sort!(res[2],rev=true) end
  res
end)
