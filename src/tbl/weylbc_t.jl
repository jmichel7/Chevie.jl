# Hand-translated part of chevie/tbl/weylbc.g
# Jean Michel & Goetz Pfeiffer (C) 1994-2012
#
# Data for types B and C
#
#  The  functions which behave  differently depending on  the isogeny type
#  take an argument after the rank, the cartanType.
#    cartanType=1 type C
#    cartanType=2 type B
#    cartanType=root(2) and rank=2 Suzuki groups.
#  other  values  happen  for  root  systems  which  occur  inside complex
#  reflection groups.

chevieset(:B, :CartanMat, function(n,ct=2)
  a=chevieget(:A, :CartanMat)(n)//one(ct)
  if n>=2
    a[1][2]=-ct
    a[2][1]=-2//ct
  end
  improve_type(a)
end)

chevieset(:B, :ReflectionName, function(r,option,cartantype=2)
  if cartantype==2
    if haskey(option,:TeX)
      if haskey(option,:Au) "D_8"
      else string("B",GAPENV.TeXIndex(r))
      end
    elseif haskey(option,:Au) "D8"
    elseif haskey(option,:arg) string("\"B\",",r)
    else string("B",r)
    end
  elseif cartantype==1
    if haskey(option,:TeX) string("C",GAPENV.TeXIndex(r))
    elseif haskey(option,:arg) string("\"C\",", r)
    else string("C", r)
    end
  elseif cartantype==root(2)
    if haskey(option,:TeX) string("B^{\\hbox{sym}}",GAPENV.TeXIndex(r))
    elseif haskey(option,:arg) string("\"Bsym\",",r)
    else string("Bsym",r)
    end
  elseif haskey(option,:TeX)
    string("B^?_",GAPENV.TeXIndex(r),"(",xrepr(cartantype;option...),")")
  elseif haskey(option,:arg) string("\"B?\",", r, ",", cartantype)
  else string("B?",r,"(",xrepr(cartantype),")")
  end
end)

chevieset(:B,:GeneratingRoots, function(l,type_)
  rts=map(i->fill(0,l),1:l)
  for i in 1:l-1 rts[i][i:i+1]=[1,-1] end
  rts[l][l]=2//type_
  rts[l:-1:1]
end)

chevieset(:B,:ParabolicRepresentatives,(l,s)->
  chevieget(:imp, :ParabolicRepresentatives)(2, 1, l, s))

chevieset(:B,:ReflectionDegrees,n->2:2:2n)

chevieset(:B,:Size,(r,arg...)->2^r*factorial(r))

chevieset(:B,:NrConjugacyClasses,(r,arg...)->npartition_tuples(r,2))

chevieset(:B,:WeightInfo,function(n,type_)
  M=Array(Int.(I(n)))
  if type_ == 2
    Dict{Symbol, Any}(:minusculeWeights => [1], :minusculeCoweights => [n], :decompositions => [[1]], :moduli => [2], :chosenAdaptedBasis => M)
  else
    for i in 1:n-1
      M[i,n]=-mod(n+i-1,2)
    end
    Dict{Symbol, Any}(:minusculeWeights => [n], :minusculeCoweights => [1], :decompositions => [[1]], :moduli => [2], :chosenAdaptedBasis => M)
  end
end)

# returns a very good representative in the sense of [Geck-Michel]
chevieset(:B,:WordClass, function(pi)
  w=Int[]
  i=1
  for l in reverse(pi[2])
    append!(w,i:-1:2)
    append!(w,1:i+l-1)
    i+=l
  end
  for l in pi[1]
    r=mod(l, 2)
    append!(w,i.+vcat(1:2:l-1-r,2:2:l+r-2))
    i+=l
  end
  w
end)

chevieset(:B,:ClassInfo, function(n)
  res=chevieget(:imp,:ClassInfo)(2,1,n)
  res[:classtext]=map(chevieget(:B,:WordClass),res[:classparams])
  res[:classes]=div.(res[:centralizers][1],res[:centralizers])
  res
end)

chevieset(:B,:LowestPowerFakeDegree,function(p)
  pp=symbol_partition_tuple(p,1)
  m=length(pp[2])
  res=dot(pp[1],m:-1:0)
  if !isempty(pp[2]) res+=dot(pp[2],m-1:-1:0) end
  2res+sum(pp[2])-div(m*(m-1)*(4m+1),6)
end)

chevieset(:B,:CharInfo,n->chevieget(:imp, :CharInfo)(2,1,n))

chevieset(:B,:PoincarePolynomial,function(n,para)
  q1=-para[1][1]//para[1][2]
  q2=-para[2][1]//para[2][2]
  prod(i->(q2^i*q1+1)*sum(k->q2^k,0:i),0:n-1)
end)

chevieset(:B,:CharTable,n->chevieget(:imp,:CharTable)(2,1,n))

chevieset(:B,:HeckeCharTable,(n,para,root)->chevieget(:imp,:HeckeCharTable)(2,1,n,para,root))

chevieset(:B,:SchurElement,(arg...)->
  chevieget(:imp, :SchurElement)(2,1,arg[1],arg[2],arg[3],[]))

chevieset(:B,:FactorizedSchurElement,(arg...)->
  chevieget(:imp, :FactorizedSchurElement)(2, 1, arg[1], arg[2], arg[3], []))

chevieset(:B,:HeckeRepresentation,(arg...)->
  chevieget(:imp, :HeckeRepresentation)(2, 1, arg[1], arg[2], [], arg[4]))

chevieset(:B,:Representation,(n,i)->
  chevieget(:imp, :Representation)(2, 1, n, i))

chevieset(:B,:FakeDegree,(n,c,q)->fegsymbol(symbol_partition_tuple(c,1))(q))

chevieset(:B, :DecompositionMatrix, function (l, p)
  decS(i)= MatrixDecompositionMatrix(DecompositionMatrix(Specht(p, p), i))
  pp=map(partitions, 0:l)
  pt=partition_tuples(l, 2)
  if p==2
   [[1:length(pt), map(pt)do p
     p = LittlewoodRichardsonRule(p[1], p[2])
     map(x->x in p,pp[l+1]) end *decS(l)]]
  else
    dd=vcat([[[1]],[[1]]],decS.(2:l))
    map(0:l)do i
     [map(x->Position(pt, x), cartesian(pp[i+1], pp[l+1-i])),
      map(x->prod.(cartesian(x)), cartesian(dd[i+1],dd[l+1-i]))]
    end
  end
end)

chevieset(:B, :UnipotentCharacters,function(rank,typ=2)
  uc=Dict{Symbol, Any}(:harishChandra=>[],:charSymbols=>[])
  for d in (0:div(-1+GAPENV.RootInt(1+4rank,2),2)).*2 .+1
    r=div(d^2-1,4)
    s=Dict{Symbol, Any}()
    s[:relativeType]=TypeIrred(;series=:B,indices=1+r:rank,rank=rank-r)
    s[:levi]=1:r
    s[:eigenvalue]=(-1)^div(d+1,4)
    s[:parameterExponents]=vcat([d],fill(1,max(0,rank-1-r)))
    s[:cuspidalName]=string("B",GAPENV.TeXIndex(r))
    push!(uc[:harishChandra], s)
    symbols = BDSymbols(rank, d)
    s[:charNumbers]=(1:length(symbols)).+length(uc[:charSymbols])
    FixRelativeType(s)
    append!(uc[:charSymbols],symbols)
  end
  uc[:harishChandra][1][:cuspidalName] = ""
  uc[:a]=valuation_gendeg_symbol.(uc[:charSymbols])
  uc[:A]=degree_gendeg_symbol.(uc[:charSymbols])
  uc[:families]=FamiliesClassical(uc[:charSymbols])
  if typ==1 uc[:harishChandra][1][:relativeType][:cartanType]=1 end
  uc
end)

chevieset(:B, :Ennola, function(n)
  uc=chevieget(:B,:UnipotentCharacters)(n)
  l=uc[:charSymbols]
  SPerm(map(1:length(l))do i
    s=EnnolaSymbol(l[i])
    if length(s[1])<length(s[2]) s=s[[2,1]] end
    findfirst(==(s),l)*(-1)^uc[:A][i]
  end)
end)

chevieset(:B,:Invariants,function(n,type_)
  m=fill(1,n);m[1]=2//type_
  m=Diagonal(m)*toM(chevieget(:imp,:GeneratingRoots)(2,1,n))
  map(f->(arg...)->f(transpose(collect(arg))*m...),CHEVIE[:imp][:Invariants](2,1,n))
end)

# References:
# [Lu]   G.Lusztig,   Character   sheaves   on   disconnected   groups,  II
# Representation Theory 8 (2004) 72--124
#
# [GM]  M.Geck and G.Malle, On the existence of a unipotent support for the
# irreducible  characters of  a finite  group of  Lie type,  Trans. AMS 352
# (1999) 429--456
# 
# [S]   N.Spaltenstein,  Classes  unipotentes  et  sous-groupes  de  Borel,
# Springer LNM 946 (1982)
#
#  The function is called with type=1 for type C and type=2 for type B
#
chevieset(:B,:UnipotentClasses,function(r,char,ctype)
  function part2dynkin(part)
    p=sort(reduce(vcat,map(d->1-d:2:d-1, part)))
    p=p[div(3+length(p),2):end]
    res= ctype==1 ? [2*p[1]] : [p[1]]
    append!(res,p[2:end]-p[1:end-1])
  end
  function addSpringer1(s,cc)
    ss=first(x for x in uc[:springerSeries] 
                     if x[:defect]==defectsymbol(s.symbol))
    if s.sp==[Int[], Int[]] p=1
    elseif s.sp==[[1], Int[]] p=2
    elseif s.sp==[Int[], [1]] p=1
    else p=findfirst(==([s.sp]),charinfo(ss[:relgroup]).charparams)
    end
    ss[:locsys][p]=[length(uc[:classes]), 
    findfirst(==(map(x->x ? [1,1] : [2], s.Au)),charinfo(cc[:Au]).charparams)]
  end
  if ctype==root(2)  #treat 2B2 as B2; make sure char=2
    ctype=2
    char=2
  end
  ss=char==2 ? XSP(4,2,r) : ctype==1 ? XSP(2,1,r) : XSP(2,0,r)
  l=union(map(c->map(x->[defectsymbol(x.symbol), sum(sum,x.sp)], c), ss)...)
  sort!(l,by=x->[abs(x[1]),-sign(x[1])])
  uc=Dict{Symbol, Any}(:classes=>[])
  uc[:springerSeries]=map(l)do d
    res=Dict(:relgroup=>coxgroup(:C,d[2]),:defect=>d[1],:levi=>1:r-d[2])
    res[:locsys]=fill([0, 0],nconjugacy_classes(res[:relgroup]))
    if char==2 res[:Z]=[1]
    elseif ctype==1 res[:Z]=[(-1)^(r-d[2])]
    elseif conductor(root(2*(r-d[2])+1))==1 res[:Z]=[1]
    else res[:Z]=[-1]
    end
    res
  end
  function symbol2para(S)
    c=sort(reduce(vcat,S))
    i=1
    part=Int[]
    if char==2 # Invert [GM] 2.7
      ex=Int[]
      while i<=length(c)
        l=2*(c[i]-2*(i-1))
        if i==length(c) || c[i+1]>c[i]+1
          push!(part, l)
          i+=1
        elseif c[i]==c[i+1]
          append!(part, [l-2, l-2])
          push!(ex,l-2)
          i+=2
        elseif c[i]+1==c[i+1]
          append!(part, [l-1, l-1])
          i+=2
        end
      end
      [reverse(filter(y->y!=0,sort(part))), ex]
    else  # Invert [GM] 2.6 and 2.10
      d=mod(ctype,2)
      while i<=length(c)
        l=2*(c[i]-(i-1))-d
        if i==length(c) || c[i+1]>c[i]
          push!(part, l+1)
          i+=1
        else
          append!(part,[l,l])
          i+=2
        end
      end
      reverse(filter(!iszero,sort(part)))
    end
  end
  if char==2 ctype=1 end
  for cl in ss
    cc = Dict{Symbol, Any}(:parameter => symbol2para(cl[1].symbol))
    v=cl[1].Au
    cc[:Au]=isempty(v) ? coxgroup() : prod(map(x->coxgroup(:A,1),v))
    if char!=2
      cc[:dynkin] = part2dynkin(cc[:parameter])
      cc[:name] = joindigits(cc[:parameter])
    else
      ctype = 1
      cc[:dimBu] = cl[1].dimBu
      cc[:name] = join(map(function(x)
        res=joindigits(fill(0,max(0,x[2])).+x[1],"[]")
        if x[1] in cc[:parameter][2] return string("(", res, ")") end
        res
      end, reverse(tally(cc[:parameter][1]))), "")
    end
    cc[:red]=coxgroup()
    if char==2 j=cc[:parameter][1]
    else j=cc[:parameter]
    end
    for j in tally(j)
      if mod(j[1],2)==mod(ctype,2) cc[:red]*=coxgroup(:C, div(j[2],2))
      elseif mod(j[2],2)!=0
        if j[2]>1 cc[:red]*=coxgroup(:B, div(j[2]-1,2)) end
      elseif j[2]>2 cc[:red]*=coxgroup(:D, div(j[2],2))
      else cc[:red]*=torus(1)
      end
    end
    push!(uc[:classes], cc)
    for s in cl addSpringer1(s,cc) end
  end
  uc[:orderClasses] = hasse(toM(map(x->
    map(function(y)
      if char!=2 return dominates(y[:parameter], x[:parameter]) end
      # cf. [S] 2.10 page 24
      m=max(x[:parameter][1][1], y[:parameter][1][1])
      f=x-> map(i->sum(filter(<(i),x))+i*count(>=(i),x) ,1:m)
      fx=f(x[:parameter][1])
      fy=f(y[:parameter][1])
      for i in 1:m
        if fx[i]<fy[i] return false
        elseif fx[i]==fy[i] && i in y[:parameter][2]
          if i in setdiff(x[:parameter][1],x[:parameter][2]) return false end
          if i<m && mod(fx[i+1]-fy[i+1],2)==1 return false end
        end
      end
      return true
    end, uc[:classes]), uc[:classes])))
  if char!=2 && ctype==2
    function LuSpin(p) # cf [Lu] 14.2
      sort!(p)
      a=Int[]
      b=Int[]
      d=[0, 1, 0, -1][map(x->1+mod(x,4), p)]
      i=1
      while i<=length(p)
        l=p[i]
        t=sum(d[1:i-1])
        if 1==mod(l, 4)
          push!(a, div(l-1, 4)-t)
          i+=1
        elseif 3==mod(l, 4)
          push!(b, div(l-3, 4)+t)
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
    function addSpringer(f, i, s, k)
      ss=first(x for x in uc[:springerSeries] if f(x))
      if s in [[Int[],[1]],[Int[],Int[]]] p=1
      elseif s==[[1],Int[]] p=2
      else p = findfirst(==([s]),charinfo(ss[:relgroup]).charparams)
      end
      ss[:locsys][p] = [i, k]
    end
    function trspringer(i, old, new)
      for ss in uc[:springerSeries]
        for c in ss[:locsys]
          if c[1] == i
            p=findfirst(==(c[2]),old)
            if !isnothing(p) c[2]=new[p] end
          end
        end
      end
    end
    d = 0
    while 4d^2-3d<=r
      i=4d^2-3d
      if mod(r-d,2)==0
          l=vcat(1:i,i+2:2:r)
          ss=Dict{Symbol, Any}(:relgroup=>coxgroup(:B,div(r-i,2)),
                               :levi => l, :Z => [-1])
          ss[:locsys]=fill([0,0],nconjugacy_classes(ss[:relgroup]))
          push!(uc[:springerSeries],ss)
          i=4d^2+3d
          if i<=r && d!=0
            l=vcat(1:i,i+2:2:r)
            ss= Dict{Symbol, Any}(:relgroup=>coxgroup(:B,div(r-i,2)),
                                  :levi => l, :Z => [-1])
            ss[:locsys]=fill([0,0],nconjugacy_classes(ss[:relgroup]))
            push!(uc[:springerSeries], ss)
          end
      end
      d+=1
    end
    l=filter(i->all(c->mod(c[1],2)==0 || c[2]==1,
       tally(uc[:classes][i][:parameter])),eachindex(uc[:classes])) 
    for i in l
      cl=uc[:classes][i]
      s=LuSpin(cl[:parameter])
      if length(cl[:Au]) == 1
        cl[:Au] = coxgroup(:A, 1)
        trspringer(i, [1], [2])
        d = 1
      elseif length(cl[:Au]) == 4
        cl[:Au] = coxgroup(:B, 2)
        trspringer(i, [1, 2, 3, 4], [1, 3, 5, 4])
        d = 2
      else
        error("Au non-commutative of order ",Size(cl[:Au])*2," not implemented")
      end
      addSpringer(ss->ss[:Z]==[-1] && rank(ss[:relgroup])==sum(sum,s),i,s,d)
    end
  end
  for ss in uc[:springerSeries] if !all(isone,ss[:Z]) ss[:hc]=0 end end
  return uc
end)

chevieset(:B, :ClassParameter, function(n, w)
  cycletype(prod(i->i==1 ? SPerm(-1) : SPerm(i-1,i),w;init=SPerm()),n)
end)
