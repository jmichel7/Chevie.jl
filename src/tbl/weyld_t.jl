# Hand-translated part of chevie/tbl/weyld.g
# (C) Jean Michel, Frank Luebeck, Goetz Pfeiffer 1994-2001
# Data for type D

chevieset(:D,:CartanMat,n->toL(cartan(:D,n)))

#  If l is even then some of the classes  and restrictions split into two
#  classes or  characters, respectively. Their  values  are given by  the
#  character  values    for W(B_l) and   those  for  the  symmetric group
#  S_(l/2). This is described in [Pfeiffer, G., Character Tables of  Weyl
#  Groups in GAP]. 
chevieset(:D,:HeckeCharTable,function(n,para,root)
function chard(n,q)
  if n%2==0
    n1=div(n,2)-1
    AHk=chevieget(:A,:HeckeCharTable)(n1,fill([q^2,-1],n1),[])[:irreducibles]
    pA=partitions(n1+1)
    Airr(x,y)=AHk[findfirst(==(x),pA)][findfirst(==(y),pA)]
  end
  BHk=chevieget(:imp,:HeckeCharTable)(2,1,n,vcat([[1,-1]],fill([q,-1],n)),[])
  pB=chevieget(:B,:CharInfo)(n)[:charparams]
  Birr(x,y)=BHk[:irreducibles][findfirst(==(x),pB)][findfirst(==(y),pB)]
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
  [[value(lambda,mu) for mu in chevieget(:D,:ClassInfo)(n)[:classparams]]
   for lambda in chevieget(:D,:CharInfo)(n)[:charparams]] 
end
   u=-para[1][1]//para[1][2]
   tbl=Dict{Symbol,Any}(:name=>"H(D_$n)")
   tbl[:identifier]=tbl[:name]
   tbl[:parameter]=fill(u,n)
   tbl[:irreducibles]=chard(n,u)
   tbl[:size]=prod(chevieget(:D,:ReflectionDegrees)(n))
#  tbl[:irredinfo]=List(CHEVIE.R("CharInfo","D")(n).charparams,p->
#     rec(charparam:=p,charname:=PartitionTupleToString(p)));
   merge!(tbl,chevieget(:D,:ClassInfo)(n))
   CHEVIE[:compat][:AdjustHeckeCharTable](tbl,para);
   tbl
  end)
chevieset(:D,:CharTable,n->chevieget(:D,:HeckeCharTable)(n,fill([1,-1],n),[]))

chevieset(:D,:CycPolPoincarePolynomial,n->CycPol(Pol()^n-1)*
          prod(i->CycPol(Pol()^2i-1),1:n-1)//CycPol(Pol()-1)^n)

chevieset(:D, :SchurElement, function (n, phi, para, sqrtparam)
  Value(chevieget(:D, :CycPolPoincarePolynomial)(n)//
        chevieget(:D, :CycPolGenericDegree)(phi), -para[1][1]//para[1][2])
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
    else p=findfirst(==([s.sp]),CharParams(ss[:relgroup]))
    end
    ss[:locsys][p]=[i,findfirst(==(map(x->x ? [1,1] : [2],s.Au)),
                                                  CharParams(cc[:Au]))]
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
                         sum(sum,fullsymbol(x.sp))],c),ss))
  sort!(l, by=x->[abs(x[1]), -sign(x[1])])
  uc = Dict{Symbol, Any}(:classes => [], :springerSeries => map(function(d)
      res = Dict{Symbol, Any}(:defect=>d[1], :levi=>1:n-d[2])
      if mod(n-d[2],4)==0 || char==2 res[:Z]=mod(n,2)==0  ? [1, 1] : [1]
      else                           res[:Z]=mod(n,2)==0  ? [-1, -1] : [-1]
      end
      res[:relgroup]=coxgroup(d[1]==0 ? :D : :B, d[2])
      res[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(res[:relgroup])]
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
  if char == 2 # see [Spaltenstein] 2.10 page 24
    uc[:orderClasses] = hasse(Poset(toM(map(uc[:classes])do x
                                    map(uc[:classes])do y
      m = max(x[:parameter][1][1], y[:parameter][1][1])
      f = x-> map(i-> sum(filter(<(i),x))+i*count(>=(i),x) , 1:m)
      fx = f(x[:parameter][1])
      fy = f(y[:parameter][1])
      for i in 1:m
        if fx[i] < fy[i] return false
        elseif fx[i] == fy[i] && i in (y[:parameter])[2]
          if i in setdiff(x[:parameter][1], x[:parameter][2]) return false end
          if i<m && mod(fx[i+1]-fy[i+1],2)==1 return false end
        end
      end
      if x[:parameter] == y[:parameter] && x != y return false end
      return true
    end
    end)))
  else
    uc[:orderClasses]=hasse(Poset(toM(map(eachindex(uc[:classes]))do i
                                  map(eachindex(uc[:classes]))do j
      dominates(uc[:classes][j][:parameter], uc[:classes][i][:parameter]) && 
      (uc[:classes][j][:parameter]!=uc[:classes][i][:parameter] || i==j)
                     end
                    end)))
  end
  if char != 2
    d = 0
    while 4*d^2-d <= n
     i = 4*d^2-d
     if mod(n-d,2)==0
       l=vcat(1:i,i+2:2:n)
       s=Dict(:relgroup=>coxgroup(:B, div(n-i,2)),:levi=>l)
       if mod(n, 2) == 0 s[:Z]=[1, -1] else s[:Z] = [E(4)] end
       s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
       push!(uc[:springerSeries], s)
       if d==0 l=vcat([1], 4:2:n) end
       s=Dict(:relgroup => coxgroup(:B, div(n-i,2)),:levi=>l)
       if mod(n,2)==0 s[:Z] = [-1, 1]
       else s[:Z] = [-(E(4))]
       end
       s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
       push!(uc[:springerSeries], s)
       i=4d^2+d
       if d != 0 && i <= n
         l = vcat(1:i, i+2:2:n)
         s = Dict(:relgroup=>coxgroup(:B, div(n-i,2)),:levi=>l)
         if mod(n,2)==0 s[:Z]=[1, -1] else s[:Z]=[E(4)] end
         s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
         push!(uc[:springerSeries], s)
         s=Dict(:relgroup=>coxgroup(:B,div(n-i,2)),:levi=>l)
         if mod(n, 2) == 0 s[:Z] = [1, 1] else s[:Z] = [-(E(4))] end
         s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
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
    else p = findfirst(==([s]),CharParams(ss[:relgroup]))
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
         addSpringer1(ss->ss[:Z] in [[1, -1], [E(4)]] && 
                     rank(ss[:relgroup]) == sum(sum,s) , i, s, k[1])
     end
     if !('+' in cl[:name])
         addSpringer1(ss->ss[:Z] in [[-1, 1], [-(E(4))]] && 
                     rank(ss[:relgroup]) == sum(sum,s), i, s, k[2])
     end
  end
end
  for ss in uc[:springerSeries] if !all(isone,ss[:Z]) ss[:hc]=0 end end
  uc
end)

chevieset(:D, :gensMODA,Dict(4=>[[perm"(1,2)(7,8)", perm"(3,4)(5,6)", perm"(2,3)(6,7)", perm"(3,5)(4,6)"],[[4],[nothing, nothing, 2]], [[2], [1, nothing, 1]]],
6=>[[perm"(1,2)(8,11)(12,14)(15,17)(16,18)(19,21)(22,25)(31,32)", 
     perm"(3,4)(5,6)(7,9)(10,13)(20,23)(24,26)(27,28)(29,30)", 
     perm"(2,3)(6,8)(9,12)(13,16)(17,20)(21,24)(25,27)(30,31)", 
     perm"(3,5)(4,6)(12,15)(14,17)(16,19)(18,21)(27,29)(28,30)", 
     perm"(5,7)(6,9)(8,12)(11,14)(19,22)(21,25)(24,27)(26,28)", 
     perm"(7,10)(9,13)(12,16)(14,18)(15,19)(17,21)(20,24)(23,26)"], 
    [[16], [4, nothing, 6], [1, nothing, nothing, nothing, 5]], 
    [[12], [2, nothing, 6], [nothing, 2, nothing, nothing, 4]]],
8=>[[perm"(1,2)(8,11)(12,15)(16,20)(17,21)(22,26)(23,27)(28,33)(29,34)(30,35)(36,41)(37,42)(43,50)(44,51)(52,59)(60,68)(61,69)(70,77)(78,85)(79,86)(87,92)(88,93)(94,99)(95,100)(96,101)(102,106)(103,107)(108,112)(109,113)(114,117)(118,121)(127,128)", 
     perm"(3,4)(5,6)(7,9)(10,13)(14,18)(19,24)(25,31)(32,38)(39,45)(40,46)(47,53)(48,54)(49,55)(56,62)(57,63)(58,64)(65,71)(66,72)(67,73)(74,80)(75,81)(76,82)(83,89)(84,90)(91,97)(98,104)(105,110)(111,115)(116,119)(120,122)(123,124)(125,126)", 
     perm"(2,3)(6,8)(9,12)(13,17)(18,23)(20,25)(24,30)(26,32)(33,39)(34,40)(41,48)(42,49)(50,57)(51,58)(53,61)(59,67)(62,70)(68,76)(71,78)(72,79)(80,87)(81,88)(89,95)(90,96)(97,103)(99,105)(104,109)(106,111)(112,116)(117,120)(121,123)(126,127)", 
     perm"(3,5)(4,6)(12,16)(15,20)(17,22)(21,26)(23,29)(27,34)(30,37)(35,42)(39,47)(45,53)(48,56)(54,62)(57,65)(58,66)(63,71)(64,72)(67,75)(73,81)(76,84)(82,90)(87,94)(92,99)(95,102)(100,106)(103,108)(107,112)(109,114)(113,117)(123,125)(124,126)", 
     perm"(5,7)(6,9)(8,12)(11,15)(22,28)(26,33)(29,36)(32,39)(34,41)(37,44)(38,45)(40,48)(42,51)(46,54)(49,58)(55,64)(65,74)(71,80)(75,83)(78,87)(81,89)(84,91)(85,92)(88,95)(90,97)(93,100)(96,103)(101,107)(114,118)(117,121)(120,123)(122,124)", 
     perm"(7,10)(9,13)(12,17)(15,21)(16,22)(20,26)(25,32)(31,38)(36,43)(41,50)(44,52)(48,57)(51,59)(54,63)(56,65)(58,67)(62,71)(64,73)(66,75)(70,78)(72,81)(77,85)(79,88)(86,93)(91,98)(97,104)(103,109)(107,113)(108,114)(112,117)(116,120)(119,122)", 
     perm"(10,14)(13,18)(17,23)(21,27)(22,29)(26,34)(28,36)(32,40)(33,41)(38,46)(39,48)(45,54)(47,56)(52,60)(53,62)(59,68)(61,70)(67,76)(69,77)(73,82)(75,84)(81,90)(83,91)(88,96)(89,97)(93,101)(95,103)(100,107)(102,108)(106,112)(111,116)(115,119)", 
     perm"(14,19)(18,24)(23,30)(27,35)(29,37)(34,42)(36,44)(40,49)(41,51)(43,52)(46,55)(48,58)(50,59)(54,64)(56,66)(57,67)(62,72)(63,73)(65,75)(70,79)(71,81)(74,83)(77,86)(78,88)(80,89)(85,93)(87,95)(92,100)(94,102)(99,106)(105,111)(110,115)"], 
    [[64], [16, nothing, 24], [nothing, nothing, 32], [4, nothing, nothing, nothing, 20], [nothing, nothing, nothing, nothing, nothing, nothing, 16]], 
    [[56], [12, nothing, 24], [6, nothing, 28], [2, 4, nothing, nothing, 18], [1, nothing, 3, nothing, nothing, nothing, 14]]]))

chevieset(:D, :ClassParameter, function (n, w)
  x=prod(i->i==1 ? SPerm(1,-2) : SPerm(i-1,i),w;init=SPerm())
  res=cycletype(x,n)
  if isempty(res[2]) && all(iseven, res[1])
  # for classes with '+' or '-' we use the cycle type for the
  # permutation representation on the cosets of the parabolic
  # subgroup 2:n (gens(Dₙ) in this representation are 
  # stored for n<=8 in the Dict 'chevieget(:D,:gensMODA)')
    u=get!(chevieget(:D,:gensMODA),n)do
      H=coxgroup("D", n)
      R=reflection_subgroup(H, 2:n)
      D=vcat(reduced(R,H)...)
      Dgens=map(s->Perm(reduced.(Ref(R),D.*s),D),gens(H))
      ci=classinfo(H)
      H=ci[:classtext][filter(i->ci[:classnames][i][end] in "+-", 
                                 1:length(ci[:classnames]))]
      H=map(a->CycleStructurePerm(prod(Dgens[a])), H)
      u=[Dgens,H[1:2:length(H)], H[2:2:length(H)]]
    end
    tmp=CycleStructurePerm(prod(u[1][w]))
    if tmp in u[2] && !(tmp in u[3]) res=[res[1],'+']
    elseif !(tmp in u[2]) && tmp in u[3] res=[res[1],'-']
    end
  end
  sort!(res[1],rev=true)
  if res[2] isa Vector sort!(res[2],rev=true) end
  res
end)
