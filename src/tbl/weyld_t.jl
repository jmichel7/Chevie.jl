# Hand-translated part of chevie/tbl/weyld.g
# (C) Jean Michel & Goetz Pfeiffer 1994-2001
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
    uc[:orderClasses] = hasse(Poset(map(uc[:classes])do x
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
    end))
  else
    uc[:orderClasses]=hasse(Poset(map(eachindex(uc[:classes]))do i
                                  map(eachindex(uc[:classes]))do j
      dominates(uc[:classes][j][:parameter], uc[:classes][i][:parameter]) && 
      (uc[:classes][j][:parameter]!=uc[:classes][i][:parameter] || i==j)
                     end
                     end))
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
        append!(a,j+div(l+mod(l,4),4)-t)
        append!(b,j+div(l-mod(l,4),4)+t)
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
  uc
end)
