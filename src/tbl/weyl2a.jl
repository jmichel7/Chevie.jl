# Hand-translated part of chevie/tbl/weyl2a.g
# Frank Luebeck & Jean Michel (C) 1992--2001
# Data for the coset W(A_r).F where F induces -w0.

chevieset("2A",:WordsClassRepresentatives,function(n,part=partitions(n+1))
  function redw(n, w)local l
    l=Int[]
    while true
      i=findfirst(j->j^w>(j+1)^w,1:n)
      if isnothing(i) return l end
      push!(l,i)
      w=Perm(i,i+1)*w
    end
    l
  end
  function guesslongest(p) # longest element in class, see [Geck-Kim]
    p=vcat(filter(iseven, p), filter(i->i!=1 && isodd(i), p))
    x=Perm()
    off=0
    for i in p
      x*=prod(j->Perm(l[off+1],l[off+j]), 2:i)
      off+=i
    end
    x
  end
  l=Int[]
  w0=Perm()
  for p in 1:div(n+1,2)
    append!(l,[p,n-p+2])
    w0*=Perm(p,n-p+2)
  end
  if iseven(n) push!(l,div(n,2)+1) end
  map(p->redw(n,guesslongest(p)*w0),part)
end)

chevieset("2A", :ClassInfo, function (n,)
  res=chevieget(:A,:ClassInfo)(n)
  res[:classtext]=chevieget("2A",:WordsClassRepresentatives)(n,res[:classparams])
  delete!(res, :orders)
  res
end)

chevieset("2A",:NrConjugacyClasses,n->npartitions(n+1))

chevieset("2A",:CharInfo,n->chevieget(:A,:CharInfo)(n))

chevieset("2A",:PhiFactors,n->map(x->(-1)^x,2:n+1))

chevieset("2A",:CharTable,function(r)
  tbl=copy(chevieget(:A,:CharTable)(r))
  tbl[:identifier]="W(^2A_$r)"
  A=chevieget(:A,:b).(chevieget(:A,:CharInfo)(r)[:charparams])
  tbl[:irreducibles]=Diagonal((-1).^A)*tbl[:irreducibles]
  merge!(tbl,chevieget("2A",:ClassInfo)(r))
end)

chevieset("2A",:FakeDegree,function(n,p,q)
  res=chevieget(:A, :FakeDegree)(n, p, Pol())
  (-1)^valuation(res)*res(-q)
end)

chevieset("2A",:HeckeCharTable, function(r, param, rootparam)
  q=-param[1][1]//param[1][2]
  v=rootparam[1]; if ismissing(v) v=root(q) end
  tbl=Dict{Symbol,Any}(:identifier=>"H(^2A_$r)")
  merge!(tbl,chevieget("2A",:ClassInfo)(r))
  W=coxgroup(:A,r)
# If q_E is the square root which deforms to 1 of the eigenvalue of T_{w_0}
# on E which deforms to 1, then we have:
#  Ẽ(T_wφ)=τ(E(T_{w^-1w_0}))q_E (trivial extension)
#  Ẽ(T_wφ)=(-1)^a_E τ(E(T_{w^-1w_0}))q_E (preferred extension)
# where τ is q->q^-1
# Here we take the preferred extension
  H=hecke(W, v^-2)
  ct=CharTable(H)
  tbl[:charnames]=ct.charnames
  tbl[:centralizers]=ct.centralizers
  T=Tbasis(H)
  cl=map(x->T(W(x...)*longest(W)), tbl[:classtext])
  tbl[:irreducibles]=transpose(toM(char_values.(cl)))
  charparams=chevieget(:A,:CharInfo)(r)[:charparams]
  a=chevieget(:A,:b).(charparams)
  qE=central_monomials(hecke(W,v))
  tbl[:irreducibles]=Diagonal((-1).^a .* qE)*tbl[:irreducibles]
  AdjustHeckeCharTable(tbl, param)
end)

chevieset("2A", :HeckeRepresentation, function (n, param, sqrtparam, i)
  W=coxgroup(:A,n)
  H=hecke(W,-param[1][1]//param[1][2])
  p=partitions(n+1)[i]
  gens=Spechtmodel(H,p)
  (gens=gens,F=prod(gens[word(W,longest(W))]).//
    root(central_monomials(H)[i])*(-1)^chevieget(:A,:b)(p))
end)

chevieset("2A", :Representation, function(n,i)
  chevieget("2A",:HeckeRepresentation)(n,map(x->[1,-1],1:n),1:1,i)
end)

# partition associated to 2-core d and partitions-pair p
function PartitionTwoCoreQuotient(d,p)
  x=Symbol_partition_tuple(reverse(p),-d).S
  partβ(sort(unique(vcat(x[1].*2,x[2].*2 .+1))))
end

chevieset("2A", :ClassParameter, function (n, w)
  x=prod(i->Perm(i,i+1),w;init=Perm())*prod(i->Perm(i,n+2-i),1:div(n+1,2))
  cycletype(x,1:n+1)
end)

# [LuB, 4.4, 4.16, 4.19]
chevieset("2A", :UnipotentCharacters, function (l,)
  uc=chevieget(:A,:UnipotentCharacters)(l)
  uc[:charSymbols] = map(i->[i],partitions(l+1))
  uc[:almostHarishChandra] = uc[:harishChandra]
  uc[:almostHarishChandra][1][:relativeType]=
  TypeIrred(;orbit=[TypeIrred(;series=:A,indices=1:l,rank=l)],
    twist=prod(i->Perm(i,l+1-i),1:div(l,2)),rank=l)
  uc[:harishChandra]=[]
  d=0
  while d*(d+1)//2<=l+1
    k=l+1-div(d*(d+1),2)
    if iseven(k)
      r=div(k, 2)
      s=Dict{Symbol, Any}(:levi=>r+1:l-r,:relativeType=> 
        TypeIrred(;series=:B,indices=r:-1:1,rank=r),
        :eigenvalue=>(-1)^div(prod(d.+(-1:2)),8))
      # for the eigenvalue see Lusztig CBMS proof of 3.34 (ii)
      if d==0 
       s[:relativeType]=TypeIrred(;series=:B,indices=r:-1:1,rank=r,cartanType=1)
      end
      if r!=0 s[:parameterExponents]=vcat(2d+1,fill(2,r-1))
      else s[:parameterExponents]=[]
      end
      if k<l s[:cuspidalName]="{}^2A"*stringind(rio(TeX=true),l-k)
      else s[:cuspidalName]=""
      end
      # see Fong/Srinivasan for this map
      s[:charNumbers]=map(a->findfirst(==([PartitionTwoCoreQuotient(d, a)]),
              uc[:charSymbols]),partition_tuples(r,2))
      FixRelativeType(s)
      push!(uc[:harishChandra],s)
    end
    d+=1
  end
  # for delta see Lusztig's book page 124 line 7
  for i in 1:length(uc[:families])
    if isodd(uc[:a][i]+uc[:A][i])
      uc[:families][i]=Family("C'1",uc[:families][i][:charNumbers])
    end
  end
  uc
end)

chevieset("2A", :UnipotentClasses, function (r, p)
  uc=copy(chevieget(:A, :UnipotentClasses)(r, p))
  uc[:classes]=map(uc[:classes])do c
    t=refltype(c[:red])
    if isempty(t) m=reflrep(c[:red],one(c[:red]))
    else m=reflrep(c[:red],Perm(vcat(map(x->reverse(indices(x)),t)...)))
    end
    p=tally(c[:parameter])
    for i in 1:length(p)-1 m[end+1-i,end+1-i]=-1 end
    c=copy(c)
    c[:red]=spets(c[:red],m)
    c
  end
  uc
end)

