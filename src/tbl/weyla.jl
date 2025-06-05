# Hand-translated from chevie/tbl/weyla.g
# Jean Michel & Goetz Pfeiffer (C) 1994--2001
# Data for W(A_n)

chevieset(:A, :CartanMat, n->Weyl.cartanmats[:A](n))

chevieset(:A, :ReflectionDegrees, n->2:n+1)

chevieset(:A, :AuName, function(r)
  ["Z_2", "S_3", "S_4", "S_5"][r]
end)

# roots for GL_n
chevieset(:A, :simpleroots, function(n)
  r=fill(0,n,n+1)
  for i in 1:n r[i,i:i+1].=[1,-1] end
  r
end)

chevieset(:A, :ParabolicRepresentatives,function(n,s)
  filter(x->!(n+1 in x),chevieget(:imp, :ParabolicRepresentatives)(1,1,n+1,s))
end)

# very good representative for partition pi in the sense of [Geck-Michel]
chevieset(:A,:WordClass,function(pi)
  w=Int[]
  i=0
  for l in pi
    r=mod(l,2)
    append!(w,i.+vcat(1:2:l-1-r,2:2:l+r-2))
    i+=l
  end
  w
end)

chevieset(:A,:centralizer,function(n,partition)
  res=k=1;last=0
  for p in partition
    res*=p
    if p==last k+=1;res*=k
    else k=1
    end
    last=p
  end
  res
end)

chevieset(:A,:ClassInfo,function(n)
  res=chevieget(:imp,:ClassInfo)(1,1,n+1)
  pp=partitions(n+1)
  res[:classparams]=pp
  res[:classtext]=map(chevieget(:A,:WordClass),pp)
  res
end)

chevieset(:A,:NrConjugacyClasses,n->npartitions(n+1))

chevieset(:A,:CharTable,function(n)
  ct=chevieget(:imp,:CharTable)(1,1,n+1)
  ct[:charnames]=chevieget(:A,:CharInfo)(n)[:charnames]
  ct
end)

chevieset(:A,:HeckeCharTable,
function(n,para,root)
  if n==1 Dict(:irreducibles=>[1 para[1][2];1 para[1][1]],
               :charnames=>["11","2"],
               :classnames=>["11","2"],
               :centralizers=>[2,2],
               :identifier=>"H(A_1)")
  else  chevieget(:imp,:HeckeCharTable)(1,1,n+1,para,root)
  end
end)

chevieset(:A,:FakeDegree,(n,p,q)->fakedegree(CharSymbol([βset(p)]))(q))

# partition for a coxeterword w
chevieset(:A, :ClassParameter, function(n,w)
  cycletype(prod(i->Perm(i,i+1),w;init=Perm()),1:n+1)
end)

chevieset(:A, :WeightInfo, function(n)
  M=Matrix(1I,n,n)
  for i in 1:n-1 M[i,n]=-i end
  Dict{Symbol, Any}(:minusculeWeights => 1:n, 
    :decompositions=>map(i->[n+1-i],1:n),:moduli=>[n+1],:chosenAdaptedBasis=>M)
end)

chevieset(:A,:b,p->dot(p,0:length(p)-1))

chevieset(:A,:B,p->binomial(sum(p),2)-sum(i->binomial(i,2),p))

chevieset(:A, :CharInfo, function(n)
  pp=partitions(n+1)
  res=Dict{Symbol, Any}(:charparams=>pp,
    :charnames=>joindigits.(pp),
    :extRefl=>map(i->findfirst(==(pushfirst!(fill(1,i),n+1-i)),pp),0:n),
    :b=>chevieget(:A,:b).(pp),
    :B=>chevieget(:A,:B).(pp))
  res[:a]=res[:b]
  res[:A]=res[:B]
  res
end)

chevieset(:A, :PoincarePolynomial,function(n,param)
  prod(i->sum(k->(-param[1][1]//param[1][2])^k,0:i),1:n)
end)

# Reference: Carter II, 13.5, p. 446.
chevieset(:A, :SchurElement, function (n, alpha, param, sqrtparam)
  q=-param[1][1]//param[1][2]
  lambda=βset(alpha)
  res=q^binomial(length(lambda), 3)
  for i in lambda
    for j in 0:i-1
      if j in lambda res//=q^j
      else res*=sum(e->q^e,0:i-j-1)
      end
    end
  end
  res
end)

chevieset(:A, :FactorizedSchurElement,function(n,alpha,para,arg...)
  chevieget(:imp, :FactorizedSchurElement)(1,1,n+1,[alpha],para,[])
end)

chevieset(:A, :HeckeRepresentation, function (n, param, sqrtparam, i)
  H=hecke(coxgroup(:A,n),-param[1][1]//param[1][2])
  -param[1][2]*Spechtmodel(H,partitions(n+1)[i])
end)

chevieset(:A, :Representation, function (n, i)
  chevieget(:A, :HeckeRepresentation)(n,fill([1,-1],n),fill(1,n),i)
end)

chevieset(:A, :FakeDegree, function (n, p, q)
  improve_type(exactdiv(chevieget(:A, :PoincarePolynomial)(sum(p)-1,[[q,-1]]),
                        chevieget(:A, :SchurElement)(sum(p)-1,p,[[q,-1]],[])))
end)

chevieset(:A, :KLeftCellRepresentatives,function(n)
  pp=partitions(n+1)
  function f(i)
    i=prod(j->Perm(j,j+1),word(W,i);init=Perm())
    p=length.(robinson_schensted((1:n+1).^i)[1])
    [findfirst(==(p),pp)]
  end
  W=coxgroup(:A,n)
  l=filter(x->isone(x^2),elements(W))
  map(x->Dict{Symbol, Any}(:duflo=>(1:n).^x,:reps=>[],:character=>f(x)),l)
end)

chevieset(:A, :DecompositionMatrix, function(n,p)
  [[1:npartitions(n+1),MatrixDecompositionMatrix(DecompositionMatrix(Specht(p,
    p),n+1))]]
end)

chevieset(:A, :UnipotentCharacters,function(n)
  ci=chevieget(:A,:CharInfo)(n)
  pp=ci[:charparams]
  Dict{Symbol,Any}(:harishChandra =>[Dict{Symbol, Any}(:levi=>Int[],
    :relativeType=>TypeIrred(series=:A,indices=1:n,rank=n), 
    :parameterExponents=>fill(1,n),:cuspidalName=>"",:eigenvalue=>1,
    :charNumbers=>1:length(pp))],
  :families=>map(i->Family("C1",[i]), 1:length(pp)), 
  :charParams=>pp,:charSymbols=>map(x->CharSymbol([βset(x)]),pp),
  :a=>ci[:a],:A=>ci[:A])
end)

chevieset(:A, :Ennola, n->n==1 ? SPerm([-1, 2]) : SPerm())

chevieset(:A, :Invariants, function(n)
  m=chevieget(:A,:simpleroots)(n)
  map(2:n+1)do i
    function(arg...)
      v=sum(arg.*eachrow(m))
      sum(a->prod(v[a]),arrangements(1:n+1,i))
    end
  end
end)

chevieset(:A, :UnipotentClasses, function (n, p)
  uc=Dict{Symbol, Any}(
    :classes=>map(p->Dict{Symbol,Any}(:parameter=>p),partitions(n+1)),
    :springerSeries=>vcat(map(d->map(i->
      Dict{Symbol,Any}(:relgroup =>coxgroup(:A,div(n+1,d)-1),:Z=>[E(d,i)],
        :levi=>filter(i->mod(i,d)!=0,1:n+1),:locsys=>[]),
            prime_residues(d)),divisors(n+1))...))
  ss(z)=uc[:springerSeries][findfirst(x->x[:Z]==[z],uc[:springerSeries])]
  function partition2parab(p)
    res=Int[]
    c=1
    for pa in p
      for i in 1:pa-1
        push!(res,c)
        c+=1
      end
      c+=1
    end
    res
  end
  for i in 1:length(uc[:classes])
    cl = uc[:classes][i]
    p = cl[:parameter]
    d = gcd(p)
    cl[:name] = joindigits(p)
    cl[:Au] = crg(d, 1, 1)
    cl[:balacarter]=vcat(map(i->sum(p[1:i-1]).+(1:p[i]-1),1:length(p))...)
    p=vcat(map(x->1-x:2:x-1,p)...)
    sort!(p)
    cl[:dynkin]=map(i->p[i+1]-p[i],1:length(p)-1)
    cl[:red]=coxgroup()
    p=tally(cl[:parameter])
    for j in p cl[:red]*=coxgroup(:A,j[2]-1) end
    cl[:red]*=torus(length(p)-1)
    cl[:AuAction]=ExtendedReflectionGroup(cl[:red],
                                         reflrep(cl[:red],one(cl[:red])))
    cl[:rep]=partition2parab(cl[:parameter])
    if d==2
      push!(ss(1)[:locsys], [i, 2])
      push!(ss(-1)[:locsys], [i, 1])
    else
      for j in 0:d-1 push!(ss(E(d,j))[:locsys],[i,j+1]) end
    end
  end
  for ss in uc[:springerSeries][2:end] ss[:hc]=0 end
  uc[:orderClasses]=hasse(Poset((x,y)->dominates(y,x),
                                map(x->x[:parameter],uc[:classes])))
  uc
end)
