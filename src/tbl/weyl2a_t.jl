# Hand-translated part of chevie/tbl/weyl2a.g
# Frank Luebeck & Jean Michel (C) 1992--2001
# Data for the coset W(A_r).F where F induces -w0.

chevieset(Symbol("2A"),:CharTable,function(r)
  tbl = copy(chevieget(:A, :CharTable)(r))
  tbl[:identifier] = "W(^2A_$r)"
  A=chevieget(:A,:LowestPowerFakeDegree).(chevieget(:A,:CharInfo)(r)[:charparams])
  tbl[:irreducibles]= (-1).^A .* tbl[:irreducibles]
  merge!(tbl, chevieget(Symbol("2A"), :ClassInfo)(r))
end)

chevieset(Symbol("2A"),:FakeDegree,function(n,p,q)
  res=chevieget(:A, :FakeDegree)(n, p, Pol())
  (-1)^valuation(res)*res(-q)
end)

chevieset(Symbol("2A"),:HeckeCharTable, function(r, param, rootparam)
  q=-param[1][1]//param[1][2]
  v=rootparam[1]
  if isnothing(v) v=root(q,2) end
  tbl=Dict{Symbol,Any}(:identifier=>"H(^2A_$r)")
  merge!(tbl,chevieget(Symbol("2A"),:ClassInfo)(r))
  W=coxgroup(:A,r)
# If q_E is the square root which deforms to 1 of the eigenvalue of T_{w_0}
# on E which deforms to 1, then we have:
#  Ẽ(T_wφ)=τ(E(T_{w^-1w_0}))q_E (trivial extension)
#  Ẽ(T_wφ)=(-1)^a_E τ(E(T_{w^-1w_0}))q_E (preferred extension)
# where τ is q->q^-1
  H=hecke(W, v^-2)
  ct=CharTable(H)
  tbl[:charnames]=ct.charnames
  tbl[:centralizers]=ct.centralizers
  T=Tbasis(H)
  cl=map(x->T(W(x...)*longest(W)), tbl[:classtext])
  tbl[:irreducibles]=toL(transpose(toM(char_values.(cl))))
  charparams=chevieget(:A,:CharInfo)(r)[:charparams]
  A=chevieget(:A,:LowestPowerFakeDegree).(charparams)
  qE=central_monomials(hecke(W,v))
  A=(-1).^A .* qE
  tbl[:irreducibles]=map((x,y)->x.*y,A,tbl[:irreducibles])
  CHEVIE[:compat][:AdjustHeckeCharTable](tbl, param)
  tbl
end)

chevieset(Symbol("2A"), :ClassParameter, function (n, w)
  x=prod(i->Perm(i,i+1),w;init=Perm())*prod(i->Perm(i,n+2-i),1:div(n+1,2))
  cycletype(x;domain=1:n+1,trivial=true)
end)

chevieset(Symbol("2A"), :UnipotentClasses, function (r, p)
  uc=deepcopy(chevieget(:A, :UnipotentClasses)(r, p))
  for c in uc[:classes]
    t=parent(c[:red])
    if length(refltype(t))>1 error() end
    if length(refltype(t))==0 || rank(t)==1 WF=spets(parent(c[:red]))
    else WF=spets(parent(c[:red]),prod(i->Perm(i,rank(t)+1-i),1:div(rank(t),2)))
    end
    t=twistings(WF,inclusiongens(c[:red]))
    m=map(x->refleigen(x,position_class(x,x())),t)
    m=map(x->count(y->y==E(2),x),m)
    p=findmax(m)[2]
    c[:red]=spets(c[:red],t[p].F)
  end
  uc
end)

