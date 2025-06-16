# Programs to get the list of ud (Unipotent Degrees) + Frobenius eigenvalues
# as outlined in ``Split Spetses...''
using Chevie
include(replace(pathof(Chevie),"src/Chevie.jl"=>"test/Tests.jl")) # for EigenAndDegHecke

"`asranges(l)` string describing as list of ranges sorted integer list l"
function asranges(l)
  if isempty(l) return "[]" end
  res=[]
  br=l[1]
  er=br
  for i in 2:length(l)
    if l[i]!=er+1 push!(res,br:er);br=l[i];er=br
    else er+=1
    end
  end
  push!(res,br:er)
  join(map(x->x.start<x.stop ? string(x) : string(x.start),res),",")
end

## returns the list of [a,A] not fully accounted by the unipotent
## characters of indices l
## Test for completion is \sum_l|ud|^2=\sum_l|fd|^2
## The fds are read from harishChandra[1] but could be from the knowledge
## of Rouquier families
#UnaccountedAa:=function(W,l)local ud,uw,aA,h,q,fd,remain;
#  uw:=UnipotentCharacters(W);ud:=CycPolUnipotentDegrees(W);
#  aA:=TransposedMat([uw.a,uw.A]);
#  h:=uw.harishChandra[1].charNumbers;
#  aA:=List(Set(aA{h}),x->rec(a:=x[1],A:=x[2],
#     fd:=Filtered([1..Length(h)],i->aA[h[i]]=x),ud:=Filtered(l,i->aA[i]=x)));
#  q:=X(Cyclotomics);q.name:="q";fd:=FakeDegrees(W,q);
#  remain:=Filtered(aA,function(r)local p;
#    p:=Sum(fd{r.fd},x->x^2)-Sum(ud{r.ud},x->Value(x*ComplexConjugate(x),q));
#    return p<>0*p;end);
#  return List(remain,x->[x.a,x.A]);
#end;

"`eW(W)`: coxnum(W)*rank(W)"
function eW(W) 
  if W isa Spets eW(Group(W)) end
  nhyp(W)+nref(W)
end

"""
`RationalMinPolParam(p::parameter_list)`

A *parameter list* is a list of pairs (f::Root1,d::Rational)
representing a list of parameters f*q^d for an Hecke algebra.

whether `p` gives for prod(x-param) a polynomial in x,q.
That is for each non integral d we must have sum(corresponding f)=0.
"""
RationalMinPolParam(p)=all(v->sum(first.(filter(x->x[2]==v,p)))==0,
                           filter(!isinteger,unique(last.(p))))

"""
`CycPolCyclicSchur(p::parameter_list,i)`

A *parameter list* is a list of pairs (f::Root1,d::Rational)
representing a list of parameters f*q^d for an Hecke algebra.

returns CycPol(i-th schur elem) of Hecke(crg(length(p),1,1),p)
"""
function CycPolCyclicSchur(p,i)
  p=collectby(last,map(x->(p[i][1]//x[1],p[i][2]-x[2]),deleteat!(copy(p),i)))
  p=map(x->(last(x[1]),first.(x)),p)
  prod(p)do l
    if isinteger(l[1])
      return prod(x->CycPol(1-Pol([x],Integer(l[1]))),l[2])
    end
    v=prod(x->1-x*Mvp(:x)^l[1],l[2])
    if any(!isinteger,degree.(monomials(v))) zero(CycPol)
    else CycPol(v)
    end
  end
end

"""
`solveit(val,pos,other,positions)`
arrangements of val completed by an arrangement of other with in the
resulting arrangement val[i] positioned at an element of pos[i].
this is recursive and positions just serves for recursion.
"""
function solveit(val,pos,other,positions=0:length(val)+length(other))
  if isempty(val) 
   return map(a->collect(zip(a,positions)),arrangements(other,length(other)))
  end
  p=pos[1]
  pp=filter(i->pos[i]==p,eachindex(val))
  v=val[pp]
  res=Vector{NTuple{2,Int}}[]
  for a in arrangements(other,length(p)-length(pp))
    v1=vcat(v,a)
    head=map(x->collect(zip(x,p)),arrangements(v1,length(v1)))
    tail=solveit(deleteat!(copy(val),pp),
                 deleteat!(copy(pos),pp),msetdiff(other,a),setdiff(positions,p))
    append!(res,map(x->vcat(x...),cartesian(head,tail)))
  end
  res
end

"""
`RouquierFamilies(W)`
Returns (partition of `1:|Irr(W)|` in Rouquier families,
         list of corresponding `a+A`)
this could be obtained in theory from Maria's program with no
knowledge of `UnipotentCharacters`.
"""
function RouquierFamilies(W)
  uc=UnipotentCharacters(W)
  ff=map(f->indexin(charnumbers(f),uc.harishChandra[1][:charNumbers]),
         uc.families)
  ff=map(x->Int.(filter(!isnothing,x)),ff)
  fd=fakedegrees(W)
  (ff,map(x->minimum(valuation.(fd[x]))+maximum(degree.(fd[x])),ff))
end

"""
`decomposeassum(l,card,max)`
Let  l  be  a  multiset  (list  with  repetitions)  of length prod(card) of
rationals in [0,sum(max)]. If it is possible to write
l=sum.(cartesian(S₁,S₂))  where Sᵢ is a list of length card[i] of rationals
with extremas [0,max[i]] then returns S₁, S₂ else gives error.
"""
function decomposeassum(l,card,maxs)
  sort!(l)
  t=cartesian(map(i->filter(x->extrema(x)==(0,maxs[i]),
                             combinations(l,card[i])),1:2)...)
  t=filter(i->tally(sum.(cartesian(i...)))==tally(l),t)
  if length(t)!=1 error(length(t)," decompositions as sum") end
  only(t)
end

principal_series(W,d;opt...)=
  Series(W,only(cuspidal_data(W, d, length(relative_degrees(W, d))))...;opt...)

"""
`degparam(W,ζ,knownud=CycPol[])`  where  `ζ`  is  a  regular eigenvalue and
`knownud` are known unipotent degrees.

returns degrees of parameters for `H_W(ζ)` for each Hplane orbit in `W(ζ)`,
each  list encoded as `[(d,n),..]`  with `d∈ℚ` and `n`  a majoration of the
multiplicity of `d`

A better majoration of the multiplicity is obtained if `knownud` is nonempty.
We use formula 5.22 in "Split Spetses":
``∑_{χ∈ F}|Feg(χ)(ζ)|^2=∑_{ψ∈  W_d|ρ_ψ∈ F}ψ(1)^2``
and that `ρ(ζ)=0` unless `ρ=ρ_ψ` and `ψ(1)` else.

We need that `W(ζ)` has `<=2` Hplane orbits. An exact list is returned if 2
orbits  Hplanes  (assuming  all  `ρ_ψ`  for  ψ  linear known) or if W(ζ) is
cyclic.
"""
function degparam(W,ζ::Root1,knownud=CycPol[])
  aA=map(x->(eW(W)-degree(x)-valuation(x),x(ζ)^2),knownud)
  lin=filter(x->x[2]==1,aA)
  s=principal_series(W,ζ) #we use only s.levi and relative_group(s)
  if length(s.levi)>1 error(ζ," is not regular") end
  ho=hyperplane_orbits(relative_group(s))
  if length(ho)>1
    if length(lin)!=prod(x->x.order,ho)
      error("unimplemented: only ",length(lin)," out of ",
      prod(x->x.order,ho)," linear characters known and >=2 Hplane orbits\n")
      return nothing
    end # else all linear chars known
    maxs=map(o->eW(reflection_subgroup(W,s.WGL.reflists[o.s]))*o.N_s,ho)
    maxs=decomposeassum(first.(lin),map(x->x.order,ho),maxs)
    return map(i->maxs[i]//ho[i].order//ho[i].N_s,eachindex(ho))
  end
# Print("maxdeg=",eW(s.WGL)/ho[1].order,"\n");
  ff,faA=RouquierFamilies(W)
  fd=fakedegrees(W,Pol())
  n=map(x->sum(y->y(ζ)*y(ζ^-1),fd[x]),ff)
  n=collect(zip((eW(W).-faA)//eW(s.WGL),n))
  n=filter(x->last(x)!=0 && first(x)>=0,n)
  n=map(y->(y[1][1],Int(sum(last,y))),collectby(first,n))
  if iscyclic(s) return n end #cyclic case
  aA=filter(x->last(x)>1 && first(x)>=0,aA)
  aA=map(y->(y[1][1],sum(last,y)),collectby(first,aA))
  aA=map(x->(first(x)//eW(s.WGL),last(x)),aA)
  for v in aA
    if last(v)>1 
      p=findfirst(x->x[1]==first(v),n)
      n[p]=(n[p][1],n[p][2]-last(v))
    end
  end
  sort(filter(x->x[2]!=0 && x[1]>=0,n))
end

"subtract multiset s from multiset l"
function msetdiff(l,s)
  l=copy(l)
  for e in s
    p=findfirst(==(e),l)
    if isnothing(p) error("$s is not a sub-multi-set of $l\n") end
    l=deleteat!(l,p)
  end
  l
end

"""
`FindRelativeHecke(W,ζ,known)`
returns the Hecke algebra of ζ-principal series of W: ζ should be a regular
eigenvalue such that W_ζ has <=2 hplane orbits 
known: list of indices of already known (deg,eig) of unipotent characters
in case W_ζ has 2 hplane orbits we need known to include indices of all 
unipotent chars corresponding to linear chars of W_ζ
"""
function FindRelativeHecke(W,ζ,known)
  s=principal_series(W,ζ)
  if length(s.levi)>1 error(ζ," is not regular\n") end
  xprintln("find Hecke parameters for ",ζ,"-series W_G(L)=",relative_group(s))
  ho=hyperplane_orbits(s.WGL)
  uw=UnipotentCharacters(W)
  ud=CycPoldegrees(uw)
  sch=map(x->sign(Int(x(ζ)))*degree(s)//x,ud)
  # known parameters <-> linear chars of WGL corresp to known => ser.ind
  sort!(known)
  lin=filter(i->abs(Int(ud[i](ζ)))==1,known)
  aA=uw.a[lin]+uw.A[lin]
  ser=(ind=lin,aA=aA,eigen=eigen(uw)[lin],
       degparam=(eW(W).-aA)//(ho[1].order*ho[1].N_s),
       specialization=map(x->Int[],lin))
  p=degparam(W,ζ,ud[known]) # list of [degparam, majoration of multiplicity]
  if isnothing(p) return p end
  if iscyclic(s)
    println("\t",length(ser.ind)," out of ",s.e," schur elements known")
  # In the cyclic case, for a parameter x=ζₑⁱ(ζ⁻¹ q)^mᵢ,
  # the Frobenius eigenvalue Fr is ζ^{(e_W-em_i)*d+i}, 
  # so Fr+em_i*d^2-e_W*d^2=i*d (mod 1) or if d=num/den then 
  # Fr*num'*den+(em_i-e_W)*d=i (mod den) where num'=num^(-1) (mod den)
  # ser.specialization contains the possible values of i for a param
    if isone(ζ) 
      ser.specialization.=map(x->(0:s.e-1).+(x%1),ser.degparam*s.e.-eW(W))
    else
      d=ζ.r
      ser.specialization.=
        map(x->(0:denominator(d):s.e-denominator(d)).+mod(x,denominator(d)),
            map(x->x.r,ser.eigen)*denominator(d)*
        invmod(numerator(d),denominator(d)).+ser.degparam*s.e*d.-eW(W)*d)
    end
    poss=findfirst(==(findfirst(==(CycPol(1,0)),ud)),ser.ind)
    ser.specialization[poss]=[0]
    poss=setdiff(eachindex(ser.specialization),[poss])
    ser.specialization[poss]=map(x->filter(!iszero,x),ser.specialization[poss])
#   tbl:=TransposedMat([ser.eigen,ser.degparam,ser.specialization]);
#   SortBy(tbl,x->x[3]);
#   Print("Known part of ",d,"-series\n",FormatTable(tbl,
#      rec(rowLabels:=ser.ind,columnLabels:=["eig","degparam","specialize"])));
    p=vcat(map(x->fill(x[1],x[2]),p)...)
    p=msetdiff(p,ser.degparam)
    if isempty(p) # all are known
      p=map(i->findfirst(==(i),sch[ser.ind]),sch[ser.ind])
      poss=dSeries.FitParameter(map(x->subs(x,Pol([ζ],1)),sch[ser.ind]),ser.degparam)
      poss=vcat(map(x->collect(zip(x...)),poss)...)
      sort!(poss,by=last)
      poss=map(y->ser.degparam[y[1]],poss)
      p=argmax(poss)
      if p!=1 poss=poss[sortperm(map(x->x.r,E(s.e)^(1-p)*E.(s.e,0:s.e-1)))] end
      poss=[poss]
    else poss=solveit(ser.degparam,ser.specialization,p,0:s.e-1)
      poss=map(poss)do v
        res=fill(0//1,maximum(last.(v))+1)
        for x in v res[x[2]+1]=x[1] end
        return res 
      end
    end
   # get pairs [root of unity part of ith param, m_i]
    poss=filter(x->maximum(x)==x[1],poss)
    if length(poss)>1 println("\t",length(poss)," parameter lists to test") end
    cnt=0
    poss=map(poss)do l
      cnt+=1
      if cnt%1000==0 print(".") end
      map(i->(E(s.e,i-1)//ζ^l[i],l[i]),eachindex(l))
    end
    poss=filter(RationalMinPolParam,poss)
    if length(poss)>1
      println("\t",length(poss)," give an algebra defined over ℂ [q]")
    end
    cnt=0
    poss=filter(poss)do x
      cnt+=1
      if cnt%100==0 print(".") end
      CycPolCyclicSchur(x,1) in sch
    end
    if length(poss)>1 error("more than one poss") end
    H=hecke(crg(s.e,1,1),[map(x->x[1]*Mvp(:q)^x[2],only(poss))])
    xprintln("found ",H)
    p=CycPol.(schur_elements(H))
    if length(p)==2 
      p=reverse(p)
    end
    if any(i->!(sch[ser.ind[i]] in p[ser.specialization[i].+1]),
           eachindex(ser.ind))
      cmpvec(sch[ser.ind],p[first.(ser.specialization).+1])
      error("expected not in schur")
    end
  elseif length(ho)==1 # one orbit of Hplanes
    p=vcat(map(x->fill(x[1],min(x[2],ho[1].order-length(ser.eigen))),p)...)
    p=filter(x->sum(x)+sum(ser.degparam)==nref(W)//ho[1].N_s,
             combinations(p,ho[1].order-length(ser.eigen)))
    p=map(x->vcat(ser.degparam,x),p)
  # Print("p=",p,"\n");
  # here each line of p is possible degparam for other parameters
    res=[]
    pp=charinfo(s.WGL).charparams[vcat([charinfo(s.WGL).positionId],ho[1].det_s)]
    for c in p
      pars=map(a->map(i->E(ho[1].order,i)//ζ^a[i+1]*Mvp(:q)^a[i+1],
                      0:ho[1].order-1),arrangements(c,ho[1].order))
      print("\t",length(pars)," to try:")
      for par in pars
        print(".")
        H=hecke(s.WGL,[par])
        ss=map(x->HeckeAlgebras.schur_element(H,x),pp)
        if all(x->all(isinteger,degree.(monomials(x))),ss)
          ss=CycPol.(ss)
          if all(i->i in sch,ss) push!(res,par) end
        end
      end
    end
    println(length(res)," selected")
    res=filter(res)do par 
      H=hecke(s.WGL,[par])
      ss=CycPol.(schur_elements(H)) 
      all(i->iszero(i) || i in ss,sch)
    end
    # Print("possible params:",res,"\n");
    # test that res[1] has highest degree param monic
    res=filter(res)do p
      p=filter(x->degree(x)==maximum(degree.(p)),p)
      if length(p)>1 error("more that one param of max. degree") end
      only(p)(q=ζ)==1
    end
    if length(res)!=1 error("possibilities:",res) end
    H=hecke(s.WGL,res)
  else # for Length(ho)>1 we assume p contains correct list of degparam
    p=reverse.(p)
    pars=cartesian(map(i->map(a->vcat([p[i][1]],a),
           arrangements(p[i][2:ho[i].order],ho[i].order-1)),eachindex(ho))...)
    print("\t",length(pars)," to try:")
    pars=map(x->map(y->Int.(y),x),pars)
    pars=map(x->map(j->map(i->E(ho[j].order,i)/ζ^x[j][i+1]*Mvp(:q)^x[j][i+1],
                           0:ho[j].order-1),eachindex(ho)),pars)
    print("\t",length(pars)," to try:")
    res=[]
    pp=[charinfo(s.WGL).charparams[charinfo(s.WGL).positionId]]
    for par in pars
      print(".")
      i=fill([],ngens(s.WGL))
      for j in eachindex(par) 
        i[ho[j].s]=par[j]
      end
      par=i
      H=hecke(s.WGL,i)
      ss=map(x->HeckeAlgebras.schur_element(H,x),pp);
      if all(x->all(isinteger,degree.(monomials(x))),ss)
	ss=CycPol.(ss)
        if all(i->i in sch,ss) push!(res,par) end
      end
    end
    if length(res)>1 error("possibilities:",res) end
    H=hecke(s.WGL,only(res))
  end
  H
end

"""
`ennola_twist(H::HeckeAlgebra,ξ::Root1,ζ::Root1)`
Given a `ζ`-cyclotomic algebra `H` with parameters monomials in `q`,
Ennola-twist it by `ξ`, getting a `ζξ`-cyclotomic algebra
"""
function ennola_twist(H::HeckeAlgebra,z::Root1,d::Root1)
  zeta=d*z
  q=Mvp(:q)
  z=map(y->map(x->x(;q=q//z),y),H.para)
  H=hecke(H.W,map(p->sort(p,by=x->Root1(scalar(x(q=zeta))).r),z))
# xprintln("by ",z,"-twisting:",H)
  H
end

# positions-with-sign in l where o or -o appears
positionssgn(l,o)=vcat(findall(==(o),l),-findall(==(-o),l))

# find positions in Uch(W) of elts of list il of (deg=+-xx,eig=xx,frac=xx)
# s is a Serie; fill s.charnumbers and s.ambig in consequence
function finddeinW(W,s,il)
  l=copy(il)
  uw=UnipotentCharacters(W)
  ud=CycPoldegrees(uw)
  res=fill(0,length(il))
  eigs=eigen(uw)
  pos=eachindex(l)
  while length(pos)>0
    r=l[pos[1]]
    m=filter(i->l[i].deg in (r.deg,-r.deg),pos)
    if haskey(r,:eig) m=filter(i->l[i].eig==r.eig,m) end
    p=abs.(positionssgn(ud,r.deg))
    if haskey(r,:eig) p=filter(y->eigs[y]==r.eig,p) end
    if length(p)!=length(m)
      println("W has ",length(p)," with deg=",r.deg," eig=",r.eig,
            " when ",length(m)," occur in given list")
      p=abs.(positionssgn(ud,r.deg))
      println("they are in positions ",p," with eig ",eigs[p])
      println("this error is for ",findfirst(==(r),il),"-th in given list")
      error()
    else
      if length(m)>1
        if !haskey(s,:ambig) s.ambig=[] end
        push!(s.ambig,[m,p])
      end
      res[m]=p
    end
    pos=setdiff(pos,m)
  end
  s.charNumbers=res
end

"""
getunpdeg(W;ennola=true,principal=Int[]) 
if given principal=list of indices of 1-series to use, 
  ie indices in uc.harishChandra
Find (unipotent degrees, eigenvalues) for W
see "split spetses" chapter 5
"""
function getunpdeg(W;ennola=true,principal=Int[])
  function ser(d,H)
    s=Series(spets(W),
    subspets(spets(W),Int[],classreps(W)[position_regular_class(W,d)]),1,d;
     NC=true)
    relative_group(s)
    s.Hecke=H
    e=Tests.EigenAndDegHecke(s)
    s.eigen=map(x->x.eig,e)
    s.frac=map(x->x.frac,e)
    finddeinW(W,s,e)
    s
  end
  chars=[]
  uc=UnipotentCharacters(W)
  function process(s)
    nbu=length(uc)
    n=setdiff(s.charNumbers,chars)
    d=setdiff(n,chars)
    if !isempty(d)
      println("By ",s.d,"-series found ",asranges(d))
    end
    p=findfirst(ss->ss.d==s.d && ss.cuspidal==s.cuspidal,sers)
    if isnothing(p) push!(sers,s)
    elseif sers[p].Hecke.para!=s.Hecke.para error()
    end
    chars=union(chars,n)
    if length(chars)==nbu
      if !isempty(d) println("****** Found all!!!") end
      return true
    elseif length(d)>0 println("  --- still to find: ",
      asranges(setdiff(1:nbu,chars)))
      return false
    end
  end
  z=gcd(degrees(refltype(W)[1]))
  H=hecke(W,Mvp(:q))
  sers=[]
  if ennola
    for k in 0:z-1 process(ser(E(z,k),ennola_twist(H,E(z,k),E(1)))) end
  else process(ser(E(1),H))
  end
  for k in principal
    h=uc.harishChandra[k]
    if !isempty(h.relativeType) && !isempty(h.levi)
      H=subspets(spets(W),h.levi)
      s=Series(spets(W),H,FindCuspidalInLevi(h.cuspidalName,H),0;NC=true)
      relative_group(s)
      s.Hecke=UnipotentCharactersOps.RelativeHecke(uc,k,Mvp(:q))
      e=Tests.EigenAndDegHecke(s)
      s.eigen=map(x->x.eig,e)
      s.frac=map(x->x.frac,e)
      finddeinW(W,s,e)
      println("Using ",s)
      process(s)
    end
  end
  cyclic=filter(x->length(relative_degrees(W,x))==1,regular_eigenvalues(W))
# todo=vcat(map(d->map(k->E(order(d),k),prime_residues(order(d))),reverse(cyclic))...)
  todo=reverse(cyclic)
  for d in todo
    H=FindRelativeHecke(W,d,chars)
    for k in 0:z-1
      process(ser(d*E(z,k),ennola_twist(H,E(z,k),d)))
      todo=setdiff(todo,[d*E(z,k)])
    end
  end
  todo=reverse(setdiff(regular_eigenvalues(W),union(cyclic,E.(z,1:z))))
  if !isempty(todo)
  todo=union(map(d->map(k->E(order(d),k),prime_residues(order(d))),todo)...)
  end
  for d in todo
    H=FindRelativeHecke(W,d,chars)
    if H==false Print("**** Could not compute ",d,"-series!!!!\n")
    else process(ser(d,H));todo=setdiff(todo,[d]) end
  end
  [sers,asranges(sort(chars))]
end

## find pairs of families with same a and A
function dbfams(W)
  uc=UnipotentCharacters(W)
  aA=map(f->(uc.a[f.charNumbers[1]],uc.A[f.charNumbers[1]]),uc.families)
  db=map(filter(x->x[2]>1,tally(aA)))do x
    filter(i->aA[i]==x[1],eachindex(aA))
  end
  for l in db
    for f in l
     print("family ",f,"#",length(uc.families[f])," a=",aA[f][1]," A=",aA[f][2],"   ")
    end
    println()
  end
end

## given a Hecke algebra [where parameters are monomials in q], 
## find the minimal field containing coeffs of prod(x-p_i) where p_i 
## are the parameters normalized so that highest degree has coeff 1.
function Nf.NF(H::HeckeAlgebra)
  l=map(x->x//first(coefficients(x[1])),H.para)
  l=unique(map(x->prod(Mvp(:x).-x),l))
# xprintln(improve_type(l))
  unique(map(x->Nf.NF(collect(coefficients(x))),l))
end

function testg(arg)
  l=getunpdeg(crg(arg...))[1]
# l=Filtered(l,x->x<>false)
  println("-----------------------------------------")
  xdisplay(map(x->(x,NF(hecke(x))),l))
end

## s is a series returned by getunpdeg; show qEigen
#CheckFrac:=function(W)local ct,r,s,uc,i,v;
#  r:=getunpdeg(W)[1];uc:=UnipotentCharacters(W);
#  for s in r do
#    ct:=CharTable(s.Hecke).irreducibles*Mvp("q")^0;
#    ct:=List(ct,x->Set(Flat(List(x,y->List(y.elm,z->z.coeff)))));
#    ct:=List(ct,x->Union(Set(List(x,Mod1)),[0]));
#    ct:=List(ct,x->Mod1(1/Lcm(List(x,Denominator))));
#    s.irch:=ct;
#  od;
#  s:=[];
#  for i in [1..Size(uc)] do
#    v:=Filtered(r,s->i in s.charNumbers);
#    v:=List(v,function(s)local p;
#      p:=Position(s.charNumbers,i);
#      return [s.d,s.irch[p],s.frac[p]];end);
#    if Length(Set(List(v,x->x[2])))>1 or Length(Set(List(v,x->x[3])))>1 then
#      Error();fi;
#    if Length(v)>0 then
#    Add(s,rec(no:=i,sers:=Set(List(v,x->x[1])),irch:=v[1][2],frac:=v[1][3],
#      fam:=PositionProperty(uc.families,x->i in x.charNumbers)));
#    fi;
#  od;
#  if ForAny(s,function(x)local h;
#    h:=First(uc.harishChandra,h->x.no in h.charNumbers);
#    if x.frac<>h.qEigen then Error("for ",x.no," computed qEig=",
#      x.frac," while stored qEig=",h.qEigen,"\n");fi;
#    return x.frac<>h.qEigen;end) then Error("oopos");fi;
#  s:=Filtered(s,x->x.irch<>0 or x.frac<>0);
#  Print(FormatTable(List(s,x->[x.irch,x.frac,x.fam,x.sers]),
#     rec(rowLabels:=List(s,x->x.no), rowsLabel:="no",
#    columnLabels:=["irr.ch.","frac.eig.","family","series"])),"\n");
#end;

#findzeta:=function(W)local ct,r,s,uc,i,v;
#  r:=getunpdeg(W)[1];uc:=UnipotentCharacters(W);
#  s:=[];
#  for i in [1..Size(uc)] do
#    v:=Filtered(r,s->i in s.charNumbers);
#    s[i]:=List(v,s->Root1(r=s.d)^(uc.a[i]+uc.A[i]
#         -Sum(ReflectionDegrees(W)+ReflectionCoDegrees(W))));
#  od;
#  return s;
#end;

#checkseries:=function(s)local e,n,W;
#  W:=s.spets;
#  if IsBound(s.eigen) then
#    e:=Eigenvalues(UnipotentCharacters(W)){s.charNumbers};
#    if s.eigen<>e then
#      n:=CharNames(UnipotentCharacters(W)){s.charNumbers};
#      ChevieErr(s," actual eigen differ from predicted eigen");
#      CHEVIE.Check.EqTables(
#                 rec(rowLabels:=["actual"],columnLabels:=n,scalar:=[e]),
#                 rec(rowLabels:=["pres"],columnLabels:=n,scalar:=[s.eigen]));
#    fi;
#  else Print("no eigen");
#  fi;
#end;
