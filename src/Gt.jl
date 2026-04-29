module Gt
using ..Chevie
export RationalUnipotentClasses, closed_subsystems, closed_subsystemsPoset,
 closed_subsystems_reps, ClassTypes

struct RationalUnipotentClasses
  WF
  p::Int
  l::Vector{@NamedTuple{card::CycPol{Cyc{Rational{Int}}}, 
                        class::UnipotentClass, 
                        classno::Int, AuNo::Int}}
end

function RationalUnipotentClasses(WF,p::Int=0)
  u=UnipotentClasses(WF, p)
  t=XTable(u;classes=true)
  RationalUnipotentClasses(WF,p,
    map(i->(card=CycPol(t.cardClass[i]), class=u.classes[t.classes[i][1]], 
    classno=t.classes[i][1], AuNo=t.classes[i][2]), eachindex(t.classes)))
end

Base.show(io::IO,::MIME"text/latex",c::RationalUnipotentClasses)=
  print(io,TeXs(c))

Base.show(io::IO,cl::RationalUnipotentClasses)=
  print(io,"RationalUnipotentClasses(",cl.WF,",",cl.p,")")

function Base.show(io::IO,::MIME"text/plain",cl::RationalUnipotentClasses)
  printTeX(io,"Rational unipotent classes of "*TeXs(cl.WF)*" in char.",cl.p,"\n")
  ow=CycPol(generic_order(cl.WF))
  cards=map(x->xrepr(io,ow/x.card),cl.l)
  cl=map(x->"{"*fromTeX(io,name(x.class,TeX=true)*"}_{"*classnames(x.class.Au)[x.AuNo]*"}"),cl.l)
  showtable(io,reshape(cards,:,1);rows_label="classes",row_labels=cl,
                                                  col_labels=["centralizer"])
end

"""
`sumr(W)` where `W` is a Weyl group

if `N=nref(W)` build the `2N×2N` array `W.sumr` such that `sumr[i,j]==k` if
`root(W,i)+root(W,j)=root(W,k)`  and `==0 if `root(W,i)+root(W,j)` is not a
root.
"""
function sumr(W)
  get!(W,:sumr)do
    function possum(i,j)
      s=roots(W,i)+roots(W,j)
      p=findfirst(==(s),roots(W))
      isnothing(p) ? 0 : p
    end
    [possum(i,j) for i in 1:2nref(W),  j in 1:2nref(W)]
  end::Matrix{Int}
end

# closure of the subset roots(W,l)
function closure(W,l)
  psum=sumr(W)
  lset=Set(l)
  res=Int[]
  l=copy(l)
  for r in l
    for s in res
      u=psum[r,s]
      if !iszero(u) && !(u in lset) 
        push!(l,u) 
        push!(lset,u)
      end
    end
    push!(res,r)
  end
  sort!(l)
end

"""
`closed_subsystemsPoset(W)` 

`W`  should  be  a  Weyl  group.  The  function returns the Poset of closed
subsystems  of the root system of `W`. Each closed subsystem is represented
by  the list of indices of its simple roots.  If `W` is the Weyl group of a
reductive  group  `𝐆  `,  then  closed  subsystem  correspond  to reductive
subgroups of maximal rank. And all such groups are obtained this way, apart
from  some exceptions  in characteristics  2 and  3, see [mt11; Proposition
13.4](@cite).

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> closed_subsystemsPoset(W)
1 2<1 4<4<∅
1 2<1 5<1<∅
1 2<2 6<6<∅
1 2<3 5<5<∅
1 4<1
1 5<6
1 5<5
2 6<2<∅
3 5<3<∅
```
"""
function closed_subsystemsPoset(W)
  get!(W, :closedsubsets)do
  l = [Int[]]
  ldict=Dict{Vector{Int},Int}()
  new = [1]
  covers = Tuple{Int,Int}[]
  for w in new
    for f in setdiff(1:nref(W), l[w])
      n=closure(W,vcat(l[w],[f,f+nref(W)]))
      if haskey(ldict,n)
        push!(covers,(ldict[n], w))
      else
        push!(l, n)
        ldict[n]=length(l)
        push!(covers,(length(l),w))
        push!(new, length(l))
      end
    end
  end
  P=Poset(CPoset(covers),l)
  P.show_element=function(io,x,n)
    e=Weyl.simpleroots_subsystem(W,x.elements[n])
    print(io,isempty(e) ? "∅" : join(e," "))
  end
  P
  end
end

function closed_subsystemsPoset(W::CoxeterCoset)
  get!(W, :closedsubsets)do
    P=closed_subsystems(Group(W))
    P=induced(P,filter(x->onsets(x,W.phi)==x,P.elements))
    P.show_element=function(io,x,n)
      e=Weyl.simpleroots_subsystem(Group(W),x.elements[n])
      print(io,isempty(e) ? "∅" : join(e," "))
    end
    P
  end
end

"""
`closed_subsystems(W)` 

`W`  should be a Weyl group. The function returns the list, sorted by rank,
of  closed subsystems of the  root system of `W`.  Each closed subsystem is
represented  by the list of indices of its simple roots. If `W` is the Weyl
group  of  a  reductive  group  `𝐆  `,  then closed subsystem correspond to
reductive  subgroups of maximal rank. And all such groups are obtained this
way,  apart from  some exceptions  in characteristics  2 and  3, see [mt11;
Proposition 13.4](@cite).

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> closed_subsystems(W)
12-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [3]
 [4]
 [5]
 [6]
 [1, 2]
 [1, 4]
 [1, 5]
 [2, 6]
 [3, 5]
```
"""
function closed_subsystems(W;verbose=false) # code from Meinolf Geck
  mat=sumr(W)
  systems=[Int[]]
  prev=map(i->[i],1:nref(W))
  append!(systems,prev)
  if verbose print("#I 1 ",length(prev)," ") end
  for _ in 1:ngens(W)-1
    new=Vector{Int}[]
    for s in prev, r in s[end]+1:nref(W)
      if all(i->iszero(mat[r,i+nref(W)]),s) push!(new,vcat(s,[r])) end
    end
    append!(systems,new)
    if verbose print(length(new)," ") end
    prev=new
  end
  if verbose println("total ",length(systems)) end
  systems
end

"""
`closed_subsystems_reps(W)` 

`W`  should be  a Weyl  group. The  function returns representatives of the
`W`-orbits  of closed subsystems of the root system of `W`, sorted by rank.
Each  closed subsystem is represented by the  list of indices of its simple
roots.  If `W`  is the  Weyl group  of a  reductive group `𝐆 `, then closed
subsystem  correspond to reductive subgroups of  maximal rank. And all such
groups are obtained this way, apart from some exceptions in characteristics
2 and 3, see [mt11; Proposition 13.4](@cite).

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> map(I->reflection_subgroup(W,I),closed_subsystems_reps(W))
6-element Vector{FiniteCoxeterSubGroup{Perm{Int16},Int64}}:
 G₂₍₎=Φ₁²
 G₂₍₁₎=A₁Φ₁
 G₂₍₂₎=Ã₁Φ₁
 G₂
 G₂₍₁₄₎=A₁×Ã₁
 G₂₍₁₅₎=A₂
```
"""
function closed_subsystems_reps(W;verbose=false)# code from Meinolf Geck
  function subsorbits(W,sys::Vector{Vector{Int}}) # W-orbits representatives on subsystems list sys
    rest=eachindex(sys)
    orb=Int[]
    while !isempty(rest)
      push!(orb,rest[1])
      r=filter(i->transporting_elt(W,sys[rest[1]],sys[i],onsets)!==nothing,rest)
      rest=setdiff(rest,r)
    end
    sys[orb]::Vector{Vector{Int}}
  end
  mat=sumr(W)
  systems=[Int[]]
  prev=map(x->[x],unique(simple_reps(W,1:ngens(W))))
  append!(systems,prev)
  if verbose print("#I 1 ",length(prev)," ") end
  for _ in 1:ngens(W)-1
    new=Vector{Int}[]
    for s in prev, r in s[end]+1:nref(W)
      if all(i->iszero(mat[r,i+nref(W)]),s) push!(new,vcat(s,[r])) end
    end
    new=subsorbits(W,new)
    append!(systems,new)
    if verbose print(length(new)," ") end
    prev=new
  end
  if verbose println("total ",length(systems)) end
  systems
end

@GapObj struct ClassType
  CGs
  cent::CycPol
  unip
end

@GapObj mutable struct ClassTypes
  p::Int
  WF
  ss::Vector{ClassType}
end

"""
`ClassTypes(G[,p])`

`G`  should be a root  datum or a twisted  root datum representing a finite
reductive  group ``𝐆 ^F`` and  `p` should be a  prime. The function returns
the class types of `G` in characteristic `p` (in good characteristic if `p`
is  omitted). Two elements  of ``𝐆 ^F``  have the same  class type if their
centralizers  are  conjugate.  If  `su`  is  the Jordan decomposition of an
element  `x`, the class type of `x` is  determined by the class type of its
semisimple part `s` and the unipotent class of `u` in ``C_𝐆 (s)``.

The   function  `ClassTypes`  is  presently  only  implemented  for  simply
connected  groups, where  ``C_𝐆 (s)``  is connected.  This section is a bit
experimental and may change in the future.

`ClassTypes`  returns a  `struct` which  contains a  list of classtypes for
semisimple  elements,  which  are  represented  by  `subspets`  and contain
additionnaly information on the unipotent classes of ``C_𝐆 (s)``.

Let us give some examples:

```julia-repl
julia> t=ClassTypes(rootdatum(:sl,3))
ClassTypes(A₂,good characteristic)
┌──────────┬─────────┐
│C_G(s)    │ |C_G(s)|│
├──────────┼─────────┤
│A₂₍₎=Φ₁²  │      Φ₁²│
│A₂₍₎=Φ₁Φ₂ │     Φ₁Φ₂│
│A₂₍₎=Φ₃   │       Φ₃│
│A₂₍₁₎=A₁Φ₁│   qΦ₁²Φ₂│
│A₂        │q³Φ₁²Φ₂Φ₃│
└──────────┴─────────┘
```
By   default,  only  information  about  semisimple  centralizer  types  is
returned:   the type, and its generic order.

```julia-rep1
julia> xdisplay(t;unip=true)
ClassTypes(A₂,good characteristic)
┌──────────┬───────────────┐
│C_G(s)    │    u |C_G(su)|│
├──────────┼───────────────┤
│A₂₍₎=Φ₁²  │    1       Φ₁²│
│A₂₍₎=Φ₁Φ₂ │    1      Φ₁Φ₂│
│A₂₍₎=Φ₃   │    1        Φ₃│
│A₂₍₁₎=A₁Φ₁│   11    qΦ₁²Φ₂│
│          │    2       qΦ₁│
│A₂        │  111 q³Φ₁²Φ₂Φ₃│
│          │   21      q³Φ₁│
│          │    3       3q²│
│          │ 3_ζ₃       3q²│
│          │3_ζ₃²       3q²│
└──────────┴───────────────┘
```
Here  we  have  displayed  information  on  unipotent  classes,  with their
centralizer.

```julia-rep1
julia> xdisplay(t;nClasses=true)
ClassTypes(A₂,good characteristic)
┌──────────┬─────────────────────────┐
│C_G(s)    │       nClasses  |C_G(s)|│
├──────────┼─────────────────────────┤
│A₂₍₎=Φ₁²  │(q²-5q+2q₃+4)/6       Φ₁²│
│A₂₍₎=Φ₁Φ₂ │       (q²-q)/2      Φ₁Φ₂│
│A₂₍₎=Φ₃   │  (q²+q-q₃+1)/3        Φ₃│
│A₂₍₁₎=A₁Φ₁│       (q-q₃-1)    qΦ₁²Φ₂│
│A₂        │             q₃ q³Φ₁²Φ₂Φ₃│
└──────────┴─────────────────────────┘
```
Here  we have added information on how many semisimple conjugacy classes of
`𝐆  ^F` have a given type. The  answer in general involves variables of the
form `qₐ` which represent `gcd(q-1,a)`.

Finally an example in bad characteristic:

```julia-rep1
julia> t=ClassTypes(coxgroup(:G,2),2);xdisplay(t;nClasses=true)
ClassTypes(G₂,char. 2)
ClassTypes(G₂,char. 2)
┌──────────┬──────────────────────────────┐
│C_G(s)    │         nClasses     |C_G(s)|│
├──────────┼──────────────────────────────┤
│G₂₍₎=Φ₁²  │(q²-8q+2q₃+10)/12          Φ₁²│
│G₂₍₎=Φ₁Φ₂ │        (q²-2q)/4         Φ₁Φ₂│
│G₂₍₎=Φ₁Φ₂ │        (q²-2q)/4         Φ₁Φ₂│
│G₂₍₎=Φ₆   │    (q²-q-q₃+1)/6           Φ₆│
│G₂₍₎=Φ₃   │    (q²+q-q₃+1)/6           Φ₃│
│G₂₍₎=Φ₂²  │ (q²-4q+2q₃-2)/12          Φ₂²│
│G₂₍₁₎=A₁Φ₁│       (q-q₃-1)/2       qΦ₁²Φ₂│
│G₂₍₁₎=A₁Φ₂│       (q-q₃+1)/2       qΦ₁Φ₂²│
│G₂₍₂₎=Ã₁Φ₁│          (q-2)/2       qΦ₁²Φ₂│
│G₂₍₂₎=Ã₁Φ₂│              q/2       qΦ₁Φ₂²│
│G₂        │                1 q⁶Φ₁²Φ₂²Φ₃Φ₆│
│G₂₍₁₅₎=A₂ │         (q₃-1)/2    q³Φ₁²Φ₂Φ₃│
│G₂₍₁₅₎=²A₂│         (q₃-1)/2    q³Φ₁Φ₂²Φ₆│
└──────────┴──────────────────────────────┘
```
We  notice that if `q` is  a power of `2` such  that `q≡2 (mod 3)`, so that
`q₃=1`,  some class types do not exist. We can see what happens by giving a
specific value to `q₃`:

```julia-rep1
julia> xdisplay(t(;q_3=1);nClasses=true)
ClassTypes(G₂,char. 2) q₃=1
┌──────────┬──────────────────────────┐
│C_G(s)    │     nClasses     |C_G(s)|│
├──────────┼──────────────────────────┤
│G₂₍₎=Φ₁²  │(q²-8q+12)/12          Φ₁²│
│G₂₍₎=Φ₁Φ₂ │    (q²-2q)/4         Φ₁Φ₂│
│G₂₍₎=Φ₁Φ₂ │    (q²-2q)/4         Φ₁Φ₂│
│G₂₍₎=Φ₆   │     (q²-q)/6           Φ₆│
│G₂₍₎=Φ₃   │     (q²+q)/6           Φ₃│
│G₂₍₎=Φ₂²  │   (q²-4q)/12          Φ₂²│
│G₂₍₁₎=A₁Φ₁│      (q-2)/2       qΦ₁²Φ₂│
│G₂₍₁₎=A₁Φ₂│          q/2       qΦ₁Φ₂²│
│G₂₍₂₎=Ã₁Φ₁│      (q-2)/2       qΦ₁²Φ₂│
│G₂₍₂₎=Ã₁Φ₂│          q/2       qΦ₁Φ₂²│
│G₂        │            1 q⁶Φ₁²Φ₂²Φ₃Φ₆│
└──────────┴──────────────────────────┘
```
"""
function ClassTypes(W,p=0)
  if W isa Spets WF=W;W=Group(W)
  else WF=spets(W)
  end
  l=vcat(twistings.(Ref(WF),sscentralizer_reps(Group(WF), p))...)
  ClassTypes(p,WF,map(x->ClassType(x,CycPol(generic_order(x,Pol(:q))),
    RationalUnipotentClasses(x,p).l,Dict{Symbol,Any}()),l),Dict{Symbol,Any}())
end

bracket_if_needed(v)=occursin(r"[-+*/]",v[nextind(v,0,2):end]) ? "("*v*")" : v

Base.show(io::IO, ::MIME"text/latex", r::ClassTypes)=print(io,TeXs(r))
 
function Base.show(io::IO,r::ClassTypes)
  res=string("ClassTypes(\$",TeX(io,r.WF),"\$")
  if r.p==0 res*=",good characteristic)"
  else res*=string(",char. ",r.p,")")
  end
  if haskey(r,:specialized)
    res*=string(" \$",join(map(x->string(x[1],"=",x[2]),collect(r.specialized))," "),"\$")
  end
  printTeX(io,res*"\n")
  function nc(p)
    p=Mvp(p)
    d=lcm(denominator.(values(p.d)))
    p=bracket_if_needed(xrepr(io,improve_type(p*d)))
    d==1 ? p : string(p,"/",d)
  end
  classes=get(io,:nClasses,false)
  col_labels=String[]
  if classes push!(col_labels, "nClasses") end
  if get(io,:unip,false)
    row_labels=[]
    push!(col_labels,"u")
    push!(col_labels,"|C_G(su)|")
    t=[]
    for x in r.ss
      u=RationalUnipotentClasses(x.CGs, r.p).l
      for c in u
        v=String[]
        if isone(c.card)
          if classes push!(v, nc(nconjugacy_classes(x,r.WF,r.p))) end
          push!(row_labels,TeX(io,x.CGs))
        else
          push!(row_labels," ")
          if classes push!(v,"") end
        end
        push!(v,Ucl.nameclass(c.class;io.dict...,class=c.AuNo))
        push!(v,xrepr(io,x.cent//c.card))
        push!(t,v)
      end
    end
    t=toM(t)
  else
    row_labels=map(x->TeX(io,x.CGs),r.ss)
    t=[]
    if classes push!(t,nc.(nconjugacy_classes(r))) end
    push!(t,map(x->xrepr(io,x.cent),r.ss))
    push!(col_labels,"|C_G(s)|")
    t=permutedims(toM(t))
  end
  showtable(io,t;col_labels,row_labels,rows_label="C_G(s)")
end

function (C::ClassTypes)(;arg...)
  C=copyGapObj(C)
  C.specialized=arg
  C.ss=map(function (x)
           res=copyGapObj(x)
           res.nClasses=nconjugacy_classes(x,C.WF,C.p)(;arg...)
           res
          end, C.ss)
  filter!(r->!iszero(nconjugacy_classes(r,C.WF,C.p)),C.ss)
  C
end

function nconjugacy_classes(C::ClassTypes)
  W=Group(C.WF)
  if length(fundamental_group(W))>1
    print("# Nr classes each type implemented only for simply connected groups")
    return
  end
  filter!(r->nconjugacy_classes(r,C.WF,C.p)!=0,C.ss)
  nconjugacy_classes.(C.ss,Ref(C.WF),C.p)
end

# See Fleischmann-Janiszczak AAECC 1996 definition 2.1
function nconjugacy_classes(r::ClassType,WF,p)
  get!(r,:nClasses)do
    HF=r.CGs
    H=Group(HF)
    W=Group(WF)
    P=closed_subsystemsPoset(W)
    elts=copy(P.elements)
    elts=filter(x->onsets(x,HF.phi)==x && issubset(inclusiongens(H),x),elts)
    P=induced(P,elts)
# here P poset of HF.phi-stable closed subsets containing roots(H)
#   InfoChevie("# ", HF, "==>", P, "")
    l=map(x->spets(reflection_subgroup(W, x), HF.phi), P.elements)
    b=collect(keys(factor(PermRoot.BadNumber(W))))
    l=map(l)do RF
      R=Group(RF)
      res=reduce(*,map(p->Mvp(:q)-p[2],filter(y->y[1]==1,degrees(RF)));init=Mvp(1))
      if semisimplerank(R)==0 return res end
      d=filter(x->!(x in [0,1,p]),smith(simpleroots(R)))
      if isempty(d) return res end
      res*prod(d)do i
                   if i==2 && (2 in b && p!=2) return 2 end
                   if mod(p,i)==1 return i end
                   Mvp(Symbol("q_",i))
                 end
    end
    mu=moebius(P)
    n=stabilizer(W,sort(inclusiongens(H)),onsets)
    n=sum(mu.*l) // length(centralizer(n, HF.phi))
#   InfoChevie("==>", n, "\n")
    n
  end
end
end
