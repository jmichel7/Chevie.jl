module Gt
using ..Gapjm
export RationalUnipotentClasses, ClosedSubsets, ClassTypes

function RationalUnipotentClasses(WF, p)
  u=UnipotentClasses(WF, p)
  t=Ucl.XTable(u;classes=true)
  map(i->Dict{Symbol, Any}(:card => CycPol(t.cardClass[i]), 
                           :class => u.classes[t.classes[i][1]], 
                           :classno => t.classes[i][1], 
                           :AuNo => t.classes[i][2]), 
           eachindex(t.classes))
end

# returns the Poset of closed subsystems of the root system of W
function ClosedSubsets(W)
  gets(W, :closedsubsets)do
  function possum(i,j)
    p=findfirst(==(roots(W,i)+roots(W,j)),roots(W))
    isnothing(p) ? 0 : p
  end
  psum=[possum(i,j) for i in 1:2nref(W),  j in 1:2nref(W)]
  closure=function(l,new)
    nnew = new
    while true
      new = nnew
      nnew = Int[]
      for i in new, j in l
        if psum[i,j]!=0 push!(nnew, psum[i,j]) end
      end
      l = union(l, new)
      nnew = setdiff(nnew, l)
      if isempty(nnew) break end
    end
    return sort(l)
  end
  l = [Int[]]
  new = [1]
  covered = [Int[]]
  for w in new
    for f in setdiff(1:nref(W), l[w])
      n = closure(l[w], [f, f + nref(W)])
      p = findfirst(==(n),l)
      if isnothing(p)
          push!(l, n)
          push!(covered, [w])
          push!(new, length(l))
      else push!(covered[p], w)
      end
    end
  end
  covered=unique.(covered)
  P=Poset(incidence(Poset(covered)))
  P.prop[:elements]=l
  P.prop[:labels]=map(x->join(x," "),l)
  P
  end
end

struct ClassType
  CGs
  cent::CycPol
  unip
  prop::Dict{Symbol,Any}
end

struct ClassTypes
  p::Int
  WF
  ss::Vector{ClassType}
  prop::Dict{Symbol,Any}
end

"""
`ClassTypes(G[,p])`

`G`  should be a root  datum or a twisted  root datum representing a finite
reductive  group `𝐆 ^F` and `p` should be a prime. The function returns the
class  types of `G` in characteristic `p` (in good characteristic if `p` is
omitted).  Two  elements  of  `𝐆  ^F`  have  the  same  class type if their
centralizers  are  conjugate.  If  `su`  is  the Jordan decomposition of an
element  `x`, the class type of `x` is  determined by the class type of its
semisimple part `s` and the unipotent class of `u` in `C_𝐆 (s)`.

The   function  `ClassTypes`  is  presently  only  implemented  for  simply
connected  groups, where  `C_𝐆 (s)`  is  connected. This  section is  a bit
experimental and may change in the future.

`ClassTypes`  returns a  `struct` which  contains a  list of classtypes for
semisimple  elements,  which  are  represented  by  `subspets`  and contain
additionnaly information on the unipotent classes of `C_𝐆 (s)`.

Let us give some examples:

```julia-repl
julia> t=ClassTypes(rootdatum(:sl,3))
ClassTypes(A₂,good characteristic)
    C_G(s)│ |C_G(s)|
──────────┼──────────
A₂₍₎=Φ₁²  │      Φ₁²
A₂₍₎=Φ₁Φ₂ │     Φ₁Φ₂
A₂₍₎=Φ₃   │       Φ₃
A₂₍₁₎=A₁Φ₁│   qΦ₁²Φ₂
A₂        │q³Φ₁²Φ₂Φ₃
```
By   default,  only  information  about  semisimple  centralizer  types  is
returned:   the type, and its generic order.

```julia-rep1
julia> xprint(t;unip=true)
ClassTypes(A₂,good characteristic)
    C_G(s)│      u |C_G(su)|
──────────┼──────────────────
A₂₍₎=.Φ₁² │              Φ₁²
A₂₍₎=.Φ₁Φ₂│             Φ₁Φ₂
A₂₍₎=.Φ₃  │               Φ₃
A₂₍₁₎=A₁Φ₁│     11    qΦ₁²Φ₂
          │      2       qΦ₁
A₂        │    111 q³Φ₁²Φ₂Φ₃
          │     21      q³Φ₁
          │      3       3q²
          │ 3_{ζ₃}       3q²
          │3_{ζ₃²}       3q²
```
Here  we  have  displayed  information  on  unipotent  classes,  with their
centralizer.

```julia-rep1
julia> xprint(t;nClasses=true)
ClassTypes(A₂,good characteristic)
    C_G(s)│       nClasses  |C_G(s)|
──────────┼──────────────────────────
A₂₍₎=.Φ₁² │(q²-5q+2q₃+4)/6       Φ₁²
A₂₍₎=.Φ₁Φ₂│       (q²-q)/2      Φ₁Φ₂
A₂₍₎=.Φ₃  │  (q²+q-q₃+1)/3        Φ₃
A₂₍₁₎=A₁Φ₁│       (q-q₃-1)    qΦ₁²Φ₂
A₂        │             q₃ q³Φ₁²Φ₂Φ₃
```
Here  we have added information on how many semisimple conjugacy classes of
`𝐆  ^F` have a given type. The  answer in general involves variables of the
form `qₐ` which represent `gcd(q-1,a)`.

Finally an example in bad characteristic:

```julia-rep1
julia> t=ClassTypes(coxgroup(:G,2),2);xprint(t;nClasses=true)
ClassTypes(G₂,char. 2)
    C_G(s)│         nClasses     |C_G(s)|
──────────┼───────────────────────────────
G₂₍₎=.Φ₁² │(q²-8q+2q₃+10)/12          Φ₁²
G₂₍₎=.Φ₁Φ₂│        (q²-2q)/4         Φ₁Φ₂
G₂₍₎=.Φ₁Φ₂│        (q²-2q)/4         Φ₁Φ₂
G₂₍₎=.Φ₆  │    (q²-q-q₃+1)/6           Φ₆
G₂₍₎=.Φ₃  │    (q²+q-q₃+1)/6           Φ₃
G₂₍₎=.Φ₂² │ (q²-4q+2q₃-2)/12          Φ₂²
G₂₍₁₎=A₁Φ₂│       (q-q₃+1)/2       qΦ₁Φ₂²
G₂₍₁₎=A₁Φ₁│       (q-q₃-1)/2       qΦ₁²Φ₂
G₂₍₂₎=Ã₁Φ₂│              q/2       qΦ₁Φ₂²
G₂₍₂₎=Ã₁Φ₁│          (q-2)/2       qΦ₁²Φ₂
G₂        │                1 q⁶Φ₁²Φ₂²Φ₃Φ₆
G₂₍₁₅₎=²A₂│         (q₃-1)/2    q³Φ₁Φ₂²Φ₆
G₂₍₁₅₎=A₂ │         (q₃-1)/2    q³Φ₁²Φ₂Φ₃
```
We  notice that if `q` is  a power of `2` such  that `q≡2 (mod 3)`, so that
`q₃=1`,  some class types do not exist. We can see what happens by giving a
specific value to `q₃`:

```julia-rep1
julia> xprint(t(;q_3=1);nClasses=true)
ClassTypes(G₂,char. 2) q_3=1
    C_G(s)│     nClasses     |C_G(s)|
──────────┼───────────────────────────
G₂₍₎=.Φ₁² │(q²-8q+12)/12          Φ₁²
G₂₍₎=.Φ₁Φ₂│    (q²-2q)/4         Φ₁Φ₂
G₂₍₎=.Φ₁Φ₂│    (q²-2q)/4         Φ₁Φ₂
G₂₍₎=.Φ₆  │     (q²-q)/6           Φ₆
G₂₍₎=.Φ₃  │     (q²+q)/6           Φ₃
G₂₍₎=.Φ₂² │   (q²-4q)/12          Φ₂²
G₂₍₁₎=A₁Φ₂│          q/2       qΦ₁Φ₂²
G₂₍₁₎=A₁Φ₁│      (q-2)/2       qΦ₁²Φ₂
G₂₍₂₎=Ã₁Φ₂│          q/2       qΦ₁Φ₂²
G₂₍₂₎=Ã₁Φ₁│      (q-2)/2       qΦ₁²Φ₂
G₂        │            1 q⁶Φ₁²Φ₂²Φ₃Φ₆
```
"""
function ClassTypes(W,p=0)
  if W isa Spets WF=W;W=Group(W)
  else WF=spets(W)
  end
  l=vcat(twistings.(Ref(WF),SScentralizer_representatives(Group(WF), p))...)
  ClassTypes(p,WF,map(x->ClassType(x,CycPol(generic_order(x,Pol(:q))),
    RationalUnipotentClasses(x,p),Dict{Symbol,Any}()),l),Dict{Symbol,Any}())
end

bracket_if_needed(v)=occursin(r"[-+*/]",v[nextind(v,0,2):end]) ? "("*v*")" : v

function Base.show(io::IO,r::ClassTypes)
  res=string("ClassTypes(",sprint(show,r.WF;context=io))
  if r.p==0 res*=",good characteristic)"
  else res*=string(",char. ",r.p,")")
  end
  if haskey(r.prop,:specialized)
    res*=string(" ",join(map(x->string(x[1],"=",x[2]),
                             collect(r.prop[:specialized]))," "))
  end
  println(io,res)
  function nc(p)
    local d
    p=Mvp(p)
    d=lcm(map(denominator,values(p.d)))
    p=bracket_if_needed(sprint(show,improve_type(p*d);context=io))
    d==1 ? p : string(p,"/",d)
  end
  classes=get(io,:nClasses,false)
  columnLabels=String[]
  if classes push!(columnLabels, "nClasses") end
  if get(io,:unip,false)
    rowLabels=[]
    push!(columnLabels,"u")
    push!(columnLabels,fromTeX(io,"|C_G(su)|"))
    t=[]
    for x in r.ss
      u=RationalUnipotentClasses(x.CGs, r.p)
      for c in u
        v=String[]
        if isone(c[:card])
          if classes push!(v, nc(nconjugacy_classes(x,r.WF,r.p))) end
          push!(rowLabels,sprint(show,x.CGs;context=io))
        else
          push!(rowLabels," ")
          if classes push!(v,"") end
        end
        push!(v,Ucl.nameclass(merge(c[:class].prop,Dict(:name=>c[:class].name)),
                              merge(io.dict,Dict(:class=>c[:AuNo]))))
        push!(v,sprint(show,x.cent//c[:card],context=io))
        push!(t,v)
      end
    end
    t=toM(t)
  else
    rowLabels=map(x->sprint(show,x.CGs;context=io),r.ss)
    t=[]
    if classes push!(t,nc.(nconjugacy_classes(r))) end
    push!(t,map(x->sprint(show,x.cent;context=io),r.ss))
    push!(columnLabels,fromTeX(io,"|C_G(s)|"))
    t=permutedims(toM(t))
  end
  format(io,t;col_labels=columnLabels,row_labels=rowLabels,
         rows_label=fromTeX(io,"C_G(s)"))
end

function (C::ClassTypes)(;arg...)
  C=deepcopy(C)
  C.prop[:specialized]=arg
  C.ss.=map(function (x)
           res=deepcopy(x)
           res.prop[:nClasses]=nconjugacy_classes(x,C.WF,C.p)(;arg...)
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

# the Moebius function on the intervals x:maximum(P)
function moebius(P::Poset)
  o=linear_extension(P)
  mu=zeros(Int,length(P))
  mu[o[length(P)]]=1
  less=i->filter(j->j!=i && incidence(P)[i,j], 1:length(P))
  for i in o[length(P)-1:-1:1] mu[i]=-sum(mu[less(i)]) end
  mu
end

# See Fleischmann-Janiszczak AAECC 1996 definition 2.1
function nconjugacy_classes(r::ClassType,WF,p)
  gets(r,:nClasses)do
    HF=r.CGs
    H=Group(HF)
    W=Group(WF)
    P=deepcopy(ClosedSubsets(W))
    elts=P.prop[:elements]
    o=filter(i->sort(elts[i].^HF.phi)==elts[i],1:length(P))
    o=filter(i->issubset(inclusiongens(H),elts[i]),o)
    P=restricted(P,o)
    P.prop[:elements]=elts[o]
# here P poset of HF.phi-stable closed subsets containing roots(H)
    InfoChevie("# ", HF, "==>", P, "")
    l=map(x->spets(reflection_subgroup(W, x), HF.phi), P.prop[:elements])
    b=collect(keys(factor(PermRoot.BadNumber(W))))
    l=map(l)do RF
      R=Group(RF)
      res=reduce(*,map(p->Mvp(:q)-p[2],filter(y->y[1]==1,degrees(RF)));init=Mvp(1))
      if semisimplerank(R)==0 return res end
      d=SmithNormalFormIntegerMat(simpleroots(R))
      d=filter(x->!(x in [0,1,p]),toM(d))
      if isempty(d) return res end
      res*prod(d)do i
                   if i==2 && (2 in b && p!=2) return 2 end
                   if mod(p,i)==1 return i end
                   Mvp(Symbol("q_",i))
                 end
    end
    mu=moebius(P)
    n=stabilizer(W,sort(inclusiongens(H)))
    n = sum(mu.*l) // length(centralizer(n, HF.phi))
    InfoChevie("==>", n, "\n")
    n
  end
end
end
