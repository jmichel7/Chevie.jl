module HasType

export reflection_name, diagram,
  charname, codegrees, ComplexReflectionGroup,
  field, getchev, Cartesian, ExtendedCox, traces_words_mats

using ..Gapjm
#-----------------------------------------------------------------------
function field(t::TypeIrred)
  if haskey(t,:orbit)
    orderphi=order(t.twist)
    t=t.orbit[1]
  else
    orderphi=1
  end
  s=t.series
  if s in [:A,:B,:D] 
     if orderphi==1 return (s,PermRoot.rank(t))
     elseif orderphi==2 
       if s==:B return (Symbol("2I"),4)
       else return (Symbol(2,s),PermRoot.rank(t))
       end
     elseif orderphi==3 return (Symbol("3D4"),)
     end
  elseif s in [:E,:F,:G]
    if orderphi==1 return (Symbol(s,PermRoot.rank(t)),) 
    elseif s==:G return (Symbol("2I"),6)
    else return (Symbol(orderphi,s,PermRoot.rank(t)),) 
    end
  elseif s==:ST 
    if haskey(t,:ST)
      if orderphi!=1 return (Symbol(orderphi,"G",t.ST),)
      elseif 4<=t.ST<=22 return (:G4_22,t.ST)
      else return (Symbol(string("G",t.ST)),)
      end
    elseif orderphi!=1
      return (:timp, t.p, t.q, t.rank)
    else
      return (:imp, t.p, t.q, t.rank)
    end
  elseif s==:I return (orderphi==1 ? :I : Symbol("2I"),t.bond)
  else return (Symbol(string(s,PermRoot.rank(t))),) 
  end
end

const needcartantype=Set([:PrintDiagram,
                          :ReflectionName,
                          :UnipotentClasses,
                          :WeightInfo])

function getchev(t::TypeIrred,f::Symbol,extra...)
  d=field(t)
# println("d=$d f=$f extra=$extra")
  o=chevieget(d[1],f)
  if o isa Function
#   o(vcat(collect(d)[2:end],collect(extra))...)
    if haskey(t,:cartantype) && f in needcartantype
#     println("args=",(d[2:end]...,extra...,t.cartantype))
      o(d[2:end]...,extra...,t.cartantype)
     else o(d[2:end]...,extra...)
    end
  else o
  end
end

function getchev(W,f::Symbol,extra...)
  [getchev(ti,f::Symbol,extra...) for ti in refltype(W)]
end

#-----------------------------------------------------------------------
struct Unknown
end

Base.:+(a,b::Unknown)=b
Base.:+(b::Unknown,a)=b
Base.:*(a,b::Unknown)=b
Base.:*(b::Unknown,a)=b
Base.zero(a::Unknown)=0
#-----------------------------------------------------------------------
function Cartesian(a::AbstractVector...)
  reverse.(vec(collect.(Iterators.product(reverse(a)...))))
end

charname(t::TypeIrred,p;TeX=false,opt...)=getchev(t,:CharName,p,
                           TeX ? Dict(:TeX=>true) : Dict())

charname(W,x;TeX=false,opt...)=join(map((t,p)->charname(t,p;TeX=TeX,opt...),
                           refltype(W),x),",")

function PositionCartesian(l,ind)
  res=prod=1
  for i in length(l):-1:1
    res+=(ind[i]-1)*prod
    prod*=l[i]
  end
  res
end

"""
`ComplexReflectionGroup(STnumber)`

`ComplexReflectionGroup(p,q,r)`

The  first form of `ComplexReflectionGroup`  returns the complex reflection
group which has Shephard-Todd number `STnumber`, see cite{ST54}. The second
form returns the imprimitive complex reflection group `G(p,q,r)`.

```julia-repl
julia> G=ComplexReflectionGroup(4)
G₄

julia> degrees(G)
2-element Array{Int64,1}:
 4
 6

julia> length(G)
24

julia> fakedegrees(G,Pol(:q))
7-element Array{Pol{Int64},1}:
 1       
 q⁴      
 q⁸      
 q⁷+q⁵   
 q⁵+q³   
 q³+q    
 q⁶+q⁴+q²

julia> ComplexReflectionGroup(2,1,6)
B₆
```
"""
function ComplexReflectionGroup(i::Int)
  if i==23     return coxgroup(:H,3)
  elseif i==28 return coxgroup(:F,4)
  elseif i==30 return coxgroup(:H,4)
  elseif i==35 return coxgroup(:E,6)
  elseif i==36 return coxgroup(:E,7)
  elseif i==37 return coxgroup(:E,8)
  end
  t=TypeIrred(Dict(:series=>:ST,:ST=>i))
  r=getchev(t,:GeneratingRoots)
  cr=getchev(t,:GeneratingCoRoots)
  if isnothing(cr)
   cr=map(coroot,r,map(x->E(;r=x),getchev(t,:EigenvaluesGeneratingReflections)))
  end
  PRG(r,cr)
end

function ComplexReflectionGroup(p,q,r)
  if !iszero(p%q) || p<=0 || r<=0 || (r==1 && q!=1) 
   error("ComplexReflectionGroup(p,q,r) must satisfy: q|p, r>0, and r=1 => q=1")
  end
  if p==1 return coxgroup(:A,r-1)
  elseif p==2 
    if q==2 return coxgroup(:D,r)
    else return coxgroup(:B,r) end
  elseif p==q && r==2 return coxgroup(:I,2,r)
  end
  t=TypeIrred(Dict(:series=>:ST,:p=>p,:q=>q,:rank=>r))
  r=getchev(t,:GeneratingRoots)
  cr=map(coroot,r,map(x->E(;r=x),getchev(t,:EigenvaluesGeneratingReflections)))
  PRG(r,map(x->convert.(Cyc{Rational{Int}},x),cr))
end

"""
`degrees(W)`

returns  a list  holding the  degrees of  `W` as  a reflection group on the
vector  space `V` on which  it acts. These are  the degrees `d₁,…,dₙ` where
`n`  is the dimension of  `V` of the basic  invariants of `W` in `SV`. They
reflect various properties of `W`; in particular, their product is the size
of `W`.

```julia-repl
julia> W=ComplexReflectionGroup(30)
H₄

julia> degrees(W)
4-element Array{Int64,1}:
  2
 12
 20
 30

julia> length(W)
14400
```
"""
function Gapjm.degrees(W::Group)
  gets(W,:degrees)do W
    vcat(fill(1,rank(W)-semisimplerank(W)),degrees.(refltype(W))...)
  end
end

"""
'degrees(WF::Spets)'

Let  `W` be  the group  of the  reflection coset  `WF`, and  let `V` be the
vector  space  of  dimension  'rank(W)'  on  which `W` acts as a reflection
group.  Let  `f₁,…,fₙ`  be  the  basic  invariants  of `W` on the symmetric
algebra  `SV` of `V`;  they can be  chosen so they  are eigenvectors of the
matrix  `WF.F`. The corresponding  eigenvalues are called  the *factors* of
`F` acting on `V`; they characterize the coset --- they are equal to 1 only
for  the trivial  coset. The  *generalized degrees*  of `WF`  are the pairs
formed of the reflection degrees and the corresponding factor.

```julia-repl
julia> W=coxgroup(:E,6)
E₆

julia> WF=spets(W)
E₆

julia> phi=W(6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1);

julia> HF=subspets(WF,2:5,phi)
³D₄Φ₃

julia> Diagram(HF)
ϕ acts as (1,2,4) on the component below
  O 2
  ￨
O—O—O
1 3 4

julia> degrees(HF)
6-element Array{Tuple{Int64,Cyc{Int64}},1}:
 (1, ζ₃) 
 (1, ζ₃²)
 (2, 1)  
 (4, ζ₃) 
 (6, 1)  
 (4, ζ₃²)
```
"""
function Gapjm.degrees(W::Spets)
  gets(W,:degrees)do W
    vcat(map(x->(1,x),E.(roots(torusfactors(W)))),degrees.(refltype(W))...)
  end
end

function Gapjm.degrees(t::TypeIrred)
  if !haskey(t,:orbit) return getchev(t,:ReflectionDegrees) end
  d=getchev(t.orbit[1],:ReflectionDegrees)
# Let  t.scalar=[s_1,..,s_r],  where  r=length(t.orbit)  and  let  p be the
# PhiFactor   of  t.twist  associated  to  the  reflection  degree  d_i  of
# t.orbit[1].   If   G0   is   the   Spets  described  by  t.orbit[1],  and
# G1:=Ennola(Product(t.scalar),G0)  then G is isomorphic  to the descent of
# scalars  of G1. According to spets 1.5, a Phifactor of Ennola(zeta,G0) is
# \zeta^{d_i}  times  that  of  G0;  and  by  spets  1.5  or [Digne-Michel,
# parabolic A.1] those of an a-descent of scalars are
# \zeta_a^j\zeta_i^{1/a} (all the a-th roots of \zeta_i).
  if order(t.twist)>1 
   f=getchev(t,:PhiFactors)
   if isnothing(f) return f end
  else f=fill(1,length(d))
  end
  if haskey(t,:scalar)
    p=prod(t.scalar) 
    f=[f[i]*p^d[i] for i in eachindex(d)]
  end
  f=collect(zip(d,f))
  a=length(t.orbit)
  if a==1 return f end
  vcat(map(f)do (d,e) map(x->(d,x),root(e,a).*E.(a,0:a-1)) end...)
end

function codegrees(t::TypeIrred)
  if !haskey(t,:orbit)
    d=getchev(t,:ReflectionCoDegrees)
    if !isnothing(d) return d
    else
      d=getchev(t,:ReflectionDegrees)
      return reverse(maximum(d).-d)
    end
  end
  d=getchev(t,:ReflectionCoDegrees)
  if isnothing(d)
    d=getchev(t.orbit[1],:ReflectionDegrees)
    a=argmax(d)
    d=reverse(d[a].-d)
    if order(t.twist)==1
      f=fill(1,length(d))
    else
      f=getchev(t,:PhiFactors)
      if isnothing(f) return f end
      f=reverse(map(x->f[a]//x,f))
    end
    d=zip(d,f)
  elseif order(t.twist)==1
    d=zip(d,fill(1,length(d)))
  end
  if haskey(t,:scalar)
    f=prod(t.scalar) 
    d=[(deg,eps*f^deg) for (deg,eps) in d]
  end
  a=length(t.orbit)
  if a==1 return d end
  vcat(map(d)do (d,e) map(x->(d,x),root(e,a).*E.(a,0:a-1)) end...)
end

function codegrees(W::Group)
  vcat(fill(0,rank(W)-semisimplerank(W)),collect.(codegrees.(refltype(W)))...)
end

function codegrees(W::Spets)
  vcat(map(x->(-1,x),E.(inv.(roots(torusfactors(W))))),
         collect.(codegrees.(refltype(W)))...)
end

function diagram(W)
  for t in refltype(W)
    getchev(t,:PrintDiagram,t[:indices],
               getchev(t,:ReflectionName,Dict()))
  end
end

nr_conjugacy_classes(W)=prod(getchev(W,:NrConjugacyClasses))

function reflection_name(io::IO,W)
  opt=IOContext(io,:TeX=>true).dict
  r=join(getchev(W,:ReflectionName,opt),"×")
  fromTeX(io,r)
end

#----------------------------------------------------------------------
# correct translations of GAP3 functions

include("../tools/gap3support.jl")
include("cheviesupport.jl")

function pad(s::String, i::Int)
  if i>0 return lpad(s,i)
  else return rpad(s,-i)
  end
end

pad(s::String)=s

PrintToSring(s,v...)=sprint(show,v...)

function Replace(s,p...)
# print("Replace s=$s p=$p")
  for (src,tgt) in (p[i]=>p[i+1] for i in 1:2:length(p))
    i=0
    while i+length(src)<=length(s)
     if src==s[i+(1:length(src))]
        if tgt isa String
          s=s[1:i]*tgt*s[i+length(src)+1:end]
        else
          s=vcat(s[1:i],tgt,s[i+length(src)+1:end])
        end
        i+=length(tgt)
      else
        i+=1
      end
    end
  end
# println("=>",s)
  s
end

struct ExtendedCox{T<:FiniteCoxeterGroup}
  group::T
  F0s::Vector{Matrix{Int}}
  phis::Vector{<:Perm}
end

function ExtendedCox(W::FiniteCoxeterGroup,F0s::Vector{Matrix{Int}})
  ExtendedCox(W,F0s,isempty(F0s) ? Perm{Int}[] : map(F->PermX(W.G,F),F0s))
end

function Base.:*(a::ExtendedCox,b::ExtendedCox)
  if isempty(gens(a.group)) return b
  elseif isempty(gens(b.group)) return a
  end
  id(r)=one(fill(0,r,r))
  ExtendedCox(a.group*b.group,vcat(
                   map(m->cat(m,id(rank(b.group)),dims=(1,2)),a.F0s),
                   map(m->cat(id(rank(b.group)),m,dims=(1,2)),b.F0s)))
end

Base.show(io::IO,W::ExtendedCox)=print(io,"Extended($(W.group),$(W.F0s))")

ExtendedReflectionGroup(W,mats::Vector{Matrix{Int}})=ExtendedCox(W,mats)
ExtendedReflectionGroup(W,mats::Matrix{Int})=ExtendedCox(W,[mats])
ExtendedReflectionGroup(W,mats::Vector{Vector{Int}})=ExtendedCox(W,[toM(mats)])

function ExtendedReflectionGroup(W,mats::Vector{Vector{Vector{Int}}})
  if isempty(mats)  ExtendedCox(W,empty([fill(0,0,0)]))
  elseif isempty(mats[1]) ExtendedCox(W,fill(fill(0,0,0),length(mats)))
  else ExtendedCox(W,toM.(mats))
  end
end

ExtendedReflectionGroup(W,p::Vector{Perm{Int}})=ExtendedCox(W,matX.(Ref(W),p))
ExtendedReflectionGroup(W,p::Perm)=ExtendedCox(W,[matX(W,p)])

function ExtendedReflectionGroup(W,mats::Vector{Any})
  if isempty(mats) ExtendedCox(W,empty([fill(0,0,0)]))
  else error("not empty")
  end
end

reflection_name(io::IO,W::ExtendedCox)=reflection_name(io,W.group)

CharTableSymmetric=Dict(:centralizers=>[
     function(n,pp) res=k=1;last=0
        for p in pp
          res*=p
          if p==last k+=1;res*=k
          else k=1
          end
          last=p
        end
        res
     end])

ChevieIndeterminate(a::Vector{<:Number})=Pol([1],0)
ChevieIndeterminate(a::Vector{<:Pol})=Mvp(:x)

"""
 EvalPolRoot(pol, x, n, p) compute pol(p*x^(1/n))
  
  The point of this routine is to avoid unnecessary root extractions
  during evaluation (e.g., if pol has no terms of odd degree and n=2,
  then no root extraction is necessary).
"""
function EvalPolRoot(pol::Pol,x,n,p)
  if isempty(pol.c) return 0 end
  P=vcat(fill(0,mod(pol.v,n)),pol.c)
  P=map(i->Pol(P[i:n:length(P)],div(pol.v-mod(pol.v,n),n))(x*p^n),1:n)
  j=findlast(!iszero,P)
  if isnothing(j) return 0 end
  pol=Pol(P[1:j],0)
  l=pol.v-1+filter(i->!iszero(pol.c[i]),eachindex(pol.c))
  push!(l,n)
  r=gcd(l...)
  pol=Pol(pol.c[1:r:length(pol.c)],div(pol.v,r))
  pol(GetRoot(x,div(n,r))*p^r)
end

function VcycSchurElement(arg...)
  n = length(arg[1])
  if length(arg) == 3
      data = arg[3]
      para = arg[1][data[:order]]
  else
      para = copy(arg[1])
  end
  monomial(mon)=prod(map(^,para,Int.(mon[1:n])))
  r = arg[2]
  if haskey(r, :rootUnity) && haskey(r,:root) error("cannot have both") end
  if haskey(r, :coeff) res = r[:coeff] else res = 1 end
  if haskey(r, :factor) res*=monomial(r[:factor]) end
  if haskey(r, :rootUnity)
    term=function(v)
     mon,pol=v
     tt=monomial(mon)
     if length(mon)==n+1 tt*=(r[:rootUnity]^data[:rootUnityPower])^mon[n+1] end
     cyclotomic_polynomial(pol)(tt)
    end
  elseif haskey(r, :root)
    term=function(v)
      mon,pol=v
     if length(mon)==n return Pol([cyclotomic_polynomial(pol)(monomial(mon))],0)
     else return cyclotomic_polynomial(pol)(Pol([monomial(mon)],mon[n+1]))
     end
    end
  else term=v->cyclotomic_polynomial(v[2])(monomial(v[1]))
  end
  res*=prod(term,r[:vcyc])
  if !haskey(r, :root) return res end
  den=lcm(denominator.(r[:root])...)
  root=monomial(den*r[:root])
  if haskey(r, :rootCoeff) root*=r[:rootCoeff] end
  EvalPolRoot(res, root, den, data[:rootPower])
end

"""
`CycPol(v::AbstractVector)`

This  form is a fast  and efficient way of  specifying a `CycPol` with only
positive multiplicities: `v` should be a vector. The first element is taken
as  a  the  `.coeff`  of  the  `CycPol`,  the  second  as the `.valuation`.
Subsequent  elements are rationals `i//d`  representing `(q-E(d)^i)` or are
integers `d` representing `Φ_d(q)`.

```julia-repl
julia> CycPol([3,-5,6,3//7])
3q⁻⁵Φ₆(q-ζ₇³)
```
"""
function CycPols.CycPol(v::AbstractVector)
  coeff=v[1]
  valuation=convert(Int,v[2])
  vv=Pair{Root1,Int}[]
  v1=convert.(Rational{Int},v[3:end])
  for i in v1
    if denominator(i)==1
      k=convert(Int,i)
      for j in prime_residues(k) push!(vv,Root1(;r=j//k)=>1) end
    else
      push!(vv,Root1(;r=i)=>1)
    end
  end
  CycPol(coeff,valuation,ModuleElt(vv))
end

"""
`traces_words_mats(mats,words)`
    
given  a list `mats`  of matrices and  a list `words`  of words returns the
list of traces of the corresponding products of the matrices

```julia-repl
julia> W=coxgroup(:F,4)
F₄

julia> r=classinfo(W)[:classtext];

julia> R=representation(W,17)
4-element Array{Array{Int64,2},1}:
 [-1 -1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] 
 [1 0 0 0; -1 -1 -1 0; 0 0 1 0; 0 0 0 1]
 [1 0 0 0; 0 1 0 0; 0 -2 -1 -1; 0 0 0 1]
 [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 -1 -1] 

julia> traces_words_mats(R,r)==CharTable(W).irr[17,:]
true
```
"""
function traces_words_mats(mats,words)
  mats=improve_type.(mats)
  dens=map(x->1,mats)
  if all(m->all(x->x isa Rational,m),mats) 
    dens=map(m->lcm(denominator.(m)),mats)
    mats=map((m,d)->Int.(m.*d),mats,dens)
  end
  words=convert.(Vector{Int},words)
  tr(m)=sum(i->m[i,i],axes(m,1))
  trace(w)=tr(prods[w])//prod(dens[w])
  prods=Dict{Vector{Int},eltype(mats)}(Int[]=>mats[1]^0)
  for i in eachindex(mats) prods[[i]]=mats[i] end
  res=map(words)do w
    i=0
    while haskey(prods,w[1:i]) 
      if i==length(w) return trace(w) end
      i+=1 
    end
    while i<=length(w) 
      prods[w[1:i]]=prods[w[1:i-1]]*mats[w[i]]
      i+=1
    end
#   println(prod(dens[w]))
    trace(w)
  end
end

function ImprimitiveCuspidalName(S)
  r=RankSymbol(convert(Vector{Vector{Int}},S))
  d=length(S)
  s=joindigits(map(length,S))
  if r==0 return "" end
  if sum(length,S)%d==1 # G(d,1,r)
    if r==1 return d==3 ? "Z_3" : "Z_{$d}^{$s}"
    else return "G_{$d,1,$r}^{$s}"
    end
  else # G(d,d,r)
    if r==2
      if d==4 return "B_2"
      elseif d==6 
        p=Dict("212010"=>"-1","221001"=>"1",
               "211200"=>"\\zeta^2","220110"=>"\\zeta_3")
        return "G_2[$(p[s])]"
      else p=CHEVIE.R("SymbolToParameter","I")(S);
	return "I_2($d)",FormatGAP(p)
      end
      elseif r==3 && d==3 
        return "G_{3,3,3}[\\zeta_3"* (s=="300" ? "" : "^2")*"]"
      elseif r==3 && d==4 
        return "G_{4,4,3}["* (s=="3010" ? "" : "-")*"\\zeta_4]"
    else return "G_{$d,$d,$r}^{$s}"
    end
  end
end

function BDSymbols(n,d)
  n-=div(d^2,4)
  if n<0 return Vector{Vector{Int}}[] end
  if d>0 return map(x->symbol_partition_tuple(x,d),partition_tuples(n,2)) end
   return map(chevieget(:D,:symbolcharparam),
              chevieget(:imp,:CharInfo)(2,2,n)[:charparams])
end

const src=[ 
#  "compat3", 
"cmp4_22", "cmplxg24", "cmplxg25", "cmplxg26", 
"cmplxg27", "cmplxg29", "cmplxg31", "cmplxg32", "cmplxg33", "cmplxg34", 
"cmplximp", "coxh3", "coxh4", "coxi", 
"weyla", "weylbc", "weyld", "weyl2a", 
"weyl2d",
"cox2i", "weyl2e6", "weyl2f4", "weyl3d4",
"weyle6", "weyle7", "weyle8", "weylf4", "weylg2", "exceptio"]

for f in src
  println("reading tbl/$f.jl")
  include("tbl/$f.jl")
end
include("table2.jl")

end
