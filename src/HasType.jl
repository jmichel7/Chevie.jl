module HasType

export charinfo, classinfo, reflection_name, diagram, chartable,
  representation, fakedegrees, unipotent_characters, unipotent_classes,
  schur_elements, charname, codegrees, ComplexReflectionGroup,
  chevieget, field, getchev, Cartesian

using Gapjm
#-----------------------------------------------------------------------
const chevie=Dict()

# extensions to get closer to GAP semantics
Base.:*(a::Array,b::Pol)=a .* Ref(b)
Base.:*(a::Pol,b::Array)=Ref(a) .* b
Base.:*(a::AbstractVector,b::AbstractVector)=sum(a.*b)
Base.:-(a::AbstractVector,b::Number)=a .- b
Base.:+(a::Integer,b::AbstractVector)=a .+ b
Base.:+(a::AbstractVector,b::Number)=a .+ b
Base.:+(a::AbstractVector,b::Pol)=a .+ Ref(b)
#Base.:*(a::Mvp, b::Array)=Ref(a).*b
#include("mvp.jl")
#Base.:*(b::Array,a::Mvp)=b.*Ref(a)
Base.getindex(s::String,a::Vector{Any})=getindex(s,Int.(a))
Cycs.:^(a::Cyc,b::Rational)=a^Int(b)
Base.:^(m::AbstractMatrix,n::AbstractMatrix)=inv(n*E(1))*m*n
Base.:^(m::Vector{<:Vector{<:Number}},n::Matrix{<:Number})=inv(n)*hcat(m...)*n
Base.isless(a::Array,b::Number)=true
Base.isless(b::Number,a::Array)=false
Base.getindex(a::Symbol,i::Int)=string(a)[i]
Base.length(a::Symbol)=length(string(a))
Base.copy(x::Char)=x

function chevieget(t::Symbol,w::Symbol)
  if haskey(chevie[t],w) return chevie[t][w] end
  println("chevie[$t] has no $w")
end

function chevieset(t::Symbol,w::Symbol,o::Any)
  if !haskey(chevie,t) chevie[t]=Dict{Symbol,Any}() end
  chevie[t][w]=o
end

function chevieset(t::Vector{String},w::Symbol,f::Function)
  for s in t 
    println("set $s $w")
    chevieset(Symbol(s),w,f(Symbol(s))) 
  end
end

function field(t::TypeIrred)
  if haskey(t,:orbit)
    orderphi=order(t[:twist])
    t=t[:orbit][1]
  else
    orderphi=1
  end
  s=t[:series]
  if s in [:A,:B,:D] 
     if orderphi==1 return (s,length(t[:indices]))
     elseif orderphi==2 return (Symbol(2,s),length(t[:indices]))
     elseif orderphi==3 return (Symbol("3D4"),)
     end
  elseif s in [:E,:F,:G]
    if orderphi==1 return (Symbol(s,length(t[:indices])),) 
    else return (Symbol(orderphi,s,length(t[:indices])),) 
    end
  elseif s==:ST 
    if haskey(t,:ST)
      if orderphi!=1 return (Symbol(orderphi,"G",t[:ST]),)
      elseif 4<=t[:ST]<=22 return (:G4_22,t[:ST])
      else return (Symbol(string("G",t[:ST])),)
      end
    elseif orderphi!=1
      return (:timp, t[:p], t[:q], t[:rank])
    else
      return (:imp, t[:p], t[:q], t[:rank])
    end
  elseif s==:I 
    return (orderphi==1 ? :I : :(2I),t[:bond])
  else return (Symbol(string(s,length(t[:indices]))),) 
  end
end

const needcartantype=Set([:PrintDiagram,:ReflectionName,:UnipotentClasses])

function getchev(t::TypeIrred,f::Symbol,extra...)
  d=field(t)
# println("d[1]=$(d[1]) f=$f")
  o=chevieget(d[1],f)
  if o isa Function
#   o(vcat(collect(d)[2:end],collect(extra))...)
    if haskey(t,:cartantype) && f in needcartantype
      o(d[2:end]...,extra...,t[:cartantype])
    else o(d[2:end]...,extra...)
    end
  else o
  end
end

function getchev(W,f::Symbol,extra...)
  [getchev(ti,f::Symbol,extra...) for ti in refltype(W)]
end

#-----------------------------------------------------------------------

function Cartesian(a::AbstractVector...)
  reverse.(vec(collect.(Iterators.product(reverse(a)...))))
end

braid_relations(W)=impl1(getchev(W,:BraidRelations))

function codegrees(W)
  vcat(map(refltype(W)) do t
    cd=getchev(t,:ReflectionCoDegrees)
    if isnothing(cd)
      cd=getchev(t,:ReflectionDegrees)
      maximum(cd).-cd
    else cd
    end
  end...)
end

impl1(l)=length(l)==1 ? l[1] : error("implemented only for irreducible groups")

charname(W,x;TeX=false,opt...)=join(map((t,p)->getchev(t,:CharName,p,
                           TeX ? Dict(:TeX=>true) : Dict()),refltype(W),x),",")

cartfields(p,f)=Cartesian(getindex.(p,f)...)

function PositionCartesian(l,ind)
  res=prod=1
  for i in length(l):-1:1
    res+=(ind[i]-1)*prod
    prod*=l[i]
  end
  res
end

allhaskey(v::Vector{<:Dict},k)=all(d->haskey(d,k),v)

function charinfo(t::TypeIrred)
  c=deepcopy(getchev(t,:CharInfo))
  c[:positionId]=c[:extRefl][1]
  c[:positionDet]=c[:extRefl][end]
  c[:charnames]=map(c[:charparams]) do p
     getchev(t,:CharName,p,Dict(:TeX=>true))
  end
  c
end

function charinfo(W)::Dict{Symbol,Any}
  gets(W,:charinfo) do W
    p=charinfo.(refltype(W))
    if isempty(p) return Dict(:a=>[0],:A=>[0],:b=>[0],:B=>[0],:positionId=>1,
      :positionDet=>1,:charnames=>["Id"],:extRefl=>[1],:charparams=>[[]])
    end
    if length(p)==1 res=copy(p[1]) else res=Dict{Symbol, Any}() end
    res[:charparams]=cartfields(p,:charparams)
    if length(p)==1 return res end
    res[:charnames]=map(l->join(l,","),cartfields(p,:charnames))
    for f in [:positionId, :positionDet]
      if allhaskey(p,f)
       res[f]=PositionCartesian(map(x->length(x[:charparams]),p),getindex.(p,f))
      end
    end
    for f in [:b, :B, :a, :A]
      if allhaskey(p,f) res[f]=Int.(map(sum,cartfields(p,f))) end
    end
    if any(x->haskey(x, :opdam),p)
      res[:opdam]=map(x->haskey(x,:opdam) ? x[:opdam] : Perm(), p)
      gt=Cartesian(map(x->1:length(x[:charparams]), p))
      res[:opdam]=PermListList(gt, map(t->map((x,i)->x^i,t,res[:opdam]),gt))
    end
    res
  end
end

function classinfo(t::TypeIrred)
  cl=deepcopy(getchev(t,:ClassInfo))
  inds=t[:indices]
  cl[:classtext]=map(x->inds[x],cl[:classtext])
  cl[:classes]=Int.(cl[:classes])
  cl
end

function classinfo(W)::Dict{Symbol,Any}
  gets(W,:classinfo) do W
    tmp=map(classinfo,refltype(W))
    if isempty(tmp) return Dict(:classtext=>[Int[]],:classnames=>[""],
                      :classparams=>[Int[]],:orders=>[1],:classes=>[1])
    end
    if any(isnothing, tmp) return nothing end
    if length(tmp)==1 res=copy(tmp[1]) else res=Dict{Symbol, Any}() end
    res[:classtext]=map(x->vcat(x...),cartfields(tmp,:classtext))
    res[:classnames]=map(join,cartfields(tmp,:classnames))
    if allhaskey(tmp, :classparam)
      res[:classparams]=cartfields(tmp,:classparams)
    end
    if allhaskey(tmp,:orders)
      res[:orders]=map(lcm, cartfields(tmp,:orders))
    end
    if allhaskey(tmp,:classes)
      res[:classes]=map(prod, cartfields(tmp,:classes))
    end
    res
  end
end

class_rep(W::FiniteCoxeterGroup)=map(x->W(x...),classinfo(W)[:classtext])

function chartable(t::TypeIrred)
  ct=getchev(t,:CharTable)
  if haskey(ct,:irredinfo) names=getindex.(ct[:irredinfo],:charname)
  else                     names=charinfo(t)[:charnames]
  end
  if !haskey(ct,:classnames) merge!(ct,classinfo(t)) end
  CharTable(permutedims(Cyc{Int}.(hcat(ct[:irreducibles]...))),names,
   ct[:classnames],Int.(ct[:centralizers]),ct[:identifier])
end

function chartable(W)::CharTable
  gets(W,:chartable) do W
    ctt=chartable.(refltype(W))
    if isempty(ctt) 
      return CharTable(hcat(1),["Id"],["1"],[1],"$W")
    end
    charnames=join.(Cartesian(getfield.(ctt,:charnames)...),",")
    classnames=join.(Cartesian(getfield.(ctt,:classnames)...),",")
    centralizers=prod.(Cartesian(getfield.(ctt,:centralizers)...))
    identifier=join(getfield.(ctt,:identifier),"×")
    if length(ctt)==1 irr=ctt[1].irr 
    else irr=kron(getfield.(ctt,:irr)...)
    end
    CharTable(irr,charnames,classnames,centralizers,identifier)
  end
end

function chartable(H::HeckeAlgebra{C})where C
  W=H.W
  ct=impl1(getchev(W,:HeckeCharTable,H.para,
       haskey(H.prop,:rootpara) ? rootpara(H) : fill(nothing,length(H.para))))
  if haskey(ct,:irredinfo) names=getindex.(ct[:irredinfo],:charname)
  else                     names=charinfo(W)[:charnames]
  end
  CharTable(Matrix(permutedims(hcat(
       map(ch->convert.(C,ch),ct[:irreducibles])...))),names,
     ct[:classnames],map(Int,ct[:centralizers]),ct[:identifier])
end

function representation(H::HeckeAlgebra,i::Int)
  ct=impl1(getchev(H.W,:HeckeRepresentation,H.para,
    haskey(H.prop,:rootpara) ? rootpara(H) : fill(nothing,length(H.para)),i))
  map(x->hcat(x...),ct)
end

function schur_elements(H::HeckeAlgebra)
  W=H.W
  map(p->getchev(W,:SchurElement,p,H.para,
     haskey(H.prop,:rootpara) ? rootpara(H) : fill(nothing,length(H.para)))[1],
      first.(charinfo(W)[:charparams]))
end

function ComplexReflectionGroup(i::Int)
  if i in [23,28,30,35,36,37]
    if i==23     return coxgroup(:H,3)
    elseif i==28 return coxgroup(:F,4)
    elseif i==30 return coxgroup(:H,4)
    elseif i==35 return coxgroup(:E,6)
    elseif i==36 return coxgroup(:E,7)
    elseif i==37 return coxgroup(:E,8)
    end
    m=getchev(t,:CartanMat)
    n=one(hcat(m...))
    return PRG(map(i->n[i,:],axes(n,1)),m)
  end
  t=TypeIrred(Dict(:series=>:ST,:ST=>i))
  r=getchev(t,:GeneratingRoots)
  cr=getchev(t,:GeneratingCoRoots)
  if cr===nothing
    e=getchev(t,:EigenvaluesGeneratingReflections)
    cr=map((x,y)->coroot(x,y),r,map(x->E(denominator(x),numerator(x)),e))
  end
  PRG(r,cr)
end

function ComplexReflectionGroup(p,q,r)
 if p==1 if r==1 return coxgroup()
   else return coxgroup(:A,r-1)
   end
  elseif p==2 
   if q==2 if r==1 return coxgroup()
           elseif r==2 return coxgroup(:A,1)*coxgroup(:A,1)
           elseif r==3 return coxgroup(:A,3)
           else return coxgroup(:D,r)
           end
   elseif r==1 return coxgroup(:A,1)
   else return coxgroup(:B,r) end
  elseif p==q && r==2 return coxgroup(:I,2,r)
  end
 t=TypeIrred(Dict(:series=>:ST,:p=>p,:q=>q,:rank=>r))
  r=getchev(t,:GeneratingRoots)
  cr=getchev(t,:EigenvaluesGeneratingReflections)
  cr=map((x,y)->coroot(x,y),r,map(x->E(Root1(numerator(x),denominator(x))),cr))
  cr=map(x->convert.(Cyc{Rational{Int}},x),cr)
  PRG(r,cr)
end

Gapjm.degrees(W)=vcat(getchev(W,:ReflectionDegrees)...)

function diagram(W)
  for t in refltype(W)
    getchev(t,:PrintDiagram,t[:indices],
               getchev(t,:ReflectionName,Dict()))
  end
end

function fakedegree(W,p,q)
  typ=refltype(W)
  if isempty(typ) return one(q) end
  prod(map((t,p)->getchev(t,:FakeDegree,p,q),typ,p))
end

function fakedegrees(W,q)
  map(p->fakedegree(W,p,q),charinfo(W)[:charparams])
end

function WeightToAdjointFundamentalGroupElement(W,i)
  t=First(refltype(W),t->i in t[:indices])
  l=W.rootInclusion[t[:indices]]
  b=longest(W,l)*longest(W,Difference(l,[W.rootInclusion[i]]))
  Add(l,W.rootInclusion[Maximum(Filtered(1:length(W.roots),
    i->ForAll(1:semisimpleRank(W),j->j in t[:indices] || W.roots[i][j]==0)))])
  RestrictedPerm(b,l)
end

# returns a record containing minuscule coweights, decompositions
# (in terms of generators of the fundamental group)
function WeightInfo(W)
  l=map(refltype(W)) do t
    r=getchev(t,:WeightInfo)
    g=filter(i->sum(r[:decompositions][i])==1,
          1:length(r[:minusculeCoweights])) # generators of fundamental group
    println("g=$g")
    r[:ww]=map(x->WeightToAdjointFundamentalGroupElement(W,x),
               t[:indices][r[:minusculeCoweights[g]]])
    C=mod1.(inv(Rational.(toM(cartan(t)))))
    r[:csi]=zeros(Int,length(g),semisimplerank(W))
    r[:csi][:,t[:indices]]=C[r[:minusculeCoweights[g]]]
    r[:minusculeWeights]=t[:indices][r[:minusculeWeights]]
    r[:minusculeCoweights]=t[:indices][r[:minusculeCoweights]]
    r
  end
  res=Dict(minusculeWeights=>Cartesian(List(l,
                                        x->vcat(x[:minusculeWeights],[0]))),
    minusculeCoweights=>Cartesian(List(l,
                                       x->vcat(x[:minusculeCoweights],[0]))),
    decompositions=>List(Cartesian(List(l,x->vcat(x[:decompositions],
      [0*x.moduli]))),Concatenation),
    moduli=>Concatenation(List(l,x->x[:moduli])))
# centre of simply connected group: the generating minuscule coweights
# mod the root lattice
  res[:CenterSimplyConnected]=vcat(getindex.(l,:csi)...)
  res[:AdjointFundamentalGroup]=vcat(getindex.(l,:ww)...)
  n=length(res[:decompositions])-1
  res[:minusculeWeights]=List(res.minusculeWeights[1:n],x->filter(y->y!=0,x))
  res[:minusculeCoweights]=List(res.minusculeCoweights[1:n],x->filter(y->y!=0,x))
  res[:decompositions]=res[:decompositions][1:n]
  res
end

nr_conjugacy_classes(W)=prod(getchev(W,:NrConjugacyClasses))

PrintToSring(s,v...)=sprint(show,v...)

reflection_name(W,opt=Dict())=join(getchev(W,:ReflectionName,opt),"×")

function representation(W,i::Int)
  map(x->hcat(x...),impl1(getchev(W,:Representation,i)))
end

include("uch.jl")
include("ucl.jl")
#----------------------------------------------------------------------
# correct translations of GAP3 functions

Append(a::Vector,b::AbstractVector)=vcat(a,b)
Append(a::String,b::String)=a*b
Append(a::String,b::Vector{Char})=a*String(b)

function Collected(v)
  d=groupby(v,v)
  sort([[k,length(v)] for (k,v) in d])
end

function CollectBy(v,f)
  d=groupby(f,v)
  [d[k] for k in sort(collect(keys(d)))]
end

Inherit(a,b)=merge!(a,b)
function Inherit(a,b,c)
  for k in c a[Symbol(k)]=b[Symbol(k)] end
  a
end

function pad(s::String, i::Int)
  if i>0 return lpad(s,i)
  else return rpad(s,-i)
  end
end

pad(s::String)=s

function Position(a::Vector,b)
  x=findfirst(isequal(b),a)
  isnothing(x) ? false : x
end

function Position(a::String,b::String)
  x=findfirst(b,a)
  isnothing(x) ? false : x.start
end

function Position(a::String,b::Char)
  x=findfirst(isequal(b),a)
  isnothing(x) ? false : x
end

function PositionProperty(a::Vector,b::Function)
  r=findfirst(b,a)
  if isnothing(r) return false end
  r
end

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

AbsInt=abs
ApplyFunc(f,x)=f(x...)
Arrangements=arrangements
BetaSet=βSet
Binomial=binomial
CartanMat(s,a...)=cartan(Symbol(s),a...)
CharParams(W)=charinfo(W)[:charparams]
Concatenation(a::String...)=prod(a)
Concatenation(a::AbstractVector{<:AbstractVector})=vcat(collect.(a)...)
Concatenation(a::Vector,b::Tuple)=vcat(a,collect(b))
Concatenation(b...)=vcat(b...)
Combinations=combinations
Copy=deepcopy
CycPolFakeDegreeSymbol=fegsymbol
DefectSymbol=defectsymbol
DiagonalMat(v...)=cat(map(m->m isa Array ? m : hcat(m),v)...,dims=(1,2))
DiagonalMat(v::Vector{<:Number})=DiagonalMat(v...)
DiagonalOfMat(m)=[m[i,i] for i in axes(m,1)]
DivisorsInt=divisors
Dominates=dominates
Drop(a::Vector,i::Int)=deleteat!(copy(a),i)
EltWord(W,x)=W(x...)
Factors(n)=vcat([fill(k,v) for (k,v) in factor(n)]...)
Filtered(l,f)=isempty(l) ? l : filter(f,l)
First(a,b)=a[findfirst(b,a)]
Flat(v)=collect(Iterators.flatten(v))
ForAll(l,f)=all(f,l)
FullSymbol=fullsymbol
Hasse=hasse
HighestPowerFakeDegreeSymbol=degree_feg_symbol
HighestPowerGenericDegreeSymbol=degree_gendeg_symbol
IdentityMat(n)=map(i->one(rand(Int,n,n))[i,:],1:n)
function Ignore() end
InfoChevie2=print
IntListToString=joindigits
IsList(l)=l isa Vector
IsInt(l)=l isa Int ||(l isa Rational && denominator(l)==1)
Join(x,y)=join(x,y)
Join(x)=join(x,",")
Lcm(a...)=Lcm(collect(a))
Lcm(a::Vector)=lcm(Int.(a))
LowestPowerFakeDegreeSymbol=valuation_feg_symbol
LowestPowerGenericDegreeSymbol=valuation_gendeg_symbol
MatXPerm=matX
Minimum(v::AbstractVector)=minimum(v)
Minimum(a::Number,x...)=min(a,x...)
NrConjugacyClasses(W)=length(classinfo(W)[:classtext])
OnMatrices(a::Vector{<:Vector},b::Perm)=Permuted(map(x->Permuted(x,b),a),b)
OrderPerm=order
PartBeta=partβ
Partitions=partitions
PartitionTuples=partition_tuples

function PartitionTupleToString(n,a=Dict())
  if n[end] isa Vector return join(map(join,n),".") end
  r=repr(E(n[end-1],n[end]),context=:limit=>true)
  if r=="1" r="+" end
  if r=="-1" r="-" end
  join(map(join,n[1:end-2]),".")*r
end

PermListList(l1,l2)=Perm(sortperm(l2))^-1*Perm(sortperm(l1))
Permuted(a,b)=[a[i^b] for i in eachindex(a)]
PrimeResidues=prime_residues
Product(v)=isempty(v) ? 1 : prod(v)
Product(v,f)=isempty(v) ? 1 : prod(f,v)
Rank=rank
RankSymbol=ranksymbol
RecFields=keys
ReflectionSubgroup(W,I::AbstractVector)=reflection_subgroup(W,convert(Vector{Int},I))
RootInt(a,b)=floor(Int,a^(1/b))
RootsCartan=roots
Rotations(a)=circshift.(Ref(a),0:length(a)-1)
gapSet(v)=unique(sort(v))
SemisimpleRank(W)=semisimplerank(W)
ShiftBeta=shiftβ
SignInt=sign
Sort=sort!
SortBy(x,f)=sort!(x,by=f)
SPrint=string
StringSymbol=stringsymbol
StringToDigits(s)=map(y->Position("01234567890", y), collect(s)).-1
Sublist(a::Vector, b::AbstractVector)=a[b]
Sum(v::AbstractVector)=sum(v)
Sum(v::AbstractVector,f)=isempty(v) ? 0 : sum(f,v)
SymbolPartitionTuple=symbol_partition_tuple
SymbolsDefect(a,b,c,d)=symbols(a,b,d)
Torus(i::Int)=torus(i)
Base.union(v::Vector)=union(v...)
Value(p,v)=p(v)
function CoxeterGroup(S::String,s...)
 if length(s)==1 return coxgroup(Symbol(S),Int(s[1])) end
 coxgroup(Symbol(S),Int(s[1]))*coxgroup(Symbol(s[2]),Int(s[3]))
end
CoxeterGroup()=coxgroup()

struct ExtendedCox{T<:FiniteCoxeterGroup}
  group::T
  F0s::Vector{Matrix{Int}}
end

ExtendedReflectionGroup(W,mats::Vector{Matrix{Int}})=ExtendedCox(W,mats)
ExtendedReflectionGroup(W,mats::Matrix{Int})=ExtendedCox(W,[mats])
ExtendedReflectionGroup(W,mats::Vector{Vector{Int}})=
  ExtendedCox(W,[hcat(mats...)])

function ExtendedReflectionGroup(W,mats::Vector{Vector{Vector{Int}}})
  if isempty(mats)  ExtendedCox(W,empty([fill(0,0,0)]))
  elseif isempty(mats[1]) ExtendedCox(W,fill(fill(0,0,0),length(mats)))
  else ExtendedCox(W,map(x->hcat(x...),mats))
  end
end

ExtendedReflectionGroup(W,p::Vector{Perm{Int}})=ExtendedCox(W,matX.(Ref(W),p))
ExtendedReflectionGroup(W,p::Perm)=ExtendedCox(W,[matX(W,p)])

function ExtendedReflectionGroup(W,mats::Vector{Any})
  if isempty(mats) ExtendedCox(W,empty([fill(0,0,0)]))
  else error("not empty")
  end
end

reflection_name(W::ExtendedCox,a...)=reflection_name(W.group,a...)

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
  j=findlast(x->!iszero(x),P)
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
  monomial=mon->prod(map(^,para,mon[1:n]))
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

function CycPols.CycPol(v::AbstractVector)
  coeff=v[1]
  valuation=convert(Int,v[2])
  vv=Pair{Root1,Int}[]
  v1=convert.(Rational{Int},v[3:end])
  for i in v1
    if denominator(i)==1
       k=convert(Int,i)
       for j in prime_residues(k) push!(vv,Root1(j//k)=>1) end
    else
     push!(vv,Root1(i)=>1)
    end
  end
  CycPol(coeff,valuation,ModuleElt(vv))
end

function CharRepresentationWords(mats,words)
  mats=map(x->hcat(x...),mats)
  map(words)do w
    if isempty(w) return size(mats[1],1) end
    m=prod(mats[w])
    sum(i->m[i,i],axes(m,1))
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
  if d>0 return map(x->symbol_partition_tuple(x,d),PartitionTuples(n,2)) end
   return map(chevieget(:D,:symbolcharparam),
              chevieget(:imp,:CharInfo)(2,2,n)[:charparams])
end

CHEVIE=Dict{Symbol,Any}(:compat=>Dict(:MakeCharacterTable=>x->x,
                           :AdjustHeckeCharTable=>(x,y)->x,
        :ChangeIdentifier=>function(tbl,n)tbl[:identifier]=n end))
include("families.jl")

#-------------------------------------------------------------------------
#  dummy translations of GAP3 functions
Format(x)=string(x)
FormatTeX(x)=repr(x,context=:TeX=>true)
FormatGAP(x)=repr(x)
Format(x,opt)= haskey(opt,:TeX) ? FormatTeX(x) : string(x)

function ReadChv(s) end
Group(v...)=v
ComplexConjugate(v)=v
GetRoot(x::Cyc,n::Number=2,msg::String="")=root(x,n)
GetRoot(x::Integer,n::Number=2,msg::String="")=root(x,n)
GetRoot(x::Pol,n::Number=2,msg::String="")=root(x,n)
function GetRoot(x,n::Number=2,msg::String="")
  error("GetRoot($x,$n) not implemented")
end

Unbind(x)=x
#-------------------------------------------------------------------------

include("tables.jl")
chevie[:D][:CharTable]=n->chevie[:imp][:CharTable](2,2,n)
chevie[:B][:CharTable]=n->chevie[:imp][:CharTable](2,1,n)
chevie[:A][:CharTable]=function(n)
  ct=chevie[:imp][:CharTable](1,1,n+1)
  ct[:irredinfo]=map(x->Dict(:charname=>joindigits(x)),chevie[:A][:CharInfo](n)[:charparams])
  ct
end
chevie[:A][:HeckeCharTable]=(n,para,root)->chevie[:imp][:HeckeCharTable](1,1,n+1,para,root)
chevie[:A][:FakeDegree]=(n,p,q)->fegsymbol([βSet(p)])(q)
chevie[:D][:HeckeCharTable]=(n,para,root)->chevie[:imp][:HeckeCharTable](2,2,n,para,root)
chevie[:imp][:PowerMaps]=(p,q,r)->[]
chevie[:imp][:GeneratingRoots]=function(p,q,r)
  if q==1 roots=[vcat([1],fill(0,r-1))]
  else
    if q!=p roots=[vcat([1],fill(0,r-1))*E(1)] end
    v=vcat([-E(p),1],fill(0,r-2))
    if r==2 && q>1 && q%2==1 v*=E(p) end
    if q == p roots = [v] else push!(roots, v) end
  end
  for i=2:r
    v=fill(0,r)
    v[[i-1,i]]=[-1,1]
    push!(roots, v)
  end
  return roots
end
chevie[:B][:UnipotentClasses]=function(r,char,ctype)
  part2dynkin=function(part)
    p=sort(vcat(map(d->1-d:2:d-1, part)...))
    p=p[div(3+length(p),2):end]
    if ctype==1 res=[2*p[1]]
    else res=[p[1]]
    end
    append!(res,p[2:end]-p[1:end-1])
  end
  addSpringer1=function(s,cc)
    ss=First(uc[:springerSeries],x->x[:defect]==DefectSymbol(s[:symbol]))
    if s[:sp] == [[], []] p = 1
    elseif s[:sp] == [[1], []] p = 2
    elseif s[:sp] == [[], [1]] p = 1
    else p = Position(CharParams(ss[:relgroup]), [s[:sp]])
    end
    ss[:locsys][p] = [length(uc[:classes]), Position(CharParams(cc[:Au]),
      map(x->x ? [1, 1] : [2], s[:Au]))]
  end
  if ctype == ER(2)
      ctype = 2
      char = 2
  end
  if char == 2 ss = XSP(4, 2, r)
  elseif ctype == 1 ss = XSP(2, 1, r)
  else ss = XSP(2, 0, r)
  end
  l = union(map(c->map(x->[DefectSymbol(x[:symbol]), Sum(x[:sp], Sum)], c), ss))
  sort!(l,by=x->[AbsInt(x[1]),-SignInt(x[1])])
  uc = Dict{Symbol, Any}(:classes => [], :springerSeries => map(function(d)
    res = Dict{Symbol, Any}(:relgroup => CoxeterGroup("C", d[2]), :defect => d[1], :levi => 1:r - d[2])
    res[:locsys] = fill([0, 0],NrConjugacyClasses(res[:relgroup]))
    if char == 2 res[:Z] = [1]
    elseif ctype == 1 res[:Z] = [(-1) ^ (r - d[2])]
    elseif IsInt(ER(2 * (r - d[2]) + 1)) res[:Z] = [1]
    else res[:Z] = [-1]
    end
    return res
  end, l))
  if char != 2
    symbol2para = function(S)
      c=sort(vcat(S...))
      i=1
      part=Int[]
      d=mod(ctype,2)
      while i <= length(c)
        if i == length(c) || c[i + 1] - c[i] > 0
            push!(part, (2 * (c[i] - (i - 1)) + 1) - d)
            i = i + 1
        else
            l = 2 * (c[i] - (i - 1)) - d
            part = Append(part, [l, l])
            i = i + 2
        end
      end
      reverse(filter(y->!iszero(y),sort(part)))
    end
  else
    symbol2para = function (S,)
      c=sort(vcat(S...))
      i = 1
      part = Int[]
      ex = Int[]
      while i <= length(c)
          if i == length(c) || c[i + 1] - c[i] > 1
              push!(part, 2 * (c[i] - 2 * (i - 1)))
              i = i + 1
          elseif c[i] == c[i + 1]
              l = 2 * (c[i] - 2 * (i - 1)) - 2
              part = Append(part, [l, l])
              push!(ex, l)
              i = i + 2
          elseif c[i] + 1 == c[i + 1]
              l = 2 * (c[i] - 2 * (i - 1)) - 1
              part = Append(part, [l, l])
              i = i + 2
          end
      end
      [reverse(filter(y->y!=0,sort(part))), ex]
    end
  end
  if char == 2 ctype = 1 end
  for cl = ss
      cc = Dict{Symbol, Any}(:parameter => symbol2para((cl[1])[:symbol]))
      cc[:Au] = CoxeterGroup(Concatenation(map(x->["A",1], cl[1][:Au]))...)
      if char != 2
          cc[:dynkin] = part2dynkin(cc[:parameter])
          cc[:name] = joindigits(cc[:parameter])
      else
          ctype = 1
          cc[:dimBu] = (cl[1])[:dimBu]
          cc[:name] = Join(map(function (x,)
              res = joindigits(fill(0, max(0, (1 + x[2]) - 1)) + x[1], "[]")
              if x[1] in cc[:parameter][2] return string("(", res, ")") end
              return res
          end, reverse(Collected(cc[:parameter][1]))), "")
      end
      cc[:red] = coxgroup()
      if char == 2 j = cc[:parameter][1]
      else j = cc[:parameter]
      end
      for j in Collected(j)
          if mod(j[1], 2) == mod(ctype, 2)
            cc[:red] = cc[:red] * CoxeterGroup("C", div(j[2],2))
          elseif mod(j[2], 2) != 0
            if j[2]>1 cc[:red]*=CoxeterGroup("B", div(j[2] - 1,2)) end
           elseif j[2]>2 cc[:red]*=CoxeterGroup("D", div(j[2], 2))
          else cc[:red]*=Torus(1)
          end
      end
      push!(uc[:classes], cc)
      for s in cl addSpringer1(s,cc) end
  end
  uc[:orderClasses] = Hasse(Poset(map((x->begin
      map(function (y,)
        if char != 2 return Dominates(y[:parameter], x[:parameter]) end
        m = maximum(((x[:parameter])[1])[1], ((y[:parameter])[1])[1])
        f = x-> map(i->Sum(Filtered(x, z->z<i)) + i*count(z->z>=i,x) ,1:m)
        fx = f(x[:parameter][1])
        fy = f(y[:parameter][1])
        for i in 1:m
          if fx[i] < fy[i] return false
          elseif fx[i] == fy[i] && i in (y[:parameter])[2]
            if i in Difference(x[:parameter][1], x[:parameter][2]) return false end
            if i < m && mod(fx[i + 1] - fy[i + 1], 2) == 1 return false end
          end
        end
        return true
      end, uc[:classes])
    end), uc[:classes])))
  if char != 2 && ctype == 2
    LuSpin=function(p)
      sort!(p)
      a = []
      b = []
      d = [0, 1, 0, -1]
      d = d[map((x->begin 1 + mod(x, 4) end), p)]
      i = 1
      while i <= length(p)
          l = p[i]
          t = Sum(d[1:i - 1])
          if 1 == mod(l, 4)
              push!(a, div(l - 1, 4) - t)
              i = i + 1
          elseif 3 == mod(l, 4)
              push!(b, div(l - 3, 4) + t)
              i = i + 1
          else
              j = i
              while i <= length(p) && p[i] == l i = i + 1 end
              j = fill(0, max(0, (1 + div(i - j, 2)) - 1))
              a = Append(a, (j + div(l + mod(l, 4), 4)) - t)
              b = Append(b, j + div(l - mod(l, 4), 4) + t)
          end
      end
      a = Filtered(a, (x->begin x != 0 end))
      a = Vector{Int}(reverse(a))
      b = Filtered(b, (x->begin x != 0 end))
      b = Vector{Int}(reverse(b))
      if Sum(d) >= 1 return [a, b]
      else return [b, a]
      end
    end
    addSpringer = function (f, i, s, k)
      ss = First(uc[:springerSeries], f)
      if s in [[[], [1]], [[], []]] p = 1
      elseif s == [[1], []] p = 2
      else p = Position(CharParams(ss[:relgroup]), [s])
      end
      ss[:locsys][p] = [i, k]
    end
    trspringer = function (i, old, new)
        for ss in uc[:springerSeries]
            for c in ss[:locsys]
                if c[1] == i
                    p = Position(old, c[2])
                    if p != false c[2] = new[p] end
                end
            end
        end
    end
    d = 0
    while 4d^2-3d<=r
      i=4d^2-3d
      if mod(r-d,2)==0
          l = Concatenation(1:i, i + 2:(i + 4) - (i + 2):r)
          ss=Dict{Symbol, Any}(:relgroup=>coxgroup(:B,div(r-i,2)),
                               :levi => l, :Z => [-1])
          ss[:locsys]=fill([0,0],NrConjugacyClasses(ss[:relgroup]))
          push!(uc[:springerSeries],ss)
          i = 4 * d ^ 2 + 3d
          if i <= r && d != 0
            l = vcat(1:i,i+2:2:r)
            ss= Dict{Symbol, Any}(:relgroup=>coxgroup(:B,div(r-i,2)),
                                  :levi => l, :Z => [-1])
            ss[:locsys]=fill([0,0],NrConjugacyClasses(ss[:relgroup]))
            push!(uc[:springerSeries], ss)
          end
      end
      d+=1
    end
    l = Filtered(eachindex(uc[:classes]), i->
     ForAll(Collected(uc[:classes][i][:parameter]), c->
                                mod(c[1],2)==0 || c[2]==1) )
    for i = l
      cl = (uc[:classes])[i]
      s = LuSpin(cl[:parameter])
      if length(cl[:Au]) == 1
          cl[:Au] = CoxeterGroup("A", 1)
          trspringer(i, [1], [2])
          d = 1
      elseif length(cl[:Au]) == 4
          cl[:Au] = CoxeterGroup("B", 2)
          trspringer(i, [1, 2, 3, 4], [1, 3, 5, 4])
          d = 2
      else
        error("Au non-commutative of order ",Size(cl[:Au])*2,"  !  implemented")
      end
      addSpringer(ss->ss[:Z]==[-1] && rank(ss[:relgroup])==sum(sum,s),i,s,d)
    end
  end
  return uc
end

chevie[:D][:UnipotentClasses]=function(n,char)
  addSpringer = function (s, i, cc)
     ss = First(uc[:springerSeries], x->x[:defect] == DefectSymbol(s[:symbol]))
     if s[:sp] in [[[], [1]], [[], []]] p = 1
     elseif s[:sp] == [[1], []] p = 2
     else p = Position(CharParams(ss[:relgroup]), [s[:sp]])
     end
     ss[:locsys][p] = [i, Position(CharParams(cc[:Au]), 
                                   map(x->x ? [1,1] : [2], s[:Au]))]
  end
  function partition2DR(part)
    p=sort(vcat(map(x->1-x:2:x-1, part)...))
    p=p[1+div(length(p),2):end]
    vcat([p[1]+p[2]], map(i->p[i+1]-p[i], 1:length(p)-1))
  end
  if char==2
    ss = XSP(4, 0, n, 1)
    symbol2partition = function (S)
      c=sort(vcat(S...))
      i = 1
      part = Int[]
      ex = Int[]
      while i <= length(c)
          l = 2 * (c[i] - 2 * (i - 1))
          if i == length(c) || c[i + 1] - c[i] > 1
              push!(part, l + 2)
              i+=1
          elseif c[i + 1] - c[i] > 0
              part = Append(part, [l+1, l+1])
              i+=2
          else
              part = Append(part, [l, l])
              i+=2
              push!(ex, l)
          end
      end
      part = filter(y->!iszero(y),sort(part)) 
      [reverse(part), ex]
    end
  else
    ss = XSP(2, 0, n, 1)
    function symbol2partition(S)
      c=sort(vcat(S...))
      i = 1
      part = Int[]
      while i <= length(c)
        l = 2 * (c[i] - (i - 1))
        if i == length(c) || c[i + 1] - c[i] > 0
          push!(part, l+1)
          i+=1
        else
          part = append!(part, [l, l])
          i+=2
        end
      end
      part = filter(y->!iszero(y),sort(part)) 
      reverse(part)
    end
  end
  l = union(map(c->map(x->
        [DefectSymbol(x[:symbol]), Sum(FullSymbol(x[:sp]), Sum)], c), ss))
  SortBy(l, x->[AbsInt(x[1]), -SignInt(x[1])])
  uc = Dict{Symbol, Any}(:classes => [], :springerSeries => map(function(d)
      res = Dict{Symbol, Any}(:defect=>d[1], :levi=>1:n - d[2])
      if mod(n - d[2], 4) == 0 || char == 2
          res[:Z]= mod(n, 2)==0  ? [1, 1] : [1]
      else
          res[:Z]= mod(n, 2)==0  ? [-1, -1] : [-1]
      end
      if d[1]==0 res[:relgroup]=coxgroup(:D, d[2])
      else res[:relgroup] = coxgroup(:B, d[2])
      end
      res[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(res[:relgroup])]
      res
  end, l))
  for cl in ss
    cc = Dict{Symbol, Any}(:parameter => symbol2partition(cl[1][:symbol]))
    if char == 2
      cc[:dimBu] = (cl[1])[:dimBu]
      cc[:name] = Join(map(function(x)
                            res=joindigits(fill(x[1], max(0, x[2])), "[]")
             if x[1] in cc[:parameter][2]
                 return SPrint("(", res, ")")
             end
             return res
         end, reverse(Collected(cc[:parameter][1]))), "")
    else
      cc[:dynkin] = partition2DR(cc[:parameter])
      cc[:name] = joindigits(cc[:parameter])
    end
    cc[:Au] = isempty(cl[1][:Au]) ? coxgroup() : 
       prod(coxgroup(:A,1) for i in eachindex(cl[1][:Au]))
    CharNames(cc[:Au])
    if char != 2
     cc[:red] = coxgroup()
     j = cc[:parameter]
     for j = Collected(j)
       if mod(j[1], 2) == 0
         cc[:red]*=coxgroup(:C, div(j[2],2))
       elseif mod(j[2], 2) != 0
         if j[2]>1 cc[:red]*=coxgroup(:B, div(j[2]-1,2)) end
       elseif j[2]>2 cc[:red]*=coxgroup(:D, div(j[2],2))
       else cc[:red]*=torus(1)
       end
     end
   end
   if !IsList(cl[1][:sp][2]) cl[1][:sp][3]=1-mod(div(n,2),2) end
   push!(uc[:classes], cc)
   for s in cl addSpringer(s, length(uc[:classes]), cc) end
   if !IsList(cl[1][:sp][2])
      cl[1][:sp][3]=1-cl[1][:sp][3]
      cc[:name]*="+"
      cc = Copy(cc)
      cc[:name]=replace(cc[:name],r".$"=>"-")
      if haskey(cc, :dynkin) cc[:dynkin][[1, 2]] = cc[:dynkin][[2, 1]] end
      push!(uc[:classes], cc)
      for s = cl addSpringer(s, length(uc[:classes]), cc) end
    end
  end
  if char == 2
    uc[:orderClasses] = Hasse(Poset(map((x->begin
    map(function (y,)
      m = maximum(((x[:parameter])[1])[1], ((y[:parameter])[1])[1])
      f = x-> map(i-> Sum(Filtered(x, z->z<i))+i*count(z->z>=i, x) , 1:m)
      fx = f(x[:parameter][1])
      fy = f(y[:parameter][1])
      for i = 1:m
          if fx[i] < fy[i] return false
          elseif fx[i] == fy[i] && i in (y[:parameter])[2]
              if i in Difference(x[:parameter][1], x[:parameter][2])
                  return false
              end
              if i<m && mod(fx[i+1]-fy[i+1],2)==1 return false end
          end
      end
      if x[:parameter] == y[:parameter] && x != y return false end
      return true
    end, uc[:classes])
    end), uc[:classes])))
  else
    uc[:orderClasses] = Hasse(Poset(map((i->begin map((j->begin
      Dominates(uc[:classes][j][:parameter], uc[:classes][i][:parameter]) && 
      (uc[:classes][j][:parameter]!=uc[:classes][i][:parameter] || i == j)
                                 end), 1:length(uc[:classes]))
                     end), 1:length(uc[:classes]))))
  end
  if char != 2
    d = 0
    while 4*d^2-d <= n
     i = 4*d^2-d
     if mod(n-d,2)==0
       l=Concatenation(1:i,i+2:2:n)
       s=Dict(:relgroup=>coxgroup(:B, div(n-i,2)),:levi=>l)
       if mod(n, 2) == 0 s[:Z]=[1, -1] else s[:Z] = [E(4)] end
       s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
       push!(uc[:springerSeries], s)
       if d==0 l=Concatenation([1], 4:2:n) end
       s=Dict(:relgroup => coxgroup(:B, div(n-i,2)),:levi=>l)
       if mod(n,2)==0 s[:Z] = [-1, 1]
       else s[:Z] = [-(E(4))]
       end
       s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
       push!(uc[:springerSeries], s)
       i=4d^2+d
       if d != 0 && i <= n
         l = vcat(1:i, i+2:2:n)
         s = Dict(:relgroup=>coxgroup(:B, div(n - i,2)),:levi=>l)
         if mod(n,2)==0 s[:Z]=[1, -1] else s[:Z]=[E(4)] end
         s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
         push!(uc[:springerSeries], s)
         s=Dict(:relgroup=>oxgroup("B",div(n-i,2)),:levi=>l)
         if mod(n, 2) == 0 s[:Z] = [1, 1] else s[:Z] = [-(E(4))] end
         s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
         push!(uc[:springerSeries], s)
       end
     end
     d+=1
  end
  function LuSpin(p)
    sort!(p)
    a = Int[]
    b = Int[]
    d = [0, 1, 0, -1]
    d = d[map((x->begin 1 + mod(x, 4) end), p)]
    i = 1
    while i <= length(p)
        l = p[i]
        t = Sum(d[1:i - 1])
        if 1 == mod(l, 4)
          push!(a, div(l - 1,4) - t)
          i+=1
        elseif 3 == mod(l, 4)
          push!(b, div(l - 3,4) + t)
          i+=1
        else
            j = i
            while i <= length(p) && p[i]==l i+=1 end
            j = fill(0, max(0, div(i - j,2)) )
            a = Append(a, j + div(l + mod(l, 4), 4) - t)
            b = Append(b, j + div(l - mod(l, 4), 4) + t)
        end
    end
    a = Filtered(a, (x->begin x != 0 end))
    a = reverse(a)
    b = Filtered(b, (x->begin x != 0 end))
    b = reverse(b)
    if Sum(d) >= 1 return [a, b]
    else return [b, a]
    end
  end
  addSpringer1= function (f, i, s, k)
    ss = First(uc[:springerSeries], f)
    if s in [[[], [1]], [[], []]] p = 1
    elseif s == [[1], []] p = 2
    else p = Position(CharParams(ss[:relgroup]), [s])
    end
    ss[:locsys][p] = [i, k]
  end
  trspringer = function (i, new)
    for ss in uc[:springerSeries]
      for c in ss[:locsys] if c[1] == i c[2] = new[c[2]] end end
    end
  end
  l=Filtered(1:length(uc[:classes]), (i->begin
                 ForAll(Collected(((uc[:classes])[i])[:parameter]), c->
                             mod(c[1], 2) == 0 || c[2] == 1) end))
  for i in l
     cl = (uc[:classes])[i]
     s = LuSpin(cl[:parameter])
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
                     rank(ss[:relgroup]) == Sum(s, Sum) , i, s, k[1])
     end
     if !('+' in cl[:name])
         addSpringer1(ss->ss[:Z] in [[-1, 1], [-(E(4))]] && 
                     rank(ss[:relgroup]) == Sum(s, Sum), i, s, k[2])
     end
  end
end
  uc
end

end
