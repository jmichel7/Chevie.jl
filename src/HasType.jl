module HasType

export reflection_name, diagram, UnipotentCharacters, 
  UnipotentClasses, schur_elements, charname, codegrees, ComplexReflectionGroup,
  chevieget, field, getchev, weightinfo, Cartesian

using ..Gapjm
#-----------------------------------------------------------------------
const chevie=Dict()

# extensions to get closer to GAP semantics
Base.:*(a::Array,b::Pol)=a .* Ref(b)
Base.:*(a::Pol,b::Array)=Ref(a) .* b
Base.:*(a::AbstractVector{<:Number},b::AbstractVector{<:Number})=sum(a.*b)
Base.:*(a::AbstractVector{Pol},b::AbstractVector{Pol})=sum(a.*b)
Base.:*(a::AbstractVector,b::AbstractVector)=toL(toM(a)*toM(b))
Base.:-(a::AbstractVector,b::Number)=a .- b
Base.:+(a::Integer,b::AbstractVector)=a .+ b
Base.:+(a::AbstractArray,b::Number)=a .+ b
Base.:+(a::AbstractArray,b::Pol)=a .+ Ref(b)
#Base.:*(a::Mvp, b::Array)=Ref(a).*b
#include("mvp.jl")
#Base.:*(b::Array,a::Mvp)=b.*Ref(a)
Base.getindex(s::String,a::Vector{Any})=getindex(s,Int.(a))
Cycs.:^(a::Cyc,b::Rational)=a^Int(b)
Base.:^(m::AbstractMatrix,n::AbstractMatrix)=inv(n*E(1))*m*n
Base.:^(m::Vector,n::Vector)=toL(inv(toM(n)*E(1))*toM(m)*toM(n))
Base.:(//)(m::Vector,n::Vector)=toL(toM(m)*inv(toM(n)*E(1)))
Base.:^(m::Vector{<:Vector{<:Number}},n::Matrix{<:Number})=inv(n)*toM(m)*n
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
#     println("args=",(d[2:end]...,extra...,t[:cartantype]))
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

impl1(l)=length(l)==1 ? l[1] : error("implemented only for irreducible groups")

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

charname(W,x;TeX=false,opt...)=join(map((t,p)->getchev(t,:CharName,p,
                           TeX ? Dict(:TeX=>true) : Dict()),refltype(W),x),",")

function PositionCartesian(l,ind)
  res=prod=1
  for i in length(l):-1:1
    res+=(ind[i]-1)*prod
    prod*=l[i]
  end
  res
end

function PermGroups.class_reps(W::PermRootGroup)
  gets(W,:classreps)do W
    map(x->W(x...),classinfo(W)[:classtext])
  end
end
PermGroups.class_reps(W::FiniteCoxeterGroup)=class_reps(W.G)


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
    n=one(toM(m))
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

function WeightToAdjointFundamentalGroupElement(W,i)
  t=First(refltype(W),t->i in t[:indices])
  l=copy(t[:indices])
  b=longest(W,l)*longest(W,setdiff(l,[i]))
  push!(l,maximum(findall(
    i->all(j->j in t[:indices] || W.rootdec[i][j]==0,1:semisimplerank(W)),
  eachindex(W.rootdec))))
  restricted(b,inclusion.(Ref(W),l))
end

# returns a record containing minuscule coweights, decompositions
# (in terms of generators of the fundamental group)
function weightinfo(W)
  l=map(refltype(W)) do t
    r=getchev(t,:WeightInfo)
    if isnothing(r)
      r=Dict{Symbol,Any}(:moduli=>Int[],:decompositions=>Vector{Vector{Int}}[],
           :minusculeWeights=>Vector{Int}[])
    end
    if !haskey(r,:minusculeCoweights)
      r[:minusculeCoweights]=r[:minusculeWeights]
    end
    if isempty(r[:moduli]) g=Int[]
      r[:ww]=Perm{Int}[]
    else g=filter(i->sum(r[:decompositions][i])==1,
          eachindex(r[:minusculeCoweights])) # generators of fundamental group
      r[:ww]=map(x->WeightToAdjointFundamentalGroupElement(W,x),
               t[:indices][r[:minusculeCoweights][g]])
    end
    r[:csi]=zeros(Rational{Int},length(g),semisimplerank(W))
    if !isempty(r[:moduli]) 
      C=mod1.(inv(Rational.(cartan(t.prop))))
      r[:csi][:,t[:indices]]=C[r[:minusculeCoweights][g],:]
      r[:minusculeWeights]=t[:indices][r[:minusculeWeights]]
      r[:minusculeCoweights]=t[:indices][r[:minusculeCoweights]]
    end
    r[:csi]=toL(r[:csi])
    r
  end
  res=Dict(:minusculeWeights=>Cartesian(map(
                                        x->vcat(x[:minusculeWeights],[0]),l)),
    :minusculeCoweights=>Cartesian(map(
                                     x->vcat(x[:minusculeCoweights],[0]),l)),
    :decompositions=>map(vcat,Cartesian(map(x->vcat(x[:decompositions],
                                 [0 .*x[:moduli]]),l))),
    :moduli=>vcat(map(x->x[:moduli],l)...))
# centre of simply connected group: the generating minuscule coweights
# mod the root lattice
  res[:CenterSimplyConnected]=vcat(getindex.(l,:csi)...)
  res[:AdjointFundamentalGroup]=vcat(getindex.(l,:ww)...)
  n=length(res[:decompositions])-1
  res[:minusculeWeights]=map(x->filter(y->y!=0,x),res[:minusculeWeights][1:n])
  res[:minusculeCoweights]=map(x->filter(y->y!=0,x),res[:minusculeCoweights][1:n])
  res[:decompositions]=res[:decompositions][1:n]
  res
end

nr_conjugacy_classes(W)=prod(getchev(W,:NrConjugacyClasses))

PrintToSring(s,v...)=sprint(show,v...)

reflection_name(W,opt=Dict())=join(getchev(W,:ReflectionName,opt),"×")

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
function DiagonalMat(v...)
  R=cat(map(m->m isa Array ? m : hcat(m),v)...,dims=(1,2))
  for i in axes(R,1), j in axes(R,2)
    if i!=j R[i,j]=zero(R[1,1]) end
  end
  R
end
DiagonalMat(v::Vector{<:Number})=DiagonalMat(v...)
DiagonalOfMat(m)=[m[i,i] for i in axes(m,1)]
DivisorsInt=divisors
Dominates=dominates
Drop(a::Vector,i::Int)=deleteat!(copy(a),i)
EltWord(W,x)=W(x...)
ExteriorPower(m,i)=toL(exterior_power(toM(m),i))
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
NullMat(i,j=i)=[zeros(Int,j) for k in 1:i]
function Ignore() end
InfoChevie2=print
IntListToString=joindigits
IsList(l)=l isa Vector
IsInt(l)=l isa Int ||(l isa Rational && denominator(l)==1)
IsString(l)=l isa String
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
  phis::Vector{<:Perm}
end

function ExtendedCox(W::FiniteCoxeterGroup,F0s::Vector{Matrix{Int}})
  phis=map(F->Perm(F,roots(parent(W.G)),action=(r,m)->permutedims(m)*r),F0s)
  ExtendedCox(W,F0s,phis)
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
  mats=toM.(mats)
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
Group(a::Perm...)=Group(collect(a))
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
include("table2.jl")

end
