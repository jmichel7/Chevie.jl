module HasType

export charinfo, classinfo, reflection_name, diagram, chartable,
  representation, fakedegrees, unipotent_characters, 
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

charname(W,x;TeX=false)=join(map((t,p)->getchev(t,:CharName,p,
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
 t=TypeIrred(Dict(:series=>:ST,:p=>p,:q=>q,:rank=>r))
  r=getchev(t,:GeneratingRoots)
  cr=getchev(t,:EigenvaluesGeneratingReflections)
  cr=map((x,y)->coroot(x,y),r,map(x->E(Root1(numerator(x),denominator(x))),cr))
  cr=map(x->convert.(Cyc{Rational{Int}},x),cr)
  PRG(r,cr)
end

degrees(W)=vcat(getchev(W,:ReflectionDegrees)...)

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

nr_conjugacy_classes(W)=prod(getchev(W,:NrConjugacyClasses))

PrintToSring(s,v...)=sprint(show,v...)

reflection_name(W)=join(getchev(W,:ReflectionName,Dict()),"×")

function representation(W,i::Int)
  map(x->hcat(x...),impl1(getchev(W,:Representation,i)))
end

UnipotentClassOps=Dict(:Name=>x->x)

unipotent_classes(W,p=0)=impl1(getchev(W,:UnipotentClasses,p))

include("uch.jl")
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

ApplyFunc(f,x)=f(x...)
Arrangements=arrangements
BetaSet=βSet
Binomial=binomial
CartanMat(s,a...)=cartan(Symbol(s),a...)
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
Drop(a::Vector,i::Int)=deleteat!(copy(a),i)
Factors(n)=vcat([fill(k,v) for (k,v) in factor(n)]...)
Filtered(l,f)=isempty(l) ? l : filter(f,l)
First(a,b)=a[findfirst(b,a)]
Flat(v)=collect(Iterators.flatten(v))
ForAll(l,f)=all(f,l)
FullSymbol=fullsymbol
HighestPowerFakeDegreeSymbol=degree_feg_symbol
HighestPowerGenericDegreeSymbol=degree_gendeg_symbol
IdentityMat(n)=map(i->one(rand(Int,n,n))[i,:],1:n)
function Ignore() end
InfoChevie2=print
IntListToString(l)=any(x->x>10,l) ? join(l,",") : join(l)
IsList(l)=l isa Vector
IsInt(l)=l isa Int ||(l isa Rational && denominator(l)==1)
Lcm(a...)=Lcm(collect(a))
Lcm(a::Vector)=lcm(Int.(a))
LowestPowerFakeDegreeSymbol=valuation_feg_symbol
LowestPowerGenericDegreeSymbol=valuation_gendeg_symbol
Minimum(v::AbstractVector)=minimum(v)
Minimum(a::Number,x...)=min(a,x...)
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
Product(v)=isempty(v) ? 1 : prod(v)
Product(v,f)=isempty(v) ? 1 : prod(f,v)
RankSymbol=ranksymbol
RecFields=keys
RootInt(a,b)=floor(Int,a^(1/b))
Rotations(a)=circshift.(Ref(a),0:length(a)-1)
gapSet(v)=unique(sort(v))
SemisimpleRank(W)=length(independent_roots(W))
Rank=rank
ShiftBeta=shiftβ
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
Value(p,v)=p(v)
CoxeterGroup(s::String,n)=coxgroup(Symbol(s),Int(n))
CoxeterGroup()=coxgroup()

ExtendedReflectionGroup(W,mats)=1

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
  term=v->cyclotomic_polynomial(v[2])(monomial(v[1]))
  r = arg[2]
  if haskey(r, :rootUnity) && haskey(r,:root) eror("cannot have both") end
  if haskey(r, :coeff) res = r[:coeff] else res = 1 end
  if haskey(r, :factor) res*=monomial(r[:factor]) end
  if haskey(r, :rootUnity)
    term=v->cyclotomic_polynomial(v[2])(prod(map(^,para,v[1][1:n]))*
        r[:rootUnity]^data[:rootUnityPower])
  end
  if haskey(r, :root)
    term=function(v)
      mon,pol=v
     if length(mon)==n return Pol([cyclotomic_polynomial(pol)(monomial(mon))],0)
     else return cyclotomic_polynomial(pol)(Pol([monomial(mon)],mon[n+1]))
     end
    end
    res*=prod(term,r[:vcyc])
    den=lcm(denominator.(r[:root])...)
    root=monomial(den*r[:root])
    if haskey(r, :rootCoeff) root*=r[:rootCoeff] end
    EvalPolRoot(res, root, den, data[:rootPower])
  else res*=prod(term,r[:vcyc])
  end
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
  s=IntListToString(map(length,S))
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
Torus(i::Int)=i
RootsCartan=x->x

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
  ct[:irredinfo]=map(x->Dict(:charname=>IntListToString(x)),chevie[:A][:CharInfo](n)[:charparams])
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
end
