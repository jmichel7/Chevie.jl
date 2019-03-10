module HasType

export charinfo, classinfo, reflection_name, diagram, chartable,
  representation, fakedegrees, unipotent_characters, 
  schur_elements, charname, codegrees, ComplexReflectionGroup,
  getclassified, chevieget, field, getfromtype, Cartesian

using Gapjm

Cartesian(a::AbstractVector...)=reverse.(vec(collect.(Iterators.product(reverse(a)...))))

braid_relations(W)=impl1(getclassified(W,:BraidRelations))

function codegrees(W)
  vcat(map(refltype(W)) do t
    cd=getfromtype(t,:ReflectionCoDegrees)
    if cd===nothing
      cd=getfromtype(t,:ReflectionDegrees)
      maximum(cd).-cd
    else cd
    end
  end...)
end

impl1(l)=length(l)==1 ? l[1] : error("implemented only for irreducible groups")

charname(W,x;TeX=false)=join(map((t,p)->getfromtype(t,:CharName,p,
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

function charinfo(W)
  p=map(refltype(W)) do t
    c=copy(getfromtype(t,:CharInfo))
    c[:positionId]=c[:extRefl][1]
    c[:positionDet]=c[:extRefl][end]
    c
  end
  if length(p)==1 res=copy(p[1]) else res=Dict{Symbol, Any}() end
  res[:charparams]=cartfields(p,:charparams)
  res[:charnames]=charname.(Ref(W),res[:charparams])
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
  return res
end

function classinfo(W)
  tmp = getclassified(W,:ClassInfo)
  if isempty(tmp) return Dict(:classtext=>[Int[]],:classnames=>[""],
                    :classparams=>[Int[]],:orders=>[1],:classes=>[1])
  end
  if any(isnothing, tmp) return nothing end
  if length(tmp)==1 res=copy(tmp[1]) else res=Dict{Symbol, Any}() end
  res[:classtext]=map(x->vcat(x...),Cartesian(map((i,t)->
                                             getindex.(Ref(i),t[:classtext]),
                        map(x->x[:indices],refltype(W)),tmp)...))
  res[:classnames]=map(join,cartfields(tmp,:classnames))
  if allhaskey(tmp, :classparam)
    res[:classparams]=cartfields(tmp,:classparams)
  end
  if allhaskey(tmp,:orders)
    res[:orders]=map(lcm, cartfields(tmp,:orders))
  end
  if allhaskey(tmp,:classes)
    res[:classes]=Int.(map(prod, cartfields(tmp,:classes)))
  end
  return res
end

function chartable(ct::Dict)
  CharTable(permutedims(Cyc{Int}.(hcat(ct[:irreducibles]...))),
     map(x->x[:charname],ct[:irredinfo]),
        ct[:classnames],map(Int,ct[:centralizers]),ct[:identifier])
end

function chartable(W)
  ctt=map(refltype(W)) do t
    ct=getfromtype(t,:CharTable)
    if haskey(ct,:irredinfo)
      names=getindex.(ct[:irredinfo],:charname)
    else
      names=map(p->getfromtype(t,:CharName,p,Dict(:TeX=>true)),
                getfromtype(t,:CharInfo)[:charparams])
    end
    if !haskey(ct,:classnames)
      merge!(ct,getfromtype(t,:ClassInfo))
    end
    CharTable(permutedims(Cyc{Int}.(hcat(ct[:irreducibles]...))),names,
    ct[:classnames],Int.(ct[:centralizers]),ct[:identifier])
  end
  charnames=join.(Cartesian(getfield.(ctt,:charnames)...),",")
  classnames=join.(Cartesian(getfield.(ctt,:classnames)...),",")
  centralizers=prod.(Cartesian(getfield.(ctt,:centralizers)...))
  identifier=join(getfield.(ctt,:identifier),"×")
  irr=length(ctt)==1 ? ctt[1].irr : kron(getfield.(ctt,:irr)...)
  CharTable(irr,charnames,classnames,centralizers,identifier)
end

function chartable(H::HeckeAlgebra{C})where C
  W=H.W
  ct=impl1(getclassified(W,:HeckeCharTable,H.para,H.sqpara))
  CharTable(Matrix(permutedims(hcat(
                map(ch->convert.(C,ch),ct[:irreducibles])...))),
     map(x->getclassified(W,:CharName,x,Dict(:TeX=>true)),
         charinfo(W)[:charparams]),
     ct[:classnames],map(Int,ct[:centralizers]),ct[:identifier])
end

function ComplexReflectionGroup(i::Int)
  t=Dict(:series=>:ST,:ST=>i)
  r=getfromtype(t,:GeneratingRoots)
  cr=getfromtype(t,:GeneratingCoRoots)
  if cr===nothing
    e=getfromtype(t,:EigenvaluesGeneratingReflections)
    cr=map((x,y)->coroot(x,y),r,E.(Root1.(e)))
  end
  PermRootGroup(r,cr)
end

function ComplexReflectionGroup(p,q,r)
  t=Dict(:series=>:ST,:p=>p,:q=>q,:rank=>r)
  r=getfromtype(t,:GeneratingRoots)
  cr=getfromtype(t,:EigenvaluesGeneratingReflections)
  cr=map((x,y)->coroot(x,y),r,map(x->E(Root1(x)),cr))
  PermRootGroup(r,cr)
end

degrees(W)=vcat(getclassified(W,:ReflectionDegrees)...)

function diagram(W)
  for t in refltype(W)
    getfromtype(t,:PrintDiagram,t[:indices],
               getfromtype(t,:ReflectionName,Dict()))
  end
end

function fakedegree(W,p,q)
  prod(map((t,p)->getfromtype(t,:FakeDegree,p,q),refltype(W),p))
end

function fakedegrees(W,q)
  map(p->fakedegree(W,p,q),charinfo(W)[:charparams])
end

nr_conjugacy_classes(W)=prod(getclassified(W,:NrConjugacyClasses))

reflection_name(W)=join(getclassified(W,:ReflectionName,Dict()),"×")

function representation(W,i::Int)
  map(x->hcat(x...),impl1(getclassified(W,:Representation,i)))
end

function representation(H::HeckeAlgebra,i::Int)
  ct=impl1(getclassified(H.W,:HeckeRepresentation,H.para,H.sqpara,i))
  map(x->hcat(x...),ct)
end

function schur_elements(H::HeckeAlgebra)
  W=H.W
  map(p->getclassified(W,:SchurElement,p,H.para,H.sqpara),
    charinfo(W)[:charparams])
end

unipotent_characters(W)=impl1(getclassified(W,:UnipotentCharacters))

unipotent_classes(W,p=0)=impl1(getclassified(W,:UnipotentClasses,p))

const chevie=Dict()

Base.:*(a::Array,b::Pol)=a .* Ref(b)
Base.:*(a::Pol,b::Array)=Ref(a) .* b
Base.:*(a::AbstractVector,b::AbstractVector)=sum(a.*b)
Base.:-(a::AbstractVector,b::Int)=a .- b
Base.:+(a::Integer,b::AbstractVector)=a .+ b
Base.:+(a::AbstractVector,b::Number)=a .+ b
#Base.:*(a::Mvp, b::Array)=Ref(a).*b
#include("mvp.jl")
#Base.:*(b::Array,a::Mvp)=b.*Ref(a)
Base.getindex(s::String,a::Vector{Any})=getindex.(s,a)
Cycs.:^(a::Cyc,b::Rational)=a^Int(b)

function chevieget(t::Symbol,w::Symbol)
 if haskey(chevie[t],w) return chevie[t][w] end
 println("chevie[$t] has no $w")
end

Cycs.E(a,b=1)=Cycs.E(Int(a),Int(b))

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

function field(t)
  s=t[:series]
  if s in [:A,:B,:D] return (s,length(t[:indices]))
  elseif s==:ST 
    if haskey(t,:ST)
      if 4<=t[:ST]<=22 return (:G4_22,t[:ST])
      else return (Symbol(string("G",t[:ST])),)
      end
    else
      return (:imp, t[:p], t[:q], t[:rank])
    end
  else return (Symbol(string(s,length(t[:indices]))),) 
  end
end

const needcartantype=Set([:PrintDiagram,:ReflectionName,:UnipotentClasses])

function getfromtype(t,f::Symbol,extra...)
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

function getclassified(W,f::Symbol,extra...)
  [getfromtype(ti,f::Symbol,extra...) for ti in refltype(W)]
end

#----------------------------------------------------------------------
# correct translations of GAP3 functions

function pad(s::String, i::Int)
  if i>0 return lpad(s,i)
  else return rpad(s,-i)
  end
end

pad(s::String)=s

function Replace(s,p...)
  p=Dict(p[i]=>p[i+1] for i in 1:2:length(p))
  p=[haskey(p,s[i:i]) ? p[s[i:i]] : s[i:i] for i in 1:length(s)]
  vcat(p...)
end

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

Sublist(a::Vector, b::AbstractVector)=a[b]

Append(a::Vector,b::AbstractVector)=vcat(a,b)
Append(a::String,b::String)=a*b
Append(a::String,b::Vector{Char})=a*String(b)

Drop(a::Vector,i::Int)=deleteat!(copy(a),i)

Base.getindex(a::Symbol,i::Int)=string(a)[i]
Base.length(a::Symbol)=length(string(a))

IsList(l)=l isa Vector
IsInt(l)=l isa Int ||(l isa Rational && denominator(l)==1)

ForAll(l,f)=all(f,l)

PartitionTuples=partition_tuples

function PartitionTupleToString(n,a=Dict())
  if n[end] isa Vector return join(map(join,n),".") end
  r=repr(E(n[end-1],n[end]),context=:limit=>true)
  if r=="1" r="+" end
  if r=="-1" r="-" end
  join(map(join,n[1:end-2]),".")*r
end

Product(v)=isempty(v) ? 1 : prod(v)
Product(v,f)=isempty(v) ? 1 : prod(f,v)

IntListToString(l)=any(x->x>10,l) ? join(l,",") : join(l)

Lcm(a...)=Lcm(collect(a))
Lcm(a::Vector)=lcm(Int.(a))

function Collected(v)
  d=groupby(v,v)
  sort([[k,length(v)] for (k,v) in d])
end

SortBy(x,f)=sort(x,by=f)

function Ignore() end

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

Base.copy(x::Char)=x

function ShiftBeta(beta,n)
  if n>=0 return [0:n-1;beta .+ n]
  elseif beta[1:-n]!=[0:-n-1]
        Error( "Cannot shift ", beta, " by ", n, "\n" );
  else return beta[1-n:end].+n
  end
end

BetaSet(alpha)=isempty(alpha) ? alpha : reverse(alpha) .+(0:length(alpha)-1)

function SymbolPartitionTuple(p,s)
  if IsInt(p[end])
    l= length(p) - 2
    e= l*p[end-1]
  else e=l=length(p)
  end
  if IsInt(s)
    if s<0  s=[0;fill(-s,l-1)]
    else    s=[s;zeros(Int,l-1)]
    end
  else s=s[1:l]
  end
  s= map(length, p[1:l]) .- s
  s= maximum(s) .- s
  p= copy(p)
  p[1:l]=map(i->ShiftBeta(BetaSet(p[i]),s[i]),1:l)
  p
end

function LowestPowerGenericDegreeSymbol(p)
  p=FullSymbol(p)
  e=length(p)
  p=sort!(vcat(p...))
  m=length(p)
  sum(p.*(m-1:-1:0))-div(m*(m-e)*(2*m-3-e),12*e)
end

function HighestPowerGenericDegreeSymbol(p)
  p=FullSymbol(p)
  r=RankSymbol(p)
  e=length(p)
  p=sort!(vcat(p...))
  m=length(p)
  if mod(m,e)==1 r=div(e*r*(r+1),2)
  else           r+= div(e*r*(r-1),2)
  end
  r+sum(p.*(0:m-1))-sum(x->div(e*x*(x+1),2),p)-div(m*(m-e)*(2*m-3-e),12*e)
end

function DefectSymbol(s)
  s=FullSymbol(s)
  length(s[1])-length(s[2])
end

function HighestPowerFakeDegreeSymbol(s,p=length(FullSymbol(s)))
  s=FullSymbol(s)
  e=length(s)
  if e==1 d=0
  else d=DefectSymbol(s)
  end
  if !(d in [0,1]) return -1 end
  r=RankSymbol(s)
  if d==1 res=div(e*r*(r+1),2)
  else    res=div(e*r*(r-1),2)+r
  end
  res+=e*sum(S->sum(S.*(0:length(S)-1))-sum(l->div(l*(l+1),2),S),
               filter(x->!isempty(x),s))
  gamma=i->sum(mod.(i+(0:e-1),e).*map(sum,s))
  if d==1 res+=gamma(0)
  else res+=maximum(gamma.(0:e-1))
  end
  res-sum(map(x->div(x*(x-1),2),e*(1:div(sum(length,s),e)-1).+mod(d,2))) 
end

function LowestPowerFakeDegreeSymbol(s)
  s=FullSymbol(s)
  if length(s)==1 d=0
  else d=DefectSymbol(s)
  end
  if !(d in [0,1]) return -1 end
  e=length(s)
  res=e*sum(S->sum(S.*(length(S)-1:-1:0)),filter(x->!isempty(x),s))
  gamma=i->sum(mod.(i+(0:e-1),e).*map(sum,s))
  if d==1 res+=gamma(0)
  else res+=minimum(gamma.(0:e-1))
  end
  res-sum(map(x->div(x*(x-1),2),e*(1:div(sum(length,s),e)-1).+mod(d,2))) 
end

function FullSymbol(S)
  if isempty(S) || S[end] isa AbstractVector return S end
  m=S[end-1]
  vcat(map(i->map(copy,S[1:end-2]),1:m)...)
end

function RankSymbol(s)
  s=FullSymbol(s)
  ss=sum(length,s)
  if isempty(s) return 0 end
  e=length(s)
  sum(sum,s)-div((ss-1)*(ss-e+1),2*e)
end

Arrangements=arrangements
Partitions=partitions
# e-symbols of rank r, Malle-defect def and content=c (mod e)
# The list is returned sorted by HC series (principal series first)
# SymbolsDefect(d,r,0,1) gives symbols of unipotent characters of G(d,1,r)
# SymbolsDefect(e,r,0,0) gives symbols of unipotent characters of G(e,e,r)
SymbolsDefect=function(e,r,def,c)
  local IsReducedSymbol, S
  function defShape(s)local e
    e=length(s)
    (binomial(e,2)*div(sum(s),e)-sum(i->i*s[i+1],0:e-1))%e
  end
 
  shapesSymbols=function(r,e,c)local f,res,m,new
    f=function(lim2,sum,nb,max)local res,a
      if nb==1   
        if sum==0   return [[sum]]
        else return [] end 
      end
      res=[]
      a=div(sum,nb-1)
      while a<=max  &&  binomial(a,2)<=lim2  &&  a<=sum  
        append!(res,map(x->vcat([a],x),f(lim2-binomial(a,2),sum-a,nb-1,a)))
        a+=1 
      end
      return res
    end

    res=[]
    m=0
    while true
      new=f(r+div((m*e+c-1)*(m*e+c-e+1),2*e),c+m*e,e,c+m*e)
      append!(res,new)
      m+=1
      if length(new)==0 break end
    end
    res=vcat(map(x->Arrangements(x,e),res)...)
    filter(s->defShape(s)==def  &&  
      all(x->defShape(x)!=def  ||  x<=s,map(i->circshift(s,i),1:length(s)-1)),res)
  end

  function IsReducedSymbol(s)
    ForAll(Rotations(s){[2..Length(s)]},x->s==x || LessSymbols(x,s))
  end

  println("shapesSymbols(r,e,c)=$(shapesSymbols(r,e,c))")
  S=vcat(map(shapesSymbols(r,e,c)) do s
    map(x->SymbolPartitionTuple(x,s),
       PartitionTuples(r-RankSymbol(map(x->0:x-1,s)),e)) end...)
  c!=0 ? S : filter(IsReducedSymbol,S)
end

Filtered(l,f)=isempty(l) ? l : filter(f,l)

Sum(v::AbstractVector)=sum(v)
Sum(v::AbstractVector,f)=isempty(v) ? 0 : sum(f,v)

Factors(n)=vcat([fill(k,v) for (k,v) in factor(n)]...)

RecFields=keys
Sort=sort!
InfoChevie2=print

function VcycSchurElement(arg...)
  local r, data, i, para, res, n, monomial, den, root
  n = length(arg[1])
  println(arg)
  if length(arg) == 3
      data = arg[3]
      para = arg[1][data[:order]]
  else
      para = copy(arg[1])
  end
  monomial = v->Product(1:length(v), i->para[i] ^ v[i])
  r = arg[2]
  if haskey(r, :coeff) res = r[:coeff] else res = 1 end
  if haskey(r, :factor) res = res * monomial(r[:factor]) end
  if haskey(r, :root)
      para = para + 0 * Product(para)
      para[n + 1] = ChevieIndeterminate(para)
  elseif haskey(r, :rootUnity)
      para[n + 1] = r[:rootUnity] ^ data[:rootUnityPower]
  end
  res = res * Product(r[:vcyc], x->
            Value(CyclotomicPolynomial(Cyclotomics, x[2]), monomial(x[1])))
  if haskey(r, :root)
      den = Lcm(map(denominator, r[:root]))
      root = monomial(den * r[:root])
      if haskey(r, :rootCoeff) root = root * r[:rootCoeff] end
      return EvalPolRoot(res, root, den, data[:rootPower])
  else
      return res
  end
end

function CycPols.CycPol(v::AbstractVector)
  coeff=v[1]
  valuation=convert(Int,v[2])
  vv=Pair{Rational{Int},Int}[]
  v1=convert.(Rational{Int},v[3:end])
  for i in v1
    if denominator(i)==1
       k=convert(Int,i)
       for j in prime_residues(k) push!(vv,j//k=>1) end
    else
      push!(vv,i=>1)
    end
  end
  CycPol(vv,valuation,coeff)
end

Minimum(v::AbstractVector)=minimum(v)
Minimum(a::Number,x...)=min(a,x...)
Value(p,v)=p(v)
SPrint=string

Inherit(a,b)=merge!(a,b)
function Inherit(a,b,c)
  for k in c a[Symbol(k)]=b[Symbol(k)] end
  a
end

function CharRepresentationWords(mats,words)
  mats=map(x->hcat(x...),mats)
  map(words)do w
    if isempty(w) return size(mats[1],1) end
    m=prod(mats[w])
    sum(i->m[i,i],axes(m,1))
  end
end

function CycPolFakeDegreeSymbol(arg...)
  local s, p, res, delta, theta, r, d, e, q, rot
  s = FullSymbol(arg[1])
  q = Pol([1],1)
  e = length(s)
  r = RankSymbol(s)
  if e == 0 return CycPol(1) end
  if length(arg) == 2 p = E(e, arg[2])
  else p = 1
  end
  delta=S->Product(collect(combinations(S,2)),x->CycPol(q^(e*x[2])-q^(e*x[1])))
  theta=S->Product(Filtered(S,x->x>0),l->Product(1:l,h->CycPol(q^(e*h)-1)))
  if length(s) == 1 d = 0
  else d = DefectSymbol(s)
  end
  if d == 1 res = theta([r])
  elseif d == 0 res = theta([r - 1]) * CycPol(q ^ r - p)
  else res = CycPol(0q)
  end
  res = (res * Product(s, (S-> delta(S) // theta(S)))) //
   CycPol(q ^ Sum(e * (1:div(Sum(s, length), e) - 1) + d % 2, (x->
                              (x * (x - 1)) // 2)))
  if d == 1
      res = res * CycPol(q ^ ((0:e - 1) * map(Sum, s)))
  else
      rot = Rotations(s)
      res*=CycPol(map(j->p^j,0:e-1)*
         map(s->q^((0:e-1)*map(Sum,s)), rot)) // count(i->i==s, rot)
      if e == 2 && p == -1
          res = -res
      end
  end
  if r == 2 && (e > 2 && p == E(e))
      res = CycPol(Value(res, E(2e) * q) // E(2e, Degree(res)))
  end
  return res
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
       return "G_{3,3,3}["* (s=="300" ? "\\zeta_3" : "\\zeta_3^2")*"]"
      elseif r==3 && d==4 
       return "G_{4,4,3}["* (s=="3010" ? "\\zeta_4" : "-\\zeta_4")*"]"
    else return "G_{$d,$d,$r}^{$s}"
    end
  end
end

WeylGroup(s::String,n)=WeylGroup(Symbol(s),Int(n))
#-------------------------------------------------------------------------
#  dummy translations of GAP3 functions
CHEVIE=Dict{Symbol,Any}(:compat=>Dict(:MakeCharacterTable=>x->x,
                           :AdjustHeckeCharTable=>(x,y)->x,
        :ChangeIdentifier=>function(tbl,n)tbl[:identifier]=n end))

FamilyOps=Dict()
CHEVIE[:families]=Dict(:C1=>
        Dict(:group=>"C1", :name=>"C_1", :explanation=>"trivial",
  :charLabels=>[""], :fourierMat=>[[1]], :eigenvalues=>[1],
  :mellin=>[[1]],:mellinLabels=>[""]),
  :C2=>Dict(:group=>"C2", :name=>"C_2",
  :explanation=>"DrinfeldDouble(Z/2)",
  :charLabels=>["(1,1)", "(g_2,1)", "(1,\\varepsilon)", "(g_2,\\varepsilon)"],
  :fourierMat=>1/2*[[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]],
  :eigenvalues=>[1,1,1,-1],
  :perm=>(),
  :mellin=>[[1,1,0,0],[1,-1,0,0],[0,0,1,1],[0,0,1,-1]],
  :mellinLabels=>["(1,1)","(1,g2)","(g2,1)","(g2,g2)"]),
  :S3=>Dict(:group=>"S3", :name=>"D(S_3)",
  :explanation=>"Drinfeld double of S3, Lusztig's version",
  :charLabels=>[ "(1,1)", "(g_2,1)", "(g_3,1)", "(1,\\rho)", "(1,\\varepsilon)",
		"(g_2,\\varepsilon)", "(g_3,\\zeta_3)", "(g_3,\\zeta_3^2)"],
  :fourierMat=>[[1, 3, 2, 2,1, 3, 2, 2],[3, 3, 0, 0,-3,-3, 0, 0],
		[2, 0, 4,-2,2, 0,-2,-2],[2, 0,-2, 4, 2, 0,-2,-2],
		[1,-3, 2, 2,1,-3, 2, 2],[3,-3, 0, 0,-3, 3, 0, 0],
		[2, 0,-2,-2,2, 0, 4,-2],[2, 0,-2,-2, 2, 0,-2, 4]]//6,
  :eigenvalues=>[1,1,1,1,1,-1,E(3),E(3,2)],
  :perm=>Perm(7,8),
  :lusztig=>true, # does not satisfy (ST)^3=1 but (SPT)^3=1
  :mellin=>[[1,0,0,2,1,0,0,0],[0,1,0,0,0,1,0,0],[0,0,1,0,0,0,1,1],[1,0,0,-1,1,0,
   0,0],[1,0,0,0,-1,0,0,0],[0,1,0,0,0,-1,0,0],[0,0,1,0,0,0,E(3),E(3,2)],
   [0,0,1,0,0,0,E(3,2),E(3)]],
  :mellinLabels=>["(1,1)","(g2,1)","(g3,1)","(1,g3)","(1,g2)","(g2,g2)",
                 "(g3,g3)","(g3,g3^2)"]),
  :X=>function(p)
    ss=combinations(0:p-1,2)
    Dict(:name=>"R_{\\BZ/$p}^{\\wedge 2}",
         :explanation=>"DoubleTaft($p)",
         :charSymbols=>ss,
         :charLabels=>map(s->repr(E(p)^s[1],context=:TeX=>true)*
      "\\!\\wedge\\!"*repr(E(p)^s[2],context=:TeX=>true),ss),
    :eigenvalues=>map(s->E(p)^Product(s),ss),
    :fourierMat=>map(i->map(j->(E(p)^(i*reverse(j))-E(p)^(i*j))/p,ss),ss),
    :special=>1,:cospecial=>p-1)
   end)

function SubFamily(f,ind,scal,label)
  ind=filter(i->ind(f,i),1:length(f[:eigenvalues]))
  res=Dict(:fourierMat=>getindex.(f[:fourierMat][ind],Ref(ind))*scal,
           :eigenvalues=>f[:eigenvalues][ind],
           :charLabels=>f[:charLabels][ind],
   :operations=>FamilyOps)
  res[:name]="$(f[:name])_{[$label]}"
  if haskey(f,:charSymbols) res[:charSymbols]=f[:charSymbols][ind] end
  if haskey(f,:group) 
    res[:group]=f[:group] 
  end
  if f[:special] in ind 
   res[:special]=findfirst(isequal(ind),f[:special]) 
  end
  res
end

function SubFamilyij(f,i,j,scal)
  g=SubFamily(f,(f,k)->sum(f[:charSymbols][k])%j==i,scal,join([i,j]))
  g[:explanation]="subfamily(sum(charsymbols)mod $j=$i of $(f[:explanation]))"
  g
end

CHEVIE[:families][:X5]=SubFamilyij(CHEVIE[:families][:X](6),1,3,1-E(3))
CHEVIE[:families][:X5][:cospecial]=5

UnipotentClassOps=Dict(:Name=>x->x)

Format(x)=string(x)
Format(x,opt)=string(x)
Torus(i::Int)=i
RootsCartan=x->x
function Family(s::String,v::Vector,d::Dict=Dict{Symbol,Any}())
  f=CHEVIE[:families][Symbol(s)]
  f[:char_numbers]=v
  merge(f,d)
end
function Family(f::Dict{Symbol,Any},v::Vector,d::Dict=Dict{Symbol,Any}())
  f[:char_numbers]=v
  merge(f,d)
end

function ReadChv(s::String)
end
DiagonalMat(v...)=v
Group(v...)=v
ComplexConjugate(v)=v
function GetRoot(x,n::Number,msg::String="")
   println("GetRoot($x,$n)")
   x
end
Unbind(x)=x

include("tables.jl")
chevie[:D][:CharTable]=n->chevie[:imp][:CharTable](2,2,n)
chevie[:B][:CharTable]=n->chevie[:imp][:CharTable](2,1,n)
chevie[:A][:CharTable]=function(n)
  ct=chevie[:imp][:CharTable](1,1,n+1)
  ct[:irredinfo]=map(x->Dict(:charname=>IntListToString(x)),chevie[:A][:CharInfo](n)[:charparams])
  ct
end
chevie[:A][:HeckeCharTable]=(n,para,root)->chevie[:imp][:HeckeCharTable](1,1,n+1,para,root)
chevie[:D][:HeckeCharTable]=(n,para,root)->chevie[:imp][:HeckeCharTable](2,2,n,para,root)
chevie[:imp][:PowerMaps]=(p,q,r)->[]
end
