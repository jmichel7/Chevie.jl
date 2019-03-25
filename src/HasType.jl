module HasType

export charinfo, classinfo, reflection_name, diagram, chartable,
  representation, fakedegrees, unipotent_characters, 
  schur_elements, charname, codegrees, ComplexReflectionGroup,
  getclassified, chevieget, field, getfromtype, Cartesian

using Gapjm

function Cartesian(a::AbstractVector...)
  reverse.(vec(collect.(Iterators.product(reverse(a)...))))
end

braid_relations(W)=impl1(getclassified(W,:BraidRelations))

function codegrees(W)
  vcat(map(refltype(W)) do t
    cd=getfromtype(t,:ReflectionCoDegrees)
    if isnothing(cd)
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
  gets(W,:charinfo) do W
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
  res
  end
end

function classinfo(W)
  gets(W,:classinfo) do W
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
    res
  end
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
  if i in [23,28,30,36,37,38]
    if i==23     return coxgroup(:H,3)
    elseif i==28 return coxgroup(:F,4)
    elseif i==30 return coxgroup(:H,4)
    elseif i==36 return coxgroup(:E,6)
    elseif i==37 return coxgroup(:E,7)
    elseif i==38 return coxgroup(:E,8)
    end
    m=getfromtype(t,:CartanMat)
    n=one(hcat(m...))
    return PermRootGroup(map(i->n[i,:],axes(n,1)),m)
  end
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

PrintToSring(s,v...)=sprint(show,v...)

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

UnipotentClassOps=Dict(:Name=>x->x)

unipotent_characters(W)=impl1(getclassified(W,:UnipotentCharacters))

unipotent_classes(W,p=0)=impl1(getclassified(W,:UnipotentClasses,p))

#-----------------------------------------------------------------------
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
Base.getindex(s::String,a::Vector{Any})=getindex(s,Int.(a))
Cycs.:^(a::Cyc,b::Rational)=a^Int(b)

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

Binomial=binomial
IdentityMat(n)=map(i->one(rand(Int,n,n))[i,:],1:n)

function pad(s::String, i::Int)
  if i>0 return lpad(s,i)
  else return rpad(s,-i)
  end
end

pad(s::String)=s

function Replace(s,p...)
# println("Replace s=$s p=$p")
  r=[p[i]=>p[i+1] for i in 1:2:length(p)]
  for (src,tgt) in r
    i=1
    while i+length(src)-1<=length(s)
      if src==s[i:i+length(src)-1]
        if tgt isa String
          s=s[1:i-1]*tgt*s[i+length(src):end]
        else
          s=vcat(s[1:i-1],tgt,s[i+length(src):end])
        end
      end
      i+=1
    end
  end
  s
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

Combinations=combinations
Drop(a::Vector,i::Int)=deleteat!(copy(a),i)

Base.getindex(a::Symbol,i::Int)=string(a)[i]
Base.length(a::Symbol)=length(string(a))

IsList(l)=l isa Vector
IsInt(l)=l isa Int ||(l isa Rational && denominator(l)==1)

Flat(v)=collect(Iterators.flatten(v))

ForAll(l,f)=all(f,l)

PartitionTuples=partition_tuples
Arrangements=arrangements
Partitions=partitions

function PartitionTupleToString(n,a=Dict())
  if n[end] isa Vector return join(map(join,n),".") end
  r=repr(E(n[end-1],n[end]),context=:limit=>true)
  if r=="1" r="+" end
  if r=="-1" r="-" end
  join(map(join,n[1:end-2]),".")*r
end

include("symbols.jl")

Product(v)=isempty(v) ? 1 : prod(v)
Product(v,f)=isempty(v) ? 1 : prod(f,v)

IntListToString(l)=any(x->x>10,l) ? join(l,",") : join(l)

Lcm(a...)=Lcm(collect(a))
Lcm(a::Vector)=lcm(Int.(a))

function Collected(v)
  d=groupby(v,v)
  sort([[k,length(v)] for (k,v) in d])
end

function CollectBy(v,f)
  d=groupby(f,v)
  [d[k] for k in sort(collect(keys(d)))]
end

SortBy(x,f)=sort!(x,by=f)

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
Filtered(l,f)=isempty(l) ? l : filter(f,l)

gapSet(v)=unique(sort(v))
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

FamilyImprimitive = function (S,)
local e, Scoll, ct, d, m, ll, eps, equiv, nrSymbols, epsreps, trace, roots, i, j, mat, frobs, symbs, newsigns, schon, orb, mult, res, IsReducedSymbol
  println("S=$S")
  e = length(S)
  Scoll = Collected(vcat(S...))
  ct = vcat(map(x->fill(x[1],x[2]), Scoll)...)
  d = length(ct) % e
  if !(d in [0, 1]) error("Length(", IntListToString(ct), ") should be 0 or 1  %  ", e, " !\n")
        end
  m = div(length(ct) - d, e)
  j = (m * binomial(e, 2)) % e
  ll = Cartesian(map(i->0:e-1, Scoll)...)
  ll = filter(x->sum(x)%e==j,ll)
  ll = map(c->map((x,y)->filter(c->sum(c)%e==y,
                   collect(combinations(0:e-1,x[2]))),Scoll,c), ll)
  nrSymbols=sum(x->prod(length,x),ll)
  ll = vcat(map(x->Cartesian(x...), ll)...)
  eps = l->(-1)^sum(i->count(j->l[i]<j,l[i+1:end]),1:length(l))
  equiv = map(x->
      Dict(:globaleps=>length(x)==1 ? 1 :
     (-1)^sum(i->sum(j->sum(y->count(k->y<k,j),x[i]),x[i+1:end]),1:length(x)-1),
     :aa=>map(y->map(x->(l=x,eps=eps(x)),arrangements(y, length(y))),x)), ll)
  epsreps = map(x->eps(vcat(x...)), ll)
  roots = map(i->E(e,i),0:e-1)
  mat = map(i->i[:globaleps]*map(k->epsreps[k]*
  prod(l->sum(j->j.eps*roots[1+mod(-sum(map((a,b)->a*b,j.l,ll[k][l])),e)],
              i[:aa][l]),
          1:length(i[:aa])), 1:nrSymbols), equiv)
  mat = ((-1)^(m*(e-1))*mat)//(E(4,binomial(e-1,2))*ER(e)^e)^m
  frobs = E(12,-(e^2-1)*m)*map(i->E(2e,-(sum(j->j*j,i))-e*sum(sum,i)),ll)
  symbs = map(function (l)local sy, j
              sy = map(j->Int[], 1:e)
              map((v,c)->begin push!(sy[v + 1], c)
                      return 1 end, vcat(l...), ct)
              return sy
          end, ll)
  newsigns = (-1) ^ (binomial(e, 2) * binomial(m, 2)) * map(i->
               (-1)^((0:e - 1)*map(x->binomial(length(x), 2), i)),symbs)
  mat = map((s,l)->s * map((x, y)->x*y, newsigns,l),newsigns,mat)
  if d == 0
  IsReducedSymbol(s)=all(x->s==x || LessSymbols(x, s),Rotations(s)[2:length(s)])
    schon = map(IsReducedSymbol, symbs)
    mult = []
    for i = 1:nrSymbols
        if schon[i]
            orb = gapSet(Rotations(symbs[i]))
            push!(mult, e // length(orb))
            for j = Filtered(i + 1:nrSymbols, (j->symbs[j] in orb))
                schon[j] = false
            end
        end
    end
    frobs = vcat(map((m,f)->fill(f,m), mult, ListBlist(frobs, schon))...)
    symbs = vcat(map(function (m, s)
                    if m==1 return [s]
                    else return map(j->vcat(s[1:e//m], [m, j]), 0:m-1)
                    end
                end, mult, ListBlist(symbs, schon))...)
    mat = vcat(map(function (m, l)return map((i->begin
      vcat(map((n,c)->fill(e*c//m//n,n),mult,ListBlist(l, schon))...)
                         end), 1:m)end, mult, ListBlist(mat, schon))...)
    mult=vcat(map(m->fill(m,m),mult)...)
    nrSymbols=length(symbs)
    for i=1:nrSymbols
      for j=1:nrSymbols
        if FullSymbol(symbs[i])==FullSymbol(symbs[j])
            mat[i][j]-=1//mult[i]
            if symbs[i]==symbs[j] mat[i][j]+=1 end
        end
      end
    end
    if (mat*DiagonalMat(frobs))^3!=mat^0
        print("** WARNING: (S*T)^3!=1\n")
    end
  end
  res=Dict{Symbol,Any}(:symbols=>symbs,
    :fourierMat=>mat,
    :eigenvalues=>frobs,
    :name=>IntListToString(ct),
    :explanation=>"classical family",
    :special=>1,
    :operations=>FamilyOps)
  res[:charLabels] = map(string, 1:length(res[:symbols]))
  res[:size] = length(res[:symbols])
  res
end
MakeFamilyImprimitive = function (S, uc)
  f=x->Position(uc[:charSymbols],x)
  if length(S)==1 return Family("C1", map(f, S)) end
  r = FamilyImprimitive(FullSymbol(S[1]))
  r[:charNumbers] = map(f, r[:symbols])
  r[:special] = PositionProperty(r[:charNumbers],(x->uc[:a][x] == uc[:b][x]))
  r[:cospecial] = PositionProperty(r[:charNumbers],(x->uc[:A][x] == uc[:B][x]))
# if length(blocks(r[:fourierMat])) > 1 error() end
  r
end
WeylGroup(s::String,n)=coxgroup(Symbol(s),Int(n))
#-------------------------------------------------------------------------
#  dummy translations of GAP3 functions
CHEVIE=Dict{Symbol,Any}(:compat=>Dict(:MakeCharacterTable=>x->x,
                           :AdjustHeckeCharTable=>(x,y)->x,
        :ChangeIdentifier=>function(tbl,n)tbl[:identifier]=n end))

DiagonalMat(v...)=cat(map(m->m isa Array ? m : hcat(m),v)...,dims=(1,2))
DiagonalOfMat(m)=[m[i,i] for i in axes(m,1)]

include("families.jl")
Format(x)=string(x)
Format(x,opt)=string(x)
FormatTeX(x)=repr(x,context=:TeX=>true)
Torus(i::Int)=i
RootsCartan=x->x

function ReadChv(s::String)
end
Group(v...)=v
ComplexConjugate(v)=v
function GetRoot(x,n::Number=2,msg::String="")
   println("GetRoot($x,$n) returns $x")
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
