# pirating extensions to get closer to GAP semantics
Base.:*(a::AbstractVector{<:Number},b::AbstractVector{<:Number})=transpose(a)*b
Base.:*(a::AbstractVector,b::AbstractVector{<:AbstractVector})=toL(toM(a)*toM(b))
Base.:*(a::AbstractVector{<:Number},b::AbstractVector{<:AbstractVector})=toL(transpose(a)*toM(b))[1]
Base.:*(a::AbstractVector,b::AbstractVector)=sum(a.*b)
Base.:*(a::Tuple,b::AbstractVector)=toL(transpose(collect(a))*toM(b))[1]
Base.:+(a::AbstractArray,b::Number)=a .+ b
Base.:+(a::Integer,b::AbstractVector)=a .+ b
Base.:-(a::AbstractVector,b::Number)=a .- b
Base.:-(a::Integer,b::AbstractVector)=a .- b
Base.:^(m::AbstractMatrix,n::AbstractMatrix)=inv(n*1//1)*m*n
Base.:^(m::AbstractVector{<:AbstractVector{<:Number}},n::Matrix{<:Number})=inv(n*1//1)*toM(m)*n
Base.:^(m::AbstractVector{<:AbstractVector},n::Integer)=toL(toM(m)^n)
Base.:^(m::AbstractVector,n::AbstractVector)=toL(inv(toM(n)*1//1)*toM(m)*toM(n))
Base.getindex(s::String,a::Vector{Any})=getindex(s,Int.(a))
Base.getindex(a::Symbol,i::Int)=string(a)[i]
Base.isless(a::Array,b::Number)=true # to sort [[1,2], 2, 0]
Base.isless(b::Number,a::Array)=false
Base.copy(x::Char)=x
Base.:(//)(m::Vector,n::Vector)=toL(toM(m)*inv(toM(n)*E(1)))
Base.inv(m::Vector)=toL(inv(toM(m)*E(1)//1))
Base.length(a::Symbol)=length(string(a))
Base.union(v::Vector)=union(v...)

# other extensions
Base.:*(a::AbstractArray,b::Pol)=a .* b
Base.:*(a::Pol,b::AbstractArray)=a .* b
Base.:*(a::AbstractArray,b::Frac)=a .* b
Base.:*(a::Frac,b::AbstractArray)=a .* b
Base.:*(a::AbstractArray,b::Mvp)=a .* b
Base.:*(a::Mvp,b::AbstractArray)=a .* b
Base.:*(W1::Spets,W2::FiniteCoxeterGroup)=Cosets.extprod(W1,spets(W2))
Base.:*(W1::FiniteCoxeterGroup,W2::Spets)=Cosets.extprod(spets(W1),W2)
Base.:+(a::AbstractArray,b::Pol)=a .+ b
Base.:/(a::AbstractArray,b::Pol)=a ./ b
Base.://(a::AbstractArray,b::Pol)=a .// b
Base.:+(a::AbstractArray,b::Mvp)=a .+ b
Base.:/(a::AbstractArray,b::Mvp)=a ./ b
Base.://(a::AbstractArray,b::Mvp)=a .// b
CyclotomicNumbers.:^(a::Cyc,b::Rational)=a^Int(b)
Base.:^(f::Family,n::Int)=galois(f,n)
Posets.Poset(m::Vector{Vector{Bool}})=Poset(toM(m))
PermRoot.roots(C::Vector{<:Vector})=roots(toM(C))

# translations of GAP3 functions for the Chevie library
ApplyWord(w,gens)=isempty(w) ? one(gens[1]) : prod(i->i>0 ? gens[i] : inv(gens[-i]),w)
CartanMat(s,a...)=cartan(Symbol(s),a...)
CharParams(W)=charinfo(W).charparams
CharRepresentationWords(mats,words)=traces_words_mats(toM.(mats),words)
CollectBy(v,f)=collectby(f,v)
ConcatenationString(s...)=prod(s)
CoxeterGroup()=coxgroup()
CyclePermInt(p::Perm,i::Integer)=orbit(p,i)
function CycleStructurePerm(p::Perm)
  if isone(p) return Int[] end
  t=cycletype(p)
  res=[0,nothing]
  resize!(res,t[1]-1)
  res.=nothing
  for (k,v) in tally(t) res[k-1]=v end
  res
end
function DiagonalMat(v...)
  arg=map(v)do m
    if m isa Array 0 .+toM(m)
    else 0 .+hcat(m)
    end
  end
  R=cat(arg...;dims=(1,2))
  u=reshape(R,(prod(size(R),)))
  for i in eachindex(u) if !isassigned(u,i) u[i]=zero(u[1]) end end
  toL(R)
end
DiagonalMat(v::Vector)=DiagonalMat(v...)
Drop(a,i::Int)=deleteat!(collect(a),i) # a AbstractVector ot Tuple
EltWord(W,x)=W(x...)
ExteriorPower(m,i)=toL(exterior_power(toM(m),i))
Factors(n)=reduce(vcat,fill(k,v) for (k,v) in factor(n))
HeckeCentralMonomials=central_monomials
function Ignore() end
Inherit(a,b)=merge!(a,b)
function Inherit(a,b,c)
  for k in c a[Symbol(k)]=b[Symbol(k)] end
  a
end
Intersection(x,y)=sort(intersect(x,y))
Intersection(x::Vector)=sort(intersect(x...))
Join(x,y)=join(x,y)
Join(x)=join(x,",")
KroneckerProduct(a,b)=toL(kron(toM(a),toM(b)))
LongestCoxeterWord(W)=word(W,longest(W))
NrConjugacyClasses(W)=length(classinfo(W)[:classtext])
OnMatrices(a::Vector{<:Vector},b::Perm)=(a.^b)^b
import Primes
phi=Primes.totient
ReflectionSubgroup(W,I::AbstractVector)=reflection_subgroup(W,convert(Vector{Int},I))
function RootInt(n,k=2)
  res=floor(Int,n^(1/k))
  if (res+1)^k<=n res+1 else res end
end
Rotations(a)=circshift.(Ref(a),length(a):-1:1)
SchurFunctor(m,p)=toL(schur_functor(toM(m),p))
SymmetricDifference(x,y)=sort(symdiff(x,y))
SymmetricPower(m,n)=SchurFunctor(m,[n])
function SortParallel(a,b)
  v=sortperm(a)
  b.=b[v]
  a.=a[v]
end
StringToDigits(s)=map(y->Position("01234567890", y), collect(s)).-1
SymbolsDefect(a,b,c,d)=symbols(a,b,d)
function TeXBracket(s)
  s=string(s)
  length(s)==1  ? s : "{"*s*"}"
end
Weyl.torus(m::Vector{<:Vector})=torus(toM(m))
Value(p,v)=p(v)
ValuePol(v,c)=isempty(v) ? 0 : evalpoly(c,v)
function CoxeterGroup(S::String,s...)
 res=coxgroup(Symbol(S),Int(s[1])) 
 for i in 2:2:length(s) res*=coxgroup(Symbol(s[i]),Int(s[i+1])) end
 res
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

ChevieIndeterminate(a::Vector{<:Number})=one(Pol)
ChevieIndeterminate(a::Vector{<:Pol})=Mvp(:x)

"""
`CycPol(v::AbstractVector)`

This  form is an  compact way unsed  in the Chevie  library of specifying a
`CycPol`  with only  positive multiplicities:  `v` should  be a vector. The
first  element is taken as the `.coeff`  of the `CycPol`, the second as the
`.valuation`.   Subsequent  elements  are   rationals  `i//d`  representing
`(q-E(d)^i)` or are integers `d` representing `Φ_d(q)`.

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
      for j in prime_residues(k) push!(vv,E(k,j)=>1) end
    else
      push!(vv,Root1(;r=i)=>1)
    end
  end
  CycPol(coeff,valuation,ModuleElt(vv))
end
#-------------------------------------------------------------------------
#  dummy translations of GAP3 functions
Format(x)=string(x)
FormatTeX(x)=repr(x,context=:TeX=>true)
FormatGAP(x)=replace(repr(x)," "=>"")
Format(x,opt)=repr(x;context=IOContext(stdout,opt...))

function ReadChv(s) end
Groups.Group(a::Perm...)=Group(collect(a))
GetRoot(x::Cyc,n::Number=2,msg...)=root(x,n)
GetRoot(x::Integer,n::Number=2,msg...)=root(x,n)
GetRoot(x::Pol,n::Number=2,msg...)=root(x,n)
GetRoot(x::Rational,n::Number=2,msg...)=root(x,n)
GetRoot(x::Mvp,n::Number=2,msg...)=root(x,n)
GetRoot(x::Root1,n::Number=2,msg...)=root(x,n)
function GetRoot(x,n::Number=2,msg...)
  error("GetRoot($x,$n) not implemented")
end

Unbind(x)=x
#-------------------------------------------------------------------------
Cosets.spets(W::FiniteCoxeterGroup,v::Vector)=spets(W,toM(v))

FamilyOps=Dict()

FactorizedSchurElementsOps=Dict{Symbol,Any}(
:Simplify=>r->HeckeAlgebras.Simplify(HeckeAlgebras.FactSchur(r[:factor],
          map(x->(pol=x[:pol], monomial=Mvp(x[:monomial])), r[:vcyc]))))

function exceptioCharName(para)
  res=string("\\phi_{",para[1],",",para[2],"}")
  if length(para)==3 res*="'"^para[3] end
  res
end
