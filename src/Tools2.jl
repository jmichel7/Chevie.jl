module Tools2
export eigmat, traces_words_mats, Unknown

using PuiseuxPolynomials
using LaurentPolynomials
using CyclotomicNumbers
using Primes: Primes
using PermGroups: Group
using ..Util: bracket_if_needed, format_coefficient
using ..Chars: CharTable
using ..Combinat: Combinat, collectby
using ..GLinearAlgebra: solutionmat, echelon, charpoly
using ..Tools: improve_type
using ..FFields: FFE
using ..CycPols: CycPol
using ..Gapjm: Gapjm, root, gap, Cyc, conductor

#----------------------------------------------------------------------
"""
`blocks(G::Group,p::Integer)`

Let  `p` be a prime. This function returns the partition of the irreducible
characters  of `G`  in `p`-blocks,  represented by  the list  of indices of
irreducibles characters in each block.

```julia-repl
julia> W=CoxSym(5)
𝔖 ₅

julia> blocks(W,2)
2-element Vector{Vector{Int64}}:
 [1, 3, 4, 5, 7]
 [2, 6]

julia> blocks(W,3)
3-element Vector{Vector{Int64}}:
 [1, 5, 6]
 [2, 3, 7]
 [4]

julia> blocks(W,7)
7-element Vector{Vector{Int64}}:
 [1]
 [2]
 [3]
 [4]
 [5]
 [6]
 [7]
```
"""
function Combinat.blocks(G::Group,p::Integer)
  T=CharTable(G)
  l=length(T.charnames)
  classes=map(c->div(T.centralizers[1],c),T.centralizers)
  v=map(chi->map(j->FFE{p}(classes[j]*chi[j]//chi[1]),1:l),eachrow(T.irr))
  sort(collectby(improve_type(v),1:l))
end

"""
`factor(p::Mvp)`

`p`  should be of degree <=2 thus represents a quadratic form. The function
returns  a list  of two  linear forms  of which  `p` is the product if such
exist, otherwise it returns [p].

```julia-repl
julia> @Mvp x,y

julia> factor(x^2-y^2+x+3y-2)
2-element Vector{Mvp{Int64, Int64}}:
 x-y+2
 x+y-1

julia> factor(x^2+x+1)
2-element Vector{Mvp{Cyc{Int64}, Int64}}:
 x-ζ₃
 x-ζ₃²

julia> factor(x*y-1)
1-element Vector{Mvp{Int64, Int64}}:
 xy-1
```
"""
function Primes.factor(p::Mvp{T,N})where {T,N}
  v=variables(p)
  r=length(v)+1
  m=zeros(T,r,r)*1//1
  for (e,t) in p.d
    n=map(x->findfirst(==(x),v),keys(e.d))
    c=collect(values(e.d))
    if c==[1,1] m[n[1],n[2]]=m[n[2],n[1]]=t*1//2
    elseif c==[2] m[n[1],n[1]]=t
    elseif c==[1] m[n[1],r]=m[r,n[1]]=t*1//2
    elseif isempty(c) m[r,r]=t
    else error("# only implemented for degree <=2")
    end
  end
  if size(m,1)==2 t=one(m)
  else n=copy(m)
    m=m[echelon(m)[2],:] # independent lines
    if size(m,1)>2 return [p] end
    t=transpose(solutionmat(m,n))
    m=solutionmat(t,m)
  end
  v=t*vcat(Mvp.(v),[Mvp(1)])
  if size(m,1)==1 return [v[1],v[1]*m[1,1]] end
  b=m[1,2]+m[2,1]
  if m[1,1]==0 return [v[2],b*v[1]+m[2,2]*v[2]] end
  b/=m[1,1]
  d=root(b^2-4*m[2,2]/m[1,1])
  if isnothing(d) 
    println("root failed")
    return p 
  end
  improve_type([v[1]+v[2]*1//2*(b-d),m[1,1]*(v[1]+v[2]*1//2*(b+d))])
end

" eigenvalues as Cycs of a matrix of finite order"
function eigmat(m)
  l=[fill(e,m) for (e,m) in CycPol(Pol(charpoly(m))).v]
  if isempty(l) Root1[] else vcat(l...) end
end

"`CycPol(x::Mvp)` converts univariate `Mvp` `x` to a `CycPol`"
function CycPol(x::Mvp)
  if !isinteger(degree(x)) x
  else CycPol(Pol(x))
  end
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
4-element Vector{Matrix{Int64}}:
 [-1 -1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
 [1 0 0 0; -1 -1 -1 0; 0 0 1 0; 0 0 0 1]
 [1 0 0 0; 0 1 0 0; 0 -2 -1 -1; 0 0 0 1]
 [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 -1 -1]

julia> traces_words_mats(R,r)==CharTable(W).irr[17,:]
true
```
"""
function traces_words_mats(mats,words)
  mats=improve_type(mats)
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
#   println(prod(dens[w])
    trace(w)
  end
end
#-----------------------------------------------------------------------
struct Unknown end

Base.:+(a,b::Unknown)=b
Base.:+(b::Unknown,a)=b
Base.:+(b::Unknown,a::Unknown)=b
Base.:*(a,b::Unknown)=iszero(a) ? a : b
Base.:*(b::Unknown,a)=iszero(a) ? a : b
Base.:*(b::Unknown,a::Unknown)=b
Base.zero(a::Unknown)=0
Base.show(io::IO,a::Unknown)=print(io,get(io,:limit,false) ? "?" : "Unknown()")
Base.isless(a::Unknown,b::Number)=false
Base.isless(b::Number,a::Unknown)=true
Base.isless(b::Unknown,a::Unknown)=false
Base.:(//)(a::Unknown,b)=a
Base.isreal(a::Unknown)=false
Base.isinteger(a::Unknown)=false
Base.conj(a::Unknown)=a
Base.broadcastable(a::Unknown)=Ref(a)
CyclotomicNumbers.galois(a::Unknown,i)=a
#-------------------- define function to backport to gap3 some values
Gapjm.gap(p::Rational)=string(numerator(p),"/",denominator(p))

Gapjm.gap(p::Integer)=string(p)

function Gapjm.gap(p::Cyc)
  res=join(map(pairs(p)) do (deg,v)
    den=denominator(v)
    v=numerator(v) 
    if deg==0 t=string(v)
    else 
      v=format_coefficient(string(v))
      t=v in ["","-"] ? v : v*"*"
      r=(deg==1 ? "E($(conductor(p)))" : "E($(conductor(p)))^$deg")
      t*=r
    end
    if t[1]!='-' t="+"*t end
    if !isone(den) t*="/$den" end
    t
  end)
  if res=="" res="0"
  elseif res[1]=='+' res=res[2:end] 
  end
  res
end

Gapjm.gap(v::AbstractVector)="["*join(gap.(v),",")*"]"

Gapjm.gap(m::AbstractMatrix)=gap(toL(m))

function Gapjm.gap(p::Pol)
  if iszero(p) return "0*q" end
  l=filter(x->x[2]!=0,collect(pairs(p.c)))
  join(map(((i,c),)->"($(gap(c)))*q^$(i+p.v-1)",l),"+")
end

Gapjm.gap(m::Monomial)=join([string(v)*(p==1 ? "" : string("^",p)) for (v,p) in m.d],"*")

function Gapjm.gap(p::Mvp)
  res=join(map(p.d)do (k,c)
    c=bracket_if_needed(gap(c))
    kk=gap(k)
    if isempty(c) kk
    elseif isempty(kk) c
    else c*"*"*kk
    end
  end,"+")
  isempty(res) ? "0" : res
end

Gapjm.gap(f::Float64)="evalf(\"$f\")"
Gapjm.gap(f::Complex{Float64})="Complex("*gap(real(f))*","*gap(imag(f))*")"

function Gapjm.gap(p::CycPol) # export positive CycPols to GAP3
  if any(<(0),values(p.v)) error("non-positive") end
  res=string("CycPol([",gap(p.coeff),",",p.valuation,",")
  res*join(map(x->join(map(gap,fill(x[1].r,x[2])),","),pairs(p.v)),",")*"])"
end

end
