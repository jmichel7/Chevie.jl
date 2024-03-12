# this contains simple functions but which need
# several of the self-contained structural packages of Chevie
module Tools
export abelian_gens, abelian_invariants, improve_type
using ModuleElts
using PuiseuxPolynomials
using Combinat: tally
using PermGroups: Group, Groups, gens, word, PermGroup, elements,
         minimal_words, isabelian
using MatInt: smith_transforms
using CyclotomicNumbers: Cyc, conductor, Root1
using ..FiniteFields: FiniteFields, FFE, Z
using ..Modulo: Mod
using ..Chevie: Chevie, order

#------------------ improve_type
best_eltype(m)=reduce(promote_type,best_type.(m))
best_eltype(p::Pol)=iszero(p) ? Int : best_eltype(p.c)
best_eltype(p::Mvp)=iszero(p) ? Int : best_eltype(values(p.d))
best_type(x)=typeof(x)
best_type(x::Cyc{Rational{T}}) where T=iszero(x) ? Int : conductor(x)==1 ? 
  best_type(Rational(x)) : denominator(x)==1 ?  Cyc{T} : typeof(x)
best_type(x::Cyc{T}) where T<:Integer=conductor(x)==1 ? T : typeof(x)
best_type(x::Rational)= denominator(x)==1 ? typeof(numerator(x)) : typeof(x)
best_type(m::AbstractArray{T,N}) where {T,N}=isempty(m) ? typeof(m) : Array{best_eltype(m),N}
best_type(p::Pol)=Pol{best_eltype(p)}
best_type(x::BigInt)= typemin(Int)<=x<=typemax(Int) ? Int : BigInt
best_type(x::Root1)= order(x)<=2 ? Int : Root1
function best_type(q::Frac)
  if isone(q.den) return best_type(q.num) end
  if q.den isa Mvp && length(q.den)==1 return best_type(Mvp(q;Rational=true)) end
  if q.den isa Pol && length(q.den.c)==1 return best_type(Pol(q)) end
  Frac{promote_type(best_type(q.num), best_type(q.den))}
end
best_type(q::Frac{Pol{Rational{T}}}) where T=Frac{Pol{T}}
best_type(q::Frac{Mvp{Rational{T}}}) where T=Frac{Mvp{T}}
function best_type(p::Mvp{T,N}) where {T,N}
  if iszero(p) return  Mvp{Int,Int} end
  n=all(m->all(isinteger,powers(m)),monomials(p)) ? Int : N
  Mvp{best_eltype(p),n}
end

function improve_type(m)
  convert(best_type(m),m)
end

#------------- extend to Cyc methods of LaurentPolynomials
Base.gcd(p::Pol{<:Cyc{<:Integer}},q::Pol{<:Cyc{<:Integer}})=srgcd(p,q)

Base.numerator(p::Mvp{<:Cyc{<:Rational{<:T}},N}) where{T,N} =convert(Mvp{Cyc{T},N},p*denominator(p))

function Base.convert(::Type{Frac{Pol{Cyc{T}}}},a::Frac{<:Pol{<:Cyc{Rational{T}}}}) where T<:Integer
  n=numerator(a)
  d=denominator(a)
  Frac(numerator(n)*denominator(d),numerator(d)*denominator(n))
end

function Base.convert(::Type{Frac{Pol{T}}},a::Frac{<:Pol{<:Cyc{Rational{T}}}}) where T<:Integer
  n=convert(Pol{Rational{T}},numerator(a))
  d=convert(Pol{Rational{T}},denominator(a))
  Frac(numerator(n)*denominator(d),numerator(d)*denominator(n))
end

#function LaurentPolynomials.Frac(a::Mvp{<:Cyc{<:Rational},Int},b::Mvp{<:Cyc{<:Rational},Int};k...)
#  Frac(numerator(a)*denominator(b),numerator(b)*denominator(a);k...)
#end

"""
`valuation(c::Union{Integer,Rational{<:Integer},p::Integer)` 
`p`-adic  valuation of `c` (largest  power of `p` which  divides `c`; for a
`Rational`, valuation of the numerator minus that of the denominator).
```julia-repl
julia> valuation.(24,(2,3,5))
(3, 1, 0)
```
"""
function LaurentPolynomials.valuation(c::Integer,p::Integer)
  if iszero(c) error("first argument should not be zero") end
  e=0
  while true
    (c,r)=divrem(c,p)
    if !iszero(r) return e end
    e+=1
  end
end

LaurentPolynomials.valuation(c::Rational{<:Integer},p::Integer)=
  valuation(numerator(c),p)-valuation(denominator(c),p)

"""
`FFE{p}(z::Cyc)`  where `z` is  a `p`-integral cyclotomic  number (that is,
`z`  times some number prime  to `p` is a  cyclotomic integer), returns the
reduction  of `z` mod.  `p`, an element  of some extension  `𝔽_{pʳ}` of the
prime field `𝔽ₚ`.

```julia_repl
julia> FFE{3}(E(7))
Z₇₂₉¹⁰⁴
```
"""
function FiniteFields.FFE{p}(c::Cyc)where p
  x=coefficients(c) 
  n=conductor(c)
  pp=p^valuation(n,p)
  np=div(n,pp)
  r=order(Mod(p,np)) # order of p mod np
  zeta=Z(p^r)^(div(p^r-1,np)*gcdx(pp,p^r-1)[2]) #n-th root of unity
  if !isone(zeta^n) error() end
  sum(i->zeta^(i-1)*x[i],1:n)
end

"""
`abelian_gens(A)`
    
`A`  should be an abelian group or the list of its generators. Such a group
has  a unique decomposition up to isomorphism as a product of cyclic groups
`C(n₁)×…×C(nₖ)`  where `C(nᵢ)`  is a  cyclic group  of order  `nᵢ` and `nᵢ`
divides  `nᵢ₊₁`. The function returns a list  of generators for each of the
`C(nᵢ)`.

```julia-repl
julia> abelian_gens([Perm(1,2),Perm(3,4,5),Perm(6,7)])
2-element Vector{Perm{Int16}}:
 (6,7)
 (1,2)(3,5,4)(6,7)
```
"""
abelian_gens(l::Vector)=isempty(l) ? l : abelian_gens(Group(l))
function abelian_gens(G::Group) # thanks to Klaus Lux for the algorithm
  l=filter(!isone,gens(G))
  if isempty(l) return l end
  rels=fill(0,length(l),length(l))
  for i in eachindex(l)
    H=i==1 ? Group([one(l[1])]) : Group(l[1:i-1])
    d=findfirst(o->l[i]^o in H,1:order(l[i]))
    rel=fill(0,length(l))
    for p in tally(word(H,l[i]^d)) rel[p[1]]=p[2] end
    rel[i]=-d
    rels[i,:]=rel
  end
  rels=Integer.(inv(Rational.(smith_transforms(rels).coltrans)))
  filter(!isone,map(r->prod(l.^r),eachrow(rels)))
end

"""
`abelian_invariants(G::Group )`
    
`G`  should be an abelian group. Such a group has a unique decomposition up
to  isomorphism as a product of cyclic groups `C(n₁)×…×C(nₖ)` where `C(nᵢ)`
is  a cyclic  group of  order `nᵢ`  and `nᵢ`  divides `nᵢ₊₁`.  The function
returns the list `n₁,…,nₖ`.

```julia-repl
julia> abelian_invariants(Group(Perm(1,2),Perm(3,4,5),Perm(6,7)))
2-element Vector{Int64}:
 2
 6
```
"""
function abelian_invariants(G::Group)
  if !isabelian(G) 
    error("abelian_invariants is only implemented for abelian groups")
  end
  order.(abelian_gens(G))
end

function Base.intersect(G::PermGroup, H::PermGroup)
  if min(length(G),length(H))<=104000 
    invoke(Base.intersect,Tuple{Group,Group},G,H) # horrible implementation
  else
    m=Base.get_extension(Chevie,:Gap4)
    if isnothing(m)
       error("*** too large intersect($G,$H)\n",
           "need Gap4.intersect; do 'using GAP' to make extension available") 
    end
    m.intersect(G,H)
  end
end

end
