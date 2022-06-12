# this contains simple functions but which need
# several of the self-contained structural packages of Gapjm
module Tools
export abelian_gens, improve_type
using UsingMerge
using ModuleElts
using LaurentPolynomials
using PuiseuxPolynomials

using ..Combinat: Combinat, tally
using ..FFields: FFields, FFE, Mod, Z
using ..Groups: Group, gens, word
using ..Gapjm: Gapjm, Cyc, conductor, order
using ..MatInt: smith_transforms

#------------------ improve_type
best_eltype(m)=reduce(promote_type,best_type.(m))
best_eltype(p::Pol)=iszero(p) ? Int : best_eltype(p.c)
best_eltype(p::Mvp)=iszero(p) ? Int : best_eltype(values(p.d))
best_type(x)=typeof(x)
best_type(x::Cyc{Rational{T}}) where T=iszero(x) ? Int : conductor(x)==1 ? 
  best_type(Rational(x)) : denominator(x)==1 ?  Cyc{T} : typeof(x)
best_type(x::Cyc{T}) where T<:Integer=conductor(x)==1 ? T : typeof(x)
best_type(x::Rational)= denominator(x)==1 ? typeof(numerator(x)) : typeof(x)
best_type(m::Array{T,N}) where {T,N}=isempty(m) ? typeof(m) : Array{best_eltype(m),N}
best_type(p::Pol)=Pol{best_eltype(p)}
function best_type(q::Frac)
  if isone(q.den) return best_type(q.num) end
  Frac{promote_type(best_type(q.num), best_type(q.den))}
end
function best_type(p::Mvp{T,N}) where {T,N}
  if iszero(p) return  Mvp{Int,Int} end
  n=all(m->all(isinteger,values(m.d)),keys(p.d)) ? Int : N
  Mvp{best_eltype(p),n}
end
  
improve_type(m)=convert(best_type(m),m)

#------------------
Base.gcd(p::Pol{<:Cyc{<:Integer}},q::Pol{<:Cyc{<:Integer}})=srgcd(p,q)

Base.numerator(p::Mvp{<:Cyc{<:Rational{<:T}},N}) where{T,N} =convert(Mvp{Cyc{T},N},p*denominator(p))

function Frac(a::Mvp{<:Cyc{<:Rational},Int},b::Mvp{<:Cyc{<:Rational},Int};k...)
  Frac(numerator(a)*denominator(b),numerator(b)*denominator(a);k...)
end

LaurentPolynomials.exactdiv(m::ModuleElt,b)=merge(exactdiv,m,b)

LaurentPolynomials.exactdiv(c::Cyc{<:Integer},b::Integer)=Cyc(
                    conductor(c),exactdiv(c.d,b))

function LaurentPolynomials.exactdiv(a::Cyc{<:Integer},b::Cyc{<:Integer})
  res=a//b
  if denominator(res)>1 error(b," does not exactly divide ",a) end
  numerator(res)
end

"valuation(c::Integer,p::Integer) p-adic valuation of c"
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
function FFields.FFE{p}(c::Cyc)where p
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
has  a unique decomposition  up to isomorphism  of the form `C₁×…×Cₙ` where
the  order of `Cᵢ` divides the order of `Cᵢ₊₁`. The function returns a list
of generators for each of the `Cᵢ`.

```julia-repl
julia> abelian_gens([Perm(1,2),Perm(3,4,5),Perm(6,7)])
2-element Vector{Perm{Int16}}:
 (1,2)(6,7)
 (3,5,4)(6,7)
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

end
