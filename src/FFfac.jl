# code for factoring a polynomial over a finite field
module FFfac

using Primes: Primes
using LinearAlgebra: exactdiv

using LaurentPolynomials: Pol, degree, shift, derivative
using ..FFields: FFE, GF

"""
`factors_same_degree(f::Pol{FFE{p}},d,F::GF)where p`

find  the irreducible factors of `f` assumed to be a square free product of
monic  irreducible  factors  of  degree  `d`  over  a  finite field `F` of
characteristic `p`.
"""
function factors_same_degree(f::Pol{FFE{p}}, d,F)where p
  if degree(f)==d return [f] end
  g=copy(f)
  while true
    g=Pol(rand(F,2d))
    k=maximum(degree.(f.c))
    if p==2 # take g+g^2+g^(2^2)+ ... +g^(2^(k*d-1)) for GF(2^k)
      h=g
      for i in 1:k*d-1
        g=powermod(g,2,f)
        h+=g
      end
    else # take g^((p^(k*d)-1)/2)-1
      h=powermod(g, div(p^(k*d)-1,2),f)-1
    end
    g=gcd(f,h)
    if degree(f)>degree(g)>0 break end
  end
  vcat(factors_same_degree(exactdiv(f,g),d,F),factors_same_degree(g,d,F))
end

"""
`factors_squarefree(f::Pol{FFE{p}},F::GF)where p`

find  the irreducible  factors of  `f` assumed  to be  a square  free monic
polynomial   with  nonzero  constant  term  over  a  finite  field  `F`  of
characteristic `p`.
"""
function factors_squarefree(f::Pol{FFE{p}},F)where p
  facs=Pol{FFE{p}}[]
  deg=0
  k=degree(F)
  q=Pol([one(FFE{p})],1)
  pow=powermod(q, p^k, f)
  # in the following pow=q^(p^(k(deg+1)))
  while 2*(deg+1)<=degree(f) # while f could still have two irreducible factors
    deg+=1
    cyc=pow-q
    pow=powermod(pow,p^k,f)
    g=gcd(f,cyc)
    if degree(g)>0
      append!(facs, factors_same_degree(g,deg,F))
      f=exactdiv(f,g)
    end
  end
  if degree(f)>0 push!(facs, f) end
  facs
end

# n-th root of a pol which is assumed to be an n-th power
# not exported since does not check that f is a power
function root(f::Pol{FFE{p}},n::Integer)where p
  d=maximum(degree.(f.c))
  z=Z(p^d)
  r=map(0:div(degree(f),n)) do i
    e=f[i*n]
    iszero(e) ? zero(z) : z^(div(log(e),n))
  end
  Pol(r,div(f.v,n))
end

"""
`factor(f::Pol{FFE{p}}[, F])`

Given  `f` a polynomial  over a finite  field of characteristic `p`, factor
`f`,  by default over the  field of its coefficients,  or if specified over
the field `F`.

```julia-repl
julia> @Pol q
Pol{Int64}: q

julia> f=(q^4-1)*Z(3)^0
Pol{FFE{3}}: q⁴-1

julia> factor(f)
3-element Vector{Pol{FFE{3}}}:
 q²+1
 q+1
 q-1

julia> factor(f,GF(9))
4-element Vector{Pol{FFE{3}}}:
 q+1
 q-1
 q+Z₉²
 q+Z₉⁶
```
"""
function Primes.factor(f::Pol{FFE{p}},F=GF(p^maximum(degree.(f.c))))where p
  facs=Pol{FFE{p}}[]
  # make the polynomial unitary, remember the leading coefficient for later
  l=f[end]
  f=f/l
  append!(facs,map(x->Pol(), 1:f.v))
  f=shift(f,-f.v)
  if degree(f)==1 facs=[f]
  elseif degree(f)>=2
    d=derivative(f)
    if iszero(d) # f is the p-th power of another polynomial
      h=factor(root(f,p),F)
      facs=vcat(fill(h,p)...)
    else
      g=gcd(f,d)
      facs=factors_squarefree(exactdiv(f,g),F)
      for h in facs
        while true
          g1,r=divrem(g,h)
          if !iszero(r) break end
          g=g1
          push!(facs, copy(h))
        end
      end
      if degree(g)>1 append!(facs, factor(g,F)) end
    end
  end
  sort!(facs)
  facs[1]*=l
  facs
end

end
