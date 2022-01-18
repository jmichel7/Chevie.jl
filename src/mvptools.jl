# this should just be "tools". It contains simple functions but which need
# several of the self-contained structural packages of Gapjm
export abelian_gens

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
function Util.factor(p::Mvp{T,N})where {T,N}
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
    t=permutedims(solutionmat(m,n))
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
function Combinat.blocks(G,p)
  T=CharTable(G)
  l=length(T.charnames)
  classes=map(c->div(T.centralizers[1],c),T.centralizers)
  v=map(chi->map(j->FFE{p}(classes[j]*chi[j]//chi[1]),1:l),eachrow(T.irr))
  sort(collectby(improve_type(v),1:l))
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
 (6,7)
 (1,2)(3,5,4)(6,7)
```
"""
abelian_gens(l::Vector)=isempty(l) ? l : abelian_gens(Group(l))
function abelian_gens(G::Group)
# thanks to Klaus Lux for the algorithm
  l=filter(!isone,gens(G))
  if isempty(l) return l end
  rels=Vector{Int}[]
  for i in eachindex(l)
    H=i==1 ? Group([one(l[1])]) : Group(l[1:i-1])
    d=findfirst(o->l[i]^o in H,1:order(l[i]))
    rel=fill(0,length(l))
    for p in tally(word(H,l[i]^d)) rel[p[1]]=p[2] end
    rel[i]=-d
    push!(rels,rel)
  end
  rels=inv(toM(MatInt.SmithNormalFormIntegerMatTransforms(rels)[:coltrans])//1)
  rels=Int.(rels)
  filter(!isone,map(r->prod(l.^r),eachrow(rels)))
end

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
  l=filter(x->x[2]!=0,collect(enumerate(p.c)))
  join(map(((i,c),)->"($(gap(c)))*q^$(i+p.v-1)",l),"+")
end

Gapjm.gap(m::Monomial)=join([string(v)*(p==1 ? "" : string("^",p)) for (v,p) in m.d],"*")

function Gapjm.gap(p::Mvp)
  res=join(map(p.d)do (k,c)
    c=Util.bracket_if_needed(gap(c))
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
