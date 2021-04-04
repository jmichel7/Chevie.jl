# this should just be "tools". It contains simple functions but which need
# several of the self-contained structural packages of Gapjm
export pblocks, abelian_gens

function Cycs.root(x::Pol,n::Number=2)
  n=Int(n)
  if length(x.c)>1 || !iszero(x.v%n)
    error("root($(repr(x;context=:limit=>true)),$n) not implemented") 
  end
  if isempty(x.c) return x end
  Pol([root(x.c[1],n)],div(x.v,n);check=false)
end

# next function is twice the speed of p(E(x))
function (p::Pol{T})(x::Root1) where T
  res=zero(T)
  for (i,c) in enumerate(p.c)
    res+=c*x^(valuation(p)+i-1)
  end
  res
end
  
"""
`Mvp(p)` converts  the `Pol`  `p` to  an  `Mvp`. 

```julia-repl
julia> Pol(:q)
Pol{Int64}: q

julia> Mvp(q^2+q)
Mvp{Int64}: q²+q
```
"""
Mvp(x::Pol)=convert(Mvp,x)

Base.convert(::Type{Mvp{T,N}},p::Pol) where{T,N}=
                     p(Mvp(convert(Monomial{N},Pols.varname[])=>one(T)))
Base.convert(::Type{Mvp},p::Pol)=p(Mvp(Pols.varname[]))

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
  m=zeros(T,r,r)//1
  for (e,t) in p.d
    n=map(x->findfirst(==(x),v),keys(e.d))
    c=values(e.d)
    if c==[1,1] m[n[1],n[2]]=m[n[2],n[1]]=t//2
    elseif c==[2] m[n[1],n[1]]=t
    elseif c==[1] m[n[1],r]=m[r,n[1]]=t//2
    elseif isempty(c) m[r,r]=t
    else error("# only implemented for degree <=2")
    end
  end
  if size(m,1)==2 t=one(m)
  else n=copy(m)
    _,ind=echelon(m)
    m=m[ind,:]
    if size(m,1)>2 return [p] end
    t=permutedims(solutionmat(m,n))
    m=solutionmat(t,m)
  end
  v=t*vcat(Mvp.(v),[Mvp(1)])
  if size(m,1)==1 return [v[1],v[1]*m[1,1]] end
  b=m[1,2]+m[2,1]
  if m[1,1]==0 return [v[2],b*v[1]+m[2,2]*v[2]] end
  b//=m[1,1]
  d=root(b^2-4*m[2,2]//m[1,1])
  if isnothing(d) 
    println("root failed")
    return p 
  end
  improve_type([v[1]+v[2]//2*(b-d),m[1,1]*(v[1]+v[2]//2*(b+d))])
end

"""
`mod(z::Cyc,p::Integer)`

`p`  should be a  prime and `z`  a cyclotomic number  which is `p`-integral
(that  is, `z` times some number prime to `p` is a cyclotomic integer). The
function  returns  the  reduction  of  `z`  mod.  `p`,  an  element of some
extension  `F_(p^r)` of the prime field  `F_p`.

```julia_repl
julia> mod(E(7),3)
Z₇₂₉¹⁰⁴
```
"""
function Base.mod(c::Cyc,p)
  x=coefficients(c) 
  n=conductor(c)
  np=MatInt.prime_part(n,p)
  pp=div(n,np)
  r=order(Mod{np}(p)) # order of p mod np
  zeta=Z(p^r)^(div(p^r-1,np)*gcdx(pp,p^r-1)[2]) #n-th root of unity
  if !isone(zeta^n) error() end
  sum(i->zeta^(i-1)*x[i],1:n)
end;

export pblocks
"""
`pblocks(G,p)`

Let  `p` be a prime. This function returns the partition of the irreducible
characters  of `G`  in `p`-blocks,  represented by  the list  of indices of
irreducibles characters in each block.

```julia-repl
julia> W=CoxSym(5)
𝔖 ₅

julia> pblocks(W,2)
2-element Vector{Vector{Int64}}:
 [1, 3, 4, 5, 7]
 [2, 6]

julia> pblocks(W,3)
3-element Vector{Vector{Int64}}:
 [1, 5, 6]
 [2, 3, 7]
 [4]

julia> pblocks(W,7)
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
function pblocks(G,p)
  T=CharTable(G)
  l=length(T.charnames)
  classes=map(c->div(T.centralizers[1],c),T.centralizers)
  v=map(chi->map(j->mod(classes[j]*chi[j]//chi[1],p),1:l),eachrow(T.irr))
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
  res=join(map(p.d) do (deg,v)
    den=denominator(v)
    v=numerator(v) 
    if deg==0 t=string(v)
    else 
      v=format_coefficient(string(v))
      t=v in ["","-"] ? v : v*"*"
      r=(deg==1 ? "E($(p.n))" : "E($(p.n))^$deg")
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
  res=""
  join(map(p.d)do (k,c)
    c=Util.bracket_if_needed(gap(c))
    if isempty(c) gap(k)
    else c*"*"*gap(k)
    end
  end,"+")
end
