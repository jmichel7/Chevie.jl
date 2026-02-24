"""
This  module  deals  with  symmetric  functions.  For details look at the books
[mac15](@cite) or [ze81](@cite).

The  algebra `R` of symmetric functions over the integers is isomorphic to
the  sum of the Grothendieck groups ``⊕_{n≥0}R[𝔖ₙ]``. The multiplication in
the algebra is given by ``\\hbox{Ind}_{𝔖ₙ×𝔖ₘ}^{𝔖_{n+m}}``.

This module implements the following bases of `R`:
  - `p` (power sums)
  - `h` (complete  symmetric  functions)
  - `e` (elementary  symmetric functions) 
  - `m` (monomial symmetric functions)
  - `s` (Schur functions)
Note  that `p`  is only  a basis  over the  rationals (we  accept arbitrary
coefficients).  Each  of  these  bases  in  degree  `n`  is  indexed by the
partitions of `n`.

Interpreted  as class functions  on `𝔖ₙ`, the  basis `s` corresponds to the
irreducible  characters,  the  basis  `p`  to the normalized characteristic
functions of the classes, the basis `h` to the induced of the identity from
Young  subgroups (standard parabolic  subgroups), and the  basis `e` to the
induced of the sign from Young subgroups.

A typical session would begin by defining these bases:
```julia-rep1
julia> using .SymFuncs: p,h,s,e,m
```
This  is  necessary  since  is  would  not  be a good policy to pollute the
namespace by exporting all these 1-letter names.

Then make elements in one of these bases. The following forms are equivalent:
```julia-rep1
julia> p(Partition(3,2,1))
p₃₂₁

julia> p([3,2,1])
p₃₂₁

julia> p(3,2,1)
p₃₂₁
```
The functions can be used to convert between bases:
```julia-rep1
julia> s(p(3,2,1))
-s₁₁₁₁₁₁-s₂₂₂+s₃₁₁₁+s₃₃-s₄₁₁+s₆

julia> h(p(3,2,1))
-h₁₁₁₁₁₁+5h₂₁₁₁₁-6h₂₂₁₁-3h₃₁₁₁+6h₃₂₁

julia> h(p(3,2,1))[3,2,1] # you can find a coefficient by indexing
6//1
```
The following operations are defined on symmetric functions, in addition to
`+` and `-`:
```julia-rep1
julia> s(2,1)*s(2,1) # product
s₂₂₁₁+s₂₂₂+s₃₁₁₁+2s₃₂₁+s₃₃+s₄₁₁+s₄₂

julia> s(2,1)⊗s(2,1) # inner product
s₁₁₁+s₂₁+s₃

julia> scalar_product(p(1,1,1),s(2,1))
2//1
```
One can mix bases in these operations. The basis of the left argument wins:
```julia-rep1
julia> s(2,1)+p(3)
s₁₁₁+s₃
```
The plethysm is implemented with two possible notations:

```julia-rep1
julia> plethysm(p(2,1),s(2,1))
(1//9)p₂₂₂₁₁₁+(-1//9)p₃₂₂₂+(-1//9)p₆₁₁₁+(1//9)p₆₃

julia> p(2,1)[s(2,1)]  # the same thing
(1//9)p₂₂₂₁₁₁+(-1//9)p₃₂₂₂+(-1//9)p₆₁₁₁+(1//9)p₆₃

julia> @Mvp u,v

julia> p(2,1)[u*p(2)+v*p(3)] # plethysm acts on Mvp coefficients
u³p₄₂+u²vp₄₃+uv²p₆₂+v³p₆₃
```

finally one can convert a symmetric function to a symmetric polynomial.
The number of variables by default is the highest degree in the function.
```julia-rep1
julia> Mvp(p(2)+p(3))
Mvp{Int64}: x₁³+x₁²+x₂³+x₂²+x₃³+x₃²

julia> Mvp(p(2),[:u,:v,:w]) # one can choose the variables used and their number
Mvp{Int64}: u²+v²+w²
```

#### Wreath Symmetric functions

We also implement symmetric functions for the group `Gₑ,₁,ₙ`, isomorphic to
`μₑ≀𝔖ₙ=(μₑⁿ)⋊𝔖ₙ`,  where `μₑ`  is the  group of  `e`-th roots of unity. The
algebra `Rₑ` of wreath symmetric function is the sum of Grothendieck groups
`⊕ₙR(Gₑ,₁,ₙ)`, with product given by
``\\hbox{Ind}_{Gₑ,₁,ₘ×Gₑ,₁,ₙ}^{Gₑ,₁,ₙ₊ₘ}``.

The  algebra `Rₑ` is  isomorphic to a  tensor product of  `e` copies of the
algebra  `R` of symmetric functions,  indexed by the irreducible characters
of  `μₑ`. The copy indexed by  ``γ∈\\hbox{Irr(μₑ)}`` has as basis in degree
`n` the irreducible characters appearing in
``\\hbox{Ind}_{μₑⁿ}^{Gₑ,₁,ₙ}γ^{⊗n}``.

Each  of the  bases of  `R` gives  thus by  tensor product  a basis of `Rₑ`
indexed  by  `e`-tuples  of  partitions.  The  basis  `s`  still represents
irreducible  characters,  but  the  basis  `p` does not represent conjugacy
classes.

The  conjugacy class of  an element `(ζ₁,…,ζₙ).σ∈(μₑⁿ)⋊𝔖ₙ`  is obtained has
follows:   decompose  `σ`  in  cycles,  and   record  the  product  of  the
corresponding `ζᵢ`. This distributes the cycles in `e` classes, thus builds
an  `e`-tuple of partitions. The  normalized characteristic function of the
corresponding  class is represented  by a basis  `π` indexed by `e`-tuples.
Then  the character table  of `Gₑ,₁,ₙ` is  encoded in the decomposition
of the basis `π` in the basis `s`.

To use the basis `π`, we recommend that you do
```julia-rep1
julia> Pi=SymFuncs.π;
```end
in order not to destroy the constant `π`.

The following are equivalent
```julia-rep1
julia> p(2,1)⊠p(1)
p₂₁.₁

julia> p([[2,1],[1]])
p₂₁.₁

julia> p(PartitionTuple([2,1],[1]))
p₂₁.₁

julia> p([2,1],[1])
p₂₁.₁
```
for  the basis `Pi` only the last 3  methods are allowed since it cannot be
build as an external tensor.
```julia-rep1
julia> s(Pi([1],[1])) # The values of the characters on the class 1.1
-s.₁₁-s.₂+s₁₁.+s₂.

julia> Pi(p([2],[1]))
(-1//4)π.₂₁+(1//4)π₁.₂+(-1//4)π₂.₁+(1//4)π₂₁.

julia> Pi(p([2],[1]))[[2],[1]] # you can get a coefficient by indexing
-1//4
```

The  function `Mvp`  converts a  wreath symmetric  function to  a symmetric
polynomial  using a different set of variables for each factor of the tensor
product:
```julia-rep1
julia> Mvp(p(3)⊠p(2))
Mvp{Int64}: x₁³y₁²+x₁³y₂²+x₂³y₁²+x₂³y₂²+x₃³y₁²+x₃³y₂²
```

For the scalar product the bases `p` and `Pi` are orthogonal and the basis `s`
is orthonormal

```julia-rep1
julia> scalar_product(Pi([1],[1]),Pi([1],[1]))
4

julia> scalar_product(p([1],[1]),p([1],[1]))
1
```
"""
module SymFuncs
using ..Chevie
using ..Symbols: z
export plethysm, ⊠

const partitionscache=Dict{Int,Vector{Partition}}()

function cpartitions(n)
  get!(partitionscache,n)do
    Partition.(partitions(n))
  end
end

const mncache=Dict{Pair{Vector{Int},Vector{Int}},Int}((Int[]=>Int[])=>1,
                                                      ([1]=>[1])=>1)
# Murnaghan-Nakayama rule: b a betaset p a partition
function mn(b::AbstractVector,p::AbstractVector)
  get!(mncache,b=>p)do
    e=p[1]
    sum(enumerate(b))do (i,j)
      if j<e return 0 end
      r=searchsorted(b,j-e)
      if !isempty(r) return 0 end
      sgn=isodd(i-r.start) ? -1 : 1
      b1=copy(b)
      for k in i:-1:r.start+1 b1[k]=b1[k-1] end; b1[r.start]=j-e
      sgn*mn(shiftβ(b1),@view p[2:end])
    end
  end
end

function charvalue(p::Partition,q::Partition)
  if rank(p)!=rank(q) error() end
  mn(βset(p),q.l)
end

function H(lambda,s)
  n=sum(lambda)
  t=UnipotentValues(UnipotentClasses(rootdatum(:gl,n));classes=true)
  i=findfirst(==(lambda),partitions(n))
  sum(t.scalar[:,i].*s.(partitions(n)))
end

struct SymFunc{b,C}
  d::ModuleElt{Partition,C}
  SymFunc{b}(m)where b =new{b,valtype(m)}(m)
end

function Base.show(io::IO, h::SymFunc{b})where b
  function showbasis(io::IO,w)
    res=string(b)
    if hasdecor(io) && !get(io,:naive,false) res*="_{"*xrepr(io,w)*"}"
    else  res*="("*join(w.l,",")*")"
    end
    fromTeX(io,res)
  end
  show(rio(io,limit=true,showbasis=showbasis,naive=!hasdecor(io)),improve_type(h.d))
end

clone(h::SymFunc{b},d) where b=SymFunc{b}(d)

Base.:+(a::SymFunc{ba},b::SymFunc) where ba=clone(a,a.d+basis(Val(ba),b).d)
Base.:+(a::SymFunc, b::Union{Number,Pol,Mvp})=a+one(a)*b
Base.:-(a::SymFunc)=clone(a,-a.d)
Base.:-(a::SymFunc{ba},b::SymFunc) where ba=clone(a,a.d-basis(Val(ba),b).d)
Base.:*(a::SymFunc, b::Union{Number,Pol,Mvp,Frac})=clone(a,a.d*b)
Base.:(//)(a::SymFunc, b::Union{Number,Pol,Mvp})=clone(a,a.d//b)
Base.:*(b::Union{Number,Pol,Mvp,Frac}, a::SymFunc)=a*b

Base.:^(a::SymFunc, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                      Base.power_by_squaring(inv(a),-n)
Base.one(a::SymFunc)=clone(a,ModuleElt(Partition()=>1))
Base.zero(a::SymFunc)=SymFunc{basis(a)}(zero(a.d))
Base.zero(b::Symbol,a::SymFunc)=SymFunc{b}(zero(a.d))
s(x...)=basis(Val(:s),x...)
p(x...)=basis(Val(:p),x...)
h(x...)=basis(Val(:h),x...)
e(x...)=basis(Val(:e),x...)
m(x...)=basis(Val(:m),x...)

Algebras.basis(h::SymFunc{b})where b=b

function Algebras.basis(::Val{b},p::Partition)where b
  if b==:π error("basis π takes only PartitionTuples") end
  SymFunc{b}(ModuleElt(p=>1))
end

function Algebras.basis(::Val{b},w1::Vector{<:Integer},
     w::Vararg{<:Vector{<:Integer}})where b
  if length(w)==0 basis(Val(b),Partition(w1))
  else basis(Val(b),PartitionTuple(pushfirst!(collect(w),w1)))
  end
end

function Algebras.basis(::Val{b},w::Vararg{<:Integer})where b
  basis(Val(b),convert(Vector{Int},collect(w)))
end

Algebras.basis(::Val{b},h::SymFunc{b}) where b =h

function Algebras.basis(::Val{b1},h::SymFunc{b}) where {b,b1}
  if b==:p error(":p to ",b1," not implemented") end
  basis(Val(b1),p(h))
end

function Algebras.basis(::Val{:p},h::SymFunc{:s})
   sum(h.d)do (p,c)
       SymFunc{:p}(ModuleElt(l=>charvalue(p,l)*c//z(l)
            for l in cpartitions(rank(p))))
   end
end
Algebras.basis(::Val{:p},h::SymFunc{:h})=sum(((μ,c),)->
  prod(i->p(basis(Val(:s),i)),μ.l;init=p())*c,h.d;init=zero(:p,h))
Algebras.basis(::Val{:p},h::SymFunc{:e})=sum(((μ,c),)->
  prod(i->p(basis(Val(:s),fill(1,i))),μ.l;init=p())*c,h.d;init=zero(:p,h))

Algebras.basis(::Val{:s},h::SymFunc{:p})=
  sum(h.d;init=zero(:s,h))do (p,c)
    SymFunc{:s}(ModuleElt(l=>charvalue(l,p)*c for l in
                     cpartitions(rank(p))))
  end

function Pto(n,b) # p(n) to basis b∈(:e,:h)
  SymFunc{b}(ModuleElt(map(cpartitions(n))do p
    p=>n*(-1)^(length(p)-1)*factorial(length(p)-1)//
    prod(i->factorial(i[2]),tally(p.l))
   end))*((b==:e && iseven(n)) ? -1 : 1)
end
function Algebras.basis(::Val{:h},f::SymFunc{:p})
  sum(((p,c),)->prod(n->Pto(n,:h),p.l;init=h())*c,f.d;init=zero(:h,f))
end
function Algebras.basis(::Val{:e},f::SymFunc{:p})
  sum(((p,c),)->prod(n->Pto(n,:e),p.l;init=e())*c,f.d;init=zero(:e,f))
end
function Algebras.basis(::Val{:p},f::SymFunc{:m})
  sum(f.d)do (μ,c)
    SymFunc{:p}(ModuleElt(l=>h(p(l))[μ]//z(l) for l in cpartitions(rank(μ))))*c
  end
end
function Algebras.basis(::Val{:m},f::SymFunc{:p})
  sum(f.d)do (μ,c)
    SymFunc{:m}(ModuleElt(l=>p(h(l))[μ] for l in cpartitions(rank(μ))))*z(μ)*c
  end
end

Base.getindex(a::SymFunc,p::Partition)=a.d[p]
Base.getindex(a::SymFunc,x::Vararg{Int})=a.d[Partition(collect(x))]

function Base.:*(f::SymFunc{b},g::SymFunc)where b
  g=basis(Val(b),g)
  if b in (:p,:h,:e)
    SymFunc{b}(ModuleElt([union(p,q)=>c*c1 for (p,c) in f.d for (q,c1) in g.d]))
  else basis(Val(b),p(f)*p(g))
  end
end

"`a⊗b` is the inner product of the symmetric functions `a` and `b`"
function FinitePosets.:⊗(f::SymFunc{b},g::SymFunc)where b
  m=ModuleElts.merge2((x,y)->x*y,p(f).d,p(g).d)
  r=SymFunc{:p}(ModuleElt(p=>c*z(p) for (p,c) in m))
  basis(Val(b),r)
end

function PuiseuxPolynomials.Mvp(a::SymFunc,
  vars=[Symbol("x",stringind(rio(),j)) for j in 1:rank(last(a.d.d)[1])])
  improve_type(sum(p(a).d)do (l,c)
                prod(i->sum(Mvp.(vars).^i),l.l)*c
  end)
end

"`scalar_product(a,b)` the scalar product of the symmetric functions `a` and `b`"
function Chars.scalar_product(a::SymFunc{ab},b::SymFunc)where ab
  b=basis(Val(ab),b)
  if ab==:s
    sum(values(ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d)))
  elseif ab==:p
    sum(((p,c),)->z(p)*c,ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d))
  else scalar_product(p(a),p(b))
  end
end

# raises all variables to rth power
pow(m::Monomial,r)=Monomial(m.d*r)
pow(p::Mvp,r)=Mvp(ModuleElt(pow(m,r)=>c for (m,c) in p.d))
pow(p::Frac{<:Mvp},r)=Frac(pow(p.num,r),pow(p.den,r))
pow(p,r)=p

"""
`plethysm(a::SymFunc,b::SymFunc)` or `a[b]`.

Plethysm  is an operation which associates  to every symmetric function `f`
a function `g↦f[g]` on the  ring of symmetric functions which
is defined by the rules
  - `(f₁+f₂)[g]=f₁[g]+f₂[g]`
  - `f₁*f₂[g]=f₁[g]*f₂[g]`
  - `pₙ[g₁+g₂]=pₙ[g₁]+[g₂]`
  - `pₙ[g₁*g₂]=pₙ[g₁]*[g₂]`
  - `pₙ[pₘ]=pₙₘ`
In addition, for `Mvp` or `Frac{Mvp}` coefficients, `g↦pₙ[g]` is no more an
algebra  homomorphism  but  acts  semi-linearly.  If  `a`  is an `Mvp` or a
`Frac{Mvp}`  we have `pₙ[a]=SymFuncs.pow(a,n)` where `Symfuncs.pow(a,n)` is
the operation which raises all variables of `a` to the `n`-th power.
```julia-rep1
julia> @Mvp u,v

julia> p(3)[(u+v^2)p()]
(u³+v⁶)p₋

julia> p(3,2)[(u+v^2)p()]
(u⁵+u³v⁴+u²v⁶+v¹⁰)p₋
```
We  implement also the plethysm for hyperoctaedral groups, which associates
to  each element `f` of `R₂` a function `R₂⊗R₂→R₂:(f₀,f₁)↦ f[f₀,f₁])` which
is defined by the rules:
  - `(f+f')[f₀,f₁]=f[f₀,f₁]+f'[f₀,f₁]`
  - `(f*f')[f₀,f₁]=f[f₀,f₁]*f'[f₀,f₁]`
  - `(pₙ⊠1)[f₀,f₁]=(pₙ⊠1)[f₀,0]`
  - `(pₙ⊠1)[f₀+f'₀,0]=(pₙ⊠1)[f₀,0]+(pₙ⊠1)[f'₀,0]`
  - `(pₙ⊠1)[a⊠b,0]=pₙ[a]⊠pₙ[b]`
  - `(1⊠pₙ)[f₀,f₁]=(1⊠pₙ)[0,f₁]`
  - `(1⊠pₙ)[0,f₁+f'₁]=(1⊠pₙ)[0,f₁]+(1⊠pₙ)[0,f'₁]`
  - `(1⊠pₙ)[0,a⊠b]=pₙ[a]⊠pₙ[b]`
"""
function plethysm(a::SymFunc{ab},b::SymFunc)where ab
  a=p(a);b=p(b)
  basis(Val(ab),sum(a.d;init=zero(a))do (P,c)
    c*prod(P.l;init=p())do i
      SymFunc{:p}(ModuleElt(Partition(i*p1.l)=>pow(c1,i) for (p1,c1) in b.d))
    end
  end)
end

Base.getindex(a::SymFunc,b::SymFunc)=plethysm(a,b)

struct WreathSymFunc{b,C}
  d::ModuleElt{PartitionTuple,C}
  WreathSymFunc{b}(m)where b =new{b,valtype(m)}(m)
end

function Base.show(io::IO, h::WreathSymFunc{b})where b
  function showbasis(io::IO,w)
    res=string(b)
    s=xrepr(io,w)
    if hasdecor(io) && !get(io,:naive,false) res*="_{"*s*"}"
    else  res*="("*join(w.P,",")*")"
    end
    fromTeX(io,res)
  end
  show(rio(io,limit=true,showbasis=showbasis,naive=!hasdecor(io)),improve_type(h.d))
end

"""
`π(λ::PartitionTuple)`  Normalized characteristic function of a class

The   basis  `π`  represents  the  characteristic  function  of  the  class
parameterized  by the `e`-tuple  of partitions `λ` times the cardinality of
the  centralizer in `Gₑ,₁,ₙ` of an element  of that class (where `n` is the
`rank` of `λ`). This cardinality can be obtained by `Symbols.z(λ)`.

It  can be  related to  the basis  `p` as  follows: both  bases satisfy the
property  that  the  product  does  the  union  of  the  partition  tuples:
`p(λ)*p(μ)==p(union(λ,μ))`  and similarly for `π`. And if `[n]ᵢ` represents
the `PartitionTuple` which is empty at all places excepted at the place `i`
where it is the partition `[n]` we have

``p(π([n]ᵢ))=∑_{j=1}^{j=e}ζ_e^{(1-j)*(i-1)}p([n]ⱼ)``

and

``π(p([n]ᵢ))=(∑_{j=1}^{j=e}ζ_e^{(j-1)*(i-1)}π([n]ⱼ))/e``

```julia_repl
julia> p(Pi(Int[],Int[],[3],Int[]))
-ζ₄p...₃-p..₃.+ζ₄p.₃..+p₃...
```
"""
π(x...)=basis(Val(:π),x...)

Algebras.basis(h::WreathSymFunc{b})where b=b
clone(h::WreathSymFunc{b},d) where b=WreathSymFunc{b}(d)

Base.:+(a::WreathSymFunc{ba},b::WreathSymFunc) where ba=clone(a,a.d+basis(Val(ba),b).d)
Base.:+(a::WreathSymFunc, b::Union{Number,Pol,Mvp})=a+one(a)*b
Base.:-(a::WreathSymFunc)=clone(a,-a.d)
Base.:-(a::WreathSymFunc{ba},b::WreathSymFunc) where ba=clone(a,a.d-basis(Val(ba),b).d)
Base.:*(a::WreathSymFunc, b::Union{Number,Pol,Mvp,Frac})=clone(a,a.d*b)
Base.:(//)(a::WreathSymFunc, b::Union{Number,Pol,Mvp})=clone(a,a.d//b)
Base.:*(b::Union{Number,Pol,Mvp,Frac}, a::WreathSymFunc)=a*b
Base.zero(a::WreathSymFunc)=WreathSymFunc{basis(a)}(zero(a.d))
Base.zero(b::Symbol,a::WreathSymFunc)=WreathSymFunc{b}(zero(a.d))

⊠(a,a1)=tens(a,a1)

function tens(aa::SymFunc...)
  if length(aa)<2 error() end
  b=basis(aa[1])
  for a in aa[2:end] a=basis(Val(b),a) end
  m=ModuleElt(map(cartesian(map(x->pairs(x.d),aa)...))do pp
               PartitionTuple(first.(pp))=>prod(last.(pp))
              end)
  WreathSymFunc{b}(m)
end

function tens(a::WreathSymFunc,b::SymFunc)
  A=b.A
  b=basis(A,Val(basis(a)),b)
  m=ModuleElt(map(cartesian(pairs(a.d),pairs(b.d)))do ((P,c),(p,c1))
               PartitionTuple(vcat(P.P,[p.l]))=>c*c1
              end)
  WreathSymFunc{basis(a)}(m)
end

Algebras.basis(::Val{b},p::PartitionTuple)where b=
  WreathSymFunc{b}(ModuleElt(p=>1))
Algebras.basis(::Val{b},w::Vector{<:Vector{<:Integer}})where b=
  basis(Val(b),PartitionTuple(w))
  
function Algebras.basis(::Val{b},h::WreathSymFunc{b1})where {b,b1}
  if b==b1 return h end
  if b==:π || b1==:π error("uuuu") end
  sum(h.d;init=zero(b,h))do (P,c)
   tens(map(p->basis(Val(b),basis(Val(b1),p)),P.P)...)*c
  end
end

LaurentPolynomials.degree(p::WreathSymFunc)=degree(first(first(p.d)))
    
Algebras.basis(::Val{:π},h::WreathSymFunc{:π})=h

function Algebras.basis(::Val{:π},h::WreathSymFunc{b})where b
  function toπ(n,j,e)
    if n==0 WreathSymFunc{:π}(ModuleElt(PartitionTuple(fill(Int[],e))=>1))
    elseif e==2 WreathSymFunc{:π}(
        ModuleElt(PartitionTuple([n],Int[])=>1//2,
                  PartitionTuple(Int[],[n])=>j==1 ? 1//2 : -1//2))
    else
      WreathSymFunc{:π}(ModuleElt(
       map(1:e)do i
          v=fill(Int[],e);v[i]=[n]
          PartitionTuple(v)=>E(e,(i-1)*(j-1))//e
       end))
    end
  end
  h=basis(Val(:p),h)
  e=degree(h)
  sum(h.d)do(P,c)
    prod(eachindex(P.P))do j
      prod(eachindex(P.P[j]);init=toπ(0,j,e))do i
        toπ(P.P[j][i],j,e)
      end
    end*c
  end
end
  
function Algebras.basis(::Val{b},h::WreathSymFunc{:π})where b
  function fromπ(n,i,e)
    if n==0 WreathSymFunc{:p}(ModuleElt(PartitionTuple(fill(Int[],e))=>1))
    elseif e==2 
       WreathSymFunc{:p}(
        ModuleElt(PartitionTuple([n],Int[])=>1,
                  PartitionTuple(Int[],[n])=>i==1 ? 1 : -1))
    else
        WreathSymFunc{:p}(ModuleElt(
         map(1:e)do j
            v=fill(Int[],e);v[j]=[n]
            PartitionTuple(v)=>E(e,-(i-1)*(j-1))
         end))
    end
  end
  e=degree(h)
  basis(Val(b),sum(h.d)do(P,c)
    prod(eachindex(P.P))do j
      prod(eachindex(P.P[j]);init=fromπ(0,j,e))do i
        fromπ(P.P[j][i],j,e)
      end
    end*c
  end)
end

function Base.:*(f::WreathSymFunc{b},g::WreathSymFunc)where b
  g=basis(Val(b),g)
  if b in (:p,:π,:h,:e)
    WreathSymFunc{b}(ModuleElt(
      [union(p,q)=>c*c1 for (p,c) in f.d for (q,c1) in g.d]))
  else basis(Val(b),p(f)*p(g))
  end
end

function PuiseuxPolynomials.Mvp(a::WreathSymFunc,vars=nothing)
  e=length(a.d.d[1][1].P)
  a=p(a)
  if isnothing(vars)
    varnames="xyztuvw"
    vars=map(i->[Symbol(varnames[i],stringind(rio(),j)) 
               for j in 1:maximum(map(k->sum(first(k).P[i]),a.d.d))],1:e)
  end
  improve_type(sum(a.d)do (P,c)
           prod(i->prod(k->sum(Mvp.(vars[i]).^k),P.P[i];init=Mvp(1)),1:e)*c
  end)
end

function Chars.scalar_product(a::WreathSymFunc{ab},b::WreathSymFunc)where ab
  b=basis(Val(ab),b)
  res=if ab==:s
    sum(values(ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d)))
  elseif ab==:p
    sum(((p,c),)->prod(z,p.P)*c,ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d))
  elseif ab==:π
    sum(((p,c),)->z(p)*c,ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d))
  else scalar_product(p(a),p(b))
  end
  improve_type(res)
end

function FinitePosets.:⊗(f::WreathSymFunc{b},g::WreathSymFunc)where b
  m=ModuleElts.merge2((x,y)->x*y,π(f).d,π(g).d)
  r=WreathSymFunc{:π}(ModuleElt(p=>c*z(p) for (p,c) in m))
  basis(Val(b),r)
end

function plethysm(g::WreathSymFunc{b},f0::WreathSymFunc,f1::WreathSymFunc)where b
  g=p(g);f0=p(f0);f1=p(f1)
  basis(Val(b),(sum(g.d;init=zero(g))do (P,c)
    c*prod(1:2)do i
      prod(P.P[i];init=p(Int[],Int[]))do n
       sum(i==1 ? f0.d : f1.d;init=zero(g))do (P1,c1)
          p(n)[p(P1.P[1])*c1]⊠p(n)[p(P1.P[2])]
        end
      end
    end
  end))
end

Base.getindex(g::WreathSymFunc,f0::WreathSymFunc,f1::WreathSymFunc)=plethysm(g,f0,f1)

Base.getindex(a::WreathSymFunc,p::PartitionTuple)=a.d[p]
Base.getindex(a::WreathSymFunc,x::Vararg{<:Vector{<:Int}})=a.d[PartitionTuple(collect(x))]

#@test mytest("SymFuncs.jl","SymFuncs.p(Partition(3,2,1))","p₃₂₁")
#@test mytest("SymFuncs.jl","SymFuncs.p([3,2,1])","p₃₂₁")
#@test mytest("SymFuncs.jl","SymFuncs.p(3,2,1)","p₃₂₁")
#@test mytest("SymFuncs.jl","SymFuncs.s(SymFuncs.p(3,2,1))","-s₁₁₁₁₁₁-s₂₂₂+s₃₁₁₁+s₃₃-s₄₁₁+s₆")
#@test mytest("SymFuncs.jl","SymFuncs.h(SymFuncs.p(3,2,1))","-h₁₁₁₁₁₁+5h₂₁₁₁₁-6h₂₂₁₁-3h₃₁₁₁+6h₃₂₁")
#@test mytest("SymFuncs.jl","SymFuncs.h(SymFuncs.p(3,2,1))[3,2,1]","6//1")
#@test mytest("SymFuncs.jl","SymFuncs.s(2,1)*SymFuncs.s(2,1)","s₂₂₁₁+s₂₂₂+s₃₁₁₁+2s₃₂₁+s₃₃+s₄₁₁+s₄₂")
#@test mytest("SymFuncs.jl","SymFuncs.s(2,1)⊗SymFuncs.s(2,1)","s₁₁₁+s₂₁+s₃")
#@test mytest("SymFuncs.jl","scalar_product(SymFuncs.p(1,1,1),SymFuncs.s(2,1))","2//1")
#@test mytest("SymFuncs.jl","SymFuncs.s(2,1)+SymFuncs.p(3)","s₁₁₁+s₃")
#@test mytest("SymFuncs.jl","plethysm(SymFuncs.p(2,1),SymFuncs.s(2,1))","(1//9)p₂₂₂₁₁₁+(-1//9)p₃₂₂₂+(-1//9)p₆₁₁₁+(1//9)p₆₃")
#@test mytest("SymFuncs.jl","SymFuncs.p(2,1)[SymFuncs.s(2,1)]","(1//9)p₂₂₂₁₁₁+(-1//9)p₃₂₂₂+(-1//9)p₆₁₁₁+(1//9)p₆₃")
#@test mytest("SymFuncs.jl","@Mvp u,v","nothing")
#@test mytest("SymFuncs.jl","SymFuncs.p(2,1)[u*SymFuncs.p(2)+v*SymFuncs.p(3)]","u³p₄₂+u²vp₄₃+uv²p₆₂+v³p₆₃")
#@test mytest("SymFuncs.jl","Mvp(SymFuncs.p(2)+SymFuncs.p(3))","Mvp{Int64}: x₁³+x₁²+x₂³+x₂²+x₃³+x₃²")
#@test mytest("SymFuncs.jl","Mvp(SymFuncs.p(2),[:u,:v,:w])","Mvp{Int64}: u²+v²+w²")
#@test mytest("SymFuncs.jl","Pi=SymFuncs.π;","")
#@test mytest("SymFuncs.jl","SymFuncs.p(2,1)⊠SymFuncs.p(1)","p₂₁.₁")
#@test mytest("SymFuncs.jl","SymFuncs.p([[2,1],[1]])","p₂₁.₁")
#@test mytest("SymFuncs.jl","SymFuncs.p(PartitionTuple([2,1],[1]))","p₂₁.₁")
#@test mytest("SymFuncs.jl","SymFuncs.p([2,1],[1])","p₂₁.₁")
#@test mytest("SymFuncs.jl","SymFuncs.s(Pi([1],[1]))","-s.₁₁-s.₂+s₁₁.+s₂.")
#@test mytest("SymFuncs.jl","Pi(SymFuncs.p([2],[1]))","(-1//4)π.₂₁+(1//4)π₁.₂+(-1//4)π₂.₁+(1//4)π₂₁.")
#@test mytest("SymFuncs.jl","Pi(SymFuncs.p([2],[1]))[[2],[1]]","-1//4")
#@test mytest("SymFuncs.jl","Mvp(SymFuncs.p(3)⊠SymFuncs.p(2))","Mvp{Int64}: x₁³y₁²+x₁³y₂²+x₂³y₁²+x₂³y₂²+x₃³y₁²+x₃³y₂²")
#@test mytest("SymFuncs.jl","scalar_product(Pi([1],[1]),Pi([1],[1]))","4")
#@test mytest("SymFuncs.jl","scalar_product(SymFuncs.p([1],[1]),SymFuncs.p([1],[1]))","1")
#@test mytest("SymFuncs.jl","@Mvp u,v","nothing")
#@test mytest("SymFuncs.jl","SymFuncs.p(3)[(u+v^2)SymFuncs.p()]","(u³+v⁶)p₋")
#@test mytest("SymFuncs.jl","SymFuncs.p(3,2)[(u+v^2)SymFuncs.p()]","(u⁵+u³v⁴+u²v⁶+v¹⁰)p₋")

end
