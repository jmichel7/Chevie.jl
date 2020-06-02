module Algebras
using ..Gapjm
using LinearAlgebra: LinearAlgebra
export fusion_algebra, involution, duality, basis, dim, idempotents
abstract type FiniteDimAlgebra end

struct FusionAlgebra<:FiniteDimAlgebra
  fourier::Matrix
  special::Int
  involution::SPerm{Int16}
  duality::SPerm{Int16}
  multable::Vector{Vector{Vector{Pair}}}
  prop::Dict{Symbol,Any}
end

struct AlgebraElt{T,TA}
  A::TA
  d::ModuleElt{Int,T}
end

"""
`FusionAlgebra(f::Family)`
`FusionAlgebra(S,special=1)`

All  the Fourier matrices `S` in CHEVIE are unitary, that is `S⁻¹=conj(S)`,
and have a *special* line `s` (the line of index `s=f.special` for a family
`f`)  such that no entry `S_(s,i)` is  equal to `0`. Further, they have the
property that the sums `C_(i,j,k)=sum_l S_(i,l)
S_(j,l)conj(S_(k,l))/S_(s,l)`  take integral  values. Finally,  `S` has the
property  that complex conjugation does a permutation with signs `σ` of the
lines of `S`.

It  follows that we can define a `Z`-algebra `A` as follows: it has a basis
`b_i`  indexed by the lines of `S`, and has a multiplication defined by the
fact that the coefficient of `b_ib_j` on `b_k` is equal to `C_(i,j,k)`.

`A`  is  commutative,  and  has  as  unit  the  element  `b_s`;  the  basis
`sigma(b_i)` is dual to `b_i` for the linear form
`(b_i,b_j)=C_(i,j,sigma(s))`.

```julia-repl
julia> W=ComplexReflectionGroup(4)
G₄

julia> uc=UnipotentCharacters(W);f=uc.families[4];

julia> A=Algebras.fusion_algebra(fourier(f),1)
Fusion Algebra dim.5

julia> b=basis(A)
5-element Array{Gapjm.Algebras.AlgebraElt{Int64,Gapjm.Algebras.FusionAlgebra},1}:
 B₁
 B₂
 B₃
 B₄
 B₅

julia> b*permutedims(b)
5×5 Array{Gapjm.Algebras.AlgebraElt{Int64,Gapjm.Algebras.FusionAlgebra},2}:
 B₁  B₂      B₃      B₄        B₅
 B₂  -B₄+B₅  B₁+B₄   B₂-B₃     B₃
 B₃  B₁+B₄   -B₄+B₅  -B₂+B₃    B₂
 B₄  B₂-B₃   -B₂+B₃  B₁+B₄-B₅  -B₄
 B₅  B₃      B₂      -B₄       B₁

julia> CharTable(A)
CharTable(Fusion Algebra dim.5)
 │1    2    3  4  5
─┼──────────────────
1│1  √-3 -√-3  2 -1
2│1    1    1  .  1
3│1   -1   -1  .  1
4│1    .    . -1 -1
5│1 -√-3  √-3  2 -1
```
"""
function fusion_algebra(S::Matrix,special::Int=1;opt...)
# zero=AlgebraElt(A,zero(ModuleElt{Int,T}))
# one=AlgebraElt(A,ModuleElt(special=>1))
  involution=SPerm(collect(eachrow(S)),collect(eachrow(conj.(S))))
  if isnothing(involution) error("complex conjugacy is not SPerm(rows)") end
  if order(involution)>2 error("complex conjugacy is of order 4") end
  irr=mapslices(x->x.//x[special],S;dims=2)
  duality=SPerm(collect(eachcol(^(irr,Perm(involution)))),
                 collect(eachcol(irr)))
  if isnothing(duality) error("the matrix does not have the * involution") end
  if order(duality)>2 error("duality is not an involution") end
  s=mapslices(x->x.//conj(x[special]),conj.(S);dims=1)
  d=size(S,1)
  multable=map(i->map(j->filter(x->x[2]!=0,
            map((x,y)->x=>y,1:d,s*(S[i,:].*S[j,:]))),1:i),1:d)
  if all(r->all(c->all(p->p[2]>=0,c),r),multable)
      InfoChevie( "# positive structure constants\n" );
  end
  if !all(r->all(c->all(p->isinteger(p[2]),c),r),multable)
      error("structure constants are not integral")
  else multable=map(r->map(c->[k=>Int(i) for (k,i) in c],r),multable)
  end
  A=FusionAlgebra(S,special,involution,duality,multable,Dict{Symbol,Any}())
  d=map(ratio,eachrow(irr),eachcol(S)) # d=inv.(S[special,:]) ?
  if nothing in d  error() end
  A.prop[:cDim]=d[special]^2
  A.prop[:qDim]=d[special].//d
  A.prop[:irr]=irr
  A.prop[:charnames]=haskey(opt,:charnames) ? opt[:charnames] : string.(1:dim(A))
  A.prop[:classnames]=haskey(opt,:classnames) ? opt[:classnames] : string.(1:dim(A))
  A
end

function fusion_algebra(f::Family)
  fusion_algebra(fourier(f),get(f.prop,:special,1);charnames=f[:charLabels],
                                             classnames=f[:charLabels])
end

dim(A)=size(A.fourier,1)

Base.show(io::IO,A::FusionAlgebra)=print(io,"Fusion Algebra dim.",dim(A))

function idempotents(A::FusionAlgebra)
  gets(A,:idempotents)do
    LinearAlgebra.Diagonal(A.fourier[A.special,:])*
      conj.(permutedims(A.fourier))*basis(A)
  end
end

function basis(A::FusionAlgebra)
  gets(A,:basis)do
    map(i->AlgebraElt(A,ModuleElt(i=>1)),1:dim(A))
  end
end

function Base.show(io::IO, h::AlgebraElt)
  function showbasis(io::IO,e)
    repl=get(io,:limit,false)
    TeX=get(io,:TeX,false)
    res="B"
    if repl || TeX
      res*="_{"*string(e)*"}"
    else res*="("*string(e)*")"
    end
    fromTeX(io,res)
  end
  show(IOContext(io,:showbasis=>showbasis),h.d)
end

Base.:+(a::AlgebraElt,b::AlgebraElt)=AlgebraElt(a.A,a.d+b.d)
Base.:-(a::AlgebraElt)=AlgebraElt(a.A,-a.d)
Base.:-(a::AlgebraElt, b::AlgebraElt)=a+(-b)
Base.:*(a::AlgebraElt, b)=AlgebraElt(a.A,a.d*b)
Base.:*(b,a::AlgebraElt)=AlgebraElt(a.A,a.d*b)
Base.zero(a::AlgebraElt)=AlgebraElt(a.A,zero(a.d))
Base.:^(a::AlgebraElt, n::Integer)=Base.power_by_squaring(a,n)

function Cycs.coefficients(a::AlgebraElt{T})where T
  v=fill(T(0),dim(a.A))
  for (i,c) in a.d v[i]=c end
  v
end

function Base.:*(a::AlgebraElt{T1}, b::AlgebraElt{T2})where {T1,T2}
  res=Pair{Int,promote_type(T1,T2)}[]
  for (i,c) in a.d
    for (i1,c1) in b.d
      l=(i>=i1) ? a.A.multable[i][i1] : a.A.multable[i1][i]
      append!(res,[k=>c*c1*c2 for (k,c2) in l])
    end
  end
  AlgebraElt(a.A,ModuleElt(res;check=true))
end

function Chars.CharTable(A::FusionAlgebra)
  irr=toM(map(e->map(b->ratio(coefficients(b*e),coefficients(e)), basis(A)),
              idempotents(A)))
  if irr!=A.prop[:irr] error() end
  labels=string.(1:dim(A))
  centralizers=fill(dim(A),dim(A))
  CharTable(irr,A.prop[:charnames],A.prop[:classnames],centralizers,
            string(A),Dict{Symbol,Any}())
end

function involution(e::AlgebraElt{FusionAlgebra})
  p=Perm(e.A.involution)
  s=signs(e.A.involution)
  AlgebraElt(e.A,ModuleElt([Int(b^p)=>c*s[b] for (b,c) in e.d]))
end

function duality(e::AlgebraElt{FusionAlgebra})
  p=Perm(e.A.duality)
  s=signs(e.A.duality)
  AlgebraElt(e.A,ModuleElt([Int(b^p)=>c*s[b] for (b,c) in e.d]))
end
end
