"""
This   package  introduces  finite  fields   using  the  GAP  syntax.  This
compatibility   with  GAP  is  the  motivation  not  to  use  the  existing
`GaloisFields`.  The  speed  is  comparable  with  `GaloisFields`, slightly
slower  for prime fields and faster for composite fields. Lke GAP3, we only
implement  fields  of  order  less  than  2^16. This package comes with the
module  `Modulo` implementing modular arithmetic without restriction on the
modulus (the modulus can be a `BigInt`).

This only dependency of this package is `Primes`.

The Galois field with `p^n` elements is obtained as `GF(p^n)`. All elements
of  Galois fields of characteristic `p`  have the same type, the parametric
type   `FFE{p}`.  The  function   `Z(p^n)`  returns  a   generator  of  the
multiplicative group of `GF(p^n)`. Other elements of `GF(p^n)` are obtained
as  powers of `Z(p^n)`, except `0`, obtained as `0*Z(p^n)`. Elements of the
prime  field can  also be  obtained as  `FFE{p}(n)` (which  is the  same as
`n*Z(p)^0`).

```julia-repl
julia> a=Z(64)
FFE{2}: Z₆₄

julia> a^9 # automatic conversion to smaller fields
FFE{2}: Z₈

julia> a^21
FFE{2}: Z₄

julia> a+1
FFE{2}: Z₆₄⁵⁶
```
Elements  of the prime field can be converted to `Mod(,p)` or to integers:

```julia-repl
julia> a=Z(19)+3
FFE{19}: 5

julia> Mod(a)
Mod{UInt64}: 5₁₉

julia> Int(a)
5

julia> order(a) # order as element of the multiplicative group
9
```

The field, `p`, `n` and `p^n` can be obtained back from an `FFE{p}` as well
as which power of `Z(p^n)` is considered

```julia-repl
julia> a=Z(8)^5
FFE{2}: Z₈⁵

julia> F=field(a)
GF(2^3)

julia> char(F)
2

julia> char(a)
2

julia> degree(F)
3

julia> degree(a)
3

julia> length(F)
8

julia> log(a)
5

julia> elements(F)
8-element Vector{FFE{2}}:
   0
   1
  Z₈
 Z₈²
 Z₈³
 Z₈⁴
 Z₈⁵
 Z₈⁶
```

A  `p`-integral integer or  rational or a  `Mod(,p)` can be  converted to a
prime field element using `FFE{p}` as a constructor.

```julia-repl
julia> FFE{19}(2)
FFE{19}: 2

julia> FFE{19}(5//3)
FFE{19}: 8

julia> FFE{19}(Mod(2,19))
FFE{19}: 2
```
```julia-rep1
julia> m=rand(GF(49),4,4)
4×4 Matrix{FFE{7}}:
 Z₄₉²⁴  Z₄₉¹⁸   Z₄₉⁹  Z₄₉⁴²
 Z₄₉²²  Z₄₉⁴¹  Z₄₉⁴⁶  Z₄₉²⁴
 Z₄₉¹⁵  Z₄₉¹⁹  Z₄₉⁴⁰   Z₄₉³
 Z₄₉²⁰  Z₄₉²⁹  Z₄₉³⁶  Z₄₉²⁰

julia> inv(m)
4×4 Matrix{FFE{7}}:
 Z₄₉³⁷   Z₄₉⁵  Z₄₉³⁶      1
 Z₄₉¹⁰    Z₄₉   Z₄₉⁶  Z₄₉⁴⁷
 Z₄₉³⁰  Z₄₉³⁸    Z₄₉     -2
 Z₄₉¹⁵   Z₄₉²      1  Z₄₉²⁸

julia> inv(m)*m
4×4 Matrix{FFE{7}}:
 1  0  0  0
 0  1  0  0
 0  0  1  0
 0  0  0  1
```
"""
module FiniteFields
include("Modulo.jl")
using .Modulo
export Modulo, Mod, order
using Primes: factor, eachfactor, totient
# next can be used from Primes when exported there
divisors(n::Integer)=vec(map(prod,Iterators.product((p.^(0:m) 
                             for (p,m) in eachfactor(n))...)))::Vector{Int}

export GF, FFE, Z, field, degree, char, elements

const conway_polynomials=Dict{Tuple{Int,Int},Vector{Int}}(
(2, 1) => [1],
(2, 2) => [1, 1],
(2, 3) => [1, 1, 0],
(2, 4) => [1, 1, 0, 0],
(2, 5) => [1, 0, 1, 0, 0],
(2, 6) => [1, 1, 0, 1, 1, 0],
(2, 7) => [1, 1, 0, 0, 0, 0, 0],
(2, 8) => [1, 0, 1, 1, 1, 0, 0, 0],
(2, 9) => [1, 0, 0, 0, 1, 0, 0, 0, 0],
(2, 10) => [1, 1, 1, 1, 0, 1, 1, 0, 0, 0],
(2, 11) => [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
(2, 12) => [1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0],
(2, 13) => [1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
(2, 14) => [1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
(2, 15) => [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
(3, 1) => [1],
(3, 2) => [2, 2],
(3, 3) => [1, 2, 0],
(3, 4) => [2, 0, 0, 2],
(3, 5) => [1, 2, 0, 0, 0],
(3, 6) => [2, 2, 1, 0, 2, 0],
(3, 7) => [1, 0, 2, 0, 0, 0, 0],
(3, 8) => [2, 2, 2, 0, 1, 2, 0, 0],
(3, 9) => [1, 1, 2, 2, 0, 0, 0, 0, 0],
(3, 10) => [2, 1, 0, 0, 2, 2, 2, 0, 0, 0],
(5, 1) => [3],
(5, 2) => [2, 4],
(5, 3) => [3, 3, 0],
(5, 4) => [2, 4, 4, 0],
(5, 5) => [3, 4, 0, 0, 0],
(5, 6) => [2, 0, 1, 4, 1, 0],
(7, 1) => [4],
(7, 2) => [3, 6],
(7, 3) => [4, 0, 6],
(7, 4) => [3, 4, 5, 0],
(7, 5) => [4, 1, 0, 0, 0],
(11, 1) => [9],
(11, 2) => [2, 7],
(11, 3) => [9, 2, 0],
(11, 4) => [2, 10, 8, 0],
(13, 1) => [11],
(13, 2) => [2, 12],
(13, 3) => [11, 2, 0],
(13, 4) => [2, 12, 3, 0],
(17, 1) => [14], (17, 2) => [3, 16], (17, 3) => [14, 1, 0],
(19, 1) => [17], (19, 2) => [2, 18], (19, 3) => [17, 4, 0],
(23, 1) => [18], (23, 2) => [5, 21], (23, 3) => [18, 2, 0],
(29, 1) => [27], (29, 2) => [2, 24], (29, 3) => [27, 2, 0],
(31, 1) => [28], (31, 2) => [3, 29], (31, 3) => [28, 1, 0],
(37, 1) => [35], (37, 2) => [2, 33], (37, 3) => [35, 6, 0],
(41, 1) => [35], (41, 2) => [6, 38],
(43, 1) => [40], (43, 2) => [3, 42],
(47, 1) => [42], (47, 2) => [5, 45],
(53, 1) => [51], (53, 2) => [2, 49],
(59, 1) => [57], (59, 2) => [2, 58],
(61, 1) => [59], (61, 2) => [2, 60],
(67, 1) => [65], (67, 2) => [2, 63],
(71, 1) => [64], (71, 2) => [7, 69],
(73, 1) => [68], (73, 2) => [5, 70],
(79, 1) => [76], (79, 2) => [3, 78],
(83, 1) => [81], (83, 2) => [2, 82],
(89, 1) => [86], (89, 2) => [3, 82],
(97, 1) => [92], (97, 2) => [5, 96],
(101, 1) => [99], (101, 2) => [2, 97],
(103, 1) => [98], (103, 2) => [5, 102],
(107, 1) => [105], (107, 2) => [2, 103],
(109, 1) => [103], (109, 2) => [6, 108],
(113, 1) => [110], (113, 2) => [3, 101],
(127, 1) => [124], (127, 2) => [3, 126],
(131, 1) => [129], (131, 2) => [2, 127],
(137, 1) => [134], (137, 2) => [3, 131],
(139, 1) => [137], (139, 2) => [2, 138],
(149, 1) => [147], (149, 2) => [2, 145],
(151, 1) => [145], (151, 2) => [6, 149],
(157, 1) => [152], (157, 2) => [5, 152],
(163, 1) => [161], (163, 2) => [2, 159],
(167, 1) => [162], (167, 2) => [5, 166],
(173, 1) => [171], (173, 2) => [2, 169],
(179, 1) => [177], (179, 2) => [2, 172],
(181, 1) => [179], (181, 2) => [2, 177],
(191, 1) => [172], (191, 2) => [19, 190],
(193, 1) => [188], (193, 2) => [5, 192],
(197, 1) => [195], (197, 2) => [2, 192],
(199, 1) => [196], (199, 2) => [3, 193],
(211, 1) => [209], (211, 2) => [2, 207],
(223, 1) => [220], (223, 2) => [3, 221],
(227, 1) => [225], (227, 2) => [2, 220],
(229, 1) => [223], (229, 2) => [6, 228],
(233, 1) => [230], (233, 2) => [3, 232],
(239, 1) => [232], (239, 2) => [7, 237],
(241, 1) => [234], (241, 2) => [7, 238],
(251, 1) => [245], (251, 2) => [6, 242],
)

function primitiveroot(p::Integer)
  phi=totient(p) # the Euler φ
  1+findfirst(x->powermod(x,phi,p)==1 && 
     all(d->powermod(x,div(phi,d),p)!=1,keys(factor(phi))),2:p-1)::Int
end

function conway_polynomial(p,n=1)
  get!(conway_polynomials,(p,n)) do
    if n==1 return [p-primitiveroot(p)] end
    error("missing conway($p,$n)") 
  end
end

struct GF # finite field with q=p^d elements
  p::Int # characteristic
  d::Int # degree
  q::Int # cardinality
  conway::Vector{UInt16}# -[coefficients of the Conway polynomial except leading]
  zech::Vector{UInt16} # prime field: log
  dict::Vector{UInt16} # prime field: dict[i+1]=Z(q)^i others: degree
end

char(x)=0
@inline char(F::GF)=F.p
@inline degree(F::GF)=F.d
@inline Base.length(F::GF)=F.q

"""
`FFE{p}` is the type of the elements of a finite field of characteristic `p`.
"""
struct FFE{p}<:Number
  i::UInt16  # represents Z(field.q)^i
  Fi::UInt16 # where field=FFvec[Fi]
end

@inline prime(x::FFE)=iszero(x.Fi)
@inline field(x::FFE{p}) where p=@inbounds FFvec[prime(x) ? iFF(p) : x.Fi]
@inline char(x::FFE{p}) where p=p
@inline char(::Type{FFE{p}}) where p=p
# the Zech logarithm of x=Z(q)^i is defined by 1+x=Z(q)^(zech(x,i))
@inline zech(x::FFE,i)=@inbounds field(x).zech[1+i]
@inline sizefield(x::FFE{p}) where p=prime(x) ? p : @inbounds FFvec[x.Fi].q
@inline degree(x::FFE{p}) where p=prime(x) ? 1 : @inbounds FFvec[x.Fi].d

const FFi=Dict{Int,Int}()
const FFvec=GF[]

function iFF(q) # get index of 𝔽_q in FFvec from q (using FFi Dict)
  get!(FFi,q) do
    if q>2^16 error(q,">2^16") end
    l=factor(q)
    if length(l)>1 error(q," should be a prime power") end
    p,n=first(l)
    pol=-Mod.(conway_polynomial(p,n),p)
    if n==1
      val=pol[1].val
#     println(val)
      dic=Vector{UInt16}(undef,p)
      c=1
      for i in 0:p-1
        dic[i+1]=c
        c=mod(c*val,p)
      end
      zz=Vector{UInt16}(undef,p-1)
      c=1
      for i in 0:p-2
        zz[c]=i
        c=mod(c*val,p)
      end
    else
      l=sort(divisors(n),rev=true)
      dic=fill(UInt16(n),q)
      for i in 2:length(l)
        dic[1:div(q-1,p^l[i]-1):q].=l[i]
      end
  #   xprint("conway=",pol)
      zz=map(i->fill(Mod(0,p),n),1:q)
      zz[1][1]=Mod(1,p)
      for i in 2:q-1
        for j in 2:n zz[i][j]=zz[i-1][j-1] end
        zz[i].+= zz[i-1][n].*pol
      end
      z=collect(pairs(reverse.(zz)))
      sort!(z,by=x->x[2])
  #   for i in 1:q xprint(i,"->",(z[i][1]-1,z[i][2])) end
      zz=Vector{UInt16}(undef,q)
      for i in 1:q
        if z[i][2][end]==Mod(p-1,p) zz[z[i][1]]=z[i-p+1][1]-1
        else zz[z[i][1]]=z[i+1][1]-1
        end
      end
  #   for i in 0:q-1 xprint(i,"=>",zz[i+1]) end
    end
    push!(FFvec,GF(p,n,q,pol,zz,dic))
    length(FFvec)
  end
end

@inbounds GF(q)=FFvec[iFF(q)]

Base.show(io::IO,F::GF)=print(io,"GF(",char(F),"^",degree(F),")")

@inline clone(x::FFE{p},i) where p=FFE{p}(i,x.Fi)

Base.iszero(x::FFE{p}) where p=x.Fi==0 && x.i==0

Base.:^(a::FFE{p},n::Integer) where p=prime(a) ? 
  clone(a,powermod(a.i,n,p)) : lower(clone(a,mod(widen(a.i)*n,sizefield(a)-1)))

Base.rand(F::GF)=FFE{F.p}(rand(0:F.q-1),F.d==1 ? 0 : iFF(F.q))
Base.rand(F::GF,n::Int...)=[rand(F) for i in fill(0,n...)]
Base.rand(::Type{FFE{p}},n::Integer...) where p=FFE{p}.(rand(0:p-1,n...),0)

"""
`Z(p^d)`

returns  a  generator  of  the  multiplicative  group  of  the finite field
`𝔽_{pᵈ}`,  where    `p`  must  be  prime  and `pᵈ` smaller than `2¹⁵`. This
multiplicative  group  is  cyclic  thus  `Z(pᵈ)ᵃ`  runs  over it for `a` in
`0:pᵈ-1`.  The zero  of the  field is  `0*Z(p)` (the  same as `0*Z(pᵈ)`; we
automatically lower an element to the smallest field which contains it).

The  various generators returned by `Z` for finite fields of characteristic
`p`  are compatible. That  is, if the  field `𝔽_{pⁿ}` is  a subfield of the
field `𝔽_{pᵐ}`, that is, `n` divides `m`, then
`Z(pⁿ)=Z(pᵐ)^div(pᵐ-1,pⁿ-1)`.  This is  achieved by  choosing `Z(p)` as the
smallest  primitive root  modulo `p`  and `Z(pⁿ)`  as a  root of the `n`-th
Conway polynomial of characteristic `p`. Those polynomials where defined by
J.H.~Conway and computed by R.A.~Parker.

```julia-repl
julia> z=Z(16)
FFE{2}: Z₁₆

julia> z^5
FFE{2}: Z₄
```
"""
function Z(q)
  Fi=iFF(q)
  F=FFvec[Fi]
  if F.d==1 FFE{F.p}(F.dict[2],0)
  else FFE{F.p}(1,Fi)
  end
end

function elements(F::GF)
  l=Z(length(F)).^(0:length(F)-2)
  pushfirst!(l,0*Z(length(F)))
  l
end

Base.promote_rule(::Type{FFE{p}},::Type{<:Integer}) where p=FFE{p}
Base.promote_rule(::Type{FFE{p}},::Type{<:Rational{<:Integer}}) where p=FFE{p}

Base.copy(a::FFE{p}) where p=clone(a,a.i)
Base.one(::Type{FFE{p}}) where p=FFE{p}(1,0)
Base.one(x::FFE{p}) where p=FFE{p}(1,0)
Base.zero(::Type{FFE{p}}) where p=FFE{p}(0,0)
Base.zero(x::FFE{p}) where p=FFE{p}(0,0)
Base.abs(a::FFE)=a
Base.conj(a::FFE)=a
Base.gcd(x::FFE,y::FFE)=one(x)
Base.gcd(x::AbstractVector{<:FFE})=one(x[1])

function Base.cmp(x::FFE, y::FFE)
  l=cmp(sizefield(x),sizefield(y))
  if !iszero(l) return l end
  cmp(x.i,y.i)
end

Base.isless(x::FFE, y::FFE)=cmp(x,y)==-1
Base.:(==)(x::FFE, y::FFE)=x.i==y.i && x.Fi==y.Fi

function Base.show(io::IO, ::MIME"text/plain", x::FFE{p})where p
  if !haskey(io,:typeinfo) print(io,typeof(x),": ") end
  show(io,x)
end

function Base.show(io::IO, x::FFE{p})where p
  sup=Dict(zip("0123456789","⁰¹²³⁴⁵⁶⁷⁸⁹"))
  sub=Dict(zip("0123456789,()","₀₁₂₃₄₅₆₇₈₉‚₍₎"))
  if !get(io,:limit,false)
    if iszero(x) print(io,"0*Z(",char(x),")")
    else print(io,"Z(",char(x))
      if degree(x)!=1 print(io,"^",degree(x)) end
      print(io,")")
      if log(x)!=1 print(io,"^",log(x)) end
    end
    return
  end
  if sizefield(x)==p
    print(io,Int(x))
  else 
    print(io,"Z",map(x->sub[x],repr(sizefield(x))))
    if x.i!=1
      print(io,map(x->sup[x],repr(Int(x.i))))
    end
  end
end

# convert i∈ 0:p-1 to a prime FFE
function primeFFE(p,i) # assumes i is already modp and p prime
  Fi=iFF(p)
  if FFvec[Fi].d!=1 error("$p should be prime") end
  FFE{p}(i,0)
end

"""
`FFE{p}(i)`  for `i` an integer or a fraction with denominator prime to `p`
returns the reduction mod `p` of `i`, an element of the prime field `𝔽ₚ`.
"""
FFE{p}(i::Integer) where p=primeFFE(p,mod(i,p))
FFE{p}(i::Mod) where p=primeFFE(Int(i.n),i.val)

function FFE{p}(i::Rational{<:Integer})where p
  primeFFE(p,mod(invmod(denominator(i),p)*numerator(i),p))
end

FFE{p}(x::FFE{p}) where p=x

function (::Type{Mod})(x::FFE{p})where p
  if !prime(x) throw(InexactError(:convert,Mod,x)) end
  Mod(x.i,p)
end

function (::Type{T})(x::FFE{p})where {T<:Integer,p}
  if !prime(x) throw(InexactError(:convert,T,x)) end
  T(Mod(x.i,p))
end

function Base.log(x::FFE{p})where p
  if prime(x)
    if iszero(x) error(x," must be invertible") end
    return Int(GF(p).zech[x.i])
  end
  Int(x.i)
end

@inline function promote_field(iF::Integer,x::FFE{p})where {p}
  q=sizefield(x)
  @inbounds nq=FFvec[iF].q
  FFE{p}(log(x)*div(nq-1,q-1),iF)
end
  
function promote_field(x::FFE{p},y::FFE{p})where p
  dx=degree(x);dy=degree(y)
  d=lcm(dx,dy)
  if d==dx 
    (x,promote_field(x.Fi,y))
  elseif d==dy 
    (promote_field(y.Fi,x),y)
  else iF=iFF(p^d)
    (promote_field(iF,x),promote_field(iF,y))
  end
end

function Base.:+(x::FFE{p},y::FFE{p}) where p
  if iszero(x) return y end
  if iszero(y) return x end
  if x.Fi!=y.Fi x,y=promote_field(x,y) end
  if prime(x) 
    res=widen(x.i)+y.i
    if res>=p res-=p end
    return clone(x,res) end
  if x.i>y.i x,y=y,x end
  res=zech(x,y.i-x.i)
  q1=sizefield(x)-1
  if res==q1 return zero(x) end
  res+=widen(x.i)
  if res>=q1 res-=q1 end
  lower(clone(x,res))
end

function Base.:-(x::FFE{p})where p
  if prime(x) 
    if iszero(x) return x end
    return clone(x,p-x.i)
  end
  if iseven(p) return x end
  q1=sizefield(x)-1
  res=x.i+div(q1,2)
  if res>=q1 res-=q1 end
  clone(x,res)
end

Base.:-(x::FFE,y::FFE)=x+(-y)

function Base.:*(x::FFE{p},y::FFE{p}) where p
  if iszero(x) return x end
  if iszero(y) return y end
  if x.Fi!=y.Fi x,y=promote_field(x,y) end
  if prime(x) return clone(x,rem(widemul(x.i,y.i),p)) end
  res=x.i+Int(y.i)
  q1=sizefield(x)-1
  if res>=q1 res-=q1 end
  lower(clone(x,res))
end

function Base.inv(x::FFE{p}) where p
  if prime(x) return clone(x,invmod(x.i,p)) end
  if x.i==0 return x end
  clone(x,sizefield(x)-1-x.i)
end

Base.:/(x::FFE,y::FFE)=x*inv(y)
Base.://(x::FFE,y::FFE)=x*inv(y)

function lower(x::FFE{p}) where p # x is never prime here
  F=field(x)
@inbounds actual_degree=F.dict[1+x.i]
  if actual_degree==F.d return x end
  newi=div(x.i,div(sizefield(x)-1,p^actual_degree-1))
  if actual_degree==1 FFE{p}(GF(p).dict[newi+1],0)
  else FFE{p}(newi,iFF(p^actual_degree))
  end
end

function Modulo.order(x::FFE{p}) where p
  if iszero(x) error(x," must be invertible") end
  div(sizefield(x)-1,gcd(log(x),sizefield(x)-1))
end

end
