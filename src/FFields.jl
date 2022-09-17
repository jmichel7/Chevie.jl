"""
This module introduces modular arithmetic and finite fields.

The  integer `x` mod. `n` is constructed  by the function `Mod(x,n)`. If `n
isa  Int`  its  type  is  `Mod{UInt64}`.  If  `n  isa  BigInt`  its type is
`Mod{BigInt}`.  Since `n` is not encoded in  the type, the elements `0` and
`1`  mod.  `n`  cannot  be  constructed  from  the  type, which causes some
problems  for some  Julia functions.  For some  prime moduli  `p`, the type
`FFE{p}` below does not have such limitations.

Example:
```julia-repl
julia> Mod(5,19)
Mod{UInt64}: 5₁₉

julia> a=Mod(5,19)
Mod{UInt64}: 5₁₉

julia> a^2
Mod{UInt64}: 6₁₉

julia> inv(a)
Mod{UInt64}: 4₁₉

julia> a*inv(a)
Mod{UInt64}: 1₁₉

julia> a+2
Mod{UInt64}: 7₁₉

julia> a*2
Mod{UInt64}: -9₁₉

julia> a+1//2
Mod{UInt64}: -4₁₉

julia> Integer(a) # get back an integer from a
5

julia> order(a) # multiplicative order of a
9
```

The Galois field with `p^n` elements is obtained as `GF(p^n)`. To work with
elements  of  this  field,  (as  in  GAP)  the  function `Z(p^n)` returns a
generator  of  the  multiplicative  group  of  `GF(p^n)`. Other elements of
`GF(p^n)` are obtained as a power of `Z(p^n)` or as `0*Z(p^n)`.

```julia-repl
julia> a=Z(64)
FFE{2}: Z₆₄

julia> a^9
FFE{2}: Z₈

julia> a^21
FFE{2}: Z₄

julia> a+1
FFE{2}: Z₆₄⁵⁶
```

Elements  of the prime field can be converted back to integers `Mod(,p)` or
to integers:

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

The  field, `p`, `n` and `p^n` can be  obtained back as well as which power
of `Z(p^n)` is considered

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

The type of an element of `GF(p^n)` is `FFE{p}`. A `p`-integral number or a
`Mod(,p)` can be converted to the prime field using this type as constructor.

```julia-repl
julia> FFE{19}(2)
FFE{19}: 2

julia> FFE{19}(Mod(2,19))
FFE{19}: 2
```
"""
module FFields
using ..Combinat: divisors
using Primes: factor
export Mod, GF, FFE, Z, field, order, degree, char, elements

struct Mod{T}<:Number
  val::T
  n::T
  global Mod_(a::T,n::T) where T=new{T}(a,n)
end

function Mod(a::Integer,n)
  if n isa BigInt return Mod_(mod(a,n),n) end
  Mod_(unsigned(mod(a,n)),unsigned(n))
end

Mod(i::Rational{<:Integer},p)=Mod(numerator(i),p)/Mod(denominator(i),p)
Base.promote(a::Integer,b::Mod)=(Mod(a,b.n),b)
Base.promote(b::Mod,a::Integer)=(b,Mod(a,b.n))
Base.promote(a::Rational{<:Integer},b::Mod)=(Mod(a,b.n),b)
Base.promote(b::Mod,a::Rational{<:Integer})=(b,Mod(a,b.n))
Base.promote(a,b::Mod)=(a*one(b),one(a)*b)
Base.zero(n::Mod)=Mod(0,n.n)
Base.one(n::Mod)=Mod(1,n.n)
Base.one(::Type{<:Mod})=Mod_(unsigned(1),unsigned(0))
Base.isone(n::Mod)=isone(n.val)
Base.:(==)(x::Mod,y::Mod)=x.n==y.n && x.val==y.val
Base.:(==)(x::Mod,y::Integer)=false
Base.:+(x::Mod, y::Mod)=x.n!=y.n ? error("moduli") : Mod(x.val+y.val,x.n)
Base.:*(x::Mod, y::Mod)=x.n!=y.n ? error("moduli") : Mod(x.val*y.val,x.n)
Base.:-(x::Mod)=Mod(-signed(Integer(x.val)),x.n)
Base.:-(x::Mod, y::Mod)=x+(-y)
Base.:/(x::Mod,y::Mod)=x*inv(y)
Base.:/(x::Mod,y::Number)=x*inv(Mod(y,x.n))
Base.inv(x::Mod)=Mod(invmod(x.val,x.n),x.n)
Base.:^(x::Mod,m::Integer)=m>=0 ? Base.power_by_squaring(x,m) :
                                  Base.power_by_squaring(inv(x),-m)
Base.cmp(x::Mod,y::Mod)=cmp(x.val,y.val)
Base.gcd(x::Mod,y::Mod)=Mod(gcd(x.val,y.val),x.n)
Base.isless(x::Mod,y::Mod)=cmp(x,y)==-1
Base.abs(x::Mod)=x      # needed for inv(Matrix) to work
Base.conj(x::Mod)=x     # needed for inv(Matrix) to work
(::Type{T})(x::Mod) where T<:Integer=2x.val<=x.n ? signed(T(x.val)) : signed(T(x.val))-signed(T(x.n))

function Base.show(io::IO, ::MIME"text/plain", x::Mod)
  if !haskey(io,:typeinfo) print(io,typeof(x),": ") end
  show(io,x)
end

const sub=Dict(zip("-0123456789,+()=aehijklmnoprstuvxβγρφχ.",
                   "₋₀₁₂₃₄₅₆₇₈₉‚₊₍₎₌ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵦᵧᵨᵩᵪ̣."))

function Base.show(io::IO, m::Mod) where p 
  if get(io,:limit,false) print(io,Integer(m),map(x->sub[x],string(m.n)))
  else print(io,"Mod(",Integer(m),",",m.n,")")
  end
end

function order(x::Mod)
  d=gcd(x.val,x.n)
  if d!=1 error("$x should be invertible") end
  o=1
  res=x
  while true
   if res.val<=1 return o end
   o+=1
   res*=x
  end
end

#-------------------------------------------------------------------------
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

function conway_polynomial(p,n)
  get!(conway_polynomials,(p,n)) do
    if n==1 
      x=2
      while x<=p-1
        if order(Mod(x,p))==p-1 return [p-x] end
        x+=1
      end
    end
    error("missing conway($p,$n)") 
  end
end

struct GF # finite field with q=p^n elements
  p::Int16 # characteristic
  d::Int16 # degree
  q::Int16 # cardinality
  conway::Vector
  zech::Vector{Int16} 
  dict::Vector{Int16}
end

@inline char(F::GF)=Int(F.p)
@inline degree(F::GF)=F.d
@inline Base.length(F::GF)=F.q

"""
`FFE{p}` is the type of the elements of a finite field of characteristic `p`.
"""
struct FFE{p}<:Number
  i::Int16  # the number is Z(field.q)^i
  Fi::Int16 # the field is FFvec[Fi]
end

@inbounds @inline field(x::FFE)=FFvec[x.Fi]
@inline char(x::FFE{p}) where p=p

#Base.:(==)(::Type{FFE{p}},::Type{FFE{q}}) where{p,q}=p==q

function Base.promote_rule(::Type{FFE{p}},::Type{FFE{q}})where {p,q}
  if p!=q error("cannot mix characteristics $p and $q") end
  FFE{p}
end

const FFi=Dict{Int,Int}()
const FFvec=GF[]

function iFF(q)
  get!(FFi,q) do
    l=factor(q)
    if length(l)>1 error(q," should be a prime power") end
    p,n=first(l)
    pol=-Mod.(conway_polynomial(p,n),p)
    if n==1
      dic=Vector{Int16}(undef,q)
      val=pol[1].val
#     println(val)
      dic[1]=-1
      dic[2]=0
      c=1
      for i in 1:p-2
        c=mod(c*val,p)
        dic[c+1]=i
      end
    else
      l=reverse(divisors(n))
      dic=fill(Int16(n),q)
      for i in 2:length(l)
        dic[1:div(q-1,p^l[i]-1):q].=l[i]
      end
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
    zz=Vector{Int16}(undef,q)
    for i in 1:q
      if z[i][2][end]==Mod(p-1,p) zz[z[i][1]]=z[i-p+1][1]-1
      else zz[z[i][1]]=z[i+1][1]-1
      end
    end
#   for i in 0:q-1 xprint(i,"=>",zz[i+1]) end
    push!(FFvec,GF(p,n,q,pol,zz,dic))
    return length(FFvec)
  end
end

GF(q)=FFvec[iFF(q)]

Base.show(io::IO,F::GF)=print(io,"GF(",char(F),"^",degree(F),")")

@inline sizefield(x::FFE)=field(x).q
@inline degree(x::FFE)=degree(field(x))
# the Zech logarithm of x=Z(q)^i is defined by 1+x=Z(q)^(zech(x,i))
@inline zech(x::FFE,i)=@inbounds field(x).zech[1+i]
@inline dict(x::FFE,i)=@inbounds field(x).dict[1+i]
@inline clone(x::FFE{p},i) where p=FFE{p}(i,x.Fi)

Base.iszero(x::FFE{p}) where p=x.i==-1

Base.:^(a::FFE{p},n::Integer) where p=iszero(a) ? a : 
                 lower(clone(a,mod(a.i*n,sizefield(a)-1)))

Base.rand(F::GF)=Z(F.p,F.d)^rand(0:F.q-2)
Base.rand(F::GF,n::Int)=map(i->rand(F),1:n)

"""
`Z(p^d)` or `Z(p,d)`

This  returns a generator  of the multiplicative  group of the finite field
`𝔽_{pᵈ}`,  where   `p`  must be  prime and  `p^d` smaller  than `2¹⁵`. This
multiplicative  group is  cyclic thus  `Z(p^d)^a` runs  over it  for `a` in
`0:p^d-1`.  The zero of the  field is `0*Z(p)` (the  same as `0*Z(p^d)`; we
automatically lower an element to the smallest field which contains it).

The  various generators returned by `Z` for finite fields of characteristic
`p`  are compatible. That  is, if the  field `𝔽_{pⁿ}` is  a subfield of the
field `𝔽_{pᵐ}`, that is, `n` divides `m`, then
`Z(p^n)=Z(p^m)^div(p^m-1,p^n-1)`.  This is  achieved by  choosing `Z(p)` as
the smallest primitive root modulo `p` and `Z(p^n)` as a root of the `n`-th
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
  if q==2 return FFE{2}(0,iFF(2)) end
  Fi=iFF(q)
  FFE{Int(char(FFvec[Fi]))}(1,Fi)
end

function Z(p,n)
  if p==2 && n==1 return FFE{2}(0,iFF(2)) end
  FFE{Int(p)}(1,iFF(p^n))
end

function elements(F::GF)
  l=Z(length(F)).^(0:length(F)-2)
  pushfirst!(l,0*Z(length(F)))
  l
end

Base.promote_rule(::Type{FFE{p}},::Type{<:Integer}) where p=FFE{p}
Base.promote_rule(::Type{FFE{p}},::Type{<:Rational{<:Integer}}) where p=FFE{p}

Base.copy(a::FFE{p}) where p=clone(a,a.i)
Base.one(::Type{FFE{p}}) where p=FFE{p}(0,iFF(p))
Base.one(x::FFE{p}) where p=sizefield(x)==p ? clone(x,0) : one(FFE{p})
Base.zero(::Type{FFE{p}}) where p=FFE{p}(-1,iFF(p))
Base.zero(x::FFE{p}) where p=sizefield(x)==p ? clone(x,-1) : zero(FFE{p})
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
    if iszero(x) print(io,"0*Z(",sizefield(x),")")
    else print(io,"Z(",sizefield(x),")^",x.i)
    end
    return
  end
  if sizefield(x)==p
    if iszero(x) print(io,"0")
    else n=(Mod(field(x).conway[1].val,p)^x.i).val
      if 2n<=p print(io,n)
      else print(io,Int(n)-p)
      end
    end
  else print(io,"Z",map(x->sub[x],repr(sizefield(x))))
    if x.i!=1
      print(io,map(x->sup[x],repr(x.i)))
    end
  end
end

function Zi(p,i) # assume i is already modp
  Fi=iFF(p)
  FFE{Int(p)}(FFvec[Fi].dict[1+i],Fi)
end

"""
`FFE{p}(i)`  for `i` an integer or a fraction with denominator prime to `p`
returns the reduction mod `p` of `i`, an element of the prime field `𝔽ₚ`.
"""
FFE{p}(i::Integer) where p=Zi(p,mod(i,p))
FFE{p}(i::Mod) where p=Zi(i.n,i.val)

function FFE{p}(i::Rational{<:Integer})where p
  i=mod(invmod(denominator(i),p)*numerator(i),p)
  Zi(p,i)
end

FFE{p}(x::FFE{q}) where {p,q}=p!=q ? error("cannot mix chars $p and $q") : x

function (::Type{Mod})(x::FFE{p})where {p}
  if sizefield(x)!=p throw(InexactError(:convert,Mod,x)) end
  iszero(x) ? Mod(0,p) : Mod(field(x).conway[1].val,p)^x.i
end

function (::Type{T})(x::FFE{p})where {T<:Integer,p}
  if sizefield(x)!=p throw(InexactError(:convert,T,x)) end
  iszero(x) ? T(0) : T(Mod(field(x).conway[1].val,p)^x.i)
end

function promote_field(r::Integer,x::FFE{p})where {p}
  if r==sizefield(x) return x end
  if r<sizefield(x) error("cannot promote to smaller field") end
  Fi=iFF(r)
@inbounds  F=FFvec[Fi]
  if char(F)!=p error("cannot mix characteristics $p and ",char(F)) end
  if mod(degree(F),degree(x))!=0 error("not an extension field") end
  FFE{p}(x.i*div(r-1,sizefield(x)-1),Fi)
end
  
function Base.promote(x::FFE{p},y::FFE{q})where {p,q}
  if p!=q error("cannot mix characteristics $p and $q") end
  if x.Fi==y.Fi return (x,y) end
  nq=p^lcm(degree(x),degree(y))
  (promote_field(nq,x),promote_field(nq,y))
end

function Base.:+(x::FFE{p},y::FFE{p}) where p
  if iszero(x) return y
  elseif iszero(y) return x
  end
  x,y=promote(x,y)
  if x.i>y.i x,y=y,x end
  @inbounds res=zech(x,y.i-x.i)
  if res==sizefield(x)-1 return zero(x) end
  res+=x.i
  if res>=sizefield(x)-1 res-=Int16(sizefield(x)-1) end
  lower(clone(x,res))
end

function Base.:-(x::FFE{p})where p
  if iseven(p) || iszero(x) return x end
  res=x.i+div(sizefield(x)-1,2)
  if res>=sizefield(x)-1 res-=Int16(sizefield(x)-1) end
  clone(x,res)
end

Base.:-(x::FFE,y::FFE)=x+(-y)

function Base.:*(x::FFE{p},y::FFE{p}) where p
  if iszero(x) return x
  elseif iszero(y) return y
  end
  (x,y)=promote(x,y)
  res=x.i+y.i
  if res>=sizefield(x)-1 res-=Int16(sizefield(x)-1) end
  lower(clone(x,res))
end

function Base.inv(x::FFE{p}) where p
  if iszero(x) error("inv(0)") end
  if x.i==0 return x end
  clone(x,sizefield(x)-1-x.i)
end

Base.:/(x::FFE,y::FFE)=x*inv(y)
Base.://(x::FFE,y::FFE)=x*inv(y)

function lower(x::FFE{p}) where p
  if sizefield(x)==p return x end
  l=dict(x,x.i)
  if l==degree(x) return x end
  FFE{p}(div(x.i,div(sizefield(x)-1,p^l-1)),iFF(p^l))
end

function order(x::FFE{p}) where p
  if iszero(x) error(x," must be invertible") end
  div(sizefield(x)-1,gcd(x.i,sizefield(x)-1))
end

function Base.log(x::FFE)
  if iszero(x) error(x," must be invertible") end
  x.i
end

# random elements of the prime field
Base.rand(::Type{FFE{p}},n::Integer...) where p=FFE{p}.(rand(1:p,n...))

end
