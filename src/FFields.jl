"""
This module introduces modular arithmetic and finite fields.

The ring of integers mod. `n` is given by the type `Mod{n}`.

Example:
```julia-repl
julia> a=Mod{19}(5)
Mod{19}(5)

julia> a^2
Mod{19}(6)

julia> inv(a)
Mod{19}(4)

julia> a*inv(a)
Mod{19}(1)

julia> a+2
Mod{19}(7)

julia> a*2
Mod{19}(10)

julia> a+1//2
Mod{19}(15)

julia> Int(a) # get back an integer from a
5

julia> order(a) # multiplicative order of a
9
```

The finite field with `p^n` elements is obtained as `FF(p^n)`. To work with
elements  of this  field, the  function `Z(p^n)`  returns a genertor of the
multiplicative group of `FF(p^n)` (this is a particular generator, obtained
as a root of the `n`-th Conway polynomial of characteristic `p`). All other
elements  of  `FF(p^n)`  can  be  obtained  as  a  power  of `Z(p^n)` or as
`0*Z(p^n)`.

```julia-repl
julia> a=Z(64)
Z₆₄

julia> a^9
Z₈

julia> a^21
Z₄

julia> a+1
Z₆₄⁵⁶
```

Elements of the prime field can be converted back to integers `Mod{p}`

```julia-repl
julia> a=Z(19)+3
Z₁₉¹⁶

julia> Mod{19}(a)
Mod{19}(5)

julia> order(a) # order as element of the multiplicative group
9
```

The  field, `p`, `n` and `p^n` can be  obtained back as well as which power
of `Z(p^n)` is considered

```julia-repl
julia> a=Z(8)^5
Z₈⁵

julia> F=field(a)
FF(2^3)

julia> F.p
2

julia> F.n
3

julia> F.q
8

julia> log(a)
5
```

The type of an element of `FF(p^n)` is `FFE{p}`. An integer or a `Mod{p}` can
be converted to the prime field using this type as constructor.

```julia-repl
julia> FFE{19}(2)
Z₁₉

julia> FFE{19}(Mod{19}(2))
Z₁₉
```
"""
module FFields
using ..Util: factor, divisors
export Mod, FF, FFE, Z, field, order

const T=UInt8 # this allows moduli up to 255

struct Mod{p} <: Number
   val::T
   function Mod{p}(a::Integer) where {p}
   	new(T(mod(a,p)))
   end
end

Base.promote_rule(::Type{Mod{p}},::Type{<:Integer}) where p=Mod{p}
Base.promote_rule(::Type{Mod{p}},::Type{<:Rational{<:Integer}}) where p=Mod{p}
Mod{p}(i::Rational{<:Integer}) where p=Mod{p}(invmod(denominator(i),p)*numerator(i))
Base.zero(::Type{Mod{p}}) where p = Mod{p}(0)
Base.:(==)(x::Mod,y::Mod) = x.val==y.val
Base.:+(x::Mod{p}, y::Mod{p}) where p = Mod{p}(Int(x.val)+y.val)
Base.:*(x::Mod{p}, y::Mod{p}) where p = Mod{p}(Int(x.val)*y.val)
Base.:-(x::Mod{p}, y::Mod{p}) where p = Mod{p}(Int(x.val)-y.val)
Base.:-(x::Mod{p}) where {p} = Mod{p}(-Int(x.val))
Base.:/(x::Mod{p}, y::Mod{p}) where p = x * inv(y)
Base.inv(x::Mod{p}) where p = Mod{p}(invmod(x.val,p))
Base.:^(x::Mod{p},m::Integer) where p=(m>=0) ? Base.power_by_squaring(x,m) :
                                     Base.power_by_squaring(inv(x),-m)
Base.cmp(x::Mod{p}, y::Mod{p}) where p=cmp(x.val,y.val)
Base.isless(x::Mod{p}, y::Mod{p}) where p=cmp(x,y)==-1
Base.abs(x::Mod)=x      # needed for inv(Matrix) to work
Base.conj(x::Mod)=x     # needed for inv(Matrix) to work
(::Type{T})(x::Mod) where T<:Integer=T(x.val)

function Base.show(io::IO, m::Mod{p}) where p 
  if get(io,:limit,false) && get(io,:typeinfo,Any)==typeof(m) print(io,m.val)
  else print(io,typeof(m),"(",m.val,")")
  end
end

function order(x::Mod{n}) where n
  if n==1 return 1 end
  d=gcd(x.val,n)
  if d!=1 error("$x should be invertible") end
  o=1
  res=x
  while true
   if isone(res) return o end
   o+=1
   res*=x
  end
end

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
    error("missing conway($p,$n)") 
  end
end

struct FF
  p::Int16
  n::Int16
  q::Int16
  conway::Vector
  zech::Vector{Int16}
  dict::Vector{Int16}
end

struct FFE{p}<:Number
  i::Int16
  Fi::Int16
end

printc(x...)=println(IOContext(stdout,:limit=>true),x...)

const FFi=Dict{Int,Int}()
const FFvec=FF[]

function iFF(q)
  get!(FFi,q) do
    l=collect(factor(q))
    if length(l)>1 error(q," should be a prime power") end
    p,n=l[1]
    pol=-Mod{p}.(conway_polynomial(p,n))
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
  # printc("conway=",pol)
    zz=map(i->fill(Mod{p}(0),n),1:q)
    zz[1][1]=Mod{p}(1)
    for i in 2:q-1
      for j in 2:n zz[i][j]=zz[i-1][j-1] end
      zz[i].+= zz[i-1][n].*pol
    end
    z=collect(enumerate(reverse.(zz)))
    sort!(z,by=x->x[2])
#   for i in 1:q printc(i,"->",(z[i][1]-1,z[i][2])) end
    zz=Vector{Int16}(undef,q)
    for i in 1:q
      if z[i][2][end]==Mod{p}(p-1) zz[z[i][1]]=z[i-p+1][1]-1
      else zz[z[i][1]]=z[i+1][1]-1
      end
    end
#   for i in 0:q-1 printc(i,"=>",zz[i+1]) end
    push!(FFvec,FF(p,n,q,pol,zz,dic))
    return length(FFvec)
  end
end

FF(q)=FFvec[iFF(q)]

Base.show(io::IO,F::FF)=print(io,"FF(",F.p,"^",F.n,")")

@inbounds @inline field(x::FFE)=FFvec[x.Fi]
@inline q(x::FFE)=field(x).q
@inline n(x::FFE)=field(x).n
@inbounds @inline zech(x::FFE,i)=field(x).zech[1+i]
@inbounds @inline dict(x::FFE,i)=field(x).dict[1+i]
@inline clone(x::FFE{p},i) where p=FFE{p}(i,x.Fi)

Base.iszero(x::FFE{p}) where p=x.i==-1

Base.:^(a::FFE{p},n::Integer) where p=iszero(a) ? a : 
                 lower(clone(a,mod(a.i*n,q(a)-1)))

function Z(q)
  if q==2 return FFE{2}(0,iFF(2)) end
  Fi=iFF(q)
  FFE{FFvec[Fi].p}(1,Fi)
end

Base.promote_rule(::Type{FFE{p}},::Type{<:Integer}) where p=FFE{p}
Base.promote_rule(::Type{FFE{p}},::Type{<:Rational{<:Integer}}) where p=FFE{p}

Base.copy(a::FFE{p}) where p=clone(a,a.i)
Base.one(::Type{FFE{p}}) where p=FFE{p}(0,iFF(p))
Base.one(x::FFE{p}) where p=q(x)==p ? clone(x,0) : one(FFE{p})
Base.zero(::Type{FFE{p}}) where p=FFE{p}(-1,iFF(p))
Base.zero(x::FFE{p}) where p=q(x)==p ? clone(x,-1) : zero(FFE{p})
Base.abs(a::FFE)=a
Base.conj(a::FFE)=a

function Base.cmp(x::FFE{p}, y::FFE{p}) where {p}
  l=cmp(q(x),q(y))
  if !iszero(l) return l end
  cmp(x.i,y.i)
end

Base.isless(x::FFE, y::FFE)=cmp(x,y)==-1
Base.:(==)(x::FFE, y::FFE)=cmp(x,y)==0

function Base.show(io::IO, x::FFE{p})where p
  sup=Dict(zip("0123456789","⁰¹²³⁴⁵⁶⁷⁸⁹"))
  sub=Dict(zip("0123456789,()","₀₁₂₃₄₅₆₇₈₉‚₍₎"))
  repl=get(io,:limit,false)
  function printz()
    if repl print(io,"Z",map(x->sub[x],repr(q(x))))
    else print(io,"Z(",q(x),")")
    end
  end
  if iszero(x) 
     if repl print(io,"0",map(x->sub[x],repr(p)))
     else print(io,"0*");printz()
     end
  elseif isone(x) 
     if repl print(io,"1",map(x->sub[x],repr(p)))
     else printz();print(io,"^",x.i)
     end
  else  printz()
    if x.i!=1
      if repl print(io,map(x->sup[x],repr(x.i))) else print(io,"^",x.i) end
    end
  end
end

function (::Type{FFE{p}})(i::Integer)where p
  Fi=iFF(p)
  FFE{p}(FFvec[Fi].dict[1+mod(i,p)],Fi)
end

function (::Type{FFE{p}})(i::Mod{p})where p
  Fi=iFF(p)
  FFE{p}(FFvec[Fi].dict[1+i.val],Fi)
end

function (::Type{Mod{r}})(x::FFE{p})where {r,p}
  if r!=p || q(x)!=p throw(InexactError(:convert,Mod{r},x)) end
  iszero(x) ? Mod{p}(0) : Mod{p}(field(x).conway[1].val)^x.i
end

function (::Type{FFE{p}})(i::Rational{<:Integer})where p
  Fi=iFF(p)
  i=mod(invmod(denominator(i),p)*numerator(i),p)
  FFE{p}(FFvec[Fi].dict[1+i],Fi)
end

function promote_field(r::Integer,x::FFE{p})where {p}
  if r==q(x) return x end
  if r<q(x) error("cannot promote to smaller field") end
  Fi=iFF(r)
@inbounds  F=FFvec[Fi]
  if F.p!=p error("different characteristic!") end
  if mod(F.n,n(x))!=0 error("not an extension field") end
  FFE{p}(x.i*div(r-1,q(x)-1),Fi)
end
  
function promote_field(x::FFE{p},y::FFE{p})where {p}
  if x.Fi==y.Fi return (x,y) end
  nq=p^lcm(n(x),n(y))
  (promote_field(nq,x),promote_field(nq,y))
end

function Base.:+(x::FFE{p},y::FFE{p}) where {p}
  if iszero(x) return y
  elseif iszero(y) return x
  end
  x,y=promote_field(x,y)
  if x.i>y.i x,y=y,x end
  @inbounds res=zech(x,y.i-x.i)
  if res==q(x)-1 return zero(x) end
  res+=x.i
  if res>=q(x)-1 res-=Int16(q(x)-1) end
  lower(clone(x,res))
end

function Base.:-(x::FFE{p})where p
  if iseven(p) || iszero(x) return x end
  res=x.i+div(q(x)-1,2)
  if res>=q(x)-1 res-=Int16(q(x)-1) end
  clone(x,res)
end

Base.:-(x::FFE,y::FFE)=x+(-y)

function Base.:*(x::FFE{p},y::FFE{p}) where {p}
  if iszero(x) return x
  elseif iszero(y) return y
  end
  (x,y)=promote_field(x,y)
  res=x.i+y.i
  if res>=q(x)-1 res-=Int16(q(x)-1) end
  lower(clone(x,res))
end

function Base.inv(x::FFE{p}) where p
  if iszero(x) error("inv(0)") end
  if x.i==0 return x end
  clone(x,q(x)-1-x.i)
end

Base.:/(x::FFE,y::FFE)=x*inv(y)
Base.://(x::FFE,y::FFE)=x*inv(y)

function lower(x::FFE{p}) where p
  if q(x)==p return x end
  l=dict(x,x.i)
  if l==n(x) return x end
  FFE{p}(div(x.i,div(q(x)-1,p^l-1)),iFF(p^l))
end

function order(x::FFE{p}) where p
  if iszero(x) error(x," must be invertible") end
  div(q(x)-1,gcd(x.i,q(x)-1))
end

function Base.log(x::FFE)
  if iszero(x) error(x," must be invertible") end
  x.i
end

end
