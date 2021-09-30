"""
An   implementation  of  univariate  Laurent  polynomials,  and  univariate
rational fractions.

A Pol contains two fields: a vector of coefficients, and the valuation.

# Examples
```julia-repl
julia> Pol(:q) # define string used for printing and returns Pol([1],1)
Pol{Int64}: q

julia> @Pol q # same as q=Pol(:q)
Pol{Int64}: q

julia> Pol([1,2]) # valuation is 0 if not specified
Pol{Int64}: 2q+1

julia> 2q+1       # same result
Pol{Int64}: 2q+1

julia> Pol()   # omitting all arguments gives Pol([1],1)
Pol{Int64}: q

julia> p=Pol([1,2],-1) # here the valuation is specified to be -1
Pol{Int64}: 2+q⁻¹

julia> valuation(p),degree(p)
(-1, 0)

julia> derivative(p)
Pol{Int64}: -q⁻²

julia> p=(q+1)^2
Pol{Int64}: q²+2q+1

julia> valuation(p),degree(p)
(0, 2)

julia> p(1//2) # value of p at 1//2
9//4

julia> p[0], p[1], p[-1] # indexing gives the coefficients
(1, 2, 0)

julia> p[-1:1]
3-element Vector{Int64}:
 0
 1
 2

julia> divrem(q^3+1,2q+1) # changes coefficients to field elements
(0.5q²-0.25q+0.125, 0.875)

julia> divrem(q^3+1,2q+1//1) # case of field elements
((1//2)q²+(-1//4)q+1//8, 7//8)

julia> Pols.pseudodiv(q^3+1,2q+1) # pseudo-division keeps the ring
(4q²-2q+1, 7)

julia> m=[q+1 q+2;q-2 q-3]
2×2 Matrix{Pol{Int64}}:
 q+1  q+2
 q-2  q-3

julia> inv(RatFrac.(m))
2×2 Matrix{RatFrac{Int64}}:
 (-q+3)/(2q-1)  (-q-2)/(-2q+1)
 (q-2)/(2q-1)   (q+1)/(-2q+1)
```

see also the individual documentation of divrem, gcd.
"""
module Pols
using ..Util: Util, format_coefficient, printTeX, exactdiv
export degree, valuation, Pol, derivative, shift, positive_part, negative_part,
       bar, derivative, srgcd, RatFrac, @Pol

const varname=Ref(:x)

struct Pol{T}
  c::Vector{T}
  v::Int
  # beware that c is not necessarily copied
  function Pol(c::AbstractVector{T},v::Integer=0;check=true)where T
    if check # normalize c so there are no leading or trailing zeroes
      b=findfirst(!iszero,c)
      if b===nothing return new{T}(empty(c),0) end
      e=findlast(!iszero,c)
      if b!=1 || e!=length(c) return new{T}(view(c,b:e),v+b-1) end
    end
    new{T}(c,v)
  end
end

Pol(a::Number)=convert(Pol,a)
Pol()=Pol([1],1;check=false)

function Pol(t::Symbol)
  varname[]=t
  Pol()
end

"@Pol q is equivalent to q=Pol(:q)"
macro Pol(t)
  if !(t isa Symbol) error("usage: @Pol <variable name>") end
  Base.eval(Main,:($t=Pol($(Core.QuoteNode(t)))))
end

Base.broadcastable(p::Pol)=Ref(p)

degree(a::Number)=0 # convenient
degree(p::Pol)=length(p.c)-1+p.v
valuation(p::Pol)=p.v
Base.lastindex(p::Pol)=degree(p)
Base.getindex(p::Pol{T},i::Integer) where T=i in p.v:lastindex(p) ?
    p.c[i-p.v+1] : zero(T)

Base.getindex(p::Pol,i::AbstractVector{<:Integer})=getindex.(Ref(p),i)

Base.copy(p::Pol)=Pol(copy(p.c),p.v;check=false)
Base.convert(::Type{Pol{T}},a::Number) where T=iszero(a) ? zero(Pol{T}) :
        Pol([T(a)];check=false)
Base.convert(::Type{Pol},a::Number)=convert(Pol{typeof(a)},a)
(::Type{Pol{T}})(a) where T=convert(Pol{T},a)
Base.convert(::Type{Pol{T}},p::Pol{T1}) where {T,T1}= T==T1 ? p :
     Pol(convert(Vector{T},p.c),p.v;check=false)

function Base.promote_rule(a::Type{Pol{T1}},b::Type{Pol{T2}})where {T1,T2}
  Pol{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Pol{T1}},b::Type{T2})where {T1,T2<:Number}
  Pol{promote_type(T1,T2)}
end

Base.isinteger(p::Pol)=iszero(p) || (iszero(p.v) && isone(length(p.c)) &&
                                     isinteger(p.c[1]))
function Base.convert(::Type{T},p::Pol) where T<:Number
  if iszero(p) return T(0) end
  if !iszero(degree(p)) || !iszero(valuation(p))
    throw(InexactError(:convert,T,p))
  end
  convert(T,p.c[1])
end

(::Type{T})(p::Pol) where T<:Number=convert(T,p)

Base.cmp(a::Pol,b::Pol)=cmp((a.c,a.v),(b.c,b.v))
Base.isless(a::Pol,b::Pol)=cmp(a,b)==-1
Base.hash(a::Pol, h::UInt)=hash(a.v,hash(a.c,h))

(p::Pol{T})(x) where T=iszero(p) ? zero(T) : evalpoly(x,p.c)*x^p.v

# efficient p↦ qˢ p
shift(p::Pol,s)=Pol(p.c,p.v+s;check=false)

Base.denominator(p::Pol)=iszero(p) ? 1 : lcm(denominator.(p.c))

function positive_part(p::Pol)
  if p.v>=0 return p end
  Pol(view(p.c,1-p.v:length(p.c)),0)
end

function negative_part(p::Pol)
  if degree(p)<=0 return p end
  Pol(view(p.c,1:1-p.v),p.v)
end

# q↦ q⁻¹ on p
bar(p::Pol)=Pol(reverse(p.c),-degree(p);check=false)

Base.:(==)(a::Pol, b::Pol)= a.c==b.c && a.v==b.v
Base.:(==)(a::Pol,b)= a==Pol(b)
Base.:(==)(b,a::Pol)= a==Pol(b)

Base.one(a::Pol{T}) where T=Pol([iszero(a) ? one(T) : one(a.c[1])];check=false)
Base.one(::Type{Pol{T}}) where T=Pol([one(T)];check=false)
Base.one(::Type{Pol})=one(Pol{Int})
Base.isone(a::Pol)=iszero(a.v) && length(a.c)==1 && isone(a.c[1])
Base.zero(::Type{Pol{T}}) where T=Pol(T[];check=false)
Base.zero(::Type{Pol})=zero(Pol{Int})
Base.zero(a::Pol)=Pol(empty(a.c);check=false)
Base.iszero(a::Pol)=isempty(a.c)
# next 3 stuff to make inv using LU work (abs is stupid)
Base.conj(p::Pol)=Pol(conj.(p.c),p.v;check=false)
Base.abs(p::Pol)=p
Base.adjoint(a::Pol)=conj(a)

function Base.show(io::IO, ::MIME"text/html", a::Pol)
  print(io, "\$")
  show(IOContext(io,:TeX=>true),a)
  print(io, "\$")
end

function Base.show(io::IO, ::MIME"text/plain", a::Pol)
  if !haskey(io,:typeinfo) print(io,typeof(a),": ") end
  show(io,a)
end

function Base.show(io::IO,p::Pol)
  if !get(io,:limit,false) && !get(io,:TeX,false)
    if length(p.c)==1 && isone(p.c[1]) && p.v==1 print(io,"Pol()")
    else print(io,"Pol(",p.c,",",p.v,")")
    end
  elseif iszero(p) print(io,"0")
  else
    mon=string(get(io,:varname,varname[]))
    for deg in degree(p):-1:valuation(p)
      c=p[deg]
      if iszero(c) continue end
      c=repr(c; context=IOContext(io,:typeinfo=>typeof(c)))
      if !iszero(deg)
        c=format_coefficient(c)*mon
        if !isone(deg) c*="^{$deg}" end
      end
      if c[1]!='-' && deg!=degree(p) c="+"*c end
      printTeX(io,c)
    end
  end
end

function Base.:*(a::Pol,b::Pol)
  if iszero(a) || iszero(b) return zero(a) end
  z=zero(a.c[1]+b.c[1])
  res=fill(z,length(a.c)+length(b.c)-1)
  for i in eachindex(a.c), j in eachindex(b.c)
@inbounds res[i+j-1]+=a.c[i]*b.c[j]
  end
  Pol(res,a.v+b.v;check=false)
end

Base.:*(a::Pol, b::Number)=Pol(a.c.*b,a.v)
Base.:*(a::Pol{T}, b::T) where T=Pol(a.c.*b,a.v)
Base.:*(b::Number, a::Pol)=a*b
Base.:*(b::T, a::Pol{T}) where T=a*b

Base.:^(a::Pol, n::Real)=length(a.c)==1 ? Pol([a.c[1]^n],n*a.v) :
         n>=0 ? Base.power_by_squaring(a,Int(n)) :
                Base.power_by_squaring(inv(a),Int(-n))

function Base.:+(a::Pol,b::Pol)
  d=b.v-a.v
  if d<0 return b+a end
  if iszero(a) return b elseif iszero(b) return a end
  z=zero(a.c[1]+b.c[1])
  res=fill(z,max(length(a.c),d+length(b.c)))
@inbounds view(res,eachindex(a.c)).=a.c
@inbounds view(res,d.+eachindex(b.c)).+=b.c
  Pol(res,a.v)
end

Base.:+(a::Pol, b::Number)=a+Pol(b)
Base.:+(b::Number, a::Pol)=Pol(b)+a
Base.:-(a::Pol)=Pol(-a.c,a.v;check=false)
Base.:-(a::Pol, b::Pol)=a+(-b)
Base.:-(a::Pol, b::Number)=a-Pol(b)
Base.:-(b::Number, a::Pol)=Pol(b)-a
Base.div(a::Pol,b::Int)=Pol(div.(a.c,b),a.v;check=false)

derivative(a::Pol)=Pol([(i+a.v-1)*v for (i,v) in enumerate(a.c)],a.v-1)

bestinv(x)=isone(x) ? x : isone(-x) ? x : inv(x)

"""
`divrem(a::Pol, b::Pol)`

computes  `(q,r)` such  that `a=q*b+r`  and `degree(r)<degree(b)`. When the
leading  coefficient  of  b  is  ±1  does  not  use  inverse.
For true polynomials (errors if the valuation of `a` or of `b` is negative).
"""
function Base.divrem(a::Pol, b::Pol)
  if iszero(b) throw(DivideError) end
  if degree(b)>degree(a) return (zero(a),a) end
  d=inv(b.c[end])
  z=zero(a.c[1]+b.c[1]+d)
  r=fill(z,1+degree(a))
  view(r,a.v+1:length(r)).=a.c
  q=fill(z,length(r)-degree(b))
  for i in length(r):-1:degree(b)+1
    if iszero(r[i]) c=zero(d)
    else c=r[i]*d
         view(r,i-length(b.c)+1:i) .-= c .* b.c
    end
    q[i-degree(b)]=c
  end
  Pol(q),Pol(r)
end

function Util.exactdiv(a::Pol,b::Pol)
  if isone(b) || iszero(a) return a end
  if iszero(b) throw(DivideError) end
  d=a.v-b.v
  if !iszero(a.v) a=shift(a,-a.v) end
  if !iszero(b.v) b=shift(b,-b.v) end
  if degree(b)>degree(a) error(b," does not divide exactly ",a) end
  z=zero(a.c[1]+b.c[1])
  r=fill(z,1+degree(a))
  view(r,a.v+1:length(r)).=a.c
  q=fill(z,length(r)-degree(b))
  for i in length(r):-1:degree(b)+1
    c=exactdiv(r[i],b.c[end])
    view(r,i-length(b.c)+1:i) .-= c .* b.c
    q[i-length(b.c)+1]=c
  end
  if !iszero(r) error(b," does not divide exactly ",a) end
  res=Pol(q)
  !iszero(d) ? shift(res,d) : res
end

"""
`pseudodiv(a::Pol, b::Pol)`

pseudo-division  of `a` by `b`.  If `d` is the  leading coefficient of `b`,
computes   `(q,r)`   such   that   `d^(degree(a)+1-degree(b))a=q*b+r`   and
`degree(r)<degree(b)`. Does not do division so works over any ring.
For true polynomials (errors if the valuation of `a` or of `b` is negative).

See Knuth AOCP2 4.6.1 Algorithm R
"""
function pseudodiv(a::Pol, b::Pol)
  if iszero(b) throw(DivideError) end
  d=b.c[end]
  if degree(a)<degree(b) return (Pol(0),d^(degree(a)+1-degree(b))*a) end
  z=zero(promote_type(eltype(a.c),eltype(b.c)))
  r=fill(z,1+degree(a))
  view(r,a.v+1:length(r)).=a.c
  q=fill(z,length(r)-degree(b))
  for i in length(r):-1:degree(b)+1
    c=r[i]
    r.*=d
    q.*=d
    if !iszero(c)
      for j in eachindex(b.c) r[j+i-length(b.c)]-=c*b.c[j] end
    end
    q[i-degree(b)]=c
  end
  Pol(q),Pol(r)
end

"""
`srgcd(a::Pol,b::Pol)`

sub-resultant gcd: gcd of polynomials over a unique factorization domain

See Knuth AOCP2 4.6.1 Algorithm C
"""
function srgcd(a::Pol,b::Pol)
  if degree(b)>degree(a) a,b=b,a end
  if iszero(b) return a end
  ca=gcd(a.c);if !isone(ca) a=Pol(exactdiv.(a.c,ca),a.v;check=false) end
  cb=gcd(b.c);if !isone(cb) b=Pol(exactdiv.(b.c,cb),b.v;check=false) end
  d=gcd(ca,cb)
  g=1
  h=1
  while true
    δ=degree(a)-degree(b)
    q,r=pseudodiv(a,b)
    if iszero(r)
      cb=gcd(b.c);if !isone(cb) b=Pol(exactdiv.(b.c,cb),b.v;check=false) end
      return isone(d) ? b : Pol(b.c .*d,b.v;check=false)
    elseif degree(r)==0
      return Pol([d];check=false)
    end
    a=b
    gh=g*h^δ
    b=isone(gh) ? r : Pol(exactdiv.(r.c,gh),r.v;check=false)
    g=a[end]
    if δ>0 h=exactdiv(g^δ,h^(δ-1)) end
  end
end

Base.gcd(p::Pol{<:Integer},q::Pol{<:Integer})=srgcd(p,q)
Base.gcd(v::Array{<:Pol})=reduce(gcd,v)
Base.lcm(p::Pol,q::Pol)=exactdiv(p*q,gcd(p,q))
Base.lcm(m::Array{<:Pol})=reduce(lcm,m)

Base.div(a::Pol, b::Pol)=divrem(a,b)[1]
Base.:%(a::Pol, b::Pol)=divrem(a,b)[2]
Base.:%(a::Pol{<:Integer}, b::Pol{<:Integer})=isone(b[end]) ?
                                      pseudodiv(a,b)[2] : divrem(a,b)[2]

"""
`gcd(p::Pol,  q::Pol)` computes the  `gcd` of the  polynomials. It uses the
subresultant algorithms for the `gcd` of integer polynomials.

```julia-repl
julia> gcd(2q+2,q^2-1)
Pol{Int64}: q+1

julia> gcd(q+1//1,q^2-1//1)
Pol{Rational{Int64}}: (1//1)q+1//1
```
"""
function Base.gcd(p::Pol,q::Pol)
  if degree(q)>degree(p)
    p,q=q,p
  end
  p,q=promote(p,q)
  while !iszero(q)
    q=q/q.c[end]
    (q,p)=(divrem(p,q)[2],q)
  end
  p*inv(p.c[end])
end

"""
  `gcdx(a,b)` works for polynomials over a field.
  returns `d,u,v` such that `d=ua+vb` and `d=gcd(a,b)`.
```julia-repl
julia> gcdx(q^3-1//1,q^2-1//1)
((1//1)q-1//1, 1//1, (-1//1)q)
```
"""
function Base.gcdx(a::Pol, b::Pol)
  a,b=promote(a, b)
  # a0, b0=a, b
  s0, s1=one(a), zero(a)
  t0, t1=s1, s0
  # The loop invariant is: s0*a0 + t0*b0 == a
  x,y=a,b
  while y != 0
    q,q1=divrem(x, y)
    x, y=y, q1
    s0, s1=s1, s0 - q*s1
    t0, t1=t1, t0 - q*t1
  end
  (x, s0, t0)./x[end]
end

#powermod for Pols over a field
function Base.powermod(p::Pol, x::Integer, q::Pol)
  x==0 && return one(q)
  b=p%q
  t=prevpow(2, x)
  r=one(q)
  while true
    if x>=t
     r=(r*b)%q
      x-=t
    end
    t >>>= 1
    t<=0 && break
    r=(r*r)%q
  end
  r
end

# random polynomial of degree d
Base.rand(::Type{Pol{T}},d::Integer) where T=Pol(rand(T,d+1))

# Interpolation: find Pol of smallest degree taking values y at points x
function Pol(x::AbstractVector,y::AbstractVector)
  t=collect(y).*1//1 # make sure coeffs are in a field
  a=map(eachindex(x))do i
    for k in i-1:-1:1
      if x[i]==x[k] error("interpolating points must be distinct") end
      t[k]=(t[k+1]-t[k])/(x[i]-x[k])
    end
    t[1]
  end
  p=Pol([a[end]])
  for i in length(x)-1:-1:1
    p=p*(Pol()-x[i])+a[i]
  end
  p
end

#---------------------- RatFrac-------------------------------------
struct RatFrac{T}
  num::Pol{T}
  den::Pol{T}
  function RatFrac(a::Pol{T1},b::Pol{T2};check=true,prime=false)where {T1,T2}
    T=promote_type(T1,T2)
    if check
      if iszero(a) return new{T}(a,one(a))
      elseif iszero(b) error("zero denominator")
      end
      v=a.v-b.v
      a=shift(a,max(v,0)-a.v)
      b=shift(b,-min(v,0)-b.v)
      if !prime
        d=gcd(a,b)
        a=exactdiv(a,d)
        b=exactdiv(b,d)
      end
      if isone(-b) a,b=(-a,-b) end
    end
    new{T}(a,b)
  end
end

function Base.convert(::Type{RatFrac{T}},p::RatFrac{T1}) where {T,T1}
  RatFrac(convert(Pol{T},p.num),convert(Pol{T},p.den);check=false)
end

function Base.convert(::Type{Pol{T}},p::RatFrac{T1}) where {T,T1}
  if length(p.den.c)==1
    return convert(Pol{T},Pol(p.num.c .//p.den.c[1],p.num.v-p.den.v))
  end
  error("cannot convert ",p," to Pol{",T,"}")
end

function Base.convert(::Type{RatFrac{T}},p::Pol{T1}) where {T,T1}
  RatFrac(convert(Pol{T},p),Pol(T(1));prime=true)
end

function Base.convert(::Type{RatFrac{T}},p::Number) where {T}
  RatFrac(convert(Pol{T},p),Pol(T(1));check=false)
end

function Base.promote_rule(a::Type{Pol{T1}},b::Type{RatFrac{T2}})where {T1,T2}
  RatFrac{promote_type(T1,T2)}
end

(::Type{RatFrac{T}})(a) where T=convert(RatFrac{T},a)

Base.broadcastable(p::RatFrac)=Ref(p)
RatFrac(a::Number)=RatFrac(Pol(a),Pol(1);check=false)
RatFrac(a::Pol)=RatFrac(a,Pol(1);prime=true)
Base.copy(a::RatFrac)=RatFrac(a.num,a.den;check=false)
Base.one(a::RatFrac)=RatFrac(one(a.num),one(a.den);check=false)
Base.one(::Type{RatFrac{T}}) where T =RatFrac(one(Pol{T}),one(Pol{T});check=false)
Base.one(::Type{RatFrac}) where T =one(RatFrac{Int})
Base.zero(::Type{RatFrac}) where T =zero(RatFrac{Int})
Base.zero(::Type{RatFrac{T}}) where T =RatFrac(zero(Pol{T}),one(Pol{T});check=false)
Base.zero(a::RatFrac)=RatFrac(zero(a.num),one(a.num);check=false)
Base.iszero(a::RatFrac)=iszero(a.num)
# next 3 stuff to make inv using LU work (abs is stupid)
Base.abs(p::RatFrac)=p
Base.conj(p::RatFrac)=RatFrac(conj(p.num),conj(p.den);check=false)
Base.adjoint(a::RatFrac)=conj(a)
Base.cmp(a::RatFrac,b::RatFrac)=cmp([a.num,a.den],[b.num,b.den])
Base.isless(a::RatFrac,b::RatFrac)=cmp(a,b)==-1

function Base.show(io::IO, ::MIME"text/plain", a::RatFrac)
  if !haskey(io,:typeinfo) print(io,typeof(a),": ") end
  show(io,a)
end

function Base.show(io::IO,a::RatFrac)
  n=sprint(show,a.num; context=io)
  if  get(io, :limit,true) && isone(a.den)
    print(io,n)
  else
    print(io,Util.bracket_if_needed(n))
    n=sprint(show,a.den; context=io)
    print(io,"/",Util.bracket_if_needed(n))
  end
end

Base.inv(a::RatFrac)=RatFrac(a.den,a.num;check=false)

function Base.://(a::Pol,b::Pol)
  if b.c==[1] return shift(a,-b.v)
  elseif b.c==[-1] return shift(-a,-b.v)
  elseif length(b.c)==1 return Pol(a.c.//b.c[1],a.v-b.v)
  end
  RatFrac(a,b)
end

function Base.inv(p::Pol)
  if length(p.c)==1 return Pol([bestinv(p.c[1])],-p.v) end
  RatFrac(Pol(1),p;prime=true)
end

Base.:/(p::Pol,q::Pol)=p*inv(q)

Base.://(a::RatFrac,b::RatFrac)=a*inv(b)
Base.://(a::RatFrac,b::Union{Number,Pol})=RatFrac(a.num,a.den*b)
Base.:/(a::RatFrac,b::Union{Number,Pol})=a//b
Base.://(a::Union{Number,Pol},b::RatFrac)=a*inv(b)
Base.:/(a::Union{Number,Pol},b::RatFrac)=a//b
Base.:/(a::RatFrac,b::RatFrac)=a*inv(b)
Base.://(p::Number,q::Pol)=RatFrac(Pol(p),q;prime=true)
Base.:/(p::Number,q::Pol)=p*inv(q)
Base.:/(p::Pol,q::Number)=Pol(p.c ./q,p.v)
Base.://(p::Pol,q::Number)=iszero(p) ? p : Pol(p.c//q,p.v)

Base.:*(a::RatFrac,b::RatFrac)=RatFrac(a.num*b.num,a.den*b.den)

Base.:*(a::RatFrac,b::Pol)=RatFrac(a.num*b,a.den)
Base.:*(b::Pol,a::RatFrac)=RatFrac(a.num*b,a.den)
Base.:*(a::RatFrac,b::T) where T =RatFrac(a.num*b,a.den;check=false)
Base.:*(b::T,a::RatFrac) where T =a*b

Base.:^(a::RatFrac, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                              Base.power_by_squaring(inv(a),-n)
Base.:+(a::RatFrac,b::RatFrac)=RatFrac(a.num*b.den+a.den*b.num,a.den*b.den)
Base.:+(a::RatFrac,b::Number)=a+RatFrac(b)
Base.:+(b::Number,a::RatFrac)=a+RatFrac(b)
Base.:-(a::RatFrac)=RatFrac(-a.num,a.den;check=false)
Base.:-(a::RatFrac,b::RatFrac)=a+(-b)
Base.:-(a::RatFrac,b)=a-RatFrac(b)
Base.:-(b,a::RatFrac)=RatFrac(b)-a

(p::RatFrac)(x)=p.num(x)//p.den(x)

end
