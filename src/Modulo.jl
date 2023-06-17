"""
This module introduces modular arithmetic.

The  integer `x` mod. `n` is constructed  by the function `Mod(x,n)`. If `n
isa  Int` the result is of type `Mod{UInt64}`. If `n isa BigInt` the result
is  of  type  `Mod{BigInt}`.  Since  `n`  is  not  encoded in the type, the
elements  `0` and `1` mod.  `n` cannot be constructed  from the type, which
causes  some problems for some Julia functionality (for instance `inv` on a
matrix does not work). For prime moduli `p`, the type `FFE{p}` in `FiniteFields`
does not have such limitations.

Example:
```julia-repl
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
"""
module Modulo
export Mod, order

struct Mod{T}<:Number
  val::T
  n::T
  global Mod_(a::T,n::T) where T=new{T}(a,n)
end

Mod(a::Integer,n)=Mod_(typeof(unsigned(n))(mod(a,n)),unsigned(n))

Mod(a::Integer,n::BigInt)=Mod_(mod(a,n),n)

# can reduce further to modulus dividing a.n
function Mod(a::Mod,n::Integer)
# in Fact.jl we use this to go from mod. n to mod. n^2
# if !iszero(a.n%n) error(n," does not divide ",a.n) end
  Mod(a.val,n)
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
Base.isone(n::Mod)=isone(n.val) || isone(n.n)
Base.:(==)(x::Mod,y::Mod)=x.n==y.n && x.val==y.val
Base.:(==)(x::Mod,y::Integer)=x==Mod(y,x.n)
function Base.:+(x::Mod, y::Mod)
  if x.n!=y.n  error("moduli") end
  res=x.val+y.val
  if res>x.n res-=x.n end
  Mod(res,x.n)
end

Base.:*(x::Mod, y::Mod)=x.n!=y.n ? error("moduli") : Mod(widemul(x.val,y.val),x.n)
Base.:-(x::Mod)=iszero(x) ? x : Mod(x.n-x.val,x.n)
Base.:-(x::Mod, y::Mod)=x+(-y)
Base.:/(x::Mod,y::Mod)=x*inv(y)
Base.:/(x::Mod,y::Number)=x*inv(Mod(y,x.n))
Base.inv(x::Mod)=Mod(invmod(x.val,x.n),x.n)
Base.:^(x::Mod,m::Integer)=m>=0 ? Base.power_by_squaring(x,m) :
                                  Base.power_by_squaring(inv(x),-m)
Base.cmp(x::Mod,y::Mod)=cmp(x.val,y.val)
Base.gcd(x::Mod,y::Mod)=Mod(gcd(x.val,y.val),x.n)
Base.isless(x::Mod,y::Mod)=cmp(x,y)==-1
Base.abs(x::Mod)=x   # needed for inv(Matrix) to work
Base.conj(x::Mod)=x  # needed for inv(Matrix) to work
(::Type{T})(x::Mod) where T<:Signed=x.val<=x.n-x.val ? T(x.val) : -T(x.n-x.val)

(::Type{T})(x::Mod) where T<:Unsigned=T(x.val)

Integer(x::Mod)=Signed(x)

function Base.show(io::IO, ::MIME"text/plain", x::Mod)
  if !haskey(io,:typeinfo) print(io,typeof(x),": ") end
  show(io,x)
end

const sub=Dict(zip("-0123456789,+()=aehijklmnoprstuvxβγρφχ.",
                   "₋₀₁₂₃₄₅₆₇₈₉‚₊₍₎₌ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵦᵧᵨᵩᵪ̣."))

function Base.show(io::IO, m::Mod)
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

end
