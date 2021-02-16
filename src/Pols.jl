"""
 An implementation of univariate Laurent polynomials.
 A Pol contains two fields: a vector of coefficients, and the valuation.

# Examples
```julia-repl
julia> Pol(:q) # define string used for printing and set variable q
Pol{Int64}: q

julia> Pol([1,2]) # valuation is 0 if not specified
Pol{Int64}: 2q+1

julia> 2q+1       # same result
Pol{Int64}: 2q+1

julia> p=Pol([1,2],-1) # here the valuation is specified to be -1
Pol{Int64}: 2+q⁻¹

julia> Pol()   # omitting all arguments gives Pol([1],1)
Pol{Int64}: q

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

julia> divrem(q^3+1,q+2)  # keeps the ring if leading coeff of divisor is ±1
(q²-2q+4, -7)
```

see also the individual documentation of divrem, gcd.
"""
module Pols
using ..Util: format_coefficient, printTeX
export degree, valuation
export Pol, shift, positive_part, negative_part, bar, derivative

const varname=Ref(:x)

struct Pol{T}
  c::Vector{T}
  v::Int
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
  Base.eval(Main,:($t=Pol()))
end

Base.broadcastable(p::Pol)=Ref(p)

degree(p::Pol)=length(p.c)-1+p.v
valuation(p::Pol)=p.v
Base.lastindex(p::Pol)=degree(p)
Base.getindex(p::Pol{T},i::Integer) where T=i in p.v:lastindex(p) ? 
    p.c[i-p.v+1] : zero(T)

Base.getindex(p::Pol,i::AbstractVector{<:Integer})=getindex.(Ref(p),i)

Base.copy(p::Pol)=Pol(p.c,p.v;check=false)
Base.convert(::Type{Pol{T}},a::Number) where T=iszero(a) ? zero(Pol{T}) :
        Pol([T(a)];check=false)
Base.convert(::Type{Pol},a::Number)=convert(Pol{typeof(a)},a)
(::Type{Pol{T}})(a) where T=convert(Pol{T},a)
Base.convert(::Type{Pol{T}},p::Pol{T1}) where {T,T1}= T==T1 ? p :
        Pol(convert.(T,p.c),p.v;check=false)

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

Base.cmp(a::Pol,b::Pol)=cmp([a.c,a.v],[b.c,b.v])
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

Base.one(a::Pol{T}) where T=Pol([one(T)];check=false)
Base.one(::Type{Pol{T}}) where T=Pol([one(T)];check=false)
Base.one(::Type{Pol})=one(Pol{Int})
Base.zero(::Type{Pol{T}}) where T=Pol(T[];check=false)
Base.zero(::Type{Pol})=zero(Pol{Int})
Base.zero(a::Pol)=Pol(empty(a.c);check=false)
Base.iszero(a::Pol)=isempty(a.c)
# next 3 stuff to make inv using LU work (abs is stupid)
Base.abs(p::Pol)=p
Base.conj(p::Pol)=Pol(conj.(p.c),p.v;check=false)
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
    if p.c==[1] && p.v==1 print(io,"Pol()")
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

function Base.:*(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  if iszero(a) || iszero(b) return zero(a) end
  res=fill(zero(promote_type(T1,T2)),length(a.c)+length(b.c)-1)
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

function Base.:+(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  d=b.v-a.v
  if d<0 return b+a end
  res=fill(zero(promote_type(T1,T2)),max(length(a.c),d+length(b.c)))
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

computes `(p,q)` such that `a=p*b+q`
When the leading coefficient of b is ±1 does not use inverse
"""
function Base.divrem(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  if iszero(b) throw(DivideError) end
  d=bestinv(b.c[end])
  T=promote_type(T1,T2,typeof(d))
  v=convert(Vector{T},copy(a.c))
  res=reverse(map(length(a.c):-1:length(b.c)) do i
    if iszero(v[i]) c=zero(d)
    else c=v[i]*d
      view(v,i-length(b.c)+1:i) .-= c .* b.c
    end
    c
  end)
  Pol(res,a.v-b.v;check=false),Pol(v,a.v)
end

Base.div(a::Pol, b::Pol)=divrem(a,b)[1]
Base.:%(a::Pol, b::Pol)=divrem(a,b)[2]

Base.://(p::Pol,q::Pol)=isone(q.c[end]^2) ? p/q : p/(q//1)
Base.:/(p::Pol,q::T) where T=Pol(p.c/q,p.v;check=false)
function Base.:/(p::Pol,q::Pol)
  if q.c==[1] return shift(p,-q.v)
  elseif q.c==[-1] return shift(-p,-q.v)
  end
  r=divrem(p,q)
  if iszero(r[2]) return r[1] end
  error("divrem=$r division $p//$q not implemented")
end

Base.://(p::Pol,q::T) where T=iszero(p) ? p : Pol(p.c//q,p.v;check=false)
Base.://(p::T,q::Pol) where T=Pol(p)//q

"""
  gcd(p::Pol, q::Pol)
the  coefficients of  p and  q must  be elements  of a  field for gcd to be
type-stable

# Examples
```julia-repl
julia> gcd(2q+2,q^2-1)
Pol{Float64}: 1.0q+1.0

julia> gcd(q+1//1,q^2-1//1)
Pol{Rational{Int64}}: (1//1)q+1//1
```
"""
function Base.gcd(p::Pol,q::Pol)
  while !iszero(q)
    q=q*bestinv(q.c[end])
    (q,p)=(divrem(p,q)[2],q)
  end
  return p*bestinv(p.c[end])
end

function Base.inv(p::Pol)
  if length(p.c)!=1 throw(InexactError(:inv,typeof(p),p)) end
  Pol([bestinv(p.c[1])],-p.v;check=false)
end
end
