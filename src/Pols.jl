"""
 An implementation of univariate Laurent polynomials.
 A Pol contains two fields: its vector of coefficients, and its valuation.

# Examples
```julia-repl
julia> Pol(:q) # define string used for printing and set variable q
Pol{Int64}: q

julia> Pol([1,2]) # valuation is 0 if not specified
Pol{Int64}: 2q+1

julia> p=Pol([1,2],-1)
Pol{Int64}: 2+q⁻¹

julia> valuation(p)
-1

julia> p=(q+1)^2
Pol{Int64}: q²+2q+1

julia> degree(p)
2

julia> p(1//2) # a Pol is a callable object, where the call evaluates the Pol
9//4

julia> p[0], p[1], p[-1] # indexing gives the coefficients
(1, 2, 0)

julia> divrem(q^3+1,2q+1) # changes coefficients to field elements
(0.5q²-0.25q+0.125, 0.875)

julia> divrem(q^3+1,q+2)  # keeps the ring, but needs leading coeff ±1
(q²-2q+4, -7)

julia> cyclotomic_polynomial(24) # the 24-th cyclotomic polynomial
Pol{Int64}: q⁸-q⁴+1

```

see also the individual documentation of divrem, gcd.
"""
module Pols
using ..Util: format_coefficient, printTeX, divisors
export degree, valuation
export Pol, cyclotomic_polynomial, shift, positive_part, negative_part, bar

const varname=Ref(:x)

struct Pol{T}
  c::Vector{T}
  v::Int
  function Pol(c::AbstractVector{T},v::Integer=0;check=true)where T
    if check
      b=findfirst(!iszero,c)
      if isnothing(b) return new{T}(empty(c),0) end
      e=findlast(!iszero,c)
      if b!=1 || e!=length(c) return new{T}(c[b:e],v+b-1) end
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

Base.lastindex(p::Pol)=length(p.c)+p.v-1

Base.getindex(p::Pol{T},i) where T=i in p.v:p.v+length(p.c)-1 ? 
    p.c[i-p.v+1] : zero(T)

Base.copy(p::Pol)=Pol(p.c,p.v;check=false)
Base.convert(::Type{Pol{T}},a::Number) where T=iszero(a) ? zero(Pol{T}) :
        Pol([T(a)];check=false)
Base.convert(::Type{Pol},a::Number)=iszero(a) ? zero(Pol{typeof(a)}) :
        Pol([a];check=false)
(::Type{Pol{T}})(a::Number) where T=convert(Pol{T},a)
Base.convert(::Type{Pol{T}},p::Pol{T1}) where {T,T1}= T==T1 ? p :
        Pol(convert.(T,p.c),p.v;check=false)
(::Type{Pol{T}})(p::Pol) where T=convert(Pol{T},p)

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
  if !isone(length(p.c)) || !iszero(p.v) throw(InexactError(:convert,T,p)) end
  convert(T,p.c[1]) 
end

(::Type{T})(p::Pol) where T<:Number=convert(T,p)

Base.cmp(a::Pol,b::Pol)=cmp([a.c,a.v],[b.c,b.v])
Base.isless(a::Pol,b::Pol)=cmp(a,b)==-1
Base.hash(a::Pol, h::UInt)=hash(a.v,hash(a.c,h))

degree(p::Pol)=length(p.c)-1+p.v

valuation(p::Pol)=p.v

(p::Pol{T})(x) where T=iszero(p) ? zero(T) : evalpoly(x,p.c)*x^p.v

# efficient p↦ qˢ p
shift(p::Pol,s)=Pol(p.c,p.v+s;check=false)

# degree ≥0
function positive_part(p::Pol)
  v=max(0,-p.v)
  if v==0 return p end
  Pol(p.c[1+v:end],p.v+v)
end

# degree ≤0
function negative_part(p::Pol)
  if degree(p)<=0 return p end
  Pol(p.c[1:1-p.v],p.v)
end

# q↦ q⁻¹ on p
bar(p::Pol)=Pol(reverse(p.c),-degree(p);check=false)

Base.:(==)(a::Pol, b::Pol)= a.c==b.c && a.v==b.v
Base.:(==)(a::Pol,b)= a==Pol(b)
Base.:(==)(b,a::Pol)= a==Pol(b)

Base.one(a::Pol)=Pol([one(eltype(a.c))];check=false)
Base.one(::Type{Pol{T}}) where T=Pol([one(T)];check=false)
Base.one(::Type{Pol})=Pol([1];check=false)
Base.zero(::Type{Pol{T}}) where T=Pol(T[];check=false)
Base.zero(::Type{Pol})=Pol(Int[];check=false)
Base.zero(a::Pol)=Pol(empty(a.c);check=false)
Base.iszero(a::Pol)=length(a.c)==0
#Base.transpose(a::Pol)=a # next 3 stupid stuff to make inv using LU work
Base.adjoint(a::Pol)=a
Base.abs(p::Pol)=p
Base.conj(p::Pol)=Pol(conj.(p.c),p.v;check=false)

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
  if get(io,:limit,false) || get(io,:TeX,false)
    s=join(map(reverse(collect(enumerate(p.c))))do (i,c)
      if iszero(c) return "" end
      c=sprint(show,c; context=io)
      deg=i+p.v-1
      if !iszero(deg) 
        mon=String(varname[])
        if deg!=1 mon*="^{$deg}" end
        c=format_coefficient(c)*mon
      end
      if isempty(c) || c[1]!='-' c="+"*c end
      c
    end)
    if s=="" s="0"
    elseif  s[1]=='+' s=s[2:end]
    end
    printTeX(io,s)
  else
    print(io,"Pol(",p.c,",",p.v,")")
  end
end

function Base.:*(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  if iszero(a) || iszero(b) return zero(a) end
  res=fill(zero(T1)*zero(T2),length(a.c)+length(b.c)-1)
  for i in eachindex(a.c), j in eachindex(b.c)
@inbounds res[i+j-1]+=a.c[i]*b.c[j]
  end
  Pol(res,a.v+b.v;check=false)
end

Base.:*(a::Pol, b::Number)=Pol(a.c.*b,a.v)
Base.:*(a::Pol{T}, b::T) where T=Pol(a.c.*b,a.v)
Base.:*(b::Number, a::Pol)=a*b
Base.:*(b::T, a::Pol{T}) where T=a*b

Base.:^(a::Pol, n::Real)= n>=0 ? Base.power_by_squaring(a,Int(n)) :
                                 Base.power_by_squaring(inv(a),-Int(n))

function Base.:+(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  d=b.v-a.v
  if d<0 return b+a end
  res=fill(zero(T1)+zero(T2),max(length(a.c),d+length(b.c)))
@inbounds view(res,eachindex(a.c)).=a.c
@inbounds view(res,d+eachindex(b.c)).+=b.c
  Pol(res,a.v)
end

Base.:+(a::Pol, b::Number)=a+Pol(b)
Base.:+(b::Number, a::Pol)=Pol(b)+a
Base.:-(a::Pol)=Pol(-a.c,a.v;check=false)
Base.:-(a::Pol, b::Pol)=a+(-b)
Base.:-(a::Pol, b::Number)=a-Pol(b)
Base.:-(b::Number, a::Pol)=Pol(b)-a
Base.div(a::Pol,b::Int)=Pol(div.(a.c,b),a.v;check=false)

"""
`divrem(a::Pol, b::Pol)`

computes `(p,q)` such that `a=p*b+q`
When the leading coefficient of b is ±1 does not change type
"""
function Base.divrem(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  if iszero(b) throw(DivideError) end
  c=b.c[end]
  if isone(c^2) d=c else d=inv(c) end
  T=promote_type(T1,T2,typeof(d))
  v=convert(Vector{T},copy(a.c))
  res=T[]
  for i=length(a.c):-1:length(b.c)
    if iszero(v[i]) c=zero(d)
    else c=v[i]*d
         v[i-length(b.c)+1:i] .-= c .* b.c
    end
    pushfirst!(res,c)
  end
  Pol(res,a.v-b.v;check=false),Pol(v,a.v)
end

Base.div(a::Pol, b::Pol)=divrem(a,b)[1]

Base.:/(p::Pol,q::Pol)=p//q
Base.:/(p::Pol,q::T) where T=Pol(p.c/q,p.v;check=false)
function Base.://(p::Pol,q::Pol)
  if q.c==[1] return shift(p,-q.v)
  elseif q.c==[-1] return shift(-p,-q.v)
  end
  r=divrem(p,q)
  if iszero(r[2]) return r[1] end
  error("r=$r division $p//$q not implemented")
end

Base.://(p::Pol,q::T) where T=Pol(p.c//q,p.v;check=false)
Base.://(p::T,q::Pol) where T=Pol(p)//q

"""
  gcd(p::Pol, q::Pol)
the  coefficients of  p and  q must  be elements  of a  field for gcd to be
type-stable

# Examples
```julia-repl
julia> gcd(q+1,q^2-1)
Pol{Float64}: 1.0q+1.0

julia> gcd(q+1//1,q^2-1//1)
Pol{Rational{Int64}}: (1//1)q+1//1
```
"""
function Base.gcd(p::Pol,q::Pol)
  while !iszero(q)
    q=q/q.c[end]
    (q,p)=(divrem(p,q)[2],q)
  end
  return p/p.c[end]
end

function Base.inv(p::Pol)
  if length(p.c)!=1 throw(InexactError(:inv,Int,p)) end
  if p.c[1]^2==1 return Pol([p.c[1]],-p.v;check=false) end
  Pol([inv(p.c[1])],-p.v;check=false)
end

const cyclotomic_polynomial_dict=Dict(1=>Pol([-1,1],0;check=false))
"""
`cyclotomic_polynomial(n)`
 
returns the `n`-th cyclotomic polynomial.
 
```julia-repl
julia> cyclotomic_polynomial(5)
Pol{Int64}: q⁴+q³+q²+q+1
```
 
The  computed  cyclotomic  polynomials  are  cached  in  the global `Dict ̀
`Pols.cyclotomic_polynomial_dict`
"""
function cyclotomic_polynomial(n::Integer)
  get!(cyclotomic_polynomial_dict,n) do
    v=fill(0,n+1);v[1]=-1;v[n+1]=1;res=Pol(v,0;check=false)
    for d in divisors(n)
      if d!=n
        res,foo=divrem(res,cyclotomic_polynomial(d))
      end
    end
    res
  end
end
end
