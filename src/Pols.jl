"""
 An implementation of univariate Laurent polynomials.
 A Pol contains two fields: its vector of coefficients, and its valuation.

# Examples
```julia-repl
julia> Pol(:q) # define string used for printing and set variable q
Pol{Int64}: q

julia> Pol([1,2],0) # coefficients should have no leading or trailing zeroes.
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

julia> divrem(q^3+1,q+2) # changes coefficients to field elements
(1.0q²-2.0q+4.0, -7.0)

julia> divrem1(q^3+1,q+2) # keeps the ring, but needs leading coeff divides
(q²-2q+4, -7)

julia> cyclotomic_polynomial(24) # the 24-th cyclotomic polynomial
Pol{Int64}: q⁸-q⁴+1

```

see also the individual documentation of divrem, divrem1, gcd.
"""
module Pols
using ..Util: bracket_if_needed, fromTeX, divisors
using ..Cycs: Cyc, Root1
import Gapjm: root, degree, valuation
# to use as a stand-alone module comment above line and uncomment next
# export root, degree, valuation
export Pol, cyclotomic_polynomial, divrem1, shift, positive_part,
  negative_part, bar, isunit, improve_type

const varname=Ref(:x)

struct Pol{T}
  c::Vector{T}
  v::Int
end

Pol(a::Number)=convert(Pol,a)

function Pol(t::Symbol)
  varname[]=t
  Base.eval(Main,:($t=Pol([1],1)))
end

Base.broadcastable(p::Pol)=Ref(p)

Base.lastindex(p::Pol)=length(p.c)+p.v-1

Base.getindex(p::Pol{T},i) where T=i in p.v:p.v+length(p.c)-1 ? 
    p.c[i-p.v+1] : zero(T)

function Polstrip(v::AbstractVector,val=0)
  b=findfirst(!iszero,v)
  if isnothing(b) return zero(Pol{eltype(v)}) end
  Pol(v[b:findlast(!iszero,v)],val+b-1)
end

Base.copy(p::Pol)=Pol(p.c,p.v)
Base.convert(::Type{Pol{T}},a::Number) where T=iszero(a) ? zero(Pol{T}) : Pol([T(a)],0)
Base.convert(::Type{Pol},a::Number)=iszero(a) ? zero(Pol{typeof(a)}) : Pol([a],0)
(::Type{Pol{T}})(a::Number) where T=convert(Pol{T},a)
Base.convert(::Type{Pol{T}},p::Pol{T1}) where {T,T1}= T==T1 ? p : Pol(convert.(T,p.c),p.v)
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
shift(p::Pol,s)=Pol(p.c,p.v+s)

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
bar(p::Pol)=Pol(reverse(p.c),-degree(p))

Base.:(==)(a::Pol, b::Pol)= a.c==b.c && a.v==b.v

Base.one(a::Pol)=Pol([one(eltype(a.c))],0)
Base.one(::Type{Pol{T}}) where T=Pol([one(T)],0)
Base.one(::Type{Pol})=Pol([1],0)
Base.zero(::Type{Pol{T}}) where T=Pol(T[],0)
Base.zero(::Type{Pol})=Pol(Int[],0)
Base.zero(a::Pol)=Pol(empty(a.c),0)
Base.iszero(a::Pol)=length(a.c)==0
Base.transpose(a::Pol)=a # next 3 stupid stuff to make inv using LU work
Base.adjoint(a::Pol)=a
Base.abs(p::Pol)=p
Base.conj(p::Pol)=Pol(conj.(p.c),p.v)

function improve_type(m::Array)
  if !isempty(m) m=convert.(reduce(promote_type,typeof.(m)),m) end
  if all(isinteger,m) m=Int.(m)
  elseif first(m) isa Cyc && all(x->x.n==1,m) m=Rational.(m)
  elseif first(m) isa Pol && all(x->all(isinteger,x.c),m) m=convert.(Pol{Int},m)
  end
  m
end

function improve_type(mm::Vector{<:Array})
  mm=improve_type.(mm)
  T=reduce(promote_type,eltype.(mm))
  map(x->convert.(T,x),mm)
end

function Base.show(io::IO, ::MIME"text/html", a::Pol)
  print(io, "\$")
  show(IOContext(io,:TeX=>true),a)
  print(io, "\$")
end

function Base.show(io::IO, ::MIME"text/plain", a::Pol)
  print(io,typeof(a),": ")
  show(io,a)
end

function Base.show(io::IO,p::Pol)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if repl||TeX
    s=join(map(reverse(eachindex(p.c)))do i
      c=p.c[i]
      if iszero(c) return "" end
      c=sprint(show,c; context=io)
      deg=i+p.v-1
      if !iszero(deg) 
        mon=String(varname[])
        if deg!=1
           if repl || TeX mon*=(1<deg<10 ? "^$deg" : "^{$deg}")
           else mon*= "^$deg"
           end
        end
        if c=="1" c="" elseif c=="-1" c="-" end
        c=bracket_if_needed(c) 
        c*=mon
      end
      if c[1]!='-' c="+"*c end
      c
    end)
    if s=="" s="0"
    elseif  s[1]=='+' s=s[2:end]
    end
    print(io, fromTeX(io,s))
  else
    print(io,"Pol(",p.c,",",p.v,")")
  end
end

#function Base.:*(a::Pol, b::Pol)
#  if iszero(a) || iszero(b) return zero(a) end
#  res=map(1:length(a.c)+length(b.c)-1)do i
#@inbounds sum(j->a.c[j]*b.c[i+1-j],i+1-min(i,length(b.c)):min(i,length(a.c)))
#  end
#  Pol(res,a.v+b.v)
#end

# stupid code is better
function Base.:*(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  if iszero(a) || iszero(b) return zero(a) end
  res=fill(zero(T1)*zero(T2),length(a.c)+length(b.c)-1)
  for i in eachindex(a.c), j in eachindex(b.c)
@inbounds res[i+j-1]+=a.c[i]*b.c[j]
  end
  Pol(res,a.v+b.v)
end

Base.:*(a::Pol, b::Number)=Polstrip(a.c.*b,a.v)
Base.:*(b::Number, a::Pol)=a*b

Base.:^(a::Pol, n::Real)= n>=0 ? Base.power_by_squaring(a,Int(n)) :
                                 Base.power_by_squaring(inv(a),-Int(n))

function Base.:+(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  d=b.v-a.v
  if d<0 return b+a end
  res=fill(zero(T1)+zero(T2),max(length(a.c),d+length(b.c)))
@inbounds res[eachindex(a.c)].=a.c
@inbounds res[d.+eachindex(b.c)].+=b.c
  Polstrip(res,a.v)
end

Base.:+(a::Pol, b::Number)=a+Pol(b)
Base.:+(b::Number, a::Pol)=Pol(b)+a
Base.:-(a::Pol)=Pol(-a.c,a.v)
Base.:-(a::Pol, b::Pol)=a+(-b)
Base.:-(a::Pol, b::Number)=a-Pol(b)
Base.:-(b::Number, a::Pol)=Pol(b)-a
Base.div(a::Pol,b::Int)=Pol(div.(a.c,b),a.v)

"""
computes (p,q) such that a=p*b+q
"""
function Base.divrem(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  if iszero(b) throw(DivideError) end
  T=promote_type(T1,T2)
  d=inv(T(b.c[end]))
  T=promote_type(T1,T2,typeof(d))
  v=T.(a.c)
  res=T[]
  for i=length(a.c):-1:length(b.c)
    if iszero(v[i]) c=zero(d)
    else c=v[i]*d
         v[i-length(b.c)+1:i] .-= c .* b.c
    end
    pushfirst!(res,c)
  end
  Pol(res,a.v-b.v),Polstrip(v,a.v)
end

Base.div(a::Pol, b::Pol)=divrem1(a,b)[1]

"""
divrem when b divides in ring: does not change type
"""
function divrem1(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  if iszero(b) throw(DivideError) end
  d=b.c[end]
  T=promote_type(T1,T2)
  if !(T<:Integer) return divrem(a,b) end
  v=T.(a.c)
  res=T[]
  for i=length(a.c):-1:length(b.c)
    if iszero(v[i]) pushfirst!(res,0)
    else r=divrem(v[i],d)
      if !iszero(r[2]) throw(InexactError) end
      v[i-length(b.c)+1:i] .-= r[1] .* b.c
      pushfirst!(res,r[1])
    end
  end
  Pol(res,a.v-b.v),Polstrip(v,a.v)
end

Base.:/(p::Pol,q::Pol)=p//q
Base.:/(p::Pol,q::T) where T=Pol(p.c/q,p.v)
function Base.://(p::Pol,q::Pol)
  if q.c==[1] return shift(p,-q.v)
  elseif q.c==[-1] return shift(-p,-q.v)
  end
  r=divrem1(p,q)
  if iszero(r[2]) return r[1] end
  error("r=$r division $p//$q not implemented")
end

Base.://(p::Pol,q::T) where T=Pol(p.c//q,p.v)
Base.://(p::T,q::Pol) where T=Pol(p)//q

"""
  gcd(p::Pol, q::Pol)
  the coefficients of p and q must be elements of a field for
  gcd to be type-stable

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

isunit(p::Pol)=length(p.c)==1 && p.c[1]^2==1

function Base.inv(p::Pol)
  if length(p.c)>1 throw(InexactError(:inv,Int,p)) end
  if p.c[1]^2==1 return Pol([p.c[1]],-p.v) end
  r=Root1(p.c[1])
  if isnothing(r) throw(InexactError(:inv,Int,p)) end
  Pol([inv(p.c[1])],-p.v)
end

const cyclotomic_polynomial_dict=Dict(1=>Pol([-1,1],0))
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
    v=fill(0,n+1);v[1]=-1;v[n+1]=1;res=Pol(v,0)
    for d in divisors(n)
      if d!=n
        res,foo=divrem1(res,cyclotomic_polynomial(d))
      end
    end
    res
  end
end

function root(x::Pol,n::Number=2)
  n=Int(n)
  if length(x.c)>1 || !iszero(x.v%n)
    error("root($(repr(x;context=:limit=>true)),$n) not implemented") 
  end
  if isempty(x.c) return x end
  Pol([root(x.c[1],n)],div(x.v,n))
end
end
