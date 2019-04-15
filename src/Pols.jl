"""
 An implementation of univariate Laurent polynomials.
 A Pol contains two fields: its vector of coefficients, and its valuation.

# Examples
```julia-repl
julia> Pol(:q) # define string used for printing and set variable q
q

julia> Pol([1,2],0) # coefficients should have no leading or trailing zeroes.
2q+1

julia> p=Pol([1,2],-1)
2+q⁻¹

julia> valuation(p)
-1

julia> p=(q+1)^2
q²+2q+1

julia> degree(p)
2

julia> p(1//2) # a Pol is a callable object, where the call evaluates the Pol
9//4

julia> divrem(q^3+1,q+2) # changes coefficients to field elements
(1.0q²-2.0q+4.0, -7.0)

julia> divrem1(q^3+1,q+2) # keeps the ring, but needs second argument unitary
(q²-2q+4, -7)

julia> cyclotomic_polynomial(24) # the 24-th cyclotomic polynomial
q⁸-q⁴+1

```

see also the individual documentation of gcd.
"""
module Pols
export Pol, valuation, cyclotomic_polynomial, divrem1, shift, positive_part

using ..Gapjm # for degree

const var=Ref(:x)
varname(a::Symbol)=(var[]=a)
varname()=var[]

struct Pol{T}
  c::Vector{T}
  v::Int
end

function Polstrip(v::AbstractVector,val=0)
  b=findfirst(x->!iszero(x),v)
  if isnothing(b) return Pol(eltype(v)[],0) end
  l=findlast(x->!iszero(x),v)
  Pol(v[b:l],val+b-1)
end

function Pol(a)
  if iszero(a) Pol(typeof(a)[],0) end
  Pol([a],0)
end

function Pol(t::Symbol)
  varname(t)
  Base.eval(Main,:($t=Pol([1],1)))
end

Base.copy(p::Pol)=Pol(p.c,p.v)
Base.convert(::Type{Pol{T}},a::Number) where T=Pol([T(a)],0)
Base.convert(::Type{Pol{T}},p::Pol{T1}) where {T,T1}= T==T1 ? p : Pol(convert.(T,p.c),p.v)

Base.cmp(a::Pol,b::Pol)=cmp([a.c,a.v],[b.c,b.v])
Base.isless(a::Pol,b::Pol)=cmp(a,b)==-1

Gapjm.degree(p::Pol)=length(p.c)-1+p.v

valuation(p::Pol)=p.v

(p::Pol)(x)=horner(x,p.c)*x^p.v

shift(p::Pol,s)=Pol(p.c,p.v+s)

function positive_part(p::Pol)
  v=max(0,-p.v)
  if v==0 return p end
  Pol(p.c[1+v:end],p.v+v)
end

Base.:(==)(a::Pol, b::Pol)= a.c==b.c && a.v==b.v

Base.one(a::Pol)=Pol([one(eltype(a.c))],0)
Base.one(::Type{Pol{T}}) where T=Pol([one(T)],0)
Base.zero(::Type{Pol{T}}) where T=Pol(T[],0)
Base.zero(a::Pol)=Pol(empty(a.c),0)
Base.iszero(a::Pol)=length(a.c)==0
Base.transpose(a::Pol)=a

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
        mon=String(var[])
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
  else
    s="Pol($(p.c),$(p.v))"
  end
  print(io, (repl && !TeX) ? TeXstrip(s) : s)
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
  res=fill(zero(promote_type(T1,T2)),length(a.c)+length(b.c)-1)
  for i in eachindex(a.c), j in eachindex(b.c)
    res[i+j-1]+=a.c[i]*b.c[j]
  end
  Pol(res,a.v+b.v)
end

Base.:*(a::Pol, b::T) where T=iszero(b) ? zero(a) : Pol(a.c.*b,a.v)
Base.:*(b::T, a::Pol) where T=iszero(b) ? zero(a) : Pol(a.c.*b,a.v)

Base.:^(a::Pol, n::Real)= n>=0 ? Base.power_by_squaring(a,Int(n)) :
                                 Base.power_by_squaring(inv(a),-Int(n))

function Base.:+(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  d=b.v-a.v
  if d<0 return b+a end
  T=promote_type(T1,T2)
  res=zeros(T,max(length(a.c),d+length(b.c)))
@inbounds  res[eachindex(a.c)].=a.c
@inbounds  res[d.+eachindex(b.c)].+=b.c
  Polstrip(res,a.v)
end

Base.:+(a::Pol, b::T) where T=a+Pol(b)
Base.:+(b::T, a::Pol) where T=Pol(b)+a

Base.:-(a::Pol)=Pol(-a.c,a.v)
Base.:-(a::Pol, b::Pol)=a+(-b)
Base.:-(a::Pol, b::T) where T=a-Pol(b)
Base.:-(b::T, a::Pol) where T=Pol(b)-a

Base.div(a::Pol,b::Int)=Pol(div.(a.c,b),a.v)

"""
computes (p,q) such that a=p*b+q
"""
function Base.divrem(a::Pol, b::Pol)
  d=inv(b.c[end])
  T=typeof(a.c[end]*d)
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

"""
divrem when b unitary: does not change type
"""
function divrem1(a::Pol{T1}, b::Pol{T2})where {T1,T2}
  d=b.c[end]
  if d^2!=1 throw(InexactError) end
  T=promote_type(T1,T2)
  v=T.(a.c)
  res=T[]
  for i=length(a.c):-1:length(b.c)
    if iszero(v[i]) c=zero(d)
    else c=v[i]*d
         v[i-length(b.c)+1:i] .-= c .* b.c
    end
    pushfirst!(res,convert(T,c))
  end
  Pol(res,a.v-b.v),Polstrip(v,a.v)
end

Base.:/(p::Pol,q::T) where T=Pol(p.c/q,p.v)
function Base.://(p::Pol,q::Pol)
  if q.c==[1] return shift(p,q.v)
  elseif q.c==[-1] return shift(-p,q.v)
  end
  r=divrem1(p,q)
  if r[2]==1 return r[1] end
  error("division $p//$q not implemented")
end

Base.://(p::Pol,q::T) where T=Pol(p.c//q,p.v)

"""
  gcd(p::Pol, q::Pol)
  the coefficients of p and q must be elements of a field for
  gcd to be type-stable

# Examples
```julia-repl
julia> gcd(q+1,q^2-1)
1.0q+1.0

julia> gcd(q+1//1,q^2-1//1)
(1//1)q+1//1
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
  if length(p.c)>1 || !(p.c[1]^2==1) Throw(InexactError()) end
  Pol([p.c[1]],-p.v)
end

const cyclotomic_polynomial_dict=Dict(1=>Pol([-1,1],0))
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
end
