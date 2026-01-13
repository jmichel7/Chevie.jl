"""
`module Trunc`: truncated Laurent series

This module depends on the package `LaurentPolynomials`.

The main function to construct a truncated series is:
```julia-repl
julia> @Pol x
Pol{Int64}: x

julia> p=(x+x^2)^5
Pol{Int64}: x¹⁰+5x⁹+10x⁸+10x⁷+5x⁶+x⁵

julia> tp=Trunc(p,4)
Trunc(4): x⁵+5x⁶+10x⁷+10x⁸
```
The result here is a truncated series with 4 terms.
You can do various operations with series
```julia-repl
julia> inv(tp)
Trunc(4): x⁻⁵-5x⁻⁴+15x⁻³-35x⁻²

julia> inv(tp)*tp
Trunc(4): 1+0x+0x²+0x³

julia> tp*2
Trunc(4): 2x⁵+10x⁶+20x⁷+20x⁸
```
By default the variable name to print `Trunc`s is the same as for `Pol`s.
You can change that:
```julia-repl
julia> Truncs.varname=:q
:q

julia> tp
Trunc(4): q⁵+5q⁶+10q⁷+10q⁸

julia> length(tp) # the number of terms
4

julia> tp+1  # "loss of precision"
Trunc(4): 1+0q+0q²+0q³

julia> tp+Pol()^5 # ok here
Trunc(4): 2q⁵+5q⁶+10q⁷+10q⁸
```
As for `Pol`s, indexing gets the term of a given degree
```julia-repl
julia> tp[1]
0

julia> tp[6]
5
```
You can convert directly to `Trunc` a rational fraction:

```julia-repl
julia> a=(3Pol()^-2+1)/p
Frac{Pol{Int64}}: (x²+3)/(x¹²+5x¹¹+10x¹⁰+10x⁹+5x⁸+x⁷)

julia> Trunc(a,4)
Trunc(4): 3q⁻⁷-15q⁻⁶+46q⁻⁵-110q⁻⁴

julia> Trunc(a,4)*tp
Trunc(4): 3q⁻²+0q⁻¹+1+0q
```

We provide also convenience functions for `Mvp`s.
```julia-repl
julia> @Mvp x,y

julia> den=x-y
Mvp{Int64}: x-y

julia> frac=(x+y)/den
Frac{Mvp{Int64, Int64}}: (x+y)/(x-y)

julia> Truncs.varname=:x
:x

julia> tden=Trunc(den,4,:x) # convenience for Trunc(Pol(den,:x),4)
Trunc(4): -y+x+0x²+0x³

julia> tfrac=Trunc(frac,4,:x)
Trunc(4): -1-2y⁻¹x-2y⁻²x²-2y⁻³x³

julia> tden*tfrac
Trunc(4): y+x+0x²+0x³
```
"""
module Truncs
struct Trunc{T}
  c::Vector{T} # vector of length l: serie with l terms
  v::Int
end

using LaurentPolynomials, PuiseuxPolynomials

export Trunc

varname::Symbol=LaurentPolynomials.varname

Base.length(p::Trunc)=length(p.c)

function Base.:*(p::Trunc,q::Trunc)
  l=min(length(p),length(q))
  Trunc(map(j->sum(k->p[k+p.v]*q[j-k+q.v],0:j),0:l-1),p.v+q.v)
end
Base.:*(p::Trunc,a::Number)=Trunc(a*p.c,p.v)
Base.:*(a::Number,p::Trunc)=p*a
Base.:*(p::Trunc{T},a::T) where T =Trunc(a*p.c,p.v)
Base.:*(a::T,p::Trunc{T}) where T =p*a

function Base.:+(p::Trunc,q::Trunc)
  if p.v>q.v p,q=q,p end
  l=min(length(p),length(q)+q.v-p.v)
  Trunc(map(i->p[i]+q[i],p.v:p.v+l-1),p.v)
end
Base.:+(p::Trunc,q::Pol)=p+Trunc(q,length(p))
Base.:+(q::Pol,p::Trunc)=p+q
Base.:+(p::Trunc,q::Number)=p+Pol(q)
Base.:+(q::Number,p::Trunc)=p+q
Base.:-(p::Trunc,q::Pol)=p-Trunc(q,length(p))
Base.:-(q::Pol,p::Trunc)=-p+q
Base.:-(p::Trunc,q::Number)=p-Pol(q)
Base.:-(q::Number,p::Trunc)=-p+q

Base.:-(p::Trunc,q::Trunc)=p+(-q)

"""
`Trunc(p::Pol,i::Integer)` 

returns the truncated Laurent series for `p` with `i` terms
```julia-repl
julia> p=Pol([2,0,1],-1)
Pol{Int64}: x+2x⁻¹

julia> Trunc(p,4)
Trunc(4): 2x⁻¹+0+x+0x²
```
"""
Trunc(p::Pol,i::Integer)=Trunc(p[p.v:p.v+i-1],p.v)
LaurentPolynomials.Pol(p::Trunc)=Pol(p.c,p.v)

@inbounds Base.getindex(p::Trunc,i::Integer)=
  i in p.v:p.v+length(p.c)-1 ? p.c[i-p.v+1] : zero(p.c[1])

Base.:^(a::Trunc, n::Integer)=Base.power_by_squaring(a,n)
Base.:-(p::Trunc)=Trunc(-p.c,p.v)
Base.one(p::Trunc)=Trunc(vcat(1,fill(0,length(p)-1)),0)
Base.copy(p::Trunc)=Trunc(copy(p.c),p.v)

function Base.inv(q::Trunc)
  v=q.v
  c=1//q[v]
  q1=Trunc(-c*q.c[2:end],1)
  Trunc((c*(Trunc(Pol(1),length(q))+sum(j->q1^j,1:length(q)))).c,-v)
end

Base.:*(a::Frac{<:Mvp}, b::Trunc)=Trunc(a.*b.c,b.v)
Base.:*(b::Trunc,a::Frac{<:Mvp})=Trunc(a.*b.c,b.v)

Trunc(Q::Frac{<:Pol},i)=Trunc(numerator(Q),i)*inv(Trunc(denominator(Q),i))

Trunc(p::Mvp,i,var::Symbol)=Trunc(Pol(p,var),i)

function Trunc(Q::Frac{<:Mvp},i,var::Symbol)
  Trunc(numerator(Q),i,var)*inv(Trunc(denominator(Q),i,var))
end

Base.:(//)(p::Trunc,q::Trunc)=p*inv(q)

function Base.show(io::IO, ::MIME"text/plain", a::Trunc)
  if !haskey(io,:typeinfo) 
    print(io,"Trunc($(length(a))): ")
    io=IOContext(io,:typeinfo=>typeof(a))
  end
  show(io,a)
end

using LaurentPolynomials: format_coefficient, stringexp
function Base.show(io::IO,p::Trunc{T})where T
  if !get(io,:limit,false) && !get(io,:TeX,false) && !get(io,:naive,false)
    if ismonomial(p) && isone(p.c[1]) && p.v==1 && T==Int print(io,"Trunc()")
    else print(io,"Trunc(",p.c)
      if !iszero(p.v) print(io,",",p.v) end
      print(io,")")
    end
  else
    for deg in p.v:p.v+length(p)-1
      c=p[deg]
      c=repr(c; context=IOContext(io,:typeinfo=>typeof(c)))
      if !iszero(deg)
        c=format_coefficient(c,prod=get(io,:prod,false))*string(varname)
        c*=stringexp(get(io,:naive,false) ? stdout : io,deg)
      end
      if c[1]!='-' && deg!=p.v c="+"*c end
      print(io,c)
    end
  end
end
end
