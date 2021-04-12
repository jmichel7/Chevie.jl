"""
The  type `Mvp`  ("multivariate polynomials")  implemented here is "Puiseux
polynomials",  that  is  linear  combinations  of  monomials  of  the  type
`x₁^{a₁}…  xₙ^{aₙ}` where `xᵢ`  are variables and  `aᵢ` are exponents which
can  be  arbitrary  rational  numbers  (we  need  Puiseux  polynomial  with
cyclotomic  coefficients as splitting fields of cyclotomic Hecke algebras).
Some  functions described below work  only with polynomials where variables
are  raised to integral powers;  we will refer to  such objects as "Laurent
polynomials"; some functions require further that variables are raised only
to  positive powers: we refer then to "true polynomials". We implement also
multivariate  rational  fractions  (type  `Mvrf`)  where  the numerator and
denominator are true polynomials.

`@Mvp x₁,…,xₙ`

assigns   to  the  Julia  names   `xᵢ`  indeterminates  suitable  to  build
multivariate  polynomials  or  rational  fractions.  `Mvp(:x1)` creates the
indeterminate `x1` without assigning a Julia variable.

```julia-repl
julia> @Mvp x,y

julia> (x+y^-1)^3
Mvp{Int64}: x³+3x²y⁻¹+3xy⁻²+y⁻³
```

`Mvp(x::Number)`   returns  the  constant   multivariate  polynomial  whose
constant term is `x`.

```julia-repl
julia> degree(Mvp(1))
0
```

Only  monomials and one-term `Mvp`s can  be raised to a non-integral power;
the  `Mvp` with one  term `cm` which  is `c` times  the monomial `m` can be
raised  to a fractional power of denominator `d` if and only if `root(c,d)`
is  defined (this is  equivalent to `c^{1//d}`  but one may  want to define
`root`  differently;  for  instance,  in  my  other package `Cycs` I define
square  roots of rationals as cyclotomics;  I also have implement in `Cycs`
arbitrary roots of roots of unity).

```julia-repl
julia> (4x)^(1//2)
Mvp{Int64,Rational{Int64}}: 2x½

julia> (2.0x)^(1//2)
Mvp{Float64,Rational{Int64}}: 1.4142135623730951x½

julia> root(2.0x)
Mvp{Float64,Rational{Int64}}: 1.4142135623730951x½
```

Using the package `Cycs`:

```julia-repl
julia> (2x)^(1//2)
Mvp{Cyc{Int64},Rational{Int64}}: √2x½

julia> (E(3)*x)^(2//3)
Mvp{Cyc{Int64},Rational{Int64}}: ζ₉²x⅔
```

One can divide an `Mvp` by another when the division is exact

```julia-repl
julia> exactdiv(x^2-y^2,x-y)
Mvp{Int64}: x+y

julia> (x+y)/(2x^2)   # or by a monomial
Mvp{Rational{Int64}}: (1//2)x⁻¹+(1//2)x⁻²y

julia> (x+y)/(x-y)    # otherwise one gets a rational fraction.
Mvrf{Int64}: (x+y)/(x-y)
```

Raising  a non-monomial  Laurent polynomial  to a  negative power returns a
rational fraction.

```julia-repl
julia> (x+1)^-2
Mvrf{Int64}: 1/(x²+2x+1)
```

One can compute the value `Mvp` or a `Ratfrac` when setting some variables 
by the call syntax:

```julia-repl
julia> p=x+y
Mvp{Int64}: x+y

julia> p(x=2)
Mvp{Int64}: y+2

julia> p(x=2,y=x)
Mvp{Int64}: x+2

julia> ((x+y)/(x-y))(x=y+1)
Mvp{Int64}: 2y+1
```

Note  that the value of  an `Mvp` is always  an `Mvp`, for consistency. The
function  `scal` converts  a constant  `Mvp` to  that constant (and returns
`nothing` if the argument is not constant:

```julia-repl
julia> p(x=1,y=2)
Mvp{Int64}: 3

julia> scal(p(x=1,y=2))
3

julia> v=(x^(1//2))(x=2.0)
Mvp{Float64}: 1.4142135623730951

julia> scal(v)
1.4142135623730951
```

see the functions `coefficient`, `coefficients, `Pol` to take apart `Mvp`s.
`Mvrf` are dissected using `numerator` and `denominator`.

Despite  the degree of generality of our  polynomials, the speed is not too
shabby. For the Fateman test f(f+1) where f=(1+x+y+z+t)^15, we take 4sec.
According to the Nemo paper, Sagemath takes 10sec and Nemo takes 1.6sec.
"""
module Mvps
# benchmark: (x+y+z)^3     2.3μs 48 alloc
using ..ModuleElts: ModuleElt, ModuleElts
using ..Util: Util, exactdiv, fromTeX, ordinal, printTeX, bracket_if_needed
using ..Pols: Pols, Pol, srgcd, positive_part, negative_part, bar, derivative,
              valuation, degree

#import Gapjm: coefficients, valuation
import ..Cycs: root
using ..Cycs: Cyc
# to use as a stand-alone module comment above line, uncomment next, and
# define root for the coefficients you want (at least root(1,n)=1)
#export root
export coefficients, coefficient
export Mvp, Monomial, @Mvp, variables, value, scal, laurent_denominator, Mvrf
#------------------ Monomials ---------------------------------------------
struct Monomial{T}
  d::ModuleElt{Symbol,T}   
end

Monomial(a::Pair...)=Monomial(ModuleElt(a...;check=false)) # check about check
Monomial()=one(Monomial{Int})
Monomial(v::Symbol...)=Monomial(ModuleElt(a=>1 for a in v))

Base.convert(::Type{Monomial{T}},v::Symbol) where T=Monomial(v=>T(1))
Base.convert(::Type{Monomial{T}},m::Monomial{N}) where {T,N}= 
  Monomial(convert(ModuleElt{Symbol,T},m.d))

# we need a promote rule to mix fractional momonials with others in an Mvp
Base.promote_rule(::Type{Monomial{N}},::Type{Monomial{M}}) where {N,M}=
  Monomial{promote_type(N,M)}

Base.:*(a::Monomial, b::Monomial)=Monomial(a.d+b.d)
Base.isone(a::Monomial)=iszero(a.d)
#Base.iszero(a::Monomial)=false
Base.one(::Type{Monomial{N}}) where N=Monomial(zero(ModuleElt{Symbol,N}))
Base.one(m::Monomial)=Monomial(zero(m.d))
Base.inv(a::Monomial)=Monomial(-a.d)
Util.exactdiv(a::Monomial, b::Monomial)=a*inv(b)
Base.:/(a::Monomial, b::Monomial)=a*inv(b)
Base.://(a::Monomial, b::Monomial)=a*inv(b)
Base.:^(x::Monomial,p)=Monomial(x.d*p)
#Base.getindex(a::Monomial,k)=getindex(a.d,k)

function Base.show(io::IO,m::Monomial)
  replorTeX=get(io,:TeX,false) || get(io,:limit,false)
  if isone(m) return end
  start=true
  for (v,d) in m.d
    if !(start || replorTeX) print(io,"*") end
    print(io,replorTeX ? fromTeX(io,string(v)) : string(v))
    if !isone(d) 
      if isone(denominator(d)) d=numerator(d) end
      if replorTeX 
        if d isa Integer printTeX(io,"^{$d}") 
        else printTeX(io,"^{\\frac{$(numerator(d))}{$(denominator(d))}}") 
        end
      elseif d isa Integer print(io,"^$d")
      else print(io,"^($d)")
      end
    end
    start=false
  end
end

function Base.show(io::IO, ::MIME"text/plain", m::Monomial)
  if !haskey(io,:typeinfo) print(io,typeof(m),":") end
  show(io,m)
end

"""
`isless(a::Monomial,b::Monomial)`

For  our implementation of `Mvp`s to  work, `isless` must define a monomial
order  (that is, for any monomial `m` we have `a<b => a*m<b*m`). By default
we define `a<b` if the first variable in `a/b` occurs to a positive power.
"""
function Base.isless(a::Monomial, b::Monomial)
  for ((va,pa),(vb,pb)) in zip(a.d,b.d)
    if va!=vb return va<vb ? pa>0 : pb<0 end
    if pa!=pb return isless(pb,pa) end
  end
  la=length(a.d)
  lb=length(b.d)
  @inbounds la>lb ? last(a.d.d[lb+1])>0 : lb>la ? last(b.d.d[la+1])<0 : false
end

Base.:(==)(a::Monomial, b::Monomial)=a.d==b.d

Base.gcd(m::Monomial,n::Monomial)=Monomial(ModuleElts.merge2(min,m.d,n.d))
Base.gcd(v::AbstractArray{<:Monomial})=reduce(gcd,v)

Base.lcm(m::Monomial,n::Monomial)=Monomial(ModuleElts.merge2(max,m.d,n.d))
Base.lcm(v::AbstractArray{<:Monomial})=reduce(lcm,v)

Base.hash(a::Monomial, h::UInt)=hash(a.d,h)

Pols.degree(m::Monomial)=sum(values(m.d))
Pols.degree(m::Monomial,var::Symbol)=m.d[var]

function root(m::Monomial,n::Integer=2)
 if all(x->iszero(x%n),last.(m.d)) Monomial((k=>div(v,n) for (k,v) in m.d)...)
 else Monomial((k=>v//n for (k,v) in m.d)...)
 end
end

#------------------------------------------------------------------------------
struct Mvp{T,N} # N=type of exponents T=type of coeffs
  d::ModuleElt{Monomial{N},T}
end

Mvp(a::Pair...;c...)=Mvp(ModuleElt(a...;c...))
Mvp(;c...)=zero(Mvp{Int,Int}) # for some calls to map() to work

macro Mvp(t) # @Mvp x,y,z defines variables to be Mvp
  if t isa Expr
    for v in t.args
      Base.eval(Main,:($v=Mvp($(Core.QuoteNode(Symbol(v))))))
    end
  elseif t isa Symbol
    Base.eval(Main,:($t=Mvp($(Core.QuoteNode(t)))))
  end
end

Base.broadcastable(p::Mvp)=Ref(p)
Base.cmp(a::Mvp,b::Mvp)=cmp(a.d,b.d)
Base.isless(a::Mvp,b::Mvp)=cmp(a,b)==-1
Base.hash(a::Mvp,h::UInt)=hash(a.d,h)

function Base.show(io::IO, ::MIME"text/html", a::Mvp)
  print(IOContext(io,:TeX=>true),"\$",a,"\$")
end

function Base.show(io::IO, ::MIME"text/plain", a::Mvp{T,N}) where{T,N}
  if !haskey(io,:typeinfo) print(io,N==Int ? "Mvp{$T}" : "Mvp{$T,$N}",": ") end
  show(io,a)
end

# necessary if called when already showing a ModuleElt
Base.show(io::IO, x::Mvp)=show(IOContext(io,:showbasis=>nothing),x.d)

Base.zero(p::Mvp)=Mvp(zero(p.d))
Base.zero(::Type{Mvp{T,N}}) where {T,N}=Mvp(zero(ModuleElt{Monomial{N},T}))
Base.one(::Type{Mvp{T,N}}) where {T,N}=Mvp(one(Monomial{N})=>one(T))
Base.one(::Type{Mvp{T}}) where T=one(Mvp{T,Int})
Base.one(::Type{Mvp})=Mvp(1)
Base.one(p::Mvp{T,N}) where {T,N}=one(Mvp{T,N})
Base.copy(p::Mvp)=Mvp(p.d)
Base.iszero(p::Mvp)=iszero(p.d)
Base.convert(::Type{Mvp},a::Number)=convert(Mvp{typeof(a),Int},a)
Base.convert(::Type{Mvp{T,N}},a::Number) where {T,N}=iszero(a) ? 
 zero(Mvp{T,N}) : Mvp(one(Monomial{N})=>convert(T,a))
(::Type{Mvp{T,N}})(a::Number) where {T,N}=convert(Mvp{T,N},a)
(::Type{Mvp{T,N}})(a::Mvp) where {T,N}=convert(Mvp{T,N},a)
Base.convert(::Type{Mvp{T,N}},a::Mvp{T1,N1}) where {T,T1,N,N1}=
  Mvp(convert(ModuleElt{Monomial{N},T},a.d))
Base.convert(::Type{Mvp},x::Mvp)=x
Base.convert(::Type{Mvp},v::Symbol)=Mvp(Monomial(v)=>1)
Mvp(x::Symbol)=convert(Mvp,x)
Mvp(x::Number)=convert(Mvp,x)
Mvp(x::Mvp)=convert(Mvp,x)
# stupid stuff to make LU work
Base.adjoint(a::Mvp)=conj(a)
Base.abs(a::Mvp)=a

Base.:(==)(a::Mvp, b::Mvp)=a.d==b.d
Base.:(==)(a::Mvp,x::Number)=a==Mvp(x)
Base.:(==)(x::Number,a::Mvp)=a==Mvp(x)

function Base.convert(::Type{T},a::Mvp) where T<:Number
  if iszero(a) return zero(T) end
  if length(a.d)>1 || !isone(first(first(a.d)))
      throw(InexactError(:convert,T,a)) 
  end
  convert(T,last(first(a.d)))
end
(::Type{T})(a::Mvp) where T<: Number=convert(T,a)

Base.isinteger(p::Mvp)=iszero(p) || (isone(length(p.d)) &&
             isone(first(first(p.d))) && isinteger(last(first(p.d))))

# we need a promote rule to handle Vectors of Mvps of different types
Base.promote_rule(::Type{Mvp{T1,N1}},::Type{Mvp{T2,N2}}) where {T1,T2,N1,N2} =
  Mvp{promote_type(T1,T2),promote_type(N1,N2)}

# and this rule to handle Vectors mixing Mvps with numbers
Base.promote_rule(::Type{Mvp{T1,N}},::Type{T2}) where {T1,N,T2<:Number} =
  Mvp{promote_type(T1,T2),N}

Base.:+(a::Mvp, b::Mvp)=Mvp(a.d+b.d) # ModuleElts.+ takes care of promotion
Base.:+(a::Number, b::Mvp)=Mvp(a)+b
Base.:+(a::Mvp, b::Number)=b+a

Base.:-(a::Mvp)=Mvp(-a.d)
Base.:-(a::Mvp, b::Mvp)=a+(-b) # not Mvp(a.d-b.d) because 0-x != x (merge)
Base.:-(a::Mvp, b::Number)=a-Mvp(b)
Base.:-(b::Number, a::Mvp)=Mvp(b)-a

Base.:*(a::Number, b::Mvp)=Mvp(b.d*a)
Base.:*(b::Mvp, a::Number)=a*b
# we use we have a monomial order so there is no order check in next line
Base.:*(a::Monomial, b::Mvp)=Mvp(ModuleElt(m*a=>c for (m,c) in b.d;check=false))
Base.:*(b::Mvp,a::Monomial)=a*b
function Base.:*(a::Mvp, b::Mvp)
  if length(a.d)>length(b.d) a,b=(b,a) end
  if iszero(a) return a end
  let b=b # needed !!!!
    sum(b*m*c for (m,c) in a.d)
  end
end

Base.:(//)(a::Mvp, b::Number)=Mvp(ModuleElt(m=>c//b for (m,c) in a.d;check=false))
Base.:(/)(a::Mvp, b::Number)=Mvp(ModuleElt(m=>c/b for (m,c) in a.d;check=false))

"""
`conj(p::Mvp)` acts on the coefficients of `p`

```julia-repl
julia> @Mvp x;conj(im*x+1)
Mvp{Complex{Int64}}: (0 - 1im)x+1 + 0im
```
"""
Base.conj(a::Mvp)=Mvp(ModuleElt(m=>conj(c) for (m,c) in a.d;check=false))

function Base.:^(x::Mvp, p::Union{Integer,Rational})
  if isinteger(p) p=Int(p) end
  if iszero(p) return one(x)
  elseif iszero(x) return x
  elseif !isinteger(p) return root(x,denominator(p))^numerator(p) 
  elseif isone(p) return x 
  elseif length(x.d)==1
    (m,c)=first(x.d)
    return Mvp(m^p=>c^p)
  else return p>=0 ? Base.power_by_squaring(x,p) :
                     Base.power_by_squaring(inv(x),-p)
  end
end

"""
The `degree` of a monomial is the sum of  the exponents of the variables.
The `degree` of an `Mvp` is the largest degree of a monomial.

```julia-repl
julia> a=x^2+x*y
Mvp{Int64}: x²+xy

julia> degree(a)
2
```

With  second argument a  variable name, `degree`  returns the degree of the
polynomial in that variable.

```julia-repl
julia> degree(a,:y)
1

julia> degree(a,:x)
2
```

"""
Pols.degree(m::Mvp)=iszero(m) ? 0 : maximum(degree.(keys(m.d)))
Pols.degree(m::Mvp,v::Symbol)=iszero(m) ? 0 : maximum(degree.(keys(m.d),v))

"""
The `valuation` of an `Mvp` is the minimal degree of a monomial.


```julia-repl
julia> @Mvp x,y; a=x^2+x*y
Mvp{Int64}: x²+xy

julia> valuation(a)
2
```

With  second argument a variable name, `valuation` returns the valuation of
the polynomial in that variable.

```julia-repl
julia> valuation(a,:y)
0

julia> valuation(a,:x)
1
```

"""
Pols.valuation(m::Mvp)=iszero(m) ? 0 : minimum(degree.(keys(m.d)))
Pols.valuation(m::Mvp,v::Symbol)=iszero(m) ? 0 : minimum(degree.(keys(m.d),v))

"""
`coefficient(p::Mvp,m::Monomial)`

The coefficient of the polynomial `p` on the monomial `m`.

```julia-repl
julia> @Mvp x,y; p=(x-y)^3
Mvp{Int64}: x³-3x²y+3xy²-y³

julia> coefficient(p,Monomial(:x,:x,:y)) # coefficient on x²y
-3

julia> coefficient(p,Monomial()) # constant coefficient
0
```
"""
coefficient(p::Mvp,m::Monomial)=p.d[m]

"""
  `coefficients(p::Mvp, var::Symbol)` 

returns  a Dict with keys the degree  in `var` and values the corresponding
coefficient of `p` with respect to `var`.

```julia-repl
julia> p=(x+y+inv(y))^4
Mvp{Int64}: x⁴+4x³y+4x³y⁻¹+6x²y²+12x²+6x²y⁻²+4xy³+12xy+12xy⁻¹+4xy⁻³+y⁴+4y²+6+4y⁻²+y⁻⁴

julia> coefficients(p,:x)
Dict{Int64, Mvp{Int64, Int64}} with 5 entries:
  0 => y⁴+4y²+6+4y⁻²+y⁻⁴
  4 => 1
  2 => 6y²+12+6y⁻²
  3 => 4y+4y⁻¹
  1 => 4y³+12y+12y⁻¹+4y⁻³

julia> coefficients(p,:y)
Dict{Int64, Mvp{Int64, Int64}} with 9 entries:
  0  => x⁴+12x²+6
  4  => 1
  -1 => 4x³+12x
  2  => 6x²+4
  -3 => 4x
  -2 => 6x²+4
  -4 => 1
  3  => 4x
  1  => 4x³+12x
```

The  same  caveat  is  applicable  to  `coefficients` as to evaluating: the
values  are always `Mvp`s. To get a list of scalars for the coefficients of
a  univariate polynomial represented  as a `Mvp`,  one should use `scal` on
the values of `coefficients`.
"""
function coefficients(p::Mvp{T,N},v::Symbol)where {T,N}
  if iszero(p) return Dict{Int,typeof(p)}() end
  d=Dict{N,typeof(p.d.d)}()
  for (m,c) in p.d
    found=false
    for (i,(v1,deg)) in enumerate(m.d)
      if v1==v 
        found=true
        d[deg]=push!(get!(d,deg,empty(p.d.d)),
                        Monomial(deleteat!(collect(m.d),i)...)=>c)
      end
    end
    if !found  d[0]=push!(get!(d,0,empty(p.d.d)),m=>c) end
  end
  Dict(dg=>Mvp(c...;check=false) for (dg,c) in d) # c... is sorted by defn of monomial order
end

"""
  `coefficient(p::Mvp, var::Symbol, d)` 

returns  the coefficient of degree `d` in the variable `var` in the `Mvp` `p`.

```julia-repl
julia> @Mvp x,y; p=(x+y^(1//2)+1)^3
Mvp{Int64,Rational{Int64}}: x³+3x²y½+3x²+3xy+6xy½+3x+y³⁄₂+3y+3y½+1

julia> coefficient(p,:y,1//2)
Mvp{Int64,Rational{Int64}}: 3x²+6x+3

julia> coefficient(p,:x,1)
Mvp{Int64,Rational{Int64}}: 3y+6y½+3
```
"""
function coefficient(p::Mvp,var::Symbol,d)
  res=empty(p.d.d)
  for (m,c) in p.d
    found=false
    for (i,(v1,deg)) in enumerate(m.d)
      if v1==var && deg==d
        found=true
        push!(res,Monomial(ModuleElt(deleteat!(copy(m.d.d),i);check=false))=>c)
        break
      elseif v1>var break
      end
    end
    if d==0 && !found push!(res,m=>c) end
  end
  Mvp(ModuleElt(res;check=false))
end

"""
  `Pol(p::Mvp)`

  converts the one-variable `Mvp` `p` to a polynomial. It is an error if `p`
  has more than one variable.

```julia-repl
julia> @Mvp x; Pol(:q); Pol(x^2+x)
Pol{Int64}: q²+q
```
"""
function Pol(x::Mvp{T})where T
  l=variables(x)
  if isempty(l) return Pol(scal(x)) end
  if length(l)>1 error("cannot convert $(length(l))-variate Mvp to Pol") end
  val=valuation(x)
  p=zeros(T,degree(x)-val+1)
  for (mon,coeff) in x.d p[degree(mon)-val+1]=coeff end
  Pol(p,val)
end

"""
  `Pol(p::Mvp,v::Symbol)`

returns  a polynomial whose coefficients are  the coefficients of the `Mvp`
`p`  with respect to the variable `v`.  The variable `v` should appear only
with integral powers in `p`.

```julia-repl
julia> p=(x+y^(1//2))^3
Mvp{Int64,Rational{Int64}}: x³+3x²y½+3xy+y³⁄₂

julia> Pol(:q); Pol(p,:x)
Pol{Mvp{Int64, Rational{Int64}}}: q³+3y½q²+3yq+y³⁄₂
```
"""
function Pols.Pol(p::Mvp{T,N},var::Symbol)where{T,N}
  d=degree(p,var)
  v=Int(valuation(p,var))
  res=[Pair{Monomial{N},T}[] for i in v:d]
  for (m,c) in p.d
    found=false
    for (i,(v1,deg)) in enumerate(m.d)
      if v1==var 
        found=true
        push!(res[Int(deg)-v+1],
             Monomial(ModuleElt(deleteat!(copy(m.d.d),i);check=false))=>c)
        break
      elseif v1>var break
      end
    end
    if !found push!(res[-v+1],m=>c) end
  end
  Pol(map(x->Mvp(ModuleElt(x;check=false)),res),v)
end

"""
`variables(p::Mvp)`
`variables(v::Array{Mvp})`

returns the list of variables of all `p` as a sorted list of `Symbol`s.

```julia-repl
julia> @Mvp x,y,z

julia> variables([x+y+1,z])
3-element Vector{Symbol}:
 :x
 :y
 :z
```
"""
variables(pp::AbstractArray{<:Mvp})=sort(unique!([k1 for p in pp for (k,v) in p.d for (k1,v1) in k.d]))

variables(p::Mvp)=sort(unique!([k1 for (k,v) in p.d for (k1,v1) in k.d]))

"""
`scal(p::Mvp)`

If `p` is a scalar, return that scalar, otherwise return `nothing`.

```julia-repl
julia> p=Mvp(:x)+1
Mvp{Int64}: x+1

julia> w=p(x=4)
Mvp{Int64}: 5

julia> scal(w)
5

julia> typeof(scal(w))
Int64
```
if  `p` is an array, then apply `scal` to its elements and return `nothing`
if it contains any `Mvp` which is not a scalar.
"""
function scal(p::Mvp{T})where T
  if iszero(p) return zero(T) end
  if length(p.d)==1 
    (m,c)=first(p.d)
    if isone(m) return c end
  end
end

function scal(m::AbstractArray{<:Mvp})
  p=scal.(m)
  if !any(isnothing,p) return p end
end

"""
    `value(p::Mvp,:x₁=>v₁,:x₂=>v₂,...)`
     ̀(p::Mvp)(x₁=v₁,…,xₙ=vₙ)`

returns  the value of  `p` when doing  the simultaneous substitution of the
variable `:x1` by `v1`, of `x2` by `v2`, …

```julia-repl
julia> p=-2+7x^5*inv(y)
Mvp{Int64}: 7x⁵y⁻¹-2

julia> p(x=2)
Mvp{Int64}: -2+224y⁻¹

julia> p(y=1)
Mvp{Int64}: 7x⁵-2

julia> p(x=2,y=1)
Mvp{Int64}: 222
```

One should pay attention to the fact that the last value is not an integer,
but  a constant `Mvp` (for consistency). See the function `scal` for how to
convert such constants to their base ring.

```julia-repl
julia> p(x=y)
Mvp{Int64}: 7y⁴-2

julia> p(x=y,y=x)
Mvp{Int64}: -2+7x⁻¹y⁵
```
Evaluating an `Mvp` which is a Puiseux polynomial may cause calls to `root`

```julia-repl
julia> p=x^(1//2)*y^(1//3)
Mvp{Int64,Rational{Int64}}: x½y⅓

julia> p(;x=y)
Mvp{Int64,Rational{Int64}}: y⁵⁄₆


julia> p(;x=4)
Mvp{Int64,Rational{Int64}}: 2y⅓

julia> p(;y=2.0)
Mvp{Float64,Rational{Int64}}: 1.2599210498948732x½
```
"""
function value(p::Mvp,k::Pair...)
  res1=badi=nothing
  if length(k)==1 (key,val)=k[1] 
  # special code when length(k)==1 and !(key∈ variables(p)) is 10x faster
  elseif isempty(k) return p
  else vv=Dict(k) 
  end 
  for (i,(m,c1)) in enumerate(p.d)
    res=badj=nothing
    for (j,(v,c)) in enumerate(m.d)
      if length(k)==1 && key!=v continue end
      if length(k)>1 
        if !haskey(vv,v) continue 
        else val=vv[v] 
        end
      end
      if badj==nothing 
        badj=Int[]
        res=Mvp(c1) 
      end
      if isinteger(c) res*=val^Int(c)
      else res*=root(val,denominator(c))^numerator(c)
      end
      push!(badj,j)
    end
    if badj!==nothing
      res*=Monomial(deleteat!(collect(m.d),badj)...)
      badj=nothing
      if badi===nothing res1=res
        badi=Int[]
      else res1+=res
      end
      push!(badi,i)
    end
 #  println("badi=$badi m=$m c=$c res1=$res1")
  end
  if badi!==nothing Mvp(deleteat!(collect(p.d),badi)...;check=false)+res1
  else p
  end
end

(p::Mvp)(;arg...)=value(p,arg...)

function root(p::Mvp,n::Real=2)
  if iszero(p) return p end
  n=Int(n)
  if length(p.d)>1 
   throw(DomainError(p,"$(ordinal(n)) root of non-monomial not implemented")) 
  end
  (m,c)=first(p.d)
  Mvp(root(m,n)=>root(c,n))
end

"""
`Base.:^(p,m;vars=variables(p))`

Implements  the action of  a matrix on  `Mvp`s. `vars` should  be a list of
symbols   representing  variables.   The  polynomial   `p`  is  changed  by
simultaneous  substitution in it of  `vᵢ` by `(v×m)ᵢ` where  `v` is the row
vector  of  the  `Mvp(vᵢ)`.  If  `vars`  is  omitted,  it  is  taken  to be
`variables(p)`.

```julia-repl
julia> @Mvp x,y

julia> (x+y)^[1 2;3 1]
Mvp{Int64}: 3x+4y
```
"""
Base.:^(p::Mvp,m::AbstractMatrix;vars=variables(p))=p(;map(Pair,vars,permutedims(Mvp.(vars))*m)...)

Pols.positive_part(p::Mvp)=
  Mvp(ModuleElt(m=>c for (m,c) in p.d if all(x->last(x)>0,m.d.d);check=false))

Pols.negative_part(p::Mvp)=
  Mvp(ModuleElt(m=>c for (m,c) in p.d if all(x->last(x)<0,m.d.d);check=false))

Pols.bar(p::Mvp)=Mvp(ModuleElt(inv(m)=>c for (m,c) in p.d))

"""
The  function 'Derivative(p,v)' returns the  derivative of 'p' with respect
to  the variable given by the string 'v'; if 'v' is not given, with respect
to the first variable in alphabetical order.

```julia-repl
julia> @Mvp x,y;p=7x^5*y^-1-2
Mvp{Int64}: 7x⁵y⁻¹-2

julia> derivative(p,:x)
Mvp{Int64}: 35x⁴y⁻¹

julia> derivative(p,:y)
Mvp{Int64}: -7x⁵y⁻²

julia> derivative(p)
Mvp{Int64}: 35x⁴y⁻¹

julia> p=x^(1//2)*y^(1//3)
Mvp{Int64,Rational{Int64}}: x½y⅓

julia> derivative(p,:x)
Mvp{Rational{Int64},Rational{Int64}}: (1//2)x⁻½y⅓

julia> derivative(p,:y)
Mvp{Rational{Int64},Rational{Int64}}: (1//3)x½y⁻²⁄₃

julia> derivative(p,:z)
Mvp{Rational{Int64},Rational{Int64}}: 0
```
"""
function Pols.derivative(p::Mvp,v=first(variables(p)))
  # check needed because 0 could appear in coeffs
  Mvp(ModuleElt(m*Monomial(v)^-1=>c*degree(m,v) for (m,c) in p.d))
end

function Util.exactdiv(p::Mvp,q::Mvp)
  if iszero(q) error("cannot divide by 0")
  elseif length(q.d)==1 
    m,c=first(q.d)
    return Mvp(ModuleElt(inv(m)*m1=>exactdiv(c1,c) for (m1,c1) in p.d;check=false))
  elseif length(p.d)<2 return iszero(p) ? p : nothing
  end 
  var=first(first(p.d)[1].d)[1]
  res=zero(p)
  mq=degree(q,var)
  cq=coefficient(q,var,mq)
  while length(p.d)!=0
#   if length(p.d)<length(q.d) return nothing end
    mp=degree(p,var)
    t=exactdiv(coefficient(p,var,mp),cq)
    if t===nothing return nothing end
    if mp!=mq t=Monomial(var=>mp-mq)*t end
    res+=t
    p-=t*q
  end
  res
end

Util.exactdiv(p::Mvp,q::Number)=exactdiv(p,Mvp(q))
 
"""
Returns  the Gcd  of the  'Mvp' arguments.  The arguments  must be  true
polynomials.

```julia-repl
julia> gcd(x^2-y^2,(x+y)^2)
Mvp{Int64}: x+y
```
"""
function Base.gcd(a::Mvp,b::Mvp)
  vars=variables([a,b])
  if isempty(vars) return Mvp(1) end
  v=first(vars)
  if length(vars)==1 gcd(Pol(a),Pol(b))(Mvp(v))
  else               srgcd(Pol(a,v),Pol(b,v))(Mvp(v))
  end
end

Base.gcd(v::AbstractArray{<:Mvp})=reduce(gcd,v;init=Mvp(0))

"""
`laurent_denominator(p1,p2,…)`

returns  the unique true monomial  `m` of minimal degree  such that for all
the Laurent polynomials `p1,p2,…` the product `m*pᵢ` is a true polynomial.

```julia-repl
julia> laurent_denominator(x^-1,y^-2+x^4)
Monomial{Int64}:xy²
```
"""
laurent_denominator(pp::Mvp...)=inv(gcd([m for p in pp for (m,c) in p.d]))

"""
`lcm(p1,p2,...)`

Returns  the Lcm  of the  `Mvp` arguments.  The arguments  must be  true
polynomials.

```julia-repl
julia> lcm(x^2-y^2,(x+y)^2)
Mvp{Int64}: x³+x²y-xy²-y³
```
"""
Base.lcm(a::Mvp,b::Mvp)=exactdiv(a*b,gcd(a,b))

Base.eltype(p::Mvp{T,N}) where{T,N} =T

Base.denominator(p::Mvp)=lcm(denominator.(values(p.d)))
Base.numerator(p::Mvp{<:Rational{T},N}) where{T,N} =convert(Mvp{T,N},p*denominator(p))
Base.numerator(p::Mvp{<:Cyc{<:Rational{<:T}},N}) where{T,N} =convert(Mvp{Cyc{T},N},p*denominator(p))
#----------------------------- Mvrf -----------------------------------
struct Mvrf{T}
  num::Mvp{T,Int}
  den::Mvp{T,Int}
  function Mvrf(a::Mvp{T,Int},b::Mvp{T,Int};check=true,check1=true)where{T}
    if iszero(b) error("division by 0") end
    if a==0 return new{T}(a,Mvp(T(1))) end
    if check
      d=laurent_denominator(a,b)
      a*=d
      b*=d
      if check1
        d=gcd(a,b)
        a=exactdiv(a,d)
        b=exactdiv(b,d)
        if T<:Rational{<:Integer} || T<:Cyc{<:Rational{<:Integer}}
          d=lcm(denominator(a),denominator(b))
          a=numerator(a*d)
          b=numerator(b*d)
        end
        T1=promote_type(eltype(a),eltype(b))
        return new{T1}(a,b)
      end
    end
    new{T}(a,b)
  end
end

Base.numerator(p::Mvrf)=p.num
Base.denominator(p::Mvrf)=p.den

function Mvrf(a::Mvp{T1,Int},b::Mvp{T2,Int};check=true,check1=true)where{T1,T2}
  T=promote_type(T1,T2)
  Mvrf(convert(Mvp{T,Int},a),convert(Mvp{T,Int},b);check,check1)
end

function Base.convert(::Type{Mvrf{T}},p::Mvrf{T1}) where {T,T1}
  Mvrf(convert(Mvp{T,Int},p.num),convert(Mvp{T,Int},p.den);check=false)
end

function Base.convert(::Type{Mvp{T,Int}},p::Mvrf{T1}) where {T,T1}
  if isone(length(p.den.d)) 
    (m,c)=first(p.den.d)
    return convert(Mvp{T,Int},p.num*inv(m)*inv(c)) 
  end
  error("cannot convert ",p," to Mvp{",T,"}")
end

function Base.convert(::Type{Mvrf{T}},p::Mvp{T1,N}) where {T,T1,N}
  Mvrf(convert(Mvp{T,Int},p),Mvp(T(1));check1=false)
end

function Base.convert(::Type{Mvrf{T}},p::Number) where {T}
  Mvrf(convert(Mvp{T,Int},p),Mvp(T(1));check=false)
end

function Base.promote_rule(a::Type{Mvp{T1,Int}},b::Type{Mvrf{T2}})where {T1,T2}
  Mvrf{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Mvrf{T1}},b::Type{Mvrf{T2}})where {T1,T2}
  Mvrf{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{T1},b::Type{Mvrf{T2}})where {T1<:Integer,T2}
  Mvrf{promote_type(T1,T2)}
end

(::Type{Mvrf{T}})(a) where T=convert(Mvrf{T},a)

Mvp(a::Mvrf{T}) where T =convert(Mvp{T,Int},a)

Base.broadcastable(p::Mvrf)=Ref(p)
Mvrf(a::Number)=Mvrf(Mvp(a),Mvp(1);check=false)
Mvrf(a::Mvp)=Mvrf(a,Mvp(1);check1=false)
Base.copy(a::Mvrf)=Mvrf(a.num,a.den;check=false)
Base.one(a::Mvrf)=Mvrf(one(a.num),one(a.den);check=false)
Base.one(::Type{Mvrf{T}}) where T =Mvrf(one(Mvp{T,Int}),one(Mvp{T,Int});check=false)
Base.one(::Type{Mvrf})=Mvrf(one(Mvp{Int,Int}),one(Mvp{Int,Int});check=false)
Base.zero(::Type{Mvrf{T}}) where T =Mvrf(zero(Mvp{T,Int}),one(Mvp{T,Int});check=false)
Base.zero(::Type{Mvrf})=Mvrf(zero(Mvp{Int,Int}),one(Mvp{Int,Int});check=false)
# next 3 stuff to make inv using LU work (abs is stupid)
Base.abs(p::Mvrf)=p
Base.conj(p::Mvrf)=Mvrf(conj(p.num),conj(p.den);check=false)
Base.adjoint(a::Mvrf)=conj(a)
Base.cmp(a::Mvrf,b::Mvrf)=cmp([a.num,a.den],[b.num,b.den])
Base.isless(a::Mvrf,b::Mvrf)=cmp(a,b)==-1

function Base.show(io::IO, ::MIME"text/plain", a::Mvrf)
  if !haskey(io,:typeinfo) print(io,typeof(a),": ") end
  show(io,a)
end

function Base.show(io::IO,a::Mvrf)
  if a.den==-1 a=Mvrf(-a.num,-a.den;check=false) end
  n=sprint(show,a.num; context=io)
  if  get(io, :limit,true) && a.den==one(a.den)
    print(io,n)
  else
    print(io,Util.bracket_if_needed(n))
    n=sprint(show,a.den; context=io)
    print(io,"/",Util.bracket_if_needed(n))
  end
end

Base.inv(a::Mvrf)=Mvrf(a.den,a.num;check=false)

function Base.://(a::Mvp,b::Mvp)
  if iszero(a) return a end
  if length(b.d)==1
    (m,c)=first(b.d)
    return c^2==1 ? Mvp([m1/m=>c1*c for (m1,c1) in a.d]...) : 
                    Mvp([m1/m=>c1//c for (m1,c1) in a.d]...)
  end
  Mvrf(a,b)
end

Base.:/(p::Mvp,q::Mvp)=p//q

function Base.inv(p::Mvp)
  if length(p.d)==1 return Mvp(1)//p end
  Mvrf(Mvp(1),p;check1=false)
end

Base.://(a::Mvrf,b::Mvrf)=a*inv(b)
Base.:/(a::Mvrf,b::Mvrf)=a*inv(b)
Base.:/(a::Mvrf,b::Union{Number,Mvp})=Mvrf(a.num,a.den*b)
Base.://(a::Mvrf,b::Union{Number,Mvp})=Mvrf(a.num,a.den*b)
Base.:/(a::Union{Number,Mvp},b::Mvrf)=a*inv(b)
Base.://(p::Number,q::Mvp)=Mvp(p)//q
Base.:/(p::Number,q::Mvp)=p//q
Base.:/(p::Mvp,q)=Mvp(p.d/q,p.v)
Base.://(p::Mvp,q)=iszero(p) ? p : Mvp(p.d//q,p.v)

Base.:*(a::Mvrf,b::Mvrf)=Mvrf(a.num*b.num,a.den*b.den)

Base.:*(a::Mvrf,b::Mvp)=Mvrf(a.num*b,a.den)
Base.:*(b::Mvp,a::Mvrf)=Mvrf(a.num*b,a.den)
Base.:*(a::Mvrf,b::T) where T =Mvrf(a.num*b,a.den;check=false)
Base.:*(b::T,a::Mvrf) where T =a*b

Base.:^(a::Mvrf, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) : 
                              Base.power_by_squaring(inv(a),-n)
Base.:+(a::Mvrf,b::Mvrf)=Mvrf(a.num*b.den+a.den*b.num,a.den*b.den)
Base.:+(a::Mvrf,b::Union{Number,Mvp})=a+Mvrf(b)
Base.:+(b::Union{Number,Mvp},a::Mvrf)=a+Mvrf(b)
Base.:-(a::Mvrf)=Mvrf(-a.num,a.den;check=false)
Base.:-(a::Mvrf,b::Mvrf)=a+(-b)
Base.:-(a::Mvrf,b)=a-Mvrf(b)
Base.:-(b,a::Mvrf)=Mvrf(b)-a
Base.zero(a::Mvrf)=Mvrf(zero(a.num),one(a.num);check=false)

value(p::Mvrf,k::Pair...)=value(p.num,k...)/value(p.den,k...)
(p::Mvrf)(;arg...)=value(p,arg...)

#julia1.6.0-rc2> @btime Mvps.fateman(15)
# 4.040 s (15219390 allocations: 5.10 GiB)
function fateman(n)
  f=(1+Mvp(:x)+Mvp(:y)+Mvp(:z)+Mvp(:t))*one(n)
  f=f^n
  length((f*(f+1)).d)
end

end
