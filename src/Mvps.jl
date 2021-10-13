"""
The  type  `Mvp`  ("multivariate  polynomials")  in  this  file  implements
"Puiseux polynomials", that is linear combinations of monomials of the type
`x₁^{a₁}…  xₙ^{aₙ}` where `xᵢ`  are variables and  `aᵢ` are exponents which
can   be  arbitrary  rational  numbers  (we  use  Puiseux  polynomial  with
cyclotomic  coefficients as splitting fields of cyclotomic Hecke algebras).
Some  functions described below work  only with polynomials where variables
are  raised to integral powers;  we will refer to  such objects as "Laurent
polynomials"; some functions require further that variables are raised only
to positive powers: we refer then to "true polynomials".

This  file implements  also multivariate  rational fractions  (type `Mvrf`)
where the numerator and denominator are true polynomials.

`@Mvp x₁,…,xₙ`

assigns  to each  Julia name  `xᵢ` an  `Mvp` representing  an indeterminate
suitable   to  build   multivariate  polynomials   or  rational  fractions.
`Mvp(:x₁)` creates the same `Mvp` without assigning it to a Julia variable.

```julia-repl
julia> @Mvp x,y

julia> (x+y^-1)^3
Mvp{Int64}: x³+3x²y⁻¹+3xy⁻²+y⁻³

julia> x+Mvp(:z)
Mvp{Int64}: x+z
```

`Mvp(x::Number)`  returns the multivariate polynomial  with only a constant
term, equal to `x`.

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

One  can evaluate an  `Mvp` or a  `Ratfrac` when setting  the value of some
variables by using the function call syntax:

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

Note  that  an  `Mvp`  always  evaluates  to an `Mvp`, for consistency. The
function  `scalar` converts a constant `Mvp`  to that constant (and returns
`nothing` if the argument is not constant:

```julia-repl
julia> p(x=1,y=2)
Mvp{Int64}: 3

julia> scalar(p(x=1,y=2))
3

julia> v=(x^(1//2))(x=2.0)
Mvp{Float64}: 1.4142135623730951

julia> scalar(v)
1.4142135623730951
```

see the functions `coefficient`, `coefficients, `Pol` to take apart `Mvp`s.
`Mvrf` are dissected using `numerator` and `denominator`.

Despite  the degree of generality of our  polynomials, the speed is not too
shabby. For the Fateman test f(f+1) where f=(1+x+y+z+t)^15, we take 4sec.
According to the Nemo paper, Sagemath takes 10sec and Nemo takes 1.6sec.
"""
module Mvps
# benchmar: (x+y+z)^3     2.3μs 48 alloc
using ModuleElts: ModuleElt, ModuleElts
using ..Util: ordinal, printTeX, stringexp
using ..Pols: Pols, Pol, srgcd, positive_part, negative_part, bar, derivative,
              valuation, degree, scalar, exactdiv

#import Gapjm: coefficients, valuation
import ..Cycs: root
using ..Cycs: Cyc
# to use as a stand-alone module comment above line, uncomment next, and
# define root for the coefficients you want (at least root(1,n)=1)
#export root
export coefficients, coefficient
export Mvp, Monomial, @Mvp, variables, value, laurent_denominator, Mvrf
#------------------ Monomials ---------------------------------------------
struct Monomial{T}
  d::ModuleElt{Symbol,T}   
end

Monomial(a::Union{Pair,Symbol}...)=Monomial(ModuleElt(map(x->x isa Symbol ?
          x=>1 : x,a)...)) # to be avoided when preformance needed
Monomial()=one(Monomial{Int})

Base.convert(::Type{Monomial{T}},v::Symbol) where T=Monomial(v=>T(1))
Base.convert(::Type{Monomial{T}},m::Monomial{N}) where {T,N}= 
  Monomial(convert(ModuleElt{Symbol,T},m.d))

# we need a promote rule to mix fractional momonials with others in an Mvp
Base.promote_rule(::Type{Monomial{N}},::Type{Monomial{M}}) where {N,M}=
  Monomial{promote_type(N,M)}

Base.:*(a::Monomial, b::Monomial)=Monomial(a.d+b.d)
Base.isone(a::Monomial)=iszero(a.d)
#Base.iszero(a::Monomial)=false
Base.one(::Type{Monomial{T}}) where T=Monomial(zero(ModuleElt{Symbol,T}))
Base.one(m::Monomial)=Monomial(zero(m.d))
Base.inv(a::Monomial)=Monomial(-a.d)
Pols.exactdiv(a::Monomial, b::Monomial)=a*inv(b)
Base.:/(a::Monomial, b::Monomial)=a*inv(b)
Base.://(a::Monomial, b::Monomial)=a*inv(b)
Base.:^(x::Monomial,p)=Monomial(x.d*p)
#Base.getindex(a::Monomial,k)=getindex(a.d,k)

function Base.show(io::IO,m::Monomial)
  replorTeX=get(io,:TeX,false) || get(io,:limit,false)
  if !replorTeX
    print(io,"Monomial(")
    join(io,map(m.d)do (s,c)
      if c==1 repr(s)
      elseif c==2 join([repr(s),repr(s)],",")
      else repr(s=>c)
      end
    end,",")
    print(io,")")
    return
  end
  if isone(m) return end
  start=true
  for (v,d) in m.d
    if !(start || replorTeX) print(io,"*") end
    if replorTeX printTeX(io,string(v)) else print(io,string(v)) end
    if !isone(d) 
      if isone(denominator(d)) d=numerator(d) end
      if replorTeX print(io,stringexp(io,d))
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
`lexless(a::Monomial, b::Monomial)`
The  "lex" ordering,  where `a<b`  if the  first variable  in `a/b`
occurs to a positive power.
"""
function lexless(a::Monomial, b::Monomial)
  for ((va,pa),(vb,pb)) in zip(a.d,b.d)
    if va!=vb return va<vb ? pa>0 : pb<0 end
    if pa!=pb return isless(pb,pa) end
  end
  la=length(a.d)
  lb=length(b.d)
  @inbounds la>lb ? last(a.d.d[lb+1])>0 : lb>la ? last(b.d.d[la+1])<0 : false
end

"""
`deglexless(a::Monomial, b::Monomial)`
The "deglex" ordering, where `a<b̀` if `degree(a)<degree(b)` or the degrees
are equal but `lexless(a,b)`.
"""
function deglexless(a::Monomial, b::Monomial)
  da=degree(a);db=degree(b)
  if da!=db return isless(da,db) end
  lexless(a,b)
end

"""
`isless(a::Monomial,b::Monomial)`

For  our implementation of `Mvp`s to  work, `isless` must define a monomial
order (that is, for monomials `m,a,b` we have `a<b => a*m<b*m`). By default
we  use the  "lex" ordering.
"""
@inline Base.isless(a::Monomial, b::Monomial)=lexless(a,b)

Base.:(==)(a::Monomial, b::Monomial)=a.d==b.d

#monomial d of greatest degree in each variable such that (m/d,n/d) positive
Base.gcd(m::Monomial,n::Monomial)=Monomial(ModuleElts.merge2(min,m.d,n.d))
Base.gcd(v::AbstractArray{<:Monomial})=reduce(gcd,v)

#monomial l of smallest degree in each variable such that (l/m,l/n) positive
Base.lcm(m::Monomial,n::Monomial)=Monomial(ModuleElts.merge2(max,m.d,n.d))
Base.lcm(v::AbstractArray{<:Monomial})=reduce(lcm,v)

Base.hash(a::Monomial, h::UInt)=hash(a.d,h)

Pols.degree(m::Monomial)=sum(values(m.d);init=0)
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

function Base.show(io::IO, x::Mvp)
  if get(io,:limit,false) || get(io,:TeX,false)
    show(IOContext(io,:showbasis=>nothing),x.d)
    # :showbasis=>nothing necessary if called when already showing a ModuleElt
  else
    print(io,"Mvp(")
    join(io,repr.(x.d),",")
    print(io,")")
  end
end

monom(x::Mvp)=length(x.d)==1 # x is a non-zero monomial
Base.zero(p::Mvp)=Mvp(zero(p.d))
Base.zero(::Type{Mvp{T,N}}) where {T,N}=Mvp(zero(ModuleElt{Monomial{N},T}))
Base.one(::Type{Mvp{T,N}}) where {T,N}=Mvp(one(Monomial{N})=>one(T))
Base.one(::Type{Mvp{T}}) where T=one(Mvp{T,Int})
Base.one(::Type{Mvp})=Mvp(1)
Base.one(p::Mvp{T,N}) where {T,N}=one(Mvp{T,N})
Base.isone(x::Mvp)=monom(x) && isone(first(first(x.d)))&&isone(last(first(x.d)))
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
  if !monom(a) || !isone(first(first(a.d)))
      throw(InexactError(:convert,T,a)) 
  end
  convert(T,last(first(a.d)))
end
(::Type{T})(a::Mvp) where T<: Number=convert(T,a)

Base.isinteger(p::Mvp)=iszero(p) || (monom(p) &&
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
Base.:-(a::Mvp, b::Mvp)=Mvp(a.d-b.d)
Base.:-(a::Mvp, b::Number)=a-Mvp(b)
Base.:-(b::Number, a::Mvp)=Mvp(b)-a

Base.:*(a::Number, b::Mvp)=Mvp(b.d*a)
Base.:*(b::Mvp, a::Number)=a*b
# we have a monomial order so there is no ordering check in next line
Base.:*(a::Monomial, b::Mvp)=Mvp(ModuleElt(m*a=>c for (m,c) in b.d;check=false))
Base.:*(b::Mvp,a::Monomial)=a*b
function Base.:*(a::Mvp, b::Mvp)
  if length(a.d)>length(b.d) a,b=(b,a) end
  if iszero(a) return a 
  elseif iszero(b) return b
  elseif isone(a) return b
  elseif isone(b) return a
  end
  let b=b # needed !!!!
    sum(b*m*c for (m,c) in a.d)
  end
end

Base.:/(p::Mvp,q::Number)=Mvp(p.d/q)
Base.://(p::Mvp,q::Number)=Mvp(p.d//q)
Base.div(a::Mvp,b::Number)=Mvp(merge(div,a.d,b))
Pols.exactdiv(a::Mvp,b::Number)=Mvp(merge(exactdiv,a.d,b;check=false))

"""
`conj(p::Mvp)` acts on the coefficients of `p`

```julia-repl
julia> @Mvp x;conj(im*x+1)
Mvp{Complex{Int64}}: (0 - 1im)x+1 + 0im
```
"""
Base.conj(a::Mvp)=Mvp(merge(conj,a.d;check=false))

function Base.:^(x::Mvp, p::Union{Integer,Rational})
  if isinteger(p) p=Int(p) end
  if iszero(p) return one(x)
  elseif iszero(x) return x
  elseif !isinteger(p) return root(x,denominator(p))^numerator(p) 
  elseif isone(p) return x 
  elseif monom(x)
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
Pols.degree(p::Mvp)=maximum(degree,keys(p.d);init=0)
Pols.degree(m::Mvp,v::Symbol)=maximum(x->degree(x,v),keys(m.d);init=0)

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
a  univariate polynomial represented as a `Mvp`, one should use `scalar` on
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
  `Pol(p::Mvp{T}) where T`

  converts the one-variable `Mvp{T}` `p` to a `Pol{T}`. It is an error if `p`
  has more than one variable.

```julia-repl
julia> @Mvp x; @Pol q; Pol(x^2+x)
Pol{Int64}: q²+q
```
"""
function Pol(x::Mvp{T})where T
  l=variables(x)
  if isempty(l) return Pol(scalar(x)) end
  if length(l)>1 error("cannot convert $(length(l))-variate Mvp to Pol") end
  val=Int(valuation(x))
  p=zeros(T,Int(degree(x))-val+1)
  for (mon,coeff) in x.d p[Int(degree(mon))-val+1]=coeff end
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
`Mvp(p)` converts  the `Pol`  `p` to  an  `Mvp`. 

```julia-repl
julia> @Pol q
Pol{Int64}: q

julia> Mvp(q^2+q)
Mvp{Int64}: q²+q
```
"""
Mvp(x::Pol)=convert(Mvp,x)

Base.convert(::Type{Mvp{T,N}},p::Pol) where{T,N}=
                     p(Mvp(convert(Monomial{N},Pols.varname[])=>one(T)))
Base.convert(::Type{Mvp},p::Pol)=p(Mvp(Pols.varname[]))

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
`scalar(p::Mvp)`

If `p` is a scalar, return that scalar, otherwise return `nothing`.

```julia-repl
julia> p=Mvp(:x)+1
Mvp{Int64}: x+1

julia> w=p(x=4)
Mvp{Int64}: 5

julia> scalar(w)
5

julia> typeof(scalar(w))
Int64
```
if  `p` is an array, then apply `scalar` to its elements and return `nothing`
if it contains any `Mvp` which is not a scalar.
"""
function Pols.scalar(p::Mvp{T})where T
  if iszero(p) return zero(T) end
  if monom(p)
    (m,c)=first(p.d)
    if isone(m) return c end
  end
end

function Pols.scalar(m::AbstractArray{<:Mvp})
  p=scalar.(m)
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
but  a constant `Mvp` (for consistency).  See the function `scalar` for how
to convert such constants to their base ring.

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
Mvp{Int64,Rational{Int64}}: y⅚

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
  if !monom(p)
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
Mvp{Rational{Int64},Rational{Int64}}: (1//3)x½y⁻⅔

julia> derivative(p,:z)
Mvp{Rational{Int64},Rational{Int64}}: 0
```
"""
function Pols.derivative(p::Mvp,v=first(variables(p)))
  # check needed because 0 could appear in coeffs
  Mvp(ModuleElt(m*Monomial(v=>-1)=>c*degree(m,v) for (m,c) in p.d))
end

# returns p/q when the division is exact, nothing otherwise
# Arguments must be true polynomials
function Pols.exactdiv(p::Mvp,q::Mvp)
  if iszero(q) error("cannot divide by 0")
  elseif iszero(p) || isone(q) return p
  elseif monom(q)
    m,c=first(q.d)
    return Mvp(ModuleElt(inv(m)*m1=>exactdiv(c1,c) for (m1,c1) in p.d;check=false))
  elseif monom(p) return nothing
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

Base.gcd(a::AbstractFloat,b::AbstractFloat)=one(a)*one(b)

"""
`gcd(p::Mvp,  q::Mvp)`  computes  the  `gcd`  of  the  'Mvp' arguments. The
arguments must be true polynomials.

```julia-repl
julia> gcd(x^2-y^2,(x+y)^2)
Mvp{Int64}: -x-y
```
"""
function Base.gcd(a::Mvp,b::Mvp)
  if iszero(a) return b
  elseif iszero(b) return a
  elseif monom(a)
    (m,c)=first(a.d)
    return Mvp(gcd(m,reduce(gcd,keys(b.d)))=>gcd(c,reduce(gcd,values(b.d))))
  elseif monom(b) return gcd(b,a)
  end
  va=variables(a)
  vb=variables(b)
  vars=intersect(va,vb)
  if isempty(vars) 
    if iszero(a) return b
    elseif iszero(b) return a
    else return Mvp(gcd(reduce(gcd,values(a.d)),reduce(gcd,values(b.d))))
    end
  end
  v=first(vars)
  if length(union(va,vb))==1 srgcd(Pol(a),Pol(b))(Mvp(v))
  else              srgcd(Pol(a,v),Pol(b,v))(Mvp(v))
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
Mvp{Int64}: -x³-x²y+xy²+y³
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
  # to speedup the constructor it can be told:
  # pol: we know a and b are true polynomials
  # prime: we know a and b are prime to each other
  function Mvrf(a::Mvp{T,Int},b::Mvp{T,Int};pol=false,prime=false)where{T}
    if iszero(b) error("division by 0") end
    if iszero(a) return new{T}(a,one(Mvp{T,Int})) end
    if !pol
      d=laurent_denominator(a,b)
      if !isone(d) a*=d;b*=d end
    end
    if !prime && !isone(b)
      d=gcd(a,b)
      if !isone(d) a=exactdiv(a,d);b=exactdiv(b,d) end
    end
    return new{T}(a,b)
  end
end

function Mvrf(a::Mvp{T1,Int},b::Mvp{T2,Int};pol=false,prime=false)where{T1,T2}
  T=promote_type(T1,T2)
  Mvrf(convert(Mvp{T,Int},a),convert(Mvp{T,Int},b);pol,prime)
end

function Mvrf(a::Mvp{T,Int},b::Mvp{T,Int})where T<:Rational
  Mvrf(numerator(a)*denominator(b),numerator(b)*denominator(a))
end

function Mvrf(a::Mvp{Cyc{T},Int},b::Mvp{Cyc{T},Int})where T<:Rational
  Mvrf(numerator(a)*denominator(b),numerator(b)*denominator(a))
end

Base.numerator(p::Mvrf)=p.num
Base.denominator(p::Mvrf)=p.den
Base.:(==)(a::Mvrf,b::Mvrf)=a.num==b.num && a.den==b.den

function Base.convert(::Type{Mvrf{T}},p::Mvrf{T1}) where {T,T1}
  Mvrf(convert(Mvp{T,Int},p.num),convert(Mvp{T,Int},p.den);pol=true,prime=true)
end

function Base.convert(::Type{Mvp{T,Int}},p::Mvrf{T1}) where {T,T1}
  if isone(p.den) return convert(Mvp{T,Int},p.num)
  elseif monom(p.den) 
    (m,c)=first(p.den.d)
    return convert(Mvp{T,Int},p.num*inv(m)*inv(c)) 
  end
  error("cannot convert ",p," to Mvp{",T,"}")
end

function Base.convert(::Type{Mvrf{T}},p::Mvp{T1,N}) where {T,T1,N}
  Mvrf(convert(Mvp{T,Int},p),Mvp(T(1));prime=true)
end

function Base.convert(::Type{Mvrf{T}},p::Number) where {T}
  Mvrf(convert(Mvp{T,Int},p),Mvp(T(1));pol=true)
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
Mvrf(a::Number)=Mvrf(Mvp(a),Mvp(1);pol=true)
Mvrf(a::Mvp)=Mvrf(a,Mvp(1);prime=true)
Base.copy(a::Mvrf)=Mvrf(a.num,a.den;pol=true,prime=true)
Base.one(a::Mvrf)=Mvrf(one(a.num),one(a.den);pol=true,prime=true)
Base.one(::Type{Mvrf{T}}) where T =Mvrf(one(Mvp{T,Int}),one(Mvp{T,Int});pol=true,prime=true)
Base.zero(::Type{Mvrf{T}}) where T =Mvrf(zero(Mvp{T,Int}),one(Mvp{T,Int});pol=true,prime=true)
Base.zero(a::Mvrf)=Mvrf(zero(a.num),one(a.num);pol=true,prime=true)
# next 3 stuff to make inv using LU work (abs is stupid)
Base.abs(p::Mvrf)=p
Base.conj(p::Mvrf)=Mvrf(conj(p.num),conj(p.den);pol=true,prime=true)
Base.adjoint(a::Mvrf)=conj(a)
Base.cmp(a::Mvrf,b::Mvrf)=cmp([a.num,a.den],[b.num,b.den])
Base.isless(a::Mvrf,b::Mvrf)=cmp(a,b)==-1

function Base.show(io::IO, ::MIME"text/plain", a::Mvrf)
  if !haskey(io,:typeinfo) print(io,typeof(a),": ") end
  show(io,a)
end

function Base.show(io::IO,a::Mvrf)
  if a.den==-1 a=Mvrf(-a.num,-a.den;pol=true) end
  n=sprint(show,a.num; context=io)
  if  get(io, :limit,true) && a.den==one(a.den)
    print(io,n)
  else
    print(io,Pols.bracket_if_needed(n))
    n=sprint(show,a.den; context=io)
    print(io,"/",Pols.bracket_if_needed(n))
  end
end

Base.inv(a::Mvrf)=Mvrf(a.den,a.num;pol=true,prime=true)

function Base.inv(p::Mvp)
  if monom(p)
    (m,c)=first(p.d)
    return Mvp(inv(m)=>c^2==1 ? c : 1/c) 
  end
  Mvrf(Mvp(1),p;pol=false,prime=true)
end

function Base.://(a::Mvp,b::Mvp)
  if iszero(a) return a end
  if monom(b) 
    (m,c)=first(b.d)
    return Mvp(ModuleElt(m1/m=>c^2==1 ? c1*c : c1//c for (m1,c1) in a.d))
  end
  Mvrf(a,b)
end

Base.:/(p::Mvp,q::Mvp)=p//q

Base.://(a::Mvrf,b::Mvrf)=a*inv(b)
Base.://(a::Mvrf,b::Number)=Mvrf(a.num,a.den*b;pol=true,prime=true)
Base.://(a::Mvrf,b::Mvp)=Mvrf(a.num,a.den*b)
Base.:/(a::Mvrf,b::Union{Mvrf,Number,Mvp})=a//b
Base.:/(a::Union{Number,Mvp},b::Mvrf)=a*inv(b)
Base.://(p::Number,q::Mvp)=Mvp(p)//q
Base.:/(p::Number,q::Mvp)=Mvp(p)/q

Base.:*(a::Mvrf,b::Mvrf)=Mvrf(a.num*b.num,a.den*b.den;pol=true)

Base.:*(a::Mvrf,b::Mvp)=Mvrf(a.num*b,a.den)
Base.:*(b::Mvp,a::Mvrf)=a*b
Base.:*(a::Mvrf,b::T) where T =Mvrf(a.num*b,a.den;pol=true)
Base.:*(b::T,a::Mvrf) where T =a*b

Base.:^(a::Mvrf, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) : 
                              Base.power_by_squaring(inv(a),-n)
Base.:+(a::Mvrf,b::Mvrf)=Mvrf(a.num*b.den+a.den*b.num,a.den*b.den;pol=true)
Base.:+(a::Mvrf,b::Union{Number,Mvp})=a+Mvrf(b)
Base.:+(b::Union{Number,Mvp},a::Mvrf)=a+Mvrf(b)
Base.:-(a::Mvrf)=Mvrf(-a.num,a.den;pol=true,prime=true)
Base.:-(a::Mvrf,b::Mvrf)=a+(-b)
Base.:-(a::Mvrf,b)=a-Mvrf(b)
Base.:-(b,a::Mvrf)=Mvrf(b)-a

value(p::Mvrf,k::Pair...)=value(p.num,k...)/value(p.den,k...)
(p::Mvrf)(;arg...)=value(p,arg...)

#julia1.6.3> @btime Mvps.fateman(15)
# 4.040 s (15219390 allocations: 5.10 GiB)
function fateman(n)
  f=(1+sum(Mvp.((:x,:y,:z,:t))))*one(n)
  f=f^n
  length((f*(f+1)).d)
end

end
