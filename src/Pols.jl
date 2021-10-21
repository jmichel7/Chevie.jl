"""
This  package,  which  depends  on  no other package, implements univariate
Laurent  polynomials (type `Pol{T}`), and  univariate rational fractions (type
`Frac{Pol{T}}`).

The  initial motivation was to have a simple way to port GAP polynomials to
Julia. The reasons for still having my own package are multiple:

  - I need  to have  a simple  and flexible  interface, which  I hope this
    provides.
  - I need my polynomials to behave  well when coefficients are in a ring,
    in which case I use pseudo-division and subresultant gcd.
  - For my use my polynomials are several times faster than those in the
    package `Polynomials`.
  - I need to  work with coefficients  of types `T`  where elements have a
    `zero`  method but `T` itself does not  have one (because `T` does not
    contain the necessary information. An example is `Mod{BigInt}`).
  - a final justification is that `Pols` is used by `Mvps`.

Here  Laurent  polynomials  have  the  parametric  type  `Pol{T}`. They are
constructed by giving a vector of coefficients of type `T`, and a valuation
(an `Int`). We call true polynomials those whose valuation is `≥0`.

There  is a  current variable  name (a  `Symbol`) used to print polynomials
nicely  at  the  repl  or  in  IJulia  or  Pluto. This name can be changed
globally,  or just changed for printing a given polynomial. But polynomials
do not record individually with which symbol they should be printed.

# Examples
```julia-repl
julia> Pol(:q) # define symbol used for printing and returns Pol([1],1)
Pol{Int64}: q

julia> @Pol q # same as q=Pol(:q)
Pol{Int64}: q

julia> Pol([1,2]) # valuation is taken to be 0 if omitted
Pol{Int64}: 2q+1

julia> 2q+1       # same polynomial
Pol{Int64}: 2q+1

julia> Pol()   # omitting all arguments gives Pol([1],1)
Pol{Int64}: q

julia> p=Pol([1,2,1],-1) # here the valuation is specified to be -1
Pol{Int64}: q+2+q⁻¹

julia> q+2+q^-1 # same polynomial
Pol{Int64}: q+2+q⁻¹
```

```julia-rep1
julia> print(p) # if not nice printing give an output which can be read back
Pol([1, 2, 1],-1)

# change the variable for printing just this time
julia> print(IOContext(stdout,:limit=>true,:varname=>"x"),p)
x+2+x⁻¹

julia> print(IOContext(stdout,:TeX=>true),p) # TeXable output (used in Pluto, IJulia)
q+2+q^{-1}
```

A  polynomial can be  taken apart with  the functions `valuation`, `degree`
and `getindex`. An index `p[i]` gives the coefficient of degree `i` of `p`.

```julia-repl
julia> valuation(p),degree(p)
(-1, 1)

julia> p[0], p[1], p[-1], p[10]
(2, 1, 1, 0)

julia> p[valuation(p):degree(p)]
3-element Vector{Int64}:
 1
 2
 1

julia> p[begin:end]  # the same as the above line
3-element Vector{Int64}:
 1
 2
 1

julia> coefficients(p)  # the same again
3-element Vector{Int64}:
 1
 2
 1
```

A  polynomial  is  a  *scalar*  if  the  valuation  and degree are `0`. The
function  `scalar` returns the constant coefficient  if the polynomial is a
scalar, and `nothing` otherwise.

```julia-repl
julia> Pol(1)
Pol{Int64}: 1

julia> convert(Pol{Int},1) # the same thing
Pol{Int64}: 1

julia> scalar(Pol(1))
1

julia> convert(Int,Pol(1)) # the same thing
1

julia> Int(Pol(1))         # the same thing
1

julia> scalar(q+1) # nothing
```

`Pol{T}` in arrays are promoted to the same type `T` (when the `T` involved
have a promotion) and a number is promoted to a polynomial. 

Usual  arithmetic (`+`, `-`,  `*`, `^`, `/`,  `//`, `one`, `isone`, `zero`,
`iszero`,  `==`) works. Elements  of type `<:Number`  or of type  `T` for a
`Pol{T}`   are  considered  as   scalars  for  scalar   operations  on  the
coefficients.

```julia-repl
julia> derivative(p)
Pol{Int64}: 1-q⁻²

julia> p=(q+1)^2
Pol{Int64}: q²+2q+1

julia> p=(q+1)^2
Pol{Int64}: q²+2q+1

julia> p/2
Pol{Float64}: 0.5q²+1.0q+0.5

julia> p//2
Pol{Rational{Int64}}: (1//2)q²+(1//1)q+1//2

julia> p(1//2) # value of p at 1//2
9//4

julia> p(0.5)
2.25

# interpolation: find p taking values [2.0,1.0,3.0] at [1,2,3]
julia> Pol([1,2,3],[2.0,1.0,3.0])  
Pol{Float64}: 1.5q²-5.5q+6.0
```

Polynomials  are scalars  for broadcasting.  They can  be sorted (they have
`cmp`   and  `isless`  functions  which   compare  the  valuation  and  the
coefficients), be keys in a `Dict` (they have a `hash` function).

The  functions  `divrem`,  `div`,  `%`,  `gcd`,  `gcdx`,  `lcm`, `powermod`
operate  between  true  polynomials  over  a  field,  using  the polynomial
division.  Over a  ring it  is better  to use  `pseudodiv` and  `srgcd`. By
default `gcd` between integer polynomials uses `srgcd`. `exactdiv` does the
division  (over a  field or  a ring)  when it  is exact,  otherwise returns
`nothing`.

```julia-repl
julia> divrem(q^3+1,2q+1) # changes coefficients to field elements
(0.5q²-0.25q+0.125, 0.875)

julia> divrem(q^3+1,2q+1//1) # case of field elements
((1//2)q²+(-1//4)q+1//8, 7//8)

julia> Pols.pseudodiv(q^3+1,2q+1) # pseudo-division keeps the ring
(4q²-2q+1, 7)

julia> (4q^2-2q+1)*(2q+1)+7 # but we get a multiple of the polynomial
Pol{Int64}: 8q³+8

julia> exactdiv(q+1,2) # nothing

julia> exactdiv(q+1,2.0)
Pol{Float64}: 0.5q+0.5
```

Finally,   `Pol`s  have   methods  `conj`,   `adjoint`  which   operate  on
coefficients,    a    `derivative`    methods,   methods   `positive_part`,
`negative_part`  and `bar` (useful for Kazhdan-Lusztig theory) and a method
`randpol` to produce random polynomials.

Rational   fractions  allow  to  invert  polynomials.  `Frac{Pol{T}}`s  are
normalized so that the numerator and denominator are true polynomials prime
to  each other. They  have the arithmetic  operations `+`, `-`  , `*`, `/`,
`//`,  `^`,  `inv`,  `one`,  `isone`,  `zero`,  `iszero` (which can operate
between a `Pol` or a `Number` and a `Frac{Pol{T}}`).


```julia-repl
julia> a=1/(q+1)
Frac{Pol{Int64}}: 1/(q+1)

julia> Pol(2/a) # convert back to `Pol`
Pol{Int64}: 2q+2

julia> numerator(a)
Pol{Int64}: 1

julia> denominator(a)
Pol{Int64}: q+1

julia> m=[q+1 q+2;q-2 q-3]
2×2 Matrix{Pol{Int64}}:
 q+1  q+2
 q-2  q-3

julia> n=inv(Frac.(m))
2×2 Matrix{Frac{Pol{Int64}}}:
 (-q+3)/(2q-1)  (-q-2)/(-2q+1)
 (q-2)/(2q-1)   (q+1)/(-2q+1)

julia> map(x->x(1),n) # evaluate at 1
2×2 Matrix{Float64}:
  2.0   3.0
 -1.0  -2.0

julia> map(x->x(1;Rational=true),n) # evaluate at 1 using //
2×2 Matrix{Rational{Int64}}:
  2//1   3//1
 -1//1  -2//1
```

Rational fractions are also scalars for broadcasting and can be sorted
(have `cmp` and `isless` methods).
"""
module Pols
export degree, valuation, Pol, derivative, shift, positive_part, negative_part,
       bar, derivative, srgcd, Frac, @Pol, scalar, coefficients, exactdiv, root

const varname=Ref(:x)

struct Pol{T}
  c::Vector{T}
  v::Int
  # Unexported inner constructor that bypasses all checks
  global Pol_(c::AbstractVector{T},v::Integer) where T=new{T}(c,v)
end

"""
  `Pol(c::AbstractVector,v::Integer=0;check=true,copy=true)`

  Make a polynomial of valuation `v` with coefficients `c`.

  Unless  `check` is `false`  normalize the result  by making sure that `c`
  has  no leading or  trailing zeroes (do  not set `check=false` unless you
  are sure this is already the case).

  Unless  `copy=false` the  contents of  `c` are  copied (you  can gain one
  allocation by setting `copy=false` if you know the contents can be shared)
"""
function Pol(c::AbstractVector{T},v::Integer=0;check=true,copy=true)where T
  if check # normalize c so there are no leading or trailing zeroes
    b=findfirst(!iszero,c)
    b===nothing && return Pol_(T[],0)
    e=findlast(!iszero,c)
    if b!=1 || e!=length(c) || copy return Pol_(view(c,b:e),v+b-1) end
  end
  Pol_(copy ? Base.copy(c) : c,v)
end

Pol()=Pol_([1],1)
Pol(a::T) where T<:Number=iszero(a) ? Pol_(T[],0) : Pol_([a],0)

"""
 `Pol(t::Symbol)`

 Sets the name of the variable for printing `Pol`s to `t`, and returns
 the polynomial of degree 1 equal to that variable.
"""
function Pol(t::Symbol)
  varname[]=t
  Pol()
end

"""
 `@Pol q`

 is equivalent to `q=Pol(:q)` excepted it creates `q` in the global scope of
 the current module, since it uses `eval`.
"""
macro Pol(t)
  if !(t isa Symbol) error("usage: @Pol <variable name>") end
  Base.eval(Main,:($t=Pol($(Core.QuoteNode(t)))))
end

Base.broadcastable(p::Pol)=Ref(p)
Base.copy(p::Pol)=Pol_(copy(p.c),p.v)

degree(a::Number)=0 # convenient
degree(p::Pol)=length(p.c)-1+p.v
valuation(p::Pol)=p.v
coefficients(p::Pol)=p.c
Base.lastindex(p::Pol)=degree(p)
Base.firstindex(p::Pol)=valuation(p)
@inbounds Base.getindex(p::Pol{T},i::Integer) where T=
i in firstindex(p):lastindex(p) ? p.c[i-p.v+1] : iszero(p) ? zero(T) : zero(p.c[1])

Base.getindex(p::Pol,i::AbstractVector{<:Integer})=getindex.(Ref(p),i)

Base.convert(::Type{Pol{T}},a::Number) where T=iszero(a) ? Pol_(T[],0) :
                                                           Pol_([T(a)],0)

(::Type{Pol{T}})(a) where T=convert(Pol{T},a)

# like convert for vectors does not always make a copy
Base.convert(::Type{Pol{T}},p::Pol{T1}) where {T,T1}=Pol_(Vector{T}(p.c),p.v)

function Base.promote_rule(a::Type{Pol{T1}},b::Type{Pol{T2}})where {T1,T2}
  Pol{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Pol{T1}},b::Type{T2})where {T1,T2<:Number}
  Pol{promote_type(T1,T2)}
end

Base.isinteger(p::Pol)=isinteger(scalar(p))

function scalar(p::Pol{T})where T
  if iszero(p) T(0)
  elseif iszero(valuation(p)) && iszero(degree(p)) p[0]
  end
end

function Base.convert(::Type{T},p::Pol) where T<:Number
  res=scalar(p)
  if res===nothing throw(InexactError(:convert,T,p)) end
  convert(T,res)
end

(::Type{T})(p::Pol) where T<:Number=convert(T,p)

Base.cmp(a::Pol,b::Pol)=cmp((a.c,a.v),(b.c,b.v))
Base.isless(a::Pol,b::Pol)=cmp(a,b)==-1
Base.hash(a::Pol, h::UInt)=hash(a.v,hash(a.c,h))

(p::Pol{T})(x) where T=iszero(p) ? zero(T) : evalpoly(x,p.c)*x^p.v

# efficient p↦ qˢ p
shift(p::Pol{T},s) where T=Pol_(p.c,p.v+s)

positive_part(p::Pol)=p.v>=0 ? copy(p) : Pol(view(p.c,1-p.v:length(p.c)),0)

negative_part(p::Pol)=degree(p)<=0 ? copy(p) : Pol(view(p.c,1:1-p.v),p.v)

# q↦ q⁻¹ on p
bar(p::Pol)=Pol_(reverse(p.c),-degree(p))

Base.:(==)(a::Pol, b::Pol)= a.c==b.c && a.v==b.v
Base.:(==)(a::Pol,b)= (iszero(a) && iszero(b))||(b!==nothing && scalar(a)==b)
Base.:(==)(b,a::Pol)= a==b

Base.one(a::Pol{T}) where T=Pol_([iszero(a) ? one(T) : one(a.c[1])],0)
Base.one(::Type{Pol{T}}) where T=Pol_([one(T)],0)
Base.one(::Type{Pol})=one(Pol{Int})
Base.isone(a::Pol)=!iszero(a) && scalar(a)==1
Base.zero(::Type{Pol{T}}) where T=Pol_(T[],0)
Base.zero(::Type{Pol})=zero(Pol{Int})
Base.zero(a::Pol{T}) where T=Pol_(T[],0)
Base.iszero(a::Pol)=isempty(a.c)
# next 3 stuff to make inv using LU work (abs is stupid)
Base.conj(p::Pol{T}) where T=Pol_(conj.(p.c),p.v)
Base.abs(p::Pol)=p
Base.adjoint(a::Pol)=conj(a)

ismonomial(p::Pol)=length(p.c)==1

function Base.show(io::IO, ::MIME"text/html", a::Pol)
  print(io, "\$")
  show(IOContext(io,:TeX=>true),a)
  print(io, "\$")
end

function Base.show(io::IO, ::MIME"text/plain", a::Pol)
  if !haskey(io,:typeinfo) print(io,typeof(a),": ") end
  show(io,a)
end

# 3 next methods copied from Util.jl in order to have a self-contained file.

# determines which coefficients should be bracketed for unambiguous display
function bracket_if_needed(c::String)
  if match(r"^[-+]?([^-+*/]|√-|{-)*(\(.*\))?$",c)!==nothing c
  else "("*c*")" 
  end
end

function format_coefficient(c::String)
  if c=="1" ""
  elseif c=="-1" "-"
  else bracket_if_needed(c)
  end
end

const supvec=['⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹']

function stringexp(io::IO,n::Integer)
  if isone(n) ""
  elseif get(io,:TeX,false) 
    n in 0:9 ? "^"*string(n) : "^{"*string(n)*"}"
  elseif get(io,:limit,false)
    res=Char[]
    if n<0 push!(res,'⁻'); n=-n end
    for i in reverse(digits(n)) push!(res,supvec[i+1]) end
    String(res)
  else "^"*string(n)
  end
end

function Base.show(io::IO,p::Pol)
  if !get(io,:limit,false) && !get(io,:TeX,false)
    if ismonomial(p) && isone(p.c[1]) && p.v==1 print(io,"Pol()")
    else print(io,"Pol(",p.c)
      if !iszero(p.v) print(io,",",p.v) end
      print(io,")")
    end
  elseif iszero(p) print(io,"0")
  else
    var=string(get(io,:varname,varname[]))
    for deg in degree(p):-1:valuation(p)
      c=p[deg]
      if iszero(c) continue end
      c=repr(c; context=IOContext(io,:typeinfo=>typeof(c)))
      if !iszero(deg)
        c=format_coefficient(c)*var*stringexp(io,deg)
      end
      if c[1]!='-' && deg!=degree(p) c="+"*c end
      print(io,c)
    end
  end
end

function Base.:*(a::Pol{T1},b::Pol{T2})where {T1,T2}
  T=promote_type(T1,T2)
  if iszero(a) || iszero(b) return zero(Pol{T}) end
  # below not zero(T) for types T like Mod{T1} which have no zero method
  res=fill(zero(a.c[1]*b.c[1]),length(a.c)+length(b.c)-1)
  for i in eachindex(a.c), j in eachindex(b.c)
@inbounds res[i+j-1]+=a.c[i]*b.c[j]
  end
  Pol_(res,a.v+b.v)
end

Base.:*(a::Pol, b::Number)=iszero(b) ? zero(a) : Pol_(a.c.*b,a.v)
Base.:*(a::Pol{T}, b::T) where T=Pol(a.c.*b,a.v;copy=false)
Base.:*(b::Number, a::Pol)=a*b
Base.:*(b::T, a::Pol{T}) where T=a*b

Base.:^(a::Pol, n::Real)=isinteger(n) ? a^Int(n) : root(a,1//n)

Base.:^(a::Pol, n::Integer)=ismonomial(a) ? Pol_([a.c[1]^n],n*a.v) :
         n>=0 ? Base.power_by_squaring(a,n) :
                Base.power_by_squaring(inv(a),-n)

function Base.:+(a::Pol,b::Pol)
  if iszero(a) return b elseif iszero(b) return a end
  z=zero(a.c[1]+b.c[1])
  d=b.v-a.v
  if d>=0
    res=fill(z,max(length(a.c),d+length(b.c)))
  @inbounds view(res,eachindex(a.c)).=a.c
  @inbounds view(res,d.+eachindex(b.c)).+=b.c
    Pol(res,a.v;copy=false)
  else
    res=fill(z,max(length(a.c)-d,length(b.c)))
  @inbounds view(res,eachindex(a.c).-d).=a.c
  @inbounds view(res,eachindex(b.c)).+=b.c
    Pol(res,b.v;copy=false)
  end
end

function Base.:-(a::Pol,b::Pol)
  if iszero(a) return -b elseif iszero(b) return a end
  z=zero(a.c[1]+b.c[1])
  d=b.v-a.v
  if d>=0
    res=fill(z,max(length(a.c),d+length(b.c)))
  @inbounds view(res,eachindex(a.c)).=a.c
  @inbounds view(res,d.+eachindex(b.c)).-=b.c
    Pol(res,a.v;copy=false)
  else
    res=fill(z,max(length(a.c)-d,length(b.c)))
  @inbounds view(res,eachindex(a.c).-d).=a.c
  @inbounds view(res,eachindex(b.c)).-=b.c
    Pol(res,b.v;copy=false)
  end
end

Base.:+(a::Pol, b::Number)=a+Pol(b)
Base.:+(b::Number, a::Pol)=Pol(b)+a
Base.:-(a::Pol)=Pol_(-a.c,a.v)
Base.:-(b::Number, a::Pol)=Pol(b)-a
Base.:-(a::Pol,b::Number)=a-Pol(b)
Base.div(a::Pol,b::Number)=Pol(div.(a.c,b),a.v;copy=false)

exactdiv(a,b)=a/b  # generic version for fields
function exactdiv(a::Integer,b::Integer) # define for integral domains
  (d,r)=divrem(a,b)
  !iszero(r) ? nothing : d
end

function coeffexactdiv(a::Pol,b)
  if isone(b) return a end
  c=exactdiv.(a.c,b)
  if !any(isnothing,c) Pol_(c,a.v) end
end

Base.:/(p::Pol,q::Number)=Pol_(p.c./q,p.v)
Base.://(p::Pol,q::Number)=Pol_(p.c.//q,p.v)

derivative(a::Pol)=Pol([(i+a.v-1)*v for (i,v) in enumerate(a.c)],a.v-1,copy=false)

"""
`divrem(a::Pol, b::Pol)`

`a` and `b` should be true polynomials (nonnegative valuation).
Computes  `(q,r)` such  that `a=q*b+r`  and `degree(r)<degree(b)`.
Type stable if the coefficients of `b` are in a field.
"""
function Base.divrem(a::Pol, b::Pol)
  if iszero(b) throw(DivideError) end
  if degree(b)>degree(a) return (zero(a),a) end
  if a.v<0 || b.v<0 error("arguments should be true polynomials") end
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

function exactdiv(a::Pol,b::Pol)
  if isone(b) || iszero(a) return a end
  if iszero(b) throw(DivideError) end
  d=a.v-b.v
  if !iszero(a.v) a=shift(a,-a.v) end
  if !iszero(b.v) b=shift(b,-b.v) end
  if degree(b)>degree(a) return nothing end
  z=zero(a.c[1]+b.c[1])
  r=fill(z,1+degree(a))
  view(r,a.v+1:length(r)).=a.c
  q=fill(z,length(r)-degree(b))
  for i in length(r):-1:degree(b)+1
    c=exactdiv(r[i],b.c[end])
    if isnothing(c) return nothing end
    view(r,i-length(b.c)+1:i) .-= c .* b.c
    q[i-length(b.c)+1]=c
  end
  if !iszero(r) return nothing end
  res=Pol(q,d)
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
 if isone(b) || iszero(a) return a,zero(a) end
  if iszero(b) throw(DivideError) end
  d=b.c[end]
  if degree(a)<degree(b) return (zero(a),d^(degree(a)+1-degree(b))*a) end
  if a.v<0 || b.v<0 error("arguments should be true polynomials") end
  z=zero(a.c[1]+b.c[1])
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
  ca=gcd(a.c);a=coeffexactdiv(a,ca)
  cb=gcd(b.c);b=coeffexactdiv(b,cb)
  d=gcd(ca,cb)
  g=1
  h=1
  while true
    δ=degree(a)-degree(b)
    q,r=pseudodiv(a,b)
    if iszero(r)
      cb=gcd(b.c);b=coeffexactdiv(b,cb)
      return isone(d) ? b : Pol_(b.c .*d,b.v)
    elseif degree(r)==0
      return Pol_([d],0)
    end
    a=b
    gh=g*h^δ
    b=coeffexactdiv(r,gh)
    g=a[end]
    if δ>0 h=coeffexactdiv(g^δ,h^(δ-1)) end
  end
end

Base.gcd(p::Pol{<:Integer},q::Pol{<:Integer})=srgcd(p,q)
Base.gcd(v::AbstractArray{<:Pol})=reduce(gcd,v)
Base.lcm(p::Pol,q::Pol)=exactdiv(p*q,gcd(p,q))
Base.lcm(m::AbstractArray{<:Pol})=reduce(lcm,m)

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
  `gcdx(a::Pol,b::Pol)` 

for  polynomials  over  a  field  returns `d,u,v`  such  that `d=ua+vb` and
`d=gcd(a,b)`.

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

"""
`powermod(p::Pol, x::Integer, q::Pol)` computes ``p^x \\pmod m``.
```julia-repl
julia> powermod(q-1,3,q^2+q+1)
Pol{Int64}: 6q+3
```
"""
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
randpol(T,d::Integer)=Pol(rand(T,d+1))

"""
`Pol(x::AbstractVector,y::AbstractVector)`

Interpolation:  find a `Pol` (of  nonnegative valuation) of smallest degree
taking  values `y` at points  `x`. The values `y`  should be in a field for
the function to be type stable.

```julia-repl
julia> p=Pol([1,1,1])
Pol{Int64}: q²+q+1

julia> vals=p.(1:5)
5-element Vector{Int64}:
  3
  7
 13
 21
 31

julia> Pol(1:5,vals*1//1)
Pol{Rational{Int64}}: (1//1)q²+(1//1)q+1//1

julia> Pol(1:5,vals*1.0)
Pol{Float64}: 1.0q²+1.0q+1.0
```
"""
function Pol(pts::AbstractVector,vals::AbstractVector)
  vals=copy(vals)
  a=map(eachindex(pts))do i
    for k in i-1:-1:1
      if pts[i]==pts[k] error("interpolating points must be distinct") end
      vals[k]=(vals[k+1]-vals[k])/(pts[i]-pts[k])
    end
    vals[1]
  end
  p=Pol_([a[end]],0)
  for i in length(pts)-1:-1:1
    p=p*(Pol()-pts[i])+a[i]
  end
  p
end

Base.denominator(p::Pol)=iszero(p) ? 1 : lcm(denominator.(p.c))
Base.numerator(p::Pol{<:Rational{T}}) where T=convert(Pol{T},p*denominator(p))

function root(x::Pol,n::Union{Integer,Rational{<:Integer}}=2)
  n=Int(n)
  if iszero(x) return x end
  if !ismonomial(x) || !iszero(x.v%n) error("root($x,$n) not implemented") end
  c=x.c[1]
  Pol_([isone(c) ? c : root(c,n)],div(x.v,n))
end

root(x::AbstractFloat,n=2)=x^(1/n)
#---------------------- Frac-------------------------------------
## cannot use Rational{T} since it forces T<:Integer
struct Frac{T}
  num::T
  den::T
  # Unexported inner constructor that bypasses all checks
  global Frac_(num::T,den::T) where T=new{T}(num,den)
end

Base.numerator(a::Frac)=a.num
Base.denominator(a::Frac)=a.den

function Base.convert(::Type{Frac{T}},p::Frac{T1}) where {T,T1}
  Frac_(convert(T,p.num),convert(T,p.den))
end

(::Type{Frac{T}})(a::Frac{T}) where T=a
(::Type{Frac{T}})(a::Number) where T=convert(Frac{T},a)

Base.broadcastable(p::Frac)=Ref(p)

Base.copy(a::Frac)=Frac_(a.num,a.den)
Base.one(a::Frac)=Frac_(one(a.num),one(a.den))
Base.one(::Type{Frac{T}}) where T =Frac_(one(T),one(T))
Base.isone(a::Frac)=a.num==a.den
Base.zero(::Type{Frac{T}}) where T =Frac_(zero(T),one(T))
Base.zero(a::Frac)=Frac_(zero(a.num),one(a.num))
Base.iszero(a::Frac)=iszero(a.num)
# next 3 methods are to make inv using LU work (abs is stupid)
Base.abs(p::Frac)=p
Base.conj(p::Frac)=Frac_(conj(p.num),conj(p.den))
Base.adjoint(a::Frac)=conj(a)
Base.:(==)(a::Frac,b::Frac)=a.num==b.num && a.den==b.den
Base.cmp(a::Frac,b::Frac)=cmp((a.num,a.den),(b.num,b.den))
Base.isless(a::Frac,b::Frac)=cmp(a,b)==-1

function Base.show(io::IO, ::MIME"text/plain", a::Frac)
  if !haskey(io,:typeinfo) 
    print(io,typeof(a),": ") 
    io=IOContext(io,:typeinfo=>typeof(a))
  end
  show(io,a)
end

function Base.show(io::IO,a::Frac)
  if !get(io, :limit, false) && !get(io, :TeX, false)
    print(io,"Frac(",a.num,",",a.den,")")
    return
  end
  if haskey(io,:typeinfo) && isone(a.den) 
    print(io,a.num)
    return 
  end
  n=sprint(show,a.num; context=io)
  print(io,bracket_if_needed(n))
  n=sprint(show,a.den; context=io)
  print(io,"/",bracket_if_needed(n))
end

Base.inv(a::Frac)=Frac_(a.den,a.num)

Base.:*(a::Frac,b::Frac)=Frac(a.num*b.num,a.den*b.den)
Base.://(a::Frac,b::Frac)=Frac(a.num*b.den,a.den*b.num)
Base.:/(a::Frac,b::Frac)=a//b

Base.:^(a::Frac, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                              Base.power_by_squaring(inv(a),-n)
Base.:+(a::Frac,b::Frac)=Frac(a.num*b.den+a.den*b.num,a.den*b.den)
Base.:-(a::Frac)=Frac_(-a.num,a.den)
Base.:-(a::Frac,b::Frac)=Frac(a.num*b.den-a.den*b.num,a.den*b.den)
Base.:+(a::Frac{<:T},b::Union{Number,T}) where T=+(promote(a,b)...)
Base.:+(b::Union{Number,T},a::Frac{<:T}) where T=+(promote(a,b)...)
Base.:-(a::Frac{<:T},b::Union{Number,T}) where T=-(promote(a,b)...)
Base.:-(b::Union{Number,T},a::Frac{<:T}) where T=-(promote(b,a)...)
Base.:*(a::Frac{<:T},b::T) where T=Frac(a.num*b,a.den)
Base.:*(b::T,a::Frac{<:T}) where T=Frac(a.num*b,a.den)

#----------------------------------------------------------------------

# make both pols, one of valuation 0
function make_positive(a::Pol,b::Pol)
  v=a.v-b.v
  shift(a,max(v,0)-a.v),shift(b,-min(v,0)-b.v)
end
  
Frac(a::Pol,b::Pol;prime=false)=Frac(promote(a,b)...;prime)

"""
`Frac(a::Pol,b::Pol;prime=false)

Polynomials  `a` and `b` are promoted to same coefficient type, and checked
for  being true polynomials (otherwise they are both multiplied by the same
power  of  the  variable  so  they  become  true  polynomials),  and unless
`prime=true` they are checked for having a non-trivial `gcd`.
"""
function Frac(a::T,b::T;prime=false)::Frac{T} where T<:Pol
  if iszero(a) return Frac_(a,one(a))
  elseif iszero(b) error("division by 0")
  end
  a,b=make_positive(a,b)
  if !(prime || ismonomial(b) || ismonomial(a))
    d=gcd(a,b)
    a,b=exactdiv(a,d),exactdiv(b,d)
  end
  if scalar(b)==-1 a,b=(-a,-b) end
  Frac_(a,b)
end

function Frac(a::Pol{<:Rational},b::Pol{<:Rational};k...)
  Frac(numerator(a)*denominator(b),numerator(b)*denominator(a);k...)
end

function Pol(p::Frac{<:Pol})
  if ismonomial(p.den)==1
    if isone(p.den.c[1]^2) return Pol(p.num.c .*p.den.c[1],p.num.v-p.den.v)
    else return Pol(p.num.c .//p.den.c[1],p.num.v-p.den.v)
    end
  end
  error("cannot convert ",p," to Pol")
end

Base.convert(::Type{Pol{T}},p::Frac) where {T}=convert(Pol{T},Pol(p))

function Base.convert(::Type{Frac{T}},p::Pol) where T
  Frac(convert(T,p),one(T);prime=true)
end

function Base.convert(::Type{Frac{T}},p::Number) where {T<:Pol}
  Frac_(convert(T,p),convert(T,1))
end

function Base.promote_rule(a::Type{T1},b::Type{Frac{T2}})where {T1<:Pol,T2<:Pol}
  Frac{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{T1},b::Type{Frac{T2}})where {T1<:Number,T2<:Pol}
  Frac{promote_type(T1,T2)}
end

function Frac(a::Pol)
  if a.v>0 return Frac_(a,one(a)) end
  Frac(a,one(a);prime=true)
end

function Base.://(a::Pol,b::Pol)
  if ismonomial(b) return Pol(a.c.//b.c[1],a.v-b.v) end
  Frac(a,b)
end

bestinv(x)=isone(x) ? x : isone(-x) ? x : inv(x)

function Base.inv(p::Pol)
  if ismonomial(p) return Pol_([bestinv(p.c[1])],-p.v) end
  Frac_(make_positive(one(p),p)...)
end

Base.:/(p::Pol,q::Pol)=p*inv(q)

Base.://(a::Frac{<:Pol},b::Pol)=Frac(a.num,a.den*b)
Base.://(a::Frac{<:Pol},b::Number)=Frac(a.num,a.den*b;prime=true)
Base.:/(a::Frac{<:Pol},b::Union{Number,Pol})=a//b
Base.://(a::Union{Number,Pol},b::Frac{<:Pol})=a*inv(b)
Base.:/(a::Union{Number,Pol},b::Frac{<:Pol})=a//b
Base.://(p::Number,q::Pol)=Frac(Pol(p),q;prime=true)
Base.:/(p::Number,q::Pol)=p*inv(q)

Base.:*(a::Frac{<:Pol},b::Number)=Frac(a.num*b,a.den;prime=true)
Base.:*(b::Number,a::Frac{<:Pol})=a*b

(p::Frac{<:Pol})(x;Rational=false)=Rational ? p.num(x)//p.den(x) : p.num(x)/p.den(x)
# @btime inv(Frac.([q+1 q+2;q-2 q-3])) setup=(q=Pol())
# 32.556 μs (938 allocations: 57.22 KiB)
end
