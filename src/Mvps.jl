"""
What   is  implemented  here  is  "Puiseux  polynomials",  that  is  linear
combinations  of monomials  of the  type `x₁^{a₁}…  xₙ^{aₙ}` where `xᵢ` are
variables  and `aᵢ` are  exponents which can  be arbitrary rational numbers
(we   need  them  with  cyclotomic  coefficients  as  splitting  fields  of
cyclotomic  Hecke  algebras).  Some  functions  described  below need their
argument  to involve  only variables  to integral  powers; we will refer to
such  objects as "Laurent polynomials"; some functions require further that
variables  are  raised  only  to  positive  powers:  we refer then to "true
polynomials".

`@Mvp x₁,…,xₙ`

declares   that  `xᵢ`are  indeterminates  suitable  to  build  multivariate
polynomials.

```julia-repl
julia> @Mvp x,y

julia> (x+y)^3
Mvp{Int64}: x³+3x²y+3xy²+y³
```

`Mvp(x::Number)`   returns  the  constant   multivariate  polynomial  whose
constant term is `x`.

```julia-repl
julia> degree(Mvp(1))
0
```

One can divide an `Mvp` by another when the division is exact
(this is equivalent to `ExactDiv`, see below).

```julia-repl
julia> (x^2-y^2)//(x-y)
Mvp{Int64}: x+y
```

Only monomials can be raised to a non-integral power; they can be raised to
a  fractional power of denominator `b` only if `root(x,b)` is defined where
`x`   is  their  leading  coefficient.  For  an  `Mvp`  `m`,  the  function
`root(m,n)`  is equivalent  to `m^(1//n)`.  Raising a  non-monomial Laurent
polynomial to a negative power returns a rational fraction.

```julia-repl
julia> (4x)^(1//2)
Mvp{Int64,Rational{Int64}}: 2x½

julia> (2.0x)^(1//2)
Mvp{Float64,Rational{Int64}}: 1.4142135623730951x½

julia> root(2.0x)
Mvp{Float64,Rational{Int64}}: 1.4142135623730951x½
```

Despite the degree of generality of our polynomials, we are not too shabby.
For  the Fateman test  f(f+1) where f=(1+x+y+z+t)^n  with n=15, we take 4.5
sec. According to the Nemo paper, Sagemath takes 10sec and Nemo takes 1.6.
"""
module Mvps
# benchmark: (x+y+z)^3     2.7μs 141 alloc
using ..ModuleElts: ModuleElt, ModuleElts
using ..Util: fromTeX, ordinal, printTeX

#import Gapjm: degree, coefficients, valuation
#import ..Pols: positive_part, negative_part, bar
import ..Cycs: root
# to use as a stand-alone module comment above line, uncomment next, and
# define root for the coefficients you want (at least root(1,n)=1)
#export root
export degree, coefficients, valuation
export positive_part, negative_part, bar
export Mvp, Monomial, @Mvp, variables, value, scal, derivative,
  laurent_denominator
#------------------ Monomials ---------------------------------------------
struct Monomial{T}
  d::ModuleElt{Symbol,T}   
end

Monomial(a::Pair...)=Monomial(ModuleElt(a...))
Monomial()=one(Monomial{Int})
Monomial(v::Symbol)=convert(Monomial{Int},v)

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
Base.div(a::Monomial, b::Monomial)=a*inv(b)
Base.:/(a::Monomial, b::Monomial)=a*inv(b)
Base.://(a::Monomial, b::Monomial)=a*inv(b)
Base.:^(x::Monomial,p)=Monomial(x.d*p)
Base.getindex(a::Monomial,k)=getindex(a.d,k)

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

# cmp must define a monomial order (a<b => a*m<b*m for all monomials a, b, m)
# We take a<b if the first variable in a/b occurs to a negative power
function Base.cmp(a::Monomial, b::Monomial)
  for ((ca,pa),(cb,pb)) in zip(a.d,b.d)
    if ca!=cb return ca<cb ? -sign(pa) : sign(pb) end
    if pa!=pb return cmp(pb,pa) end
  end
  la=length(a.d)
  lb=length(b.d)
  if    la==lb return 0
  elseif la>lb return -sign(last(a.d.d[lb+1]))
  else         return sign(last(b.d.d[la+1]))
  end
end

Base.isless(a::Monomial, b::Monomial)=cmp(a,b)<0
Base.:(==)(a::Monomial, b::Monomial)=a.d==b.d

function Base.gcd(l::Monomial...)
  if length(l)!=2 reduce(gcd,l;init=Monomial())
  else Monomial(merge(min,l[1].d,l[2].d))
  end
end

function Base.lcm(l::Monomial...)
  if length(l)!=2 reduce(lcm,l;init=Monomial())
  else Monomial(merge(max,l[1].d,l[2].d))
  end
end

Base.hash(a::Monomial, h::UInt)=hash(a.d,h)

degree(m::Monomial)=sum(values(m.d))
degree(m::Monomial,var::Symbol)=m[var]

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
Mvp()=Mvp(ModuleElt(Pair{Monomial{Int},Int}[]))

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

# necessary to clean up a ModuleElt upstream
Base.show(io::IO, x::Mvp)=show(IOContext(io,:showbasis=>nothing),x.d)

Base.zero(p::Mvp)=Mvp(zero(p.d))
Base.zero(::Type{Mvp{T,N}}) where {T,N}=Mvp(ModuleElt(Pair{Monomial{N},T}[]))
Base.one(::Type{Mvp{T,N}}) where {T,N}=Mvp(one(Monomial{N})=>one(T))
Base.one(::Type{Mvp{T}}) where T=one(Mvp{T,Int})
Base.one(::Type{Mvp})=Mvp(1)
Base.one(p::Mvp{T,N}) where {T,N}=one(Mvp{T,N})
Base.copy(p::Mvp)=Mvp(p.d)
Base.iszero(p::Mvp)=iszero(p.d)
Base.convert(::Type{Mvp},a::Number)=convert(Mvp{typeof(a),Int},a)
Base.convert(::Type{Mvp{T,N}},a::Number) where {T,N}=iszero(a) ? 
 zero(Mvp{T,N}) : Mvp(one(Monomial{N})=>convert(T,a))
Base.convert(::Type{Mvp{T,N}},a::Mvp{T1,N1}) where {T,T1,N,N1}=
  Mvp(convert(ModuleElt{Monomial{N},T},a.d))
Base.convert(::Type{Mvp},x::Mvp)=x
Base.convert(::Type{Mvp},v::Symbol)=Mvp(Monomial(v)=>1)
Mvp(x::Symbol)=convert(Mvp,x)
Mvp(x::Number)=convert(Mvp,x)
Mvp(x::Mvp)=convert(Mvp,x)

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
Base.:*(a::Monomial, b::Mvp)=Mvp(ModuleElt(m*a=>c for (m,c) in b.d))
Base.:*(b::Mvp,a::Monomial)=a*b
function Base.:*(a::Mvp, b::Mvp)
  if length(a.d)>length(b.d) a,b=(b,a) end
  if iszero(a) return a end
  let b=b # needed !!!!
    sum(b*m*c for (m,c) in a.d)
  end
end

Base.:(//)(a::Mvp, b::Number)=Mvp(ModuleElt(m=>c//b for (m,c) in a.d))
Base.:(/)(a::Mvp, b::Number)=Mvp(ModuleElt(m=>c/b for (m,c) in a.d))

"""
`conj(p::Mvp)` acts on the coefficients of `p`

```julia-repl
julia> @Mvp x;conj(im*x+1)
Mvp{Complex{Int64}}: (0 - 1im)x+1 + 0im
```
"""
Base.conj(a::Mvp)=Mvp(ModuleElt(m=>conj(c) for (m,c) in a.d))

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
degree(m::Mvp)=iszero(m) ? 0 : maximum(degree.(keys(m.d)))
degree(m::Mvp,v::Symbol)=iszero(m) ? 0 : maximum(degree.(keys(m.d),v))

"""
The `valuation` of an `Mvp` is the minimal degree of a monomial.


```julia-repl
julia> a=x^2+x*y
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
valuation(m::Mvp)=iszero(m) ? 0 : minimum(degree.(keys(m.d)))
valuation(m::Mvp,v::Symbol)=iszero(m) ? 0 : minimum(degree.(keys(m.d),v))

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
function coefficients(p::Mvp,v::Symbol)
  if iszero(p) return Dict{Int,typeof(p)}() end
  d=Dict{Int,typeof(p.d.d)}()
  for (m,c) in p.d
#   print("$m=>$c d=$d\n")
    found=false
    for (i,(v1,deg)) in enumerate(m.d)
      if v1==v 
        found=true
        d[deg]=push!(get(d,deg,empty(p.d.d)),
                        Monomial(deleteat!(collect(m.d),i)...)=>c)
      end
    end
    if !found  d[0]=push!(get(d,0,empty(p.d.d)),m=>c) end
  end
  Dict(dg=>Mvp(c...) for (dg,c) in d) # c... is sorted by defn of monomial order
end

"""
`variables(p::Mvp...)`

returns the list of variables of all `p` as a sorted list of `Symbol`s.

```julia-repl
julia> @Mvp x,y,z

julia> variables(x+y+1,z)
3-element Vector{Symbol}:
 :x
 :y
 :z
```
"""
variables(pp::Mvp...)=unique!(sort([k1 for p in pp for (k,v) in p.d for (k1,v1) in k.d]))

"""
`scal(p::Mvp)`

If  `p`  is a  scalar,  return that  scalar,
otherwise return  `nothing`. 

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
`value(p::Mvp,:x1=>v1,:x2=>v2,...)`

gives  the value  of `p`  when doing  the simultaneous  substitution of the
variable  `:x1`  by  `v1`,  of  `x2`  by  `v2`,  … This can also be written
`p(;x1=v1,x2=v2,...)`.

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
  if badi!==nothing Mvp(deleteat!(collect(p.d),badi)...)+res1
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

function Base.inv(p::Mvp)
  if length(p.d)>1 || iszero(p) throw(InexactError(:inv,Mvp,p)) end
  (m,c)=first(p.d)
  if c==1 || c==-1 return Mvp(inv(m)=>c) end
  Mvp(inv(m)=>1//c)
end

function Base.://(p::Mvp,q::Mvp)
  res=ExactDiv(p,q)
  if isnothing(res) return p*inv(q) end
  res
end
Base.:/(p::Mvp,q::Mvp)=p//q
Base.://(p::Number,q::Mvp)=Mvp(p)//q

function ExactDiv(p::Mvp,q::Mvp)
# println("p=$p q=$q")
  if iszero(q) error("cannot divide by 0")
  elseif length(q.d)==1 return p*inv(q)
  elseif length(p.d)<2 return iszero(p) ? p : nothing
  end 
  var=first(first(p.d)[1].d)[1]
  res=zero(p)
  cq=coefficients(q,var)
  mq=maximum(keys(cq))
# println("var=$var mq=$mq cq=$cq q=$q")
  while length(p.d)!=0
    if length(p.d)<length(q.d) return nothing end
    cp=coefficients(p,var)
    mp=maximum(keys(cp))
#   println("mp=$mp cp=$cp")
    t=ExactDiv(cp[mp],cq[mq])
    if t===nothing return nothing end
    if mp!=mq t=Monomial(var=>mp-mq)*t end
    res+=t
#   print("t=$t res=$res p=$p=>")
    p-=t*q
#   println("$p")
  end
  res
end

function Base.gcd(a::Mvp,b::Mvp)
  if any(x->length(x.d)==1,[a,b]) 
    return Mvp([gcd(reduce(vcat,map(x->x[1],a.d),map(x->x[1],b.d)))=>1])
  end
  listvar=(variables(a),variables(b))
  v=symdiff(listvar...)
  listvar=union(listvar...)
  nvar=length(listvar)
# if nvar==2 Print("v=",v," listvar=",listvar,"\n");end
  v=length(v)>0 ? v[1] : listvar[1]
  coef=[coefficients(a,v),coefficients(b,v)]
# if nvar==2 Print("coef=",coef,"\n");end

  function Vecmod(p, q)
    lq=length(q)
    if lq==1 return [] end
    qlq=q[lq]
    lp=length(p)
    p=copy(p)
    f=scal(qlq)
    q=q[1:lq-1]
    if f!=false q=q/f end
    while lp>=lq
      plp=p[lp]
      lp-=1 
      if f==false p=p*qlq end
      p[lp-lq+2:lp]=p[lp-lq+2:lp]-plp*q
      while lp>0 && p[lp]==0*p[lp] lp-=1 end
    end
    p=p[1:lp]
    if lp==0 return p end
    plp=scal(p[lp])
    if plp==false return p/ApplyFunc(MvpGcd,p) end
    p/plp
  end

  function VecGcd(p,q)
    while length(q)>0 (q,p)=(Vecmod(p,q),q) end
    p/p[end]
  end
  
  if nvar==1 return Mvp(v,VecGcd(coef[1],coef[2])) end # faster
  cont=map(x->MvpGcd(x...),coef)
  for i in 1:2
     if scal(cont[i]==false) coef[i]=List(coef[i],x->ExactDiv(x,cont[i]))
     else coef[i]=List(coef[i],Mvp)
          cont[i]=Mvp(cont[i])
     end
  end
  reseuc=VecGcd(coef[1],coef[2])
  l=MvpLcm(map(x->RatFrac(x).den,reseuc)...)
  reseuc=map(p->IsRatFrac(p) ? p.num*ExactDiv(l,p.den) : p*l,reseuc)
# if not ForAll(reseuc,IsMvp) Error() end
# if not IsMvp(MvpGcd(cont[1],cont[2])) Error() end
  MvpGcd(cont[1],cont[2])* ValuePol(reseuc,Mvp(v))
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

positive_part(p::Mvp)=
  Mvp(ModuleElt([m=>c for (m,c) in p.d if all(>(0),values(m.d))]))

negative_part(p::Mvp)=
  Mvp(ModuleElt([m=>c for (m,c) in p.d if all(<(0),values(m.d))]))

bar(p::Mvp)=Mvp((inv(m)=>c for (m,c) in p.d)...)

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
function derivative(p::Mvp,v=first(variables(p)))
  Mvp(ModuleElt([m*Monomial(v)^-1=>c*m[v] for (m,c) in p.d];check=true))
end

"""
`laurent_denominator(p1,p2,…)`

returns  the unique monomial  `m` of minimal  degree such that  for all the
Laurent  polynomial  arguments  `p1,p2,…`  the  product  `m*pᵢ`  is  a true
polynomial.

```julia-repl
julia> laurent_denominator(x^-1,y^-2+x^4)
Monomial{Int64}:xy²
```
"""
function laurent_denominator(p::Mvp...)
  res=Monomial()
  for v in variables(p...)
    val=minimum(valuation.(p,v))
    if val<0 res*=Monomial(v)^-val end
  end
  return res
end


"""
```julia-repl
julia> p=E(3)x+E(5)
Mvp{Cyc{Int64}}: ζ₃x+ζ₅

julia> q=convert(Mvp{ComplexF64,Int},p)
Mvp{ComplexF64}: (-0.4999999999999998 + 0.8660254037844387im)x+0.30901699437494745 + 0.9510565162951535im

julia> conj(q)
Mvp{ComplexF64}: (-0.4999999999999998 - 0.8660254037844387im)x+0.30901699437494745 - 0.9510565162951535im
```
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The next functions have been provided by Gwenaëlle Genet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'MvpGcd( <p1>, <p2>, ...)'

Returns  the Gcd  of the  'Mvp' arguments.  The arguments  must be  true
polynomials.

|    gap> MvpGcd(x^2-y^2,(x+y)^2);
    x+y|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'MvpLcm( <p1>, <p2>, ...)'

Returns  the Lcm  of the  'Mvp' arguments.  The arguments  must be  true
polynomials.

|    gap> MvpLcm(x^2-y^2,(x+y)^2);
    xy^2-x^2y-x^3+y^3|
"""

#julia1.6.0-rc1> @btime Mvps.fateman(15)
# 4.523 s (15219388 allocations: 5.10 GiB)
function fateman(n)
  f=(1+Mvp(:x)+Mvp(:y)+Mvp(:z)+Mvp(:t))*one(n)
  f=f^n
  length((f*(f+1)).d)
end

end
