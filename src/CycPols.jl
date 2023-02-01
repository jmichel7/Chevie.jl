"""
This package deals with products of Cyclotomic polynomials.

Cyclotomic  numbers, and cyclotomic polynomials  over the rationals or some
cyclotomic  field, are important in the theories of finite reductive groups
and  Spetses. In particular Schur elements of cyclotomic Hecke algebras are
products of cyclotomic polynomials.

The  type  `CycPol`  represents  the  product  of  a `coeff` (a constant, a
polynomial or a rational fraction in one variable) with a rational fraction
in  one variable with all poles or zeroes equal to 0 or roots of unity. The
advantages  of  representing  as  `CycPol`  such  objects are: nice display
(factorized),  less storage, fast  multiplication, division and evaluation.
The drawback is that addition and subtraction are not implemented!

This    package   depends   on   the   packages   `Primes`,   `ModuleElts`,
`CyclotomicNumbers`, `LaurentPolynomials` and `Combinat`.

The  method  `CycPol(a::Pol)`  converts  `a`  to  a `CycPol` by finding the
largest  cyclotomic polynomial  dividing, leaving  a `Pol` `coefficient` if
some roots of the polynomial are not roots of unity.

```julia-repl
julia> @Pol q
Pol{Int64}: q

julia> p=CycPol(q^25-q^24-2q^23-q^2+q+2) # a `Pol` coefficient remains
(q-2)Φ₁Φ₂Φ₂₃

julia> p(q) # evaluate CycPol p at q
Pol{Int64}: q²⁵-q²⁴-2q²³-q²+q+2

julia> p*inv(CycPol(q^2+q+1)) # `*`, `inv`, `/` and `//` are defined
(q-2)Φ₁Φ₂Φ₃⁻¹Φ₂₃

julia> -p  # one can multiply by a scalar
(-q+2)Φ₁Φ₂Φ₂₃

julia> valuation(p)
0

julia> degree(p)
25

julia> lcm(p,CycPol(q^3-1)) # lcm is fast between CycPols
(q-2)Φ₁Φ₂Φ₃Φ₂₃
```

```julia-rep1
julia> print(p)
CycPol(Pol([-2, 1]),0,(1,0),(2,0),(23,0)) # a format which can be read in Julia
```
Evaluating  a `CycPol` at some `Pol` value  gives in general a `Pol`. There
are  exceptions  where  we  can  keep  the  value a `CycPol`: evaluating at
`Pol()^n`  (that is `q^n`) or at  `Pol([E(n,k)],1)` (that is `qζₙᵏ`). Then
`subs` gives that evaluation:

```julia-repl
julia> subs(p,Pol()^-1) # evaluate as a CycPol at q⁻¹
(2-q⁻¹)q⁻²⁴Φ₁Φ₂Φ₂₃

julia> subs(p,Pol([E(2)],1)) # or at -q
(-q-2)Φ₁Φ₂Φ₄₆
```
The variable name used when printing a `CycPol` is the same as for `Pol`s.

When  showing  a  `CycPol`,  some  factors  over  extension  fields  of the
cyclotomic polynomial `Φₙ` are given a special name. If `n` has a primitive
root  `ξ`, `ϕ′ₙ` is the product of the  `(q-ζ)` where `ζ` runs over the odd
powers  of `ξ`, and `ϕ″ₙ` is the  product for the even powers. Some further
factors are recognized for small `n`. 
```julia-repl
julia> CycPol(q^6-E(4))
Φ″₈Φ⁽¹³⁾₂₄
```
The  function `show_factors` gives the  complete list of recognized factors
for a given `n`:
```julia-repl
julia> CycPols.show_factors(24)
15-element Vector{Tuple{CycPol{Int64}, Pol}}:
 (Φ₂₄, q⁸-q⁴+1)
 (Φ′₂₄, q⁴+ζ₃²)
 (Φ″₂₄, q⁴+ζ₃)
 (Φ‴₂₄, q⁴-√2q³+q²-√2q+1)
 (Φ⁗₂₄, q⁴+√2q³+q²+√2q+1)
 (Φ⁽⁵⁾₂₄, q⁴-√6q³+3q²-√6q+1)
 (Φ⁽⁶⁾₂₄, q⁴+√6q³+3q²+√6q+1)
 (Φ⁽⁷⁾₂₄, q⁴+√-2q³-q²-√-2q+1)
 (Φ⁽⁸⁾₂₄, q⁴-√-2q³-q²+√-2q+1)
 (Φ⁽⁹⁾₂₄, q²+ζ₃²√-2q-ζ₃)
 (Φ⁽¹⁰⁾₂₄, q²-ζ₃²√-2q-ζ₃)
 (Φ⁽¹¹⁾₂₄, q²+ζ₃√-2q-ζ₃²)
 (Φ⁽¹²⁾₂₄, q²-ζ₃√-2q-ζ₃²)
 (Φ⁽¹³⁾₂₄, q⁴-ζ₄q²-1)
 (Φ⁽¹⁴⁾₂₄, q⁴+ζ₄q²-1)
```
Such a factor can be obtained directly as:

```julia-repl
julia> CycPol(;conductor=24,no=7)
Φ⁽⁷⁾₂₄

julia> CycPol(;conductor=24,no=7)(q)
Pol{Cyc{Int64}}: q⁴+√-2q³-q²-√-2q+1
```
This package also defines the function `cylotomic_polynomial`:
```julia-repl
julia> p=cyclotomic_polynomial(24)
Pol{Int64}: q⁸-q⁴+1

julia> CycPol(p) # same as CycPol(;conductor=24,no=0)
Φ₂₄
```
"""
module CycPols

export CycPol, cyclotomic_polynomial, subs

using Primes: primes, factor, eachfactor, totient #Euler φ
using ModuleElts: ModuleElts, ModuleElt
using LaurentPolynomials: stringexp, bracket_if_needed
using CyclotomicNumbers: CyclotomicNumbers, Root1, E, conductor, Cyc, order
using LaurentPolynomials: Pol, LaurentPolynomials, degree, valuation,
                          coefficients, pseudodiv, exactdiv, Frac
using Combinat: primitiveroot, divisors, collectby

Base.numerator(p::Pol{<:Integer})=p  # to put in LaurentPolynomials
Base.numerator(p::Pol{Cyc{Rational{T}}}) where T<:Integer =
  Pol{Cyc{T}}(p*denominator(p))
Base.numerator(p::Pol{Cyc{T}}) where T<:Integer =p
CyclotomicNumbers.conductor(x::Pol)=lcm(conductor.(coefficients(x)))

using CyclotomicNumbers: prime_residues, stringind, format_coefficient

function stringprime(io::IO,n)
  if iszero(n) return "" end
  if get(io,:TeX,false) return "'"^n end
  n<5 ? ["′","″","‴","⁗"][n] : "⁽"*stringexp(io,n)*"⁾"
end
  
# The  computed  cyclotomic  polynomials  are  cached 
const cyclotomic_polynomial_dict=Dict(1=>Pol([-1,1]))
"""
`cyclotomic_polynomial(n)`
 
returns the `n`-th cyclotomic polynomial.
 
```julia-repl
julia> cyclotomic_polynomial(5)
Pol{Int64}: q⁴+q³+q²+q+1

julia> cyclotomic_polynomial(24)
Pol{Int64}: q⁸-q⁴+1
```
"""
function cyclotomic_polynomial(n::Integer)
  get!(cyclotomic_polynomial_dict,n) do
    res=Pol(fill(1,n),0;check=false)
    for d in divisors(n)
      if d!=1 && d!=n
        res,_=pseudodiv(res,cyclotomic_polynomial(d))
      end
    end
    res
  end::Pol{Int}
end

"""
`CycPol`s are internally a `struct` with fields:

`.coeff`:  a coefficient, usually a `Cyc` or a `Pol`. The `Pol` case allows
   to represent as `CycPol`s arbitrary `Pol`s which is useful sometimes.

`.valuation`: an `Int`.

`.v`: a ModuleElt{Root1,Int} where pairs `ζ=>m` give multiplicity `m` of `ζ`.

So `CycPol(coeff,val,v)` represents `coeff*q^val*prod((q-ζ)^m for (ζ,m) in v)`.
"""
struct CycPol{T}
  coeff::T
  valuation::Int
  v::ModuleElt{Root1,Int}
end

# CycPols are scalars for broadcasting
Base.broadcastable(p::CycPol)=Ref(p)

Base.convert(::Type{CycPol{T1}},p::CycPol{T}) where {T1,T}=T==T1 ? p :
                               CycPol(T1(p.coeff),p.valuation,p.v)

Base.hash(a::CycPol, h::UInt)=hash(a.coeff,hash(a.valuation,hash(a.v,h)))

function Base.cmp(a::CycPol,b::CycPol)
  res=cmp(a.valuation,b.valuation)
  if !iszero(res) return res end
  res=cmp(a.v,b.v)
  if !iszero(res) return res end
  cmp(a.coeff,b.coeff)
end

Base.isless(a::CycPol,b::CycPol)=cmp(a,b)==-1

Base.:(==)(a::CycPol,b::CycPol)=cmp(a,b)==0

# see if check should be false
#CycPol(c,val::Int,v::Pair{Rational{Int},Int}...;check=false)=CycPol(c,val,
#  ModuleElt(Pair{Root1,Int}[Root1(;r=r)=>m for (r,m) in v];check)) 

CycPol(c,val::Int,tt::Tuple...;check=false)=CycPol(c,val,sum(x->v(x...),tt;
                                   init=zero(ModuleElt{Root1,Int})))

Base.one(::Type{CycPol})=CycPol(1,0)
Base.one(p::CycPol)=CycPol(one(p.coeff),0)
Base.isone(p::CycPol)=isone(p.coeff) && iszero(p.valuation) && iszero(p.v)
Base.zero(::Type{CycPol{T}}) where T=CycPol(zero(T),0)
Base.zero(::Type{CycPol})=zero(CycPol{Int})
Base.zero(a::CycPol)=CycPol(zero(a.coeff),0)
Base.iszero(a::CycPol)=iszero(a.coeff)
Base.copy(a::CycPol)=CycPol(a.coeff,a.valuation,a.v)

LaurentPolynomials.degree(a::CycPol)=reduce(+,values(a.v);init=0)+a.valuation+degree(a.coeff)
LaurentPolynomials.valuation(a::CycPol)=a.valuation
LaurentPolynomials.valuation(a::CycPol,d::Root1)=reduce(+,c for (r,c) in a.v if r==d;init=0)

function Base.:*(a::CycPol,b::CycPol)
  if iszero(a) || iszero(b) return zero(a) end
  CycPol(a.coeff*b.coeff,a.valuation+b.valuation,a.v+b.v)
end

Base.conj(a::CycPol)=CycPol(conj(a.coeff),a.valuation,ModuleElt(
           inv(r)=>m for (r,m) in a.v))
Base.transpose(a::CycPol)=a
Base.:*(a::CycPol,b::Number)=iszero(b) ? zero(a) : CycPol(a.coeff*b,a.valuation,a.v)
Base.:*(b::Number,a::CycPol)=a*b
Base.:-(a::CycPol)=CycPol(-a.coeff,a.valuation,a.v)
Base.inv(a::CycPol)=CycPol(LaurentPolynomials.bestinv(a.coeff), -a.valuation, -a.v)
Base.:^(a::CycPol, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                      Base.power_by_squaring(inv(a),-n)
Base.://(a::CycPol,b::CycPol)=CycPol(a.coeff//b.coeff, a.valuation-b.valuation, a.v-b.v)
Base.://(a::CycPol,b::Number)=CycPol(a.coeff//b,a.valuation,a.v)
Base.:/(a::CycPol,b::CycPol)=CycPol(a.coeff/b.coeff, a.valuation-b.valuation, a.v-b.v)
Base.:/(a::CycPol,b::Number)=CycPol(a.coeff/b,a.valuation,a.v)
Base.:div(a::CycPol,b::Number)=CycPol(div(a.coeff,b),a.valuation,a.v)

function Base.lcm(a::CycPol,b::CycPol) # forgets .coeff
  if (b.coeff isa Pol && !iszero(degree(b))) 
    error(b,".coeff should be scalar")
  end
  CycPol(a.coeff,max(a.valuation,b.valuation),ModuleElts.merge2(max,a.v,b.v))
end

Base.lcm(v::CycPol...)=reduce(lcm,collect(v);init=one(CycPol))

Base.lcm(v::AbstractArray{<:CycPol})=reduce(lcm,v;init=one(CycPol))

const dec_dict=Dict(1=>[[1]],2=>[[1]],
  8=>[[1,3,5,7],[1,5],[3,7],[1,7],[3,5],[1,3],[5,7]],
 12=>[[1,5,7,11],[1,5],[7,11],[1,7],[5,11],[1,11],[5,7]],
 15=>[[1,2,4,7,8,11,13,14],[1,4,11,14],[2,7,8,13],[1,4,7,13],[2,8,11,14],
      [1,4],[7,13],[11,14],[2,8]],
 16=>[[1,3,5,7,9,11,13,15],[1,7,9,15],[3,5,11,13],[1,5,9,13],[3,7,11,15]],
 20=>[[1,3,7,9,11,13,17,19],[1,9,11,19],[3,7,13,17],[1,9,13,17],[3,7,11,19]],
 21=>[[1,2,4,5,8,10,11,13,16,17,19,20],[1,4,10,13,16,19],[2,5,8,11,17,20]],
 24=>[[1,5,7,11,13,17,19,23],[1,7,13,19],[5,11,17,23],[1,7,17,23],[5,11,13,19],
      [1,5,19,23],[7,11,13,17],[1,11,17,19],[5,7,13,23],
      [7,13],[1,19],[5,23],[11,17],[1,5,13,17],[7,11,19,23]],
 30=>[[1,7,11,13,17,19,23,29],[1,11,19,29],[7,13,17,23],[1,7,13,19],
      [11,17,23,29],[11,29],[17,23],[1,19],[7,13]],
 36=>[[1,5,7,11,13,17,19,23,25,29,31,35],[1,7,13,19,25,31],[5,11,17,23,29,35],
      [7,11,19,23,31,35],[1,5,13,17,25,29]],
42=>[[1,5,11,13,17,19,23,25,29,31,37,41],[1,13,19,25,31,37],[5,11,17,23,29,41]])

# returns list of subsets of primitive_roots(d) wich have a "name" Φ_d^(i)
function dec(d::Int)
  get!(dec_dict,d) do
    dd=[prime_residues(d)]
    if (r=primitiveroot(d))!==nothing
      for a in 0:1
        push!(dd,sort(powermod.(r,(0:2:totient(d)-2).+a,d)))
      end
    end
    dd
  end::Vector{Vector{Int}}
end

v(conductor=1,no=0,mul=1)=no<0 ? ModuleElt(E(conductor,-no)=>mul) : 
 ModuleElt(map(i->E(conductor,i)=>mul,dec(conductor)[no+1]);check=false)

CycPol(;conductor=1,no=0)=CycPol(1,0,v(conductor,no))
  
function show_factors(d)
 map(eachindex(CycPols.dec(d))) do i
   p=CycPol(;conductor=d,no=i-1)
   (p,p(Pol()))
 end
end

pr()=vcat(map(show_factors,sort(collect(keys(dec_dict))))...)

# decompose the .v of a CycPol in subsets Φ^i (used for printing and value)
function decompose(v::Vector{Pair{Root1,Int}})
  rr=@NamedTuple{conductor::Int, no::Int, mul::Int}[]
  for t in collectby(x->order(first(x)),v)
    c=order(first(t[1]))
    if c==1 
      push!(rr,(conductor=c,no=0,mul=last(t[1])))
      continue
    end
    res=@NamedTuple{conductor::Int, no::Int, mul::Int}[]
    v=fill(0,c)
    @views v[exponent.(first.(t))].=last.(t)
    for (i,r) in enumerate(dec(c))
      if (n=minimum(@view v[r]))>0 || (n=maximum(@view v[r]))<0 
        @views v[r].-=n 
        push!(res,(conductor=c,no=i-1,mul=n))
      end
    end
    for i in 1:c  
      if v[i]!=0 push!(res,(conductor=c,no=-i,mul=v[i])) end 
    end
    append!(rr,res)
  end
  rr
end

function Base.show(io::IO, ::MIME"text/html", a::CycPol)
  print(io, "\$")
  show(IOContext(io,:TeX=>true),a)
  print(io, "\$")
end

function Base.show(io::IO, ::MIME"text/plain", a::CycPol)
  if isempty(a.v) && !haskey(io,:typeinfo) print(io,typeof(a),": ") end
  show(io,a)
end

function Base.show(io::IO,a::CycPol)
  if !(get(io,:limit,false) || get(io,:TeX,false))
    print(io,"CycPol(",a.coeff,",",a.valuation)
#   for (r,m) in a.v print(io,",",r.r,"=>",m) end
    for e in decompose(a.v.d) 
      print(io,",(",e.conductor,",",e.no)
      if e.mul!=1 print(io,",",e.mul) end
      print(io,")")
    end
    print(io,")")
    return
  end
  if iszero(a.valuation) && isempty(a.v) 
    if isone(denominator(a.coeff)) print(io,numerator(a.coeff))
    else print(io,a.coeff)
    end
    return
  end
  print(io,format_coefficient(repr(numerator(a.coeff); context=io))) 
  v=LaurentPolynomials.varname[]
  if a.valuation==1 print(io,v)
  elseif a.valuation!=0 print(io,v,stringexp(io,a.valuation)) end
  for e in decompose(a.v.d)
#   println(e)
    if e.no>=0  
      if get(io,:expand,false)
        print(io,"(",prod(i->Pol()-E(e.conductor,i),dec(e.conductor)[e.no+1]),")")
      else print(io,get(io,:TeX,false) ? "\\Phi" : "Φ")
        print(io,stringprime(io,e.no))
        print(io,stringind(io,e.conductor))
      end
    else print(io,"(",v,"-",E(e[1],-e.no),")")
    end
    if e.mul!=1 print(io,stringexp(io,e.mul)) end
  end
  den=denominator(a.coeff)
  if !isone(den) print(io,"/",bracket_if_needed(repr(den;context=io))) end
end

# fields to test first: all n such that totient(n)<=12 except 11,13,22,26
const tested=[1,2,4,3,6,8,12,5,10,9,18,7,14,24,16,20,15,30,36,28,21,42]

# list of i such that φᵢ/φ_(i∩ conductor))≤d, so a polynomial of
# degree ≤d with coeffs in Q(ζ_conductor) could have roots power of ζᵢ
function bounds(conductor::Int,d::Int)::Vector{Int}
  if d==0 return Int[] end
  t=Vector{Int}[];t1=Vector{Int}[]
  for (p,m) in eachfactor(conductor)
    tp=[1];pw=p;while pw<=d push!(tp,pw);pw*=p end # tp=={powers of p<=d}
    push!(t,tp)
    push!(t1,tp*p^m)
  end
  for p in setdiff(primes(d+1),keys(factor(conductor)))
    tp=[1,p-1];pw=p*(p-1)
    while pw<=d push!(tp,pw);pw*=p end
    push!(t,tp)
    push!(t1,p.^(eachindex(tp).-1))
  end
  produ=function(l,d)
    p=filter(x->l[1][x]<=d,eachindex(l[1]))
    if length(l)==1 return map(x->[x],p) end
    return reduce(vcat,map(i->vcat.([i],produ(l[2:end],div(d,l[1][i]))),p))
  end
  p=[prod(map((i,j)->j[i],l,t1)) for l in produ(t,d)]
  p=union(map(divisors,p)...)
  p=setdiff(p,tested)
  sort(p,by=x->x/length(divisors(x)))
end

# next function is twice the speed of p(Cyc(x))
(p::Pol)(x::Root1)=transpose(p.c)*x.^(p.v:degree(p))
  
"""
`CycPol(p::Pol)` Converts `Pol` `p` to `CycPol`
    
```julia-repl
julia> @Pol q;CycPol(3*q^3-3q)
3qΦ₁Φ₂
```
"""
function CycPol(p::Pol{T};trace=false)where T
 # lot of code to be as efficient as possible in all cases
  if iszero(p) return zero(CycPol{T})
  elseif length(p.c)==1 # p==ax^s
    return CycPol(p.c[1],valuation(p))
  elseif 2==count(!iszero,p.c) # p==ax^s+bx^t
    a=Root1(-p[begin]//p[end])
    if a===nothing return CycPol(Pol(coefficients(p)),valuation(p)) end
    d=degree(p)-valuation(p)
    vcyc=(Root1(;r=(a.r+i)//d)=>1 for i in 0:d-1)
    return CycPol(p[end],valuation(p),ModuleElt(vcyc))
  end
  val=valuation(p)
  p=Pol(coefficients(p))
  coeff=p[end]
  p=coeff^2==1 ? p*coeff : p//coeff
  if denominator(p)==1 p=numerator(p) end
  vcyc=Pair{Root1,Int}[]

  # find factors Phi_i
  testcyc=function(c)
    if trace print("C$c") end
    found=false
    while true
      np,r=pseudodiv(p,cyclotomic_polynomial(c))
      if iszero(r) 
        append!(vcyc,[E(c,j)=>1 for j in (c==1 ? [0] : prime_residues(c))])
        p=(np[begin] isa Cyc) ? np : Pol(coefficients(np)) # why?
        if trace print("(d°$(degree(p))c$(conductor(p.c)))") end
        found=true
      else break
      end
    end
    return found
  end

  # find other primitive i-th roots of unity
  testall=function(i)
    if eltype(coefficients(p))<:Integer return false end # cannot have partial product
    found=false
    to_test=prime_residues(i)
    while true 
      to_test=filter(r->iszero(p(E(i,r))),to_test)
      if isempty(to_test) return found end
      found=true
      p=exactdiv(p,prod(r->Pol([-E(i,r),1],0),to_test))
      append!(vcyc,E.(i,to_test).=>1)
      if trace print("[d°$(degree(p))c$(conductor(p.c)),$i:",join(to_test,","),"]") end
      if degree(p)<div(totient(i),totient(gcd(i,conductor(p.c)))) return found end
    end
  end

  # first try commonly occuring fields
  for i in tested 
    if degree(p)>=totient(i) testcyc(i) end
    if degree(p)>0 testall(i) end
    if degree(p)==0 return CycPol(coeff,val,ModuleElt(vcyc)) end
  end
  
  # if not finished do a general search.
  # p is in Q(zeta_conductor)[x] so can only have a root in mu_i for i below
  cond=(p.c[1] isa Cyc) ? conductor(p.c) : 1
  try_=bounds(cond,degree(p))
# println("try_=$try_\n")
  i=1
  while i<=length(try_) 
#   push!(tested,try_[i])
    found=if cond==1 # All factors are Phi_i
         testcyc(try_[i])
    else testall(try_[i])
    end
    if found  
      cond=(p.c[1] isa Cyc) ? conductor(p.c) : 1
      try_=bounds(cond,degree(p))
      i=1
#	  print("tested==",tested,"\n")
    else i+=1
    end
  end
# println("now p=$p val=$val")
  
  CycPol(degree(p)==0 ? coeff : p*coeff,val,ModuleElt(vcyc))
end

CycPol(p::Frac)=CycPol(numerator(p))//CycPol(denominator(p))

function (p::CycPol)(x)
  res=p.valuation<0 ? (x*1//1)^p.valuation : x^p.valuation
  l=decompose(p.v.d)
  for e in l
    if iszero(res) return res end
    if e.no==0 res*=(cyclotomic_polynomial(e.conductor)(x))^e.mul end
  end
  for e in l
    if iszero(res) return res end
    if e.no>0 
      res*=prod(x-E(e.conductor,j) for j in dec(e.conductor)[e.no+1])^e.mul 
    end
  end
  pp=one(x)
  co=0
  pp=one(x)
  for e in l
    if iszero(res) return res end
    if e.no<0 
      if co==e.conductor pp*=(x-E(e.conductor,-e.no))^e.mul 
      else co=e.conductor
        res*=pp
        pp=(x-E(e.conductor,-e.no))^e.mul 
      end
    end
  end
  res*=pp
  if p.coeff isa Pol res*p.coeff(x) 
  elseif p.coeff isa Frac res*p.coeff(x;Rational=true)
  else res*p.coeff end
end

"""
`subs(p::CycPol,v::Pol)`

a fast routine to compute `CycPol(p(v))` but works for only two types
of polynomials:

  - `v=Pol([e],1)` for `e` a `Root1`, that is the value at `qe` for `e=ζₙᵏ`
  - `v=Pol([1],n)` that is the value at `qⁿ`
"""
function subs(p::CycPol,v::Pol{Root1})
  if degree(v)!=1 || valuation(v)!=1 error(v," should be Pol([Root1],1)") end
  e=v[1]
  coeff=p.coeff*e^degree(p)
  if coeff isa Pol coeff=(coeff*e^-degree(coeff))(v) end
  re=inv(Root1(e))
  CycPol(coeff,valuation(p),ModuleElt(r*re=>m for (r,m) in p.v))
end

function subs(p::CycPol,v::Pol{Int})
  if degree(v)!=valuation(v) || coefficients(v)!=[1] error(v," should be Pol()^n") end
  n=valuation(v)
  if n==0 return CycPol(p(1),0) end
  n=Int(n)
  val=n*valuation(p)
  if p.coeff isa Pol coeff=p.coeff(v) else coeff=p.coeff end
  if n>0
    vcyc=vcat((map(i->Root1(;r=i)=>pow,((0:n-1).+r.r)/n) for (r,pow) in p.v)...)
  else
    val+=n*sum(values(p.v))
    coeff*=(-1)^sum(values(p.v))*Cyc(prod(r^p for (r,p) in p.v))
    vcyc=vcat((map(i->Root1(;r=i)=>pow,((0:-n-1).-r.r)/-n) for (r,pow) in p.v)...)
  end
  CycPol(coeff,val,ModuleElt(vcyc))
end
    
# 281st generic degree of G34
const p=CycPol(E(3)//6,19,(1,0,3),(2,0,6),(4,0,2),(5,0),(6,0,4),(7,0),(8,0),
               (9,1),(10,0),(14,0),(15,4),(18,0),(21,0),(24,1),(30,0),(42,0))
#=
julia> @btime u=CycPols.p(Pol()) # gap 1.25 ms
#1.8.5 254.158 μs (5626 allocations: 416.19 KiB)
julia> @btime CycPol(u) # gap 8.2ms
#1.8.5 4.749 ms (92895 allocations: 7.35 MiB)
julia> @btime u(1)  # gap 40μs
#1.8.5 24.669 μs (553 allocations: 41.91 KiB)
julia> @btime CycPols.p(1) # gap 142μs
#1.8.5 5.140 μs (101 allocations: 19.88 KiB)
=#

const p1=Pol([1,0,-1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1,0,1],0)

const p2=CycPol(-4E(3),-129,(3,0),(6,0),(8,1),(8,3),(9,2),(12,2),(16,2),
(16,-1),(16,-9),(18,2,2),(21,2),(27,2),(32,-7),(32,-15),(32,-23),(32,-31),(39,
-1),(39,-4),(39,-7),(39,-10),(39,-16),(39,-19),(39,-22),(39,-25),(39,-28),(39,
-31),(39,-34),(39,-37),(42,1),(48,-11),(48,-19),(48,-35),(48,-43),(60,-7),(60,
-19),(60,-31),(60,-43),(78,-5),(78,-11),(78,-17),(78,-23),(78,-29),(78,-35),
(78,-41),(78,-47),(78,-53),(78,-59),(78,-71),(78,-77),(80,-7),(80,-23),(80,
-31),(80,-39),(80,-47),(80,-63),(80,-71),(80,-79),(88,-3),(88,-19),(88,-27),
(88,-35),(88,-43),(88,-51),(88,-59),(88,-67),(88,-75),(88,-83),(90,-1),(90,
-7),(90,-13),(90,-19),(90,-31),(90,-37),(90,-43),(90,-49),(90,-61),(90,-67),
(90,-73),(90,-79),(96,-5),(96,-13),(96,-29),(96,-37),(96,-53),(96,-61),(96,
-77),(96,-85),(104,-1),(104,-9),(104,-17),(104,-25),(104,-33),(104,-41),(104,
-49),(104,-57),(104,-73),(104,-81),(104,-89),(104,-97),(144,-1),(144,-17),
(144,-25),(144,-41),(144,-49),(144,-65),(144,-73),(144,-89),(144,-97),(144,
-113),(144,-121),(144,-137),(152,-5),(152,-13),(152,-21),(152,-29),(152,-37),
(152,-45),(152,-53),(152,-61),(152,-69),(152,-77),(152,-85),(152,-93),(152,
-101),(152,-109),(152,-117),(152,-125),(152,-141),(152,-149),(204,-11),(204,
-23),(204,-35),(204,-47),(204,-59),(204,-71),(204,-83),(204,-95),(204,-107),
(204,-131),(204,-143),(204,-155),(204,-167),(204,-179),(204,-191),(204,-203))
#=
julia> @btime u=CycPols.p2(Pol()) # gap 9.6ms
#1.8.5 14.701 ms (127154 allocations: 16.82 MiB)
julia> @btime CycPol(u) # gap 1.33s
#1.8.5 1.295 s (9661765 allocations: 1.83 GiB)
julia> @btime u(1)  # gap  43μs
#1.8.5 88.881 μs (905 allocations: 118.00 KiB)
julia> @btime CycPols.p2(1) # gap 1.1ms
#1.8.5 6.840 ms (3514 allocations: 7.13 MiB)
=#
end
