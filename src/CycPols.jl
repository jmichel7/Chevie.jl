"""
Cyclotomic  numbers, and cyclotomic polynomials  over the rationals or some
cyclotomic field, are important in reductive groups or Spetses. This module
deals  with them: the type `CycPol`  represents the product of a polynomial
with  a rational fraction in one variable with all poles or zeroes equal to
0  or  roots  of  unity.  The  advantages  of representing as `CycPol` such
objects    are:   nice   display   (factorized),   less   storage,   faster
multiplication,  division and evaluation. The drawback is that addition and
subtraction are not implemented!

```julia-repl
julia> @Pol q
Pol{Int64}: q

julia> p=CycPol(q^25-q^24-2q^23-q^2+q+2)
(q-2)Φ₁Φ₂Φ₂₃

julia> p(q) # a CycPol is a callable object, this call evaluates p at q
Pol{Int64}: q²⁵-q²⁴-2q²³-q²+q+2

julia> p*inv(CycPol(q^2+q+1))
(q-2)Φ₁Φ₂Φ₃⁻¹Φ₂₃

```
The variable name in a `CycPol` is set by default to the same as for `LaurentPolynomials`.

`CycPol`s are internally a `struct` with fields:

`.coeff`:  a coefficient, usually a cyclotomic number or a polynomial.

`.valuation`: an `Int`.

`.v`: a list of pairs `r=>m` of a root of unity `r` and a multiplicity `m`.
Here `r` is a `Root1`, internally a fraction `n//e` with `n<e` representing
`Cyc(r)=E(e,n)`.

So `CycPol(c,val,v)` represents `c*q^val*prod((q-Cyc(r))^m for (r,m) in v)`.

When   showing,  some  factors  of   the  cyclotomic  polynomial  `Φₙ`  are
represented.  If `n` has a primitive root  `ξ`, `ϕ′ₙ` is the product of the
`(q-ζ)` where `ζ` runs over the odd powers of `ξ`, and `ϕ″ₙ` is the product
for the even powers. The function `show_factors` gives the complete list of
recognized factors:

```julia-rep1
julia> CycPols.show_factors(24)
Φ₂₄=q⁸-q⁴+1
Φ′₂₄=q⁴+ζ₃²
Φ″₂₄=q⁴+ζ₃
Φ‴₂₄=q⁴-√2q³+q²-√2q+1
Φ⁗₂₄=q⁴+√2q³+q²+√2q+1
Φ⁽⁵⁾₂₄=q⁴-√6q³+3q²-√6q+1
Φ⁽⁶⁾₂₄=q⁴+√6q³+3q²+√6q+1
Φ⁽⁷⁾₂₄=q⁴+(√-2)q³-q²+(-√-2)q+1
Φ⁽⁸⁾₂₄=q⁴+(-√-2)q³-q²+(√-2)q+1
Φ⁽⁹⁾₂₄=q⁴-ζ₄q²-1
Φ⁽¹⁰⁾₂₄=q⁴+ζ₄q²-1
Φ⁽¹¹⁾₂₄=q²+(ζ₂₄+ζ₂₄¹⁹)q-ζ₃
Φ⁽¹²⁾₂₄=q²+(-ζ₂₄-ζ₂₄¹⁹)q-ζ₃
Φ⁽¹³⁾₂₄=q²+(ζ₂₄¹¹+ζ₂₄¹⁷)q-ζ₃²
Φ⁽¹⁴⁾₂₄=q²+(-ζ₂₄¹¹-ζ₂₄¹⁷)q-ζ₃²
```
"""
module CycPols

export CycPol,descent_of_scalars,ennola_twist, cyclotomic_polynomial,eigmat
# to use as a stand-alone module uncomment the next line
# export roots
import ..Gapjm: roots, gap

using ModuleElts: ModuleElt
using LaurentPolynomials: Pol, LaurentPolynomials, degree, valuation
using PuiseuxPolynomials: Mvp
using ..Cyclotomics: Root1, E, conductor, Cyc
using ..GLinearAlgebra: charpoly
using ..Combinat: collectby
using ..Util: prime_residues, primitiveroot, phi, divisors, factor, 
              stringexp, stringprime, format_coefficient, xprintln, stringind
using ..PermRoot: improve_type
import Primes: Primes, primes

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
function cyclotomic_polynomial(n::Integer)::Pol{Int}
  get!(cyclotomic_polynomial_dict,n) do
    res=Pol(fill(1,n),0;check=false)
    for d in divisors(n)
      if d!=1 && d!=n
        res,_=LaurentPolynomials.pseudodiv(res,cyclotomic_polynomial(d))
      end
    end
    res
  end
end

function lindivrem(a::Pol,z::Root1) # divrem(a,Pol()-z) where z root of 1
  res=Vector{promote_type(eltype(a.c),Cyc{Int})}(a.c[2:end])
  for i in length(res):-1:2
    res[i-1]+=res[i]*z
  end
  if iszero(res[1]*z+a.c[1]) return Pol(res,a.v;check=false) end
  # else return nothing since remainder is not zero
end

struct CycPol{T}
  coeff::T
  valuation::Int
  v::ModuleElt{Root1,Int}
end

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
CycPol(c,val::Int,v::Pair{Rational{Int},Int}...)=CycPol(c,val,
  ModuleElt(Pair{Root1,Int}[Root1(;r=r)=>m for (r,m) in v];check=false)) 

# all the Root1 roots of c
function roots(c::CycPol)
  function f(e,m)
    if m<0 error("should be a true polynomial") end
    fill(e,m)
  end
  if isempty(c.v) return Root1[] end
  vcat([f(e,m) for (e,m) in c.v]...)
end

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
LaurentPolynomials.valuation(a::CycPol,d::Integer)=valuation(a,E(d))
LaurentPolynomials.valuation(a::CycPol,d::Rational)=valuation(a,Root1(;r=d))

function Base.:*(a::CycPol,b::CycPol)
  if iszero(a) || iszero(b) return zero(a) end
  CycPol(a.coeff*b.coeff,a.valuation+b.valuation,a.v+b.v)
end

Base.conj(a::CycPol)=CycPol(conj(a.coeff),a.valuation,ModuleElt(
           inv(r)=>m for (r,m) in a.v))
Base.:*(a::CycPol,b::Number)=iszero(b) ? zero(a) : CycPol(a.coeff*b,a.valuation,a.v)
Base.:*(b::Number,a::CycPol)=a*b
Base.:-(a::CycPol)=CycPol(-a.coeff,a.valuation,a.v)
Base.inv(a::CycPol)=CycPol(LaurentPolynomials.bestinv(a.coeff), -a.valuation, -a.v)
Base.:^(a::CycPol, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                      Base.power_by_squaring(inv(a),-n)
Base.://(a::CycPol,b::CycPol)=a*inv(b)
Base.://(a::CycPol,b::Number)=CycPol(a.coeff//b,a.valuation,a.v)
Base.:div(a::CycPol,b::Number)=CycPol(div(a.coeff,b),a.valuation,a.v)

Base.lcm(a::CycPol,b::CycPol)=CycPol(1,max(a.valuation,b.valuation),
                                       ModuleElts.merge2(max,a.v,b.v))
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

# returns list of subsets of primitive_roots(d) wich have a `name` Φ_d^i
function dec(d::Int)
  get!(dec_dict,d) do
    dd=[prime_residues(d)]
    if (r=primitiveroot(d))!==nothing
      for a in 0:1
        push!(dd,sort(powermod.(r,(0:2:phi(d)-2).+a,d)))
      end
    end
    dd
  end
end

CycPol(;cond=1,no=1)=CycPol(1,0,map(i->i//cond=>1,dec(cond)[no])...)
  
function show_factors(d)
  for i in eachindex(CycPols.dec(d))
    p=CycPol(;cond=d,no=i)
    xprintln(p,"=",p(Pol()))
  end
end

pr()=for d in sort(collect(keys(dec_dict))) show_factors(d) end

function segment(v::Vector{Pair{Root1,Int}})
  res=typeof(v)[]
  n=empty(v)
  c=0
  for p in v
    c1=denominator(p[1])
    if c!=c1 
      if !iszero(c) push!(res,n) end
      c=c1
      n=empty(v)
    end
    push!(n,p)
  end
  if !iszero(c) push!(res,n) end
  res
end

# decompose the .v of a CycPol in subsets Φ^i (for printing)
function decompose(v::Vector{Pair{Root1,Int}})
  rr=Pair{NamedTuple{(:conductor, :no),Tuple{Int,Int}},Int}[]
  for t in segment(v)
    c=denominator(t[1][1])
    t=[(numerator(e),p) for (e,p) in t]
    if c==1 
      push!(rr,(conductor=c,no=1)=>t[1][2])
      continue
    end
    res=[]
    v=fill(0,c)
    v[first.(t)]=last.(t)
    dd=dec(c)
    for (i,r) in enumerate(dd)
      if (n=minimum(v[r]))>0 v[r].-=n
      elseif (n=maximum(v[r]))<0 v[r].-=n 
      else n=0 
      end
      if n!=0 push!(res,(conductor=c,no=i)=>n) end
    end
    for i in 1:c  
      if v[i]!=0 push!(res,(conductor=c,no=-i)=>v[i]) end 
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

function Base.show(io::IO,a::CycPol)
  if !(get(io,:limit,false) || get(io,:TeX,false))
    print(io,"CycPol(",a.coeff,",",a.valuation)
    for (r,m) in a.v print(io,",",r.r,"=>",m) end
    print(io,")")
    return
  end
  den=denominator(a.coeff)
  c=improve_type(a.coeff*den)
  s=repr(c; context=IOContext(io,:varname=>:q))
  if iszero(a.valuation) && isempty(a.v) print(io,s)
  else
    s=format_coefficient(s)
    print(io,s) 
    if a.valuation==1 print(io,"q")
    elseif a.valuation!=0 print(io,"q",stringexp(io,a.valuation)) end
    for (e,pow) in decompose(a.v.d)
  #   println(e)
      if e.no>0  print(io,get(io,:TeX,false) ? "\\Phi" : "Φ")
        print(io,stringprime(io,e.no-1))
        print(io,stringind(io,e.conductor))
      else print(io,"(",Pol([-E(e[1],-e.no),1],0),")")
      end
      if pow!=1 print(io,stringexp(io,pow)) end
    end
  end
  if !isone(den) print(io,"/",den) end
end

# fields to test first: all n such that phi(n)<=12 except 11,13,22,26
const tested=[1,2,4,3,6,8,12,5,10,9,18,24,16,20,7,14,15,30,36,28,21,42]

# list of i such that φᵢ/φ_(i∩ conductor))≤d, so a polynomial of
# degree ≤d with coeffs in Q(ζ_conductor) could have roots power of ζᵢ
function bounds(conductor::Int,d::Int)::Vector{Int}
  if d==0 return Int[] end
  f=factor(conductor)
  t=Vector{Int}[];t1=Vector{Int}[]
  local p
  for (p,m) in f
   tp=[1];pw=p;while pw<=d push!(tp,pw);pw*=p end # tp=={powers of p<=d}
    push!(t,tp)
    push!(t1,tp*p^m)
  end
  for p in setdiff(Primes.primes(d+1),keys(f))
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
(p::Pol)(x::Root1)=sum(map(*,p.c,x.^(p.v:degree(p))))
  
"""
`CycPol(p::Pol)`
    
Converts a polynomial to `CycPol`
    
```julia-repl
julia> CycPol(3*q^3-3)
3Φ₁Φ₃
```
"""
function CycPol(p::Pol{T};trace=false)where T
 # lot of code to be as efficient as possible in all cases
  if iszero(p) return zero(CycPol{T})
  elseif length(p.c)==1 # p==ax^s
    return CycPol(p.c[1],valuation(p))
  elseif 2==count(!iszero,p.c) # p==ax^s+bx^t
    a=Root1(-p.c[1]//p.c[end])
    if a===nothing return CycPol(Pol(p.c,0),valuation(p)) end
    d=length(p.c)-1
    vcyc=(Root1(;r=(a.r+i)//d)=>1 for i in 0:d-1)
    return CycPol(p.c[end],valuation(p),ModuleElt(vcyc))
  end
  val=valuation(p)
  coeff=p.c[end]
  p=improve_type(coeff^2==1 ? Pol(p.c.*coeff,0) : Pol(p.c.//coeff,0))
  vcyc=Pair{Root1,Int}[]

  # find factors Phi_i
  testcyc=function(c)
    if trace print("C$c") end
    found=false
    while true
      np,r=LaurentPolynomials.pseudodiv(p,cyclotomic_polynomial(c))
      if iszero(r) 
        append!(vcyc,[E(c,j)=>1 for j in (c==1 ? [0] : prime_residues(c))])
        p=(np.c[1] isa Cyc) ? np : Pol(np.c,0)
        if trace print("(d°$(degree(p))c$(conductor(p.c)))") end
        found=true
      else break
      end
    end
    return found
  end

  # find other primitive i-th roots of unity
  testall=function(i)
    if eltype(p.c)<:Integer return false end # cannot have partial product
    found=false
    to_test=prime_residues(i)
    while true 
      to_test=filter(r->iszero(p(E(i,r))),to_test)
      if isempty(to_test) return found end
      found=true
      p=LaurentPolynomials.pseudodiv(p,prod(r->Pol([-E(i,r),1],0),to_test))[1]
      append!(vcyc,E.(i,to_test).=>1)
      if trace print("[d°$(degree(p))c$(conductor(p.c))]") end
      if degree(p)<div(phi(i),phi(gcd(i,conductor(p.c)))) return found end
    end
#   for r in to_test
#     while true 
#       p1=lindivrem(p,E(i,r))
#       if !isnothing(p1)
#         p=p1
#         found=true
#         push!(vcyc,E(i,r)=>1)
#         if trace print("(d°$(degree(p)) c$(conductor(p.c)) e$i.$r)") end
#         if degree(p)<div(phi(i),phi(gcd(i,conductor(p.c)))) return found end
#       else break
#       end
#     end
#   end
  end

  # first try commonly occuring fields
  for i in tested 
    if degree(p)>=phi(i) testcyc(i) end
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
  
  CycPol(degree(p)==0 ? coeff : improve_type(p*coeff),val,ModuleElt(vcyc))
end

function CycPol(x::Mvp)
  if !isinteger(degree(x)) CycPol(x,0)
  else CycPol(Pol(x))
  end
end

function (p::CycPol)(x)
  res=x^p.valuation
  if !isempty(p.v)
    v=decompose(p.v.d)
    for ((cond,no),mul) in filter(x->x[1].no==1,v)
      res*=(cyclotomic_polynomial(cond)(x))^mul
    end
    for ((cond,no),mul) in filter(x->x[1].no>1,v)
      res*=prod(x-E(cond,j) for j in dec(cond)[no])^mul
    end
    for v in collectby(x->x[1].conductor,filter(x->x[1].no<0,v))
      res*=prod((x-E(cond,-no))^mul for ((cond,no),mul) in v)
    end
  end
  if p.coeff isa Pol res*p.coeff(x) else res*p.coeff end
end

# Fast routine for CycPol(Value(p,q*e)) for a root of unity e
function ennola_twist(p::CycPol,e)
  if p.coeff isa Pol coeff=p.coeff(Pol()*e) else coeff=p.coeff end
  coeff=improve_type(p.coeff*e^degree(p))
  re=Root1(e)
  vcyc=[r*inv(re)=>m for (r,m) in p.v]
  CycPol(coeff,p.valuation,ModuleElt(vcyc))
end

# Fast routine for  CycPol(Value(p,q^n))
function descent_of_scalars(p::CycPol,n)
  if n==0 return CycPol(value(p,1)) end
  n=Int(n)
  valuation=n*p.valuation
  if n>0
    coeff=p.coeff
    vcyc=vcat((map(i->Root1(;r=i)=>pow,((0:n-1).+r.r)/n) for (r,pow) in p.v)...)
  else
    valuation+=n*sum(values(p.v))
    coeff=p.coeff*(-1)^sum(values(p.v))*Cyc(prod(r^p for (r,p) in p.v))
    vcyc=vcat((map(i->Root1(;r=i)=>pow,((0:-n-1).-r.r)/-n) for (r,pow) in p.v)...)
  end
  CycPol(coeff,valuation,ModuleElt(vcyc))
end
    
# export positive CycPols to GAP3
function gap(p::CycPol)
  if any(<(0),values(p.v)) error("non-positive") end
  res=string("[",gap(p.coeff),",",p.valuation,",")
  res*join(map(x->join(map(gap,fill(x[1].r,x[2])),","),pairs(p.v)),",")*"]"
end

" eigenvalues as Cyclotomics of a matrix of finite order"
eigmat(m)=roots(CycPol(Pol(charpoly(m))))

# 281st generic degree of G34; p(Pol()) 2ms back to CycPol 14ms (gap3 10ms)
const p=CycPol(E(3)//6,19,0//1=>3, 1//2=>6, 1//4=>2, 3//4=>2,
1//5=>1, 2//5=>1, 3//5=>1, 4//5=>1, 1//6=>4, 5//6=>4, 1//7=>1, 2//7=>1,
3//7=>1, 4//7=>1, 5//7=>1, 6//7=>1, 1//8=>1, 3//8=>1, 5//8=>1, 7//8=>1,
1//9=>1, 4//9=>1, 7//9=>1, 1//10=>1, 3//10=>1, 7//10=>1, 9//10=>1, 1//14=>1,
3//14=>1, 5//14=>1, 9//14=>1, 11//14=>1, 13//14=>1, 2//15=>1, 8//15=>1,
11//15=>1, 14//15=>1, 1//18=>1, 5//18=>1, 7//18=>1, 11//18=>1, 13//18=>1,
17//18=>1, 1//21=>1, 2//21=>1, 4//21=>1, 5//21=>1, 8//21=>1, 10//21=>1,
11//21=>1, 13//21=>1, 16//21=>1, 17//21=>1, 19//21=>1, 20//21=>1, 1//24=>1,
7//24=>1, 13//24=>1, 19//24=>1, 1//30=>1, 7//30=>1, 11//30=>1, 13//30=>1,
17//30=>1, 19//30=>1, 23//30=>1, 29//30=>1, 1//42=>1, 5//42=>1, 11//42=>1,
13//42=>1, 17//42=>1, 19//42=>1, 23//42=>1, 25//42=>1, 29//42=>1, 31//42=>1,
37//42=>1, 41//42=>1)

# a worse polynomial; p2(Pol()) 28ms (gap3 13ms) back to CycPol 1.8sec(gap3 1.8)
const p2=CycPol(-4E(3),-129,1//3=>1,2//3=>1,1//6=>1,5//6=>1,1//8=>2,5//8=>1,7//8=>1,
2//9=>1,5//9=>1,8//9=>1,7//12=>1,11//12=>1,1//16=>1,3//16=>1,5//16=>1,
9//16=>1,11//16=>1,13//16=>1,5//18=>2,11//18=>2,17//18=>2,2//21=>1,5//21=>1,
8//21=>1,11//21=>1,17//21=>1,20//21=>1,2//27=>1,5//27=>1,8//27=>1,11//27=>1,
14//27=>1,17//27=>1,20//27=>1,23//27=>1,26//27=>1,7//32=>1,15//32=>1,
23//32=>1,31//32=>1,1//39=>1,4//39=>1,7//39=>1,10//39=>1,16//39=>1,19//39=>1,
22//39=>1,25//39=>1,28//39=>1,31//39=>1,34//39=>1,37//39=>1,1//42=>1,
13//42=>1,19//42=>1,25//42=>1,31//42=>1,37//42=>1,11//48=>1,19//48=>1,
35//48=>1,43//48=>1,7//60=>1,19//60=>1,31//60=>1,43//60=>1,5//78=>1,11//78=>1,
17//78=>1,23//78=>1,29//78=>1,35//78=>1,41//78=>1,47//78=>1,53//78=>1,
59//78=>1,71//78=>1,77//78=>1,7//80=>1,23//80=>1,31//80=>1,39//80=>1,
47//80=>1,63//80=>1,71//80=>1,79//80=>1,3//88=>1,19//88=>1,27//88=>1,
35//88=>1,43//88=>1,51//88=>1,59//88=>1,67//88=>1,75//88=>1,83//88=>1,
1//90=>1,7//90=>1,13//90=>1,19//90=>1,31//90=>1,37//90=>1,43//90=>1,49//90=>1,
61//90=>1,67//90=>1,73//90=>1,79//90=>1,5//96=>1,13//96=>1,29//96=>1,
37//96=>1,53//96=>1,61//96=>1,77//96=>1,85//96=>1,1//104=>1,9//104=>1,
17//104=>1,25//104=>1,33//104=>1,41//104=>1,49//104=>1,57//104=>1,73//104=>1,
81//104=>1,89//104=>1,97//104=>1,1//144=>1,17//144=>1,25//144=>1,41//144=>1,
49//144=>1,65//144=>1,73//144=>1,89//144=>1,97//144=>1,113//144=>1,
121//144=>1,137//144=>1,5//152=>1,13//152=>1,21//152=>1,29//152=>1,37//152=>1,
45//152=>1,53//152=>1,61//152=>1,69//152=>1,77//152=>1,85//152=>1,93//152=>1,
101//152=>1,109//152=>1,117//152=>1,125//152=>1,141//152=>1,149//152=>1,
11//204=>1,23//204=>1,35//204=>1,47//204=>1,59//204=>1,71//204=>1,83//204=>1,
95//204=>1,107//204=>1,131//204=>1,143//204=>1,155//204=>1,167//204=>1,
179//204=>1,191//204=>1,203//204=>1)

const p1=Pol([1,0,-1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1,0,1],0)
#=
  benchmark: 
julia> @btime u=CycPols.p(Pol()) # gap 2ms
  868.139 μs (14931 allocations: 1.03 MiB)
    1.827 ms (90128 allocations: 4.62 MiB)
  763.055 μs (18850 allocations: 1.15 MiB)
julia> @btime CycPol(u) # gap 12ms
  44.270 ms (500873 allocations: 56.74 MiB)
  31.011 ms (666125 allocations: 41.54 MiB)
  14.205 ms (371324 allocations: 23.16 MiB)
  14.055 ms (237536 allocations: 17.07 MiB)
julia> @btime u(1)  # gap 57μs
  107.098 μs (2559 allocations: 196.03 KiB)
   69.927 μs (2819 allocations: 149.95 KiB)
   77.146 μs (1875 allocations: 126.95 KiB)
=#

end
