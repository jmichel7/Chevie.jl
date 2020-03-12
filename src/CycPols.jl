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
julia> Pol(:q)
Pol{Int64}: q

julia> p=CycPol(q^25-q^24-2q^23-q^2+q+2)
(q-2)Φ₁Φ₂Φ₂₃

julia> p(q) # a CycPol is a callable object, this call evaluates p at q
Pol{Cyc{Int64}}: q²⁵-q²⁴-2q²³-q²+q+2

julia> p*inv(CycPol(q^2+q+1))
(q-2)Φ₁Φ₂Φ₃⁻¹Φ₂₃

```
The variable name in a `CycPol` is set by default to the same as for `Pols`.

`CycPol`s are internally a `struct` with fields:

`.coeff`:  a coefficient, usually a cyclotomic number or a polynomial.

`.valuation`: an `Int`.

`.v`: a list of pairs `r=>m` of a root of unity `r` and a multiplicity `m`.
Here `r` is a `Root1`, internally a fraction `n//e` with `n<e` representing
`E(r)=E(e,n)`.

So `CycPol(c,val,v)` represents `c*q^val*prod((q-E(r))^m for (r,m) in v)`.

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

export CycPol
# to use as a stand-alone module uncomment the next line
# export roots, degree

using ..ModuleElts: ModuleElt, norm!
using ..Cycs: Root1, E, conductor, Cyc
using ..Pols
using ..Util: fromTeX, prime_residues, primitiveroot, phi
using ..Gapjm
import Primes

struct CycPol{T}
  coeff::T
  valuation::Int
  v::ModuleElt{Root1,Int}
end

CycPol(c,val::Int,v::Pair{Rational{Int},Int}...)=CycPol(c,val,
  ModuleElt(Pair{Root1,Int}[Root1(;r=r)=>m for (r,m) in v]))

function roots(c::CycPol)
  function f(e,m)
    if m<0 error("should be a true polynomial") end
    fill(e,m)
  end
  vcat([f(e,m) for (e,m) in c.v]...)
end

# 281st generic degree of G34
const p=CycPol((1//6)E(3),19,0//1=>3, 1//2=>6, 1//4=>2, 3//4=>2,
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

const p1=Pol([1,0,-1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1,0,1],0)
#=
  benchmark: 
julia> @btime u=CycPols.p(q) # gap 2ms
  868.139 μs (14931 allocations: 1.03 MiB)
    1.827 ms (90128 allocations: 4.62 MiB)
julia> @btime CycPol(u) # gap 12ms
  44.270 ms (500873 allocations: 56.74 MiB)
  31.011 ms (666125 allocations: 41.54 MiB)
julia> @btime u(1)  # gap 57μs
  107.098 μs (2559 allocations: 196.03 KiB)
   69.927 μs (2819 allocations: 149.95 KiB)
=#

Base.one(::Type{CycPol})=CycPol(1,0)
Base.isone(p::CycPol)=isone(p.coeff) && iszero(p.valuation) && iszero(p.v)
Base.zero(::Type{CycPol{T}}) where T=CycPol(zero(T),0)
Base.zero(::Type{CycPol})=zero(CycPol{Int})
Base.zero(a::CycPol)=CycPol(zero(a.coeff),0)

Gapjm.degree(a::CycPol)=sum(last,a.v)+a.valuation+degree(a.coeff)
Gapjm.valuation(a::CycPol)=a.valuation
Gapjm.valuation(a::CycPol,d::Root1)=reduce(+,c for (r,c) in a.v if r==d;init=0)
Gapjm.valuation(a::CycPol,d::Integer)=valuation(a,Root1(1,d))
Gapjm.valuation(a::CycPol,d::Rational)=valuation(a,Root1(;r=d))

function Base.:*(a::CycPol,b::CycPol)
  CycPol(a.coeff*b.coeff,a.valuation+b.valuation,a.v+b.v)
end
Base.:*(a::CycPol,b::Number)=iszero(b) ? zero(a) : CycPol(a.coeff*b,a.valuation,a.v)
Base.:*(b::Number,a::CycPol)=a*b

Base.inv(a::CycPol)=CycPol(a.coeff^2==1 ? a.coeff : inv(a.coeff), -a.valuation,
                                         -a.v)

Base.://(a::CycPol,b::CycPol)=a*inv(b)
Base.://(a::CycPol,b::Number)=CycPol(a.coeff//b,a.valuation,a.v)
Base.:div(a::CycPol,b::Number)=CycPol(div(a.coeff,b),a.valuation,a.v)

const dec_dict=Dict(1=>[[1]],2=>[[1]],
  8=>[[1,3,5,7],[1,5],[3,7],[1,7],[3,5],[1,3],[5,7]],
 12=>[[1,5,7,11],[1,5],[7,11],[1,7],[5,11],[1,11],[5,7]],
 15=>[[1,2,4,7,8,11,13,14],[1,4,11,14],[2,7,8,13],[1,4,7,13],[2,8,11,14],
      [1,4],[7,13],[11,14],[2,8]],
 16=>[[1,3,5,7,9,11,13,15],[1,7,9,15],[3,5,11,13],[1,5,9,13],[3,7,11,15]],
 20=>[[1,3,7,9,11,13,17,19],[1,9,11,19],[3,7,13,17],[1,9,13,17],[3,7,11,19]],
 21=>[[1,2,4,5,8,10,11,13,16,17,19,20],[1,4,10,13,16,19],[2,5,8,11,17,20]],
 24=>[[1,5,7,11,13,17,19,23],[1,7,13,19],[5,11,17,23],[1,7,17,23],[5,11,13,19],
      [1,5,19,23],[7,11,13,17],[1,11,17,19],[5,7,13,23],[1,5,13,17],
      [7,11,19,23],[7,13],[1,19],[5,23],[11,17]],
 30=>[[1,7,11,13,17,19,23,29],[1,11,19,29],[7,13,17,23],[1,7,13,19],
      [11,17,23,29],[11,29],[17,23],[1,19],[7,13]],
 36=>[[1,5,7,11,13,17,19,23,25,29,31,35],[1,7,13,19,25,31],[5,11,17,23,29,35],
      [7,11,19,23,31,35],[1,5,13,17,25,29]],
42=>[[1,5,11,13,17,19,23,25,29,31,37,41],[1,13,19,25,31,37],[5,11,17,23,29,41]])


# returns list of subsets of primitive_roots(d) wich have a `name` Φ^i
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

function show_factors(d)
  for i in eachindex(CycPols.dec(d))
    p=CycPol(;cond=d,no=i)
    println(IOContext(stdout,:limit=>true),p,"=",p(Pol(:q)))
  end
end
pr()=for d in sort(collect(keys(dec_dict))) show_factors(d) end

CycPol(;cond=1,no=1)=CycPol(1,0,map(i->i//cond=>1,dec(cond)[no])...)
  
function segment(v::ModuleElt{Root1,Int})
  res=typeof(v.d)[]
  n=empty(v.d)
  c=0
  for p in v.d
    c1=conductor(p[1])
    if c!=c1 
      if !iszero(c) push!(res,n) end
      c=c1
      n=empty(v.d)
    end
    push!(n,p)
  end
  if !iszero(c) push!(res,n) end
  res
end

# decompose the .v of a CycPol in subsets Φ^i (for printing)
function decompose(v::ModuleElt{Root1,Int})
  rr=Pair{NamedTuple{(:cond, :no),Tuple{Int,Int}},Int}[]
  for t in segment(v)
    c=conductor(t[1][1])
    t=[(exponent(e),p) for (e,p) in t]
    if c==1 
      push!(rr,(cond=c,no=1)=>t[1][2])
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
      if n!=0 push!(res,(cond=c,no=i)=>n) end
    end
    for i in 1:c  
      if v[i]!=0 push!(res,(cond=c,no=-i)=>v[i]) end 
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
  nov=iszero(a.valuation) && isempty(a.v)
  if a.coeff!=1 || nov
    c=a.coeff
    if c isa Rational && isone(denominator(c)) c=numerator(c) end
    s=sprint(show,c; context=io)
    if s=="-1" && !nov s="-" end
    if occursin(r"[+\-*/]",s[nextind(s,1):end]) && !nov s="($s)" end
    print(io,s) 
  end
  if a.valuation==1 print(io,"q")
  elseif a.valuation!=0 print(io,fromTeX(io,"q^{$(a.valuation)}")) end
  for (e,pow) in decompose(a.v)
#   println(e)
    if e.no>0  print(io,fromTeX(io,"\\Phi"*"'"^(e.no-1)*"_{$(e.cond)}"))
    else print(io,"(",Pol([-E(e[1])^-e.no,1],0),")")
    end
    if pow!=1 print(io,fromTeX(io,"^{$pow}")) end
  end
end

# fields to test first: all n such that phi(n)<=12 except 11,13,22,26
const tested=[1,2,4,3,6,8,12,5,10,9,18,24,16,20,7,14,15,30,36,28,21,42]

# list of i such that phi_i/phi(gcd(i,conductor))<=d
# so a pol of degree <=d with coeffs in Q(ζ_conductor) could have roots ζ_i^h
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

const trace=false
"""
`CycPol(p::Pol)`
    
Converts a polynomial to `CycPol`
    
```julia-repl
julia> CycPol(3*q^3-3)
3Φ₁Φ₃
```

Special code makes the conversion fast if `p` has not more than two nonzero
coefficients.
"""
function CycPol(p::Pol{T})where T
 # lot of code to be as efficient as possible in all cases
  if iszero(p) return zero(CycPol{T})
  elseif length(p.c)==1 # p==ax^s
    return CycPol(p.c[1],Pols.valuation(p))
  elseif 2==count(!iszero,p.c) # p==ax^s+bx^t
    a=Root1(-p.c[1]//p.c[end])
    if a===nothing return CycPol(Pol(p.c,0),Pols.valuation(p)) end
    d=length(p.c)-1
    vcyc=[Root1(;r=u)=>1 for u in (a.r .+(0:d-1))//d]
    return CycPol(p.c[end],Pols.valuation(p),ModuleElt(sort(vcyc)))
  end
  val=Pols.valuation(p)
  p=Pol(p.c,0)
  vcyc=zero(ModuleElt{Root1,Int})

  # find factors Phi_i
  testcyc=function(c)
    if trace print("C$c") end
    found=false
    while true
      np,r=Pols.divrem1(p,cyclotomic_polynomial(c))
      if iszero(r) 
        append!(vcyc,[Root1(j,c)=>1 for j in (c==1 ? [1] : prime_residues(c))])
        p=(np.c[1] isa Cyc) ? np : Pol(np.c,0)
        if trace print("(d°$(degree(p)) c$(conductor(p.c)))") end
        found=true
      else break
      end
    end
    return found
  end

  # find other primitive i-th roots of unity
  testall=function(i)
    found=false
    while true 
      l=filter(r->iszero(p(E(i,r))),prime_residues(i))
      if isempty(l) return found end
      found=true
      p=Pols.divrem1(p,prod(r->Pol([-E(i,r),1],0),l))[1]
      append!(vcyc,Root1.(l,i).=>1)
      if trace print("(d°$(degree(p)) c$(conductor(p.c)) e$i.$l)") end
      if degree(p)<div(phi(i),phi(gcd(i,conductor(p.c)))) return found end
    end
  end

  # first try commonly occuring fields
  for i in tested 
    if degree(p)>=phi(i) testcyc(i) end
    testall(i)
    if degree(p)==0 return CycPol(p.c[end],val,norm!(vcyc)) end
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
  
  CycPol(degree(p)==0 ? p.c[1] : p,val,norm!(vcyc))
end

function (p::CycPol)(x)
  res=x^p.valuation
  if !isempty(p.v)
    for v in segment(p.v)
      pp=prod((x-E(r))^m for (r,m) in v)
#     if all(isone,conductor.(pp.c)) pp=Pol(convert.(Int,pp.c),pp.v) end
      res*=pp
      if iszero(res) return res end
    end
  end
  res*p.coeff
end

end
