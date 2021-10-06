"""
Cyclotomic  numbers means complex numbers which are sums of rationals times
roots of unity.

They  are a very important feature of GAP, since character values of finite
groups are cyclotomics.

They  have a normal form given by writing them in the Zumbroich basis. This
form  allows to find  the smallest Cyclotomic  field which contains a given
number,   and  decide   in  particular   if  a   cyclotomic  is  zero.  Let
ζₙ=exp(2iπ/n).  The Zumbroich basis is a  particular subset of size φ(n) of
1,ζₙ,ζₙ²,…,ζₙⁿ⁻¹ which forms a basis of ℚ (ζₙ).

I  started  this  file  by  porting  Christian  Stump's Sage code, which is
simpler to understand than GAP's code. The reference for the algorithms is

T. Breuer, Integral bases for subfields of cyclotomic fields AAECC 8 (1997)

As  does  GAP,  I  lower  automatically  numbers  after  each  computation;
currently  the code about 50% slower than the C code in GAP since it is not
as  much optimized. GAP also  converts a Cyclotomic which  is rational to a
Rational,  a Rational which is integral to  an Int, a BigInt which is small
to  a small Int, etc... This is tremendously useful but needs a new type of
number to be added to Julia, which I am not competent enough to try.

The main way to build a Cyclotomic number is to use the function `E(n,k=1)`
which constructs ζₙᵏ.

# Examples
```julia-repl
julia> E(3)+E(4)
Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹

julia> E(3,2)
Cyc{Int64}: ζ₃²

julia> 1+E(3,2)
Cyc{Int64}: -ζ₃

julia> a=E(4)-E(4)
Cyc{Int64}: 0

julia> conductor(a) # a has been lowered to ℚ (ζ₁)=ℚ 
1

julia> typeof(convert(Int,a))
Int64

julia> convert(Int,E(4))
ERROR: InexactError: convert(Int64, E(4))

julia> inv(1+E(4)) # inverses often need Rational coefficients
Cyc{Rational{Int64}}: (1-ζ₄)/2

julia> inv(E(5)+E(5,4)) # but not always
Cyc{Int64}: -ζ₅²-ζ₅³

julia> Cyc(1//2+im) # one can convert Gaussian integers or rationals
Cyc{Rational{Int64}}: (1+2ζ₄)/2

julia> conj(1+E(4)) # complex conjugate
Cyc{Int64}: 1-ζ₄

julia> real(E(5))  # real part
Cyc{Rational{Int64}}: (-1+√5)/4

julia> imag(E(5))  # imaginary part
Cyc{Rational{Int64}}: (ζ₅-ζ₅⁴)/2

julia> c=E(9)   # an effect of the Zumbroich basis
Cyc{Int64}: -ζ₉⁴-ζ₉⁷

julia> Root1(c) # but you can decide whether a Cyc is a root of unity
Root1: ζ₉

julia> Root1(1+E(4)) # it returns nothing for a non-root

julia> Root1(4,1)
Root1: ζ₄

julia> c=Root1(;r=1//4)*Root1(3,1) # faster computation for roots of unity
Root1: ζ₁₂⁷

julia> E(c) # convert back to Cyc
Cyc{Int64}: ζ₁₂⁷

julia> c=Complex{Float64}(E(3))  # convert to float is sometimes useful
-0.4999999999999999 + 0.8660254037844387im
```

`Cyc`s have methods `copy, hash, ==, cmp, isless` (total order) so they can
be  keys in hashes or  elements of sets. Cyclotomics  which are integers or
rationals compare correctly to integers or rationals:

```julia-repl
julia> -1<Cyc(0)<1
true
```
For more information see the methods conductor, coefficients, denominator,
ER, Quadratic, galois, root. 

Finally, a benchmark:

```benchmark
julia> function testmat(p) 
         ss=[[i,j] for i in 0:p-1 for j in i+1:p-1]
         [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
       end
testmat (generic function with 1 method)

julia> @btime Cycs.testmat(12)^2;
  346.079 ms (4367402 allocations: 366.17 MiB)
```
The equivalent in GAP:

```
testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
end; 
```
testmat(12)^2 takes 0.35s in GAP3, 0.29s in GAP4
"""
module Cycs
#import Gapjm: coefficients, root
# to use as a stand-alone module comment above line and uncomment next
export coefficients, root, E, ER, Cyc, conductor, galois, Root1, Quadratic

using ..Util: fromTeX, printTeX, format_coefficient, factor, prime_residues, 
              phi, bracket_if_needed, xprint
using ..Combinat: constant

const use_list=false # I tried two different implementations. 
                     # The ModuleElt is twice the speed of the list one
if use_list
struct Cyc{T <: Real}<: Number   # a cyclotomic number
  n::Int       # conductor
  d::Vector{T} # the i-th element is the coefficient on zumbroich_basis[i]
end
else
using ..ModuleElts
const MM=ModuleElt # HModuleElt is twice slower
struct Cyc{T <: Real}<: Number   # a cyclotomic number
  n::Int              # conductor
  d::MM{Int,T} # list of pairs: i=>coeff on zumbroich_basis[i]
end
end

"""
   `conductor(c::Cyc)`
   `conductor(v::AbstractVector)`

returns the smallest positive integer  n  uch that `c∈ ℚ (ζₙ)` (resp. all
elements of `v` are in `ℚ (ζₙ)`).

```julia-repl
julia> conductor(E(9))
9

julia> conductor([E(3),1//2,E(4)])
12
```
"""
conductor(c::Cyc)=c.n
conductor(a::Array)=lcm(conductor.(a))
conductor(i::Integer)=1
conductor(i::Rational)=1

const zumbroich_basis_dict=Dict(1=>[0]) # The Zumbroich basis is memoized
"""
  zumbroich_basis(n::Int) 

  returns  the Zumbroich basis of  ℚ (ζₙ) as the  vector of i in 0:n-1 such
  that `ζₙⁱ` is in the basis
"""
function zumbroich_basis(n::Int)::Vector{Int}
  get!(zumbroich_basis_dict,n) do
  if n==1 return [0] end
  function J(k::Int, p::Int) # see [Breuer] Rem. 1 p. 283
    if k==0 if p==2 return 0:0 else return 1:p-1 end
    elseif p==2 return 0:1
    else return div(1-p,2):div(p-1,2)
    end
  end
  nfact=factor(n)
  res=
  let J=J
    [[div(n*i,p^k) for i in J(k-1,p)] for (p,np) in nfact for k in 1:np]
  end
  v=sum.(vec(collect(Iterators.product(res...))))
  sort(v.%n)
  end
end

"""
`coefficients(c::Cyc)`

for  a cyclotomic `c` of conductor `n`,  returns a vector `v` of length `n`
such that `c==∑ᵢ vᵢ₋₁ ζⁱ`.

```julia-repl
julia> coefficients(E(9))
9-element Vector{Int64}:
  0
  0
  0
  0
 -1
  0
  0
 -1
  0
```
"""
function coefficients(c::Cyc{T})where T
  res=zeros(T,conductor(c))
if use_list
  for (p,i) in enumerate(zumbroich_basis(length(res))) res[i+1]=c.d[p] end
else
  for (i,v) in c.d res[i+1]=v end
end
  res
end
  
"""
`denominator(c::Cyc{Rational})`

returns the smallest `d` such that `d*c` has integral coefficients (thus is
an algebraic integer).
"""
Base.denominator(c::Cyc)=lcm(denominator.(values(c.d)))

Base.numerator(c::Cyc{<:Rational{T}}) where T =Cyc{T}(c*denominator(c))

const Elist_dict=Dict((1,0)=>(true=>[0])) # to memoize Elist
"""
  Elist(n,i)  
  
  expresses  ζₙⁱ  in  zumbroich_basis(n):  it  is  a  sum  of some ζₙʲ with
  coefficients all 1 or all -1. The result is a Pair sgn=>inds where sgn is
  true  if coefficients are all 1 and false otherwise, and inds is the list
  of i in 0:n-1 such that ζₙⁱ occurs with a non-zero coefficient.
"""
function Elist(n::Int,i1::Int=1)
  i=mod(i1,n)
  get!(Elist_dict,(n,i)) do
    mp=Int[]
    j=i
    for (p,np) in factor(n)
      f=p^np
      m=div(n,f)
      cnt=mod(j*invmod(m,f),f)
      j-=cnt*m
      if p==2
        if 1==div(cnt,div(f,p)) push!(mp,p) end
      else
        tmp=zeros(Int,np)
        for k in 1:np
          f=div(f,p)
          tmp[k],cnt=divrem(cnt,f)
        end
        for k in np-1:-1:1
          if tmp[k+1]>div(p-1,2)
            tmp[k]+=1
            if k==1 && tmp[k]==p tmp[k]=0 end
          end
        end
        if tmp[1]==0 push!(mp,p) end
      end
    end
if use_list
    z=zumbroich_basis(n)
    if isempty(mp) return true=> Vector{Int}(indexin([i],z)) end
    v=vec(sum.(Iterators.product((div(n,p)*(1:p-1) for p in mp)...)))
    iseven(length(mp))=>Vector{Int}(indexin((i .+ v).%n,z))
else
    if isempty(mp) return true=> [i] end
    v=vec(sum.(Iterators.product((div(n,p)*(1:p-1) for p in mp)...)))
    iseven(length(mp))=>sort((i .+ v).%n)
end
  end
end

const E_dict=Dict((1,0)=>Cyc(1, use_list ? [1] : MM(0=>1)))
"""
  E(n::Integer,k::Integer=1) is exp(2i k π/n)
"""
function E(n1,i1=1)
  n=Int(n1)
  i=mod(Int(i1),n)
  get!(E_dict,(n,i)) do
    s,l=Elist(n,i) #::Pair{Bool,Vector{Int}}
if use_list
  v=zeros(Int,length(zumbroich_basis(n)))
  v[l].=ifelse(s,1,-1)
  lower(Cyc(n,v))
else
  lower(Cyc(n,MM(l.=>ifelse(s,1,-1);check=false)))
end
  end
end

E(;r)=E(denominator(r),numerator(r))

if use_list
Base.zero(c::Cyc)=Cyc(1,eltype(c.d)[0])
Base.zero(::Type{Cyc{T}}) where T=Cyc(1,T[0])
Base.iszero(c::Cyc)=c.n==1 && iszero(c.d[1])
else
Base.zero(c::Cyc)=Cyc(1,zero(c.d))
Base.zero(::Type{Cyc{T}}) where T=Cyc(1,zero(MM{Int,T}))
Base.iszero(c::Cyc)=iszero(c.d)
end
#Base.zero(m::Array{Cyc})=zero.(m)
Base.one(c::Cyc)=E(1,0)

function Cyc(c::Complex{T}) where T
  if iszero(imag(c)) return Cyc(real(c)) end
if use_list
  if iszero(real(c)) return Cyc(4,[0,imag(c)])
  else return Cyc(4,[real(c),imag(c)])
  end
else
  if iszero(real(c)) return Cyc(4,MM(1=>imag(c)))
  else return Cyc(4,MM(0=>real(c),1=>imag(c);check=false))
  end
end
end

if use_list
Cyc(i::Real)=Cyc(1,[i])
else
Cyc(i::Real)=iszero(i) ? zero(Cyc{typeof(i)}) : Cyc(1,MM(0=>i))
end
Cyc{T}(i::Real) where T<:Real=Cyc(T(i))

Cyc{T}(c::Complex{T1}) where {T,T1}=T(real(c))+E(4)*T(imag(c))

function Cyc{T}(c::Cyc{T1}) where {T,T1}
  if T==T1 return c end
if use_list
  Cyc(c.n,T.(c.d))
else
  Cyc(c.n,convert(MM{Int,T},c.d))
end
end

if use_list
  num(c::Cyc)=c.d[1]
else
  num(c::Cyc{T}) where T =iszero(c) ? zero(T) : last(first(c.d))
end

function (::Type{T})(c::Cyc)where T<:Union{Integer,Rational}
  if c.n==1 return T(num(c)) end
  throw(InexactError(:convert,T,c))
end

function Base.convert(::Type{T},c::Cyc;check=true)where T<:AbstractFloat
  if check && !isreal(c) throw(InexactError(:convert,T,c)) end
  real(convert(Complex{T},c))
end

function Complex{T}(c::Cyc)where T<:AbstractFloat
if use_list
  sum(x->x[2]*cispi(2*T(x[1])/c.n), zip(zumbroich_basis(c.n),c.d))
else
  iszero(c) ? Complex{T}(0.0) : sum(v*cispi(2*T(k)/c.n) for (k,v) in c.d)
end
end

function complexgaussian(c::Cyc{T})where T
if use_list
  c.d[1]+im*c.d[2]
else
  r=i=zero(T)
  for (e,c) in c.d
    if e==0 r=c
    else i=c
    end
  end
  Complex(r,i)
end
end
  
Base.complex(c::Cyc)=c.n==1 ? num(c) : c.n==4 ? complexgaussian(c) : 
                                                           Complex{Float64}(c)

function Complex{T}(c::Cyc)where T<:Union{Integer,Rational}
  if c.n==1 return Complex{T}(num(c)) end
  if c.n==4 return Complex{T}(complexgaussian(c)) end
  throw(InexactError(:convert,Complex{T},c))
end

Base.isinteger(c::Cyc)=c.n==1 && isinteger(num(c))

Base.isreal(c::Cyc)=c.n==1 || c==conj(c)

function Base.real(c::Cyc{T}) where T<:Real
  if c.n==1 return num(c) end
  (c+conj(c))/2
end

function Base.imag(c::Cyc{T}) where T<:Real
  if c.n==1 return 0 end
  (c-conj(c))/2
end

if false
# l is a list of pairs i=>c representing E(n,i)*c
function sumroots(n::Int,l)
  res=typeof(first(l))[]
# res=eltype(l)[]
  for (i,c) in l 
    (s,v)=Elist(n,i)
    if !s c=-c end
    for k in v push!(res,k=>c) end
  end
  Cyc(n,MM(res;check=length(l)>1))
end
else # 10% slower for small fields, faster for big ones
function sumroots(n::Int,l)
  res=fill(zero(last(first(l))),n)
  for (i,c) in l 
    (s,v)=Elist(n,i)
    if !s c=-c end
    for k in v res[k+1]+=c end
  end
  zb=zumbroich_basis(n)
  Cyc(n,MM(filter(x->!iszero(last(x)),map(i->i=>res[i+1],zb));check=false))
# Cyc(n,MM(i=>res[i+1] for i in zb if res[i+1]!=0;check=false))
end
end

function sumroot(res,n,i,c)
  (s,v)=Elist(n,i)
  if !s c=-c end
  for k in v push!(res,k=>c) end
end

function raise(n::Int,c::Cyc) # write c in Q(ζ_n) if c.n divides n
  if n==c.n return c end
  m=div(n,c.n)
if use_list
  res=zeros(eltype(c.d),length(zumbroich_basis(n)))
  z=zumbroich_basis(c.n)
  for (i,v) in enumerate(c.d)
    if iszero(v) continue end
    @inbounds (b,k)=Elist(n,z[i]*m)
    @inbounds res[k].+=b ? v : -v
  end
  Cyc(n,res)
else
  sumroots(n,i*m=>u for (i,u) in c.d)
end
end

function promote_conductor(a::Cyc,b::Cyc)
  if a.n==b.n return (a, b) end
  l=lcm(a.n,b.n)
  (raise(l,a),raise(l,b))
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{T2})where {T1,T2<:AbstractFloat}
  Complex{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Complex{T2}})where {T1,T2<:AbstractFloat}
  Complex{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Complex{T2}})where {T1,T2<:Integer}
  Cyc{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Complex{T2}})where {T1,T2<:Rational{<:Integer}}
  Cyc{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{T2})where {T1,T2<:Integer}
  Cyc{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{T2})where {T1,T2<:Rational{<:Integer}}
  Cyc{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Cyc{T2}})where {T1,T2}
  Cyc{promote_type(T1,T2)}
end

# total order is necessary to put Cycs in a sorted list
# for c.n==1  a<b is as expected
function Base.cmp(a::Cyc,b::Cyc)
  t=cmp(a.n,b.n)
  if !iszero(t) return t end
  a.n==1 ? cmp(num(a),num(b)) : cmp(a.d,b.d)
end

Base.:(==)(a::Cyc,b::Cyc)=a.n==b.n && a.d==b.d
Base.isless(a::Cyc,b::Cyc)=cmp(a,b)==-1
Base.isless(c::Cyc,d::Real)=c<Cyc(d)
Base.isless(d::Real,c::Cyc)=Cyc(d)<c

# hash is necessary to put Cycs as keys of a Dict or make a Set
Base.hash(a::Cyc, h::UInt)=hash(a.d, hash(a.n, h))

function Base.show(io::IO, ::MIME"text/html", a::Cyc)
  print(io, "\$")
  show(IOContext(io,:TeX=>true),a)
  print(io, "\$")
end

function Base.show(io::IO, ::MIME"text/plain", r::Cyc)
  if !haskey(io,:typeinfo) print(io,typeof(r),": ") end
  show(io,r)
end

function normal_show(io::IO,p::Cyc{T})where T
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if T<:Rational{<:Integer}
    den=denominator(p)
    p=Cyc{typeof(den)}(p*den)
  else
    den=1
  end
if use_list
  it=zip(zumbroich_basis(p.n),p.d)
else
  it=p.d
end
  res=join( map(it) do (deg,v)
if use_list
    if iszero(v) return "" end
end
    if deg==0 t=string(v)
    else 
      t=format_coefficient(string(v))
      if repl || TeX
        r="\\zeta"* (p.n==1 ? "" : p.n<10 ? "_$(p.n)" : "_{$(p.n)}")
        if deg>=1 r*= deg==1 ? "" : deg<10 ? "^$deg" : "^{$deg}" end
      else
        r=(deg==1 ? "E($(p.n))" : "E($(p.n),$deg)")
      end
      t*=r
    end
    if t[1]!='-' t="+"*t end
    t
  end)
  if res[1]=='+' res=res[2:end] end
  if !isone(den) 
    res=bracket_if_needed(res)
    res*=fromTeX(io,TeX ? "/{$den}" : repl ? "/$den" : "//$den" ) 
  end
  fromTeX(io,res)
end

function Base.show(io::IO, p::Cyc{T})where T
  quadratic=get(io,:quadratic,true)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if iszero(p)
    if repl||TeX print(io,"0")
    else print(io,"zero(Cyc{",T,"})")
    end
    return
  end
  rqq=[normal_show(io,p)]
  if quadratic && (T<:Integer || T<:Rational{<:Integer})
    q=Quadratic(p)
    if !isnothing(q)
      push!(rqq,repr(q;context=io))
    end
    for test in [1-E(4),1+E(4),E(3),E(3,2),1-E(3),1-E(3,2),1+E(3),1+E(3,2)]
      q=Quadratic(p/test)
      if !isnothing(q)
        rq=repr(q;context=io)
        rq=format_coefficient(rq;allow_frac=true)
        t=format_coefficient(normal_show(io,test))
        if !isempty(rq) && rq[1]=='-' rq="-"*t*rq[2:end] else rq=t*rq end
        push!(rqq,rq)
      end
    end
  end
  print(io,rqq[argmin(length.(rqq))])
end

Base.gcd(v::Vector{<:Cyc})=one(v[1])
Base.gcd(a::Cyc,b::Cyc)=one(a)
Base.gcd(a::Cyc,b::Number)=one(a)
Base.gcd(b::Number,a::Cyc)=one(a)

function Base.:+(x::Cyc,y::Cyc)
  a,b=promote(x,y)
  if iszero(a) return b
  elseif iszero(b) return a
  end
if use_list
  a,b=promote_conductor(a,b)
  lower(Cyc(a.n,a.d+b.d))
else
  n=lcm(a.n,b.n)
  na=div(n,a.n)
  nb=div(n,b.n)
  res=eltype(a.d)[]
  for (i,va) in a.d sumroot(res,n,na*i,va) end
  for (i,vb) in b.d sumroot(res,n,nb*i,vb) end
  lower(Cyc(n,MM(res)))
end
end

Base.:-(a::Cyc)=Cyc(a.n,-a.d)
Base.:-(a::Cyc,b::Cyc)=a+(-b)

if use_list
Base.div(c::Cyc,a::Real)=Cyc(c.n,div.(c.d,a))
Base.://(c::Cyc,a::Real)=Cyc(c.n,c.d.//a)
else
function Base.div(c::Cyc,a::Real)
  n=merge(div,c.d,a)
  Cyc(iszero(n) ? 1 : c.n,n)
end
Base.://(c::Cyc,a::Real)=Cyc(c.n,c.d//a)
Base.:*(c::Cyc,a::Real)=Cyc(iszero(a) ? 1 : c.n,c.d*a)
end
Base.://(a::Cyc,c::Cyc)=a*inv(c)
Base.://(a::Real,c::Cyc)=a*inv(c)
Base.:/(c::Cyc,a::Real)=c//a
Base.:/(a::Cyc,c::Cyc)=a//c
Base.:/(a::Real,c::Cyc)=a//c
Base.:*(a::Real,c::Cyc)=c*a

function Base.:*(a::Cyc,b::Cyc)
  a,b=promote(a,b)
  if iszero(a) return a end
  if iszero(b) return b end
if use_list
  a,b=promote_conductor(a,b)
  zb=zumbroich_basis(a.n)
  res=zero(a.d)
  for i in eachindex(a.d), j in eachindex(b.d)
@inbounds  c=a.d[i]*b.d[j]
    if iszero(c) continue end
@inbounds  (v,k)=Elist(a.n,zb[i]+zb[j])
@inbounds  res[k].+=v ? c : -c
  end
  lower(Cyc(a.n,res))
else
# a,b=promote_conductor(a,b)
# let ad=a.d,bd=b.d
#   res=sumroots(a.n,[i+j=>va*vb for (i,va) in ad, (j,vb) in bd])
# end
# lower(res)
#---------------------------------------
  if a.n==1 return Cyc(b.n,b.d*num(a))
  elseif b.n==1 return Cyc(a.n,a.d*num(b))
  end
  n=lcm(a.n,b.n)
  na=div(n,a.n)
  nb=div(n,b.n)
  let ad=a.d,bd=b.d,na=na,nb=nb
    lower(sumroots(n,na*i+nb*j=>va*vb for (i,va) in ad, (j,vb) in bd))
  end
# res=eltype(a.d)[]
# for (i,va) in a.d, (j,vb) in b.d sumroot(res,n,na*i+nb*j,va*vb) end
# lower(Cyc(n,MM(res)))
end
end

if use_list
function Cyc(m::Int,z::Vector{Int},c::Vector{T})where T
  zb=zumbroich_basis(m)
  res=zeros(T,length(zb))
  res[indexin(z,zb)]=c
  Cyc(m,res)
end
end

function lower(c::Cyc{T})where T # write c in smallest Q(ζ_n) where it sits
  n=c.n
# println("lowering $(c.n):$(c.d)")
  if n==1 return c end
if use_list
  nz=findall(!iszero,c.d)
  if length(nz)==0 return Cyc(c.d[1]) end
  zb=zumbroich_basis(n)
else
  if iszero(c) return zero(Cyc{T}) end
end
  for (p,np) in factor(n)
    m=div(n,p)
if use_list
    let m=m, zb=zb
    if np>1
      if all(z->z%p==0,zb[nz])
        return lower(Cyc(m,div.(zb[nz],p),c.d[nz]))
      end
    elseif count(!iszero,c.d)%(p-1)==0
      cnt=zeros(Int,m)
      for (i,v) in enumerate(zumbroich_basis(c.n)) 
        if !iszero(c.d[i]) cnt[1+(v%m)]+=1 end 
      end
      if all(x->iszero(x) || x==p-1,cnt) 
        u=findall(!iszero,cnt).-1
        kk=@. div(u+m*mod(-u,p)*invmod(m,p),p)%m
        if p==2 return lower(Cyc(m,kk,c.d[indexin((kk*p).%n,zb)]))
        elseif all(k->constant(c.d[indexin((m*(1:p-1).+k*p).%n,zb)]),kk)
          return lower(Cyc(m,kk,-c.d[indexin((m.+kk*p).%n,zb)]))
        end
      end
    end
    end
else
    if np>1 
      if all(k->first(k)%p==0,c.d) 
        return lower(Cyc(m,MM(div(k,p)=>v for (k,v) in c.d;check=false)))
      end
    elseif iszero(length(c.d)%(p-1))
      cnt=zeros(Int,m)
      for (k,v) in c.d cnt[1+(k%m)]+=1 end
      if all(x->iszero(x) || x==p-1,cnt) 
        u=findall(!iszero,cnt).-1
        kk=@. div(u+m*mod(-u,p)*invmod(m,p),p)%m
        if p==2  
          return lower(Cyc(m,MM(k=>c.d[(k*p)%n] for k in kk)))
        elseif all(k->constant(map(i->c.d[(m*i+k*p)%n],1:p-1)),kk)
          return lower(Cyc(m,MM(k=>-c.d[(m+k*p)%n] for k in kk)))
        end
      end
    end
end
  end
  c
end

galois(c::Rational,n::Int)=c

galois(c::Integer,n::Int)=c
"""
  galois(c::Cyc,n::Int) applies to c the galois automorphism
  of Q(ζ_conductor(c)) raising all roots of unity to the n-th power.
  n should be prime to conductor(c).
# Examples
```julia-repl
julia> galois(1+E(4),-1) # galois(c,-1) is the same as conj(c)
Cyc{Int64}: 1-ζ₄

julia> galois(ER(5),2)==-ER(5)
true
```
"""
function galois(c::Cyc,n::Int)
  if gcd(n,c.n)!=1 error("$n should be prime to conductor($c)=$(c.n)") end
if use_list
  zb=zumbroich_basis(c.n)
  nz=findall(!iszero,c.d)
  if isempty(nz) c
  else let zb=zb
     sum(t->c.d[t]*E(c.n,zb[t]*n),nz)
       end
  end
else
  if iszero(c) return c end
  let cn=c.n, n=n
  sumroots(c.n,(e*n)%cn=>p for (e,p) in c.d)
  end
end
end

Base.conj(c::Cyc)=galois(c,-1)

function Base.inv(c::Cyc)
  if c.n==1
    r=num(c)
    if r==1 || r==-1 return Cyc(r) else return Cyc(1//r) end
  end
  l=setdiff(unique(map(i->galois(c,i),prime_residues(c.n))),[c])
  r=prod(l)
  n=num(c*r)
  n==1 ? r : (n==-1 ? -r : r//n)
end

Base.:^(a::Cyc, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                   Base.power_by_squaring(inv(a),-n)

const ER_dict=Dict(1=>Cyc(1),-1=>E(4))
"""
  ER(n::Int) computes as a Cyc the square root of the integer n.
# Examples
```julia-repl
julia> ER(-1)
Cyc{Int64}: ζ₄

julia> ER(3)
Cyc{Int64}: √3
```
"""
function ER(n::Integer)
  get!(ER_dict,n) do 
  for (p,d) in factor(n)
    h=p^div(d,2)
    if h>1 return h*ER(div(n,h^2)) end
  end
  if n==0       return 0
  elseif n<0    return E(4)*ER(-n)
  elseif n%4==1 return sum(k->E(n,k^2),1:n)
  elseif n%4==2 return (E(8)-E(8,3))*ER(div(n,2))
  elseif n%4==3 return E(4,3)*sum(k->E(n,k^2),1:n)
  else          return 2*ER(div(n,4))
  end
  end
end 

Base.abs(c::Cyc)=c*conj(c)

#------------------------ type Root1 ----------------------------------
struct Root1 <: Number # E(c,n)
  r::Rational{Int}
  Root1(c::Int,n::Int)=new(mod(n,c)//c)
end

Cycs.conductor(a::Root1)=denominator(a.r)
Base.exponent(a::Root1)=numerator(a.r)

Root1(;r=0)=Root1(denominator(r),numerator(r)) # does a mod1
E(a::Root1)=E(;r=a.r)

function Root1(c::Real)
  if c==1 Root1(1,0)
  elseif c==-1 Root1(2,1)
  else nothing
  end
end

Base.broadcastable(r::Root1)=Ref(r)

function Base.show(io::IO, ::MIME"text/plain", r::Root1)
  if !haskey(io,:typeinfo) print(io,typeof(r),": ") end
  show(io,r)
end

function Base.show(io::IO, r::Root1)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  d=exponent(r)
  c=conductor(r)
  if repl || TeX
    if c==1 print(io,"1")
    elseif c==2 print(io,"-1")
    else r="\\zeta"* (c==1 ? "" : c<10 ? "_$(c)" : "_{$(c)}")
      if d>=1 r*=(d==1 ? "" : d<10 ? "^$d" : "^{$d}") end
      printTeX(io,r)
    end
  else
    print(io,"Root1($c,$d)")
  end
end

"""
`Root1(c)`
    
`c` should be a cyclotomic number (a `Cyc`), or a `Real`. `Root1` returns a
`Root1` object containing the rational `e/n` with `0≤e<n` (that is, `e/n∈ ℚ
/ℤ`) if `c==E(n,e)`, and `nothing` if `c` is not a root of unity.

```julia-repl
julia> r=Root1(-E(9,2)-E(9,5))
Root1: ζ₉⁸

julia> conductor(r)
9

julia> exponent(r)
8

julia> E(r)
Cyc{Int64}: -ζ₉²-ζ₉⁵

julia> Root1(-E(9,4)-E(9,5)) # nothing
```
""" 
function Root1(c::Cyc)
if use_list
  if !(all(x->x==0 || x==1,c.d) ||all(x->x==0 || x==-1,c.d))
    return nothing
  end
else
  if !(all(x->last(x)==1,c.d) || all(x->last(x)==-1,c.d))
    return nothing
  end
end
  for i in prime_residues(c.n)
    if c==E(c.n,i) return Root1(c.n,i) end
    if -c==E(c.n,i) 
      if c.n%2==0 return Root1(c.n,div(c.n,2)+i)
      else return Root1(2*c.n,c.n+2*i)
      end
    end
  end
  return nothing
end

function Base.cmp(a::Root1,b::Root1)
  r=cmp(conductor(a),conductor(b))
  if !iszero(r) return r end
  cmp(exponent(a),exponent(b))
end

Base.isless(a::Root1,b::Root1)=cmp(a,b)==-1
#Base.:(==)(a::Root1,b::Root1)=iszero(cmp(a,b))
Base.one(a::Root1)=Root1(1,0)
Base.isone(a::Root1)=iszero(a.r)
Base.:*(a::Root1,b::Root1)=Root1(;r=a.r+b.r)
Base.:^(a::Root1,n::Integer)=Root1(;r=n*a.r)
function Base.:^(a::Root1,r::Rational) # "canonical" way to extract roots
  n=denominator(r)
  d=conductor(a)
  j=1
  while true
    k=gcd(n,d)
    n=div(n,k)
    j*=k
    if k==1 break end
  end
  Root1(j*d,numerator(r)*exponent(a)*gcdx(n,d)[2])
end
Base.inv(a::Root1)=Root1(;r=-a.r)
Base.conj(a::Root1)=inv(a)
Base.:/(a::Root1,b::Root1)=a*inv(b)

if use_list
Base.:*(a::Cyc,b::Root1)=a*E(b)
else
function Base.:*(a::Cyc,b::Root1)
  n=lcm(conductor(a),conductor(b))
  na=div(n,conductor(a))
  nb=div(n,conductor(b))
  res=eltype(a.d)[]
  for (i,va) in a.d sumroot(res,n,na*i+nb*exponent(b),va) end
  lower(Cyc(n,MM(res)))
end
end

Base.:*(b::Root1,a::Cyc)=a*b
Base.:*(a::Union{Integer,Rational},b::Root1)=a*E(b)
Base.:*(b::Root1,a::Integer)=E(b)*a

#------------------- end of Root1 ----------------------------------------

struct Quadratic
  a
  b
  root
  den
end

"""
  `Quadratic(c::Cyc)` 
  
determines  if `c`  lives in  a quadratic  extension of  `ℚ`. It  returns a
`Quadratic`  struct with fields `a`, `b`, `root`, `den` representing `c` as
`(a  + b ER(root))//den`  if such a  representation is possibe or `nothing`
otherwise

# Examples
```julia-repl
julia> Quadratic(1+E(3))
(1+√-3)/2

julia> Quadratic(1+E(5))

```
"""
function Quadratic(cyc::Cyc{T})where T
  l1=coefficients(cyc)
  den=lcm(denominator.(l1))
  cyc*=den
  l=numerator.(coefficients(cyc))
  if cyc.n==1 return Quadratic(l[1],0,1,den) end

  f=factor(cyc.n)
  v2=get(f,2,0)

  if v2>3 || (v2==2 && any(p->p[1]!=2 && p[2]!=1,f)) ||
     (v2<2 && any(x->x!=1,values(f)))
    return nothing
  end

  f=keys(f)
  if v2==0
    root=cyc.n
    if root%4==3 root=-root end
    gal=Set(galois.(cyc,prime_residues(cyc.n)))
    if length(gal)!=2 return nothing end
    a=numerator(convert(T,sum(gal)))      # trace of 'cyc' over the rationals
    if length(f)%2==0 b=2*l[2]-a
    else b=2*l[2]+a
    end
    if a&1==0 && b&1==0 a>>=1; b>>=1; d=1
    else d=2
    end
  elseif v2==2
    root=cyc.n>>2
    if root==1 a=l[1];b=-l[2]
    else
      a=l[5]
      if length(f)%2==0 a=-a end
      b=-l[root+5]
    end
    if root%4==1 root=-root; b=-b end
    d=1
  else		# v2 = 3
    root=cyc.n>>2
    if root==2
      a=l[1];b=l[2]
      if b==l[4] root=-2 end
    else
      a=l[9]
      if length(f)%2==0 a=-a end
      b=l[(root>>1)+9]
      if b!=-l[3*(root>>1)-7] root=-root
      elseif (root>>1)%4==3 b=-b
      end
    end
    d=1
  end
  if d*cyc!=a+b*ER(root) return nothing end
  return Quadratic(a,b,root,den*d)
end

function Base.show(io::IO,q::Quadratic)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  rq=string(q.a)
  if q.b!=0 
    if iszero(q.a) rq=""
    elseif q.b>0 rq*="+" end
    rq*=q.b==1 ? "" : q.b==-1 ? "-" : string(q.b)
    r=string(q.root)
    rq*=TeX ? "\\sqrt{$r}" : repl ? "√$r" : "ER($r)"
    if !iszero(q.a) && q.den!=1 rq="("*rq*")" end
  end
  print(io,rq)
  if q.den!=1 && rq!="0" print(io,(repl||TeX) ? "/" : "//",q.den) end
end

const inforoot=Ref(true)
function proot(x,n,r)
  if inforoot[] 
    xprint("root(",x)
    if n!=2 xprint(",",n) end
    xprint(")=",r,"\n")
  end
end
const Irootdict=Dict{Tuple{Int,Int},Any}()

"""
`root(x,n=2)`

computes  the `n`-th root of `x` when we know  how to do it. We know how to
compute  `n`-th  roots  for  roots  of  unity, square roots of integers and
`n`-th  roots of  integers wich  are perfect  `n`-th powers  of integers or
square roots of integers.

```julia-repl
julia> root(-1)
Cyc{Int64}: ζ₄

julia> root(E(4))
Cyc{Int64}: ζ₈

julia> root(27,6)
Cyc{Int64}: √3
```
"""
function root(x::Integer,n=2)
  if isone(n) || isone(x) return x end
  if !(n isa Int) n=Int(n) end
  get!(Irootdict,(n,x)) do
    if x==1 || (x==-1 && n%2==1) return x end
    l=factor(x)
    if any(y->(2y)%n!=0,values(l)) error("root($x,$n) not implemented") end
    a=prod(p^div(pow,n) for (p,pow) in l)
    b=[p for (p,pow) in l if pow%n!=0]
    res=isempty(b) ? a : a*ER(prod(b))
    proot(x,n,res)
    res
  end
end

root(x::AbstractFloat,n)=x^(1//n);
root(x::Rational{<:Integer},n::Number=2)=root(numerator(x),n)//root(denominator(x),n)

const Crootdict=Dict{Tuple{Int,Cyc},Any}()
function root(x::Cyc,n=2)
  if isone(n) || isone(x) return x end
  if !(n isa Int) n=Int(n) end
  get!(Crootdict,(n,x)) do
  r=Root1(x)
  if isnothing(r) 
    if conductor(x)>1 return nothing end
    return root(num(x),n)
  end
  d=conductor(r)
  j=1
  n1=n
  while true
    k=gcd(n1,d)
    n1=div(n1,k)
    j*=k
    if k==1 break end
  end
  res=E(j*d,exponent(r)*gcdx(n1,d)[2])
  proot(x,n,res)
  res
  end
end

# 347.534 ms (4367402 allocations: 366.17 MiB) in 1.5.3
# 565.431 ms (5861810 allocations: 775.28 MiB) in 1.5.3, ModuleElts Dict
function testmat(p) # testmat(12)^2 takes 0.27s in 1.0
  ss=[[i,j] for i in 0:p-1 for j in i+1:p-1]
  [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
end
end
