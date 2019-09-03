"""
Cyclotomic  numbers means complex numbers which are sums of rationals times
roots of unity.

They are a very important feature of GAP, since entries of character tables
of finite groups are cyclotomics.

They  have a normal form given by the Zumbroich basis, which allows to find
the  smallest Cyclotomic field which contains a given number, and decide in
particular if a cyclotomic is zero. Let ζ_n:=e^{2iπ/n}. The Zumbroich basis
of Q(ζ_n) is a particular subset of 1,ζ,ζ^2,...,ζ^{n-1} which forms a basis
of Q(ζ_n) with good properties.

I  ported here Christian Stump's Sage  code, which is simpler to understand
than GAP's code. The reference for the algorithms is

T. Breuer, Integral bases for subfields of cyclotomic fields AAECC 8 (1997)

As in GAP, I lower automatically numbers after each computation; this makes
this code about twice slower than GAP since lower is not as much optimized.
GAP  also converts a Cyclotomic which is rational to a Rational, a Rational
which is integral to an Int, etc... This is tremendously useful but needs a
new  type of  number to  be added  to Julia,  which requires more competent
people than me.

The main way to build a Cyclotomic number is to use the function `E(n,k=1)`
which constructs ζ_n^k.

# Examples
```julia-repl
julia> E(3)+E(4)
ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹

julia> E(3,2)
ζ₃²

julia> 1+E(3,2)
-ζ₃

julia> a=E(4)-E(4)
0

julia> conductor(a) # a is lowered to Q(ζ_1)=Q
1

julia> typeof(convert(Int,a))
Int64

julia> convert(Int,E(4))
ERROR: InexactError: convert(Int64, E(4))

julia> c=inv(1+E(4)) # inverses need Rationals
1/2-ζ₄/2

julia> typeof(c)
Cyc{Rational{Int64}}

julia> typeof(1+E(4))
Cyc{Int64}

julia> Cyc(1+im) # one can convert Gaussian integers or rationals
1+ζ₄

julia> 1//(1+E(4))
1/2-ζ₄/2

julia> typeof(Cyc(1//2)) # another way of building a Cyc
Cyc{Rational{Int64}}

julia> conj(1+E(4))
1-ζ₄

julia> c=E(9)   # an effect of the Zumbroich basis
-ζ₉⁴-ζ₉⁷

julia> Root1(c) # but you can decide whether a Cyc is a root of unity
Root1(1//9)

julia> c=Complex(E(3))   # convert to float is probably not very useful
-0.4999999999999998 + 0.8660254037844387im
```

For more information see ER, Quadratic, galois. 

Finally, a benchmark:

```benchmark
julia> function testmat(p) 
         ss=vcat([[[i,j] for j in i+1:p-1] for i in 0:p-1]...)
         [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
       end
testmat (generic function with 1 method)

julia> @btime testmat(12)^2;
  472.964 ms (8324504 allocations: 707.18 MiB)
```
The equivalent in GAP:

```
testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
end; 
```

for testmat(12) takes 0.4s in GAP3, 0.3s in GAP4
"""
module Cycs
export E, ER, Cyc, conductor, galois, Root1, Quadratic

using Gapjm
import ..Util: ModuleElt, norm!
import ..Util: TeXstrip, bracket_if_needed, constant
import ..Util: factor, prime_residues, phi

const use_list=false
if use_list
struct Cyc{T <: Real}<: Number   # a cyclotomic number
  n::Int
  d::Vector{T} # coefficients on the Zumbroich basis
               # the i-th element is the coefficient on zumbroich_basis[i]
end
else
struct Cyc{T <: Real}<: Number   # a cyclotomic number
  n::Int
  d::ModuleElt{Int,T} # coefficients on the Zumbroich basis
end
end

conductor(c::Cyc)=c.n
conductor(a::Array)=lcm(conductor.(a))
conductor(i::Int)=1

const zumbroich_basis_dict=Dict(1=>[0])
"""
  zumbroich_basis(n::Int) returns the Zumbroich basis of Q(ζ_n)
  as the vector of i in 0:n-1 such that ``ζ_n^i`` is in the basis
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

#=
  Elist(n,i)  expresses ζ_n^i  in zumbroich_basis(n):  it is  a sum of some
  ζ_n^j  with coefficients all 1 or all  -1. The result is a Pair sgn=>inds
  where sgn is true if coefficients are all 1 and false otherwise, and inds
  is  the  list  of  i  in  0:n-1  such  that  ζ_n^i occurs with a non-zero
  coefficient.
=#
const Elist_dict=Dict((1,0)=>(true=>[0]))
function Elist(n::Int,i1::Int=1)::Pair{Bool,Vector{Int}}
  i=mod(i1,n)
  get!(Elist_dict,(n,i)) do
if use_list
    z=zumbroich_basis(n)
end
    nfact=factor(n)
    mp=Int[]
    j=i
    for (p,np) in nfact
      factor=p^np
      m=div(n,factor)
      cnt=mod(j*invmod(m,factor),factor)
      j-=cnt*m
      if p==2
        if 1==div(cnt,div(factor,p)) push!(mp,p) end
      else
        tmp=zeros(Int,np)
        for k in 1:np
          factor=div(factor,p)
          tmp[k],cnt=divrem(cnt,factor)
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
    if isempty(mp) 
      return true=> use_list ? Vector{Int}(indexin([i],z)) : [i] 
    end
    v=vec(sum.(Iterators.product((div(n,p)*(1:p-1) for p in mp)...)))
if use_list
    iseven(length(mp))=>Vector{Int}(indexin((i .+ v).%n,z))
else
    iseven(length(mp))=>sort((i .+ v).%n)
end
  end
end

#=
  E(n::Integer,k::Integer=1) is exp(2i k π/n)
=#
const E_dict=Dict((1,0)=>Cyc(1, use_list ? [1] : ModuleElt(0=>1)))
function E(n1::Integer,i1::Integer=1)::Cyc{Int}
  n=Int(n1)
  i=mod(Int(i1),n)
  get!(E_dict,(n,i)) do
    s,l=Elist(n,i) #::Pair{Bool,Vector{Int}}
if use_list
  v=zeros(Int,length(zumbroich_basis(n)))
  v[l].=ifelse(s,1,-1)
  lower(Cyc(n,v))
else
  lower(Cyc(n,ModuleElt(l.=>ifelse(s,1,-1))))
end
  end
end

E(a,b=1)=Cycs.E(Int(a),Int(b))

if use_list
Base.zero(c::Cyc)=Cyc(1,eltype(c.d)[0])
Base.zero(::Type{Cyc{T}}) where T=Cyc(1,T[0])
Base.iszero(c::Cyc)=c.n==1 && iszero(c.d[1])
else
Base.zero(c::Cyc)=Cyc(1,zero(c.d))
Base.zero(::Type{Cyc{T}}) where T=Cyc(1,zero(ModuleElt{Int,T}))
Base.iszero(c::Cyc)=iszero(c.d)
end
Base.one(c::Cyc)=E(1,0)

function Cyc(c::Complex{T}) where T
  if iszero(imag(c)) return Cyc(real(c)) end
if use_list
  if iszero(real(c)) return Cyc(4,[0,imag(c)])
  else return Cyc(4,[real(c),imag(c)])
  end
else
  if iszero(real(c)) return Cyc(4,ModuleElt(1=>imag(c)))
  else return Cyc(4,ModuleElt(0=>real(c),1=>imag(c)))
  end
end
end

if use_list
Cyc(i::T) where T<:Real=Cyc(1,[i])
else
Cyc(i::T) where T<:Real=iszero(i) ? zero(Cyc{T}) : Cyc(1,ModuleElt(0=>i))
end
Base.convert(::Type{Cyc{T}},i::Real) where T<:Real=Cyc(convert(T,i))
Cyc{T}(i::T1) where {T<:Real,T1<:Real}=convert(Cyc{T},i)

function Base.convert(::Type{Cyc{T}},c::Cyc{T1}) where {T,T1}
  if T==T1 return c end
if use_list
  Cyc(c.n,convert.(T,c.d))
else
  Cyc(c.n,ModuleElt(n=>convert(T,v) for (n,v) in c.d))
end
end

Cyc{T}(c::Cyc{T1}) where {T,T1}=convert(Cyc{T},c)

if use_list
 num(c::Cyc)=c.d[1]
else
 num(c::Cyc)=first(c.d)[2]
end

function Base.convert(::Type{T},c::Cyc)::T where T<:Real
  if c.n!=1 throw(InexactError(:convert,T,c)) end
if !use_list
  if iszero(c) return zero(T) end
end
  convert(T,num(c))
end

Int(c::Cyc)=convert(Int,c)

Base.isinteger(c::Cyc)=isone(conductor(c)) && (iszero(c) || isinteger(num(c)))

function promote_conductor(a::Cyc,b::Cyc)
  if a.n==b.n return (a, b) end
  l=lcm(a.n,b.n)
  (raise(l,a),raise(l,b))
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{T2})where {T1<:Real,T2<:Real}
  Cyc{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Cyc{T2}})where {T1<:Real,T2<:Real}
  Cyc{promote_type(T1,T2)}
end

" isless is necessary to put Cycs in a sorted list"
function Base.cmp(a::Cyc,b::Cyc)
  t=cmp(a.n,b.n)
  if !iszero(t) return t end
  cmp(a.d,b.d)
end

Base.isless(a::Cyc,b::Cyc)=cmp(a,b)==-1
Base.:(==)(a::Cyc,b::Cyc)=cmp(a,b)==0

" hash is necessary to put Cycs as keys of a Dict"
function Base.hash(a::Cyc, h::UInt)
   b = 0x595dee0e71d271d0%UInt
if use_list
   for i in a.d
     b = xor(b,xor(hash(i, h),h))
     b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
else
   for (k,v) in a.d
     b = xor(b,xor(hash(k, h),h))
     b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
     b = xor(b,xor(hash(v, h),h))
     b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
end
   return b
end

Base.show(io::IO,p::Cyc)=format(io,p;repl=get(io,:limit,false))

function Util.format(io::IO, p::Cyc; opt...)
  opt=Dict(opt)
  quadratic=get(opt,:quadratic,true)
  repl=get(opt,:repl,false)
  TeX=get(opt,:TeX,false)
  rq=""
  if quadratic
    q=Quadratic(p)
    if !isnothing(q)
      rq=sprint((io,x)->format(io,x;opt...),q;context=io)
    end
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
    den=denominator(v)
    v=numerator(v) 
    if deg==0 t=string(v)
    else 
      if repl || TeX
        t= v==1 ? "" : v==-1 ? "-" : bracket_if_needed(string(v))
        r="\\zeta"* (p.n==1 ? "" : p.n<10 ? "_$(p.n)" : "_{$(p.n)}")
        if deg>=1 r*= (deg==1 ? "" : deg<10 ? "^$deg" : "^{$deg}") end
      else
        v=bracket_if_needed(string(v))
        t= v=="1" ? "" : v=="-1" ? "-" : v
        r=(deg==1 ? "E($(p.n))" : "E($(p.n),$deg)")
      end
      t*=r
    end
    if t[1]!='-' t="+"*t end
    if !isone(den) t*="/$den" end
    t
  end)
  if res=="" res="0"
  elseif res[1]=='+' res=res[2:end] 
  end
  res=(repl && !TeX) ? TeXstrip(res) : res
  print(io, (quadratic && rq!="" && length(rq)<length(res)) ? rq : res)
end

function Base.:+(x::Cyc,y::Cyc)
# if iszero(x) return y
# elseif iszero(y) return x
# end
  a,b=promote(x,y)
  a,b=promote_conductor(a,b)
if use_list
  lower(Cyc(a.n,a.d.+b.d))
else
  lower(Cyc(a.n,a.d+b.d))
end
end

Base.:+(a::Cyc,b::Number)=a+Cyc(b)
Base.:+(b::Number,a::Cyc)=a+Cyc(b)

Base.:-(a::Cyc)=Cyc(a.n,-a.d)
Base.:-(a::Cyc,b::Cyc)=a+(-b)
Base.:-(b::Number,a::Cyc)=Cyc(b)-a
Base.:-(b::Cyc,a::Number)=b-Cyc(a)

Base.:*(a::Real,c::Cyc)= iszero(a) ? zero(c) : Cyc(c.n,c.d*a)
Base.:*(c::Cyc,a::Real)=a*c

if use_list
Base.div(c::Cyc,a::Real)=Cyc(c.n,div.(c.d,a))
Base.://(c::Cyc,a::Real)=Cyc(c.n,c.d.//a)
else
Base.div(c::Cyc,a::Real)=Cyc(c.n,ModuleElt(k=>div(v,a) for (k,v) in c.d))
Base.://(c::Cyc,a::Real)=Cyc(c.n,ModuleElt(k=>v//a for (k,v) in c.d))
end
Base.://(a::Cyc,c::Cyc)=a*inv(c)
Base.://(a::Real,c::Cyc)=a*inv(c)
Base.:/(c::Cyc,a::Real)=c//a
Base.:/(a::Cyc,c::Cyc)=a//c
Base.:/(a::Real,c::Cyc)=a//c

function addroot!(l::ModuleElt,n::Int,p)
  (i,c)=p
  (s,v)=Elist(n,i)
  if !s c=-c end
if ModuleElts.usedict
  for k in v if haskey(l,k) l.d[k]+=c else push!(l,k=>c) end end
else
  for k in v push!(l,k=>c) end
end
end

function sumroots(n::Int,l::Vector{Pair{K,V}})where {K,V}
  res=zero(ModuleElt{K,V})
  for p in l addroot!(res,n,p) end
  Cyc(n,norm!(res))
end

function Base.:*(a::Cyc,b::Cyc)
  a,b=promote(a,b)
  if iszero(a) return a end
  if iszero(b) return b end
  if a.n==1 return num(a)*b end
  if b.n==1 return num(b)*a end
  a,b=promote_conductor(a,b)
if use_list
  zb=zumbroich_basis(a.n)
  res=zero(a.d)
  for i in eachindex(a.d), j in eachindex(b.d)
@inbounds  c=a.d[i]*b.d[j]
    if iszero(c) continue end
    (v,k)=Elist(a.n,zb[i]+zb[j])
@inbounds  res[k].+=v ? c : -c
  end
  lower(Cyc(a.n,res))
else
#  lower(sumroots(a.n,[i+j=>va*vb for (i,va) in a.d for (j,vb) in b.d]))
   res=zero(a.d)
   for (i,va) in a.d, (j,vb) in b.d addroot!(res,a.n,i+j=>va*vb) end
   lower(Cyc(a.n,norm!(res)))
end
end

function raise(n::Int,c::Cyc) # write c in Q(ζ_n) if c.n divides n
  if n==c.n return c end
  m=div(n,c.n)
if use_list
  res=zeros(eltype(c.d),length(zumbroich_basis(n)))
  z=zumbroich_basis(c.n).*m
  for i in eachindex(c.d)
    if iszero(c.d[i]) continue end
    (v,k)=Elist(n,z[i])
    @inbounds res[k].+=v ? c.d[i] : -c.d[i]
  end
  Cyc(n,res)
else
# sumroots(n,[i*m=>u for (i,u) in c.d])
  res=zero(c.d)
  for (i,u) in c.d addroot!(res,n,i*m=>u) end
  Cyc(n,norm!(res))
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
  nz=findall(x->!iszero(x),c.d)
  if length(nz)==0 return Cyc(c.d[1]) end
  zb=zumbroich_basis(n)
else
  if iszero(c) return zero(Cyc{T}) end
end
  for (p,np) in factor(n)
    m=div(n,p)
if use_list
    let p=p, m=m, zb=zb
    if np>1
      if all(z->z%p==0,zb[nz])
        return lower(Cyc(m,div.(zb[nz],p),c.d[nz]))
      end
    elseif count(x->!iszero(x),c.d)%(p-1)==0
      cnt=zeros(Int,m)
      for (i,v) in enumerate(zumbroich_basis(c.n)) 
        if !iszero(c.d[i]) cnt[1+(v%m)]+=1 end 
      end
      if all(x->iszero(x) || x==p-1,cnt) 
        u=findall(x->!iszero(x),cnt).-1
        kk=[div(k+m*mod(-k,p)*invmod(m,p),p)%m for k in u]
        if p==2 return lower(Cyc(m,kk,c.d[indexin((kk*p).%n,zb)]))
        elseif all(k->constant(c.d[indexin((m*(1:p-1).+k*p).%n,zb)]),kk)
         return lower(Cyc(m,kk,-c.d[indexin((m.+kk*p).%n,zb)]))
        end
      end
    end
    end
else
    if np>1 
      let p=p
      if all(k->k[1]%p==0,c.d) 
        return lower(Cyc(m,ModuleElt(div(k,p)=>v for (k,v) in c.d)))
      end
      end
    elseif iszero(length(c.d)%(p-1))
      cnt=zeros(Int,m)
      for (k,v) in c.d cnt[1+(k%m)]+=1 end
      let p=p, m=m, n=n
      if all(x->iszero(x) || x==p-1,cnt) 
        u=findall(x->!iszero(x),cnt).-1
        kk=sort(Int.([div(k+m*mod(-k,p)*invmod(m,p),p)%m for k in u]))
        if p==2 return lower(Cyc(m,ModuleElt(k=>c.d[(k*p)%n] for k in kk)))
        elseif all(k->constant(map(i->c.d[(m*i+k*p)%n],1:p-1)),kk)
          return lower(Cyc(m,ModuleElt(k=>-c.d[(m+k*p)%n] for k in kk)))
        end
      end
      end
    end
end
  end
  c
end

"""
  galois(c::Cyc,n::Int) applies to c the galois automorphism
  of Q(ζ_conductor(c)) raising all roots of unity to the n-th power.
  n should be prime to c.
# Examples
```julia-repl
julia> galois(1+E(4),-1) # galois(c,-1) is the same as conj(c)
1-ζ₄

julia> galois(ER(5),2)==-ER(5)
true
```
"""
function galois(c::Cyc,n::Int)
  if gcd(n,c.n)!=1 error("$n should be prime to conductor($c)=$(c.n)") end
if use_list
  zb=zumbroich_basis(c.n)
  nz=findall(x->!iszero(x),c.d)
  if isempty(nz) c
  else let zb=zb
     sum(t->c.d[t]*E(c.n,zb[t]*n),nz)
       end
  end
else
  sumroots(c.n,[(e*n)%c.n=>p for (e,p) in c.d])
end
end

Base.conj(c::Cyc)=galois(c,-1)

function Base.inv(c::Cyc)
  c=lower(c)
  if c.n==1
    r=Real(c)
    if r^2==1 return Cyc(r) else return Cyc(1//r) end
  else
    r=prod(i->galois(c,i),prime_residues(c.n)[2:end])
    n=Real(c*r)
  end
  n==1 ? r : (n==-1 ? -r : r//n)
end

Base.:^(a::Cyc, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                                    Base.power_by_squaring(inv(a),-n)

const ER_dict=Dict(1=>Cyc(1),-1=>E(4))
"""
  ER(n::Int) computes as a Cyc the square root of the integer n.
# Examples
```julia-repl
julia> ER(-1)
ζ₄

julia> ER(3)
√3
```
"""
function ER(n::Int)
  get!(ER_dict,n) do 
  if n==0       return 0
  elseif n<0    return E(4)*ER(-n)
  elseif n%4==1 return sum(k->E(n,k^2),1:n)
  elseif n%4==2 return (E(8)-E(8,3))*ER(div(n,2))
  elseif n%4==3 return E(4,3)*sum(k->E(n,k^2),1:n)
  else          return 2*ER(div(n,4))
  end
  end
end 

if use_list
Base.Complex(c::Cyc)=sum(x->x[2]*exp(2*x[1]*im*pi/c.n),
                          zip(zumbroich_basis(c.n),c.d))
else
Base.Complex(c::Cyc)=iszero(c) ? Complex(0.0) : 
             sum(v*exp(2*k*im*pi/c.n) for (k,v) in c.d)
end

function Base.Real(c::Cyc{T}) where T
if !use_list
  if iszero(c) return zero(T) end
end
  if c.n==1 return num(c) end
  if c!=conj(c) error("$c is not real") end
  return real(Complex(c))
end

Base.:<(c::Cyc,d::Real)=Real(c)<d

Base.abs(c::Cyc)=c*conj(c)

struct Root1 # E(c,n)
  r::Rational{Int}
  Root1(n::Int,c::Int)=new(mod(n,c)//c)
end

Root1(r::Rational)=Root1(numerator(r),denominator(r))

function Root1(c::Cyc)
if use_list
  if !(all(x->x==0 || x==1,c.d) ||all(x->x==0 || x==-1,c.d))
    return nothing
  end
else
  if !(all(x->x[2]==1,c.d) || all(x->x[2]==-1,c.d))
    return nothing
  end
end
  for i in 0:c.n-1
    if c==E(c.n,i) return Root1(i,c.n) end
    if -c==E(c.n,i) 
      if c.n%2==0 return Root1(div(c.n,2)+i,c.n)
      else return Root1(c.n+2*i,2*c.n)
      end
    end
  end
  return nothing
end

function Root1(c::Real)
  if c==1 Root1(0,1)
  elseif c==-1 Root1(1,2)
  else nothing
  end
end

Cycs.conductor(a::Root1)=a.r.den
Base.exponent(a::Root1)=a.r.num

function Base.cmp(a::Root1,b::Root1)
  r=cmp(conductor(a),conductor(b))
  if !iszero(r) return r end
  cmp(exponent(a),exponent(b))
end

Base.isless(a::Root1,b::Root1)=cmp(a,b)==-1

function Base.:*(a::Root1,b::Root1)
  r=a.r+b.r
  Root1(numerator(r),denominator(r))
end

E(a::Root1)=E(conductor(a),exponent(a))
E(a::Rational{<:Integer})=E(denominator(a),numerator(a))


struct Quadratic{T<:Integer}
  a::T
  b::T
  root::T
  den::T
end

"""
  `Quadratic(c::Cyc)` 
  
  determines  if `c` lives  in a quadratic  extension of `ℚ`.  It returns a
  `Quadratic`  object representing `c` as `(a  + b ER(root))//d` or nothing
  if no such tuple exists

# Examples
```julia-repl
julia> Quadratic(1+E(3))
(1+√-3)/2

julia> Quadratic(1+E(5))

```
"""
function Quadratic(cyc::Cyc)
if use_list
  den=lcm(denominator.(cyc.d))::Int
  zb=zumbroich_basis(cyc.n)
  l=fill(0,cyc.n)
  l[1 .+ zb]=Int.(cyc.d*den)
else
  den=lcm(denominator.(map(x->x[2],cyc.d)))::Int
  if den!=1 cyc=Cyc(cyc.n,ModuleElt([k=>Int(v*den) for (k,v) in cyc.d])) end
  l=fill(0,cyc.n)
  for (k,v) in cyc.d l[k+1]=v end
end
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
    a=convert(Int,sum(gal))          # trace of 'cyc' over the rationals
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

function Util.format(io::IO,q::Quadratic;repl=false,TeX=false,opt...)
  rq=string(q.a)
  if q.b!=0 
    if iszero(q.a) rq=""
    elseif q.b>0 rq*="+" end
    rq*=q.b==1 ? "" : (q.b==-1 ? "-" : string(q.b))
    rq*=repl ? "√$(q.root)" : (TeX ? "\\sqrt{$(q.root)}" : "ER($(q.root))")
    if !iszero(q.a) && q.den!=1 rq="("*rq*")" end
  end
  if q.den!=1 && rq!="0" rq*="/$(q.den)" end
  print(io,rq)
end

Base.show(io::IO,q::Quadratic)=format(io,q;repl=get(io,:limit,false))

function Gapjm.root(x::Cyc,n::Number=2)
  r=Root1(x)
  n1=Int(n)
  if isnothing(r) 
    if conductor(x)>1 return nothing end
    x=Real(x)
    if denominator(x)>1 return nothing end
    return root(Int(x),n)
  end
  d=denominator(r.r)
  j=1
  while true
    k=gcd(n1,d)
    n1=div(n1,k)
    j*=k
    if k==1 break end
  end
  res=E(j*d,numerator(r.r)*gcd_repr(n1,d)[1])
  println("root($x,$n)=$res")
  res
end

function testmat(p) # testmat(12)^2 takes 0.27s in 1.0
  ss=vcat([[[i,j] for j in i+1:p-1] for i in 0:p-1]...)
  [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
end
end
