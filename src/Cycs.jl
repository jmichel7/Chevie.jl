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
1//2+(-1//2)ζ₄

julia> typeof(c)
Cyc{Rational{Int64}}

julia> typeof(1+E(4))
Cyc{Int64}

julia> Cyc(1+im) # one can convert Gaussian integers or rationals
1+ζ₄

julia> 1//(1+E(4))
1//2+(-1//2)ζ₄

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

julia> Cyc(c) # even less useful
-0.4999999999999998+0.8660254037844387ζ₄
```

For more information see ER, quadratic, galois. 

Finally, a benchmark:

```benchmark
julia> function testmat(p) 
         ss=vcat([[[i,j] for j in i+1:p-1] for i in 0:p-1]...)
         [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
       end
testmat (generic function with 1 method)

julia> @btime testmat(12)^2;
  261.577 ms (3012634 allocations: 290.07 MiB)
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
export E, ER, Cyc, conductor, galois, Root1, quadratic

using ..Util

struct Cyc{T <: Real}<: Number   # a cyclotomic number
  n::Int
  d::SortedPairs{Int,T} # coefficients on the Zumbroich basis
end

conductor(c::Cyc)=c.n
conductor(a::Array)=lcm(conductor.(a))
conductor(i::Int)=1

"""
  zumbroich_basis(n::Int) returns the Zumbroich basis of Q(ζ_n)
  as the vector of i in 0:n-1 such that ``ζ_n^i`` is in the basis
"""
function zumbroich_basis(n::Int)::Vector{Int}
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
  v=(cartesian(res...)*fill(1,length(res)))::Vector{Int}
  sort(v.%n)
end

#=
  Elist(n,i)  expresses ζ_n^i  in zumbroich_basis(n):  it is  a sum of some
  ζ_n^j  with coefficients all 1 or all  -1. The result is a Pair sgn=>inds
  where sgn is true if coefficients are all 1 and false otherwise, and inds
  is  the  list  of  i  in  0:n-1  such  that  ζ_n^i occurs with a non-zero
  coefficient.
=#
const ElDict=Dict((1,0)=>(true=>[0]))
function Elist(n::Int,i::Int=1)::Pair{Bool,Vector{Int}}
  get!(ElDict,(n,i)) do
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
    if isempty(mp) return true=>[i%n] end
    v=cartesian((1:p-1 for p in mp)...)*div.(n,mp)
    iseven(length(mp))=>sort((i .+ v).%n)
  end
end

#=
  E(n::Integer,k::Integer=1) is exp(2i k π/n)
=#
const EDict=Dict((1,0)=>Cyc(1,[0=>1]))
function E(n::Integer,i::Integer=1)::Cyc{Int}
  get!(EDict,(n,i%n)) do
    s,l=Elist(n,i%n) #::Pair{Bool,Vector{Int}}
    Cyc(n,l.=>ifelse(s,1,-1))
  end
end

function Cyc(c::Complex{T}) where T
  if iszero(imag(c)) return Cyc(real(c)) end
  if iszero(real(c)) return Cyc(4,[1=>imag(c)])
  else return Cyc(4,[0=>real(c),1=>imag(c)])
  end
end

function Base.convert(::Type{Cyc{T}},i::Real)where T<:Real
  Cyc(1,iszero(i) ? Pair{Int,T}[] : [0=>convert(T,i)])
end

Cyc(i::T) where T<:Real=Cyc(1,iszero(i) ? Pair{Int,T}[] : [0=>i])
Cyc{T}(i::T1) where {T<:Real,T1<:Real}=convert(Cyc{T},i)

function Base.convert(::Type{Cyc{T}},c::Cyc{T1}) where {T,T1}
  if T==T1 return c end
  Cyc(c.n,convert.(Pair{Int,T},c.d))
end

Cyc{T}(c::Cyc{T1}) where {T,T1}=convert(Cyc{T},c)

function Base.convert(::Type{T},c::Cyc)::T where T<:Real
  if c.n!=1 throw(InexactError(:convert,T,c)) end
  if iszero(c) return zero(T) end
  convert(T,c.d[1][2])
end

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

Base.zero(c::Cyc)=Cyc(c.n,empty(c.d))
Base.zero(::Type{Cyc{T}},n::Int=1) where T=Cyc(n,Pair{Int,T}[])
Base.iszero(c::Cyc)=isempty(c.d)
Base.one(c::Cyc)=E(1,0)

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
   for i in a.d
     b = xor(b,xor(hash(i[1], h),h))
     b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
     b = xor(b,xor(hash(i[2], h),h))
     b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

function Base.show(io::IO, p::Cyc)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if get(io,:quadratic,false) && !(nothing==(q=quadratic(p)))
    if q.b==0 print(io,q.a)
    else
      if !iszero(q.a) 
        if q.den!=1 print(io,"(") end
        print(io,q.a)
        if q.b>0 print(io,"+") end
      end
      print(io,q.b==1 ? "" : q.b==-1 ? "-" : q.b)
      print(io,"√$(q.root)")
      if !iszero(q.a) && q.den!=1 print(io,")") end
    end
    if q.den!=1 && !iszero(q.a) print(io,"//$(q.den)") end
    return
  end
  res=join(map(p.d) do (deg,v)
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
    t
  end)
  if res=="" res="0"
  elseif res[1]=='+' res=res[2:end] 
  end
  print(io, (repl && !TeX) ? TeXstrip(res) : res)
end

function Base.:+(x::Cyc,y::Cyc)
  a,b=promote(x,y)
  a,b=promote_conductor(a,b)
  lower(Cyc(a.n,mergesum(a.d,b.d)))
end

Base.:+(a::Cyc,b::Number)=a+Cyc(b)
Base.:+(b::Number,a::Cyc)=a+Cyc(b)

Base.:-(a::Cyc)=Cyc(a.n,[k=>-v for (k,v) in a.d])
Base.:-(a::Cyc,b::Cyc)=a+(-b)
Base.:-(b::Number,a::Cyc)=Cyc(b)-a
Base.:-(b::Cyc,a::Number)=b-Cyc(a)

Base.:*(a::Real,c::Cyc)=iszero(a) ? zero(c) : Cyc(c.n,[k=>a*v for (k,v) in c.d])
Base.:*(c::Cyc,a::Real)=a*c

Base.div(c::Cyc,a::Real)=Cyc(c.n,[k=>Div(v,a) for (k,v) in c.d])

Base.://(c::Cyc,a::Real)=Cyc(c.n,[k=>v//a for (k,v) in c.d])
Base.://(a::Cyc,c::Cyc)=a*inv(c)
Base.://(a::Real,c::Cyc)=a*inv(c)

#function sumroots(n::Int,l::SortedPairs{Int,T})where T<:Real
function sumroots(n::Int,l)
  res=empty(l)
  sizehint!(res,phi(n)*length(l))
# res=Pair{Int,T}[]
  for (i,b) in l
    (s,v)=Elist(n,mod(i,n))#::Pair{Bool,Vector{Int}}
    if !s b=-b end
    for j in v 
      push!(res,j=>b)#::Pair{Int,T}) 
    end
#   append!(res,v.=>b)
  end
  lower(Cyc(n,norm!(res)))
end

function Base.:*(a::Cyc,b::Cyc)
  if iszero(a) return a end
  if iszero(b) return b end
  if a.n==1 return a.d[1][2]*b end
  if b.n==1 return b.d[1][2]*a end
  a,b=promote(a,b)
  a,b=promote_conductor(a,b)
# sumroots(a.n,[i+j=>va*vb for (i,va) in a.d for (j,vb) in b.d])
  res=empty(a.d)
  sizehint!(res,length(a.d)*length(b.d)*phi(a.n))
  for (i,va) in a.d for (j,vb) in b.d
    (s,v)=Elist(a.n,mod(i+j,a.n))#::Pair{Bool,Vector{Int}}
    u=va*vb
    if !s u=-u end
    for k in v push!(res,k=>u) end
  end end
  lower(Cyc(a.n,norm!(res)))
end

function raise(n::Int,c::Cyc{T})where T # write c in Q(zeta_n) if c.n divides n
  if n==c.n return c end
  if iszero(c) return zero(Cyc{T},n) end
  m=div(n,c.n)
# let m=m, c=c
# sumroots(n,[(i*m=>e) for (i,e) in c.d])
# end
  res=empty(c.d)
  sizehint!(res,phi(n)*length(c.d))
  for (i,u) in c.d
    (s,v)=Elist(n,mod(i*m,n))#::Pair{Bool,Vector{Int}}
    if !s u=-u end
    for k in v push!(res,k=>u) end
  end
  Cyc(n,norm!(res))
end

function cntmod(m::Int,p::Int,c::Cyc)
  cnt=zeros(Int,m)
  for (k,v) in c.d cnt[1+(k%m)]+=1 end
  all(x->iszero(x) || x==p-1,cnt)
end

function lower(c::Cyc{T})where T # write c in smallest Q(ζ_n) where it sits
  n=c.n
# println("lowering $(c.n):$(c.d)")
  if n==1 return c end
  if iszero(c) return zero(Cyc{T}) end
  for (p,np) in factor(n)
    m=div(n,p)
    let p=p, m=m
    if np>1 
      if all(k->k[1]%p==0,c.d) 
        return lower(Cyc(m,[div(k,p)=>v for (k,v) in c.d]))
      end
    elseif cntmod(m,p,c) # if all(v->length(v)==p-1,values(res))
      res=groupby(x->x%m,first.(c.d))
      kk=sort(Int.([div(k+m*mod(-k,p)*invmod(m,p),p)%m for k in keys(res)]))
      if p==2 return lower(Cyc(m,[k=>getvalue(c.d,(k*p)%n) for k in kk]))
      elseif all(k->constant(map(i->getvalue(c.d,(m*i+k*p)%n),1:p-1)),kk)
        return lower(Cyc(m,[k=>-getvalue(c.d,(m+k*p)%n) for k in kk]))
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
  sumroots(c.n,[(e*n)%c.n=>p for (e,p) in c.d])
end

Base.conj(c::Cyc)=galois(c,-1)

function Base.inv(c::Cyc)
  if c.n==1
    r=Real(c)
    if r^2==1 return r else return Cyc(inv(r)) end
  else
    r=prod(i->galois(c,i),prime_residues(c.n)[2:end])
    n=Real(c*r)
  end
  n==1 ? r : (n==-1 ? -r : r//n)
end

Base.:^(a::Cyc, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                                    Base.power_by_squaring(inv(a),-n)

"""
  ER(n::Int) computes as a Cyc the square root of the integer n.
# Examples
```julia-repl
julia> ER(-3)
ζ₃-ζ₃²

julia> ER(3)
-ζ₁₂⁷+ζ₁₂¹¹
```
"""
function ER(n::Int)
  if n==0       return 0
  elseif n<0    return E(4)*ER(-n)
  elseif n%4==1 return sum(k->E(n,k^2),1:n)
  elseif n%4==2 return (E(8)-E(8,3))*ER(div(n,2))
  elseif n%4==3 return -E(4)*sum(k->E(n,k^2),1:n)
  else          return 2*ER(div(n,4))
  end
end 

Base.Complex(c::Cyc)=iszero(c) ? Complex(0.0) : sum(v*exp(2*k*im*pi/c.n) for (k,v) in c.d)

function Base.Real(c::Cyc{T}) where T
  if iszero(c) return zero(T) end
  if c.n==1 return c.d[1][2] end
  if c!=conj(c) error("$c is not real") end
  return real(Complex(c))
end

struct Root1
 r::Rational{Int}
 Root1(r::Rational)=new(mod(r.num,r.den)//r.den)
end

function Root1(c::Cyc)
  if !(all(x->x[2]==1,c.d) || all(x->x[2]==-1,c.d))
    return nothing
  end
  for i in 0:c.n-1
    if c==E(c.n,i) return Root1(i//c.n) end
    if -c==E(c.n,i) return Root1(1//2+i//c.n) end
  end
  return nothing
end

function Root1(c::Real)
  if c==1 Root1(0//1)
  elseif c==-1 Root1(1//2)
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

Cycs.E(a::Root1)=E(conductor(a),exponent(a))


"""
  quadratic(c::Cyc) determines if c lives in a quadratic extension of Q
  it  returns a named  tuple (a=a,b=b,root=root,d=d) of  integers such that
  c=(a + b ER(root))//d or nothing if no such tuple exists

# Examples
```julia-repl
julia> quadratic(1+E(3))
(a = 1, b = 1, root = -3, den = 2)

julia> quadratic(1+E(5))

```
"""
quadratic=function(cyc::Cyc)
  den=lcm(denominator.(map(x->x[2],cyc.d)))::Int
  if den!=1 cyc=Cyc(cyc.n,[k=>Int(v*den) for (k,v) in cyc.d]) end
  l=fill(0,cyc.n)
  for (k,v) in cyc.d l[k+1]=v end
  if cyc.n==1 return (a=l[1],b=0,root=1,den=den) end

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
    gal=
    let cyc=cyc
      Set(map(i->galois(cyc,i),prime_residues(cyc.n)))
    end
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
  return (a=a,b=b,root=root,den=den*d)
end

function testmat(p) # testmat(12)^2 takes 0.27s in 1.0
  ss=vcat([[[i,j] for j in i+1:p-1] for i in 0:p-1]...)
  [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
end
end
