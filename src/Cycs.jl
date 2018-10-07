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

Contrary   to  GAP,  I  do  not  lower  automatically  numbers  after  each
computation since this is too expensive with this code. GAP also converts a
Cyclotomic which is rational to a Rational, a Rational which is integral to
an  Int, etc... This is tremendously useful  but needs a new type of number
to be added to Julia, which requires more competent people than me.

The main way to build a Cyclotomic number is to use the function `E(n,k=1)`
which constructs ζ_n^k.

# Examples
```julia-repl
julia> E(3)+E(4)
E(12)^4-E(12)^7-E(12)^11

julia> E(3,2)
E(3)^2

julia> 1+E(3,2)
-E(3)

julia> a=E(4)-E(4)
0

julia> conductor(a) # a not lowered so still in Q(ζ_4)
4

julia> conductor(lower(a)) # but now a is lowered to Q(ζ_1)=Q
1

julia> typeof(convert(Int,a))
Int64

julia> convert(Int,E(4))
ERROR: InexactError: convert(Int64, E(4))

julia> inv(1+E(4)) # inverses need Rationals
(1/2)+(-1/2)E(4)

julia> typeof(ans)
Cyc{Rational{Int64}}

julia> typeof(1+E(4))
Cyc{Int64}

ulia> Cyc(1+im) # one can convert Gaussian integers or rationals
1+E(4)

julia> 1//(1+E(4))
(1/2)+(-1/2)E(4)

julia> typeof(Cyc(1//2)) # another way of building a Cyc
Cyc{Rational{Int64}}

julia> conj(1+E(4))
1-E(4)

julia> E(9)     # an effect of the Zumbroich basis
-E(9)^4-E(9)^7

julia> AsRootOfUnity(ans) # but you can decide if a Cyc is a root of unity
1/9

julia> Complex(E(3))   # convert to float is probably not very useful
-0.4999999999999998 + 0.8660254037844387im

julia> Cyc(ans) # even less useful
-0.4999999999999998+0.8660254037844387E(4)
```

For more information look at the documentation of ER, quadratic, galois. 
"""
module Cycs
export E, ER, Cyc, conductor, lower, galois, AsRootOfUnity, quadratic
using Memoize,Util

struct Cyc{T <: Real}<: Number  # a cyclotomic number
  n::Int
  c::Vector{T} # coefficients on the Zumbroich basis
               # the i-th element is the coefficient on zumbasis[i]
end

conductor(c::Cyc)=c.n
conductor(a::Array)=lcm(conductor.(a))
conductor(i::Int)=1

#=
  zumbasis(n::Int) returns the Zumbroich basis of Q(ζ_n) 
  as the vector of i in 0:n-1 such that ``ζ_n^i`` is in the basis
=#
@memoize Dict function zumbasis(n::Int)::Vector{Int}
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
  Elist(n,i) returns the coefficients of E(n)^i in zumbasis(n).
  These  coeffs are all  1 or all  -1. Thus the  result is a Pair sgn=>inds
  where  sgn is true if  coeffs are all 1  and false otherwise, and inds is
  the list of indices in zumbasis of the nonzero coeffs.
=#
@memoize Dict function Elist(n::Int,i::Int=1)::Pair{Bool,Vector{Int}}
  function modp(n::Int,i::Int)::Vector{Int}
    nfact=factor(n)
    res=Int[]
    for (p,np) in nfact
      factor=p^np
      m=div(n,factor)
      cnt=mod(i*invmod(m,factor),factor)
      i-=cnt*m
      if p==2
        if 1==div(cnt,div(factor,p)) push!(res,p) end
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
        if tmp[1]==0 push!(res,p) end
      end
    end
    res
  end
  if n==1 return true=>[1] end
  z=zumbasis(n)::Vector{Int}
  i=mod(i,n)
  mp=modp(n,i)::Vector{Int}
  if isempty(mp) 
    return true=>Vector{Int}(indexin([i],z))
  end
  v=cartesian((1:p-1 for p in mp)...)*div.(n,mp)
  u=(i .+ v).%n
  (length(mp)%2==0)=>Vector{Int}(indexin(u,z))
end

Cyc(i::Real)=Cyc(1,[i])

function Cyc(c::Complex{<:Real})
   Cyc(4,[real(c),imag(c)])
end

function Cyc(n::Int,ind::Vector{Int},coeff)
  zb=zumbasis(n)
  c=Cyc(n,zeros(eltype(coeff),length(zb)))
  c.c[indexin(ind.%n,zb)]=coeff
  c
end

function Base.convert(::Type{T},c::Cyc)where T<:Real
  if iszero(c) return 0 end
  c=lower(c)
  if c.n!=1 throw(InexactError(:convert,T,c)) end
  convert(T,c.c[1])
end

function Base.convert(::Type{Cyc{T}},a::Real)where T<:Real
   Cyc(T(a))
end

function Base.promote_type(::Type{Cyc{T}},::Type{Cyc{T}})where T <: Real
   Cyc{T}
end

function Base.promote_type(::Type{Cyc{T}},::Type{Cyc{T1}})where T where T1 <: Real
   Cyc{promote_type(T,T1)}
end

function Base.promote_type(::Type{Cyc{T}},::Type{T1})where T where T1 <: Real
   Cyc{promote_type(T,T1)}
end

function Base.promote_type(::Type{T1},::Type{Cyc{T}})where T where T1 <: Real
   Cyc{promote_type(T,T1)}
end

function Base.convert(::Type{Cyc{T}},a::Cyc)where T<:Real
 if eltype(a.c)==T return a end
 Cyc(a.n,Vector{T}(a.c))
end

Base.zero(c::Cyc)=Cyc(c.n,fill!(similar(c.c), zero(eltype(c.c))))
Base.iszero(c::Cyc)=all(iszero,c.c)

Base.:(==)(a::Cyc,b::Cyc)=(a.n==b.n) && (a.c==b.c)

function Base.hash(a::Cyc, h::UInt)
   b = 0x595dee0e71d271d0%UInt
   for i in a.c
     b = xor(b,xor(hash(i, h),h))
     b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

function Base.show(io::IO,p::Cyc)
  p=lower(p)
  e="E($(p.n))"
  l=map(p.c,zumbasis(p.n))do v,deg
    if v==0 return "" end
    s="$v"
    if occursin(r"(//|/|\+|\-|\*)",s[2:end]) s="($s)" end
    if deg==0 t=s
    else t=(v==1 ? "" : (v==-1 ? "-" : s))*e*(deg==1 ? "" : "^$deg")
    end
    (t[1]!='-' ? "+" : "")*t
  end::Array{String,1}
  s=join(l)
  if s=="" print(io,"$(p.c[1])")
  elseif  s[1]=='+' print(io,s[2:end])
  else print(io,s) end
end

"""
  E(n::Integer,k::Integer=1) is exp(2i k π/n)
"""
function E(n::Integer,i::Integer=1)
  (s,t)=Elist(n,mod(i,n))
  c=Cyc(n,zero(zumbasis(n)))
  c.c[t].=ifelse(s,1,-1)
  c
end

function Base.promote(a::Cyc,b::Cyc)
  if a.n==b.n 
    ac,bc=promote(a.c,b.c)
    if typeof(ac)!=typeof(a.c) a=Cyc(a.n,ac) end
    if typeof(bc)!=typeof(b.c) b=Cyc(b.n,bc) end
    return (a,b) end
  l=lcm(a.n,b.n)
  return promote(raise(l,a),raise(l,b))
end

function Base.:+(a::Cyc,b::Cyc)
  (a,b)=promote(a,b)
  Cyc(a.n,a.c+b.c)
end

function Base.isless(a::Cyc,b::Cyc)
  (a,b)=promote(a,b)
  lexless(a.c,b.c)
end

Base.:-(a::Cyc)=Cyc(a.n,-a.c)
Base.:-(a::Cyc,b::Cyc)=a+(-b)

Base.:*(a::Real,c::Cyc)=Cyc(c.n,a*c.c)
Base.:*(c::Cyc,a::Real)=Cyc(c.n,a*c.c)

function addExb!(a::Cyc,i::Int,b)  # a+=b*E(a.n,i)
  s,t=Elist(a.n,i)
  if !s b=-b end
  @inbounds  a.c[t].+=b
end

function Base.:*(a::Cyc,b::Cyc)
  a,b=promote(a,b)
  res=zero(a)
  zb=zumbasis(a.n)
  for i in eachindex(a.c), j in eachindex(b.c) 
    if !iszero(a.c[i]) && !iszero(b.c[j])
      addExb!(res,zb[i]+zb[j],a.c[i]*b.c[j]) 
#    res+=a.c[i]*b.c[j]*E(res.n,zb[i]+zb[j])
    end
  end
  res
end

Base.div(c::Cyc,a::Real)=Cyc(c.n,div.(c.c,a))

function raise(n::Int,c::Cyc) # write c in Q(ζ_n) if c.n divides n
  if n==c.n return c end
  if n%c.n !=0 error("raise to $n not multiple of $(c.n)") end
  zn=zumbasis(n)
  res=Cyc(n,zeros(eltype(c.c),length(zn)))
  z=zumbasis(c.n)*div(n,c.n)
  for i in eachindex(c.c) 
@inbounds if !iszero(c.c[i]) addExb!(res,z[i],c.c[i]) end
  end
  res
end

function lower(c::Cyc) # write c in smallest Q(ζ_n) where it sits
  n=c.n
  nz=findall(x->!iszero(x),c.c)
  if length(nz)==0 return Cyc(c.c[1]) end
  zb=zumbasis(n)
  for (p,np) in factor(n)
    m=div(n,p)
    let p=p, m=m, zb=zb
    if np>1
      if all(z->z%p==0,zb[nz]) 
        return lower(Cyc(m,div.(zb[nz],p),c.c[nz]))
      end
    else
      res=groupby(x->x%m,zb[nz])
      if all(v->length(v)==p-1,values(res))
        kk=[div(k+m*mod(-k,p)*invmod(m,p),p)%m for k in keys(res)]
        if p==2 return lower(Cyc(m,kk,c.c[indexin((kk*p).%n,zb)]))
        elseif all(k->constant(c.c[indexin((m*(1:p-1).+k*p).%n,zb)]),kk)
          return lower(Cyc(m,kk,-c.c[indexin((m.+kk*p).%n,zb)]))
        end
      end
    end
    end
  end
  c
end

# write elts of v in smallest Q(ζ_n) where all elements sit
function lower(v::Array)
  v=lower.(v)
  raise(conductor(v),v)
end

function raise(n::Int,v::Array)
  raise.(n,v)
end

"""
  galois(c::Cyc,n::Int) applies to c the galois automorphism
  of Q(ζ_conductor(c)) raising all roots of unity to the n-th power.
  n should be prime to c.
# Examples
```julia-repl
julia> galois(1+E(4),-1) # galois(c,-1) is the same as conj(c)
1-E(4)

julia> galois(ER(5),2)==-ER(5)
true
```
"""
function galois(c::Cyc,n::Int)
  if gcd(n,c.n)!=1 error("$n should be prime to conductor($c)=$(c.n)") end
  zb=zumbasis(c.n)
  nz=findall(x->!iszero(x),c.c)
  if isempty(nz) c 
  else let zb=zb
     sum(t->c.c[t]*E(c.n,zb[t]*n),nz)
       end
  end
end

Base.conj(c::Cyc)=galois(c,-1)

function Base.inv(c::Cyc)
  if c.n==1 
    n=convert(Int,c)
    r=Cyc(1)
  else
    p=prime_residues(c.n)[2:end]
    r=prod(i->galois(c,i),p)
    n=convert(Int,c*r)
  end
  n==1 ? r : (n==-1 ? -r : r//n)
end

Base.://(a::Cyc,c::Cyc)=a*inv(c)
Base.://(c::Cyc,a::Real)=Cyc(c.n,c.c//a)
Base.://(a::Real,c::Cyc)=a*inv(c)

"""
  ER(n::Int) computes as a Cyc the square root of the integer n.
# Examples
```julia-repl
julia> ER(-3)
E(3)-E(3)^2

julia> ER(3)
-E(12)^7+E(12)^11
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

Base.Complex(c::Cyc)=sum(x->x[2]*exp(2*x[1]*im*pi/c.n),zip(zumbasis(c.n),c.c))

function AsRootOfUnity(c::Cyc)
  if !(all(x->x==0 || x==1,c.c) ||all(x->x==0 || x==-1,c.c))
    return nothing
  end
  for i in 0:c.n-1
    if c==E(c.n,i) return i//c.n end
  end
  return nothing
end

function AsRootOfUnity(c::Real)
  if c==1 0//1
  elseif c==-1 1//2
  else nothing
  end
end

"""
  quadratic(c::Cyc) determines if c lives in a quadratic extension of Q
  it returns tuple (a,b,root,d) of integers such that  c=(a + b ER(root))//d 
  or false if no such tuple exists
# Examples
```julia-repl
julia> quadratic(1+E(3))
(1, 1, -3, 2)

julia> quadratic(1+E(5))
false
```
"""
quadratic=function(cyc::Cyc)
  cyc=lower(cyc)
  den=lcm(denominator.(cyc.c))::Int
  zb=zumbasis(cyc.n)
  l=fill(0,cyc.n)
  l[1 .+ zb]=Int.(cyc.c*den)
  if cyc.n==1 return (l[1],0,1,den) end

  f=factor(cyc.n)
  v2=get(f,2,0)

  if v2>3 || (v2==2 && any(p->p[1]!=2 && p[2]!=1,f)) ||
     (v2<2 && any(x->x!=1,values(f))) 
    return false
  end

  f=keys(f)
  if v2==0
    root=cyc.n
    if root%4==3 root=-root end
    gal=
    let cyc=cyc
      Set(map(i->lower(galois(cyc,i)),prime_residues(cyc.n)))
    end
    if length(gal)!=2 return false end
    a=convert(Int,sum(gal))                # trace of 'cyc' over the rationals
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
  
  if den*d*cyc!=lower(a+b*ER(root)) return false end
  return (a,b,root,den*d)
end

function testmat(p) # testmat(12)^2 takes 0.31s in 1.0
  ss=vcat([[[i,j] for j in i+1:p-1] for i in 0:p-1]...)
  [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
end
# testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
#  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
#end; in GAP3 takes 0.4s for testmat(12)^2, 0.3s in GAP4
end
