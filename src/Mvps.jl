module Mvps
using Gapjm
export Mvp, @Mvp, variables
# benchmark: (x+y+z)^3     2.7μs 141 alloc
#------------------ Monomials ---------------------------------------------
PowType=Int # could be int8 to save space if limiting degree
struct Monomial
  d::ModuleElt{Symbol,PowType}   
end
Monomial(a::Pair...)=Monomial(ModuleElt{Symbol,PowType}(collect(a)))

const fractional=Dict{Symbol,PowType}() 
# if fractional(:x)==y then x interpreted as x^(1/y)

function Base.convert(::Type{Monomial},v::Symbol)
  if haskey(fractional,v) pow=fractional[v]
  else pow=1
  end
  Monomial(v=>PowType(pow))
end
Monomial(v::Symbol)=convert(Monomial,v)

Base.:*(a::Monomial, b::Monomial)=Monomial(a.d+b.d)
Base.isone(a::Monomial)=iszero(a.d)
Base.iszero(a::Monomial)=false
Base.one(::Type{Monomial})=Monomial()
Base.one(m::Monomial)=Monomial(zero(m.d))
Base.inv(a::Monomial)=Monomial(-a.d)
Base.div(a::Monomial, b::Monomial)=a*inv(b)
Base.:^(x::Monomial, p)= iszero(p) ? one(x) : Monomial(x.d*p)
@inline function ModuleElts.drop(a::Monomial,k)
   u=ModuleElts.drop(a.d,k)
   if isnothing(u) return u end
   Monomial(u[1]),u[2]
end
@inline Base.getindex(a::Monomial,k)=getindex(a.d,k)

function Base.show(io::IO, m::Monomial)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if isone(m) print(io,"1") 
  else
    for (v,d) in m.d
      print(io,v)
      if haskey(fractional,v) d//=fractional[v] end
      if !isone(d)
        if repl || TeX  
          e=1<d<10 ? "^$d" : "^{$d}"
          print(io, TeX ? e : TeXstrip(e)) 
        else print(io,"^$d") 
        end
      end
    end
  end
end

function Base.show(io::IO, ::MIME"text/plain", m::Monomial)
  print(io,"Monomial:")
  show(io,m)
end

# cmp must define a monomial order (a<b => a m< b m for all monomials a, b, m)
# We take the sign of the power of the first variable in a/b
function Base.cmp(a::Monomial, b::Monomial)
  for (ea,eb) in zip(a.d,b.d)
    c=cmp(ea[1],eb[1])
    if c<0 return -sign(ea[2])
    elseif c>0 return sign(eb[2])
    end
    c=cmp(ea[2],eb[2])
    if c!=0 return -c end
  end
  if length(a.d)>length(b.d) return -sign(a.d.d[length(b.d)+1][2])
  elseif length(a.d)<length(b.d) return sign(b.d.d[length(a.d)+1][2])
  else return 0
  end
end

Base.isless(a::Monomial, b::Monomial)=cmp(a,b)==-1
Base.:(==)(a::Monomial, b::Monomial)=cmp(a,b)==0
#Base.:(==)(a::Monomial, b::Monomial)=a.d==b.d

# gcd(m_1,...,m_k)= largest m such that m_i>=m
function Base.gcd(l::Monomial...)
  if length(l)==0 return one(Monomial)
  elseif length(l)==1 return l[1]
  elseif length(l)>2 return reduce(gcd,l)
  else (a,b)=l
  end
  res=empty(a.d)
  ai=bi=1
  la=length(a.d)
  lb=length(b.d)
  while ai<=la && bi<=lb
    if a.d[ai][1]>b.d[bi][1] bi+=1
    elseif a.d[ai][1]<b.d[bi][1] ai+=1
    else 
      if    a.d[ai][2]<=b.d[bi][2] push!(res,a.d[ai]) 
      elseif a.d[ai][2]>b.d[bi][2] push!(res,b.d[bi]) 
      end
      ai+=1
      bi+=1
    end
  end
  Monomial(res)
end

function Base.hash(a::Monomial, h::UInt)
   b = 0x595dee0e71d271d0%UInt
   for (s,p) in a.d
     b = xor(b,xor(hash(s, h),h))
     b = xor(b,xor(hash(p, h),h))
     b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

Gapjm.degree(m::Monomial)=isone(m) ? 0 : sum(last.(m.d))

function Gapjm.degree(m::Monomial,var::Symbol)
  for (v,d) in m.d
    if v>var return 0 end
    if v==var return d end
  end
  return 0
end
#------------------------------------------------------------------------------
struct Mvp{T} # N=type of exponents T=type of coeffs
  d::ModuleElt{Monomial,T}
end

Mvp(a::Pair...)=Mvp(ModuleElt(a...))
Mvp(v::Symbol)=Mvp(Monomial(v)=>1)
Mvp{T}(n::T) where T=Mvp(one(Monomial)=>n)

macro Mvp(t) # Mvp x,y,z defines variables to be Mvp
  if t isa Expr
    for v in t.args
      Base.eval(Main,:($v=Mvp($(Core.QuoteNode(Symbol(v))))))
    end
  elseif t isa Symbol
    Base.eval(Main,:($t=Mvp($(Core.QuoteNode(t)))))
  end
end

Base.show(io::IO, x::Mvp)=ModuleElts.helpshow(io,x.d;showbasis=show)

Base.zero(p::Mvp)=Mvp(zero(p.d))
Base.zero(::Type{Mvp})=Mvp(ModuleElt{Monomial,Int}())
Base.zero(::Type{Mvp{T}}) where T=Mvp(ModuleElt{Monomial,T}())
Base.one(::Type{Mvp{T}}) where T=Mvp(one(T))
Base.one(p::Mvp{T}) where T=Mvp(one(T))
Base.copy(p::Mvp)=Mvp(p.d)
Base.iszero(p::Mvp)=length(p.d)==0
Base.transpose(p::Mvp)=p
Base.convert(::Type{Mvp},a::Number)=iszero(a) ? zero(Mvp{typeof(a)}) : 
                                          Mvp(one(Monomial)=>a)
Base.convert(::Type{Mvp{T}},a::Number) where T=Mvp(one(Monomial)=>T(a))
Mvp(a::Number)=convert(Mvp,a)
Base.:(==)(a::Mvp, b::Mvp)=a.d==b.d

function Base.promote(a::Mvp{T1},b::Mvp{T2}) where {T1,T2}
  T=promote_type(T1,T2)
  let T=T, a=a, b=b
    if T!=T1 a=Mvp(ModuleElt{Monomial,T}(m=>T(c) for (m,c) in a.d)) end
    if T!=T2 b=Mvp(ModuleElt{Monomial,T}(m=>T(c) for (m,c) in b.d)) end
    a,b
  end
end

function Base.:+(a::Mvp, b::Mvp)
  a,b=promote(a,b)
  Mvp(a.d+b.d)
end

Base.:+(a::Number, b::Mvp)=Mvp(a)+b
Base.:+(a::Mvp, b::Number)=b+a
Base.:-(a::Mvp)=Mvp(-a.d)
Base.:-(a::Mvp, b::Mvp)=a+(-b)
Base.:-(a::Mvp, b::Number)=a-Mvp(b)
Base.:-(b::Number, a::Mvp)=Mvp(a)-b
Base.:*(a::Monomial, b::Mvp)=Mvp(ModuleElt(m*a=>c for (m,c) in b.d))
Base.:*(a::Mvp, b::Mvp)=sum(Mvp(ModuleElt(m*m1=>c*c1 for (m1,c1) in b.d)) for (m,c) in a.d)
Base.:*(a, b::Mvp)=Mvp(b.d*a)
Base.:*(b::Mvp, a)=a*b

function Base.inv(x::Mvp)
  if length(x.d)!=1 error("can only take inverse of monomial") end
  (m,c)=first(x.d)
  Mvp(inv(m)=>(c^2==1 ? c : 1//c))
end

function Base.:^(x::Mvp, p::Int)
  if iszero(x) return 0 end
  if p<0
    x=inv(x)
    p=-p
  end
  if length(x.d)==1
    (m,c)=first(x.d)
    return Mvp(m^p=>c^p)
  end
  Base.power_by_squaring(x,p)
end

Gapjm.degree(m::Mvp)=maximum(map(degree,keys(m.d)))
Pols.valuation(m::Mvp)=minimum(map(degree,keys(m.d)))

function coefficients(p::Mvp,v::Symbol)
  if iszero(p) return Dict() end
  d=Dict{PowType,typeof(p.d)}()
  for (m,c) in p.d
#   print("$m=>$c d=$d\n")
    u=ModuleElts.drop(m,v)
    if isnothing(u)
      d[0]=push!(get(d,0,zero(p.d)),m=>c)
    else 
      (m1,deg)=u
      d[deg]=push!(get(d,deg,zero(p.d)),m1=>c)
    end
  end
  Dict(dg=>Mvp(c) for (dg,c) in d)
end

variables(p::Mvp)=sort(union(map(c->map(x->x[1],c[1].d),p.d)...))

function scal(p::Mvp{T})where T
  if iszero(p) return zero(T) end
  if length(p.d)!=1 return nothing end
  (m,c)=first(p.d)
  if isone(m) return c end
  return nothing
end

function (p::Mvp)(;arg...)
  res=p
  for s in keys(arg)
    res1=zero(res.d)
    for (m,c) in res.d
      u=ModuleElts.drop(m,s)
      if isnothing(u) push!(res1,m=>c)
      else 
        (m1,deg)=u
        append!(res1,(Mvp(m1=>c)*arg[s]^deg).d)
      end
    end
    res=Mvp(norm!(res1))
  end
  res
end

function Base.://(p::Mvp,q::Mvp)
  res=ExactDiv(p,q)
  if isnothing(res) error("rational fractions of Mvp not yet implemented") end
  res
end

function ExactDiv(p::Mvp,q::Mvp)
# println("p=$p q=$q")
  if length(q.d)==1 return p*inv(q)
  elseif length(p.d)<2 return iszero(p) ? p : nothing
  elseif iszero(q) error("cannot divide by 0")
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
    if isnothing(t) return nothing end
    if mp!=mq t*=Monomial(var=>mp-mq) end
    res+=t
#   print("t=$t res=$res p=$p=>")
    p-=t*q
#   println("$p")
  end
  res
end


function Base.gcd(a::Mvp,b::Mvp)
  if any(x->length(x.d)==1,[a,b]) 
    return Mvp([gcd(vcat(map(x->x[1],a.d),map(x->x[1],b.d))...)=>1])
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
    plp=ScalMvp(p[lp])
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
     if ScalMvp(cont[i]==false) coef[i]=List(coef[i],x->ExactDiv(x,cont[i]))
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
end
