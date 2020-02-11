"""
What   is  implemented  here  is  "Puiseux  polynomials",  that  is  linear
combinations  of monomials  of the  type `x₁^{a₁}…  xₙ^{aₙ}` where `xᵢ` are
variables  and `aᵢ` are exponents which  can be arbitrary rational numbers.
Some  functions  described  below  need  their  argument  to  involve  only
variables  to integral  powers; we  will refer  to such objects as "Laurent
polynomials"; some functions require further that variables are raised only
to positive powers: we refer then to "true polynomials".

`@Mvp x₁,…,xₙ`

declares   that  `xᵢ`are  indeterminates  suitable  to  build  multivariate
polynomials.

```julia-repl
julia> @Mvp x,y

julia> (x+y)^3
Mvp{Int64}: x³+3x²y+3xy²+y³
```
`Mvp(p)` converts  the `Pol`  `p` to  an  `Mvp`. 

```julia-repl
julia> Pol(:q)
Pol{Int64}: q

julia> Mvp(q^2+q)
Mvp{Int64}: q²+q
```

`Mvp(x::Number)`

Returns the  constant multivariate polynomial whose constant term is `x`.

```julia-repl
julia> degree(Mvp(1))
0
```

One can divide an `Mvp` by another when the divison is exact
(this is equivalent to `ExactDiv`, see below).

```julia-repl
julia> (x^2-y^2)//(x-y)
Mvp{Int64}: x+y
```

Only monomials can be raised to a non-integral power; they can be raised
to  a fractional  power of  denominator  'b' only  if 'GetRoot(x,b)'  is
defined  where 'x'  is  their  leading coefficient.  For  an 'Mvp'  <m>,
the  function  'GetRoot(m,n)' is  equivalent  to  'm^(1/n)'. Raising  a
non-monomial Laurent polynomial  to a negative power  returns a rational
fraction.

|    gap> (2*x)^(1/2);
    ER(2)x^(1/2)
    gap> (evalf(2)*x)^(1/2);
    1.4142135624x^(1/2)
    gap> GetRoot(evalf(2)*x,2);
    1.4142135624x^(1/2)|

"""
module Mvps
# to use as a stand-alone module uncomment the next line
# export degree, root, coefficients, valuation
export Mvp, Monomial, @Mvp, variables, value
# benchmark: (x+y+z)^3     2.7μs 141 alloc
using ..ModuleElts: ModuleElt, ModuleElts
using ..Util: fromTeX
using ..Cycs: Cyc
#------------------ Monomials ---------------------------------------------
PowType=Int # could be int8 to save space if limiting degree
struct Monomial
  d::ModuleElt{Symbol,PowType}   
end
Monomial(a::Pair...)=Monomial(ModuleElt(convert.(Pair{Symbol,PowType},collect(a))))

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
Base.one(::Type{Monomial})=Monomial(ModuleElt(Pair{Symbol,PowType}[]))
Base.one(m::Monomial)=Monomial(zero(m.d))
Base.inv(a::Monomial)=Monomial(-a.d)
Base.div(a::Monomial, b::Monomial)=a*inv(b)
Base.:^(x::Monomial, p)= iszero(p) ? one(x) : Monomial(x.d*p)
@inline Base.getindex(a::Monomial,k)=getindex(a.d,k)

function Base.show(io::IO,m::Monomial)
  if isone(m) return end
  for (v,d) in m.d
    print(io,v)
    if haskey(fractional,v) d//=fractional[v] end
    if !isone(d) print(io,fromTeX(io,"^{$d}")) end
  end
end

function Base.show(io::IO, ::MIME"text/plain", m::Monomial)
  print(io,typeof(m),":")
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

Base.hash(a::Monomial, h::UInt)=hash(a.d,h)

degree(m::Monomial)=isone(m) ? 0 : sum(last.(m.d))

function degree(m::Monomial,var::Symbol)
  for (v,d) in m.d
    if v>var return 0 end
    if v==var return d end
  end
  return 0
end

function root(m::Monomial,n::Integer=2)
 if !all(x->iszero(x%n),last.(m.d)) throw(InexactError(:root,Monomial,n)) end
 Monomial((k=>div(v,n) for (k,v) in m.d)...)
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

Base.broadcastable(p::Mvp)=Ref(p)
Base.cmp(a::Mvp,b::Mvp)=cmp(a.d,b.d)
Base.isless(a::Mvp,b::Mvp)=cmp(a,b)==-1
Base.hash(a::Mvp,h::UInt)=hash(a.d,h)

function Base.show(io::IO, ::MIME"text/plain", a::Mvp)
  print(io,typeof(a),": ")
  show(io,a)
end

# necessary to clean up a ModuleElt upstream
Base.show(io::IO, x::Mvp)=show(IOContext(io,:showbasis=>nothing),x.d)

Base.zero(p::Mvp)=Mvp(zero(p.d))
Base.zero(::Type{Mvp})=Mvp(zero(ModuleElt{Monomial,Int}))
Base.zero(::Type{Mvp{T}}) where T=Mvp(zero(ModuleElt{Monomial,T}))
Base.one(::Type{Mvp{T}}) where T=Mvp(one(T))
Base.one(::Type{Mvp})=Mvp(1)
Base.one(p::Mvp{T}) where T=Mvp(one(T))
Base.copy(p::Mvp)=Mvp(p.d)
Base.iszero(p::Mvp)=length(p.d)==0
Base.adjoint(a::Mvp)=a
Base.transpose(p::Mvp)=p
Base.abs(p::Mvp)=p
Base.convert(::Type{Mvp},a::Number)=iszero(a) ? zero(Mvp{typeof(a)}) : 
                                          Mvp(one(Monomial)=>a)
Base.convert(::Type{Mvp{T}},a::Number) where T=iszero(a) ? zero(Mvp{T}) : 
                                          Mvp(one(Monomial)=>T(a))
(::Type{Mvp{T}})(a::Number) where T=convert(Mvp{T},a)
(::Type{Mvp{T}})(a::Mvp) where T=convert(Mvp{T},a)
Base.convert(::Type{Mvp{T}},a::Mvp{T1}) where {T,T1}=T==T1 ? a : iszero(a) ? zero(Mvp{T}) : Mvp([k=>T(v) for (k,v) in a.d]...)
Mvp(a::Number)=convert(Mvp,a)
Base.:(==)(a::Mvp, b::Mvp)=a.d==b.d

function Base.convert(::Type{T},a::Mvp) where T<:Number
  if iszero(a) return zero(T) end
  if length(a.d)>1 || !isone(a.d.d[1].first)
      throw(InexactError(:convert,T,a)) 
  end
  convert(T,a.d.d[1].second)
end
(::Type{T})(a::Mvp) where T<: Number=convert(T,a)

Base.isinteger(p::Mvp)=iszero(p) || (isone(length(p.d.d)) &&
                           isone(p.d.d[1].first) && isinteger(p.d.d[1].second))

function Base.promote(a::Mvp{T1},b::Mvp{T2}) where {T1,T2}
  T=promote_type(T1,T2)
  let T=T, a=a, b=b
   if T!=T1 a=Mvp(ModuleElt(Pair{Monomial,T}[m=>T(c) for (m,c) in a.d])) end
   if T!=T2 b=Mvp(ModuleElt(Pair{Monomial,T}[m=>T(c) for (m,c) in b.d])) end
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
Base.:-(b::Number, a::Mvp)=Mvp(b)-a
Base.:*(a::Monomial, b::Mvp)=Mvp(ModuleElt(m*a=>c for (m,c) in b.d))
Base.:*(a::Mvp, b::Mvp)=iszero(a) ? zero(b) : sum(Mvp(ModuleElt(m*m1=>c*c1 for (m1,c1) in b.d)) for (m,c) in a.d)
Base.:*(a::Number, b::Mvp)=Mvp(b.d*a)
Base.:*(b::Mvp, a::Number)=a*b
Base.:(//)(a::Mvp, b)=Mvp(ModuleElt(m=>c//b for (m,c) in a.d))
Base.conj(a::Mvp)=Mvp(map(x->(x[1]=>conj(x[2])),a.d)...)

function Base.:^(x::Mvp, p::Real)
  if iszero(x) return x end
  p=Int(p)
  if iszero(p) return one(x) end
  if isone(p) return x end
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

"""
The `degree` of a monomial is the sum of  the exponent of the variables.
The `degree` of an `Mvp` is the largest degree of a monomial.

```julia-repl
julia> a=x^2+x*y
Mvp{Int64}: x²+xy

julia> degree(a)
2
```

There  exists also a form of `degree`  taking as second argument a variable
name, which returns the degree of the polynomial in that variable.

```julia-repl
julia> degree(a,:y)
1

julia> degree(a,:x)
2
```

"""
degree(m::Mvp)=maximum(degree.(keys(m.d)))
degree(m::Mvp,v::Symbol)=maximum(map(x->degree(x,v),keys(m.d)))

"""
The `valuation` of an `Mvp` is the minimal degree of a monomial.


```julia-repl
julia> a=x^2+x*y
Mvp{Int64}: x²+xy

julia> valuation(a)
2
```

There  exists  also  a  form  of  `valuation`  taking  as second argument a
variable  name,  which  returns  the  valuation  of  the polynomial in that
variable.

```julia-repl
julia> valuation(a,:y)
0

julia> valuation(a,:x)
1
```

"""
valuation(m::Mvp)=minimum(map(degree,keys(m.d)))
valuation(m::Mvp,v::Symbol)=minimum(map(x->degree(x,v),keys(m.d)))

"""
  `coefficients(p::Mvp, var::Symbol)` 

returns as a Dict the list of coefficients of `p` with respect to `var`.

```julia-repl
julia> p=(x+y+inv(y))^4
Mvp{Int64}: x⁴+4x³y+4x³y⁻¹+6x²y²+12x²+6x²y⁻²+4xy³+12xy+12xy⁻¹+4xy⁻³+y⁴+4y²+6+4y⁻²+y⁻⁴

julia> coefficients(p,:x)
Dict{Int64,Mvp{Int64}} with 5 entries:
  0 => y⁴+4y²+6+4y⁻²+y⁻⁴
  4 => 1
  2 => 6y²+12+6y⁻²
  3 => 4y+4y⁻¹
  1 => 4y³+12y+12y⁻¹+4y⁻³

julia> coefficients(p,:y)
Dict{Int64,Mvp{Int64}} with 9 entries:
  0  => x⁴+12x²+6
  4  => 1
  -4 => 1
  -3 => 4x
  2  => 6x²+4
  -2 => 6x²+4
  -1 => 4x³+12x
  3  => 4x
  1  => 4x³+12x
```

The  same caveat is  applicable to 'coefficients'  as to values: the values
are  always `Mvp`s.  To get  a list  of scalars  for univariate polynomials
represented as `Mvp`s, one should use `ScalMvp`.
"""
function coefficients(p::Mvp,v::Symbol)
  if iszero(p) return Dict{PowType,typeof(p.d)}() end
  d=Dict{PowType,typeof(p.d)}()
  for (m,c) in p.d
#   print("$m=>$c d=$d\n")
    u=ModuleElts.drop(m.d,v)
    if isnothing(u)
      d[0]=push!(get(d,0,zero(p.d)),m=>c)
    else 
      (m1,deg)=u
      d[deg]=push!(get(d,deg,zero(p.d)),Monomial(m1)=>c)
    end
  end
  Dict(dg=>Mvp(c) for (dg,c) in d)
end

"""
`variables(p)`

returns the list of variables of the `Mvp` as a sorted list of strings.

```julia-repl
julia> variables(x+x^4+y)
2-element Array{Symbol,1}:
 :x
 :y
```
"""
function variables(p::Mvp)
  l=map(c->map(x->x[1],c[1].d),p.d)
  if length(l)==1 l[1]
  elseif  length(l)==0 Symbol[]
  else sort(union(l...))
  end
end

function scal(p::Mvp{T})where T
  if iszero(p) return zero(T) end
  if length(p.d)!=1 return nothing end
  (m,c)=first(p.d)
  if isone(m) return c end
  return nothing
end

"""
Value of an `Mvp`

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
but  a constant `Mvp`  (for consistency). See  the function 'ScalMvp' below
for how to convert such constants to their base ring.

```julia-repl
julia> p(x=y)
Mvp{Int64}: 7y⁴-2

julia> p(x=y,y=x)
Mvp{Int64}: 7x⁴-2
```
    gap> Value(p,["x",y,"y",x]);
    -2+7x^-1y^5|

Evaluating an  'Mvp' which is  a Puiseux  polynomial may cause  calls to
'GetRoot'

|    gap> p:=x^(1/2)*y^(1/3);
    x^(1/2)y^(1/3)
    gap> Value(p,["x",y]);
    y^(5/6)
    gap>  Value(p,["x",2]);
    ER(2)y^(1/3)
    gap>  Value(p,["y",2]);
    Error, : unable to compute 3-th root of 2
     in
    GetRoot( values[i], d[i] ) called from
    f.operations.Value( f, x ) called from
    Value( p, [ "y", 2 ] ) called from
    main loop
    brk>|
"""
# generic version:
#function value(p::Mvp,vv::Pair{Symbol,<:Mvp})
#  (s,v)=vv
#  res1=Tuple{typeof(v),Monomial}[]
#  for (m,c) in p.d
#    u=ModuleElts.drop(m.d,s)
#    if !isnothing(u)
#      (m1,deg)=u
#      push!(res1,(c*(Monomial(m1)*v^deg),m))
#    end
#  end
#  if isempty(res1) return p end
#  newd=p.d
#  newp=zero(p)
#  for (p,m) in res1
#    newd=ModuleElts.drop(newd,m)[1]
#    newp+=p
#  end
#  Mvp(newd)+newp
#end
function value(p::Mvp,vv::Pair)
  (s,v)=vv
  if !(v isa Mvp) v=Mvp(v) end
  res1=Tuple{typeof(v),Int}[]
  for i in eachindex(p.d.d)
    (m,c)=p.d.d[i]
    u=ModuleElts.drop(m.d,s)
    if !isnothing(u)
      push!(res1,(c*(Monomial(first(u))*v^last(u)),i))
    end
  end
  if isempty(res1) return p end
  newd=deleteat!(copy(p.d.d),map(last,res1))
  Mvp(ModuleElt(newd))+
     Mvp(ModuleElt(vcat(map(x->first(x).d.d,res1)...);check=true))
end

function (p::Mvp)(;arg...)
  for vv in arg p=value(p,vv) end
  p
end

function root(p::Mvp,n::Real=2)
  n=Int(n)
  if length(p.d)>1 throw(InexactError(:root,n,p)) end
  p=p.d.d[1]
  Mvp(root(first(p),n)=>root(last(p),n))
end

function Base.inv(p::Mvp)
  if length(p.d)>1 throw(InexactError(:inv,Mvp,p)) end
  (m,c)=first(p.d)
  if c^2==1 return Mvp(inv(m)=>c) end
  if (c isa Cyc) || (c isa Rational) return Mvp(inv(m)=>inv(c)) end
  throw(InexactError(:inv,typeof(c),p))
end

Base.:/(p::Mvp,q::Mvp)=p* inv(q)

function Base.://(p::Mvp,q::Mvp)
  res=ExactDiv(p,q)
  if isnothing(res) return p*inv(q) end
  res
end

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
    if isnothing(t) return nothing end
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

"""
The  function 'Derivative(p,v)' returns the  derivative of 'p' with respect
to  the variable given by the string 'v'; if 'v' is not given, with respect
to the first variable in alphabetical order.

|    gap>  p:=7*x^5*y^-1-2;
    -2+7x^5y^-1
    gap> Derivative(p,"x");
    35x^4y^-1
    gap> Derivative(p,"y");
    -7x^5y^-2
    gap> Derivative(p);
    35x^4y^-1
    gap>  p:=x^(1/2)*y^(1/3);
    x^(1/2)y^(1/3)
    gap>  Derivative(p,"x");
    1/2x^(-1/2)y^(1/3)
    gap>  Derivative(p,"y");
    1/3x^(1/2)y^(-2/3)
    gap>  Derivative(p,"z");
    0|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Finally we  mention the  functions 'ComplexConjugate' and  'evalf' which
are defined using  for coefficients the 'Complex'  and 'Decimal' numbers
of the CHEVIE package.

|    gap> p:=E(3)*x+E(5);
    E5+E3x
    gap> evalf(p);
    0.3090169944+0.9510565163I+(-0.5+0.8660254038I)x
    gap> p:=E(3)*x+E(5);          
    E5+E3x
    gap> ComplexConjugate(p);
    E5^4+E3^2x
    gap> evalf(p);
    0.3090169944+0.9510565163I+(-0.5+0.8660254038I)x
    gap> ComplexConjugate(last);
    0.3090169944-0.9510565163I+(-0.5-0.8660254038I)x|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'ScalMvp( <p> )'

If  <p> is  an  'Mvp' then  if  <p>  is a  scalar,  return that  scalar,
otherwise return  'false'. Or  if <p>  is a  list, then  apply 'ScalMvp'
recursively to  it (but return false  if it contains any  'Mvp' which is
not a scalar). Else assume <p> is already a scalar and thus return <p>.

|    gap> v:=[Mvp("x"),Mvp("y")];        
    [ x, y ]
    gap> ScalMvp(v);
    false
    gap> w:=List(v,p->Value(p,["x",2,"y",3]));
    [ 2, 3 ]
    gap> Gcd(w);
    Error, sorry, the elements of <arg> lie in no common ring domain in
    Domain( arg[1] ) called from
    DefaultRing( ns ) called from
    Gcd( w ) called from
    main loop
    brk> 
    gap> Gcd(ScalMvp(w));
    1|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'LaurentDenominator( <p1>, <p2>, ... )'

Returns the unique monomial 'm' of minimal degree such that for all the
Laurent polynomial arguments <p1>, <p2>, etc... the product `m* p_i` is
a true polynomial.

|    gap> LaurentDenominator(x^-1,y^-2+x^4);
    xy^2|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'OnPolynomials( <m>, <p> [,<vars>] )'

Implements  the action of  a matrix on  'Mvp's. <vars> should  be a list of
strings representing variables. If `v`'=List(vars,Mvp)', the polynomial `p`
is  changed  by  replacing  in  it  `v_i`  by `(v×m)_i`. If <vars> is
omitted, it is taken to be 'Variables(p)'.

|    gap> OnPolynomials([[1,2],[3,1]],x+y);    
    3x+4y|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'FactorizeQuadraticForm( <p>)'

<p>  should be an 'Mvp' of degree  2 which represents a quadratic form. The
function  returns a list of two linear forms of which <p> is the product if
such  forms exist, otherwise it returns 'false' (it returns [Mvp(1),<p>] if
<p> is of degree 1).

|    gap> FactorizeQuadraticForm(x^2-y^2+x+3*y-2);
    [ -1+x+y, 2+x-y ]
    gap> FactorizeQuadraticForm(x^2+x+1);        
    [ -E3+x, -E3^2+x ]
    gap> FactorizeQuadraticForm(x*y+1);  
    false|

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
end
