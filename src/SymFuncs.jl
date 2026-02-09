"""
This  module  deals  with  symmetric  functions.  The  main  object  is the
(truncated) algebra of symmetric functions, and the bases `p` (power sums),
`h`  (complete symmetric  functions), `e`  (elementary symmetric functions)
and  `s` (Schur functions). Each of these bases in degree `n` is indexed by
the partitions of `n`.

A typical session would begin by defining the symmetric functions and these
bases:
```julia-repl
julia> A=SymFuncAlgebra(10) # truncated in degree >10
SymFuncAlgebra(10)

julia> p=pbasis(A);h=hbasis(A);e=ebasis(A);s=sbasis(A);
```
Then making elements in one of these bases. The following forms are equivalent:
```julia-repl
julia> p(Partition(3,2,1))
p₃₂₁

julia> p([3,2,1])
p₃₂₁

julia> p(3,2,1)
p₃₂₁
```
The functions can be used to convert between bases:
```julia-repl
julia> s(p(3,2,1))
-s₁₁₁₁₁₁-s₂₂₂+s₃₁₁₁+s₃₃-s₄₁₁+s₆

julia> h(p(3,2,1))
-h₁₁₁₁₁₁+5h₂₁₁₁₁-6h₂₂₁₁-3h₃₁₁₁+6h₃₂₁
```
The following operations are defined on symmetric functions, in addition to
`+` and `-`:
```julia-repl
julia> s(2,1)*s(2,1) # product
s₂₂₁₁+s₂₂₂+s₃₁₁₁+2s₃₂₁+s₃₃+s₄₁₁+s₄₂

julia> s(2,1)⊗s(2,1) # inner product
s₁₁₁+s₂₁+s₃

julia> scalar_product(p(1,1,1),s(2,1))
2//1
```
One can mix bases in these operations. The basis of the left argument wins:
```julia-repl
julia> s(2,1)+p(3)
s₁₁₁+s₃
```
The plethysm is implemented with two possible notations:

```julia-repl
julia> plethysm(p(2,1),s(2,1))
(1//9)p₂₂₂₁₁₁+(-1//9)p₃₂₂₂+(-1//9)p₆₁₁₁+(1//9)p₆₃

julia> p(2,1)[s(2,1)]  # the same thing
(1//9)p₂₂₂₁₁₁+(-1//9)p₃₂₂₂+(-1//9)p₆₁₁₁+(1//9)p₆₃

julia> @Mvp u,v

julia> p(2,1)[u*p(2)+v*p(3)] # plethysm acts on Mvp coefficients
u³p₄₂+u²vp₄₃+uv²p₆₂+v³p₆₃
```

finally one can convert a symmetric function to a symmetric polynomial.
The number of variables is the highest degree in the function.
```julia-repl
julia> Mvp(p(2)+p(3))
Mvp{Int64}: x₁³+x₁²+x₂³+x₂²+x₃³+x₃²

julia> Mvp(p(2),[u,v]) # one can choose the variables used
Mvp{Int64}: u²+v²
```
"""
module SymFuncs
using ..Chevie
export SymFuncAlgebra, sbasis, pbasis, hbasis, ebasis, plethysm

struct SymFuncAlgebra
  n::Int
  charvalues::Dict{Pair{Partition,Partition},Int}
  centralizers::Dict{Partition,Int}
  partitions::Dict{Int,Vector{Partition}}
end

function charvalue(A::SymFuncAlgebra,p::Pair{Partition,Partition})
  l=size(p[1])
  if l!=size(p[2]) || l>A.n return 0 end
  get!(A.charvalues,p)do 
    ct=CharTable(coxsym(l))
    pp=Partitions(A,l)
    for x in eachindex(pp), y in eachindex(pp) 
      A.charvalues[pp[x]=>pp[y]]=ct.irr[x,y]
    end
    for x in eachindex(pp) A.centralizers[pp[x]]=ct.centralizers[x] end
    A.charvalues[p]
  end
end

function Groups.centralizer(A::SymFuncAlgebra,p::Partition)
  if size(p)>A.n return 0 end
  get!(A.centralizers,p)do 
    charvalue(A,p=>p)
    A.centralizers[p]
  end
end

function Combinat.Partitions(A::SymFuncAlgebra,n)
  if n>A.n return Partition[] end
  get!(A.partitions,n)do
    Partition.(partitions(n))
  end
end

function Base.show(io::IO,A::SymFuncAlgebra)
  print(io,"SymFuncAlgebra(",A.n,")")
end

function H(lambda,s)
  n=sum(lambda)
  t=UnipotentValues(UnipotentClasses(rootdatum(:gl,n));classes=true)
  i=findfirst(==(lambda),partitions(n))
  sum(t.scalar[:,i].*s.(partitions(n)))
end

"`SymFuncAlgebra(n)` the algebra of symmetric functions truncated in `degree>n`"
SymFuncAlgebra(n)=SymFuncAlgebra(n,
   Dict{Pair{Partition,Partition},Int}((Partition()=>Partition())=>1),
   Dict{Partition,Int}(Partition()=>1),
   Dict{Int,Vector{Partition}}())
                                 
struct SymFunc{b,C,TS}# b=:s,:p or :h depending on the basis
  d::ModuleElt{Partition,C}
  A::TS
  SymFunc{b}(m,A)where b =new{b,valtype(m),typeof(A)}(m,A)
end

function Base.show(io::IO, h::SymFunc{b})where b
  function showbasis(io::IO,w)
    res=string(b)
    if hasdecor(io) && !get(io,:naive,false) res*="_{"*xrepr(io,w)*"}"
    else  res*="("*xrepr(io,w)*")"
    end
    fromTeX(io,res)
  end
  show(rio(io,limit=true,showbasis=showbasis),improve_type(h.d))
end

clone(h::SymFunc{b},d) where b=SymFunc{b}(d,h.A)

Base.:+(a::SymFunc{ba},b::SymFunc) where ba=clone(a,a.d+basis(Val(ba),b).d)
Base.:+(a::SymFunc, b::Union{Number,Pol,Mvp})=a+one(a)*b
Base.:-(a::SymFunc)=clone(a,-a.d)
Base.:-(a::SymFunc{ba},b::SymFunc) where ba=clone(a,a.d-basis(Val(ba),b).d)
Base.:*(a::SymFunc, b::Union{Number,Pol,Mvp,Frac})=clone(a,a.d*b)
Base.:(//)(a::SymFunc, b::Union{Number,Pol,Mvp})=clone(a,a.d//b)
Base.:*(b::Union{Number,Pol,Mvp,Frac}, a::SymFunc)=a*b

Base.:^(a::SymFunc, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                      Base.power_by_squaring(inv(a),-n)
Base.one(a::SymFunc)=clone(a,ModuleElt(Partition()=>1))
Base.zero(a::SymFunc)=clone(a,zero(a.d))
Base.zero(b::Symbol,a::SymFunc)=SymFunc{b}(zero(a.d),a.A)
sbasis(H::SymFuncAlgebra)=(x...)->basis(H,Val(:s),x...)
pbasis(H::SymFuncAlgebra)=(x...)->basis(H,Val(:p),x...)
hbasis(H::SymFuncAlgebra)=(x...)->basis(H,Val(:h),x...)
ebasis(H::SymFuncAlgebra)=(x...)->basis(H,Val(:e),x...)
sbasis(h::SymFunc)=basis(h.A,Val(:s),h)
pbasis(h::SymFunc)=basis(h.A,Val(:p),h)
hbasis(h::SymFunc)=basis(h.A,Val(:h),h)
ebasis(h::SymFunc)=basis(h.A,Val(:e),h)

Algebras.basis(h::SymFunc{b})where b=b
Algebras.basis(H::SymFuncAlgebra,::Val{b},p::Partition)where b=
  SymFunc{b}(ModuleElt(p=>1),H)
Algebras.basis(H::SymFuncAlgebra,::Val{b},w::Vector{<:Integer})where b=
  basis(H,Val(b),Partition(w))
Algebras.basis(H::SymFuncAlgebra,::Val{b},w::Vararg{Integer})where b=
  basis(H,Val(b),collect(w))
Algebras.basis(::Val{b},h::SymFunc)where b=basis(h.A,Val(b),h)

Algebras.basis(H::SymFuncAlgebra,::Val{b},h::SymFunc{b}) where b =h
function Algebras.basis(H::SymFuncAlgebra,::Val{b1},h::SymFunc{b}) where {b,b1}
  if b==:p error(":p to ",b1," not implemented") end
  basis(h.A,Val(b1),pbasis(h))
end

function Algebras.basis(H::SymFuncAlgebra,::Val{:p},h::SymFunc{:s})
   sum(h.d)do (p,c)
       SymFunc{:p}(ModuleElt(l=>charvalue(H,p=>l)*c//centralizer(H,l)
            for l in Partitions(H,size(p))),H)
   end
end
Algebras.basis(H::SymFuncAlgebra,::Val{:p},h::SymFunc{:h})=sum(((p,c),)->
 prod(i->pbasis(basis(H,Val(:s),i)),p.l;init=pbasis(H)())*c,h.d;init=zero(:p,h))
Algebras.basis(H::SymFuncAlgebra,::Val{:p},h::SymFunc{:e})=sum(((p,c),)->
   prod(i->pbasis(basis(H,Val(:s),fill(1,i))),p.l;init=pbasis(H)())*c,
                                                         h.d;init=zero(:p,h))

Algebras.basis(H::SymFuncAlgebra,::Val{:s},h::SymFunc{:p})=
  sum(h.d;init=zero(:s,h))do (p,c)
    SymFunc{:s}(ModuleElt(l=>charvalue(H,l=>p)*c for l in
                     Partitions(H,size(p))),H)
  end

function Pto(H,n,b) # p(n) to basis b∈(:e,:h)
  SymFunc{b}(ModuleElt(map(Partitions(H,n))do p
    p=>n*(-1)^(length(p)-1)*factorial(length(p)-1)//
    prod(i->factorial(i[2]),tally(p.l))
   end),H)*((b==:e && iseven(n)) ? -1 : 1)
end
function Algebras.basis(H::SymFuncAlgebra,::Val{:h},h::SymFunc{:p})
  sum(((p,c),)->prod(n->Pto(H,n,:h),p.l;init=hbasis(H)())*c,h.d,init=zero(:h,h))
end
function Algebras.basis(H::SymFuncAlgebra,::Val{:e},f::SymFunc{:p})
  sum(((p,c),)->prod(n->Pto(H,n,:e),p.l;init=ebasis(H)())*c,f.d,init=zero(:e,f))
end

Base.getindex(a::SymFunc,p::Partition)=a.d[p]
Base.getindex(a::SymFunc,x::Vararg{Int})=a.d[Partition(collect(x))]

function Base.:*(a::SymFunc{ba},b::SymFunc)where ba
  b=basis(Val(ba),b)
  if ba in (:p,:h,:e)
    clone(a,ModuleElt([union(p,q)=>c*c1 for (p,c) in a.d for (q,c1) in b.d
                       if size(p)+size(q)<=a.A.n]))
  else basis(Val(ba),pbasis(a)*pbasis(b))
  end
end

"`a⊗b` the inner product of the symmetric functions `a` and `b`"
function FinitePosets.:⊗(a::SymFunc{ba},b::SymFunc)where ba
  H=a.A
  a1=pbasis(a)
  b1=pbasis(b)
  m=ModuleElts.merge2((x,y)->x*y,a1.d,b1.d)
  m=ModuleElt(p=>c*centralizer(H,p) for (p,c) in m)
  basis(Val(ba),SymFunc{:p}(m,H))
end

function PuiseuxPolynomials.Mvp(a::SymFunc,
   vars=Mvp.(Symbol("x",stringind(rio(),j)) for j in 1:size(last(a.d.d)[1])))
  sum(pbasis(a).d)do (l,c)
    prod(i->sum(vars.^i),l.l)*c
  end
end

"`scalar_product(a,b)` the scalar product of the symmetric functions `a` and `b`"
function Chars.scalar_product(a::SymFunc{ab},b::SymFunc)where ab
  b=basis(Val(ab),b)
  if ab==:s
    sum(values(ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d)))
  elseif ab==:p
    sum(((p,c),)->centralizer(a.A,p)*c,
        ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d))
  else scalar_product(pbasis(a),pbasis(b))
  end
end

# raises all variables to rth power
powr(m::Monomial,r)=Monomial(m.d*r)
powr(p::Mvp,r)=Mvp(ModuleElt(powr(m,r)=>c for (m,c) in p.d))
powr(p::Frac{<:Mvp},r)=Frac(powr(p.num,r),powr(p.den,r))
powr(p,r)=p

function plethysm(a::SymFunc{ab},b::SymFunc)where ab
  a=pbasis(a)
  b=pbasis(b)
  basis(Val(ab),sum(a.d;init=zero(a))do (p,c)
    c*prod(p.l)do i
     SymFunc{:p}(ModuleElt(Partition(i*p1.l)=>powr(c1,i) for (p1,c1) in b.d 
                           if i*size(p1)<=a.A.n),a.A)
     end
    end)
end

Base.getindex(a::SymFunc,b::SymFunc)=plethysm(a,b)

end
