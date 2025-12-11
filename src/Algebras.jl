"""
This  is a port of the GAP3 package Algebras by Cédric Bonnafé. 

Finite  dimensional algebras are of the abstract type `FiniteDimAlgebra{C}`
where the parameter `C` is the type of the coefficients.

There  are few generic methods for such an  algebra `A `. By default, it is
assumed  there is  a method  `dim` (which  can also be spelled `dimension`)
giving  the  dimension,  and  that  `A`  has  some canonical basis index by
`1:dim(A)`.  Elements of  `A`, of  type `AlgebraElt{TA}`  where `TA` is the
type  of `A`, are internally `ModuleElt`s containing pairs `i=>c` where `i`
is  the integer index in `1:dim(A)` specifiying  an element of the basis of
`A`  and  `c`  is  the  coefficient  of  that  basis  element. The function
`basis(A,i)` returns the `AlgebraElt` with only non-zero coefficient on `i`
and  equal to `1`. It is assumed  that `basis(A,1)` is the identity of `A`.
The function `basis(A)` returns `basis.(A,1:dim(A))`. Addition, subtraction
and  multiplication by scalars (anything of type `<:Number` or of type `C`)
of `AlgebraElt`s works.

If  `A` has  a method  `basismult(A,i,j)` returning `basis(A,i)*basis(A,j)`
then  the multiplication of elements is implemented using this function. If
`A`  has  a  field  `.multable`  containing  a  matrix  `dim(A)×dim(A)`  of
`AlgebraElt`s,  then `basismult` by  default uses this  table understood as
giving the multiplication table of `basis(A)`.

Another  view of algebra elements is as the vector of their coefficients on
each  basis  element.  If  `a`  is an `AlgebraElt` then `v=coefficients(a)`
returns  this  vector  and  `AlgebraElt(A,v)`  recovers  `a`.  The function
`Matrix(a)`  will  transform  `a`  to  a  matrix  whose  `i`-th  row is the
coefficients  of `a` times the `i`-th basis  element (such a matrix acts on
the left on the vector space underlying `A`).

Generic  functions  are  `left_ideal`  and  `twosided_ideal`  which take as
argument `A` and a list of `AlgebraElt` and return the generated ideal.

`radical`   is  not  defined  generically,  but,  if  defined,  it  defines
automatically `radicalpower` and `loewylength`.

Some of the implemented algebras in this module are 
the [`PolynomialQuotientAlgebra`](@ref) of a polynomial, 
the [`GrothendieckRing`](@ref)  and the [`GroupAlgebra`](@ref)  of a group, 
and the [`ZeroHecke`](@ref) and [`SolomonAlgebra`](@ref) of a Coxeter group.
"""
module Algebras
using ..Chevie
export FiniteDimAlgebra, coefftype, AlgebraElt, basis, idempotents,
       involution, isassociative, SubAlgebra, TwoSidedIdeal, 
       loewylength, radicalpower, centralidempotents,
       PolynomialQuotientAlgebra, GrothendieckRing, GroupAlgebra,
       Quaternions, SolomonAlgebra, ZeroHecke
# finite dimensional algebra over ring whose elements have type C
abstract type FiniteDimAlgebra{C} end

"""
`coefftype(A::FiniteDimAlgebra{C}) where C=C`
returns the type of the coefficients of the algebra `A`.
"""
coefftype(A::FiniteDimAlgebra{C}) where C=C
coefftype(T::Type{<:FiniteDimAlgebra{C}}) where C=C

InfoAlgebra=print

"""
`ppart(g,p)`

Given  an element `g` of  a group and an  integer `p` we can uniquely write
`g=g₁g₂` where the order of `g₂` is prime to `p` and the order of `g₁` is a
product of primes which divide `p`. this function returns `g₁`.
```julia-repl
julia> Algebras.ppart(Perm(1,2,3,4,5,6),2)
(1,4)(2,5)(3,6)

julia> Algebras.ppart(Perm(1,2,3,4,5,6),3)
(1,5,3)(2,6,4)
```
"""
function ppart(g,k)
  n=order(g)
  np=prod(p^m for (p,m) in factor(n) if k%p==0;init=1)
  g^(gcdx(np,div(n,np))[3]*div(n,np))
end

"""
`pprimesections(G::Group,p::Integer)`

This  function returns a partition of  the indices of the conjugacy classes
of  `G`. Indices `i,j` are in the  same part if the elements `x/ppart(x,p)`
for `x` in `classreps(G)[i],classreps(G)[j]` are in the same class.
```julia-repl
julia> G=symmetric_group(5)
Group((1,2),(2,3),(3,4),(4,5))

julia> map(x->classreps(G)[x],Algebras.pprimesections(G,2))
3-element Vector{Vector{Perm{Int16}}}:
 [(), (4,5), (2,3)(4,5), (2,3,4,5)]
 [(3,4,5), (1,2)(3,4,5)]
 [(1,2,3,4,5)]

julia> map(x->classreps(G)[x],Algebras.pprimesections(G,3))
5-element Vector{Vector{Perm{Int16}}}:
 [(), (3,4,5)]
 [(4,5), (1,2)(3,4,5)]
 [(2,3)(4,5)]
 [(2,3,4,5)]
 [(1,2,3,4,5)]
```
"""
function pprimesections(G::Group,k::Integer)
  pp=map(x->x/ppart(x,k),classreps(G))
  collectby(i->position_class(G,pp[i]),eachindex(pp))
end

# Algebra element of algebra A over ring elements of type T
# d coefficients on basis indexed by number
# basis[1] should be identity element
struct AlgebraElt{TA<:FiniteDimAlgebra,C}
  A::TA
  d::ModuleElt{Int,C}
end

Base.broadcastable(h::AlgebraElt)=Ref(h)

AlgebraElt(A::FiniteDimAlgebra,v::AbstractVector)=
  AlgebraElt(A,ModuleElt(Pair.(1:dim(A),v)))

Weyl.dim(A::FiniteDimAlgebra)=error("not implemented in general")

idempotents(A::FiniteDimAlgebra)=error("not implemented in general")

struct AlgebraHom{TA<:FiniteDimAlgebra,TB<:FiniteDimAlgebra}
  source::TA
  target::TB
  images::Vector{<:AlgebraElt{TB}}
end

(H::AlgebraHom)(h::AlgebraElt)=sum(H.images[k]*c for (k,c) in h.d)

"""
`gens(A::FiniteDimAlgebra)`

The generators of the algebra `A`. By default `basis(A)` but should
be specialized when a more efficient answer is known.
"""
Groups.gens(A::FiniteDimAlgebra)=basis(A)

"`isabelian(A::FiniteDimAlgebra)` whether `A` is commutative"
function Groups.isabelian(A::FiniteDimAlgebra)
  get!(A,:isabelian)do
    all(x*y==y*x for x in gens(A), y in gens(A))
  end::Bool
end

"`isassociative(A::FiniteDimAlgebra)` whether `A` is associative"
function isassociative(A::FiniteDimAlgebra)
  get!(A, :isassociative)do
    l=basis(A)
    all((x*y)*z==x*(y*z) for x in l, y in l, z in l)
  end::Bool
end

function basis(A::FiniteDimAlgebra{T};force=false)where T
  get!(A,:basis)do
    if !applicable(dim,A) error(A," has no dim method") end
    if !force && dim(A)>10000
      error(A," is of dimension ",dim(A),"\nif you really want to do that ",
            "call basis(A;force=true)")
    end
    map(i->basis(A,i),1:dim(A))
  end
end

basis(A::FiniteDimAlgebra{T},i::Integer) where T =
  AlgebraElt(A,ModuleElt(Int(i)=>T(1)))

Base.one(A::FiniteDimAlgebra)=basis(A,1)

function Base.show(io::IO, h::AlgebraElt)
  if !hasdecor(io)
    print(io,"AlgebraElt(",h.A,",",h.d,")")
    return 
  end
  if haskey(h.A,:showbasis) showbasis=h.A.showbasis
  else showbasis(io::IO,e)=fromTeX(io,"B"*"_{"*string(e)*"}")
  end
  show(IOContext(io,:showbasis=>showbasis),h.d)
end

Base.:+(a::AlgebraElt,b::AlgebraElt)=AlgebraElt(a.A,a.d+b.d)
Base.:+(a::AlgebraElt,b::Number)=a+one(a.A)*b
Base.:-(a::AlgebraElt)=AlgebraElt(a.A,-a.d)
Base.:-(a::AlgebraElt, b::AlgebraElt)=a+(-b)
Base.:*(a::AlgebraElt, b::Vector)=a.*b #for tbl
Base.:*(a::AlgebraElt, b)=AlgebraElt(a.A,a.d*b)
Base.:*(b,a::AlgebraElt)=AlgebraElt(a.A,b*a.d)
Base.://(a::AlgebraElt, b)=AlgebraElt(a.A,a.d//b)
Base.zero(a::AlgebraElt)=AlgebraElt(a.A,zero(a.d))
Base.:^(a::AlgebraElt, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                          Base.power_by_squaring(inv(a),-n)
Base.isless(a::AlgebraElt,b::AlgebraElt)=isless(a.d,b.d)
Base.:(==)(a::AlgebraElt,b::AlgebraElt)=a.d==b.d
Base.hash(a::AlgebraElt,k::UInt)=hash(a.d,k)
Base.one(a::AlgebraElt)=one(a.A)

function LaurentPolynomials.coefficients(a::AlgebraElt)
  v=fill(zero(valtype(a.d)),dim(a.A))
  for (i,c) in a.d v[i]=c end
  v
end

Base.Matrix(a::AlgebraElt)=toM(map(coefficients,a.*basis(a.A)))

using LinearAlgebra
function Base.inv(a::AlgebraElt)
  m=Matrix(a)
  d=LinearAlgebra.det_bareiss(m)
  if iszero(d) error(a," is not invertible") end
  mi=inv(m//1)
  if eltype(m)<:Integer && (isone(d) || d==-1) mi=Integer.(mi) end
  AlgebraElt(a.A,mi[1,:])
end

GenLinearAlgebra.charpoly(a::AlgebraElt)=Pol(charpoly(Matrix(a)))

function minpoly(a::AlgebraElt)
  p=charpoly(a)
  exactdiv(p,gcd(p,derivative(p)//1))
end

# generic fallback
basismult(A::FiniteDimAlgebra,i::Integer,j::Integer)=A.multable[i,j]

function Base.:*(a::AlgebraElt, b::AlgebraElt)
  res=Pair{Int,promote_type(valtype(a.d),valtype(b.d))}[]
  for (i,c) in a.d, (i1,c1) in b.d
    append!(res,[k=>c*c1*c2 for (k,c2) in basismult(a.A,i,i1)])
  end
  AlgebraElt(a.A,ModuleElt(res))
end

function GenLinearAlgebra.rowspace(v::AbstractVector{<:ModuleElt})
  v=filter(!iszero,v)
  if isempty(v) return v end
  kk=unique_sorted!(sort(vcat(map(x->collect(keys(x)),v)...)))
  V=fill(zero(first(values(v[1]))),length(v),length(kk))
  for i in eachindex(v), (k,c) in v[i]
     V[i,searchsortedfirst(kk,k)]=c
  end
  map(r->ModuleElt(Pair.(kk,r)),eachrow(rowspace(V)))
end

function GenLinearAlgebra.rowspace(v::AbstractVector{<:AlgebraElt})
  if isempty(v) return v end
  A=v[1].A
  map(x->AlgebraElt(A,x),rowspace(map(x->x.d,v)))
end

function saturate_left(gens,v)
  l=0
  newl=1
  while newl>l
    l=newl
    v=rowspace(vcat(v,map(g->g.*v,gens)...))
    newl=length(v)
  end
  v
end

function saturate_right(gens,v)
  l=0
  newl=1
  while newl>l
    l=newl
    v=rowspace(vcat(v,map(g->v.*g,gens)...))
    newl=length(v)
  end
  v
end

#----------------------------- LeftIdeal --------------------------------
@GapObj struct LeftIdeal{T<:FiniteDimAlgebra}
  parent::T
  basis
end

function left_ideal(A,v::AbstractVector)
  LeftIdeal(A,saturate_left(gens(A),rowspace(v)),Dict{Symbol,Any}())
end

function Base.show(io::IO,I::LeftIdeal)
 print(io,"LeftIdeal(",I.parent,",[")
 join(io,I.basis,",")
 print(io,"])")
end

basis(I::LeftIdeal)=I.basis
Weyl.dim(I::LeftIdeal)=length(basis(I))
#----------------------------- TwoSidedIdeal ----------------------------
@GapObj struct TwoSidedIdeal{T<:FiniteDimAlgebra}
  parent::T
  basis
end

function twosided_ideal(A,v::AbstractVector)
  v=rowspace(v)
  v=saturate_left(gens(A),v)
  if !isabelian(A) v=saturate_right(gens(A),v) end
  TwoSidedIdeal(A,v,Dict{Symbol,Any}())
end

function Base.show(io::IO,I::TwoSidedIdeal)
 print(io,"TwoSidedIdeal(",I.parent,",[")
 join(io,I.basis,",")
 print(io,"])")
end

Groups.gens(I::TwoSidedIdeal)=basis(I)
basis(I::TwoSidedIdeal)=I.basis
Weyl.dim(I::TwoSidedIdeal)=length(basis(I))
#------------------------------------------------------------------------

# radicalpower(A,n) computes radical(A)^n. This is again a two-sided ideal
function radicalpower(A::FiniteDimAlgebra,n)
  if n==0 return A
  elseif n==1 return radical(A)
  elseif haskey(A,:radicalpowers) 
   if length(A.radicalpowers)>=n return A.radicalpowers[n] end
  else A.radicalpowers=[radical(A)]
  end
  res=vcat(map(b->b.*basis(radicalpower(A,n-1)),basis(radical(A)))...)
  ideal=TwoSidedIdeal(A,rowspace(res),Dict{Symbol,Any}())
# if dim(ideal)==0
#   ideal.operations.LeftTrace=x->0*A.field.one;
#   ideal.operations.RightTrace=x->0*A.field.one;
# else 
#   ideal.operations.LeftTrace=x->TraceMat(List(base.vectors,v->
#     Coefficients(base,Matrix(x*AlgebraElt(A,v)))));
#   ideal.operations.RightTrace=x->TraceMat(List(base.vectors,v->
#     Coefficients(base,Matrix(A.AlgebraElt(A,v)*x))));
# end
  push!(A.radicalpowers,ideal)
  ideal
end

# loewylength(A) is the smallest d such that radicalpower(A,d)=0
function loewylength(A::FiniteDimAlgebra)
  get!(A,:loewylength)do
    d=1
    if length(gens(radical(A)))==1 
      gen=gens(radical(A))[1] 
      res=gen
      while !iszero(res) d+=1; res*=gen end
    else while dim(radicalpower(A,d))>0 d+=1 end
    end
    d
  end
end 

function LeftIndecomposableProjectives(A::FiniteDimAlgebra)
  get!(A,:leftprojectives)do
    map(i->left_ideal(A,[i]),idempotents(A))
  end
end

function centralidempotents(A::FiniteDimAlgebra)
  get!(A,:centralidempotents) do
    idem=idempotents(A)
    mat=[length(rowspace(i*basis(j))) for i in idem, 
                           j in LeftIndecomposableProjectives(A)]
    A.blocks=diagblocks(mat)
    perm=vcat(A.blocks...)
    A.cartan=mat[perm,perm]
    map(i->sum(idem[i]),A.blocks)
  end
end

function PermRoot.cartan(A::FiniteDimAlgebra)
  if !haskey(A,:cartan)
    centralidempotents(A) # computes A.blocks and A.cartan
  end
  if applicable(CharTable,A) && hasproperty(CharTable(A),:rows)
    rows=CharTable(A).rows 
  else rows=string.(1:length(idempotents(A)))
  end
  rows=joindigits.(rows[vcat(A.blocks...)])
  m=map(e->iszero(e) ? "." : repr(e),A.cartan)
  Format.Table(m,row_labels=rows)
end

#------------------------- SubAlgebra -------------------------
@GapObj struct SubAlgebra{T}<:FiniteDimAlgebra{T}
  parent::FiniteDimAlgebra{T}
  inclusion::Vector
  space::Matrix{T}
end

Weyl.dim(A::SubAlgebra)=length(A.inclusion)

function SubAlgebra(A::FiniteDimAlgebra,l)
  inclusion=saturate_left(vcat([one(A)],l),vcat([one(A)],l))
  space=improve_type(toM(coefficients.(inclusion)))
  res=SubAlgebra(A,inclusion,space,Dict{Symbol,Any}())
  res.gens=restriction.(Ref(res),l)
  res
end

Groups.gens(A::SubAlgebra)=A.gens

function Base.show(io::IO,A::SubAlgebra)
  print(io,"SubAlgebra(",A.parent,",",inclusion.(Ref(A),gens(A)),")")
end
  
#  SUB.injection:=AlgebraHomomorphismByLinearity(SUB,A,SUB.basisinclusion);
function PermRoot.inclusion(A::SubAlgebra,h::AlgebraElt)
  sum(c*A.inclusion[i] for (i,c) in h.d)
end

#  SUB.belongs:=v->A.EltToVector(v) in subspace;
# restriction of basis element b to A
function PermRoot.restriction(A::SubAlgebra,b)
  res=solutionmat(A.space,coefficients(b))
  if !isnothing(res) AlgebraElt(A,improve_type(res)) end
end

function basismult(A::SubAlgebra,i::Integer,j::Integer)
  restriction(A,A.inclusion[i]*A.inclusion[j]).d.d
end
#----------------------- Quaternions -----------------------------------
@GapObj struct Quaternions{T}<:FiniteDimAlgebra{T}
  multable::Matrix{Vector{Pair{Int,T}}}
end

"""
`Quaternions()` The quaternion algebra.
```julia-repl
julia> b=basis(Quaternions())
4-element Vector{AlgebraElt{Quaternions{Int64}, Int64}}:
 1
 i
 j
 k

julia> b*permutedims(b)
4×4 Matrix{AlgebraElt{Quaternions{Int64}, Int64}}:
 1  i   j   k
 i  -1  k   -j
 j  -k  -1  i
 k  j   -i  -1
```
"""
function Quaternions(a::T=-1,b::T=-1)where T
  multable=[[1=>1] [2=>1]  [3=>1] [4=>1];
            [2=>1] [1=>a]  [4=>1] [3=>-1];
            [3=>1] [4=>-1] [1=>b] [2=>-b];
            [4=>1] [3=>-a] [2=>b] [1=>-a*b]]
  multable=map(x->[x],multable)
  A=Quaternions(multable,Dict{Symbol,Any}())
  A.a=a;A.b=b
  A.isabelian=false
  A.showbasis=(io::IO,i)->["","i","j","k"][i]
  A
end

Weyl.dim(A::Quaternions)=4

Base.show(io::IO,A::Quaternions{T}) where T=print(io,"Quaternions(",A.a,",",A.b,")")

Groups.gens(A::Quaternions)=basis.(Ref(A),[2,3])

involution(h::AlgebraElt{<:Quaternions})=
  AlgebraElt(h.A,coefficients(h).*[1,-1,-1,-1])

LinearAlgebra.norm(h::AlgebraElt{<:Quaternions})=h*involution(h)
#----------------------- 0-Hecke -----------------------------------
#@GapObj struct ZeroHecke{TE,TW<:CoxeterGroup{TE}}<:FiniteDimAlgebra{Int}
@GapObj struct ZeroHecke{TE,TW<:CoxeterGroup{TE},T}<:FiniteDimAlgebra{T}
  W::TW
  e::Vector{TE}
  u::T
end

Base.show(io::IO,A::ZeroHecke)=
   print(io,hasdecor(io) ? "0-Hecke" : "ZeroHecke","(",A.W,")")

Weyl.dim(A::ZeroHecke)=length(A.W)

function basismult(A::ZeroHecke,i::Integer,j::Integer)
  res=A.e[j]
  for k in reverse(word(A.W,A.e[i]))
    if !isleftdescent(A.W,res,k) res=A.W(k)*res end
  end
  [findfirst(==(res),A.e)=>1]
end
  
function Tbasis(A::ZeroHecke)
  return function(v::Int...)
    prod(i->gens(A)[i],v;init=one(A))
  end
end

Groups.gens(A::ZeroHecke)=basis.(Ref(A),2:ngens(A.W)+1)

"""
`ZeroHecke(W::CoxeterGroup,T=Int)` create the Hecke algebra of `W` with
parameters 0 and 1.
```julia-repl
julia> A=ZeroHecke(coxsym(3))
0-Hecke(𝔖 ₃)

julia> radical(A)
TwoSidedIdeal(0-Hecke(𝔖 ₃),[(1//1)T₁₂+(-1//1)T₁₂₁,(1//1)T₂₁+(-1//1)T₁₂₁])
```
"""
function ZeroHecke(W::CoxeterGroup,T=Int)
  e=elements(W);# we use they are stored by increasing CoxeterLength
  A=ZeroHecke(W,e,one(T),Dict{Symbol,Any}())
  A.isabelian=false
  A.showbasis=(io::IO,i)->fromTeX(io,string("T_{",joindigits(word(A.W,A.e[i])),"}"))
  A
end

function PermRoot.radical(A::ZeroHecke)
  get!(A,:radical)do
    c=groupby(x->unique!(sort!(word(A.W,A.e[x]))),eachindex(A.e))
    l=filter(x->length(x)>1,collect(values(c)))
    twosided_ideal(A,vcat(map(i->map(j->basis(A,i[1])-basis(A,j),i[2:end]),l)...))
  end
end

#    subsets:=Combinations([1..W.semisimpleRank]));
#  A.tbasis:=function(arg) local mot;
#    if Length(arg)>1 then mot:=arg;
#    elif IsList(arg[1]) then mot:=arg[1];
#    elif IsInt(arg[1]) then mot:=Digits(arg[1]);
#    else mot:=CoxeterWord(W,arg[1]);
#    fi;
#    if mot=[] then return A.basis[findfirst(isempty,A.mots)];
#    else return prod(i->A.basis[findfirst(==([i]),A.mots)],mot);
#    fi;
#  end;
#  A.radical:=TwoSidedIdeal(A,vcat(A.radical...));
#  A.radicalpowers:=[A.radical];
#  A.embedding:=function(g) return A.basis[findfirst(==(g),e)];end;

#------------------------  T[q]/p(q)----------------------------------------
## The function A.class sends a polynomial to its image in A
## An element of A is printed as "Class(polynomial)"
@GapObj struct PolynomialQuotientAlgebra{T}<:FiniteDimAlgebra{T}
  p::Pol{T}
  var::Symbol
  multable::Matrix{Vector{Pair{Int,T}}}
end

Weyl.dim(A::PolynomialQuotientAlgebra)=degree(A.p)

Groups.gens(A::PolynomialQuotientAlgebra)=[basis(A,2)]

Base.show(io::IO,A::PolynomialQuotientAlgebra)=
    print(io,eltype(A.p.c),"[",A.var,"]/",A.p)

"""
`PolynomialQuotientAlgebra(p::Pol)`

The quotient of the polynomial algebra with coefficients of type `eltype(p)`
by the polynomial `p`.

```julia-repl
julia> A=PolynomialQuotientAlgebra(Pol()^2+1)
Int64[x]/x²+1

julia> basis(A)
2-element Vector{AlgebraElt{PolynomialQuotientAlgebra{Int64}, Int64}}:
 ⋅1
 x
```
"""
function PolynomialQuotientAlgebra(p::Pol{T})where T
#  A.string:=r->string("Class(",sum(i->i[1]*q^(i[2]-1),r.coefficients),")");
  d=degree(p)
  function class(r::Pol)
    r*=one(p)
    r=isone(p[end]) ? pseudodiv(r,p)[2] : rem(r//1,p)
    filter(x->x[2]!=0,map(Pair,1:degree(r)+1,r[0:end]))
  end
  basismult(i,j)=class(Pol([1],i+j-2))
  A=PolynomialQuotientAlgebra(p,LaurentPolynomials.varname,
           [basismult(i,j) for i in 1:d, j in 1:d],Dict{Symbol,Any}())
  A.showbasis=(io::IO,i)->fromTeX(io,i==1 ? "⋅1" :
                           string(A.var)*(i==2 ? "" : string("^{",i-1,"}")))
  A.isAbelian=true
  A.class=class
  A.fp=sort(collect(factor(p)),by=x->degree(x[1]))
  if T<:FFE A.char=char(T)
  else A.char=0
  end
  A
end

function idempotents(A::PolynomialQuotientAlgebra)
  map(A.fp)do (pp,mul)
    a=pp^mul
    b=exactdiv(A.p,a)
    bezout=gcdx(a*1//1,b*1//1)
    res=bezout[2]*a+bezout[3]*b
    res=res[0]
    res=(bezout[3]*b)//res
    AlgebraElt(A,ModuleElt(A.class(res)))
  end
end

PermRoot.cartan(A::PolynomialQuotientAlgebra)=Diagonal(last.(A.fp))

function Chars.CharTable(A::PolynomialQuotientAlgebra)
  factors=first.(A.fp)
  irr=fill(0,dim(A),length(factors))
  for k in eachindex(factors)
    coefs=factors[k][0:end]
    deg=degree(factors[k])
    mat=fill(zero(coefs[1]),deg,deg)
    for i in 1:deg-1  mat[i+1,i]=1 end
    for i in 1:deg  mat[i,deg]=-coefs[i] end
    for j in 1:dim(A) irr[j,k]=tr(mat^(j-1)) end
  end
  Format.Table(irr,col_labels=factors,row_labels=basis(A))
end

function PermRoot.radical(A::PolynomialQuotientAlgebra)
  get!(A,:radical) do
    if A.char==0
      der=derivative(A.p)
      res=gcd(A.p//gcd(der,A.p),der)
      A.radical=twosided_ideal(A,[AlgebraElt(A,ModuleElt(A.class(res)))])
      A.radicalpowers=[A.radical]
      return A.radical
    end
    res=prod(first.(filter(i->i[2]>1,A.fp)),init=Pol([1],0))
    A.radical=twosided_ideal(A,[AlgebraElt(A,ModuleElt(A.class(res)))])
    A.radicalpowers=[A.radical]
    A.radical
  end
end

@GapObj struct GrothendieckRing{T}<:FiniteDimAlgebra{T}
  G::Group
  ct::CharTable
  multable::Matrix{Vector{Pair{Int,T}}}
  positionid::Int
end

Weyl.dim(A::GrothendieckRing)=size(A.ct.irr,1)

function Chars.CharTable(A::GrothendieckRing{T})where T
  if char(T)==0 rows=1:dim(A)
    irr=permutedims(A.ct.irr)
  else cl=classreps(A.G)
    rows=filter(i->gcd(order(cl[i]),char(T))==1,eachindex(cl))
    irr=T.(permutedims(A.ct.irr)[rows,:])
  end
  CharTable(irr,A.ct.classnames[rows],A.ct.charnames,A.ct.centralizers[rows],
            A.ct.order,Dict{Symbol,Any}(:name=>string(A)))
end

Base.show(io::IO,A::GrothendieckRing{T}) where T=
  print(io,"GrothendieckRing(",A.G,",",T,")")

function PermRoot.radical(A::GrothendieckRing{T})where T
  get!(A,:radical)do
    if char(T)==0 || length(A.G)%char(T)!=0
      b=empty(basis(A))
    else
      p=char(T); while p<dim(A) p*=char(T) end
      mat=toM(coefficients.(basis(A).^p))
      mat=rowspace(lnullspace(mat))
      b=map(x->AlgebraElt(A,x),eachrow(mat))
    end
    twosided_ideal(A,b)
  end
end

## La fonction Degre envoie un element de l'anneau de 
## Grothendieck sur son degre (virtuel a  priori bien sur)
function LaurentPolynomials.degree(h::AlgebraElt)
  sum(c*Int(h.A.ct.irr[k,1]) for (k,c) in h.d)
end

function PermRoot.cartan(A::GrothendieckRing{T})where T
  if char(T)==0 Diagonal(fill(1,dim(A)))
  else Diagonal(length.(pprimesections(A.G,char(T))))
  end
end

Base.one(A::GrothendieckRing)=basis(A,A.positionid)

"""
`GrothendieckRing(G::Group,T=Int)`

constructs the Grothendieck ring of the group `G` with coefficients of type
`T`.  The only condition  is to be  able to compute  the character table of
`G`.
```julia-repl
julia> A=GrothendieckRing(coxsym(3))
GrothendieckRing(𝔖 ₃,Int64)

julia> basis(A).*permutedims(basis(A))
3×3 Matrix{AlgebraElt{GrothendieckRing{Int64}, Int64}}:
 [3]    [21]            [111]
 [21]   [111]+[21]+[3]  [21]
 [111]  [21]            [3]

julia> A=GrothendieckRing(coxsym(3),FFE{2})
GrothendieckRing(𝔖 ₃,FFE{2})

julia> radical(A)
TwoSidedIdeal(GrothendieckRing(𝔖 ₃,FFE{2}),[[111]+[3]])
```
"""
function GrothendieckRing(G::Group,T=Int)
  ct=CharTable(G)
  d=size(ct.irr,1)
  basismult(i,j)=filter(x->!iszero(x[2]),map(Pair,1:size(ct.irr,1),
         T.(decompose(ct,ct.irr[i,:].*ct.irr[j,:]))))
  multable=[basismult(i,j) for i in 1:d, j in 1:d]
  positionid=findfirst(i->all(isone,ct.irr[i,:]),1:d)
  A=GrothendieckRing{T}(G,ct,multable,positionid,Dict{Symbol,Any}())
  A.isabelian=true
  A.showbasis=(io::IO,i)->"["*fromTeX(io,ct.charnames[i])*"]"
  A
end;

function idempotents(A::GrothendieckRing{T}) where T
  e=map(i->sum((conj.(A.ct.irr[:,i]).//A.ct.centralizers[i]).*basis(A)),1:dim(A))
  if char(T)==0 return e end
  e=map(i->sum(e[i]),pprimesections(A.G,char(T)))
  map(x->AlgebraElt(A,ModuleElt(k=>T(c) for (k,c) in x.d)),e)
end

#GrothendieckRingOps.PrincipalIdempotent:=function(A)local G,e,c,res;
#  G:=A.group;
#  e:=Idempotents(GrothendieckRing(G));
#  c:=position_class(G,G.identity);
#  res:=sum(e{First(pprimesections(G,A.field.char),i->c in i)});
#  return AlgebraElement(A,List(res.coefficients,i-> 
#    [A.field.one*Numerator(i[1])/Denominator(i[1]),i[2]]));
#end;

#GrothendieckRingOps.PrincipalLoewyLength:=function(A)local e;
#  e:=GrothendieckRingOps.PrincipalIdempotent(A);
#  LoewyLength(A);
#  e:=List(A.radicalpowers,i->LeftIdeal(A,i.basis*e));
#  A.principalloewyseries:=Filtered(e,i->i.dimension>0);
#  return Length(A.principalloewyseries)+1;
#end;

#RestrictionHomomorphism:=function(A,B)local k;
#  k:=induction_table(B.group,A.group).scalar;
#  return AlgebraHomomorphismByLinearity(A,B,k*B.basis);
#end;

function AdamsOperation(A::GrothendieckRing,n)
   r=map(x->position_class(A.G,x^n),classreps(A.G))
   v=map(x->decompose(A.ct,x[r]),eachrow(A.ct.irr))
   AlgebraHom(A,A,map(x->AlgebraElt(A,x),v))
end

#----------------- Group Algebras -----------------------------------------
@GapObj struct GroupAlgebra{T,E,TG<:Group{E}}<:FiniteDimAlgebra{T}
  group::TG
  elts::Vector{E}
  one::Pair{Int,T}
end

#GroupAlgebraOps.CharacterDegrees:=function(A) 
#  if A.field.char=0 then 
#    A.characterDegrees:=CharacterDegrees(A.group);
#  else Error("Cannot compute CharacterDegrees in positive characteristic");
#  fi;
#  return A.characterdegrees;
#end;

function Chars.CharTable(A::GroupAlgebra)
  if coefftype(A)<:FFE
    error("Cannot compute CharTable in positive characteristic")
  end
  A.charTable=CharTable(A.group)
#  conj:=List(Elements(A.group),i->position_class(A.group,i));
#  A.charTable.basistraces:=List(A.charTable.irr,chi->chi{conj});
  return A.charTable
end

function PermRoot.radical(A::GroupAlgebra{T})where T
  if char(T)!=0 error("not implemented in positive characteristic") end
  TwoSidedIdeal(A,empty(basis(A)),Dict{Symbol,Any}())
end

Base.show(io::IO,A::GroupAlgebra{T}) where T =print(io,"Algebra(",A.group,",",T,")")

embedding(A::GroupAlgebra{T,E},g::E) where{T,E}=findfirst(==(g),A.elts)

"""
`GroupAlgebra(G,T=Int)` group algebra of `G` with coefficients `T`
```julia-repl
julia> W=Group(Perm(1,2))
Group((1,2))

julia> A=GroupAlgebra(W)
Algebra(Group((1,2)),Int64)

julia> basis(A)
2-element Vector{AlgebraElt{GroupAlgebra{Int64, Perm{Int16}, PermGroups.PG{Int16}}, Int64}}:
 e
 e₁
```
"""
function GroupAlgebra(G,T=Int)
  A=GroupAlgebra(G,elements(G),1=>T(1),Dict{Symbol,Any}())
  A.gens=map(g->basis(A,embedding(A,g)),gens(G))
  if A.elts[1]!=one(G) error(one(G)," should be first in ",elements(G)) end
  A.showbasis=function(io::IO,e)
    fromTeX(io,"e"*"_{"*joindigits(word(A.group,A.elts[e]))*"}")
  end
  if T<:FFE
    A.char=char(T)
  else
    A.char=0
    A.radical=TwoSidedIdeal(A,AlgebraElt[],Dict{Symbol,Any}())
    A.radicalpowers=[A.radical];
  end
  A
end

Weyl.dim(A::GroupAlgebra)=length(A.elts)

Groups.gens(A::GroupAlgebra)=A.gens

basismult(A::GroupAlgebra{T},i::Integer,j::Integer) where T=
  [embedding(A,A.elts[i]*A.elts[j])=>T(1)]

function centralidempotents(A::GroupAlgebra{T}) where T
  p=map(i->position_class(A.group,i),A.elts)
  res=map(eachrow(CharTable(A.group).irr))do chi
    ModuleElt(map(Pair,1:dim(A),chi[1].*conj.(chi[p]).//dim(A)))
  end
  if char(T)==0 return map(x->AlgebraElt(A,x),res) end
  pid=map(b->sum(res[b]),blocks(A.group,char(T)))
  map(i->AlgebraElt(A,ModuleElt(convert(ModuleElt{Int,T},i).d)),pid)
end

# augmentation morphism
augmentation(r)=sum(values(r.d))

################################################################################
## La fonction SolomonAlgebra(W,Q) construit l'algebre de Solomon 
## d'un groupe de Coxeter fini W sur le corps Q (s'il est omis, 
## la fonction prend pour Q le corps des rationnels). Le resultat 
## est une algebre A munie en plus des champs suivants (ici, r est 
## le rang semi-simple de W) :
##   A.group  = W
##   A.xbase  = fonction qui envoie la partie I de [1..r] sur x_I
##   A.ybase  = fonction qui envoie la partie I de [1..r] sur y_I
##   A.radical = radical de A : 
##        A.radical.base = base de A.radical (x_I-x_J)_{I et J conjugues}
##        A.radical.dimension = dimension de A.radical
##   A.radicalpowers = [A.radical]
##   A.injection = injection dans l'algebre de groupe (disponible 
##                 si |W| <= 2000 : jusqu'a D5 donc...)
##
## De plus
##   A.base       = base x_I
##   A.parameters = parties de [1..r] (comme nombres : ex. 123 pour [1,2,3])
##   A.basisname   = "X" (par defaut)
##   
## La fonction SolomonAlgebra munit W du record W.solomon contenant 
## les champs suivants :
##   W.solomon.mackey = constantes de structure (coeff de x_I * x_J sur x_K)
##   W.solomon.conjugacy = classes de conjugaison de parties de [1..r] 
##                         reperees par leur position dans A.parameters
##   W.solomon.domain = A
##   W.solomon.subsets = parties de [1..r]

#SolomonAlgebraOps.CharacterDegrees:=function(A) 
#  if A.field.char=0 then 
#    if A.type="Solomon algebra" then 
#      A.characterDegrees:=List(A.group.solomon.conjugacy, i-> 1);
#    elif A.type="Generalized Solomon algebra" then 
#      A.characterDegrees:=List(A.group.generalizedsolomon.conjugacy, i-> 1);
#    else
#      Error("<A> must be a Solomon algebra");
#    fi;
#  else Error("Cannot compute character degrees of Solomon algebra in positive characteristic");
#  fi;
#  return A.characterDegrees;
#end;

@GapObj struct SolomonAlgebra{T}<:FiniteDimAlgebra{T}
  W::FiniteCoxeterGroup
  multable::Matrix{Vector{Pair{Int,T}}}
end

function Base.show(io::IO,A::SolomonAlgebra{T})where T
  print(io,"SolomonAlgebra(",A.W,",",T,")")
end

@GapObj struct SolomonCharTable{TA}
  A::TA
  irr::Matrix{Int}
  rows::Vector{Vector{Int}}
end

function Base.show(io::IO,ct::SolomonCharTable)
  print(io,"CharTable(",ct.A,")")
  if !hasdecor(io) return end
  println(io)
  showtable(io,ct.irr,row_labels=joindigits.(ct.rows),dotzero=true)
end

function Chars.CharTable(A::SolomonAlgebra)
  W=A.W
  conj=first.(W.solomon_conjugacy)
  orb=map(i->W.solomon_subsets[i],conj)
  irr=[W.solomon_mackey[j,i][i] for i in conj, j in conj]
  SolomonCharTable(A,irr,orb,Dict{Symbol,Any}())
#  res.columns[Length(res.columns)]:=["0"];
#  res.basistraces:=List([1..Length(res.irreducibles)],function(ii)local r;
#    r:=vcat(List([1..Length(W.solomon.conjugacy)],
#      i->List(W.solomon.conjugacy[i], j-> [j,res.irreducibles[ii][i]])));
#    Sort(r); return List(r,i->i[2]);end);
#  if A.field.char=0 then return res;fi;
#  irr:=CyclotomicModP(res.irreducibles,A.field.char);
#  inc:=collectby(irr,eachindex(irr));Sort(inc);
#  res.rows:=List(inc, i-> res.rows[i[1]]);
#  res.irreducibles:=List(inc, i-> irr[i[1]]);
#  res.basistraces:=CyclotomicModP(res.basistraces,A.field.char);
#  res.parameters:=List(inc, i-> A.parameters[W.solomon.conjugacy[i[1]][1]]);
#  res.invertiblematrix:=List(inc,i-> List(inc,j->irr[i[1]][j[1]]));
#  res.incidence:=inc;
#  return res;
end

function injection(A::SolomonAlgebra{T})where T
  W=A.W
  B=GroupAlgebra(W,T)
# if A.type="Solomon algebra" then 
  X=map(W.solomon_subsets)do i
   sum(j->basis(B,embedding(B,inv(j))),vcat(reduced(reflection_subgroup(W,i),W)...))
  end
# else
#   X:=List(W.generalizedsolomon.subgroups, i-> 
#     sum(i-> B.embedding(i^-1),ReducedRightCosetRepresentatives(W, i)));
# fi;
  AlgebraHom(A,B,X)
end

Weyl.dim(A::SolomonAlgebra)=2^ngens(A.W)

const LIM=10000
"""
`SolomonAlgebra(W,K)`

Let  `(W,S)`  be  a  finite  Coxeter  group.  For  `w∈W`,  let `R(w)={s∈S |
l(ws)>l(w)}`.  If `I` is a subset  of `S`, we set ``Y_I=\\{w∈W|R(w)=I\\}``,
``X_I=\\{w∈W|R(w)⊃ I\\}``.

Note  that ``X_I`` is the set  of minimal length left coset representatives
of ``W/W_I``. 

Let  ``y_I=∑_{w∈Y_I}w``, ``x_I=∑_{w∈X_I}w``. They are elements of the group
algebra `ℤW` of `W` over `ℤ`.

Now,  let  ``Σ(W)  =  ⊕_{I  ⊂  S}  ℤ  y_I  =  ⊕_{I  ⊂ S} ℤ x_I``. This is a
sub-`ℤ`-module of `ℤW`. In fact, Solomon proved that it is a sub-algebra of
`ℤW`  called  the  *Solomon  descent  algebra*.  Now,  let  `K(W)`  be  the
Grothendieck  ring of  `W` and  let `θ:Σ(W)→  K(W)` be  the map  defined by
``θ(x_I)  = Ind_{W_I}^W 1``. Solomon proved that this is an homomorphism of
algebras. We call it the *Solomon homomorphism*.

`SolomonAlgebra`  returns the Solomon descent algebra of the finite Coxeter
group `(W,S)` over `K`. If `S=[s₁,…,sᵣ]`, the element ``x_I`` corresponding
to  the  subset  `I=[s₁,s₂,s₄]`  of  `S`  is  printed  as `X₁₂₄`. Note that
`A:=SolomonAlgebra(W,K)` is endowed with the following fields:

`A.W`: the group `W`

`A.basis`: the basis `(x_I)_{I ⊂ S}`.

`A.xbasis`: the  function   sending   `I`   to   ``x_I``  (for  instance
                     `A.xbasis(1,2,4)=X₁₂₄`)

`A.ybasis`: the function sending the subset `I` to `y_I`.

`A.injection`:  the injection of `A` in the group algebra of `W`, obtained
by calling `SolomonAlgebraOps.injection(A)`.

Note that `SolomonAlgebra(W,K)` endows `W` with the field `W.solomon` which
is a record containing the following fields:

`W.solomon_subsets`: the set of subsets of `S`

`W.solomon_conjugacy`:  conjugacy classes  of parabolic  subgroups of `W` (a
conjugacy   class  is  represented  by  the   list  of  the  positions,  in
`W.solomon.subsets`, of the subsets `I` of `S` such that `W_I` lies in this
conjugacy class).

`W.solomon_mackey`:  essentially  the  structure  constants  of  the Solomon
algebra over the rationals.

```julia-repl
julia> W=coxgroup(:B,4)
B₄

julia> A=SolomonAlgebra(W)
SolomonAlgebra(B₄,Int64)

julia> X=A.xbasis; X(1,2,3)*X(2,4)
2X₂+2X₄

julia> W.solomon_subsets
16-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [1, 2, 3]
 [1, 2, 4]
 [1, 3, 4]
 [2, 3, 4]
 [1, 2]
 [1, 3]
 [1, 4]
 [2, 3]
 [2, 4]
 [3, 4]
 [1]
 [2]
 [3]
 [4]
 []

julia> W.solomon_conjugacy
12-element Vector{Vector{Int64}}:
 [1]
 [2]
 [3]
 [4]
 [5]
 [6]
 [7, 8]
 [9, 11]
 [10]
 [12]
 [13, 14, 15]
 [16]

julia> Algebras.injection(A)(X(1,2,3))
e+e₄+e₃₄+e₂₃₄+e₁₂₃₄+e₂₁₂₃₄+e₃₂₁₂₃₄+e₄₃₂₁₂₃₄
```
"""
function SolomonAlgebra(W::FiniteCoxeterGroup,T=Int)
  r=ngens(W)
  d=2^r
  n=length(W)
#  if r <=7 then A.generators:=List([2..r+1],i->A.basis[i]);fi;
  I=sort(combinations(1:r),by=x->-length(x))
  if !haskey(W,:solomon_mackey)
    if n>LIM InfoAlgebra("# structure constants...") end;c=n
    Idict=Dict(p=>i for (i,p) in enumerate(I))
    subsets=map(x->map(j->Idict[j],combinations(x)),I)
    mackey=[Pair{Int,Int}[] for i in 1:d,j in 1:d]
    for w in elements(W)
      leftascents=Idict[setdiff(1:r,leftdescents(W,w))]
      rightascents=Idict[setdiff(1:r,leftdescents(W,inv(w)))]
 #    l1=map(i->action(W,i,w),1:r)
      l=map(i->(a=findfirst(==(W(i)^w),gens(W));
                (isnothing(a) || isleftdescent(W,w,i)) ? 0 : a),1:r)
 #    if filter(<=(r),l1)!=filter(!iszero,l) error(l1,l2) end
      for i in subsets[leftascents] 
        II=sort!(l[I[i]])
 	for j in subsets[rightascents]
          push!(mackey[i,j],Idict[intersect_sorted(II,I[j])]=>1)
        end
      end
      if c%LIM==0 InfoAlgebra(c,".") end; c-=1
    end
    if n>LIM InfoAlgebra("\n") end
    W.solomon_mackey=ModuleElt.(mackey)
    W.solomon_subsets=I
  end
  if !haskey(W,:solomon_conjugacy)
    if n>LIM InfoAlgebra("# orbits of parabolic subgroups...") end
    orb=1:d
    conj=Vector{Int}[]
    while !isempty(orb)
      i=orb[1]
      j=filter(k->!iszero(W.solomon_mackey[i,k][k]) && 
              length(I[i])==length(I[k]),orb)
      push!(conj,j)
      orb=setdiff(orb,j)
    end
    W.solomon_conjugacy=conj
   if n>LIM InfoAlgebra("\n") end
  end
  basismult(i,j)=[k=>T(c) for (k,c) in W.solomon_mackey[i,j].d]
  multable=[basismult(i,j) for i in 1:d, j in 1:d]
  A=SolomonAlgebra(W,multable,Dict{Symbol,Any}())
  A.xbasis=function(a...)
    a=sort(collect(a))
    basis(A,findfirst(==(a),W.solomon_subsets))
  end
  if iszero(char(T))
    A.xprimebasis=function(a...)
      a=sort(collect(a))
      sum(i->(-1//2)^(length(a)-length(i))*A.xbasis(i...),combinations(a))
    end
  end
  A.ybasis=function(a...)
    a=sort(collect(a))
    res=zero(one(A))
    for i in filter(i->issubset(a,I[i]),eachindex(I))
      res+=(-1)^(length(I[i])-length(a))*basis(A,i)
    end
    res
  end
  A.showbasis=function(io::IO,i)
    fromTeX(io,string("X_{",joindigits(A.W.solomon_subsets[i]),"}"))
  end
  if char(T)==0 
    A.characterDegrees=fill(1,length(W.solomon_conjugacy))
  end
  A
end

function PermRoot.radical(A::SolomonAlgebra{T})where T
  get!(A,:radical) do
  if iszero(char(T))
    if length(A.W)>LIM InfoAlgebra("# Computing radical...") end
      k=vcat(map(i->map(j->basis(A,i[1])-basis(A,j),i[2:end]),
                  A.W.solomon_conjugacy)...)
      radical=TwoSidedIdeal(A,k,Dict{Symbol,Any}())
#   if dim(radical)==0 
#     A.radical.operations.LeftTrace=x->0*A.field.one
#     A.radical.operations.RightTrace=x->0*A.field.one
#   else 
#     subspace=Subspace(A.vectorspace,Matrix(k))
#     base=CanonicalBasis(subspace)
#     A.radical.operations.LeftTrace=x->TraceMat(List(base.vectors,v->
#       Coefficients(base,Matrix(x*AlgebraElt(A,v)))));
#     A.radical.operations.RightTrace=x->TraceMat(List(base.vectors,v->
#       Coefficients(base,Matrix(AlgebraElt(A,v)*x))));
#   end
    A.radicalpowers=[radical]
    if length(A.W)>LIM InfoAlgebra("done!\n") end
    radical
  end
  end
end

function permutation_character(W,S)
  t=induction_table(S,W).scalar
  i=findfirst(x->all(==(1),x),eachrow(CharTable(S).irr))
  permutedims(CharTable(W).irr)*t[:,i]
end

#SolomonHomomorphism:=function(x) local A,i,W,T,irr,mat,WI,p;
#  A:=x.domain; W:=A.group; T:=CharTable(W); irr:=T.irreducibles;
#  mat:=List(irr, i-> 0);
#  for i in x.coefficients do 
#    if A.type="Solomon algebra" then 
#      WI:=reflection_subgroup(W,W.solomon.subsets[i[2]]);
#    elif A.type="Generalized Solomon algebra" then 
#      WI:=W.generalizedsolomon.subgroups[i[2]];
#    fi;
##   p:=PermutationCharacter(W,WI);
#    p:=PermutationCharacterByInductiontable(W,WI); # faster for large groups
#    mat:=mat+i[1]*MatScalarProducts(T,irr,[p])[1];
#  od;
#  return GrothendieckRing(W,A.field).VectorToElt(mat);
#end;

#SolomonAlgebraOps.xprimePrint:=function(r) local res,i,A,f,xprime;
#  A:=r.domain;
#  A.xprimebasisname:="X'";
#  res:=Matrix(r);
#  xprime:=Matrix(List(A.group.solomon.subsets,A.xprimebasis));
#  res:=AlgebraElt(A,res*xprime^-1);
#  res:=res.coefficients;
#  if IsBound(A.xprimeprint) then A.xprimeprint(r);
#  else 
#    f:=function(coef) 
#    if coef[1]=A.field.one then 
#      Print(A.xprimebasisname,"(",A.parameters[coef[2]],")");
#    elif coef[1]=-A.field.one then 
#      Print(" - ",A.xprimebasisname,"(",A.parameters[coef[2]],")");
#    else Print(coef[1],"*",A.xprimebasisname,"(",A.parameters[coef[2]],")");
#    fi;
#    end;
#  if Length(res)=0 then Print("0*",A.xprimebasisname,"(",
#     A.parameters[A.one.coefficients[1][2]],")");
#  else 
#    f(res[1]);
#    for i in [2..Length(res)] do 
#      if res[i][1] < 0 or (res[i][1]=-A.field.one and A.field.char <> 2) 
#         then Print(" - ");f([-res[i][1],res[i][2]]);
#      else Print(" + ");f(res[i]);
#      fi;
#    od;
#  fi;
#  fi;
#end;
#
#ProjectionMatrix:=function(quotient) local res;
#  quotient.zero:=quotient.operations.Zero(quotient);
#  quotient.generators:=quotient.operations.Generators(quotient);
#  Basis(quotient);
#  res:=List(Basis(quotient.factorNum).vectors, 
#    i-> quotient.semiEchelonBasis.operations.Coefficients(Basis(quotient),
#    VectorSpaceProjection(quotient,i)));
#  return res;
#end;
#
#SolomonAlgebraOps.Radical:=function(A)local B,space,quotient,mat;
#  if A.field.char=0 then return A.radical;fi;
#  B:=GrothendieckRing(A.group,GF(A.field.char));
#  space:=B.vectorspace;
#  quotient:=space/Subspace(space,Matrix(Radical(B).basis));
#  mat:=List(A.basis,i->Matrix(SolomonHomomorphism(i)));
#  mat:=mat*ProjectionMatrix(quotient);
#  A.radical:=TwoSidedIdeal(A,AlgebraElt(A,NullspaceMat(mat)));
#  A.radicalpowers:=[A.radical];
#  return A.radical;
#end;

function idempotents(A::SolomonAlgebra{T})where T
  get!(A,:idempotents)do
  W=A.W
  if char(T)==0
    M=CharTable(A).irr
#   if A.type="Solomon algebra" then 
    I=W.solomon_subsets[first.(W.solomon_conjugacy)]
#   elseif A.type="Generalized Solomon algebra" then 
#     I:=List(W.generalizedsolomon.conjugacy, i-> 
#       W.generalizedsolomon.signedcompositions[i[1]]);
#   end
  else 
#    #if A.type="Solomon algebra" then 
    M=CharTable(A).invertiblematrix
    I=CharTable(A).parameters
#    #fi;
  end
  r=length(I)
  e=permutedims(inv(M*1//1))*map(i->A.xbasis(i...),I)
  f=map(i->sum(e[i:r]),1:r)
  n=loewylength(A)
  v=Pol(T[1],1)
  pol=sum(j->binomial(2*n,j)*v^(2*n-j)*(1-v)^j,0:n)
  InfoAlgebra("# Computations to do: ",r)
#  f[1]=one(A)
  for i in 2:r
    f[i]=pol(f[i-1]*f[i]*f[i-1])
    InfoAlgebra(".",r-i)
  end
  InfoAlgebra("\n")
  e[r]=f[r]
  for i in 1:r-1 e[i]=f[i]-f[i+1] end
  e
  end
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Section{Generalized Solomon algebras}
#
#In this subsection, we refer to the paper cite{BH05}.
#
#If  `n` is a non-zero natural number, we  denote by `W_n` the Weyl group of
#type  `B_n` and by `W_{-n}` the Weyl group of type `A_{n-1}` (isomorphic to
#the  symmetric  group  of  degree  `n`).  If `C=[i_1,...,i_r]` is a *signed
#composition* of `n`, we denote by `W_C` the subgroup of `W_n` equal to `W_C
# = W_{i_1} x ... x W_{i_r}`. This is a subgroup generated by reflections (it
#is  not in general a  parabolic subgroup of `W_n`).  Let `X_C = {x ∈ W_C
#|  l(xw) ≥ l(x) ∀ w ∈  W_C}`. Note that `X_C` is the set of
#minimal   length  left   coset  representatives   of  `W_n/W_C`.  Now,  let
#`x_C=∑_{w ∈ X_C} w`. We now define `Σ'(W_n) = ⊕_C ℤ
#`,  where `C` runs over the  signed compositions of `n`. By cite{BH05},
#this  is a subalgebra of  `ZW_n`. Now, let `Y_C`  be the set of elements of
#`X_C`  which are not in  any other `X_D` and  let `y_C=∑_{w ∈ Y_C} w`.
#Then  `Σ'(W_n) =  ⊕_C ℤ  y_C`. Moreover,  the linear map
#`θ' :  Σ^'(W_n)    →   K(W_n)`    defined   by
#`θ'(x_C) = Ind_{W_C}^{W_n} 1` is a *surjective homomorphism* of
#algebras (see cite{BH05}). We still call it the *Solomon homomorphism*.
 
#SolomonAlgebraOps.Print:=function(A) 
#    Print("GeneralizedSolomonAlgebra(",A.group,",",A.field,")");
#end;

end
