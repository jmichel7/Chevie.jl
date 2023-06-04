module Algebras
using ..Gapjm
export FiniteDimAlgebra, AlgebraElt, basis, idempotents,
       involution, isassociative, SubAlgebra, TwoSidedIdeal, 
       loewylength, radicalpower, centralidempotents,
       PolynomialQuotientAlgebra, 
       Quaternions, 
       SolomonAlgebra, 
       ZeroHecke
abstract type FiniteDimAlgebra{T} end

InfoAlgebra=print

"""
`ppart(g,p)`

Given  an element `g` of  a group and an  integer `p` we can uniquely write
`g=g₁g₂` where the order of `g₂` is prime to `p` and the order of `g₁` is a
product of primes which divide `p`. this function returns `g₁`.
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

julia> Algebras.pprimesections(G,2)
3-element Vector{Vector{Int64}}:
 [1, 2, 4, 5]
 [3, 6]
 [7]

julia> Algebras.pprimesections(G,3)
5-element Vector{Vector{Int64}}:
 [1, 3]
 [2, 6]
 [4]
 [5]
 [7]
```
"""
function pprimesections(G::Group,k::Integer)
  pp=map(x->x/ppart(x,k),classreps(G))
  collectby(i->position_class(G,pp[i]),eachindex(pp))
end

# properties of Algebras: have .multable, dim

# Algebra element of algebra A
# d coefficients on basis indexed by number
# basis[1] should be identity element
struct AlgebraElt{TA<:FiniteDimAlgebra,T}
  A::TA
  d::ModuleElt{Int,T}
end

Base.broadcastable(h::AlgebraElt)=Ref(h)

AlgebraElt(A::FiniteDimAlgebra,v::AbstractVector)=
  AlgebraElt(A,ModuleElt(Pair.(1:dim(A),v)))

Weyl.dim(A::FiniteDimAlgebra)=error("not implemented in general")

idempotents(A::FiniteDimAlgebra)=error("not implemented in general")

struct AlgebraHom
  source::FiniteDimAlgebra
  target::FiniteDimAlgebra
  images::Vector
end

(H::AlgebraHom)(h::AlgebraElt)=sum(H.images[k]*c for (k,c) in h.d)

"`isabelian(A::FiniteDimAlgebra)` whether `A` is commutative"
function Groups.isabelian(A::FiniteDimAlgebra)
  get!(A,:isabelian)do
    l=applicable(gens,A) ? gens(A) : basis(A)
    all(x*y==y*x for x in l, y in l)
  end::Bool
end

function isassociative(A::FiniteDimAlgebra)
  get!(A, :isassociative)do
    l=basis(A)
    all((x*y)*z==x*(y*z) for x in l, y in l, z in l)
  end::Bool
end

function basis(A::FiniteDimAlgebra{T})where T
  get!(A,:basis)do
    map(i->AlgebraElt(A,ModuleElt(i=>T(1))),1:dim(A))
  end
end
basis(A::FiniteDimAlgebra,i::Integer)=basis(A)[i]

Base.one(A::FiniteDimAlgebra)=basis(A,1)

function Base.show(io::IO, h::AlgebraElt)
  if !get(io,:limit,false) && !get(io,:TeX,false)
    print(io,"AlgebraElt(",h.A,",",h.d,")")
    return 
  end
  if haskey(h.A,:showbasis) showbasis=h.A.showbasis
  else showbasis=function(io::IO,e)
      fromTeX(io,"B"*"_{"*string(e)*"}")
     end
  end
  show(IOContext(io,:showbasis=>showbasis),h.d)
end

Base.:+(a::AlgebraElt,b::AlgebraElt)=AlgebraElt(a.A,a.d+b.d)
Base.:+(a::AlgebraElt,b::Number)=a+one(a.A)*b
Base.:-(a::AlgebraElt)=AlgebraElt(a.A,-a.d)
Base.:-(a::AlgebraElt, b::AlgebraElt)=a+(-b)
Base.:*(a::AlgebraElt, b::Vector)=a.*b #for tbl
Base.:*(a::AlgebraElt, b)=AlgebraElt(a.A,a.d*b)
Base.:*(b,a::AlgebraElt)=AlgebraElt(a.A,a.d*b)
Base.://(a::AlgebraElt, b)=AlgebraElt(a.A,a.d//b)
Base.zero(a::AlgebraElt)=AlgebraElt(a.A,zero(a.d))
Base.:^(a::AlgebraElt, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                          Base.power_by_squaring(inv(a),-n)
Base.isless(a::AlgebraElt,b::AlgebraElt)=isless(a.d,b.d)
Base.:(==)(a::AlgebraElt,b::AlgebraElt)=a.d==b.d
Base.hash(a::AlgebraElt,k::UInt)=hash(a.d,k)
Base.one(a::AlgebraElt)=one(a.A)

function LaurentPolynomials.coefficients(a::AlgebraElt{A,T})where {A,T}
  v=fill(T(0),dim(a.A))
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

basismult(A::FiniteDimAlgebra,i::Integer,j::Integer)=A.multable[i,j]

function Base.:*(a::AlgebraElt{A,T1}, b::AlgebraElt{A,T2})where {A<:FiniteDimAlgebra,T1,T2}
  res=Pair{Int,promote_type(T1,T2)}[]
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
  v=rowspace(v)
  v=saturate_left(applicable(gens,A) ? gens(A) : basis(A),v)
  LeftIdeal(A,v,Dict{Symbol,Any}())
end

function Base.show(io::IO,I::LeftIdeal)
 print(io,"LeftIdeal(",I.parent,",[")
 join(io,I.basis,",")
 print(io,"])")
end

basis(I::LeftIdeal)=I.basis
Groups.gens(I::LeftIdeal)=basis(I)
Weyl.dim(I::LeftIdeal)=length(basis(I))
#----------------------------- TwoSidedIdeal ----------------------------
@GapObj struct TwoSidedIdeal{T<:FiniteDimAlgebra}
  parent::T
  basis
end

function twosided_ideal(A,v::AbstractVector)
  v=rowspace(v)
  v=saturate_left(applicable(gens,A) ? gens(A) : basis(A),v)
  if !haskey(A,:isabelian) && A.isabelian
    v=saturate_right(applicable(gens,A) ? gens(A) : basis(A),v)
  end
  TwoSidedIdeal(A,v,Dict{Symbol,Any}())
end

function Base.show(io::IO,I::TwoSidedIdeal)
 print(io,"TwoSidedIdeal(",I.parent,",[")
 join(io,I.basis,",")
 print(io,"])")
end

basis(I::TwoSidedIdeal)=I.basis
Groups.gens(I::TwoSidedIdeal)=basis(I)
Weyl.dim(I::TwoSidedIdeal)=length(basis(I))
#------------------------------------------------------------------------

# radicalpower(A,n) computes radical(A)^n. This is again a two-sided ideal
function radicalpower(A::FiniteDimAlgebra,n)
  if n==0 return A
  elseif n==1 return radical(A)
  elseif haskey(A,:radicalpowers) && length(A.radicalpowers)>=n
    return A.radicalpowers[n]
  end
  res=vcat(map(b->b.*basis(radicalpower(A,n-1)),basis(radical(A)))...)
  ideal=TwoSidedIdeal(A,rowspace(res),Dict{Symbol,Any}())
# if dim(ideal)==0
#   ideal.operations.LeftTrace=x->0*A.field.one;
#   ideal.operations.RightTrace=x->0*A.field.one;
# else 
#   ideal.operations.LeftTrace=x->TraceMat(List(base.vectors,v->
#     Coefficients(base,A.EltToVector(x*A.VectorToElt(v)))));
#   ideal.operations.RightTrace=x->TraceMat(List(base.vectors,v->
#     Coefficients(base,A.EltToVector(A.VectorToElt(v)*x))));
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
  showtable(rio(),m,row_labels=rows)
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

function Quaternions(a::T=-1,b::T=-1)where T
  multable=[[1=>1] [2=>1]  [3=>1] [4=>1];
            [2=>1] [1=>a]  [4=>1] [3=>-1];
            [3=>1] [4=>-1] [1=>b] [2=>-b];
            [4=>1] [3=>-a] [2=>b] [1=>-a*b]]
  A=Quaternions(multable,Dict{Symbol,Any}())
  A.a=a;A.b=b
  A.isabelian=false
  A.showbasis=(io::IO,i)->["","i","j","k"][i]
  A
end

Weyl.dim(A::Quaternions)=4

Base.show(io::IO,A::Quaternions{T}) where T=print(io,"Quaternions(",A.a,",",A.b,")")

Groups.gens(A::Quaternions)=[basis(A,2),basis(A,3)]

involution(h::AlgebraElt{<:Quaternions})=
  AlgebraElt(h.A,coefficients(h).*[1,-1,-1,-1])

LinearAlgebra.norm(h::AlgebraElt{<:Quaternions})=h*involution(h)
#----------------------- 0-Hecke -----------------------------------
@GapObj struct ZeroHecke{TE,TW<:CoxeterGroup{TE},T}<:FiniteDimAlgebra{T}
  W::TW
  e::Vector{TE}
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

Groups.gens(A::ZeroHecke)=basis(A)[2:ngens(A.W)+1]

# 0-Hecke algebra of W with parameters 0,1
function ZeroHecke(W::CoxeterGroup)
  e=elements(W);# we use they are stored by increasing CoxeterLength
  A=ZeroHecke(W,e,Dict{Symbol,Any}())
  A.isabelian=false
  A.showbasis=(io::IO,i)->fromTeX(io,string("T_{",joindigits(word(A.W,A.e[i])),"}"))
  A
end

function PermRoot.radical(A::ZeroHecke)
  get!(A,:radical)do
    c=groupby(x->unique!(sort!(word(A.W,A.e[x]))),eachindex(A.e))
    l=filter(x->length(x)>1,collect(values(c)))
    TwoSidedIdeal(A,vcat(map(i->map(j->basis(A,i[1])-basis(A,j),i[2:end]),l)...))
  end
end

#    subsets:=Combinations([1..W.semisimpleRank]));
#  A.tbasis:=function(arg) local mot;
#    if Length(arg)>1 then mot:=arg;
#    elif IsList(arg[1]) then mot:=arg[1];
#    elif IsInt(arg[1]) then mot:=Digits(arg[1]);
#    else mot:=CoxeterWord(W,arg[1]);
#    fi;
#    if mot=[] then return A.basis[Position(A.mots,[])];
#    else return Product(mot,i->A.basis[Position(A.mots,[i])]);
#    fi;
#  end;
#  A.radical:=TwoSidedIdeal(A,Concatenation(A.radical));
#  A.radicalpowers:=[A.radical];
#  A.embedding:=function(g) return A.basis[Position(e,g)];end;

#------------------------  T[q]/p(q)----------------------------------------
## The function A.class sends a polynomial to its image in A
## An element of A is printed as "Class(polynomial)"
@GapObj struct PolynomialQuotientAlgebra{T}<:FiniteDimAlgebra{T}
  p::Pol{T}
  fp::Vector{Pair{Pol{T},Int}}
  var::Symbol
  multable::Matrix{Vector{Pair{Int,T}}}
end

Weyl.dim(A::PolynomialQuotientAlgebra)=degree(A.p)

Groups.gens(A::PolynomialQuotientAlgebra)=basis(A,2)

Base.show(io::IO,A::PolynomialQuotientAlgebra)=
    print(io,eltype(A.p.c),"[",A.var,"]/",A.p)

function PolynomialQuotientAlgebra(p::Pol{T})where T
#  A.string:=r->SPrint("Class(",Sum(r.coefficients, i-> i[1]*q^(i[2]-1)),")");
  d=degree(p)
  function class(r::Pol)
    r*=one(p)
    r=isone(p[end]) ? pseudodiv(r,p)[2] : rem(r//1,p)
    filter(x->x[2]!=0,map(Pair,1:degree(r)+1,r[0:end]))
  end
  basismult(i,j)=class(Pol([1],i+j-2))
  A=PolynomialQuotientAlgebra(p,sort(tally(factor(p)),by=x->degree(x[1])),
                                     LaurentPolynomials.varname[],
           [basismult(i,j) for i in 1:d, j in 1:d],Dict{Symbol,Any}())
  A.showbasis=(io::IO,i)->fromTeX(io,i==1 ? "⋅1" :
                           string(A.var)*(i==2 ? "" : string("^{",i-1,"}")))
  A.isAbelian=true
  A.class=class
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
  showtable(irr,col_labels=map(x->repr(x;context=rio()),factors),
            row_labels=map(i->A.showbasis(rio(),i),1:dim(A)))
end

function PermRoot.radical(A::PolynomialQuotientAlgebra)
  res=prod(first.(filter(i->i[2]>1,A.fp)),init=Pol([1],0))
  A.radical=TwoSidedIdeal(A,[AlgebraElt(A,ModuleElt(A.class(res)))])
  A.radicalpowers=[A.radical]
  A.radical
end

######################################################################
## GrothendieckRing(G,Q) renvoie l'anneau de Grothendieck A=Q Irr(G) 
##(si Q est omis, on prend pour Q le corps des rationnels)
##
##   A.basisname  = "χ" (par defaut)
## 
## La fonction Degre envoie un element de l'anneau de 
## Grothendieck sur son degre (virtuel a  priori bien sur)

@GapObj struct GrothendieckRing{T}<:FiniteDimAlgebra{T}
  G::Group
  ct::CharTable
  multable::Matrix{Vector{Pair{Int,T}}}
  positionid::Int
end

Weyl.dim(A::GrothendieckRing)=size(A.ct.irr,1)

function CharTable(A::GrothendieckRing{T})where T
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
  if char(T)==0 || length(A.G)%char(T)!=0
    b=empty(basis(A))
  else
    p=char(T); while p<dim(A) p*=char(T) end
    mat=toM(coefficients.(basis(A).^p))
    mat=rowspace(lnullspace(mat))
    b=map(x->AlgebraElt(A,x),eachrow(mat))
  end
  radical=TwoSidedIdeal(A,b)
  A.radicalpowers=[radical]
  radical
end

function LaurentPolynomials.degree(h::AlgebraElt)
  sum(c*Int(h.A.ct.irr[k,1]) for (k,c) in h.d)
end

function PermRoot.cartan(A::GrothendieckRing{T})where T
  if char(T)==0 Diagonal(fill(1,dim(A)))
  else Diagonal(length.(pprimesections(A.G,char(T))))
  end
end

Base.one(A::GrothendieckRing)=basis(A,A.positionid)

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
#  c:=PositionClass(G,G.identity);
#  res:=Sum(e{First(pprimesections(G,A.field.char),i->c in i)});
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
#  k:=InductionTable(B.group,A.group).scalar;
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
#
#GroupAlgebraOps.CharTable:=function(A)local conj;
#  if A.field.char>0 then 
#    Error("Cannot compute CharTable in positive characteristic");
#  fi;
#  A.charTable:=CharTable(A.group);
#  conj:=List(Elements(A.group),i->PositionClass(A.group,i));
#  A.charTable.basistraces:=List(A.charTable.irreducibles,chi->chi{conj});
#  A.charTable.Print:=function(table)Display(table);end;
#  return A.charTable; 
#end;

function PermRoot.radical(A::GroupAlgebra{T})where T
  if char(T)!=0 error("not implemented in positive characteristic") end
  TwoSidedIdeal(A,empty(basis(A)),Dict{Symbol,Any}())
end

Base.show(io::IO,A::GroupAlgebra{T}) where T =print(io,"Algebra(",A.group,",",T,")")

embedding(A::GroupAlgebra{T,E},g::E) where{T,E}=findfirst(==(g),A.elts)

"""
`GroupAlgebra(G,T=Int)` group algebra of `G` with coefficients `T`
"""
function GroupAlgebra(G,T=Int)
  A=GroupAlgebra(G,elements(G),1=>T(1),Dict{Symbol,Any}())
  A.gens=map(g->basis(A,embedding(A,g)),gens(G))
  if A.elts[1]!=one(G) error(one(G)," should be first in ",elements(G)) end
  A.showbasis=function(io::IO,e)
    fromTeX(io,"e"*"_{"*joindigits(word(A.group,A.elts[e]))*"}")
  end
#  if Q.char=0 then 
#    A.radical:=TwoSidedIdeal(A,[]);
#    A.radicalpowers:=[A.radical];
#  fi;
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
##   A.basename   = "X" (par defaut)
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
  rows=map(w->joindigits(w),ct.rows)
  irr=map(e->iszero(e) ? "." : repr(e),ct.irr)
  showtable(io,irr,row_labels=rows)
end

function Chars.CharTable(A::SolomonAlgebra)
  W=A.W
  conj=first.(W.solomon_conjugacy)
  orb=map(i->reflection_subgroup(W,W.solomon_subsets[i]),conj)
  cox=map(g->prod(gens(g);init=one(g)),orb)
  irr=[W.solomon_mackey[j,i][i] for i in conj, j in conj]
  SolomonCharTable(A,irr,map(w->word(W,w),cox),Dict{Symbol,Any}())
#  res.columns[Length(res.columns)]:=["0"];
#  res.basistraces:=List([1..Length(res.irreducibles)],function(ii)local r;
#    r:=Concatenation(List([1..Length(W.solomon.conjugacy)],
#      i->List(W.solomon.conjugacy[i], j-> [j,res.irreducibles[ii][i]])));
#    Sort(r); return List(r,i->i[2]);end);
#  if A.field.char=0 then return res;fi;
#  irr:=CyclotomicModP(res.irreducibles,A.field.char);
#  inc:=CollectBy([1..Length(irr)],irr);Sort(inc);
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
#     Sum(ReducedRightCosetRepresentatives(W, i),i-> B.embedding(i^-1)));
# fi;
  AlgebraHom(A,B,X)
end

Weyl.dim(A::SolomonAlgebra)=2^ngens(A.W)

const LIM=10000
"""
'SolomonAlgebra(W,K)'

Let  `(W,S)` be a  finite Coxeter group.  If `w` is  an element of `W`, let
`R(w)={s  ∈ S | l(ws) >  l(w)}`. If `I` is a  subset of `S`, we set
`Y_I={w ∈ W | R(w)=I}`, `X_I={w ∈ W | R(w) ⊃ I}`.

Note  that `X_I` is the set of minimal length left coset representatives of
`W/W_I`. Now, let `y_I=∑_{w ∈ Y_I} w`, `x_I=∑_{w ∈ X_I} w`.

They  are elements  of the  group algebra  `ℤ W`  of `W` over `Z`. Now, let
``Σ(W)  = ⊕_{I ⊂ S} ℤ y_I = ⊕_{I ⊂ S} ℤ x_I``. This is a sub-`ℤ`-module of
`ℤW`.  In fact, Solomon proved  that it is a  sub-algebra of `ℤW`. Now, let
`K(W)`  be the Grothendieck ring  of `W` and let  `θ:Σ(W)→ K(W)` be the map
defined  by  `θ(x_I)  =  Ind_{W_I}^W  1`.  Solomon  proved  that this is an
homomorphism of algebras. We call it the *Solomon homomorphism*.

returns  the Solomon  descent algebra  of the  finite Coxeter group `(W,S)`
over  `K`.  If  `S=[s₁,…,sᵣ]`,  the  element `x_I` corresponding to the
subset   `I=[s₁,s₂,s₄]`  of  `S`  is  printed  as  |X(124)|.  Note  that
'A:=SolomonAlgebra(W,K)' is endowed with the following fields:

'A.W': the group `W`

'A.basis': the basis `(x_I)_{I ⊂ S}`.

'A.xbasis':  the function sending the subset `I` (written as a number: for
instance `124` for `[s_1,s_2,s_4]`) to `x_I`.

'A.ybasis': the function sending the subset `I` to `y_I`.

'A.injection':  the injection of `A` in the group algebra of `W`, obtained
by calling 'SolomonAlgebraOps.injection(A)'.

Note that 'SolomonAlgebra(W,K)' endows `W` with the field `W.solomon` which
is a record containing the following fields:

'W.solomon_subsets': the set of subsets of `S`

'W.solomon_conjugacy':  conjugacy classes  of parabolic  subgroups of `W` (a
conjugacy   class  is  represented  by  the   list  of  the  positions,  in
'W.solomon.subsets', of the subsets `I` of `S` such that `W_I` lies in this
conjugacy class).

'W.solomon_mackey':  essentially  the  structure  constants  of  the Solomon
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
```

gap> injection(A)(X(123));
e(1) + e(2) + e(3) + e(8) + e(19) + e(45) + e(161) + e(361)
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
    for w in W
      leftascents=Idict[setdiff(1:r,leftdescents(W,w))]
      rightascents=Idict[setdiff(1:r,leftdescents(W,inv(w)))]
      l=map(i->action(W,i,w),1:r)
      for i in subsets[leftascents] 
        II=sort!(l[I[i]])
 	for j in subsets[rightascents]
          push!(mackey[i,j],Idict[Combinat.intersect_sorted(II,I[j])]=>1)
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
#     subspace=Subspace(A.vectorspace,A.EltToVector(k))
#     base=CanonicalBasis(subspace)
#     A.radical.operations.LeftTrace=x->TraceMat(List(base.vectors,v->
#       Coefficients(base,A.EltToVector(x*A.VectorToElt(v)))));
#     A.radical.operations.RightTrace=x->TraceMat(List(base.vectors,v->
#       Coefficients(base,A.EltToVector(A.VectorToElt(v)*x))));
#   end
    A.radicalpowers=[radical]
    if length(A.W)>LIM InfoAlgebra("done!\n") end
    radical
  end
  end
end

function permutation_character(W,S)
  t=InductionTable(S,W).scalar
  i=findfirst(x->all(==(1),x),eachrow(CharTable(S).irr))
  permutedims(CharTable(W).irr)*t[:,i]
end

#SolomonHomomorphism:=function(x) local A,i,W,T,irr,mat,WI,p;
#  A:=x.domain; W:=A.group; T:=CharTable(W); irr:=T.irreducibles;
#  mat:=List(irr, i-> 0);
#  for i in x.coefficients do 
#    if A.type="Solomon algebra" then 
#      WI:=ReflectionSubgroup(W,W.solomon.subsets[i[2]]);
#    elif A.type="Generalized Solomon algebra" then 
#      WI:=W.generalizedsolomon.subgroups[i[2]];
#    fi;
##   p:=PermutationCharacter(W,WI);
#    p:=PermutationCharacterByInductiontable(W,WI); # faster for large groups
#    mat:=mat+i[1]*MatScalarProducts(T,irr,[p])[1];
#  od;
#  return GrothendieckRing(W,A.field).VectorToElt(mat);
#end;
#
#SolomonAlgebraOps.xprimePrint:=function(r) local res,i,A,f,xprime;
#  A:=r.domain;
#  A.xprimebasisname:="X'";
#  res:=A.EltToVector(r);
#  xprime:=A.EltToVector(List(A.group.solomon.subsets,A.xprimebasis));
#  res:=A.VectorToElt(res*xprime^-1);
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
#  quotient:=space/Subspace(space,B.EltToVector(Radical(B).basis));
#  mat:=List(A.basis,i->B.EltToVector(SolomonHomomorphism(i)));
#  mat:=mat*ProjectionMatrix(quotient);
#  A.radical:=TwoSidedIdeal(A,A.VectorToElt(NullspaceMat(mat)));
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
