module Algebras
using ..Gapjm
export FiniteDimAlgebra, AlgebraElt, basis, dim, idempotents,
       Quaternions, involution, ZeroHecke
abstract type FiniteDimAlgebra end

# properties of Algebras: have .multable, dim

# Algebra element of algebra A
# d coefficients on basis indexed by number
# basis[1] should be identity element
struct AlgebraElt{TA,T}
  A::TA
  d::ModuleElt{Int,T}
end

AlgebraElt(A::FiniteDimAlgebra,v::AbstractVector)=
     AlgebraElt(A,ModuleElt(Pair.(1:dim(A),v)))

dim(A::FiniteDimAlgebra)=error("not implemented in general")
idempotents(A::FiniteDimAlgebra)=error("not implemented in general")

"`isabelian(A::FiniteDimAlgebra)` whether `A` is commutative"
function Groups.isabelian(A::FiniteDimAlgebra)
  get!(A,:isabelian)do
    l=haskey(A,:generators) ? A.generators : basis(A)
    all(x*y==y*x for x in l, y in l)
  end::Bool
end

function isassociative(A::FiniteDimAlgebra)
  get!(A,:isassociative)do
    l=basis(A)
    all((x*y)*z==x*(y*z) for x in l, y in l, z in l)
  end::Bool
end

function basis(A::FiniteDimAlgebra)
  get!(A,:basis)do
    map(i->AlgebraElt(A,ModuleElt(i=>1)),1:dim(A))
  end
end
basis(A::FiniteDimAlgebra,i::Integer)=basis(A)[i]

one(A::FiniteDimAlgebra)=basis(A,1)

function Base.show(io::IO, h::AlgebraElt)
  function showbasis(io::IO,e)
    repl=get(io,:limit,false)
    TeX=get(io,:TeX,false)
    res="B"
    if repl || TeX
      res*="_{"*string(e)*"}"
    else res*="("*string(e)*")"
    end
    fromTeX(io,res)
  end
  show(IOContext(io,:showbasis=>showbasis),h.d)
end

Base.:+(a::AlgebraElt,b::AlgebraElt)=AlgebraElt(a.A,a.d+b.d)
Base.:-(a::AlgebraElt)=AlgebraElt(a.A,-a.d)
Base.:-(a::AlgebraElt, b::AlgebraElt)=a+(-b)
Base.:*(a::AlgebraElt, b::Vector)=Ref(a).*b #for tbl
Base.:*(a::AlgebraElt, b)=AlgebraElt(a.A,a.d*b)
Base.:*(b,a::AlgebraElt)=AlgebraElt(a.A,a.d*b)
Base.zero(a::AlgebraElt)=AlgebraElt(a.A,zero(a.d))
Base.:^(a::AlgebraElt, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                          Base.power_by_squaring(inv(a),-n)
Base.isless(a::AlgebraElt,b::AlgebraElt)=isless(a.d,b.d)
Base.:(==)(a::AlgebraElt,b::AlgebraElt)=a.d==b.d
Base.hash(a::AlgebraElt,k::UInt)=hash(a.d,k)
Base.one(a::AlgebraElt)=one(a.A)

function Gapjm.coefficients(a::AlgebraElt{A,T})where {A,T}
  v=fill(T(0),dim(a.A))
  for (i,c) in a.d v[i]=c end
  v
end

Base.Matrix(a::AlgebraElt)=toM(map(coefficients,Ref(a).*basis(a.A)))

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

basismult(A::FiniteDimAlgebra,i::Integer,j::Integer)=A.multable[i][j]

function Base.:*(a::AlgebraElt{A,T1}, b::AlgebraElt{A,T2})where {A,T1,T2}
  res=Pair{Int,promote_type(T1,T2)}[]
  for (i,c) in a.d
    for (i1,c1) in b.d
      l=(i>=i1 || !isabelian(a.A)) ? basismult(a.A,i,i1) : basismult(a.A,i1,i)
      append!(res,[k=>c*c1*c2 for (k,c2) in l])
    end
  end
  AlgebraElt(a.A,ModuleElt(res))
end

#----------------------- Quaternions -----------------------------------
@GapObj struct Quaternions{T}<:FiniteDimAlgebra
  multable::Vector{Vector{Vector{Pair{Int,T}}}}
end

function Quaternions(a::T=-1,b::T=-1)where T
  multable=
   [[[1=>1],[2=>1], [3=>1],[4=>1]],
    [[2=>1],[1=>a], [4=>1],[3=>-1]],
    [[3=>1],[4=>-1],[1=>b],[2=>-b]],
    [[4=>1],[3=>-a],[2=>b],[1=>-a*b]]]
  A=Quaternions(multable,Dict{Symbol,Any}())
  A.a=a;A.b=b
  A.isabelian=false
  A
end

dim(A::Quaternions)=4

Base.show(io::IO,A::Quaternions{T}) where T=print(io,"Quaternions(",A.a,",",A.b,")")

function Base.show(io::IO,h::AlgebraElt{<:Quaternions})
  showbasis(io::IO,e)=["","i","j","k"][e]
  show(IOContext(io,:showbasis=>showbasis),h.d)
end

involution(h::AlgebraElt{<:Quaternions})=
  AlgebraElt(h.A,coefficients(h).*[1,-1,-1,-1])

LinearAlgebra.norm(h::AlgebraElt{<:Quaternions})=h*involution(h)
#----------------------- 0-Hecke -----------------------------------
@GapObj struct ZeroHecke{TE,TW<:CoxeterGroup{TE}}<:FiniteDimAlgebra
  W::TW
  e::Vector{TE}
end

Base.show(io::IO,A::ZeroHecke)=print(io,"0-Hecke(",A.W,")")

dim(A::ZeroHecke)=length(A.W)

function basismult(A::ZeroHecke,i::Integer,j::Integer)
  res=A.e[j]
  for k in reverse(word(A.W,A.e[i]))
    if !isleftdescent(A.W,res,k) res=A.W(k)*res end
  end
  [findfirst(==(res),A.e)=>1]
end
  
function Tbasis(A::ZeroHecke)
  return function(v::Int...)
    prod(i->A.generators[i],v;init=basis(A,1))
  end
end

function Base.show(io::IO,h::AlgebraElt{<:ZeroHecke})
  showbasis(io::IO,i)=fromTeX(io,string("T_{",joindigits(word(h.A.W,h.A.e[i])),"}"))
  show(IOContext(io,:showbasis=>showbasis),h.d)
end

# 0-Hecke algebra of W with parameters 0,1
function ZeroHecke(W::CoxeterGroup)
  e=elements(W);# we use they are stored by increasing CoxeterLength
  A=ZeroHecke(W,e,Dict{Symbol,Any}())
  A.generators=basis(A)[2:ngens(W)+1]
  A.isabelian=false
  A
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
#  c:=List(CollectBy(List([1..d],i->[i,A.mots[i]]),x->Set(x[2])),
#     x->List(x,x->x[1]));
#  A.radical:=List(c,i->List([2..Length(i)],j->A.basis[i[1]]-A.basis[i[j]]));
#  A.radical:=TwoSidedIdeal(A,Concatenation(A.radical));
#  A.radicalpowers:=[A.radical];
#  A.embedding:=function(g) return A.basis[Position(e,g)];end;
end
