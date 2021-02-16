module Algebras
using ..Gapjm
export FiniteDimAlgebra, AlgebraElt, basis, dim, idempotents, iscommutative
abstract type FiniteDimAlgebra end

struct AlgebraElt{T,TA}
  A::TA
  d::ModuleElt{Int,T}
end

dim(A::FiniteDimAlgebra)=error("not implemented in general")
idempotents(A::FiniteDimAlgebra)=error("not implemented in general")
iscommutative(A::FiniteDimAlgebra)=error("not implemented in general")

function basis(A::FiniteDimAlgebra)
  get!(A,:basis)do
    map(i->AlgebraElt(A,ModuleElt(i=>1)),1:dim(A))
  end
end

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

function Cycs.coefficients(a::AlgebraElt{T})where T
  v=fill(T(0),dim(a.A))
  for (i,c) in a.d v[i]=c end
  v
end

function Base.:*(a::AlgebraElt{T1}, b::AlgebraElt{T2})where {T1,T2}
  res=Pair{Int,promote_type(T1,T2)}[]
  for (i,c) in a.d
    for (i1,c1) in b.d
     l=(i>=i1 || !iscommutative(a.A)) ? a.A.multable[i][i1] : a.A.multable[i1][i]
      append!(res,[k=>c*c1*c2 for (k,c2) in l])
    end
  end
  AlgebraElt(a.A,ModuleElt(res;check=true))
end

end
