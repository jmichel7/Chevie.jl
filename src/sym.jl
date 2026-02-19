using .SymFuncs: SymFunc

struct PartitionTuple
  P::Vector{Vector{Int}}
  function PartitionTuple(P::Vector{<:Vector{<:Integer}})
    for p in P
      if !all(i->p[i]≥p[i+1],1:length(p)-1)
        error(p," not a partition: parts should be non-increasing")
       elseif length(p)>0 && p[end]<=0
        error(p," not a partition: parts should be positive")
      end
    end
    new(P)
  end
end

Base.:(==)(a::PartitionTuple,b::PartitionTuple)=a.P==b.P

PartitionTuple(v::Vector{Partition})=PartitionTuple(map(x->x.l,v))

PartitionTuple(v::Vararg{Vector{<:Integer}})=PartitionTuple(collect(v))

Base.isless(p1::PartitionTuple,p2::PartitionTuple)=p1.P<p2.P

function Base.show(io::IO,p::PartitionTuple)
  join(io,map(x->isempty(x) ? "" : joindigits(x),p.P),".")
end

struct WreathSymFunc{b,C}# b=:s,:p or :h depending on the basis
  d::ModuleElt{PartitionTuple,C}
  A::SymFuncAlgebra
  WreathSymFunc{b}(m,A)where b =new{b,valtype(m)}(m,A)
end

function Base.show(io::IO, h::WreathSymFunc{b})where b
  function showbasis(io::IO,w)
    res=string(b)
    s=xrepr(io,w)
    if hasdecor(io) && !get(io,:naive,false) res*="_{"*s*"}"
    else  res*="("*s*")"
    end
    fromTeX(io,res)
  end
  show(rio(io,limit=true,showbasis=showbasis),improve_type(h.d))
end

Algebras.basis(H::SymFuncAlgebra,::Val{b},p::PartitionTuple)where b=
  WreathSymFunc{b}(ModuleElt(p=>1),H)
Algebras.basis(H::SymFuncAlgebra,::Val{b},w::Vector{<:Vector{<:Integer}})where b=
  basis(H,Val(b),PartitionTuple(w))
Pibasis(H::SymFuncAlgebra)=(x...)->basis(H,Val(:π),x...)
Algebras.basis(h::WreathSymFunc{b})where b=b
Algebras.basis(H::SymFuncAlgebra,::Val{b},h::WreathSymFunc)where b=basis(Val(b),h)
clone(h::WreathSymFunc{b},d) where b=WreathSymFunc{b}(d,h.A)

Base.:+(a::WreathSymFunc{ba},b::WreathSymFunc) where ba=clone(a,a.d+basis(Val(ba),b).d)
Base.:+(a::WreathSymFunc, b::Union{Number,Pol,Mvp})=a+one(a)*b
Base.:-(a::WreathSymFunc)=clone(a,-a.d)
Base.:-(a::WreathSymFunc{ba},b::WreathSymFunc) where ba=clone(a,a.d-basis(Val(ba),b).d)
Base.:*(a::WreathSymFunc, b::Union{Number,Pol,Mvp,Frac})=clone(a,a.d*b)
Base.:(//)(a::WreathSymFunc, b::Union{Number,Pol,Mvp})=clone(a,a.d//b)
Base.:*(b::Union{Number,Pol,Mvp,Frac}, a::WreathSymFunc)=a*b

⊠(a::SymFunc{b},a1::SymFunc{b1})where {b,b1}=tens(a,a1)

function tens(aa::SymFunc...)
  if length(aa)<2 error() end
  A=aa[1].A
  b=basis(aa[1])
  for a in aa[2:end] a=basis(A,Val(b),a) end
  m=ModuleElt(map(cartesian(map(x->pairs(x.d),aa)...))do pp
               PartitionTuple(first.(pp))=>prod(last.(pp))
              end)
  WreathSymFunc{b}(m,A)
end

function Algebras.basis(::Val{b},h::WreathSymFunc{b1})where {b,b1}
  if b==b1 return h end
  if b==:π || b1==:π error("uuuu") end
  sum(h.d)do (P,c)
   tens(map(p->basis(Val(b),basis(A,Val(b1),p)),P.P)...)*c
  end
end

LaurentPolynomials.degree(P::PartitionTuple)=length(P.P)
LaurentPolynomials.degree(p::WreathSymFunc)=degree(first(first(p.d)))
    
Algebras.basis(::Val{:π},h::WreathSymFunc{:π})=h

function Algebras.basis(::Val{:π},h::WreathSymFunc{b})where b
  function toπ(A,n,j,e)
    if n==0 WreathSymFunc{:π}(ModuleElt(PartitionTuple(fill(Int[],e))=>1),A)
    elseif e==2 WreathSymFunc{:π}(
        ModuleElt(PartitionTuple([n],Int[])=>1//2,
                  PartitionTuple(Int[],[n])=>j==1 ? 1//2 : -1//2),A)
    else
      WreathSymFunc{:π}(ModuleElt(
       map(1:e)do j
          v=fill(Int[],e);v[j]=[n]
          PartitionTuple(v)=>E(e,j-1)//e
       end),A)
    end
  end
  h=basis(Val(:p),h)
  e=degree(h)
  sum(h.d)do(P,c)
    prod(eachindex(P.P))do j
      prod(eachindex(P.P[j]);init=toπ(A,0,j,e))do i
        n=P.P[j][i]
        toπ(A,P.P[j][i],j,e)
      end
    end*c
  end
end
  
function Algebras.basis(::Val{b},h::WreathSymFunc{:π})where b
  function fromπ(A,n,i,e)
    if n==0 WreathSymFunc{:p}(ModuleElt(PartitionTuple(fill(Int[],e))=>1),A)
    elseif e==2 
       WreathSymFunc{:p}(
        ModuleElt(PartitionTuple([n],Int[])=>1,
                  PartitionTuple(Int[],[n])=>i==1 ? 1 : -1),A)
    else
        WreathSymFunc{:p}(ModuleElt(
         map(1:e)do j
            v=fill(Int[],e);v[j]=[n]
            PartitionTuple(v)=>E(e,j-1)//e
         end),A)
    end
  end
  e=degree(h)
  basis(Val(b),sum(h.d)do(P,c)
    prod(eachindex(P.P))do j
      prod(eachindex(P.P[j]);init=fromπ(A,0,j,e))do i
        fromπ(A,P.P[j][i],j,e)
      end
    end*c
  end)
end

function Base.:*(a::WreathSymFunc{b},a1::WreathSymFunc)where {b}
  a1=basis(Val(b),a1)
  A=a.A
  sum(cartesian(pairs(a.d),pairs(a1.d)))do ((P1,c1),(P2,c2))
   tens(map((p1,p2)->basis(A,Val(b),p1)*basis(A,Val(b),p2),P1.P,P2.P)...)*c1*c2
     end
end

function Base.:*(a::WreathSymFunc{:π},a1::WreathSymFunc{:π})
  WreathSymFunc{:π}(ModuleElt(
    map(cartesian(pairs(a.d),pairs(a1.d)))do ((P1,c1),(P2,c2))
     PartitionTuple(map(union,Partition.(P1.P),Partition.(P2.P)))=>c1*c2
    end),a.A)
end
