using Chevie

struct SymFunAlgebra
  n::Int
  charvalues::Dict{Pair{Partition,Partition},Int}
  centralizers::Dict{Partition,Int}
  partitions::Dict{Int,Vector{Partition}}
end

function charvalue(A::SymFunAlgebra,p::Pair{Partition,Partition})
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

function Groups.centralizer(A::SymFunAlgebra,p::Partition)
  if size(p)>A.n return 0 end
  get!(A.centralizers,p)do 
    charvalue(A,p=>p)
    A.centralizers[p]
  end
end

function Combinat.Partitions(A::SymFunAlgebra,n)
  if n>A.n return Partition[] end
  get!(A.partitions,n)do
    Partition.(partitions(n))
  end
end

function Base.show(io::IO,A::SymFunAlgebra)
  print(io,"SymFunAlgebra(",A.n,")")
end

SymFunAlgebra(n)=SymFunAlgebra(n,
   Dict{Pair{Partition,Partition},Int}((Partition()=>Partition())=>1),
   Dict{Partition,Int}(Partition()=>1),
   Dict{Int,Vector{Partition}}())
                                 
struct SymFun{b,TS,C}# b=:S,:P or :h depending on the basis
  d::ModuleElt{Partition,C}
  A::TS
  SymFun{b}(m,A)where b =new{b,typeof(A),valtype(m)}(m,A)
end

function Base.show(io::IO, h::SymFun{b})where b
  function showbasis(io::IO,w)
    res=string(b)
    if hasdecor(io) && !get(io,:naive,false) res*="_{"*xrepr(io,w)*"}"
    else  res*="("*xrepr(io,w)*")"
    end
    fromTeX(io,res)
  end
  show(rio(io,limit=true,showbasis=showbasis),improve_type(h.d))
end

clone(h::SymFun{b},d) where b=SymFun{b}(d,h.A)

Base.:+(a::SymFun{ba},b::SymFun) where ba=clone(a,a.d+basis(Val(ba),b).d)
Base.:-(a::SymFun)=clone(a,-a.d)
Base.:-(a::SymFun{ba},b::SymFun) where ba=clone(a,a.d-basis(Val(ba),b).d)
Base.:*(a::SymFun, b::Union{Number,Pol,Mvp})=clone(a,a.d*b)
Base.:(//)(a::SymFun, b::Union{Number,Pol,Mvp})=clone(a,a.d//b)
Base.:*(b::Union{Number,Pol,Mvp}, a::SymFun)=a*b

Base.:^(a::SymFun, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                      Base.power_by_squaring(inv(a),-n)
Base.one(a::SymFun{b})where b=clone(a,ModuleElt(Partition()=>1))
Sbasis(H::SymFunAlgebra)=(x...)->basis(H,Val(:S),x...)
Pbasis(H::SymFunAlgebra)=(x...)->basis(H,Val(:P),x...)
hbasis(H::SymFunAlgebra)=(x...)->basis(H,Val(:h),x...)

Algebras.basis(h::SymFun{b})where b=b
Algebras.basis(H::SymFunAlgebra,::Val{b},p::Partition)where b=
  SymFun{b}(ModuleElt(p=>1),H)
Algebras.basis(H::SymFunAlgebra,::Val{b},w::Vector{<:Integer})where b=
  basis(H,Val(b),Partition(w))
Algebras.basis(H::SymFunAlgebra,::Val{b},w::Vararg{Integer})where b=
  basis(H,Val(b),collect(w))
Algebras.basis(::Val{b},h::SymFun)where b=basis(h.A,Val(b),h)

function Algebras.basis(H::SymFunAlgebra,::Val{:P},h::SymFun{:S})
   sum(h.d)do (p,c)
       SymFun{:P}(ModuleElt(l=>charvalue(H,p=>l)*c//centralizer(H,l)
            for l in Partitions(H,size(p))),H)
   end
end
Algebras.basis(H::SymFunAlgebra,::Val{:S},h::SymFun{:S})=h
Algebras.basis(H::SymFunAlgebra,::Val{:S},h::SymFun{:P})=
  sum(h.d)do (p,c)
    SymFun{:S}(ModuleElt(l=>charvalue(H,l=>p)*c for l in
                     Partitions(H,size(p))),H)
  end
Algebras.basis(H::SymFunAlgebra,::Val{:P},h::SymFun{:P})=h
Algebras.basis(H::SymFunAlgebra,::Val{:S},h::SymFun{:h})=basis(Val(:S),basis(Val(:P),h))
Algebras.basis(H::SymFunAlgebra,::Val{:P},h::SymFun{:h})=
sum(((p,c),)->reduce(⊗,map(i->basis(Val(:P),basis(H,Val(:S),i)),p.l);init=Pbasis(H)())*c,h.d)
Algebras.basis(H::SymFunAlgebra,::Val{:h},h::SymFun{:h})=h
function Algebras.basis(H::SymFunAlgebra,::Val{:h},h::SymFun{:P})
  function Ptoh(n)
    SymFun{:h}(ModuleElt(map(Partitions(H,n))do p
      p=>n*(-1)^(length(p)-1)*factorial(length(p)-1)//
      prod(i->factorial(i[2]),tally(p.l))
     end),H)
  end
  sum(((p,c),)->reduce(⊗,map(Ptoh,p.l);init=hbasis(H)())*c,h.d)
end
Algebras.basis(H::SymFunAlgebra,::Val{:h},h::SymFun{:S})=basis(Val(:h),basis(Val(:P),h))

Base.getindex(a::SymFun,p::Partition)=a.d[p]
Base.getindex(a::SymFun,x::Vararg{Int})=a.d[Partition(collect(x))]

function FinitePosets.:⊗(a::SymFun{ba},b::SymFun)where ba
  b=basis(Val(ba),b)
  if ba in (:P,:h)
    clone(a,ModuleElt([p+q=>c*c1 for (p,c) in a.d for (q,c1) in b.d]))
  else
    basis(Val(ba),basis(Val(:P),a)⊗basis(Val(:P),b))
  end
end

function Base.:*(a::SymFun{ba},b::SymFun)where ba
  H=a.A
  a1=basis(Val(:P),a)
  b1=basis(Val(:P),b)
  m=ModuleElts.merge2((x,y)->x*y,a1.d,b1.d)
  m=ModuleElt(p=>c*centralizer(H,p) for (p,c) in m)
  basis(Val(ba),SymFun{:P}(m,H))
end

function PuiseuxPolynomials.Mvp(a::SymFun,n=size(last(a.d.d)[1]))
  f=basis(Val(:P),a)
  sum(f.d)do (l,c)
    prod(i->sum(j->Mvp(Symbol("x",stringind(rio(),j)))^i,1:n),l.l)*c
  end
end

function scalarproduct(a::SymFun{ab},b::SymFun)where ab
  b=basis(Val(ab),b)
  if ab==:S
    sum(values(ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d)))
  elseif ab==:P
    sum(((p,c),)->centralizer(a.A,p)*c,
        ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d))
  else scalarproduct(basis(Val(:P),a),basis(Val(:P),b))
  end
end
