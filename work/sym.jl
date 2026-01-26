using Chevie
struct Partition
  l::Vector{Int}
  function Partition(l::Vector{<:Integer})
   if length(l)>0 && (!all(i->l[i]≥l[i+1],1:length(l)-1) || l[end]<=0)
      error(l," is not a partition")
    end
    new(l)
  end
end

Base.hash(p::Partition,k::UInt)=hash(p.l,k)
Base.:(==)(a::Partition,b::Partition)=a.l==b.l
Base.:*(a::Partition,b::Partition)=Partition(sort!(vcat(a.l,b.l),rev=true))

struct SymFunAlgebra
  n::Int
  charvalues::Dict{Pair{Partition,Partition},Int}
  centralizers::Dict{Partition,Int}
  partitions::Dict{Int,Vector{Partition}}
end

pos(l,p)=searchsorted(l,p)[1]

function charvalue(A::SymFunAlgebra,p::Pair{Partition,Partition})
  chi,c=p
  l=length(chi)
  if l!=length(c) || l>A.n return 0 end
  get!(A.charvalues,chi=>c)do 
    ct=CharTable(coxsym(l))
    pp=Partitions(A,l)
    for x in pp, y in pp A.charvalues[x=>y]=ct.irr[pos(pp,x),pos(pp,y)] end
    for x in pp A.centralizers[x]=ct.centralizers[pos(pp,x)] end
    A.charvalues[p]
  end
end

function centralizer(A::SymFunAlgebra,p::Partition)
  l=length(p)
  if l>A.n return 0 end
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
   Dict{Pair{Partition,Partition},Int}(),Dict{Partition,Int}(),
   Dict{Int,Vector{Partition}}())
                                 
function Base.show(io::IO, p::Partition)
  if hasdecor(io)
    print(io,isempty(p.l) ? "." : joindigits(p.l))
  else print(io,p.l)
  end
end

Partition(a...)=Partition(collect(a))
Base.length(p::Partition)=sum(p.l)
Base.isless(a::Partition,b::Partition)=length(a)<length(b) || a.l<b.l

struct SymFun{TS,C}
  d::ModuleElt{Partition,C}
  A::TS
  b::Symbol # :S,:P or :h depending on the basis
end

function Base.show(io::IO, h::SymFun)
  function showbasis(io::IO,w)
    res=string(h.b)
    if hasdecor(io) && !get(io,:naive,false) res*="_{"*xrepr(io,w)*"}"
    else  res*="("*xrepr(io,w)*")"
    end
    fromTeX(io,res)
  end
  show(rio(io,limit=true,showbasis=showbasis),improve_type(h.d))
end

clone(h::SymFun,d)=SymFun(d,h.A,h.b)

Base.:+(a::SymFun,b::SymFun)=clone(a,a.d+basis(a.b,b).d)
Base.:-(a::SymFun)=clone(a,-a.d)
Base.:-(a::SymFun, b::SymFun)=clone(a,a.d-basis(a.b,b).d)

Base.:*(a::SymFun, b::Union{Number,Pol,Mvp})=clone(a,a.d*b)
Base.:(//)(a::SymFun, b::Union{Number,Pol,Mvp})=clone(a,a.d//b)
Base.:*(b::Union{Number,Pol,Mvp}, a::SymFun)=a*b

Base.:^(a::SymFun, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                        Base.power_by_squaring(inv(a),-n)

Sbasis(H::SymFunAlgebra)=(x...)->basis(H,:S,x...)
Pbasis(H::SymFunAlgebra)=(x...)->basis(H,:P,x...)
hbasis(H::SymFunAlgebra)=(x...)->basis(H,:h,x...)

basis(H::SymFunAlgebra,t::Symbol,w::Vararg{Integer})=basis(H,t,collect(w))
basis(H::SymFunAlgebra,t::Symbol,w::Vector{<:Integer})=basis(H,t,Partition(w))
basis(H::SymFunAlgebra,t::Symbol,p::Partition)=SymFun(ModuleElt(p=>1),H,t)

function basis(H::SymFunAlgebra,t::Symbol,h::SymFun)
  if h.b==t return h end
  if h.b==:S
    if t==:P 
      sum(h.d)do (p,c)
        SymFun(ModuleElt(l=>charvalue(H,p=>l)*c//centralizer(H,l)
            for l in Partitions(H,length(p))),H,:P)
      end
    else error(":S=>:h not done")
    end
  elseif h.b==:P
    if t==:S 
      sum(h.d)do (p,c)
        SymFun(ModuleElt(l=>charvalue(H,l=>p)*c for l in
                         Partitions(H,length(p))),H,:S)
      end
    else error(":P=>:h not done")
    end
  elseif h.b==:h
    if t==:S basis(:S,basis(:P,h))
    else
      sum(h.d)do (p,c)
        prod(i->basis(:P,Sbasis(h.A)(i)),p.l)*c
      end
    end
  end
end

basis(t::Symbol,h::SymFun)=basis(h.A,t,h)

Base.getindex(a::SymFun,p::Partition)=a.d[p]
Base.getindex(a::SymFun,x::Vararg{Int})=a.d[Partition(collect(x))]

function Base.:*(a::SymFun,b::SymFun)
  H=a.A
  a1=basis(:P,a)
  b1=basis(:P,b)
  basis(a.b,SymFun(ModuleElt([p*q=>c*c1 for (p,c) in a1.d for (q,c1) in b1.d]),
                     H,:P))
end

function PuiseuxPolynomials.Mvp(a::SymFun,n=length(last(a.d.d)[1]))
  f=basis(:P,a)
  sum(a.d)do (l,c)
    prod(i->sum(j->Mvp(Symbol("x",stringind(rio(),j)))^i,1:n),l.l)*c
  end
end

function scalarproduct(a::SymFun,b::SymFun)
  b=basis(a.b,b)
  if a.b==:S
    sum(values(ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d)))
  elseif a.b==:P
   sum(((p,c),)->centralizer(a.A,p)*c,
        ModuleElts.merge2((x,y)->x*conj(y),a.d,b.d))
  else scalarproduct(basis(:P,a),basis(:P,b))
  end
end
