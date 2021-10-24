@GapObj struct CPMonoid{T}<:GarsideMonoid{T}
  e::Int
  n::Int
  k::Int
  δ::T
  orderδ::Int
  stringδ::String
  atoms::Vector{T}
  W::Group{T}
  one::T
end

Garside.IntervalStyle(M::CPMonoid)=Garside.Interval()

function Base.show(io::IO,M::CPMonoid)
  print(io,"CorranPicantinMonoid(",M.e,",",M.e,",",M.n,")")
  if M.k!=1 print("[",M.k,"]") end
end

# for the next function, see Proposition 3.11 in [1]
function Garside.isleftdescent(M::CPMonoid,simp,i)
  m=mat(M,simp)
  if i<=M.e
    c1=findfirst(!iszero,@view m[1,:])
    c2=findfirst(!iszero,@view m[2,:])
    if c1<c2 return m[2,c2]!=1
    else return m[1,c1]==E(M.e,1-i)
    end
  end
  i+=2-M.e
  c1=findfirst(!iszero,@view m[i-1,:])
  c2=findfirst(!iszero,@view m[i,:])
  if c1<c2 m[i,c2]!=1
  else m[i-1,c1]==1
  end
end

mat(M::CPMonoid,s)=reflrep(M.W,s)
showmat(M::CPMonoid,s::Perm)=display(mat(M,s))
showmat(M::CPMonoid,s::GarsideElt)=for t in s.elm showmat(M,t) end

"""
© July 2017 --- Jean Michel and Georges Neaime

The  function `CorranPicantinMonoid(e,n,k=1)` computes  the interval monoid
generalizing   the  Corran-Picantin  monoid,  as   defined  in  G.  Neaime,
http://arxiv.org/abs/1707.06864.

It  returns the  monoid with  δ the  diagonal matrix whose diagonal entries
except  the first are  equal to E(e)^k;  this monoid is the Corran-Picantin
monoid for G(e,e,n) when Gcd(k,e)=1
"""
function CorranPicantinMonoid(e,n,k=1;revMonoid=nothing)
  k%=e
  if k==0 error("k must not be divisible by e") end
  r=fill(Cyc(0),n,n)
  r[1,1:2]=[-E(e),1]
  for i in 2:n r[i,i-1:i]=[-1,1] end
  cr=copy(r)
  cr[1,1]=-E(e,-1)
  W=PRG(r,cr)
  atoms=gens(W)[[2,1]]
  for i in 2:e-1 push!(atoms,atoms[i-1]^atoms[i]) end
  append!(atoms,gens(W)[3:end])
  δmat=reflrep(W,one(W))*E(e)^k
  δmat[1,1]=E(e)^(k*(1-n))
  δ=PermX(W,δmat)
  M=CPMonoid(e,n,k,δ,order(δ),"δ",atoms,W,one(W),Dict{Symbol,Any}())
  if k!=e-k 
    if revMonoid===nothing 
      M.revMonoid=CorranPicantinMonoid(e,n,e-k;revMonoid=M)
    else M.revMonoid=revMonoid
    end
  end
  M
end
