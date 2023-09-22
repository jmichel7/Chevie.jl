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

Base.Matrix(M::CPMonoid,s)=reflrep(M.W,s)
showmat(M::CPMonoid,s)=display(Matrix(M,s))
function showmat(s::GarsideElt)
  for i in 1:s.pd showmat(s.M,s.M.δ) end
  for t in s.elm showmat(s.M,t) end
end

"""
`CorranPicantinMonoid(e,n,k=1)`

returns the interval monoid defined by G. Neaime,
http://arxiv.org/abs/1707.06864,   which  generalizes  the  Corran-Picantin
monoid for `G(e,e,n)`.

In  this monoid `δ` has `image`  the element of `G(e,e,n)` corresponding to
the  diagonal matrix whose  diagonal entries except  the first are equal to
`E(e)^k`;  this  monoid  is  isomorphic  to  the Corran-Picantin monoid for
`G(e,e,n)` when `gcd(k,e)=1`.

```julia-repl
julia> C=CorranPicantinMonoid(3,3)
CorranPicantinMonoid(3,3,3)

julia> word(C(C.δ))
6-element Vector{Int64}:
 1
 3
 4
 1
 3
 4

julia> Matrix(C,C.δ)
3×3 Matrix{Cyc{Int64}}:
 ζ₃   0   0
  0  ζ₃   0
  0   0  ζ₃

julia> b=C(1,2,3,4)^3
1.2.341.2.341.2.34

julia> Matrix(C,b[3])
3×3 Matrix{Cyc{Int64}}:
 0    0  ζ₃
 0  ζ₃²   0
 1    0   0
```
© July 2017 --- Jean Michel and Georges Neaime
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

export CorranPicantinMonoid
