"""
A  signed permutation of `1:n` is  a permutation of the set `-n,…,-1,1,…,n`
which  preserves the  pairs `(-i,i)`.  It is  represented internally as the
images of `1:n`. It is printed as a product of signed cycles.

# Examples
```julia-repl
julia> SPerm([-2,-1,-3])
SPerm{Int64}: (1,-2)(3,-3)

julia> p=SPerm(-1)
(1,-1)

julia> q=SPerm(1,2)
(1,2)

julia> sort(elements(Group([p,q])))
8-element Vector{SPerm{Int16}}:
 (1,-2)
 (1,-2,-1,2)
 (1,-1)(2,-2)
 (1,-1)
 (2,-2)
 ()
 (1,2,-1,-2)
 (1,2)
```

The  complete type of signed permutations is `SPerm{T}` where `T<:Integer`,
where `Vector{T}` is the type of the vector which holds the image of `1:n`.
This  can used to save space or time when possible. If `T` is not specified
we  take it to be `Int16` since this is a good compromise between speed and
compactness.

SPerms  have methods `copy, hash, ==, isless`  (total order) so they can be
keys in hashes or elements of sets; two `SPerms` are equal if they move the
same points to the same images. For instance,
```julia-repl
julia> SPerm([-2,-1,-3])==SPerm([-2,-1,-3,4])
true
```
SPerms are considered as scalars for broadcasting.
"""
module SPerms
using ..Gapjm
export SPerm, CoxHyperoctaedral, sstab_onmats, SPerm_onmats, @sperm_str, signs

"""
`struct SPerm`

An  `SPerm` represents a signed permutation of `1:n`, that is a permutation
of  the  set  `-n,…,-1,1,…,n`  which  preserves  the  pairs `(-i,i)`. It is
implemented  by a `struct SPerm`  with one field `d`,  a vector holding the
images of `1:n`.
"""
struct SPerm{T<:Integer}
   d::Vector{T}
end

"""
SPerm{T}(x::Integer...)where T<:Integer

returns   a   signed   cycle.   For  instance  `SPerm{Int8}(1,2,-1,2)`  and
`SPerm({Int8}[1,-2])`  define  the  same  signed  permutation. If not given
`{T}` is taken to be `{Int16}`.
"""
function SPerm{T}(x::Integer...)where T<:Integer
  if isempty(x) return SPerm(T[]) end
  d=T.(1:max(abs.(x)...))
  for i in 1:length(x)-1
   d[abs(x[i])]=sign(x[i])*x[i+1]
  end
  if length(x)==1 d[abs(x[1])]=x[1]
  else d[abs(x[end])]=sign(x[end])*x[1]
  end
  SPerm(d)
end

const Idef=Int16 # you can change the default type T for SPerm here

SPerm(x::Integer...)=SPerm{Idef}(x...)

"""
   @sperm"..."

 make a `SPerm` from a string; allows style `sperm"(1,-2)(5,-6,7)(-4,9)"`
"""
macro sperm_str(s::String)
  start=1
  res=SPerm()
  while true
    m=match(r"\((\s*-?\d+\s*,)+\s*-?\d+\)",s[start:end])
    if isnothing(m) break end
    start+=m.match.ncodeunits
    res*=SPerm(Meta.parse(m.match).args...)
  end
  res::SPerm
end

Base.convert(::Type{SPerm{T}},p::SPerm{T1}) where {T,T1}=T==T1 ? p : SPerm(T.(p.d))

@GapObj struct SPermGroup{T}<:Group{SPerm{T}}
  gens::Vector{SPerm{T}}
end

function Base.show(io::IO,G::SPermGroup)
  print(io,"Group([")
  join(io,gens(G),',')
  print(io,"])")
end

function Groups.Group(a::AbstractVector{SPerm{T}}) where T
  SPermGroup(filter(!isone,a),Dict{Symbol,Any}())
end

SPermGroup()=Group(SPerm{Idef}[])

Base.one(p::SPerm)=SPerm(empty(p.d))
Base.one(::Type{SPerm{T}}) where T=SPerm(T[])
Base.copy(p::SPerm)=SPerm(copy(p.d))

Base.vec(a::SPerm)=a.d

"`Perm(p::SPerm)` returns the underlying Perm of an SPerm"
Perms.Perm(p::SPerm)=Perm(abs.(p.d))
SPerm(p::Perm)=SPerm(vec(p))
signs(p::SPerm)=sign.(p.d)

# SPerms are scalars for broadcasting"
Base.broadcastable(p::SPerm)=Ref(p)

# hash is needed for using SPerms in Sets/Dicts
function Base.hash(a::SPerm, h::UInt)
  for (i,v) in enumerate(a.d)
    if v!=i h=hash(v,h) end
  end
  h
end

function Base.promote_rule(a::Type{SPerm{T1}},b::Type{SPerm{T2}})where {T1,T2}
  SPerm{promote_type(T1,T2)}
end

extend!(a::SPerm,n::Integer)=if length(a.d)<n append!(a.d,length(a.d)+1:n) end

"""
 `promote_degree(a::SPerm, b::SPerm)` promotes `a` and `b` to the same type,
 then extends `a` and `b` to the same degree
"""
function promote_degree(a::SPerm,b::SPerm)
  a,b=promote(a,b)
  extend!(a,length(b.d))
  extend!(b,length(a.d))
  (a,b)
end

# total order is needed to use SPerms in sorted lists
function Base.isless(a::SPerm, b::SPerm)
  a,b=promote_degree(a,b)
  for (ai,bi) in zip(a.d,b.d) ai!=bi && return ai<bi end
  false
end

function Base.:(==)(a::SPerm, b::SPerm)
  a,b=promote_degree(a,b)
  a.d==b.d
end

"""
  orbit(a::SPerm,i::Integer) returns the orbit of a on i
"""
function Groups.orbit(a::SPerm{T},i::Integer)where T
  if abs(i)>length(a.d) return T[i] end
  res=T[]
  sizehint!(res,length(a.d))
  j=i
  while true
    push!(res,j)
@inbounds j=sign(j)*a.d[abs(j)]
    if j==i return res end
  end
end

function Perms.cycles(p::SPerm)
  cycles=Vector{eltype(p.d)}[]
  to_visit=trues(length(p.d))
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    cyc=orbit(p,i)
@inbounds to_visit[abs.(cyc)].=false
    if length(cyc)>1 push!(cycles,cyc) end
  end
  cycles
end

"""
order(a)

order of the signed permutation a
"""
Gapjm.order(a::SPerm) = lcm(length.(cycles(a)))

function Base.show(io::IO, a::SPerm)
  replorTeX=get(io,:limit,false)||get(io,:TeX,false)
  if !replorTeX print(io,"sperm\"") end
  cyc=cycles(a)
  if isempty(cyc) print(io,"()")
  else for c in cyc print(io,"(",join(c,","),")") end
  end
  if !replorTeX print(io,"\"") end
end

function Base.show(io::IO, ::MIME"text/plain", p::SPerm{T})where T
  if T!=Idef && !haskey(io,:typeinfo) print(io,typeof(p),": ") end
  show(io,p)
end

function Base.:*(a::SPerm, b::SPerm)
  a,b=promote_degree(a,b)
  r=similar(a.d)
  for (i,v) in enumerate(a.d)
@inbounds if v<0 r[i]=-b.d[-v] else r[i]=b.d[v] end
  end
  SPerm(r)
end

function Base.inv(a::SPerm)
  r=similar(a.d)
  for (i,v) in enumerate(a.d)
 @inbounds if v<0 r[-v]=-i  else r[v]=i end
  end
  SPerm(r)
end

# less allocations than inv(a)*b
function Base.:\(a::SPerm, b::SPerm)
  a,b=promote_degree(a,b)
  r=similar(a.d)
  @inbounds for (i,v) in enumerate(a.d) 
     if v<0 r[-v]=-b.d[i] else r[v]=b.d[i] end
  end
  SPerm(r)
end

Base.:/(a::SPerm,b::SPerm)=a*inv(b)
Base.:^(a::SPerm, b::SPerm)=inv(b)*a*b
Base.:^(a::SPerm, b::Perm)=a^SPerm(b)

@inline function Base.:^(n::Integer, a::SPerm{T}) where T
  if abs(n)>length(a.d) return T(n) end
@inbounds n<0 ? -a.d[-n] : a.d[n]
end

Base.:^(a::SPerm, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                     Base.power_by_squaring(inv(a),-n)

"""
`Base.:^(l::AbstractVector,p::SPerm)`

returns `l` permuted by `p`, a vector `r` such that `r[abs(i^p)]=l[i]*sign(i^p)`.

# Examples
```julia-repl
julia> p=SPerm([-2,-1,-3])
SPerm{Int64}: (1,-2)(3,-3)

julia> [20,30,40]^p
3-element Vector{Int64}:
 -30
 -20
 -40
```
"""
function Base.:^(l::AbstractVector,a::SPerm)
  res=copy(l)
  for i in eachindex(l)
    v=i^a
    if v>0 res[v]=l[i] else res[-v]=-l[i] end
  end
  res
end

function Base.:^(m::AbstractMatrix,a::SPerm;dims=1)
  if dims==1 hcat(map(c->c^a,eachcol(m))...)
  elseif dims==2 permutedims(hcat(map(c->c^a,eachrow(m))...))
  elseif dims==(1,2) hcat(map(c->c^a,eachcol(m))^a...)
  end
end

"""
`SPerm{T}(l::AbstractVector,l1::AbstractVector)`

return  a `SPerm{T}`  `p` such  that `l1^p==l`  if such `p` exists; returns
nothing  otherwise. If not given `{T}` is  taken to be `{Int16}`. Needs the
objects in `l` and `l1` to be sortable.

```julia-repl
julia> p=SPerm([20,30,40],[-40,-20,-30])
(1,-2,3,-1,2,-3)

julia> [20,30,40]^p
3-element Vector{Int64}:
 -40
 -20
 -30
```
"""
function SPerm{T}(a::AbstractVector,b::AbstractVector)where T<:Integer
  pair(x)=x<-x ? (x,-x) : (-x,x)
  p=Perm(pair.(a),pair.(b))
  if isnothing(p) return p end
  res=eachindex(a)^p
  for i in eachindex(a)
    if b[i^(p^-1)]!=a[i] res[i]=-res[i] end
  end
  SPerm{T}(res)
end

SPerm(l::AbstractVector,l1::AbstractVector)=SPerm{Idef}(l,l1)

function SPerm{T}(l::AbstractMatrix,l1::AbstractMatrix;dims=1)where T<:Integer
  if     dims==1 SPerm{T}(collect(eachrow(l)),collect(eachrow(l1)))
  elseif dims==2 SPerm{T}(collect(eachcol(l)),collect(eachcol(l1)))
  end
end

SPerm(l::AbstractMatrix,l1::AbstractMatrix;dims=1)=SPerm{Idef}(l,l1,dims=dims)

"""
`Matrix(a::SPerm)` is the permutation matrix for a
# Examples
```julia-repl
julia> Matrix(SPerm([-2,-1,-3]))
3×3 Matrix{Int64}:
  0  -1   0
 -1   0   0
  0   0  -1
```
"""
function Base.Matrix(a::SPerm,n=length(a.d))
  res=zeros(Int,n,n)
  for (i,v) in enumerate(a.d) res[i,abs(v)]=sign(v) end
  res
end

##underlying Signs of a SignedPerm
#Signs:=p->sign.(p.d)^Perm(p)
#
## We have the properties p=SPerm(Perm(p),Signs(p)) and
## if N=onmats(M,p) then
##    M=onmats(N,Perm(p))^Diagonal(Signs(p)))
#
## Transforms matrix to SPerm
#SignedPerm:=function(m::AbstractMatrix)
#  n=size(m,1)
#  res=[]
#  for i in 1:n
#    nz=filter(x->m[i,x]!=0,1:n)
#    if length(nz)!=1 return false end
#    nz=nz[1]
#    if m[i,nz]==1 res[i]=nz
#    elif m[i,nz]==-1 res[i]=-nz
#    else return false
#    end
#  end
#end
## Transforms perm+signs to a signed perm,
#SPerm(p,n)=SPerm(map(*,eachindex(n),n)^inv(p))
#
#
#------------ Example II: HyperOctaedral groups as Coxeter groups

@GapObj struct CoxHyperoctaedral{T} <: CoxeterGroup{SPerm{T}}
  G::SPermGroup{T}
  n::Int
end

Base.iterate(W::CoxHyperoctaedral,r...)=iterate(W.G,r...)
Groups.gens(W::CoxHyperoctaedral)=gens(W.G)

"""
`CoxHyperoctaedral(n)`  The Hyperoctaedral  group on  ±1,…,±n as  a Coxeter
group  of type  B, with  generators (1,-1)  and (i,i+1)(-i,-i-1); it is the
group of all signed permutations of `1:n`.
```julia-repl
julia> elements(CoxHyperoctaedral(2))
8-element Vector{SPerm{Int8}}:
 ()
 (1,2)
 (1,-1)
 (1,2,-1,-2)
 (1,-2,-1,2)
 (2,-2)
 (1,-2)
 (1,-1)(2,-2)
```
"""
function CoxHyperoctaedral(n::Int)
  gens=[SPerm{Int8}(1,-1)]
  for i in 2:n push!(gens,SPerm{Int8}(i-1,i)) end
  CoxHyperoctaedral{Int8}(Group(gens),n,Dict{Symbol,Any}())
end

function Base.show(io::IO, W::CoxHyperoctaedral)
  print(io,"CoxHyperoctaedral($(W.n))")
end

PermRoot.refltype(W::CoxHyperoctaedral)=[TypeIrred(Dict(:series=>:B,
                                        :indices=>collect(1:W.n)))]

CoxGroups.nref(W::CoxHyperoctaedral)=ngens(W)^2

Gapjm.degrees(W::CoxHyperoctaedral)=2*(1:ngens(W))

Base.length(W::CoxHyperoctaedral)=prod(degrees(W))

function CoxGroups.isleftdescent(W::CoxHyperoctaedral,w,i::Int)
  if i==1 return i^w<0 end
  i^w<(i-1)^w
end

function PermRoot.reflections(W::CoxHyperoctaedral)
  get!(W,:reflections)do
    refs=vcat(gens(W),map(i->SPerm{Int8}(i,-i),2:W.n))
    for i in 2:W.n-1 append!(refs,map(j->SPerm{Int8}(j,j+i),1:W.n-i)) end
    for i in 1:W.n-1 append!(refs,map(j->SPerm{Int8}(j,-j-i),1:W.n-i)) end
    refs
  end
end

PermRoot.reflection(W::CoxHyperoctaedral,i)=reflections(W)[i]

function Perms.reflength(W::CoxHyperoctaedral,w)
  sym=nsym=0
  for x in cycles(w)
    if sort(x)==sort(-x) sym+=length(x)
    else nsym+=length(x)-1
    end
  end
  div(sym,2)+nsym
end

" Only parabolics defined are I=1:m for m≤n"
function PermRoot.reflection_subgroup(W::CoxHyperoctaedral,I::AbstractVector{Int})
  if I!=1:length(I) error(I," should be 1:n for some n") end
  CoxHyperoctaedral(Group(gens(W)[I]),length(I),Dict{Symbol,Any}())
end

#--------------------- action on matrices -----------------------------------
# to find orbits under SPerms transform objects to pairs
pair(x)=x<-x ? (x,-x) : (-x,x)

# duplicate lines and cols of M so group(dup(...)) operates
function dup(M::AbstractMatrix)
  res=zeros(eltype(M),size(M).*2)
  for i in axes(M,1), j in axes(M,2)
    res[[2*i-1,2*i],[2*j-1,2*j]]=[M[i,j] -M[i,j];-M[i,j] M[i,j]]
  end
  res
end

dedup(M::AbstractMatrix)=M[1:2:size(M,1),1:2:size(M,2)]

# transform SPerm on -n:n to Perm acting on  1:2n
dup(p::SPerm)=isone(p) ? Perm() :
    Perm{Idef}(vcat(map(i->i>0 ? [2i-1,2i] : [-2i,-2i-1],p.d)...))

dedup(p::Perm)=SPerm{Idef}(map(i->iseven(i) ? -div(i,2) : div(i+1,2),
                         p.d[1:2:length(p.d)-1]))

dup(g::CoxHyperoctaedral)=Group(dup.(gens(g)))

Base.length(g::SPermGroup)=length(Group(dup.(gens(g))))

function invblocks(m,extra=nothing)
  if isnothing(extra) extra=zeros(Int,size(m,1)) end
  blk1=[collect(axes(m,1))]
  while true
    blk=blk1
    blk1=vcat(map(I->collectby(map(i->
            (tally(pair.(m[i,I])),m[i,i],extra[i]),I),I), blk)...)
    if blk==blk1 return blk end
  end
end

"""
`sstab_onmats([G,]M[,l])`

If `onmats(m,p)=^(M,p;dims=(1,2))` (simultaneous signed conjugation of rows
and  columns, or conjugating by the  matrix of the signed permutation `p`),
and  the argument `G`  is given (which  should be an  `SPermGroup`) this is
just  a fast implementation of  `centralizer(G,M;action=onmats)`. If `G` is
omitted  it is taken to be `CoxHyperoctaedral(size(M,1))`. The program uses
sophisticated  algorithms, and can  handle matrices up  to 80×80. If `l` is
given the return group should also centralize `l` (for the action ^)

```julia-repl
julia> uc=UnipotentCharacters(ComplexReflectionGroup(6));

julia> g=sstab_onmats(fourier(uc.families[2]))
Group([(1,18)(3,-6)(8,-21)(10,-16)(11,22)(13,15),(1,-15)(2,-19)(3,-11)(6,22)(7,-12)(13,-18),(2,19)(4,-14)(5,20)(7,12),(1,-11)(2,-19)(3,-15)(5,-20)(6,13)(8,10)(16,21)(17,-17)(18,-22),(1,-22)(2,-19)(3,-13)(5,-20)(6,15)(8,-16)(10,-21)(11,-18)(17,-17),(1,6)(2,-19)(3,-18)(4,14)(8,16)(9,-9)(10,21)(11,-13)(15,-22),(1,13)(3,22)(4,14)(5,-20)(6,-11)(8,21)(9,-9)(10,16)(15,18)(17,-17)])
julia> length(g)
32
```
"""
function sstab_onmats(M,extra=nothing)
  k=size(M,1)
  if M!=permutedims(M) error("M should be symmetric") end
  if isnothing(extra) extra=fill(1,size(M,1)) end
  blocks=sort(invblocks(M),by=length)
  gen=SPerm{Idef}[]
  I=Int[]
  for r in blocks
    if length(r)>5 InfoChevie("#IS Large Block:",r,"\n") end
    gr=stab_onmats(dup(CoxHyperoctaedral(length(r))),dup(M[r,r]))
    p=SPerm(mappingPerm(1:length(r),r).d)
    append!(gen,map(x->dedup(x)^p,gens(gr)))
    append!(I,r)
    p=SPerm(mappingPerm(I,eachindex(I)).d)
    gen=gen.^p
    gen=dedup.(gens(stab_onmats(Group(dup.(gen)),dup(M[I,I]))))
    gen=gen.^inv(p)
  end
  return Group(gen)
end

"""
`SPerm_onmats(M,N[,m,n])`

`M`  and `N` should be symmetric  matrices. `SPerm_onmats` returns a signed
permutation `p` such that `onmats(M,p)=N` if such a permutation exists, and
`nothing`  otherwise. If  in addition  vectors `m`  and `n`  are given, the
signed permutation `p` should also satisfy `m^p==n`.

This  routine is  useful to  identify two  objects which are isomorphic but
with  different  labelings.  It  is  used  in   Chevie  to identify Lusztig
Fourier  transform matrices  with standard  (classified) data.  The program
uses  sophisticated  algorithms,  and  can  often  handle  matrices  up  to
80×80.

Efficient version of
`transporting_elt(CoxHyperoctaedral(size(M,1)),M,N;action=onmats)`

```julia-repl
julia> f=SubFamilyij(chevieget(:families,:X)(12),1,3,(3+root(-3))/2);

julia> M=fourier(conj(f));

julia> uc=UnipotentCharacters(ComplexReflectionGroup(6));

julia> N=fourier(uc.families[2]);

julia> p=SPerm_onmats(M,N)
(1,3)(2,19,-2,-19)(4,-14,-4,14)(5,-5)(6,-18)(7,-7)(8,10)(11,15,-11,-15)(12,-12)(13,22)(16,21,-16,-21)

julia> ^(M,p;dims=(1,2))==N
true
```
"""
function SPerm_onmats(M,N,extra1=nothing,extra2=nothing)
  onmats(M,p)=^(M,p;dims=(1,2))
  if M!=permutedims(M) error("M should be symmetric") end
  if N!=permutedims(N) error("N should be symmetric") end
  if isnothing(extra1)
    extra1=fill(1,size(M,1))
    extra2=fill(1,size(M,1))
  end
  function ind(I,J)local iM,iN,p,n;
    invM=map(i->(tally(pair.(M[i,I])),M[i,i],extra1[i]),I)
    invN=map(i->(tally(pair.(M[i,J])),M[i,i],extra2[i]),J)
    if tally(invM)!=tally(invN) InfoChevie("content differs");return false end
    iM=collectby(invM,I)
    iN=collectby(invN,J)
    if length(iM)==1
      if length(I)>6 InfoChevie("large block:",length(I),"\n")
        p=transporting_elt(Group(reflections(
          CoxHyperoctaedral(length(I)))), M[I,I],N[J,J];action=onmats,
                    dist=(M,N)->count(i->M[i]!=N[i],eachindex(M)))
      else p=transporting_elt(CoxHyperoctaedral(length(I)),
         M[I,I],N[J,J],action=onmats)
      end
      if isnothing(p) InfoChevie("could not match block");return nothing end
      return [[I,J,p]]
    else p=map(ind,iM,iN)
      if false in p return false else return vcat(p...) end
    end
  end
  l=ind(axes(M,1),axes(N,1))
  if l==false return false end
  I=Int[];J=Int[];g=SPermGroup();tr=SPerm()
  for r in l
#   Print("r=",r,"\n");
    n=length(r[1])
#   q=mappingPerm(eachindex(I),eachindex(I))
    q=SPerm()
    p=SPerm(mappingPerm(1:n,(1:n)+length(I)).d)
    append!(I,r[1]);append!(J,r[2])
#   Print("#I=",Length(I),"\c");
    if comm(r[3]^p,tr^q)!=SPerm() error("noncomm") end
    tr=tr^q*r[3]^p
    h=gens(sstab_onmats(M[r[1],r[1]])).^p
    g=Group(vcat(gens(g).^q,h))
#   Print(" #g=",Size(g),"\c");
    e=transporting_elt(g,M[I,I],onmats(N[J,J],tr^-1),action=onmats)
    if e==false return false
    elseif e^-1*e^tr!=SPerm() print("*** tr does not commute to e\n")
         tr=e*tr
    end
    g=stab_onmats(Group(dup.(gens(g))),dup(M[I,I]))
#   Print(" #stab=",Size(g),"\n");
    g=Group(dedup.(gens(g)))
  end
  # transporter of a ps from [1..Length(I)] to I
  trans=I->SPerm(mappingPerm(eachindex(I),I).d)
  trans(I)^-1*tr*trans(J)
end
end
