"""
This module is a port of the GAP permutations type.

A  permutation here is a permutation of the set `1:n` and is represented as
a  vector of `n` integers representing the images of `1:n`. The integer `n`
is called the *degree* of the permutation, even if it is not moved.

The  permutations in this module  follow the GAP design:  it is possible to
multiply, or to store in the same group, permutations of different degrees.
Slightly  faster is the MAGMA design where any permutation has to belong to
a  group and the degree is determined by that group. Then multiplication of
permutations  in a given group is slightly faster, but it is more difficult
to multiply permutations coming from different groups, like a group and one
of its subgroups.

The  permutation whose cycle decomposition is `(1,2,3)(4,5)` can be written
`Perm(1,2,3)*Perm(4,5)`   or   `perm"(1,2,3)(4,5)"`.   It   is  represented
internally  as `[2,3,1,5,4]`; note that `[2,3,1,5,4,6]` represents the same
permutation.

As  in  GAP  `i^p`  applies  `p`  to  the  integer  `i`,  while `p^q` means
`p^-1*q*p`.

The  complete  type  of  our  permutations is `Perm{T}` where `T<:Integer`,
where `Vector{T}` is the type of the vector which holds the image of `1:n`.
This   can  used  to  save  space  or  time  when  possible.  For  instance
`Perm{UInt8}`  uses less  space than  `Perm{Int}` and  can be used for Weyl
groups of rank <=8 since they have at most 240 roots.

# Examples
```julia-repl
julia> p=Perm(1,2)*Perm(2,3)
(1,3,2)

julia> Perm{Int8}(p)
Int8(1,3,2)

julia> 1^p
3

julia> Matrix(p)
3×3 Array{Int64,2}:
 0  0  1
 1  0  0
 0  1  0

julia> p^Perm(3,10)
(1,10,2)

julia> inv(p)
(1,2,3)

julia> one(p)
()

julia> order(p)
3

julia> degree.((Perm(1,2),Perm(2,3)))
(2, 3)

julia> largest_moved_point(Perm(1,2)*Perm(2,3)^2)
2

julia> smallest_moved_point(Perm(2,3))
2
```

```not-in-tests
julia> rand(Perm,10)
(1,8,4,2,9,7,5,10,3,6)

```

Operations on permutations are `*, /, inv, \` and `mul!`
Perms  have methods `copy, hash,  ==, cmp, isless` (total order)  so they can be
keys in hashes or elements of sets.

other functions are: `cycles, cycletype, sign, rand`. 
See individual documentations.

GAP→ Julia dictionary
```
     PermList(v)                      →  Perm(v) 
     Permuted(v,p)                    →  v[p.d]
     ListPerm(p)                      →  p.d
     PermListList(l1,l2)              →  Perm(l1,l2)
     OnTuples(l,p)                    →  l.^p or (faster) p.d[l]
     RestrictedPerm                   →  restricted
```
"""
module Perms

using Gapjm 

export Perm, largest_moved_point, cycles, cycletype, order, sign,
  @perm_str, smallest_moved_point

struct Perm{T<:Integer}
   d::Vector{T}
end

# Gap's Permuted(a,p) is a[p.d], Gap's ListPerm(p) is p.d
#---------------- Constructors ---------------------------------------
"for example  Perm{Int8}(1,2,3) constructs a cycle"
function Perm{T}(x::Int...)where T<:Integer
  if isempty(x) return Perm(T[]) end
  d=T.(1:max(x...))
  for i in 1:length(x)-1
    d[x[i]]=x[i+1]
  end
  d[x[end]]=x[1]
  Perm(d)
end

"Perm(1,2,3)==Perm{Int}(1,2,3)"
Perm(x::Int...)=Perm{Int}(x...)

"Perm{Int8}(Perm(1,2,3))==Perm{Int8}(1,2,3)"
Perm{T}(p::Perm) where T<:Integer=Perm(T.(p.d))

"""
allows GAP-style perm"(1,2)(5,6,7)(4,9)"
"""
macro perm_str(s::String)
  start=1
  res=Perm()
  while true
    m=match(r"\((\s*\d+\s*,)+\s*\d+\)",s[start:end])
    if isnothing(m) break end
    start+=m.match.ncodeunits
    res*=Perm(Meta.parse(m.match).args...)
  end
  res::Perm
end

"""
just for fun: Perm[1 2;5 6 7;4 9]=perm"(1,2)(5,6,7)(4,9)"
"""
function Base.typed_hvcat(::Type{Perm},a::Tuple{Vararg{Int64,N} where N},
  b::Vararg{Number,N} where N)
  res=Perm()
  for i in a
    res*=Perm(Iterators.take(b,i)...)
    b=Iterators.drop(b,i)
  end
  res
end

" find permutation mapping l to l1 if exists (like GAP's PermListList)"
function Perm{T}(l::AbstractVector,l1::AbstractVector)where T<:Integer
  s=sortperm(l)
  s1=sortperm(l1)
  if !all(i->l[s[i]]==l1[s1[i]],eachindex(l)) error("not permuted") end
  Perm{T}(s1)\Perm{T}(s)
end

Perm(l::AbstractVector,l1::AbstractVector)=Perm{Int}(l,l1)
#---------------------------------------------------------------------
Base.one(p::Perm)=Perm(empty(p.d))
Base.one(::Type{Perm{T}}) where T=Perm(T[])
Base.copy(p::Perm)=Perm(copy(p.d))

import ..Gapjm.degree
@inline degree(a::Perm)=length(a.d)
#Base.vec(a::Perm)=a.d

" hash is needed for using Perms in Sets/Dicts"
function Base.hash(a::Perm, h::UInt)
  b = 0x595dee0e71d271d0%UInt
  for (i,v) in enumerate(a.d)
    if v!=i
      b = xor(b,xor(hash(v, h),h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
    end
  end
  b
end

" permutations are scalars for broadcasting"
Base.broadcastable(p::Perm)=Ref(p)

" total order is needed to use Perms in sorted lists"
function Base.cmp(a::Perm, b::Perm)
  da=length(a.d)
  db=length(b.d)
  for i in 1:min(da,db)
@inbounds if a.d[i]<b.d[i] return -1 end
@inbounds if a.d[i]>b.d[i] return  1 end
  end
  if     da<db for i in (da+1:db) b.d[i]==i || return -1 end
  elseif da>db for i in (db+1:da) a.d[i]==i || return  1 end
  end
  0
end

Base.isless(a::Perm, b::Perm)=cmp(a,b)==-1

Base.:(==)(a::Perm, b::Perm)= cmp(a,b)==0

Base.rand(::Type{Perm},i::Integer)=Perm(sortperm(rand(1:i,i)))
Base.rand(::Type{Perm{T}},i::Integer) where T=Perm(T.(sortperm(rand(1:i,i))))

"Matrix(a::Perm) is the permutation matrix for a"
Base.Matrix(a::Perm,n=length(a.d))=Int[j==i^a for i in 1:n, j in 1:n]

" largest_moved_point(a::Perm) is the largest integer moved by a"
largest_moved_point(a::Perm)=findlast(x->a.d[x]!=x,eachindex(a.d))

" smallest_moved_point(a::Perm) is the smallest integer moved by a"
smallest_moved_point(a::Perm)=findfirst(x->a.d[x]!=x,eachindex(a.d))

#------------------ operations on permutations --------------------------

" `promote(a::Perm, b::Perm)` promotes `a` and `b` to the same degree"
function Base.promote(a::Perm,b::Perm)
  da=length(a.d)
  db=length(b.d)
  if da<db
    resize!(a.d,db)
@inbounds    a.d[da+1:db]=da+1:db
  elseif db<da
    resize!(b.d,da)
@inbounds    b.d[db+1:da]=db+1:da
  end
  (a,b)
end

function Base.:*(a::Perm, b::Perm)
  a,b=promote(a,b)
  r=similar(a.d)
@inbounds for (i,v) in enumerate(a.d) r[i]=b.d[v] end
  Perm(r)
end

# this is a*=b without allocation
function mul!(a::Perm, b::Perm)
  a,b=promote(a,b)
@inbounds for (i,v) in enumerate(a.d) a.d[i]=b.d[v] end
  a
end

function Base.inv(a::Perm)
  r=similar(a.d)
@inbounds for (i,v) in enumerate(a.d) r[v]=i end
  Perm(r)
end

# I do not know how to do this one faster
Base.:/(a::Perm, b::Perm)=a*inv(b)

# less allocations than inv(a)*b
function Base.:\(a::Perm, b::Perm)
  a,b=promote(a,b)
  r=similar(a.d)
@inbounds for (i,v) in enumerate(a.d) r[v]=b.d[i] end
  Perm(r)
end

# less allocations than inv(a)*b*a
function Base.:^(a::Perm, b::Perm)
  a,b=promote(a,b)
  r=similar(a.d)
@inbounds for (i,v) in enumerate(a.d) r[b.d[i]]=b.d[v] end
  Perm(r)
end

@inline Base.:^(n::Integer, a::Perm{T}) where T=
   @inbounds if n>length(a.d) T(n) else a.d[n] end

Base.:^(a::Perm, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                               Base.power_by_squaring(inv(a),-n)

#---------------------- cycle decomposition -------------------------
function cycle(a::Perm{T},i::Integer,check=false)where T
  res=T[]
  j=i
  while true
    if check && j in res error("point $j occurs twice") end
    push!(res,j)
    if (j=a.d[j])==i return res end
  end
end
  
"""
  cycles(a::Perm) returns the cycles of a
# Example
```julia-repl
julia> cycles(Perm(1,2)*Perm(4,5))
3-element Array{Array{Int64,1},1}:
 [1, 2]
 [3]
 [4, 5]
```
"""
function cycles(a::Perm{T};domain=1:length(a.d),check=false)where T
  to_visit=falses(length(a.d))
  to_visit[domain].=true
  cycles=Vector{T}[]
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    cyc=cycle(a,i,check)
    to_visit[cyc].=false
    push!(cycles,cyc)
  end
  cycles
end

function Base.show(io::IO, a::Perm{T}) where T
  if T!=Int print(io,T) end
  cyc=filter(c->length(c)>1,cycles(a,check=true))
  if isempty(cyc) print(io,"()")
  else for c in cyc print(io,"(",join(c,","),")") end
  end
end

order(a::Perm) = lcm(length.(cycles(a)))

"""
  cycletype(a::Perm) describes the partition of degree(a) associated to the
  conjugacy class of a in the symmetric group, with ones removed. It is
  represented as a Dict of cyclesize=>multiplicity
# Example
```julia-repl
julia> cycletype(Perm(1,2)*Perm(3,4))
1-element Array{Pair{Tuple{Int64,Int64},Int64},1}:
 (2, 1) => 2
```
"""
function cycletype(a::Perm;domain=ones(Int16,length(a.d)))
  res=Dict{Tuple{Int,Int},Int}()
  to_visit=domain[1:length(a.d)] # this makes a copy
  for i in eachindex(to_visit)
    if iszero(to_visit[i]) continue end
    l=0
    j=i
    color=to_visit[i]
    while true
      to_visit[j]=0
      l+=1
      if (j=a.d[j])==i break end
    end
    if l>1 
      if haskey(res,(l,color)) res[(l,color)]+=1
      else res[(l,color)]=1
      end
    end
  end
  sort(collect(res),by=x->x[1])
end

function nrcycles(a::Perm)
  to_visit=trues(length(a.d))
  nr=0
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    nr+=1
    j=i
    while true
      to_visit[j]=false
      j=a.d[j]
      if j==i break end
    end
  end
  nr
end

" sign(a::Perm) is the signature of  the permutation a"
Base.sign(a::Perm)=(-1)^(length(a.d)-nrcycles(a)) # nr of even cycles

# l should be a union of cycles of p
# returns p restricted to l
function Gapjm.restricted(a::Perm{T},l::AbstractVector{<:Integer})where T
  res=one(a)
  while !isempty(l)
    c=cycle(a,l[1])
    l=setdiff(l,c)
    res*=Perm{T}(Int.(c)...)
  end
  res
end

end
