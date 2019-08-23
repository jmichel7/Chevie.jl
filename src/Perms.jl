"""
This module is a port of the GAP permutations type.

A  permutation here is a permutation of the set `1:n` and is represented as
a  vector of `n` integers representing the images of `1:n`. The integer `n`
is called the degree of the permutation, even if it is not moved. We follow
the  GAP design: it is possible to multiply, or to store in the same group,
permutations of different degrees; this is implemented by promoting both to
the   higher  degree.  Slightly  faster  is  the  MAGMA  design  where  any
permutation  has to belong to a group  and the degree is determined by that
group. Then multiplication of permutations within a given group is slightly
faster,  but  it  is  more  difficult  to multiply permutations coming from
different  groups, like a group and one  of its subgroups. The degree is an
implementation  detail so usually  it should not  be used. One sould rather
use the function `largest_moved_point`.

A  permutation  can  be  defined  by  the  list  of  images  of `1:n`, like
`Perm([2,3,1,5,4])`.  Usually it is rather  given by its cycle decomposion:
the  permutation whose cycle decomposition is `(1,2,3)(4,5)` can be written
`Perm(1,2,3)*Perm(4,5)`  or  `perm"(1,2,3)(4,5)"`.  The  list  of images of
`1:n`  is gotten back from the permutation by the function `vec`; note that
since  equal  permutations  may  have  different  degrees,  they  may  have
different `vec`.

The  complete  type  of  our  permutations is `Perm{T}` where `T<:Integer`,
where `Vector{T}` is the type of the vector which holds the image of `1:n`.
This   can  used  to  save  space  or  time  when  possible.  For  instance
`Perm{UInt8}`  uses less  space than  `Perm{Int}` and  can be used for Weyl
groups of rank <=8 since they have at most 240 roots.

# Examples
```julia-repl
julia> a=Perm(1,2,3)
(1,2,3)

julia> vec(a)
3-element Array{Int64,1}:
 2
 3
 1

julia> a==Perm(vec(a))
true

julia> b=Perm(1,2,3,4)
(1,2,3,4)

julia> a*b     # product
(1,3,2,4)

julia> inv(a)  # inverse
(1,3,2)

julia> a/b     # quotient  a*inv(b)
(3,4)

julia> a\\b     # left quotient inv(a)*b
(1,4)

julia> a^b     # conjugation inv(b)*a*b
(2,3,4)

julia> b^2
(1,3)(2,4)

julia> 1^a     # apply a to point 1
2

julia> one(a)
()

julia> sign(a)
1

julia> order(a)
3

julia> largest_moved_point(a)
3

julia> smallest_moved_point(a)
1

julia> Perm{Int8}(a) # convert to Perm{Int8}
Int8(1,2,3)

julia> Matrix(b)
4×4 Array{Int64,2}:
 0  1  0  0
 0  0  1  0
 0  0  0  1
 1  0  0  0
```

```not-in-tests
julia> rand(Perm,10)
(1,8,4,2,9,7,5,10,3,6)

```

Perms  have methods `copy, hash, ==, cmp, isless` (total order) so they can
be  keys in hashes or elements of  sets; two permutations are equal if they
move   the  same  points.  Permutations   are  considered  as  scalars  for
broadcasting.

other functions are: 
`cycles, cycletype, orbit, orbits, permuted, rand, restricted, sign`. 
See individual documentations.

GAP→ Julia dictionary
```
     PermList(v)                      →  Perm(v) 
     Permuted(v,p)                    →  permuted(v,p)
     ListPerm(p)                      →  vec(p)
     PermListList(l1,l2)              →  Perm(l1,l2)
     OnTuples(l,p)                    →  l.^p
     RestrictedPerm(p,d)              →  restricted(p,d)
```
"""
module Perms

using Gapjm, ..Groups
import ..Gapjm.degree
import ..Gapjm.restricted

export Perm, largest_moved_point, cycles, cycletype, order, sign,
  @perm_str, smallest_moved_point, reflength, permuted

struct Perm{T<:Integer}
   d::Vector{T}
end

Base.vec(a::Perm)=a.d

#---------------- Constructors ---------------------------------------
"for example  Perm{Int8}(1,2,3) constructs the cycle (1,2,3)"
function Perm{T}(x::Integer...)where T<:Integer
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
allow GAP-style perm"(1,2)(5,6,7)(4,9)"
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

"""
  Perm{T}(l::AbstractVector,l1::AbstractVector)
  return permutation p such that Permuted(l1,p)==l if such p exists
"""
function Perm{T}(l::AbstractVector,l1::AbstractVector)where T<:Integer
  s=sortperm(l)
  s1=sortperm(l1)
  if !all(i->l[s[i]]==l1[s1[i]],eachindex(l)) error("not permuted") end
  Perm{T}(s1)\Perm{T}(s)
end

Perm(l::AbstractVector,l1::AbstractVector)=Perm{Int}(l,l1)

"""
assume l is union of orbits under group elt g; return permutation of l by g
needs objects in l sortable
"""
Perm{T}(g,l::AbstractVector;action::Function=^) where T<:Integer=Perm{T}(l,action.(l,Ref(g)))

Perm(g,l::AbstractVector;action::Function=^)=Perm{Int}(g,l,action=action)

#---------------------------------------------------------------------
Base.one(p::Perm)=Perm(empty(p.d))
Base.one(::Type{Perm{T}}) where T=Perm(T[])
Base.copy(p::Perm)=Perm(copy(p.d))

@inline degree(a::Perm)=length(a.d)

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

"`Matrix(a::Perm)` is the permutation matrix for a"
Base.Matrix(a::Perm,n=length(a.d))=Int[j==i^a for i in 1:n, j in 1:n]

" `largest_moved_point(a::Perm)` is the largest integer moved by a"
largest_moved_point(a::Perm)=findlast(x->a.d[x]!=x,eachindex(a.d))

" `smallest_moved_point(a::Perm)` is the smallest integer moved by a"
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

#---------------------- cycles -------------------------
"""
  orbit(a::Perm,i::Integer) returns the orbit of a on i
"""
function Groups.orbit(a::Perm{T},i::Integer,check=false)where T
  res=T[]
  j=i
  while true
    if check && j in res error("point $j occurs twice") end
    push!(res,j)
    if (j=a.d[j])==i return res end
  end
end
  
"""
orbits(a::Perm,d::Vector=1:length(vec(a))) returns the orbits of a on domain d
# Example
```julia-repl
julia> orbits(Perm(1,2)*Perm(4,5),1:5)
3-element Array{Array{Int64,1},1}:
 [1, 2]
 [3]
 [4, 5]
```
"""
function Groups.orbits(a::Perm,domain=1:length(a.d);trivial=true,check=false)
  to_visit=falses(length(a.d))
  to_visit[domain].=true
  cycles=Vector{eltype(a.d)}[]
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    cyc=orbit(a,i,check)
    to_visit[cyc].=false
    if length(cyc)>1 || trivial push!(cycles,cyc) end
  end
  cycles
end

"""
  cycles(a::Perm) returns the non-trivial cycles of a
# Example
```julia-repl
julia> cycles(Perm(1,2)*Perm(4,5))
3-element Array{Array{Int64,1},1}:
 [1, 2]
 [4, 5]
```
"""
cycles(a::Perm;check=false)=orbits(a;trivial=false)

function Base.show(io::IO, a::Perm{T}) where T
  if T!=Int print(io,T) end
  cyc=orbits(a;trivial=false,check=true)
  if isempty(cyc) print(io,"()")
  else for c in cyc print(io,"(",join(c,","),")") end
  end
end

" order(a) is the order of the permutation a"
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

" reflength(a::Perm) minimum number of transpositions of which a is product"
function reflength(a::Perm)
  to_visit=trues(length(a.d))
  l=0
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    j=i
    while true
      to_visit[j]=false
      if (j=a.d[j])==i break end
      l+=1
    end
  end
  l
end

" sign(a::Perm) is the signature of  the permutation a"
Base.sign(a::Perm)=(-1)^reflength(a)
#---------------------- other -------------------------

"""
   permuted(l,a) returns l permuted by a as a new list r,
   that is r[i^a]==l[i]
"""
function permuted(l::AbstractVector,a::Perm)
  res=copy(l)
  res[eachindex(l).^a].=l
  res
end

"""
   restricted(a::Perm{T},l::AbstractVector{<:Integer})

l should be a union of cycles of p; returns p restricted to l
"""
function restricted(a::Perm{T},l::AbstractVector{<:Integer})where T
  o=orbits(a,l;trivial=false)
  isempty(o) ? Perm{T}() : prod(c->Perm{T}(c...),o)
end

end
