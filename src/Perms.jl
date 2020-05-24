"""
This module is a port of the GAP permutations type.

They  are permutations  of the  set `1:n`  represented as  a vector  of `n`
integers  holding the images of `1:n`. The integer `n` is called the degree
of  the permutation, even if it is not  moved. We follow the GAP design: it
is  possible to multiply,  or to store  in the same  group, permutations of
different  degrees; this  is implemented  by promoting  both to  the higher
degree. Slightly different is the MAGMA design where any permutation has to
belong  to  a  group  and  the  degree  is  determined  by that group; then
multiplication of permutations within a given group is slightly faster, but
it is more difficult to multiply permutations coming from different groups,
like  a group  and one  of its  subgroups. The  degree is an implementation
detail so usually it should not be used. One should rather use the function
`largest_moved_point`.

The default constructor for a permutation uses the list of images of `1:n`,
like  `Perm([2,3,1,5,4])`.  Often  it  is  more  convenient  to  use  cycle
decompositions:    the   above   permutation    has   cycle   decomposition
`(1,2,3)(4,5)`    thus   can   be    written   `Perm(1,2,3)*Perm(4,5)`   or
`perm"(1,2,3)(4,5)"`  (this last form  can parse any  GAP permutation). The
list  of images  of `1:n`  can be  gotten back  from the permutation by the
function  `vec`;  note  that  permutations  may  be equal even if they have
different  degrees  (if  they  move  the  same points), then they will have
different `vec`.

The  complete type of a permutation  is `Perm{T}` where `T<:Integer`, where
`Vector{T}`  is the type of the vector which holds the image of `1:n`. This
can  be used to save space or  time. For instance `Perm{UInt8}` can be used
for  Weyl groups of rank≤8 since they permute  at most 240 roots. If `T` is
not  specified we  take it  to be  `Int16` since  this is a good compromise
between   speed,  compactness  and  possible  size  of  `n`.  One  can  mix
permutations  of different  types `T`;  they are  promoted to the wider one
when multiplying.

# Examples of operations with permutations
```julia-repl
julia> a=Perm(1,2,3)
(1,2,3)

julia> vec(a)
3-element Array{Int16,1}:
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

julia> b^2     # square
(1,3)(2,4)

julia> 1^a     # image by a of point 1
2

julia> one(a)
()

julia> sign(a) # sigature of permutation
1

julia> order(a)
3

julia> largest_moved_point(a)
3

julia> smallest_moved_point(a)
1

julia> Perm{Int8}(a) # convert to Perm{Int8}
Perm{Int8}: (1,2,3)

julia> Matrix(b)
4×4 Array{Bool,2}:
 0  1  0  0
 0  0  1  0
 0  0  0  1
 1  0  0  0
```
```julia-rep1
julia> rand(Perm,10)
(1,8,4,2,9,7,5,10,3,6)
```

`Perm`s have methods `copy`, `hash`, `==`, so they can be keys in hashes or
elements  of sets; two permutations are equal  if they move the same points
to  the same images. They have methods `cmp`, `isless` (lexicographic order
on   moved  points)  so  they  can  be  sorted.  `Perm`s  are  scalars  for
broadcasting.
"""
module Perms

#import Gapjm: degree, restricted, order
#import ..Groups: orbit, orbits
# to use as a stand-alone module comment above 2 lines and uncomment next
export degree, restricted, orbit, orbits, order

export Perm, largest_moved_point, cycles, cycletype, support,
  @perm_str, smallest_moved_point, reflength, mappingPerm

"""
`struct Perm{T<:Integer}`

A  Perm represents a permutation  of the set `1:n`  and is implemented by a
`struct` with one field, a `Vector{T}` holding the images of `1:n`.
"""
struct Perm{T<:Integer}
   d::Vector{T}
end

Base.vec(a::Perm)=a.d

#---------------- Constructors ---------------------------------------
"""
   `Perm{T}(x::Integer...)where T<:Integer`

   returns  a cycle.  For example  `Perm{Int8}(1,2,3)` constructs the cycle
   `(1,2,3)` as a `Perm{Int8}`. If omitted `{T}` is taken as to be `Int16`.
"""
function Perm{T}(x::Integer...)where T<:Integer
  if isempty(x) return Perm(T[]) end
  d=T.(1:max(x...))
  for i in 1:length(x)-1
    d[x[i]]=x[i+1]
  end
  d[x[end]]=x[1]
  Perm(d)
end

Perm(x::Integer...)=Perm{Int16}(x...)

"""
   `Perm{T}(p::Perm) where T<:Integer`

   change the type of `p` to `Perm{T}`
   for example `Perm{Int8}(Perm(1,2,3))==Perm{Int8}(1,2,3)`
"""
Perm{T}(p::Perm) where T<:Integer=convert(Perm{T},p)

"""
   @perm"..."

 make a `Perm` from a string; allows GAP-style `perm"(1,2)(5,6,7)(4,9)"`
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
function Base.typed_hvcat(::Type{Perm},a::Tuple{Vararg{Int,N} where N},
  b::Vararg{Number,N} where N)
  res=Perm()
  for i in a
    res*=Perm(Iterators.take(b,i)...)
    b=Iterators.drop(b,i)
  end
  res
end

"""
  `Perm{T}(l::AbstractVector,l1::AbstractVector)`

returns `p`, a `Perm{T}`, such that `l1^p==l` if such a `p` exists; returns
`nothing` otherwise. If not given `{T}` is taken to be `{Int16}`. Needs the
elements of `l` and `l1` to be sortable.

```julia-repl
julia> Perm([0,2,4],[4,0,2])
(1,3,2)
```
"""
function Perm{T}(l::AbstractVector,l1::AbstractVector)where T<:Integer
  p=sortperm(l)
  p1=sortperm(l1)
  if l[p]==l1[p1] Perm{T}(p1)\Perm{T}(p) end
end

Perm(l::AbstractVector,l1::AbstractVector)=Perm{Int16}(l,l1)

function Perm{T}(l::AbstractMatrix,l1::AbstractMatrix;dims=1)where T<:Integer
  if     dims==1 Perm{T}(collect(eachrow(l)),collect(eachrow(l1)))
  elseif dims==2 Perm{T}(collect(eachcol(l)),collect(eachcol(l1)))
  end
end

Perm(l::AbstractMatrix,l1::AbstractMatrix;dims=1)=Perm{Int16}(l,l1,dims=dims)

#---------------------------------------------------------------------
Base.one(p::Perm)=Perm(empty(p.d))
Base.one(::Type{Perm{T}}) where T=Perm(T[])
Base.copy(p::Perm)=Perm(copy(p.d))

@inline degree(a::Perm)=length(a.d)

# hash is needed for using Perms in Sets/Dicts
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

# Perms are scalars for broadcasting"
Base.broadcastable(p::Perm)=Ref(p)

# total order is needed to use Perms in sorted lists
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

Base.rand(::Type{Perm},i::Integer)=Perm(Int16.(sortperm(rand(1:i,i))))
Base.rand(::Type{Perm{T}},i::Integer) where T=Perm(T.(sortperm(rand(1:i,i))))

"""
`Matrix(a::Perm,n=degree(a))` 
the  permutation  matrix  for  `a`  operating  on `n` points. If given, `n`
should be larger than `largest_moved_point(a)`.

```julia-repl
julia> Matrix(Perm(2,3,4),5)
5×5 Array{Bool,2}:
 1  0  0  0  0
 0  0  1  0  0
 0  0  0  1  0
 0  1  0  0  0
 0  0  0  0  1
```
"""
Base.Matrix(a::Perm,n=length(a.d))=[j==i^a for i in 1:n, j in 1:n]

" `largest_moved_point(a::Perm)` is the largest integer moved by a"
largest_moved_point(a::Perm)=findlast(x->a.d[x]!=x,eachindex(a.d))

" `smallest_moved_point(a::Perm)` is the smallest integer moved by a"
smallest_moved_point(a::Perm)=findfirst(x->a.d[x]!=x,eachindex(a.d))

support(a::Perm)=findall(x->a.d[x]!=x,eachindex(a.d))
#------------------ operations on permutations --------------------------

Base.convert(::Type{Perm{T}},p::Perm{T1}) where {T,T1}=T==T1 ? p : Perm(T.(p.d))

function Base.promote_rule(a::Type{Perm{T1}},b::Type{Perm{T2}})where {T1,T2}
  Perm{promote_type(T1,T2)}
end

extend!(a::Perm,n::Integer)=if length(a.d)<n append!(a.d,length(a.d)+1:n) end

# `promote_length(a::Perm, b::Perm)` extends `a` and `b` to the same degree"
function promote_length(a::Perm,b::Perm)
  a,b=promote(a,b)
  extend!(a,length(b.d))
  extend!(b,length(a.d))
  (a,b)
end

function Base.:*(a::Perm, b::Perm)
  a,b=promote_length(a,b)
  r=similar(a.d)
@inbounds for (i,v) in enumerate(a.d) r[i]=b.d[v] end
  Perm(r)
end

# this is a*=b without allocation
function mul!(a::Perm, b::Perm)
  a,b=promote_length(a,b)
@inbounds for (i,v) in enumerate(a.d) a.d[i]=b.d[v] end
  a
end

function Base.inv(a::Perm)
  r=similar(a.d)
@inbounds for (i,v) in enumerate(a.d) r[v]=i end
  Perm(r)
end

# less allocations than inv(a)*b
function Base.:\(a::Perm, b::Perm)
  a,b=promote_length(a,b)
  r=similar(a.d)
@inbounds for (i,v) in enumerate(a.d) r[v]=b.d[i] end
  Perm(r)
end

# less allocations than inv(a)*b*a
function Base.:^(a::Perm, b::Perm)
  a,b=promote_length(a,b)
  r=similar(a.d)
@inbounds for (i,v) in enumerate(a.d) r[b.d[i]]=b.d[v] end
  Perm(r)
end

# I do not know how to do this one faster
Base.:/(a::Perm, b::Perm)=a*inv(b)

@inline Base.:^(n::Integer, a::Perm{T}) where T=
  if n>length(a.d) T(n) 
  else
@inbounds a.d[n] 
  end

Base.:^(a::Perm, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                               Base.power_by_squaring(inv(a),-n)

"""
`Base.:^(l::AbstractVector,p::Perm)` 

   returns `l` permuted by `p`, a vector `r` such that `r[i^p]==l[i]`

# Examples
```julia-repl
julia> [5,4,6,1,7,5]^Perm(1,3,5,6,4)
6-element Array{Int64,1}:
 1
 4
 5
 5
 6
 7
```
"""
function Base.:^(l::AbstractVector,a::Perm)
  res=collect(l)
  res[eachindex(l).^a].=l
  res
end

"""
Base.:^(m::AbstractMatrix,p::Perm;dims=1)

Applies the permutation `p` on the lines, columns or both of the matrix `m`
depending on the value of `dims`

```julia-repl
julia> m=[3*i+j for i in 0:2,j in 1:3]
3×3 Array{Int64,2}:
 1  2  3
 4  5  6
 7  8  9

julia> p=Perm(1,2,3)
(1,2,3)

julia> m^p
3×3 Array{Int64,2}:
 7  8  9
 1  2  3
 4  5  6

julia> ^(m,p;dims=2)
3×3 Array{Int64,2}:
 3  1  2
 6  4  5
 9  7  8

julia> ^(m,p;dims=(1,2))
3×3 Array{Int64,2}:
 9  7  8
 3  1  2
 6  4  5
```
"""
function Base.:^(m::AbstractMatrix,a::Perm;dims=1)
  if dims==2 m[:,axes(m,2)^a]
  elseif dims==1 m[axes(m,1)^a,:]
  elseif dims==(1,2) m[axes(m,1)^a,axes(m,2)^a]
  end
end

#---------------------- cycles -------------------------

# takes 20% more time than GAP CyclePermInt for rand(Perm,1000)
"""
  orbit(a::Perm,i::Integer) returns the orbit of a on i
"""
function orbit(a::Perm{T},i::Integer,check=false)where T
  if i>length(a.d) return T[i] end
  res=T[]
  sizehint!(res,length(a.d))
  j=i
  while true
    if check && j in res error("point $j occurs twice") end
    push!(res,j)
@inbounds j=a.d[j]
    if j==i return res end
  end
end
  
"""
`orbits(a::Perm,d::Vector=1:degree(a))` 

returns the orbits of a on domain d

# Example
```julia-repl
julia> orbits(Perm(1,2)*Perm(4,5),1:5)
3-element Array{Array{Int16,1},1}:
 [1, 2]
 [3]
 [4, 5]
```
"""
function orbits(a::Perm,domain=1:length(a.d);trivial=true,check=false)
  cycles=Vector{eltype(a.d)}[]
  if isempty(a.d) return cycles end
  to_visit=falses(max(length(a.d),maximum(domain)))
  to_visit[domain].=true
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    cyc=orbit(a,i,check)
    to_visit[cyc].=false
    if length(cyc)>1 || trivial push!(cycles,cyc) end
  end
  cycles
end

# 15 times faster than GAP Cycles for rand(Perm,1000)
"""
  cycles(a::Perm) returns the non-trivial cycles of a
# Example
```julia-repl
julia> cycles(Perm(1,2)*Perm(4,5))
2-element Array{Array{Int16,1},1}:
 [1, 2]
 [4, 5]
```
"""
cycles(a::Perm;check=false)=orbits(a;trivial=false,check=check)

function Base.show(io::IO, a::Perm)
  cyc=orbits(a;trivial=false,check=true)
  if !get(io,:limit,false) print(io,"perm\"") end
  if isempty(cyc) print(io,"()")
  else for c in cyc print(io,"(",join(c,","),")") end
  end
  if !get(io,:limit,false) print(io,"\"") end
end

function Base.show(io::IO, ::MIME"text/plain", p::Perm{T})where T
 if T!=Int16 print(io,typeof(p),": ") end
  show(io,p)
end

"""
`order(a::Perm)` is the order of the permutation a
"""
order(a::Perm)=lcm(length.(cycles(a)))

"""
  cycletype(a::Perm) describes the partition of degree(a) associated to the
  conjugacy class of a in the symmetric group, with ones removed. It is
  represented as a list of pairs cyclesize=>multiplicity
# Example
```julia-repl
julia> cycletype(Perm(1,2)*Perm(3,4))
1-element Array{Pair{Int64,Int64},1}:
 2 => 2
```
"""
function cycletype(a::Perm;domain=nothing)
  res=Dict{Tuple{Int,Int},Int}()
  if isnothing(domain) to_visit=ones(length(a.d))
  else to_visit=domain[1:length(a.d)] # this makes a copy
  end
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
  if isnothing(domain) sort!(map(x->x[1][1]=>x[2],collect(res)))
  else sort!(collect(res))
  end
end

" reflength(a::Perm) minimum number of transpositions of which a is product"
function reflength(a::Perm)
  to_visit=ones(Bool,length(a.d))
  l=0
  for i in eachindex(to_visit)
@inbounds if !to_visit[i] continue end
    j=i
    while true
@inbounds to_visit[j]=false
@inbounds j=a.d[j]
      if j==i break end
      l+=1
    end
  end
  l
end

" sign(a::Perm) is the signature of  the permutation a"
Base.sign(a::Perm)=(-1)^reflength(a)
#---------------------- other -------------------------
"""
   restricted(a::Perm{T},l::AbstractVector{<:Integer})

l should be a union of cycles of p; returns p restricted to l

```julia-repl
julia> restricted(Perm(1,2)*Perm(3,4),3:4)
(3,4)
```
"""
function restricted(a::Perm{T},l::AbstractVector{<:Integer})where T
  o=orbits(a,l;trivial=false)
  isempty(o) ? Perm{T}() : prod(c->Perm{T}(c...),o)
end

"""
`mappingPerm(a,b)`

given two lists of positive integers without repetition `a` and `b`, this
function finds a permutation `p` such that `a.^p==b`.

```julia-repl
julia> mappingPerm([1,2,5,3],[2,3,4,6])
(1,2,3,6,5,4)
```
"""
function mappingPerm(::Type{T},a,b)where T
  r=1:max(maximum(a),maximum(b))
  a=vcat(a,setdiff(r,a))
  b=vcat(b,setdiff(r,b))
  Perm{T}(a)\Perm{T}(b)
end
mappingPerm(a,b)=mappingPerm(Int16,a,b)

end
