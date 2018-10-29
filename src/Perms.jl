"""
This module is a port of some GAP functionality on permutations.

A  permutation here is a permutation of the set 1:n and is represented as a
list  of n integers representing the images of 1:n. The integer n is called
the *degree* of the permutation.

Permutations  in  this  module  follow  the  GAP  design: it is possible to
multiply, or to store in the same group, permutations of different degrees.
A  slightly faster  design is  the MAGMA  one where  any permutation has to
belong  to  a  group  and  the  degree  is  determined by that group. There
multiplication of permutations in a given group is a faster, but it is more
difficult  to multiply  permutations coming  from different  groups, like a
group and one of its subgroups.

The  GAP permutation  (1,2,3)(4,5) can  be written Perm(1,2,3)*Perm(4,5) or
perm"(1,2,3)(4,5)".  It is represented internally as [2,3,1,5,4]; note that
[2,3,1,5,4,6] represents the same permutation.

As in GAP i^p applies p to integer i, while p^q means p^-1*q&ast;p.

Another  Perm  constructor  is  Perm{T}(p) which converts the
perm  p to a permutation on integers of type T; for instance Perm{UInt8} is
more  efficient that Perm{Int} and can be  used for Weyl groups of rank <=8
since they have at most 240 roots.

# Examples
```julia-repl
julia> p=Perm(1,2)*Perm(2,3)
(1,3,2)

julia> Perm{Int8}(p)
{Int8}(1,3,2)

julia> 1^p
3

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

julia> Matrix(p)
3×3 Array{Float64,2}:
 0.0  0.0  1.0
 1.0  0.0  0.0
 0.0  1.0  0.0

julia> Matrix{Int}(p)
3×3 Array{Int64,2}:
 0  0  1
 1  0  0
 0  1  0
```

Perms  have methods copy, hash,  ==, cmp, isless (total order)  so they can be
keys in hashes or elements of sets.

other functions are: cycles, cycletype, sign. See individual documentation.
"""
module Perms
export Perm, largest_moved_point, cycles, cycletype, order, sign,
  @perm_str, smallest_moved_point
using ..Gapjm # for degree

struct Perm{T<:Integer}
   d::Vector{T}
end

function Perm{T}(x::Int...)where T<:Integer
  if isempty(x) return Perm(T[]) end
  d=Vector{T}(collect(1:max(x...)))
  for i in 1:length(x)-1
    d[x[i]]=x[i+1]
  end
  d[x[end]]=x[1]
  Perm(d)
end

Perm(x::Int...)=Perm{Int}(x...)

function Perm{T}(p::Perm)where T<:Integer
  Perm(Vector{T}(p.d))
end

# just for fun, to provide Perm[1 2;5 6 7;4 9]=perm"(1,2)(5,6,7)(4,9)"
function Base.typed_hvcat(::Type{Perm},a::Tuple{Vararg{Int64,N} where N},
  b::Vararg{Number,N} where N)
  res=Perm()
  for i in a
    res*=Perm(Iterators.take(b,i)...)
    b=Iterators.drop(b,i)
  end
  res
end

macro perm_str(s::String)
  start=1
  res=Perm()
  while true
    m=match(r"\((\s*\d+\s*,)+\s*\d+\)",s[start:end])
    if m===nothing break end
    start+=m.match.ncodeunits
    res*=Perm(eval(Meta.parse(m.match))...)
  end
  res::Perm
end

Base.one(p::Perm{T}) where T=Perm(T[])
Base.one(::Type{Perm{T}}) where T=Perm(T[])

Gapjm.degree(a::Perm)= length(a.d)

@inline Base.:^(n::Integer, a::Perm{T}) where T=if n>degree(a) T(n) else a.d[n] end
Base.:^(a::Perm, b::Perm)=inv(b)*a*b
Base.:^(a::Perm, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                               Base.power_by_squaring(inv(a),-n)

Base.copy(p::Perm)=Perm(copy(p.d))

# hash is needed for using permutations in Sets/Dicts
function Base.hash(a::Perm, h::UInt)
   b = 0x595dee0e71d271d0%UInt
   for i in eachindex(a.d)
     if i^a!=i
       b = xor(b,xor(hash(i^a, h),h))
       b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
     end
   end
   return b
end

# total order is needed to use Perms in sorted lists
function Base.cmp(a::Perm, b::Perm)
  da=degree(a)
  db=degree(b)
  for i in 1:min(da,db)
@inbounds if a.d[i]<b.d[i] return -1 end
@inbounds if a.d[i]>b.d[i] return  1 end
  end
  if   da<db for i in (da+1:db) b.d[i]==i || return -1 end
  else       for i in (db+1:da) a.d[i]==i || return  1 end
  end
  0
end

Base.isless(a::Perm, b::Perm)=cmp(a,b)==-1

Base.:(==)(a::Perm, b::Perm)= cmp(a,b)==0

" largest_moved_point(a::Perm) is the largest integer moved by a"
largest_moved_point(a::Perm)=findlast(x->x^a!=x,a.d)
" smallest_moved_point(a::Perm) is the largest integer moved by a"
smallest_moved_point(a::Perm)=findfirst(k->k^a!=k,1:degree(a))

function Base.:*(b::Perm, a::Perm=one(b))
  if degree(b)<=degree(a)
    d=copy(a.d)
@inbounds @simd for i in 1:degree(b) d[i]=a.d[b.d[i]] end
  else
@inbounds d=[i^a for i in b.d]
  end
  Perm(d)
end

Base.inv(a::Perm)=Perm(invperm(a.d))

"""
  cycles(a::Perm) returns the non-trivial cycles of a
# Example
```julia-repl
julia> cycles(Perm(1,2)*Perm(4,5))
3-element Array{Array{Int64,1},1}:
 [1, 2]
 [3]
 [4, 5]
```
"""
function cycles(a::Perm{T},check::Bool=false)where T
  to_visit=trues(degree(a))
  cycles=Vector{Vector{T}}()
  for i in eachindex(to_visit)
   if !to_visit[i] continue end
    cycle=Vector{T}()
    j=i
    while true
      if check && j in cycle error("point $j occurs twice") end
      to_visit[j]=false
      push!(cycle,j)
      if (j^=a)==i break end
    end
    push!(cycles,cycle)
  end
  cycles
end

function Base.show(io::IO, a::Perm{T}) where T
  if T==Int t="" else t="{$T}" end
  cyc=(c for c in cycles(a) if length(c)>1)
  if isempty(cyc) print(io,t,"()")
  else print(io,t,join("("*join(c,",")*")" for c in cyc))
  end
end

order(a::Perm) = lcm(length.(cycles(a)))

"""
  cycletype(a::Perm) is the partition of degree(a) associated to the
  conjugacy class of a in the symmetric group, with ones removed
# Example
```julia-repl
julia> cycletype(Perm(1,2)*Perm(3,4))
2-element Array{Int64,1}:
 2
 2

```
"""
cycletype(a::Perm) = sort([length(c) for c in cycles(a)], rev=true)

" sign(a::Perm) is the signature of  the permutation a"
function Base.sign(a::Perm)
  parity = degree(a)
  to_visit = trues(parity)
  k = 1
  while (k=findnext(to_visit, k))!==nothing
     parity -= 1
     next=k
     while true
       to_visit[next]=false
       if (next^=a)==k break end
     end
  end
  (-1)^parity
end

using LinearAlgebra
"Matrix(a::Perm) is the permutation matrix for a"
Base.Matrix(a::Perm)=Matrix(1.0I,degree(a),degree(a))[a.d,:]
Base.Matrix{T}(a::Perm) where T= Matrix(one(T)I,degree(a),degree(a))[a.d,:]

end
