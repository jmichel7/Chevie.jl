"""
This package implements permutations and some functions of them. It depends
only  on  the  package  `Combinat` (which itself depends on `Primes`). 

This  package  follows  the  design  of  permutations  in the GAP language.
`Perm`s  are permutations  of the  set `1:n`,  represented internally  as a
vector  of `n`  integers holding  the images  of `1:n`.  The integer `n` is
called  the degree  of the  permutation. In  this package,  as in  GAP (and
contrary  to the philosophy of Magma or the package `Permutations.jl`), two
permutations of different  degrees  can  be  multiplied (the result has the
larger  degree). Two permutations  are equal if  and only if  they move the
same points in the same way, so two permutations of different degree can be
equal; the degree is thus an implementation detail so usually it should not
be used. One should rather use the function `largest_moved_point`.

This  design makes it  easy to multiply  permutations coming from different
groups, like a group and one of its subgroups. It has a negligible overhead
compared to the design where the degree is fixed.

The default constructor for a permutation uses the list of images of `1:n`,
like  `Perm([2,3,1,5,4])`.  Often  it  is  more  convenient  to  use  cycle
decompositions:    the   above   permutation    has   cycle   decomposition
`(1,2,3)(4,5)`    thus   can   be    written   `Perm(1,2,3)*Perm(4,5)`   or
`perm"(1,2,3)(4,5)"`  (this last form  can parse a  permutation coming from
GAP  or the default printing at the REPL).  The list of images of `1:n` can
be  recovered from the  permutation by the  function `vec`; note that equal
permutations with different degrees will have different `vec`.

The  complete type of a permutation  is `Perm{T}` where `T<:Integer`, where
`Vector{T}`  is the type of the vector which holds the image of `1:n`. This
can  be used to save space or  time. For instance `Perm{UInt8}` can be used
for  Weyl groups of rank≤8 since they permute  at most 240 roots. If `T` is
not  specified we  take it  to be  `Int16` since  this is a good compromise
between   speed,  compactness  and  possible  size  of  `n`.  One  can  mix
permutations of different integer types; they are promoted to the wider one
when multiplying.

# Examples of operations with permutations
```julia-repl
julia> a=Perm(1,2,3)
(1,2,3)

julia> vec(a)
3-element Vector{Int16}:
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

julia> one(a)  # trivial permutation
()

julia> sign(a) # signature of permutation
1

julia> order(a) # order (least trivial power) of permutation
3

julia> largest_moved_point(a)
3

julia> smallest_moved_point(a)
1

julia> Perm{Int8}(a) # convert a to Perm{Int8}
Perm{Int8}: (1,2,3)

julia> Matrix(b)  # permutation matrix of b
4×4 Matrix{Bool}:
 0  1  0  0
 0  0  1  0
 0  0  0  1
 1  0  0  0
```
```julia-rep1
julia> randPerm(10) # random permutation of 1:10
(1,8,4,2,9,7,5,10,3,6)
```

`Perm`s have methods `copy`, `hash`, `==`, so they can be keys in hashes or
elements  of sets; two permutations are equal  if they move the same points
to  the same images. They have methods `cmp`, `isless` (lexicographic order
on   moved  points)  so  they  can  be  sorted.  `Perm`s  are  scalars  for
broadcasting.

Other   methods   on   permutations   are  `cycles,  cycletype,  reflength,
mappingPerm, sortPerm, Perm_rowcol`.

No  method is given in  this package to enumerate  `Perm`s; you can use the
method   `arrangements`  from   `Combinat`  or   iterate  the  elements  of
`symmetric_group` with `PermGroups`.
"""
module Perms

export restricted, orbit, orbits, order, Perm, largest_moved_point, cycles,
  cycletype, support, @perm_str, smallest_moved_point, reflength,
  mappingPerm, sortPerm, Perm_rowcol, randPerm

using Combinat: tally, collectby, arrangements

"""
`struct Perm{T<:Integer}`

A  Perm represents a permutation  of the set `1:n`  and is implemented by a
`struct` with one field, a `Vector{T}` holding the images of `1:n`.

```julia-repl
julia> p=Perm(Int16[1,3,2,4])
(2,3)

julia> vec(p)
4-element Vector{Int16}:
 1
 3
 2
 4
```
"""
struct Perm{T<:Integer}
  d::Vector{T}
end
Perm(v)=Perm(collect(v))

const Idef=Int16 # you can change the default type T for Perm here

Base.vec(a::Perm)=a.d

#---------------- Constructors ---------------------------------------
"""
   `Perm{T}(x::Integer...)where T<:Integer`

   returns  a cycle.  For example  `Perm{Int8}(1,2,3)` constructs the cycle
   `(1,2,3)` as a `Perm{Int8}`. If omitted `{T}` is taken to be `{$Idef}`.
"""
function Perm{T}(x::Vararg{<:Integer,N};degree=0)where {T<:Integer,N}
  if isempty(x) return Perm(T[]) end
  d=T.(1:(degree==0 ? maximum(x) : degree))
  for i in 1:length(x)-1
    d[x[i]]=x[i+1]
  end
  d[x[end]]=x[1]
  Perm(d)
end

Perm(x::Integer...;degree=0)=Perm{Idef}(x...;degree)

Base.convert(::Type{Perm{T}},p::Perm{T1}) where {T,T1}=T==T1 ? p : Perm(T.(p.d))

"""
   `Perm{T}(p::Perm) where T<:Integer`

   change the type of `p` to `Perm{T}`
   for example `Perm{Int8}(Perm(1,2,3))==Perm{Int8}(1,2,3)`
"""
Perm{T}(p::Perm) where {T<:Integer}=convert(Perm{T},p)

"""
  @perm"..."

make a `Perm` from a string; allows GAP-style `perm"(1,2)(5,6,7)(4,9)"`
"""
macro perm_str(s::String)
  start=1
  s=replace(s,"\\n"=>"\n")
  res=Perm()
  while match(r"^\s*$"s,s[start:end])===nothing
    m=match(r"^\s*\((\s*\d+\s*,)+\s*\d+\)"s,s[start:end])
    if m===nothing error("malformed permutation: ",s) end
    start+=m.match.ncodeunits
    res*=Perm(Meta.parse(replace(m.match,r"\s*"=>"")).args...)
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
`Matrix(a::Perm,n=length(a.d))` 
the  permutation matrix  for `a`  operating on  `n` points (by default, the
degree of `a`). If given, `n` should be larger than `largest_moved_point(a)`.

```julia-repl
julia> Matrix(Perm(2,3,4),5)
5×5 Matrix{Bool}:
 1  0  0  0  0
 0  0  1  0  0
 0  0  0  1  0
 0  1  0  0  0
 0  0  0  0  1
```
"""
Base.Matrix(a::Perm,n=length(a.d))=[j==i^a for i in 1:n, j in 1:n]

"""
`Perm{T}(m::AbstractMatrix)`
If  `m` is a  permutation matrix, returns  the corresponding permutation of
type `T`. If omitted, `T` is taken to be `$Idef`.

```julia-repl
julia> m=[0 1 0;0 0 1;1 0 0]
3×3 Matrix{Int64}:
 0  1  0
 0  0  1
 1  0  0

julia> Perm(m)
(1,2,3)
```
"""
Perm(m::AbstractMatrix{<:Integer})=Perm{Idef}(m)

function Perm{T}(m::AbstractMatrix{<:Integer}) where T<:Integer
  if size(m,1)!=size(m,2) error("matrix should be square") end
  if any(x->count(!iszero,x)!=1,eachrow(m)) 
    error("not a permutation matrix")
  end
  l=map(x->findfirst(!iszero,x),eachrow(m))
  if !isperm(l) error("not a permutation matrix") end
  Perm{T}(l)
end

"""
  `Perm{T}(l::AbstractVector,l1::AbstractVector)`

returns `p`, a `Perm{T}`, such that `l1^p==l` if such a `p` exists; returns
`nothing` otherwise. If not given `{T}` is taken to be `{$Idef}`. Needs the
elements of `l` and `l1` to be sortable.

```julia-repl
julia> Perm([0,2,4],[4,0,2])
(1,3,2)
```
"""
function Perm{T}(l::AbstractVector,l1::AbstractVector)where T<:Integer
  p=sortperm(l)
  p1=sortperm(l1)
@inbounds if view(l,p)==view(l1,p1) Perm{T}(p1)\Perm{T}(p) end
end

Perm(l::AbstractVector,l1::AbstractVector)=Perm{Idef}(l,l1)

"""
  `Perm{T}(m::AbstractMatrix,m1::AbstractMatrix;dims=1)`

returns  `p`, a `Perm{T}`, which  permutes the rows of  `m1` (the coluns of
`m1`  if `dims=2`)  to bring  them to  those of  `m`, if such a `p` exists;
returns  `nothing` otherwise. If not given  `{T}` is taken to be `{$Idef}`.
Needs the elements of `m` and `m1` to be sortable.

```julia-repl
julia> Perm([0 1 0;0 0 1;1 0 0],[1 0 0;0 1 0;0 0 1];dims=1)
(1,3,2)

julia> Perm([0 1 0;0 0 1;1 0 0],[1 0 0;0 1 0;0 0 1];dims=2)
(1,2,3)
```
"""
function Perm{T}(l::AbstractMatrix,l1::AbstractMatrix;dims=1)where T<:Integer
  if     dims==1 Perm{T}(collect(eachrow(l)),collect(eachrow(l1)))
  elseif dims==2 Perm{T}(collect(eachcol(l)),collect(eachcol(l1)))
  end
end

Perm(l::AbstractMatrix,l1::AbstractMatrix;dims=1)=Perm{Idef}(l,l1,dims=dims)

#---------------------------------------------------------------------
Base.one(p::Perm)=Perm(empty(p.d))
Base.one(::Type{Perm{T}}) where T=Perm(T[])
Base.copy(p::Perm)=Perm(copy(p.d))

# hash is needed for using Perms in Sets or as keys in Dicts
function Base.hash(a::Perm, h::UInt)
  for (i,v) in pairs(a.d)
    if v!=i h=hash(v,h) end
  end
  h
end

# Perms are scalars for broadcasting
Base.broadcastable(p::Perm)=Ref(p)

Base.typeinfo_implicit(::Type{Perm{T}}) where T=T==Idef

function Base.promote_rule(a::Type{Perm{T1}},b::Type{Perm{T2}})where {T1,T2}
  Perm{promote_type(T1,T2)}
end

function extend!(a::Perm{T},n::Integer)where T
  if length(a.d)<n append!(a.d,T(length(a.d)+1):T(n)) end
end

"""
 `promote_degree(a::Perm, b::Perm)` promotes `a` and `b` to the same type,
 then extends `a` and `b` to the same degree
"""
function promote_degree(a::Perm,b::Perm)
  a,b=promote(a,b)
  extend!(a,length(b.d))
  extend!(b,length(a.d))
  (a,b)
end

# permutations need to be totally ordered to use them in sorted lists
function Base.isless(a::Perm, b::Perm)
  a,b=promote_degree(a,b)
  isless(a.d,b.d)
end

function Base.:(==)(a::Perm, b::Perm)
  a,b=promote_degree(a,b)
  a.d==b.d
end

" `largest_moved_point(a::Perm)` is the largest integer moved by a"
largest_moved_point(a::Perm)=findlast(x->a.d[x]!=x,eachindex(a.d))

" `smallest_moved_point(a::Perm)` is the smallest integer moved by a"
smallest_moved_point(a::Perm)=findfirst(x->a.d[x]!=x,eachindex(a.d))

" `support(a::Perm)` is the set of all points moved by `a`"
support(a::Perm)=eachindex(a.d)[a.d.!=eachindex(a.d)]
# 5 times faster than findall(x->a.d[x]!=x,eachindex(a.d))

" for convenience: `sortPerm(a)=Perm(sortperm(a))`"
sortPerm(::Type{T},a::AbstractVector;k...) where T=Perm{T}(sortperm(a;k...))
sortPerm(a::AbstractVector;k...)=sortPerm(Idef,a;k...)

"""
`randPerm([T,]n::Integer)` a random permutation of `1:n` of type `T`.
If omitted `T` is taken to be `$Idef`
"""
randPerm(::Type{T},n::Integer) where T =sortPerm(T,rand(1:n,n))
randPerm(n::Integer)=sortPerm(Idef,rand(1:n,n))

#------------------ operations on permutations --------------------------

function Base.:*(a::Perm, b::Perm)
  a,b=promote_degree(a,b)
  r=similar(a.d)
@inbounds for (i,v) in pairs(a.d) r[i]=b.d[v] end
  Perm(r)
end

# this is a*=b without allocation
function mul!(a::Perm, b::Perm)
  a,b=promote_degree(a,b)
@inbounds for (i,v) in pairs(a.d) a.d[i]=b.d[v] end
  a
end

function Base.inv(a::Perm)
  r=similar(a.d)
@inbounds for (i,v) in pairs(a.d) r[v]=i end
  Perm(r)
end

# less allocations than inv(a)*b
function Base.:\(a::Perm, b::Perm)
  a,b=promote_degree(a,b)
  r=similar(a.d)
@inbounds for (i,v) in pairs(a.d) r[v]=b.d[i] end
  Perm(r)
end

# less allocations than inv(b)*a*b
function Base.:^(a::Perm, b::Perm)
  a,b=promote_degree(a,b)
  r=similar(a.d)
@inbounds for (i,v) in pairs(a.d) r[b.d[i]]=b.d[v] end
  Perm(r)
end

# I do not know how to do this one faster
Base.:/(a::Perm, b::Perm)=a*inv(b)

@inline Base.:^(n::Integer, a::Perm{T}) where T=
  n>length(a.d) ? T(n) : @inbounds a.d[n] 

Base.:^(a::Perm, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                    Base.power_by_squaring(inv(a),-n)

"""
`Base.:^(l::AbstractVector,p::Perm)` 

returns `l` permuted by `p`, a vector `r` such that `r[i^p]==l[i]`

# Examples
```julia-repl
julia> [5,4,6,1,7,5]^Perm(1,3,5,6,4)
6-element Vector{Int64}:
 1
 4
 5
 5
 6
 7
```

note that we follow here the convention for the GAP function `Permuted`, but
this has the consequence that `sort(a)==a^inv(Perm(sortperm(a)))`.
"""
function Base.:^(l::AbstractVector,a::Perm)
  res=similar(l)
@inbounds for i in eachindex(l) res[i^a]=l[i] end
  res
end

"""
`Base.:^(m::AbstractMatrix,p::Perm;dims=1)`

Applies the permutation `p` on the lines, columns or both of the matrix `m`
depending on the value of `dims`

```julia-repl
julia> m=[3*i+j for i in 0:2,j in 1:3]
3×3 Matrix{Int64}:
 1  2  3
 4  5  6
 7  8  9

julia> p=Perm(1,2,3)
(1,2,3)

julia> m^p
3×3 Matrix{Int64}:
 7  8  9
 1  2  3
 4  5  6

julia> ^(m,p;dims=2)
3×3 Matrix{Int64}:
 3  1  2
 6  4  5
 9  7  8

julia> ^(m,p;dims=(1,2))
3×3 Matrix{Int64}:
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

"""
`Base.:^(m::AbstractMatrix,p::Tuple{Perm,Perm})`

given  a tuple `(p1,p2)` of  `Perm`s, applies `p1` to  the lines of `m` and
`p2` to the columns of `m`.

```julia-repl
julia> m=[1 2 3;4 5 6;7 8 9]
3×3 Matrix{Int64}:
 1  2  3
 4  5  6
 7  8  9

julia> m^(Perm(1,2),Perm(2,3))
3×3 Matrix{Int64}:
 4  6  5
 1  3  2
 7  9  8
```
"""
Base.:^(m::AbstractMatrix,p::Tuple{Perm,Perm})=m[axes(m,1)^p[1],axes(m,2)^p[2]]
#---------------------- cycles -------------------------

# 20% slower than GAP CyclePermInt for randPerm(1000)
"""
  orbit(a::Perm,i::Integer) returns the orbit of a on i
"""
function orbit(a::Perm,i::Integer)
  res=empty(a.d)
  sizehint!(res,length(a.d))
  j=i
  while true
    push!(res,j)
    j^=a
    if j==i return res end
  end
end
  
"""
`orbits(a::Perm,d::Vector=1:length(a.d))` 

returns the orbits of `a` on domain `d`

# Example
```julia-repl
julia> orbits(Perm(1,2)*Perm(4,5),1:5)
3-element Vector{Vector{$Idef}}:
 [1, 2]
 [3]
 [4, 5]
```
"""
function orbits(a::Perm,domain=1:length(a.d);trivial=true)
  cycles=Vector{eltype(a.d)}[]
  if isempty(a.d) return cycles end
  to_visit=falses(max(length(a.d),maximum(domain)))
@inbounds  to_visit[domain].=true
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    cyc=orbit(a,i)
@inbounds  to_visit[cyc].=false
    if length(cyc)>1 || trivial push!(cycles,cyc) end
  end
  cycles
end

# 15 times faster than GAP Cycles for randPerm(1000)
"""
  `cycles(a::Perm)` returns the non-trivial cycles of `a`
# Example
```julia-repl
julia> cycles(Perm(1,2)*Perm(4,5))
2-element Vector{Vector{$Idef}}:
 [1, 2]
 [4, 5]
```
"""
cycles(a::Perm)=orbits(a;trivial=false)

function Base.show(io::IO, a::Perm)
  if !isperm(a.d) error("malformed permutation") end
  replorTeX=get(io,:limit,false)||get(io,:TeX,false)
  if !replorTeX print(io,"perm\"") end
  cyc=cycles(a)
  if isempty(cyc) print(io,"()")
  else for c in cyc print(io,"(",join(c,","),")") end
  end
  if !replorTeX print(io,"\"") end
end

function Base.show(io::IO, ::MIME"text/plain", p::Perm{T})where T
  if T!=Idef && !haskey(io,:typeinfo) print(io,typeof(p),": ") end
  show(io,p)
end

"""
`cycletype(a::Perm;domain=1:length(a.d),trivial=false)`

returns  the  partition  of  `maximum(domain)`  associated to the conjugacy
class of `a` in the symmetric group of `domain`, with ones removed (thus it
does  not  depend  on  `domain`  but  just  on  the  moved  points)  unless
`trivial=true`.

# Example
```julia-repl
julia> cycletype(Perm(1,2)*Perm(4,5))
2-element Vector{Int64}:
 2
 2

julia> cycletype(Perm(1,2)*Perm(4,5);trivial=true)
3-element Vector{Int64}:
 2
 2
 1

julia> cycletype(Perm(1,2)*Perm(4,5);trivial=true,domain=1:6)
4-element Vector{Int64}:
 2
 2
 1
 1
```
"""
function cycletype(a::Perm;domain=1:length(a.d),trivial=false)
  lengths=Int[]
  if isempty(domain) return lengths end
  to_visit=falses(max(length(a.d),maximum(domain)))
@inbounds  to_visit[domain].=true
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    j=i
    l=0
    while true
      l+=1
      to_visit[j]=false
      j^=a
      if j==i break end
    end
    if l>1 || trivial push!(lengths,l) end
  end
  sort(lengths,rev=true)
end

"""
`order(a::Perm)` is the order of the permutation a
"""
order(a::Perm)=lcm(cycletype(a))

"""
`reflength(a::Perm)`

is   the  "reflection   length"  of   `a`,  that   is,  minimum  number  of
transpositions of which `a` is the product
"""
function reflength(a::Perm)
  to_visit=trues(length(a.d))
  l=0
  for i in eachindex(to_visit)
@inbounds if !to_visit[i] continue end
    j=i
    while true
@inbounds to_visit[j]=false
@inbounds j=Int(a.d[j])
      if j==i break end
      l+=1
    end
  end
  l
end

" `sign(a::Perm)` is the signature of  the permutation `a`"
Base.sign(a::Perm)=(-1)^reflength(a)

"""
`restricted(a::Perm,l::AbstractVector{<:Integer})`

`l` should be a union of cycles of `p`; returns `p` restricted to `l`

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
`mappingPerm(a)`

given  a list  of positive  integers without  repetition `a`, this function
finds  a permutation  `p` such  that `sort(a).^p==a`.  This can  be used to
translate between arrangements and `Perm`s.

```julia-repl
julia> p=mappingPerm([6,7,5])
(5,6,7)

julia> (5:7).^p
3-element Vector{$Idef}:
 6
 7
 5
```
"""
function mappingPerm(::Type{T},a)where T
  r=collect(1:maximum(a))
  r[sort(a)]=a
  Perm{T}(r)
end
mappingPerm(a)=mappingPerm(Idef,a)

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
mappingPerm(a,b)=mappingPerm(Idef,a,b)

"""
`Perm_rowcol(m1::AbstractMatrix, m2::AbstractMatrix)`

whether `m1` is conjugate to `m2` by row/col permutations.

`m1` and `m2` should be rectangular matrices of the same size. The function
returns a pair of permutations `(p1,p2)` such that `m1^(p1,p2)==m2` if such
permutations exist, `nothing` otherwise.

The entries of `m1` and `m2` must be sortable.

```julia-repl
julia> a=[1 1 1 -1 -1; 2 0 -2 0 0; 1 -1 1 -1 1; 1 1 1 1 1; 1 -1 1 1 -1]
5×5 Matrix{Int64}:
 1   1   1  -1  -1
 2   0  -2   0   0
 1  -1   1  -1   1
 1   1   1   1   1
 1  -1   1   1  -1

julia> b=[1 -1 -1 1 1; 1 1 -1 -1 1; 1 -1 1 -1 1; 2 0 0 0 -2; 1 1 1 1 1]
5×5 Matrix{Int64}:
 1  -1  -1   1   1
 1   1  -1  -1   1
 1  -1   1  -1   1
 2   0   0   0  -2
 1   1   1   1   1

julia> p1,p2=Perm_rowcol(a,b)
((1,2,4,5,3), (3,5,4))

julia> a^(p1,p2)==b
true
```
"""
function Perm_rowcol(m1::AbstractMatrix, m2::AbstractMatrix;debug=false)
  if size(m1)!=size(m2) error("not same dimensions") end
  if isempty(m1) return [Perm(), Perm()] end
  dist(m,n)=count(i->m[i]!=n[i],eachindex(m))
  dist(m,n,dim,l)=dim==1 ? dist(m[l,:],n[l,:]) : dist(m[:,l],n[:,l])
  mm=[m1,m2]
  if debug print("# ", dist(m1, m2), "") end
  rcperm=[Perm(), Perm()],[Perm(), Perm()]
  crg=Vector{Int}[],Vector{Int}[]
  crg1=[axes(m1,1)],[axes(m1,2)]
  while true
    crg=crg1
    crg1=Vector{Int}[],Vector{Int}[]
    for dim in 1:2
      for g in crg[dim]
        invars=map(1:2) do i
          invar=map(j->map(k->tally(dim==1 ? mm[i][j,k] : mm[i][k,j]),
                           crg[3-dim]), g)
          p=mappingPerm(vcat(collectby(invar,g)...), g)
          rcperm[dim][i]*=p
          mm[i]=^(mm[i],p,dims=dim)
          sort!(invar)
        end
        if invars[1]!=invars[2] return nothing end
        append!(crg1[dim], collectby(invars[1],g))
      end
    end
    if debug print("==>",dist(mm[1],mm[2])) end
    if crg==crg1 break end
  end
  function best(l,dim)
    if length(l)==1 return false end
    d=dist(mm[1], mm[2], dim, l)
#   if debug print("l=",l,"\n") end
    for e in mappingPerm.(arrangements(l,length(l)))
      m=dist(^(mm[1], e;dims=dim), mm[2], dim, l)
      if m<d
        if debug print("\n",("rows","cols")[dim],l,":$d->",m) end
        rcperm[dim][1]*=e
        mm[1]=^(mm[1],e;dims=dim)
        return true
      end
    end
    return false
  end
  while true
    s=false
    for dim in 1:2 for g in crg[dim] s=s || best(g,dim) end end
    if !s break end
  end
  if debug print("\n") end
  if !iszero(dist(mm...)) error("Perm_rowcol failed") end
  (rcperm[1][1]/rcperm[1][2],rcperm[2][1]/rcperm[2][2])
end

end
