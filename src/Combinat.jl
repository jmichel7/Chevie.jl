module Combinat
export combinations, arrangements, partitions, npartitions, partition_tuples,
  conjugate_partition, dominates, compositions, submultisets, cartesian,
  npartition_tuples, NrArrangements, groupby, constant, tally, collectby

"""
  group items of list l according to the corresponding values in list v

    julia> groupby([31,28,31,30,31,30,31,31,30,31,30,31],
           [:Jan,:Feb,:Mar,:Apr,:May,:Jun,:Jul,:Aug,:Sep,:Oct,:Nov,:Dec])
    Dict{Int64,Array{Symbol,1}} with 3 entries:
      31 => Symbol[:Jan, :Mar, :May, :Jul, :Aug, :Oct, :Dec]
      28 => Symbol[:Feb]
      30 => Symbol[:Apr, :Jun, :Sep, :Nov]

"""
function groupby(v::AbstractArray{K},l::AbstractArray{V})where {K,V}
  res=Dict{K,Vector{V}}()
  for (k,val) in zip(v,l) push!(get!(res,k,V[]),val) end
  res
end

"""
  group items of list l according to the values taken by function f on them

    julia> groupby(iseven,1:10)
    Dict{Bool,Array{Int64,1}} with 2 entries:
      false => [1, 3, 5, 7, 9]
      true  => [2, 4, 6, 8, 10]

Note:in this version l is required to be non-empty since I do not know how to
access the return type of a function
"""
function groupby(f,l::AbstractArray)
  res=Dict(f(l[1])=>[l[1]]) # l should be nonempty
  for val in l[2:end]
    push!(get!(res,f(val),empty(l)),val)
  end
  res
end

"count how many times each element of v occurs and return a list of (elt,count)"
tally(v)=sort([(k,length(v)) for (k,v) in groupby(v,v)])

"group the elements of v in packets where f takes the same value"
function collectby(f,v)
  d=groupby(f,v)
  [d[k] for k in sort(collect(keys(d)))]
end

" whether all elements in list a are equal"
function constant(a::AbstractArray)
   all(i->a[i]==a[1],2:length(a))
end

# faster than unique! for sorted vectors
function unique_sorted(v::Vector)
  i=1
@inbounds  for j in 2:length(v)
    if v[j]==v[i]
    else i+=1; v[i]=v[j]
    end
  end
  resize!(v,i)
end

# reverse to get the same order as GAP
function cartesian(a::AbstractVector...)
  reverse.(vec(collect.(Iterators.product(reverse(a)...))))
end

function combinations_sorted(mset::AbstractVector,k)
  if iszero(k) return [eltype(mset)[]] end
  res=Vector{eltype(mset)}[]
  for (i,e) in enumerate(mset)
    append!(res,map(x->vcat([e],x),combinations_sorted(mset[i+1:end],k-1)))
  end
  res
end

"""
`combinations(mset[,k])`

'combinations'  returns  all  combinations  of  the  multiset `mset` (a not
necessarily  sorted list with  possible repetitions). If  a second argument
`k` is given, it returns the combinations with `k` elements.

A  *combination*  is  an  unordered  selection  without  repetitions and is
represented  by a sorted sublist of `mset`.  If `mset` is a proper set, the
set of all combinations is just the *powerset* of `mset`.

```julia-repl
julia> combinations([1,2,2,3])
12-element Array{Array{Int64,1},1}:
 []
 [1]
 [2]
 [3]
 [1, 2]
 [1, 3]
 [2, 2]
 [2, 3]
 [1, 2, 2]
 [1, 2, 3]
 [2, 2, 3]
 [1, 2, 2, 3]
```
"""
combinations(mset,k)=unique(combinations_sorted(sort(mset),k))
combinations(mset)=isempty(mset) ? [Int[]] :
   union(combinations.(Ref(mset),0:length(mset)))

function ArrangementsK(mset,blist,k)
  if iszero(k) return [eltype(mset)[]] end
  combs=Vector{eltype(mset)}[]
  for i in eachindex(mset)
    if blist[i] && (i==length(mset) || mset[i+1]!=mset[i] || !blist[i+1])
      blist1=copy(blist)
      blist1[i]=false
      append!(combs,pushfirst!.(ArrangementsK(mset,blist1,k-1),Ref(mset[i])))
    end
  end
  combs
end

"""
`arrangements(mset[,k])`

`arrangements`  returns  the  arrangements  of  the  multiset `mset` (a not
necessarily  sorted list with  possible repetitions). If  a second argument
`k` is given, it returns arrangements with `k` elements.

An  *arrangement* of  `mset` is  a selection  without repetitions  and thus
lists a subset `mset`, but in arbitrary order.

As  an example of arrangements  of a multiset, think  of the game Scrabble.
Suppose  you have the six  characters of the word  'settle' and you have to
make a four letter word. Then the possibilities are given by

```julia-repl
julia> length(arrangements(collect("settle"),4))
102

julia> length(arrangements(collect("settle")))
523
```
The  result  returned  by  'arrangements'  is  sorted,  which means in this
example  that the possibilities are listed in the same order as they appear
in the dictionary.
"""
arrangements(mset,k)=ArrangementsK(sort(mset),fill(true,length(mset)),k)
arrangements(mset)=isempty(mset) ? [Int[]] :
   union(arrangements.(Ref(mset),0:length(mset)))
NrArrangements(a...)=length(arrangements(a...))

# partitions of n of first (greatest) part <=m
function partitions_less(n,m)
  if m==1 return [fill(1,n)] end
  if iszero(n) return [Int[]] end
  res=Vector{Int}[]
  for i in 1:min(m,n)
    append!(res,map(x->vcat([i],x),partitions_less(n-i,i)))
  end
  res
end

"""
`partitions(n)`

`partitions` returns the set of all partitions of the positive integer `n`.

A  *partition*  is  a  decomposition  `n=p₁+p₂+…+pₖ`  in integers such that
`p₁≥p₂≥…≥pₖ>0`,  and is represented by  the list `p=[p₁,p₂,…,pₖ]`. We write
`p⊢n`. There are approximately `exp(π√(2n/3))/(4√3 n)` such partitions.

```julia-repl
julia> partitions(7)
15-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1, 1, 1, 1]
 [2, 1, 1, 1, 1, 1]
 [2, 2, 1, 1, 1]
 [2, 2, 2, 1]
 [3, 1, 1, 1, 1]
 [3, 2, 1, 1]
 [3, 2, 2]
 [3, 3, 1]
 [4, 1, 1, 1]
 [4, 2, 1]
 [4, 3]
 [5, 1, 1]
 [5, 2]
 [6, 1]
 [7]
```
"""
partitions(n)=partitions_less(n,n)
# partitions of n of first (greatest) part <=m with k parts
function partitions_less(n,m,k)
# if m==1 return [fill(1,n)] end
  res=Vector{Int}[]
  if n<k return res end
  if k==1 return m<n ? res : [[n]] end
  for i in 1:min(m,n)
    append!(res,map(x->vcat([i],x),partitions_less(n-i,i,k-1)))
  end
  res
end
partitions(n,k)=partitions_less(n,n,k)

function npartitions(n)
  s=1
  p=fill(s,n+1)
  for m in 1:n
    s=0
    k=1
    l=1
    while 0<=m-(l+k)
      s-=(-1)^k*(p[m-l+1]+p[m-(l+k)+1])
      k+=1
      l+=3*k-2
    end
    if l<=m s-=(-1)^k*p[m-l+1] end
    p[m+1]=s
  end
  s
end

function npartitions(n,k)
  if n==k return 1
  elseif n<k || k==0 return 0
  end
  p=fill(1,n)
  for l  in 2:k
    for m  in l+1:n-l+1 p[m]+=p[m-l] end
  end
  p[n-k+1]
end

if false
function partition_tuples(n,r)::Vector{Vector{Vector{Int}}}
  if r==1 return iszero(n) ? [[Int[]]] : map(x->[x],partitions(n)) end
  res=Vector{Vector{Int}}[]
  for i in  n:-1:1
    for p1 in partitions(i), p2 in partition_tuples(n-i,r-1)
      push!(res,vcat([p1],p2))
    end
  end
  for p2 in partition_tuples(n,r-1)
    push!(res,vcat([Int[]],p2))
  end
  res
end
else # bad implementation but which is ordered as GAP3; needed for
# compatibility with CHEVIE data library
"""
`partition_tuples(n,r)`

the `r`-tuples of  partitions that together partition `n`.

```julia-repl
julia> partition_tuples(3,2)
10-element Array{Array{Array{Int64,1},1},1}:
 [[1, 1, 1], []]
 [[1, 1], [1]]
 [[1], [1, 1]]
 [[], [1, 1, 1]]
 [[2, 1], []]
 [[1], [2]]
 [[2], [1]]
 [[], [2, 1]]
 [[3], []]
 [[], [3]]
```
"""
function partition_tuples(n, r)
   if n==0 return [fill(Int[],r)] end
   empty=(tup=[Int[] for i in 1:r], pos=fill(1,n-1))
   pm=[typeof(empty)[] for i in 1:n-1]
   for m in 1:div(n,2)
      for i in 1:r
         t1=map(copy,empty)
         t1.tup[i]=[m]
         t1.pos[m]=i
         push!(pm[m],t1)
      end
      for k in m+1:n-m
         for t in pm[k-m]
            for i in t.pos[m]:r
               t1=map(copy,t)
               t1.tup[i]=vcat([m],t1.tup[i])
               t1.pos[m]= i
               push!(pm[k], t1)
            end
         end
      end
   end
   res= Vector{Vector{Int}}[]
   for k in 1:n-1
      for t in pm[n-k]
         for i in t.pos[k]:r
            s=copy(t.tup)
            s[i]=vcat([k],s[i])
            push!(res,s)
         end
      end
   end
   for i in 1:r
      s=copy(empty.tup)
      s[i]=[n]
      push!(res,s)
   end
   res
end
end

function npartition_tuples(n,k)
  res=0
  for l in 1:k
    r=binomial(k,l)
    res+=r*sum(a->NrArrangements(a,l)*prod(npartitions.(a)),partitions(n,l))
  end
  res
end

"""
`conjugate_partition(λ)`

returns  the  conjugate  partition  of  the  partition  `λ`,  that  is, the
partition having the transposed of the Young diagram of `λ`.

```julia-repl
julia> conjugate_partition([4,2,1])
4-element Array{Int64,1}:
 3
 2
 1
 1

julia> conjugate_partition([6])
6-element Array{Int64,1}:
 1
 1
 1
 1
 1
 1
```
"""
function conjugate_partition(p)
  res=zeros(eltype(p),maximum(p))
  for i in p, j in 1:i res[j]+=1 end
  res
end

"""
`compositions(n[,k])`

This  function returns the compositions of  `n` (the compositions of length
`k`  if a second argument `k` is given), where a composition of the integer
`n` is a decomposition `n=p₁+…+pₖ` in positive integers, represented as the
list `[p₁,…,pₖ]`.

```julia-repl
julia> compositions(4)
8-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1]
 [2, 1, 1]
 [1, 2, 1]
 [3, 1]
 [1, 1, 2]
 [2, 2]
 [1, 3]
 [4]

julia> compositions(4,2)
3-element Array{Array{Int64,1},1}:
 [3, 1]
 [2, 2]
 [1, 3]
```
"""
function compositions(n)
  if iszero(n) return [Int[]] end
  vcat(map(i->map(c->push!(c,i),compositions(n-i)),1:n)...)
end

function compositions(n,k)
  if isone(k) return [[n]] end
  vcat(map(i->map(c->push!(c,i),compositions(n-i,k-1)),1:n-1)...)
end

"""
`submultisets(set,k)`

`submultisets`  returns the set of all  multisets of length `k` of elements
of the set `set`.

An  *multiset* of length `k` is a  selection with repetitions of length `k`
from `set` and is represented by a vector of length `k` containing elements
from `set`. There are `binomial(|set|+k-1,k)` such sub-multisets.

```julia-repl
julia> Combinat.submultisets(1:4,3)
20-element Array{Array{Int64,1},1}:
 [1, 1, 1]
 [1, 1, 2]
 [1, 1, 3]
 [1, 1, 4]
 [1, 2, 2]
 [1, 2, 3]
 [1, 2, 4]
 [1, 3, 3]
 [1, 3, 4]
 [1, 4, 4]
 [2, 2, 2]
 [2, 2, 3]
 [2, 2, 4]
 [2, 3, 3]
 [2, 3, 4]
 [2, 4, 4]
 [3, 3, 3]
 [3, 3, 4]
 [3, 4, 4]
 [4, 4, 4]
```
"""
function submultisets(A,i)
  if i==1 return map(x->[x],A) end
  vcat(map(j->map(v->pushfirst!(v,A[j]),submultisets(A[j:end],i-1)),
           eachindex(A))...)
end

end
