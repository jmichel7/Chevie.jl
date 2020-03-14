module Combinat
export combinations, arrangements, partitions, NrPartitions, partition_tuples,
  conjugate_partition, evalpoly, dominates, compositions, submultisets,
  NrPartitionTuples, NrArrangements

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

'combinations'  returns the set of all combinations of the multiset `mset`.
If  a second argument `k` is given,  it returns the set of all combinations
of the multiset `mset` with `k` elements.

A *combination* of `mset` is an unordered selection without repetitions and
is  represented by a sorted  sublist of `mset`. If  `mset` is a proper set,
the set of all combinations is just the *powerset* of `mset`.

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

`arrangements` returns the set of arrangements of the multiset `mset`. It a
second argument `k` is given, it returns all arrangements with `k` elements
of the multiset `mset`.

An  *arrangement* of `mset`  is an ordered selection  without repetitions
and is represented by a list that contains only elements from `mset`, but
maybe  in a different  order.

As an example of arrangements of a multiset, think  of the game Scrabble.
Suppose you have the six characters of the word 'settle'  and you have to
make a four letter word.  Then the possibilities are given by

``julia-repl
julia> length(arrangements(collect("settle"),4))
102

julia> length(arrangements(collect("settle")))
523
```

Note that the fact that the list returned by 'arrangements' is a proper set
means  in this example that the possibilities  are listed in the same order
as they appear in the dictionary.
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

A  *partition*  is  a  sum  `n=p₁+p₂+…+  pₖ`  of  positive  integers and is
represented  by the  list `p=[p₁,p₂,…,pₖ]`,  in nonincreasing  order, i.e.,
`p₁≥p₂≥…≥pₖ`.  We write `p⊢n`.  There are approximately `exp(π√(2n/3))/(4√3
n)` such partitions.

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
NrPartitions(n...)=length(partitions(n...))

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

returns  the list of  all `r`-tuples of  partitions that together partition
`n`.

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

function NrPartitionTuples(n,k)
  res=0
  for l in 1:k
    r=binomial(k,l)
    res+=r*sum(a->NrArrangements(a,l)*prod(NrPartitions.(a)),partitions(n,l))
  end
  res
end

"""
`conjugate_partition(pi)`

returns the conjugate partition of the partition `pi`.

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

The  *conjugate  partition*  of  a  partition  `pi`  is  defined  to be the
partition belonging to the transposed of the Young diagram of `pi`.
"""
function conjugate_partition(p)
  res=zeros(eltype(p),maximum(p))
  for i in p, j in 1:i res[j]+=1 end
  res
end

# horner scheme
function evalpoly(x,p::Vector)
  value=zero(x)
  for i in length(p):-1:1
    value=x*value+p[i]
  end
  value
end

"""
`dominates(mu,nu)`

The  dominance  ordering  is  an  important partial order in representation
theory.  `dominates(mu,nu)` returns  `true` if  either `mu==nu`  or for all
`i≥1` we have `sumⱼ₌₁ⁱ muⱼ≥sumⱼ₌₁ⁱ nuⱼ`, and `false` otherwise.

```julia-repl
julia> dominates([5,4],[4,4,1])
true
```
"""
dominates(mu,nu)=all(i->i>length(nu) || sum(mu[1:i])>=sum(nu[1:i]),eachindex(mu))

"""
`compositions(n[,i])`

returns  the list of compositions of the integer `n` (the compositions with
`i` parts if a second argument `i` is given).

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

`submultisets` returns the  set of all  multisets length `k`
of elements of the set `set`.

An *multiset* of length `k` is a  selection with
repetitions  of length `k` from `set` and  is represented by a vector
of length `k` containing  elements  from  `set`.   There  are  
`binomial(|set|+k-1,k)` such sub-multisets.

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
