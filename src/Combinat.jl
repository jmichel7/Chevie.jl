module Combinat
export combinations, arrangements, partitions, npartitions, partition_tuples,
  conjugate_partition, compositions, submultisets, cartesian,
  npartition_tuples, narrangements, groupby, constant, tally, collectby,
  partitions_set, npartitions_set, bell, stirling2, unique_sorted!

"""
`groupby(v,l)`

group  elements of collection `l` according  to the corresponding values in
the coillection `v`

```julia-rep1
julia> groupby([31,28,31,30,31,30,31,31,30,31,30,31],
  [:Jan,:Feb,:Mar,:Apr,:May,:Jun,:Jul,:Aug,:Sep,:Oct,:Nov,:Dec])
Dict{Int64,Array{Symbol,1}} with 3 entries:
  31 => Symbol[:Jan, :Mar, :May, :Jul, :Aug, :Oct, :Dec]
  28 => Symbol[:Feb]
  30 => Symbol[:Apr, :Jun, :Sep, :Nov]
```
"""
function groupby(v,l)
  res=Dict{eltype(v),Vector{eltype(l)}}()
  for (k,val) in zip(v,l) push!(get!(res,k,empty(l)),val) end
  res
end

"""
`groupby(f::Function,l)`

group  elements of collection `l` according to the values taken by function
`f` on them

```julia-repl
julia> groupby(iseven,1:10)
Dict{Bool,Array{Int64,1}} with 2 entries:
  false => [1, 3, 5, 7, 9]
  true  => [2, 4, 6, 8, 10]
```
Note: `l` is required to be non-empty since I do not know how to access the
return type of a function
"""
function groupby(f::Function,l)
  if isempty(l) return Dict{Any,eltype(l)}() end
  res=Dict(f(l[1])=>[l[1]])
  for val in l[2:end]
    push!(get!(res,f(val),empty(l)),val)
  end
  res
end

"""
`tally(v)` 

count  how many times  each element of  collection `v` occurs  and return a
`Vector` of `elt=>count` (a variation on StatsBase.countmap)
"""
function tally(v)
  res=Pair{eltype(v),Int}[]
  if isempty(v) return res end
  if !(v isa AbstractArray) v=collect(v) end
  v=sort(v)
  c=1
@inbounds  for j in 2:length(v)
    if v[j]==v[j-1] c+=1
    else push!(res,v[j-1]=>c)
      c=1
    end
  end
  push!(res,v[end]=>c)
end

"""
`collectby(f,v)`

group  the elements of `v` in packets  (`Vector`s) where `f` takes the same
value. The resulting vector of vectors is sorted by the values of `f`.
"""
function collectby(f,v)
  d=groupby(f,v)
  [d[k] for k in sort(collect(keys(d)))]
end

" `constant(a)` whether all elements in collection `a` are equal"
constant(a)=isempty(a) || all(==(first(a)),a)

# faster than unique! for sorted vectors
function unique_sorted!(v::Vector)
  i=1
@inbounds  for j in 2:length(v)
    if v[j]<v[i] error("not sorted")
    elseif v[j]==v[i]
    else i+=1; v[i]=v[j]
    end
  end
  resize!(v,i)
end

# reverse to get the same order as GAP
function cartesian(a::AbstractVector...)
  reverse.(vec(collect.(Iterators.product(reverse(a)...))))
end

# k<=length(mset) && k>0
function combinations_sorted(sorted_mset,k)
  if k==1 return map(x->[x],sorted_mset) end
  reduce(vcat,map(enumerate(sorted_mset[1:end-k+1]))do (i,e)
    map(x->pushfirst!(x,e),combinations_sorted(sorted_mset[i+1:end],k-1))
  end)
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
function combinations(mset,k)
  if k>length(mset) return Vector{eltype(mset)}[] end
  if iszero(k) return [empty(mset)] end
  unique_sorted!(combinations_sorted(sort(mset),k))
end
combinations(mset)=isempty(mset) ? [Int[]] :
  vcat(combinations.(Ref(mset),0:length(mset))...)


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
function arrangements(mset,k)
  blist=trues(length(mset))
  function arr(mset,k)
    if iszero(k) return [empty(mset)] end
    combs=Vector{eltype(mset)}[]
    for i in eachindex(blist)
      if blist[i] && (i==length(blist) || mset[i+1]!=mset[i] || !blist[i+1])
        blist[i]=false
        append!(combs,pushfirst!.(arr(mset,k-1)::Vector{Vector{eltype(mset)}},
                                   Ref(mset[i])))
        blist[i]=true
      end
    end
    combs
  end
  arr(sort(mset),k)
end

arrangements(mset)=vcat(arrangements.(Ref(mset),0:length(mset))...)

function narrangements(mset,k)
  blist=trues(length(mset))
  function arr(mset,k)
    if iszero(k) return 1 end
    combs=0
@inbounds for i in eachindex(blist)
      if blist[i] && (i==length(blist) || mset[i+1]!=mset[i] || !blist[i+1])
        blist[i]=false
        combs+=arr(mset,k-1)
        blist[i]=true
      end
    end
    combs
  end
  arr(sort(mset),k)
end
narrangements(mset)=sum(narrangements.(Ref(mset),0:length(mset)))

@inbounds function partitions_less(v,n)
  l=v[end]
  k=length(l)
  for i in 1:(isempty(l) ? n : min(l[end],n))
    if i>1 push!(v,v[end][1:k]) end
    push!(v[end],i)
    if n>i partitions_less(v,n-i) end
  end
end

"""
`partitions(n)`

`partitions` returns the set of all partitions of the positive integer `n`.

A  *partition*  is  a  decomposition  `n=p₁+p₂+…+pₖ`  in integers such that
`p₁≥p₂≥…≥pₖ>0`,  and is represented by  the vector `p=[p₁,p₂,…,pₖ]`. We write
`p⊢n`. 

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
function partitions(n)
  v=[Int[]]
  partitions_less(v,n)
  v
end

# partitions of n of first (greatest) part <=m with k parts
function partitions_less(n,m,k)
# if m==1 return [fill(1,n)] end
  res=Vector{Int}[]
  if n<k return res end
  if k==1 return m<n ? res : [[n]] end
  for i in 1:min(m,n)
    append!(res,pushfirst!.(partitions_less(n-i,i,k-1),i))
  end
  res
end
"""
`partitions(n,k)`

The set of all partitions of the positive integer `n` with `k` parts.

```julia-repl
julia> partitions(7,3)
4-element Array{Array{Int64,1},1}:
 [3, 2, 2]
 [3, 3, 1]
 [4, 2, 1]
 [5, 1, 1]
```
"""
partitions(n,k)=partitions_less(n,n,k)

"""
`npartitions(n)` The number of partitions of `n`.
There are approximately `exp(π√(2n/3))/(4√3 n)` such partitions.
"""
function npartitions(n)
  s=one(n)
  p=fill(s,n+1)
  for m in 1:n
    s=zero(n)
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

" npartitions(n,k) number of partitions of n with k parts"
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

function partition_tuples_inner(v,n,r,l)
  k=length(v[end])
  start=true
  for i in (r==1 ? (n:n) : (0:n)), p in l[i+1]
    if !start push!(v,v[end][1:k]) end
    start=false
    push!(v[end],p)
    if r>1 partition_tuples_inner(v,n-i,r-1,l) end
  end
end

# the `r`-tuples of  partitions that together partition `n`.
function partition_tuples2(n,r)
  l=partitions.(0:n)
  v=[Vector{Int}[]]
  if iszero(r) if n>0 error("non-sensical 0-tuples of sum $n") end
  else partition_tuples_inner(v,n,r,l) end
  v
end

# bad implementation but which is ordered as GAP3; needed for
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

"`npartition_tuples(n,k)` number of `r`-tuples of  partitions of `n`."
function npartition_tuples(n,k)
  res=0
  for l in 1:min(n,k)
    r=binomial(k,l)
    res+=r*sum(a->narrangements(a,l)*prod(npartitions.(a)),partitions(n,l))
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
  if isempty(p) return p end
  res=zeros(eltype(p),maximum(p))
  for i in p, j in 1:i res[j]+=1 end
  res
end

"""
`compositions(n[,k];start=1)`

This  function returns the compositions of  `n` (the compositions of length
`k`  if a second argument `k` is given), where a composition of the integer
`n` is a decomposition `n=p₁+…+pₖ` in integers `≥start`, represented as the
vector `[p₁,…,pₖ]`. Unless `k` is given, `start` must be `>0`.

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

julia> compositions(4,2;start=0)
5-element Array{Array{Int64,1},1}:
 [4, 0]
 [3, 1]
 [2, 2]
 [1, 3]
 [0, 4]
```
"""
function compositions(n;start=1)
  if iszero(n) return [Int[]] end
  if start<=0 error("start must be ≥1") end
  vcat(map(i->map(c->push!(c,i),compositions(n-i;start)),start:n)...)
end

function compositions(n,k;start=1)
  if isone(k) return [[n]] end
  vcat(map(i->map(c->push!(c,i),compositions(n-i,k-1;start)),start:n-start)...)
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

"""
`partitions_set(set)`

the set of all unordered partitions of the set `set`.

An *unordered partition* of `set` is  a set of pairwise disjoint nonempty
sets with union `set`  and is represented by  a sorted list of such sets.

```julia-repl
julia> partitions_set(1:3)
5-element Array{Array{Array{Int64,1},1},1}:
 [[1], [2], [3]]
 [[1], [2, 3]]
 [[1, 2], [3]]
 [[1, 2, 3]]
 [[1, 3], [2]]
```

Note  that `partitions_set` does not currently support multisets and that
there is currently no ordered counterpart.
"""
function partitions_set(set)
  if unique(set)!=set error("partitions_set: ",set," must be a set") end
  if isempty(set) return [empty(set)] end
  m=1 .!=eachindex(set)
"""
`pp(part)`
all  partitions  of  `set`  that  begin  with `part[1:end-1]` and where the
`length(part)`-th set begins with `part[end]`.
 
To find that first we consider the set `part[end]` to be complete and add a
new  set to `part`, which must start with the smallest element of `set` not
yet  taken,  because  we  require  the  returned  partitions  to  be sorted
lexicographically.  `m` is  a boolean  list that  contains `true` for every
element of `set` not yet taken.

Second  we find all elements of `set`  that can be added to `part[end]` and
call  `pp`  recursively  for  each  candidate.  If  `o`  is the position of
`part[end][end]`  in `mset`, these candidates are set[l]` for `l` for which
`o<=l` and `m[l]` is `true`.
"""
  function pp(part)
    local l,o
    l=findfirst(m)
    if l===nothing return [part] end
    m[l]=false
    npart=copy.(part)
    push!(npart,[set[l]])
    parts=pp(npart)
    m[l]=true
    part=copy(part)
    part[end]=copy(part[end])
    o=findfirst(==(part[end][end]),set)
    push!(part[end],set[o])
    for l in o:length(m)
      if m[l]
        m[l]=false
        part[end][end]=set[l]
        append!(parts,pp(part))
        m[l]=true
      end
    end
    parts
  end
  pp([[set[1]]])
end

"""
`partitions_set(set,k)`
the set of  all unordered partitions of the  set `set` into  `k` pairwise
disjoint nonempty sets.

```julia-repl
julia> partitions_set(1:4,2)
7-element Array{Array{Array{Int64,1},1},1}:
 [[1], [2, 3, 4]]
 [[1, 2], [3, 4]]
 [[1, 2, 3], [4]]
 [[1, 2, 4], [3]]
 [[1, 3], [2, 4]]
 [[1, 3, 4], [2]]
 [[1, 4], [2, 3]]
```
"""
function partitions_set(set,k)
  if unique(set)!=set error("partitions_set: ",set," must be a set") end
  if isempty(set)
    if k==0  return [empty(set)]
    else return typeof(set)[]
    end
  end
  m=1 .!=eachindex(set)
"""
`pp(k,part)`
set  of  all  partitions  of  the  set  `set`  of  length  `n`,  that  have
`k+length(part)-1`  subsets, that begin with  `part[1:end-1]` and where the
`length(part)`-th set begins with `part[end]`.

To do so it does two things. It finds all elements of `mset` that can go at
`part[end][end+1]` and calls itself recursively for each candidate. And, if
`k`  is larger than 1, it considers  the set `part[end]` to be complete and
starts  a new set `part[end+1]`, which must start with the smallest element
of  `mset` not yet taken, because we  require the returned partitions to be
sorted  lexicographically. `m` is  a boolean list  that contains `true` for
every   element  of   `set`  not   yet  taken.   `o`  is  the  position  of
`part[end][end]` in `set`, so the candidates for `part[end][end+1]` are the
`set[l]` for which `o<l` and `m[l]` is `true`.
"""
  function pp(k,part)
    local l,o
    l=findfirst(m)
    parts=typeof(part)[]
    if k==1
      part=copy.(part)
      for l in k:length(set)
        if m[l] push!(part[end],set[l]) end
      end
      return [part]
    end
    if l===nothing return parts end
    m[l]=false
    npart=copy(part)
    push!(npart,[set[l]])
    parts=pp(k-1,npart)
    m[l]=true
    part=copy(part)
    part[end]=copy(part[end])
    o=findfirst(==(part[end][end]),set)
    push!(part[end],set[o])
    for l in o:length(set)
      if m[l]
        m[l]=false
        part[end][end]=set[l]
        append!(parts,pp(k,part));
        m[l]=true
      end
    end
    parts
  end
  pp(k,[[set[1]]])
end

"""
'bell(n)'

The  Bell numbers are  defined by `bell(0)=1`  and ``bell(n+1)=∑_{k=0}^n {n
\\choose  k}bell(k)``, or by the fact  that `bell(n)/n!` is the coefficient
of `x^n` in the formal series `e^(e^x-1)`.
    
```julia-repl
julia> bell.(0:6)
7-element Array{Int64,1}:
   1
   1
   2
   5
  15
  52
 203

julia> bell(14)
190899322
```julia-repl
"""
function bell(n)
  bell_=[1]
  for i in 1:n-1
    push!(bell_,bell_[1])
    for k in 0:i-1 bell_[i-k]+=bell_[i-k+1] end
  end
  bell_[1]
end

"""
`npartitions_set(set)`

the number of  unordered partitions of the set `set`.

```julia-repl
julia> npartitions_set(1:6)
203
```
"""
npartitions_set(set)=bell(length(set))

"""
`stirling2(n,k)`

the   *Stirling  number   of  the   second  kind*.   They  are  defined  by
`stirling2(0,0)=1`,  `stirling2(n,0)=stirling2(0,k)=0`  if  `n,  k!=0`  and
`stirling2(n,k)=k   stirling2(n-1,k)+stirling2(n-1,k-1)`,   and   also   as
coefficients of the generating function 
``x^n=\\sum_{k=0}^{n}stirling2(n,k) k!{x\\choose k}``.

```julia-repl
julia> stirling2.(4,0:4)
5-element Array{Int64,1}:  # Knuth calls this the trademark of stirling2
 0
 1
 7
 6
 1

julia> [stirling2(i,j) for i in 0:6, j in 0:6]
7×7 Array{Int64,2}:
 1  0   0   0   0   0  0 # Note the similarity with Pascal's triangle
 0  1   0   0   0   0  0
 0  1   1   0   0   0  0
 0  1   3   1   0   0  0
 0  1   7   6   1   0  0
 0  1  15  25  10   1  0
 0  1  31  90  65  15  1

julia> stirling2(50,big(10))
26154716515862881292012777396577993781727011
```
"""
function stirling2( n, k )
  if n<k return 0
  elseif n==k return 1
  elseif n<0 && k<0 return stirling1(-k,-n)
  elseif k<=0 return 0
  end
  bin,sti,fib=1,0,1
  for i in 1:k
    bin=div(bin*(k-i+1),i)
    sti=bin*i^n-sti
    fib*=i
  end
  div(sti,fib)
end

"""
`npartitions_set(set,k)`
the  number  of  unordered  partitions  of  the set `set` into `k` pairwise
disjoint  nonempty sets.

```julia-repl
julia> npartitions_set(1:10,3)
9330
```
"""
npartitions_set(set,k)=stirling2(length(set),k)
end
