"""
Posets  are represented as objects where at  least one of the two following
fields is present:

  - `incidence`:  a  boolean  matrix  such that `incidence[i][j]==true` iff
   `i<=j` in the poset.

  - `hasse`:  a list representing  the Hasse diagram  of the poset: the i-th
    entry is the list of indices of elements which are immediate successors
    (covers) of the i-th element, that is the list of j such that `i<j` and
    such that there is no k such that i<k<j.

There   are  thus  two   constructors,  `Poset(I::Matrix{Bool})`  from  the
incidence  matrix `I`, and `Poset(H::Vector{<:Vector{<:Integer}})` from the
Hasse diagram.

`Poset`s  are  printed  as  a  list  of covering chains. Elements which are
equivalent for the `Poset` are printed together separated by commas.

```julia-repl
julia> p=Poset(coxgroup(:A,2))
.<1,2<21,12<121

julia> length(p) # the number of elements of the `Poset`
6
```

note in the above example that the poset has been printed with labels which
are  the words in `coxgroup(:A,2)`. If no labels have been set for `p`, the
labels  shown are `1:length(p)`. The labels should have same length as  `p`
but can be any objects.

```julia-repl
julia> p.labels="abcdef"; p
a<b,c<d,e<f

julia> hasse(p)
6-element Vector{Vector{Int64}}:
 [2, 3]
 [4, 5]
 [4, 5]
 [6]   
 [6]   
 []    

julia> incidence(p)
6×6 Matrix{Bool}:
 1  1  1  1  1  1
 0  1  0  1  1  1
 0  0  1  1  1  1
 0  0  0  1  0  1
 0  0  0  0  1  1
 0  0  0  0  0  1
```
"""
module Posets 
# this module has only the 2 dependencies below which could be copied (27 lines)
using ..Combinat: collectby
using ..Util: @GapObj
export Poset, linear_extension, hasse, incidence, partition,
transitive_closure, is_join_lattice, is_meet_lattice
export restricted #import ..Gapjm: restricted to not use using_merge

@GapObj struct Poset end

"""
`transitive_closure(M)`

`M`  should be a  square boolean matrix  representing a relation; returns a
boolean  matrix representing the  transitive closure of  this relation. The
transitive  closure is computed  by the Floyd-Warshall  algorithm, which is
quite fast even for large matrices.

```julia-repl
julia> m=[j-i in [0,1] for i in 1:5, j in 1:5]
5×5 Matrix{Bool}:
 1  1  0  0  0
 0  1  1  0  0
 0  0  1  1  0
 0  0  0  1  1
 0  0  0  0  1

julia>transitive_closure(m)
5×5 Matrix{Bool}:
 1  1  1  1  1
 0  1  1  1  1
 0  0  1  1  1
 0  0  0  1  1
 0  0  0  0  1
```
"""
function transitive_closure(m)
  m=copy(m)
  for k in axes(m,1), i in axes(m,1)
   if m[i,k] m[i,:].|=m[k,:] end
  end
  m
end

"""
`Poset(m::Matrix{Bool})`

Creates a poset from an incidence matrix `m`, that is `m[i,j]==true` if and
only if `i≤j` in the poset,

```julia-repl
julia> Poset(Bool[1 1 1 1 1;0 1 0 1 1;0 0 1 1 1;0 0 0 1 0;0 0 0 0 1])
1<2,3<4,5
```
"""
Poset(m::Matrix{Bool})=Poset(Dict{Symbol,Any}(:incidence=>m))

"""
`Poset(h::Vector{<:Vector{<:Integer}})`

Creates a poset from a Hasse diagram given as a `Vector` whose `i`-th entry
is  the list of indices of elements which are immediate successors (covers)
of the `i`-th element, that is `h[i]` is the list of `j` such that `i<j` in
the poset and such that there is no `k` such that `i<k<j`.

```julia-repl
julia> Poset([[2,3],[4,5],[4,5],Int[],Int[]])
1<2,3<4,5
```
"""
Poset(m::Vector{<:Vector{<:Integer}})=Poset(Dict{Symbol,Any}(:hasse=>m))
Base.length(p::Poset)=haskey(p,:hasse) ? length(hasse(p)) : size(incidence(p),1)

function label(io::IO,p::Poset,n)
  haskey(p,:labels) ? p.labels[n] : (haskey(p,:label) ? p.label(io,n) : string(n))
end
  
function Base.show(io::IO,x::Poset)
  s=hasse(x)
  pp=partition(x)
  labels=map(y->join(map(n->label(io,x,n),y),","),pp)
  p=Poset(map(x->unique!(sort(convert(Vector{Int},
                        map(y->findfirst(z->y in z,pp),s[x[1]])))),pp))
  sep=get(io,:Symbol,false)
  TeX=get(io,:TeX,false)
  if sep==false sep=TeX ? "{<}" : "<" end
  ch=map(x->join(labels[x],sep), chains(p)) 
  if TeX print(io,join(map(x->"\$\$$x\$\$\n",ch)))
  else join(io,ch,"\n")
  end
end

"""
`linear_extension(P)`

returns  a  linear  extension  of  the  poset  `P`,  that  is  a vector `l`
containing  a permutation of the integers  `1:length(P)` such that if `i<j`
in `P`, then `findfirst(==(i),l)<findfirst(==(j),l)`. This is also called a
topological sort of `P`.

```julia-repl
julia> p=Poset([j%i==0 for i in 1:6, j in 1:6])
1<5
1<2<4
1<3<6
2<6

julia> linear_extension(p)
6-element Vector{Int64}:
 1
 2
 3
 5
 4
 6
```
"""
function linear_extension(P::Poset)
  ord=hasse(P)
  n=zeros(length(ord))
  for v in ord for x in v n[x] += 1 end end
  Q=filter(x->iszero(n[x]),1:length(n))
  res=Int[]
  while length(Q) > 0
    push!(res, Q[1])
    for x in ord[Q[1]]
      n[x]-=1
      if iszero(n[x]) push!(Q, x) end
    end
    popfirst!(Q)
  end
  if sum(n)>0 error("cycle") end
  res
end

hasse(m::Matrix{Bool})=map(x->filter(y->x[y]==2,1:length(x)),eachrow(m*m))

"""
`hasse(P)`

returns the Hasse diagram of the poset `P`.

```julia-repl
julia> p=Poset([j%i==0 for i in 1:5, j in 1:5])
1<3,5
1<2<4

julia> hasse(p)
5-element Vector{Vector{Int64}}:
 [2, 3, 5]
 [4]      
 []       
 []       
 []       
```
"""
hasse(p::Poset)=get!(()->hasse(incidence(p)),p,:hasse)

"""
`incidence(P)`

returns the incidence matrix of the poset `P`.

```julia-repl
julia> p=Poset(push!([[i+1] for i in 1:5],Int[]))
1<2<3<4<5<6

julia> incidence(p)
6×6 Matrix{Bool}:
 1  1  1  1  1  1
 0  1  1  1  1  1
 0  0  1  1  1  1
 0  0  0  1  1  1
 0  0  0  0  1  1
 0  0  0  0  0  1
```
"""
function incidence(p::Poset)::Matrix{Bool}
  get!(p,:incidence)do
    n=linear_extension(p)
    incidence=one(Matrix{Bool}(undef,length(n),length(n)))
    for i in length(n)-1:-1:1
      for x in hasse(p)[n[i]] incidence[n[i],:].|= incidence[x,:] end
    end
    incidence
  end
end

"A (greedy: the first is longest possible) list of covering chains of P."
function chains(P::Poset)
  ch=Vector{Int}[]
  h=hasse(P)
  for i in linear_extension(P), j in h[i]
    p=findfirst(c->i==c[end],ch)
    if p===nothing push!(ch,[i,j])
    else push!(ch[p],j)
    end
  end
  ch
end

"""
`reverse(P)`

the opposed poset to `P`.

```julia-repl
julia> p=Poset([i==j || i%4<j%4 for i in 1:8, j in 1:8])
4,8<1,5<2,6<3,7

julia> reverse(p)
3,7<2,6<1,5<4,8
```
"""
function Base.reverse(p::Poset)
  res=deepcopy(p)
  if haskey(p,:incidence) res.incidence=permutedims(incidence(p)) end
  if haskey(p,:hasse)
    res.hasse=map(empty,hasse(p))
    for i in 1:length(p), j in hasse(p)[i] push!(hasse(res)[j], i) end
  end
  return res
end

"""
`partition(P::Poset)`

returns  the  partition  of  `1:length(P)`  determined  by  the equivalence
relation  associated to `P`; that  is, `i` and `j`  are in the same part of
the  partition if the `k` such that `i<k` and `j<k` are the same as well as
the `k` such that `k<i` and `k<j`.

```julia-repl
julia> p=Poset([i==j || i%4<j%4 for i in 1:8, j in 1:8])
4,8<1,5<2,6<3,7

julia> partition(p)
4-element Vector{Vector{Int64}}:
 [4, 8]
 [2, 6]
 [3, 7]
 [1, 5]
```
"""
function partition(p::Poset)
  if haskey(p,:hasse)
    l=reverse(p)
    collectby(i->(hasse(l)[i], hasse(p)[i]),1:length(p))
  else
    I=incidence(p)
    ind=1:length(p)
    l=map(i->[map(j->j!=i && I[i,j], ind), map(j->j!=i && I[j,i], ind)], ind)
    map(x->filter(i->l[i] == x,ind), unique!(sort(l)))
  end
end

"""
`restricted(P::Poset,inds::AbstractVector{<:Integer})`

returns  the sub-poset of `P` determined by `inds`, which must be a sublist
of `1:length(P)`. The indices in this sub-poset will be renumbered to their
position in `inds`. If `P` has labels, they will be transmitted to the
sub-poset.

```julia-repl
julia> p=Poset([i==j || i%4<j%4 for i in 1:8, j in 1:8])
4,8<1,5<2,6<3,7

julia> restricted(p,2:6)
3<4<1,5<2

julia> p.labels=1:8; restricted(p,2:6)
4<5<2,6<3
```
"""
function restricted(p::Poset,ind::AbstractVector{<:Integer})
  res=Poset(copy(p.prop))
  if length(ind) == length(p) && sort(ind) == 1:length(p)
    if haskey(res,:hasse)
      res.hasse=Vector{Int}.(map(x->map(y->findfirst(==(y),ind),x),res.hasse[ind]))
    end
    if haskey(res, :incidence) res.incidence=incidence(res)[ind,ind] end
  else
    inc=incidence(p)
    res.incidence=[i!=j && ind[i]==ind[j] ? false : inc[ind[i],ind[j]]
                            for i in eachindex(ind), j in eachindex(ind)]
    delete!(res.prop, :hasse)
  end
  if haskey(p, :label) res.label=(io,x, n)->p.label(io, x, ind[n])
  elseif haskey(p, :labels) res.labels=p.labels[ind]
  end
  res
end

function checkl(ord::Matrix{Bool})
  subl=Set(Vector{Bool}[])
  n=size(ord,1)
  for i in 1:n, j in 1:i-1
    if !ord[j,i] || ord[i,j]
      l=ord[i,:].&ord[j,:]
      if !(l in subl)
        if !any(y->l==l.&y,ord[l,:])
          for k in (1:n)[l]
            ll=copy(ord[k,:])
            ll[k]=false
            l.&=.!ll
          end
          println("# $i  &&  $j have bounds ", (1:n)[l])
          return false
        end
        push!(subl, l)
      end
    end
  end
  true
end

"""
`is_join_lattice(P::Poset)`

returns  `true` if `P` is  a join semilattice, that  is any two elements of
`P` have a unique smallest upper bound; returns `false` otherwise.

```julia-repl
julia> p=Poset([i==j || i%4<j%4 for i in 1:8, j in 1:8])
4,8<1,5<2,6<3,7

julia> is_join_lattice(p)
false
```
"""
is_join_lattice(P::Poset)=checkl(incidence(P))

"""
`is_meet_lattice(P)`

returns  `true` if `P` is  a meet semilattice, that  is any two elements of
`P` have a unique highest lower bound; returns `false` otherwise.

```julia-repl
julia> p=Poset([i==j || i%4<j%4 for i in 1:8, j in 1:8])
4,8<1,5<2,6<3,7

julia> is_meet_lattice(p)
false
```
"""
is_meet_lattice(P::Poset)=checkl(permutedims(incidence(P)))

end
