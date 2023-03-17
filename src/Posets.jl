"""
This  package deals  with finite  posets. They  always have  internally the
following fields:

  - `elements` The list of elements of the poset.
  - `hasse`:  a list representing  the Hasse diagram  of the poset: the
    `i`-th  entry  is  the  list  of  indices  of elements which cover (are
    immediate  successors of) the  `i`-th element, that  is the list of `j`
    such that `i<j` and there is no `k` such that `i<k<j`.

The following field is computed on demand (for certain algorithms), since
it is quite large for large posets:

  - `incidence`:  a  boolean  matrix  such that `incidence[i,j]==true` iff
    `i<=j` in the poset. This is sometimes called the ζ-matrix of `P`.

There are several ways of defining a poset.

  - `Poset(h::Vector{<:Vector{<:Integer}},elements=1:length(h))`

by entering the Hasse diagram. By default the elements are `1:length(h)`.
```julia-repl
julia> p=Poset([[2,3],[4],[4],Int[]])
1<2,3<4
```
by  default posets are shown as a covering list of maximal chains. Elements
which  are  equivalent  for  the  poset  are  printed together separated by
commas.

```julia-repl
julia> summary(p) # useful for big posets
"Poset{Int64} with 4 elements"

julia> length(p) # the number of elements of the Poset
4

julia> incidence(p) # asking for the incidence matrix
4×4 Matrix{Bool}:
 1  1  1  1
 0  1  0  1
 0  0  1  1
 0  0  0  1
```

A convenient constructor is:
  - `Poset(isless::Function,elements)`

this  constructs a poset from the incidence matrix computed by applying the
`isless`  function  to  each  pair  of  elements. For `isless` one can give
either  a function implementing  `<` or a  function implementing `≤` (it is
`or`-ed with `=` in any case).

```julia-repl
julia> l=vec(collect(Iterators.product(1:2,1:2)))
4-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (2, 1)
 (1, 2)
 (2, 2)

julia> P=Poset((x,y)->all(map(<=,x,y)),l)
(1, 1)<(2, 1),(1, 2)<(2, 2)

julia> eltype(P)
Tuple{Int64, Int64}
```

A poset can be constructed from an incidence matrix
  - `Poset(m::Matrix{Bool},elements=1:size(m,1))`

The last example could also be entered as
```julia-repl
julia> P=Poset([all(map(<=,x,y)) for x in l, y in l],l)
(1, 1)<(2, 1),(1, 2)<(2, 2)
```

Finally  a  poset  can  be  specified  by  a list of tuples specifing cover
relations. The transitive closure of these relations is computed, giving an
incidence  matrix from which the poset is built. The elements of the poset,
if  not specified by a separate list, will be all the elements which appear
in the tuples.

```julia-repl
julia> Poset([(:a,:b),(:c,:d)])
a<b
c<d
```
By  default printing shows the elements.  The printing can be controlled by
`IO`  properties. The package gives a  convenient constructor to build rich
`IO`:
```julia-rep1
julia> rio(io::IO=stdout;p...)=IOContext(io,:limit=>true,p...)

julia> print(rio(indices=true),P) # print indices instead of elements
1<2,3<4
```

ultimate  flexibility  on  printing  is  obtained  by  setting the function
`show_element`  which takes as arguments an  `IO`, the poset, and the index
of the element to print:
```julia-repl
julia> P.show_element=(io,p,n)->join(io,p.elements[n],".");

julia> P
1.1<2.1,1.2<2.2
```

```julia-rep1
julia> print(P) # a form which can be input back in Julia
Poset([[2, 3], [4], [4], Int64[]],[(1, 1), (2, 1), (1, 2), (2, 2)])
```

The  function `≤(p,a,b)` gives  the interval between  `a` and `b` in
the poset `p`.

```julia-repl
julia> <=(P,(2,1),(2,2)) # elements between (2,1) and (2,2)
2-element Vector{Tuple{Int64, Int64}}:
 (2, 1)
 (2, 2)

julia> <=(P,(1,2)) # elements below (1,2)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (1, 2)

julia> >=(P,(1,2)) # elements above (1,2)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (2, 2)

julia> <=(P,(1,2),(2,1)) # nothing in between
Tuple{Int64, Int64}[]
```
see the on-line help on `linear_extension, hasse, incidence, partition, 
covering_chains, transitive_closure, isjoinlattice, ismeetlattice, moebius,
moebiusmatrix,reverse, induced, minimum, maximum` for more information
"""
module Posets 
# this module has only 1 dependency below which could be copied.
using ..Combinat: collectby  # 13 lines
export Poset, linear_extension, hasse, incidence, partition, covering_chains,
transitive_closure, isjoinlattice, ismeetlattice, moebius, moebiusmatrix
export induced #needs using_merge to use with Gapjm

struct Poset{T}
  hasse::Vector{Vector{Int}}
  elements::Vector{T}
  prop::Dict{Symbol,Any}
end

Base.getproperty(o::Poset,s::Symbol)=hasfield(Poset,s) ? getfield(o,s) : 
         getfield(o,:prop)[s]
Base.setproperty!(o::Poset,s::Symbol,v)=getfield(o,:prop)[s]=v
Base.haskey(o::Poset,s::Symbol)=haskey(getfield(o,:prop),s)
Base.get!(f::Function,o::Poset,s::Symbol)=get!(f,getfield(o,:prop),s)

"""
`transitive_closure(M)`
`transitive_closure!(M)`

`M`   should  be   a  square   boolean  matrix   representing  a  relation;
`transitive_closure`  returns a boolean  matrix representing the transitive
closure  of this  relation; `transistive_closure!`  modifies `M`  in place,
doing   no  allocations.  The   transitive  closure  is   computed  by  the
Floyd-Warshall algorithm, which is quite fast even for large matrices.

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
transitive_closure(m)=transitive_closure!(copy(m))

function transitive_closure!(m)
  for k in axes(m,2), i in axes(m,2)
    if m[k,i]
      @views m[:,i].|=m[:,k]
    end
  end
  m
end

"""
`Poset(m::Matrix{Bool},elements=1:size(m,1))`

Creates a poset from an incidence matrix `m`, that is `m[i,j]==true` if and
only if `i≤j` in the poset,

```julia-repl
julia> Poset(Bool[1 1 1 1 1;0 1 0 1 1;0 0 1 1 1;0 0 0 1 0;0 0 0 0 1])
1<2,3<4,5
```
"""
function Poset(m::Matrix{Bool},e=collect(1:size(m,1)))
  p=Poset(hasse(m),e,Dict{Symbol,Any}(:incidence=>m))
  p.show_element=(io,x,n)->print(io,x.elements[n])
  p
end

"""
`Poset(h::Vector{<:Vector{<:Integer}},elements=1:length(h))`

Creates a poset from a Hasse diagram given as a `Vector` whose `i`-th entry
is  the list of indices of elements which are immediate successors (covers)
of the `i`-th element, that is `h[i]` is the list of `j` such that `i<j` in
the poset and such that there is no `k` such that `i<k<j`.

```julia-repl
julia> Poset([[2,3],[4,5],[4,5],Int[],Int[]])
1<2,3<4,5
```
"""
function Poset(h::Vector{<:Vector{<:Integer}},e=collect(1:length(h)))
  p=Poset(h,e,Dict{Symbol,Any}())
  p.show_element=(io,x,n)->print(io,x.elements[n])
  p
end

Poset(f::Function,e)=Poset([(f(x,y)|| x==y) for x in e, y in e],collect(e))

function Poset(h::Vector{Tuple{T,T}},e::AbstractVector{T})where T
  inc=zeros(Bool,length(e),length(e))
  for i in eachindex(e) inc[i,i]=true end
  for (i,j) in h inc[findfirst(==(i),e),findfirst(==(j),e)]=true end
  transitive_closure!(inc)
  Poset(inc,collect(e))
end
  
Poset(h::Vector{Tuple{T,T}}) where T=Poset(h,unique(vcat(collect.(h)...)))

Base.length(p::Poset)::Int=length(elements(p))
Base.eltype(p::Poset{T}) where T=T

elements(p::Poset)=p.elements

function Base.show(io::IO,x::Poset)
  s=hasse(x)
  if !(get(io,:TeX,false) || get(io,:limit,false))
    print(io,"Poset(",s)
    if x.elements!=1:length(x) print(io,",",x.elements) end
    print(io,")")
    return 
  end
  pp=partition(x)
  sh=get(io,:indices,false) ? (io,x,n)->print(io,n) : x.show_element
  e=elements(x)
  labels=map(y->join(map(n->sprint(sh,x,n;context=io),y),","),pp)
  p=Poset(map(x->unique!(sort(map(y->Int(findfirst(z->y in z,pp)),s[x[1]]))),pp))
  sep=get(io,:Symbol,false)
  TeX=get(io,:TeX,false)
  if sep==false sep=TeX ? "{<}" : "<" end
  ch=map(x->join(labels[x],sep), covering_chains(p)) 
  if TeX print(io,join(map(x->"\$\$$x\$\$\n",ch)))
  else join(io,ch,"\n")
  end
end

Base.summary(io::IO,p::Poset)=print(io,typeof(p)," with ",length(p)," elements")

"""
`linear_extension(P)`

returns  a  linear  extension  of  the  poset  `P`,  that  is  a vector `l`
containing  a permutation of the integers  `1:length(P)` such that if `i<j`
in `P`, then `findfirst(==(i),l)<findfirst(==(j),l)`. This is also called a
topological sort of `P`.

```julia-repl
julia> p=Poset((i,j)->j%i==0,1:6)
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
function linear_extension(P::Poset)::Vector{Int}
  get!(P,:linear_extension)do
    n=zeros(length(P))
    for v in hasse(P), x in v n[x]+=1 end
    Q=filter(x->iszero(n[x]),1:length(n))
    res=Int[]
    while !isempty(Q)
      i=popfirst!(Q)
      push!(res, i)
      for x in hasse(P)[i]
        n[x]-=1
        if iszero(n[x]) push!(Q, x) end
      end
    end
    if !iszero(n) error("cycle") end
    res
   end::Vector{Int}
end

hasse(m::Matrix{Bool})=map(x->filter(y->x[y]==2,1:length(x)),eachrow(m*m))

"""
`hasse(P::Poset)`

the Hasse diagram of `P`.

```julia-repl
julia> p=Poset((i,j)->j%i==0,1:5)
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
hasse(p::Poset)=p.hasse

"""
`incidence(P::Poset)`

returns the incidence matrix (also called the ζ matrix) of `P`.

```julia-repl
julia> p=Poset([i==6 ? Int[] : [i+1] for i in 1:6])
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
function incidence(p::Poset)
  get!(p,:incidence)do
    inc=zeros(Bool,length(p),length(p))
    for (i,s) in enumerate(hasse(p))
      inc[i,i]=true
      for j in s inc[i,j]=true end
    end
    transitive_closure!(inc)
  end::Matrix{Bool}
end

coxetermatrix(p::Poset)=-moebiusmatrix(p)*transpose(incidence(p))

"""
`covering_chains(P::Poset)`

A (greedy: the first is longest possible) list of covering chains for P.
"""
function covering_chains(P::Poset)
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
`reverse(P::Poset)`

the opposed poset to `P`.

```julia-repl
julia> p=Poset((i,j)->i%4<j%4,1:8)
4,8<1,5<2,6<3,7

julia> reverse(p)
3,7<2,6<1,5<4,8
```
"""
function Base.reverse(p::Poset)
  h=hasse(p)
  resh=map(empty,h)
  for i in 1:length(p), j in h[i] push!(resh[j], i) end
  res=Poset(resh,p.elements,copy(p.prop))
  if haskey(p,:incidence) res.incidence=permutedims(incidence(p)) end
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
  l=hasse(reverse(p))
  return collectby(i->(l[i], hasse(p)[i]),1:length(p))
  if false
    I=incidence(p)
    n=.!(one(I))
    collectby(map(i->(I[i,:].&n[i,:],I[:,i].&n[i,:]),1:length(p)),1:length(p))
  end
end

"""
`induced(P::Poset,S)`

returns  the subposet  induced by  `P` on  `S`, a sublist of `elements(P)`.
Note  that  the  sublist  `S`  does  not  have  to  be in the same order as
`elements(P)`, so this can be just used to renumber the elements of `P`.

```julia-repl
julia> p=Poset((i,j)->i%4<j%4,1:8)
4,8<1,5<2,6<3,7

julia> induced(p,2:6)
4<5<2,6<3
```
"""
function induced(p::Poset,S::AbstractVector)
  ind=indexin(S,elements(p))
  if length(ind)==length(p) && sort(ind)==1:length(p)
    resh=Vector{Int}.(map(x->map(y->findfirst(==(y),ind),x),hasse(p)[ind]))
    res=Poset(resh,p.elements[ind],copy(p.prop))
    if haskey(res, :incidence) res.incidence=incidence(res)[ind,ind] end
  else
    inc=incidence(p)
    inc=[i!=j && ind[i]==ind[j] ? false : inc[ind[i],ind[j]]
               for i in eachindex(ind), j in eachindex(ind)]
    res=Poset(hasse(inc),collect(S),copy(p.prop))
    res.incidence=inc
  end
  res
end

function checkl(ord::AbstractMatrix{Bool})
  subl=Set(Vector{Bool}[])
  n=size(ord,1)
  for i in 1:n, j in 1:i-1
    if !ord[j,i] || ord[i,j]
      @views l=ord[i,:].&ord[j,:]
      if !(l in subl)
        if !any(y->l==l.&y,eachrow(ord[l,:]))
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
`isjoinlattice(P::Poset)`

returns  `true` if `P` is  a join semilattice, that  is any two elements of
`P` have a unique smallest upper bound; returns `false` otherwise.

```julia-repl
julia> p=Poset((i,j)->j%i==0,1:8)
1<5,7
1<2<4<8
1<3<6
2<6

julia> isjoinlattice(p)
false
```
"""
isjoinlattice(P::Poset)=checkl(incidence(P))

"""
`ismeetlattice(P)`

returns  `true` if `P` is  a meet semilattice, that  is any two elements of
`P` have a unique highest lower bound; returns `false` otherwise.

```julia-repl
julia> p=Poset((i,j)->j%i==0,1:8)
1<5,7
1<2<4<8
1<3<6
2<6

julia> ismeetlattice(p)
true
```
"""
ismeetlattice(P::Poset)=checkl(transpose(incidence(P)))

"""
`moebius(P,y=maximum(P))`

the vector of values `μ(x,y)` of the Moebius function of `P` for `x` varying.
Here is an example giving the ususal Moebius function on integers.
```julia_repl
julia> p=Poset((i,j)->i%j==0,1:8)
5,7<1
6<2<1
6<3<1
8<4<2

julia> moebius(p)
8-element Vector{Int64}:
  1
 -1
 -1
  0
 -1
  1
 -1
  0
```
"""
function moebius(P::Poset,y=0)
  o=linear_extension(P)
  if y==0 y=length(o)
  else y=findfirst(==(y),o)
  end
  mu=zeros(Int,length(P))
  mu[o[y]]=1
  I=incidence(P)
  for i in y-1:-1:1 
    mu[o[i]]=0
    for j in i+1:y
      if I[o[i],o[j]]
        mu[o[i]]-=mu[o[j]]
      end
    end
  end
  mu
end

function unitriangularinv(b::Matrix)
  a=Int.(b)
  for k in axes(b,1), i in k+1:size(b,1) 
    a[k,i]=-b[k,i]-sum(j->b[j,i]*a[k,j],k+1:i-1;init=0)
  end
  a
end

"`moebiusmatrix(P::Poset)` the matrix of the Moebius function `μ(x,y)`"
function moebiusmatrix(P::Poset)
  o=linear_extension(P)
  r=unitriangularinv(incidence(P)[o,o])
  o=invperm(o)
  r[o,o]
end

function Base.argmax(p::Poset)
  m=incidence(p)
  maxs=findall(i->all(@view m[:,i]),axes(m,2))
  if length(maxs)==1 return only(maxs) end
end

"`maximum(P::Poset)` the maximum element of `P` if it has one, otherwise `nothing`"
function Base.maximum(p::Poset)
  m=argmax(p)
  if !isnothing(m) return elements(p)[m] end
end

function Base.argmin(p::Poset)
  m=incidence(p)
  mins=findall(i->all(@view m[i,:]),axes(m,1))
  if length(mins)==1 return only(mins) end
end

"`minimum(P::Poset)` the minimum element of `P` if it has one, otherwise `nothing`"
function Base.minimum(p::Poset)
  m=argmin(p)
  if !isnothing(m) return elements(p)[m] end
end

index(P,a)=findfirst(==(a),elements(P))
Base.:≤(p::Poset,a)=p.elements[incidence(p)[:,index(p,a)]]
Base.:≤(p::Poset,a,b)=p.elements[incidence(p)[index(p,a),:].&incidence(p)[:,index(p,b)]]
Base.:≥(p::Poset,a)=p.elements[incidence(p)[index(p,a),:]]

function Base.getindex(p::Poset,r::UnitRange)
  I=1:length(p)
  if r.start==0 
    if r.stop==1+length(p) I
    else I[incidence(p)[:,r.stop]]
    end
  else 
    if r.stop==1+length(p) I[incidence(p)[r.start,:]]
    else I[incidence(p)[r.start,:].&incidence(p)[:,r.stop]]
    end
  end
end

Base.lastindex(p::Poset)=1+length(p)
Base.firstindex(p::Poset)=0

end
