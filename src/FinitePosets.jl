"""
This  package deals with  finite posets. 

There  are two kinds of  posets. A "canonical poset"  or `CPoset` is on the
elements  `1:n`  where  `n=length(P)`.  A  `Poset`  is  on  a given list of
elements  which may be of any type. A `Poset` contains internally a`CPoset`
working  on  the  indices  of  the  elements,  which is more efficient than
working  with  the  elements  themselves.  Similarly  for  efficiency  many
functions  work with the  `CPoset` by transforming  their input to indices
and their output to elements.

A  `CPoset` has the field:

  - `hasse`:  a list representing  the Hasse diagram  of the poset: the `i`-th  entry is the list of elements which cover (are immediate  successors of) `i`, that  is the list of `j` such that `i<j` and there is no `k` such that `i<k<j`.

The following is cached when computed to speed up subsequent computations:

  - `incidence`: a  boolean matrix  such that `incidence[i,j]==true` iff `i<=j`. This is sometimes called the ζ-matrix of the poset.

There are several ways of defining a poset.  By entering the Hasse diagram:
```julia-repl
julia> p=CPoset([[2,3],[4],[4],Int[]])
1<2,3<4
```
As  seen above, `p` is shown as a list of covering maximal chains; elements
which  are  equivalent  for  the  poset  are  printed together separated by
commas.

```julia-repl
julia> length(p) # the number of elements of `p`
4

julia> incidence(p)
4×4 Matrix{Bool}:
 1  1  1  1
 0  1  0  1
 0  0  1  1
 0  0  0  1

julia> linear_extension(p) # a total order compatible with p
4-element Vector{Int64}:
 1
 2
 3
 4
```

A `Poset` is constructed from a `CPoset` and a list of elements

```julia-repl
julia> P=Poset(p,[:a,:b,:c,:d])
a<b,c<d

julia> P.C # the CPoset attached to P
1<2,3<4
```

A  convenient  constructor  for  `Poset`s  takes  a  function  representing
`isless`  for the poset and  the list of elements  and constructs the poset
from  the incidence matrix, computed by  applying the function to each pair
of  elements. For `isless` one can  give either a function implementing `<`
or a function implementing `≤` (it is `or`-ed with `==` in any case).
```julia-repl
julia> l=vec(collect(Iterators.product(1:2,1:2)))
4-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (2, 1)
 (1, 2)
 (2, 2)

julia> P=Poset((x,y)->all(map(<=,x,y)),l)
(1, 1)<(2, 1),(1, 2)<(2, 2)

julia> eltype(P) # the type of the elements of P
Tuple{Int64, Int64}

julia> summary(P) # useful for big posets
"Poset{Tuple{Int64, Int64}} of length 4"
```
A  poset  can  also  be  constructed  from  an incidence matrix so the last
example could also be entered as
```julia-repl
julia> P=Poset(CPoset([all(map(<=,x,y)) for x in l, y in l]),l)
(1, 1)<(2, 1),(1, 2)<(2, 2)
```
Flexibility  on  printing  a  `Poset`  is  obtained by setting the function
`show_element`  which takes as arguments an  `IO`, the poset, and the index
of the element to print:
```julia-repl
julia> P.show_element=(io,p,n)->join(io,p.elements[n],".");

julia> P
1.1<2.1,1.2<2.2

julia> delete!(P,:show_element); # back to default
```
The above fancy printing applies only when printing at the REPL or in pluto
or  Jupyter. The default printing  gives a form which  can be input back in
Julia
```julia-rep1
julia> print(P) 
Poset(CPoset([[2, 3], [4], [4], Int64[]]),[(1, 1), (2, 1), (1, 2), (2, 2)])
```
A poset can be specified by a list of tuples specifing order relations. The
transitive  closure  of  these  relations  is computed, giving an incidence
matrix  from which the  poset is built.  The elements of  the poset, if not
specified separately, will be all the elements which appear in the tuples.
```julia-repl
julia> Poset([(:a,:b),(:c,:d)])
a<b
c<d

julia> CPoset([(1,3),(2,5)]) # the CPoset is on 1:maximum(entries)
4
1<3
2<5
```
To get the order relation `≤` of the poset `p` between elements
`i` and `j` just call `≤(p,i,j)`. 
```julia-repl
julia> ≤(P,(1,1),(2,1))
true

julia> ≤(P.C,1,2) # the same
true
```

Intervals in a poset can be computed with strict or not bounds.
```julia-repl
julia> interval(P,≤,(1,2)) # elements below (1,2)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (1, 2)

julia> interval(P,≥,(1,2)) # elements above (1,2)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (2, 2)

julia> interval(P,<,(1,2)) # elements strictly below (1,2)
1-element Vector{Tuple{Int64, Int64}}:
 (1, 1)

julia> interval(P,≥,(2,1),≤,(2,2)) # elements between (2,1) and (2,2)
2-element Vector{Tuple{Int64, Int64}}:
 (2, 1)
 (2, 2)

julia> interval(P,>,(1,1),<,(2,2)) # elements strictly between
2-element Vector{Tuple{Int64, Int64}}:
 (2, 1)
 (1, 2)
julia> interval(P.C,>,1,<,4) # in terms of indices
2-element Vector{Int64}:
 2
 3
```
A sample of other functions available on posets:
```julia-repl
julia> maximal_chains(P)
2-element Vector{Vector{Tuple{Int64, Int64}}}:
 [(1, 1), (2, 1), (2, 2)]
 [(1, 1), (1, 2), (2, 2)]

julia> height(P) # the length of a maximal chain
3

julia> moebiusmatrix(P)
4×4 Matrix{Int64}:
 1  -1  -1   1
 0   1   0  -1
 0   0   1  -1
 0   0   0   1

julia> minima(P)
1-element Vector{Tuple{Int64, Int64}}:
 (1, 1)

julia> maxima(P)
1-element Vector{Tuple{Int64, Int64}}:
 (2, 2)

julia> Q=CPoset(:chain,3)
1<2<3

julia> P1=Poset(Q) # transformed to a Poset with elements 1:3
1<2<3

julia> P⊕ P1 # the ordinal sum
(1, 1)<(2, 1),(1, 2)<(2, 2)<1<2<3

julia> P1*P1
(1, 1)<(2, 1)<(3, 1)<(3, 2)<(3, 3)
(1, 1)<(1, 2)<(2, 2)<(3, 2)
(2, 1)<(2, 2)
(1, 2)<(1, 3)<(2, 3)<(3, 3)
(2, 2)<(2, 3)

julia> P1⊗ P1 # the ordinal product
(1, 1)<(1, 2)<(1, 3)<(2, 1)<(2, 2)<(2, 3)<(3, 1)<(3, 2)<(3, 3)
```

Finally `showpic(p)` where `p` is a `CPoset` or a `Poset` gives a graphical
display of the poset (on Linux) provided you have the command `dot` and the
command  `display` of  `imagemagick` installed.  It may  work on  MacOs and
Windows but I did not test it.

see the on-line help on
`⊕,
⊗, 
+,
*,
chains,
chainpoly,
covering_chains,
coxetermatrix, 
dual, 
hasse, 
height,
incidence, 
induced, 
interval, 
isjoinlattice, 
ismeetlattice, 
linear_extension, 
maxima,
maximal_chains, 
minima,
moebius, 
moebiusmatrix, 
partition, 
showpic,
transitive_closure`
for more information
"""
module FinitePosets 
# this module has only one dependency.
using Combinat: collectby, combinations, tally, partitions, dominates 
export CPoset, Poset,
⊕,
⊗, 
chains,
chainpoly,
covering_chains,
coxetermatrix, 
dual, 
hasse, 
height,
incidence, 
induced, 
interval, 
isjoinlattice, 
ismeetlattice, 
linear_extension, 
maxima,
maximal_chains, 
minima,
moebius, 
moebiusmatrix, 
partition, 
showpic,
transitive_closure

"""
`transitive_closure(M)`
`transitive_closure!(M)`

`M`   should  be   a  square   boolean  matrix   representing  a  relation;
`transitive_closure`  returns a boolean  matrix representing the transitive
closure  of  this  relation;  `transitive_closure!`  modifies `M` in place,
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
`hasse(m::Matrix{Bool})`
Given  an  incidence  matrix  for  a  poset returns the corresponding Hasse
diagram.
```julia-repl
julia> m=incidence(CPoset(:diamond,5))
5×5 Matrix{Bool}:
 1  1  1  1  1
 0  1  0  0  1
 0  0  1  0  1
 0  0  0  1  1
 0  0  0  0  1

julia> hasse(m)
5-element Vector{Vector{Int64}}:
 [2, 3, 4]
 [5]
 [5]
 [5]
 []
```
"""
hasse(m::Matrix{Bool})=map(x->filter(y->x[y]==2,1:length(x)),eachrow(m*m))

abstract type AbstractPoset{T} end

struct CPoset<:AbstractPoset{Int}
  hasse::Vector{Vector{Int}}
  prop::Dict{Symbol,Any}
end

Base.getproperty(o::CPoset,s::Symbol)=hasfield(CPoset,s) ? getfield(o,s) : 
         getfield(o,:prop)[s]
Base.setproperty!(o::CPoset,s::Symbol,v)=getfield(o,:prop)[s]=v
Base.haskey(o::CPoset,s::Symbol)=haskey(getfield(o,:prop),s)
Base.get!(f::Function,o::CPoset,s::Symbol)=get!(f,getfield(o,:prop),s)
Base.delete!(p::CPoset,s::Symbol)=delete!(p.prop,s)

Base.length(p::AbstractPoset)=length(hasse(p))

Base.eltype(p::AbstractPoset{T}) where T=T

struct Poset{T}<:AbstractPoset{T}
  C::CPoset
  elements::Vector{T}
  prop::Dict{Symbol,Any}
end

Base.getproperty(o::Poset,s::Symbol)=hasfield(Poset,s) ? getfield(o,s) : 
         getfield(o,:prop)[s]
Base.setproperty!(o::Poset,s::Symbol,v)=getfield(o,:prop)[s]=v
Base.haskey(o::Poset,s::Symbol)=haskey(getfield(o,:prop),s)
Base.get!(f::Function,o::Poset,s::Symbol)=get!(f,getfield(o,:prop),s)
Base.delete!(p::Poset,s::Symbol)=delete!(p.prop,s)

"""
`Poset(p::CPoset,e::AbstractVector=1:length(p))`

creates a `Poset` with order specified by `p` and elements `e`.
```julia-repl
julia> Poset(CPoset([[2,3],[4],[4],Int[]]),[:a,:b,:c,:d])
a<b,c<d
```
with no second argument transforms a `CPoset` into a `Poset`.
"""
Poset(p::CPoset,e::AbstractVector{T}) where T=Poset(p,e isa Vector ? e : collect(e),Dict{Symbol,Any}())
Poset(P::CPoset)=Poset(P,1:length(P))
CPoset(P::Poset)=P.C

hasse(P::Poset)=hasse(P.C)
linear_extension(P::Poset)=linear_extension(P.C)
incidence(P::Poset)=incidence(P.C)
partition(P::Poset)=partition(P.C)
moebiusmatrix(P::Poset)=moebiusmatrix(P.C)
coxetermatrix(P::Poset)=coxetermatrix(P.C)
chainpoly(P::Poset)=chainpoly(P.C)
height(P::Poset)=height(P.C)
isjoinlattice(P::Poset)=isjoinlattice(P.C)
ismeetlattice(P::Poset)=ismeetlattice(P.C)
"""
`CPoset(m::Matrix{Bool})`

Creates a poset from an incidence matrix `m`, that is `m[i,j]==true` if and
only if `i≤j` in the poset,

```julia-repl
julia> CPoset(Bool[1 1 1 1 1;0 1 0 1 1;0 0 1 1 1;0 0 0 1 0;0 0 0 0 1])
1<2,3<4,5
```
"""
CPoset(m::Matrix{Bool})=CPoset(hasse(m),Dict{Symbol,Any}(:incidence=>m))

"""
`CPoset(h::Vector{<:Vector{<:Integer}})`

Creates a poset from a Hasse diagram given as a `Vector` whose `i`-th entry
is  the  list  of  indices  which  are immediate successors (covers) of the
`i`-th  element, that is `h[i]`  is the list of  `j` such that `i<j` in the
poset and such that there is no `k` such that `i<k<j`.
```julia-repl
julia> CPoset([[2,3],[4,5],[4,5],Int[],Int[]])
1<2,3<4,5
```
"""
CPoset(h::Vector{<:Vector{<:Integer}})=CPoset(h,Dict{Symbol,Any}())

"""
`Poset(f::Function,e::AbstractVector)`

creates a `Poset` with elements `e` and order between two elements given by
function `f`.
```julia-repl
julia> Poset((x,y)->all(x.≤y),vec(collect(Iterators.product(1:2,1:3))))
(1, 1)<(2, 1)<(2, 2)<(2, 3)
(1, 1)<(1, 2)<(2, 2)
(1, 2)<(1, 3)<(2, 3)
```
"""
function Poset(f::Function,e::AbstractVector)
  Poset(CPoset([(f(x,y)|| x==y) for x in e, y in e]),collect(e))
end
"""
`CPoset(f::Function,n::integer)`

creates the `Poset` on `1:n` with order given by function `f`.
```julia-repl
julia> CPoset((x,y)->y%x==0,8)  # the divisibility poset
1<5,7
1<2<4<8
1<3<6
2<6
```
"""
CPoset(f::Function,n::Integer)=CPoset([(f(x,y)|| x==y) for x in 1:n, y in 1:n])

"""
`CPoset(covers::Vector{Tuple{Int,Int}})`

creates a poset representing the transitive closure of the given relations.
The  poset is on `1:n` where `n` is the maximum number which appears in the
relations.
```julia-repl
julia> CPoset([(6,2),(5,1)])
3,4
5<1
6<2
```
"""
function CPoset(covers::Vector{Tuple{Int,Int}})
  n=maximum(maximum.(covers))
  inc=zeros(Bool,n,n)
  for i in 1:n inc[i,i]=true end
  for (i,j) in covers inc[i,j]=true end
  transitive_closure!(inc)
  CPoset(inc);
end

"""
`Poset(covers::Vector{Tuple{T,T}}) where T`

creates a poset representing the transitive closure of the given relations.
The  poset is on the elements which appear in the relations.
```julia-repl
julia> Poset([(:a,:b),(:d,:c)])
a<b
d<c
```
"""
function Poset(covers::Vector{Tuple{T,T}})where T
  e=sort(unique(collect(Iterators.flatten(covers))))
  P=CPoset(map(x->map(i->findfirst(==(i),e),x),covers))
  Poset(P,e)
end
  
"""
  - `CPoset(:chain,n)`  a chain on `1:n`
  - `CPoset(:antichain,n)`  an antichain on `1:n`
"""
CPoset(s::Symbol,arg...)=CPoset(Val(s),arg...)
Poset(s::Symbol,arg...)=Poset(Val(s),arg...)
CPoset(::Val{:chain},n::Int)=CPoset(push!(map(i->[i],2:n),Int[]))
"""
  - `Poset(:chain,e)`  a chain with elements `e`
  - `Poset(:antichain,e)`  an antichain with elements `e`
  - `Poset(:powerset,n::Integer)`  the powerset of the set `1:n` with inclusion
```julia-repl
julia> p=Poset(:powerset,3);p.show_element=(io,p,n)->join(io,p.elements[n]);

julia> p
<1<12<123
<2<12
<3<13<123
1<13
2<23<123
3<23
```
  - `Poset(:powerset,e)`  the powerset of the set `e` with inclusion
  - `Poset(:partitionsdominance,n)`  the poset of partitions of `n` with dominance order
```julia-repl
julia> Poset(:partitionsdominance,5)
[1, 1, 1, 1, 1]<[2, 1, 1, 1]<[2, 2, 1]<[3, 1, 1]<[3, 2]<[4, 1]<[5]
```
"""
Poset(::Val{:chain},v::AbstractVector)=Poset(CPoset(:chain,length(v)),v)
CPoset(::Val{:antichain},n::Int)=CPoset([Int[] for i in 1:n])
Poset(::Val{:antichain},v::AbstractVector)=Poset(CPoset(:antichain,length(v)),v)
Poset(::Val{:powerset},n::Int)=Poset(issubset,combinations(1:n))
Poset(::Val{:powerset},v::AbstractVector)=Poset(issubset,combinations(v))
Poset(::Val{:partitionsdominance},n::Int)=dual(Poset(dominates,partitions(n)))
function CPoset(::Val{:diamond},n::Int)
  h=[[n] for i in 1:n]; h[n]=Int[]; h[1]=2:n-1
  CPoset(h)
end

function Base.show(io::IO,x::AbstractPoset)
  s=hasse(x)
  if !(get(io,:TeX,false) || get(io,:limit,false))
    if x isa CPoset print(io,"CPoset(",s,")")
    else print(io,"Poset(",CPoset(x),",",x.elements,")")
    end
    return 
  end
  if haskey(x,:name) println(io,x.name) end
  pp=partition(x)
  if x isa Poset
    sh=haskey(x,:show_element) ? x.show_element : (io,x,n)->print(io,x.elements[n])
    labels=map(y->join(map(n->sprint(sh,x,n;context=io),y),","),pp)
  else
    labels=join.(pp,",")
  end
  p=CPoset(map(x->unique!(sort(map(y->Int(findfirst(z->y in z,pp)),s[x[1]]))),pp))
  sep=get(io,:Symbol,false)
  TeX=get(io,:TeX,false)
  if sep==false sep=TeX ? "{<}" : "<" end
  ch=map(x->join(labels[x],sep), covering_chains(p)) 
  if TeX print(io,join(map(x->"\$\$$x\$\$\n",ch)))
  else join(io,ch,"\n")
  end
end

Base.summary(io::IO,p::AbstractPoset)=print(io,typeof(p)," of length ",length(p))

"""
`linear_extension(P::CPoset)`

returns a linear extension of the `CPoset`, that is a vector `l` containing
a permutation of the integers `1:length(P)` such that if `i<j` in `P` (that
is  `incidence(P)[i,j]` is `true`), then `i` is  before `j` in `l`. This is
also called a topological sort of `P`.

```julia-repl
julia> p=CPoset((i,j)->j%i==0,6) # divisibility poset on 1:6
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
`linear_extension(P::Poset)` returns a linear extension of `P.C`.
"""
function linear_extension(P::CPoset)::Vector{Int}
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

"""
`hasse(P::CPoset)`

the Hasse diagram of `P`.

```julia-repl
julia> p=CPoset((i,j)->j%i==0,5)
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
`hasse(P::Poset)` returns `hasse(P.C)`.
"""
hasse(p::CPoset)=p.hasse

"""
`incidence(P::CPoset)`

returns the incidence matrix (also called the ζ matrix) of `P`.

```julia-repl
julia> p=CPoset([i==6 ? Int[] : [i+1] for i in 1:6])
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
`incidence(P::Poset)` returns `incidence(P.C)`.
"""
function incidence(p::CPoset)
  get!(p,:incidence)do
    inc=zeros(Bool,length(p),length(p))
    for (i,s) in enumerate(hasse(p))
      inc[i,i]=true
      for j in s inc[i,j]=true end
    end
    transitive_closure!(inc)
  end::Matrix{Bool}
end

"""
`P+Q` returns the sum of two `CPoset`s or of two `Poset`s.
```julia-repl
julia> CPoset(:chain,2)+CPoset(:chain,3)
1<2
3<4<5

julia> Poset(:chain,[1,2])+Poset(:chain,[:a,:b,:c])
1<2
a<b<c
```
"""
Base.:+(P::CPoset,Q::CPoset)=CPoset(vcat(hasse(P),map(x->x.+length(P),hasse(Q))))
Base.:+(P::Poset,Q::Poset)=Poset(P.C+Q.C,vcat(P.elements,Q.elements))

"""
`P⊕ Q` returns the ordinal sum of two `CPoset`s or of two `Poset`s.
```julia-repl
julia> CPoset(:chain,2)⊕ CPoset(:chain,3)
1<2<3<4<5

julia> Poset(:chain,[1,2])⊕ Poset(:chain,[:a,:b,:c])
1<2<a<b<c
```
"""
function ⊕(P::CPoset,Q::CPoset)
  m=cat(incidence(P),incidence(Q),dims=(1,2))
  m[1:length(P),length(P).+(1:length(Q))].=true
  CPoset(m)
end
⊕(P::Poset,Q::Poset)=Poset(P.C⊕Q.C,vcat(P.elements,Q.elements))

"""
`P*Q` returns the product of two `CPoset`s or of two `Poset`s.
```julia-repl
julia> CPoset(:chain,2)*CPoset(:chain,3)
1<2<3<6
1<4<5<6
2<5

julia> Poset(:chain,[1,2])*Poset(:chain,[:a,:b,:c])
(1, :a)<(2, :a)<(1, :b)<(2, :c)
(1, :a)<(2, :b)<(1, :c)<(2, :c)
(2, :a)<(1, :c)
```
"""
Base.:*(P::CPoset,Q::CPoset)=CPoset(kron(incidence(P),incidence(Q)))
Base.:*(P::Poset,Q::Poset)=Poset(P.C*Q.C,
                vec(collect(Iterators.product(P.elements,Q.elements))))

covers(P::CPoset)=vcat([map(j->(i,j),s) for (i,s) in enumerate(hasse(P))]...)

"""
`dot(p)` gives a rendering of the Hasse diagram of the
`Poset` or `CPoset` in the graphical language `dot`.
"""
function dot(P::AbstractPoset)
  if P isa CPoset q=string
  else  io=IOContext(stdout,:limit=>true)
    q(i)=haskey(P,:show_element) ? 
    sprint(P.show_element,P,i;context=io) :
    sprint(show,P.elements[i];context=io)
  end
  res="digraph {\n"
  for i in 1:length(P) res*="\""*q(i)*"\";" end
  res*="\n"
  c= P isa CPoset ? covers(P) : covers(P.C)
  for (i,j) in c res*="\""*q(j)*"\"->\""*q(i)*"\"[dir=back];" end
  res*"\n}"
end

"""
`showpic(p)` display a graphical representation of the Hasse diagram of the
`Poset` or `CPoset` using the commands `dot` and `display`.
"""
function showpic(P::AbstractPoset)
  open(pipeline(`dot -Tpng`,`display`),"w")do f
    print(f,dot(P))
  end
end

"""
`P⊗ Q` returns the ordinal product of two `CPoset`s or of two `Poset`s.
```julia-repl
julia> CPoset(:chain,2)⊗ CPoset(:chain,3)
1<3<5<2<4<6

julia> Poset(:chain,[1,2])⊗ Poset(:chain,[:a,:b,:c])
(1, :a)<(1, :b)<(1, :c)<(2, :a)<(2, :b)<(2, :c)
```
"""
function ⊗(P::CPoset,Q::CPoset)
  a=covers(P)
  b=covers(Q)
  res=Tuple{Int,Int}[]
  e=LinearIndices((length(P),length(Q)))
  for i in 1:length(P)
    for (j,k) in covers(Q) push!(res,(e[i,j],e[i,k])) end
  end
  for (i,j) in covers(P)
    for k in 1:length(Q), l in 1:length(Q) push!(res,(e[i,k],e[j,l])) end
  end
  CPoset(res)
end
⊗(P::Poset,Q::Poset)=Poset(P.C⊗ Q.C,
              vec(collect(Iterators.product(P.elements,Q.elements))))

"""
`coxetermatrix(p)` the Coxeter matrix of the `Poset` or `CPoset`, defined
as `-m*transpose(inv(m))` where `m` is the ζ or incidence matrix.
```julia-repl
julia> coxetermatrix(CPoset(:diamond,5))
5×5 Matrix{Int64}:
  0  -1  -1  -1  -2
  0   0   1   1   1
  0   1   0   1   1
  0   1   1   0   1
 -1  -1  -1  -1  -1
```
"""
coxetermatrix(p::CPoset)=-moebiusmatrix(p)*transpose(incidence(p))

"""
`covering_chains(P::CPoset)`

A (greedy: the first is longest possible) list of covering chains for P.
"""
function covering_chains(P::CPoset)
  ch=Vector{Int}[]
  h=hasse(P)
  for i in linear_extension(P)
    p=findfirst(c->i==c[end],ch)
    if p===nothing
      if isempty(h[i]) push!(ch,[i])
      else for j in h[i] push!(ch,[i,j]) end
      end
    else
      for j in h[i][2:end] push!(ch,[i,j]) end
      if !isempty(h[i]) push!(ch[p],h[i][1]) end
    end
  end
  ch
end

"""
`dual(P)`

the dual poset to the `Poset` or `CPoset` (the order relation is reversed).

```julia-repl
julia> p=CPoset((i,j)->i%4<j%4,8)
4,8<1,5<2,6<3,7

julia> dual(p)
3,7<2,6<1,5<4,8
```
"""
function dual(p::CPoset)
  h=hasse(p)
  resh=map(empty,h)
  for i in 1:length(p), j in h[i] push!(resh[j], i) end
  res=CPoset(resh)
  if haskey(p,:incidence) res.incidence=permutedims(incidence(p)) end
  return res
end

dual(p::Poset)=Poset(dual(p.C),p.elements)

"""
`partition(P::CPoset)`

returns  the partition of `1:length(P)` induced by the equivalence relation
associated  to  `P`;  that  is,  `i`  and  `j`  are in the same part of the
partition  if the `k` such that `i<k` and `j<k` are the same as well as the
`k` such that `k<i` and `k<j`.

```julia-repl
julia> p=CPoset([i==j || i%4<j%4 for i in 1:8, j in 1:8])
4,8<1,5<2,6<3,7

julia> partition(p)
4-element Vector{Vector{Int64}}:
 [4, 8]
 [2, 6]
 [3, 7]
 [1, 5]
```
`partition(P::Poset)` returns `partition(P.C)`
"""
function partition(p::CPoset)
  l=hasse(dual(p))
  return collectby(i->(l[i], hasse(p)[i]),1:length(p))
  if false
    I=incidence(p)
    n=.!(one(I))
    collectby(map(i->(I[i,:].&n[i,:],I[:,i].&n[i,:]),1:length(p)),1:length(p))
  end
end

"""
`induced(P,S)`

returns the subposet induced by `P` on `S`, a sublist of `P.elements` if `P
isa  Poset` or a subset  of `1:length(P)` if `P  isa CPoset`. Note that the
sublist  `S` does not have to be in the same order as `P.elements`, so this
can be just used to renumber the elements of `P`.

```julia-repl
julia> p=CPoset((i,j)->i%4<j%4,8)
4,8<1,5<2,6<3,7

julia> induced(p,2:6) # indices are renumbered
3<4<1,5<2

julia> induced(Poset(p),2:6) # elements are kept
4<5<2,6<3
```
"""
function induced(p::CPoset,ind::AbstractVector{<:Integer})
  if length(ind)==length(p) && sort(ind)==1:length(p)
    resh=Vector{Int}.(map(x->map(y->findfirst(==(y),ind),x),hasse(p)[ind]))
    res=CPoset(resh)
    if haskey(p, :incidence) res.incidence=incidence(p)[ind,ind] end
  else
    inc=incidence(p)
    inc=[i!=j && ind[i]==ind[j] ? false : inc[ind[i],ind[j]]
               for i in eachindex(ind), j in eachindex(ind)]
    res=CPoset(hasse(inc))
    res.incidence=inc
  end
  res
end

function induced(p::Poset{T},S::AbstractVector{T})where T
  ind=Int.(indexin(S,p.elements))
  Poset(induced(p.C,ind),p.elements[ind])
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
`isjoinlattice(P::CPoset)`

returns  `true` if `P` is  a join semilattice, that  is any two elements of
`P` have a unique smallest upper bound; returns `false` otherwise.

```julia-repl
julia> p=CPoset((i,j)->j%i==0,8)
1<5,7
1<2<4<8
1<3<6
2<6

julia> isjoinlattice(p)
false
```
`isjoinlattice(P::Poset)` returns `isjoinlattice(P.C)`
"""
isjoinlattice(P::CPoset)=checkl(incidence(P))

"""
`ismeetlattice(P)`

returns  `true` if `P` is  a meet semilattice, that  is any two elements of
`P` have a unique highest lower bound; returns `false` otherwise.

```julia-repl
julia> p=CPoset((i,j)->j%i==0,8)
1<5,7
1<2<4<8
1<3<6
2<6

julia> ismeetlattice(p)
true
```
`ismeetlattice(P::Poset)` returns `ismeetlattice(P.C)`
"""
ismeetlattice(P::CPoset)=checkl(transpose(incidence(P)))

"""
`moebius(P::CPoset,y=first(maxima(P)))`

the vector of values `μ(x,y)` of the Moebius function of `P` for `x` varying.
Here is an example giving the ususal Moebius function on integers.
```julia_repl
julia> p=CPoset((i,j)->i%j==0,1:8)
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
function moebius(P::CPoset,y=0)
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
moebius(P::Poset,y=0)=moebius(P.C,y)

function unitriangularinv(b::Matrix)
  a=Int.(b)
  for k in axes(b,1), i in k+1:size(b,1) 
    a[k,i]=-b[k,i]-sum(j->b[j,i]*a[k,j],k+1:i-1;init=0)
  end
  a
end

"`moebiusmatrix(P::CPoset)` 
the  matrix  of  the  Moebius  function  `μ(x,y)`  (the inverse of the ζ or
incidence matrix)
```julia-repl
julia> moebiusmatrix(CPoset(:diamond,5))
5×5 Matrix{Int64}:
 1  -1  -1  -1   2
 0   1   0   0  -1
 0   0   1   0  -1
 0   0   0   1  -1
 0   0   0   0   1
```
`moebiusmatrix(P::Poset)` returns `moebiusmatrix(P.C)`
"
function moebiusmatrix(P::CPoset)
  o=linear_extension(P)
  r=unitriangularinv(incidence(P)[o,o])
  o=invperm(o)
  r[o,o]
end

"""
`minima(P)`
the minimal elements of the `Poset` or `CPoset`
```julia-repl
julia> p=CPoset([[3],[3],[4,5],Int[],Int[]])
1,2<3<4,5

julia> minima(p)
2-element Vector{Int64}:
 1
 2
```
"""
function minima(p::CPoset)
  m=incidence(p)
  findall(i->count(@view m[:,i])==1,1:length(p))
end
minima(p::Poset)=p.elements[minima(p.C)]

"""
`maxima(P)`
the maximal elements of the `Poset` or `CPoset`
```julia-repl
julia> p=CPoset([[3],[3],[4,5],Int[],Int[]])
1,2<3<4,5

julia> maxima(p)
2-element Vector{Int64}:
 4
 5
```
"""
function maxima(p::CPoset)
  m=incidence(p)
  findall(i->count(@view m[i,:])==1,1:length(p))
end
maxima(p::Poset)=p.elements[maxima(p.C)]

index(p::Poset{T},x::T) where T=findfirst(==(x),p.elements)
Base.:≤(p::CPoset,a,b)=incidence(p)[a,b]
Base.:≤(p::Poset,a,b)=≤(p.C,index(p,a),index(p,b))

"""
`interval(P,f::Function,a)`
`interval(P,f::Function,a,g::Function,b)`

returns  an interval in the `Poset` or  `CPoset` given by `P`. The function
`f` must be one of the comparison functions `≤, <, ≥, >`. In the first form
it returns the interval between `a` and one end (or the other, depending on
the  comparison function). In  the second form  it returns the
intersection of the intervals `interval(P,f,a)` and `interval(P,g,b)`.
```julia-repl
julia> l=vec(collect(Iterators.product(1:2,1:2)))
4-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (2, 1)
 (1, 2)
 (2, 2)

julia> P=Poset((x,y)->all(map(<=,x,y)),l)
(1, 1)<(2, 1),(1, 2)<(2, 2)

julia> interval(P,≤,(1,2)) # elements below (1,2)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (1, 2)

julia> interval(P,≥,(1,2)) # elements above (1,2)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (2, 2)

julia> interval(P,<,(1,2)) # elements strictly below (1,2)
1-element Vector{Tuple{Int64, Int64}}:
 (1, 1)

julia> interval(P,≥,(2,1),≤,(2,2)) # elements between (2,1) and (2,2)
2-element Vector{Tuple{Int64, Int64}}:
 (2, 1)
 (2, 2)

julia> interval(P,>,(1,1),<,(2,2)) # elements strictly between
2-element Vector{Tuple{Int64, Int64}}:
 (2, 1)
 (1, 2)
julia> interval(P.C,>,1,<,4) # in terms of indices
2-element Vector{Int64}:
 2
 3
```
"""
function interval(p::CPoset,f::Function,i)
  s=Symbol(f)
  if s==:<=  
    (1:length(p))[incidence(p)[:,i]]
  elseif s==:< 
    [j for j in 1:length(p) if incidence(p)[j,i] && i!=j]
  elseif s==:>= 
    (1:length(p))[incidence(p)[i,:]]
  elseif s==:> 
    [j for j in 1:length(p) if incidence(p)[i,j] && i!=j]
  else error(s," unknown")
  end
end
interval(p::Poset,f::Function,a)=p.elements[interval(p.C,f,index(p,a))]

interval(p::CPoset,f::Function,a,g::Function,b)=
  intersect(interval(p,f,a),interval(p,g,b))
interval(p::Poset,f::Function,a,g::Function,b)=
  p.elements[intersect(interval(p.C,f,index(p,a)),interval(p.C,g,index(p,b)))]

"""
`maximal_chains(P)` the maximal chains of the `Poset` or `CPoset`.
```julia-repl
julia> p=Poset([(:a,:b),(:a,:c),(:b,:d),(:c,:d)])
a<b,c<d

julia> maximal_chains(p)
2-element Vector{Vector{Symbol}}:
 [:a, :b, :d]
 [:a, :c, :d]

julia> maximal_chains(p.C)
2-element Vector{Vector{Int64}}:
 [1, 2, 4]
 [1, 3, 4]
```
"""
function maximal_chains(P::CPoset)
  get!(P,:maxchains)do
    p=hasse(P)
    function mc(o)local res, i, f
      f=o[1]
      if isempty(o) return Int[] 
      elseif length(o)==1 return [[f]]
      end
      res=filter(x->!(x[1] in p[f]),mc(o[2:end]))
      if isempty(p[f]) return vcat([[f]],res) end
      for i in p[f] 
        append!(res,map(c->vcat([f],c),
                        filter(x->x[1]==i,mc(o[findfirst(==(i),o):end]))))
      end
      res
    end
    mc(linear_extension(P))
  end
end

maximal_chains(P::Poset)=map(x->P.elements[x],maximal_chains(P.C))

"""
`chains(P)` the chains of the `Poset` or `CPoset`.
```julia-repl
julia> chains(CPoset(:chain,3))
8-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [3]
 [1, 2]
 [1, 3]
 [2, 3]
 [1, 2, 3]
```
"""
function chains(P::CPoset)
  l=maximal_chains(P)
  unique(vcat(combinations.(l)...))
end
chains(P::Poset)=map(x->P.elements[x],chains(P.C))

"""
`chainspoly(P)` the chain polynomial of the `Poset` or `CPoset`, returned
as the list of its coefficients.
```julia-repl
julia> chainpoly(Poset(:powerset,3))
5-element Vector{Int64}:
  1
  8
 19
 18
  6
```
"""
function chainpoly(P::CPoset)
  l=tally(length.(chains(P)))
  v=fill(0,l[end][1]+1)
  for (l,m) in l v[l+1]=m end
  v
end

"""
`height(P)` the height the `Poset` or `CPoset` (the longest length of a chain).
"""
height(P::CPoset)=maximum(length.(maximal_chains(P)))
end
