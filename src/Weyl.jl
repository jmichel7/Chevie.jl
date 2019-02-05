"""
*Weyl groups* are the finite real reflection groups (finite Coxeter groups)
which  can be defined  over the rational  numbers (that is, finite rational
reflection  groups). We implement them as a subgroup of the permutations of
a  root system. Root systems play an  important role in mathematics as they
classify semi-simple Lie algebras and algebraic groups.

Let  us give precise definitions. Let `V`  be a real vector space, `Vⱽ` its
dual  and let `(,)`  be the natural  pairing between `Vⱽ`  and `V`. A *root
system*  is a finite set  of vectors `R` which  generate `V` (the *roots*),
together  with  a  map  `r↦  rⱽ`  from  `R`  to  a subset `Rⱽ` of `Vⱽ` (the
*coroots*) such that:

-  For any `r∈  R`, we have  `(rⱽ,r)=2` so that  the formula `x↦ x-(rⱽ,x)r`
defines a reflection `s_r:V→ V` with root `r` and coroot `rⱽ`.
- The reflection `s_r` stabilizes `R`.

We  will only  consider *reduced*  root systems,  i.e., such  that the only
elements of `R` colinear with `r∈ R` are `r` and `-r`; we also ask that the
root  system be *crystallographic*, that is `(rⱽ,s)` is an integer, for any
`s∈ R,rⱽ∈ Rⱽ`.

The  subgroup `W=W(R)` of  `GL(V)` generated by  the reflections `s_r` is a
Weyl  group --- the representation `V` of  `W` is defined over the rational
numbers.  All other finite-dimensional (complex)  representations of a Weyl
group  can also be realized  over the rational numbers.  Weyl groups can be
characterized  amongst finite Coxeter  groups by the  fact that all numbers
`m(s,t)` in the Coxeter matrix are in `{2,3,4,6}`.

If  we identify  `V` with  `Vⱽ` by  choosing a  `W`-invariant bilinear form
`(.;.)`;  then we have `rⱽ=2r/(r;r)`. A root system `R` is *irreducible* if
it is not the union of two orthogonal subsets. If `R` is reducible then the
corresponding  Weyl  group  is  the  direct  product  of  the  Weyl  groups
associated   with  the  irreducible  components  of  `R`.  The  irreducible
crystallographic  root  systems  are  classified  by  the following list of
*Dynkin  diagrams*. We show the labeling of the nodes given by the function
'diagram' described below.

```
A_n O—O—O—…—O   B_n O⇐O—O—…—O  C_n O⇒O—O—…—O  D_n  O 2
    1 2 3 … n       1 2 3 … n      1 2 3 … n       ￨
                                                 O—O—…—O
                                                 1 3 … n

G_2 O⇛O   F_4 O—O⇒O—O  E_6   O 2   E_7   O 2     E_8    O 2
    1 2       1 2 3 4        ￨           ￨              ￨
                         O—O—O—O—O   O—O—O—O—O—O    O—O—O—O—O—O—O
                         1 3 4 5 6   1 3 4 5 6 7    1 3 4 5 6 7 8
```

These diagrams encode the presentation of the Coxeter group `W` as follows:
the vertices represent the generating reflections; an edge is drawn between
`s`  and `t` if the order `m(s,t)` of `st` is greater than `2`; the edge is
single  if  `m(s,t)=3`,  double  if  `m(s,t)=4`,  triple if `m(s,t)=6`. The
arrows  indicate the relative root lengths when `W` has more than one orbit
on  `R`, as explained below; we  get the *Coxeter Diagram*, which describes
the  underlying Weyl group, if  we ignore the arrows:  we see that the root
systems `B_n` and `C_n` correspond to the same Coxeter group.

Let us now describe how the root systems are encoded in these diagrams. Let
`R`  be a root system in `V`. Then we can choose a linear form on `V` which
vanishes  on no element of `R`. According to  the sign of the value of this
linear  form on a root  `r ∈ R` we  call `r` *positive* or *negative*. Then
there  exists a unique subset `Π` of  the positive roots, called the set of
*simple  roots*, such that  any positive root  is a linear combination with
non-negative  coefficients of roots  in `Π`. It  can be shown  that any two
sets of simple roots (corresponding to different choices of linear forms as
above)  can be transformed into  each other by a  unique element of `W(R)`.
Hence, since the pairing between `V` and `Vⱽ` is `W`-invariant, if `Π` is a
set  of simple roots and if we define  the *Cartan matrix* as being the `n`
times  `n` matrix `C={rⱽ(s)}_{rs}`, for `r,s∈Π` this matrix is unique up to
simultaneous  permutation of rows and columns.  It is precisely this matrix
which is encoded in a Dynkin diagram, as follows.

The  indices for the rows of `C` label the nodes of the diagram. The edges,
for  `r ≠ s`, are  given as follows. If  `C_{rs}` and `C_{sr}` are integers
such  that `|C_{rs}| ≥  |C_{sr}|` the vertices  are connected by `|C_{rs}|`
lines,  and if `|C_{rs}|>1`  then we put  an additional arrow  on the lines
pointing towards the node with label `s`. In all other cases, we simply put
a  single line  equipped with  the unique  integer `p_{rs}  ≥ 1`  such that
`C_{rs}C_{sr}=cos^2 (π/p_{sr})`.

It  is now important to note that, conversely, the whole root system can be
recovered from the simple roots. The reflections in `W(R)` corresponding to
the  simple roots are  called *simple* reflections.  They are precisely the
generators  for which the Coxeter diagram encodes the defining relations of
`W(R)`. Each root is in the orbit of a simple root, so that `R` is obtained
as  the orbit of the  simple roots under the  group generated by the simple
reflections.  The  restriction  of  the  simple  reflections  to  `V_R`  is
determined  by the Cartan matrix, so `R` is determined by the Cartan matrix
and the set of simple roots.

The  Cartan  matrix  corresponding  to  one  of  the above irreducible root
systems  (with the specified labeling) is  returned by the command 'cartan'
which takes as input a symbol giving the type (that is ':A', ':B', …, ':I')
and a positive integer giving the rank. This function returns a matrix with
entries  in `bbZ`.  Given two  Cartan matrices  `c1` and `c2`, their matrix
direct sum (corresponding to the orthogonal direct sum of the root systems)
can be produced by the `cat(c1,c2,dims=[1,2])`.

The function 'WeylGroup' takes as input some data which determine the roots
and  the coroots and  produces a `struct`  containing information about the
root  system `R` and about `W(R)`. If we label the positive roots by '1:N',
and  the  negative  roots  by  'N+1:2N',  then  each  simple  reflection is
represented by the permutation of '1:2N' which it induces on the roots.

The function 'WeylGroup' has 2 methods; in one of them, the argument is the
Cartan  matrix of the root system. This  constructs a root system where the
simple  roots are the canonical basis of `V`, and the matrix of the coroots
expressed in the dual basis of `Vⱽ` is then equal to the Cartan matrix.

If one only wants to work with Cartan matrices with a labeling as specified
by  the  above  list,  the  function  call  can  be  simplified. Instead of
'WeylGroup(CartanMat(:D,4))' the following is also possible.

```
julia> W=WeylGroup(:D,4)
WeylGroup(:D,4)

julia> W.cartan
4×4 Array{Int8,2}:
  2   0  -1   0
  0   2  -1   0
 -1  -1   2  -1
  0   0  -1   2
```

Also,  the Weyl group struct associated to a direct sum of irreducible root
systems can be obtained as a product

```
julia> W=WeylGroup(:A,2)*WeylGroup(:B,2)
WeylGroup(:A,<1,2>)xWeylGroup(:B,<3,4>)

julia> W.cartan
4×4 Array{Int8,2}:
  2  -1   0   0
 -1   2   0   0
  0   0   2  -2
  0   0  -1   2
```

The  same `struct`  is constructed  by applying  'WeylGroup' to  the matrix
'cat(cartan(:A,2), cartan(:B,2),dims=[1,2])'.

The elements of a Weyl group are permutations of the roots:
```julia-repl
julia> W=WeylGroup(:D,4)
WeylGroup(:D,4)

julia> p=element(W,[1,3,2,1,3])
{UInt8}(1,14,13,2)(3,17,8,18)(4,12)(5,20,6,15)(7,10,11,9)(16,24)(19,22,23,21)

julia> word(W,p)
5-element Array{Int64,1}:
 1
 3
 1
 2
 3

```
This module is mostly a port of the basic functions on Weyl groups in CHEVIE.
The dictionary from CHEVIE is as follows:
```
     CartanMat("A",5)                       →  cartan(:A,5) 
     CoxeterGroup("A",5)                    →  WeylGroup(:A,5) 
     Size(W)                                →  length(W) 
     ForEachElement(W,f)                    →  for w in W f(w) end 
     ReflectionDegrees(W)                   →  degrees(W) 
     IsLeftDescending(W,w,i)                →  isleftdescent(W,w,i) 
     ReflectionSubgroup : only standard parabolics now
     TwoTree(m)                             →  twotree(m) 
     FiniteCoxeterTypeFromCartanMat(m)      →  type_cartan(m) 
     RootsCartan(m)                         →  roots(m) 
     PrintDiagram(W)                        →  diagram(W) 
     Inversions                             →  inversions 
     Reflection                             →  reflection 
     W.orbitRepresentative[i]               →  simple_representative(W,i) 
```
finally, a benchmark on julia 1.0.2
```benchmark
julia> @btime length(elements(WeylGroup(:E,7)))
  807.009 ms (8608502 allocations: 787.93 MiB)
```
GAP3 for the same computation takes 2.2s
"""
module Weyl

export cartan, WeylGroup, diagram, inversions, simple_conjugating_element, 
  two_tree, CharTable

using Gapjm, LinearAlgebra
#------------------------ Cartan matrices ----------------------------------
"""
    `cartan(type, rank)`

Cartan matrix for a Weyl group:

```julia-repl
julia> cartan(:A,4)
4×4 Array{Int8,2}:
  2  -1   0   0
 -1   2  -1   0
  0  -1   2  -1
  0   0  -1   2
```
"""
function cartan(t::Symbol,r::Int=0)
  if t==:A return Array{Int8,2}(SymTridiagonal(fill(2,r),fill(-1,r-1))) end
  m=cartan(:A,r) 
  if t==:B m[1,2]=-2 
  elseif t==:C m[2,1]=-2 
  elseif t==:D m[1:3,1:3]=[2 0 -1; 0 2 -1;-1 -1 2]
  elseif t==:E m[1:4,1:4]=[2 0 -1 0; 0 2 0 -1;-1 0 2 -1;0 -1 -1 2]
  elseif t==:F m[3,2]=-2 
  elseif t==:G m[2,1]=-3
  end
  m
end

"""
    two_tree(m)

 Given  a square  matrix m  with zeroes  (or falses,  for a boolean matrix)
 symmetric  with respect to the diagonal, let  G be the graph with vertices
 axes(m)[1] and an edge between i and j iff !iszero(m[i,j]).
 If G  is a line this function returns it as a Vector{Int}. 
 If  G  is  a  tree  with  one  vertex  c of valence 3 the function returns
 (c,b1,b2,b3)  where b1,b2,b3 are  the branches from  this vertex sorted by
 increasing length.
 Otherwise the function returns `nothing`
```julia-repl
julia> CoxGroups.two_tree(cartan(:A,4))
4-element Array{Int64,1}:
 1
 2
 3
 4

julia> CoxGroups.two_tree(cartan(:E,8))
(4, [2], [3, 1], [5, 6, 7, 8])
```
"""
two_tree=function ( m )
  branch= function ( x )
    while true
      x=findfirst(i->m[x,i]!=0 && !(i in line),1:r)
      if !isnothing(x) push!(line,x) else break end
    end
  end
  r=size(m,1)
  line=[1]
  branch(1)
  l=length(line)
  branch(1)
  line=vcat(line[end:-1:l+1],line[1:l])
  l=length(line)
  if any(i->any(j->m[line[j],line[i]]!=0,1:i-2),1:l) return nothing end
  if l==r return line end
  p = findfirst(x->any(i->!(i in line)&&(m[x,i]!=0),1:r),line)
  branch(line[p])
  if length( line ) != r  return nothing end
  (line[p],sort([line[p-1:-1:1],line[p+1:l],line[l+1:r]], by=length)...)
end

"""
    type_cartan(C)

 return a list of (series=s,indices=[i1,..,in]) for a Cartan matrix
"""
function type_cartan(m::Matrix{<:Integer})
" (series,rank) for an irreducible Cartan matrix"
function type_irred_cartan(m::Matrix{<:Integer})
  rank=size(m,1)
  l=two_tree(m)
  if isnothing(l) return nothing end
  if l isa Tuple # types D,E
    (vertex,b1,b2,b3)=l
    if length(b2)==1 series=:D 
      indices=vcat(b1,b2,[vertex],b3) 
    else series=:E 
      indices=vcat([b2[2],b1[1],b2[1],vertex],b3) 
    end 
  else  # types A,B,C,F,G
    n=m[l,l] 
    co=i->n[i,i+1]*n[i+1,i] 
    function rev()
      l=l[end:-1:1] 
      n=m[l,l]
    end
    if rank==1 series=:A 
    elseif rank==2 
      if co(1)==1 series=:A 
      elseif co(1)==2 series=:B  
        if n[1,2]==-1 rev() end # B2 preferred to C2
      elseif co(1)==3 series=:G  
        if n[1,2]!=-1 rev() end 
      end
    else
      if co(rank-1)!=1 rev() end 
      if co(1)==1
        if co(2)==1 series=:A 
        else series=:F
          if n[2,3]==-2 rev() end 
        end 
      elseif n[1,2]==-2 series=:B
      else series=:C
      end  
    end 
    indices=l 
  end 
  if cartan(series,rank)!=m[indices,indices] return nothing end  # countercheck
  (series,indices)
end
  map(blocks(m))do I
    (t,indices)=type_irred_cartan(m[I,I])
    (series=t,indices=I[indices])
  end
end

function diagram(;series::Symbol,indices::AbstractVector{Int})
  ind=repr.(indices)
  l=length.(ind)
  bar(n)="\u2014"^n
  rdarrow(n)="\u21D0"^(n-1)*" "
  ldarrow(n)="\u21D2"^(n-1)*" "
  tarrow(n)="\u21DB"^(n-1)*" "
  vbar="\UFFE8" # "\u2503"
  node="O"
  if series==:A
    println(join(map(l->node*bar(l),l[1:end-1])),node)
    println(join(ind," "))
  elseif series==:B
    println(node,rdarrow(max(l[1],2)),join(map(l->node*bar(l),l[2:end-1])),node)
    println(ind[1]," "^max(3-l[1],1),join(ind[2:end]," "))
  elseif series==:C
    println(node,ldarrow(max(l[1],2)),join(map(l->node*bar(l),l[2:end-1])),node)
    println(ind[1]," "^max(3-l[1],1),join(ind[2:end]," "))
  elseif series==:D
    println(" "^l[1]," O $(ind[2])\n"," "^l[1]," ",vbar)
    println(node,bar(l[1]),map(l->node*bar(l),l[3:end-1])...,node)
    println(ind[1]," ",join(ind[3:end]," "))
  elseif series==:E
    dec=2+l[1]+l[3]
    println(" "^dec,"O $(ind[2])\n"," "^dec,vbar)
    println(node,bar(l[1]),node,bar(l[3]),
              join(map(l->node*bar(l),l[4:end-1])),node)
    println(join(vcat(ind[1],ind[3:end])," "))
  elseif series==:F
    println(node,bar(l[1]),node,ldarrow(max(l[2],2)),node,bar(l[3]),node)
    println(ind[1]," ",ind[2]," "^max(3-l[2],1),ind[3]," ",ind[4])
  elseif series==:G
    println(node,tarrow(max(l[1],2)),node)
    println(ind[1]," "^max(3-l[1],1),ind[2])
  end
end

function diagram(v::Vector{NamedTuple{(:series,:indices),
                                      Tuple{Symbol,Vector{Int}}}})
  for t in v diagram(;t...) end
end

"""
    roots(C)

 return the set of positive roots defined by integral Cartan matrix C
"""
function roots(C::Matrix{<:Integer})
  l=size(C,1)
  Pi=one(C)
  Pi=[Pi[i,:] for i in 1:l]
  R=typeof(Pi)[]
  old=Pi
  while length(old)!=0
    old=union(old)
    push!(R,old)
    new=eltype(Pi)[]
    for r in old, j in 1:l
      p =C[j,:]'*r
      if p<0 || (length(R)>p+1 && r[j]>p && (r-(p+1)*Pi[j] in R[end-p-1]))
        push!(new,r+Pi[j])
      end
    end
    old=new
  end
  vcat(R...)
end
#--------------- Weyl groups --------------------------------------------
struct WeylGroup <: CoxeterGroup{Perm{UInt8}}
  G::PermGroup{UInt8}
  matgens::Vector{Matrix{Int8}}
  roots::Vector{Vector{Int8}}
  parentN::Int
  N::Int
  cartan::Array{Int8,2}
  weyltype::Vector{NamedTuple{(:series,:indices),Tuple{Symbol,Vector{Int}}}}
  prop::Dict{Symbol,Any}
end

" the matrix of the reflection of given root and coroot"
function CoxGroups.reflection(root::Vector,coroot::Vector)
  root,coroot=promote(root,coroot)
  m=[i*j for i in coroot, j in root]
  one(m)-m
end

" Weyl group from Cartan matrix"
function WeylGroup(C::Matrix{T}) where T<:Integer
  r=roots(C)
  N=length(r)
  r=vcat(r,-r)
  # the reflections defined by Cartan matrix C
  matgens=[reflection(one(C)[i,:],C[i,:]) for i in axes(C,1)]
  """
    the permutations of the roots r effected by the matrices mm
  """
  function perms(r,mm)
    s=Perm(Vector{UInt8}(sortperm(r)))
    map(mm)do m
      inv(Perm(Vector{UInt8}(sortperm(map(v->(v'*m)',r)))))*s
    end
  end
  gens=perms(r,matgens)
  rank=size(C,1)
  t=type_cartan(C)
  name=join(map(t)do (series,indices)
             n="W($series"*TeXstrip("_{$(length(indices))}")
    if indices!=eachindex(indices)
      ind=any(i->i>10,indices) ? join(indices,",") : join(indices)
      n*=TeXstrip("_{($ind)}")
    end
    n*")"
  end,"×")
  WeylGroup(PermGroup(gens),matgens,r,N,N,C,t,Dict(:name=>name))
end

"Weyl group from type"
WeylGroup(t::Symbol,r::Int=0)=WeylGroup(cartan(t,r))

(W::WeylGroup)(x...)=element(W.G,collect(x))
(W::WeylGroup)(x::AbstractVector)=element(W.G,x)

"diagram of Weyl group"
diagram(W::WeylGroup)=diagram(W.weyltype)

"number of reflections of W"
CoxGroups.nref(W::WeylGroup)=W.N

function matX(W::WeylGroup,w)
  vcat(permutedims(hcat(W.roots[(1:coxrank(W)).^Ref(w)]...)))
end

function cartancoeff(W::WeylGroup, i, j )
  v=findfirst(x->!iszero(x),W.roots[i])
  j=W.roots[j]-W.roots[j^reflection(W,i)]
  div(j[v],W.roots[i][v])
end

function CoxGroups.ReflectionSubgroup(W::WeylGroup,I::AbstractVector{Int})
  if all(i->i in 1:coxrank(W),I)
    roots=filter(r->all(i->i[1] in I || iszero(i[2]),enumerate(r)),W.roots)
    n=any(i->i>=10,I) ? join(I,",") : join(I)
    n=name(W)*TeXstrip("_{$n}")
    WeylGroup(PermGroup(coxgens(W)[I]),W.matgens[I],roots,W.parentN,
      div(length(roots),2),
      W.cartan[I,I],type_cartan(W.cartan[I,I]),Dict(:name=>n))
  else
    G=PermGroup(map(i->reflection(W,i),I))
    o=vcat(orbits(G,I)...)
    sort!(o)
    N=div(length(o),2)
    I=filter(o[1:N]) do i
      cnt=0
      r=reflection(W,i)
      for j in o[1:N]; if j!=i && j^r>W.parentN return false end end
      return true
    end
    n=any(i->i>=10,I) ? join(I,",") : join(I)
    n=name(W)*TeXstrip("_{$n}")
    gens=reflection.(Ref(W),I)
    C=[cartancoeff(W,i,j) for i in I, j in I]
    WeylGroup(PermGroup(gens),matX.(Ref(W),gens),W.roots[o],W.parentN,
      N,C,type_cartan(C),Dict(:name=>n,:inclusion=>o))
  end
end

function Base.:*(W::WeylGroup...)
  WeylGroup(cat(map(x->x.cartan,W)...,dims=[1,2]))
end

CoxGroups.isleftdescent(W::WeylGroup,w,i::Int)=i^w>W.parentN

inclusion(W::WeylGroup)::Vector{Int}=W.N==W.parentN ? (1:W.N) :
  W.prop[:inclusion][1:W.N]

inversions(W::WeylGroup,w)=[i for i in inclusion(W) if isleftdescent(W,w,i)]

Base.length(W::WeylGroup,w)=count(i->isleftdescent(W,w,i),1:W.N)

CoxGroups.name(W::WeylGroup)::String=W.prop[:name]

"""
  The reflection degrees of W
"""
function Gapjm.degrees(W::WeylGroup)
  l=sort(map(length,values(groupby(sum,W.roots))),rev=true)
  reverse(1 .+div.(conjugate_partition(l),2))
end

Base.length(W::WeylGroup)=prod(degrees(W))

function root_representatives(W::WeylGroup)
  reps=fill(0,length(W.roots))
  repelts=fill(one(W),length(W.roots))
  for i in eachindex(coxgens(W))
    if iszero(reps[i])
      d=orbit_and_representative(W.G,i)
      for (n,e) in d 
        reps[n]=i
        repelts[n]=e
      end
    end
  end
  W.prop[:rootreps]=reps
  W.prop[:repelms]=repelts
  W.prop[:reflections]=map((i,p)->coxgens(W)[i]^p,reps,repelts)
end

"for each root index of simple representative"
function CoxGroups.simple_representative(W::WeylGroup,i)::Int
  getp(root_representatives,W,:rootreps)[i]
end
  
"for each root element conjugative representative to root"
function simple_conjugating_element(W::WeylGroup,i)::Perm{UInt8}
  getp(root_representatives,W,:repelms)[i]
end

function CoxGroups.reflection(W::WeylGroup,i::Integer)::Perm{UInt8}
  getp(root_representatives,W,:reflections)[i]
end

struct CharTable{T}
  irr::Matrix{T}
  charnames::Vector{String}
  classnames::Vector{String}
  centralizers::Vector{Int}
  identifier::String
end

function Base.show(io::IO,ct::CharTable)
  println(io,"CharTable(",ct.identifier,")")
  format(io,ct.irr,row_labels=map(TeXstrip,ct.charnames),
                column_labels=map(TeXstrip,ct.classnames))
end

end
