"""
Finite Coxeter groups are the finite complex reflection groups which can be
defined on a real vector space `V`.

*Weyl  groups* are the  finite Coxeter groups  which can be  defined over a
rational vector space `V` (thus over the integers).

Like   finite  complex   reflection  groups,   finite  Coxeter  groups  are
implemented  as groups  of permutations  of a  set of roots. The particular
*crystallographic*  root systems for Weyl groups  play an important role in
mathematics as they classify semi-simple Lie algebras and algebraic groups.

Let  us give precise definitions.  Let `V` be a  real vector space and `Vⱽ`
its  dual. A *root system* is a finite set of vectors `R⊂ V` (the *roots*),
together  with  a  map  `r↦  rⱽ`  from  `R`  to  a subset `Rⱽ` of `Vⱽ` (the
*coroots*) such that:

  - For any `r∈ R`,  we have `rⱽ(r)=2`, so  that the formula `x↦ x-rⱽ(x)r`
    defines  a  reflection  `sᵣ:V→  V`  with  root  `r`  and coroot `rⱽ`.
  - The reflection `sᵣ` stabilizes `R`.

The  subgroup `W=W(R)` of `GL(V)` generated by the reflections `sᵣ` for `r∈
R`  is a finite Coxeter  group. We require *reduced*  root systems, that is
such  that the only elements of `R` colinear  with `r∈ R` are `r` and `-r`;
for Weyl groups, we also require that the root system be
*crystallographic*, that is `rⱽ(s)` is an integer, for any `s∈ R,rⱽ∈ Rⱽ`.

If  we identify  `V` with  `Vⱽ` by  choosing a  `W`-invariant bilinear form
`(.;.)`,  then we have `rⱽ=2r/(r;r)`. A root system `R` is *irreducible* if
`R`   is  not  the  union  of  two  orthogonal  subsets;  equivalently  the
representation  of `W` on the subspace  generated by `R` is irreducible. If
`R` is reducible then the corresponding Coxeter group is the direct product
of the Coxeter groups associated with the irreducible components of `R`.

Let  us  now  describe  how  a  root  system  `R` and a presentation of the
corresponding  `W` are encoded in  a Cartan matrix or  a Dynkin diagram. We
can  choose a linear  form on `V`  which does not  vanish on any element of
`R`.  Depending on the sign of the value of this linear form on a root `r ∈
R`  we call `r` *positive* or *negative*. Then there exists a unique subset
`Π`  of the positive roots, the *simple roots*, such that any positive root
is a linear combination with non-negative coefficients of the roots in `Π`.
Any  two sets of simple roots (corresponding to different choices of linear
forms) can be transformed into each other by a unique element of `W(R)`. If
`S`  is  the  set  of  reflections  with  respect to the simple roots, then
`(W,S)`  is  a  Coxeter  system.  These  generating  reflections are called
*Coxeter generators* or *simple reflections*.

Since the pairing between `V` and `Vⱽ` is `W`-invariant, if `Π` is a set of
simple  roots and if we  define the *Cartan matrix*  as being the `n` times
`n`  matrix `C={rⱽ(r')}`  for `r,r'∈Π`,  this matrix  is independent of the
chosen  linear form  up to  simultaneous permutation  of rows  and columns.
Since the action of `sᵣ` on `r'` for `r,r'∈Π` is given by
`sᵣ(r')=r'-C(r,r')r`,   the   Cartan   matrix   determines  the  reflection
representation of `W`.

For  a crystallographic root system the Cartan matrix has integral entries,
and  in the basis `Π` (completed by a basis of the orthogonal), `sᵣ` has an
integral  matrix.  All  finite-dimensional  (complex)  representations of a
finite  Coxeter  group  can  be  realized  over  the field generated by the
entries of the Cartan matrix.

The  Cartan matrix is encoded  in a *Dynkin diagram*,  a tree with weighted
edges  and  an  orientation  on  edges  of  even weight >2, as follows. The
vertices are indexed by the simple reflections; an edge is drawn between
`s` and `t` if the order `mₛₜ` of `st` is greater than `2` and is given the
weight  `mₛₜ`. These  weights are  encoded by  drawing the  edge single for
weight  3, double for weight 4 and triple for weight 6. The arrows indicate
the relative root lengths (going from the longer to the shorter root) which
may  differ between different orbits of  `W` on `R`. Alternately the Dynkin
diagram  can be obtained  from the Cartan  matrix as follows:  if `Cᵣₛ` and
`Cₛᵣ`  are integers  such that  `|Cₛᵣ|≥|Cᵣₛ|=1` there  is an edge of weight
`|Cₛᵣ|`  from `r` to `s`  with an arrow pointing  to `s` if `|Cₛᵣ|>1`. Note
that  the Cartan matrices  we consider here  are not necessarily symmetric,
contrary  to the  Cartan matrices  we considered  describing the reflection
representation  of a general Coxeter  group; being symmetric corresponds to
all roots being taken of the same length.

The  irreducible  crystallographic  root  systems  are  classified  by  the
following  list of Dynkin diagrams. The labeling  of the nodes is the order
of the generators and is shown by the function `diagram`.
```
Aₙ O—O—O—…—O   Bₙ O⇐ O—O—…—O  Cₙ O⇒ O—O—…—O  Dₙ O 2
   1 2 3 … n      1  2 3 … n     1  2 3 … n     ￨
                                              O—O—…—O
                                              1 3 … n

G₂ O⇛ O  F₄ O—O⇒O—O    E₆  O 2    E₇  O 2      E₈  O 2
   1  2     1 2 3 4        ￨          ￨            ￨
                       O—O—O—O—O  O—O—O—O—O—O  O—O—O—O—O—O—O
                       1 3 4 5 6  1 3 4 5 6 7  1 3 4 5 6 7 8
```
We get the *Coxeter diagram*, which describes the underlying Weyl group, if
we  ignore  the  arrows:  we  see  that  the  root  systems `B_n` and `C_n`
correspond to the same Coxeter group (the Coxeter diagram is defined by the
*Coxeter  matrix*). Weyl  groups can  also be  characterized as  the finite
Coxeter groups such that all off-diagonal entries of the Coxeter matrix are
in `{2,3,4,6}`.

Here  are the Coxeter diagrams for the  finite Coxeter groups which are not
crystallographic (`I₂(e)` is not crystallographic if `e∉ {2,3,4,6}`).
```
       e        5         5
I₂(e) O—O   H₃ O—O—O  H₄ O—O—O—O
      1 2      1 2 3     1 2 3 4
```
The function `cartan` gives the cartan matrix for an irreducible root system

```julia-repl
julia> cartan(:D,4)
4×4 Matrix{Int64}:
  2   0  -1   0
  0   2  -1   0
 -1  -1   2  -1
  0   0  -1   2

julia> cartan(:I,2,5) # for type I₂(e) give e as 3rd argument
2×2 Matrix{Cyc{Int64}}:
       2  ζ₅²+ζ₅³
 ζ₅²+ζ₅³        2
```

Given   two  Cartan  matrices  `c1`  and  `c2`,  their  matrix  direct  sum
(corresponding  to the  orthogonal direct  sum of  the root systems) can be
obtained by `cat(c1,c2,dims=[1,2])`.

The  whole  root  system  can  be  recovered  from the simple roots and the
corresponding  coroots, since each root  is in the orbit  of a simple root.
The  restriction of the simple reflections to the span of `R` is determined
by the Cartan matrix, so `R` is determined by the Cartan matrix and the set
of simple roots.

The  function  `rootdatum`  takes  as  arguments  a  matrix `r` whose lines
represents the list of simple roots and another matrix `cr` whose lines are
the  corresponding  coroots  and  produces  a `FiniteCoxeterGroup`. Such an
object  is a permutation group containing  in addition the description of a
root  system. `cr*transpose(r)` should be a  Cartan matrix. Each element of
the  coxeter  group  is  represented  as  the permutation it induces on the
roots,  coded as a permutation of `1:2N`  where we label the positive roots
by `1:N`, and the negative roots by `N+1:2N`.

If a single matrix argument is given to `rootdatum` it is taken as `cr` and
`r`  is taken  to be  the identity  matrix; we  get thus  a particular root
system  where the  roots are  the canonical  basis of `V`. For convenience,
`rootdatum(cartan(t...))` can be simplified to `coxgroup(t...)`.

```julia-repl
julia> W=coxgroup(:D,4) # same as rootdatum(cartan(:D,4))
D₄

julia> cartan(W)
4×4 Matrix{Int64}:
  2   0  -1   0
  0   2  -1   0
 -1  -1   2  -1
  0   0  -1   2
```

Also,  the `FiniteCoxeterGroup` associated to a direct sum of irreducible root
systems can be obtained as

```julia-repl
julia> W=coxgroup(:A,2)*coxgroup(:B,2)
A₂×B₂

julia> cartan(W) # same as cat(cartan(:A,2), cartan(:B,2),dims=[1,2])
4×4 Matrix{Int64}:
  2  -1   0   0
 -1   2   0   0
  0   0   2  -2
  0   0  -1   2
```
The elements of a Weyl group are permutations of the roots:
```julia-repl
julia> W=coxgroup(:D,4)
D₄

julia> p=W(1,3,2,1,3) # permutes the 24 roots
(1,14,13,2)(3,17,8,18)(4,12)(5,20,6,15)(7,10,11,9)(16,24)(19,22,23,21)

julia> word(W,p)
5-element Vector{Int64}:
 1
 3
 1
 2
 3
```

finally, a benchmark on julia 1.0.2
```benchmark
julia> @btime length(elements(coxgroup(:E,7)))
  531.385 ms (5945569 allocations: 1.08 GiB)
```
GAP3 for the same computation takes 2.2s
"""
module Weyl

export two_tree, rootdatum, torus, istorus, derived_datum,
 dim, dimension, with_inversions, standard_parabolic, describe_involution,
 relative_group, rootlengths, highest_short_root, badprimes,
 ComplexReflectionGroup

using ..Chevie
#------------------------ Cartan matrices ----------------------------------

"""
`cartan(M::AbstractMatrix)` Cartan matrix from Coxeter matrix

The  argument should be the Coxeter matrix  `M` for a Coxeter group `W` and
the   result  is  the  Cartan  Matrix   `C`  for  the  standard  reflection
representation  of `W`. We have `C[s,t]=-2cos(π/M[s,t])`, where `M[s,s]==1`
and  by  convention  `π/M[s,t]==0`  if  `M[s,t]==∞`,  which we represent by
`M[s,t]==0`.  Since  `M`  is  symmetric,  the  resulting  `C` is symmetric,
meaning  that all roots  in the constructed  reflection representation have
same length.

```julia-repl
julia> cartan([1 3;3 1])
2×2 Matrix{Cyc{Int64}}:
  2  -1
 -1   2
```
"""
PermRoot.cartan(m::AbstractMatrix)=map(c->iszero(c) ? -2 : -E(2*c,-1)-E(2*c),m)

const cartanmats=Dict{Symbol,Function}()

cartanmats[:a]=cartanmats[:A]=function(r)
  m=fill(0,r,r)
  for i in 1:r
    m[i,i]=2;if i<r m[i,i+1]=-1 end;if i>1 m[i,i-1]=-1 end
  end
  m
end
cartanmats[:b]=cartanmats[:B]=function(r,cartanType=2)
  m=cartanmats[:A](r)
  if r>1 
    if cartanType==2 m[1,2]=-2 
    else m=typeof(cartanType//1).(m);m[1,2]=-cartanType;m[2,1]=-2//cartanType
    end
  end
  m
end
cartanmats[:bsym]=cartanmats[:Bsym]=function(r)
  m=Cyc.(cartanmats[:A](r))
  m[1:2,1:2]=cartanmats[:Isym](4)
  m
end
cartanmats[:c]=cartanmats[:C]=function(r)
  m=cartanmats[:A](r)
  if r>1 m[2,1]=-2 end
  m
end
cartanmats[:d]=cartanmats[:D]=function(r)
  m=cartanmats[:A](r)
  u=min(r,3); m[1:u,1:u]=[2 0 -1; 0 2 -1;-1 -1 2][1:u,1:u]
  m
end
cartanmats[:e]=cartanmats[:E]=function(r)
  if r>8 error("type :E is defined only for rank ≤8") end
  m=cartanmats[:A](r)
  u=1:min(r,4); m[u,u]=[2 0 -1 0; 0 2 0 -1;-1 0 2 -1;0 -1 -1 2][u,u]
  m
end
cartanmats[:f]=cartanmats[:F]=function(r,cartanType=1)
  if r!=4 error("type :F is defined only for rank 4") end
  m=cartanmats[:A](r)
  if cartanType==1 m[3,2]=-2 
  else m=typeof(cartanType//1).(m);m[3,2]=-2//cartanType;m[2,3]=-cartanType
  end
  m
end
cartanmats[:fsym]=cartanmats[:Fsym]=function(r)
  if r!=4 error("type :Fsym is defined only for rank 4") end
  m=Cyc.(cartanmats[:A](r))
  m[3,2]=m[2,3]=-root(2)
  m
end
cartanmats[:g]=cartanmats[:G]=function(r,cartanType=1)
  if r!=2 error("type :G is defined only for rank 2") end
  m=cartanmats[:A](r)
  if cartanType==1 m[2,1]=-3
  else m=typeof(cartanType//1).(m);m[1,2]=-cartanType;m[2,1]=-3//cartanType
  end
  m
end
cartanmats[:gsym]=cartanmats[:Gsym]=r->cartanmats[:Isym](6)
cartanmats[:h]=cartanmats[:H]=function(r)
  if r>4 error("type :H is defined only for rank ≤4") end
  m=Cyc.(cartanmats[:A](r))
  m[1:2,1:2]=cartanmats[:Isym](5)
  m
end
cartanmats[:isym]=cartanmats[:Isym]=function(bond)
  c=improve_type(-E(2*bond)-E(2*bond,-1))
  [2 c;c 2]
end
cartanmats[:i]=cartanmats[:I]=function(r,bond,cartanType=1)
  if r!=2 error("$r should be 2") end
  if isodd(bond) cartanmats[:Isym](bond)
  elseif bond==2 
    [2 0;0 2]
  else 
    improve_type([2 -cartanType;(-2-E(bond)-E(bond,-1))//cartanType 2])
  end
end
"""
`cartan(type, rank [,bond])`

the  Cartan matrix for a  finite Coxeter group described  by type and rank.
The  recognized types are `:A, :B, :Bsym, :C, :D, :E, :F, :Fsym, :G, :Gsym,
:I,  :H`. For type `:I` a third  argument must be given describing the bond
between the two generators. The `sym` types correspond to
(non-crystallographic)  root systems where all  roots have the same length;
they  afford automorphisms that  the crystallographic root  system does not
afford, which allow to define the "very twisted" Chevalley groups.

```julia-repl
julia> cartan(:F,4)
4×4 Matrix{Int64}:
  2  -1   0   0
 -1   2  -1   0
  0  -2   2  -1
  0   0  -1   2

julia> cartan(:I,2,5)
2×2 Matrix{Cyc{Int64}}:
       2  ζ₅²+ζ₅³
 ζ₅²+ζ₅³        2

julia> cartan(:Bsym,2)
2×2 Matrix{Cyc{Int64}}:
   2  -√2
 -√2    2
```
"""
function PermRoot.cartan(t::Symbol,r::Integer,args::Number...)
  if haskey(cartanmats,t) cartanmats[t](r,args...)
  else error("Unknown Cartan type $(repr(t)). Known types are:\n",
             join(sort(repr.(collect(keys(cartanmats)))),", "))
  end
end

"""
`coxeter_matrix(type, rank [,bond])` or `coxmat`

Like `cartan`, the function `coxmat` can be defined from the type and rank
of a finite Coxeter group.
"""
CoxGroups.coxeter_matrix(t::Symbol,r::Integer,b::Integer=0)=coxeter_matrix(cartan(t,r,b))

"""
`two_tree(m)`

Given  a  square  matrix  `m`  with  zeroes  symmetric  with respect to the
diagonal,  let  `G`  be  the  graph  with vertices `axes(m)[1]` and an edge
between `i` and `j` iff `!iszero(m[i,j])`.

If  `G` is a line this function returns  it as a `Vector{Int}`. If `G` is a
tree with one vertex `c` of valence `3` the function returns `(c,b1,b2,b3)`
where  `b1,b2,b3` are  the branches  from this  vertex sorted by increasing
length. Otherwise the function returns `nothing`.

This function is used when recognizing the type of a Cartan matrix.
```julia-repl
julia> two_tree(cartan(:A,4))
4-element Vector{Int64}:
 1
 2
 3
 4

julia> two_tree(cartan(:E,8))
(4, [2], [3, 1], [5, 6, 7, 8])
```
"""
two_tree=function(m::AbstractMatrix)
  function branch(x)
    while true
      x=findfirst(i->m[x,i]!=0 && !(i in line),axes(m,2))
      if !isnothing(x) push!(line,x) else break end
    end
  end
  line=[1]
  branch(1)
  l=length(line)
  branch(1)
  line=vcat(line[end:-1:l+1],line[1:l])
  l=length(line)
  if any(i->any(j->m[line[j],line[i]]!=0,1:i-2),1:l) return nothing end
  r=size(m,1)
  if l==r return line end
  p = findfirst(x->any(i->!(i in line)&&(m[x,i]!=0),1:r),line)
  branch(line[p])
  if length(line)!=r return nothing end
  (line[p],sort([line[p-1:-1:1],line[p+1:l],line[l+1:r]], by=length)...)
end

" TypeIrred or nothing for an irreducible Cartan matrix"
function type_fincox_cartan(m::AbstractMatrix)
  rank=size(m,1)
  if !all(==(2),diag(m)) return nothing end
  s=two_tree(m)
  if isnothing(s) return nothing end
  t=TypeIrred(Dict{Symbol,Any}(:rank=>rank))
  if s isa Tuple # types D,E
    (vertex,b1,b2,b3)=s
    if length(b2)==1 t.series=:D
      t.indices=[b1;b2;vertex;b3]::Vector{Int}
    else t.series=:E
      t.indices=[b2[2];b1[1];b2[1];vertex;b3]::Vector{Int}
    end
  else  # types A,B,C,F,G,H,I
    l=i->m[s[i],s[i+1]]
    r=i->m[s[i+1],s[i]]
    if rank==1 t.series=:A
    elseif rank==2
      bond=l(1)*r(1)
      if bond==1 t.series=:A
      elseif bond==2 t.series=:B
        if l(1)==-1 reverse!(s) end # B2 preferred to C2
        t.cartanType=improve_type(-l(1))
      elseif bond==3 t.series=:G
        if r(1)==-1 reverse!(s) end
        t.cartanType=improve_type(-l(1))
      else n=conductor(bond)
        if r(1)==-1 reverse!(s) end
        if bond==2+E(n)+E(n,-1) bond=n else bond=2n end
        t.series=:I
        if bond%2==0 t.cartanType=-l(1) end
        t.bond=bond
      end
    else
      if l(rank-1)*r(rank-1)!=1 reverse!(s) end
      if l(1)*r(1)==1
        if l(2)*r(2)==1 t.series=:A
        else t.series=:F
          if r(2)==-1 reverse!(s) end
          t.cartanType=improve_type(-l(2))
        end
      else n=conductor(l(1)*r(1))
        if n==5 t.series=:H
        else t.series=:B
          t.cartanType=improve_type(-l(1))
        end
      end
    end
    t.indices=s::Vector{Int}
  end
  return t
# println("t=$t")
# println("indices=",t.indices]," cartan=",cartan(t)," m=$m")
  if cartan(t,permute=true)==m return t end  # countercheck
end

"""
    type_cartan(C)

 return a list of (series=s,indices=[i1,..,in]) for a Cartan matrix
"""
function type_cartan(m::AbstractMatrix)
  map(diagblocks(m)) do I
    t=type_fincox_cartan(m[I,I])
    if isnothing(t) return nothing end
    t[:indices].=I[t[:indices]]
    t
  end
end

"""
`roots(C::AbstractMatrix)`

returns the set of positive roots defined by the Cartan matrix `C`, which
should be the Cartan matrix of a finite Coxeter group.

For  an integer Cartan matrix, the returned  roots are sorted by height and
reverse lexicographically for a given height.
"""
function PermRoot.roots(C::AbstractMatrix)
  o=one(C)
  R=[o[i,:] for i in axes(C,1)] # fast way to get rows of one(C)
  set=Set(R)
  for a in R
    c=C*a
    for i in axes(C,1)
      if a[i]!=1 || c[i]!=2
        v=copy(a)
        v[i]-=c[i]
        if !(v in set) 
          push!(R,v) 
          push!(set,v)
        end
      end
    end
  end
  if all(isinteger,C) sort!(R,by=x->(sum(x),-x)) end
  # important that roots are sorted as in CHEVIE for e.g. KLeftCells to work
  R
end

#-------Finite Coxeter groups --- T=type of elements----T1=type of roots------
const ComplexReflectionGroup=Union{PermRootGroup,FiniteCoxeterGroup}
@GapObj struct FCG{T,T1}<:FiniteCoxeterGroup{Perm{T}}
  G::PRG{T1,T}
  rootdec::Vector{Vector{T1}}
  N::Int
end

@GapObj struct FCSG{T,T1} <: FiniteCoxeterGroup{Perm{T}}
  G::PRSG{T1,T}
  rootdec::Vector{Vector{T1}}
  N::Int
  parent::FCG{T,T1}
end

const FC=Union{FCG,FCSG}

"""
`inversions(W::FiniteCoxeterGroup, w::AbstractVector{<:Integer})`

Given  a word `w=[s₁,…,sₙ]` (a vector of integers) representing the element
`W(w...)`,  returns the inversions of  `w`, that is the  list of indices of
the reflections of `W` given by `W(s₁), W(s₁,s₂,s₁), …,
W(s₁,s₂,…,sₙ,sₙ₋₁,…,s₁)`.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> inversions(W,[2,1,2])
3-element Vector{Int64}:
 2
 4
 1
```
"""
CoxGroups.inversions(W::FC,w::AbstractVector{<:Integer})=
   map(i->action(W,w[i],W(w[i-1:-1:1]...)),eachindex(w))

"""
`with_inversions(W,N)`

`W`  should be  a finite  Coxeter group  and `N`  a subset  of `1:nref(W)`.
Returns  the  element  `w`  of  `W` such that `N==inversions(W,w)`. Returns
`nothing` if no such element exists.

```julia-repl
julia> W=coxgroup(:A,2)
A₂

julia> map(N->with_inversions(W,N),combinations(1:nref(W)))
8-element Vector{Union{Nothing, Perm{Int16}}}:
 ()
 (1,4)(2,3)(5,6)
 (1,3)(2,5)(4,6)
 nothing
 nothing
 (1,6,2)(3,5,4)
 (1,2,6)(3,4,5)
 (1,5)(2,4)(3,6)
```
"""
function with_inversions(W::FC,N)
  w=copy(one(W))
  n=copy(N)
  while !isempty(n)
    p=findfirst(<=(ngens(W)),n)
    if p===nothing return nothing end
    r=refls(W,n[p])
    deleteat!(n,p)
    n.=action.(Ref(W),n,r)
    Perms.mul!(w,r)
  end
  w
end

"""
`standard_parabolic(W,H)`

Given  a reflection subgroup `H` or the indices of its simple roots returns
`nothing` if `H` is not parabolic, otherwise returns `w` such that `H^w` is
a standard parabolic subgroup of `W`.

```julia-repl
julia> W=coxgroup(:E,6)
E₆

julia> R=reflection_subgroup(W,[20,30,19,22])
E₆₍₁₉‚₁‚₉‚₂₀₎=A₄₍₃₁₂₄₎Φ₁²

julia> p=standard_parabolic(W,R)
(1,4,49,12,10)(2,54,62,3,19)(5,17,43,60,9)(6,21,34,36,20)(7,24,45,41,53)(8,65,50,15,22)(11,32,31,27,28)(13,48,46,37,40)(14,51,58,44,29)(16,23,35,33,30)(18,26,39,55,38)(42,57,70,72,56)(47,68,67,63,64)(52,59,71,69,66)

julia> p==standard_parabolic(W,[19,1,9,20]) # can give inclusiongens
true

julia> reflection_subgroup(W,[20,30,19,22].^p) # same as R^p
E₆₍₂₄₅₆₎=A₄Φ₁²

julia> R=reflection_subgroup(W,[1,2,3,5,6,35])
E₆₍₁‚₃‚₂‚₃₅‚₅‚₆₎=A₂₍₁₃₎×A₂₍₂₆₎×A₂₍₄₅₎

julia> standard_parabolic(W,R)
```
"""
function standard_parabolic(W::FC,I::AbstractVector{<:Integer})
  if isempty(I) return Perm() end
  b=W.rootdec[I]
  heads=map(x->findfirst(!iszero,x),toL(rowspace(toM(b))))
  b=vcat(W.rootdec[setdiff(1:ngens(W),heads)],b)
  # complete basis of <roots(W,I)> with part of S to make basis
  l=map(eachrow(toM(W.rootdec[1:W.N])*inv(toM(b).*1//1)))do v
   for x in v
     if (x isa Rational && x<0) || (isreal(x) && real(x)<0) return true
     elseif (x isa Rational && x>0) || (isreal(x) && real(x)>0) return false
     end
   end end
  N=(1:W.N)[l]
# find negative roots for associated order and make order standard
  w=with_inversions(W,N)
  if issubset(action.(Ref(W),I,w),eachindex(gens(W))) return w
  else return nothing
  end
end

standard_parabolic(W::FC,H::FC)=
  if !all(isinteger,cartan(W)) standard_parabolic(W.G,H.G)
  else standard_parabolic(W,inclusiongens(H,W))
  end

"""
`badprimes(W)`

Let  `W`  be  a  Weyl  group.  A  prime  is  *bad*  for `W` if it divides a
coefficient  of some  root on  the simple  roots. The  function `badprimes`
returns the list of primes which are bad for `W`.

```julia-repl
julia> W=coxgroup(:E,8)
E₈

julia> badprimes(W)
3-element Vector{Int64}:
 2
 3
 5
```
"""
function badprimes(W::FC)
  if isempty(W.rootdec) return Int[] end
  collect(setdiff(vcat(collect.(keys.(factor.(Int.(W.rootdec[nref(W)]))))...),[0,-1]))
end

"""
`describe_involution(W,w)`

Given  an  involution  `w`  of  a  Coxeter  group  `W`,  by  a  theorem  of
[Richardson1982](biblio.htm#rich82)  there is  a unique  parabolic subgroup
`P`  of `W` such that `P` is finite  and `w` is the longest element of `P`,
and is central in `P`. The function returns `I` such that
`P==reflection_subgroup(W,I)` and `w==longest(reflection_subgroup(W,I))`.

```julia-repl
julia> W=coxgroup(:A,2)
A₂

julia> w=longest(W)
(1,5)(2,4)(3,6)

julia> describe_involution(W,w)
1-element Vector{Int64}:
 3

julia> w==longest(reflection_subgroup(W,[3]))
true
```
For now only works for finite Coxeter groups.
"""
describe_involution(W,w)=simpleroots_subsystem(W::FC,
                                        filter(i->action(W,i,w)==i+W.N,1:W.N))

Base.length(W::FiniteCoxeterGroup,w)=count(i->isleftdescent(W,w,i),1:nref(W))

dimension(W::FiniteCoxeterGroup)=2*nref(W)+rank(W)
const dim=dimension
Base.length(W::FiniteCoxeterGroup)=prod(degrees(W))
@inline Base.parent(W::FiniteCoxeterGroup)=W
Base.in(w,W::FiniteCoxeterGroup)=w in W.G
PermRoot.nhyp(W::FiniteCoxeterGroup)=nref(W)

Base.:(==)(W::FiniteCoxeterGroup,W1::FiniteCoxeterGroup)=W.G==W1.G

# Weyl groups should  be  derived  from PermRoot and CoxeterGroup.
# Since  such inheritance  is impossible  in Julia,  we have  to choose. We
# derive  it  from  CoxeterGroup  and  represent  the  other inheritance by
# composition.  Thus, implementations  of FiniteCoxeterGroups  have a field
# .G, a PermRootGroup, to which we forward the following methods.
@forward FC.G PermRoot.action, PermRoot.cartan, 
 PermRoot.coroots, PermRoot.coxnum, PermRoot.inclusion, 
 PermRoot.inclusiongens, PermRoot.independent_roots, PermRoot.invariants, 
 PermRoot.invariant_form, PermRoot.YMatrix, PermRoot.PermX, PermRoot.PermY, 
 PermRoot.rank, PermRoot.roots, PermRoot.reflection_character,
 PermRoot.refleigen, PermRoot.reflrep, PermRoot.refltype, PermRoot.restriction,
 PermRoot.simplecoroots,
 PermRoot.simple_conjugating, PermRoot.simple_reps, PermRoot.simpleroots,
 PermRoot.unique_refls, PermRoot.torus_order, PermRoot.baseX, PermRoot.invbaseX,
 PermRoot.central_action, 
 Perms.reflection_length, Perms.last_moved

 PermRoot.refls(W::FC)=refls(W.G)

if false # slower now with changes in PermGroups
function Groups.transporting_elt(W::FC,H1::FC,H2::FC)
  if parent(W)!=parent(H1) || parent(W)!=parent(H1) error("not same parent") end
  if isomorphism_type(H1)!=isomorphism_type(H2) return end
  PH1=parabolic_closure(W,inclusiongens(H1,W))
  p1=standard_parabolic(W,PH1);PH1=sort(PH1.^p1);H1=H1^p1;
  PH2=parabolic_closure(W,inclusiongens(H2,W))
  p2=standard_parabolic(W,PH2);PH2=sort(PH2.^p2);H2=H2^p2;
  C=CoxGroups.parabolic_category(W,PH1)
  i=findfirst(==(PH2),C.obj)
  if isnothing(i) return end
  p=Garside.from(C,1,i)
  if isnothing(p) return end
  H1=H1^p
  if H1==H2 return p1*p*inv(p2) end
  q=transporting_elt(reflection_subgroup(W,PH2),
     sort(inclusiongens(H1)),sort(inclusiongens(H2)),onsets)
  if isnothing(q) return end
  p1*p*q*inv(p2)
end
else
Groups.transporting_elt(H::FC,R::FC,G::FC)=
  transporting_elt(H.G,sort(inclusiongens(R)),sort(inclusiongens(G)),onsets)
end

#--------------- FCG -----------------------------------------
function Base.show(io::IO, W::FCG)
  if haskey(W,:callname)
    if hasdecor(io) printTeX(io,W.TeXcallname)
    else print(io,W.callname)
    end
  elseif !hasdecor(io) && semisimplerank(W)==0
    print(io,"torus(",rank(W),")")  
  else show(io,W.G)
  end
end

function Base.show(io::IO,t::Type{FCG{T,T1}})where {T,T1}
  print(io,"FiniteCoxeterGroup{Perm{",T,"},",T1,"}")
end

@inline CoxGroups.nref(W::FCG)=W.N

# several functions on finite coxeter groups depend on the fact that
# isleftdescent is defined for all i in 1:nref(W), for example length and
# inversions
CoxGroups.isleftdescent(W::FCG,w,i::Integer)=i^w>W.N
CoxGroups.isrightdescent(W::FCG,w,i::Integer)=preimage(i,w)>W.N

"""
`coxeter_group(type,rank[,bond];sc=false)` (or `coxgroup`)

If `C=cartan(type,rank[,bond])`, this is equivalent to `rootdatum(C)`.
If `sc=true` this is equivalent to `rootdatum(permutedims(C),one(C))`.

The  resulting object `W`, a  `FiniteCoxeterGroup`, has an additional entry
compared to a `PermRootGroup`.

  - `W.rootdec`:  the root vectors, given  as linear combinations of simple
    roots.  The first `nref(W)` roots are  positive, the next `nref(W)` are
    the corresponding negative roots. Moreover, the first
    `semisimplerank(W)`  roots are the simple roots. The positive roots are
    ordered by increasing height.

and `roots(W)` is ordered is the same way as `W.rootdec`.

For  how to  get various  information on  the root  system and  the Coxeter
group,   see  the  functions   `nref,  coroots,  rootlengths,  simple_reps,
simple_conjugating,  reflrep,  simpleroots,  simplecoroots,  PermX, cartan,
inclusion, restriction, action, rank, semisimplerank`

In terms of root data, this function returns the adjoint root datum of Weyl
group  `W`.  If  `sc=true`  it  returns  the  simply  connected root datum.
```julia_repl
julia> W=coxgroup(:G,2)
G₂

julia> cartan(W)
2×2 Matrix{Int64}:
  2  -1
 -3   2

julia> W.rootdec
12-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [1, 1]
 [1, 2]
 [1, 3]
 [2, 3]
 [-1, 0]
 [0, -1]
 [-1, -1]
 [-1, -2]
 [-1, -3]
 [-2, -3]

julia> reflrep(W)
2-element Vector{Matrix{Int64}}:
 [-1 0; 1 1]
 [1 3; 0 -1]
```
"""
CoxGroups.coxeter_group(t::Symbol,r::Int=0,b::Int...;sc=false)=iszero(r) ? coxgroup() :
sc ? rootdatum(permutedims(cartan(t,r,b...)),Matrix{Int}(I,r,r)) : rootdatum(cartan(t,r,b...))

"""
`rootdatum(C::AbstractMatrix)`  adjoint root datum  from Cartan matrix `C`.
The adjoint group is also the default returned for
`coxeter_group(type,rank)`. The following methods all define `pgl₃`.
```julia-repl
julia> rootdatum(cartan(:A,3))==coxgroup(:A,3)
true

julia> rootdatum(:pgl,3)
pgl₃
```
"""
rootdatum(C::AbstractMatrix)=isempty(C) ? rootdatum(C,C) : rootdatum(one(C),C)

"""
`rootdatum(R::AbstractMatrix,CR::AbstractMatrix)`

root  datum from `R` whose  rows are the simple  roots on a basis of `X(T)`
and `CR` whose rows are the simple coroots on a basis of `Y(T)`. The following
methods all define `gl₃`.

```julia-repl
julia> rootdatum(:gl,3)==rootdatum("gl",3)
true

julia> rootdatum([1 -1 0;0 1 -1],[1 -1 0;0 1 -1])
A₂Φ₁
```
"""
function rootdatum(rr::AbstractMatrix,cr::AbstractMatrix)
  C=cr*transpose(rr) # Cartan matrix
  rootdec=roots(C) # difference with PermRootGroup is order of roots here
  N=length(rootdec)
  if isempty(rr) G=PRG(size(rr,2)) # a torus
  else
  r=Ref(transpose(rr)).*rootdec
  r=vcat(r,-r)
  rootdec=vcat(rootdec,-rootdec)
  mats=map(reflectionMatrix,eachrow(rr),eachrow(cr)) # reflrep
  # permutations of the roots effected by mats
# gens=map(M->sortPerm(Ref(transpose(M)).*r),mats).\sortPerm(r)
  gens=map(M->Perm(Vector{Perms.Idef}(indexin(Ref(transpose(M)).*r,r))),mats)
  coroots=Vector{Vector{eltype(cr)}}(undef,length(r))
  coroots[axes(cr,1)].=eachrow(cr)
  G=PRG(gens,one(gens[1]),mats,r,coroots,
        Dict{Symbol,Any}(:cartan=>C,:refltype=>type_cartan(C)))
  end
  FCG(G,rootdec,N,Dict{Symbol,Any}())
end


"""
`torus(rank::Integer)`

This  function returns the object corresponding to the notion of a torus of
dimension  `rank`, a Coxeter  group of semisimple  rank 0 and given `rank`.
This  corresponds to a split torus; the extension to Coxeter cosets is more
useful.

```julia-repl
julia> torus(3)
Φ₁³
```
"""
torus(i::Integer)=FCG(PRG(i),Vector{Int}[],0,Dict{Symbol,Any}())

"`istorus(W)` whether `W` is a torus"
istorus(W)=isempty(roots(W))

PermRoot.radical(W::FC)=torus(rank(W)-semisimplerank(W))

"`coxeter_group()` or `coxgroup()` the trivial Coxeter group"
CoxGroups.coxeter_group()=torus(0)

function uniinv(m) # inverse of unimodular matrix
  ch,c=GenLinearAlgebra.charpolyandcomatrix(m)
  (-1)^size(m,1)*ch[1]*c
end

"""
`derived_datum(G)`

The  derived datum of `(X,Φ,Y,Φ^∨)` is `(X/Φ^{∨⟂}, Φ, Y∩ ℚ Φ^∨, Φ^∨)` where
`⟂`  denotes the orthogonal. This function  starts with a root datum object
`G` and returns the root datum object corresponding to the derived datum.
"""
function derived_datum(G::FC)
  s=semisimplerank(G)
# Find a basis of Y where the last rank(G)-s columns of simplecoroots are 0
# Then  in  this  basis  Y∩  ℚ  Φ^∨  is  the ℤ-space spanned by the first s
# columns.  (gv^-t,gv) is a  (co)-isomorphism of root  data G->G1 and in G1
# projecting on the first s columns gives the derived datum.
  gv=col_hermite_transforms(simplecoroots(G)).coltrans
  g=uniinv(permutedims(gv))  # so gv = g^-t.
  rootdatum(simpleroots(G)*g[:,1:s],simplecoroots(G)*gv[:,1:s])
end

#reflection representation in the basis of rootdec
#function reflrep(W::FCG,w)
#  vcat(transpose(hcat(roots.(Ref(W),(1:ngens(W)).^w)...)))
#end

"""
`rootlengths(W::FiniteCoxeterGroup)`
 the vector  of the (squared)  length of the roots of `W`.
 The  shortest roots in an irreducible subsystem are given the length 1. In
 a Weyl group the others then have length 2 (or 3 in type `G₂`). The matrix
 of the `W`-invariant bilinear form is given by
 `Diagonal(rootlengths(W)[1:ngens(W)])*cartan(W)`.
"""
function rootlengths(W::FiniteCoxeterGroup)
  get!(W,:rootlengths)do
    C=cartan(W)
    lengths=fill(eltype(C)(1),2*nref(W))
    for t in refltype(W)
      I=t.indices
      if length(I)>1 && C[I[1],I[2]]!=C[I[2],I[1]]
        lengths[I[2]]=-C[I[1],I[2]]
        lengths[I[1]]=-C[I[2],I[1]]
      elseif length(I)>2 && C[I[2],I[3]]!=C[I[3],I[2]]
        lengths[I[3]]=-C[I[2],I[3]]
        lengths[I[1]]=-C[I[3],I[2]]
      end
    end
    for i in eachindex(lengths) lengths[i]=lengths[simple_reps(W,i)] end
    lengths
  end::Vector{eltype(cartan(W))}
end

rootlengths(W,i)=rootlengths(W)[i]
"""
`highest_short_root(W)`

It  is  an  error  if  `W`  is  not an irreducible Coxeter group. Otherwise
`HighestShortRoot`  returns the index  of the unique  short root of maximal
height  of `W`. If all roots have the same length then this is the index of
the unique root of maximal height, equal to `nref(W)`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> highest_short_root(W)
4
```
"""
function highest_short_root(W)
  if length(refltype(W))!=1 error("W should be irreducible and non-trivial") end
  findlast(==(1),rootlengths(W,1:nref(W)))
end

PermRoot.invariant_form(W::FiniteCoxeterGroup)=
  Diagonal(rootlengths(W,1:ngens(W)))*cartan(W)

function Base.:*(W1::FC,W2::FC)
  r=cat(simpleroots(parent(W1)),simpleroots(parent(W2)),dims=[1,2])
  cr=cat(simplecoroots(parent(W1)),simplecoroots(parent(W2)),dims=[1,2])
  res=rootdatum(r,cr)
  if haskey(W1,:callname) && haskey(W2,:callname)
    res.callname=W1.callname*"*"*W2.callname
    res.TeXcallname=W1.TeXcallname*"\\times "*W2.TeXcallname
  end
  if W1==parent(W1) && W2==parent(W2) return res end
  if !issubset(inclusiongens(W1),1:ngens(parent(W1)))
    error("not implemented") end
  if !issubset(inclusiongens(W2),1:ngens(parent(W2)))
    error("not implemented") end
  reflection_subgroup(res,vcat(inclusiongens(W1),
                               inclusiongens(W2).+ngens(parent(W1))))
end

Base.:*(W1::FC,W2::PermRootGroup)=W1.G*W2
Base.:*(W1::PermRootGroup,W2::FC)=W1*W2.G

function Base.:^(W::FC,p::Perm)
  WW=parent(W)
  if !(p in WW) error("can only conjugate in parent") end
  reflection_subgroup(WW,inclusiongens(W).^p)
end

#--------------- FCSG -----------------------------------------
Base.show(io::IO, W::FCSG)=show(io,W.G)

function Base.show(io::IO,t::Type{FCSG{T,T1}})where {T,T1}
  print(io,"FiniteCoxeterSubGroup{Perm{$T},$T1}")
end

CoxGroups.nref(W::FCSG)=W.N

Base.parent(W::FCSG)=W.parent

# if I are the positive roots or all roots of a subsystem find the simple ones
function simpleroots_subsystem(W,I)
  N=nref(W)
  filter(I) do i
    if i>N return false end
    r=refls(W,i)
    all(j->j>N || j==i || !isleftdescent(W,r,j),I)
  end
end

"""
`reflection_subgroup(W::FiniteCoxeterGroup,I)`

For `I⊆1:nref(W)`, the subgroup `H` of `W` generated by `refls(W,I)`.

A   theorem  found  independently  by  [Deodhar1989](biblio.htm#Deo89)  and
[Dyer1990](biblio.htm#Dye90)  is that  a subgroup  `H` of  a Coxeter system
`(W,S)`  generated by reflections  has a canonical  Coxeter generating set,
formed  of the `t  ∈ Ref(H)` such  `length(W,tt')>length(W,t)` for any `t'∈
Ref(H)`  different  from  `t`.  This  is  used  by `reflection_subgroup` to
determine the Coxeter system of `H`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> diagram(W)
O⇛ O G₂
1  2

julia> H=reflection_subgroup(W,[2,6])
G₂₍₂₆₎=Ã₁×A₁

julia> diagram(H)
O Ã₁
1
O A₁
2
```

The  notation `G₂₍₂₆₎`  means that  `roots(W,[2,6])` is  a system of simple
roots for `H`.

If  `H` is a  standard parabolic subgroup  of a Coxeter  group `W` then the
length  function on  `H` (with  respect to  its set  of generators)  is the
restriction  of the length function on `W`. This need not no longer be true
for  arbitrary reflection subgroups of  `W`. 
```julia-repl
julia> elH=word.(Ref(H),elements(H))
4-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [1, 2]

julia> elW=word.(Ref(W),elements(H))
4-element Vector{Vector{Int64}}:
 []
 [2]
 [1, 2, 1, 2, 1]
 [1, 2, 1, 2, 1, 2]

julia> map(w->H(w...),elH)==map(w->W(w...),elW)
true
```
We  implement finite  reflection groups  as permutation  groups on a set of
roots.  Consequently,  a  reflection  subgroup  `H⊆  W`  is  a  permutation
subgroup, thus its elements are represented as permutations of the roots of
the parent group.
"""
function PermRoot.reflection_subgroup(W::FCG,I::AbstractVector{<:Integer})
# contrary to GAP3, I is indices in W and not parent(W)
  if !all(<=(ngens(W)),I) 
    inclusion=sort!(vcat(orbits(refls(W,I),I)...))
    inclusion=inclusion[1:div(length(inclusion),2)]
    I=simpleroots_subsystem(W,inclusion)
  end
  if !haskey(W,:reflection_subgroups) 
    W.reflection_subgroups=Dict{Vector{Int},FCSG}()
  end
  get!(W.reflection_subgroups,I)do
  C=cartan(W,I)
  rootdec=isempty(C) ? empty(W.rootdec) : roots(C)
  rootdec=vcat(rootdec,-rootdec)
  if isempty(rootdec) inclusion=Int[]
  else m=transpose(toM(W.rootdec[I]))
    inclusion=map(rootdec)do r
      findfirst(==(m*r),W.rootdec)
    end
  end
  restriction=zeros(Int,2*W.N)
  restriction[inclusion]=1:length(inclusion)
  refltypes=map(type_cartan(C)) do t
    if (t.series in [:A,:D]) && rootlengths(W,inclusion[t.indices[1]])==1
      p=simple_reps(W)[inclusion[t.indices[1]]]
      for s in refltype(W)
        if p in s.indices && s.series in [:B,:C,:F,:G]
          t.short=true
        end
      end
    end
    t
  end
  prop=Dict{Symbol,Any}(:cartan=>C,:refltype=>refltypes)
  if isempty(inclusion) prop[:rank]=PermRoot.rank(W) end
  gens=isempty(I) ? eltype(W)[] : refls(W,I)
  G=PRSG(gens,one(W.G),inclusion,restriction,W.G,prop)
  FCSG(G,rootdec,div(length(inclusion),2),W,Dict{Symbol,Any}())
  end
end

PermRoot.reflection_subgroup(W::FCSG,I::AbstractVector{<:Integer})=
  reflection_subgroup(W.parent,inclusion(W)[I])

CoxGroups.isleftdescent(W::FCSG,w,i::Integer)=inclusion(W,i)^w>nref(parent(W))
# next is 25% slower
#CoxGroups.isleftdescent(W::FCSG,w,i::Integer)=action(W,i,w)>nref(W)
CoxGroups.isrightdescent(W::FCSG,w,i::Integer)=preimage(inclusion(W,i),w)>nref(parent(W))

function rootlengths(W::FCSG)
  get!(W,:rootlengths)do
    rootlengths(parent(W))[inclusion(W)]
  end
end

"""
`relative_group(W::FiniteCoxeterGroup,J)`

`J`  should be a if *distinguished* subset of `S==eachindex(gens(W))`, that
is if for `s∈ S-J` we set ``v(s,J)=w₀^{J∪ s}w₀ᴶ`` then `J` is stable by all
`v(s,J)`.  Then ``N_W(W_J)=W_J⋊ N₁``  where `N₁` is  the group generated by
the  `v(s,J)`,  which  form  a  Coxeter  system for `N₁`. Equivalently `N₁`
consists   of  the   `J`-reduced  elements   of  `N_W(W_J)`.  The  quotient
`R=N_W(W_J)/W_J` has a natural reflection representation on ``X(ZL_J/ZG)``,
using  that by [Lusztig1976](biblio.htm#Lus76), the  images of the roots of
`W`  in  ``X(ZL_J)``  form  a  root  system.  The function returns `R` as a
reflection  group on ``X(ZL_J/ZG)``, with  some extra attributes reflecting
its origin
  - `R.relative_indices=setdiff(S,J)` in a certain order
  - `R.toparent=` the list of `v(s,J)` corresponding to `.relative_indices`;
     defines an isomorphism `R→ N₁`.
  - `R.fromparent`  is  a  function  mapping elements of `N₁` to `R`. The 
    inverse mapping to `.toparent`.
"""
function relative_group(W::FC,J::Vector{<:Integer})
  S=eachindex(gens(W))
  if !issubset(J,S)
    error("implemented only for standard parabolic subgroups")
  end
  I=setdiff(S,J) # keeps the order of S
  sort!(I,by=x->findfirst(t->x in indices(t),refltype(W)))
  # order I compatibly with components of W
  vI=map(i->longest(W,vcat([i],J))*longest(W,J),I)
  if any(g->sort(action.(Ref(W),J,g))!=sort(J),vI)
    error("$J is not distinguished in $W")
  end
  qr=i->W.rootdec[i][I]
  res=isempty(I) ? coxgroup() :
    rootdatum(improve_type([ratio(qr(j)-qr(action(W,j,vI[ni])),qr(i))
                            for (ni,i) in pairs(I), j in I]))
  res.relative_indices=I
  res.toparent=vI
  res.fromparent=
    function(w)c=Perm()
      while true
        r=findfirst(x->isleftdescent(W,w,x),I)
        if isnothing(r) return c end
 	w=vI[r]*w
        c*=res(r)
      end
    end
  res
end

function Groups.position_class(W::FiniteCoxeterGroup,w)
  l=PermGroups.positions_class(W.G,w)
  if length(l)==1 return only(l) end
  ww=word(W,w)
  iw=map(refltype(W))do t
    v=map(i->findfirst(==(i),t.indices),filter(i->i in t.indices,ww))
    getchev(t,:ClassParameter,v)
  end
  if any(isnothing,iw) return position_class(W.G,w) end
  findfirst(==(iw),map(x->x.param,conjugacy_classes(W)))
end

end
