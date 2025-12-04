"""
# Reductive algebraic groups and root data

Let  us fix an  algebraically closed field  `K` and let  `𝐆` be a connected
reductive  algebraic group over `K`. Let `𝐓` be a maximal torus of `𝐆`, let
`X(𝐓)`  be the  character group  of `𝐓`  (resp. `Y(𝐓)`  the dual lattice of
one-parameter  subgroups  of  `𝐓`)  and  `Φ`  (resp  `Φ^`) the roots (resp.
coroots) of `𝐆` with respect to `𝐓`.

Then  `𝐆` is  determined up  to isomorphism  by the  *root datum* `(X(𝐓),Φ,
Y(𝐓),Φ^)`.  In algebraic terms, this consists  in giving a free `ℤ`-lattice
`X=X(𝐓)` of dimension the *rank* of `𝐓` (which is also called the *rank* of
`𝐆`),  and a root system `Φ ⊂ X`,  and similarly giving the dual roots `Φ^⊂
Y=Y(𝐓)`.

This  can be  constructed by  giving integral  matrices as arguments to our
function  [`rootdatum`](@ref). We  assume that  the canonical  basis of the
vector  space `V` on which `W` acts is  a `ℤ`-basis of `X`, and the rows of
the first argument [`simpleroots`](@ref)`(W)` of `rootdatum` are the simple
roots  `Π` expressed in this basis of `X`. Similarly the rows of the second
argument  [`simplecoroots`](@ref)`(W)` describes the  simple coroots in the
basis  of `Y` dual to the chosen  basis of `X`. The duality pairing between
`X`  and `Y` is the canonical one,  that is the pairing between vectors `x∈
X` and `y∈ Y` is given by `transpose(x)*y`. Thus, we must have the relation
`simplecoroots(W)*permutedims(simpleroots(W))=cartan(W)`.

The  roots need not generate  `V`, so the matrices  need not be square. For
example, the root datum of the linear group of rank 3 can be obtained as:

```julia-repl
julia> W=rootdatum([-1 1 0;0 -1 1],[-1 1 0;0 -1 1])
A₂Φ₁

julia> reflrep(W,W(1))
3×3 Matrix{Int64}:
 0  1  0
 1  0  0
 0  0  1
```
here  the  Weyl  group  is  the  symmetric  group  on  3  letters acting by
permutation of the basis of `X`. The dimension of `X`,
`size(simpleroots(W),2)`,  is the [`rank`](@ref)`(W)`] and the dimension of
the  subspace  generated  by  the  roots,  `size(simpleroots(W),1)`, is the
[`semisimplerank`](@ref)`(W)`.  In  the  example  the  rank  is  3  and the
semisimple rank is 2.

The  default representation  obtained as  `W=coxgroup(:A,2)` corresponds to
the  adjoint algebraic group (the group in its isogeny class with a trivial
centre).  In this case  `Φ` is a  basis of `X`,  so `simpleroots(W)` is the
identity  matrix and `simplecoroots(W)` is the Cartan matrix `cartan(W)` of
the root system.

The   form  `coxgroup(:A,2,sc=true)`   constructs  the   semisimple  simply
connected  algebraic  group,  where  `simpleroots(W)`  is the transposed of
`cartan(W)` and `simplecoroots(W)` is the identity matrix.

An  extreme form of root data can  be specified more conveniently: when `W`
is  the trivial `coxgroup()` (in this case `𝐆 ` is a torus), the root datum
has  no roots, thus  is entirely determined  by the rank  `r`. The function
[`torus`](@ref)`(r)`  constructs  such  a  root  datum  (it  could  be also
obtained by giving two `0×r` matrices to `rootdatum`).

Finally,  the `rootdatum` function also understands some familiar names for
the  algebraic groups for which it gives the results that could be obtained
by giving the appropriate `simpleroots(W)` and `simplecoroots(W)` matrices:
```julia-repl
julia> rootdatum(:gl,3)   # same as the previous example
gl₃
```
The types of root data which are understood are
 `:gl, :pgl, :sl, :slmod, :tgl :sp, :csp, :psp, :so, :pso, :cso, :halfspin, 
  :gpin, :spin, :E6, :E6sc, :CE6, :E7, :E7sc, :CE7, :E8, :F4, :G2`.

The  group `𝐆` is *semisimple* if the rank is equal to the semisimple rank.
In  this case, things  are constrained: the  lattice `X` having an integral
pairing  withe the coroots,  is in the  dual lattice of  the lattice of the
coroots.  This dual lattice is called  the *weight lattice*. The dual basis
to  the simple  coroots is  called the  *fundamental weights*. Similarly we
have the *coweight lattice* and the *fundamental coweights*. In general the
lattice  `X` is  an intermediate  lattice between  the root  and the weight
lattices.

It  follows  that  the  finite  abelian  group  obtained by quotienting the
coweight  lattice by the  coroot lattice describes  all possibilities. This
group   is  called  the  *fundamental  group*   of  the  root  system;  the
*fundamental  group* of `𝐆` is  the quotient of `Y`  by the coroot lattice.
The  fundamental  group  of  the  root  system  is also called the *adjoint
fundamental  group* since for an adjoint group `Y` is the coweight lattice.
The  fundamental  group  of  `𝐆`  is  a subgroup of the adjoint fundamental
group.   It  turns  out  that  the  non-trivial  elements  of  the  adjoint
fundamental  group are in bijection with the *minuscule coweights*, that is
the  coweights whose pairing with every root is in `-1,0,1` (this bijection
is  compatible  with  the  group  structure,  taking  coweights  modulo the
coroots).  The  function  [`intermediate_group`](@ref)  can  construct  any
semisimple group using this.

The fundamental group has another incarnation which can be more convenient.

  - The *extended simple roots* is `Π̃=Π∪{-α₀}` where `α₀` is the highest root.

Via  the theory  of the  affine Weyl  group (see `affine`), the fundamental
group  is  isomorphic  to  the  subgroup  of  `W` which stabilizes `Π̃`. We
represent  `Π̃` by the  indices of the  roots it contains  and the function
[`fundamental_group`](@ref)  returns it as a group of permutations of these
indices.
"""
module Rootdata
using ..Chevie
export torus, istorus, derived_datum, SubTorus, weightinfo, fundamental_group,
       intermediate_group, isomorphism_type, weights, coweights
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
torus(i::Integer)=Weyl.FCG(PRG(i),Vector{Int}[],0,Dict{Symbol,Any}())

"`istorus(W)` whether `W` is a torus"
istorus(W)=isempty(roots(W))

PermRoot.radical(W::Weyl.FC)=torus(rank(W)-semisimplerank(W))

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
function derived_datum(G::Weyl.FC)
  s=semisimplerank(G)
# Find a basis of Y where the last rank(G)-s columns of simplecoroots are 0
# Then  in  this  basis  Y∩  ℚ  Φ^∨  is  the ℤ-space spanned by the first s
# columns.  (gv^-t,gv) is a  (co)-isomorphism of root  data G->G1 and in G1
# projecting on the first s columns gives the derived datum.
  gv=col_hermite_transforms(simplecoroots(G)).coltrans
  g=uniinv(permutedims(gv))  # so gv = g^-t.
  rootdatum(simpleroots(G)*g[:,1:s],simplecoroots(G)*gv[:,1:s])
end

struct SubTorus{T<:Integer}
  gens::Matrix{T}
  complement::Matrix{T}
  group
end

"""
`SubTorus(W,Y::Matrix{<:Integer})`

The  function  returns  the  subtorus  𝐒  of  the  maximal torus `𝐓` of the
reductive  group represented by the Weyl group  `W` such that `Y(𝐒)` is the
(pure)  sublattice of  `Y(𝐓)` generated  by the  (integral) vectors  `Y`. A
basis  of `Y(𝐒)` in  Hermite normal form  (for easy memebership testing) is
computed  and stored in the field `S.gens` of the returned SubTorus struct.
A  basis of `Y(𝐓)` adapted to `Y(𝐒)` is  also computed. This means a set of
integral   vectors,  stored  in  `S.complement`,   is  computed  such  that
`M=vcat(S.gens,S.complement)`   is   a   basis   of   `Y(𝐓)`  (equivalently
`M∈GL(Z^{rank(W)})`.  An  error  is  raised  if  `Y` does not define a pure
sublattice.

```julia-repl
julia> W=coxgroup(:A,4)
A₄

julia> SubTorus(W,[1 2 3 4;2 3 4 1;3 4 1 1])
SubTorus(A₄,[1 0 3 -13; 0 1 2 7; 0 0 4 -3])

julia> SubTorus(W,[1 2 3 4;2 3 4 1;3 4 1 2])
ERROR: not a pure sublattice
Stacktrace:
 [1] error(::String) at ./error.jl:33
 [2] Chevie.Weyl.SubTorus(::FiniteCoxeterGroup{Perm{Int16},Int64}, ::Matrix{Int64}) at /home/jmichel/julia/Chevie.jl/src/Weyl.jl:1082
 [3] top-level scope at REPL[25]:1
```
"""
function SubTorus(W,V::Matrix{<:Integer}=reflrep(W,one(W)))
  C=complementInt(V)
  if any(!=(1),C.moduli) error("not a pure sublattice") end
  SubTorus(baseInt(C.sub),C.complement,W)
end

Base.show(io::IO,T::SubTorus)=print(io,"SubTorus(",T.group,",",T.gens,")")

Chevie.rank(T::SubTorus)=size(T.gens,1)

#----------------------- Fundamental group --------------------------
function WeightToAdjointFundamentalGroupElement(W,l::Vector;full=false)
  if isempty(l) return Perm();end
  prod(x->WeightToAdjointFundamentalGroupElement(W,x;full),l)
end

function WeightToAdjointFundamentalGroupElement(W,i;full=false)
  t=refltype(W)[findfirst(t->i in t.indices,refltype(W))]
  l=copy(t.indices)
  b=longest(W,l)*longest(W,setdiff(l,[i]))
  push!(l,maximum(findall(
    i->all(j->j in t.indices || W.rootdec[i][j]==0,1:semisimplerank(W)),
  eachindex(W.rootdec))))
  full ? b : restricted(b,inclusion.(Ref(W),l))
end

"""
`weightinfo(W)`

`W`  is a  Coxeter group  record describing  an algebraic  group `𝐆 `, or a
`IypeIrred`. The function is independent of the isogeny type of `𝐆`(so just
depends  on `refltype(W)`, that is  on the root system).  It returns a dict
with the following keys:

  - `:minusculeWeights` the minuscule weights, described as their position 
    in  the  list  of  fundamental  weights.  For non-irreducible groups, a
    weight  is the sum of at most one weight in each irreducible component.
    It  is represented as  the list of  its components. For consistency, in
    the  case  of  an  irreducible  system,  a  weight  is represented as a
    one-element list.

  - `:minusculeCoweights` the minuscule coweights, represented in the same
    manner as the minuscule weights

  - `:decompositions` for each coweight, its decomposition in terms of the
    generators  of the adjoint fundamental group  (given by the list of the
    exponents  of the generators). Together with  the next field it enables
    to work out the group structure of the adjoint fundamental group.

  - `:moduli` the list of orders of the generators of the fundamental group.

  - `:AdjointFundamentalGroup` the list of generators of the adjoint
    fundamental group, given as permutations of the extended root basis.

  - `:CenterSimplyConnected` A list of semisimple elements generating the
    center of the universal covering of  𝐆

  - `:chosenAdaptedBasis` A basis  of the  weight lattice  adapted to the
    root lattice. In the basis of the fundamental weights, the root lattice
    is  given  by  `C=transpose(cartan(W))`.  The  property  specifying the
    `:chosenAdaptedBasis` is that the Hermite normal form of
    `C*chosenAdaptedBasis`  is almost in Smith  normal form (it is diagonal
    but  the diagonal entries may be  permuted compared to the Smith normal
    form; the non-trivial entries are in the positions corresponding to the
    generators of the fundamental group as indicated by `:decompositions`).

```julia-repl
julia> weightinfo(coxgroup(:A,2)*coxgroup(:B,2))
Dict{Symbol, Array} with 8 entries:
  :moduli                  => [3, 2]
  :minusculeWeights        => [[1, 3], [1], [2, 3], [2], [3]]
  :decompositions          => [[2, 1], [2, 0], [1, 1], [1, 0], [0, 1]]
  :highestroot             => [5, 7]
  :chosenAdaptedBasis      => [1 -1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
  :minusculeCoweights      => [[1, 4], [1], [2, 4], [2], [4]]
  :CenterSimplyConnected   => Vector{Rational{Int64}}[[1//3, 2//3, 0, 0], [0, 0…
  :AdjointFundamentalGroup => [(1,12,2), (4,14)]
```
"""
function weightinfo(t::TypeIrred)
  r=chevieget(t,:WeightInfo)
  if isnothing(r)
    r=Dict{Symbol,Any}(:moduli=>Int[],:decompositions=>Vector{Vector{Int}}[],
         :minusculeWeights=>Vector{Int}[])
  end
  if !haskey(r,:minusculeCoweights)
    r[:minusculeCoweights]=r[:minusculeWeights]
  end
  if !haskey(r,:chosenAdaptedBasis)
    r[:chosenAdaptedBasis]=Matrix{Int}(I,rank(t),rank(t))
  end
  r
end

function weightinfo(W)
  get!(W,:weightinfo)do
  M=Matrix{Int}(I,semisimplerank(W),semisimplerank(W))
  if isempty(refltype(W)) return Dict(:minusculeWeights=>Vector{Int}[],
         :minusculeCoweights=>Vector{Int}[],
         :decompositions=>Vector{Vector{Int}}[],
         :moduli=>Int[],
         :highestroot=>Int[],
         :CenterSimplyConnected=>Vector{Rational{Int}}[],
         :AdjointFundamentalGroup=>eltype(W)[]
        )
  end
  l=map(refltype(W)) do t
    r=weightinfo(t)
    if isempty(r[:moduli]) g=Int[]
      r[:ww]=eltype(W)[]
    else g=filter(i->sum(r[:decompositions][i])==1,
          eachindex(r[:minusculeCoweights])) # generators of fundamental group
      r[:ww]=map(x->WeightToAdjointFundamentalGroupElement(W,x),
               t.indices[r[:minusculeCoweights][g]])
    end
    r[:csi]=zeros(Rational{Int},length(g),semisimplerank(W))
    if !isempty(r[:moduli])
      C=modZ.(inv(Rational.(cartan(t))))
      r[:csi][:,t.indices]=C[r[:minusculeCoweights][g],:]
      r[:minusculeWeights]=t.indices[r[:minusculeWeights]]
      r[:minusculeCoweights]=t.indices[r[:minusculeCoweights]]
      M[t.indices,t.indices]=r[:chosenAdaptedBasis]
    end
    r[:csi]=Array.(toL(r[:csi]))
    r[:highestroot]=last(findmax(x->sum(x[indices(t)]),W.rootdec))
    r
  end
  res=Dict(:minusculeWeights=>cartesian(map(
                                   x->vcat(x[:minusculeWeights],[0]),l)...),
    :minusculeCoweights=>cartesian(map(
                                   x->vcat(x[:minusculeCoweights],[0]),l)...),
    :decompositions=>splat(vcat).(tcartesian(map(x->vcat(x[:decompositions],
                                 [0 .*x[:moduli]]),l)...)),
    :moduli=>reduce(vcat,map(x->x[:moduli],l)),
# center of simply connected group: the generating minuscule coweights
# mod the root lattice
    :CenterSimplyConnected=>reduce(vcat,getindex.(l,:csi)),
    :AdjointFundamentalGroup=>reduce(vcat,getindex.(l,:ww)),
    :chosenAdaptedBasis=>M,
    :highestroot=>map(x->x[:highestroot],l))
  n=1:length(res[:decompositions])-1
  res[:minusculeWeights]=map(x->filter(!iszero,x),res[:minusculeWeights][n])
  res[:minusculeCoweights]=map(x->filter(!iszero,x),res[:minusculeCoweights][n])
  res[:decompositions]=res[:decompositions][n]
  res
  end
end

" `weights(W)` fundamental weights in the basis of X(T)"
weights(W)=permutedims(inv(Rational.(cartan(W))))*simpleroots(W)

" `coweights(W)` fundamental coweights in the basis of Y(T)"
coweights(W)=inv(Rational.(cartan(W)))*simplecoroots(W)

"""
`fundamental_group(W;full=false)`

This  function returns the fundamental group of the algebraic group defined
by  the root datum  `W`, that is  the subgroup of  `W` which stabilises the
extended root basis (the set of simple roots enriched by the lowest root of
each  irreducible component). Unless the keyword `full=true` is given, this
is  returned as a group of permutations of  the indices of the roots in the
extended root basis.

The  definition  we  take  of  the  fundamental group of a (not necessarily
semisimple)  reductive group is (P∩ Y(𝐓))/Q where P is the coweight lattice
(the  dual lattice  in Y(𝐓)⊗  ℚ of  the root  lattice) and  Q is the coroot
latice.  The  bijection  between  elements  of  P/Q and permutations of the
extended  root basis is explained in  the context of non-irreducible groups
for example in [bon05; §3.B](@cite).
```julia-repl
julia> W=coxgroup(:A,3) # the adjoint group
A₃

julia> fundamental_group(W) # the 12th root is the lowest one
Group((1,2,3,12))

julia> W=rootdatum(:sl,4) # the semisimple simply connected group
sl₄

julia> fundamental_group(W)
Group(Perm{Int16}[])
```
"""
function fundamental_group(W;full=false)
  if istorus(W) return Group(Perm()) end
  e=weightinfo(W)[:minusculeCoweights]
  e=filter(x->all(isinteger,sum(coweights(W)[x,:];dims=1)),e) # minusc. coweights in Y
  if isempty(e) return Group(Perm()) end
  e=map(x->WeightToAdjointFundamentalGroupElement(W,x;full),e)
  Group(abelian_gens(e))
end

"""
`intermediate_group(W, indices...)`

returns a `rootdatum` which represents a semisimple algebraic group between
the   adjoint  group,  obtained  by  a  call  like  `rootdatum(:pgl,4)`  or
`coxgroup(:A,3)`, and the simply connected semi-simple group, obtained by a
call  like  `rootdatum(:sl,4)`  or  `coxgroup(:A,3;sc=true)`. This group is
specified  by a  subset of  the *minuscule  coweights*, the coweights whose
scalar  product with each root is in  `-1,0,1` (a coweight is an element of
the  *coweight lattice*  `Pᵛ`, the  lattice dual  to the root lattice). The
non-trivial  characters of the  (algebraic) center of  a semi-simple simply
connected  algebraic group are  in bijection with  the minuscule coweights,
and also in bijection with `Pᵛ/Qᵛ` where `Qᵛ` is the coroot lattice. If `W`
is  irreducible,  the  minuscule  coweights  are  part  of the basis of the
coweight  lattice given by  the *fundamental coweights*,  which is the dual
basis  of the simple roots. They can  therefore be specified by an index in
the  Dynkin diagram.  The `indices`  of minuscule  coweights in  the dynkin
diagram are indices where the coefficient of the highest root on the simple
roots  is  `1`.  The  constructed  intermediate  group  has  lattice `Y(𝐓)`
generated  by the root lattice and the given minuscule coweights. If `W` is
not  irreducible, a minuscule  coweight is a  sum of minuscule coweights in
different  components;  an  element  of  `indices`  is  thus itself a list,
interpreted as representing the sum of the corresponding coweights.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> fundamental_group(intermediate_group(W)) # simply connected
Group(Perm{Int16}[])

julia> fundamental_group(intermediate_group(W,2)) # intermediate
Group((1,3)(2,12))

julia> fundamental_group(intermediate_group(W,1)) # adjoint
Group((1,2,3,12))
```
"""
function intermediate_group(W,I...)
  C=cartan(W)
  w=inv(Rational.(C)); # w = coweights in terms of coroots
  R=one(w)
  for v in I
    if v isa Integer R=vcat(R,w[[v],:])
    else R=vcat(R,sum(w[v,:],dims=1))
    end
  end
  d=lcm(denominator.(R))
  R=baseInt(Integer.(d*R))//d
  rootdatum(Integer.(transpose(R*C)),Integer.(R^-1))
end

isomorphism_type(t::TypeIrred)=repr(t;context=(:TeX=>true))

function isomorphism_type(W;torus=false,TeX=false,limit=false)
  t=reverse(tally(map(x->fromTeX(isomorphism_type(x);TeX,limit),refltype(W))))
  l=map(x-> x[2]==1 ? x[1] : string(x[2],x[1]),t)
  d=rank(W)-semisimplerank(W)
  if d>0 && torus push!(l,fromTeX("T_{$d}";TeX,limit)) end
  join(l,"+")
end

end
