"""
Let  `R` be a  root system in  the real vector  space `V`. We say that `F₀∈
GL(V)`  is an  *automorphism of  `R`* if  it permutes  `R` and is of finite
order  (finite  order  is  automatic  if  `R` generates `V`). It follows by
cite[chap.  VI, S1.1, lemme 1]{Bou68} that  the dual `F₀*∈ GL(V*)` permutes
the  coroots  `R*⊂  V*`;  thus  `F₀`  normalizes  the  reflection group `W`
associated  to `R`, that is `w↦ F₀wF₀⁻¹` is an automorphism of `W`. Thus we
get a reflection coset `WF₀`, which we call a *Coxeter coset*.

The  motivation for introducing Coxeter  cosets comes from automorphisms of
algebraic  reductive groups, giving rise to non-split reductive groups over
finite fields. Let `𝐆` be a connected reductive algebraic group `𝐆` over an
algebraic  closure `𝔽̄_q` of a finite field `𝔽_q`, defined over `𝔽_q`; this
corresponds  to a  Frobenius endomorphism  `F` so  that the finite group of
rational  points `𝐆(𝔽_q)` identifies to the  subgroup `𝐆^F` of fixed points
under `F`.

Let `𝐓` be a maximal torus of `𝐆`, and `Φ` (resp. `Φ*`) be the roots (resp.
coroots)  of `𝐆` with respect  to `𝐓` in the  character group `X(𝐓)` (resp.
the  group of one-parameter subgroups `Y(𝐓)`). Then `𝐆` is determined up to
isomorphism  by `(X(𝐓),Φ,Y(𝐓),Φ*)`; this corresponds  to give a root system
in   the  vector  space  `V=ℚ ⊗ X(𝐓)`   and  a  rational  reflection  group
`W=N_𝐆(𝐓)/𝐓` acting on it.

If  `𝐓` is `F`-stable the Frobenius endomorphism `F` acts also naturally on
`X(T)`  and defines thus  an endomorphism of  `V`, which is  of the form `q
F₀`, where `F₀∈ GL(V)` is of finite order and normalizes `W`. We get thus a
Coxeter  coset `WF₀⊂GL(V)`.  The data  `(X(𝐓), Φ,  Y(𝐓), Φ*,  F₀)`, and the
integer   `q`  completely  determine  up   to  isomorphism  the  associated
*reductive finite group* `𝐆^F`. Thus these data is a way of representing in
the  essential  information  which  determines  a  finite  reductive group.
Indeed, all properties of Chevalley groups can be computed from that datum:
symbols  representing characters, conjugacy classes,  and finally the whole
character table of `𝐆^F`.

It  turns out that  many interesting objects  attached to this datum depend
only on `(V,W, F₀)`: the order of the maximal tori, the *fake degrees*, the
order  of `𝐆^F`, symbols representing unipotent characters, Deligne-Lusztig
induction  in  terms  of  *almost  characters*, the Fourier matrix relating
characters  and almost  characters, etc…  (see, e.g.,  cite{BMM93}). It is
thus  possible to extend their  construction to non-crystallographic groups
(or  even to more general complex  reflection groups, see "Spets"); this is
why  we did  not include  a root  system in  the definition of a reflection
coset. However, unipotent conjugacy classes for instance depend on the root
system.

We assume now that `𝐓` is contained in an `F`-stable Borel subgroup of `𝐆`.
This  defines an order  on the roots,  and there is  a unique element `ϕ∈ W
F₀`,  the  *reduced  element*  of  the  coset,  which  preserves the set of
positive  roots.  It  thus  defines  a  *diagram  automorphism*, that is an
automorphism  of the Coxeter system `(W,S)`.  This element is stored in the
component  '.phi' of the coset record. It may be defined without mentioning
the  roots,  as  follows:  `(W,F₀(S))`  is  another  Coxeter  system,  thus
conjugate to `S` by a unique element of `W`, thus there is a unique element
`ϕ∈ WF₀` which stabilizes `S` (a proof follows from cite[Theoreme 1, chap.
V,  S  3]{Bou68}).  We  consider  thus  cosets  of the form `Wϕ` where `ϕ`
stabilizes  `S`. The coset  `W ϕ` is  completely defined by the permutation
'.phi'  when `𝐆` is semi-simple --- equivalently when `Φ` generates `V`; in
this case we just need to specify 'phi' to define the coset.

There is a slight generalisation of the above setup, covering in particular
the  case of the Ree  and Suzuki groups. We  consider `𝐆^F` where `F` not a
Frobenius  endomorphism, but  an isogeny  such that  some power  `F^n` is a
Frobenius endomorphism. Then `F` still defines an endomorphism of `V` which
normalizes  `W`; we define a real number `q` such that `F^n` is attached to
an  `𝔽_{qⁿ}`-structure. Then we still have `F=q F₀` where `F₀` is of finite
order  but `q` is no more an integer.  Thus `F₀∈ GL(V⊗ ℝ)` but `F₀∉ GL(V)`.
For  instance, for the  Ree and Suzuki  groups, `F₀` is  an automorphism of
order  `2` of `W`, which is of type `G₂`, `B₂` or `F₄`, and `q=√2` for `B₂`
and  `F₄` and `q=√3`  for `G₂` This  can be constructed  starting from root
systems  for `G₂`, `B₂` or  `F₄` where all the  roots have the same length.
This kind of root system is *not* crystallographic. Such
non-crystallographic  root systems exist for all finite Coxeter groups such
as  the dihedral groups, `H₃` and `H₄`. We will call here *Weyl cosets* the
cosets  corresponding to rational forms  of algebraic groups, which include
thus some non-rational roots systems for `B₂`, `G₂` and `F₄`.

Conjugacy  classes and irreducible characters of Coxeter cosets are defined
as  for  general  reflection  cosets.  For  irreducible  characters of Weyl
cosets,  we choose (following Lusztig) for each `ϕ`-stable character of `W`
a  particular extension to a character of  `W⋊ ⟨ϕ⟩`, which we will call the
*preferred extension*. The *character table* of the coset `Wϕ` is the table
of  the restrictions to  `Wϕ` of the  preferred extensions. The question of
finding the conjugacy classes and character table of a Coxeter coset can be
reduced to the case of irreducible root systems `R`.

  *   The automorphism `ϕ` permutes the irreducible components of `W`, and
    `Wϕ`  is a direct  product of cosets  where `ϕ` permutes cyclically the
    irreducible components of `W`. The preferred extension is defined to be
    the  direct  product  of  the  preferred  extension  in  each  of these
    situations.

*  Assume now that `Wϕ` is a  descent of scalars, that is the decomposition
    in irreducible components `W=W₁× ⋯ × Wₖ` is cyclically permuted by `ϕ`.
    Then there are natural bijections from the `ϕ`-conjugacy classes of `W`
    to  the `ϕᵏ`-conjugacy classes  of `W₁` as  well as from the `ϕ`-stable
    characters  of `W` to the `ϕᵏ`-stable  characters of `W₁`, which reduce
    the  definition of preferred  extensions on `Wϕ`  to the definition for
    `W₁ϕᵏ`.

  *   Assume now  that `W`  is the  Coxeter group  of an  irreducible root
    system.   `ϕ`  permutes  the  simple   roots,  hence  induces  a  graph
    automorphism  on  the  corresponding  Dynkin  diagram.  If  `ϕ=1`  then
    conjugacy  classes and  characters coincide  with those  of the Coxeter
    group `W`.

The  nontrivial cases for crystallographic roots  systems are (the order of
`ϕ`  is written as left exponent to  the type): `²Aₙ`, `²Dₙ`, `³D₄`, `²E₆`.
For  non-crystallographic root  systems where  all the  roots have the same
length the additional cases `²B₂`, `²G₂`, `²F₄` and `²I₂(k)` arise.

  *   In  case  `³D₄`  the  group  `W⋊ ⟨ϕ⟩`  can be embedded into the
    Coxeter  group of type `F₄`, which induces a labeling for the conjugacy
    classes of the coset. The preferred extension is chosen as the (single)
    extension with rational values.

  *   In case  `²Dₙ` the  group `W⋊ ⟨ϕ⟩`  is isomorphic  to a Coxeter
    group of type `Bₙ`. This induces a canonical labeling for the conjugacy
    classes  of the coset and allows to define the preferred extension in a
    combinatorial  way  using  the  labels  (pairs  of  partitions) for the
    characters of the Coxeter group of type `Bₙ`.

  *  In the remaining crystallographic cases `ϕ` identifies to `-w₀` where
    `w₀`  is the longest element of `W`.  So, there is a canonical labeling
    of  the conjugacy classes and characters of  the coset by those of `W`.
    The  preferred extensions  are defined  by describing  the signs of the
    character values on `-w₀`.

The  most general  construction of  a Coxeter  coset is  by starting from a
Coxeter   datum   specified   by   the   matrices   of   'simpleRoots'  and
'simpleCoroots',  and  giving  in  addition  the  matrix 'F0Mat' of the map
`F₀:V→ V` (see the commands  'CoxeterCoset' and 'CoxeterSubCoset'). As for
Coxeter  groups,  the  elements  of  `Wϕ`  are  uniquely  determined by the
permutation  they  induce  on  the  set  of  roots  `R`.  We consider these
permutations as 'Elements' of the Coxeter coset.

Coxeter  cosets are implemented by a struct which points to a Coxeter datum
record  and  has  additional  fields  holding 'F0Mat' and the corresponding
element  'phi'. Functions on the coset (for example, 'classinfo') are about
properties  of  the  group  coset  `W  ϕ`  ;  however, most definitions for
elements of untwisted Coxeter groups apply without change to elements in `W
ϕ`:  e.g., if we define the length of  an element `wϕ∈ Wϕ` as the number of
positive  roots it sends to negative ones, it  is the same as the length of
`w`,  i.e., `ϕ` is of length `0`, since `ϕ` has been chosen to preserve the
set of positive roots. Similarly, the 'Coxeter word' describing `wϕ` is the
same as the one for `w`, etc…

We associate to a Coxeter coset `Wϕ` a *twisted Dynkin diagram*, consisting
of  the Dynkin diagram of `W` and  the graph automorphism induced by `ϕ` on
this  diagram (this specifies the  group `W⋊ ⟨F⟩`, mentioned above, up
to  isomorphism). See the  functions 'ReflectionType', 'ReflectionName' and
'Diagram' for Coxeter cosets.

Below  is an example showing first how to *not* define, then how to define,
the Weyl coset for a Suzuki group:

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> spets(W,Perm(1,2))
ERROR: matrix F must preserve the roots
Stacktrace:
 [1] error(::String) at ./error.jl:33
 [2] spets(::Gapjm.Weyl.FCG{Int16,Int64,PRG{Int64,Int16}}, ::Array{Int64,2}) at /home/jmichel/julia/Gapjm/src/Cosets.jl:241 (repeats 2 times)
 [3] top-level scope at REPL[19]:1

julia> W=coxgroup(:Bsym,2)
Bsym₂

julia> spets(W,Perm(1,2))
²Bsym₂

julia> CharTable(W)
CharTable(H(G(2,1,2)))
   │11. 1.1 .11 2. .2
───┼──────────────────
11.│  1   1   1 -1 -1
1.1│  2   .  -2  .  .
.11│  1  -1   1 -1  1
2. │  1   1   1  1  1
.2 │  1  -1   1  1 -1
```

A *subcoset* `Hwϕ` of `Wϕ` is given by a reflection subgroup `H` of `W` and
an  element `w` of `W`  such that `wϕ` induces  an automorphism of the root
system of `H`. For algebraic groups, this corresponds to a rational form of
a  reductive subgroup of maximal rank.  For example, if `Wϕ` corresponds to
the  algebraic group `𝐆` and  `H` is the trivial  subgroup, the coset `Hwϕ`
corresponds to a maximal torus `𝐓_w` of type `w`.

```julia-repl
julia> W=coxgroup(:Bsym,2)
Bsym₂

julia> WF=spets(W,Perm(1,2))
²Bsym₂

julia> subspets(WF,Int[],W(1))
Bsym₂₍₎=.Φ‴₈
```

A subgroup `H` which is a parabolic subgroup corresponds to a rational form
of  a Levi  subgroup of  `𝐆`. The  command 'twistings'  gives all rational
forms of such a Levi.

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> twistings(W,[1])
2-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterSubGroup{Perm{Int16},Int64}},1}:
 B₂₍₁₎=Ã₁Φ₁
 B₂₍₁₎=Ã₁Φ₂

julia> twistings(W,[2])
2-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterSubGroup{Perm{Int16},Int64}},1}:
 B₂₍₂₎=A₁Φ₂
 B₂₍₂₎=A₁Φ₁
```

Notice how we distinguish between subgroups generated by short roots and by
long  roots. A general  `H` corresponds to  a reductive subgroup of maximal
rank.  Here we consider the subgroup generated  by the long roots in `B₂`,
which  corresponds to a  subgroup of type  `SL₂× SL₂` in `SP₄`, and
show its possible rational forms.

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> twistings(W,[2,4])
2-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterSubGroup{Perm{Int16},Int64}},1}:
 B₂₍₂₄₎=(A₁A₁)
 B₂₍₂₄₎=A₁×A₁
```
"""
module Cosets

using ..Gapjm
export twistings, spets, Frobenius, Spets, torusfactors, subspets, 
  relative_coset, generic_sign, PhiOnDiscriminant

abstract type Spets{TW}<:Coset{TW} end

abstract type CoxeterCoset{TW}<:Spets{TW} end

function twisting_elements(W::FiniteCoxeterGroup,J::AbstractVector{<:Integer})
  if isempty(J) C=W
  elseif all(J.<=coxrank(W))
    C=Group(collect(endomorphisms(CoxGroups.parabolic_category(W,J),1)))
   else C=centralizer(W,sort(J);action=(J,w)->sort(action.(Ref(W),J,w)))
  end
  classreps(C)
end

function twisting_elements(WF::CoxeterCoset,J::AbstractVector{<:Integer})
  if isempty(J) return classreps(WF)./WF.phi end
  W=Group(WF)
  h=transporting_elt(W,sort(action.(Ref(W),J,WF.phi)),sort(J),
                             action=(x,p)->sort(action.(Ref(W),x,p)))
  if isnothing(h)
    println( "\n# no subspets for ", J )
    return Perm[]
  end
  W_L=centralizer(W,sort(collect(J)),action=(x,p)->sort(action.(Ref(W),x,p)))
  e=classreps(Group(vcat(gens(W_L),[WF.phi*h])))
  return filter(x->WF.phi*h*inv(x) in W_L,e).*inv(WF.phi)
end

Groups.Group(WF::Spets)=WF.W
PermRoot.inclusion(WF::Spets,a...)=inclusion(WF.W,a...)
PermRoot.restriction(WF::Spets,a...)=restriction(WF.W,a...)
PermRoot.semisimplerank(WF::Spets)=semisimplerank(WF.W)

"""
`twistings(W,I)`

`W`  should be a  Coxeter group.
                 
The  function returns the list, up  to `W`-conjugacy, of Coxeter sub-cosets
of  `W` whose  Coxeter group  is `reflection_subgroup(W,I)`  --- In term of
algebraic groups, it corresponds to representatives of the possible twisted
forms of the corresponding reductive subgroup of maximal rank `L`.

`W`  could also be a coset `Wϕ`; then the subgroup `L` must be conjugate to
`ϕ(L)`  for  a  rational  form  to  exist.  If `ϕ` normalizes `L`, then the
rational forms are classified by the the `ϕ`-classes of `N_W(L)/L`.

```julia-repl
julia> W=coxgroup(:E,6)
E₆

julia> WF=spets(W,Perm(1,6)*Perm(3,5))
²E₆

julia> twistings(W,2:5)
3-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterSubGroup{Perm{Int16},Int64}},1}:
 E₆₍₂₃₄₅₎=²D₄Φ₁Φ₂
 E₆₍₂₅₄₃₎=³D₄₍₁₄₃₂₎Φ₃
 E₆₍₂₃₄₅₎=D₄Φ₁²

julia> twistings(WF,2:5)
3-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterSubGroup{Perm{Int16},Int64}},1}:
 E₆₍₂₅₄₃₎=³D₄₍₁₄₃₂₎Φ₆
 E₆₍₂₃₄₅₎=²D₄Φ₁Φ₂
 E₆₍₂₃₄₅₎=D₄Φ₂²
```
"""
twistings(W,J::AbstractVector{<:Integer})=
  subspets.(Ref(W),Ref(J),twisting_elements(W,J))

"""  
`twistings(W)`

`W`  should be a Coxeter group which is not a proper reflection subgroup of
another  reflection group.  The function  returns all  'spets' representing
twisted forms of algebraic groups of type `W`.

```julia-repl
julia> twistings(coxgroup(:A,3)*coxgroup(:A,3))
8-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterGroup{Perm{Int16},Int64}},1}:
 (A₃A₃)
 ²(A₃A₃)
 ²A₃×A₃
 ²A₃×²A₃
 ²(A₃A₃)₍₁₂₃₆₅₄₎
 (A₃A₃)₍₁₂₃₆₅₄₎
 A₃×A₃
 A₃×²A₃

julia> twistings(coxgroup(:D,4))
6-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterGroup{Perm{Int16},Int64}},1}:
 ³D₄₍₁₄₃₂₎
 ²D₄₍₁₄₃₂₎
 ³D₄
 ²D₄
 ²D₄₍₂₄₃₁₎
 D₄

julia> W=rootdatum(:so,8)
D₄

julia> twistings(W)
2-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterGroup{Perm{Int16},Int64}},1}:
 ²D₄
 D₄
 
```
"""  
function twistings(W::FiniteCoxeterGroup)
  if W!=parent(W)
    error(W," must not be a proper subgroup of another reflection group")
  end
  gen=empty(gens(W))
  for (n,t) in groupby(repr,refltype(W))
    for i in 1:length(t)-1
      push!(gen,prod(Perm.(inclusion(W,t[i].indices),
                           inclusion(W,t[i+1].indices))))
    end
    J=inclusion(W,t[1].indices)
    rk=length(J)
    if t[1].series==:A 
      if rk>1  
        push!(gen,prod(i->Perm(J[i],J[rk+1-i]),1:div(rk,2)))
      end
    elseif t[1].series==:D push!(gen,Perm(J[1],J[2]))
      if rk==4  push!(gen,Perm(J[1],J[4])) end
    elseif t[1].series==:E && rk==6 
      push!(gen,Perm(J[1],J[6])*Perm(J[3],J[5]))
    end
  end
  l=filter(y->all(isinteger,matY(W.G,y)),elements(Group(gen)))
  spets.(Ref(W),l)
end

function Groups.position_class(G::Spets,g)
  if isone(G.phi) return position_class(G.W,g) end
  findfirst(c->g in c,conjugacy_classes(G))
end

#-------------- finite coxeter cosets ---------------------------------
struct FCC{T,TW<:FiniteCoxeterGroup{Perm{T}}}<:CoxeterCoset{TW}
  phi::Perm{T}
  F::Matrix
  W::TW
  prop::Dict{Symbol,Any}
end

#function Base.show(io::IO,::MIME"text/plain",t::Type{FCC{T,TW}})where {T,TW}
#  print(io,"spets{",TW,"}")
#end

spets(W::FiniteCoxeterGroup,w::Perm=Perm())=spets(W,reflrep(W,w))
Groups.Coset(W::FiniteCoxeterGroup,w::Perm=one(W))=spets(W,w)
spets(phi,F::Matrix,W::FiniteCoxeterGroup,P::Dict{Symbol,Any})=FCC(phi,F,W,P)

Base.parent(W::Spets)=gets(()->W,W,:parent)

"""
`spets(W::FiniteCoxeterGroup,F::Matrix=I(rank(W)))`

This  function returns a  Coxeter coset. `F`  must be an invertible matrix,
representing  an  automorphism  of  the  vector  space  `V` of dimension of
dimension  `rank(W)` which  induces an  automorphism of  the root system of
`parent(W)`.

The returned struct has in particular the following fields:

`.W`: the Coxeter group `W`

`.F`: the matrix acting on `V` which represents the unique element `phi` in
`WF` which preserves the positive roots.

'.phi': the permutation of the roots of `W` induced by `.F`
(also the element of smallest length in the Coset  `W .phi`).

In the first example we create a Coxeter coset corresponding to the general
unitary group `GU_3(q)` over the finite field `FF(q)`.

```julia-repl
julia> W=rootdatum(:gl,3)
A₂

julia> gu3=spets(W,-reflrep(W,W()))
²A₂Φ₂

julia> F4=coxgroup(:F,4);D4=reflection_subgroup(F4,[1,2,16,48])
F₄₍₉‚₂‚₁‚₁₆₎=D₄₍₃₂₁₄₎
```
"""
function spets(W::FiniteCoxeterGroup{Perm{T}},F::Matrix) where{T}
  perm=PermX(W.G,F)
  if isnothing(perm) error("matrix F must preserve the roots") end
  phi=reduced(W,perm)
  FCC(phi,F*reflrep(W,perm\phi),W,Dict{Symbol,Any}())
end

function torusfactors(WF::CoxeterCoset)
  M=baseX(WF.W.G)
  M*=WF.F*inv(M*E(1)//1)
  M=improve_type(M)
  r=semisimplerank(WF.W)
  CycPol(charpoly(M[r+1:end,r+1:end]))
end

"""
`torus(m::Matrix)`

`m`  should be an integral matrix of finite order. The function returns the
coset `T` of the trivial Coxeter group such that `T.F==m`. This corresponds
to  an algebraic torus `𝐓 ` of rank `size(m,1)`, with an isogeny which acts
by `m` on `X(𝐓)`.

```julia-repl
julia> torus([0 -1;1 -1])
.Φ₃
```
"""
Weyl.torus(m::Matrix)=spets(torus(size(m,1)),m)

"""
`torus(W,i)`

This  returns the torus twisted by a representative of the `i`-th conjugacy
class of `W`. This is the same as `twistings(W,Int[])[i]`.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> twistings(W,Int[])
5-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterSubGroup{Perm{Int16},Int64}},1}:
 A₃₍₎=.Φ₁³
 A₃₍₎=.Φ₁²Φ₂
 A₃₍₎=.Φ₁Φ₂²
 A₃₍₎=.Φ₁Φ₃
 A₃₍₎=.Φ₂Φ₄

julia> torus(W,2)
A₃₍₎=.Φ₁²Φ₂

julia> WF=spets(W,Perm(1,3))
²A₃

julia> twistings(WF,Int[])
5-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterSubGroup{Perm{Int16},Int64}},1}:
 A₃₍₎=.Φ₂³
 A₃₍₎=.Φ₁Φ₂²
 A₃₍₎=.Φ₁²Φ₂
 A₃₍₎=.Φ₂Φ₆
 A₃₍₎=.Φ₁Φ₄

julia> torus(WF,2)
A₃₍₎=.Φ₁Φ₂²
```
"""
Weyl.torus(W::Spets,i)=subspets(W,Int[],W.phi\classreps(W)[i])
Weyl.torus(W,i)=subspets(W,Int[],classreps(W)[i])

function Groups.conjugacy_classes(WF::Spets)
  gets(WF,:classes)do
    map(x->orbit(Group(WF),x),classreps(WF))
  end
end

function PermRoot.generic_order(WF::Spets,q)
  if rank(Group(WF))==0 return one(q) end
  generic_sign(WF)*q^sum(x->x[1]+1,codegrees(WF))*
      prod(p->1-q^p[1]*conj(p[2]),degrees(WF))
end

function PhiOnDiscriminant(WF)
  tt=refltype(WF)
  isempty(tt) ? 1 : prod(t->haskey(t,:scalar) ?  
    prod(t.scalar)^sum(degrees(t.orbit[1]).+codegrees(t.orbit[1])) : 1, tt)
end

generic_sign(WF)=(-1)^rank(Group(WF))*
   prod(last.(degrees(WF)))*conj(PhiOnDiscriminant(WF))

function PermRoot.refltype(WF::CoxeterCoset)::Vector{TypeIrred}
  gets(WF,:refltype)do
    t=refltype(WF.W)
    c=map(x->PermRoot.indices(x),t)
    phires=Perm(inclusion(WF.W.G),inclusion(WF.W.G).^WF.phi)
    map(orbits(Perm(sort.(c),map(i->sort(i.^phires),c))))do c
      o=deepcopy(t[c])
      J=PermRoot.indices(o[1])
      twist=Perm(J,J.^(phires^length(c)))
      if o[1].series==:D && length(J)==4
        if order(twist)==2
          rf=reduce(vcat,cycles(twist))
          getfield(o[1],:prop)[:indices]=J[vcat(rf,[3],setdiff([1,2,4],rf))]
        elseif order(twist)==3 && J[1]^twist!=J[2]
          getfield(o[1],:prop)[:indices]=J[[1,4,3,2]]
        end
      end
      for i in 2:length(c) 
        getfield(o[i],:prop)[:indices]=PermRoot.indices(o[i-1]).^phires
      end
      TypeIrred(Dict(:orbit=>o,:twist=>twist))
    end
  end
end

function PermRoot.parabolic_representatives(WF::CoxeterCoset,s)
  W=Group(WF)
  res=Vector{Int}[]
  for I in parabolic_representatives(W,s)
    if sort(I.^WF.phi)==sort(I) push!(res,I)
    else 
      c=filter(x->sort(x.^WF.phi)==sort(x),standard_parabolic_class(W,I))
      if !isempty(c) push!(res,c[1]) end
    end
  end
  res
end

PermRoot.Diagram(W::Spets)=PermRoot.Diagram(refltype(W))

function Base.show(io::IO, WF::Spets)
  W=Group(WF)
  if !get(io,:limit,false) && !get(io,:TeX,false)
    print(io,"spets(",W,",",WF.phi,")") 
    return
  end
  if isdefined(W,:parent)
    I=inclusiongens(W)
    n=I[PermRoot.indices(refltype(WF))]
    if n!=eachindex(gens(W.parent))
      print(io,W.parent,fromTeX(io,"_{"*joindigits(n;always=true)*"}="))
    end
  end
  PermRoot.showtypes(io,refltype(WF))
  t=torusfactors(WF)
  if !isone(t) show(io,t) end
end

PermRoot.reflrep(WF::Spets,w)=WF.F*reflrep(Group(WF),w)
  
function PermGroups.classreps(W::Spets)
  gets(W,:classreps)do
    map(x->W(x...),classinfo(W)[:classtext])
  end
end

function PermRoot.refleigen(W::Spets)
  gets(W,:refleigen)do
    map(W.phi.\classreps(W)) do x
      p=CycPol(charpoly(reflrep(W,x)))
      vcat(map(r->fill(r[1],r[2]),p.v.d)...)
    end
  end
end
  
PermRoot.refleigen(W::Spets,i)=refleigen(W)[i]

function Frobenius(WF::CoxeterCoset)
  f(w,i=1)=Frobenius(w,WF.phi^i)
end

Frobenius(w::Perm,phi)=w^phi

function twisting_elements(WF::Spets,J::AbstractVector{<:Integer})
  if isempty(J) return classreps(WF)./WF.phi end
  W=Group(WF)
  L=reflection_subgroup(W,J)
  N=normalizer(W,L)
  W_L=N/L
  if !isone(WF.phi) error( "not implemented for twisted parent Spets" ) end
  if length(W_L)>=10
    H=Group(map(x->GetRelativeAction(W,L,x.phi),gens(W_L)))
    if length(H)==length(W_L)
      h=Hom(H,W_L,gens(W_L))
      e=classreps(H)
      return map(x->h(x).phi,e)
    end
  end
  return map(x->x.phi,classreps(W_L))
end

function twisting_elements(W::PermRootGroup,J::AbstractVector{<:Integer})
  if isempty(J) classreps(W) end
  L=reflection_subgroup(W,J)
  s=unique!(sort(reflections(L)))
  C=centralizer(W,s;action=(p,g)->sort(p.^g))
  W_L=C/L
  map(x->reduced(L,x.phi),classreps(W_L))
end

function relative_coset(WF::CoxeterCoset,J)
# Print("CoxeterCosetOps.RelativeCoset ",WF,J," called \n");
  res=relative_group(Group(WF),J)
  p=isempty(J) ? WF.phi : Perm(res.prop[:parentMap],res.prop[:parentMap].^WF.phi)
  spets(res,p)
end
#-------------- subcoset ---------------------------------
"""
`subspets(WF,I,w=one(Group(WF)))`

Returns   the   reflection   subcoset   of   the   coset  `WF`  with  group
`reflection_subgroup(Group(WF),I)`  and torsion `w*WF.phi`.  `w` must be an
element  of `Group(WF)` such that  'w*WF.phi' normalizes the subroot system
generated by `I`.

```julia-repl
julia> WF=spets(coxgroup(:F,4))
F₄

julia> w=transporting_elt(Group(WF),[1,2,9,16],[1,9,16,2],action=(s,g)->s.^g);

julia> LF=subspets(WF,[1,2,9,16],w)
F₄₍₉‚₁₆‚₁‚₂₎=³D₄₍₃₄₁₂₎

julia> Diagram(LF)
ϕ acts as (1,4,2) on the component below
  O 4
  ￨
O—O—O
3 1 2
```
"""
function subspets(WF::Spets,I::AbstractVector{<:Integer},w=one(Group(WF)))
# printc("WF=",WF," I=",I," w=",w,"\n")
  RF=WF
  WF=parent(WF)
  phi=RF.phi/WF.phi
  W=Group(WF)
  if !(w in W) error(w," should be in ",W) end
  phi=w*phi
  R=reflection_subgroup(W,I)
  if (W isa CoxeterGroup) && 
     (sort(action.(Ref(R),1:2nref(R),phi*WF.phi))!=1:2nref(R))
    error("w*WF.phi does not normalize subsystem")
  end
  phi=reduced(R,phi*WF.phi)
  spets(phi,reflrep(W,phi/WF.phi)*WF.F,R,Dict{Symbol,Any}(:parent=>WF))
end

subspets(W::Group,I::AbstractVector{<:Integer},w=one(W))=subspets(spets(W),I,w)
#-------------- spets ---------------------------------
struct PRC{T,T1,TW<:PermRootGroup{T,T1}}<:Spets{TW}
  phi::Perm{T1}
  F::Matrix
  W::TW
  prop::Dict{Symbol,Any}
end

function spets(W::PermRootGroup{T,T1},w::Perm{T1}=one(W))where {T,T1}
  w=reduced(W,w)
  F=reflrep(W,w)
# println("w=$w\nF=$F")
  res=PRC(w,F,W,Dict{Symbol,Any}())
  refltype(res) # changes phi
  res
end

function spets(W::PermRootGroup,F::Matrix)
  w=PermX(W,F)
  if isnothing(w) error("matrix F must preserve the roots") end
  spets(W,w)
end

function relative_coset(WF::Spets,J=Int[],a...)
# printc("relative_coset(",WF,",",J,")\n");
  W=Group(WF)
  R=relative_group(W,J,a...)
  res=spets(R,GetRelativeAction(W,reflection_subgroup(W,J),WF.phi))
  if isempty(J)
    Group(res).prop[:MappingFromNormalizer]=R.prop[:MappingFromNormalizer]
# else Group(res).MappingFromNormalizer:=function(x)Error("MappingFromNormalizer failed");return false;end;
  end
  res
end

#function Base.show(io::IO,::MIME"text/plain",t::Type{PRC{T,T1,TW}})where {T,T1,TW}
#  print(io,"Spets{",TW,"}")
#end

Groups.Coset(W::PermRootGroup,w::Perm=one(W))=spets(W,w)

spets(phi,F::Matrix,W::PermRootGroup,P::Dict{Symbol,Any})=PRC(phi,F,W,P)

Groups.Group(W::PRC)=W.W

function PermRoot.refltype(WF::PRC)
  gets(WF,:refltype)do
    W=Group(WF)
    t=refltype(W)
    if isone(WF.phi) 
      return map(x->TypeIrred(Dict(:orbit=>[x],:twist=>Perm())),t) 
    end
    subgens=map(x->gens(reflection_subgroup(W,x.indices)),t)
    c=Perm(map(x->sort(x.^WF.phi),subgens),map(sort,subgens))
    c=orbits(inv(c))
    roots=parent(W).roots
    function scals(rr,img)
      map(ratio,roots[inclusion(W)[rr].^WF.phi],roots[inclusion(W)[img]])
    end
    map(c) do orb
      to=TypeIrred(Dict{Symbol,Any}(:orbit=>map(copy,t[orb])))
      scalar=Cyc{Rational{Int}}[]
      for i in eachindex(orb)
        if i==length(orb) next=1 else next=i+1 end
        u=Perm(subgens[orb[next]],subgens[orb[i]].^WF.phi)
        tn=t[orb[next]]
        ti=t[orb[i]]
        if i!=length(orb)  tn.indices^=u
          scal=scals(ti.indices,tn.indices)
        else to.twist=u
          scal=scals(ti.indices,tn.indices^inv(u))
        end
        if any(isnothing,scal) || !constant(scal)
          error("no element of coset acts as scalar on orbits")
          return nothing
        end
        scal=Root1(scal[1])
        sub=reflection_subgroup(W,ti.indices)
        zg=Groups.centre(sub)
        z=length(zg)
 #      println("zg=$zg")
        if z>1 # simplify scalars using centre
          e=collect(elements(zg))
          zg=e[findfirst(x->order(x)==z,e)]
          i=inclusion(sub)[1]
          v=Root1(ratio(roots[i.^zg],roots[i]))
          zg^=invmod(exponent(v),conductor(v)) # distinguished
          v=scal.*Root1.(0:z-1,z)
          m=argmin(conductor.(v))
          Perms.mul!(WF.phi,(zg^(m-1))^WF.phi)
          scal=v[m]
        end
        # simplify again by -1 in types 2A(>1), 2D(odd), 2E6
        if mod(conductor(scal),4)==2 && 
          (ti.series in [:A,:D] || (ti.series==:E && ti.rank==6))
          sb=coxgroup(ti.series,ti.rank)
          w0=sub(word(sb,longest(sb))...)
          Perms.mul!(WF.phi,w0)
          u=Perm(subgens[orb[next]],subgens[orb[i]].^WF.phi)
          #print("l=$l i=$i scal before:$scal")
          if i!=length(orb)
            tn.indices^=u
            subgens[orb[next]]=gens(reflection_subgroup(W,tn.indices))
            scal=scals(ti.indices,tn.indices)
          else to.twist=u
            scal=scals(ti.indices,tn.indices^inv(u))
          end
          #println(" scal after:$scal")
          if !constant(scal) error()
          else scal=Root1(scal[1])
          end
        end
        push!(scalar,E(scal))
      end
      to.scalar=scalar 
      to
    end
  end
end

function torusfactors(WF::PRC)
  W=Group(WF)
  M=baseX(W)
  M*=WF.F*inv(M)
  r=semisimplerank(W)
  CycPol(charpoly(M[r+1:end,r+1:end]))
end
#--------------------- Root data ---------------------------------
Weyl.rootdatum(t::String,r::Int...)=rootdatum(Symbol(t),r...)
" root datum from type "
function Weyl.rootdatum(t::Symbol,r::Int...)
   if haskey(rootdata,t) return rootdata[t](r...) end
   error("Unknown root datum $(repr(t)). Known types are:\n",
              join(sort(collect(keys(rootdata))),", "))
end

id(r)=one(fill(0,r,r))

const  rootdata=Dict{Symbol,Function}()
rootdata[:gl]=function(r)
  if r==1 return torus(1) end
  R=id(r)
  R=R[1:end-1,:]-R[2:end,:]
  rootdatum(R,R)
end
rootdata[:sl]=r->rootdatum(cartan(:A,r-1),id(r-1))
rootdata[:tgl]=function(n, k)
  if gcd(n,k)!=1 error(k," should be prime to ",n) end
  X=hcat(cartan(:A,n-1),fill(0,n-1))
  # We intertwine the last weight with the torus
  lat=vcat(X,permutedims(vcat(fill(0,n-2), [k,1])))
  # Get a basis for the sublattice
  lat=MatInt.HermiteNormalFormIntegerMat(toL(lat))
  # Find the roots in a basis of the sublattice
  Y=id(n)[1:n-1,:]
  return rootdatum(toM(map(x->MatInt.SolutionIntMat(lat,x),toL(X))),
                   Y*permutedims(toM(lat)))
end
rootdata[:pgl]=r->coxgroup(:A,r-1)
rootdata[:sp]=function(r)
  R=id(div(r,2))
  for i in 2:div(r,2) R[i,i-1]=-1 end
  R1=copy(R)
  R1[1,:].*=2
  rootdatum(R1,R)
end
rootdata[:csp]=function(r)
  r=div(r,2)
  R=id(r+1)
  R=R[1:r,:]-R[2:end,:] 
  cR=copy(R) 
  R[1,1:2]=[0,-1] 
  cR[1,1:2]=[1,-2]
  rootdatum(cR,R)
end
rootdata[:psp]=r->coxgroup(:C,div(r,2))
rootdata[:so]=function(r)
  R=id(div(r,2)) 
  for i in 2:div(r,2) R[i,i-1]=-1 end
  if isodd(r) R1=copy(R)
    R1[1,1]=2
    rootdatum(R,R1)
  else R[1,2]=1 
    rootdatum(R,R)
  end
end
rootdata[:pso]=function(r)
  r1=div(r,2)
  isodd(r) ? coxgroup(:B,r1) : coxgroup(:D,r1)
end
rootdata[:spin]=function(r)
  r1=div(r,2)
  isodd(r) ? rootdatum(cartan(:C,r1),id(r1)) : rootdatum(cartan(:D,r1),id(r1))
end
rootdata[:halfspin]=function(r)
  if !iszero(r%4) error("halfspin groups only exist in dimension 4r") end
  r=div(r,2)
  R=id(r)
  R[r,:]=vcat([-div(r,2),1-div(r,2)],2-r:1:-2,[2])
  rootdatum(R,Int.(cartan(:D,r)*permutedims(inv(Rational.(R)))))
end
rootdata[:gpin]=function(r)
  d=div(r,2)
  R=id(d+1)
  R=R[1:d,:]-R[2:end,:]
  cR=copy(R)
  if isodd(r)
    R[1,1]=0
    cR[1,2]=-2
  else
    R[1,3]=-1
    cR[1,1:3]=[0,-1,-1]
    R=hcat(fill(0,d),R)
    cR=hcat(fill(0,d),cR)
    cR[1,1]=1
  end
  rootdatum(R,cR)
end
rootdata[:E6]=()->coxgroup(:E,6)
rootdata[:E7]=()->coxgroup(:E,7)
rootdata[:E8]=()->coxgroup(:E,8)
rootdata[:E6sc]=()->rootdatum(cartan(:E,6),id(6))
rootdata[:E7sc]=()->rootdatum(cartan(:E,7),id(7))
rootdata[:F4]=()->coxgroup(:F,4)
rootdata[:G2]=()->coxgroup(:G,2)
rootdata[:u]=r->spets(rootdatum(:gl,r),reverse(-id(r),dims=1))
rootdata[:su]=r->r==2 ? spets(rootdatum(:sl,r)) :
    spets(rootdatum(:sl,r),prod(i->Perm(i,r-i),1:div(r-1,2)))
rootdata[:psu]=function(r)g=rootdatum(:pgl,r)
  r==2 ? g : spets(g,prod(i->Perm(i,r-i),1:div(r-1,2))) end
rootdata[Symbol("3gpin8")]=()->spets(rootdatum(:gpin,8),[1 1 1 0 0 0;
 -2 0 -1 -1 -1 -1;-1 0 -1 0 0 -1;-1 0 -1 0 -1 0;
 -1 0 -1 -1 0 0;-1 -1 0 0 0 0])
rootdata[Symbol("gpin-")]=function(r)
  d=r/2
  F=id(d+2)
  F[1,1:3]=[1,-1,1]
  F[1:d+2,2]=fill(-1,d+2)
  F[2:3,2:3]=-id(2)
  spets(rootdatum(:gpin,r),F)
end
rootdata[Symbol("so-")]=r->spets(rootdatum(:so,r),Perm(1,2))
rootdata[Symbol("pso-")]=r->spets(rootdatum(:pso,r),Perm(1,2))
rootdata[Symbol("spin-")]=r->spets(rootdatum(:spin,r),Perm(1,2))
rootdata[Symbol("2I")]=e->spets(coxgroup(:Isym,2,e),Perm(1,2))
rootdata[:suzuki]=()->spets(coxgroup(:Bsym,2),Perm(1,2))
rootdata[Symbol("2B2")]=()->rootdatum(:suzuki)
rootdata[:ree]=()->spets(coxgroup(:Gsym,2),Perm(1,2))
rootdata[Symbol("2G2")]=()->rootdatum(:ree)
rootdata[:triality]=()->spets(coxgroup(:D,4),Perm(1,2,4))
rootdata[Symbol("3D4")]=()->rootdatum(:triality)
rootdata[Symbol("3D4sc")]=()->spets(rootdatum(:spin,8),Perm(1,2,4))
rootdata[Symbol("2E6")]=()->spets(coxgroup(:E,6),Perm(1,6)*Perm(3,5))
rootdata[Symbol("2E6sc")]=()->spets(rootdatum(:E6sc),Perm(1,6)*Perm(3,5))
rootdata[Symbol("2F4")]=()->spets(coxgroup(:Fsym,4),Perm(1,4)*Perm(2,3))
end
