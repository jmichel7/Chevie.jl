"""
Let  `R` be a  root system in  the real vector  space `V`. We say that `F₀∈
GL(V)`  is an  *automorphism of  `R`* if  it permutes  `R` and is of finite
order  (finite  order  is  automatic  if  `R` generates `V`). It follows by
[chap.  VI,  §1.1,  lemme  1  Bourbaki1968](biblio.htm#Bou68) that the dual
`F₀*∈  GL(V*)`  permutes  the  coroots  `R*⊂  V*`; thus `F₀` normalizes the
reflection  group  `W`  associated  to  `R`,  that  is  `w↦  F₀wF₀⁻¹` is an
automorphism  of `W`. Thus we get a reflection coset `WF₀`, which we call a
*Coxeter coset*.

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
characters and almost characters, etc… (see, e.g.,
[Broue-Malle-Michel1993](biblio.htm#BMM93)).  It is thus possible to extend
their  construction to non-crystallographic groups (or even to more general
complex  reflection groups,  see [`spets`](@ref));  this is  why we did not
include  a root  system in  the definition  of a reflection coset. However,
unipotent conjugacy classes for instance depend on the root system.

We assume now that `𝐓` is contained in an `F`-stable Borel subgroup of `𝐆`.
This  defines an order  on the roots,  and there is  a unique element `ϕ∈ W
F₀`,  the  *reduced  element*  of  the  coset,  which  preserves the set of
positive  roots.  It  thus  defines  a  *diagram  automorphism*, that is an
automorphism  of the Coxeter system `(W,S)`.  This element is stored in the
component  `.phi` of the coset record. It may be defined without mentioning
the  roots,  as  follows:  `(W,F₀(S))`  is  another  Coxeter  system,  thus
conjugate to `S` by a unique element of `W`, thus there is a unique element
`ϕ∈  WF₀` which stabilizes `S` (a proof  follows from [Theoreme 1, chap. V,
§3  Bourbaki1968](biblio.htm#Bou68)). We  consider thus  cosets of the form
`Wϕ` where `ϕ` stabilizes `S`. The coset `W ϕ` is completely defined by the
permutation  `.phi`  when  `𝐆`  is  semi-simple  ---  equivalently when `Φ`
generates  `V`; in this  case we just  need to specify  `phi` to define the
coset.

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

## Spets

We  now extend the above notions  to general complex reflection groups. Let
`W⊂  GL(V)` be a complex reflection group  on the vector space `V`. Let `ϕ`
be  an element  of `GL(V)`  which normalizes  `W`. Then  the coset  `Wϕ` is
called a reflection coset.

A reference for these cosets is [Broue-Malle-Michel 1999](biblio.htm#BMM99).
When `W` is a so-called *Spetsial* group, they are the basic object for the
construction  of  a  *Spetses*,  which  is  an object attached to a complex
reflection  group from which one can derive combinatorially some attributes
shared with finite reductive groups, like unipotent degrees, etc….

We  say that  a reflection  coset is  irreducible if  `W` is irreducible. A
general  coset is a direct  product of *descents of  scalars*, which is the
case  where `ϕ`  is transitive  on the  irreducible components  of `W`. The
irreducible    cosets   have   been   classified   in   [Broue-Malle-Michel
1999](biblio.htm#BMM99):  up to multiplication of `ϕ` by a scalar, there is
usually only one or two possible cosets for a given irreducible group.

We  deal only  with *finite  order* cosets,  that is,  we assume there is a
(minimal) integer `δ` such that `(Wϕ)^δ=W`. Then the group generated by `W`
and `ϕ` is finite, of order `δ|W|`.

A  subset `C`  of a  `Wϕ` is  called a  *conjugacy class*  if one of the
following equivalent conditions is fulfilled:

  * `C` is the orbit of an element in `Wϕ` under the conjugation action of `W`.

  * `C` is a conjugacy class of `⟨W,ϕ⟩` contained in `Wϕ`.

  * The set `{w∈ W|wϕ∈ C}` is  a `ϕ`-conjugacy  class of `W` (two elements
`v,w∈  W` are called `ϕ`-conjugate, if and only if there exists `x∈ W` with
`v=xwϕ(x⁻¹)`).

An irreducible character of `⟨W,ϕ⟩` has some non-zero values on `Wϕ` if and
only if its restriction to `W` is irreducible. Further, two characters `χ₁`
and  `χ₂`  which  have  same  irreducible  restriction  to  `W` differ by a
character  of  the  cyclic  group  `⟨ϕ⟩`  (which identifies to the quotient
`⟨W,ϕ⟩/W`). A set containing one extension to `⟨W,ϕ⟩` of each `ϕ`-invariant
character  of `W` is called a *set  of irreducible characters of `Wϕ`*. Two
such  characters  are  orthogonal  for  the  scalar  product  on  the class
functions on `Wϕ` given by ``⟨χ,ψ⟩:=|W|¹∑_{w∈ W}χ(wϕ)\\overline{ψ(wϕ)}.``
For rational groups (Weyl groups), Lusztig has defined a choice of a set of
irreducible  characters for  `Wϕ` (called  the *preferred extensions*), but
for  more  general  reflection  cosets  we  have made some rather arbitrary
choices,  which  however  have  the  property  that their values lie in the
smallest possible field.

The  *character  table*  of  `Wϕ`  is  the  table  of  values  of  a set of
irreducible characters on the conjugacy classes.

A *subcoset* `Lwϕ` of `Wϕ` is given by a reflection subgroup `L` of `W` and
an element `w` of `W` such that `wϕ` normalizes `L`.

We  then have a natural notion of  *restriction* of class functions on `Wϕ`
to  class  functions  on  `Lwϕ`  as  well  as  of  *induction* in the other
direction.  These  maps  are  adjoint  with  respect  to the scalar product
defined above (see [Broue-Malle-Michel 1999](biblio.htm#BMM99)).

In  this package the most general construction  of a reflection coset is by
starting  from a reflection datum, and giving in addition the matrix `F` of
the  map `ϕ:V→ V`  (see the command  `spets`). However, at present, general
cosets are only implemented for groups represented as permutation groups on
a  set of roots, and  it is required that  the automorphism given preserves
this  set up to  a scalar (it  is allowed that  these scalars depend on the
pair  of an  irreducible component  and its  image). It  is also allowed to
specify  `ϕ` by the permutation it induces on the roots; in this case it is
assumed  that `ϕ` acts  trivially on the  orthogonal of the  roots, but the
roots  could be those of a parent group, generating a larger space. Thus in
any  case we have  a permutation representation  of `⟨W,ϕ⟩` and we consider
the coset to be a set of permutations.

Reflection  cosets  are  implemented  in  by  a  `struct` which points to a
reflection  group  record  and  has  additional  fields holding `F` and the
corresponding  permutation `phi`. In the general case, on each component of
`W`  which is  a descent  of scalars,  `F` will  permute the components and
differ  by a scalar on each  component from an automorphism which preserves
the  roots. In this case, we have  a permutation `phi` and a `scalar` which
is stored for that component.

The  most common situation where cosets  with non-trivial `phi` arise is as
sub-cosets  of reflection groups. Here is an "exotic" example, see the next
chapter for more classical examples involving Coxeter groups.

```julia-repl
julia> W=complex_reflection_group(14)
G₁₄

julia> R=reflection_subgroup(W,[2,4])
G₁₄₍₂₄₎=G(-ζ₃²√-2)₅

julia> RF=spets(R,W(1))
G₁₄₍₂₄₎=²G(-ζ₃²√-2)₅
```

```julia-rep1
julia> diagram(RF)
ϕ acts as (1,2) on the component below
③ ══③ G₅
1   2
```

```julia-repl
julia> degrees(RF)
2-element Vector{Tuple{Int64, Cyc{Int64}}}:
 (6, 1)
 (12, -1)
```

The  last line shows for each  reflection degree the corresponding *factor*
of  the coset, which is  the scalar by which  `ϕ` acts on the corresponding
fundamental reflection invariant. The factors characterize the coset.

A  spets by default is  printed in an abbreviated  form which describes its
type,  as above (`G₅` twisted by 2, with a Cartan matrix which differs from
the  standard one by  a factor of  `√6`). The function  `repr` gives a form
which could be input back in Julia. With the same data as above we have:

```julia-rep1
julia> print(RF)
spets(reflection_subgroup(complex_reflection_group(14),[2, 4]),perm"(1,3)(2,4)(5,9)(6,10)(7,11)(8,12)(13,21)(14,22)(15,23)(16,24)(17,25)(18,26)(19,27)(20,28)(29,41)(30,42)(31,43)(32,44)(33,45)(34,46)(35,47)(36,48)(37,49)(38,50)(39,51)(40,52)(53,71)(54,72)(55,73)(56,74)(57,75)(58,76)(59,77)(60,78)(62,79)(64,80)(65,81)(66,82)(67,69)(68,70)(83,100)(84,101)(85,102)(87,103)(89,99)(90,97)(91,98)(92,96)(93,104)(94,95)(105,113)(106,114)(109,111)(110,112)(115,118)(116,117)(119,120)")
```

Conjugacy  classes and irreducible characters of Coxeter cosets are defined
as  for  general  reflection  cosets.  For  irreducible  characters of Weyl
cosets,  we choose (following Lusztig) for each `ϕ`-stable character of `W`
a  particular extension to a character of  `W⋊ ⟨ϕ⟩`, which we will call the
*preferred extension*. The *character table* of the coset `Wϕ` is the table
of  the restrictions to  `Wϕ` of the  preferred extensions. The question of
finding the conjugacy classes and character table of a Coxeter coset can be
reduced to the case of irreducible root systems `R`.

  * The automorphism `ϕ` permutes the irreducible components of `W`, and
    `Wϕ`  is a direct  product of cosets  where `ϕ` permutes cyclically the
    irreducible components of `W`. The preferred extension is defined to be
    the  direct  product  of  the  preferred  extension  in  each  of these
    situations.

  * Assume now that `Wϕ` is a descent of scalars, that is the decomposition
    in irreducible components `W=W₁× ⋯ × Wₖ` is cyclically permuted by `ϕ`.
    Then there are natural bijections from the `ϕ`-conjugacy classes of `W`
    to  the `ϕᵏ`-conjugacy classes  of `W₁` as  well as from the `ϕ`-stable
    characters  of `W` to the `ϕᵏ`-stable  characters of `W₁`, which reduce
    the  definition of preferred  extensions on `Wϕ`  to the definition for
    `W₁ϕᵏ`.

  * Assume now  that `W`  is the  Coxeter group  of an  irreducible root
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

  * In case  `²Dₙ` the  group `W⋊ ⟨ϕ⟩`  is isomorphic  to a Coxeter
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
Coxeter   datum   specified   by   the   matrices   of   `simpleRoots`  and
`simplecoroots`,  and  giving  in  addition  the  matrix `F0Mat` of the map
`F₀:V→ V` (see the commands  `CoxeterCoset` and `CoxeterSubCoset`). As for
Coxeter  groups,  the  elements  of  `Wϕ`  are  uniquely  determined by the
permutation  they  induce  on  the  set  of  roots  `R`.  We consider these
permutations as `elements` of the Coxeter coset.

Coxeter  cosets are implemented by a struct which points to a Coxeter datum
record  and  has  additional  fields  holding `F0Mat` and the corresponding
element  `phi`. Functions on the coset (for example, `classinfo`) are about
properties  of  the  group  coset  `W  ϕ`  ;  however, most definitions for
elements of untwisted Coxeter groups apply without change to elements in `W
ϕ`:  e.g., if we define the length of  an element `wϕ∈ Wϕ` as the number of
positive  roots it sends to negative ones, it  is the same as the length of
`w`,  i.e., `ϕ` is of length `0`, since `ϕ` has been chosen to preserve the
set of positive roots. Similarly, the Coxeter `word` describing `wϕ` is the
same as the one for `w`, etc…

We associate to a Coxeter coset `Wϕ` a *twisted Dynkin diagram*, consisting
of  the Dynkin diagram of `W` and  the graph automorphism induced by `ϕ` on
this  diagram (this  specifies the  group `W⋊  ⟨F⟩`, mentioned above, up to
isomorphism).  See  the  functions  `refltype`,  `print`  and `diagram` for
Coxeter cosets.

Below  is an example showing first how to *not* define, then how to define,
the Weyl coset for a Suzuki group:

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> spets(W,Perm(1,2))
ERROR: matrix F must preserve the roots
Stacktrace:
 [1] error(::String) at ./error.jl:33
 [2] spets(::Chevie.Weyl.FCG{Int16,Int64,PRG{Int64,Int16}}, ::Matrix{Int64}) at /home/jmichel/julia/Chevie/src/Cosets.jl:241 (repeats 2 times)
 [3] top-level scope at REPL[19]:1

julia> W=coxgroup(:Bsym,2)
Bsym₂

julia> WF=spets(W,Perm(1,2))
²Bsym₂

julia> CharTable(WF)
CharTable(²Bsym₂)
┌───┬─────────┐
│   │    1 121│
├───┼─────────┤
│2. │1   1   1│
│.11│1  -1  -1│
│1.1│. -√2  √2│
└───┴─────────┘
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
²Bsym₂₍₎=Φ‴₈
```

A subgroup `H` which is a parabolic subgroup corresponds to a rational form
of  a Levi  subgroup of  `𝐆`. The  command `twistings`  gives all rational
forms of such a Levi.

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> twistings(W,[1])
2-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 B₂₍₁₎=Ã₁Φ₁
 B₂₍₁₎=Ã₁Φ₂

julia> twistings(W,[2])
2-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 B₂₍₂₎=A₁Φ₁
 B₂₍₂₎=A₁Φ₂
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
2-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 B₂₍₂₄₎=A₁×A₁
 B₂₍₂₄₎=(A₁A₁)
```
"""
module Cosets

using ..Chevie
export twistings, spets, Frobenius, Spets, subspets,
  relative_coset, generic_sign, PhiOnDiscriminant, graph_automorphisms,
  CoxeterCoset, twisted_power

Base.rand(G::Coset)=rand(Group(G))*G.phi # put into Groups

abstract type Spets{T,TW}<:NormalCoset{T,TW} end

abstract type CoxeterCoset{T,TW}<:Spets{T,TW} end

extprod(W1::Spets,W2::Spets)=spets(W1.W*W2.W,cat(W1.F,W2.F,dims=(1,2)))
extprod(W1::Spets,W2::FiniteCoxeterGroup)=extprod(W1,spets(W2))
extprod(W1::FiniteCoxeterGroup,W2::Spets)=extprod(spets(W1),W2)
extprod(W1::FiniteCoxeterGroup,W2::FiniteCoxeterGroup)=W1*W2

function graph_automorphisms(W::FiniteCoxeterGroup,J::AbstractVector{<:Integer})
  if isempty(J) return W end
  H=reflection_subgroup(W,J)
  if !isparabolic(W,H)
    return stabilizer(W,sort(inclusiongens(H)),onsets)
  end
  if issubset(inclusiongens(H),1:ngens(W)) p=one(W)
  else p=standard_parabolic(W,H); H=H^p
  end
  J=inclusiongens(H)
  C=collect(endomorphisms(CoxGroups.parabolic_category(W,J),1))
  Group(C.^inv(p))
end

twisting_elements(W::FiniteCoxeterGroup,J::AbstractVector{<:Integer})=
  classreps(graph_automorphisms(W,J))

function twisting_elements(WF::CoxeterCoset,J::AbstractVector{<:Integer})
  if isempty(J) return classreps(WF)./WF.phi end
  if isone(WF.phi) return twisting_elements(Group(WF),J) end
  W=Group(WF)
  iJ=inclusion(W,J)
  h=transporting_elt(W,onsets(iJ,WF.phi),sort(iJ),onsets)
  if isnothing(h)
    println("\n# no subspets for ",J)
    return eltype(W)[]
  end
  W_L=graph_automorphisms(W,iJ)
  classreps(NormalCoset(W_L,WF.phi*h))./WF.phi
end

Groups.Group(WF::Spets)=WF.W
@forward Spets.W PermRoot.inclusion, PermRoot.restriction, 
  PermRoot.inclusiongens, PermRoot.semisimplerank, Symbols.rank

Weyl.dimension(WF::CoxeterCoset)=dimension(WF.W)

CoxGroups.word(WF::CoxeterCoset,w)=word(WF.W,w/WF.phi)

"""
`twistings(W,I)`

`W`  should be a  complex reflection group.

The  function returns  the list,  up to  `W`-conjugacy, of  subspets of `W`
whose  group is `reflection_subgroup(W,I)` --- In  the case of Weyl groups,
it  corresponds to  representatives of  the possible  twisted forms  of the
reductive subgroup of maximal rank `L` corresponding to
`reflection_subgroup(W,I)`.   These  are  classified   by  the  classes  of
`N_W(L)/L`.

`W`  could also be a coset `Wϕ`; then the subgroup `L` must be conjugate to
`ϕ(L)`  for a  rational form  to exist.  If `wϕ`  normalizes `L`,  then the
rational forms are classified by the the `wϕ`-classes of `N_W(L)/L`.

```julia-repl
julia> W=coxgroup(:E,6)
E₆

julia> WF=spets(W,Perm(1,6)*Perm(3,5))
²E₆

julia> twistings(W,2:5)
3-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 E₆₍₂₃₄₅₎=D₄Φ₁²
 E₆₍₂₃₄₅₎=³D₄Φ₃
 E₆₍₂₃₄₅₎=²D₄Φ₁Φ₂


julia> twistings(WF,2:5)
3-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 ²E₆₍₂₃₄₅₎=²D₄₍₁₄₃₂₎Φ₁Φ₂
 ²E₆₍₂₃₄₅₎=³D₄₍₁₄₃₂₎Φ₆
 ²E₆₍₂₃₄₅₎=D₄Φ₂²
```
"""
function twistings(W,J::AbstractVector{<:Integer})
  subspets.(Ref(W),Ref(J),twisting_elements(W,J))
end

twistings(W,R::FiniteCoxeterGroup)=twistings(W,inclusiongens(R,W))

"""
`graph_automorphisms(t::Vector{TypeIrred})`

Given  the `refltype` of a  finite Coxeter group, returns  the group of all
Graph automorphisms of `t` as a group of permutations of `indices(t)`.

```julia-repl
julia> W=coxgroup(:D,4)
D₄

julia> graph_automorphisms(refltype(W*W))
Group((1,5)(2,6)(3,7)(4,8),(1,2),(1,4))
```
"""
function graph_automorphisms(t::Vector{TypeIrred})
  gen=empty([Perm()])
  for (n,t) in groupby(repr,t)
    for i in 1:length(t)-1
      push!(gen,prod(Perm.(t[i].indices,t[i+1].indices)))
    end
    J=t[1].indices
    rk=length(J)
    if t[1].series==:A
    if rk>1
    push!(gen,prod(i->Perm(J[i],J[rk+1-i]),1:div(rk,2)))
    end
    elseif t[1].series==:B && t[1].cartanType==root(2) push!(gen,Perm(J[1],J[2]))
    elseif t[1].series==:D push!(gen,Perm(J[1],J[2]))
    if rk==4  push!(gen,Perm(J[1],J[4])) end
    elseif t[1].series==:E && rk==6
    push!(gen,Perm(J[1],J[6])*Perm(J[3],J[5]))
    elseif t[1].series==:F && t[1].cartanType==root(2)
    push!(gen,Perm(J[1],J[4])*Perm(J[2],J[3]))
    elseif t[1].series==:G && t[1].cartanType==root(3) push!(gen,Perm(J[1],J[2]))
    end
  end
  Group(gen)
end

graph_automorphisms(W::FiniteCoxeterGroup)=
  Group(map(g->PermX(W,reflrep(W,g)),gens(graph_automorphisms(refltype(W)))))

"""
`twistings(W)`

`W`  should be a Coxeter group which is not a proper reflection subgroup of
another reflection group (so that `inclusion(W)==eachindex(roots(W))`). The
function returns all `spets` representing twisted forms of algebraic groups
of type `W`.

```julia-repl
julia> twistings(coxgroup(:A,3)*coxgroup(:A,3))
8-element Vector{Spets{FiniteCoxeterGroup{Perm{Int16},Int64}}}:
 A₃×A₃
 A₃×²A₃
 ²A₃×A₃
 ²A₃×²A₃
 (A₃A₃)
 ²(A₃A₃)
 ²(A₃A₃)₍₁₂₃₆₅₄₎
 (A₃A₃)₍₁₂₃₆₅₄₎

julia> twistings(coxgroup(:D,4))
6-element Vector{Spets{FiniteCoxeterGroup{Perm{Int16},Int64}}}:
 D₄
 ²D₄₍₂₄₃₁₎
 ²D₄
 ³D₄
 ²D₄₍₁₄₃₂₎
 ³D₄₍₁₄₃₂₎

julia> W=rootdatum(:so,8)
so₈

julia> twistings(W)
2-element Vector{Spets{FiniteCoxeterGroup{Perm{Int16},Int64}}}:
 D₄
 ²D₄

```
"""
function twistings(W::FiniteCoxeterGroup)
  if W!=parent(W)
    error(W," must not be a proper subgroup of another reflection group")
  end
  l=elements(graph_automorphisms(refltype(W)))
  l=filter(y->all(isinteger,YMatrix(W.G,y)),l)
  spets.(Ref(W),l)
end

#function Groups.position_class(G::Coset,g)
#  if isone(G.phi) return position_class(G.W,g) end
#  findfirst(c->g in c,conjugacy_classes(G))
#end

twisted_power(x,n,F::Function)=iszero(n) ? one(x) : x*F(twisted_power(x,n-1,F))
#-------------- finite coxeter cosets ---------------------------------
@GapObj struct FCC{T,TW<:FiniteCoxeterGroup{T}}<:CoxeterCoset{T,TW}
  phi::T
  F::Matrix
  W::TW
end

function Base.show(io::IO,t::Type{FCC{T,TW}})where {T,TW}
#function Base.show(io::IO,::MIME"text/plain",t::Type{FCC{T,TW}})where {T,TW}
  if get(io,:limit,true) print(io,"Spets{",TW,"}") 
  else print(io,t) 
  end
end

Groups.NormalCoset(W::FiniteCoxeterGroup,w::Perm=one(W))=spets(W,w)
spets(phi,F::Matrix,W::FiniteCoxeterGroup,P::Dict{Symbol,Any})=FCC(phi,F,W,P)

Base.parent(W::Spets)=get!(()->W,W,:parent)

"""
`spets(W::ComplexReflectionGroup, F::Matrix=I(rank(W)))`

This  function returns a or  a `CoxeterCoset` or a  `Spets`. `F` must be an
invertible  matrix, representing an automorphism of the vector space `V` of
dimension  of dimension `rank(W)` which for  a finite Coxeter group induces
an  automorphism of the root  system of `parent(W)`, or  for a more general
complex reflection group just stabilizes `W`.

The returned struct has in particular the following fields:

`.W`: the group `W`

`.F`: the matrix acting on `V` which represents the unique element `phi` in
`WF` which preserves the positive roots (for finite Coxeter groups) or some
"canonical" representative of the coset for more general complex reflection
groups.

`.phi`:  a `Perm`, the permutation of the roots of `W` induced by `.F` (for
general  complex reflection groups this may be a permutation up to scalars)
(also  for Coxeter groups the element of smallest length in the NormalCoset
`W .phi`).

In the first example we create a Coxeter coset corresponding to the general
unitary group `GU_3(q)` over the finite field `FF(q)`.

```julia-repl
julia> W=rootdatum(:gl,3)
gl₃

julia> gu3=spets(W,-reflrep(W,W()))
²A₂Φ₂

julia> F4=coxgroup(:F,4);D4=reflection_subgroup(F4,[1,2,16,48])
F₄₍₁‚₂‚₉‚₁₆₎=D₄₍₃₂₁₄₎

julia> spets(D4,[1 0 0 0;0 1 2 0;0 0 0 1;0 0 -1 -1])
F₄₍₉‚₁₆‚₁‚₂₎=³D₄₍₃₄₁₂₎
```
`spets(W::ComplexReflectionGroup,p::Perm)`

In  this version `F` is  defined by the permutation  of the simple roots it
does.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> spets(W,Perm(1,3))
²A₃
```
"""
function spets(W::FiniteCoxeterGroup{Perm{T}},F::AbstractMatrix) where{T}
  perm=PermX(W.G,F)
  if isnothing(perm) error("matrix F must preserve the roots") end
  phi=reduced(W,perm)
  FCC(phi,F*reflrep(W,perm\phi),W,Dict{Symbol,Any}())
end

function spets(W::FiniteCoxeterGroup{Perm{T}}) where{T}
  get!(W,:trivialspets)do
    FCC(Perm{T}(),reflrep(W,one(W)),W,Dict{Symbol,Any}())
  end
end

spets(W::FiniteCoxeterGroup,w::Perm)=spets(W,reflrep(W,w))

PermRoot.radical(WF::CoxeterCoset)=torus(central_action(Group(WF),WF.F))

"""
`torus(m::AbstractMatrix)`

`m`  should be a matrix of finite order. The function returns the coset `T`
of  the  trivial  group  such  that  `T.F==m`.  When  `m`  is  integral his
corresponds to an algebraic torus `𝐓 ` of rank `size(m,1)`, with an isogeny
which acts by `m` on `X(𝐓)`.

```julia-repl
julia> torus([0 -1;1 -1])
Φ₃
```
"""
Weyl.torus(m::AbstractMatrix)=spets(torus(size(m,1)),m)

"""
`torus(W,i)`

where  `W` is  a `Spets`  or a  `ComplexReflectionGroup`. This  returns the
torus  twisted by  a representative  of the  `i`-th conjugacy class of `W`.
This is the same as `twistings(W,Int[])[i]`.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> twistings(W,Int[])
5-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 A₃₍₎=Φ₁³
 A₃₍₎=Φ₁²Φ₂
 A₃₍₎=Φ₁Φ₂²
 A₃₍₎=Φ₁Φ₃
 A₃₍₎=Φ₂Φ₄

julia> torus(W,2)
A₃₍₎=Φ₁²Φ₂

julia> WF=spets(W,Perm(1,3))
²A₃

julia> twistings(WF,Int[])
5-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 ²A₃₍₎=Φ₂³
 ²A₃₍₎=Φ₁Φ₂²
 ²A₃₍₎=Φ₁²Φ₂
 ²A₃₍₎=Φ₂Φ₆
 ²A₃₍₎=Φ₁Φ₄

julia> torus(WF,2)
²A₃₍₎=Φ₁Φ₂²
```
"""
Weyl.torus(W::Spets,i::Integer)=subspets(W,Int[],W.phi\classreps(W)[i])
Weyl.torus(W::ComplexReflectionGroup,i::Integer)=subspets(W,Int[],classreps(W)[i])

function Groups.elements(C::ConjugacyClass{T,TW})where {T,TW<:Spets}
  get!(C,:elements)do
    orbit(Group(C.G),C.representative)
  end
end

function PermRoot.generic_order(WF::Spets,q;nc=false)
  if rank(Group(WF))==0 return one(q) end
  if nc 
    prod(x->conj(x[2]),degrees(WF))^2*q^sum(x->x[1]-1,degrees(WF))*
    prod(((d,e),)->q^d-e,degrees(WF))
  else generic_sign(WF)*q^sum(x->x[1]+1,codegrees(WF))*
    prod(((d,e),)->1-q^d*conj(e),degrees(WF))
  end
end

function PhiOnDiscriminant(WF)
  tt=refltype(WF)
  isempty(tt) ? 1 : prod(t->haskey(t,:scalar) ?
    prod(t.scalar)^sum(degrees(t.orbit[1]).+codegrees(t.orbit[1])) : 1, tt)
end

generic_sign(WF)=(-1)^rank(Group(WF))*
   prod(last.(degrees(WF)))*conj(PhiOnDiscriminant(WF))

function PermRoot.refltype(WF::CoxeterCoset)::Vector{TypeIrred}
  get!(WF,:refltype)do
    W=Group(WF)
    t=refltype(W)
    c=indices.(t)
    phires=Perm(action.(Ref(W),eachindex(roots(W)),WF.phi))
    map(orbits(Perm(sort.(c),map(i->sort(i.^phires),c)),eachindex(c)))do c
      o=copy.(t[c])
      J=indices(o[1])
      twist=Perm(J,J.^(phires^length(c)))
      if o[1].series==:D && length(J)==4
        if order(twist)==2
          rf=reduce(vcat,cycles(twist))
          o[1].indices=J[vcat(rf,[3],setdiff([1,2,4],rf))]
        elseif twist==Perm(1,4,2)
          o[1].indices=J[[1,4,3,2]]
          twist=Perm(1,2,4)
        end
      end
      for i in 2:length(c) o[i].indices=indices(o[i-1]).^phires end
      TypeIrred(;orbit=o,twist)
    end
  end
end

function PermRoot.parabolic_reps(WF::CoxeterCoset,s)
  W=Group(WF)
  res=Vector{Int}[]
  for I in parabolic_reps(W,s)
   if sort(action.(Ref(W),I,WF.phi))==sort(I) push!(res,I)
    else
      c=filter(x->sort(action.(Ref(W),x,WF.phi))==sort(x),
                 standard_parabolic_class(W,I))
      if !isempty(c) push!(res,c[1]) end
    end
  end
  res
end

PermRoot.parabolic_reps(WF::Spets,s)=if !isone(WF.phi) error() else
  parabolic_reps(Group(WF),s) end

function Base.show(io::IO, WF::Spets)
  if haskey(WF,:callname) 
    if hasdecor(io) printTeX(io,WF.TeXcallname) 
    else print(io,WF.callname)
    end
    return 
  end
  W=Group(WF)
  if !hasdecor(io)
    print(io,"spets(",W,",",WF.phi,")")
    return
  end
  if get(io,:parent,true)
    if haskey(WF,:parent)
#     n=inclusion(WF,indices(refltype(WF)))
      n=inclusiongens(WF)
      if n!=eachindex(gens(Group(WF.parent)))
        printTeX(io,"{",WF.parent,"}_{"*joindigits(n;always=true)*"}=")
      end
    elseif isdefined(W,:parent)
      n=inclusion(W,indices(refltype(WF)))
      if n!=eachindex(gens(W.parent))
        printTeX(io,"{",W.parent,"}_{"*joindigits(n;always=true)*"}=")
      end
    end
  end
  show(io,refltype(WF))
  t=CycPol(1,0,ModuleElt(map(x->x.r,torusfactors(WF)).=>1))
  if !isone(t) show(io,t) 
  elseif isempty(refltype(WF)) print(io,".")
  end
end

PermRoot.reflrep(WF::Spets,w)=WF.F*reflrep(Group(WF),WF.phi\w)

"""
    Frobenius(WF)(x,i=1)

If   `WF`  is  a  Coxeter  coset  associated  to  the  Coxeter  group  `W`,
`Frobenius(WF)`  returns a  function `F`  such that  `x↦ F(x,i=1)` does the
automorphism induced by `WF.phi^i` on the object `x`.

```julia-repl
julia> W=coxgroup(:D,4)
D₄

julia> WF=spets(W,Perm(1,2,4))
³D₄

julia> u=unichar(W,2)
[D₄]:<11->

julia> F=Frobenius(WF);F(u)
[D₄]:<.211>

julia> F(u,-1)
[D₄]:<11+>

julia> F(1)
4
```
"""
function Frobenius(WF::CoxeterCoset)
  (w,i=1)->Frobenius(w,WF.phi^i)
end

Frobenius(w::Perm,phi)=w^(inv(phi))
Frobenius(w::Integer,phi)=w^(inv(phi))

function twisting_elements(WF::Spets,J::AbstractVector{<:Integer})
  if isempty(J) return classreps(WF)./WF.phi end
  if !isone(WF.phi) error( "not implemented for twisted parent Spets" ) end
  twisting_elements(Group(WF),J)
# W=Group(WF)
# L=reflection_subgroup(W,J)
# N=normalizer(W,L)
# W_L=N/L
# if length(W_L)>=10
#   H=Group(map(x->central_action(L,reflrep(L,x.phi)),gens(W_L)))
#   if length(H)==length(W_L)
#     h=Hom(H,W_L,gens(W_L))
#     e=classreps(H)
#     return map(x->h(x).phi,e)
#   end
# end
# sort(map(x->x.phi,classreps(W_L)))
end

function twisting_elements(W::PermRootGroup,J::AbstractVector{<:Integer})
  if isempty(J) return classreps(W) end
  L=reflection_subgroup(W,J)
  N=normalizer(W,L)
#  W_L=C/N
  W_L=Group(gens(N))/Group(gens(L))
  map(x->x.phi,classreps(W_L))
end

function relative_coset(WF::CoxeterCoset,J)
# Print("CoxeterCosetOps.RelativeCoset ",WF,J," called \n");
  res=relative_group(Group(WF),J)
  spets(res,Perm(res.toparent,res.toparent.^WF.phi))
end
#-------------- subcoset ---------------------------------
"""
`subspets(WF,I,w=one(Group(WF)))`

Returns   the   reflection   subcoset   of   the   coset  `WF`  with  group
`reflection_subgroup(Group(WF),I)`  and torsion `w*WF.phi`.  `w` must be an
element  of `Group(WF)` such that  `w*WF.phi` normalizes the subroot system
generated by `I`.

```julia-repl
julia> WF=spets(coxgroup(:F,4))
F₄

julia> w=transporting_elt(Group(WF),[1,2,9,16],[1,9,16,2],ontuples);

julia> LF=subspets(WF,[1,2,9,16],w)
F₄₍₁‚₂‚₉‚₁₆₎=³D₄₍₃₄₁₂₎

julia> diagram(LF)
ϕ acts as (2,3,4) on the component below
  O 4
  ￨
O—O—O D₄
3 1 2
```
"""
function subspets(WF::Spets,I::AbstractVector{<:Integer},w=one(Group(WF));NC=false)
# xprintln("WF=",WF," I=",I," w=",w)
  W=Group(WF)
  if !(w in W) error(w," should be in ",W) end
  phi=w*WF.phi
  R=reflection_subgroup(W,I)
  if (W isa CoxeterGroup) &&
     (sort(action.(Ref(R),1:2nref(R),phi))!=1:2nref(R))
    error("w*WF.phi does not normalize subsystem")
  end
  phi=reduced(R,phi)
  if !(phi isa Perm) phi=phi.phi end
  if isone(phi) RF=spets(R)
  else RF=spets(phi,reflrep(WF,phi),R,Dict{Symbol,Any}())
  end
  if R!==W RF.parent=WF end
  if !NC refltype(RF) end
  RF
end

PermRoot.reflection_subgroup(WF::Spets,I::AbstractVector{<:Integer})=subspets(WF,I)

subspets(W::ComplexReflectionGroup,I::AbstractVector{<:Integer},w=one(W))=subspets(spets(W),I,w)

#-------------- spets ---------------------------------
@GapObj mutable struct PRC{T,TW<:PermRootGroup}<:Spets{T,TW}
  phi::T
  F::Matrix
  W::TW
end

function Base.show(io::IO,t::Type{PRC{T,TW}})where {T,TW}
  print(io,"Spets{",TW,"}")
end

function spets(W::PermRootGroup)
  get!(W,:trivialspets)do
    PRC(one(W),reflrep(W,one(W)),W,Dict{Symbol,Any}())
  end
end

function spets(W::PermRootGroup,w::Perm;NC=false)
  if isone(w) return spets(W) end
  w=reduced(W,w)
  if !(w isa Perm) w=w.phi end
  F=reflrep(W,w)
# println("w=$w\nF=$F")
  res=PRC(w,F,W,Dict{Symbol,Any}())
  if !NC refltype(res) end # changes phi and F
  res
end

function spets(W::PermRootGroup,F::AbstractMatrix;NC=false)
  w=PermX(W,F)
  if !isnothing(w) return spets(W,w;NC) end
# if W isa PRSG error("that's all for subgroups") end
# check if there exists a permutation perm and for each W-orbit of roots O 
# a scalar l_O such that (roots(W,O)*WF.F0Mat=l_O*roots(W,O.^perm)
  t=collectby(j->simple_reps(W)[j],eachindex(gens(W)))
  s=map(t)do inds
    scal=map(inds)do y
      r=transpose(F)*roots(W,y)
      scals=ratio.(Ref(r),roots(W))
      l=filter(i->!isnothing(scals[i]),eachindex(roots(W)))
      (ind=l,scal=scals[l])
    end
    if length(unique(map(x->unique(sort(x.scal)),scal)))>1 error("theory") end
    l=intersect(map(x->x.scal,scal)...)
    if isempty(l) return false end
    l=Root1.(l)
    if nothing in l
      error("only reflection cosets of finite order implemented")
      return false
    end
    # choose simplest scal
    l=minimum(map(x->[order(x),exponent(x)],l))
    l=E(l[1],l[2])
    [l,map(x->x.ind[findfirst(==(l),x.scal)],scal)]
  end
  if false in s
    ChevieErr("Spets(",W,",F=",F,
    " must normalize set of roots of parent up to scalars.\n")
    return false
  end
  scalars=map(x->x[1],s)[simple_reps(W)[eachindex(gens(W))]]
  perm=fill(0,length(roots(W)))
  for i in eachindex(t) perm[t[i]]=s[i][2] end
  while true
    i=filter(j->!iszero(perm[j]),eachindex(perm))
    if length(i)==length(perm) break end
    for j in eachindex(gens(W))
      perm[i.^refls(W,j)].=perm[i].^refls(W,perm[j])
    end
  end
  if length(unique(perm))<length(perm) return false end
  w=Perm(eltype(one(W)).(perm))
  if isnothing(w) error("matrix F must preserve the roots") end
  res=PRC(w,F,W,Dict{Symbol,Any}(:scalars=>scalars))
end

"""
spets(s::String)
builds a few of the exceptional spets
```julia-repl
julia> spets("3G422")
³G₄‚₂‚₂

julia> spets("2G5")
²G(ζ₃²√-2)₅

julia> spets("3G333")
G₃‚₃‚₃₍₁‚₂‚₃‚₄₄₎=³G₃‚₃‚₃₍₁‚₂‚₃‚₄₄₎

julia> spets("3pG333")
G₃‚₃‚₃₍₁‚₂‚₃‚₄₄₎=³G₃‚₃‚₃₍₁‚₂‚₃‚₄₄₎

julia> spets("4G333")
G₃‚₃‚₃₍₂‚₁₂‚₁₁‚₁₆‚₅₃‚₁₀‚₄₃‚₃₆₎=⁴G₃‚₃‚₃₍₁‚₂‚₃‚₃₂‚₁₆‚₃₆‚₃₀‚₁₀₎
```
"""
function spets(s::String;NC=false)
  if s=="3G422" 
    W=PRG([2 (root(3)-1)E(3);2 (-1+root(3))E(3,2);2 root(3)-1].//1,
     [(3+root(3))//2 root(3)E(3,2);(3+root(3))//2 root(3)E(3);(3+root(3))//2 root(3)].//3)
    return spets(W,reflrep(W,Perm(1,2,3));NC)
  elseif s=="2G5"  # reflection_subgroup(G14,[10,52])
    W=PRG([(-E(3)-2E(3,2))*(-3+root(6))//3 E(3);
           (-E(3)-2E(3,2))*(3-root(6))//3 E(3)],
          [E(3)//2 (-E(3)-2E(3,2))*(-3-root(6))//6;
           -E(3)//2 (-E(3)-2E(3,2))*(-3-root(6))//6])
    return spets(W,[-1 0;0 1];NC)
  elseif s=="3G333" 
    W=crg(3,3,3);return spets(reflection_subgroup(W,1:3),reflrep(W,Perm(1,2,44));NC)
  elseif s=="3pG333" 
    W=crg(3,3,3);return spets(reflection_subgroup(W,1:3),reflrep(W,Perm(1,44,2));NC)
  elseif s=="4G333" 
    W=crg(3,3,3)
    return spets(reflection_subgroup(W,1:3),
    perm"(1,44,32,37)(2,12,16,53)(3,50,30,15)
    ( 4,49,39,19)( 5, 9, 6,48)( 7,41,54,25)( 8,33,18,51)(10,43,36,11)
    (13,47,28,42)(14,52,22,17)(21,46,38,31)(24,27,26,40)";NC)
  elseif s=="3G422" 
    W=PRG([[2,(-1+root(3))*E(3)],[2,(-1+root(3))*E(3,2)],[2,(-1+root(3))]],
[[(3+root(3))//2,root(3)*E(3,2)],[(3+root(3))//2,root(3)*E(3)],[(3+root(3))//2,root(3)]]//3)
    return spets(W,[1 0;0 E(3)];NC)
  else error("argument should be 2G5, 3G422, 3G333, 3pG333 or 4G333")
  end
end

function relative_coset(WF::Spets,J=Int[],a...)
# xprintln("relative_coset(",WF,",",J,")");
  W=Group(WF)
  R=relative_group(W,J,a...)
  L=reflection_subgroup(W,J)
  res=spets(R,central_action(L,reflrep(L,WF.phi)))
  if isempty(J)
    Group(res).fromparent=R.fromparent
# else Group(res).fromparent:=function(x)Error("fromparent failed");return false;end;
  end
  res
end

Groups.NormalCoset(W::PermRootGroup,w::Perm=one(W))=spets(W,w)

spets(phi,F::Matrix,W::PermRootGroup,P::Dict{Symbol,Any})=PRC(phi,F,W,P)

Groups.Group(W::PRC)=W.W

isG333(t::TypeIrred)=haskey(t,:p) && (t.p,t.q,t.rank)==(3,3,3)

function PermRoot.refltype(WF::PRC)
  get!(WF,:refltype)do
    W=Group(WF)
    t=deepcopy(refltype(W))
    if isone(WF.phi)
      return map(x->TypeIrred(orbit=[x],twist=Perm()),t)
    end
    subgens=map(x->refls(W,x.indices),t)
    c=Perm(map(x->sort(x.^WF.phi),subgens),map(sort,subgens))
    if isnothing(c) && any(isG333,t)
      for a in t
        if isG333(a)
          a.subgroup=reflection_subgroup(W,a.indices;NC=true)
          c=chevieget(:timp,:ReducedInRightCoset)(a.subgroup,WF.phi)
          WF.phi=c.phi
          c=restriction(W,c.gen)
          a.indices=c
          a.subgroup=reflection_subgroup(W,c;NC=true)
          WF.W=reflection_subgroup(W,indices(t);NC=true)
          W=WF.W
          if order(WF.phi) in [3,6] G333=[1,2,3,44]
          elseif order(WF.phi)==4  G333=[1,2,3,32,16,36,30,10]
          else G333=1:3
          end
          a.subgens=refls(a.subgroup,G333)
          a.indices=inclusion(a.subgroup,W,G333)
          WF.W.refltype=t
        end
      end
  #   subgens=map(x->gens(reflection_subgroup(W,x.indices)),t)
      subgens=map(x->refls(W,x.indices),t)
      c=Perm(map(x->sort(x.^WF.phi),subgens),map(sort,subgens))
    end
    c=orbits(inv(c),eachindex(t))
    prr=roots(parent(W))
    function scals(rr,img)
      map(ratio,prr[inclusion(W,rr).^WF.phi],prr[inclusion(W,img)])
    end
    typ=map(c) do orb
      to=TypeIrred(orbit=map(copy,t[orb]))
      scalar=Cyc{Rational{Int}}[]
      for i in eachindex(orb)
        next=i==length(orb) ? 1 : i+1
        u=Perm(subgens[orb[next]],subgens[orb[i]].^WF.phi)
        tn=t[orb[next]]
        ti=t[orb[i]]
        if i!=length(orb)  tn.indices=invpermute(tn.indices,u)
          scal=scals(ti.indices,tn.indices)
        else to.twist=u
          scal=scals(ti.indices,invpermute(tn.indices,inv(u)))
        end
        if any(isnothing,scal) || !allequal(scal)
          if ti.series==:B && ti.rank==2
            scal=[root(prod(scal))]
          else println("W=$W\nphi=",WF.phi)
            error("scal=$scal :no element of coset acts as scalar on orbits")
            return nothing
          end
        end
        scal=Root1(scal[1])
        if !isone(scal)
        sub=reflection_subgroup(W,ti.indices;NC=true)
        zg=center(sub)
        z=length(zg)
 #      println("zg=$zg")
        if z>1 # simplify scalars using center
          ee=sort(elements(zg))
          zg=ee[findfirst(x->order(x)==z,ee)]
          i=inclusion(sub)[1]
          v=Root1(ratio(prr[i.^zg],prr[i]))
          zg^=invmod(exponent(v),order(v)) # distinguished
          v=scal.*E.(z,0:z-1)
          m=argmin(order.(v))
          Perms.mul!(WF.phi,(zg^(m-1))^WF.phi)
          WF.F=reflrep(Group(WF),WF.phi)
          scal=v[m]
        end
        # simplify again by -1 in types 2A(>1), 2D(odd), 2E6
        if mod(order(scal),4)==2 &&
          (ti.series in [:A,:D] || (ti.series==:E && ti.rank==6))
          sb=coxgroup(ti.series,ti.rank)
          w0=sub(word(sb,longest(sb))...)
          Perms.mul!(WF.phi,w0)
          u=Perm(subgens[orb[next]],subgens[orb[i]].^WF.phi)
          #print("l=$l i=$i scal before:$scal")
          if i!=length(orb)
            tn.indices^=u
            subgens[orb[next]]=refls(W,tn.indices)
            scal=scals(ti.indices,tn.indices)
          else to.twist=u
            scal=scals(ti.indices,invpermute(tn.indices,inv(u)))
          end
          #println(" scal after:$scal")
          if !allequal(scal) error(" scal after:$scal")
          else scal=Root1(scal[1])
          end
        end
        end
        push!(scalar,Cyc(scal))
      end
      to.scalar=scalar
      to
    end
  # some adjustment such that in type ^2D_4 the first two simple
  # roots are permuted by res.twist, and in type ^3D_4 the permutation 
  # of the simple roots is 1 -> 2 -> 4 -> 1:
  # in type 4G333 restore the indices to be in 1,2,3
    typ=map(typ)do a
      i=a.orbit[1].indices
      b=a.orbit[1]
      if b.series==:D && b.rank==4
        if order(a.twist)==2 # ^2D_4
          rf=filter(x->refls(W,i[x])!=refls(W,i[x])^WF.phi,1:4)
          rf=vcat(rf,[3],setdiff([1,2,4],rf))
          for j in a.orbit j.indices=j.indices[rf] end
        elseif order(a.twist)==3 # ^3D_4
          if refls(W,i[1])^WF.phi!=refls(W,i[2])
           for j in a.orbit j.indices=j.indices[[1,4,3,2]] end
          end
        end
      elseif isG333(b) # needless?
        if order(a.twist)==4
          for j in a.orbit j.indices[2]=inclusion(j.subgroup,W,2) end
        end
      end
      a
    end
  end
end

#--------------------- Root data ---------------------------------
Weyl.rootdatum(t::String,r::Int...)=rootdatum(Symbol(t),r...)

"""
`rootdatum(type::String or Symbol[,dimension or bond::Integer])` root datum
from type. The known types are

`2B2, 2E6, 2E6sc, 2F4, 2G2, 2I2, 3D4, 3D4sc, 3gpin8, CE6, CE7, E6, E6sc, E7,
E7sc, E8, F4, G2, cso, csp, gl, gpin, gpin-, halfspin, pgl, pso, pso-, psp,
psu,  ree, sl, slmod, so, so-, sp,  spin, spin-, su, suzuki, tgl, triality,
u`
"""
function Weyl.rootdatum(t::Symbol,r::Int...)
   if haskey(rootdata,t) 
     res=rootdata[t](r...) 
     res.callname="rootdatum("*repr(t)
     if !isempty(r) res.callname*=","*join(r,",") end
     res.callname*=")"
     n=string(t)
     sc=endswith(n,"sc")
     if sc n=n[1:end-2] end
     pre=match(r"^[0-9]+",n)
     if !isnothing(pre) res.TeXcallname="{}^{"*pre.match*"}"
     else res.TeXcallname=""
     end
     mid=match(r"[a-zA-Z-]+",n)
     if !isnothing(mid) res.TeXcallname*=mid.match end
     las=match(r"[0-9]+$",n)
     if !isnothing(las) res.TeXcallname*="_{"*las.match*"}" end
     if !isempty(r) 
       if mid.match=="I" res.TeXcallname*="("*join(r,",")*")"
       else res.TeXcallname*="_{"*join(r,",")*"}"
       end
     end
     if sc res.TeXcallname*="sc" end
     return res
   end
   error("unknown root datum $t. Known types are:\n",
         sprint(cut,join(sort(collect(keys(rootdata))),", ")))
end

id(r)=Matrix{Int}(I,r,r)


const  rootdata=Dict{Symbol,Function}()
rootdata[:gl]=function(r)
  if r==1 return torus(1) end
  R=id(r)
  R=R[1:end-1,:]-R[2:end,:]
  rootdatum(R,R)
end
rootdata[:sl]=r->rootdatum(cartan(:A,r-1),id(r-1))
rootdata[:slmod]=function(dim,q)
  if dim%q!=0 error(q," should divide ",dim) end
  intermediate_group(rootdata[:sl](dim),[q])
end
rootdata[:tgl]=function(n, k)
  if gcd(n,k)!=1 error(k," should be prime to ",n) end
  X=hcat(cartan(:A,n-1),fill(0,n-1))
  # We intertwine the last weight with the torus
  lat=vcat(X,transpose(vcat(fill(0,n-2), [k,1])))
  # Get a basis for the sublattice
  lat=hermite(lat)
  # Find the roots in a basis of the sublattice
  Y=id(n)[1:n-1,:]
  rootdatum(toM(map(x->solutionmatInt(lat,x),toL(X))),Y*transpose(lat))
end
rootdata[:pgl]=r->coxgroup(:A,r-1)
rootdata[:sp]=function(r)
  r=div(r,2)
  R=id(r)
  for i in 2:r R[i,i-1]=-1 end
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
rootdata[:cso]=function(dim)
  n=div(dim,2)
  R=fill(0,n,n+1)
  cR=fill(0,n,n+1)
  R[1,1:3]=[1,-1,-1]
  cR[1,2:3]=[-1,-1]
  for i in 2:n
    R[i,i:i+1]=[1,-1]
    cR[i,i:i+1]=[1,-1]
  end
  rootdatum(R,cR)
end
rootdata[:so]=function(r)
  r1=div(r,2)
  R=id(r1)
  for i in 2:r1 R[i,i-1]=-1 end
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
  R=cartan(:D,r)
  R[:,r].=0;R[1,r]=1
  rootdatum(Int.(transpose(inv(Rational.(R))*cartan(:D,r))),R)
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
rootdata[:CE6]=()->rootdatum([2 0 -1 0 0 0 1;0 2 0 -1 0 0 0;-1 0 2 -1 0 0 -1;
       0 -1 -1 2 -1 0 0;0 0 0 -1 2 -1 1;0 0 0 0 -1 2 -1],id(7)[1:6,:])
rootdata[Symbol("2CE6")]=()->spets(rootdatum(:CE6),[0 0 0 0 0 1 0;0 1 0 0 0 0 0;
     0 0 0 0 1 0 0;0 0 0 1 0 0 0;0 0 1 0 0 0 0;1 0 0 0 0 0 0;0 0 0 0 0 0 -1])
rootdata[:CE7]=()->rootdatum([0 2 -1 0 0 0 0 0;1 0 0 -1 0 0 0 0;
   0 -1 2 -1 0 0 0 0;-1 0 -1 2 -1 0 0 0;1 0 0 -1 2 -1 0 0;
   -1 0 0 0 -1 2 -1 0;1 0 0 0 0 -1 2 0],
   [0 1 0 0 0 0 0 0;2 0 0 0 -1 0 -1 -1;0 0 1 0 0 0 0 0;0 0 0 1 0 0 0 0;
    0 0 0 0 1 0 0 0;0 0 0 0 0 1 0 0;0 0 0 0 0 0 1 0])
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
  d=div(r,2)
  F=id(d+2)
  F[1,1:3]=[1,-1,1]
  F[1:d+2,2]=fill(-1,d+2)
  F[2:3,2:3]=-id(2)
  spets(rootdatum(:gpin,r),F)
end
rootdata[Symbol("so-")]=r->spets(rootdatum(:so,r),Perm(1,2))
rootdata[Symbol("pso-")]=r->spets(rootdatum(:pso,r),Perm(1,2))
rootdata[Symbol("spin-")]=r->spets(rootdatum(:spin,r),Perm(1,2))
rootdata[Symbol("2I2")]=e->spets(coxgroup(:Isym,2,e),Perm(1,2))
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
