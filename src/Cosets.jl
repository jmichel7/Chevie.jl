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
.Φ‴₈
```

A subgroup `H` which is a parabolic subgroup corresponds to a rational form
of  a Levi  subgroup of  `𝐆`. The  command 'Twistings'  gives all rational
forms of such a Levi.

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> twistings(W,[1])
2-element Array{Gapjm.Cosets.FCC{Int16,Gapjm.Weyl.FCSG{Int16,Int64,PRSG{Int64,Int16}}},1}:
 Ã₁Φ₁
 Ã₁Φ₂

julia> twistings(W,[2])
2-element Array{Gapjm.Cosets.FCC{Int16,Gapjm.Weyl.FCSG{Int16,Int64,PRSG{Int64,Int16}}},1}:
 A₁Φ₂
 A₁Φ₁
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
2-element Array{Gapjm.Cosets.FCC{Int16,Gapjm.Weyl.FCSG{Int16,Int64,PRSG{Int64,Int16}}},1}:
 (A₁A₁)
 A₁×A₁ 
```
"""
module Cosets

using ..Gapjm
export twistings, spets, Frobenius, Spets, torusfactors, subspets

abstract type Coset{TW} end

abstract type Spets{TW}<:Coset{TW} end

abstract type CoxeterCoset{TW}<:Spets{TW} end

function twisting_elements(W::FiniteCoxeterGroup,J::AbstractVector{<:Integer})
  if isempty(J) C=W
  elseif W isa CoxeterGroup && all(x->x in 1:coxrank(W),J)
    C=Group(collect(endomorphisms(CoxGroups.parabolic_category(W,J),1)))
  else C=centralizer(W,sort(J);action=(J,w)->sort(J.^w))
  end
  class_reps(C)
end

twistings(W::FiniteCoxeterGroup,J::AbstractVector{<:Integer})=
  spets.(Ref(reflection_subgroup(W,J)),twisting_elements(W,J))

#-------------- finite coxeter cosets ---------------------------------
struct FCC{T,TW<:FiniteCoxeterGroup{Perm{T}}}<:CoxeterCoset{TW}
  phi::Perm{T}
  F::Matrix
  W::TW
  prop::Dict{Symbol,Any}
end

spets(W::FiniteCoxeterGroup,w::Perm=Perm())=spets(W,matX(W,w))
Base.parent(W::FCC)=W

function spets(W::FiniteCoxeterGroup{Perm{T}},F::Matrix) where{T}
  perm=Perm{T}(F,roots(parent(W.G)),action=(r,m)->permutedims(m)*r)
  if isnothing(perm) error("matrix F must preserve the roots") end
  phi=reduced(W,perm)
  FCC(phi,F,W,Dict{Symbol,Any}())
end

function torusfactors(WF::CoxeterCoset)
  M=PermRoot.baseX(WF.W.G)
  M*=WF.F*inv(convert.(Cyc{Rational},M))
  r=length(gens(WF.W))
  M=M[r+1:end,r+1:end]
  if isempty(M) return CycPol(Pol([1],0)) end
  M=Ref(Pol([1],1)).*one(M)-M
  CycPol(GLinearAlgebra.Det(M))
end

function PermRoot.refltype(WF::CoxeterCoset)::Vector{TypeIrred}
  gets(WF,:refltype)do WF
    t=refltype(WF.W)
    c=map(x->PermRoot.indices(x),t)
    phires=Perm(WF.phi,inclusion(WF.W.G))
    map(orbits(Perm(sort.(c),map(i->sort(i.^phires),c))))do c
      o=deepcopy(t[c])
      J=PermRoot.indices(o[1])
      twist=Perm{Int16}(phires^length(c),J)
      if o[1][:series]==:D && length(J)==4
        if order(twist)==2
          rf=reduce(vcat,cycles(twist))
          o[1].prop[:indices]=J[vcat(rf,[3],setdiff([1,2,4],rf))]
        elseif order(twist)==3 && J[1]^twist!=J[2]
          o[1].prop[:indices]=J[[1,4,3,2]]
        end
      end
      for i in 2:length(c) 
        o[i].prop[:indices]=PermRoot.indices(o[i-1]).^phires
      end
      TypeIrred(Dict(:orbit=>o,:twist=>twist))
    end
  end
end

function Base.show(io::IO, W::CoxeterCoset)
   PermRoot.showtypes(io,refltype(W))
   t=torusfactors(W)
   if !isone(t) show(io,t) end
end

function Frobenius(WF::CoxeterCoset)
  f(w,i=1)=Frobenius(w,WF.phi^i)
end

Frobenius(w::Perm,phi)=w^phi

#-------------- finite coxeter subcoset ---------------------------------
struct FCSC{T,TW<:FiniteCoxeterGroup{Perm{T}}}<:CoxeterCoset{TW}
  phi::Perm{T}
  F::Matrix
  W::TW
  prop::Dict{Symbol,Any}
end

function subspets(WF,I::AbstractVector{Int},w=one(WF.W))
  RF=WF
  WF=parent(WF)
  phi=RF.phi/WF.phi
  W=WF.W
  if !(w in W) error(w," should be in ",W) end
  phi=w*phi
  R=reflection_subgroup(W,I)
  tmp=Set(inclusion(R))
  if Set(tmp.^(phi*WF.phi))!=tmp 
    error("w*WF.phi does not normalize subsystem")
  end
  phi=reduced(R,phi*WF.phi)
  FCSC(phi,matX(W,phi/WF.phi)*WF.F,R,Dict{Symbol,Any}(:parent=>WF))
end

end
