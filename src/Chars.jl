"""
Characters and conjugacy classes of complex reflection groups.

The  [`CharTable`](@ref)  of  a  finite  complex  reflection  group  `W` is
computed  using  the  decomposition  of  `W`  in  irreducible  groups,  see
[`refltype`](@ref).  For an irreducible group the character table is either
computed  using recursive formulas for the  infinite series, or read from a
library  file  for  the  exceptional  types.  Thus, character tables can be
obtained  quickly  even  for  very  large  groups  like `W(E₈)` or ``𝔖₁₅``.
Similar remarks apply for conjugacy classes.

The  conjugacy  classes  and  irreducible  characters of irreducible finite
complex  reflection groups  are labelled  by certain combinatorial objects;
these  labels are shown in  the tables we give.  For the classes, these are
partitions or partition tuples for the infinite series, or, for exceptional
Coxeter  groups,  Carter's  admissible  diagrams  [car72](@cite); for other
primitive  complex reflection groups we just use words in the generators to
specify  representatives of the  classes. For the  characters, these labels
are  again partitions or partition tuples  for the infinite series, and for
the  exceptional groups they are pairs of two integers `(d,b)` where `d` is
the  degree of the character and `b` is the smallest symmetric power of the
reflection  representation containing  the character  as a constituent (the
valuation of the [fakedegree](@ref)); further a `'` or `''` may be added to
disambiguate  pairs with  the same  `(d,b)`. These  labels are shown by the
functions   [`classinfo`](@ref)  and   [`charinfo`](@ref).  Displaying  the
character table also shows the labelings for classes and characters.

A  typical example is `coxgroup(:A,n)`,  the symmetric group ``𝔖ₙ₊₁`` where
classes  and characters are parameterized by  partitions of `n+1` (the same
table is shown for `coxsym(n+1)`).

```julia-repl
julia> CharTable(coxgroup(:A,3))
CharTable(A₃)
┌────┬─────────────────┐
│    │1111 211 22 31  4│
├────┼─────────────────┤
│1111│   1  -1  1  1 -1│
│211 │   3  -1 -1  .  1│
│22  │   2   .  2 -1  .│
│31  │   3   1 -1  . -1│
│4   │   1   1  1  1  1│
└────┴─────────────────┘

julia> ct=CharTable(coxgroup(:G,2))
CharTable(G₂)
┌─────┬────────────────────┐
│     │A₀ Ã₁ A₁ G₂ A₂ A₁+Ã₁│
├─────┼────────────────────┤
│φ₁‚₀ │ 1  1  1  1  1     1│
│φ₁‚₆ │ 1 -1 -1  1  1     1│
│φ′₁‚₃│ 1  1 -1 -1  1    -1│
│φ″₁‚₃│ 1 -1  1 -1  1    -1│
│φ₂‚₁ │ 2  .  .  1 -1    -2│
│φ₂‚₂ │ 2  .  . -1 -1     2│
└─────┴────────────────────┘

julia> charnames(ct,TeX=true)
6-element Vector{String}:
 "\\phi_{1,0}"
 "\\phi_{1,6}"
 "\\phi_{1,3}'"
 "\\phi_{1,3}''"
 "\\phi_{2,1}"
 "\\phi_{2,2}"

julia> classnames(ct,TeX=true)
6-element Vector{String}:
 "A_0"
 "\\tilde A_1"
 "A_1"
 "G_2"
 "A_2"
 "A_1+\\tilde A_1"
```

Reflection  groups  have  `fakedegrees`](@ref),  whose valuation and degree
give  two  invariants  `b,B`  for  each  irreducible  character of `W`. For
spetsial  groups (which include  finite Coxeter groups),  the valuation and
degree  of the generic degrees of the  spetsial Hecke algebra give two more
invariants  `a,A`; for Coxeter groups  see [car85; chap.11](@cite) for more
details.  These integers are used in the operations of truncated induction,
see [`j_induction_table`](@ref) and [`J_induction_table`](@ref).

Iwahori-Hecke  algebras and  cyclotomic Hecke  algebras also have character
tables, see [`HeckeAlgebras`](@ref).

We  now describe for each type our conventions for labeling the classes and
characters.

-  Type `Aₙ` (`n≥0`). 
We  have  ``W(Aₙ)≅𝔖ₙ₊₁``;  its  classes  and  characters  are  labelled  by
partitions  of `n+1`. The partition labelling a  class is the cycle type of
the  elements in that  class; the representative  word in the generators in
`.classtext`  is the concatenation of the words corresponding to each part,
where  the  word  for  a  part  `i`  is  the  product  of `i-1` consecutive
generators  (starting  one  higher  than  the  last  generator used for the
previous  parts). The partition labelling a character describes the type of
the  Young  subgroup  such  that  the  trivial  character induced from this
subgroup  contains that character with multiplicity `1` and such that every
other character occurring in this induced character has a higher `a`-value.
Thus,  the sign  character is  labelled by  the partition  `(1ⁿ⁺¹)` and the
trivial character by the partition `(n+1)`. The character of the reflection
representation of `W` is labelled by `(n,1)`.

-  Type  `Bₙ`  (`n≥2`).  
The group `W=W(Bₙ)` is isomorphic to the wreath product of the cyclic group
of  order  `2`  with  the  symmetric  group ``𝔖ₙ``. Hence  the  classes and
characters are parameterized by pairs of partitions such that the total sum
of  their parts equals `n`. The pair corresponding to a class describes the
signed  cycle type for the elements in that class, as in [car72](@cite). We
use  the convention that if `(λ,μ)` is  such a pair then `λ` corresponds to
the  positive  and  `μ`  to  the  negative  cycles.  Thus,  `((1ⁿ),-)`  and
`(-,(1ⁿ))`  label  respectively  the  trivial  class  and  the class of the
longest element.

The  pair  corresponding  to  an  irreducible  character  is determined via
Clifford  theory using the  semidirect product decomposition ``W(Bₙ)=N⋊𝔖ₙ``
where `N` is the standard `n`-dimensional `𝔽₂ⁿ`-vector space. For `a,b ≥ 0`
such that `n=a+b` let ``η_{a,b}`` be the irreducible character of `N` which
takes  value `1` on the first `a`  standard basis vectors and value `-1` on
the  last `b` standard basis  vectors of `N`. Then  the inertia subgroup of
``η_{a,b}``   has  the  form  ``T_{a,b}=N⋊(𝔖_a×𝔖_b)``  and  we  can  extend
``η_{a,b}``   trivially  to   an  irreducible   character  ``η̃_{a,b}``  of
``T_{a,b}``. Let `α` and `β` be partitions of `a` and `b`, respectively. We
take  the  tensor  product  of  the corresponding irreducible characters of
``𝔖_a``  and  ``𝔖_b``  and  regard  this  as  an  irreducible  character of
``T_{a,b}``.  Multiplying this character with  ``η̃_{a,b}`` and inducing to
`W(Bₙ)`  yields an irreducible character  ``χ= χ_{(α,β)}`` of `W(Bₙ)`. This
defines  the  correspondence  between  irreducible  characters and pairs of
partitions as above.

For example, the pair `((n),-)` labels the trivial character and `(-,(1ⁿ))`
labels  the  sign  character.  The  character  of  the  natural  reflection
representation is labeled by `((n-1),(1))`.

-  Type  `Dₙ` (`n≥4`). 
The  group `W(Dₙ)`  can be  embedded as  a subgroup  of index  `2` into the
Coxeter  group `W(Bₙ)`. The intersection of a class of `W(Bₙ)` with `W(Dₙ)`
is  either empty or a single class in `W(Dₙ)` or splits up into two classes
in `W(Dₙ)`. This also leads to a parameterization of the classes of `W(Dₙ)`
by  pairs of partitions `(λ,μ)` as before  but where the number of parts of
`μ`  is even and where there  are two classes of this  type if `μ` is empty
and all parts of `λ` are even. In the latter case we denote the two classes
in  `W(Dₙ)` by `(λ,+)`  and `(λ,-)`, where  we use the  convention that the
class  labeled by `(λ,+)` contains a representative which can be written as
a word in `{s₁,s₃,…,sₙ}` and `(λ,-)` contains a representative which can be
written as a word in `{s₂,s₃, …,sₙ}`.

By  Clifford theory the restriction of  an irreducible character of `W(Bₙ)`
to  `W(Dₙ)`  is  either  irreducible  or  splits  up  into  two irreducible
components.  Let `(α,β)` be  a pair of  partitions with total  sum of parts
equal to `n`. If `α!=β` then the restrictions of the irreducible characters
of  `W(Bₙ)` labeled  by `(α,β)`  and `(β,α)`  are irreducible and equal. If
`α=β`  then the restriction of the character labeled by `(α,α)` splits into
two  irreducible components  which we  denote by  `(α,+)` and `(α,-)`. Note
that  this can only happen if `n` is  even. In order to fix the notation we
use  a result of [ste89](@cite) which describes the value of the difference
of  these two  characters on  a class  of the  form `(λ,+)` in terms of the
character  values of  the symmetric  group ``𝔖_{n/2}``.  Recall that  it is
implicit  in the notation `(λ,+)` that all  parts of `λ` are even. Let `λ'`
be  the partition of `n/2` obtained by  dividing each part by `2`. Then the
value  of ``χ_{(α,-)}-χ_{(α,+)}``  on an  element in  the class  `(λ,+)` is
given  by  `2^{k(λ)}`  times  the  value  of  the  irreducible character of
``𝔖_{n/2}``  labeled by `α` on the class  of cycle type `λ'`. (Here, `k(λ)`
denotes the number of non-zero parts of `λ`.)

The  labels for the trivial, the  sign and the natural reflection character
are the same as for `W(Bₙ)`, since these characters are restrictions of the
corresponding characters of `W(Bₙ)`.

-  The groups `G(d,1,n)`.
They are isomorphic to the wreath product ``μ_d≀𝔖ₙ`` of the cyclic group of
order  `d` with the  symmetric group. Hence  the classes and characters are
parameterized  by `d`-tuples of partitions such that the total sum of their
parts  equals `n`. The words chosen  as representatives of the classes are,
when `d>2`, computed in a slightly different way than for `Bₙ`, in order to
agree  with the words on which Ram  and Halverson compute the characters of
the  Hecke algebra. First the parts of the `d` partitions are merged in one
big  partition and sorted in  increasing order. Then, to  a part `i` coming
from  the `j`-th partition is  associated the word `(l+1…1… l+1)ʲ⁻¹l+2…l+i`
where `l` is the highest generator used to express the previous part.

The  `d`-tuple corresponding to an  irreducible character is determined via
Clifford  theory in  a similar  way than  for the  `Bₙ` case.  The identity
character  has the first  partition with one  part equal `n`  and the other
ones  empty. The character of the  reflection representations has the first
two  partitions with one part  equal respectively to `n-1`  and to `1`, and
the other partitions empty.

-  The groups `G(de,e,n)`.
They  are normal  subgroups of  index `e`  in `G(de,1,n)`.  The quotient is
cyclic,  generated by the image `g`  of the first generator of `G(de,1,n)`.
The  classes are parameterized as the  classes of `G(de,e,n)` with an extra
information for a component of a class which splits.

According  to [hu85](@cite), a class `C`  of `G(de,1,n)` parameterized by a
`de`-partition  ``(S₀,…,S_{de-1})`` is in `G(de,e,n)` if `e` divides ``∑ᵢ i
∑_{p∈ Sᵢ}p``. It splits in `d` classes for the largest `d` dividing `e` and
all  parts of all `Sᵢ` and  such that `Sᵢ` is empty  if `d` does not divide
`i`.  If `w∈C` then `[gⁱ w g⁻ⁱ for i in 0:d-1]` are representatives of the
classes  of `G(de,e,n)` which meet `C`.  They are labelled by appending the
integer `i` to the label for `C`.

The  characters are described by Clifford theory. We make `g` act on labels
for  characters of `G(de,1,n)`  . The action  of `g` permutes circularly by
`d`  the partitions in the `de`-tuple.  A character has same restriction to
`G(de,e,n)`  as its transform by `g`.  The number of irreducible components
of its restriction is equal to the order `k` of its stabilizer under powers
of  `g`.  We  encode  a  character  of  `G(de,e,n)`  by first, choosing the
smallest  for lexicographical order label  of a character whose restriction
contains  it; then this label is periodic with a motive repeated `k` times;
we  represent the  character by  one of  these motives,  to which we append
`E(k)ⁱ` for `i in 0:k-1` to describe which component of the restriction we
choose.

-  Types `G₂` and `F₄`. 
A character is labeled by a pair `(d,b)` where `d` denotes the degree and `b`
the  corresponding `b`-invariant. If there  are several characters with the
same pair `(d,b)` we attach a prime to them, as in [car85](@cite).

The  matrices of  character values  and the  orderings and labelings of the
irreducible characters are exactly the same as in [car85;
p.412--413](@cite):  in type `G₂` the character  `φ₁,₃'` takes the value -1
on  the reflection associated  to the long  simple root; in  type `F₄`, the
characters  `φ₁,₁₂'`, `φ₂,₄'`,  `φ₄,₇'`, `φ₈,₉'`  and `φ₉,₆'`  occur in the
induced  of the  identity from  the `A₂`  corresponding to the short simple
roots;  the pairs (`φ₂,₁₆'`, `φ₂,₄″`) and (`φ₈,₃'`, `φ₈,₉″`) are related by
tensoring  by  sign;  and  finally  `φ₆,₆″`  is  the exterior square of the
reflection  representation. Note, however, that we put the long root at the
left  of the Dynkin  diagrams to be  in accordance with  the conventions in
[lus85; (4.8) and (4.10)](@cite).

The  classes are labeled by  Carter's admissible diagrams [car72](@cite).

-  Types  `E₆,E₇,E₈`. 
The  character tables are obtained by  specialization of those of the Hecke
algebra.   The  classes   are  labeled   by  Carter's  admissible  diagrams
[car72](@cite).  A  character  is  labeled  by  the  pair `(d,b)` where `d`
denotes  the degree and  `b` is the  corresponding `b`-invariant. For these
types, this gives a unique labeling of the characters.

- Non-crystallographic  types `I₂(m)`, `H₃`, `H₄`. 
In these cases we do not have canonical labelings for the classes. We label
them by reduced expressions.

Each  character for  type `H₃`  is uniquely  determined by the pair `(d,b)`
where  `d` is the degree and  `b` the corresponding `b`-invariant. For type
`H₄`  there are just  two characters (those  of degree `30`)  for which the
corresponding  pairs are  the same.  These two  characters are nevertheless
distinguished  by  their  fake  degrees:  the  character `φ₃₀,₁₀'` has fake
degree  `q¹⁰+q¹²+` higher terms, while `φ₃₀,₁₀″` has fake degree `q¹⁰+q¹⁴+`
higher  terms. The characters in the table for type `H₄` are ordered in the
same way as in [al82](@cite).

Finally,  the characters  of degree `2`  for type  `I₂(m)` are  ordered as
follows.  The matrix representations affording the characters of degree `2`
are given by: `` ρ_j : s₁s₂ ↦
\\begin{pmatrix}\\zeta_m^j&0\\\\0&\\zeta_m^{-j}\\end{pmatrix},
 s₁↦\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix},``
where  ``1 ≤ j ≤  ⌊(m-1)/2⌋``. The reflection representation
is  `ρ₁`. The  characters in  the table  are ordered by listing
first the characters of degree 1 and then `ρ₁,ρ₂,…`.

## Primitive complex reflection groups `G₄` to `G₃₄`.
The  groups `G₂₃=H₃`, `G₂₈=F₄`, `G₃₀=H₄` are exceptional Coxeter groups and
have  been  explained  above.  Similarly  for  the  other groups labels for
characters  consist primarily  of the  pair `(d,b)`  where `d`  denotes the
degree  and `b` is the corresponding  `b`-invariant. This is sufficient for
`G₄`,  `G₁₂`, `G₂₂` and `G₂₄`. For other  groups there are pairs or triples
of  characters which  have the  same `(d,b)`  value. We  disambiguate these
according to the conventions of [mal00](@cite) for `G₂₇, G₂₉, G₃₁, G₃₃` and
`G₃₄`:

-  For `G₂₇`:
The  fake degree  of `φ₃,₅'`  (resp. `φ₃,₂₀'`,  `φ₈,₉″`) has smaller degree
that  of  `φ₃,₅″`  (resp.  `φ₃,₂₀″`,  `φ₈,₉'`). The characters `φ₅,₁₅'` and
`φ₅,₆'` occur with multiplicity 1 in the induced from the trivial character
of  the parabolic subgroup  of type `A₂`  generated by the  first and third
generator  (it is asserted  mistakenly in [mal00](@cite)  that `φ₅,₆″` does
not occur in this induced; it occurs with multiplicity 2).

-  For `G₂₉`:
The  character  `φ₆,₁₀‴`  is  the  exterior  square  of `φ₄,₁`; its complex
conjugate  is `φ₆,₁₀⁗`. The  character `φ₁₅,₄″` occurs  in `φ₄,₁⊗φ₄,₃`; the
character  `φ₁₅,₁₂″`  is  tensored  by  the  sign  character from `φ₁₅,₄″`.
Finally  `φ₆,₁₀'` occurs in  the induced from  the trivial character of the
standard parabolic subgroup of type `A₃` generated by the first, second and
fourth generators.

-  For `G₃₁`:
The  characters `φ₁₅,₈'`, `φ₁₅,₂₀'` and `φ₄₅,₈″` occur in `φ₄,₁⊗φ₂₀,₇`; the
character   `φ₂₀,₁₃'`  is  complex  conjugate  of  `φ₂₀,₇`;  the  character
`φ₄₅,₁₂'`  is tensored by sign of `φ₄₅,₈'`. The two terms of maximal degree
of  the fakedegree of `φ₃₀,₁₀'` are  `q⁵⁰+q⁴⁶` while for `φ₃₀,₁₀″` they are
`q⁵⁰+2q⁴⁶`.

-  For `G₃₃`:
The two terms of maximal degree of the fakedegree of `φ₁₀,₈'` are `q²⁸+q²⁶`
while  for `φ₁₀,₈″` they are `q²⁸+q²⁴`. The  terms of maximal degree of the
fakedegree   of  `φ₄₀,₅'`  are  `q³¹+q²⁹`   while  for  `φ₄₀,₅″`  they  are
`q³¹+2q²⁹`.  The character  `φ₁₀,₁₇'` is  tensored by  sign of `φ₁₀,₈'` and
`φ₄₀,₁₄'` is tensored by sign of `φ₄₀,₅'`.

-  For `G₃₄`:
The  character `φ₂₀,₃₃'` occurs in `φ₆,₁⊗φ₁₅,₁₄`. The character `φ₇₀,₉'` is
rational.  The character  `φ₇₀,₉″` occurs  in `φ₆,₁⊗φ₁₅,₁₄`.  The character
`φ₇₀,₄₅'`   is  rational.The   character  `φ₇₀,₄₅″`   is  tensored  by  the
determinant  character of  `φ₇₀,₉″`. The  character `φ₅₆₀,₁₈'` is rational.
The character `φ₅₆₀,₁₈‴` occurs in `φ₆,₁⊗φ₃₃₆,₁₇`. The character `φ₂₈₀,₁₂'`
occurs    in   `φ₆,₁⊗φ₃₃₆,₁₇`.   The   character   `φ₂₈₀,₃₀″`   occurs   in
`φ₆,₁⊗φ₃₃₆,₁₇`.  The  character  `φ₅₄₀,₂₁'`  occurs  in `φ₆,₁⊗φ₁₀₅,₂₀`. The
character  `φ₁₀₅,₈'` is  complex conjugate  of `φ₁₀₅,₄`,  and `φ₈₄₀,₁₃'` is
complex  conjugate  of  `φ₈₄₀,₁₁`.  The  character  `φ₈₄₀,₂₃'`  is  complex
conjugate  of  `φ₈₄₀,₁₉`.  Finally  `φ₁₂₀,₂₁'`  occurs  in induced from the
trivial character of the standard parabolic subgroup of type `A₅` generated
by the generators of `G₃₄` with the third one omitted.

For  the groups `G₅` and `G₇` we  adopt the following conventions. For `G₅`
they are compatible with those of [mr03](@cite) and [bmm14](@cite).

-  For `G₅`:
We  let `W=complex_reflection_group(5)`,  so the  generators are  `W(1)` and
`W(2)`.

The  character `φ₁,₄'` (resp. `φ₁,₁₂'`, `φ₂,₃'`) takes the value `1` (resp.
`ζ₃`,  `-ζ₃`)  on  `W(1)`.  The  character  `φ₁,₈″` is complex conjugate to
`φ₁,₁₆`,  and the character  `φ₁,₈'` is complex  conjugate to `φ₁,₄'` . The
character  `φ₂,₅″` is complex conjugate to  `φ₂,₁`; `φ₂,₅'` takes the value
`-1` on `W(1)`. The character `φ₂,₇'` is complex conjugate to `φ₂,₅'`.

-  For `G₇`:
We  let `W=complex_reflection_group(7)`,  so the  generators are
`W(1)`, `W(2)` and `W(3)`.

The  characters  `φ₁,₄'`  and  `φ₁,₁₀'`  take  the value `1` on `W(2)`. The
character  `φ₁,₈″` is complex  conjugate to `φ₁,₁₆`  and `φ₁,₈'` is complex
conjugate  to `φ₁,₄'`. The characters `φ₁,₁₂'`  and `φ₁,₁₈'` take the value
`ζ₃`  on `W(2)`. The character `φ₁,₁₄″` is complex conjugate to `φ₁,₂₂` and
`φ₁,₁₄'`  is complex conjugate to `φ₁,₁₀'`. The character `φ₂,₃'` takes the
value  `-ζ₃` on  `W(2)` and  `φ₂,₁₃'` takes  the value  `-1` on `W(2)`. The
characters  `φ₂,₁₁″`, `φ₂,₅″`, `φ₂,₇‴` and  `φ₂,₁` are Galois conjugate, as
well  as  the  characters  `φ₂,₇'`,  `φ₂,₁₃'`,  `φ₂,₁₁'`  and  `φ₂,₅'`. The
character  `φ₂,₉'` is complex  conjugate to `φ₂,₁₅`  and `φ₂,₉‴` is complex
conjugate to `φ₂,₃'`.

Finally,  for the remaining groups `G₆, G₈`  to `G₁₁, G₁₃` to `G₂₁`, `G₂₅`,
`G₂₆`,  `G₃₂` and `G₃₃` there are only  pairs of characters with same value
`(d,b)`.  We give labels uniformly to these characters by applying in order
the following rules :

-  If the two characters have  different fake degrees, label `φ_{d,b}'` the
   one  whose  fake  degree  is  minimal  for  the  lexicographic  order of
   polynomials (starting with the highest term).

-  For the not yet labeled pairs, if only one of the two characters has the
   property   that  in  its   Galois  orbit  at   least  one  character  is
   distinguished by its `(d,b)`-invariant, label it `φ_{d,b}'`.

-  For the not yet labeled pairs,  if the minimum of the `(d,b)`-value (for
   the  lexicographic  order  `(d,b)`)  in  the  Galois  orbits  of the two
   character  is different, label `φ_{d,b}'` the character with the minimal
   minimum.

-  We define now a new invariant `t` for a character `φ`: consider all the
   pairs of irreducible characters `χ` and `ψ` uniquely determined by their
   `(d,b)`-invariant such that `φ` occurs with non-zero multiplicity `m` in
   `χ⊗ψ`.  We define  `t(φ)` to  be the  minimal (for  lexicographic order)
   possible list `(d(χ),b(χ),d(ψ),b(ψ),m)`.

For  the not  yet labeled  pairs, if  the t-invariants are different, label
`φ_{d,b}'` the character with the minimal `t`-invariant.

After  applying  the  last  rule  all  the  pairs  will be labelled for the
considered  groups. The labelling obtained  is compatible for `G₂₅`, `G₂₆`,
`G₃₂`  and  `G₃₃`  with  that  of  [mal00](@cite)  and  for  `G₈` with that
described in [mr03](@cite).

We  should  emphasize  that  for  all  groups  with  a  few exceptions, the
parameters  for characters do  not depend on  any non-canonical choice. The
exceptions  are `G(de,e,n)` with `e>1`, and `G₅`, `G₇`, `G₂₇`, `G₂₈`, `G₂₉`
and  `G₃₄`, groups  which admit  outer automorphisms  preserving the set of
reflections,  so choices  of a  particular value  on a particular generator
must be made for characters which are not invariant by these automorphisms.

-  Labels  for the classes. 
For  the exceptional complex reflection groups,  the labels for the classes
represent  the decomposition of a representative  of the class as a product
of  generators,  with  the  additional  conventions that `z` represents the
generator  of the  center and  for well-generated  groups `c`  represents a
Coxeter element (a product of the generators which is a regular element for
the highest reflection degree).
"""
module Chars

using ..Chevie

export charinfo, classinfo, fakedegrees, CharTable, representation,
  WGraphToRepresentation, DualWGraph, WGraph2Representation, charnames,
  representations, InductionTable, induction_table, classes, j_induction_table,
  J_induction_table, decompose, on_chars, detPerm, conjPerm,
  classnames, decomposition_matrix, eigen, schur_functor, charnumbers

"""
`schur_functor(mat,λ)`

`mat`  should be  a square  matrix and  `λ` a  partition. The result is the
Schur  functor  of  the  matrix  `mat`  corresponding to partition `λ`; for
example,   if  `λ==[n]`  it  returns  the   n-th  symmetric  power  and  if
`λ==[1,1,1]` it returns the 3rd exterior power. The current algorithm (from
Littlewood)  is rather inefficient so it is  quite slow for partitions of n
where `n>6`.

```julia-repl
julia> m=cartan(:A,3)
3×3 Matrix{Int64}:
  2  -1   0
 -1   2  -1
  0  -1   2

julia> schur_functor(m,[2,2])
6×6 Matrix{Rational{Int64}}:
   9   -6    4  3//2   -2    1
 -12   16  -16  -4      8   -4
   4   -8   16   2     -8    4
  12  -16   16  10    -16   12
  -4    8  -16  -4     16  -12
   1   -2    4  3//2   -6    9
```
"""
function schur_functor(A::AbstractMatrix,λ)
  n=sum(λ)
  Sn=coxsym(n)
  p=findfirst(==(λ),partitions(n))
# r=representation(Sn,p)
# does not work with the above repr but works with the below one. Fishy!
  r=chevieget(:imp,:Representation)(1,1,n,p)[2:end]
  rep(x)=isone(x) ? one(r[1]) : prod(r[word(Sn,x)])
  f(μ)=prod(factorial,last.(tally(μ)))
  basis=multisets(axes(A,1),n)
  basis=filter(b->!isempty(semistandard_tableaux(λ,b)),basis)
  if isempty(basis) return A[1:0,1:0] end
  M=sum(x->kron(rep(x),
    toM(map(basis)do i
      i=invpermute(i,x)
      map(j->prod(k->A[i[k],j[k]],1:n),basis)//f(i)
     end)),elements(Sn))
  #print(size(M),"=>")
  M=M[filter(i->!iszero(@view M[i,:]),axes(M,1)),
      filter(i->!iszero(@view M[:,i]),axes(M,2))]
  m=sort.(collectby(i->M[:,i],axes(M,2)))
  m=sort(m)
  M=M[:,first.(m)]
  M=map(x->sum(M[x,:],dims=1)[1,:],m)
  M=improve_type(toM(M))
  #println(size(M))
  #M
end

"`dim_schur(n,λ)` the dimension of `schur_functor(M,λ)` for `M` of size (n,n)"
function dim_schur(n,λ)
  l(i)=i>length(λ) ? 0 : λ[i]
  Integer(prod(1+(l(i)-l(j))//(j-i) for i in 2:n for j in 1:i-1))
end

"""
`fakedegree(W, φ, q=Pol())`

returns the fake degree (see [`fakedegrees`](@ref) for a definition) of the
character  of parameter  φ (see  [`charinfo`](@ref)`(W).charparams`) of the
reflection group `W`, evaluated at `q` .

```julia-repl
julia> fakedegree(coxgroup(:A,2),[[2,1]],Pol(:q))
Pol{Int64}: q²+q
```
"""
function Symbols.fakedegree(W,p,q=Pol())
  typ=refltype(W)
  if isempty(typ) return one(q) end
# prod(map((t,p)->fakedegree(t,p,q),typ,p))
  prod(fakedegree.(typ,p,q))
end

function Symbols.fakedegree(t::TypeIrred,p,q=Pol())
  if haskey(t,:scalar) q=prod(s->q*conj(s),t.scalar)
  elseif haskey(t,:orbit) q=q^length(t.orbit)
  end
  chevieget(t,:FakeDegree,p,q)
end

"""
`fakedegrees(W, q=Pol())`

returns  a list holding the fake degrees of the reflection group `W` on the
vector  space `V`, evaluated at `q`. These are the graded multiplicities of
the  irreducible characters of `W` in the quotient `SV/I` where `SV` is the
symmetric  algebra of `V` and `I` is the ideal generated by the homogeneous
invariants  of  positive  degree  in  `SV`.  The  ordering  of  the  result
corresponds to the ordering of the characters in `charinfo(W)`.

```julia-repl
julia> fakedegrees(coxgroup(:A,2),Pol(:q))
3-element Vector{Pol{Int64}}:
 q³
 q²+q
 1
```
"""
function fakedegrees(W,q=Pol();recompute=false)
  if !recompute
    res=get!(W,:fakedegrees)do
      improve_type(map(p->fakedegree(W,p,Pol()),charinfo(W).charparams))
    end
    if !any(isnothing,res) && !all(iszero,res)
      if q==Pol() return res
      else return map(x->x(q),res)
      end
    end
  end
  # recompute from general principles
  InfoChevie("# recomputing fakedegrees for ",W,"\n")
  qq=Pol()
  P=generic_order(W,qq)
  P=shift(P,-valuation(P))
  ct=CharTable(W)
  P=ct.irr*map(enumerate(ct.centralizers))do (i,c)
    exactdiv(P,
        improve_type(prod(l->(qq*conj(l)-1),refleigen(W,i);init=one(qq))))//c
  end
  if q!=qq P=map(x->x(q),P) end
  P=improve_type(P)
  if W isa Spets P.*=(-1)^rank(W)*generic_sign(W) end
  P
end

@GapObj struct CharInfo
end

charnumbers(d::Dict)=d[:charNumbers]

function charinfo(t::TypeIrred)
  c=CharInfo(copy(chevieget(t,:charinfo)))
  c.positionId=c.extRefl[1]
  c.positionDet=c.extRefl[end]
  if !haskey(c,:charnames) error("charnames(",t,") missing") end
  if !haskey(c,:b) c.b=chevieget(t,:b) end
  if !haskey(c,:B) c.B=chevieget(t,:B) end
  if !haskey(c,:a) c.a=chevieget(t,:a) end
  if !haskey(c,:A) c.A=chevieget(t,:A) end
  if isnothing(c.a)
    uc=chevieget(t,:UnipotentCharacters)
    if !isnothing(uc) && uc!=false
#     printstyled(rio(),"using UniChars\n";color=:red)
      if haskey(uc,:almostHarishChandra)
        c.a=uc[:a][charnumbers(uc[:almostHarishChandra][1])]
        c.A=uc[:A][charnumbers(uc[:almostHarishChandra][1])]
      else
        c.a=uc[:a][charnumbers(uc[:harishChandra][1])]
        c.A=uc[:A][charnumbers(uc[:harishChandra][1])]
      end
    elseif !haskey(t,:orbit)
      para=ordergens(t)
      if !isnothing(para)
        para=map(x->vcat([Mvp(:x)],map(j->E(x,j),1:x-1)),para)
        s=map(p->chevieget(t,:SchurElement,p,para,fill(missing,t.rank+1)),c.charparams)
        c.a=-valuation.(s)
        c.a.+=valuation(s[c.positionId])
        c.A=-degree.(s)
        c.A.+=degree(s[c.positionId])
      end
    end
  end
  for f in [:a,:A,:b,:B]
    if haskey(c,f)
      if isnothing(c[f]) || c[f]==false delete!(c.prop,f)
      else c.prop[f]=improve_type(c[f]) end
    end
  end
  if haskey(t,:orbit)
    if !haskey(c,:charRestrictions)
      c.charRestrictions=eachindex(c.charparams)
      c.nrGroupClasses=length(c.charparams) # assume orbit twist trivial
    end
    for f in [:a,:A,:b,:B]
      if haskey(c,f) c.prop[f]*=length(t.orbit) end
    end
  end
  c
end

cartfields(p,f)=cartesian(getproperty.(p,f)...)
tcartfields(p,f)=tcartesian(getproperty.(p,f)...)

"""
`charinfo(W)`

returns   information  about  the  irreducible  characters  of  the  finite
reflection group or Spets `W`. The result is an object with various entries
describing  properties of  the irreducible  characters of  `W`. This object
prints  at the  Repl or  in Pluto  or Jupyter  as a table synthesizing most
information.

A  field  not  printed  is  `.charparams`:  it  contains parameters for the
irreducible  characters.  A  parameter  is  a  list  with one item for each
irreducible  component of `W` (as given  by `refltype`). For an irreducible
`W` see the explanation in [`Chars`](@ref) for what are the parameters.

```julia-repl
julia> charinfo(coxgroup(:G,2)).charparams
6-element Vector{Vector{Vector{Int64}}}:
 [[1, 0]]
 [[1, 6]]
 [[1, 3, 1]]
 [[1, 3, 2]]
 [[2, 1]]
 [[2, 2]]
```

```julia-repl
julia> charinfo(coxgroup(:G,2))
┌──┬──────────────────────────────────────────────────────────┐
│n0│ name ext b B a A spaltenstein lusztig              symbol│
├──┼──────────────────────────────────────────────────────────┤
│1 │ φ₁‚₀  Id 0 0 0 0            1       1       (0,0,0,0,0,2)│
│2 │ φ₁‚₆ det 6 6 6 6            ε       ε (01,01,01,01,01,12)│
│3 │φ′₁‚₃     3 3 1 5           εₗ      ε′            (0,0,1+)│
│4 │φ″₁‚₃     3 3 1 5          ε_c      ε″            (0,0,1-)│
│5 │ φ₂‚₁  Λ¹ 1 5 1 5           θ′      θ′       (0,0,0,0,1,1)│
│6 │ φ₂‚₂     2 4 1 5           θ″      θ″       (0,0,0,1,0,1)│
└──┴──────────────────────────────────────────────────────────┘
```
In  the table printed  at the Repl,  the columns reflect  various fields of
`charinfo`.  The  column  `name`  reflects  the  field `.charnames`, a name
computed  from `.charparams`. This  is the same  as `charnames(io,W)` where
here `io` being the Repl has the property `:limit` true.

The   column   `ext`   shows   the   exterior   powers  of  the  reflection
representation.  It corresponds  to the  field `.extrefl`  which is present
only  if `W`  is irreducible.  Otherwise, only  two items  are shown in the
column:  `Id` corresponds to the field  `.positionId` and shows the trivial
character.  `det`  corresponds  to  the  field `.positionDet` and shows the
determinant  character (for Coxeter groups the sign character). When `W` is
irreducible,  the characters marked  `Λⁱ` are the  `i`-th exterior power of
the  reflection  representation.  They  are  irreducible  by  a  theorem of
Steinberg.

The  column  `b`  shows  the  field  `.b`  listing  for  each character the
valuation  of the fake degree, and the column `B` shows the field `.B`, the
degree of the fake degree.

The  columns `a` and  `A` only appear  for Spetsial groups. They correspond
then  to the fields  `.a` and `.A`,  and contain respectively the valuation
and the degree of the generic degree of the character (in the one-parameter
Hecke algebra `hecke(W,Pol())` for `W`).

For  irreducible  groups,  the  table  shows  sometimes additional columns,
corresponding to a field of the same name.

for  `F₄`,  the  column  `kondo`  gives  the  labeling of the characters in
[kon65](@cite). The columns `spaltenstein` and `lusztig` show variants used
in [spalt85](@cite) and [lus85; (4.10)](@cite).

for `E₆, E₇, E₈` the column `frame` gives the labeling of the characters in
[frame51,frame70](@cite), also used in [lus85; (4.11) to (4.13)](@cite).

for  `G₂` the  column `spaltenstein`  gives the  labeling of the characters
used  by [spalt85](@cite) and the column  `lusztig` shows a variant used in
[lus85](@cite).

for  `H₃` and `H₄` the  column `gp`  gives the  labeling of the characters
used  in [gp2000](@cite).

for `G(de,e,2)` even `e` and `d>1`, the column `malle` gives the parameters
for the characters used in [mal96](@cite).

If  `W`  is  irreducible  spetsial  and  imprimitive,  the  column `symbol`
(corresponding  to the field  `.charSymbols`) shows the  symbol attached to
the corresponding unipotent caracter.

If  `W  isa  Spets`,  the  column  `restr.`  (corresponding  to  the  field
`.charRestrictions`)  gives  the  index  of  the corresponding character of
`Group(W)`.

Finally,  the  field  `.hgal`  contains  the  permutation of the characters
resulting  from a Galois  action on the  characters of `H=hecke(W,Pol()^e)`
where  `e` is the order of  the center of `W`. `H`  splits by taking `v` an
`e`-th root of `Pol()`, and `.hgal` records the permutation effected by the
Galois action `v->E(e)*v` (`charinfo` does not have the key `:hgal` if this
permutation   is  trivial).  `.hgal*conj`,  where  `conj`  is  the  complex
conjugaison, is the Opdam involution.
```julia-repl
julia> charinfo(complex_reflection_group(24))
┌──┬─────────────────────┐
│n0│ name ext  b  B  a  A│
├──┼─────────────────────┤
│1 │ φ₁‚₀  Id  0  0  0  0│
│2 │φ₁‚₂₁ det 21 21 21 21│
│3 │ φ₃‚₈      8 18  8 20│
│4 │ φ₃‚₁  Λ¹  1 11  1 13│
│5 │φ₃‚₁₀  Λ² 10 20  8 20│
│6 │ φ₃‚₃      3 13  1 13│
│7 │ φ₆‚₂      2 12  1 13│
│8 │ φ₆‚₉      9 19  8 20│
│9 │ φ₇‚₆      6 18  6 18│
│10│ φ₇‚₃      3 15  3 15│
│11│ φ₈‚₄      4 16  4 17│
│12│ φ₈‚₅      5 17  4 17│
└──┴─────────────────────┘
hgal=(11,12)
```
"""
function charinfo(W)
  get!(W,:charinfo)do
    p=charinfo.(refltype(W))
    if isempty(p)
      res=CharInfo(Dict(:a=>[0],:A=>[0],:b=>[0],:B=>[0],:positionId=>1,
        :positionDet=>1,:charnames=>["Id"],:extRefl=>[1],:charparams=>[[]]))
      if W isa Spets
        res.charRestrictions=[1]
        res.nrGroupClasses=1
      end
      return res
    end
    if length(p)==1 res=p[1]
    else res=CharInfo(Dict{Symbol, Any}())
    end
    res.charparams=cartfields(p,:charparams)
    if W isa Spets
      gt=map(x->sort(x.indices),refltype(Group(W)))
      t=refltype(W)
      n=fill(0,length(gt))
      for i in eachindex(t), f in t[i].orbit
        n[findfirst(==(sort(f.indices)),gt)]=p[i].nrGroupClasses
      end
      res.charRestrictions=map(cartfields(p,:charRestrictions))do y
        m=fill(0,length(gt))
        for i in eachindex(t), f in t[i].orbit
          m[findfirst(==(sort(f.indices)),gt)]=y[i]
        end
        cart2lin(n,m)
      end
      res.nrGroupClasses=prod(i->p[i].nrGroupClasses^length(t[i].orbit),
                                                          eachindex(t))
    end
    if length(p)>1 res.charnames=join.(cartfields(p,:charnames),",") end
    for f in [:positionId, :positionDet]
      if all(d->haskey(d,f),p)
        res.prop[f]=cart2lin(map(x->length(x.charparams),p),getindex.(p,f))
      end
    end
    for f in [:b, :B, :a, :A]
      if all(d->haskey(d,f),p) 
        res.prop[f]=improve_type(sum.(tcartfields(p,f)))
      end
    end
    if !haskey(res,:b) || !haskey(res,:B)
       P=fakedegrees(W;recompute=true)
       res.b=valuation.(P)
       res.B=degree.(P)
    end
    if any(x->haskey(x,:hgal),p)
      res.hgal=map(x->haskey(x,:hgal) ? x.hgal : Perm(), p)
      gt=cartesian(map(x->1:length(x.charparams), p)...)
      res.hgal=Perm(gt, map(t->t.^res.hgal,gt))
    end
    res
  end::CharInfo
end

charinfo(W::Weyl.FC)=charinfo(W.G)

using LaTeXStrings
function Base.show(io::IO, ::MIME"text/html", ci::CharInfo)
  print(io,latexstring(xrepr(MIME("text/plain"),ci,TeX=true)))
end

function Base.show(io::IO,t::CharInfo)
  print(io,"CharInfo(",t.prop,")")
end

function Base.show(io::IO, ::MIME"text/plain", ci::CharInfo)
  if !hasdecor(io)
    print(io,"CharInfo(",ci.prop,")")
    return
  end
  n=length(ci.charnames)
  t=fromTeX.(Ref(io),ci.charnames);cl=String["name"]
  ext=fill("",n)
  if haskey(ci,:extRefl)
    ext[ci.extRefl[[1,end]]]=["Id","det"]
    for (i,j) in pairs(ci.extRefl[2:end-1])
      ext[j]=fromTeX(io,"\\Lambda^$i")
    end
  else ext[ci.positionId]="Id";ext[ci.positionDet]="det"
  end
  t=hcat(t,ext); push!(cl,"ext")
  for key in [:b,:B,:a,:A,:kondo,:spaltenstein,:frame,:malle,:lusztig,:gp]
    if haskey(ci,key) && ci[key]!=false
      t=hcat(t,fromTeX.(Ref(io),string.(ci[key])))
      push!(cl,string(key))
    end
  end
  if haskey(ci,:charRestrictions)
    t=hcat(t,string.(ci.charRestrictions))
    push!(cl,"restr.")
  end
  if haskey(ci,:charSymbols)
    t=hcat(t,repr.(ci.charSymbols,context=io));push!(cl,"symbol")
  end
  showtable(io,string.(t);row_labels=string.(1:n),col_labels=cl,rows_label="n0")
  if haskey(ci,:hgal) println(io,"hgal=",ci.hgal) end
end

"""
`detPerm(W)`

return  the permutation of the characters of the reflection group `W` which
is effected when tensoring by the determinant character (for Coxeter groups
this is the sign character).

```julia-repl
julia> W=coxgroup(:D,4)
D₄

julia> detPerm(W)
(1,8)(2,9)(3,11)(4,13)(7,12)
```
"""
function detPerm(W)
  get!(W,:detPerm)do
    t=CharTable(W).irr
    Perm(t,t.*transpose(t[charinfo(W).positionDet,:]);dims=1)
  end::Perm{Perms.Idef}
end

"""
`conjPerm(W)`

return  the permutation of the characters of the group `W` which
is effected when taking the complex conjugate of the character table.

```julia-repl
julia> W=complex_reflection_group(4)
G₄

julia> conjPerm(W)
(2,3)(5,6)
```
"""
function conjPerm(W)
  get!(W,:conjPerm)do
    t=CharTable(W).irr
    Perm(t,conj.(t);dims=1)
  end::Perm{Perms.Idef}
end

@GapObj struct ClassInfo
end

function classinfo(t::TypeIrred)
  cl=copy(convert(Dict{Symbol,Any},chevieget(t,:classinfo)))
  if haskey(t,:orbit)
     l=length(t.orbit)
     t=t.orbit[1]
     if l>1 && haskey(cl,:classes)
       cl[:classes]=cl[:classes].*prod(degrees(t))^(l-1)
     end
  end
  inds=t.indices
  cl[:classtext]=map(x->inds[x],cl[:classtext])
  if haskey(cl,:classes) cl[:classes]=Int.(cl[:classes]) end
  if haskey(cl,:centralizers) cl[:centralizers]=Int.(cl[:centralizers]) end
  cl[:classnames]=String.(cl[:classnames])
  if !haskey(cl,:classparams) cl[:classparams]=cl[:classtext] end
  ClassInfo(cl)
end

Groups.nconjugacy_classes(t::TypeIrred)=chevieget(t,:nconjugacy_classes)

"""
`classinfo(W)`

returns  information about the  conjugacy classes of  the finite reflection
group or Spets `W`. The result is an object with various entries describing
properties  of the conjugacy classes of `W`. This object prints at the Repl
or in Pluto or Jupyter as a table synthesizing most information.

A field not shown in the table is `.classparams`, containing parameters for
the  conjugacy classes. Each parameter  is a vector which  has one item for
each  irreducible  component  of  `W`.  For  what are the parameters for an
irreducible `W`, see the description for each type in [`Chars`](@ref).

```julia-repl
julia> classinfo(coxgroup(:A,2))
┌──┬──────────────────────┐
│n0│name length order word│
├──┼──────────────────────┤
│1 │ 111      1     1    .│
│2 │  21      3     2    1│
│3 │   3      2     3   12│
└──┴──────────────────────┘
```
The table contains the columns:
  - `name`, corresponding to the field `.classnames`:  strings describing the
    conjugacy classes, made out of the information in `:classparams`.
  - `length`, corresponding to the field `.classes`, is the number of elements
    in the conjugacy class.
  - `order`, corresponding to the field `.orders`, is the order of elements
    in the conjugacy class.
  - `word`, corresponding to the field `.classtext`, describes a word in
    the  generators for the  representatives of each  conjugacy class. Each
    word is a list of integers where the generator `W(i)` is represented by
    the  integer  `i`.  For  finite  Coxeter  groups,  it  is  the  same as
    `word.(Ref(W),classreps(W))`,   and  each  such  representative  is  of
    minimal  length in its conjugacy class and  is a "very good" element in
    the sense of [gm97](@cite).
"""
function classinfo(W)
  get!(W,:classinfo)do
    tmp=map(classinfo,refltype(W))
    if isempty(tmp) return ClassInfo(Dict(:classtext=>[Int[]],:classnames=>[""],
      :classparams=>[Int[]],:orders=>[1],:classes=>[1],:powermaps=>Any[]))
    end
    if any(isnothing, tmp) return nothing end
    if length(tmp)==1 res=copy(tmp[1].prop)
    else res=Dict{Symbol, Any}() 
      res[:classtext]=map(x->reduce(vcat,x),tcartfields(tmp,:classtext))
      res[:classnames]=map(x->join(x,","),tcartfields(tmp,:classnames))
    end
    if all(haskey.(tmp,:classparams))
      res[:classparams]=cartfields(tmp,:classparams)
    end
    res=ClassInfo(res)
    if length(tmp)>1
      if all(haskey.(tmp,:orders)) prodpowmaps(tmp,res) end
      if all(haskey.(tmp,:classes))
        res.classes=map(prod, tcartfields(tmp,:classes))
      end
    end
    res
  end::ClassInfo
end

using Primes: primes
function prodpowmaps(rec::Vector,res)
  if length(rec)==1
    if haskey(rec[1],:orders) res.orders=rec[1].orders end
    if haskey(rec[1],:powermaps) res.powermaps=rec[1].powermaps end
    return 
  end
  if all(haskey.(rec,:orders))
    res.orders=map(x->lcm(x...),tcartfields(rec,:orders))
  end
  if all(haskey.(rec,:powermaps))
    classes=cartesian(eachindex.(getproperty.(rec,:orders))...)
    maxo=maximum(res[:orders])
    l=primes(maxo)
    pm=Vector{Any}(fill(nothing,maxo))
    for p in l 
      powmap=fill(0,length(classes))
      for (i,c) in enumerate(classes)
        powmap[i]=findfirst(==(map(eachindex(c))do j
          power(rec[j][:powermaps],c[j],mod1(p,rec[j][:orders][c[j]]))
        end),classes)
      end
      pm[p]=powmap
    end
    res.powermaps=pm
  end
end

@forward Weyl.FC.G classinfo, charinfo

function Base.show(io::IO, ::MIME"text/html", ci::ClassInfo)
  print(io,latexstring(xrepr(MIME("text/plain"),ci,TeX=true)))
end

function Base.show(io::IO,t::ClassInfo)
  print(io,"ClassInfo(",t.prop,")")
end

function limitto(s,sz)
  if length(s)<sz || sz<=1 return s end
  s[1:prevind(s,sz-1)]*"…"
end

function showperiodic(io::IO,v::Vector{Int})
  if isempty(v) return "." end
  i=prev=1
  res=""
  while i<=length(v)
    for p in 2:div(length(v)-i+1,2)
#     @show i,p
      if (@view v[i:i+p-1])!=(@view v[i+p:i+2p-1]) continue end
      c=2
      while i+(c+1)*p-1<=length(v) && (@view v[i:i+p-1])==(@view v[i+c*p:i+(c+1)*p-1]) c+=1 end
      if c==2 && p==2 continue end
      res*=string(joindigits(@view v[prev:i-1]),"(",joindigits(@view v[i:i+p-1]),")",stringexp(io,c))
      i+=c*p-1;prev=i+1
#     @show i,c,p,res,prev
      break
    end
    i+=1
  end
  res*=joindigits(@view v[prev:end])
end

function Base.show(io::IO, ::MIME"text/plain", ci::ClassInfo)
  if !hasdecor(io)
    print(io,"ClassInfo(",ci.prop,")")
    return
  end
  n=length(ci.classnames)
  sz=length(string(n))+1
  t=Any[fromTeX.(Ref(io),ci.classnames)];cl=String["name"]
  sz+=max(maximum(length.(t[end])),4)+1
  for (key,tkey) in [(:classes,:length),(:orders,:order)]
    if haskey(ci,key) && ci[key]!=false
      push!(t,string.(ci[key]))
      push!(cl,string(tkey))
      sz+=max(maximum(length.(t[end])),length(cl[end]))+1
    end
  end
  if haskey(ci,:refleigen)
    m=toM(ci.refleigen)
    for i in axes(m,2)
      push!(t,xrepr.(io,m[:,i]))
      push!(cl,i==size(m,2) ? "eig" : "")
      sz+=max(maximum(length.(t[end])),length(cl[end]))+1
    end
  end
  push!(t,limitto.(showperiodic.(io,ci.classtext),displaysize(io)[2]-sz-2))
  push!(cl,"word")
  showtable(io,permutedims(string.(toM(t)));row_labels=string.(1:n),col_labels=cl,rows_label="n0")
end

const Hastype=Union{PermRootGroup,Spets,CoxSym,CoxHyp}

function Groups.conjugacy_classes(W::TW)where TW<:Hastype
  get!(W,:classes) do
    c=classinfo(W)
    map(1:length(c.classtext))do i
      C=ConjugacyClass(W,W(c.classtext[i]...),
        Dict{Symbol,Any}(:length=>c.classes[i],
                         :word=>c.classtext[i],
                         :name=>c.classnames[i],
                         :param=>c.classparams[i]))
      if haskey(c,:orders) C.order=c.orders[i] end
      if haskey(c,:malle) C.malle=c.malle[i] end
      if haskey(c,:centralizers) C.centralizers=c.centralizers[i] end
      if haskey(c,:indexclasses) C.index=c.indexclasses[i] end
      C
    end
  end
end

function Base.show(io::IO,C::ConjugacyClass{T,TW})where{T,TW<:Hastype}
  if get(io,:limit,false) || get(io,:TeX,false)
    print(io,"conjugacy_class(",C.G,",",fromTeX(io,C.name),")")
  else
    print(io,"conjugacy_class(",C.G,",",C.representative,")")
  end
end

function Groups.position_class(W::Hastype,w;verbose=false)
  l=PermGroups.positions_class(W,w)
  if length(l)==1 return only(l)
  elseif verbose println("ambiguity: ",l) end
# p=eigmat(reflrep(W,w))
# l=filter(x->eigen(conjugacy_classes(W)[x])==p,l)
# if length(l)==1 return only(l) end
  p=findfirst(c->w in c,conjugacy_classes(W)[l])
  if isnothing(p)
    @show w
    error()
  end
  l[p]
end

function Perms.reflength(C::ConjugacyClass{T,TW})where{T,TW<:Hastype}
  getp(eigen,C,:reflength)
end

Base.length(C::ConjugacyClass{T,TW}) where{T,TW<:Hastype}=C.length

Groups.word(C::ConjugacyClass{T,TW}) where{T,TW<:Hastype}=C.word

"""
`eigen(C::ConjugacyClass)` the eigenvalues on the reflection representation
of the class `C` of a reflection group or spets.
"""
function LinearAlgebra.eigen(C::ConjugacyClass{T,TW})where{T,TW<:Hastype}
  get!(C,:eigen)do
    t=refltype(C.G)
    if !any(x->haskey(x,:orbit) && (length(x.orbit)>1 || order(x.twist)>1 ||
       (haskey(x,:scalar) && !all(isone,x.scalar))),t)
      l=refleigen(C.G)
      for (i,C1) in enumerate(conjugacy_classes(C.G))
        C1.eigen=l[i]
        C1.reflength=count(!isone,C1.eigen)
      end
    else
      C.eigen=eigmat(reflrep(C.G,C.representative)) # eigmat is sorted
      C.reflength=count(!isone,C.eigen)
    end
    C.eigen
  end::Vector{Root1}
end

#--------------- CharTables -----------------------------------------
@GapObj struct CharTable{T}
  irr::Matrix{T}
  charnames::Vector{String}
  classnames::Vector{String}
  centralizers::Vector{Int}
  order::Int
end

@doc """
 CharTable is a structure to hold character tables of groups and Hecke
 algebras
""" CharTable

function Base.show(io::IO, ::MIME"text/html", ct::CharTable)
  print(io,latexstring(xrepr(MIME("text/plain"),ct,TeX=true)))
end

function Base.show(io::IO,t::CharTable)
  if hasdecor(io) || !haskey(t,:repr)
    printTeX(io,"CharTable(\$",haskey(t,:name) ? t.name : "?","\$)")
  else  print(io,t.repr)
  end
end

function Base.show(io::IO, ::MIME"text/plain", ct::CharTable)
  println(io,ct)
  showtable(io,ct.irr,row_labels=ct.charnames,col_labels=ct.classnames,dotzero=true)
end

function CharTable(t::TypeIrred;opt...)
  ct=chevieget(t,:CharTable)
  irr=collect(ct[:irreducibles]) # to make sure it is a abstractmatrix
  if !(irr isa Matrix) error("should not happen for $t") end
# irr=improve_type(irr)
  dict=Dict{Symbol,Any}(:name=>xrepr(t;TeX=true))
  for k in keys(ct)
    if !(k in (:irreducibles,:size,:name,:order,:identifier))
      dict[k]=ct[k]
    end
  end
  CharTable(irr,charnames(t;opt...,TeX=true),classnames(t;opt...,TeX=true),
            improve_type(ct[:centralizers]),ct[:size],dict)
end

function Base.prod(ctt::Vector{<:CharTable})
  if isempty(ctt)
    return CharTable(hcat(1),["Id"],["."],[1],1,Dict{Symbol,Any}(:name=>"."))
  elseif length(ctt)==1 return ctt[1]
  end
  charnames=join.(tcartfields(ctt,:charnames),",")
  classnames=join.(tcartfields(ctt,:classnames),",")
  centralizers=prod.(tcartfields(ctt,:centralizers))
  order=prod(getfield.(ctt,:order))
  irr=kron(getfield.(ctt,:irr)...)
  res=CharTable(irr,charnames,classnames,centralizers,order,Dict{Symbol,Any}(:name=>"x"))
  prodpowmaps(ctt,res)
  res
end

"""
`CharTable(WF::Spets)`

This  function returns the character table of the reflection coset `WF`. We
call *characters* of the coset `WF=W.ϕ` of the group `W` the restriction to
`W.ϕ`  of a set containing one extension of each `ϕ`-invariant character of
W  to the semidirect  product of W  with the cyclic  group generated by `ϕ`
(for  Coxeter  cosets  we  choose,  following  Lusztig,  in  each  case one
extension, called the preferred extension.)

```julia-repl
julia> W=spets(coxgroup(:D,4),Perm(1,2,4))
³D₄

julia> CharTable(W)
CharTable(³D₄)
┌─────┬─────────────────────────────────┐
│     │C₃ Ã₂ C₃+A₁ Ã₂+A₁ F₄ Ã₂+A₂ F₄(a₁)│
├─────┼─────────────────────────────────┤
│.4   │ 1  1     1     1  1     1      1│
│.1111│-1  1     1    -1  1     1      1│
│.22  │ .  2     2     . -1    -1     -1│
│11.2 │ .  .     .     . -1     3      3│
│1.3  │ 1  1    -1    -1  .    -2      2│
│1.111│-1  1    -1     1  .    -2      2│
│1.21 │ .  2    -2     .  .     2     -2│
└─────┴─────────────────────────────────┘
```
"""
function CharTable(W::Union{Hastype,FiniteCoxeterGroup};opt...)
  get!(W,:chartable)do
    t=refltype(W)
    ct=isempty(t) ? prod(CharTable[]) : prod(CharTable.(t;opt...))
    ct.name=xrepr(W;TeX=true)
    ct.repr=string("CharTable(",W,")")
    ct
  end::CharTable
end

function CharTable(W::Group)
  get!(W,:chartable)do
    error("to have CharTable(s) for arbitrary PermGroup(s) do 'using GAP'")
  end
end

function classes(ct::CharTable)
  get!(ct,:classes)do
    div.(ct.order,ct.centralizers)
  end::Vector{Int}
end

function scalarproduct(ct::CharTable,c1::AbstractVector,c2::AbstractVector;exact=true)
  v=c2'*(c1.*classes(ct))
  exact ? exactdiv.(v,ct.order) : improve_type(v//ct.order)
end

"""
`decompose(ct::CharTable,c::Vector;exact=true)`

decompose  class function `c` (given by its values on conjugacy classes) on
irreducible  characters as  given by  `CharTable` `ct`.  By default  `c` is
expected to be a virtual character so the result will be an integer vector.
If  `c` is not a virtual character  give the keyword `exact=false` to get a
correct result.
"""
function decompose(ct::CharTable,c::AbstractVector;exact=true)
  v=conj(ct.irr)*Diagonal(classes(ct))*c
  exact ? Integer.(exactdiv.(v,ct.order)) : improve_type(v//ct.order)
end

function  power(powermap,class,order)
  for (p,i) in factor(order)
    if p>length(powermap) error("missing $p-powermap")
    else pm=powermap[p]
      if isnothing(pm) error("missing $p-powermap") end
      for j in 1:i class=pm[class] end
    end
  end
  class
end
	
# returns the restriction of char to the cyclic group generated by an
# element of class (supposed of order n)
function char_powers_class(powermap::AbstractVector,char::AbstractVector,
  class::Integer,n::Integer)
  if n==1 return [char[class]] end
  MISSING=1+char[1]
  powers=fill(MISSING,n)
  powers[n]=char[1]
  for i in 1:n-1
    if powers[i]!=MISSING continue end
    # necessarily `i` divides `n`, since `i/Gcd(n,i)` is invertible
    # mod `n`, and thus powering with `i` is Galois conjugate to
    # powering with `Gcd(n,i)`
    powers[i]=char[power(powermap,class,i)]
    for j in prime_residues( div(n,i) )
      powers[mod(i*j,n)]=galois(powers[i],j) # Note: mod(i*j,n) cannot be 0.
    end
  end
  powers
end

"""
`schur_functor(ct::CharTable,φ::AbstractVector,λ)`

`schur_functor(ct::CharTable,i::Integer,λ)`

`ct`  should be a character table for some group `G`. In the first form `φ`
should  be a vector representing a character of `G`, and in the second form
`i`  denotes the `i`-th  character `φ` of  `G`. `λ` is  a partition of some
integer  `n`, representing an irreducible character  of ``𝔖 ₙ``. The function
returns  the Schur functor `φ[λ]` of φ defined by λ . This is the character
defined by the formula

``φ[λ](g)=∑_μ λ(μ)∏_{k∈1:n}φ(gᵏ)^{nₖ(μ)}``

where  `μ` runs over the partitions of  `n`, representing the classes of ``𝔖
ₙ``,  where `λ(μ)` is the  value of the character  `λ` at the class `μ`, and
`nₖ(μ)` is the number of parts of `μ` equal to `k`.

Particular  examples are the exterior  powers, corresponding to `λ=1ⁿ`, and
the symmetric powers, corresponding `λ=n`.

`ct`  should  have  `.powermaps`  containing  powermaps  for  at  least the
prime divisors of the order of `G`.
"""
function schur_functor(ct::CharTable,ch::AbstractVector,λ::AbstractVector)
  n=sum(λ)
  if any(i->λ[i-1]<λ[i],2:length(λ)) error("λ shoud be decreasing") end
  Sn=coxsym(n)
  pp=partitions(n)
  Csn=CharTable(Sn)
  chλ=Csn.irr[findfirst(==(λ),pp),:]
  map(eachindex(ch))do i
    pow=char_powers_class(ct.powermaps,ch,i,ct.orders[i])
    sum(enumerate(pp))do (j,μ)
      chλ[j]*prod(1:n)do k
        pow[mod1(k,ct.orders[i])]^count(==(k),μ)
      end//Csn.centralizers[j]
    end
  end
end
      
schur_functor(ct::CharTable,i::Integer,λ)=schur_functor(ct,ct.irr[i,:],λ)

"""
`eigen(ct::CharTable,char::AbstractVector,class)``

`eigen(ct::CharTable,char::Integer,class)``

In   the  first   form,  the   eigenvalues  of   `class`-th  class  on  the
representation affording the character whose values are given by the vector
`char`.

In   the  second   form,  the   eigenvalues  of   `class`-th  class  on  the
representation affording the `char`-th character of `ct`.

The character table `ct` must have fields `.order` containing the orders of
elements  in a conjugacy class and `.powermaps` containing powermaps for at
least the divisors of the order of the group.
"""
function LinearAlgebra.eigen(ct::CharTable,char::AbstractVector,class::Integer)
  n=ct.orders[class]
  powers=char_powers_class(ct.powermaps,char,class,n)
# compute the scalar products of the characters given by `E(n)->E(n)^i`
# with the restriction of char to the cyclic group generated by class
  eigen=Root1[]
  for i in 1:n
    append!(eigen,fill(E(n,i),div(Int(sum(j->E(n,-i*j)*powers[j],1:n)),n)))
  end
  sort(eigen,by=x->x.r)
end

LinearAlgebra.eigen(ct::CharTable,char::Integer,class::Integer)=
  eigen(ct,ct.irr[char,:],class)

"""
`on_chars(G,aut)`

`aut`  is an automorphism of  the group `G` (for  a permutation group, this
could  be  given  as  a  permutation  normalizing  `G`).  The result is the
permutation of the indices of the irreducible characters induced by `aut`.
```julia-repl
julia> WF=rootdatum("3D4")
³D₄

julia> on_chars(Group(WF),WF.phi)
(1,2,7)(8,9,12)
```
"""
function on_chars(W,aut)
  ct=CharTable(W).irr
  inv(Perm(ct,invpermute(ct,on_classes(W, aut),dims=2),dims=1))
end

"""
`representation(W,i)`

returns,   for  the  `i`-th  irreducible   representation  of  the  complex
reflection  group or Spets `W`, a list of matrices images of the generating
reflections  of `W` in a model of the representation (for Spets, the result
is  a `NamedTuple` with fields `gens`,  a representation of `Group(W)`, and
`F`,  the matrix for `W.phi` in the representation). This function is based
on  the  classification,  and  is  not  yet fully implemented for `G₃₄`; 60
representations   are  missing  out   of  169,  that   is,  about  half  of
representations  of  dimension  ≥140.  `nothing`  is returned for a missing
representation.

```julia-repl
julia> representation(complex_reflection_group(24),3)
3-element Vector{Matrix{Cyc{Int64}}}:
 [1 0 0; -1 -1 0; -1 0 -1]
 [-1 0 -1; 0 -1 (1-√-7)/2; 0 0 1]
 [-1 -1 0; 0 1 0; 0 (1+√-7)/2 -1]
```
"""
function representation(W::Union{Hastype,FiniteCoxeterGroup},i::Integer)
  tt=refltype(W)
  if isempty(tt) return Matrix{Int}[] end
  if length(tt)==1 return chevieget(tt[1],:Representation,i) end
  v=lin2cart(chevieget.(tt,:nconjugacy_classes),i)
  mm=chevieget.(tt,:Representation,v)
  if any(isnothing,mm) || any(==(false),mm) return nothing end
  if W isa Spets
    FF=map(x->x.F,mm)
    F=length(FF)==1 ? FF[1] : kron(FF...)
    mm=map(x->x.gens,mm)
  end
  n=length(tt)
  if n==1 reps=mm[1]
  else reps=vcat(map(1:n) do i
      map(mm[i]) do m
        kron(map(j->j==i ? m : one(mm[j][1]),1:n)...)
      end
    end...)
  end
  (W isa Spets) ? (gens=reps,F=F) : reps
end

"""
`representations(W)`

returns  the list  of representations  of the  complex reflection  group or
Spets `W` (see `representation`).

```julia-repl
julia> representations(coxgroup(:B,2))
5-element Vector{Vector{Matrix{Int64}}}:
 [[1;;], [-1;;]]
 [[1 0; -1 -1], [1 2; 0 -1]]
 [[-1;;], [-1;;]]
 [[1;;], [1;;]]
 [[-1;;], [1;;]]
```
"""
representations(W::Union{Hastype,FiniteCoxeterGroup})=representation.(Ref(W),1:nconjugacy_classes(W))

using SparseArrays
"""
`WGraphToRepresentation(coxrank::Integer,graph,v)`

We  store some  representations of  one-parameter Iwahori-Hecke algebras as
`W`-graphs.  For a Coxeter system `(W,S)`, a  `W`-graph is defined by a set
of  vertices `C`  with a  function `I`  which attaches  to `x∈  C` a subset
`I(x)⊂ S`, and *edge labels* which to `(x,y)∈ C^2` attach `μ(x,y)∈ K` where
`K` is the field of definition of `W`; this defines a representation of the
Hecke  algebra  with  parameters  `v`  and  `-v⁻¹`  on  a  space with basis
``{e_y}_{y∈ C}`` by:

``Tₛ(e_y)=-e_y`` if `s∈ I(y)` and otherwise
``Tₛ(e_y)=v^2 e_y+∑_{x∣s∈ I(x)} vμ(x,y)eₓ``.

The  `W`-graphs are  stored in  a compact  format to  save space.  They are
represented as a pair.
  - The  first element is a list describing `C`. Its  elements are either a
    vector   `I(x)`  of  indices  in  `eachindex(S)`,  or  an  integer  `n`
    specifying to repeat the previous element `n` more times.

  - The  second element is a list which  specifies `μ`.

We   first   describe   the   `μ`-list   for   symmetric  `W`-graphs  (when
`μ(x,y)=μ(y,x)`).  There is one  element of the  `μ`-list for each non-zero
value `m` taken by `μ`, which consists of a pair whose first element is `m`
and  whose second element is a list of  lists; if `l` is one of these lists
each  pair `[l[1],l[i]]`  represents an  edge (`x=l[1]`,`y=l[i]`) such that
`μ(x,y)=μ(y,x)=m`.  For non-symmetric `W`-graphs, the first element of each
pair  in the `μ`-list  is a pair  `[m1,m2]` and each  edge `[x,y]` obtained
from  the lists in the second element  has to be interpreted as `μ(x,y)=m1`
and `μ(y,x)=m2`.

```julia-repl
julia> W=coxgroup(:H,3)
H₃

julia> g=Wgraph(W,3)
2-element Vector{Vector{Vector{Any}}}:
 [[2], [1, 2], [1, 3], [1, 3], [2, 3]]
 [[-1, [[1, 3], [2, 4], [3, 5], [4, 5]]]]

julia> WGraphToRepresentation(3,g,Pol(:x))
3-element Vector{Matrix{Pol{Int64}}}:
 [x² 0 … 0 0; 0 -1 … 0 0; … ; 0 0 … -1 -x; 0 0 … 0 x²]
 [-1 0 … 0 0; 0 -1 … -x 0; … ; 0 0 … x² 0; 0 0 … -x -1]
 [x² 0 … 0 0; 0 x² … 0 0; … ; 0 -x … -1 0; 0 0 … 0 -1]
```
"""
function WGraphToRepresentation(rk::Integer,gr::Vector,v)
# Jean Michel june/december 2003 from  code/data of Geck, Marin, Alvis,
# Naruse, Howlett,Yin)
  I=Vector{Int}[]
  for S in gr[1]
    if S isa Integer append!(I,fill(I[end],S))
    else push!(I,S)
    end
  end
  dim=length(I)
  T=Int
  function prom(a)
    if a isa Vector
      for u in a prom(u) end
    else T=promote_type(T,typeof(a))
    end
  end
  prom(gr[2])
  v=improve_type(v)
  T=promote_type(T,typeof(v))
  S=map(i->spzeros(T,dim,dim),1:rk)
  for j in 1:dim
    for i in 1:rk
      if i in I[j] S[i][j,j]=-one(v)
      else         S[i][j,j]=v^2
      end
    end
  end
  for i in gr[2]
    if i[1] isa Vector mu=i[1] else mu=[i[1],i[1]] end
    for l in i[2]
      x=l[1]
      for y in l[2:end]
        for j in I[y] if !(j in I[x]) S[j][y,x]=mu[2]*v end  end
        for j in I[x] if !(j in I[y]) S[j][x,y]=mu[1]*v end  end
      end
    end
  end
  density=sum(nnz.(S))/(rk*dim^2)
  density>0.2 ? Array.(S) : S
end

# the next function returns the dual W-graph of gr (for an Hecke algebra of
# rank rk). A dual W-graph corresponds to a Curtis Dual representation.
function DualWGraph(rk,gr)
  [map(x->x isa Integer ? x : setdiff(1:rk,x),gr[1]),
   map(((x,y),)->x isa Vector ? [-reverse(x),y] : [-x,y],gr[2])]
end

############################################################################
# How to interpret W-graphs for complex reflection groups with one orbit of
# reflections, for hecke(W,[vars]).

function WGraph2Representation(a,vars)
# println("a=$a vars=$vars")
  nodes=a[1]
  varno=function(n,j) # a function from 1:ngens(W) to eachindex(vars)
    if n[1] isa Vector p=findfirst(x->j in x,n)
      if isnothing(p) p=length(vars) end
    else p=(j in n) ? 1 : 2
    end
    p
  end
  flat(l)=any(x->x isa Vector,l) ? flat(reduce(vcat,l)) : l
  rk=maximum(Int.(flat(nodes))) # ngens(W)
  dim=length(nodes)
  R=map(j->map(k->vars[varno(nodes[k],j)],1:dim),1:rk)
  T=eltype(flat(a[2]))
  R=Array.(Diagonal.(map(x->T.(x),R)))
# println("R=$(typeof(R))$R")
  for r in a[2]
#   println("r=$r")
    for k in [3,4]
      i1=Int(r[k-2]);i2=Int(r[5-k])
      if r[k] isa Vector
        for j in 2:2:length(r[k]) R[Int(r[k][j-1])][i1,i2]=r[k][j] end
      else
        for i in 1:rk
          if varno(nodes[i1],i)<varno(nodes[i2],i)
            R[i][i1,i2]=r[k]
          end
        end
      end
    end
  end
# println("R=$(typeof(R))$R")
  R
end

function charnames(io::IO,c::Union{CharInfo,CharTable})
  cn=c.charnames
  for k in [:spaltenstein, :frame, :malle, :kondo, :gp, :lusztig]
    if get(io,k,false) && haskey(c,k) cn=string.(c[k]) end
  end
  fromTeX.(Ref(io),cn)
end

charnames(io::IO,W::Weyl.FC)=charnames(io,W.G)

"""
`charnames(W;options...)` or `charnames(io::IO,W)`

returns  the list of character names for  the reflection group or Spets `W`
(`W`  could also  be a  `CharTable` or  `CharInfo`). The  options may imply
alternative  names in certain cases, or  a different formatting of names in
general.  They are specified by `IO` attributes  if giving an `IO` as first
argument.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> charnames(W;limit=true)
6-element Vector{String}:
 "φ₁‚₀"
 "φ₁‚₆"
 "φ′₁‚₃"
 "φ″₁‚₃"
 "φ₂‚₁"
 "φ₂‚₂"

julia> charnames(W;TeX=true)
6-element Vector{String}:
 "\\phi_{1,0}"
 "\\phi_{1,6}"
 "\\phi_{1,3}'"
 "\\phi_{1,3}''"
 "\\phi_{2,1}"
 "\\phi_{2,2}"

julia> charnames(W;spaltenstein=true,limit=true)
6-element Vector{String}:
 "1"
 "ε"
 "εₗ"
 "ε_c"
 "θ′"
 "θ″"

julia> charnames(W;spaltenstein=true,TeX=true)
6-element Vector{String}:
 "1"
 "\\varepsilon"
 "\\varepsilon_l"
 "\\varepsilon_c"
 "\\theta'"
 "\\theta''"
```
The  last two commands  show the character  names used in [spalt85](@cite).
Other  keywords  available  are  `frame,  kondo,  malle,  gp, lusztig`, see
[`charinfo`](@ref) for a description.
"""
function charnames(io::IO,W::Union{Group,Coset})
  if applicable(refltype,W) charnames(io,charinfo(W))
  else fromTeX.(Ref(io),CharTable(W).charnames)
  end
end

charnames(W;opt...)=charnames(IOContext(stdout,opt...),W)

charnames(t::TypeIrred;opt...)=charnames(IOContext(stdout,opt...),t)

function charnames(io::IO,t::TypeIrred)
  ci=chevieget(t,:charinfo)
  cn=ci[:charnames]
  for k in [:spaltenstein, :frame, :malle, :kondo, :gp, :lusztig]
    if get(io,k,false) && haskey(ci,k) cn=string.(ci[k]) end
  end
  fromTeX.(Ref(io),cn)
end

function charnames(io::IO,tt::Vector{TypeIrred})
  if isempty(tt) return ["Id"] end
  cn=map(t->charnames(io,t),tt)
  if length(cn)==1 cn[1]
  else map(x->join(x,fromTeX(io,"\\otimes ")),tcartesian(cn...))
  end
end

classnames(io::IO,ct::Union{ClassInfo,CharTable})=fromTeX.(Ref(io),ct.classnames)

"""
`classnames(W;options...)` or `classnames(io::IO,W)`

returns  the list of class names for the reflection group `W` (`W` may also
be  a `CharTable` or  the result of  `classinfo`). The optional options are
IOContext  attributes which can give alternative names in certain cases, or
a different formatting of names in general. They can be specified by giving
an IO as argument.
```julia-repl
julia> classnames(coxgroup(:G,2),TeX=true)
6-element Vector{String}:
 "A_0"
 "\\tilde A_1"
 "A_1"
 "G_2"
 "A_2"
 "A_1+\\tilde A_1"
```
"""
function classnames(io::IO,W)
  if applicable(refltype,W)
    c=classinfo(W)
    cn=c.classnames
    for k in [:malle]
      if get(io,k,false) && haskey(c,k) cn=string.(c[k]) end
    end
  else
    cn=CharTable(W).classnames
  end
  fromTeX.(Ref(io),cn)
end

classnames(W;opt...)=classnames(IOContext(stdout,opt...),W)

function classnames(t::TypeIrred;opt...)
  c=classinfo(t)
  cn=c.classnames
  for k in [:malle]
    if get(opt,k,false) && haskey(c,k) cn=string.(c[k]) end
  end
  cn
end

@GapObj struct InductionTable{T}
  scalar::Matrix{T}
  gcharnames::Vector{String}
  ucharnames::Vector{String}
  identifier::String
end

"""
   `induction_table(u,g)`

returns   an  object  describing  the   decomposition  of  the  irreducible
characters  of the subgroup  `u` induced to  the group `g`.  At the repl or
IJulia  or Pluto,  a table  is displayed  where the  rows correspond to the
characters  of the parent group, and the  columns to those of the subgroup.
The  returned  object  has  a  field  `scalar`  which  is  a  `Matrix{Int}`
containing  the  induction  table,  and  the  other fields contain labeling
information taken from the character tables of `u` and `g` when it exists.

```julia-rep1
julia> g=Group([Perm(1,2),Perm(2,3),Perm(3,4)])
Group([(1,2),(2,3),(3,4)])

julia> u=Group( [ Perm(1,2), Perm(3,4) ])
Group([(1,2),(3,4)])

julia> induction_table(u,g)  #     needs "using GAP"
Induction table from Group((1,2),(3,4)) to Group((1,2),(2,3),(3,4))
┌───┬───────────────┐
│   │X.1 X.2 X.3 X.4│
├───┼───────────────┤
│X.1│  .   1   .   .│
│X.2│  .   1   1   1│
│X.3│  1   1   .   .│
│X.4│  1   .   1   1│
│X.5│  1   .   .   .│
└───┴───────────────┘
```

```julia-repl
julia> g=coxgroup(:G,2)
G₂

julia> u=reflection_subgroup(g,[1,6])
G₂₍₁₅₎=A₂

julia> t=induction_table(u,g)
Induction table from G₂₍₁₅₎=A₂ to G₂
┌─────┬────────┐
│     │111 21 3│
├─────┼────────┤
│φ₁‚₀ │  .  . 1│
│φ₁‚₆ │  1  . .│
│φ′₁‚₃│  1  . .│
│φ″₁‚₃│  .  . 1│
│φ₂‚₁ │  .  1 .│
│φ₂‚₂ │  .  1 .│
└─────┴────────┘
```

`IO` attributes can be transmitted to the table format method

```julia-rep1
julia> xdisplay(t;rows=[5],cols=[3,2])
Induction table from G₂₍₁₅₎=A₂ to G₂
┌─────┬────┐
│     │3 21│
├─────┼────┤
│φ₂‚₁ │.  1│
└─────┴────┘
```

It is also possible to TeX induction tables with `xdisplay(t;TeX=true)`.

`induction_table` also works for spets (reflection cosets).
"""
function induction_table(u,g)
  tu=CharTable(u)
  tg=CharTable(g)
  f=fusion_conjugacy_classes(u,g)
  function myexactdiv(a,b)
    q,r=divrem(a,b)
    if !iszero(r) error("$b does not exactly divide $a") end
    q
  end
  function myexactdiv(a::Cyc,b)
    q=Cyc(conductor(a),exactdiv(a.d,b))
    if q*b!=a error("$b does not exactly divide $a") end
    q
  end
  cl=myexactdiv.(length(u),tu.centralizers)
  sc=myexactdiv.(tu.irr*Diagonal(cl)*tg.irr[:,f]',length(u))
  InductionTable(permutedims(sc),tg.charnames,tu.charnames,
  "Induction table from \$"*TeXs(u)*"\$ to \$"*TeXs(g)*"\$",
  Dict{Symbol,Any}(:repr=>string("induction_table(",u,",",g,")")))
end

function Base.show(io::IO, ::MIME"text/html", t::InductionTable)
  print(io,latexstring(xrepr(MIME("text/plain"),t,TeX=true)))
end

function Base.show(io::IO,t::InductionTable)
  if !hasdecor(io) && haskey(t,:repr) print(io,t.repr)
  else printTeX(io,t.identifier)
  end
end

function Base.show(io::IO,::MIME"text/plain",t::InductionTable)
  println(io,t)
  showtable(io,t.scalar;row_labels=t.gcharnames,col_labels=t.ucharnames,dotzero=true)
end

"""
`j_induction_table(H, W)`

computes  the decomposition  into irreducible  characters of the reflection
group  `W`  of  the  `j`-induced  of  the  irreducible  characters  of  the
reflection  subgroup  `H`.  The  `j`-induced  of  `φ`  is  the  sum  of the
irreducible  components of the induced of  `φ` which have same `b`-function
(see `charinfo`) as `φ`. What is returned is an `InductionTable` struct.

```julia-repl
julia> W=coxgroup(:D,4)
D₄

julia> H=reflection_subgroup(W,[1,3])
D₄₍₁₃₎=A₂Φ₁²

julia> j_induction_table(H,W)
j-induction table from D₄₍₁₃₎=A₂Φ₁² to D₄
┌─────┬────────┐
│     │111 21 3│
├─────┼────────┤
│11+  │  .  . .│
│11-  │  .  . .│
│1.111│  .  . .│
│.1111│  .  . .│
│11.2 │  .  . .│
│1.21 │  1  . .│
│.211 │  .  . .│
│2+   │  .  . .│
│2-   │  .  . .│
│.22  │  .  . .│
│1.3  │  .  1 .│
│.31  │  .  . .│
│.4   │  .  . 1│
└─────┴────────┘
```
"""
function j_induction_table(u,g)
  tbl=induction_table(u,g)
  bu=charinfo(u).b
  bg=charinfo(g).b
  t=copy(tbl.scalar)
  for (i,bi) in pairs(bu), (j,bj) in pairs(bg)
    if bi!=bj t[j,i]=0 end
  end
  InductionTable(t,tbl.gcharnames,tbl.ucharnames,
  "j-induction table from \$"*TeXs(u)*"\$ to \$"*TeXs(g)*"\$",
  Dict{Symbol,Any}(:repr=>string("j_induction_table(",u,",",g,")")))
end

"""
`J_induction_table(H, W)`

computes  the decomposition  into irreducible  characters of the reflection
group  `W`  of  the  `J`-induced  of  the  irreducible  characters  of  the
reflection  subgroup  `H`.  The  `J`-induced  of  `φ`  is  the  sum  of the
irreducible  components of the induced of  `φ` which have same `a`-function
(see `charinfo`) as `φ`. What is returned is an `InductionTable` struct.

```julia-repl
julia> W=coxgroup(:D,4)
D₄

julia> H=reflection_subgroup(W,[1,3])
D₄₍₁₃₎=A₂Φ₁²

julia> J_induction_table(H,W)
J-induction table from D₄₍₁₃₎=A₂Φ₁² to D₄
┌─────┬────────┐
│     │111 21 3│
├─────┼────────┤
│11+  │  .  . .│
│11-  │  .  . .│
│1.111│  .  . .│
│.1111│  .  . .│
│11.2 │  1  . .│
│1.21 │  1  . .│
│.211 │  .  . .│
│2+   │  .  . .│
│2-   │  .  . .│
│.22  │  .  . .│
│1.3  │  .  1 .│
│.31  │  .  . .│
│.4   │  .  . 1│
└─────┴────────┘
```
"""
function J_induction_table(u,g)
  tbl=induction_table(u,g)
  bu=charinfo(u).a
  bg=charinfo(g).a
  t=copy(tbl.scalar)
  for (i,bi) in pairs(bu), (j,bj) in pairs(bg)
    if bi!=bj t[j,i]=0 end
  end
  InductionTable(t,tbl.gcharnames,tbl.ucharnames,
  "J-induction table from \$"*TeXs(u)*"\$ to \$"*TeXs(g)*"\$",
  Dict{Symbol,Any}(:repr=>string("J_induction_table(",u,",",g,")")))
end

"""
`discriminant(W)`

returns  the  discriminant  of  the  complex  reflection  group  `W`,  as a
polynomial in the fundamental invariants. The discriminant is the invariant
obtained  by  taking  the  product  of  the  linear  forms  describing  the
reflecting   hyperplanes  of  `W`,   each  raised  to   the  order  of  the
corresponding  reflection. The discriminant  is returned as  a function `f`
such  that  the  discriminant  in  the  variables  `a₁,…,aₙ` is obtained by
calling `f(a₁,…,aₙ)`. For the moment, this function is implemented only for
the exceptional complex reflection groups `G₄` to `G₃₃`.

```julia-repl
julia> W=complex_reflection_group(4);@Mvp x,y

julia> discriminant(W)(x,y)
Mvp{Int64}: x³-y²
```
"""
function LaurentPolynomials.discriminant(W::Group)
  t=refltype(W)
  if isempty(t) return ()->Mvp(1)
  elseif length(t)==1 return chevieget(t[1],:Discriminant)
  else error("not implemented for non-irreducible ",W)
  end
end

function decomposition_matrix(t::TypeIrred,p)
  m=chevieget(t,:DecompositionMatrix,p)
  if m==false
    error("decomposition_matrix(",t,",",p,") not implemented")
  end
  n=chevieget(t,:nconjugacy_classes)
  append!(m,map(i->[[i],[[1]]],setdiff(1:n,union(first.(m)))))
  res=cat(map(x->toM(x[2]),m)...;dims=(1,2))
  res[sortperm(vcat(first.(m)...)),:]
end

"""
`decomposition_matrix(W,p)`

This  provides an interface to some  decomposition matrices for Weyl groups
available  in the Chevie  library: those for  `W(E₆), W(E₇), W(E₈)` for `p`
dividing `|W|`.
"""
function decomposition_matrix(W,p)
  m=decomposition_matrix.(refltype(W),p)
  map(x->prod.(cartesian(x)),cartesian(toL.(m)...))
end

end
