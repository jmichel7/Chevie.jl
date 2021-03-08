"""
The  `CharTable` of a finite complex reflection group `W` is computed using
the  decomposition of `W` in irreducible  groups (see `refltype`). For each
irreducible  group the character  table is either  computed using recursive
formulas  for the infinite series,  or read into the  system from a library
file  for the  exceptional types.  Thus, character  tables can  be obtained
quickly  even for very large groups  (e.g., E₈). Similar remarks apply for
conjugacy classes.

The  conjugacy  classes  and  irreducible  characters of irreducible finite
complex reflection groups have canonical labelings by certain combinatorial
objects;  these labelings are used in the  tables we give. For the classes,
these  are partitions or partition tuples  for the infinite series, or, for
exceptional  Coxeter  groups,  Carter's  admissible  diagrams
[Carter1972](biblio.htm#Car72); for
other  primitive  complex  reflection  groups  we  just  use  words  in the
generators  to specify  the classes.  For the  characters, these  are again
partitions  or partition tuples for the infinite series, and for the others
they  are pairs  of two  integers `(d,e)`  where `d`  is the  degree of the
character  and  `e`  is  the  smallest  symmetric  power  of the reflection
representation  containing  the  given  character  as  a  constituent  (the
`b`-invariant  of the character). This information is obtained by using the
functions `classinfo` and `charinfo`. When you display the character table,
the canonical labelings for classes and characters are those displayed.

A  typical example  is `coxgroup(:A,n)`,  the symmetric  group `𝔖ₙ₊₁` where
classes and characters are parameterized by partitions of `n+1`.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> CharTable(W)
CharTable(A₃)
    │1111 211 22 31  4
────┼──────────────────
1111│   1  -1  1  1 -1
211 │   3  -1 -1  .  1
22  │   2   .  2 -1  .
31  │   3   1 -1  . -1
4   │   1   1  1  1  1

julia> W=coxgroup(:G,2)
G₂

julia> ct=CharTable(W)
CharTable(G₂)
     │A₀ Ã₁ A₁ G₂ A₂ A₁+Ã₁
─────┼─────────────────────
φ₁‚₀ │ 1  1  1  1  1     1
φ₁‚₆ │ 1 -1 -1  1  1     1
φ′₁‚₃│ 1  1 -1 -1  1    -1
φ″₁‚₃│ 1 -1  1 -1  1    -1
φ₂‚₁ │ 2  .  .  1 -1    -2
φ₂‚₂ │ 2  .  . -1 -1     2

julia> ct.charnames
6-element Vector{String}:
 "\\phi_{1,0}"
 "\\phi_{1,6}"
 "\\phi_{1,3}'"
 "\\phi_{1,3}''"
 "\\phi_{2,1}"
 "\\phi_{2,2}"

julia> ct.classnames
6-element Vector{String}:
 "A_0"
 "\\tilde A_1"
 "A_1"
 "G_2"
 "A_2"
 "A_1+\\tilde A_1"
```

Recall  that our groups acts a reflection group on the vector space `V`, so
have  fake degrees  (see "fakeDegree").  The valuation  and degree of these
give  two  integers  `b,B`  for  each  irreducible  character  of  `W` (see
`charinf(W)[:b]`  and  `charinfo(W)[:B]`).  For  finite Coxeter groups, the
valuation  and degree of  the generic degrees  of the one-parameter generic
Hecke  algebra  give  two  more  integers  `a,A` (see `charinfo(W)[:a]` and
`charinfo(W)[:A]`,   and  [Carter1985,  Ch.11](biblio.htm#Car85]  for  more
details). These will also be used in the operations of truncated inductions
explained in the chapter "Reflection subgroups".

Iwahori-Hecke  algebras and  cyclotomic Hecke  algebras also have character
tables, see the corresponding chapters.

We  now describe for each type our conventions for labeling the classes and
characters.

Type  `Aₙ` (`n≥0`). In this  case we have  `W ≅ 𝔖ₙ₊₁`. The classes and
characters  are labeled by partitions of `n+1`. The partition corresponding
to  a class describes  the cycle type  for the elements  in that class; the
representative   in  '.classtext'   is  the   concatenation  of  the  words
corresponding  to each part, and to a part `i` is associated the product of
`i-1`  consecutive generators (starting one  higher that the last generator
used  for the previous  parts). The partition  corresponding to a character
describes  the type of  the Young subgroup  such that the trivial character
induced  from this  subgroup contains  that character with multiplicity `1`
and such that every other character occurring in this induced character has
a  higher `a`-value. Thus, the sign  character corresponds to the partition
`(1ⁿ⁺¹)`  and  the  trivial  character  to  the  partition  `(n+1)`. The
character of the reflection representation of `W` is labeled by `(n,1)`.

Type  `Bₙ`  (`n≥2`).  In  this  case  `W=W(Bₙ)` is isomorphic to the wreath
product  of the cyclic  group of order  `2` with the  symmetric group `𝔖ₙ`.
Hence  the classes and characters are  parameterized by pairs of partitions
such  that the total sum of their  parts equals `n`. The pair corresponding
to  a class describes the signed cycle type for the elements in that class,
as in [Carter1972](biblio.htm#Car72). We use the convention that if `(λ,μ)`
is such a pair then `λ` corresponds to the positive and `μ` to the negative
cycles.  Thus, `(1ⁿ,-)` and `(-,1ⁿ)` label  the trivial class and the class
containing  the longest element, respectively. The pair corresponding to an
irreducible character is determined via Clifford theory, as follows.

We  have a semidirect product decomposition `W(Bₙ)=N ⋊ 𝔖ₙ` where `N` is the
standard  `n`-dimensional  `𝔽₂ⁿ`-vector  space.  For  `a,b  ≥  0` such that
`n=a+b` let `η_{a,b}` be the irreducible character of `N` which takes value
`1`  on the first `a` standard basis vectors and value `-1` on the next `b`
standard  basis vectors of `N`. Then  the inertia subgroup of `η_{a,b}` has
the  form `T_{a,b}=N.(𝔖_a × 𝔖_b)` and  we can extend `η_{a,b}` trivially to
an  irreducible  character  `η̃_{a,b}`  of  `T_{a,b}`.  Let  `α` and `β` be
partitions  of `a` and `b`, respectively. We take the tensor product of the
corresponding  irreducible characters of `𝔖_a` and `𝔖_b` and regard this as
an  irreducible  character  of  `T_{a,b}`.  Multiplying this character with
`η̃_{a,b}`  and  inducing  to  `W(Bₙ)`  yields an irreducible character `χ=
χ_{(α,β)}`  of `W(Bₙ)`. This defines the correspondence between irreducible
characters and pairs of partitions as above.

For example, the pair `((n),-)` labels the trivial character and `(-,(1ⁿ))`
labels  the  sign  character.  The  character  of  the  natural  reflection
representation is labeled by `((n-1),(1))`.

Type  `Dₙ` (`n≥4`). In this case `W=W(Dₙ)` can be embedded as a subgroup of
index  `2` into the Coxeter  group `W(Bₙ)`. The intersection  of a class of
`W(Bₙ)` with `W(Dₙ)` is either empty or a single class in `W(Dₙ)` or splits
up  into two classes in  `W(Dₙ)`. This also leads  to a parameterization of
the  classes of `W(Dₙ)` by pairs of  partitions `(λ,μ)` as before but where
the  number of parts of `μ` is even and where there are two classes of this
type  if `μ` is empty and all parts of  `λ` are even. In the latter case we
denote  the two classes in `W(Dₙ)` by `(λ,+)` and `(λ,-)`, where we use the
convention  that  the  class  labeled  by `(λ,+)` contains a representative
which  can be written  as a word  in `{s₁,s₃,…,sₙ}` and  `(λ,-)` contains a
representative which can be written as a word in `{s₂,s₃, …,sₙ}`.

By  Clifford theory the restriction of  an irreducible character of `W(Bₙ)`
to  `W(Dₙ)`  is  either  irreducible  or  splits  up  into  two irreducible
components.  Let `(α,β)` be  a pair of  partitions with total  sum of parts
equal to `n`. If `α!=β` then the restrictions of the irreducible characters
of  `W(Bₙ)` labeled  by `(α,β)`  and `(β,α)`  are irreducible and equal. If
`α=β`  then the restriction of the character labeled by `(α,α)` splits into
two  irreducible components  which we  denote by  `(α,+)` and `(α,-)`. Note
that  this can only happen if `n` is  even. In order to fix the notation we
use  a  result  of  [Stembridge1989](biblio.htm#Ste89)  which describes the
value  of the  difference of  these two  characters on  a class of the form
`(λ,+)`  in terms of the character values of the symmetric group `𝔖_{n/2}`.
Recall  that it is implicit  in the notation `(λ,+)`  that all parts of `λ`
are even. Let `λ'` be the partition of `n/2` obtained by dividing each part
by  `2`. Then the value of `χ_{(α,-)}-χ_{(α,+)}` on an element in the class
`(λ,+)` is given by `2^{k(λ)}` times the value of the irreducible character
of  `𝔖_{n/2}` labeled by `α` on the class of cycle type `λ'`. (Here, `k(λ)`
denotes the number of non-zero parts of `λ`.)

The  labels for the trivial, the  sign and the natural reflection character
are the same as for `W(Bₙ)`, since these characters are restrictions of the
corresponding characters of `W(Bₙ)`.

The groups `G(d,1,n)`.
They  are isomorphic to the wreath product of the cyclic group of order `d`
with  the  symmetric  group  `𝔖ₙ`.  Hence  the  classes  and characters are
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

The groups `G(de,e,n)`.
They  are normal  subgroups of  index `e`  in `G(de,1,n)`.  The quotient is
cyclic,  generated by the image `g`  of the first generator of `G(de,1,n)`.
The  classes are parameterized as the  classes of `G(de,e,n)` with an extra
information for a component of a class which splits.

According  to  [Hugues1985](biblio.htm#Hu85),  a  class  `C` of `G(de,1,n)`
parameterized  by a `de`-partition ``(S₀,…,S_{de-1})`` is in `G(de,e,n)` if
`e`  divides ``∑ᵢ i ∑_{p∈ Sᵢ}p``. It  splits in `d` classes for the largest
`d`  dividing `e` and all parts of all  `Sᵢ` and such that `Sᵢ` is empty if
`d`  does not divide `i`. If `w` is in `C` then 'gⁱ w g⁻ⁱ' for 'i in 0:d-1'
are  representatives of the classes of `G(de,e,n)` which meet `C`. They are
described by appending the integer `i` to the label for `C`.

The  characters are described by Clifford theory. We make `g` act on labels
for  characters of `G(de,1,n)`  . The action  of `g` permutes circularly by
`d`  the partitions in the `de`-tuple.  A character has same restriction to
`G(de,e,n)`  as its transform by `g`.  The number of irreducible components
of its restriction is equal to the order `k` of its stabilizer under powers
of  `g`.  We  encode  a  character  of  `G(de,e,n)`  by first, choosing the
smallest  for lexicographical order label  of a character whose restriction
contains  it; then this label is periodic with a motive repeated `k` times;
we  represent the  character by  one of  these motives,  to which we append
`E(k)ⁱ` for 'i in 0:k-1' to describe which component of the restriction we
choose.

Types `G₂` and `F₄`. The matrices of character values and the orderings and
labelings  of  the  irreducible  characters  are  exactly  the  same  as in
[Carter1985,  p.412/413](biblio.htm#Car85):  in  type  `G₂`  the  character
`φ₁,₃'`  takes the value -1 on the reflection associated to the long simple
root;  in type `F₄`, the characters `φ₁,₁₂'`, `φ₂,₄'`, `φ₄,₇'`, `φ₈,₉'` and
`φ₉,₆'` occur in the induced of the identity from the `A₂` corresponding to
the  short  simple  roots;  the  pairs  (`φ₂,₁₆'`,  `φ₂,₄″`)  and (`φ₈,₃'`,
`φ₈,₉″`)  are  related  by  tensoring  by  sign; and finally `φ₆,₆″` is the
exterior  square of the  reflection representation. Note,  however, that we
put  the long root at  the left of the  Dynkin diagrams to be in accordance
with the conventions in [Lusztig1985, (4.8) and (4.10)](biblio.htm#Lus85).

The classes are labeled by Carter's admissible diagrams
[Carter1972](biblio.htm#Car72).  A character  is labeled  by a pair `(d,b)`
where  `d` denotes the  degree and `b`  the corresponding `b`-invariant. If
there  are several characters with the same  pair `(d,b)` we attach a prime
to them, as in [Carter1985](biblio.htm#Car85).

Types  `E₆,E₇,E₈`. The character  tables are obtained  by specialization of
those  of the Hecke algebra. The classes are labeled by Carter's admissible
diagrams [Carter1972](biblio.htm#Car72). A character is labeled by the pair
`(d,b)`  where  `d`  denotes  the  degree  and  `b`  is  the  corresponding
`b`-invariant.  For  these  types,  this  gives  a  unique  labeling of the
characters.

Non-crystallographic  types `I₂(m)`, `H₃`, `H₄`. In these cases we do not
have  canonical  labelings  for  the  classes.  We  label  them  by reduced
expressions.

Each  character for  type `H₃`  is uniquely  determined by the pair `(d,b)`
where  `d` is the degree and  `b` the corresponding `b`-invariant. For type
`H₄`  there are just  two characters (those  of degree `30`)  for which the
corresponding  pairs are  the same.  These two  characters are nevertheless
distinguished  by  their  fake  degrees:  the  character `φ₃₀,₁₀'` has fake
degree  `q¹⁰+q¹²+` higher terms, while `φ₃₀,₁₀″` has fake degree `q¹⁰+q¹⁴+`
higher  terms. The characters in the table for type `H₄` are ordered in the
same way as in [Alvis and Lusztig1982](biblio.htm#AL82).

Finally,  the characters  of degree `2`  for type  `I₂(m)` are  ordered as
follows.  The matrix representations affording the characters of degree `2`
are given by:
`` ρ_j : s₁s₂ ↦
\\begin{pmatrix}E(m)^j&0\\\\0&E(m)^{-j}\\end{pmatrix},
 s₁↦\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix},``
where  `1 ≤ j ≤  ⌊(m-1)/2⌋`. The reflection representation
is  `ρ₁`. The  characters in  the table  are ordered by listing
first the characters of degree 1 and then `ρ₁,ρ₂,…`.

Primitive complex reflection groups `G₄` to `G₃₄`.
The  groups `G₂₃=H₃`, `G₂₈=F₄`, `G₃₀=H₄` are exceptional Coxeter groups and
have  been  explained  above.  Similarly  for  the  other groups labels for
characters  consist primarily  of the  pair `(d,b)`  where `d`  denotes the
degree  and `b` is the corresponding  `b`-invariant. This is sufficient for
`G₄`,  `G₁₂`, `G₂₂` and `G₂₄`. For other  groups there are pairs or triples
of  characters which  have the  same `(d,b)`  value. We  disambiguate these
according  to  the  conventions  of [Malle2000](biblio.htm#Mal00) for `G₂₇,
G₂₉, G₃₁, G₃₃` and `G₃₄`:

-  For `G₂₇`:
The  fake degree  of `φ₃,₅'`  (resp. `φ₃,₂₀'`,  `φ₈,₉″`) has smaller degree
that  of  `φ₃,₅″`  (resp.  `φ₃,₂₀″`,  `φ₈,₉'`). The characters `φ₅,₁₅'` and
`φ₅,₆'` occur with multiplicity 1 in the induced from the trivial character
of  the parabolic subgroup  of type `A₂`  generated by the  first and third
generator  (it is asserted mistakenly in [Malle2000](biblio.htm#Mal00) that
`φ₅,₆″` does not occur in this induced; it occurs with multiplicity 2).

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
The  terms of  maximal degree  of the  fakedegree of `φ₁₀,₈'` are `q²⁸+q²⁶`
while  for `φ₁₀,₈'` they are `q²⁸+q²⁴`. The  terms of maximal degree of the
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
they are compatible with those of [MalleRouquier2003](biblio.htm#MR03) and
[BroueMalleMichel2014](biblio.htm#BMM14).

-  For `G₅`:
We  let `W=ComplexReflectionGroup(5)`,  so the  generators are  `W(1)` and
`W(2)`.

The  character `φ₁,₄'` (resp. `φ₁,₁₂'`, `φ₂,₃'`) takes the value `1` (resp.
`ζ₃`,  `-ζ₃`)  on  `W(1)`.  The  character  `φ₁,₈″` is complex conjugate to
`φ₁,₁₆`,  and the character  `φ₁,₈'` is complex  conjugate to `φ₁,₄'` . The
character  `φ₂,₅″` is complex  conjugate to `φ₂,₁`;  `φ₂,₅'` take the value
`-1` on `W(1)`. The character `φ₂,₇'` is complex conjugate to `φ₂,₅'`.

-  For `G₇`:
We  let `W=ComplexReflectionGroup(7)`,  so the  generators are
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

-  We define now a new invariant  for characters: consider all the pairs of
   irreducible   characters  `χ`  and  `ψ`  uniquely  determined  by  their
   `(d,b)`-invariant such that `φ` occurs with non-zero multiplicity `m` in
   `χ⊗ψ`.  We define  `t(φ)` to  be the  minimal (for  lexicographic order)
   possible list `(d(χ),b(χ),d(ψ),b(ψ),m)`.

For  the not  yet labeled  pairs, if  the t-invariants are different, label
`φ_{d,b}'` the character with the minimal `t`-invariant.

After  applying  the  last  rule  all  the  pairs  will be labelled for the
considered  groups. The labelling obtained  is compatible for `G₂₅`, `G₂₆`,
`G₃₂`  and `G₃₃`  with that  of [Malle2000](biblio.htm#Mal00)  and for `G₈`
with that described in [MalleRouquier2003](biblio.htm#MR03).

We  should  emphasize  that  for  all  groups  with  a  few exceptions, the
parameters  for characters do  not depend on  any non-canonical choice. The
exceptions  are `G(de,e,n)` with `e>1`, and `G₅`, `G₇`, `G₂₇`, `G₂₈`, `G₂₉`
and  `G₃₄`, groups  which admit  outer automorphisms  preserving the set of
reflections,  so choices  of a  particular value  on a particular generator
must be made for characters which are not invariant by these automorphisms.

Labels  for the classes. For the exceptional complex reflection groups, the
labels  for the classes represent the  decomposition of a representative of
the  class as a product of generators, with the additional conventions that
'z'  represents the generator  of the center  and for well-generated groups
'c'  represents a Coxeter element  (a product of the  generators which is a
regular element for the highest reflection degree).
"""
module Chars

using ..Gapjm

export charinfo, classinfo, fakedegree, fakedegrees, CharTable, representation,
  WGraphToRepresentation, DualWGraph, WGraph2Representation, charnames,
  representations, InductionTable, classes, jInductionTable, JInductionTable,
  decompose, on_chars, detPerm

"""
`fakedegree(W, φ, q)`

returns  the  fake degree  of  the  character  of parameter  φ  (see
:CharParams)  of  the  reflection  group `W`,  evaluated  at  `q`  (see
"fakeDegrees" for a definition of the fake degrees).

```julia-repl
julia> Chars.fakedegree(coxgroup(:A,2),[[2,1]],Pol(:q))
Pol{Cyc{Int64}}: q²+q
```
"""
function fakedegree(W,p,q)
  typ=refltype(W)
  if isempty(typ) return one(q) end
  prod(map((t,p)->fakedegree(t,p,q),typ,p))
end

function fakedegree(t::TypeIrred,p,q)
  if haskey(t,:scalar) q=prod(s->q*conj(s),t.scalar)
  elseif haskey(t,:orbit) q=q^length(t.orbit)
  end
  getchev(t,:FakeDegree,p,q)
end

"""
`fakedegrees(W , q)`

returns  a list holding the fake degrees of the reflection group `W` on the
vector  space `V`, evaluated at `q`. These are the graded multiplicities of
the  irreducible characters of `W` in the quotient `SV/I` where `SV` is the
symmetric  algebra of `V` and `I` is the ideal generated by the homogeneous
invariants  of  positive  degree  in  `SV`.  The  ordering  of  the  result
corresponds to the ordering of the characters in 'CharTable(W)'.

```julia-repl
julia> fakedegrees(coxgroup(:A,2),Pol(:q))
3-element Vector{Pol{Int64}}:
 q³
 q²+q
 1
```
"""
function fakedegrees(W,q;recompute=false)
  if !recompute
    res=improve_type(map(p->fakedegree(W,p,q),charinfo(W)[:charparams]))
    if !any(isnothing,res) && !any(iszero,res) return res end
  end
  # need general routine
  InfoChevie("# recomputing fakedegrees for ",W,"\n")
  qq=Pol()
  P=generic_order(W,qq)
  P=shift(P,-valuation(P))
  ct=CharTable(W)
  if rank(W)==0 P=[one(qq)]
  else
    P=ct.irr*map(1:length(ct.centralizers))do i
      P/(ct.centralizers[i]*prod(l->(qq*E(inv(l))-1),refleigen(W,i)))
    end
  end
  charinfo(W)[:B]=degree.(P)
  charinfo(W)[:b]=valuation.(P)
  P=map(x->x(q),P)
  P=improve_type(P)
  W isa Spets ?  P.*(-1)^rank(W)*generic_sign(W) : P
end

function charinfo(t::TypeIrred)
  c=deepcopy(getchev(t,:CharInfo))
  c[:positionId]=c[:extRefl][1]
  c[:positionDet]=c[:extRefl][end]
  c[:charnames]=map(c[:charparams]) do p
     getchev(t,:CharName,p,Dict(:TeX=>true))
  end
  if !haskey(c,:b) c[:b]=getchev(t,:LowestPowerFakeDegrees) end
  if !haskey(c,:B) c[:B]=getchev(t,:HighestPowerFakeDegrees) end
  if !haskey(c,:a) c[:a]=getchev(t,:LowestPowerGenericDegrees) end
  if !haskey(c,:A) c[:A]=getchev(t,:HighestPowerGenericDegrees) end
  if isnothing(c[:a])
    uc=getchev(t,:UnipotentCharacters)
    if uc!=false && uc!==nothing
      c[:a]=uc[:a][uc[:harishChandra][1][:charNumbers]]
      c[:A]=uc[:A][uc[:harishChandra][1][:charNumbers]]
    end
  end
  for f in [:a,:A,:b,:B]
    if isnothing(c[f]) delete!(c,f) else c[f]=Int.(c[f]) end
  end
  if haskey(t,:orbit)
    if !haskey(c,:charRestrictions)
      c[:charRestrictions]=eachindex(c[:charparams])
      c[:nrGroupClasses]=length(c[:charparams]) # assume ortit twist trivial
    end
    for f in [:a,:A,:b,:B]
      if haskey(c,f) c[f]*=length(t.orbit) end
    end
  end
  c
end

cartfields(p,f)=cartesian(getindex.(p,f)...)

"""
`charinfo(W)`

returns   information  about  the  irreducible  characters  of  the  finite
reflection group `W`. The result is a Dict with the following entries:

`:charparams`:  contains  parameters  for  the  irreducible  characters  as
described in the introduction. The parameters are tuples with one component
for  each irreducible  component of  `W` (as  given by  `refltype`). For an
irreducible   component  which  is  an  imprimitive  reflection  group  the
component  of the `charparam` is a tuple of partitions (partitions for type
`:A`,  double partitions  for type  `:B`), and  for a primitive irreducible
group it is a pair `(d,e)` where `d` is the degree of the character and `e`
is  the  smallest  symmetric  power  of  the  character  of  the reflection
representation  which  contains  the  given  character  as  a component. In
addition,  there is an ordinal number if more than one character shares the
first two invariants.

```julia-repl
julia> charinfo(coxgroup(:G,2))[:charparams]
6-element Vector{Vector{Vector{Int64}}}:
 [[1, 0]]
 [[1, 6]]
 [[1, 3, 1]]
 [[1, 3, 2]]
 [[2, 1]]
 [[2, 2]]
```

`:charnames`:  strings describing the  irreducible characters, computed from
the `charparams`. This is the same as `charnames(W)`.

`:positionId`:  the position of the trivial character in the character table
of `W`.

```julia-repl
julia> charinfo(coxgroup(:D,4))[:positionId]
13
```

`:positionDet`:  Contains the position  of the determinant  character in the
character   table  of  `W`. For Coxeter groups this is the sign character.

```julia-repl
julia> charinfo(coxgroup(:D,4))[:positionDet]
4
```

`:extRefl`: Only present if `W` is irreducible, in which case the reflection
representation  of `W` and all its exterior powers are irreducible. It then
contains   the  position   of  the   exterior  powers   of  the  reflection
representation in the character table.

```julia-repl
julia> charinfo(coxgroup(:D,4))[:extRefl]
5-element Vector{Int64}:
 13
 11
  5
  3
  4
```

`:b`:   contains  a  list  holding  the  `b`-function  for  all  irreducible
characters  of `W`, that is,  for each character `χ`,  the valuation of the
fake  degree of `χ`. The ordering of the result corresponds to the ordering
of  the  characters  in  `CharTable(W)`.  The  advantage  of  this function
compared  to calling `fakeDegrees` is that one  does not have to provide an
indeterminate,  and that  it may  be much  faster to  compute than the fake
degrees.

```julia-repl
julia> charinfo(coxgroup(:D,4))[:b]
13-element Vector{Int64}:
  6
  6
  7
 12
  4
  3
  6
  2
  2
  4
  1
  2
  0
```

`:B`:   contains  a  list  holding  the  `B`-function  for  all  irreducible
characters  of `W`, that is, for each character `χ`, the degree of the fake
degree  of `χ`. The ordering  of the result corresponds  to the ordering of
the  characters in `CharTable(W)`. The  advantage of this function compared
to  calling  `fakeDegrees`  is  that  one  does  not  have  to  provide  an
indeterminate,  and that  it may  be much  faster to  compute than the fake
degrees.

```julia-repl
julia> charinfo(coxgroup(:D,4))[:B]
13-element Vector{Int64}:
 10
 10
 11
 12
  8
  9
 10
  6
  6
  8
  5
  6
  0
```

`:a`:  Only  filled  for  Spetsial  groups.  Contains  a  list  holding  the
`a`-function  for  all  irreducible  characters  of  the  Coxeter  group or
Spetsial  reflection  group  `W`,  that  is,  for  each  character `χ`, the
valuation  of the generic degree of `χ` (in the one-parameter Hecke algebra
`hecke(W,Pol(:q))`  corresponding  to  `W`).  The  ordering  of  the result
corresponds to the ordering of the characters in `CharTable(W)`.

```julia-repl
julia> charinfo(coxgroup(:D,4))[:a]
13-element Vector{Int64}:
  6
  6
  7
 12
  3
  3
  6
  2
  2
  3
  1
  2
  0
```

`:A`:  Only  filled  for  Spetsial  groups.  Contains  a  list  holding  the
`A`-function  for  all  irreducible  characters  of  the  Coxeter  group or
Spetsial  reflection group `W`, that is, for each character `χ`, the degree
of   the  generic  degree  of  `χ`  (in  the  one-parameter  Hecke  algebra
`hecke(W,Pol(:q))`  corresponding  to  `W`).  The  ordering  of  the result
corresponds to the ordering of the characters in `CharTable(W)`.

```julia-repl
julia> charinfo(coxgroup(:D,4))[:A]
13-element Vector{Int64}:
 10
 10
 11
 12
  9
  9
 10
  6
  6
  9
  5
  6
  0
```

`:opdam`:  Contains the permutation of  the characters obtained by composing
the  Opdam  involution  with  complex  conjugation. This permutation has an
interpretation as a Galois action on the characters of
`H=hecke(W,Pol(:x))`:  if `H` splits  by taking `v`  an `e`-th root of `x`,
`.opdam` records the permutation effected by the Galois action `v->E(e)*v`.

```julia-repl
julia> charinfo(ComplexReflectionGroup(22))[:opdam]
(3,5)(4,6)(11,13)(12,14)(17,18)
```

```julia-repl
julia> charinfo(coxgroup(:A,2))
Dict{Symbol, Any} with 9 entries:
  :a           => [3, 1, 0]
  :b           => [3, 1, 0]
  :positionId  => 3
  :charnames   => ["111", "21", "3"]
  :A           => [3, 2, 0]
  :B           => [3, 2, 0]
  :extRefl     => [3, 2, 1]
  :charparams  => [[[1, 1, 1]], [[2, 1]], [[3]]]
  :positionDet => 1
```

For  irreducible groups, the returned  record contains sometimes additional
information:

for  `F₄`: the entry `:kondo` gives the labeling of the characters given by
Kondo, also used in [Lusztig1985, (4.10)](biblio.htm#Lus85).

for  `E₆, E₇, E₈`: the entry `:frame`  gives the labeling of the characters
given   by  Frame,   also  used   in  [Lusztig1985,   (4.11),  (4.12),  and
(4.13)](biblio.htm#Lus85).

for  `G₂`: the  entry `:spaltenstein`  gives the  labeling of the characters
given by Spaltenstein.

```julia-repl
julia> charinfo(coxgroup(:G,2))[:spaltenstein]
6-element Vector{String}:
 "1"
 "\\varepsilon"
 "\\varepsilon_l"
 "\\varepsilon_c"
 "\\theta'"
 "\\theta''"
```

for  `G(de,e,2)`  even  `e`  and  `d>1`:  the  entry  `:malle`  gives  the
parameters for the characters used by Malle in [Malle1996](biblio.htm#Mal96).
"""
function charinfo(W)::Dict{Symbol,Any}
  get!(W,:charinfo)do
    p=charinfo.(refltype(W))
    if isempty(p)
      res=Dict(:a=>[0],:A=>[0],:b=>[0],:B=>[0],:positionId=>1,
      :positionDet=>1,:charnames=>["Id"],:extRefl=>[1],:charparams=>[[]])
      if W isa Spets
        res[:charRestrictions]=[1]
        res[:nrGroupClasses]=1
      end
      return res
    end
    if length(p)==1 res=copy(p[1]) else res=Dict{Symbol, Any}() end
    res[:charparams]=cartfields(p,:charparams)
    if W isa Spets
      gt=map(x->sort(x.indices),refltype(Group(W)))
      t=refltype(W)
      n=fill(0,length(gt))
      for i in eachindex(t), f in t[i].orbit
        n[findfirst(==(sort(f.indices)),gt)]=p[i][:nrGroupClasses]
      end
      res[:charRestrictions]=
      map(cartesian(getindex.(p,:charRestrictions)...))do y
        m=fill(0,length(gt))
        for i in eachindex(t), f in t[i].orbit
          m[findfirst(==(sort(f.indices)),gt)]=y[i]
        end
        return HasType.PositionCartesian(n,m)
      end
      res[:nrGroupClasses]=prod(i->p[i][:nrGroupClasses]^length(t[i].orbit),
                                                          eachindex(t))
    end
    if length(p)==1 return res end
    res[:charnames]=map(l->join(l,","),cartfields(p,:charnames))
    for f in [:positionId, :positionDet]
     if all(d->haskey(d,f),p)
       res[f]=HasType.PositionCartesian(map(x->length(x[:charparams]),p),getindex.(p,f))
      end
    end
    for f in [:b, :B, :a, :A]
      if all(d->haskey(d,f),p) res[f]=Int.(map(sum,cartfields(p,f))) end
    end
    if any(x->haskey(x, :opdam),p)
      res[:opdam]=map(x->haskey(x,:opdam) ? x[:opdam] : Perm(), p)
      gt=cartesian(map(x->1:length(x[:charparams]), p))
      res[:opdam]=PermListList(gt, map(t->map((x,i)->x^i,t,res[:opdam]),gt))
    end
    res
  end
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
  t=permutedims(CharTable(W).irr)
  Perm(t,t.*t[:,charinfo(W)[:positionDet]];dims=2)
end

function classinfo(t::TypeIrred)
  cl=deepcopy(getchev(t,:ClassInfo))
  if haskey(t,:orbit)
     l=length(t.orbit)
     t=t.orbit[1]
     if l>1 && haskey(cl,:classes)
       cl[:classes].*=prod(degrees(t))^(l-1)
     end
  end
  inds=t.indices
  cl[:classtext]=map(x->inds[x],cl[:classtext])
  if haskey(cl,:classes) cl[:classes]=Int.(cl[:classes]) end
  cl
end

Groups.nconjugacy_classes(t::TypeIrred)=getchev(t,:NrConjugacyClasses)

"""
`classinfo(W)`

returns  information about the  conjugacy classes of  the finite reflection
group `W`. The result is a Dict with three entries:

`:classtext`:  contains words in  the generators describing representatives
of  each  conjugacy  class.  Each  word  is  a  list  of integers where the
generator  `sᵢ`  is  represented  by  the  integer  `i`. For finite Coxeter
groups, it is the same as
`map(x->word(W,representative(x)),conjugacyclasses(W))`,   and   each  such
representative  is of minimal length in its  conjugacy class and is a "very
good" element in the sense of [GeckMichel1997](biblio.htm#GM97).

`:classparams`:  The  elements  of  this  list  are  tuples  which have one
component  for each irreducible component of  `W`. These components for the
infinite  series,  contain  partitions  or  partition tuples describing the
class  (see  the  introduction).  For  the  exceptional Coxeter groups they
contain  Carter's admissible  diagrams, see [Carter1972](biblio.htm#Car72).
For  exceptional complex reflection groups they contain in general the same
information as in classtext.

`:classnames`:  Contains strings describing the conjugacy classes, made out
of the information in `:classparams`.

```julia-repl
julia> classinfo(coxgroup(:A,2))
Dict{Symbol, Any} with 5 entries:
  :classes     => [1, 3, 2]
  :orders      => [1, 2, 3]
  :classtext   => [Int64[], [1], [1, 2]]
  :classnames  => ["111", "21", "3"]
  :classparams => [[1, 1, 1], [2, 1], [3]]
```

See also the introduction of this section.
"""
function classinfo(W)::Dict{Symbol,Any}
  get!(W,:classinfo)do
    tmp=map(classinfo,refltype(W))
    if isempty(tmp) return Dict(:classtext=>[Int[]],:classnames=>[""],
                      :classparams=>[Int[]],:orders=>[1],:classes=>[1])
    end
    if any(isnothing, tmp) return nothing end
    if length(tmp)==1 res=copy(tmp[1]) else res=Dict{Symbol, Any}() end
    res[:classtext]=map(x->reduce(vcat,x),cartfields(tmp,:classtext))
    res[:classnames]=map(join,cartfields(tmp,:classnames))
    if all(haskey.(tmp,:classparam))
      res[:classparams]=cartfields(tmp,:classparams)
    end
    if all(haskey.(tmp,:orders))
      res[:orders]=map(lcm, cartfields(tmp,:orders))
    end
    if all(haskey.(tmp,:classes))
      res[:classes]=map(prod, cartfields(tmp,:classes))
    end
    res
  end
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
  show(IOContext(io,:TeX=>true), "text/plain",ct)
end

function Base.show(io::IO,t::CharTable)
  if get(io,:TeX,false) || get(io,:limit,false) || !haskey(t,:repr) 
    printTeX(io,"CharTable(\$",t.name,"\$)")
  else  print(io,t.repr)
  end
end

function Base.show(io::IO, ::MIME"text/plain", ct::CharTable)
  println(io,ct)
  irr=map(ct.irr)do e
    if iszero(e) "." else repr(e; context=io) end
  end
  showtable(io,irr,row_labels=map(s->fromTeX(io,s),ct.charnames),
                col_labels=map(s->fromTeX(io,s),ct.classnames))
end

function CharTable(t::TypeIrred)
  ct=getchev(t,:CharTable)
  if haskey(ct,:irredinfo) names=getindex.(ct[:irredinfo],:charname)
  else                     names=charinfo(t)[:charnames]
  end
  if !haskey(ct,:classnames) merge!(ct,classinfo(t)) end
  irr=toM(ct[:irreducibles])
  if all(isinteger,irr) irr=Int.(irr)
  else irr=Cyc{Int}.(irr)
  end
  CharTable(irr,names,String.(ct[:classnames]),Int.(ct[:centralizers]),
            ct[:size],Dict{Symbol,Any}(:name=>repr(t;context=:TeX=>true)))
end

function Base.prod(ctt::Vector{<:CharTable})
  if isempty(ctt)
    return CharTable(hcat(1),["Id"],["1"],[1],1,Dict{Symbol,Any}(:name=>"."))
  end
  charnames=join.(cartesian(getfield.(ctt,:charnames)...),",")
  classnames=join.(cartesian(getfield.(ctt,:classnames)...),",")
  centralizers=prod.(cartesian(getfield.(ctt,:centralizers)...))
  order=prod(getfield.(ctt,:order))
  if length(ctt)==1 irr=ctt[1].irr
  else irr=kron(getfield.(ctt,:irr)...)
  end
  CharTable(irr,charnames,classnames,centralizers,order,Dict{Symbol,Any}(:name=>"x"))
end

"""
`CharTable(WF::Spets)`

This  function returns the character table of the reflection coset `WF`. We
call "characters" of the coset `WF=W.ϕ` of the group `W` the restriction to
`W.ϕ`  of a set containing one extension of each `ϕ`-invariant character of
W  to the semidirect  product of W  with the cyclic  group generated by `ϕ`
(for  Coxeter  cosets  we  choose,  following  Lusztig,  in  each  case one
extension, called the preferred extension.)

```julia-repl
julia> W=spets(coxgroup(:D,4),Perm(1,2,4))
³D₄

julia> CharTable(W)
CharTable(³D₄)
     │C₃ Ã₂ C₃+A₁ Ã₂+A₁ F₄ Ã₂+A₂ F₄(a₁)
─────┼──────────────────────────────────
.4   │ 1  1     1     1  1     1      1
.1111│-1  1     1    -1  1     1      1
.22  │ .  2     2     . -1    -1     -1
11.2 │ .  .     .     . -1     3      3
1.3  │ 1  1    -1    -1  .    -2      2
1.111│-1  1    -1     1  .    -2      2
1.21 │ .  2    -2     .  .     2     -2
```
"""
function CharTable(W::Union{PermRootGroup,FiniteCoxeterGroup,Spets})::CharTable
  get!(W,:chartable)do
    t=refltype(W)
    ct=isempty(t) ? 
      CharTable(fill(1,1,1),["Id"],["."],[1],1,Dict{Symbol,Any}()) :
      prod(CharTable.(t))
    ct.name=repr(W;context=:TeX=>true)
    ct.repr="CharTable($W)"
    ct
  end
end

function classes(ct::CharTable)
  get!(ct,:classes)do
    div.(ct.order,ct.centralizers)
  end
end

function scalarproduct(ct::CharTable,c1::AbstractVector,c2::AbstractVector)
  div(c2'*(c1.*classes(ct)),ct.centralizers[1])
end

"""
`decompose(ct::CharTable,c::Vector)` 

decompose character `c` (given by its values on conjugacy classes) 
on irreducible characters as given by `CharTable` `ct`
"""
function decompose(ct::CharTable,c::Vector)
  map(i->scalarproduct(ct,view(ct.irr,i,:),c),axes(ct.irr,1))
end

"""
`on_chars(G,aut)`

`aut`  is an automorphism of  the group `G` (for  a permutation group, this
could  be  given  as  a  permutation  normalizing  `G`).  The result is the
permutation  `1:nconjugacy_classes(G)` (representing indices of irreducible
characters) induced by `aut`.

```julia-repl
julia> WF=rootdatum("3D4")
³D₄

julia> on_chars(Group(WF),WF.phi)
(1,2,7)(8,9,12)
```
"""
function on_chars(W,aut)
  ct=CharTable(W).irr
  inv(Perm(ct,^(ct,on_classes(W, aut),dims=2),dims=1))
end

CharTable(W::CoxSym)=CharTable(refltype(W)[1])
"""
`representation(W,i)`

returns a list holding, for the `i`-th irreducible character of the complex
reflection  group  `W`,  a  list  of  matrices  images  of  the  generating
reflections  of `W`  in a  model of  the corresponding representation. This
function  is based on the classification,  and is not yet fully implemented
for   `G₃₄`;  78  representations   are  missing  out   of  169,  that  is,
representations  of dimension ≥140, except half of those of dimensions 315,
420 and 840.

```julia-repl
julia> representation(ComplexReflectionGroup(24),3)
3-element Vector{Matrix{Cyc{Int64}}}:
 [1 0 0; -1 -1 0; -1 0 -1]
 [-1 0 -1; 0 -1 (1-√-7)/2; 0 0 1]
 [-1 -1 0; 0 1 0; 0 (1+√-7)/2 -1]
```
"""
function representation(W::Union{Group,Spets},i::Int)
  dims=Tuple(getchev(W,:NrConjugacyClasses))
  if isempty(dims) return Matrix{Int}[] end
  tt=refltype(W)
  inds=reverse(Tuple(CartesianIndices(reverse(dims))[i]))
  mm=map((t,j)->getchev(t,:Representation,j),tt,inds)
  if any(isnothing,mm) error("no representation for ",W) end
  if W isa Spets
    FF=map(x->x[:F],mm)
    if !(FF[1] isa Matrix) FF=map(toM,FF) end
    F=length(FF)==1 ? FF[1] : kron(FF...)
    mm=map(x->x[:gens],mm)
  end
  if !(mm[1][1] isa Matrix) mm=map(x->toM.(x),mm) end
  mm=improve_type.(mm)
  n=length(tt)
  if n==1 reps=mm[1] 
  else reps=vcat(map(1:n) do i
           map(mm[i]) do m
             kron(map(j->j==i ? m : mm[j][1]^0,1:n)...)
           end
         end...)
  end
  if !(W isa Spets) return reps end
  (gens=reps,F=F)
end

"""
`representations(W)`

returns the representations of `W` (see `representation`).

```julia-repl
julia> representations(coxgroup(:B,2))
5-element Vector{Vector{Matrix{Int64}}}:
 [[1], [-1]]
 [[1 0; -1 -1], [1 2; 0 -1]]
 [[-1], [-1]]
 [[1], [1]]
 [[-1], [1]]
```
"""
representations(W::Union{Group,Spets})=representation.(Ref(W),1:nconjugacy_classes(W))

"""
`WGraphToRepresentation(coxrank::Integer,graph,v)`

(Jean Michel june/december 2003 from  code/data of Geck, Marin, Alvis,
Naruse, Howlett,Yin)
We  store  some  representations  of  some  Hecke  algebras  with  only one
parameter  `v` as `W`-graphs. For a  Coxeter system `(W,S)` where `coxrank`
is  the length of `S`, a `W`-graph is  defined by a set of vertices `C`; to
`x∈  C` is  attached `I(x)⊂  S` and  to `(x,y)∈  C^2` is attached an "edge"
`μ(x,y)`  in the field of definition  of `W`; this defines a representation
of  the Hecke algebra with  single rootparameter `v` on  a space with basis
``{e_y}_{y∈ C}`` by:

``T_s(e_y)=\\begin{cases}-e_y& if s∈ I(y)\\\\
       v^2 e_y+∑_{x∣s∈ I(x)} vμ(x,y)e_x&otherwise\\end{cases}``

The  `W`-graphs are  stored in  a compact  format to  save space.  They are
represented as a pair.
  - The  first element is a list describing `C`.
Its  elements are  either a  set `I(x)`,  or an  integer `n`  specifying to
repeat the previous element `n` more times.

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
"""
function WGraphToRepresentation(rk::Integer,gr::Vector,v)
  V=Vector{Int}[]
  for S in gr[1]
    if S isa Integer append!(V,map(i->V[end],1:S))
    else push!(V,S)
    end
  end
  dim=length(V)
  T=Int
  function prom(a)
    if a isa Vector for u in a prom(u) end
    else T=promote_type(T,typeof(a))
    end
  end
  prom(gr[2])
  S=map(i->one(fill(T(0),dim,dim))*v^2,1:rk)
  for j in 1:dim for i in V[j] S[i][j,j]=-one(v) end end
  for i in gr[2]
    if i[1] isa Vector mu=i[1] else mu=[i[1],i[1]] end
    for l in i[2]
      x=l[1]
      for y in l[2:end]
        for j in setdiff(V[y],V[x]) S[j][y,x]=mu[2]*v end
        for j in setdiff(V[x],V[y]) S[j][x,y]=mu[1]*v end
      end
    end
  end
  map(toL,S)
end

############################################################################
# How to interpret W-graphs for complex reflection groups with one orbit of
# reflections, for hecke(W,[vars]).

function WGraph2Representation(a,vars)
# println("a=$a vars=$vars")
  nodes=a[1]
  pos=function(n,j)
    if n[1] isa Vector p=findfirst(x->j in x,n)
      if isnothing(p) p=length(vars) end
    elseif j in n p=1
    else p=2 end
    p
  end
  flat(l)=l[1] isa Vector ? flat(reduce(vcat,l)) : l
  rk=maximum(Int.(flat(nodes))) # number of generators
  dim=length(nodes)
  R=map(j->map(k->vars[pos(nodes[k],j)],1:dim),1:rk)
  R=map(x->toM(HasType.DiagonalMat(x...)),R)
  R=map(x->x.+0*E(1)//1,R)
# println("R=$(typeof(R))$R")
  for r in a[2]
#   println("r=$r")
    for k in [3,4]
    if HasType.IsList(r[k])
      for j in 2:2:length(r[k]) R[Int(r[k][j-1])][r[k-2],r[5-k]]=r[k][j] end
    else
      r2=Int.(r[1:2])
      j=filter(i->pos(nodes[r2[k-2]],i)<pos(nodes[r2[5-k]],i),1:rk)
      for i in j R[i][r2[k-2],r2[5-k]]=r[k] end
    end
    end
  end
# println("R=$(typeof(R))$R")
  toL.(R)
end

# the next function returns the dual W-graph of gr (for an Hecke algebra of
# rank rk). A dual W-graph corresponds to a Curtis Dual representation.
function DualWGraph(rk,gr)
  [map(x->x isa Integer ? x : setdiff(1:rk,x),gr[1]),
   map(((x,y),)->x isa Vector ? [-reverse(x),y] : [-x,y],gr[2])]
end

"""
`charnames(W;options...)`
`charnames(io::IO,W)`

returns  the  list  of  character  names  for the reflection group `W`. The
optional  options are IOContext attributes which can give alternative names
in  certain cases, or a different formatting  of names in general. They can
be specified by giving an IO as argument.

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

The  last two  commands show  the character  names used by Spaltenstein and
Lusztig when describing the Springer correspondence.
"""
function charnames(io::IO,W)
  if hasmethod(refltype,(typeof(W),))
    c=charinfo(W)
    cn=c[:charnames]
    for k in [:spaltenstein, :frame, :malle, :kondo]
      if get(io,k,false) && haskey(c,k) cn=c[k] end
    end
  else
    cn=CharTable(W).charnames
  end
  fromTeX.(Ref(io),cn)
end

charnames(W;opt...)=charnames(IOContext(stdout,opt...),W)

@GapObj struct InductionTable{T}
  scalar::Matrix{T}
  gcharnames::Vector{String}
  ucharnames::Vector{String}
  identifier::String
end

"""
   `InductionTable(u,g)`

returns   an  object  describing  the   decomposition  of  the  irreducible
characters  of the subgroup  `u` induced to  the group `g`.  In the default
show method, the rows correspond to the characters of the parent group, and
the  columns  to  those  of  the  subgroup.  The  return object has a field
`scalar`  which is a `Matrix{Int}` containing  the induction table, and the
other  fields contain labeling information  taken from the character tables
of `u` and `g` when it exists.

```julia-rep1        needs Gap4
julia> g=Group([Perm(1,2),Perm(2,3),Perm(3,4)])
Group([(1,2),(2,3),(3,4)])

julia> u=Group( [ Perm(1,2), Perm(3,4) ])
Group([(1,2),(3,4)])

julia> InductionTable(u,g)
Induction Table from Group([perm"(1,2)",perm"(3,4)"]) to Group([perm"(1,2)",perm"(2,3)",perm"(3,4)"])
   │X.1 X.2 X.3 X.4
───┼────────────────
X.1│  .   1   .   .
X.2│  .   1   1   1
X.3│  1   1   .   .
X.4│  1   .   1   1
X.5│  1   .   .   .
```

```julia-repl
julia> g=coxgroup(:G,2)
G₂

julia> u=reflection_subgroup(g,[1,6])
G₂₍₁₅₎=A₂

julia> t=InductionTable(u,g)
Induction Table from G₂₍₁₅₎=A₂ to G₂
     │111 21 3
─────┼─────────
φ₁‚₀ │  .  . 1
φ₁‚₆ │  1  . .
φ′₁‚₃│  1  . .
φ″₁‚₃│  .  . 1
φ₂‚₁ │  .  1 .
φ₂‚₂ │  .  1 .
```

using an `IOContext` allows to transmit attributes to the table format method

```julia-rep1
julia> xprint(t;rows=[5],cols=[3,2])
Induction Table
    │3 21
────┼─────
φ₂‚₁│.  1
```

It is also possible to TeX induction tables

##  This function also works for Spets (Reflection Cosets)
"""
function InductionTable(u,g)
  tu=CharTable(u)
  tg=CharTable(g)
  f=fusion_conjugacy_classes(u,g)
  cl=div.(length(u),tu.centralizers)
  scal(c,c1)=sum(map(*,conj.(c),c1,cl))//length(u)
  lu=repr(u;context=:TeX=>true)
  lg=repr(g;context=:TeX=>true)
  sc=[scal(tg.irr[j,f],tu.irr[i,:]) for j in axes(tg.irr,1),i in axes(tu.irr,1)]
  InductionTable(improve_type(sc),tg.charnames,tu.charnames,
  "Induction Table from \$$lu\$ to \$$lg\$",
  Dict{Symbol,Any}(:repr=>"InductionTable($(repr(u)),$(repr(g)))"))
end

function Base.show(io::IO, ::MIME"text/html", t::InductionTable)
  show(IOContext(io,:TeX=>true),"text/plain",t)
end

function Base.show(io::IO,t::InductionTable)
  if !get(io,:TeX,false) && !get(io,:limit,false) && haskey(t,:repr) 
    print(io,t.repr)
  else printTeX(io,t.identifier)
  end
end

function Base.show(io::IO,::MIME"text/plain",t::InductionTable)
  println(io,t)
  column_labels=fromTeX.(Ref(io),t.ucharnames)
  row_labels=fromTeX.(Ref(io),t.gcharnames)
  scal=map(t.scalar)do e
    if iszero(e) "." else repr(e; context=io) end
  end
  showtable(io,scal,row_labels=row_labels,col_labels=column_labels)
end

"""
`jInductionTable(H, W)`

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

julia> jInductionTable(H,W)
j-Induction Table from D₄₍₁₃₎=A₂Φ₁² to D₄
     │111 21 3
─────┼─────────
11+  │  .  . .
11-  │  .  . .
1.111│  .  . .
.1111│  .  . .
11.2 │  .  . .
1.21 │  1  . .
.211 │  .  . .
2+   │  .  . .
2-   │  .  . .
.22  │  .  . .
1.3  │  .  1 .
.31  │  .  . .
.4   │  .  . 1
```
"""
function jInductionTable(u,g)
  tbl=InductionTable(u,g)
  bu=charinfo(u)[:b]
  bg=charinfo(g)[:b]
  t=copy(tbl.scalar)
  for (i,bi) in enumerate(bu), (j,bj) in enumerate(bg)
    if bi!=bj t[j,i]=0 end
  end
  lu=repr(u;context=:TeX=>true)
  lg=repr(g;context=:TeX=>true)
  InductionTable(t,tbl.gcharnames,tbl.ucharnames,
  "j-Induction Table from \$$lu\$ to \$$lg\$",
  Dict{Symbol,Any}(:repr=>"jInductionTable($(repr(u)),$(repr(g)))"))
end

"""
`JInductionTable(H, W)`

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

julia> JInductionTable(H,W)
J-Induction Table from D₄₍₁₃₎=A₂Φ₁² to D₄
     │111 21 3
─────┼─────────
11+  │  .  . .
11-  │  .  . .
1.111│  .  . .
.1111│  .  . .
11.2 │  1  . .
1.21 │  1  . .
.211 │  .  . .
2+   │  .  . .
2-   │  .  . .
.22  │  .  . .
1.3  │  .  1 .
.31  │  .  . .
.4   │  .  . 1
```
"""
function JInductionTable(u,g)
  tbl=InductionTable(u,g)
  bu=charinfo(u)[:a]
  bg=charinfo(g)[:a]
  t=copy(tbl.scalar)
  for (i,bi) in enumerate(bu), (j,bj) in enumerate(bg)
    if bi!=bj t[j,i]=0 end
  end
  lu=repr(u;context=:TeX=>true)
  lg=repr(g;context=:TeX=>true)
  InductionTable(t,tbl.gcharnames,tbl.ucharnames,
  "J-Induction Table from \$$lu\$ to \$$lg\$",
  Dict{Symbol,Any}(:repr=>"JInductionTable($(repr(u)),$(repr(g)))"))
end
end
