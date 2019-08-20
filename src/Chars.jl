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
exceptional  Coxeter  groups,  Carter's  admissible  diagrams [@Car72] (for
other  primitive  complex  reflection  groups  we  just  use  words  in the
generators  to specify  the classes).  For the  characters, these are again
partitions  or partition tuples for the infinite series, and for the others
they  are pairs  of two  integers `(d,e)`  where `d`  is the  degree of the
character  and  `e`  is  the  smallest  symmetric  power  of the reflection
representation  containing  the  given  character  as  a  constituent  (the
`b`-invariant  of the character). This information is obtained by using the
functions `classinfo` and `charinfo`. When you display the character table,
the canonical labelings for classes and characters are those displayed.

A   typical   example   is   `coxgroup(:A,n)`,  the  symmetric  group
`𝔖ₙ₊₁`  where classes and characters  are parameterized by partitions
of `n+1`.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> CharTable(W)
CharTable(H(G(1,1,4)))
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
CharTable(W(G2))
     │A₀ Ã₁ A₁ G₂ A₂ A₁+Ã₁
─────┼─────────────────────
φ₁‚₀ │ 1  1  1  1  1     1
φ₁‚₆ │ 1 -1 -1  1  1     1
φ′₁‚₃│ 1  1 -1 -1  1    -1
φ″₁‚₃│ 1 -1  1 -1  1    -1
φ₂‚₁ │ 2  .  .  1 -1    -2
φ₂‚₂ │ 2  .  . -1 -1     2


julia> ct.charnames
6-element Array{String,1}:
 "\\phi_{1,0}"  
 "\\phi_{1,6}"  
 "\\phi_{1,3}'" 
 "\\phi_{1,3}''"
 "\\phi_{2,1}"  
 "\\phi_{2,2}"  

julia> ct.classnames
6-element Array{String,1}:
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
"LowestPowerfakeDegrees" and "HighestPowerfakeDegrees"). For finite Coxeter
groups,   the  valuation  and   degree  of  the   generic  degrees  of  the
one-parameter  generic Hecke algebra give two  more integers `a,A` (see the
functions  "LowestPowerGenericDegrees",  "HighestPowerGenericDegrees",  and
[@Car85, Ch.11]   for  more details).  These will  also be  used in the
operations  of truncated  inductions explained  in the  chapter "Reflection
subgroups".

Iwahori-Hecke  algebras and  cyclotomic Hecke  algebras also have character
tables, see the corresponding chapters.

We  now describe for each type our conventions for labeling the classes and
characters.

{Type  `Aₙ` (`n ≥  0`)}. In this  case we have  `W ≅ 𝔖ₙ₊₁`. The classes and
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
`(1^{n+1})`  and  the  trivial  character  to  the  partition  `(n+1)`. The
character of the reflection representation of `W` is labeled by `(n,1)`.

{Type  `Bₙ` (`n ≥ 2`)}. In this  case `W=W(Bₙ)` is isomorphic to the wreath
product  of the cyclic  group of order  `2` with the  symmetric group `𝔖ₙ`.
Hence  the classes and characters are  parameterized by pairs of partitions
such  that the total sum of their  parts equals `n`. The pair corresponding
to  a class describes the signed cycle type for the elements in that class,
as  in [@Car72]. We use the convention that  if `(λ,μ)` is such a pair then
`λ`  corresponds  to  the  positive  and  `μ` to the negative cycles. Thus,
`(1^n,-)`  and `(-,1^n)` label  the trivial class  and the class containing
the longest element, respectively. The pair corresponding to an irreducible
character is determined via Clifford theory, as follows.

We  have a semidirect product decomposition `W(Bₙ)=N ⋊ 𝔖ₙ` where `N` is the
standard  `n`-dimensional  `𝔽₂^n`-vector  space.  For  `a,b  ≥ 0` such that
`n=a+b` let `η_{a,b}` be the irreducible character of `N` which takes value
`1`  on the first `a` standard basis vectors and value `-1` on the next `b`
standard  basis vectors of `N`. Then  the inertia subgroup of `η_{a,b}` has
the  form `T_{a,b}:=N.(𝔖_a × 𝔖_b)` and we can extend `η_{a,b}` trivially to
an  irreducible  character  `η̃_{a,b}`  of  `T_{a,b}`.  Let  `α` and `β` be
partitions  of `a` and `b`, respectively. We take the tensor product of the
corresponding  irreducible characters of `𝔖_a` and `𝔖_b` and regard this as
an  irreducible  character  of  `T_{a,b}`.  Multiplying this character with
`η̃_{a,b}`  and  inducing  to  `W(Bₙ)`  yields an irreducible character `χ=
χ_{(α,β)}`  of `W(Bₙ)`. This defines the correspondence between irreducible
characters and pairs of partitions as above.

For   example,  the  pair  `((n),-)`   labels  the  trivial  character  and
`(-,(1^n))`  labels  the  sign  character.  The  character  of  the natural
reflection representation is labeled by `((n-1),(1))`.

{Type `Dₙ` (`n ≥ 4`)}. In this case `W=W(Dₙ)` can be embedded as a subgroup
of index `2` into the Coxeter group `W(Bₙ)`. The intersection of a class of
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
equal  to  `n`.  If  `α  !=  β`  then  the  restrictions of the irreducible
characters  of `W(Bₙ)` labeled by `(α,β)`  and `(β, α)` are irreducible and
equal.  If `α=β` then  the restriction of  the character labeled by `(α,α)`
splits  into  two  irreducible  components  which  we denote by `(α,+)` and
`(α,-)`. Note that this can only happen if `n` is even. In order to fix the
notation  we use  a result  of [@Ste89]  which describes  the value  of the
difference  of these two characters on a class of the form `(λ,+)` in terms
of the character values of the symmetric group `𝔖ₙ_{/2}`. Recall that it is
implicit  in the notation `(λ,+)` that all  parts of `λ` are even. Let `λ'`
be  the partition of `n/2` obtained by  dividing each part by `2`. Then the
value  of `χ_{(α,-)}-χ_{(α,+)}` on an element in the class `(λ,+)` is given
by  `2^{k(λ)}` times  the value  of the  irreducible character of `𝔖ₙ_{/2}`
labeled  by `α` on the class of  cycle type `λ'`. (Here, `k(λ)` denotes the
number of non-zero parts of `λ`.)

The  labels for the trivial, the  sign and the natural reflection character
are the same as for `W(Bₙ)`, since these characters are restrictions of the
corresponding characters of `W(Bₙ)`.

{The groups `G(d,1,n)`}.
They  are isomorphic to the wreath product of the cyclic group of order `d`
with  the  symmetric  group  `𝔖ₙ`.  Hence  the  classes  and characters are
parameterized  by `d`-tuples of partitions such that the total sum of their
parts  equals `n`. The words chosen  as representatives of the classes are,
when `d>2`, computed in a slightly different way than for `Bₙ`, in order to
agree  with the words on which Ram  and Halverson compute the characters of
the  Hecke algebra. First the parts of the `d` partitions are merged in one
big  partition and sorted in  increasing order. Then, to  a part `i` coming
from  the `j`-th partition  is associated the  word `(l+1…1… l+1)^{j-1}l+2…
l+i` where `l` is the highest generator used to express the previous part.

The  `d`-tuple corresponding to an  irreducible character is determined via
Clifford  theory in  a similar  way than  for the  `Bₙ` case.  The identity
character  has the first  partition with one  part equal `n`  and the other
ones  empty. The character of the  reflection representations has the first
two  partitions with one part  equal respectively to `n-1`  and to `1`, and
the other partitions empty.

{The groups `G(de,e,n)`}.
They  are normal  subgroups of  index `e`  in `G(de,1,n)`.  The quotient is
cyclic,  generated by the image `g`  of the first generator of `G(de,1,n)`.
The  classes are parameterized as the  classes of `G(de,e,n)` with an extra
information for a component of a class which splits.

According  to  [@Hu85],  a  class  `C`  of  `G(de,1,n)`  parameterized by a
`de`-partition  `(S₀,…,S_{de-1})` is in  `G(de,e,n)` if `e`  divides `∑_i i
∑_{p∈ S_i}p`. It splits in `d` classes for the largest `d` dividing `e` and
all  parts of all `S_i` and such that `S_i` is empty if `d` does not divide
`i`.  If `w` is in `C` then 'gⁱ w g⁻ⁱ' for 'i in 0:d-1' are representatives
of  the  classes  of  `G(de,e,n)`  which  meet  `C`.  They are described by
appending the integer `i` to the label for `C`.

The  characters are described by Clifford theory. We make `g` act on labels
for  characters of `G(de,1,n)`  . The action  of `g` permutes circularly by
`d`  the partitions in the `de`-tuple.  A character has same restriction to
`G(de,e,n)`  as its transform by `g`.  The number of irreducible components
of its restriction is equal to the order `k` of its stabilizer under powers
of  `g`.  We  encode  a  character  of  `G(de,e,n)`  by first, choosing the
smallest  for lexicographical order label  of a character whose restriction
contains  it; then this label is periodic with a motive repeated `k` times;
we  represent the  character by  one of  these motives,  to which we append
`E(k)^i` for 'i in 0:k-1' to describe which component of the restriction we
choose.

{Types  `G₂` and `F₄`}. The matrices  of character values and the orderings
and  labelings of  the irreducible  characters are  exactly the  same as in
[@Car85,  p.412/413]: in type `G₂` the character `φ₁,₃'` takes the value -1
on  the reflection associated  to the long  simple root; in  type `F₄`, the
characters  `φ₁,₁₂'`, `φ₂,₄'`, `φ₄,₇'`, `φ₈,₉'` and `φ₉,_{6}'` occur in the
induced  of the  identity from  the `A₂`  corresponding to the short simple
roots;  the  pairs  (`φ₂,₁_{6}'`,  `φ₂,₄″`)  and  (`φ₈,₃'`, `φ₈,₉″`) are
related by tensoring by sign; and finally `φ₆,_{6}″` is the exterior square
of  the reflection representation. Note, however, that we put the long root
at the left of the Dynkin diagrams to be in accordance with the conventions
in [@Lus85, (4.8) and (4.10)].

The  classes  are  labeled  by  Carter's  admissible  diagrams  [@Car72]. A
character is labeled by a pair `(d,b)` where `d` denotes the degree and `b`
the  corresponding `b`-invariant. If there  are several characters with the
same pair `(d,b)` we attach a prime to them, as in [@Car85].

{Types  `E₆,E₇,E₈`}. The character tables are obtained by specialization of
those  of the Hecke algebra. The classes are labeled by Carter's admissible
diagrams  [@Car72]. A  character is  labeled by  the pair `(d,b)` where `d`
denotes  the degree and  `b` is the  corresponding `b`-invariant. For these
types, this gives a unique labeling of the characters.

{Non-crystallographic  types `I₂(m)`, `H₃`, `H₄`}. In these cases we do not
have  canonical  labelings  for  the  classes.  We  label  them  by reduced
expressions.

Each  character for  type `H₃`  is uniquely  determined by the pair `(d,b)`
where  `d` is the degree and  `b` the corresponding `b`-invariant. For type
`H₄`  there are just  two characters (those  of degree `30`)  for which the
corresponding  pairs are  the same.  These two  characters are nevertheless
distinguished  by  their  fake  degrees:  the  character `φ₃₀,₁₀'` has fake
degree  `q^{10}+q^{12}+`  higher  terms,  while  `φ₃₀,₁₀″`  has fake degree
`q^{10}+q^{14}+`  higher terms. The  characters in the  table for type `H₄`
are ordered in the same way as in [@AL82].

Finally,  the characters  of degree `2`  for type  `I₂(m)` are  ordered as
follows.  The matrix representations affording the characters of degree `2`
are given by:
` ρ_j : s₁s₂ ↦
(\\begin{array}{cc}E(m)^j&0\\0&E(m)^{-j}\\end{array}),
 s₁↦(\\begin{array}{cc}0&1\\1&0\\end{array}),`
where  `1 ≤ j ≤  ⌊(m-1)/2⌋`. The reflection representation
is  `ρ₁`. The  characters in  the table  are ordered by listing
first the characters of degree 1 and then `ρ₁,ρ₂,…`.

{Primitive complex reflection groups `G₄` to `G₃₄`}.
The  groups `G₂₃=H₃`, `G₂₈=F₄`, `G₃₀=H₄` are exceptional Coxeter groups and
have  been  explained  above.  Similarly  for  the  other groups labels for
characters  consist primarily  of the  pair `(d,b)`  where `d`  denotes the
degree  and `b` is the corresponding  `b`-invariant. This is sufficient for
`G₄`,  `G₁₂`, `G₂₂` and `G₂₄`. For other  groups there are pairs or triples
of  characters which  have the  same `(d,b)`  value. We  disambiguate these
according  to  the  conventions  of  [@Mal00]  for `G₂₇, G₂₉, G₃₁, G₃₃` and
`G₃₄`:

-  For `G₂₇`:
The  fake degree  of `φ₃,₅'`  (resp. `φ₃,₂₀'`,  `φ₈,₉″`) has smaller degree
that  of `φ₃,₅″` (resp. `φ₃,₂₀″`,  `φ₈,₉'`). The characters `φ₅,_{15}'` and
`φ₅,_{6}'`  occur  with  multiplicity  1  in  the  induced from the trivial
character of the parabolic subgroup of type `A₂` generated by the first and
third generator (it is asserted mistakenly in [@Mal00] that `φ₅,_{6}″` does
not occur in this induced; it occurs with multiplicity 2).

-  For `G₂₉`:
The  character  `φ₆,₁₀‴`  is  the  exterior  square  of `φ₄,₁`; its complex
conjugate is `φ₆,₁₀⁗`. The character `φ₁₅,₄″` occurs in `φ₄,₁⊗φ₄,₃`; the
character  `φ₁₅,_{12}″` is  tensored by  the sign  character from `φ₁₅,₄″`.
Finally  `φ₆,₁₀'` occurs in  the induced from  the trivial character of the
standard parabolic subgroup of type `A₃` generated by the first, second and
fourth generators.

-  For `G₃₁`:
The   characters  `φ₁₅,_{8}'`,   `φ₁₅,_{20}'`  and   `φ₄₅,_{8}″`  occur  in
`φ₄,₁⊗φ₂₀,_{7}`;   the  character  `φ₂₀,₁₃'`  is  complex  conjugate  of
`φ₂₀,_{7}`;  the character `φ₄₅,_{12}'` is tensored by sign of `φ₄₅,_{8}'`.
The  two  terms  of  maximal  degree  of  the  fakedegree  of `φ₃₀,₁₀'` are
`q^{50}+q^{46}` while for `φ₃₀,₁₀″` they are `q^{50}+2q^{46}`.

-  For `G₃₃`:
The   terms  of   maximal  degree   of  the   fakedegree  of  `φ₁₀,₈'`  are
`q^{28}+q^{26}`  while for `φ₁₀,₈'` they  are `q^{28}+q^{24}`. The terms of
maximal  degree of the fakedegree of `φ₄₀,₅'` are `q^{31}+q^{29}` while for
`φ₄₀,₅″`  they are `q^{31}+2q^{29}`. The character `φ₁₀,₁₇'` is tensored by
sign of `φ₁₀,₈'` and `φ₄₀,₁₄'` is tensored by sign of `φ₄₀,₅'`.

-  For `G₃₄`:
The  character  `φ₂₀,₃₃'`  occurs  in  `φ₆,₁⊗φ₁₅,_{14}`.  The  character
`φ₇₀,₉'`  is rational.  The character  `φ₇₀,₉″` occurs in `φ₆,₁⊗φ₁₅,_{14}`.
The  character  `φ₇₀,₄_{5}'`  is  rational.The  character  `φ₇₀,₄_{5}″`  is
tensored   by  the   determinant  character   of  `φ₇₀,₉″`.  The  character
`φ₅_{60,18}'`   is   rational.   The   character  `φ₅_{60,18}‴`  occurs  in
`φ₆,₁⊗φ₃₃₆,₁_{7}`. The character `φ₂₈₀_{,12}'` occurs in `φ₆,₁⊗φ₃₃₆,₁_{7}`.
The  character  `φ₂₈₀_{,30}″`  occurs  in  `φ₆,₁⊗φ₃₃₆,₁_{7}`. The character
`φ₅₄₀_{,21}'`  occurs in  `φ₆,₁⊗φ₁₀₅,_{20}`. The  character `φ₁₀₅,_{8}'` is
complex  conjugate of `φ₁₀₅,₄`,  and `φ₈₄₀_{,13}'` is  complex conjugate of
`φ₈₄₀_{,11}`.   The  character   `φ₈₄₀_{,23}'`  is   complex  conjugate  of
`φ₈₄₀_{,19}`.  Finally  `φ₁₂₀,_{21}'`  occurs  in  induced from the trivial
character  of the standard parabolic subgroup of type `A₅` generated by the
generators of `G₃₄` with the third one omitted.

For the groups `G₅` and `G₇` we adopt the following conventions.
For `G₅` they are compatible with those of [@MR03] and [@BMM14].

-  For `G₅`:
We  let |W:=ComplexReflectionGroup(5)|,  so the  generators are
|W.1| and |W.2|.

The  character  `φ₁,₄'`  (resp.  `φ₁,₁₂'`,  `φ₂,₃'`) takes the value `1`
(resp.  |E(3)|,  |-E(3)|)  on  |W.1|.  The  character  `φ₁,₈″`  is  complex
conjugate  to `φ₁,₁₆`,  and the  character `φ₁,₈'`  is complex conjugate to
`φ₁,₄'`  . The  character `φ₂,₅″`  is complex  conjugate to `φ₂,₁`; `φ₂,₅'`
take the value `-1` on |W.1|. The character `φ₂,₇'` is complex conjugate to
`φ₂,₅'`.

-  For `G₇`:
We  let |W:=ComplexReflectionGroup(7)|,  so the  generators are
|W.1|, |W.2| and |W.3|.

The  characters  `φ₁,₄'`  and  `φ₁,₁₀'`  take  the  value `1` on |W.2|. The
character  `φ₁,₈″` is complex  conjugate to `φ₁,₁₆`  and `φ₁,₈'` is complex
conjugate  to `φ₁,₄'`. The characters `φ₁,₁₂'`  and `φ₁,₁₈'` take the value
|E(3)| on |W.2|. The character `φ₁,₁₄″` is complex conjugate to `φ₁,₂₂` and
`φ₁,₁₄'`  is complex conjugate to  `φ₁,₁₀'`. The character `φ₂,₃'` takes
the  value |-E(3)| on |W.2| and `φ₂,₁₃'`  takes the value `-1` on |W.2|.
The  characters `φ₂,₁₁″`, `φ₂,₅″`, `φ₂,₇‴` and `φ₂,₁` are Galois conjugate,
as  well as the characters `φ₂,₇'`,  `φ₂,₁₃'`, `φ₂,₁₁'` and `φ₂,₅'`. The
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
`G₃₂`  and `G₃₃` with that of [@Mal00]  and for `G₈` with that described in
[@MR03].

We  should  emphasize  that  for  all  groups  with  a  few exceptions, the
parameters  for characters do  not depend on  any non-canonical choice. The
exceptions  are `G(de,e,n)` with `e>1`, and `G₅`, `G₇`, `G₂₇`, `G₂₈`, `G₂₉`
and  `G₃₄`, groups  which admit  outer automorphisms  preserving the set of
reflections,  so choices  of a  particular value  on a particular generator
must be made for characters which are not invariant by these automorphisms.

{Labels  for  the  classes.}  For  the  exceptional  complex
reflection  groups, the labels for  the classes represent the decomposition
of  a representative  of the  class as  a product  of generators,  with the
additional  conventions that 'z' represents the generator of the center and
for  well-generated groups 'c'  represents a Coxeter  element (a product of
the  generators  which  is  a  regular  element  for the highest reflection
degree).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'CharNames( `W` [,<options>] )'

returns  the  list  of  character  names  for the reflection group `W`. The
optional  <options> is a record which can give alternative names in certain
cases, or a different formatting of names in general.

|    gap> W:=coxgroup("G",2);
    coxgroup("G",2)
    gap> CharNames(W);
    [ "phi{1,0}", "phi{1,6}", "phi{1,3}'", "phi{1,3}''", "phi{2,1}",
      "phi{2,2}" ]
    gap> CharNames(W,rec(TeX:=true));
    [ "φ₁,₀", "φ₁,₆", "φ₁,₃'", "φ₁,₃''",
      "φ₂,₁", "φ₂,₂" ]
    gap> CharNames(W,rec(spaltenstein:=true));
    [ "1", "eps", "epsl", "epsc", "theta'", "theta''" ]
    gap> CharNames(W,rec(spaltenstein:=true,TeX:=true));
    [ "1", "\\varepsilon", "\\varepsilon_l", "\\varepsilon_c",
      "\\theta'", "\\theta''" ]|

The  last two  commands show  the character  names used by Spaltenstein and
Lusztig when describing the Springer correspondence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'CharParams( `W` )'

this  function returns the list of parameters for irreducible characters of
`W`:  partitions for type 'A', double partitions for type 'B', etc`…`
as  described in the introduction. For exceptional groups they are pairs or
triples,  beginning with the  dimension, the valuation  of the fake degree,
and  an ordinal  number if  more than  one character  shares the  first two
invariants.

|    gap> CharParams(coxgroup("G",2));
    [ [ [ 1, 0 ] ], [ [ 1, 6 ] ], [ [ 1, 3, 1 ] ], [ [ 1, 3, 2 ] ],
      [ [ 2, 1 ] ], [ [ 2, 2 ] ] ]|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'LowestPowerfakeDegrees( `W` )'

return a list holding the `b`-function for all irreducible characters of
`W`,  that is,  for each  character `χ`,  the valuation  of the  fake
degree of `χ`. The ordering of the result corresponds to the ordering
of  the characters  in 'CharTable(W)'.  The advantage  of this  function
compared to calling  'fakeDegrees' is that one does not  have to provide
an indeterminate,  and that it  may be much  faster to compute  than the
fake degrees.

|    gap> LowestPowerfakeDegrees( coxgroup( :D, 4 ) );
    [ 6, 6, 7, 12, 4, 3, 6, 2, 2, 4, 1, 2, 0 ]|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'HighestPowerfakeDegrees( `W` )'

returns a list  holding the `B`-function for  all irreducible characters
of  `W`, that  is, for  each character  `χ`, the  degree of  the fake
degree of `χ`. The ordering of the result corresponds to the ordering
of  the characters  in 'CharTable(W)'.  The advantage  of this  function
compared to calling  'fakeDegrees' is that one does not  have to provide
an indeterminate,  and that it  may be much  faster to compute  than the
fake degrees.

|    gap> HighestPowerfakeDegrees( coxgroup( :D, 4 ) );
    [ 10, 10, 11, 12, 8, 9, 10, 6, 6, 8, 5, 6, 0 ]|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'Representations( `W`[, `l`])'

returns  a  list  holding,  for  each  irreducible character of the complex
reflection  group  `W`,  a  list  of  matrices  images  of  the  generating
reflections  of `W`  in a  model of  the corresponding representation. This
function  is based on the classification,  and is not yet fully implemented
for  `G₃₄`;  88  representations  are  missing  out  of  169,  that is 4
representations of dim. 105, 3 of dim. 315, 6 of dim. 420, 4 of dim.840 and
those  of dim. 120, 140, 189, 280, 384,  504, 540, 560, 630, 720, 729, 756,
896, 945, 1260 and 1280.

If  there is a  second argument, it  can be a  list of indices (or a single
integer) and only the representations with these indices (or that index) in
the list of all representations are returned.

|    gap> Representations(coxgroup("B",2));
    [ [ [ [ 1 ] ], [ [ -1 ] ] ],
      [ [ [ 1, 0 ], [ -1, -1 ] ], [ [ 1, 2 ], [ 0, -1 ] ] ],
      [ [ [ -1 ] ], [ [ -1 ] ] ], [ [ [ 1 ] ], [ [ 1 ] ] ],
      [ [ [ -1 ] ], [ [ 1 ] ] ] ]
    gap> Representations(ComplexReflectionGroup(4),7);
    [ [ [ E(3)^2, 0, 0 ], [ 2*E(3)^2, E(3), 0 ], [ E(3), 1, 1 ] ],
      [ [ 1, -1, E(3) ], [ 0, E(3), -2*E(3)^2 ], [ 0, 0, E(3)^2 ] ] ]|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'LowestPowerGenericDegrees( `W` )'

returns  a list holding the `a`-function  for all irreducible characters of
the  Coxeter  group  or  Spetsial  reflection  group `W`, that is, for each
character  `χ`, the  valuation of  the generic  degree of `χ` (in the
one-parameter  Hecke  algebra  'Hecke(W,X(Cyclotomics))'  corresponding  to
`W`).  The  ordering  of  the  result  corresponds  to  the ordering of the
characters in 'CharTable(W)'.

|    gap> LowestPowerGenericDegrees( coxgroup( :D, 4 ) );
    [ 6, 6, 7, 12, 3, 3, 6, 2, 2, 3, 1, 2, 0 ]|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'HighestPowerGenericDegrees( `W` )'

returns  a list holding the `A`-function  for all irreducible characters of
the  Coxeter  group  or  Spetsial  reflection  group `W`, that is, for each
character  `χ`,  the  degree  of  the  generic  degree of `χ` (in the
one-parameter  Hecke  algebra  'Hecke(W,X(Cyclotomics))'  corresponding  to
`W`).  The  ordering  of  the  result  corresponds  to  the ordering of the
characters in 'CharTable(W)'.

|    gap> HighestPowerGenericDegrees( coxgroup( :D, 4 ) );
    [ 10, 10, 11, 12, 9, 9, 10, 6, 6, 9, 5, 6, 0 ]|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'PositionDet( `W` )'

return  the position of the determinant character in the character table of
the group `W` (for Coxeter groups this is the sign character).

|    gap> W := coxgroup( :D, 4 );;
    gap> PositionDet( W );
    4|

See also 'ChevieCharInfo' ("ChevieCharInfo").

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'DetPerm( `W` )'

return  the permutation of the characters of the reflection group `W` which
is effected when tensoring by the determinant character (for Coxeter groups
this is the sign character).

|    gap> W := coxgroup( :D, 4 );;
    gap> DetPerm( W );
    [ 8, 9, 11, 13, 5, 6, 12, 1, 2, 10, 3, 7, 4 ]|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
module Chars

using Gapjm

export charinfo, classinfo, fakedegrees, CharTable

fakedegree(t::TypeIrred,p,q)=getchev(t,:FakeDegree,p,q)

"""
`fakeDegree(W, φ, q)`

returns  the  fake degree  of  the  character  of parameter  φ  (see
:CharParams)  of  the  reflection  group `W`,  evaluated  at  `q`  (see
"fakeDegrees" for a definition of the fake degrees).

```julia-repl
julia> Chars.fakedegree(coxgroup(:A,2),[[2,1]],Pol(:q))
q²+q
```
"""
function fakedegree(W,p,q)
  typ=refltype(W)
  if isempty(typ) return one(q) end
  prod(map((t,p)->fakedegree(t,p,q),typ,p))
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
3-element Array{Pol{Int64},1}:
 q³  
 q²+q
 1   
```
"""
function fakedegrees(W,q)
  map(p->fakedegree(W,p,q),charinfo(W)[:charparams])
end

function charinfo(t::TypeIrred)
  c=deepcopy(getchev(t,:CharInfo))
  c[:positionId]=c[:extRefl][1]
  c[:positionDet]=c[:extRefl][end]
  c[:charnames]=map(c[:charparams]) do p
     getchev(t,:CharName,p,Dict(:TeX=>true))
  end
  c
end

cartfields(p,f)=Cartesian(getindex.(p,f)...)

"""
'ChevieCharInfo( `W` )'

returns  information  about the  irreducible  characters  of the  finite
reflection group `W`.  The result is a record with the following
components:

'charparams': contains parameters for the irreducible characters as
     described  in the  introduction or  returned by 'CharParams(`W`)'. The
     parameters   are  tuples  with  one  component  for  each  irreducible
     irreducible  component of `W`  (as given by  'ReflectionType'). For an
     irreducible  component which  is an  imprimitive reflection  group the
     component  of  the  'charparam'  is  a  tuple of partitions, and for a
     primitive  irreducible group  it is  a pair  `(d,e)` where  `d` is the
     degree of the character and `e` is the smallest symmetric power of the
     character  of the  reflection representation  which contains the given
     character as a component.

'charnames': strings describing the irreducible characters, computed from
     the 'charparams'. This is the same as 'CharNames(`W`)'.

'positionId': the  position of the  trivial character  in the  character
     table  of `W`  (which is  also returned by the function 'PositionId').

'positionDet':  Contains the position of the determinant character in the
   character  table  of  `W`  (which  is  also  returned  by  the  function
   'PositionDet'). For Coxeter groups this is the sign character.

'extRefl':  Only  present  if  `W`  is  irreducible,  in  which  case the
     reflection  representation  of  `W`  and  all  its exterior powers are
     irreducible.  It then contains the position  of the exterior powers of
     the reflection representation in the character table.

'b': contains  the    result   of    'LowestPowerFakeDegrees(`W`)'.

'B': contains  the    result   of    'HighestPowerFakeDegrees(`W`)'.

'a': Only filled for Spetsial groups. Contains the result of
   'LowestPowerGenericDegrees(`W`)'.

'A': Only filled for Spetsial groups. Contains the result of
   'HighestPowerGenericDegrees(`W`)'.

'opdam': Contains the permutation of the characters obtained by composing
   the  Opdam involution with complex  conjugation. This permutation has an
   interpretation  as a Galois action on the characters of 'H:=Hecke(W,x)'
   where  'x:=Indeterminate(Cyclotomics)': if 'H' splits by taking `v` an
   `e`-th  root of  `x`, '.opdam'  records the  permutation effected by the
   Galois action 'v->E(e)*v'.

```julia-repl
julia> charinfo(ComplexReflectionGroup(22))
Dict{Symbol,Any} with 8 entries:
  :positionDet => 2
  :b           => [0, 30, 11, 13, 1, 7, 2, 6, 12, 16, 3, 6, 9, 8, 4, 10, 7, 5]
  :positionId  => 1
  :charnames   => ["\\phi_{1,0}", "\\phi_{1,30}", "\\phi_{2,11}", "\\phi_{2,13}…
  :extRefl     => [1, 5, 2]
  :indexchars  => [1, 16, 65, 32, 35, 62, 91, 92, 101, 102, 151, 166, 181, 196,…
  :charparams  => Array{Array{Int64,1},1}[[[1, 0]], [[1, 30]], [[2, 11]], [[2, …
  :opdam       => (3,5)(4,6)(11,13)(12,14)(17,18)

julia> charinfo(coxgroup(:G,2))
Dict{Symbol,Any} with 10 entries:
  :a            => [0, 6, 1, 1, 1, 1]
  :b            => [0, 6, 3, 3, 1, 2]
  :positionId   => 1
  :charnames    => ["\\phi_{1,0}", "\\phi_{1,6}", "\\phi_{1,3}'", "\\phi_{1,3}'…
  :A            => [0, 6, 5, 5, 5, 5]
  :B            => [0, 6, 3, 3, 5, 4]
  :spaltenstein => ["1", "\\varepsilon", "\\varepsilon_l", "\\varepsilon_c",…
  :extRefl      => [1, 5, 2]
  :charparams   => Array{Array{Int64,1},1}[[[1, 0]], [[1, 6]], [[1, 3, 1]], [[1…
  :positionDet  => 2
```

For  irreducible groups, the returned  record contains sometimes additional
information:

for  `F₄`: the field 'kondo' gives  the labeling of the characters given
  by Kondo, also used in [@Lus85, (4.10)].

for  `E₆,  E₇,  E₈`:  the  field  'frame'  gives  the  labeling of the
  characters  given  by  Frame,  also  used  in  [@Lus85, (4.11),  (4.12), and
  (4.13)].

for  `G₂`: the field 'spaltenstein' gives the labeling of the characters
  given by Spaltenstein.

for  `G(de,e,2)`  even  `e`  and  `d>1`:  the  field  'malle'  gives  the
parameters for the characters used by Malle in [@Mal96].
"""
function charinfo(W)::Dict{Symbol,Any}
  gets(W,:charinfo) do W
    p=charinfo.(refltype(W))
    if isempty(p) return Dict(:a=>[0],:A=>[0],:b=>[0],:B=>[0],:positionId=>1,
      :positionDet=>1,:charnames=>["Id"],:extRefl=>[1],:charparams=>[[]])
    end
    if length(p)==1 res=copy(p[1]) else res=Dict{Symbol, Any}() end
    res[:charparams]=cartfields(p,:charparams)
    if length(p)==1 return res end
    res[:charnames]=map(l->join(l,","),cartfields(p,:charnames))
    for f in [:positionId, :positionDet]
     if all(d->haskey(d,f),p)
       res[f]=PositionCartesian(map(x->length(x[:charparams]),p),getindex.(p,f))
      end
    end
    for f in [:b, :B, :a, :A]
      if all(d->haskey(d,f),p) res[f]=Int.(map(sum,cartfields(p,f))) end
    end
    if any(x->haskey(x, :opdam),p)
      res[:opdam]=map(x->haskey(x,:opdam) ? x[:opdam] : Perm(), p)
      gt=Cartesian(map(x->1:length(x[:charparams]), p))
      res[:opdam]=PermListList(gt, map(t->map((x,i)->x^i,t,res[:opdam]),gt))
    end
    res
  end
end

"""
'ChevieClassInfo( `W` )'

returns  information about the  conjugacy classes of  the finite reflection
group `W`. The result is a record with three components:

'classtext': contains words in the generators  describing representatives
     of  each conjugacy class.  Each word is  a list of  integers where the
     generator  `s_i` is represented by the integer `i`. For finite Coxeter
     groups, it is the same as
     'List(ConjugacyClasses(W),x->CoxeterWord(W,Representative(x)))',   and
     each  such representative is of minimal  length in its conjugacy class
     and is a "very good" element in the sense of [@GM97].

'classparams':   The elements of this list are tuples which have one
     component  for each irreducible component of `W`. These components for
     the infinite series, contain partitions or partition tuples describing
     the  class (see the introduction).  For the exceptional Coxeter groups
     they  contain  Carter's  admissible  diagrams,  see [@Car72]. For
     exceptional complex reflection groups they contain in general the same
     information as in classtext.

'classnames': Contains strings describing  the conjugacy classes, made
     out of the information in 'classparams'.

```julia-repl
julia> classinfo(coxgroup(:D,4))
Dict{Symbol,Any} with 6 entries:
  :classes      => [1, 6, 1, 12, 24, 12, 6, 6, 12, 32, 32, 24, 24]
  :orders       => [1, 2, 2, 2, 4, 2, 2, 2, 4, 3, 6, 4, 4]
  :classtext    => Array{Int64,1}[[], [1, 2], [1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 3,…
  :centralizers => Any[192//1, 32//1, 192//1, 16//1, 8//1, 16//1, 32//1, 32//1,…
  :classnames   => ["1111.", "11.11", ".1111", "211.", "1.21", "2.11", "22.+", …
  :classparams  => Array{T,1} where T[Array{Int64,1}[[1, 1, 1, 1], []], Array{I…
```

See also the introduction of this section.
"""
function classinfo(t::TypeIrred)
  cl=deepcopy(getchev(t,:ClassInfo))
  inds=t[:indices]
  cl[:classtext]=map(x->inds[x],cl[:classtext])
  cl[:classes]=Int.(cl[:classes])
  cl
end

function classinfo(W)::Dict{Symbol,Any}
  gets(W,:classinfo) do W
    tmp=map(classinfo,refltype(W))
    if isempty(tmp) return Dict(:classtext=>[Int[]],:classnames=>[""],
                      :classparams=>[Int[]],:orders=>[1],:classes=>[1])
    end
    if any(isnothing, tmp) return nothing end
    if length(tmp)==1 res=copy(tmp[1]) else res=Dict{Symbol, Any}() end
    res[:classtext]=map(x->vcat(x...),cartfields(tmp,:classtext))
    res[:classnames]=map(join,cartfields(tmp,:classnames))
    if all(x->haskey(x,:classparam),tmp)
      res[:classparams]=cartfields(tmp,:classparams)
    end
    if all(x->haskey(x,:orders),tmp)
      res[:orders]=map(lcm, cartfields(tmp,:orders))
    end
    if all(x->haskey(x,:classes),tmp)
      res[:classes]=map(prod, cartfields(tmp,:classes))
    end
    res
  end
end

#--------------- CharTables -----------------------------------------
"""
 CharTable is a structure to hold character tables of groups and Hecke
 algebras
"""
struct CharTable{T}
  irr::Matrix{T}
  charnames::Vector{String}
  classnames::Vector{String}
  centralizers::Vector{Int}
  identifier::String
end

function Base.show(io::IO,ct::CharTable)
  println(io,"CharTable(",ct.identifier,")")
  irr=map(ct.irr)do e
   if iszero(e) "." else sprint(show,e; context=io) end
  end
  format(io,irr,row_labels=TeXstrip.(ct.charnames),
                column_labels=TeXstrip.(ct.classnames))
end

function CharTable(t::TypeIrred)
  ct=getchev(t,:CharTable)
  if haskey(ct,:irredinfo) names=getindex.(ct[:irredinfo],:charname)
  else                     names=charinfo(t)[:charnames]
  end
  if !haskey(ct,:classnames) merge!(ct,classinfo(t)) end
  irr=toM(ct[:irreducibles])
  if eltype(irr)==Rational{Int} irr=Int.(irr)
  elseif eltype(irr)==Cyc{Rational{Int}} irr=Cyc{Int}.(irr)
  end
  CharTable(irr,names,ct[:classnames],Int.(ct[:centralizers]),ct[:identifier])
end

function CharTable(W)::CharTable
  gets(W,:chartable) do W
    ctt=CharTable.(refltype(W))
    if isempty(ctt) 
      return CharTable(hcat(1),["Id"],["1"],[1],"$W")
    end
    charnames=join.(Cartesian(getfield.(ctt,:charnames)...),",")
    classnames=join.(Cartesian(getfield.(ctt,:classnames)...),",")
    centralizers=prod.(Cartesian(getfield.(ctt,:centralizers)...))
    identifier=join(getfield.(ctt,:identifier),"×")
    if length(ctt)==1 irr=ctt[1].irr 
    else irr=kron(getfield.(ctt,:irr)...)
    end
    CharTable(irr,charnames,classnames,centralizers,identifier)
  end
end

function CharTable(H::HeckeAlgebra{C})where C
  W=H.W
  ct=impl1(getchev(W,:HeckeCharTable,H.para,
       haskey(H.prop,:rootpara) ? rootpara(H) : fill(nothing,length(H.para))))
  if haskey(ct,:irredinfo) names=getindex.(ct[:irredinfo],:charname)
  else                     names=charinfo(W)[:charnames]
  end
  CharTable(Matrix(convert.(C,toM(ct[:irreducibles]))),names,
     ct[:classnames],map(Int,ct[:centralizers]),ct[:identifier])
end
end
