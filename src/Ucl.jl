"""
This  module gives information  about the unipotent  conjugacy classes of a
connected  reductive group  `𝐆` over  an algebraically  closed field `k` of
characteristic   `p`,  and  various  associated  invariants,  such  as  the
generalised  Springer correspondence. The unipotent  classes depend on `p`;
their  classification differs  when `p`  is not  *good* (that  is, when `p`
divides  one  of  the  coefficients  of  the  highest  root  of  one of the
irreducible  components). In good characteristic, the unipotent classes are
in bijection with nilpotent orbits on the Lie algebra `𝔤` of `𝐆`.

For each unipotent conjugacy class, we give the following information for a
representative `u` of that class:

  - the centralizer ``C_𝐆 (u)``, that we describe by the reductive part of
    ``C_𝐆  (u)^0``, by the  group of components  ``A(u):=C_𝐆 (u)/C_𝐆 (u)^0``,
    and by the dimension of its unipotent radical.

  - in good characteristic, the  Dynkin-Richardson  diagram.

  - the generalised  Springer correspondence,  attaching a  character of the
    Weyl group or a relative Weyl group to each character of `A(u)`.

The  Dynkin-Richarson diagram is  attached to a  nilpotent element ``e`` of
the  Lie algebra  ``𝔤``. By  the Jacobson-Morozov  theorem there  exists an
``𝔰𝔩₂`` subalgebra of ``𝔤`` containing ``e`` as the element
``\\begin{pmatrix}0&1\\\\0&0  \\end{pmatrix}``.  Let  ``𝐒``  be  the  torus
``\\begin{pmatrix}h&0\\\\0&h^{-1} \\end{pmatrix}`` of ``SL₂`` and let ``𝐓``
be  a  maximal  torus  containing  ``𝐒``,  so  that ``𝐒`` is the image of a
one-parameter   subgroup  ``σ∈  Y(𝐓)``.  Consider  the  root  decomposition
``𝔤=∑_{α∈Φ}𝔤_α``  given by ``𝐓`` and the  root system `Φ`; then ``α↦⟨σ,α⟩``
defines  a linear form  on ``Φ``, determined  by its value  on simple roots
`Π`.  It is possible to choose a  system of simple roots such that ``⟨σ,α⟩≥
0``  for ``α∈Π``,  and then  ``⟨σ,α⟩∈{0,1,2}`` for  any ``α∈Π``. The Dynkin
diagram  of  ``Π``  decorated  by  these  values  ``0,1,2``  is  called the
Dynkin-Richardson  diagram  of  ``e``,  and  in  good  characteristic  is a
complete  invariant  of  the  ``𝔤``-orbit  of  ``e``. The Dynkin-Richardson
diagrams thus classify unipotent classes of `𝐆` in good characteristic.

Another classification of unipotent classes was given by Bala and Carter. A
standard  parabolic subgroup `𝐏` of `𝐆` associated with the subset `I⊂Π` of
the simple roots is *distinguished* if the linear form `σ` taking the value
`2` on `α∈ I` and `0` on other simple roots satisfies
`2n₀+semisimplerank(𝐆)=n₂`,  where `nᵢ` is the number of roots in `Φ` where
`σ`  takes the value `i`.  Given a distinguished parabolic  `𝐏`, there is a
unique  unipotent class which is dense in the unipotent radical of `𝐏`. For
this  class, the linear form described  by the Dynkin-Richardson diagram is
equal  to  `σ`.  Such  unipotent  classes  are  called *distinguished*. The
Bala-Carter theorem says that every unipotent class is distinguished in the
smallest  Levi subgroup `𝐋` which  contains it, and that  such pairs of `𝐋`
and the distinguished parabolic `𝐏` of `𝐋` taken up to `𝐆`-conjugacy are in
bijection with unipotent classes of `𝐆`. The function `induced_linear_form`
allows to find the Dynkin-Richardson diagram in `𝐆` of a unipotent class of
`𝐋` and thus recover these diagrams from the Bala-Carter classification.

Let  ``ℬ`` be  the variety  of all  Borel subgroups  and let  ``ℬᵤ`` be the
subvariety  of Borel subgroups  containing the unipotent  element `u`. Then
``dim C_𝐆(u)=rank 𝐆 + 2 dim ℬ_u`` and in good characteristic this dimension
can  be computed from the linear  form `σ` defined by the Dynkin-Richardson
diagram:  the dimension of the class of `u` is the number of roots `α` such
that ``⟨σ,α⟩∉{0,1}``.

We   now  describe  the  Springer  correspondence.  Indecomposable  locally
constant  ``𝐆``-equivariant  sheaves  on  a  unipotent  class ``C``, called
*local  systems*, are  parameterised by  irreducible characters of ``A(u)``
for  `u∈ C`. The *ordinary* Springer  correspondence is a bijection between
irreducible  characters of the Weyl  group and a large  subset of the local
systems  containing all trivial  local systems (those  parameterised by the
trivial  character  of  ``A(u)``  for  each  ``u``).  More  generally,  the
*generalised*  Springer correspondence  associates to  each local  system a
(unique  up to ``𝐆``-conjugacy) *cuspidal datum*,  a Levi subgroup ``𝐋`` of
``𝐆``  and a *cuspidal* local system on  an unipotent class of ``𝐋``, and a
character  of the relative Weyl group ``W_𝐆 (𝐋):=N_𝐆 (𝐋)/𝐋``. There are few
cuspidal  local  systems  (at  most  one  in  each  dimension for classical
groups).  The ordinary  Springer correspondence  is the  special case where
``𝐋``  is a maximal torus, the cuspidal  local system is the trivial system
on the identity element, and ``W_𝐆 (𝐋)=W``, the Weyl group.

The  Springer correspondence gives information on the character values of a
finite  reductive  group  as  follows:  assume  that ``k`` is the algebraic
closure  of a finite field ``𝔽_q`` and that ``F`` is the Frobenius attached
to  an ``𝔽_q``-structure of  ``𝐆``. Let ``C``  be an ``F``-stable unipotent
class  and let ``u∈ C^F``; we call ``C`` the *geometric class* of ``u`` and
the ``𝐆^F``-classes within ``C^F`` are parameterised by the ``F``-conjugacy
classes  of ``A(u)``,  denoted by  ``H¹(F,A(u))`` (most  of the time we can
find ``u∈ C`` such that ``F`` acts trivially on ``A(u)`` and ``H¹(F,A(u))``
parametrises  the conjugacy classes). To an ``F``-stable character ``φ`` of
``A(u)``  we associate  the *characteristic  function* of the corresponding
local  system (actually  associated with  an extension  ``φ̃`` of  ``φ`` to
``A(u).F``);  it is  a class  function ``Y_{u,φ}``  on ``𝐆^F`` which can be
normalised  so  that:  ``Y_{u,φ}(u₁)=φ̃(cF)``  if  ``u₁``  is geometrically
conjugate   to  ``u``  and  its   ``𝐆^F``-class  is  parameterised  by  the
``F``-conjugacy  class ``cF`` of  ``A(u)``, otherwise ``Y_{u,φ}(u₁)=0``. If
the  pair  ``u,φ``  corresponds  via  the  Springer  correspondence  to the
character  ``χ`` of  ``W_𝐆(𝐋)``, then  ``Y_{u,φ}`` is  also denoted ``Yᵪ``.
There  is another important class of functions indexed by local systems: to
a  local  system  on  class  ``C``  is  attached an intersection cohomology
complex,  which is  a complex  of sheaves  supported on the closure ``C̄``.
Such a complex of sheaves is associated with its *characteristic function*,
a class function of ``𝐆^F`` obtained by taking the alternating trace of the
Frobenius acting on the stalks of the cohomology sheaves. If ``Y_ψ`` is the
characteristic  function of a local  system, the characteristic function of
the  corresponding intersection  cohomology complex  is denoted by ``X_ψ``.
This  function is supported on ``C̄``, and Lusztig has shown that ``X_ψ=∑_χ
P_{ψ,χ}  Yᵪ`` where ``P_{ψ,χ}`` are integer polynomials in ``q`` and ``Yᵪ``
are attached to local systems on classes lying in ``C̄``.

Lusztig   and  Shoji  have  given  an   algorithm  to  compute  the  matrix
``P_{ψ,χ}``, which is implemented in Chevie. The relation to the characters
of   ``𝐆(𝔽_q)``,   considering   for   simplicity   the  ordinary  Springer
correspondence,  is that the  restriction to the  unipotent elements of the
almost  character ``R_χ`` is equal to  ``q^{bᵪ} Xᵪ``, where ``bᵪ`` is ``dim
ℬᵤ``  for an element `u` of  the class `C` such that  the support of `χ` is
``C̄``.  The restrictions of the Deligne-Lusztig characters ``R_w`` for `w∈
W`  on the  unipotents are  called the  *Green functions*  and can  also be
computed  by Chevie.  The values  of all  unipotent characters on unipotent
elements  can also be  computed in principle  by applying Lusztig's Fourier
transform  matrix (see the  section on the  Fourier matrix) but  there is a
difficulty  in that the  ``Xᵪ`` must first  be multiplied by  some roots of
unity  which are not known  in all cases (and  when known may depend on the
congruence class of ``q`` modulo some small prime).

Finally,  we describe  how unipotent  classes of  `𝐆` are  parameterised in
various   quasisimple  groups.   In  classical   types,  the   classes  are
parametrised by partitions corresponding to the Jordan form in the standard
representation. Thus,

  - for `Aₙ` we have partitions of `n+1`.
  - for  `B_n` we have partitions of `2n+1`  where even parts occur an even
    number  of times. In characteristic  2, types B and  C are isogenous so
    have the same classification; thus see the next paragraph.
  - for  `C_n` we  have partitions  of `2n`  where odd  parts occur an even
    number  of times. In characteristic 2,  there are `2ᵏ` classes attached
    to  a partition where  `k` is the  number of even  parts which occur an
    even number of times.
  - for  `D_n` we have partitions of `2n`  where even parts occur an even
    number  of times,  excepted there  are two  classes when  all parts are
    even.  In characteristic 2, we have  partitions of `2n` where odd parts
    occur  an  even  number  of  times,  excepted  there are `2ᵏ+δ` classes
    attached  to a partition  where `k` is  the number of  even parts which
    occur an even number of times, and `δ` is 2 when all parts are even and
    0 otherwise.
In  exceptional  groups,  the  names  of  the  classes are derived from the
Bala-Carter  classification. The name of a class parametrised by `(𝐋,𝐏)` is
of  the form `l(p)`  where `l` is  the name of  `𝐋` and `(p)` is present if
there  is more than one distinguished  parabolic in `𝐋` and describes which
one  it is. Before the  Bala-Carter classification was universally adopted,
Shoji  and  Mizuno  used  a  different  scheme  where a class was sometimes
parametrised  by a reductive subgroup of maximal rank which was not a Levi.
These  older labels  can be  obtained instead  by giving  the `IO` property
`:shoji=>true`  or `:mizuno=>true`. In a  bad characteristic `p`, there are
extra  classes.  Each  of  them  is  associated  with  a  class `c` in good
characteristic and is called `(c)ₚ`.

We illustrate the above descriptions on some examples:

```julia-repl
julia> UnipotentClasses(rootdatum(:sl,4))
UnipotentClasses(sl₄)
1111<211<22<31<4
┌────┬───────────────────────────────────────────────────────────────┐
│u   │D-R dℬ ᵤ B-C   C_𝐆(u) A₃(Φ₁³) A₁(A₁×A₁Φ₁)/-1 .(A₃)/ζ₄ .(A₃)/ζ₄³│
├────┼───────────────────────────────────────────────────────────────┤
│4   │222    0 222    q³.Z₄     1:4           -1:2    ζ₄:Id    ζ₄³:Id│
│31  │202    1 22.    q⁴.Φ₁   Id:31                                  │
│22  │020    2 2.2 q⁴.A₁.Z₂    2:22          11:11                   │
│211 │101    3 2..  q⁵.A₁Φ₁  Id:211                                  │
│1111│000    6 ...       A₃ Id:1111                                  │
└────┴───────────────────────────────────────────────────────────────┘
```
The first column of the table gives the name of the unipotent class, here a
partition  describing  the  Jordan  form.  The  partial  order on unipotent
classes  given by the Zariski closure is given before the table. The column
`D-R`,   which   is   only   shown   in   good  characteristic,  gives  the
Dynkin-Richardson  diagram  for  each  class;  the  column  `dBu` gives the
dimension  of the variety  ``ℬ ᵤ``. The  column `B-C` gives the Bala-Carter
classification of ``u``, that is in the case of ``sl₄`` it shows ``u`` as a
regular  unipotent  in  a  Levi  subgroup  by  giving the Dynkin-Richardson
diagram  of a regular unipotent (all  2's) for the entries corresponding to
the  Levi and `.` for the entries  which not corresponding to the Levi. The
column `C(u)` describes the group ``C_𝐆(u)``: a power ``qᵈ`` describes that
the  unipotent  radical  of  ``C_𝐆(u)``  has  dimension  ``d`` (thus ``qᵈ``
rational  points); then follows a description  of the reductive part of the
neutral  component of ``C_𝐆(u)``, given by the name of its root datum. Then
if  ``C_𝐆(u)`` is not connected, the description of ``A(u)`` is given using
a  different vocabulary: a cyclic group of order  4 is given as `Z4`, and a
symmetric group on 3 points would be given as `S3`.

For  example, the first  class `4` has  ``C_𝐆(u)^0`` unipotent of dimension
`3` and ``A(u)`` equal to `Z4`, the cyclic group of order 4. The class `22`
has  ``C_G(u)`` with unipotent radical of  dimension `4`, reductive part of
type  `A1` and  ``A(u)`` is  `Z2`, the  cyclic group  of order 2. The other
classes have ``C_𝐆(u)`` connected. For the class `31` the reductive part of
``C_G(u)`` is a torus of rank 1.

Then  there is a column  for each *Springer series*,  giving for each class
the  pairs  'a:b'  where  'a'  is  the  name  of  the character of ``A(u)``
describing  the local system involved and 'b'  is the name of the character
of  the (relative) Weyl group corresponding by the Springer correspondence.
At  the top of  the column is  the name of  the relative Weyl group, and in
brackets  the name  of the  Levi affording  a cuspidal  local system; then,
separated  by a `/`,  is a description  of the central character associated
with  the Springer series  (omitted if this  central character is trivial):
all  local systems in a given Springer  series have the same restriction to
the  centre of ``𝐆``. To see what the picture becomes for another algebraic
group  in the same isogeny class, for example the adjoint group, one simply
discards the Springer series whose central character becomes trivial on the
centre  of ``𝐆``; and each group ``A(u)``  must be quotiented by the common
kernel  of  the  remaining  characters.  Here  is the table for the adjoint
group:

```julia-repl
julia> UnipotentClasses(coxgroup(:A,3))
UnipotentClasses(A₃)
1111<211<22<31<4
┌────┬────────────────────────────┐
│u   │D-R dℬ ᵤ B-C  C_𝐆(u) A₃(Φ₁³)│
├────┼────────────────────────────┤
│4   │222    0 222      q³    Id:4│
│31  │202    1 22.   q⁴.Φ₁   Id:31│
│22  │020    2 2.2   q⁴.A₁   Id:22│
│211 │101    3 2.. q⁵.A₁Φ₁  Id:211│
│1111│000    6 ...      A₃ Id:1111│
└────┴────────────────────────────┘
```
Here is another example:

```julia-repl
julia> UnipotentClasses(coxgroup(:G,2))
UnipotentClasses(G₂)
1<A₁<Ã₁<G₂(a₁)<G₂
┌──────┬──────────────────────────────────────────┐
│u     │D-R dℬ ᵤ B-C C_𝐆(u)         G₂(Φ₁²)  .(G₂)│
├──────┼──────────────────────────────────────────┤
│G₂    │ 22    0  22     q²         Id:φ₁‚₀       │
│G₂(a₁)│ 20    1  20  q⁴.S₃ 21:φ′₁‚₃ 3:φ₂‚₁ 111:Id│
│Ã₁    │ 01    2  .2  q³.A₁         Id:φ₂‚₂       │
│A₁    │ 10    3  2.  q⁵.A₁        Id:φ″₁‚₃       │
│1     │ 00    6  ..     G₂         Id:φ₁‚₆       │
└──────┴──────────────────────────────────────────┘
```
which illustrates that on class `G₂(a₁)` there are two local systems in the
principal series of the Springer correspondence, and another cuspidal local
system. It also illustrates how we represent the Bala-Carter classification
in  general. If a class is attached to `(𝐋,𝐏)` then the simple roots in the
complement  of `𝐋` have a `.`.  Those in `𝐋` have a  `0` or a `2`, the `2`s
characterising `𝐏`. So, from the `B-C` column, we see that that `G₂(a₁)` is
not  in a proper Levi, in which case the Bala-Carter diagram coincides with
the Dynkin-Richardson diagram.

Characteristics 2 and 3 are not good for `G2`. To get the unipotent classes
and  the Springer  correspondence in  a bad  characteristic, give  a second
argument to the function `UnipotentClasses`:

```julia-repl
julia> UnipotentClasses(coxgroup(:G,2),3)
UnipotentClasses(G₂,3)
1<A₁,(Ã₁)₃<Ã₁<G₂(a₁)<G₂
┌──────┬───────────────────────────────────────────┐
│u     │dℬ ᵤ B-C C_𝐆(u)  G₂(Φ₁²) .(G₂) .(G₂)  .(G₂)│
├──────┼───────────────────────────────────────────┤
│G₂    │   0  22  q².Z₃   1:φ₁‚₀       ζ₃:Id ζ₃²:Id│
│G₂(a₁)│   1  20  q⁴.Z₂   2:φ₂‚₁ 11:Id             │
│Ã₁    │   2  .2     q⁶  Id:φ₂‚₂                   │
│A₁    │   3  2.  q⁵.A₁ Id:φ″₁‚₃                   │
│(Ã₁)₃ │   3  ??  q⁵.A₁ Id:φ′₁‚₃                   │
│1     │   6  ..     G₂  Id:φ₁‚₆                   │
└──────┴───────────────────────────────────────────┘
```

The  function `ICCTable` gives the  transition matrix between the functions
``Xᵪ``  and ``Y_ψ``.

```julia-repl
julia> uc=UnipotentClasses(coxgroup(:G,2));
julia> t=ICCTable(uc;q=Pol(:q))
Coefficients of Xᵪ on Yᵩ for series L=G₂₍₎=Φ₁² W_G(L)=G₂
┌──────┬─────────────────────────────┐
│X\\Y   │G₂ G₂(a₁)⁽²¹⁾ G₂(a₁) Ã₁ A₁  1│
├──────┼─────────────────────────────┤
│Xφ₁‚₀ │ 1          .      1  1  1  1│
│Xφ′₁‚₃│ .          1      .  1  . q²│
│Xφ₂‚₁ │ .          .      1  1  1 Φ₈│
│Xφ₂‚₂ │ .          .      .  1  1 Φ₄│
│Xφ″₁‚₃│ .          .      .  .  1  1│
│Xφ₁‚₆ │ .          .      .  .  .  1│
└──────┴─────────────────────────────┘
```

An example which illustrates how to get the `shoji` names of classes
```julia-rep1
julia> uc=UnipotentClasses(coxgroup(:F,4));

julia> uc.classes[10:end]
7-element Vector{UnipotentClass}:
 UnipotentClass(C₃(a₁))
 UnipotentClass(F₄(a₃))
 UnipotentClass(C₃)
 UnipotentClass(B₃)
 UnipotentClass(F₄(a₂))
 UnipotentClass(F₄(a₁))
 UnipotentClass(F₄)

julia> xdisplay(uc.classes[10:end],shoji=true)
7-element Vector{UnipotentClass}:
 UnipotentClass(A₁+B₂)
 UnipotentClass(A₃+Ã₁)
 UnipotentClass(C₃)
 UnipotentClass(B₃)
 UnipotentClass(C₃+A₁)
 UnipotentClass(B₄)
 UnipotentClass(F₄)
```

In  the `ICCTable` the  row labels and  column labels show  the two ways of
indexing  local systems: the row labels  give the character of the relative
Weyl  group and the column labels give the  class and the name of the local
system as a character of `A(u)`: for example, `G2(a1)` is the trivial local
system  of the  class `G2(a1)`,  while `G2(a1)(21)`  is the local system on
this class corresponding to the 2-dimensional character of ``A(u)=A₂``.

The  data on unipotent classes for  arbitrary reductive groups are obtained
as follows. We may reduce to the case where the radical is trivial, that is
`𝐆`  is semisimple. Then `𝐆` is the quotient of a simply connected group by
a  finite central group. The data  for quasi-simple simply connected groups
have  been entered  by hand  for each  type and  a general simply connected
group  is  a  direct  product  of  quasi-simple  ones.  In  such a group, a
character  of `A(Z)`, the group of components of the centre, is attached to
each  Springer series.  For an  isogenous reductive  group `𝐆'` with centre
`Z'`  the group  `A(Z')` is  a quotient  of the  group `A(Z)`. The Springer
series  for `𝐆'` are those such  that the corresponding character of `A(Z)`
factors through `A(Z')` (to compute `A(Z')` see
[`algebraic_center`](@ref)).  The  geometric  unipotent  classes of `𝐆` and
`𝐆'`  are in bijection. For  `u` a unipotent element  of `𝐆'` (which we can
also consider as a unipotent element of `𝐆`) the group `A₁=A(u)` in `𝐆'` is
a quotient of `A=A(u)` in `𝐆` which we can compute as follows: the Springer
correspondence  for `𝐆'` tells us which  characters of `A` survive in `𝐆'`.
Then `A'` is the quotient of `A` by the common kernel of these characters.
"""
module Ucl

using ..Chevie

export UnipotentClasses, UnipotentClass, ICCTable, XTable,
 GreenTable, UnipotentValues, induced_linear_form, special_pieces, name,
 distinguished_parabolics, representative

@GapObj mutable struct UnipotentClass
  name::String
  parameter::Any
  dimBu::Int
end

@doc """
A `struct UnipotentClass` representing the class `C` of a unipotent element
`u`  of the reductive  group `𝐆` with  Weyl group `W`,  contains always the
following information
  * `.name`  The name of `C`
  * `.parameter` A parameter describing `C`. Sometimes the same as `.name`; a partition describing the Jordan form, for classical groups.
  * `.dimBu` The dimension of the variety of Borel subgroups containing `u`.

For  some  classes  in  types  `E₆,  E₇,  E₈`  there  is  a field `.mizuno`
containing the names given by Mizuno for some classes, and for some classes
in type `F₄` a field `.shoji` containing the names given by Shoji.

A  `UnipotentClass` contains also some of  the following information (all of
it for some types and some characteristics but sometimes much less)
  * `.dynkin` the Dynkin-Richardson diagram of `C` (a vector giving a weight 0, 1 or 2 to the simple roots).
  *  `.dimred` the dimension of the reductive part of `C_G(u)`.
  *  `.red` a `CoxeterCoset` recording the type of the reductive part of `C_G(u)`, with the twisting induced by the Frobenius if any.
  *  `.Au` the group `A_G(u)=C_G(u)/C^0_G(u)`.
  *  `.balacarter` encodes the Bala-Carter classification of `C`, which says that `u` is distinguished in a Levi `L` (the Richardson class in a parabolic `P_L`) as a vector listing the indices of the simple roots in `L`, with those not in `P_L` negated.
  *  `.rep` a list of indices for roots such that if `U=UnipotentGroup(W)` then `prod(U,u.rep)` is a representative of `C` (which can be obtained also by `representative(W,u)`).
  *  `.dimunip` the dimension of the unipotent part of `C_G(u)`.
  *  `.AuAction` an `ExtendedCoxeterGroup` recording the action of `A_G(u)` on `red`.
""" UnipotentClass

@GapObj struct UnipotentClasses
  classes::Vector{UnipotentClass}
  p::Int
  orderclasses::Poset
  springerseries::Vector{Dict}
end

function nameclass(u;opt...)
# println("u=$u, opt=$opt")
  if haskey(opt,:mizuno) && haskey(u,:mizuno) n=u[:mizuno]
  elseif haskey(opt,:shoji) && haskey(u,:shoji) n=u[:shoji]
  else n=u isa Dict ? u[:name] : u.name
  end
  TeX=haskey(opt,:TeX)
  n=fromTeX(n;opt...)
  if haskey(opt,:locsys) && opt[:locsys]!=charinfo(u[:Au]).positionId
    cl="("*charnames(u[:Au];opt...)[opt[:locsys]]*")"
    if haskey(opt,:noclass_nontrivial) && opt[:noclass_nontrivial]
      n="{}^{$cl}"
    else n*="^{$cl}"
    end
    n=fromTeX(n;opt...)
  elseif haskey(opt,:class) && opt[:class]!=position_class(u[:Au],one(u[:Au]))
    cl=conjugacy_classes(u[:Au])[opt[:class]].name
    n=TeX ? "\\mbox{\$$n\$}_{($cl)}" : fromTeX("$(n)_{$cl}";opt...)
  end
  n
end

name(io::IO,u::UnipotentClass)=nameclass(u;IOContext(io).dict...)

name(u;opt...)=name(IOContext(stdout,opt...),u)

function Base.show(io::IO,u::UnipotentClass)
  print(io,"UnipotentClass(",name(io,u),")")
end

representative(W,u::UnipotentClass)=prod(UnipotentGroup(W),u.rep)

"""
`induced_linear_form(W, H, h)`

This routine can be used to find the unipotent class in the reductive group
with  Weyl group `W` which contains a  given unipotent class of a reductive
subgroup of maximum rank represented by the reflection subgroup `H` of `W`.

The  argument `h` is a linear form on  the roots of `H`, given by its value
on  the simple roots; this  linear form is extended  to the roots of `W` by
`0`  on the orthogonal of the roots  of `H`; and finally the resulting form
is  conjugated by an element of `W` so that it takes positive values on the
simple roots. If the initial form describes a Dynkin-Richardson diagram for
`H`, the result will describe a Dynkin-Richardson diagram for `W`.

```julia-repl
julia> W=coxgroup(:F,4)
F₄

julia> H=reflection_subgroup(W,[1,3])
F₄₍₁₃₎=A₁×Ã₁Φ₁²

julia> induced_linear_form(W,H,[2,2])
4-element Vector{Int64}:
 0
 1
 0
 0

julia> uc=UnipotentClasses(W);

julia> uc.classes[4].dynkin
4-element Vector{Int64}:
 0
 1
 0
 0

julia> uc.classes[4]
UnipotentClass(A₁+Ã₁)
```
The  example above shows that the class containing the regular class of the
Levi subgroup of type `A₁×Ã₁` is the class `A₁+Ã₁`.
"""
function induced_linear_form(W,K,h)
  if semisimplerank(K)==0 return fill(0,semisimplerank(W)) end
  h=vcat(h,fill(0,rank(W)-semisimplerank(K)))
  h=Int.(invbaseX(K)*h)
  r=roots(W)
  v=dot.(view(r,1:W.N),Ref(h))
  w=inv(with_inversions(W,findall(<(0),v)))
  dot.(view(r,action.(Ref(W),eachindex(gens(W)),w)),Ref(h))
end

"""
`distinguished_parabolics(W)`

the  list of distinguished (in the  sense of Richardson) standard parabolic
subgroups `W_I` of `W`, each given by the list `I` of indices in the simple
reflections. The distinguished unipotent conjugacy classes of the reductive
group  `𝐆` with Weyl group `W` consist of  the dense unipotent orbit in the
unipotent  radical of a  parabolic subgroup `𝐏_I`  of `𝐆` associated with a
distinguished  `W_I`.  Their  Dynkin-Richardson  diagram  contains  a  0 at
indices  `I` and a 2 in other entries. If `𝐏_I=𝐋𝐔` is a Levi decomposition,
`𝐏_I` is distinguished iff `dim 𝐋∩𝐆'=dim 𝐔/𝐔'`.

```julia-repl
julia> W=coxgroup(:F,4)
F₄

julia> distinguished_parabolics(W)
4-element Vector{Vector{Int64}}:
 []
 [3]
 [1, 3]
 [1, 3, 4]
```
"""
function distinguished_parabolics(W)
  filter(combinations(eachindex(gens(W)))) do J
    p=fill(1,ngens(W))
    p[J].=0
    p=dot.(view(W.rootdec,1:nref(W)),Ref(p))
    2*count(iszero,p)+ngens(W)==count(isone,p)
  end
end

function BalaCarterLabels(W)
  vcat(map(parabolic_reps(W)) do J
    L=reflection_subgroup(W,J)
    map(distinguished_parabolics(L))do D
      w=fill(2,length(J));w[D].=0
      u=copy(J);u[D]=-u[D]
      [induced_linear_form(W,L,w),u]
    end
  end...)
end

"""
`QuotientAu(Au,chars)` 

`chars` is a list of indices of characters of `Au`.
If  `k` is the common kernel of `chars`, `QuotientAu` returns a `NamedTuple`
`(Au=>Au/k,
  chars=>`index of chars in `Irr(Au/k),
  gens=>words` in `gens(Au)` preimages of `gens(Au/k))`

Since we have problems with quotient groups, we are forced to program an ad
hoc  solution which  works only  for `Au`  actually occuring  for unipotent
classes of a reductive group `G`.
"""
function QuotientAu(Au,chars)
  AbGens=function(g)
    res=empty(gens(g))
    l=gens(g)
    while !isempty(l)
      sort!(l,by=order,rev=true)
      if isempty(res) t=[l[1]] else t=Ref(l[1]).*elements(Group(res)) end
      p=findfirst(x->order(x)<order(l[1]),t)
      if !isnothing(p)
        t=t[p]
	if order(t)>1 l[1]=t
	else l=deleteat!(l,1)
        end
      else push!(res,l[1])
      end
    end
    res
  end
  # q=Au/k,  ww=words in q images of gens(Au)
  finish=function(q,ww)
    h=Hom(Au,q,map(x->q(x...),ww))
    fusion=map(c->position_class(q,h(c)),classreps(Au))
    ctu=CharTable(Au).irr
    cth=CharTable(q).irr
    ch(c)=map(j->ctu[c,findfirst(==(j),fusion)],1:nconjugacy_classes(q))
    return (Au=q,
      chars=map(c->findfirst(i->cth[i,:]==ch(c),axes(cth,1)),chars),
      gens=map(x->word(Au,elements(Au)[findfirst(y->h(y)==x,elements(Au))]),gens(q)))
  end
  Z=n->crg(n,1,1)
# println("Au=$Au chars=$chars")
  ct=transpose(CharTable(Au).irr[chars,:])
  cl=filter(i->(@view ct[i,:])==(@view ct[1,:]),axes(ct,1))
# println("ct=$ct cl=$cl")
  if length(cl)==1 return (Au=Au,chars=chars,
                              gens=map(x->[x],eachindex(gens(Au)))) end
  ct=transpose(toM(unique!(sort(toL(ct)))))
# println("ct=$ct")
# k=Subgroup(Au,filter(x->position_class(Au,x) in cl,elements(Au)))
  k=Group(filter(x->position_class(Au,x) in cl,elements(Au)))
  if length(k)==length(Au) return (Au=coxgroup(),chars=[1],gens=[])
  end
# println("Au=$Au k=$k")
  if semisimplerank(Au)==1 return finish(Z(div(length(Au),length(k))),[[1]])
  elseif isabelian(Au/k)
    q=Au/k
    q=Group(AbGens(q))
    h=Hom(Au,q,map(x->NormalCoset(k,x),gens(Au)))
    f=map(x->word(q,h(x)),gens(Au))
#  Print(prod(map(x->Z(order(x)),gens(q)))," ",f,"\n");
    return finish(prod(map(x->Z(order(x)),gens(q))),f)
  else
    p=findfirst(t->all(x->x in reflection_subgroup(Au,t.indices),
                elements(k)),refltype(Au))
    if !isnothing(p)
      p=refltype(Au)[p].indices
      if length(k)==length(reflection_subgroup(Au,p))
       return finish(reflection_subgroup(Au,setdiff(eachindex(gens(Au)),p)),
                     map(i->i in p ? [] : [i],eachindex(gens(Au))))
      elseif length(p)==1
        t=copy(refltype(Au))
        p=findfirst(t->t.indices==p,refltype(Au))
	t[p].p/=length(k)
        return finish(ReflectionGroup(t...),map(x->[x],eachindex(gens(Au))))
      end
     elseif xrepr(rio(),Au)=="A₁×B₂" && length(k)==2 && longest(Au) in k
      return finish(coxgroup(:B,2),[[1,2,1,2],[1],[2]])
     elseif xrepr(rio(),Au)=="A₂×A₁×A₁" && length(k)==2 &&
      longest(reflection_subgroup(Au,[3,4])) in k
      return finish(coxgroup(:A,2)*coxgroup(:A,1),[[1],[2],[3],[3]])
    end
  end
  # Print(" Au=",xrepr(rio(),Au)," sub=",map(e.Get,gens(k)),"\n");
  error("not implemented ",xrepr(rio(),Au),chars)
# q:=Au/k; f:=FusionConjugacyClasses(Au,q); Print(" quot=",q," fusion=",f,"\n");
# return rec(Au:=Au,chars:=chars);
end

# When some Springer series have been suppressed/weeded out, we  quotient
# the Au's by the common  kernel of the remaining characters of the Au's.
function AdjustAu!(classes,springerseries)
  for (i, u) in pairs(classes)
    l=map(s->filter(k->s[:locsys][k][1]==i,eachindex(s[:locsys])),
          springerseries)
#   println(springerseries)
    chars=reduce(vcat,map(j->last.(springerseries[j][:locsys][l[j]]),
                   eachindex(l)))
    f=QuotientAu(u.Au,chars)
#   if Size(Au)<>Size(f.Au) then
#     Print("class ",i,"=",classes[i].name," ",[Au,chars],"=>",f,"\n");
#   fi;
#   print(u.name,":");ds(u.AuAction)
    u.Au=f.Au
    if haskey(u,:AuAction)
      R=u.AuAction.group
      if rank(R)==0
        u.AuAction=ExtendedCox(R,[fill(0,0,0) for x in f.gens])
      else
       if isempty(f.gens) F0s=[reflrep(R,R())]
       else F0s=map(x->prod(u.AuAction.F0s[x]),f.gens)
       end
       u.AuAction=ExtendedCox(R,F0s)
      end
#     u.AuAction.phis=map(x->prod(u.AuAction.phis[x]),f.gens)
    end
    k=1
    for j in eachindex(l)
      springerseries[j][:locsys]=copy(springerseries[j][:locsys])
      for s in l[j]
        springerseries[j][:locsys][s][2]=f.chars[k]
        k+=1
      end
    end
  end
end

# unipotent classes for the simply connected quasisimple group of type t
function UnipotentClasses(t::TypeIrred,p=0)
  uc=getchev(t,:UnipotentClasses,p)
  if uc===nothing error("no UnipotentClasses for ",t) end
  uc[:orderClasses]=map(uc[:orderClasses])do v
    if v isa String && isempty(v) return Int[] end
    Vector{Int}(v)
  end
  rank=PermRoot.rank(t)
  c=haskey(t,:orbit) ? cartan(t.orbit[1]) : cartan(t)
  rr=toM(roots(c))
  classes=map(uc[:classes])do u # fill omitted fields
    name=u[:name]
    parameter= haskey(u,:parameter) ? u[:parameter] : u[:name]
    dimBu= haskey(u,:dimBu)  ? u[:dimBu] : -1
    delete!.(Ref(u),[:name,:parameter,:dimBu])
    cl=UnipotentClass(name,parameter,dimBu,u)
    if haskey(cl,:dynkin)
      weights=rr*cl.dynkin
      n0=count(iszero,weights)
      if cl.dimBu==-1 cl.dimBu=n0+div(count(isone,weights),2)
      elseif cl.dimBu!=n0+div(count(isone,weights),2) error("theory")
      end
      n0=2*n0-count(==(2),weights)
      cl.dimunip=2*cl.dimBu-n0
      cl.dimred=n0+rank
    elseif haskey(cl,:red)
      cl.dimred=dimension(cl.red)
      cl.dimunip=2*cl.dimBu+rank-cl.dimred
    elseif haskey(cl,:dimred)
      cl.dimunip=2*cl.dimBu+rank-cl.dimred
    end
    cl
  end
  springerseries=uc[:springerSeries]
  for s in springerseries
    if isempty(s[:levi]) s[:levi]=Int[] end
#   s[:levi]=indices(t)[s[:levi]]
    s[:locsys]=Vector{Int}.(s[:locsys])
  end
  orderclasses=Poset(CPoset(uc[:orderClasses]),classes)
  delete!.(Ref(uc),[:classes,:orderClasses,:springerSeries])
# uc[:spets]=t
  UnipotentClasses(classes,p,orderclasses,springerseries,uc)
end

Base.length(uc::UnipotentClasses)=length(uc.classes)

"""
`UnipotentClasses(W[,p])`

`W` should be a `FiniteCoxeterGroup` record for a Weyl group or `rootdatum`
describing  a reductive algebraic group `𝐆`.  The function returns a record
containing information about the unipotent classes of `𝐆` in characteristic
`p`  (if omitted, `p`  is assumed to  be any good  characteristic for `𝐆`).
This contains the following fields:

`group`: a pointer to `W`

`p`: the characteristic of the field for which the unipotent classes were
computed. It is `0` for any good characteristic.

`orderclasses`:  a list describing the Hasse diagram of the partial order
induced   on   unipotent   classes   by   the  closure  relation.  That  is
`.orderclasses[i]`  is the list of `j` such that ``C̄ⱼ⊋ C̄ᵢ``  and  there  is
no  class  ``Cₖ``  such  that ``C̄ⱼ⊋ C̄ₖ⊋ C̄ᵢ``.

`classes`:  a list of records holding information for each unipotent class.
See  the  help  for  [`UnipotentClass`](@ref)  for  a  description of these
records.

`springerseries`:  a list of records, each  of which describes a Springer
series  of `𝐆`.

The  records  describing  individual  Springer  series  have  the following
fields:

`levi`:the  indices of the  reflections corresponding to  the Levi subgroup
`𝐋`  where  lives  the  cuspidal  local  system `ι` from which the Springer
series is induced.

`relgroup`: The relative Weyl group ``N_𝐆(𝐋,ι)/𝐋``. The first series is the
principal series for which `.levi=[]` and `.relgroup=W`.

`locsys`:  a  list  of  length  `nconjugacy_classes(.relgroup)`, holding in
`i`-th  position a  pair describing  which local  system corresponds to the
`i`-th  character of  ``N_𝐆(𝐋,ι)``. The  first element  of the  pair is the
index  of the concerned unipotent class `u`, and the second is the index of
the corresponding character of `A(u)`.

`Z`:  the central character associated  to the Springer series, specified
by its value on the generators of the center.

```julia-repl
julia> W=rootdatum(:sl,4)
sl₄

julia> uc=UnipotentClasses(W);

julia> uc.classes
5-element Vector{UnipotentClass}:
 UnipotentClass(1111)
 UnipotentClass(211)
 UnipotentClass(22)
 UnipotentClass(31)
 UnipotentClass(4)
```
The  `show`  function  for  unipotent  classes  accepts  all the options of
`showtable`  and  of  `charnames`.  Giving  the  option  `mizuno`  (resp.
`shoji`)  uses  the  names  given  by  Mizuno  (resp.  Shoji) for unipotent
classes.  Moreover,  there  is  also  an  option  `fourier` which gives the
Springer  correspondence tensored with the  sign character of each relative
Weyl  group,  which  is  the  correspondence obtained via a Fourier-Deligne
transform  (here  we  assume  that  `p`  is  very  good, so that there is a
nondegenerate  invariant bilinear form on the Lie algebra, and also one can
identify nilpotent orbits with unipotent classes).

Here is how to display the non-cuspidal part of the Springer correspondence
of  the unipotent  classes of  `E₆` using  the notations  of Mizuno for the
classes  and those  of Frame  for the  characters of  the Weyl group and of
Spaltenstein  for the characters  of `G₂` (this  is convenient for checking
our data with the original paper of Spaltenstein):

```julia-rep1
julia> uc=UnipotentClasses(rootdatum(:E6sc));

julia> xdisplay(uc;cols=[5,6,7],spaltenstein=true,frame=true,mizuno=true,
      order=false)
┌──────┬─────────────────────────────────────────────────────┐
│u     │             E₆(Φ₁⁶) G₂(A₂×A₂Φ₁²)/ζ₃ G₂(A₂×A₂Φ₁²)/ζ₃²│
├──────┼─────────────────────────────────────────────────────┤
│E₆    │                1:1ₚ            ζ₃:1            ζ₃²:1│
│E₆(a₁)│                1:6ₚ           ζ₃:εₗ           ζ₃²:εₗ│
│D₅    │              Id:20ₚ                                 │
│A₅+A₁ │        -1:15ₚ 1:30ₚ           ζ₃:θ′           ζ₃²:θ′│
│A₅    │              1:15_q           ζ₃:θ″           ζ₃²:θ″│
│D₅(a₁)│              Id:64ₚ                                 │
│A₄+A₁ │              Id:60ₚ                                 │
│D₄    │              Id:24ₚ                                 │
│A₄    │              Id:81ₚ                                 │
│D₄(a₁)│111:20ₛ 3:80ₛ 21:90ₛ                                 │
│A₃+A₁ │              Id:60ₛ                                 │
│2A₂+A₁│               1:10ₛ          ζ₃:ε_c          ζ₃²:ε_c│
│A₃    │             Id:81ₚ′                                 │
│A₂+2A₁│             Id:60ₚ′                                 │
│2A₂   │              1:24ₚ′            ζ₃:ε            ζ₃²:ε│
│A₂+A₁ │             Id:64ₚ′                                 │
│A₂    │      11:15ₚ′ 2:30ₚ′                                 │
│3A₁   │            Id:15_q′                                 │
│2A₁   │             Id:20ₚ′                                 │
│A₁    │              Id:6ₚ′                                 │
│1     │              Id:1ₚ′                                 │
└──────┴─────────────────────────────────────────────────────┘
```
"""
function UnipotentClasses(W::Union{FiniteCoxeterGroup,CoxeterCoset},p=0)
# println("UnipotentClasses(",W,")")
  if !(p in (W isa Spets ? badprimes(Group(W)) : badprimes(W))) p=0 end
  get!(W,Symbol("unipotentclasses",p))do
  spetscase=W isa Spets
  if spetscase
    WF=W
    W=Group(WF)
    t=refltype(W)
    l=map(x->findfirst(y->any(z->sort(x.indices)==sort(z.indices),
                                       y.orbit),refltype(WF)),t)
    uc=UnipotentClasses.(refltype(WF)[l],p)
   else WF=spets(W)
    t=refltype(W)
    uc=UnipotentClasses.(t,p)
  end
  if isempty(t)
    classes=[UnipotentClass("1",[],0,
        Dict(:Au=>coxgroup(),:dynkin=>[],:balacarter=>[],
             :dimunip=>0,:red=>torus(rank(W))))]
    uc=[UnipotentClasses(classes,p,Poset(CPoset([Int[]])),
      [Dict(:Z=>Int[],:levi=>Int[],:locsys=>[[1,1]],:relgroup=>coxgroup())],
      Dict{Symbol,Any}(:spets=>W))]
    l=Vector{Int}[]
  else
    classes=map(cartesian(map(x->x.classes,uc)...)) do v
      l=indices.(t)
      if length(v)==1 && issorted(l[1]) u=deepcopy(v[1])
      else
        u=UnipotentClass(join(map(x->x.name,v),","),map(x->x.parameter,v),
                         sum(map(x->x.dimBu,v)),Dict{Symbol,Any}())
        u.Au=prod(x->x.Au,v)
        if all(x->haskey(x,:dimred),v)
          u.dimred=sum(x->x.dimred,v)+rank(W)-semisimplerank(W)
        end
        if all(x->haskey(x,:dimunip),v)
          u.dimunip=sum(x->x.dimunip,v) end
        if all(x->haskey(x,:red),v)
#         println(map(x->x.red,v))
          u.red=reduce(Cosets.extprod,map(x->x.red,v)) end
        if all(x->haskey(x,:AuAction),v)
          u.AuAction=prod(x->x.AuAction,v) end
        if all(x->haskey(x,:dynkin),v)
          u.dynkin=zeros(Int,sum(length,l))
          for i in 1:length(l) u.dynkin[l[i]]=v[i].dynkin end
        end
      end
      if all(x->haskey(x,:balacarter),v)
        u.balacarter=reduce(vcat,map(j->j>0 ? x[j] : -x[-j],v[i].balacarter)
                                                   for (i,x) in pairs(l))
      end
      if haskey(u, :red)
        if !(u.red isa Spets) u.red=spets(u.red) end
        if rank(W)>semisimplerank(W) && haskey(u, :red)
          u.red=Cosets.extprod(u.red,radical(WF))
          if haskey(u,:AuAction)
            T=radical(W)
            u.AuAction=ExtendedCox(u.AuAction.group*T,
               map(x->cat(x,reflrep(T,T()),dims=(1,2)),u.AuAction.F0s))
          end
        end
      end
      u
    end
  end
  if !(p in badprimes(W)) && !haskey(classes[1],:balacarter)
    bc=BalaCarterLabels(W)
#   println("W=$W bc=$bc")
#   println(map(u->u.dynkin,classes))
    for u in classes
      pp=findfirst(p->p[1]==u.dynkin,bc)
      if isnothing(pp) error("not found:",u.dynkin) end
      u.balacarter=bc[findfirst(p->p[1]==u.dynkin,bc)][2]
    end
  end
  ll=length.(uc)
  orderclasses=map(cartesian(map(x->1:x,ll)...)) do v
    o=cartesian(map(j->vcat(hasse(uc[j].orderclasses)[v[j]],[v[j]]),1:length(v))...)
    o=map(x->cart2lin(ll,x),o)
    setdiff(o,[cart2lin(ll,v)])
  end
  springerseries=map(cartesian(map(x->x.springerseries,uc)...)) do v
#   if isempty(v) return Dict(:Z=>[],:levi=>[],:locsys=>[[1,1]])
    if length(v)==1
      if !isempty(l)
        v[1]=deepcopy(v[1]); v[1][:levi]=l[1][v[1][:levi]];
      end
      return v[1]
    end
    s=Dict{Symbol,Any}(:levi=>reduce(vcat,map(i->l[i][v[i][:levi]],eachindex(v))))
    s[:Z]=reduce(vcat,getindex.(v,:Z))
    s[:locsys]=map(cartesian(getindex.(v,:locsys)...)) do v
      v=collect.(zip(v...))
      u=map(i->nconjugacy_classes(uc[i].classes[v[1][i]].Au),
              eachindex(v[1]))
      [cart2lin(ll,v[1]),cart2lin(u,v[2])]
    end
    if all(haskey.(v,:parameter)) s[:parameter]=getindex.(v,:parameter) end
    s[:relgroup]=prod(getindex.(v,:relgroup))
    if length(v)==1
      for k in setdiff(keys(v[1]),[:levi,:Z,:locsys,:parameter])
        s[k]=v[1][k]
      end
    end
    s
  end
  if length(uc)==1 prop=uc[1].prop else prop=Dict{Symbol,Any}() end
  prop[:spets]=spetscase ? WF : W
  if spetscase
    springerseries=filter(x->sort(inclusion(W,x[:levi]).^WF.phi)==
                             sort(inclusion(W,x[:levi])),springerseries)
  end
# To deal with a general group intermediate between Gad and Gsc, we discard
# the  Springer series  corresponding to  a central  character which is not
# trivial on the fundamental group (seen as a subgroup of ZGsc)
# algebraic_center(W).descAZ returns the generators of the fundamental group
# of  the  algebraic  group  W  as  words  in  generators  of  the absolute
# fundamental group.
  if !all(x->unique(x[:Z])==[1],springerseries)
    springerseries=filter(s->all(y->prod(s[:Z][y])==1,
             algebraic_center(W).descAZ),springerseries)
    AdjustAu!(classes,springerseries)
  end
  if spetscase
    g=Group(weightinfo(W)[:AdjointFundamentalGroup])
    permZ=map(x->word(g,x),gens(g).^WF.phi)
    springerseries=filter(x->map(i->prod(x[:Z][i]),permZ)==x[:Z],springerseries)
  end
# println(springerseries[1])
  s=springerseries[1]
  if spetscase
    s[:relgroup]=relative_coset(WF,s[:levi])
    s[:locsys]=s[:locsys][charinfo(s[:relgroup]).charRestrictions]
  end
  l=filter(i->any(y->i==y[1],s[:locsys]),1:length(classes))
  s[:locsys]=map(((c,s),)->[findfirst(==(c),l),s],s[:locsys])
  # for now only springerseries[1] properly twisted
  for s in springerseries[2:end]
    if spetscase
      s[:relgroup]=relative_coset(WF,s[:levi])
      s[:locsys]=s[:locsys][charinfo(s[:relgroup]).charRestrictions]
    end
    s[:locsys]=map(((c,s),)->[findfirst(==(c),l),s],s[:locsys])
  end
  classes=classes[l]
  AdjustAu!(classes,springerseries)
  orderclasses=Poset(induced(CPoset(orderclasses),l),classes)
  orderclasses.show_element=(io,x,n)->print(io,name(io,x.elements[n]))
  ucl=UnipotentClasses(classes,p,orderclasses,springerseries,prop)
  ucl
  end
end

FinitePosets.Poset(uc::UnipotentClasses)=uc.orderclasses

function reflection_name(io::IO,W)
  r=join(getchev(W,:ReflectionName,IOContext(io,:TeX=>true).dict),"×")
  fromTeX(io,r)
end

function showcentralizer(io::IO,u)
  c=""
  function AuName(u)
    if length(u.Au)==1 return "" end
    res=haskey(u,:AuAction) ||
        (haskey(u,:dimred) && iszero(u.dimred)) ? "." : "?"
        res*=reflection_name(IOContext(io,:Au=>true),u.Au)
  end
  if haskey(u,:dimunip)
    if u.dimunip>0 c*=sprint(show,Mvp(:q)^u.dimunip;context=io) end
  else c*="q^?" end
  if haskey(u,:AuAction)
    if rank(u.red)>0
      c*="."
      if all(isone,u.AuAction.F0s)
        c*=xrepr(io,u.red)*AuName(u)
      elseif length(u.Au)==1 ||
         length(u.Au)==length(Group(u.AuAction.phis...))
        if length(u.Au)==1 || isone(u.red.phi)
          c*=xrepr(io,u.AuAction)
        else
          c*="["*xrepr(io,u.red)*"]"*xrepr(io,u.AuAction)
        end
      else
        c*=xrepr(io,u.AuAction)*AuName(u)
      end
    else
      c*=AuName(u)
    end
  elseif haskey(u,:red)
    n=xrepr(io,u.red)
    if n!="." c*="."*n end
    c*=AuName(u)
  else
    if haskey(u,:dimred) if u.dimred>0 c*="[red dim ",u.dimred,"]" end
    else c*="[red??]"
    end
    c*=AuName(u)
  end
  c=replace(c,r"^\."=>"")
  c=replace(c,r"\.\."=>".")
  replace(c,"()"=>"")
  c
end

function Base.show(io::IO, ::MIME"text/html", uc::UnipotentClasses)
  show(IOContext(io,:TeX=>true),"text/plain",uc)
end

function Base.show(io::IO,uc::UnipotentClasses)
  TeX=get(io,:TeX,false)
  print(io,TeX ? "\$\\mbox{UnipotentClasses}" : "UnipotentClasses")
  print(io,"(",uc.spets)
  if uc.p!=0 print(io,",",uc.p) end
  print(io,")")
  if TeX print(io,"\$") end
end

function Base.show(io::IO,::MIME"text/plain",uc::UnipotentClasses)
  println(io,uc)
  TeX=get(io,:TeX,false)
  if get(io,:order,true) println(io,uc.orderclasses) end
  sp = map(copy, uc.springerseries)
  if get(io,:fourier,false)
    for p in sp p[:locsys] = p[:locsys][detPerm(p[:relgroup])] end
  end
  WF=uc.spets
  if WF isa Spets W=Group(WF)
  else W=WF;WF=spets(W) end
  if uc.p!=0 || !any(x->haskey(x,:balacarter),uc.classes)
    io=IOContext(io,:balacarter=>false)
  end
  tbl = map(uc.classes)do u
    res= iszero(uc.p) ? [joindigits(u.dynkin)] : String[]
    push!(res, string(u.dimBu))
    if get(io,:balaCarter,true)
      if haskey(u, :balacarter)
        b=fill('.',ngens(W))
        for i in u.balacarter if i>0 b[i]='2' else b[-i]='0' end end
      else
        b=fill('?',ngens(W))
      end
      push!(res, String(b))
    end
    if get(io,:centralizer,true)
      push!(res,showcentralizer(io,u))
    end
    if get(io,:springer,true)
      i=findfirst(==(u),uc.classes)
      cc(ss)=map(function (i)
                c1 = charnames(io,u.Au)[ss[:locsys][i][2]]
                c2 = charnames(io,ss[:relgroup])[i]
                (c1=="") ? c2 : c1*":"*c2
           end, findall(y->y[1]==i,ss[:locsys]))
      append!(res, map(ss->join(cc(ss), TeX ? "\\kern 0.8em " : " "),sp))
    end
    res
  end
  col_labels= String[]
  if iszero(uc.p) push!(col_labels,"\\mbox{D-R}") end
  push!(col_labels, TeX ? "\\dim{\\cal B}_u" : "dℬ ᵤ")
  if get(io,:balaCarter,true) push!(col_labels, "\\mbox{B-C}") end
  if get(io,:centralizer,true) push!(col_labels, "C_{\\bf G}(u)") end
  if get(io,:springer,true)
    append!(col_labels,
      map(function (ss,)
        res=string(xrepr(io,ss[:relgroup]),"(",
          xrepr(io,subspets(WF,ss[:levi]);parent=false),")")
        if !all(isone,ss[:Z])
          res*=string("/", join(map(q->xrepr(io,q),ss[:Z]),","))
        end
        return res
    end, sp))
  end
  row_labels=name.(Ref(io),uc.classes)
  if get(io,:rows,false)==false
    io=IOContext(io,:rows=>sortperm(map(x->x.dimBu, uc.classes)))
  end
  showtable(io,toM(tbl);rows_label="u",col_labels=col_labels,row_labels=row_labels)
end

# decompose tensor product of characters (given as their indices in CharTable)
function decompose_tensor(W,c::Int...)
  ct=CharTable(W)
# println("eltype=",eltype(irr))
  decompose(ct,vec(prod(view(ct.irr,collect(c),:),dims=1)))
end

@GapObj struct ICCTable end

"""
`ICCTable(uc,seriesNo=1;q=Pol())`

This function gives the table of decompositions of the functions ``X_ι`` in
terms  of the functions ``Y_ι``. Here `ι` is a `𝐆`-equivariant local system
on  the  class  `C`  of  a  unipotent  element  `u`. Such a local system is
parametrised  by the pair  `(u,ϕ)` of `u`  and a character  of the group of
components   `A(u)`   of   ``C_𝐆   (u)``.   The  function  ``Y_ι``  is  the
characteristic   function  of  this   local  system  and   ``X_ι``  is  the
characteristic   function  of  the  corresponding  intersection  cohomology
complex  on `C̄`. The  Springer correspondence says  that the local systems
can  also be  indexed by  characters of  a relative  Weyl group.  Since the
coefficient of `Xᵪ` on `Yᵩ` is `0` if `χ` and `φ` are not characters of the
same  relative Weyl group (are not in  the same Springer series), the table
is  for one  Springer series,  specified by  the argument  'seriesNo' (this
defaults  to 'seriesNo=1' which is the principal series). The decomposition
multiplicities  are graded,  and are  given as  polynomials in one variable
(specified by the argument `q`; if not given `Pol()` is assumed).

```julia-repl
julia> uc=UnipotentClasses(coxgroup(:A,3));t=ICCTable(uc)
Coefficients of Xᵪ on Yᵩ for series L=A₃₍₎=Φ₁³ W_G(L)=A₃
┌─────┬────────────────┐
│X\\Y  │4 31 22 211 1111│
├─────┼────────────────┤
│X4   │1  1  1   1    1│
│X31  │.  1  1  Φ₂   Φ₃│
│X22  │.  .  1   1   Φ₄│
│X211 │.  .  .   1   Φ₃│
│X1111│.  .  .   .    1│
└─────┴────────────────┘
```
In  the  above  the  multiplicities  are  given  as  products of cyclotomic
polynomials  to display them  more compactly. However  the format of such a
table can be controlled more precisely.

For  instance,  one  can  ask  to  not  display  the entries as products of
cyclotomic polynomials and not display the zeroes as '.'

```julia-rep1
julia> xdisplay(t;cycpol=false,dotzero=false)
Coefficients of Xᵪ on Yᵩ for A3
┌─────┬──────────────────┐
│X\\Y  │4 31 22 211   1111│
├─────┼──────────────────┤
│X4   │1  1  1   1      1│
│X31  │0  1  1 q+1 q²+q+1│
│X22  │0  0  1   1   q²+1│
│X211 │0  0  0   1 q²+q+1│
│X1111│0  0  0   0      1│
└─────┴──────────────────┘
```

Since  `show`  uses  the  function  `showtable`,  all  the  options of this
function  are  also  available.  We  can  use  this to restrict the entries
displayed  to a  given sublist  of the  rows and  columns (here the indices
correspond  to the number  in Chevie of  the corresponding character of the
relative Weyl group of the given Springer series):

```julia-rep1
julia> uc=UnipotentClasses(coxgroup(:F,4));
julia> t=ICCTable(uc);
julia> sh=[13,24,22,18,14,9,11,19];
julia> xdisplay(t,rows=sh,cols=sh)
Coefficients of Xᵪ on Yᵩ for series L=F₄₍₎=Φ₁⁴ W_G(L)=F₄
┌───────┬────────────────────────────────────────────┐
│X\\Y    │A₁+Ã₁ A₂ Ã₂ A₂+Ã₁ Ã₂+A₁ B₂⁽¹¹⁾ B₂ C₃(a₁)⁽¹¹⁾│
├───────┼────────────────────────────────────────────┤
│Xφ₉‚₁₀ │    1  .  .     .     .      .  .          .│
│Xφ″₈‚₉ │    1  1  .     .     .      .  .          .│
│Xφ′₈‚₉ │    1  .  1     .     .      .  .          .│
│Xφ″₄‚₇ │    1  1  .     1     .      .  .          .│
│Xφ′₆‚₆ │   Φ₄  1  1     1     1      .  .          .│
│Xφ₄‚₈  │   q²  .  .     .     .      1  .          .│
│Xφ″₉‚₆ │   Φ₄ Φ₄  .     1     .      .  1          .│
│Xφ′₄‚₇ │   q²  . Φ₄     .     1      .  .          1│
└───────┴────────────────────────────────────────────┘
```
The  `ìo` option `rowlocsys=true`  will display local  systems also for the
row labels.

The   function  'ICCTable'  returns  an   object  with  various  pieces  of
information which can help further computations.

`.scalar`:  this contains the table of  multiplicities `Pᵪᵩ` of the `Xᵪ` on
the  `Yᵩ`.  One  should  pay  attention  that  by default, the table is not
displayed  in the same order as the  stored |.scalar|, which is in order in
Chevie  of  the  characters  in  the  relative  Weyl  group;  the  table is
transposed,  then lines  and rows  are sorted  by `dimBu,class  no,index of
character in A(u)` while displayed.

`.group`: The group `W`.

`.relgroup`: The relative Weyl group for the Springer series.

`.series`: The index of the Springer series given for `W`.

`.dimBu`: The list of ``dimℬᵤ`` for each local system `(u,φ)` in the series.

`:L`:  The matrix of (unnormalised) scalar products of the functions ``Yᵩ``
with  themselves,  that  is  the  ``(φ,χ)``  entry  is ``∑_{g∈𝐆(𝔽_q)} Yᵩ(g)
Yᵪ(g)``. This is thus a symmetric, block-diagonal matrix where the diagonal
blocks  correspond to geometric unipotent conjugacy classes. This matrix is
obtained as a by-product of Lusztig's algorithm to compute ``Pᵩᵪ``.
"""
function ICCTable(uc::UnipotentClasses,i=1;q=Pol())
  W=uc.spets # W=Group(uc.spets)
  if W isa Spets W=W.W end
  ss=uc.springerseries[i]
  res=ICCTable(Dict(:spets=>uc.spets,:relgroup=>ss[:relgroup],
                    :series=>i,:q=>q,:p=>uc.p))
# We are going to solve the equation in "unipotent support", page 151
# ᵗPΛP=Ω  where $Λ_{i,j}$ is  $∑_{g∈ G^F} Yᵢ(g)Ȳⱼ(g)$ and $Ω_{i,j}$ is equal
# to # $|Z^0(G^F)|q^{-semisimple rank L}|G^F|/P(W_G(L))
#  q^{-bᵢ-bⱼ}FakeDegree(χᵢ⊗χⱼ⊗sgn)$
# where $P(W_G(L))$ is the Poincare polynomial $∏ᵢ(q^{dᵢ}-1)$
# where $dᵢ$ are the reflection degrees of $W_G(L)$
# res[:scalar] is the matrix $P$
  R=ss[:relgroup]
  ct=CharTable(R)
  k=charinfo(R).positionDet
# Partition on characters of ss.relgroup induced by poset of unipotent classes
  res.dimBu=map(x->uc.classes[x[1]].dimBu,ss[:locsys])
  res.blocks=collectby(-res.dimBu,eachindex(ss[:locsys]))
  subst=!(q isa Pol)
  if subst var=q; q=Pol() end
  f=fakedegrees(R,q)
  n=length(f)
  # matrix of q^{-bᵢ-bⱼ}*fakedegree(χᵢ ⊗ χⱼ ⊗ sgn)
  tbl=bigcell_decomposition([q^(-res.dimBu[i]-res.dimBu[j])*
                            sum(map(*,f,decompose_tensor(R,i,j,k)))
     for i in 1:n,j in 1:n], res.blocks) # //1 needed in D7
  if subst tbl=map(y->map(x->x(var),y),tbl); q=var end
  res.scalar=tbl[1]
  res.locsys=ss[:locsys]
# res.L=tbl[2]*GenericOrder(W,q)/prod(ReflectionDegrees(R),d->q^d-1)/
#   q^(W.semisimplerank-R.semisimplerank);
  res.L=tbl[2]*q^(nref(W)+semisimplerank(R)-semisimplerank(W))
  res.uc=uc
  res.levi=reflection_subgroup(W,ss[:levi])
  if haskey(ss,:parameter) res.parameter=ss[:parameter]
  else res.parameter=(1:length(ss[:locsys])).+100*(i-1)
  end
  res.scalar=improve_type(res.scalar)
  res.L=improve_type(res.L)
  res
end

function Base.show(io::IO, ::MIME"text/html", x::ICCTable)
  show(IOContext(io,:TeX=>true),"text/plain",x)
end

Base.show(io::IO,x::ICCTable)=print(io,"ICCTable(",x.uc,",",x.series,")")

function Base.show(io::IO,::MIME"text/plain",x::ICCTable)
 printTeX(io,"Coefficients of \$X_\\chi\$ on \$Y_\\phi\$ for series \$L=",x.levi,"\$ \$W_G(L)=",x.relgroup,"\$\n")
  if get(io,:cols,false)==false && get(io,:rows,false)==false
    rows=collect(eachindex(x.dimBu))
    sort!(rows,by=i->[x.dimBu[i],x.locsys[i]])
    io=IOContext(io,:rows=>rows,:cols=>rows)
  end
  tbl=get(io,:cycpol,true) ? map(CycPol,x.scalar) : x.scalar
  col_labels=map(((c,s),)->name(IOContext(io,:locsys=>s),x.uc.classes[c]),
                  x.locsys)
  if get(io,:rowlocsys,false) row_labels=col_labels
  else row_labels=map(x->get(io,:TeX,false) ? "X_{$x}" : "X$x",charnames(io,x.relgroup))
  end
  showtable(io,transpose(tbl);row_labels=row_labels,col_labels=col_labels,
    rows_label="X\\Y", dotzero=get(io,:dotzero,true))
end

@GapObj struct XTable end

"""
`XTable(uc;classes=false)`

This  function presents  in a  different way  the information obtained from
`ICCTable`. Let ``X̃_{u,ϕ}=q^{1/2(codim C-dim Z(𝐋 ))}X_{u,ϕ}`` where `C` is
the  class of `u` and `Z(𝐋 )` is the center of Levi subgroup on which lives
the cuspidal local system attached to the local system `(u,ϕ)`.

Then  `XTable(uc)` gives the decomposition of the functions ``X̃_{u,ϕ}`` on
local   systems.  `t=XTable(uc,classes==true)`  gives  the  values  of  the
functions   ``X̃_{u,ϕ}``   on   unipotent   classes.   A   side  effect  of
`classes=true`  is  to  compute  the  cardinal  of  the unipotent conjugacy
classes,  available in `t.cardClass`; in this case displaying `t` will show
the  cardinal  of  the  centralizers  of  unipotent  elements, available in
`t.centClass`.


```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> XTable(UnipotentClasses(W))
Values of character sheaves X̃ᵪ on local systems φ
┌──────────┬───────────────────────────────────────────┐
│X̃ᵪ|φ      │   1 A₁ Ã₁ G₂(a₁)⁽¹¹¹⁾ G₂(a₁)⁽²¹⁾ G₂(a₁) G₂│
├──────────┼───────────────────────────────────────────┤
│X_φ₁‚₀^G₂ │   1  1  1           .          .      1  1│
│X_φ₁‚₆^G₂ │  q⁶  .  .           .          .      .  .│
│X_φ′₁‚₃^G₂│  q³  .  q           .          q      .  .│
│X_φ″₁‚₃^G₂│  q³ q³  .           .          .      .  .│
│X_φ₂‚₁^G₂ │ qΦ₈  q  q           .          .      q  .│
│X_φ₂‚₂^G₂ │q²Φ₄ q² q²           .          .      .  .│
│X_Id^.    │   .  .  .          q²          .      .  .│
└──────────┴───────────────────────────────────────────┘
```

The functions `X̃` in the first column are decorated by putting as an
exponent the relative groups ``W_𝐆 (𝐋)``.

```julia-repl
julia> t=XTable(UnipotentClasses(W);classes=true)
Values of character sheaves X̃ᵪ on unipotent classes
┌──────────┬─────────────────────────────────────────────────────────┐
│X̃ᵪ|class  │           1     A₁     Ã₁ G₂(a₁) G₂(a₁)₍₂₁₎ G₂(a₁)₍₃₎ G₂│
├──────────┼─────────────────────────────────────────────────────────┤
│X_φ₁‚₀^G₂ │           1      1      1      1          1         1  1│
│X_φ₁‚₆^G₂ │          q⁶      .      .      .          .         .  .│
│X_φ′₁‚₃^G₂│          q³      .      q     2q          .        -q  .│
│X_φ″₁‚₃^G₂│          q³     q³      .      .          .         .  .│
│X_φ₂‚₁^G₂ │         qΦ₈      q      q      q          q         q  .│
│X_φ₂‚₂^G₂ │        q²Φ₄     q²     q²      .          .         .  .│
│X_Id^.    │           .      .      .     q²        -q²        q²  .│
├──────────┼─────────────────────────────────────────────────────────┤
│|C_𝐆(u)|  │q⁶Φ₁²Φ₂²Φ₃Φ₆ q⁶Φ₁Φ₂ q⁴Φ₁Φ₂    6q⁴        2q⁴       3q⁴ q²│
└──────────┴─────────────────────────────────────────────────────────┘

julia> XTable(UnipotentClasses(W,2))
Values of character sheaves X̃ᵪ on local systems φ
┌──────────┬──────────────────────────────────────────────────┐
│X̃ᵪ|φ      │   1 A₁ Ã₁ G₂(a₁)⁽¹¹¹⁾ G₂(a₁)⁽²¹⁾ G₂(a₁) G₂⁽¹¹⁾ G₂│
├──────────┼──────────────────────────────────────────────────┤
│X_φ₁‚₀^G₂ │   1  1  1           .          .      1      .  1│
│X_φ₁‚₆^G₂ │  q⁶  .  .           .          .      .      .  .│
│X_φ′₁‚₃^G₂│  q³  .  q           .          q      .      .  .│
│X_φ″₁‚₃^G₂│  q³ q³  .           .          .      .      .  .│
│X_φ₂‚₁^G₂ │ qΦ₈  q  q           .          .      q      .  .│
│X_φ₂‚₂^G₂ │q²Φ₄ q² q²           .          .      .      .  .│
│X_Id^.    │   .  .  .          q²          .      .      .  .│
│X_Id^.    │   .  .  .           .          .      .      q  .│
└──────────┴──────────────────────────────────────────────────┘

julia> XTable(UnipotentClasses(rootdatum(:sl,4)))
Values of character sheaves X̃ᵪ on local systems φ
┌────────┬────────────────────────────────────────────┐
│X̃ᵪ|φ    │1111 211 22⁽¹¹⁾ 22 31 4 4^(ζ₄) 4⁽⁻¹⁾ 4^(ζ₄³)│
├────────┼────────────────────────────────────────────┤
│X₁₁₁₁^A₃│  q⁶   .      .  .  . .      .     .       .│
│X₂₁₁^A₃ │q³Φ₃  q³      .  .  . .      .     .       .│
│X₂₂^A₃  │q²Φ₄  q²      . q²  . .      .     .       .│
│X₃₁^A₃  │ qΦ₃ qΦ₂      .  q  q .      .     .       .│
│X₄^A₃   │   1   1      .  1  1 1      .     .       .│
│X₁₁^A₁  │   .   .     q³  .  . .      .     .       .│
│X₂^A₁   │   .   .     q²  .  . .      .     q       .│
│X_Id^.  │   .   .      .  .  . .   q³⁄₂     .       .│
│X_Id^.  │   .   .      .  .  . .      .     .    q³⁄₂│
└────────┴────────────────────────────────────────────┘
```
"""
function XTable(uc::UnipotentClasses;q=Mvp(:q),classes=false)
# println("here uc=",uc)
  pieces=map(i->ICCTable(uc,i),eachindex(uc.springerseries))
# Note that c_ι=βᵤ+(rkss L_\CI)/2
  greenpieces=map((x,y)->map(x->x(q),x.scalar)*Diagonal(q.^x.dimBu)*
                  q^(length(y[:levi])//2),pieces,uc.springerseries)
  l=vcat(getproperty.(pieces,:locsys)...)
  p=inv(sortPerm(l))
  res=XTable(Dict(
    :scalar=>transpose(invpermute(cat(greenpieces...,dims=(1,2)),p)),
    :uc=>uc,
    :Y=>invpermute(cat(getproperty.(pieces,:L)...,dims=(1,2)),p,dims=(1,2)),
    :relgroups=>getindex.(uc.springerseries,:relgroup),
    :q=>q,
    :class=>classes))
  res.Y=map(x->x(q),res.Y)
  if classes
    res.scalar*=E(1)
    res.cardClass=zeros(eltype(res.scalar),length(l))*1//1
    res.centClass=zeros(eltype(res.scalar),length(l))*1//1
    res.classes=invpermute(l,p)
    for i in eachindex(uc.classes)
      Au=uc.classes[i].Au
      b=filter(j->res.classes[j][1]==i,eachindex(res.classes))
 #    println("i=",i," b=",b," Au=",Au)
      res.scalar[:,b]*=CharTable(Au).irr
      res.cardClass[b]=res.Y[[b[charinfo(Au).positionId]],b]*CharTable(Au).irr
      res.cardClass[b]=map((x,y)->x*y//length(Au),
                           res.cardClass[b],length.(conjugacy_classes(Au)))
      res.centClass[b]=generic_order(uc.spets)(q).//res.cardClass[b]
    end
    res.scalar=improve_type(res.scalar)
  else
    res.locsys=invpermute(l,p)
  end
  res
end
#@test haskey(t,:centClass)
#@test haskey(t,:cardClass)

function Base.show(io::IO, ::MIME"text/html", x::XTable)
  show(IOContext(io,:TeX=>true),"text/plain",x)
end

Base.show(io::IO,x::XTable)=print(io,"XTable(",x.uc,",q=",x.q,",classes=$(x.class))")

function Base.show(io::IO,::MIME"text/plain",x::XTable)
  printTeX(io,"Values of character sheaves \$\\tilde X_\\chi\$ on")
  row_labels=vcat(map(g->map(n->"X_{"*n*"}^{"*TeX(io,g)*"}",
                             charnames(TeX(io),g)),x.relgroups)...)
  rows_label="\\tilde X_\\chi|"
  tbl=x.scalar
  if x.class
    rows_label*="class"
    print(io," unipotent classes\n")
    col_labels=map(p->name(TeX(io;class=p[2]),x.uc.classes[p[1]]),x.classes)
    tbl=vcat(tbl,permutedims(x.centClass))
    push!(row_labels,"|C_{𝐆}(u)|")
    row_seps=[-1,0,size(tbl,1)-1,size(tbl,1)]
  else
    rows_label*=fromTeX(io,"\\phi")
    printTeX(io," local systems \$\\phi\$\n")
    col_labels=map(p->name(TeX(io,locsys=p[2]),x.uc.classes[p[1]]),x.locsys)
    row_seps=[-1,0,size(tbl,1)]
  end
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  dotzero=get(io,:dotzero,true)
  showtable(io,tbl;row_labels,col_labels,rows_label,dotzero,row_seps)
end

@GapObj struct GreenTable end

"""
`GreenTable(uc;classes=false)`

Keeping the same notations as in the description of 'XTable', this function
returns a table of the functions ``Q_{wF}``, attached to elements ``wF∈ W_𝐆
(𝐋)⋅F`` where ``W_𝐆 (𝐋)`` are the relative weyl groups attached to cuspidal
local  systems.  These  functions  are  defined  by ``Q_{wF}=∑_{u,ϕ} ϕ̃(wF)
X̃_{u,ϕ}``. An point to note is that in the principal Springer series, when
`𝐓`  is  a  maximal  torus,  the  function  ``Q_{wF}``  coincides  with the
Deligne-Lusztig  character ``R^𝐆  _{𝐓_W}(1)``. As  for 'XTable', by default
the  table  gives  the  values  of  the  functions  on  local  systems.  If
`classes=true`  is  given,  then  it  gives  the  values  of  the functions
``Q_{wF}`` on conjugacy classes.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> GreenTable(UnipotentClasses(W))
Values of Green functions Q_wF on local systems φ
┌──────────┬──────────────────────────────────────────────────────────┐
│Qᴵ_wF|φ   │        1     A₁       Ã₁ G₂(a₁)⁽¹¹¹⁾ G₂(a₁)⁽²¹⁾ G₂(a₁) G₂│
├──────────┼──────────────────────────────────────────────────────────┤
│Q_A₀^G₂   │  Φ₂²Φ₃Φ₆   Φ₂Φ₃ (2q+1)Φ₂           .          q   2q+1  1│
│Q_Ã₁^G₂   │-Φ₁Φ₂Φ₃Φ₆  -Φ₁Φ₃       Φ₂           .          q      1  1│
│Q_A₁^G₂   │-Φ₁Φ₂Φ₃Φ₆   Φ₂Φ₆      -Φ₁           .         -q      1  1│
│Q_G₂^G₂   │ Φ₁²Φ₂²Φ₃ -Φ₁Φ₂²    -Φ₁Φ₂           .         -q     Φ₂  1│
│Q_A₂^G₂   │ Φ₁²Φ₂²Φ₆  Φ₁²Φ₂    -Φ₁Φ₂           .          q    -Φ₁  1│
│Q_A₁+Ã₁^G₂│  Φ₁²Φ₃Φ₆  -Φ₁Φ₆ (2q-1)Φ₁           .         -q  -2q+1  1│
│Q_^.      │        .      .        .          q²          .      .  .│
└──────────┴──────────────────────────────────────────────────────────┘
```

The  functions ``Q_{wF}`` depend only on the conjugacy class of `wF`, so in
the  first column the indices of 'Q' are the names of the conjugacy classes
of ``W_𝐆(𝐋)``. The exponents are the names of the groups ``W_𝐆(𝐋)``.

```julia-repl
julia> GreenTable(UnipotentClasses(W);classes=true)
Values of Green functions Q_wF on unipotent classes
┌───────────┬────────────────────────────────────────────────────────┐
│Qᴵ_wF|class│        1     A₁       Ã₁ G₂(a₁) G₂(a₁)₍₂₁₎ G₂(a₁)₍₃₎ G₂│
├───────────┼────────────────────────────────────────────────────────┤
│Q_A₀^G₂    │  Φ₂²Φ₃Φ₆   Φ₂Φ₃ (2q+1)Φ₂   4q+1       2q+1        Φ₂  1│
│Q_Ã₁^G₂    │-Φ₁Φ₂Φ₃Φ₆  -Φ₁Φ₃       Φ₂   2q+1          1       -Φ₁  1│
│Q_A₁^G₂    │-Φ₁Φ₂Φ₃Φ₆   Φ₂Φ₆      -Φ₁  -2q+1          1        Φ₂  1│
│Q_G₂^G₂    │ Φ₁²Φ₂²Φ₃ -Φ₁Φ₂²    -Φ₁Φ₂    -Φ₁         Φ₂      2q+1  1│
│Q_A₂^G₂    │ Φ₁²Φ₂²Φ₆  Φ₁²Φ₂    -Φ₁Φ₂     Φ₂        -Φ₁     -2q+1  1│
│Q_A₁+Ã₁^G₂ │  Φ₁²Φ₃Φ₆  -Φ₁Φ₆ (2q-1)Φ₁  -4q+1      -2q+1       -Φ₁  1│
│Q_^.       │        .      .        .     q²        -q²        q²  .│
└───────────┴────────────────────────────────────────────────────────┘

julia> GreenTable(UnipotentClasses(rootdatum(:sl,4)))
Values of Green functions Q_wF on local systems φ
┌────────┬──────────────────────────────────────────────────────────────────┐
│Qᴵ_wF|φ │     1111          211 22⁽¹¹⁾       22   31 4 4^(ζ₄) 4⁽⁻¹⁾ 4^(ζ₄³)│
├────────┼──────────────────────────────────────────────────────────────────┤
│Q₁₁₁₁^A₃│  Φ₂²Φ₃Φ₄ (3q²+2q+1)Φ₂      . (2q+1)Φ₂ 3q+1 1      .     .       .│
│Q₂₁₁^A₃ │-Φ₁Φ₂Φ₃Φ₄   -q³+q²+q+1      .       Φ₂   Φ₂ 1      .     .       .│
│Q₂₂^A₃  │  Φ₁²Φ₃Φ₄        -Φ₁Φ₄      .  2q²-q+1  -Φ₁ 1      .     .       .│
│Q₃₁^A₃  │ Φ₁²Φ₂²Φ₄        -Φ₁Φ₂      .    -Φ₁Φ₂    1 1      .     .       .│
│Q₄^A₃   │ -Φ₁³Φ₂Φ₃        Φ₁²Φ₂      .      -Φ₁  -Φ₁ 1      .     .       .│
│Q₁₁^A₁  │        .            .   q²Φ₂        .    . .      .     q       .│
│Q₂^A₁   │        .            .  -q²Φ₁        .    . .      .     q       .│
│Q_^.    │        .            .      .        .    . .   q³⁄₂     .       .│
│Q_^.    │        .            .      .        .    . .      .     .    q³⁄₂│
└────────┴──────────────────────────────────────────────────────────────────┘
```
"""
function GreenTable(uc::UnipotentClasses;q=Mvp(:q),classes=false)
  t=GreenTable(XTable(uc;classes,q).prop)
  m=cat(map(g->transpose(CharTable(g).irr),t.relgroups)...;dims=(1,2))
  t.scalar=m*t.scalar
  t.indices=Vector{Int}[]
  i=0
  for g in t.relgroups
    n=nconjugacy_classes(g)
    push!(t.indices,(1:n).+i)
    i+=n
  end
  t
end

function Base.show(io::IO, ::MIME"text/html", x::GreenTable)
  show(IOContext(io,:TeX=>true),"text/plain",x)
end

Base.show(io::IO,x::GreenTable)=print(io,"GreenTable(",x.uc,",q=",x.q,")")

function Base.show(io::IO,::MIME"text/plain",x::GreenTable)
  printTeX(io,"Values of Green functions \$Q_{wF}\$ on")
  row_labels=vcat(map(x.relgroups) do g
     map(n->string("Q_{",n,"}^{",TeX(io,g),"}"),classnames(io,g))
    end...)
  rows_label="Q^I_{wF}|"
  if x.class
    println(io," unipotent classes")
    rows_label*="class"
    col_labels=map(p->name(TeX(io;class=p[2]),x.uc.classes[p[1]]),x.classes)
  else
    rows_label*="\\varphi"
    printTeX(io," local systems \$\\varphi\$\n")
    col_labels=map(p->name(TeX(io;locsys=p[2]),x.uc.classes[p[1]]),x.locsys)
  end
  tbl=x.scalar
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  showtable(io,tbl;row_labels,col_labels,rows_label,dotzero=get(io,:dotzero,true))
end

@GapObj struct ValuesTable end

"""
`UnipotentValues(uc,classes=false)`

This  function returns  a table  of the  values of  unipotent characters on
local systems (by default) or on unipotent classes (if `classes=true`).

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> UnipotentValues(UnipotentClasses(W);classes=true)
Values of unipotent characters for G₂ on unipotent classes
┌───────┬─────────────────────────────────────────────────────────────┐
│       │        1          A₁     Ã₁   G₂(a₁) G₂(a₁)₍₂₁₎ G₂(a₁)₍₃₎ G₂│
├───────┼─────────────────────────────────────────────────────────────┤
│φ₁‚₀   │        1           1      1        1          1         1  1│
│φ₁‚₆   │       q⁶           .      .        .          .         .  .│
│φ′₁‚₃  │  qΦ₃Φ₆/3    -qΦ₁Φ₂/3      q (q+5)q/3     -qΦ₁/3     qΦ₁/3  .│
│φ″₁‚₃  │  qΦ₃Φ₆/3  (2q²+1)q/3      .    qΦ₁/3     -qΦ₁/3  (q+2)q/3  .│
│φ₂‚₁   │ qΦ₂²Φ₃/6 (2q+1)qΦ₂/6  qΦ₂/2 (q+5)q/6     -qΦ₁/6     qΦ₁/6  .│
│φ₂‚₂   │ qΦ₂²Φ₆/2       qΦ₂/2  qΦ₂/2   -qΦ₁/2      qΦ₂/2    -qΦ₁/2  .│
│G₂[-1] │ qΦ₁²Φ₃/2      -qΦ₁/2 -qΦ₁/2   -qΦ₁/2      qΦ₂/2    -qΦ₁/2  .│
│G₂[1]  │ qΦ₁²Φ₆/6 (2q-1)qΦ₁/6 -qΦ₁/2 (q+5)q/6     -qΦ₁/6     qΦ₁/6  .│
│G₂[ζ₃] │qΦ₁²Φ₂²/3    -qΦ₁Φ₂/3      .    qΦ₁/3     -qΦ₁/3  (q+2)q/3  .│
│G₂[ζ₃²]│qΦ₁²Φ₂²/3    -qΦ₁Φ₂/3      .    qΦ₁/3     -qΦ₁/3  (q+2)q/3  .│
└───────┴─────────────────────────────────────────────────────────────┘

julia> UnipotentValues(UnipotentClasses(W,3);classes=true)
Values of unipotent characters for G₂ on unipotent classes
┌───────┬─────────────────────────────────────────────────────────────────────┐
│       │        1          A₁         Ã₁ G₂(a₁) G₂(a₁)₍₂₎    G₂       G₂_(ζ₃)│
├───────┼─────────────────────────────────────────────────────────────────────┤
│φ₁‚₀   │        1           1          1      1         1     1             1│
│φ₁‚₆   │       q⁶           .          .      .         .     .             .│
│φ′₁‚₃  │  qΦ₃Φ₆/3    -qΦ₁Φ₂/3        q/3  qΦ₂/3    -qΦ₁/3 -2q/3           q/3│
│φ″₁‚₃  │  qΦ₃Φ₆/3  (2q²+1)q/3        q/3  qΦ₂/3    -qΦ₁/3 -2q/3           q/3│
│φ₂‚₁   │ qΦ₂²Φ₃/6 (2q+1)qΦ₂/6  (3q+1)q/6  qΦ₂/6    -qΦ₁/6  2q/3          -q/3│
│φ₂‚₂   │ qΦ₂²Φ₆/2       qΦ₂/2      qΦ₂/2 -qΦ₁/2     qΦ₂/2     .             .│
│G₂[-1] │ qΦ₁²Φ₃/2      -qΦ₁/2     -qΦ₁/2 -qΦ₁/2     qΦ₂/2     .             .│
│G₂[1]  │ qΦ₁²Φ₆/6 (2q-1)qΦ₁/6 (-3q+1)q/6  qΦ₂/6    -qΦ₁/6  2q/3          -q/3│
│G₂[ζ₃] │qΦ₁²Φ₂²/3    -qΦ₁Φ₂/3        q/3  qΦ₂/3    -qΦ₁/3   q/3 (-ζ₃+2ζ₃²)q/3│
│G₂[ζ₃²]│qΦ₁²Φ₂²/3    -qΦ₁Φ₂/3        q/3  qΦ₂/3    -qΦ₁/3   q/3  (2ζ₃-ζ₃²)q/3│
└───────┴─────────────────────────────────────────────────────────────────────┘

┌───────┬─────────────────────────┐
│       │     G₂_(ζ₃²)       (Ã₁)₃│
├───────┼─────────────────────────┤
│φ₁‚₀   │            1           1│
│φ₁‚₆   │            .           .│
│φ′₁‚₃  │          q/3  (2q²+1)q/3│
│φ″₁‚₃  │          q/3    -qΦ₁Φ₂/3│
│φ₂‚₁   │         -q/3 (2q+1)qΦ₂/6│
│φ₂‚₂   │            .       qΦ₂/2│
│G₂[-1] │            .      -qΦ₁/2│
│G₂[1]  │         -q/3 (2q-1)qΦ₁/6│
│G₂[ζ₃] │ (2ζ₃-ζ₃²)q/3    -qΦ₁Φ₂/3│
│G₂[ζ₃²]│(-ζ₃+2ζ₃²)q/3    -qΦ₁Φ₂/3│
└───────┴─────────────────────────┘
```
"""
function UnipotentValues(uc;q=Mvp(:q),classes=false)
  t=ValuesTable(XTable(uc;classes,q).prop)
  uw=UnipotentCharacters(uc.spets)
  f=toL(fourier(uw))
  m=Vector{eltype(f[1])}[]
  for (i,ss) in pairs(uc.springerseries)
   if i==1 append!(m,f[charnumbers(uw.harishChandra[1])])
    elseif !haskey(ss,:hc) error("not implemented")
    elseif ss[:hc]==0 append!(m,map(i->zero(f[1]),eachindex(ss[:locsys])))
    else append!(m,f[charnumbers(uw.harishChandra[ss[:hc]])])
    end
  end
  t.scalar=transpose(toM(m))*t.scalar
  t
end

Base.show(io::IO,x::ValuesTable)=print(io,"UnipotentValues(",x.uc,",q=",x.q,")")

function Base.show(io::IO,::MIME"text/plain",x::ValuesTable)
  printTeX(io,"Values of unipotent characters for \$",x.uc.spets,"\$ on ")
  if x.class
    println(io,"unipotent classes")
    col_labels=map(p->name(TeX(io;class=p[2]),x.uc.classes[p[1]]),x.classes)
  else
    col_labels=map(p->name(TeX(io;locsys=p[2]),x.uc.classes[p[1]]),x.locsys)
    println(io,"local systems")
  end
  row_labels=charnames(TeX(io),UnipotentCharacters(x.uc.spets))
  tbl=improve_type(x.scalar)
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  showtable(io,tbl;row_labels,col_labels,dotzero=get(io,:dotzero,true))
end

@GapObj struct TwoVarGreenTable end

# two-variable green functions table.
# for now only implemented when W is split.
function TwoVarGreen(W,L)
  if !(W isa Spets) W=spets(W) end
  if !(L isa Spets) L=spets(L) end
  uG=UnipotentClasses(W)
  uL=UnipotentClasses(L)
  tG=GreenTable(uG;classes=true)
  tL=GreenTable(uL;classes=true)
  q=Mvp(:q)
  mm=map(eachindex(uL.springerseries))do i
    s=uL.springerseries[i]
    p=findfirst(S->S[:levi]==inclusion(L,s[:levi]) &&
                (isempty(S[:Z]) || S[:Z][1]==s[:Z][1]),uG.springerseries)
    if isnothing(p) error("not found ",s[:levi]) end
    RG=relative_coset(W,inclusion(L,s[:levi]))
    RLF=relative_coset(L,s[:levi])
    RL=Group(RLF)
    l=map(x->findfirst(==(x),Group(RG).relative_indices),
      inclusion(L,RL.relative_indices))
    if nothing in l error("not implemented") end
    RLF=subspets(RG,convert(Vector{Int},l),Group(RG).fromparent(L.phi))
    RL=Group(RLF)
    f=fusion_conjugacy_classes(RLF,RG)
    cl=classreps(RLF)
    d=map(cl)do w
      pw=word(RLF,w)
      if isempty(pw) pw=L.phi
      else pw=prod(Group(RG).toparent[pw])*L.phi
      end
      Lo=subspets(L,Int.(s[:levi]),pw/L.phi)
      r=map(last,filter(x->isone(first(x)),degrees(Lo)))
      prod(x->q-x,r)//length(centralizer(RL,w))
    end
    transpose(tL.scalar[tL.indices[i],:])*
    Diagonal(d)*conj(tG.scalar[tG.indices[p][f],:])
  end
  oL=generic_order(L,q)
  mm=improve_type(mm)
  mm=toM(map((x,y)->exactdiv(x*y,oL),eachrow(sum(mm)),tL.cardClass))
  res=TwoVarGreenTable(Dict(:W=>W,:L=>L,:scalar=>mm,:uL=>uL,:uG=>uG))
  res.classL=tL.classes
  res.classG=tG.classes
  res
end

function Base.show(io::IO, ::MIME"text/html", x::TwoVarGreenTable)
  show(IOContext(io,:TeX=>true),"text/plain",x)
end

Base.show(io::IO,x::TwoVarGreenTable)=print(io,"TwoVarGreenTable(",x.W,",",x.L,")")

function Base.show(io::IO,::MIME"text/plain",x::TwoVarGreenTable)
  print(io,"Values of two-variable Green functions for \$",x.W,"\$ and \$",
        x.L,"\$\n")
  rowLabels=map(x.classL)do p
    name(IOContext(io,:class=>p[2]),x.uL.classes[p[1]])
  end
  columnLabels=map(x.classG)do p
    name(IOContext(io,:class=>p[2]),x.uG.classes[p[1]])
  end
  tbl=improve_type(x.scalar)
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  showtable(io,tbl,row_labels=rowLabels,col_labels=columnLabels,rows_label="v/u")
end

"""
'special_pieces(<uc>)'

The special pieces form a partition of the unipotent variety of a reductive
group   `𝐆`   which   was   first   defined   in   [Spaltenstein1982  chap.
III](biblio.htm#spalt82)  as the fibres  of `d^2`, where  `d` is a "duality
map". Another definition is as the set of classes in the Zariski closure of
a  special class  and not  in the  Zariski closure  of any  smaller special
class,  where a  special class  is the  support of  the image  of a special
character by the Springer correspondence.

Each  piece is a union of unipotent  conjugacy classes so is represented in
Chevie  as a  list of  class numbers.  Thus the  list of  special pieces is
returned  as  a  list  of  lists  of  class  numbers. The list is sorted by
increasing  piece dimension, while each piece is sorted by decreasing class
dimension, so that the special class is listed first.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> special_pieces(UnipotentClasses(W))
3-element Vector{Vector{Int64}}:
 [1]
 [4, 3, 2]
 [5]

julia> special_pieces(UnipotentClasses(W,3))
3-element Vector{Vector{Int64}}:
 [1]
 [4, 3, 2, 6]
 [5]
```

The   example  above  shows  that  the  special  pieces  are  different  in
characteristic 3.
"""
function special_pieces(uc)
  W=uc.spets
  ch=charinfo(W)
  specialch=findall(iszero,ch.a-ch.b) # special characters of W
  specialc=first.(uc.springerseries[1][:locsys][specialch])
  sort!(specialc,by=c->-uc.classes[c].dimBu)
  m=transpose(incidence(uc.orderclasses))
  map(eachindex(specialc))do i
    p=m[specialc[i],:]
    for j in 1:i-1 p.&=.!m[specialc[j],:] end
    p=eachindex(p)[p]
    sort!(p,by=c->uc.classes[c].dimBu)
    p
  end
end
"""

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%'UnipotentValues(<W>,<w>)'
"""

end
