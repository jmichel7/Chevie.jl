"""
Let  ``𝐆 `` be a connected reductive group over an algebraic closure of the
finite  field ``𝔽_q``,  defined over  ``𝔽_q``, with corresponding Frobenius
automorphism  ``F``. We  want to  study the  irreducible characters  of ``𝐆
^F``.  More generally we consider ``𝐆 ^F`` where ``F`` is a Frobenius root,
an  isogeny of  ``𝐆 ``  such that  a power  is a Frobenius (this covers the
Suzuki and Ree groups).

If  ``𝐓`` is an  ``F``-stable maximal torus  of ``𝐆``, and  ``𝐁`` is a (not
necessarily  ``F``-stable) Borel  subgroup containing  ``𝐓``, we define the
*Deligne-Lusztig*  variety ``X_𝐁=\\{g𝐁∈𝐆/𝐁 ∣ g𝐁∩F(g𝐁)≠∅ \\}``. This variety
affords  a  natural  action  of  ``𝐆^F``  on the left, so the corresponding
*Deligne-Lusztig  virtual module*  given by  the `ℓ`-adic  cohomology with
compact  support  ``H^*_c(X_𝐁):=∑ᵢ(-1)ⁱHⁱ_c(X_𝐁,ℚ̄_ℓ)``  has  an  action of
``𝐆^F``  on  the  right.  The  (virtual)  character  of  this module is the
*Deligne-Lusztig* character ``R_𝐓^𝐆(1)``; the notation reflects the theorem
that  this character does not depend on the choice of ``𝐁``. This character
can  be parameterized by an ``F``-conjugacy class of ``W``: if ``𝐓₀⊂𝐁₀`` is
an ``F``-stable pair, there is an unique ``w∈ W=N_𝐆 (𝐓₀)/𝐓₀`` such that the
triple  ``(𝐓,𝐁,F)``  is  ``𝐆``-conjugate  to  ``(𝐓₀,𝐁₀,wF)``.  We will thus
denote  ``R_w``  for  ``R_𝐓^𝐆(1)``;  this  character  depends  only  on the
``F``-class of ``w``.

The  *unipotent characters* of ``𝐆^F`` are the irreducible constituents of
the  ``R_w``. In a similar way that the Jordan decomposition shows that the
unipotent classes are a building block for describing the conjugacy classes
of  a  reductive  group,  Lusztig  has  defined  a  Jordan decomposition of
characters  where  the  unipotent  characters  are  the building block. The
unipotent  characters are parameterized by  combinatorial data that Lusztig
has  defined just  from the  coset ``Wφ``,  where `φ`  is the  finite order
automorphism  of ``X(𝐓₀)``  such that  ``F=qφ``. Thus,  from our viewpoint,
unipotent  characters  are  objects  combinatorially  attached to a Coxeter
coset.

A  subset  of  the  unipotent  characters, the *principal series* unipotent
characters,  can  be  described  in  a  more  elementary  way. They are the
constituents  of  ``R₁``,  or  equivalently  the  characters of the virtual
module   ``H^*_c(X_{𝐁₀})``,  where  ``X_{𝐁₀}``   is  the  discrete  variety
``(𝐆/𝐁₀)^F``;   this   virtual   module   reduces   to  the  actual  module
``ℚ̄_ℓ[(𝐆/𝐁₀)^F]``.  Thus the  Deligne-Lusztig induction ``R₁=R_{𝐓₀}^𝐆(1)``
reduces  to Harish-Chandra induction, defined  as follows: let ``𝐏=𝐔⋊𝐋`` be
an ``F``-stable Levi decomposition of an ``F``-stable parabolic subgroup of
``𝐆``.  Then the *Harish-Chandra* induced ``R_𝐋^𝐆`` of a character ``χ`` of
``𝐋^F``  is the character ``Ind_{𝐏^F}^{𝐆 ^F}χ̃``,  where ``χ̃`` is the lift
to  ``𝐏^F``  of  ``χ``  via  the  quotient  ``𝐏^F/𝐔^F=𝐋^F``; Harish-Chandra
induction  is a  particular case  of *Lusztig  induction*, which is defined
when  ``𝐏``  is  not  ``F``-stable  using  the  variety  ``X_𝐔=\\{ g𝐔∈𝐆/𝐔 ∣
g𝐔∩F(g𝐔)≠∅\\}``, and gives for an ``𝐋^F``-module a virtual ``𝐆 ^F``-module.
Like  ordinary induction, these  functors have adjoint  functors going from
representations    of   ``𝐆^F``    to   representations    (resp.   virtual
representations)   of  ``𝐋^F``  called  Harish-Chandra  restriction  (resp.
Lusztig restriction).

The  commuting algebra of ``𝐆^F``-endomorphisms of ``R₁=R_{𝐓₀}^𝐆(1)`` is an
Iwahori-Hecke  algebra for ``W^φ``, with parameters some powers of `q`; the
parameters  are  all  equal  to  `q`  when ``W^φ=W``. Thus principal series
unipotent characters are parametrized by characters of ``W^φ``.

To   understand  the  decomposition  of  more  general  ``R_w``,  and  thus
parameterize unipotent characters, is is useful to introduce another set of
class  functions which are  parameterized by irreducible  characters of the
coset  ``Wφ``.  If  ``χ``  is  such  a  character, we define the associated
*almost  character* by:  ``Rᵪ=|W|⁻¹∑_{w∈ W}χ(wφ)  R_w``. The  name reflects
that these class function are close to irreducible characters. They satisfy
``⟨Rᵪ,  R_ψ⟩_{𝐆^F}=δ_{χ,ψ}``;  for  the  linear  and unitary group they are
actually  unipotent characters (up to sign in the latter case). They are in
general the sum (with rational coefficients) of a small number of unipotent
characters in the same *Lusztig family*, see [`Families`](@ref). The degree
of  ``Rᵪ``  is  a  polynomial  in  ``q``  equal  to  the fake degree of the
character ``χ`` of ``Wφ`` (see [`fakedegree`](@ref)).

We   now  describe  the  parameterization   of  unipotent  characters  when
``W^φ=W``,  in  which  case  the  coset  ``Wφ``  identifies with ``W`` (the
general  situation is  similar but  a bit  more difficult to describe). The
(rectangular) matrix of scalar products ``⟨ρ, Rᵪ⟩_{𝐆 ^F}``, when characters
of  ``W``  and  unipotent  characters  are  arranged in the right order, is
block-diagonal   with  rather  small  blocks   which  are  called  *Lusztig
families*.

For  the characters  of ``W``  a family  `𝓕` corresponds  to a block of the
Hecke  algebra  over  a  ring  called  the  Rouquier  ring.  To `𝓕` Lusztig
associates  a small group ``Γ`` (equal  to ``(ℤ/2)ⁿ`` for rather small `n`,
or ``𝔖ᵢ`` for ``i≤5``) such that the unipotent characters in the family are
parameterized  by the  characters of  the Drinfed  double of  `Γ``, that is
pairs  ``(x,θ)`` taken up to ``Γ``-conjugacy, where ``x∈Γ`` and ``θ`` is an
irreducible   character  of  ``C_Γ(x)``.  Further,   the  elements  of  `𝓕`
themselves are parameterized by a subset of such pairs, and Lusztig defines
a  pairing  between  such  pairs  which  computes  the scalar product ``⟨ρ,
Rᵪ⟩_{𝐆^F}``,  called  the  *Lusztig  Fourier  matrix*. For more details see
[`drinfeld_double`](@ref).

A  second parameterization  of unipotent  character is  via *Harish-Chandra
series*.  A character is called *cuspidal* if all its proper Harish-Chandra
restrictions  vanish. There are few  cuspidal unipotent characters (none in
``GLₙ``  for  ``n>1``,  and  at  most  one  in other classical groups). The
``𝐆^F``-endomorphism algebra of an Harish-Chandra induced
``R_{𝐋^F}^{𝐆^F}λ``,  where  ``λ``  is  a  cuspidal unipotent character is a
Hecke algebra associated to the group ``W_{𝐆^F}(𝐋^F):=N_{𝐆^F}(𝐋)/𝐋``, which
turns  out  to  be  a  Coxeter  group.  Thus another parameterization is by
triples  ``(𝐋,λ,φ)``,  where  ``λ``  is  a  cuspidal unipotent character of
``𝐋^F``  and  ``φ``  is  an  irreducible  character of the *relative group*
``W_{𝐆^F}(𝐋^F)``.  Such characters are said to belong to the Harish-Chandra
series determined by ``(𝐋,λ)``.

A  final  piece  of  information  attached  to  unipotent characters is the
*eigenvalues of Frobenius*. Let ``Fᵟ`` be the smallest power of the isogeny
``F``  which  is  a  split  Frobenius  (that  is, ``Fᵟ`` is a Frobenius and
``φᵟ=1``). Then ``Fᵟ`` acts naturally on Deligne-Lusztig varieties and thus
on  the  corresponding  virtual  modules,  and  commutes  to  the action of
``𝐆^F``;  thus for  a given  unipotent character  ``ρ``, a submodule of the
virtual  module which  affords ``ρ``  affords a  single eigenvalue ``μ`` of
``Fᵟ``.  [lu78;  3.9](@cite)  and  [dm85;  II,  2.3](@cite)  show that this
eigenvalue  is of the form ``qᵃᵟλᵨ`` where ``2a∈ℤ`` and ``λᵨ`` is a root of
unity, and where the parity of `2a` and ``λᵨ`` depend only on ``ρ`` and not
the  considered module. This ``λᵨ`` is  called the *eigenvalue of Frobenius
attached  to ``ρ``*. Unipotent characters in the Harish-Chandra series of a
pair ``(𝐋,λ)`` have the same eigenvalue of Frobenius as ``λ``.

`Chevie`   contains  tables  of  all  this  information,  and  can  compute
Harish-Chandra  and Lusztig  induction of  unipotent characters  and almost
characters. We illustrate this on some examples:

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> uc=UnipotentCharacters(W)
UnipotentCharacters(G₂)
┌───────┬────────────────────────────────────────────────────┐
│γ      │n₀    Deg(γ)  Feg              Symbol Fr(γ)    label│
├───────┼────────────────────────────────────────────────────┤
│φ₁‚₀   │ 1         1    1       (0,0,0,0,0,2)     1         │
│φ₁‚₆   │ 2        q⁶   q⁶ (01,01,01,01,01,12)     1         │
│φ′₁‚₃  │ 3   qΦ₃Φ₆/3   q³            (0,0,1+)     1    (1,ρ)│
│φ″₁‚₃  │ 4   qΦ₃Φ₆/3   q³            (0,0,1-)     1   (g₃,1)│
│φ₂‚₁   │ 5  qΦ₂²Φ₃/6  qΦ₈       (0,0,0,0,1,1)     1    (1,1)│
│φ₂‚₂   │ 6  qΦ₂²Φ₆/2 q²Φ₄       (0,0,0,1,0,1)     1   (g₂,1)│
│G₂[-1] │ 7  qΦ₁²Φ₃/2    0       (01,0,01,,0,)    -1   (g₂,ε)│
│G₂[1]  │ 8  qΦ₁²Φ₆/6    0       (01,01,0,,,0)     1    (1,ε)│
│G₂[ζ₃] │ 9 qΦ₁²Φ₂²/3    0       (01,0,0,01,,)    ζ₃  (g₃,ζ₃)│
│G₂[ζ₃²]│10 qΦ₁²Φ₂²/3    0       (01,01,,0,0,)   ζ₃² (g₃,ζ₃²)│
└───────┴────────────────────────────────────────────────────┘
```
The  first column gives  the name of  the unipotent character, derived from
its  Harish-Chandra  classification;  the  first  6  characters  are in the
principal  series  so  are  named  by  characters  of  `W`.  The last 4 are
cuspidal,  and named by the corresponding eigenvalue of Frobenius, which is
displayed  in the fourth  column. For classical  groups, the Harish-Chandra
data can be synthesized combinatorially to give a *symbol*.

The  first two characters are  each in a Lusztig  family by themselves. The
last  eight are in a family associated to the group `Γ=𝔖₃`: the last column
shows  the parameters  `(x,θ)`. The  third column  shows the  degree of the
unipotent characters, which is transformed by the Lusztig Fourier matrix of
the  fourth  column,  which  gives  the  degree of the corresponding almost
character,  or equivalently the fake  degree of the corresponding character
of `W` (extended by `0` outside the principal series).

One  can get  more information  on the  Lusztig Fourier  matrix of  the big
family by asking

```julia-repl
julia> uc.families[1]
Family(D(𝔖 ₃)) Drinfeld double of 𝔖 ₃, Lusztig′s version
┌────────┬───────────────────────────────────────────────────────┐
│label   │no eigen                                               │
├────────┼───────────────────────────────────────────────────────┤
│(1,1)*  │ 5     1 1//6  1//2  1//3  1//3  1//6  1//2  1//3  1//3│
│(g₂,1)  │ 6     1 1//2  1//2     .     . -1//2 -1//2     .     .│
│(g₃,1)  │ 4     1 1//3     .  2//3 -1//3  1//3     . -1//3 -1//3│
│(1,ρ)   │ 3     1 1//3     . -1//3  2//3  1//3     . -1//3 -1//3│
│(1,ε)-e │ 8     1 1//6 -1//2  1//3  1//3  1//6 -1//2  1//3  1//3│
│(g₂,ε)  │ 7    -1 1//2 -1//2     .     . -1//2  1//2     .     .│
│(g₃,ζ₃) │ 9    ζ₃ 1//3     . -1//3 -1//3  1//3     .  2//3 -1//3│
│(g₃,ζ₃²)│10   ζ₃² 1//3     . -1//3 -1//3  1//3     . -1//3  2//3│
└────────┴───────────────────────────────────────────────────────┘
```
We  can also do computations with  individual unipotent characters. Here we
construct  the Coxeter torus, and then the identity character of this torus
as a unipotent character.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> T=spets(reflection_subgroup(W,Int[]),W(1,2))
G₂₍₎=Φ₆

julia> u=unipotent_character(T,1)
[G₂₍₎=Φ₆]:<Id>
```

To construct `T` one could equivalently do
```julia-repl
julia> T=torus(W,position_regular_class(W,6))
G₂₍₎=Φ₆
```
Here  are two ways to construct the Deligne-Lusztig character associated to
the Coxeter torus:
```julia-repl
julia> lusztig_induce(W,u)
[G₂]:<φ₁‚₀>+<φ₁‚₆>-<φ₂‚₁>+<G₂[-1]>+<G₂[ζ₃]>+<G₂[ζ₃²]>

julia> v=deligne_lusztig_character(W,[1,2])
[G₂]:<φ₁‚₀>+<φ₁‚₆>-<φ₂‚₁>+<G₂[-1]>+<G₂[ζ₃]>+<G₂[ζ₃²]>

julia> degree(v)
Pol{Int64}: q⁶+q⁵-q⁴-2q³-q²+q+1

julia> scalar_product(v,v)
6
```
The  last two lines ask  for the degree of  the virtual character `v`, then
for the scalar product of `v` with itself.

Finally  we mention  that Chevie  can also  provide unipotent characters of
Spetses, as defined in [bmm14](@cite). An example:

```julia-repl
julia> UnipotentCharacters(complex_reflection_group(4))
UnipotentCharacters(G₄)
┌─────┬────────────────────────────────────────┐
│γ    │n₀            Deg(γ)    Feg Fr(γ)  label│
├─────┼────────────────────────────────────────┤
│φ₁‚₀ │ 1                 1      1     1       │
│φ₁‚₄ │ 2  -√-3q⁴Φ″₃Φ₄Φ″₆/6     q⁴     1   1∧ζ₆│
│φ₁‚₈ │ 3   √-3q⁴Φ′₃Φ₄Φ′₆/6     q⁸     1 -1∧ζ₃²│
│φ₂‚₅ │ 4         q⁴Φ₂²Φ₆/2   q⁵Φ₄     1  1∧ζ₃²│
│φ₂‚₃ │ 5 -ζ₃√-3qΦ″₃Φ₄Φ′₆/3   q³Φ₄     1  1∧ζ₃²│
│φ₂‚₁ │ 6 ζ₃²√-3qΦ′₃Φ₄Φ″₆/3    qΦ₄     1   1∧ζ₃│
│φ₃‚₂ │ 7            q²Φ₃Φ₆ q²Φ₃Φ₆     1       │
│Z₃:2 │ 8     -√-3qΦ₁Φ₂Φ₄/3      0   ζ₃² ζ₃∧ζ₃²│
│Z₃:11│ 9    -√-3q⁴Φ₁Φ₂Φ₄/3      0   ζ₃² ζ₃∧ζ₆⁵│
│G₄   │10        -q⁴Φ₁²Φ₃/2      0    -1  ζ₆∧-1│
└─────┴────────────────────────────────────────┘
```
"""
module Uch

using ..Chevie

export UnipotentCharacters, UniChar, unichar, 
       unipotent_character, almostchar, almost_character, dlchar, dlCharTable,
       deligne_lusztigCharTable,
       deligne_lusztig_character, dllefschetz, deligne_lusztig_lefschetz,
       lusztig_induce, lusztig_restrict, 
       hc_induce, hc_restrict, 
       cuspidal, cuspidal_data,
       CycPoldegrees, on_unipotents, almostcharnames

@GapObj struct UnipotentCharacters
  harishChandra::Vector{Dict{Symbol,Any}}
  almostHarishChandra::Vector{Dict{Symbol,Any}}
  families::Vector{Family}
end

function params(sers)
  chh=map(ser->charinfo(ser[:relativeType]),sers)
  res=fill(Pair("foo",[]),sum(x->length(x.charnames),chh))
  for (ch,ser) in zip(chh,sers)
#   t=ser[:relativeType]
#   t.rank=haskey(t,:orbit) ? t.orbit[1].rank : t.rank
    res[charnumbers(ser)]=Pair.(ser[:cuspidalName],ch.charparams)
  end
  res
end

function SerNames(io::IO,sers)
  res=fill("",sum(x->length(charnumbers(x)),sers))
  for ser in sers
    tt=ser[:relativeType]
    if !(tt isa Vector) tt=[tt] end
    n=fromTeX(io,ser[:cuspidalName])
    if isempty(tt) nn=[n]
    else nn=charnames(io,tt)
      if !isempty(ser[:levi]) nn=string.(n,":",nn) end
    end
    res[charnumbers(ser)]=nn
  end
  res
end

"""
charnames(uc;options...) charnames(io::IO,uc)

returns  the list of character names for the unipotent characters uc. The
optional  options  are  IOContext  attributes  which can give alternative
names  in certain cases,  or a different  formatting of names in general.
They can be specified by giving an IO as argument.
"""
function Chars.charnames(io::IO,uc::UnipotentCharacters)
  if get(io,:cyclicparam,false) && haskey(uc,:cyclicparam)
    map(uc.cyclicparam)do x
      if length(x[1])==1 "Id"
      else fromTeX(io,string("\\rho_{",x[1][1],",",x[1][2],"}"))
      end
    end
  else SerNames(io,uc.harishChandra)
  end
end

almostcharnames(io::IO,uc::UnipotentCharacters)=SerNames(io,uc.almostHarishChandra)

function UnipotentCharacters(t::TypeIrred)
  uc=chevieget(t,:UnipotentCharacters)
  if uc===nothing
    InfoChevie("# $t is not spetsial!!")
    return
  end
  uc=copy(uc)
  uc[:harishChandra]=copy.(uc[:harishChandra])
  if haskey(uc,:almostHarishChandra)
    uc[:almostHarishChandra]=copy.(uc[:almostHarishChandra])
  end
  uc[:charParams]=params(uc[:harishChandra])
  if !haskey(uc,:charSymbols) uc[:charSymbols]=uc[:charParams] end
  # adjust things for descent of scalars
  # we would like to adjust indices so they fit with those stored in t
  # but we cannot when indices mention non-generating reflections!
  a=length(t.orbit)
  if a>1
    if haskey(uc,:a) uc[:a]*=a end
    if haskey(uc,:A) uc[:A]*=a end
    for s in uc[:harishChandra]
      s[:parameterExponents]*=a
#     s[:eigenvalue]^=a
      if s[:cuspidalName]=="" s[:cuspidalName]="Id" end
#     s[:cuspidalName]=join(fill(s[:cuspidalName],a),"\\otimes ")
    end
  end

  if !haskey(uc,:almostHarishChandra)
    uc[:almostHarishChandra]=map(uc[:harishChandra])do s
      res=Dict{Symbol,Any}()
      for f in [:levi, :cuspidalName, :eigenvalue, :charNumbers] res[f]=s[f] end
      res[:relativeType]=TypeIrred(;orbit=
      map(eachindex(t.orbit))do i
          r=copy(s[:relativeType])
          r.indices=r.indices.+(i-1)*rank(t.orbit[1])
          r
        end, twist=Perm())
      if haskey(s[:relativeType],:twist) && s[:relativeType][:twist]!=Perm()
        error()
      end
      if !isone(t.twist)
        a=t.orbit[1].indices[s.relativeType[:indices]]
        res[:relativeType][:twist]=prod(Perm.(a,a.^t.twist))
      end
      res[:levi]=vcat(map(eachindex(t.orbit))do i
       res[:levi].+(i-1)*rank(t.orbit[1])
      end...)
      res
    end
  else
    for s in uc[:almostHarishChandra]
      if !haskey(s[:relativeType],:orbit)
        s[:relativeType]=TypeIrred(;orbit=[s[:relativeType]],twist=Perm())
      end
    end
  end
  if !haskey(uc,:almostCharSymbols) uc[:almostCharSymbols]=uc[:charSymbols] end
  uc[:almostCharParams]=params(uc[:almostHarishChandra])
  uc[:type]=t
  uch=UnipotentCharacters(uc[:harishChandra],uc[:almostHarishChandra],
                          Family.(copy(uc[:families])),uc)
  delete!(uc,:families)
  delete!(uc,:harishChandra)
  delete!(uc,:almostHarishChandra)
  uch
end

function UnipotentCharacters(W::Group)
  get!(W,:UnipotentCharacters) do
    UnipotentCharacters(spets(W))
  end
end

function CartesianSeries(sers)
  ser=Dict{Symbol,Any}()
  ser[:levi]=reduce(vcat,getindex.(sers,:levi))
  ser[:relativeType]=filter(x->rank(x)!=0,getindex.(sers,:relativeType))
  if haskey(sers[1],:eigenvalue)
    ser[:eigenvalue]=prod(getindex.(sers,:eigenvalue))
  end
  if any(x->haskey(x,:qEigen),sers)
    ser[:qEigen]=sum(sers)do x
      haskey(x,:qEigen) ? x[:qEigen] : 0
    end
  else
    ser[:qEigen]=0
  end
  if all(haskey.(sers,:parameterExponents))
    ser[:parameterExponents]=vcat(getindex.(sers,:parameterExponents)...)
  end
  ser[:charNumbers]=tcartesian(charnumbers.(sers)...)
  ser[:cuspidalName]=join(map(x->x[:cuspidalName]=="" ? "Id" :
                                 x[:cuspidalName], sers),"\\otimes ")
  ser
end

"""
`UnipotentCharacters(W)`

`W`  should be a Coxeter group, a  Coxeter Coset or a Spetses. The function
gives  back a record containing  information about the unipotent characters
of the associated algebraic group (or Spetses). This contains the following
fields:

`.harishChandra`:  information  about  Harish-Chandra  series  of  unipotent
characters.  This is itself a list of records, one for each pair `(𝐋,λ)` of
a  Levi  of  an  `F`-stable  parabolic  subgroup  and  a cuspidal unipotent
character of ``𝐋^F``. These records themselves have the following fields:

`:levi`: a list 'l' such that `𝐋` corresponds to 'reflection_subgroup(W,l)'.

`:cuspidalName`: the name of the unipotent cuspidal character `lambda`.

`:eigenvalue`: the eigenvalue of Frobenius for `λ`.

`:relativeType`: the reflection type of ``W_𝐆 (𝐋)``;

`:parameterExponents`:  the ``𝐆 ^F``-endomorphism  algebra of ``R_𝐋^𝐆 (λ)``
is  a  Hecke  algebra  for  ``W_𝐆  (𝐋)``  with  some parameters of the form
``q^{a_s}``. This holds the list of exponents ``a_s``.

`:charNumbers`:  the  indices  of  the  unipotent  characters indexed by the
irreducible characters of ``W_𝐆 (𝐋)``.

`.almostHarishChandra`:   information   about   Harish-Chandra   series  of
unipotent  character sheaves.  This is  identical to  ̀harishChandra` for a
split  reductive group,  and reflects  the situation  for the corresponding
split group for a nonsplit group.

`.families`:  information  about  Lusztig  families of unipotent characters.
This  is itself a list  of records, one for  each family. These records are
described in the section about families below.

the following information is computed on demand from
`uc=UnipotentCharacters(W)`:

`spets(uc)`: the reductive group `W`.

```julia-repl
julia> W=coxgroup(:Bsym,2)
Bsym₂

julia> WF=spets(W,Perm(1,2))
²Bsym₂

julia> uc=UnipotentCharacters(WF)
UnipotentCharacters(²Bsym₂)
┌────────┬─────────────────────────────────────────────────────┐
│γ       │n₀ almostch    Deg(γ)   Feg        Symbol Fr(γ) label│
├────────┼─────────────────────────────────────────────────────┤
│2       │ 1       2.         1     1     (02,,0,0)     1      │
│11      │ 2      .11        q⁴    q⁴ (012,1,01,01)     1      │
│²B₂[1,3]│ 3      1.1 √2qΦ₁Φ₂/2 qΦ₁Φ₂     (01,,1,0)   ζ₈³     1│
│²B₂[1,5]│ 4       B₂ √2qΦ₁Φ₂/2     0     (01,,0,1)   ζ₈⁵     2│
└────────┴─────────────────────────────────────────────────────┘

julia> uc.families
3-element Vector{Family}:
 Family(C₁,[1])
 Family(C₁,[2])
 Family(?4,[3, 4])

julia> uc.families[3]
Family(?4) 
┌─────┬───────────────────┐
│label│no eigen   1*     2│
├─────┼───────────────────┤
│1*   │ 3   ζ₈³ √2/2 -√2/2│
│2    │ 4   ζ₈⁵ √2/2  √2/2│
└─────┴───────────────────┘
```

`charnames(uc)`:  the list of names of the unipotent characters.  Using
   appropriate keywords, one can control the display in various ways.

```julia-repl
julia> uc=UnipotentCharacters(coxgroup(:G,2));

julia> charnames(uc;limit=true)
10-element Vector{String}:
 "φ₁‚₀"
 "φ₁‚₆"
 "φ′₁‚₃"
 "φ″₁‚₃"
 "φ₂‚₁"
 "φ₂‚₂"
 "G₂[-1]"
 "G₂[1]"
 "G₂[ζ₃]"
 "G₂[ζ₃²]"

julia> charnames(uc;TeX=true)
10-element Vector{String}:
 "\\phi_{1,0}"
 "\\phi_{1,6}"
 "\\phi_{1,3}'"
 "\\phi_{1,3}''"
 "\\phi_{2,1}"
 "\\phi_{2,2}"
 "G_2[-1]"
 "G_2[1]"
 "G_2[\\zeta_3]"
 "G_2[\\zeta_3^2]"
```

One  can control  the display  of unipotent  characters in  various ways by
`IOContext` properties. In the display, the row labels are the names of the
unipotent characters. The possible columns are numbered as follows:

  1. The index of the character in the list of unipotent characters.
  2. The degree of the unipotent character.
  3. The degree of the corresponding almost character.
  4. for imprimitive groups, the symbol attached to the unipotent character.
  5. The eigenvalue of Frobenius attached to the unipotent character.
  6. The parameter the character has in its Lusztig family.

Which  columns  are  displayed  can  be  controlled by the property `:cols`
(default [2,3,5,6] and 4 when applicable).

In  addition if  ':byfamily=true', the  characters are  displayed family by
family  instead  of  in  index  order.  Finally,  the properties `rows` and
`columnrepartition`  of  `format`  can  be  set,  giving more tuning of the
table.

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> uc=UnipotentCharacters(W)
UnipotentCharacters(B₂)
┌───┬──────────────────────────────────┐
│γ  │n₀ Deg(γ) Feg   Symbol Fr(γ) label│
├───┼──────────────────────────────────┤
│11.│ 1  qΦ₄/2  q²   (12,0)     1   +,-│
│1.1│ 2 qΦ₂²/2 qΦ₄   (02,1)     1   +,+│
│.11│ 3     q⁴  q⁴ (012,12)     1      │
│2. │ 4      1   1     (2,)     1      │
│.2 │ 5  qΦ₄/2  q²   (01,2)     1   -,+│
│B₂ │ 6 qΦ₁²/2   0   (012,)    -1   -,-│
└───┴──────────────────────────────────┘
```

```julia-rep1
julia> xdisplay(uc;byfamily=true)
┌────┬──────────────────────────────────┐
│γ   │n₀ Deg(γ) Feg   Symbol Fr(γ) label│
├────┼──────────────────────────────────┤
│11. │ 1  qΦ₄/2  q²   (12,0)     1   +,-│
│1.1ˢ│ 2 qΦ₂²/2 qΦ₄   (02,1)     1   +,+│
│.2  │ 5  qΦ₄/2  q²   (01,2)     1   -,+│
│B₂  │ 6 qΦ₁²/2   0   (012,)    -1   -,-│
├────┼──────────────────────────────────┤
│2.  │ 4      1   1     (2,)     1      │
├────┼──────────────────────────────────┤
│.11 │ 3     q⁴  q⁴ (012,12)     1      │
└────┴──────────────────────────────────┘

julia> xdisplay(uc;cols=[1,4])
UnipotentCharacters(B₂)
┌───┬───────────┐
│γ  │n₀   Symbol│
├───┼───────────┤
│11.│ 1   (12,0)│
│1.1│ 2   (02,1)│
│.11│ 3 (012,12)│
│2. │ 4     (2,)│
│.2 │ 5   (01,2)│
│B₂ │ 6   (012,)│
└───┴───────────┘
```
"""
function UnipotentCharacters(WF::Spets)
  get!(WF,:UnipotentCharacters) do
  tt=refltype(WF)
  if isempty(tt) # UnipotentCharacters(coxgroup())
    return UnipotentCharacters(
      [Dict(:relativeType=>TypeIrred[],
	    :levi=>Int[], :parameterExponents=>Int[],
	    :cuspidalName=>"Id", :eigenvalue=>1, :charNumbers =>[ 1 ])],
      [Dict(:relativeType=>TypeIrred[],
	    :levi=>Int[], :parameterExponents=>Int[],
	    :cuspidalName=>"Id", :eigenvalue=>1, :charNumbers =>[ 1 ])],
     [Family("C1",[1])],
     Dict(:charParams=>[["",[1]]], :charSymbols=>[[CharSymbol([[1]])]],
      :size=>1, :a =>[0], :A =>[0], :spets=>WF))
  end

  W=WF.W
  simp=map(tt) do t
# adjust indices of Levis, almostLevis, relativetypes so they agree with
# Parent(Group(WF))
    uc=UnipotentCharacters(t)
    if isnothing(uc) return end
 #  H=map(x->reflection_subgroup(W,x.indices[1:x.rank]),t.orbit)
    i=indices(t)
    H=reflection_subgroup(W,sort(i))
    p=mappingPerm(i)^mappingPerm(sort(i),eachindex(i))
#   @show t,i,H,p
    for s in uc.harishChandra
      s[:levi]=inclusion(H,W,s[:levi].^p)
      s[:relativeType]=copy(s[:relativeType])
      s[:relativeType].indices=inclusion(H,W,s[:relativeType].indices.^p)
    end
    for s in uc.almostHarishChandra
      s[:levi]=inclusion(H,W,s[:levi].^p)
      s[:relativeType]=copy(s[:relativeType])
      s[:relativeType].orbit=vcat(
        map(s[:relativeType].orbit)do r
	  r=copy(r)
          r.indices=inclusion(H,W,r.indices.^p)
	  r
        end...)
      s[:relativeType].twist^=mappingPerm(1:length(inclusion(H)),inclusion(H))
    end

    for f in uc.families
      f.fourierMat=fourier(f)
      if !haskey(f,:charLabels) f.charLabels=string.(1:length(f)) end
    end
    uc
  end

  # "Kronecker product" of records in simp:
  r=simp[1]
  if isnothing(r) return end
  f=keys(r.prop)
  extra=Dict{Symbol,Any}()
  for a in f
    if a==:type continue end
    if length(simp)==1
      extra[a]=map(x->[x],getproperty(r,a))
    elseif all(x->haskey(x,a),simp)
      extra[a]=cartesian(getproperty.(simp,a)...)
    end
  end

  extra[:size]=length(extra[:charParams])

  hh=CartesianSeries.(cartesian(map(x->x.harishChandra,simp)...))
  ah=CartesianSeries.(cartesian(map(x->x.almostHarishChandra,simp)...))
  # finally the new 'charNumbers' lists
  ll=length.(simp)
  for s in hh s[:charNumbers]=cart2lin.(Ref(ll),s[:charNumbers]) end
  for s in ah s[:charNumbers]=cart2lin.(Ref(ll),s[:charNumbers]) end

  if length(tt)==1
    ff=r.families
  else
    ff=Family.(prod.(cartesian(map(x->x.families,simp)...)))
    for f in ff f.charNumbers=cart2lin.(Ref(ll),f.charNumbers) end
  end

  for a in [:a, :A]
    if haskey(extra,a) extra[a]=sum.(extra[a]) end
  end

  extra[:spets]=WF
  UnipotentCharacters(hh,ah,ff,extra)
  end
end

function Chars.fakedegrees(uc::UnipotentCharacters,q=Pol())
  f=fakedegrees(spets(uc),q)
  fd=fill(zero(f[1]),length(uc))
  fd[uc.almostHarishChandra[1][:charNumbers]]=f
  fd
end

Base.show(io::IO,::MIME"text/latex",uc::UnipotentCharacters)=print(io,TeXs(uc))

function Base.show(io::IO,uc::UnipotentCharacters)
  print(io,"UnipotentCharacters(",spets(uc),")")
end

function Base.show(io::IO,::MIME"text/plain",uc::UnipotentCharacters)
  TeX=get(io,:TeX,false)
  print(io,TeX ? "\$\\mbox{UnipotentCharacters}" : "UnipotentCharacters")
  print(io,"(",spets(uc),")")
  println(io,TeX ? "\$" : "")
  col_labels=["n_0"]
  m=hcat(repr.(1:length(uc)))
  row_labels=charnames(io,uc)
  almost=almostcharnames(io,uc)
  if almost!=row_labels
    m=hcat(m,almost)
    push!(col_labels,"almostch")
  end
  io=IOContext(io,:varname=>:q)
  cycpol=get(io,:cycpol,true)
  m=hcat(m,xrepr.(io,cycpol ? CycPoldegrees(uc) : degrees(uc)))
  push!(col_labels,"\\mbox{Deg}(\\gamma)")
  feg=fakedegrees(uc)
  m=hcat(m,xrepr.(io,cycpol ? CycPol.(feg) : feg))
  push!(col_labels,"\\mbox{Feg}")
  if haskey(uc,:charSymbols) && (uc.charSymbols!=uc.charParams)
    m=hcat(m,xrepr.(io,first.(uc.charSymbols)))
    push!(col_labels,"\\mbox{Symbol}")
  end
  m=hcat(m,xrepr.(io,eigen(uc)))
  push!(col_labels,"\\mbox{Fr}(\\gamma)")
  m=hcat(m,fromTeX.(io,labels(uc)))
  push!(col_labels,"\\mbox{label}")
  if get(io,:signs,false)
    m=hcat(m,xrepr.(io,signs(uc)))
    push!(col_labels,"\\mbox{signs}")
  end
  if get(io,:byfamily,false)
    rows=vcat(map(x->x[:charNumbers],uc.families)...)
    row_seps=vcat([-1,0],rows[cumsum(length.(uc.families))])
    for f in uc.families
      if special(f)==1 && cospecial(f)==special(f) continue end
      row_labels[f.charNumbers[special(f)]]*="^{s}"
      if cospecial(f)==special(f) continue end
      row_labels[f.charNumbers[cospecial(f)]]*="^{c}"
    end
  else
    rows=get(io,:rows,1:length(uc))
    row_seps=get(io,:row_seps,[-1,0,length(uc)])
  end
  showtable(io,m;row_labels,rows,rows_label="\\gamma",row_seps,col_labels)
end

Cosets.spets(uc::UnipotentCharacters)=uc.spets

Base.length(uc::UnipotentCharacters)=length(uc.charParams)

"`fourier(uc::UnipotentCharacters)` the Lusztig Fourier matrix for `uc`."
function Families.fourier(uc::UnipotentCharacters)
  get!(uc,:fourier)do
    l=length(uc)
    T=reduce(promote_type,eltype.(improve_type(fourier.(uc.families))))
    i=fill(T(0),l,l)
    for f in uc.families
     i[f.charNumbers,f.charNumbers]=fourier(f)
    end
    i
  end
end

function Families.qeigen(uc::UnipotentCharacters)
  get!(uc,:qeigen)do
    res=zeros(Rational{Int},length(uc))
    for f in uc.harishChandra
      if haskey(f,:qEigen)
        res[f[:charNumbers]]=fill(f[:qEigen],length(f[:charNumbers]))
      end
    end
    res
  end
end

"""
`degrees(uc::UnipotentCharacters,q=Pol())`

Returns  the  list  of  degrees  of  the unipotent characters of the finite
reductive group (or Spetses) with Weyl group (or Spetsial reflection group)
`W`, evaluated at `q`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> uc=UnipotentCharacters(W);

julia> degrees(uc)
10-element Vector{Pol{Rational{Int64}}}:
 1
 q⁶
 (1//3)q⁵+(1//3)q³+(1//3)q
 (1//3)q⁵+(1//3)q³+(1//3)q
 (1//6)q⁵+(1//2)q⁴+(2//3)q³+(1//2)q²+(1//6)q
 (1//2)q⁵+(1//2)q⁴+(1//2)q²+(1//2)q
 (1//2)q⁵+(-1//2)q⁴+(-1//2)q²+(1//2)q
 (1//6)q⁵+(-1//2)q⁴+(2//3)q³+(-1//2)q²+(1//6)q
 (1//3)q⁵+(-2//3)q³+(1//3)q
 (1//3)q⁵+(-2//3)q³+(1//3)q
```
"""
function Chevie.degrees(uc::UnipotentCharacters,q=Pol())
  if !haskey(uc,:degrees) uc.degrees=Dict{Any,Any}() end
  d=uc.degrees
  if haskey(d,q) return d[q] end
  v=uc.almostHarishChandra[1][:charNumbers]
  d[q]=improve_type(fourier(uc)[v,:]')*fakedegrees(spets(uc),q)
end

"""
`CycPoldegrees(uc::UnipotentCharacters)`

Taking  advantage that  the degrees  of unipotent  characters of the finite
reductive group (or Spetses) with Weyl group (or Spetsial reflection group)
`W`  are products  of cyclotomic  polynomials, this  function returns these
degrees as a list of `CycPol`s. It is faster than  `CycPol.(degrees(uc))`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> CycPoldegrees(UnipotentCharacters(W))
10-element Vector{CycPol{Rational{Int64}}}:
 1
 q⁶
 qΦ₃Φ₆/3
 qΦ₃Φ₆/3
 qΦ₂²Φ₃/6
 qΦ₂²Φ₆/2
 qΦ₁²Φ₃/2
 qΦ₁²Φ₆/6
 qΦ₁²Φ₂²/3
 qΦ₁²Φ₂²/3
```
"""
function CycPoldegrees(uc::UnipotentCharacters)
  get!(uc,:cycpoldegrees) do
    CycPol.(degrees(uc))
  end
end

"`eigen(uc::UnipotentCharacters)` the eigenvalues of Frobenius for `uc`"
function LinearAlgebra.eigen(uc::UnipotentCharacters)
  get!(uc,:eigen)do
    eig=fill(E(1),sum(length,uc.families))
    for f in uc.families eig[f.charNumbers]=eigen(f) end
    eig
  end::Vector{Root1}
end

"`signs(uc::UnipotentCharacters)` the signs in families of  `uc`"
function SignedPerms.signs(uc::UnipotentCharacters)
  get!(uc,:signs)do
    sgn=fill(1,sum(length,uc.families))
    for f in uc.families 
      if haskey(f,:signs) sgn[f.charNumbers]=signs(f) end
    end
    sgn
  end::Vector{Int}
end

function labels(uc::UnipotentCharacters)::Vector{String}
  get!(uc,:labels)do
    lab=fill("",length(uc))
    for f in uc.families lab[f.charNumbers]=f.charLabels
    end
    lab
  end
end

#-------------------------- UniChars -------------------------------
@GapObj struct UniChar{T,C}
  group::T
  v::Vector{C}
end

UniChar(W,r)=UniChar(W,r,Dict{Symbol,Any}())

"""
`unipotent_character(W,l)` or `unichar(W,l)`

Constructs  an object representing the unipotent character specified by `l`
of  the algebraic  group associated  to the  Coxeter group or Coxeter coset
specified  by `W`. There are 3 possibilities  for `l`: if it is an integer,
the  `l`-th unipotent character of `W` is  returned. If it is a string, the
unipotent  character of `W` whose name is  `l` is returned (where the names
are as given by `charnames(UnipotentCharacters(W))`). Finally, `l` can be a
list  of length the number of  unipotent characters of `W`, which specifies
the coefficient to give to each unipotent character.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> u=unichar(W,7)
[G₂]:<G₂[-1]>

julia> v=unichar(W,"G2[E3]")
[G₂]:<G₂[ζ₃]>

julia> w=unichar(W,[1,0,0,-1,0,0,2,0,0,1])
[G₂]:<φ₁‚₀>-<φ″₁‚₃>+2<G₂[-1]>+<G₂[ζ₃²]>

julia> unichar(W,fourier(UnipotentCharacters(W))[3,:])
[G₂]:2//3<φ′₁‚₃>-1//3<φ″₁‚₃>+1//3<φ₂‚₁>+1//3<G₂[1]>-1//3<G₂[ζ₃]>-1//3<G₂[ζ₃²]>
```
The  last line shows  the almost character  associated to the 3rd unipotent
character of `W`.

some limited arithmetic is available on unipotent characters:

```julia-repl
julia> coefficients(u) # so that u==unichar(W,coefficients(u))
10-element Vector{Int64}:
 0
 0
 0
 0
 0
 0
 1
 0
 0
 0

julia> w-2u
[G₂]:<φ₁‚₀>-<φ″₁‚₃>+<G₂[ζ₃²]>

julia> scalar_product(w,w)
7

julia> degree(w)
Pol{Int64}: q⁵-q⁴-q³-q²+q+1
```
"""
function unipotent_character(W,v::Int)
  r=zeros(Int,length(UnipotentCharacters(W)))
  r[v] = 1
  UniChar(W,r)
end

const unichar=unipotent_character

unichar(W,v::String)=
  unichar(W,findfirst(==(v),charnames(UnipotentCharacters(W))))

unichar(W,v::AbstractVector)=UniChar(W,collect(v))

LaurentPolynomials.coefficients(r::UniChar)=r.v

"""
`Base.show(io::IO,w::UniChar)`

The  formatting  of  unipotent  characters  is  affected  by  IO property
:compact .  If `true` (the default) they are printed in a compact form.
Otherwise, they are printed one unipotent character per line:

```julia-rep1
julia> xdisplay(w;compact=false)
[G₂]:
<φ₁‚₀>    1
<φ₁‚₆>    0
<φ′₁‚₃>   0
<φ″₁‚₃>   -1
<φ₂‚₁>    0
<φ₂‚₂>    0
<G₂[-1]>  2
<G₂[1]>   0
<G₂[ζ₃]>  0
<G₂[ζ₃²]> 1
```
"""
function Base.show(io::IO,r::UniChar)
  if !hasdecor(io)
    print(io,"unichar(",r.group,",",coefficients(r),")")
    return
  end
  if haskey(r,:name) print(io,"[",r.group,":",r.name,"]:")
  else print(io,"[",r.group,"]:")
  end
  res=""
  s=get(io,:names,String[])
  if isempty(s) s=charnames(io,UnipotentCharacters(r.group)) end
  m=maximum(length.(s))+3
  for (i,c) in pairs(r.v)
    n = "<"*s[i]*">"
    if get(io,:compact,true)
      if !iszero(c)
        if isone(c) res*= "+"
        elseif isone(-c) res*="-"
        else
          c=xrepr(io,c)
          if occursin(r".[+-]",c) c = "("* c* ")" end
          if !(c[1] in "+-") res*="+" end
          res*=c
        end
        res*=n
      end
     elseif !iszero(c) || !get(io,:nozero,false)
      res*="\n"*rpad(n,m)*xrepr(io,c)
    end
  end
  if isempty(res) res="0" end
  if res[1]=='+' res=res[2:end] end
  print(io,res)
end

Chars.scalar_product(u1::UniChar,u2::UniChar)=sum(u1.v .* conj.(u2.v))
Base.:+(u1::UniChar,u2::UniChar)=UniChar(u1.group,u1.v+u2.v)
Base.:-(u1::UniChar,u2::UniChar)=UniChar(u1.group,u1.v-u2.v)
Base.:*(u1::UniChar,a)=UniChar(u1.group,u1.v .* a)
Base.:*(a,u1::UniChar)=u1*a
Base.:(==)(u1::UniChar,u2::UniChar)=u1.group==u2.group && u1.v==u2.v

LaurentPolynomials.degree(u::UniChar,q=Pol(:q))=improve_type(sum(u.v .*
                                     degrees(UnipotentCharacters(u.group),q)))

"""
`lusztig_induce(W,u)`

`u`  should be a unipotent character of a parabolic subcoset of the Coxeter
coset  `W`. It represents  a unipotent character  `λ` of a  Levi `𝐋` of the
algebraic  group  `𝐆`  attached  to  `W`.  The  program returns the Lusztig
induced ``R_𝐋^𝐆(λ)``.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> WF=spets(W)
G₂

julia> T=subspets(WF,Int[],W(1))
G₂₍₎=Φ₁Φ₂

julia> u=unichar(T,1)
[G₂₍₎=Φ₁Φ₂]:<Id>

julia> lusztig_induce(WF,u)
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> dlchar(W,W(1))
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>
```
"""
function lusztig_induce(WF, u)
  t=lusztig_induction_table(u.group, WF)
  if !isnothing(t) unichar(WF, improve_type(t.scalar*u.v)) end
end

"""
`lusztig_restrict(R,u)`

`u`  should be a unipotent character of a parent Coxeter coset `W` of which
`R` is a parabolic subcoset. It represents a unipotent character `γ` of the
algebraic  group `𝐆` attached to `W`,  while `R` represents a Levi subgroup
`L`. The program returns the Lusztig restriction ``*R_𝐋^𝐆(γ)``.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> WF=spets(W)
G₂

julia> T=subspets(WF,Int[],W(1))
G₂₍₎=Φ₁Φ₂

julia> u=dlchar(W,W(1))
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> lusztig_restrict(T,u)
[G₂₍₎=Φ₁Φ₂]:4<Id>

julia> T=subspets(WF,Int[],W(2))
G₂₍₎=Φ₁Φ₂

julia> lusztig_restrict(T,u)
[G₂₍₎=Φ₁Φ₂]:0
```
"""
lusztig_restrict(HF,u)=unichar(HF,improve_type(transpose(
                            lusztig_induction_table(HF,u.group).scalar)*u.v))

harish_chandra_induce(WF,u)=
   unichar(WF,improve_type(hc_induction_table(u.group,WF).scalar*u.v))
const hc_induce=harish_chandra_induce

harish_chandra_restrict(HF,u)=unichar(HF,improve_type(u.v*hc_induction_table(HF,u.group).scalar))
const hc_restrict=harish_chandra_restrict

"""
`deligne_lusztigCharTable(W)` or `dlCharTable(W)`

for  each conjugacy class of `W`, gives the decomposition of `R_{T_w}^G` in
unipotent characters.

```julia-repl
julia> dlCharTable(W)
6×10 Matrix{Int64}:
 1   1   1   1   2   2   0   0   0   0
 1  -1   1  -1   0   0   0   0   0   0
 1  -1  -1   1   0   0   0   0   0   0
 1   1   0   0  -1   0   1   0   1   1
 1   1   0   0   0  -1   0   1  -1  -1
 1   1  -1  -1   0   0  -2  -2   0   0
```
"""
function deligne_lusztigCharTable(W)
  get!(W,:rwTable)do
    uc=UnipotentCharacters(W)
    improve_type(CharTable(W).irr'*fourier(uc)[charnumbers(uc.almostHarishChandra[1]),:])
  end
end

const dlCharTable=deligne_lusztigCharTable

"""
`deligne_lusztig_character(W,w)` or `dlchar(W,w)`

This  function returns the Deligne-Lusztig character  ``R_𝐓 ^𝐆 (1)`` of the
algebraic  group `𝐆 ` associated to the Coxeter group or Coxeter coset `W`.
The  torus  `𝐓`  can  be  specified  in  3  ways:  if `w` is an integer, it
represents the `w`-th conjugacy class (or `phi`-conjugacy class for a coset
`Wϕ`)  of `W`. Otherwise  `w` can be  a word or  an element of  `W`, and it
represents the class (or `ϕ`-class) of `w`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> dlchar(W,3)
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> dlchar(W,W(1))
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> dlchar(W,[1])
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> dlchar(W,[1,2])
[G₂]:<φ₁‚₀>+<φ₁‚₆>-<φ₂‚₁>+<G₂[-1]>+<G₂[ζ₃]>+<G₂[ζ₃²]>
```
"""
deligne_lusztig_character(W,i::Int)=unichar(W,dlCharTable(W)[i,:])
const dlchar=deligne_lusztig_character

dlchar(W,w::Perm)=dlchar(W,position_class(W,w))

dlchar(W,w::Vector{Int})=dlchar(W,W(w...))

"""
`almost_character(W,i)` or `almostchar(W,i)`

This  function  returns  the  `i`-th  almost  unipotent  character  of  the
algebraic  group 𝐆 associated to the Coxeter group or Coxeter coset `W`. If
`φ` is the `i`-th irreducible character of `W`, the `i`-th almost character
is  ``R_φ=W⁻¹∑_{w∈ W}  φ(w) R_{𝐓_w}^𝐆  (1)`` where  ``𝐓_w`` is  the maximal
torus  associated  to  the  conjugacy  class  (or `ϕ`-conjugacy class for a
coset) of `w`.

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> almostchar(W,3)
[B₂]:<.11>

julia> almostchar(W,1)
[B₂]:1//2<11.>+1//2<1.1>-1//2<.2>-1//2<B₂>
```
"""
almost_character=function(W,i)
  ct=CharTable(W)
  dl=dlchar.(Ref(W),1:length(ct.charnames))
  sum(ct.irr[i,:] .* classes(ct).//length(W).*dl)
end

const almostchar=almost_character

"""
`deligne_lusztig_lefschetz(h,m=0)` or `dllefschetz(h,m=0)`

Here `h` is an element of a Hecke algebra associated to a Coxeter group `W`
or  Coxeter coset `Wϕ` which itself is  associated to an algebraic group `𝐆
`.  By [dm85](@cite), for ``g∈ 𝐆^F``, the number of fixed points of `Fᵐ` on
the Deligne-Lusztig variety associated to the element `wϕ∈Wϕ`, have for `m`
divisible   by   a   sufficently   large   integer   `d`,  the  form  ``∑_φ
φ_{(qᵐ)}(T_wϕ)R_φ(g)``  where `φ`  runs over  the irreducible characters of
``Wϕ``,  where  ``R_φ``  is  the  corresponding almost character, and where
``φ_{(qᵐ)}``  is a  character value  of the  Hecke algebra  ``H(Wϕ,qᵐ)`` of
``Wϕ``  with  parameter  `qᵐ`.  This  expression  is  called the *Lefschetz
character*  of  the  Deligne-Lusztig  variety.  If  we  consider `qᵐ` as an
indeterminate  `x`, it can  be seen as  a sum of  unipotent characters with
coefficients  character values of the  generic Hecke algebra ``H(Wϕ,x)``. A
more complicated formula involving the eigenvalues of Frobenius attached to
unipotent characters applies for `m` not prime to `d`. The function returns
this formula when a second parameter `m≠0` is given.

The  function 'dllefschetz' takes  as argument a  Hecke element and returns
the  corresponding Lefschetz character. This is defined on the whole of the
Hecke  algebra by linearity.  The Lefschetz character  of various varieties
related   to   Deligne-Lusztig   varieties,   like   their  completions  or
desingularisation,  can be  obtained by  taking the  Lefschetz character at
various elements of the Hecke algebra.

```julia-repl
julia> W=coxgroup(:A,2)
A₂

julia> H=hecke(W,Pol(:q))
hecke(A₂,q)

julia> T=Tbasis(H);

julia> dllefschetz(T(1,2))
[A₂]:<111>-q<21>+q²<3>

julia> dllefschetz((T(1)+T())*(T(2)+T()))
[A₂]:q<21>+(q²+2q+1)<3>
```

The   last  line  shows  the   Lefschetz  character  of  the  Samelson-Bott
desingularisation of the Coxeter element Deligne-Lusztig variety.

We now show an example with a coset (corresponding to the unitary group).

```julia-repl
julia> H=hecke(spets(W,Perm(1,2)),Pol(:q)^2)
hecke(²A₂,q²)

julia> T=Tbasis(H);dllefschetz(T(1))
[²A₂]:-<11>-q<²A₂>+q²<2>
```
Finally,  there is a second form `dllefschetz(H::HeckeAlgebra,w,i=0)` where
the  arguments are a Hecke algebra and an  element of `w`. This may be used
for  Spetses where we know the column of the `CharTable` of `H` for `w` but
not other columns of the spetsial Hecke algebra charcater table.
"""
function deligne_lusztig_lefschetz(h,i=0)
  W=h.H.W
  uc=UnipotentCharacters(W)
  uniform=charnumbers(uc.almostHarishChandra[1])
  unichar(W,improve_type((char_values(h)'*fourier(uc)[uniform,:])[1,:].*eigen(uc).^i))
end

const dllefschetz=deligne_lusztig_lefschetz

function dllefschetz(H::HeckeAlgebra,w,i=0)
  W=H.W
  uc=UnipotentCharacters(W)
  uniform=charnumbers(uc.almostHarishChandra[1])
  unichar(W,improve_type((char_values(H,w)'*fourier(uc)[uniform,:])[1,:].*eigen(uc).^i))
end

function dllefschetzTable(H,i=0)
  WF=H.W
  t=CharTable(H).irr
  uc=UnipotentCharacters(WF)
  improve_type(t'*fourier(uc)[charnumbers(uc.almostHarishChandra[1]),:]*
               Diagonal(Cyc.(eigen(uc)))^i)
end

"""
`on_unipotents(W,aut)`

`W`  is  a  reflection  group  or  reflection  coset  representing a finite
reductive group ``𝐆 ^F``, and `aut` is an automorphism of ``𝐆 ^F`` (for `W`
a  permutation group, this can be given as a permutation of the roots). The
function  returns the permutation  of the unipotent  characters of ``𝐆 ^F``
induced  by `aut`. This makes sense  for Spetsial complex reflection groups
and is implemented for them.

```julia-repl
julia> WF=rootdatum("3D4")
³D₄

julia> on_unipotents(Group(WF),WF.phi)
(1,7,2)(8,12,9)
```
"""
function on_unipotents(W,aut)
  uc=UnipotentCharacters(W)
  t=dlCharTable(W)
  t=vcat(t,transpose(eigen(uc)))
  l=fill(0,length(uc))
  n=charnumbers(uc.harishChandra[1])
  l[n]=1:length(n)
  t=vcat(t,transpose(l))
  if length(unique(eachcol(t)))<size(t,2)
    error("Rw + eigen + principal series cannot disambiguate\n")
  end
  t1=invpermute(t[1:end-1,:],on_classes(W, aut),dims=1)
  l[n]=l[n].^inv(on_chars(W,aut))
  t1=vcat(t1,transpose(l))
  Perm(t,t1,dims=2)
end

function Cosets.Frobenius(x::UniChar, phi)
  W=x.group
  unichar(W,invpermute(x.v,inv(on_unipotents(W,phi))))
end

cuspidal(uc::UnipotentCharacters,d::Integer)=cuspidal(uc,E(d))
cuspidal(uc::UnipotentCharacters,d::Rational)=cuspidal(uc,Root1(;r=d))
"""
`cuspidal(uc::UnipotentCharacters[,e])`

A  unipotent character `γ` of a  finite reductive group `𝐆` is `e`-cuspidal
if  its  Lusztig  restriction  to  any  proper `e`-split Levi is zero. When
`e==1`  (the default when  `e` is omitted)  we recover the  usual notion of
cuspidal character. Equivalently the `Φₑ`-part of the generic degree of `γ`
is equal to the `Φₑ`-part of the generic order of the adjoint group of `𝐆`.
This  makes  sense  for  any  Spetsial  complex  reflection  group  and  is
implemented for them.

The  function returns the list of indices of unipotent characters which are
`e`-cuspidal.

```julia-repl
julia> W=coxgroup(:D,4)
D₄

julia> cuspidal(UnipotentCharacters(W))
1-element Vector{Int64}:
 14

julia> cuspidal(UnipotentCharacters(W),6)
8-element Vector{Int64}:
  1
  2
  6
  7
  8
  9
 10
 12

julia> cuspidal(UnipotentCharacters(complex_reflection_group(4)),3)
4-element Vector{Int64}:
  3
  6
  7
 10
```
"""
function cuspidal(uc::UnipotentCharacters,d=E(1))
  if length(uc)==1 return [1] end
  ud=CycPoldegrees(uc)
  ad=count(!isone,relative_degrees(spets(uc),d))
  filter(i->ad==valuation(ud[i],d),eachindex(ud))
end

"""
`cuspidal_data(W[,d[,ad]];proper=false,all=false)`

returns  named tuples `(levi=LF,cuspidal=λ,d=d)` where  `LF` is a `d`-split
Levi  (with `d`-center  of dimension  `ad` if  `ad` is  given) and `λ` is a
`d`-cuspidal  character of  `LF`. If  `d=1` this  returns ordinary cuspidal
characters.  The  character  `λ`  is  given  as  its  index  in the list of
unipotent  characters. If `d` was given as  an integer, it is returned as a
`Root1` representing `E(d)`.

For a  given `d`, the pairs `(levi=LF,cuspidal=λ)` are called `d`-cuspidal
pairs.

If  the keyword `proper=true`  is given, only  the tuples where `LF!=W` (or
equivalently `ad>0`) are returned.

If  `d` is omitted, data  for all `d` orders  of eigenvalues of elements of
`W`  is returned. If in addition  the keyword argument `all=true` is given,
data  for all eigenvalues of elements of `W`  is returned --- a pair may be
cuspidal for some `E(d)^i` while not cuspidal for `E(d)`.

```julia-repl
julia> cuspidal_data(coxgroup(:F,4),1)
9-element Vector{@NamedTuple{levi::Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}, cuspidal::Int64, d::Root1}}:
 (levi = F₄, cuspidal = 31, d = 1)
 (levi = F₄, cuspidal = 32, d = 1)
 (levi = F₄, cuspidal = 33, d = 1)
 (levi = F₄, cuspidal = 34, d = 1)
 (levi = F₄, cuspidal = 35, d = 1)
 (levi = F₄, cuspidal = 36, d = 1)
 (levi = F₄, cuspidal = 37, d = 1)
 (levi = F₄₍₂₃₎=B₂₍₂₁₎Φ₁², cuspidal = 6, d = 1)
 (levi = F₄₍₎=Φ₁⁴, cuspidal = 1, d = 1)

julia> cuspidal_data(complex_reflection_group(4),3)
5-element Vector{@NamedTuple{levi::Spets{PRSG{Cyc{Rational{Int64}}, Int16}}, cuspidal::Int64, d::Root1}}:
 (levi = G₄, cuspidal = 3, d = ζ₃)
 (levi = G₄, cuspidal = 6, d = ζ₃)
 (levi = G₄, cuspidal = 7, d = ζ₃)
 (levi = G₄, cuspidal = 10, d = ζ₃)
 (levi = G₄₍₎=Φ₁Φ′₃, cuspidal = 1, d = ζ₃)

julia> cuspidal_data(crg(4),E(3,2))
5-element Vector{@NamedTuple{levi::Spets{PRSG{Cyc{Rational{Int64}}, Int16}}, cuspidal::Int64, d::Root1}}:
 (levi = G₄, cuspidal = 2, d = ζ₃²)
 (levi = G₄, cuspidal = 5, d = ζ₃²)
 (levi = G₄, cuspidal = 7, d = ζ₃²)
 (levi = G₄, cuspidal = 10, d = ζ₃²)
 (levi = G₄₍₎=Φ₁Φ″₃, cuspidal = 1, d = ζ₃²)
```
"""
cuspidal_data(W,d::Integer,ad)=cuspidal_data(W,E(d),ad)
cuspidal_data(W,d::Rational,ad)=cuspidal_data(W,Root1(;r=d),ad)
cuspidal_data(W,d::Root1,ad)=[(levi=L,cuspidal=char,d=d)
                        for L in split_levis(W, d, ad)
                        for char in cuspidal(UnipotentCharacters(L),d)]

cuspidal_data(W,d;proper=false)=[p for ad in
         (proper ? 1 : 0):length(relative_degrees(W,d))
         for p in cuspidal_data(W,d,ad)]

cuspidal_data(W;proper=false,all=false)=[p for d in
  sort(unique(all ? vcat(refleigen(W)...) : order.(vcat(refleigen(W)...))))
  for p in cuspidal_data(W,d;proper)]

function relative_hecke(uc::UnipotentCharacters,i,q)
  hw=uc.harishChandra[i]
  return hecke(reflection_group(hw[:relativeType]),
    map(hw[:parameterExponents])do i
    if i isa Vector return map(j->E(length(i),j-1)*q^i[j],1:length(i))
#   elif i<0  then return -q^-i; # JM 14/2/2018 I think obsolete
    else return q^i
    end end)
end

end
