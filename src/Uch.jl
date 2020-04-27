"""
Let  `𝐆 ` be a connected reductive group defined over the algebraic closure
of  a finite field `𝔽_q`, with corresponding Frobenius automorphism `F`, or
more  generally  let  `F`  be  an  isogeny  of  `𝐆 ` such that a power is a
Frobenius (this covers the Suzuki and Ree groups).

If  `𝐓`  is  an  `F`-stable  maximal  torus  of  `𝐆  `,  and  `𝐁` is a (not
necessarily  `F`-stable)  Borel  subgroup  containing  `𝐓`,  we  define the
*Deligne-Lusztig*  variety `X_𝐁={g𝐁 ∈ 𝐆 /𝐁 ∣ g𝐁 ∩ F(g𝐁 )≠∅ }`. This variety
has  a  natural  action  of  `𝐆  ^F`  on  the  left,  so  the corresponding
*Deligne-Lusztig  virtual  module*  `∑ᵢ  (-1)ⁱ  Hⁱ_c(X_𝐁,ℚ̄  _ℓ)` also. The
character of this virtual module is the *Deligne-Lusztig* character `R_𝐓 ^𝐆
(1)`; the notation reflects the fact that one can prove that this character
does  not  depend  on  the  choice  of  `𝐁`.  Actually,  this  character is
parameterized by an `F`-conjugacy class of `W`: if `𝐓₀⊂𝐁₀` is an `F`-stable
pair,  there is an unique `w∈ W=N_𝐆 (𝐓₀)/𝐓₀` such that the triple `(𝐓,𝐁,F)`
is  `𝐆 `-conjugate to `(𝐓₀,𝐁₀,wF)`. In this case we denote `R_w` for `R_𝐓^𝐆
(1)`; it depends only on the `F`-class of `w`.

The  *unipotent characters* of  `𝐆 ^F` are  the irreducible constituents of
the `R_w`. In a similar way that the unipotent classes are a building block
for  describing the conjugacy  classes of a  reductive group, the unipotent
characters  are  a  building  block  for  the  irreducible  characters of a
reductive  group.  They  can  be  parameterized  by combinatorial data that
Lusztig  has attached just to the coset `Wφ`, where `φ` is the finite order
automorphism  of  `X(𝐓₀)`  such  that  `F=qφ`.  Thus, from the viewpoint of
Chevie, they are objects combinatorially attached to a Coxeter coset.

A  subset  of  the  unipotent  characters, the *principal series* unipotent
characters,   can  be  described  in  an   elementary  way.  They  are  the
constituents  of `R₁`, or equivalently the characters of the virtual module
defined  by the cohomology of `X_{𝐁 ₀}`,  which is the discrete variety `(𝐆
/𝐁₀)^F`;  the virtual  module reduces  to the  actual module `ℚ̄ _ℓ[(𝐆 /𝐁₀)
^F]`.   Thus  the   Deligne-Lusztig  induction   `R_𝐓₀^𝐆  (1)`  reduces  to
Harish-Chandra  induction,  defined  as  follows:  let  `𝐏  =𝐔  ⋊ 𝐋 ` be an
`F`-stable  Levi decomposition of an `F`-stable parabolic subgroup of `𝐆 `.
Then  the *Harish-Chandra* induced `R_𝐋^𝐆 ` of  a character `χ` of `𝐋^F` is
the  character `Ind_{𝐏^F}^{𝐆 ^F}χ̃`, where `χ̃` is the lift to `𝐏^F` of `χ`
via  the quotient `𝐏^F/𝐔 ^F=𝐋^F`;  Harish-Chandra induction is a particular
case  of *Lusztig induction*,  which is defined  when `𝐏` is not `F`-stable
using  the variety `X_𝐔  ={ g𝐔 ∈𝐆  /𝐔 ∣ g𝐔  ∩ F(g𝐔 )≠∅}`,  and gives for an
`𝐋^F`-module  a  virtual  `𝐆  ^F`-module.  Like  ordinary  induction, these
functors  have adjoint  functors going  from representations  of `𝐆  ^F` to
representations   (resp.   virtual   representations)   of   `𝐋^F`   called
Harish-Chandra restriction (resp. Lusztig restriction).

The  commuting  algebra  of  `𝐆^F`-endomorphisms  of  `R_{𝐓₀}^𝐆(1)`  is  an
Iwahori-Hecke  algebra for `W^φ`, with parameters  which are some powers of
`q`;  they  are  all  equal  to  `q`  when  `W^φ=W`.  Thus principal series
unipotent characters correspond to characters of `W^φ`.

To  understand the  decomposition of  Deligne-Lusztig characters,  and thus
unipotent  characters,  is  is  useful  to  introduce  another set of class
functions  which are parameterized  by irreducible characters  of the coset
`Wφ`.  If  `χ`  is  such  a  character,  we  define  the associated *almost
character* by: `Rᵪ=|W|⁻¹∑_{w∈ W}χ(wφ) R_w`. The reason to the name is that
these  class  function  are  close  to irreducible characters: they satisfy
`⟨Rᵪ, R_ψ⟩_{𝐆^F}=δ_{χ,ψ}`;  for  the  linear  and  unitary group they are
actually  unipotent characters (up to sign in the latter case). They are in
general  sum (with  rational coefficients)  of a  small number of unipotent
characters  in  the  same  *Lusztig  family*  (see  "Families  of unipotent
characters").  The degree of `Rᵪ` is a polynomial in `q` equal to the fake
degree  of  the  character  `χ`  of  `Wφ`  (see  "Functions  for Reflection
cosets").

We  now describe the parameterization of unipotent characters when `W^φ=W`,
thus  when the coset `Wφ` identifies with `W` (the situation is similar but
a  bit more difficult to describe  in general). The (rectangular) matrix of
scalar  products  `⟨ρ, Rᵪ⟩_{𝐆 ^F}`,  when  characters of `W` and unipotent
characters  are arranged in the right  order, is block-diagonal with rather
small blocks which are called *Lusztig families*.

For  the characters of `W` a family `𝓕` corresponds to a block of the Hecke
algebra  over a ring called the Rouquier  ring. To `𝓕` Lusztig associates a
small  group `Γ` (not bigger  than `(ℤ/2)^n`, or `𝔖ᵢ`  for `i≤5`) such that
the  unipotent  characters  in  the  family  are parameterized by the pairs
`(x,θ)`  taken up to  `Γ`-conjugacy, where `x∈Γ`  and `θ` is an irreducible
character  of  `C_Γ(x)`.  Further,  the  elements  of  `𝓕`  themselves  are
parameterized  by a  subset of  such pairs,  and Lusztig  defines a pairing
between  such pairs which computes the scalar product `⟨ρ, Rᵪ⟩_{𝐆^F}`. For
more details see "DrinfeldDouble".

A  second parameterization  of unipotent  character is  via *Harish-Chandra
series*.  A character is called *cuspidal* if all its proper Harish-Chandra
restrictions  vanish. There are few  cuspidal unipotent characters (none in
linear   groups,  and  at   most  one  in   other  classical  groups).  The
`𝐆^F`-endomorphism  algebra of an  Harish-Chandra induced `R_{𝐋^F}^{𝐆^F}λ`,
where `λ` is a cuspidal unipotent character turns out to be a Hecke algebra
associated to the group `W_{𝐆^F}(𝐋^F):=N_{𝐆^F}(𝐋)/𝐋`, which turns out to be
a  Coxeter group.  Thus another  parameterization is  by triples `(𝐋,λ,φ)`,
where  `λ`  is  a  cuspidal  unipotent  character  of  `𝐋^F`  and `φ` is an
irreducible   character  of  the   *relative  group*  `W_{𝐆^F}(𝐋^F)`.  Such
characters  are said to  belong to the  Harish-Chandra series determined by
`(𝐋,λ)`.

A  final  piece  of  information  attached  to  unipotent characters is the
*eigenvalues  of Frobenius*. Let `F^δ` be the smallest power of the isogeny
`F` which is a split Frobenius (that is, `F^δ` is a Frobenius and `φ^δ=1`).
Then  `F^δ` acts  naturally on  Deligne-Lusztig varieties  and thus  on the
corresponding  virtual modules, and  commutes to the  action of `𝐆^F`; thus
for  a given  unipotent character  `ρ`, a  submodule of  the virtual module
which  affords `ρ`  affords a  single eigenvalue  `μ` of  `F^δ`. Results of
Lusztig  and  Digne-Michel  show  that  this  eigenvalue  is  of  the  form
`q^{aδ}λ_ρ` where `2a∈ℤ` and `λ_ρ` is a root of unity which depends only on
`ρ`  and not the considered module. This  `λ_ρ` is called the eigenvalue of
Frobenius  attached  to  `ρ`.  Unipotent  characters  in the Harish-Chandra
series of a pair `(𝐋,λ)` have the same eigenvalue of Frobenius as `λ`.

Chevie   contains  tables  of  all   this  information,   and  can  compute
Harish-Chandra  and Lusztig  induction of  unipotent characters  and almost
characters. We illustrate the information on some examples:

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> uc=UnipotentCharacters(W)
UnipotentCharacters(G₂)
      γ│   Deg(γ)  Feg Fr(γ)    label
───────┼──────────────────────────────
φ₁‚₀   │        1    1     1         
φ₁‚₆   │       q⁶   q⁶     1         
φ′₁‚₃  │  qΦ₃Φ₆/3   q³     1    (1,ρ)
φ″₁‚₃  │  qΦ₃Φ₆/3   q³     1   (g₃,1)
φ₂‚₁   │ qΦ₂²Φ₃/6  qΦ₈     1    (1,1)
φ₂‚₂   │ qΦ₂²Φ₆/2 q²Φ₄     1   (g₂,1)
G₂[-1] │ qΦ₁²Φ₃/2    0    -1   (g₂,ε)
G₂[1]  │ qΦ₁²Φ₆/6    0     1    (1,ε)
G₂[ζ₃] │qΦ₁²Φ₂²/3    0    ζ₃  (g₃,ζ₃)
G₂[ζ₃²]│qΦ₁²Φ₂²/3    0   ζ₃² (g₃,ζ₃²)
```

The first column gives the name of the unipotent character; the first 6 are
in  the  principal  series  so  are  named  according  to the corresponding
characters  of `W`. The last 4 are cuspidal, and named by the corresponding
eigenvalue  of  Frobenius,  which  is  displayed  in  the fourth column. In
general   the   names   of   the   unipotent  characters  come  from  their
parameterization  by  Harish-Chandra  series;  in  addition,  for classical
groups, they are associated to *symbols*.

The first two characters are each in a family by themselves. The last eight
are  in a family associated to the group `Γ=𝔖_3`: the last column shows the
parameters  `(x,θ)`. The  second column  shows the  degree of the unipotent
characters, which is transformed by the Lusztig Fourier matrix of the third
column,  which gives the  degree of the  corresponding almost character, or
equivalently the fake degree of the corresponding character of `W`.

One  can get  more information  on the  Lusztig Fourier  matrix of  the big
family by asking

```julia-repl
julia> uc.families[1]
Family(D(S₃):[5, 6, 4, 3, 8, 7, 9, 10])
   label│eigen                                               
────────┼─────────────────────────────────────────────────────
(1,1)   │    1 1//6  1//2  1//3  1//3  1//6  1//2  1//3  1//3
(g₂,1)  │    1 1//2  1//2  0//1  0//1 -1//2 -1//2  0//1  0//1
(g₃,1)  │    1 1//3  0//1  2//3 -1//3  1//3  0//1 -1//3 -1//3
(1,ρ)   │    1 1//3  0//1 -1//3  2//3  1//3  0//1 -1//3 -1//3
(1,ε)   │    1 1//6 -1//2  1//3  1//3  1//6 -1//2  1//3  1//3
(g₂,ε)  │   -1 1//2 -1//2  0//1  0//1 -1//2  1//2  0//1  0//1
(g₃,ζ₃) │   ζ₃ 1//3  0//1 -1//3 -1//3  1//3  0//1  2//3 -1//3
(g₃,ζ₃²)│  ζ₃² 1//3  0//1 -1//3 -1//3  1//3  0//1 -1//3  2//3
```

One  can  do  computations  with  individual  unipotent characters. Here we
construct  the Coxeter torus, and then the identity character of this torus
as a unipotent character.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> T=spets(reflection_subgroup(W,Int[]),W(1,2))
.Φ₆

julia> u=UniChar(T,1)
[.Φ₆]:<.>
```

Then  here  are  two  ways  to  construct  the  Deligne-Lusztig  character
associated to the Coxeter torus:

```julia-repl
julia> LusztigInduce(W,u)
[G₂]:<φ₁‚₀>+<φ₁‚₆>-<φ₂‚₁>+<G₂[-1]>+<G₂[ζ₃]>+<G₂[ζ₃²]>

julia> v=DLChar(W,[1,2])
[G₂]:<φ₁‚₀>+<φ₁‚₆>-<φ₂‚₁>+<G₂[-1]>+<G₂[ζ₃]>+<G₂[ζ₃²]>

julia> degree(v)
Pol{Cyc{Rational{Int64}}}: q⁶+q⁵-q⁴-2q³-q²+q+1

julia> v*v
Cyc{Rational{Int64}}: 6
```

The  last two lines ask for the degree  of `v`, then for the scalar product
of `v` with itself.

Finally  we mention  that Chevie  can also  provide unipotent characters of
Spetses, as defined in [@BMM14]. An example:

```julia-repl
julia> UnipotentCharacters(ComplexReflectionGroup(4))
UnipotentCharacters(G₄)
    γ│               Deg(γ)    Feg Fr(γ)   label
─────┼───────────────────────────────────────────
φ₁‚₀ │                    1      1     1        
φ₁‚₄ │   (-√-3)q⁴Φ″₃Φ₄Φ″₆/6     q⁴     1  1∧-ζ₃²
φ₁‚₈ │    (√-3)q⁴Φ′₃Φ₄Φ′₆/6     q⁸     1  -1∧ζ₃²
φ₂‚₅ │            q⁴Φ₂²Φ₆/2   q⁵Φ₄     1   1∧ζ₃²
φ₂‚₃ │(-ζ₃-2ζ₃²)qΦ″₃Φ₄Φ′₆/3   q³Φ₄     1   1∧ζ₃²
φ₂‚₁ │(-2ζ₃-ζ₃²)qΦ′₃Φ₄Φ″₆/3    qΦ₄     1    1∧ζ₃
φ₃‚₂ │               q²Φ₃Φ₆ q²Φ₃Φ₆     1        
Z₃:2 │      (-√-3)qΦ₁Φ₂Φ₄/3      0   ζ₃²  ζ₃∧ζ₃²
Z₃:11│     (-√-3)q⁴Φ₁Φ₂Φ₄/3      0   ζ₃²  ζ₃∧-ζ₃
G₄   │           -q⁴Φ₁²Φ₃/2      0    -1 -ζ₃²∧-1
```
"""
module Uch

using Gapjm

export UnipotentCharacters, FixRelativeType, fourierinverse, UniChar,
AlmostChar, DLChar, DLLefschetz, LusztigInduce, LusztigRestrict

struct UnipotentCharacters
  harishChandra::Vector{Dict{Symbol,Any}}
  almostHarishChandra::Vector{Dict{Symbol,Any}}
  families::Vector{Family}
  prop::Dict{Symbol,Any}
end

function params_and_names(sers)
  function maketype(s)
    if s isa TypeIrred return s end
    if haskey(s,:orbit) 
      s[:orbit]=maketype.(s[:orbit])
    else s[:series]=Symbol(s[:series])
#     if s[:rank]==0 return Dict(:charnames=>[""],:charparams=>[[]]) end
    end
    TypeIrred(convert(Dict{Symbol,Any},s))
  end
  for ser in sers ser[:relativeType]=maketype(ser[:relativeType]) end
  chh=map(ser->charinfo(ser[:relativeType]),sers)
  l=sum(x->length(x[:charnames]),chh)
  res=Dict{Symbol,Any}()
  res[:charParams]=fill([],l)
  res[:TeXCharNames]=fill("",l)
  for (i,ser) in enumerate(sers)
    t=ser[:relativeType]
    n=ser[:cuspidalName]
    ch=chh[i]
    res[:charParams][ser[:charNumbers]]=map(x->[n,x],ch[:charparams])
    res[:TeXCharNames][ser[:charNumbers]]=map(ch[:charnames])do x
      rk=haskey(t,:orbit) ? t.orbit[1].rank : t.rank
      s=n
      if length(s)>0 && rk>0 s*=":" end
      if rk>0 s*=x end
      s
      end
  end
  res
end

function UnipotentCharacters(t::TypeIrred) 
  uc=copy(getchev(t,:UnipotentCharacters))
  if uc==false 
    println("Warning: $t is not a Spets!!")
    return false 
  end
  merge!(uc,params_and_names(uc[:harishChandra]))
  if !haskey(uc,:charSymbols) uc[:charSymbols]=uc[:charParams] end
  # adjust things for descent of scalars
  # we would like to adjust indices so they fit with those stored in t
  # but we cannot when indices mention non-generating reflections!
  a=length(t.orbit)
  if a>1
    if haskey(uc,:a) uc[:a].*=a end
    if haskey(uc,:A) uc[:A].*=a end
    for s in uc[:harishChandra]
      s[:parameterExponents].*=a
      s[:eigenvalue]^=a
      s[:cuspidalName]=join(map(i->s[:cuspidalName],1:a),"\\otimes ")
    end
  end

  if !haskey(uc,:almostHarishChandra)
    uc[:almostHarishChandra]=map(uc[:harishChandra])do s
    res=Dict{Symbol,Any}()
    for f in [:levi, :cuspidalName, :eigenvalue, :charNumbers] res[f]=s[f] end
    res[:relativeType]=TypeIrred(Dict(:orbit=>[copy(s[:relativeType])],:twist=>Perm()))
    if !isone(t.twist)
      a=t.orbit[1].indices[s.relativeType[:indices]]
      res[:relativeType][:twist]=prod(map(Perm,a,a.^t.twist))
    end
    res
    end
  else
    for s in uc[:almostHarishChandra]
      if !haskey(s[:relativeType],:orbit)
        s[:relativeType]=Dict(:orbit=>[s[:relativeType]],:twist=>Perm())
      end
    end
  end
  if !haskey(uc,:almostCharSymbols) uc[:almostCharSymbols]=uc[:charSymbols] end
  a=params_and_names(uc[:almostHarishChandra])
  uc[:almostCharParams]=a[:charParams]
  uc[:almostTeXCharNames]=a[:TeXCharNames]
  uc[:group]=t
  uch=UnipotentCharacters(uc[:harishChandra],uc[:almostHarishChandra],
                          Family.(copy(uc[:families])),uc)
  delete!(uc,:families)
  delete!(uc,:harishChandra)
  delete!(uc,:almostHarishChandra)
  uch
end

UnipotentCharacters(W::Group)=UnipotentCharacters(spets(W))

"""
`UnipotentCharacters(W)`

`W`  should be a Coxeter group, a  Coxeter Coset or a Spetses. The function
gives  back a record containing  information about the unipotent characters
of the associated algebraic group (or Spetses). This contains the following
fields:

`:group`: a pointer to `W`

`:charNames`:  the list of names of the unipotent characters.

`:charSymbols`: the list of symbols associated to unipotent characters,
for classical groups.

`:harishChandra`:  information  about  Harish-Chandra  series  of  unipotent
characters.  This is itself a list of records, one for each pair `(𝐋,λ)` of
a  Levi  of  an  `F`-stable  parabolic  subgroup  and  a cuspidal unipotent
character of `𝐋^F`. These records themselves have the following fields:

`:levi`: a list 'l' such that `𝐋` corresponds to 'ReflectionSubgroup(W,l)'.

`:cuspidalName`: the name of the unipotent cuspidal character `lambda`.

`:eigenvalue`: the eigenvalue of Frobenius for `λ`.

`:relativeType`: the reflection type of `W_𝐆 (𝐋)`;

`:parameterExponents`:  the  `𝐆 ^F`-endomorphism  algebra  of `R_𝐋^𝐆 (λ)` is a
Hecke algebra for `W_𝐆 (𝐋)` with some parameters of the form `q^{a_s}`. This
holds the list of exponents `a_s`.

`:charNumbers`:  the  indices  of  the  unipotent  characters indexed by the
irreducible characters of `W_𝐆 (𝐋)`.

`:families`:  information  about  Lusztig  families of unipotent characters.
This  is itself a list  of records, one for  each family. These records are
described in the section about families below.

```julia-repl
julia> W=coxgroup(:Bsym,2)
Bsym₂

julia> WF=spets(W,Perm(1,2))
²Bsym₂

julia> uc=UnipotentCharacters(WF)
UnipotentCharacters(²Bsym₂)
       γ│   Deg(γ)   Feg Fr(γ) label
────────┼────────────────────────────
2       │        1     1     1      
11      │       q⁴    q⁴     1      
²B₂[1,3]│√2qΦ₁Φ₂/2 qΦ₁Φ₂   ζ₈³     1
²B₂[1,5]│√2qΦ₁Φ₂/2     0   ζ₈⁵     2

julia> uc.families
3-element Array{Family,1}:
 Family(C₁:[1]) 
 Family(C₁:[2]) 
 Family(?4:3:4)

julia> uc.families[3]
Family(?4:3:4)
label│eigen    1     2
─────┼─────────────────
1    │  ζ₈³ √2/2 -√2/2
2    │  -ζ₈ √2/2  √2/2
```

`:charnames`:  returns  the  names  of  the  unipotent characters. Using the
version  with an additional  option record as  the second argument, one can
control the display in various ways.

```julia-repl
julia> uc=UnipotentCharacters(coxgroup(:G,2));

julia> charnames(uc;limit=true)
10-element Array{String,1}:
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
10-element Array{String,1}:
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

`:Display`:  One can control the display  of unipotent characters in various
ways.  In the record controlling 'Display', a field 'items' specifies which
columns are displayed. The possible values are

`:n0`:  The index of the character in the list of unipotent characters.

`:Name`:   The name of the unipotent character.

`:Degree`:  The degree of the unipotent character.

`:FakeDegree`: The degree of the corresponding almost character.

`:Eigenvalue`:  The eigenvalue of Frobenius attached to the unipotent
character.

`:Symbol`: for classical groups, the symbol attached to the unipotent
character.

`:Family`: The parameter the character has in its Lusztig family.

`:Signs`: The signs attached to the character in the Fourier transform.

The default value is
 'items:=[:Name,:Degree,:FakeDegree,:Eigenvalue,:Family]`

This  can be changed by setting the variable 'UnipotentCharactersOps.items`
which holds this default value. In addition if the field 'byFamily' is set,
the  characters are displayed  family by family  instead of in index order.
Finally,  the field 'chars' can be  set, indicating which characters are to
be displayed in which order.

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> uc=UnipotentCharacters(W)
UnipotentCharacters(B₂)
  γ│Deg(γ) Feg Fr(γ) label
───┼───────────────────────
11.│ qΦ₄/2  q²     1   +,-
1.1│qΦ₂²/2 qΦ₄     1   +,+
.11│    q⁴  q⁴     1      
2. │     1   1     1      
.2 │ qΦ₄/2  q²     1   -,+
B₂ │qΦ₁²/2   0    -1   -,-
```

    gap> Display(uc,rec(byFamily:=true));
    Unipotent characters for B2
    Name |  Degree FakeDegree Eigenvalue Label
    ___________________________________________
    *.11 |     q^4        q^4          1
    ___________________________________________
    11.  |  1/2qP4        q^2          1   +,-
    *1.1 |1/2qP2^2        qP4          1   +,+
    .2   |  1/2qP4        q^2          1   -,+
    B2   |1/2qP1^2          0         -1   -,-
    ___________________________________________
    *2.  |'|'|       1          1          1
    gap> Display(uc,items=[:n0,:Name,:Symbol]));
    Unipotent characters for B2
    n0 |Name   Symbol
    __________________
    1  | 11.   (12,0)
    2  | 1.1   (02,1)
    3  | .11 (012,12)
    4  |  2.     (2,)
    5  |  .2   (01,2)
    6  |  B2   (012,)|
"""
function UnipotentCharacters(WF::Spets) 
  gets(WF,:UnipotentCharacters) do
  function CartesianSeries(sers)
    ser=Dict{Symbol,Any}()
    ser[:levi]=reduce(vcat,getindex.(sers,:levi))
    ser[:relativeType]=filter(x->rank(x)!=0,getindex.(sers,:relativeType))
    if haskey(sers[1],:eigenvalue)
      ser[:eigenvalue]=prod(getindex.(sers,:eigenvalue))
    end
    if any(x->haskey(x,:qEigen),sers)
      ser[:qEigen]=sum(sers)do x
       if !haskey(x,:qEigen) return 0
       elseif x[:qEigen]==false return false
       else return x[:qEigen]
       end end
    else 
      ser[:qEigen]=0
    end
    if all(haskey.(sers,:parameterExponents))
      ser[:parameterExponents]=vcat(getindex.(sers,:parameterExponents))
    end
    ser[:charNumbers]=Cartesian(getindex.(sers,:charNumbers)...)
    ser[:cuspidalName]=join(getindex.(sers,:cuspidalName),"\\otimes ")
    ser
  end

  tt=refltype(WF)
  if isempty(tt) # UnipotentCharacters(coxgroup())
    return UnipotentCharacters(
      [Dict(:relativeType=>Dict[], 
	    :levi=>Int[], :parameterExponents=>Int[],
	    :cuspidalName=>"", :eigenvalue=>1, :charNumbers =>[ 1 ])],
      [Dict(:relativeType=>Dict[], 
	    :levi=>Int[], :parameterExponents=>Int[],
	    :cuspidalName=>"", :eigenvalue=>1, :charNumbers =>[ 1 ])],
     [Family("C1",[1])],
     Dict( :charParams => [ [ "", [ 1 ] ] ],
      :TeXCharNames => [ "." ],
      :almostTeXCharNames => [ "." ],
      :charSymbols => [ [ "", [ 1 ] ] ],
      :size=>1,
      :a => [ 0 ],
      :A => [ 0 ],
      :group=>WF))
  end

  W=WF.W
  simp=map(tt) do t
# adjust indices of Levis, almostLevis, relativetypes so they agree with
# Parent(Group(WF))
    uc=UnipotentCharacters(t)
    if uc==false return false end
    H=map(x->reflection_subgroup(W,x.indices),t.orbit)
    inc=vcat(map(x->x.indices,t.orbit)...)
    for s in uc.harishChandra
      s[:levi]=restriction(W,vcat(map(R->inclusion(R,s[:levi]),H)...))
      s[:relativeType].indices=restriction(W,inclusion(H[1],s[:relativeType].indices))
    end
    for s in uc.almostHarishChandra
      s[:levi]=restriction(W,vcat(map(R->inclusion(R,s[:levi]),H)...))
      s[:relativeType].orbit=vcat(map(x->
        map(s[:relativeType].orbit)do r
	  r=copy(r)
          r.indices=restriction(W,inclusion(x,r.indices))
	  r
        end,H)...)
      s[:relativeType].twist^=prod(map(Perm,1:length(inclusion(H[1])),inclusion(H[1])))
    end

    for f in uc.families
      f[:fourierMat]=fourier(f)
      if !haskey(f,:charLabels) 
        f[:charLabels]=string.(1:length(f[:eigenvalues]))
      end
    end
    uc
  end

  # "Kronecker product" of records in simp:
  r=simp[1]
  f=keys(r.prop)
  res=Dict{Symbol,Any}()
  for a in f
    if a!=:group 
    if length(simp)==1 
      res[a]=map(x->[x],r.prop[a])
    elseif all(x->haskey(x.prop,a),simp)
      res[a]=Cartesian(map(x->x.prop[a],simp)...)
    end
    end
  end
  
  for a in [:TeXCharNames,:almostTeXCharNames]
    res[a]=join.(res[a],"\\otimes ")
  end

  res[:size]=length(res[:TeXCharNames])
  
  # finally the new 'charNumbers' lists
  tmp=Cartesian(map(a->1:length(a.prop[:TeXCharNames]),simp)...)

  hh=CartesianSeries.(Cartesian(map(x->x.harishChandra,simp)...))
  ah=CartesianSeries.(Cartesian(map(x->x.almostHarishChandra,simp)...))
  for s in hh
    s[:charNumbers]=map(y->findfirst(isequal(y),tmp),s[:charNumbers])
  end
  for s in ah
    s[:charNumbers]=map(y->findfirst(isequal(y),tmp),s[:charNumbers])
  end

  if length(tt)==1
    ff=r.families
  else 
    ff=Family.(prod.(Cartesian(map(x->x.families,simp)...)))
    for f in ff
      f[:charNumbers]=map(y->findfirst(isequal(y),tmp),f[:charNumbers])
    end
  end
  
  for a in ["a", "A"]
    if haskey(res,a) res[a]=sum.(res[a]) end
  end


  res[:group]=WF
  UnipotentCharacters(hh,ah,ff,res)
  end
end

function Base.show(io::IO, ::MIME"text/html", uc::UnipotentCharacters)
  print(io, "\$")
  show(IOContext(io,:TeX=>true),uc)
  print(io, "\$")
end

Chars.charnames(io::IO,uc::UnipotentCharacters)=
   fromTeX.(Ref(io),uc.prop[:TeXCharNames])

function Base.show(io::IO,uc::UnipotentCharacters)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if !TeX print(io,"UnipotentCharacters(",uc.prop[:group],")") end
  if !repl && !TeX return end
  println(io,"")
  strip(x)=fromTeX(io,x)
  m=hcat(sprint.(show,CycPol.(degrees(uc)); context=io),
         sprint.(show,CycPol.(fakedegrees(uc)); context=io),
         sprint.(show,Root1.(eigen(uc)); context=io),
         strip.(labels(uc)))
  format(io,m,row_labels=charnames(io,uc),
         rows_label=strip("\\gamma"),
         col_labels=strip.(["Deg(\\gamma)","Feg","Fr(\\gamma)","label"]))
end

Groups.Group(uc::UnipotentCharacters)=uc.prop[:group]

Base.length(uc::UnipotentCharacters)=length(uc.prop[:TeXCharNames])

function Chars.fakedegrees(uc::UnipotentCharacters,q=Pol([1],1))
  if !haskey(uc.prop,:fakedegrees)
    uc.prop[:fakedegrees]=Dict{Any,Any}()
  end
  d=uc.prop[:fakedegrees]
  if haskey(d,q) return d[q] end
  f=fakedegrees(Group(uc),q)
  if isa(q,Pol) f=improve_type(f) end
  fd=fill(zero(f[1]),length(uc))
  fd[uc.almostHarishChandra[1][:charNumbers]]=f
  d[q]=fd
end

# FourierInverse times the vector of fake degrees is the vector of unip degrees
function fourierinverse(uc::UnipotentCharacters)
  gets(uc,:fourierinverse)do
     l=length(uc)
     T=reduce(promote_type,map(f->eltype(f[:fourierMat]),uc.families))
#    println(map(f->eltype(f[:fourierMat]),uc.families),"=> T=$T")
     i=fill(T(0),l,l)
     for f in uc.families
       i[f[:charNumbers],f[:charNumbers]]=f[:fourierMat]'
     end
     i
  end
end

function Families.fourier(uc::UnipotentCharacters)
  gets(uc,:fourier)do
     l=length(uc)
     i=fill(0*E(1)//1,l,l)
     for f in uc.families
       i[f[:charNumbers],f[:charNumbers]]=f[:fourierMat]
     end
     i
  end
end

"""
`degrees(uc::UnipotentCharacters,q=Pol([1],1))`

Returns  the  list  of  degrees  of  the unipotent characters of the finite
reductive group (or Spetses) with Weyl group (or Spetsial reflection group)
`W`, evaluated at `q`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> uc=UnipotentCharacters(W);

julia> degrees(uc)
10-element Array{Pol{Rational{Int64}},1}:
 1//1                                         
 (1//1)q⁶                                     
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
function Gapjm.degrees(uc::UnipotentCharacters,q=Pol([1],1))
  if !haskey(uc.prop,:degrees)
    uc.prop[:degrees]=Dict{Any,Any}()
  end
  d=uc.prop[:degrees]
  if haskey(d,q) return d[q] end
  d[q]=fourierinverse(uc)*fakedegrees(uc,q)
end

function eigen(uc::UnipotentCharacters)
  gets(uc,:eigen)do
    eig=fill(E(1),length(uc))
    for f in uc.families eig[f[:charNumbers]]=f[:eigenvalues] end
    eig
  end
end

function labels(uc::UnipotentCharacters)::Vector{String}
  gets(uc,:labels)do
    lab=fill("",length(uc))
    for f in uc.families lab[f[:charNumbers]]=f[:charLabels]
    end
    lab
  end
end

"""
fix illegal relativeTypes B1 and C2 which appear in HC or almost HC
series of classical groups
"""
function FixRelativeType(t)
  d=t[:relativeType]
  if d[:series]=="B" 
    if d[:rank]==1
      d[:series]="A"
      t[:charNumbers]=collect(t[:charNumbers]) # map B1->A1
      reverse!(view(t[:charNumbers],1:2)) # map B1->A1
    elseif d[:rank]==2 && haskey(d,:cartanType) && d[:cartanType]==1
      d[:cartanType]=2
      d[:indices]=reverse(collect(d[:indices]))
      reverse!(view(t[:charNumbers],[1,5])) # map C2->B2
      if haskey(t,:parameterExponents) reverse!(t[:parameterExponents]) end
    end
  end
end

"""
`CycPolUnipotentDegrees(W)`

Taking  advantage that  the degrees  of unipotent  characters of the finite
reductive group (or Spetses) with Weyl group (or Spetsial reflection group)
`W`  are products  of cyclotomic  polynomials, this  function returns these
degrees as a list of `CycPol`s.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> Uch.CycPolUnipotentDegrees(W)
10-element Array{CycPol{Rational{Int64}},1}:
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
CycPolUnipotentDegrees(W)=CycPol.(degrees(UnipotentCharacters(W)))
# slow implementation

#-------------------------- UniChars -------------------------------
struct UniChar{T,T1}
  group::T
  v::T1
  prop::Dict{Symbol,Any}
end

UniChar(W,v::AbstractVector)=UniChar(W,v,Dict{Symbol,Any}())

"""
`UniChar(W,l)`

Constructs  an object representing the unipotent character specified by `l`
of  the algebraic  group associated  to the  Coxeter group or Coxeter coset
specified  by `W`. There are 3 possibilities  for `l`: if it is an integer,
the  `l`-th unipotent character of `W` is  returned. If it is a string, the
unipotent  character of `W` whose name is `l` is returned. Finally, `l` can
be  a  list  of  length  the  number  of unipotent characters of `W`, which
specifies the coefficient to give to each.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> UniChar(W,7)
[G₂]:<G₂[-1]>

julia> UniChar(W,"G2[E3]")
[G₂]:<G₂[ζ₃]>

julia> UniChar(W,[1,0,0,-1,0,0,2,0,0,1])
[G₂]:<φ₁‚₀>-<φ″₁‚₃>+2<G₂[-1]>+<G₂[ζ₃²]>
```
"""
function UniChar(W,v::Int)
  r=zeros(Int,length(UnipotentCharacters(W)))
  r[v] = 1
  UniChar(W,r)
end

function UniChar(W,v::String)
  n=charnames(stdout,UnipotentCharacters(W))
  UniChar(W,findfirst(==(v),n))
end

const short=Ref(true)

function Base.show(io::IO,r::UniChar)
  res=""
  s=charnames(io,UnipotentCharacters(r.group))
  m=maximum(length.(s))+3
  for i = 1:length(r.v)
    n = "<"*s[i]*">"
    c = sprint(show,r.v[i];context=io)
    if short[]
      if c != "0"
        if c == "1" res*= "+"
        elseif c == "-1" res*="-"
        else
          if occursin(r".[+-]",c) c = "("* c* ")" end
          if !(c[1] in "+-") res*="+" end
          res*=c
        end
        res*=n
      end
    elseif c!="0" || !get(io,:nozero,false)
      res *= "\n"* rpad(n,m)* c
    end
  end
  if length(res) == 0 res = "0" end
  if res[1] == '+' res = res[2:end] end
  if haskey(r.prop, :name)
    res="DLvar["*sprint(show,r.group; context=io)*","*r[:name],"]:",res
  else
    res="["*sprint(show,r.group; context=io)*"]:"* res
  end
  print(io,res)
end

Base.:+(u1::UniChar,u2::UniChar)=UniChar(u1.group,u1.v+u2.v)
Base.:-(u1::UniChar,u2::UniChar)=UniChar(u1.group,u1.v-u2.v)
Base.:*(u1::UniChar,u2::UniChar)=sum(u1.v .* conj.(u2.v))
Base.:*(u1::UniChar,a)=UniChar(u1.group,u1.v .* a)
Base.:*(a,u1::UniChar)=u1*a

Gapjm.degree(u::UniChar,q=Pol(:q))=sum(u.v .*
                             degrees(UnipotentCharacters(u.group),q))

"""
`LusztigInduce(W,u)`

`u`  should be a unipotent character of a parabolic subcoset of the Coxeter
coset  `W`. It represents  a unipotent character  `λ` of a  Levi `𝐋` of the
algebraic  group  `𝐆`  attached  to  `W`.  The  program returns the Lusztig
induced `R_𝐋^𝐆(λ)`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> WF=spets(W)
G₂

julia> T=subspets(WF,Int[],W(1))
.Φ₁Φ₂

julia> u=UniChar(T,1)
[.Φ₁Φ₂]:<.>

julia> LusztigInduce(WF,u)
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> DLChar(W,W(1))
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>
```
"""
function LusztigInduce(WF, u)
  t=LusztigInductionTable(u.group, WF)
  if !isnothing(t) UniChar(WF, t.scalar*u.v) end
end

"""
`LusztigRestrict(R,u)`

`u`  should be a unipotent character of a parent Coxeter coset `W` of which
`R` is a parabolic subcoset. It represents a unipotent character `γ` of the
algebraic  group `𝐆` attached to `W`,  while `R` represents a Levi subgroup
`L`. The program returns the Lusztig restriction `*R_𝐋^𝐆(γ)`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> WF=spets(W)
G₂

julia> T=subspets(WF,Int[],W(1))
.Φ₁Φ₂

julia> u=DLChar(W,W(1))
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> Uch.LusztigRestrict(T,u)
[.Φ₁Φ₂]:4<.>

julia> T=subspets(WF,Int[],W(2))
.Φ₁Φ₂

julia> Uch.LusztigRestrict(T,u)
[.Φ₁Φ₂]:0
```
"""
LusztigRestrict(HF,u)=UniChar(HF,permutedims(LusztigInductionTable(HF,
                                                          u.group).scalar)*u.v)

HCInduce(WF,u)=UniChar(WF,HCInductionTable(u.group,WF).scalar*u.v)

HCRestrict(HF,u)=UniChar(HF,u.v*HCInductionTable(HF,u.group).scalar)

function DLCharTable(W)
  gets(W,:rwTable)do
    uc=UnipotentCharacters(W)
    CharTable(W).irr'*fourier(uc)[uc.almostHarishChandra[1][:charNumbers],:]
  end
end

"""
`DLChar(W,w)`

This  function returns  the Deligne-Lusztig  character `R_𝐓  ^𝐆 (1)` of the
algebraic  group `𝐆 ` associated to the Coxeter group or Coxeter coset `W`.
The  torus  `𝐓`  can  be  specified  in  3  ways:  if `w` is an integer, it
represents the `w`-th conjugacy class (or `phi`-conjugacy class for a coset
`Wϕ`)  of `W`. Otherwise  `w` can be  a word or  an element of  `W`, and it
represents the class (or `ϕ`-class) of `w`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> DLChar(W,3)
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> DLChar(W,W(1))
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> DLChar(W,[1])
[G₂]:<φ₁‚₀>-<φ₁‚₆>-<φ′₁‚₃>+<φ″₁‚₃>

julia> DLChar(W,[1,2])
[G₂]:<φ₁‚₀>+<φ₁‚₆>-<φ₂‚₁>+<G₂[-1]>+<G₂[ζ₃]>+<G₂[ζ₃²]>
```
"""
DLChar(W,i::Int)=UniChar(W,DLCharTable(W)[i,:])

DLChar(W,w::Perm)=DLChar(W,position_class(W,w))

DLChar(W,w::Vector{Int})=DLChar(W,W(w...))

"""
`AlmostChar(W,i)`

This  function  returns  the  `i`-th  almost  unipotent  character  of  the
algebraic  group 𝐆 associated to the Coxeter group or Coxeter coset `W`. If
`φ` is the `i`-th irreducible character of `W`, the `i`-th almost character
is  `R_φ=W⁻¹∑_w∈  W  φ(w)  R_𝐓_w^𝐆  (1)`  where  `𝐓_w` is the maximal torus
associated  to the conjugacy class (or  `ϕ`-conjugacy class for a coset) of
`w`.

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> AlmostChar(W,3)
[B₂]:<.11>

julia> AlmostChar(W,1)
[B₂]:1/2<11.>+1/2<1.1>-1/2<.2>-1/2<B₂>
```
"""
AlmostChar=function(W,i)
  ct=CharTable(W)
  dl=DLChar.(Ref(W),1:length(ct.charnames))
  sum(ct.irr[i,:] .* classes(ct).//length(W).*dl)
end

"""
`DLLefschetz(h)`

Here `h` is an element of a Hecke algebra associated to a Coxeter group <W>
which  itself  is  associated  to  an  algebraic  group `𝐆 `. By results of
Digne-Michel,  for `g∈  𝐆 ^F`,  the number  of fixed  points of `Fᵐ` on the
Deligne-Lusztig variety associated to the element `wϕ` of the Coxeter coset
`Wϕ`, have for `m` sufficiently divisible, the form `∑_φ φ_(qᵐ)(T_wϕ)R_φ(g)`
where  `φ` runs over the irreducible characters of `Wϕ`, where `R_φ` is the
corresponding  almost character, and where `φ_(qᵐ)` is a character value of
the  Hecke algebra `ℋ (Wϕ,qᵐ)` of `Wϕ` with parameter `qᵐ`. This expression
is  called the *Lefschetz character* of  the Deligne-Lusztig variety. If we
consider `qᵐ` as an indeterminate `x`, it can be seen as a sum of unipotent
characters  with coefficients character values of the generic Hecke algebra
`ℋ (Wϕ,x)`.

The  function 'DLLefschetz' takes  as argument a  Hecke element and returns
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

julia> DLLefschetz(T(1,2))
[A₂]:<111>-q<21>+q²<3>

julia> DLLefschetz((T(1)+T())*(T(2)+T()))
[A₂]:q<21>+(q²+2q+1)<3>
```

The   last  line  shows  the   Lefschetz  character  of  the  Samelson-Bott
desingularisation of the Coxeter element Deligne-Lusztig variety.

We now show an example with a coset (corresponding to the unitary group).

gap> H:=Hecke(CoxeterCoset(W,(1,2)),q^2);
Hecke(2A2,q^2)
gap> T:=Basis(H,"T");
function ( arg ) ... end
gap> DeligneLusztigLefschetz(T(1));
[2A2]=-<11>-q<2A2>+q^2<2>
"""
DLLefschetz=function(h,i=0)
# if haskey(h, :coset) W = ReflectionCoset(h[:coset])
# else 
  W=h.H.W
# end
  uc=UnipotentCharacters(W)
  UniChar(W, fourier(uc)[:,uc.harishChandra[1][:charNumbers]]*
          conj.(char_values(h)).*Uch.eigen(uc).^i)
end

DLLefschetzTable=function(H)
# if haskey(H, :spets) WF = ReflectionCoset(H)
# else 
  WF=H.W
# end
  t=CharTable(H).irr
  uc=UnipotentCharacters(WF)
  return t'*fourier(uc)[uc.harishChandra[1][:charNumbers],:]
end

Frobenius=function(WF, x::UniChar, i)
  W=x.group
  p=Perm(map(x->position_class(W,x^WF.phi), class_reps(W)))
  uc=UnipotentCharacters(W)
  t=vcat(DLCharTable(W), permutedims(eigen(uc)))
  pt=t^p
  p = map(findall(i->x==t[:;i], axes(t,2)), eachcol(pt))
  if any(x->length(x)>1,p) error("Rw + eigen cannot disambiguate\n") end
  p=Perm(map(x->x[1], p))
  UniChar(W,x.v^(p^-i))
end

"""
`+`: Adds the specified characters.

`-`: Subtracts the specified characters

`*`:  Multiplies  a  character  by  a  scalar,  or  if  given two unipotent
characters returns their scalar product.

We go on from examples of the previous section:

|    gap> u+v;
    [G2]=<G2[-1]>+<G2[E3]>
    gap> w-2*u;
    [G2]=<phi{1,0}>-<phi{1,3}''>+<G2[E3^2]>
    gap> w*w;
    7|

`:Degree`: returns the degree of the unipotent character.

|    gap> Degree(w);
    q^5 - q^4 - q^3 - q^2 + q + 1
    gap> Degree(u+v);
    (5/6)*q^5 + (-1/2)*q^4 + (-2/3)*q^3 + (-1/2)*q^2 + (5/6)*q|

`:String' and 'Print`: the formatting of unipotent characters is affected by
the  variable 'CHEVIE.PrintUniChars'. It is a  record; if the field 'short'
is  bound (the default)  they are printed  in a compact  form. If the field
`:long' is bound, they are printed one character per line:

|    gap> CHEVIE.PrintUniChars:=rec(long:=true);
    rec(
      long := true )
    gap> w;
    [G2]=
    <phi{1,0}>   1
    <phi{1,6}>   0
    <phi{1,3}'>  0
    <phi{1,3}''> -1
    <phi{2,1}>   0
    <phi{2,2}>   0
    <G2[-1]>     2
    <G2[1]>      0
    <G2[E3]>     0
    <G2[E3^2]>   1
    gap> CHEVIE.PrintUniChars:=rec(short:=true);;|

`Frobenius(  <WF> )`: If 'WF' is a  Coxeter coset associated to the Coxeter
group  `W`, the function 'Frobenius(WF)' returns  a function which does the
corresponding automorphism on the unipotent characters

|    gap> W:=CoxeterGroup("D",4);WF:=CoxeterCoset(W,(1,2,4));
    CoxeterGroup("D",4)
    3D4
    gap> u:=UnipotentCharacter(W,2);
    [D4]=<11->
    gap> Frobenius(WF)(u);
    [D4]=<.211>
    gap> Frobenius(WF)(u,-1);
    [D4]=<11+>|

"""
end
