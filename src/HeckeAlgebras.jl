"""
This   module  implements  Hecke  algebras  associated  to  finite  complex
reflection  groups and arbitrary Coxeter  groups (these algebras are called
Iwahori-Hecke  algebras  in  this  last  case),  and  also  implements  the
character  tables, Schur elements and representations of Hecke algebras for
finite  groups. For Iwahori-Hecke algebras  and `G(d,1,1)` this module also
implements  the standard `T` basis;  see the module `KL`for Kazhdan-Lusztig
bases.

Let  `(W,S)` be  a Coxeter  system and  let `mₛₜ`  be the order of `st` for
`s,t∈ S`. Let `R` be a commutative ring with 1 and for `s∈ S` let `uₛ₀,uₛ₁∈
R` be elements which depend only on the conjugacy class of `s` in `W` (this
is  the  same  as  requiring  that  `uₛᵢ=uₜᵢ`  whenever  `mₛₜ` is odd). The
Iwahori-Hecke   algebra  of  `W`  over  `R`  with  parameters  `uₛᵢ`  is  a
deformation  of the group algebra of `W` over `R` defined as follows: it is
the  unitary  associative  `R`-algebra  generated  by  elements  `Tₛ, s∈ S`
subject to the relations:

``(Tₛ-uₛ₀)(Tₛ-uₛ₁)=0`` for all `s∈ S` (the quadratic relations)

``TₛTₜTₛ…= TₜTₛTₜ…`` with `mₛₜ` factors on each side (the braid relations)

If  `uₛ₀=1` and  `uₛ₁=-1` for  all `s`  then the quadratic relations become
`Tₛ²=1` and the deformation of the group algebra is trivial.

Since  the generators `Tₛ`  satisfy the braid  relations, `H` is  in fact a
quotient  of the group algebra of the  braid group associated with `W`. The
braid relations also imply that for any reduced expression `s_1⋯ s_m` of `w
∈  W` the product `Tₛ_1⋯ Tₛ_m` has the same value, that we denote `T_w`. We
have `T_1=1`; if one of the `uₛᵢ` is invertible, the `{T_w}_{w∈ W}` form an
`R`-basis  of the Iwahori-Hecke algebra  which specializes to the canonical
basis of the group algebra `R[W]` for `uₛ₀↦1` and `uₛ₁↦-1`.

When  one  of  the  `uₛᵢ`  is  invertible,  the  structure  constants  (the
decomposion  of  a  product  `T_vT_w`)  in  the `T_w` basis are obtained as
follows.  Choose a reduced expression for `v`,  say `v=s_1 ⋯ s_k` and apply
inductively the formula:

``T_sT_w=T_{sw}``               if `l(sw)=l(w)+1`

``T_sT_w=-uₛ₀uₛ₁T_{sw}+(uₛ₀+uₛ₁)T_w`` if `l(sw)=l(w)-1`.

If  one of `uₛ₀` or `uₛ₁` is invertible  in `R`, for example `uₛ₁`, then by
changing  the generators  to `T′ₛ=-Tₛ/uₛ₁`,  and setting `qₛ=-uₛ₀/uₛ₁`, the
braid  relations do no change  (since when `mₛₜ` is  odd we have `uₛᵢ=uₜᵢ`)
but  the quadratic relations become `(T′ₛ-qₛ)(T′ₛ+1)=0`. This normalisation
is  the most common form considered  in the literature. Another common form
in  the context of  Kazhdan-Lusztig theory, is  `uₛ₀=√qₛ` and `uₛ₁=-√qₛ⁻¹`.
The  form provided, with two parameters per generator, is often useful, for
instance  when constructing  the Jones  polynomial. If  for all `s` we have
`uₛ₀=q`,   `uₛ₁=-1`   then   we   call   the   corresponding   algebra  the
"one-parameter" or "Spetsial" Iwahori-Hecke algebra associated with `W`.

For  some  Iwahori-Hecke  algebras  the  character  table,  and  in general
Kazhdan-Lusztig  bases, require  a square  root of  `-uₛ₀uₛ₁`. These square
roots  can be specified  with the keyword  `rootpara` when constructing the
algebra;  after  this  the  function  `rootpara(H)`  will return the chosen
roots. If not specified, we try to extract roots automatically when needed;
`rootpara(H)`  informs  on  the  choices  made. Note that some mathematical
results  require an explicit choice of one  of the two possible roots which
cannot be automatically made thus require a keyword initialisation.

There  is a universal choice  for `R` and `uₛᵢ`:  Let `uₛᵢ:s∈ S,i∈[0,1]` be
indeterminates   such  that  `uₛᵢ=uₜᵢ`  whenever  `mₛₜ`  is  odd,  and  let
`A=ℤ[uₛᵢ]` be the corresponding polynomial ring. Then the Hecke algebra `H`
of  `W` over `A` with parameters `uₛᵢ` is called the *generic Iwahori-Hecke
algebra*  of  `W`.  Any  Hecke  algebra  `H₁`  with parameters `vₛᵢ` can be
obtained  by  specialization  from  `H`,  since  there  is  a  unique  ring
homomorphism  `f:A → R` such that `f(uₛᵢ)=vₛᵢ` for all `i`. Then via `f` we
can identify `H₁` to ``R⊗ _A H``.

Certain invariants of the irreducible characters of the one-parameter Hecke
algebra  play a special role in the representation theory of the underlying
finite  Coxeter  groups,  namely  the  `a`-  and  `A`-invariants. For basic
properties   of  Iwahori-Hecke   algebras  and   their  relevance   to  the
representation theory of finite groups of Lie type, see for example [Curtis
and Reiner 1987](biblio.htm#CR87) Sections~67 and 68.

In  the  following  example,  we  compute  the multiplication table for the
`0`-Iwahori--Hecke algebra associated with the Coxeter group of type `A_2`.

```julia-repl
julia> W=coxgroup(:A,2)
A₂

julia> H=hecke(W,0)            # One-parameter algebra with `q=0`
hecke(A₂,0)

julia> T=Tbasis(H);            # Create the `T` basis

julia> b=T.(elements(W))       # the basis
6-element Vector{HeckeTElt{HeckeAlgebra{Int64, Perm{Int16}, FiniteCoxeterGroup{Perm{Int16},Int64}}, Int64, Perm{Int16}}}:
 T.
 T₁
 T₂
 T₁₂
 T₂₁
 T₁₂₁

julia> b*permutedims(b)       # multiplication table
6×6 Matrix{HeckeTElt{HeckeAlgebra{Int64, Perm{Int16}, FiniteCoxeterGroup{Perm{Int16},Int64}}, Int64, Perm{Int16}}}:
 T.    T₁     T₂     T₁₂    T₂₁    T₁₂₁
 T₁    -T₁    T₁₂    -T₁₂   T₁₂₁   -T₁₂₁
 T₂    T₂₁    -T₂    T₁₂₁   -T₂₁   -T₁₂₁
 T₁₂   T₁₂₁   -T₁₂   -T₁₂₁  -T₁₂₁  T₁₂₁
 T₂₁   -T₂₁   T₁₂₁   -T₁₂₁  -T₁₂₁  T₁₂₁
 T₁₂₁  -T₁₂₁  -T₁₂₁  T₁₂₁   T₁₂₁   -T₁₂₁
```
Thus,  we work  with algebras  with arbitrary  parameters. We will see that
this also works on the level of characters and representations.

For  general complex reflection  groups, the picture  is similar. The Hecke
algebras  are deformations  of the  group algebras,  generalizing those for
real reflection groups.

The  definition is  as a  quotient of  the algebra  of the  braid group. We
assume  now that `W` is  a *finite* reflection group  in the complex vector
space  `V`. The *braid group* associated  is the fundamental group `Π₁` of
the  space  ``(V-\\bigcup_{H\\in\\mathcal H}  H)/W``, where ``\\mathcal H``
is  the set of  reflecting hyperplanes of  `W`. This group  is generated by
*braid reflections*, elements which by the natural map from the braid group
to  the reflection  group project  to distinguished  reflections. The braid
reflections   which  project  to  a  given  `W`-orbit  of  reflections  are
conjugate.  Let `𝐬` be a representative of  such a conjugacy class of braid
reflections,  let `e`  be the  order of  the image  of `𝐬`  in `W`, and let
``u_{𝐬,0},…,u_{𝐬,e-1}`` be indeterminates. The generic Hecke algebra of `W`
is  the  ``ℤ[u_{𝐬,i}^{±  1}]_{𝐬,i}``-algebra  quotient  of  the braid group
algebra  by the relations ``(𝐬-u_{𝐬,0})…(𝐬-u_{𝐬,e-1})=0``, and an arbitrary
Hecke  algebra for `W` is an algebra  obtained from this generic algebra by
specializing some of the parameters.

The  generic Hecke algebras are explicitely  described by a presentation of
the  braid group. The braid group can be presented by homogeneous relations
in   the  braid   reflections,  called   *braid  relations*,  described  in
[Broué-Malle-Rouquier     1998](biblio.htm#BMR98)     and    [Bessis-Michel
2003](biblio.htm#BM03)  (some  of  which  were  obtained  using the VKCURVE
GAP3-package,  also ported to Julia).  Furthermore, these presentations are
such  that the  reflection group  is presented  by the same relations, plus
relations  describing the order  of the generating  reflections, called the
*order  relations*. Thus  the Hecke  algebra has  a presentation similar to
that of `W`, with the same braid relations but the order relations replaced
by a deformed version.

If  `S⊂ W`  is the  set of  distinguished reflections  of `W` which lift to
generating  braid reflections in the braid  group, for each conjugacy class
of  an  `s`  of  order  `e`  we take indeterminates `uₛ₀,…,uₛₑ₋₁`. Then the
generic  Hecke algebra is the ``ℤ[uₛᵢ^{±1}]ₛᵢ``-algebra `H` with generators
`T_s`  for each `s∈  S` presented by  the braid relations  and the deformed
order relations ``(T_s-u_{s,0})…(T_s-u_{s,e-1})=0``.

Ariki,  Koike and Malle have computed the  character table of some of these
algebras,  including  those  for  all  2-dimensional reflection groups, see
[Broué-Malle 1993](biblio.htm#BM93) and [Malle 1996](biblio.htm#Mal96); our
data  has  models  of  all  representation  and  character  tables for real
reflection  groups; it  contains the  same for  imprimitive groups  and for
primitive groups of dimension 2 and 3 (these last representations have been
computed  in [Malle-Michel 2010](biblio.htm#MM10)) and contains also models
and  character tables computed  by Michel for  `G₂₉` and `G₃₃`; it contains
also  partial lists of representations and partial character tables for the
remaining  groups `G₃₁,G₃₂`  and `G₃₄`,  computed by  Malle and  Michel for
`G₃₂` and by Michel for the other two algebras.

The quotient of the Hecke algebra obtained by the specialisation ``u_{𝐬,i}↦
ζₑⁱ``  is isomorphic to the group algebra of `W`. It was conjectured for 20
years  that over a splitting ring the Hecke algebra is itself isomorphic to
the  group algebra of `W` over the  same ring. This was called the freeness
conjecture since the main problem is to show that the Hecke algebra is free
of dimension `|W|`. This has finally been proved in 2019 thanks to the work
of  many  people  including  Marin,  Pfeiffer,  Chavli  and  Tsuchioka  for
exceptional  groups. Along the way  it has been proven  that there exists a
set  `{b_w}_{w∈ W}` of  elements of the  Braid group such  that `b_1=1` and
`b_w` maps to `w` by the natural quotient map, such that their images `T_w`
form a basis of the Hecke algebra.

It  is  conjectured  that  such  a  basis  `T_w`  can  be  chosen such that
additionnaly  the  linear  form  `t`  defined  by  `t(T_w)=0` if `w≠ 1` and
`t(1)=1` is a symmetrizing form for the symmetric algebra `H`. This is well
known  for all real reflection groups  and has been proved in [Malle-Mathas
1998](biblio.htm#MM98)   for   imprimitive   reflection   groups   and   in
[Malle-Michel 2010](biblio.htm#MM10) for some primitive groups of dimension
2  and  3.  Chlouveraki  and  Chavli  have handled some other 2-dimensional
cases.  For  each  irreducible  character  `φ`  of `H` we define the *Schur
element*  `Sᵩ` associated to `φ` by the  condition that for any element `T`
of  `H` we have `t(T)=∑ᵩ φ(T)/Sᵩ`. It  can be shown that the Schur elements
are  Laurent polynomials, and they  do not depend on  the choice of a basis
having  the  above  property.  Malle  has  computed  these  Schur elements,
assuming the above conjecture; they are in the Chevie data.

See the function `hecke` for various ways of specifying the parameters of a
Hecke   algebra.  Look  also  at   the  docstrings  of  `central_monomials,
char_values,     class_polynomials,    schur_elements,    isrepresentation,
factorized_schur_elements`,  and  at  the  methods  for  Hecke  algebras of
`CharTable, representations, reflrep`.

Taking  apart  Hecke  elements  is  done  with  the  functions  `getindex`,
`setindex!`, `keys`, `values`, `iterate`.

```julia-repl
julia> H=hecke(W,Pol(:q))
hecke(A₂,q)

julia> T=Tbasis(H);

julia> h=T(1,2)^2
qT₂₁+(q-1)T₁₂₁

julia> length(h) # h has 2 terms
2

julia> h[W(2,1)] # coefficient of W(2,1)
Pol{Int64}: q

julia> collect(h) # pairs perm=>coeff
2-element Vector{Any}:
  (1,2,6)(3,4,5) => q
 (1,5)(2,4)(3,6) => q-1

julia> collect(values(h)) # the coefficients
2-element Vector{Pol{Int64}}:
 q
 q-1

julia> collect(keys(h)) # the corresponding Perms
2-element Vector{Perm{Int16}}:
 (1,2,6)(3,4,5)
 (1,5)(2,4)(3,6)

julia> h[W(2,1)]=Pol(3)
Pol{Int64}: 3

julia> h
3T₂₁+(q-1)T₁₂₁
```
finally, benchmarks on julia 1.8
```benchmark
julia> function test_w0(n)
         W=coxgroup(:A,n)
         Tbasis(hecke(W,Pol(:q)))(longest(W))^2
       end
test_w0 (generic function with 1 method)

julia> @btime test_w0(7);
   97.210 ms (1776476 allocations: 127.52 MiB)
```
in GAP3 the analogous function takes 920ms
```
test_w0:=function(n)local W,T,H;
  W:=CoxeterGroup("A",n);H:=Hecke(W,X(Rationals));T:=Basis(H,"T");
  return T(LongestCoxeterWord(W))^2;
end;
```
"""
module HeckeAlgebras
using ..Chevie
export HeckeElt, Tbasis, central_monomials, hecke, HeckeAlgebra, HeckeTElt,
  rootpara, equalpara, class_polynomials, char_values, schur_elements,
  schur_element, isrepresentation, alt, HeckeCoset,
  FactorizedSchurElements, factorized_schur_element, VFactorSchurElement,
  factorized_schur_elements

@GapObj struct HeckeAlgebra{C,T,TW<:Group{T}}<:FiniteDimAlgebra{C}
  W::TW
  para::Vector{Vector{C}}
end

Algebras.dim(H::HeckeAlgebra)=length(H.W)
Algebras.basis(H::HeckeAlgebra,i::Integer)=Tbasis(H)(elements(H.W)[i])

"""
`hecke( W [, parameter];rootpara=nothing)`

Hecke  algebra for the complex reflection group or Coxeter group `W`. If no
`parameter` is given, `1` is assumed which gives the group algebra of `W`.

The  following forms are accepted for  `parameter`: if `parameter` is not a
vector or a tuple, it is replaced by the vector `fill(parameter,ngens(W))`.
If it is a vector with one entry, it is replaced with
`fill(parameter[1],ngens(W))`.  If `parameter`  is a  vector with more than
one  entry, it should  have length `ngens(W)`,  each entry representing the
parameters   for   the   corresponding   generator   of  `W`,  and  entries
corresponding  to  the  same  `W`-orbit  of generators should be identical.
Finally, if `parameter` is a `Tuple`, the tuple should have as many entries
as  there are hyperplane  orbits in `W`  and each entry  will represent the
parameters for the corresponding conjugacy class of braid reflections.

An  entry in  `parameter` for  a reflection  of order  `e` can  be either a
single scalar value or a `Vector` of length 'e'. If it is a `Vector`, it is
interpreted as the list `[u₀,…,u_(e-1)]` of parameters for that reflection.
If  it is not  a vector, let  `q` be its  value; it is  then interpreted as
specifying  the  list  of  parameters  for  the Spetsial algebra, which are
`[q,ζ_e,…,ζ_{e-1}]`  (thus the  list `[q,-1]`  of the one-parameter algebra
for Coxeter groups).

When  printing an Hecke algebra the parameter list is abbreviated using the
same conventions.

Computing characters or representations of Hecke algebra needs sometimes to
extract  roots of the  parameters. These roots  are extracted automatically
(when  possible). For Coxeter groups it  is possible to give explicit roots
by  giving  a  keyword  argument  `rootpara`:  if  it is a vector it should
contain at the `i`-th position a square root of
`-parameter[i][1]*parameter[i][2]`;   if  a   scalar  it   is  replaced  by
`fill(rootpara,ngens(W))`.

# Example
```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> @Pol q
Pol{Int64}: q

julia> H=hecke(W,q)
hecke(B₂,q)

julia> H.para
2-element Vector{Vector{Pol{Int64}}}:
 [q, -1]
 [q, -1]

julia> H=hecke(W,q^2,rootpara=-q)
hecke(B₂,q²,rootpara=-q)

julia> H=hecke(W,q^2)
hecke(B₂,q²)

julia> rootpara(H)
2-element Vector{Pol{Int64}}:
 q
 q

julia> H
hecke(B₂,q²,rootpara=q)

julia> H=hecke(W,[q^2,q^4],rootpara=[q,q^2])
hecke(B₂,Pol{Int64}[q², q⁴],rootpara=Pol{Int64}[q, q²])

julia> H.para,rootpara(H)
(Vector{Pol{Int64}}[[q², -1], [q⁴, -1]], Pol{Int64}[q, q²])

julia> H=hecke(W,9,rootpara=3)
hecke(B₂,9,rootpara=3)

julia> H.para,rootpara(H)
([[9, -1], [9, -1]], [3, 3])

julia> @Mvp x,y,z,t

julia> H=hecke(W,[[x,y]])
hecke(B₂,Vector{Mvp{Int64, Int64}}[[x, y]])

julia> rootpara(H);H
hecke(B₂,Vector{Mvp{Int64, Int64}}[[x, y]],rootpara=ζ₄x½y½)

julia> H=hecke(W,[[x,y],[z,t]])
hecke(B₂,Vector{Mvp{Int64, Int64}}[[x, y], [z, t]])

julia> rootpara(H);H
hecke(B₂,Vector{Mvp{Int64, Int64}}[[x, y], [z, t]],rootpara=Mvp{Cyc{Int64}, Rational{Int64}}[ζ₄x½y½, ζ₄t½z½])

julia> hecke(coxgroup(:F,4),(q,q^2)).para
4-element Vector{Vector{Pol{Int64}}}:
 [q, -1]
 [q, -1]
 [q², -1]
 [q², -1]

julia> hecke(complex_reflection_group(3,1,2),q).para # spetsial parameters
2-element Vector{Vector{Pol{Cyc{Int64}}}}:
 [q, ζ₃, ζ₃²]
 [q, -1]
```
"""
function hecke(W::Group,para::Vector{<:Vector{C}};rootpara::Vector=C[])where C
  if applicable(simple_reps,W)
  para=map(eachindex(gens(W)))do i
    j=simple_reps(W)[i]
    if i<=length(para)
     if j<i && para[i]!=para[j] error("one should have  para[$i]==para[$j]") end
      return para[i]
    elseif length(para)==1 return para[1]
    elseif j<i return para[j]
    else error("parameters should be given for first reflection in a class")
    end
  end
  end
  d=Dict{Symbol,Any}(:equal=>allequal(para))
  if !isempty(rootpara) d[:rootpara]=rootpara end
  p=findfirst(i->length(para[i])!=order(W(i)),eachindex(para))
  if !isnothing(p) error("parameter no. $p should be of length ",order(W(p))) end
  HeckeAlgebra(W,para,d)
end

function hecke(W::Group,p::Vector;rootpara::Vector=Any[])
  oo=ordergens(W)
  para=map(p,oo)do p, o
    if p isa Vector return p end
    all(==(2),oo) ? [p,-one(p)] : vcat([p],E.(o,1:o-1))
  end
  if isempty(para)
   return HeckeAlgebra(W,Vector{Int}[],Dict{Symbol,Any}(:rootpara=>rootpara))
  end
  hecke(W,para;rootpara=convert(Vector{eltype(para[1])},rootpara))
end

function hecke(W::Group,p::C=1;rootpara::C=zero(p))where C
  if ngens(W)==0 para=Vector{C}[]
  elseif all(==(2),ordergens(W)) para=[[p,-one(p)] for o in ordergens(W)]
  else para=map(o->vcat([p],E.(o,1:o-1)),ordergens(W))
  end
  H=HeckeAlgebra(W,para,Dict{Symbol,Any}(:equal=>true))
  if !iszero(rootpara) H.rootpara=fill(rootpara,ngens(W)) end
  H
end

function hecke(W::Group,p::Tuple;rootpara=nothing)
  if length(p)==1
    para=fill(p[1],ngens(W))
    rootpara= isnothing(rootpara) ? eltype(p[1])[] : rtpara
  else
    s=simple_reps(W,1:ngens(W))
    C=sort(unique(s))
    para=fill(first(p),ngens(W))
    if !isnothing(rootpara) rtpara=fill(first(rootpara),ngens(W)) end
    for i in 1:ngens(W)
      para[i]=p[findfirst(==(s[i]),C)]
      if !isnothing(rootpara)
        rtpara[i]=rootpara[findfirst(==(s[i]),C)]
      end
    end
    rootpara= isnothing(rootpara) ? eltype(p[1])[] : rtpara
  end
  hecke(W,para;rootpara=rootpara)
end

function rootpara(H::HeckeAlgebra)
  get!(H,:rootpara)do
    map(eachindex(H.para)) do i
       if isone(-prod(H.para[i])) return -prod(H.para[i]) end
       return root(-prod(H.para[i]))
    end
  end
end

equalpara(H::HeckeAlgebra)::Bool=H.equal

function simplify_para(para)
  tr(p)=all(i->p[i]==E(length(p),i-1),2:length(p)) ? p[1] : p
  if isempty(para) para
  elseif allequal(tr.(para))
    p=tr(para[1])
    p isa Vector ? [p] : p
  else map(tr,para)
  end
end

function Base.show(io::IO, H::HeckeAlgebra)
  if isempty(H.para) print(io,"hecke(",H.W,")"); return end
  print(io,"hecke(",H.W,",",simplify_para(H.para))
  if haskey(H,:rootpara)
    rp=rootpara(H)
    if !isempty(rp) && allequal(rp) print(io,",rootpara=",rp[1])
    else print(io,",rootpara=",rp)
    end
  end
  print(io,")")
end

"""
`CharTable(H::HeckeAlgebra or HeckeCoset)`

returns  the `CharTable` of the Hecke algebra `H`. For the primitive groups
`G₃₁,  G₃₂,  G₃₄`  there  are  `Unknown()` entries corresponding to missing
representations   (see  [`representation`](@ref)).   The  columns   of  the
`CharTable`  are labelled  by `classnames(H.W)`  and contain  the character
values for the corresponding element given by `classreps(H.W)`.

```julia-repl
julia> H=hecke(crg(4),Pol())
hecke(G₄,q)

julia> CharTable(H)
CharTable(hecke(G₄,q))
┌────┬──────────────────────────────────────┐
│    │.    z 212   12    z12     1        1z│
├────┼──────────────────────────────────────┤
│φ₁‚₀│1   q⁶  q³   q²     q⁸     q        q⁷│
│φ₁‚₄│1    1   1  ζ₃²    ζ₃²    ζ₃        ζ₃│
│φ₁‚₈│1    1   1   ζ₃     ζ₃   ζ₃²       ζ₃²│
│φ₂‚₅│2   -2   .    1     -1    -1         1│
│φ₂‚₃│2 -2q³   . ζ₃²q -ζ₃²q⁴ q+ζ₃² -q⁴-ζ₃²q³│
│φ₂‚₁│2 -2q³   .  ζ₃q  -ζ₃q⁴  q+ζ₃  -q⁴-ζ₃q³│
│φ₃‚₂│3  3q²  -q    .      .   q-1     q³-q²│
└────┴──────────────────────────────────────┘
```
"""
function Chars.CharTable(H::HeckeAlgebra;opt...)
  get!(H,:chartable)do
    W=H.W
    if isempty(refltype(W)) ct=prod(CharTable[])
    else
      ct=prod(map(refltype(W))do t
        ct=getchev(t,:HeckeCharTable,H.para[t.indices], haskey(H,:rootpara) ?
                 rootpara(H)[t.indices] : fill(nothing,length(H.para)))
        CharTable(improve_type(toM(ct[:irreducibles])),
                  charnames(t;opt...,TeX=true),
             string.(classnames(t;opt...,TeX=true)),Int.(ct[:centralizers]),
             Int(ct[:centralizers][1]),Dict{Symbol,Any}())
      end)
    end
    ct.name=xrepr(H;TeX=true)
    ct.group=H
    ct
  end::CharTable
end

using SparseArrays
"""
`representation(H::HeckeAlgebra or HeckeCoset,i)`

returns,  for the `i`-th irreducible representation of the Hecke algebra or
Hecke  coset `H`, a list  of matrices images of  the generators of `H` in a
model of the representation (for Hecke cosets, the result is a `NamedTuple`
with fields `gens`, a representation of `hecke(H)`, and `F`, the matrix for
the automorphism of `H` in the representation).

This  function  is  based  on  the  classification,  and  is  not yet fully
implemented for the Hecke algebras of the groups `G₃₁`, `G₃₂` and `G₃₄`: we
have 50 representations out of 59 for type `G₃₁`, 30 representations out of
102  for  type  `G₃₂`  and  38  representations  out of 169 for type `G₃₄`;
`nothing` is returned for a missing representation.

```julia-repl
julia> W=crg(24)
G₂₄

julia> H=hecke(W,Pol(:q))
hecke(G₂₄,q)

julia> representation(H,3)
3-element Vector{Matrix{Pol{Cyc{Int64}}}}:
 [q 0 0; -q -1 0; -q 0 -1]
 [-1 0 -1; 0 -1 ((1-√-7)/2)q; 0 0 q]
 [-1 -1 0; 0 q 0; 0 (1+√-7)/2 -1]
```

The  models  implemented  for  imprimitive  types `G(de,e,n)` for `n>2` and
`de>1` (this includes Coxeter type `Dₙ`), excepted for `G(2,2,4), G(3,3,3),
G(3,3,4), G(3,3,5)` and `G(4,4,3)`, involve rational fractions.

```julia-repl
julia> H=hecke(coxgroup(:D,5),Pol())
hecke(D₅,q)

julia> representation(H,7)
5-element Vector{Matrix{Frac{Pol{Int64}}}}:
 [q 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]
 [q 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]
 [1/(-q-1) q/(q+1) 0 0; (q²+q+1)/(q+1) q²/(q+1) 0 0; 0 0 -1 0; 0 0 0 -1]
 [-1 0 0 0; 0 1/(-q²-q-1) (-q²-q)/(-q²-q-1) 0; 0 (q³+q²+q+1)/(q²+q+1) q³/(q²+q+1) 0; 0 0 0 -1]
 [-1 0 0 0; 0 -1 0 0; 0 0 1/(-q³-q²-q-1) (-q³-q²-q)/(-q³-q²-q-1); 0 0 (q⁴+q³+q²+q+1)/(q³+q²+q+1) q⁴/(q³+q²+q+1)]
```
"""
function Chars.representation(H::HeckeAlgebra,i::Integer)
  dims=getchev(H.W,:NrConjugacyClasses)
  if isempty(dims) return Matrix{Int}[] end
  tt=refltype(H.W)
  rp=haskey(H,:rootpara) ? rootpara(H) : fill(nothing,length(H.para))
  mm=map((t,j)->getchev(t,:HeckeRepresentation,H.para[t.indices],rp,j),tt,
                                                    lin2cart(dims,i))
  if any(==(false),mm) return nothing end
  if !(mm[1][1] isa AbstractMatrix) mm=map(x->toM.(x),mm) end
  if !all(m->m isa Vector{<:SparseMatrixCSC},mm) mm=improve_type.(mm) end
  n=length(tt)
  if n==1 return mm[1] end
  vcat(map(1:n) do i
     map(mm[i]) do m
       kron(map(j->j==i ? m : mm[j][1]^0,1:n)...)
    end
  end...)
end

"""
`isrepresentation(H::HeckeAlgebra,r)`

returns `true` or `false`, according to whether a given set `r` of elements
corresponding  to  the  standard  generators  of the reflection group `H.W`
defines a representation of the Hecke algebra `H` or not.

```julia-repl
julia> H=hecke(coxgroup(:F,4))
hecke(F₄,1)

julia> isrepresentation(H,reflrep(H))
true

julia> isrepresentation(H,Tbasis(H).(1:4))
true
```
"""
function isrepresentation(H::HeckeAlgebra,t;verbose=false)
  W=H.W
  res=true
#   bug in sparsearrays: q*one(m) is of type Any for SparseMatrix m
  myone(m,q)=m isa AbstractMatrix ? Diagonal(fill(q,size(m,1))) : q*one(m)
  for i in eachindex(gens(W))
    if !iszero(prod(q->t[i]-myone(t[i],q+0),H.para[i]))
      if !verbose return false end
      println("Error in ",ordinal(i)," parameter relation");
      res=false
    end
  end
  for (l,r) in braid_relations(W)
    if !iszero(prod(t[l])-prod(t[r]))
      if !verbose return false end
      println("Error in relation ",l,"=",r)
      res=false
    end
  end
  res
end

"""
`reflection_representation(H::HeckeAlgebra)` or `reflrep(H)`

returns  a  list  of  matrices  for  the  generators  of `H` which give the
reflection representation of the Iwahori-Hecke algebra `H`.

```julia-repl
julia> W=coxgroup(:B,2);H=hecke(W,Pol(:q))
hecke(B₂,q)

julia> reflrep(H)
2-element Vector{Matrix{Pol{Int64}}}:
 [-1 0; -q q]
 [q -2; 0 -1]

julia> H=hecke(coxgroup(:H,3))
hecke(H₃,1)

julia> reflrep(H)
3-element Vector{Matrix{Cyc{Int64}}}:
 [-1 0 0; -1 1 0; 0 0 1]
 [1 (-3-√5)/2 0; 0 -1 0; 0 -1 1]
 [1 0 0; 0 1 -1; 0 0 -1]
```
"""
function PermRoot.reflection_representation(H::HeckeAlgebra)
  W=H.W
  if !equalpara(H) || !(W isa CoxeterGroup)
        error("Reflexion representation of Cyclotomic Hecke algebras or\n",
          "Hecke algebras with unequal parameters not implemented")
  end
  q=-1*H.para[1][1]//H.para[1][2]
  r=ngens(W)
  C=fill(q*E(1),r,r)
  CM=coxmat(W)
  for i  in eachindex(gens(W))
    for j  in 1:i-1
      m=CM[i,j]
      if m!=0 m=E(m)+E(m,-1) else m=2 end
      C[i,j]=2+m
      if m==-2 C[j,i]=0 else C[j,i]=q end
    end
    C[i,i]=q+1
  end
  improve_type(map(eachindex(gens(W)))do i
    a=fill(0*q*E(1),r,r)
    for j  in eachindex(gens(W))
      a[j,j]=q
      a[j,i]-=C[i,j]
    end
    -H.para[1][2]*a
  end)
end

"""
`WGraphToRepresentation(H::HeckeAlgebra,gr::Vector)`

`H`  should be  a one-parameter  Hecke algebra  for a  finite Coxeter group
where  `rootpara`  is  defined.  The  function  returns the matrices of the
representation  of `H` defined by the W-graph `gr`.

```julia-repl
julia> W=coxgroup(:H,3)
H₃

julia> H=hecke(W,Pol(:q)^2)
hecke(H₃,q²)

julia> g=Wgraph(W,3)
2-element Vector{Vector{Vector{Any}}}:
 [[2], [1, 2], [1, 3], [1, 3], [2, 3]]
 [[-1, [[1, 3], [2, 4], [3, 5], [4, 5]]]]

julia> WGraphToRepresentation(H,g)
3-element Vector{Matrix{Pol{Int64}}}:
 [q² 0 … 0 0; 0 -1 … 0 0; … ; 0 0 … -1 q; 0 0 … 0 q²]
 [-1 0 … 0 0; 0 -1 … q 0; … ; 0 0 … q² 0; 0 0 … q -1]
 [q² 0 … 0 0; 0 q² … 0 0; … ; 0 q … -1 0; 0 0 … 0 -1]
```
"""
function Chars.WGraphToRepresentation(H::HeckeAlgebra,gr::Vector)
  if !equalpara(H)
    error("cell representations for unequal parameters not yet implemented")
  end
  S=-H.para[1][2]*WGraphToRepresentation(length(H.para),gr,
                                         rootpara(H)[1]//H.para[1][2])
  if !isrepresentation(H,S;verbose=true) error() end
  improve_type(S)
end

"""
`central_monomials(H)`

Let  `H`  be  an  Hecke  algebra  for  the finite reflection group `W`. The
function  returns the scalars by which the image  in `H` of `π` acts on the
irreducible representations of `H`.

When  `W` is irreducible,  `π` is the  generator of the  center of the pure
braid  group.  In  general,  it  is  the  product of such elements for each
irreducible  component. When  `W` is  a Coxeter  group, the  image of  π in
`H` is ``T_{w_0}^2``.

```julia-repl
julia> H=hecke(coxgroup(:H,3),Pol(:q))
hecke(H₃,q)

julia> central_monomials(H)
10-element Vector{Pol{Cyc{Int64}}}:
 1
 q³⁰
 q¹²
 q¹⁸
 q¹⁰
 q¹⁰
 q²⁰
 q²⁰
 q¹⁵
 q¹⁵
```
"""
function central_monomials(H::HeckeAlgebra)
# Cf. BrMi, 4.16 for the formula used
  W=H.W
  v=hyperplane_orbits(W)
  map(eachrow(CharTable(W).irr)) do irr
    dim=Int(irr[1])
    prod(v)do C
      q=H.para[restriction(W)[C.s]]
      m=Int.(map(0:C.order-1)do j
       (dim+sum(l->irr[C.cl_s[l]]*E(C.order,-j*l),1:C.order-1))//C.order
      end)
      E.(dim,-C.N_s*sum(m.*(0:C.order-1)))*
          prod(j->q[j]^Int(C.N_s*C.order*m[j]//irr[1]),1:C.order)
    end
  end
end

#--------------------------------------------------------------------------
# TH= typeof Algebra P=typeof(keys) [Perms] C typeof(coeffs)
abstract type HeckeElt{TH, C, P} end 

Base.zero(h::HeckeElt)=clone(h,zero(h.d))
Base.iszero(h::HeckeElt)=iszero(h.d)
Base.:(==)(a::HeckeElt,b::HeckeElt)=a.H===b.H && a.d==b.d
Base.copy(h::HeckeElt)=clone(h,h.d)
Base.getindex(a::HeckeElt{TH,C,P},p::P) where {TH,C,P}=a.d[p]
Base.setindex!(a::HeckeElt{TH,C,P},c::C,p::P) where {TH,C,P}=setindex!(a.d,c,p)
@inline Base.iterate(a::HeckeElt,s...)=iterate(a.d,s...)
Base.length(a::HeckeElt)=length(a.d)
Base.keys(a::HeckeElt)=keys(a.d)
Base.values(a::HeckeElt)=values(a.d)

# HeckeElts are scalars for broadcasting
Base.broadcastable(h::HeckeElt)=Ref(h)

function Base.show(io::IO, h::HeckeElt)
  function showbasis(io::IO,e)
    w=word(h.H.W,e)
    res=basisname(h)
    if hasdecor(io) res*=isempty(w) ? "." : "_"*joindigits(w,"{}";always=true)
    else            res*="("*join(w,",")*")"
    end
    fromTeX(io,res)
  end
  show(IOContext(io,:showbasis=>showbasis),h.d)
end


Base.:+(a::HeckeElt, b::HeckeElt)=clone(a,a.d+b.d)
Base.:-(a::HeckeElt)=clone(a,-a.d)
Base.:-(a::HeckeElt, b::HeckeElt)=clone(a,a.d-b.d)

Base.:*(a::HeckeElt, b::Union{Number,Pol,Mvp})=clone(a,a.d*b)
Base.:*(b::Union{Number,Pol,Mvp}, a::HeckeElt)=a*b

Base.:^(a::HeckeElt, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                        Base.power_by_squaring(inv(a),-n)
#--------------------------------------------------------------------------
const MM=ModuleElt # HModuleElt is 3 times slower

struct HeckeTElt{TH,C,P}<:HeckeElt{TH,C,P}
  d::MM{P,C}
  H::TH
end

function Groups.gens(H::HeckeAlgebra{C,T,TW})where {C,T,TW}
  get!(H,:gens)do
    map(i->Tbasis(H,H.W(i)),1:ngens(H.W))
  end::Vector{HeckeTElt{HeckeAlgebra{C,T,TW},C,T}}
end

clone(h::HeckeTElt,d)=HeckeTElt(d,h.H) # d could be different type from h.d
basisname(h::HeckeTElt)="T"

function Base.one(H::HeckeAlgebra{C,T,TW})where {C,T,TW}
  get!(H,:one)do
    Tbasis(H,one(H.W))
  end::HeckeTElt{HeckeAlgebra{C,T,TW},C,T}
end

Base.one(h::HeckeTElt)=one(h.H)
Base.zero(H::HeckeAlgebra{C,T}) where {C,T} =HeckeTElt(zero(MM{T,C}),H)

"""
`Tbasis(H::HeckeAlgebra)`
The  `T` basis of  `H`. It is  defined currently for Iwahori-Hecke algebras
and  for Hecke algebras of cyclic  complex reflection groups `G(d,1,1)`. It
returns  a function, say `T`,  which can take an  argument of the following
forms

  - `T(i::Integer)`: the generator `T_s` where `s=H.W(i)`.
  - `T(i₁,…,iᵣ)`: the product `T(i₁)…T(iᵣ)`
  - `T([i₁,…,iᵣ])`: same as `T(i₁,…,iᵣ)`
  - `T(w)` where `w∈ H.W`: returns `T_w`

```julia-repl
julia> H=hecke(coxgroup(:A,2),Pol(:q))
hecke(A₂,q)

julia> T=Tbasis(H);T(longest(H.W))^2
q³T.+(q³-2q²+q)T₂₁+(q³-q²)T₂+(q³-q²)T₁+(q³-2q²+2q-1)T₁₂₁+(q³-2q²+q)T₁₂

julia> W=crg(3,1,1)
G₃‚₁‚₁

julia> H=hecke(crg(3,1,1),Pol(:q))
hecke(G₃‚₁‚₁,q)

julia> T=Tbasis(H);T(1)^3
(q-1)T.+(q-1)T₁+qT₁₁
```
"""
Tbasis(H::HeckeAlgebra)=(x...)->x==() ? one(H) : Tbasis(H,x...)
function Tbasis(H::HeckeAlgebra{C,T,TW},w::Vector{<:Integer}) where {C,T,TW<:Group{T}}
  ww=H.W(w...)::T
  if length(w)==1 && w[1]>0 Tbasis(H,ww)
  elseif all(>(0),w) && H.W isa CoxeterGroup && length(H.W,ww)==length(w) 
    Tbasis(H,ww)
  else prod(i->i>0 ? Tbasis(H,H.W(i)) : inv(Tbasis(H,H.W(-i))),w)
  end
end
Tbasis(H::HeckeAlgebra,w::Vararg{Integer})=Tbasis(H,collect(w))
Tbasis(H::HeckeAlgebra,h::HeckeTElt)=h
Tbasis(H::HeckeAlgebra,h::HeckeElt)=Tbasis(h)
Tbasis(H::HeckeAlgebra,w)=HeckeTElt(MM(w=>one(coefftype(H));check=false),H)

# for each parameter p relation T^length(p)=lower terms
function polynomial_relations(H::HeckeAlgebra{C})where C
  get!(H,:polrel)do
    map(p->Pol([1],length(p))-prod(x->Pol([-x,1]),p),H.para)
  end::Vector{Pol{C}}
end

function innermul(W::CoxeterGroup,a,b)
  sum(a.d) do (ea,pa)
    h=b.d*pa
    for i in reverse(word(W,ea))
      s=W(i)
      up=empty(h.d)
      down=empty(h.d)
      for (e,p)  in h
        if isleftdescent(W,e,i) push!(down,e=>p) else push!(up,s*e=>p) end
      end
      h=MM(up)
      if isempty(down) continue end
      ss=sum(a.H.para[i])
      if !iszero(ss) h+=MM(down;check=false)*ss end
      p=-prod(a.H.para[i])
      if !iszero(p) h+=MM(s*e=>c*p for (e,c) in down) end
    end
    HeckeTElt(h,a.H)
  end
end

function innermul(W::PermRootGroup,a,b)
  if length(refltype(W))>1 || !iscyclic(W)
   error("T basis implemented only for Coxeter groups and G(d,1,1)")
  end
  sum(a.d) do (ea,pa)
    h=b.d*pa
    for i in reverse(word(W,ea))
      new=zero(h)
      for (eb,pb) in h
        lb=length(word(W,eb))
        if 1+lb<length(W) push!(new.d,W(1)*eb=>pb)
        else
          p=polynomial_relations(a.H)[1]
          append!(new.d,W(1)^(i+p.v+1)=>pb*c for (i,c) in pairs(p.c))
        end
      end
      h=new
    end
    HeckeTElt(MM(h.d),a.H)
  end
end
# implement also inverse see cmplximp.g/H.inverse

function Base.:*(a::HeckeTElt, b::HeckeTElt)
  if iszero(a) return a end
  if iszero(b) return b end
  W=a.H.W
  innermul(W,a,b)  # function barrier needed for performance
end

function Base.inv(a::HeckeTElt)
  H=a.H
  W=H.W
  if !(W isa CoxeterGroup) error("only implemented for Coxeter groups") end
  if length(a)!=1 error("can only invert single T(w)") end
  w,coeff=first(a)
  l=word(W,w)
  if isempty(l) return inv(coeff)*one(H) end
  if length(l)==1 
    i=only(l)
    s=sum(H.para[i]);p=inv(coeff)*inv(prod(H.para[i]))
    return p*clone(a,ModuleElt(W()=>s,W(i)=>-one(s)))
  end
  d=div(length(l),2)
  T=Tbasis(H)
  inv(T(W(l[d+1:end]...)))*inv(T(W(l[1:d]...)))
end

Cosets.Frobenius(x::HeckeElt,phi)=
     clone(x,MM(Frobenius(k,phi)=>v for (k,v) in x))

"""
`alt(a::HeckeTElt)`

`a` should be an element of an Iwahori-Hecke algebra `H`. the involution on
`H`   defined  by  `x↦  bar(x)`   on  coefficients  and  `Tₛ↦  uₛ,₀uₛ,₁Tₛ`.
Essentially it corresponds to tensoring with the sign representation.

```julia-repl
julia> W=coxgroup(:G,2);H=hecke(W,Pol(:q))
hecke(G₂,q)

julia> T=Tbasis(H);h=T(1,2)*T(2,1)
q²T.+(q²-q)T₁+(q-1)T₁₂₁

julia> alt(h)
q⁻²T.+(q⁻²-q⁻³)T₁+(q⁻³-q⁻⁴)T₁₂₁
```
"""
function alt(a::HeckeTElt)
  clone(a,MM(isone(w) ? w=>bar(c) : w=>prod(prod(inv.(a.H.para[i]))
                for i in word(a.H.W,w))* bar(c) for (w,c) in a;check=false))
end

"""
`α(a::HeckeTElt)`

the anti-involution on the Hecke algebra defined by ``T_w↦ T_{inv(w)}``.
"""
Garside.α(h::HeckeTElt)=clone(h,MM(inv(p)=>c for (p,c) in h))

"""
`class_polynomials(h::HeckeElt)`

returns  the  class  polynomials  of  the  element `h` of the Iwahori-Hecke
algebra or coset given by `h.H` with respect to the `T` basis for a set `R`
of  representatives  of  minimal  length  in  the  conjugacy classes of the
Coxeter group or coset `H.W`. Such minimal length representatives are given
by  `classreps(H.W)`. The vector `p` of  these polynomials has the property
that  if `X` is the  matrix of the values  of the irreducible characters of
`H`  on `T_w` (for `w∈ R`), then the product `X*p` is the list of values of
the irreducible characters on `h`.

```julia-repl
julia> W=coxsym(4)
𝔖 ₄

julia> H=hecke(W,Pol(:q))
hecke(𝔖 ₄,q)

julia> h=Tbasis(H,longest(W))
T₁₂₁₃₂₁

julia> p=class_polynomials(h)
5-element Vector{Pol{Int64}}:
 0
 0
 q²
 q³-2q²+q
 q³-q²+q-1
```
The class polynomials were introduced in [Geck-Pfeiffer1993](biblio.htm#GP93).
"""
function class_polynomials(h::HeckeElt)
  H=h.H
  WF=H.W
  if H isa HeckeCoset
    W=Group(WF)
    para=H.H.para
  else W=WF
    para=H.para
  end
  minl=length.(word.(conjugacy_classes(WF)))
  h=Tbasis(H,h)
# Since  vF is not of minimal length in its class there exists wF conjugate
# by   cyclic  shift  to  vF  and  a  generating  reflection  s  such  that
# l(swFs)=l(vF)-2. Return T_sws.T_s^2
  function orb(w)
    orbit=[w]
    for w in orbit
      for s in leftdescents(W,w)
        sw=W(s)*w
        sws=sw*W(s)
        if isrightdescent(W,sw,s)
          q1,q2=para[s]
          return (elm=[sws,sw],coeff=[-q1*q2,q1+q2])
        elseif !(sws in orbit) push!(orbit,sws)
        end
      end
    end
    error("Geck-Kim-Pfeiffer theory")
  end

  elm,coeff=first(h.d)
  min=fill(zero(coeff),length(minl))
  while length(h.d)>0
    elms=typeof(elm)[]
    coeffs=typeof(coeff)[]
    l=[length(W,elm) for (elm,coeff) in h.d]
    maxl=maximum(l)
    for (elm,coeff) in h.d
      if length(W,elm)<maxl
        push!(elms,elm)
        push!(coeffs,coeff)
      else
        p=position_class(WF,elm)
        if minl[p]==maxl min[p]+=coeff
        else o=orb(elm)
          append!(elms,o.elm)
          append!(coeffs,o.coeff.*coeff) end
      end
    end
    h=clone(h,MM(Pair.(elms,coeffs)))
  end
  return min
end

"""
`char_values(h::HeckeTElt)`

`h` is an element of an Iwahori-Hecke algebra `H`. The function returns the
values  of the irreducible characters of `H`  on `h` (the method used is to
convert to the `T` basis, and then use `class_polynomials`).

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> H=hecke(W,q^2;rootpara=q)
hecke(B₂,q²,rootpara=q)

julia> char_values(Cpbasis(H)(1,2,1))
5-element Vector{Pol{Int64}}:
 -q-q⁻¹
 q+q⁻¹
 0
 q³+2q+2q⁻¹+q⁻³
 0
```
"""
char_values(h::HeckeElt,ch=CharTable(h.H).irr)=ch*class_polynomials(h)

"""
`char_values(H::HeckeAlgebra,v::Vector{<:Integer})`

For an Iwahori-Hecke algebra this computes the character values of `H` on
the `Tbasis(H)(v)`.

For  `H` the Hecke algebra  of a complex reflection  group `W` this routine
computes  character values on a  lift of the element  of `W` defined by the
word `v` in `gens(W)`.

For  complex reflection  groups the  character table  of the  generic Hecke
algebra  of  `W`  has  been  computed  (not  entirely for 3 exceptions, see
`representation`)  in the  sense that,  if `s₁,…,sₙ`  are generators of the
braid  group lifting  the Broué-Malle-Rouquier-Bessis-Michel  generators of
`W`,  there is at least one element `v`  in each conjugacy class of `W` and
one  expression in the generators for it  such that the character values of
the  image `Tᵥ`  in the  Hecke algebra  of the  lift to the braid group are
known.  Such an expression in the generators  will be called a *known* word
(the  list of known words  is obtained by `word.(conjugacy_classes(W))`. If
the  word `v` is known, the computation is quick using the character table.
If  not,  the  function  computes  the  trace  of  `Tᵥ` in each irreducible
representation.   The   values   returned   are   `Unknown()`  for  missing
representations (see `representation`).
```julia-repl
julia> W=crg(4)
G₄

julia> H=hecke(W,Pol(:q))
hecke(G₄,q)

julia> char_values(H,[2,1,2])
7-element Vector{Pol{Cyc{Int64}}}:
 q³
 1
 1
 0
 0
 0
 -q
```
"""
function char_values(H::HeckeAlgebra,w::Vector{<:Integer})
  W=H.W
  if W isa CoxeterGroup return char_values(Tbasis(H)(w)) end
  p=findfirst(==(w),word.(conjugacy_classes(W)))
  if !isnothing(p) return CharTable(H).irr[:,p] end
  map(representations(H))do r
    if isnothing(r) return nothing end
    first(traces_words_mats(r,[w]))
  end
end

function schur_element(H::HeckeAlgebra,p)
  t=map((t,phi)->getchev(t,:SchurElement,phi,H.para[t.indices],
      haskey(H,:rootpara) ?  H.rootpara[t.indices] :
      fill(nothing,length(H.para))),
      refltype(H.W),p)
  if any(==(false),t) return nothing end
  prod(t)
end

"""
`schur_elements(H)`

returns the list of Schur elements for the Hecke algebra `H`

```julia-repl
julia> H=hecke(complex_reflection_group(4),Pol(:q))
hecke(G₄,q)

julia> s=schur_elements(H)
7-element Vector{Pol{Cyc{Int64}}}:
 q⁸+2q⁷+3q⁶+4q⁵+4q⁴+4q³+3q²+2q+1
 2√-3+(6+4√-3)q⁻¹+12q⁻²+(6-4√-3)q⁻³-2√-3q⁻⁴
 -2√-3+(6-4√-3)q⁻¹+12q⁻²+(6+4√-3)q⁻³+2√-3q⁻⁴
 2+2q⁻¹+4q⁻²+2q⁻³+2q⁻⁴
 ζ₃²√-3q³+(3-√-3)q²+3q+3+√-3-ζ₃√-3q⁻¹
 -ζ₃√-3q³+(3+√-3)q²+3q+3-√-3+ζ₃²√-3q⁻¹
 q²+2q+2+2q⁻¹+q⁻²

julia> CycPol.(s)
7-element Vector{CycPol{Cyc{Int64}}}:
 Φ₂²Φ₃Φ₄Φ₆
 2√-3q⁻⁴Φ₂²Φ′₃Φ′₆
 -2√-3q⁻⁴Φ₂²Φ″₃Φ″₆
 2q⁻⁴Φ₃Φ₄
 ζ₃²√-3q⁻¹Φ₂²Φ′₃Φ″₆
 -ζ₃√-3q⁻¹Φ₂²Φ″₃Φ′₆
 q⁻²Φ₂²Φ₄
```
"""
schur_elements(H::HeckeAlgebra)=
  improve_type(map(p->schur_element(H,p), charinfo(H.W).charparams))

#@test (H=hecke(coxgroup(:I,2,8),[Mvp(:x)^2,Mvp(:y)^2]);transpose(CharTable(H).irr)*inv.(schur_elements(H))==[1,0,0,0,0,0,0])

#----------------------- Factorized Schur elements
"""
A  `FactSchur` representing  a Schur  element of  the form  `M∏ᵩφ(Mᵩ)` (see
[`factorized_schur_element`](@ref))  is  a  `struct`  with a field `factor`
which  holds the  monomial `M`,  and a  field `vcyc`  which holds a list of
`NamedTuples`  describing each  factor `Mᵩ`  in the  product. An element of
`vcyc`  representing a  term `φ(Mᵩ)`  is itself  a `NamedTuple` with fields
`monomial` holding `Mᵩ` (as an `Mvp` with a single term), and a field `pol`
holding a `CycPol` (see `CycPol`) representing `φ`.

A  few operations are implemented for  `FactSchur`, like `*, lcm`. They can
be  evaluated  partially  or  completely  keeping  as  much as possible the
factored form.

```julia-repl
julia> @Mvp x,y; W=crg(4); H=hecke(W,[[1,x,y]])
hecke(G₄,Vector{Mvp{Int64, Int64}}[[1, x, y]])

julia> p=factorized_schur_element(H,[[2,5]])
-x⁻¹y(xy+1)(x-1)Φ₆(xy⁻¹)(y-1)

julia> q=p(;x=E(3)) # partial evaluation
ζ₃²√-3y⁻¹Φ₁Φ₂Φ′₆²(y)

julia> q(;y=2//1)
-9√-3/2
```

In contrast, the next operation expands `p` to an `Mvp`:

```julia-repl
julia> HeckeAlgebras.expand(p)
Mvp{Cyc{Rational{Int64}},Rational{Int64}}: -x³y+x³+x²y²-2x²+x²y⁻¹-xy³+2xy-xy⁻¹+y³-2y²+1+x⁻¹y²-x⁻¹y
```
"""
mutable struct FactSchur
  factor::Mvp{Cyc{Rational{Int}},Rational{Int}}
  vcyc::Vector{NamedTuple{(:pol,:monomial),
         Tuple{CycPol{Int},Mvp{Cyc{Rational{Int}},Rational{Int}}}}}
end

Base.getindex(x::FactSchur,s::Symbol)=s==:factor ? x.factor : x.vcyc
Base.setindex!(x::FactSchur,v,s::Symbol)=
  if s==:factor x.factor=v else x.vcyc=v end

function Base.show(io::IO,x::FactSchur)
 if !(get(io,:TeX,false)||get(io,:limit,false))
   print(io,"FactSchur(",x.factor,",",x.vcyc,")")
   return
 end
 v=map(x.vcyc) do l
    if get(io,:Maple,false) || degree(l.pol)==1
      "("*xrepr(io,l.pol(l.monomial))*")"
    else
      xrepr(io,l.pol)*"("*xrepr(io,l.monomial)*")"
    end
  end
  if get(io,:GAP,false) || get(io,:Maple,false) v=join(v,"*")
  else v=join(v,"")
  end
  c=xrepr(io,x.factor)
  if length(v)>0 && degree(x.factor)==0 c=format_coefficient(c) end
  print(io,c,v)
end

expand(x::FactSchur)=x.factor*prod(v->v.pol(v.monomial),x.vcyc;init=1)

Base.:*(a::FactSchur, b::Number)=FactSchur(a.factor*b,copy(a.vcyc))
Base.:*(b::Number,a::FactSchur)=FactSchur(a.factor*b,copy(a.vcyc))

function Base.:*(a::FactSchur, b::FactSchur)
  vcyc=copy(a.vcyc)
  for t in b.vcyc
    p = findfirst(x->x.monomial == t.monomial,vcyc)
    if p===nothing push!(vcyc, t)
    else vcyc[p]=(pol=vcyc[p].pol*t.pol,monomial=vcyc[p].monomial)
    end
  end
  FactSchur(a.factor*b.factor,filter(x->degree(x.pol)>0,vcyc))
end

Base.://(a::FactSchur, b::Number)=FactSchur(a.factor//b,copy(a.vcyc))

function Base.://(a::Number, b::FactSchur)
  FactSchur(a//b.factor,map(t->(pol=1//t.pol,monomial=t.monomial),a.vcyc))
end

function Base.://(a::FactSchur,b::FactSchur)
  vcyc=copy(a.vcyc)
  for t in b.vcyc
    p=findfirst(x->x.monomial==t.monomial,vcyc)
    if p===nothing push!(vcyc,(pol=1//t.pol,monomial=t.monomial))
    else vcyc[p]=(pol=vcyc[p].pol//t.pol,monomial=vcyc[p].monomial)
    end
  end
  FactSchur(a.factor//b.factor,filter(x->degree(x.pol)>0,vcyc))
end

function (x::FactSchur)(y...;z...)
  simplify(FactSchur(x.factor(y...;z...),
        map(p->(pol=p.pol,monomial=p.monomial(y...;z...)) , x.vcyc)))
end

function simplify(res::FactSchur)
  R=Rational{Int}
  T=Cyc{R}
  evcyc=NamedTuple{(:pol,:monomial,:power),Tuple{CycPol{T},Mvp{T,R},R}}[]
  factor=res.factor
  for (pol,monomial) in res.vcyc
    k=scalar(monomial)
    if k!==nothing
      factor*=pol(k)
      continue
    end
    k=collect(powers(first(monomials(monomial))))
    if k[1]<0
      pol=subs(pol,Pol()^-1)
      k=-k
      monomial=1//monomial
      factor*=pol.coeff*monomial^pol.valuation
      pol=CycPol(1,0,pol.v)
    end
    c=first(coefficients(monomial))
    n=c*conj(c)
    if isinteger(n)
      n=root(n)//c
      if n!=1
        monomial*=n
        pol=subs(pol,Pol([inv(Root1(n))],1))
        factor*=pol.coeff
        if isone(pol.coeff^2) pol*=pol.coeff
        else pol//=pol.coeff
        end
      end
    end
    power=abs(gcd(numerator.(k)))//lcm(denominator.(k))
    monomial=monomial^(1//power)
    push!(evcyc,(pol=pol,monomial=monomial,power=power))
  end
  if isempty(evcyc)
    FactSchur(factor,NamedTuple{(:pol,:monomial),Tuple{CycPol{T},Mvp{T,R}}}[])
  else
    vcyc=map(collectby(x->x.monomial,evcyc))do fil
      D=lcm(map(x->denominator(x.power), fil))
      P=prod(x->subs(x.pol,Pol()^(D*x.power)),fil)
      p=P(Pol())
      p=improve_type(p)
      f=filter(i->p.c[i]!=0,eachindex(p.c)).-1
      f=gcd(gcd(f),D)
      if f>1
        p=Pol(p.c[1:f:length(p.c)],0)
        D//=f
      end
      (pol=CycPol(p), monomial=fil[1].monomial^(1//D))
    end
    FactSchur(factor,vcyc)
  end
end

function Base.lcm(l::FactSchur...)
  v=collectby(x->x.monomial,vcat(map(x->x.vcyc,l)...))
  FactSchur(1,map(x->(pol=lcm(map(y->y.pol,x)),monomial=x[1].monomial),v))
end

# para=vcat along Hplanes(parameters of Hecke algebra)
# r=rec(coeff,vcyc=vector([vector{Int}(length(para))],n0 cyclotomic pol),
#       possibly: root, rootCoeff, rootUnity
# data=rec(order: perm same length as para,:name,
#          possibly: rootPower, :rootUnityPower)
# u is always nothing
function VFactorSchurElement(para,r,data=nothing,u=nothing)
  if data!==nothing para=para[data[:order]] end
  function monomial(v)
    res=prod(map((x,p)->(x*1//1)^p,para,v))
    if length(v)>length(para) res*=rt^v[end] end
    res
  end
  factor=haskey(r,:coeff) ? r[:coeff] : 1
  if haskey(r, :factor) factor*= monomial(r[:factor]) end
  if haskey(r, :root)
    den=lcm(denominator.(r[:root]))
    rt=monomial(r[:root]*den)
    if haskey(r, :rootCoeff) rt*=r[:rootCoeff] end
    rt=root(rt,den)
    if !isnothing(data) rt*=data[:rootPower] end
  elseif haskey(r, :rootUnity)
    rt=r[:rootUnity]^data[:rootUnityPower]
  end
  vcyc=[(pol=CycPol([1,0,p]),monomial=Mvp(monomial(v))) for (v,p) in r[:vcyc]]
  if factor==0 || isempty(vcyc) return factor end
  return simplify(FactSchur(factor,vcyc))
end

"""
`factorized_schur_element(H,phi)`

returns  the factorized `schur_element`  (see `factorized_schur_elements`) of
the  Hecke algebra  `H` for  the irreducible  character of `H` of parameter
`phi` (see `charinfo(W).charparams`)

```julia-repl
julia> W=complex_reflection_group(4)
G₄

julia> @Mvp x,y; H=hecke(W,[[1,x,y]])
hecke(G₄,Vector{Mvp{Int64, Int64}}[[1, x, y]])

julia> factorized_schur_element(H,[[2,5]])
-x⁻¹y(xy+1)(x-1)Φ₆(xy⁻¹)(y-1)
```
"""
function factorized_schur_element(H::HeckeAlgebra,phi)
  t=map(refltype(H.W),phi)do t,psi
     getchev(t,:FactorizedSchurElement,psi,H.para[t.indices],
    haskey(H,:rootpara) ?  H.rootpara[t.indices] : nothing)
  end
  if false in t return false
  else return prod(t)
  end
end

"""
`factorized_schur_elements(H)`

Let  `H` be  a Hecke  algebra for  the complex  reflection group `W`, whose
parameters are all (Laurent) monomials in some variables `x₁,…,xₙ`, and let
K  be the field of definition of `W`. Then Maria Chlouveraki has shown that
the  Schur elements of `H` take the  particular form `M ∏ᵩ φ(Mᵩ)` where `φ`
runs  over  a  list  of  K-cyclotomic  polynomials,  and  `M`  and `Mᵩ` are
(Laurent)  monomials (in possibly some  fractional powers) of the variables
`xᵢ`.  The  function  `factorized_schur_elements`  returns a data structure
(see `HeckeAlgebras.FactSchur`) which shows this factorization.

```julia-repl
julia> W=complex_reflection_group(4)
G₄

julia> @Mvp x,y; H=hecke(W,[[1,x,y]])
hecke(G₄,Vector{Mvp{Int64, Int64}}[[1, x, y]])

julia> factorized_schur_elements(H)
7-element Vector{Chevie.HeckeAlgebras.FactSchur}:
 x⁻⁴y⁻⁴(xy+1)Φ₁Φ₆(x)Φ₁Φ₆(y)
 (x²y⁻¹+1)Φ₁Φ₆(x)Φ₁Φ₆(xy⁻¹)
 -x⁻⁴y⁵Φ₁Φ₆(xy⁻¹)(xy⁻²+1)Φ₁Φ₆(y)
 -x⁻¹y(xy+1)(x-1)Φ₆(xy⁻¹)(y-1)
 -x⁻⁴y(x²y⁻¹+1)(x-1)(xy⁻¹-1)Φ₆(y)
 x⁻¹y⁻¹Φ₆(x)(xy⁻¹-1)(xy⁻²+1)(y-1)
 x⁻²y(x²y⁻¹+1)(xy+1)(xy⁻²+1)
```
"""
factorized_schur_elements(H::HeckeAlgebra)=
    map(p->factorized_schur_element(H,p),charinfo(H.W).charparams)

const FactorizedSchurElements=factorized_schur_elements
#---------------------- Hecke Cosets
@doc """
`HeckeCoset`s  are  `Hϕ`  where  `H`  is  an  Iwahori-Hecke algebra of some
Coxeter  group `W` on which the automorphism `ϕ` of some Spets `Wϕ` acts by
`ϕ(T_w)=T_{ϕ(w)}`.  For Weyl groups, this corresponds  to the action of the
Frobenius  automorphism  on  the  commuting  algebra  of the induced of the
trivial  representation from the  rational points of  some `F`-stable Borel
subgroup to `𝐆 ^F`.

```julia-repl
julia> WF=rootdatum(:u,3)
u₃

julia> HF=hecke(WF,Pol(:v)^2;rootpara=Pol())
hecke(u₃,v²,rootpara=v)

julia> CharTable(HF)
CharTable(hecke(u₃,v²,rootpara=v))
┌───┬──────────┐
│   │ 111 21  3│
├───┼──────────┤
│111│  -1  1 -1│
│21 │-2v³  .  v│
│3  │  v⁶  1 v²│
└───┴──────────┘
```
Thanks  to the work of Xuhua He and Sian Nie, 'class_polynomials' also make
sense for these cosets. This is used to compute such character tables.
""" HeckeCoset
@GapObj struct HeckeCoset{TH<:HeckeAlgebra,TW<:Spets}
  H::TH
  W::TW
end

"`hecke(HF::HeckeCoset)` returns the underlying Hecke algebra"
hecke(HF::HeckeCoset)=HF.H

"""
`hecke(WF::Spets, H)`

`hecke(WF::Spets, params)`

Construct  a `HeckeCoset`  from a  Coxeter coset  `WF` and an Hecke algebra
associated to `Group(WF)`. The second form is equivalent to
`Hecke(WF,Hecke(Group(WF),params))`. See the doc for `HeckeCoset`.
"""
hecke(WF::Spets,H::HeckeAlgebra)=HeckeCoset(H,WF,Dict{Symbol,Any}())
hecke(WF::Spets,a...;b...)=HeckeCoset(hecke(Group(WF),a...;b...),WF,Dict{Symbol,Any}())

function Base.show(io::IO, H::HeckeCoset)
  print(io,"hecke(",H.W,",")
  tr(p)= p[2]==-one(p[2]) ? p[1] : p
  if allequal(H.H.para) print(io,tr(H.H.para[1]))
  else print(io,map(tr,H.H.para))
  end
  if haskey(H.H,:rootpara)
    rp=rootpara(H.H)
    if allequal(rp) print(io,",rootpara=",rp[1])
    else print(io,",rootpara=",rp)
    end
  end
  print(io,")")
end

function Chars.CharTable(H::HeckeCoset;opt...)
  get!(H,:chartable)do
    W=H.W
    cts=map(refltype(W))do t
      inds=t.orbit[1].indices
      ct=getchev(t,:HeckeCharTable,H.H.para[inds], haskey(H.H,:rootpara) ?
               rootpara(H.H)[inds] : fill(nothing,length(H.H.para)))
      if haskey(ct,:irredinfo) names=getindex.(ct[:irredinfo],:charname)
      else                     names=charnames(t;opt...,TeX=true)
      end
      CharTable(improve_type(toM(ct[:irreducibles])),names,
         ct[:classnames],Int.(ct[:centralizers]),length(W),
         Dict{Symbol,Any}(:name=>ct[:identifier]))
    end
    ct=prod(cts)
    ct.name=xrepr(H;TeX=true)
    ct.group=H
    ct
  end::CharTable
end

function Chars.representation(H::HeckeCoset,i::Int)
  dims=getchev(H.W,:NrConjugacyClasses)
  if isempty(dims) return (gens=Matrix{Int}[],F=fill(0,0,0)) end
  tt=refltype(H.W)
  rp=haskey(H.H,:rootpara) ? rootpara(H.H) : fill(nothing,length(H.H.para))
  mm=map(tt,lin2cart(dims,i)) do t,j
    r=getchev(t,:HeckeRepresentation,H.H.para[t.orbit[1].indices],rp,j)
    if r==false return nothing
    elseif r isa Vector
      if !(r[1] isa Matrix) r=toM.(r) end
      r=improve_type(r)
      (gens=r,F=one(r[1]))
    else
      if !(r[:gens][1] isa Matrix) r1=toM.(r[:gens]) else r1=r[:gens] end
      if !(r[:F] isa Matrix) F=toM(r[:F]) else F=r[:F] end
      (gens=improve_type(r1),F=improve_type(F))
    end
  end
  if any(isnothing,mm) return nothing end
  n=length(tt)
  if n==1 return (gens=mm[1].gens,F=mm[1].F) end
  (gens=vcat(map(1:n) do i
     map(mm[i].gens) do m
       kron(map(j->j==i ? m : one(mm[j].gens[1]),1:n)...)
     end
   end...), F=kron(getindex.(mm,:F)...))
end

"""
`representations(H)`

returns  the list  of representations  of the  Hecke algebra or Hecke coset
`H` (see `representation`).

```julia-repl
julia> WF=rootdatum("2B2")
²B₂

julia> H=hecke(WF,Pol(:x)^2;rootpara=Pol())
hecke(²B₂,x²,rootpara=x)

julia> representations(H)
3-element Vector{NamedTuple{(:gens, :F)}}:
 (gens = Matrix{Pol{Int64}}[[x²;;], [x²;;]], F = [1;;])
 (gens = [[-1;;], [-1;;]], F = [1;;])
 (gens = Matrix{Pol{Cyc{Int64}}}[[-1 0; √2x x²], [x² √2x; 0 -1]], F = [0 -1; -1 0])
```
"""
Chars.representations(H::Union{HeckeAlgebra,HeckeCoset})=representation.(Ref(H),1:nconjugacy_classes(H.W))

function isrepresentation(H::HeckeCoset,r)
  all(i->iszero(r.gens[i]*r.F-r.F*r.gens[action(Group(H.W),i,H.W.phi)]),
      eachindex(r.gens)) && isrepresentation(H.H,r.gens)
end

struct HeckeTCElt{TH<:HeckeCoset,C,P}<:HeckeElt{TH,C,P}
  d::MM{P,C}
  H::TH
end

basisname(h::HeckeTCElt)="T"
clone(h::HeckeTCElt,d)=HeckeTCElt(d,h.H)

Tbasis(H::HeckeCoset)=(x...)->isempty(x) ? one(H) : Tbasis(H,x...)
Tbasis(H::HeckeCoset,w::Vararg{Integer})=Tbasis(H,H.W(w...))
Tbasis(H::HeckeCoset,w::Vector{<:Integer})=Tbasis(H,H.W(w...))
Tbasis(H::HeckeCoset,h::HeckeTCElt)=h
Tbasis(H::HeckeCoset,h::HeckeElt)=Tbasis(h)
Tbasis(H::HeckeCoset,w)=HeckeTCElt(MM(w=>one(coefftype(H.H))),H)

end
