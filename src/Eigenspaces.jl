"""
Eigenspaces and `d`-Harish-Chandra series

Let  `Wϕ` be a  reflection coset on  a vector space  `V`; that is `Φ∈GL(V)`
normalizes  the reflection group  `W`. Let `Lwϕ`  be a reflection subcoset;
that  is `L` is a  parabolic subgroup of `W`  (the fixator of a subspace of
`V`)  and  `w∈  W`  is  such  that  `wΦ`  normalizes `L`. There are several
interesting  cases where the *relative group* ``N_W(Lwϕ)/L``, or a subgroup
of it normalizing some further data attached to `L`, is itself a reflection
group.

A first example is the case where `ϕ=1` and `w=1`, `W` is the Weyl group of
a   finite  reductive   group  ``𝐆^F``   and  the   Levi  subgroup  ``𝐋^F``
corresponding  to `L` has a cuspidal unipotent character. Then ``N_W(L)/L``
is  a  Coxeter  group  acting  on  the  space  `X(Z𝐋)⊗ℝ`.  A  combinatorial
characterization of such parabolic subgroups of Coxeter groups is that they
are  normalized by the  longest element of  larger parabolic subgroups (see
[5.7.1 Lusztig1976](biblio.htm#Lus76)).

A  second  example  is  when  `L`  is  trivial  and  `wϕ` is a *`ζ`-regular
element*,  that is  the `ζ`-eigenspace  ``V_ζ`` of  `wϕ` contains  a vector
outside  all the reflecting hyperplanes of `W`. Then ``N_W(Lwϕ)/L=C_W(wϕ)``
is a reflection group in its action on ``V_ζ``.

A similar but more general example is when ``V_ζ`` is the `ζ`-eigenspace of
some  element of  the reflection  coset `Wϕ`,  and is  of maximal dimension
among  such `ζ`-eigenspaces. Then the set of  elements of `Wϕ` which act by
`ζ`  on  ``V_ζ``  is  a  certain  subcoset  `Lwϕ`,  and ``N_W(Lwϕ)/L`` is a
reflection group in its action on ``V_ζ`` (see [2.5
Lehrer-Springer1999](biblio.htm#LS99)).

Finally,  a  still  more  general  example,  but which only occurs for Weyl
groups  or  Spetsial  reflection  groups,  is  when `𝐋` is a `ζ`-split Levi
subgroup  (which means that  the corresponding subcoset  `Lwϕ` is formed of
all  the elements which act by `ζ` on  some subspace `V_ζ` of `V`), and `λ`
is  a  `d`-cuspidal  unipotent  character  of  `𝐋`  (which  means  that the
multiplicity  of `ζ`  as a  root of  the degree  of `λ`  is the same as the
multiplicity  of `ζ` as a root of the generic order of the semi-simple part
of  `𝐆`); then ``N_W(Lwϕ,λ)/L`` is a complex reflection group in its action
on `V_ζ`.

Further,  in the above cases the relative group describes the decomposition
of a Lusztig induction.

When  ``𝐆^F`` is  a finite  reductive group,  and `λ`  a cuspidal unipotent
character  of  the  Levi  subgroup  ``𝐋^F``,  then the ``𝐆^F``-endomorphism
algebra  of  the  Harish-Chandra  induced  representation ``R_𝐋^𝐆(λ)`` is a
Hecke algebra attached to the group ``N_W(L)/L``, thus the dimension of the
characters  of this group describe the multiplicities in the Harish-Chandra
induced.

Similarly, when `𝐋` is a `ζ`-split Levi subgroup, and `λ` is a `d`-cuspidal
unipotent  character of  `𝐋` then  (conjecturally) the ``𝐆^F``-endomorphism
algebra  of the Lusztig induced ``R_𝐋^𝐆(λ)``  is a cyclotomic Hecke algebra
for  to the  group ``N_W(Lwϕ,λ)/L``.  The constituents  of ``R_𝐋^𝐆(λ)`` are
called  a  `ζ`-Harish-Chandra  series.  In  the  case of rational groups or
cosets,  corresponding to finite  reductive groups, the  conjugacy class of
`Lwϕ`  depends  only  on  the  order  `d`  of  `ζ`,  so  one  also talks of
`d`-Harish-Chandra  series. These series correspond to `ℓ`-blocks where `l`
is  a prime divisor of `Φ_d(q)` which  does not divide any other cyclotomic
factor of the order of ``𝐆^F``.

The functions described in this module allow to explore these situations.
"""
module Eigenspaces
export relative_degrees, regular_eigenvalues,
  position_regular_class, eigenspace_projector, relative_root, split_levis

using ..Gapjm
"""
`relative_degrees(WF,ζ::Root1=1)`

Let  `WF` be a reflection group or a  reflection coset and `ζ` be a root of
unity.  Then if ``V_ζ`` is the `ζ`-eigenspace  of some element of `WF`, and
is of maximal dimension among such `ζ`-eigenspaces (and if `WF ` is a coset
`W`  is the group of `WF`) then ``N_W(V_ζ)/C_W(V_ζ)`` is a reflection group
in  its  action  on  ``V_ζ``.  The  function `relative_degrees` returns the
reflection  degrees of this complex reflection group, which are a subset of
those of `W`. These degrees are computed by an invariant-theoretic formula:
if   `(d₁,ε₁),…,(dₙ,εₙ)`  are   the  generalized   degrees  of   `WF`  (see
[`degrees`](@ref)) they are the `dᵢ` such that `ζ^{dᵢ}=εᵢ`.

The  eigenvalue `ζ` can also  be specified by giving  an integer `d` (which
then  specifies  `ζ=E(d)`)  or  a  fraction  `a//b`  which  then  specifies
`ζ=E(b,a)`. If omitted, `ζ` is taken to be `1`.
```julia-repl
julia> W=coxgroup(:E,8)
E₈

julia> relative_degrees(W,4) # the degrees of G₃₂
4-element Vector{Int64}:
  8
 12
 20
 24
```
"""
relative_degrees(W,d::Integer)=relative_degrees(W,E(d))
relative_degrees(W,d::Rational)=relative_degrees(W,Root1(;r=d))
relative_degrees(W,d::Root1)=filter(x->isone(d^x),degrees(W))
relative_degrees(W)=relative_degrees(W,E(1))
relative_degrees(W::Spets,z::Root1)=[d for (d,f) in degrees(W) if Cyc(z)^d==f]

"""
`regular_eigenvalues(W)`

Let `W` be a reflection group or a reflection coset. A root of unity `ζ` is
a *regular eigenvalue* for `W` if some element of `W` has a `ζ`-eigenvector
which   lies   outside   of   the   reflecting  hyperplanes.  The  function
returns the list of regular eigenvalues for `W`.
```julia-repl
julia> regular_eigenvalues(coxgroup(:G,2))
6-element Vector{Root1}:
   1
  -1
  ζ₃
 ζ₃²
  ζ₆
 ζ₆⁵

julia> W=complex_reflection_group(6)
G₆

julia> L=twistings(W,[2])[2]
G₆₍₂₎=G₃‚₁‚₁[ζ₄]Φ′₄

julia> regular_eigenvalues(L)
3-element Vector{Root1}:
    ζ₄
  ζ₁₂⁷
 ζ₁₂¹¹
```
"""
function regular_eigenvalues(W)
  d=degrees(W)
  c=codegrees(W)
  l=filter(x->count(iszero.(d.%x))==count(iszero.(c.%x)),
                  sort(union(divisors.(d)...)))
  sort(vcat(map(x->E.(x,prime_residues(x)),l)...),by=x->(order(x),exponent(x)))
end

function regular_eigenvalues(W::Spets)
  d=degrees(W)
  c=codegrees(W)
  l=union(map(p->divisors(order(Root1(p[2]))*p[1]),d)...)
  res=Root1[]
  for n in l
    p=prime_residues(n)
    p1=filter(i->count(p->E(n,i*p[1])==p[2],d)==count(p->E(n,i*p[1])==p[2],c),p)
    append!(res,E.(n,p1))
  end
  res
end

"""
`position_regular_class(WF,ζ::Root1=1)`

Let  `WF` be a reflection group or a  reflection coset and `ζ` be a root of
unity  such that some element of `WF` has a non-trivial `ζ`-eigenspace. The
function   returns  the   index  of   a  conjugacy   class  of  `WF`  whose
`ζ`-eigenspace  is  maximal  (amongst  all  `ζ `-eigenspaces of elements of
`W`).  If not element of `WF` has a non-trivial `ζ`-eigenspace the function
returns `nothing`.

The  eigenvalue `ζ` can also  be specified by giving  an integer `d` (which
then  specifies  `ζ=E(d)`)  or  a  fraction  `a//b`  which  then  specifies
`ζ=E(b,a)`. If omitted, `ζ` is taken to be `1`.
```julia-repl
julia> W=coxgroup(:E,8)
E₈

julia> position_regular_class(W,30)
65

julia> W=complex_reflection_group(6)
G₆

julia> L=twistings(W,[2])[2]
G₆₍₂₎=G₃‚₁‚₁[ζ₄]Φ′₄

julia> position_regular_class(L,7//12)
2
```
"""
position_regular_class(W,d::Integer)=position_regular_class(W,E(d,1))
position_regular_class(W,d::Rational)=position_regular_class(W,Root1(;r=d))
function position_regular_class(W,d::Root1=E(1))
  drank=length(relative_degrees(W,d))
  if drank==0 return nothing end
  return findfirst(x->count(==(d),x)==drank,refleigen(W))
end

"""
`eigenspace_projector(WF,w,ζ::Root1=1)`

Let `WF` be a reflection group or a reflection coset, let `w` be an element
of  `WF` and let  `ζ` be a  root of unity.  The function returns the unique
`w`-invariant projector on the `ζ`-eigenspace of `w`.

The  eigenvalue `ζ` can also  be specified by giving  an integer `d` (which
then  specifies  `ζ=E(d)`)  or  a  fraction  `a//b`  which  then  specifies
`ζ=E(b,a)`. If omitted, `ζ` is taken to be `1`.
```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> w=W(1:3...)
(1,12,3,2)(4,11,10,5)(6,9,8,7)

julia> p=eigenspace_projector(W,w,1//4)
3×3 Matrix{Cyc{Rational{Int64}}}:
  (1+ζ₄)/4   ζ₄/2  (-1+ζ₄)/4
  (1-ζ₄)/4   1//2   (1+ζ₄)/4
 (-1-ζ₄)/4  -ζ₄/2   (1-ζ₄)/4

julia> GenLinearAlgebra.rank(p)
1
```
"""
eigenspace_projector(W,w,d::Integer=1)=eigenspace_projector(W,w,E(d))
eigenspace_projector(W,w,d::Rational)=eigenspace_projector(W,w,Root1(;r=d))
function eigenspace_projector(W,w,d::Root1)
  c=refleigen(W)[position_class(W,w)]
  c=Cyc.(filter(x->x!=d,c))
  f=reflrep(W,w)
  if length(c)==0 one(f)
  else prod(x->f-one(f)*x,c)//prod(x->Cyc(d)-x,c)
  end
end

"""
`relative_root(W,L,i)`

return 'action' of refls(W,i)∉ L on X(Z_L)
a namedtuple with root, coroot, eigenvalue, index=i and parentMap
"""
function relative_root(W,L,i)
# xprintln("W=",W," i=",i," L=",L)
  N=normalizer(reflection_subgroup(W,vcat(inclusiongens(L,W),[i])),L)
  if length(N)==length(L) error("N(L)==L") end
  F=N/L
# xprintln(abelian_gens(elements(F)))
  if  !iscyclic(F)  error("in theory N/L expected to be cyclic") end
  d=length(F)
# println("d=$d ",order.(elements(F)))
  for rc in filter(x->order(x)==d,elements(F))
    m=central_action(L,reflrep(L,rc.phi))
    r=reflection(m)
    if r===nothing error("not a reflection --- This should not happen") end
#   println("rc=$rc")
    if r.eig==E(d)
      rc=filter(c->central_action(L,reflrep(L,c))==m,classreps(N))
#     println("rc=$rc")
      c=unique!(sort(map(x->position_class(W,x),rc)))
      m=maximum(map(x->count(isone,x),refleigen(W)[c]))
      m=filter(x->count(isone,refleigen(W)[x])==m,c)
      m=filter(x->position_class(W,x) in m,rc)
      if any(x->order(x)==d,m)  m=filter(x->order(x)==d,m) end
      return (root=r.root,coroot=r.coroot,eigenvalue=r.eig,index=i,
              parentMap=m[1])
    end
  end
  error("no root found for reflection along root ",i,"\n")
end

"""
`split_levis(WF,ζ::Root1=1[,ad])`

Let  `WF`  be  a  reflection  group  or  a  reflection  coset.  If `W` is a
reflection group it is treated as the trivial coset 'Spets(W)'.
A  *Levi*  is  a  subcoset  of  the  form `W₁F₁` where `W₁` is a *parabolic
subgroup* of `W`, that is the centralizer of some subspace of `V`, and `F₁∈
WF`.

Let `ζ` be a root of unity. `split_levis` returns a list of representatives
of  conjugacy classes  of `ζ`-split  Levis of  `W`. A  `ζ`-split Levi  is a
subcoset  of `WF` formed  of all the  elements which act  by `ζ` on a given
subspace  `V_ζ`. If the additional argument  `ad` is given, it returns only
those subcosets such that the common `ζ`-eigenspace of their elements is of
dimension  `ad`. These notions make sense  and thus are implemented for any
complex reflection group.

In  terms of algebraic groups, an `F`-stable Levi subgroup of the reductive
group  `𝐆  `  is  `ζ`-split  if  and  only  if it is the centralizer of the
`Φ`-part  of its center, where `Φ` is a cyclotomic polynomial with `ζ` as a
root. When `ζ=1`, we get the notion of a *split* Levi, which is the same as
a Levi sugroup of an `F`-stable parabolic subgroup of `𝐆 `.

The  eigenvalue `ζ` can also  be specified by giving  an integer `d` (which
then  specifies  `ζ=E(d)`)  or  a  fraction  `a//b`  which  then  specifies
`ζ=E(b,a)`. If omitted, `ζ` is taken to be `1`.
```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> split_levis(W,4)
2-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 A₃
 A₃₍₎=Φ₂Φ₄

julia> W=spets(coxgroup(:D,4),Perm(1,2,4))
³D₄

julia> split_levis(W,3)
3-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 ³D₄
 ³D₄₍₁₃₎=A₂Φ₃
 ³D₄₍₎=Φ₃²

julia> W=coxgroup(:E,8)
E₈

julia> split_levis(W,4,2)
3-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 E₈₍₃₂₄₅₎=D₄₍₁₃₂₄₎Φ₄²
 E₈₍₅₇₂₃₎=(A₁A₁)×(A₁A₁)Φ₄²
 E₈₍₃₁₅₆₎=²(A₂A₂)₍₁₄₂₃₎Φ₄²

julia> split_levis(complex_reflection_group(5))
4-element Vector{Spets{PRSG{Cyc{Rational{Int64}}, Int16}}}:
 G₅
 G₅₍₁₎=G₃‚₁‚₁Φ₁
 G₅₍₂₎=G₃‚₁‚₁Φ₁
 G₅₍₎=Φ₁²
```
"""
split_levis(W,d=E(1))=[L for ad in 0:length(relative_degrees(W,d))
                             for L in split_levis(W,d,ad)]
split_levis(W,d::Integer,ad)=split_levis(W,E(d),ad)
split_levis(W,d::Rational,ad)=split_levis(W,Root1(;r=d),ad)
function split_levis(WF,d::Root1,ad)
  if WF isa Spets W=WF.W
  else W=WF; WF=spets(W)
  end
  refs=unique_refls(W)
  eig=refleigen(WF)
  cl=filter(j->count(==(d),eig[j])==ad,1:length(eig))
  res=typeof(subspets(WF,Int[]))[]
  while length(cl)>0
    w=classreps(WF)[cl[1]]
    if rank(W)==0 V=fill(0,0,0)
    else m=reflrep(WF,w)
      V=lnullspace(m-Cyc(d)*one(m))
    end
    I=filter(r->iszero(V*coroots(W,r)),refs)
#   println("I=$I\nphi=",w/WF.phi)
    HF=subspets(WF, I, w/WF.phi)
    if isnothing(HF)
      error("subspets($WF,",I,",class#",cl[1],") failed")
    else
      f=intersect(fusion_conjugacy_classes(HF, WF), cl)
      if isempty(f)
        error("fusion is wrong")
        return nothing
      end
#   println("class=",cl[1]," cl=",cl," f=",f)
      cl=setdiff(cl,f)
      H=HF.W
      P=standard_parabolic(W, H)
      if isnothing(P) P=Perm() end
      J=inclusiongens(H).^P
      if P!=Perm() || J!=inclusiongens(H)
         HF=subspets(WF,J,HF.phi^P/WF.phi)
      end
      push!(res, HF)
    end
  end
  return res
end

"""
`standard_parabolic(W::PermRootGroup, H)`

Let  `H` be a  reflection subgroup of  `W`. Returns an  element `w∈ W` such
that  `H^w` is a standard  parabolic subgroup of `W`  (or nothing if `H` is
not parabolic or not conjugate to a standard parabolic).
"""
function Weyl.standard_parabolic(W::PermRootGroup, H)
  hr=inclusiongens(H,W)
# println("hr=",hr)
  if issubset(hr,eachindex(gens(W))) return Perm() end
  I=combinations(eachindex(gens(W)),length(hr))
  # I=map(x->inclusion(W,x),parabolic_reps(W))
  I=filter(l->length(reflection_subgroup(W,l))==length(H),I)
  try_=function(a)
    H1=H
    w=Perm()
    C=W
    for i in 1:length(a)
      if i>1 C=centralizer(C, gens(W)[a[i-1]]) end
      t=transversal(C,inclusion(H1,i))
      w1=get(t,inclusion(W,a[i]))do
        t=transversal(C,H1(i))
        get(t,W(a[i]),nothing)
      end
      if isnothing(w1) return nothing end
      w=w*w1
      H1=H1^w1
    end
    return w
  end
  for l in I
    for a in arrangements(l, length(l))
      w=try_(a)
      if !isnothing(w) return w end
    end
  end
  InfoChevie("# ****** ",H," is not conjugate to a standard parabolic\n")
  return nothing
end

function Weyl.relative_group(W,J::Vector{<:Integer},indices=false)
# println("relative_group: W=$W J=$J indices=$indices")
  res = Dict{Symbol, Any}(:callarg => joindigits(J))
  if indices!=false res[:callarg]*=","*joindigits(indices) end
# res = CHEVIE[:GetCached](W, "RelativeGroups", res, x->x[:callarg])
# if length(keys(res))>1 return res end
  L=reflection_subgroup(W, J)
  if length(J)==0
    W.MappingFromNormalizer=w->PermX(W,central_action(L, reflrep(L,w)))
    W.relativeIndices=eachindex(gens(W))
    return W
  end
  if indices==false
    indices=filter(x->!(refls(W, x) in L),eachindex(gens(W)))
  end
  if sort(inclusion(L))==sort(inclusion(W))
    res[:relativeIndices]=Int[]
    res[:MappingFromNormalizer]=w->Perm()
    res[:rank]=0
    return PRG(Perm{Int16}[],one(W),Matrix{Int}[],Vector{Int}[],Vector{Int}[],res)
  end
  res[:roots] = []
  res[:simpleCoroots] = []
  res[:relativeIndices] = []
  res[:parentMap] = []
  for R in filter(!isnothing,map(r->relative_root(W,L,r),indices))
    p=findfirst(==(R.root),res[:roots])
    if p===nothing || p!=findfirst(==(R.coroot),res[:simpleCoroots])
      push!(res[:roots], R.root)
      push!(res[:simpleCoroots], R.coroot)
      push!(res[:relativeIndices], R.index)
      push!(res[:parentMap], R.parentMap)
    end
  end
  N=length(normalizer(W,L))//length(L)
  res[:roots]=improve_type(res[:roots])
  res[:simpleCoroots]=improve_type(res[:simpleCoroots])
# println(res[:roots],res[:simpleCoroots])
  R=PRG(res[:roots], res[:simpleCoroots])
  if N==length(R)
    merge!(R.prop,res)
    R.MappingFromNormalizer=w->PermX(R, central_action(L, reflrep(L,w)))
#   InfoChevie("W_",W,"(L_",res[:callarg],")==",R,"<",join(R.prop[:relativeIndices]),">")
    printTeX(rio(),"W_",W,"(L_{",res[:callarg],"})==",R,
                  "<",joindigits(R.relativeIndices),">\n")
    return R
  end
  printTeX(rio(),"# warning: W_",W,"(L_{",res[:callarg],"})",":size ",N,
   " not generated by ",joindigits(indices,"!=")," ==>size ",length(R),"\n")
  print("#           trying all other roots\n")
  indices=filter(r->!(refls(W,r) in gens(W)) && !(refls(W,r) in L), 
                  eachindex(roots(W)))
  for r in indices
    l=relative_root(W, L, r)
    if l.root in res[:roots] continue end
    nr=vcat(res[:roots],[l.root])
    ncr=vcat(res[:simpleCoroots],[l.coroot])
    R=PRG(nr,ncr)
    if N==length(R)
      res[:roots]=nr
      res[:simpleCoroots]=ncr
      push!(res[:relativeIndices], l.index)
      push!(res[:parentMap], l.parentMap)
      R.prop=res
      R.MappingFromNormalizer=w->PermX(R, central_action(L, reflrep(L,w)))
      return R
    end
  end
  error("relgroup not found")
end

end
