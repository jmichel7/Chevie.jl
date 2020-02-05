"""
Eigenspaces and `d`-Harish-Chandra series

Let `Wϕ` be a reflection coset on a vector space `V` and `Lwϕ` a reflection
subcoset  where `L` is a  parabolic subgroup (the fixator  of a subspace of
`V`).  There  are  several  interesting  cases  where  the *relative group*
`N_W(Lwϕ)/L`, or a subgroup of it normalizing some further data attached to
`L`, is itself a reflection group.

A first example is the case where `ϕ=1` and `w=1`, `W` is the Weyl group of
a finite reductive group `𝐆^F` and the Levi subgroup `𝐋^F` corresponding to
`L`  has a cuspidal unipotent character. Then `N_W(L)/L` is a Coxeter group
acting  on the  space `X(Z𝐋)⊗ℝ`.  A combinatorial  characterization of such
parabolic  subgroups of Coxeter  groups is that  they are normalized by the
longest element of larger parabolic subgroups (see cite[5.7.1]{Lus76}).

A  second  example  is  when  `L`  is  trivial  and  `wϕ` is a *`ζ`-regular
element*,  that  is  the  `ζ`-eigenspace  `V_ζ`  of  `wϕ` contains a vector
outside all the reflecting hyperplanes of `W`. Then `N_W(Lwϕ)/L=C_W(wϕ)` is
a reflection group in its action on `V_ζ`.

A  similar but more general example is  when `V_ζ` is the `ζ`-eigenspace of
some  element of  the reflection  coset `Wϕ`,  and is  of maximal dimension
among such possible `ζ`-eigenspaces. Then the set of elements of `Wϕ` which
act  by `ζ`  on `V_ζ`  is a  certain subcoset  `Lwϕ`, and `N_W(Lwϕ)/L` is a
reflection group in its action on `V_ζ` (see cite[2.5]{LS99}).

Finally,  a  still  more  general  example,  but which only occurs for Weyl
groups  or  Spetsial  reflection  groups,  is  when `𝐋` is a `ζ`-split Levi
subgroup  (which means that  the corresponding subcoset  `Lwϕ` is formed of
all  the elements which act  by `ζ` on some  subspace `V_ζ` of `V`), and `λ
`d`-cuspidal  unipotent character of `𝐋` (which means that the multiplicity
of  `ζ` as a root of  the degree of `λ` is  the same as the multiplicity of
`ζ`  as a root of  the generic order of  the semi-simple part of `𝐆`); then
`N_W(Lwϕ,λ)/L` is a complex reflection group in its action on `V_ζ`.

Further,  in the above cases the relative group describes the decomposition
of a Lusztig induction.

When  `𝐆^F`  is  a  finite  reductive  group,  and `λ` a cuspidal unipotent
character  of the Levi subgroup  `𝐋^F`, then the `𝐆^F`-endomorphism algebra
of  the Harish-Chandra induced representation `R_𝐋^𝐆(λ)` is a Hecke algebra
attached  to the group `N_W(L)/L`, thus  the dimension of the characters of
this group describe the multiplicities in the Harish-Chandra induced.

Similarly, when `𝐋` is a `ζ`-split Levi subgroup, and `λ` is a `d`-cuspidal
unipotent  character  of  `𝐋`  then  (conjecturally) the `𝐆^F`-endomorphism
algebra of the Lusztig induced `R_𝐋^𝐆(λ)` is a cyclotomic Hecke algebra for
to  the group `N_W(Lwϕ,λ)/L`.  The constituents of  `R_𝐋^𝐆(λ)` are called a
`ζ`-Harish-Chandra  series.  In  the  case  of  rational  groups or cosets,
corresponding  to  finite  reductive  groups,  the conjugacy class of `Lwϕ`
depends only on the order `d` of `ζ`, so one also talks of
`d`-Harish-Chandra  series. These series correspond to `ℓ`-blocks where `l`
is  a prime divisor of `Φ_d(q)` which  does not divide any other cyclotomic
factor of the order of `𝐆^F`.

The functions described in this module allow to explore these situations.
"""
module Eigenspaces
export relative_degrees, regular_eigenvalues, regular_eigenvalues,
  PositionRegularClass, EigenspaceProjector, GetRelativeAction,
  GetRelativeRoot, SplitLevis, RelativeGroup

using Gapjm
"""
`relative_degrees(WF,d=0)`

Let  `WF` be a reflection group or a reflection coset. Here `d` specifies a
root  of unity `ζ`: either `d` is an integer and specifies `ζ=E(d)` or is a
fraction  smaller `a/b` with `0<a<b`  and specifies `ζ=E(b,a)`. If omitted,
`d`   is  taken  to  be  `0`,  specifying  `ζ=1`.  Then  if  `V_ζ`  is  the
`ζ`-eigenspace  of some element of `WF`,  and is of maximal dimension among
such   possible  `ζ`-eigenspaces,  and  `W`  is  the  group  of  `WF`  then
`N_W(V_ζ)/C_W(V_ζ)`  is  a  reflection  group  in  its action on `V_ζ`. The
function  `RelativeDegrees` returns the reflection  degrees of this complex
reflection group, which are a subset of those of `W`.

These  degrees are obtained quickly by invariant-theoretic computations: if
`(d₁,ε₁),…,(dₙ,εₙ)`  are the generalized degrees of  `WF` they are the `dᵢ`
such that `ζ^{dᵢ}=εᵢ`.

```julia-repl
julia> W=coxgroup(:E,8)
E₈

julia> relative_degrees(W,4)
4-element Array{Int64,1}:
  8
 12
 20
 24
```
"""
function relative_degrees(W,d=0//1)
  if d isa Integer && !iszero(d) d=1//d end
  if W isa Spets 
    zeta=E(d//1)
    [d for (d,f) in degrees(W) if zeta^d==f]
  else filter(x->iszero(mod1(d*x)),degrees(W))
  end
end

"""
`regular_eigenvalues(W)`

Let `W` be a reflection group or a reflection coset. A root of unity `ζ` is
a *regular eigenvalue* for `W` if some element of `W` has a `ζ`-eigenvector
which   lies   outside   of   the   reflecting  hyperplanes.  The  function
'RelativeDegree' returns a list describing the regular eigenvalues for `W`.
If  all the primitive  `n`-th roots of  unity are regular eigenvalues, then
`n`  is put on the result list.  Otherwise the fractions `a/n` are added to
the  list for each `a` such that  `E(n,a)` is a primitive `n`-root of unity
and a regular eigenvalue for `W`.

```julia_repl
julia> W=coxgroup(:E,8)
E₈

julia> regular_eigenvalues(W)
13-element Array{Int64,1}:
  1
  2
  4
  8
  3
  6
 12
  5
 10
 20
 24
 15
 30

julia> W=ComplexReflectionGroup(6)
G₆
```
    gap> L:=Twistings(W,[2])[4];
    Z3[I]<2>.(q-I)
    gap> RegularEigenvalues(L);
    [ 1/4, 7/12, 11/12 ]|
"""
function regular_eigenvalues(W)
  d=degrees(W)
  c=codegrees(W)
  if !(W isa Spets)
    return filter(x->
     count(i->mod(i,x)==0,d)==count(i->mod(i,x)==0,c),union(divisors.(d)...))
  end
  l=union(map(p->divisors(denominator(Root1(p[2]).r)*p[1]),d)...)
  res=Rational{Int}[]
  for n in l
    p=prime_residues(n)
    p1=filter(i->count(p->E(n,i*p[1])==p[2],d)==count(p->E(n,i*p[1])==p[2],c),p)
    if p1==p push!(res,n)
    else res=append!(res,p1//n)
    end
  end
  res
end

"""
`PositionRegularClass(WF,d=0)`

Let  `WF` be a reflection group or a reflection coset. Here `d` specifies a
root  of unity `ζ`: either `d` is an integer and specifies `ζ=E(d)' or is a
fraction  smaller `a/b` with `0<a<b`  and specifies `ζ=E(b,a)'. If omitted,
`d`  is taken to be `0`, specifying `ζ=1`. The root `ζ` should be a regular
eigenvalue  for `WF` (see "regular_eigenvalues").  The function returns the
index of the conjugacy class of `WF` which has a `ζ`-regular eigenvector.

```julia-repl
julia> W=coxgroup(:E,8)
E₈

julia> PositionRegularClass(W,30)
65

julia> W=ComplexReflectionGroup(6)
G₆
```
    gap> L:=Twistings(W,[2])[4];
    Z3[I]<2>.(q-I)
    gap> PositionRegularClass(L,7/12);
    2|
"""
function PositionRegularClass(W,d=0)
  if d isa Int && d!=0 d=mod1(1//d) end
  drank=length(relative_degrees(W,d))
  if drank==0 return nothing end
  return findfirst(x->count(==(d),x)==drank,refleigen(W))
end

"""
`EigenspaceProjector(WF,w,d)`

Let  `WF` be a reflection group or a reflection coset. Here `d` specifies a
root  of unity `ζ`: either `d` is an integer and specifies `ζ=E(d)' or is a
fraction  smaller `a/b` with `0<a<b` and specifies `ζ=E(b,a)'. The function
returns the unique <w>-invariant projector on the `ζ`-eigenspace of `w`.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> w=W(1:3...)
(1,12,3,2)(4,11,10,5)(6,9,8,7)

julia> EigenspaceProjector(W,w,1//4)
3×3 Array{Cyc{Rational{Int64}},2}:
  1/4+ζ₄/4   ζ₄/2  -1/4+ζ₄/4
  1/4-ζ₄/4    1/2   1/4+ζ₄/4
 -1/4-ζ₄/4  -ζ₄/2   1/4-ζ₄/4

```
    gap> RankMat(last);
    1|
"""
function EigenspaceProjector(WF, w, d)
  if d isa Int && d!=0 d=mod1(1//d) end
  c=Root1.(refleigen(WF)[position_class(WF,w)])
  c=E.(filter(x->x!=Root1(d),c))
  f=matX(WF,w)
  if length(c)==0 f^0
  else prod(x->f-f^0*x,c)//prod(x->E(Root1(d))-x,c)
  end
end

function GetRelativeAction(W,L,w)
  m=matX(parent(W), w)
  if size(m,2)==0 return m end
  m=m^inv(BaseX(L))
  m[semisimplerank(L)+1:end,semisimplerank(L)+1:end]
end

function GetRelativeRoot(W,L,i)
  J=inclusion(L)[L.generatingReflections]
  N=Normalizer(ReflectionSubgroup(W,Concatenation(J,[i])),L)
  F=N/L
  if  !IsCyclic(F)  Error("in theory N/L expected to be cyclic") end
  d=Size(F)
  for rc in Filtered(Elements(F),x->Order(F,x)==d)  
    m=GetRelativeAction(W,L,rc.element.representative)
    r=AsReflection(m)
    if r==false  Error("I thought this does  !  happen") end
    if r.eigenvalue==E(d)   
      r.index=i
      rc=Filtered(List(ConjugacyClasses(N),Representative),
          c->GetRelativeAction(W,L,c)==m)
      c=Set(List(rc,x->PositionClass(W,x)))
      m=Maximum(List(refleigen(W){c},x->Number(x,y->y==0)))
      m=filter(x->count(==(0),refleigen(W,x))==m,c)
      m=filter(x->PositionClass(W,x) in m,rc)
      if any(x->order(x)==d,m)  m=filter(x->order(x)==d,m) end
      r.parentMap=m[1]
      return r 
    end 
  end
  Error("no root found for reflection",i,"\n")
end

"""
`SplitLevis(WF,d=0,ad=-1)`

Let  `WF`  be  a  reflection  group  or  a  reflection  coset.  If `W` is a
reflection group it is treated as the trivial coset 'Spets(W)'.

Here  `d`  specifies  a  root  of  unity  `ζ`: either `d` is an integer and
specifies  `ζ=E(d)`  or  is  a  fraction  `a/b`  with `0<a<b` and specifies
`ζ=E(b,a)`. If omitted, `d` is taken to be `0`, specifying `ζ=1`.

A  *Levi*  is  a  subcoset  of  the  form `W₁F₁` where `W₁` is a *parabolic
subgroup* of `W`, that is the centralizer of some subspace of `V`.

The  function returns  a list  of representatives  of conjugacy  classes of
`d`-split  Levis of `W`. A  `d`-split Levi is a  subcoset of `WF` formed of
all  the  elements  which  act  by  `ζ`  on  a given subspace `V_ζ`. If the
additional  argument `ad`  is given,  it returns  only those subcosets such
that the common `ζ`-eigenspace of their elements is of dimension `ad`.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> SplitLevis(W,4)
2-element Array{Any,1}:
 A₃   
 .Φ₂Φ₄

julia> W=spets(coxgroup(:D,4),Perm(1,2,4))
³D₄

julia> SplitLevis(W,3)
3-element Array{Any,1}:
 ³D₄ 
 A₂Φ₃
 ³D₄ 

julia> W=coxgroup(:E,8)
E₈
```

|   gap> SplitLevis(W,4,2); # too difficult fusion classes
    [ D4<3,2,4,5>.(q^2+1)^2, (A1xA1)<5,7>x(A1xA1)<2,3>.(q^2+1)^2,
      2(A2xA2)<3,1,5,6>.(q^2+1)^2 ]|
"""
function SplitLevis(WF,d=0,ad=-1)
  if WF isa Spets W=WF.W
  else W=WF; WF=spets(W)
  end
  if d isa Int && d!=0 d=mod1(1 // d) end
  if ad==-1 return vcat(map(ad->SplitLevis(WF, d,
                                   ad),0:length(relative_degrees(WF,d)))...)
  end
  refs=eachindex(reflections(W)[1:nref(W)]) #hum
#  refs=filter(i->Position(reflections(W),reflection(W,i))==i,
#              eachindex(roots(W)))
  mats=map(i->matX(W,reflection(W,i)),refs)
  eig=refleigen(WF)
  cl=filter(j->count(==(d),eig[j])==ad,1:length(eig))
  res=[]
  while length(cl)>0
    w=class_reps(WF)[cl[1]]
    if rank(W)==0 V=fill(0,0,0)
    else m=matX(WF,w)
      V=permutedims(GLinearAlgebra.nullspace(permutedims(m-E(denominator(d),numerator(d))*one(m))))
    end
    I=refs[map(m->V==V*m, mats)]
    HF=subspets(WF, inclusion(W)[I], w/WF.phi)
    if isnothing(HF)
      error("subspets($WF,",inclusion(W)[I],",class#",cl[1],") failed")
    else
      f=intersect(fusion_conjugacy_classes(HF, WF), cl)
      if length(f)==0
          error("fusion is wrong")
          return false
      end
      cl=setdiff(cl,f)
      H=HF.W
      P=standard_parabolic(W, H)
      if P==false P=Perm() end
      P=[P, inclusion(H)[eachindex(gens(H))].^P]
      if P[1]!=Perm() || P[2]!=inclusion(H)[eachindex(gens(H))]
         HF=subspets(WF,P[2],HF.phi^P[1]/WF.phi)
      end
      push!(res, HF)
    end
  end
  return res
end

function PermRoot.standard_parabolic(W::PermRootGroup, H)
  wr=inclusion(W,eachindex(gens(W)))
  hr=inclusion(H,eachindex(gens(H)))
  if issubset(hr,wr) return Perm() end
  I=combinations(wr,length(hr))
  I=filter(l->length(reflection_subgroup(W,l))==length(H),I)
  try_=function(a)
    H1=H
    w=Perm()
    for i in 1:length(a)
      C=centralizer(W, reflection_subgroup(W,a[1:i-1]))
      t=transversal(C,inclusion(H1)[i])
      if haskey(t,a[i]) w1=t[a[i]]
      else
        t=transversal(C,H1(i))
        r=reflection(W,restriction(W)[a[i]])
        if haskey(t,r) w1=t[r] else return nothing end
      end
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
  InfoChevie("\n# ****** ", hr, " not conjugate to a standard parabolic\n")
  return nothing
end

function RelativeGroup(W,J,indices=false)
  res = Dict{Symbol, Any}(:callarg => IntListToString(J))
  if indices!=false
      res[:callarg]=Append(res[:callarg], ",")
      res[:callarg]=Append(res[:callarg], IntListToString(arg[3]))
  end
  res = CHEVIE[:GetCached](W, "RelativeGroups", res, x->x[:callarg])
  if length(keys(res))>1 return res end
  L = ReflectionSubgroup(W, J)
  if length(J)==0
    W[:MappingFromNormalizer]=w->PermMatX(W,GetRelativeAction(W, L, w))
    W[:relativeIndices] = W[:generatingReflections]
    return W
  end
  if indices==false indices=Filtered(inclusion(W)[W[:generatingReflections]],
                         x->!(Reflection(Parent(W), x)) in L)
  end
  if gapSet(L[:rootInclusion]) == gapSet(W[:rootInclusion])
    Inherit(res, PermRootGroup([]))
    res[:relativeIndices] = []
    res[:MappingFromNormalizer] = (w->begin Perm() end)
    return res
  end
  res[:MappingFromNormalizer] = w->PermMatX(res, GetRelativeAction(W, L, w))
  res[:roots] = []
  res[:simpleCoroots] = []
  res[:relativeIndices] = []
  res[:parentMap] = []
  for R = Filtered(map(r->GetRelativeRoot(W, L, r), indices), x->x!=false)
    if Position(res[:roots],R[:root])==false || 
        Position(res[:roots],R[:root])!=Position(res[:simpleCoroots],R[:coroot])
      push!(res[:roots], R[:root])
      push!(res[:simpleCoroots], R[:coroot])
      push!(res[:relativeIndices], R[:index])
      push!(res[:parentMap], R[:parentMap])
    end
  end
  N = Size(Normalizer(W, L)) // Size(L)
  R = PermRootGroup(res[:roots], res[:simpleCoroots])
  relname = SPrint("W_", ReflectionName(W), "(L_", res[:callarg], ")")
  if N == Size(R)
    Inherit(res, R)
    InfoChevie(relname,"==",ReflectionName(res),"<",Join(res[:relativeIndices]),">\n")
    return res
  end
  print("# warning: ", relname, ":size ", N, " not generated by ", 
    IntListToString(indices, "!="), " ==>size ", Size(R), "\n#", 
    "           trying all other roots\n")
  indices = Filtered(W[:rootInclusion], r->!Reflection(Parent(W), r) in
                     W[:generators] && !(Reflection(Parent(W), r)) in L )
  for r = indices
    l = GetRelativeRoot(W, L, r)
    if l != false && !(l[:root]) in res[:roots]
      R=PermRootGroup(vcat(res[:roots],[l[:root]]),vcat(res[:simpleCoroots],[l[:coroot]]))
      if N == Size(R)
          push!(res[:roots], l[:root])
          push!(res[:simpleCoroots], l[:coroot])
          push!(res[:relativeIndices], l[:index])
          push!(res[:parentMap], l[:parentMap])
          Inherit(res, R)
          return res
      end
    end
  end
  error("relgroup  !  found")
end

end
