"""
# Quasi-Semisimple elements of non-connected reductive groups

We also use Coxeter cosets to represented non-connected reductive groups of
the  form `𝐆 ⋊ σ` where  `𝐆 ` is a connected  reductive group and `σ` is an
algebraic automorphism of `𝐆 `; more specifically to represent the coset `𝐆
.σ`.  We may always choose `σ∈𝐆⋅σ` *quasi-semisimple*, which means that `σ`
preserves  a pair `𝐓 ⊂ 𝐁` of a maximal  torus and a Borel subgroup of `𝐆 `,
and  further *quasi-central*, which means that the Weyl group of ``C_𝐆(σ)``
is ``W^σ``. Then `σ` defines an automorphism `F₀` of the root datum `(X(𝐓),
Φ, Y(𝐓), Φᵛ)`, thus a Coxeter coset. We refer to [ss](@cite) for details.

We  have  extended  the  functions  for  semi-simple  elements to work with
quasi-semisimple  elements `tσ∈  𝐓 ⋅σ`.  Here, as  in [ss](@cite), `σ` is a
quasi-central  automorphism uniquely  defined by  a diagram automorphism of
`(W,S)`, taking `σ` symplectic in type `A₂ₙ`.

Here are some examples:

```julia-repl
julia> WF=rootdatum(:u,6)
u₆
```

We  can  see  `WF`  as  the  coset  `GL₆⋅σ`  where  `σ`  is the composed of
transpose, inverse and the longest element of `W`.

```julia-repl
julia> l=quasi_isolated_reps(WF)
4-element Vector{SemisimpleElement{Root1}}:
 <1,1,1,1,1,1>
 <ζ₄,ζ₄,ζ₄,ζ₄³,ζ₄³,ζ₄³>
 <ζ₄,ζ₄,1,1,ζ₄³,ζ₄³>
 <ζ₄,1,1,1,1,ζ₄³>
```

we  define an element `tσ∈ 𝐓 ⋅σ` to  be quasi-isolated if the Weyl group of
`C_𝐆  (tσ)`  is  not  in  any  proper  parabolic  subgroup of ``W^σ``. This
generalizes  the  definition  for  connected  groups.  The  above shows the
elements  `t`  where  `tσ`  runs  over  representatives  of  quasi-isolated
quasi-semisimple  classes of  `𝐆 ⋅σ`.  The given  representatives have been
chosen `σ`-stable.

```julia-repl
julia> centralizer.(Ref(WF),l)
4-element Vector{ExtendedCox{Perm{Int16}, FiniteCoxeterGroup{Perm{Int16},Rational{Int64}}}}:
 C₃₍₃₂₁₎
 ²A₃₍₃₁₂₎
 (A₁A₁)₍₁₃₎×A₁₍₂₎
 B₂Φ₁
```
in  the above example,  the groups `C_𝐆(tσ)`  are computed and displayed as
extended  Coxeter groups (following the same convention as for centralisers
in connected reductive groups).

We  define an element  `tσ∈ 𝐓⋅σ` to  be isolated if  the Weyl group of `C_𝐆
(tσ)⁰` is not in any proper parabolic subgroup of ``W^σ``. This generalizes
the definition for connected groups.

```julia-repl
julia> isisolated.(Ref(WF),l)
4-element BitVector:
 1
 1
 1
 0
```
"""
module Sscoset

using ..Chevie

"""
`isspecial(WF,c)` where `c` is an orbit of `WF.phi` on `roots(Group(WF))`
return true iff `c` is special in the sense of [ss](@cite]
"""
function isspecial(WF::Spets,c)
  if mod(length(c),2)!=0 return false end
  W=Group(WF)
  roots(W,c[1])+roots(W,c[1+div(length(c),2)]) in roots(W)
end

function Chevie.order(m::Matrix)
  o=1
  while m^o!=one(m) o+=1 end
  o
end

"""
`relativedatum(WF)`

computes `X_σ, X^σ, Y_σ, Y^σ, R(σ)` see [ss; 1.1 to 1.7](@cite)
"""
function relativedatum(WF)
  get!(WF,:Rs)do
    W=Group(WF)
    n=order(WF.F) # matrix of sigma on X
    π=sum(i->WF.F^i,0:n-1)//n
    Xₛ=baseInt(Int.(n*π))//n # basis of X_σ
    Yₛ=baseInt(Int.(n*transpose(π)))//n # basis of Y_σ
    Xˢ=(Xₛ*transpose(Yₛ))^-1*Xₛ # basis of X^σ
    Yˢ=(Yₛ*transpose(Xₛ))^-1*Yₛ # basis of Y^σ
    cc=restriction.(Ref(W),orbits(WF.phi, inclusiongens(W)))
    Phis=map(cc)do c
      res=sum(roots(W,c))
      if isspecial(WF,c) res=2res end
      res
    end
    cPhis=map(c->sum(coroots(W,c))//length(c), cc)
    WF.pi=π
    WF.X_s=Xₛ
    WF.Y_s=Yₛ
    WF.Xs=Xˢ
    WF.Ys=Yˢ
    rootdatum(improve_type(toM(map(x->solutionmat(Xˢ, x), Phis))),
              improve_type(toM(map(x->solutionmat(Yₛ, x), cPhis))))
  end
end

"`Cso(WF)` computes the constants ``C_{σ,𝒪}`` of [ss; page 1305](@cite)"
function Cso(WF)
  get!(WF,:Cso)do
    W=Group(WF)
    res=fill(1,nref(W))
    for o in refltype(WF)
      if order(o.twist)==2 && o.orbit[1].series==:A && mod(rank(o.orbit[1]),2)==0
        for p in o.orbit
         res[p.indices[div(rank(p),2).+[0,1]]]=[-1,-1]
        end
      end
    end
    function C(i)
      for j in 1:semisimplerank(W)
        if W.rootdec[i][j]>0
          p=findfirst(==(roots(W,i)-roots(W,j)),roots(W))
          if !isnothing(p) return res[j]*res[p] end
        end
      end
    end
    for i in semisimplerank(W)+1:nref(W) res[i]=C(i) end
    append!(res,res)
  end
end

"""
`centralizer(WF::Spets,t::SemisimpleElement{Root1})`  

`WF`  should be  a Coxeter  coset representing  an algebraic  coset `𝐆 ⋅σ`,
where  `𝐆 ` is a connected reductive group (represented by `W:=Group(WF)`),
and  `σ`  is  a  quasi-central  automorphism  of  `𝐆 ` defined by `WF`. The
element `t` should be a semisimple element of `𝐆 `. The function returns an
extended reflection group describing ``C_𝐆(tσ)``, with the reflection group
part  representing ``C_𝐆(tσ)⁰``,  and the  diagram automorphism  part being
those induced by ``C_𝐆(tσ)/C_𝐆(tσ)⁰`` on ``C_𝐆(tσ)⁰``.

```julia-repl
julia> WF=rootdatum(:u,6)
u₆

julia> s=ss(Group(WF),[1//4,0,0,0,0,3//4])
SemisimpleElement{Root1}: <ζ₄,1,1,1,1,ζ₄³>

julia> centralizer(WF,s)
B₂Φ₁

julia> centralizer(WF,one(s))
C₃₍₃₂₁₎
```
"""
function Groups.centralizer(WF::Spets,t::SemisimpleElement{Root1})
  W=Group(WF)
  Rs=relativedatum(WF)
  refC=centralizer(Rs,ss(Rs,solutionmat(WF.Y_s,WF.pi*map(x->x.r,t.v))))
  Rs=map(restriction.(Ref(W),orbits(WF.phi,inclusion(W,1:nref(W)))))do c
    res=(c,sum(roots(W,c))//length(c),sum(coroots(W,c)))
    if isspecial(WF,c) res[3]=2res[3] end
    res
  end
  labels=joindigits.(orbits(WF.phi, inclusion(W,1:nref(W))))
  good=map(p->Cyc(prod(i->t^roots(W,i),p[1]))*Cso(WF)[p[1][1]]==1,Rs)
  Rs=Rs[good]
  cRs=map(x->x[3], Rs)
  cRs=map(x->solutionmat(WF.Ys, x), cRs)
  scRs=sum.(Iterators.product(cRs, cRs))
  cRs=filter!(x->!(x in scRs),cRs)
  Rs=map(x->x[2], Rs)
  Rs=map(x->solutionmat(WF.X_s, x), Rs)
  sRs=sum.(Iterators.product(Rs, Rs))
  good=map(x->!(x in sRs), Rs)
  Rs=Rs[good]
  labels=labels[good]
  if length(Rs)>0 C=rootdatum(toM(Rs), toM(cRs))
  else C=torus(size(WF.Xs,1))
  end
# C[:operations][:ReflectionFromName] = function (W, x)
#         return findfirst(==(x),inclusion(W))
#     end
  p=solutionmat(WF.X_s,WF.Xs)
# transfer matrix on X^σ to X_σ
  if isempty(refC.F0s) return ExtendedReflectionGroup(C) end
  ExtendedReflectionGroup(C,map(x->Int.(inv(p)*x*p), refC.F0s))
end

"""
`quasi_isolated_reps(WF::Spets,p=0)`

`WF`  should be  a Coxeter  coset representing  an algebraic  coset `𝐆 ⋅σ`,
where  `𝐆 ` is a connected  reductive group (represented by `W=Group(WF)`),
and  `σ`  is  a  quasi-central  automorphism  of  `𝐆 ` defined by `WF`. The
function returns a list of semisimple elements of `𝐆 ` such that `tσ`, when
`t`  runs over this  list, are representatives  of the conjugacy classes of
quasi-isolated quasisemisimple elements of `𝐆 ⋅σ` (an element `tσ∈ 𝐓 ⋅σ` is
quasi-isolated  if  the  Weyl  group  of  `C_𝐆  (tσ)`  is not in any proper
parabolic  subgroup of `W^σ`). If a second  argument `p` is given, it lists
only those representatives which exist in characteristic `p`.

```julia-repl
julia> WF=rootdatum("2E6sc")
²E₆sc

julia> quasi_isolated_reps(WF)
5-element Vector{SemisimpleElement{Root1}}:
 <1,1,1,1,1,1>
 <1,-1,ζ₄,1,ζ₄,1>
 <1,1,1,-1,1,1>
 <1,ζ₃²,1,ζ₃,1,1>
 <1,ζ₄³,1,-1,1,1>

julia> quasi_isolated_reps(WF,2)
2-element Vector{SemisimpleElement{Root1}}:
 <1,1,1,1,1,1>
 <1,ζ₃²,1,ζ₃,1,1>

julia> quasi_isolated_reps(WF,3)
4-element Vector{SemisimpleElement{Root1}}:
 <1,1,1,1,1,1>
 <1,-1,ζ₄,1,ζ₄,1>
 <1,1,1,-1,1,1>
 <1,ζ₄³,1,-1,1,1>
```
"""
function Semisimple.quasi_isolated_reps(WF::Spets,p=0)
  map(x->ss(Group(WF),transpose(WF.Y_s)*map(x->x.r,x.v)), 
      quasi_isolated_reps(relativedatum(WF), p))
end

"""
`isisolated(WF::Spets,t::SemisimpleElement{Root1})`

`WF`  should be  a Coxeter  coset representing  an algebraic  coset `𝐆 ⋅σ`,
where  `𝐆 ` is a connected  reductive group (represented by `W=Group(WF)`),
and  `σ`  is  a  quasi-central  automorphism  of  `𝐆 ` defined by `WF`. The
element  `t` should be a semisimple element of `𝐆 `. The function returns a
boolean describing whether `tσ` is isolated, that is whether the Weyl group
of `C_𝐆 (tσ)⁰` is not in any proper parabolic subgroup of `W^σ`.

```julia-repl
julia> WF=rootdatum(:u,6)
u₆

julia> l=quasi_isolated_reps(WF)
4-element Vector{SemisimpleElement{Root1}}:
 <1,1,1,1,1,1>
 <ζ₄,ζ₄,ζ₄,ζ₄³,ζ₄³,ζ₄³>
 <ζ₄,ζ₄,1,1,ζ₄³,ζ₄³>
 <ζ₄,1,1,1,1,ζ₄³>

julia> isisolated.(Ref(WF),l)
4-element BitVector:
 1
 1
 1
 0
```
"""
function Semisimple.isisolated(WF::Spets,t::SemisimpleElement{Root1})
  Rs=relativedatum(WF)
  t=ss(Rs,solutionmat(WF.Y_s,WF.pi*map(x->x.r,t.v)))
  isisolated(Rs, t)
end

end
