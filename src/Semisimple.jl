"""
## Semisimple elements

We  construct semi-simple elements in two ways. The first way is for finite
order  elements of `𝐓`, which over an algebraically closed field `K` are in
bijection  with elements  of `Y⊗  ℚ /ℤ`  whose denominator  is prime to the
characteristic of `K`. These are represented as a vector of `Rational`s `r`
such  that `0≤r<1`. The function `ss`  constructs such a semisimple element
from a vector of `Rational`s.

More generally a torus `𝐓` over a field `K` is isomorphic to `(Kˣ)^n` where
`n`  is the dimension  of `𝐓`, so  a vector of  elements of `Kˣ`  is a more
general representation which is produced by the function
`SemisimpleElement`;  in  this  setting  the  result  of  `ss` is naturally
interpreted  as a  `Vector{Root1}`, so  it can  also be obtained by calling
`SemisimpleElement` which such a vector.

```julia-repl
julia> G=rootdatum(:sl,4)
sl₄

julia> ss(G,[1//3,1//4,3//4,2//3])
SemisimpleElement{Root1}: <ζ₃,ζ₄,ζ₄³,ζ₃²>

julia> SemisimpleElement(G,[E(3),E(4),E(4,3),E(3,2)])
SemisimpleElement{Root1}: <ζ₃,ζ₄,ζ₄³,ζ₃²>

julia> L=reflection_subgroup(G,[1,3])
A₃₍₁₃₎=A₁×A₁Φ₁

julia> C=algebraic_center(L)
(Z0 = SubTorus(A₃₍₁₃₎=A₁×A₁Φ₁,[1 2 1]), AZ = Group(SemisimpleElement{Root1}[<1,1,-1>]), descAZ = [[1, 2]], ZD = Group(SemisimpleElement{Root1}[<-1,1,1>, <1,1,-1>]))

julia> T=torsion_subgroup(C.Z0,3)
Group(SemisimpleElement{Root1}[<ζ₃,ζ₃²,ζ₃>])

julia> e=sort(elements(T))
3-element Vector{SemisimpleElement{Root1}}:
 <1,1,1>
 <ζ₃,ζ₃²,ζ₃>
 <ζ₃²,ζ₃,ζ₃²>
```
In  the above, the Levi subgroup  `L` of `SL₄` consisting of block-diagonal
matrices  of shape  `2×2` is  constructed. The  function `algebraic_center`
returns  a named tuple with : the  neutral component `Z⁰` of the center `Z`
of `L`, represented by a basis of `Y(Z⁰)`, a complement subtorus `S` of `𝐓`
to  `Z⁰`  represented  similarly  by  a  basis  of  `Y(S)`, and semi-simple
elements  representing the classes of `Z` modulo  `Z⁰` , chosen in `S`. The
classes  `Z/Z⁰` also biject to the fundamental  group as given by the field
`.descAZ`,  see [`algebraic_center`](@ref) for  an explanation. Finally the
semi-simple elements of order 3 in `Z⁰` are computed.

```julia-repl
julia> e[3]^G(2)
SemisimpleElement{Root1}: <ζ₃²,1,ζ₃²>

julia> orbit(G,e[3])
6-element Vector{SemisimpleElement{Root1}}:
 <ζ₃²,ζ₃,ζ₃²>
 <ζ₃²,1,ζ₃²>
 <ζ₃,1,ζ₃²>
 <ζ₃²,1,ζ₃>
 <ζ₃,1,ζ₃>
 <ζ₃,ζ₃²,ζ₃>
```

Here  is the same  computation as above  performed with semisimple elements
whose  coefficients are in the  finite field `FF(4)`, representing elements
of `sl₄(𝔽₄)`.
```julia-repl
julia> G=rootdatum(:sl,4)
sl₄

julia> s=SemisimpleElement(G,Z(4).^[1,2,1])
SemisimpleElement{FFE{2}}: <Z₄,Z₄²,Z₄>

julia> s^G(2)
SemisimpleElement{FFE{2}}: <Z₄,1,Z₄>

julia> orbit(G,s)
6-element Vector{SemisimpleElement{FFE{2}}}:
 <Z₄,Z₄²,Z₄>
 <Z₄,1,Z₄>
 <Z₄²,1,Z₄>
 <Z₄,1,Z₄²>
 <Z₄²,1,Z₄²>
 <Z₄²,Z₄,Z₄²>
```
We can compute the centralizer ``C_𝐆 (s)`` of a semisimple element in `𝐆 `:
```julia-repl
julia> G=coxgroup(:A,3)
A₃

julia> s=ss(G,[0,1//2,0])
SemisimpleElement{Root1}: <1,-1,1>

julia> centralizer(G,s)
A₃₍₁₃₎=(A₁A₁)Φ₂
```
The  result is an  `ExtendedReflectionGroup`; the reflection  group part is
the Weyl group of ``C_𝐆 ⁰(s)`` and the extended part are representatives of
``C_𝐆  (s)``  modulo  ``C_𝐆⁰(s)``  taken  as  diagram  automorphisms of the
reflection  part.  Here  it  is  printed  as  a  coset  ``C_𝐆 ⁰(s)ϕ`` which
generates ``C_𝐆 (s)``.
"""
module Semisimple
using ..Chevie
export algebraic_center, isisolated, SemisimpleElement, ss, torsion_subgroup,
     quasi_isolated_reps, structure_rational_points_connected_centre, 
     semisimple_centralizer_representatives, sscentralizer_reps, affine
export ExtendedCox, ExtendedReflectionGroup
#----------------- Extended Coxeter groups-------------------------------
struct ExtendedCox{T,TW<:FiniteCoxeterGroup{T}}<:Group{T}
  group::TW
  F0s::Vector{Matrix{Int}}
  phis::Vector{T}
end

Groups.gens(W::ExtendedCox)=vcat(gens(W.group),W.phis)

function ExtendedCox(W::FiniteCoxeterGroup{T},F0s::Vector{<:AbstractMatrix})where T
  if isempty(F0s) return ExtendedCox(W,[reflrep(W,W())],[one(W.G)]) end
  ExtendedCox(W,F0s,isempty(F0s) ? T[] : map(F->PermX(W.G,F),F0s))
end

function Base.:*(a::ExtendedCox,b::ExtendedCox)
  id(r)=Matrix{Int}(I,r,r)
  ExtendedCox(a.group*b.group,improve_type(vcat(
                   map(m->cat(m,id(rank(b.group)),dims=(1,2)),a.F0s),
                   map(m->cat(id(rank(a.group)),m,dims=(1,2)),b.F0s))))
end

function Base.show(io::IO,W::ExtendedCox)
  if !get(io,:limit,false) && !get(io,:TeX,false)
     print(io,"ExtendedCox(",W.group,",",W.F0s,",",W.phis,")")
     return
  end
  if isempty(W.phis) print(io,"Extended(",W.group,")")
  elseif length(W.phis)==1 print(io,spets(W.group,W.phis[1]))
  elseif all(x->isone(x^2),W.phis) && length(Group(W.phis))==6
    print(io,W.group,fromTeX(io,"\\rtimes\\mathfrak S_3"))
  else
    ff=map(x->restricted(x,inclusiongens(W.group)),W.phis)
    if all(!isone,ff) || rank(W.group)==0
         print(io,"Extended(",W.group,",");join(io,ff,",");print(io,")")
    else print(io,"<");join(io,spets.(Ref(W.group),W.phis),",");print(io,">")
    end
  end
end

ExtendedReflectionGroup(W,mats::AbstractVector{<:AbstractMatrix{<:Integer}})=ExtendedCox(W,mats)
ExtendedReflectionGroup(W,mats::AbstractMatrix{<:Integer})=ExtendedCox(W,[mats])
ExtendedReflectionGroup(W)=ExtendedReflectionGroup(W,AbstractMatrix{<:Integer}[])

ExtendedReflectionGroup(W,p::Vector{<:Perm})=ExtendedCox(W,
       isempty(p) ? Matrix{Int}[] : reflrep.(Ref(W),p))
ExtendedReflectionGroup(W,p::Perm)=ExtendedCox(W,[reflrep(W,p)])

function ExtendedReflectionGroup(W,mats::Vector{Any})
  if isempty(mats) ExtendedCox(W,empty([fill(0,0,0)]))
  else error("not empty")
  end
end

#----------------------------------------------------------------------------

struct SemisimpleElement{T}
  W::FiniteCoxeterGroup
  v::Vector{T}
end

Base.:*(a::SemisimpleElement,b::SemisimpleElement)=SemisimpleElement(a.W,
                                                                a.v .* b.v)

Base.inv(a::SemisimpleElement)=SemisimpleElement(a.W,inv.(a.v))
Base.:/(a::SemisimpleElement,b::SemisimpleElement)=a*inv(b)
Base.one(a::SemisimpleElement)=SemisimpleElement(a.W,one.(a.v))
Base.isone(a::SemisimpleElement)=all(isone,a.v)
Base.cmp(a::SemisimpleElement,b::SemisimpleElement)=cmp(a.v,b.v)
Base.isless(a::SemisimpleElement,b::SemisimpleElement)=cmp(a,b)==-1

ss(W::FiniteCoxeterGroup,v::AbstractVector{<:Number})=
                      SemisimpleElement(W,map(x->Root1(;r=x),v))

ss(W::FiniteCoxeterGroup)=SemisimpleElement(W,fill(E(1),rank(W)))

Base.:^(a::SemisimpleElement,n::Integer)=SemisimpleElement(a.W,a.v .^n)

Base.:^(a::SemisimpleElement,m::AbstractMatrix)=SemisimpleElement(a.W,
                                 map(v->prod(a.v.^v),eachcol(m)))

Base.:^(a::SemisimpleElement,p::Perm)=a^YMatrix(parent(a.W.G),inv(p))

# scalar product with a root
Base.:^(a::SemisimpleElement,alpha::Vector{<:Number})=prod(a.v.^Int.(alpha))

function Base.show(io::IO, ::MIME"text/plain", r::SemisimpleElement)
  if !haskey(io,:typeinfo) print(io,typeof(r),": ") end
  show(io,r)
end

function Base.show(io::IO,a::SemisimpleElement)
  if hasdecor(io)
    print(io,"<")
    join(io,a.v,",")
    print(io,">")
  else
    print(io,"SemisimpleElement(",a.W,",[")
    join(io,a.v,",")
    print(io,"])")
  end
end

# hash is needed for using SemisimpleElement in Sets/Dicts
Base.hash(a::SemisimpleElement, h::UInt)=hash(a.v, h)
Base.:(==)(a::SemisimpleElement, b::SemisimpleElement)=a.v==b.v

Chevie.order(a::SemisimpleElement{Root1})=lcm(order.(a.v))

function Base.in(s::SemisimpleElement{Root1},T::SubTorus)
  x=solutionmat(vcat(T.gens,T.complement),map(x->x.r,s.v))
  all(isinteger,x[size(T.gens,1)+1:end])
end

"""
`torsion_subgroup(S::SubTorus,n)`

This  function  returns  the  subgroup  of  semi-simple  elements  of order
dividing `n` in the subtorus `S`.

```julia-repl
julia> G=rootdatum(:sl,4)
sl₄

julia> L=reflection_subgroup(G,[1,3])
A₃₍₁₃₎=A₁×A₁Φ₁

julia> C=algebraic_center(L)
(Z0 = SubTorus(A₃₍₁₃₎=A₁×A₁Φ₁,[1 2 1]), AZ = Group(SemisimpleElement{Root1}[<1,1,-1>]), descAZ = [[1, 2]], ZD = Group(SemisimpleElement{Root1}[<-1,1,1>, <1,1,-1>]))

julia> T=torsion_subgroup(C.Z0,3)
Group(SemisimpleElement{Root1}[<ζ₃,ζ₃²,ζ₃>])

julia> sort(elements(T))
3-element Vector{SemisimpleElement{Root1}}:
 <1,1,1>
 <ζ₃,ζ₃²,ζ₃>
 <ζ₃²,ζ₃,ζ₃²>
```
"""
torsion_subgroup(T::SubTorus,n)=Group(map(x->ss(T.group,x//n),eachrow(T.gens)))

# returns (Tso,s-stable representatives of T/Tso) for automorphism s of T
# here m is the matrix of s on Y(T)
# use ss 1.2(1): Ker(1+m+m^2+...)/Im(m-Id)
function FixedPoints(T::SubTorus,m)
# @show T,m
  n=map(z->solutionmat(T.gens,transpose(m)*z),eachrow(T.gens)) #action on subtorus
  if nothing in n || !all(isinteger,toM(n))
    error(m," does not stabilize ",T)
  end
  n=Int.(toM(n))
  fix=lnullspaceInt(n-n^0) # pure sublattice Y(Tso)
  o=order(n)
  Y1=lnullspaceInt(sum(i->n^i,0:o-1)) # pure sublattice Y(T1) where
  # T=T1.Tso almost direct product, thus spaces Y.(1-s) and Y1.(1-s) coincide
  n=baseInt(n-n^0) # basis of Im(1-s)
  m=map(v->solutionmat(n,v),eachrow(Y1)) # basis of Im[(1-s)^{-1} restricted to Y1]
  # generates elements y of Y1⊗ ℚ such that (1-s)y\in Y1
  (SubTorus(T.group,fix*T.gens),
   abelian_gens(map(v->ss(T.group,transpose(Y1*T.gens)*v),m)))
end

"""
`algebraic_center(W)`

`W`  should  be  a  Weyl  group,  or  an extended Weyl group. This function
returns  a description  of the  center `Z` of  the algebraic  group `𝐆 `
defined by `W` as a named tuple with the following fields:

`Z0`: the neutral component `Z⁰` of `Z` as a subtorus of `𝐓`.

`AZ`: representatives in `Z` of `A(Z):=Z/Z⁰` given as a group of semisimple
elements.

`ZD`:  center of the derived subgroup of `𝐆` given as a group of semisimple
elements.

`descAZ`:  if `W`  is not  an extended  Weyl group,  describes `A(Z)`  as a
quotient  of the center  `pi` of the  simply connected covering  of `𝐆` (an
incarnation of the fundamental group). It contains a list of elements given
as  words  in  the  generators  of  `pi`  which  generate the kernel of the
quotient map.
```julia_repl
julia> G=rootdatum(:sl,4)
sl₄

julia> L=reflection_subgroup(G,[1,3])
A₃₍₁₃₎=A₁×A₁

ulia> C=algebraic_center(L)
(Z0 = SubTorus(A₃₍₁₃₎=A₁×A₁Φ₁,[1 2 1]), AZ = Group(SemisimpleElement{Root1}[<1,1,-1>]), descAZ = [[1, 2]], ZD = Group(SemisimpleElement{Root1}[<-1,1,1>, <1,1,-1>]))

julia> G=coxgroup(:A,3)
A₃

julia> s=ss(G,[0,1//2,0])
SemisimpleElement{Root1}: <1,-1,1>

julia> C=centralizer(G,s)
A₃₍₁₃₎=(A₁A₁)Φ₂

julia> algebraic_center(C)
(Z0 = SubTorus(A₃₍₁₃₎=A₁×A₁Φ₁,Matrix{Int64}(undef, 0, 3)), AZ = Group(SemisimpleElement{Root1}[<1,-1,1>]))
```
"""
function algebraic_center(W)
#   [implemented only for connected groups 18/1/2010]
#   [I added something hopefully correct in general. JM 22/3/2010]
#   [introduced subtori JM 2017 and corrected AZ computation]
  extended=W isa ExtendedCox
  if extended
    F0s=W.F0s
    W=W.group
  end
  if istorus(W) Z0=reflrep(W,one(W))
  else 
    Z0=lnullspaceInt(transpose(simpleroots(W)))
  end
  Z0=SubTorus(W,Z0)
  if isempty(Z0.complement) AZ=Vector{Rational{Int}}[]
  else
    m=Z0.complement
    AZ=toL(inv(Rational.(m*permutedims(simpleroots(W))))*m)
  end
  AZ=ss.(Ref(W),AZ)
  if extended # compute fixed space of F0s in Y(T)
    for m in F0s
      AZ=filter(s->s/ss(W,m*map(x->x.r,s.v)) in Z0,AZ)
      if rank(Z0)>0 Z0=FixedPoints(Z0,transpose(m))
        append!(AZ,Z0[2])
        Z0=Z0[1]
      end
    end
  end
  AZ=Group(abelian_gens(AZ),ss(W))
  if extended && length(F0s)>0 return (;Z0,AZ) end
  descAZ=ss.(Ref(W),weightinfo(W)[:CenterSimplyConnected])
  if isempty(descAZ) return (;Z0,AZ,descAZ) end
  descAZ=Group(descAZ)
  ZD=Group(map(s->ss(W,permutedims(simplecoroots(W))*map(x->x.r,s.v)),gens(descAZ)),ss(W))
  toAZ=function(s)
    s=vec(transpose(map(x->x.r,s.v))*simplecoroots(W))
    s=transpose(s)*inv(Rational.(vcat(Z0.complement,Z0.gens)))
    ss(W,vec(transpose(vec(s)[1:semisimplerank(W)])*Z0.complement))
  end
  ssl=toAZ.(gens(descAZ))
  #println("AZ=$descAZ")
  #println("res=",res)
  #println("gens(AZ)=",gens(descAZ))
  #println("ss=$ss")
  descAZ=if isempty(gens(AZ)) map(x->[x],eachindex(gens(descAZ)))
         elseif gens(descAZ)==ssl Vector{Int}[]
         else # map of root data Y(Wsc)->Y(W)
           h=Hom(descAZ,AZ,ssl)
#          println("h=$h")
           map(x->word(descAZ,x),gens(kernel(h)))
         end
  (;Z0,AZ,descAZ,ZD)
end


#------------------------- Affine Weyl groups ----------------------------
@GapObj struct Affine{T,TW}<:CoxeterGroup{Matrix{T}}
  W::TW
  G::CoxGroups.MatCox{T}
end

"""
A  *generalized Cartan matrix* `C`  is a square integer  matrix of size `n`
such  that `cᵢᵢ=2`, `cᵢⱼ≤0` if `i≠j`, and `cᵢⱼ==0` if and only if `cⱼᵢ==0`.
We  say  that  `C`  is  *indecomposable*  if  it  does  not admit any block
decomposition.

Let  `C` be a generalized  Cartan matrix. For `I`  a subset of `{1,…,n}` we
denote  by `C_I` the square  submatrix with indices `i,j`  taken in `I`. If
`v`  is a real vector of length `n`, we write `v>0` if for all `i∈ {1,…,n}`
we  have `vᵢ>0`. It can be shown that `C` is a Cartan matrix if and only if
for  all sets  `I`, we  have `det  C_I>0`; or  equivalently, if and only if
there  exists  `v>0`  such  that  `C.v>0`.  `C` is called an *affine Cartan
matrix*  if for all proper subsets `I` we have `det C_I>0`, but `det C==0`;
or equivalently if there exists `v>0` such that `C.v==0`.

Given  an  irreducible  Weyl  group  `W`  with  Cartan  matrix  `C`, we can
construct  a generalized  Cartan matrix  `C̃` as  follows. Let  `α₀` be the
opposed of the highest root. Then the matrix
``\\left(\\begin{array}{cc}C&C.α₀\\\\  α₀.C&2\\end{array}\\right)``
is  an  affine  Cartan  matrix.  The  affine  Cartan  matrices which can be
obtained  in this way  are those we  are interested in,  which give rise to
affine Weyl groups.

Let `d=n-rank(C)`. A *realization* of a generalized Cartan matrix is a pair
`V,Vᵛ`  of vector spaces of dimension `n+d` together with vectors `α₁,…,αₙ∈
V`  (the *simple roots*), `αᵛ₁,…,αᵛₙ∈ Vᵛ` (the *simple coroots*), such that
`(αᵛᵢ,  αⱼ)=c_{i,j}`.  Up  to  isomorphism,  a  realization  is obtained as
follows: write
``C=\\left(\\begin{array}{c}C_1\\\\C_2\\end{array}\\right)``
where  `C₁` is  of same  rank as  `C`. Then  take `αᵢ`  to be the first `n`
vectors  in a basis of `V`, and take `αᵛⱼ` to be given in the dual basis by
the rows of the matrix
``\\left(\\begin{array}{cc}C₁&0\\\\ C_2&\\hbox{Id}_d\\\\ \\end{array}\\right).``
To  `C` we associate a reflection group  in the space `V`, generated by the
*fundamental  reflections*  `rᵢ`  given  by  `rᵢ(v)=v-(αᵛᵢ,v)αᵢ`. This is a
Coxeter  group, called the *affine Weyl group* `\tilde W` associated to `W`
when we start with the affine Cartan matrix associated to a Weyl group `W`.

The  affine Weyl  group is  infinite; it  has one additional generator `s₀`
(the  reflection with respect to `α₀`) compared  to `W`. We can not use `0`
as  a label  by default  for a  generator of  a Coxeter  group (because the
default  labels are used as indices, and indices start at 1 in Julia) so we
label it as `n+1` where `n` is the numbers of generators of `W`.

```julia-repl
julia> W=affine(coxgroup(:A,4))
Ã₄

julia> diagram(W)
       ————5————
      /         \\
Ã₄   1———2———3———4
```
"""
function affine(W)
  t=refltype(W)
  if length(t)!=1 || !(t[1].series in Symbol.('A':'G'))
    error("affine needs an irreducible Weyl group")
  else
    t=copy.(t)
    t[1].series=Symbol(string(t[1].series),Char(0x00303))
    t[1].indices=vcat(t[1].indices,[length(t[1].indices)+1])
  end
  ex=vcat(1:semisimplerank(W),2*nref(W))
  C=improve_type([cartan(W.G,i,j) for i in ex, j in ex])
  res=Affine(W,coxgroup(C),Dict{Symbol,Any}())
  res.refltype=t
  res
end

# to make Hecke elements of an affine group work need the following
Base.isless(A::Matrix,B::Matrix)=isless(vec(A),vec(B))

Base.show(io::IO,W::Affine)=print(io,refltype(W))

PermRoot.refltype(W::Affine)=W.refltype

@forward Affine.G Base.length, Base.one,
 Groups.gens, Groups.ngens, #PermGroups.reduced,
 Groups.word, Garside.BraidMonoid,
 KL.KLPol, FinitePosets.Poset, CoxGroups.isleftdescent,
 CoxGroups.bruhatless, CoxGroups.coxmat,
 CoxGroups.leftdescents, PermRoot.semisimplerank

## Given an affine Weyl group W, x a vector in the basis of
## simple roots of W.W and w in W, returns the image of x under w.
function AffineRootAction(W,w,x)
  y=permutedims(vcat(x,[0,1]))*w
  y[eachindex(x)]-y[length(x)+1].*roots(W.W)[nref(W.W)]
end

Base.isfinite(W::Affine)=false

function Perms.reflength(W::Affine,w)
  W0=W.W
  r=W0.semisimpleRank
  Id=vcat(Matrix(1I,r,r),fill(0,r)')
  mov=map(v->AffineRootAction(W,w,v)-v,eachrow(Id))
  l=push!(map(i->refls(W0,i),eachindex(gens(W0))),refls(W0,nref(W0)))
  p=reflength(W0,prod(l[word(W,w)]))
  dimw=minimum(length.(filter(x->RankMat(vcat(mov,roots(W0,x)))==length(x),
                      ParabolicSubgroups(W0))))
  2*dimw-p
end

#-----------------------------------------------------
# returns w,s1 such that s1 is in the fundamental alcove of the affine
# Weyl group and s1^w==s
function to_alcove(s::SemisimpleElement{Root1})
  W=s.W
  w=one(W)
  i=0
  l=map(x->x.r,s.v)
  while i<=ngens(W)
    if i==0 
      for j in weightinfo(W)[:highestroot]
        if dot(l,roots(W,j))>1
          l-=(dot(l,roots(W,j))-1)*coroots(W,j)
          w*=refls(W,j)
          continue
        end
      end
    elseif dot(l,roots(W,i))<0
      l-=dot(l,roots(W,i))*coroots(W,i)
      w*=W(i)
      i=0
      continue
    end
    i+=1
  end
  w,ss(W,l)
end

"""
`centralizer(W,s::SemisimpleElement)`

`W`  should  be  a  Weyl  group  or  an extended reflection group and `s` a
semisimple  element of the  algebraic group `G`  corresponding to `W`. This
function  returns the  Weyl group  of ``C_G(s)``,  which describes  it. The
stabilizer  is an extended reflection group, with the reflection group part
equal to the Weyl group of ``C_{G⁰}(s)``, and the diagram automorphism part
being those induced by ``C_G(s)``.

```julia-repl
julia> G=coxgroup(:A,3)
A₃
julia> s=ss(G,[0,1//2,0])
SemisimpleElement{Root1}: <1,-1,1>
julia> centralizer(G,s)
A₃₍₁₃₎=(A₁A₁)Φ₂
```
"""
function Groups.centralizer(W::Group,s::SemisimpleElement)
  # for the computation of A_G(s) see Bonnafé, "Quasi-isolated elements in
  # reductive groups", comm. in algebra 33 (2005) proposition 3.14
  if W isa ExtendedCox
    totalW=Group(vcat(gens(fundamental_group(W.group;full=true)),W.phis))
    W=W.group
  else totalW=fundamental_group(W;full=true)
  end
  p=filter(i->isone(s^roots(W,i)),1:nref(W))
  W0s=reflection_subgroup(W,p)
  w,s0=to_alcove(s)
  l=map(x->reduced(W0s,x),filter(w->s0==s0^w,elements(totalW)).^inv(w))
  N=Group(abelian_gens(l))
  if rank(W)!=semisimplerank(W)
    N=Group(reflrep.(Ref(W),isempty(gens(N)) ? [one(W)] : gens(N)))
  end
  ExtendedReflectionGroup(W0s,gens(N))
end

"""
`quasi_isolated_reps(W,p=0)`

`W`  should be a Weyl  group corresponding to an  algebraic group 𝐆 over an
algebraically  closed field  of characteristic  0. This  function returns a
list  of  semisimple  elements  for  𝐆,  which  are  representatives of the
𝐆-orbits  of quasi-isolated  semisimple elements.  It follows the algorithm
given  in  [bon05](@cite).  If  a  second  argument  `p` is given, it gives
representatives   of   those   quasi-isolated   elements   which  exist  in
characteristic `p`.

```julia-repl
julia> W=coxgroup(:E,6);l=quasi_isolated_reps(W)
5-element Vector{SemisimpleElement{Root1}}:
 <1,1,1,1,1,1>
 <1,-1,1,1,1,1>
 <1,1,1,ζ₃,1,1>
 <ζ₃,1,1,1,1,ζ₃>
 <1,ζ₆,ζ₆,1,ζ₆,1>

julia> map(s->isisolated(W,s),l)
5-element Vector{Bool}:
 1
 1
 1
 0
 0

julia> W=rootdatum(:E6sc);l=quasi_isolated_reps(W)
7-element Vector{SemisimpleElement{Root1}}:
 <1,1,1,1,1,1>
 <-1,1,1,-1,1,-1>
 <ζ₃,1,ζ₃²,1,ζ₃,ζ₃²>
 <ζ₃²,1,ζ₃,1,ζ₃,ζ₃²>
 <ζ₃²,1,ζ₃,1,ζ₃²,ζ₃>
 <ζ₆⁵,1,ζ₃²,1,ζ₃,ζ₃²>
 <ζ₃²,1,ζ₃,1,ζ₃²,ζ₆⁵>

julia> map(s->isisolated(W,s),l)
7-element Vector{Bool}:
 1
 1
 1
 1
 1
 1
 1

julia> Semisimple.quasi_isolated_reps(W,3)
2-element Vector{SemisimpleElement{Root1}}:
 <1,1,1,1,1,1>
 <-1,1,1,-1,1,-1>
```
"""
function quasi_isolated_reps(W::FiniteCoxeterGroup,p=0)
##  This function follows Theorem 4.6 in
##  C.Bonnafe, ``Quasi-Isolated Elements in Reductive Groups''
##  Comm. in Algebra 33 (2005), 2315--2337
##  after one fixes the following bug: at the beginning of section 4.B
##  ``the stabilizer of `Ω ∩ Δ̃ᵢ` in `𝓐 _G` acts transitively on `Ω ∩ Δ̃ᵢ`''
##  should be
##  ``the stabilizer of `Ω` in `𝓐 _G` acts transitively on `Ω ∩ Δ̃ᵢ`''
  if istorus(W) return [ss(W,fill(0//1,rank(W)))] end
  H=fundamental_group(W)
  w=Vector{Vector{Rational{Int}}}[]
  ind=Vector{Int}[]
  iso=coweights(W)
  l=map(refltype(W))do t
    n=t.indices; # n is Δₜ
    # next line  uses that negative roots are listed by decreasing height!
    r=findlast(rr->sum(rr[n])!=0,W.rootdec)
    d=inclusion(W,vcat(n,[r])) # d is Δ̃ₜ
    push!(ind,d)
    push!(w,toL(vcat(iso[n,:].//(-W.rootdec[r][n]),0*iso[1:1,:])))
    pp=vcat(map(i->combinations(d,i),1:length(H))...)
    filter(P->length(orbits(stabilizer(H,P,onsets),P))==1,pp) #possible sets Ωₜ
  end
  res=splat(vcat).(tcartesian(l...))
  res=filter(res)do P
    S=stabilizer(H,P,onsets)
    all(I->length(orbits(S,intersect(P,I)))==1,ind)
  end
  res=map(x->x[1],orbits(H,map(x->unique!(sort(x)),res),
          (s,g)->unique!(sort(s.^g)))) # possible sets Ω
  if p!=0
    res=filter(res)do P
      all(map(ind,w)do I,W
        J=intersect(P,I)
        length(J)%p!=0 && all(v->lcm(denominator.(v))%p!=0,
                              W[map(x->findfirst(==(x),I),J)])
      end)
    end
  end
  res=map(res)do P
      sum(map(ind,w)do I,p
      J=intersect(P,I)
      sum(p[map(x->findfirst(==(x),I),J)])//length(J)
     end)
  end
  res=sort(unique!(map(s->ss(W,s),res)),by=x->(order(x),x))
  Z0=algebraic_center(W).Z0
  if rank(Z0)>0
    res=res[filter(i->!any(j->res[i]/res[j] in Z0,1:i-1),eachindex(res))]
  end
  res
end

isisolated(W,s)=rank(algebraic_center(centralizer(W,s).group).Z0)==
    rank(W)-semisimplerank(W)

"""
`structure_rational_points_connected_centre(W,q)`

`W`  should be  a Coxeter  group or  a Coxeter  coset representing a finite
reductive  group ``𝐆 ^F``, and `q` should  be the prime power associated to
the  isogeny `F`. The function returns the abelian invariants of the finite
abelian group ``Z⁰𝐆 ^F`` where `Z⁰𝐆 ` is the connected center of `𝐆 `.

In  the following example one determines the structure of `𝐓(𝔽₃)` where `𝐓`
runs over all the maximal tori of `SL`₄.

```julia-repl
julia> l=twistings(rootdatum(:sl,4),Int[])
5-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 A₃₍₎=Φ₁³
 A₃₍₎=Φ₁²Φ₂
 A₃₍₎=Φ₁Φ₂²
 A₃₍₎=Φ₁Φ₃
 A₃₍₎=Φ₂Φ₄

julia> structure_rational_points_connected_centre.(l,3)
5-element Vector{Vector{Int64}}:
 [2, 2, 2]
 [2, 8]
 [4, 8]
 [26]
 [40]
```
"""
function structure_rational_points_connected_centre(MF,q)
  if MF isa Spets M=Group(MF)
  else M=MF;MF=Spets(M)
  end
  W=parent(M)
  Z0=algebraic_center(M).Z0
  Phi=YMatrix(W.G,MF.phi)
  Z0F=Z0.gens*(Phi*q-I)
  Z0F=map(x->solutionmatInt(Z0.gens,x),eachrow(Z0F))
  Z0F=smith(toM(Z0F))
  filter(!isone,map(i->Z0F[i,i],axes(Z0F,1)))
end

"""
`semisimple_centralizer_representatives(W [,p])` or `sscentralizer_reps`

`W`  should be a Weyl group corresponding  to an algebraic group `𝐆 `. This
function  returns a list describing representatives  `𝐇 ` of `𝐆 `-orbits of
reductive  subgroups  of  `𝐆 `  which  are  the  identity component of the
centralizer of a semisimple element. Each group `𝐇 ` is specified by a list
`h`   of  reflection  indices  in  `W`   such  that  `𝐇  `  corresponds  to
`reflection_subgroup(W,h)`.  If a  second argument  `p` is  given, only the
list of the centralizers which occur in characteristic `p` is returned.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> sscentralizer_reps(W)
6-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [1, 2]
 [1, 5]
 [2, 6]

julia> reflection_subgroup.(Ref(W),sscentralizer_reps(W))
6-element Vector{FiniteCoxeterSubGroup{Perm{Int16},Int64}}:
 G₂₍₎=Φ₁²
 G₂₍₁₎=A₁Φ₁
 G₂₍₂₎=Ã₁Φ₁
 G₂
 G₂₍₁₅₎=A₂
 G₂₍₂₆₎=Ã₁×A₁

julia> sscentralizer_reps(W,2)
5-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [1, 2]
 [1, 5]
```
"""
function semisimple_centralizer_representatives(W,p=0)
# W-orbits of subsets of Π∪ {-α₀}
  l=map(refltype(W))do t
    H=reflection_subgroup(W,t.indices)
    cent=reflection_subgroup.(Ref(H),parabolic_reps(H))
    npara=length(cent)
    ED=vcat(1:rank(t),[nref(H)])
    for J in combinations(ED)
      if length(J)==length(ED) continue end
      R=reflection_subgroup(H,J)
      if !isnothing(standard_parabolic(H,R)) continue end
      u=findall(G->isomorphism_type(R)==isomorphism_type(G),cent[npara+1:end])
#     if length(u)>0 xprintln("comparing ",R," to ",cent[npara+1:end][u]) end
      if all(G->isnothing(transporting_elt(H,R,G)),cent[npara.+u])
        push!(cent,R)
      end
    end
    cent=inclusiongens.(cent)
    if p==0 return cent end
    filter(I->all(x->x==0 || x%p!=0, smith(toM(W.rootdec[I]))),cent)
  end
  if isempty(l) return [Int[]] end
  splat(vcat).(tcartesian(l...))
end

const sscentralizer_reps=semisimple_centralizer_representatives

end
