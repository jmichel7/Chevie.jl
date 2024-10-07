"""
# Algebraic groups and semi-simple elements

Let  us fix an  algebraically closed field  `K` and let  `𝐆` be a connected
reductive  algebraic group over `K`. Let `𝐓` be a maximal torus of `𝐆`, let
`X(𝐓)`  be the  character group  of `𝐓`  (resp. `Y(𝐓)`  the dual lattice of
one-parameter  subgroups  of  `𝐓`)  and  `Φ`  (resp  `Φ^`) the roots (resp.
coroots) of `𝐆` with respect to `𝐓`.

Then  `𝐆` is  determined up  to isomorphism  by the  *root datum* `(X(𝐓),Φ,
Y(𝐓),Φ^)`.  In algebraic terms, this consists  in giving a free `ℤ`-lattice
`X=X(𝐓)` of dimension the *rank* of `𝐓` (which is also called the *rank* of
`𝐆`),  and a root system `Φ ⊂ X`,  and similarly giving the dual roots `Φ^⊂
Y=Y(𝐓)`.

This  is obtained  by a  slight generalisation  of our  setup for a Coxeter
group `W`. This time we assume that the canonical basis of the vector space
`V`  on which `W`  acts is a  `ℤ`-basis of `X`,  and `Φ` is  specified by a
matrix  `simpleroots(W)` whose rows are the  simple roots expressed in this
basis  of `X`. Similarly  `Φ^` is described  by a matrix `simplecoroots(W)`
whose  rows are the simple  coroots in the basis  of `Y` dual to the chosen
basis of `X`. The duality pairing between `X` and `Y` is the canonical one,
that  is  the  pairing  between  vectors  `x∈  X`  and  `y∈  Y` is given by
`transpose(x)*y`. Thus, we must have the relation
`simplecoroots(W)*permutedims(simpleroots(W))=cartan(W)`.

We  get this by using the function `rootdatum`, whose arguments are the two
matrices `simpleroots(W)` and `simplecoroots(W)` described above. The roots
need not generate `V`, so the matrices need not be square. For example, the
root datum of the linear group of rank 3 can be obtained as:

```julia-repl
julia> W=rootdatum([-1 1 0;0 -1 1],[-1 1 0;0 -1 1])
A₂Φ₁

julia> reflrep(W,W(1))
3×3 Matrix{Int64}:
 0  1  0
 1  0  0
 0  0  1
```
here  the  Weyl  group  is  the  symmetric  group  on  3  letters acting by
permutation of the basis of `X`. The dimension of `X`,
`size(simpleroots(W),2)`,  is the *rank* and  the dimension of the subspace
generated   by   the   roots,   `size(simpleroots(W),1)`,   is  called  the
*semi-simple rank*. In the example the rank is 3 and the semisimple rank is
2.

The  default form  `W=coxgroup(:A,2)` defines  the adjoint  algebraic group
(the group in its isogeny class with a trivial centre). In this case `Φ` is
a   basis  of  `X`,   so  `simpleroots(W)`  is   the  identity  matrix  and
`simplecoroots(W)` is the Cartan matrix `cartan(W)` of the root system.

The   form  `coxgroup(:A,2,sc=true)`   constructs  the   semisimple  simply
connected  algebraic  group,  where  `simpleroots(W)`  is the transposed of
`cartan(W)` and `simplecoroots(W)` is the identity matrix.

An  extreme form of root data can  be specified more conveniently: when `W`
is  the trivial `coxgroup()` (in this case `𝐆 ` is a torus), the root datum
has  no roots, thus  is entirely determined  by the rank  `r`. The function
`torus(r)`  constructs  such  a  root  datum  (it could be also obtained by
giving two `0×r` matrices to `rootdatum`).

Finally,  the `rootdatum` function also understands some familiar names for
the  algebraic groups for which it gives the results that could be obtained
by giving the appropriate `simpleroots(W)` and `simplecoroots(W)` matrices:
```julia-repl
julia> rootdatum(:gl,3)   # same as the previous example
gl₃
```
The types of root data which are understood are
 `:gl, :pgl, :sl, :slmod, :tgl :sp, :csp, :psp, :so, :pso, :cso, :halfspin, 
  :gpin, :spin, :E6, :E6sc, :CE6, :E7, :E7sc, :CE7, :E8, :F4, :G2`.

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
export algebraic_center, SubTorus, weightinfo, fundamental_group, isisolated,
SemisimpleElement, ss, torsion_subgroup, quasi_isolated_reps,
structure_rational_points_connected_centre, 
semisimple_centralizer_representatives,sscentralizer_reps, 
intermediate_group,
isomorphism_type, weights, coweights, affine
export ExtendedCox, ExtendedReflectionGroup
#----------------- Extended Coxeter groups-------------------------------
struct ExtendedCox{T,TW<:FiniteCoxeterGroup{T}}<:Group{T}
  group::TW
  F0s::Vector{Matrix{Int}}
  phis::Vector{T}
end

Groups.gens(W::ExtendedCox)=vcat(gens(W.group),W.phis)

function ExtendedCox(W::FiniteCoxeterGroup{T},F0s::Vector{<:Matrix})where T
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

ExtendedReflectionGroup(W,mats::AbstractVector{Matrix{Int}})=ExtendedCox(W,mats)
ExtendedReflectionGroup(W,mats::Matrix{Int})=ExtendedCox(W,[mats])
ExtendedReflectionGroup(W,mats::AbstractVector{<:AbstractVector{Int}})=ExtendedCox(W,[toM(mats)])
ExtendedReflectionGroup(W)=ExtendedReflectionGroup(W,Matrix{Int}[])

function ExtendedReflectionGroup(W,mats::Vector{Vector{Vector{Int}}})
  if isempty(mats)  ExtendedCox(W,empty([fill(0,0,0)]))
  elseif isempty(mats[1]) ExtendedCox(W,fill(fill(0,0,0),length(mats)))
  else ExtendedCox(W,toM.(mats))
  end
end

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
SemisimpleElement(W,map(x->Root1(;r=Rational{Int}(x)),v))

ss(W::FiniteCoxeterGroup)=SemisimpleElement(W,fill(E(1),rank(W)))

Base.:^(a::SemisimpleElement,n::Integer)=SemisimpleElement(a.W,a.v .^n)

Base.:^(a::SemisimpleElement,m::AbstractMatrix)=SemisimpleElement(a.W,
                                 map(v->prod(a.v .^v),eachcol(m)))

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

struct SubTorus{T<:Integer}
  gens::Matrix{T}
  complement::Matrix{T}
  group
end

"""
`SubTorus(W,Y::Matrix{<:Integer})`

The  function  returns  the  subtorus  𝐒  of  the  maximal torus `𝐓` of the
reductive  group represented by the Weyl group  `W` such that `Y(𝐒)` is the
(pure)  sublattice of  `Y(𝐓)` generated  by the  (integral) vectors  `Y`. A
basis  of `Y(𝐒)` in  Hermite normal form  (for easy memebership testing) is
computed  and stored in the field `S.gens` of the returned SubTorus struct.
A  basis of `Y(𝐓)` adapted to `Y(𝐒)` is  also computed. This means a set of
integral   vectors,  stored  in  `S.complement`,   is  computed  such  that
`M=vcat(S.gens,S.complement)`   is   a   basis   of   `Y(𝐓)`  (equivalently
`M∈GL(Z^{rank(W)})`.  An  error  is  raised  if  `Y` does not define a pure
sublattice.

```julia-repl
julia> W=coxgroup(:A,4)
A₄

julia> SubTorus(W,[1 2 3 4;2 3 4 1;3 4 1 1])
SubTorus(A₄,[1 0 3 -13; 0 1 2 7; 0 0 4 -3])

julia> SubTorus(W,[1 2 3 4;2 3 4 1;3 4 1 2])
ERROR: not a pure sublattice
Stacktrace:
 [1] error(::String) at ./error.jl:33
 [2] Chevie.Weyl.SubTorus(::FiniteCoxeterGroup{Perm{Int16},Int64}, ::Matrix{Int64}) at /home/jmichel/julia/Chevie.jl/src/Weyl.jl:1082
 [3] top-level scope at REPL[25]:1
```
"""
function SubTorus(W,V::Matrix{<:Integer}=reflrep(W,one(W)))
  C=complementInt(V)
  if any(!=(1),C.moduli) error("not a pure sublattice") end
  SubTorus(baseInt(C.sub),C.complement,W)
end

Base.show(io::IO,T::SubTorus)=print(io,"SubTorus(",T.group,",",T.gens,")")

Chevie.rank(T::SubTorus)=size(T.gens,1)

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
  ssl=map(toAZ,gens(descAZ))
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

function WeightToAdjointFundamentalGroupElement(W,l::Vector;full=false)
  if isempty(l) return Perm();end
  prod(x->WeightToAdjointFundamentalGroupElement(W,x;full),l)
end

function WeightToAdjointFundamentalGroupElement(W,i;full=false)
  t=refltype(W)[findfirst(t->i in t.indices,refltype(W))]
  l=copy(t.indices)
  b=longest(W,l)*longest(W,setdiff(l,[i]))
  push!(l,maximum(findall(
    i->all(j->j in t.indices || W.rootdec[i][j]==0,1:semisimplerank(W)),
  eachindex(W.rootdec))))
  full ? b : restricted(b,inclusion.(Ref(W),l))
end

"""
`weightinfo(W)`

`W`  is a  Coxeter group  record describing  an algebraic  group `𝐆 `, or a
`IypeIrred`. The function is independent of the isogeny type of `𝐆`(so just
depends  on `refltype(W)`, that is  on the root system).  It returns a dict
with the following keys:

  - `:minusculeWeights` the minuscule weights, described as their position 
    in  the  list  of  fundamental  weights.  For non-irreducible groups, a
    weight  is the sum of at most one weight in each irreducible component.
    It  is represented as  the list of  its components. For consistency, in
    the  case  of  an  irreducible  system,  a  weight  is represented as a
    one-element list.

  - `:minusculeCoweights` the minuscule coweights, represented in the same
    manner as the minuscule weights

  - `:decompositions` for each coweight, its decomposition in terms of the
    generators  of the adjoint fundamental group  (given by the list of the
    exponents  of the generators). Together with  the next field it enables
    to work out the group structure of the adjoint fundamental group.

  - `:moduli` the list of orders of the generators of the fundamental group.

  - `:AdjointFundamentalGroup` the list of generators of the adjoint
    fundamental group, given as permutations of the extended root basis.

  - `:CenterSimplyConnected` A list of semisimple elements generating the
    center of the universal covering of  𝐆

  - `:chosenAdaptedBasis` A basis  of the  weight lattice  adapted to the
    root lattice. In the basis of the fundamental weights, the root lattice
    is  given  by  `C=transpose(cartan(W))`.  The  property  specifying the
    `:chosenAdaptedBasis` is that the Hermite normal form of
    `C*chosenAdaptedBasis`  is almost in Smith  normal form (it is diagonal
    but  the diagonal entries may be  permuted compared to the Smith normal
    form; the non-trivial entries are in the positions corresponding to the
    generators of the fundamental group as indicated by `:decompositions`).

```julia-repl
julia> weightinfo(coxgroup(:A,2)*coxgroup(:B,2))
Dict{Symbol, Array} with 8 entries:
  :moduli                  => [3, 2]
  :minusculeWeights        => [[1, 3], [1], [2, 3], [2], [3]]
  :decompositions          => [[2, 1], [2, 0], [1, 1], [1, 0], [0, 1]]
  :highestroot             => [5, 7]
  :chosenAdaptedBasis      => [1 -1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
  :minusculeCoweights      => [[1, 4], [1], [2, 4], [2], [4]]
  :CenterSimplyConnected   => Vector{Rational{Int64}}[[1//3, 2//3, 0, 0], [0, 0…
  :AdjointFundamentalGroup => [(1,12,2), (4,14)]
```
"""
function weightinfo(t::TypeIrred)
  r=getchev(t,:WeightInfo)
  if isnothing(r)
    r=Dict{Symbol,Any}(:moduli=>Int[],:decompositions=>Vector{Vector{Int}}[],
         :minusculeWeights=>Vector{Int}[])
  end
  if !haskey(r,:minusculeCoweights)
    r[:minusculeCoweights]=r[:minusculeWeights]
  end
  if !haskey(r,:chosenAdaptedBasis)
   r[:chosenAdaptedBasis]=Matrix{Int}(I,rank(t),rank(t))
  else
   r[:chosenAdaptedBasis]=toM(r[:chosenAdaptedBasis])
  end
  r
end

function weightinfo(W)
  M=Matrix{Int}(I,semisimplerank(W),semisimplerank(W))
  if isempty(refltype(W)) return Dict(:minusculeWeights=>Vector{Int}[],
         :minusculeCoweights=>Vector{Int}[],
         :decompositions=>Vector{Vector{Int}}[],
         :moduli=>Int[],
         :highestroot=>Int[],
         :CenterSimplyConnected=>Vector{Rational{Int}}[],
         :AdjointFundamentalGroup=>eltype(W)[]
        )
  end
  l=map(refltype(W)) do t
    r=weightinfo(t)
    if isempty(r[:moduli]) g=Int[]
      r[:ww]=eltype(W)[]
    else g=filter(i->sum(r[:decompositions][i])==1,
          eachindex(r[:minusculeCoweights])) # generators of fundamental group
      r[:ww]=map(x->WeightToAdjointFundamentalGroupElement(W,x),
               t.indices[r[:minusculeCoweights][g]])
    end
    r[:csi]=zeros(Rational{Int},length(g),semisimplerank(W))
    if !isempty(r[:moduli])
      C=modZ.(inv(Rational.(cartan(t))))
      r[:csi][:,t.indices]=C[r[:minusculeCoweights][g],:]
      r[:minusculeWeights]=t.indices[r[:minusculeWeights]]
      r[:minusculeCoweights]=t.indices[r[:minusculeCoweights]]
      M[t.indices,t.indices]=r[:chosenAdaptedBasis]
    end
    r[:csi]=Array.(toL(r[:csi]))
    r[:highestroot]=last(findmax(x->sum(x[indices(t)]),W.rootdec))
    r
  end
  res=Dict(:minusculeWeights=>cartesian(map(
                                   x->vcat(x[:minusculeWeights],[0]),l)...),
    :minusculeCoweights=>cartesian(map(
                                   x->vcat(x[:minusculeCoweights],[0]),l)...),
    :decompositions=>map(x->vcat(x...),cartesian(map(x->vcat(x[:decompositions],
                                 [0 .*x[:moduli]]),l)...)),
    :moduli=>reduce(vcat,map(x->x[:moduli],l)),
# center of simply connected group: the generating minuscule coweights
# mod the root lattice
    :CenterSimplyConnected=>reduce(vcat,getindex.(l,:csi)),
    :AdjointFundamentalGroup=>reduce(vcat,getindex.(l,:ww)),
    :chosenAdaptedBasis=>M,
    :highestroot=>map(x->x[:highestroot],l))
  n=1:length(res[:decompositions])-1
  res[:minusculeWeights]=map(x->filter(!iszero,x),res[:minusculeWeights][n])
  res[:minusculeCoweights]=map(x->filter(!iszero,x),res[:minusculeCoweights][n])
  res[:decompositions]=res[:decompositions][n]
  res
end

" `weights(W)` simple weights in the basis of X(T)"
weights(W)=permutedims(inv(Rational.(cartan(W))))*simpleroots(W)

" `coweights(W)` simple coweights in the basis of Y(T)"
coweights(W)=inv(Rational.(cartan(W)))*simplecoroots(W)

"""
`fundamental_group(W)`

This  function returns the fundamental group of the algebraic group defined
by  the Coxeter  group struct  `W`. This  group is  returned as  a group of
diagram  automorphisms  of  the  corresponding  affine Weyl group, which is
represented  as a group of permutations of the set of simple roots enriched
by the lowest root of each irreducible component. The definition we take of
the  fundamental group of a (not necessarily semisimple) reductive group is
(P∩ Y(𝐓))/Q where P is the coweight lattice (the dual lattice in Y(𝐓)⊗ ℚ of
the  root  lattice)  and  Q  is  the  coroot  latice. The bijection between
elements  of P/Q and  diagram automorphisms is  explained in the context of
non-irreducible groups for example in [§3.B Bonnafé2005](biblio.htm#Bon05).

```julia-repl
julia> W=coxgroup(:A,3) # the adjoint group
A₃

julia> fundamental_group(W) # the 12th root is the lowest one
Group((1,12,3,2))

julia> W=rootdatum(:sl,4) # the semisimple simply connected group
sl₄

julia> fundamental_group(W)
Group(Perm{Int16}[])
```
"""
function fundamental_group(W;full=false)
  if istorus(W) return Group(Perm()) end
  e=weightinfo(W)[:minusculeCoweights]
  e=filter(x->all(isinteger,sum(coweights(W)[x,:];dims=1)),e) # minusc. coweights in Y
  if isempty(e) return Group(Perm()) end
  e=map(x->WeightToAdjointFundamentalGroupElement(W,x;full),e)
  Group(abelian_gens(e))
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
   t=deepcopy(t)
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
  Id=vcat(Matrix{Int}(I,r,r),fill(0,r)')
  mov=map(v->AffineRootAction(W,w,v)-v,eachrow(Id))
  l=push!(map(i->refls(W0,i),eachindex(gens(W0))),refls(W0,nref(W0)))
  p=reflength(W0,prod(l[word(W,w)]))
  dimw=Minimum(List(Filtered(ParabolicSubgroups(W0),
      x->RankMat(Concatenation(mov,roots(W0,x)))==Length(x)),Length))
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
given  in  [Bonnafe2005](biblio.htm#Bon05).  If  a  second  argument `p` is
given,  it  gives  representatives  of  those quasi-isolated elements which
exist in characteristic `p`.

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
  res=map(x->vcat(x...),cartesian(l...))
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
`intermediate_group(W, indices)`

This  computes  a  `rootdatum`  representing  a  semisimple algebraic group
intermediate  between  the  adjoint  group  ---  obtained  by  a  call like
`rootdatum(:pgl,4)`---  and  the  simply  connected  semi-simple  group ---
obtained  by  a  call  like  `rootdatum(:sl,4)`.  The group is specified by
specifying  a subset  of the  *minuscule weights*,  which are weights whose
scalar  product  with  every  coroot  is  in  `-1,0,1` (the weights are the
elements  of the *weight  lattice*, the lattice  in `X(𝐓)⊗ ℚ  ` dual to the
coroot  lattice). The non-trivial characters of the (algebraic) center of a
semi-simple  simply  connected  algebraic  group  are in bijection with the
minuscule  weights; this set is  also in bijection with  `P/Q` where `P` is
the  weight lattice and `Q` is the root lattice. If `W` is irreducible, the
minuscule  weights are part of the basis of the weight lattice given by the
*fundamental  weights*, which is the dual basis of the simple coroots. They
can  thus be specified by an `index` in the Dynkin diagram. The constructed
group  has lattice `X(𝐓)` generated by the  sum of the root lattice and the
weights  with the given  `indices`. If `W`  is not irreducible, a minuscule
weight is a sum of minuscule weights in different components. An element of
`indices` is thus itself a list, interpreted as representing the sum of the
corresponding weights.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> fundamental_group(intermediate_group(W,Int[])) # adjoint
Group((1,12,3,2))

julia> fundamental_group(intermediate_group(W,Int[2])) # intermediate
Group((1,3)(2,12))
```
"""
function intermediate_group(W,I)
  C=cartan(W)
  w=transpose(inv(C*1//1)); # w = weights in terms of roots
  R=one(w)
  for v in I
    if isinteger(v) R=vcat(R,transpose(w[v,:]))
    else R=vcat(R,transpose(sum(eachrow(w[v,:]))))
    end
  end
  d=lcm(denominator.(R))
  R=baseInt(Int.(d*R))//d
  rootdatum(Int.(R^-1),Int.(C*transpose(R)))
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
  map(x->vcat(x...),cartesian(l...))
end

const sscentralizer_reps=semisimple_centralizer_representatives

function isomorphism_type(t::TypeIrred;TeX=false,limit=false)
  if !limit && !TeX context=(:TeX=>true,:limit=>false)
  else context=(:TeX=>TeX,:limit=>limit)
  end
  t=repr(t;context)
  if !limit && !TeX
    t=Format.TeXstrip(t)
    t=replace(t,"^"=>"")
  end
  t
end

function isomorphism_type(W;torus=false,TeX=false,limit=false)
  if !limit && !TeX context=(:TeX=>true,:limit=>false)
  else context=(:TeX=>TeX,:limit=>limit)
  end
  t=reverse(tally(map(x->isomorphism_type(x;TeX,limit),refltype(W))))
  t=join(map(x-> x[2]==1 ? x[1] : string(x[2],x[1]),t),"+")
  d=rank(W)-semisimplerank(W)
  if d>0 && torus
    if t!="" t*="+" end
    t*="T"*string(d)
  end
  t
end

end
