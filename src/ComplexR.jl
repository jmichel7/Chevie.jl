module ComplexR
using ..Chevie
export complex_reflection_group, crg, diagram, degrees, codegrees,
  Reflection, reflections, isdistinguished, 
  hyperplane, hyperplane_orbits, hyperplane_orbit, simple_rep,
  reflection_group, torusfactors

@forward Weyl.FC.G hyperplane_orbits,  codegrees,  degrees

# roots for the adjoint group
Chevie.roots(t::TypeIrred)=
 t.series==:ST ? getchev(t,:GeneratingRoots) : collect(eachrow(one(cartan(t))))

# coroots for the adjoint group
function PermRoot.coroots(t::TypeIrred)
  if t.series==:ST
    cr=getchev(t,:GeneratingCoRoots)
    if isnothing(cr)
      r=getchev(t,:GeneratingRoots)
      cr=coroot.(r,E.(ordergens(t)))
    end
    return map(x->convert.(Cyc{Rational{Int}},x),cr)
  end
  toL(cartan(t))
end

"""
`complex_reflection_group(STnumber)` or `crg(STnumber)`

`complex_reflection_group(p,q,r)` or `crg(p,q,r)`

The  first form of `complex_reflection_group`  returns the complex reflection
group which has Shephard-Todd number `STnumber`, see
[Shephard-Todd1954](biblio.htm#ST54).   The   second   form   returns   the
imprimitive complex reflection group `G(p,q,r)`.

```julia-repl
julia> G=complex_reflection_group(4)
G₄

julia> degrees(G)
2-element Vector{Int64}:
 4
 6

julia> length(G)
24

julia> W*coxgroup(:A,2) # how to make a non-irreducible group
G₄×A₂

julia> complex_reflection_group(1,1,3) # another way to enter A₂
A₂

julia> crg(4) # there is also a short alias
G₄
```
"""
function complex_reflection_group(i::Int)
  if i==23     coxgroup(:H,3)
  elseif i==28 coxgroup(:F,4)
  elseif i==30 coxgroup(:H,4)
  elseif i==35 coxgroup(:E,6)
  elseif i==36 coxgroup(:E,7)
  elseif i==37 coxgroup(:E,8)
  else t=TypeIrred(Dict(:series=>:ST,:ST=>i))
    PRG(roots(t),coroots(t))
  end
end

function complex_reflection_group(p,q,r)
  if !iszero(p%q) || p<=0 || r<=0 || (r==1 && q!=1)
   error("complex_reflection_group(p,q,r) must satisfy: q|p, r>0, and r=1 => q=1")
  end
  if p==1 return coxgroup(:A,r-1)
  elseif p==2
    if q==2 return coxgroup(:D,r)
    else return coxgroup(:B,r) end
  elseif p==q && r==2 return coxgroup(:I,2,p)
  end
  t=TypeIrred(Dict(:series=>:ST,:p=>p,:q=>q,:rank=>r))
  PRG(roots(t),coroots(t))
end

const crg=complex_reflection_group

# converts a type back to a group
function reflection_group(t::TypeIrred)
  if haskey(t,:orbit)
    W=reflection_group(t.orbit)
    if length(t.orbit)>1
      spets(W,reflrep(W,Perm(vcat(circshift(map(x->x.indices,
                                        refltype(W)),-1)...))*t.twist))
    else spets(W,t.twist)
    end
  elseif t.series==:ST PRG(roots(t),coroots(t))
  else C=cartan(t)
    all(isreal,C) ? rootdatum(C) : PRG(one(C),C)
  end
end

function reflection_group(l::Vector{TypeIrred})
  if isempty(l) return coxgroup() end
  prod(reflection_group.(l))
end

"""
`degrees(W::ComplexReflectionGroup)`

returns  a list  holding the  degrees of  `W` as  a reflection group on the
vector  space `V` on which  it acts. These are  the degrees `d₁,…,dₙ` where
`n`  is the dimension of  `V` of the basic  invariants of `W` in `SV`. They
reflect  various properties  of `W`;  in particular,  their product  is the
cardinality of `W`.

```julia-repl
julia> W=complex_reflection_group(30)
H₄

julia> degrees(W)
4-element Vector{Int64}:
  2
 12
 20
 30

julia> length(W)
14400
```
"""
function degrees(W::ComplexReflectionGroup)
  get!(W,:degrees)do
    vcat(fill(1,rank(W)-semisimplerank(W)),degrees.(refltype(W))...)
  end::Vector{Int}
end

torusfactors(WF::Spets)=eigmat(central_action(Group(WF),WF.F))

"""
`degrees(WF::Spets)`

Let  `W` be  the group  of the  reflection coset  `WF`, and  let `V` be the
vector  space  of  dimension  `rank(W)`  on  which `W` acts as a reflection
group.  Let  `f₁,…,fₙ`  be  the  basic  invariants  of `W` on the symmetric
algebra  `SV` of `V`;  they can be  chosen so they  are eigenvectors of the
matrix  `WF.F`. The corresponding  eigenvalues are called  the *factors* of
`F` acting on `V`; they characterize the coset --- they are equal to 1 only
for  the trivial  coset. The  *generalized degrees*  of `WF`  are the pairs
formed of the reflection degrees and the corresponding factor.

```julia-repl
julia> W=coxgroup(:E,6)
E₆

julia> WF=spets(W)
E₆

julia> phi=W(6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1);

julia> HF=subspets(WF,2:5,phi)
E₆₍₂₃₄₅₎=³D₄Φ₃

julia> diagram(HF)
ϕ acts as (1,2,4) on the component below
  O 2
  ￨
O—O—O D₄
1 3 4

julia> degrees(HF)
6-element Vector{Tuple{Int64, Cyc{Int64}}}:
 (1, ζ₃)
 (1, ζ₃²)
 (2, 1)
 (4, ζ₃)
 (6, 1)
 (4, ζ₃²)
```
"""
function degrees(W::Spets)
  get!(W,:degrees)do
    a=torusfactors(W) 
    if isempty(refltype(W)) b=Tuple{Int,Cyc{Int}}[]
    else b=vcat(degrees.(refltype(W))...) 
    end
    # separate/recombine for promotions to work
    collect(zip(vcat(fill(1,length(a)),first.(b)),vcat(a,last.(b))))
  end::Vector{Tuple{Int,Cyc{Int}}}
end

function degrees(t::TypeIrred)
  if !haskey(t,:orbit) return getchev(t,:ReflectionDegrees) end
  d=getchev(t.orbit[1],:ReflectionDegrees)
# Let t.scalar=[s₁,…,sᵣ], where r=length(t.orbit) and ζ=prod(t.scalar). let
# p  be the PhiFactor of t.twist associated  to the reflection degree dᵢ of
# t.orbit[1].   If   G0   is   the   Spets  described  by  t.orbit[1],  and
# G1:=Ennola(ζ,G0)  then G is  isomorphic to the  descent of scalars of G1.
# According to spets 1.5, a Phifactor of Ennola(ζ,G0) is ζ^dᵢ times that of
# G0;  and  by  spets  1.5  or  [Digne-Michel,  parabolic  A.1] those of an
# a-descent of scalars are ζₐʲζᵢ^{1/a} (all the a-th roots of ζᵢ).
  if order(t.twist)>1
   f=getchev(t,:PhiFactors)
   if isnothing(f) return nothing end
  else f=fill(1,length(d))
  end
  if haskey(t,:scalar)
    p=improve_type(prod(t.scalar))
    f=[f[i]*p^d[i] for i in eachindex(d)]
  end
  f=collect(zip(d,Cyc.(f)))
  a=length(t.orbit)
  if a==1 return f end
  vcat(map(f)do (d,e) map(x->(d,x),root(e,a).*E.(a,0:a-1)) end...)
end

function codegrees(t::TypeIrred)
  if !haskey(t,:orbit)
    d=getchev(t,:ReflectionCoDegrees)
    if !isnothing(d) return d
    else
      d=getchev(t,:ReflectionDegrees)
      return reverse(maximum(d).-d)
    end
  end
  d=getchev(t,:ReflectionCoDegrees)
  if isnothing(d)
    d=getchev(t.orbit[1],:ReflectionDegrees)
    a=argmax(d)
    d=reverse(d[a].-d)
    if order(t.twist)==1
      f=fill(1,length(d))
    else
      f=getchev(t,:PhiFactors)
      if isnothing(f) return f end
      f=Cyc.(improve_type(reverse(map(x->f[a]//x,f))))
    end
    d=collect(zip(d,Cyc.(f)))
  elseif order(t.twist)==1
    d=collect(zip(d,fill(Cyc(1),length(d))))
  end
  if haskey(t,:scalar)
    f=Cyc.(improve_type(prod(t.scalar)))
    d=[(deg,eps*f^deg) for (deg,eps) in d]
  end
  a=length(t.orbit)
  if a==1 return d end
  vcat(map(d)do (d,e) map(x->(d,x),root(e,a).*E.(a,0:a-1)) end...)
end

"""
`codegrees(W::ComplexReflectionGroup)`

returns  the vector of codegrees of `W`  as a reflection group on the space
`V`  of `reflrep(W)`.  These are  one less  than the  degrees of  the basic
derivations of ` W` on `SV⊗ V^vee`.

```julia-repl
julia> W=complex_reflection_group(4)
G₄

julia> codegrees(W)
2-element Vector{Int64}:
 0
 2
```
"""
function codegrees(W::ComplexReflectionGroup)
  vcat(fill(-1,rank(W)-semisimplerank(W)),collect.(codegrees.(refltype(W)))...)
end

function codegrees(W::Spets)
  get!(W,:codegrees)do
    a=inv.(torusfactors(W))
    if isempty(refltype(W)) b=Tuple{Int,Cyc{Int}}[]
    # separate/recombine for promotions to work
    else b=vcat(codegrees.(refltype(W))...) 
    end
    collect(zip(vcat(fill(-1,length(a)),first.(b)),vcat(a,last.(b))))
  end::Vector{Tuple{Int,Cyc{Int}}}
end

"""
`hyperplane_orbits(W::ComplexReflectionGroup)`

returns  a  list  of  named  tuples,  one  for each hyperplane orbit of the
reflection  group `W`. If `o` is the named tuple for such an orbit, and `s`
is  the first  element of  `gens(W)` whose  hyperplane is  in the orbit, it
contains the following fields

 `.s`:     index of `s` in `gens(W)`

 `.order`: order of s

 `.cl_s`:  for i in `1:order-1`, `position_class(W,W(s)^i)`

 `.N_s`:    size of hyperplane orbit

 `.det_s`:  for i in `1:order-1`, position in `CharTable(W)` of `detₛⁱ`, where
   `detₛ` is the linear character taking the value `det(reflrep(W,s))` on `s`
   and `1` on non-conjugate reflections.

```julia-repl
julia> W=coxgroup(:B,2)
B₂

julia> hyperplane_orbits(W)
2-element Vector{@NamedTuple{s::Int64, cl_s::Vector{Int64}, order::Int64, N_s::Int64, det_s::Vector{Int64}}}:
 (s = 1, cl_s = [2], order = 2, N_s = 2, det_s = [5])
 (s = 2, cl_s = [4], order = 2, N_s = 2, det_s = [1])
```
"""
function hyperplane_orbits(W::ComplexReflectionGroup)
  T=@NamedTuple{s::Int,cl_s::Vector{Int},order::Int,N_s::Int,det_s::Vector{Int}}
  get!(W,:hyperplane_orbits)do
  sr=simple_reps(W)
  if isempty(sr) return T[] end
  rr=refls(W)
  cr=classreps(W)
  orb=unique(sort(sr))
  class=map(orb)do s
    map(1:ordergens(W)[s]-1)do o # return position_class(W,refls(W,s)^o)
      for i in eachindex(sr)
        if sr[i]==s
          p=findfirst(==(rr[i]^o),cr)
          if p!==nothing return p end
        end
      end
      error("not found")
    end
  end
  if !allunique(class)  # when there are too many roots in rootsystem
    l=unique(indexin(class,class));class=class[l];orb=orb[l]
  end
  chars=CharTable(W).irr
  pairs=zip(orb,class)
  res=map(pairs) do (s,c)
    ord=ordergens(W)[s]
    dets=map(1:ord-1) do j
      findfirst(i->chars[i,1]==1 && chars[i,c[1]]==E(ord,j) &&
         all(p->chars[i,p[2][1]]==1 || p[1]==s,pairs),axes(chars,1))
    end
    (s=s,cl_s=c,order=ord,N_s=div(length(W),CharTable(W).centralizers[c[1]]),det_s=dets)
  end
  end::Vector{T}
end

#------------------------- Reflection(s) --------------------------------
@GapObj struct Reflection{TW<:ComplexReflectionGroup}
  W::TW
  rootno::Int
  eigen::Root1
  word::Vector{Int}
end

@doc """
`Reflection`  is a `struct`  representing a reflection  `r` in a reflection
group `W`. The fields are `W`, the index of a root for `r`, the non-trivial
eigenvalue of `r`, and a word for `r` in the generating reflections of `W`.
```julia-repl
julia> W=crg(8);

julia> r=reflections(W)[7]
Reflection(G₈,1,-1)

julia> r==Reflection(W,1,-1) # specify r with root number and eigenvalue
true

julia> Reflection(W,1) # specify with root number gets distinguished reflection
Reflection(G₈,1,ζ₄)

julia> r.eigen # the non-trival eigenvalue, as a Root1
Root1: -1

julia> root(r)
2-element Vector{Cyc{Rational{Int64}}}:
  0
 ζ₄

julia> coroot(r)
2-element Vector{Cyc{Int64}}:
    0
 -2ζ₄

julia> Matrix(r)
2×2 Matrix{Cyc{Rational{Int64}}}:
 1   0
 0  -1

julia> hyperplane(r) # the fixed hyperplane, as a rowspace
1×2 Matrix{Cyc{Rational{Int64}}}:
 1  0

julia> hyperplane(r)*Matrix(r)==hyperplane(r)
true

julia> isdistinguished(r) # r is not distinguished
false

julia> exponent(r) # which power of a distinguished reflection it is
2

julia> Perm(r)
(1,8)(2,9)(3,16)(4,15)(5,17)(6,18)(7,19)(10,22)(11,21)(12,23)

julia> hyperplane_orbit(r) # r is in the first hyperplane orbit
1

julia> position_class(r) # the index of the conjugacy class of r in W 
15

julia> simple_rep(r) # smallest root index affording a conjugate reflection
1

julia> word(r) # a word in the generators for r
2-element Vector{Int64}:
 1
 1
```
""" Reflection

function Base.show(io::IO,r::Reflection)
  print(io,"Reflection(",r.W,",",r.rootno,",",r.eigen,")")
end

LaurentPolynomials.root(r::Reflection)=roots(r.W,r.rootno)

function PermRoot.coroot(r::Reflection)
  get!(r,:coroot)do
    rr=root(r)//1
    cr=coroots(r.W,r.rootno)
    improve_type(cr*(1-r.eigen)/(transpose(cr)*rr))
  end
end

function Base.Matrix(r::Reflection)
  get!(r,:matrix)do
    reflectionMatrix(root(r),coroot(r))
  end
end

simple_rep(r::Reflection)=simple_reps(r.W)[r.rootno]

Base.exponent(r::Reflection)=Int(r.eigen.r*ordergens(r.W)[simple_rep(r)])

isdistinguished(r::Reflection)=exponent(r)==1

Perms.Perm(r::Reflection)=refls(r.W,r.rootno)^exponent(r)

function hyperplane_orbit(r::Reflection)
  findfirst(x->x.s==simple_rep(r),hyperplane_orbits(r.W))
end

function hyperplane(r::Reflection)
  get!(r,:hyperplane)do
    Matrix(rowspace(lnullspace(hcat(coroot(r)))))
  end
end

function Groups.position_class(r::Reflection)
  get!(r,:position_class)do
    hyperplane_orbits(r.W)[hyperplane_orbit(r)].cl_s[exponent(r)]
  end::Int
end

# invert word w in gens(w) to a positive word in gens(w)
function invert_word(W,w)
  if isempty(w) return w end
  lastg=0
  mul=0
  seq=Pair{Int,Int}[]
  for i in length(w):-1:1
    if w[i]==lastg mul+=1
    else
      if lastg!=0 push!(seq,lastg=>mul) end
      lastg=w[i]; mul=1
    end
  end
  if lastg!=0 push!(seq,lastg=>mul) end
  o=ordergens(W)
  [i for (i,mul) in seq for y in 1+mul:o[i]]
end

function Reflection(W::ComplexReflectionGroup,i::Integer,eig)
  if !(eig isa Root1) eig=Root1(eig) end
  p=findfirst(r->refls(W,r.rootno)==refls(W,i) && r.eigen==eig,reflections(W))
  reflections(W)[p]
end

function Reflection(W::ComplexReflectionGroup,i::Integer)
  Reflection(W,i,E(order(refls(W,i))))
end

"""
`reflections(W)`  the list of  all reflections of  the reflection group `W`
(including  the  non-distinguished  ones),  given as a `Vector{Reflection}`
(see   [`Reflection`](@ref)).   `reflections(W)[1:nhyp(W)]`   contains  the
distinguished reflections.

```julia-repl
julia> W=crg(4)
G₄

julia> reflections(W)
8-element Vector{Reflection{PRG{Cyc{Rational{Int64}}, Int16}}}:
 Reflection(G₄,1,ζ₃)
 Reflection(G₄,2,ζ₃)
 Reflection(G₄,4,ζ₃)
 Reflection(G₄,5,ζ₃)
 Reflection(G₄,1,ζ₃²)
 Reflection(G₄,2,ζ₃²)
 Reflection(G₄,4,ζ₃²)
 Reflection(G₄,5,ζ₃²)
```
"""
function reflections(W::ComplexReflectionGroup)
  get!(W,:reflections)do
    sreps=sort(unique(simple_reps(W)))
    pnts=refls(W,sreps)
    if W isa PermRootGroup
      dd=map(x->Groups.words_transversal(gens(W),x),pnts)
    end
    res=map(i->Reflection{typeof(W)}[],1:maximum(ordergens(W))-1)
    for i in unique_refls(W)
      e=ordergens(W)[simple_reps(W)[i]]
      if W isa CoxeterGroup w=word(W,refls(W,i))
      else
        rep=simple_reps(W)[i]
        w=dd[findfirst(==(rep),sreps)][refls(W,i)]
        w=vcat(invert_word(W,w),[rep],w)
      end
      for j in 1:e-1
        push!(res[j],Reflection(W,i,E(e,j),repeat(w,j),Dict{Symbol,Any}()))
      end
    end
    vcat(res...)
   end::Vector{Reflection{typeof(W)}}
end
  
Groups.word(r::Reflection)=r.word

end
