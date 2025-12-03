#------------------------- Reflection(s) --------------------------------
export Reflection, reflections, isdistinguished, hyperplane, hyperplane_orbit,
  simple_rep

@GapObj struct Reflection{TW}
  W::TW
  rootno::Int
  eigen::Root1
  word::Vector{Int}
end

@doc """
`Reflection` is a `struct` representing a reflection in a reflection group.
```julia-repl
julia> W=crg(8);

julia> r=reflections(W)[7] # shows (r.W,r.rootno,r.eigen)
Reflection(G₈,1,-1)

julia> r.rootno # r is a reflection for the first root
1

julia> r.eigen # the non-trival eigenvalue, as a Root1
Root1: -1

julia> r.W # the group of which r is a reflection
G₈

julia> r==Reflection(W,1,-1) # specify r with .rootno and .eigen
true

julia> Reflection(W,1) # specify with .rootno gets the distinguished reflection
Reflection(G₈,1,ζ₄)

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

julia> word(r) # a word in the generators of r.W for r
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

Perms.order(r::Reflection)=order(r.eigen) 

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

function Reflection(W::PermRootGroup,i::Integer,eig)
  if !(eig isa Root1) eig=Root1(eig) end
  p=findfirst(r->refls(W,r.rootno)==refls(W,i) && r.eigen==eig,reflections(W))
  reflections(W)[p]
end

function Reflection(W::PermRootGroup,i::Integer)
  Reflection(W,i,E(order(refls(W,i))))
end

"""
`reflections(W)`   a  `Vector{Reflection}`   of  all   reflections  of  the
reflection   group   `W`   (including   the   non-distinguished  ones;  see
[`Reflection`](@ref)).  `reflections(W)[1:nhyp(W)]`  are  the distinguished
reflections.

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
function reflections(W::PermRootGroup)
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
