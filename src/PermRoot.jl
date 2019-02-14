module PermRoot

export PermRootGroup, PermRootSubGroup, ReflectionSubGroup, simple_representatives,
simple_conjugating_element, reflections, reflection, refltype, diagram

using Gapjm

struct PermRootGroup{T,T1<:Integer}
  matgens::Vector{Matrix{T}}
  roots::Vector{Vector{T}}
  coroots::Vector{Vector{T}}
  cartan::Array{T,2}
  G::PermGroup{T1}
  prop::Dict{Symbol,Any}
end

function refltype end

Gapjm.gens(W::PermRootGroup)=gens(W.G)
Base.one(W::PermRootGroup)=one(W.G)
PermGroups.orbit_and_representative(W::PermRootGroup,i)=orbit_and_representative(W.G,i)
Base.length(W::PermRootGroup)=length(W.G)
PermGroups.element(W::PermRootGroup,x...)=element(W.G,x...)

function PermRootGroup(r::Vector{Vector{T}},cr::Vector{Vector{T1}}) where{T,T1}
  matgens=map(reflection,r,cr)

  # the following section is quite subtle: it has the (essential -- this is
  # what  allows  to  construct  reflexion  subgroups  in a consistent way)
  # property  that the order of the  constructed roots (thus the generating
  # permutations) depends only on the Cartan matrix of g, not on the actual
  # root values.

  println("# roots: ")
  roots=map(x->convert.(eltype(matgens[1]),x),r)
  refls=map(x->Int[],roots)
  newroots=true
  while newroots
    newroots=false
    for j in eachindex(matgens)
      lr=length(roots)
      for y in Ref(permutedims(matgens[j])).*roots[length(refls[j])+1:end]
        p=findfirst(isequal(y),roots[1:lr]) 
	if isnothing(p)
          push!(roots,y)
#         println("j=$j roots[$(length(refls[j])+1)...] ",length(roots),":",y)
          newroots=true
          push!(refls[j],length(roots))
        else push!(refls[j],p)
	end
      end
    end
    println(" ",length(roots))
  end
  roots=map(x->convert.(eltype(matgens[1]),x),roots)
  PermRootGroup(matgens,roots,cr,PermGroup(map(Perm,refls)),Dict{Symbol,Any}())
# for r in Orbits(G,G.rootInclusion) do
#   G.orbitRepresentative{G.rootRestriction{r}}:=[1..Length(r)]*0+Minimum(r);
# od
end

function root_representatives(W::PermRootGroup)
  reps=fill(0,length(W.roots))
  repelts=fill(one(W),length(W.roots))
  for i in eachindex(gens(W))
    if iszero(reps[i])
      d=orbit_and_representative(W,i)
      for (n,e) in d 
        reps[n]=i
        repelts[n]=e
      end
    end
  end
  W.prop[:rootreps]=reps
  W.prop[:repelms]=repelts
  W.prop[:reflections]=map((i,p)->gens(W)[i]^p,reps,repelts)
end

"for each root index of simple representative"
function simple_representatives(W::PermRootGroup{T})::Vector{T} where T
  getp(root_representatives,W,:rootreps)
end
  
"for each root element conjugative representative to root"
function simple_conjugating_element(W::PermRootGroup{T,T1},i)::Perm{T1} where{T,T1}
  getp(simple_conjugating_element,W.G,:repelms)[i]
end

function reflections(W::PermRootGroup{T,T1})::Vector{Perm{T1}} where{T,T1}
  getp(root_representatives,W,:reflections)
end

" the matrix of the reflection of given root and coroot"
function reflection(root::Vector,coroot::Vector)
  root,coroot=promote(root,coroot)
  m=[i*j for i in coroot, j in root]
  one(m)-m
end

reflection(W::PermRootGroup,i)=reflections(W)[i]

function independent_roots(W::PermRootGroup)::Vector{Int}
  gets(W,:indeproots) do W
    echelon(permutedims(hcat(W.roots...)))[2]
  end
end

function baseX(W::PermRootGroup{T})::Matrix{T} where T
  gets(W,:baseX) do W
    ir=independent_roots(W)
    res=permutedims(hcat(W.roots[ir]...))
    res=hcat(res,NullSpace(hcat(W.coroots[ir]...)))
  end
end

" as Chevie's MatXPerm"
function matX(W::PermRootGroup,w)
  X=baseX(W)
  inv(X)*vcat(permutedims(hcat(W.roots[independent_roots(W).^w]...)),
            X[length(ir)+1:end,:]);
end

function Base.show(io::IO, W::PermRootGroup)
  print(io,"PermRootGroup($(length(W.roots)) roots)")
end

function cartan_coeff(W::PermRootGroup,i,j)
  v=findfirst(x->!iszero(x),W.roots[i])
  r=W.roots[j]-W.roots[j^reflection(W,i)]
  return r[v]/W.roots[i][v];
end

function cartanmat(W::PermRootGroup)
  [cartan_coeff(W,i,j) for i in eachindex(gens(W)), j in eachindex(gens(W))]
end

struct PermRootSubGroup{T}
  G::PermGroup{T}
  inclusion::Vector{Int}
  restriction::Vector{Int}
  parent::PermRootGroup{T}
  cartan::Array{T,2}
  type::Vector{NamedTuple{(:series,:indices),Tuple{Symbol,Vector{Int}}}}
  prop::Dict{Symbol,Any}
end

cartancoeff(W::PermRootSubGroup,i,j)=cartancoeff(W.parent,W.inclusion[i],W.inclusion[j])
Gapjm.gens(W::PermRootSubGroup)=gens(W.G)
Base.one(W::PermRootSubGroup)=one(W.G)

ReflectionSubGroup(W::PermRootSubGroup,I::AbstractVector{Int})=
  ReflectionSubGroup(W.parent,W.inclusion[I])

function root_representatives(W::PermRootSubGroup)
  reps=fill(0,2*W.N)
  repelts=fill(one(W),2*W.N)
  for i in eachindex(gens(W))
    if iszero(reps[i])
      d=orbit_and_representative(W.G,W.inclusion[i])
      for (n,e) in d 
        reps[W.restriction[n]]=i
        repelts[W.restriction[n]]=e
      end
    end
  end
  W.prop[:rootreps]=reps
  W.prop[:repelms]=repelts
  W.prop[:reflections]=map((i,p)->gens(W)[i]^p,reps,repelts)
end

end
