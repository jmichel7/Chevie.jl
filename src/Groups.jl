"""
This module is a port of some GAP functionality on groups.

The only field of a Group G at the start is gens, the list of generators of
G.  To  mimic  GAP  records  where  attributes/properties  of an object are
computed  on demand when asked for, other fields are computed on demand and
stored in the field prop of the Group, which starts as Dict{Symbol,Any}()

# Examples
```julia-repl
julia> G=Group([Perm(i,i+1) for i in 1:2])
Group([(1,2),(2,3)])

julia> gens(G)
2-element Array{Perm{Int64},1}:
 (1,2)
 (2,3)

julia> nbgens(G)
2

# orbit of point 1 under G
julia> orbit(G,1) 
3-element Array{Int64,1}:
 2
 3
 1

# transversal of point 1: Dict  t with keys orbit of 1 and for i in orbit
# t[i] is a representatives moving 1 to i
julia> transversal(G,1)
Dict{Int64,Perm{Int64}} with 3 entries:
  2 => (1,2)
  3 => (1,3,2)
  1 => ()

# orbit functions can take any action of G as keyword argument
julia> transversal(G,[1,2],action=(x,y)->x.^Ref(y))
Dict{Array{Int64,1},Perm{Int64}} with 6 entries:
  [1, 3] => (2,3)
  [1, 2] => ()
  [2, 3] => (1,2,3)
  [3, 2] => (1,3)
  [2, 1] => (1,2)
  [3, 1] => (1,3,2)

julia> minimal_words(G)  # minimal word in gens for each element of G
Dict{Perm{Int64},Array{Int64,1}} with 6 entries:
  ()      => Int64[]
  (2,3)   => [2]
  (1,3,2) => [1, 2]
  (1,3)   => [1, 2, 1]
  (1,2)   => [1]
  (1,2,3) => [2, 1]
```

"""

module Groups
using ..Gapjm # for gens, minimal_words
export Group, minimal_words, element, gens, nbgens, class_reps,
  conjugacy_classes, orbit, transversal, orbits, Hom

#--------------general groups and functions for "black box groups" -------
abstract type Group{T} end # T is the type of elements of G

Base.one(G::Group{T}) where T=isempty(gens(G)) ? one(T) : one(gens(G)[1])
gens(G::Group)=G.gens
nbgens(G::Group)=length(gens(G))
@inline gen(W,i)=i>0 ? gens(W)[i] : inv(gens(W)[-i])

" element of W corresponding to a sequence of generators and their inverses"
element(W::Group,w...)=isempty(w) ? one(W) : length(w)==1 ? gen(W,w[1]) : 
    prod( gen(W,i) for i in w)

" orbit(G,p) is the orbit of p under Group G"
function orbit(G::Group,p;action::Function=^)
  new=Set([p])
  res=empty(new)
  while !isempty(new)
    union!(res,new)
    n=empty(res)
    for p in new, s in gens(G)
      w=action(p,s)
      if !(w in res) push!(n,w) end
    end
    new=n
  end
  collect(res)
end

"returns Dict x=>g for x in orbit(G,p) giving g such that x=action(p,g)"
function transversal(G::Group,p;action::Function=^)
  new=[p]
  res=Dict(p=>one(G))
  while !isempty(new)
    old=copy(new)
    empty!(new)
    for p in old, s in gens(G)
      w=action(p,s)
      if !haskey(res,w)
        push!(new,w)
        res[w]=res[p]*s
      end
    end
  end
  res
end

function orbits(G::Group,v::AbstractVector=1:degree(G);action::Function=^)
  res=Vector{eltype(v)}[]
  while !isempty(v)
    o=orbit(G,v[1],action=action)
    push!(res,o)
    v=setdiff(v,o)
  end
  res
end

" dict giving for each element of G a minimal word in the generators"
function minimal_words(G::Group)
  gets(G,:words)do G
  words=Dict(one(G)=>Int[])
  for i in eachindex(gens(G))
    rw = [one(G)=>Int[]]
    j=1
    nwords=deepcopy(words)
    while j<=length(rw)
      for k in 1:i
        e=rw[j][1]*gens(G)[k]
        if !haskey(nwords,e)
          we=vcat(rw[j][2],[k])
          push!(rw,e=>we)
          for (e1,w1) in words nwords[e1*e]=vcat(w1,we) end
        end
      end
      j+=1
    end
    words = nwords
  end
  words
  end
end

Gapjm.word(G,w)=minimal_words(G)[w]

Gapjm.elements(G::Group)=collect(keys(minimal_words(G)))
Gapjm.length(G::Group)=length(elements(G))

function conjugacy_classes(G::Group{T})::Vector{Vector{T}} where T
  gets(G,:classes) do G
    if length(G)>1000
      println("length(G)=",length(G),": not supposed to do it the hard way")
    end
    cl=orbits(G,collect(G))
    G.prop[:classreps]=first.(cl)
    cl
  end
end

function class_reps(G::Group{T})::Vector{T} where T
  getp(conjugacy_classes,G,:classreps)
end

# hom from source to target sending gens(source) to images
struct Hom{T}
  source::Group{T}
  target::Group{T}
  images::Vector{T}
end

function Base.show(io::IO,h::Hom)
  print(io,"Hom(",h.source,"→ ",h.target,";",gens(h.source),"↦ ",h.images)
end

function Gapjm.kernel(h::Hom)
  if all(isone,h.images) return h.source end
  if length(h.source)==length(Group(h.images)) 
    return Group(empty(gens(h.source)))
  end
  error("not implemented: kernel(",h,")")
end

# h(w) is the image of w by h
(h::Hom)(w)=isone(w) ? one(h.target) : prod(h.images[word(h.source,w)])

end
