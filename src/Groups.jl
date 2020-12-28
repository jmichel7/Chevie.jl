"""
This module gives some basic functionality on groups.

`Group`  is  an  abstract  type,  but  the  following is assumed of all its
concrete implementations:

  - The function `gens(G)` returns the list of generators of the group `G`. 
  - The function `one(G)` returns the identity element of `G`.
  
Further,  a group  is the  kind of  object where  attributes/properties are
computed on demand when asked for; such attributes when computed are stored
in the field `.prop` of `G`, which should be of type `Dict{Symbol,Any}()`.

# Examples
```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3)])
Group([(1,2),(1,2,3)])

julia> gens(G)
2-element Array{Perm{Int16},1}:
 (1,2)  
 (1,2,3)

julia> ngens(G)
2

julia> minimal_words(G)
Dict{Perm{Int16},Array{Int64,1}} with 6 entries:
  ()      => Int64[]
  (1,3)   => [1, 2]
  (1,2,3) => [2]
  (1,2)   => [1]
  (2,3)   => [2, 1]
  (1,3,2) => [1, 2, 1]

julia> G.prop
Dict{Symbol,Any} with 1 entry:
  :words => Dict(()=>[],(1,3)=>[1, 2],(1,2,3)=>[2],(1,2)=>[1],(2,3)=>[2, 1],(1,…
```
"""
module Groups
export Group, minimal_words, gens, ngens, classreps, centralizer,
  conjugacy_classes, orbit, transversal, orbits, Hom, isabelian,
  position_class, fusion_conjugacy_classes, Coset, transporting_elt,
  centre, normalizer, stabilizer, abelian_generators,iscyclic,comm,
  nconjugacy_classes

using ..Util: gets, InfoChevie
#import Gapjm: word, elements, kernel, order
# to use as a stand-alone module comment above line and uncomment next line
export word, elements, kernel, order
#--------------general groups and functions for "black box groups" -------
abstract type Group{T} end # T is the type of elements of G

function normalizer end

Base.one(G::Group{T}) where T=one(T)

gens(G::Group)=G.gens
ngens(G::Group)=length(gens(G))

"""
`(G::Group)(i...)`

A  Group  may  be  used  as  a  function  wich  takes  integer arguments in
`eachindex(gens(W))`.  This constructs  the element  of `G`  product of the
generators  with the specified  indices. An argument  can also be negative,
then the inverse of the corresponding generator is used.

```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3)])
Group([(1,2),(1,2,3)])

julia> G(2,1,-2) # returns gens(G)[2]*gens(G)[1]/gens(G)[2]
(1,3)
```
"""
(W::Group)(w::Integer...)=isempty(w) ? one(W) : prod((i>0 ? gens(W)[i] : inv(gens(W)[-i]))
                                            for i in w)

comm(a,b)=inv(a)*inv(b)*a*b

"""
`orbit(gens::vector,p;action::Function=^)`

`orbit(G::Group,p;action::Function=^)`

the  orbit of point  `p` under repeated  action of generators `gens`. Point
`p`  should be hashable. The default action  of a group element is `^`. For
example  if `g` is a permutation and `p`  an integer, `p^g` is the image of
`p`  by `g`; if `h` and `g` are group elements, then `h^g` is the conjugate
`inv(g)*h*g`.  If a group  is given instead  of generators, the orbit under
`gens(G)` is returned.

```julia-repl
julia> orbit([Perm(1,2),Perm(2,3)],1) 
3-element Array{Int64,1}:
 1
 2
 3

julia> orbit([Perm(1,2),Perm(2,3)],[1,3];action=(v,g)->v.^g)
6-element Array{Array{Int64,1},1}:
 [1, 3]
 [2, 3]
 [1, 2]
 [3, 2]
 [2, 1]
 [3, 1]

julia> orbit([Perm(1,2),Perm(2,3)],[1,3];action=(v,g)->sort(v.^g))
3-element Array{Array{Int64,1},1}:
 [1, 3]
 [2, 3]
 [1, 2]
```
"""
function orbit(gens::Vector,pnt;action::F=^) where F<:Function
  set=Set([pnt])
  orb=[pnt]
  for pnt in orb, gen in gens
    img=action(pnt,gen)
    if !(img in set)
      push!(orb,img)
      push!(set,img)
    end
  end
  orb
end

orbit(G::Group,pnt;action::F=^) where F<:Function=orbit(gens(G),pnt;action)

"""
`transversal(G::Group,p;action::Function=^)`

returns  a `Dict` with entries `x=>g`  where `x` runs over `orbit(G,p)` and
where `g` is such that `x=action(p,g)`

```julia-repl
julia> G=Group([Perm(1,2),Perm(2,3)]);
julia> transversal(G,1)
Dict{Int64,Perm{Int16}} with 3 entries:
  2 => (1,2)
  3 => (1,3,2)
  1 => ()
```
orbit functions can take any action of `G` as keyword argument

```julia-repl
julia> transversal(G,[1,2],action=(x,y)->x.^y)
Dict{Array{Int64,1},Perm{Int16}} with 6 entries:
  [1, 3] => (2,3)
  [2, 1] => (1,2)
  [1, 2] => ()
  [3, 1] => (1,3,2)
  [2, 3] => (1,2,3)
  [3, 2] => (1,3)
```
"""
function transversal(G::Group,pnt;action::F=^) where F<:Function
  trans=Dict(pnt=>one(G))
  orb=[pnt]
  for pnt in orb, gen in gens(G)
    img=action(pnt,gen)
    if !haskey(trans,img)
      push!(orb,img)
      trans[img]=trans[pnt]*gen
    end
  end
  trans
end

function orbits(gens::Vector,v;action::F=^,
                trivial=true) where F<:Function
  res=Vector{eltype(v)}[]
  while !isempty(v)
    o=orbit(gens,first(v);action)
    if length(o)>1 || trivial push!(res,o) end
    v=setdiff(v,o)
  end
  res
end

"""
`orbits(gens::Vector,v;action=^)`

`orbits(G,v;action=^)`
    
the  orbits on `v`  of the repeated  action of `gens`;  the elements of `v`
should  be hashable. If a  group is given instead  of generators, the orbit
under `gens(G)` is returned.

```julia-repl
julia> G=Group([Perm(1,2),Perm(2,3)]);
julia> orbits(G,1:4)
2-element Array{Array{Int64,1},1}:
 [1, 2, 3]
 [4]
```
"""
orbits(G::Group,v;action::F=^,trivial=true) where F<:Function=
           orbits(gens(G),v;action,trivial)

"""
`centralizer(G::Group,p;action=^)`

computes  the centralizer of  `p` in group  `G` (by default  for the action
`action(g,p)=g^p`)

```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3)]);
julia> centralizer(G,1)
Group([(2,3)])
```
"""
function centralizer(G::Group,p;action::Function=^)
# compute Schreier generators
  t=transversal(G,p;action=action)
  C=[wx*s/t[action(x,s)] for (x,wx) in t for s in gens(G)]
  Group(unique!(sort(C)))
end

centralizer(G::Group,H::Group)=centralizer(G,gens(H);action=(x,s)->x.^s)

"""
`stabilizer(G::Group,s)`

Assume that `s` is a set, represented as a sorted list without repetitions.
The  action  of  `g∈  G`  on  sets  is  given  by  `(g,p)->sort(p.^g)`. The
*stabilizer* of `s` in `G` is the centralizer of `s` for that action.

```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3,4)])
Group([(1,2),(1,2,3,4)])

julia> centralizer(G,[1,2];action=(s,g)->sort(s.^g))
Group([(3,4),(1,2),(1,2)(3,4)])

julia> stabilizer(G,[1,2])
Group([(3,4),(1,2),(1,2)(3,4)])
```
"""
stabilizer(G::Group,p)=centralizer(G,p;action=(s,g)->sort(s.^g))

function centre(G::Group)
  gets(G,:centre) do
    centralizer(G,G)
  end
end

"""
    `minimal_words(G)`
  returns a Dict giving for each element of `G` a minimal positive word in 
  the generators representing it.

```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3)]);
julia> minimal_words(G)
Dict{Perm{Int16},Array{Int64,1}} with 6 entries:
  ()      => Int64[]
  (1,3)   => [1, 2]
  (1,2,3) => [2]
  (1,2)   => [1]
  (2,3)   => [2, 1]
  (1,3,2) => [1, 2, 1]
```
"""
function minimal_words(G::Group)
  gets(G,:words)do
    words=Dict(one(G)=>Int[])
    for i in eachindex(gens(G))
      nwords=copy(words)
      rw = [one(G)=>Int[]]
      while !isempty(rw)          # look at group generated by gens(W)[1:i]
        p=popfirst!(rw)
        for k in 1:i
          e=first(p)*gens(G)[k]
          if !haskey(nwords,e)
            we=vcat(last(p),[k])
            push!(rw,e=>we)
            for (e1,w1) in words nwords[e1*e]=vcat(w1,we) end
          end
        end
      end
      words = nwords
    end
    words
  end
end

function minimal2(G::Group)
  gets(G,:words)do
    words=Dict(one(G)=>0)
    for i in eachindex(gens(G))
      nwords=copy(words)
      rw=[one(G)=>0]
      while !isempty(rw)          # look at group generated by gens(W)[1:i]
        p=popfirst!(rw)
        for k in 1:i
          e=first(p)*gens(G)[k]
          if !haskey(nwords,e)
            push!(rw,e=>k)
            for (e1,w1) in words nwords[e1*e]=k end
          end
        end
      end
      words = nwords
    end
    words
  end
end

"word(G::Group,w): a word in  gens(G) representing element w of G"
word(G::Group,w)=minimal_words(G)[w]

function word2(W::Group,w)
  res=Int[]
  d=minimal2(W)
  while !isone(w)
    i=d[w]
    w/=gens(W)[i]
    pushfirst!(res,i)
  end
  res
end

"elements(G::Group): the list of elements of G"
function elements(G::Group)
  collect(keys(minimal_words(G)))
end

Base.in(w,G::Group)=haskey(minimal_words(G),w)

"length(G::Group): the number of elements of G"
Base.length(G::Group)=length(minimal_words(G))

function conjugacy_classes(G::Group{T})::Vector{Vector{T}} where T
  gets(G,:classes) do
    if haskey(G.prop,:classreps)
      map(x->orbit(G,x),classreps(G))
    elseif length(G)>10000
      error("length(G)=",length(G),": should call Gap4")
    else
      orbits(G,elements(G))
    end
  end
end

function position_class(G::Group,g)
  findfirst(c->g in c,conjugacy_classes(G))
end

function fusion_conjugacy_classes(H::Group,G::Group)
  map(x->position_class(G,x),classreps(H))
end

"classreps(G::Group): representatives of conjugacy classes of G"
function classreps(G::Group{T})::Vector{T} where T
  gets(G,:classreps) do
    if length(G)>10000 Gap4.classreps(G)
    else first.(conjugacy_classes(G))
    end
  end
end

nconjugacy_classes(G::Group)=length(classreps(G))

# hom from source to target sending gens(source) to images
struct Hom{T,T1}
  source::Group{T}
  target::Group{T1}
  images::Vector{T1}
end

function Base.show(io::IO,h::Hom)
  print(io,"Hom(",h.source,"→ ",h.target,";",gens(h.source),"↦ ",h.images)
end

function kernel(h::Hom)
  if all(isone,h.images) return h.source
  elseif length(h.source)==length(Group(h.images)) 
    return Group(empty(gens(h.source)),one(h.source))
  elseif length(h.source)<1000
    return Group(filter(x->isone(h(x)),elements(h.source)))
  else error("not implemented: kernel(",h,")")
  end
end

# h(w) is the image of w by h
(h::Hom)(w)=isone(w) ? one(h.target) : prod(h.images[word(h.source,w)])

isabelian(W::Group)=all(x*y==y*x for x in gens(W), y in gens(W))

rand(W::Group)=W(rand(eachindex(gens(W)),20)...)

"""
`transporting_elt(G,p,q;action=^,dist=nothing)`   

returns  an  element  `g∈ G`  such  that  `p^g==q` (or `action(p,g)==q` if
`action`  is given) if such a `g`  exists and nothing otherwise. The set of
possible `g` forms a right coset of the centralizer of p in G.

```julia-repl
julia> g=Group(perm"(1,2,3)(6,7)",perm"(3,4,5)(7,8)")
Group([(1,2,3)(6,7),(3,4,5)(7,8)])

julia> transporting_elt(g,1,5)
(1,5,4,3,2)

julia> transporting_elt(g,1,6)

julia> transporting_elt(g,[1,2,3,4],[2,3,4,5];action=(s,g)->sort(s.^g))
(1,2,3,4,5)(6,7,8)

julia> transporting_elt(g,[1,2,3,4],[3,4,5,2];action=(s,g)->s.^g)
```
"""
function transporting_elt(W::Group,x,y;action::Function=^,dist=nothing)
  if isnothing(dist)
    t=transversal(W,x;action=action)
    if haskey(t,y) return t[y] else return nothing end
  end
# InfoChevie("#I  group:", length(W), " too big - trying random walk\n")
  p=one(W)
  InfoChevie(dist(x, y), " ")
  x1=x
  while true
    prev=dist(x1, y)
    if prev==0
      InfoChevie("\n")
      return p
    end
    dmin=minimum(map(g->(dist(action(x1, g), y),g),gens(W)))
    if dmin[1]<prev 
      InfoChevie("->",dmin)
      p*=dmin[2]
      x1=action(x1,dmin[2])
    else
      InfoChevie("\n[ Info: stalled -- restarting at a random element of",
                 "W(size $(length(W)))\n")
      p*=rand(W)
      x1=action(x,p)
    end
  end
end

# horrible implementation
function Base.intersect(G::Group, H::Group)
  if min(length(G),length(H))>104000 error("too large intersect($G,$H)") end
  if length(G)<length(H) Group(filter(x->x in H,elements(G)))
  else Group(filter(x->x in G,elements(H)))
  end
end

#------------------- "abstract" concrete groups -------------------------------
struct Groupof{T}<:Group{T}
  gens::Vector{T}
  one::T
  prop::Dict{Symbol,Any}
end

"""
`Group(l::AbstractVector{T}[,one]) where T`

A  group may be constructed  from a list of  `l` elements of the same type.
These  elements must respond  to the function  `*` and `inv`.  If it is not
possible to compute `one` from the elements (because they do not respond to
`one`,  or  `l`  is  empty  and  `T`  does  not respond to `one`), then the
identity element of the group must be given as a second argument.

```julia-repl
julia> G=Group([[-1 -1;1 0]])
Gapjm.Groups.Groupof{Array{Int64,2}}([[-1 -1; 1 0]], [1 0; 0 1], Dict{Symbol,Any}())

julia> elements(G)
3-element Array{Array{Int64,2},1}:
 [-1 -1; 1 0]
 [0 1; -1 -1]
 [1 0; 0 1]
```
"""
function Group(a::AbstractVector{T}) where T
  Groupof(filter(!isone,a),!isempty(a) ? one(a[1]) : one(T),Dict{Symbol,Any}())
end

# for the case one(a::T) is defined but not one(T)
function Group(a::AbstractVector{T},one) where T
  Groupof(filter(!=(one),a),one,Dict{Symbol,Any}())
end

Base.one(G::Groupof)=G.one

#------------------- cosets ----------------------------------------
# for now normal cosets (phi normalizes G)
abstract type Coset{TW<:Group} end

struct Cosetof{T,TW<:Group{T}}<:Coset{TW}
  phi::T
  G::TW
  prop::Dict{Symbol,Any}
end

Base.isone(a::Coset)=isone(a.phi)

Coset(G::Group,phi)=Cosetof(phi,G,Dict{Symbol,Any}())
Coset(G::Group)=Cosetof(one(G),G,Dict{Symbol,Any}())

Group(W::Coset)=W.G

(W::Coset)(x...)=Group(W)(x...)*W.phi

Base.cmp(a::Coset, b::Coset)=cmp(a.phi,b.phi)

Base.isless(a::Coset, b::Coset)=cmp(a,b)==-1

Base.:(==)(a::Coset, b::Coset)= (Group(a)==Group(b)) && cmp(a,b)==0

Base.hash(a::Coset, h::UInt)=hash(a.phi,h)

Base.copy(C::Coset)=Coset(Group(C),C.phi)

Base.one(C::Coset)=Coset(Group(C),one(C.phi))

Base.inv(C::Coset)=Coset(Group(C),inv(C.phi))

Base.:*(a::Coset,b::Coset)=Coset(Group(a),a.phi*b.phi)

Base.:/(a::Coset,b::Coset)=a*inv(b)

Base.length(C::Coset)=length(Group(C))

Base.:^(a::Coset, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                     Base.power_by_squaring(inv(a),-n)

Base.:^(a::Coset, b::Coset)= inv(b)*a*b

order(a::Coset)=findfirst(i->a.phi^i in Group(a),1:order(a.phi))

Base.show(io::IO,C::Coset)=print(io,Group(C),".",C.phi)

struct CosetGroup{TW}<:Group{Coset{TW}}
  gens::Vector{Coset{TW}}
  one::Coset{TW}
  prop::Dict{Symbol,Any}
end

Groups.Group(g::Vector{Coset},o=one(g[1]))=CosetGroup(filter(!isone,g),o,Dict{Symbol,Any}())

Base.one(G::CosetGroup)=G.one

Base.:/(W::Group,H::Group)=Group(unique(map(x->Coset(H,x),gens(W))),Coset(H))

elements(C::Coset)=C.phi .* elements(Group(C))

Base.in(w,C::Coset)=C.phi\w in Group(C)

function classreps(G::Coset{Group{T}})::Vector{T} where T
  gets(G,:classreps) do
    first.(conjugacy_classes(G))
  end
end

nconjugacy_classes(G::Coset)=length(classreps(G))

function conjugacy_classes(G::Coset{<:Group{T}})::Vector{Vector{T}} where T
  gets(G,:classes) do
    if length(G)>10000
      error("length(G)=",length(G),": should call Gap4")
    else
      orbits(Group(G),elements(G))
    end
  end
end

function position_class(G::Coset,g)
  findfirst(c->g in c,conjugacy_classes(G))
end

function fusion_conjugacy_classes(H::Coset,G::Coset)
  map(x->position_class(G,x),classreps(H))
end

#---------------------- miscellaneous functions ----------------------
function abelian_generators(l::Array)
  res=empty(l)
  l=filter(!isone,l)
  while !isempty(l)
    o=order.(l)
    push!(res,l[findfirst(isequal(maximum(o)),o)])
    l=setdiff(l,elements(Group(res)))
  end
  res
end

iscyclic(W::Group)=isabelian(W) && lcm(order.(gens(W)))==length(W)
end
