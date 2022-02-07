"""
This module gives some basic functionality on groups.

`Group`  is  an  abstract  type,  but  the  following is assumed of all its
concrete implementations:

  - The function `gens(G)` returns the list of generators of the group `G`. 
  - The function `one(G)` returns the identity element of `G`.
  
Further,  a group  is a  `GapObj`, a  kind of  object where  properties are
computed on demand when asked for; such properties when computed are stored
in the field `.prop` of `G`, which is of type `Dict{Symbol, Any}()`.

# Examples
```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3)])
Group([(1,2),(1,2,3)])

julia> gens(G)
2-element Vector{Perm{Int16}}:
 (1,2)  
 (1,2,3)

julia> ngens(G)
2

julia> minimal_words(G)
Dict{Perm{Int16}, Vector{Int64}} with 6 entries:
  ()      => []
  (1,2)   => [1]
  (1,3)   => [1, 2]
  (1,2,3) => [2]
  (2,3)   => [2, 1]
  (1,3,2) => [2, 2]
  

julia> G.prop
Dict{Symbol, Any} with 1 entry:
  :minwords => Dict(()=>[], (1,2)=>[1], (1,3)=>[1, 2], (1,2,3)=>[2], (2,3)=>[2,…
  
```
"""
module Groups
export Group, centralizer, centre, classreps, comm,
  conjugacy_class, conjugacy_classes, Coset, fusion_conjugacy_classes,
  gens, Hom, isabelian, iscyclic, words, minimal_words, nconjugacy_classes,
  ngens, normalizer, orbit, orbits, position_class, stabilizer,
  transporting_elt, transversal

using ..Util: InfoChevie, @GapObj
#import Gapjm: word, elements, kernel, order
# to use as a stand-alone module comment above line and uncomment next line
export word, elements, kernel, order
#--------------general groups and functions for "black box groups" -------
abstract type Group{T} end # T is the type of elements of G

function normalizer end

"`one(G::Group)` returns the identity element of `G`."
Base.one(G::Group{T}) where T=one(T)

"`gens(G::Group)` returns the `Vector` of generators of `G`."
gens(G::Group)=G.gens
"`ngens(G::Group)` returns the number of generators of `G`."
ngens(G::Group)=length(gens(G))

"""
`(G::Group)(i...)`

A Group used as a function takes integer arguments in `eachindex(gens(W))`.
This  constructs  the  element  of  `G`  product of the generators with the
specified  indices. An argument  can also be  negative, then the inverse of
the corresponding generator is used.

```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3)])
Group([(1,2),(1,2,3)])

julia> G(2,1,-2) # returns gens(G)[2]*gens(G)[1]/gens(G)[2]
(1,3)
```
"""
function (W::Group)(w::Vararg{Integer,N}where N)
  isempty(w) ? one(W) : prod(W(i) for i in w)
end

(W::Group)(i::Integer)=i>0 ? gens(W)[i] : inv(gens(W)[-i])

"`comm(a,b)` the commutator of `a` and `b`"
comm(a,b)=inv(a)*inv(b)*a*b

"""
`orbit(gens::AbstractVector,p;action::Function=^)`

`orbit(G::Group,p;action::Function=^)`

the  orbit of point  `p` under repeated  action of generators `gens`. Point
`p`  should be hashable. The default action  of a group element is `^`. For
example  if `g` is a permutation and `p`  an integer, `p^g` is the image of
`p`  by `g`; if `h` and `g` are group elements, then `h^g` is the conjugate
`inv(g)*h*g`.  If a group  is given instead  of generators, the orbit under
`gens(G)` is returned.

```julia-repl
julia> orbit([Perm(1,2),Perm(2,3)],1) 
3-element Vector{Int64}:
 1
 2
 3

julia> orbit([Perm(1,2),Perm(2,3)],[1,3];action=(v,g)->v.^g) # Gap "OnTuples"
6-element Vector{Vector{Int64}}:
 [1, 3]
 [2, 3]
 [1, 2]
 [3, 2]
 [2, 1]
 [3, 1]

julia> orbit([Perm(1,2),Perm(2,3)],[1,3];action=(v,g)->sort(v.^g)) # "OnSets"
3-element Vector{Vector{Int64}}:
 [1, 3]
 [2, 3]
 [1, 2]
```
"""
function orbit(gens::AbstractVector,pnt;action::F=^) where F<:Function
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

returns a `Dict` with entries `x=>g` where `x` runs over
`orbit(G,p;action)` and where `g∈ G` is such that `x=action(p,g)`

```julia-repl
julia> G=Group([Perm(1,2),Perm(2,3)]);
julia> transversal(G,1)
Dict{Int64, Perm{Int16}} with 3 entries:
  2 => (1,2)
  3 => (1,3,2)
  1 => ()
```

orbit functions can take any action of `G` as keyword argument

```julia-repl
julia> transversal(G,[1,2],action=(x,y)->x.^y)
Dict{Vector{Int64}, Perm{Int16}} with 6 entries:
  [2, 3] => (1,2,3)
  [2, 1] => (1,2)
  [1, 3] => (2,3)
  [3, 1] => (1,3,2)
  [1, 2] => ()
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

function orbits(gens::AbstractVector,v::AbstractVector;action=^,trivial=true)
  res=Vector{eltype(v)}[]
  while !isempty(v)
    o=orbit(gens,first(v);action)
    if length(o)>1 || trivial push!(res,o) end
    v=setdiff(v,o)
  end
  res
end

"""
`orbits(gens::Vector,v;action=^,trivial=true)`

`orbits(G,v;action=^,trivial=true)`
    
the  orbits on `v`  of the repeated  action of `gens`;  the elements of `v`
should  be hashable. If a  group is given instead  of generators, the orbit
under  `gens(G)` is returned. If `trivial=false` the one-element orbits are
not returned.

```julia-repl
julia> G=Group([Perm(1,2),Perm(2,3)]);
julia> orbits(G,1:4)
2-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [4]
```
"""
orbits(G::Group,v;action=^,trivial=true)=orbits(gens(G),v;action,trivial)

"""
`centralizer(G::Group,p;action=^)`

computes  the centralizer of `p` in `G`.

```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3)]);
julia> centralizer(G,1)
Group([(2,3)])
```
"""
function centralizer(G::Group,p;action::Function=^)
  t=transversal(G,p;action)
  if length(t)==1 return G end
  C=[wx*s/t[action(x,s)] for (x,wx) in t for s in gens(G)] #Schreier generators
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
  get!(G,:centre) do
    centralizer(G,G)
  end
end

"""
`minimal_words(G)`

returns  a `Dict` giving for each element of `G` a minimal positive word in
the generators representing it.

```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3)]);
julia> minimal_words(G)
Dict{Perm{Int16}, Vector{Int64}} with 6 entries:
  ()      => []
  (1,2)   => [1]
  (1,3)   => [1, 2]
  (1,2,3) => [2]
  (2,3)   => [2, 1]
  (1,3,2) => [2, 2]
```
"""
function minimal_words(G::Group)
  get!(G,:minwords)do
    words=Dict(one(G)=>Int[])
    l=words
    while !isempty(l)
      newl=Dict{eltype(G),Vector{Int}}()
      for (i,g) in enumerate(gens(G)), (h,w) in l
        e=h*g
        if !haskey(words,e) newl[e]=vcat(w,i) end
      end
      l=newl
      merge!(words,l)
    end
    words
  end
end

"""
`words(G)`

returns  a `Dict`  giving for  each element  of `G`  a positive word in the
generators representing it. It is faster than `minimal_words` but the words
are not guaranteed minimal.

```julia-repl
julia> G=Group([Perm(1,2),Perm(1,2,3)]);
julia> words(G)
Dict{Perm{Int16}, Vector{Int64}} with 6 entries:
  ()      => []
  (1,2)   => [1]
  (1,3)   => [1, 2]
  (1,2,3) => [2]
  (2,3)   => [2, 1]
  (1,3,2) => [1, 2, 1]
```
"""
function words(G::Group)
  get!(G,:words)do
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

# returns Dict: (word w, n0 generator such that l(w/G(n0))<l(w)
function words2(G::Group) # faster than words but longer to retrieve word
  get!(G,:words2)do
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

"word(G::Group,w): a minimal word in  gens(G) representing element w of G"
word(G::Group,w)=minimal_words(G)[w]

function word2(W::Group,w)
  res=Int[]
  d=words2(W)
  while !isone(w)
    i=d[w]
    w/=gens(W)[i]
    pushfirst!(res,i)
  end
  res
end

"elements(G::Group): the list of elements of G"
function elements(G::Group)
  collect(keys(words2(G)))
end

Base.in(w,G::Group)=haskey(words2(G),w)

"length(G::Group): the number of elements of G"
Base.length(G::Group)=length(words2(G))

"`conjugacy_classes(G::Group)` conjugacy classes of `G` (as a `Vector{Vector}`)"
function conjugacy_classes(G::Group{T})::Vector{Vector{T}} where T
  get!(G,:classes) do
    if haskey(G,:classreps) 
      [conjugacy_class(G,i) for i in 1:nconjugacy_classes(G)]
    else orbits(G,elements(G))
    end
  end
end

"`conjugacy_classes(G::Group,i::Integer)` the `i`-th class of `G`"
function conjugacy_class(G::Group{T},i::Int)::Vector{T} where T
  if !haskey(G,:classes) 
    G.classes=Vector{Vector{T}}(undef,nconjugacy_classes(G))
  end
  if !isassigned(G.classes,i) G.classes[i]=orbit(G,classreps(G)[i]) end
  G.classes[i]
end

"`position_class(G::Group,g)` index of conjugacy class to which `g` belongs"
function position_class(G::Group,g)
  findfirst(c->g in c,conjugacy_classes(G))
end

"""
'fusion_conjugacy_classes(H::Group,G::Group)`

A `Vector{Int}` telling for each conjugacy class of subgroup `H` of which class
of `G` is is a subset
"""
function fusion_conjugacy_classes(H::Group,G::Group)
  map(x->position_class(G,x),classreps(H))
end

"""
`classreps(G::Group)` 

representatives of conjugacy classes of `G`. Fills `G.classreps`. If this
field is filled it is used by  `conjugacy_classes`.
"""
function classreps(G::Group{T})::Vector{T} where T
  get!(G,:classreps) do
    if length(G)>10000 error("length(G)=",length(G),": should call Gap4")
    else minimum.(conjugacy_classes(G))
    end
  end
end

"`nconjugacy_classes(G::Group)` the number of conjugacy classes of `G`"
nconjugacy_classes(G::Group)=length(classreps(G))

# hom from source to target sending gens(source) to images
struct Hom{T,T1}
  source::Group{T}
  target::Group{T1}
  images::Vector{T1}
end

function Base.show(io::IO,h::Hom)
  if h.source==h.target
    print(io,"Aut(",h.source,";",gens(h.source),"↦ ",h.images)
  else
    print(io,"Hom(",h.source,"→ ",h.target,";",gens(h.source),"↦ ",h.images)
  end
end

"`kernel(h::Hom)` the kernel of the homomorphism `h`"
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
(h::Hom)(w)=isone(w) ? one(h.target) : prod(
  (i>0 ? h.images[i] : inv(h.images[-i])) for i in word(h.source,w))

"`isabelian(G::Group)` whether `G` is abelian"
isabelian(W::Group)=all(x*y==y*x for x in gens(W), y in gens(W))

"`iscyclic(G::Group)` whether `G` is cyclic"
iscyclic(W::Group)=isabelian(W) && lcm(order.(gens(W)))==length(W)

Base.rand(W::Group)=W(rand(eachindex(gens(W)),20)...)

"""
`transporting_elt(G,p,q;action=^,dist=nothing)`   

returns  an  element  `g∈  G`  such  that  `p^g==q` (or `action(p,g)==q` if
`action` is given), if such a `g` exists, and nothing otherwise. The set of
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
    if x==y return one(W) end
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

# suppress unneeded generators
function weedgens(G::Group)
  l=length(G)
  gen=empty(gens(G))
  if l==1 return Group(gen) end
  lr=1
  for g in gens(G)
    N=Group(vcat(gen,[g]))
    ln=length(N)
    if ln==l return N end
    if ln>lr
      lr=ln
      gen=vcat(gen,[g])
    end
  end
end

# horrible implementation
function Base.intersect(G::Group, H::Group)
  if min(length(G),length(H))>104000 error("too large intersect($G,$H)") end
  if length(G)<length(H) res=Group(filter(x->x in H,elements(G)))
  else res=Group(filter(x->x in G,elements(H)))
  end
  weedgens(res)
end

#------------------- "abstract" concrete groups -------------------------------
@GapObj struct Groupof{T}<:Group{T}
  gens::Vector{T}
  one::T
end

"""
`Group(l::AbstractVector{T}[,one]) where T`

A  group may be constructed  from a list of  `l` elements of the same type.
These  elements must respond to  the functions `*` and  `inv`. If it is not
possible  to compute  `one` from  `l` (because  ``l[1]` does not respond to
`one`,  or  `l`  is  empty  and  `T`  does  not respond to `one`), then the
identity element of the group must be given as a second argument.

```julia-repl
julia> G=Group([[-1 -1;1 0]])
Gapjm.Groups.Groupof{Matrix{Int64}}([[-1 -1; 1 0]], [1 0; 0 1], Dict{Symbol, Any}())

julia> elements(G)
3-element Vector{Matrix{Int64}}:
 [0 1; -1 -1]
 [1 0; 0 1]
 [-1 -1; 1 0]
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
# for now only normal cosets (phi normalizes G)
abstract type Coset{TW<:Group} end

@GapObj struct Cosetof{T,TW<:Group{T}}<:Coset{TW}
  phi::T
  G::TW
end

Base.isone(a::Coset)=a.phi in Group(a)

Coset(G::Group,phi=one(G))=Cosetof(phi,G,Dict{Symbol,Any}())

Group(W::Coset)=W.G

(W::Coset)(x...)=Group(W)(x...)*W.phi

Base.cmp(a::Coset, b::Coset)=cmp(a.phi,b.phi)

Base.isless(a::Coset, b::Coset)=cmp(a,b)==-1

Base.:(==)(a::Coset, b::Coset)= (Group(a)==Group(b)) && cmp(a,b)==0

Base.hash(a::Coset, h::UInt)=hash(a.phi,h)

Base.copy(C::Coset)=Coset(Group(C),C.phi)

Base.one(C::Coset)=Coset(Group(C))

Base.inv(C::Coset)=Coset(Group(C),inv(C.phi))

Base.:*(a::Coset,b::Coset)=Coset(Group(a),a.phi*b.phi)

Base.:/(a::Coset,b::Coset)=a*inv(b)

Base.length(C::Coset)=length(Group(C))

Base.:^(a::Coset, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                     Base.power_by_squaring(inv(a),-n)

Base.:^(a::Coset, b::Coset)= inv(b)*a*b

order(a::Coset)=findfirst(i->a.phi^i in Group(a),1:order(a.phi))

Base.show(io::IO,C::Coset)=print(io,Group(C),".",C.phi)

@GapObj struct CosetGroup{TW}<:Group{Coset{TW}}
  gens::Vector{Coset{TW}}
  one::Coset{TW}
end

Groups.Group(g::Vector{Coset},o=one(g[1]))=CosetGroup(filter(!isone,g),o,Dict{Symbol,Any}())

Base.one(G::CosetGroup)=G.one

Base.:/(W::Group,H::Group)=Group(unique(map(x->Coset(H,x),gens(W))),Coset(H))

elements(C::Coset)=C.phi .* elements(Group(C))

Base.in(w,C::Coset)=C.phi\w in Group(C)

function classreps(G::Coset{Group{T}})::Vector{T} where T
  get!(G,:classreps) do
    first.(conjugacy_classes(G))
  end
end

nconjugacy_classes(G::Coset)=length(classreps(G))

function conjugacy_classes(G::Coset{<:Group{T}})::Vector{Vector{T}} where T
  get!(G,:classes) do
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
end
