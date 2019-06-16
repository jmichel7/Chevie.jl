"""
This module is a port of some GAP functionality on groups, in particular
permutation groups.

This  codes refers to Holt "Handbook of computational group theory" chapter
4 for basic algorithms.

The only field of a Group G at the start is gens, the list of generators of
G.  To  mimic  GAP  records  where  attributes/properties  of an object are
computed  on demand when asked for, other fields are computed on demand and
stored in the field prop of the Group, which starts as Dict{Symbol,Any}()

A  PermGroup is  a group  where gens  are Perms,  which allows  for all the
algorithms like base, centralizer chain, etc...

# Examples
```julia-repl
julia> G=PermGroup([Perm(i,i+1) for i in 1:2])
PermGroup((1,2),(2,3))

# PermGroups are iterators over their elements
julia> collect(G)  
6-element Array{Perm{Int64},1}:
 (1,2)
 (1,3,2)
 ()
 (1,2,3)
 (1,3)
 (2,3)

# maximum degree of an element of G
julia> degree(G)  
3

# orbit of point 1 under G
julia> orbit(G,1) 
3-element Array{Int64,1}:
 2
 3
 1

# orbit decorated with representatives moving 1 to given point
julia> orbit_and_representative(G,1)
Dict{Int64,Perm{Int64}} with 3 entries:
  2 => (1,2)
  3 => (1,3,2)
  1 => ()

# orbit functions can take any action of G as keyword argument
julia> orbit_and_representative(G,[1,2],action=(x,y)->x.^Ref(y))
Dict{Array{Int64,1},Perm{Int64}} with 6 entries:
  [1, 3] => (2,3)
  [1, 2] => ()
  [2, 3] => (1,2,3)
  [3, 2] => (1,3)
  [2, 1] => (1,2)
  [3, 1] => (1,3,2)

julia> Perm(1,2) in G
true

julia> Perm(1,2,4) in G
false

# Elements,  appartenance test and  other function are  computed on G using
# Schreier-Sims theory, that is computing the following

# a list of points that no element of G fixes
julia> base(G) 
2-element Array{Int64,1}:
 1
 2

# the i-th element is the centralizer of base[1:i-1]
julia> centralizers(G) 
2-element Array{PermGroup{Int64},1}:
 PermGroup((1,2),(2,3))
 PermGroup((2,3))

# i-th element is orbit_and_representative of centralizer[i] on base[i]
julia> centralizer_orbits(G)
2-element Array{Dict{Int64,Perm{Int64}},1}:
 Dict(2=>(1,2),3=>(1,3,2),1=>())
 Dict(2=>(),3=>(2,3))

julia> minimal_words(G)  # minimal word in gens for each element of G
Dict{Perm{Int64},Array{Int64,1}} with 6 entries:
  ()      => Int64[]
  (2,3)   => [2]
  (1,3,2) => [1, 2]
  (1,3)   => [1, 2, 1]
  (1,2)   => [1]
  (1,2,3) => [2, 1]
```

finally, benchmarks on julia 1.0.1
```benchmark
julia> @btime length(collect(symmetric_group(8)))
  5.481 ms (270429 allocations: 12.40 MiB)

julia> @btime minimal_words(symmetric_group(8));
  10.477 ms (122062 allocations: 15.22 MiB)
```

Compare to GAP3 Elements(SymmetricGroup(8)); takes 3.8 ms
"""
module PermGroups
using ..Perms
using ..Gapjm # for degree, gens, minimal_words
export Group, PermGroup, orbit, orbit_and_representative, orbits,
  base, centralizer_orbits, centralizers, minimal_words, element,
  symmetric_group, gens, nbgens, centralizer, CharTable, class_reps

#--------------general groups and functions for "black box groups" -------
abstract type Group{T} end # T is the type of elements of G

Base.one(G::Group{T}) where T=one(T)
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
function orbit_and_representative(G::Group,p;action::Function=^)
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

function centralizer(G::Group,p;action::Function=^)
  new=[p]
  res=Dict(p=>one(G))
  C=PermGroup(empty(gens(G)))
  while !isempty(new)
    old=copy(new)
    empty!(new)
    for p in old, s in gens(G)
      w=action(p,s)
      if !haskey(res,w)
        push!(new,w)
        res[w]=res[p]*s
      else
        extend!(C,res[p]*s/res[w])
      end
    end
  end
  C
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

"""
assume l is union of orbits under group elt g; return permutation of l by g
needs objects in l sortable
"""
Perm{T}(g,l::AbstractVector;action::Function=^) where T<:Integer=Perm{T}(l,action.(l,Ref(g)))

Perm(g,l::AbstractVector;action::Function=^)=Perm{Int}(g,l,action=action)

" dict giving for each element of G a minimal word in the generators"
function minimal_words(G::Group)
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

function conjugacy_classes(G::Group{T})::Vector{Vector{T}} where T
  gets(G,:classes) do G
    if length(G)>1000
      println("length(G)=",length(G),": not supposed to do it the hard way")
    end
    cl=orbits(G,collect(G))
    G.prop[:classreps]=first.(cl)
    G.prop[:classes]=cl
  end
end

function class_reps(G::Group{T})::Vector{T} where T
  getp(conjugacy_classes,G,:classreps)
end

#--------------- CharTables -----------------------------------------
struct CharTable{T}
  irr::Matrix{T}
  charnames::Vector{String}
  classnames::Vector{String}
  centralizers::Vector{Int}
  identifier::String
end

function Base.show(io::IO,ct::CharTable)
  println(io,"CharTable(",ct.identifier,")")
  irr=map(ct.irr)do e
   if iszero(e) "." else sprint(show,e; context=io) end
  end
  format(io,irr,row_labels=TeXstrip.(ct.charnames),
                column_labels=TeXstrip.(ct.classnames))
end

#-------------------- now permutation groups -------------------------
struct PermGroup{T}<:Group{Perm{T}}
  gens::Vector{Perm{T}}
  prop::Dict{Symbol,Any}
end

function PermGroup(a::AbstractVector{Perm{T}})where T
  a=filter(x->!isone(x),a)
  PermGroup(a,Dict{Symbol,Any}())
end

function Base.show(io::IO,G::PermGroup)
  print(io,"PermGroup($(join(map(repr,gens(G)),',')))")
end

function Gapjm.degree(G::PermGroup)::Int
  gets(G,:degree)do G maximum(map(largest_moved_point,gens(G))) end
end

(W::PermGroup)(x...)=element(W,x...)

" describe the orbit of Int p under PermGroup G as a Schreier vector "
function schreier_vector(G::PermGroup,p::Integer;action::Function=^)
  res=zeros(Int,degree(G))
  res[p]=-1
  new=BitSet([p])
  while true
    n=new
    new=BitSet([])
    for p in n, i in eachindex(gens(G))
      q=action(p,gens(G)[i])
      if res[q]==0
        res[q]=i
        push!(new,q)
      end
    end
    if isempty(new) break end
  end
  res
end

" The symmetric group of degree n "
function symmetric_group(n::Int)
  PermGroup([Perm(i,i+1) for i in 1:n-1])
end

"""
 The input is
 -  g: an element of a PermGroup G
 -  B: a base (or partial base) of G
 -  Δ: Δ[i] is the orbit_and_representative of C_G(B[1:i-1]) on B[i]
 The function returns g "stripped" of its components in all C_G(B[1:i])
"""
function strip(g::Perm{T},B::Vector{T},Δ::Vector{Dict{T,Perm{T}}}) where T
  h=g
  for i in eachindex(B)
    β=B[i]^h
    if !haskey(Δ[i],β)
      return h,i
    end
    h*=inv(Δ[i][β])
  end
  h,length(B)+1
end

"""
  see Holt, 4.4.2

  This function creates in G.prop the fields base, centralizers,
  centralizer_orbits. See the description in the functions with the same name.
"""
function schreier_sims(G::PermGroup{T})where T
  B=T[]
  S=Vector{Perm{T}}[]
  for x in gens(G)
    j=1
    while j<=length(B)
      push!(S[j],x)
      if B[j]^x!=B[j] break end
      j+=1
    end
    if j>length(B)
      push!(B,smallest_moved_point(x))
      push!(S,[x])
    end
  end
  H=[PermGroup(s) for s in S]
  Δ=map(orbit_and_representative,H,B)
  rep(v)=join(map(repr,v),',')
  i=length(B)
  while i>=1
   for β in keys(Δ[i]), x in S[i]
     h=Δ[i][β]* x *inv(Δ[i][β^x])
     if !isone(h)
       y=true
       h,j=strip(h,B,Δ)
       if j<=length(B)
         y=false
       elseif !isone(h)
         y=false
         push!(B,smallest_moved_point(h))
         push!(S,Perm{T}[])
       end
       if y==false
         for l in i+1:j
           push!(S[l],h)
           if l>length(H)
            push!(H,PermGroup(S[l]))
            push!(Δ,orbit_and_representative(H[l],B[l]))
           else
           H[l]=PermGroup(S[l])
           Δ[l]=orbit_and_representative(H[l],B[l])
           end
         end
         i=j
         @goto nexti
       end
     end
   end
   i-=1
   @label nexti
  end
  G.prop[:base]=B
  G.prop[:centralizers]=H
  G.prop[:centralizer_orbits]=Δ
end

" centralizers: the i-th element is the centralizer of base[1:i-1]"
function centralizers(G::PermGroup{T})::Vector{PermGroup{T}} where T
  getp(schreier_sims,G,:centralizers)
end

"""
  The  i-th element  is  a description of  the orbit of :centralizers[i] on
  :base[i]  as a Dict where each point q is the key to a permutation p such
  that :base[i]^p=q
"""
function centralizer_orbits(G::PermGroup{T})::Vector{Dict{T,Perm{T}}} where T
  getp(schreier_sims,G,:centralizer_orbits)
end

" A list of points stabilized by no element of G "
function base(G::PermGroup{T})::Vector{T} where T
  getp(schreier_sims,G,:base)
end

" length(G::PermGroup) returns the cardinality of G "
function Base.length(G::PermGroup)::Int
  gets(G,:length)do G 
    prod(map(length,centralizer_orbits(G)))
  end
end

" Tells whether permutation g is an element of G "
function Base.in(g::Perm,G::PermGroup)
  g,i=strip(g,base(G),centralizer_orbits(G))
  isone(g)
end

" extend G by adding generator s to G"
function extend!(G::PermGroup,s::Perm)
  if !(s in G)
    push!(G.gens,s)
    schreier_sims(G)
  end
  G
end

# if l1,...,ln are the centralizer orbits the elements are the products
# of one element in each li
function Base.iterate(G::PermGroup)
  prod=one(G)
  state=map(reverse(values.(centralizer_orbits(G)))) do l
    u=iterate(l)
    if isnothing(u) return nothing end
    prod*=u[1]
    prod,u[2]
  end
  prod,reverse(state)
end

function Base.iterate(G::PermGroup,state)
  for i in eachindex(state)
    u=iterate(values(centralizer_orbits(G)[i]),state[i][2])
    if isnothing(u) continue end
    if i==length(state)
      state[i]=u
    else
      state[i]=(state[i+1][1]*u[1],u[2])
    end
    for j in i-1:-1:1
      u=iterate(values(centralizer_orbits(G)[j]))
      state[j]=(state[j+1][1]*u[1],u[2])
    end
    return state[1][1],state
  end
  return nothing
end

Base.eltype(::Type{PermGroup{T}}) where T=Perm{T}

end
