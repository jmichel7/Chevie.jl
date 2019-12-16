"""
This module is a port of some GAP functionality on permutation groups.

This code refers to Holt "Handbook of computational group theory" chapter 4
for basic algorithms.

A  PermGroup is  a group  where gens  are Perms,  which allows  for all the
algorithms like base, centralizer chain, etc...

# Examples
```julia-repl
julia> G=Group([Perm(i,i+1) for i in 1:2])
Group([(1,2),(2,3)])

# PermGroups are iterators over their elements
julia> collect(G)  
6-element Array{Perm{Int16},1}:
 (1,2)
 (1,3,2)
 ()
 (1,2,3)
 (1,3)
 (2,3)

# maximum degree of an element of G
julia> degree(G)  
3

julia> Perm(1,2) in G
true

julia> Perm(1,2,4) in G
false

# Elements,  appartenance test and  other function are  computed on G using
# Schreier-Sims theory, that is computing the following

# a list of points that no element of G fixes
julia> base(G) 
2-element Array{Int16,1}:
 1
 2

# the i-th element is the centralizer of base[1:i-1]
julia> centralizers(G) 
2-element Array{PermGroup{Int16},1}:
 Group([(1,2),(2,3)])
 Group([(2,3)])

# i-th element is transversal of centralizer[i] on base[i]
julia> transversals(G)
2-element Array{Dict{Int16,Perm{Int16}},1}:
 Dict(2 => (1,2),3 => (1,3,2),1 => ())
 Dict(2 => (),3 => (2,3))
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
export PermGroup, base, transversals, centralizers, symmetric_group, Coset

#-------------------- now permutation groups -------------------------
struct PermGroup{T}<:Group{Perm{T}}
  gens::Vector{Perm{T}}
  prop::Dict{Symbol,Any}
end

PermGroup(W::PermGroup)=W

function Groups.Group(a::AbstractVector{Perm{T}}) where T
  PermGroup(filter(!isone,a),Dict{Symbol,Any}())
end

function Base.show(io::IO,G::PermGroup)
  print(io,"Group([$(join(map(repr,gens(G)),','))])")
end

function Gapjm.degree(G::PermGroup)::Int
  gets(G,:degree)do G maximum(map(largest_moved_point,gens(G))) end
end

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
symmetric_group(n::Int)=Group([Perm(i,i+1) for i in 1:n-1])

"""
 The input is
 -  g: an element of a PermGroup G
 -  B: a base (or partial base) of G
 -  Δ: Δ[i] is the transversal of C_G(B[1:i-1]) on B[i]
 The function returns g "stripped" of its components in all C_G(B[1:i])
"""
function strip(g::Perm{T},B::Vector{T},Δ::Vector{Dict{T,Perm{T}}}) where T
  for i in eachindex(B)
    β=B[i]^g
    if !haskey(Δ[i],β) return g,i end
    g*=inv(Δ[i][β])
  end
  g,length(B)+1
end

"""
  see Holt, 4.4.2

  This function creates in G.prop the fields base, centralizers,
  transversals. See the description in the functions with the same name.
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
  H=[Group(s) for s in S]
  Δ=map(transversal,H,B)
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
            push!(H,Group(S[l]))
            push!(Δ,transversal(H[l],B[l]))
           else
           H[l]=Group(S[l])
           Δ[l]=transversal(H[l],B[l])
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
  G.prop[:transversals]=Δ
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
function transversals(G::PermGroup{T})::Vector{Dict{T,Perm{T}}} where T
  getp(schreier_sims,G,:transversals)
end

" A list of points stabilized by no element of G "
function base(G::PermGroup{T})::Vector{T} where T
  getp(schreier_sims,G,:base)
end

" Tells whether permutation g is an element of G "
function Base.in(g::Perm,G::PermGroup)
  g,i=strip(g,base(G),transversals(G))
  isone(g)
end
#-------------- iteration on product of lists of group elements ---------------
struct GroupProdIterator{T}
  iterators::Vector{T}
end
# if the iterators are i1,...,in iterate on all products i1[j1]*...*in[jn]

Base.length(I::GroupProdIterator)=prod(length.(I.iterators))

function Base.iterate(I::GroupProdIterator{T})where T
  prod=one(eltype(T))
  state=map(I.iterators) do l
    u=iterate(l)
    if isnothing(u) return nothing end
    prod*=first(u)
    prod,last(u)
  end
  prod::eltype(T),state::Vector{Tuple{eltype(T),Int}}
end

function Base.iterate(I::GroupProdIterator,state)
  for i in length(state):-1:1
    u=iterate(I.iterators[i],last(state[i]))
    if isnothing(u) continue end
    p,st=u
    state[i]= i==1 ? u : (first(state[i-1])*p,st)
    for j in i+1:length(state)
      p,st=iterate(I.iterators[j])
      state[j]=(first(state[j-1])*p,st)
    end
    return first(state[end]),state
  end
  return nothing
end

Base.eltype(::Type{GroupProdIterator{T}}) where T=eltype(T)
#------------------------- iteration for PermGroups -----------------------
function Base.length(G::PermGroup)::Int
  gets(G,:length)do G 
    prod(length.(transversals(G)))
  end
end

Base.eltype(::Type{PermGroup{T}}) where T=Perm{T}

function Base.iterate(G::PermGroup)
  I=GroupProdIterator(reverse(values.(transversals(G))))
  u=iterate(I)
  if isnothing(u) return u end
  p,st=u
  p,(I,st)
end

function Base.iterate(G::PermGroup,(I,state))
  u=iterate(I,state)
  if isnothing(u) return u end
  p,st=u
  p,(I,st)
end

# elements is much faster than collect(G), should not be
function Gapjm.elements(G::PermGroup)
  t=reverse(values.(transversals(G)))
  res=t[1]
  for i in 2:length(t)
    res=vcat(map(x->res.*x,t[i])...)
  end
  res
end
#------------------- cosets ----------------------------------------
struct Coset{T,TW<:Group{T}}
  w::T
  G::TW
end

# computes "canonical" element of W.w
function Coset(W::PermGroup,w)
  for i in eachindex(base(W))
    t=transversals(W)[i]
    (kw,e)=minimum((k^w,e) for (k,e) in t)
    w=e*w
  end
  Coset(w,W)
end

Base.cmp(a::Coset, b::Coset)=cmp(a.w,b.w)

Base.isless(a::Coset, b::Coset)=cmp(a,b)==-1

Base.:(==)(a::Coset, b::Coset)= cmp(a,b)==0

Base.hash(a::Coset, h::UInt)=hash(a.w,h)

Base.copy(C::Coset)=Coset(C.w,C.G)

Base.one(C::Coset)=Coset(one(C.w),C.G)

Base.inv(C::Coset)=Coset(C.G,inv(C.w))

Base.:*(a::Coset,b::Coset)=Coset(a.G,a.w*b.w)

Base.:^(a::Coset, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                               Base.power_by_squaring(inv(a),-n)

Perms.order(a::Coset)=findfirst(i->isone(a^i),1:order(a.w))

Base.show(io::IO,C::Coset)=print(io,C.G,".",C.w)

struct CosetGroup{T,TW}<:Group{Coset{T,TW}}
  gens::Vector{Coset{T,TW}}
  prop::Dict{Symbol,Any}
end

Groups.Group(g::Vector{Coset{T,TW}}) where {T,TW}=
  CosetGroup(filter(!isone,g),Dict{Symbol,Any}())

Base.:/(W::PermGroup,H::PermGroup)=Group(map(x->Coset(H,x),gens(W)))

end
