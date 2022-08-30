"""
This module is a port of some GAP functionality on permutation groups.

The code refers to the [Handbook of computational group theory, chapter 4]
(biblio.htm#H05) for basic algorithms.

A  PermGroup is  a group  where gens  are Perms,  which allows  for all the
algorithms like base, centralizer chain, etc...

# Examples
```julia-repl
julia> G=Group([Perm(i,i+1) for i in 1:2])
Group([(1,2), (2,3)])

# PermGroups are iterators over their elements
julia> collect(G)  
6-element Vector{Perm{Int16}}:
 (1,2)
 (1,3,2)
 ()
 (1,2,3)
 (1,3)
 (2,3)

# maximum moved point of an element of G
julia> largest_moved_point(G)  
3

julia> Perm(1,2) in G
true

julia> Perm(1,2,4) in G
false

# Elements,  appartenance test and  other function are  computed on G using
# Schreier-Sims theory, that is computing the following

# a list of points that no element of G fixes
julia> base(G) 
2-element Vector{Int16}:
 1
 2

# the i-th element is the centralizer of base[1:i-1]
julia> centralizers(G) 
2-element Vector{PermGroup{Int16}}:
 Group([(1,2), (2,3)])
 Group([(2,3)])

# i-th element is transversal of centralizer[i] on base[i]
julia> transversals(G)
2-element Vector{Dict{Int16, Perm{Int16}}}:
 Dict(2 => (1,2), 3 => (1,3,2), 1 => ())
 Dict(2 => (), 3 => (2,3))
```

finally, benchmarks on julia 1.7
```benchmark
julia> @btime collect(symmetric_group(8));
  3.720 ms (271800 allocations: 7.94 MiB)

julia> @btime words(symmetric_group(8));
  8.730 ms (202029 allocations: 12.86 MiB)
  
julia> @btime elements(symmetric_group(8));
  1.616 ms (52681 allocations: 3.77 MiB)
```
Compare to GAP3 Elements(SymmetricGroup(8)); takes 8 ms (GAP4 9 ms)
"""
module PermGroups
#using ..Groups
#using UsingMerge
#@usingmerge Perms # will do that when Perms is a package
using ..Gapjm
using ..Util: getp, InfoChevie, @GapObj, printTeX, hasdecor
using ..Combinat: tally, collectby
export PermGroup, base, transversals, centralizers, symmetric_group, reduced,
  stab_onmats, Perm_onmats, onmats, on_classes
#-------------------- now permutation groups -------------------------
abstract type PermGroup{T}<:Group{Perm{T}} end

PermGroup()=Group(Perm{Int16}[])

Base.show(io::IO,G::PermGroup)=print(io,"Group(",gens(G),")")

Base.one(G::PermGroup)=G.one # PermGroups should have fields gens and one

"`largest_moved_point(G::PermGroup)` the largest moved point by any `g∈ G`"
function Perms.largest_moved_point(G::PermGroup)::Int
  get!(G,:largest_moved)do 
    if isempty(gens(G)) return 0 end
    maximum(largest_moved_point.(gens(G)))
  end
end

function Perms.orbits(G::PermGroup)
  get!(G,:orbits)do
    orbits(G,1:largest_moved_point(G);trivial=false)
  end
end
 
"""
describe the orbit of Int p under PermGroup G as a Schreier vector v.
That is, v[p]==-1 and v[k]=i means that k^inv(G(i)) is the antecessor of k
in the orbit of p.
"""
function schreier_vector(G::PermGroup,p::Integer;action::Function=^)
  res=zeros(Int,largest_moved_point(G))
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

"""
 The input is
 -  g: a permutation
 -  B: a base (or partial base) of a PermGroup G
 -  Δ: Δ[i]==transversal(C_G(B[1:i-1]),B[i])
 The function returns g "stripped" of its components in all C_G(B[1:i]),
 that is a pair (an element which fixes B[1:i],i+1)
"""
function strip(g::Perm,B::Vector{<:Integer},Δ::Vector{Dict{T,Perm{T}}}) where T
  for i in eachindex(B)
    β=B[i]^g
    if !haskey(Δ[i],β) return g,i end
    g/=Δ[i][β]
  end
  g,length(B)+1
end

"""
  see Holt, 4.4.2

  This function creates in G.prop the fields base, centralizers,
  transversals. See the description in the functions with the same name.
"""
function schreier_sims(G::PermGroup{T})where T
  B=T[]  # base
  C=PG{T}[]  # C[i] will become C_G(B[1:i-1])
  for x in gens(G)
    j=1
    while j<=length(B)
      push!(gens(C[j]),x)
      if B[j]^x!=B[j] break end
      j+=1
    end
    if j>length(B)
      push!(B,smallest_moved_point(x))
      push!(C,Group([x]))
    end
  end
  Δ=transversal.(C,B)
  i=length(B)
  while i>=1
    for (β,wβ) in Δ[i], x in gens(C[i])
      h=wβ*x/Δ[i][β^x] # possibly new elt of C_G(B[1:i])
      if isone(h) continue end
      h,j=strip(h,B,Δ)
      if isone(h) continue end
      for l in i+1:j # now h is in C[l] for those l
        if l>length(C)
          push!(B,smallest_moved_point(h))
          push!(C,Group([h]))
          push!(Δ,transversal(C[l],B[l]))
        else
          push!(gens(C[l]),h)
          Δ[l]=transversal(C[l],B[l])
        end
      end
      i=j
      @goto nexti
    end
    i-=1
    @label nexti
  end
  G.base=B
  G.centralizers=C
  G.transversals=Δ
end

" centralizers: the i-th element is the centralizer of base[1:i-1]"
function centralizers(G::PermGroup{T})::Vector{PermGroup{T}} where T
  getp(schreier_sims,G,:centralizers)
end

"""
`transversals(G::PermGroup)`

returns a list whose `i`-th element is the transversal of
`G.centralizers[i]` on `G.base[i]`
"""
function transversals(G::PermGroup{T})::Vector{Dict{T,Perm{T}}} where T
  getp(schreier_sims,G,:transversals)
end

" `base(G::PermGroup)` A `Vector` of points stabilized by no element of `G` "
function base(G::PermGroup{T})::Vector{T} where T
  getp(schreier_sims,G,:base)
end

" `g in G` tells whether the permutation `g` is an element of `G` "
function Base.in(g::Perm,G::PermGroup)
  g,i=strip(g,base(G),transversals(G))
  isone(g)
end

function Base.intersect(G::PermGroup, H::PermGroup)
  if all(x->x in H,gens(G)) return G end
  if all(x->x in G,gens(H)) return H end
  if min(length(G),length(H))>104000 
    println("*** too large intersect($G,$H) -- calling Gap4.intersect") 
    return Gapjm.Gap4.intersect(G,H)
  end
  if length(G)<length(H) res=Group(filter(x->x in H,elements(G)))
  else res=Group(filter(x->x in G,elements(H)))
  end
  Groups.weedgens(res)
end

cycletypes(W::PermGroup,x)=map(o->cycletype(x,domain=o),orbits(W)) # first invariant

# eventually rewrite this dispatching on types
function classinv(W::PermGroup)
  get!(W,:classinv)do
    cycletypes.(Ref(W),classreps(W))
  end
end

# internal function accepting ambiguity
function positions_class(W::PermGroup,w)
  l=findall(==(cycletypes(W,w)),classinv(W))
  if length(l)==1 return l end
  if length(centre(W))>1
    central=gens(centre(W))
    l=filter(i->cycletypes.(Ref(W),classreps(W)[i].*central)==
                cycletypes.(Ref(W),w.*central),l)
    if length(l)==1 return l end
  end
  l
end

function Groups.position_class(W::PermGroup,w)
  l=positions_class(W,w)
  if length(l)==1 return only(l) end
  for i in eachindex(l) 
    if w in conjugacy_classes(W,l[i]) return l[i] end
  end
end

#-------------- iteration on product of lists of group elements ---------------
# if iterators=[i1,...,in] iterate on all products i1[j1]*...*in[jn]
struct ProdIterator{T}
  iterators::Vector{T}
end

Base.length(I::ProdIterator)=prod(length.(I.iterators))
Base.eltype(::ProdIterator{T}) where T=eltype(T)

function Base.iterate(I::ProdIterator)
  T=eltype(eltype(I.iterators))
  state=Tuple{T,Int}[]
  p=findfirst(!isempty,I.iterators)
  if p===nothing n=one(T) # the best we can do in this case
  else n=one(first(I.iterators[p]))
  end
  for it in I.iterators
    u=iterate(it)
    if u===nothing return u end
    n*=first(u)
    push!(state, (n,last(u)))
  end
  return n,state 
end

function Base.iterate(I::ProdIterator,state)
  for i in length(state):-1:1
    u=iterate(I.iterators[i],last(state[i]))
    if u===nothing continue end
    state[i]= i==1 ? u : (first(state[i-1])*first(u),last(u))
    for j in i+1:length(state)
      u=iterate(I.iterators[j])
      if u===nothing error() end
      state[j]=(first(state[j-1])*first(u),last(u))
    end
    return first(state[end]),state
  end
end

#------------------------- iteration for PermGroups -----------------------
function Base.length(G::PermGroup;type=Int)
  get!(G,:length)do
    prod(type.(length.(transversals(G))))
  end
end

Base.eltype(::Type{<:PermGroup{T}}) where T=Perm{T}

# iterating I directly is 25% faster unfortunately
function Base.iterate(G::PermGroup{T}) where T
  I=ProdIterator(reverse(values.(transversals(G))))
  u=iterate(I)
  if u===nothing return u end
  p,st=u
  p,(I,st)
end

function Base.iterate(G::PermGroup,(I,state))
  u=iterate(I,state)
  if u===nothing return u end
  p,st=u
  p,(I,st)
end

# elements is twice faster than collect(G), should not be
function Groups.elements(G::PermGroup)
  t=reverse(sort.(collect.(values.(transversals(G)))))
  if isempty(t) return [one(G)] end
  res=t[1]
  for i in 2:length(t)
    res=vcat(map(x->res.*x,t[i])...)
  end
  res
end

"""
`on_classes(G, aut)`

`aut`  is an automorphism of  the group `G` (for  a permutation group, this
could  be  given  as  a  permutation  normalizing  `G`).  The result is the
permutation of `1:nconjugacy_classes(G)` induced ny `aut`.

```julia-repl
julia> WF=rootdatum("3D4")
³D₄

julia> on_classes(Group(WF),WF.phi)
Perm{Int64}: (2,8,7)(5,13,12)
```
"""
on_classes(W, aut)=Perm(map(c->position_class(W,c^aut),classreps(W)))

#------------------------- cosets for PermGroups -----------------------

# computes "canonical" element of W.phi
function reduced(W::PermGroup,phi)
  for i in eachindex(base(W))
    t=transversals(W)[i]
    (kw,e)=minimum((k^phi,e) for (k,e) in t)
    phi=e*phi
  end
  phi
end

function Groups.Coset(W::PermGroup,phi::Perm=one(W))
  Groups.Cosetof(reduced(W,phi),W,Dict{Symbol,Any}())
end

#-------------------------- now a concrete type-------------------------
@GapObj struct PG{T}<:PermGroup{T}
  gens::Vector{Perm{T}}
  one::Perm{T}
end

function Groups.Group(a::AbstractVector{Perm{T}},one=one(Perm{T})) where T
  PG(filter(!isone,a),one,Dict{Symbol,Any}())
end

" The symmetric group of degree n "
function symmetric_group(n::Int)
  G=Group([Perm(i,i+1) for i in 1:n-1])
  G.name=string("symmetric_group(",n,")")
  G
end

#---------------- application to matrices ------------------------------
"""
`onmats(m::AbstractMatrix,g::Perm)` simultaneous action of `g` on 
 the columns and rows of `m`.
"""
onmats(m,g)=m^(g,g)

function invblocks(m,extra=nothing)
  if isnothing(extra) extra=zeros(Int,size(m,1)) end
  blk1=[collect(axes(m,1))]
  while true
    blk=blk1
    blk1=vcat(map(I->collectby(map(i->
                 (tally(m[i,I]),tally(m[I,i]),m[i,i],extra[i]),I),I), blk)...)
    if blk==blk1 return blk end
  end
end

# fast method for centralizer(g,M;action=onmats) additionally centralizing extra
function stab_onmats(g::PermGroup,M::AbstractMatrix,extra=nothing)
# if length(g)>1
#   print("g=",
#        length.(filter(x->length(x)>1,orbits(g,1:maximum(degree.(gens(g)))))),
#         " M:",size(M))
# end
  blocks=invblocks(M,extra)
  for r in blocks g=stabilizer(g,r) end
  if isempty(gens(g)) return g end
  n=vec(CartesianIndices(M))
  e=Group(map(y->Perm(n,map(x->CartesianIndex(x.I.^y),n)),gens(g)))
  for s in collectby(vec(M),1:length(M)) e=stabilizer(e,s) end
  if isempty(gens(e)) return e end
  Group(map(p->Perm(map(i->n[i^p][1], axes(M,1))), gens(e)))
end

"""
`stab_onmats([G,]M;extra=nothing)`

If  `onmats(m,p)=^(M,p;dims=(1,2))`, and  the argument  `G` is given (which
should   be  a  `PermGroup`)   this  is  just   a  fast  implementation  of
`centralizer(G,M;action=onmats)`.  If  `G`  is  omitted  it  is taken to be
`symmetric_group(size(M,1))`.  The  program  uses sophisticated algorithms,
and can handle matrices up to 80×80. If a list `extra` is given the result
centralizes also `extra`.

```julia-repl
julia> uc=UnipotentCharacters(complex_reflection_group(34));

julia> stab_onmats(fourier(uc.families[20]))
Group([(7,38), (39,44)(40,43)(41,42)])
```
"""
function stab_onmats(M,extra=nothing)
  k=length(M)
  blocks=sort(invblocks(M,extra),by=x->-length(x))
  g=PermGroup()
  I=Int[]
  for r in blocks
    if length(r)>7 InfoChevie("Large Block:$r\n")  end
    if length(r)>1
      gr=stab_onmats(symmetric_group(length(r)), M[r,r])
      g=Group(vcat(gens(g),gens(gr).^mappingPerm(eachindex(r),r)))
    end
    append!(I,r)
    p=mappingPerm(I,eachindex(I))
    g=stab_onmats(Group(gens(g).^p), M[I,I])
    g=Group(gens(g).^inv(p))
  end
  return g
end

"""
`Perm_onmats(M, N[, m ,n])` 

If  `onmats(M,p)=^(M,p;dims=(1,2))`, return `p`  such that `onmats(M,p)=N`;
so is just an efficient version of
`transporting_elt(symmetric_group(size(M,1)),M,N;action=onmats)`    If   in
addition the vectors `m` and `n` are given, `p` should satisfy `m^p=n`.

```julia-repl
julia> m=cartan(:D,12);

julia> n=^(m,Perm(1,5,2,8,12,4,7)*Perm(3,9,11,6);dims=(1,2));

julia> Perm_onmats(m,n)
(1,5,2,8,12,4,7)(3,9,11,6)
```
"""
function Perm_onmats(M,N,m=nothing,n=nothing)
  if isnothing(m) && M==N return Perm() end
  if size(M,1)!=size(M,2) || size(N,1)!=size(N,2)
    error("matrices are  not  square")
  end
  if size(M,1)!=size(N,1)
    @info "matrices do not have same dimensions"
    return nothing
  end
  sg=n->n==1 ? PermGroup() : Group(vcat(map(i->map(j->Perm(i,j),i+1:n),1:n-1)...))
  function ind(I,J)
    iM=map(i->[tally(M[i,I]), tally(M[I,i]), M[i,i]], I)
    iN=map(i->[tally(N[i,J]), tally(N[J,i]), N[i,i]], J)
    if !isnothing(m)
      iM=map(push!,iM,m[I])
      iN=map(push!,iN,n[J])
    end
    if tally(iM)!=tally(iN) return false end
    iM=collectby(iM,J)
    iN=collectby(iN,I)
    p=map(function(I,J)
      if length(I)>7
        InfoChevie("large block size $(length(I))\n")
        if length(iM)==1
          p=transporting_elt(sg(length(I)),M[I,I],N[J,J];
                action=onmats,dist=(M,N)->sum(x->count(!iszero,x),M-N))
        elseif isnothing(m) p=Perm_onmats(M[I,I],N[J,J])
        else p=Perm_onmats(M[I,I],N[J,J],m[I],n[J])
        end
      else 
        p=transporting_elt(sg(length(I)),M[I,I],N[J,J],action=onmats)
      end
      if isnothing(p) return false end
      I^=p
      p=mappingPerm(eachindex(I), I)
      return [I,J,Group(gens(stab_onmats(M[I,I])).^p)] end, iM, iN)
    if false in p return false else return p end
  end
  l=ind(axes(M,1),axes(N,1))
  if l==false return nothing end
  I=Int[]
  J=Int[]
  g=PermGroup()
  for r in l
    append!(I, r[1])
    append!(J, r[2])
    s=length(r[1])
    g=Group(vcat(gens(g), gens(r[3])))
    p=mappingPerm(I, eachindex(I))
    h=Group(gens(g).^p)
    if M[I,I]!=N[J,J]
      InfoChevie("I==$(length(I)) stab==$(length(g)) ")
      e = transporting_elt(h, M[I,I], N[J,J], action=onmats)
      if isnothing(e) return nothing else I^=e end
    end
    h=centralizer(h, M[I,I], action=onmats)
    g=Group(gens(h).^inv(p))
  end
  return mappingPerm(I,J)
end

end
