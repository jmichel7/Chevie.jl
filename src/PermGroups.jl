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
  5.995 ms (391728 allocations: 13.89 MiB)

julia> @btime minimal_words(symmetric_group(8));
  10.477 ms (122062 allocations: 15.22 MiB)
  
julia> @btime length(elements(symmetric_group(8)))
  2.136 ms (98328 allocations: 5.94 MiB)
```
Compare to GAP3 Elements(SymmetricGroup(8)); takes 3.8 ms
"""
module PermGroups
using ..Perms
using ..Groups
using ..Util: gets, getp, joindigits, InfoChevie
import ..Gapjm: degree, elements
using ..Combinat: tally, collectby
export PermGroup, base, transversals, centralizers, symmetric_group, reduced,
  stab_onmats, Perm_onmats, Perm_rowcolmat
#-------------------- now permutation groups -------------------------
abstract type PermGroup{T}<:Group{Perm{T}} end

PermGroup()=Group(Perm{Int16}[])

function Base.show(io::IO,G::PermGroup)
  print(io,"Group([");join(io,gens(G),',');print(io,"])")
# print(io,"Group(",gens(G),")")
end

Base.one(G::PermGroup{T}) where T=one(Perm{T})

function degree(G::PermGroup)::Int
  gets(G,:degree)do 
    maximum(largest_moved_point.(gens(G)))
  end
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
# if iterators=[i1,...,in] iterate on all products i1[j1]*...*in[jn]
struct ProdIterator{T}
  iterators::Vector{T}
end

Base.length(I::ProdIterator)=prod(length.(I.iterators))
Base.eltype(::Type{ProdIterator{T}}) where T=eltype(T)

function Base.iterate(I::ProdIterator)
  T=eltype(eltype(I.iterators))
  state=Tuple{T,Int}[]
  p=findfirst(!isempty,I.iterators)
  if isnothing(p) n=one(T) # the best we can do in this case
  else n=one(first(I.iterators[p]))
  end
  for it in I.iterators
    u=iterate(it)
    if isnothing(u) return u end
    n*=first(u)
    push!(state, (n,last(u)))
  end
  return n,state 
end

function Base.iterate(I::ProdIterator,state)
  for i in length(state):-1:1
    u=iterate(I.iterators[i],last(state[i]))
    if isnothing(u) continue end
    state[i]= i==1 ? u : (first(state[i-1])*first(u),last(u))
    for j in i+1:length(state)
      u=iterate(I.iterators[j])
      state[j]=(first(state[j-1])*first(u),last(u))
    end
    return first(state[end]),state
  end
end

#------------------------- iteration for PermGroups -----------------------
function Base.length(G::PermGroup)
  gets(G,:length)do
    prod(length.(transversals(G)))
  end
end

Base.eltype(::Type{<:PermGroup{T}}) where T=Perm{T}

# iterating I directly is faster, but how to do that?
function Base.iterate(G::PermGroup)
  I=ProdIterator(reverse(values.(transversals(G))))
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
function elements(G::PermGroup)
  t=reverse(values.(transversals(G)))
  if isempty(t) return [one(G)] end
  res=collect(t[1])
  for i in 2:length(t)
    res=vcat(map(x->res.*x,t[i])...)
  end
  res
end

function reduced(W::PermGroup,phi)
  for i in eachindex(base(W))
    t=transversals(W)[i]
    (kw,e)=minimum((k^phi,e) for (k,e) in t)
    phi=e*phi
  end
  phi
end

# computes "canonical" element of W.phi
function Groups.Coset(W::PermGroup,phi::Perm)
  Groups.CosetofAny(reduced(W,phi),W,Dict{Symbol,Any}())
end

#-------------------------- now the concrete type-------------------------
struct PG{T}<:PermGroup{T}
  gens::Vector{Perm{T}}
  prop::Dict{Symbol,Any}
end

function Groups.Group(a::AbstractVector{Perm{T}}) where T
  PG(filter(!isone,a),Dict{Symbol,Any}())
end

" The symmetric group of degree n "
symmetric_group(n::Int)=Group([Perm(i,i+1) for i in 1:n-1])

#---------------- application to matrices ------------------------------
onmats(m,g)=^(m,g;dims=(1,2))

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

function stab_onmats(g::Group,M::AbstractMatrix,extra=nothing)
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
`stab_onmats([G,]M[,l])`

If  `onmats(m,p)=^(M,p;dims=(1,2))`, and  the argument  `G` is given (which
should   be  a  `PermGroup`)   this  is  just   a  fast  implementation  of
`centralizer(G,M;action=onmats)`.  If  `G`  is  omitted  it  is taken to be
`symmetric_group(size(M,1))`.  The  program  uses sophisticated algorithms,
and can handle matrices up to 80×80.

```julia-repl
julia> uc=UnipotentCharacters(ComplexReflectionGroup(34));

julia> stab_onmats(fourier(uc.families[20]))
Group([(7,38),(39,44)(40,43)(41,42)])
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

If  `onmats(M,p)=^(M,p;dims=(1,2))`, return `p`  such that `onmats(M,p)=N`.
If  in  addition  the  vectors  `m`  and  `n` are given, `p` should satisfy
`m^p=n`.

Efficient version of 
`transporting_elt(symmetric_group(size(M,1)),M,N;action=onmats)`

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
  ind=function(I,J)local p
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
      local g
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
      g=stab_onmats(M[I,I])
      p=mappingPerm(eachindex(I), I)
      return [I,J,Group(gens(g).^p)] end, iM, iN)
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

"""
`Perm_rowcolmat(m1,m2)`
  whether matrix `m1` is conjugate to matrix `m2` by row/col permutations

  `m1`  and `m2` should be rectangular matrices of the same dimensions. The
  function   returns   a   pair   of   permutations   `[p1,p2]`  such  that
  `^(m1^p[1],p[2];dims=2)==m2`   if  such   permutations  exist,  `nothing`
  otherwise.
"""
function Perm_rowcolmat(m1, m2)
  if size(m1)!=size(m2) error("not same dimensions") end
  if isempty(m1) return [Perm(), Perm()] end
  dist(m,n)=count(i->m[i]!=n[i],eachindex(m))
  dist(m,n,dim,l)=dim==1 ? dist(m[l,:],n[l,:]) : dist(m[:,l],n[:,l])
  mm=[m1,m2]
  InfoChevie("# ", dist(m1, m2), "")
  rcperm=[Perm(), Perm()],[Perm(), Perm()]
  crg=Vector{Int}[],Vector{Int}[]
  crg1=[axes(m1,1)],[axes(m1,2)]
  while true
    crg=crg1
    crg1=Vector{Int}[],Vector{Int}[]
    for dim in 1:2
      for g in crg[dim]
        invars=map(1:2) do i
          invar=map(j->map(k->tally(dim==1 ? mm[i][j,k] : mm[i][k,j]),
                           crg[3-dim]), g)
          p=mappingPerm(vcat(collectby(invar,g)...), g)
          rcperm[dim][i]*=p
          mm[i]=^(mm[i],p,dims=dim)
          sort!(invar)
        end
        if invars[1]!=invars[2] return nothing end
        append!(crg1[dim], collectby(invars[1],g))
      end
    end
    InfoChevie("==>",dist(mm[1],mm[2]))
    if crg==crg1 break end
  end
  function best(l,dim)
    if length(l)==1 return false end
    d=dist(mm[1], mm[2], dim, l)
    for e in elements(Group(map(i->Perm(l[i],l[i+1]),1:length(l)-1)))
      m=dist(^(mm[1], e;dims=dim), mm[2], dim, l)
      if m<d
        InfoChevie("\n",("rows","cols")[dim],joindigits(l),":$d->",m)
        rcperm[dim][1]*=e
        mm[1]=^(mm[1],e;dims=dim)
        return true
      end
    end
    return false
  end
  while true
    s=false
    for dim in 1:2 for g in crg[dim] s=s || best(g,dim) end end
    if !s break end
  end
  InfoChevie("\n")
  if !iszero(dist(mm...)) error("RepresentativeRowColOperation failed") end
  [rcperm[1][1]/rcperm[1][2],rcperm[2][1]/rcperm[2][2]]
end

end
