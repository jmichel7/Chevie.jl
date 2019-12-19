"""
A  suitable  reference  for  the  general  theory of Coxeter groups is, for
example, Bourbaki "Lie Groups and Lie Algebras" chapter 4.

A *Coxeter group* is a group which has the presentation
`W=⟨S|(st)^m(s,t)=1`  for  `s,t∈  S⟩`  for  some  symmetric  integer matrix
`m(s,t)`  called  the  *Coxeter  matrix*,  where  `m(s,t)>1`  for `s≠t` and
`m(s,s)=1`.  It is true (but a non-trivial theorem) that in a Coxeter group
the  order of `st` is exactly `m(s,t)`, thus a Coxeter group is the same as
a  *Coxeter system*, that is a pair `(W,S)` of a group `W` and a set `S` of
involutions,  such that the group is  presented by relations describing the
order  of the product of two elements of `S`. A Coxeter group has a natural
representation, its *reflection representation*, on a real vector space `V`
of  dimension `length(S)` (the *Coxeter rank*  of W), where each element of
`S`  acts as a  reflection; the faithfulness  of this representation in the
main  argument to prove  that the order  of `st` is  exactly `m(s,t)`. Thus
Coxeter groups are real reflection groups. The converse need not be true if
the  set of reflecting  hyperplanes has bad  topological properties, but it
turns out that finite Coxeter groups are the same as finite real reflection
groups.  The possible Coxeter matrices for  finite Coxeter groups have been
completely  classified; the corresponding finite groups play a deep role in
several areas of mathematics.

Coxeter  groups  have  a  nice  solution  to the word problem. The *length*
`l(w)`  of an element  `w∈ W` is  the minimum number  of elements of `S` of
which it is a product (since the elements of `S` are involutions, we do not
need inverses). An expression of `w` of minimal length is called a *reduced
word*  for `w`. The main property of  reduced words is the *exchange lemma*
which  states that if `s₁…sₖ` is a  reduced word for `w` (thus`k=l(w)`) and
`s∈  S` is such that `l(sw)≤l(w)` then one  of the `sᵢ` in the word for `w`
can be deleted to obtain a reduced word for `sw`. Thus given `s∈ S` and `w∈
W`,  either `l(sw)=l(w)+1` or  `l(sw)=l(w)-1` and we  say in this last case
that  `s` belongs to  the *left descent  set* of `w`.  The computation of a
reduced word for an element, and other word problems, are easily done if we
know  the left descent sets. For the Coxeter groups that we implement, this
left  descent set  can be  easily determined  (see e.g. 'CoxSym' below), so
this suggests how to deal with Coxeter groups.

The type `CoxeterGroup` is an abstact type; an actual struct which implements
it must define a function

`isleftdescent(W,w,i)` which tells whether the
      `i`-th element of `S` is in the left descending set of `w`.

the other functions needed in an instance of a Coxeter group are
- `gens(W)` which returns the set `S` (the list of *Coxeter generators*)
- `nref(W)` which  returns the  number of  reflections of  `W`, if  `W` is
   finite or `nothing` if `W` is infinite

It  should be  noted that  a Coxeter group can be
*any* kind of group implementing the above functions.

A  common occurrence in code for Coxeter groups is a loop like:

`findfirst(eachindex(gens(W)),x->isleftdescent(W,w,x))`

if you provide a function `firstleftdescent(W,w)` it will be called instead
of the above loop.

Because  of the  easy solution  of the  word problem  in Coxeter  groups, a
convenient  way  to  represent  their  elements  is as words in the Coxeter
generators.  They are represented as lists of labels for the generators. By
default  these labels are  given as the  index of a  generator in `S`, so a
Coxeter  word is just  a list of  integers in `1:length(S)`. For reflection
subgroups, the labels are indices of the reflections in the parent group.

The functions 'word' and 'W(...)' will do the conversion between
Coxeter words and elements of the group.

# Examples
```julia-repl
julia> W=CoxSym(4)
𝔖 ₄

julia> p=W(1,3,2,1,3)
Perm{UInt8}: (1,4)

julia> word(W,p)
5-element Array{Int64,1}:
 1
 2
 3
 2
 1

```
We  notice that the word we started with and the one that we ended up with,
are not the same, though they represent the same element of `W`. The reason
is  that the function 'word' computes a lexicographically smallest word for
`w`.  Below  are  some  other  possible  computations with the same Coxeter
group:

```julia-repl
julia> word(W,longest(W))  # the (unique) longest element in W
6-element Array{Int64,1}:
 1
 2
 1
 3
 2
 1

julia> w0=longest(W)
Perm{UInt8}: (1,4)(2,3)
julia> length(W,w0)
6
julia> map(i->word(W,reflection(W,i)),1:nref(W))
6-element Array{Array{Int64,1},1}:
 [1]            
 [2]            
 [3]            
 [1, 2, 1]      
 [2, 3, 2]      
 [1, 2, 3, 2, 1]
julia> [length(elements(W,i)) for i in 0:nref(W)]
7-element Array{Int64,1}:
 1
 3
 5
 6
 5
 3
 1

```

The above line tells us that there is 1 element of length 0, there are 6 of
length 3, …

For  most basic functions the convention is that the input is an element of
the  group, rather than  a Coxeter word.  The reason is  that for a Coxeter
group  which  is  a  permutation  group,  using the low level functions for
permutations  is usually  much faster  than manipulating lists representing
reduced expressions.

This  file contains mostly a port of  the basic functions on Coxeter groups
in  Chevie. The only Coxeter group  constructor implemented here is CoxSym.
The file Weyl.jl defines coxgroup.

The dictionary from Chevie is as follows:
```
     CoxeterElements(W[,l])                → elements(W[,l])
     CoxeterLength(W,w)                    → length(W,w)
     CoxeterWord(W,w)                      → word(W,w)
     LongestCoxeterElement(W)              → longest(W)
     ReducedExpressions(W,w)               → reduced_words(W,w)
     FirstLeftDescending(W,w)              → firstleftdescent(W,w)
     LeftDescenTSet(W,w)                   → leftdescents(W,w)
     ReducedInRightCoset(W,w)              → reduced(W,w)
     ReducedRightCosetRepresentatives(W,H) → reduced(H,W)
     SemiSimpleRank(W)                     → coxrank(W)
     CoxeterGroupSymmetricGroup(n)         → CoxSym(n)
     CoxeterGroupHyperoctaedralGroup(n)    → CoxHyperoctaedral(n)
     ReflectionSubGroup                    only standard parabolics now
     IsLeftDescending(W,w,i)               → isleftdescent(W,w,i)
     ReflectionDegrees(W)                  → degrees(W)
     ReflectionLength(W,w)                 → reflength(W,w)
     W.N                                   → nref(W)
```
"""
module CoxGroups

export bruhatless, CoxeterGroup, coxrank, firstleftdescent, leftdescents, 
  longest, reduced, braid_relations, coxetermat, parabolic_representatives,
  CoxSym, reduced_words

export isleftdescent, nref # 'virtual' methods (exist only for concrete types)

using Gapjm
#-------------------------- Coxeter groups
abstract type CoxeterGroup{T}<:Group{T} end 

function firstleftdescent(W::CoxeterGroup,w)
  findfirst(i->isleftdescent(W,w,i),eachindex(gens(W)))
end

function leftdescents(W::CoxeterGroup,w)
  # 3 times faster than filter
  [i for i in eachindex(gens(W)) if isleftdescent(W,w,i)]
end

isrightdescent(W::CoxeterGroup,w,i)=isleftdescent(W,inv(w),i)

"""
  word(W::CoxeterGroup,w)

returns  a reduced word in the standard generators of the Coxeter group `W`
for  the  element  `w`  (represented  as  the  vector  of the corresponding
generator indices).

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> w=perm"(1,11)(3,10)(4,9)(5,7)(6,12)"
(1,11)(3,10)(4,9)(5,7)(6,12)

julia> w in W
true

julia> word(W,w)
5-element Array{Int64,1}:
 1
 2
 3
 2
 1
```
The  result  of   `word`  is  the  lexicographically  smallest reduced word
for~`w` (for the ordering of the Coxeter generators given by `gens(W)`).
"""
function Gapjm.word(W::CoxeterGroup,w)
  ww=Int[]
  while w!=one(W)
    i=firstleftdescent(W,w)
    push!(ww,i)
    w=W(i)*w
  end
  ww
end

Base.length(W::CoxeterGroup,w)=length(word(W,w))
Base.one(W::CoxeterGroup)=one(W.G)
Base.eltype(W::CoxeterGroup)=eltype(W.G)
Gapjm.gens(W::CoxeterGroup)=gens(W.G)
coxrank(W::CoxeterGroup)=length(gens(W))
function nref end

"""
The longest element of reflection_subgroup(W,I) --- never ends if infinite
"""
function longest(W::CoxeterGroup,I::AbstractVector{<:Integer}=eachindex(gens(W)))
  w=one(W)
  i=1
  while i<=length(I)
    if isleftdescent(W,w,I[i]) i+=1
    else w=W(I[i])*w
      i=1
    end
  end
  w
end

"""
reduced(W,w)
  The unique element in the coset W.w which stabilises the positive roots of W
```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> H=reflection_subgroup(W,[2,6])
G₂₍₂₆₎

julia> word.(Ref(W),Set(reduced.(Ref(H),elements(W))))
3-element Array{Array{Int64,1},1}:
 []    
 [1, 2]
 [1]   
```
"""
function reduced(W::CoxeterGroup,w)
  while true
    i=firstleftdescent(W, w)
    if isnothing(i) return w end
    w = W(i) * w
  end
end

"""
reduced(H,W)
  The elements in W which are H-reduced
```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> H=reflection_subgroup(W,[2,6])
G₂₍₂₆₎

julia> [word(W,w) for S in reduced(H,W) for w in S]
3-element Array{Array{Int64,1},1}:
 []    
 [1]   
 [1, 2]
```
"""
function reduced(H::CoxeterGroup,W::CoxeterGroup)
  res=[Set([one(W)])]
  while true
    new=reduced(H,W,res[end])
    if isempty(new) break
    else push!(res,new)
    end
  end
  vcat(res)
end

"""
reduced(H,W,S)
  The elements in W which are H-reduced of length i from the set S of length i-1
"""
function reduced(H::CoxeterGroup,W::CoxeterGroup,S)
  res=empty(S)
  for w in S
    for i in eachindex(gens(W))
      if !isrightdescent(W,w,i)
        w1=w*W(i)
        if w1==reduced(H,w1) push!(res,w1) end
      end
    end
  end
  res
end

function Gapjm.elements(W::CoxeterGroup{T}, l::Int)::Vector{T} where T
  elts=gets(W->Dict(0=>[one(W)]),W,:elements)::Dict{Int,Vector{T}}
  if haskey(elts,l) return elts[l] end
  if coxrank(W)==1 return l>1 ? T[] : gens(W) end
  H=gets(W->reflection_subgroup(W,1:coxrank(W)-1),W,:maxpara)::CoxeterGroup{T}
  rc=gets(W->[Set([one(W)])],W,:rc)::Vector{Set{T}}
  while length(rc)<=l
    new=reduced(H,W,rc[end])
    if isempty(new) break
    else push!(rc,new)
    end
  end
# println("l=$l W=$W H=$H rc=$rc")
  elts[l]=T[]
  for i in max(0,l+1-length(rc)):l
    for x in rc[1+l-i] append!(elts[l],elements(H,i).*Ref(x)) end
  end
# N=nref(W)
# if !isnothing(N) && N-l>l elts[N-l]=elts[l].*longest(W) end
  elts[l]
end

function Gapjm.elements(W::CoxeterGroup)
  reduce(vcat,map(i->elements(W,i),0:nref(W)))
end

const Wtype=Vector{Int8}
function Gapjm.words(W::CoxeterGroup{T}, l::Int)where T
 ww=gets(W->Dict(0=>[Wtype([])]),W,:words)::Dict{Int,Vector{Wtype}}
  if haskey(ww,l) return ww[l] end
  if coxrank(W)==1
    if l!=1 return Vector{Int}[]
    else return [[1]]
    end
  end
  H=gets(W->reflection_subgroup(W,1:coxrank(W)-1),W,:maxpara)::CoxeterGroup{T}
  rc=gets(W->[[Wtype([])]],W,:rcwords)::Vector{Vector{Wtype}}
  while length(rc)<=l
    new=reduced(H,W,Set((x->W(x...)).(rc[end])))
    if isempty(new) break
    else push!(rc,Wtype.(word.(Ref(W),new)))
    end
  end
  ww[l]=Wtype([])
  for i in max(0,l+1-length(rc)):l
    e=words(H,i)
    for x in rc[1+l-i], w in e push!(ww[l],vcat(w,x)) end
#   somewhat slower variant
#   for x in rc[1+l-i] append!(ww[l],vcat.(e,Ref(x))) end
  end
  ww[l]
end

function Gapjm.words(W::CoxeterGroup)
  reduce(vcat,map(i->words(W,i),0:nref(W)))
end

"""
   `bruhatless(W, x, y)`  whether x≤y in the Bruhat order, for x, y ∈ W.
"""
function bruhatless(W::CoxeterGroup,x,y)
  if x==one(W) return true end
  d=length(W,y)-length(W,x)
  while d>0
    i=firstleftdescent(W,y)
    s=W(i)
    if isleftdescent(W,x,i)
      if x==s return true end
      x=s*x 
    else d-=1
    end
    y=s*y
  end
  return x==y
end

"""
   `bruhatless(W, y)`  all x≤y in the Bruhat order, for y ∈ W.
"""
function bruhatless(W::CoxeterGroup,w)
  if w==one(W) return [[w]] end
  i=firstleftdescent(W,w)
  s=W(i)
  res=bruhatless(W,s*w)
  for j in 1:length(res)-1
    res[j+1]=union(res[j+1],s.*filter(x->!isleftdescent(W,x,i),res[j]))
  end
  push!(res,s.*filter(x->!isleftdescent(W,x,i),res[end]))
end

function reduced_words(W::CoxeterGroup,w)
  l=leftdescents(W,w)
  if isempty(l) return [Int[]] end
  reduce(vcat,map(x->vcat.(Ref([x]),reduced_words(W,W(x)*w)),l))
end

"diagram of finite Coxeter group"
PermRoot.Diagram(W::CoxeterGroup)=Diagram(refltype(W))

# return all subsets of S which are W-conjugate to I
function standard_parabolic_class(W,I::Vector{Int})
   new=[I]
   res=empty(new)
   while !isempty(new)
     n=vcat(map(new) do K
       J=map(setdiff(1:coxrank(W),K)) do i
	 w0=longest(W,push!(copy(K),i)::Vector{Int})
         if isnothing(w0) one(w) else w0 end 
         # isnothing==parabolic W_I1 infinite
       end
       map(J) do w
         H=gens(W)[K].^w
         findall(x->x in H,gens(W))
       end
       end...)
       res=union(res,new)
       new=setdiff(n,res)
   end
   res
end

function parabolic_category(W,I::AbstractVector{<:Integer})
  Category(collect(sort(I));action=(J,e)->sort(J.^e))do J
    map(setdiff(1:coxrank(W),J)) do i
      longest(W,J)*longest(W,push!(copy(J),i))
    end
  end
end
    
# representatives of parabolic classes
function parabolic_representatives(W,s)
  l=collect(combinations(1:coxrank(W),s))
  orbits=[]
  while !isempty(l) 
    o=parabolic_category(W,l[1]).obj
    push!(orbits,o)
    l=setdiff(l,o)
  end
  first.(orbits)
end

parabolic_representatives(W)=union(parabolic_representatives.(Ref(W),
          0:semisimplerank(W))...)

function coxetermat end

function braid_relations(W::CoxeterGroup)
  p(i,j,b)=map(k->iszero(k%2) ? j : i,1:b)
  m=coxetermat(W)
  reduce(vcat,map(i->map(j->[p(i,j,m[i,j]),p(j,i,m[i,j])],1:i-1),1:size(m,1)))
end
#--------------------- CoxSymmetricGroup ---------------------------------
struct CoxSym{T} <: CoxeterGroup{Perm{T}}
  G::PermGroup{T}
  n::Int
  prop::Dict{Symbol,Any}
end

Base.iterate(W::CoxSym,r...)=iterate(W.G,r...)

"""
  `Coxsym(n)` The symmetric group on `n` letters as a Coxeter group
"""
function CoxSym(n::Int)
  CoxSym{UInt8}(Group([Perm{UInt8}(i,i+1) for i in 1:n-1]),n,Dict{Symbol,Any}())
end

function Base.show(io::IO, W::CoxSym)
  name="\\frakS _{$(W.n)}"
  if get(io,:TeX,false) print(io,name)
  elseif get(io,:limit,false) print(io,TeXstrip(name))
  else print(io,"CoxSym($(W.n))")
  end
end

PermRoot.refltype(W::CoxSym)=[TypeIrred(Dict(:series=>:A,
                                        :indices=>collect(1:W.n-1)))]

Perms.reflength(W::CoxSym,a)=reflength(a)
  
nref(W::CoxSym)=div(W.n*(W.n-1),2)

isleftdescent(W::CoxSym,w,i::Int)=i^w>(i+1)^w

Gapjm.degrees(W::CoxSym)=2:length(gens(W))+1

Base.length(W::CoxSym)=prod(degrees(W))

function Base.length(W::CoxSym,w)
  l=0
  for k in 1:W.n-1 for i in 1:W.n-k
    l+=i^w>(i+k)^w
  end end
  l
end
# 3 times longer on 1.3
function length2(W::CoxSym,w)
  count(i^w>(i+k)^w for k in 1:W.n-1 for i in 1:W.n-k)
end

" Only parabolics defined are I=1:m for m≤n"
function PermRoot.reflection_subgroup(W::CoxSym,I::AbstractVector{Int})
  if length(I)>0 n=maximum(I) 
    if I!=1:n error(I," should be 1:n for some n") end
  else n=0 end
  CoxSym(Group(gens(W)[I]),n+1,Dict{Symbol,Any}())
end
  
PermRoot.simple_representatives(W::CoxSym)=fill(1,nref(W))

function PermRoot.reflection(W::CoxSym{T},i::Int)where T
  ref=gets(W,:reflections)do W
    [Perm{T}(i,i+k) for k in 1:W.n-1 for i in 1:W.n-k]
  end::Vector{Perm{T}}
  ref[i]
end

PermRoot.reflections(W::CoxSym)=reflection.(Ref(W),1:nref(W))

end
