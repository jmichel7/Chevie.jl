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
left  descent set  can be  easily determined  (see e.g. 'coxsym' below), so
this suggests how to deal with Coxeter groups.

The type `CoxeterGroup` is an abstact type; an actual struct which implements
it must define a function

`isleftdescent(W,w,i)` which tells whether the
      `i`-th element of `S` is in the left descending set of `w`.

the other functions needed in an instance of a Coxeter group are
- `coxgens(W)` which returns the set `S` (the list of *Coxeter generators*)
-  `nref(W)` which  returns the  number of  reflections of  `W`, if  `W` is
   finite or `nothing` if `W` is infinite

It  should be  noted that  a Coxeter group can be
*any* kind of group implementing the above functions.

A  common occurrence in code for Coxeter groups is a loop like:

`findfirst(eachindex(coxgens(W)),x->isleftdescent(W,w,x))`

if you provide a function `firstleftdescent(W,w)` it will be called instead
of the above loop.

Because  of the  easy solution  of the  word problem  in Coxeter  groups, a
convenient  way  to  represent  their  elements  is as words in the Coxeter
generators.  They are represented as lists of labels for the generators. By
default  these labels are  given as the  index of a  generator in `S`, so a
Coxeter  word is just  a list of  integers in `1:length(S)`. For reflection
subgroups, the labels are indices of the reflections in the parent group.

The functions 'word' and 'eltword' will do the conversion between
Coxeter words and elements of the group.

# Examples
```julia-repl
julia> W=coxsym(4)
coxsym(4)

julia> p=eltword(W,[1,3,2,1,3])
{UInt8}(1,4)

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
{UInt8}(1,4)(2,3)
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
in  CHEVIE. The only Coxeter group  constructor implemented here is coxsym.
The file Weyl.jl defines WeylGroup.

The dictionary from CHEVIE is as follows:
```
     CoxeterElements(W[,l])                → elements(W[,l])
     CoxeterLength(W,w)                    → length(W,w)
     CoxeterWord(W,w)                      → word(W,w)
     LongestCoxeterElement(W)              → longest(W)
     FirstLeftDescending(W,w)              → firstleftdescent(W,w)
     ReducedInRightCoset(W,w)              → reduced(W,w)
     ReducedRightCosetRepresentatives(W,H) → reduced(H,W)
     SemiSimpleRank(W)                     → coxrank(W)
     CoxeterGroupSymmetricGroup(n)         → coxsym(n)
     ReflectionSubgroup                    only standard parabolics now
     IsLeftDescending(W,w,i)               → isleftdescent(W,w,i)
     ReflectionDegrees(W)                  → degrees(W)
     ReflectionLength(W,w)                 → reflength(W,w)
     W.N                                   → nref(W)
```
"""
module CoxGroups

export CharTable, coxrank, CoxeterGroup, coxgens, coxsym, firstleftdescent, 
  longest, bruhatless, reduced, name

export isleftdescent, nref, ReflectionSubgroup, reflection,
  simple_representative # 'virtual' methods

using Gapjm
#-------------------------- Coxeter groups
abstract type CoxeterGroup{T}<:Group{T} end 

function firstleftdescent(W::CoxeterGroup,w)
  findfirst(i->isleftdescent(W,w,i),eachindex(coxgens(W)))
end

"""
  The Coxeter word for element w of W
"""
function Gapjm.word(W::CoxeterGroup,w)
  ww=Int[]
  while w!=one(W)
    i=firstleftdescent(W,w)
    push!(ww,i)
    w=coxgens(W)[i]*w
  end
  ww
end

Base.one(W::CoxeterGroup)=one(W.G)
Base.eltype(W::CoxeterGroup)=eltype(W.G)
coxgens(W::CoxeterGroup)=gens(W.G)
Gapjm.gens(W::CoxeterGroup)=coxgens(W)
coxrank(W::CoxeterGroup)=length(coxgens(W))
Base.show(io::IO, W::CoxeterGroup)=print(io,name(W))
function nref end

"""
The longest element of W --- never ends if W is infinite
"""
function longest(W::CoxeterGroup,I::AbstractVector{<:Integer}=
                  eachindex(coxgens(W)))
  w=one(W)
  i=1
  while i<=length(I)
    if isleftdescent(W,w,I[i]) i+=1
    else w=coxgens(W)[I[i]]*w
      i=1
    end
  end
  w
end

"""
reduced(W,w)
  The unique element in the coset W.w which stabilises the positive roots of W
"""
function reduced(W::CoxeterGroup,w)
  while true
    i=firstleftdescent(W, w)
    if isnothing(i) return w end
    w = coxgens(W)[i] * w
  end
end

"""
reduced(H,W)
  The elements in W which are H-reduced
"""
function reduced(H::CoxeterGroup,W::CoxeterGroup)
  res=[Set([one(W)])]
  while true
    new=vcat(map(collect(res[end])) do w
      l=filter(i->!isleftdescent(W,inv(w),i),eachindex(coxgens(W)))
      l=filter(x->x==reduced(H,x),Ref(w).*coxgens(W)[l])
    end...)
    new=Set(new)
    if isempty(new) break
    else push!(res,new)
    end
  end
  vcat(res)
end

function Gapjm.elements(W::CoxeterGroup, l::Int)
  elts=gets(W,:elements)do W Dict(0=>[one(W)]) end
  if haskey(elts,l) return Vector{typeof(one(W))}(elts[l]) end
  elts[l]=typeof(one(W))[]
  H=gets(W,:maxpara)do W 
    ReflectionSubgroup(W,collect(1:length(coxgens(W))-1))
  end
  rc=gets(W,:rc)do W reduced(H,W) end
  N=nref(W)
  for i in max(0,l+1-length(rc)):min(l,N+1-length(rc))
    e=elements(H,i)
    for x in rc[1+l-i], w in e push!(elts[l],w*x) end
#   2 times slower variant
#   for x in rc[1+l-i] append!(elts[l],e.*x) end
  end
  if N-l>l elts[N-l]=elts[l].*Ref(longest(W)) end
  Vector{typeof(one(W))}(elts[l])
end

function Gapjm.elements(W::CoxeterGroup)
  vcat(map(i->elements(W,i),0:nref(W))...)
end

function Gapjm.words(W::CoxeterGroup, l::Int)
  ww=gets(W,:words)do W Dict(0=>[Int[]]) end
  if haskey(ww,l) return Vector{Vector{Int}}(ww[l]) end
  ww[l]=Vector{Int}[]
  H=gets(W,:maxpara)do W ReflectionSubgroup(W,collect(1:length(coxgens(W))-1))
  end
  rc=gets(W->[[word(W,x) for x in y] for y in reduced(H,W)],W,:rcwords)
  N=nref(W)
  for i in max(0,l+1-length(rc)):min(l,N+1-length(rc))
    e=words(H,i)
    for x in rc[1+l-i], w in e push!(ww[l],vcat(w,x)) end
#   somewhat slower variant
#   for x in rc[1+l-i] append!(ww[l],vcat.(e,(x,))) end
  end
# if N-l>l ww[N-l]=ww[l].*longest(W) end
  Vector{Vector{Int}}(ww[l])
end

function Gapjm.words(W::CoxeterGroup)
  vcat(map(i->words(W,i),0:nref(W))...)
end

"""
   `bruhatless(W, x, y)`  whether x≤y in the Bruhat order, for x, y ∈ W.
"""
function bruhatless(W::CoxeterGroup,x,y)
  if x==one(W) return true end
  d=length(W,y)-length(W,x)
  while d>0
    i=firstleftdescent(W,y)
    s=coxgens(W)[i]
    if isleftdescent(W,x,i)
      if x==s return true end
      x=s*x 
    else d-=1
    end
    y=s*y
  end
  return x==y
end

#--------------------- CoxSymmetricGroup ---------------------------------
struct CoxSymmetricGroup{T} <: CoxeterGroup{Perm{T}}
  G::PermGroup{T}
  n::Int
  prop::Dict{Symbol,Any}
end

"The symmetric group on n letters as a Coxeter group"
function coxsym(n::Int)
  gens=map(i->Perm{UInt8}(i,i+1),1:n-1)
  CoxSymmetricGroup{UInt8}(PermGroup(gens),n,Dict{Symbol,Any}(:name=>"coxsym($n)"))
end

name(W::CoxSymmetricGroup)::String=W.prop[:name]
coxgens(W::CoxSymmetricGroup)=W.G.gens
nref(W::CoxSymmetricGroup)=div(W.n*(W.n-1),2)

function isleftdescent(W::CoxSymmetricGroup,w,i::Int)
  i^w>(i+1)^w
end

Gapjm.degrees(W::CoxSymmetricGroup)=2:length(coxgens(W))+1

Base.length(W::CoxSymmetricGroup)=prod(degrees(W))

function Base.length(W::CoxSymmetricGroup,w)
  l=0
  for j in 1 .+eachindex(coxgens(W)) for i in 1:j-1 
    l+=i^w>j^w
  end end
  l
end
# 2.5 times longer on 1.0.2
function length2(W::CoxSymmetricGroup,w)
  count(i^w>j^w for j in 1 .+eachindex(coxgens(W)) for i in 1:j-1) 
end

function reflength(W::CoxSymmetricGroup,w::Perm{T})where T
 sum(x->x-1,length.(cycles(w)))
end

" Only parabolics defined are I=1:m for m≤n"
function ReflectionSubgroup(W::CoxSymmetricGroup,I::Vector{Int})
  if length(I)>0 n=maximum(I) 
    if I!=1:n error(I," should be 1:n for some n") end
  else n=0 end
  CoxSymmetricGroup(PermGroup(coxgens(W)[I]),n+1,Dict{Symbol,Any}(:name=>name(W)*"_$I",))
end
  
simple_representative(W::CoxSymmetricGroup,i)=1

function reflection(W::CoxSymmetricGroup,i::Int)
  (gets(W,:reflections)do W
    [Perm(i,i+k) for k in 1:W.n-1 for i in 1:W.n-k]
  end)[i]
end

end
