"""
A  suitable  reference  for  the  general  theory of Coxeter groups is, for
example, Bourbaki "Lie Groups and Lie Algebras" chapter 4.

A *Coxeter group* is a group which has the presentation
``W=⟨S|(st)^{m(s,t)}=1\\text{  for  }s,t∈  S⟩``  for some symmetric integer
matrix `m(s,t)` called the *Coxeter matrix*, where `m(s,t)>1` for `s≠t` and
`m(s,s)=1`;  `m(s,t)=∞` is allowed meaning there is no relation between `s`
and `t`. It is true (but a non-trivial theorem) that in a Coxeter group the
order  of `st` is exactly  `m(s,t)`, thus a Coxeter  group is the same as a
*Coxeter  system*, that is a pair `(W,S)` of a group `W` and a set `S⊂W` of
involutions,  such  that  the  group  is  presented  by  generators `S` and
relations describing the order of the product of two elements of `S`.

A   Coxeter   group   has   a   natural   representation,  its  *reflection
representation*, on a real vector space `V` of dimension `length(S)` (which
is  the  *Coxeter  rank*  of  W),  where  each  element  of  `S`  acts as a
reflection; the faithfulness of this representation in the main argument to
prove  that the order  of `st` is  exactly `m(s,t)`. This representation is
defined as follows on a space `V` with basis `{eₛ}` for `s∈ S`. The *cartan
matrix*  associated to the  Coxeter matrix `m(s,t)`  is the matrix `C` with
entries  `C(s,t)=-2cos(π/m(s,t))`; we  set `C(s,t)=-2`  if `m(s,t)=∞`. Then
the action of `s∈ S` on `V` is given by `s(eₜ)=eₜ-C(s,t)eₛ`.

Thus, Coxeter groups are  real reflection groups.  The converse need not be
true  if the set of reflecting  hyperplanes has bad topological properties,
but  it turns out  that finite Coxeter  groups are the  same as finite real
reflection  groups. The possible Coxeter matrices for finite Coxeter groups
have  been  completely  classified,  see  [`Weyl`](@ref); the corresponding
finite groups play a deep role in several areas of mathematics.

Coxeter  groups  have  a  nice  solution  to the word problem. The *length*
`l(w)`  of an element  `w∈ W` is  the minimum number  of elements of `S` of
which it is a product (since the elements of `S` are involutions, we do not
need inverses). An expression of `w` of minimum length is called a *reduced
word*  for `w`. The main property of  reduced words is the *exchange lemma*
which  states that if `s₁…sₖ` is a reduced word for `w` (thus `k=l(w)`) and
`s∈  S` is such that `l(sw)≤l(w)` then one  of the `sᵢ` in the word for `w`
can be deleted to obtain a reduced word for `sw`. Thus given `s∈ S` and `w∈
W`,  either `l(sw)=l(w)+1` or `l(sw)=l(w)-1` and  in the latter case we say
that `s` belongs to the *left descent set* of `w`. Computing a reduced word
for  an  element,  and  other  word  problems,  are  easy if we know how to
multiply elements and know left descent sets. In each of the Coxeter groups
that we implement, the left descent set is easy to compute (see for example
[`CoxSym`](@ref)  below), so this suggests how  to deal with Coxeter groups
generically:

The  type  `CoxeterGroup`  is  an  abstract  type;  an  actual struct which
implements it must define a function

`isleftdescent(W,w,i)` which tells whether the `i`-th element of `S` is in
   the left descent set of `w`.

the other functions needed in an instance of a Coxeter group are
- `gens(W)` which returns the set `S` (the list of *Coxeter generators*)
- `nref(W)` which  returns the  number of  reflections of  `W`, if  `W` is
   finite or `nothing` if `W` is infinite.

It  should  be  noted  that  a  Coxeter  group  can  be *any* kind of group
implementing the above functions.

Because  of the  easy solution  of the  word problem  in Coxeter  groups, a
convenient  way  to  represent  their  elements  is as words in the Coxeter
generators,  that  is  lists  of  integers  in `1:length(S)`. The functions
'word'  and 'W(...)' do the conversion between Coxeter words and elements
of the group.

# Examples
```julia-repl
julia> W=CoxSym(4)
𝔖 ₄

julia> p=W(1,3,2,1,3)
(1,4)

julia> word(W,p)
5-element Vector{Int64}:
 1
 2
 3
 2
 1
```
We  notice that the word we started with and the one that we ended up with,
are  not the same, even though they  represent the same element of `W`. The
reason  is that there are several reduced  words for an element of `W`. The
function 'word' calculates a lexicographically smallest word for `w`. Below
are some other possible computations using the same Coxeter group:

```julia-repl
julia> word(W,longest(W))  # the (unique) longest element in W
6-element Vector{Int64}:
 1
 2
 1
 3
 2
 1

julia> w0=longest(W)
(1,4)(2,3)

julia> length(W,w0)
6
julia> map(w->word(W,w),refls(W,1:nref(W)))
6-element Vector{Vector{Int64}}:
 [1]
 [2]
 [3]
 [1, 2, 1]
 [2, 3, 2]
 [1, 2, 3, 2, 1]
julia> [length(elements(W,i)) for i in 0:nref(W)]
7-element Vector{Int64}:
 1
 3
 5
 6
 5
 3
 1
```

The  last list tells us that there is 1 element of length 0, there are 6 of
length 3, …

For  most basic functions the convention is that the input is an element of
the  group, rather than a  Coxeter word. The reason  for this is that for a
Coxeter  group which is a permutation  group, using the low level functions
for   permutations  is   usually  much   faster  than   manipulating  lists
representing reduced expressions.

The only Coxeter group constructors implemented in this module are `CoxSym`
and  `coxgroup`; the last constructor takes  a Cartan matrix and builds the
corersponding  Coxeter group as  a matrix group.  The module [`Weyl`](@ref)
defines  other methods for `coxgroup` building  a finite Coxeter group as a
permutation group, given its type.
"""
module CoxGroups

export bruhatless, CoxeterGroup, firstleftdescent, leftdescents,
  longest, braid_relations, coxmat, coxeter_matrix,
  CoxSym, standard_parabolic_class,
  inversions, degrees, FiniteCoxeterGroup,
  coxgroup, coxeter_group

export isleftdescent # 'virtual' methods (exist only for concrete types)

using ..Gapjm
#-------------------------- Coxeter groups
abstract type CoxeterGroup{T}<:Group{T} end

"""
`firstleftdescent(W,w)`

returns the index in `gens(W)` of the first element of the left descent set
of  `w` --- that is, the first  `i` such that if `s=W(i)` then `l(sw)<l(w).
It returns `nothing` for `one(W)`.

```julia-repl
julia> W=CoxSym(3)
𝔖 ₃

julia> firstleftdescent(W,Perm(2,3))
2
```
"""
function firstleftdescent(W::CoxeterGroup{T},w::T)where T
  findfirst(i->isleftdescent(W,w,i),eachindex(gens(W)::Vector{T}))
end

"""
`leftdescents(W,w)`

The  left descents of the element `w` of the Coxeter group `W`, that is the
set of `i` such that `length(W,W(i)*w)<length(W,w)`.

```julia-repl
julia> W=CoxSym(3)
𝔖 ₃

julia> leftdescents(W,Perm(1,3))
2-element Vector{Int64}:
 1
 2
```
"""
function leftdescents(W::CoxeterGroup{T},w)where T
  filter(i->isleftdescent(W,w,i),eachindex(gens(W)::Vector{T}))
end

isrightdescent(W::CoxeterGroup,w,i)=isleftdescent(W,inv(w),i)

"""
`word(W::CoxeterGroup,w)`

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
5-element Vector{Int64}:
 1
 2
 3
 2
 1
```
The  result  of   `word`  is  the  lexicographically  smallest reduced word
for `w` (for the ordering of the Coxeter generators given by `gens(W)`).
"""
function Groups.word(W::CoxeterGroup{T},w)where T
  ww=Int[]
  while true
    i=firstleftdescent(W,w)
    if i===nothing return ww end
    push!(ww,i)
    w=(gens(W)::Vector{T})[i]*w
  end
end

"""
`length(W::CoxeterGroup ,w)`

returns the length of a reduced expression in the Coxeter generators of the
element `w` of `W`.

```julia-repl
julia> W=CoxSym(4)
𝔖 ₄

julia> p=W(1,2,3,1,2,3)
(1,3)(2,4)

julia> length(W,p)
4

julia> word(W,p)
4-element Vector{Int64}:
 2
 1
 3
 2
```
"""
Base.length(W::CoxeterGroup,w)=length(word(W,w))
PermRoot.semisimplerank(W::CoxeterGroup)=ngens(W)

"""
`reduced(W,w)`

The  unique element of  minimal length in  the coset W.w.  This makes sense
when  `isleftdescent(W,u)` makes sense for `u∈  W.w` which happens when `w`
is  a Coxeter automorphism of `W` or  when `w` lives in a Coxeter overgroup
of `W`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> H=reflection_subgroup(W,[2,6])
G₂₍₂₆₎=Ã₁×A₁

julia> word.(Ref(W),unique(reduced.(Ref(H),elements(W))))
3-element Vector{Vector{Int64}}:
 []
 [1]
 [1, 2]
```
"""
function PermGroups.reduced(W::CoxeterGroup,w)
  while true
    i=firstleftdescent(W,w)
    if i===nothing return w end
    w=W(i)*w
  end
end

"""
`reduced(H::CoxeterGroup,W::CoxeterGroup,i=nref(W))`

The  elements `w∈ W` which are `H`-reduced, and of length `≤i`
(by default all of them), grouped by length.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> H=reflection_subgroup(W,[2,6])
G₂₍₂₆₎=Ã₁×A₁

julia> [word(W,w) for S in reduced(H,W) for w in S]
3-element Vector{Vector{Int64}}:
 []
 [1]
 [1, 2]
```
"""
function PermGroups.reduced(H::CoxeterGroup,W::CoxeterGroup,i::Integer=nref(W))
  res=[[one(W)]]
  while length(res)<=i
    new=reducedfrom(H,W,res[end])
    if isempty(new) break
    else push!(res,new)
    end
  end
  vcat(res)
end

#reduced(H,W,S)
#  The elements in `W` which are `H`-reduced of length `i` given the set `S`
#  of those of length `i-1`
function reducedfrom(H::CoxeterGroup,W::CoxeterGroup,S)
  res=empty(S)
  for w in S
    for i in eachindex(gens(W))
      if !isrightdescent(W,w,i)
        w1=w*W(i)
        if w1==reduced(H,w1) push!(res,w1) end
      end
    end
  end
  unique(res)
end

maxpara(W::CoxeterGroup)=2:ngens(W)

"""
`elements(W::CoxeterGroup[,l])`

When  `l` is  not given  this works  only if  `W` is finite; it returns the
elements of `W` sorted by increasing Coxeter length. If the second argument
is an integer `l`, the elements of `W` of Coxeter length `l` are returned.

```julia_repl
julia> W=coxgroup(:G,2)
G₂

julia> e=elements(W,6)
1-element Vector{Perm{Int16}}:
 (1,7)(2,8)(3,9)(4,10)(5,11)(6,12)

julia> e[1]==longest(W)
true
```
"""
function Groups.elements(W::CoxeterGroup{T}, l::Int)::Vector{T} where T
  elts=get!(()->OrderedDict(0=>[one(W)]),W,:elements)::OrderedDict{Int,Vector{T}}
  get!(elts,l)do
  if ngens(W)==1 return l>1 ? T[] : gens(W) end
# if l==1 return elts[1]=gens(W) end
  H=get!(()->reflection_subgroup(W,maxpara(W)),W,:maxpara)::CoxeterGroup{T}
  rc=get!(()->[[one(W)]],W,:rc)::Vector{Vector{T}}
  while length(rc)<=l
    new=reducedfrom(H,W,rc[end])
    if isempty(new) break
    else push!(rc,new)
    end
  end
# println("l=$l W=$W H=$H rc=$rc")
  v=T[]
  for i in max(0,l+1-length(rc)):l, x in rc[1+l-i]
    append!(v,elements(H,i).*Ref(x))
  end
# if applicable(nref,W) # this does not reduce much time
#   N=nref(W)
#   if N-l>l && !haskey(elts,N-l) elts[N-l]=elts[l].*longest(W) end
# end
  v
  end
end

function Groups.elements(W::CoxeterGroup,I::AbstractVector{<:Integer})
  reduce(vcat,map(i->elements(W,i),I))
end

Groups.elements(W::CoxeterGroup)=elements(W,0:nref(W))

const Wtype=Vector{Int8}
function Groups.words(W::CoxeterGroup{T}, l::Integer)where T
  ww=get!(()->OrderedDict(0=>[Wtype([])]),W,:words)::OrderedDict{Int,Vector{Wtype}}
  if haskey(ww,l) return ww[l] end
  if ngens(W)==1
    if l!=1 return Vector{Int}[]
    else return [[1]]
    end
  end
  H=get!(()->reflection_subgroup(W,maxpara(W)),W,:maxpara)::CoxeterGroup{T}
  rc=get!(()->[[Wtype([])]],W,:rcwords)::Vector{Vector{Wtype}}
  while length(rc)<=l
    new=reducedfrom(H,W,(x->W(x...)).(rc[end]))
    if isempty(new) break
    else push!(rc,Wtype.(word.(Ref(W),new)))
    end
  end
  ww[l]=Wtype([])
  d=inclusiongens(H)
  for i in max(0,l+1-length(rc)):l
    e=words(H,i)
    for x in rc[1+l-i], w in e push!(ww[l],vcat(d[w],x)) end
#   somewhat slower variant
#   for x in rc[1+l-i] append!(ww[l],vcat.(e,Ref(x))) end
  end
  ww[l]
end

"""
`words(W::CoxeterGroup[,l::Integer])`

With  one argument this works only if `W` is finite; it returns the reduced
Coxeter  words  of  elements  of  `W`  by  increasing length. If the second
argument  is an integer `l`, only the  elements of length `l` are returned;
this works for infinite Coxeter groups.

```julia_repl
julia> W=coxgroup(:G,2)
G₂

julia> e=elements(W,6)
1-element Vector{Perm{Int16}}:
 (1,7)(2,8)(3,9)(4,10)(5,11)(6,12)

julia> e[1]==longest(W)
true
```
"""
function Groups.words(W::CoxeterGroup)
  reduce(vcat,map(i->words(W,i),0:nref(W)))
end

"""
`bruhatless(W, x, y)`

whether `x≤y` in the Bruhat order, for `x,y∈ W`. We have `x≤y` if a reduced
expression  for `x`  can be  extracted from  one for  `w`). See  [(5.9) and
(5.10) Humphreys1990](biblio.htm#Hum90) for properties of the Bruhat order.

```julia-repl
julia> W=coxgroup(:H,3)
H₃

julia> w=W(1,2,1,3);

julia> b=filter(x->bruhatless(W,x,w),elements(W));

julia> word.(Ref(W),b)
12-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [3]
 [1, 2]
 [2, 1]
 [1, 3]
 [2, 3]
 [1, 2, 1]
 [1, 2, 3]
 [2, 1, 3]
 [1, 2, 1, 3]
```
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
`bruhatless(W, y)`

returns  a vector  whose `i`-th  element is  the vector  of elements of `W`
smaller for the Bruhat order than `w` and of Coxeter length `i-1`. Thus the
first  element  of  the  returned  list  contains  only  `one(W)`  and  the
`length(W,w)`-th element contains only `w`.

```julia-repl
julia> W=CoxSym(3)
𝔖 ₃

julia> bruhatless(W,Perm(1,3))
4-element Vector{Vector{Perm{Int16}}}:
 [()]
 [(1,2), (2,3)]
 [(1,2,3), (1,3,2)]
 [(1,3)]
```

see also [`Poset`](@ref) for Coxeter groups.
"""
function bruhatless(W::CoxeterGroup,w)
  if w==one(W) return [[w]] end
  i=firstleftdescent(W,w)
  s=W(i)
  res=bruhatless(W,s*w)
  for j in 1:length(res)-1
    res[j+1]=union(res[j+1],Ref(s).*filter(x->!isleftdescent(W,x,i),res[j]))
  end
  push!(res,Ref(s).*filter(x->!isleftdescent(W,x,i),res[end]))
end

"""
`Poset(W::CoxeterGroup,w=longest(W))`

returns  as a poset the Bruhat interval `[1,w]`of `W`. If `w` is not given,
the whole Bruhat Poset of `W` is returned (`W` must then be finite).

```julia-repl
julia> W=CoxSym(3)
𝔖 ₃

julia> Poset(W)
.<1,2<21,12<121
```
The  above  poset  is  constructed  efficiently  by  constructing the Hasse
diagram, but it could be constructed naively as follows:
```julia-repl
julia> p=Poset((x,y)->bruhatless(W,x,y),elements(W))
()<(1,2),(2,3)<(1,3,2),(1,2,3)<(1,3)
```
The  output is not so nice, showing permutations instead of words. This can
be fixed by defining:
```julia-repl
julia> p.show_element=(io,x,n)->join(io,word(W,x.elements[n]));

julia> p
<1,2<12,21<121

julia> W=CoxSym(4)
𝔖 ₄

julia> Poset(W,W(1,3))
.<3,1<13
```
"""
function FinitePosets.Poset(W::CoxeterGroup,w=longest(W))
  if w==one(W)
    p=Poset(CPoset([Int[]]),[w],Dict(:action=>map(x->[0],gens(W)), :W=>W))
  else
  s=firstleftdescent(W,w)
  p=Poset(W,W(s)*w)
  l=length(p)
  # action: the Cayley graph: for generator i, action[i][w]=sw
  # where w and sw are represented by their index in :elements
  new=filter(k->iszero(p.action[s][k]),1:l)
  append!(p.elements,Ref(W(s)).*p.elements[new])
  append!(hasse(p),map(x->Int[],new))
  p.action[s][new]=l.+(1:length(new))
  for i in eachindex(gens(W))
    append!(p.action[i],i==s ? new : fill(0,length(new)))
  end
  for i in 1:length(new) push!(hasse(p)[new[i]],l+i) end
  for i in 1:l
    j=p.action[s][i]
    if j>i
      for h in p.action[s][hasse(p)[i]]
        if h>l push!(hasse(p)[j],h)
          k=findfirst(==(p.elements[h]/p.elements[j]),gens(W))
          if k!==nothing p.action[k][j]=h;p.action[k][h]=j end
        end
      end
    end
  end
  end
  p.show_element=(io,x,n)->(e=x.elements[n];print(io,isone(e) ? "." :
                                                 joindigits(word(W,e))))
  p
end

"""
`words(W::CoxeterGroup,w)`

returns  the list  of all  reduced expressions  of the  element `w`  of the
Coxeter group `W`.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> words(W,longest(W))
16-element Vector{Vector{Int64}}:
 [1, 2, 1, 3, 2, 1]
 [1, 2, 3, 1, 2, 1]
 [1, 2, 3, 2, 1, 2]
 [1, 3, 2, 1, 3, 2]
 [1, 3, 2, 3, 1, 2]
 [2, 1, 2, 3, 2, 1]
 [2, 1, 3, 2, 1, 3]
 [2, 1, 3, 2, 3, 1]
 [2, 3, 1, 2, 1, 3]
 [2, 3, 1, 2, 3, 1]
 [2, 3, 2, 1, 2, 3]
 [3, 1, 2, 1, 3, 2]
 [3, 1, 2, 3, 1, 2]
 [3, 2, 1, 2, 3, 2]
 [3, 2, 1, 3, 2, 3]
 [3, 2, 3, 1, 2, 3]
```
"""
function Groups.words(W::CoxeterGroup{T},w::T)where T
  l=leftdescents(W,w)
  if isempty(l) return [Int[]] end
  reduce(vcat,map(x->vcat.(Ref([x]),words(W,W(x)*w)),l))
end

"""
`inversions(W,w)`

Returns  the inversions of the element `w` of the finite Coxeter group `W`,
that  is, the list of the  indices of reflections `r` of `W` such that
`l(rw)<l(w)` where `l` is the Coxeter length.

```julia-repl
julia> W=coxgroup(:A,3)
A₃

julia> inversions(W,W(1,2,1))
3-element Vector{Int64}:
 1
 2
 4
```
"""
inversions(W::CoxeterGroup,w)=filter(x->isleftdescent(W,w,x),1:nref(W))
# assumes isleftdescent works for all reflections

function parabolic_category(W,I::AbstractVector{<:Integer})
   Category(collect(sort(I));action=(J,e)->sort(action.(Ref(W),J,e)))do J
    map(setdiff(1:ngens(W),J)) do i
      longest(W,J)*longest(W,push!(copy(J),i))
    end
  end
end

"""
`standard_parabolic_class(W,I)`

`I`  should be a  subset of `eachindex(gens(W))`.  The function returns the
list of such subsets `W`-conjugate to `I`.

```julia-repl
julia> CoxGroups.standard_parabolic_class(coxgroup(:E,8),[7,8])
7-element Vector{Vector{Int64}}:
 [7, 8]
 [6, 7]
 [5, 6]
 [4, 5]
 [2, 4]
 [3, 4]
 [1, 3]
```
"""
standard_parabolic_class(W,I::Vector{Int})=parabolic_category(W,I).obj

# representatives of parabolic classes
function PermRoot.parabolic_reps(W::CoxeterGroup,s)
  l=collect(combinations(1:ngens(W),s))
  orbits=[]
  while !isempty(l)
    o=standard_parabolic_class(W,l[1])
    push!(orbits,o)
    l=setdiff(l,o)
  end
  first.(orbits)
end

"""
`coxeter_matrix(m::AbstractMatrix)` or `coxmat`

returns  the  Coxeter  matrix  of  the  Coxeter group defined by the cartan
matrix `m`

```julia-repl
julia> C=cartan(:H,3)
3×3 Matrix{Cyc{Int64}}:
       2  ζ₅²+ζ₅³   0
 ζ₅²+ζ₅³        2  -1
       0       -1   2

julia> coxmat(C)
3×3 Matrix{Int64}:
 1  5  2
 5  1  3
 2  3  1
```
"""
function coxmat(m::AbstractMatrix)
  function find(c)
    if c in 0:4 return [2,3,4,6,0][Int(c)+1] end
    x=conductor(c)
    if c==2+E(x)+E(x,-1) return x
    elseif c==2+E(2x)+E(2x,-1) return 2x
    else error("not a Cartan matrix of a Coxeter group")
    end
  end
  res=Int.([i==j for i in axes(m,1), j in axes(m,2)])
  for i in 2:size(m,1), j in 1:i-1
    res[i,j]=res[j,i]=find(m[i,j]*m[j,i])
  end
  res
end

"""
`coxeter_matrix(W)` or `coxmat`

returns the Coxeter matrix of the Coxeter group `W`, that is the matrix `m`
whose  entry `m[i,j]` contains the order of `W(i)*W(j)` where `W(i)` is the
`i`-th  Coxeter generator of  `W`. An infinite  order is represented by the
entry `0`.

```julia-repl
julia> W=CoxSym(4)
𝔖 ₄

julia> coxmat(W)
3×3 Matrix{Int64}:
 1  3  2
 3  1  3
 2  3  1
```
"""
coxmat(W::CoxeterGroup)=coxmat(cartan(W))

const coxeter_matrix=coxmat

"""
`braid_relations(W)`

this  function returns the  relations which present  the braid group of the
reflection group `W`. These are homogeneous (both sides of the same length)
relations  between generators in bijection  with the generating reflections
of  `W`. A presentation  of `W` is  obtained by adding relations specifying
the order of the generators.

```julia-repl
julia> W=complex_reflection_group(29)
G₂₉

julia> braid_relations(W)
7-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 1], [2, 1, 2]]
 [[2, 4, 2], [4, 2, 4]]
 [[3, 4, 3], [4, 3, 4]]
 [[2, 3, 2, 3], [3, 2, 3, 2]]
 [[1, 3], [3, 1]]
 [[1, 4], [4, 1]]
 [[4, 3, 2, 4, 3, 2], [3, 2, 4, 3, 2, 4]]
```

each  relation  is  represented  as  a  pair  of lists, specifying that the
product  of the  generators according  to the  indices on  the left side is
equal  to the product according to the  indices on the right side. See also
`diagram`.
"""
function braid_relations(t::TypeIrred)
  r=if t.series==:ST getchev(t,:BraidRelations)
  else
    m=coxmat(cartan(t))
    p(i,j)=map(k->iszero(k%2) ? j : i,1:m[i,j])
    vcat(map(i->map(j->[p(i,j),p(j,i)],1:i-1),axes(m,1))...)
  end
  haskey(t,:indices) ? map(x->map(y->t.indices[y],x),r) : r
end

braid_relations(W::Group)=vcat(braid_relations.(refltype(W))...)

# construct Coxeter generators and rootinclusion of reflection subgroup
# defined by J by Dyer's method.
function dyer(W,J)
  refs=refls(W,J)
  refs=vcat(orbits(Group(refs),refs)...)
  inc=Int.(indexin(refs,refls(W)))
  gens=filter(t->count(s->isrightdescent(W,refls(W,t),s),inc)==1,inc)
  gens=sort(gens)
  (gens=gens,rootinclusion=vcat(gens,sort(setdiff(inc,gens))))
end

"""
`longest(W)`

If  `W` is  finite, returns  the unique  element of  maximal length  of the
Coxeter group `W`. May loop infinitely otherwise.

```julia-repl
julia> longest(CoxSym(4))
(1,4)(2,3)
```

`longest(W,I)`

returns  the longest element of the  parabolic subgroup of `W` generated by
the generating reflections of indices in `I`.

```julia-repl
julia> longest(CoxSym(4))
(1,4)(2,3)
```
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

abstract type FiniteCoxeterGroup{T} <: CoxeterGroup{T} end

#--------------------- CoxSym ---------------------------------
@GapObj struct CoxSym{T} <: FiniteCoxeterGroup{Perm{T}}
  G::PermGroup{T}
  inversions::Vector{Tuple{Int,Int}}
  refls::Vector{Perm{T}}
  d::UnitRange{Int}
end

#@forward CoxSym.G Base.iterate, Groups.gens, Base.one, PermGroups.orbits,
#  PermGroups.orbit

"""
  `Coxsym(n::Integer)` or `CoxSym(m:n)`

The  symmetric group on the  letters `1:n` (or if  a `m≤n` is given, on the
letters  `m:n`)  as  a  Coxeter  group.  The  coxeter  generators  are  the
`Perm(i,i+1)` for `i` in `m:n-1`.
```julia-repl
julia> W=CoxSym(3)
𝔖 ₃

julia> gens(W)
2-element Vector{Perm{Int16}}:
 (1,2)
 (2,3)

julia> e=elements(W)
6-element Vector{Perm{Int16}}:
 ()
 (1,2)
 (2,3)
 (1,3,2)
 (1,2,3)
 (1,3)

julia> length.(Ref(W),e) # length in the genrators of the elements
6-element Vector{Int64}:
 0
 1
 1
 2
 2
 3
```
"""
CoxSym(n::Int)=CoxSym(1:n)

function CoxSym(d::UnitRange)
  gens=map(i->Perm(i,i+1;degree=d.stop),d.start:d.stop-1)
  inversions=[(i,i+k) for k in 1:length(d)-1 for i in d.start:d.stop-k]
  refs=[Perm(r...;degree=d.stop) for r in inversions]
  append!(refs,refs)
  CoxSym(Group(gens),inversions,refs,d,Dict{Symbol,Any}())
end

function Base.show(io::IO, W::CoxSym)
  if W.d.start==1 n=string(W.d.stop) else n=string(W.d) end
  if hasdecor(io) printTeX(io,"\\mathfrak S_{$n}")
  else print(io,"CoxSym(n)")
  end
end

PermRoot.action(W::CoxSym,i,p)=i^p

PermRoot.refltype(W::CoxSym)=get!(W,:refltype)do
  [TypeIrred(Dict(:series=>:A,:indices=>collect(1:length(W.d)-1)))]
end

PermRoot.inclusiongens(W::CoxSym)=W.d
Groups.classreps(W::CoxSym)=map(x->W(x...),classinfo(W).classtext)
Perms.reflength(W::CoxSym,a)=reflength(a)
PermRoot.nref(W::CoxSym)=length(W.inversions)
function PermRoot.simple_reps(W::CoxSym)
  get!(W,:simple_reps)do
    W.unique_refls=collect(1:nref(W))
    fill(1,length(W.refls))
  end::Vector{Int}
end
PermRoot.simple_reps(W::CoxSym,i)=simple_reps(W)[i]
PermRoot.refls(W::CoxSym)=W.refls
PermRoot.refls(W::CoxSym,i)=W.refls[i]
PermRoot.rank(W::CoxSym)=ngens(W)

"""
`isleftdescent(W::CoxeterGroup,w,i)`

returns  `true`  if  and  only  if  the `i`-th generating reflection of the
Coxeter  group `W` is  in the left  descent set of  the element `w` of `W`,
that is iff `length(W,W(i)*w)<length(W,w)`.

```julia-repl
julia> W=CoxSym(3)
𝔖 ₃

julia> isleftdescent(W,Perm(1,2),1)
true
```
"""
function isleftdescent(W::CoxSym,w,i::Integer)
 j,k=W.inversions[i]
 j^w>k^w
end

degrees(W::CoxSym)=2:ngens(W)+1
"""
`cartan(W::CoxeterGroup)`  The Cartan matrix of `W`.
"""
PermRoot.cartan(W::CoxSym)=cartan(:A,ngens(W))

# for reflection_subgroups note the difference with Chevie:
# leftdescents, rightdescents, word  use indices in W and not in parent(W)
"""
`reflection_subgroup(W::CoxSym,I)`

The only reflection subgroups defined for `CoxSym(n)` are for `I=1:m` for `m≤n`
"""
function PermRoot.reflection_subgroup(W::CoxSym,I::AbstractVector{Int})
  if sort(I)!=minimum(I):maximum(I) error(I," should be a:b for some a,b") end
  if W.d.start+maximum(I)>W.d.stop error("range too large") end
  CoxSym(W.d.start.+(minimum(I)-1:maximum(I)))
end
#------------------------ MatCox ------------------------------

@GapObj struct MatCox{T}<:CoxeterGroup{Matrix{T}}
  gens::Vector{Matrix{T}}
end

Base.one(W::MatCox)=one(W(1))
PermRoot.cartan(W::MatCox)=W.cartan
isleftdescent(W::MatCox,w,i::Int)=real(sum(w[i,:]))<0

"""
`coxeter_group(m)` or `coxgroup(m)`

`m`  should be a square  matrix of real cyclotomic  numbers. It returns the
Coxeter  group  whose  Cartan  matrix  is  `m`.  This is a matrix group `W`
constructed  as  follows.  Let  `V`  be  a  real  vector space of dimension
`size(m,1)`, let `eᵢ` be the canonical basis of `V`. Then `W` is the matrix
group generated by the reflections `sᵢ(eⱼ)=eⱼ-mᵢⱼ eᵢ`.

```julia-repl
julia> W=coxgroup([2 -2;-2 2])
coxeter_group([2 -2; -2 2])
```

Above is a way to construct the affine Weyl group  `Ã₁`.
"""
function coxeter_group(C::Matrix{T})where T
  I=one(C)
  MatCox(reflectionmat.(eachrow(I),eachrow(C)),Dict{Symbol,Any}(:cartan=>C))
end

const coxgroup=coxeter_group

maxpara(W::MatCox)=1:ngens(W)-1

function PermRoot.reflection_subgroup(W::MatCox,I::AbstractVector{Int})
  if length(I)>0 n=maximum(I)
    if I!=1:n error(I," should be 1:n for some n") end
  else n=0 end
  MatCox(gens(W)[I],Dict{Symbol,Any}())
end

function Base.show(io::IO,W::MatCox)
  print(io,"coxeter_group(",cartan(W),")")
end

end
