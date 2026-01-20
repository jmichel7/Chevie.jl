"""
A  suitable  reference  for  the  general  theory of Coxeter groups is, for
example, [bou68; chapter 4].

A *Coxeter group* is a group which has the presentation
``W=⟨S|(st)^{m(s,t)}=1\\text{  for  }s,t∈  S⟩``  for some symmetric integer
matrix  `m(s,t)` called the  [`coxeter_matrix`](@ref), where `m(s,t)>1` for
`s≠t`  and `m(s,s)=1`; `m(s,t)=∞`  is allowed meaning  there is no relation
between  `s` and  `t`. It  is true  (but a  non-trivial theorem)  that in a
Coxeter  group the order of `st` is  exactly `m(s,t)`, thus a Coxeter group
is  the same as a *Coxeter  system*, that is a pair  `(W,S)` of a group `W`
and  a  set  `S⊂W`  of  involutions,  such  that  the group is presented by
generators  `S` and  relations describing  the order  of the product of two
elements of `S`.

A   Coxeter   group   has   a   natural   representation,  its  *reflection
representation*, on a real vector space `V` of dimension `length(S)` (which
is  the  *Coxeter  rank*  of  W),  where  each  element  of  `S`  acts as a
reflection;  the faithfulness of this representation (a theorem of Tits) is
the main argument to prove that the order of `st` is exactly `m(s,t)`. This
representation  is defined as follows on a  space `V` with basis `{eₛ}` for
`s∈  S`. The *Cartan  matrix* associated to  the Coxeter matrix `m(s,t)` is
the matrix `C` with entries `C(s,t)=-2cos(π/m(s,t))`; we set `C(s,t)=-2` if
`m(s,t)=∞`. Then the action of `s∈ S` on `V` is given by `s(eₜ)=eₜ-C(s,t)eₛ`.

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
[`coxeter_symmetric_group`](@ref) below), so this suggests how to deal with
Coxeter groups generically:

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
julia> W=coxsym(4)
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

The  only  Coxeter  group  constructors  implemented  in  this  module  are
[`coxeter_symmetric_group`](@ref)  and  [`coxeter_group`](@ref);  the  last
constructor  takes  a  Cartan  matrix  and builds the corresponding Coxeter
group  as a matrix  group. The module  [`Weyl`](@ref) defines other methods
for  `coxgroup` building  a finite  Coxeter group  as a  permutation group,
given its type.
"""
module CoxGroups

export bruhatless, bruhatPoset, CoxeterGroup, firstleftdescent, leftdescents,
  isleftdescent, isrightdescent,
  longest, braid_relations, coxmat, 
  coxeter_symmetric_group, coxsym, CoxSym,
  coxeter_hyperoctaedral_group, coxhyp, CoxHyp,
  coxeter_group, coxgroup,
  standard_parabolic_class, inversions, degrees, FiniteCoxeterGroup

using ..Chevie
#-------------------------- Coxeter groups
abstract type CoxeterGroup{T}<:Group{T} end

"""
`firstleftdescent(W,w)`

returns the index in `gens(W)` of the first element of the left descent set
of  `w` --- that is, the first  `i` such that if `s=W(i)` then `l(sw)<l(w).
It returns `nothing` for `one(W)`.

```julia-repl
julia> W=coxsym(3)
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
julia> W=coxsym(3)
𝔖 ₃

julia> leftdescents(W,Perm(1,3))
2-element Vector{Int64}:
 1
 2
```
"""
function leftdescents(W::CoxeterGroup,w)
  filter(i->isleftdescent(W,w,i),eachindex(gens(W)))
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
julia> W=coxsym(4)
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
function Groups.elements(W::CoxeterGroup{T}, l::Integer)::Vector{T} where T
  elts=get!(()->OrderedDict(0=>[one(W)]),W,:elements)::OrderedDict{Int,Vector{T}}
  get!(elts,l)do
  if ngens(W)==1 return l>1 ? T[] : gens(W) end
  if isfinite(W)
    H=get!(()->reflection_subgroup(W,maxpara(W)),W,:maxpara)::CoxeterGroup{T}
    rc=get!(()->[[one(W)]],W,:rc)::Vector{Vector{T}}
    while length(rc)<=l
      new=reducedfrom(H,W,rc[end])
      if isempty(new) break
      else push!(rc,new)
      end
    end
  # @show l,W,H,rc
    v=T[]
    for i in max(0,l+1-length(rc)):l, x in rc[1+l-i]
      append!(v,elements(H,i).*Ref(x))
    end
  # if applicable(nref,W) # this does not reduce much time
  #   N=nref(W)
  #   if N-l>l && !haskey(elts,N-l) elts[N-l]=elts[l].*longest(W) end
  # end
    v
  else
    if l==1 return gens(W) end
    v=elements(W,l-1)
    res=Set(empty(v))
    for e in v, i in 1:ngens(W)
      if !isleftdescent(W,e,i) push!(res,W(i)*e) end
    end
    collect(res)
  end
  end
end

function Groups.elements(W::CoxeterGroup,I::AbstractVector{<:Integer})
  reduce(vcat,elements.(Ref(W),I))
end

function Groups.elements(W::CoxeterGroup)
  if isfinite(W) return elements(W,0:nref(W))
  else error(W," is infinite")
  end
end

PermRoot.nref(W::CoxeterGroup)=length(W,longest(W))

const wordtype=Vector{Int8}
function Groups.words(W::CoxeterGroup{T}, l::Integer)where T
  if !isfinite(W) return word.(Ref(W),elements(W,l)) end
  ww=get!(W,:words)do
    OrderedDict(0=>[wordtype([])])
  end::OrderedDict{Int,Vector{wordtype}}
  if haskey(ww,l) return ww[l] end
  if ngens(W)==1
    if l!=1 return Vector{Int}[]
    else return [[1]]
    end
  end
  H=get!(()->reflection_subgroup(W,maxpara(W))::CoxeterGroup{T},W,:maxpara)
  rc=get!(()->[[wordtype([])]]::Vector{Vector{wordtype}},W,:rcwords)
  while length(rc)<=l
    new=reducedfrom(H,W,(x->W(x...)).(rc[end]))
    if isempty(new) break
    else push!(rc,wordtype.(word.(Ref(W),new)))
    end
  end
  ww[l]=wordtype([])
  d=maxpara(W) # inclusiongens(H,W)
  for i in max(0,l+1-length(rc)):l
    e=words(H,i)
    for x in rc[1+l-i], w in e push!(ww[l],vcat(d[w],x)) end
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
  reduce(vcat,words.(Ref(W),0:nref(W)))
end

"""
`bruhatless(W, x, y)`

whether `x≤y` in the Bruhat order, for `x,y∈ W`. We have `x≤y` if a reduced
expression for `x` can be extracted from one for `w`. See [hum90; (5.9) and
(5.10)](@cite) for properties of the Bruhat order.

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
  lx=length(W,x)
  d=length(W,y)-lx
  while d>0
    if iszero(lx) return true end
    i=firstleftdescent(W,y)
    s=W(i)
    if isleftdescent(W,x,i) x=s*x;lx-=1
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
julia> W=coxsym(3)
𝔖 ₃

julia> bruhatless(W,Perm(1,3))
4-element Vector{Vector{Perm{Int16}}}:
 [()]
 [(1,2), (2,3)]
 [(1,2,3), (1,3,2)]
 [(1,3)]
```

see also [`bruhatPoset`](@ref) for Coxeter groups.
"""
function bruhatless(W::CoxeterGroup,w)
  if isone(w) return [[w]] end
  i=firstleftdescent(W,w)
  s=W(i)
  res=bruhatless(W,s*w)
  for j in 1:length(res)-1
    res[j+1]=union(res[j+1],Ref(s).*filter(x->!isleftdescent(W,x,i),res[j]))
  end
  push!(res,Ref(s).*filter(x->!isleftdescent(W,x,i),res[end]))
end

"""
`bruhatPoset(W::CoxeterGroup,w=longest(W))`

returns  as a poset the Bruhat interval `[1,w]`of `W`. If `w` is not given,
the whole Bruhat Poset of `W` is returned (`W` must then be finite).

```julia-repl
julia> W=coxsym(3)
𝔖 ₃

julia> bruhatPoset(W)
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

julia> W=coxsym(4)
𝔖 ₄

julia> bruhatPoset(W,W(1,3))
.<3,1<13
```
"""
function bruhatPoset(W::CoxeterGroup,w=longest(W))
  if isone(w)
    p=Poset(CPoset([Int[]]),[w],Dict(:action=>map(x->[0],gens(W)), :W=>W))
  else
  s=firstleftdescent(W,w)
  p=bruhatPoset(W,W(s)*w)
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
  if isone(w) return [Int[]] end
  reduce(vcat,map(x->pushfirst!.(words(W,W(x)*w),x),leftdescents(W,w)))
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
   Category(sort(I);action=(J,e)->sort!(action.(Ref(W),J,e)))do J
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
  l=combinations(1:ngens(W),s)
  res=Vector{Int}[]
  while !isempty(l)
    o=standard_parabolic_class(W,first(l))
    push!(res,first(o))
    l=setdiff(l,o)
  end
  res
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
function FinitePosets.coxeter_matrix(m::AbstractMatrix)
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
julia> W=coxsym(4)
𝔖 ₄

julia> coxmat(W)
3×3 Matrix{Int64}:
 1  3  2
 3  1  3
 2  3  1
```
"""
FinitePosets.coxeter_matrix(W::CoxeterGroup)=coxeter_matrix(cartan(W))

const coxmat=FinitePosets.coxeter_matrix

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
7-element Vector{Tuple{Vector{Int64}, Vector{Int64}}}:
 ([1, 2, 1], [2, 1, 2])
 ([2, 4, 2], [4, 2, 4])
 ([3, 4, 3], [4, 3, 4])
 ([2, 3, 2, 3], [3, 2, 3, 2])
 ([1, 3], [3, 1])
 ([1, 4], [4, 1])
 ([4, 3, 2, 4, 3, 2], [3, 2, 4, 3, 2, 4])
```

each  relation  is  represented  as  a  pair  of lists, specifying that the
product  of the  generators according  to the  indices on  the left side is
equal  to the product according to the  indices on the right side. See also
`diagram`.
"""
function braid_relations(t::TypeIrred)
  r=if t.series==:ST Tuple.(chevieget(t,:BraidRelations))
  else
    m=coxmat(cartan(t))
    p(i,j)=map(k->iseven(k) ? j : i,1:m[i,j])
    vcat(map(i->map(j->(p(i,j),p(j,i)),1:i-1),axes(m,1))...)
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
julia> longest(coxsym(4))
(1,4)(2,3)
```

`longest(W,I)`

returns  the longest element of the  parabolic subgroup of `W` generated by
the generating reflections of indices in `I`.

```julia-repl
julia> longest(coxsym(4))
(1,4)(2,3)
```
"""
function longest(W::CoxeterGroup,I::AbstractVector{<:Integer}=eachindex(gens(W)))
  if length(I)==ngens(W) && !isfinite(W) error(W," must be finite") end
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

# FiniteCoxeterGroup  should  be  derived  from PermGroup and CoxeterGroup.
# Since  such inheritance  is impossible  in Julia,  we have  to choose. We
# derive  it  from  CoxeterGroup  and  represent  the  other inheritance by
# composition.  Thus, implementations  of FiniteCoxeterGroups  have a field
# .G, a PermGroup, to which we forward the following methods.
@forward FiniteCoxeterGroup.G Base.iterate, Base.one,
 Groups.gens, Groups.conjugacy_classes,
 PermGroups.classreps

Groups.transporting_elt(W::FiniteCoxeterGroup,x,y,F::Function)=
  transporting_elt(W.G,x,y,F)
Base.hash(W::FiniteCoxeterGroup,i::UInt64)=hash(W.G,i)
PermGroups.orbit(W::FiniteCoxeterGroup,p,F::Function)=orbit(W.G,p,F)
PermGroups.orbits(W::FiniteCoxeterGroup,v,F::Function=^)=orbits(W.G,v,F)
PermGroups.orbits(W::FiniteCoxeterGroup)=orbits(W.G)
PermGroups.stabilizer(W::FiniteCoxeterGroup,p,F::Function)=stabilizer(W.G,p,F)
PermGroups.stabilizer(W::FiniteCoxeterGroup,p,::typeof(ontuples))=stabilizer(W.G,p,ontuples)

Base.isfinite(W::FiniteCoxeterGroup)=true
#--------------------- CoxSym ---------------------------------
@GapObj struct CoxSym{T} <: FiniteCoxeterGroup{Perm{T}}
  G::PermGroup{T}
  inversions::Vector{Tuple{Int,Int}}
  refls::Vector{Perm{T}}
  d::Vector{Int} # permuted ints
  n::Int  # degree of elements
end

"""
`coxeter_symmetric_group(n::Integer)` or `coxeter_symmetric_group(v::AbstractVector{<:Integer})` or
`coxsym(n)` or `coxsym(v::AbstractVector{<:Integer})`

The  symmetric group  on the  numbers `1:n`  (or if  a `v` is given, on the
numbers  in `v`, assumed to be positive and increasing) as a Coxeter group.
The Coxeter generators are the `Perm(v[i],v[i+1])` for `i` in `1:length(v)-1`.
```julia-repl
julia> W=coxsym(3)
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

julia> length.(Ref(W),e) # length in the generators of the elements
6-element Vector{Int64}:
 0
 1
 1
 2
 2
 3
```
"""
coxeter_symmetric_group(n::Int)=coxeter_symmetric_group(1:n)

const coxsym=coxeter_symmetric_group

function coxsym(d::AbstractVector{<:Integer},n=maximum(d))
  if !(issorted(d) && d[1]>=1) error(d," should be increasing and positive") end
  gens=map(i->Perm(d[i],d[i+1];degree=n),1:length(d)-1)
  inversions=[(d[i],d[i+k]) for k in 1:length(d)-1 for i in 1:length(d)-k]
  refs=[Perm(r...;degree=n) for r in inversions]
  nref=length(refs)
  append!(refs,refs)
  G=Group(gens,one(prod(gens)))
  G.classreps=map(x->G(chevieget(:A,:WordClass)(x)...),partitions(length(d)))
  CoxSym(G,inversions,refs,collect(d),n,
   Dict{Symbol,Any}(:unique_refls=>collect(1:nref),:simple_reps=>fill(1,2nref)))
end

function Base.show(io::IO, W::CoxSym)
  if W.d==W.d[1]:W.d[end] 
    n= W.d[1]==1 ? string(W.d[end]) : string(W.d[1]:W.d[end])
  else n= hasdecor(io) ? joindigits(W.d) : string(W.d) end
  if hasdecor(io) printTeX(io,"\\mathfrak S_{$n}")
  else print(io,"coxsym($n)")
  end
end

PermRoot.action(W::CoxSym,i,p)=i^p

PermRoot.refltype(W::CoxSym)=get!(W,:refltype)do
  [TypeIrred(;series=:A,indices=collect(1:length(W.d)-1))]
end

PermRoot.inclusiongens(W::CoxSym)=W.d[1:end-1]
Perms.reflength(W::CoxSym,a)=reflength(a)
PermRoot.nref(W::CoxSym)=length(W.inversions)
PermRoot.refls(W::CoxSym)=W.refls
PermRoot.simple_reps(W::CoxSym)=W.simple_reps
Symbols.rank(W::CoxSym)=ngens(W)+1
PermRoot.reflrep(W::CoxSym,w::Perm)=Matrix(w,W.n)
PermRoot.reflrep(W::CoxSym,i::Integer)=Matrix(W(i),W.n)
PermRoot.reflrep(W::CoxSym)=reflrep.(Ref(W),1:ngens(W))

"""
`isleftdescent(W::CoxeterGroup,w,i::Integer)`

returns  `true` iff the  `i`-th generating reflection  of the Coxeter group
`W`  is in  the left  descent set  of the  element `w`  of `W`, that is iff
`length(W,W(i)*w)<length(W,w)`.
```julia-repl
julia> W=coxsym(3)
𝔖 ₃

julia> isleftdescent(W,Perm(1,2),1)
true
```
"""
isleftdescent(W::CoxSym,w,i::Integer)= >(W.inversions[i].^w...)

"""
`isrightdescent(W::CoxeterGroup,w,i::Integer)`

returns  `true` iff the  `i`-th generating reflection  of the Coxeter group
`W`  is in the  right descent set  of the element  `w` of `W`,  that is iff
`length(W,w*W(i))<length(W,w)`.
```julia-repl
julia> W=coxsym(3)
𝔖 ₃

julia> isrightdescent(W,Perm(1,2),1)
true
```
"""
isrightdescent(W::CoxSym,w,i::Integer)= >(preimage.(W.inversions[i],w)...)

"""
`cartan(W::CoxeterGroup)`  The Cartan matrix of `W`.
"""
PermRoot.cartan(W::CoxSym)=cartan(:A,ngens(W))

# for reflection_subgroups note the difference with Chevie:
# leftdescents, rightdescents, word  use indices in W and not in parent(W)
"""
`reflection_subgroup(W::CoxSym,I)`

The  only reflection subgroups defined for  `coxsym(n)` are for `I=a:b` for
`1≤a≤b≤n`
"""
function PermRoot.reflection_subgroup(W::CoxSym,I::AbstractVector{Int})
  if !(issorted(I) && I[end]<=length(W.d)-1 && I==I[1]:I[end])
    error(I," should be a condecutive subset of ",1:length(W.d))
  end
  coxsym(W.d[I[1]:I[end]+1],W.n)
end
#------------ Example II: HyperOctaedral groups as Coxeter groups

@GapObj struct CoxHyp{T} <: FiniteCoxeterGroup{SPerm{T}}
  G::SPermGroup{T}
  roots::Vector{Tuple{T,T}}
  n::Int
end

Groups.gens(W::CoxHyp)=gens(W.G)

"""
`coxeter_hyperoctaedral_group(n)`  or `coxhyp(n)`

The  Hyperoctaedral group (the group of all signed permutations of ±1,…,±n)
as   a  Coxeter   group  of   type  `Bₙ`,   with  generators  `(1,-1)`  and
`(i,i+1)(-i,-i-1)`.

```julia-repl
julia> elements(coxhyp(2))
8-element Vector{SPerm{Int8}}:
 ()
 (1,2)
 (1,-1)
 (1,2,-1,-2)
 (1,-2,-1,2)
 (2,-2)
 (1,-2)
 (1,-1)(2,-2)
```
"""
function coxeter_hyperoctaedral_group(n::Int)
  roots=Tuple{Int8,Int8}[(1,-1)]
  append!(roots,map(i->(i-1,i),2:n))
  append!(roots,map(i->(i,-i),2:n))
  for i in 2:n-1 append!(roots,map(j->(j,j+i),1:n-i)) end
  for i in 1:n-1 append!(roots,map(j->(j,-j-i),1:n-i)) end
  append!(roots,roots)
  G=hyperoctaedral_group(n)
  G.classreps=map(x->G(chevieget(:B,:WordClass)(x)...),partition_tuples(n,2))
  CoxHyp{Int8}(G,roots,n,Dict{Symbol,Any}())
end

const coxhyp=coxeter_hyperoctaedral_group

Base.show(io::IO, W::CoxHyp)=print(io,"coxhyp($(W.n))")

PermRoot.action(W::CoxHyp,i,p)=abs(i^p)

PermRoot.refltype(W::CoxHyp)=[TypeIrred(;series=:B,indices=collect(1:W.n))]

CoxGroups.nref(W::CoxHyp)=div(length(W.roots),2)

# use that "roots" correspond to roots of so_{2n+1}. 
# If eᵢ are the basis vectors (i,-i) <-> eᵢ, (i,j) <-> eⱼ-eᵢ, (i,-j) <-> eᵢ+eⱼ
function CoxGroups.isleftdescent(W::CoxHyp,w,i)
  a,b=W.roots[i]
  aw=a^w
  if a==-b return aw<0 end
  bw=b^w
  if b>0 aw=-aw else bw=-bw end
  if aw>0 && bw>0 return false end
  if aw<0 && bw<0 return true end
  max(aw,bw)<-min(aw,bw)
end

Symbols.rank(W::CoxHyp)=W.n
CoxGroups.maxpara(W::CoxHyp)=1:W.n-1
PermRoot.inclusiongens(W::CoxHyp)=1:W.n

function PermRoot.refls(W::CoxHyp{T})where T
  get!(W,:refls)do
    refs=map(t->SPerm{T}(t...),W.roots)
    W.simple_reps=map(x->x[1]==-x[2] ? 1 : 2,W.roots)
    W.unique_refls=collect(1:nref(W))
    refs
  end::Vector{SPerm{T}}
end

PermRoot.simple_reps(W::CoxHyp)=getp(refls,W,:simple_reps)
PermRoot.unique_refls(W::CoxHyp)=getp(refls,W,:unique_refls)

function Perms.reflength(W::CoxHyp,w)
  sym=nsym=0
  for x in cycles(w)
    if sort(x)==sort(-x) sym+=length(x)
    else nsym+=length(x)-1
    end
  end
  div(sym,2)+nsym
end

PermRoot.reflrep(W::CoxHyp,w::SPerm)=Matrix(w,W.n)
PermRoot.reflrep(W::CoxHyp,i::Integer)=
  reflrep(W,i<=ngens(W) ? gens(W)[i] : refls(W)[i])
PermRoot.reflrep(W::CoxHyp)=reflrep.(Ref(W),1:ngens(W))

"""
`reflection_subgroup(W::CoxHyp,I)`
is defined only for `I=1:m` for `m≤ngens(W)`.
"""
function PermRoot.reflection_subgroup(W::CoxHyp,I::AbstractVector{Int})
  if I!=1:length(I) error(I," should be 1:n for some n") end
  CoxHyp(Group(gens(W)[I]),filter(t->abs(t[2]!=W.n),W.roots),
                                    length(I),Dict{Symbol,Any}())
end
PermRoot.cartan(W::CoxHyp)=cartan(:B,W.n)
#------------------------ MatCox ------------------------------

@GapObj struct MatCox{T}<:CoxeterGroup{Matrix{T}}
  gens::Vector{Matrix{T}}
end

Base.one(W::MatCox)=one(W(1))
PermRoot.cartan(W::MatCox)=W.cartan
isleftdescent(W::MatCox,w,i::Int)=real(sum(w[i,:]))<0

"""
`coxeter_group(C)` or `coxgroup(C)`

`C`  should be  a square  matrix of  integers, rationals or real cyclotomic
numbers.  The function returns the Coxeter group whose Cartan matrix is `C`
as  a matrix group with generators the  following reflections: let `V` be a
real  vector space of  dimension `size(C,1)` and  let `eᵢ` be the canonical
basis of `V`. Then the generators are the reflections `sᵢ(eⱼ)=eⱼ-Cᵢⱼ eᵢ`.

```julia-repl
julia> W=coxgroup([2 -2;-2 2]) # the affine Weyl group  `Ã₁`.
coxeter_group([2 -2; -2 2])

julia> gens(W) # the matrix generators
2-element Vector{Matrix{Int64}}:
 [-1 0; 2 1]
 [1 2; 0 -1]
```
"""
function coxeter_group(C::Matrix{T})where T
  I=one(C)
  MatCox(reflectionMatrix.(eachrow(I),eachrow(C)),Dict{Symbol,Any}(:cartan=>C))
end

function Base.isfinite(W::MatCox)
  get!(W,:isfinite)do
    !any(isnothing,Weyl.fincox_refltype(cartan(W)))
  end::Bool
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
