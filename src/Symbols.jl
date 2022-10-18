"""
The  combinatorial objects  in this  module are  *partitions*, *β-sets* and
*symbols*.

A  partition is a  non-increasing list of  positive integers `p₁≥p₂≥…pₙ>0`,
represented as a `Vector{Int}`.

A  *β-set* is a strictly increasing `Vector` of nonnegative integers, up to
*shift*,  the  equivalence  relation  generated  by the *elementary shifts*
`[b₁,…,bₙ]∼[0,1+b₁,…,1+bₙ]`.  An equivalence  class has  exactly one member
which does not contain `0`: it is called a normalized β-set.

To  a  partition  `p₁≥p₂≥…pₙ>0`  is  associated  a  β-set, whose normalized
representative   is   `pₙ,pₙ₋₁+1,…,p₁+n-1`.   Conversely,   to  each  β-set
`b₁<b₂<…<bₙ` is associated the partition `bₙ-n+1≥…≥b₂-1≥b₁`, which may have
some trailing zeros if starting from a non-normalized representative.

`2`-symbols  where introduced by [Lusztig1977] and more general `e`-symbols
by  [Malle1995](biblio.htm#Mal95). An `e`-symbol  is a vector `S=[S₁,…,Sₑ]`
of   β-sets,  taken  modulo  the  equivalence  relation  generated  by  the
simultaneous  elementary shift of all β-sets, and by cyclic permutations of
`S`;  in the particular case where `e=2`,  `S` is thus an unordered pair of
β-sets.  `S` is a *normalized symbol* if  `0` is not in the intersection of
the   `Sᵢ`;  equivalent   normalized  symbols   are  equivalent  by  cyclic
permutation.  The *content* of `S`  is `mod(c,e)` where `c=sum(length(S))`;
it  is an invariant  of the symbol,  as well as  the *rank*, defined for an
`e`-symbol  as `sum(sum,S)-div((c-1)*(c-e+1),2*e)`. Invariant  by shift but
not cyclic permutation is the *shape* `s-minimum(s)` where `s=map(length,S)`.

When  `e=2` we  choose representatives  of the  symbols `[S₁,S₂]` such that
`length(S₁)≥length(S₂)`,  so the shape is `[d,0]` for some `d≥0` called the
*defect*  of the  symbol; the  content is  equal to `mod(d,2)`. For symbols
`[S₁,S₂]` with `length(S₁)==length(S₂)` we choose representatives such that
`P₁≤P₂`  lexicographically where  `P₁,P₂` are  the partitions associated to
`S₁,S₂`.

Partitions  and  pairs  of  partitions  parametrize  characters of the Weyl
groups  of classical types, and tuples of partitions parametrize characters
of  imprimitive complex reflection  groups. 2-Symbols parametrize unipotent
characters  of  classical  Chevalley  groups,  and more general `e`-symbols
parametrize   unipotent  characters  of   Spetses  associated  to  spetsial
imprimitive  complex reflection groups. The rank of a symbol is equal tothe
semi-simple rank of the corresponding Chevalley group or Spets.

Symbols of rank `n` and defect `0` parametrize characters of the Weyl group
of  type  `Dₙ`,  and  symbols  of  rank  `n`  and  defect  divisible by `4`
parameterize  unipotent characters of split  orthogonal groups of dimension
`2n`.  Symbols of  rank `n`  and defect`≡2  (mod 4)` parameterize unipotent
characters  of non-split  orthogonal groups  of dimension  `2n`. Symbols of
rank  `n` and defect `1`  parametrize characters of the  Weyl group of type
`Bₙ`,  and  symbols  of  rank  `n`  and  odd  defect  parametrize unipotent
characters  of symplectic groups of dimension  `2n` or orthogonal groups of
dimension `2n+1`.

`e`-symbols  of rank `n` and  content `1` parameterize unipotent characters
of  `G(e,1,n)`. Those of  content `0` parameterize  unipotent characters of
`G(e,e,n)`.  The  principal  series  (in  bijection  with characters of the
reflection  group)  is  parametrized  by  symbols  of shape `[1,0,…,0]` for
`G(e,1,n)` and `[0,…,0]` for `G(e,e,n)`.
"""
module Symbols
using ..Util: joindigits
using ..Combinat: arrangements, partition_tuples, allequal, collectby
if VERSION<=v"1.7.6"
using ..Combinat: allequal
end
using ..CyclotomicNumbers: E
using ..CycPols: CycPol, subs
using LaurentPolynomials
export shiftβ, βset, partβ, symbol_partition_tuple,
valuation_gendeg_symbol,      degree_gendeg_symbol,      degree_fegsymbol,
valuation_fegsymbol,   defectsymbol,   fullsymbol,   ranksymbol,  symbols,
fegsymbol, stringsymbol, XSP, PartitionTupleToString, gendeg_symbol, 
EnnolaSymbol

"""
`PartitionTupleToString(tuple)`

converts  the partition tuple `tuple` to  a string where the partitions are
separated by a dot.

```julia-repl
julia> d=partition_tuples(3,2)
10-element Vector{Vector{Vector{Int64}}}:
 [[1, 1, 1], []]
 [[1, 1], [1]]
 [[1], [1, 1]]
 [[], [1, 1, 1]]
 [[2, 1], []]
 [[1], [2]]
 [[2], [1]]
 [[], [2, 1]]
 [[3], []]
 [[], [3]]

julia> PartitionTupleToString.(d)
10-element Vector{String}:
 "111."
 "11.1"
 "1.11"
 ".111"
 "21."
 "1.2"
 "2.1"
 ".21"
 "3."
 ".3"
```
"""
function PartitionTupleToString(n,a=Dict())
  if n[end] isa Vector return join(joindigits.(n),".") end
  r=repr(E(n[end-1],n[end]),context=:limit=>true)
  if r=="1" r="+" end
  if r=="-1" r="-" end
  join(joindigits.(n[1:end-2]),".")*r
end

"""
`shiftβ( β, n)` shift the β-set `β` by `n`

```julia-repl
julia> shiftβ([2,3],2)
4-element Vector{Int64}:
 0
 1
 4
 5

julia> shiftβ([0,1,4,5],-2)
2-element Vector{Int64}:
 2
 3

```
"""
function shiftβ(β,n=1)
  if n>=0 return [0:n-1;β .+ n]
  elseif β[1:-n]!= 0:-n-1 error("Cannot shift $β by $n\n")
  else return β[1-n:end].+n
  end
end

"""
`βset(p)` normalized β-set of partition `p`

```julia-repl
julia> βset([3,3,1])
3-element Vector{Int64}:
 1
 4
 5
```
"""
function βset(p,s=0)
# shorter code if we don't care about allocations:
#  p=vcat(p,fill(0,s))
#  reverse(p).+(0:length(p)-1)
  if !iszero(s)
    q=fill(0,length(p)+s)
    @inbounds for i in eachindex(p) q[i]=p[i] end
  else q=p
  end
  if isempty(q) return q end
  p=collect(0:length(q)-1)
  @inbounds for i in eachindex(p) p[i]+=q[1+length(q)-i] end
  p
end

"""
`partβ(β)` partition defined by β-set `β`

```julia-repl
julia> partβ([0,4,5])
2-element Vector{Int64}:
 3
 3
```
"""
partβ(β)=filter(!iszero,reverse(β.-(0:length(β)-1)))

"""
`symbol_partition_tuple(p, s)` symbol of shape `s` for partition tuple `p`.

In  the general case, `s` is a `Vector{Int}`  of same length as `p` and the
`i`-th  element of the result is the β-set for `pᵢ` shifted to be of length
`sᵢ` (the minimal integer which makes this possible is added to `s`).

When  `s` is  a positive  integer it  is interpreted  as `[s,0,0,…]`  and a
negative  integer is interpreted  as `[0,-s,-s,…]` so  when `p` is a double
partition  one gets the  symbol of defect  `s` associated to  `p`; as other
uses  the  unipotent  symbol  for  a  character  of the principal series of
`G(e,1,r)`   parameterized   by   an   `e`-tuple   `p`   of  partitions  is
`symbol_partition_tuple(p,1)` and for `G(e,e,r)` the similar computation is
`symbol_partition_tuple(p,0)`  (the function handles coded periodic `p` for
`G(e,e,r)`).

```julia-repl
julia> symbol_partition_tuple([[1,2],[1]],1)
2-element Vector{Vector{Int64}}:
 [2, 2]
 [1]

julia> symbol_partition_tuple([[1,2],[1]],0)
2-element Vector{Vector{Int64}}:
 [2, 2]
 [0, 2]

julia> symbol_partition_tuple([[1,2],[1]],-1)
2-element Vector{Vector{Int64}}:
 [2, 2]
 [0, 1, 3]
```
"""
function symbol_partition_tuple(p,S)
  if p[end] isa Number
    l=1:length(p)-2
    return vcat(symbol_partition_tuple(p[l],S isa Integer ? S : S[l]),
                         p[end-1:end])
  end
  if S isa Integer
    if S<0 s=fill(-S,length(p));s[1]=0
    else   s=fill(0,length(p));s[1]=S
    end
  else s=copy(S)
  end
  s.-=length.(p)
  s.-=minimum(s)
  βset.(p,s)
end

function fullsymbol(S)::Vector{Vector{Int}}
  if isempty(S) || S[end] isa AbstractVector return S end
  reduce(vcat,map(i->map(copy,S[1:end-2]),1:S[end-1]))
end

"""
`ranksymbol(S)` rank of symbol `S`.

```julia-repl
julia> ranksymbol([[1,5,6],[1,2]])
11
```
"""
function ranksymbol(s)
  if isempty(s) return 0 end
  s=fullsymbol(s)
  ss=sum(length,s)
  e=length(s)
  sum(sum,s)-div((ss-1)*(ss-e+1),2*e)
end

"""
`valuation_gendeg_symbol(s)`

the   valuation  of   the  generic   degree  of   the  unipotent  character
parameterized by the `e`-symbol `s`.

```julia-repl
julia> valuation_gendeg_symbol([[1,5,6],[1,2]])
13
```
"""
function valuation_gendeg_symbol(p)
  p=fullsymbol(p)
  e=length(p)
  p=sort!(reduce(vcat,p))
  m=length(p)
  transpose(p)*(m-1:-1:0)-div(m*(m-e)*(2*m-3-e),12*e)
end

"""
`degree_gendeg_symbol(s)`

the  degree of the generic degree  of the unipotent character parameterized
by the `e`-symbol `s`.

```julia-repl
julia> degree_gendeg_symbol([[1,5,6],[1,2]])
91
```
"""
function degree_gendeg_symbol(p)
  p=fullsymbol(p)
  r=ranksymbol(p)
  e=length(p)
  p=sort!(reduce(vcat,p))
  m=length(p)
  if mod(m,e)==1 r=div(e*r*(r+1),2)
  else           r+= div(e*r*(r-1),2)
  end
  r+transpose(p)*(0:m-1)-sum(x->div(e*x*(x+1),2),p)-div(m*(m-e)*(2*m-3-e),12*e)
end

"""
`defectsymbol(s)'

For an `e`-symbol `[S₁,S₂,…,Sₑ]` returns `length(S₁)-length(S₂)`.

```julia-repl
julia> defectsymbol([[1,5,6],[1,2]])
1
```
"""
function defectsymbol(S)
 S[end] isa AbstractVector && length(S)>1 ? length(S[1])-length(S[2]) : 0
end

"""
`degree_fegsymbol(s)`

the  degree  of  the  fake  degree  of  the  character parameterized by the
`e`-symbol `s`.

```julia-repl
julia> degree_fegsymbol([[1,5,6],[1,2]])
88
```
"""
function degree_fegsymbol(s)
  s=fullsymbol(s)
  e=length(s)
  d=defectsymbol(s)
  if !(d in [0,1]) return -1 end
  r=ranksymbol(s)
  if d==1 res=div(e*r*(r+1),2)
  else    res=div(e*r*(r-1),2)+r
  end
  res+=e*sum(S->transpose(S)*(0:length(S)-1)-sum(l->div(l*(l+1),2),S),
               filter(!isempty,s))
  gamma=i->sum(mod.(i.+(0:e-1),e).*map(sum,s))
  if d==1 res+=gamma(0)
  else res+=maximum(gamma.(0:e-1))
  end
  res-sum(map(x->div(x*(x-1),2),e*(1:div(sum(length,s),e)-1).+mod(d,2)))
end

"""
`valuation_fegsymbol(s)`

the  valuation of  the fake  degree of  the character  parameterized by the
`e`-symbol `s`.

```julia-repl
julia> valuation_fegsymbol([[1,5,6],[1,2]])
16
```
"""
function valuation_fegsymbol(s)
  s=fullsymbol(s)
  d=defectsymbol(s)
  if !(d in [0,1]) return -1 end
  e=length(s)
  res=e*sum(S->transpose(S)*(length(S)-1:-1:0),filter(!isempty,s))
  gamma=i->sum(mod.(i.+(0:e-1),e).*map(sum,s))
  if d==1 res+=gamma(0)
  else res+=minimum(gamma.(0:e-1))
  end
  res-sum(map(x->div(x*(x-1),2),e*(1:div(sum(length,s),e)-1).+mod(d,2)))
end

" stringsymbol(S) string for symbol S"
function stringsymbol(io,S)
  if isempty(S) return "()" end
  if S[end] isa AbstractVector "("*join(joindigits.(S),",")*")"
  else
    v=E(S[end-1],S[end])
    n= v==1 ? "+" : v==-1 ? "-" : repr(v;context=io)
    "("*join(joindigits.(S[1:end-2]),",")*"$n)"
  end
end

stringsymbol(S)=stringsymbol(stdout,S)

"""
  `lesssymbols(x,y)`  `<` for symbols

  A symbol is smaller than another if the shape is lexicographically smaller,
  or the shape is the same and the the symbol is lexicographically bigger
"""
function lesssymbols(x,y)
  c=cmp(length.(x),length.(y))
  c==-1 || (c==0 && x>y)
end

"""
`defshape(s::Vector{Int})`  

Malle-defect of symbols of shape s. This is an invariant by shift but not 
under cyclic permutations.
"""
function defshape(s)
  e=length(s)
  mod(binomial(e,2)*div(sum(s),e)-sum(i->(i-1)*s[i],1:e),e)
end

"All shapes for e-symbols of rank r, content c, Malle-defect def"
function shapesSymbols(e,r,c=1,def=0)
  if e==1 return [[0]] end
  function f(lim2,len,nb,max) local res,a # possible decreasing shapes
    if nb==1
      if len==0   return [[len]]
      else return Vector{Int}[] end
    end
    res=Vector{Int}[]
    a=div(len,nb-1)
    while a<=max  &&  binomial(a,2)<=lim2  &&  a<=len
      append!(res,map(x->pushfirst!(x,a),f(lim2-binomial(a,2),len-a,nb-1,a)))
      a+=1
    end
    return res
  end

  res=Vector{Int}[]
  m=0
  while true
    new=f(r+div((m*e+c-1)*(m*e+c-e+1),2*e),c+m*e,e,c+m*e)
    if length(new)==0 break end
    append!(res,new)
    m+=1
  end
  res=reduce(vcat,map(x->arrangements(x,e),res))
  # for symbols of content 1 only one circshift of the shape has defshape=0
  filter(s->defshape(s)==def  &&
    all(x->defshape(x)!=def  ||  x<=s,map(i->circshift(s,i),1:length(s)-1)),res)
end

function symbolsshape(r,s)
  S=map(x->symbol_partition_tuple(x,s),
       partition_tuples(r-ranksymbol(map(x->0:x-1,s)),length(s)))
  if !iszero(sum(s)%length(s)) return S end
  S=filter(S)do s
    all(1:length(s)-1)do i
      x=circshift(s,i)
      x==s || lesssymbols(x,s)
    end
  end
  if !iszero(defshape(s)) return S end
  res=[]
  for s in S
    p=findfirst(i->s==circshift(s,1-i),2:length(s))
    if p===nothing push!(res,s)
    else append!(res,map(i->vcat(map(copy,s[1:p]),[div(length(s),p),i]),
                      0:div(length(s),p)-1))
    end
  end
  res
end

"""
`symbols(e,r,c=1,def=0)` `e`-symbols of rank `r`, content `c` and Malle-defect
`def`

An `e`-symbol is a symbol of length `e`.
The content of an `e`-symbol `S` is `sum(length,S)%e`.
The symbols for unipotent  characters of:
  - `G(d,1,r)` are `symbols(d,r)`
  - `G(e,e,r)` are `symbols(e,r,0)`.
  - `G(e,e,r).s₁ᵗ` where `s₁` is the first generator of `G(e,1,r)` and `t|e`
    are `symbols(e,r,0,t)`

```julia-repl
julia> stringsymbol.(symbols(3,2,1))
14-element Vector{String}:
 "(12,0,0)"
 "(02,1,0)"
 "(02,0,1)"
 "(012,12,01)"
 "(01,1,1)"
 "(012,01,12)"
 "(2,,)"
 "(01,2,0)"
 "(01,0,2)"
 "(1,012,012)"
 "(,02,01)"
 "(,01,02)"
 "(0,,012)"
 "(0,012,)"

julia> stringsymbol.(symbols(3,3,0))
12-element Vector{String}:
 "(1+)"
 "(1E(3))"
 "(1E(3,2))"
 "(01,12,02)"
 "(01,02,12)"
 "(012,012,123)"
 "(0,1,2)"
 "(0,2,1)"
 "(01,01,13)"
 "(0,0,3)"
 "(012,,)"
 "(012,012,)"
```
"""
symbols(e,r,c=1,def=0)=vcat(map(s->symbolsshape(r,s),shapesSymbols(e,r,c,def))...)

# See Mal95, 2.11 and 5.7
"""
`fegsymbol(S,p=0)`

Let  `s=[S₁,…,Sₑ]` be an `e`-symbol  given as a `Vector{Vector{Int}}`. This
function  returns as a `CycPol` the fake  degree of the character of symbol
`S`.

```julia-repl
julia> fegsymbol([[1,5,6],[1,2]])
q¹⁶Φ₅Φ₇Φ₈Φ₉Φ₁₀Φ₁₁Φ₁₄Φ₁₆Φ₁₈Φ₂₀Φ₂₂
```

When  given a  second argument  `p` dividing  `e`, and  a first argument of
shape  `(0,…,0)` representing the restriction  to `G(e,e,r)`, works for the
coset `G(e,e,r).s₁ᵖ`.
"""
function fegsymbol(s,p=0)
  if isempty(s) return one(CycPol) end
  s=fullsymbol(s)
  e=length(s)
  r=ranksymbol(s)
  ep=E(e,p)
  function delta(S)
    l=length(S)
    if l<2 return one(CycPol) end
    prod(CycPol(Pol([1],e*S[j])-Pol([1],e*S[i])) for i in 1:l for j in i+1:l)
  end
  function theta(S)
    if iszero(sum(S)) return one(CycPol) end
    prod(CycPol(Pol([1],e*h)-Pol(1)) for l in S for h in 1:l)
  end
  if !allequal(length.(s[2:end])) return zero(CycPol) end
  d=defectsymbol(s)
  if d==1 res=theta([r])
  elseif d==0 res=theta([r-1])*CycPol(Pol([1],r)-Pol(ep))
  else return zero(CycPol)
  end
  res*=prod(S->delta(S)//theta(S),s)
  res//=CycPol(1,sum([div(x*(x-1),2) for x in e.*(1:div(sum(length,s),e)-1).+d%2]))
  if d==1 res*=CycPol(1,sum((0:e-1).*map(sum,s)))
  else
    rot=circshift.(Ref(s),e:-1:1)
    u=map(j->ep^j,0:e-1).*map(s->Pol([1],sum((0:e-1).*map(sum,s))),rot)
    res*=CycPol(sum(u))
    res=div(res,count(==(s), rot))
    if e==2 && ep==-1 res=-res end
  end
  if r==2 && (e>2 && ep==E(e))
    res=subs(res,Pol([E(2e)],1))//E(2e,degree(res))
  end
  return res
end

"""
`gendeg_symbol(s)`

Let  `s=[S₁,…,Sₑ]` be an `e`-symbol  given as a `Vector{Vector{Int}}`. This
function  returns  as  a  `CycPol`  the  generic  degree  of  the unipotent
character parameterized by `s`.

```julia-repl
julia> gendeg_symbol([[1,2],[1,5,6]])
q¹³Φ₅Φ₆Φ₇Φ₈²Φ₉Φ₁₀Φ₁₁Φ₁₄Φ₁₆Φ₁₈Φ₂₀Φ₂₂/2
```
Works for symbols for:

    G(e,1,r) (c==1, d==0)
    G(e,e,r) (c==0, d==0)
   ²G(e,e,r) (c==0, d==1) (e,r even. This includes ²Dₙ, ²B₂, ²G₂)

here  `c` is the content  of the symbol and  `d` the Malle-defect, see [3.9
and 6.4 Malle1995](biblio.htm#Mal95).
"""
function gendeg_symbol(S)
  S=fullsymbol(S)
  r=ranksymbol(S)
  e=length(S)
  sh=length.(S)
  if e==0 return one(CycPol) end
  m=div(sum(sh),e)
  d=sum(sh)% e
  defect=(binomial(e,2)*m-transpose(sh)*(0:e-1))%e
  function theta(S)
    if iszero(sum(S)) return one(CycPol) end
    prod(CycPol(Pol([1],e*h)-Pol(1)) for l in S for h in 1:l)
  end

  # initialize with the q'-part of the group order
  if d==1 res=theta([r])
  elseif d==0 res=theta([r-1])*CycPol(Pol([1],r)-E(e,defect))
  end

  res*=(-1)^(sum((0:e-1).*binomial.(sh,2)))*
    prod(i->prod(j->reduce(*,map(l->reduce(*,map(m->CycPol(l-m),
    filter(m->i<j || degree(m)<degree(l),
           map(l->Pol([E(e,j)],l),S[j+1])));init=one(CycPol)),
           map(l->Pol([E(e,i)],l),S[i+1]));init=one(CycPol)),i:e-1),0:e-1)//
     (prod(theta,S)*(E(4)^binomial(e-1,2)*root(e)^e)^m
      *CycPol(1,sum(map(x->binomial(x,2),e.*(1:m-1).+d))))

  # Dudas' sign change
  if r==1 m=findfirst(isempty,S)
    if !isnothing(m) res*=(-1)^m end
  elseif r==2 && e==3 && S in
    [[[1],[0,1,2],[0,1,2]],[Int[],[0,2],[0,1]],[Int[],[0,1],[0,2]]] res=-res
  elseif r==3 && e==3 && d==0 && S in
    [[[0,1,2],Int[],Int[]],[[0,1,2],[0,1,2],Int[]]] res=-res
  elseif r==4 && e==3 && d==0 && S in
    [[[0,1,3],Int[],Int[]],[[0,1,2,3],[1],[0]],[[0,1,2,3],[0],[1]],
    [[0,1,3],[0,1,2],Int[]],[[0,1,2],[0,1,3],Int[]],[[0,1,2,3],[0,1,2,3],[1]]]
    res=-res
  elseif r==5 && e==3 && d==0 && S in
  [[[0,2,3],Int[],Int[]],[[0,1,2,4],[1],[0]],[[0,1,2,4],[0],[1]],
   [[0,1,2,3,4],[1,2],[0,1]],[[0,1,2,3],[1],[1]],[[0,1,2,3,4],[0,1],[1,2]],
   [[0,1,4],Int[],Int[]],[[0,1,2,3],[2],[0]],[[0,1,2,3],[0],[2]],
   [[0,2,3],[0,1,2],Int[]],[[0,1,3],[0,1,3],Int[]],[[0,1,2,4],[0,1,2,3],[1]],
   [[0,1,2],[0,2,3],Int[]],[[0,1,2,3],[0,1,2,4],[1]],[[0,1,2,3,4],[0,1,2,3,4],[1,2]],
   [[0,1,4],[0,1,2],Int[]],[[0,1,2],[0,1,4],Int[]],[[0,1,2,3],[0,1,2,3],[2]]]
    res=-res
  end

  if d==0 res*=findfirst(i->circshift(S,-i)==S,1:e) end
  if defect==0 || r!=2 || e<=2 return res
  else return E(e)^-1*subs(res,Pol([E(2*e)],1)) # 2I(e)
  end
end

# Ennola of e-symbol s (of content 1 or 0)
# The order of Ennola (order of center of group) is computed automatically:
# it is e for content 1 and gcd(e,rank(s)) for content 0.
function EnnolaSymbol(s)
  repeat=!(s[end] isa Vector)
  if repeat
    times=s[end-1]
    ind=s[end]
    s=fullsymbol(s)
  end
  e=length(s)
  z=sum(length,s)%e==1 ? e : gcd(e,ranksymbol(s))
  if isone(z) return s end
  res=map(x->Int[],s)
  for i in eachindex(s), k in s[i]
    push!(res[1+(i+k*div(e,z))%e],k)
  end
  res=map(sort,res)
  if repeat
    if div(e,times)!=findfirst(i->circshift(res,-i)==res,1:e)
      error("period changed!")
    end
    if e%2==0 ind=(ind+div(e,2)-1)%times end
    # works for types D,I and G(3,3,3) but ???
    res=[res[1:div(e,times)]...,times,ind]
  end
  res
end

function xsp(rho,s,n,d)
  nrsd=rho*div(d^2,4)-s*div(d-mod(d,2),2)
  if n<nrsd return Vector{Vector{Int}}[] end
  return map(partition_tuples(n-nrsd,2)) do S
    S = symbol_partition_tuple(S, d)
    S = map(x->isempty(x) ? x : x+(0:length(x)-1)*(rho-1),S)
    [S[1],S[2].+s]
  end
end

"""
`XSP(ρ,s,n,even=false)`

returns  the union of the Lusztig-Spaltenstein ``X̃^{ρ-s,s}_{n,d}`` for all
`d` even when `even=true`, all `d` odd otherwise; these symbols parametrize
local  systems  on  unipotent  conjugacy  classes  for classical groups. In
[Lusztig2004](biblio.htm#Lus04),  13.2 the notation  is ``{}^ρ X^s_{n,d}``.
The result is a list of lists, each one corresponding to a similarity class
(which  correspond to a given conjugacy  class for the support). If `s==0`,
only positive defects are considered.

  - `XSP(2,1,n)` gives L-S symbols for Sp₂ₙ
  - `XSP(4,2,n)` gives L-S symbols for Sp₂ₙ in char.2
  - `XSP(2,0,n)` gives L-S symbols for SO₂ₙ₊₁ [defect odd]
  - `XSP(2,0,n,true)` gives L-S symbols for SO₂ₙ [defect even]
  - `XSP(4,0,n,true)` gives L-S symbols for SO₂ₙ in char 2

each item is a `NamedTuple` giving some information on the local system.
It has fields

  - `symbol` the Lusztig-Spaltenstein symbol
  - `dimBu` for the support `u` of the local system
  - `Au` describes the character  of `A(u)` for the  local system as a list:
    `true`->sgn, `false`->Id
  - `sp`  parameter (double partition) of the generalized Springer
     correspondent (a character of the relative Weyl group)
"""
function XSP(rho,s,n,even=false)
  d=Int(!Bool(even))
  res = Vector{Vector{Int}}[]
  while true
    S=xsp(rho, s, n, d)
    if iszero(d) S=unique!(sort!(sort.(S))) end
    append!(res,S)
    d+=2
    if isempty(S) break end
  end
  if !iszero(s) && !iszero(mod(d,2))
    d=-1
    while true
      S=xsp(rho, s, n, d)
      append!(res, S)
      d-=2
      if isempty(S) break end
    end
  end
  map(collectby(x->[sort(unique!(union(x...))),sort(unique!(intersect(x...)))],res)) do f
    ii=Vector{Int}[]
    d=sort(symdiff(f[1]...))
    if length(d)>0
      i=[d[1]]
      for j in d[2:end]
        if j-i[end]<rho push!(i, j)
        else push!(ii, i)
          i=[j]
        end
      end
      push!(ii,i)
      ii=filter(x->x[1]>=s,ii)
    end
    u=findfirst(f) do x
      if !(defectsymbol(x) in [0,1]) return false end
      l=zeros(Int,1+sum(length,x))
      l[1:2:2*length(x[1])-1]=x[1]
      l[2:2:2*length(x[2])]=x[2]
      issorted(l[1:sum(length,x)])
    end
    if isnothing(u) error("f=$f") end
    dist=f[u]
    n=sum(length,dist)
    d=defectsymbol(dist)
    m=div(n,2)
    i=sort(reduce(vcat,dist),rev=true)
    n=sum(i.*(0:n-1))-div(rho*m*(m-1)*((4m-5)+6d),6)-s*m*(m+d-1)
    return map(f) do S
      function rr(x,s)
        isempty(x) ? x : reverse(filter(!iszero,x.-(0:length(x)-1).*rho.-s))
      end
      sp=[rr(S[1],0), rr(S[2],s)]
      if defectsymbol(S) == 0
        if sp>reverse(sp) reverse!(sp) end
        if sp[1]==sp[2] sp=[sp[1],2,0] end
      elseif defectsymbol(S)<0 reverse!(sp)
      end
      Au=map(i->intersect(S[1],i)!=intersect(dist[1],i), ii)
      if iszero(s) && S[1]!=S[2]
        if Au[end] Au=.!(Au) end
        Au=Au[1:end-1]
      end
      (symbol=S,Au,dimBu=n,sp)
  end
end
end

showxsp(r)=println("(symbol=",PartitionTupleToString(r.symbol),
    ", sp=",PartitionTupleToString(r.sp),", dimBu=",r.dimBu,", Au=",r.Au,")")

end
