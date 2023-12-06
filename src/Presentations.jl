"""
This  is a  port of  some GAP3/VkCurve  functionality on *presentations* of
*finitely presented groups*.

We  have defined just enough functionality  on finitely presented groups so
that  presentations  can  be  translated  to  finitely presented groups and
vice-versa. The focus is on presentations, the goal being to simplify them.

The  elements of finitely presented groups are `AbsWord` or abstract words,
representing elements of a free group. In order to speed up the algorithms,
the   relators  in  a  presentation   are  not  represented  internally  by
`AbsWord`s, but by lists of positive or negative generator numbers which we
call *Tietze words*.

```julia-repl
julia> @AbsWord a,b # same as a=AbsWord(:a);b=AbsWord(:b)

julia> F=FpGroup([a,b])
FreeGroup(a,b)

julia> G=F/[a^2,b^7,comm(a,a^b),comm(a,a^(b^2))*inv(b^a)]
FreeGroup(a,b)/[a²,b⁷,a⁻¹b⁻¹a⁻¹bab⁻¹ab,a⁻¹b⁻²a⁻¹b²ab⁻²ab²a⁻¹b⁻¹a]

julia> P=Presentation(G) # by default give a summary
Presentation: 2 generators, 4 relators, total length 30

julia> relators(P)
4-element Vector{AbsWord}:
 a²
 b⁷
 ab⁻¹abab⁻¹ab
 b⁻²ab²ab⁻²ab²ab⁻¹
```

```julia-rep1
julia> showgens(P)
1. a 10 occurrences involution
2. b 20 occurrences

julia> dump(P) # here in relators inverses are represented by capitalizing
# F relator
1:3 aa
2:0 bbbbbbb
3:0 aBabaBab
4:0 abbaBBabbaBBB
gens=AbsWord[a, b] involutions:AbsWord[a] modified=false numredunds=0

julia> display_balanced(P)
1: a=A
2: bbbbbbb=1
3: aBab=BAbA
4: BBabbaBBabbaB=1
```

The  functions `Presentation`  and `FpGroup`  create a  presentation from a
finitely presented group and vice versa.

for  more  information  look  at  the  help  strings  of `AbsWord, FpGroup,
Presentation,     relators,    display_balanced,    simplify,    conjugate,
tryconjugate`. 

A  minimal thing to add to this package so it would be a reasonable package
for finitely preented groups is the Coxeter-Todd algorithm.
"""
module Presentations
## Changing Presentations
#
#The  functions `AddGenerator`, `AddRelator`, `RemoveRelator` can be used to
#change a presentation. In general, they will change the isomorphism type of
#the  group defined  by the  presentation, hence,  though they are sometimes
#used  as subroutines by Tietze transformations functions like `Substitute`,
#they do *not* perform Tietze transformations themselves.
#
## Tietze Transformations
#
#The  functions described in this section can be used to modify a
#group presentation by Tietze transformations.
#
#In  general, the aim of such modifications  will be to *simplify* the given
#presentation,  i.e., to reduce  the number of  generators and the number of
#relators  without increasing too much the  sum of all relator lengths which
#we  will  call  the  *total  length*  of the presentation. Depending on the
#concrete presentation under investigation one may end up with a nice, short
#presentation   or  with  a   very  huge  one.
#
#There  is no way  to find the  shortest presentation which  can be obtained
#from  a given  one. Therefore,  what we  offer are  some lower-level Tietze
#transformation   functions  and,  in  addition,  a  heuristic  higher-level
#function   (which  of  course   cannot  be  the   optimal  choice  for  all
#presentations).
#
#The design of these functions follows closely the concept of the ANU Tietze
#transformation  program designed by George Havas cite{Hav69} which has been
#available  from Canberra since 1977 in a stand-alone version implemented by
#Peter   Kenne  and  James  Richardson  and   later  on  revised  by  Edmund
#F.~Robertson (see cite{HKRR84}, cite{Rob88}).
#
#The  higher-level  function  is  `simplify`.  The lower-level functions are
#`Eliminate`, `Search`, `SearchEqual`, and `FindCyclicJoins`.
#
#Some  of  these  functions  may  eliminate  generators,  but  they  do *not*
#introduce  new generators. However,  sometimes you will  need to substitute
#certain  words as  new generators  in order  to improve  your presentation.
#Therefore  there are the  functions `Substitute` and `SubstituteCyclicJoins`
#which introduce new generators.
#
#Finally  the functions `tracing`  and `images` can
#be used to determine and to display the images or preimages of the involved
#generators under the isomorphism which is defined by the sequence of Tietze
#transformations which are applied to a presentation.
#
#The  functions, `show_pairs`, and `PrintOptions`,  can be useful. There are
#also  the  *Tietze  options*:  parameters  which  essentially influence the
#performance  of the  functions mentioned  above; they  are not specified as
#arguments of function calls. Instead, they are stored in the presentation.
using ..Util: xprintln, fromTeX, rio
function stringind(io::IO,n::Integer)
  if get(io,:TeX,false) 
    n in 0:9 ? "_"*string(n) : "_{"*string(n)*"}"
  elseif get(io,:limit,false)
    if n<0 res=['₋']; n=-n else res=Char[] end
    for i in reverse(digits(n)) push!(res,Char(0x2080+i)) end
    String(res)
  else "_"*string(n)
  end
end

using PermGroups
using Combinat: tally
export AbsWord, @AbsWord, Presentation, FpGroup, Go, GoGo, conjugate, 
       tryconjugate, simplify, relators, display_balanced, tracing, images,
       showgens

plural(n,w)=string(n)*" "*w*(n==1 ? "" : "s")

#------------------ Abstract Words ----------------------------------
"""
An  `AbsWord` represents an  element of the  free group on some generators.
The  generators  are  indexed  by  `Symbols`.  For  example  the  `Absword`
representing `a³b⁻²a` is represented internally as
`[:a => 3, :b => -2, :a => 1]`. The mulitiplcation follows the group rule:
```julia-repl
julia> w=AbsWord([:a => 3, :b => -2, :a => 1])
a³b⁻²a

julia> w*AbsWord([:a=>-1,:b=>1])
a³b⁻¹

```
A positive `AbsWord` may be obtained by giving `Symbols` as arguments
```julia-repl
julia> AbsWord(:b,:a,:a,:b)
ba²b
```
"""
struct AbsWord 
  d::Vector{Pair{Symbol,Int}}
  function AbsWord(v::Vector{Pair{Symbol,Int}};check=true)
    if check && length(v)>0
      ri=1
      for i in 2:length(v)
        if ri==0 || first(v[i])!=first(v[ri])
          ri+=1
          v[ri]=v[i]
        else c=last(v[ri])+last(v[i])
          if iszero(c) && ri>0 ri-=1 
          else v[ri]=first(v[ri])=>c
          end
        end
      end
      resize!(v,ri)
    end
    new(v)
  end
end

AbsWord(x::Symbol...)=AbsWord([s=>1 for s in x])

"`@AbsWord x,y` is the same as `x=AbsWord(:x);y=AbsWord(y)`"
macro AbsWord(t) 
  if t isa Expr
    for v in t.args
      Base.eval(Main,:($v=AbsWord($(Core.QuoteNode(Symbol(v))))))
    end
  elseif t isa Symbol
    Base.eval(Main,:($t=AbsWord($(Core.QuoteNode(t)))))
  end
end

"AbsWord(s::String) defines an `AbsWord` from a line of display_balanced"
function AbsWord(s::AbstractString)
  ss=split(s,"=")
  if ss[end]=="1" pop!(ss) end
  if length(ss)==2 return AbsWord(ss[1])*inv(AbsWord(ss[2])) end
  s=ss[1]
  res=Pair{Symbol,Int}[]
  for c in collect(s)
    if islowercase(c) push!(res,Symbol(c)=>1)
    else push!(res,Symbol(lowercase(c))=>-1)
    end
  end
  AbsWord(res)
end

function Base.show(io::IO,a::AbsWord)
  if isone(a) print(io,".") end
  for (s,c) in a.d
    print(io,string(s))
    if c!=1 print(io,fromTeX(io,"^{$c}")) end
  end
end

Base.one(::Type{AbsWord})=AbsWord(Pair{Symbol,Int}[])
Base.isone(a::AbsWord)=isempty(a.d)
Base.:*(a::AbsWord,b::AbsWord)=AbsWord(vcat(a.d,b.d))
Base.inv(a::AbsWord)=AbsWord([k=>-v for (k,v) in reverse(a.d)];check=false)
Base.:^(a::AbsWord, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                       Base.power_by_squaring(inv(a),-n)
Base.:^(a::AbsWord, b::AbsWord)=inv(b)*a*b
Base.:/(a::AbsWord, b::AbsWord)=a*inv(b)
Base.:\(a::AbsWord, b::AbsWord)=inv(a)*b
Base.length(a::AbsWord)=sum(x->abs(last(x)),a.d;init=0)
Base.copy(a::AbsWord)=AbsWord(copy(a.d))
Base.:(==)(a::AbsWord,b::AbsWord)=a.d==b.d

"returns symbol for an AbsWord of length 1"
function mon(w::AbsWord)
  if length(w.d)!=1 || last(w.d[1])!=1 error("not a generator") end
  first(w.d[1])
end

function Base.getindex(w::AbsWord,i::Integer)
  for (s,n) in w.d
    if i>abs(n) i-=abs(n)
    elseif n<0 return s=>-1
    else return s=>1
    end
  end
  error("index out of bounds")
end

Base.getindex(w::AbsWord,i::AbstractVector{<:Integer})=AbsWord(getindex.(Ref(w),i))
#-------------------------- FpGroups -----------------------------
struct FpGroup <: Group{AbsWord}
  gens::Vector{AbsWord}
  rels::Vector{AbsWord}
end

FpGroup(gens::Vector{AbsWord})=FpGroup(gens,AbsWord[])

FpGroup(s::Symbol...)=FpGroup(AbsWord.(collect(s)),AbsWord[])

function Base.show(io::IO,G::FpGroup)
  print(io,"FreeGroup(",join(gens(G),","),")")
  if length(G.rels)>0 print(io,"/[");join(io,G.rels,",");print(io,"]") end
end

function Base.:/(G::FpGroup,rel::Vector{AbsWord})
  append!(G.rels,rel)
  G
end

# Tietze word from AbsWord
function Groups.word(G::FpGroup,w::AbsWord)
  res=Vector{Int}(undef,length(w))
  x=1
  for (s,m) in w.d
    p=findfirst(x->x.d[1][1]==s,gens(G))
    if isnothing(p) error(w," is not a word for ",G) end
    for i in 1:abs(m) 
      res[x]=m<0 ? -p : p
      x+=1
    end
  end
  res
end

function Groups.elements(F::FpGroup,i)
  l1=map(x->x.d[1],gens(F))
  l=map(x->x[1]=>-x[2],l1)
  l=vcat(l,l1)
  res=AbsWord.(collect.(Iterators.product(fill(l,i)...)))
  filter(x->length(x)==i,res)
end

#--------------------------  Tietze words and structs -------------
"""
`TietzeWord(word::AbsWord, generators::Vector{AbsWord})`

Let  `generators` be a  list of abstract  generators and `word` an abstract
word   in  these   generators.  The   function  `TietzeWord`   returns  the
corresponding (reduced) Tietze word.

```julia-repl
julia> F=FpGroup(:a,:b,:c)
FreeGroup(a,b,c)

julia> Presentations.TietzeWord(comm(F(1),F(2))*inv(F(3)^2*F(2)),gens(F))
5-element Vector{Int64}:
 -1
 -2
  1
 -3
 -3
```
"""
function TietzeWord(w::AbsWord,gens::Vector{AbsWord})
  ss=mon.(gens)
  res=Int[]
  for (s,c) in w.d
    p=findfirst(==(s),ss)
    if p===nothing error(s," is not in ",gens) end
    if c>0 append!(res,fill(p,c))
    else append!(res,fill(-p,-c))
    end
  end
  res
end

@GapObj mutable struct Presentation
  generators::Vector{AbsWord} # copy of initial gens
  inverses::Vector{Int}
  relators::Vector{Vector{Int}}
  flags::Vector{Int}
  modified::Bool
  numredunds::Int
end
# possible fields: tietze, nextFree, eliminationsLimit,
# expandLimit, generatorsLimit, lengthLimit, loopLimit, debug,
# saveLimit, searchSimultaneous, protected, status

" `relators(P::Presentation)` relators of `P` as `AbsWord`s. "
relators(P::Presentation)=AbsWord.(P.relators,Ref(P.generators))

"""
`FpGroup(P::Presentation)`

returns the finitely presented group defined  by the presentation `P`.
"""
function FpGroup(P::Presentation)
  debug(P,3,"# converting presentation to FpGroup")
  if P.numredunds>0 RemoveGenerators(P) end
  sort!(P)
  G=FpGroup(deepcopy(P.generators),relators(P))
  if haskey(P, :imagesOldGens)
    G.imagesOldGens=deepcopy(P.imagesOldGens)
  end
  if haskey(P, :preImagesNewGens)
    G.preImagesNewGens=deepcopy(P.preImagesNewGens)
  end
  G
end

Presentation(s::Vector{Symbol},r,pl=1)=Presentation(AbsWord.(s),r,pl)

function Presentation(gens::Vector{AbsWord},grels::Vector{AbsWord},debug::Int=1)
  numgens=length(gens)
  rels=map(w->reduceword(TietzeWord(w,gens);cyclically=true),grels)
  P=Presentation(deepcopy(gens),numgens:-1:-numgens,rels,fill(1,length(rels)),
                 false, 0, Dict{Symbol,Any}())
  P.nextFree=numgens+1
  P.eliminationsLimit=100
  P.expandLimit=150
  P.generatorsLimit=0
  P.lengthLimit=0 # infinity
  P.loopLimit=0 # infinity
  P.debug=debug
  P.saveLimit=10
  P.searchSimultaneous=20
  PrintStatus(P,2)
  P.protected=length(P.generators);HandleLength1Or2Relators(P);P.protected=0
  sort!(P)
  PrintStatus(P,2)
  P
end

function Presentation(v::Vector{Vector{Int}})
  gens=sort(unique(abs(x) for r in v for x in r))
  if length(gens)!=maximum(gens) error("hummm!") end
  vars=length(gens)<=26 ? Symbol.('a':'a'+length(gens)-1) :
      map(i->Symbol("x",stringind(rio(),i)),1:length(gens))
  Presentation(AbsWord.(vars),map(r->AbsWord(r,vars),v))
end
  
"""
`Presentation( G::FpGroup[, debug=1])`

returns  the  presentation  corresponding  to  the given finitely presented
group `G`.

The  optional `debug` parameter  can be used  to restrict or  to extend the
amount  of output  provided by  Tietze transformation  functions when being
applied  to the created  presentation. The default  value 1 is designed for
interactive  use and implies  explicit messages to  be displayed by most of
these functions. A `debug` value of 0 will suppress these messages, whereas
a `debug` value of 2 will enforce some additional output.
"""
Presentation(G::FpGroup,printlevel::Integer=1)=Presentation(G.gens,G.rels)

# takes as input the output of display_balanced
function Presentation(s::String;level=1)
  s=replace(s,r"\n\s*="=>"=")
  l=split(s,"\n")
  l=map(s->replace(s,r"^ *[0-9]*: *"=>""),l)
  l=filter(x->match(r"^ *$",x)===nothing,l)
  rels=AbsWord.(l)
  l=map(x->first.(x.d),rels)
  atoms=length(l)==1 ? unique(l[1]) : union(map(x->first.(x.d),rels)...)
  sort!(atoms)
  Presentation(AbsWord.(atoms),rels,level)
end

Base.getindex(T::Presentation,i)=T.inverses[length(T.generators)+1-i]
Base.setindex!(T::Presentation,j,i)=T.inverses[length(T.generators)+1-i]=j

# reduce cyclically a Tietze word
function reduceword(w::Vector{Int};cyclically=false)
  i=1;
  res=Int[]
  for j in eachindex(w)
    if isempty(res) || res[end]!=-w[j] if w[j]!=0 push!(res,w[j]) end
    else pop!(res)
    end
  end
  b=1;e=length(res)
  if cyclically while b<e && res[b]==-res[e] b+=1;e-=1 end end
  res[b:e]
end
   
# reduce cyclically a Tietze word using the table of inverses
function reduceword(w::Vector{Int},T::Presentation)
  i=1;
  res=Int[]
  for j in eachindex(w)
    if T[w[j]]==0 continue end
    if isempty(res) || T[res[end]]!=T[-w[j]] push!(res,T[w[j]])
    else pop!(res)
    end
  end
  b=1;e=length(res)
  while b<e && T[res[b]]==T[-res[e]] b+=1;e-=1 end
  res[b:e]
end

# sorts by length relators and remove duplicates and empty ones
function Base.sort!(T::Presentation)
  if isempty(T.relators) return end
  p=sortperm(T.relators,by=length)
  p=unique(i->T.relators[i],p)
  if isempty(T.relators[p[1]]) p=p[2:end] end
  if p==eachindex(T.relators) return end
  if T.debug>=3 println("#sort! the relators") end
  numrels=length(p)
  T.relators[1:numrels].=T.relators[p];resize!(T.relators,numrels)
  T.flags[1:numrels].=T.flags[p];resize!(T.flags,numrels)
end

function debug(P,i,arg...)
  if P.debug>=i println("#",arg...) end
end

"""
`PrintStatus(P::Presentation)`

updates  and prints  the current  status of  a presentation  `P`, i.e., the
number  of generators, the number of relators,  and the total length of all
relators.

The  printing is suppressed if  none of the three  values has changed since
the last call.
"""
function PrintStatus(P::Presentation,level=0)
  if P.debug<level return end
  status=[length(P.generators)-P.numredunds, length(P.relators), 
          sum(length,P.relators;init=0)]
  if !haskey(P,:status) || status!=P.status
    P.status=status
    show(stdout,P)
    println()
  end
end

function Base.show(io::IO, p::Presentation)
  print(io,"Presentation")
  if haskey(p,:name) print(io," ",p.name) end
  print(io,": ",plural(length(p.generators)-p.numredunds,"generator"),
        ", ",plural(length(p.relators),"relator"),", total length ",
        sum(length,p.relators;init=0))
end

function Base.dump(T::Presentation)
  l=filter(i->T[-i]!=-T[i],eachindex(T.generators))
  l1=filter(i->T[-i]==T[i],l)
  println("# F relator")
  for i in eachindex(T.relators)
    println("$i:",T.flags[i]," ",alphab(T.relators[i],T.generators))
  end
  print("gens=",T.generators)
  if length(l1)>0 print(" involutions:",T.generators[l1])end
  println(" modified=",T.modified," numredunds=",T.numredunds)
  if l1!=l display(stack(map(i->[T.generators[i],T[i],T[-i]],l))) end
end

function display_balanced(T,dumb=false)
  f(i,w)=vcat(w,w)[i:i+div(length(w),2)-1]
  used=Set{Int}();for x in T.relators, y in x push!(used,abs(y)) end
  if length(T.generators)>length(used)
   print("There are ",length(T.generators)-length(used)," free generators\n")
  end
  for (i,w) in enumerate(T.relators)
    lw=length(w)
    if isodd(lw) || dumb println(i, ": ", alphab(w,T.generators), "=1")
    elseif lw>0
      m=argmax(map(i->count(>(0), f(i,w)), 1:lw))
      println("$i: ",alphab(f(m,w),T.generators),"=",alphab(-reverse(f(m+div(lw,2),w)),T.generators))
    end
  end
end

"""
`showgens(P,list=eachindex(P.generators))`

prints  the generators of `P` with the total number of their occurrences in
the  relators, and notes involutions. A  second `list` argument prints only
those generators.
"""
function showgens(T::Presentation,list=eachindex(T.generators))
  gens=T.generators
  if isempty(gens) println("#I  there are no generators");return end
  occur=Occurrences(T)
  if list isa Integer list=[list] end
  for i in list
    print(i, ". ", gens[i], " ",plural(occur[1][i],"occurrence"))
    if T[-i]>0 print(" involution") end
    println()
  end
end

"""
`show_pairs(P,n=10)`

shows  the  `n`  most  frequently  occurring squarefree relator subwords of
length 2 with their number of occurrences.`n=0` is interpreted as infinity.

This  is  useful  information  in  the context of the `Substitute` command.
"""
function show_pairs(P::Presentation,n::Integer=10)
  for (m,(num,i,j,k)) in enumerate(MostFrequentPairs(P,n))
    geni=P.generators[i];if k>1 geni=inv(geni) end
    genj=P.generators[j];if isodd(k) genj=inv(genj) end
    xprintln("#I $m. ",plural(num,"occurrence")," of ",geni*genj)
  end
end

"""
`AbsWord(word, generators)`

Tranforms Tietze word `word` to an absword, using given generators `gens`.

```julia-repl
julia> AbsWord([-1,-2,1,-3,-3],AbsWord.([:a,:b,:c]))
a⁻¹b⁻¹ac⁻²
```
"""
AbsWord(tz::Vector{Int},gens::Vector{AbsWord})=AbsWord(tz,mon.(gens))
AbsWord(tz::Vector{Int},ss::Vector{Symbol})=AbsWord(map(i->ss[abs(i)]=>sign(i),tz))

(P::Presentation)(l::Int...)=AbsWord(collect(l),P.generators)

"""
`AddGenerator( P[, generator])`

`AddGenerator` adds a new generator to the list of generators.

If  you don't specify a second  argument, then `AddGenerator` will define a
new  abstract generator `_xi` where `i` is the least positive integer which
has not yet been used as a generator number. Though this new generator will
be  printed as `_xi`, you  will have to use  the external variable `P.i` if
you want to access it.

If  you specify  a second  argument, then  `generator` must  be an abstract
generator which does not yet occur in the presentation. `AddGenerator` will
add it to the presentation.
"""
function AddGenerator(P::Presentation,gen=nothing)
  Check(P)
  if gen===nothing
    gen=NewGenerator(P)
    debug(P,1,"AddGenerator  new generator is ", gen)
  else
    if !(gen isa AbsWord && length(gen)==1)
      error("second argument must be an abstract generator")
    end
    if gen in P.generators || inv(gen) in P.generators
      println("#I generator ",gen," is already in the presentation")
      return
    end
    push!(P.generators,gen)
    P.inverses=vcat([length(P.generators)],P.inverses, [-length(P.generators)])
  end
  tracing(P,false)
  P.modified=true
end

"""
NewGenerator(P::Presentation)

defines a new abstract generator and adds it to `P`.

Let i be the smallest positive integer for which the generator `_xi` is not
a  generator of  `P`. A  new abstract  generator `_xi`  is defined and then
added to `P.generators`.

Warning:  `NewGenerator`  is  an  internal  subroutine  of  the  Tietze
routines.  You should not call it.  Instead, you should call the function
`AddGenerator`, if needed.
"""
function NewGenerator(P::Presentation)
  new=P.nextFree
  while true
    gen=AbsWord(Symbol("_x",new))
    if gen in P.generators new+=1;continue end
  end
  P.nextFree=new+1
  push!(P.generators,gen)
  P.inverses=vcat([length(P.generators)],P.inverses,[-length(P.generators)])
  gen
end

"""
`AddRelator(P::Presentation, word::AbsWord)`

adds  the word `word` to the list of  relators.  `word` must
be a word in the generators of the given presentation.
"""
function AddRelator(P::Presentation, word::AbsWord)
  Check(P)
  debug(P,3,"AddRelator adding ", word)
  rel=reduceword(TietzeWord(word,P.generators);cyclically=true)
  if !isempty(rel)
    push!(P.relators,rel)
    push!(P.flags,1)
    P.modified=true
  end
  tracing(P,false)
end

"""
HandleLength1Or2Relators(presentation) . . . handle short relators

`HandleLength1Or2Relators`  searches  for  relators  of  length  1 or 2 and
performs suitable Tietze transformations for each of them:

Generators occurring in relators of length 1 are eliminated.

Generators  occurring  in  square  relators  of  length  2 are marked to be
involutions.

If  a  relator  of  length  2  involves  two different generators, then the
generator  with the larger  number is substituted  by the other  one in all
relators and finally eliminated from the set of generators.
"""
function HandleLength1Or2Relators(T::Presentation)
  debug(T,3,"Handle length 1 or 2 relators")
  gens=T.generators
  rels=T.relators
  done=false
  while !done
    done=true
    i=0
    while i<length(T.relators)
      i+=1
      lg=length(rels[i])
      if 0<lg<=2 && T.flags[i]<=2
        rep1=T[rels[i][1]]
        if lg==1
          rep1=abs(rep1)
          if rep1>T.protected
            if T[rep1]==rep1 T.numredunds+=1 end
            T[rep1]=T[-rep1]=0
            debug(T,2,"Handle12 eliminating ",gens[rep1]," redund=",T.numredunds)
            UpdateGeneratorImages(T, rep1, Int[])
            done=false
          end
        else
          rep2=T[rels[i][2]]
          if abs(rep2)<abs(rep1) rep1,rep2=rep2,rep1 end
          if rep1<0 rep1,rep2=(-rep1,-rep2) end
          if rep1==0 # already eliminated
            rep2=abs(rep2)
            if rep2>T.protected
              if T[rep2]==rep2 T.numredunds+=1 end
              T[rep2]=T[-rep2]=0
              debug(T,2,"Handle12 eliminating ",gens[rep2]," redund=",T.numredunds)
              UpdateGeneratorImages(T, rep2, Int[])
              done=false
            end
          elseif rep1!=-rep2  # otherwise not reduced
            if rep1!=rep2 # not an involution
              if T[rep2]==T[-rep2] && T[-rep1]<0 # rep2^2=1 => rep1^2=1
                push!(rels,[rep1, rep1])
                push!(T.flags,1)
              end
              if abs(rep2)>T.protected
                if T[rep2]==rep2 T.numredunds+=1 end
                T[rep2]=T[-rep1];T[-rep2]=rep1
                if haskey(T,:imagesOldGens) || T.debug>=2
                  if rep2>0 rep1=T[-rep1] end
                  debug(T,2,"Handle12 eliminating ",gens[abs(rep2)],"=",
                      alphab([rep1],T.generators)," redund=",T.numredunds)
                  UpdateGeneratorImages(T,abs(rep2),[rep1])
                end
                done=false
              end
            elseif T[-rep1]<0  # an involution not yet detected
              rels[i]=[rep1,rep1]
              T.flags[i]=3
              T[-rep1]=rep1
              done=false
            end
          end
        end
      end
    end
    if !done
      for i in eachindex(T.generators)
        if T[i]!=i T[i],T[-i]=T[T[i]],T[-T[i]] end
      end
# the next loop is FunTzReplaceGens
      for i in eachindex(T.relators)
       if length(T.relators[i])==2 && T.flags[i]==3 && abs(T[T.relators[i][1]])== abs(T.relators[i][1])continue end#don't remove involutions
        T.relators[i]=reduceword(T.relators[i],T)
        T.flags[i]=1
      end
    end
  end
  if T.numredunds>0 RemoveGenerators(T) end
end

"""
RelsViaCosetTable(G,cosets)  . . . . . . . construct relators for the
RelsViaCosetTable(G,cosets,ggens)  . . . . . . . . . given concrete
RelsViaCosetTable(G,cosets,F,words,H,F1)  . . . . . . . group

`RelsViaCosetTable`  constructs a defining set of relators  for the given
concrete group using John Cannon's relations finding algrotithm.

It is a  subroutine  of function  `PresentationViaCosetTable`.  Hence its
input  and  output  are   specifically  designed  for  this  purpose.  In
particular, it does not check the arguments.

G,            # given group
cosets,       # right cosets of G with respect to H
"""
function RelsViaCosetTable(G,cosets,arg...)
  #  F,        # given free group
  #  words,    # given words for the generators of H
  #  H,        # subgroup, if specified
  #  F1,       # f.p. group isomorphic to H
  #  F2,       # f.p. group isomorphic to G
  #  ng1,      # position number of identity element in G
  #  nh1,      # position number of identity element in H
  #  perms,    # permutations induced by the gens on the cosets
  #  stage,    # 1 or 2
  #  table,    # columns in the table for gens
  #  rels,     # representatives of the relators
  #  relsGen,  # relators sorted by start generator
  #  subgroup, # rows for the subgroup gens
  #  rels1,    # list of relators
  #  app,      # arguments list for `MakeConsequences`
  #  index,    # index of the table
  #  col,      # generator col in auxiliary table
  #  perm,     # permutations induced by a generator on the cosets
  #  gens,     # abstract gens in which the relators are written
  #  gens2,    # the above abstract gens and their inverses
  #  ggens,    # concrete generators of G
  #  ngens,    # number of generators of G
  #  ngens2,   # twice the above number
  #  order,    # order of a generator
  #  actcos,   # part 1 of Schreier vector of G by H
  #  actgen,   # part 2 of Schreier vector of G by H
  #  tab0,     # auxiliary table in parallel to table <table>
  #  cosRange, # range from 1 to index (= number of cosets)
  #  genRange, # range of the odd integers from 1 to 2*ngens-1
  #  geners,   # order in which the table cols are worked off
  #  next,     # local coset number
  #  left1,    # part 1 of Schreier vector of H by trivial group
  #  right1,   # part 2 of Schreier vector of H by trivial group
  #  n,        # number of subgroup element
  #  words2,   # words for the generators of H and their inverses
  #  h         # subgroup element
  if length(arg)==1 ggens=arg[1]
  else ggens=G[:generators]
  end
  ngens=length(ggens)
  ngens2=ngens * 2
  if cosets[1] in G ng1=PositionSorted(cosets, cosets[1] ^ 0)
  else ng1=1
  end
  index=length(cosets)
  tab0=[]
  table=[]
  subgroup=[]
  cosRange=1:index
  genRange=map((i->begin 2i-1 end), 1:ngens)
  if length(arg)<2
      stage=1
      F2=FreeGroup(ngens)
      rels=[]
  else
      stage=2
      F=arg[1]
      words=arg[2]
      F2=Group(F[:generators], IdWord)
      if haskey(F, :namesGenerators)
          F2[:namesGenerators]=F[:namesGenerators]
      end
      H=arg[3]
      nh1=PositionSorted(H[:elements], (H[:elements])[1] ^ 0)
      F1=arg[4]
      rels=map(rel->MappedWord(rel, F1[:generators], words), F1[:relators])
      left1=F1[:actcos]
      right1=F1[:actgen]
      words2=[]
      for i in 1:length(F1[:generators])
          push!(words2, words[i])
          push!(words2, words[i] ^ -1)
      end
  end
  gens=F2[:generators]
  gens2=[]
  perms=map((gen->begin Permutation(gen, cosets, OnRight) end), ggens)
  for i in 1:ngens
      push!(gens2, gens[i])
      push!(gens2, gens[i] ^ -1)
      perm=perms[i]
      col=OnTuples(cosRange, perm)
      gen=0*cosRange
      push!(tab0, col)
      push!(table, gen)
      order=Order(G, ggens[i])
      if order==2 push!(rels, gens[i] ^ 2)
      else col=OnTuples(cosRange, perm ^ -1)
          gen=0*cosRange
      end
      push!(tab0, col)
      push!(table, gen)
  end
  cosets=0*cosRange
  actcos=0*cosRange
  actgen=0*cosRange
  cosets[1]=ng1
  actcos[ng1]=ng1
  j=1
  i=0
  while i<index
      i+=1
      c=cosets[i]
      g=0
      while g<ngens2
          g+=1
          next=(tab0[g])[c]
          if next>0 && actcos[next]==0
              g1=(g+2 * mod(g, 2))-1
              table[g][c]=next
              table[g1][next]=c
              tab0[g][c]=0
              tab0[g1][next]=0
              actcos[next]=c
              actgen[next]=g
              j=j+1
              cosets[j]=next
              if j==index
                  g=ngens2
                  i=index
              end
          end
      end
  end
  rels=RelatorRepresentatives(rels)
  app=fill(0,11)
  app[1]=table
  app[5]=subgroup
  if stage==2
      relsGen=RelsSortedByStartGen(gens, rels, table)
      app[4]=relsGen
      for g in genRange
          gen0=tab0[g]
          for c=cosRange
              if gen0[c]==0
                  app[10]=g
                  app[11]=c
                  n=MakeConsequences(app)
              end
          end
      end
  end
  geners=1:ngens2
  for i in cosets
      for j=geners
          if table[j][i]==0
              g=(j+2 * mod(j, 2))-1
              c=tab0[j][i]
              table[j][i]=c
              table[g][c]=i
              tab0[j][i]=0
              tab0[g][c]=0
              rel=IdWord
              while c != ng1
                  g=actgen[c]
                  rel=rel * gens2[g] ^ -1
                  c=actcos[c]
              end
              rel=rel ^ -1 * gens2[j] ^ -1
              c=i
              while c != ng1
                  g=actgen[c]
                  rel=rel * gens2[g] ^ -1
                  c=actcos[c]
              end
              if stage==2
                  h=MappedWord(rel, gens, ggens)
                  n=PositionSorted(H[:elements], h)
                  while n != nh1
                      g=right1[n]
                      rel=rel * words2[g] ^ -1
                      n=left1[n]
                  end
              end
              rels1=RelatorRepresentatives([rel])
              if length(rels1)>0
                  rel=rels1[1]
                  if !rel in rels
                      push!(rels, rel)
                  end
              end
              relsGen=RelsSortedByStartGen(gens, rels, table)
              app[4]=relsGen
              for g in genRange
                  gen=table[g]
                  gen0=tab0[g]
                  inv0=tab0[g+1]
                  for c=cosRange
                      if gen[c]>0
                          gen0[c]=0
                          inv0[gen[c]]=0
                      end
                  end
              end
              for g in genRange
                  gen0=tab0[g]
                  for c=cosRange
                      if gen0[c]==0
                          app[10]=g
                          app[11]=c
                          n=MakeConsequences(app)
                      end
                  end
              end
          end
      end
  end
  F2=F2 // rels
  F2[:actcos]=actcos
  F2[:actgen]=actgen
  return F2
end

function Ignore end
if !@isdefined InfoFpGroup1 
  InfoFpGroup1=Ignore 
end
if !@isdefined InfoFpGroup2 
  InfoFpGroup2=Ignore 
end

"""
`PresentationViaCosetTable( G[, F, words] )`

`PresentationViaCosetTable`  constructs a presentation  for the given group
`G`.  The method  being used  is John  Cannon's relations finding algorithm
which has been described in cite{Can73} or in cite{Neu82}.

In  its first form,  if only the  group `G` has  been specified, it applies
Cannon's  single stage  algorithm which,  by plain  element multiplication,
computes a coset table of `G` with respect to its trivial subgroup and then
uses coset enumeration methods to find a defining set of relators for `G`.

|    gap> G := GeneralLinearGroup( 2, 7 );
    GL(2,7)
    gap> G.generators;
    [ [ [ Z(7), 0*Z(7) ], [ 0*Z(7), Z(7)^0 ] ],
      [ [ Z(7)^3, Z(7)^0 ], [ Z(7)^3, 0*Z(7) ] ] ]
    gap> Size( G );
    2016
    gap> P := PresentationViaCosetTable( G );
    << presentation with 2 gens and 5 rels of total length 46 >>
    gap> relators( P );
    &I  1. f.2^3
    &I  2. f.1^6
    &I  3. f.1*f.2*f.1*f.2*f.1*f.2*f.1*f.2*f.1*f.2*f.1*f.2
    &I  4. f.1*f.2*f.1^-1*f.2*f.1*f.2^-1*f.1^-1*f.2*f.1*f.2*f.1^-1*f.2^-1
    &I  5. f.1^2*f.2*f.1*f.2*f.1*f.2^-1*f.1^-1*f.2^-1*f.1^3*f.2^-1 |
    gap> DisplayPresentation(P);
    1: bbb=1
    2: aaa=AAA
    3: bababa=ABABAB
    4: baBAba=aBAbaB
    5: aababaBABaaaB=1

The  second form  allows to  call Cannon's  two stage algorithm which first
applies  the single stage  algorithm to an  appropriate subgroup `H` of `G`
and  then uses the resulting relators of `H`  and a coset table of `G` with
respect  to `H` to find relators of  `G`. In this case the second argument,
`F`,  is assumed to be  a free group with  the same number of generators as
`G`, and `words` is expected to be a list of words in the generators of `F`
which, when being evaluated in the corresponding generators of `G`, provide
subgroup generators for `H`.

|    gap> M12 := MathieuGroup( 12 );;
    gap> M12.generators;
    [ ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11), ( 3, 7,11, 8)( 4,10, 5, 6),
      ( 1,12)( 2,11)( 3, 6)( 4, 8)( 5, 9)( 7,10) ]
    gap> F := FreeGroup( "a", "b", "c" );
    Group( a, b, c )
    gap> words := [ F.1, F.2 ];
    [ a, b ]
    gap> P := PresentationViaCosetTable( M12, F, words );
    << presentation with 3 gens and 10 rels of total length 97 >>
    gap> DisplayPresentation(P);
    1: c=C
    2: bb=BB
    3: cac=ACA
    4: aBBaBBaBB=1
    5: aaaaaaaaaaa=1
    6: aaBaab=bAbbaa
    7: bAbABa=AbaBaB
    8: aabaaBBAbABAB=1
    9: aaBABacbcabab=1
    10: aaabaabAAcabAca=1
    gap> G := FpGroup(P);
    Group( a, b, c )
    gap> G.relators;
    [ c^2, b^4, a*c*a*c*a*c, a*b^-2*a*b^-2*a*b^-2, a^11,
      a^2*b*a^-2*b^-2*a*b^-1*a^2*b^-1,
      a*b*a^-1*b*a^-1*b^-1*a*b*a^-1*b*a^-1*b^-1,
      a^2*b*a^2*b^-2*a^-1*b*a^-1*b^-1*a^-1*b^-1,
      a^2*b^-1*a^-1*b^-1*a*c*b*c*a*b*a*b, a^3*b*a^2*b*a^-2*c*a*b*a^-1*c*a
     ] |

Before  it is returned,  the resulting presentation  is being simplified by
appropriate  calls  of  the  function  `simplify`  (see "Tietze
Transformations"),  but without  allowing it  to eliminate  any generators.
This  restriction guarantees  that we  get a  bijection between the list of
generators of `G` and the list of generators in the presentation. Hence, if
the  generators  of  `G`  are  redundant  and  if  you  don't  care for the
bijection, it may be convenient to apply the function
`simplify` again.

julia> P=Presentation("
1: aaa=1
2: AFB=1
3: bbb=1
4: ccc=1
5: CEd=1
6: ddd=1
7: af=ba
8: ac=ca
9: bc=cb
10: ad=da
11: bd=db
12: dc=cE
")
Presentation: 6 generators, 12 relators, total length 42
    gap> simplify( P );
    &I  there are 4 generators and 10 relators of total length 36 |

PresentationViaCosetTable(G[,F,words]) . . . . .  
construct a presentation for the given group

`PresentationViaCosetTable`   constructs   a  presentation  for  a  given
concrete  group.  It applies  John  Cannon's  relations finding algorithm
which has been described in

  John J. Cannon:  Construction of  defining relators  for  finte groups.
  Discrete Math. 5 (1973), 105-129, and in

  Joachim Neubueser: An elementary introduction to coset table methods in
  computational  group theory.  Groups-St. Andrews 1981,  edited by C. M.
  Campbell and E. F. Robertson, pp. 1-45.  London Math. Soc. Lecture Note
  Series no. 71, Cambridge Univ. Press, Cambridge, 1982.

If only a group  `G`  has been  specified,  the single stage algorithm is
applied.

If the  two stage algorithm  is to  be used,  `PresentationViaCosetTable`
expects a subgroup `H` of `G` to be described by two additional arguments
`F`  and  `words`,  where  `F`  is a  free group  with the same number of
generators as  `G`,  and  `words` is a list of words in the generators of
`F`  which supply  a list of generators of  `H`  if they are evaluated as
words in the corresponding generators of `G`.
"""
function PresentationViaCosetTable(G,arg...)
  #  F,          # given free group
  #  fgens,      # generators of F
  #  fwords,     # given words in the generators of F
  #  words,      # tidied up words for the subgroup generators
  #  H,          # subgroup
  #  elts,       # elements of G or H
  #  cosets,     # right cosets of G with respect to H
  #  F1,         # f.p. group isomorphic to H
  #  F2,         # f.p. group isomorphic to G
  #  P,          # resulting presentation
  #  ggens,      # concrete generators of G
  #  hgens,      # concrete generators of H
  #  ngens,      # number of generators of G
  if !(IsRec(G) && (haskey(G, :isGroup) && G[:isGroup]))
      error("first argument must be a group")
  end
  ggens=G[:generators]
  ngens=length(ggens)
  if length(arg)==0
    InfoFpGroup1("#I  calling the single stage relations finding algorithm\n")
    elts=elements(G)
    F2=RelsViaCosetTable(G, elts)
  else
    F=arg[1]
    if !(IsRec(F) && (haskey(F, :isFpGroup) && F[:isFpGroup]))
        error("second argument must be a free group")
    end
    fgens=F[:generators]
    if length(fgens) != ngens
        error("the given groups have different number of generators")
    end
    fwords=arg[2]
    hgens=map(w->MappedWord(w, fgens, ggens), fwords)
    H=Subgroup(G, hgens)
    if Size(H)==1 || Size(H)==Size(G)
      InfoFpGroup1("#I  calling the single stage relations finding algorithm\n")
      elts=elements(G)
      F2=RelsViaCosetTable(G, elts)
    else
      InfoFpGroup1("#I  ","calling the two stage relations finding algorithm\n",
                   "#I  using a subgroup of size ",Size(H),"  &&  index ",
                   Size(G) // Size(H), "\n")
      words=[]
      hgens=[]
      for w=fwords
          h=MappedWord(w, fgens, ggens)
          if h != G.identity && (!h in hgens && !(h ^ -1) in hgens)
              push!(words, w)
              push!(hgens, h)
          end
      end
      hgens=map(w->MappedWord(w, fgens, ggens), words)
      elts=elements(H)
      F1=RelsViaCosetTable(H, elts, hgens)
      cosets=RightCosets(G, H)
      F2=RelsViaCosetTable(G, cosets, F, words, H, F1)
    end
  end
  P=Presentation(F2, 0)
  P.generatorsLimit=ngens
  GoGo(P)
  P.generatorsLimit=0
  P.debug=1
  return P
end

"""
`RemoveRelator(P,n)`

`RemoveRelator`  removes the `n`th relator and  then resorts the  list of
relators in the given presentation `P`.
"""
function RemoveRelator(T::Presentation, n)
  if !(n in eachindex(T.relators)) error("relator number out of range") end
  rels=T.relators
  debug(T,3,"RemoveRelator removing $n-th relator")
  leng=length(T.relators[n])
  if leng==2 && rels[n][1]==rels[n][2]
    num=rels[n][1]
    if num<0 num=-num end
    if T[-num]==num T[-num]=-num end
  end
  rels[n]=Int[]
  sort!(T)
  PrintStatus(T,2)
  tracing(T,false)
end

"""
`simplify(G)`

`simplify`  applies  Tietze  transformations  to a  copy of  the
presentation of the given finitely presented group `G` in order to reduce
it with respect to  the number of generators, the number of relators, and
the relator lengths.

`simplify` returns the resulting finitely presented group (which
is isomorphic to `G`).

```julia-repl
julia> @AbsWord a,b,c,d,e,f

julia> F=FpGroup([a,b,c,d,e,f])
FreeGroup(a,b,c,d,e,f)

julia> G=F/[a^2,b^2,d*f^-1,e^2,f^2,a*b^-1*c,a*e*c^-1,b*d^-1*c,c*d*e^-1,a*f*c^-2,c^4]
FreeGroup(a,b,c,d,e,f)/[a²,b²,df⁻¹,e²,f²,ab⁻¹c,aec⁻¹,bd⁻¹c,cde⁻¹,afc⁻²,c⁴]

julia> simplify(G)
FreeGroup(a,c)/[a²,ac⁻¹ac⁻¹,c⁴]
```

In fact, the call

```julia-rep1
julia> simplify(G)
```

is an abbreviation of the call sequence

```julia-rep1
julia> P=Presentation(G,0);simplify(P);FpGroup(P)
```

which applies  a rather simple-minded strategy of  Tietze transformations
to the intermediate presentation `P`.
If for  some  concrete group the resulting presentation  is unsatisfying,
then  you  should  try  a  more  sophisticated,  interactive  use of  the
available Tietze transformation functions  (see "Tietze Transformations").
"""
function simplify(G::FpGroup)
  T=Presentation(G, 0)
  GoGo(T)
  FpGroup(T)
end

"""
`Check(presentation)`

checks some fields of the  given presentation to be consistent.
"""
function Check(P::Presentation)
  if length(P.inverses)!=2*length(P.generators)+1
    error("inconsistent inverses")
  end
  if length(P.flags)!=length(P.relators) error("inconsistent relators") end
end

"""
`Eliminate(P[, gen])`
`Eliminate(P[, n])`

`Eliminate`  tries to eliminate one or  more generators from a presentation
`P` via Tietze transformations.

Any  relator  which  contains  some  generator  just  once  can  be used to
substitute  that generator by  a word in  the remaining generators. If such
generators  and relators  exist, then  `Eliminate` chooses  a generator for
which  the  product  of  its  number  of  occurrences and the length of the
substituting  word is minimal,  and then it  eliminates this generator from
the  presentation, provided that the resulting total length of the relators
does  not exceed  the associated  Tietze option parameter `P.spaceLimit`.
The  default value of `P.spaceLimit` is  `infinity`, but you may alter it
appropriately (see Tietze options below).

If you specify a generator `gen` as  second argument, then  `Eliminate`
only tries to eliminate that generator.

If you  specify an  integer  `n`  as second argument,  then `Eliminate`
tries  to   eliminate   up  to   `n`  generators.  Note  that  the  calls
`Eliminate( `P` )` and `Eliminate( `P`, 1 )` are equivalent.
"""
function Eliminate(T::Presentation,n=nothing)
  debug(T,3,"Eliminate n=$n")
  Check(T)
  tietze=T
  gennum=0
  if n===nothing n=1
  elseif !(n isa Integer)
    gennum=findfirst(==(n),tietze.generators)
    if gennum===nothing
      error("usage: Eliminate(presentation[,nbgens or gen])")
    end
  end
  if n==1
    if gennum==0 EliminateGen1(T)
    else EliminateGen(T, gennum)
    end
    PrintStatus(T)
    if tietze.numredunds>0 RemoveGenerators(T) end
    HandleLength1Or2Relators(T)
    sort!(T)
    PrintStatus(T)
  else
    eliminations=T.eliminationsLimit;T.eliminationsLimit=n
    numgenslimit=T.generatorsLimit;T.generatorsLimit=length(tietze.generators)-n
    EliminateGens(T)
    T.eliminationsLimit=eliminations
    T.generatorsLimit=numgenslimit
  end
end

# substitute gen by w
function SubstituteGen(T::Presentation, gen, w)
  iw=-reverse(w)
  for i in eachindex(T.relators)
    res=Int[]
    occ=0
    for j in T.relators[i]
      if j==gen append!(res,w);occ+=1
      elseif j==-gen append!(res,iw);occ+=1
      else push!(res,j)
      end
    end
    if occ>0 T.flags[i]=1 end
    if length(res)==2 && allequal(res) T.relators[i]=res
    else T.relators[i]=reduceword(res,T)
    end
  end
  UpdateGeneratorImages(T, abs(gen), gen<0 ? iw : w)
end

"""
EliminateGen( presentation,n) . . . eliminates the nth generator

`EliminateGen` eliminates the Tietze generator tietze.generators[n]
if possible, i. e. if that generator can be isolated  in some appropriate
Tietze relator.  However,  the elimination  will not be  performed if the
resulting total length of the relators cannot be guaranteed to not exceed
the parameter P.lengthLimit (if !=0).
"""
function EliminateGen(P::Presentation, num::Int)
  spacelimit=sum(length,P.relators)+P.lengthLimit
  P.modified=false
  if !(num in eachindex(P.generators))
    error("EliminateGen: second argument is not a valid generator number")
  end
  occTotal,occRelNum,occMult=Occurrences(P, num)
  if occMult==1
    length=length(P.relators[occRelNum])
    space=(occTotal-1)*(length-1)-length
    if P.lengthLimit!=0 && sum(length,P.relators)+space>spacelimit return end
    gen=num
    rel=P.relators[occRelNum]
    pos=findfirst(==(gen),rel)
    if pos===nothing
      gen=-gen
      pos=findfirst(==(gen),rel)
    end
    word=vcat(rel[pos+1:end],rel[1:pos-1])
    debug(P,2,"EliminateGen eliminating ", P.generators[num], "=",
         let a=gen>0 ? -reverse(word) : word;  alphab(word,P.generators) end)
    SubstituteGen(P,-gen, word)
    P[num]=0
    P.numredunds+=1
    P.modified=true
  end
end

# returns 3 vectors:
# 1: nb. occurences of each generator
# 2: n0 of a relator of minimal length where occurs with minimal>0 multiplicity
# 3: the corresponding multiplicity
function Occurrences(T::Presentation)
  occ=fill(0,length(T.generators))
  min=fill(0,length(T.generators))
  mult=fill(0,length(T.generators))
  for j in eachindex(T.relators)
    lococc=fill(0,length(T.generators))
    for i in T.relators[j] lococc[abs(i)]+=1 end
    occ+=lococc
    for i in eachindex(T.generators)
      if lococc[i]>0
      if mult[i]==0 || lococc[i]<mult[i] || 
        (lococc[i]==mult[i] && length(T.relators[j])<length(T.relators[min[i]]))
          min[i]=j
          mult[i]=lococc[i]
        end
      end
    end
  end
  [occ,min,mult]
end

Occurrences(T::Presentation,gen)=getindex.(Occurrences(T),gen)

"""
EliminateGen1(presentation)  . . . . . . .  eliminates a generator

`EliminateGen1`  tries to  eliminate a  Tietze generator:  If there are
Tietze generators which occur just once in certain Tietze relators,  then
one of them is chosen  for which the product of the length of its minimal
defining word  and the  number of its  occurrences  is minimal.  However,
the elimination  will not be performed  if the resulting  total length of
the  relators   cannot  be  guaranteed   to  not  exceed   the  parameter
P.lengthLimit (if !=0).
"""
function EliminateGen1(T::Presentation)
  occTotals,occRelNo,occMultiplicities=Occurrences(T)
  modified=false
  num=0
  space=0
  for i in T.protected+1:length(T.generators)
    if occMultiplicities[i]==1
      len=length(T.relators[occRelNo[i]])
      ispace=(occTotals[i]-1)*(len-1)-len
      if num==0 || ispace<=space
        num=i
        space=ispace
      end
    end
  end
  if num>0 && (T.lengthLimit==0 || space<=T.lengthLimit)
    gen=num
    occRelNo=occRelNo[num]
    rel=T.relators[occRelNo]
    len=length(T.relators[occRelNo])
    pos=findfirst(==(gen),rel)
    if pos===nothing
      gen=-gen
      pos=findfirst(==(gen),rel)
    end
#   display_balanced(T)
#   @show rel
    word=vcat(rel[pos+1:end],rel[1:pos-1])
    debug(T,2,"EliminateGen1 eliminating ",T.generators[num],"=",
      let a=gen>0 ? -reverse(word) : word;  alphab(word,T.generators) end)
    SubstituteGen(T, -gen, word)
    T[num]=0
    T.numredunds+=1
    modified=true
  end
  T.modified=modified
end

"""
EliminateGens(P::presentation)

`EliminateGens`  repeatedly eliminates generators  from the presentation of
the given group until one of the following conditions is violated:

(1) The  current  number of  generators  is greater than P.generatorsLimit.
(2) The number of generators eliminated so far is less than P.eliminationsLimit.
(3) The  total length of the relators  has not yet grown  to a percentage
    greater than the parameter P.expandLimit.
(4) The  next  elimination  will  not  extend the total length to a value
    greater than the parameter P.lengthLimit.

the function will not eliminate any protected generators.
"""
function EliminateGens(P::Presentation)
  debug(P,3,"EliminateGens eliminating generators")
  redundantsLimit=5
  bound=sum(length,P.relators)*P.expandLimit//100
  modified=false
  P.modified=true
  num=0
  while P.modified && num<P.eliminationsLimit && sum(length,P.relators)<=bound &&
   length(P.generators)-P.numredunds>P.generatorsLimit
    EliminateGen1(P)
#   print("a ");PrintStatus(P);@show P.flags
    if P.numredunds==redundantsLimit 
      RemoveGenerators(P)
#     print("b ");PrintStatus(P);@show P.flags
    end
    modified=modified || P.modified
    num+=1
  end
  P.modified=modified
  if P.numredunds>0 
    RemoveGenerators(P)
#   print("c ");PrintStatus(P);@show P.flags
  end
  if modified
    HandleLength1Or2Relators(P)
#   print("d ");PrintStatus(P);@show P.flags
    sort!(P)
    PrintStatus(P,2)
  end
end

"""
`FindCyclicJoins(P)`

`FindCyclicJoins` performs Tietze transformations on a presentation `P`. It
searches  for pairs of  generators which generate  the same cyclic subgroup
and eliminates one of the two generators of each such pair it finds.

More  precisely: `FindCyclicJoins` searches for pairs of generators `a` and
`b`  such that (possibly after inverting  or conjugating some relators) the
set  of  relators  contains  the  commutator  `[a,b]`, a power `a^n`, and a
product  of the form `a^s  b^t` with `s` prime  to `n`. For each such pair,
`FindCyclicJoins` uses the Euclidian algorithm to express `a` as a power of
`b`, and then it eliminates `a`.
"""
function FindCyclicJoins(P::Presentation)
  debug(P,3,"FindCyclicJoins searching for cyclic joins")
  T=P
  T.modified=false
  newstart=true
  while newstart
    exponents=find_exponents(P)
    if iszero(exponents) return end
    newstart=false
    gens=T.generators
    numgens=length(T.generators)
    rels=T.relators
    numrels=length(rels)
    flags=T.flags
    i=0
    while i<numrels
      i+=1
      rel=rels[i]
      if length(T.relators[i])!=4 || rel[1]!=T[-rel[3]] || rel[2]!=T[-rel[4]] continue end
      num=abs.(rel[1:2])
      exp=exponents[num]
      if iszero(exp) continue end
      fac=[0, 0]
      e=[0, 0]
      j=0
      while j<numrels
        j+=1
        if length(T.relators[j])==0 || j==i continue end
        if !all(x->abs(x) in num,rels[j]) continue end
        rel=rels[j]
        e.=0
        for x in rel, i in 1:2
          if x==num[i] e[i]+=1 elseif x==-num[i] e[i]-=1 end
        end
        powers=count(!iszero,e)
        if powers<2 continue end
        for k in 1:2
          fac[k]=num[k]
          if exp[k]>0
            e[k]=mod(e[k], exp[k])
            if 2e[k]>exp[k]
              e[k]=exp[k]-e[k]
              fac[k]=-fac[k]
            end
          elseif e[k]<0
            e[k]=-e[k]
            fac[k]=-fac[k]
          end
          if fac[k]<0 fac[k]=T[fac[k]] end
        end
        for k in 1:2
          if e[k]>0 && e[3-k]==0
            exp[k]=gcd(e[k], exp[k])
            if exp[k] != exponents[num[k]]
              exponents[num[k]]=exp[k]
              e[k]=exp[k]
            end
          end
        end
        rel=vcat(fill(fac[1],e[1]), fill(fac[2],e[2]))
        if rels[j]!=rel
          rels[j]=rel
          flags[j]=1
          T.modified=true
          debug(P,3,"FindCyclicJoins rels[$j] reduced to ",rels[j])
        end
        if e[1]==1 n=num[1]
        elseif e[2]==1 n=num[2]
        else
          n=0
          for k in 1:2
            if n==0 && e[k]>1 && gcd(e[k], exp[k])==1
              gen=gens[num]
              if fac[1]<0 gen[1]=inv(gen[1]) end
              if fac[2]<0 gen[2]=inv(gen[2]) end
              word=gen[k]*gen[3-k]^(e[3-k]*gcdx(e[k], exp[k])[2])
              AddRelator(P, word)
              numrels=length(T.relators)
              n=num[k]
            end
          end
        end
        if n!=0 && P.generatorsLimit<numgens
          Eliminate(P)
          T.modified=true
          j=numrels
          i=numrels
          if 1<numgens newstart=true end
        end
      end
    end
  end
  if T.modified
    HandleLength1Or2Relators(P)
    sort!(P)
    PrintStatus(P,1)
  end
end

"""
`Go(P::Presentation[,silent])`

`Go`  performs Tietze transformations on a  presentation `P`. It is perhaps
the  most convenient of the  interactive Tietze transformation commands. It
offers  a  kind  of  default  strategy  which,  in  general, saves you from
explicitly calling the lower-level commands it involves.

Roughly  speaking, `Go` consists of a  loop over a procedure which involves
two  phases:  In  the  *search  phase*  it calls `Search` and `SearchEqual`
described  below which  try to  reduce the  relator lengths by substituting
common  subwords  of  relators,  in  the  *elimination  phase* it calls the
command  `Eliminate` described below  (or, more precisely,  a subroutine of
`Eliminate`  in order to save some  administrative overhead) which tries to
eliminate  generators  that  can  be  expressed  as  words in the remaining
generators.

If  `Go`  succeeds  in  reducing  the  number  of generators, the number of
relators,  or the total  length of all  relators, then it  displays the new
status  before returning (provided that you did  not set the print level to
zero).  However, it does not  provide any output if  all these three values
have  remained unchanged,  even if  the `SearchEqual`  command involved has
changed  the  presentation  such  that  another  call of `Go` might provide
further  progress. Hence, in such a case  it makes sense to repeat the call
of  the command for  several times (or  to call instead  the `GoGo` command
which we will describe next).

As an example  we  compute  a presentation of a  subgroup of index 408 in
`PSL(2,17)`.

|    gap> F2 := FreeGroup( "a", "b" );;
    gap> G := F2 / [ F2.1^9, F2.2^2, (F2.1*F2.2)^4, (F2.1^2*F2.2)^3 ];;
    gap> a := G.1;;  b := G.2;;
    gap> H := Subgroup( G, [ (a*b)^2, (a^-1*b)^2 ] );;
    gap> Index( G, H );
    408
    gap> P := PresentationSubgroup( G, H );
    << presentation with 8 gens and 36 rels of total length 111 >>
    gap> DisplayPresentation(P);
    1: a=A
    2: g=G
    3: c=C
    4: f=F
    5: h=H
    6: d=D
    7: fbc=1
    8: dBa=1
    9: gdb=1
    10: ceb=1
    11: eab=1
    12: caB=1
    13: BBB=1
    14: cBd=1
    15: acb=1
    16: afB=1
    17: bbb=1
    18: bac=1
    19: cba=1
    20: aBc=1
    21: haB=1
    22: dab=1
    23: Bec=1
    24: Bad=1
    25: Bca=1
    26: bbb=1
    27: afB=1
    28: agb=1
    29: agb=1
    30: ab=BA
    31: ab=bC
    32: ab=BA
    33: ac=Ce
    34: ac=bb
    35: cB=bC
    36: aca=CAC
    gap> P.primaryGeneratorWords;
    [ b, a*b*a ]
    gap> P.protected := 2;;
    gap> P.debug := 2;;
    gap> simplify( P );
    &I  eliminating _x7 = _x5
    &I  eliminating _x5 = _x4
    &I  eliminating _x18 = _x3
    &I  eliminating _x8 = _x3
    &I  there are 4 generators and 8 relators of total length 21
    &I  there are 4 generators and 7 relators of total length 18
    &I  eliminating _x4 = _x3^-1*_x2^-1
    &I  eliminating _x3 = _x2*_x1^-1
    &I  there are 2 generators and 4 relators of total length 14
    &I  there are 2 generators and 4 relators of total length 13
    &I  there are 2 generators and 3 relators of total length 9
    gap> relators( P );
    &I  1. _x1^2
    &I  2. _x2^3
    &I  3. _x2*_x1*_x2*_x1 |

Note that the number  of loops over the two phases as well  as the number
of  subword   searches  or  generator  eliminations  in  each  phase  are
determined by a set of option parameters which may  heavily influence the
resulting presentation and the computing time (see Tietze options below).

`Go` is just another name  for the  `simplify` command.  It
has  been introduced for  the convenience  of those {GAP} users  who are
used to  that  name from the *go* option of the ANU Tietze transformation
stand-alone program or from the *go* command in SPAS.

If "silent" is specified as true, then the printing of the status line by
`Go` in case of T.debug=1 is suppressed.
"""
function Go(P::Presentation,silent=false)
  printstatus=P.debug==1 && !silent
  Check(P)
  Search(P)
  count=0
  while (P.loopLimit==0 || count<P.loopLimit) && !all(isempty,P.relators)
#   println("Go loop")
    SearchEqual(P)
#   @show P.flags
    if P.modified Search(P) end
#   @show P.flags
    EliminateGens(P)
#   print("before search");dump(P)
    if P.modified
      Search(P)
      count+=1
    end
#   print("after search");dump(P)
    if printstatus print("#");PrintStatus(P) end
    if !P.modified break end
  end
  if !all(isempty,P.relators)
    FindCyclicJoins(P)
    if P.modified Search(P) end
    if printstatus print("#");PrintStatus(P) end
  end
  P
end

"""
`GoGo(P)`

`GoGo` performs Tietze transformations on a presentation `P`. It repeatedly
calls  the  `Go`  command  until  neither  the number of generators nor the
number of relators nor the total length of all relators have changed during
five consecutive calls of `Go`.
"""
function GoGo(T::Presentation)
  numgens=length(T.generators)
  numrels=length(T.relators)
  total=sum(length,T.relators)
  silentGo=T.debug==1
  count=0
  while count<5
    Go(T, silentGo)
    count+=1
    if silentGo && (length(T.generators)<numgens || length(T.relators)<numrels)
        PrintStatus(T)
    end
    if length(T.generators)<numgens || length(T.relators)<numrels || 
      sum(length,T.relators;init=0)<total
      numgens=length(T.generators)
      numrels=length(T.relators)
      total=sum(length,T.relators;init=0)
      count=0
    end
  end
  if silentGo PrintStatus(T) end
  T
end

"""
OccurrencesPairs(T::Presentation,gen)
counts occurrences of pairs gen,geni or their inverses (where 0<gen<geni)
circularly in the relators or their inverses
returns a 4×numgens matrix; each line is zero until gen+1
the first line  counts the occurrences of (gen,geni)
the second             the occurrences of (gen,-geni)
the third              the occurrences of (-gen,geni)
the last               the occurrences of (-gen,-geni)
"""
function OccurrencesPairs(T::Presentation,gen)
  res=fill(0,4,length(T.generators))
  for rel in T.relators, i in eachindex(rel)
    if rel[i]==gen add=0
    elseif rel[i]==-gen add=2
    else continue
    end
    geni=rel[i==end ? 1 : i+1]
    if geni>gen res[1+add,geni]+=1
    elseif geni<-gen res[2+add,-geni]+=1
    end
    geni=rel[i==1 ? end : i-1]
    if geni>gen res[4-add,geni]+=1
    elseif geni<-gen res[3-add,-geni]+=1
    end
  end
  res
end

"""
MostFrequentPairs(presentation, n) ....  occurrences of pairs

`MostFrequentPairs`  returns a list  describing the  n  most frequently
occurring relator subwords of the form  g1*g2,  where  g1  and  g2 are
different generators or their inverses. n=0 interpreted as infinity
"""
function MostFrequentPairs(T::Presentation, nmax::Int)
  if nmax<0 error("second argument must be≥0") end
  gens=T.generators
  numgens=length(T.generators)
  pairs=Vector{Int}[]
  n=0
  if nmax==1
    max=0
    for i in 1:numgens-1
      occlist=OccurrencesPairs(T, i)
      for j in i+1:numgens, k in 4:-1:1
        if occlist[k,j]>=max
          max=occlist[k,j]
          if isempty(pairs) push!(pairs,[max, i, j, k-1])
          else pairs[1]=[max, i, j, k-1]
          end
          n=1
        end
      end
    end
  else
    for i in 1:numgens-1
      occlist=OccurrencesPairs(T, i)
      for j=i+1:numgens, k in 1:4
        if occlist[k,j]>0
          n+=1
          if length(pairs)==n-1 push!(pairs,[occlist[k,j], i, j, k-1])
          else pairs[n]=[occlist[k,j], i, j, k-1]
          end
        end
      end
      if nmax!=0 && n>nmax 
        sort!(pairs,rev=true)
        resize!(pairs,nmax)
        n=nmax
      end
    end
    sort!(pairs,rev=true)
  end
  pairs
end

"""
`PrintOptions(P)`

Several   of  the  Tietze  transformation   commands  described  above  are
controlled  by certain parameters, the *Tietze options*, which often have a
tremendous  influence on  their performance  and results.  However, in each
application  of  the  commands,  an  appropriate  choice  of  these  option
parameters  will depend  on the  concrete presentation under investigation.
Therefore  we have implemented the  Tietze options in such  a way that they
are  associated to the presentation: Each presentation keeps its own set of
Tietze option parameters in the form of ordinary fields. In particular, you
may  alter the value of any of these Tietze options by just assigning a new
value to the respective component.

`PrintOptions`  prints the Tietze  option fields of  the  specified
presentation `P`.

The Tietze options have the following meaning.

`protected`: 
  The  first `P.protected` generators  in a presentation  `P` are protected
  from  being eliminated by the Tietze transformations functions. There are
  only two exceptions: The option `P.protected` is ignored by the functions
  `Eliminate(P,gen)` and `Substitute(P,n,eliminate)`    because   they
  explicitly  specify the generator to be  eliminated. The default value of
  `protected` is 0.

`eliminationsLimit`: 
  Whenever  the  elimination  phase  of  the  `Go` command is entered for a
  presentation  `P`, then  it will  eliminate at most `P.eliminationsLimit`
  generators (except for further ones which have turned out to be trivial).
  Hence  you may use the `eliminationsLimit` parameter as a break criterion
  for  the  `Go`  command.  Note,  however,  that  it  is  ignored  by  the
  `Eliminate` command. The default value of `eliminationsLimit` is 100.

`expandLimit`: 
  Whenever the routine for eliminating more than 1 generators is called for
  a presentation `P` by the `Eliminate` command or the elimination phase of
  the  `Go` command, then it saves the  given total length of the relators,
  and  subsequently it  checks the  current total  length against its value
  before  each elimination. If the total  length has increased to more than
  `P.expandLimit`  per cent of its original value, then the routine returns
  instead   of  eliminating  another  generator.  Hence  you  may  use  the
  `expandLimit`  parameter as a  break criterion for  the `Go` command. The
  default value of `expandLimit` is 150.

`generatorsLimit`: 
  Whenever  the  elimination  phase  of  the  `Go` command is entered for a
  presentation  `P` with  `n` generators,  then it  will eliminate  at most
  `n-P.generatorsLimit` generators (except for generators which turn out to
  be trivial). Hence you may use the `generatorsLimit` parameter as a break
  criterion for the `Go` command. The default `generatorsLimit` is 0.

`lengthLimit`: 
  The  Tietze transformation commands will never eliminate a generator of a
  presentation  `P`,  if  they  cannot  exclude  the  possibility  that the
  resulting   total   length   of   the   relators  exceeds  the  value  of
  `P.lengthLimit`. The default value of `lengthLimit` is `infinity`.

`loopLimit`: 
  Whenever  the `Go` command is called for a presentation `P`, then it will
  loop over at most `P.loopLimit` of its basic steps. Hence you may use the
  `loopLimit`  parameter as  a break  criterion for  the `Go`  command. The
  default value of `loopLimit` is `infinity`.

`debug`: 
  Whenever Tietze transformation commands are called for a presentation `P`
  with  `P.debug` `= 0`, they will not  provide any output except for error
  messages. If `P.debug` `= 1`, they will display some reasonable amount of
  output  which allows you to watch the  progress of the computation and to
  decide  about your next commands.  In the case `P.debug`  `= 2`, you will
  get  a much more generous amount of  output. Finally, if `P.debug` `= 3`,
  various  messages on internal details will be added. The default value of
  `debug` is 1.

`saveLimit`: 
  Whenever  the  `Search`  command  has  finished  its  main  loop over all
  relators  of a presentation `P`, then  it checks whether during this loop
  the   total  length  of  the  relators  has  been  reduced  by  at  least
  `P.saveLimit`  per cent. If  this is the  case, then `Search` repeats its
  procedure  instead  of  returning.  Hence  you  may  use  the `saveLimit`
  parameter  as  a  break  criterion  for  the  `Search`  command  and,  in
  particular,  for the search phase of  the `Go` command. The default value
  of `saveLimit` is 10.

`searchSimultaneous`: 
  Whenever  the  `Search`  or  the  `SearchEqual`  command  is called for a
  presentation `P`, then it is allowed to handle up to `searchSimultaneous`
  short  relators simultaneously (see  for the description  of the `Search`
  command  for  more  details).  The  choice  of this parameter may heavily
  influence  the performance as well as the  result of the `Search` and the
  `SearchEqual`  commands and  hence also  of the  search phase of the `Go`
  command. The default value of `searchSimultaneous` is 20.

As soon as a presentation has been defined, you may alter any of its Tietze
option  parameters  at  any  time  by  just  assigning  a  new value to the
respective component.

To  demonstrate the  effect of  the `eliminationsLimit`  parameter, we will
give  an example in which we  handle a subgroup of index  240 in a group of
order  40320  given  by  a  presentation  due  to  B.~H.  Neumann. First we
construct  a presentation  of the  subgroup, and  then we  apply to  it the
`GoGo`  command for  different values  of the `eliminationsLimit` parameter
(including  the  default  value  100).  In  fact, we also alter the `debug`
parameter,  but this is only done in  order to suppress most of the output.
In  all cases  the resulting  presentations cannot  be improved any more by
applying the `GoGo` command again, i.e., they are the best results which we
can get without substituting new generators.
```julia-rep1
    gap> F3 := FreeGroup( "a", "b", "c" );;
    gap> G := F3 / [ F3.1^3, F3.2^3, F3.3^3, (F3.1*F3.2)^5,
            (F3.1^-1*F3.2)^5, (F3.1*F3.3)^4, (F3.1*F3.3^-1)^4,
            F3.1*F3.2^-1*F3.1*F3.2*F3.3^-1*F3.1*F3.3*F3.1*F3.3^-1,
            (F3.2*F3.3)^3, (F3.2^-1*F3.3)^4 ];;
    gap> H := Subgroup( G, [ G.1, G.3 ] );;
    gap> P := PresentationSubgroup( G, H );
    << presentation with 224 gens and 593 rels of total length 2769 >>
    julia> for i in [ 28, 29, 30, 94, 100 ]
    Pi=deepcopy(P)
    Pi.eliminationsLimit=i
    println("# eliminationsLimit set to ",i)
      Pi.debug=0
      Presentations.GoGo(Pi)
      print(Pi)
    end
    #  eliminationsLimit set to 28
    #  there are 2 generators and 95 relators of total length 10817
    #  eliminationsLimit set to 29
    #  there are 2 generators and 5 relators of total length 35
    #  eliminationsLimit set to 30
    #  there are 3 generators and 98 relators of total length 2928
    #  eliminationsLimit set to 94
    #  there are 4 generators and 78 relators of total length 1667
    #  eliminationsLimit set to 100
    #  there are 3 generators and 90 relators of total length 3289 |
```

Similarly,  we demonstrate  the influence  of the  `saveLimit` parameter by
just  continuing the  preceding example  for some  different values  of the
`saveLimit`  parameter  (including  its  default  value  10),  but  without
changing  the `eliminationsLimit`  parameter which  keeps its default value
100.

```julia-rep1
    julia> for i in [ 9, 10, 11, 12, 15 ]
       Pi=deepcopy(P)
       Pi.saveLimit=i
       println("# saveLimit set to ",i)
       Pi.debug=0
       GoGo(Pi)
       print(Pi)
    end
    &I  saveLimit set to 9
    &I  there are 3 generators and 97 relators of total length 5545
    &I  saveLimit set to 10
    &I  there are 3 generators and 90 relators of total length 3289
    &I  saveLimit set to 11
    &I  there are 3 generators and 103 relators of total length 3936
    &I  saveLimit set to 12
    &I  there are 2 generators and 4 relators of total length 21
    &I  saveLimit set to 15
    &I  there are 3 generators and 143 relators of total length 18326 |
```
"""
function PrintOptions(P::Presentation)
  names=[:eliminationsLimit,:expandLimit,:generatorsLimit,:lengthLimit,
         :loopLimit,:debug,:saveLimit,:searchSimultaneous]
  lst=sort(collect(intersect(keys(P.prop),names)))
  len=maximum(length.(string.(lst)))
  foreach(lst)do nam
    println(rpad(nam,len),"= ",getproperty(P,nam))
  end
end

"""
`RemoveGenerators(P::presentation)`

deletes  the redundant  Tietze generators  and renumbers  the non-redundant
ones  accordingly. The  redundant generators  are assumed  to be  marked in
`P.inverses` by an entry `P[i]!=i`.
"""
function RemoveGenerators(P::Presentation)
  debug(P,3,"RemoveGenerators renumbering the Tietze generators")
  redunds=P.numredunds
  if redunds==0 return end
  gens=P.generators
  numgens=length(P.generators)
  j=0
# @show j,P.inverses, P.numredunds
  for i in 1:numgens
    if P[i]==i
      j+=1
      if j<i
        P[i]=j
        if P[-i]>0 P[-i]=j else P[-i]=-j end
      end
    else
      P[i]=P[-i]=0
    end
  end
# @show j,P.inverses
  if j!=numgens-redunds
    error("j=$j numgens=$numgens redunds=$redunds. You should never get here.\n")
  end
# @show j, P.inverses, numgens
  for rel in P.relators, i in eachindex(rel) rel[i]=P[rel[i]] end
  if haskey(P, :imagesOldGens)
    for i in 1:length(P.imagesOldGens)
      image=P.imagesOldGens[i]
      newim=Int[]
      for j=1:length(image) push!(newim, P[image[j]]) end
      P.imagesOldGens[i]=reduceword(newim)
    end
  end
  for i in 1:numgens
    j=P[i]
    if j<i && j>0
      gens[j]=gens[i]
      P[j]=j
      P[-j]=P[-i]
      if haskey(P, :preImagesNewGens)
          P.preImagesNewGens[j]=P.preImagesNewGens[i]
      end
    end
  end
  j=numgens
  numgens1=numgens+1
  numgens-=redunds
  P.inverses=P.inverses[numgens1-numgens:numgens1+numgens]
  resize!(P.generators,numgens)
  if haskey(P, :preImagesNewGens) resize!(P.preImagesNewGens,numgens) end
  P.numredunds=0
end

# hash of cyclic subword of relator w starting at i and of length l
# i in 1:length(w) or i negative: backward subword
function hashword(w,i,l,T)
  res=UInt(0)
  if i<0
    for j in -i:-1:max(1,-i+1-l) res=hash(T[-w[j]],res) end
    for j in length(w):-1:length(w)-l-i+1 res=hash(T[-w[j]],res) end
  else
    for j in i:min(i+l-1,length(w)) res=hash(T[w[j]],res) end
    for j in 1:i+l-1-length(w) res=hash(T[w[j]],res) end
  end
  res
end

# finds intersection of 2 sorted lists; returns list of pairs (i1,i2)
# such that by(a[i1])==by(b[i2])
function commonsorted(a,b;by=x->x)
  res=Tuple{Int,Int}[]
  ai=bi=1
  ri=0
@inbounds while ai<=length(a) && bi<=length(b)
    c=cmp(by(a[ai]),by(b[bi]))
    if     c>0 bi+=1
    elseif c<0 ai+=1
    else push!(res,(ai,bi)); ai+=1; bi+=1
    end
  end
  res
end
  
# perhaps w1[abs(i1)]*sign(i1)=w2[abs(i2)]*sign(i2), extend the match
# as much as possible; returns (new i1, new i2, length of found match)
function maxmatch(w1,i1,w2,i2,T)
  l=0;lw=length(w1)
  prev(i,w)=i>0 ? (i==1 ? length(w) : i-1) : (i==-length(w) ? -1 : i-1)
  next(i,w)=i>0 ? (i==length(w) ? 1 : i+1) : (i==-1 ? -length(w) : i+1)
  eq(i1,i2)=T[w1[abs(i1)]*sign(i1)]==T[w2[abs(i2)]*sign(i2)]
  n1=i1;n2=i2
  while l<lw
    if !eq(n1,n2) break end
    j1=n1;j2=n2;l+=1
    n1=next(j1,w1);n2=next(j2,w2)
  end
  if l>0
  while l<lw
    n1=prev(i1,w1);n2=prev(i2,w2)
    if !eq(n1,n2) break end
    i1=n1;i2=n2;l+=1
  end
  end
  (i1,i2,l)
end
   
# cyclic subword of relator w starting at i and of length l; takes i mod
# length(w) but remembers sign (positive: forward word, negative: -backward)
function subword(w,i,l)
  neg=sign(i)
  i=mod1(abs(i),length(w))
  if neg<0 vcat(-w[i:-1:max(1,i+1-l)],-w[length(w):-1:length(w)-l+i+1])
  else vcat(w[i:min(i+l-1,length(w))],w[1:i+l-1-length(w)])
  end
end

# doreplace(w1,w2,pos1,pos2,l) replace in relator w2 substring of length l
# starting at pos2 with complement in w1 of identical substring of length l
# starting at pos1 in w1.
# return new w2
function doreplace(w1,w2,pos1,pos2,l)
  l1=length(w1)
  compl=subword(w1,pos1>0 ? mod(1-pos1,l1)-l1 : mod1(-pos1+1,l1),l1-l)
  vcat(compl,subword(w2,pos2+l,length(w2)-l))
# circshift(vcat(compl,subword(w2,pos2+l,length(w2)-l)),pos2-1)
end

alphab(rel)=String(map(i->i>0 ? Char(i+96) : uppercase(Char(-i+96)),rel))

function alphab(rel,gens)
  prod(i->i>0 ? string(gens[i]) : uppercase(string(gens[-i])),rel;init="")
end

# search from i-th relator, hashing subwords until j-th
function SearchC(T::Presentation,i,j,equal=false)
# @show i,j,equal
  rels=T.relators
# println("SearchC(",join(alphab.(rels[i:j],Ref(T.generators)),","),",",equal,")")
  len=length(rels[i])
  if !equal l=div(len,2)+1;lmin=len-len&1;lmax=lmin+1
  elseif isodd(len) error("searchequal and odd length")
  else l=div(len,2);lmin=lmax=len
  end
# @show l,lmin,lmax
  hh=NamedTuple{(:hash, :pos, :relno), Tuple{UInt, Int, Int}}[]
# if T.inverses!=T.numgens:-1:-T.numgens error("numgens=",T.numgens," invs=",
#                    T.inverses) end
  altered=Int[]
  for k in i:length(rels)
    rel=rels[k]
    if length(rel)<lmin continue end # sort should be called before SearchC
    rel=reduceword(rel,T)
    while true
      hk=map(u->(hash=hashword(rel,u,l,T),pos=u),eachindex(rel))
      sort!(hk,by=x->x.hash)
      pp=commonsorted(hh,hk;by=x->x.hash)
      if isempty(pp) break end
      hk=map(pp)do (p1,p2)
        (pos1,pos2,w)=(hh[p1].pos,hk[p2].pos,hh[p1].relno)
        (match=maxmatch(rels[w],pos1,rel,pos2,T),relno1=w,relno2=k)
      end
#     @show hk
      (pos1,pos2,l1),w,k=hk[argmax(map(x->(x[1][3],-x[1][2]),hk))]
#     if l1<l println("!!!!hash",hk) end
      if l1<l error("!")
      else
#       print("$k:",alphab(rel),"+$w:(",alphab(rels[w]),")$pos1:$pos2($l1)")
        rel=doreplace(rels[w],rel,pos1,pos2,l1)
#       println("=>",alphab(rel))
      end
      rel=reduceword(rel,T)
#     println("=>",alphab(rel))
#     if T.flags[k]==1 T.flags[k]=2 end
#     if equal T.flags[k]=1 end
      if equal || length(rel)<l break end
    end
    if k<=j && length(rel) in lmin:lmax
      for u in 1:length(rel) 
        push!(hh,(hash=hashword(rel,u,l,T),pos=u,relno=k))
        push!(hh,(hash=hashword(rel,-u,l,T),pos=-u,relno=k))
      end
      sort!(hh,by=x->x.hash)
    end
    if rel!=rels[k] 
      push!(altered,k)
      rels[k]=rel
    end
  end
  altered
end

"""
`Search(P)`

`Search`  performs Tietze transformations on a presentation `P`. It tries
to  reduce the relator lengths by  substituting common subwords of relators
by shorter words.

The idea is to find pairs of relators `r₁` length `l₁` and `r₂` length `l₂`
with  `l₁≤l₂` such that `r₁` and `r₂` coincide (possibly after inverting or
conjugating  one of them)  on some maximal  subword `w`, of length `>l₁/2`,
and  then to substitute each copy of  `w` in `r₂` by the inverse complement
of `w` in `r₁`.

Two  of the Tietze  option parameters which  are listed at  the end of this
section  may  strongly  influence  the  performance  and the results of the
`Search`   command.   These   are   the   parameters   `P.saveLimit`  and
`P.searchSimultaneous`. The first of them has the following effect.

When  Search  has  finished  its  main  loop  over all relators, then, in
general,  there are relators which have changed and hence should be handled
again in another run through the whole procedure. However, experience shows
that  it really does  not pay to  continue this way  until no more relators
change.  Therefore,  `Search`  starts  a  new  loop only if the loop just
finished  has  reduced  the  total  length  of  the  relators  by  at least
`P.saveLimit` per cent.

The default value of `P.saveLimit` is 10.

To  understand the effect of  the parameter `P.searchSimultaneous`, we have
to look in more detail at how `Search` proceeds.

First,  it  sorts  the  list  of  relators  by  increasing lengths. Then it
performs  a loop  over this  list. In  each step  of this loop, the current
relator  is treated  as *short  relator* `r₁`,  and a  subroutine is called
which  loops over the succeeding relators, treating them as *long relators*
`r₂` and performing the respective comparisons and substitutions.

As  this  subroutine  performs  a  very  expensive  process,  it  has  been
implemented as SearchC. For the given relator `r₁` of length `l₁`, it first
determines  the *minimal  match length*  `l=⌊l₁/2⌋+1`. Then  it builds up a
hash  list for all  subwords of length  `l` occurring in  the conjugates of
`r₁`  or  `r₁⁻¹`,  and  finally  it  loops  over all long relators `r₂` and
compares the hash values of their subwords of length `l` against this list.
A  true comparison is  only done if a hash match has been found.

To improve the efficiency of this process we allow the subroutine to handle
several  short  relators  simultaneously  provided  that they have the same
minimal  match  length.  If,  for  example,  it  handles `n` short relators
simultaneously,  then you save `n-1` loops over the long relators `r₂`, but
you pay for it by additional fruitless subword comparisons. In general, you
will  not get the best performance  by always choosing the maximal possible
number of short relators to be handled simultaneously. In fact, the optimal
choice  of  the  number  will  depend  on  the  concrete presentation under
investigation.   You  can  use   the  parameter  `P.searchSimultaneous`  to
prescribe  an upper bound  for the number  of short relators  to be handled
simultaneously.

The default value of `P.searchSimultaneous` is 20.
"""
function Search(P::Presentation)
  debug(P,3,"Search subwords")
  Check(P)
  T=P
  rels=T.relators
  T.modified=false
  save=P.saveLimit//100
  while sum(length,T.relators;init=0)>0
    sort!(P)
    modified=false
    oldtotal=sum(length,T.relators)
    flag=0
    for i in eachindex(T.relators)
      if T.flags[i]>1 || isempty(rels[i]) continue end
      len=length(rels[i])
      lmax=iseven(len) ? len+1 : len
      if flag<T.flags[i] flag=T.flags[i] end
      simultan=1
      j=i
      lastj=0
      k=i+1
      while k<=length(T.relators) && length(rels[k])<=lmax && simultan<P.searchSimultaneous
        if T.flags[k]<=1
          lastj=j
          j=k
          simultan+=1
        end
        k+=1
      end
      while k<=length(T.relators) && (T.flags[k]>1 || (flag==0 && T.flags[k]==0))
        k+=1
      end
      if k>length(T.relators) j=lastj end
      if i<=j
#       @show T.flags
        altered=SearchC(T,i,j)
        if !isempty(altered)
          debug(P,3,"SearchC($i,$j): ",length(altered)," altered")
          modified=true
          T.flags[altered].=2
#         @show T.flags,altered
        end
        i=j
      end
    end
    for i in eachindex(T.flags) if T.flags[i] in 1:2 T.flags[i]-=1 end end
    if modified
     if sum(length,T.relators)<oldtotal
        T.modified=true
#       print("in search avant handle");display_balanced(P)
        HandleLength1Or2Relators(P)
#       print("in search apres handle");display_balanced(P)
        sort!(P)
        PrintStatus(P,2)
      end
    end
    if !modified break end
  end
end

"""
`SearchEqual( P )`

`SearchEqual`  performs Tietze transformations on a  presentation  `P`.
It tries to alter relators by substituting common subwords of relators by
subwords of equal length.

The  idea is  to find  pairs of  relators `r₁`  length `l₁` and `r₂` length
`l₂`, such that `l₁` is even, `l₁≤l₂`, and `r₁` and `r₂` coincide (possibly
after  inverting or conjugating one of them) in some maximal subword `w` of
length  `≥l₁/2`. Let `l` be the length  of `w`. Then, if `l>l₁/2`, the pair
is  handled  as  in  `Search`.  Otherwise,  if `l=l₁/2`, then `SearchEqual`
substitutes  each copy of `w`  in `r₂` by the  inverse complement of `w` in
`r₁`.

The  Tietze  option   parameter   `P.searchSimultaneous`  is   used  by
`SearchEqual` in the same way as described for `Search`.

However, `SearchEqual` does  not use the  parameter `P.saveLimit`:
The loop over the relators is executed exactly once.
"""
function SearchEqual(P::Presentation)
  debug(P,3,"SearchEqual length subwords")
  Check(P)
  T=P
  simultanlimit=P.searchSimultaneous
  sort!(P)
  rels=T.relators
  modified=false
  oldtotal=sum(length,T.relators)
  i=1
  while i<length(T.relators)
    leng=length(rels[i])
    if leng>3 && iseven(leng)
      simultan=1
      j=i
      lastj=0
      k=i+1
      while k<=length(T.relators) && (length(rels[k])<=leng && simultan<simultanlimit)
        if length(rels[k])==leng
          lastj=j
          j=k
          simultan+=1
        end
        k+=1
      end
      while k<=length(T.relators) && length(rels[k])<leng k+=1 end
      if k>length(T.relators) j=lastj end
      if i<=j
#       display_balanced(P,true)
        altered=SearchC(T, i, j, true)
        if !isempty(altered)
          debug(P,3,"SearchCeq($i,$j) altered=",altered)
          T.flags[altered].=1
          modified=true
        end
        i=j
      end
    end
    i+=1
  end
  if modified
   if sum(length,T.relators)<oldtotal
      T.modified=true
      HandleLength1Or2Relators(P)
      sort!(P)
      PrintStatus(P,2)
    end
  end
end

"""
`Substitute(P::Presentation,n=1,eliminate=0)`

Basically,  it substitutes a squarefree word of length 2 as a new generator
and  then eliminates a generator from  the extended generator list. We will
describe this process in more detail.

`n`  is expected to be a positive integer, and `eliminate` is expected to be
0, 1, or 2.

`Substitute`   first  determines  the  `n`  most  frequently  occurring
squarefree  relator subwords  of  length  2 and sorts  them by decreasing
numbers of occurrences.  Let `ab` be the `n`th word in that list, and let
`i` be the smallest  positive integer  which has not  yet  been used as a
generator number.  Then `Substitute` defines a new  generator `P.i`
(see  `AddGenerator` for  details), adds it to  the presentation together
with a new relator `P.i^{-1}ab`, and  replaces all occurrences of `ab` in
the given relators by `P.i`.

Finally,  it eliminates some generator  from  the  extended presentation.
The  choice  of that  generator  depends  on  the  actual  value  of  the
`eliminate` parameter:

If `eliminate` is zero, then the generator to be eliminated is  chosen as
by the  `Eliminate` command.  This means  that in this case it may well
happen that  it is the generator `P.i` just introduced which  is  now
deleted  again  so  that  you do  not  get  any  remarkable  progress  in
transforming  your  presentation.   On  the other  hand,  this  procedure
guaranties that the total length of the relators will not be increased by
a call of `Substitute` with `eliminate = 0`.

Otherwise,  if `eliminate` is 1 or 2, then  `Substitute` eliminates the
respective factor of the  substituted word `ab`, i.e., `a` for `eliminate
=  1` or  `b` for `eliminate = 2`.  In this case, it may well happen that
the  total  length  of  the  relators  increases,  but  sometimes such an
intermediate  extension  is  the  only  way to  finally  reduce  a  given
presentation.

In order to decide which arguments might be appropriate for the next call
of `Substitute`, often it  is helpful to print out a  list of the  most
frequently occurring squarefree  relator subwords of  length 2.  You  may
use the `show_pairs` command described below to do this.

As an example we handle a subgroup of index 266 in the Janko group `J₁`.

|    gap> F2 := FreeGroup( "a", "b" );;
    gap> J1 := F2 / [ F2.1^2, F2.2^3, (F2.1*F2.2)^7,
         Comm(F2.1,F2.2)^10, Comm(F2.1,F2.2^-1*(F2.1*F2.2)^2)^6 ];;
    gap> a := J1.1;;  b := J1.2;;
    gap> H := Subgroup ( J1, [ a, b^(a*b*(a*b^-1)^2) ] );;
    gap> P := PresentationSubgroup( J1, H );
    << presentation with 23 gens and 82 rels of total length 530 >>|
    P=Presentation("
1: a=A
2: r=R
3: b=B
4: g=G
5: p=P
6: q=Q
7: l=L
8: j=J
9: ied=1
10: ccc=1
11: dca=1
12: iii=1
13: eee=1
14: tao=1
15: fac=1
16: www=1
17: vnl=1
18: mHj=1
19: utf=1
20: wup=1
21: Mkg=1
22: snq=1
23: tSf=1
24: vvv=1
25: gc=mN
26: rg=kS
27: lo=oP
28: gb=QG
29: Hf=fi
30: bh=dO
31: gj=LG
32: lb=Vf
33: gr=RG
34: gCbeb=1
35: weedU=1
36: rhnng=1
37: uuuuu=1
38: kop=hAB
39: gal=feB
40: gsf=sfB
41: agf=LBE
42: ahn=iBG
43: hab=koP
44: wUDbgbf=1
45: wUDbgbf=1
46: jCabeDb=1
47: wuOnfDU=1
48: kHbdiib=1
49: jbdEbac=1
50: gbca=HBAc
51: qlaf=fIUP
52: babc=EcBA
53: gbca=HBAc
54: abac=cBAB
55: KgHm=MhGk
56: rbed=fdFT
57: kMha=RBdF
58: babc=EcBA
59: Fafp=PFAf
60: qlaf=fIUP
61: rlns=hKRh
62: grgs=sLQL
63: rlns=hKRh
64: rbed=fdFT
65: FHbEdFacb=1
66: sfOfIFaHg=1
67: vlgbaccsl=1
68: vbIeebrgs=1
69: bdIFa=kfWuF
70: bdIFa=kfWuF
71: bgbdi=diiAI
72: pUtHm=UtHmG
73: ababab=BABABA
74: lgrkgl=VGRGsn
75: aTFaft=TFAftA
76: lnqNln=nQNLnQ
77: gqgbdp=RGQGBd
78: lgrkgl=VGRGsn
79: jrglql=mGMJRG
80: bhgMjmg=dFLfDBh
81: kMjmKrkMj=RkMJmKRkM
82: jmgMjmgMjm=mGMJmGMJmG
|    gap> GoGo( P );
    &I  there are 3 generators and 47 relators of total length 1368
    &I  there are 2 generators and 46 relators of total length 3773
    &I  there are 2 generators and 46 relators of total length 2570
    gap> GoGo( P );
    &I  there are 2 generators and 46 relators of total length 2568
    gap> GoGo( P );
    gap> & We do not get any more progress without substituting a new
    gap> & generator
    gap> Substitute( P );
    &I  substituting new generator _x28 defined by _x6*_x23^-1
    &I  eliminating _x28 = _x6*_x23^-1
    gap> & GAP cannot substitute a new generator without extending the
    gap> & total length, so we have to explicitly ask for it
    gap> show_pairs( P );
    &I  1.  504  occurrences of  _x6 * _x23^-1
    &I  2.  504  occurrences of  _x6^-1 * _x23
    &I  3.  448  occurrences of  _x6 * _x23
    &I  4.  448  occurrences of  _x6^-1 * _x23^-1
    gap> Substitute( P, 2, 1 );
    &I  substituting new generator _x29 defined by _x6^-1*_x23
    &I  eliminating _x6 = _x23*_x29^-1
    &I  there are 2 generators and 46 relators of total length 2867
    gap> GoGo( P );
    &I  there are 2 generators and 45 relators of total length 2417
    &I  there are 2 generators and 45 relators of total length 2122
    gap> Substitute( P, 1, 2 );
    &I  substituting new generator _x30 defined by _x23*_x29^-1
    &I  eliminating _x29 = _x30^-1*_x23
    &I  there are 2 generators and 45 relators of total length 2192
    gap> GoGo( P );
    &I  there are 2 generators and 42 relators of total length 1637
    &I  there are 2 generators and 40 relators of total length 1286
    &I  there are 2 generators and 36 relators of total length 807
    &I  there are 2 generators and 32 relators of total length 625
    &I  there are 2 generators and 22 relators of total length 369
    &I  there are 2 generators and 18 relators of total length 213
    &I  there are 2 generators and 13 relators of total length 141
    &I  there are 2 generators and 12 relators of total length 121
    &I  there are 2 generators and 10 relators of total length 101
    gap> show_pairs( P );
    &I  1.  19  occurrences of  _x23 * _x30^-1
    &I  2.  19  occurrences of  _x23^-1 * _x30
    &I  3.  14  occurrences of  _x23 * _x30
    &I  4.  14  occurrences of  _x23^-1 * _x30^-1
    gap> & If we save a copy of the current presentation, then later we
    gap> & will be able to restart the computation from the current state
    gap> P1 := Copy( P );;
    gap> & Just for demonstration, let's make an inconvenient choice
    gap> Substitute( P, 3, 1 );
    &I  substituting new generator _x31 defined by _x23*_x30
    &I  eliminating _x23 = _x31*_x30^-1
    &I  there are 2 generators and 10 relators of total length 122
    gap> GoGo( P );
    &I  there are 2 generators and 9 relators of total length 105
    gap> & The presentation is worse than the one we have saved, so let's
    gap> & restart from that one again
    gap> P := Copy( P1 );
    << presentation with 2 gens and 10 rels of total length 101 >>
    gap> Substitute( P, 2, 1);
    &I  substituting new generator _x31 defined by _x23^-1*_x30
    &I  eliminating _x23 = _x30*_x31^-1
    &I  there are 2 generators and 10 relators of total length 107
    gap> GoGo( P );
    &I  there are 2 generators and 9 relators of total length 84
    &I  there are 2 generators and 8 relators of total length 75
    gap> Substitute( P, 2, 1);
    &I  substituting new generator _x32 defined by _x30^-1*_x31
    &I  eliminating _x30 = _x31*_x32^-1
    &I  there are 2 generators and 8 relators of total length 71
    gap> GoGo( P );
    &I  there are 2 generators and 7 relators of total length 56
    &I  there are 2 generators and 5 relators of total length 36
    gap> relators(P);
    &I  1. _x32^5
    &I  2. _x31^5
    &I  3. _x31^-1*_x32^-1*_x31^-1*_x32^-1*_x31^-1*_x32^-1
    &I  4. _x31*_x32*_x31^-1*_x32*_x31^-1*_x32*_x31*_x32^-2
    &I  5. _x31^-1*_x32^2*_x31*_x32^-1*_x31^2*_x32^-1*_x31*_x32^2 |

As shown in the preceding example, you can use the `Copy` command to save
a copy of a presentation and to restart from it  again if you want
to try an alternative strategy.  
"""
function Substitute(P::Presentation,n::Int=1,elim::Int=0)
  if !(0<=elim<=2) error("last argument must be 0, 1, or 2") end
  T=P
  if length(T.generators)<2 return end
  pairs=MostFrequentPairs(P, n)
  if length(pairs)<n
    debug(P,1,"Substitute nbpairs $n is out of range")
    return
  end
  m,i,j,k=pairs[n]
  if k>1 i=T[-i] end
  if isodd(k) j=T[-j] end
  pair=[i, j]
  gen=NewGenerator(P)
  word=AbsWord(pair,T.generators)
  debug(P,1,"Substitute new generator ",gen," defined by ", word)
  AddRelator(P,inv(gen)*word)
  UpdateGeneratorImages(P, pair)
  printlevel=P.debug;P.debug=0;Search(P);
  P.debug=2
  if elim>0 EliminateGen(P, abs(pair[elim])) else EliminateGen1(P) end
  P.debug=printlevel
  if T.modified Search(P) end
  if T.numredunds>0 RemoveGenerators(P) end
  PrintStatus(P,1)
end

"""
`find_exponents(presentation)`  tries  to  find  exponents  for  the Tietze
generators  and  returns them  in  a  list  parallel  to  the  list  of the
generators.
"""
function find_exponents(T::Presentation)
  debug(T,3,"find generator exponents")
  rels=T.relators
  exponents=fill(0,length(T.generators))
  for i in eachindex(rels)
    if isempty(rels[i]) || !allequal(rels[i]) continue end
    if rels[i][1]<0 rels[i]=-rels[i] end
    num=rels[i][1]
    exponents[num]=exp=gcd(exponents[num], length(rels[i]))
    if exp<length(rels[i])
      rels[i]=fill(num,exp)
      T.flags[i]=1
      T.modified=true
    end
  end
  exponents
end

"""
`SubstituteCyclicJoins(presentation)`

`SubstituteCyclicJoins`    performs   Tietze   transformations   on   a
presentation `P`.  It tries to find pairs of generators `a` and `b`, say,
for which among  the  relators (possibly after  inverting  or conjugating
some of them) there are the commutator `[a,b]` and powers `a^m` and `b^n`
with mutually  prime exponents  `m` and  `n`.   For each  such  pair,  it
substitutes the product `ab`  as  a new generator, and then it eliminates
the generators `a` and `b`.
"""
function SubstituteCyclicJoins(P::Presentation)
  debug(P,3,"SubstituteCyclicJoins")
  exponents=find_exponents(P)
  P.modified=false
  if iszero(exponents) return end
  sort!(P)
  gens=P.generators
  rels=P.relators
  i=1
  while i<=length(P.relators) && length(P.relators[i])<=4
    rel=rels[i]
    if length(P.relators[i])==4 && rel[1]==P[-rel[3]] && rel[2]==P[-rel[4]]
      num1=abs(rel[1])
      exp1=exponents[num1]
      num2=abs(rel[2])
      exp2=exponents[num2]
      if exp1>0 && exp2>0 && gcd(exp1, exp2)==1
        gen=NewGenerator(P)
        gen2=gens[num2]
        debug(P,1,"I substituting new generator ",gen," defined by ",
                gens[num1]*gen2)
        AddRelator(P, gens[num1]/gen^exp2)
        AddRelator(P, gen2/gen^exp1)
        UpdateGeneratorImages(P,[num1, num2])
        debug=P.debug;P.debug=2;
        Eliminate(P,gens[num1]);Eliminate(P,gen2);
        P.debug=debug
        P.modified=true
        i=0
      end
    end
    i+=1
  end
  if P.modified
    HandleLength1Or2Relators(P)
    sort!(P)
    PrintStatus(P,1)
  end
end

"""
`Substitute(P::Presentation,word[,string])`

`word`  is be either an abstract word or a Tietze word in the generators of
`P`.  It substitutes the given word as a new generator of `P`. This is done
as follows.

First,  `Substitute` creates a new abstract generator, `g` say, and adds it
to the presentation `P`, then it adds a new relator `g^{-1}⋅word` to `P`.

If  a  string  `string`  has  been  specified  as  third  argument, the new
generator  `g` will be named  by `string`, otherwise it  will get a default
name  `_xi` as  described with  the function  `AddGenerator` (see "Changing
Presentations").

More precisely: If, for instance, |word| is an abstract word, a call

|    Substitute( P, word );|

is more or less equivalent to

|    AddGenerator( P );
    g := P.generators[Length( P.generators )];
    AddRelator( P, g^-1 * word );|

whereas a call
|    Substitute( P, word, string );|
is more or less equivalent to

|    g := AbstractGenerator( string );
    AddGenerator( P, g );
    AddRelator( P, g^-1 * word );|

The   essential  difference   is,   that  `Substitute`,  as  a   Tietze
transformation  of `P`, saves and updates  the lists  of generator images
and preimages if  they are being  traced under the Tietze transformations
applied to `P` (see  the function `tracing` below), whereas
a call  of the  function `AddGenerator` (which   does not perform  Tietze
transformations) will delete these lists and hence terminate the tracing.

Example.

|    gap> G := PerfectGroup( 960, 1 );
    PerfectGroup(960,1)
    gap> P := Presentation( G );
    << presentation with 6 gens and 21 rels of total length 84 >>
```julia-rep1
P=Presentation("
1: a=A
2: e=E
3: f=F
4: c=C
5: d=D
6: bbb=1
7: fc=CF
8: dc=CD
9: fe=EF
10: ec=CE
11: da=AF
12: ed=DE
13: fd=DF
14: ca=AE
15: ea=AC
16: fa=AD
17: fb=bE
18: Bcbfd=1
19: Bebfe=1
20: Bdbfedc=1
21: babab=ABABA
")
julia> P.generators
6-element Vector{AbsWord}:
 a
 b
 c
 d
 e
 f
julia> Presentations.GoGo(P)
Presentation: 3 generators, 10 relators, total length 81
Presentation: 3 generators, 10 relators, total length 80
Presentation: 3 generators, 10 relators, total length 80

julia> showgens(P)
1. a 31 occurrences involution
2. b 26 occurrences
3. d 23 occurrences involution

julia> a,b=P.generators[1:2]
2-element Vector{AbsWord}:
 a
 b

julia> Presentations.Substitute( P, a*b, :ab )
#started tracing generator images
#Substitute new generator ab defined by ab
#started tracing generator images
Presentation: 4 generators, 11 relators, total length 83

julia> Presentations.Go(P)
#Presentation: 3 generators, 10 relators, total length 74
#Presentation: 3 generators, 10 relators, total length 74

julia> Presentations.showgens(P)
1. a 25 occurrences involution
2. d 23 occurrences involution
3. ab 26 occurrences
```
"""
function Substitute(P::Presentation,word::Union{Vector{Int},AbsWord},arg...)
  Check(P)
  gens=P.generators
  if word isa Vector
    tzword=reduceword(word)
    word=AbsWord(tzword, gens)
  else
    tzword=TietzeWord(word, gens)
  end
  images=0
  if haskey(P,:imagesOldGens)
    images=P.imagesOldGens
    delete!(P,:imagesOldGens)
  end
  if length(arg)==1
    gen=AbsWord(arg[1])
    AddGenerator(P, gen)
  else
    AddGenerator(P)
    gen=P.generators[end]
  end
  debug(P,1,"Substitute new generator ", gen," defined by ", word)
  AddRelator(P,inv(gen)*word)
  if images isa Vector
    P.imagesOldGens=images
    UpdateGeneratorImages(P, tzword)
  end
  PrintStatus(P,1)
end

"""
`tracing(P::Presentation,stop=false)`

a  sequence of  Tietze transformations  applied to  a presentation `P₁` and
ending  up with a presentation `P₂`, defines an isomorphism `φ` between the
groups defined by `P₁` and `P₂`, respectively. To know `φ` (resp. `φ⁻¹`) we
need  to  know  `φ(old  generators)`  or  `φ⁻¹(new generators)`. The Tietze
transformations  functions  can  trace  these  images.  This is not done by
default  since the involved words may grow to tremendous length; it will be
done  on request  by calling  `tracing(P)`. A  call `tracing(P,false)` will
stop tracing.

`tracing` initializes three fields of `P`:

`P.oldGenerators`:  is  initialized  to  a  copy of `P.generators`.

`P.imagesOldGens`:  is the list `φ(old generators)`  as Tietze words in the
new generators. Its `i`-th entry is initialized to the Tietze word |[i]|.

`P.preImagesNewGens`: The list `φ⁻¹(new generators)` as Tietze words in the
old generators. Its `i`-th entry is initialized by the Tietze word |[i]|.

The  existence of  these fields  will cause  the Tietze  transformations to
update the lists `imagesOldGens` and `preImagesNewGens`.

There are a few restrictions concerning the tracing of generator images:

The  functions `AddGenerator`, `AddRelator`, and `RemoveRelator` may change
the  isomorphism type  of the  presentation. Therefore,  if any  of them is
called, it will call `tracing(P,false)`.

You  can reinitialize  tracing at  any later  state by  calling `tracing()`
again: if the above fields do already exist when `tracing` is being called,
they will be initialized again.
"""
function tracing(P::Presentation,stop=false)
  if stop 
    if !haskey(P, :imagesOldGens) return end
    delete!(P, :imagesOldGens)
    delete!(P, :preImagesNewGens)
    delete!(P, :oldGenerators)
    debug(P,1,"terminated tracing generator images")
  else
    P.oldGenerators=deepcopy(P.generators)
    P.imagesOldGens=map(i->[i],eachindex(P.generators))
    P.preImagesNewGens=map(i->[i],eachindex(P.generators))
    debug(P,1,"started tracing generator images")
  end
end

"""
`UpdateGeneratorImages(T::Presentation, word )`

assumes  that a new  generator defined by  the Tietze word  `word` has just
been  added to the presentation.  It converts `word` from  a Tietze word in
the  new generators to  a Tietze word  in the old  generators and adds that
word to the list of preimages.
"""
function UpdateGeneratorImages(T::Presentation, word)
  if !haskey(T,:preImagesNewGens) return end
  preim=T.preImagesNewGens
  newim=Int[]
  for num in word append!(newim, num>0 ? preim[num] : -reverse(preim[-num])) end
  push!(preim, reduceword(newim))
end

"""
`UpdateGeneratorImages(T::Presentation, n::Integer, word )`

assumes  that  the  `n`-th  generator  has  just  been replaced by the word
`word`.  It  updates  the  images  of  the old generators by replacing each
occurrence of the `n`-th generator by the given Tietze word `word`. This is
all if `word` does not involve `n` itself. If `word` involves `n` once then
the  preimages are updated  in turn. If  `word` involves `n`  more than one
time then the preimages are no more maintained.
"""
function UpdateGeneratorImages(T::Presentation, n::Integer, word)
  if n<0 error("i did not know it was possible") end
  invert(w)=-reverse(w)
  if !haskey(T,:imagesOldGens) return end
  invword=invert(word)
  for (i,image) in enumerate(T.imagesOldGens)
    newim=Int[]
    for (j,x) in enumerate(image)
      if x==n append!(newim, word)
      elseif x==-n append!(newim, invword)
      else push!(newim, x)
      end
    end
    T.imagesOldGens[i]=reduceword(newim)
  end
  if !haskey(T,:preImagesNewGens) return end
  p=findall(x->x==n || x==-n, word)
  if length(p)!=1 
   if length(p)!=0 delete!(T,preImagesNewGens) end
    return
  end
  p=only(p)
  n1=word[p]
  word=vcat(invert(word[1:p-1]),n,invert(word[p+1:end]))
  word=vcat(map(x->x>0 ? T.preImagesNewGens[x] : invert(T.preImagesNewGens[-x]),
                word)...)
  T.preImagesNewGens[abs(n1)]=reduceword(n1>0 ? word : invert(word))
# images(T)
end

"""
`images(P::Presentation)`

If  `P` is a presentation in which generator images and preimages are being
traced  through all Tietze transformations  applied to `P`, `images` prints
the  preimages of the current generators as words in the old generators and
the images of the old generators as words in the current generators.

```julia-repl
julia> P=Presentation("
1: a=A
2: e=E
3: f=F
4: c=C
5: d=D
6: bbb=1
7: fc=CF
8: dc=CD
9: fe=EF
10: ec=CE
11: da=AF
12: ed=DE
13: fd=DF
14: ca=AE
15: ea=AC
16: fa=AD
17: fb=bE
18: Bcbfd=1
19: Bebfe=1
20: Bdbfedc=1
21: babab=ABABA
") 
Presentation: 6 generators, 21 relators, total length 84
julia> tracing(P)
julia> Presentations.Go(P)
Presentation: 3 generators, 10 relators, total length 81

julia> P.imagesOldGens
6-element Vector{Vector{Int64}}:
 [1]
 [2]
 [1, -2, 1, 3, 1, 2, 1]
 [3]
 [-2, 1, 3, 1, 2]
 [1, 3, 1]

julia> P.preImagesNewGens
3-element Vector{Vector{Int64}}:
 [1]
 [2]
 [4]
```

```julia-rep1
julia> Presentations.images(P)
current generators in terms of the old ones:
  a=a
  b=b
  d=d
old generators in terms of the current ones:
  a=a
  b=b
  c=ab⁻¹adaba
  d=d
  e=b⁻¹adab
  f=ada
```
"""
function images(P::Presentation)
  Check(P)
  if !haskey(P, :imagesOldGens)
    println("#I presentation is not tracing generator images"); return
  end
  if haskey(P,:preImagesNewGens)
    println("current generators in terms of the old ones:")
    for i in eachindex(P.generators) xprintln("  ",P.generators[i],"=",
           AbsWord(P.preImagesNewGens[i],P.oldGenerators)) end
  end
  println("old generators in terms of the current ones:")
  for i in eachindex(P.oldGenerators) xprintln("  ",P.oldGenerators[i],"=",
         AbsWord(P.imagesOldGens[i],P.generators))
  end
end

"""
RenumberGenerators( presentation, sort list )

renumbers the generators of the given presentation
according to the given sort list.
"""
function RenumberGenerators(P::Presentation, L::AbstractVector)
  debug(P,3,"RenumberGenerators")
  RemoveGenerators(P)
  numgens=length(P.generators)
  if sort(L)!=1:numgens
    error("<L> must determine a permutation of the generator numbers")
  end
  invsort=invperm(L)
  gens=P.generators
  rels=P.relators
  oldgens=deepcopy(gens)
  oldinvs=copy(P.inverses)
  new=copy(P.inverses)
  println("gens==",gens,"\ninvs==",P.inverses,"\nrels==",rels)
  for i in 1:numgens
    j=L[i]
    gens[i]=oldgens[j]
    if oldinvs[numgens+1+j]==j P[-i]=i else P[-i]=-i end
    new[numgens+1+i]=invsort[i]
    new[numgens+1-i]=-invsort[i]
  end
  println("gens==",gens,"\ninvs==",P.inverses,"\nnew==",new)
  for rel in P.relators
    for i in eachindex(rel) rel[i]=new[numgens+1+rel[i]] end
  end
  println("rels==", rels)
  if haskey(P, :imagesOldGens)
    images(P)
    for i in 1:length(P.imagesOldGens)
      P.imagesOldGens[i]=map(j->new[numgens+1+j],P.imagesOldGens[i])
    end
    permute!(P.preImagesNewGens, L)
    println("renumbered to")
    images(P)
  end
end

"""
`simplify(p [,tries])`

simplify  the  presentation  `p`.  We  have  found heuristics which make it
somewhat  efficient, but the algorithm depends  on random numbers so is not
reproducible.  The main  idea is  to rotate  relators between  calls to the
basic  `Presentations.Go` function. By default 100 such rotations are tried
(unless  the  presentation  is  so  small  that  less rotations exhaust all
possible  ones), but the actual number tried  can be controlled by giving a
second  parameter `tries` to the function. Another useful tool to deal with
presentations is `tryconjugate`.

```julia-rep1
julia> display_balanced(p)
1: ab=ba
2: dbd=bdb
3: bcb=cbc
4: cac=aca
5: adca=cadc
6: dcdc=cdcd
7: adad=dada
8: Dbdcbd=cDbdcb
9: adcDad=dcDadc
10: dcdadc=adcdad
11: dcabdcbda=adbcbadcb
12: caCbdcbad=bdcbadBcb
13: cbDadcbad=bDadcbadc
14: cdAbCadBc=bdcAbCdBa
15: cdCbdcabdc=bdcbadcdaD
16: DDBcccbdcAb=cAbCdcBCddc
17: CdaBdbAdcbCad=abdcAbDadBCbb
18: bdbcabdcAADAdBDa=cbadcbDadcBDABDb
19: CbdbadcDbbdCbDDadcBCDAdBCDbdaDCDbdcbadcBCDAdBCDBBdacDbdccb=abdbcabdcAdcbCDDBCDABDABDbbdcbDadcbCDAdBCabDACbdBadcaDbAdd

julia> simplify(p)
Presentation: 4 generators, 18 relators, total length 304
Presentation: 4 generators, 18 relators, total length 284
Presentation: 4 generators, 17 relators, total length 264
Presentation: 4 generators, 16 relators, total length 256
Presentation: 4 generators, 15 relators, total length 244
Presentation: 4 generators, 15 relators, total length 240
Presentation: 4 generators, 15 relators, total length 226
Presentation: 4 generators, 15 relators, total length 196
Presentation: 4 generators, 15 relators, total length 178
Presentation: 4 generators, 15 relators, total length 172
Presentation: 4 generators, 14 relators, total length 158

julia> display_balanced(p)
1: ab=ba
2: dbd=bdb
3: bcb=cbc
4: cac=aca
5: adAc=cadA
6: dcdc=cdcd
7: adad=dada
8: CdBcbd=bCdBcb
9: adcDad=dcDadc
10: dcdadc=adcdad
11: cbdcbdc=dcbdcbd
12: dcbadcbda=adcbcadcb
13: cbCDadcab=DadcbadcD
14: caDCbdBcADbda=bDBaDbADcbadc
```
"""
function simplify(P::Presentation,lim=100)
  if isempty(P.relators) return end
  rot(tt,i)=tt[i]=circshift(tt[i],-1)
  function test()
    if prod(filter(!iszero,big.(length.(P.relators))))<lim
      tt=P.relators
      v=fill(0,length(tt))
      if length(v)==0 return false end
      while true
        before=[length(tt),sum(length,tt)]
        j=length(v)
        while v[j]==length(tt[j])-1
          rot(tt,j)
          v[j]=0
          j-=1
          if j==0 return false end
        end
        rot(tt,j)
        v[j]+=1
        GoGo(P)
        tt=P.relators
        if [length(tt),sum(length,tt)]<before return true end
      end
    else
      for k in 1:lim
        tt=P.relators
        before=[length(tt),sum(length,tt)]
        rot(tt,rand(eachindex(tt)))
        GoGo(P)
        tt=P.relators
        if [length(tt),sum(length,tt;init=0)]<before return true end
      end
    end
    return false
  end
  while test() end
end

"""
`conjugate(p,conjugation)`

This program modifies a presentation by conjugating a generator by another.
The  conjugation to  apply is  described by  a length-3  string of the same
style  as  the  result  of  `display_balanced`,  that is `"abA"` means
replace  the second generator by its  conjugate by the first, and  `"Aba"`
means replace it by its conjugate by the inverse of the first.

```julia-rep1
julia> display_balanced(P)
1: dabcd=abcda
2: dabcdb=cabcda
3: bcdabcd=dabcdbc

julia> display_balanced(conjugate(P,"Cdc"))
<< presentation with 4 generators, 3 relators of total length 36>>
1: dcabdc=cabdca
2: abdcab=cabdca
3: bdcabd=cabdca
```
"""
function conjugate(p::Presentation, s)
  if !(s isa Vector{<:Integer})
    l=string.(mon.(p.generators))
    if !all(l)do ss
      if length(ss)!=1 return false end
      if !('a'<=ss[1]<='z') return false end
      true
    end
    error("not all gens in a-z") end
    minmaj=vcat(map(x->x[1],l),map(x->uppercase(x[1]),l))
    s=map(c->findfirst(==(c),minmaj),collect(s))
    s=map(x->x>length(l) ? length(l)-x : x,s)
  end
  if length(s)!=3 || s[1]!=-s[3] error("should be conjugate like abA") end
  p=deepcopy(p);SubstituteGen(p,s[2],s)
  simplify(p, 100)
  p
end

"""
`tryconjugate(p[,goal[,printlevel]])`

This program tries to simplify group presentations by applying conjugations
to  the  generators.  The  algorithm  depends  on  random  numbers,  and on
tree-searching,  so is  not reproducible.  By default  the program stops as
soon  as a shorter presentation is found.  Sometimes this does not give the
desired  presentation.  One  can  give  a  second argument `goal`, then the
program  will only stop when  a presentation of length  less than `goal` is
found.  Finally, a third  argument can be  given and then all presentations
the  programs runs  over which  are of  length less  than or  equal to this
argument are displayed. Due to the non-deterministic nature of the program,
it  may be useful to  run it several times  on the same input. Upon failure
(to improve the presentation), the program returns `p`.

```julia-rep1
julia> display_balanced(p)
1: ba=ab
2: dbd=bdb
3: cac=aca
4: bcb=cbc
5: dAca=Acad
6: dcdc=cdcd
7: adad=dada
8: dcDbdc=bdcbdB
9: dcdadc=adcdad
10: adcDad=dcDadc
11: BcccbdcAb=dcbACdddc
julia> p=tryconjugate(p)
Presentation: 4 generators, 11 relators, total length 100
dcD=> Presentation: 4 generators, 10 relators, total length 90
# dcD gives Presentation: 4 generators, 10 relators, total length 90
Presentation: 4 generators, 10 relators, total length 90

julia> p=tryconjugate(p)
Dcd=> Presentation: 4 generators, 10 relators, total length 88
# Dcd gives Presentation: 4 generators, 10 relators, total length 88
Presentation: 4 generators, 10 relators, total length 88

julia> p=tryconjugate(p)
dcD=> Presentation: 4 generators, 10 relators, total length 90
Dbd=> Presentation: 4 generators, 10 relators, total length 96
Aca=> Presentation: 4 generators, 9 relators, total length 84
Presentation: 4 generators, 8 relators, total length 76
# Aca gives Presentation: 4 generators, 8 relators, total length 76
Presentation: 4 generators, 8 relators, total length 76

julia> p=tryconjugate(p)
Bcb=> Presentation: 4 generators, 8 relators, total length 70
# Bcb gives Presentation: 4 generators, 8 relators, total length 70
Presentation: 4 generators, 8 relators, total length 70

julia> p=tryconjugate(p)
Cac=> Presentation: 4 generators, 8 relators, total length 64
# Cac gives Presentation: 4 generators, 8 relators, total length 64
Presentation: 4 generators, 8 relators, total length 64

julia> p=tryconjugate(p)
caC=> Presentation: 4 generators, 8 relators, total length 58
# caC gives Presentation: 4 generators, 8 relators, total length 58
Presentation: 4 generators, 8 relators, total length 58

julia> p=tryconjugate(p)
Cac=> Presentation: 4 generators, 8 relators, total length 64
Cbc=> Presentation: 4 generators, 7 relators, total length 50
# Cbc gives Presentation: 4 generators, 7 relators, total length 50
Presentation: 4 generators, 7 relators, total length 50

julia> p=tryconjugate(p)
cdC=> Presentation: 4 generators, 7 relators, total length 56
Dcd=> Presentation: 4 generators, 7 relators, total length 54
Cac=> Presentation: 4 generators, 7 relators, total length 48
# Cac gives Presentation: 4 generators, 7 relators, total length 48
Presentation: 4 generators, 7 relators, total length 48

julia> p=tryconjugate(p)
caC=> Presentation: 4 generators, 7 relators, total length 50
Cdc=> Presentation: 4 generators, 7 relators, total length 50
Dbd=> Dcd=> Presentation: 4 generators, 7 relators, total length 60
Bab=> Aba=> Aca=> Presentation: 4 generators, 7 relators, total length 46
# Aca gives Presentation: 4 generators, 7 relators, total length 46
Presentation: 4 generators, 7 relators, total length 46

julia> display_balanced(p)
1: db=bd
2: ba=ab
3: cac=aca
4: ada=dad
5: bcb=cbc
6: cdcd=dcdc
7: AdCacd=cAdCac
```
"""
function tryconjugate(p::Presentation,tp=[0,0];info=[0,0])
  perf(p)=[length(p.relators), sum(length,p.relators)]
  function triples(p)
    res=Vector{Int}[]
    for v in p.relators, i in 1:length(v)-2
      if v[i]==-v[i+2]
        if v[i+1]>0 push!(res,v[i:i+2])
        else        push!(res,-reverse(v[i:i+2]))
        end
      end
    end
    first.(sort(tally(res),by=x->-x[2]))
  end
  function expand(p)
    res=map(x->[x, conjugate(p, x)], triples(p))
    sort!(res,by=x->perf(x[2]))
  end
  p=deepcopy(p)
  GoGo(p)
  if iszero(tp) tp=perf(p) end
  p1=[]
  for c in triples(p)
    print(alphab(c),"=> ")
    n=deepcopy(p);SubstituteGen(n,c[2],[-c[1],c[2],-c[3]]);GoGo(n)
    if perf(n)<=info
      println("# ",alphab(c)," gives rels/len ", perf(n))
      display_balanced(n)
    end
    if perf(n)<tp
      println("# ",alphab(c)," gives ",n)
      return n
    end
    push!(p1,(c,n))
  end
  sort!(p1,by=x->perf(x[2]))
  for p2 in p1, c in triples(p2[2])
    n=deepcopy(p2[2]);SubstituteGen(n,c[2],[-c[1],c[2],-c[3]]);GoGo(n)
    print(alphab(p2[1]),"+",alphab(c),"=>")
    if perf(n)<=info
      println("# ",alphab(p2[1]),"->",alphab(c)," gives rels/len ", perf(n))
      display_balanced(n)
    end
    if perf(n)<tp
      println("# ",alphab(p2[1]),"->",alphab(c)," gives ",n)
      return n
    end
  end
  println("\n# could not shrink ",p)
  p
end

using ..Gapjm:gap, Gapjm

function Gapjm.gap(p::Presentation)
  t=p
  s="F:=FreeGroup($(length(t.generators)));\n"
  s*="F.relators:=["
  s*=join(map(t.relators)do r
      join(map(r)do i
        i>0 ? string("F.",i) : string("F.",-i,"^-1")
      end,"*")
          end,",")*"];\n"
  s*"PresentationFpGroup(F);\n"
end

p1=Presentation(
[[1, 1], [7, 7], [3, 3], [6, 6], [8, 8], [4, 4], [6, 2, 3], [4, -2, 1], [7, 4,
 2], [3, 5, 2], [5, 1, 2], [3, 1, -2], [-2, -2, -2], [3, -2, 4], [1, 3, 2],
 [1, 6, -2], [2, 2, 2], [2, 1, 3], [3, 2, 1], [1, -2, 3], [8, 1, -2], [4, 1,
 2], [-2, 5, 3], [-2, 1, 4], [-2, 3, 1], [1, 7, 2], [1, 2, 1, 2], [1, 2, 3,
 -2], [1, 3, -5, 3], [1, 3, -2, -2], [3, -2, 3, -2], [1, 3, 1, 3, 1, 3]])
p2=Presentation(
[[1, 1], [5, 5], [6, 6], [3, 3], [4, 4], [2, 2, 2], [6, 3, 6, 3], [4, 3, 4,
 3], [6, 5, 6, 5], [5, 3, 5, 3], [4, 1, 6, 1], [5, 4, 5, 4], [6, 4, 6, 4], [3,
 1, 5, 1], [5, 1, 3, 1], [6, 1, 4, 1], [6, 2, 5, -2], [-2, 3, 2, 6, 4], [-2,
 5, 2, 6, 5], [-2, 4, 2, 6, 5, 4, 3], [2, 1, 2, 1, 2, 1, 2, 1, 2, 1]])
p3=Presentation(
[[1, 1], [18, 18], [2, 2], [7, 7], [16, 16], [17, 17], [12, 12], [10, 10], [9,
 5, 4], [3, 3, 3], [4, 3, 1], [9, 9, 9], [5, 5, 5], [20, 1, 15], [6, 1, 3],
 [23, 23, 23], [22, 14, 12], [13, -8, 10], [21, 20, 6], [23, 21, 16], [-13,
 11, 7], [19, 14, 17], [20, -19, 6], [22, 22, 22], [7, 3, 14, -13], [18, 7,
 19, -11], [12, 15, 16, -15], [7, 2, 7, 17], [-8, 6, -9, -6], [2, 8, 15, -4],
 [7, 10, 7, 12], [12, 2, -6, 22], [7, 18, 7, 18], [7, -3, 2, 5, 2], [23, 5, 5,
 4, -21], [18, 8, 14, 14, 7], [21, 21, 21, 21, 21], [11, 15, 16, 2, 1, -8],
 [7, 1, 12, 2, -5, -6], [7, 19, 6, 2, -6, -19], [1, 7, 6, 5, 2, 12], [1, 8,
 14, 7, 2, -9], [8, 1, 2, 16, -15, -11], [23, -21, -4, 2, 7, 2, 6], [10, -3,
 1, 2, 5, -4, 2], [23, 21, -15, 14, 6, -4, -21], [11, -8, 2, 4, 9, 9, 2], [10,
 2, 4, -5, 2, 1, 3], [7, 2, 3, 1, -3, 1, 2, 8], [17, 12, 1, 6, 16, 21, 9, -6],
 [2, 1, 2, 3, 1, 2, -3, 5], [1, 2, 1, 3, 2, 1, 2, -3], [-11, 7, -8, 13, -11,
 7, -8, 13], [18, 2, 5, 4, 20, 6, -4, -6], [11, -13, 8, 1, 6, -4, 2, 18], [-6,
 1, 6, 16, -6, 1, 6, 16], [18, 12, 14, 19, -8, 18, 11, -8], [7, 18, 7, 19, 12,
 17, 12, -19], [-6, -8, 2, -5, 4, -6, 1, 3, 2], [19, 6, -15, 6, -9, -6, 1, -8,
 7], [22, 12, 7, 2, 1, 3, 3, 19, 12], [22, 2, -9, 5, 5, 2, 18, 7, 19], [2, 4,
 -9, -6, 1, 6, -21, 23, -6, -11], [2, 7, 2, 4, 9, 9, 1, -9, -9, -4], [16, -21,
 20, -8, 13, 7, -13, 8, -20, 21], [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2], [12,
 7, 18, 11, 7, 12, -14, -19, 7, 18, 7, 22], [1, -20, -6, 1, 6, 20, 1, -20, -6,
 1, 6, 20], [12, 14, 17, -14, 12, 14, 17, -14, 12, 14, 17, -14], [7, 17, 7, 2,
 4, 16, -4, 2, 7, 17, 7, 18], [10, 18, 7, 12, 17, 12, 7, 18, 10, 13, 7, -13],
 [2, 8, 7, -13, 10, 13, 7, -8, 2, 4, -6, 12, 6, -4], [11, -13, 10, 13, -11,
 18, 11, -13, 10, 13, -11, 18, 11, -13, 10, 13, -11, 18], [10, 13, 7, -13, 10,
 13, 7, -13, 10, 13, 7, -13, 10, 13, 7, -13, 10, 13, 7, -13]])
p3simp=Presentation(
[[2, 2, 2], [2, -1, 2, -1, 2, -1], [1, 1, -2, -1, -2, 1, 1, -2, -1, -2], [1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 1, -2, -1, -1, -2, 1, 1, 2, -1, -1, -1,
 2], [1, 1, -2, -1, -1, -2, -1, 2, 1, -2, -1, 2, 1]])
p4=Presentation(
[[2, 3, -2, -3], [5, 3, 5, -3, -5, -3], [3, 4, 3, -4, -3, -4], [4, 2, 4, -2,
 -4, -2], [2, 5, 4, 2, -4, -5, -2, -4], [5, 4, 5, 4, -5, -4, -5, -4], [2, 5,
 2, 5, -2, -5, -2, -5], [5, 4, 5, 2, 5, 4, 1, -5, -2, -5, -4, -5, -2], [-5, 3,
 5, 4, 3, 5, 1, 1, -3, -4, -5, -3, 5, -4], [2, 5, 4, -5, 2, 5, 1, 1, -4, -5,
 -2, 5, -4, -5], [5, 4, 2, 3, 5, 4, 3, 5, 2, -3, -4, -5, -2, -3, -4, -3, -5,
 -2], [4, 2, -4, 3, 5, 4, 3, 2, 5, -3, -4, 3, -5, -2, -3, -4, -5, -3], [4, 3,
 -5, 2, 5, 4, 3, 2, 5, -4, -5, -2, -3, -4, -5, -2, 5, -3], [4, 5, -2, 3, -4,
 2, 5, -3, 4, -2, 3, -5, 4, -3, 2, -4, -5, -3], [4, 5, -4, 3, 5, 4, 2, 3, 5,
 4, 5, -2, -5, -4, -5, -2, -3, -4, -5, -3], [-5, -5, -3, 4, 4, 4, 3, 5, 4, -2,
 3, -4, -5, -5, 4, 3, -4, -5, 4, -3, 2, -4], [-4, 5, 2, -3, 5, 3, -2, 5, 4, 3,
 -4, 2, 5, -3, -3, 4, 3, -5, -2, 5, -3, 2, -4, -5, -3, -2], [3, 5, 3, 4, 2, 3,
 5, 4, -2, -2, -5, -2, 5, -3, -5, 2, -3, 5, 3, 2, 5, 3, -4, -5, -2, 5, -3, -4,
 -5, -2, -3, -4], [-4, 3, 5, 3, 2, 5, 4, -5, 3, 3, 5, -4, 3, -5, -5, 2, 5, 4,
 -3, -4, -5, -2, 5, -3, -4, -5, 3, 5, 2, -5, -4, -5, 3, 5, 4, 3, 2, 5, 4, -3,
 -4, -5, -2, 5, -3, -4, -5, -3, -3, 5, 2, 4, -5, 3, 5, 4, 4, 3, -5, -5, 2, -3,
 5, -2, -4, -5, -2, 3, -5, -3, 4, 2, 5, -3, -2, 4, 3, -5, 2, 5, 4, -3, -4, -5,
 -2, 5, -3, -4, -5, -3, -3, 5, 3, 2, 5, 3, 2, 5, 4, 3, 5, 5, 4, -3, -4, -5, 2,
 -4, -5, -3, -2, -4, -3, -5, -3, -2]])
p5=Presentation(
[[3, -1], [2, -3], [2, 3, -1, -3], [1, 4, -1, -4], [4, 1, -4, -1], [1, 3, -2,
 -3], [1, 2, -3, -2], [4, 3, -4, -3], [3, 4, -3, -4], [2, 3, 4, 1], [3, 2, 3,
 -2, -3, -2], [-1, 4, 1, 2, -1, -4, 1, 4, -2, -4], [-1, 4, 1, 3, -1, -4, 1, 4,
 -3, -4]])
p6=Presentation(
[[1, 2, -1, -2], [4, 2, 4, -2, -4, -2], [2, 3, 2, -3, -2, -3], [3, 1, 3, -1,
 -3, -1], [1, 4, 3, 1, -3, -4, -1, -3], [4, 3, 4, 3, -4, -3, -4, -3], [1, 4,
 1, 4, -1, -4, -1, -4], [-4, 2, 4, 3, 2, 4, -2, -3, -4, -2, 4, -3], [1, 4, 3,
 -4, 1, 4, -3, -4, -1, 4, -3, -4], [4, 3, 4, 1, 4, 3, -4, -1, -4, -3, -4, -1],
 [4, 3, 1, 2, 4, 3, 2, 4, 1, -2, -3, -4, -1, -2, -3, -2, -4, -1], [3, 1, -3,
 2, 4, 3, 2, 1, 4, -2, -3, 2, -4, -1, -2, -3, -4, -2], [3, 2, -4, 1, 4, 3, 2,
 1, 4, -3, -4, -1, -2, -3, -4, -1, 4, -2], [3, 4, -1, 2, -3, 1, 4, -2, 3, -1,
 2, -4, 3, -2, 1, -3, -4, -2], [3, 4, -3, 2, 4, 3, 1, 2, 4, 3, 4, -1, -4, -3,
 -4, -1, -2, -3, -4, -2], [-4, -4, -2, 3, 3, 3, 2, 4, 3, -1, 2, -3, -4, -4, 3,
 2, -3, -4, 3, -2, 1, -3], [-3, 4, 1, -2, 4, 2, -1, 4, 3, 2, -3, 1, 4, -2, -2,
 3, 2, -4, -1, 4, -2, 1, -3, -4, -2, -1], [2, 4, 2, 3, 1, 2, 4, 3, -1, -1, -4,
 -1, 4, -2, -4, 1, -2, 4, 2, 1, 4, 2, -3, -4, -1, 4, -2, -3, -4, -1, -2, -3],
 [-3, 2, 4, 2, 1, 4, 3, -4, 2, 2, 4, -3, 2, -4, -4, 1, 4, 3, -2, -3, -4, -1,
 4, -2, -3, -4, 2, 4, 1, -4, -3, -4, 2, 4, 3, 2, 1, 4, 3, -2, -3, -4, -1, 4,
 -2, -3, -4, -2, -2, 4, 1, 3, -4, 2, 4, 3, 3, 2, -4, -4, 1, -2, 4, -1, -3, -4,
 -1, 2, -4, -2, 3, 1, 4, -2, -1, 3, 2, -4, 1, 4, 3, -2, -3, -4, -1, 4, -2, -3,
 -4, -2, -2, 4, 2, 1, 4, 2, 1, 4, 3, 2, 4, 4, 3, -2, -3, -4, 1, -3, -4, -2,
 -1, -3, -2, -4, -2, -1]])
end
