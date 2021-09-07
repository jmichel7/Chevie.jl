"""
This is a port of some GAP3/VkCurve functionality on presentations and
finitely presented groups.

*finitely  presented groups*  are distinguished  from *group presentations*
which  are objects of their own. The  reason is that when a presentation is
changed  (e.g. simplified)  by Tietze  transformations, new  generators and
relators are introduced; thus all words in a finitely presented group would
also have to be changed if such a Tietze transformation were applied to the
group. Therefore, it is better to work separately with the presentation and
reflect only carefully the changes on the group.

In order to speed up the algorithms, the relators in a presentation are not
represented  by `AbsWord`s, but by lists  of positive or negative generator
numbers which we call *Tietze words*.

We  now  describe  the  available  functions;  for more information look at
individual  docstrings. The functions `Presentation` and `FpGroup` create a
presentation from a finitely presented group or, vice versa.

Since  it could be large, by default a presentation is printed as a summary
of  the number of generators, the number  of relators, and the total length
of  all  relators.  We  illustrate  below  functions  displaying  more of a
presentation.

```julia-repl
julia> @AbsWord a,b

julia> F=FpGroup([a,b])
FreeGroup(a,b)

julia> G=F/[a^2,b^7,comm(a,a^b),comm(a,a^(b^2))*inv(b^a)]
FreeGroup(a,b)/[a²,b⁷,a⁻¹b⁻¹a⁻¹bab⁻¹ab,a⁻¹b⁻²a⁻¹b²ab⁻²ab²a⁻¹b⁻¹a]

julia> P=Presentation(G)
<< presentation with 2 gens and 4 rels of total length 30 >>

julia> relators(P)
4-element Vector{AbsWord}:
 a²
 b⁷
 ab⁻¹abab⁻¹ab
 b⁻²ab²ab⁻²ab²ab⁻¹
```

```julia-rep1
julia> Presentations.PrintGenerators(P)
#I 1. a 10 occurrences involution
#I 2. b 20 occurrences

julia> Presentations.display_balanced(P)
1: a=A
2: bbbbbbb=1
3: aBab=BAbA
4: BBabbaBBabbaB=1
```

for more information look at the help strings of 
AbsWord, Presentation, FpGroup, Go, GoGo, conjugate, 
       tryconjugate, simplify, relators, display_balanced
"""
module Presentations
using ..Gapjm
export AbsWord, @AbsWord, Presentation, FpGroup, Go, GoGo, conjugate, 
       tryconjugate, simplify, relators, display_balanced

"""
# Changing Presentations

The  functions `AddGenerator`, `AddRelator`, `RemoveRelator` can be used to
change a presentation. In general, they will change the isomorphism type of
the  group defined  by the  presentation, hence,  though they are sometimes
used  as subroutines by Tietze transformations functions like `Substitute`,
they do *not* perform Tietze transformations themselves.

# Group Presentations

The function `PresentationViaCosetTable`  can be used  to compute a
presentation for a concrete (e.,g. permutation or matrix) group.

# Tietze Transformations

The  functions described in this section can be used to modify a
group presentation by Tietze transformations.

In  general, the aim of such modifications  will be to *simplify* the given
presentation,  i.e., to reduce  the number of  generators and the number of
relators  without increasing too much the  sum of all relator lengths which
we  will  call  the  *total  length*  of the presentation. Depending on the
concrete presentation under investigation one may end up with a nice, short
presentation   or  with  a   very  huge  one.

There  is no way  to find the  shortest presentation which  can be obtained
from  a given  one. Therefore,  what we  offer are  some lower-level Tietze
transformation   functions  and,  in  addition,  a  heuristic  higher-level
function   (which  of  course   cannot  be  the   optimal  choice  for  all
presentations).

The design of these functions follows closely the concept of the ANU Tietze
transformation  program designed by George Havas cite{Hav69} which has been
available  from Canberra since 1977 in a stand-alone version implemented by
Peter   Kenne  and  James  Richardson  and   later  on  revised  by  Edmund
F.~Robertson (see cite{HKRR84}, cite{Rob88}).

The  higher-level  function  is  `simplify`.  The lower-level functions are
`Eliminate`, `Search`, `SearchEqual`, and `FindCyclicJoins`.

Some  of  these  functions  may  eliminate  generators,  but  they  do *not*
introduce  new generators. However,  sometimes you will  need to substitute
certain  words as  new generators  in order  to improve  your presentation.
Therefore  there are the  functions `Substitute` and `SubstituteCyclicJoins`
which introduce new generators.

Finally  the functions `InitGeneratorImages`  and `PrintGeneratorImages` can
be used to determine and to display the images or preimages of the involved
generators under the isomorphism which is defined by the sequence of Tietze
transformations which are applied to a presentation.

The  functions, `PrintPairs`, and `PrintOptions`,  can be useful. There are
also  the  *Tietze  options*:  parameters  which  essentially influence the
performance  of the  functions mentioned  above; they  are not specified as
arguments of function calls. Instead, they are stored in the presentation.
"""
plural(n,w)=string(n)*" "*w*(n==1 ? "" : "s")

#------------------ Abstract Words ----------------------------------
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

"A positive `AbsWord` is obtained by giving `Symbols` as arguments"
AbsWord(x::Symbol...)=AbsWord([s=>1 for s in x])

"`@AbsWord x,y,z` defines the variables `x,y,z` to be `AbsWord`s"
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
  if length(w.d)!=1 || last(w.d[1])!=1 error("not generator") end
  first(w.d[1])
end

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

function CoxGroups.word(G::FpGroup,w::AbsWord)
  res=Int[]
  for (s,m) in w.d
   p=findfirst(==(AbsWord(s)),gens(G))
   if isnothing(p) error(w," is not a word for ",G) end
   append!(res,fill(m<0 ? -p : p,abs(m)))
  end
  res
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
TietzeWord(w::AbsWord,gens::Vector{AbsWord})=TietzeWord(w,mon.(gens))

# hope this reflects faithfully the c function Relator in GAP3
function TietzeWord(w::AbsWord,ss::Vector{Symbol})
  res=Int[]
  for (s,c) in w.d
    p=findfirst(==(s),ss)
    if p===nothing error(s," is not in ",gens) end
    if c>0 append!(res,fill(p,c))
    else append!(res,fill(-p,-c))
    end
  end
  reduceword(res)
end

# reduce cyclically a Tietze word
function reduceword(w::Vector{Int})
  i=1;
  res=Int[]
  for j in eachindex(w)
    if isempty(res) || res[end]!=-w[j] if w[j]!=0 push!(res,w[j]) end
    else pop!(res)
    end
  end
  b=1;e=length(res)
  while b<=e && res[b]==-res[e] b+=1;e-=1 end
  res[b:e]
end
   
mutable struct TietzeStruct
  numgens::Int   # number of generators
  numrels::Int   # number of relators
  total::Int     # total length of relators
  generators::Vector{AbsWord} # copy of initial gens
  inverses::Vector{Int}
  relators::Vector{Vector{Int}}
  lengths::Vector{Int} # lengths of the relators
  flags::Vector{Int}
  modified::Bool
  numredunds::Int
  status::Vector{Int}
end

Base.getindex(T::TietzeStruct,i)=T.inverses[T.numgens+1-i]
Base.setindex!(T::TietzeStruct,j,i)=T.inverses[T.numgens+1-i]=j

# sorts by length relators after removing the empty ones
function Base.sort!(T::TietzeStruct)
  p=sortperm(T.relators,by=length)
  T.relators=T.relators[p]
  T.lengths=T.lengths[p]
  T.flags=T.flags[p]
  p=unique(i->T.relators[i],1:T.numrels)
  if isempty(T.relators[p[1]]) p=p[2:end] end
  T.relators=T.relators[p]
  T.lengths=T.lengths[p]
  T.flags=T.flags[p]
  T.numrels=length(T.lengths)
  T.total=sum(T.lengths)
end

@GapObj struct Presentation
end
# possible fields: tietze, operations, generators, components,
# nextFree, identity, eliminationsLimit, expandLimit, generatorsLimit,
# lengthLimit, loopLimit, printLevel, saveLimit, searchSimultaneous,
# protected

"""
PrintStatus(P [, norepeat=false]) . . .  print status line

`PrintStatus`  prints the current status of a presentation `P`,
i.e., the number of generators,  the  number of relators, and  the  total
length of all relators.

If  `norepeat`  is true,  then the printing is
suppressed if none of the three values has changed since the last call.
"""
function PrintStatus(P::Presentation;norepeat=false)
  T=P.tietze
  status=[T.numgens-T.numredunds, T.numrels, T.total]
  if !(status==T.status && norepeat)
    T.status=status
    print("Presentation:")
    if haskey(P, :name) print(" ",P.name)end
    print(" ",plural(status[1],"generator"))
    print(", ",plural(status[2],"relator"))
    print(", total length ",status[3])
    println()
  end
end

function Base.show(io::IO, T::Presentation)
  t=T.tietze
  print(io,"Presentation: ",plural(t.numgens-t.numredunds,"generator"),
        ", ",plural(t.numrels,"relator"),", total length ",t.total)
end

function Base.sort!(P::Presentation)
  if P.printLevel>=3 print("#I  sorting the relators\n") end
  sort!(P.tietze)
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

"""
`AddGenerator( P[, generator])`

`AddGenerator` adds a new generator to the list of generators.

If you don't specify a second argument, then  `AddGenerator` will define
a  new  abstract  generator  `_xi`  and save  it  in a  new  component
`P.i` of  the  given  presentation  where  `i`  is the  least
positive  integer which  has  not  yet  been used as a  generator number.
Though this  new generator will be printed as `_xi`,  you will have to
use the external variable `P.i` if you want to access it.

If you  specify a second  argument, then `generator` must be  an abstract
generator which does  not yet occur in the presentation.   `AddGenerator`
will add it to the presentation and save  it in a new component `P.i`
in the same way as described for `_xi` above.
"""
function AddGenerator(P::Presentation,gen=nothing)
  Check(P)
  T=P.tietze
  gens=T.generators
  if gen===nothing
    gen=NewGenerator(P)
    if P.printLevel>=1 println("#I AddGenerator: new generator is ", gen) end
  else
    if !(gen isa AbsWord && length(gen)==1)
      error("second argument must be an abstract generator")
    end
    if gen in gens || inv(gen) in gens
      println("#I generator ",gen," is already in the presentation")
      return
    end
    new=P.nextFree
    while Symbol(new) in keys(P.prop) new+=1 end
    P.nextFree=new+1
    push!(gens,gen)
    P[Symbol(string(new))]=gen
    push!(P.components,new)
    T.numgens+=1
    T.inverses=vcat([T.numgens],T.inverses, [-T.numgens])
  end
  if haskey(P, :imagesOldGens) StopTracingGeneratorImages(P) end
  T.modified=true
end

"""
`AddRelator(P, word)`

`AddRelator` adds  the word `word` to the list of  relators.  `word` must
be a word in the generators of the given presentation.
"""
function AddRelator(P::Presentation, word)
  Check(P)
  if P.printLevel>=3 print("#I  adding relator ", word, "\n") end
  T=P.tietze
  rel=TietzeWord(word,T.generators)
  if !isempty(rel)
    push!(T.relators,rel)
    push!(T.lengths,length(rel))
    push!(T.flags,1)
    T.numrels+=1
    T.total+=T.lengths[end]
    T.modified=true
  end
  if haskey(P,:imagesOldGens) StopTracingGeneratorImages(P) end
end

# reduce cyclically a Tietze word using the table of inverses
function reduceword(w::Vector{Int},T::TietzeStruct)
  i=1;
  res=Int[]
  for j in eachindex(w)
    if isempty(res) ||  (T[res[end]]!=-T[w[j]]  && T[w[j]]!=0) 
         push!(res,T[w[j]])
    elseif T[w[j]]!=0 pop!(res)
    end
  end
  b=1;e=length(res)
  while b<=e && res[b]==-res[e] b+=1;e-=1 end
  res[b:e]
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
function HandleLength1Or2Relators(P::Presentation)
  if P.printLevel>=3 print("#I  handling short relators\n") end
  T=P.tietze
  protected=P.protected
  tracingImages=haskey(P,:imagesOldGens)
  gens=T.generators
  rels=T.relators
  lengths=T.lengths
  flags=T.flags
  numgens=T.numgens
  numrels=T.numrels
  redunds=T.numredunds
  done=false
  while !done
    done=true
    i=0
    while i<numrels
      i+=1
      lg=lengths[i]
      if 0<lg<=2 && flags[i]<=2
        rep1=rels[i][1]
        if T[rep1]!=rep1 rep1=T[rep1] end
        if lg==1
          rep1=abs(rep1)
          if rep1>protected
            T[rep1]=T[-rep1]=0
            if P.printLevel>=2 println("#I  Handle12:eliminating ",gens[rep1]) end
            if tracingImages UpdateGeneratorImages(P, rep1, Int[]) end
            redunds+=1
            done=false
          end
        else
          rep2=rels[i][2]
          if T[rep2]!=rep2 rep2=T[rep2] end # an equivalence already declared
          if abs(rep2)<abs(rep1) rep1,rep2=rep2,rep1 end
          if rep1<0 rep1,rep2=(-rep1,-rep2) end
          if rep1==0 # already eliminated
            rep2=abs(rep2)
            if rep2>protected
              T[rep2]=T[-rep2]=0
              if P.printLevel>=2 println("#I Handle12:eliminating ",gens[rep2]) end
              if tracingImages UpdateGeneratorImages(P, rep2, Int[]) end
              redunds+=1
              done=false
            end
          elseif rep1!=-rep2  # otherwise not reduced
            if rep1!=rep2 # not an involution
              if T[rep2]==T[-rep2] && T[-rep1]<0 # rep2^2=1 => rep1^2=1
                numrels+=1
                push!(rels,[rep1, rep1])
                push!(lengths,2)
                push!(flags,1)
                T.numrels=numrels
                T.total+=2
              end
              if abs(rep2)>protected
                T[rep2]=T[-rep1]
                T[-rep2]=rep1
                if tracingImages || P.printLevel>=2
                  if rep2>0 rep1=T[-rep1] end
                  if P.printLevel>=2
                    println("#I Handle12:eliminating ",gens[abs(rep2)],"=",
                                        AbsWord([rep1],T.generators))
                  end
                  if tracingImages
                    UpdateGeneratorImages(P,abs(rep2),[rep1])
                  end
                end
                redunds=redunds+1
                done=false
              end
            elseif T[-rep1]<0  # an involution not yet detected
              rels[i]=[rep1,rep1]
              flags[i]=3
              T[-rep1]=rep1
              done=false
            end
          end
        end
      end
    end
    if !done
      for i in 1:numgens
        if T[i]!=i T[i],T[-i]=T[T[i]],T[-T[i]] end
      end
#     @show T.inverses
# the next loop should be FunTzReplaceGens
      for i in eachindex(T.relators)
        T.total-=T.lengths[i]
        T.relators[i]=reduceword(T.relators[i],T)
        T.lengths[i]=length(T.relators[i])
        T.flags[i]=1
        T.total+=T.lengths[i]
      end
    end
  end
  T.numredunds=redunds
  if redunds>0 RemoveGenerators(P) end
end

"""
`FpGroup(P::Presentation)`

returns the finitely presented group defined  by the presentation `P`.
"""
function FpGroup(P::Presentation)
  if P.printLevel>=3 
    println("#I  converting the Tietze presentation to a group")
  end
  Check(P)
  T=P.tietze
  if T.numredunds>0 RemoveGenerators(P) end
  sort!(T)
  gens=deepcopy(T.generators)
  G=FpGroup(gens,map(i->AbsWord(i,gens),T.relators))
  if haskey(P, :imagesOldGens)
    G.imagesOldGens=deepcopy(P.imagesOldGens)
    G.preImagesNewGens=deepcopy(P.preImagesNewGens)
  end
  G
end

Presentation(s::Vector{Symbol},r,pl=1)=Presentation(AbsWord.(s),r,pl)

function Presentation(ggens::Vector{AbsWord},grels::Vector{AbsWord},printlevel::Int=1)
  numgens=length(ggens)
  numrels=length(grels)
  rels=Vector{Int}[]
  lengths=Int[]
  total=0
  for i in grels
    rel=TietzeWord(i,ggens)
    push!(rels,rel)
    l=length(rel)
    push!(lengths,l)
    total+=l
  end
  T=TietzeStruct(numgens, numrels, total, deepcopy(ggens),
    numgens:-1:-numgens, rels, lengths, fill(1,numrels), false, 0, [0, 0, -1])
  P=Presentation(Dict{Symbol,Any}())
  P.generators=T.generators
  P.tietze=T
  P.components=collect(1:numgens)
  for i in 1:numgens
    setproperty!(P,Symbol(i),ggens[i])
  end
  P.nextFree=numgens+1
  P.identity=one(AbsWord)
  P.eliminationsLimit=100
  P.expandLimit=150
  P.generatorsLimit=0
  P.factor=100
  P.lengthLimit=100000000 # infinity
  P.loopLimit=100000000 # infinity
  P.printLevel=printlevel
  P.saveLimit=10
  P.searchSimultaneous=20
  if P.printLevel>=2 PrintStatus(P;norepeat= true) end
  P.protected=T.numgens
  HandleLength1Or2Relators(P)
  P.protected=0
  sort!(P)
  if P.printLevel>=2 PrintStatus(P;norepeat= true) end
  return P
end

"""
`Presentation( G::FpGroup[, printlevel])`

`Presentation` returns a presentation containing a copy of the presentation
of the given finitely presented group `G` on the same set of generators.

The  optional `printlevel` parameter  can be used  to restrict or to extend
the  amount of output provided by Tietze transformation functions when being
applied  to the created  presentation. The default  value 1 is designed for
interactive  use and implies  explicit messages to  be displayed by most of
these  functions. A  `printlevel` value  of 0  will suppress these messages,
whereas a `printlevel` value of 2 will enforce some additional output.
"""
Presentation(G::FpGroup,printlevel::Integer=1)=Presentation(G.gens,G.rels)

# takes as input the output of display_balanced
function Presentation(s::String)
  s=replace(s,r"\n\s*="=>"=")
  l=split(s,"\n")
  l=map(s->replace(s,r"^ *[0-9]*: *"=>""),l)
  l=filter(x->match(r"^ *$",x)===nothing,l)
  rels=AbsWord.(l)
  l=map(x->first.(x.d),rels)
  atoms=length(l)==1 ? unique(l[1]) : union(map(x->first.(x.d),rels)...)
  sort!(atoms)
  Presentation(AbsWord.(atoms),rels)
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
  local   F,            # given free group
          words,        # given words for the generators of H
          H,            # subgroup, if specified
          F1,           # f.p. group isomorphic to H
          F2,           # f.p. group isomorphic to G
          ng1,          # position number of identity element in G
          nh1,          # position number of identity element in H
          perms,        # permutations induced by the gens on the cosets
          stage,        # 1 or 2
          table,        # columns in the table for gens
          rels,         # representatives of the relators
          relsGen,      # relators sorted by start generator
          subgroup,     # rows for the subgroup gens
          i, j,         # loop variables
          gen,          # loop variables for generator
          gen0, inv0,   # loop variables for generator cols
          g, g1,        # loop variables for generator cols
          c,            # loop variable for coset
          rel,          # loop variables for relator
          rels1,        # list of relators
          app,          # arguments list for `MakeConsequences`
          index,        # index of the table
          col,          # generator col in auxiliary table
          perm,         # permutations induced by a generator on the cosets
          gens,         # abstract gens in which the relators are written
          gens2,        # the above abstract gens and their inverses
          ggens,        # concrete generators of G
          ngens,        # number of generators of G
          ngens2,       # twice the above number
          order,        # order of a generator
          actcos,       # part 1 of Schreier vector of G by H
          actgen,       # part 2 of Schreier vector of G by H
          tab0,         # auxiliary table in parallel to table <table>
          cosRange,     # range from 1 to index (= number of cosets)
          genRange,     # range of the odd integers from 1 to 2*ngens-1
          geners,       # order in which the table cols are worked off
          next,         # local coset number
          w,            # loop variable
          left1,        # part 1 of Schreier vector of H by trivial group
          right1,       # part 2 of Schreier vector of H by trivial group
          n,            # number of subgroup element
          words2,       # words for the generators of H and their inverses
          h             # subgroup element
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

|    gap> H := Group(
    >  [ (2,5,3), (2,7,5), (1,8,4), (1,8,6), (4,8,6), (3,5,7) ], () );;
    gap> P := PresentationViaCosetTable( H );
    << presentation with 6 gens and 12 rels of total length 42 >>
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
  local   F,          # given free group
          fgens,      # generators of F
          fwords,     # given words in the generators of F
          words,      # tidied up words for the subgroup generators
          H,          # subgroup
          elts,       # elements of G or H
          cosets,     # right cosets of G with respect to H
          F1,         # f.p. group isomorphic to H
          F2,         # f.p. group isomorphic to G
          P,          # resulting presentation
          ggens,      # concrete generators of G
          hgens,      # concrete generators of H
          ngens,      # number of generators of G
          h, w        # loop variables
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
  P.printLevel=1
  return P
end

"""
`RemoveRelator(P,n)`

`RemoveRelator`  removes the `n`th relator and  then resorts the  list of
relators in the given presentation `P`.
"""
function RemoveRelator(P::Presentation, n)
  Check(P)
  T=P.tietze
  numrels=T.numrels
  if !(n in 1:numrels) error("relator number out of range") end
  rels=T.relators
  lengths=T.lengths
  if P.printLevel>=3 print("#I  removing the ", n, "th relator\n") end
  leng=lengths[n]
  if leng==2 && rels[n][1]==rels[n][2]
    num=rels[n][1]
    if num<0 num=-num end
    if T[-num]==num T[-num]=-num end
  end
  rels[n]=Int[]
  lengths[n]=0
  T.total=T.total-leng
  sort!(P)
  if P.printLevel>=2 PrintStatus(P;norepeat= true) end
  if haskey(P, :imagesOldGens) StopTracingGeneratorImages(P) end
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

#F  Check(presentation)  . . . .  check presentation components
##
##  `Check`  checks some components of the  given presentation to be
##  consistent.
function Check(P::Presentation)
  T=P.tietze
  if !(P.generators===T.generators) || length(T.generators)!=T.numgens 
    error("inconsistent generator lists")
  end
  if length(T.inverses)!=2*T.numgens+1
    error("inconsistent generator inverses")
  end
  if length(T.relators)!=T.numrels || length(T.lengths)!=T.numrels ||
     length(T.flags)!=T.numrels error("inconsistent relators")
  end
  if any(i->length(T.relators[i])!=T.lengths[i],1:T.numrels)
    error("inconsistent lengths")
  end
  if sum(length,T.relators)!=T.total error("inconsistent total") end
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
  Check(T)
  tietze=T.tietze
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
    if tietze.numredunds>0 RemoveGenerators(T) end
    HandleLength1Or2Relators(T)
    sort!(T)
    if T.printLevel>=1 PrintStatus(T;norepeat= true) end
  else
    eliminations=T.eliminationsLimit;T.eliminationsLimit=n
    numgenslimit=T.generatorsLimit;T.generatorsLimit=tietze.numgens-n
    EliminateGens(T)
    T.eliminationsLimit=eliminations
    T.generatorsLimit=numgenslimit
  end
end

# substitute gen by w
function SubstituteGen(T::TietzeStruct, gen, w)
  iw=-reverse(w)
  for i in eachindex(T.relators)
    res=Int[]
    for j in T.relators[i]
      if j==gen append!(res,w)
      elseif j==-gen append!(res,iw)
      else push!(res,j)
      end
    end
    T.relators[i]=reduceword(res,T)
    T.total-=T.lengths[i]
    T.lengths[i]=length(T.relators[i])
    T.total+=T.lengths[i]
  end
end

"""
EliminateGen( presentation,n) . . . eliminates the nth generator

`EliminateGen` eliminates the Tietze generator tietze.generators[n]
if possible, i. e. if that generator can be isolated  in some appropriate
Tietze relator.  However,  the elimination  will not be  performed if the
resulting total length of the relators cannot be guaranteed to not exceed
the parameter P.lengthLimit.
"""
function EliminateGen(P::Presentation, num::Int)
  T=P.tietze
# spacelimit=P.lengthLimit
  spacelimit=P.total+P.factor
  T.modified=false
  gens=T.generators
  numgens=T.numgens
  rels=T.relators
  numrels=T.numrels
  lengths=T.lengths
  if !(0<num<=numgens)
    error("EliminateGen: second argument is not a valid generator number")
  end
  occur=Occurrences(T, num)
  occTotal=occur[1][1]
  if occTotal>0 && occur[3][1]==1
    occRelNum=occur[2][1]
    length=lengths[occRelNum]
    space=(occTotal-1)*(length-1)-length
    if T.total+space>spacelimit return end
    gen=num
    rel=rels[occRelNum]
    length=lengths[occRelNum]
    pos=findfirst(==(gen),rel)
    if pos===nothing
      gen=-gen
      pos=findfirst(==(gen),rel)
    end
    word=vcat(rel[pos+1:end],rel[1:pos-1])
    if P.printLevel >= 2
      print("#I EliminateGen: eliminating ", gens[num], "=")
      if gen>0 println(inv(AbsWord(word,T.generators)))
      else println(AbsWord(word,T.generators))
      end
    end
    SubstituteGen(T,-gen, word)
    if haskey(P, :imagesOldGens)
      if gen>0 word=-reverse(word) end
      UpdateGeneratorImages(P, num, word)
    end
    T[num]=0
    T.numredunds+=1
    T.modified=true
  end
end

# returns 3 vectors:
# 1: nb. occurences of each generator
# 2: in which relator occur with minimal multiplicity
# 3: the corresponding multiplicity
function Occurrences(T::TietzeStruct,gen=0)
  occ=fill(0,T.numgens)
  min=fill(0,T.numgens)
  mult=fill(0,T.numgens)
  for j in 1:T.numrels
    lococc=fill(0,T.numgens)
    for i in T.relators[j] lococc[abs(i)]+=1 end
    occ+=lococc
    for i in 1:T.numgens
      if lococc[i]>0
        if mult[i]==0 || lococc[i]<mult[i]
          min[i]=j
          mult[i]=lococc[i]
        end
      end
    end
  end
  gen==0 ? [occ,min,mult] : [[occ[gen]],[min[gen]],[mult[gen]]]
end

"""
EliminateGen1(presentation)  . . . . . . .  eliminates a generator

`EliminateGen1`  tries to  eliminate a  Tietze generator:  If there are
Tietze generators which occur just once in certain Tietze relators,  then
one of them is chosen  for which the product of the length of its minimal
defining word  and the  number of its  occurrences  is minimal.  However,
the elimination  will not be performed  if the resulting  total length of
the  relators   cannot  be  guaranteed   to  not  exceed   the  parameter
P.lengthLimit.
"""
function EliminateGen1(P::Presentation)
  T=P.tietze
  protected=P.protected
# spacelimit=P.lengthLimit
  spacelimit=T.total+P.factor
  gens=T.generators
  numgens=T.numgens
  rels=T.relators
  numrels=T.numrels
  lengths=T.lengths
  occur=Occurrences(T)
  occTotals=occur[1]
  occRelNums=occur[2]
  occMultiplicities=occur[3]
  modified=false
  num=0
  space=0
  for i in protected+1:numgens
    if occMultiplicities[i]==1
      total=occTotals[i]
      length=lengths[occRelNums[i]]
      ispace=(total-1)*(length-1)-length
      if num==0 || ispace<=space
        num=i
        space=ispace
      end
    end
  end
  if num>0 && T.total+space<=spacelimit
    gen=num
    occRelNum=occRelNums[num]
    rel=rels[occRelNum]
    length=lengths[occRelNum]
    pos=findfirst(==(gen),rel)
    if pos===nothing
      gen=-gen
      pos=findfirst(==(gen),rel)
    end
    word=vcat(rel[pos+1:length],rel[1:pos-1])
    if P.printLevel>=2
      print("#I EliminateGen1: eliminating ", gens[num], "==")
      if gen>0 println(inv(AbsWord(word,T.generators)))
      else     println(AbsWord(word,T.generators))
      end
    end
    SubstituteGen(T, -gen, word)
    if haskey(P, :imagesOldGens)
      if gen>0 word=-reverse(word) end
      UpdateGeneratorImages(P, num, word)
    end
    T[num]=0
    T.numredunds+=1
    modified=true
  end
  T.modified=modified
end

"""
EliminateGens(presentation)  .  Eliminates generators

`EliminateGens`  repeatedly eliminates generators from the presentation
of the given group until at least one  of  the  following  conditions  is
violated:

(1) The  current  number of  generators  is greater  than  the  parameter
    T.generatorsLimit.
(2) The   number   of   generators   eliminated   so  far  is  less  than
    the parameter T.eliminationsLimit.
(3) The  total length of the relators  has not yet grown  to a percentage
    greater than the parameter T.expandLimit.
(4) The  next  elimination  will  not  extend the total length to a value
    greater than the parameter T.lengthLimit.

the function will not eliminate any protected generators.
"""
function EliminateGens(T::Presentation)
  Check(T)
  if T.printLevel>=3 print("#I  eliminating generators\n") end
  tietze=T.tietze
  redundantsLimit=5
  maxnum=T.eliminationsLimit
  bound=tietze.total * T.expandLimit // 100
  modified=false
  tietze.modified=true
  num=0
  while tietze.modified && (num<maxnum && (tietze.total <=
       bound && tietze.numgens -tietze.numredunds>T.generatorsLimit))
    EliminateGen1(T)
    if tietze.numredunds==redundantsLimit RemoveGenerators(T) end
    modified=modified || tietze.modified
    num+=1
  end
  tietze.modified=modified
  if tietze.numredunds>0 RemoveGenerators(T)
  end
  if modified
    HandleLength1Or2Relators(T)
    sort!(T)
    if T.printLevel>=2 PrintStatus(T;norepeat= true) end
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
  if P.printLevel>=3 print("#I  searching for cyclic joins\n") end
  Check(P)
  T=P.tietze
  T.modified=false
  newstart=true
  while newstart
    exponents=GeneratorExponents(P)
    if sum(exponents)==0 return end
    newstart=false
    gens=T.generators
    numgens=T.numgens
    rels=T.relators
    numrels=T.numrels
    lengths=T.lengths
    flags=T.flags
    i=0
    while i<numrels
      i+=1
      rel=rels[i]
      if lengths[i]==4 && (rel[1]==T[-rel[3]] && rel[2]==T[-rel[4]])
        num=[abs(rel[1]), abs(rel[2])]
        exp=[exponents[num[1]], exponents[num[2]]]
        fac=[0, 0]
        e=[0, 0]
        if exp[1]>0 || exp[2]>0
          j=0
          while j<numrels
            j=j+1
            if lengths[j]>0 && j != i
              rel=rels[j]
              length=lengths[j]
              e[1]=0
              e[2]=0
              powers=0
              prev=0
              l=0
              while l<length
                l+=1
                next=rel[l]
                if next != prev
                    powers+=1
                    prev=next
                end
                if next==num[1] e[1]+=1
                elseif next==num[2] e[2]+=1
                elseif next==-num[1] e[1]-=1
                elseif next==-num[2] e[2]-=1
                else l=length+1
                end
              end
              if l==length && powers>1
                for k=[1, 2]
                    fac[k]=num[k]
                    if exp[k]>0
                        e[k]=mod(e[k], exp[k])
                        if e[k]>exp[k] // 2
                            e[k]=exp[k]-e[k]
                            fac[k]=-(fac[k])
                        end
                    elseif e[k]<0
                        e[k]=-(e[k])
                        fac[k]=-(fac[k])
                    end
                    if fac[k]<0 fac[k]=T[fac[k]] end
                end
                for k=[1, 2]
                    if e[k]>0 && e[3-k]==0
                        exp[k]=gcd(e[k], exp[k])
                        if exp[k] != exponents[num[k]]
                            exponents[num[k]]=exp[k]
                            e[k]=exp[k]
                        end
                    end
                end
                if e[1]+e[2]<length || powers>2
                  rel=[]
                  if e[1]>0 rel=Concatenation(rel, fac[1]+fill(0,e[1])) end
                  if e[2]>0 rel=Concatenation(rel, fac[2]+fill(0,e[2])) end
                  rels[j]=rel
                  lengths[j]=e[1]+e[2]
                  T.total=(T.total-length)+lengths[j]
                  flags[j]=1
                  T.modified=true
                  if P.printLevel >= 3
                      print("#I  rels[", j, "] reduced to ", rels[j], "\n")
                  end
                end
                if e[1]==1 n=num[1]
                elseif e[2]==1 n=num[2]
                else
                  n=0
                  for k=[1, 2]
                    if n==0 && (e[k]>1 && gcd(e[k], exp[k])==1)
                      ggt=Gcdex(e[k], exp[k])
                      gen=[gens[num[1]], gens[num[2]]]
                      if fac[1]<0 gen[1]=gen[1]^-1 end
                      if fac[2]<0 gen[2]=gen[2]^-1 end
                      word=gen[k] * gen[3-k]^(e[3-k] * ggt[:coeff1])
                      AddRelator(P, word)
                      numrels=T.numrels
                      n=num[k]
                    end
                  end
                end
                if n != 0 && P.generatorsLimit<numgens
                  Eliminate(P)
                  T.modified=true
                  j=numrels
                  i=numrels
                  if 1<numgens newstart=true end
                end
              end
            end
          end
        end
      end
    end
  end
  if T.modified
      HandleLength1Or2Relators(P)
      sort!(P)
      if P.printLevel>=1 PrintStatus(P;norepeat= true) end
  end
end

"""
GeneratorExponents(presentation) . . . list of generator exponents

`GeneratorExponents`  tries to find exponents for the Tietze generators
and return them in a list parallel to the list of the generators.
"""
function GeneratorExponents(T::Presentation)
  if T.printLevel>=3 println("#I  trying to find generator exponents") end
  tietze=T.tietze
  numgens=tietze.numgens
  rels=tietze.relators
  numrels=tietze.numrels
  lengths=tietze.lengths
  flags=tietze.flags
  exponents=fill(0,numgens)
  for i in 1:numrels
      if lengths[i]>0
          rel=rels[i]
          length=lengths[i]
          num1=rel[1]
          j=2
          while j <= length && rel[j]==num1 j+=1 end
          if j>length
              num=abs(num1)
              if exponents[num]==0 exp=length
              else exp=gcd(exponents[num], length)
              end
              exponents[num]=exp
              if exp<length
                  rels[i]=num+fill(0,exp)
                  lengths[i]=exp
                  tietze.total=(tietze.total-length)+exp
                  flags[i]=1
                  tietze.modified=true
              elseif num1<0
                  rels[i]=-rel
              end
          end
      end
  end
  return exponents
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
    gap> P.primaryGeneratorWords;
    [ b, a*b*a ]
    gap> P.protected := 2;;
    gap> P.printLevel := 2;;
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
`Go` in case of T.printLevel=1 is suppressed.
"""
function Go(P::Presentation,silent=false)
  printstatus=P.printLevel==1 && !silent
  Check(P)
  T=P.tietze
  Search(P)
  looplimit=P.loopLimit
  count=0
  while count<looplimit && T.total>0
    SearchEqual(P)
    if T.modified Search(P) end
    EliminateGens(P)
    if T.modified
      Search(P)
      count+=1
    else count=looplimit
    end
    if printstatus PrintStatus(P;norepeat= true) end
  end
  if T.total>0
    FindCyclicJoins(P)
    if T.modified Search(P) end
    if printstatus PrintStatus(P;norepeat= true) end
  end
end

simplify(P::Presentation)=GoGo(P)

"""
`GoGo(P)`

`GoGo` performs Tietze transformations on a presentation `P`. It repeatedly
calls  the  `Go`  command  until  neither  the number of generators nor the
number of relators nor the total length of all relators have changed during
five consecutive calls of `Go`.
"""
function GoGo(T::Presentation)
  tietze=T.tietze
  numgens=tietze.numgens
  numrels=tietze.numrels
  total=tietze.total
  silentGo=T.printLevel==1
  count=0
  while count<5
    Go(T, silentGo)
    count+=1
    if silentGo && (tietze.numgens<numgens || tietze.numrels<numrels)
        PrintStatus(T;norepeat= true)
    end
    if tietze.numgens<numgens || (tietze.numrels<numrels ||
                                    tietze.total<total)
      numgens=tietze.numgens
      numrels=tietze.numrels
      total=tietze.total
      count=0
    end
  end
  if silentGo PrintStatus(T;norepeat= true) end
end

"""
`InitGeneratorImages(P)`
Any  sequence of  Tietze transformations  applied to  a presentation
`P`,  starting  from  an  ``old''  presentation  `P₁`  and ending up with a
``new''  presentation `P₂`,  defines an  isomorphism, `φ`  say, between the
groups defined by `P₁` and `P₂`, respectively. Sometimes it is desirable to
know  the  images  of  the  old  generators  or  the  preimages  of the new
generators  under `φ`. The {GAP}  Tietze transformations functions are able
to  trace these images. This is not automatically done because the involved
words  may grow to tremendous length, but it will be done if you explicitly
request for it by calling the function `InitGeneratorImages`.

`InitGeneratorImages` initializes three components of `P`:

`P.oldGenerators`: 
        This is the  list of the old generators.  It is initialized  by a
        copy of the current list of generators, `P.generators`.

`P.imagesOldGens`: 
        This  will be  the list  of the images  of the  old generators as
        Tietze words in the new generators. For each generator `g_i`, the
        `i`-th entry of the list is initialized by the Tietze word |[i]|.

`P.preImagesNewGens`: 
        This will be the  list of the  preimages of the new generators as
        Tietze words in the old generators. For each generator `g_i`, the
        `i`-th entry of the list is initialized by the Tietze word |[i]|.

This means, that  `P₁` is  defined  to be the  current presentation  and
`φ` to be the identity on `P₁`. From now  on, the existence of the
component  `P.imagesOldGens`   will  cause the  Tietze  transformations
functions to update  the lists of images and  preimages whenever they are
called.

You can  reinitialize the tracing of  the  generator images at  any later
state by just calling the function `InitGeneratorImages` again. For, if
the   above components do   already exist when `InitGeneratorImages` is
being called, they will first be deleted and then initialized again.

There    are  a few restrictions    concerning  the tracing  of generator
images:

In   general,    the    functions    `AddGenerator`,  `AddRelator`,   and
`RemoveRelator`  described  in section  "Changing  Presentations"  do not
perform Tietze transformations as they may change the isomorphism type of
the presentation.  Therefore, if any of them is called for a presentation
in which generator images and preimages are  being traced, it will delete
these lists.
"""
#F  InitGeneratorImages( T )  . . . . . . . initialize the generator images
#F                                               under Tietze transformations
##
##  `InitGeneratorImages`  expects  T  to  be  a  presentation.  It
##  defines the current generators to be the "old" generators and initializes
##  two components T.imagesOldGens and T.preImagesNewGens which, at any time,
##  will contain a list of the images  of the old generators  as Tietze words
##  in the new ones  or a list of the  preimages of the new generators in the
##  old ones,  respectively,  under the  isomorphism  defined  by all  Tietze
##  transformations which then will have been applied to T.
##
function InitGeneratorImages(T::Presentation)
  Check(T)
  numgens=length(T.generators)
  T.oldGenerators=deepcopy(T.generators)
  T.imagesOldGens=map(i->[i],1:numgens)
  T.preImagesNewGens=map(i->[i],1:numgens)
end

"""
MostFrequentPairs( presentation, n ) . . . .  occurrences of pairs

`MostFrequentPairs`  returns a list  describing the  n  most frequently
occurring relator subwords of the form  g1*g2,  where  g1  and  g2 are
different generators or their inverses. n=0 interpreted as infinity
"""
function MostFrequentPairs(P::Presentation, nmax::Int)
  if nmax<=0 error("second argument must be≥0") end
# counts occurrences of pairs gen,geni or their inverses (where 0<gen<geni)
# circularly in the relators or their inverses
# returns a (4,numgens) matrix; each line is zero until gen+1
# the first line  counts the occurrences of (gen,geni)
# the second             the occurrences of (gen,-geni)
# the third              the occurrences of (-gen,geni)
# the last               the occurrences of (-gen,-geni)
function OccurrencesPairs(T::TietzeStruct,gen)
  res=fill(0,4,T.numgens)
  for rel in T.relators, i in 1:length(rel)
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
  T=P.tietze
  gens=T.generators
  numgens=T.numgens
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
      if n>nmax
        sort!(pairs,rev=true)
        pairs=pairs[1:nmax]
        n=nmax
      end
    end
    sort!(pairs,rev=true)
  end
  return pairs
end

"""
NewGenerator(presentation) . . . . . . . . .  adds a new generator

`NewGenerator`  defines a  new  abstract generator  and adds it  to the
given presentation.

Let  i  be the smallest positive integer  which has not yet been used  as
a generator number  and for which no component  T.i  exists so far in the
given  presentation  T,  say.  A new abstract generator  _xi  is defined
and then added as component  T.i  to the given presentation.

Warning:  `NewGenerator`  is  an  internal  subroutine  of  the  Tietze
routines.  You should not call it.  Instead, you should call the function
`AddGenerator`, if needed.
"""
function NewGenerator(P::Presentation)
  T=P.tietze
  numgens=T.numgens
  new=P.nextFree
  recfields=keys(P.prop)
  while Symbol(new) in recfields new+=1 end
  P.nextFree=new+1
  string=Symbol("_x",new)
  gen=AbsWord(string)
  P.prop[string]=gen
  numgens=numgens+1
  push!(T.generators,gen)
  push!(P.components,new)
  T.numgens=numgens
  T.inverses=vcat([numgens], T.inverses, [-numgens])
  gen
end

"""
`Print(P[,list])`

`Print` provides a  kind of  *fast print out* for a presentation `P`.

Remember that in order to speed up the Tietze transformation routines, each
relator  in  a  presentation  `P`  is  internally  represented by a list of
positive  or negative  generator numbers,  i.e., each  factor of the proper
{GAP}  word  is  represented  by  the  position number of the corresponding
generator  with  respect  to  the  current  list  of  generators, or by the
respective  negative number,  if the  factor is  the inverse of a generator
which  is  not  known  to  be  an  involution.  In contrast to the commands
`relators`  and `PrintPresentation` described  above, `Print` does not
convert these lists back to the corresponding {GAP} words.

`Print`  prints the current  list of generators,  and then for each relator
its  length  and  its  internal  representation  as  a  list of positive or
negative generator numbers.

If a list `list` has been specified as second argument, then it is expected
to  be a list of the position numbers of the relators to be printed. `list`
need  not be  sorted and  may contain  duplicate elements. The relators are
printed  in  the  order  in  which  and  as often as their numbers occur in
`list`.  Position  numbers  out  of  range  (with  respect  to  the list of
relators) will be ignored.
"""
function Print(P::Presentation,list=1:P.tietze.numrels)
  Check(P)
  T=P.tietze
  if isempty(T.generators) print("#I  there are no generators\n")
  else print("#I  generators: ", T.generators, "\n")
  end
  if T.numrels==0 print("#I  there are no relators\n");return end
  print("#I  relators:\n")
  for i in list
    if i in 1:T.numrels 
      println("#I  ",i,".  ",T.lengths[i],"  ",T.relators[i])
    end
  end
end

"""
`PrintGenerators(P[,list])`

prints the current list of generators of a presentation
`P`,  providing  for  each  generator  its  name,  the  total number of its
occurrences  in the  relators, and,  if that  generator is  known to  be an
involution, an appropriate message.

a second argument `list` of generator indices prints only those generators.
"""
function PrintGenerators(P::Presentation,list::AbstractVector{Int}=1:P.tietze.numgens)
  Check(P)
  T=P.tietze
  gens=T.generators
  numgens=T.numgens
  if numgens==0 print("#I  there are no generators\n");return end
  if isempty(list) list=[0] end
  leng=length(list)
  _min=max(minimum(list),1)
  _max=min(maximum(list),numgens)
  if _min==_max
    occur=Occurrences(T, _max)
    num=occur[1][1]
    print("#I ", _max,". ",gens[_max]," ",plural(num,"occurrence"))
    if T[-_max]>0 print(" involution") end
    println()
  elseif _min<_max
    occur=Occurrences(T)
    for i in list
      if 1<=i<=numgens
        num=occur[1][i]
        print("#I ", i, ". ", gens[i], " ",plural(num,"occurrence"))
        if T[-i]>0 print(" involution") end
        println()
      end
    end
  end
end

"""
`PrintGeneratorImages(P)`

If `P` is   a presentation in  which  generator images and preimages  are
being  traced  through  all  Tietze  transformations   applied  to `P`,
`PrintGeneratorImages` prints the preimages  of the current  generators
as  Tietze   words in the   old  generators  and  the  images of  the old
generators as Tietze words in the current generators.

|    gap> G := PerfectGroup( 960, 1 );
    PerfectGroup(960,1)
    gap> P := Presentation( G );
    << presentation with 6 gens and 21 rels of total length 84 >>
    gap> InitGeneratorImages( P );
    gap> Go( P );
    &I  there are 3 generators and 11 relators of total length 96
    &I  there are 3 generators and 10 relators of total length 81
    gap> PrintGeneratorImages( P );
    &I  preimages of current generators as Tietze words in the old ones:
    &I  1. [ 1 ]
    &I  2. [ 2 ]
    &I  3. [ 4 ]
    &I  images of old generators as Tietze words in the current ones:
    &I  1. [ 1 ]
    &I  2. [ 2 ]
    &I  3. [ 1, -2, 1, 3, 1, 2, 1 ]
    &I  4. [ 3 ]
    &I  5. [ -2, 1, 3, 1, 2 ]
    &I  6. [ 1, 3, 1 ]
    gap> & Print the old generators as words in the new generators.
    gap> gens := P.generators;
    [ a, b, t ]
    gap> oldgens := P.oldGenerators;
    [ a, b, s, t, u, v ]
    gap> for i in [ 1 .. Length( oldgens ) ] do
    >  Print( oldgens[i], " = ",
    >  AbsWord( P.imagesOldGens[i], gens ), "\n" );
    >  od;
    a = a
    b = b
    s = a*b^-1*a*t*a*b*a
    t = t
    u = b^-1*a*t*a*b
    v = a*t*a |
"""
function PrintGeneratorImages(T::Presentation)
  Check(T)
  if haskey(T, :imagesOldGens)
    images=T.preImagesNewGens
    println("#I  preimages of current generators as Tietze words in the old ones:")
    for i in 1:length(images)
      print("#I  ", i, ". ", images[i], "\n")
    end
    images=T.imagesOldGens
    println("#I  images of old generators as Tietze words in the current ones:")
    for i in 1:length(images) print("#I  ",i,". ",images[i], "\n") end
  else
    print("#I  generator images are not available\n")
  end
end

"""
`PrintOptions(P)`

Several  of  the  Tietze  transformation  commands  described  above  are
controlled by  certain parameters, the *Tietze options*, which often have
a  tremendous  influence on their performance and  results.   However, in
each application of the  commands, an appropriate choice of these  option
parameters  will depend on the concrete presentation under investigation.
Therefore we have implemented the Tietze options in such  a way that they
are associated to the presentation:  Each  presentation
keeps its own  set  of  Tietze option parameters in the form of  ordinary
components.   In  particular, you may alter the  value  of any  of
these  Tietze  options by  just assigning  a  new value to the respective
component.

`PrintOptions`  prints the Tietze  option components of  the  specified
presentation `P`.

The Tietze options have the following meaning.

`protected`: 
        The first  `P.protected`  generators in a presentation  `P` are
        protected from  being  eliminated  by the  Tietze transformations
        functions.   There  are  only   two  exceptions:   The  option
        `P.protected`     is     ignored      by     the      functions
        `Eliminate(P,gen)` and  `Substitute(P,n,eliminate)`
        because they explicitly specify  the generator  to be eliminated.
        The default value of `protected` is 0.

`eliminationsLimit`: 
        Whenever the elimination  phase of  the `Go` command is entered
        for  a  presentation  `P`,   then  it   will  eliminate  at  most
        `P.eliminationsLimit` generators (except for further ones which
        have  turned  out  to   be  trivial).  Hence  you  may  use   the
        `eliminationsLimit` parameter as a break criterion for the `Go`
        command.  Note, however, that it  is ignored by the `Eliminate`
        command. The default value of `eliminationsLimit` is 100.

`expandLimit`: 
        Whenever the routine for eliminating  more than  1  generators is
        called for a presentation `P` by the `Eliminate` command or the
        elimination phase of the `Go` command, then it saves the  given
        total  length of the  relators,  and  subsequently it checks  the
        current total length against  its value before  each elimination.
        If the total length  has increased to more than `P.expandLimit`
        per cent of its original value, then the  routine returns instead
        of   eliminating  another  generator.   Hence  you  may  use  the
        `expandLimit` parameter  as a  break  criterion  for  the  `Go`
        command. The default value of `expandLimit` is 150.

`generatorsLimit`: 
        Whenever the elimination  phase of the  `Go` command is entered
        for  a  presentation  `P`  with  `n`  generators,  then  it  will
        eliminate at most `n-P.generatorsLimit`  generators (except
        for generators which turn out to  be trivial).  Hence you may use
        the  `generatorsLimit`  parameter  as  a break  criterion for the
        `Go` command. The default value of `generatorsLimit` is 0.

`lengthLimit`: 
        The  Tietze  transformation  commands  will  never  eliminate   a
        generator of a presentation  `P`,  if  they  cannot  exclude  the
        possibility  that  the  resulting total  length  of  the relators
        exceeds  the value of  `P.lengthLimit`.  The default  value  of
        `lengthLimit` is `infinity`.

`loopLimit`: 
        Whenever the `Go`  command  is called  for a  presentation `P`,
        then  it  will  loop  over  at  most `P.loopLimit` of its basic
        steps.  Hence you  may  use the  `loopLimit` parameter as a break
        criterion  for   the  `Go`   command.  The   default  value  of
        `loopLimit` is `infinity`.

`printLevel`: 
        Whenever   Tietze  transformation  commands  are  called  for   a
        presentation `P` with  `P.printLevel`  `=  0`,  they  will  not
        provide any output except for error messages. If `P.printLevel`
        `=  1`, they will display some  reasonable amount of output which
        allows you to watch the progress of the computation and to decide
        about your next commands. In the case `P.printLevel` `= 2`, you
        will  get  a much more  generous amount  of output.  Finally,  if
        `P.printLevel` `= 3`, various messages on internal details will
        be added. The default value of `printLevel` is 1.

`saveLimit`: 
        Whenever the  `Search` command has finished its main loop  over
        all relators of a presentation `P`, then it checks whether during
        this loop the total length of the relators has been reduced by at
        least  `P.saveLimit`  per  cent.  If  this is  the  case,  then
        `Search` repeats its procedure instead of returning.  Hence you
        may use the `saveLimit` parameter  as a  break  criterion for the
        `Search` command  and, in  particular,  for the search phase of
        the `Go` command. The default value of `saveLimit` is 10.

`searchSimultaneous`: 
        Whenever the `Search`  or the `SearchEqual` command is called
        for  a  presentation  `P`,  then it is  allowed  to handle  up to
        `P.searchSimultaneously` short relators simultaneously (see for
        the description of the `Search` command for more details).  The
        choice of this parameter may heavily influence the performance as
        well  as  the result of the  `Search`  and  the `SearchEqual`
        commands and  hence also  of  the  search  phase  of  the  `Go`
        command. The default value of `searchSimultaneous` is 20.

As soon as  a presentation has been defined, you  may
alter any of its Tietze option parameters at any time by just assigning a
new value to the respective component.

To demonstrate  the  effect of the `eliminationsLimit` parameter, we will
give  an example in which we handle a subgroup of index 240 in a group of
order 40320  given by  a presentation  due  to B.~H.  Neumann.   First we
construct a presentation of the subgroup, and  then  we  apply  to it the
`GoGo`   command  for  different   values  of  the  `eliminationsLimit`
parameter (including the default value 100).  In  fact, we also alter the
`printLevel` parameter, but this  is only done in order to suppress  most
of  the output.   In all  cases the  resulting  presentations  cannot  be
improved any more by applying the `GoGo`  command again, i.e., they are
the best results which we can get without substituting new generators.

|    gap> F3 := FreeGroup( "a", "b", "c" );;
    gap> G := F3 / [ F3.1^3, F3.2^3, F3.3^3, (F3.1*F3.2)^5,
            (F3.1^-1*F3.2)^5, (F3.1*F3.3)^4, (F3.1*F3.3^-1)^4,
            F3.1*F3.2^-1*F3.1*F3.2*F3.3^-1*F3.1*F3.3*F3.1*F3.3^-1,
            (F3.2*F3.3)^3, (F3.2^-1*F3.3)^4 ];;
    gap> a := G.1;;  b := G.2;;  c := G.3;;
    gap> H := Subgroup( G, [ a, c ] );;
    gap> P := PresentationSubgroup( G, H );
    << presentation with 224 gens and 593 rels of total length 2769 >>
    gap> for i in [ 28, 29, 30, 94, 100 ] do
    >       Pi := Copy( P );
    >       Pi.eliminationsLimit := i;
    >       Print( "&I  eliminationsLimit set to ", i, "\n" );
    >       Pi.printLevel := 0;
    >       GoGo( Pi );
    >       print( Pi );
    >    od;
    &I  eliminationsLimit set to 28
    &I  there are 2 generators and 95 relators of total length 10817
    &I  eliminationsLimit set to 29
    &I  there are 2 generators and 5 relators of total length 35
    &I  eliminationsLimit set to 30
    &I  there are 3 generators and 98 relators of total length 2928
    &I  eliminationsLimit set to 94
    &I  there are 4 generators and 78 relators of total length 1667
    &I  eliminationsLimit set to 100
    &I  there are 3 generators and 90 relators of total length 3289 |

Similarly,  we  demonstrate the influence of the `saveLimit` parameter by
just continuing the  preceding example  for some  different values of the
`saveLimit` parameter  (including  its  default  value  10),  but without
changing the `eliminationsLimit`  parameter which keeps its default value
100.

|    gap> for i in [ 9, 10, 11, 12, 15 ] do
    >       Pi := Copy( P );
    >       Pi.saveLimit := i;
    >       Print( "&I  saveLimit set to ", i, "\n" );
    >       Pi.printLevel := 0;
    >       GoGo( Pi );
    >       print( Pi );
    >    od;
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
"""
OptionNames=[:eliminationsLimit,:expandLimit,:generatorsLimit,:lengthLimit,
               :loopLimit,:printLevel,:saveLimit,:searchSimultaneous]
function PrintOptions(T::Presentation)
  len=0
  for nam in keys(T.prop)
    if nam in OptionNames
      if len<length(string(nam)) len=length(string(nam)) end
      lst=nam
    end
  end
  print("rec(\n",
    join(map(filter(nam->nam in OptionNames,collect(keys(T.prop))))do nam
    string("  ",nam," "^(len+1-length(string(nam))),"= ",getproperty(T,nam))
  end,",\n")," )")
end

"""
`PrintPairs(P[,n=10])`

`PrintPairs` determines  in  the  given presentation  `P`  the `n` most
frequently occurring squarefree relator subwords  of length  2 and prints
them together with their numbers of  occurrences.
A value `n=0` is interpreted as `infinity`.

This  list is a  useful piece of information  in the context of using the
`Substitute` command described above.
"""
function PrintPairs(T::Presentation,n=10)
  Check(T)
  if !(n isa Int) || n<0 error("second argument must be a ≥0 integer") end
  tietze=T.tietze
  gens=tietze.generators
  pairs=MostFrequentPairs(T, n)
  n=length(pairs)
  for m=1:n
    num=pairs[m][1]
    k=pairs[m][4]
    geni=gens[pairs[m][2]]
    if k>1 geni=geni ^ -1 end
    genj=gens[pairs[m][3]]
    if mod(k, 2)==1 genj=genj ^ -1 end
    if num==1
      print("#I  ",m,".  ", num, "  occurrence  of  ", geni, " * ", genj, "\n")
    elseif num>1
      print("#I  ",m,".  ", num, "  occurrences of  ", geni, " * ", genj, "\n")
    end
  end
end

"""
`PrintPresentation(P)`

`PrintPresentation` prints the current lists of generators and relators
and the current state of a presentation `P`.  In fact, the command
"""
function PrintPresentation(T::Presentation)
  Check(T)
  print("#I  generators:\n")
  PrintGenerators(T)
  print("#I  relators:\n")
  println(relators(T))
  print(T)
end

"""
`relators(P)`

prints the current list of relators  of  a presentation `P`.

If  a  list  `list`  has been specified  as  second argument, then it  is
expected  to  be a list  of the position numbers  of the  relators  to be
printed.  `list` need  not be sorted  and may contain duplicate elements.
The relators are printed as Tietze words in  the order  in  which (and as
often as)  their numbers occur in  `list`.  Position numbers out of range
(with respect to the list of relators) will be ignored.
"""
function relators(P::Presentation)
  Check(P)
  T=P.tietze
  map(i->AbsWord(T.relators[i],T.generators),1:T.numrels)
end

"""
RemoveGenerators(P::presentation) . . . . Remove redundant generators

`RemoveGenerators`  deletes the  redundant Tietze  generators and renumbers
the non-redundant ones accordingly. The redundant generators are assumed to
be marked in T.inverses list by an entry T[i]!=i.
"""
function RemoveGenerators(P::Presentation)
  if P.printLevel>=3 print("#I  renumbering the Tietze generators\n") end
  T=P.tietze
  redunds=T.numredunds
  if redunds==0 return end
  tracingImages=haskey(P, :imagesOldGens)
  if tracingImages preimages=P.preImagesNewGens end
  comps=P.components
  gens=T.generators
  numgens=T.numgens
  j=0
# @show j,T.inverses, T.numgens,T.numredunds
  for i in 1:numgens
    if T[i]==i
      j+=1
      if j<i
        comps[j]=comps[i]
        T[i]=j
        if T[-i]>0 T[-i]=j else T[-i]=-j end
      end
    else
      delete!(P.prop,Symbol(comps[i]))
      T[i]=0
      T[-i]=0
    end
  end
# @show j,T.inverses
  if j!=numgens-redunds
    error("This is a bug.  You should never get here.\n", 
          "Please send a copy of your job to the GAP administrators.\n")
  end
# @show j, T.inverses, numgens
  for rel in T.relators, i in eachindex(rel) rel[i]=T[rel[i]] end
  if tracingImages
    for i in 1:length(P.imagesOldGens)
      image=P.imagesOldGens[i]
      newim=[]
      for j=1:length(image) push!(newim, T[image[j]]) end
      P.imagesOldGens[i]=reduceword(newim)
    end
  end
  for i in 1:numgens
    j=T[i]
    if j<i && j>0
      gens[j]=gens[i]
      T[j]=j
      T[-j]=T[-i]
      if tracingImages preimages[j]=preimages[i] end
    end
  end
  j=numgens
  numgens1=numgens+1
  numgens-=redunds
  T.inverses=T.inverses[numgens1-numgens:numgens1+numgens]
  resize!(T.generators,numgens)
  resize!(P.components,numgens)
  while j>numgens
    if tracingImages Unbind(preimages[j]) end
    j-=1
  end
  T.numgens=numgens
  T.numredunds=0
end

#myhash(x::Int,s::UInt)=x+(s<<7)+s
myhash(x::Int,s::UInt)=hash(x,s)

# hash of cyclic subword of relator w starting at i and of length l
# works for i≠0, -length(w)≤i≤length(w)
function hashword(w,i,l,T)
  res=UInt(0)
  if i<0
    for j in -i:-1:max(1,-i+1-l) res=myhash(T[-w[j]],res) end
    for j in length(w):-1:length(w)-l-i+1 res=myhash(T[-w[j]],res) end
  else
    for j in i:min(i+l-1,length(w)) res=myhash(T[w[j]],res) end
    for j in 1:i+l-1-length(w) res=myhash(T[w[j]],res) end
  end
  res
end

modpos(i,l)=1+mod(i-1,l)
modneg(i,l)=mod(i,l)-l

# cyclic subword of relator w starting at i and of length l; takes i mod
# length(w) but remembers sign (positive: forward word, negative: -backward)
function subword(w,i,l)
  neg=sign(i)
  i=modpos(abs(i),length(w))
  if neg<0 vcat(-w[i:-1:max(1,i+1-l)],-w[length(w):-1:length(w)-l+i+1])
  else vcat(w[i:min(i+l-1,length(w))],w[1:i+l-1-length(w)])
  end
end

# finds intersection of 2 sorted lists; returns list of pairs (i1,i2)
# such that a[i1]==b[i2]
function intersectsorted(a,b;by=x->x)
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
   
# doreplace(w1,w2,pos1,pos2,l) replace in relator w2 substring of length l
# starting at pos2 with complement in w1 of identical substring of length l
# starting at pos1 in w1.
# return new w2
function doreplace(w1,w2,pos1,pos2,l)
  l1=length(w1)
  compl=subword(w1,pos1>0 ? modneg(1-pos1,l1) : modpos(-pos1+1,l1),l1-l)
  vcat(compl,subword(w2,pos2+l,length(w2)-l))
end

pr(w)=String(map(i->i>0 ? Char(i+96) : uppercase(Char(-i+96)),w))

function SearchC(T::TietzeStruct,i,j,equal=false)
# @show i,j,equal
  rels=T.relators
# println("SearchC(",join(pr.(rels[i:j]),","),",",equal,")")
  len=length(rels[i])
  if !equal l=div(len,2)+1;lmin=len-(len%2);lmax=lmin+1
  elseif isodd(len) error("searchequal and odd length")
  else l=div(len,2);lmin=len-(len%2);lmax=lmin
  end
# @show l,lmin,lmax
  hh=Tuple{UInt,Int,Int}[]
# if T.inverses!=T.numgens:-1:-T.numgens error("numgens=",T.numgens," invs=",
#                    T.inverses) end
  altered=Int[]
  for k in i:length(rels)
    rel=rels[k]
    if length(rel)<lmin continue end # sort should be called before SearchC
    rel=reduceword(rel,T)
    hk=Vector{Tuple{UInt,Int}}(undef,length(rel))
    for u in 1:length(rel) hk[u]=(hashword(rel,u,l,T),u) end
    sort!(hk,by=first)
    pp=intersectsorted(hh,hk;by=first)
    hk=map(((p1,p2),)->(hh[p1][2],hk[p2][2],hh[p1][3]),pp)
    if !isempty(hk) 
      hk=map(hk)do (pos1,pos2,w)
        (maxmatch(rels[w],pos1,rel,pos2,T),w,k)
      end
#     @show hk
      (pos1,pos2,l1),w,k=hk[argmax(map(x->x[1][3],hk))]
#     if l1<l println("!!!!hash",hk)
      if l1<l print("!")
      else
#       print("rels[$k]=",pr(rel),"=>($w)")
        rel=doreplace(rels[w],rel,pos1,pos2,l1)
#       println(pr(rel))
      end
    end
    rel=reduceword(rel)
    if rel!=rels[k] push!(altered,k) end
    if k<=j && length(rel) in lmin:lmax
      for u in 1:length(rel) push!(hh,(hashword(rel,u,l,T),u,k)) end
      for u in 1:length(rel) push!(hh,(hashword(rel,-u,l,T),-u,k)) end
      sort!(hh,by=first)
    end
    rels[k]=rel
    T.total-=T.lengths[k]
    T.lengths[k]=length(rel)
    T.total+=T.lengths[k]
  end
# if !isempty(altered) 
#   @show altered
# end
  length(altered)
end

"""
`Search(P)`

`Search`  performs Tietze transformations on a presentation `P`. It tries
to  reduce the relator lengths by  substituting common subwords of relators
by shorter words.

The  idea is  to find  pairs of  relators `r₁`  and `r₂` of length `l₁` and
`l₂`,  respectively,  such  that  `l₁  ≤  l₂`  and  `r₁`  and `r₂` coincide
(possibly  after  inverting  or  conjugating  one  of them) in some maximal
subword  `w`, say,  of length  greater than  `l₁/2`, and then to substitute
each copy of `w` in `r₂` by the inverse complement of `w` in `r₁`.

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
implemented  as a C routine in the {GAP} kernel. For the given relator `r₁`
of  length  `l₁`,  say,  it  first  determines  the  *minimal match length*
`l=⌊l₁/2⌋+1`.  Then it builds up a hash list for all subwords of length `l`
occurring in the conjugates of `r₁` or `r₁^{-1}`, and finally it loops over
all  long relators `r₂` and  compares the hash values  of their subwords of
length  `l` against this list. A comparison  of subwords which is much more
expensive is only done if a hash match has been found.

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
  if P.printLevel>=3 print("#I  searching subwords\n") end
  Check(P)
  T=P.tietze
  rels=T.relators
  lengths=T.lengths
  flags=T.flags
  T.modified=false
  save=P.saveLimit/100
  if T.total<=0 return end
  while true
    sort!(P)
    numrels=T.numrels
    modified=false
    oldtotal=T.total
    flag=0
    i=1
    while i<numrels
      if flags[i]<=1 && lengths[i]>0
        leng=lengths[i]
        lmax=iseven(leng) ? leng+1 : leng
        if flag<flags[i] flag=flags[i] end
        simultan=1
        j=i
        lastj=0
        k=i+1
        while k<=numrels && lengths[k]<=lmax && simultan<P.searchSimultaneous
          if flags[k]<=1 && (lengths[k]==leng || lengths[k]==lmax)
            lastj=j
            j=k
            simultan+=1
          end
          k+=1
        end
        while k<=numrels && 
              (lengths[k]<leng || flags[k]>1 || flag==0 && flags[k]==0)
          k+=1
        end
        if k>numrels j=lastj end
        if i<=j
          modified=modified || Presentations.SearchC(T,i,j)>0
          i=j
        end
      end
      i+=1
    end
    for i in 1:numrels
      if flags[i]==1 || flags[i]==2 flags[i]=flags[i]-1 end
    end
    if modified
      if T.total<oldtotal
        T.modified=true
        HandleLength1Or2Relators(P)
        sort!(P)
        if P.printLevel>=2 PrintStatus(P;norepeat= true) end
      end
    end
    if T.total>=oldtotal || T.total<=0 || (oldtotal-T.total)/oldtotal<save
       break
    end
  end
end

"""
`SearchEqual( P )`

`SearchEqual`  performs Tietze transformations on a  presentation  `P`.
It tries to alter relators by substituting common subwords of relators by
subwords of equal length.

The  idea is  to find  pairs of  relators `r₁`  and `r₂` of length `l₁` and
`l₂`,  respectively, such that `l₁`  is even, `l₁ ≤  l₂`, and `r₁` and `r₂`
coincide  (possibly after  inverting or  conjugating one  of them)  in some
maximal  subword `w`, say, of length at least `l₁/2`. Let `l` be the length
of `w`. Then, if `l>l₁/2`, the pair is handled as in `Search`. Otherwise,
if `l = l₁/2`, then `SearchEqual` substitutes each copy of `w` in `r₂` by
the inverse complement of `w` in `r₁`.

The  Tietze  option   parameter   `P.searchSimultaneous`  is   used  by
`SearchEqual` in the same way as described for `Search`.

However, `SearchEqual` does  not use the  parameter `P.saveLimit`:
The loop over the relators is executed exactly once.
"""
function SearchEqual(T::Presentation)
  if T.printLevel>=3 print("#I  searching subwords of equal length\n") end
  Check(T)
  tietze=T.tietze
  simultanlimit=T.searchSimultaneous
  sort!(T)
  rels=tietze.relators
  lengths=tietze.lengths
  numrels=tietze.numrels
  modified=false
  oldtotal=tietze.total
  i=1
  while i<numrels
    leng=lengths[i]
    if leng>3 && iseven(leng)
      simultan=1
      j=i
      lastj=0
      k=i+1
      while k <= numrels && (lengths[k] <= leng && simultan<simultanlimit)
          if lengths[k]==leng
              lastj=j
              j=k
              simultan=simultan+1
          end
          k+=1
      end
      while k<=numrels && lengths[k]<leng k+=1 end
      if k>numrels j=lastj end
      if i<=j
        altered=SearchC(tietze, i, j, true)
        modified=modified || altered>0
        i=j
      end
    end
    i+=1
  end
  if modified
    if tietze.total<oldtotal
      tietze.modified=true
      HandleLength1Or2Relators(T)
      sort!(T)
      if T.printLevel>=2 PrintStatus(T;norepeat= true) end
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
use the `PrintPairs` command described below to do this.

As an example we handle a subgroup of index 266 in the Janko group `J₁`.

|    gap> F2 := FreeGroup( "a", "b" );;
    gap> J1 := F2 / [ F2.1^2, F2.2^3, (F2.1*F2.2)^7,
         Comm(F2.1,F2.2)^10, Comm(F2.1,F2.2^-1*(F2.1*F2.2)^2)^6 ];;
    gap> a := J1.1;;  b := J1.2;;
    gap> H := Subgroup ( J1, [ a, b^(a*b*(a*b^-1)^2) ] );;
    gap> P := PresentationSubgroup( J1, H );
    << presentation with 23 gens and 82 rels of total length 530 >>
    gap> GoGo( P );
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
    gap> PrintPairs( P );
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
    gap> PrintPairs( P );
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
  Check(P)
  if !(0<=elim<=2) error("last argument must be 0, 1, or 2") end
  T=P.tietze
  numgens=T.numgens
  if numgens<2 return end
  printlevel=P.printLevel
  pairs=MostFrequentPairs(P, n)
  if length(pairs)<n
    if printlevel>=1 println("#I  Substitute: nbpairs is out of range") end
    return
  end
  i=pairs[n][2]
  j=pairs[n][3]
  k=pairs[n][4]
  if k>1 i=T[-i] end
  if isodd(k) j=T[-j] end
  pair=[i, j]
  gen=NewGenerator(P)
  word=AbsWord(pair,T.generators)
  if P.printLevel >= 1
      print("#I  substituting new generator ", gen, " defined by ", word, "\n")
  end
  AddRelator(P,inv(gen)*word)
  if haskey(P, :imagesOldGens) UpdateGeneratorImages(P, 0, pair) end
  P.printLevel=0;Search(P);P.printLevel=printlevel
  if printlevel==1 P.printLevel=2 end
  if elim>0 EliminateGen(P, abs(pair[elim]))
  else EliminateGen1(P)
  end
  P.printLevel=printlevel
  if T.modified Search(P) end
  if T.numredunds>0 RemoveGenerators(P) end
  if P.printLevel >= 1 PrintStatus(P;norepeat= true) end
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
  Check(P)
  if P.printLevel>=3 print("#I  substituting cyclic joins\n") end
  exponents=GeneratorExponents(P)
  T=P.tietze
  T.modified=false
  if Sum(exponents)==0 return end
  sort!(P)
  gens=T.generators
  numgens=T.numgens
  rels=T.relators
  numrels=T.numrels
  lengths=T.lengths
  printlevel=P.printLevel
  i=1
  while i <= numrels && lengths[i] <= 4
    rel=rels[i]
    if lengths[i]==4 && (rel[1]==T[-rel[3]] && rel[2]==T[-rel[4]])
      num1=abs(rel[1])
      exp1=exponents[num1]
      num2=abs(rel[2])
      exp2=exponents[num2]
      if exp1>0 && (exp2>0 && gcd(exp1, exp2)==1)
        gen=NewGenerator(P)
        numgens=T.numgens
        if printlevel >= 1
          print("#I  substituting new generator ", gen, " defined by ",
                gens[num1] * gens[num2], "\n")
        end
        AddRelator(P, gens[num1]*gen^-exp2)
        AddRelator(P, gens[num2]*gen^-exp1)
        if haskey(P, :imagesOldGens)
            UpdateGeneratorImages(P, 0, [num1, num2])
        end
        gen2=gens[num2]
        if printlevel==1 P.printLevel=2 end
        Eliminate(P, gens[num1])
        num2=findfirst(==(gen2),gens)
        Eliminate(P, gens[num2])
        P.printLevel=printlevel
        numgens=T.numgens
        numrels=T.numrels
        T.modified=true
        i=0
      end
    end
    i+=1
  end
  if T.modified
    HandleLength1Or2Relators(P)
    sort!(P)
    if printlevel>=1 PrintStatus(P;norepeat= true) end
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
applied to `P` (see  the function `InitGeneratorImages` below), whereas
a call  of the  function `AddGenerator` (which   does not perform  Tietze
transformations) will delete these lists and hence terminate the tracing.

Example.

|    gap> G := PerfectGroup( 960, 1 );
    PerfectGroup(960,1)
    gap> P := Presentation( G );
    << presentation with 6 gens and 21 rels of total length 84 >>
    gap> P.generators;
    [ a, b, s, t, u, v ]
    gap> GoGo( P );
    &I  there are 3 generators and 10 relators of total length 81
    &I  there are 3 generators and 10 relators of total length 80
    gap> PrintGenerators( P );
    &I  1.  a   31 occurrences   involution
    &I  2.  b   26 occurrences
    &I  3.  t   23 occurrences   involution
    gap> a := P.generators[1];;
    gap> b := P.generators[2];;
    gap> Substitute( P, a*b, "ab" );
    &I  substituting new generator ab defined by a*b
    &I  there are 4 generators and 11 relators of total length 83
    gap> Go(P);
    &I  there are 3 generators and 10 relators of total length 74
    gap> PrintGenerators( P );
    &I  1.  a   23 occurrences   involution
    &I  2.  t   23 occurrences   involution
    &I  3.  ab   28 occurrences |
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
  if narg==1
    gen=AbstractGenerator(arg[1])
    AddGenerator(P, gen)
  else
    AddGenerator(P)
    gen=(P.generators)[length(P.generators)]
  end
  if P.printLevel>=1
    print("#I  substituting new generator ", gen, " defined by ", word, "\n")
  end
  AddRelator(P,inv(gen)*word)
  if images isa Vector
    P.imagesOldGens=images
    UpdateGeneratorImages(P, 0, tzword)
  end
  if P.printLevel>=1 PrintStatus(P;norepeat= true) end
end

"""
UpdateGeneratorImages( T, n, word )  . . . . update the generator images
                                            after a Tietze transformation

`UpdateGeneratorImages`  assumes  that  it  is  called  by  a function that
performs  Tietze transformations to  a presentation `T`  in which images of
the  old generators are being traced as  Tietze words in the new generators
as  well as  preimages of  the new  generators as  Tietze words  in the old
generators.

If  `n` is zero, it assumes that a new generator defined by the Tietze word
`word`  has just been added to the  presentation. It converts `word` from a
Tietze  word in the new  generators to a Tietze  word in the old generators
and adds that word to the list of preimages.

If  `n` is greater than zero, it assumes that the `n`-th generator has just
been  eliminated from  the presentation.  It updates  the images of the old
generators  by replacing  each occurrence  of the  `n`-th generator  by the
given Tietze word `word`.

Note:  `UpdateGeneratorImages` is  considered to  be an  internal function.
Hence it does not check the arguments.
"""
function UpdateGeneratorImages(T::Presentation, n, word)
  if n==0
    newim=Int[]
    for num in word
      if num>0 append!(newim, T.preImagesNewGens[num])
      else     append!(newim, -reverse(T.preImagesNewGens[-num]))
      end
    end
    push!(T.preImagesNewGens, reduceword(newim))
  elseif n>0
    invword=-reverse(word)
    oldnumgens=length(T.imagesOldGens)
    for i in 1:oldnumgens
      image=T.imagesOldGens[i]
      newim=[]
      for j=1:length(image)
        if image[j]==n append!(newim, word)
        elseif image[j]==-n append!(newim, invword)
        else push!(newim, image[j])
        end
      end
      T.imagesOldGens[i]=reduceword(newim)
    end
  end
end

"""
Terminates   the  tracing  of  generator   images,  i.e.,  it  deletes  the
corresponding components of `T`.
"""
function StopTracingGeneratorImages(T::Presentation)
  delete!(T, :imagesOldGens)
  delete!(T, :preImagesNewGens)
  if haskey(T, :oldGenerators) delete!(T, :oldGenerators) end
  if T.printLevel>=1 println("#I  terminated the tracing of generator images")
  end
end

"""
RenumberGenerators( presentation, sort list ) . . . . Renumber the
             generators of T according to L

`RenumberGenerators` renumbers the generators of the given presentation
T according to the given sort list L.
"""
function RenumberGenerators(P::Presentation, L)
  if P.printLevel>=3 println("#I  renumbering the Tietze generators") end
  RemoveGenerators(P)
  T=P.tietze
  numgens=T.numgens
  if !(L isa Vector) error("second argument must be a list") end
  if length(L) != numgens || gapSet(L) != 1:numgens
      error("<L> must determine a permutation of the generator numbers")
  end
  perm=PermList(L) ^ -1
  invsort=ListPerm(perm)
  comps=P.components
  gens=T.generators
  rels=T.relators
  numrels=T.numrels
  oldcomps=deepcopy(comps)
  oldgens=deepcopy(gens)
  oldinvs=copy(T.inverses)
  new=copy(T.inverses)
  println("gens==",gens,"\ncomps==",comps,"\ninvs==",T.inverses,"\nrels==",rels)
  for i in 1:numgens
    j=L[i]
    gens[i]=oldgens[j]
    comps[i]=oldcomps[j]
    if oldinvs[numgens+1+j]==j T[-i]=i else T[-i]=-i end
    new[numgens+1+i]=invsort[i]
    new[numgens+1-i]=-invsort[i]
  end
  println("gens==",gens,"\ncomps==",comps,"\ninvs==",T.inverses,"\nnew==",new)
  for i in 1:numrels
      rels[i]=map(j->new[numgens+1+j], rels[i])
  end
  print("rels==", rels, "\n")
  if haskey(P, :imagesOldGens)
    println("preimages==", P.preImagesNewGens)
    println("images==", P.imagesOldGens)
    for i in 1:length(P.imagesOldGens)
      P.imagesOldGens[i]=map(j->new[numgens+1+j],P.imagesOldGens[i])
    end
    P.preImagesNewGens=Permuted(P.preImagesNewGens, perm)
    println("preimages==", P.preImagesNewGens)
    println("images==", P.imagesOldGens)
  end
end

"""
`Shrink(p [,tries])`

This is our  own program to simplify group presentations.  We have found
heuristics which make it somewhat  more efficient than GAP' s programs
`simplify` and  `GoGo`, but  the algorithm depends  on random
numbers so  is not  reproducible. The  main idea  is to  rotate relators
between calls  to GAP  functions. By default  1000 such  rotations are
tried (unless the  presentation is so small that  less rotations exhaust
all possible  ones), but the  actual number  tried can be  controlled by
giving a second  parameter `tries` to the function.  Another useful tool
to deal  with presentations  is `tryconjugate`  described in
the utility functions.

|    gap> display_balanced(p);
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
    19: CbdbadcDbbdCbDDadcBCDAdBCDbdaDCDbdcbadcBCDAdBCDBBdacDbdccb
        =abdbcabdcAdcbCDDBCDABDABDbbdcbDadcbCDAdBCabDACbdBadcaDbAdd
    gap> Shrink(p);   
    #I  there are 4 generators and 19 relators of total length 332
    #I  there are 4 generators and 17 relators of total length 300
    #I  there are 4 generators and 17 relators of total length 282
    #I  there are 4 generators and 17 relators of total length 278
    #I  there are 4 generators and 16 relators of total length 254
    #I  there are 4 generators and 15 relators of total length 250
    #I  there are 4 generators and 15 relators of total length 248
    #I  there are 4 generators and 15 relators of total length 246
    #I  there are 4 generators and 14 relators of total length 216
    #I  there are 4 generators and 13 relators of total length 210
    #I  there are 4 generators and 13 relators of total length 202
    #I  there are 4 generators and 13 relators of total length 194
    #I  there are 4 generators and 12 relators of total length 174
    #I  there are 4 generators and 12 relators of total length 170
    #I  there are 4 generators and 12 relators of total length 164
    #I  there are 4 generators and 12 relators of total length 162
    #I  there are 4 generators and 12 relators of total length 148
    #I  there are 4 generators and 12 relators of total length 134
    #I  there are 4 generators and 12 relators of total length 130
    #I  there are 4 generators and 12 relators of total length 126
    #I  there are 4 generators and 12 relators of total length 124
    #I  there are 4 generators and 12 relators of total length 118
    #I  there are 4 generators and 12 relators of total length 116
    #I  there are 4 generators and 11 relators of total length 100
    gap> display_balanced(p);
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
    11: BcccbdcAb=dcbACdddc|
"""
function Shrink(g::Presentation,lim=1000)
  T=g.tietze
  if isempty(T.relators) return end
  rot(i)=T.relators[i]=circshift(T.relators[i],-1)
  function test()
    if prod(filter(!iszero,T.lengths))<lim
      v=fill(0,T.numrels)
      if length(v)==0 return false end
      while true
        before=[T.numrels,T.total]
        j=length(v)
        while v[j]==T.lengths[j]-1
          rot(j)
          v[j]=0
          j-=1
          if j==0 return false end
        end
        rot(j)
        v[j]+=1
        GoGo(g)
        if [T.numrels,T.total]<before return true end
      end
    else
      for i in 1:lim
        i=rand(1:T.total)
        j=1
        while i>T.lengths[j]
          i-=T.lengths[j]
          j+=1
        end
        rot(j)
        before=[T.numrels,T.total]
        GoGo(g)
        if [T.numrels,T.total]<before return true end
      end
    end
    return false
  end
  while test() end
end

function invertcase(s::String)
  String(map(collect(s))do x
    if isuppercase(x) lowercase(x) else uppercase(x) end
  end)
end

function display_balanced(g,dumb=false)
  f(i,w)=(w*w)[i:i+div(length(w),2)-1]
  minusc="abcdefghijklmnopqrstuvwxyz"
  used=Set{Int}()
  l=map(g.tietze.relators)do x
    map(x)do y
      push!(used, abs(y))
      y<0 ? uppercase(minusc[-y]) : minusc[y]
    end 
  end
  l=String.(l)
  if g.tietze.numgens>length(used)
    print("There are ",g.tietze.numgens-length(used)," free generators\n")
  end
  for i in 1:length(l)
    w=l[i]
    lw=length(w)
    if isodd(lw) || dumb print(i, ": ", w, "=1\n")
    elseif lw>0
      m=argmax(map(i->count(islowercase, f(i,w)), 1:lw))
      print(i,": ",f(m,w),"=",invertcase(reverse(f(m+div(lw,2),w))), "\n")
    end
  end
end

"""
`conjugate(p,conjugation)`

This program modifies a presentation by conjugating a generator by another.
The  conjugation to  apply is  described by  a length-3  string of the same
style  as  the  result  of  `display_balanced`,  that is `"abA"` means
replace  the second generator by its  conjugate by the first, and  `"Aba"`
means replace it by its conjugate by the inverse of the first.

```julia-rep1
julia> Presentations.display_balanced(P)
1: dabcd=abcda
2: dabcdb=cabcda
3: bcdabcd=dabcdbc

julia> Presentations.display_balanced(conjugate(P,"Cdc"))
<< presentation with 4 generators, 3 relators of total length 36>>
1: dcabdc=cabdca
2: abdcab=cabdca
3: bdcabd=cabdca
```
"""
function conjugate(p::Presentation, s)
  minmaj="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
  l=map(l->findfirst(==(l),minmaj),collect(s))
  l=map(x->x>26 ? 26-x : x,l)
  if length(l)!=3 || l[1]!=-l[3] error("should be conjugate like abA") end
  p=deepcopy(p);SubstituteGen(p.tietze,l[2],l)
# Shrink(p, 100)
  GoGo(p)
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

|    gap> display_balanced(p);
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
    gap> p:=tryconjugate(p); 
    #I  there are 4 generators and 11 relators of total length 100
    #I  there are 4 generators and 11 relators of total length 120
    #I  there are 4 generators and 10 relators of total length 100
    #I  there are 4 generators and 11 relators of total length 132
    #I  there are 4 generators and 11 relators of total length 114
    #I  there are 4 generators and 11 relators of total length 110
    #I  there are 4 generators and 11 relators of total length 104
    #I  there are 4 generators and 11 relators of total length 114
    #I  there are 4 generators and 11 relators of total length 110
    #I  there are 4 generators and 11 relators of total length 104
    #I  there are 4 generators and 8 relators of total length 76
    #I  there are 4 generators and 8 relators of total length 74
    #I  there are 4 generators and 8 relators of total length 72
    #I  there are 4 generators and 8 relators of total length 70
    #I  there are 4 generators and 7 relators of total length 52
    # d->adA gives length 52
    << presentation with 4 gens and 7 rels of total length 52 >>
    gap> display_balanced(p); 
    1: ba=ab
    2: dc=cd
    3: aca=cac
    4: dbd=bdb
    5: bcb=cbc
    6: adad=dada
    7: aBcADbdac=dBCacbdaB
    gap> tryconjugate(p,48);
    #I  there are 4 generators and 7 relators of total length 54
    #I  there are 4 generators and 7 relators of total length 54
    #I  there are 4 generators and 7 relators of total length 60
    #I  there are 4 generators and 7 relators of total length 60
    #I  there are 4 generators and 7 relators of total length 48
    # d->bdB gives length 48
    << presentation with 4 gens and 7 rels of total length 48 >>
    gap> display_balanced(last);
    1: ba=ab
    2: bcb=cbc
    3: cac=aca
    4: dbd=bdb
    5: cdc=dcd
    6: adad=dada
    7: dAbcBa=bAcBad|
"""
function tryconjugate(p::Presentation,tp=[0,0];info=[0,0])
  perf(p)=[p.tietze.numrels, p.tietze.total]
  function triples(p) local res, v, i, r
    res=Vector{Int}[]
    for v in p.tietze.relators, i in 1:length(v)-2
      if v[i]==-v[i+2]
        if v[i+1]>0 push!(res,v[i:i+2])
        else        push!(res,-reverse(v[i:i+2]))
        end
      end
    end
    first.(sort(tally(res),by=x->-x[2]))
  end
  function expand(p)local res
    res=map(x->[x, conjugate(p, x)], triples(p))
    sort!(res,by=x->perf(x[2]))
  end
  p=deepcopy(p)
  GoGo(p)
  if iszero(tp) tp=perf(p) end
  p1=[]
  for c in triples(p)
    print(pr(c),"=> ")
    n=deepcopy(p);SubstituteGen(n.tietze,c[2],[-c[1],c[2],-c[3]]);GoGo(n)
    if perf(n)<=info
      println("# ",pr(c)," gives rels/len ", perf(n))
      display_balanced(n)
    end
    if perf(n)<tp
      println("# ",pr(c)," gives ",n)
      return n
    end
    push!(p1,(c,n))
  end
  sort!(p1,by=x->perf(x[2]))
  for p2 in p1, c in triples(p2[2])
    n=deepcopy(p2[2]);SubstituteGen(n.tietze,c[2],[-c[1],c[2],-c[3]]);GoGo(n)
    print(pr(p2[1]),"+",pr(c),"=>")
    if perf(n)<=info
      println("# ",pr(p2[1]),"->",pr(c)," gives rels/len ", perf(n))
      display_balanced(n)
    end
    if perf(n)<tp
      println("# ",pr(p2[1]),"->",pr(c)," gives ",n)
      return n
    end
  end
  println("\n# could not shrink ",p)
  p
end

end
