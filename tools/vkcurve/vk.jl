"""
This is a port to Julia of the GAP3 package VKcurve written by David Bessis
and Jean Michel in 2002.

The  main function of the VKcurve package computes the fundamental group of
the   complement  of   a  complex   algebraic  curve   in  `ℂ²`,  using  an
implementation  of the Van  Kampen method (see  for example

D. Cheniot. "Une démonstration du théorème de   Zariski sur les sections hyperplanes d'une hypersurface projective et du   théorème de Van Kampen sur le groupe fondamental du   complémentaire d'une courbe projective plane."
Compositio Math., 27:141--158, 1973.

for  a clear and modernized account of this method). Here is an example for
curves given by the zeroes of two-variable polynomials in `x` and `y`.

```julia-repl
julia> using VKcurve

julia> @Mvp x,y

julia> fundamental_group(x^2-y^3)
Presentation: 2 generators, 1 relator, total length 6
1: bab=aba

julia> fundamental_group((x+y)*(x-y)*(x+2*y))
Presentation: 3 generators, 2 relators, total length 12
1: abc=bca
2: cab=abc
```
Here  we define the variables and then  give the curves as argument. Though
approximate  calculations are used  at various places,  they are controlled
and  the final result is exact  (technically speaking, the computations use
`Rational{BigInt}` or `Complex{Rational{BigInt}}` since the precision given
by  floats  in  unsufficient.  It  might  be  possible  to use intervals of
bigfloats  to make faster  computations, but it  would make the programming
more difficult).

The output is a  `struct` which  contains lots  of information  about the
computation, including a presentation of the computed fundamental group,
which is what is displayed by default when printing it.

Our  motivation  for  writing  this  package  in  2002 was to find explicit
presentations  for  generalized  braid  groups  attached to certain complex
reflection  groups. Though presentations  were known for  almost all cases,
six  exceptional cases were missing (in the notations of Shephard and Todd,
these  cases are `G₂₄`, `G₂₇`, `G₂₉`, `G₃₁`,  `G₃₃` and `G₃₄`). Since the a
priori  existence  of  nice  presentations  for  braid groups was proved in

D. Bessis. "Zariski theorems and diagrams for braid groups.",
Invent. Math. 145:487--507, 2001

it  was upsetting not to  know them explicitly. In  the absence of any good
grip on the geometry of these six examples, brute force was a way to get an
answer. Using VKcurve, we have obtained presentations for all of them (they
have been confirmed by less computational methods since).

If  you  are  not  interested  in  the  details  of  the  algorithm, and if
'fundamental_group' as in the above examples gives you satisfactory answers
in a reasonable time, then you do not need to read this manual any further.

To  implement the algorithms, we needed  to write auxiliary facilities, for
instance  find  `Complex{Rational}`  approximations  of  zeros  of  complex
polynomials,  or work with piecewise linear  braids, which may be useful on
their own. These various facilities are documented in this manual.

Before  discussing  our  actual  implementation,  let  us  give an informal
summary  of the mathematical  background. Our strategy  is adapted from the
one  originally described in the 1930's by Van Kampen. Let `C` be an affine
algebraic  curve, given as the set of  zeros in `ℂ^2` of a non-zero reduced
polynomial  `P(x,y)`.  The  problem  is  to  compute  a presentation of the
fundamental  group of  `ℂ^2-C`. Consider  `P` as  a polynomial in `x`, with
coefficients  in the ring  of polynomials in  `y`, that is `P=α₀(y)xⁿ+α₁(y)
xⁿ⁻¹+…+αₙ₋₁(y)x+αₙ(y)`, where the `αᵢ∈ℂ[y]`. Let `Δ(y)` be the discriminant
of  `P` or, in other words, the resultant  of `P` and `∂P/∂x`. Since `P` is
reduced,  `Δ` is non-zero. Let `y₁,…,y_d` be the roots of the corresponding
reduced polynomial `Δ_{red}`. For a generic value of `y`, the polynomial in
`x`  given by  `P(x,y)` has  `n` distinct  roots. When  `y=yⱼ`, with `j` in
`1,…,d`,  we  are  in  exactly  one  of  the  following  situations: either
`P(x,yⱼ)=0`  (we then say that  `yⱼ` is bad), or  `P(x,yⱼ)` has a number of
roots  in  `x`  strictly  smaller  than  `n`.  Fix  `y₀` in `ℂ-{y₁,…,y_d}`.
Consider  the projection `p:  ℂ²→ ℂ, (x,y)↦  y`. It restricts  to a locally
trivial  fibration with base space `B=ℂ-{y₁,…,y_d}` and fibers homeomorphic
to  the complex plane with  `n` points removed. We  denote by `E` the total
space `p⁻¹(B)` and by `F` the fiber over `y₀`. The fundamental group of `F`
is  isomorphic to the free group on `n` generators. Let `γ₁,…,γ_d` be loops
in  the  pointed  space  `(B,y₀)`  representing  a  generating  system  for
`π₁(B,y₀)`.  By trivializing  the pullback  of `p`  along `γᵢ`,  one gets a
(well-defined  up to  isotopy) homeomorphism  of `F`,  and a (well-defined)
automorphism `φᵢ` of the fundamental group of `F`, identified with the free
group `Fₙ` by the choice of a generating system `f₁,…,fₙ`. An effective way
of  computing `φᵢ` is by following the solutions in `x` of `P(x,y)=0`, when
`y`  moves along `φᵢ`. This defines a loop in the space of configuration of
`n`  points in a plane, hence an element  `bᵢ` of the braid group `Bₙ` (via
an  identification of `Bₙ` with the fundamental group of this configuration
space).  Let `φ` be the Hurwitz action of  `Bₙ` on `Fₙ`. All choices can be
made in such a way that `φᵢ=φ(bᵢ)`. The theorem of Van Kampen asserts that,
if  there are  no bad  roots of  the discriminant,  a presentation  for the
fundamental  group of  `ℂ^2-C` is  `⟨f₁,…,fₙ∣∀i,j,φᵢ(fⱼ)=fⱼ⟩`. A variant of
the  above presentation  (see 'VKquotient')  can be  used to  deal with bad
roots of the discriminant.

This algorithm is implemented in the following way.

  - As input, we have a polynomial `P`. The polynomial is reduced if it was
    not.

  - The discriminant `Δ` of `P`  with respect to `x` is computed.  It is a
    polynomial in `y`.

  - The  roots  of  `Δ`  are  approximated,  via  the  following procedure.
    First,  we reduce `Δ` and get  `Δ_{red}` (generating the radical of the
    ideal  generated  by  `Δ`).  The  roots  `{y₁,…,y_d}`  of `Δ_{red}` are
    separated by 'separate_roots' (which uses Newton's method).

  - Loops around  these roots  are computed  by 'loops_around_punctures'.
    This  function first computes  some sort of  honeycomb, consisting of a
    set  `S` of  affine segments,  isolating the  `yᵢ`. Since  it makes the
    computation  of the monodromy  more effective, each  inner segment is a
    fragment of the mediatrix of two roots of `Δ`. Then a vertex of one the
    segments  is chosen as a basepoint, and  the function returns a list of
    lists  of  oriented  segments  in  `S`:  each list of segment encodes a
    piecewise linear loop `γᵢ` circling one of `yᵢ`.

  - For each  segment in  `S`, we  compute the  monodromy braid  obtained by
    following  the solutions in `x` of  `P(x,y)=0` when `y` moves along the
    segment. By default, this monodromy braid is computed by
    `follow_monodromy`. The strategy is to compute a piecewise-linear braid
    approximating  the actual monodromy geometric braid. The approximations
    are controlled. The piecewise-linear braid is constructed step-by-step,
    by  computations of linear pieces. As soon as new piece is constructed,
    it  is converted  into an  element of  `Bₙ` and  multiplied; therefore,
    though  the braid may consist of a  huge number of pieces, the function
    `follow_monodromy`  works  with  constant  memory.  The  packages  also
    contains  a  variant  function  `approx_follow_monodromy`,  which  runs
    faster, but without guarantee on the result (see below).

  - The monodromy braids `bᵢ` corresponding to the loops `γᵢ` are obtained
    by  multiplying  the  corresponding  monodromy  braids of segments. The
    action  of these elements of `Bₙ` on the free group `Fₙ` is computed by
    'hurwitz'  and the resulting  presentation of the  fundamental group is
    computed  by 'VKquotient'. It happens for  some large problems that the
    whole  fundamental group  process fails  here, because  the braids `bᵢ`
    obtained  are  too  long  and  the  computation  of  the action on `Fₙ`
    requires thus too much memory. We have been able to solve such problems
    when  they occur  by calling  on the  `bᵢ` at  this stage  our function
    'shrink'  which  finds  smaller  generators  for  the  subgroup of `Bₙ`
    generated  by the `bᵢ`  (see the description  in `Gapjm.Garside`). This
    function is called automatically at this stage if 'VK.shrinkBraid'
    is set to 'true' (the default for this variable is 'false').

  - Finally, the  presentation is  simplified by  'simplify'. This function
    is  a heuristic  adaptation and  refinement of  the basic functions for
    simplifying presentations. It is non-deterministic.

From  the algorithmic point of view, memory should not be an issue, but the
procedure  may  take  a  lot  of  CPU  time  (the  critical  part being the
computation  of the monodromy braids  by 'follow_monodromy'). For instance,
an  empirical study with  the curves `x²-yⁿ`  suggests that the needed time
grows exponentially with `n`. Two solutions are offered to deal with curves
for which the computation time becomes unreasonable.

A  global  variable  `VK.approx_monodromy`  controls  which  monodromy
function  is used.  The default  value of  this variable  is `false`, which
means  that `follow_monodromy` will be used. If  the variable is set by the
user  to `true`  then the  function `approx_follow_monodromy`  will be used
instead.  This  function  runs  faster  than  `follow_monodromy`,  but  the
approximations  are no longer  controlled. Therefore presentations obtained
while  `VK.approx_monodromy`  is  set  to  'true'  are  not certified.
However,  though  it  is  likely  that  there  exists  examples  for  which
`approx_follow_monodromy` actually returns incorrect answers, we still have
not seen one.
"""
#module VKcurve
using Gapjm
"""
`nearest_pair(v::Vector{<:Complex})`

returns  a pair whose first element is the minimum distance (in the complex
plane)  between two elements  of `v`, and  the second is  a pair of indices
`[i,j]` such that `v[i],v[j]` achieves this minimum.

julia> nearest_pair([1+im,0,1])
1=>[1,3]
"""
function nearest_pair(v)
  l=combinations(eachindex(v),2)
  m,c=findmin(((x,y),)->abs(v[x]-v[y]),l)
  m=>l[c]
end

"`dist_seg(z,a,b)` distance (in the complex plane) of `z` to segment `[a,b]` "
function dist_seg(z,a,b)
  b-=a
  z-=a
  r=abs(b)
  z*=r/b
  rz,iz=reim(z)
  rz<0 ? abs(z) : rz>r ? abs(z-r) : iz>0 ? iz : -iz
end

#---------------------- global functions --------------------------
@GapObj struct VKopt end
const VK=VKopt(Dict(
:name=>"vkcurve",
:version=>"2.0",
:date=>[2009,3],
:homepage=>"http://webusers.imj-prg.fr/~jean.michel/vkcurve.html",
:copyright=>
"(C) David Bessis, Jean Michel -- compute Pi₁ of hypersurface complements",
:approx_monodromy=>false,
:showallnewton=>false,
:NewtonLim=>800,
:AdaptivityFactor=>10,
:shrinkBraid=>false,
:mvp=>1))

@GapObj struct VKrec{T}
  curve::Mvp{T,Int}
  ismonic::Bool
end

function Base.show(io::IO,r::VKrec)
  if haskey(r,:presentation) display_balanced(r.presentation)
  else xdisplay(r.prop)
  end
end

# Loops(r)
#  r should be a record with the fields
#    .roots   -- roots of the curve discriminant
#    .ismonic -- tells if the discriminant is monic in x
#  the function computes the following fields describing loops around the
#  .roots around from a basepoint:
#
#  .loops --  a list of loops, each described by a list of indices in r.segments
#             (a negative index tells to follow the reverse segment)
#  .segments -- oriented segments represented as a pair of indices in r.points.
#  .points -- a list of points stored a complex decimal numbers.
#  .basepoint -- holds the chosen basepoint
#
function Loops(r)
# here we  have loops around the  'true' roots and around  the 'extra'
# roots setdiff(r.roots,r.trueroots). We get rid of the extra loops
# and the associated segments and points, first saving the basepoint.
# (its location is known now, and maybe not later? J.M.)
  p=r.loops[1][1]
  r.basepoint=p<0 ? last(r.segments[-p]) : first(r.segments[p])
  if !r.ismonic
    r.loops=r.loops[sort(indexin(r.trueroots,r.roots))]
    segmentNumbers=sort(union(map(x->abs.(x),r.loops)...))
    r.segments=r.segments[segmentNumbers]
    r.loops=map(x->map(y->y<0 ? -findfirst(==(-y),segmentNumbers) :
                    findfirst(==(y),segmentNumbers), x),r.loops)
    uniquePoints=sort(union(r.segments...))
    r.points=r.points[uniquePoints]
    r.segments=map(x->indexin(x,uniquePoints), r.segments)
    r.basepoint=findfirst(==(r.basepoint),uniquePoints)
  end
  if VK.showSegments
    println("# There are ",length(r.segments)," segments in ",
            length(r.loops)," loops")
  end
  if VK.showWorst
    l=map(enumerate(r.segments))do (i,s)
      m,ixm=findmin(dist_seg.(r.roots, r.points[s[1]], r.points[s[2]]))
      (m,i,ixm)
    end
    sort!(l)
    print("worst segments:\n")
    for i in 1:min(5,length(l))
      d,s,s1=l[i]
      println("segment ",s,"==",r.segments[s]," dist to ",s1,"-th root is ",approx(d))
    end
  end
  if length(r.roots)>1
    d,p=nearest_pair(r.roots) # find the minimum distance between two roots
    if VK.showRoots
      println("\nMinimum distance==",d," between roots ",join(p," and ")," of discriminant")
    end
    r.dispersal=d
  else
    r.dispersal=1/1000
  end
# and round points to m/100
end

"""
`Discy(r)`

`r`  should be a record with field `r.curve`, a quadratfrei `Mvp` in `x,y`.
The discriminant of this curve with respect to `x` (a polynomial in `y`) is
computed. First, the curve is split in
  - `r.curveVerticalPart`: gcd(coefficients(r.curve,:x)) (an Mvp in y).
  - `r.nonVerticalPart`:   r.curve/r.curveVerticalPart
Then `r.discy=discriminant(Pol(r.nonVerticalPart,:x))` is computed. 
Its   quadratfrei  part  is  computed,  stripped  of  factors  common  with
`r.curveVerticalPart` and then factored (if possible which for now means it
is  a  polynomial  over  the  rationals),  and  its  factors  are stored in
`r.discyFactored` as a list of `Mvp` in `x`.
"""
function Discy(r)
  r.curveVerticalPart=gcd(Pol.(values(coefficients(r.curve,:x))))
  if VK.showRoots && degree(r.curveVerticalPart)>0
    println("Curve has ",degree(r.curveVerticalPart)," linear factors in y")
  end
  r.nonVerticalPart=exactdiv(r.curve,r.curveVerticalPart(Mvp(:y)))
  # discriminant wrto x of r.nonVerticalPart
  d=Pol(discriminant(Pol(r.nonVerticalPart,:x)))
  if iszero(d)
    error("Discriminant is 0 but ", r.curve," should be quadratfrei")
  end
  if VK.showRoots print("Discriminant has ",degree(d)," roots, ") end
  d=exactdiv(d,gcd(d,derivative(d)))
  if VK.showRoots println(" of which ", degree(d), " are distinct") end
  common=gcd(d,r.curveVerticalPart)
  if VK.showRoots && degree(common)>0
    println(" and of which ",degree(common)," are roots of linear factors")
  end
  d=exactdiv(d,common)
  d//=d[end]
  r.discy=d
  if eltype(coefficients(d))<:Union{Complex{<:Rational},Rational}
    r.discyFactored=factor(r.discy)
    if r.discyFactored isa Tuple r.discyFactored=r.discyFactored[1] end
    r.discyFactored=collect(keys(r.discyFactored)) # no multiplicities
  else
    r.discyFactored=[d]
  end
end

function SearchHorizontal(r) # Searching for a good horizontal
  height=9
  while true
    height+=1
    section=Pol(r.curve(x=height))
    section=exactdiv(section,gcd(derivative(section), section))
    if degree(section)==degree(r.curve,:y) &&
       degree(gcd(r.discy,section))==0 break end
  end
  section=exactdiv(section,gcd(section,r.curveVerticalPart))
  println("Curve is not monic in x -- Trivializing along horizontal line x == ", height)
  r1=VKrec(r.curve*(Mvp(:x)-height),false,Dict{Symbol,Any}()) # is it monic?
  r1.height=height
  r1.input=r.curve
  if haskey(r,:name) r1.name=r.name end
  Discy(r1)
# set trueroots  to roots  of Discy  which do  not correspond only to
# intersections of the curve with the chosen horizontal
  r1.trueroots=r.roots
  r1.verticallines=r.roots[1:degree(r.curveVerticalPart)]
  r1.roots=copy(r1.trueroots)
  append!(r1.roots,separate_roots(section,1000))
  r1
end

function TrivialCase(r)
  r.presentation=Presentation(FpGroup(degree(r.curve,:x)))
  r
end

"""
`fundamental_group(curve::Mvp; printlevel=0)`

`curve` should be an `Mvp` in `x` and `y` representing an equation `f(x,y)`
for  a  curve  in  `ℂ²`.  The  coefficients should be rationals or gaussian
rationals.  The result is  a record with  a certain number  of fields which
record steps in the computation described in this introduction:

```julia-repl
julia> @Mvp x,y; @Pol y

julia> r=fundamental_group(x^2-y^3)
Presentation: 2 generators, 1 relator, total length 6
1: bab=aba

julia> propertynames(r)
(:curve, :ismonic, :prop, :rawPresentation, :B, :basepoint, :dispersal, :monodromy, :discyFactored, :segments, :braids, :roots, :nonVerticalPart, :discy, :zeros, :curveVerticalPart, :points, :loops, :presentation)

julia> r.curve # the given equation
Mvp{Int64}: x²-y³

julia> r.discy # its discriminant wrt x
Pol{Rational{Int64}}: (1//1)y

julia> r.roots  # roots of the discriminant
1-element Vector{Rational{Int64}}:
 0//1

julia> r.points # for points, segments and loops see loops_around_punctures
4-element Vector{Complex{Rational{Int64}}}:
  0//1 - 1//1*im
 -1//1 + 0//1*im
  1//1 + 0//1*im
  0//1 + 1//1*im

julia> r.segments
4-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 3]
 [2, 4]
 [3, 4]

julia> r.loops
1-element Vector{Vector{Int64}}:
 [4, -3, -1, 2]

julia> r.zeros # zeroes of curve(y=pt) when pt runs over r.points
4-element Vector{Vector{Complex{Rational{BigInt}}}}:
 [2378//3363 + 2378//3363*im, -2378//3363 - 2378//3363*im]
 [0//1 + 1//1*im, 0//1 - 1//1*im]
 [1//1 + 0//1*im, -1//1 + 0//1*im]
 [-2378//3363 + 2378//3363*im, 2378//3363 - 2378//3363*im]

julia> r.monodromy # monodromy around each r.segment
4-element Vector{GarsideElt{Perm{Int16}, BraidMonoid{Perm{Int16}, CoxSym{Int16}}}}:
 (Δ)⁻¹
 Δ
 .
 Δ

julia> r.braids # monodromy around each r.loop
1-element Vector{GarsideElt{Perm{Int16}, BraidMonoid{Perm{Int16}, CoxSym{Int16}}}}:
 Δ³

julia> display_balanced(r.presentation) # printing of r by default
1: bab=aba
```
The second optional argument triggers  the display of information on the
progress of the  computation. It is recommended to  set the `printlevel`
at 1 or 2  when the computation seems to take a  long time without doing
anything. `printlevel` set  at 0 is the default and  prints nothing; set
at 1 it shows which segment is  currently active, and set at 2 it traces
the computation inside each segment.
```julia-repl
julia> fundamental_group(x^2-y^3,printlevel=1);
# There are 4 segments in 1 loops
# The following braid was computed by follow_monodromy in 8 steps.
monodromy[1]=B(-1)
#segment 1/4 Time==0.00495sec
# The following braid was computed by follow_monodromy in 8 steps.
monodromy[2]=B(1)
#segment 2/4 Time==0.00481sec
# The following braid was computed by follow_monodromy in 8 steps.
monodromy[3]=B()
#segment 3/4 Time==0.00403sec
# The following braid was computed by follow_monodromy in 8 steps.
monodromy[4]=B(1)
#segment 4/4 Time==0.00392sec
# Computing monodromy braids
loops=[
r.B(1,1,1),
]
Presentation: 2 generators, 1 relator, total length 6
```
"""
function fundamental_group(curve::Mvp;printlevel=0,abort=false)
  VK.showSingularProj=VK.showBraiding=VK.showLoops=
  VK.showAction=VK.showInsideSegments=VK.showWorst=
  VK.showZeros=VK.showNewton=VK.showRoots=printlevel>=2
  VK.showSegments=VK.showgetbraid=printlevel>=1
  if !issubset(variables(curve),[:x,:y])
    error(curve," should be an Mvp in x,y")
  end
  d=gcd(curve, derivative(curve,:x))
  if degree(d,:x)>0
    xprintln("**** Warning: curve is not quadratfrei: dividing by ", d)
    curve=exactdiv(curve,d)
  end
#  record with fields .curve and .ismonic if the curve is monic in x
  r=VKrec(curve,degree(Pol(curve,:x)[end])==0,Dict{Symbol,Any}())
#  we should make coefficients(curve) in  Complex{Rational}
  Discy(r);
  if VK.showRoots println("Computing roots of discriminant...") end
  r.roots=vcat(map(p->separate_roots(p,1000),r.discyFactored)...)
  if degree(r.curveVerticalPart)>0
    prepend!(r.roots,separate_roots(r.curveVerticalPart, 1000))
  end
  if isempty(r.roots) return TrivialCase(r) end
  if !r.ismonic r=SearchHorizontal(r) end
  loops=convert_loops(loops_around_punctures(r.roots))
  merge!(r.prop,pairs(loops))
  Loops(r)
  #------------------compute r.zeros[i]=zeros of r.curve(y=r.points[i])
  if VK.showRoots
    println("Computing zeros of curve at the ",length(r.points)," segment extremities...")
  end
  mins=Tuple{Float64,Int}[]
  r.zeros=map(1:length(r.points))do i
    if VK.showZeros print("<",i,"/",length(r.points),">") end
    zz=separate_roots(Pol(r.curve(y=r.points[i])), 10^4)
    if length(zz)>1
      m=nearest_pair(zz); push!(mins,(m[1],i))
      if VK.showZeros println(" d==",approx(m[1])," for ",m[2]) end
    end
    zz
  end
  if VK.showWorst && length(r.zeros[1])>1
    sort!(mins)
    println("worst points:")
    for m in mins[1:min(5,length(mins))]
      println("mindist(zeros[",m[2],"])==",approx(m[1]))
    end
  end
  if isempty(r.zeros[1]) return TrivialCase(r) end
#---------------------- zeros computed -------------------------------
  r.B=BraidMonoid(coxsym(length(r.zeros[1])))
  r.monodromy=fill(r.B(),length(r.segments))
  for segno in eachindex(r.segments)
    tm=time()
    r.monodromy[segno]=(VK.approx_monodromy ?
                    approx_follow_monodromy : follow_monodromy)(r, segno)
    tm=time()-tm
    if VK.showSegments
      println("monodromy[$segno]=",r.monodromy[segno])
      println("#segment $segno/",length(r.segments)," Time==",approx(tm),"sec")
    end
  end
#---------------------- braids --------------------------
  if VK.showgetbraid println("# Computing monodromy braids") end
  r.braids=empty([r.B()])
  if VK.showgetbraid println("loops=[") end
  for i in eachindex(r.loops)
    l=filter(s->!(r.monodromy[s]!==nothing),abs.(r.loops[i]))
    if length(l)>0 
      if VK.showgetbraid println("# loop[$i] missing segments ",l) end
    else
      bb=prod(s->s<0 ? r.monodromy[-s]^-1 : r.monodromy[s],r.loops[i])
      if VK.showgetbraid println("r.",bb,",") end
      push!(r.braids, bb)
    end
  end
  if VK.showgetbraid println("]") end
#---------------------- end braids --------------------------
  if r.ismonic 
    if VK.shrinkBraid
      r.rawBraids=r.braids
      r.braids=shrink(r.braids)
    end
    F=VKquotient(r.braids)
  else F=DBVKquotient(r)
  end
  r.presentation=Presentation(F)
  r.rawPresentation=Presentation(F)
  simplify(r.presentation)
  r
end

#---------------------- simp ----------------------------------
# simplest fraction approximating t closer than prec
function simp(t0::Real,prec=10^-15)
  t=t0
  a=BigInt[]
  k=BigInt[1,0]
  h=BigInt[0,1]
  while abs(h[end]//k[end]-t0)>prec
   n=floor(BigInt,t)
   push!(a,n)
   push!(h,n*h[end]+h[end-1])
   push!(k,n*k[end]+k[end-1])
   t=1/(t-n)
  end
  h[end]//k[end]
end

simp(t::Complex,prec=10^-15)=simp(real(t),prec)+im*simp(imag(t),prec)

#---------------------- root-finding ----------------------------------
"""
`NewtonRoot(p::Pol,initial_guess,precision;showall=false,show=false,lim=800)`

Here   `p`  is   a  polynomial   with  `Rational`   or  `Complex{Rational}`
coefficients.  The function  computes an  approximation to  a root  of `p`,
guaranteed of distance closer than `precision` to an actual root. The first
approximation  used is `initial`.  A possibility is  that the Newton method
starting  from `initial` does not converge  (the number of iterations after
which  this is decided  is controlled by  `lim`); then the function returns
`nothing`.  Otherwise the function returns a pair: the approximation found,
and an upper bound on the distance between that approximation and an actual
root.  The point of returning  an upper bound is  that it is usually better
than the asked-for `precision`. For the precision estimate a good reference
is

J.  Hubbard, D. Schleicher,  and S. Sutherland.  "How to find  all roots of
complex polynomials by Newton's method.", Invent. Math. 146:1--33, 2001.
```julia-repl
julia> p=Pol([1,0,1])
Pol{Int64}: x²+1

julia> NewtonRoot(p,1+im,10^-7)
2-element Vector{ComplexF64}:
  8.463737877036721e-23 + 1.0im
 1.8468154345305913e-11 + 0.0im

julia> NewtonRoot(p,1,10^-7;show=true)
****** Non-Convergent Newton after 800 iterations ******
p=x²+1 initial=-1.0 prec=1.0000000000000004e-7
```
"""
function NewtonRoot(p::Pol,z,precision;showall=VK.showallnewton,
                          show=VK.showNewton,lim=VK.NewtonLim)
  deriv=derivative(p)
  for cnt in 1:lim
    a=p(z)
    b=deriv(z)
    c=iszero(b) ? a : a/b
    err=abs(c)
    if iszero(err) err=(precision/100)/(degree(p)+1) end
    if err>precision err=precision end
    z=simp(z-c,precision/100/(degree(p)+1)/2)
    if showall println(cnt,": ",z) end
    if err<=(precision/100)/(degree(p)+1)
      if show print(cnt,":") end
      return (z,err)
    end
  end
  if show
    println("\n****** Non-Convergent Newton after ", lim," iterations ******")
    @show p,z,precision
    return nothing
  end
end

"""
'separate_roots_initial_guess(p::Pol, v, safety)'

Here  `p` is a complex  polynomial, and `v` is  a list of approximations to
roots  of `p` which should lie in different attraction basins for Newton' s
method.  The  result  is  a  list  `l`  of  complex  rationals representing
approximations  to the  roots of  `p` (each  element of  `l` is the root in
whose attraction basin the corresponding element of `v` lies), such that if
`d`  is the minimum distance  between two elements of  `l`, then there is a
root  of `p` within radius  `d/(2*safety)` of any element  of `l`. When the
elements  of  `v`  do  not  lie  in  different  attraction basins (which is
necessarily the case if `p` has multiple roots), `nothing` is returned.

```julia-repl
julia> p=Pol([1,0,1])
Pol{Int64}: x²+1

julia> separate_roots_initial_guess(p,[1+im,1-im],10^5)
2-element Vector{ComplexF64}:
 8.463737877036721e-23 + 1.0im
 8.463737877036721e-23 - 1.0im

julia> separate_roots_initial_guess(p,[1+im,2+im],1000)
    # 1+im and 2+im in same attraction basins
```
"""
function separate_roots_initial_guess(p, v, safety)
  if degree(p)==1 return [-p[0]/p[1]] end
  radv=nearest_pair(v)[1]/safety/2
  res=map(e->NewtonRoot(p,e,radv),v)
  if !any(isnothing,res) && nearest_pair(first.(res))[1]/safety/2>maximum(last.(res))
    return first.(res)
  end
  print("dispersal required=",nearest_pair(first.(res))[1]/safety/2)
  println(" obtained=",approx(maximum(last.(res))))
  println(join(v[nearest_pair(first.(res))[2]]," and ")," lie in the same attraction basin")
end

"""
'separate_roots(<p>, <safety>)'

Here  `p` is  a complex  polynomial. The  result is  a list  `l` of complex
numbers  representing approximations to the roots  of `p`, such that if `d`
is  the minimum distance between two elements  of `l`, then there is a root
of  `p` within  radius `d/(2*safety)`  of any  element of  `l`. This is not
possible when `p` has multiple roots, in which case `nothing` is returned.

```julia-repl
julia> @Pol q
Pol{Int64}: q

julia> separate_roots(q^2+1,100)
2-element Vector{ComplexF64}:
  2.3541814200656927e-43 + 1.0im
 -2.3541814200656927e-43 - 1.0im

julia> separate_roots((q-1)^2,100)

julia> separate_roots(q^3-1,100)
3-element Vector{ComplexF64}:
 -0.5 - 0.8660254037844386im
  1.0 - 1.232595164407831e-32im
 -0.5 + 0.8660254037844387im
```
"""
function separate_roots(p,safety)
  subtractroot(p,r)=divrem(p,Pol([-r,1]))[1]
  if p isa Mvp p=Pol(p) end
  if degree(p)<1 return empty(p.c)
  elseif degree(p)==1 return [-p[0]/p[1]]
  end
  p//=p[end]
# e=complex(E(7))
  e=big(simp(complex(E(7))))
  v=nothing
  cnt = 0
  while isnothing(v) && cnt<2*(degree(p)+1)
    if VK.showNewton && cnt>0
      println("****** ", cnt, " guess failed for p degree ", degree(p))
    end
    v=NewtonRoot(p,e,(1/safety)*10.0^(-degree(p)-4))
    e*=simp(complex(5//4*E(2*(degree(p)+1))))
#   e*=complex(5//4*E(2*(degree(p)+1)))
    cnt+=1
  end
  if cnt>=2*(degree(p)+1) error("no good initial guess") end
  v=[v[1]]
  append!(v,separate_roots(subtractroot(p,v[1]), safety))
  safety==0 ? v : separate_roots_initial_guess(p, v, safety)
end

"""
`FindRoots(p::Pol, <approx>)`

`p`  should have rational or  `Complex{Rational} coefficients. The function
returns  'Complex' rational  approximations to  the roots  of `p` which are
better  than  `approx`  (a  positive  rational).  Contrary to the functions
`separate_roots`,  etc... described in the  previous chapter, this function
handles  quite  well  polynomials  with  multiple  roots.  We  rely  on the
algorithms explained in detail in cite{HSS01}.

```julia-repl
julia> FindRoots((Pol()-1)^5,1/5000)
5-element Vector{Complex{Rational{BigInt}}}:
 1//1 + 0//1*im
 1//1 + 0//1*im
 1//1 + 0//1*im
 1//1 + 0//1*im
 1//1 - 1//4837347*im

julia> FindRoots(Pol()^3-1,10^-5)
3-element Vector{Complex{Rational{BigInt}}}:
 -1//2 - 16296//18817*im
  1//1 + 0//1*im
 -1//2 + 16296//18817*im

julia> round.(Complex{Float64}.(ans.^3);sigdigits=3)
3-element Vector{ComplexF64}:
 1.0 - 1.83e-9im
 1.0 + 0.0im
 1.0 + 1.83e-9im
```
"""
function FindRoots(p,prec)
  subtractroot(p,r)=divrem(p,Pol([-r,1]))[1]
  if degree(p)<1 return empty(p.c)
  elseif degree(p)==1 return [-p[0]//p[1]]
  end
  e=big(simp(complex(E(7))))
  v=nothing
  while isnothing(v)
    v=NewtonRoot(p,e,10.0^(-degree(p)-1))
    e*=simp(complex(E(degree(p)+1)))
  end
   v=vcat([v[1]],FindRoots(subtractroot(p,v[1]),prec))
   map(e->NewtonRoot(p,e,prec)[1],v)
end

#------------------ Loops --------------------------------------------
# sorts a list of points trigonometrically around a center
# starting from -im+ε and going anticlockwise
function cycorder(list, center)
  right=empty(list)
  left=empty(list)
  top=empty(list)
  bottom=empty(list)
  for y in list
    if real(y)>real(center) push!(right, y)
    elseif real(y)<real(center) push!(left, y)
    elseif imag(y)>imag(center) push!(top, y)
    else push!(bottom, y)
    end
  end
  sort!(right,by=x->imag(x-center)/real(x-center))
  sort!(left,by=x->imag(x-center)/real(x-center))
  vcat(right, top, left, bottom)
end
function cycorder2(list,center) # slightly slower
  angles=map(x->iszero(x-center) ? pi/2 : angle(-im*(x-center)),list)
  list[sortperm(angles)]
end

# Input: (l::Vector{Complex},center::Complex)
# Output: sublist of l in cycorder of "neighbours" of center,
# y is neighbour of center iff y≠center and no z∈l, z∉(y,center) is in the
# disk of diameter [y,center]
function neighbours(l, center)
  cycorder(filter(l)do y
    if y==center return false end
    for z in l
      if z==y || z==center continue end
      if abs2(y-z)+abs2(z-center)<=abs2(y-center) return false end
    end
    return true
  end,center)
end

# value at z of an equation of the line (x,y)
function lineq(x, y, z)
  if real(x)==real(y)
    if imag(x)==imag(y) error("Undefined line\n")
    else return real(z)-real(x)
    end
  else
    return (imag(y)-imag(x))*(real(z)-real(x))/(real(y)-real(x))+imag(x)-imag(z)
  end
end

# mediatrix of segment (x,y) of length abs2(x-y) on each side of segment
function mediatrix(x, y)
  if x==y error("Undefined mediatrix") end
  (x+y)/2 .+[im,-im].*(x-y)
end

crossing(v1,v2)=crossing(v1...,v2...)

# Computes the intersection of lines (x1,x2) and (y1,y2)
# returns nothing if the lines are parallel or a pair is a single
function crossing(x1,x2,y1,y2)
  if x1==x2 || y1==y2 return nothing end
  if !(real(x1)==real(x2))
    λx=(imag(x1)-imag(x2))/(real(x1)-real(x2))
    μx=-λx*real(x1)+imag(x1)
    if !(real(y1)==real(y2))
      λy=(imag(y1)-imag(y2))/(real(y1)-real(y2))
      μy=-λy*real(y1)+imag(y1)
      if λx==λy return nothing end
      resr=(μy-μx)/(λx-λy)
      res=resr+(λx*resr+μx)*im
      return res
    else
      E3=simp(complex(E(3)))
      res=crossing(E3*x1, E3*x2, E3*y1, E3*y2)
      if isnothing(res) return nothing end
      return res/E3
    end
  else
    res=crossing(im*x1, im*x2, im*y1, im*y2)
    if isnothing(res) return nothing end
    return res/im
  end
end

function detectsleftcrossing(c, w, y, z)
  res=fill(false,length(c)-1)
  a,b=mediatrix(y, z)
  for k in 1:length(c)-1
    if lineq(a, b, c[k])*lineq(a, b, c[k+1])<=0
      x=crossing(a, b, c[k], c[k+1])
      if !isnothing(x) res[k]=imag((z-y)/(w[k]-y))>=0 end
    end
  end
  res
end

# eliminates trivial segments and contracts pairs [a,b],[b,a]
function Garside.shrink(l)local k
  k=findfirst(i->l[i]==l[i+1],1:length(l)-1)
  if !isnothing(k) return shrink(vcat(l[1:k],l[k+2:end])) end
  k=findfirst(i->l[i]==l[i+2],1:length(l)-2)
  if !isnothing(k) return shrink(vcat(l[1:k],l[k+3:end])) end
  l
end

"""
`convert_loops(ll)`

The  input is a list  of loops, each a  list of complex numbers representing
the vertices of the loop.

The output is a named tuple with fields
  - `points`: a list of complex  numbers.
  - `segments`:  a list of oriented segments, each of them  encoded by the
    list of the positions in 'points' of  its two endpoints.
  - `loops`: a list of loops. Each loops is a list  of integers representing
    a  piecewise  linear  loop,  obtained  by  concatenating the `segments`
    indexed  by the  integers, where  a negative  integer is  used when the
    opposed orientation of the segment is taken.
"""
function convert_loops(ll)
  points=unique(vcat(ll...))
  points=points[filter(i->!any(==(points[i]),points[1:i-1]),eachindex(points))]
  points=sort(points,by=x->(imag(x),real(x)))
  np(p)=findfirst(==(p),points)
  loops=map(l->np.(l),ll)
  loops=shrink.(loops)
  loops=map(l->map(i->l[i-1:i],2:length(l)),loops)
  segments=sort(unique(sort.(vcat(loops...))))
  loops=map(loops)do l
    map(l)do seg
     seg[1]<seg[2] ? findfirst(==(seg),segments) :
                    -findfirst(==(reverse(seg)),segments)
    end
  end
  (;points, segments, loops)
end

function Box(l)
  minr,maxr=extrema(real.(l))
  mini,maxi=extrema(imag.(l))
  [Complex(minr-2, mini-2), Complex(minr-2, maxi+2),
   Complex(maxr+2, maxi+2), Complex(maxr+2, mini-2), 
   Complex((maxr+minr)/2, mini-2-(maxr-minr)/2),
   Complex(minr-2-(maxi-mini)/2, (maxi+mini)/2),
   Complex((maxr+minr)/2, maxi+2+(maxr-minr)/2),
   Complex(maxr+2+(maxi-mini)/2, (maxi+mini)/2)]
end

"""
'loops_around_punctures(points)'

`points`  should be complex numbers. The function computes piecewise-linear
loops representing generators of the fundamental group of `ℂ -{points}`.

```julia-repl
julia> loops_around_punctures([0])
1-element Vector{Vector{Complex{Int64}}}:
 [1 + 0im, 0 + 1im, -1 + 0im, 0 - 1im, 1 + 0im]
```
# Guarantees on loops_around_punctures:
# For a set Z of zeroes and z in Z, let R(z):=1/2 dist(z,Z-z).
# The  input of  loops_around_punctures is  a set  Z of approximate zeroes of
# r.discy such that for any z one of the zeroes is closer than R(z)/S where
# S is a global constant of the program (in practice we may take S=100).
# Let  d=inf_{z in  Z}(R(z)); we  return points  with denominator  10^-k or
# 10^-k<d/S'  (in practive we take S'=100) and  such that the distance of a
# segment to a zero of r.discy is guaranteed >= d-d/S'-d/S
"""
function loops_around_punctures(originalroots)
# tol=first(nearest_pair(originalroots))
  roots=originalroots*(1+0im)
  n=length(roots)
  if n==1 return [roots[1].+[1,im,-1,-im,1]] end
  average=sum(roots)/n
  sort!(roots, by=x->abs2(x-average))
  ys=map(x->(neighbours=Int[],friends=Int[], lovers=Int[],
             cycorder=empty(roots),circle=empty(roots),
             witness=empty(roots),path=empty(roots),handle=empty(roots),
             loop=empty(roots)),roots)
  err=filter(x->==(roots[x]...),combinations(eachindex(roots),2))
  if length(err)>0 error("roots too close ",err) end
  iy(y)=findfirst(==(y),roots)
  sy(y)=ys[iy(y)]
  for (yi,y) in enumerate(ys)
    append!(y.neighbours,iy.(neighbours(roots, roots[yi])))
    push!(y.friends,yi)
  end
  if VK.showLoops println("neighbours computed") end
  for (yi,y) in enumerate(ys)
    for z in y.neighbours
      if !(z in y.friends)
        push!(y.lovers, z)
        push!(ys[z].lovers, yi)
        newfriends=vcat(y.friends, ys[z].friends)
        for t in y.friends 
          empty!(ys[t].friends);append!(ys[t].friends,newfriends)
        end
        for t in ys[z].friends 
          empty!(ys[t].friends);append!(ys[t].friends,newfriends)
        end
      end
    end
  end
  for (yi,y) in enumerate(ys) 
    sort!(y.neighbours,by=z->abs2(roots[yi]-roots[z]))
  end
# To avoid trouble with points on the border of the convex hull,
# we make a box around all the points;
  box=Box(roots)
  for (yi,y) in enumerate(ys) 
    append!(y.cycorder,cycorder(vcat(deleteat!(copy(roots),yi),box),roots[yi]))
    n1=roots[y.neighbours[1]]
    k=findfirst(==(n1),y.cycorder)
    y.cycorder.=circshift(y.cycorder,1-k)
    push!(y.cycorder, y.cycorder[1])
    push!(y.circle,(roots[yi]+n1)/2)
    push!(y.witness,n1)
    for z in y.cycorder[2:end]
      cut=detectsleftcrossing(y.circle, y.witness, roots[yi], z)
      if any(cut)
        k=findfirst(cut)
        resize!(y.circle,k)
        resize!(y.witness,k)
      end
      k=length(y.circle)
      newcirc=crossing(mediatrix(roots[yi],y.witness[k]),mediatrix(roots[yi], z))
      if !isnothing(newcirc)
        push!(y.circle, newcirc)
        push!(y.witness, z)
      end
      if iy(z) in y.lovers
        push!(y.circle, (roots[yi]+z)/2)
        push!(y.witness, z)
      end
    end
  end
  if VK.showLoops println("circles computed") end

  function boundpaths(path, i) # y must be an element of ys
    if !isempty(ys[i].path) return end
    append!(ys[i].path,path);push!(ys[i].path,roots[i])
    for z in ys[i].lovers boundpaths(ys[i].path, z) end
  end

  boundpaths(empty(roots), 1)
  for (yi,y) in enumerate(ys) 
    k=length(y.path)
    if k>1
      circleorigin=(roots[yi]+y.path[k-1])/2
      k=findfirst(==(circleorigin),y.circle)
      y.circle.=circshift(y.circle,1-k)
    end
  end
  for y in ys
    k=length(y.path)
    append!(y.handle,vcat(map(1:k-1)do i
      l=sy(y.path[i]).circle
      l[1:findfirst(==((y.path[i]+y.path[i+1])/2),l)]
     end...))
    append!(y.loop,vcat(y.handle, y.circle, reverse(y.handle)))
  end
  ys=ys[sort(eachindex(ys),by=i->findfirst(==(roots[i]),originalroots))]
  map(ys)do y
#   loop=map(x->round(x;sigdigits=8),y.loop)
    y.loop
  end
end

#-------------------- ApproxMonodromy ----------------------------
# for each point of a find closest point in b
# Complain if the result is not a bijection between a and b of if
# the distance between an a and the corresponding b is bigger than 1/10
# of minimum distance between two b's
function fit(a, b)
  dm=map(p->findmin(abs.(b.-p)),a)
  monodromyError=maximum(first.(dm))
# println("# Monodromy error==",monodromyError)
  if monodromyError>nearest_pair(b)[1]/10 error("monodromy error too big") end
  pos=last.(dm)
  if sort(pos)!=1:length(pos) error("monodromy cannot find perm") end
  b[pos]
end

# Decimal Log of Norm of polynomial d evaluated at point p
function normdisc(d, p)
  p=abs(prod(map(f->f(p), d)))
  if log10(p)==0 return round(Float64(-log10(1/p));digits=3)
  else return round(Float64(log10(p));digits=3)
  end
end

# keep only 3 significant digits of x
approx(x::Real)=round(Float64(x);sigdigits=3)
approx(x::Complex)=round(Complex{Float64}(x);sigdigits=3)

"""
`approx_follow_monodromy(<r>,<segno>,<pr>)`

This function  computes an approximation  of the monodromy braid  of the
solution in `x`  of an equation `P(x,y)=0` along  a segment `[y₀,y₁]`.
It is called  by 'fundamental_group', once for each of  the segments. The
first  argument is  a  global record,  similar to  the  one produced  by
'fundamental_group'  (see the  documentation of  this function)  but only
containing intermediate information. The second argument is the position
of the segment in 'r.segments'. The  third argument is a print function,
determined  by the  printlevel set  by the  user (typically,  by calling
'fundamental_group' with a second argument).

Contrary to 'follow_monodromy',  'approx_follow_monodromy' does not control
the approximations; it just uses a  heuristic for how much to move along
the segment  between linear braid  computations, and this  heuristic may
possibly fail. However,  we have not yet found an  example for which the
result is actually incorrect, and thus the existence is justified by the
fact that  for some difficult  computations, it is sometimes  many times
faster  than 'follow_monodromy'.  We illustrate  its typical  output when
<printlevel> is 2.

|   VK.approx_monodromy=true
julia-rep1```
julia> fundamental_group((x+3*y)*(x+y-1)*(x-y);printlevel=2)

  ....

546 ***rejected
447<15/16>mindist=2.55 step=0.5 total=0 logdisc=0.55 ***rejected
435<15/16>mindist=2.55 step=0.25 total=0 logdisc=0.455 ***rejected
334<15/16>mindist=2.55 step=0.125 total=0 logdisc=0.412 ***rejected
334<15/16>mindist=2.55 step=0.0625 total=0 logdisc=0.393
334<15/16>mindist=2.55 step=0.0625 total=0.0625 logdisc=0.412
334<15/16>mindist=2.56 step=0.0625 total=0.125 logdisc=0.433
334<15/16>mindist=2.57 step=0.0625 total=0.1875 logdisc=0.455
334<15/16>mindist=2.58 step=0.0625 total=0.25 logdisc=0.477
======================================
==    Nontrivial braiding B(2)      ==
======================================
334<15/16>mindist=2.6 step=0.0625 total=0.3125 logdisc=0.501
334<15/16>mindist=2.63 step=0.0625 total=0.375 logdisc=0.525
334<15/16>mindist=2.66 step=0.0625 total=0.4375 logdisc=0.55
334<15/16>mindist=2.69 step=0.0625 total=0.5 logdisc=0.576
334<15/16>mindist=2.72 step=0.0625 total=0.5625 logdisc=0.602
334<15/16>mindist=2.76 step=0.0625 total=0.625 logdisc=0.628
334<15/16>mindist=2.8 step=0.0625 total=0.6875 logdisc=0.655
334<15/16>mindist=2.85 step=0.0625 total=0.75 logdisc=0.682
334<15/16>mindist=2.9 step=0.0625 total=0.8125 logdisc=0.709
334<15/16>mindist=2.95 step=0.0625 total=0.875 logdisc=0.736
334<15/16>mindist=3.01 step=0.0625 total=0.9375 logdisc=0.764
# Minimal distance==2.55
# Minimal step==0.0625==-0.0521 + 0.0104im
# Adaptivity==10
monodromy[15]=[2]

# segment 15/16 Time==0.002741098403930664sec
```

Here at each  step the following information is  displayed: first, how
many iterations of  the Newton method were necessary to  compute each of
the 3  roots of the current  polynomial `f(x,y₀)` if we  are looking at
the point `y₀` of the segment.  Then, which segment we are dealing with
(here the  15th of  16 in  all). Then the  minimum distance  between two
roots of  `f(x,y₀)` (used in our  heuristic). Then the current  step in
fractions of the length of the segment  we are looking at, and the total
fraction of the segment we have  done. Finally, the decimal logarithm of
the absolute  value of the discriminant  at the current point  (used in
the heuristic). Finally, an indication if the heuristic predicts that we
should  halve the  step  ('***rejected')  or that  we  may double  it
('***up').

The function returns an element of the ambient braid group 'r.B'.
"""
function approx_follow_monodromy(r,segno)
  if VK.showInsideSegments ipr=print
  else ipr=function(x...)end
  end
  p,q=r.segments[segno]
  res=r.B()
  prevzeros=r.zeros[p]
  n=length(prevzeros)
  if n==1 return r.B() end
  mindm=nearest_pair(prevzeros)[1]
  p=r.points[p]
  v=r.points[q]-p
  prev=p
  step=1//1
  minstep=step
  total=0//1
  nextzeros=nothing
  while true
    next=prev+step*v
    P=Pol(r.curve(y=next))
    nextzeros=separate_roots_initial_guess(P, prevzeros, 100)
    if isnothing(nextzeros) ||
      (iszero(maximum(abs.(nextzeros-prevzeros))) && step>1//16)
      rejected=true
    else
      dm=map(i->minimum(abs.(prevzeros[i]-prevzeros[j] for j in 1:n if j!=i)),
                                                                          1:n)
      mdm=minimum(dm)
      if step<1 ipr("<$segno/",length(r.segments),">mindist=",approx(mdm),
         " step=$step total=$total logdisc=",normdisc(r.discyFactored,next))
      end
      dn=abs.(prevzeros-nextzeros)
      rejected=any(dm.<VK.AdaptivityFactor.*dn)
      if !rejected && mdm<mindm mindm=mdm end
    end
    if rejected
      step/=2
      ipr(" ***rejected\n")
      if step<minstep minstep=step end
    else
      total+=step
      if all(dm.>2 .*VK.AdaptivityFactor.*dn) && total+step!=1
        step*=2
        ipr(" ***up")
      end
      ipr("\n")
      if total != 1
        res*=LBraidToWord(prevzeros, nextzeros, r.B)
        prevzeros=nextzeros
      end
      prev=next
    end
    if total+step>1 step=1-total end
    if total==1 break end
  end
  if VK.showSegments
    println("# Minimal distance=", approx(mindm))
    println("# Minimal step=", minstep, "=", approx(v*minstep))
    println("# Adaptivity=", VK.AdaptivityFactor)
  end
  res*LBraidToWord(prevzeros,fit(nextzeros,r.zeros[q]),r.B)
end
#-------------------- Monodromy ----------------------------
# ceil(-log2(p)) ifor 0<p<1
function Intlog2(p)
  k=0
  q=p
  while q<1
    q=2q
    k+=1
  end
  k
end

# computes the lower approximation of the rational a by a
# rational with denominator 2^k
function binlowevalf(a, time)
  k=Intlog2(a-time)+3
  b=floor(Int,a*2^k)
  a>=0 ? b//2^k : (b-1)//2^k
end

# one step of the Newton method
function mynewton(p,z)
  a=p(z)
  b=derivative(p)(z)
  if iszero(b) c=a
    println("NewtonError")
  else c=a/b
  end
  err=degree(p)*abs(c)
  if err==0 prec=1
  else prec=max(0,ceil(Int,-log10(err)))+2
  end
  simp(z-c,(1//10)^(prec+1))
end

# for each point of a find closest point in b
function myfit(a, b)
  d=length(a)
  dist=fill(zero(real(eltype(a))),d,d)
  for k in 1:d, l in k+1:d
    dist[k,k]=dist[l,k]=dist[k,l]=abs2(a[k]-a[l])
  end
  dist[d,d]=dist[d,d-1]
  R=map(k->minimum(dist[k,:])*1//4,1:d)
  map(k->only(filter(i->abs2(i-a[k])<R[k], b)),1:d)
end

# sets coeff of degree i of p to x
Base.setindex!(p::Pol{T},x::T,i::Integer) where T=p.c[i+1-p.v]=x

# Sturm(pp,time)
# if polynomial pp is positive  at time<1
# returns some rational number t such that
#    time<t<=1  and  pp  is positive on [time,t]
# otherwise returns 0
# [third input and second output is an adaptive factor to
#  accelerate the computation]
function Sturm(pp::Pol, tm, adapt::Integer)
  q=Pol()
  pol=pp((1-q)*tm+q)
  if pol[0]<=0
    print("*****",Float32(pol[0]))
    return [0, 0]
  end
  k=1
  while k<degree(pol) && pol[k]>=0 k+=1 end
  while k<degree(pol) && pol[k]<=0 k+=1 end
  for i in k:degree(pol) if pol[i]>0 pol[i]=zero(pol[i]) end end
  t=big(1//2)^adapt
  m=adapt
  while pol(t)<=0
    t//=2
    m+=1
  end
  if VK.showInsideSegments print(m) end
  if m==adapt && adapt>0
    if pol(3t//2)>0
      if pol(2t)>0 res=[(1-2t)*tm+2t, adapt-1]
      else res=[(1-3t//2)*tm+3t//2, adapt-1]
      end
    else res=[(1-t)*tm+t, adapt]
    end
  else res=[(1-t)*tm+t, m]
  end
  res
end

Base.real(p::Pol)=Pol(real.(p.c),p.v)
Base.imag(p::Pol)=Pol(imag.(p.c),p.v)
Base.abs2(p::Pol)=real(p)^2+imag(p)^2

# fraction 0≤tm≤1 in hexa width l
function formattm(tm,l)
  if iszero(tm) return rpad("0",l) end
  d=denominator(tm)
  n=numerator(tm)
  m=2^Int(4-mod1(log2(d),4))
# @show n,d,m
  d*=m
  n*=m
  rpad("0."*"0"^max(0,floor(Int,(log2(d)-log2(n+1))/4))*string(n,base=16),l)
end 

"""
`follow_monodromy(r,segno)`
This  function computes the  monodromy braid of  the solution in  `x` of an
equation   `P(x,y)=0`  along   a  segment   `[y₀,y₁]`.  It   is  called  by
`fundamental_group`, once for each of the segments. The first argument is a
global  record, similar to the one produced by `fundamental_group` (see the
documentation   of   this   function)   but  only  containing  intermediate
information.  The  second  argument  is  the  position  of  the  segment in
`r.segments`.

The function returns an element of the ambient braid group `r.B`.

This function has no reason to be called directly by the user, so we do not
illustrate  its behavior. Instead,  we explain what  is displayed on screen
when the user sets the printlevel to `2`.

What is quoted below is an excerpt of what is printed during the execution of
```julia_repl
julia> fundamental_group((x+3*y)*(x+y-1)*(x-y),printlevel=2)
......
<1/16>    1 time=0           ?2?1?3
<1/16>    2 time=0.2         R2. ?3
<1/16>    3 time=0.48        R2. ?2
<1/16>    4 time=0.74        ?2R1?2
<1/16>    5 time=0.94        R1. ?2
======================================
==    Nontrivial braiding B(2)      ==
======================================
<1/16>    6 time=0.bc        R1. ?1
<1/16>    7 time=0.d8        . ?0. 
<1/16>    8 time=0.dc        ?1R0?1
# The following braid was computed by follow_monodromy in 8 steps.
monodromy[1]=B(2)
#segment 1/16 Time==0.0536sec
```

`follow_monodromy` computes  its results by subdividing  the segment into
smaller  subsegments  on which  the  approximations  are controlled.  It
starts at one  end and moves subsegment after subsegment.  A new line is
displayed at each step.

The  first column  indicates which  segment is  studied. In  the example
above, the function  is computing the monodromy along  the first segment
(out  of  `16`).  This  gives  a  rough  indication  of  the  time  left
before  completion of  the total  procedure.  The second  column is  the
number of  iterations so  far (number of  subsegments). In  our example,
`follow_monodromy`  had to  cut the  segment into  `8` subsegments.  Each
subsegment has its own length. The cumulative length at a given step, as
a  fraction of  the  total length  of the  segment,  is displayed  after
`time=`.  This  gives  a  rough  indication  of  the  time  left  before
completion  of the  computation of  the monodromy  of this  segment. The
segment is completed when this fraction reaches `1`.

The last column has to do with the piecewise-linear approximation of the
geometric monodromy  braid. It is  subdivided into sub-columns  for each
string. In  the example above,  there are  three strings. At  each step,
some strings are fixed (they are  indicated by `. ` in the corresponding
column). A symbol like `R5` or `?3` indicates that the string is moving.
The exact meaning of the symbol has to do with the complexity of certain
sub-computations.

As  some strings  are moving,  it  happens that  their real  projections
cross. When such a crossing occurs, it is detected and the corresponding
element of `Bₙ` is displayed on screen (`Nontrivial braiding =`...) The
monodromy braid is the product of these elements of `Bₙ`, multiplied in
the order in which they occur.
"""
function follow_monodromy(r,seg)
# Exact computation of the monodromy braid along a segment
# r: global VK record
# seg: segment number
# sprint: Print function (to screen, to file, or none)
  iPrint=VK.showInsideSegments ? print : function(arg...) end
  p=r.curve
  dpdx=derivative(r.curve,:x)
  a,b=r.segments[seg]
  v=r.zeros[a]
  B=r.B
  res=B()
  # If there is only one string, the braid is trivial
  if length(v)==1 return res end
  d=length(r.zeros[1])
  t=Mvp(:t)
  tm=big(0)
  pt=p(;y=r.points[b]*t+r.points[a]*(1-t))
  dpdxt=dpdx(;y=r.points[b]*t+r.points[a]*(1-t))
  RR=fill(big(0.0),d)
  adapt=fill(0,d)
  protected=fill(0//1,d)
  protp=map(i->zero(Pol(real(v[1]))),1:d)
  protdpdx=map(i->zero(Pol(real(v[1]))),1:d)
  steps=0
  dist=fill(big(0.0),d,d)
  while true
    steps+=1
#   if steps>540 error() end
    iPrint("<$seg/",length(r.segments),">",lpad(steps,5))
    iPrint(" time=",formattm(tm,9),"   ")
    for k in 1:d, l in k+1:d
      dist[k,k]=dist[l,k]=dist[k,l]=abs2((v[k]-v[l])*big(1))
    end
    dist[d,d]=dist[d,d-1]
    for k in 1:d
      Rk=minimum(dist[k,:])/4
      z=v[k]
      if protected[k]>tm && Rk>=RR[k]
        iPrint(". ")
      elseif protected[k]>tm
        if adapt[k]+2<maximum(adapt) Rk/=2 end
        iPrint("R")
        s,adapt[k]=Sturm(Rk*protdpdx[k]-protp[k], tm, adapt[k])
        if s>tm protected[k]=binlowevalf(s,tm)
        else iPrint("How bizarre...")
#         @show Rk,protdpdx[k],protp[k]
#         @show Rk*protdpdx[k]-protp[k], tm, adapt[k]
        end
        RR[k]=Rk
      else
        iPrint("?")
        cptz=Pol(pt(;x=z))
#       @show pt, z,cptz
        protp[k]=d^2*abs2(cptz)
        cdpdxtz=Pol(dpdxt(;x=z))
        protdpdx[k]=abs2(cdpdxtz)
        s,adapt[k]=Sturm(Rk*protdpdx[k]-protp[k], tm, adapt[k])
        if s>tm protected[k]=binlowevalf(s,tm)
        else error("Something's wrong...s=",s,"≤time=",tm)
#         @show R[k],protdpdx[k],protp[k]
#         @show R[k]*protdpdx[k]-protp[k], tm, adapt[k]
        end
        RR[k]=Rk
      end
    end
    allowed=minimum(protected)
    tm=allowed
    py=Pol(p(;y=r.points[a]*(1-tm)+r.points[b]*tm))
    iPrint("\n")
    newv=map(1:d)do k
      if protected[k]>allowed v[k]
      else mynewton(py,v[k])
      end
    end
    res*=LBraidToWord(v, newv, B)
    v=newv
    if tm==1 break end
  end
  if VK.showSegments
    println("# The following braid was computed by follow_monodromy in $steps steps.")
  end
  res*LBraidToWord(v, myfit(v, r.zeros[b]), B)
end

#------------------- Compute PLBraid ----------------------------------
# Deals with "star" linear braids, those with associated permutation w₀
function starbraid(y, offset, B)
  n=length(y)
  if n==1 return B() end
  k=argmin(y)
  B((k:n-1).+offset...)*starbraid(deleteat!(copy(y),k),offset,B)/
  B((n+1-k:n-1).+offset...)
end

# In case two points have the same real projection, we use
# a "lexicographical" desingularization by "infinitesimal rotation"
function desingularized(v1, v2)
  n=length(v1)
  tan=1
  for k in 1:n, l in k+1:n
    rv,iv=reim(v1[k]-v1[l])
    if !iszero(iv*rv) tan=min(tan,abs(rv/iv)) end
    rv,iv=reim(v2[k]-v2[l])
    if !iszero(iv*rv) tan=min(tan,abs(rv/iv)) end
  end
  [v1, v2].*(1-im*tan/2)
end

"""
`LBraidToWord(v1,v2,B)`

This function converts  the linear braid joining the points in `v1` to the
corresponding ones in `v2` into an element of the braid group.

```julia-repl
julia> B=BraidMonoid(coxsym(3))
BraidMonoid(𝔖 ₃)

julia> LBraidToWord([1+im,2+im,3+im],[2+im,1+2im,4-6im],B)
1
```

The lists `v1` and `v2` must have the same length, say `n`. Then `B` should
be  `BraidMonoid(coxsym(n))`, the braid group  on `n` strings. The elements
of  `v1` (resp. `v2`)  should be `n`  distinct complex rational numbers. We
use the Brieskorn basepoint, namely the contractible set `C+iV_ℝ` where `C`
is  a real chamber; therefore the endpoints  need not be equal. The strings
defined  by `v1` and `v2` should be  non-crossing. When the numbers in `v1`
(resp.  `v2`)  have  distinct  real  parts,  the  real picture of the braid
defines a unique element of `B`. When some real parts are equal, we apply a
lexicographical  desingularization, corresponding to a rotation of `v1` and
`v2` by an arbitrary small positive angle.
"""
function LBraidToWord(v1, v2, B)
# two printlevel control fields: VK.showSingularProj
#				 VK.showBraiding
# Deals with linear braids
# 1) singular real projections are identified
# 2) calls starbraid for each
  n=length(v1)
  x1,y1=reim(v1)
  x2,y2=reim(v2)
  if length(Set(x1))<length(x1) || length(Set(x2))<length(x2)
    if VK.showSingularProj
      println("WARNING: singular projection(resolved)")
    end
    return LBraidToWord(desingularized(v1, v2)..., B)
  end
  q=sortPerm(x1)
  crit=empty(x1)
  for i in 1:n-1, j in i+1:n
    iq=i^q;jq=j^q
    if x2[iq]>x2[jq]
      push!(crit,(x1[iq]-x1[jq])/((x2[jq]-x1[jq]+x1[iq])-x2[iq]))
    end
  end
  tcrit=unique(sort(crit))
  res=B()
  u=0
  for t in tcrit
    xt=map(k->x1[k]+t*(x2[k]-x1[k]),1:n)
    yt=map(k->y1[k]+t*(y2[k]-y1[k]),1:n)
    ut=(u+t)/2
    xut=map(k->x1[k]+ut*(x2[k]-x1[k]),1:n)
    put=inv(sortPerm(xut))
    xt=permute(xt,put)
    yt=permute(yt,put)
    xcrit=unique(sort(xt))
    for x in xcrit
      posx=findfirst(==(x),xt)
      nx=count(==(x),xt)
      res*=starbraid(yt[posx:posx+nx-1], posx-1, B)
    end
    u=t
  end
  if VK.showBraiding
   if !isempty(tcrit)
      if VK.showInsideSegments
        println("======================================")
        println("==    Nontrivial braiding ",rpad(res,10),"==")
        println("======================================")
#       print("v1:=",gap(v1),";\n")
#       print("v2:=",gap(v2),";\n")
      else println("==    Nontrivial braiding ",rpad(res,10),"==")
      end
    end
  end
  return res
end
#----------------------- Presentation -------------------------------
"""
`VKquotient(braids)`

The  input `braids` is a list `b₁,…,bᵣ`, living in the braid group
on `n` strings. Each `bᵢ` defines by Hurwitz action an automorphism `φᵢ` of
the free group `Fₙ`. The function returns the group defined by the abstract
presentation: ``< f₁,…,fₙ ∣ ∀ i,j φᵢ(fⱼ)=fⱼ >``

```julia-repl
julia> B=BraidMonoid(coxsym(3))
BraidMonoid(𝔖 ₃)

julia> g=VKquotient([B(1,1,1),B(2)])
FreeGroup(a,b,c)/[b⁻¹a⁻¹baba⁻¹,b⁻¹a⁻¹b⁻¹aba,.,.,cb⁻¹,c⁻¹b]

julia> p=Presentation(g)
Presentation: 3 generators, 4 relators, total length 16

julia> display_balanced(p)
1: c=b
2: b=c
3: bab=aba
4: aba=bab

julia> simplify(p)
Presentation: 2 generators, 1 relator, total length 6
Presentation: 2 generators, 1 relator, total length 6

julia> display_balanced(p)
1: bab=aba
```
"""
function VKquotient(braids)
  F=FpGroup(Symbol.('a'.+(0:ngens(braids[1].M.W)))...)
  F/reduce(vcat,map(b->hurwitz(gens(F),b).*inv.(gens(F)),braids))
end

# A variant of VKquotient
# See arXiv:math.GR/0301327 for more mathematical details.
# Input: global VK record
# Output: the quotient, encoded as an FpGroup
function DBVKquotient(r)
  # get the true monodromy braids and the Hurwitz action basic data
  n=ngens(r.braids[1].M.W)+1
  F=FpGroup(Symbol.('a'.+(0:n+length(r.verticallines)-1))...)
# above the basepoint for the loops, locate the position of the string
# corresponding to the trivializing horizontal line
  bzero=r.zeros[r.basepoint]
  height=bzero[argmin(abs2.(bzero.-r.height))]
  fbase=F(count(z->reim(z)<=reim(height),bzero))
  rels=AbsWord[]
  auts=map(b->Hom(F,F,hurwitz(gens(F),b)),r.braids)
  for (i,aut) in enumerate(auts)
# Find an element conjugator such that aut(fbase)^inv(conjugator)=fbase
    ifbase=aut(fbase)
    conjugator=one(F)
    while length(ifbase)>1
      x=ifbase[1:1]
      ifbase=ifbase^x
      conjugator*=x
    end
# Replacing aut by  correctaut:= Conj(conjugator)*aut
    conj=Hom(F, F, gens(F).^Ref(conjugator))
    correctaut=x->conj(aut(x))
    g=i>length(r.verticallines) ? one(F) : F(i+n)
    append!(rels, map(f->correctaut(f)*g*inv(g*f), gens(F)[1:n]))
  end
  push!(rels, fbase)
  F/rels
end

gg(x)="Complex(evalf(\"$(real(x))\"),evalf(\"$(imag(x))\"))"
data=Dict()

x=Mvp(:x);y=Mvp(:y);z=Mvp(:z);t=Mvp(:t)
data[23]=discriminant(crg(23))(x,y,z)(;x=1,z=x)
data[24]=discriminant(crg(24))(x,y,z)(;x=1,z=x)
data[27]=discriminant(crg(27))(x,y,z)(;x=1,z=x)
data[29]=discriminant(crg(29))(x,y,z,t)(;t=y+1,z=x)
data[31]=discriminant(crg(31))(x,y,z,t)(;t=x+1,z=y)
data[34]=
95864732434895657396628326400//164799823*x*y^3-598949723065092000//
1478996726054382501274179923886253687929138281*x*y^7-
67840632073999787861633181671139840000*x^2-
7622790471072621273612030528032173587500421120000*y^2-273861000//
27158981660831329*x^2*y^4+37130333513291749382400//7130353846013*x^3*y-
13608525//50841945969352380915996169*x^4*y^2-2606867429323404078970327//
1675017448527954334139901265590107596081497211494411528*x^6+
3269025273548225517660538475128200000//390195840687434028022928452202401489*y^6
#---------plotting
function _segs(r,v)
  rp=Float64[]
  ip=Float64[]
  for seg in v
    s=r[:points][seg<0 ? reverse(r[:segments][-seg]) : r[:segments][seg]]
    append!(rp,real.(s))
    append!(ip,imag.(s))
  end
  rp,ip
end

function _loops(r,v)
  rp=Float64[]
  ip=Float64[]
  for l in r[:loops][v]
    nrp,nip=_segs(r,l)
    append!(rp,nrp)
    append!(ip,nip)
  end
  rp,ip
end

function plotloops(r,v)
 colors=[:black,:green,:blue,:red,:yellow,:cyan,:magenta]
 plt=lineplot(_loops(r,v[1])...,color=:black)
 for i in 2:length(v)
   lineplot!(plt,_loops(r,v[i])...,color=colors[i])
 end
 plt
end

#end
