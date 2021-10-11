#------------------- utilities -------------------------------
"""
`Dispersal(v)`
`v`  is a list of complex numbers. The result is a pair whose first element
is the minimum distance (in the complex plane) between two elements of `v`,
and  the  second  is  a  pair  of  indices `[i,j]` such that `v[i]`, `v[j]`
achieves this minimum distance.

julia> Dispersal([1+im,0,1])
(1, [1,3])
"""
function Dispersal(v)
  l=combinations(eachindex(v),2)
  m,c=findmin(map(x->abs(v[x[1]]-v[x[2]]),l))
  (m,l[c])
end

# distance of z to segment [a,b]
function DistSeg(z,a,b)
  b-=a
  z-=a
  r=abs(b)
  z*=r/b
  real(z)<0 ? abs(z) : real(z)>r ? abs(z-b) : imag(z)>0 ? imag(z) : -imag(z)
end

"""
`discy2(p::Mvp)`

The input should be an 'Mvp' in 'x' and 'y', with rational coefficients.
The  function  returns the  discriminant  of  `p`  with respect  to  'x'
(an  `Mvp` in  `y`);  it uses  interpolation to  reduce  the problem  to
discriminants of univariate polynomials,  and works reasonably fast (not
hundreds of times slower than MAPLE...).

|    gap> discy2(x+y^2+x^3+y^3);
    4+27y^4+54y^5+27y^6|
"""
function discy2(p)
  n=2*(1+degree(p,:x))*(1+degree(p,:y))
  v=map(1:n)do i
    q=scal.(Pol(p(y=i),:x).c)
    if isempty(q) return 0
    elseif length(q)==1 return q[1]
    else return GLinearAlgebra.det(resultant(q,map(*,q[2:end],1:length(q)-1)))
    end
  end
  Pol(1:n,v)
end

Chars.discriminant(p::Pol)=GLinearAlgebra.det(resultant(p[0:end],
                                                 derivative(p)[0:end]))

discy(p)=Pol(discriminant(Pol(p,:x)))

"""
`resultant(v,w)`

`v` and  `w` are vectors  representing coefficients of  two polynomials.
The function returns  Sylvester matrix for these  two polynomials (whose
determinant  is  the resultant  of  the  two  polynomials). It  is  used
internally by Discy.

```julia-repl
julia> @Mvp x,y

julia> p=x+y^2+x^3+y^3
Mvp{Int64}: x³+x+y³+y²

julia> c=Pol(p,:x).c
4-element Vector{Mvp{Int64, Int64}}:
 y³+y²
 1
 0
 1

julia> d=Pol(derivative(p,:x),:x).c
3-element Vector{Mvp{Int64, Int64}}:
 1
 0
 3

julia> resultant(c,d)
5×5 Matrix{Mvp{Int64, Int64}}:
 1  0  1  y³+y²  0
 0  1  0  1      y³+y²
 3  0  1  0      0
 0  3  0  1      0
 0  0  3  0      1

julia> GLinearAlgebra.det(m)
Mvp{Int64}: 27y⁶+54y⁵+27y⁴+4
```
"""
function resultant(v,w)
  v=reverse(v)
  w=reverse(w)
  l=max(0,length(v)+length(w)-2)
  m=fill(zero(v[1]),l,l)
  for i in 1:length(w)-1 m[i,i:i+length(v)-1]=v end
  for i in 1:length(v)-1 m[i+length(w)-1,i:i+length(w)-1]=w end
  m
end
#---------------------- global functions --------------------------
VKCURVE=Dict(
:name=>"vkcurve",
:version=>"2.0",
:date=>[2009,3],
:homepage=>"http://webusers.imj-prg.fr/~jean.michel/vkcurve.html",
:copyright=>
"(C) David Bessis, Jean Michel -- compute Pi_1 of hypersurface complements",
:monodromyApprox=>false, ##########################
:showSingularProj=>false,
:showBraiding=>false,
:showLoops=>false,
:showAction=>false,
:showSegments=>false,
:showInsideSegments=>false,
:showWorst=>false,
:showZeros=>false,
:showNewton=>false,
:showgetbraid=>false,
:showRoots=>false,
:showallnewton=>false,
:NewtonLim=>800,
:AdaptivityFactor=>10,
:shrinkBraid=>false,
:mvp=>1)

function SetPrintLevel(printlevel)
  VKCURVE[:showSingularProj]=VKCURVE[:showBraiding]=VKCURVE[:showLoops]=
  VKCURVE[:showAction]=VKCURVE[:showInsideSegments]=VKCURVE[:showWorst]=
  VKCURVE[:showZeros]=VKCURVE[:showNewton]=VKCURVE[:showRoots]=printlevel>=2
  VKCURVE[:showSegments]=VKCURVE[:showgetbraid]=printlevel>=1
end

SetPrintLevel(2) # while developping

@GapObj struct VK{T}
  curve::Mvp{T,Int}
  ismonic::Bool
end

function Base.show(io::IO,r::VK)
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
  r.basepoint=r.loops[1][1]<0 ? r.segments[-r.loops[1][1]][2] :
                                r.segments[r.loops[1][1]][1]
  r.loops=r.loops[indexin(r.ismonic ? r.roots : r.trueroots,r.roots)]
  segmentNumbers=sort(union(map(x->abs.(x),r.loops)))
  r.segments=r.segments[segmentNumbers]
  r.loops=map(x->map(y->y<0 ? -findfirst(==(-y),segmentNumbers) :
                     findfirst(==(y),segmentNumbers), x),r.loops)
  uniquePoints=sort(union(r.segments))
  r.points=r.points[uniquePoints]
  r.segments=map(x->indexin(x,uniquePoints), r.segments)
  r.basepoint=findfirst(==(r.basepoint),uniquePoints)
  if VKCURVE[:showSegments]
    println("# There are ",length(r.segments)," segments in ",
            length(r.loops)," loops")
  end
  if VKCURVE[:showWorst]
    l=map(enumerate(r.segments))do (i,s)
      m,ixm=findmin(DistSeg.(r.roots, r.points[s[1]], r.points[s[2]]))
      (m,i,ixm)
    end
    sort!(l)
    print("worst segments:\n")
    for i in 1:min(5,length(l))
      d,s,s1=l[i]
      println("segment ",s,"==",r.segments[s]," dist to ",s1,"-th root is ",d)
    end
  end
# find the minimum distance m between two roots
  if length(r.roots)>1
    m=Dispersal(r.roots)
    if VKCURVE[:showRoots]
      print("\nMinimum distance==",m[1]," between roots ",m[2][1]," and ",m[2][2]," of discriminant\n")
    end
    r.dispersal=m[1]
  else
    r.dispersal=1/1000
  end
# and round points to m/100
end

# Zeros(r)  r should have fields
#   r.points
#   r.curve
# It computes r.zeros: r.zeros[i]=zeros of r.curve(y=r.points[i])
function Zeros(r)
  if VKCURVE[:showRoots]
    println("Computing zeros of curve at the ", length(r.points), " segment extremities...")
  end
  mins=Tuple{Float64,Int}[]
  r.zeros=map(1:length(r.points))do i
    if VKCURVE[:showZeros] print("<",i,"/",length(r.points),">") end
    zz=SeparateRoots(Pol(r.curve(y=r.points[i])), 10^4)
    if length(zz)>1
      m=Dispersal(zz)
      if VKCURVE[:showZeros] println(" d==", m) end
      push!(mins,(m[1],i))
    end
    zz
  end
  if VKCURVE[:showWorst] && length(r.zeros[1])>1
    sort!(mins)
    println("worst points:")
    for i in 1:min(5,length(mins))
      println(mins[i][2],": mindist(zeros)==",mins[i][1])
    end
  end
end

# PrepareCurve(curve)
#   curve should be an Mvp in x and y.
#   This  function  makes  sure  the  curve  is  quadratfrei  and makes its
#   coefficients  in  CF(4)=GaussianRationals (the coeffs could
#   be  decimals or complex  decimals) returns a  record with fields .curve
#   and .ismonic if the curve is monic in x
function PrepareCurve(curve::Mvp)
  if !issubset(variables(curve),[:x,:y])
    error(curve," should be an Mvp in x,y")
  end
  d=gcd(curve, derivative(curve,:x))
  if degree(d,:x)>0
    xprintln("**** Warning: curve is not quadratfrei: dividing by ", d)
    curve=exactdiv(curve,d)
  end
  VK(curve,degree(Pol(curve,:x)[end])==0,Dict{Symbol,Any}())
end

complexpol(p::Pol)=Pol(Complex{Float64}.(p.c),p.v)

# VKCURVE.Discy(r)
#  r should be a record with field r.curve, a quadratfrei Mvp in x,y.
#  The discriminant of this curve with respect to x (a polynomial in y)
#  is computed.
#  First, the curve is split in
#    r.curveVerticalPart  -- the Gcd of the coeffs in x (an Mvp in y).
#    r.nonVerticalPart    -- curve/curveVerticalPart
#  Then disc=discriminant of r.nonVerticalPart is computed. Its quadratfrei
#  part  is computed, stripped  of factors common  with d and then factored
#  (if possible which in GAP3 means it is a polynomial over the rationals),
#  disc  is  stored  in  r.discy  as  a  GAP polynomial, and its factors in
#  r.discyFactored as a list of Mvp in x.
#  Some of  these computations may be  too costly for GAP,  in which case
#  one has a better hope to complete them in MAPLE.
function Discy(r)
  r.curveVerticalPart=gcd(Pol.(values(coefficients(r.curve,:x))))
  if VKCURVE[:showRoots] && degree(r.curveVerticalPart)>0
    println("Curve has ",degree(r.curveVerticalPart)," linear factors in y")
  end
  r.nonVerticalPart=exactdiv(r.curve,r.curveVerticalPart(Mvp(:y)))
  d=discy(r.nonVerticalPart)
  if iszero(d)
    error("Discriminant is 0 but ", r.curve," should be quadratfrei")
  end
  if VKCURVE[:showRoots] print("Discriminant has ",degree(d)," roots, ") end
  d=exactdiv(d,gcd(d,derivative(d)))
  if VKCURVE[:showRoots] println(" of which ", degree(d), " are distinct") end
  common=gcd(d,r.curveVerticalPart)
  if VKCURVE[:showRoots] && degree(common)>0
    println(" and of which ",degree(common)," are roots of linear factors")
  end
  d=exactdiv(d,common)
  d//=d[end]
  r.discy=d
# r.discyFactored=[d]
  r.discyFactored=factor(r.discy)
end

complexmvp(p::Mvp)=Mvp(ModuleElt(m=>Complex{Float64}(c) for (m,c) in p.d))

function GetDiscyRoots(r)
  if VKCURVE[:showRoots] println("Computing roots of discriminant...") end
  if degree(r.curveVerticalPart)==0 r.roots=Complex{Float64}[]
  else r.roots=SeparateRoots(r.curveVerticalPart, 1000)
  end
  append!(r.roots,vcat(map(p->SeparateRoots(p,1000),r.discyFactored)...))
end

function Braids(r)
 pr=VKCURVE[:showgetbraid] ? println : function(x...)end
  pr("# Computing monodromy braids")
  r.braids=[]
  for i in eachindex(r.loops)
    l=filter(s->!(r.monodromy[s]!==nothing),abs.(r.loops[i]))
    if length(l)>0 pr("# loop[$i] missing segments ",l)
    else
      bb=prod(s->s<0 ? r.monodromy[-s]^-1 : r.monodromy[s],r.loops[i])
      pr("# loop[",i,"]==",bb)
      push!(r.braids, bb)
    end
  end
  if VKCURVE[:shrinkBraid]
    r[:rawBraids]=r.braids
    r.braids=shrink(r.braids)
  end
end

function Segment(r,segno)
  if haskey(r,:name)
    PrintTo(string(r.name,".",segno), "")
    pr=function(arg...,)
      AppendTo(string(r.name,".",segno),arg...)
      println(arg...)
    end
  elseif VKCURVE[:showSegments] pr=print
  else pr=function(x...)end
  end
  tm=time()
  r.monodromy[segno]=(VKCURVE[:monodromyApprox] ?
                  ApproxFollowMonodromy : FollowMonodromy)(r, segno, pr)
  tm=time()-tm
  if haskey(r, :name) pr(r.name, ".") end
  pr("monodromy[", segno,"]=",word(r.monodromy[segno]), "\n\n")
  pr("# segment ",segno,"/",length(r.segments)," Time==",tm,"sec\n")
end

function TrivialCase(r)
  r.presentation=Presentation(FpGroup(degree(r.curve,:x)))
  r
end

# VKCURVE.Segments(name[,range])
# builds  files  describing  monodromy  along  segments  of r.segments (all
# segments by default, those in range if given).
# If  r  has  a  .name  field  then  computed segments are written on files
# beginning with r.name, otherwise are added to r.monodromy
# if r does not have the right fields then one tries to read (r.name).tmp
function Segments(r,range=1:length(r.segments))
  read(r.name*".tmp")
  r.B=BraidMonoid(CoxSym(length(r.zeros[1])))
  if !haskey(r,:monodromy) r.monodromy=fill(r.B(),length(r.segments)) end
  for segno in range Segment(r, segno) end
end

function Finish(r)
  Braids(r)
  if r.ismonic F=VKQuotient(r.braids)
  else F=DBVKQuotient(r)
  end
  r.presentation=Presentation(F)
  r.rawPresentation=Presentation(F)
  return r
  simplify(r.presentation)
end

function SearchHorizontal(r)
  # Searching for a good horizontal
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
  r1=VK(r.curve*(Mvp(:x)-height),false,Dict{Symbol,Any}()) # is it monic?
  r1.height=height
  r1.input=r.curve
  if haskey(r,:name) r1.name=r.name end
  Discy(r1)
# set trueroots  to roots  of Discy  which do  not only  correspond to
# intersections of the curve with the chosen horizontal
  r1.trueroots=r.roots
  r1.verticallines=r.roots[1:degree(r.curveVerticalPart)]
  r1.roots=copy(r1.trueroots)
  append!(r1.roots,SeparateRoots(section,1000))
  r1
end

# curve --  an Mvp in x and y describing a curve in complex^2
function FundamentalGroup(c::Mvp;printlevel=0)
  SetPrintLevel(printlevel)
  r=PrepareCurve(c)
  Discy(r)
  GetDiscyRoots(r)
  if isempty(r.roots) return TrivialCase(r) end
  if !r.ismonic r=SearchHorizontal(r) end
  loops=convert_loops(LoopsAroundPunctures(r.roots))
  merge!(r.prop,pairs(loops))
  Loops(r)
  Zeros(r)
  if isempty(r.zeros[1]) return TrivialCase(r) end
  r.B=BraidMonoid(CoxSym(length(r.zeros[1])))
  r.monodromy=fill(r.B(),length(r.segments))
# return r
  for i in eachindex(r.segments) Segment(r, i) end
  Finish(r)
end

# Guarantees on LoopsAroundPunctures:
# For a set Z of zeroes and z in Z, let R(z):=1/2 dist(z,Z-z).
# The  input of  LoopsAroundPunctures is  a set  Z of approximate zeroes of
# r.discy such that for any z one of the zeroes is closer than R(z)/S where
# S is a global constant of the program (in practice we may take S=100).
# Let  d=inf_{z in  Z}(R(z)); we  return points  with denominator  10^-k or
# 10^-k<d/S'  (in practive we take S'=100) and  such that the distance of a
# segment to a zero of r.discy is guaranteed >= d-d/S'-d/S

# curve --  an Mvp in x and y "monic in x" describing a curve P(x,y)
function PrepareFundamentalGroup(curve, name)
  r=PrepareCurve(curve)
  Discy(r)
  r.name=name
  GetDiscyRoots(r)
  if isempty(r.roots) return TrivialCase(r) end
  if !r.ismonic r=SearchHorizontal(r) end
  merge!(r.prop,pairs(convert_loops(LoopsAroundPunctures(r.roots))))
  Loops(r)
  Zeros(r)
  if isempty(r.zeros[1]) return TrivialCase(r) end
  fname=r.name*".tmp"
  for f in [:curve,:discy,:ismonic,:roots,:loops,:segments,:points,:zeros]
    cut(string(r.name,".$f=",getproperty(r,f),"\n\n");before="+*",file=fname)
  end
  println("     ----------------------------------")
  println("Data saved in ", fname)
  println("You can now compute segments 1 to ", length(r.segments))
  println("in different GAP sessions by doing in each of them:")
  println("    ", r.name, "=rec(name=\"", r.name, "\")\n")
  println("    VKCURVE.Segments(",r.name,",[1..",length(r.segments),"])\n")
  println("(or some other range depending on the session)")
  println("Then when all files ",r.name,".xx have been computed finish by")
  println("    ",r.name,"=rec(name=\"", r.name, "\")\n")
  println("    FinishFundamentalGroup(", r.name, ")\n")
end

function FinishFundamentalGroup(r)
  read(r.name*".tmp")
  if !haskey(r, :monodromy) r.monodromy=[] end
  r[:B]=BraidMonoid(CoxSym(length(r.zeros[1])))
  for i in 1:length(r.segments)
    read(string(r.name,".",i))
    if isnothing(r.monodromy[i]) println("***** ", i, " missing ****") end
    r.monodromy[i]=r[:B](r.monodromy[i]...)
  end
  Finish(r)
end

#---------------------- root-finding ----------------------------------
"""
`NewtonRoot(p::Pol,initial,precision;showall=false,show=false,lim=800)`

Here `p` is a complex polynomial. The function computes an approximation to a
root of `p`, guaranteed of distance closer than `precision` to
an  actual root. The first approximation used is `initial`. If `initial` is
in  the  attraction  basin  of  a  root  of  `p`,  the  one approximated. A
possibility  is that  the Newton  method starting  from `initial`  does not
converge  (the  number  of  iterations  after  which  this  is  decided  is
controlled  by  `lim`);  then  the  function returns `nothing`.
Otherwise  the function returns  a pair: the  approximation found, and an
upper  bound of the distance between that approximation and an actual root.
The point of returning
an upper bound is that it is usually better than the asked-for `precision`.
For the precision estimate a good reference is cite{HSS01}.

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
function NewtonRoot(p::Pol,z,precision;showall=VKCURVE[:showallnewton],
                          show=VKCURVE[:showNewton],lim=VKCURVE[:NewtonLim])
  deriv=derivative(p)
  for cnt in 1:lim
    a=p(z)
    b=deriv(z)
    c=iszero(b) ? a : a/b
    err=abs(c)
    if iszero(err) err=(precision/100)/(degree(p)+1) end
    if err>precision err=precision end
    z-=c
    if showall println(cnt,": ",z) end
    if err<=(precision/100)/(degree(p)+1)
      if show print(cnt) end
      return (z,err)
    end
  end
  if show
    println("\n****** Non-Convergent Newton after ", lim," iterations ******")
    @show p,z,precision
  end
end

"""
'SeparateRootsInitialGuess(p, v, safety)'

Here  `p` is a complex  polynomial, and `v` is  a list of approximations to
roots  of `p` which should lie in different attraction basins for Newton' s
method.  The  result  is  a  list  `l`  of  complex  rationals representing
approximations  to the  roots of  `p` (each  element of  `l` is the root in
whose attraction basin the corresponding element of `v` lies), such that if
`d`  is the minimum distance  between two elements of  `l`, then there is a
root  of `p` within radius  `d/(2*safety)` of any element  of `l`. When the
elements  of  `v`  do  not  lie  in  different  attraction basins (which is
necessarily the case if `p` has multiple roots), 'false' is returned.

```julia-repl
julia> p=Pol([1,0,1])
Pol{Int64}: x²+1

julia> SeparateRootsInitialGuess(p,[1+im,1-im],10^5)
2-element Vector{ComplexF64}:
 8.463737877036721e-23 + 1.0im
 8.463737877036721e-23 - 1.0im

julia> SeparateRootsInitialGuess(p,[1+im,2+im],1000)
    # 1+im and 2+im in same attraction basins
```
"""
function SeparateRootsInitialGuess(p, v, safety)
  if degree(p)==1 return [-p[0]/p[1]] end
  radv=Dispersal(v)[1]/safety/2
  res=map(e->NewtonRoot(p,e,radv),v)
  if !any(isnothing,res) && Dispersal(v)[1]/2>=maximum(last.(res))
    return first.(res)
  end
  @show p,v,safety
  println("dispersal required=",Dispersal(first.(res))[1]/safety/2)
  println("obtained=",maximum(last.(res)))
  error()
end

"""
'SeparateRoots(<p>, <safety>)'

Here  `p` is  a complex  polynomial. The  result is  a list  `l` of complex
numbers  representing approximations to the roots  of `p`, such that if `d`
is  the minimum distance between two elements  of `l`, then there is a root
of  `p` within  radius `d/(2*safety)`  of any  element of  `l`. This is not
possible when `p` has multiple roots, in which case `nothing` is returned.

```julia-repl
julia> @Pol q
Pol{Int64}: q

julia> SeparateRoots(q^2+1,100)
2-element Vector{ComplexF64}:
  2.3541814200656927e-43 + 1.0im
 -2.3541814200656927e-43 - 1.0im

julia> SeparateRoots((q-1)^2,100)

julia> SeparateRoots(q^3-1,100)
3-element Vector{ComplexF64}:
 -0.5 - 0.8660254037844386im
  1.0 - 1.232595164407831e-32im
 -0.5 + 0.8660254037844387im
```
"""
function SeparateRoots(p,safety)
  subtractroot(p,r)=divrem(p,Pol([-r,1]))[1]
  if p isa Mvp p=Pol(p) end
  if degree(p)<1 return empty(p.c)
  elseif degree(p)==1 return [-p[0]/p[1]]
  end
  p/=p[end]
  e=Complex{Float64}(E(7))
  v=nothing
  cnt = 0
  while isnothing(v) && cnt<2*(degree(p)+1)
    if VKCURVE[:showNewton] && cnt>0
      println("****** ", cnt, " guess failed for p degree ", degree(p))
    end
    v=NewtonRoot(p,e,safety*10.0^-(degree(p)+4))
    e*=5/4*Complex{Float64}(E(2*(degree(p)+1)))
    cnt+=1
  end
  if cnt>=2*(degree(p)+1) error("no good initial guess") end
  v=[v[1]]
  append!(v,SeparateRoots(subtractroot(p,v[1]), safety))
  safety==0 ? v : SeparateRootsInitialGuess(p, v, safety)
end

"""
'FindRoots(<p>, <approx>)'

<p>  should be a univariate 'Mvp'  with cyclotomic or 'Complex' rational or
decimal  coefficients or  a list  of cyclotomics  or 'Complex' rationals or
decimals  which represents  the coefficients  of a  complex polynomial. The
function  returns  'Complex'  rational  approximations  to the roots of <p>
which  are  better  than  <approx>  (a  positive rational). Contrary to the
functions  'SeparateRoots', etc... described in  the previous chapter, this
function handles quite well polynomials with multiple roots. We rely on the
algorithms explained in detail in cite{HSS01}.

julia> FindRoots((Pol()-1.0)^5,1/5000)
5-element Vector{ComplexF64}:
 1.0009973168670234 - 6.753516026574732e-9im
 1.0004572203796391 - 0.00033436001135679156im
  0.999029224439348 + 5.075797152907413e-12im
  0.999723670720127 - 0.0008487747754577878im
 1.0007950023584915 - 0.0005779801979057327im

julia> FindRoots(Pol()^3-1,10^-5)
3-element Vector{ComplexF64}:
 -0.5 - 0.8660254037844386im
  1.0 - 1.232595164407831e-32im
 -0.5 + 0.8660254037844387im

julia> map(x->x^3,ans)
3-element Vector{ComplexF64}:
 0.9999999999999998 - 1.1102230246251565e-16im
                1.0 - 3.697785493223493e-32im
 1.0000000000000002 - 1.1102230246251565e-16im
"""
function FindRoots(p,prec)
  subtractroot(p,r)=divrem(p,Pol([-r,1]))[1]
  if degree(p)<1 return empty(p.c)
  elseif degree(p)==1 return [-p[0]/p[1]]
  end
  e=Complex{Float64}(E(7))
  v=nothing
  while isnothing(v)
    v=NewtonRoot(p,e,10.0^(-degree(p)-1))
    e*=Complex{Float64}(E(degree(p)+1))
  end
   v=vcat([v[1]],FindRoots(subtractroot(p,v[1]),prec))
   map(e->NewtonRoot(p,e,prec)[1],v)
end

#------------------ Loops --------------------------------------------
# Ordonne une liste de points trigonometriquement autour d'un centre
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

# Input: (list of complex numbers,complex number)
# Output: sublist of "neighbours" of the second input,
#   x and y are neighbours iff no z is in the disk of diameter [x,y]
function neighbours(l, center)
  function isneighbour(y)
    d=abs2(y-center)
    for z in l
      if z!=y && abs2(y-z)+abs2(z-center)<=d return false end
    end
    return true
  end
  l=filter(y->y!=center,l)
  l=filter(isneighbour, l)
  cycorder(l, center)
end

# value at z of an equation of the line (x,y)
function lineq(x, y, z)
  if real(x)≈real(y)
    if imag(x)≈imag(y) error("Undefined line\n")
    else return real(z)-real(x)
    end
  else
    return (imag(y)-imag(x))*(real(z)-real(x))/(real(y)-real(x))+imag(x)-imag(z)
  end
end

function mediatrix(x, y)
  if x≈y error("Undefined mediatrix") end
  (x+y)/2 .+[im,-im]*(x-y)
end

crossing(v1,v2)=crossing(v1...,v2...)

# Computes the intersecting point of two lines
# returns false if the lines are parallel or a pair is a single
function crossing(x1,x2,y1,y2)
  if x1≈x2 || y1≈y2 return nothing end
  if !(real(x1)≈real(x2))
    lambdax=(imag(x1)-imag(x2))/(real(x1)-real(x2))
    mux=-lambdax*real(x1)+imag(x1)
    if !(real(y1)≈real(y2))
      lambday=(imag(y1)-imag(y2))/(real(y1)-real(y2))
      muy=-lambday*real(y1)+imag(y1)
      if 1+lambdax≈1+lambday return nothing end
      resr=(muy-mux)/(lambdax-lambday)
      resi=lambdax*resr+mux
      res=Complex(resr, resi)
    else
      E3=Complex{Float64}(E(3))
      res=crossing(E3*x1, E3*x2, E3*y1, E3*y2)
      if isnothing(res) return nothing end
      res/=E3
    end
  else
    res=crossing(im*x1, im*x2, im*y1, im*y2)
    if isnothing(res) return nothing end
    res/=im
  end
  res
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

# y must be an element of ys
function boundpaths(ys, sy, path, y)
  if !haskey(y, :path)
    y[:path]=vcat(path, [y[:y]])
    for z in y[:lovers] boundpaths(ys, sy, y[:path], sy(z)) end
  end
end

# eliminates trivial segments and contracts pairs [a,b],[b,a]
function shrink(l)
  k=findfirst(i->l[i]==l[i+1],1:length(l)-1)
  if !isnothing(k) return shrink(vcat(l[1:k],l[k+2:end])) end
  k=findfirst(i->l[i]==l[i+2],1:length(l)-2)
  if !isnothing(k) return shrink(vcat(l[1:k],l[k+3:end])) end
  l
end

# converts old format (list of loops) to the new format :
#  .points   : set of all endpoints of all segments
#  .segments : set of all segments used, where endpoints are indexed
#              as in points
#     .loops : list of sequence of numbers of used segments
"""
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
  points=sort(unique(vcat(ll...)),by=x->(imag(x),real(x)))
  np(p)=findfirst(==(p),points)
  loops=map(l->map(i->np.(l[i-1:i]),2:length(l)),ll)
  segments=sort(unique(sort.(vcat(loops...))))
  loops=map(loops)do l
    map(l)do seg
     seg[1]<seg[2] ? findfirst(==(seg),segments) :
                    -findfirst(==(reverse(seg)),segments)
    end
  end
  (;points, segments, loops)
end

"""
'LoopsAroundPunctures(points)'

`points`  should be complex numbers. The function computes piecewise-linear
loops representing generators of the fundamental group of `ℂ -{points}`.

```julia-repl
julia> LoopsAroundPunctures([0])
1-element Vector{Vector{Complex{Int64}}}:
 [1 + 0im, 0 + 1im, -1 + 0im, 0 - 1im, 1 + 0im]
```
"""
function LoopsAroundPunctures(originalroots)
  roots=originalroots
  n=length(roots)
  if n==1 return [roots[1].+[1,im,-1,-im,1]] end
  average=sum(roots)/n
  sort!(roots, by=x->abs2(x-average))
  ys=map(x->Dict{Symbol, Any}(:y=>x), roots)
  sy(y)=ys[findfirst(==(y),roots)]
  for y in ys
    y[:neighbours]=neighbours(roots, y[:y])
    y[:friends]=[y[:y]]
    y[:lovers]=empty(y[:friends])
  end
  if VKCURVE[:showLoops] println("neighbours computed") end
  for y in ys
    for z in y[:neighbours]
      if !(z in y[:friends])
        push!(y[:lovers], z)
        push!(sy(z)[:lovers], y[:y])
        newfriends=vcat(y[:friends], sy(z)[:friends])
        for t in y[:friends] sy(t)[:friends]=newfriends end
        for t in sy(z)[:friends] sy(t)[:friends]=newfriends end
      end
    end
  end
  for y in ys sort!(y[:neighbours],by=z->abs2(y[:y]-z)) end
  minr,maxr=extrema(map(y->real(y[:y]), ys))
  mini,maxi=extrema(map(y->imag(y[:y]), ys))
# To avoid trouble with points on the border of the convex hull,
# we make a box around all the points;
  box=[Complex(minr-2, mini-2), Complex(minr-2, maxi+2),
       Complex(maxr+2, mini-2), Complex(maxr+2, maxi+2),
       Complex((maxr+minr)/2, mini-(maxr-minr)/2-2),
       Complex((maxr+minr)/2, maxi+(maxr-minr)/2+2),
       Complex(minr-(maxi-mini)/2-2, (maxi+mini)/2),
       Complex(maxr+(maxi-mini)/2+2, (maxi+mini)/2)]
  for y in ys
    y[:cycorder]=cycorder(vcat(filter(z->z!=y[:y],roots),box),y[:y])
    k=findfirst(==(y[:neighbours][1]),y[:cycorder])
    y[:cycorder]=circshift(y[:cycorder],1-k)
    push!(y[:cycorder], y[:cycorder][1])
    y[:circle]=Complex{Float64}[(y[:y]+y[:neighbours][1])/2]
    y[:witness]=Complex{Float64}[y[:neighbours][1]]
    for z in y[:cycorder][2:end]
      cut=detectsleftcrossing(y[:circle], y[:witness], y[:y], z)
      if true in cut
        k=findfirst(cut)
        y[:circle]=y[:circle][1:k]
        y[:witness]=y[:witness][1:k]
      end
      k=length(y[:circle])
      newcirc=crossing(mediatrix(y[:y],y[:witness][k]),mediatrix(y[:y], z))
      if !isnothing(newcirc)
        push!(y[:circle], newcirc)
        push!(y[:witness], z)
      end
      if z in y[:lovers]
        push!(y[:circle], (y[:y]+z)/2)
        push!(y[:witness], z)
      end
    end
  end
  if VKCURVE[:showLoops] println("circles computed") end
  boundpaths(ys, sy, [], ys[1])
  for y in ys
    k=length(y[:path])
    if k>1
      circleorigin=(y[:y]+y[:path][k-1])/2
      k=findfirst(≈(circleorigin),y[:circle])
      y[:circle]=circshift(y[:circle],1-k)
    end
  end
  for y in ys
    k=length(y[:path])
    y[:handle]=Vector{Complex{Float64}}(vcat(map(1:k-1)do i
      l=sy(y[:path][i])[:circle]
      l[1:findfirst(≈((y[:path][i]+y[:path][i+1])/2),l)]
     end...))
    y[:loop]=vcat(y[:handle], y[:circle], reverse(y[:handle]))
  end
  sort!(ys,by=y->findfirst(==(y[:y]),originalroots))
  for y in ys
    y[:loop]=map(x->round(x;sigdigits=8),y[:loop])
    y[:loop]=shrink(y[:loop])
  end
  map(y->y[:loop], ys)
end

function segs(r,v)
  rp=Float64[]
  ip=Float64[]
  for seg in v
    s=r[:points][seg<0 ? reverse(r[:segments][-seg]) : r[:segments][seg]]
    append!(rp,real.(s))
    append!(ip,imag.(s))
  end
  rp,ip
end

function loops(r,v)
  rp=Float64[]
  ip=Float64[]
  for l in r[:loops][v]
    nrp,nip=segs(r,l)
    append!(rp,nrp)
    append!(ip,nip)
  end
  rp,ip
end
#-------------------- ApproxMonodromy ----------------------------
# for each point of a find closest point in b
# Complain if the result is not a bijection between a and b of if
# the distance between an a and the correponding b is bigger than 1/10
# of minimum distance between two b's
function fit(a, b)
  dm=map(p->findmin(abs.(b.-p)),a)
  monodromyError=maximum(first.(dm))
# println("# Monodromy error==",monodromyError)
  if monodromyError>Dispersal(b)[1]/10 error("monodromy error too big") end
  pos=last.(dm)
  if sort(pos)!=1:length(pos) error("monodromy cannot find perm") end
  return b[pos]
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
'ApproxFollowMonodromy(<r>,<segno>,<pr>)'

This function  computes an approximation  of the monodromy braid  of the
solution in `x`  of an equation `P(x,y)=0` along  a segment `[y_0,y_1]`.
It is called  by 'FundamentalGroup', once for each of  the segments. The
first  argument is  a  global record,  similar to  the  one produced  by
'FundamentalGroup'  (see the  documentation of  this function)  but only
containing intermediate information. The second argument is the position
of the segment in 'r.segments'. The  third argument is a print function,
determined  by the  printlevel set  by the  user (typically,  by calling
'FundamentalGroup' with a second argument).

Contrary to 'FollowMonodromy',  'ApproxFollowMonodromy' does not control
the approximations; it just uses a  heuristic for how much to move along
the segment  between linear braid  computations, and this  heuristic may
possibly fail. However,  we have not yet found an  example for which the
result is actually incorrect, and thus the existence is justified by the
fact that  for some difficult  computations, it is sometimes  many times
faster  than 'FollowMonodromy'.  We illustrate  its typical  output when
<printlevel> is 2.

|   VKCURVE.monodromyApprox:=true;
julia-rep1```
julia> FundamentalGroup((x+3*y)*(x+y-1)*(x-y);printlevel=2)

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
the 3  roots of the current  polynomial `f(x,y_0)` if we  are looking at
the point `y_0` of the segment.  Then, which segment we are dealing with
(here the  15th of  16 in  all). Then the  minimum distance  between two
roots of  `f(x,y_0)` (used in our  heuristic). Then the current  step in
fractions of the length of the segment  we are looking at, and the total
fraction of the segment we have  done. Finally, the decimal logarithm of
the absolute  value of the discriminant  at the current point  (used in
the heuristic). Finally, an indication if the heuristic predicts that we
should  halve the  step  ('***rejected')  or that  we  may double  it
('***up').

The function returns an element of the ambient braid group 'r.B'.
"""
function ApproxFollowMonodromy(r,segno,pr)
  if VKCURVE[:showInsideSegments] ipr=print
  else ipr=function(x...)end
  end
  p,q=r.segments[segno]
  res=r.B()
  prevzeros=r.zeros[p]
  n=length(prevzeros)
  if n==1 return r.B() end
  mindm=Dispersal(prevzeros)[1]
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
    nextzeros=SeparateRootsInitialGuess(P, prevzeros, 100)
    if isnothing(nextzeros) ||
       (1+maximum(abs.(nextzeros-prevzeros))≈1 && step>1//16)
      rejected=true
    else
      dm=map(i->minimum(abs.(prevzeros[i]-prevzeros[j] for j in 1:n if j!=i)),
                                                                          1:n)
      mdm=minimum(dm)
      if step<1 ipr("<$segno/",length(r.segments),">mindist=",approx(mdm),
         " step=$step total=$total logdisc=",normdisc(r.discyFactored,next))
      end
      dn=abs.(prevzeros-nextzeros)
      rejected=any(dm.<VKCURVE[:AdaptivityFactor].*dn)
      if !rejected && mdm<mindm mindm=mdm end
    end
    if rejected
      step/=2
      ipr(" ***rejected\n")
      if step<minstep minstep=step end
    else
      total+=step
      if all(dm.>2 .*VKCURVE[:AdaptivityFactor] .*dn) && total+step!=1
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
  pr("# Minimal distance=", approx(mindm), "\n")
  pr("# Minimal step=", minstep, "=", approx(v*minstep), "\n")
  pr("# Adaptivity=", VKCURVE[:AdaptivityFactor], "\n")
  res*LBraidToWord(prevzeros,fit(nextzeros,r.zeros[q]),r.B)
end
#-------------------- Monodromy ----------------------------
# ceil(-log2(p)) for 0<p<1
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

# truncated iteration of the Newton method
function mynewton(p,z)
  a=p(z)
  b=derivative(p)(z)
  if iszero(b) c=a
    print("NewtonError\n")
  else c=a/b
  end
  err=degree(p)*abs(c)
  if err==0 prec=1
  else prec=max(0,ceil(Int,-log10(err)))+2
  end
  round(z-c;digits=prec)
end

# for each point of a find closest point in b
function myfit(a, b)
  d=length(a)
  dist=fill(0.0,d,d)
  for k in 1:d, l in k+1:d
    dist[k,k]=dist[l,k]=dist[k,l]=abs2(a[k]-a[l])
  end
  dist[d,d]=dist[d,d-1]
  R=map(k->minimum(dist[k,:])/4,1:d)
  map(k->only(filter(i->abs2(i-a[k])<R[k], b)),1:d)
end

Base.setindex!(p::Pol{T},x::T,i::Integer) where T=p.c[i+1-p.v]=x

# Sturm(pp,time)
# if polynomial pp is positive  at time
# returns some rational number t such that
#    time<t<=1  and  pp  is positive on [time,t]
# otherwise returns 0
# [third input and second output is an adaptive factor to
#  accelerate the computation]
function Sturm(pp::Pol, time, adapt::Integer;pr=print)
  q=Pol()
  pol=pp((1-q)*time+q)
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
  pr(m)
  if m==adapt && adapt>0
    if pol(3t//2)>0
      if pol(2t)>0 res=[(1-2t)*time+2t, adapt-1]
      else res=[(1-3t//2)*time+3t//2, adapt-1]
      end
    else res=[(1-t)*time+t, adapt]
    end
  else res=[(1-t)*time+t, m]
  end
  res
end

Base.real(p::Pol)=Pol(real.(p.c),p.v)
Base.imag(p::Pol)=Pol(imag.(p.c),p.v)

"""
'FollowMonodromy(<r>,<segno>,<print>)'
This function computes the monodromy braid  of the solution in `x` of an
equation  `P(x,y)=0`  along  a  segment `[y_0,y_1]`.  It  is  called  by
'FundamentalGroup', once for each of the segments. The first argument is
a global record, similar to  the one produced by 'FundamentalGroup' (see
the  documentation of  this function)  but only  containing intermediate
information.  The second  argument is  the  position of  the segment  in
'r.segments'. The third argument is  a print function, determined by the
printlevel  set by  the user  (typically, by  calling 'FundamentalGroup'
with a second argument).

The function returns an element of the ambient braid group 'r.B'.

This function has no reason to be  called directly by the user, so we do
not illustrate its  behavior. Instead, we explain what  is displayed on
screen when the user sets the printlevel to `2`.

What is quoted below is an excerpt of what is displayed on screen
during the execution of
|    gap>  FundamentalGroup((x+3*y)*(x+y-1)*(x-y),2);

    <1/16>    1 time=          0   ?2?1?3
    <1/16>    2 time=      0.125   R2. ?3
    <1/16>    3 time=    0.28125   R2. ?2
    <1/16>    4 time=   0.453125   ?2R1?2
    <1/16>    5 time=   0.578125   R1. ?2
    ======================================
    =    Nontrivial braiding = 2         =
    ======================================
    <1/16>    6 time=   0.734375   R1. ?1
    <1/16>    7 time=    0.84375   . ?0.
    <1/16>    8 time=   0.859375   ?1R0?1
    # The following braid was computed by FollowMonodromy in 8 steps.
    monodromy[1]:=B(2);
    # segment 1/16 Time=0.1sec|

'FollowMonodromy' computes  its results by subdividing  the segment into
smaller  subsegments  on which  the  approximations  are controlled.  It
starts at one  end and moves subsegment after subsegment.  A new line is
displayed at each step.

The  first column  indicates which  segment is  studied. In  the example
above, the function  is computing the monodromy along  the first segment
(out  of  `16`).  This  gives  a  rough  indication  of  the  time  left
before  completion of  the total  procedure.  The second  column is  the
number of  iterations so  far (number of  subsegments). In  our example,
'FollowMonodromy'  had to  cut the  segment into  `8` subsegments.  Each
subsegment has its own length. The cumulative length at a given step, as
a  fraction of  the  total length  of the  segment,  is displayed  after
'time='.  This  gives  a  rough  indication  of  the  time  left  before
completion  of the  computation of  the monodromy  of this  segment. The
segment is completed when this fraction reaches `1`.

The last column has to do with the piecewise-linear approximation of the
geometric monodromy  braid. It is  subdivided into sub-columns  for each
string. In  the example above,  there are  three strings. At  each step,
some strings are fixed (they are  indicated by '. ' in the corresponding
column). A symbol like 'R5' or '?3' indicates that the string is moving.
The exact meaning of the symbol has to do with the complexity of certain
sub-computations.

As  some strings  are moving,  it  happens that  their real  projections
cross. When such a crossing occurs, it is detected and the corresponding
element of `B_n` is displayed on screen ('Nontrivial braiding ='...) The
monodromy braid is the product of these elements of `B_n`, multiplied in
the order in which they occur.
"""
# Exact computation of the monodromy braid along a segment
# r: global VKCURVE record
# seg: segment number
# sprint: Print function (to screen, to file, or none)
function FollowMonodromy(r,seg,sPrint)
  if VKCURVE[:showInsideSegments] iPrint=print
  else iPrint=function(arg...) end
  end
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
  time=big(0)
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
    iPrint(" time=",lpad(time,11),"   ")
    for k in 1:d, l in k+1:d
      dist[k,k]=dist[l,k]=dist[k,l]=abs2(big(v[k]-v[l]))
    end
    dist[d,d]=dist[d,d-1]
    for k in 1:d
      Rk=minimum(dist[k,:])/4
      z=v[k]
      if protected[k]>time && Rk>=RR[k]
        iPrint(". ")
      elseif protected[k]>time
        if adapt[k]+2<maximum(adapt) Rk/=2 end
        iPrint("R")
        s,adapt[k]=Sturm(Rk*protdpdx[k]-protp[k], time, adapt[k])
        if s>time protected[k]=binlowevalf(s,time)
        else iPrint("How bizarre...")
#         @show Rk,protdpdx[k],protp[k]
#         @show Rk*protdpdx[k]-protp[k], time, adapt[k]
        end
        RR[k]=Rk
      else
        iPrint("?")
        cptz=Pol(pt(;x=z))
#       @show pt
#       @show z
#       @show cptz
        protp[k]=d^2*(real(cptz)^2+imag(cptz)^2)
        cdpdxtz=Pol(dpdxt(;x=z))
        protdpdx[k]=real(cdpdxtz)^2+imag(cdpdxtz)^2
        s,adapt[k]=Sturm(Rk*protdpdx[k]-protp[k], time, adapt[k])
        if s>time protected[k]=binlowevalf(s,time)
        else error("Something's wrong...")
#         @show R[k],protdpdx[k],protp[k]
#         @show R[k]*protdpdx[k]-protp[k], time, adapt[k]
        end
        RR[k]=Rk
      end
    end
    allowed=minimum(protected)
    time=allowed
    py=Pol(p(;y=r.points[a]*(1-time)+r.points[b]*time))
    iPrint("\n")
    newv=map(1:d)do k
      if protected[k]>allowed v[k]
      else mynewton(py,v[k])
      end
    end
    res*=LBraidToWord(v, newv, B)
    v=newv
    if time==1 break end
  end
  sPrint("# The following braid was computed by FollowMonodromy in $steps steps.\n")
  res*LBraidToWord(v, myfit(v, r.zeros[b]), B)
end

#------------------- Compute PLBraid ----------------------------------
# Deals with "star" linear braids, those with associated permutation w_0
function starbraid(y, offset, B)
  n=length(y)
  if n==1 return B() end
  k=argmin(y)
  B((k:n-1)+offset...)*starbraid(y[setdiff(1:n,[k])],offset,B)/
  B((n+1-k:n-1)+offset...)
end

# In case two points have the same real projection, we use
# a "lexicographical" desingularization by "infinitesimal rotation"
function desingularized(v1, v2)
  n=length(v1)
  tan=1
  for k in 1:n, l in k+1:n
    rv,iv=reim(v1[k]-v1[l])
    if abs(iv*rv)>10^-8 tan=min(tan,abs(rv/iv)) end
    rv,iv=reim(v2[k]-v2[l])
    if abs(iv*rv)>10^-8 tan=min(tan,abs(rv/iv)) end
  end
  [v1, v2].*Complex(1,-tan/2)
end

function setapprox(l)
  if isempty(l) return l end
  l=sort(l)
  prev=1
  for next in 2:length(l)
    if isapprox(l[prev],l[next];rtol=10^-8,atol=10^-10) continue end
    prev+=1
    l[prev]=l[next]
  end
  l[1:prev]
end

"""
'LBraidToWord(v1,v2,B)'

This function converts  the linear braid joining the points in `v1` to the
corresponding ones in `v2` into an element of the braid group.

```julia-repl
julia> B=BraidMonoid(CoxSym(3))
BraidMonoid(𝔖 ₃)

julia> LBraidToWord([1+im,2+im,3+im],[2+im,1+2im,4-6im],B)
1
```

The lists `v1` and `v2` must have the same length, say `n`. Then `B` should
be  `BraidMonoid(CoxSym(n))`, the braid group  on `n` strings. The elements
of  `v1` (resp. `v2`)  should be `n`  distinct complex rational numbers. We
use the Brieskorn basepoint, namely the contractible set `C+iV_ℝ` where `C`
is  a real chamber; therefore the endpoints  need not be equal. The strings
defined  by `v1` and `v2` should be  non-crossing. When the numbers in `v1`
(resp.  `v2`)  have  distinct  real  parts,  the  real picture of the braid
defines a unique element of `B`. When some real parts are equal, we apply a
lexicographical  desingularization, corresponding to a rotation of `v1` and
`v2` by an arbitrary small positive angle.
"""
# two printlevel control fields: VKCURVE.showSingularProj
#				 VKCURVE.showBraiding
# Deals with linear braids
# 1) singular real projections are identified
# 2) calls starbraid for each
function LBraidToWord(v1, v2, B)
  n=length(v1)
  x1,y1=reim(v1)
  x2,y2=reim(v2)
  if length(setapprox(x1))<length(x1) || length(setapprox(x2))<length(x2)
    if VKCURVE[:showSingularProj]
      println("WARNING: singular projection(resolved)")
    end
    return LBraidToWord(desingularized(v1, v2)..., B)
  end
  q=sortPerm(x1)
  crit=Float64[]
  for i in 1:n-1, j in i+1:n
    if x2[i^q]>x2[j^q]
      push!(crit,(x1[i^q]-x1[j^q])/((x2[j^q]-x1[j^q]+x1[i^q])-x2[i^q]))
    end
  end
  tcrit=setapprox(crit)
  res=B()
  u=0
  for t in tcrit
    xt=map(k->x1[k]+t*(x2[k]-x1[k]),1:n)
    yt=map(k->y1[k]+t*(y2[k]-y1[k]),1:n)
    ut=(u+t)/2
    xut=map(k->x1[k]+ut*(x2[k]-x1[k]),1:n)
    put=inv(sortPerm(xut))
    xt=xt^put
    yt=yt^put
    xcrit=setapprox(xt)
    for x in xcrit
      posx=findfirst(≈(x;atol=10^-10),xt)
      nx=count(≈(x;atol=10^-10),xt)
      res*=starbraid(yt[posx:posx+nx-1], posx-1, B)
    end
    u=t
  end
  if VKCURVE[:showBraiding]
   if !isempty(tcrit)
      if VKCURVE[:showInsideSegments]
        println("======================================")
        println("==    Nontrivial braiding ",rpad(res,10),"==")
        println("======================================")
#       print("v1:=[",join(gg.(v1),","),"];\n")
#       print("v2:=[",join(gg.(v2),","),"];\n")
      else println("==    Nontrivial braiding ",rpad(res,10),"==")
      end
    end
  end
  return res
end
#----------------------- Presentation -------------------------------
# Hurwitz action of the braid b on the list l of group elements
function hurwitz(b,l)
  l=copy(l)
  for s in word(b)
    if s>0 l[s:s+1]=[l[s+1],l[s]^l[s+1]]
    else l[-s:-s+1]=[l[-s]*l[-s+1]*inv(l[-s]),l[-s]]
    end
  end
  return l
end

# Implements the action of B_n on F_n
#   Input:  element of B_n, ambient free group F_m
#           (with m >= n; when m>n, Hurwitz
#            action is extended trivially to the extra generators)
#   Output: automorphism of F_n
"""
'BnActsOnFn(<braid b>,<Free group F>)'

This function  implements the Hurwitz action  of the braid group  on `n`
strings  on  the  free  group  on `n`  generators,  where  the  standard
generator  `σ_i` of  `B_n`  fixes  the generators  `f_1,…,f_n`,
except `f_i` which is mapped to  `f_{i+1}` and `f_{i+1}` which is mapped
to `f_{i+1}^{-1}f_if_{i+1}`.

```julia-repl
julia> B=BraidMonoid(CoxSym(3))
BraidMonoid(𝔖 ₃)

julia> b=B(1)
1

julia> BnActsOnFn(b,FpGroup(:a,:b,:c))
Aut(FreeGroup(a,b,c);AbsWord[a, b, c]↦ AbsWord[b, b⁻¹ab, c]

julia> BnActsOnFn(b^2,FpGroup(:a,:b,:c))
Aut(FreeGroup(a,b,c);AbsWord[a, b, c]↦ AbsWord[b⁻¹ab, b⁻¹a⁻¹bab, c]
```
The second input is the free group on `n` generators. The first input is
an  element  of  the  braid  group  on  `n`  strings,  in  its  CHEVIE
implementation.
"""
BnActsOnFn(b,F)=Hom(F,F,hurwitz(b,gens(F)))

"""
'VKQuotient(braids)'

The input `braid` is a list of braids `b₁,…,b_d`, living in the braid group
on `n` strings. Each `bᵢ` defines by Hurwitz action an automorphism `φᵢ` of
the  free group `Fₙ`. The function return the group defined by the abstract
presentation: ``< f₁,…,fₙ ∣ ∀ i,j φᵢ(fⱼ)=fⱼ > ``

|    gap> B:=Braid(CoxeterGroupSymmetricGroup(3));
    function ( arg ) ... end
    gap> b1:=B(1)^3; b2:=B(2);
    1.1.1
    2
    gap> g:=VKQuotient([b1,b2]);
    Group( f.1, f.2, f.3 )
    gap>  last.relators;
    [ f.2^-1*f.1^-1*f.2*f.1*f.2*f.1^-1, IdWord,
      f.2^-1*f.1^-1*f.2^-1*f.1*f.2*f.1, f.3*f.2^-1, IdWord, f.3^-1*f.2 ]
    gap> p:=PresentationFpGroup(g);display_balanced(p);
    << presentation with 3 gens and 4 rels of total length 16 >>
    1: c=b
    2: b=c
    3: bab=aba
    4: aba=bab
    gap> SimplifyPresentation(p);display_balanced(p);
    #I  there are 2 generators and 1 relator of total length 6
    1: bab=aba|
"""
function VKQuotient(braids)
  # get the true monodromy braids and the Hurwitz action basic data
  n=braids[1].M.W.n
  F=FpGroup(Symbol.('a'.+(0:n-1))...)
  rels=AbsWord[]
  f=gens(F)
  for b in braids append!(rels,map((a,b)->a*inv(b),hurwitz(b,f),f)) end
  F/rels
end

# A variant of the previous function.
# See arXiv:math.GR/0301327 for more mathematical details.
# Input: global VKCURVE record
# Output: the quotient, encoded as an FpGroup
############################################################
# Printing controlled by VKCURVE.showAction
function DBVKQuotient(r)
  # get the true monodromy braids and the Hurwitz action basic data
  n=r.braids[1].M.W.n
  F=FpGroup(Symbol.('a'.+(0:n+length(r.verticallines)-1))...)
# above the basepoint for the loops, locate the position of the string
# corresponding to the trivializing horizontal line
  bzero=r.zeros[r.basepoint]
  dist=abs2.(bzero-r.height)
  height=bzero[argmin(dist)]
  basestring=count(z->reim(z)<=reim(height),bzero)
  fbase=F(basestring)
  rels=AbsWord[]
  auts=map(b->BnActsOnFn(b, F),r.braids)
  for (i,aut) in enumerate(auts)
# Find an element conjugator such that aut(fbase)^inv(conjugator)=fbase
    ifbase=aut(fbase)
    conjugator=F.gens[1]*inv(F.gens[1])
    choices=vcat(map(f->[f,inv(f)], gens(F)[1:n])...)
    while true
      k=0
      while true
        k+=1
        if length(choices[k]*ifbase)<length(ifbase) break end
      end
      ifbase=choices[k]*ifbase*inv(choices[k])
      conjugator=choices[k]*conjugator
      if length(ifbase)==1 break end
    end
# Replacing aut by  correctaut:= Conj(conjugator)*aut
    conj=Hom(F, F, gens(F).^Ref(inv(conjugator)))
    correctaut=x->conj(aut(x))
    if i>length(r.verticallines)
      append!(rels, map(f->correctaut(f)*inv(f), gens(F)[1:n]))
    else
      g=F(i+n)
      append!(rels, map(f->correctaut(f)*g*inv(g*f), gens(F)[1:n]))
    end
  end
  push!(rels, fbase)
  F/rels
end

gg(x)="Complex(evalf(\"$(real(x))\"),evalf(\"$(imag(x))\"))"
VKCURVE[:showInsideSegments]=true
VKCURVE[:showBraiding]=true
VKCURVE[:showNewton]=true

data=Dict()

@Mvp x,y,z,t
d=discriminant(ComplexReflectionGroup(24))(x,y,z)
data[24]=d(;x=1,z=x)
d=discriminant(ComplexReflectionGroup(27))(x,y,z)
data[27]=d(;x=1,z=x)
d=discriminant(ComplexReflectionGroup(23))(x,y,z)
data[23]=d(;x=1,z=x)
d=discriminant(ComplexReflectionGroup(29))(x,y,z,t)
data[29]=d(;t=y+1,z=x)
d=discriminant(ComplexReflectionGroup(31))(x,y,z,t)
data[31]=d(;t=x+1,z=y)
data[34]=
95864732434895657396628326400//164799823*x*y^3-598949723065092000//
1478996726054382501274179923886253687929138281*x*y^7-
67840632073999787861633181671139840000*x^2-
7622790471072621273612030528032173587500421120000*y^2-273861000//
27158981660831329*x^2*y^4+37130333513291749382400//7130353846013*x^3*y-
13608525//50841945969352380915996169*x^4*y^2-2606867429323404078970327//
1675017448527954334139901265590107596081497211494411528*x^6+
3269025273548225517660538475128200000//390195840687434028022928452202401489*y^6
