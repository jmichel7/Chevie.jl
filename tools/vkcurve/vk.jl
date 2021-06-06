#------------------- utilities -------------------------------
# alternative to ComplexOps.Norm(x)
BigNorm(x)=abs(real(x))+abs(imag(x))

# another alternative
SmallNorm(x)=max(abs(real(x)),abs(imag(x)))

# For a non-zero rational x, returns k such that 10^k<|x|<=10^(k+1) 
function DecimalLog(x)
  if iszero(x) error("trying to take decimal log of 0") end
  d = div(denominator(x), numerator(x))
  if iszero(d)
    d=div(numerator(x), denominator(x))
    return length(string(abs(d)))-1
  else return -length(string(abs(d)))
  end
end

# a complex rational approximation of a (a Cyc or ??)
function ComplexRational(a)
  if a isa Cyc
      a = Complex(a)
      if !(IsRat(a[:r]))
          a[:r] = evalf(a[:r])
      end
      if !(IsRat(a[:i]))
          a[:i] = evalf(a[:i])
      end
  end
  return Complex(Rational(a[:r]), Rational(a[:i]))
end

"""
`Dispersal(v)`
`v`  is a list of complex numbers  representing points in the real plane.
The  result is a pair  whose first element is  the minimum distance between
two  elements of `v`, and the second is a pair of indices `[i,j]` such that
`v[i]`, `v[j]` achieves this minimum distance.

julia> Dispersal([1+im,0,1])
(1, [1,3])
"""
function Dispersal(v)
  l=combinations(eachindex(v),2)
  p=findmin(map(x->BigNorm(v[x[1]]-v[x[2]]),l))
  (p[1],l[p[2]])
end

# distance of z to segment [a,b]
function DistSeg(z,a,b)
  b-=a
  z-=a
  r=abs(b)
  z*=r/b
  real(z)<0 ? abs(z) : real(z)>r ? abs(z-b) : imag(z)>0 ? imag(z) : -imag(z)
end

InterpolatedPolynomial=function(x,y)
  t=copy(y)*1//1 # make sure coeffs are in a field
  a=map(1:length(x))do i
    for k in i-1:-1:1
      t[k]=(t[k+1]-t[k])/(x[i]-x[k])
      if isnothing(t[k]) error("cannot interpolate polynomial") end
    end
    t[1]
  end
  p=a[length(x)]
  for i in length(x)-1:-1:1
    p=p*(Pol()-x[i])+a[i]
  end
  p
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
  InterpolatedPolynomial(1:n,v)
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
:monodromyApprox=>true, ##########################
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
  if haskey(r,:presentation) Presentations.DisplayPresentation(r.presentation)
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
  merge!(r.prop,pairs(LoopsAroundPunctures(r.roots)))
# here we  have loops around the  'true' roots and around  the 'extra'
# roots Difference(r.roots,r.trueroots). We get rid of the extra loops
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
# m=-DecimalLog(Rational(r.dispersal)//100)
# r.points=map(y->evalf(y, m),r.points)
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
    zz=SeparateRoots(Pol(complexmvp(r.curve)(y=r.points[i])), 1000)
    if length(zz)>1
      m=Dispersal(zz)
#     m[1]=evalf(m[1], -(DecimalLog(m[1] // 1000)))
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
  r.discyFactored=[d]
  # r.discyFactored=Factors(r.discy)
end

complexmvp(p::Mvp)=Mvp(ModuleElt(m=>Complex{Float64}(c) for (m,c) in p.d))

function GetDiscyRoots(r)
  if VKCURVE[:showRoots] println("Computing roots of discriminant...") end
  if degree(r.curveVerticalPart)==0 r.roots=Complex{Float64}[]
  else r.roots=SeparateRoots(complexpol(r.curveVerticalPart), 1000)
  end
  append!(r.roots,vcat(map(p->SeparateRoots(complexpol(p),1000),
                           r.discyFactored)...))
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
# simplify(r.presentation)
  r
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
  append!(r1.roots,SeparateRoots(complexpol(section),1000))
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
  Loops(r)
  Zeros(r)
  if isempty(r.zeros[1]) return TrivialCase(r) end
  r.B=BraidMonoid(CoxSym(length(r.zeros[1])))
  r.monodromy=fill(r.B(),length(r.segments))
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
  Loops(r)
  Zeros(r)
  if isempty(r.zeros[1]) return TrivialCase(r) end
  fname=r.name*".tmp"
  for f in [:curve,:discy,:ismonic,:roots,:loops,:segments,:points,:zeros]
    cut(string(r.name,".$f=",r[f],"\n\n");before=>"+*",file=fname)
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
    read(Strin(r.name,".",i))
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
  deriv = derivative(p)
  for cnt in 1:lim
    a=p(z)
    b=deriv(z)
    c=iszero(b) ? a : a/b
    err=abs(BigNorm(c))
    if iszero(err) err=(precision/100)/(degree(p)+1) end
    if err>precision err=precision end
#   pr=maximum(0,-1-DecimalLog(err))
#   z=ComplexRational(evalf(z-c,pr+1))
#   if showall println(cnt,": ",evalf(z, pr)) end
    z-=c
    if showall println(cnt,": ",z) end
#   if 10^-pr*(length(p)-1)<=precision
    if err*degree(p)<=precision
      if show print(cnt) end
#     return [z,10^-pr]
      return (z,err)
    end
  end
  if show
    println("\n****** Non-Convergent Newton after ", lim," iterations ******")
    xprintln("p=",p," initial=",z," prec=",precision)
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
  v=map(e->NewtonRoot(p,e,radv),v)
  if !any(isnothing,v) && Dispersal(first.(v))[1]/safety/2>=maximum(last.(v))
    return first.(v)
  end
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
function convert(ll)
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
  (points=points, segments, loops)
end

"""
'LoopsAroundPunctures(points)'

`points`  should be complex numbers. The function computes piecewise-linear
loops representing generators of the fundamental group of `ℂ -{points}`.

```julia-repl
julia> LoopsAroundPunctures([0])
(points = Complex{Int64}[0 - 1im, -1 + 0im, 1 + 0im, 0 + 1im], segments = [[1, 2], [1, 3], [2, 4], [3, 4]], loops = [[4, -3, -1, 2]])
```

The output is a named tuple with fields 
  - `points`: a list of complex  numbers. 
  - `segments`:  a list of oriented segments, each of them  encoded by the
    list of the positions in 'points' of  its two endpoints. 
  - `loops`: a list of loops. Each loops is a list  of integers representing
    a  piecewise  linear  loop,  obtained  by  concatenating the `segments`
    indexed  by the  integers, where  a negative  integer is  used when the
    opposed orientation of the segment is taken.
"""
function LoopsAroundPunctures(originalroots)
  roots=originalroots
  n=length(roots)
  average=sum(roots)/n
  sort!(roots, by=x->abs2(x-average))
  if n==1 return convert([roots[1].+[1,im,-1,-im,1]]) end
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
  rs=map(y->real(y[:y]), ys)
  is=map(y->imag(y[:y]), ys)
  minr=real(ys[argmin(rs)][:y])
  maxr=real(ys[argmax(rs)][:y])
  mini=imag(ys[argmin(is)][:y])
  maxi=imag(ys[argmax(is)][:y])
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
  @show ys[1]
  convert(map(y->y[:loop], ys))
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
Gapjm.gap(f::Float64)="evalf(\"$f\")"
Gapjm.gap(f::Complex{Float64})="Complex("*gap(real(f))*","*gap(imag(f))*")"
# for each point of a find closest point in b
# Complain if the result is not a bijection between a and b of if
# the distance between an a and the correponding b is bigger than 1/10
# of minimum distance between two b's
function fit(a, b)
# a=map(ComplexRational, a)
# b=map(ComplexRational, b)
  dm=map(p->minimum(BigNorm.(b.-p)),a)
  monodromyError=maximum(dm)
  println("# Monodromy error==",monodromyError)
  if maximum(dm)>Dispersal(b)[1]/10 error("monodromy error too big") end
  dm=map(i->findfirst(z->BigNorm(z-a[i])==dm[i],b),1:length(dm))
  if sort(dm)!=1:length(dm) error("monodromy cannot find perm") end
  return b[dm]
end

# Decimal Log of Norm of polynomial d evaluated at point p
function normdisc(d, p)
  p=SmallNorm(prod(map(f->f(p), d)))
  if log10(p)==0 return round(-log10(1/p);digits=3)
  else return round(log10(p);digits=3)
  end
end

# keep only 3 significant digits of x
approx=x->round(x;digits=3)

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
    FundamentalGroup((x+3*y)*(x+y-1)*(x-y),2);|

  ....

|    5.3.6. ***rejected
    4.3.6.<15/16>mindist=3 step=1/2 total=0 logdisc=1 ***rejected
    3.3.4.<15/16>mindist=3 step=1/4 total=0 logdisc=1 ***rejected
    3.3.4.<15/16>mindist=3 step=1/8 total=0 logdisc=1 ***rejected
    3.3.3.<15/16>mindist=3 step=1/16 total=0 logdisc=1
    3.2.3.<15/16>mindist=2.92 step=1/16 total=1/16 logdisc=1
    3.3.3.<15/16>mindist=2.83 step=1/16 total=1/8 logdisc=1
    3.2.3.<15/16>mindist=2.75 step=1/16 total=3/16 logdisc=1
    3.3.3.<15/16>mindist=2.67 step=1/16 total=1/4 logdisc=1
    ======================================
    =    Nontrivial braiding = 2         =
    ======================================
    3.2.3.<15/16>mindist=2.63 step=1/16 total=5/16 logdisc=1
    3.2.3.<15/16>mindist=2.75 step=1/16 total=3/8 logdisc=1
    3.3.3.<15/16>mindist=2.88 step=1/16 total=7/16 logdisc=1
    3.2.3.<15/16>mindist=3 step=1/16 total=1/2 logdisc=1
    3.3.3.<15/16>mindist=3.13 step=1/16 total=9/16 logdisc=1
    3.2.3.<15/16>mindist=3.25 step=1/16 total=5/8 logdisc=1
    3.3.3.<15/16>mindist=3.38 step=1/16 total=11/16 logdisc=1
    3.2.3.<15/16>mindist=3.5 step=1/16 total=3/4 logdisc=1
    3.2.3.<15/16>mindist=3.63 step=1/16 total=13/16 logdisc=1
    3.2.3.<15/16>mindist=3.75 step=1/16 total=7/8 logdisc=1
    3.2.3.<15/16>mindist=3.88 step=1/16 total=15/16 logdisc=1 ***up
    # Monodromy error=0
    # Minimal distance=2.625
    # Minimal step=1/16=-0.05208125+0.01041875I
    # Adaptivity=10
    monodromy[15]:=B(2);
    # segment 15/16 Time=0.2sec|

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
  step=1
  minstep=step
  total=0
  nextzeros=nothing
  while true
    next=prev+step*v
    P=Pol(complexmvp(r.curve)(y=next))
    nextzeros=SeparateRootsInitialGuess(P, prevzeros, 100)
    if isnothing(nextzeros) || 
       (1+maximum(BigNorm.(nextzeros-prevzeros))≈1) && step>1//16)
      rejected=true
    else
      dm=map(i->minimum(BigNorm.(prevzeros[i].-prevzeros[filter(j->j!=i,1:n)])),1:n)
      mdm=minimum(dm)
      if step<1 ipr("<$segno/",length(r.segments),">mindist==",approx(mdm),
         " step==$step total==$total logdisc==",normdisc(r.discyFactored,next))
      end
      dn=map(i->BigNorm(prevzeros[i]-nextzeros[i]),1:n)
      rejected=any(i->dm[i]<VKCURVE[:AdaptivityFactor]*dn[i],1:n)
      if !rejected && mdm<mindm mindm=mdm end
    end
    if rejected
      step/=2
      ipr(" ***rejected\n")
      if step<minstep minstep=step end
    else
      total+=step
      if all(i->dm[i]>2*VKCURVE[:AdaptivityFactor]*dn[i],1:n) && total+step!=1
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
  res*=LBraidToWord(prevzeros,fit(nextzeros,r.zeros[q]),r.B)
  pr("# Minimal distance==", approx(mindm), "\n")
  pr("# Minimal step==", minstep, "==", approx(v*minstep), "\n")
  pr("# Adaptivity==", VKCURVE[:AdaptivityFactor], "\n")
  res
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
  for k in 1:n
    for l in k+1:n
      rv=(real(v1[k])-real(v1[l]))
      iv=(imag(v1[k])-imag(v1[l]))
      if abs(iv*rv)>10^-8 tan=min(tan,abs(rv/iv)) end
      rv=(real(v2[k])-real(v2[l]))
      iv=(imag(v2[k])-imag(v2[l]))
      if abs(iv*rv)>10^-8 tan=min(tan,abs(rv/iv)) end
    end
  end
  [v1, v2]*Complex(1,-tan/2)
end

"""
'LBraidToWord(v1,v2,B)'

This function converts  the linear braid given by `v1`  and `v2` into an
element of the braid group `B`.

|    gap> B:=Braid(CoxeterGroupSymmetricGroup(3)); 
    function ( arg ) ... end
    gap> i:=Complex(0,1);
    I
    gap> LBraidToWord([1+i,2+i,3+i],[2+i,1+2*i,4-6*i],B);
    1|

The  list `v1` and `v2` must have the same length, say `n`. The braid group
`B`   should  be  the   braid  group  on   `n`  strings,  in  its  CHEVIE
implementation.  The elements of  `v1` (resp. `v2`)  should be `n` distinct
complex  rational  numbers.  We  use  the  Brieskorn  basepoint, namely the
contractible  set  `C+iV_ℝ`  where  `C`  is  a real chamber; therefore the
endpoints need not be equal (hence, if the path is indeed a loop, the final
endpoint must be given). The linear braid considered is the one with affine
strings  connecting each point in `v1`  to the corresponding point in `v2`.
These strings should be non-crossing. When the numbers in `v1` (resp. `v2`)
have  distinct real parts, the  real picture of the  braid defines a unique
element  of `B`. When some real parts are equal, we apply a lexicographical
desingularization,  corresponding  to  a  rotation  of  `v1` and `v2` by an
arbitrary small positive angle.
"""
# Computes, from a piecewise linear braid, the corresponding element of B_n.
# Convention:
#    we  use the  Brieskorn basepoint,  namely the  contractible set C+iV_R
#    where  C is a real chamber; therefore  the endpoints need not be equal
#    (hence,  if the  path is  indeed a  loop, the  final endpoint  must be
#    given).
#    To get rid of singular projections, a lexicographical
#    desingularization is applied.
##########################################
# Input: pair of n-tuples of complex rational numbers, ambient braid function
#         (the ambient braid group should be B_n)
# Output: the corresponding element of B_n
#########################################
# two printlevel control fields: VKCURVE.showSingularProj
#				 VKCURVE.showBraiding
#########################################

# Deals with linear braids
# 1) singular real projections are identified
# 2) calls starbraid for each
function LBraidToWord(v1, v2, B)
  n=length(v1)
  x1=real.(v1)
  y1=imag.(v1)
  x2=real.(v2)
  y2=imag.(v2)
  c=combinations(1:n,2)
  if any(x->isapprox(x1[x[1]],x1[x[2]];rtol=10^-8),c) || 
     any(x->isapprox(x2[x[1]],x2[x[2]];rtol=10^-8),c)
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
  tcrit=sort(unique(crit))
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
    xcrit=sort(unique(xt))
    for x in xcrit
      posx=findfirst(==(x),xt)
      nx=count(==(x),xt)
      res*=starbraid(yt[posx:(posx+nx)-1], posx-1, B)
    end
    u=t
  end
  if VKCURVE[:showBraiding]
   if !isempty(tcrit)
      if VKCURVE[:showInsideSegments]
        println("======================================")
        println("==    Nontrivial braiding == ",lpad(res,10),"==")
        println("======================================")
      else println("==    Nontrivial braiding == ",lpad(res,10),"==")
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

|    gap> B:=Braid(CoxeterGroupSymmetricGroup(3));
    function ( arg ) ... end
    gap> b:=B(1);
    1
    gap> BnActsOnFn(b,FreeGroup(3));
    GroupHomomorphismByImages( Group( f.1, f.2, f.3 ), Group( f.1, f.2, f.3 ), 
    [ f.1, f.2, f.3 ], [ f.2, f.2^-1*f.1*f.2, f.3 ] )
    gap> BnActsOnFn(b^2,FreeGroup(3));
    GroupHomomorphismByImages( Group( f.1, f.2, f.3 ), Group( f.1, f.2, f.3 ), 
    [ f.1, f.2, f.3 ], [ f.2^-1*f.1*f.2, f.2^-1*f.1^-1*f.2*f.1*f.2, f.3 ] )|

The second input is the free group on `n` generators. The first input is
an  element  of  the  braid  group  on  `n`  strings,  in  its  CHEVIE
implementation.
"""
BnActsOnFn(b,F)=GroupHomomorphismByImages(F,F,gens(F),hurwitz(b,gens(F)))

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
    gap> p:=PresentationFpGroup(g);DisplayPresentation(p);
    << presentation with 3 gens and 4 rels of total length 16 >>
    1: c=b
    2: b=c
    3: bab=aba
    4: aba=bab
    gap> SimplifyPresentation(p);DisplayPresentation(p);
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
  zero=r.zeros[r.basepoint]
  dist=abs2.(zero-r.height)
  height=zero[argmin(dist)]
  basestring=count(z->(real(z),imag(z))<(real(height),imag(height)),zero)
  @show basestring,r.braids
  fbase=gens(F)[basestring]
  rels=[]
  auts=map(b->BnActsOnFn(b, F),r.braids)
  for aut in auts
# Find an element conjugator such that aut(fbase)^inv(conjugator)=fbase
    ifbase=Image(aut, fbase)
    conjugator=F.gens[1]/F.gens[1]
    choices=vcat(map(f->[f,f^-1], gens(f)[1:n])...)
    while true
      k=0
      while true
        k+=1
        if LengthWord(choices[k]*ifbase) < LengthWord(ifbase) break end
      end
      ifbase=(choices[k]*ifbase) // choices[k]
      conjugator=choices[k]*conjugator
      if LengthWord(ifbase)==1 break end
    end
# Replacing aut by  correctaut:= Conj(conjugator)*aut
    conj=GroupHomomorphismByImages(F, F, gens(F), gens(F).^inv(conjugator))
    conj[:isMapping]=true
    correctaut=CompositionMapping(conj, aut)
    if Position(auts, aut) > length(r.verticallines)
      rels=Append(rels, map(f->Image(correctaut, f)/ f, gens(F)[1:n]))
    else
      g=gens(F)[Position(auts, aut)+n]
      append!(rels, map(f->(Image(correctaut, f)*g)/ (g*f), gens(F)[1:n]))
    end
  end
  push!(rels, fbase)
  F/rels
end
