VKCURVE=Dict(
:name=>"vkcurve",
:version=>"1.2",
:date=>[2009,3],
:homepage=>"http://webusers.imj-prg.fr/~jean.michel/vkcurve.html",
:copyright=>
"(C) David Bessis, Jean Michel -- compute Pi_1 of hypersurface complements",
:monodromyApprox=>false,
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
:showallnewton=>false,# for NewtonRoot
:NewtonLim=>800,      # for NewtonRoot
:AdaptivityFactor=>10, # for ApproxFollowMonodromy
:shrinkBraid=>false,
:mvp=>1)

function SetPrintLevel(printlevel)
  VKCURVE[:showSingularProj]=printlevel>=2
  VKCURVE[:showBraiding]=printlevel>=2
  VKCURVE[:showLoops]=printlevel>=2
  VKCURVE[:showAction]=printlevel>=2
  VKCURVE[:showSegments]=printlevel>=1
  VKCURVE[:showInsideSegments]=printlevel>=2
  VKCURVE[:showWorst]=printlevel>=2
  VKCURVE[:showZeros]=printlevel>=2
  VKCURVE[:showNewton]=printlevel>=2
  VKCURVE[:showgetbraid]=printlevel>=1
  VKCURVE[:showRoots]=printlevel>=2
end

@GapObj struct VK{T}
  curve::Mvp{T,Int}
  ismonic::Bool
end

function Base.show(io::IO,r::VK)
  if haskey(r,:presentation) DisplayPresentationr.presentation
  else print(io,r.prop)
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
  Inherit(r, LoopsAroundPunctures(r.roots))
# here we  have loops around the  'true' roots and around  the 'extra'
# roots Difference(r.roots,r.trueroots). We get rid of the extra loops
# and the associated segments and points, first saving the basepoint.
# (its location is known now, and maybe not later? J.M.) 
  r.basepoint=r.loops[1][1]<0 ? r.segments[-r.loops[1][1]][2] :
                                r.segments[r.loops[1][1]][1]
  r.loops=r.loops[map(z->Position(r.roots,z),r.ismonic ? r.roots : r.trueroots)]
  segmentNumbers=union(map(x->abs.(x),r.loops))
  r.loops=map(x->map(y->y<0 ? -Position(segmentNumbers, -y) :
                               Position(segmentNumbers, y), x),r.loops)
  r.segments = r.segments[segmentNumbers]
  uniquePoints = union(r.segments)
  r.segments = map(x->map(y->Position(uniquePoints, y), x), r.segments)
  r.points = r.points[uniquePoints]
  r.basepoint = Position(uniquePoints, r.basepoint)
  if VKCURVE[:showSegments]
    println("# There are ",length(r.segments)," segments in ",length(r.loops)," loops")
  end
  if VKCURVE[:showWorst]
   l=map(enumerate(r.segments))do i,s
      m=findmin(map(z->DistSeg(z, r.points[s[1]], r.points[s[2]]),r.roots))
      [m[1], i, m[2]]
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
    r.dispersal=1//1000
  end
# and round points to m/100
  m = -DecimalLog(Rational(r.dispersal)//100)
  r.points = map(y->evalf(y, m),r.points)
end

# VKCURVE.Zeros(r)
#   r should be a record with fields
#   .points
#   .curve
#
# It computes r.zeros as the zeros of the .curve at each of the .points
function Zeros(r)
  if VKCURVE[:showRoots]
    println("Computing zeros of curve at the ", length(r.points), " segment extremities...")
  end
  r.zeros = []
  min = []
  for i = 1:length(r.points)
    if VKCURVE[:showZeros]
        print("<", i, "/", length(r.points), ">")
    end
    r.zeros[i]=SeparateRoots(r.curve(y=r.points[i]), 1000)
    if length(r.zeros[i]) > 1
        m = Dispersal(r.zeros[i])
        m[1] = evalf(m[1], -(DecimalLog(m[1] // 1000)))
        if VKCURVE[:showZeros] println(" d==", m) end
        min[i] = [m[1], i]
    end
  end
  if VKCURVE[:showWorst] && length(r.zeros[1]) > 1
    sort!(min)
    println("worst points:")
    for i in 1:min(5,length(min))
        println(min[i][2],": mindist(zeros)==",min[i][1])
    end
  end
end

# PrepareCurve(curve)
#   curve should be an Mvp in x and y.
#   This  function  makes  sure  the  curve  is  quadratfrei  and makes its
#   coefficients  in  CF(4)=GaussianRationals (the coeffs could
#   be  decimals or complex  decimals) returns a  record with fields .curve
#   and .ismonic if the curve is monic in x
function PrepareCurve(curve)
  d=gcd(curve, derivative(curve,:x))
  if degree(d,:x)>0
    xprintln("**** Warning: curve is not quadratfrei: dividing by ", d)
    curve=exactdiv(curve,d)
  end
  if eltype(curve)<:Complex
    curve=Mvp(ModuleElt(m=>Cyc(c) for (m,c) in curve.d))
  end
  VK(curve,degree(Pol(curve,:x)[end])==0)
end

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
  r[:curveVerticalPart] = gcd(values(coefficients(r.curve,:x))...)
  if VKCURVE[:showRoots] && degree(r[:curveVerticalPart]) > 0
    println("Curve has ",degree(r[:curveVerticalPart])," linear factors in y")
  end
  r[:nonVerticalPart] = exactdiv(r.curve,r[:curveVerticalPart])
  d=discy(r[:nonVerticalPart])
  if iszero(d)
    error("Discriminant is 0 but ", r.curve, " should be quadratfrei")
  end
  if VKCURVE[:showRoots] print("Discriminant has ", degree(d), " roots, ") end
  d = exactdiv(d,gcd(d, derivative(d)))
  if VKCURVE[:showRoots] println(" of which ", degree(d), " are distinct") end
  common = gcd(d(Mvp(:y)), r[:curveVerticalPart])
  if VKCURVE[:showRoots] && degree(common) > 0
      println(" and of which ", degree(common), " are roots of linear factors")
  end
  d=exactdiv(d,Pol(common))
  d/=d[end]
  r.discy = d
  r[:discyFactored] = [d]
  # r[:discyFactored] = Factors(r.discy)
end

function GetDiscyRoots(r)
  if VKCURVE[:showRoots] println("Computing roots of discriminant...") end
  if degree(r[:curveVerticalPart])==0 r.roots=[]
  else r.roots = SeparateRoots(r[:curveVerticalPart], 1000)
  end
  append!(r.roots, vcat(map(p->SeparateRoots(p, 1000), r[:discyFactored])))
end

function Braids(r)
  pr=VKCURVE[:showgetbraid] ? println : Ignore
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

function Segment(r, segno)
  if haskey(r, :name)
    PrintTo(SPrint(r.name, ".", segno), "")
    pr=function(arg...,)
      AppendTo(string(r.name,".",segno),arg...)
      println(arg...)
    end
  elseif VKCURVE[:showSegments] pr = print
  else pr = Ignore
  end
  tm=Runtime()
  r.monodromy[segno]=(VKCURVE[:monodromyApprox] ?
                  ApproxFollowMonodromy : FollowMonodromy)(r, segno, pr)
  tm=Runtime()-tm
  if haskey(r, :name) pr(r.name, ".") end
  pr("monodromy[", segno,"]=",word(r.monodromy[segno]), "\n\n")
  pr("# segment ",segno,"/",length(r.segments)," Time==",evalf(tm//1000,1),
      "sec\n")
end

function TrivialCase(r)
  r.presentation=PresentationFpGroup(FreeGroup(degree(r.curve,:x)))
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
  r[:B]=BraidMonoid(CoxSym(length(r.zeros[1])))
  if !haskey(r,:monodromy) r.monodromy=[] end
  for segno in range VKCURVE[:Segment](r, segno) end
end

function Finish(r)
  VKCURVE[:Braids](r)
  if r.ismonic r[:rawPresentation] = VKQuotient(r.braids)
  else r[:rawPresentation] = DBVKQuotient(r)
  end
  r.presentation = PresentationFpGroup(r[:rawPresentation])
  ShrinkPresentation(r.presentation)
  r[:rawPresentation] = PresentationFpGroup(r[:rawPresentation])
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
  section=section//gcd(section, r[:curveVerticalPart])
  println("Curve is not monic in x -- Trivializing along horizontal line x == ", height)
  r1=Dict{Symbol,Any}(:height=>height,:input=>r.curve,:curve=>r.curve*
                      (Mvp(:x)-height))
  if haskey(r,:name) r1[:name] = r.name end
  Discy(r1)
# set trueroots  to roots  of Discy  which do  not only  correspond to
# intersections of the curve with the chosen horizontal
  r1[:trueroots] = r.roots
  r1[:verticallines] = r.roots[1:degree(r[:curveVerticalPart])]
  r1[:roots] = copy(r1[:trueroots])
  r1[:ismonic] = false
  append!(r1[:roots],SeparateRoots(section, 1000))
  return r1
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
  r[:B] = BraidMonoid(CoxSym(length(r.zeros[1])))
  r.monodromy = []
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
  r = VKCURVE[:PrepareCurve](curve)
  VKCURVE[:Discy](r)
  r.name = name
  VKCURVE[:GetDiscyRoots](r)
  if isempty(r.roots) return VKCURVE[:TrivialCase](r) end
  if !r.ismonic r=VKCURVE[:SearchHorizontal](r) end
  VKCURVE[:Loops](r)
  VKCURVE[:Zeros](r)
  if isempty(r.zeros[1]) return VKCURVE[:TrivialCase](r) end
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
  r[:B] = BraidMonoid(CoxSym(length(r.zeros[1])))
  for i in 1:length(r.segments)
    read(Strin(r.name,".",i))
    if isnothing(r.monodromy[i]) println("***** ", i, " missing ****") end
    r.monodromy[i] = r[:B](r.monodromy[i]...)
  end
  return VKCURVE[:Finish](r)
end
