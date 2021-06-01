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
  iszero=x->x+1≈1
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
       (iszero(maximum(BigNorm.(nextzeros-prevzeros))) && step>1//16)
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
