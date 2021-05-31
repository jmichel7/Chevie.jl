##############################################################################
##
#A  segtobrd.g       VKCURVE package         Jean Michel
##
#Y  Copyright (C) 2001 - 2002  University Paris VII, France.
##
##  This file holds the implementation of ApproxFollowMonodromy
## 
#############################################################################
# arguments: <Global record>,<segment No>,<print function>
ApproxFollowMonodromy:=function(r,segno,pr)local p,q,v,step,prev,next,res,
  prevzeros,nextzeros,dm,total,dn,i,mindm,minstep,monodromyError,P,ratio,n,
  rejected,mdm,iszero,fit,normdisc,approx,ipr;

  # for each point of a find closest point in b
  # Complain if the result is not a bijection between a and b of if
  # the distance between an a and the correponding b is bigger than 1/10
  # of minimum distance between two b's
  fit:=function(a,b)local dm,monodromyError;
    a:=List(a,ComplexRational);
    b:=List(b,ComplexRational);
    dm:=List(a,p->Minimum(List(b,z->BigNorm(z-p))));
    monodromyError:=Maximum(dm);
    pr("# Monodromy error=",evalf(monodromyError,10),"\n");
    if Maximum(dm)>Dispersal(b)[1]/10 then Error("monodromy error too big");fi;
    dm:=List([1..Length(dm)],i->PositionProperty(b,z->BigNorm(z-a[i])=dm[i]));
    if Set(dm)<>[1..Length(dm)] then Error("monodromy cannot find perm");fi;
    return b{dm};
  end;

  # Decimal Log of Norm of polynomial d evaluated at point p
  normdisc:=function(d,p)
    p:=SmallNorm(Product(List(d,f->ScalMvp(Value(f,["x",p])))));
  # p:=Rational(p);
    if DecimalLog(p)=0 then return -DecimalLog(1/p);
    else return DecimalLog(p);
    fi;
  end;

  # keep only 3 significant digits of x
  approx:=x->evalf(x,-DecimalLog(Rational(x))+2);

  if VKCURVE.showInsideSegments then ipr:=Print;else ipr:=Ignore;fi;
  iszero:=x->x=0 or (IsDecimal(x) and x.mantissa<10);
  p:=r.segments[segno][1];
  q:=r.segments[segno][2];
  res:=r.B();
  prevzeros:=r.zeros[p];
  n:=Length(prevzeros);
  if n=1 then return r.B();fi;
  mindm:=Dispersal(prevzeros)[1];
  p:=ComplexRational(r.points[p]);
  v:=ComplexRational(r.points[q])-p;prev:=p;
  step:=1;minstep:=step;
  total:=0;
  repeat 
    next:=prev+step*v;
   P:=ScalMvp(Coefficients(Value(r.curve,["y",ComplexRational(next)]),"x"));
    nextzeros:=SeparateRootsInitialGuess(P,prevzeros,100);
    if false=nextzeros 
      or (iszero(Maximum(List(nextzeros-prevzeros,BigNorm))) and step>1/16)
    then rejected:=true;
    else
      dm:=List([1..n],i->Minimum(List(prevzeros[i]
	-prevzeros{Filtered([1..n],j->j<>i)},BigNorm)));
      mdm:=Minimum(dm);
      if step<1 then ipr("<",segno,"/",Length(r.segments),">",
	      "mindist=",approx(mdm)," step=",step," total=",total,
	      " logdisc=",normdisc(r.discyFactored,next));
      fi;
      dn:=List([1..n],i->BigNorm(prevzeros[i]-nextzeros[i]));
      rejected:=ForAny([1..n],i->dm[i]<VKCURVE.AdaptivityFactor*dn[i]);
      if not rejected and mdm<mindm then mindm:=mdm;fi;
    fi;
    if rejected then 
      step:=step/2;ipr(" ***rejected\n");
      if step<minstep then minstep:=step;fi;
    else
      total:=total+step;
      if ForAll([1..n],i->dm[i]>2*VKCURVE.AdaptivityFactor*dn[i]) 
	and total+step<>1 then 
	step:=2*step;ipr(" ***up");
      fi;
      ipr("\n");
      if total<>1 then 
	res:=res*LBraidToWord(prevzeros,nextzeros,r.B);
	prevzeros:=nextzeros;
      fi;
      prev:=next;
    fi;
    if total+step>1 then step:=1-total;fi;
  until total=1;
  res:=res*LBraidToWord(prevzeros,fit(nextzeros,r.zeros[q]),r.B);
  pr("# Minimal distance=",evalf(evalf(mindm),10),"\n");
  pr("# Minimal step=",minstep,"=",evalf(v*minstep,10),"\n");
  pr("# Adaptivity=",VKCURVE.AdaptivityFactor,"\n");
  return res;
end;
