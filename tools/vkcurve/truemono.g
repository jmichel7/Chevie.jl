##############################################################################
##
#A  truemono.g       VKCURVE package         David bessis
##
#Y  Copyright (C) 2001 - 2002  University Paris VII, France.
##
##  This file holds the implementation of FollowMonodromy
## 
#############################################################################
# Exact computation of the monodromy braid along a segment
# arg[1]: global VKCURVE record
# arg[2]: segment number
# arg[3]: Print function (to screen, to file, or none)
# inner print level controlled by VKCURVE.showinsidesegments

FollowMonodromy:=function(arg) local globalrec,a,b,
      p,px,dpdx,py,d,y0,y1,B,v,s,t,time,increment,res,qr,
      z,i,R,dist,pos,allowed,newv,
      protected,protp,protdpdx,prec,
      pt,dpdxt,ptz,dpdxtz,cptz,cdpdxtz,ptzr,ptzi,dpdxtzr,dpdxtzi,
      debut,fin,py,dpy,k,l,phi,RR,steps,x,
      square,file,iPrint,sPrint,
      mynorm,lowevalf,binlowevalf,myfit,mybinlog,mynewton,Sturm,
      adapt,stu;

  if VKCURVE.showInsideSegments then iPrint:=Print;
			        else iPrint:= function(arg); end;
  fi;

  mynorm:=x->x.r^2+x.i^2;

  # computes the lower approximation of the rational a by a
  # decimal with denominator 10^k
  lowevalf:=function(a,k)local b;
    b:=EuclideanQuotient(Numerator(a)*10^k,Denominator(a));
    if a >= 0 then return  b/10^k;
    else return (b-1)/10^k;
    fi;
  end;

  # computes the lower approximation of the rational a by a
  # rational with denominator 2^k
  binlowevalf:=function(a,k)local b;
    b:=EuclideanQuotient(Numerator(a)*2^k,Denominator(a));
    if a >= 0 then return b/2^k;
    else return (b-1)/2^k;
    fi;
  end;
    
  # computes an integral approximation of the binary logarithm
  # of the positive rational p
  mybinlog:=function(p) local q,k; k:=0; q:=p;
    while q < 1 do q:=2*q; k:=k+1; od;
    return k;
  end;

  # truncated iteration of the Newton method
  mynewton:=function(p,deriv,z)local a,b,c,d,err,prec,z;
    d:=Length(p)-1; a:=ValuePol(p,z); b:=ValuePol(deriv,z);
    if b=0*b then c:=a; Print("NewtonError\n");else c:=a/b;fi;
    err:=d*BigNorm(c);
    if err = 0 then prec:=1; else prec:=Maximum(0,-DecimalLog(err))+2; fi;
    return ComplexRational(evalf(z-c,prec));
  end;

  # for each point of a find closest point in b
  myfit:=function(a,b) local R,dist,res,k,l,d,z,i;
    d:=Length(a);
    dist:=List([1..d],i->[1..d]);
    for k in [1..d] do 
       for l in [k+1..d] do
	 dist[k][l]:=mynorm(v[k]-v[l]);
	 dist[l][k]:=dist[k][l];
	 dist[k][k]:=dist[k][l];
       od;
    od;
    dist[d][d]:=dist[d][d-1];
    R:=List([1..d],k->Minimum(dist[k])/4);
    res:=[1..d];
    for k in [1..d] do
	z:=Filtered(b,i->mynorm(i-a[k])<R[k]);
	if Length(z) <> 1 then Error("something's wrong"); fi;
	res[k]:=z[1];
    od;
  return res;
  end;

  # computes square of a one variable polynomial, encoded in a list l
  square:=function(l) local k,res;
    res:=l*l[1];
    for k in [2..Length(l)] do
      Add(res,0);
      res:=res + Concatenation([1..k-1]*0,l[k]*l); od;
    return res;
  end;

  # pp is polynomial
  # if pp is positive  at time
  #   the function returns some rational number t such that
  #    time < t <= 1  and  pp  is positive on [time,t]
  # otherwise returns 0
  # [third input and second output is an adaptive factor to 
  #  accelerate the computation]
  Sturm:=function(pp,time,adapt) local pol,q,i,j,l,t,m,k;
    q:=X(Rationals);
    pol:=ValuePol(pp,(1-q)*time+q);
    pol:=pol.coefficients;
    if pol[1]<=0 then Print("*****"); return [0,0]; fi;
    l:=Length(pol);
    k:=2;
    #Print("\n");
    #for i in [1..l] do if pol[i]>0 then Print("+"); elif pol[i]=0 then
    #  Print("0"); else Print("-"); fi; od; Print("\n");
    while (k<l) and (pol[k]>=0) do k:=k+1; od;
    while (k<l) and (pol[k]<=0) do k:=k+1; od;
    for i in [k..l] do
	if pol[i]>0 then pol[i]:=0; fi;
    od;
    #for i in [1..l] do if pol[i]>0 then Print("+"); elif pol[i]=0 then
    #  Print("0"); else Print("-"); fi; od; Print("\n");
    t:=1/2^adapt;
    m:=adapt;
    while ValuePol(pol,t) <= 0
    do t:=t/2; m:=m+1;
    od;
    iPrint(m);
    if (m=adapt) and (adapt>0) then
	   if ValuePol(pol,3*t/2) > 0
	   then 
		      if ValuePol(pol,2*t)>0
		      then return [(1-2*t)*time+2*t,adapt-1]; 
		      else return [(1-3*t/2)*time+3*t/2,adapt-1]; fi;
	   else return [(1-t)*time+t,adapt];
	   fi;
    else return [(1-t)*time+t,m];
    fi;
  end;

  globalrec:=arg[1];
  x:=Mvp("x");
  #We encode the polynomial as the vector p of its coefficients in y,
  # (initial conversion tricks used only to compute dp/dx)
  if IsBound(globalrec.monop) 
  then p:=globalrec.monop; dpdx:=globalrec.monodpdx;
  else iPrint("Initializing monodromy data\n");
    p:=ShallowCopy(globalrec.curve); p:=Coefficients(p,"y");
    p:=List(p,i->ScalMvp(Coefficients(i,"x")));
    dpdx:=List(p,Derivative); p:=List(p,i->ValuePol(i,x));
    p:=p*Mvp(1); dpdx:=List(dpdx,i->ValuePol(i,x));
    dpdx:=dpdx*Mvp(1); globalrec.monop:=p;
    globalrec.monodpdx:=dpdx;
  fi;
  a:=globalrec.segments[arg[2]][1];
  b:=globalrec.segments[arg[2]][2];
  B:=globalrec.B;
  sPrint:=arg[3];
  debut:=ComplexRational(globalrec.points[a]);
  fin:=ComplexRational(globalrec.points[b]);
  v:=List(globalrec.zeros[a],ComplexRational);
  res:=B();
  # If there is only one string, the braid is trivial
  if Length(v)=1 then return res; fi;
  d:=Length(globalrec.zeros[1]);
  t:=Mvp("t"); 
  time:=0; y0:=debut; y1:=fin;
  pt:=ValuePol(p,y1*t+y0*(1-t));
  dpdxt:=ValuePol(dpdx,y1*t+y0*(1-t));
  pt:=Coefficients(pt,"x");
  dpdxt:=Coefficients(dpdxt,"x");
  RR:=List([1..d],k->0);
  adapt:=[1..d]*0;
  protected:=[1..d]*0;
  protp:=[1..d];
  protdpdx:=[1..d];
  steps:=0;
  dist:=List([1..d],k->List([1..d],l->0));
  repeat
    steps:=steps+1;
    iPrint("<",arg[2],"/",Length(globalrec.segments),">");
    iPrint(String(steps,5));
    iPrint(" time=",String(evalf(time),11),"   ");
    for k in [1..d] do 
       for l in [k+1..d] do
	 dist[k][l]:=mynorm(v[k]-v[l]);
	 dist[l][k]:=dist[k][l];
	 dist[k][k]:=dist[k][l];
       od;
    od;
    dist[d][d]:=dist[d][d-1];
    R:=List([1..d],k->Minimum(dist[k])/4);
    for k in [1..d] do
       z:=v[k];
       if (protected[k] > time) and (R[k] >= RR[k]) then
	      iPrint(". ");
       elif (protected[k] > time) then
	     if adapt[k] + 2 < Maximum(adapt) then R[k]:=R[k]/2; fi;
	     iPrint("R");
	     phi:=R[k]*protdpdx[k]-protp[k];
	     stu:=Sturm(phi,time,adapt[k]);
	     s:=stu[1]; adapt[k]:=stu[2];
	   if s>time then protected[k]:=binlowevalf(s,mybinlog(s-time)+3);
		else iPrint("How bizarre...");
	   fi;
	      RR[k]:=R[k];
       else    
	     # if adapt[k] + 2 < Maximum(adapt) then R[k]:=R[k]/2; fi;
	      iPrint("?");
	      ptz:=ValuePol(pt,z);
	      if IsComplex(ptz) then cptz:=[ptz];
	      else cptz:=ScalMvp(Coefficients(ptz,"t"));fi;
	      if Length(cptz)=0 then cptz:=[Complex(0)];fi;
	      cptz:=cptz+Complex(0);
	      dpdxtz:=ValuePol(dpdxt,z);
	      if IsComplex(dpdxtz) then cdpdxtz:=[dpdxtz];
	      else cdpdxtz:=ScalMvp(Coefficients(dpdxtz,"t"));fi;
	      if Length(cdpdxtz)=0 then cdpdxtz:=[Complex(0)];fi;
	      cdpdxtz:=cdpdxtz+Complex(0);
	      ptzr:=List(cptz,i->i.r);
	      ptzi:=List(cptz,i->i.i);
	      dpdxtzr:=List(cdpdxtz,i->i.r);
	      dpdxtzi:=List(cdpdxtz,i->i.i);
	      protp[k]:=d^2*(square(ptzr)+square(ptzi));
	      protdpdx[k]:=square(dpdxtzr)+square(dpdxtzi);
	      if Length(protdpdx[k]) < Length(protp[k])
		   then
		      for i in [1..Length(protp[k])-Length(protdpdx[k])]
		      do Add(protdpdx[k],0); od;
	      elif Length(protdpdx[k]) > Length(protp[k])
		   then 
		      for i in [1..Length(protdpdx[k])-Length(protp[k])]
		      do Add(protp[k],0); od;
	      fi;
	      phi:=R[k]*protdpdx[k]-protp[k];
	      stu:=Sturm(phi,time,adapt[k]);
	      s:=stu[1]; adapt[k]:=stu[2];
		   if s > time then
		     protected[k]:=
		     binlowevalf(s,mybinlog(s-time)+3);
		     else Print("Something's wrong...");
		   fi;
	      RR[k]:=R[k];
	fi;
    od;
    allowed:=Minimum(protected);
    time:=allowed;
    y0:=debut*(1-time)+fin*time;
    py:=ScalMvp(Coefficients(ValuePol(p,y0),"x"));
    dpy:=Derivative(py);
    newv:=[1..d];
    iPrint("\n");
    for k in [1..d] do
	if protected[k] > allowed
		  then newv[k]:=v[k];
		      # iPrint(".");
		  else newv[k]:=mynewton(py,dpy,v[k]);
		    #   Print(" ",(mynorm((newv[k]-v[k])/
		    #      (mynewton(py,dpy,newv[k])-newv[k])),0)," ");
		      # iPrint("n");
	fi;
    od;
    res:=res*LBraidToWord(v,newv,B);
    v:=newv;
  until time=1;
  res:=res*LBraidToWord(v,myfit(v,globalrec.zeros[b]),B);
  sPrint("# The following braid was computed by FollowMonodromy in ",
  steps," steps.\n");
  return res;
end;
