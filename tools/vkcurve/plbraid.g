##############################################################################
##
#A  plbraid.g 2.0      VKCURVE package         David Bessis
##
#Y  Copyright (C) 2001 - 2003  David Bessis.
##
##  This file holds the implementation of LBraidToWord
## 
#############################################################################
# Computes,  from a  piecewise linear  braid, the  corresponding element of
# B_n.
# Convention:
#    we  use the  Brieskorn basepoint,  namely the  contractible set C+iV_R
#    where  C is a real chamber; therefore  the endpoints need not be equal
#    (hence,  if the  path is  indeed a  loop, the  final endpoint  must be
#    given).
#    To get rid of singular projections, a lexicographical
#    desingularization is applied.
####i######################################
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
LBraidToWord:=function(v1,v2,B)
     local des,k,x1,x2,y1,y2,p,q,inv,i,j,tcrit,crit,t,u,ut,
     xut,put,res,xt,yt,x,xcrit,posx,nx,xx,n,starbraid,desingularized;

# Deals with "star" linear braids, those with associated permutation being w_0
  starbraid:=function(y,offset,B)local n,k;
    n:=Length(y);
    if n=1 then return B();fi;
    k:=Position(y,Minimum(y));
    return B([k..n-1]+offset)*
      starbraid(y{Difference([1..n],[k])},offset,B)/B([n+1-k..n-1]+offset);
  end;

  # In case two points have the same real projection, we use
  # a "lexicographical" desingularization by "infinitesimal rotation"
  desingularized:=function(v1,v2) local n,k,l,tan;
    n:=Length(v1);
    tan:=1;
    for k in [1..n] do
      for l in [k+1..n] do
	if (v1[k].i-v1[l].i)*(v1[k].r-v1[l].r) <>0 then
	  tan:=Minimum(tan,AbsInt((v1[k].r-v1[l].r)/(v1[k].i-v1[l].i)));
	fi;
	if (v2[k].i-v2[l].i)*(v2[k].r-v2[l].r) <>0  then
	  tan:=Minimum(tan,AbsInt((v2[k].r-v2[l].r)/(v2[k].i-v2[l].i)));
	fi;
      od;
    od;
    return [v1,v2]*Complex(1,-tan/2);
  end;

  n:=Length(v1);
  x1:=List([1..n],k->v1[k].r);
  y1:=List([1..n],k->v1[k].i);
  x2:=List([1..n],k->v2[k].r);
  y2:=List([1..n],k->v2[k].i);
  if (Length(Set(x1)) < n) or (Length(Set(x2)) < n) then 
  if VKCURVE.showSingularProj
    then   Print("WARNING: singular projection (resolved)\n"); fi;
       des:=desingularized(v1,v2);
       return LBraidToWord(des[1],des[2],B); fi;
  p:=SortingPerm(x1);
  q:=p^-1;
  inv:=[];
  for i in [1..n-1] do
    for j in [i+1..n] do
      if x2[i^q] > x2[j^q] then Add(inv,[i,j]); fi;
    od;
  od;
  crit:=List(inv,k->
      (x1[k[1]^q]-x1[k[2]^q])/(x2[k[2]^q]-x1[k[2]^q]+x1[k[1]^q]-x2[k[1]^q]));
  tcrit:=Set(crit);
  #Print(crit); Print("\n");
  #Print(tcrit); Print("\n");
  res:=B( );
  u:=0;
  for t in tcrit do
     xt:=List([1..n], k->x1[k]+t*(x2[k]-x1[k]));
     yt:=List([1..n], k->y1[k]+t*(y2[k]-y1[k]));
     ut:=(u+t)/2;
     xut:=List([1..n], k->x1[k]+ut*(x2[k]-x1[k]));
     put:=SortingPerm(xut);
     xt:=Permuted(xt,put);
     yt:=Permuted(yt,put);
     xcrit:=Set(xt);
  #Print(xt); Print("\n");
  #Print(xcrit); Print("\n");
     for x in xcrit do
	 posx:=Position(xt,x);
	 nx:=Length(Filtered(xt,xx->x=xx));
	 res:=res*starbraid(yt{[posx..posx+nx-1]},posx-1,B);
	 od;
     u:=t; 
  od;
  if VKCURVE.showBraiding then 
      if tcrit <> [] then
	  if VKCURVE.showInsideSegments then
	     Print("======================================\n");
	     Print("=    Nontrivial braiding = ");
	     Print(String(res,-10));Print("=\n");
	     Print("======================================\n");
	  else
	     Print("=    Nontrivial braiding = ");
	     Print(String(res,-10));Print("=\n");
	  fi;
      fi;
  fi;
  return res;
end;
