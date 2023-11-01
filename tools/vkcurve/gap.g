p23:=395507812500*x^3-100195312500*x^2*y+1191796875*x^2+131835937500*x*y^3
+1371093750*x*y^2-74250000*x*y+438750*x+175781250000*y^5-28369140625*y^4
+1162187500*y^3-22233750*y^2+213700*y-829;

p24:=-1024*x^3+5632*x^2*y+18*x*y^4-4352*x*y^2-229376*x
-(27/3136)*y^7-67*y^5-5504*y^3-114688*y;

p27:=11664*x^3+5832*x^2*y^2+7776*x^2*y+648*x^2-1944*x*y^5-5508*x*y^4-198*x*y^3
+954*x*y^2+198*x*y+18*x-1404*y^7-3078*y^6-3271*y^5-1094*y^4-204*y^3
-20*y^2-1*y;

p29:=-4*x^8*y^2-512*x^8*y+4*x^7*y^3+516*x^7*y^2-7680*x^7*y-8192*x^7-1*x^6*y^4
-122*x^6*y^3+12543*x^6*y^2+458624*x^6*y+8192*x^6-8*x^5*y^4-6792*x^5*y^3
-453248*x^5*y^2+6629376*x^5*y+7075840*x^5+3*x^4*y^5+1410*x^4*y^4
+112899*x^4*y^3-7151232*x^4*y^2-14635520*x^4*y-7372800*x^4+4*x^3*y^5
+6532*x^3*y^4+2067328*x^3*y^3+6163200*x^3*y^2+6150400*x^3*y+2048000*x^3
-3*x^2*y^6-1286*x^2*y^5-117763*x^2*y^4-231680*x^2*y^3-115200*x^2*y^2
-32000*x*y^5-5216000*x*y^4-15456000*x*y^3-15392000*x*y^2-5120000*x*y
+1*y^7+2*y^6+1*y^5-51200000*y^4-204800000*y^3-307200000*y^2-204800000*y
-51200000;

p31:=1*x^9-3*x^8*y^2+3*x^8-3*x^7*y^2-7773*x^7-2*x^6*y^4+23974*x^6*y^2
-31103*x^6-1290*x^5*y^4+32394*x^5*y^2+15069888*x^5+61566*x^4*y^4
-63016062*x^4*y^2+75551616*x^4+4*x^3*y^6+48198241*x^3*y^4-189034562*x^3*y^2
+151157664*x^3+32397*x^2*y^6+13167363*x^2*y^4-189019008*x^2*y^2
+151165440*x^2+43756197*x*y^6-69961317*x*y^4-63016704*x*y^2+75582720*x
-2*y^8+15746416200*y^6-34991999*y^4-7776*y^2+15116544;

#FundamentalGroup(c,2);

mediatrix:=function(x,y) local mid;
  if x = y then Print("Undefined mediatrix"); return false; fi;
  return (x+y)/2+[E(4),-E(4)]*(x-y);
end;

# value at z of an equation of the line (x,y)
lineq:=function(x,y,z)
  if x.r=y.r then  
    if x.i=y.i then Print("Undefined line\n"); return false;
    else return z.r-x.r; fi;
  else return (y.i-x.i)*(z.r-x.r)/(y.r-x.r)+x.i-z.i;
  fi;
end;

# Computes the intersecting point of two lines, each given by either a pair
# of points or a vector; returns false if the lines are parallel or a pair
# is a single
crossing:=function(arg) local x1,x2,y1,y2,lambdax,mux,lambday,muy,res,resr,resi;
if Length(arg)=4 then x1:=arg[1]; x2:=arg[2]; y1:=arg[3]; y2:=arg[4];
    else  x1:=arg[1][1]; x2:=arg[1][2]; y1:=arg[2][1]; y2:=arg[2][2];
fi;
if (x1=x2) or (y1=y2) then return false; fi; #Undefined line
if x1.r <> x2.r then
       lambdax:=(x1.i-x2.i)/(x1.r-x2.r);
       mux:=-lambdax*x1.r+x1.i;
       if y1.r <> y2.r then
           lambday:=(y1.i-y2.i)/(y1.r-y2.r);
           muy:=-lambday*y1.r+y1.i;
	   if lambdax=lambday then return false; fi;
	   resr:=(muy-mux)/(lambdax-lambday);
	   resi:=lambdax*resr+mux;
	   res:=Complex(resr,resi);
       else res:=crossing(E(3)*x1,E(3)*x2,E(3)*y1,E(3)*y2);
            if res = false then return false; fi;
            res:=res/E(3);
       fi;
else res:=crossing(E(4)*x1,E(4)*x2,E(4)*y1,E(4)*y2);
     if res = false then return false; fi;
     res:=res/E(4);
fi;
return res;
end;


detectsleftcrossing:=function(c,w,y,z) local res,med,a,b,x,k,xx;
res:=[1..Length(c)-1]; med:=mediatrix(y,z); a:=med[1]; b:=med[2];
for k in [1..Length(c)-1] do
    if lineq(a,b,c[k])*lineq(a,b,c[k+1])<= 0
       then
       x:=crossing(a,b,c[k],c[k+1]);
       if x = false then res[k]:=false;
          else xx:=(z-y)/(w[k]-y);
          if xx.i >= 0 then res[k]:=true;
                           else res[k]:=false;
          fi;
       fi;
       else res[k]:=false;
    fi;
od;
return res;
end;

c:=[Complex(-124999999971025551/125000000000000000,
          -13589067290477/2000000000000000000), 
  Complex(-39975786053546411338970619/40000000000000000000000000,
-32875264541918234242255377709133/4853238318027500000000000000000000000) ];
w:=[ Complex(-62499999971025551/62500000000000000,
   -13589067290477/1000000000000000000), 
  Complex(-19975786053546411338970619/20000000000000000000000000,0)];
yy:=-1;
z:=Complex(-62499999971025551/62500000000000000,
   -13589067290477/1000000000000000000);
med:=mediatrix(yy,z); x1:=med[1]; x2:=med[2];
y1:=c[1];y2:=c[2];
Print("x1=",evalf(x1,13),"\n");
Print("x2=",evalf(x2,13),"\n");
Print("y1=",evalf(y1,13),"\n");
Print("y2=",evalf(y2,13),"\n");
Print(crossing(x1,x2,y1,y2),"\n");
