# programmes accompagnant le bilan de Luminy 2017
using Chevie
include("getunpdeg.jl")
Chevie.CheckIndexChars=true

KnownSeries(W)=get!(()->getunpdeg(W)[1],W,:series)

"""
`hxy(W,nx,ny;var=false,sh=false)` where `W` is a crg and `nx`, `ny`∈ℚ/ℤ

returns a `UnipotentCharacter`, the `(q,t)-Trace(Dₓ|H^*_c(X(y)))` 
for `y=π^ny` and `x=π^nx` using formula

``∑_{φ∈Irr(H_y)}φ(x)ε(φ)exp(2iπnₓn_yαᵩ)(q^{nₓ}t^{n_y})^(N+N^*-αᵩ) ρᵩ``
where `αᵩ=aᵩ+Aᵩ`

if `var=true` returns a polynomial in variables `u_{y,φ}`;
if `sh=true`  returns a polynomial in variables  `su_{y,φ}`
"""
function hxy(W,nx,ny;var=false,sh=false)
  ey=Root1(;r=ny)
  ex=Root1(;r=nx)
  s=KnownSeries(W)
  s=s[findfirst(x->x.d==ey,s)] # utilise H^*(X(y))=H^*(X(yπ))
  if !(ex in regular_eigenvalues(s.WGL))
    InfoChevie(nx," is not a regular number for ",s.WGL,"\n")
    return false
  end
  q=Mvp(:q);t=Mvp(:t)
  uc=UnipotentCharacters(W)
  ud=degrees(uc,q)[s.charNumbers]
  signes=map(x->sign(Int.(scalar(x(q=ey)))),ud)
  ct=CharTable(s.WGL).irr
  chix=ct[:,position_regular_class(s.WGL,ex)]
  aA=map(P->degree(P)+valuation(P),ud)
  NN=sum(degrees(W))+sum(codegrees(W))
  if var
    if !haskey(W,:hxy) W.hxy=Dict(:eigen=>Dict()) end
    u=rho->Mvp(Symbol("u",modZ(ny),".",rho))
    oy=ct[:,position_regular_class(s.WGL,ey)].//ct[:,1]
    eig=map((o,a)->o*Root1(;r=a*ny^2),oy,aA)
    for i in eachindex(ud)
      W.hxy[:eigen][Symbol("u",modZ(ny),".",i)]=eig[i]
    end
  elseif sh u=rho->Mvp(Symbol("su",ny,".",rho))
  else u=rho->unichar(W,s.charNumbers[rho])
  end
  sum(i->chix[i]*ey^(aA[i]*nx)*
      q^((NN-aA[i])*nx)*t^((NN-aA[i])*ny)*signes[i]*u(i),eachindex(ud))
end

"""
Shintani on a `UnipotentCharacter` with coefficients in `ℂ[q^±1,t^±1]`
(the map `q↦qt`)
"""
function Shintani(u)
  W=u.group
  uc=UnipotentCharacters(W)
  S=fourier(uc)
  O=Diagonal(Cyc.(eigen(uc)))
  Sh=O*inv(S)*(any(haskey.(uc.families,:lusztig)) ? conj(O) : O)
  UniChar(W,permutedims(Sh)*map(x->x(q=Mvp(:q)*Mvp(:t)),u.v))
end

function testΩ(W,x,y)
  u=hxy(W,x,y)
  Ou=sl2action(u,[1 0;1 1])
  Ou-hxy(W,x+y,y)
end

function testS(W,x,y)
  u=hxy(W,x,y)
  Ou=sl2action(u,[0 -1;1 0])
  Ou-hxy(W,y,-x)
end

# action de SL2: S=[0 -1;1 0], Ω=[1 0;1 1]
function sl2action(u::UniChar,m::Matrix)
  W=u.group
  uc=UnipotentCharacters(W)
  O=Diagonal(Cyc.(eigen(uc)))
  S=fourier(uc)
  uv=u.v
  for l in redsl2(m)
    if l==0 uv=S*map(x->x(q=inv(Mvp(:t)),t=Mvp(:q)),uv )
    elseif l>0 uv=O^l*map(x->x(q=Mvp(:q),t=Mvp(:q)^l*Mvp(:t)),uv)
    else uv=conj(O)^-l*map(x->x(q=Mvp(:q),t=Mvp(:q)^-l*Mvp(:t)),uv)
    end
  end
  UniChar(W,uv)
end

function testm(W,x,y,m)
  u=hxy(W,x,y)
  vars=[Mvp(:q)^m[1,1]*Mvp(:t)^m[1,2],Mvp(:q)^m[2,1]*Mvp(:t)^m[2,2]]
  u=sl2action(u,m)
  u=UniChar(W,value.(u.v,:q=>vars[1],:t=>vars[2]))
  u-hxy(W,m'*[x,y]...)
end

testallm(m,i)=map(l->testm(l...,m),testdata[i])

#UnipCharOps.Projection:=function(u,L)
#  return UnipotentCharacter(u.group,List([1..Length(u.v)],
#     function(i) if i in L then return u.v[i]; else return 0; fi; end));
#end;
#
#UnambigFamilies:=function(W)local N,F;
#  N:=Filtered(KnownSeries(W),s->IsBound(s.ambig));
#  N:=Union(List(N,s->Union(List(s.ambig,x->x[2]))));
#  F:=UnipotentCharacters(W).families;
#  F:=Filtered(F,f->Intersection(N,f.charNumbers)=[]);
#  return Union(List(F,f->f.charNumbers));
#end;
#
#AdmissiblePairs:=function(W)local D,res;
#  D:=List(KnownSeries(W),s->s.d);
#  return List(Cartesian(D,D),i->[modZ(i[2]-i[1]),i[1]]);
#end;
#
#g:=function(W,i,j)
#  local Un,a,b;
#  Un:=UnambigFamilies(W);
## Print("unambig:",Un,"\n");
#  a:=hxy(W,i,j);
#  if a=false then return false;fi;
#  b:=hxy(W,i,i+j);
#  if b=false then return false;fi;
## return Projection(b-Shintani(a),Un);
#  return b-Shintani(a);
#end;
#  
## Coefficients of an Mvp in q,t and other variables for each monomial in q,t
#explode2:=function(p)
#  p:=Zip(p.elm,p.coeff,function(e,c)local pt,ct,pq,cq,left;
#    pt:=Position(e.elm,"t");if pt=false then ct:=0;else ct:=e.coeff[pt];fi;
#    pq:=Position(e.elm,"q");if pq=false then cq:=0;else cq:=e.coeff[pq];fi;
#    left:=Filtered([1..Length(e.elm)],i->i<>pt and i<>pq);
#    return [c,rec(elm:=e.elm{left},coeff:=e.coeff{left}),[ct,cq]];
#    end);
#  p:=CollectBy(p,x->x[3]);
#  return List(p,x->Mvp(List(x,x->x[2]),List(x,x->x[1])));
#end;
#
## Cut a unipotent character according to each monomial in q,t
#explode:=function(u)local mm,l;
#  l:=u.v*Mvp("x")^0;
#  mm:=Union(List(l,x->x.elm));
#  return List(mm,m->List(l,function(P)local p;
#    p:=Position(P.elm,m);
#    if p=false then return 0;
#    else return P.coeff[p];fi;end));
#end;
#
## Cut a polynomial in u_{y,φ} according to each eigenvalue of u_{y,φ}
#explodeEigen:=function(W,p)local e,eig;
#  e:=Mvp(p[1])-p[2];
#  eig:=List(Variables(e),v->W.hxy.eigen.(v));
#  eig:=CollectBy(Variables(e),eig);
#  return List(eig,function(c)local s;
#     s:=Filtered([1..Length(e.elm)],i->e.elm[i].elm[1] in c);
#     return Mvp(e.elm{s},e.coeff{s});
#     end);
#end;

# utilise Ω^-1*Sh*Ω^-1=Sh*Ω^-1*Sh
#usebraid:=function(W,p)local og,rm,lm,a,b,vv,v,l;
#  l:=W.hxy.sols;
#  og:=W.hxy.eigen.(p[1]{[2..Length(p[1])]});
#  rm:=Mvp(p[2].elm,Zip(p[2].coeff,p[2].elm,
#    function(c,e)return c*og^-1*W.hxy.eigen.(e.elm[1])^-1;end));
#  a:=List(p[2].elm,e->rec(coeff:=[1],elm:=[ConcatenationString("s",e.elm[1])]));
#  b:=Zip(p[2].coeff,p[2].elm,
#    function(c,e)return c*W.hxy.eigen.(e.elm[1])^-1;end);
#  lm:=Mvp(a,b);
#  rm:=rm-lm;
#  for v in Filtered(Variables(rm),x->x[1]='s') do
#    vv:=First(l,x->x[1]=v);
#    rm:=Value(rm,vv);
#  od;
#  return rm;
#end;
#
#poltovec:=function(W,p)local v,res,i;
#  v:=RecFields(W.hxy.eigen);
#  res:=List(v,x->0);
#  for i in [1..Length(p.elm)] do
#    res[Position(v,p.elm[i].elm[1])]:=p.coeff[i];
#  od;
#  return res;
#end;
#  
#shmatrix:=function(W)local v,m,l;
#  l:=W.hxy.sols;l:=Filtered(l,x->x[1][1]='s');
#  v:=RecFields(W.hxy.eigen);
#  if Length(l)<>Length(v) then Error("not all variables s... found");fi;
#  m:=List(l,e->poltovec(W,e[2]));
#  SortParallel(List(l,x->Position(v,x[1]{[2..Length(x[1])]})),m);
#  return m;
#end;

#omegamatrix:=function(W)
#  return DiagonalMat(List(RecFields(W.hxy.eigen),x->W.hxy.eigen.(x)));
#end;
#
#testsh:=function(W)local l,a,b,p,v,addsol;
#  addsol:=function(v)local p;
#    Add(W.hxy.sols,v);for p in W.hxy.sols do p[2]:=Value(p[2],v);od;
#  end;
#  l:=[];
#  for p in AdmissiblePairs(W) do
#    a:=hxy(W,p[1],p[2],rec(sh:=true));
#    if a<>false then
#      a:=Value(a,["q",Mvp("q")*Mvp("t")]);
#      b:=hxy(W,p[1],p[1]+p[2],rec(var:=true));
#      if b<>false then Add(l,b-a);fi;
#    fi;
#  od;
#  W.hxy.sols:=[];
#  l:=Concatenation(List(l,explode2));
#  while Length(l)>0 do
#    v:=solve(l[1]);addsol(v);l:=Value(l,v);l:=sortbylg(l);
#  od;
#  l:=Filtered(W.hxy.sols,x->x[1] in RecFields(W.hxy.eigen) and 
#    Length(x[2].coeff)>1);
#  l:=Concatenation(List(l,p->explodeEigen(W,p)));
#  while Length(l)>0 do
#    v:=solve(l[1]);addsol(v);l:=Value(l,v);l:=sortbylg(l);
#  od;
#  p:=PositionProperty(W.hxy.sols,x->ForAny(Variables(x[2]),
#    v->v[1]='s'));
#  if p<>false then 
#    Error("one of the s.. variables not resolved",W.hxy.sols[p],"\n");
#  fi;
#  for p in Filtered(W.hxy.sols,x->x[1][1]='s') do
#    v:=usebraid(W,p); if v<>0 then v:=solve(v);addsol(v);fi;
#  od;
#  W.hxy.sols:=Set(W.hxy.sols);
#  return W.hxy.sols;
#end;
#
#checkennola:=function(W,x,y)local z,dep,arr,e,ew,uc;
#  z:=OrderCenter(W.type[1]);
#  uc:=UnipotentCharacters(W);
#  dep:=hxy(W,x,y);
#  if dep=false then return false;fi;
#  if IsBound(W.Nhyp) then ew:=W.N+W.Nhyp;
#  else ew:=2*W.N;
#  fi;
#  arr:=hxy(W,x,y+1/z);
#  if arr=false then return false;fi;
#  arr:=E(z)^(-ew*x)*arr;
#  arr.v:=List(arr.v,x->Value(x,["q",Mvp("q")*E(z)]));
#  arr.v:=Zip(arr.v,uc.a+uc.A,function(p,n)return
#     p*Mvp("t")^((n-ew)/z);end);
#  e:=Ennola(W)[1].ls;
## e:=e^-1;
#  arr.v:=Permuted(arr.v,e);
#  Print("dep=");Display(dep,rec(cyclicparam:=true));
#  Display(hxy(W,x,y+1/z),rec(cyclicparam:=true));
#  Print("arr=");Display(arr,rec(cyclicparam:=true));
#  return dep-arr;
#end;
#
#checkennolascalar:=function(W,y)local z,s,uc,cn,aA,om,e,f,ew;
#  z:=OrderCenter(W.type[1]);
#  s:=PrincipalSeries(W,y);Hecke(s);
#  if IsBound(W.Nhyp) then ew:=W.N+W.Nhyp;
#  else ew:=2*W.N;
#  fi;
#  uc:=UnipotentCharacters(W);
#  cn:=s.charNumbers;
#  aA:=uc.A+uc.a;aA:=aA{cn};
#  om:=List(CharTable(s.WGL).irreducibles,
#    x->x[PositionRegularClass(s.WGL,modZ(-1/z))]/x[1]);
#  e:=Ennola(W)[1].allscal{cn};
#  f:=Zip(aA,om,function(x,o)return E(z)^(y*(x-3*ew))*o;end);
#  Print("aA=",aA,"om=",List(om,AsRootOfUnity),"\n");
#  if e=f then return "Ok!";fi;
## e:=TransposedMat([e,f]);
## f:=Filtered([1..Length(e)],x->e[x][1]<>e[x][2]);
## e:=TransposedMat(e{f});
## Add(e,f);
#  return [e,f,Hecke(s),y,CharNames(uc,rec(cyclicparam:=true)){s.charNumbers}];
#end;

testdata=Vector{Tuple{PRG,Rational{Int},Rational{Int}}}[]

push!(testdata,map(x->(crg(5,1,1),x...),
                   [(0,0), (0,1//5), (1//5,0), (1//5,1//5)]))

push!(testdata,map(x->(crg(4),x...),
[ ( 0, 0 ), ( 0, 1//2 ), ( 0, 1//3 ), ( 0, 1//4 ), ( 0, 1//6 ), ( 1//2, 0 ), 
  ( 1//2, 1//2 ), ( 1//2, 1//3 ), ( 1//2, 1//4 ), ( 1//2, 1//6 ), ( 1//3, 0 ), 
  ( 1//3, 1//2 ), ( 1//3, 1//3 ), ( 1//3, 1//6 ), ( 1//4, 0 ), ( 1//4, 1//2 ), 
  ( 1//4, 1//4 ), ( 1//6, 0 ), ( 1//6, 1//2 ), ( 1//6, 1//3 ),( 1//6, 1//6 )]))

#W:=ComplexReflectionGroup(24);
#Append(testdata,List(
#[ [ 0, 0 ], [ 0, 1/2 ], [ 0, 1/3 ], [ 0, 1/6 ], [ 0, 1/7 ], [ 0, 1/14 ], 
#  [ 1/2, 0 ], [ 1/2, 1/2 ], [ 1/2, 1/3 ], [ 1/2, 1/6 ], [ 1/2, 1/7 ], 
#  [ 1/2, 1/14 ], [ 1/3, 0 ], [ 1/3, 1/2 ], [ 1/3, 1/3 ], [ 1/3, 1/6 ], 
#  [ 1/6, 0 ], [ 1/6, 1/2 ], [ 1/6, 1/3 ], [ 1/6, 1/6 ], [ 1/7, 0 ], 
#  [ 1/7, 1/2 ], [ 1/7, 1/7 ], [ 1/7, 1/14 ], [ 1/14, 0 ], [ 1/14, 1/2 ], 
#  [ 1/14, 1/7 ], [ 1/14, 1/14 ] ],
# p->Concatenation([W],p)));
#end;
#add();

"""
`redsl2(m)` takes any m in `sl2(ℤ)` and returns a word in `Ω=[1 0;1 1]` and 
`Fo=[0 -1;1 0]` of which it is the product. 
In the word `Ω^i` is encoded as `i` and `Fo` as `0`.
"""
function redsl2(m)
  pushnz!(a,i)=if i!=0 push!(a,i) end
  res=Int[]
  while true
    if m[2,1]==0
      if m[1,1]==1 push!(res,0);m=-m
      else append!(res,[0,0,0]) end
      pushnz!(res,m[1,2]); append!(res,[0,0,0])
      return res
    end
    a=m[1,1]
    c=m[2,1]
    if a==0
      if m[1,2]>0 append!(res,[0,0]);m=-m end
      pushnz!(res,-m[2,2]); append!(res,[0])
      return res
    elseif abs(a)<=abs(c)
      alpha=div(c-mod(c,a),a)
      pushnz!(res,alpha)
      m=[1 0;-alpha 1]*m
    else
      alpha=div(a-mod(a,c),c)
      append!(res,[0,0,0]); pushnz!(res,-alpha)
      m=[0 -1;1 -alpha]*m
    end
  end
end

checksl2(s)=prod(i->i==0 ? [0 -1;1 0] : [1 0;i 1],s) # inverse of redsl2

# decompose a UniChars with coeffs Mvp(q,t) in list monomial(q,t)=>Unichar
function decomp(u)
  l=pairs.(u.v)
  monoms=unique(vcat(map(x->first.(x),l)...))
  toM(map(monoms)do m
    map(u.v)do p
      coefficient(p,m)
    end
   end)
end

# truncate u to ith family
function trunc(u::UniChar,i)
  r=zero(u.v)
  W=u.group
  l=UnipotentCharacters(W).families[i].charNumbers
  r[l]=u.v[l]
  UniChar(W,r)
end
