# Hand-translated part of chevie/tbl/exceptio.jl
# (C) Jean Michel 1999-2017
# data common to several (but not all) types of reflection groups

# an addition
chevieset(["A","B","D"],:EigenvaluesGeneratingReflections,t->r->fill(1//2,r))

chevieset(["G25","G26","G29","G31","G32","G34"],:CartanMat,
  function(t)
    r=chevieget(t,:GeneratingRoots)
    eig=map(x->Root1(;r=x),chevieget(t,:EigenvaluesGeneratingReflections))
    toL(toM(coroot.(r,eig))*transpose(toM(r)))
  end)

chevieset(["E7", "E8", "H3", "H4"], :Invariants, t->
  function()
    C=chevieget(t, :CartanMat)
    r=roots(C)*C
    map(d->function(arg...)sum(a->sum(arg.*a)^d,r) end, 
        chevieget(t, :ReflectionDegrees))
  end
)

chevieset(["G24","G27","G29","G33","G34","E6","E7","E8","H3","H4"], 
  :FactorizedSchurElement, t->function(ch,para,arg...)
   c=chevieget(t,:CycPolSchurElements)[findfirst(==(ch),chevieget(t,:CharInfo)()[:charparams])]
   q=-para[1][1]//para[1][2]
   res=HeckeAlgebras.FactSchur(Mvp(c[1]*q^Int(c[2])), 
                 map(v->(pol=CycPol([1,0,v]),monomial=q),c[3:length(c)]))
   HeckeAlgebras.simplify(res)
 end
)

function ExpandRep(r, d, l) # decompress representation of r gens of dim d
  T=reduce(promote_type,typeof.(first.(l)))
  m=map(i->map(j->fill(zero(T),d),1:d), 1:r)
  for v in l
    for k in @view v[2:end]
      q,i=divrem(Int(k),d^2)
      q1,r1=divrem(i,d)
      m[q+1][q1+1][r1+1]=v[1]
    end
  end
  return m
end

"""
 EvalPolRoot(pol, x, n, p) compute pol(p*x^(1/n))
  
  The point of this routine is to avoid unnecessary root extractions
  during evaluation (e.g., if pol has no terms of odd degree and n=2,
  then no root extraction is necessary).

 this was in lib/util.g but is used only here
"""
function EvalPolRoot(pol::Pol,x,n,p)
# println("pol=",pol,"\nx=",x,"\nn=",n,"\np=",p)
  if isempty(pol.c) return 0 end
  P=vcat(fill(0,mod(pol.v,n)),pol.c)
  P=map(i->Pol(P[i:n:length(P)],div(pol.v-mod(pol.v,n),n))(x*p^n),1:n)
  j=findlast(!iszero,P)
  if isnothing(j) return 0 end
  pol=Pol(P[1:j],0)
  l=pol.v-1+filter(i->!iszero(pol.c[i]),eachindex(pol.c))
  r=gcd(n,l...)
  pol=Pol(pol.c[1:r:length(pol.c)],div(pol.v,r))
  pol(GetRoot(x,div(n,r))*p^r)
end

#  VcycSchurElement(para,r(schur model)[,data(schur data)])
#
#  This function computes the Schur elements for G4-22,  G25-26, G28, G32
#  according to the data computed by M. Chlouveraki.
#  para is the list of parameters of the algebra.
#  schur model describes the shape of the Schur element: it has the fields
#   .factor=(possibly fractional) vecmonomial
#   .coeff= a constant
#   [nothing] or [.root=vecmonomial] or [.rootUnity]  
#   vcyc= a list of pairs [vecmonomial, cyclotomic polynomial index]
#   rootCoeff=  a constant by which multiply .root before taking root
#  where vecmonomial=vector of powers for elts of para (plus possibly
#     the power to which to raise .root or .rootUnity)
#  schur data describes the Schur element in its Galois orbit : it has fields
#   order: in which order to take the variables
#   rootPower: by which E(root)^i multiply .root
function VcycSchurElement(para,r,data=nothing)
# println("para=",para,"\nr=",r,"\ndata=",data)
  n=length(para)
  if !isnothing(data) para=para[data[:order]] else para = copy(para) end
  monomial(mon)=prod(map(^,para//1,Int.(mon[1:n])))
  if haskey(r, :rootUnity) && haskey(r,:root) error("cannot have both") end
  if haskey(r, :coeff) res = r[:coeff] else res = 1 end
  if haskey(r, :factor) res*=monomial(r[:factor]) 
     res=Pol([res],0)
  end
  function term(v)
    mon,pol=v
    if haskey(r,:rootUnity)
      tt=monomial(mon)
      if length(mon)==n+1 tt*=(r[:rootUnity]^data[:rootUnityPower])^mon[n+1] end
      Pol([cyclotomic_polynomial(pol)(tt)],0)
    elseif haskey(r, :root)
     if length(mon)==n return Pol([cyclotomic_polynomial(pol)(monomial(mon))],0)
     else return cyclotomic_polynomial(pol)(Pol([monomial(mon)],mon[n+1]))
     end
    else 
     Pol([cyclotomic_polynomial(pol)(monomial(mon))],0)
    end
  end
  res*=prod(term.(r[:vcyc]))
  if !haskey(r, :root) return res.c[1] end
  den=lcm(denominator.(r[:root])...)
  root=monomial(den*r[:root])
  if haskey(r, :rootCoeff) root*=r[:rootCoeff] end
  EvalPolRoot(res, root, den, data[:rootPower])
end

"""
`BDSymbols(n,d)`
    
returns  2-symbols of defect `d` and rank `n` (for Weyl types B,C,D,2D). If
`d==0`  the symbols with  equal entries are  returned twice, represented as
the  first entry, followed by the repetition factor 2 and an ordinal number
0 or 1, so that `BDSymbols(n, 0)` is a set of parameters for the characters
of the Weyl group of type `Dₙ`.

```julia-repl
julia> GAPENV.BDSymbols(2,1)
5-element Vector{Vector{Vector{Int64}}}:
 [[1, 2], [0]]
 [[0, 2], [1]]
 [[0, 1, 2], [1, 2]]
 [[2], []]
 [[0, 1], [2]]

julia> GAPENV.BDSymbols(4,0)
13-element Vector{Vector{T} where T}:
 Any[[1, 2], 2, 0]
 Any[[1, 2], 2, 1]
 [[0, 1, 3], [1, 2, 3]]
 [[0, 1, 2, 3], [1, 2, 3, 4]]
 [[1, 2], [0, 3]]
 [[0, 2], [1, 3]]
 [[0, 1, 2], [1, 2, 4]]
 Any[[2], 2, 0]
 Any[[2], 2, 1]
 [[0, 1], [2, 3]]
 [[1], [3]]
 [[0, 1], [1, 4]]
 [[0], [4]]
```
"""
function BDSymbols(n,d)
  n-=div(d^2,4)
  if n<0 return Vector{Vector{Int}}[] end
  if d>0 return map(x->symbol_partition_tuple(x,d),partition_tuples(n,2)) end
   return map(chevieget(:D,:symbolcharparam),
              chevieget(:imp,:CharInfo)(2,2,n)[:charparams])
end

bar="\u2014"
rdarrow(n)="\u21D0"^(n-1)*" "
ldarrow(n)="\u21D2"^(n-1)*" "
tarrow(n)="\u21DB"^(n-1)*" "
dbar="\u2550"
#dbar="="
tbar(n)="\u2261"^(n-1)*" "
vbar="\UFFE8" # "\u2503"
c2="\u2461" # "\u2B55"
c3="\u2462"
c4="\u2463"
c5="\u2464"
node="O"

function getlind(d)
  t=d.t
  indices=t.indices
  if isnothing(indices) ind=fill("?",rank(t))
  else ind=repr.(indices)
  end
  length.(ind),ind
end

function Base.show(io::IO,d::Diagram,::Val{:A})
  l,ind=getlind(d)
  join(io,node.*bar.^l[1:end-1]);print(io,node);println(io," ",d.t)
  join(io,ind," ")
end

function Base.show(io::IO,d::Diagram,::Val{:B})
  l,ind=getlind(d)
  c=d.t.cartanType
  print(io,node)
  if c==2 l1=max(l[1],2);print(io,rdarrow(l1))
  elseif c==1 l1=max(l[1],2);print(io,ldarrow(l1))
  elseif c==root(2) l1=max(l[1],2);print(io,dbar^l1)
  else xprint(io,"=",c,"=");
     l1=length(repr(c;context=rio()))+2
  end
  join(io,node.*bar.^l[2:end-1]);println(io,node," ",d.t)
  print(io,rpad(ind[1],l1+1));join(io,ind[2:end]," ")
end

function Base.show(io::IO,d::Diagram,::Val{:D})
  l,ind=getlind(d)
  println(io," "^(l[1]+1),node," $(ind[2])")
  println(io," "^(l[1]+1),vbar)
  print(io,node,map(l->bar^l*node,l[[1;3:end-1]])...)
  println(io," ",d.t)
  join(io,ind[[1;3:end]]," ")
end

function Base.show(io::IO,d::Diagram,::Val{:E})
  l,ind=getlind(d)
  println(io," "^(2+l[1]+l[3]),node," $(ind[2])")
  println(io," "^(2+l[1]+l[3]),vbar)
  print(io,node,map(l->bar^l*node,l[[1;3:end-1]])...)
  println(io," ",d.t)
  join(io,ind[[1;3:end]]," ")
end

function Base.show(io::IO,d::Diagram,::Val{:F})
  l,ind=getlind(d)
  print(io,node,bar^l[1],node)
  if d.t.cartanType==1 l1=max(l[2],2);print(io,ldarrow(l1))
  else l1=max(l[2],1);print(io,dbar^l1)
  end
  println(io,node,bar^l[3],node," ",d.t)
  print(io,ind[1]," ",rpad(ind[2],l1+1),ind[3]," ",ind[4])
end

function Base.show(io::IO,d::Diagram,::Val{:G})
  l,ind=getlind(d)
  if d.t.cartanType==1 print(io,node,tarrow(max(l[1],2)))
  else print(io,node,tbar(max(l[1],2)))
  end
  println(io,node," ",d.t)
  print(io,ind[1]," "^max(3-l[1],1),ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:H})
  l,ind=getlind(d)
  println(io," "^l[1],"₅")
  println(io,map(i->node*bar^l[i],1:length(l)-1)...,node," ",d.t)
  join(io,ind," ")
end

function Base.show(io::IO,d::Diagram,::Val{:I})
  l,ind=getlind(d)
  println(io," "^l[1],d.t.bond)
  println(io,node,bar^l[1],node)
  join(io,ind," ")
end

function Base.show(io::IO,d::Diagram,::Val{:G4})
  l,ind=getlind(d)
  println(io,c3," ",bar^2,c3," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G5})
  l,ind=getlind(d)
  println(io,c3," ",dbar^2,c3," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G6})
  l,ind=getlind(d)
  println(io,c2," ",tbar(2),c3," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G7})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^2,c3," ",ind[2]," ",d.t);
  println(io," ","/3\\")
  println(io,c2," ",bar^2,c3)
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3),"=", f(2,3,1),"=",f(3,1,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G8})
  l,ind=getlind(d)
  println(io,c4," ",bar^2,c4," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G9})
  l,ind=getlind(d)
  println(io,c2," ",tbar(2),c4," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G10})
  l,ind=getlind(d)
  println(io,c3," ",dbar^2,c4," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G11})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^2,c3," ",ind[2]," ",d.t);println(io," /3\\");
  println(io,c2," ",bar^2,c4)
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3),"=",f(2,3,1),"=",f(3,1,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G12})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",c2," ",ind[2]," ",d.t);
  println(io," /4\\");println(io,c2," ",bar^2,c2)
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3,1),"=",f(2,3,1,2),"=",f(3,1,2,3))
end

function Base.show(io::IO,d::Diagram,::Val{:G13})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",c2," ",ind[1]," ",d.t)
  println(io," / \\");println(io,c2," ",bar^2,c2)
  print(io,ind[3]," "^3,ind[2]," ",f(2,3,1,2),"=",f(3,1,2,3)," ",f(1,2,3,1,2),"=")
  print(io,f(3,1,2,3,1))
end

function Base.show(io::IO,d::Diagram,::Val{:G14})
  l,ind=getlind(d)
  println(io,"  ₈");println(io,c2," ",bar,c3," ",d.t)
  print(io,ind[1]," "^2,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G15})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",c2," ",ind[1]," ",d.t);
  println(io," /5");println(io,c2," ",ind[3]);
  println(io," \\");print(io,"  ",c3," ",ind[2]," ")
  print(io,f(1,2,3),"=",f(3,1,2)," ",f(2,3,1,2,1),"=",f(3,1,2,1,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G16})
  l,ind=getlind(d)
  println(io,c5," ",bar^2,c5," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G17})
  l,ind=getlind(d)
  println(io,c2," ",tbar(2),c5," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G18})
  l,ind=getlind(d)
  println(io,c3," ",dbar^2,c5," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G19})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^2,c3," ",ind[2]," ",d.t);println(io," /3\\");
  println(io,c2," ",bar^2,c5)
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3),"=",f(2,3,1),"=",f(3,1,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G20})
  l,ind=getlind(d)
  println(io,"  ₅");println(io,c3," ",bar,c3," ",d.t)
  print(io,ind[1]," "^2,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G21})
  l,ind=getlind(d)
  println(io,"  ₁₀");
  println(io,c2," ",bar^2,c3," ",d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G22})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",c2," ",ind[2]," ",d.t);println(io," /5\\")
  println(io,c2," ",bar^2,c2)
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3,1,2),"=",f(2,3,1,2,3),"=",f(3,1,2,3,1))
end

function Base.show(io::IO,d::Diagram,::Val{:G24})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",c2," ",ind[1]," ",d.t)
  println(io," / \\")
  print(io,c2," ",dbar^2,c2)
  println(io,"  ",f(2,3,1,2,3,1,2,3,1),"==",f(3,2,3,1,2,3,1,2,3))
  print(io,ind[2]," "^3,ind[3])
end

function Base.show(io::IO,d::Diagram,::Val{:G25})
  l,ind=getlind(d)
  println(io,c3," ",bar^2,c3," ",bar^2,c3," ",d.t)
  print(io,ind[1]," "^3,ind[2]," "^3,ind[3])
end

function Base.show(io::IO,d::Diagram,::Val{:G26})
  l,ind=getlind(d)
  println(io,c2," ",dbar^2,c3," ",bar^2,c3," ",d.t)
  print(io,ind[1]," "^3,ind[2]," "^3,ind[3])
end

function Base.show(io::IO,d::Diagram,::Val{:G27})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",c2," ",ind[1]," ",d.t)
  println(io, " / \\")
  print(io,c2," ",dbar^2,c2)
  println(io,"  ",f(3,2,3,1,2,3,1,2,3,1,2,3),"==",f(2,3,1,2,3,1,2,3,1,2,3,2))
  print(io,ind[2]," "^3,ind[3])
end

function Base.show(io::IO,d::Diagram,::Val{:G29})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^6,c2," ",f(4)," "^2,d.t)
  println(io,"     /\u2016\\")
  print(io,c2," ",bar^2,c2," ",dbar^2,c2)
  println(io,"  ",f(4, 3, 2, 4, 3, 2),"==",f(3, 2, 4, 3, 2, 4))
  print(io,f(1)," "^3,f(2)," "^3,f(3))
end

function Base.show(io::IO,d::Diagram,::Val{:G31})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,ind[4], bar^3, ind[2], bar^3, ind[5]," "^3,d.t)
  println(io," \\ /3\\ /")
  print(io,"  ",ind[1],bar^3,ind[3],"     i.e. A₅ on ")
  println(io,f(1,4,2,5,3)," plus ",f(1,2,3),"==",f(2,3,1),"==",f(3,1,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G32})
  l,ind=getlind(d)
  println(io,c3," ",bar^2,c3," ",bar^2,c3," ",bar^2,c3," ",d.t)
  print(io,ind[1]," "^3,ind[2]," "^3,ind[3]," "^3,ind[4])
end

function Base.show(io::IO,d::Diagram,::Val{:G33})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^6, ind[3]," "^7,d.t)
  println(io," "^5,"/^\\")
  print(io,ind[1],bar^3,ind[2],bar^3,ind[4],bar^3,ind[5])
  println(io," ",f(4,2,3,4,2,3),"==",f(3,4,2,3,4,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G34})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^6,ind[3]," "^11,d.t)
  println(io," "^5,"/^\\")
  print(io,ind[1],bar^3,ind[2],bar^3,ind[4],bar^3,ind[5],bar^3,ind[6])
  println(io," ",f(4,2,3,4,2,3),"==",f(3,4,2,3,4,2))
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("A",Char(0x00303))})
  v=d.t.indices
  r=length(v)-1
  if r==1 println(io,d.t.series,"₁  ",v[1]," ∞ ",v[2])
  else 
    n=string(d.t.series,stringind(io,r),"   ")
    s=string(join(v[1:r],bar^3))
    o=length(s)-4
    println(io," "^length(n)," ",bar^div(o,2),v[r+1],bar^div(o,2))
    println(io," "^length(n),"/"," "^o,"\\")
    println(io,n,s)
  end
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("D",Char(0x00303))})
  v=d.t.indices
  r=length(v)-1
  lr=length(string(r))
  println(io,d.t.series,stringind(io,r)," ",v[1]," "^(4*(r-3)-1),v[r+1])
  println(io," "^(lr+3),"\\"," "^(1+4*(r-4)),"/")
  print(io," "^(lr+4),v[3])
  for i in 4:r-1 print(io,bar^3,v[i]) end
  println(io)
  println(io," "^(lr+3),"/"," "^(1+4*(r-4)),"\\")
  println(io," "^(lr+2),v[2]," "^(4*(r-3)-1),v[r])
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("E",Char(0x00303))})
  v=d.t.indices
  if length(v)==7
    println(io,"            ",v[7],"\n            |\n",
          "            ",v[2],"\n            |\n",
          d.t.series,"₆  ",join(v[[1,3,4,5,6]],bar^3))
  elseif length(v)==8
    println(io,"                ",v[2],"\n                |")
    println(io,d.t.series,"₇  ",join(v[[8,1,3,4,5,6,7]],bar^3))
  elseif length(v)==9
    println("            ",v[2],"\n            |")
    println(io,d.t.series,"₈  ",join(v[[1,3,4,5,6,7,8,9]],bar^3))
  end
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("F",Char(0x00303))})
  v=d.t.indices
  println(io,d.t.series,"₄  ",v[5],bar^3,v[1],bar^3,v[2],ldarrow(2),v[3],bar^3,v[4])
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("G",Char(0x00303))})
  v=d.t.indices
  println(io,d.t.series,"₂  ",v[3],bar^3,v[1],tarrow(2),v[2])
end
