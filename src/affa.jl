"""
# affa.g                             (C) François Digne (Amiens university)
#
# This file contains:
#   -functions for computing with periodic permutations of the integers
#     The main functions are:
#     PPerm to create a periodic permutation.
#     Operations *, /, ^
#     functions Cycles and CycleType
#
#   -the function CoxeterGroupAtildeGroup which creates the Coxeter group
#    of type ~A_n as a group of periodic permutations.
#
#   -the function DualBraidMonoid for the above kind of groups.
# It is based on the papers
# [Digne] Presentations duales pour les groupes de tresses de type affine A
#         Comment. Math. Helv. 81 (2006) 23--47
# [Shi] The Kazhdan-Lusztig cells in certain affine Weyl groups 
#       Springer LNM 1179 (1986) 
#--------------------------------------------------------------------------
"""
#module AffA

"""
Periodic permutation f of the integers
- of period n, which means f(i+n)=f(i)+n
- with no shift, which means sum(f.(1:n))=sum(1:n)
Stored as [f(1),...,f(n)]
"""
struct PPerm{T<:Integer}
  d::Vector{T}
  function PPerm(d::AbstractVector{T};check=false)where T
    if check isvalid(d) end
    new{T}(d)
  end
end

function isvalid(d)
  n=length(d)
  if Set(mod.(d,n))!=Set(0:n-1) error(d,": images must be distinct mod ",n) end
  if sum(d)!=sum(1:n) error(d,": sum of shifts must be 0") end
end

Base.one(p::PPerm)=PPerm(1:length(p.d))
Base.isone(p::PPerm)=p.d==eachindex(p.d)
Base.copy(p::PPerm)=PPerm(copy(p.d))
Base.broadcastable(p::PPerm)=Ref(p)

"""
`PPerm(n,(i₁,…,iₖ)[=>d=0])`  where  `i₁,…,iₖ`  differ  mod  `n` represents the
permutation `i₁↦ i₂↦ …↦ iₖ↦ i₁+d*n`
"""
function PPerm(n::Int,cc...)
  if isempty(cc) return PPerm(1:n) end
  cc=map(cc) do cyc
    if cyc isa Int 
      ((cyc,),0)
    elseif cyc isa Tuple 
      (cyc,0)
    elseif cyc isa Pair 
      if first(cyc) isa Int
        ((first(cyc),),last(cyc))
      else cyc
      end
    end
  end
  u=collect(Iterators.flatten(map(x->mod.(x[1],n),cc)))
  if length(unique(u))!=length(u)
    error(cc," : the cycles must be disjoint mod ",n)
  end
  p=prod(cc)do (cyc,d)
    perm=collect(1:n)
    cyc=cyc.-1
    for i in eachindex(cyc)
      k=1+mod(cyc[i],n)
      perm[k]=cyc[1+mod(i,length(cyc))]+k-cyc[i]
    end
    perm[1+mod(cyc[end],n)]+=d*n
    PPerm(perm)
  end
  isvalid(p.d)
  p
end

function Base.:*(x::PPerm,y::PPerm)
  n=length(x.d)
  ymodn=mod1.(y.d,n)
  PPerm(x.d[ymodn].+y.d.-ymodn)
end

function Base.inv(x::PPerm)
  n=length(x.d)
  l=map(i->i-mod(i,n),x.d-1)
  ll=(1:n).^inv(Perm(x.d.-l))
  PPerm(ll.-l[ll])
end

function Base.:^(i::Int,p::PPerm)
  y=1+mod(i-1,length(p.d))
  i+p.d[y]-y
end

Base.:^(q::PPerm,p::PPerm)= inv(p)*q*p

Base.:^(a::PPerm, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                               Base.power_by_squaring(inv(a),-n)

Base.:/(a::PPerm,b::PPerm)=a*inv(b)
Base.:\(a::PPerm,b::PPerm)=inv(a)*b

#PPermOps.\=:=function(a,b)return a.perm=b.perm;end;

# Non-trivial cycles of a PPerm; each cycle i_1,..,i_k,[d] is normalized
# such that i_1 mod n is the smallest of i_j mod n and i_1 is in [1..n]
function cycles(a::PPerm)
  res=Pair{Vector{Int},Int}[]
  n=length(a.d)
  l=trues(n)
  while true
    x=findfirst(l)
    if isnothing(x) break end
    cyc=Int[]
    while true
      push!(cyc,x)
      xn=mod1(x,n)
      l[xn]=false
      x+=a.d[xn]-xn
      if mod(x-cyc[1],n)==0 break end
    end
    pp=cyc=>div(x-cyc[1],n)
    if length(cyc)>1 || pp[2]!=0 push!(res,pp) end
  end
  res
end    

function Base.show(io::IO,a::PPerm)
  n=length(a.d)
  function stringdec(d)
    if get(io,:sgn,false)
      d>0 ? "₊"^d : "₋"^(-d)
    else
      d==0 ? "" : fromTeX(io,"_{$d}")
    end
  end
  if !get(io,:limit,false) && !get(io,:TeX,false) 
    print(io,"PPerm(",a.d,")");return
  end
  c=cycles(a)
  if isempty(c) print(io,"()")
  else
    for cc in c
      cyc,d=cc
      print(io,"(",join(map(cyc)do y
        y-=1
        x=mod(y,n)
        string(1+x,stringdec(div(y-x,n)))
      end,","),")")
      print(io,stringdec(d))
    end
  end
end

##------------------------AtildeGroup----------------------------
##The following formula is from [Shi] Lemma 4.2.2
Base.length(w::PPerm)=sum(j->sum(i->abs(fld(j^w-i^w,length(w.d))),1:j-1),
                      eachindex(w.d))

isrightdescent(w::PPerm,i)= i==length(w.d) ? w.d[i]>w.d[1]+length(w.d) : 
                     w.d[i]>w.d[1+i]

CoxGroups.isleftdescent(w::PPerm,i)=isrightdescent(w^-1,i)

#PPermOps.FirstLeftDescending:=function(x)local i,IRD;
#  x:=x^-1;IRD:=x.operations.IsRightDescending;
#  for i in [1..Length(x.perm)] do if IRD(x,i) then return i;fi;od;
#  return false;
#end;

# for this function see [Digne],2.8
function Perms.reflength(w::PPerm)
  function vec(pp,i)
    pp=map(x->sum.(getindex.(Ref(pp),x)),partitions_set(eachindex(pp),i))
    if pp[1][1]<0 pp=-pp end
    return tally.(pp)
  end
  d=last.(cycles(w))
  n=length(w.d)
  p0=count(iszero,d)
  res=n+length(d)-2*p0-count(i->i==i^w,1:n)
  if p0==length(d) return res end
  pos=filter(x->x>0,d)
  neg=filter(x->x<0,d)
  m=min(length(pos),length(neg))
  r=m:-1:1
  res-2*r[findfirst(i->length(intersect(vec(pos,i),vec(neg,i)))>0,r)]
end

function refword(w::PPerm)
  n=length(w.d)
  cnt=0
  function ff(w)local s
    l=reflength(w)
    d=1
    while true
      for i in 1:n
        s=PPerm(n,(i,i+d))
        if reflength(s*w)<l return s end
      end
      d+=1
      if mod(d,n)==0 d+=1 end
    end
  end
  res=PPerm[]
  while !isone(w)
    s=ff(w)
    push!(res,s)
    w=s*w
  end
  res
end

# descent sets are encoded as a pair: a list of atoms, and a list
# of atoms of the form [1,u] representing all atoms [1,u+i*n]
# This uses lemma 2.20 of [Digne] and is valid only if there are
# 0 or 2 cycles with a non-zero shift
function dualleftdescents(a::PPerm)
  res=[PPerm[],PPerm[]]
  n=length(a.d)
  for (x,d) in cycles(a)
    if x[1]!=1 || length(x)!=1
      for j in 1:length(x)
        for k in j+1:length(x)
          push!(res[1],PPerm(n,(x[j],x[k])))
          if d!=0 push!(res[1],PPerm(n,(x[j]+n,x[k]))) end
        end
        if d!=0 push!(res[2],PPerm(n,(1,1+mod(x[j]-1,n)))) end
      end
    end
  end
  res
end

# a,b are results of dualleftdescents
function firstintersectiondualleftdescents(a,b)
  for t in a[1]
    if t in b[1] || mod1.(t.d,length(t.d)) in map(x->x.d,b[2]) 
      return t
    end
  end
  for t in b[1]
    if mod1.(t.d,length(t.d)) in map(x->x.d,a[2]) return t end
  end
  for t in a[2] if t in b[2] return t end end
end

#PPermOps.CycleType:=function(a)local res;
#  res:=List(PPermOps.Cycles(a),cyc->[Length(cyc)-1,cyc[Length(cyc)]]);
#  Sort(res);
#  return res;
#end;
##--------------------------------------------------------------------------
## The next function constructs W(~A_{n-1}) as a group of periodic
## permutations of period n.
# 
#AtildeGroupOps:=OperationsRecord("AtildeGroupOps",GroupOps);
#
#AtildeGroupOps.Print:=function(W)Print(W.name);end;
#
#AtildeGroupOps.PrintDiagram:=function(W)
#    CHEVIE.R("PrintDiagram","AffineA")(W.reflectionsLabels);end;
# 
#AtildeGroupOps.IsRightDescending:=function(W,w,i)
#  return PPermOps.IsRightDescending(w,i);
#end;
#
#AtildeGroupOps.FirstLeftDescending:=function(W,x)
#  return PPermOps.FirstLeftDescending(x);
#end;
#
## Reflections (a,b[i]) are enumerated by lexicographical order of [i,a,b-a]
## with i positive --- recall that when a>b this reflection is printed (b,a[-i])
#AtildeGroupOps.Reflection:=function(W,i) local n,p,r,ecart,pos;
#  n:=W.semisimpleRank;
#  p:=QuoInt(i-1,n*(n-1)); r:=(i-1) mod (n*(n-1));
#  ecart:= QuoInt(r,n); pos:=r mod n;
#  return PPerm([1+pos,2+pos+ecart+p*n],n);
#end;

struct Atilde{T} <: CoxeterGroup{PPerm{T}}
  gens::Vector{PPerm{T}}
  prop::Dict{Symbol,Any}
end

Base.show(io::IO,G::Atilde)=print(io,"Atilde(",G.gens,")")
  
function Atilde(n)
  if n<2 error(n," should be >=2") end
  gens=map(i->PPerm((1:n)^Perm(i,i+1)),1:n-1)
  push!(gens,PPerm(vcat([0],2:n-1,[n+1])))
  Atilde(gens,Dict{Symbol,Any}())
end

Base.length(W::Atilde,w)=length(w)
CoxGroups.isleftdescent(W::Atilde,w,i)=isleftdescent(w,i)
Perms.reflength(W::Atilde,w)=reflength(w)

PermRoot.reflection_subgroup(W::Atilde,I)=Atilde(gens(W)[I],Dict{Symbol,Any}())

struct AffaDualBraidMonoid{T,TW}<:Garside.GarsideMonoid{T}
  δ::T
  stringδ::String
  W::TW
  prop::Dict{Symbol,Any}
end

Garside.IntervalStyle(M::AffaDualBraidMonoid)=Garside.Interval()

function Garside.DualBraidMonoid(W::Atilde)
  n=length(gens(W))
  delta=PPerm(vcat([1-n],3:n,[2+n]))
  AffaDualBraidMonoid(delta,"c",W,Dict{Symbol,Any}())
end

function Garside.leftgcd(M::AffaDualBraidMonoid,a,b)
  x=one(M)
  while true
    t=firstintersectiondualleftdescents(dualleftdescents(a),dualleftdescents(b))
    if isnothing(t) return x,(a,b) end
    x*=t
    a=t^-1*a
    b=t^-1*b
  end
end

CoxGroups.word(M::AffaDualBraidMonoid,w)=refword(w)
function CoxGroups.word(io::IO,M::AffaDualBraidMonoid,w)
  join(map(x->sprint(show,x;context=io),refword(w)))
end

Garside.δad(M::AffaDualBraidMonoid,x,i::Integer)=iszero(i) ? x : x^(M.δ^i)

#  W.operations.\in:=function(e,W)return Length(e.perm)=W.rank;end;
#
## DualMonoid(W[,M])
## constructs dual monoid for tilde A_{n-1}
##
## If a second argument is given, constructs the reversed monoid
## [which allows to fill the field M.revMonoid after building the dual monoid]
#
#AtildeGroupOps.DualBraidMonoid:=function(arg)local M,n,W,delta;
#  W:=arg[1];n:=W.rank;
#  if Length(arg)=1 then
#  else delta:=arg[2].delta^-1;
#  fi;
#  M:=rec(rank:=n,
#    delta:=delta,
#    stringDelta:="c",
#    group:=W,
#    identity:=W.identity,
#    FormatSimple:=function(a,opt)local option;
#      option:=ShallowCopy(PPermOps.PrintOptions);Inherit(option,opt);
#      if IsBound(option.word) then return IntListToString(CoxeterWord(W,a));
#      else return Format(a,option);fi;end,
#    Reverse:=function(b)local res,s;
#      if Length(b.elm)=0 then return M.revMonoid.Elt([],b.pd);fi;
#      res:=[];
#      for s in List(Reversed(b.elm),y->M.revMonoid.DeltaAction(y^-1,b.pd))
#      do res:=M.revMonoid.AddToNormal(res,s);od;
#      return GarsideEltOps.Normalize(M.revMonoid.Elt(res,b.pd));
#    end,
#    RightComplementToDelta:=a->a^-1*M.delta,
#    RightAscentSet:=a->M.LeftDescentSet(M.RightComplementToDelta(a)),
#    operations:=rec(Print:=function(M) Print("DualBraidMonoid(",W,")");end)
#  );
#  M.Elt:=function(arg)local res;
#    res:=rec(elm:=arg[1],operations:=GarsideEltOps,monoid:=M);
#    if Length(arg)>1 then res.pd:=arg[2]; else res.pd:=0; fi;
#    return GarsideEltOps.Normalize(res);
#  end;
#  M.B:=function(arg)local x,p;
#    if IsList(arg[1]) then x:=[arg[1]]; else x:=arg; fi;
#    p:=PositionProperty(x,y->not Number(Cycles(y),c->c[Length(c)][1]<>0)in [0,2]);
#    if p<>false then Error(x[p]," is not a dual simple");fi;
#    x:=Concatenation(List(x,M.AtomListSimple));
#    if not ForAll(x,M.IsDualAtom) then
#      Error("not atom of dual monoid: ",First(x,y->not M.IsDualAtom(y)));fi;
#    return GarsideEltOps.Normalize(Product(List(x,y->M.Elt([y]))));
#  end;
#  M.AtomListSimple:=function(w) local s,res,v;res:=[];v:=w;
#    while v<>M.group.identity do
#     s:=Flat(M.LeftDescentSet(v))[1];
#     v:=s^-1*v;
#     Add(res,s);
#    od;
#    return res;
#  end;
#  M.IsDualAtom:=function(a)local c;
#   c:=Cycles(a);
#   return Length(c)=1 and Length(c[1])=3 and c[1][3]=[0] and
#    (c[1][1]=1 or AbsInt(c[1][1]-c[1][2])<M.rank);
#  end;
#  CompleteGarsideRecord(M,rec(interval:=true));
#  if Length(arg)=1 then M.revMonoid:=DualBraidMonoid(W,M);
#  else M.revMonoid:=arg[2];
#  fi;
#  return M;
#end;
#
#AtildeBraid:=n->DualBraidMonoid(CoxeterGroupAtildeGroup(n)).B;
#Atilde:=PPerm;
