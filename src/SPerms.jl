"""
A  signed permutation of `1:n` is  a permutation of the set `-n,…,-1,1,…,n`
which  preserves the  pairs `(-i,i)`.  It is  represented internally as the
images of `1:n`. It is printed as a product of signed cycles.

# Examples
```julia-repl
julia> SPerm([-2,-1,-3])
SPerm{Int64}: (1,-2)(3,-3)

julia> p=SPerm(-1)
(1,-1)

julia> q=SPerm(1,2)
(1,2)

julia> elements(Group([p,q]))
8-element Array{SPerm{Int16},1}:
 ()          
 (1,-1)(2,-2)
 (1,-2,-1,2) 
 (1,-2)      
 (1,2)       
 (1,2,-1,-2) 
 (2,-2)      
 (1,-1)      
```

The  complete type of signed permutations is `SPerm{T}` where `T<:Integer`,
where `Vector{T}` is the type of the vector which holds the image of `1:n`.
This  can used to save space or time when possible. If `T` is not specified
we  take it to be `Int16` since this is a good compromise between speed and
compactness.

SPerms have methods `copy, hash, ==, cmp, isless` (total order) so they can
be  keys in hashes or elements of sets; two `SPerms` are equal if they move
the same points to the same images. For instance,
```julia-repl
julia> SPerm([-2,-1,-3])==SPerm([-2,-1,-3,4])
true
```
SPerms are considered as scalars for broadcasting.
"""
module SPerms
using Gapjm
export SPerm, CoxHyperoctaedral

"""
`struct SPerm`

An  `SPerm` represents a signed permutation of `1:n`, that is a permutation
of  the  set  `-n,…,-1,1,…,n`  which  preserves  the  pairs `(-i,i)`. It is
implemented  by a `struct SPerm`  with one field `d`,  a vector holding the
images of `1:n`.
"""
struct SPerm{T<:Integer}
   d::Vector{T}
end

"""
SPerm{T}(x::Integer...)where T<:Integer

returns   a   signed   cycle.   For  instance  `SPerm{Int8}(1,2,-1,2)`  and
`SPerm({Int8}[1,-2])`  define  the  same  signed  permutation. If not given
`{T}` is taken to be `{Int16}`.
"""
function SPerm{T}(x::Integer...)where T<:Integer
  if isempty(x) return SPerm(T[]) end
  d=T.(1:max(abs.(x)...))
  for i in 1:length(x)-1
   d[abs(x[i])]=sign(x[i])*x[i+1]
  end
  if length(x)==1 d[abs(x[1])]=x[1]
  else d[abs(x[end])]=sign(x[end])*x[1]
  end
  SPerm(d)
end

SPerm(x::Integer...)=SPerm{Int16}(x...)

struct SPermGroup{T}<:Group{SPerm{T}}
  gens::Vector{SPerm{T}}
  prop::Dict{Symbol,Any}
end

function Groups.Group(a::AbstractVector{SPerm{T}}) where T
  SPermGroup(filter(!isone,a),Dict{Symbol,Any}())
end

Base.one(p::SPerm)=SPerm(empty(p.d))

Base.vec(a::SPerm)=a.d

"`Perm(p::SPerm)` returns the underlying Perm of an SPerm"
Perms.Perm(p::SPerm)=Perm(abs.(p.d))

# SPerms are scalars for broadcasting"
Base.broadcastable(p::SPerm)=Ref(p)

# hash is needed for using SPerms in Sets/Dicts
function Base.hash(a::SPerm, h::UInt)
  b = 0x595dee0e71d271d0%UInt
  for (i,v) in enumerate(a.d)
    if v!=i
      b = xor(b,xor(hash(v, h),h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
    end
  end
  b
end

# total order is needed to use SPerms in sorted lists
function Base.cmp(a::SPerm, b::SPerm)
  da=length(a.d)
  db=length(b.d)
  for i in 1:min(da,db)
@inbounds if a.d[i]<b.d[i] return -1 end
@inbounds if a.d[i]>b.d[i] return  1 end
  end
  if     da<db for i in (da+1:db) b.d[i]==i || return -1 end
  elseif da>db for i in (db+1:da) a.d[i]==i || return  1 end
  end
  0
end

Base.isless(a::SPerm, b::SPerm)=cmp(a,b)==-1

Base.:(==)(a::SPerm, b::SPerm)= cmp(a,b)==0


"""
  orbit(a::SPerm,i::Integer) returns the orbit of a on i
"""
function Groups.orbit(a::SPerm{T},i::Integer)where T
  if abs(i)>length(a.d) return T[i] end
  res=T[]
  sizehint!(res,length(a.d))
  j=i
  while true
    push!(res,j)
@inbounds j=sign(j)*a.d[abs(j)]
    if j==i return res end
  end
end

# argument SignedPerm
function Perms.cycles(p::SPerm)
  cycles=Vector{eltype(p.d)}[]
  to_visit=trues(length(p.d))
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    cyc=orbit(p,i)
    to_visit[abs.(cyc)].=false
    if length(cyc)>1 push!(cycles,cyc) end
  end
  cycles
end

" order(a) is the order of the permutation a"
order(a::SPerm) = lcm(length.(cycles(a)))

function Base.show(io::IO, a::SPerm)
  cyc=cycles(a)
  if isempty(cyc) print(io,"()")
  else for c in cyc print(io,"(",join(c,","),")") end
  end
end

function Base.show(io::IO, ::MIME"text/plain", p::SPerm{T})where T
  if T!=Int16 print(io,typeof(p),": ") end
  show(io,p)
end

" `promote(a::Perm, b::Perm)` promotes `a` and `b` to the same degree"
function Base.promote(a::SPerm,b::SPerm)
  da=length(a.d)
  db=length(b.d)
  if da<db
    resize!(a.d,db)
@inbounds    a.d[da+1:db]=da+1:db
  elseif db<da
    resize!(b.d,da)
@inbounds    b.d[db+1:da]=db+1:da
  end
  (a,b)
end

function Base.:*(a::SPerm, b::SPerm)
  a,b=promote(a,b)
  r=similar(a.d)
  for (i,v) in enumerate(a.d) 
@inbounds if v<0 r[i]=-b.d[-v] else r[i]=b.d[v] end
  end
  SPerm(r)
end

function Base.inv(a::SPerm)
  r=similar(a.d)
  for (i,v) in enumerate(a.d) 
 @inbounds if v<0 r[-v]=-i  else r[v]=i end
  end
  SPerm(r)
end

Base.:^(a::SPerm, b::SPerm)=inv(a)*b*a

@inline function Base.:^(n::Integer, a::SPerm{T}) where T
  if abs(n)>length(a.d) return T(n) end
@inbounds n<0 ? -a.d[-n] : a.d[n]
  end

Base.:^(a::SPerm, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                               Base.power_by_squaring(inv(a),-n)

"""
`Base.:^(l::AbstractVector,p::SPerm)`

returns `l` permuted by `p`, a vector `r` such that `r[abs(i^p)]=l[i]*sign(i^p)`.

# Examples
```julia-repl
julia> p=SPerm([-2,-1,-3])
SPerm{Int64}: (1,-2)(3,-3)

julia> [20,30,40]^p
3-element Array{Int64,1}:
 -30
 -20
 -40
```
"""
function Base.:^(l::AbstractVector,a::SPerm)
  res=copy(l)
  for i in eachindex(l)
    v=i^a
    if v>0 res[v]=l[i] else res[-v]=-l[i] end
  end
  res
end

"""
  `SPerm{T}(l::AbstractVector,l1::AbstractVector)`

  return a `SPerm{T}` `p` such that `l1^p==l` if such `p` exists;
  returns nothing otherwise. If not given `{T}` is taken to be `{Int16}`.
  Needs the objects in `l` and `l1` to be sortable.
"""
function SPerm{T}(a::AbstractVector,b::AbstractVector)where T<:Integer
  p=Perm(map(x->sort([x,-x]),a),map(x->sort([x,-x]),b))
  if isnothing(p) return p end
  res=collect(eachindex(a))^p
  for i in eachindex(a)
    if b[i^(p^-1)]!=a[i] res[i]=-res[i] end
  end
  SPerm{T}(res)
end

SPerm(l::AbstractVector,l1::AbstractVector)=SPerm{Int16}(l,l1)

"""
`Matrix(a::SPerm)` is the permutation matrix for a
# Examples
```julia-repl
julia> Matrix(SPerm([-2,-1,-3]))
3×3 Array{Int64,2}:
  0  -1   0
 -1   0   0
  0   0  -1
```
"""
function Base.Matrix(a::SPerm,n=length(a.d))
  res=zeros(Int,n,n)
  for (i,v) in enumerate(a.d) res[i,abs(v)]=sign(v) end
  res
end

#SignedPermOps.Comm:=function(p,q)return p^-1*q^-1*p*q;end;
#
##underlying Signs of a SignedPerm
#Signs:=p->Permuted(List(p.l,SignInt),Perm(p));
#
## We have the properties p=SignedPerm(Perm(p),Signs(p)) and
## if N=OnMatrices(M,p) then
##    M=OnMatrices(N,Perm(p))^DiagonalMat(Signs(p)))
#
## Transforms matrix or perm+signs to a signed perm,
#SignedPerm:=function(arg)local ls,n,i,nz,res;
#  ls:=arg[1];
#  if IsMat(ls) then
#    n:=Length(ls);
#    res:=[];
#    for i in [1..n] do
#      nz:=Filtered([1..n],x->ls[i][x]<>0);
#      if Length(nz)<>1 then return false;fi;
#      nz:=nz[1];
#      if ls[i][nz]=1 then res[i]:=nz;
#      elif ls[i][nz]=-1 then res[i]:=-nz;
#      else return false;
#      fi;
#    od;
#    return SignedPerm(res);
#  elif Length(arg)=2 then n:=arg[2];  # perm, signs
#    return SignedPerm(Permuted(Zip([1..Length(n)],n,
#       function(x,y)return x*y;end),ls^-1));
#  fi;
#end;
#
## duplicate lines and cols of M so HOgroup operates
#SignedPermOps.dup:=function(M)local res,i,j;
#  res:=List([1..2*Length(M)],i->[1..2*Length(M)]*0);
#  for i in [1..Length(M)] do for j in [1..Length(M)] do
#    res{[2*i-1,2*i]}{[2*j-1,2*j]}:=[[M[i][j],-M[i][j]],[-M[i][j],M[i][j]]];
#  od;od;
#  return res;
#end;
#
## SignedMatStab(M [,extra]) find permutations with signs stabilizing M
## (and such that the associated permutations in addition stabilizes extra
##  -- which could be for instance a list of eigenvalues)
#SignedMatStab:=function(arg)local blocks,stab,g,r,gens,I,p,ss,M,extra,k,gr;
#  M:=arg[1];k:=Length(M);
#  ss:=x->Set([x,-x]);
#  if M<>TransposedMat(M) then Error("M should be symmetric");fi;
#  if Length(arg)>1 then extra:=arg[2];else extra:=List(M,x->1);fi;
#  blocks:=CollectBy([1..k],function(i)local inv;
#    inv:=[Collected(List(M[i],ss)),M[i][i]];
#    if IsBound(extra) then Add(inv,extra[i]);fi;
#    return inv;end);
#  g:=Group(PermList([]));I:=[];
#  for r in blocks do
#    if Length(r)>5 then InfoChevie("#I Large Block:",r,"\n");fi;
#    gr:=MatStab(CoxeterGroupHyperoctaedralGroup(Length(r)),SignedPermOps.dup(M{r}{r}));
#    p:=MappingPermListList([1..Length(r)],r);
#    gens:=List(gr.generators,x->SignedPerm(x,Length(r))^p);
#    g:=ApplyFunc(Group,Concatenation(g.generators,gens));
#    Append(I,r); 
#    p:=MappingPermListList(I,[1..Length(I)]);
#    g:=ApplyFunc(Group,List(g.generators,x->x^p));
#    g:=Group(List(g.generators,SignedPermOps.HO),());
#    g:=MatStab(g,SignedPermOps.dup(M{I}{I}));
#    g:=Group(List(g.generators,x->SignedPerm(x,k)),SignedPerm([]));
#    g:=ApplyFunc(Group,List(g.generators,x->x^(p^-1)));
#  od;
#  return g;
#end;
#
## SignedPermMatMat(M, N[, extra1, extra2]) find p such that PsOnMatrices(M,p)=N
## [and such that Permuted(extra1,p)=extra2]
#SignedPermMatMat:=function(arg)
#  local ind,l,I,J,r,p,e,g,n,h,trans,tr,q,ss,M,N,extra1,extra2,PsOnMatrices;
#  PsOnMatrices:=function(M,p)return OnMatrices(M,SignedPerm(p,Length(M)));end;
#  M:=arg[1];N:=arg[2];
#  ss:=x->Set([x,-x]);
#  if M<>TransposedMat(M) then Error("M should be symmetric");fi;
#  if N<>TransposedMat(N) then Error("N should be symmetric");fi;
#  if Length(arg)=2 then extra1:=List(M,x->1);extra2:=List(N,x->1);fi;
#  ind:=function(I,J)local iM,iN,p,n;
#    iM:=List(I,function(i)local inv;
#      inv:=[Collected(List(M[i]{I},ss)),M[i][i]];
#      if IsBound(extra1) then Add(inv,extra1[i]);fi;
#      return inv;end);
#    iN:=List(J,function(i)local inv;
#      inv:=[Collected(List(N[i]{J},ss)),N[i][i]];
#      if IsBound(extra2) then Add(inv,extra2[i]);fi;
#      return inv;end);
#    if Collected(iM)<>Collected(iN) then 
#       InfoChevie("content differs");return false;fi;
#    iM:=CollectBy(I,iM); iN:=CollectBy(J,iN);
#    if Length(iM)=1 then
#      if Length(I)>6 then InfoChevie("large block:",Length(I),"\n");
#        p:=DistHelpedRepresentativeOperation(Group(Reflections(CoxeterGroupHyperoctaedralGroup(Length(I))),()),
#           M{I}{I},N{J}{J},PsOnMatrices,
#	   function(M,N)return Sum(M-N,x->Number(x,y->y<>0));end);
#      else p:=RepresentativeOperation(CoxeterGroupHyperoctaedralGroup(Length(I)),
#         M{I}{I},N{J}{J},PsOnMatrices);
#      fi;
#      if p=false then InfoChevie("could not match block");return false;fi;
#      return [[I,J,SignedPerm(p,Length(I))]];
#    else p:=Zip(iM,iN,ind);
#      if false in p then return false;
#      else return Concatenation(p);
#      fi;
#    fi;
#  end;
#  l:=ind([1..Length(M)],[1..Length(N)]);
#  if l=false then return false;fi;
#  I:=[];J:=[];g:=Group(SignedPerm([]));tr:=SignedPerm([]);
#  for r in l do
##   Print("r=",r,"\n");
#    n:=Length(r[1]);
#    q:=MappingPermListList([1..Length(I)],[1..Length(I)]);
#    p:=MappingPermListList([1..n],[1..n]+Length(I));
#    Append(I,r[1]);Append(J,r[2]);
##   Print("#I=",Length(I),"\c");
#    if Comm(r[3]^p,tr^q)<>SignedPerm([]) then Error("noncomm");fi;
#    tr:=tr^q*r[3]^p;
#    h:=OnTuples(SignedMatStab(M{r[1]}{r[1]}).generators,p);
#    g:=Group(Concatenation(OnTuples(g.generators,q),h),SignedPerm([]));
##   Print(" #g=",Size(g),"\c");
#    e:=RepresentativeOperation(g,M{I}{I},
#      OnMatrices(N{J}{J},tr^-1),OnMatrices);
#    if e=false then return false;
#    else if e^-1*e^tr<>SignedPerm([]) then 
#            Print("*** tr does not commute to e\n");
#         fi;
#         tr:=e*tr;
#    fi;
#    g:=MatStab(Group(List(g.generators,SignedPermOps.HO),()),SignedPermOps.dup(M{I}{I}));
##   Print(" #stab=",Size(g),"\n");
#    g:=Group(List(g.generators,x->SignedPerm(x,Length(I))),SignedPerm([]));
#  od;
#  # transporter of a ps from [1..Length(I)] to I
#  trans:=I->SignedPerm(ListPerm(MappingPermListList([1..Length(I)],I)));
#  return trans(I)^-1*tr*trans(J);
#------------ Example II: HyperOctaedral groups as Coxeter groups

struct CoxHyperoctaedral{T} <: CoxeterGroup{SPerm{T}}
  G::SPermGroup{T}
  n::Int
  prop::Dict{Symbol,Any}
end

Base.iterate(W::CoxHyperoctaedral,r...)=iterate(W.G,r...)

"""
  `CoxHyperoctaedral(n)` The Hyperoctaedral group on ±1,…,±n as a Coxeter group
  of type B, with generators (1,-1) and (i,i+1)(-i,-i-1)
```julia-repl
julia> elements(CoxHyperoctaedral(2))
8-element Array{SPerm{Int8},1}:
 ()          
 (1,2)       
 (1,-1)      
 (1,2,-1,-2) 
 (1,-2,-1,2) 
 (2,-2)      
 (1,-2)      
 (1,-1)(2,-2)
```
"""
function CoxHyperoctaedral(n::Int)
  gens=[SPerm{Int8}(1,-1)]
  for i in 2:n push!(gens,SPerm{Int8}(i-1,i)) end
  CoxHyperoctaedral{Int8}(Group(gens),n,Dict{Symbol,Any}())
end

function Base.show(io::IO, W::CoxHyperoctaedral)
  print(io,"CoxHyperoctaedral($(W.n))")
end

PermRoot.refltype(W::CoxHyperoctaedral)=[TypeIrred(Dict(:series=>:B,
                                        :indices=>collect(1:W.n)))]

CoxGroups.nref(W::CoxHyperoctaedral)=length(gens(W))^2

Gapjm.degrees(W::CoxHyperoctaedral)=2*(1:length(gens(W)))

Base.length(W::CoxHyperoctaedral)=prod(degrees(W))

function CoxGroups.isleftdescent(W::CoxHyperoctaedral,w,i::Int)
  if i^w<0 return true end
  if i==1 return false end
  return (i-1)^w>0 && (i-1)^w>i^w
end

function PermRoot.reflection(W::CoxHyperoctaedral,i::Int)
  ref=gets(W,:reflections)do W
    refs=vcat(gens(W),map(i->SPerm{Int8}(i,-i),2:W.n))
    for i in 2:W.n-1 append!(refs,map(j->SPerm{Int8}(j,j+i),1:W.n-i)) end
    for i in 1:W.n-1 append!(refs,map(j->SPerm{Int8}(j,-j-i),1:W.n-i)) end
    refs
  end
  ref[i]
end

PermRoot.reflections(W::CoxHyperoctaedral)=reflection.(Ref(W),1:nref(W))

function Perms.reflength(W::CoxHyperoctaedral,w)
  c=cycles(w)
  c1=filter(x->Set(x)==Set(-x),c)
  div(reduce(+,map(length,c1);init=0),2)+
      reduce(+,map(x->length(x)-1,setdiff(c,c1));init=0)
end

" Only parabolics defined are I=1:m for m≤n"
function PermRoot.reflection_subgroup(W::CoxHyperoctaedral,I::AbstractVector{Int})
  if I!=1:length(I) error(I," should be 1:n for some n") end
  CoxHyperoctaedral(Group(gens(W)[I]),length(I),Dict{Symbol,Any}())
end

end;
