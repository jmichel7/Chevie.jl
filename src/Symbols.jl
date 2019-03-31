module Symbols
using Gapjm
export ShiftBeta, BetaSet, PartBeta, SymbolPartitionTuple,
LowestPowerGenericDegreeSymbol, HighestPowerGenericDegreeSymbol,
HighestPowerFakeDegreeSymbol, LowestPowerFakeDegreeSymbol,
DefectSymbol, FullSymbol, RankSymbol, SymbolsDefect,
CycPolFakeDegreeSymbol

function ShiftBeta(beta,n)
  if n>=0 return [0:n-1;beta .+ n]
  elseif beta[1:-n]!=[0:-n-1]
        Error( "Cannot shift ", beta, " by ", n, "\n" );
  else return beta[1-n:end].+n
  end
end

BetaSet(alpha)=isempty(alpha) ? alpha : reverse(alpha) .+(0:length(alpha)-1)

PartBeta(l)=filter(x->!iszero(x),reverse(l.-(0:length(l)-1)))

function SymbolPartitionTuple(p,s)
  if p[end] isa Integer
    l= length(p) - 2
    e= l*p[end-1]
  else e=l=length(p)
  end
  if s isa Integer
    if s<0  s=[0;fill(-s,l-1)]
    else    s=[s;zeros(Int,l-1)]
    end
  else s=s[1:l]
  end
  s= map(length, p[1:l]) .- s
  s= maximum(s) .- s
  p= copy(p)
  p[1:l]=map(i->ShiftBeta(BetaSet(p[i]),s[i]),1:l)
  p
end

function LowestPowerGenericDegreeSymbol(p)
  p=FullSymbol(p)
  e=length(p)
  p=sort!(vcat(p...))
  m=length(p)
  sum(p.*(m-1:-1:0))-div(m*(m-e)*(2*m-3-e),12*e)
end

function HighestPowerGenericDegreeSymbol(p)
  p=FullSymbol(p)
  r=RankSymbol(p)
  e=length(p)
  p=sort!(vcat(p...))
  m=length(p)
  if mod(m,e)==1 r=div(e*r*(r+1),2)
  else           r+= div(e*r*(r-1),2)
  end
  r+sum(p.*(0:m-1))-sum(x->div(e*x*(x+1),2),p)-div(m*(m-e)*(2*m-3-e),12*e)
end

function DefectSymbol(s)
  s=FullSymbol(s)
  length(s[1])-length(s[2])
end

function HighestPowerFakeDegreeSymbol(s,p=length(FullSymbol(s)))
  s=FullSymbol(s)
  e=length(s)
  if e==1 d=0
  else d=DefectSymbol(s)
  end
  if !(d in [0,1]) return -1 end
  r=RankSymbol(s)
  if d==1 res=div(e*r*(r+1),2)
  else    res=div(e*r*(r-1),2)+r
  end
  res+=e*sum(S->sum(S.*(0:length(S)-1))-sum(l->div(l*(l+1),2),S),
               filter(x->!isempty(x),s))
  gamma=i->sum(mod.(i+(0:e-1),e).*map(sum,s))
  if d==1 res+=gamma(0)
  else res+=maximum(gamma.(0:e-1))
  end
  res-sum(map(x->div(x*(x-1),2),e*(1:div(sum(length,s),e)-1).+mod(d,2))) 
end

function LowestPowerFakeDegreeSymbol(s)
  s=FullSymbol(s)
  if length(s)==1 d=0
  else d=DefectSymbol(s)
  end
  if !(d in [0,1]) return -1 end
  e=length(s)
  res=e*sum(S->sum(S.*(length(S)-1:-1:0)),filter(x->!isempty(x),s))
  gamma=i->sum(mod.(i+(0:e-1),e).*map(sum,s))
  if d==1 res+=gamma(0)
  else res+=minimum(gamma.(0:e-1))
  end
  res-sum(map(x->div(x*(x-1),2),e*(1:div(sum(length,s),e)-1).+mod(d,2))) 
end

function FullSymbol(S)
  if isempty(S) || S[end] isa AbstractVector return S end
  m=S[end-1]
  vcat(map(i->map(copy,S[1:end-2]),1:m)...)
end

function RankSymbol(s)
  s=FullSymbol(s)
  ss=sum(length,s)
  if isempty(s) return 0 end
  e=length(s)
  sum(sum,s)-div((ss-1)*(ss-e+1),2*e)
end

# e-symbols of rank r, Malle-defect def and content=c (mod e)
# The list is returned sorted by HC series (principal series first)
# SymbolsDefect(d,r,0,1) gives symbols of unipotent characters of G(d,1,r)
# SymbolsDefect(e,r,0,0) gives symbols of unipotent characters of G(e,e,r)
SymbolsDefect=function(e,r,def,c)
  local IsReducedSymbol, S
  function defShape(s)local e
    e=length(s)
    (binomial(e,2)*div(sum(s),e)-sum(i->(i-1)*s[i],1:e))%e
  end
 
  shapesSymbols=function(r,e,c)local f,res,m,new
    f=function(lim2,sum,nb,max)local res,a
      if nb==1   
        if sum==0   return [[sum]]
        else return Vector{Int}[] end 
      end
      res=Vector{Int}[]
      a=div(sum,nb-1)
      while a<=max  &&  binomial(a,2)<=lim2  &&  a<=sum  
        append!(res,map(x->vcat([a],x),f(lim2-binomial(a,2),sum-a,nb-1,a)))
        a+=1 
      end
      return res
    end

    res=[]
    m=0
    while true
      new=f(r+div((m*e+c-1)*(m*e+c-e+1),2*e),c+m*e,e,c+m*e)
      append!(res,new)
      m+=1
      if length(new)==0 break end
    end
    res=vcat(map(x->arrangements(x,e),res)...)
    filter(s->defShape(s)==def  &&  
      all(x->defShape(x)!=def  ||  x<=s,map(i->circshift(s,i),1:length(s)-1)),res)
  end

  function IsReducedSymbol(s)
    ForAll(Rotations(s){[2..Length(s)]},x->s==x || LessSymbols(x,s))
  end

  S=vcat(map(shapesSymbols(r,e,c)) do s
    map(x->SymbolPartitionTuple(x,s),
       partition_tuples(r-RankSymbol(map(x->0:x-1,s)),e)) end...)
  c!=0 ? S : filter(IsReducedSymbol,S)
end

function CycPolFakeDegreeSymbol(arg...)
  local s, p, res, delta, theta, r, d, e, q, rot
  s = FullSymbol(arg[1])
  q = Pol([1],1)
  e = length(s)
  r = RankSymbol(s)
  if e == 0 return CycPol(1) end
  if length(arg) == 2 p = E(e, arg[2])
  else p = 1
  end
  function delta(S)
    S=collect(combinations(S,2))
    if isempty(S) return CycPol(Pol([1],0)) end
    prod(x->CycPol(q^(e*x[2])-q^(e*x[1])),S)
  end
  function theta(S)
    S=filter(x->x>0,S)
    if isempty(S) return CycPol(Pol([1],0)) end
    prod(l->prod(h->CycPol(q^(e*h)-1),1:l),S)
  end
  if length(s) == 1 d = 0
  else d = DefectSymbol(s)
  end
  if d == 1 res = theta([r])
  elseif d == 0 res = theta([r - 1]) * CycPol(q ^ r - p)
  else res = CycPol(0q)
  end
  res*=prod(S->delta(S)//theta(S),s)
  res//=CycPol(q^sum([div(x*(x-1),2) for x in e*(1:div(sum(length,s),e)-1)+d%2]))
  if d == 1
    res*=CycPol(q^sum(map((x,y)->x*sum(y),0:e-1,s)))
  else
      rot = Rotations(s)
      res*=CycPol(map(j->p^j,0:e-1)*
         map(s->q^((0:e-1)*map(Sum,s)), rot)) // count(i->i==s, rot)
      if e == 2 && p == -1
          res = -res
      end
  end
  if r == 2 && (e > 2 && p == E(e))
      res = CycPol(Value(res, E(2e) * q) // E(2e, Degree(res)))
  end
  return res
end

end
