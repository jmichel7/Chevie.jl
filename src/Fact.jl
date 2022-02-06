module Fact
import Primes: nextprime
using LaurentPolynomials: Pol, @Pol, shift, degree, derivative, exactdiv
using ..FFields: FFields, FFE, Mod
using ..Util: Util, factor
using ..Combinat: combinations, nrestrictedpartitions

#export factor

InfoPoly2=print
##  <f> must be squarefree.  We test 3 "small" and 2 "big" primes.
function FactorsModPrime(f::Pol{<:Union{Integer,Rational}})
  min=degree(f)+1 # set minimal number of factors to the degree of <f>
  lc=f[end]
    # find a suitable prime
  t=Dict{Symbol, Any}()
  p=1;deg=0;LP=Any[];P=1
  for i in 1:5
        # reset <p> to big prime after first 3 test
    if i==4 p=max(p, 1000) end
    # find a prime not dividing lc(f) and f_p squarefree
    fp=nothing
    while true
      while true
        p=Primes.nextprime(p+1)
        if mod(lc, p)!=0 && mod(f[0], p)!=0 break end
      end
      fp=Pol(FFE{p}.(f.c),f.v)/lc
      if 0==degree(gcd(fp,derivative(fp))) break end
    end
    InfoPoly2("#I  starting factorization mod p:  ", time(), "\n")
    lp=factor(fp)
    sort!(lp,by=degree)
    InfoPoly2("#I  finishing factorization mod p: ", time(), "\n")
        # if <fp> is irreducible so is <f>
    if 1==length(lp)
      InfoPoly2("#I  <f> mod ", p, " is irreducible\n")
      t[:isIrreducible]=true
      return t
    else
      InfoPoly2("#I  found ",length(lp)," factors mod ",p," of degree ",degree.(lp), "\n")
    end
        # choose a maximal prime with minimal number of factors
    if length(lp) <= min
      min=length(lp)
      P=p
      LP=lp
    end
        # compute the possible degrees
    tmp=unique(sum.(combinations(map(degree, lp))))
    if 1==i deg=tmp
    else deg=intersect(deg, tmp)
    end
        # if there is only one possible degree!=0 then <f> is irreducible
    if 2==length(deg)
        InfoPoly2("#I  <f> must be irreducible, only one degree left\n")
        t[:isIrreducible]=true
        return t
    end
  end
    # convert factors <LP> back to the integers
  for i in 1:length(LP)
    LP[i]=Pol(map(x->Int(x),LP[i].c),LP[i].v)
  end
    # return the chosen prime
  InfoPoly2("#I  choosing prime ", P, " with ", length(LP), " factors\n")
  InfoPoly2("#I  possible degrees: ", deg, "\n")
  t[:isIrreducible]=false
  t[:prime]=P
  t[:factors]=LP
  t[:degrees]=deg
  t
end

#  `f` must be square free and must have a non-zero constant term.
function factorSQF(f::Pol{<:Union{Integer,Rational}})
# find a suitable prime, if <f> is irreducible return
  t=FactorsModPrime(f)
  if t[:isIrreducible] return [f] end
  InfoPoly2("#I  using prime ", t[:prime], " for factorization\n")
# for easy combining, we want large degree factors first
  sort!(t[:factors], by=x->-degree(x))
# start Hensel
  h=SquareHensel(f, t)
# combine remaining factors
  fac=[]
# first the factors found by hensel
  if 0<length(h[:remaining])
    InfoPoly2("#I  found ", length(h[:remaining]), " remaining terms\n")
    tmp=TryCombinations(h[:remPolynomial], h[:lc], h[:remaining], h[:primePower], t[:degrees], h[:bounds], true)
    append!(fac, tmp[:irrFactors])
    append!(fac, tmp[:redFactors])
  else
    tmp=Dict{Symbol, Any}()
  end
# append the irreducible ones
  if 0<length(h[:irrFactors])
    InfoPoly2("#I  found ", length(h[:irrFactors]), " irreducibles\n")
    append!(fac, h[:irrFactors])
  end
# and try to factorize the (possible) reducible ones
  if 0<length(h[:redFactors])
    InfoPoly2("#I  found ", length(h[:redFactors]), " reducibles\n")
    if !(haskey(tmp, :stop) || haskey(h, :stop))
  # the stopping criterion has not yet been reached
      for g in h[:redFactors] fac=Append(fac, factorSQF(g)) end
    else fac=Append(fac, h[:redFactors])
    end
  end
  return fac
end

#F  ApproxRational:  approximate r with a <=s-digits denominator
function ApproxRational(r, s)
  n=numerator(r)
  d=denominator(r)
  u=ndigits(d)-s
  if u<=0 return r end
  u=big(10)^u
  div(n,u)//div(d,u)
end

#F  ApproximateRoot(r,e,f=10) . . approximate th e-th root of r
##   with a denominator of 'f' digits.
function ApproximateRoot(r,e=2,f=10)
  RootInt(x,e)=Integer(floor(x^(1/e)))
  x=big(RootInt(numerator(r), e)//RootInt(denominator(r), e))
  nf=r
  c=0
  while true
    lf=nf
    x=ApproxRational(1//e*((e-1)*x+r//x^(e-1)), f+6)
    nf=abs(x^e-r)
    if nf==0 c=6
    else
      if nf>lf lf=nf//lf else lf=lf//nf end
      if lf<2 c+=1 else c=0 end
    end
    if c>2 break end
  end
  x
end

#F  BombieriNorm(<pol>) . . . . . . . . . . . . compute weighted Norm [pol]_2
function BombieriNorm(f)
  n=degree(f)
  f*=big(1)
  ApproximateRoot(sum(i->abs(f[i])^2//binomial(n,i),0:n))
end

#F  MinimizeBombieriNorm(<pol>) . . . . . . . minimize weighted Norm [pol]_2
##                                            by shifting roots
function MinimizeBombieriNorm(f)
  if !(haskey(f, :minimization))
  # this stepwidth should be corrected
    step=1//lcm(denominator.(f.c))
    step=1
    bb="infinity"
    bf=f
    bd=0
    bn=[]
  # evaluation of norm, including storing it (avoids expensive double evals)
    bnf=function(dis) local p, g
      p=filter(i->i[1]==dis, bn)
      if p==[]
          g=Value(f, x+dis)
          p=[dis, BombieriNorm(g)]
          push!(bn, p)
          if bb>p[2]
              bf=g
              bb=p[2]
              bd=dis
          end
          return p[2]
      else return p[1][2]
      end
    end
    x=X(f[:baseRing])
    d=0
    while true
      InfoPoly2("#I Minimizing BombieriNorm, x->x+(", d, ")\n")
  # lokale Parabelann"aherung
      a=bnf(d-step)
      b=bnf(d)
      c=bnf(d+step)
      if a<b && c<b
        if a<c d-=step
        else d+=step
        end
      elseif !(a>b && c>b) && a+c!=2b
        a=-(c-a)//2//(a+c-2b)*step
  # stets aufrunden (wir wollen weg)
        a=step*Int(abs(a)//step+1)*sign(a)
        if a==0 error("sollte nicht")
        else d+=a
        end
      end
# no better can be reached
      if a>b && c>b || all(i->!isempty(filter(j->j[1]==i,bn)), [d-1, d, d+1])
          break
      end
    end
    f[:minimization]=[bf, bd]
  end
  f[:minimization]
end

"""
`factor(f::Pol{<:Union{Integer,Rational}})`

Factor over the integers a polynomial with integral coefficients, or do the
same over the rationals.

"""
function Util.factor(f::Pol{<:Union{Integer,Rational}})
  InfoPoly2("#I  starting integer factorization: ", time(), "\n")
  if iszero(f)
    InfoPoly2("#I  f is zero\n")
    return [f]
  end
  d=lcm(denominator.(f.c))
  f*=d
  f=Pol(Integer.(f.c),f.v)
  v=f.v
  f=shift(f,-f.v)
  if 0==degree(f)
    InfoPoly2("#I  f is a power of x\n")
    s=map(f->Pol(),1:v)
    s[1]*=f.c[1]*lc
    return s
  end
  if 1==degree(f)
    InfoPoly2("#I  f is a linear\n")
    s=map(f->Pol(),1:v)
    push!(s,f)
    return s
  end
  # shift the zeros of f if appropriate
  if degree(f)>20
    g=MinimizeBombieriNorm(f)
    f=g[1]
    shft=-g[2]
  else shft=0
  end
    # make <f> integral, primitive and square free
  g=gcd(f, derivative(f))
  q=exactdiv(f,g)
# @show q
  q=q*sign(q[end])
  InfoPoly2("#I  factorizing polynomial of degree ", degree(q), "\n")
    # and factorize <q>
  if degree(q)<2
    InfoPoly2("#I  <f> is a linear power\n")
    s=[q]
  else
    if q.v>0
      s=[Pol()]
      shift(q,-1)
    else s=[]
    end
    append!(s, factorSQF(q))
  end
  for r in s # find multiple factors
    if 0<degree(g) && degree(g)>=degree(r)
      q=exactdiv(g, r)
      while 0<degree(g) && !isnothing(q)
        push!(s, r)
        g=q
        if degree(g)>=degree(r) q=exactdiv(g,r)
        else q=nothing
        end
      end
    end
  end
    # reshift
  if shft!=0
    InfoPoly2("#I shifting zeros back\n")
    s=map(i->Value(i,Pol()+shft),s)
  end
  append!(s,map(f->Pol(),1:v))
  sort!(s)
  s[1]//=d
  s
end

"""
`LogInt(n,base)`

returns  Int(floor(log(base,n))) but very accurately. Assumes n>0 and base>1.

```julia-repl
julia> Fact.LogInt(1030,2)
10

julia> Fact.LogInt(1,10)
0
```
"""
function LogInt(n, base)
  if n <= 0 error("<n> must be positive") end
  if base <= 1 error("<base> must be greater than 1") end
  function log(b)
    local i
    if b>n return 0 end
    i=2log(b^2)
    if b>n return i
    else n=div(n, b)
        return i+1
    end
  end
  log(base)
end

FFields.Mod(x::Pol,p)=Pol(Mod.(x.c,p),x.v)
FFields.Mod(x::Pol{FFE{p}}) where p=Pol(Mod.(x.c),x.v)

function SquareHensel(f::Pol{<:Union{Integer,Rational}}, t)
local   p,              # prime
      l,              # factorization mod <q>
      k,              # Lift boundary
      prd,            # product of <l>
      rep,            # lifted representation of gcd(<lp>)
      fcn,            # index of true factor in <l>
      dis,            # distance of <f> and <l>
      cor,            # correction
      rcr,            # inverse corrections
      quo,            # quotient
      sum,            # temp
      aa,  bb,        # left and right subproducts
      lq1,            # factors mod <q1>
      max,            # maximum absolute coefficient of <f>
      res,            # result
      gcd,            # used in gcd representation
      i,  j,  x;      # loop
  p=big(t[:prime])
  l=map(x->Pol(Mod.(Integer.(Mod(x).c),p),x.v),t[:factors])
  lc=f[end]
  max=maximum(abs.(f.c))
# compute the factor coefficient bounds
  ofb=2*abs(lc)*OneFactorBound(f)
  InfoPoly2("#I  One factor bound==", ofb, "\n")
  bounds=2*abs(lc)*HenselBound(f)
# compute a representation of the 1 mod <p>
  InfoPoly2("#I  computing gcd representation:",time(),"\n")
  prd=Mod(f,p)/lc
  g=exactdiv(prd,l[1])
  rep=[Pol(Mod(1,p))]
  for i in 2:length(l)
    dis=exactdiv(prd,l[i])
#   xprint("g=",g)
    g,cor1,cor2=gcdx(g,dis)
#   xprintln(" dis=",dis," g=",g," cor1=",cor1," cor2=",cor2)
    rep*=cor1
    push!(rep, cor2)
  end
# for (i,u) in enumerate(rep) xprintln("rep[$i]=",rep[i]) end
  InfoPoly2("#I  representation computed:      ",time(), "\n")
  res=Dict{Symbol,Any}(:irrFactors=>[],:redFactors=>[],
                       :remaining=>[],:bounds=>bounds)
  q=p^2
# start Hensel until <q> is greater than k
  k=bounds[maximum(filter(i->2i<=degree(f), t[:degrees]))]
  InfoPoly2("#I  Hensel bound==", k, "\n")
  q1=p
  while q1<k
    InfoPoly2("#I  computing mod ", q, "\n")
#   for (i,u) in enumerate(l) xprintln("l[$i]=",l[i]) end
    for i in 1:length(l)
      l[i]=Mod(Pol(Integer.(l[i].c),l[i].v),q)
      dis=rem(Mod(f//lc,q),l[i])
      rep[i]=Mod(Pol(Integer.(rep[i].c),rep[i].v),q)
      l[i]+=rem(rep[i]*dis,l[i])
    end
# if this is not the last step update <rep> and check for factors
    if q<k
    # correct the inverses
      for i in 1:length(l)
        if length(l)==1 dis=l[1]^0
        else dis=one(l[1])
          for j in setdiff(1:length(l),[i]) dis=rem(dis*l[j],l[i]) end
        end
        rep[i]=rem(rep[i]*(Mod(2,q)-rep[i]*dis), l[i])
      end
  # try to find true factors
      if max<=q || ofb<q
        InfoPoly2("#I  searching for factors: ", time(), "\n")
        fcn=TryCombinations(f, lc, l, q, t[:degrees], bounds, false)
        InfoPoly2("#I  finishing search:      ", time(), "\n")
      else
        fcn=Dict{Symbol, Any}(:irreducibles => [], :reducibles => [])
      end
  # if we have found a true factor update everything
      if 0<length(fcn[:irreducibles])+length(fcn[:reducibles])
        append!(res[:irrFactors], fcn[:irrFactors])
        append!(res[:redFactors], fcn[:redFactors]) # reducible factors
        prd=prod(fcn[:redFactors];init=one(f))*prod(fcn[:irrFactors];init=one(f))
        f=exactdiv(f,prd)
        if haskey(fcn, :stop)
          res[:stop]=true
          return res
        end
        lc=f[end]
        ofb=2*abs(lc)*OneFactorBound(f)
        InfoPoly2("#I  new one factor bound==", ofb, "\n")
    # degree arguments or OFB arguments prove f irreducible
        if all(i->i==0 || 2i>=degree(f), t[:degrees]) || ofb<q
          push!(fcn[:irrFactors], f)
          push!(res[:irrFactors], f)
          f=f^0
        end
        if degree(f)<1
          InfoPoly2("#I  found non-trivial factorization\n")
          return res
        end
    # compute the factor coefficient bounds
        k=HenselBound(f)
        bounds=map(i->min(bounds[i], k[i]), 1:length(k))
        k=2*abs(lc)*bounds[maximum(filter(i->2i<=degree(f), t[:degrees]))]
        InfoPoly2("#I  new Hensel bound==", k, "\n")
    # remove true factors from <l> and corresponding <rep>
        prd=Mod(prd//prd[end],q)
        l=l[fcn[:remaining]]
        rep=prd.*(rep[fcn[:remaining]])
    # reduce <rep>[i] mod <l>[i]
        for i in 1:length(l) rep[i]=Mod(rem(rep[i], l[i]),q) end
  # if there was a factor, we ought to have found it
      elseif ofb<q
        push!(res[:irrFactors], f)
        InfoPoly2("#I  f irreducible, since one factor would have been found now\n")
        return res
      end
    end
    q1=q
    q=q^2
# avoid a modulus too big
    if q>k q=p^(LogInt(k,p)+1) end
  end
  res[:remPolynomial]=f
  res[:remaining]=l
  res[:primePower]=q1
  res[:lc]=lc
  return res
end

#F  TrialQuotient(<f>,<g>,<b>)  . . . . . . f/g if coeffbounds are given by b
##
function TrialQuotient(f::Pol{T},g,b)where T
# @show f,g,b
  CoeffAbs(p::Pol)=maximum(abs.(numerator.(p.c)))
  CoeffAbs(p::Rational)=abs(numerator(p))
  a=degree(f)-degree(g)
  if a==0
    # Sonderfall abfangen (es mu"s dann eh 0 rauskommen
    a=b[1]
  else
    a=b[a] 
  end
  if iszero(f) return f end
  if f.v<g.v return nothing end
  val=f.v-g.v
  q=[]
  n=length(g.c)
  m=length(f.c)-n
  if T<:Rational
   q=reverse(map(0:m)do i
      c=f[m-i+n-1]/g.c[n]
      if CoeffAbs(c)>a
        InfoPoly2("#I early1 break\n")
        return nothing
      end
      for k in 1:n
        f.c[m-i+k-f.v]-=c*g.c[k]
      end
      c
    end)
  else
    f//=1
    for i in 0:m
      c=exactdiv(f[m-i+n-1],g.c[n])
      if isnothing(c) return c end
      if CoeffAbs(c)>a
        InfoPoly2("#I early2 break\n");
        return nothing
      end
      for k in 1:n
        f.c[m-i+k-f.v]-=c*g.c[k]
#       @show f,m,i,k,b
        if CoeffAbs(f[m-i+k-1])>b
          InfoPoly2("#I early3 break\n")
          return nothing
        end
      end
      q[m-i+1]=c
    end
  end
  for i in 1:m+n
    if !iszero(f[i-1])
       InfoPoly2("#I early4 break\n")
      return nothing
    end
  end
  Pol(q,val)
end

#F  TryCombinations( <f>, ... )  . . . . . . . . . . . . . . . .  try factors
function TryCombinations(f,lc,l,p,alldegs,bounds,split;onlydegs=nothing,stopdegs=nothing)
# xprintln("lc=",lc," p=",p," alldegs=",alldegs," bounds=",bounds," split=",split)
# for u in l xprintln("***=",u) end
  res=Dict{Symbol, Any}(:irreducibles => [], :irrFactors => [], :reducibles => [], :redFactors => [], :remaining => collect(1:length(l)))
  # coefficients should be in [-p/2,p/2]
  deli=map(degree, l)
  sel=collect(1:length(l))
# create List of binomial coefficients to speed up the 'Combinations' process
  binoli=map(0:length(l)-1) do i
    map(j->binomial(i,j),0:i)
  end
  fastbinomial(i,j)=binoli[i+1][j+1]
  step=0
  act=1
  while true
   # factors of larger than half remaining degree we will find as final cofactor
    degf=degree(f)
    degs=filter(i->2i<=degf,alldegs)
    if !isnothing(onlydegs) degs=Intersection(degs, onlydegs) end
    if act in sel
    # search all combinations of length step+1 containing the act-th
    # factor, that are allowed
      good=true
      da=degs.-deli[act]
  # check, whether any combination will be of suitable degree
      cnew=sort(unique(deli[filter(i->i>act, sel)]))
      if any(i->nrestrictedpartitions(i, cnew, step)>0, da)
        InfoPoly2("#I  trying length ", step+1, " containing ", act, "\n")
        cnew=filter(i->i>act,sel)
      else
        InfoPoly2("#I  length ", step+1," containing ",act, " not feasible\n")
        cnew=[]
      end
      mind=sum(deli) # the maximum of the possible degrees. We surely will find something smaller
      lco=binomial(length(cnew), step)
      if 0==lco mind=0
      else
        InfoPoly2("#I  ",lco," combinations\n")
        i=1
        while good && i<=lco
          q=i
          d=length(cnew)
          o=0
          combi=[]
          for ii in step-1:-1:0
            j=1
            b=fastbinomial(d-1,ii)
            while q>b
              q-=b
              b*=(d-j-ii)//(d-j) # b==binomial(d-(j+1),ii)
              j+=1
            end
            o+=j
            d-=j
            push!(combi, cnew[o]);sort!(combi)
          end
  # check whether this yields a minimal degree
          d=sum(deli[combi])
          if d<mind mind=d end
          if d in da
            if !(act in combi) push!(combi, act);sort!(combi) end
  # make sure that the quotient has a chance, compute the
  # extremal coefficient of the product:
            q=Integer(prod(map(i->i.c[1],l[combi]))*lc)
            if length(combi)==2
#           @show q,l[combi]
            end
  # As  we don't know yet the gcd of all the products coefficients (to make
  # it primitive), we do a slightly weaker test: (test of leading coeffs is
  # first  in  'TrialQuotient')  this  just  should  reduce  the  number of
  # 'ProductMod'  neccessary. the absolute part  of the product must divide
  # the absolute part of f up to a divisor of <lc>
            q=f.c[1]//q*lc
            if length(combi)==2
#           @show q
            end
            prd=zero(l[1])
            if !isinteger(q)
              InfoPoly2("#I  ignoring combination ", combi, "\n")
              q=nothing
            else
              InfoPoly2("#I  testing combination ", combi, "\n")
              prd=prod(l[combi])
              cof=Integer.(prd.c*lc)
  # make the product primitive
              cof*=1//gcd(cof)
              prd=Pol(cof,prd.v)
              q=TrialQuotient(f//1, prd, bounds)
            end
            if !isnothing(q)
              f=q
              InfoPoly2("#I  found true factor of degree ", degree(prd), "\n")
              if length(combi)==1 || split q=0
              else
                q=2*lc*OneFactorBound(prd)
                if q <= p
                  InfoPoly2("#I  proven irreducible by ", "'OneFactorBound'\n")
                end
              end
              if q <= p
                append!(res[:irreducibles], combi)
                push!(res[:irrFactors], prd)
                if !isnothing(stopdegs) && degree(prd) in stopdegs
                  InfoPoly2("#I  hit stopdegree\n")
                  push!(res[:redFactors], f)
                  res[:stop]=true
                  return res
                end
              else
                push!(res[:reducibles], combi)
                push!(res[:redFactors], prd)
              end
              setdiff!(res[:remaining], combi)
              good=false
              setdiff!(sel, combi)
            end
          end
          i+=1
        end
      end
  # we can forget about the actual factor, as any longer combination
  # is too big
      if length(degs)>1 && deli[act]+mind>=maximum(degs)
        InfoPoly2("#I  factor ", act, " can be further neglected\n")
        setdiff!(sel, [act])
      end
    end
    act+=1
    if 0<length(sel) && act>maximum(sel)
      step+=1
      act=sel[1]
    end
    if 0==length(sel) || length(sel)<step break end
  end
  if split && (0<length(res[:remaining]) && f!=f^0)
    append!(res[:irreducibles], res[:remaining])
    res[:remaining]=[]
    push!(res[:irrFactors], f)
  end
  res
end

function OneFactorBound(f)
  n=degree(f)
  if n>=3
    # Single factor bound of Beauzamy, Trevisan and Wang (1993)
      return numerator(floor(10912//10000*(ApproximateRoot(2^n,2)//
             ApproximateRoot(n^3,8)*ApproximateRoot(BombieriNorm(f),2))))+1
  else
    # Mignotte's single factor bound
    d=div(n,2)
    binomial(d,div(d,2))*(1+Integer(floor(sqrt(sum(i->i^2,f.c)))))
  end
end

#F  Beauzamy's Bound for Factors Coefficients  cf. JSC 13 (1992), 463-472
function BeauzamyBound(f)
  n=big(degree(f))
  1+Integer(floor(
  # the strange number in the next line is an (upper) rational approximation
  # for 3^{3/4}/2/\sqrt(\pi)
  643038/1000000*ApproximateRoot(3^n,2)/ApproximateRoot(n,2)*BombieriNorm(f)))
end

#############################################################################
##
#F  ApproxRootBound(f) Numerical approximation of RootBound (better, but
##  may fail)
function ApproxRootBound(f::Pol{<:Rational})
  x=Pol()*one(f.c[1])
  f=shift(f,-f.v)
  # probably first test, whether polynomial should be inverted. However,
  # we expect roots larger than one.
  d=degree(f)
  f=Pol(reverse(f.c))
  app=1//2
  diff=1//4
  nkon=true
  while true
    # pol, whose roots are the 1/app of the roots of f
    tp=f(x*app)
    tp=tp.c
    tp/=tp[1]
    tp=map(i->ApproxRational(i,10),tp)
    # now check, by using the Lehmer/Schur method, whether tp has a root
    # in the unit circle, i.e. f has a root in the app-circle
    v=0
    while true
      fail=false
      p=tp
      while true
        d=length(p)
        # compute T[p]=\bar a_n p-a_0 p*, everything rational.
        pl=p
        p=p[1]*p-p[d]*reverse(p)
        p=map(i->ApproxRational(i,10),p)
        d=findlast(!iszero,p[2:end])
        if isnothing(d) d=1 else d+=1 end
        p=p[1:d]
        v=p[1]
        if v==0 fail=nkon end
        nkon=any(!zero,p[2:end])
        if v<=0 break end
      end
      if fail
        # we fail due to rounding errors
        return nothing
      else
        if v<0 # zero in the unit circle, app smaller
          app-=diff
        else # no zero in the unit circle, app larger
          app+=diff
        end
      end
      if !fail break end
    end
    diff/=2
    # until good circle found, which does not contain roots.
    if v==0 && (1-app/(app+diff))<1/40 break end
  end

  # revert last enlargement
  app-=2*diff
  1/app+1/20
end

#F  RootBound(f) . . . . bound for absolute value of (complex) roots of f
function RootBound(f)
  # valuation gives only 0 as zero, this can be neglected
  if f.v>0 f=shift(f,-f.v) end
  # normieren
  f=f//f[end]
  a=ApproxRootBound(f)
  # did the numerical part fail?
  if isnothing(a)
    b=f.baseRing
    # we use, that AbsInt works for also for rationals!
    b=Pol([2],degree(f))-Pol(abs.(f.c))
    a=ApproxRootBound(b)
    if isnothing(a)
      c=abs.(f.c)
      d=length(c)
      a=max(1,sum(c[1:d-1]))
      b=1+maximum(c)
      if b<a a=b end
      b=maximum(map(i->RootInt(d*Int(abs(c[d-i])+1//2),i)+1,1:d-1)) 
      if b<a a=b end
      if all(!iszero,c)
        b=map(i->2*ans(c[i-1]/c[i]),3:d)
        push!(b,ans(c[1]/c[2]))
        b=maximum(b)
        if b<a a=b end
      end
      b=sum(i->abs(c[i]-c[i+1]),1:d-1)+abs(c[1])
      if b<a a=b end
      b=maximum(map(i->RootInt(Int(abs(c[d-i]/binomial(d-1,i))+1/2),i)+1,
                     1:d-1))/(ApproximateRoot(2,d-1)-1)+10^-10
      if b<a a=b end
    end
  end
  a
end

#F  HenselBound(<pol>,[<minpol>,<den>]) . . . Bounds for Factor coefficients
##    if the computation takes place over an algebraic extension, then
##    minpol and denominator must be given
function HenselBound(pol,arg...)
  if length(arg)>0
    n,d=arg
    dis=discriminant(n)
    nalpha=RootBound(n)
    if !(all(IsRat, pol[:coefficients]))
    # now try to bound the roots of f accordingly. As in all estimates by
    # RootBound only the absolute value of the coefficients is used, we will
    # estimate these first, and replace f by the polynomial
    # x^n-b_{n-1}x^(n-1)-...-b_0 whose roots are certainly larger
      a=[]
      for i in pol[:coefficients]
        if IsRat(i) push!(a, abs(i))
        else push!(a, Sum(i[:coefficients], abs)*nalpha)
        end
      end
      a=-a
      a[length(a)]=-(a[length(a)])
      pol=Polynomial(Rationals, a)
    else
      pol=Polynomial(Rationals, pol[:coefficients])
    end
    n=degree(n)
  else
    n=1
  end
  rb=0
  bea=BeauzamyBound(pol)
  # compute Landau-Mignotte bound for absolute values of
  # coefficients of any factor
  w=sum(i->i^2,pol.c)
  # we want an upper bound of the root
  # As we nowhere selected a specific galois representative,
  # this bound (which is rational!) will bound all conjugactes as well.
  lm=Integer(ceil(sqrt(w)))
  lb=2^div(degree(pol), 2)*lm
  bound=map(1:degree(pol))do k
    l=2^k*lm
    if l<bea w=l
    else w=bea
    end
    if bea>big(10)^200 || n>1
      if rb==0
        rb=RootBound(pol)
        a=rb
      end
      bin=1
      for j in 1:k
        bin*=(k-j+1)//j
        w=bin*rb^j
        if w>a a=w end
      end
      if a<l w=a else w=l end
    end
    if n>1
    # algebraic Extension case
    # finally we have to bound (again) the coefficients of \alpha when
    # writing the coefficients of the factor as \sum c_i/d\alpha^i.
      w=Int(d*w*factorial(n)//RootRat(abs(dis))*nalpha^((n*(n-1))//2))+1
    end
    Integer(floor(w))+1
  end
  bound
end

q=Pol(:q)
testf=
q^17 + (28475495154853885209666581640806568193580138878464//35)*q^13 + (
-230113592530021909980533244143226208753492718798402862328238781356987464425343796365306003458195456//1225)*q^9 + (
-1331663813192583136877354316616861253307245376897389186616397648441393194005806591992969924972327944810995692943683480696873183913411848253822992384//
42875)*q^5 + (
-11476141319715558689839068201773325238493770231736218788232819936653287068633455026978626687107040865767792194345449338617565687175571484254603185906485112033995725244553550967874479529656320//343)*q
end
