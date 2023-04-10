"""
Factoring polynomials over the rationals
"""
module Fact
import Primes: Primes, nextprime, factor
using LinearAlgebra:exactdiv
using LaurentPolynomials: Pol, @Pol, shift, degree, derivative, valuation, coefficients
using ..FFields: FFields, FFE
using ..Modulo: Modulo, Mod
using ..Combinat: Combinat, combinations, npartitions
using ..Tools: improve_type

const verbose=Ref(false)
function info(x...)
 if verbose[] println(x...) end
end

const ltime=Ref(0.0)

function etime()
  newt=time()
  d=newt-ltime[]
  ltime[]=newt
  string(round(d;digits=5),"s")
end

##  f must be squarefree and f[0]!=0.  We test 3 "small" and 2 "big" primes.
function FactorsModPrime(f::Pol{<:Union{Integer,Rational}})
  min=degree(f)+1 # set minimal number of factors to the degree of <f>
  lc=f[end]
    # find a suitable prime
  t=Dict{Symbol, Any}()
  p=1;deg=0;LP=Any[];P=1
  for i in 1:5
    if i==4 p=max(p, 1000) end # reset <p> to big prime after first 3 test
    # find a prime not dividing lc(f) and f_p squarefree
    fp=nothing
    while true
      while true
        p=Primes.nextprime(p+1)
        if mod(lc, p)!=0 && mod(f[0], p)!=0 break end
      end
      fp=Pol{FFE{p}}(f)/lc
      if 0==degree(gcd(fp,derivative(fp))) break end
    end
    info("#I  starting factorization mod p:  ", etime())
    lp=factor(fp)
    sort!(lp,by=degree)
    info("#I  finishing factorization mod p: ", etime())
        # if <fp> is irreducible so is <f>
    if 1==length(lp)
      info("#I  <f> mod ", p, " is irreducible")
      t[:isIrreducible]=true
      return t
    else
      info("#I  found ",length(lp)," factors mod ",p," of degree ",degree.(lp))
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
        info("#I  <f> must be irreducible, only one degree left")
        t[:isIrreducible]=true
        return t
    end
  end
  LP=map(Pol{Int},LP) # convert factors <LP> back to the integers
    # return the chosen prime
  info("#I  choosing prime ", P, " with ", length(LP), " factors")
  info("#I  possible degrees: ", deg)
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
  info("#I  using prime ", t[:prime], " for factorization")
# for easy combining, we want large degree factors first
  sort!(t[:factors], by=x->-degree(x))
# start Hensel
  h=SquareHensel(f, t)
# combine remaining factors
  fac=[]
# first the factors found by hensel
  if 0<length(h[:remaining])
    info("#I  found ", length(h[:remaining]), " remaining terms")
    tmp=TryCombinations(h[:remPolynomial], h[:lc], h[:remaining], h[:primePower], t[:degrees], h[:bounds], true)
    append!(fac, tmp[:irrFactors])
    append!(fac, tmp[:redFactors])
  else
    tmp=Dict{Symbol, Any}()
  end
# append the irreducible ones
  if 0<length(h[:irrFactors])
    info("#I  found ", length(h[:irrFactors]), " irreducibles")
    append!(fac, h[:irrFactors])
  end
# and try to factorize the (possible) reducible ones
  if 0<length(h[:redFactors])
    info("#I  found ", length(h[:redFactors]), " reducibles")
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

#F  BombieriNorm(pol) . . . . . . . . . . . . compute weighted Norm [pol]_2
function BombieriNorm(f)
  n=degree(f)
  f*=big(1)
  ApproximateRoot(sum(i->abs(f[i])^2//binomial(n,i),0:degree(f)))
end

#F  MinimizeBombieriNorm(pol) . . . . . . . minimize weighted Norm [pol]_2
##                                            by shifting roots
function MinimizeBombieriNorm(f::Pol)
  # this stepwidth should be corrected
# step=1//denominator(f)
  step=1
  mn=(bb=big(-1),pol=f,dis=0)
  bn=Tuple{Int,Rational{BigInt}}[]
# evaluation of norm, including storing it (avoids expensive double evals)
  bnf=function(dis)
    p=filter(i->first(i)==dis, bn)
    if !isempty(p) return last(p[1]) end
    g=f(Pol()+dis) # assumes f isa Pol
    q=(dis, BombieriNorm(g))
    push!(bn, q)
    if mn.bb==-1 || mn.bb>last(q)
      mn=(bb=last(q),pol=g,dis=dis)
    end
    q[2]
  end
  d=0//1
  while true
    info("#I Minimizing BombieriNorm, x->x+(", d, ")")
# approximation of local parabola
    a=Rational{BigInt}(bnf(d-step))
    b=Rational{BigInt}(bnf(d))
    c=Rational{BigInt}(bnf(d+step))
    if a<b && c<b
      if a<c d-=step
      else d+=step
      end
    elseif !(a>b && c>b) && a+c!=2b
      a=-(c-a)//2//(a+c-2b)*step
# always round (we want that)
      a=step*(div(abs(a),step)+1)*sign(a)
      if a==0 error("theory")
      else d+=a
      end
    end
# no better can be reached
    if a>b && c>b || all(i->!iszero(count(j->first(j)==i,bn)),[d-1,d,d+1])
        break
    end
  end
  (mn.pol,-mn.dis)
end

"""
`factor(f::Pol{<:Union{Integer,Rational}})`

Factor over the integers a polynomial with integral coefficients, or do the
same over the rationals.

```julia-repl
julia> factor(Pol(:q)^24-1)
8-element Vector{Pol{Int64}}:
 q-1
 q²-q+1
 q⁴-q²+1
 q⁸-q⁴+1
 q⁴+1
 q²+1
 q+1
 q²+q+1
```
"""
function Primes.factor(f::Pol{<:Union{Integer,Rational}})
  info("#I  starting integer factorization: ", etime())
  if iszero(f)
    info("#I  f is zero")
    return [f]
  end
  d=denominator(f)
  f*=d
  f=Pol{typeof(d)}(f)
  v=valuation(f)
  f=shift(f,-v)
  if 0==degree(f)
    info("#I  f is a power of x")
    s=map(f->Pol(),1:v)
    s[1]*=f[0]
    return s
  end
  if 1==degree(f)
    info("#I  f is linear")
    s=map(f->Pol(),1:v)
    push!(s,f)
    return s
  end
  # shift the zeros of f if appropriate
  if degree(f)>20 f,shft=MinimizeBombieriNorm(f)
  else shft=0
  end
    # make <f> integral, primitive and square free
  g=gcd(f, derivative(f))
  q=exactdiv(f,g)
# @show q
  q=q*sign(q[end])
  info("#I  factorizing polynomial of degree ", degree(q))
    # and factorize <q>
  if degree(q)<2
    info("#I  <f> is a linear power")
    s=[q]
  else
    if valuation(q)>0
      s=[Pol()]
      shift(q,-1)
    else s=empty([q])
    end
    append!(s, factorSQF(q))
  end
  for r in s # find multiple factors
    if 0<degree(g) && degree(g)>=degree(r)
      q=exactdiv(g, r)
      while 0<degree(g) && !isnothing(q)
        push!(s, r)
        g=q
        if degree(g)>=degree(r) 
          if iszero(rem(g//1,r)) q=exactdiv(g,r)
          else q=nothing
          end
        else q=nothing
        end
      end
    end
  end
    # reshift
  if shft!=0
    info("#I shifting zeros back")
    s=map(i->Value(i,Pol()+shft),s)
  end
  append!(s,map(f->Pol(),1:v))
  sort!(s)
  s[1]//=d
  improve_type(s)
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
    if b>n return 0 end
    i=2log(b^2)
    if b>n return i
    else n=div(n, b)
        return i+1
    end
  end
  log(base)
end

Modulo.Mod(x::Pol,p)=Pol(Mod.(x.c,p),x.v)
Modulo.Mod(x::Pol{FFE{p}}) where p=Pol(Mod.(x.c),x.v)
Modulo.Mod(x::Mod,p)=Mod(x.val,p)

function SquareHensel(f::Pol{<:Union{Integer,Rational}}, t)
#   p,              # prime
#   l,              # factorization mod <q>
#   k,              # Lift boundary
#   prd,            # product of <l>
#   rep,            # lifted representation of gcd(<lp>)
#   fcn,            # index of true factor in <l>
#   dis,            # distance of <f> and <l>
#   cor,            # correction
#   rcr,            # inverse corrections
#   quo,            # quotient
#   aa,  bb,        # left and right subproducts
#   lq1,            # factors mod <q1>
#   max,            # maximum absolute coefficient of <f>
#   gcd,            # used in gcd representation
  p=big(t[:prime])
  l=Mod.(t[:factors],p)
  lc=f[end]
  max=maximum(abs.(coefficients(f)))
# compute the factor coefficient bounds
  ofb=2*abs(lc)*OneFactorBound(f)
  info("#I  One factor bound==", ofb)
  bounds=2*abs(lc)*HenselBound(f)
# compute a representation of the 1 mod <p>
  info("#I  computing gcd representation:",etime())
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
# for (i,u) in pairs(rep) xprintln("rep[$i]=",rep[i]) end
  info("#I  representation computed:      ",etime())
  res=Dict{Symbol,Any}(:irrFactors=>[],:redFactors=>[],
                       :remaining=>[],:bounds=>bounds)
  q=p^2
# start Hensel until <q> is greater than k
  k=bounds[maximum(filter(i->2i<=degree(f), t[:degrees]))]
  info("#I  Hensel bound==", k)
  q1=p
  while q1<k
    info("#I  computing mod ", q)
#   for (i,u) in pairs(l) xprintln("l[$i]=",l[i]) end
    for i in 1:length(l)
      l[i]=Mod(l[i],q)
      dis=rem(Mod(f//lc,q),l[i])
      rep[i]=Mod(rep[i],q)
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
        info("#I  searching for factors: ", etime())
        fcn=TryCombinations(f, lc, l, q, t[:degrees], bounds, false)
        info("#I  finishing search:      ", etime())
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
        info("#I  new one factor bound==", ofb)
    # degree arguments or OFB arguments prove f irreducible
        if all(i->i==0 || 2i>=degree(f), t[:degrees]) || ofb<q
          push!(fcn[:irrFactors], f)
          push!(res[:irrFactors], f)
          f=f^0
        end
        if degree(f)<1
          info("#I  found non-trivial factorization")
          return res
        end
    # compute the factor coefficient bounds
        k=HenselBound(f)
        bounds=map(i->min(bounds[i], k[i]), 1:length(k))
        k=2*abs(lc)*bounds[maximum(filter(i->2i<=degree(f), t[:degrees]))]
        info("#I  new Hensel bound==", k)
    # remove true factors from <l> and corresponding <rep>
        prd=Mod(prd//prd[end],q)
        l=l[fcn[:remaining]]
        rep=prd.*(rep[fcn[:remaining]])
    # reduce <rep>[i] mod <l>[i]
        for i in 1:length(l) rep[i]=Mod(rem(rep[i], l[i]),q) end
  # if there was a factor, we ought to have found it
      elseif ofb<q
        push!(res[:irrFactors], f)
        info("#I  f irreducible, since one factor would have been found now")
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
      c=f[m-i+n-1]/g[end]
      if CoeffAbs(c)>a
        info("#I early1 break")
        return nothing
      end
      f.c[(m-i-f.v).+(1:n)].-=c*g.c
      c
    end)
  else
    f//=1
    for i in 0:m
      c=exactdiv(f[m-i+n-1],g[end])
      if isnothing(c) return c end
      if CoeffAbs(c)>a
        info("#I early2 break");
        return nothing
      end
      for k in 1:n
        f.c[m-i+k-f.v]-=c*g.c[k]
#       @show f,m,i,k,b
        if CoeffAbs(f[m-i+k-1])>b
          info("#I early3 break")
          return nothing
        end
      end
      q[m-i+1]=c
    end
  end
  for i in 1:m+n
    if !iszero(f[i-1])
       info("#I early4 break")
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
      if any(i->npartitions(i, cnew, step)>0, da)
        info("#I  trying length ", step+1, " containing ", act)
        cnew=filter(i->i>act,sel)
      else
        info("#I  length ", step+1," containing ",act, " not feasible")
        cnew=[]
      end
      mind=sum(deli) # the maximum of the possible degrees. We surely will find something smaller
      lco=binomial(length(cnew), step)
      if 0==lco mind=0
      else
        info("#I  ",lco," combinations")
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
            q=Integer(prod(map(i->i[0],l[combi]))*lc)
            if length(combi)==2
#           @show q,l[combi]
            end
  # As  we don't know yet the gcd of all the products coefficients (to make
  # it primitive), we do a slightly weaker test: (test of leading coeffs is
  # first  in  'TrialQuotient')  this  just  should  reduce  the  number of
  # 'ProductMod'  neccessary. the absolute part  of the product must divide
  # the absolute part of f up to a divisor of <lc>
            q=f[0]//q*lc
            if length(combi)==2
#           @show q
            end
            prd=zero(l[1])
            if !isinteger(q)
              info("#I  ignoring combination ", combi)
              q=nothing
            else
              info("#I  testing combination ", combi)
              prd=prod(l[combi])
              cof=Integer.(prd.c*lc)
  # make the product primitive
              cof*=1//gcd(cof)
              prd=Pol(cof,prd.v)
              q=TrialQuotient(f//1, prd, bounds)
            end
            if !isnothing(q)
              f=q
              info("#I  found true factor of degree ", degree(prd))
              if length(combi)==1 || split q=0
              else
                q=2*lc*OneFactorBound(prd)
                if q <= p
                  info("#I  proven irreducible by ", "'OneFactorBound'")
                end
              end
              if q <= p
                append!(res[:irreducibles], combi)
                push!(res[:irrFactors], prd)
                if !isnothing(stopdegs) && degree(prd) in stopdegs
                  info("#I  hit stopdegree")
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
        info("#I  factor ", act, " can be further neglected")
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
  x=Pol()*one(f[0])
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
p=-1434517664964444836229883525221665654811721278967027348529102492081660883579181878372328335888380108220974024293181167327195710896946435531825398238310639004249465655569193870984309941207040000*q-1331663813192583136877354316616861253307245376897389186616397648441393194005806591992969924972327944810995692943683480696873183913411848253822992384*q^5-8053975738550766849318663545012917306372245157944100181488357347494561254887032872785710121036840960*q^9+34882481564696009381841562509988046037135670126118400*q^13+42875*q^17
end
