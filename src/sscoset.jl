##  This file contains functions dealing with semisimple elements of algebraic
##  cosets.
##
## An algebraic coset G.σ quasi-central is represented
## by a CoxeterCoset WF, where WF.F0Mat is the action of σ on X(T).
##
## A finite order quasi-semisimple element tσ is represented as t, an 
## element of Y(T)\otimes Q/Z=(Q/Z)^n, that is a list of length n of rationals r
## such that 0<=r<1.
"""
\Section{Quasi-Semisimple elements of non-connected reductive groups}

We  may  also  use  Coxeter  cosets  to represented non-connected reductive
groups  of the form $\bG\rtimesσ$ where $\bG$ is a connected reductive
group   and  $σ$  an   algebraic  automorphism  of   $\bG$,  and  more
specifically the coset $\bG.σ$. We may always choose
$σ\in\bG⋅σ$  *quasi-semisimple*,  which  means  that $σ$
preserves a pair $\bT\subset\bB$ of a maximal torus and a Borel subgroup of
$\bG$.  If $σ$  is of  finite order,  it then  defines an automorphism
$F_0$ of the root datum $(X(\bT), \Phi, Y(\bT), \Phi^\vee)$, thus a Coxeter
coset. We refer to \cite{ss} for details.

We  have  extended  the  functions  for  semi-simple  elements to work with
quasi-semisimple   elements   $tσ\in\bT⋅σ$.   Here,   as   in
\cite{ss},  $σ$ is a quasi-central  automorphism uniquely defined by a
diagram  automorphism  of  $(W,S)$,  taking  $σ$  symplectic  in  type
$A_{2n}$.  We  recall  that  a  quasi-central element is a quasi-semisimple
element such that the Weyl group of $C_\bG(σ)$ is equal to $W^σ$;
such an element always exists in the coset $\bG⋅σ$.

Here are some examples\:

|    gap> WF:=RootDatum("u",6);
    2A5.(q+1)|

The above defines the coset $\GL_6⋅σ$ where $σ$ is the composed
of transpose, inverse and the longest element of $W$.

|    gap> l:=QuasiIsolatedRepresentatives(WF);
    [ <0,0,0,0,0,0>, <1/4,0,0,0,0,3/4>, <1/4,1/4,0,0,3/4,3/4>,
      <1/4,1/4,1/4,3/4,3/4,3/4> ]|

we define an element $tσ\in\bT⋅σ$ to be quasi-isolated if the
Weyl  group of $C_\bG(tσ)$ is not  in any proper parabolic subgroup of
$W^σ$. This generalizes the definition for connected groups. The above
shows  the  elements  $t$  where  $tσ$  runs  over  representatives of
quasi-isolated  quasi-semisimple  classes  of  $\bG⋅σ$.  The given
representatives have been chosen $σ$-stable.

|    gap> List(l,s->Centralizer(WF,s));
    [ C3<3,2,1>, B2.(q+1), (A1xA1)<1,3>xA1<2>, 2A3<3,1,2> ]|

in  the above,  the groups  $C_\bG(tσ)$ are  computed and displayed as
extended  Coxeter groups (following the same convention as for centralisers
in connected reductive groups).

We  define an element $tσ\in\bT⋅σ$ to be isolated if the Weyl
group  of $C_\bG(tσ)^0$  is not  in any  proper parabolic  subgroup of
$W^σ$. This generalizes the definition for connected groups.

|    gap> List(l,s->IsIsolated(WF,s));
    [ true, false, true, true ]|
"""


# IsSpecial(WF,c) c is an orbit of WF.phi on WF.rootInclusion
# return true iff c is special in the sense of Digne-Michel
function IsSpecial(WF,c)
  if mod(length(c),2)!=0 return false end
  W=Group(WF)
  roots(W,restriction(W,c[1]))+roots(W,restriction(W,c[1+length(c)//2])) in roots(W)
end

# rootdatum R(σ)
# Computes X_σ X^σ, Y_σ, Y^σ, R(σ) of ss
# see 1.1 to 1.7 of \cite{ss}
function RelativeDatum(WF)
  gets(WF,:Rs)do
    W=Group(WF)
    n=OrderMat(WF[:F0Mat]) # matrix of sigma on X
    WF[:pi] = sum(i->WF[:F0Mat]^i,0:n-1)//n
    WF[:X_s] = BaseIntMat(n * WF[:pi]) // n # basis of X_σ
    WF[:Y_s] = BaseIntMat(n * TransposedMat(WF[:pi])) // n # basis of Y_σ
    WF[:Xs] = (WF[:X_s]*TransposedMat(WF[:Y_s]))^ -1 * WF[:X_s] # basis of X^σ
    WF[:Ys] = (WF[:Y_s]*TransposedMat(WF[:X_s]))^ -1 * WF[:Y_s] # basis of Y^σ
    cc = orbits(WF[:phi], inclusiongens(W))
    Phis = map(function(c)
      res = Sum(roots(W,restriction(W,c))) * W[:simpleRoots]
      if IsSpecial(WF, c) res = 2res end
      res
    end, cc)
    cPhis=map(c->Sum(coroots(W,restriction(W,c))*W[:simpleCoroots])
                //length(c), cc)
    CoxeterGroup(map(x->SolutionMat(WF[:Xs], x), Phis), 
                 map(x->SolutionMat(WF[:Y_s], x), cPhis))
  end
end

function Cso(WF)# compute constants C_σ,α
  gets(WF, Cso)do
    W=Group(WF)
    res=fill(1, max(0,W[:N]))
    for o in refltype(WF)
      if order(o.twist)==2 && o.orbit[1].series==:A && mod(o.orbit[1].rank,2)==0
        for p in o.orbit
          res[p[:indices][p[:rank]//2+[0, 1]]]=[-1,-1]
        end
      end
    end
    function C(i)
      for j in 1:semisimplerank(W)
        if W.rootdec[i][j]>0
          p=Position(W[:roots], W[:roots][i]-W[:roots][j])
          if p!=false return res[j]*res[p] end
        end
      end
    end
    for i in semisimplerank(W).+(1:nref(W)) res[i]=C(i) end
    append!(res,res)
  end
end

# computes the centralizer of tσ as an ExtendedReflectionGroup
"""
\Section{Centralizer for quasisemisimple elements}

'Centralizer(<WF>, <t>)'

<WF>   should  be   a  Coxeter   coset  representing   an  algebraic  coset
$\bG⋅σ$,  where $\bG$ is a  connected reductive group (represented
by  'W:=Group(WF)'), and $σ$ is a quasi-central automorphism of $\bG$
defined  by <WF>. The element <t> should  be a semisimple element of $\bG$.
The    function   returns   an   extended   reflection   group   describing
$C_\bG(tσ)$,    with   the   reflection    group   part   representing
$C_\bG^0(tσ)$,  and the diagram automorphism  part being those induced
by $C_\bG(tσ)/C_\bG(tσ)^0$ on $C_\bG(tσ)^0$.

|    gap> WF:=RootDatum("u",6);
    2A5.(q+1)
    gap> s:=SemisimpleElement(Group(WF),[1/4,0,0,0,0,3/4]);
    <1/4,0,0,0,0,3/4>
    gap> Centralizer(WF,s);
    B2.(q+1)
    gap> Centralizer(WF,s^0);
    C3<3,2,1>|
"""
function centralizer(WF::Spets,t::SemisimpleElement)
  W=Group(WF)
  RelativeDatum(WF)
  refC = Centralizer(WF[:Rs], SemisimpleElement(WF[:Rs], SolutionMat(WF[:Y_s], t[:v] * TransposedMat(WF[:pi]))))
  Rs = map(function (c,)
              local res
              res = [c, Sum((W[:roots])[(W[:rootRestriction])[c]] * W[:simpleRoots]) // length(c), Sum((W[:coroots])[(W[:rootRestriction])[c]]) * W[:simpleCoroots]]
              if IsSpecial(WF, c)
                  res[3] = 2 * res[3]
              end
              return res
          end, orbits(WF[:phi], (W[:rootInclusion])[1:W[:N]]))
  labels = map(joindigits, orbits(WF[:phi], (W[:rootInclusion])[1:W[:N]]))
  if t[:additive]
      good = map((p->begin
                      Mod1(Sum(p[1], (i->begin t ^ ((parent(W))[:roots])[i]
                                      end)) + AsRootOfUnity(Cso(WF)[(p[1])[1]])) == 0
                  end), Rs)
  else
    good=map((p->begin
                 Product(p[1], (i->begin
                                    t ^ ((parent(W))[:roots])[i]
                                   end)) * Cso(WF)[(p[1])[1]] == (t[:v])[1] ^ 0
                end), Rs)
  end
  Rs = ListBlist(Rs, good)
  labels = ListBlist(labels, good)
  cRs = map((x->begin x[3] end), Rs)
  cRs = map((x->begin SolutionMat(WF[:Ys], x) end), cRs)
  cRs = Filtered(cRs, (x->begin !x in map(Sum, cartesian(cRs, cRs)) end))
  Rs = map((x->begin x[2] end), Rs)
  Rs = map((x->begin SolutionMat(WF[:X_s], x) end), Rs)
  good = map((x->begin !x in map(Sum, cartesian(Rs, Rs)) end), Rs)
  Rs = ListBlist(Rs, good)
  labels = ListBlist(labels, good)
  if length(Rs) > 0 C = CoxeterGroup(Rs, cRs)
  else C = torus(length(WF[:Xs]))
  end
  (C[:operations])[:ReflectionFromName] = function (W, x)
          return Position(W[:rootInclusion], x)
      end
  p = map((x->begin SolutionMat(WF[:X_s], x) end), WF[:Xs])
# transfer matrix on X^σ to X_σ
  if refC[:F0s] == [] return ExtendedReflectionGroup(C) end
  return ExtendedReflectionGroup(C, ApplyFunc(Group, map((x->begin
                          x^p end), refC[:F0s])))
end

# returns representatives of quasi-isolated classes of G.σ
"""
\Section{QuasiIsolatedRepresentatives for Coxeter cosets}

'QuasiIsolatedRepresentatives(<WF>[, <p>])'

<WF>   should  be   a  Coxeter   coset  representing   an  algebraic  coset
$\bG⋅σ$,  where $\bG$ is a  connected reductive group (represented
by  'W:=Group(WF)'), and $σ$ is  a quasi-central automorphism of $\bG$
defined  by <WF>.  The function  returns a  list of  semisimple elements of
$\bG$   such  that   $tσ$,  when   $t$  runs   over  this   list,  are
representatives  of the conjugacy classes of quasi-isolated quasisemisimple
elements  of  $\bG⋅σ$  (an  element  $tσ\in\bT⋅σ$ is
quasi-isolated  if the Weyl group of  $C_\bG(tσ)$ is not in any proper
parabolic  subgroup of $W^σ$).  If a second  argument <p> is given, it
lists only those representatives which exist in characteristic <p>.

|    gap> QuasiIsolatedRepresentatives(RootDatum("2E6sc"));
    [ <0,0,0,0,0,0>, <0,0,0,1/2,0,0>, <0,1/2,1/4,0,1/4,0>,
      <0,2/3,0,1/3,0,0>, <0,3/4,0,1/2,0,0> ]
    gap> QuasiIsolatedRepresentatives(RootDatum("2E6sc"),2);
    [ <0,0,0,0,0,0>, <0,2/3,0,1/3,0,0> ]
    gap> QuasiIsolatedRepresentatives(RootDatum("2E6sc"),3);
    [ <0,0,0,0,0,0>, <0,0,0,1/2,0,0>, <0,1/2,1/4,0,1/4,0>,
      <0,3/4,0,1/2,0,0> ]|
"""
CoxeterCosetOps[:QuasiIsolatedRepresentatives] = function (arg...,)
  WF = arg[1]
  if length(arg) == 2 p = arg[2]
  else p = 0
  end
  return map((x->begin
                  SemisimpleElement(Group(WF), x[:v] * WF[:Y_s])
              end), QuasiIsolatedRepresentatives(RelativeDatum(WF), p))
end

# whether tσ is isolated
"""
\Section{IsIsolated for Coxeter cosets}

'IsIsolated(<WF>, <t>)'

<WF>   should  be   a  Coxeter   coset  representing   an  algebraic  coset
$\bG⋅σ$,  where $\bG$ is a  connected reductive group (represented
by  'W:=Group(WF)'), and $σ$ is  a quasi-central automorphism of $\bG$
defined  by <WF>. The element <t> should  be a semisimple element of $\bG$.
The  function returns a  boolean describing whether  $tσ$ is isolated,
that  is whether the Weyl group of  $C_\bG(tσ)^0$ is not in any proper
parabolic subgroup of $W^σ$.

|    gap> WF:=RootDatum("u",6);
    2A5.(q+1)
    gap> l:=QuasiIsolatedRepresentatives(WF);
    [ <0,0,0,0,0,0>, <1/4,0,0,0,0,3/4>, <1/4,1/4,0,0,3/4,3/4>,
      <1/4,1/4,1/4,3/4,3/4,3/4> ]
    gap> List(l,s->IsIsolated(WF,s));
    [ true, false, true, true ]|

"""
CoxeterCosetOps[:IsIsolated] = function (WF, t)
  t=SemisimpleElement(WF[:Rs],SolutionMat(WF[:Y_s],t[:v]*TransposedMat(WF[:pi])))
  return IsIsolated(WF[:Rs], t)
end
