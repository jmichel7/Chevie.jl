module KL
export KLPol, Cpbasis
using Gapjm

"""
 `CriticalPair(W, y, w)` returns the critical pair (z,w) associated to (y,w).

 Let  L (resp. R) be the left (resp. right) descent set. A pair of elements
 (y,w)  of W is called  critical if L(y)⊃ L(w)  and R(y)⊃ R(w). If (y,w) is
 not  critical, y can be  multiplied from the left  (resp. the right) by an
 element  of L(w) (resp. R(w))  which is not in  L(y) (resp. R(y)) until we
 get a pair (z,w) critical. The function returns z. If y<=w then y<=z<=w.

 The significance of this construction is that `KLPol(W,y,w)==KLPol(W,z,w)`
"""
function CriticalPair(W::CoxeterGroup,y,w)::typeof(y)
  Lw=filter(i->isleftdescent(W,w,i),eachindex(coxgens(W)))
  Rw=filter(i->isleftdescent(W,inv(w),i),eachindex(coxgens(W)))
  function cr(y)::typeof(y)
    for s in Lw if !isleftdescent(W,y,s) return cr(coxgens(W)[s]*y) end end
    for s in Rw if !isleftdescent(W,inv(y),s) return cr(y*coxgens(W)[s]) end end
    y
  end
  cr(y)
end

"""
  P(W,y,w) returns the Kazhdan-Lusztig polynomial P_{y,w} of W
"""
function KLPol(W::CoxeterGroup,y,w)::Pol{Int}
  if !bruhatless(W,y,w) return Pol(Int[],0) end
  y=CriticalPair(W,y,w)
  lw=length(W,w)
  if lw-length(W,y)<=2 return Pol(1) end
  d=gets(W->Dict{Tuple{Perm,Perm},Pol{Int}}(),W,:klpol)
  if haskey(d,(w,y)) return  d[(w,y)] end
  s=firstleftdescent(W,w)
  v=coxgens(W)[s]*w
  pol=KLPol(W,coxgens(W)[s]*y,v)+shift(KLPol(W,y,v),1)
  lz=lw-2
  while div(lw-lz,2)<=Pols.degree(pol)
   for z in CoxGroups.elements(W,lz)::Vector{typeof(w)}
      if div(lw-lz,2)<=Pols.degree(pol) && pol.c[div(lw-lz,2)+1]>0 && 
        isleftdescent(W,z,s) && bruhatless(W,y,z)
        let z=z, m=m=KLMue(W,z,v)
        if m!=0 
          pol-=m*shift(KLPol(W,y,z),div(lw-lz,2)) 
        end
        end
      end
    end
    lz-=2
  end
  d[(w,y)]=pol
end

"""
  KLMue(W, y, w) highest coefficient of KLPol(W,y,w)

  KLMue returns the coefficient of highest possible degree (l(w)-l(y)-1)/2
  of  KLPol(W,y,w). This is 0 unless y<=w for the Bruhat order.
"""
function KLMue(W::CoxeterGroup,y,w)
  ly=length(W,y)
  lw=length(W,w)
  if ly>=lw || !bruhatless(W,y,w) return 0 end
  if lw==ly+1 return 1 end
  if any(s->(isleftdescent(W,w,s) && !isleftdescent(W,y,s)) 
    || (isleftdescent(W,inv(w),s) && !isleftdescent(W,inv(y),s)), 
    eachindex(coxgens(W))) 
    return 0
  end
  pol=KLPol(W,y,w)
  if degree(pol)==div(lw-ly-1,2) return pol.c[div(lw-ly+1,2)]
  else return 0 end
end

function QXHalf(H::HeckeAlgebra,x::Perm)
  (-H.sqpara[1]*H.para[1][2])^length(H.W,x)
end

struct HeckeCpElt{P,C}<:HeckeElt{P,C}
  d::SortedPairs{P,C} # has better merge performance than Dict
  H::HeckeAlgebra{C,<:CoxeterGroup}
end

Hecke.clone(h::HeckeCpElt,d)=HeckeCpElt(d,h.H)

Hecke.basename(h::HeckeCpElt)="C'"

function Cpbasis(H::HeckeAlgebra{C,TW})where C where TW<:CoxeterGroup{P} where P
  function f(w::Vector{<:Integer})
    if isempty(w) return HeckeCpElt([one(H.W)=>one(C)],H) end
    HeckeCpElt([eltword(H.W,collect(w))=>one(C)],H)
  end
  f(w::Vararg{Integer})=f(collect(w))
  f(w::P)=f(word(H.W,w))
end

"""
    getCp(H,w)

  return ``C_w`` expressed in the basis T of the Hecke algebra H.

  New implementation JM and FD 3/2000 replacing A.Mathas 12/1994.
  We use the formulae:
  ``C'_w=Σ_{y≤w}P_{y,w}(q)q^{-l(w)/2}T_y``
  and if ``sw<w`` then
  ``C'_s C'_{sw}=C'w+Σ_{y<sw}μ(y,sw)C'y=Σ_{v≤w}μ_v T_v``
  where
  ``μ_v=P_{v,w}(q)q^{-l(w)/2}+Σ_{v≤y≤sw}μ(y,sw)P_{v,y}(q)q^{-l(y)/2}``
 
  It  follows that if  ``deg(μ_v)>=-l(v)`` then ``deg(μ_v)=-l(v)`` with
  leading  coefficient ``μ(v,sw)``  (this happens  exactly for ``y=v`` in
  the sum which occurs in the formula for ``μ_v``).
"""
function getCp(H::HeckeAlgebra,w::P)where P
  T=Tbasis(H)
  W=H.W
  cdict=gets(H,Symbol("C'->T")) do H 
     Dict(one(W)=>T()) end::Dict{P,HeckeTElt{P,typeof(H.sqpara[1])}}
  if haskey(cdict,w) return cdict[w] end
  if w in coxgens(W) 
    return inv(H.sqpara[1])*(T()+T(findfirst(isequal(w),coxgens(W))))
  else
    s=coxgens(W)[firstleftdescent(W,w)]
    res=getCp(H,s)*getCp(H,s*w)
    tmp=zero(T())
    for (e,coef) in res.d 
      if e!=w tmp+=positive_part(coef*QXHalf(H,e))*getCp(H,e) end
    end
    res-=tmp
  end
  cdict[w]=res
end

Hecke.Tbasis(h::HeckeCpElt)=sum(getCp(h.H,e)*c for (e,c) in h.d)
end
