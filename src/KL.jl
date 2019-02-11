#A  good example to see how long  the programs will take for computations in
#big Coxeter groups is the following:
#
#|    gap> W:=CoxeterGroup("D",5);;
#    gap> LeftCells(W);;|
#
#which  takes `10` seconds cpu time on 3Ghz computer. The computation of all
#Kazhdan-Lusztig  polynomials  for  type  `F_4`  takes  a  bit more than~`1`
#minute.  Computing the Bruhat order is a bottleneck for these computations;
#they can be speeded up by a factor of two if one does:
#
#|    gap> ReadChv("contr/brbase");
#    gap> BaseBruhat(W);;|
#
#after  which the computation  of the Bruhat  order will be  speeded up by a
#large factor.
#
#However,  Alvis'  computation  of  the  Kazhdan--Lusztig polynomials of the
#Coxeter  group of type  `H_4` in a  computer algebra system  like GAP would
#still take many hours. For such applications, it is probably more efficient
#to use a special purpose program like the one provided by F. DuCloux DuC91.
"""
This  module ports Chevie functionality for Kazhdan-Lusztig polynomials and
bases.

Let  `ℋ ` be  the Iwahori-Hecke algebra  of a Coxeter  system `(W,S)`, with
quadratic  relations `(Tₛ-uₛ₀)(Tₛ-uₛ₁)=0`  for `s∈  S`. If  `-uₛ₀uₛ₁` has a
square  root  `wₛ`,  we  can  scale  the  basis  `Tₛ`  to  get  a new basis
`tₛ=-Tₛ/wₛ`    with   quadratic    relations   `(tₛ-vₛ)(tₛ+vₛ⁻¹)=0`   where
`vₛ=wₛ/uₛ₁`.   The  most  general  case   when  Kazhdan-Lusztig  bases  and
polynomials  can be defined is when the parameters `vₛ` belong to a totally
ordered  abelian group `Γ`  for multiplication, see  Lus83. We set `Γ⁺= {γ∈
Γ∣γ>0}` and `Γ⁻={γ⁻¹∣γ∈ Γ⁺}={γ∈ Γ∣γ<0}`.

Thus  we assume `ℋ ` defined over the ring `ℤ[Γ]`, the group algebra of `Γ`
over  `ℤ`, and the quadratic  relations of `ℋ `  associate to each `s∈ S` a
`vₛ∈  Γ⁺` such that  `(tₛ-vₛ)(tₛ+vₛ⁻¹)=0`. We also  set `qₛ=vₛ²` and define
the  basis `Tₛ=vₛtₛ` with quadratic relations `(Tₛ-qₛ)(Tₛ+1)=0`; for `w∈ W`
with reduced expression `w=s₁…sₙ` we define `q_w∈ Γ⁺` by
`q_w^½=v_{s₁}…v_{sₙ}` and let `q_w=(q_w^½)²`.

We  define the bar involution on `ℋ `  by linearity: on `ℤ[Γ]` we define it
by  ``\\overline{∑_{γ∈ Γ}a_γγ}= ∑_{γ∈ Γ} a_γ γ⁻¹`` and we extend it to `ℋ `
by  ``\\overline  Tₛ=Tₛ⁻¹  ``.  Then  the  Kazhdan-Lusztig  basis `C′_w` is
defined  as  the  only  basis  of  `ℋ  `  stable  by the bar involution and
congruent to `t_w` modulo `∑_{w∈ W}Γ⁻ t_w`.

The  basis `C′_w` can be computed  as follows. We define elements `R_{x,y}`
of  `ℤ[Γ]` by  `T_y⁻¹=∑_x \\overline{R_{x,y⁻¹}}  q_x⁻¹T_x`. We  then define
inductively  the Kazhdan-Lusztig  polynomials (in  this general  context we
should  say the  Kazhdan-Lusztig elements  of `ℤ[Γ]`,  which belong  to the
subalgebra  of `ℤ[Γ]` generated by  the `qₛ`) by ``P_{x,w}=τ_{≤(q_w/q_x)^½}
(∑_{x<y≤w}R_{x,y}P_{y,w})``  where `τ`  is the  truncation: ``τ_≤ν ∑_{γ∈ Γ}
a_γγ= ∑_{γ≤ν}a_γγ``; the induction is thus on decreasing `x` for the Bruhat
order  and  starts  at  `P_{w,w}=1`.  We  have  then  ``C′_w=∑_y q_w^{-1/2}
P_{y,w}T_y``.

The  Chevie code  for the  Kazhdan-Lusztig bases  `C`, `D` and their primed
versions, has been initially written by Andrew Mathas around 1994, who also
contributed  to  the  design  of  the programs dealing with Kazhdan-Lusztig
bases. He also implemented some other bases, such as the Murphy basis which
can  be  found  in  the  Chevie  contributions  directory. The code for the
unequal  parameters  case  has  been  written  around  1999  by F.Digne and
J.Michel. The other Kazhdan-Lusztig bases are computed in terms of the `C′`
basis.

When  the `ℤ[Γ]` is a  Laurent polynomial ring the  bar operation is taking
the  inverse of  the variables,  and truncation  is keeping terms of degree
smaller or equal to that of `ν`. It is possible to use arbitrary groups `Γ`
as  long as  methods `bar`,  `positive_part` and  `negative_part` have been
defined on them. which perform the operations respectively ``∑_{γ∈ Γ} a_γγ↦
∑_{γ∈  Γ} a_γγ⁻¹``, ``∑_{γ∈ Γ}  a_γγ↦ ∑_{γ≥ 1} a_γγ``  and ``∑_{γ∈ Γ} a_γγ↦
∑_{γ≤  1} a_γγ`` on  elements 'p' of  `ℤ[Γ]`. The operations  above will be
used internally by the programs to compute Kazhdan-Lusztig bases.

finally, benchmarks on julia 1.0.2
```benchmark
julia> function test_kl(W)
         q=Pol([1],1); H=hecke(W,q^2,q)
         C=Cpbasis(H); T=Tbasis(H)
         [T(C(w)) for w in elements(W)]
       end
test_kl (generic function with 1 method)

julia> @btime test_kl(WeylGroup(:F,4));
2.265 s (22516606 allocations: 1.81 GiB)
```
Compare to GAP3 where the following function takes 11s for F4
```
test_kl:=function(W)local q,H,T,C;
  q:=X(Rationals);H:=Hecke(W,q^2,q);
  T:=Basis(H,"T");C:=Basis(H,"C'");
  List(Elements(W),e->T(C(e)));
end;
```
Another benchmark:
```benchmark
function test_kl2(W)
  el=elements(W)
  [KLPol(W,x,y) for x in el, y in el]
end

test_kl2 (generic function with 1 method)

julia>@btime test_kl2(WeylGroup(:F,4));
  8s (97455915 allocations: 6.79 GiB)
```
Compare to GAP3 where the following function takes 42s for F4
```
test_kl2:=function(W)local el;
  el:=Elements(W);
  List(el,x->List(el,y->KazhdanLusztigPolynomial(W,x,y)));
end;
```
"""
module KL
export KLPol, Cpbasis
using Gapjm

#--------- Meinolf Geck's code for KL polynomials ---------------------------
"""
 `CriticalPair(W, y, w)` returns the critical pair z≤w associated to y≤w.

Let  `ℒ` (resp.  `ℛ `)  be the  left (resp.  right) descent  set. A pair of
elements y≤w of W is called critical if `ℒ(y)⊃ ℒ(w)` and `ℛ (y)⊃ ℛ (w)`. If
y≤w is not critical, y can be multiplied from the left (resp. the right) by
an  element of  `ℒ(w)` (resp.  `ℛ (w)`)  which is  not in `ℒ (y)` (resp. `ℛ
(y)`) until we get a critical pair z≤w. The function returns z. If y≤w then
y≤z≤w.

The significance of this construction is that `KLPol(W,y,w)==KLPol(W,z,w)`

```julia-repl
julia> W=WeylGroup(:F,4)
W(F₄)

julia> w=longest(W)*coxgens(W)[1];length(W,w)
23

julia> y=element(W,1:4);length(W,y)
4

julia> cr=KL.CriticalPair(W,y,w);length(W,cr)
16

julia> Pol(:x);KLPol(W,y,w)
x³+1

julia> KLPol(W,cr,w)
x³+1
```julia-repl

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
  KLMue(W, y, w) highest coefficient of KLPol(W,y,w)

KLMue returns the coefficient of highest possible degree (l(w)-l(y)-1)/2
of  KLPol(W,y,w). This is 0 unless y≤w for the Bruhat order.
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

"""
  KLPol(W,y,w) returns the Kazhdan-Lusztig polynomial P_{y,w} of W

To  compute Kazhdan-Lusztig polynomials in  the one-parameter case it seems
that  the best  approach still  is by  using the  recursion formula  in the
original  article KL79. One can first run  a number of standard checks on a
given  pair  of  elements  to  see  if the computation of the corresponding
polynomial  can be reduced to a similar computation for elements of smaller
length. One such check involves the notion of critical pairs (cf. Alv87): a
pair  of elements `w₁,w₂∈  W` such that  `w₁≤w₂` is *critical*  if `ℒ(w₂) ⊆
ℒ(w₁)`  and `ℛ (w₂)⊆ ℛ (w₁)`, where `ℒ`  and `ℛ ` denote the left and right
descent  set, respectively.  Now if  `y≤w ∈  W` are arbitrary elements then
there   always  exists  a  critical  pair   `z≤w`  with  `y≤z≤w`  and  then
`P_{y,w}=P_{z,w}`.  Given two elements `y` and `w`, such a critical pair is
found by the function 'CriticalPair'. Whenever the polynomial corresponding
to a critical pair is computed then this pair and the polynomial are stored
in the property `:klpol` of the underlying Coxeter group.

```julia-repl
julia> W=WeylGroup(:B,3)
W(B₃)

julia> map(i->map(x->KLPol(W,one(W),x),elements(W,i)),1:W.N)
9-element Array{Array{Pol{Int64},1},1}:
 [1, 1, 1]                       
 [1, 1, 1, 1, 1]                 
 [1, 1, 1, 1, 1, 1, 1]           
 [1, 1, 1, x+1, 1, 1, 1, 1]      
 [x+1, 1, 1, x+1, x+1, 1, x+1, 1]
 [1, x+1, 1, x+1, x+1, x²+1, 1]  
 [x+1, x+1, x²+x+1, 1, 1]        
 [x²+1, x+1, 1]                  
 [1]
```
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

#---------- JM & FD code for the C' basis -------------------------------------

QXHalf(H::HeckeAlgebra,x::Perm)=H.sqpara[1]^length(H.W,x)

struct HeckeCpElt{P,C,G<:CoxeterGroup}<:HeckeElt{P,C}
  d::SortedPairs{P,C} # has better merge performance than Dict
  H::HeckeAlgebra{C,G}
end

Hecke.clone(h::HeckeCpElt,d)=HeckeCpElt(d,h.H)

Hecke.basename(h::HeckeCpElt)="C'"

function Cpbasis(H::HeckeAlgebra{C,TW})where C where TW<:CoxeterGroup{P} where P
  function f(w::Vector{<:Integer})
    if isempty(w) return HeckeCpElt([one(H.W)=>one(C)],H) end
    HeckeCpElt([element(H.W,collect(w))=>one(C)],H)
  end
  f(w::Vararg{Integer})=f(collect(w))
  f(w::P)=f(word(H.W,w))
end

"""
    getCp(H,w)

return ``C_w`` expressed in the basis T of the Hecke algebra H.

Implementation JM and FD 1999. We use the formulae:
``C'_w=Σ_{y≤w}P_{y,w}(q)q^{-l(w)/2}T_y``
and if ``sw<w`` then
``C'ₛ C'_{sw}=C'w+Σ_{y<sw}μ(y,sw)C'y=Σ_{v≤w}μᵥ Tᵥ``
where
``μᵥ=P_{v,w}(q)q^{-l(w)/2}+Σ_{v≤y≤sw}μ(y,sw)P_{v,y}(q)q^{-l(y)/2}``

It  follows that if ``deg(μᵥ)>=-l(v)``  then ``deg(μᵥ)=-l(v)`` with leading
coefficient  ``μ(v,sw)`` (this happens exactly for ``y=v`` in the sum which
occurs in the formula for ``μᵥ``).
"""
function getCp(H::HeckeAlgebra{C,G},w::P)where {P,C,G}
  T=Tbasis(H)
  W=H.W
  cdict=gets(H,Symbol("C'->T")) do H 
  Dict(one(W)=>one(H)) end::Dict{P,HeckeTElt{P,C,G}}
  if haskey(cdict,w) return cdict[w] end
  l=firstleftdescent(W,w)
  s=coxgens(W)[l]
  if w==s
    return inv(H.sqpara[l])*(one(H)-inv(H.para[l][2])*T(s))
  else
    res=getCp(H,s)*getCp(H,s*w)
    tmp=zero(H)
    for (e,coef) in res.d 
      if e!=w tmp+=positive_part(coef*QXHalf(H,e))*getCp(H,e) end
    end
    res-=tmp
  end
  cdict[w]=res
end

"""
```julia-repl
julia> W=WeylGroup(:B,3)
W(B₃)

julia> Pol(:v);H=hecke(W,v^2,v)
Hecke(W(B₃),v²,v)

julia> C=Cpbasis(H)
(::getfield(Gapjm.KL, Symbol("#f#10")){Pol{Int64},Perm{Int16},HeckeAlgebra{Pol{Int64},WeylGroup{Int16}}}) (generic function with 3 methods)

julia> T=Tbasis(H)
(::getfield(Gapjm.Hecke, Symbol("#f#20")){Pol{Int64},Perm{Int16},HeckeAlgebra{Pol{Int64},WeylGroup{Int16}}}) (generic function with 4 methods)

julia> T(C(1,2))
v⁻²T.+v⁻²T₂+v⁻²T₁+v⁻²T₁₂
```
"""
Hecke.Tbasis(h::HeckeCpElt)=sum(getCp(h.H,e)*c for (e,c) in h.d)
end
