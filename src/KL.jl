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
congruent to ``t_w`` modulo ``∑_{w∈ W}Γ⁻ t_w``.

The  basis `C′_w` can be computed  as follows. We define elements ``R_{x,y}``
of  `ℤ[Γ]` by  ``T_y⁻¹=∑_x \\overline{R_{x,y⁻¹}}  q_x⁻¹T_x``. We  then define
inductively  the Kazhdan-Lusztig  polynomials (in  this general  context we
should  say the  Kazhdan-Lusztig elements  of `ℤ[Γ]`,  which belong  to the
subalgebra  of `ℤ[Γ]` generated by  the `qₛ`) by ``P_{x,w}=τ_{≤(q_w/q_x)^½}
(∑_{x<y≤w}R_{x,y}P_{y,w})``  where `τ`  is the  truncation: ``τ_≤\\nu ∑_{γ∈ Γ}
a_γγ= ∑_{γ≤\\nu}a_γγ``; the induction is thus on decreasing `x` for the Bruhat
order  and  starts  at  ``P_{w,w}=1``.  We  have  then  ``C′_w=∑_y q_w^{-1/2}
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
as   long   as   methods   `bar`:``∑_{γ∈   Γ}   a_γγ↦  ∑_{γ∈  Γ}  a_γγ⁻¹``,
`positive_part`  : ``∑_{γ∈  Γ} a_γγ↦  ∑_{γ≥ 1}  a_γγ`` and `negative_part`:
``∑_{γ∈  Γ}  a_γγ  ↦  ∑_{γ≤  1}  a_γγ``  have been defined on `ℤ[Γ]`. These
operations   will   be   used   internally   by  the  programs  to  compute
Kazhdan-Lusztig bases.

finally, benchmarks on julia 1.0.2
```benchmark
julia> function test_kl(W)
         q=Pol([1],1); H=hecke(W,q^2,rootpara=q)
         C=Cpbasis(H); T=Tbasis(H)
         [T(C(w)) for w in elements(W)]
       end
test_kl (generic function with 1 method)

julia> @btime test_kl(coxgroup(:F,4));
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

julia>@btime test_kl2(coxgroup(:F,4));
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
export KLPol, Cpbasis, LeftCell, LeftCells
using Gapjm

#--------- Meinolf Geck's code for KL polynomials ---------------------------
"""
 `critical_pair(W, y, w)` returns the critical pair z≤w associated to y≤w.

Let  `ℒ` (resp.  `ℛ `)  be the  left (resp.  right) descent  set. A pair of
elements y≤w of W is called critical if `ℒ(y)⊃ ℒ(w)` and `ℛ (y)⊃ ℛ (w)`. If
y≤w is not critical, y can be multiplied from the left (resp. the right) by
an  element of  `ℒ(w)` (resp.  `ℛ (w)`)  which is  not in `ℒ (y)` (resp. `ℛ
(y)`) until we get a critical pair z≤w. The function returns z. If y≤w then
y≤z≤w.

The significance of this construction is that `KLPol(W,y,w)==KLPol(W,z,w)`

```julia-repl
julia> W=coxgroup(:F,4)
F₄

julia> w=longest(W)*gens(W)[1];length(W,w)
23

julia> y=W(1:4...);length(W,y)
4

julia> cr=KL.critical_pair(W,y,w);length(W,cr)
16

julia> Pol(:x);KLPol(W,y,w)
Pol{Int64}: x³+1

julia> KLPol(W,cr,w)
Pol{Int64}: x³+1
```julia-repl

"""
function critical_pair(W::CoxeterGroup,y,w)::typeof(y)
  Lw=filter(i->isleftdescent(W,w,i),eachindex(gens(W)))
  Rw=filter(i->isleftdescent(W,inv(w),i),eachindex(gens(W)))
  function cr(y)::typeof(y)
    for s in Lw if !isleftdescent(W,y,s) return cr(gens(W)[s]*y) end end
    for s in Rw if !isleftdescent(W,inv(y),s) return cr(y*gens(W)[s]) end end
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
    eachindex(gens(W))) 
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
found by the function 'critical_pair'. Whenever the polynomial corresponding
to a critical pair is computed then this pair and the polynomial are stored
in the property `:klpol` of the underlying Coxeter group.

```julia-repl
julia> W=coxgroup(:B,3)
B₃

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
  y=critical_pair(W,y,w)
  lw=length(W,w)
  if lw-length(W,y)<=2 return Pol(1) end
  d=gets(W->Dict{Tuple{Perm,Perm},Pol{Int}}(),W,:klpol)
  if haskey(d,(w,y)) return  d[(w,y)] end
  s=firstleftdescent(W,w)
  v=gens(W)[s]*w
  pol=KLPol(W,gens(W)[s]*y,v)+shift(KLPol(W,y,v),1)
  lz=lw-2
  while div(lw-lz,2)<=degree(pol)
   for z in CoxGroups.elements(W,lz)::Vector{typeof(w)}
      if div(lw-lz,2)<=degree(pol) && pol.c[div(lw-lz,2)+1]>0 && 
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

Hecke.rootpara(H::HeckeAlgebra,x::Perm)=equalpara(H) ?  rootpara(H)[1]^length(H.W,x) : prod(rootpara(H)[word(H.W,x)])

struct HeckeCpElt{P,C,G<:CoxeterGroup}<:HeckeElt{P,C}
  d::ModuleElt{P,C} # has better merge performance than Dict
  H::HeckeAlgebra{C,G}
end

Hecke.clone(h::HeckeCpElt,d)=HeckeCpElt(d,h.H)

Hecke.basename(h::HeckeCpElt)="C'"

function Cpbasis(H::HeckeAlgebra{C,TW})where C where TW<:CoxeterGroup{P} where P
  function f(w::Vector{<:Integer})
    if isempty(w) return HeckeCpElt(ModuleElt(one(H.W)=>one(C)),H) end
    HeckeCpElt(ModuleElt(H.W(w...)=>one(C)),H)
  end
  f(w::Vararg{Integer})=f(collect(w))
  f(w::P)=f(word(H.W,w))
  f(h::HeckeElt)=Cpbasis(h)
end

function Cpbasis(h::HeckeTElt)
  H=h.H
  res=HeckeCpElt(zero(h.d),H)
  while !iszero(h)
    lens=map(x->length(H.W,x[1]),h.d)
    lens=findall(isequal(maximum(lens)),lens)
    tmp=HeckeCpElt(ModuleElt(x[1]=>x[2]*rootpara(H,x[1]) for x in h.d.d[lens]),H)
    res+=tmp
    h-=Tbasis(H)(tmp)
  end
  res
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
  if equalpara(H)
    l=firstleftdescent(W,w)
    s=gens(W)[l]
    if w==s
      return inv(rootpara(H)[l])*(T(s)-H.para[l][2]*one(H))
    else
      res=getCp(H,s)*getCp(H,s*w)
      tmp=zero(H)
      for (e,coef) in res.d 
        if e!=w tmp+=positive_part(coef*rootpara(H,e))*getCp(H,e) end
      end
      res-=tmp
    end
  else
# we follow formula 2.2 in Lusztig's 'Left cells in Weyl groups'
#
# bar(P̄̄_{x,w})-P_{x,w}=∑_{x<y≤w} R_{x,y} P_{y,w}
#
#  where R_{x,y}=bar(T_{y^-1}^{-1}|T_x)
#  where T is the basis with parameters q_s,-q_s^-1
#
# thus we compute P_{x,w} by induction on l(w)-l(x) by
# P_{x,w}=\neg ∑_{x<y≤w} R_{x,y} P_{y,w}
    elm=reduce(vcat,reverse(bruhatless(W,w)))
    coeff=fill(inv(rootpara(H,w)),length(elm))# start with Lusztig  ̃T basis
    f(w)= w==one(W) ? 1 : prod(y->-H.para[y][2],word(W,w))
    for i in 2:length(elm)
      x=elm[i]
      qx=rootpara(H,x)
      z=critical_pair(W,x,w)
      if x!=z coeff[i]=f(z)*inv(f(x))*coeff[findfirst(isequal(z),elm)]
      else 
        coeff[i]=-negative_part(sum(j->
          bar(qx*inv(T(inv(elm[j]))).d[x])*coeff[j],1:i-1))*inv(qx)
      end
    end
    res=HeckeTElt(ModuleElt(map((x,y)->x=>y,elm,coeff);check=true),H)
  end
  cdict[w]=res
end

"""
```julia-repl
julia> W=coxgroup(:B,3)
B₃

julia> Pol(:v);H=hecke(W,v^2,rootpara=v)
Hecke(B₃,v²,rootpara=v)

julia> C=Cpbasis(H)
(::Gapjm.KL.var"#f#10"{Pol{Int64},Perm{Int16},HeckeAlgebra{Pol{Int64},Gapjm.Weyl.FCG{Int16,Int64,PRG{Int64,Int16}}}}) (generic function with 4 methods)

julia> T=Tbasis(H)
(::Gapjm.Hecke.var"#f#25"{Pol{Int64},Perm{Int16},HeckeAlgebra{Pol{Int64},Gapjm.Weyl.FCG{Int16,Int64,PRG{Int64,Int16}}}}) (generic function with 4 methods)

julia> T(C(1,2))
v⁻²T.+v⁻²T₂+v⁻²T₁+v⁻²T₁₂
```
"""
Hecke.Tbasis(h::HeckeCpElt)=sum(getCp(h.H,e)*c for (e,c) in h.d)

Base.:*(a::HeckeCpElt,b::HeckeCpElt)=Cpbasis(Tbasis(a)*Tbasis(b))

#----------------------------------Left cells --------------------------
struct LeftCell{G<:Group}
  group::G
  prop::Dict{Symbol,Any}
## Optional (computed) fields are:
## .reps representatives of the cell
## .duflo the Duflo involution of the cell
## .character the character of the cell
## .a the a-function of the cell
## .elements the elements of the cell
end

Base.copy(c::LeftCell)=LeftCell(c.group,copy(c.prop))

duflo(c::LeftCell)=c.prop[:duflo]

Base.hash(c::LeftCell, h::UInt)=hash(duflo(c),h)

function Base.length(c::LeftCell)
  if haskey(c.prop,:character)
    sum(first.(CharTable(c.group).irr)[c.prop[:character]])
  else length(elements(c))
  end
end

## returns character as list of irred. numbers repeated if multiplicity
function character(c::LeftCell)
  gets(c,:character) do c
    r=representation(c,hecke(c.group))
    cc=HasType.CharRepresentationWords(r,classinfo(c.group)[:classtext])
    ct=CharTable(c.group)
    cc=Chars.decompose(ct,cc)
    char=vcat(map(i->fill(i,cc[i]),1:length(cc))...)
    c.prop[:a]=charinfo(c.group)[:a][char]
    if length(Set(c.prop[:a]))>1 error() else c.prop[:a]=c.prop[:a][1] end
    char
  end
end

function Base.show(io::IO,c::LeftCell)
   print(io,"LeftCell<",c.group,": ")
   if haskey(c.prop,:duflo)
     print(io,"duflo=",joindigits(Weyl.DescribeInvolution(c.group,duflo(c))))
   end
   if haskey(c.prop,:character)
     uc=UnipotentCharacters(c.group)
     ch=c.prop[:character]
     i=findfirst(f->uc.harishChandra[1][:charNumbers][ch[1]] in f[:charNumbers],
                 uc.families)
     f=uc.families[i]
     i=f[:charNumbers][HasType.special(f)]
     i=findfirst(==(i),uc.harishChandra[1][:charNumbers])
     p=findfirst(==(i),ch)
     p=vcat([[i,1]],HasType.Collected(vcat(ch[1:p-1],ch[p+1:end])))
     print(io," character=",join(
       map(p)do v
        (v[2]!=1 ? string(v[2]) : "")*charnames(io,c.group)[v[1]]
      end,"+"))
   end
   print(io,">")
 end
 
# br is a braid_relation of W
# returns corresponding * op on w if applicable otherwise returns w
function leftstar(W,br,w)
  # st is a left or right member of a braidrelation of W
  # leftdescents(W,w) contains st[1] and not st[2]
  # returns corresponding * operation applied to w
  function leftstarNC(W,st,w)
    rst=reflections(W)[st]
    w0=prod(rst)
    i=1
    while true
      w=rst[i]*w
      i+=1
      w0*=rst[i]
      if !isleftdescent(W,w,st[i]) break end
    end
    w0*w
  end
  l,r=map(x->restriction(W)[x],br)
  if isleftdescent(W,w,l[1]) isleftdescent(W,w,r[1]) ? w : leftstarNC(W,l,w)
  else isleftdescent(W,w,r[1]) ? leftstarNC(W,r,w) : w
  end
end

# List of functions giving all possible left * images of w
leftstars(W)=map(st->(w->leftstar(W,st,w)),
                 filter(r->length(r[1])>2,braid_relations(W)))

function Gapjm.elements(c::LeftCell)
  gets(c,:elements) do c
    elements=orbit(leftstars(c.group),duflo(c);action=(x,f)->f(x))
    for w in c.prop[:reps]
      append!(elements,orbit(leftstars(c.group),w;action=(x,f)->f(x)))
    end
    Set(elements)
  end
end

Gapjm.words(c::LeftCell)=word.(Ref(c.group),elements(c))

Base.:(==)(a::LeftCell,b::LeftCell)=duflo(a)==duflo(b)

Base.in(w,c::LeftCell)=w in elements(c)

#F  KLMueMat( <W>, <list> )  . . . . (symmetrized) matrix of leading 
#F coefficients of Kazhdan-Lusztig polynomials of elements in a given list
##
function KLMueMat(W,c)
  w0c=longest(W).*c
  lc=length.(Ref(W),c)
  m=zeros(Int,length(c),length(c))
  for k in 0:maximum(lc)-minimum(lc)
    for i in eachindex(c)
      for j in filter(j->lc[i]-lc[j]==k,eachindex(c))
        if lc[i]+lc[j]>W.N m[i,j]=KLMue(W,w0c[i],w0c[j])
        else m[i,j]=KLMue(W,c[j],c[i])
        end
        m[j,i]=m[i,j]
      end
    end
  end
  return m
end

function Mu(c::LeftCell)
  gets(c,:mu) do
    KLMueMat(c.group,elements(c))
  end
end

function Chars.representation(c::LeftCell,H)
  W=H.W
  if !equalpara(H)
    error("cell representations for unequal parameters not yet implemented")
  else v=rootpara(H)[1]
  end
  WGraphToRepresentation(semisimplerank(c.group),WGraph(c),v)
end

# returns right star operation st (a BraidRelation) of LeftCell c
function RightStar(st,c)
  res=copy(c)
  W=c.group
  rs(w)=leftstar(W,st,w^-1)^-1
  if haskey(c.prop,:duflo)
   res.prop[:duflo]=rs(rs(duflo(c))^-1)
  end
  if haskey(c.prop,:reps) res.prop[:reps]=rs.(c.prop[:reps]) end
  if haskey(c.prop,:elements)
    res.prop[:elements]=rs.(c.prop[:elements])
    n=sortperm(res.prop[:elements])
    res.prop[:elements]=res.prop[:elements][n]
    if haskey(c.prop,:mu) res.prop[:mu]=c.prop[:mu][n,n] end
    if haskey(c.prop,:graph) res.prop[:orderGraph]=c.prop[:orderGraph][n] end
  end
  res
end

# Fo an irreducible type, reps contain:
# .duflo,  .reps: elements of W represented as images of simple roots
# .character: decomposition of left cell in irreducibles
function LeftCellRepresentatives(W)
  res=map(refltype(W))do t
#   R=ReflectionGroup(t)
    R=rootdatum(cartan(t.prop)) # above implemented for now like that
    rr=getchev(t,:KLeftCellRepresentatives)
    if isnothing(rr) return nothing end
    return map(rr)do r
      r=copy(r)
      function f(l)
        m=permutedims(toM(R.rootdec[l]))
        w=Perm(R.rootdec,Ref(m).*R.rootdec)
        inclusion(W)[t[:indices][word(R,w)]]
      end
      r[:duflo]=f(r[:duflo])
      if isempty(r[:reps]) r[:reps]=Vector{Int}[]
      else r[:reps]=f.(r[:reps])
      end
      push!(r[:reps],r[:duflo])
      r
    end
  end
  if isempty(res) return res end
  if nothing in res return nothing end
  n=getchev(W,:NrConjugacyClasses)
  return map(Cartesian(res...)) do l
    duflo=W(vcat(map(x->x[:duflo],l)...)...)
    reps=map(v->W(vcat(v...)...),Cartesian(map(x->x[:reps],l)...))
    reps=setdiff(reps,[duflo])
    character=map(p->HasType.PositionCartesian(n,p),Cartesian(map(x->x[:character],l)...))
    a=charinfo(W)[:a][character]
    if length(Set(a))>1 error() else a=a[1] end
    return LeftCell(W,Dict{Symbol,Any}(:duflo=>duflo,:reps=>reps,
                 :character=>character,:a=>a))
  end
end

InfoChevie(x...)=println(x...)
function OldLeftCellRepresentatives(W)
  st=map(st->(c->RightStar(st,c)),filter(r->length(r[1])>2,braid_relations(W)))
  rw=groupby(x->leftdescents(W,x^-1),elements(W))
  rw=[(rd=k,elements=Set(v)) for (k,v) in rw]
  sort!(rw;by=x->length(x.elements))
  cells0=LeftCell{typeof(W)}[]
  while length(rw)>0
    c=collect(rw[1].elements)
    InfoChevie("#I R(w)=",rw[1].rd," : #Elts=",length(c))
    mu=KLMueMat(W,c)
    n=1:length(c)
    Lleq=[x==y || (mu[x,y]!=0 && any(i->isleftdescent(W,c[x],i) && 
       !isleftdescent(W,c[y],i),eachindex(gens(W)))) for x in n, y in n]
    Lleq=transitive_closure(Lleq)
    m=Set(n[Lleq[i,:].&Lleq[:,i]] for i in n)
    x=[LeftCell(W,Dict{Symbol,Any}(:elements=>c[d],:mu=>mu[d,d])) for d in m]
    while length(x)>0
      c=x[1]
      n=sortperm(c.prop[:elements])
      c.prop[:elements]=c.prop[:elements][n]
      c.prop[:mu]=c.prop[:mu][n,n]
      i=filter(x->isone(x^2),c.prop[:elements])
      if length(i)==1 c.prop[:duflo]=i[1]
      else m=map(x->length(W,x)-2*degree(KLPol(W,one(W),x)),i)
        p=argmin(m)
        c.prop[:a]=m[p]
        c.prop[:duflo]=i[p] # Duflo involutions minimize Delta
      end
      i=filter(x->!(c.prop[:duflo] in x),
               orbits(leftstars(W),c.prop[:elements];action=(x,f)->f(x)))
      c.prop[:reps]=first.(i)
      push!(cells0,c)
      n=orbit(st,c;action=(x,f)->f(x))
      InfoChevie(", ",length(n)," new cell")
      if length(n)>1 InfoChevie("s") end
      for e in n
        rd=leftdescents(W,duflo(e))
        i=findfirst(x->x.rd==rd,rw)
        if i==1 x=filter(c->!(c.prop[:elements][1] in elements(e)),x)
        elseif i!=nothing
	  setdiff!(rw[i].elements,elements(e))
        end
      end
    end
    InfoChevie(" \n")
    rw=filter(x->length(x.elements)>0,rw)
    rw=rw[2:end]
    sort!(rw,by=x->length(x.elements))
  end
  cells0
end
  
function cellreps(W)
  gets(W,:cellreps) do W
    cc=LeftCellRepresentatives(W)
    if isnothing(cc) cc=OldLeftCellRepresentatives(W) end
    cc
  end
end
    
"""
  `LeftCells(W[,i])` left cells of `W` [in `i`-th 2-sided cell]
  for the 1-parameter Hecke algebra
"""
function LeftCells(W,i=0)
  cc=cellreps(W)
  if !iszero(i)
    uc=UnipotentCharacters(W)
    cc=filter(c->uc.harishChandra[1][:charNumbers][character(c)[1]]
                                in uc.families[i][:charNumbers],cc)
  end
  st=map(st->(c->RightStar(st,c)),filter(r->length(r[1])>2,braid_relations(W)))
  d=map(c->orbit(st,c;action=(x,f)->f(x)),cc)
  length(d)==1 ? d[1] : union(d...)
end

# gens is a list each element of which can operate on element e
# returns minimal word w such that Composition(gens[w]) applied to e 
# satifies cond
function MinimalWordProperty(e,gens::Vector,cond::Function;action::Function=^)
  if cond(e) return Int[] end
  elements=[e]
  nbLength=[1]
  cayleyGraph=[Int[]]
  bag=Set(elements)
  InfoChevie("#I ")
  while true
    new=map(1+sum(nbLength[1:end-1]):length(elements))do h
      map(g->[action(elements[h],gens[g]),[h,g]],eachindex(gens))
    end
    new=first.(values(groupby(first,vcat(new...))))
    new=filter(x->!(x[1] in bag),new)
    append!(cayleyGraph,map(x->x[2],new))
    new=first.(new)
    p=findfirst(cond,new)
    if !isnothing(p) InfoChevie("\n")
      res=Int[] 
      p=length(elements)+p
      while p!=1 push!(res,cayleyGraph[p][2])
        p=cayleyGraph[p][1] 
      end
      return res
    end
    if all(x->x in elements,new) error("no solution") end
    append!(elements,new)
    bag=union(bag,new)
#   if Length(new)>10 then 
       InfoChevie(length(new)," ")
#   fi;
    push!(nbLength,length(new))
  end
end

# LeftCell containing w
function LeftCell(W,w)
  l=cellreps(W)
  sst=filter(r->length(r[1])>2,braid_relations(W))
  word=MinimalWordProperty(w,map(st->(w->leftstar(W,st,w^-1)^-1),sst),
     w->any(c->w in c,l);action=(x,f)->f(x))
  v=w
  for g in reverse(word) v=leftstar(W,sst[g],v^-1)^-1 end
  cell=l[findfirst(c->v in c,l)]
  for g in word cell=RightStar(sst[g],cell) end
  return cell
end

# WGraph of LeftCell c
function WGraph(c::LeftCell)
  gets(c,:graph) do c
    e=elements(c)
    mu=Mu(c)
    n=length(e)
    nodes=leftdescents.(Ref(c.group),e)
    p=sortperm(nodes)
    nodes=nodes[p]
    mu=mu[p,p]
    c.prop[:orderGraph]=p
    nodes=vcat(map(HasType.Collected(nodes)) do p
        p[2]==1  ? [p[1]] : [p[1],p[2]-1]
        end...)
    graph=[nodes,[]]
    l=vcat(map(i->map(j->[mu[i,j],mu[j,i],i,j],1:i-1),1:n)...)
    l=filter(x->x[1]!=0 || x[2]!=0,l)
    if isempty(l) return graph end
    l=groupby(x->x[[1,2]],l)
    for u in values(l)
      if u[1][1]==u[1][2] value=u[1][1] else value=u[1][[1,2]] end
      w=[value,[]]
      s=groupby(first,map(x->x[[3,4]],u))
      for k in values(s) push!(w[2],vcat([k[1][1]],map(x->x[2],k))) end
      push!(graph[2],w)
    end
    graph
  end
end

end
