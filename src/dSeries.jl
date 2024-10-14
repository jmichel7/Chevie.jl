"""
`d`-Harish-Chandra   series  describe  unipotent  `l`-blocks  of  a  finite
reductive  group ``𝐆(𝔽_q)`` for ``l|Φ_d(q)`` (at least, when `l` is not too
small which means mostly not a bad prime for `𝐆`). Some of the facts stated
below  are still partly conjectural, we do not try to distinguish precisely
what has been established and what is still conjectural.

If  `(𝐋,λ)` is  a `d`-cuspidal  pair then  the constituents  of the Lusztig
induced  ``R_𝐋^𝐆(λ)`` are called a `d`-Harish-Chandra series; they form the
unipotent part of an `l`-block of ``𝐆^F``. It is conjectured (and proven in
some   cases)  that  the  ``𝐆^F``-endomorphism   algebra  of  the  `l`-adic
cohomology  of the  variety `𝐗`  which defines  the Lusztig  induction is a
`d`-cyclotomic Hecke algebra ``H_𝐆(𝐋,λ)`` for the group
``W_𝐆(𝐋,λ):=N_𝐆(𝐋,λ)/𝐋``,  which  is  a  complex  reflection group --- here
`d`-cyclotomic  means that the parameters  of ``H_𝐆(𝐋,λ)`` are monomials in
`q`  and that ``H_𝐆(𝐋,λ)``  specializes to the  algebra of ``W_𝐆(𝐋,λ)`` for
``q↦ζ_d``.

It  follows that the decomposition of the  Lusztig induction is of the form
``R_𝐋^𝐆(λ)=∑_{ϕ∈Irr(W_𝐆(𝐋,λ))}(-1)^{nᵩ} ϕ(1)γᵩ,`` where `γᵩ` is a unipotent
character   of  `𝐆^F`  attached  to  `ϕ`  and  where  `nᵩ`  is  the  degree
``H^{nᵩ}_c(𝐗)``  where  `γᵩ`  occurss;  and  further  for  any  `ϕ` we have
``R_𝐋^𝐆(λ)(1)=  (-1)^{nᵩ} γᵩ(1)Sᵩ`` where `Sᵩ` is  the Schur element of the
character  of  ``H_𝐆(𝐋,λ)``  which  deforms  to  `ϕ`. The function |Series|
allows to explore a `d`-Harish-Chandra series.

```julia-repl
julia> W=rootdatum("3D4")
³D₄

julia> l=cuspidal_data(W,3)
2-element Vector{@NamedTuple{levi::Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}, cuspidal::Int64, d::Root1}}:
 (levi = ³D₄, cuspidal = 8, d = ζ₃)
 (levi = ³D₄₍₎=Φ₃², cuspidal = 1, d = ζ₃)

julia> Series(W,l[2]...)
ζ₃-series R^³D₄_{³D₄₍₎=Φ₃²}(λ==Id)  H_G(L,λ)==hecke(G₄,Vector{Mvp{Cyc{Int64}, Int64}}[[ζ₃q², ζ₃, ζ₃q]])
┌─┬───────────────────────┐
│ │    γᵩ    φ  ε family #│
├─┼───────────────────────┤
│1│  φ₁‚₀ φ₁‚₀  1        1│
│2│  φ₁‚₆ φ₁‚₄  1        2│
│3│  φ₂‚₂ φ₁‚₈ -1        5│
│6│ φ″₁‚₃ φ₂‚₅  1        4│
│5│ φ′₁‚₃ φ₂‚₃ -1        3│
│7│  φ₂‚₁ φ₂‚₁ -1        5│
│4│³D₄[1] φ₃‚₂  1        5│
└─┴───────────────────────┘
```
Above  we explore the 3-series corresponding  to ``R_𝐓^𝐆(Id)`` where `𝐆` is
the  triality group  and `𝐓`  is the  torus of  type `(q²+q+1)²`. The group
``W_𝐆(𝐓)``  is the complex reflection group `G₄`. The displays shows in the
column   'γᵩ'  the  name  of   the  unipotent  characters  constituents  of
``R_𝐓^𝐆(Id)``,  and in the  first column the  number of these characters in
the  list  of  unipotent  characters.  In  the  column  'φ' the name of the
character  of ``W_𝐆(𝐓)`` corresponding  to the unipotent  character `γᵩ` is
shown;  in the column  'ε' we show  the sign ``(-1)^{nᵩ}``.  Finally in the
last column we show in which family of unipotent characters is `γᵩ`.

The theory of `d`-Harish-Chandra series can be generalized to spetsial complex
reflection groups using some axioms. We show below an example.

```julia-repl
julia> W=complex_reflection_group(4)
G₄

julia> l=cuspidal_data(W,3)
5-element Vector{@NamedTuple{levi::Spets{PRSG{Cyc{Rational{Int64}}, Int16}}, cuspidal::Int64, d::Root1}}:
 (levi = G₄, cuspidal = 3, d = ζ₃)
 (levi = G₄, cuspidal = 6, d = ζ₃)
 (levi = G₄, cuspidal = 7, d = ζ₃)
 (levi = G₄, cuspidal = 10, d = ζ₃)
 (levi = G₄₍₎=Φ₁Φ′₃, cuspidal = 1, d = ζ₃)

julia> Series(W,l[5]...)
ζ₃-series R^G₄_{G₄₍₎=Φ₁Φ′₃}(λ==Id)  W_G(L,λ)==Z₆
┌─┬────────────────────────────────────┐
│ │   γᵩ φ(mod 3)  ε parameter family #│
├─┼────────────────────────────────────┤
│1│ φ₁‚₀        1  1      ζ₃q²        1│
│5│ φ₂‚₃       ζ₆  1      -ζ₃q        2│
│2│ φ₁‚₄       ζ₃ -1        ζ₃        4│
│8│ Z₃:2       -1 -1     -ζ₃²q        2│
│9│Z₃:11      ζ₃² -1       ζ₃²        4│
│4│ φ₂‚₅      ζ₆⁵ -1       -ζ₃        4│
└─┴────────────────────────────────────┘
```

Above  we explore the `3`-series corresponding  to the trivial character of
the  torus of type `(q-1)(q-ζ₃)`. For cyclic groups ``W_𝐆(𝐋,λ)`` we display
the  parameters in  the table  since they  are associated  to characters of
``W_𝐆(𝐋,λ)``. Finally the mention '(mod 3)' which appears in the 'φ' column
means that in this case the axioms leave an ambiguity in the correspondence
between  unipotent  characters  `γᵩ`  and  characters  `ϕ` (as well as with
parameters):  the correspondence is known only up to a translation by 3 (in
this case, the same as a global multiplication of all `ϕ` by `-1`).

Finally,  we should note that  if the reflection group  or coset `W` is not
defined  over the integers,  what counts is  not cyclotomic polynomials but
factors  of them  over the  field of  definition of  `W`. In this case, one
should not give as argument an integer `d` representing ``ζ_d`` but specify
a  root of unity. For instance, in the above case we get a different answer
with:

```julia-repl
julia> cuspidal_data(W,E(3,2))
5-element Vector{@NamedTuple{levi::Spets{PRSG{Cyc{Rational{Int64}}, Int16}}, cuspidal::Int64, d::Root1}}:
 (levi = G₄, cuspidal = 2, d = ζ₃²)
 (levi = G₄, cuspidal = 5, d = ζ₃²)
 (levi = G₄, cuspidal = 7, d = ζ₃²)
 (levi = G₄, cuspidal = 10, d = ζ₃²)
 (levi = G₄₍₎=Φ₁Φ″₃, cuspidal = 1, d = ζ₃²)
```
"""
module dSeries
using ..Chevie
export Series, ennola

function SpetsEnnola(t::TypeIrred;sperm=true)
  if haskey(t,:orbit) z=gcd(first.(degrees(t))) else z=gcd(degrees(t)) end
  ξ=E(z,-1)
  uc=getchev(t,:UnipotentCharacters)
  ff=improve_type(uc[:families])
  fd=fill(Pol(0),length(uc[:a]))
  l=haskey(uc,:almostHarishChandra) ? uc[:almostHarishChandra] : uc[:harishChandra]
  fd[charnumbers(l[1])]=fakedegrees(reflection_group(t),Pol())

  # positions-with-sign in l where o or -o appears
  findallsgn(l,o)=vcat(findall(==(o),l),-findall(==(-o),l))

  # EnnolaBete[i]: possible action of ξ-Ennola on i-th family
  # list of possible destinations for each char, taking just degree in account
  EnnolaBete=map(ff)do f
    m=f.fourierMat;if !(m isa Matrix) m=toM(m) end
    ud=CycPol.(m'*fd[f.charNumbers])
    map(p->findallsgn(ud,subs(p,Pol([ξ],1))),ud)
  end

  for h in uc[:harishChandra] 
   fd[charnumbers(h)]=
    fakedegrees(reflection_group(Uch.maketype(h[:relativeType])), Pol())
  end
  ΩΧ=map(f->f(ξ)//f(1),fd)# List of ω_χ(ξ) for character sheaves (ρ,χ)
  eig=eigen(ff)

  function predeigen(i) # if Ennola_ξ(U_χ=ρ_i)=ρ returns deduced Frob(ρ)
    if i>nconjugacy_classes(t) error("only for principal series") end
    ΩΧ[i]^-1*E(z^2,uc[:a][i]+uc[:A][i])*eig[i]
  end

  # EnnolaBE(i,f): More sophisticated than EnnolaBete: take in account predeigen
  function EnnolaBE(f,e)
    for j in 1:length(f)
      if f.charNumbers[j] in charnumbers(uc[:harishChandra][1])
        e[j]=filter(k->eigen(f)[abs(k)]==predeigen(f.charNumbers[j]),e[j])
      end
    end
    e
  end

  l=map(eachindex(ff),ff)do i,f
    poss=EnnolaBete[i]
    if prod(length,poss)>1 poss=EnnolaBE(f,poss) end
    A=fusion_algebra(f)
    b=basis(A)
    res=Tuple{Int,SPerm{Int16}}[]
    for i in eachindex(b)
      p=SPerm(b[i]*b,b)
      if p!==nothing
        if all(j->j^p in poss[j],eachindex(b)) push!(res,(i,p)) end
        p=SPerm(-p.d)
        if all(j->j^p in poss[j],eachindex(b)) push!(res,(-i,p)) end
      end
    end
    res
  end
  if !sperm return map(x->first.(x),l) end
  l=map(x->last.(x),l)
  l=cartesian(l...)[1]
  res=fill(0,length(uc[:a]))
  for (i,p) in pairs(l)
    res[ff[i].charNumbers]=invpermute(ff[i].charNumbers,p)
  end
  SPerm(res)
end

"""
`ennola(W[,z1])`

Let  `W` be an irreducible spetsial reflection  group or coset, and `z` the
generator  of the center of `W`, viewed as  a root of unity. Let `𝔾` be the
spets  attached to  `W`. A  property checked  case-by case  is that,  for a
unipotent  character `γ` of `𝔾` with polynomial generic degree `deg γ(q)` ,
`deg  γ(zq)` is equal to `±deg γ'(q)` for another unipotent character `γ'`;
`±γ'`  is called  the Ennola  transform of  `γ`. For  `W` a Weyl group, the
spets  `𝔾` is a finite reductive group, in  which case `z=-1` if `-1` is in
`W`  and `z=1` otherwise.  The function returns  the signed permutation `e`
done by `ennola` on the unipotent degrees (as an `SPerm` of
`1:length(UnipotentCharacters(W))`).

The `SPerm` `e` is not uniquely determined by the degrees since two of them
may  be equal, but is uniquely determined by some additional axioms that we
do not recall here (they include a description of the Ennola-permutation in
terms of the fusion algebras attached to each Lusztig family).

If  a second argument  `z1` is given,  it should be  a power of `z` and the
corresponding power of `e` is returned.

```julia-repl
julia> ennola(rootdatum("3D4"))
SPerm{Int64}: (3,-4)(5,-5)(6,-6)(7,-8)

julia> ennola(complex_reflection_group(14))
SPerm{Int64}: (2,43,-14,16,41,34)(3,35,40,18,-11,42)(4,-37,25,-17,-26,-36)(5,-6,-79)(7,-7)(8,-74)(9,-73)(10,-52,13,31,-50,29)(12,53,15,32,-51,-30)(19,71,70,21,67,68,20,69,72)(22,-39,27,-33,-28,-38)(23,24,-66,-23,-24,66)(44,46,49,-44,-46,-49)(45,48,47,-45,-48,-47)(54,-63,-55,-57,62,-56)(58,-65,-59,-61,64,-60)(75,-77)(76,-76)(78,-78)

```
The  above example shows  that it may  happen that the  order of `z`-Ennola
(here 18) is greater than the order of `z` (here 6); this is related to the
presence  of irrationalities  `q⅓` in  the character  table of the spetsial
Hecke algebra of `W`.

For a non-irreducible group, `z`-ennola is defined if `z` can be considered
an element of the the centre of each irreducible component.
"""
function ennola(t::TypeIrred)
  res=getchev(t,:Ennola)
  if res!==nothing return res end
  InfoChevie("#using SpetsEnnola\n")
  SpetsEnnola(t)
end

function ennola(W)
  get!(W,:ennola)do
    tt=refltype(W)
    if isempty(tt) return SPerm() end
    if length(tt)==1 return ennola(tt[1]) end
    zz=map(tt)do t
      if haskey(t,:orbit) gcd(first.(degrees(t))) else gcd(degrees(t)) end
    end
    z=gcd(zz)
    ee=map((t,z1)->SpetsEnnola(t)^div(z1,z),tt,zz)
    l=cartesian(map(e->1:length(e.d),ee)...)
    p=Perm(l,map(x->map((y,e)->y^Perm(e),x,ee),l))
    sgn=map(x->prod(map((y,e)->sign(y^e),x,ee)),l)
    SPerm(p.d .*sgn)
  end
end

function ennola(W,xi)
  z=gcd(degrees(W))
  if !isone(xi^z) error(xi,"^$z!=1") end
  o=Int(Root1(xi).r*z)
  ennola(W)^o
end

# s is a vector of tuples. Return a vector of vectors E_1,...,E_n such that
# s=Tuple.(cartesian(E_1,...,E_n))
# Assumes all E_i but possibly one are of size 2
function factorset(s)
  if isempty(s[1]) return s[1] end
  for i in s[1]
    s1=filter(y->i in y,s)
    if 2length(s1)!=length(s) continue end
    i=intersect(s1...)
    r=setdiff(s[1], i)
    s2=filter(y->isempty(intersect(y,i)),s)
    a=findfirst(y->issubset(r,y),s2)
    if 2length(s2)!=length(s) || isnothing(a) continue end
    j=setdiff(s2[a],r)
    s1=sort(unique(map(x->setdiff(x, i), s1)))
    s2=sort(unique(map(x->setdiff(x, j), s2)))
    if length(i)==1 i=i[1] end
    if length(j)==1 j=j[1] end
    if s1==s2 return vcat([[i, j]], factorset(s1)) end
  end
  if all(x->length(x)==1,s) 
       [vcat(s...)]
  else [s]
  end
end

# l is a matrix and S a list length size(l,1).
# Find subsets P of axes(l,1) such that sum(l[P,:];dims=1)=S.
# In  addition, lv is a vector of  length size(l,1), v is a sub-multiset of
# lv and P should satisfy tally(lv[P])=tally(v).
function SubsetsSum(S, l::AbstractMatrix, v::Vector, lv::Vector;lim=10)
# println("S=$S;l=$l;v=$v;lv=$lv")
  function sievev(good::Vector{Int}, v)
    v=copy(v)
    for i in good
      p=findfirst(==(lv[i]),v)
      if isnothing(p) return empty(v) end
      deleteat!(v,p)
    end
    v
  end
  c=0
  found=0
  # s assumed to be a subset of sought P at this stage
  # S= initial S minus sum(l[s,:];dims=1)
  # t= remaining elements of axes(l,1) which could be in P
  # nonsolved= indices of nonsolved entries of S
  # v= remaining v to match
  function inner(S, s::AbstractVector{Int}, t::AbstractVector{Int}, 
    nonsolved::AbstractVector{Int}, v::Vector, factor)
    local bad, good, p, sols, res, i, sol, f, ll, solved
#   println("#solved=",size(l,1)-length(nonsolved)," ",join(s),"=>",join(t),
#         " v=",join(map(x->join(x,":"),tally(v))," "))
    c+=1
    if c%1000==0
      InfoChevie("# ",factor,": xcols:",length(nonsolved)-1," xrows:",length(t),
                 " found:", found, "\n")
    end
    t=filter(i->lv[i] in v,t)
    if isempty(t)
      if iszero(S) found+=1
         return [Int[]]
      else return Vector{Int}[]
      end
    end
    ll=map(i->(pos=i,cand=filter(j->l[j,i]!=0,t),sols=Vector{Int}[]),nonsolved)
    sort!(ll,by=x->length(x.cand))
    if length(ll[1].cand)>lim ll=[ll[1]]
    else ll=filter(x->length(x.cand)<=lim,ll)
    end
    solved=Int[]
    good=Int[]
    bad=Int[]
    for p in ll
      append!(p.sols,filter(e->sum(l[e,p.pos])==S[p.pos],combinations(p.cand)))
      if length(p.sols)==0 return Vector{Int}[]
      elseif length(p.sols)==1 push!(solved,p.pos)
      end
      if sum(length,p.sols)>0
       good=union_sorted(good,intersect(p.sols...)) #lines part of any solution
        bad=union_sorted(bad,setdiff(p.cand, union(p.sols)))#part of no solution
      else
        bad=union_sorted(bad,p.cand) #part of no solution
      end
    end
    nonsolved=setdiff(nonsolved, solved)
    if length(good)+length(bad)>0 #progress
      return map(r->vcat(good, r), inner(S-vec(sum(l[good,:];dims=1)), 
      union_sorted(s,good),setdiff(t,union_sorted(good,bad)),nonsolved,sievev(good,v),factor))
    end
    res=Vector{Int}[]
    p=argmax(x->length(x.cand),ll)
    f=length(p.sols)
    nonsolved=setdiff(nonsolved, [p.pos])
    InfoChevie("# ", factor, ": xcols:", length(nonsolved), " xrows:",
       length(t)," in comb(",length(p.cand),")==>", length(p.sols), "\n")
    for sol in p.sols
      good=sol
      bad=setdiff(p.cand, sol)
      if isempty(intersect(good, bad))
         append!(res, map(r->vcat(good,r),
       inner(S-vec(sum(l[good,:];dims=1)),union_sorted(s, good),
        setdiff(t, union_sorted(good, bad)), nonsolved, sievev(good, v),
        string(factor, ":", f))))
      end
      f-=1
    end
    return res
  end
  inner(S, Int[], axes(l,1), eachindex(S), v, "")
end

positive(p::CycPol)=all(>(0),values(p.v))

"""
`FitParameter(sch,m)` given:
sch: schur elements for H(ℤ/e) given as CycPols
m:   a list of length e of Rationals, exponents of params

finds all permutations σ of 1:e such that the parameters pₖ=E(e,σₖ-1)q^mₖ
gives  sch by the formula schᵢ=∏_{j≠i}(1-pᵢ/pⱼ). Since multiplying the pₖ
by  a  scalar  leaves  invariant  the  sch,  σ  is  known  only  modulo a
translation.

The  result is  a list  of pairs  [v1,v2] telling that globally σ(v1)=v2,
where the v1 are sort(collectby(sch,eachindex(sch)),by=minimum)
"""
function FitParameter(sch, m::AbstractVector{<:Rational{<:Integer}})
  den=lcm(denominator.(m))
  e=length(m)
  a=tally(m)
  sch=map(x->subs(x,Pol()^den),sch) # m_k replaced by den*m_k
  function poss(i)
    # for each element (m_k,c_k) of tally(m) p will hold a minimal
    # corresponding possible set of j=σ(i)-σ(l) for l such that m_l=m_k
    p=map(x->collect(0:e-1),a)
    avail=eachindex(a)
    good=Int[]
    term(j,k)=CycPol(1-E(e,j)*Pol([1],Int(den*(m[i]-first(a[k])))))
    v=sch[i]
    while true
      for k in avail
        p[k]=setdiff(p[k],vcat(p[good]...))
        if m[i]!=first(a[k]) p[k]=filter(j->positive(v//term(j,k)),p[k]) end
      end
      good=filter(i->length(p[i])==last(a[i]),avail)
      if isempty(good) return p end
      avail=setdiff(avail, good)
      v//=reduce(*,map(k->prod(j->term(j,k),p[k]),
                       filter(k->m[i]!=first(a[k]),good));init=1)
    end
  end
  pp=poss.(1:e)
  p=findfirst(i->last.(a)==length.(pp[i]),1:e)
  if isnothing(p) error("no solution\n") end
# v=map(i->map(x->sort(unique(mod.(x,e))),i.+pp[p]),0:e-1)
  v=map(0:e-1)do i
    u=map(j->j.+i,pp[p])
    map(x->sort(unique(mod.(x,e))),u)
  end
  v=map(k->filter(i->all(j->issubset(v[i][j],pp[k][j]),1:length(a)),1:e),1:e)
  G=map(x->findfirst(==(x),sch), sch)
  v=map(x->v[x], sort(unique(G)))
  G=map(x->findfirst(==(x),sort(unique(G))), G)
  bad=1:length(v)
  while true
    good=filter(i->length(v[i])==count(j->G[j]==i,1:e),bad)
    bad=setdiff(bad, good)
    for p in bad for a in good v[p]=setdiff(v[p],v[a]) end end
    if isempty(good) break end
  end
  if sum(length,v)>e error("non-unique solution\n") end
  res=collectby(G,eachindex(G))
  pp=vcat(v...)[sortperm(vcat(res...))]
  para=map(i->E(e,pp[i])*Pol([1],Int(den*m[i])),1:e)
  if map(k->prod(j->CycPol(1-exactdiv(para[k],para[j])),filter(j->j!=k,1:e)),1:e)!=sch
    error("schur elms don't match\n")
  end
  map((x,y)->[x,y],res, v)
end

#---------------------- Series -----------------------------------------
# A d-Harish-Chandra series \CE(𝐆,(𝐋,λ) is a record with fields:
#  .d           AsRootOfUnity(ζ) such that L=C_G(V_ζ)
#  .spets       𝐆
#  .levi        𝐋
#  .cuspidal    λ (index in UnipotentCharacters(.levi))
#  .principal   true iff λ=\Id
# optional fields (in s.prop)
#  .charNumbers \CE(𝐆,(𝐋,λ) (indices in UnipotentCharacters(.spets))
#  .degree      The degree of R_𝐋^𝐆(λ)=|𝐆/L|_{q'} degλ
#  .eps         For each χ in .charNumbers sign of <χ,R_𝐋^𝐆(λ)>
#  .dims        χ(1) for γ_χ
#  .WGL         W_𝐆(L,λ)
#  .Hecke       H_𝐆(L,λ)
#  .e           |W_𝐆(L,λ)|
@GapObj struct Series
  spets
  levi
  cuspidal::Int
  d::Root1
  principal::Bool
end

"""
`Series(W, L, cuspidal, d)`

If the reflection coset or group `W` corresponds to the algebraic group `𝐆`
and  `cuspidal`  to  the  `d`-cuspidal  unipotent  character  `λ`  of  `𝐋`,
constructs  the `d`-series corresponding to ``R_𝐋^𝐆(λ)``. The result `s` it
is a record with the following fields and functions:

`s.spets`: the reflection group or coset `W`.

`s.levi`: the subcoset `L`.

`s.cuspidal`: the index of `λ` in `UnipotentCharacters(L)`.

`s.d`: the value of `d` (a `Root1`).

`relative_group(s)`: the group ``W_𝐆(𝐋,λ)``.

`dSeries.RLG(s)`: the `UnipotentCharacter` given by ``R_𝐋^𝐆(λ)``.

`dSeries.eps(s)`:  for each  character `φ`  of `relative_group(s)` the sign
``(-1)^{n_φ}``  in the cohomology  of the variety  defining `RLG(s)` of the
corresponding constituent `γᵩ` of `RLG(s)`.

`degree(s)`: the generic degree of `RLG(s)`, as a `CycPol`.

`charnumbers(s)`:  the indices in  `UnipotentCharacters(W)` of the
constituents of `RLG(s)`.

`hecke(s)`: the hecke algebra ``H_𝐆(𝐋,λ)``.

The function `Series` has another form:

`Series(<W> [,<d> [,<ad>]];k...)`

where  it returns a  vector of `Series`  corresponding to the cuspidal data
described   by  the   arguments  and   the  keywords   (see  the  help  for
`cuspidal_data`).

```julia-repl
julia> W=complex_reflection_group(4)
G₄

julia> Series(W,3;proper=true)
1-element Vector{Series}:
 ζ₃-series R^G₄_{G₄₍₎=Φ₁Φ′₃}(λ==Id)  W_G(L,λ)==Z₆

julia> s=Series(W,3,1)[1]
ζ₃-series R^G₄_{G₄₍₎=Φ₁Φ′₃}(λ==Id)  W_G(L,λ)==Z₆
┌─┬────────────────────────────────────┐
│ │   γᵩ φ(mod 3)  ε parameter family #│
├─┼────────────────────────────────────┤
│1│ φ₁‚₀        1  1      ζ₃q²        1│
│5│ φ₂‚₃       ζ₆  1      -ζ₃q        2│
│2│ φ₁‚₄       ζ₃ -1        ζ₃        4│
│8│ Z₃:2       -1 -1     -ζ₃²q        2│
│9│Z₃:11      ζ₃² -1       ζ₃²        4│
│4│ φ₂‚₅      ζ₆⁵ -1       -ζ₃        4│
└─┴────────────────────────────────────┘

julia> s.spets
G₄

julia> s.levi
G₄₍₎=Φ₁Φ′₃

julia> s.cuspidal
1

julia> s.d
Root1: ζ₃

julia> hecke(s)
hecke(G₆‚₁‚₁,Vector{Mvp{Cyc{Int64}, Int64}}[[ζ₃q², -ζ₃q, ζ₃, -ζ₃²q, ζ₃², -ζ₃]])

julia> degree(s)
ζ₃Φ₁Φ₂²Φ″₃Φ₄Φ₆

julia> dSeries.RLG(s)
[G₄]:<φ₁‚₀>-<φ₁‚₄>-<φ₂‚₅>+<φ₂‚₃>-<Z₃:2>-<Z₃:11>

julia> charnumbers(s)
6-element Vector{Int64}:
 1
 5
 2
 8
 9
 4

julia> dSeries.eps(s)
6-element Vector{Int64}:
  1
  1
 -1
 -1
 -1
 -1

julia> relative_group(s)
G₆‚₁‚₁
```
"""
function Series(WF, levi, cuspidal, d;NC=false)
  if d isa Int && d>0 d=1//d end
  if !(d isa Root1) d=Root1(;r=d) end
  if !(WF isa Spets) WF=spets(WF) end
  principal=UnipotentCharacters(levi).a[cuspidal]==0 &&
            UnipotentCharacters(levi).A[cuspidal]==0
  s=Series(WF,levi,cuspidal,d,principal,Dict{Symbol,Any}())
  if !NC hecke(s) end
  s
end

Series(W,r...;NC=false,k...)=map(x->Series(W,x...;NC),cuspidal_data(W,r...;k...))

function LaurentPolynomials.degree(s::Series)
  get!(s,:degree) do
    deg=exactdiv(conj(generic_sign(s.spets))*generic_order(s.spets,Pol()),
    conj(generic_sign(s.levi))*generic_order(s.levi,Pol()))
    CycPol(shift(deg,-deg.v))*CycPoldegrees(UnipotentCharacters(s.levi))[s.cuspidal]
  end
end

function Base.show(io::IO,s::Series)
  tex=get(io,:TeX,false)
  repl=get(io,:limit,false)
  if !(repl || tex)
    print(io,"Series($(s.spets),$(s.levi),$(s.cuspidal),$(s.d))")
    return
  end
  cname=charnames(io,UnipotentCharacters(s.levi))[s.cuspidal]
  n=repl || tex ? "\\lambda" : "c"
  quad=tex ? "\\quad" : " "
  if s.spets == s.levi
    printTeX(io,s.d, "-cuspidal ",cname," of \$", s.spets,"\$")
  else
    print(io,s.d, "-series ")
    printTeX(io,"\$R^{",s.spets,"}_","{",s.levi,"}($n==",cname,")\$")
  end
  if haskey(s, :WGL)
    if iscyclic(s) && (!haskey(s,:Hecke) || e(s)>4)
      printTeX(io,"\$$quad W_G(L,$n)==Z_{$(e(s))}\$")
    elseif haskey(s, :Hecke)
      printTeX(io,"\$$quad H_G(L,$n)==",TeX(io,hecke(s)),"\$")
    elseif haskey(relative_group(s), :refltype)
      printTeX(io,"\$$quad W_G(L,$n)==",TeX(io,relative_group(s)),"\$")
    else
      printTeX(io,"\$$quad |W_G(L,$n)|==$(length(relative_group(s)))\$")
    end
  elseif haskey(s, :WGLdims)
    printTeX(io,"\$$quad |W_G(L,$n)|==$(sum(WGLdims(s).^2))\$")
  end
# if haskey(s, :translation)
#   Util.print(io,"$quad translation==",s.translation)
# end
  if !haskey(s,:charNumbers) return end
  if get(io,:typeinfo,Any)==typeof(s) return end
end

function Base.show(io::IO, ::MIME"text/plain", s::Series)
  show(io,s)
  if haskey(io,:typeinfo) return end
  if haskey(s,:charNumbers) && !isnothing(charnumbers(s)) format(io,s) end
end

function format(io::IO,s::Series)
  uw = UnipotentCharacters(s.spets)
  e = length(charnumbers(s))
  function f(texn,val)
    push!(col_labels,texn)
    push!(m,val)
  end
  m = []
  col_labels = String[]
  f("\\gamma_\\phi", charnames(io,uw)[charnumbers(s)])
  n="\\varphi"
  if haskey(s,:translation) n*="(mod $(s.translation))" end
  f(n, charnames(io,relative_group(s)))
  f("\\varepsilon", s.eps)
# if haskey(s, :eigen) f("\\hbox{eigen}", s.eigen) end
  if s.cyclic && s.e>1
    if haskey(s, :Hecke) f("\\hbox{parameter}", s.Hecke.para[1])
    elseif haskey(s, :mC) f("\\hbox{mC}",mC(s))
    end
  end
  f("\\hbox{family \\#}", map(j->findfirst(k->j in k.charNumbers,
                                    uw.families), charnumbers(s)))
  if haskey(s, :permutable) && any(x->x!=false,s.permutable)
      f("\\hbox{permutable}", s.permutable)
  end
  m=permutedims(toM(m))
  println(io)
  showtable(io,m;row_labels=charnumbers(s),col_labels)
end

ChevieErr(x...)=xprint("!!!!!!! ",x...)

#  .WGL         W_𝔾 (𝕃,λ) as a relgroup, contains toparent
#  (lifting reflections to elts of W) and reflists (lifting a generator s to
#   reflections of the parabolic 𝕄  of W such that W_𝕄 (𝕃,λ)=<s>)
#  .WGLdims     irr dims of WGL
function Weyl.relative_group(s::Series)
  get!(s,:WGL) do
  WF=s.spets
  W=Group(WF)
  L=Group(s.levi)
  if W==L
    WGL=coxgroup()
    WGL.reflists=Vector{Int}[]
    s.WGLdims=[1]
    s.e=1
    return WGL
  end
  if isone(projector(s)) #central series
    WGL=deepcopy(W)
    WGL.toparent=gens(WGL)
    WGL.reflists=map(x->[x],inclusiongens(W))
    s.WGLdims=CharTable(WGL).irr[:,1]
    s.e=length(WGL)
    return WGL
  end
  if W isa FiniteCoxeterGroup 
    N1=graph_automorphisms(W,inclusiongens(L))
    N=Group(vcat(gens(N1),gens(L)))
  else N=normalizer(W, L)
  end
  if !isone(s.levi.phi)
    if length(L) == 1
      sz=length(conjugacy_classes(WF)[position_class(WF,s.levi.phi)])
      if sz>100000 
        m=Base.get_extension(Chevie,:Gap4)
        if isnothing(m)
          error("*** class too big ($sz>100000)\n",
           "need Gap4.centralizer; do 'using GAP' to make extension available") 
        end
        @time N=m.centralizer(N,s.levi.phi)
      else
        N=centralizer(N,s.levi.phi)
      end
    elseif s.spets isa Cosets.CoxeterCoset
      N=centralizer(N1,s.levi.phi)
      N=Group(vcat(gens(N),gens(Group(s.levi))))
    else #   N:=Stabilizer(N,s.levi); # is shorter but slower...
      NF=centralizer(N, s.levi.phi)
      if length(NF)==length(L)*prod(relative_degrees(s.spets, s.d)) N=NF
      else N=centralizer(N,s.levi.phi,function(p,g)
                  w=reduced(s.levi.W,p^g)
                  w isa Perm ? w : w.phi
                 end)
      if false
        NF=Group(vcat(gens(N),[s.levi.phi]))
        NFQ=NF/L
        N=N/L
        NF=centralizer(NFQ, s.levi) # image of s.levi.phi in NFQ
        NFQ=intersect(NF, N)
        N=Group(vcat(gens(L),map(x->x.phi,gens(NFQ))))
      end
      end
    end
  end
  eig=eigen(UnipotentCharacters(s.levi))
  inds=filter(i->eig[i]==eig[s.cuspidal],eachindex(eig))
  if length(inds)>1
    ud=CycPoldegrees(UnipotentCharacters(s.levi))
    inds=filter(i->ud[i]==ud[s.cuspidal],inds)
    if length(inds)>1
      c=length(N)
      N=centralizer(N,s.cuspidal,(c,g)->c^on_unipotents(s.levi, g))
      if c!=length(N)
        ChevieErr("# WGL:",length(N)//length(L),"/",c//length(L)," fix c\n")
      end
    end
  end
  refs=filter(isdistinguished,reflections(W))
  if length(L) == 1
    WGL=N
  elseif length(N)==length(L)
    if s.levi!=s.spets
      error(s, " N==L\n")
      return nothing
    end
    WGL=coxgroup()
    WGL.reflists=[]
    s.WGLdims=[1]
    s.:e=1
    return WGL
  elseif W isa CoxeterGroup
    WGL=N/L
  else
    WGL=N/Group(gens(L)) # while problem with G333 not solved
  end
  V=lnullspace(projector(s)-I)#The E(d)-eigenspace
  m=vcat(V,lnullspace(projector(s)))^-1
  # V∩ fix(r) for r a reflection in W. The result is a rowspace
  hplane(r)=intersect_rowspace(hyperplane(r),V)
  # hplane of refl. matrix on V to itself in global space
  function Vtoglobal(x)
    h=lnullspace(x-x^0)
    rowspace(hcat(h, fill(0,size(h,1),rank(W)-size(V,1)))*m^-1)
  end
  # restriction of matrix x to V
  restrV(x)=(inv(m)*x*m)[axes(V,1),axes(V,1)]
  # for   reflections of W with same Hplane in V return Dict
  # :hom  element of W above reflection of WGL
  # :refs generators of parabolic of W fixing Hplane
  # :WH   cyclic subgroup of WGL fixing Hplane
  # :eigenvalue of refl. of WGL
  # :root       of refl. of WGL
  # :coroot     of refl. of WGL
  function getreflection(rr)
    rH = hplane(rr[1]) # Hyperplane of WGL
    if length(rH)==0 # cyclic WGL
      H=W
    else
      H=reflection_subgroup(W,restriction(W,vcat(inclusion(W,
                               map(x->x.rootno,rr)),inclusiongens(L))))
    end
    res=Dict{Symbol, Any}(:refs=>inclusiongens(H))
    if H isa CoxeterGroup H=H.G end
    H=intersect(H, N)
    if length(L)!=1 H=H/L end
  # if length(L)!=1 H=H/Group(gens(L)) end # until problem in G333
    if length(H)==1 error("H=L") end
    gen=abelian_gens(H)
    if length(gen)>1 error("H not cyclic |H|=",length(H),order.(gen)) end
    res[:hom]=only(gen)
    r=restrV(reflrep(W, length(L)==1 ? res[:hom] : res[:hom].phi))
    ref=asreflection(improve_type(r))
    n=ref.eig
    n=invmod(exponent(n), order(n))
    r^=n # distinguished reflection
    res[:hom]^=n
    merge!(res,pairs(asreflection(improve_type(r))))
    res[:WH] = H
    res
  end
  #hplanes hashable!sortable
  rrefs=collect(values(groupby(hplane,refs)))
  rrefs=filter(x->!(any(y->inclusion(W,y.rootno) in inclusion(L),x)),rrefs)
  sort!(rrefs,by=x->first(x).rootno)
  reflist=[]
  for r in rrefs
    push!(reflist, getreflection(r))
    if length(Group(map(x->x[:hom],reflist)))==length(WGL)
      rr=improve_type(map(x->x[:root],reflist))
      crr=improve_type(map(x->x[:coroot],reflist))
      WGL=PRG(rr,crr)
      reflist=Vtoglobal.(reflrep(WGL))
      reflist=map(h->rrefs[findfirst(rr->hplane(rr[1])==h,rrefs)], reflist)
      WGL.reflists = map(getreflection, reflist)
      WGL.toparent = map(x->x[:hom], WGL.reflists)
      WGL.reflists = map(x->x[:refs], WGL.reflists)
      s.WGLdims=CharTable(WGL).irr[:,1]
      s.e=length(WGL)
      return WGL
    end
  end
  error("b")
end
end

# Degree in q of the parameters (normalized so the smallest is 0)
function mC(s::Series)
  get!(s,:mC) do
  e = hyperplane_orbits(relative_group(s))
  if length(e)!=1 return else e = e[1].order end
  uc = UnipotentCharacters(s.spets)
  cn = charnumbers(s)[filter(i->s.dims[i]==1,1:length(s.dims))]
  aA = uc.:a[cn] + uc.A[cn]
  lpi(W)=nref(W)+nhyp(W)
  if s.principal
    if minimum(aA)!=0 error("id not in RLG(1)") end
    pG=lpi(Group(s.spets))
    pL=lpi(Group(s.levi))
    D0=pG-pL
    xiL = Root1(PhiOnDiscriminant(s.levi))^order(s.d)
    xiG = Root1(PhiOnDiscriminant(s.spets))^order(s.d)
    if xiL != xiG
      ChevieErr("fixing dimension of variety by xiL-xiG==", xiL - xiG, "\n")
      D0 = (D0 + xiL) - xiG
    end
    # Id in chars
    if !isinteger(D0*s.d.r)
      ChevieErr(s, "\n    ==> (l(pi_G)==", pG, ")-(l(pi_L)==", pL,
                ")+(xiL*d==", xiL, ")-(xiG*d==", xiG, ")  !  ==0( %
         d==", s.d, ") Max(aA)-Min(aA)==", maximum(aA) - minimum(aA), "\n")
    end
    if !isinteger(D0//e)
      ChevieErr("fixing dimension of variety to be divisible by e==", e, "\n")
      D0-=mod(D0, e)
    end
    (D0.-aA)//e
  else
    # (JM+Gunter 18-3-2004) in any case we normalize so that
    # the smallest mC is 0 since the above choice gives unpleasant
    # results for G27
    D0=maximum(aA)-minimum(aA) # just so that mC are positive
    (D0+minimum(aA).-aA)//e
  end
  end
end

#  .projector   .element-equivariant projector on V_ζ
function projector(s)
  get!(s,:projector)do
    ad=count(x->x==(1,Cyc(s.d)), degrees(s.levi))
    q=minimum(map(x->count(==(s.d),x),refleigen(s.levi)))
    if q!=ad error("bad start") end
    q=findall(x->count(==(s.d),x)==q,refleigen(s.levi))
# found an element w∈ Levi such that V_ζ=ζ-eigenspace of w
# .projector= w-equivariant projector on V_ζ
    s.ad=ad
    eigenspace_projector(s.spets,classreps(s.levi)[q[1]],s.d)
  end
end

Groups.iscyclic(s::Series)=get!(()->iscyclic(relative_group(s)),s,:cyclic)

WGLdims(s::Series)=getp(relative_group,s,:WGLdims)

e(s::Series)=getp(relative_group,s,:e)

function RLG(s::Series)
  get!(s,:RLG) do
  RLG=lusztig_induce(s.spets, unichar(s.levi, s.cuspidal))
  if isnothing(RLG) && isone(s.levi.phi)
    RLG=Uch.hc_induce(s.spets, unichar(s.levi,s.cuspidal))
  end
  if isnothing(RLG) ChevieErr(s, ":RLG failed\n")
  elseif degree(s)!=CycPol(degree(RLG))
    ChevieErr(s,":Deg RLG!=Sum(eps[i]*ud[i])\n")
  end
  RLG
  end
end

function canfromdeg(s::Series)
  ad=count(!isone,relative_degrees(s.levi, s.d))
  ud=CycPoldegrees(UnipotentCharacters(s.spets))
  cand=filter(i->ad==valuation(ud[i],s.d),1:length(ud))
# now use that S_φ(q)=ε_φ Deg(RLG(λ))/Deg(γ_φ)
  cand=map(c->Dict(:charNumbers=>c,:sch=>degree(s)//ud[c]),cand)
  cand=filter(c->positive(c[:sch]),cand)
  q=Mvp(:q)
  ad=CycPol(Pol()-Cyc(s.d))^ad
  f=degree(s)//ad
  if !positive(f)
    ChevieErr(s, " cuspidal is not\n")
    return nothing
  end
  v=f(Cyc(s.d))
  filter(cand)do c
    c[:span]=degree(c[:sch])-valuation(c[:sch])
    f=sum(x->x^2,WGLdims(s))*(ud[c[:charNumbers]]//ad)(Cyc(s.d))//v
    if !isinteger(f) return false end
    f=Int(f)
    if !(abs(f) in WGLdims(s)) return false end
    c[:dims]=abs(f)
    c[:eps] =sign(f)
    true
  end
end

# takes a d-series s with s.spets split; fills in s.charNumbers, s.eps, s.dims
function Chars.charnumbers(s::Series)
  get!(s,:charNumbers) do
  if s.levi==s.spets
    s.eps=[1]
    s.dims=[1]
    return [s.cuspidal]
  end
  if !haskey(s,:WGL) && isnothing(relative_group(s)) return nothing end
  rlg=RLG(s)
  uc=UnipotentCharacters(s.spets)
  if rlg!==nothing
    charNumbers=findall(!iszero,rlg.v)
    s.eps=map(sign,rlg.v[charNumbers])
    s.dims=map(abs,rlg.v[charNumbers])
    s.span=degree(degree(s))-valuation(degree(s)).+uc.a.-uc.A
    return charNumbers
  end
  cand=canfromdeg(s)
  ud=CycPoldegrees(uc)
  eig=eigen(UnipotentCharacters(s.levi))[s.cuspidal]
  eig*=map(i->E(order(s.d)^2,i),1:order(s.d)^2)
  cand=filter(c->eigen(UnipotentCharacters(s.spets))[c[:charNumbers]] in eig,cand)
  if length(cand)<length(WGLdims(s))
    ChevieErr(s,": not enough left with predicted eigenvalues in ",Root1.(eig),"\n")
    return nothing
  end
  sort!(cand,by=x->x[:dims])
  sort!(WGLdims(s))
  foo(n)=getindex.(cand,n)
  function foo(n,l)
    if l isa Vector return map(x->x[n], cand[l])
    else return cand[l][n]
    end
  end
  check=function ()
    for n in [:charNumbers, :eps, :dims, :span] setproperty!(s,n,foo(n)) end
    if !isnothing(RLG(s)) && (filter(i->RLG(s).v[i]!=0,eachindex(RLG(s).v))!=
              sort(s.charNumbers) || RLG(s).v[s.charNumbers]!=s.dims.*s.eps)
      ChevieErr(s, ":RLG does not match")
    end
    s.charNumbers
  end
  if length(cand)==length(WGLdims(s)) return check() end
  ud=foo(:dims).*CycPoldegrees(uc)[foo(:charNumbers)].*foo(:eps)
  t=maximum(degree.(ud))
  function c(p)
    pp=p(Pol())
    vcat(fill(0,pp.v),pp.c,fill(0,max(0,t-degree(p))))
  end
  v = SubsetsSum(improve_type(c(degree(s))), improve_type(toM(map(c, ud))),
                 improve_type(WGLdims(s)), improve_type(foo(:dims)))
  InfoChevie("# ", length(v), " found\n")
  if length(v)>10000
    InfoChevie("# ", length(v), " combinations sum to dimRLG\n")
  elseif length(v) == 0
    ChevieErr(s, " no combination sums to dimRLG\n")
    return nothing
  end
  if iscyclic(s)
    s.charNumbers=foo(:charNumbers)
    s.dims=foo(:dims)
    mC(s)
    if length(v)>1
      v=filter(a->foo(:span,a)==map(i->
           sum(j->abs(mC(s)[i]-mC(s)[j]),setdiff(a,[i])),a),v)
      if length(v)>10000
          InfoChevie("# ", length(v), " combinations have right span\n")
      end
    elseif length(v) == 0
      ChevieErr(s, " no combination has right span\n")
      return nothing
    end
    delete!.(Ref(s.prop),[:charNumbers,:dims,:mC])
  end
  if length(v)>1
    InfoChevie("# after span ", length(v), " combinations\n")
    InfoChevie("# Warning: using Mackey with tori for ", s, "\n")
    i=fusion_conjugacy_classes(s.levi, s.spets)
    c=CharTable(s.spets).centralizers[i].//CharTable(s.levi).centralizers
    t=dlCharTable(s.levi)[:,s.cuspidal].*c
    t=map(k->sum(t[filter(j->i[j]==k,eachindex(i))]),
                       1:nconjugacy_classes(s.spets))
    c=dlCharTable(s.spets)
    v=filter(a->c[:,foo(:charNumbers,a)]*(foo(:dims,a).*foo(:eps,a))==t,v)
  end
  if length(v)>1
    ChevieErr(s," ",join(factorset(map(x->foo(:charNumbers,x),v)),"x"),
              " chars candidates: using RLG\n")
    if isnothing(RLG(s)) return nothing end
    v=filter(l->all(i->RLG(s).v[foo(:charNumbers,i)]!=0,l),v)
    if length(v)!=1 error() end
  elseif length(v)==0
    ChevieErr(s," no candidates left\n")
    return nothing
  end
  cand=cand[v[1]]
  check()
  end
end

eps(s::Series)=getp(charnumbers,s,:eps)

COMPACTCOHOMOLOGY=true
function paramcyclic(s::Series)
  if !iscyclic(s) error("fill assumes relative_group(s) cyclic\n") end
  if isnothing(charnumbers(s)) return nothing end
  if e(s)==1 return (s.Hecke=hecke(relative_group(s))) end
  uc=UnipotentCharacters(s.spets)
  Schur=CycPoldegrees(uc)[charnumbers(s)]
  Schur=map(x->degree(s)//Schur[x]*eps(s)[x],1:e(s))
  s.eigen=eigen(uc)[charnumbers(s)]
  LFrob=eigen(UnipotentCharacters(s.levi))[s.cuspidal]
  m=degrees(Group(s.spets))
  s.delta=lcm(map(x->order(Root1(x[2])),filter(x->x[1]!=1,degrees(s.spets))))
  rr(j,i)=(i-1)//e(s)-mC(s)[j]*s.d.r
  param(j,i)=Mvp(:q)^mC(s)[j]*Root1(;r=rr(j,i))
  # parameters of Hecke algebra are map(i->param(i,i),1:e(s))
  mmp=FitParameter(subs.(Schur,Pol([s.d],1)),COMPACTCOHOMOLOGY ? mC(s) : -mC(s))
# possible perms encoded by mmp[i][1]->mmp[i][2]
  nid=charnumbers(uc.almostHarishChandra[1])[charinfo(s.spets).positionId]
  function predeigen(j, i)#eigenvalue for mC[j] and ζ_e^{i-1}
    if nid in charnumbers(s)# id should correspond to id(WGL)
      Root1(;r=s.d.r*e(s)*s.delta*
                   (rr(j,i)+s.d.r*mC(s)[findfirst(==(nid),charnumbers(s))]))
    else Root1(;r=s.delta*rr(j,i)*e(s)*s.d.r+LFrob.r)
              # as fraction predeigen_i:=delta di -delta m_i d^2e
    end
  end
  series=map(x->findfirst(y->x in charnumbers(y),uc.harishChandra),charnumbers(s))
  unique=filter(p->length(p[1])==1,mmp)
  ratio=map(p->s.eigen[p[1][1]]//Cyc(predeigen(p[1][1], p[2][1])), unique)
  if length(Set(ratio))>1
    ChevieErr(s, " eigenvalue ratios==", ratio, "\n")
    return nothing
  end
  ratio=Root1(ratio[1]).r
# now find integer translation t such that mod 1. we have t delta d=ratio
  ratio*=order(s.d^s.delta)
  if !isinteger(ratio)
    ChevieErr(s, "non-integral ratio==", ratio, "\n")
    return nothing
  end
  if ratio==0 r=0
  else r=s.d^s.delta
   r=mod(ratio*invmod(exponent(r),order(r)),order(r))
  end
  mmp=map(x->(x[1],map(y->mod(y,e(s))+1,map(j->j.+(r-1),x[2]))), mmp)
  r=fill(0,e(s))
  s.permutable=map(x->0,charnumbers(s))
  j=1
  for i in mmp
    a=arrangements(i[2], length(i[2]))
    p=findall(A->s.eigen[i[1]]==map(j->Cyc(predeigen(i[1][1],j)),A),a)
    if isempty(p)
      ChevieErr(s, "predicted eigenvalues cannot match actual\n")
      return nothing
    else
      if length(p)>1
        o=orbits(Group(map(x->Perm(a[p[1]],x),a[p])), 1:length(i[2]))
        if length(p)!=prod(x->factorial(length(x)),o) error() end
        for u in o
          s.permutable[i[1][u]]=fill(j,length(u))
          j+=1
        end
      end
      r[i[1]] = a[p[1]]
    end
  end
  p=inv(Perm(sortperm(r)))
  for i in [:mC, :charNumbers, :eigen, :span, :eps, :dims, :permutable]
    setproperty!(s,i,invpermute(getproperty(s,i),p))
  end
  s.translation=filter(t->
       s.eigen==map(i->Cyc(predeigen(i,mod(i+t,e(s)))),1:e(s)),0:e(s)-1)
  if length(s.translation)==1 delete!(s.prop,:translation)
  else
    p=s.translation[2]
    if any(x->mod(x, p)!=0,s.translation) || length(s.translation)!=e(s)//p
        error()
    end
    quality=map(t->map(i->param(mod(t+i-1,e(s))+1,i),1:e(s)), s.translation)
    quality=map(x->conductor(collect(values(prod(y->Mvp(:x)-y,x).d))),quality)
    quality=1 .+s.translation[filter(i->quality[i]==minimum(quality),eachindex(quality))]
    if isempty(quality) quality=[1] end
    m=GAPENV.Rotations(mC(s))
    m=findfirst(==(maximum(m[quality])),m)
    m=circshift(1:e(s),m-1)
    for i in [:mC,:charNumbers,:eigen,:span,:eps,:dims,:permutable]
      setproperty!(s,i,getproperty(s,i)[m])
    end
    s.translation=p
  end
  s.Hecke=hecke(relative_group(s), improve_type([map(i->param(i,i),1:e(s))]))
  s.Hecke
end

# get Hecke of series s. If scalar-twisted untwist first
function getHecke(s)
# println("s=",s.spets,"\n")
  t=refltype(s.spets)
  if length(t)==1 && haskey(t[1],:scalar) && !all(isone,t[1].scalar)
    t=copy(t[1])
    scal=t.scalar
    InfoChevie("# removing scal==", scal, "\n")
    t.scalar=map(x->1,t.scalar)
    if t.orbit[1].series==:B && t.orbit[1].rank==2 && t.twist==perm"(1,2)"
      t.orbit[1].cartanType=E(8)-E(8, 3)
    end
    g=reflection_group(t)
    if isnothing(g) return nothing end
    W=Group(s.spets)
    p=Group(g)(word(W,s.levi.phi/s.spets.phi)...)
    inds=convert(Vector{Int},map(x->findfirst(==(Group(g)(word(W,x)...)),
                   refls(Group(g))),gens(Group(s.levi))))
    l=subspets(g,inds,p)
    if length(scal)>1
      ChevieErr("scal==", scal, " unimplemented\n")
      return nothing
    end
    scal=Root1(scal[1])
# one should do an Ennola of the Levi's cuspidal by the absorbed part
# of scal by the center of the levi
    if !(s.cuspidal in cuspidal(UnipotentCharacters(l),s.d/scal))
      e = ennola(Group(l))
      c = abs(s.cuspidal^e)
      if !(c in cuspidal(UnipotentCharacters(l), s.d/scal))
        c=abs(s.cuspidal^(e^-1))
      end
    else c=s.cuspidal
    end
    s1=Series(g, l, c, (s.d/scal).r;NC=true)
    p=getHecke(s1)
    if !isnothing(p)
      p=map(x->x(q=Cyc(inv(scal))*Mvp(:q)),p)
    end
    return p
  else
    if any(u->haskey(u, :scalar) && !all(x->x==1,u.scalar), t)
      ChevieErr("scals==", map(u->u.scalar,t), " unimplemented\n")
#     return nothing
    end
    paramcyclic(s)
    haskey(s,:Hecke) ? s.Hecke.para[1] : nothing
  end
end

#  .classno     class of s.levi with ζ-eigenspace V_ζ
#  .element     representative of classno
function RelativeSeries(s)
# if !(haskey(s, :charNumbers)) charnumbers(s) end
# find simplest regular eigenvalue q of s.levi
  eig=union(map(x->x isa Int ? prime_residues(x).//x : [x],
                map(x->x.r,regular_eigenvalues(s.levi)))...)
  c=minimum(denominator.(eig))
  c=minimum(filter(x->denominator(x)==c,eig))
  c=Root1(;r=c)
  eig=refleigen(s.levi)
  q=maximum(map(x->count(==(c),x), eig))
  s.classno=findall(x->count(==(c),x)==q,eig)
  if length(s.classno)>1 error("classno==",s.classno) end
  s.element=conjugacy_classes(s.levi)[s.classno[1]].representative
  WGL=relative_group(s)
  res=map(enumerate(WGL.reflists)) do (i,r)
 # println("r=",r,"\nphi=",s.element/s.spets.phi)
   R=subspets(s.spets,restriction(s.spets,r),s.element/s.spets.phi)
    l=inclusion(s.levi)
    if !issubset(l,inclusion(R))
      r=map(r)do w
        rr=roots(Group(s.spets))
        p=findfirst(y->!isnothing(ratio(rr[w],rr[y])),l)
        !isnothing(p) ? l[p] :  w
      end
      R=subspets(s.spets,r,s.element/s.spets.phi)
      if !issubset(r,inclusion(R))
        l=sort(unique(map(l)do w
          rr=roots(Group(s.spets))
          p=findfirst(y->!isnothing(ratio(rr[w],rr[y])),inclusion(Group(R)))
          !isnothing(p) ? inclusion(R,p) : w
        end))
        if !issubset(l,inclusion(R))
          error("could not change r==",r," so that Subspets(r) contains",l,"\n")
        end
      end
    end
    #   println("R=",R," L=",subspets(R,l,s.element/R.phi))
    #   println("R=",R," l=",l," s.element/R.phi=",s.element/R.phi)
    p=Series(R,subspets(R,restriction(R,l),s.element/R.phi),s.cuspidal,s.d.r;NC=true)
    p.orbit=simple_reps(WGL)[i]
    r=reflections(WGL)[1:nhyp(WGL)]
    if gens(WGL)==WGL.toparent
      r=filter(x->Perm(x) in Group(p.spets),r)
    else
      w=map(x->prod(WGL.toparent[word(x)]), r)
#     xprint(w,"\n")
      r=r[map(w)do x
          g=Group(p.spets)
          if !(x isa Perm) return x.phi in g
          else return x in g
          end
         end]
    end
    p.WGL=reflection_subgroup(WGL, map(x->x.rootno,r))
    p.WGLdims=CharTable(relative_group(p)).irr[:,1]
    p.e=length(p.WGL)
    return p
  end
  s.relativeSpets=map(x->x.spets, res)
  p=getHecke.(res)
  if nothing in p return res end
  s.Hecke=hecke(WGL,improve_type(p))
# xprint("Hecke=",s.Hecke,"\n")
  u1=schur_elements(s.Hecke)//1
  if any(x->any(y->!all(isinteger, values(y.d)), keys(x.d)),u1)
    ChevieErr(s.Hecke, " fractional: wrong set of SchurElements")
    return res
  end
  u1=map(x->degree(s)//CycPol(x), u1)
  degcusp=CycPoldegrees(UnipotentCharacters(s.levi))[s.cuspidal]
  ud=map(x->x*sign(Int((x//degcusp)(Cyc(s.d)))),
         CycPoldegrees(UnipotentCharacters(s.spets))[charnumbers(s)])
  p=Perm(u1,ud)
# the permutation should also take in account eigenvalues
  if isnothing(p)
    xprint("u1=",u1," ud=",ud,"\n")
    println(indexin(ud,u1))
    println(indexin(u1,ud))
    ChevieErr(s.Hecke, " wrong set of SchurElements")
    return res
  end
  s.charNumbers=invpermute(s.charNumbers,p)
  if haskey(s,:span) s.span=invpermute(s.span,p) end
  aA=map(x->E(order(s.d)^2,valuation(x)+degree(x)),u1)
  p=position_regular_class(WGL,s.d)
  if p == false && length(WGL)==1 p=1 end
  o=map(x->x[p]//x[1], eachrow(CharTable(WGL).irr))
  s.predictedEigen=map(i->aA[i]*o[i], 1:length(o)) *
    eigen(UnipotentCharacters(s.levi))[s.cuspidal]
  return res
end

CHEVIE[:relativeSeries]=true
function HeckeAlgebras.hecke(s::Series)
  if haskey(s,:Hecke) return s.Hecke end
  if iscyclic(s) paramcyclic(s)
  elseif CHEVIE[:relativeSeries]
    InfoChevie("# Relative: ",s,"\n")
    RelativeSeries(s)
  end
  if haskey(s,:Hecke) return s.Hecke
  else return nothing
  end
end
end
