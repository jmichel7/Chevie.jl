"""
This  module gives information  about the unipotent  conjugacy classes of a
connected  reductive  group  over  an  algebraically  closed field `k`, and
various  invariants attached to  them. The unipotent  classes depend on the
characteristic of `k`; their classification differs when the characteristic
is  not *good*  (that is,  when it  divides one  of the coefficients of the
highest  root).  In  good  characteristic,  the  unipotent  classes  are in
bijection with nilpotent orbits on the Lie algebra.

We  give  the  following  information  for  a unipotent element `u` of each
class:

- the centralizer ``C_𝐆 (u)``, that we describe by the reductive part of
  ``C_𝐆  (u)^0``, by the  group of components  ``A(u):=C_𝐆 (u)/C_𝐆 (u)^0``,
  and by the dimension of its radical.

- in good characteristic, the  Dynkin-Richardson  diagram.

- the Springer correspondence,  attaching characters of  the Weyl group or
  relative Weyl groups to each character of `A(u)`.

The  Dynkin-Richarson diagram is attached to a nilpotent element `e` of the
Lie  algebra `𝔤`.  By the  Jacobson-Morozov theorem  there exists an `𝔰𝔩₂`
subalgebra of `𝔤` containing `e` as the element ``\\begin{pmatrix}1&0\\\\0&1
\\end{pmatrix}``. Let `𝐒` be the torus ``\\begin{pmatrix}h&0\\\\0&h^{-1}
\\end{pmatrix}`` of `SL₂` and let `𝐓` be a
maximal  torus containing `𝐒`, so that `𝐒`  is the image of a one-parameter
subgroup  `σ∈ Y(𝐓)`. Consider the root decomposition ``𝔤=∑_{α∈Σ}𝔤_α`` given
by  `𝐓`; then  `α↦⟨σ,α⟩` defines  a linear  form on  `Σ`, determined by its
value  on simple roots. It  is possible to choose  a system of simple roots
`Π`  so that `⟨σ,α⟩≥ 0` for `α∈Π`,  and then `⟨σ,α⟩∈{0,1,2}` for any `α∈Π`.
The  Dynkin diagram of `Π` decorated by  these values `0,1,2` is called the
Dynkin-Richardson  diagram of `e`, and in good characteristic is a complete
invariant of its `𝔤`-orbit.

Let  `𝓑`  be  the  variety  of  all  Borel  subgroups  and  let `𝓑ᵤ` be the
subvariety  of Borel subgroups  containing the unipotent  element `u`. Then
`dim C_𝐆(u)=rank 𝐆 + 2 dim 𝓑_u` and in good characteristic this dimension
can  be computed from  the Dynkin-Richardson diagram:  the dimension of the
class of `u` is the number of roots `α` such that `⟨σ,α⟩∉{0,1}`.

We   describe  now  the  Springer  correspondence.  Indecomposable  locally
constant  `𝐆`-equivariant  sheaves  on  `C`,  called  *local  systems*, are
parameterized  by irreducible characters of `A(u)`. The *ordinary* Springer
correspondence  is a bijection  between irreducible characters  of the Weyl
group  and a large subset  of the local systems  which contains all trivial
local  systems (those parameterized by the  trivial character of `A(u)` for
each  `u`).  More  generally,  the  *generalized*  Springer  correspondence
associates  to each local  system a (unique  up to `𝐆`-conjugacy) *cuspidal
pair*  of a Levi  subgroup `𝐋` of  `𝐆` and a  `cuspidal` local system on an
unipotent  class of `𝐋`, such that the set of local systems associated to a
given cuspidal pair is parameterized by the characters of the relative Weyl
group ``W_𝐆 (𝐋):=N_𝐆 (𝐋)/𝐋``. There are only few cuspidal pairs.

The  Springer correspondence gives information on the character values of a
finite  reductive  groups  as  follows:  assume  that  `k` is the algebraic
closure of a finite field ``𝔽_q`` and that `F` is the Frobenius attached to
an  ``𝔽_q``-structure of `𝐆`. Let `C`  be an `F`-stable unipotent class and
let  ``u∈ C^F``;  we  call  `C`  the  *geometric  class*  of  `u`  and the
``𝐆^F``-classes  inside  ``C^F``  are  parameterized  by  the `F`-conjugacy
classes  of `A(u)`, denoted `H¹(F,A(u))` (most of  the time we can find `u`
such  that `F` acts trivially  on `A(u)` and `H¹(F,A(u))`  is then just the
conjugacy  classes). To an `F`-stable character  `φ` of `A(u)` we associate
the  *characteristic function* of the  corresponding local system (actually
associated to an extension `φ̃` of `φ` to `A(u).F`); it is a class function
``Y_{u,φ}`` on ``𝐆^F`` which can be normalized so that:
``Y_{u,φ}(u₁)=φ̃(cF)``  if `u₁` is  geometrically conjugate to  `u` and its
``𝐆^F``-class  is parameterized by the  `F`-conjugacy class `cF` of `A(u)`,
otherwise ``Y_{u,φ}(u₁)=0``. If the pair `u,φ` corresponds via the Springer
correspondence to the character `χ` of ``W_𝐆(𝐋)``, then ``Y_{u,φ}`` is also
denoted  `Yᵪ`. There  is another  important class  of functions  indexed by
local  systems: to a local system on  class `C` is attached an intersection
cohomology  complex, which is a complex of sheaves supported on the closure
`C̄`.  To  such  a  complex  of  sheaves  is associated its *characteristic
function*,  a class function of ``𝐆^F``  obtained by taking the alternating
trace  of the Frobenius acting on the  stalks of the cohomology sheaves. If
``Y_ψ``   is   the   characteristic   function   of  a  local  system,  the
characteristic   function  of  the  corresponding  intersection  cohomology
complex  is denoted  by ``X_ψ``.  This function  is supported  on `C̄`, and
Lusztig  has shown that ``X_ψ=∑ᵩ P_{ψ,χ} Yᵪ`` where ``P_{ψ,χ}`` are integer
polynomials  in `q` and `Yᵪ` are attached to local systems on classes lying
in `C̄`.

Lusztig   and  Shoji  have  given  an   algorithm  to  compute  the  matrix
``P_{ψ,χ}``,   which  is  implemented  in  Chevie.  The  relationship  with
characters   of  ``𝐆(𝔽_q)``,  taking  to  simplify  the  ordinary  Springer
correspondence,  is that the  restriction to the  unipotent elements of the
almost  character ``R_χ`` is equal to ``q^{bᵪ} Xᵪ``, where `bᵪ` is `dim 𝓑ᵤ`
for  an element `u` of the class `C`  such that the support of `χ` is `C̄`.
The restriction of the Deligne-Lusztig characters ``R_w`` to the unipotents
are  called the *Green functions*  and can also be  computed by Chevie. The
values  of  all  unipotent  characters  on  unipotent  elements can also be
computed  in principle by applying  Lusztig's Fourier transform matrix (see
the  section on the Fourier  matrix) but there is  a difficulty in that the
`Xᵪ` must be first multiplied by some roots of unity which are not known in
all  cases (and when known may depend on the congruence class of `q` modulo
some small primes).

We illustrate these computations on some examples:

```julia-repl
julia> UnipotentClasses(rootdatum(:sl,4))
UnipotentClasses(A₃)
1111<211<22<31<4
   u│D-R dBu B-C          C(u) A₃(A₃₍₎=Φ₁³) A₁(A₃₍₁₃₎=A₁×A₁Φ₁)/-1 .(A₃)/ζ₄
────┼──────────────────────────────────────────────────────────────────────
4   │222   0 222         q³.Z₄          1:4                  -1:2    ζ₄:Id
31  │202   1 22.    q⁴.A₁₍₎=Φ₁        Id:31                               
22  │020   2 2.2      q⁴.A₁.Z₂         2:22                 11:11         
211 │101   3 2.. q⁵.A₂₍₁₎=A₁Φ₁       Id:211                               
1111│000   6 ...            A₃      Id:1111

   u│.(A₃)/-ζ₄
────┼──────────
4   │   -ζ₄:Id
31  │
22  │
211 │
1111│
```

The  first column in the table gives the name of the unipotent class, which
here  is  a  partition  describing  the  Jordan  form. The partial order on
unipotent  classes given by Zariski closure  is given before the table. The
column   'D-R',   displayed   only   in   good  characteristic,  gives  the
Dynkin-Richardson  diagram  for  each  class;  the  column  'dBu' gives the
dimension  of  the  variety  `𝓑ᵤ`.  The  column 'B-C' gives the Bala-Carter
classification  of `u`, that is  in the case of  `sl₄` it displays `u` as a
regular  unipotent  in  a  Levi  subgroup  by  giving the Dynkin-Richardson
diagram  of a regular  unipotent (all 2's)  at entries corresponding to the
Levi  and '.' at  entries which do  not correspond to  the Levi. The column
'C(u)'  describes the  group ``C_𝐆(u)``:  a power  `qᵈ` describes  that the
unipotent  radical  of  ``C_𝐆(u)``  has  dimension  `d` (thus `qᵈ` rational
points);  then follows a  description of the  reductive part of the neutral
component  of ``C_𝐆(u)``,  given by  the name  of its  root datum.  Then if
``C_𝐆(u)``  is  not  connected,  the  description  of `A(u)` is given using
another  vocabulary: a  cyclic group  of order  4 is  given as  'Z4', and a
symmetric group on 3 points would be given as 'S3'.

For  instance, the first class '4'  has ``C_𝐆(u)^0`` unipotent of dimension
`3`  and `A(u)` equal to 'Z4', the cyclic  group of order 4. The class '22'
has  ``C_G(u)`` with unipotent radical of  dimension `4`, reductive part of
type  'A1' and  `A(u)` is  'Z2', that  is the  cyclic group of order 2. The
other  classes have ``C_𝐆(u)`` connected. For  the class '31' the reductive
part of ``C_G(u)`` is a torus of rank 1.

Then  there is one column for each *Springer series*, giving for each class
the pairs 'a:b' where 'a' is the name of the character of `A(u)` describing
the  local system  involved and  'b' is  the name  of the  character of the
(relative)  Weyl group corresponding by the Springer correspondence. At the
top  of the column is  written the name of  the relative Weyl group, and in
brackets  the name  of the  Levi affording  a cuspidal  local system; next,
separated  by a `/` is a description of the central character associated to
the  Springer series  (omitted if  this central  character is trivial): all
local  systems  in  a  given  Springer  series have same restriction to the
center of `𝐆`. To find what the picture becomes for another algebraic group
in  the  same  isogeny  class,  for  instance the adjoint group, one simply
discards the Springer series whose central character becomes trivial on the
center  of `𝐆`; and  each group `A(u)`  has to be  quotiented by the common
kernel  of  the  remaining  characters.  Here  is the table for the adjoint
group:

```julia-repl
julia> UnipotentClasses(coxgroup(:A,3))
UnipotentClasses(A₃)
1111<211<22<31<4
   u│D-R dBu B-C          C(u) A₃(A₃₍₎=Φ₁³)
────┼───────────────────────────────────────
4   │222   0 222            q³         Id:4
31  │202   1 22.    q⁴.A₁₍₎=Φ₁        Id:31
22  │020   2 2.2         q⁴.A₁        Id:22
211 │101   3 2.. q⁵.A₂₍₁₎=A₁Φ₁       Id:211
1111│000   6 ...            A₃      Id:1111
```

Here is another example:

```julia-repl
julia> UnipotentClasses(coxgroup(:G,2))
UnipotentClasses(G₂)
1<A₁<Ã₁<G₂(a₁)<G₂
     u│D-R dBu B-C  C(u)    G₂(G₂₍₎=Φ₁²)  .(G₂)
──────┼─────────────────────────────────────────
G₂    │ 22   0  22    q²         Id:φ₁‚₀       
G₂(a₁)│ 20   1  20 q⁴.S₃ 21:φ′₁‚₃ 3:φ₂‚₁ 111:Id
Ã₁    │ 01   2  .2 q³.A₁         Id:φ₂‚₂       
A₁    │ 10   3  2. q⁵.A₁        Id:φ″₁‚₃       
1     │ 00   6  ..    G₂         Id:φ₁‚₆       
```

which illustrates that on class `G₂(a₁)` there are two local systems in the
principal  series of  the Springer  correspondence, and  a further cuspidal
local system. Also, from the 'B-C' column, we see that that class is not in
a  proper Levi,  in which  case the  Bala-Carter diagram coincides with the
Dynkin-Richardson diagram.

The  characteristics 2 and  3 are not  good for 'G2'.  To get the unipotent
classes  and the Springer correspondence in bad characteristic, one gives a
second argument to the function 'UnipotentClasses':

```julia-repl
julia> UnipotentClasses(coxgroup(:G,2),3)
UnipotentClasses(G₂)
1<A₁,(Ã₁)₃<Ã₁<G₂(a₁)<G₂
     u│dBu B-C  C(u) G₂(G₂₍₎=Φ₁²) .(G₂) .(G₂)  .(G₂)
──────┼──────────────────────────────────────────────
G₂    │  0  22 q².Z₃       1:φ₁‚₀       ζ₃:Id ζ₃²:Id
G₂(a₁)│  1  20 q⁴.Z₂       2:φ₂‚₁ 11:Id             
Ã₁    │  2  .2    q⁶      Id:φ₂‚₂                   
A₁    │  3  2. q⁵.A₁     Id:φ″₁‚₃                   
(Ã₁)₃ │  3  ?? q⁵.A₁     Id:φ′₁‚₃                   
1     │  6  ..    G₂      Id:φ₁‚₆
```

The  function 'ICCTable' gives the  transition matrix between the functions
`Xᵪ`  and `Y_ψ`.

```julia-repl
julia> uc=UnipotentClasses(coxgroup(:G,2));
julia> t=ICCTable(uc)
Coefficients of Xᵪ on Yᵩ for G₂
      │G₂ G₂(a₁)⁽²¹⁾ G₂(a₁) Ã₁ A₁  1
──────┼──────────────────────────────
Xφ₁‚₀ │ 1          0      1  1  1  1
Xφ′₁‚₃│ 0          1      0  1  0 q²
Xφ₂‚₁ │ 0          0      1  1  1 Φ₈
Xφ₂‚₂ │ 0          0      0  1  1 Φ₄
Xφ″₁‚₃│ 0          0      0  0  1  1
Xφ₁‚₆ │ 0          0      0  0  0  1
```

Here  the row labels  and the column  labels show the  two ways of indexing
local  systems: the  row labels  give the  character of the relative Weyl
group and the column labels give the class and the name of the local system
as  a character  of `A(u)`:  for instance,  'G2(a1)' is the trivial local
system  of the  class 'G2(a1)',  while 'G2(a1)(21)'  is the local system on
that class corresponding to the 2-dimensional character of `A(u)=A₂`.
"""
module Ucl

using ..Gapjm

export UnipotentClasses, UnipotentClassOps, UnipotentClassesOps, ICCTable,
 induced_linear_form, special_pieces

struct UnipotentClass
  name::String
  parameter::Any
  dimBu::Int
  prop::Dict{Symbol,Any}
end

function Base.getproperty(c::UnipotentClass, s::Symbol)
  if s===:name return getfield(c, :name)
  elseif s===:parameter return getfield(c, :parameter)
  elseif s===:dimBu return getfield(c, :dimBu)
  elseif s===:prop return getfield(c, :prop)
  else return getfield(c, :prop)[s]
  end
end

Base.setproperty!(c::UnipotentClass,s::Symbol,v)=getfield(c,:prop)[s]=v

Base.haskey(c::UnipotentClass, s::Symbol)=haskey(c.prop,s)

struct UnipotentClasses
  classes::Vector{UnipotentClass}
  p::Int
  orderclasses::Poset
  springerseries::Vector{Dict}
  prop::Dict{Symbol,Any}
end

function nameclass(u::Dict,opt=Dict{Symbol,Any}())
# println("u=$u, opt=$opt")
  if haskey(opt,:mizuno) && haskey(u,:mizuno) n=u[:mizuno]
  elseif haskey(opt,:shoji) && haskey(u,:shoji) n=u[:shoji]
  else n=u[:name]
  end
  TeX=haskey(opt,:TeX)
  n=fromTeX(n;opt...)
  if haskey(opt,:locsys) && opt[:locsys]!=charinfo(u[:Au])[:positionId]
    cl="("*charnames(u[:Au];opt...)[opt[:locsys]]*")"
    n*="^{$cl}"
    n=fromTeX(n;opt...)
  elseif haskey(opt,:class) && opt[:class]!=charinfo(u[:Au])[:positionId]
    cl=classinfo(u[:Au])[:classnames][opt[:class]]
    n=TeX ? "\\mbox{\$$n\$}_{($cl)}" : fromTeX("$(n)_{$cl}";opt...)
  end
  n
end

function name(io::IO,u::UnipotentClass)
  nameclass(merge(u.prop,Dict(:name=>u.name)),io.dict)
end

function Base.show(io::IO,u::UnipotentClass)
  print(io,"UnipotentClass(",name(io,u),")")
end

UnipotentClassOps=Dict{Symbol,Any}(:Name=>nameclass)

"""
`induced_linear_form(W, K, h)`

This routine can be used to find the Richardson-Dynkin diagram of the class
in  the algebraic group `𝐆`  which contains a given  unipotent class of a
reductive subgroup of maximum rank `𝐒` of `𝐆`.

It  takes a linear  form on the  roots of `K`,  defined by its value on the
simple  roots (these values  can define a  Dynkin-Richardson diagram); then
extends  this linear form to the roots of `𝐆` by `0` on the orthogonal of
the  roots of `K`; and finally conjugates  the resulting form by an element
of the Weyl group so that it takes positive values on the simple roots.

```julia-repl
julia> W=coxgroup(:F,4)
F₄

julia> H=reflection_subgroup(W,[1,3])
F₄₍₁₃₎=A₁×Ã₁Φ₁²

julia> Ucl.induced_linear_form(W,H,[2,2])
4-element Array{Int64,1}:
 0
 1
 0
 0

julia> uc=UnipotentClasses(W);

julia> uc.classes[4].prop
Dict{Symbol,Any} with 7 entries:
  :dynkin     => [0, 1, 0, 0]
  :dimred     => 6
  :red        => A₁×A₁
  :Au         => .
  :balacarter => [1, 3]
  :dimunip    => 18
  :AuAction   => A₁×A₁

julia> uc.classes[4]
UnipotentClass(A₁+Ã₁)
```

The  example above shows that the class containing the regular class of the
Levi subgroup of type `A₁× Ã₁` is the class |A1+~A1|.
"""
function induced_linear_form(W,K,h)
# print("W=$W K=$K h=$h");
  if semisimplerank(K)==0 return fill(0,semisimplerank(W)) end
  h=vcat(h,zeros(Int,rank(W)-semisimplerank(K)))
  h=Int.(inv(Rational.(PermRoot.baseX(K.G)))*h)
  r=parent(W).G.roots[inclusion(W)]
  v=toM(r[1:W.N])*h
  w=with_inversions(W,filter(i->v[i]<0,1:W.N))
  map(i->r[i]*h,restriction.(Ref(W),
        inclusion.(Ref(W),eachindex(gens(W))).^(w^-1)))
end

function DistinguishedParabolicSubgroups(W)
  filter(combinations(eachindex(gens(W)))) do J
    if isempty(J) return true end
    p=fill(1,semisimplerank(W))
    p[J]=fill(0,length(J))
    p=toM(W.rootdec[1:W.N])*p
    2*count(iszero,p)+semisimplerank(W)==count(isone,p)
  end
end

function BalaCarterLabels(W)
  vcat(map(parabolic_representatives(W)) do J
    L=reflection_subgroup(W,J)
    map(DistinguishedParabolicSubgroups(L))do D
      w=fill(2,length(J));w[D].=0
      u=copy(J);u[D]=-u[D]
      [induced_linear_form(W,L,w),u]
    end
  end...)
end

# QuotientAu(Au,chars): chars is a list of indices of characters of Au.
# If  k is the common kernel of chars, QuotientAu returns a
# Dict(Au=>Au/k,
#      chars=>index of chars as characters of Au/k,
#      gens=>words in Au preimages of generators of Au/k)
# Since  we have problems with quotient groups, we are forced to program an
# ad  hoc solution which works only  for Au actually occuring for unipotent
# classes of a reductive group G.
function QuotientAu(Au,chars)
  AbGens=function(g)
    res=empty(gens(g))
    l=gens(g)
    while !isempty(l)
      sort!(l,by=order,rev=true)
      if isempty(res) t=[l[1]] else t=Ref(l[1]).*elements(Group(res)) end
      p=findfirst(x->order(x)<order(l[1]),t)
      if !isnothing(p)
        t=t[p]
	if order(t)>1 l[1]=t
	else l=deleteat!(l,1)
        end
      else push!(res,l[1])
      end
    end
    res
  end
  # q=Au/k,  ww=words in q images of gens(Au)
  finish=function(q,ww)
    h=Hom(Au,q,map(x->q(x...),ww))
    fusion=map(c->position_class(q,h(c)),classreps(Au))
    ctu=CharTable(Au).irr
    cth=CharTable(q).irr
    ch(c)=map(j->ctu[c,findfirst(isequal(j),fusion)],1:HasType.NrConjugacyClasses(q))
    return Dict(:Au=>q,
      :chars=>map(c->findfirst(i->cth[i,:]==ch(c),axes(cth,1)),chars),
      :gens=>map(x->word(Au,HasType.First(elements(Au),y->h(y)==x)),gens(q)))
  end
  Z=n->ComplexReflectionGroup(n,1,1)
# println("Au=$Au chars=$chars")
  ct=permutedims(CharTable(Au).irr[chars,:])
  cl=filter(i->ct[i,:]==ct[1,:],axes(ct,1))
# println("ct=$ct cl=$cl")
  if length(cl)==1 return Dict(:Au=>Au,:chars=>chars,
                              :gens=>map(x->[x],eachindex(gens(Au)))) end
  ct=permutedims(toM(unique!(sort(toL(ct)))))
# println("ct=$ct")
# k=Subgroup(Au,filter(x->position_class(Au,x) in cl,elements(Au)))
  k=Group(filter(x->position_class(Au,x) in cl,elements(Au)))
  if length(k)==length(Au) return Dict(:Au=>coxgroup(),:chars=>[1],:gens=>[])
  end
# println("Au=$Au k=$k")
  if semisimplerank(Au)==1 return finish(Z(div(length(Au),length(k))),[[1]])
  elseif isabelian(Au/k)
    q=Au/k
    q=Group(AbGens(q))
    h=Hom(Au,q,map(x->Coset(k,x),gens(Au)))
    f=map(x->word(q,h(x)),gens(Au))
#  Print(prod(map(x->Z(order(x)),gens(q)))," ",f,"\n");
    return finish(prod(map(x->Z(order(x)),gens(q))),f)
  else
    p=findfirst(t->all(x->x in reflection_subgroup(Au,t.indices),
                elements(k)),refltype(Au))
    if !isnothing(p)
      p=refltype(Au)[p].indices
      if length(k)==length(reflection_subgroup(Au,p))
       return finish(reflection_subgroup(Au,setdiff(eachindex(gens(Au)),p)),
                     map(i->i in p ? [] : [i],eachindex(gens(Au))))
      elseif length(p)==1
        t=copy(refltype(Au))
        p=findfirst(t->t.indices==p,refltype(Au))
	t[p].p/=length(k)
        return finish(ReflectionGroup(t...),map(x->[x],eachindex(gens(Au))))
      end
    elseif ReflectionName(Au)=="A1xB2" && length(k)==2 && longest(Au) in k
      return finish(coxgroup(:B,2),[[1,2,1,2],[1],[2]])
    end
  end
# Print(" Au=",ReflectionName(Au)," sub=",map(e.Get,gens(k)),"\n");
  error("not implemented ",ReflectionName(Au),chars)
# q:=Au/k; f:=FusionConjugacyClasses(Au,q); Print(" quot=",q," fusion=",f,"\n");
# return rec(Au:=Au,chars:=chars);
end

# When some Springer series have been suppressed/weeded out, we  quotient
# the Au's by the common  kernel of the remaining characters of the Au's.
function AdjustAu!(classes,springerseries)
  for (i, u) in enumerate(classes)
    l=map(s->filter(k->s[:locsys][k][1]==i,eachindex(s[:locsys])),
          springerseries)
#   println(springerseries)
    chars=reduce(vcat,map(j->last.(springerseries[j][:locsys][l[j]]),
                   eachindex(l)))
    f=QuotientAu(u.Au,chars)
#   if Size(Au)<>Size(f.Au) then
#     Print("class ",i,"=",classes[i].name," ",[Au,chars],"=>",f,"\n");
#   fi;
#   print(u.name,":");ds(u.AuAction)
    u.Au=f[:Au]
    if haskey(u,:AuAction)
      R=u.AuAction.group
      if rank(R)==0
        u.AuAction=ExtendedCox(R,[fill(0,0,0) for x in f[:gens]])
      else
       if isempty(f[:gens]) F0s=[reflrep(R,R())]
       else F0s=map(x->prod(u.AuAction.F0s[x]),f[:gens])
       end
       u.AuAction=ExtendedCox(R,F0s)
      end
#     u.AuAction.phis=map(x->prod(u.AuAction.phis[x]),f[:gens])
    end
    k=1
    for j in eachindex(l)
      springerseries[j][:locsys]=copy(springerseries[j][:locsys])
      for s in l[j]
        springerseries[j][:locsys][s][2]=f[:chars][k]
        k+=1
      end
    end
  end
end

"""
`UnipotentClasses(W[,p])`

`W`  should  be  a  `CoxeterGroup`  record  for a Weyl group or `RootDatum`
describing a reductive algebraic group `𝐆`. The function returns a record
containing   information   about   the   unipotent   classes  of  `𝐆`  in
characteristic   `p`  (if   omitted,  `p`   is  assumed   to  be  any  good
characteristic for `𝐆`). This contains the following fields:

`group`: a pointer to `W`

`p`: the characteristic of the field for which the unipotent classes were
computed. It is `0` for any good characteristic.

`orderClasses`:  a list describing the Hasse diagram of the partial order
induced   on   unipotent   classes   by   the  closure  relation.  That  is
`.orderclasses[i]`  is the list of `j` such that `C̄ⱼ⊋ C̄ᵢ`  and  there  is
no  class  `Cₖ`  such  that `C̄ⱼ⊋ C̄ₖ⊋ C̄ᵢ`.

`classes`:  a  list  of  records  holding information for each unipotent
class (see below).

`springerSeries`:  a list of records, each  of which describes a Springer
series  of `𝐆`.

The  records  describing  individual  unipotent  classes have the following
fields:

`name`: the name of the unipotent class.

`parameter`:  a parameter  describing the  class (for  example, a partition
describing the Jordan form, for classical groups).

`Au`: the group `A(u)`.

`dynkin`:  present in good characteristic; contains the Dynkin-Richardson
diagram,  given  as  a  list  of  0,1,2  describing  the coefficient on the
corresponding simple root.

`red`:  the reductive part of ``C_𝐆(u)``.

`dimBu`:  the dimension of the variety `𝓑ᵤ`.

The  records for classes contain additional  fields for certain groups: for
instance,  the names given to classes by Mizuno in `E₆, E₇, E₈` or by Shoji
in `F₄`.

The  records  describing  individual  Springer  series  have  the following
fields:

`levi`:the  indices of the  reflections corresponding to  the Levi subgroup
`𝐋`  where  lives  the  cuspidal  local  system `ι` from which the Springer
series is induced.

`relgroup`: The relative Weyl group ``N_𝐆(𝐋,ι)/𝐋``. The first series is the
principal series for which `.levi=[]` and `.relgroup=W`.

`locsys`:  a  list  of  length  `NrConjugacyClasses(.relgroup)`, holding in
`i`-th  position a  pair describing  which local  system corresponds to the
`i`-th  character of  ``N_𝐆(𝐋,ι)``. The  first element  of the  pair is the
index  of the concerned unipotent class `u`, and the second is the index of
the corresponding character of `A(u)`.

`Z`:  the central character associated  to the Springer series, specified
by its value on the generators of the centre.

```julia-repl
julia> W=rootdatum(:sl,4)
A₃

julia> uc=UnipotentClasses(W);

julia> uc.classes
5-element Array{Gapjm.Ucl.UnipotentClass,1}:
 UnipotentClass(1111)
 UnipotentClass(211)
 UnipotentClass(22)
 UnipotentClass(31)
 UnipotentClass(4)
```

The  `show`  function  for  unipotent  classes  accepts  all the options of
`formatTable`  and  of  `charnames`.  Giving  the  option  `mizuno`  (resp.
`shoji`)  uses  the  names  given  by  Mizuno  (resp.  Shoji) for unipotent
classes.  Moreover,  there  is  also  an  option  `fourier` which gives the
correspondence  tensored  with  the  sign  character  of each relative Weyl
group, which is the correspondence obtained via a Fourier-Deligne transform
(here  we assume that  `p` is very  good, so that  there is a nondegenerate
invariant  bilinear  form  on  the  Lie  algebra, and also one can identify
nilpotent orbits with unipotent classes).

Here is how to display the non-cuspidal part of the Springer correspondence
of  the unipotent  classes of  `E₆` using  the notations  of Mizuno for the
classes  and those  of Frame  for the  characters of  the Weyl group and of
Spaltenstein  for the characters  of `G₂` (this  is convenient for checking
our data with the original paper of Spaltenstein):

```julia-rep1
julia> uc=UnipotentClasses(rootdatum(:E6sc));

julia> xprint(uc;cols=[5,6,7],spaltenstein=true,frame=true,mizuno=true,
      order=false)
UnipotentClasses(E₆)
     u│            E₆(E₆₍₎) G₂(E₆₍₁₃₅₆₎=A₂×A₂)/ζ₃ G₂(E₆₍₁₃₅₆₎=A₂×A₂)/ζ₃²
──────┼──────────────────────────────────────────────────────────────────
E₆    │                1:1ₚ                  ζ₃:1                  ζ₃²:1
E₆(a₁)│                1:6ₚ                ζ₃:ε_c                ζ₃²:ε_c
D₅    │              Id:20ₚ
A₅+A₁ │        -1:15ₚ 1:30ₚ                 ζ₃:θ′                 ζ₃²:θ′
A₅    │              1:15_q                 ζ₃:θ″                 ζ₃²:θ″
D₅(a₁)│              Id:64ₚ
A₄+A₁ │              Id:60ₚ
D₄    │              Id:24ₚ
A₄    │              Id:81ₚ
D₄(a₁)│111:20ₛ 3:80ₛ 21:90ₛ
A₃+A₁ │              Id:60ₛ
2A₂+A₁│               1:10ₛ                 ζ₃:εₗ                 ζ₃²:εₗ
A₃    │             Id:81ₚ′
A₂+2A₁│             Id:60ₚ′
2A₂   │              1:24ₚ′                  ζ₃:ε                  ζ₃²:ε
A₂+A₁ │             Id:64ₚ′
A₂    │      11:15ₚ′ 2:30ₚ′
3A₁   │            Id:15_q′
2A₁   │             Id:20ₚ′
A₁    │              Id:6ₚ′
1     │              Id:1ₚ′
```
"""
function UnipotentClasses(t::TypeIrred,p=0)
  uc=getchev(t,:UnipotentClasses,p)
  if uc===nothing error("no UnipotentClasses for ",t) end
  rank=PermRoot.rank(t)
  classes=UnipotentClass[]
  for u in uc[:classes] # fill omitted fields
    name=u[:name]
    parameter= haskey(u,:parameter) ? u[:parameter] : u[:name]
    dimBu= haskey(u,:dimBu)  ? u[:dimBu] : -1
    if haskey(u,:dynkin)
      c=haskey(t,:orbit) ? cartan(t.orbit[1]) : cartan(t)
      weights=toM(roots(c))*u[:dynkin]
      n0=count(iszero,weights)
      if dimBu==-1 dimBu=n0+div(count(isone,weights),2)
      elseif dimBu!=n0+div(count(isone,weights),2) error("theory")
      end
      n0=2*n0-count(isequal(2),weights)
      u[:dimunip]=2*dimBu-n0
      u[:dimred]=n0+rank
    elseif haskey(u,:red)
      u[:dimred]=dimension(u[:red])
      u[:dimunip]=2*dimBu+rank-u[:dimred]
    elseif haskey(u,:dimred)
      u[:dimunip]=2*dimBu+rank-u[:dimred]
    end
    delete!.(Ref(u),[:name,:parameter,:dimBu])
    push!(classes,UnipotentClass(name,parameter,dimBu,u))
  end
  springerseries=uc[:springerSeries]
  for s in springerseries
    if isempty(s[:levi]) s[:levi]=Int[] end
#   s[:levi]=PermRoot.indices(t)[s[:levi]]
    s[:locsys]=Vector{Int}.(s[:locsys])
  end
  orderclasses=map(x->isempty(x) ? Int[] : x,uc[:orderClasses])
  delete!.(Ref(uc),[:classes,:orderClasses,:springerSeries])
# uc[:spets]=t
  UnipotentClasses(classes,p,Poset(orderclasses),springerseries,uc)
end

Base.length(uc::UnipotentClasses)=length(uc.classes)

function UnipotentClasses(W,p=0)
# println("UnipotentClasses(",W,")")
  spetscase=W isa Spets
  if spetscase
    WF=W
    W=Group(WF)
    t=refltype(W)
    l=map(x->findfirst(y->any(z->sort(x.indices)==sort(z.indices),
                                       y.orbit),refltype(WF)),t)
    uc=UnipotentClasses.(refltype(WF)[l],p)
   else WF=spets(W)
    t=refltype(W)
    uc=UnipotentClasses.(t,p)
  end
  if isempty(t)
    classes=[UnipotentClass("",[],0,
        Dict(:Au=>coxgroup(),:dynkin=>[],:balacarter=>[],
             :dimunip=>0,:red=>torus(rank(W))))]
    uc=[UnipotentClasses(classes,p,Poset([Int[]]),
      [Dict(:Z=>Int[],:levi=>Int[],:locsys=>[[1,1]],:relgroup=>coxgroup())],
      Dict{Symbol,Any}(:spets=>W))]
    l=Vector{Int}[]
  else
    classes=map(cartesian(map(x->x.classes,uc)...)) do v
      l=PermRoot.indices.(t)
      if length(v)==1 && issorted(l[1]) u=deepcopy(v[1])
      else
        u=UnipotentClass(join(map(x->x.name,v),","),map(x->x.parameter,v),
                         sum(map(x->x.dimBu,v)),Dict{Symbol,Any}())
        u.Au=prod(x->x.Au,v)
        if all(x->haskey(x,:dimred),v)
          u.dimred=sum(x->x.dimred,v)+rank(W)-semisimplerank(W) 
        end
        if all(x->haskey(x,:dimunip),v)
          u.dimunip=sum(x->x.dimunip,v) end
        if all(x->haskey(x,:red),v)
#         println(map(x->x.red,v))
          u.red=reduce(Cosets.extprod,map(x->x.red,v)) end
        if all(x->haskey(x,:AuAction),v)
          u.AuAction=prod(x->x.AuAction,v) end
        if all(x->haskey(x,:dynkin),v)
          u.dynkin=zeros(Int,sum(length,l))
          for i in 1:length(l) u.dynkin[l[i]]=v[i].dynkin end
        end
      end
      if all(x->haskey(x,:balacarter),v)
        u.balacarter=reduce(vcat,[map(j->j>0 ? x[j] : -x[-j],
                    v[i].balacarter) for (i,x) in enumerate(l)])
      end
      if haskey(u, :red)
        if !(u.red isa Spets) u.red=spets(u.red) end
        if rank(W)>semisimplerank(W) && haskey(u, :red)
          u.red=Cosets.extprod(u.red,radical(WF))
          if haskey(u,:AuAction)
            T=radical(W)
            u.AuAction=ExtendedCox(u.AuAction.group*T,
               map(x->cat(x,reflrep(T,T()),dims=(1,2)),u.AuAction.F0s))
          end
        end
      end
      u
    end
  end
  if iszero(p) && !haskey(classes[1],:balacarter)
    bc=BalaCarterLabels(W)
#   println("W=$W bc=$bc")
#   println(map(u->u.dynkin,classes))
    for u in classes
      pp=findfirst(p->p[1]==u.dynkin,bc)
      if isnothing(pp) error("not found:",u.dynkin) end
      u.balacarter=bc[findfirst(p->p[1]==u.dynkin,bc)][2]
    end
  end
  ll=length.(uc)
  orderclasses=map(cartesian(map(x->1:x,ll)...)) do v
    o=cartesian(map(j->vcat(hasse(uc[j].orderclasses)[v[j]],[v[j]]),1:length(v))...)
    o=map(x->HasType.PositionCartesian(ll,x),o)
    setdiff(o,[HasType.PositionCartesian(ll,v)])
  end
  springerseries=map(cartesian(map(x->x.springerseries,uc)...)) do v
#   if isempty(v) return Dict(:Z=>[],:levi=>[],:locsys=>[[1,1]])
    if length(v)==1 
      if !isempty(l)
        v[1]=deepcopy(v[1]); v[1][:levi]=l[1][v[1][:levi]]; 
      end
      return v[1]
    end
    s=Dict{Symbol,Any}(:levi=>reduce(vcat,map(i->l[i][v[i][:levi]],eachindex(v))))
    s[:Z]=reduce(vcat,getindex.(v,:Z))
    s[:locsys]=map(cartesian(getindex.(v,:locsys)...)) do v
      v=collect.(zip(v...))
      u=map(i->HasType.NrConjugacyClasses(uc[i].classes[v[1][i]].Au),
              eachindex(v[1]))
      [HasType.PositionCartesian(ll,v[1]),HasType.PositionCartesian(u,v[2])]
    end
    if all(haskey.(v,:parameter)) s[:parameter]=getindex.(v,:parameter) end
    s[:relgroup]=prod(getindex.(v,:relgroup))
    if length(v)==1
      for k in setdiff(keys(v[1]),[:levi,:Z,:locsys,:parameter])
        s[k]=v[1][k]
      end
    end
    s
  end
  if length(uc)==1 prop=uc[1].prop else prop=Dict{Symbol,Any}() end
  prop[:spets]=spetscase ? WF : W
  if spetscase
    springerseries=filter(x->sort(inclusion(W,x[:levi]).^WF.phi)==
                             sort(inclusion(W,x[:levi])),springerseries)
  end
# To deal with a general group intermediate between Gad and Gsc, we discard
# the  Springer series  corresponding to  a central  character which is not
# trivial on the fundamental group (seen as a subgroup of ZGsc)
# algebraic_centre(W).descAZ returns the generators of the fundamental group
# of  the  algebraic  group  W  as  words  in  generators  of  the absolute
# fundamental group.
  if !all(x->Set(x[:Z])==Set([1]),springerseries)
    springerseries=filter(s->all(y->prod(s[:Z][y])==1,
             algebraic_centre(W)[:descAZ]),springerseries)
    AdjustAu!(classes,springerseries)
  end
# println(springerseries[1])
  s=springerseries[1]
  if spetscase
    s[:relgroup]=relative_coset(WF,s[:levi])
    s[:locsys]=s[:locsys][charinfo(s[:relgroup])[:charRestrictions]]
  end
  l=filter(i->any(y->i==y[1],s[:locsys]),1:length(classes))
  s[:locsys]=map(((c,s),)->[findfirst(==(c),l),s],s[:locsys])
  # for now only springerseries[1] properly twisted
  for s in springerseries[2:end]
    if spetscase
      s[:relgroup]=relative_coset(WF,s[:levi])
      s[:locsys]=s[:locsys][charinfo(s[:relgroup])[:charRestrictions]]
    end
    s[:locsys]=map(((c,s),)->[findfirst(==(c),l),s],s[:locsys])
  end
  classes=classes[l]
  AdjustAu!(classes,springerseries)
  orderclasses=Poset(hasse(restricted(Poset(orderclasses),l)))
  orderclasses.prop[:label]=(io,n)->name(io,classes[n])
  ucl=UnipotentClasses(classes,p,orderclasses,springerseries,prop)
  ucl
end

function showcentralizer(io::IO,u)
  c=""
  function AuName(u)
    if length(u.Au)==1 return "" end
    res=haskey(u,:AuAction) ||
        (haskey(u,:dimred) && iszero(u.dimred)) ? "." : "?"
        res*=reflection_name(IOContext(io,:Au=>true),u.Au)
  end
  if haskey(u,:dimunip)
    if u.dimunip>0 c*=HasType.Format(Mvp(:q)^u.dimunip,io.dict) end
  else c*="q^?" end
  if haskey(u,:AuAction)
    if rank(u.red)>0
      c*="."
      if all(isone,u.AuAction.F0s)
        c*=sprint(show,u.red;context=io)*AuName(u)
      elseif length(u.Au)==1 ||
         length(u.Au)==length(Group(u.AuAction.phis...))
        if length(u.Au)==1 || isone(u.red.phi)
          c*=sprint(show,u.AuAction;context=io)
        else
          c*="["*sprint(show,u.red;context=io)*"]"*
                 sprint(show,u.AuAction;context=io)
        end
      else
        c*=reflection_name(io,u.AuAction)*AuName(u)
      end
    else
      c*=AuName(u)
    end
  elseif haskey(u,:red)
    n=sprint(show,u.red;context=io)
    if n!="." c*="."*n end
    c*=AuName(u)
  else
    if haskey(u,:dimred) if u.dimred>0 c*="[red dim ",u.dimred,"]" end
    else c*="[red??]"
    end
    c*=AuName(u)
  end
  c=replace(c,r"^\."=>"")
  c=replace(c,r"\.\."=>".")
  replace(c,"()"=>"")
  c
end

function Base.show(io::IO, ::MIME"text/html", uc::UnipotentClasses)
  show(IOContext(io,:TeX=>true),uc)
end

function Base.show(io::IO,uc::UnipotentClasses)
  TeX=get(io,:TeX,false)
  repl=get(io,:limit,false)
  deep=get(io,:typeinfo,false)!=false
  printTeX(io,"UnipotentClasses(\$",uc.prop[:spets],"\$)")
  if !(repl||TeX) || deep return end
  print(io,"\n")
  if get(io,:order,true) println(io,uc.orderclasses) end
  sp = map(copy, uc.springerseries)
  if get(io,:fourier,false)
    for p in sp p[:locsys] = p[:locsys][DetPerm(p[:relgroup])] end
  end
  WF=uc.prop[:spets]
  if WF isa Spets W=Group(WF)
  else W=WF;WF=spets(W) end
  if uc.p!=0 || !any(x->haskey(x,:balacarter),uc.classes)
    io=IOContext(io,:balacarter=>false)
  end
  tbl = map(uc.classes)do u
    res= iszero(uc.p) ? [joindigits(u.dynkin)] : String[]
    push!(res, string(u.dimBu))
    if get(io,:balaCarter,true)
      if haskey(u, :balacarter)
        b=fill('.',coxrank(W))
        for i in u.balacarter if i>0 b[i]='2' else b[-i]='0' end end
      else
        b=fill('?',coxrank(W))
      end
      push!(res, String(b))
    end
    if get(io,:centralizer,true)
      push!(res,showcentralizer(io,u))
    end
    if get(io,:springer,true)
      i=findfirst(isequal(u),uc.classes)
      cc(ss)=map(function (i)
                c1 = charnames(io,u.Au)[ss[:locsys][i][2]]
                c2 = charnames(io,ss[:relgroup])[i]
                (c1=="") ? c2 : c1*":"*c2
           end, findall(y->y[1]==i,ss[:locsys]))
      append!(res, map(ss->join(cc(ss), TeX ? "\\kern 0.8em " : " "),sp))
    end
    res
  end
  col_labels= String[]
  if iszero(uc.p)
    push!(col_labels, TeX ? "\\mbox{Dynkin-Richardson}" : "D-R")
  end
    push!(col_labels, TeX ? "\\dim{\\cal B}_u" : "dBu")
  if get(io,:balaCarter,true)
     push!(col_labels, TeX ? "\\mbox{Bala-Carter}" : "B-C")
  end
  if get(io,:centralizer,true)
     push!(col_labels, TeX ? "C_{\\bf G}(u)" : "C(u)")
  end
  if get(io,:springer,true)
   append!(col_labels,
      map(function (ss,)
        res=string(sprint(show,ss[:relgroup];context=io),"(",
          sprint(show,subspets(WF,ss[:levi]);context=io),")")
        if !all(isone,ss[:Z])
          res*=string("/", join(map(q->sprint(show,q;context=io),ss[:Z]),","))
        end
        return res
    end, sp))
  end
  row_labels=name.(Ref(io),uc.classes)
  if get(io,:rows,false)==false
    io=IOContext(io,:rows=>sortperm(map(x->x.dimBu, uc.classes)))
  end
  format(io,toM(tbl);rows_label="u",col_labels=col_labels,row_labels=row_labels)
end

# decompose tensor product of characters (given as their indices in CharTable)
function DecomposeTensor(W,c::Int...)
  ct=CharTable(W)
# println("eltype=",eltype(irr))
  decompose(ct,vec(prod(view(ct.irr,collect(c),:),dims=1)))
end

struct ICCTable
  prop::Dict{Symbol,Any}
end

Base.getproperty(t::ICCTable,k::Symbol)=getfield(t,:prop)[k]
Base.setproperty!(t::ICCTable,k::Symbol,x)=getfield(t,:prop)[k]=x


"""
`ICCTable(uc,seriesNo=1;q=Pol())`

This function gives the table of decompositions of the functions ``X_ι`` in
terms  of the functions ``Y_ι``. Here `ι` is a `𝐆`-equivariant local system
on  the  class  `C`  of  a  unipotent  element  `u`. Such a local system is
parametrized  by the pair  `(u,ϕ)` of `u`  and a character  of the group of
components   `A(u)`   of   ``C_𝐆   (u)``.   The  function  ``Y_ι``  is  the
characteristic   function  of  this   local  system  and   ``X_ι``  is  the
characteristic   function  of  the  corresponding  intersection  cohomology
complex  on `C̄`. The  Springer correspondence says  that the local systems
can  also be  indexed by  characters of  a relative  Weyl group.  Since the
coefficient of `Xᵪ` on `Yᵩ` is `0` if `χ` and `φ` are not characters of the
same  relative Weyl group (are not in  the same Springer series), the table
is  for one  Springer series,  specified by  the argument  'seriesNo' (this
defaults  to 'seriesNo=1' which is the principal series). The decomposition
multiplicities  are graded,  and are  given as  polynomials in one variable
(specified by the argument `q`; if not given `Pol()` is assumed).

```julia-repl
julia> t=ICCTable(uc)
Coefficients of Xᵪ on Yᵩ for A₃
     │4 31 22 211 1111
─────┼─────────────────
X4   │1  1  1   1    1
X31  │0  1  1  Φ₂   Φ₃
X22  │0  0  1   1   Φ₄
X211 │0  0  0   1   Φ₃
X1111│0  0  0   0    1
```
In  the  above  the  multiplicities  are  given  as  products of cyclotomic
polynomials  to display them  more compactly. However  the format of such a
table can be controlled more precisely.

For  instance,  one  can  ask  to  not  display  the entries as products of
cyclotomic polynomials:

```julia-rep1
julia> xprint(t;cycpol=false)
Coefficients of Xᵪ on Yᵩ for A3
     │4 31 22 211   1111
─────┼───────────────────
X4   │1  1  1   1      1
X31  │0  1  1 q+1 q²+q+1
X22  │0  0  1   1   q²+1
X211 │0  0  0   1 q²+q+1
X1111│0  0  0   0      1
```

Since `show` uses the function `format` for tables, all the options of this
function  are  also  available.  We  can  use  this to restrict the entries
displayed  to a  given sublist  of the  rows and  columns (here the indices
correspond  to the number  in Chevie of  the corresponding character of the
relative Weyl group of the given Springer series):

```julia-rep1
julia> uc=UnipotentClasses(coxgroup(:F,4));
julia> t=ICCTable(uc);
julia> sh=[13,24,22,18,14,9,11,19];
julia> show(IOContext(stdout,:rows=>sh,:cols=>sh,:limit=>true),t);
Coefficients of Xᵪ on Yᵩ for F₄
      │A₁+Ã₁ A₂ Ã₂ A₂+Ã₁ Ã₂+A₁ B₂⁽¹¹⁾ B₂ C₃(a₁)⁽¹¹⁾
──────┼─────────────────────────────────────────────
Xφ₉‚₁₀│    1  0  0     0     0      0  0          0
Xφ″₈‚₉│    1  1  0     0     0      0  0          0
Xφ′₈‚₉│    1  0  1     0     0      0  0          0
Xφ″₄‚₇│    1  1  0     1     0      0  0          0
Xφ′₆‚₆│   Φ₄  1  1     1     1      0  0          0
Xφ₄‚₈ │   q²  0  0     0     0      1  0          0
Xφ″₉‚₆│   Φ₄ Φ₄  0     1     0      0  1          0
Xφ′₄‚₇│   q²  0 Φ₄     0     1      0  0          1
```

The   function  'ICCTable'  returns  an   object  with  various  pieces  of
information which can help further computations.

`.scalar`:  this contains the table of  multiplicities `Pᵪᵩ` of the `Xᵪ` on
the  `Yᵩ`.  One  should  pay  attention  that  by default, the table is not
displayed  in the same order as the  stored |.scalar|, which is in order in
Chevie  of  the  characters  in  the  relative  Weyl  group;  the  table is
transposed,  then lines  and rows  are sorted  by `dimBu,class  no,index of
character in A(u)` while displayed.

`.group`: The group `W`.

`.relgroup`: The relative Weyl group for the Springer series.

`.series`: The index of the Springer series given for `W`.

`.dimBu`: The list of `dim𝓑ᵤ` for each local system `(u,φ)` in the series.

`:L`:  The matrix of  (unnormalized) scalar products  of the functions `Yᵩ`
with themselves, that is the `(φ,χ)` entry is ``∑_{g∈𝐆(𝔽_q)} Yᵩ(g) Yᵪ(g)``.
This  is thus a symmetric, block-diagonal  matrix where the diagonal blocks
correspond  to  geometric  unipotent  conjugacy  classes.  This  matrix  is
obtained as a by-product of Lusztig's algorithm to compute `Pᵩᵪ`.
"""
function ICCTable(uc::UnipotentClasses,i=1;q=Pol())
  W=uc.prop[:spets] # W=Group(uc.spets)
  if W isa Spets W=W.W end
  ss=uc.springerseries[i]
  res=ICCTable(Dict(:spets=>uc.prop[:spets],:relgroup=>ss[:relgroup],
                    :series=>i,:q=>q,:p=>uc.p))
  if haskey(ss,:warning) println("# ",ss[:warning])
    res.warning=ss[:warning]
  end
# We are going to solve the equation in "unipotent support", page 151
# ᵗPΛP=Ω  where $Λ_{i,j}$ is  $∑_{g∈ G^F} Yᵢ(g)Ȳⱼ(g)$ and $Ω_{i,j}$ is equal
# to # $|Z^0(G^F)|q^{-semisimple rank L}|G^F|/P(W_G(L))
#  q^{-bᵢ-bⱼ}FakeDegree(χᵢ⊗χⱼ⊗sgn)$
# where $P(W_G(L))$ is the Poincare polynomial $∏ᵢ(q^{dᵢ}-1)$
# where $dᵢ$ are the reflection degrees of $W_G(L)$
# res[:scalar] is the matrix $P$
  R=ss[:relgroup]
  ct=CharTable(R)
  var=q
  q=Pol()
  f=fakedegrees(R,q)
  k=charinfo(R)[:positionDet]
  n=length(f)
# Partition on characters of ss.relgroup induced by poset of unipotent classes
  res.dimBu=map(x->uc.classes[x[1]].dimBu,ss[:locsys])
  res.blocks=HasType.CollectBy(eachindex(ss[:locsys]),-res.dimBu)
  # matrix of q^{-bᵢ-bⱼ}*fakedegree(χᵢ ⊗ χⱼ ⊗ sgn)
  tbl=bigcell_decomposition([q^(-res.dimBu[i]-res.dimBu[j])*
                            sum(map(*,f,DecomposeTensor(R,i,j,k)))
     for i in 1:n,j in 1:n], res.blocks) # //1 needed in D7
  res.scalar=tbl[1]
  res.locsys=ss[:locsys]
# res[:L]=tbl[2]*GenericOrder(W,q)/prod(ReflectionDegrees(R),d->q^d-1)/
#   q^(W.semisimplerank-R.semisimplerank);
  res.L=tbl[2]*q^(W.N+semisimplerank(R)-semisimplerank(W))
  res.uc=uc
  if haskey(ss,:parameter) res.parameter=ss[:parameter]
  else res.parameter=1:length(ss[:locsys])
  end
  if !(var isa Pol)
    res.scalar=improve_type(map(x->x(var),res.scalar))
    res.L=improve_type(map(x->x(var),res.L))
  end
  res
end

function Base.show(io::IO,x::ICCTable)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if !(repl || TeX) print(io,"ICCTable(",x.uc,",",x.series,")"); return end
  text="Coefficients of \$X_\\chi\$ on \$Y_\\phi\$ for \$"*
        sprint(show,x.relgroup;context=io)*"\$\n"
  if TeX text*="\\medskip\n\n" else text=fromTeX(io,text) end
  print(io,text)
  if get(io,:cols,false)==false && get(io,:rows,false)==false
    rows=collect(eachindex(x.dimBu))
    sort!(rows,by=i->[x.dimBu[i],x.locsys[i]])
    io=IOContext(io,:rows=>rows,:cols=>rows)
  end
  tbl=get(io,:cycpol,true) ? map(CycPol,x.scalar) : x.scalar
  col_labels=map(((c,s),)->name(IOContext(io,:locsys=>s),x.uc.classes[c]),
                  x.locsys)
  rowLabels=map(x->TeX ? "X_{$x}" : "X$x",charnames(io,x.relgroup))
  format(io,permutedims(tbl),row_labels=rowLabels,col_labels=col_labels)
end

struct XTable
  prop::Dict{Symbol,Any}
end

Base.getproperty(t::XTable,k::Symbol)=getfield(t,:prop)[k]
Base.setproperty!(t::XTable,k::Symbol,x)=getfield(t,:prop)[k]=x

function repair(u::Matrix) # can disappear with julia1.6
  for i in 1:length(u)
    if !isassigned(u,i) u[i]=0 end
  end
  u
end
  
# XTable(uc[,opt]) values of X̃ᵪ on unipotent classes or local systems
# Note that c_ι=βᵤ+(rkss L_\CI)/2
#
# Formatting: options of FormatTable + [.classes, .CycPol]
function XTable(uc::UnipotentClasses;q=Pol(),classes=false)
# println("here uc=",uc)
  pieces=map(i->ICCTable(uc,i;q=q),eachindex(uc.springerseries))
  greenpieces=map(x->x.scalar*toM(HasType.DiagonalMat(q.^x.dimBu...)),pieces)
  l=vcat(getproperty.(pieces,:locsys)...)
  p=inv(sortPerm(l))
  res=XTable(Dict(
    :scalar=>permutedims(repair(cat(greenpieces...,dims=(1,2)))^p),
    :uc=>uc,
    :Y=>^(repair(cat(getproperty.(pieces,:L)...,dims=(1,2))),p,dims=(1,2)),
    :parameter=>vcat(getproperty.(pieces,:parameter)...),
    :relgroups=>getindex.(uc.springerseries,:relgroup)))
  if classes
    res.scalar*=E(1)
    res.cardClass=zeros(eltype(res.scalar),length(l))//1
    res.classes=l^p
    for i in eachindex(uc.classes)
      Au=uc.classes[i].Au
      b=filter(j->res.classes[j][1]==i,eachindex(res.classes))
 #    println("i=",i," b=",b," Au=",Au)
      res.scalar[:,b]*=CharTable(Au).irr
      res.cardClass[b]=res.Y[[b[charinfo(Au)[:positionId]]],b]*CharTable(Au).irr
      res.cardClass[b]=map((x,y)->x*y//length(Au),
                             res.cardClass[b],classinfo(Au)[:classes])
    end
    res.scalar=improve_type(res.scalar)
  else
    res.locsys=l^p
  end
  res
end

function Base.show(io::IO,x::XTable)
  if !get(io,:limit,false) && !get(io,:TeX,false)
    print(io,"XTable(",x.uc.prop[:spets],",variable=",x.q,")")
    return
  end
  class=haskey(getfield(x,:prop),:classes)
  print(io,"Values of character sheaves X̃ᵪ on")
  rowLabels=vcat(map(g->map(n->fromTeX(io,"X^{"*sprint(show,g;context=io)
             *"}_{"*n*"}"),charnames(io,g)),x.relgroups)...)
  rowsLabel=fromTeX(io,"\\tilde X_\\chi\\"*(class ? "class" : "locsys"))
  if class
    print(io," unipotent classes\n")
    columnLabels=map(p->name(IOContext(io,:class=>p[2]),x.uc.classes[p[1]]),
                     x.classes)
  else print(io," local systems\n")
    columnLabels=map(p->name(IOContext(io,:locsys=>p[2]),x.uc.classes[p[1]]),x.locsys)
  end
  tbl=x.scalar
  if get(io,:cycpol,false) tbl=CycPol.(tbl) end
  format(io,tbl,row_labels=rowLabels,col_labels=columnLabels,rows_label=rowsLabel)
end

struct GreenTable
  prop::Dict{Symbol,Any}
end

Base.getproperty(t::GreenTable,k::Symbol)=getfield(t,:prop)[k]
Base.setproperty!(t::GreenTable,k::Symbol,x)=getfield(t,:prop)[k]=x

# GreenTable(uc;q=Pol())
# values of Green functions Q^\CI_{wF} on unipotent classes
# method: use formula DLM3 (3.1)
# Lines indexed by (\CI,wF). Columns by unip. classes or local systems
function GreenTable(uc::UnipotentClasses;q=Pol())
  t=GreenTable(getfield(XTable(uc;classes=true,q=q),:prop))
  m=cat(map(g->permutedims(CharTable(g).irr),t.relgroups)...;dims=(1,2))
  t.scalar=m*t.scalar
  t.indices=Vector{Int}[]
  i=0
  for g in t.relgroups
    push!(t.indices,(1:nconjugacy_classes(g))+i)
    i+=nconjugacy_classes(g)
  end
  t
end

function Base.show(io::IO,x::GreenTable)
  if !get(io,:limit,false) && !get(io,:TeX,false)
    print(io,"GreenTable(",x.uc.prop[:spets],",variable=",x[:q],")")
    return
  end
  print(io,"Values of Green functions Q_{wF} on unipotent classes\n")
  rowLabels=vcat(map(x.relgroups) do g
    classnames=map(x->fromTeX(io,x),classinfo(g)[:classnames])
    map(n->string("Q^{",sprint(show,g;context=io),"}_{",n,"}"),classnames)
    end...)
  rowsLabel="Q^I_{wF}\\class"
  columnLabels=map(p->name(IOContext(io,:class=>p[2]),x.uc.classes[p[1]]),
                     x.classes)
  tbl=x.scalar
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  format(io,tbl,row_labels=rowLabels,col_labels=columnLabels,rows_label=rowsLabel)
end

struct TwoVarGreenTable
  prop::Dict{Symbol,Any}
end

Base.getproperty(t::TwoVarGreenTable,k::Symbol)=getfield(t,:prop)[k]
Base.setproperty!(t::TwoVarGreenTable,k::Symbol,x)=getfield(t,:prop)[k]=x

# two-variable green functions table.
# for now only implemented when W is split.
function TwoVarGreen(W,L)
  if !(W isa Spets) W=spets(W) end
  if !(L isa Spets) L=spets(L) end
  uG=UnipotentClasses(W)
  uL=UnipotentClasses(L)
  tG=GreenTable(uG)
  tL=GreenTable(uL)
  q=Pol()
  mm=map(eachindex(uL.springerseries))do i
    s=uL.springerseries[i]
    p=findfirst(S->S[:levi]==inclusion(L,s[:levi]) &&
                S[:Z][1]==s[:Z][1],uG.springerseries)
    if isnothing(p) error("not found ",s[:levi]) end
    RG=relative_coset(W,inclusion(L,s[:levi]))
    RLF=relative_coset(L,s[:levi])
    RL=Group(RLF)
    l=map(x->findfirst(==(x),Group(RG).prop[:relativeIndices]),
      inclusion(L,RL.prop[:relativeIndices]))
    if nothing in l error("not implemented") end
    RLF=subspets(RG,convert(Vector{Int},l),
                Group(RG).prop[:MappingFromNormalizer](L.phi))
    RL=Group(RLF)
    f=fusion_conjugacy_classes(RLF,RG)
    cl=classreps(RLF)
    d=map(cl)do w
      if isempty(s[:levi]) pw=w
      elseif isone(w) pw=L.phi
      else pw=word(RLF,w)
        if isempty(pw) pw=L.phi
        else pw=prod(Group(RG).prop[:parentMap][pw])*L.phi
        end
      end
      Lo=subspets(L,Int.(s[:levi]),pw/L.phi)
      r=map(last,filter(x->isone(first(x)),degrees(Lo)))
      prod(x->q-x,r)/length(centralizer(RL,w))
    end
    q^(length(s[:levi]))*permutedims(tL.scalar[tL.indices[i],:])*
    toM(HasType.DiagonalMat(d...))*conj(tG.scalar[tG.indices[p][f],:])
  end
  oL=generic_order(L,q)
  mm=toM(map((x,y)->x*y/oL,eachrow(sum(mm)),tL.cardClass))
  res=TwoVarGreenTable(Dict(:W=>W,:L=>L,:scalar=>mm,:uL=>uL,:uG=>uG))
  res.classL=tL.classes
  res.classG=tG.classes
  res
end

function Base.show(io::IO, ::MIME"text/html", x::TwoVarGreenTable)
  show(IOContext(io,:TeX=>true),x)
end

function Base.show(io::IO,x::TwoVarGreenTable)
  if !get(io,:limit,false) && !get(io,:TeX,false)
    print(io,"TwoVarGreenTable(",x.W,",",x.L,")")
    return
  end
  print(io,"Values of two-variable Green functions for ",x.W," and ",x.L,"\n")
  rowLabels=map(x.classL)do p
    name(IOContext(io,:class=>p[2]),x.uL.classes[p[1]])
  end
  columnLabels=map(x.classG)do p
    name(IOContext(io,:class=>p[2]),x.uG.classes[p[1]])
  end
  tbl=improve_type(x.scalar)
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  format(io,tbl,row_labels=rowLabels,col_labels=columnLabels,rows_label="v/u")
end

"""
'special_pieces(<uc>)'

The  special  pieces  forme  a  partition  of  the  unipotent  variety of a
reductive  group `𝐆` which was defined  the first time in [Spaltenstein1982
chap.  III](biblio.htm#spalt82)  as  the  fibers  of  `d^2`, where `d` is a
"duality  map". Another definition is as the  set of classes in the Zariski
closure  of a special class  and not in the  Zariski closure of any smaller
special  class, where  a special  class in  the support  of the  image of a
special character by the Springer correspondence.

Each  piece is a union of unipotent  conjugacy classes so is represented in
Chevie  as a  list of  class numbers.  Thus the  list of  special pieces is
returned  as  a  list  of  lists  of  class  numbers. The list is sorted by
increasing  piece dimension, while each piece is sorted by decreasing class
dimension, so the special class is listed first.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> special_pieces(UnipotentClasses(W))
3-element Array{Array{Int64,1},1}:
 [1]
 [4, 3, 2]
 [5]

julia> special_pieces(UnipotentClasses(W,3))
3-element Array{Array{Int64,1},1}:
 [1]
 [4, 3, 2, 6]
 [5]
```

The   example  above  shows  that  the  special  pieces  are  different  in
characteristic 3.
"""
function special_pieces(uc)
  W=uc.prop[:spets]
  ch=charinfo(W)
  specialch=findall(iszero,ch[:a]-ch[:b]) # special characters of W
  specialc=first.(uc.springerseries[1][:locsys][specialch])
  sort!(specialc,by=c->-uc.classes[c].dimBu)
  m=permutedims(incidence(uc.orderclasses))
  map(eachindex(specialc))do i
    p=m[specialc[i],:]
    for j in 1:i-1 p.&=.!m[specialc[j],:] end
    p=eachindex(p)[p]
    sort!(p,by=c->uc.classes[c].dimBu)
    p
  end
end
"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%'UnipotentValues(<W>,<w>)'
"""

end
