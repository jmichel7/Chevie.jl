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

The  Dynkin-Richarson diagram is  attached to a  nilpotent element ``e`` of
the  Lie algebra  ``𝔤``. By  the Jacobson-Morozov  theorem there  exists an
``𝔰𝔩₂`` subalgebra of ``𝔤`` containing ``e`` as the element
``\\begin{pmatrix}1&0\\\\0&1  \\end{pmatrix}``.  Let  ``𝐒``  be  the  torus
``\\begin{pmatrix}h&0\\\\0&h^{-1} \\end{pmatrix}`` of ``SL₂`` and let ``𝐓``
be  a  maximal  torus  containing  ``𝐒``,  so  that ``𝐒`` is the image of a
one-parameter   subgroup  ``σ∈  Y(𝐓)``.  Consider  the  root  decomposition
``𝔤=∑_{α∈Σ}𝔤_α``  given by ``𝐓``; then ``α↦⟨σ,α⟩`` defines a linear form on
``Σ``,  determined by its value on simple roots. It is possible to choose a
system  of simple roots  ``Π`` so that  ``⟨σ,α⟩≥ 0`` for  ``α∈Π``, and then
``⟨σ,α⟩∈{0,1,2}`` for any ``α∈Π``. The Dynkin diagram of ``Π`` decorated by
these  values ``0,1,2`` is  called the Dynkin-Richardson  diagram of ``e``,
and in good characteristic is a complete invariant of its ``𝔤``-orbit.

Let  ``ℬ`` be  the variety  of all  Borel subgroups  and let  ``ℬᵤ`` be the
subvariety  of Borel subgroups  containing the unipotent  element `u`. Then
``dim C_𝐆(u)=rank 𝐆 + 2 dim ℬ_u`` and in good characteristic this dimension
can  be computed from  the Dynkin-Richardson diagram:  the dimension of the
class of `u` is the number of roots `α` such that ``⟨σ,α⟩∉{0,1}``.

We   describe  now  the  Springer  correspondence.  Indecomposable  locally
constant  ``𝐆``-equivariant sheaves  on ``C``,  called *local systems*, are
parameterized   by  irreducible  characters  of  ``A(u)``.  The  *ordinary*
Springer  correspondence is  a bijection  between irreducible characters of
the  Weyl group and a large subset  of the local systems which contains all
trivial  local  systems  (those  parameterized  by the trivial character of
``A(u)``  for  each  ``u``).  More  generally,  the  *generalized* Springer
correspondence   associates  to   each  local   system  a   (unique  up  to
``𝐆``-conjugacy)  *cuspidal pair* of  a Levi subgroup  ``𝐋`` of ``𝐆`` and a
*cuspidal*  local system on an unipotent class  of ``𝐋``, such that the set
of  local systems associated  to a given  cuspidal pair is parameterized by
the characters of the relative Weyl group ``W_𝐆 (𝐋):=N_𝐆 (𝐋)/𝐋``. There are
only few cuspidal pairs.

The  Springer correspondence gives information on the character values of a
finite  reductive groups  as follows:  assume that  ``k`` is  the algebraic
closure  of a finite field ``𝔽_q`` and that ``F`` is the Frobenius attached
to  an ``𝔽_q``-structure of  ``𝐆``. Let ``C``  be an ``F``-stable unipotent
class  and let ``u∈ C^F``; we call ``C`` the *geometric class* of ``u`` and
the ``𝐆^F``-classes inside ``C^F`` are parameterized by the ``F``-conjugacy
classes  of ``A(u)``, denoted ``H¹(F,A(u))`` (most  of the time we can find
``u`` such that ``F`` acts trivially on ``A(u)`` and ``H¹(F,A(u))`` is then
just the conjugacy classes). To an ``F``-stable character ``φ`` of ``A(u)``
we  associate  the  *characteristic  function*  of  the corresponding local
system (actually associated to an extension ``φ̃`` of ``φ`` to ``A(u).F``);
it  is a class function  ``Y_{u,φ}`` on ``𝐆^F`` which  can be normalized so
that:  ``Y_{u,φ}(u₁)=φ̃(cF)`` if ``u₁`` is geometrically conjugate to ``u``
and  its ``𝐆^F``-class is parameterized by the ``F``-conjugacy class ``cF``
of  ``A(u)``, otherwise ``Y_{u,φ}(u₁)=0``. If  the pair ``u,φ`` corresponds
via  the Springer correspondence to the character ``χ`` of ``W_𝐆(𝐋)``, then
``Y_{u,φ}``  is also  denoted ``Yᵪ``.  There is  another important class of
functions  indexed by local  systems: to a  local system on  class ``C`` is
attached  an intersection cohomology complex, which is a complex of sheaves
supported on the closure ``C̄``. To such a complex of sheaves is associated
its  *characteristic  function*,  a  class  function of ``𝐆^F`` obtained by
taking  the alternating trace of the Frobenius  acting on the stalks of the
cohomology  sheaves. If ``Y_ψ``  is the characteristic  function of a local
system,  the  characteristic  function  of  the  corresponding intersection
cohomology  complex is  denoted by  ``X_ψ``. This  function is supported on
``C̄``,  and Lusztig has shown that ``X_ψ=∑ᵩ P_{ψ,χ} Yᵪ`` where ``P_{ψ,χ}``
are  integer polynomials in ``q`` and  ``Yᵪ`` are attached to local systems
on classes lying in ``C̄``.

Lusztig   and  Shoji  have  given  an   algorithm  to  compute  the  matrix
``P_{ψ,χ}``,   which  is  implemented  in  Chevie.  The  relationship  with
characters   of  ``𝐆(𝔽_q)``,  taking  to  simplify  the  ordinary  Springer
correspondence,  is that the  restriction to the  unipotent elements of the
almost  character ``R_χ`` is equal to  ``q^{bᵪ} Xᵪ``, where ``bᵪ`` is ``dim
ℬᵤ``  for an element `u` of  the class `C` such that  the support of `χ` is
``C̄``.  The restriction of  the Deligne-Lusztig characters  ``R_w`` to the
unipotents  are called  the *Green  functions* and  can also be computed by
Chevie.  The values of  all unipotent characters  on unipotent elements can
also  be  computed  in  principle  by  applying Lusztig's Fourier transform
matrix (see the section on the Fourier matrix) but there is a difficulty in
that  the ``Xᵪ`` must be first multiplied  by some roots of unity which are
not  known in all cases (and when  known may depend on the congruence class
of ``q`` modulo some small primes).

We illustrate these computations on some examples:

```julia-repl
julia> UnipotentClasses(rootdatum(:sl,4))
UnipotentClasses(sl₄)
1111<211<22<31<4
   u│D-R dBu B-C          C(u) A₃(A₃₍₎=Φ₁³) A₁(A₃₍₁₃₎=A₁×A₁Φ₁)/-1 .(A₃)/ζ₄
────┼──────────────────────────────────────────────────────────────────────
4   │222   0 222         q³.Z₄          1:4                  -1:2    ζ₄:Id
31  │202   1 22.    q⁴.A₁₍₎=Φ₁        Id:31                               
22  │020   2 2.2      q⁴.A₁.Z₂         2:22                 11:11         
211 │101   3 2.. q⁵.A₂₍₁₎=A₁Φ₁       Id:211                               
1111│000   6 ...            A₃      Id:1111

   u│.(A₃)/ζ₄³
────┼──────────
4   │   ζ₄³:Id
31  │         
22  │         
211 │         
1111│         

```

The  first column in the table gives the name of the unipotent class, which
here  is  a  partition  describing  the  Jordan  form. The partial order on
unipotent  classes given by Zariski closure  is given before the table. The
column   `D-R`,   displayed   only   in   good  characteristic,  gives  the
Dynkin-Richardson  diagram  for  each  class;  the  column  `dBu` gives the
dimension  of the variety  ``ℬ ᵤ``. The  column `B-C` gives the Bala-Carter
classification  of ``u``, that is in the  case of ``sl₄`` it displays ``u``
as  a regular unipotent in a  Levi subgroup by giving the Dynkin-Richardson
diagram  of a regular  unipotent (all 2's)  at entries corresponding to the
Levi  and `.` at  entries which do  not correspond to  the Levi. The column
`C(u)`  describes the group  ``C_𝐆(u)``: a power  ``qᵈ`` describes that the
unipotent  radical of ``C_𝐆(u)`` has  dimension ``d`` (thus ``qᵈ`` rational
points);  then follows a  description of the  reductive part of the neutral
component  of ``C_𝐆(u)``,  given by  the name  of its  root datum.  Then if
``C_𝐆(u)``  is not  connected, the  description of  ``A(u)`` is given using
another  vocabulary: a  cyclic group  of order  4 is  given as  `Z4`, and a
symmetric group on 3 points would be given as `S3`.

For  instance, the first class `4`  has ``C_𝐆(u)^0`` unipotent of dimension
`3` and ``A(u)`` equal to `Z4`, the cyclic group of order 4. The class `22`
has  ``C_G(u)`` with unipotent radical of  dimension `4`, reductive part of
type  `A1` and ``A(u)`` is  `Z2`, that is the  cyclic group of order 2. The
other  classes have ``C_𝐆(u)`` connected. For  the class `31` the reductive
part of ``C_G(u)`` is a torus of rank 1.

Then  there is one column for each *Springer series*, giving for each class
the  pairs  'a:b'  where  'a'  is  the  name  of  the character of ``A(u)``
describing  the local system involved and 'b'  is the name of the character
of  the (relative) Weyl group corresponding by the Springer correspondence.
At  the top of the  column is written the  name of the relative Weyl group,
and  in brackets the  name of the  Levi affording a  cuspidal local system;
next,  separated  by  a  `/`  is  a  description  of  the central character
associated  to the  Springer series  (omitted if  this central character is
trivial):   all  local  systems  in  a  given  Springer  series  have  same
restriction  to the center of  ``𝐆``. To find what  the picture becomes for
another algebraic group in the same isogeny class, for instance the adjoint
group,  one  simply  discards  the  Springer series whose central character
becomes  trivial on the center of ``𝐆``;  and each group ``A(u)`` has to be
quotiented  by the common  kernel of the  remaining characters. Here is the
table for the adjoint group:

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
local system. Also, from the `B-C` column, we see that that class is not in
a  proper Levi,  in which  case the  Bala-Carter diagram coincides with the
Dynkin-Richardson diagram.

The  characteristics 2 and  3 are not  good for `G2`.  To get the unipotent
classes  and the Springer correspondence in bad characteristic, one gives a
second argument to the function `UnipotentClasses`:

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

The  function `ICCTable` gives the  transition matrix between the functions
``Xᵪ``  and ``Y_ψ``.

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
local systems: the row labels give the character of the relative Weyl group
and  the column labels give the class and the name of the local system as a
character  of `A(u)`: for instance, `G2(a1)` is the trivial local system of
the  class `G2(a1)`, while  `G2(a1)(21)` is the  local system on that class
corresponding to the 2-dimensional character of ``A(u)=A₂``.
"""
module Ucl

using ..Gapjm

export UnipotentClasses, UnipotentClassOps, ICCTable, XTable, GreenTable,
 UnipotentValues, induced_linear_form, special_pieces, name,
 distinguished_parabolics

@GapObj struct UnipotentClass
  name::String
  parameter::Any
  dimBu::Int
end

@GapObj struct UnipotentClasses
  classes::Vector{UnipotentClass}
  p::Int
  orderclasses::Poset
  springerseries::Vector{Dict}
end

function nameclass(u::Dict,opt=Dict{Symbol,Any}())
# println("u=$u, opt=$opt")
  if haskey(opt,:mizuno) && haskey(u,:mizuno) n=u[:mizuno]
  elseif haskey(opt,:shoji) && haskey(u,:shoji) n=u[:shoji]
  else n=u[:name]
  end
  TeX=haskey(opt,:TeX)
  n=fromTeX(n;opt...)
  if haskey(opt,:locsys) && opt[:locsys]!=charinfo(u[:Au]).positionId
    cl="("*charnames(u[:Au];opt...)[opt[:locsys]]*")"
    n*="^{$cl}"
    n=fromTeX(n;opt...)
  elseif haskey(opt,:class) && opt[:class]!=position_class(u[:Au],one(u[:Au]))
    cl=classinfo(u[:Au])[:classnames][opt[:class]]
    n=TeX ? "\\mbox{\$$n\$}_{($cl)}" : fromTeX("$(n)_{$cl}";opt...)
  end
  n
end

function name(io::IO,u::UnipotentClass)
 nameclass(merge(u.prop,Dict(:name=>u.name)),IOContext(io).dict)
end

name(u;opt...)=name(IOContext(stdout,opt...),u)

function Base.show(io::IO,u::UnipotentClass)
  print(io,"UnipotentClass(",name(io,u),")")
end

const UnipotentClassOps=Dict{Symbol,Any}(:Name=>nameclass)

"""
`induced_linear_form(W, K, h)`

This routine can be used to find the Richardson-Dynkin diagram of the class
in  the algebraic  group `𝐆`  which contains  a given  unipotent class of a
reductive subgroup of maximum rank `𝐒` of `𝐆`.

It  takes a linear  form on the  roots of `K`,  defined by its value on the
simple  roots (these values  can define a  Dynkin-Richardson diagram); then
extends  this linear form to  the roots of `𝐆`  by `0` on the orthogonal of
the  roots of `K`; and finally conjugates  the resulting form by an element
of the Weyl group so that it takes positive values on the simple roots.

```julia-repl
julia> W=coxgroup(:F,4)
F₄

julia> H=reflection_subgroup(W,[1,3])
F₄₍₁₃₎=A₁×Ã₁Φ₁²

julia> Ucl.induced_linear_form(W,H,[2,2])
4-element Vector{Int64}:
 0
 1
 0
 0

julia> uc=UnipotentClasses(W);

julia> uc.classes[4].prop
Dict{Symbol, Any} with 8 entries:
  :dynkin     => [0, 1, 0, 0]
  :dimred     => 6
  :red        => A₁×A₁
  :Au         => .
  :balacarter => [1, 3]
  :rep        => [1, 3]
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

"""
`distinguished_parabolics(W)` 

the  list  of  distinguished  parabolic  subgroups  of  `W` in the sense of
Richardson,  each  given  as  the  list  of  the corresponding indices. The
distinguished  unipotent  conjugacy  classes  of  `W`  consist of the dense
unipotent orbit in the unipotent radical of such a parabolic.

```julia-repl
julia> W=coxgroup(:F,4)
F₄

julia> distinguished_parabolics(W)
4-element Vector{Vector{Int64}}:
 []
 [3]
 [1, 3]
 [1, 3, 4]
```
"""
function distinguished_parabolics(W)
  filter(combinations(eachindex(gens(W)))) do J
    if isempty(J) return true end
    p=fill(1,semisimplerank(W))
    p[J]=fill(0,length(J))
    p=toM(W.rootdec[1:W.N])*p
    2*count(iszero,p)+semisimplerank(W)==count(isone,p)
  end
end

function BalaCarterLabels(W)
  vcat(map(parabolic_reps(W)) do J
    L=reflection_subgroup(W,J)
    map(distinguished_parabolics(L))do D
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
    ch(c)=map(j->ctu[c,findfirst(==(j),fusion)],1:HasType.NrConjugacyClasses(q))
    return Dict(:Au=>q,
      :chars=>map(c->findfirst(i->cth[i,:]==ch(c),axes(cth,1)),chars),
      :gens=>map(x->word(Au,HasType.First(elements(Au),y->h(y)==x)),gens(q)))
  end
  Z=n->ComplexReflectionGroup(n,1,1)
# println("Au=$Au chars=$chars")
  ct=transpose(CharTable(Au).irr[chars,:])
  cl=filter(i->ct[i,:]==ct[1,:],axes(ct,1))
# println("ct=$ct cl=$cl")
  if length(cl)==1 return Dict(:Au=>Au,:chars=>chars,
                              :gens=>map(x->[x],eachindex(gens(Au)))) end
  ct=transpose(toM(unique!(sort(toL(ct)))))
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

`orderclasses`:  a list describing the Hasse diagram of the partial order
induced   on   unipotent   classes   by   the  closure  relation.  That  is
`.orderclasses[i]`  is the list of `j` such that ``C̄ⱼ⊋ C̄ᵢ``  and  there  is
no  class  ``Cₖ``  such  that ``C̄ⱼ⊋ C̄ₖ⊋ C̄ᵢ``.

`classes`:  a  list  of  records  holding information for each unipotent
class (see below).

`springerseries`:  a list of records, each  of which describes a Springer
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
sl₄

julia> uc=UnipotentClasses(W);

julia> uc.classes
5-element Vector{Gapjm.Ucl.UnipotentClass}:
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

julia> xdisplay(uc;cols=[5,6,7],spaltenstein=true,frame=true,mizuno=true,
      order=false)
UnipotentClasses(E₆sc)
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
      n0=2*n0-count(==(2),weights)
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
#   s[:levi]=indices(t)[s[:levi]]
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
    classes=[UnipotentClass("1",[],0,
        Dict(:Au=>coxgroup(),:dynkin=>[],:balacarter=>[],
             :dimunip=>0,:red=>torus(rank(W))))]
    uc=[UnipotentClasses(classes,p,Poset([Int[]]),
      [Dict(:Z=>Int[],:levi=>Int[],:locsys=>[[1,1]],:relgroup=>coxgroup())],
      Dict{Symbol,Any}(:spets=>W))]
    l=Vector{Int}[]
  else
    classes=map(cartesian(map(x->x.classes,uc)...)) do v
      l=indices.(t)
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
  if !(p in badprimes(W)) && !haskey(classes[1],:balacarter)
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
    o=map(x->cart2lin(ll,x),o)
    setdiff(o,[cart2lin(ll,v)])
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
      [cart2lin(ll,v[1]),cart2lin(u,v[2])]
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
             algebraic_centre(W).descAZ),springerseries)
    AdjustAu!(classes,springerseries)
  end
  if spetscase
    g=Group(weightinfo(W)[:AdjointFundamentalGroup])
    permZ=map(x->word(g,x),gens(g).^WF.phi)
    springerseries=filter(x->map(i->prod(x[:Z][i]),permZ)==x[:Z],springerseries)
  end
# println(springerseries[1])
  s=springerseries[1]
  if spetscase
    s[:relgroup]=relative_coset(WF,s[:levi])
    s[:locsys]=s[:locsys][charinfo(s[:relgroup]).charRestrictions]
  end
  l=filter(i->any(y->i==y[1],s[:locsys]),1:length(classes))
  s[:locsys]=map(((c,s),)->[findfirst(==(c),l),s],s[:locsys])
  # for now only springerseries[1] properly twisted
  for s in springerseries[2:end]
    if spetscase
      s[:relgroup]=relative_coset(WF,s[:levi])
      s[:locsys]=s[:locsys][charinfo(s[:relgroup]).charRestrictions]
    end
    s[:locsys]=map(((c,s),)->[findfirst(==(c),l),s],s[:locsys])
  end
  classes=classes[l]
  AdjustAu!(classes,springerseries)
  orderclasses=Poset(hasse(restricted(Poset(orderclasses),l)),classes)
  orderclasses.show_element=function(io,x,n)
     print(io,name(io,x.elements[n]))
  end
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
        c*=repr(u.red;context=io)*AuName(u)
      elseif length(u.Au)==1 ||
         length(u.Au)==length(Group(u.AuAction.phis...))
        if length(u.Au)==1 || isone(u.red.phi)
          c*=repr(u.AuAction;context=io)
        else
          c*="["*repr(u.red;context=io)*"]"*repr(u.AuAction;context=io)
        end
      else
        c*=reflection_name(io,u.AuAction)*AuName(u)
      end
    else
      c*=AuName(u)
    end
  elseif haskey(u,:red)
    n=repr(u.red;context=io)
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
  show(IOContext(io,:TeX=>true),"text/plain",uc)
end

function Base.show(io::IO,uc::UnipotentClasses)
  print(io,"UnipotentClasses(",uc.spets,")")
end

function Base.show(io::IO,::MIME"text/plain",uc::UnipotentClasses)
  println(io,uc)
  TeX=get(io,:TeX,false)
  if get(io,:order,true) println(io,uc.orderclasses) end
  sp = map(copy, uc.springerseries)
  if get(io,:fourier,false)
    for p in sp p[:locsys] = p[:locsys][detPerm(p[:relgroup])] end
  end
  WF=uc.spets
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
        b=fill('.',ngens(W))
        for i in u.balacarter if i>0 b[i]='2' else b[-i]='0' end end
      else
        b=fill('?',ngens(W))
      end
      push!(res, String(b))
    end
    if get(io,:centralizer,true)
      push!(res,showcentralizer(io,u))
    end
    if get(io,:springer,true)
      i=findfirst(==(u),uc.classes)
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
        res=string(repr(ss[:relgroup];context=io),"(",
          repr(subspets(WF,ss[:levi]);context=io),")")
        if !all(isone,ss[:Z])
          res*=string("/", join(map(q->repr(q;context=io),ss[:Z]),","))
        end
        return res
    end, sp))
  end
  row_labels=name.(Ref(io),uc.classes)
  if get(io,:rows,false)==false
    io=IOContext(io,:rows=>sortperm(map(x->x.dimBu, uc.classes)))
  end
  showtable(io,toM(tbl);rows_label="u",col_labels=col_labels,row_labels=row_labels)
end

# decompose tensor product of characters (given as their indices in CharTable)
function decompose_tensor(W,c::Int...)
  ct=CharTable(W)
# println("eltype=",eltype(irr))
  decompose(ct,vec(prod(view(ct.irr,collect(c),:),dims=1)))
end

@GapObj struct ICCTable end

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
julia> xdisplay(t;cycpol=false)
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

`.dimBu`: The list of ``dimℬᵤ`` for each local system `(u,φ)` in the series.

`:L`:  The matrix of (unnormalized) scalar products of the functions ``Yᵩ``
with  themselves,  that  is  the  ``(φ,χ)``  entry  is ``∑_{g∈𝐆(𝔽_q)} Yᵩ(g)
Yᵪ(g)``. This is thus a symmetric, block-diagonal matrix where the diagonal
blocks  correspond to geometric unipotent conjugacy classes. This matrix is
obtained as a by-product of Lusztig's algorithm to compute ``Pᵩᵪ``.
"""
function ICCTable(uc::UnipotentClasses,i=1;q=Pol())
  W=uc.spets # W=Group(uc.spets)
  if W isa Spets W=W.W end
  ss=uc.springerseries[i]
  res=ICCTable(Dict(:spets=>uc.spets,:relgroup=>ss[:relgroup],
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
  k=charinfo(R).positionDet
# Partition on characters of ss.relgroup induced by poset of unipotent classes
  res.dimBu=map(x->uc.classes[x[1]].dimBu,ss[:locsys])
  res.blocks=HasType.CollectBy(eachindex(ss[:locsys]),-res.dimBu)
  subst=!(q isa Pol)
  if subst var=q; q=Pol() end
  f=fakedegrees(R,q)
  n=length(f)
  # matrix of q^{-bᵢ-bⱼ}*fakedegree(χᵢ ⊗ χⱼ ⊗ sgn)
  tbl=bigcell_decomposition([q^(-res.dimBu[i]-res.dimBu[j])*
                            sum(map(*,f,decompose_tensor(R,i,j,k)))
     for i in 1:n,j in 1:n], res.blocks) # //1 needed in D7
  if subst tbl=map(y->map(x->x(var),y),tbl); q=var end
  res.scalar=tbl[1]
  res.locsys=ss[:locsys]
# res.L=tbl[2]*GenericOrder(W,q)/prod(ReflectionDegrees(R),d->q^d-1)/
#   q^(W.semisimplerank-R.semisimplerank);
  res.L=tbl[2]*q^(nref(W)+semisimplerank(R)-semisimplerank(W))
  res.uc=uc
  if haskey(ss,:parameter) res.parameter=ss[:parameter]
  else res.parameter=(1:length(ss[:locsys])).+100*(i-1)
  end
  res.scalar=improve_type(res.scalar)
  res.L=improve_type(res.L)
  res
end

function Base.show(io::IO, ::MIME"text/html", x::ICCTable)
  show(IOContext(io,:TeX=>true),"text/plain",x)
end

Base.show(io::IO,x::ICCTable)=print(io,"ICCTable(",x.uc,",",x.series,")")

function Base.show(io::IO,::MIME"text/plain",x::ICCTable)
  printTeX(io,"Coefficients of \$X_\\chi\$ on \$Y_\\phi\$ for \$",x.relgroup,"\$\n")
  TeX=get(io,:TeX,false)
  if get(io,:cols,false)==false && get(io,:rows,false)==false
    rows=collect(eachindex(x.dimBu))
    sort!(rows,by=i->[x.dimBu[i],x.locsys[i]])
    io=IOContext(io,:rows=>rows,:cols=>rows)
  end
  tbl=get(io,:cycpol,true) ? map(CycPol,x.scalar) : x.scalar
  col_labels=map(((c,s),)->name(IOContext(io,:locsys=>s),x.uc.classes[c]),
                  x.locsys)
  rowLabels=map(x->TeX ? "X_{$x}" : "X$x",charnames(io,x.relgroup))
  showtable(io,transpose(tbl),row_labels=rowLabels,col_labels=col_labels)
end

@GapObj struct XTable end

"""
`XTable(uc;classes=false)`

This  function presents  in a  different way  the information obtained from
`ICCTable`.  Let  ``X̃_{u,ϕ}=q^{1/2(codim  C-dim  Z(𝐋 ))}`` where `C` is the
class  of `u` and `Z(𝐋 )` is the center of Levi subgroup on which lives the
cuspidal local system attached to the local system `(u,ϕ)`.

Then  `XTable` gives the decomposition of the functions ``X̃_{u,ϕ}`` on local
systems,  by  default.  If  `classes==true`,  it  gives  the  values of the
functions ``X̃_{u,ϕ}`` on unipotent classes.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> XTable(UnipotentClasses(W))
Values of character sheaves X̃ᵪ on local systems φ
      X̃ᵪ|φ│   1 A₁ Ã₁ G₂(a₁)⁽¹¹¹⁾ G₂(a₁)⁽²¹⁾ G₂(a₁) G₂
──────────┼────────────────────────────────────────────
X_φ₁‚₀^G₂ │   1  1  1           0          0      1  1
X_φ₁‚₆^G₂ │  q⁶  0  0           0          0      0  0
X_φ′₁‚₃^G₂│  q³  0  q           0          q      0  0
X_φ″₁‚₃^G₂│  q³ q³  0           0          0      0  0
X_φ₂‚₁^G₂ │ qΦ₈  q  q           0          0      q  0
X_φ₂‚₂^G₂ │q²Φ₄ q² q²           0          0      0  0
X_Id^.    │   0  0  0          q²          0      0  0
```

The functions `X̃` in the first column are decorated by putting as an
exponent the relative groups ``W_𝐆 (𝐋)``.

```julia-repl
julia> XTable(UnipotentClasses(W);classes=true)
Values of character sheaves X̃ᵪ on unipotent classes
  X̃ᵪ|class│   1 A₁ Ã₁ G₂(a₁) G₂(a₁)₍₂₁₎ G₂(a₁)₍₃₎ G₂
──────────┼──────────────────────────────────────────
X_φ₁‚₀^G₂ │   1  1  1      1          1         1  1
X_φ₁‚₆^G₂ │  q⁶  0  0      0          0         0  0
X_φ′₁‚₃^G₂│  q³  0  q     2q          0        -q  0
X_φ″₁‚₃^G₂│  q³ q³  0      0          0         0  0
X_φ₂‚₁^G₂ │ qΦ₈  q  q      q          q         q  0
X_φ₂‚₂^G₂ │q²Φ₄ q² q²      0          0         0  0
X_Id^.    │   0  0  0     q²        -q²        q²  0

julia> XTable(UnipotentClasses(W,2))
Values of character sheaves X̃ᵪ on local systems φ
      X̃ᵪ|φ│   1 A₁ Ã₁ G₂(a₁)⁽¹¹¹⁾ G₂(a₁)⁽²¹⁾ G₂(a₁) G₂⁽¹¹⁾ G₂
──────────┼───────────────────────────────────────────────────
X_φ₁‚₀^G₂ │   1  1  1           0          0      1      0  1
X_φ₁‚₆^G₂ │  q⁶  0  0           0          0      0      0  0
X_φ′₁‚₃^G₂│  q³  0  q           0          q      0      0  0
X_φ″₁‚₃^G₂│  q³ q³  0           0          0      0      0  0
X_φ₂‚₁^G₂ │ qΦ₈  q  q           0          0      q      0  0
X_φ₂‚₂^G₂ │q²Φ₄ q² q²           0          0      0      0  0
X_Id^.    │   0  0  0          q²          0      0      0  0
X_Id^.    │   0  0  0           0          0      0      q  0

julia> XTable(UnipotentClasses(rootdatum(:sl,4)))
Values of character sheaves X̃ᵪ on local systems φ
    X̃ᵪ|φ│1111 211 22⁽¹¹⁾ 22 31 4 4^(ζ₄) 4⁽⁻¹⁾ 4^(ζ₄³)
────────┼─────────────────────────────────────────────
X₁₁₁₁^A₃│  q⁶   0      0  0  0 0      0     0       0
X₂₁₁^A₃ │q³Φ₃  q³      0  0  0 0      0     0       0
X₂₂^A₃  │q²Φ₄  q²      0 q²  0 0      0     0       0
X₃₁^A₃  │ qΦ₃ qΦ₂      0  q  q 0      0     0       0
X₄^A₃   │   1   1      0  1  1 1      0     0       0
X₁₁^A₁  │   0   0     q³  0  0 0      0     0       0
X₂^A₁   │   0   0     q²  0  0 0      0     q       0
X_Id^.  │   0   0      0  0  0 0   q³⁄₂     0       0
X_Id^.  │   0   0      0  0  0 0      0     0    q³⁄₂
```

A  side effect  of calling  `XTable` with  `classes=true` is to compute the
cardinal of the unipotent conjugacy classes:

```julia-repl
julia> t=Ucl.XTable(UnipotentClasses(coxgroup(:G,2));classes=true);

julia> CycPol.(t.cardClass)
7-element Vector{CycPol{Cyc{Rational{Int64}}}}:
 1
 Φ₁Φ₂Φ₃Φ₆
 q²Φ₁Φ₂Φ₃Φ₆
 q²Φ₁²Φ₂²Φ₃Φ₆/6
 q²Φ₁²Φ₂²Φ₃Φ₆/2
 q²Φ₁²Φ₂²Φ₃Φ₆/3
 q⁴Φ₁²Φ₂²Φ₃Φ₆
```
"""
function XTable(uc::UnipotentClasses;q=Mvp(:q),classes=false)
# println("here uc=",uc)
  pieces=map(i->ICCTable(uc,i),eachindex(uc.springerseries))
# Note that c_ι=βᵤ+(rkss L_\CI)/2
  greenpieces=map((x,y)->map(x->x(q),x.scalar)*
                  toM(HasType.DiagonalMat(q.^x.dimBu...))*
                  q^(length(y[:levi])//2),pieces,uc.springerseries)
  l=vcat(getproperty.(pieces,:locsys)...)
  p=inv(sortPerm(l))
  res=XTable(Dict(
    :scalar=>transpose(cat(greenpieces...,dims=(1,2))^p),
    :uc=>uc,
    :Y=>^(cat(getproperty.(pieces,:L)...,dims=(1,2)),p,dims=(1,2)),
    :relgroups=>getindex.(uc.springerseries,:relgroup),
    :q=>q,
    :class=>classes))
  res.Y=map(x->x(q),res.Y)
  if classes
    res.scalar*=E(1)
    res.cardClass=zeros(eltype(res.scalar),length(l))*1//1
    res.classes=l^p
    for i in eachindex(uc.classes)
      Au=uc.classes[i].Au
      b=filter(j->res.classes[j][1]==i,eachindex(res.classes))
 #    println("i=",i," b=",b," Au=",Au)
      res.scalar[:,b]*=CharTable(Au).irr
      res.cardClass[b]=res.Y[[b[charinfo(Au).positionId]],b]*CharTable(Au).irr
      res.cardClass[b]=map((x,y)->x*y//length(Au),
                             res.cardClass[b],classinfo(Au)[:classes])
    end
    res.scalar=improve_type(res.scalar)
  else
    res.locsys=l^p
  end
  res
end

function Base.show(io::IO, ::MIME"text/html", x::XTable)
  show(IOContext(io,:TeX=>true),"text/plain",x)
end

Base.show(io::IO,x::XTable)=print(io,"XTable(",x.uc,",q=",x.q,",classes=$(x.class))")

function Base.show(io::IO,::MIME"text/plain",x::XTable)
  printTeX(io,"Values of character sheaves \$\\tilde X_\\chi\$ on")
  row_labels=vcat(map(g->map(n->"X_{"*n*"}^{"*TeX(io,g)*"}",
                             charnames(TeX(io),g)),x.relgroups)...)
  rows_label="\\tilde X_\\chi|"
  if x.class
    rows_label*="class"
    print(io," unipotent classes\n")
    col_labels=map(p->name(TeX(io;class=p[2]),x.uc.classes[p[1]]),x.classes)
  else 
    rows_label*=fromTeX(io,"\\phi")
    printTeX(io," local systems \$\\phi\$\n")
    col_labels=map(p->name(TeX(io,locsys=p[2]),x.uc.classes[p[1]]),x.locsys)
  end
  tbl=x.scalar
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  showtable(io,tbl;row_labels,col_labels,rows_label)
end

@GapObj struct GreenTable end

"""
`GreenTable(uc;classes=false)`

Keeping the same notations as in the description of 'XTable', this function
returns a table of the functions ``Q_{wF}``, attached to elements ``wF∈ W_𝐆
(𝐋)⋅F`` where ``W_𝐆 (𝐋)`` are the relative weyl groups attached to cuspidal
local  systems.  These  functions  are  defined  by ``Q_{wF}=∑_{u,ϕ} ϕ̃(wF)
X̃_{u,ϕ}``. An point to note is that in the principal Springer series, when
`𝐓`  is  a  maximal  torus,  the  function  ``Q_{wF}``  coincides  with the
Deligne-Lusztig  character ``R^𝐆  _{𝐓_W}(1)``. As  for 'XTable', by default
the  table  gives  the  values  of  the  functions  on  local  systems.  If
`classes=true`  is  given,  then  it  gives  the  values  of  the functions
``Q_{wF}`` on conjugacy classes.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> GreenTable(UnipotentClasses(W))
Values of Green functions Q_wF on local systems φ
   Qᴵ_wF|φ│        1     A₁       Ã₁ G₂(a₁)⁽¹¹¹⁾ G₂(a₁)⁽²¹⁾ G₂(a₁) G₂
──────────┼───────────────────────────────────────────────────────────
Q_A₀^G₂   │  Φ₂²Φ₃Φ₆   Φ₂Φ₃ (2q+1)Φ₂           0          q   2q+1  1
Q_Ã₁^G₂   │-Φ₁Φ₂Φ₃Φ₆  -Φ₁Φ₃       Φ₂           0          q      1  1
Q_A₁^G₂   │-Φ₁Φ₂Φ₃Φ₆   Φ₂Φ₆      -Φ₁           0         -q      1  1
Q_G₂^G₂   │ Φ₁²Φ₂²Φ₃ -Φ₁Φ₂²    -Φ₁Φ₂           0         -q     Φ₂  1
Q_A₂^G₂   │ Φ₁²Φ₂²Φ₆  Φ₁²Φ₂    -Φ₁Φ₂           0          q    -Φ₁  1
Q_A₁+Ã₁^G₂│  Φ₁²Φ₃Φ₆  -Φ₁Φ₆ (2q-1)Φ₁           0         -q  -2q+1  1
Q_^.      │        0      0        0          q²          0      0  0
```

The  functions ``Q_{wF}`` depend only on the conjugacy class of `wF`, so in
the  first column the indices of 'Q' are the names of the conjugacy classes
of ``W_𝐆(𝐋)``. The exponents are the names of the groups ``W_𝐆(𝐋)``.

```julia-repl
julia> GreenTable(UnipotentClasses(W);classes=true)
Values of Green functions Q_wF on unipotent classes
Qᴵ_wF|class│        1     A₁       Ã₁ G₂(a₁) G₂(a₁)₍₂₁₎ G₂(a₁)₍₃₎ G₂
───────────┼─────────────────────────────────────────────────────────
Q_A₀^G₂    │  Φ₂²Φ₃Φ₆   Φ₂Φ₃ (2q+1)Φ₂   4q+1       2q+1        Φ₂  1
Q_Ã₁^G₂    │-Φ₁Φ₂Φ₃Φ₆  -Φ₁Φ₃       Φ₂   2q+1          1       -Φ₁  1
Q_A₁^G₂    │-Φ₁Φ₂Φ₃Φ₆   Φ₂Φ₆      -Φ₁  -2q+1          1        Φ₂  1
Q_G₂^G₂    │ Φ₁²Φ₂²Φ₃ -Φ₁Φ₂²    -Φ₁Φ₂    -Φ₁         Φ₂      2q+1  1
Q_A₂^G₂    │ Φ₁²Φ₂²Φ₆  Φ₁²Φ₂    -Φ₁Φ₂     Φ₂        -Φ₁     -2q+1  1
Q_A₁+Ã₁^G₂ │  Φ₁²Φ₃Φ₆  -Φ₁Φ₆ (2q-1)Φ₁  -4q+1      -2q+1       -Φ₁  1
Q_^.       │        0      0        0     q²        -q²        q²  0

julia> GreenTable(UnipotentClasses(rootdatum(:sl,4)))
Values of Green functions Q_wF on local systems φ
 Qᴵ_wF|φ│     1111          211 22⁽¹¹⁾       22   31 4 4^(ζ₄) 4⁽⁻¹⁾ 4^(ζ₄³)
────────┼───────────────────────────────────────────────────────────────────
Q₁₁₁₁^A₃│  Φ₂²Φ₃Φ₄ (3q²+2q+1)Φ₂      0 (2q+1)Φ₂ 3q+1 1      0     0       0
Q₂₁₁^A₃ │-Φ₁Φ₂Φ₃Φ₄   -q³+q²+q+1      0       Φ₂   Φ₂ 1      0     0       0
Q₂₂^A₃  │  Φ₁²Φ₃Φ₄        -Φ₁Φ₄      0  2q²-q+1  -Φ₁ 1      0     0       0
Q₃₁^A₃  │ Φ₁²Φ₂²Φ₄        -Φ₁Φ₂      0    -Φ₁Φ₂    1 1      0     0       0
Q₄^A₃   │ -Φ₁³Φ₂Φ₃        Φ₁²Φ₂      0      -Φ₁  -Φ₁ 1      0     0       0
Q₁₁^A₁  │        0            0   q²Φ₂        0    0 0      0     q       0
Q₂^A₁   │        0            0  -q²Φ₁        0    0 0      0     q       0
Q_^.    │        0            0      0        0    0 0   q³⁄₂     0       0
Q_^.    │        0            0      0        0    0 0      0     0    q³⁄₂
```
"""
function GreenTable(uc::UnipotentClasses;q=Mvp(:q),classes=false)
  t=GreenTable(XTable(uc;classes,q).prop)
  m=cat(map(g->transpose(CharTable(g).irr),t.relgroups)...;dims=(1,2))
  t.scalar=m*t.scalar
  t.indices=Vector{Int}[]
  i=0
  for g in t.relgroups
    n=nconjugacy_classes(g)
    push!(t.indices,(1:n)+i)
    i+=n
  end
  t
end

function Base.show(io::IO, ::MIME"text/html", x::GreenTable)
  show(IOContext(io,:TeX=>true),"text/plain",x)
end

Base.show(io::IO,x::GreenTable)=print(io,"GreenTable(",x.uc,",q=",x.q,")")

function Base.show(io::IO,::MIME"text/plain",x::GreenTable)
  printTeX(io,"Values of Green functions \$Q_{wF}\$ on")
  row_labels=vcat(map(x.relgroups) do g
     map(n->string("Q_{",n,"}^{",TeX(io,g),"}"),classinfo(g)[:classnames])
    end...)
  rows_label="Q^I_{wF}|"
  if x.class
    println(io," unipotent classes")
    rows_label*="class"
    col_labels=map(p->name(TeX(io;class=p[2]),x.uc.classes[p[1]]),x.classes)
  else
    rows_label*="\\varphi"
    printTeX(io," local systems \$\\varphi\$\n")
    col_labels=map(p->name(TeX(io;locsys=p[2]),x.uc.classes[p[1]]),x.locsys)
  end
  tbl=x.scalar
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  showtable(io,tbl;row_labels,col_labels,rows_label)
end

@GapObj struct ValuesTable end

# values of unipotent characters
# UnipotentValues(uc[,opt]) values on unipotent classes (opt.classes bound)
# or local systems (opt.classes unbound)
"""
`UnipotentValues(uc,classes=false)`

This  function returns  a table  of the  values of  unipotent characters on
local systems (by default) or on unipotent classes (if `classes=true`).

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> UnipotentValues(UnipotentClasses(W);classes=true)
Values of unipotent characters for G₂ on unipotent classes
       │        1          A₁     Ã₁   G₂(a₁) G₂(a₁)₍₂₁₎ G₂(a₁)₍₃₎ G₂
───────┼──────────────────────────────────────────────────────────────
φ₁‚₀   │        1           1      1        1          1         1  1
φ₁‚₆   │       q⁶           0      0        0          0         0  0
φ′₁‚₃  │  qΦ₃Φ₆/3    -qΦ₁Φ₂/3      q (q+5)q/3     -qΦ₁/3     qΦ₁/3  0
φ″₁‚₃  │  qΦ₃Φ₆/3  (2q²+1)q/3      0    qΦ₁/3     -qΦ₁/3  (q+2)q/3  0
φ₂‚₁   │ qΦ₂²Φ₃/6 (2q+1)qΦ₂/6  qΦ₂/2 (q+5)q/6     -qΦ₁/6     qΦ₁/6  0
φ₂‚₂   │ qΦ₂²Φ₆/2       qΦ₂/2  qΦ₂/2   -qΦ₁/2      qΦ₂/2    -qΦ₁/2  0
G₂[-1] │ qΦ₁²Φ₃/2      -qΦ₁/2 -qΦ₁/2   -qΦ₁/2      qΦ₂/2    -qΦ₁/2  0
G₂[1]  │ qΦ₁²Φ₆/6 (2q-1)qΦ₁/6 -qΦ₁/2 (q+5)q/6     -qΦ₁/6     qΦ₁/6  0
G₂[ζ₃] │qΦ₁²Φ₂²/3    -qΦ₁Φ₂/3      0    qΦ₁/3     -qΦ₁/3  (q+2)q/3  0
G₂[ζ₃²]│qΦ₁²Φ₂²/3    -qΦ₁Φ₂/3      0    qΦ₁/3     -qΦ₁/3  (q+2)q/3  0


julia> UnipotentValues(UnipotentClasses(W,3);classes=true)
Values of unipotent characters for G₂ on unipotent classes
       │        1          A₁         Ã₁ G₂(a₁) G₂(a₁)₍₂₎    G₂       G₂_(ζ₃)
───────┼──────────────────────────────────────────────────────────────────────
φ₁‚₀   │        1           1          1      1         1     1             1
φ₁‚₆   │       q⁶           0          0      0         0     0             0
φ′₁‚₃  │  qΦ₃Φ₆/3    -qΦ₁Φ₂/3        q/3  qΦ₂/3    -qΦ₁/3 -2q/3           q/3
φ″₁‚₃  │  qΦ₃Φ₆/3  (2q²+1)q/3        q/3  qΦ₂/3    -qΦ₁/3 -2q/3           q/3
φ₂‚₁   │ qΦ₂²Φ₃/6 (2q+1)qΦ₂/6  (3q+1)q/6  qΦ₂/6    -qΦ₁/6  2q/3          -q/3
φ₂‚₂   │ qΦ₂²Φ₆/2       qΦ₂/2      qΦ₂/2 -qΦ₁/2     qΦ₂/2     0             0
G₂[-1] │ qΦ₁²Φ₃/2      -qΦ₁/2     -qΦ₁/2 -qΦ₁/2     qΦ₂/2     0             0
G₂[1]  │ qΦ₁²Φ₆/6 (2q-1)qΦ₁/6 (-3q+1)q/6  qΦ₂/6    -qΦ₁/6  2q/3          -q/3
G₂[ζ₃] │qΦ₁²Φ₂²/3    -qΦ₁Φ₂/3        q/3  qΦ₂/3    -qΦ₁/3   q/3 (-ζ₃+2ζ₃²)q/3
G₂[ζ₃²]│qΦ₁²Φ₂²/3    -qΦ₁Φ₂/3        q/3  qΦ₂/3    -qΦ₁/3   q/3  (2ζ₃-ζ₃²)q/3

       │     G₂_(ζ₃²)       (Ã₁)₃
───────┼──────────────────────────
φ₁‚₀   │            1           1
φ₁‚₆   │            0           0
φ′₁‚₃  │          q/3  (2q²+1)q/3
φ″₁‚₃  │          q/3    -qΦ₁Φ₂/3
φ₂‚₁   │         -q/3 (2q+1)qΦ₂/6
φ₂‚₂   │            0       qΦ₂/2
G₂[-1] │            0      -qΦ₁/2
G₂[1]  │         -q/3 (2q-1)qΦ₁/6
G₂[ζ₃] │ (2ζ₃-ζ₃²)q/3    -qΦ₁Φ₂/3
G₂[ζ₃²]│(-ζ₃+2ζ₃²)q/3    -qΦ₁Φ₂/3
```
"""
function UnipotentValues(uc;q=Mvp(:q),classes=false)
  t=ValuesTable(XTable(uc;classes,q).prop)
  uw=UnipotentCharacters(uc.spets)
  f=toL(fourier(uw))
  m=Vector{eltype(f[1])}[]
  for (i,ss) in enumerate(uc.springerseries)
    if i==1 append!(m,f[uw.harishChandra[1][:charNumbers]])
    elseif !haskey(ss,:hc) error("not implemented") 
    elseif ss[:hc]==0 append!(m,map(i->zero(f[1]),eachindex(ss[:locsys])))
    else append!(m,f[uw.harishChandra[ss[:hc]][:charNumbers]])
    end
  end
  t.scalar=transpose(toM(m))*t.scalar
  t
end

Base.show(io::IO,x::ValuesTable)=print(io,"UnipotentValues(",x.uc,",q=",x.q,")")

function Base.show(io::IO,::MIME"text/plain",x::ValuesTable)
  printTeX(io,"Values of unipotent characters for \$",x.uc.spets,"\$ on ")
  if x.class
    println(io,"unipotent classes")
    col_labels=map(p->name(TeX(io;class=p[2]),x.uc.classes[p[1]]),x.classes)
  else 
    col_labels=map(p->name(TeX(io;locsys=p[2]),x.uc.classes[p[1]]),x.locsys)
    println(io,"local systems")
  end
  row_labels=charnames(TeX(io),UnipotentCharacters(x.uc.spets))
  tbl=improve_type(x.scalar)
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  showtable(io,tbl;row_labels,col_labels)
end

@GapObj struct TwoVarGreenTable end

# two-variable green functions table.
# for now only implemented when W is split.
function TwoVarGreen(W,L)
  if !(W isa Spets) W=spets(W) end
  if !(L isa Spets) L=spets(L) end
  uG=UnipotentClasses(W)
  uL=UnipotentClasses(L)
  tG=GreenTable(uG)
  tL=GreenTable(uL)
  q=Mvp(:q)
  mm=map(eachindex(uL.springerseries))do i
    s=uL.springerseries[i]
    p=findfirst(S->S[:levi]==inclusion(L,s[:levi]) &&
                (isempty(S[:Z]) || S[:Z][1]==s[:Z][1]),uG.springerseries)
    if isnothing(p) error("not found ",s[:levi]) end
    RG=relative_coset(W,inclusion(L,s[:levi]))
    RLF=relative_coset(L,s[:levi])
    RL=Group(RLF)
    l=map(x->findfirst(==(x),Group(RG).relativeIndices),
      inclusion(L,RL.relativeIndices))
    if nothing in l error("not implemented") end
    RLF=subspets(RG,convert(Vector{Int},l),
                Group(RG).MappingFromNormalizer(L.phi))
    RL=Group(RLF)
    f=fusion_conjugacy_classes(RLF,RG)
    cl=classreps(RLF)
    d=map(cl)do w
      pw=word(RLF,w)
      if isempty(pw) pw=L.phi
      else pw=prod(Group(RG).parentMap[pw])*L.phi
      end
      Lo=subspets(L,Int.(s[:levi]),pw/L.phi)
      r=map(last,filter(x->isone(first(x)),degrees(Lo)))
      prod(x->q-x,r)/length(centralizer(RL,w))
    end
    transpose(tL.scalar[tL.indices[i],:])*
    toM(HasType.DiagonalMat(d...))*conj(tG.scalar[tG.indices[p][f],:])
  end
  oL=generic_order(L,q)
  mm=improve_type(mm)
  mm=toM(map((x,y)->exactdiv(x*y,oL),eachrow(sum(mm)),tL.cardClass))
  res=TwoVarGreenTable(Dict(:W=>W,:L=>L,:scalar=>mm,:uL=>uL,:uG=>uG))
  res.classL=tL.classes
  res.classG=tG.classes
  res
end

function Base.show(io::IO, ::MIME"text/html", x::TwoVarGreenTable)
  show(IOContext(io,:TeX=>true),"text/plain",x)
end

Base.show(io::IO,x::TwoVarGreenTable)=print(io,"TwoVarGreenTable(",x.W,",",x.L,")")

function Base.show(io::IO,::MIME"text/plain",x::TwoVarGreenTable)
  print(io,"Values of two-variable Green functions for \$",x.W,"\$ and \$",
        x.L,"\$\n")
  rowLabels=map(x.classL)do p
    name(IOContext(io,:class=>p[2]),x.uL.classes[p[1]])
  end
  columnLabels=map(x.classG)do p
    name(IOContext(io,:class=>p[2]),x.uG.classes[p[1]])
  end
  tbl=improve_type(x.scalar)
  if get(io,:cycpol,true) tbl=CycPol.(tbl) end
  showtable(io,tbl,row_labels=rowLabels,col_labels=columnLabels,rows_label="v/u")
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
3-element Vector{Vector{Int64}}:
 [1]
 [4, 3, 2]
 [5]

julia> special_pieces(UnipotentClasses(W,3))
3-element Vector{Vector{Int64}}:
 [1]
 [4, 3, 2, 6]
 [5]
```

The   example  above  shows  that  the  special  pieces  are  different  in
characteristic 3.
"""
function special_pieces(uc)
  W=uc.spets
  ch=charinfo(W)
  specialch=findall(iszero,ch.a-ch.b) # special characters of W
  specialc=first.(uc.springerseries[1][:locsys][specialch])
  sort!(specialc,by=c->-uc.classes[c].dimBu)
  m=transpose(incidence(uc.orderclasses))
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
