"""
This  module contains  functions for  computing with  unipotent elements of
reductive  groups; specifically to  compute with elements  of the unipotent
radical of a Borel subgroup of a connected algebraic reductive group; these
functions  were initially written by Olivier Dudas in GAP3, partly inspired
by older C code of Jean Michel.

The  unipotent radical of a  Borel subgroup is the  product in any order of
the root subgroups associated to the positive roots. We fix an order, which
gives a canonical form to display elements and to compare them.

The  computations use the Steinberg relations between root subgroups, which
come from the choice of a Chevalley basis of the Lie algebra. The reference
we  follow is [Carter1972,  chapters 4 to  6](biblio.htm#Car72b), but it is
possible  to use another  choice of Chevalley  basis, see the documentation
for `UnipotentGroup`.

We  start with  a root  datum specified  by a  Weyl group  `W` and  build a
`struct  UnipotentGroup`  which  contains  information  about  the  maximal
unipotent  subgroup  of  the  corresponding  reductive  group,  that is the
unipotent radical of the Borel subgroup determined by the positive roots.

```julia-repl
julia> W=coxgroup(:E,6)
E₆

julia> U=UnipotentGroup(W)
UnipotentGroup(E₆)
```

Now, if `α=roots(W,2)`, we construct the element `u_α(4)`
of the root subgroup `u_α`:

```julia-repl
julia> U(2=>4)
u2(4)
```

If we do not specify the coefficient we make by default `u_α(1)`, so we
have also:

```julia-repl
julia> U(2)^4
u2(4)
```

We can make more complicated elements:

```julia-repl
julia> U(2=>4)*U(4=>5)
u2(4)u4(5)

julia> U(2=>4,4=>5)
u2(4)u4(5)
```

If the roots are not in order the element is normalized:

```julia-repl
julia> U(4=>5,2=>4)
u2(4)u4(5)u8(-20)
```

It is possible to display the decomposition of the roots in simple roots
instead of their index:

```julia-rep1
julia> xdisplay(U(4=>5,2=>4);root=true)
u₀₁₀₀₀₀(4)u₀₀₀₁₀₀(5)u₀₁₀₁₀₀(-20)
```

The coefficients in the root subgroups can be elements of arbitrary rings.
Here is an example using `Mvp`'s:

```julia-repl
julia> W=coxgroup(:E,8);U=UnipotentGroup(W)
UnipotentGroup(E₈)

julia> u=U(map(i->i=>Z(2)*Mvp(Symbol("x",Char(i+0x2080))),1:8)...)
u1(x₁)u2(x₂)u3(x₃)u4(x₄)u5(x₅)u6(x₆)u7(x₇)u8(x₈)
```

```julia-rep1
julia> cut(xrepr(u^16;limit=true,root=true),before="u",width=60)
u₂₂₃₄₃₂₁₀(x₁²x₂²x₃³x₄⁴x₅³x₆²x₇)
u₁₂₃₄₃₂₁₁(x₁x₂²x₃³x₄⁴x₅³x₆²x₇x₈)
u₁₂₂₄₃₂₂₁(x₁x₂²x₃²x₄⁴x₅³x₆²x₇²x₈)
u₁₂₂₃₃₃₂₁(x₁x₂²x₃²x₄³x₅³x₆³x₇²x₈)
u₂₂₃₄₃₂₁₁(x₁²x₂²x₃³x₄⁴x₅³x₆²x₇x₈)
u₁₂₂₄₃₃₂₁(x₁x₂²x₃²x₄⁴x₅³x₆³x₇²x₈)
u₁₂₂₄₄₃₂₁(x₁x₂²x₃²x₄⁴x₅⁴x₆³x₇²x₈)
u₂₂₃₄₃₃₂₁(x₁²x₂²x₃³x₄⁴x₅³x₆³x₇²x₈)
u₁₂₃₄₄₃₂₁(x₁x₂²x₃³x₄⁴x₅⁴x₆³x₇²x₈)
u₂₂₃₄₄₃₂₁(x₁²x₂²x₃³x₄⁴x₅⁴x₆³x₇²x₈)
u₂₃₃₅₄₃₂₁(x₁²x₂³x₃³x₄⁵x₅⁴x₆³x₇²x₈)
u₂₂₄₅₄₃₂₁(x₁²x₂²x₃⁴x₄⁵x₅⁴x₆³x₇²x₈)
u₂₃₄₆₅₄₃₂(x₁²x₂³x₃⁴x₄⁶x₅⁵x₆⁴x₇³x₈²)
```

```julia-repl
julia> u^32
()
```

The  above computation shows  that in characteristic  2 the exponent of the
unipotent  group of `E₈` is 32. More precisely, squaring doubles the height
of  the involved roots, so in the above `u¹⁶` involves only roots of height
16 or more.

Various  actions are  defined on  unipotent elements.  Elements of the Weyl
group  act (through certain representatives) as long as no root subgroup is
in their inversion set:

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> U=UnipotentGroup(W);@Mvp x,y

julia> u=U(1=>x,3=>y)
u1(x)u3(y)

julia> u^W(2,1)
u4(y)u5(x)
```

```julia-repl
julia> u^W(1)
ERROR: u1(x)u3(y) should have no coefficient on root 1
```

Semisimple elements act by conjugation:

```julia-repl
julia> s=SemisimpleElement(W,[E(3),2])
SemisimpleElement{Cyc{Int64}}: <ζ₃,2>

julia> u^s
u1(ζ₃x)u3(2ζ₃y)
```

As well as unipotent elements:

```julia-repl
julia> u^U(2)
u1(x)u3(x+y)u4(-x-2y)u5(x+3y)u6(x²+3xy+3y²)
```
"""
module Urad
using ..Chevie
export UnipotentElement, UnipotentGroup, reorder, abelianpart

@GapObj struct UnipotentGroup
  W
  order::Vector{Int} # total order on roots used to normalize unipotents
  special::Vector{NamedTuple{(:r,:s,:rs,:N,:comm),
                             Tuple{Int,Int,Int,Int,Vector{NTuple{4,Int}}}}}
end

@doc """
A  `struct UnipotentGroup` represents  the unipotent  radical `𝐔` of a
Borel subgroup of a reductive group `G`.

See   [Carter1972,  section  4.2](biblio.htm#Car72b)  for  details  on  the
following.  A Chevalley basis  of the Lie  algebra of `𝐔`  is a basis `eᵣ`,
where   each  `eᵣ`  is  in  the  corresponding  root  subspace,  such  that
`[eᵣ,eₛ]=Nᵣₛ e_{r+s}` for some integer constants `Nᵣₛ`. The constants `Nᵣₛ`
for general roots are computed from the case where `r` and `s` are positive
roots whose sum is a root; such a pair `(r,s)` is called *special*.

Constants  `Cᵣₛᵢⱼ` are defined, see [Carter1972, 5.2.3](biblio.htm#Car72b),
by

``u_s(u)u_r(t)=u_r(t)u_s(u)\\prod_{i,j>0}u_{ir+js}(C_{rsij}(-t)^iu^j)``

Where  `ir+js` runs over the positive  integral combinations of `r` and `s`
which  are roots,  taken in  lexicographic order  on `(i,j)`. The constants
`Cᵣₛᵢⱼ`  are computed from the constants  `Nᵣₛ`, see [Carter1972, bottom of
page 61 and top of page 76](biblio.htm#Car72b).

The fields of `struct Unipotent Group` are:

  * `W`:         the Weyl group of `G`
  * `order`:     the total order on the roots used to normalize products of root subgroups (by default `1:nref(W)`)
  * `special:    a `NamedTuple` for each special pair of roots `(r,s)`
  - `r`  index of `r`
  - `s`  index of `s`
  - `rs` index of `r+s`
  - `N`: the constant `Nᵣₛ`
  - `comm`: stores the `Cᵣₛᵢⱼ` as the list of quadruples `(i,j,ir+js,Cᵣₛᵢⱼ)`.
```julia-repl
julia> U=UnipotentGroup(coxgroup(:G,2))
UnipotentGroup(G₂)

julia> U.special
10-element Vector{@NamedTuple{r::Int64, s::Int64, rs::Int64, N::Int64, comm::Vector{NTuple{4, Int64}}}}:
 (r = 1, s = 2, rs = 3, N = 1, comm = [(1, 1, 3, 1), (1, 2, 4, -1), (1, 3, 5, 1), (2, 3, 6, 2)])
 (r = 2, s = 3, rs = 4, N = 2, comm = [(1, 1, 4, 2), (2, 1, 5, 3), (1, 2, 6, -3)])
 (r = 2, s = 4, rs = 5, N = 3, comm = [(1, 1, 5, 3)])
 (r = 1, s = 5, rs = 6, N = 1, comm = [(1, 1, 6, 1)])
 (r = 3, s = 4, rs = 6, N = 3, comm = [(1, 1, 6, 3)])
 (r = 2, s = 1, rs = 3, N = -1, comm = [(1, 1, 3, -1), (2, 1, 4, -1), (3, 1, 5, -1), (3, 2, 6, -1)])
 (r = 3, s = 2, rs = 4, N = -2, comm = [(1, 1, 4, -2), (2, 1, 6, -3), (1, 2, 5, 3)])
 (r = 4, s = 2, rs = 5, N = -3, comm = [(1, 1, 5, -3)])
 (r = 5, s = 1, rs = 6, N = -1, comm = [(1, 1, 6, -1)])
 (r = 4, s = 3, rs = 6, N = -3, comm = [(1, 1, 6, -3)])
```
""" UnipotentGroup

struct UnipotentElement{T}
  U::UnipotentGroup
  list::Vector{Pair{Int,T}}
end

function Base.show(io::IO,u::UnipotentElement)
  if hasdecor(io)
    if isempty(u.list) print(io,"()")
    else for (r,c) in u.list
      if get(io,:root,false)
        printTeX(io,"u_{"*joindigits(u.U.W.rootdec[r])*"}")
      else print(io,"u",r)
      end
      print(io,"(",c,")")
    end
    end
  else print(io,"U(",join(repr.(u.list),", "),")")
  end
end

Base.:*(u::UnipotentElement,v::UnipotentElement)=
  UnipotentElement(u.U,reorder(u.U,vcat(u.list, v.list)))

Base.:/(u::UnipotentElement,v::UnipotentElement)=u*inv(v)

Base.inv(u::UnipotentElement)=u.U(reverse(map(((r,c),)->r=>-c,u.list))...)

Base.:^(u::UnipotentElement,v::UnipotentElement)=inv(v)*u*v

Base.:(==)(u::UnipotentElement,v::UnipotentElement)=u.U==v.U && u.list==v.list

Base.copy(u::UnipotentElement)=UnipotentElement(u.U,copy(u.list))

Base.hash(u::UnipotentElement,h::UInt)=hash(u.list,h)

Base.:^(u::UnipotentElement,v::SemisimpleElement)=
  UnipotentElement(u.U,map(((r,c),)->r=>v^roots(u.U.W,r)*c,u.list))

n⁰(W,r)=findfirst(==(r),roots(W))

# Computes the constants ηᵣₛ defined in [Carter1972, 6.4.2 and 6.4.3]
# used for conjugating a unipotent element by a reflection.
η=function(U::UnipotentGroup,a,b)
  W=U.W
  L=map(j->n⁰(W,roots(W,a)*j+roots(W,b)),-4:4)
  L=filter(j->!isnothing(L[5+j]) && L[5+j]<=nref(W),-4:4)
  p=-L[1]
  q=L[end]
  eta=(-1)^p
  for i in 0:max(p-1, q-1)
    n=N(U,a,n⁰(W,roots(W,a)*(i-p)+roots(W,b)))
    if i<p || i<q eta*=sign(n) end
  end
  eta
end

function Base.:^(u::UnipotentElement,n::Perm)
  if isone(n) return u end
  W=u.U.W
  p=findfirst(x->isleftdescent(W,n,x[1]),u.list)
  if !isnothing(p) 
    error(u," should have no coefficient on root ", u.list[p][1][1], "\n")
  end
  s=firstleftdescent(W, n)
  u.U(map(((r,c),)->Int(r^n)=>η(u.U,s,r)*c,u.list)...)^(W(s)*n)
end

Base.one(u::UnipotentElement{T}) where T =UnipotentElement(u.U, Pair{Int,T}[])

Base.:^(u::UnipotentElement, n::Integer)=n>=0 ? Base.power_by_squaring(u,n) :
                                            Base.power_by_squaring(inv(u),-n)

# Let  < be  the ordering  of the  roots of  W according to the linear form
# determined by the weights cumsum((1//2).^(1:rank(W))) on the simple roots
# It  is checked that this induces a total  order on the roots and that the
# resulting order is the default order on the roots.
function check_root_order(W)
  l=roots(W,1:nref(W))*cumsum((1//2).^(1:rank(W)))
  if length(unique(l))!=length(l)
    error("the linear form takes the same value on two roots")
  end
  i=sortPerm(l)
  if !isone(i) error("Need ", i, " to sort the roots of ", W, "\n") end
end

# Compute N_{a,b} for non-necessarily positive roots [Carter1972b, 4.2.1
# (ii)](biblio.htm#Car72b] 
# Here a negative root α is represented by -(index of -α)
N(U::UnipotentGroup,a::Integer,b::Integer;scaled=false)=
   N(U.W,U.special,U.special,a,b;scaled)

function N(W,special,Nc,a::Integer,b::Integer;scaled=false)
  ra= a<0 ? -roots(W,-a) : roots(W,a)
  rb= b<0 ? -roots(W,-b) : roots(W,b)
  c=findfirst(==(ra+rb),roots(W))
  if isnothing(c) return 0 end
  l(r)=rootlengths(W,r)
  lc=l(c)
  if c>nref(W) c=-(c-nref(W)) end
  if a>0
    if b>0 
      p=findfirst(x->x.r==a && x.s==b,special)
      res=Nc[p].N
    elseif c<0 res=N(W,special,Nc,-c,a)*lc//l(-b)
    else res=-N(W,special,Nc,-b, c)*lc//l(a)
    end
  elseif b<0 res=-N(W,special,Nc,-a,-b)
  elseif c<0 res=N(W,special,Nc,b,-c)*lc//l(-a)
  else res=-N(W,special,Nc,c,-a)*lc//l(b)
  end
  scaled ? res//lc : res
end

function Ncarter(W,special)
# the order on the positive roots must be compatible with an order
# on the ambient vector space. 1:nref(W) is such an order (see
# check_root_order)
# Compute Nᵣₛ for each special pair, using arbitraryness for extraspecial
# pairs. See formula in proof of [Carter1972b, 4.2.2](biblio.htm#Car72b)
  r=s=rs=0
  ns=div(length(special),2)
  Nc=fill((N=0,),2*ns)
  for i in 1:ns
    if i==1 || special[i-1].rs!=special[i].rs #extraspecial
      (r,s,rs)=special[i]
      Nc[i]=(N=count(j->roots(W,s)-roots(W,r)*j in roots(W),0:3),)
    else # special of sum rs
      (r1,s1,rs1)=special[i]
      l=rootlengths(W)[rs1] # length of root r
      Nc[i]=(N=l//N(W,special,Nc,r,s)*(N(W,special,Nc,s,-r1;scaled=true)*
N(W,special,Nc,r,-s1)+N(W,special,Nc,-r1,r;scaled=true)*N(W,special,Nc,s,-s1)),)
    end
    Nc[ns+i]=(N=-Nc[i].N,) #fill in now; used for later i
  end
  Nc
end

"""
`UnipotentGroup(W;chevalley=nothing,order=1:nref(W))`

`W` should be a Weyl group. This function returns a `struct UnipotentGroup`
associated  to the reductive group of Weyl group `W`. 

If  the keyword `order` is given it is  a total order on the positive roots
used to normalize unipotent elements.

By  default the structure constants `Nᵣₛ`  are computed using the method of
[Carter1972](biblio.htm#Car72b)  from  extraspecial  pairs.  Another set of
structure  constants can be specified by giving for the keyword `chevalley`
a `Dict` associating to a pair `(r,s)` of root indices some object `p` such
that `first(p)` is the corresponding `N_{r,s}`. Here is an example of using
different constants from the package `ChevLie` of Meinolf Geck.
```julia-rep1
julia> W=coxgroup(:G,2)
G₂

julia> using ChevLie: LieAlg, structconsts

julia> l=LieAlg(:g,2)
#I dim = 14
LieAlg('G2')

julia> U=UnipotentGroup(W;chevalley=structconsts(l))
#I calculating structconsts
#I calculating eps-canonical base (100/.) 
UnipotentGroup(G₂)
```
"""
function UnipotentGroup(W::FiniteCoxeterGroup;chevalley=nothing,order=1:nref(W))
  # basic sanity check
  if roots(W,nref(W).+(1:nref(W)))!=-roots(W,1:nref(W)) error() end
  # compute special pairs for the order 1:nref(W). That is pairs
  # where r<s and r+s is a root
  special=NamedTuple{(:r,:s,:rs),NTuple{3,Int}}[]
  for s in 1:nref(W), r in 1:s-1
    pos=n⁰(W,roots(W,r)+roots(W,s))
    if !isnothing(pos) push!(special,(r=r,s=s,rs=pos)) end
  end
  sort!(special,by=x->(x.rs,x.r))
  ns=length(special) # half of final length
  append!(special,map(x->(r=x.s,s=x.r,rs=x.rs),special))
  # `special`:   triples of indices of roots `(r,s,r+s)`
  # where `(r,s)` is special and r<s, ordered by `(r+s,r)`, followed
  # by the triples `(s,r,r+s)` for the same list.
  if isnothing(chevalley) N=Ncarter(W,special) 
  else  N=fill((N=0,),2*ns)
    for i in eachindex(special)
      r,s,rs=special[i]
      N[i]=(N=first(chevalley[(r,s)]),)
    end
  end
  function M(r,s,i) # M_{r,s,i} of [Carter1972, bottom of page 61]
    m=1
    for j in 1:i
      d=n⁰(W,roots(W,r)+roots(W,s))
      if isnothing(d) || d>nref(W) return 0 end
      m*=N[findfirst(x->x.r==r && x.s==s,special)].N//j
      s=d
    end
    m
  end
  comm=[NTuple{4,Int}[] for i in 1:2*ns]
  for (i,(r,s,rs)) in enumerate(special)
    for c in [[1, 1], [2, 1], [1, 2], [3, 1], [1, 3], [3, 2], [2, 3]]
# possible (i,j) such that there may exist a root is+jr.
# see [Carter1972, top of page 76] for the formulas
      if c[2]==1 C=M(r,s,c[1])
      elseif c[1]==1 C=(-1)^c[2]*M(s,r,c[2])
      elseif c[[1,2]]==[3, 2] C=M(rs,r,2)//3
      elseif c[[1,2]]==[2, 3] C=(-2*M(rs,s,2))//3
      else C=0
      end
      if C!=0 
       push!(comm[i],(c[1],c[2],n⁰(W,c[1]*roots(W,r)+c[2]*roots(W,s)),Int(C)))
      end
    end
  end
  UnipotentGroup(W,1:nref(W),
    map((s,n,c)->(r=s.r,s=s.s,rs=s.rs,N=n.N,comm=c),special,N,comm),
        Dict{Symbol,Any}())
end

"""
'reorder(U,l[,order])'

This  function  takes  a  list  of  pairs  'r=>c'  representing a unipotent
element,  where 'r'  is a  root and  'c' the corresponding coefficient, and
puts  it in  canonical form,  reordering the  terms to agree with 'U.order'
using  the commutation  relations. If  a second  argument is given, this is
used instead of 'U.order'.

```julia-repl
julia> U=UnipotentGroup(coxgroup(:G,2))
UnipotentGroup(G₂)

julia> l=reorder(U,[2=>4,1=>2])
6-element Vector{Pair{Int64, Int64}}:
 1 => 2
 2 => 4
 3 => -8
 4 => 32
 5 => -128
 6 => 512

julia> reorder(U,l,6:-1:1)
2-element Vector{Pair{Int64, Int64}}:
 2 => 4
 1 => 2
```
"""
function reorder(U::UnipotentGroup,l,order=U.order)
  i=1
  l=copy(l)
  while i<=length(l)
    if iszero(l[i][2])
      splice!(l,i)
      if i>1 i-=1 end
    elseif i<length(l) && l[i][1]==l[i+1][1]
      splice!(l,i:i+1,[l[i][1]=>l[i][2]+l[i+1][2]])
    elseif i<length(l) &&
      findfirst(==(l[i][1]),order)>findfirst(==(l[i+1][1]),order)
# The Chevalley relation is
#
# uₛ(u)uᵣ(t)=uᵣ(t)uₛ(u)∏_{i,j>0} uᵢᵣ₊ⱼₛ(Cᵣₛᵢⱼ(-t)ⁱuʲ)
#
# Here l[i]=[s,u] and l[i+1]=[r,t].
      W=U.W
      res=l[[i+1,i]]
      s,u=l[i]; r,t=l[i+1]
      c=findfirst(p->p.r==r && p.s==s,U.special)
      if !isnothing(c)
        c=map(k->k[3]=>k[4]*(-t)^k[1]*u^k[2],U.special[c].comm)
        append!(res,filter(x->!iszero(x[2]),c))
      end
      splice!(l,i:i+1,res)
      if i>1 i-=1 end
    else i+=1
    end
  end
  return l
end

"""
'U(r)'

'U(r_1=>c_1`,..,r_n=>c_n`)'

Where `U` is a `UnipotentGroup`. In the first form the function creates the
element  `uᵣ(1)`,  and  in  the  second  form  the  element  `u_{r_1}(c_1)…
u_{r_n}(c_n)`

```julia-repl
julia> U=UnipotentGroup(coxgroup(:G,2))
UnipotentGroup(G₂)

julia> U(2)
u2(1)

julia> U(1=>2,2=>4)
u1(2)u2(4)

julia> U(2=>4,1=>2)
u1(2)u2(4)u3(-8)u4(32)u5(-128)u6(512)
```
"""
(U::UnipotentGroup)(i::Integer)=UnipotentElement(U, [i=>1])

(U::UnipotentGroup)(l::Pair{Int}...)=UnipotentElement(U,reorder(U,collect(l)))

Base.show(io::IO,U::UnipotentGroup)=print(io,"UnipotentGroup(",U.W,")")

"""
'abelianpart(u::UnipotentElement)'

If  `𝐔` is the unipotent subgroup and `D(𝐔)` its derived subgroup, this
function   returns  the  projection   of  the  unipotent   element  'u'  on
`𝐔/D(𝐔)`, that is its coefficients on the simple roots.

```julia-repl
julia> U=UnipotentGroup(coxgroup(:G,2));@Mvp x,y

julia> u=U(2=>y,1=>x)
u1(x)u2(y)u3(-xy)u4(xy²)u5(-xy³)u6(2x²y³)

julia> abelianpart(u)
u1(x)u2(y)
```
"""
abelianpart(u::UnipotentElement)=UnipotentElement(u.U,filter(x->x[1] in
                            1:ngens(u.U.W),u.list))

"""
`decompose(w,u::UnipotentElement)`

`w`  should be an element of the Weyl group corresponding to `u`. If `𝐔` is
the  unipotent radical  of the  Borel subgroup  determined by  the positive
roots,  and `𝐔⁻` the  opposite unipotent radical,  this function decomposes
`u` into its component in `𝐔 ∩ ʷ𝐔⁻` and its component in `𝐔 ∩ ʷ𝐔`.
```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> U=UnipotentGroup(W);@Mvp x,y

julia> u=U(2=>y,1=>x)
u1(x)u2(y)u3(-xy)u4(xy²)u5(-xy³)u6(2x²y³)

julia> decompose(W(1),u)
2-element Vector{UnipotentElement{Mvp{Int64, Int64}}}:
 u1(x)
 u2(y)u3(-xy)u4(xy²)u5(-xy³)u6(2x²y³)

julia> decompose(W(2),u)
2-element Vector{UnipotentElement{Mvp{Int64, Int64}}}:
 u2(y)
 u1(x)
```
"""
function Chars.decompose(w,u::UnipotentElement)
  U=u.U
  W=U.W
  order=vcat(filter(i->isleftdescent(W,w,i),U.order),
             filter(i->!isleftdescent(W,w,i),U.order))
  l=reorder(U,u.list, order)
  [U(filter(i->isleftdescent(W,w,i[1]),l)...),
   U(filter(i->!isleftdescent(W,w,i[1]),l)...)]
end

end
