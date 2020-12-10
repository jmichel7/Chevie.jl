"""
This  module contains  functions for  computing with  unipotent elements of
reductive  groups;  specifically  to  compute  with  elements  of unipotent
radical  of a Borel subgroup of  a connected algebraic reductive group; the
implementation of these functions was initially written by Olivier Dudas in
GAP3.

The  unipotent radical of a  Borel subgroup is the  product in any order of
root  subgroups associated  to the  positive roots.  We fix an order, which
gives a canonical form to display elements and to compare them.

The  computations use the Steinberg relations between root subgroups, which
come from the choice of a Chevalley basis of the Lie algebra. The reference
we  follow is [Carter1972, chapters 4 to 6](biblio.htm#Car72b).

We  start with  a root  datum specified  by a  Weyl group  `W` and  build a
`struct` which contains information about the maximal unipotent subgroup of
the  corresponding reductive  group, that  is the  unipotent radical of the
Borel subgroup determined by the positive roots.

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
julia> xprint(U(4=>5,2=>4);root=true)
u₀₁₀₀₀₀(4)u₀₀₀₁₀₀(5)u₀₁₀₁₀₀(-20)
```

The coefficients in the root subgroups can be elements of arbitrary rings.
Here is an example using `Mvp`'s:

```julia-repl
julia> W=coxgroup(:E,8);U=UnipotentGroup(W)
UnipotentGroup(E₈)

julia> u=U(map(i->i=>Z(2)*Mvp(Symbol("x_",i)),1:8)...)
u1(1₂x₁)u2(1₂x₂)u3(1₂x₃)u4(1₂x₄)u5(1₂x₅)u6(1₂x₆)u7(1₂x₇)u8(1₂x₈)
```

```julia-rep1
julia> cut(sprint(show,u^16;context=rio(root=true)),before="u",width=60)
u₂₂₃₄₃₂₁₀(1₂x₁²x₂²x₃³x₄⁴x₅³x₆²x₇)
u₁₂₃₄₃₂₁₁(1₂x₁x₂²x₃³x₄⁴x₅³x₆²x₇x₈)
u₁₂₂₄₃₂₂₁(1₂x₁x₂²x₃²x₄⁴x₅³x₆²x₇²x₈)
u₁₂₂₃₃₃₂₁(1₂x₁x₂²x₃²x₄³x₅³x₆³x₇²x₈)
u₂₂₃₄₃₂₁₁(1₂x₁²x₂²x₃³x₄⁴x₅³x₆²x₇x₈)
u₁₂₂₄₃₃₂₁(1₂x₁x₂²x₃²x₄⁴x₅³x₆³x₇²x₈)
u₁₂₂₄₄₃₂₁(1₂x₁x₂²x₃²x₄⁴x₅⁴x₆³x₇²x₈)
u₂₂₃₄₃₃₂₁(1₂x₁²x₂²x₃³x₄⁴x₅³x₆³x₇²x₈)
u₁₂₃₄₄₃₂₁(1₂x₁x₂²x₃³x₄⁴x₅⁴x₆³x₇²x₈)
u₂₂₃₄₄₃₂₁(1₂x₁²x₂²x₃³x₄⁴x₅⁴x₆³x₇²x₈)
u₂₃₃₅₄₃₂₁(1₂x₁²x₂³x₃³x₄⁵x₅⁴x₆³x₇²x₈)
u₂₂₄₅₄₃₂₁(1₂x₁²x₂²x₃⁴x₄⁵x₅⁴x₆³x₇²x₈)
u₂₃₄₆₅₄₃₂(1₂x₁²x₂³x₃⁴x₄⁶x₅⁵x₆⁴x₇³x₈²)
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
using ..Gapjm
export UnipotentElement, UnipotentGroup, norm, abelianpart

"""
A  `struct UnipotentGroup`  `U` represents  the unipotent  radical `𝐔` of a
Borel subgroup of the reductive group with Weyl group `U.W`.

See   [Carter1972,  section  4.2](biblio.htm#Car72b)  for  details  on  the
following.  A Chevalley basis  of the Lie  algebra of `𝐔`  is a basis `eᵣ`,
where   each  `eᵣ`  is  in  the  corresponding  root  subspace,  such  that
`[eᵣ,eₛ]=Nᵣₛ e_{r+s}` for some integer constants `Nᵣₛ`.

To  build such  a basis,  let `<`  be a  total order  on the positive roots
induced  by a total order on the ambient vector space (the default order of
roots  of `W` in this package  is an example; we use  it in this module). A
pair `(r,s)` of roots is *special* if `0<r<s` and `r+s` is a root.

Constants  `Cᵣₛᵢⱼ` are defined, see [Carter1972, 5.2.3](biblio.htm#Car72b),
by

``u_s(u)u_r(t)=u_r(t)u_s(u)\\prod_{i,j>0}u_{ir+js}(C_{rsij}(-t)^iu^j)``

Where  `ir+js` runs over the positive  integral combinations of `r` and `s`
which are roots, taken in lexicographic order.

The fields of `struct Unipotent Group` are:

  * `W`:         the underlying Weyl group
  * `specialPairs`:   triples of indices of the roots `(r,s,r+s)`
                  where `(r,s)` is special, ordered by `(r+s,r)`, followed
                  by the triples `(s,r,r+s)` for the same list.
  * `ns`:         The number of   `specialPairs` where `r<s`.
  * `N`:          the constants `Nᵣₛ` for `specialPairs`
  * `order`:      the order on positive roots used to normalize products
  * `commutatorConstants`: stores the `Cᵣₛᵢⱼ` by storing for each special pair
                     `(r,s)` the list of quadruples `[i,j,ir+js,Cᵣₛᵢⱼ]`.
"""
struct UnipotentGroup
  W
  specialPairs::Vector{Vector{Int}}
  ns::Int
  N::Vector{Int}
  order
  commutatorConstants
end

struct UnipotentElement{T}
  U::UnipotentGroup
  list::Vector{Pair{Int,T}}
end

function Base.show(io::IO,u::UnipotentElement)
  if get(io,:limit,false) || get(io,:TeX,false)
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

Base.:*(u::UnipotentElement,v::UnipotentElement)=u.U(vcat(u.list, v.list)...)

Base.:/(u::UnipotentElement,v::UnipotentElement)=u*inv(v)

Base.inv(u::UnipotentElement)=u.U(reverse(map(p->p[1]=>-p[2],u.list))...)

Base.:^(u::UnipotentElement,v::UnipotentElement)=inv(v)*u*v

Base.:(==)(u::UnipotentElement,v::UnipotentElement)=u.U==v.U && u.list==v.list

function Base.:^(u::UnipotentElement,v::SemisimpleElement)
  UnipotentElement(u.U, map(p->p[1]=>v^roots(u.U.W,p[1])*p[2],u.list))
end

# Computes the constants ηᵣₛ defined in [Carter1972, 6.4.2 and 6.4.3]
# used for conjugating a unipotent element by a reflection.
η=function(U::UnipotentGroup,a,b)
  W=U.W
  no(r)=findfirst(==(r),roots(W))
  L=map(j->no(roots(W,a)*j+roots(W,b)),-4:4)
  L=filter(j->!isnothing(L[5+j]) && L[5+j]<=nref(W),-4:4)
  p=-L[1]
  q=L[end]
  eta=(-1)^p
  for i in 0:max(p-1, q-1)
    n=N(U,a,no(roots(W,a)*(i-p)+roots(W,b)))
    if i<p || i<q eta*=sign(n) end
  end
  eta
end

function Base.:^(u::UnipotentElement,n::Perm)
  if isone(n) return u end
  W=u.U.W
  s=firstleftdescent(W, n)
  p=filter(x->isleftdescent(W,n,x[1]),u.list)
  if !isempty(p) error(u," should have no coefficient on root ", p[1][1], "\n") end
  u.U(map(i->Int(i[1]^n)=>η(u.U,s,i[1])*i[2], u.list)...)^(W(s)*n)
end

Base.one(u::UnipotentElement)=UnipotentElement(u.U, Vector{Int}[])

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

# Compute N_{a,b} for non-necessarily positive roots [@Car72b, 4.2.1 (ii)]
function N(U::UnipotentGroup,a::Integer,b::Integer;scaled=false)
  W=U.W
  ra= a<0 ? -roots(W,-a) : roots(W,a)
  rb= b<0 ? -roots(W,-b) : roots(W,b)
  c=findfirst(==(ra+rb),roots(W))
  if isnothing(c) return 0 end
  l(r)=isnothing(r) ? 1 : rootlengths(W)[r] # length of root r
  lc=l(c)
  if c>nref(W) c=-(c-nref(W)) end
  if a>0
    if b>0 res=U.N[findfirst(==([a,b,c]),U.specialPairs)]
    elseif c<0 res=N(U,-c,a)*lc//l(-b)
    else res=-N(U,-b, c)*lc//l(a)
    end
  elseif b<0 res=-N(U,-a,-b)
  elseif c<0 res=N(U,b,-c)*lc//l(-a)
  else res=-N(U,c,-a)*lc//l(b)
  end
  scaled ? res//lc : res
end

"""
`UnipotentGroup(W)`

`W`  should be a Weyl group.  This function returns a `struct` representing
the  unipotent radical `𝐔`  of a Borel  subgroup of the  reductive group of
Weyl group `W`.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> U=UnipotentGroup(W)
UnipotentGroup(G₂)

julia> U.specialPairs
10-element Array{Array{Int64,1},1}:
 [1, 2, 3]
 [2, 3, 4]
 [2, 4, 5]
 [1, 5, 6]
 [3, 4, 6]
 [2, 1, 3]
 [3, 2, 4]
 [4, 2, 5]
 [5, 1, 6]
 [4, 3, 6]

julia> U.N
10-element Array{Int64,1}:
  1
  2
  3
  1
  3
 -1
 -2
 -3
 -1
 -3

julia> U.commutatorConstants
10-element Array{Array{Array{Int64,1},1},1}:
 [[1, 1, 3, 1], [1, 2, 4, -1], [1, 3, 5, 1], [2, 3, 6, 2]]
 [[1, 1, 4, 2], [2, 1, 5, 3], [1, 2, 6, -3]]
 [[1, 1, 5, 3]]
 [[1, 1, 6, 1]]
 [[1, 1, 6, 3]]
 [[1, 1, 3, -1], [2, 1, 4, -1], [3, 1, 5, -1], [3, 2, 6, -1]]
 [[1, 1, 4, -2], [2, 1, 6, -3], [1, 2, 5, 3]]
 [[1, 1, 5, -3]]
 [[1, 1, 6, -1]]
 [[1, 1, 6, -3]]
```
"""
function UnipotentGroup(W::FiniteCoxeterGroup)
  if roots(W,nref(W)+(1:nref(W)))!=-roots(W,1:nref(W)) error() end
  no(r)=findfirst(==(r),roots(W))
  # compute special pairs. We take `1:nref(W)` as order on the roots.
  special=Vector{Int}[]
  for s in 1:nref(W), r in 1:s-1
    pos=no(roots(W,r)+roots(W,s))
    if !isnothing(pos) push!(special,[r,s,pos]) end
  end
  sort!(special,by=x->x[[3,1]])
  ns=length(special)
  append!(special,map(x->x[[2,1,3]],special))
  # we initialize U with the order on U_α given by 1:nref(W)
  U=UnipotentGroup(W,special,ns,fill(0,2*ns),1:nref(W),fill(Vector{Int}[],2*ns))
  l(r)=isnothing(r) ? 1 : rootlengths(W)[r] # length of root r
# Compute Nᵣₛ for each special pair. See formula in proof of [@Car72b, 4.2.2]
  r=s=rs=0
  for i in 1:U.ns
   if i==1 || U.specialPairs[i-1][3]!=U.specialPairs[i][3] #extraspecial
      (r,s,rs)=U.specialPairs[i]
      U.N[i]=count(j->roots(W,s)-roots(W,r)*j in roots(W),0:3)
    else # special of sum rs
      (r1,s1,rs1)=U.specialPairs[i]
      U.N[i]=l(rs1)//N(U,r,s)*(N(U,s,-r1;scaled=true)*N(U,r,-s1)+
                  N(U,-r1,r;scaled=true)*N(U,s,-s1))
    end
    U.N[U.ns+i]=-U.N[i] #fill in now; used for next i
  end
  function M(r,s,i) # M_{r,s,i} of [Carter1972, bottom of page 61]
    m=1
    for j in 1:i
      d=no(roots(W,r)+roots(W,s))
      if isnothing(d) || d>nref(W) return 0 end
      m*=U.N[findfirst(==([r,s,d]),U.specialPairs)]//j
      s=d
    end
    m
  end
  for i in eachindex(U.specialPairs)
    (r,s,rs)=U.specialPairs[i]
    L=Vector{Int}[]
    for c in [[1, 1], [2, 1], [1, 2], [3, 1], [1, 3], [3, 2], [2, 3]]
# possible (i,j) such that there may exist a root is+jr.
# see [Carter1972, top of page 76] for the formulas
      if c[2]==1 C=M(r,s,c[1])
      elseif c[1]==1 C=(-1)^c[2]*M(s,r,c[2])
      elseif c[[1,2]]==[3, 2] C=M(rs,r,2)//3
      elseif c[[1, 2]]==[2, 3] C=(-2*M(rs,s,2))//3
      else C=0
      end
      if C!=0 push!(L,[c[1],c[2],no(c[1]*roots(W,r)+c[2]*roots(W,s)),C]) end
    end
    U.commutatorConstants[i]=L
  end
  U
end

"""
'norm(U,l[,order])'

This  function  takes  a  list  of  pairs  'r=>c'  representing a unipotent
element,  where 'r'  is a  root and  'c' the corresponding coefficient, and
puts  it in  canonical form,  reordering the  terms to agree with 'U.order'
using  the commutation  relations. If  a second  argument is given, this is
used instead of 'U.order'.

```julia-repl
julia> U=UnipotentGroup(coxgroup(:G,2))
UnipotentGroup(G₂)

julia> l=norm(U,[2=>4,1=>2])
6-element Array{Pair{Int64,Int64},1}:
 1 => 2
 2 => 4
 3 => -8
 4 => 32
 5 => -128
 6 => 512

julia> norm(U,l,6:-1:1)
2-element Array{Pair{Int64,Int64},1}:
 2 => 4
 1 => 2
```
"""
function norm(U::UnipotentGroup,l,order=U.order)
  W=U.W
  i=1
  while i<=length(l)
    if iszero(l[i][2])
      l=vcat(l[1:i-1],l[i+1:end])
      if i>1 i-=1 end
    elseif i<length(l) && l[i][1]==l[i+1][1]
      splice!(l,i:i+1,[l[i][1]=>l[i][2]+l[i+1][2]])
    elseif i<length(l) &&
      findfirst(==(l[i][1]),order)>findfirst(==(l[i+1][1]),order)
# The Chevalley relation is
#
# uₛ(u)uᵣ(t)=uᵣ(t)uₛ(u)∏_{i,j>0} u_{ir+js}(Cᵣₛᵢⱼ(-t)^iu^j)
#
# Here l[i]=[s,u] and l[i+1]=[r,t].
      res=l[vcat(1:i-1,[i+1,i])]
      s,u=l[i]; r,t=l[i+1]
      c=findfirst(==(roots(W,r)+roots(W,s)),roots(W))
      if !isnothing(c) c=findfirst(==([r,s,c]),U.specialPairs) end
      if !isnothing(c)
        c=map(k->k[3]=>k[4]*(-t)^k[1]*u^k[2],U.commutatorConstants[c])
        append!(res,filter(x->!iszero(x[2]),c))
      end
      append!(res,l[i+2:end])
      l=res
      if i>1 i-=1 end
    else i+=1
    end
  end
  return l
end

"""
'U(r)'

'U(r_1=>c_1`,..,r_n=>c_n`)'

In the first form the function creates the element `uᵣ(1)`, and in the second
form the element `u_{r_1}(c_1)… u_{r_n}(c_n)`

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

(U::UnipotentGroup)(l::Pair{Int,T}...) where T=UnipotentElement(U,
                                                norm(U,collect(l)))

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
`decompose(w,u)`

`u`  should be a unipotent element and  `w` an element of the corresponding
Weyl  group.  If  `𝐔`  is  the  unipotent  radical  of  the  Borel subgroup
determined  by the  positive roots,  and `𝐔⁻`  the unipotent radical of the
opposite  Borel, this  function decomposes  `u` into  its component in `𝐔 ∩
ʷ𝐔⁻` and its component in `𝐔 ∩ ʷ𝐔 `.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> U=UnipotentGroup(W);@Mvp x,y

julia> u=U(2=>y,1=>x)
u1(x)u2(y)u3(-xy)u4(xy²)u5(-xy³)u6(2x²y³)

julia> decompose(W(1),u)
2-element Array{UnipotentElement{Mvp{Int64,Int64}},1}:
 u1(x)
 u2(y)u3(-xy)u4(xy²)u5(-xy³)u6(2x²y³)

julia> decompose(W(2),u)
2-element Array{UnipotentElement{Mvp{Int64,Int64}},1}:
 u2(y)
 u1(x)
```
"""
function Chars.decompose(w,u::UnipotentElement)
  U=u.U
  W=U.W
  order=vcat(filter(i->isleftdescent(W,w,i),U.order),
             filter(i->!isleftdescent(W,w,i),U.order))
  l=norm(U,u.list, order)
  [U(filter(i->isleftdescent(W,w,i[1]),l)...),
   U(filter(i->!isleftdescent(W,w,i[1]),l)...)]
end

end
