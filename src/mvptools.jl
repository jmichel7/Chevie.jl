# this should just be "tools". It contains simple functions but which need
# several of the self-contained structural packages of Gapjm
"""
`Mvp(p)` converts  the `Pol`  `p` to  an  `Mvp`. 

```julia-repl
julia> Pol(:q)
Pol{Int64}: q

julia> Mvp(q^2+q)
Mvp{Int64}: q²+q
```
"""
Mvp(x::Pol)=convert(Mvp,x)
Base.convert(::Type{Mvp{T,N}},p::Pol) where{T,N}=
                     p(Mvp(convert(Monomial{N},Pols.varname[])=>one(T)))
Base.convert(::Type{Mvp},p::Pol)=p(Mvp(Pols.varname[]))

"""
`factor(p::Mvp)`

`p`  should be of degree <=2 thus represents a quadratic form. The function
returns  a list  of two  linear forms  of which  `p` is the product if such
exist, otherwise it returns [p].

```julia-repl
julia> factor(x^2-y^2+x+3y-2)
2-element Array{Mvp{Int64,Int64},1}:
 x+y-1
 x-y+2

julia> factor(x^2+x+1)
2-element Array{Mvp{Cyc{Int64},Int64},1}:
 x-ζ₃
 x-ζ₃²

julia> factor(x*y-1)
1-element Array{Mvp{Int64,Int64},1}:
 xy-1
```
"""
function Util.factor(p::Mvp{T,N})where {T,N}
  v=variables(p)
  r=length(v)+1
  m=zeros(T,r,r)//1
  for (e,t) in p.d
    n=map(x->findfirst(==(x),v),keys(e.d))
    c=values(e.d)
    if c==[1,1] m[n[1],n[2]]=m[n[2],n[1]]=t//2
    elseif c==[2] m[n[1],n[1]]=t
    elseif c==[1] m[n[1],r]=m[r,n[1]]=t//2
    elseif isempty(c) m[r,r]=t
    else error("# only implemented for degree <=2")
    end
  end
  if size(m,1)==2 t=one(m)
  else n=copy(m)
    m,ind=echelon!(m)
    m=m[ind,:]
    if size(m,1)>2 return [p] end
    t=permutedims(solutionmat(m,n))
    m=solutionmat(t,m)
  end
  v=t*vcat(Mvp.(v),[Mvp(1)])
  if size(m,1)==1 return [v[1],v[1]*m[1,1]] end
  b=m[1,2]+m[2,1]
  if m[1,1]==0 return [v[2],b*v[1]+m[2,2]*v[2]] end
  b//=m[1,1]
  d=root(b^2-4*m[2,2]//m[1,1])
  if isnothing(d) 
    println("root failed")
    return p 
  end
  improve_type([v[1]+v[2]//2*(b-d),m[1,1]*(v[1]+v[2]//2*(b+d))])
end

"""
`mod(z::Cyc,p::Integer)`

`p`  should be a  prime and `z`  a cyclotomic number  which is `p`-integral
(that  is, `z` times some number prime to `p` is a cyclotomic integer). The
function  returns  the  reduction  of  `z`  mod.  `p`,  an  element of some
extension  `F_(p^r)` of the prime field  `F_p`.

```julia_repl
julia> mod(E(7),3)
Z₇₂₉¹⁰⁴
```
"""
function Base.mod(c::Cyc,p)
  x=coefficients(c) 
  n=conductor(c)
  np=MatInt.prime_part(n,p)
  pp=div(n,np)
  r=order(Mod{np}(p)) # order of p mod np
  zeta=Z(p^r)^(div(p^r-1,np)*gcdx(pp,p^r-1)[2]) #n-th root of unity
  if !isone(zeta^n) error() end
  sum(i->zeta^(i-1)*x[i],1:n)
end;

export pblocks
"""
`pblocks(G,p)`

Let  `p` be a prime. This function returns the partition of the irreducible
characters  of `G`  in `p`-blocks,  represented by  the list  of indices of
irreducibles characters in each block.

gap> W:=SymmetricGroup(5);
Group( (1,5), (2,5), (3,5), (4,5) )
gap> PBlocks(W,2);
[ [ 1, 2, 5, 6, 7 ], [ 3, 4 ] ]
gap> PBlocks(W,3);
[ [ 1, 3, 6 ], [ 2, 4, 5 ], [ 7 ] ]
gap> PBlocks(W,7);
[ [ 1 ], [ 2 ], [ 3 ], [ 4 ], [ 5 ], [ 6 ], [ 7 ] ]
"""
function pblocks(G,p)
  T=CharTable(G)
  l=length(T.charnames)
  classes=map(c->div(T.centralizers[1],c),T.centralizers)
  v=map(chi->map(j->mod(classes[j]*chi[j]//chi[1],p),1:l),eachrow(T.irr))
  sort(collect(values(groupby(improve_type(v),1:l))))
end
