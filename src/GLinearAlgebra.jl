"""
GLinearAlgebra: linear algebra over arbitrary fields and rings

The  linear  algebra  package  in  Julia  is  not  suitable  for  a general
mathematics  package: it assumes  the field is  the Real or Complex numbers
and uses floating point to do approximate computations.
Here we are interested in functions which work over any field (or sometimes
any ring).
"""
module GLinearAlgebra
using Gapjm
export echelon, exterior_power, CoFactors, bigcell_decomposition, diagblocks, 
  ratio

"""
    `echelon!(m)`

  puts `m` in echelon form and returns:
  `m`, indices of linearly independent rows of `m`
  works in any field.
  The echelon form transforms the rows of m into a particular basis
  of the rowspace. The first non-zero element of each line is 1, and
  such an element is also the only non-zero in its column.
"""
function echelon!(m::Matrix)
  T=typeof(one(eltype(m))//1)
  if T!=eltype(m) m=convert.(T,m) end
  rk=0
  inds=collect(axes(m,1))
  for k in axes(m,2)
    j=findfirst(!iszero,m[rk+1:end,k])
    if isnothing(j) continue end
    j+=rk
    rk+=1
    row=m[j,:]
    m[j,:].=m[rk,:]
    m[rk,:].=inv(row[k]).*row
    inds[[j,rk]]=inds[[rk,j]]
    for j in axes(m,1)
      if rk!=j && !iszero(m[j,k]) m[j,:].-=m[j,k].*m[rk,:] end
    end
#   println(m)
  end
  m,inds[1:rk]
end

echelon(m::Matrix)=echelon!(copy(m))

"""
 computes rank of m
 not exported to avoid conflict with LinearAlgebra
"""
rank(m::Matrix)=length(echelon(m)[2])

"""
 computes right nullspace of m in a type_preserving way
 not exported to avoid conflict with LinearAlgebra
"""
function nullspace(m::Matrix)
  m=echelon(m)[1]
  n=size(m)[2]
  z=Int[]
  j=0
  lim=size(m,1)
  for i in axes(m,1)
    f=findfirst(!iszero,m[i,j+1:end])
    if isnothing(f)
      lim=i-1
      break
    end
    j+=f
    push!(z,j)
  end
# println("z=$z lim=$lim")
  zz=zeros(eltype(m),n,n)
  zz[z,:]=m[1:lim,:]
  nn=filter(k->iszero(zz[k,k]),1:n)
  for i in nn zz[i,i]=-one(eltype(m)) end
  zz[:,nn]
end

function det(A::Matrix)
  if size(A,1)==1 return A[1,1] end
  sum(i->A[i,1]*det(A[vcat(1:i-1,i+1:size(A,1)),2:end])*(-1)^(i-1),axes(A,1))
end

Pols.isunit(c::Cyc)=true
Pols.isunit(c::Any)=false

function Det(m)
  if isempty(m) return 0 end
  n=size(m,1)
  if n==1 return m[1,1]
  elseif n==2 return m[1,1]*m[2,2]-m[1,2]*m[2,1]
  elseif n==3
    return sum(p->prod(i->m[i,i^p],1:n)*sign(p),elements(symmetric_group(n)))
  end
  function compl(m,i,j)
    v=axes(m,1)
    m[setdiff(v,[i]),setdiff(v,[j])]
  end
  i=findfirst(i->count(!iszero,m[i,:])<=2,axes(m,1))
  if !isnothing(i)
    j=findall(!iszero,m[i,:])
    if isempty(j) return 0 end
    return sum(k->(-1)^(i+k)*m[i,k]*Det(compl(m,i,k)),j)
  end
  v=axes(m,1)
# if length(v)<=3 return det*Det(m) end
  for j in v
    i=findfirst(isunit,m[v,j])
    if isnothing(i) continue end
#   println(size(m,1),":",[j,i])
    n=copy(m)
    f=inv(m[i,j])
    for k in setdiff(v,[i]) n[k,:]-=m[k,j]*f.*m[i,:] end
    return (-1)^(i+j)*n[i,j]*Det(compl(n,i,j))
  end
  print("m=");display(iszero.(m))
  v=map(x->count(!iszero,x),toL(m))
  j=findfirst(isequal(minimum(v)),v)
  if length(m)>71 println("\n",length(m),":",v[j]) end
# if v[j]>5 return det*DeterminantMat(m) end
  sum(function(i)
#   if Length(m) mod 5=0 then Print(Length(m),":",i," \c");fi;
    if iszero(m[i,j]) return 0
    else return (-1)^(i+1)*m[i,j]*Det(compl(m,i,j))
   end end,axes(m,1))
end

# Cofactors of the square matrix M, defined by CoFactors(M)*M=Det(M)*one(M)
function CoFactors(m)
  if size(m,1)==1 return fill(1,1,1) end
  v=axes(m,1)
  permutedims([(-1)^(i+j)*Det(m[filter(x->x!=i,v),filter(x->x!=j,v)])
               for i in v, j in v])
end

"""
`bigcell_decomposition(M [, b])`

`M`  should be a square  matrix, and `b` specifies  a block structure for a
matrix  of  same  size  as  `M`  (it  is  a  `Vector`  of  `Vector`s  whose
concatenation  is `1:size(M,1)`).  If `b`  is not  given, the trivial block
structure `[[i] for i in axes(M,1)]` is assumed.

The  function  decomposes  `M`  as  a  product  `P₁ L P` where `P` is upper
block-unitriangular   (with  identity  diagonal   blocks),  `P₁`  is  lower
block-unitriangular  and `L` is block-diagonal for the block structure `b`.
If  `M` is symmetric then  `P₁` is the transposed  of `P` and the result is
the  pair  `[P,L]`;  else  the  result  is  the triple `[P₁,L,P]`. The only
condition  for  this  decomposition  of  `M`  to  be  possible  is that the
principal  minors  according  to  the  block  structure be invertible. This
routine  is used  in the  Lusztig-Shoji algorithm  for computing  the Green
functions  and the example  below is extracted  from the computation of the
Green functions for `G₂`.

```julia-repl
julia> Pol(:q)
Pol{Int64}: q

julia> M=[q^6 q^0 q^3 q^3 q^5+q q^4+q^2; q^0 q^6 q^3 q^3 q^5+q q^4+q^2; q^3 q^3 q^6 q^0 q^4+q^2 q^5+q; q^3 q^3 q^0 q^6 q^4+q^2 q^5+q; q^5+q q^5+q q^4+q^2 q^4+q^2 q^6+q^4+q^2+1 q^5+2*q^3+q; q^4+q^2 q^4+q^2 q^5+q q^5+q q^5+2*q^3+q q^6+q^4+q^2+1]
6×6 Array{Pol{Int64},2}:
 q⁶     1      q³     q³     q⁵+q        q⁴+q²
 1      q⁶     q³     q³     q⁵+q        q⁴+q²
 q³     q³     q⁶     1      q⁴+q²       q⁵+q
 q³     q³     1      q⁶     q⁴+q²       q⁵+q
 q⁵+q   q⁵+q   q⁴+q²  q⁴+q²  q⁶+q⁴+q²+1  q⁵+2q³+q
 q⁴+q²  q⁴+q²  q⁵+q   q⁵+q   q⁵+2q³+q    q⁶+q⁴+q²+1

julia> bb=[[2],[4],[6],[3,5],[1]];

julia> (P,L)=bigcell_decomposition(M,bb);

julia> P
6×6 Array{Pol{Int64},2}:
 1    0  0    0    0        0
 q⁻⁶  1  q⁻³  q⁻³  q⁻¹+q⁻⁵  q⁻²+q⁻⁴
 0    0  1    0    0        0
 q⁻³  0  0    1    q⁻²      q⁻¹
 q⁻¹  0  0    0    1        0
 q⁻²  0  q⁻¹  0    q⁻¹      1

julia> L
6×6 Array{Pol{Int64},2}:
 q⁶-q⁴-1+q⁻²  0   0            0     0            0
 0            q⁶  0            0     0            0
 0            0   q⁶-q⁴-1+q⁻²  0     0            0
 0            0   0            q⁶-1  0            0
 0            0   0            0     q⁶-q⁴-1+q⁻²  0
 0            0   0            0     0            q⁶-1

julia> M==permutedims(P)*L*P
true
```
"""
function bigcell_decomposition(M,b=map(i->[i],axes(M,1)))
  L=one(M)
  P=one(M)
  block(X,i,j)=X[b[i],b[j]]
  if M==permutedims(M)
    for j in eachindex(b)
      L[b[j],b[j]]=block(M,j,j)
      if j>1 L[b[j],b[j]]-=sum(k->block(P,j,k)*block(L,k,k)*
                              permutedims(block(P,j,k)),1:j-1) end
      cb=CoFactors(block(L,j,j))
      db=Det(block(L,j,j))
      for i in j+1:length(b)
        P[b[i],b[j]]=block(M,i,j)
        if j>1 P[b[i],b[j]]-=
          sum(k->block(P,i,k)*block(L,k,k)*permutedims(block(P,j,k)),1:j-1)
        end
        P[b[i],b[j]]*=cb
        P[b[i],b[j]]=div.(block(P,i,j),db)
      end
    end
    return permutedims(P),L
  end
  tP=one(M)
  for j in eachindex(b)
    L[b[j],b[j]]=block(M,j,j)-sum(k->block(P,j,k)*block(L,k,k)*
                                  permutedims(block(P,j,k)),1:j-1)
    cb=CoFactors(block(L,j,j))
    db=Det(block(L,j,j))
    for i in j+1:length(b)
      P[b[i],b[j]]=block(M,i,j)-sum(k->block(P,i,k)*block(L,k,k)*
                                     block(tP,k,j),1:j-1)
      P[b[i],b[j]]=div.(P[b[i],b[j]]*cb,db)
      tP[b[j],b[i]]=cb*(block(M,j,i)-sum(k->block(P,j,k)*
                                         block(L,k,k)*block(tP,k,i),1:j-1))
      tP[b[i],b[j]]=div.(tP[b[i],b[j]],db)
    end
  end
  (P,L,tP)
end

"""
`exterior_power(mat,n)`

`mat`  should be a square matrix.  The function returns the `n`-th exterior
power  of  `mat`,  in  the  basis naturally indexed by`combinations(1:r,n)`
where`r=size(mat,1)`

```julia-repl
julia> M=[1 2 3 4;2 3 4 1;3 4 1 2;4 1 2 3]
4×4 Array{Int64,2}:
 1  2  3  4
 2  3  4  1
 3  4  1  2
 4  1  2  3

julia> exterior_power(M,2)
6×6 Array{Int64,2}:
  -1   -2   -7   -1  -10  -13
  -2   -8  -10  -10  -12    2
  -7  -10  -13    1    2    1
  -1  -10    1  -13    2    7
 -10  -12    2    2    8   10
 -13    2    1    7   10   -1
```
"""
function exterior_power(A,m)
  basis=combinations(1:size(A,1),m)
  [det(A[i,j]) for i in basis, j in basis]
end

"""
`diagblocks(M::Matrix)`

`M`  should  be  a  square  matrix.  Define  a  graph  `G`  with vertices
`1:size(M,1)` and with an edge between `i`  and `j` if either `M[i,j]` or
`M[j,i]` is not zero or `false`. `diagblocks` returns a vector of vectors
`I`  such that  `I[1]`,`I[2]`, etc..  are the  vertices in each connected
component  of `G`.  In other  words, `M[I[1],I[1]]`,`M[I[2],I[2]]`,etc...
are diagonal blocks of `M`.

```julia-repl
julia> m=[0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0]
4×4 Array{Int64,2}:
 0  0  0  1
 0  0  1  0
 0  1  0  0
 1  0  0  0

julia> diagblocks(m)
2-element Array{Array{Int64,1},1}:
 [1, 4]
 [2, 3]

julia> m[[1,4],[1,4]]
2×2 Array{Int64,2}:
 0  1
 1  0
```
"""
function diagblocks(M::Matrix)::Vector{Vector{Int}}
  l=size(M,1)
  if l==0 return Vector{Int}[] end
  cc=collect(1:l) # cc[i]: in which block is i, initialized to different blocks
  for i in 1:l, j in i+1:l
    # if new relation i~j then merge components:
    if !(iszero(M[i,j]) && iszero(M[j,i])) && cc[i]!=cc[j]
      cj=cc[j]
      for k in 1:l
         if cc[k]==cj cc[k]=cc[i] end
      end
    end
  end
  sort(collect(values(groupby(cc,collect(1:l)))))
end

"""
`Transporter(l1, l2 )`

`l1`  and `l2` should be vectors of  the same length of square matrices all
of the same size. The result is a basis of the vector space of matrices `A`
such  that for any `i` we have  `A*l1[i]=l2[i]*A` --- the basis is returned
as  a vector of matrices, empty if the vector space is 0. This is useful to
find whether two representations are isomorphic.
"""
function Transporter(l1::Vector{<:Matrix}, l2::Vector{<:Matrix})
  n=size(l1[1],1)
  M=Vector{eltype(l1[1])}[]
  for i in 1:n, j in 1:n,  m in 1:length(l1)
    v=zero(l1[1])
    v[i,:]+=l1[m][:,j]
    v[:,j]-=l2[m][i,:]
    push!(M, vec(v))
  end
  M=echelon(toM(M))[1]
  M=M[filter(i->!iszero(M[i,:]),axes(M,1)),:]
  v=nullspace(M)
  if isempty(v) return v end
  map(w->reshape(w,(n,n)), eachcol(v))
end

Transporter(l1::Matrix, l2::Matrix)=Transporter([l1],[l2])

"ratio of two vectors"
function ratio(v::Vector, w::Vector)
  i=findfirst(x->!iszero(x),w)
  if isnothing(i) return nothing end
  r=v[i]//w[i]
  if v!= r.* w  return nothing end
  r
end

end
