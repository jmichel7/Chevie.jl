"""
GLinearAlgebra: linear algebra over arbitrary fields and rings

The  linear  algebra  package  in  Julia  is  not  suitable  for  a general
mathematics  package: it assumes  the field is  the Real or Complex numbers
and uses floating point to do approximate computations.
Here we are interested in functions which work over any field (or sometimes
any ring) and assume exact computations.

For more information, look at the helpstrings of
`echelon!, exterior_power, comatrix, bigcell_decomposition, 
  ratio, charpoly, solutionmat, transporter, 
  permanent, symmetric_power, diagconj_elt, lnullspace,
  intersect_rowspace, in_rowspace, rowspace, independent_rows`.
"""
module GLinearAlgebra
using CyclotomicNumbers: Cyc # case in echelon
using Combinat: combinations, multisets, tally
using LinearAlgebra: tr, exactdiv, det_bareiss
using ..Util: toM, toL
using ..Tools: improve_type
export echelon!, exterior_power, comatrix, bigcell_decomposition, 
  ratio, charpoly, solutionmat, transporter, 
  permanent, symmetric_power, diagconj_elt, lnullspace,
  intersect_rowspace, in_rowspace, rowspace, independent_rows

"""
`echelon!(m::AbstractMatrix)`

puts  `m` in echelon  form in-place and  returns: 

(`m`, indices of linearly independent  rows of original `m`)

The  echelon form transforms the  rows of m into  a particular basis of the
rowspace. The first non-zero element of each line is 1, and such an element
is also the only non-zero in its column. works in any field.
"""
function echelon!(m::AbstractMatrix)
  rk=0
  inds=collect(axes(m,1))
  for k in axes(m,2)
    j=findfirst(!iszero,@view m[rk+1:end,k])
    if j===nothing continue end
    j+=rk # m[j,k]!=0
    rk+=1
    if rk!=j Base.swaprows!(m,j,rk)
      inds[j],inds[rk]=inds[rk],inds[j]
    end
    if !isone(m[rk,k])
      @views m[rk,:].*=inv(m[rk,k])
    end
    for j in axes(m,1)
      if rk!=j && !iszero(m[j,k]) 
        @views m[j,:].-=m[j,k].*m[rk,:]
      end
    end
#   println(m)
  end
  m,inds[1:rk]
end

"""
`rowspace(m::AbstractMatrix)`

returns a canonical form of the rowspace of `m` (the vector space generated
by  the rows of `m`),  the echelon form. This  is a particular basis of the
rowspace  of `m`: the first non-zero element of each line is 1, and such an
element  is also the only non-zero  element in its column. Integer matrices
are converted to `Rational` before applying this function.

```julia-repl
julia> m=[1 2;2 4;5 6]
3×2 Matrix{Int64}:
 1  2
 2  4
 5  6

julia> GLinearAlgebra.rowspace(m)
2×2 view(::Matrix{Rational{Int64}}, 1:2, :) with eltype Rational{Int64}:
 1//1  0//1
 0//1  1//1
```
"""
function rowspace(m::AbstractMatrix;dup=true)
  if dup m=copy(m) end
  m,inds=echelon!(m)
  @view m[1:length(inds),:]
end
rowspace(m::AbstractMatrix{<:Integer})=rowspace(m*1//1;dup=false)
rowspace(m::AbstractMatrix{<:Cyc{<:Integer}})=rowspace(m*1//1;dup=false)

"""
`independent_rows(m::AbstractMatrix)`

returns  a vector of the  indices of a set  of independent rows of `m` (the
lexicographically first such set).

```julia-repl
julia> m=[1 2;2 4;5 6]
3×2 Matrix{Int64}:
 1  2
 2  4
 5  6

julia> GLinearAlgebra.independent_rows(m)
2-element Vector{Int64}:
 1
 3
```
"""
function independent_rows(m::AbstractMatrix;dup=true)
  if dup m=copy(m) end
  last(echelon!(m))
end
independent_rows(m::AbstractMatrix{<:Integer})=
    independent_rows(m*1//1;dup=false)
independent_rows(m::AbstractMatrix{<:Cyc{<:Integer}})=
    independent_rows(m*1//1;dup=false)

"""
 `GLinearAlgebra.rank(m::AbstractMatrix)`
 computes exactly the rank of `m`

```julia-rep1
julia> hilbert(r)=[1//(n+m) for n in 1:r, m in 1:r];

julia> LinearAlgebra.rank(hilbert(11))
10

julia> GLinearAlgebra.rank(hilbert(11)) # correct value
11
```
"""
rank(m::AbstractMatrix)=length(independent_rows(m))

"""
`GLinearAlgebra.nullspace(m::AbstractMatrix)`

 computes the right nullspace of `m` in a type-preserving way.
 Is not exported to avoid conflict with LinearAlgebra
"""
function nullspace(m::AbstractMatrix)
  if eltype(m)<:Integer || eltype(m)<:Cyc{<:Integer} m=m//1 
  else m=copy(m)
  end
  m=echelon!(m)[1]
  n=size(m)[2]
  z=Int[]
  j=0
  lim=size(m,1)
  for i in axes(m,1)
    f=findfirst(!iszero,m[i,j+1:end])
    if f===nothing
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

"""
`lnullspace(m::AbstractMatrix)`

The left nullspace of `m`
"""
lnullspace(m::AbstractMatrix)=transpose(nullspace(transpose(m)))

"""
`intersect_rowspace(m::AbstractMatrix,n::AbstractMatrix)` 

The intersection of the rowspaces of `m` and `n`.
"""
function intersect_rowspace(m::AbstractMatrix,n::AbstractMatrix)
  mat=rowspace([m m;n zero(n)])
  in=mat[:,size(m,2).+axes(m,2)]
  in=in[count(!iszero,eachrow(mat[:,axes(m,2)]))+1:end,:]
  in[1:count(!iszero,eachrow(in)),:]
end

"""
`in_rowspace(v::AbstractVector,m::AbstractMatrix)` 

whether `v` is in the rowspace of `m`.
"""
in_rowspace(v::AbstractVector,m::AbstractMatrix)=solutionmat(m,v)!=nothing

# see Cohen computational algebraic number theory 2.2.7
function charpolyandcomatrix(m)
  C=one(m)
  res=C
  n=size(m,1)
  a=ones(eltype(m),n+1)
  for i in 1:n
    if i==n res=(-1)^(n-1)*C end
    C*=m
    a[n+1-i]=exactdiv(-tr(C),i)
    if i!=n C+=a[n+1-i]*one(C) end
  end
  a,res
end

"""
`charpoly(M::Matrix)` 

The characteristic polynomial of `M` (as a `Vector` of coefficients).
This function works over any ring.
"""
charpoly(m)=first(charpolyandcomatrix(m))

"""
`comatrix(M::Matrix)` 

is defined by `comatrix(M)*M=det(M)*one(M)`.
This function works over any ring.
"""
comatrix(m)=last(charpolyandcomatrix(m))

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
julia> @Pol q
Pol{Int64}: q

julia> M=[q^6 q^0 q^3 q^3 q^5+q q^4+q^2; q^0 q^6 q^3 q^3 q^5+q q^4+q^2; q^3 q^3 q^6 q^0 q^4+q^2 q^5+q; q^3 q^3 q^0 q^6 q^4+q^2 q^5+q; q^5+q q^5+q q^4+q^2 q^4+q^2 q^6+q^4+q^2+1 q^5+2*q^3+q; q^4+q^2 q^4+q^2 q^5+q q^5+q q^5+2*q^3+q q^6+q^4+q^2+1]
6×6 Matrix{Pol{Int64}}:
 q⁶     1      q³     q³     q⁵+q        q⁴+q²
 1      q⁶     q³     q³     q⁵+q        q⁴+q²
 q³     q³     q⁶     1      q⁴+q²       q⁵+q
 q³     q³     1      q⁶     q⁴+q²       q⁵+q
 q⁵+q   q⁵+q   q⁴+q²  q⁴+q²  q⁶+q⁴+q²+1  q⁵+2q³+q
 q⁴+q²  q⁴+q²  q⁵+q   q⁵+q   q⁵+2q³+q    q⁶+q⁴+q²+1

julia> bb=[[2],[4],[6],[3,5],[1]];

julia> (P,L)=bigcell_decomposition(M,bb);

julia> P
6×6 Matrix{Pol{Int64}}:
 1    0  0    0    0        0
 q⁻⁶  1  q⁻³  q⁻³  q⁻¹+q⁻⁵  q⁻²+q⁻⁴
 0    0  1    0    0        0
 q⁻³  0  0    1    q⁻²      q⁻¹
 q⁻¹  0  0    0    1        0
 q⁻²  0  q⁻¹  0    q⁻¹      1

julia> L
6×6 Matrix{Pol{Int64}}:
 q⁶-q⁴-1+q⁻²  0   0            0     0            0
 0            q⁶  0            0     0            0
 0            0   q⁶-q⁴-1+q⁻²  0     0            0
 0            0   0            q⁶-1  0            0
 0            0   0            0     q⁶-q⁴-1+q⁻²  0
 0            0   0            0     0            q⁶-1

julia> M==transpose(P)*L*P
true
```
"""
function bigcell_decomposition(M,b=map(i->[i],axes(M,1)))
  L=one(M)
  P=one(M)
  block(X,i,j)=X[b[i],b[j]]
  if M==transpose(M)
    for j in eachindex(b)
      L[b[j],b[j]]=block(M,j,j)
      if j>1 L[b[j],b[j]]-=sum(k->block(P,j,k)*block(L,k,k)*
                              transpose(block(P,j,k)),1:j-1) end
      cb=comatrix(block(L,j,j))
      db=det_bareiss(block(L,j,j))
      for i in j+1:length(b)
        P[b[i],b[j]]=block(M,i,j)
        if j>1 P[b[i],b[j]]-=
          sum(k->block(P,i,k)*block(L,k,k)*transpose(block(P,j,k)),1:j-1)
        end
        P[b[i],b[j]]=improve_type((block(P,i,j)*cb).//db)
      end
    end
    return permutedims(P),L
  end
  tP=one(M)
  for j in eachindex(b)
    L[b[j],b[j]]=block(M,j,j)-sum(k->block(P,j,k)*block(L,k,k)*
                                  transpose(block(P,j,k)),1:j-1)
    cb=comatrix(block(L,j,j))
    db=det_bareiss(block(L,j,j))
    for i in j+1:length(b)
      P[b[i],b[j]]=block(M,i,j)-sum(k->block(P,i,k)*block(L,k,k)*
                                     block(tP,k,j),1:j-1)
      P[b[i],b[j]]=improve_type((block(P,i,j)*cb).//db)
      tP[b[j],b[i]]=cb*(block(M,j,i)-sum(k->block(P,j,k)*
                                         block(L,k,k)*block(tP,k,i),1:j-1))
      tP[b[i],b[j]]=improve_type(block(tP,i,j).//db)
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
4×4 Matrix{Int64}:
 1  2  3  4
 2  3  4  1
 3  4  1  2
 4  1  2  3

julia> exterior_power(M,2)
6×6 Matrix{Int64}:
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
  [det_bareiss(A[i,j]) for i in basis, j in basis]
end

"""
`permanent(m)`

returns the *permanent* of the square matrix `m`, which is defined by 
``\\sum_{p\\in\\frak S_n}\\prod_{i=1}^n m[i,p(i)]``.

Note the similarity of the definition of  the permanent to the definition
of the determinant.  In  fact the only  difference is the missing sign of
the permutation.  However the  permanent is quite unlike the determinant,
for example   it is  not  multilinear or  alternating.  It   has  however
important combinatorical properties.

```julia-repl
julia> permanent([0 1 1 1;1 0 1 1;1 1 0 1;1 1 1 0])
9 # inefficient way to compute the number of derangements of 1:4

julia> permanent([1 1 0 1 0 0 0; 0 1 1 0 1 0 0;0 0 1 1 0 1 0; 0 0 0 1 1 0 1;1 0 0 0 1 1 0;0 1 0 0 0 1 1;1 0 1 0 0 0 1])
24 # 24 permutations fit the projective plane of order 2 
```
"""
function permanent(m)
  n=size(m,1)
  function Permanent2(m,i,sum)
    if i==n+1 reduce(*,sum;init=1)
    else Permanent2(m,i+1,sum)-Permanent2(m,i+1,sum+m[i,:])
    end
  end 
  (-1)^n*Permanent2(m,1,zeros(Int,n));
end

"""
`symmetric_power(m,n)`

returns the `n`-th symmetric power of the square matrix `m`, in the basis 
naturally indexed by the `multisets` of `1:n`, where `n=size(m,1)`.

```julia-repl
julia> m=[1 2;3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> Int.(symmetric_power(m,2))
3×3 Matrix{Int64}:
 1   2   4
 6  10  16
 9  12  16
```
"""
function symmetric_power(A,m)
  f(j)=prod(factorial,last.(tally(j)))
  basis=multisets(axes(A,1),m)
  [permanent(A[i,j])//f(i) for i in basis, j in basis]
end

"""
`transporter(l1, l2 )`

`l1`  and `l2` should be vectors of  the same length of square matrices all
of the same size. The result is a basis of the vector space of matrices `A`
such  that for any `i` we have  `A*l1[i]=l2[i]*A` --- the basis is returned
as  a vector of matrices, empty if the vector space is 0. This is useful to
find whether two representations are isomorphic.
"""
function transporter(l1::Vector{<:Matrix}, l2::Vector{<:Matrix})
  n=size(l1[1],1)
  M=Vector{eltype(l1[1])}[]
  for i in 1:n, j in 1:n,  m in 1:length(l1)
    v=zero(l1[1])
    v[i,:]+=l1[m][:,j]
    v[:,j]-=l2[m][i,:]
    push!(M, vec(v))
  end
  M=rowspace(toM(M))
  M=M[filter(i->!iszero(M[i,:]),axes(M,1)),:]
  v=nullspace(M)
  if isempty(v) return v end
  map(w->reshape(w,(n,n)), eachcol(v))
end

transporter(l1::Matrix, l2::Matrix)=transporter([l1],[l2])

"""
`ratio(v::abstractVector,w::abstractVector)` 

ratio `v/w`, `nothing` if `v` is not a multiple of `w`.
"""
function ratio(v::AbstractVector, w::AbstractVector)
  i=findfirst(x->!iszero(x),w)
  if i===nothing return nothing end
  r=v[i]//w[i]
  if v==r.*w return r end
end

"""
`solutionmat(mat,v)`

returns  one solution of the  equation `transpose(mat)*x=v` or `nothing` if
no  such  solution  exists.  Similar  to  `transpose(mat)\\v` when `mat` is
invertible.

```julia-repl
julia> m=[2 -4 1;0 0 -4;1 -2 -1]
3×3 Matrix{Int64}:
 2  -4   1
 0   0  -4
 1  -2  -1

julia> x=solutionmat(m,[10,-20, -10])
3-element Vector{Rational{Int64}}:
  5//1
 15//4
  0//1

julia> m'*x
3-element Vector{Rational{Int64}}:
  10//1
 -20//1
 -10//1

julia> solutionmat(m,[10, 20, -10])
```
"""
function solutionmat(m::AbstractMatrix,v::AbstractVector)
  m=transpose(m).*1//1
  if length(v)!=size(m,1) error("dimension mismatch") end
  v=v.*1//1
  r=0
  c=1
  while c<=size(m,2) && r<size(m,1)
    s=findfirst(!iszero,m[r+1:end,c])
    if s!==nothing
      s+=r
      r+=1
      piv=inv(m[s,c])
      m[s,:],m[r,:]=m[r,:],m[s,:].*piv
      v[s],v[r]=v[r],v[s]*piv
      for s in 1:size(m,1)
        if s!=r && !iszero(m[s,c])
          v[s]-=m[s,c]*v[r]
  @views  m[s,:].-=m[s,c]*m[r,:]
        end
      end
    end
    c+=1
  end
  if any(!iszero,v[r+1:end]) return nothing end
  h=eltype(m)[]
  s=size(m,2)
  z=zero(eltype(m))
  r=1
  c=1
  while c<=s && r<=size(m,1)
    while c<=s && iszero(m[r,c])
      c+=1
      push!(h, z)
    end
    if c<=s
      push!(h, v[r])
      r+=1
      c+=1
    end
  end
  while c<=s
    push!(h,z)
    c+=1
  end
  return h
end

"""
`solutionmat(m,n::AbstractMatrix)`

return a matrix `x` such that `x*m==n`. This is interesting when `m` 
is not invertible.
"""
solutionmat(m,n::AbstractMatrix)=toM(solutionmat.(Ref(m),eachrow(n)))

"""
`diagconj_elt(M,N)`

`M` and `N` must be square matrices of the same size. This function returns
a  list `d`  such that  `N==inv(Diagonal(d))*M*Diagonal(d)` if  such a list
exists, and `nothing` otherwise.

```julia_repl
julia> M=[1 2;2 1];N=[1 4;1 1]
2×2 Matrix{Int64}:
 1  4
 1  1

julia> diagconj_elt(M,N)
2-element Vector{Rational{Int64}}:
 1//1
 2//1
```
"""
function diagconj_elt(M, N)
  d=fill(zero(eltype(M))*1//1,size(M,1));d[1]=1
  n=size(M,1)
  for i in 1:n, j in i+1:n
    if M[i,j]!=0
      if N[i,j]==0 return nothing end
      if d[i]!=0
        c=d[i]*N[i,j]//M[i,j]
        if d[j]!=0 if c!=d[j] return nothing end 
        else d[j]=c end
      end
    end
  end
  if 0 in d return nothing end
  if !all([N[i,j]==M[i,j]*d[j]//d[i] for i in 1:n for j in 1:n])
    return nothing end
  d
end

end
