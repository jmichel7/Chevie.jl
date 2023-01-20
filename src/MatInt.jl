"""
Smith and Hermite normal forms for integral matrices, ported from GAP4, and
Diaconis-Graham  normal form  for sets  of generators  of an abelian group,
ported from GAP3/Chevie.

The  code ported  from GAP4  for NormalFormIntMat  is horrible (unreadable)
like the original one.

The  best way to  make sure of  the validity of  the result is to work with
matrices of SaferIntegers.

For  the API, look at the docstrings for `smith, smith_transforms, hermite,
hermite_transforms,  col_hermite,  col_hermite_transforms, diaconis_graham,
baseInt, complementInt, lnullspaceInt, solutionmatInt, intersect_rowspaceInt`
"""
module MatInt

export complementInt, lnullspaceInt, solutionmatInt, smith, smith_transforms,
  hermite, hermite_transforms, col_hermite, col_hermite_transforms, 
  diaconis_graham, baseInt, intersect_rowspaceInt

using LinearAlgebra: I, dot

"`prime_part(N,a)`  largest factor of `N` prime to `a`"
function prime_part(N, a)
  while true
    a=gcd(a, N)
    if a==1 return N end
    N=div(N, a)
  end
end

"`rgcd(N,a)` smallest `c≥0` such that `gcd(N,a+c)==1`"
function rgcd(N, a)
  if N==1 return 0 end
  r=[mod(a-1, N)]
  d=[N]
  c=0
  while true
    for i in eachindex(r) r[i]=mod(r[i]+1,d[i]) end
    if all(>(0),r)
      g=1
      i=0
      while g==1 && i<length(r)
        i+=1
        g=gcd(r[i], d[i])
      end
      if g==1 return c end
      q=prime_part(div(d[i], g), g)
      if q>1
       push!(r,mod(r[i], q))
       push!(d, q)
      end
      r[i]=0
      d[i]=g
    end
    c+=1
  end
end

"""
`Gcdex(m,n)`

`Gcdex`  returns a  named tuple  with fields  `gcd=gcd(m,n)` and `coeff`, a
unimodular (invertible integer) 2x2 matrix such that `coeff*[m,n]=[gcd,0]`.

If `m*n!=0`, `abs(coeff[1,1])≤abs(n)/(2*gcd)` and
`abs(coeff[1,2])≤abs(m)/(2*gcd)`.   If  `m`  and  `n`  are  not  both  zero
`coeff[2,1]` is `-n/gcd` and `coeff[2,2]` is `m/gcd`.

```julia-repl
julia> MatInt.Gcdex(123,66)
(gcd = 3, coeff = [7 -13; -22 41])
# [7 -13;-22 41]*[123,66]==[3,0]
julia> MatInt.Gcdex(0,-3)
(gcd = 3, coeff = [0 -1; 1 0])
julia> MatInt.Gcdex(0,0)
(gcd = 0, coeff = [1 0; 0 1])
```
"""
function Gcdex( m, n )
  if 0<=m  f=m; fm=1
  else f=-m; fm=-1 end
  g=0<=n ? n : -n
  gm=0
  while g!=0
    q=div(f,g)
    h=g
    hm=gm
    g=f-q*g
    gm=fm-q*gm
    f=h
    fm=hm
  end
  (gcd=f, coeff= n==0 ? [fm 0; gm 1] : [fm div(f-fm*m,n); gm div(-gm*m,n)])
end

"""
`bezout(A)` 

`A` should be  a 2x2 matrix. Returns a `NamedTuple` with 2 fields
  - `.sign`  the sign of `det(A)`
  - `.rowtrans` such that `rowtrans*A=[e f;0 g]`  (Hermite normal form)
"""
function bezout(A)
  e=Gcdex(A[1,1],A[2,1])
  @views f,g=e.coeff*A[:,2]
  if iszero(g) return e end
  coeff=e.coeff
  if g<0
    @views coeff[2,:]=-coeff[2,:]
    g=-g
  end
  @views coeff[1,:].+=-div(f-mod(f,g),g).*coeff[2,:]
  (sign=sign(A[1,1]*A[2,2]-A[1,2]*A[2,1]),rowtrans=coeff)
end

"""
`mgcdex(N,a,v)` returns `c[1],c[2],…c[k]` such that
`gcd(N,a+c[1]*v[1]+…+c[n]*v[k])==gcd(N,a,v[1],v[2],…,v[k])`
"""
function mgcdex(N, a, v)
  l=length(v)
  h=N
  M=zeros(eltype(v),l)
  for i in 1:l
    g=h
    h=gcd(g, v[i])
    M[i]=div(g, h)
  end
  h=gcd(a,h)
  g=div(a,h)
  for i in l:-1:1
    b=div(v[i], h)
    d=prime_part(M[i], b)
    if d==1 c=0
    else 
      u=g//b
      c=rgcd(d, numerator(u)*invmod(denominator(u),d))
      g+=c*b
    end
    M[i]=c
  end
  M
end

## SNFofREF - fast SNF of REF matrix
function SNFofREF(R)
  n,m=size(R)
  piv=findfirst.(!iszero,eachrow(R))
  r=findfirst(isnothing,piv)
  if isnothing(r) r=length(piv)
  else
    r-=1
    piv.=@view piv[1:r]
  end
  append!(piv, setdiff(1:m, piv))
  T=zeros(eltype(R),n,m)
  for j in 1:m
    for i in 1:min(r,j) T[i,j]=R[i,piv[j]] end
  end
  si=1
  A=Vector{eltype(R)}(undef,n)
  d=2
  for k in 1:m
    if k<=r
      d*=abs(T[k,k])
      T[k,:].=mod.((@view T[k,:]), 2d)
    end
    t=min(k, r)
    for i in t-1:-1:si
      t=mgcdex(A[i], T[i,k], (T[i+1],T[k]))[1]
      if t!=0
        T[i,:].+=(@view T[i+1,:]).*t
        T[i,:].=mod.((@view T[i,:]), A[i])
      end
    end
    for i in si:min(k-1, r)
      g=gcdx(A[i], T[i,k])
      T[i,k]=0
      if g[1]!=A[i]
        b=div(A[i], g[1])
        A[i]=g[1]
        for ii in i+1:min(k-1,r)
          T[ii,:].+=mod.((@view T[i,:])*(-g[3]*div(T[ii,k],A[i])),A[ii])
          T[ii,k]*=b
          T[ii,:].=mod.((@view T[ii,:]),A[ii])
        end
        if k<=r
          t=g[3]*div(T[k,k], g[1])
          T[k,:].+=-t*@view T[i,:]
          T[k,k]*=b
        end
        T[i,:].=mod.(T[i,:], A[i])
        if A[i]==1 si=i+1 end
      end
    end
    if k<=r
      A[k]=abs(T[k,k])
      T[k,:].=mod.((@view T[k,:]), A[k])
    end
  end
  for i in 1:r T[i,i]=A[i] end
  T
end

"""
general operation for computation of various Normal Forms.

Options:
  - TRIANG Triangular Form / Smith Normal Form.
  - REDDIAG Reduce off diagonal entries.
  - ROWTRANS Row Transformations.
  - COLTRANS Col Transformations.

Compute  a Triangular, Hermite or  Smith form of the  `n × m` integer input
matrix `A`. Optionally, compute `n × n` and `m × m` unimodular transforming
matrices `Q, P` which satisfy `Q A==H` or `Q A P==S`.

Compute a Triangular, Hermite or Smith form of the n x m 
integer input matrix A.  Optionally, compute n x n / m x m
unimodular transforming matrices which satisfy Q C A==H 
or  Q C A B P==S.

Triangular / Hermite :

Let I be the min(r+1,n) x min(r+1,n) identity matrix with r=rank(A).
Then Q and C can be written using a block decomposition as

             [ Q1 |   ]  [ C1 | C2 ]
             [----+---]  [----+----]  A== H.
             [ Q2 | I ]  [    | I  ]

Smith :

  [ Q1 |   ]  [ C1 | C2 ]     [ B1 |   ]  [ P1 | P2 ]
  [----+---]  [----+----]  A  [----+---]  [----+----] ==S.
  [ Q2 | I ]  [    | I  ]     [ B2 | I ]  [ *  | I  ]

 * - possible non-zero entry in upper right corner...

The routines used are based on work by Arne Storjohann and were implemented
in GAP4 by him and R.Wainwright.

Returns a Dict with entry `:normal` containing the computed normal form and
optional  entries `:rowtrans` and/or `:coltrans`  which hold the respective
transformation matrix. Also in the dict are entries holding the sign of the
determinant if A is square, `:signdet`, and the rank of the matrix, `:rank`.

```julia-repl
julia> m=[1 15 28;4 5 6;7 8 9]
3×3 Matrix{Int64}:
 1  15  28
 4   5   6
 7   8   9

julia> MatInt.NormalFormIntMat(m,REDDIAG=true,ROWTRANS=true)
Dict{Symbol, Any} with 6 entries:
  :rowQ     => [-2 62 -35; 1 -30 17; -3 97 -55]
  :normal   => [1 0 1; 0 1 1; 0 0 3]
  :rowC     => [1 0 0; 0 1 0; 0 0 1]
  :rank     => 3
  :signdet  => 1
  :rowtrans => [-2 62 -35; 1 -30 17; -3 97 -55]

julia> r=MatInt.NormalFormIntMat(m,TRIANG=true,ROWTRANS=true,COLTRANS=true)
Dict{Symbol, Any} with 9 entries:
  :rowQ     => [-2 62 -35; 1 -30 17; -3 97 -55]
  :normal   => [1 0 0; 0 1 0; 0 0 3]
  :colQ     => [1 0 -1; 0 1 -1; 0 0 1]
  :coltrans => [1 0 -1; 0 1 -1; 0 0 1]
  :rowC     => [1 0 0; 0 1 0; 0 0 1]
  :rank     => 3
  :signdet  => 1
  :rowtrans => [-2 62 -35; 1 -30 17; -3 97 -55]
  :colC     => [1 0 0; 0 1 0; 0 0 1]

julia> r[:rowtrans]*m*r[:coltrans]
3×3 Matrix{Int64}:
 1  0  0
 0  1  0
 0  0  3
```
"""
function NormalFormIntMat(mat::AbstractMatrix; TRIANG=false, REDDIAG=false, ROWTRANS=false, COLTRANS=false)
# The gap code for INPLACE cannot work -- different memory model to julia
  sig=1
  #Embed nxm mat in a (n+2)x(m+2) larger "id" matrix
  n,m=size(mat).+(2,2)
  A=zeros(eltype(mat),n,m)
  A[2:end-1,2:end-1]=mat
  A[1,1]=1
  A[n,m]=1
  if ROWTRANS
    Q=zeros(eltype(mat),n,n)
    Q[1,1]=1
    C=one(Q)
  end
  if TRIANG && COLTRANS
    B=one(zeros(eltype(mat),m,m))
    P=copy(B)
  end
  r=0
  c2=1
  rp=Int[]
  while m>c2
    c1=c2
    push!(rp,c1)
    r+=1
    if ROWTRANS Q[r+1,r+1]=1 end
    j=c1+1
    k=0
    while j<=m
      k=r+1
      while k<=n && A[r,c1]*A[k,j]==A[k,c1]*A[r,j] k+=1 end
      if k<=n
        c2=j
        j=m
      end
      j+=1
    end
    #Smith with some transforms..
    if TRIANG && ((COLTRANS || ROWTRANS) && c2<m)
      N=gcd(@view A[r:n,c2])
      for j in Iterators.flatten((c1+1:c2-1,c2+1:m-1,c2))
        if j==c2
          b=A[r,c2]
          a=A[r,c1]
          for i in r+1:n
            if b!=1
              g=gcdx(b, A[i,c2])
              b=g[1]
              a=g[2]*a+g[3]*A[i,c1]
            end
          end
          N=0
          T=typeof(N)
          for i in r:n
            if N!=1 N=gcd(N, A[i,c1]-div(A[i,c2],b)*widen(a)) end
          end
          N=T(N)
        else
          c=mgcdex(N, A[r,j], @view A[r+1:n,j])
          c=dot(c,@view A[r+1:n,j])
          b=A[r,j]+c
          a=A[r,c1]+c
        end
        t=mgcdex(N, a, (b,))[1]
        tmp=A[r,c1]+t*A[r,j]
        if tmp==0 || tmp*A[k,c2]==(A[k,c1]+t*A[k,j])*A[r,c2]
          t+=1+mgcdex(N, a+t*b+b,(b,))[1]
        end
        if t>0
        @views  A[:,c1].+=t*A[:,j]
          if COLTRANS B[j,c1]+=t end
        end
      end
      if A[r,c1]*A[k,c1+1]==A[k,c1]*A[r,c1+1]
        @views A[:,c1+1].+=A[:,c2]
        if COLTRANS B[c2,c1+1]=1 end
      end
      c2=c1+1
    end
    c=mgcdex(abs(A[r,c1]), A[r+1,c1], @view A[r+2:n,c1])
    for i in r+2:n
      if c[i-r-1]!=0
        @views A[r+1,:].+=c[i-r-1].*A[i,:]
        if ROWTRANS
          C[r+1,i]=c[i-r-1]
          @views Q[r+1,:].+=c[i-r-1].*Q[i,:]
        end
      end
    end
    i=r+1
    while A[r,c1]*A[i,c2]==A[i,c1]*A[r,c2] i+=1 end
    if i>r+1
      c=mgcdex(abs(A[r,c1]), A[r+1,c1]+A[i,c1], (A[i],A[c1]))[1]+1
      @views A[r+1,:].+=c.*A[i,:]
      if ROWTRANS
        C[r+1,i]+=c
        @views Q[r+1,:].+=c.*Q[i,:]
      end
    end
    g=bezout(@view A[r:r+1,[c1,c2]])
    sig*=g.sign
    @views A[r:r+1,:].=g.rowtrans*A[r:r+1,:]
    if ROWTRANS @views Q[r:r+1,:].=g.rowtrans*Q[r:r+1,:] end
    for i in r+2:n
      q=div(A[i,c1], A[r,c1])
      @views A[i,:].-=q.*A[r,:]
      if ROWTRANS @views Q[i,:].-=q.*Q[r,:] end
      q=div(A[i,c2], A[r+1,c2])
      @views A[i,:].-=q.*A[r+1,:]
      if ROWTRANS @views Q[i,:].-=q.*Q[r+1,:] end
    end
  end
  push!(rp,m) # length(rp)==r+1
  if n==m && r+1<n sig=0 end
  #smith w/ NO transforms - farm the work out...
  if TRIANG && !(ROWTRANS || COLTRANS)
    A=@view A[2:end-1,2:end-1]
    R=Dict(:normal => SNFofREF(A), :rank=>r-1)
    if n==m R[:signdet]=sig end
    return R
  end
  # hermite or (smith w/ column transforms)
  if (!TRIANG && REDDIAG) || (TRIANG && COLTRANS)
    for i in r:-1:1
      for j in i+1:r+1
        q=div(A[i,rp[j]]-mod(A[i,rp[j]], A[j,rp[j]]), A[j,rp[j]])
        @views A[i,:].-=q.*A[j,:]
        if ROWTRANS @views Q[i,:].-=q.*Q[j,:] end
      end
      if TRIANG && i<r
        for j in i+1:m
          q=div(A[i,j], A[i,i])
          @views A[1:i,j].-=q.*A[1:i,i]
          if COLTRANS P[i,j]=-q end
        end
      end
    end
  end
  #Smith w/ row but not col transforms
  if TRIANG && ROWTRANS && !COLTRANS
    for i in 1:r-1
      t=A[i,i]
      A[i,:].=0
      A[i,i]=t
    end
    for j in r+1:m-1
      A[r,r]=gcd(A[r,r], A[r,j])
      A[r,j]=0
    end
  end
  #smith w/ col transforms
  if TRIANG && COLTRANS && r<m-1
    c=mgcdex(A[r,r], A[r,r+1], @view A[r,r+2:m-1])
    for j in r+2:m-1
      A[r,r+1]+=c[j-r-1]*A[r,j]
      B[j,r+1]=c[j-r-1]
      @views P[1:r,r+1].+=c[j-r-1].*P[1:r,j]
    end
    P[r+1,:].=0
    P[r+1,r+1]=1
    g=Gcdex(A[r,r], A[r,r+1])
    A[r,r]=g.gcd
    A[r,r+1]=0
    @views P[1:r+1,r:r+1]*=transpose(g.coeff)
    for j in r+2:m-1
      q=div(A[r,j], A[r,r])
      @views P[1:r+1,j].-=q.*P[1:r+1,r]
      A[r,j]=0
    end
    P[r+2:m-1,:].=0
    for i in r+2:m-1 P[i,i]=1 end
  end
  #row transforms finisher
  if ROWTRANS for i in r+2:n Q[i,i]=1 end end
  R=Dict{Symbol,Any}(:normal => A[2:end-1,2:end-1],
                     :rank => r-1,:signdet=>n==m ? sig : nothing)
  if ROWTRANS
    R[:rowC]=C[2:end-1,2:end-1]
    R[:rowQ]=Q[2:end-1,2:end-1]
    R[:rowtrans]=R[:rowQ]*R[:rowC]
  end
  if TRIANG && COLTRANS
    R[:colC]=B[2:end-1,2:end-1]
    R[:colQ]=P[2:end-1,2:end-1]
    R[:coltrans]=R[:colC]*R[:colQ]
  end
  R
end

"""
`TriangulizedIntegerMat(mat)`
Changes  `mat` to be in  upper triangular form.

```julia-repl
julia> m=[1 15 28;4 5 6;7 8 9]
3×3 Matrix{Int64}:
 1  15  28
 4   5   6
 7   8   9

julia> MatInt.TriangulizedIntegerMat(m)
3×3 Matrix{Int64}:
 1  15  28
 0   1   1
 0   0   3
```
"""
TriangulizedIntegerMat(mat)=NormalFormIntMat(mat;)[:normal]

"""
```julia-repl
julia> m=[1 15 28;4 5 6;7 8 9]
3×3 Matrix{Int64}:
 1  15  28
 4   5   6
 7   8   9

julia> n=MatInt.TriangulizedIntegerMatTransform(m)
Dict{Symbol, Any} with 6 entries:
  :rowQ     => [1 0 0; 1 -30 17; -3 97 -55]
  :normal   => [1 15 28; 0 1 1; 0 0 3]
  :rowC     => [1 0 0; 0 1 0; 0 0 1]
  :rank     => 3
  :signdet  => 1
  :rowtrans => [1 0 0; 1 -30 17; -3 97 -55]


julia> n[:rowtrans]*m==n[:normal]
true
```
"""
TriangulizedIntegerMatTransform(mat)=NormalFormIntMat(mat;ROWTRANS=true)

"""
`hermite(m::AbstractMatrix{<:Integer})`

returns  the row Hermite normal  form `H` of `m`,  a row equivalent integer
upper  triangular form. If a pivot is the  first non-zero entry on a row of
`H`,  the quadrant  below left  a pivot  is zero,  pivots are  positive and
entries  above a  pivot are  nonnegative and  smaller than the pivot. There
exists a unique invertible integer matrix `r` such that `rm==H`.

```julia-repl
julia> m=[1 15 28;4 5 6;7 8 9]
3×3 Matrix{Int64}:
 1  15  28
 4   5   6
 7   8   9

julia> hermite(m)
3×3 Matrix{Int64}:
 1  0  1
 0  1  1
 0  0  3
```
"""
function hermite(mat::AbstractMatrix{<:Integer})
  NormalFormIntMat(mat;REDDIAG=true)[:normal]
end

"""
`hermite_transforms(m::AbstractMatrix{<:Integer})`

The  row Hermite  normal form  `H` of  the `m`  is a row equivalent integer
upper  triangular form. If a pivot is the  first non-zero entry on a row of
`H`,  the quadrant  below left  a pivot  is zero,  pivots are  positive and
entries  above a  pivot are  nonnegative and  smaller than the pivot. There
exists   a  unique  invertible  integer   matrix  `r`  such  that  `rm==H`.
`hermite_transforms`  returns  a  tuple  with  components  `.normal=H`  and
`.rowtrans=r`.

```julia-repl
julia> m=[1 15 28;4 5 6;7 8 9]
3×3 Matrix{Int64}:
 1  15  28
 4   5   6
 7   8   9

julia> n=hermite_transforms(m)
(normal = [1 0 1; 0 1 1; 0 0 3], rowtrans = [-2 62 -35; 1 -30 17; -3 97 -55], rank = 3, signdet = 1)

julia> n.rowtrans*m==n.normal
true
```
"""
function hermite_transforms(mat::AbstractMatrix{<:Integer})
  res=NormalFormIntMat(mat;REDDIAG=true,ROWTRANS=true)
  (normal=res[:normal], rowtrans=res[:rowtrans], 
   rank=res[:rank], signdet=res[:signdet])
end

"""
`col_hermite(m::AbstractMatrix{<:Integer})`

returns  the column Hermite  normal form `H`  of the integer  matrix `m`, a
column  equivalent lower triangular form. If  a pivot is the first non-zero
entry on a column of `H` (the quadrant above right a pivot is zero), pivots
are  positive and entries left of a  pivot are nonnegative and smaller than
the  pivot. There exists  a unique invertible  integer matrix `Q` such that
`mQ==H`.

```julia-repl
julia> m=[1 15 28;4 5 6;7 8 9]
3×3 Matrix{Int64}:
 1  15  28
 4   5   6
 7   8   9

julia> col_hermite(m)
3×3 Matrix{Int64}:
 1  0  0
 0  1  0
 0  1  3
```
"""
function col_hermite(mat::AbstractMatrix{<:Integer})
  permutedims(NormalFormIntMat(transpose(mat);REDDIAG=true)[:normal])
end

"""
`col_hermite_transforms(m::AbstractMatrix{<:Integer})`

The  column Hermite normal form  `H` of the integer  matrix `m` is a column
equivalent lower triangular form. If a pivot is the first non-zero entry on
a  column of  `H` (the  quadrant above  right a  pivot is zero), pivots are
positive  and entries left of a pivot  are nonnegative and smaller than the
pivot.  There  exists  a  unique  invertible  integer  matrix `Q` such that
`mQ==H`. The function returns a named tuple with components `.normal=H` and
`.rowtrans=Q`.

```julia-repl
julia> m=[1 15 28;4 5 6;7 8 9]
3×3 Matrix{Int64}:
 1  15  28
 4   5   6
 7   8   9

julia> n=col_hermite_transforms(m)
(normal = [1 0 0; 0 1 0; 0 1 3], coltrans = [-1 13 -50; 2 -27 106; -1 14 -55], rank = 3, signdet = 1)

julia> m*n.coltrans==n.normal
true
```
"""
function col_hermite_transforms(mat::AbstractMatrix{<:Integer})
  res=NormalFormIntMat(transpose(mat);REDDIAG=true,ROWTRANS=true)
  (normal=permutedims(res[:normal]), coltrans=permutedims(res[:rowtrans]), 
   rank=res[:rank], signdet=res[:signdet])
end

"""
`smith(m::AbstractMatrix{<:Integer})`

computes  the Smith normal form  `S` of `m`, the  unique equivalent (in the
sense  that  there  exist  unimodular  integer  matrices  `P,  Q` such that
`Q*m*P==S`) diagonal matrix such that `Sᵢ,ᵢ` divides `Sⱼ,ⱼ` for `i≤j`.

```julia-repl
julia> m=[1 15 28 7;4 5 6 7;7 8 9 7]
3×4 Matrix{Int64}:
 1  15  28  7
 4   5   6  7
 7   8   9  7

julia> smith(m)
3×4 Matrix{Int64}:
 1  0  0  0
 0  1  0  0
 0  0  3  0
```
"""
smith(mat::AbstractMatrix{<:Integer})=NormalFormIntMat(mat,TRIANG=true)[:normal]

"""
`smith_transforms(m::AbstractMatrix{<:Integer})`

The  Smith normal form of `m` is  the unique equivalent diagonal matrix `S`
such  that `Sᵢ,ᵢ` divides `Sⱼ,ⱼ` for  `i≤j`. There exist unimodular integer
matrices  `P, Q` such that `Q*m*P==S`.  This function returns a named tuple
with `.normal=S`, `.rowtrans=Q` and `.coltrans=P`.

```julia-repl
julia> m=[1 15 28 7;4 5 6 7;7 8 9 7]
3×4 Matrix{Int64}:
 1  15  28  7
 4   5   6  7
 7   8   9  7

julia> n=smith_transforms(m)
(normal = [1 0 0 0; 0 1 0 0; 0 0 3 0], coltrans = [1 0 -1 -84; 0 1 -1 175; 0 0 1 -91; 0 0 0 1], rowtrans = [-2 62 -35; 1 -30 17; -3 97 -55], rank = 3, signdet = nothing)

julia> n.rowtrans*m*n.coltrans==n.normal
true
```
"""
function smith_transforms(mat::AbstractMatrix{<:Integer})
  res=NormalFormIntMat(mat;TRIANG=true,ROWTRANS=true,COLTRANS=true)
  (normal=res[:normal], coltrans=res[:coltrans], rowtrans=res[:rowtrans],
   rank=res[:rank], signdet=res[:signdet])
end

"""
`baseInt(m::Matrix{<:Integer})`

returns  a list of vectors that forms a  basis of the integral row space of
`m`, i.e. of the set of integral linear combinations of the rows of `m`.

```julia-repl
julia> m=[1 2 7;4 5 6;10 11 19]
3×3 Matrix{Int64}:
  1   2   7
  4   5   6
 10  11  19

julia> baseInt(m)
3×3 Matrix{Int64}:
 1  2   7
 0  3   7
 0  0  15
```
"""
function baseInt(mat::Matrix{<:Integer})
  norm=NormalFormIntMat(mat;REDDIAG=true)
  norm[:normal][1:norm[:rank],:]
end

"""
`intersect_rowspaceInt(m::Matrix{<:Integer}, n::Matrix{<:Integer})`

returns  a  matrix  whose  rows  forms  a  basis of the intersection of the
integral row spaces of `m` and `n`.

```julia-repl
julia> mat=[1 2 7;4 5 6;10 11 19]; nat=[5 7 2;4 2 5;7 1 4]
3×3 Matrix{Int64}:
 5  7  2
 4  2  5
 7  1  4

julia> intersect_rowspaceInt(mat,nat)
3×3 Matrix{Int64}:
 1  5  509
 0  6  869
 0  0  960
```
"""
function intersect_rowspaceInt(M1::Matrix{<:Integer}, M2::Matrix{<:Integer})
  M=vcat(M1, M2)
  r=TriangulizedIntegerMatTransform(M)
  T=r[:rowtrans][r[:rank]+1:size(M,1),axes(M1,1)]
  if !isempty(T) T*=M1 end
  baseInt(T)
end

"""
`complementInt(full::Matrix{<:Integer}, sub::Matrix{<:Integer})`

Let  `M` be the integral row module of  `full` and let `S`, a submodule of
`M`,  be the integral row  module of `sub`. This  function computes a free
basis for `M` that extends `S`, that is, if the dimension of `S` is `n` it
determines  a basis  `B={b₁,…,bₘ}` for  `M`, as  well as `n` integers `xᵢ`
such that the `n` vectors `sᵢ:=xᵢ⋅bᵢ` form a basis for `S`.

It returns a named tuple with the following components:
  - `complement` a matrix whose lines are `bₙ₊₁,…,bₘ`.
  - `sub` a matrix whose lines are the `sᵢ` (a basis for `S`).
  - `moduli` the factors `xᵢ`.

```julia-repl
julia> m=one(zeros(Int,3,3))
3×3 Matrix{Int64}:
 1  0  0
 0  1  0
 0  0  1

julia> n=[1 2 3;4 5 6]
2×3 Matrix{Int64}:
 1  2  3
 4  5  6

julia> complementInt(m,n)
(complement = [0 0 1], sub = [1 2 3; 0 3 6], moduli = [1, 3])
```
"""
function complementInt(full::Matrix{<:Integer}, sub::Matrix{<:Integer})
  F=baseInt(full)
  if isempty(sub) || iszero(sub) return (complement=F,sub=sub,moduli=Int[]) end
  S=intersect_rowspaceInt(F, sub)
  if S!=baseInt(sub) error(sub," must be submodule of ",full) end
  M=vcat(F,S)
  T=Integer.(inv(Rational.(TriangulizedIntegerMatTransform(M)[:rowtrans])))
  T=T[size(F,1)+1:end,axes(F,1)]
  r=smith_transforms(T)
  M=Integer.(inv(Rational.(r.coltrans))*F)
  (complement=baseInt(M[1+r.rank:end,:]), sub=r.rowtrans*T*F, 
   moduli=map(i->r.normal[i,i],1:r.rank))
end

"""
`lnullspaceInt(m::Matrix{<:Integer})

returns  a matrix whose rows form a  basis of the integral lnullspace of
`m`,  that is of elements  of the left nullspace  of `m` that have integral
entries.

```julia-repl
julia> m=[1 2 7;4 5 6;7 8 9;10 11 19;5 7 12]
5×3 Matrix{Int64}:
  1   2   7
  4   5   6
  7   8   9
 10  11  19
  5   7  12

julia> MatInt.lnullspaceInt(m)
2×5 Matrix{Int64}:
 1  18   -9  2  -6
 0  24  -13  3  -7
```
"""
function lnullspaceInt(mat)
  norm=TriangulizedIntegerMatTransform(mat)
  baseInt(norm[:rowtrans][norm[:rank]+1:size(mat,1),:])
end

"""
`solutionmatInt(mat::Matrix{<:Integer}, v::Vector{<:Integer})`

returns  a  vector  `x`  with  integer  entries  that  is a solution of the
equation `mat'*x=vec`. It returns `false` if no such vector exists.

```julia-repl
julia> mat=[1 2 7;4 5 6;7 8 9;10 11 19;5 7 12]
5×3 Matrix{Int64}:
  1   2   7
  4   5   6
  7   8   9
 10  11  19
  5   7  12

julia> solutionmat(mat,[95,115,182])
5-element Vector{Rational{Int64}}:
  47//4
 -17//2
  67//4
   0//1
   0//1

julia> solutionmatInt(mat,[95,115,182])
5-element Vector{Int64}:
  2285
 -5854
  4888
 -1299
     0
```
"""
function solutionmatInt(mat, v)
  if iszero(mat)
    if iszero(v) return fill(0,length(mat))
    else return false
    end
  end
  norm=TriangulizedIntegerMatTransform(mat)
  t=norm[:rowtrans]
  rs=norm[:normal][1:norm[:rank],:]
  M=vcat(rs, transpose(v))
  r=TriangulizedIntegerMatTransform(M)
  if r[:rank]==size(r[:normal],1) || r[:rowtrans][end,end]!=1
      return false
  end
  -transpose(t[1:r[:rank],:])*r[:rowtrans][end,1:r[:rank]]
end
    
"""
This  function returns  a Tuple of length  two, its  first entry  being the
result  of a call  to `solutionmatInt` with  same arguments, the second the
result of `lnullspaceInt` applied to the matrix `mat`. The calculation is
performed faster than if two separate calls would be used.
```julia_repl
julia> mat=[1 2 7;4 5 6;7 8 9;10 11 19;5 7 12]
julia> MatInt.SolutionNullspaceIntMat(mat,[95,115,182])
([2285, -5854, 4888, -1299, 0], [1 18 … 2 -6; 0 24 … 3 -7])
```
"""
function SolutionNullspaceIntMat(mat, v)
  if iszero(mat)
    len=size(mat,1)
    if iszero(v) return [fill(0,max(0,len)), Matrix{Int}(I,len,len)]
    else return [false, Matrix{Int}(I,len,len)]
    end
  end
  norm=TriangulizedIntegerMatTransform(mat)
  kern=norm[:rowtrans][norm[:rank]+1:size(mat,1),:]
  kern=baseInt(kern)
  t=norm[:rowtrans]
  rs=norm[:normal][1:norm[:rank],:]
  M=vcat(rs, transpose(v))
  r=TriangulizedIntegerMatTransform(M)
  if r[:rank]==size(r[:normal],1) || r[:rowtrans][end,end]!=1
      return [false, kern]
  end
  (-transpose(t[1:r[:rank],:])*r[:rowtrans][end,1:r[:rank]], kern)
end

function DeterminantIntMat(mat)
  sig=1
  n=size(mat,1)+2
  if n<22 return DeterminantMat(mat) end
  m=size(mat,2)+2
  if n!=m error("DeterminantIntMat: <mat> must be a square matrix") end
  A=fill(zero(eltype(mat)),m,n)
  A[2:end-1,2:end-1]=m
  A[1,1]=1
  A[n,m]=1
  r=0
  c2=1
  while m>c2
    r+=1
    c1=c2
    j=c1+1
    while j <= m
      k=r+1
      while k<=n && A[r,c1]*A[k,j]==A[k,c1]*A[r,j] k+=1 end
      if k<=n
        c2=j
        j=m
      end
      j+=1
    end
    c=mgcdex(abs(A[r,c1]), A[r+1,c1], A[r+2:n,c1])
    for i in r+2:n
      if c[i-r-1]!=0 A[r+1,:]+=A[i,:].*c[i-r-1] end
    end
    i=r+1
    while A[r,c1]*A[i,c2]==A[i,c1]*A[r,c2] i+=1 end
    if i>r+1
      c=mgcdex(abs(A[r,c1]), A[r+1,c1]+A[i,c1], [A[i,c1]])[1]+1
      A[r+1,:]+=A[i,:].* c
    end
    g=bezout(@view A[r:r+1,[c1,c2]])
    sig*=g.sign
    if sig==0 return 0 end
    A[r:r+1,:]=g.rowtrans*A[r:r+1,:]
    for i in r+2:n
      q=div(A[i,c1], A[r,c1])
      A[i,:]-=A[r,:].*q
      q=div(A[i,c2], A[r+1,c2])
      A[i,:]-=q.*A[r+1,:]
    end
  end
  for i in 2:r+1 sig*=A[i,i] end
  sig
end

function IntersectionLatticeSubspace(m)
  m*=lcm(denominator.(vcat(m...)))
  r=smith_transforms(m)
  for i in 1:length(r[:normal])
    if !iszero(r[:normal][i,:])
      r[:normal][i,:]//=maximum(abs.(r[:normal][i,:]))
    end
  end
  r[:rowtrans]^-1*r[:normal]*r[:coltrans]
end

"""
`diaconis_graham(m::Matrix{<:Integer}, moduli::Vector{<:Integer})`

returns the normal  form defined  in [Diaconis-Graham1999](biblio.htm#dg99)
for  the set of generators  defined by `m` of  the abelian group defined by
`moduli`.

`moduli`  should  have  positive  entries  such  that `moduli[i+1]` divides
`moduli[i]` for all `i`, representing the abelian group
`A=ℤ/moduli[1]×…×ℤ/moduli[n]`, where `n=length(moduli)`.

`m`  should have `n` columns, and each  line, with the `i`-th element taken
`mod  moduli[i]`, represents  an element  of `A`;  the set  of lines of `m`
should generate `A`.

The  function returns 'nothing'  if the lines  of `m` do  not generate `A`.
Otherwise it returns a named tuple `r` with fields

`r.normal`:  the Diaconis-Graham normal form, a matrix of same shape as `m`
where  either the first `n` lines are the identity matrix and the remaining
lines  are `0`,  or `length(m)=n`  and `.normal`  differs from the identity
matrix only in the entry `.normal[n,n]`, which is prime to `moduli[n]`.

`r.rowtrans`: unimodular matrix such that `r.normal==mod.(r.rowtrans*m,moduli')`

Here is an example:

```julia-repl
julia> r=diaconis_graham([3 0;4 1],[10,5])
(rowtrans = [-13 10; 4 -3], normal = [1 0; 0 2])

julia> r.normal==mod.(r.rowtrans*[3 0;4 1],[10,5]')
true
```
"""
function diaconis_graham(m::Matrix{<:Integer}, moduli::Vector{<:Integer})
  if moduli==[] return (rowtrans=[],normal=m) end
  if any(i->moduli[i]%moduli[i+1]!=0,1:length(moduli)-1)
    error("moduli[i+1] should divide moduli[i] for all i")
  end
  r=hermite_transforms(m)
  res=r.rowtrans
  m=r.normal
  n=length(moduli)
  if size(m,1)>0 && n!=size(m,2)
    error("length(moduli) should equal size(m,2)")
  end
  for i in 1:min(n,size(m,1)-1)
    l=m[i,i]
    if m[i,i] != 1
      if gcd(m[i,i], moduli[i])!=1 return end
      l1=invmod(l, moduli[i])
      e=Matrix{Int}(I,size(m,1),size(m,1))
      e[i:i+1,i:i+1]=[l1 l1-1;1-l*l1 (l+1)-l*l1]
      res=e*res
      m=mod.(e*m,moduli')
    end
  end
  r=hermite_transforms(m)
  m=mod.(r.normal, moduli')
  res=r.rowtrans*res
  if size(m,1)==n
    if m[n,n]>div(moduli[n], 2)
      m[n,n]=mod(-m[n,n], moduli[n])
      res[n,:]*=-1
    end
    l=m[n,n]
    if gcd(l,moduli[n])!=1 return end
    l1=invmod(l, moduli[n])
    for i in 1:n-1
      if m[i,n]!=0
        e=Matrix{Int}(I,size(m,1),size(m,1))
        e[[i,n],[i,n]]=[1 -m[i,n]*l1;0 1]
        res=e*res
        m=mod.(e*m, moduli')
      end
    end
  end
  (rowtrans=res, normal=m)
end

end
