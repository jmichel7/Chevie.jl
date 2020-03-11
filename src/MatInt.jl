module MatInt
using ..Gapjm

export ComplementIntMat, NullspaceIntMat

IdentityMat(n)=map(i->one(rand(Int,n,n))[i,:],1:n)
NullMat(i,j)=[zeros(Int,j) for k in 1:i]
MatMul(a,b)=[[sum(j->a[i][j]*b[j][k],eachindex(a[1])) 
             for k in eachindex(b[1])] for i in eachindex(a)]

# MATINTsplit(N,a) largest factor of N prime to a
function MATINTsplit(N, a)
  while a != 1
      a = gcd(a, N)
      N = div(N, a)
  end
  N
end

# MATINTrgcd(N,a) smallest nonnegative c such that gcd(N,a+c) = 1
function MATINTrgcd(N, a)
  if N == 1 return 0 end
  r = [mod(a-1, N)]
  d = [N]
  c = 0
  while true
    for i in 1:length(r) r[i] = mod(r[i] + 1, d[i]) end
    i = findfirst(<=(0),r)
    if isnothing(i)
      g=1
      i=0
      while g==1 && i<length(r)
        i+=1
        g=gcd(r[i], d[i])
      end
      if g == 1 return c end
      q = MATINTsplit(div(d[i], g), g)
      if q>1
       push!(r,mod(r[i], q))
       push!(d, q)
      end
      r[i] = 0
      d[i] = g
    end
    c+=1
  end
end

mutable struct GcdRes
  gcd::Int
  coeff1::Int
  coeff2::Int
  coeff3::Int
  coeff4::Int
end
#############################################################################
##
#F  Gcdex( <m>, <n> ) . . . . . . . . . . greatest common divisor of integers
##
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
  if n==0
    GcdRes(f, fm, 0, gm, 1)
  else
    GcdRes(f, fm, div(f-fm*m,n), gm, div(-gm*m,n))
  end
end

#  MATINTmgcdex(<N>,<a>,<v>) - Returns c[1],c[2],...c[k] such that
#   gcd(N,a+c[1]*v[1]+...+c[n]*v[k]) = gcd(N,a,v[1],v[2],...,v[k])
function MATINTmgcdex(N, a, v)
  l = length(v)
  M = zeros(Int,l)
  h = N
  for i = 1:l
    g = h
    h = gcd(g, v[i])
    M[i] = div(g, h)
  end
  h=gcd(a,h)
  g=div(a,h)
  c = zeros(Int,l)
  for i in l:-1:1
    b = div(v[i], h)
    d = MATINTsplit(M[i], b)
    if d==1 c[i] = 0
    else c[i] = MATINTrgcd(d, mod(div(g,b), d))
      g+=c[i]*b
    end
  end
  c
end

#  MATINTbezout(a,b,c,d) - returns P to transform A to hnf (PA=H);
#
#  [st][ab]=[ef] 
#  [uv][cd] [0g]
#
function MATINTbezout(a, b, c, d)
  e=Gcdex(a, c)
  f=e.coeff1*b+e.coeff2*d
  g=e.coeff3*b+e.coeff4*d
  if g < 0
    e.coeff3 = -e.coeff3
    e.coeff4 = -e.coeff4
    g = -g
  end
  if g > 0
    q = div(f-mod(f, g), g)
    e.coeff1-=q*e.coeff3
    e.coeff2-=q*e.coeff4
  end
  e
end

## SNFofREF - fast SNF of REF matrix
function SNFofREF(R, destroy)
  n = length(R)
  m = length(R[1])
  piv = map(x->findfirst(!iszero,x), R)
  r = findfirst(==(false),piv)
  if r == false r = length(piv)
  else
      r-=1
      piv=piv[1:r]
  end
  append!(piv, Difference(1:m, piv))
  if destroy
      T = R
      for i in 1:r T[i][1:m] = T[i][piv] end
  else
      T = NullMat(n, m)
      for j in 1:m
          for i in 1:min(r, j) T[i][j] = R[i][piv[j]] end
      end
  end
  si = 1
  A = []
  d = 2
  for k = 1:m
      if k <= r
          d*=abs(T[k][k])
          Apply(T[k], x->mod(x, 2d))
      end
      t = min(k, r)
      for i = t-1:-1:si
          t = MATINTmgcdex(A[i], T[i][k], [T[i+1][k]])[1]
          if t != 0
              T[i]+=T[i+1]*t
              Apply(T[i], (x->begin mod(x, A[i]) end))
          end
      end
      for i = si:min(k - 1, r)
          g = Gcdex(A[i], T[i][k])
          T[i][k] = 0
          if g.gcd != A[i]
              b = div(A[i], g.gcd)
              A[i] = g.gcd
              for ii = i + 1:min(k - 1, r)
                  T[ii]+=map(x->mod(x,A[ii]),T[i]*-g.coeff2*div(T[ii][k],A[i]))
                  T[ii][k]= b*T[ii][k]
                  Apply(T[ii],x->mod(x,A[ii]))
              end
              if k <= r
                  t = g.coeff2 * div(T[k][k], g.gcd)
                  T[k]+=-T[i]-t
                  T[k][k]*=b
              end
              Apply(T[i], x->mod(x, A[i]))
              if A[i]==1 si=i+1 end
          end
      end
      if k <= r
          A[k] = abs(T[k][k])
          Apply(T[k], x->mod(x, A[k]))
      end
  end
  for i = 1:r
      T[i][i] = A[i]
  end
  return T
end

BITLISTS_NFIM = [
  [false, false, false, false, false], 
  [true, false, false, false, false], 
  [false, true, false, false, false], 
  [true, true, false, false, false], 
  [false, false, true, false, false], 
  [true, false, true, false, false], 
  [false, true, true, false, false], 
  [true, true, true, false, false], 
  [false, false, false, true, false], 
  [true, false, false, true, false], 
  [false, true, false, true, false], 
  [true, true, false, true, false], 
  [false, false, true, true, false], 
  [true, false, true, true, false], 
  [false, true, true, true, false], 
  [true, true, true, true, false], 
  [false, false, false, false, true], 
  [true, false, false, false, true], 
  [false, true, false, false, true], 
  [true, true, false, false, true], 
  [false, false, true, false, true], 
  [true, false, true, false, true], 
  [false, true, true, false, true], 
  [true, true, true, false, true], 
  [false, false, false, true, true], 
  [true, false, false, true, true], 
  [false, true, false, true, true], 
  [true, true, false, true, true], 
  [false, false, true, true, true], 
  [true, false, true, true, true], 
  [false, true, true, true, true], 
  [true, true, true, true, true]]
###########################################################
#
# DoNFIM(<mat>,<options>)
#
# Options bit values:
#
# 1  - Triangular / Smith
# 2  - No / Yes  Reduce off diag entries
# 4  - No / Yes  Row Transforms 
# 8  - No / Yes  Col Transforms
# 16 - change original matrix in place (The rows still change) -- save memory
#
# Compute a Triangular, Hermite or Smith form of the n x m 
# integer input matrix A.  Optionally, compute n x n / m x m
# unimodular transforming matrices which satisfy Q C A = H 
# or  Q C A B P = S.
#
# Triangular / Hermite :
#
# Let I be the min(r+1,n) x min(r+1,n) identity matrix with r = rank(A).
# Then Q and C can be written using a block decomposition as
#
#             [ Q1 |   ]  [ C1 | C2 ]
#             [----+---]  [----+----]  A  =  H.
#             [ Q2 | I ]  [    | I  ]
#
# Smith :
#
#  [ Q1 |   ]  [ C1 | C2 ]     [ B1 |   ]  [ P1 | P2 ]
#  [----+---]  [----+----]  A  [----+---]  [----+----] = S.
#  [ Q2 | I ]  [    | I  ]     [ B2 | I ]  [ *  | I  ]
#
# * - possible non-zero entry in upper right corner...
#				
#
function DoNFIM(mat::AbstractVector, opt::Int)
  if isempty(mat) mat=[Int[]] end
# println("mat=$mat opt=$opt")
  opt = BITLISTS_NFIM[opt + 1]
  TRIANG = opt[1]
  REDDIAG = opt[2]
  ROWTRANS = opt[3]
  COLTRANS = opt[4]
  INPLACE = opt[5]
  sig = 1
  n=length(mat)+2
  m=length(mat[1])+2
  k=fill(0, max(0,m))
  if INPLACE
    A=mat
    for i in n-1:-1:2
      A[i]=copy(k)
      A[i][2:m-1]=A[i-1]
    end
  else
    A=[copy(k) for i in 1:n]
    for i in 2:n-1 A[i][2:m-1].=mat[i-1] end
  end
  A[1]=copy(k)
  A[n]=k
  A[1][1]=1
  A[n][m]=1
  if ROWTRANS
    C=IdentityMat(n)
    Q=NullMat(n, n)
    Q[1][1] = 1
  end
  if TRIANG && COLTRANS
      B = IdentityMat(m)
      P = IdentityMat(m)
  end
  r = 0
  c2 = 1
  rp=Int[]
  while m > c2
    c1=c2
    push!(rp,c1)
    r=r+1
    if ROWTRANS Q[r+1][r+1]=1 end
    j=c1+1
    while j <= m
        k=r+1
        while k <= n && A[r][c1]*A[k][j]==A[k][c1]*A[r][j] k=k+1 end
        if k <= n
            c2 = j
            j = m
        end
        j=j+1
    end
    if TRIANG && ((COLTRANS || ROWTRANS) && c2 < m)
        N = gcd(map(x->x[c2],A[r:n]))
        L = vcat(c1+1:c2-1, c2+1:m-1)
        push!(L, c2)
        for j = L
            if j == c2
                b = A[r][c2]
                a = A[r][c1]
                for i = r+1:n
                    if b != 1
                        g = Gcdex(b, A[i][c2])
                        b = g.gcd
                        a = g.coeff1 * a + g.coeff2 * A[i][c1]
                    end
                end
                N = 0
                for i = r:n
                  if N != 1 N = gcd(N, A[i][c1] - div(A[i][c2], b) * a) end
                end
            else
              c=MATINTmgcdex(N, A[r][j], map(x->x[j],A[r+1:n]))
              b=A[r][j] + c*map(x->x[j],A[r+1:n])
              a=A[r][c1] + c*map(x->x[c1],A[r+1:n])
            end
            t = MATINTmgcdex(N, a, [b])[1]
            tmp = A[r][c1] + t * A[r][j]
            if tmp==0 || tmp*A[k][c2]==(A[k][c1]+t*A[k][j])*A[r][c2]
                t+=1+MATINTmgcdex(N, a+t*b+b,[b])[1]
            end
            if t > 0
                for i in 1:n A[i][c1]+=t*A[i][j] end
                if COLTRANS B[j][c1]+=t end
            end
        end
        if A[r][c1]*A[k][c1+1]==A[k][c1]*A[r][c1+1]
            for i in 1:n A[i][c1+1]+=A[i][c2] end
            if COLTRANS B[c2][c1+1]=1 end
        end
        c2 = c1 + 1
    end
    c = MATINTmgcdex(abs(A[r][c1]), A[r+1][c1], map(i->A[i][c1],r+2:n))
    for i in r+2:n
      if c[i-r-1]!=0
        A[r+1]+=A[i]*c[i-r-1]
        if ROWTRANS
          C[r+1][i]=c[i-r-1]
          Q[r+1]+=Q[i]*c[i-r-1]
        end
      end
    end
    i=r+1
    while A[r][c1]*A[i][c2]==A[i][c1]*A[r][c2] i+=1 end
    if i>r+1
      c=MATINTmgcdex(abs(A[r][c1]), A[r+1][c1]+A[i][c1], [A[i][c1]])[1]+1
      A[r+1]+=A[i]*c
      if ROWTRANS
        C[r+1][i]+=c
        Q[r+1]+=Q[i]*c
      end
    end
    g = MATINTbezout(A[r][c1], A[r][c2], A[r+1][c1], A[r+1][c2])
    sig*= sign(A[r][c1]*A[r+1][c2]-A[r][c2]*A[r+1][c1])
    A[[r,r+1]]=MatMul([[g.coeff1,g.coeff2],[g.coeff3, g.coeff4]],A[[r,r+1]])
    if ROWTRANS
      Q[[r,r+1]]=MatMul([[g.coeff1,g.coeff2],[g.coeff3,g.coeff4]],Q[[r,r+1]])
    end
    for i in r+2:n
      q=div(A[i][c1], A[r][c1])
      A[i]-=A[r]*q
      if ROWTRANS Q[i]-=Q[r]*q end
      q=div(A[i][c2], A[r+1][c2])
      A[i]-=A[r+1]*q
      if ROWTRANS Q[i]-=Q[r+1]*q end
    end
  end
  push!(rp,m) # length(rp)==r+1
  if n==m && r+1<n sig=0 end
  if TRIANG && !(ROWTRANS || COLTRANS)
    for i in 2:n-1 A[i-1]=A[i][2:m-1] end
    A=A[1:n-2]
    R = Dict(:normal => SNFofREF(A, INPLACE), :rank => r - 1)
    if n==m R[:signdet] = sig end
    return R
  end
  if !TRIANG && REDDIAG || TRIANG && COLTRANS
    for i in r:-1:1
      for j in i+1:r+1
        q = div(A[i][rp[j]]-mod(A[i][rp[j]], A[j][rp[j]]), A[j][rp[j]])
        A[i]-=A[j]*q
        if ROWTRANS Q[i]-=Q[j]*q end
      end
      if TRIANG && i<r
        for j in i+1:m
          q = div(A[i][j], A[i][i])
          for k in 1:i A[k][j]-=q*A[k][i] end
          if COLTRANS P[i][j] = -q end
        end
      end
    end
  end
  if TRIANG && (ROWTRANS && !COLTRANS)
    for i in 1:r-1
      t=A[i][i]
      A[i]=map(x->0, 1:m)
      A[i][i] = t
    end
    for j in r+1:m-1
      A[r][r] = gcd(A[r][r], A[r][j])
      A[r][j] = 0
    end
  end
  if TRIANG && COLTRANS && r<m-1
    c=MATINTmgcdex(A[r][r], A[r][r+1], A[r][r+2:m-1])
    for j in r+2:m-1
        A[r][r+1]+=c[j-r-1] * A[r][j]
        B[j][r+1]=c[j-r-1]
        for i in 1:r P[i][r+1]+=c[j-r-1] * P[i][j] end
    end
    P[r+1]=zeros(Int,m)
    P[r+1][r+1]=1
    g = Gcdex(A[r][r], A[r][r + 1])
    A[r][r] = g.gcd
    A[r][r+1] = 0
    for i in 1:r+1
      t = P[i][r]
      P[i][r] = P[i][r] * g.coeff1 + P[i][r+1] * g.coeff2
      P[i][r+1] = t * g.coeff3 + P[i][r+1] * g.coeff4
    end
    for j in r+2:m-1
      q = div(A[r][j], A[r][r])
      for i in 1:r+1 P[i][j]-=q*P[i][r] end
      A[r][j]=0
    end
    for i in r+2:m-1
      P[i] = map(x->0, 1:m)
      P[i][i] = 1
    end
  end
  if ROWTRANS for i in r+2:n Q[i][i]=1 end end
  for i in 2:n-1 A[i-1]=A[i][2:m-1] end
  A=A[1:n-2]
  R = Dict{Symbol,Any}(:normal => A)
  if ROWTRANS
    R[:rowC] = map(x->x[2:n-1],C[2:n-1])
    R[:rowQ] = map(x->x[2:n-1],Q[2:n-1])
  end
  if TRIANG && COLTRANS
    R[:colC] = map(x->x[2:m-1],B[2:m-1])
    R[:colQ] = map(x->x[2:m-1],P[2:m-1])
  end
  R[:rank]=r-1
  if n==m R[:signdet]=sig end
  return R
end

function NormalFormIntMat(mat, options)
  r=DoNFIM(mat, options)
  opt=BITLISTS_NFIM[options+1]
  if opt[3] r[:rowtrans]=MatMul(r[:rowQ],r[:rowC]) end
  if opt[1] && opt[4] r[:coltrans]=MatMul(r[:colC],r[:colQ]) end
  r
end

TriangulizedIntegerMat(mat)=DoNFIM(mat,0)[:normal]
TriangulizedIntegerMatTransform(mat)=NormalFormIntMat(mat,4)
TriangulizeIntegerMat(mat)=DoNFIM(mat,16)
HermiteNormalFormIntegerMat(mat)=DoNFIM(mat,2)[:normal]
HermiteNormalFormIntegerMatTransform(mat)=NormalFormIntMat(mat, 6)
SmithNormalFormIntegerMat(mat)=DoNFIM(mat,1)[:normal]
SmithNormalFormIntegerMatTransforms(mat)=NormalFormIntMat(mat,13)
DiagonalizeIntMat(mat)=DoNFIM(mat,17)
function BaseIntMat(mat)
  norm=NormalFormIntMat(mat, 2)
  norm[:normal][1:norm[:rank]]
end
function BaseIntersectionIntMats(M1, M2)
  M = vcat(M1, M2)
  r = NormalFormIntMat(M,4)
  T = map(x->x[1:length(M1)],r[:rowtrans][r[:rank]+1:length(M)])
  if !isempty(T) T=MatMul(T,M1) end
  BaseIntMat(T)
end

function ComplementIntMat(full, sub)
  F = BaseIntMat(full)
  if length(sub) == 0 || iszero(sub)
    return Dict(:complement => F, :sub => Vector{Int}[], :moduli => Int[])
  end
  S = BaseIntersectionIntMats(F, sub)
  if S!=BaseIntMat(sub) error("sub must be submodule of full") end
  M = vcat(F, S)
  T = toL(Int.(inv(Rational.(toM(NormalFormIntMat(M,4)[:rowtrans])))))
  T = map(x->x[1:length(F)],T[length(F)+1:length(T)])
  r = NormalFormIntMat(T, 13)
  M = toL(Int.(inv(Rational.(toM(r[:coltrans])))*toM(F)))
  Dict(:complement=>BaseIntMat(M[1+r[:rank]:length(M)]),
      :sub=>MatMul(MatMul(r[:rowtrans],T),F), 
      :moduli=>map(i->r[:normal][i][i],1:r[:rank]))
end

function NullspaceIntMat(mat)
  norm=NormalFormIntMat(mat, 4)
  kern=norm[:rowtrans][norm[:rank]+1:length(mat)]
  return BaseIntMat(kern)
end

function SolutionIntMat(mat, v)
  if iszero(mat)
      if iszero(v) return fill(0, max(0, (1 + length(mat)) - 1))
      else return false
      end
  end
  norm = NormalFormIntMat(mat, 4)
  t = norm[:rowtrans]
  rs = (norm[:normal])[1:norm[:rank]]
  M = vcat(rs, [v])
  r = NormalFormIntMat(M, 4)
  if r[:rank]==length(r[:normal]) || r[:rowtrans][length(M)][length(M)]!=1
      return false
  end
  return -r[:rowtrans][length(M)][1:r[:rank]]*t[1:r[:rank]]
end
    
function SolutionNullspaceIntMat(mat, v)
  if iszero(mat)
    len = length(mat)
    if iszero(v) return [fill(0, max(0, (1 + len) - 1)), IdentityMat(len)]
    else return [false, IdentityMat(len)]
    end
  end
  norm = NormalFormIntMat(mat, 4)
  kern = (norm[:rowtrans])[norm[:rank] + 1:length(mat)]
  kern = BaseIntMat(kern)
  t = norm[:rowtrans]
  rs = (norm[:normal])[1:norm[:rank]]
  M = Concatenation(rs, [v])
  r = NormalFormIntMat(M, 4)
  if r[:rank] == length(r[:normal]) || r[:rowtrans][length(M)][length(M)]!=1
      return [false, kern]
  end
  return [-r[:rowtrans][length(M)][1:r[:rank]]*t[1:r[:rank]], kern]
end

function DeterminantIntMat(mat)
  sig = 1
  n = length(mat) + 2
  if n < 22 return DeterminantMat(mat) end
  m = length(mat[1]) + 2
  if !n == m error("DeterminantIntMat: <mat> must be a square matrix") end
  A = [map(x->0, 1:m)]
  for i = 2:n - 1
      A[i] = [0]
      A[i] = Append(A[i], mat[i - 1])
      A[i][m] = 0
  end
  A[n] = map(x->0, 1:m)
  A[1][1] = 1
  A[n][m] = 1
  r = 0
  c2 = 1
  while m > c2
    r = r + 1
    c1 = c2
    j = c1 + 1
    while j <= m
        k = r + 1
        while k <= n && (A[r])[c1] * (A[k])[j] == (A[k])[c1] * (A[r])[j]
            k = k + 1
        end
        if k <= n
            c2 = j
            j = m
        end
        j = j + 1
    end
    c = MATINTmgcdex(abs((A[r])[c1]), (A[r + 1])[c1], (A[r + 2:n])[c1])
    for i in r + 2:n
        if c[i-r-1]!=0 A[r+1]+=A[i]*c[i-r-1] end
    end
    i = r + 1
    while A[r][c1]*A[i][c2]==A[i][c1]*A[r][c2] i=i+1 end
    if i>r+1
      c=MATINTmgcdex(abs(A[r][c1]), A[r+1][c1]+A[i][c1], [A[i][c1]])[1]+1
      A[r+1]+=A[i] * c
    end
    g = MATINTbezout(A[r][c1], A[r][c2], A[r + 1][c1], A[r + 1][c2])
    sig*=sign(A[r][c1] * A[r + 1][c2] - A[r][c2] * A[r + 1][c1])
    if sig == 0 return 0 end
    A[[r,r+1]]=[[g[:coeff1], g[:coeff2]], [g[:coeff3], g[:coeff4]]]*A[[r,r+1]]
    for i in r+2:n
      q = div(A[i][c1], A[r][c1])
      A[i]-=A[r]*q
      q = div(A[i][c2], A[r+1][c2])
      A[i]-=q*A[r+1]
    end
  end
  for i in 2:r+1 sig*=A[i][i] end
  return sig
end

function IntersectionLatticeSubspace(m)
  m = m * Lcm(map(denominator, Concatenation(m)))
  r = SmithNormalFormIntegerMatTransforms(m)
  for i = 1:length(r[:normal])
    if (r[:normal])[i] != 0 * (r[:normal])[i]
      r[:normal][i]=r[:normal][i] // maximum(map(abs, r[:normal][i]))
    end
  end
  return r[:rowtrans] ^ -1 * r[:normal] * r[:coltrans]
end
#############################################################################
##
#F  DiaconisGraham(m,moduli)
#
# m  should be  an integral  matrix with  n columns where n=Length(moduli),
# such  that each  line (with  the i-th  element taken  mod moduli[i]) of m
# represents an element of the group A=Z/moduli[1] x ... x Z/moduli[n], and
# such  that the set of  lines of m generates  A. The moduli should be such
# that moduli[i+1] divides moduli[i] for all i.
#
# The function returns  a record  r with fields
# r.normal          the Diaconis-Graham normal form, a matrix of same shape
#    as  m where either the  first n lines are  the identity matrix and the
#    remaining  lines are  0, or  Length(m)=n and  .normal differs from the
#    identity  matrix only  in the  entry .normal[n][n],  which is prime to
#    moduli[n]
# r.rowtrans        a unimodular matrix such that  
#   r.normal=List(r.rowtrans*m,v->Zip(v,moduli,function(x,y)return x mod y;end))
#
function DiaconisGraham(m, moduli)
  if moduli == [] return Dict{Symbol, Any}(:rowtrans => [], :normal => m) end
  if any(i->mod(moduli[i],moduli[i+1])!=0,1:length(moduli)-1)
    error("DiaconisGraham(m,moduli): moduli[i+1] should divide moduli[i] for all i")
  end
  r = HermiteNormalFormIntegerMatTransform(m)
  res = r[:rowtrans]
  m = r[:normal]
  n = length(moduli)
  if length(m) > 0 && n != length(m[1])
    error("DiaconisGraham(m,moduli): moduli  &&  m[1] should have same length")
  end
  ZipMod = function (m, moduli)
    return map(v-> map(function(x, y) return mod(x, y) end, v, moduli), m)
  end
  for i = 1:min(n,length(m)-1)
      l = m[i][i]
      if m[i][i] != 1
          if gcd(m[i][i], moduli[i]) != 1 return false end
          l1 = invmod(l, moduli[i])
          e = IdentityMat(length(m))
          e[[i, i + 1]][[i, i+1]] = [[l1, l1-1], [1-l*l1, (l+1)-l*l1]]
          res = e*res
          m = ZipMod(e * m, moduli)
      end
  end
  r = HermiteNormalFormIntegerMatTransform(m)
  m = ZipMod(r[:normal], moduli)
  res = MatMul(r[:rowtrans],res)
  if length(m) == n
    if m[n][n]>div(moduli[n], 2)
      m[n][n]=mod(-m[n][n], moduli[n])
      res[n]*=-1
    end
    l=m[n][n]
    if gcd(l,moduli[n])!=1 return false end
    l1 = invmod(l, moduli[n])
    for i in 1:n-1
      if m[i][n]!=0
        e=IdentityMat(length(m))
        e[i][[i,n]] = [1,-m[i][n]*l1]
        e[n][[i,n]] = [0, 1]
        res = MatMul(e,res)
        m = ZipMod(MatMul(e,m), moduli)
      end
    end
  end
  return Dict{Symbol, Any}(:rowtrans => res, :normal => m)
end

end
