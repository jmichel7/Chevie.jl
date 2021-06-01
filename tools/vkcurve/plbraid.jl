# Deals with "star" linear braids, those with associated permutation w_0
function starbraid(y, offset, B)
  n=length(y)
  if n==1 return B() end
  k=argmin(y)
  B((k:n-1)+offset...)*starbraid(y[setdiff(1:n,[k])],offset,B)/
  B((n+1-k:n-1)+offset...)
end

# In case two points have the same real projection, we use
# a "lexicographical" desingularization by "infinitesimal rotation"
function desingularized(v1, v2)
  n=length(v1)
  tan=1
  for k in 1:n
    for l in k+1:n
      rv=(real(v1[k])-real(v1[l]))
      iv=(imag(v1[k])-imag(v1[l]))
      if abs(iv*rv)>10^-8 tan=min(tan,abs(rv/iv)) end
      rv=(real(v2[k])-real(v2[l]))
      iv=(imag(v2[k])-imag(v2[l]))
      if abs(iv*rv)>10^-8 tan=min(tan,abs(rv/iv)) end
    end
  end
  [v1, v2]*Complex(1,-tan/2)
end

"""
'LBraidToWord(v1,v2,B)'

This function converts  the linear braid given by `v1`  and `v2` into an
element of the braid group `B`.

|    gap> B:=Braid(CoxeterGroupSymmetricGroup(3)); 
    function ( arg ) ... end
    gap> i:=Complex(0,1);
    I
    gap> LBraidToWord([1+i,2+i,3+i],[2+i,1+2*i,4-6*i],B);
    1|

The  list `v1` and `v2` must have the same length, say `n`. The braid group
`B`   should  be  the   braid  group  on   `n`  strings,  in  its  CHEVIE
implementation.  The elements of  `v1` (resp. `v2`)  should be `n` distinct
complex  rational  numbers.  We  use  the  Brieskorn  basepoint, namely the
contractible  set  `C+iV_ℝ`  where  `C`  is  a real chamber; therefore the
endpoints need not be equal (hence, if the path is indeed a loop, the final
endpoint must be given). The linear braid considered is the one with affine
strings  connecting each point in `v1`  to the corresponding point in `v2`.
These strings should be non-crossing. When the numbers in `v1` (resp. `v2`)
have  distinct real parts, the  real picture of the  braid defines a unique
element  of `B`. When some real parts are equal, we apply a lexicographical
desingularization,  corresponding  to  a  rotation  of  `v1` and `v2` by an
arbitrary small positive angle.
"""
# Computes, from a piecewise linear braid, the corresponding element of B_n.
# Convention:
#    we  use the  Brieskorn basepoint,  namely the  contractible set C+iV_R
#    where  C is a real chamber; therefore  the endpoints need not be equal
#    (hence,  if the  path is  indeed a  loop, the  final endpoint  must be
#    given).
#    To get rid of singular projections, a lexicographical
#    desingularization is applied.
##########################################
# Input: pair of n-tuples of complex rational numbers, ambient braid function
#         (the ambient braid group should be B_n)
# Output: the corresponding element of B_n
#########################################
# two printlevel control fields: VKCURVE.showSingularProj
#				 VKCURVE.showBraiding
#########################################

# Deals with linear braids
# 1) singular real projections are identified
# 2) calls starbraid for each
function LBraidToWord(v1, v2, B)
  n=length(v1)
  x1=real.(v1)
  y1=imag.(v1)
  x2=real.(v2)
  y2=imag.(v2)
  c=combinations(1:n,2)
  if any(x->isapprox(x1[x[1]],x1[x[2]];rtol=10^-8),c) || 
     any(x->isapprox(x2[x[1]],x2[x[2]];rtol=10^-8),c)
    if VKCURVE[:showSingularProj]
      println("WARNING: singular projection(resolved)")
    end
    return LBraidToWord(desingularized(v1, v2)..., B)
  end
  q=sortPerm(x1)
  crit=Float64[]
  for i in 1:n-1, j in i+1:n
    if x2[i^q]>x2[j^q] 
      push!(crit,(x1[i^q]-x1[j^q])/((x2[j^q]-x1[j^q]+x1[i^q])-x2[i^q]))
    end
  end
  tcrit=sort(unique(crit))
  @show tcrit
  res=B()
  u=0
  for t in tcrit
    xt=map(k->x1[k]+t*(x2[k]-x1[k]),1:n)
    yt=map(k->y1[k]+t*(y2[k]-y1[k]),1:n)
    ut=(u+t)/2
    xut=map(k->x1[k]+ut*(x2[k]-x1[k]),1:n)
    put=sortPerm(xut)
    xt=xt^put
    yt=yt^put
    xcrit=sort(unique(xt))
    for x in xcrit
      @show x,xt
      posx=findfirst(==(x),xt)
      nx=count(==(x),xt)
      res*=starbraid(yt[posx:(posx+nx)-1], posx-1, B)
    end
    u=t
  end
  if VKCURVE[:showBraiding]
   if !isempty(tcrit)
      if VKCURVE[:showInsideSegments]
        println("======================================")
        println("==    Nontrivial braiding == ",lpad(res,10),"==")
        println("======================================")
      else println("==    Nontrivial braiding == ",lpad(res,10),"==")
      end
    end
  end
  return res
end
