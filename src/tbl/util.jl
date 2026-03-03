# util.jl  utility functions used in several data files
"""
fix illegal relativeTypes B1 and C2 which appear in HC or almost HC
series of classical groups
"""
function FixRelativeType(t)
  d=t[:relativeType]
  if d.series==:B
    if d.rank==1
      d.series=:A
      t[:charNumbers]=reverse(t[:charNumbers]) # map B1->A1
    elseif d.rank==2 && haskey(d,:cartanType) && d.cartanType==1
      d.cartanType=2;d.indices=reverse(d.indices)
      reverse!(view(t[:charNumbers],[1,5])) # map C2->B2
      if haskey(t,:parameterExponents) reverse!(t[:parameterExponents]) end
    end
  end
  InitChevie.field(d)
end

function expandrep(r,d,l)
  T=reduce(promote_type,typeof.(first.(l)))
  m=fill(zero(T),d*d*r)
  for v in l, k in v[2]
    m[Int(k)]=v[1]
  end
  Matrix.(eachslice(reshape(m,r,d,d),dims=1))
end

function compactrep(r) # inverse of expandrep
  r=stack(r,dims=1)
  v=collectby(r,eachindex(r))
  v=filter(x->!iszero(r[x[1]]),v)
  (size(r,1),size(r,2),map(x->(r[x[1]],x),v))
end

function exceptioCharName(para)
  res="\\phi_{"*string(para[1])*","*string(para[2])*"}"
  if length(para)==3 res*="'"^para[3] end
  res
end

function Replace(s,p...)
  for (src,tgt) in (p[i]=>p[i+1] for i in 1:2:length(p))
    res=empty(s)
    i=0
    while i+length(src)<=length(s)
      if @views src==s[i+1:i+length(src)]
        append!(res,tgt)
        i+=length(src)
      else
        push!(res,s[i+1])
        i+=1
      end
    end
    @views append!(res,s[i+1:end])
    s=res
  end
# println("=>",s)
  s
end

# make a cuspidal harish-chandra series record
function mkcuspidal(name,charnum,eigen;no=0,qeig=0,eig=true)
  cn=name
  if eig # if only one cuspidal in W omit display of eigenvalue
    if no!=0 cn*="^"*string(no) end
    if eigen==E(4) cn*="[i]"
    elseif eigen==E(4,3) cn*="[-i]"
    else cn*="["*xrepr(eigen,TeX=true)*"]"
    end
  end
  res=Dict(:relativeType=>TypeIrred(;series=:A,indices=Int[],rank=0),
    :parameterExponents=>Int[],:charNumbers=>[charnum],
    :eigenvalue=>eigen,:cuspidalName=>cn)
  res[:levi]=1:Dict{String,Int}( "G_4"=>2, "G_6"=>2, "G_8"=>2, "G_{14}"=>2,
 "H_3"=>3, "G_{24}"=>3, "G_{25}"=>3, "G_{26}"=>3, "G_{27}"=>3, "F_4"=>4,
 "G_{29}"=>4, "H_4"=>4, "G_{32}"=>4, "G_{33}"=>5, "G_{34}"=>6, "E_6"=>6, 
 "E_7"=>7, "E_8"=>8, "G_2"=>2, "{}^2F_4"=>4,"{}^2E_6"=>6,"{}^3G_{3,3,3}"=>3,
 "G_{3,3,3}"=>3,"{}^4G_{3,3,3}"=>3,"{}^3D_4"=>4)[name]
  if qeig!=0 res[:qEigen]=qeig end
  res
end

# adjust table prepared for [qₛ,-1]-Hecke algebra to more general [xₛ,yₛ]
function AdjustHeckeCharTable(tbl,param)
  param=improve_type(param)
  if !(tbl[:irreducibles] isa AbstractMatrix) 
    tbl[:irreducibles]=toM(tbl[:irreducibles])
  end
  for (i,w) in enumerate(tbl[:classwords])
    tbl[:irreducibles][:,i].*=prod(-last.(param[w]))
  end
  tbl
end

function braidrel(i,j,o)
  s(i,j)=[ifelse(iseven(k),j,i) for k in 1:o]
  (s(i,j),s(j,i))
end

"""
`CycPol(v::AbstractVector)`

This  form is  a compact  representation used  in the  Chevie library  of a
`CycPol`  with only  positive multiplicities:  `v` should  be a vector. The
first  element is taken as the `.coeff`  of the `CycPol`, the second as the
`.valuation`.   Subsequent  elements  are   rationals  `i//d`  representing
`(q-E(d)^i)` or are integers `d` representing `Φ_d(q)`.

```julia-repl
julia>Pol(:q);

julia> CycPol([3,-5,6,3//7])
3q⁻⁵Φ₆(q-ζ₇³)
```
"""
function CycPols.CycPol(v::AbstractVector)
  coeff=v[1]
  valuation=convert(Int,v[2])
  vv=Pair{Rational{Int},Int}[]
  v1=convert.(Rational{Int},v[3:end])
  for i in v1
    if denominator(i)==1
      k=convert(Int,i)
      for j in prime_residues(k) push!(vv,j//k=>1) end
    else
      push!(vv,i=>1)
    end
  end
  CycPol(coeff,valuation,ModuleElt(vv))
end

"""
 EvalPolRoot(pol, x, n, p) compute pol(p*x^(1/n))
  
  The point of this routine is to avoid unnecessary root extractions
  during evaluation (e.g., if pol has no terms of odd degree and n=2,
  then no root extraction is necessary).
"""
function EvalPolRoot(p::Pol,x,n,f)
  if iszero(p) return 0 end
  v=valuation(p)
  d=degree(p)
  r=n
  for i in v:d
   if !iszero(p[i]) r=gcd(r,i) end
  end
  Pol(p[v:r:d],div(v,r))(root(x,div(n,r))*f^r)
end

"""
`VcycSchurElement(para,r[,data])``

This function computes the Schur elements for G4-22,  G25-26, G28, G32 and
imprimitive groups according to the data computed by Maria Chlouveraki.
  - `para` is `vcat(H.para...)` for some Hecke algebra `H`.
  - `r` is a *schur model* and data is a *schur data*.
a schur model describes the shape of the Schur element: it has the fields
  - `.factor` a (possibly fractional) *vecmonomial* coefficient
  - `.coeff` a constant coefficient
  - one or none of `.root` (a vecmonomial) or `.rootUnity`  
  - `vcyc` a vector of [vecmonomial, cyclotomic polynomial index]
  - `rootCoeff`  a constant by which multiply .root before taking root

a vecmonomial is a vector of powers to apply to `para` (plus if of length
`length(para)+1`, the power to which to raise `.root` or `.rootUnity`)

a schur data describes the Schur element in its Galois orbit : it has fields
  - `order` in which order to take the variables (para)
  - `rootPower` by which `Root1` multiply `.root`

During computation the result is evaluated as a Pol in .root in order to
use EvalPolRoot at the end.
"""
function VcycSchurElement(para::Vector,r,data=nothing)
# println("para=",para,"\nr=",r,"\ndata=",data)
  if !isnothing(data) para=para[data[:order]] else para=copy(para) end
  monomial(m,p)=prod(i->(p[i]//1)^(Int(m[i])),eachindex(p))
  if haskey(r,:rootUnity) && haskey(r,:root) error("cannot have both") end
  res=haskey(r,:coeff) ? r[:coeff] : 1
  if haskey(r,:factor) res=Pol([res*monomial(r[:factor],para)]) end
  res*=prod(r[:vcyc])do (mon,pol)
    C=cyclotomic_polynomial(pol)
    n=length(para)
    if haskey(r,:rootUnity)
      tt=monomial(mon,para)
      if length(mon)==n+1 tt*=(r[:rootUnity]^data[:rootUnityPower])^mon[n+1] end
      Pol([C(tt)])
    elseif haskey(r,:root)
      if length(mon)==n Pol([C(monomial(mon,para))])
      else C(Pol([monomial(mon,para)],mon[n+1]))
      end
    else Pol([C(monomial(mon,para))])
    end
  end
  if !haskey(r,:root) return res[begin] end
  den=lcm(denominator.(r[:root]))
  root=monomial(den*r[:root],para)
  if haskey(r,:rootCoeff) root*=r[:rootCoeff] end
  EvalPolRoot(res,root,den,data[:rootPower])
end
