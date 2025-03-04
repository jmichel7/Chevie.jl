Base.:*(a::AbstractArray,b::Union{Pol,Frac,Mvp})=a.*b
Base.:*(a::Union{Pol,Frac,Mvp},b::AbstractArray)=a.*b

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

function TeXIndex(s)
  s=string(s)
  length(s)==1  ? "_"*s : "_{"*s*"}"
end

function Replace(s,p...)
# print("Replace s=$s p=$p")
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
function mkcuspidal(n,charnum,eig;no=0,qeig=0,E4=false)
  if n==28 name="F_4"
  elseif n==30 name="H_4"
  elseif n>=36 name="E_"*string(n-29)
  else name="G_{"*string(n)*"}"
  end
  if no!=0 name*="^"*string(no) end
  if E4 && eig==E(4) name*="[i]"
  elseif E4 && eig==E(4,3) name*="[-i]"
  else name*="["*xrepr(eig,TeX=true)*"]"
  end
  res=Dict(:relativeType=>TypeIrred(series=:A,indices=Int[],rank=0),
    :parameterExponents=>Int[],:charNumbers=>[charnum],
    :eigenvalue=>eig,:cuspidalName=>name)
  if n<=22 res[:levi]=1:2
  elseif n<=27 res[:levi]=1:3
  elseif n<=32 res[:levi]=1:4
  elseif n<=33 res[:levi]=1:5
  elseif n<=35 res[:levi]=1:6
  elseif n<=36 res[:levi]=1:7
  else res[:levi]=1:8
  end
  if qeig!=0 res[:qEigen]=qeig end
  res
end

"""
`CycPol(v::AbstractVector)`

This  form is an  compact way unsed  in the Chevie  library of specifying a
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
function EvalPolRoot(pol::Pol,x,n,p)
# println("pol=",pol,"\nx=",x,"\nn=",n,"\np=",p)
  if isempty(pol.c) return 0 end
  P=vcat(fill(0,mod(pol.v,n)),pol.c)
  P=map(i->Pol(P[i:n:length(P)],div(pol.v-mod(pol.v,n),n))(x*p^n),1:n)
  j=findlast(!iszero,P)
  if isnothing(j) return 0 end
  pol=Pol(P[1:j],0)
  l=(pol.v-1).+filter(i->!iszero(pol.c[i]),eachindex(pol.c))
  r=gcd(n,l...)
  pol=Pol(pol.c[1:r:length(pol.c)],div(pol.v,r))
  pol(root(x,div(n,r))*p^r)
end

"""
`VcycSchurElement(para,r[,data])``

This function computes the Schur elements for G4-22,  G25-26, G28, G32 and
imprimitive groups according to the data computed by Maria Chlouveraki.
  - `para` are vcat(the parameters of some Hecke algebra).
  - `r` is a *schur model* and data is some *schur data*.
a schur model describes the shape of the Schur element: it has the fields
  - `.factor` a (possibly fractional) *vecmonomial* coefficient
  - `.coeff` a constant coefficient
  - one or none of `.root` (a vecmonomial) or `.rootUnity`  
  - `vcyc` a vector of [vecmonomial, cyclotomic polynomial index]
  - `rootCoeff`  a constant by which multiply .root before taking root

a vecmonomial is a vector of powers for variables in `para` (plus if of length
`length(para)+1`, the power to which to raise `.root` or `.rootUnity`)

a schur data describes the Schur element in its Galois orbit : it has fields
  - `order` in which order to take the variables (para)
  - `rootPower` by which `Root1` multiply `.root`

During computation the result is evaluated as a Pol in .root in order to
use EvalPolRoot at the end.
"""
function VcycSchurElement(para,r,data=nothing)
# println("para=",para,"\nr=",r,"\ndata=",data)
  n=length(para)
  if !isnothing(data) para=para[data[:order]] else para = copy(para) end
  monomial(mon)=prod(map(^,para//1,Int.(mon[1:n])))
  if haskey(r, :rootUnity) && haskey(r,:root) error("cannot have both") end
  if haskey(r, :coeff) res = r[:coeff] else res = 1 end
  if haskey(r, :factor) res*=monomial(r[:factor]) 
     res=Pol([res],0)
  end
  function term(v)
    mon,pol=v
    if haskey(r,:rootUnity)
      tt=monomial(mon)
      if length(mon)==n+1 tt*=(r[:rootUnity]^data[:rootUnityPower])^mon[n+1] end
      Pol([cyclotomic_polynomial(pol)(tt)],0)
    elseif haskey(r, :root)
      if length(mon)==n Pol([cyclotomic_polynomial(pol)(monomial(mon))],0)
      else cyclotomic_polynomial(pol)(Pol([monomial(mon)],mon[n+1]))
      end
    else 
      Pol([cyclotomic_polynomial(pol)(monomial(mon))],0)
    end
  end
  res*=prod(term.(r[:vcyc]))
  if !haskey(r, :root) return res.c[1] end
  den=lcm(denominator.(r[:root])...)
  root=monomial(den*r[:root])
  if haskey(r, :rootCoeff) root*=r[:rootCoeff] end
  EvalPolRoot(res, root, den, data[:rootPower])
end
