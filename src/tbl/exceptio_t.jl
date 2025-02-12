# Hand-translated part of chevie/tbl/exceptio.jl
# (C) Jean Michel 1999-2017
# data common to several (but not all) types of reflection groups

# an addition
chevieset([:A,:B,:D],:EigenvaluesGeneratingReflections,t->r->fill(1//2,r))

chevieset([:G25,:G26,:G29,:G31,:G32,:G34],:CartanMat,
  function(t)
    r=chevieget(t,:GeneratingRoots)
    eig=map(x->Root1(;r=x),chevieget(t,:EigenvaluesGeneratingReflections))
    toL(toM(coroot.(r,eig))*transpose(toM(r)))
  end)

chevieset([:E7, :E8, :H3, :H4], :Invariants, t->
  function()
    C=toM(chevieget(t, :CartanMat))
    r=toM(roots(C))*C
    map(d->function(arg...)sum(a->sum(arg.*a)^d,eachrow(r)) end, 
        chevieget(t, :ReflectionDegrees))
  end
)

chevieset([:G24,:G27,:G29,:G33,:G34,:E6,:E7,:E8,:H3,:H4], 
  :FactorizedSchurElement, t->function(ch,para,arg...)
   c=chevieget(t,:CycPolSchurElements)[findfirst(==(ch),chevieget(t,:CharInfo)()[:charparams])]
   q=-para[1][1]//para[1][2]
   res=HeckeAlgebras.FactSchur(Mvp(c[1]*q^Int(c[2])), 
                 map(v->(pol=CycPol([1,0,v]),monomial=q),c[3:length(c)]))
   HeckeAlgebras.simplify(res)
 end
)

chevieset(["G24","G25","G26","G27","G29","G31","G32","G33","G34","H3","H4","2E6","2F4","3D4","E6","E7","E8","F4","G2"],:IrredInfo,function(t)
  ci=chevieget(t,:CharInfo)()
  map((x,y)->Dict{Symbol,Any}(:charparam=>x,:charname=>y),ci[:charparams],ci[:charnames])
end)

chevieset(["G24", "G25", "G27", "G29", "G31", "G32", "G33", "G34", "E6", "E7", "E8", "2E6", "2F4", "3D4", "H3", "H4"], :ReflectionName, t->
function(option)
  i=["G24","G25","G27","G29","G31","G32","G33","G34","E6","E7","E8","2E6","2F4","3D4","H3","H4"]
  o=["G_{24}","G_{25}","G_{27}","G_{29}","G_{31}","G_{32}","G_{33}","G_{34}","E_6","E_7","E_8","{}^2E_6","{}^2F_4","{}^3D_4","H_3","H_4"]
  haskey(option,:TeX) ? o[findfirst(==(t),i)] : t
end)

chevieset(["D", "2A", "2D"],:ReflectionName,t->
  function (r, option)
    i=["D", "2A", "2D"]
    o=["D", "{}^2A", "{}^2D"]
    if haskey(option, :arg) string(t, ",", r)
    elseif haskey(option,:TeX)
     string(o[findfirst(==(t),i)],"_",r<10 ? string(r) : string("{",r,"}"))
    else string(t, r)
    end
end)

chevieset(["3D4", "E6", "2E6", "E7", "E8", "F4", "2F4", "G2", "H3", "H4"], :CharTable, t-> function ()
  rank=t[end]-'0'
  res=chevieget(t, :HeckeCharTable)(map(x->[1,-1],1:rank),map(x->1,1:rank))
  CHEVIE[:compat][:ChangeIdentifier](res,string("W(",t,")"))
  res
end)

chevieset([:G24,:G27,:G29,:G33,:G34,:H3,:H4,:E6,:E7,:E8],:PoincarePolynomial,t->
function (q)
  prod(x->sum(y->(-q[1][1]//q[1][2])^y,0:x-1),chevieget(t,:ReflectionDegrees))
end)

function ExpandRep(r, d, l) # decompress representation of r gens of dim d
  T=reduce(promote_type,typeof.(first.(l)))
  m=map(i->map(j->fill(zero(T),d),1:d), 1:r)
  for v in l
    for k in @view v[2:end]
      q,i=divrem(Int(k),d^2)
      q1,r1=divrem(i,d)
      m[q+1][q1+1][r1+1]=v[1]
    end
  end
  return m
end
export ExpandRep

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
export VcycSchurElement
