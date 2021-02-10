# Hand-translated part of chevie/tbl/exceptio.jl
# (C) Jean Michel 1999-2017
# data common to several (but not all) types of reflection groups

# an addition
chevieset(["A","B","D"],:EigenvaluesGeneratingReflections,(t->r->fill(1//2,r)))

chevieset(["G25","G26","G29","G31","G32","G34"],:CartanMat,
  function(t)
    r=chevieget(t,:GeneratingRoots)
    eig=map(x->E(;r=x),chevieget(t,:EigenvaluesGeneratingReflections))
    toL(toM(coroot.(r,eig))*permutedims(toM(r)))
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

"""
 EvalPolRoot(pol, x, n, p) compute pol(p*x^(1/n))
  
  The point of this routine is to avoid unnecessary root extractions
  during evaluation (e.g., if pol has no terms of odd degree and n=2,
  then no root extraction is necessary).

 this was in lib/util.g but is used only here
"""
function EvalPolRoot(pol::Pol,x,n,p)
# println("pol=",pol,"\nx=",x,"\nn=",n,"\np=",p)
  if isempty(pol.c) return 0 end
  P=vcat(fill(0,mod(pol.v,n)),pol.c)
  P=map(i->Pol(P[i:n:length(P)],div(pol.v-mod(pol.v,n),n))(x*p^n),1:n)
  j=findlast(!iszero,P)
  if isnothing(j) return 0 end
  pol=Pol(P[1:j],0)
  l=pol.v-1+filter(i->!iszero(pol.c[i]),eachindex(pol.c))
  r=gcd(n,l...)
  pol=Pol(pol.c[1:r:length(pol.c)],div(pol.v,r))
  pol(GetRoot(x,div(n,r))*p^r)
end

#  VcycSchurElement(para,r(schur model)[,data(schur data)])
#
#  This function computes the Schur elements for G4-22,  G25-26, G28, G32
#  according to the data computed by M. Chlouveraki.
#  para is the list of parameters of the algebra.
#  schur model describes the shape of the Schur element: it has the fields
#   .factor=(possibly fractional) vecmonomial
#   .coeff= a constant
#   [nothing] or [.root=vecmonomial] or [.rootUnity]  
#   vcyc= a list of pairs [vecmonomial, cyclotomic polynomial index]
#   rootCoeff=  a constant by which multiply .root before taking root
#  vecmonomial=vector of powers for elts of para (plus possibly
#     the power to which to raise root or rootUnity)
#  schur data describes the Schur element in its Galois orbit : it has fields
#   order: in which order to take the variables
#   rootPower: by which E(root)^i multiply .root
function VcycSchurElement(para,r,data=nothing)
# println("para=",para,"\nr=",r,"\ndata=",data)
  n=length(para)
  if !isnothing(data) para=para[data[:order]] else para = copy(para) end
  monomial(mon)=prod(map(^,para,Int.(mon[1:n])))
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
     if length(mon)==n return Pol([cyclotomic_polynomial(pol)(monomial(mon))],0)
     else return cyclotomic_polynomial(pol)(Pol([monomial(mon)],mon[n+1]))
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

# symbols for types B,C,D,2D: 2-line symbols of defect d and rank n
function BDSymbols(n,d)
  n-=div(d^2,4)
  if n<0 return Vector{Vector{Int}}[] end
  if d>0 return map(x->symbol_partition_tuple(x,d),partition_tuples(n,2)) end
   return map(chevieget(:D,:symbolcharparam),
              chevieget(:imp,:CharInfo)(2,2,n)[:charparams])
end
