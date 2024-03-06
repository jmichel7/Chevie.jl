function findfractionalpowers(W)
  if isempty(gens(W)) return [0] end
  uc=UnipotentCharacters(W)
  uc.fractions=Union{Rational{Int},Nothing}[]
  append!(uc.fractions,fill(nothing,length(uc)))
  reasons=map(x->[], 1:length(uc))
  for h in uc.harishChandra
    if length(h[:levi])!=ngens(W)
      L=reflection_subgroup(W, h[:levi])
      cusp=Lusztig.FindCuspidalInLevi(h[:cuspidalName], L)
      n=findfractionalpowers(L)
      if !isnothing(n[cusp])
        uc.fractions[h[:charNumbers]].=n[cusp]
        for i in h[:charNumbers] push!(reasons[i], joindigits(h[:levi])) end
      end
    end
  end
  sers=Series(W;proper=true)
  sers=unique(filter(x->haskey(x, :mC) && iscyclic(x),sers))
  while !isempty(sers)
    for s in sers
      p = findfirst(i->uc.fractions[i]!==nothing,charnumbers(s))
      if isnothing(p) print("reticent series", s, "\n")
      else
        frac=dSeries.mC(s)[p]*dSeries.e(s)*s.d.r
        fix=modZ(uc.fractions[charnumbers(s)[p]]-frac)
        if fix!=0 println(" **** Badly normalized series ",s," adding ",fix) end
        print(rio(),"\n", s, "==>")
        for i in 1:dSeries.e(s)
          cn=charnumbers(s)[i]
          print(cn,".c")
          frac = modZ(dSeries.mC(s)[i]*dSeries.e(s)*s.d.r+fix)
          if isnothing(uc.fractions[cn]) push!(reasons[cn], s.d.r)
            uc.fractions[cn]=frac
          elseif uc.fractions[cn]!=frac
           println(rio(),"Failed! ",cn, "==", charnames(uc)[cn],
                 " in ", s, "\n     where mC modZ==", frac, "\n
                 conflicts with ", reasons[cn])
          else push!(reasons[cn],s.d.r)
          end
        end
        sers = setdiff(sers, [s])
      end
    end
  end
  for i=uc.harishChandra
    if length(i[:charNumbers])>1
      println(uc.fractions[i[:charNumbers]],length.(reasons[i[:charNumbers]]))
    else
      println(uc.fractions[i[:charNumbers][1]],reasons[i[:charNumbers][1]])
    end
    p=sort(unique(uc.fractions[i[:charNumbers]]))
    if haskey(i,:qEigen)
      if p!=[i[:qEigen]]
        if p==[nothing]
          println("!!!!!!!!!!!! for ", i[:charNumbers], 
                " qEigen should be nothing is ", i[:qEigen])
        else
          error("for HCseries $i of ",W," qEigen should be ",p)
          uc.fractions[i[:charNumbers]]=i[:qEigen].+i[:charNumbers]
        end
      end
    elseif p!=[0] println("!!!!!!!!!!!!!! qEigen unbound should be ",p)
    end
  end
  uc.fractions
end

# series of d-regular element w
function PrincipalSeries(W, d)
  s=cuspidal_data(W, d, length(relative_degrees(W, d)))
  if length(s)!=1 error(" not one ", d, "-Sylow") end
  Series(W,s[1]...)
end

function CheckMaximalPairs(W)
  res=map(d->PrincipalSeries(W, d), EigenspaceNumbers(W))
  for s in res
    e=elements(s.levi)
    any(function(x)
      local Z
      Z=centralizer(Group(s.spets), x)
      Z=Group(vcat(gens(Z), gens(Group(s.levi))))
      if length(Z)!=length(relative_group(s))*length(s.levi) error() end
      return true
    end, e)
  end
  res
end

Chevie.coefftype(p::Mvp{T1,T2}) where{T1,T2} = T2

# find the specialization of Mvp parameters leading to group algebra
# assume each parameter is the power of a single variable
function SpecializationToGroup(H)
  vcat(unique(vcat(map(H.para)do p
    e=length(p)
    res=Pair{Symbol,coefftype(p[1])}[]
    for i in 1:e
      v=scalar(p[i])
      if !isnothing(v)
        if v!=E(e,i-1) error(v," is not ",E(e,i-1)) end
        continue
      end
      if length(p[i])!=1 || length(first(monomials(p[i])))!=1
       error(repr(p[i],context=rio())," is not a power of a variable")
      end
      push!(res,variables(p[i])[1]=>root(E(e,i-1),degree(p[i])))
    end
    res
  end))...)
end

# get central characters without having to know Chartable
# for algebras above group algebra
function HeckeCentralCharacters(H)
  W=H.W
  z=length(center(W))
  v=map(m->root(m,z),central_monomials(H))
  v1=CharTable(W).irr[:,position_regular_class(W,z)]
  v2=CharTable(W).irr[:,position_class(W,one(W))]
  map((x,y,z)->x*y//z//value(x,SpecializationToGroup(H)...),v,v1,v2)
end

# for an algebra with single parameter q, find for each character φ 
# an e such that φ takes values in ℂ[q^{1/e}]
function gete(H)
  function fp(p)
    if coefftype(p)==Int || iszero(p) return 1 end
    lcm(denominator.(vcat(map(x->modZ.(powers(x)),monomials(p))...)))
  end
  ct=CharTable(H)
  if !isnothing(ct) map(x->lcm(fp.(x)),eachrow(ct.irr))
  else fp.(HeckeCentralCharacters( H ))
  end
end

function Checkzegen(W)
  l1=Series(W;proper=true)
  l=filter(s->haskey(s,:Hecke),l1)
  println(rio(),"with no Hecke:",setdiff(l1,l))
  println("    center==",gcd(degrees(W)))
  uc=UnipotentCharacters(W)
  aA=uc.a+uc.A
  srat(x)=denominator(x)==1 ? numerator(x) : x
  for s in l
    xprintln(s)
    e=gete(hecke(s))
    z=gcd(degrees(relative_group(s)))
    ucL=UnipotentCharacters(s.levi)
    aAL=(ucL.A.+ucL.a)[s.cuspidal]
    showtable(srat.(toM([modZ.((aA[charnumbers(s)].-aAL)//z),
                         modZ.(aA[charnumbers(s)]//z), modZ.(1 .//e)])), 
     row_labels=["(aA-aAL)/z", "aA/z", "local Hecke irr"])
  end
end

# minimal polynomials for Cyclotomic hecke algebras of series s
function ratser(s)
  map(x->prod(y->Mvp(:T)-y,x),collectby(degree,hecke(s).para[1]))
end

function Checkratser(W)
  ss=Series(W;proper=true)
  F=NF(cartan(Group(ss[1].spets))...)
  xprintln("F=",F)
  for s in ss
    if isnothing(hecke(s)) continue end
    p=ratser(s)
    F1=NF(F.gens...,s.d)
    p=filter(x->any(y->!(y in F1),coefficients(x)),p)
    if !isempty(p) printstyled(rio(),s.d,": ",p;color=:red) end
  end
end

function CheckPiGPiL(n)
  W=crg(n)
  s=Series(W;proper=true)
  s=filter(ser->iscyclic(ser) && ser.principal,s)
  D0(s)=sum(degrees(Group(s.spets))+codegrees(Group(s.spets)))-
        sum(degrees(Group(s.levi))+codegrees(Group(s.levi)))
  map(x->D0(x)//dSeries.e(x),s)
end

function Checkdovere(n)
  W=crg(n)
  s=Series(W;proper=true)
  s=filter(ser->iscyclic(ser) && ser.principal,s)
  filter(x->x.e>2*order(x.d),s)
end

function CheckxiL(n)
  W=crg(n)
  l=Series(W;proper=true)
  l=filter(s->iscyclic(s) && s.principal,l)
  filter(x->!iszero(modZ(Root1(PhiOnDiscriminant(x.levi)).r*x.d.r)),l)
end

# check formula for product parameters
function CheckCoN(i)
  W=crg(i)
  l=Series(W;proper=true)
  l=filter(x->length(x.levi)==1,l)
  l=filter(x->haskey(x,:Hecke),l)
  map(l)do s
    m=det_bareiss(reflrep(W, s.element))
    sg=(-1)^sum(degrees(relative_group(s)).-1)
    if length(hyperplane_orbits(relative_group(s)))>1 false
    else m*sg*prod(s.Hecke.para[1])^hyperplane_orbits(relative_group(s))[1].N_s
    end
  end
end

CheckLuCox(s)=[prod(x->1-x,filter(x->x!=1,s.Hecke.para[1])),degree(s)(Mvp(:q))]

# data for linear char:
# returns d, e_W/d (mod d), ω_c (mod d)
function datafor(W)
  l=map(d->PrincipalSeries(W, d),regular_eigenvalues(W))
  e=sum(degrees(W)+codegrees(W))
  map(s->[s.d, gcd(e*s.d,conductor(s.d)) // gcd(gcd(map(x->x[:N_s], 
          hyperplane_orbits(relative_group(s)))), conductor(s.d))], l)
end

# description of centralizers of regular elements
function cent(W)
  l=map(d->PrincipalSeries(W, d),regular_eigenvalues(W))
  res=map(s->(r=s.d, g=repr(relative_group(s);context=:limit=>true)), l)
  l=sort(unique(map(x->x.g,res)))
  res=map(x->[map(z->z[1],filter(y->y.g==x,res)), x],l)
  for x in sort(res) xprintln(x[1]," ",x[2]) end
end

function EigenspaceNumbers(W)
  d=degrees(W)
  if W isa Spets 
    e=map(p->lcm(p[1],conductor(Root1(p[2]))),d)
    e=sort(union(divisors.(e)))
    filter(x->any(p->E(x,p[1])==p[2],d),e)
  else sort(union(divisors.(d)...))
  end
end

# nrconjecture: si kd est le multiple maximum de d tq
# !isempty(relative_degrees(W,kd)), alors
# length(relative_degrees(W,d))=length(relative_wegrees(W,kd))*Phi(kd)/Phi(d)
function nrconjecture(W)
  e=setdiff(EigenspaceNumbers(W), regular_eigenvalues(W))
  e=filter(x->!(x[1] in r),e)
  for d in e
    f = filter(x->mod(x,d)==0,r)
    if length(f)>0 println(rio(),"**** failed: ",W, d, f) end
    p = First(reverse(e), (i->begin mod(i[1], d) == 0 end))
    if p!=d
      if length(RelativeDegrees(W,d))*Primes.totient(d)!=p[2]*Primes.totient(p[1])
        println(rio(),"**** failed: ",W, p, d)
      end
    end
  end
end

# An element which has a maximal ζ_d-eigenspace lives in the minimal
# d-split levi of same d-rank as W
function minimaldparabolic(W,d)
  l=map(i->setdiff(eachindex(gens(W)),[i]),eachindex(gens(W)))
  l=filter(x->length(relative_degrees(reflection_subgroup(W, x), d)) 
            ==length(relative_degrees(W, d)),l)
  map(x->reflection_subgroup(W, x),l)
end
