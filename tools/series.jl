function findfractionalpowers(W)
  if length(gens(W))==0 return [0] end
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
  divs=sort(unique(conductor.(refleigen(W))))
  res = map(d->PrincipalSeries(W, d), divs)
  for s in res
    e = elements(s.levi)
    any(function(x)
      local Z
      Z = centralizer(Group(s.spets), x)
      Z = Group(vcat(gens(Z), gens(Group(s.levi))))
      if length(Z) != length(relative_group(s)) * length(s.levi) error() end
      return true
    end, e)
  end
  res
end

function Checkzegen(W)
  l=Series(W;proper=true)
  println(rio(),"with no Hecke:",filter(s->!haskey(s,:Hecke),l))
  println("    center==",gcd(degrees(W)))
  l=filter(s->haskey(s,:Hecke),l)
  uc=UnipotentCharacters(W)
  aA=uc.a+uc.A
  for s in l
    e=gete(s[:Hecke])
    z=OrderCenter(relative_group(s))
    ucL=UnipotentCharacters(s.levi)
    aAL=(ucL.A.+ucL.a)[s.cuspidal]
    println(s)
    print(FormatTable([modZ.((aA[charNumbers(s)]-aAL)//z),
                       modZ.(aA[charNumbers(s)]//z), modZ.(1 .//e)], 
   Dict{Symbol, Any}(:rowLabels => ["(aA-aAL)/z", "aA/z", "local Hecke irr"])))
  end
end

# checks minimal polynomials for Cyclotomic hecke algebras
function CheckRatCyc(s)
  if !haskey(s,:Hecke) error(" Hecke not bound") end
  l=collectby(degree,hecke(s).para[1])
  fact=map(x->prod(y->Mvp(:T)-y,x),l)
  F=NF(cartan(Group(s.spets))...)
  fact,F
end

function CheckRatSer(arg...)
  W=crg(arg...)
  s=Series(W;proper=true)
  s=filter(x->x.spets!=x.levi,s)
  join(map(CheckRatCyc,s), "\\hfill\\break\n")
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
  map(Hecke,s)
  s=filter(ser->iscyclic(ser) && ser.principal,s)
  filter(x->e(x)//x.d>2,s)
end

function CheckxiL(n)
  W=crg(n)
  l=Series(W;proper=true)
  map(Hecke,l)
  l=filter(s->iscyclic(s) && s.principal,l)
  filter(x->!isone(Root1(PhiOnDiscriminant(x.levi))^x.d),l)
end

# check formula for product parameters
function CheckCoN(i)
  W=crg(i)
  l=Series(W;proper=true)
  l=filter(x->length(x.levi)==1,l)
  map(Hecke, l)
  l=filter(x->haskey(x.prop,:Hecke),l)
  map(l)do s
   m=DeterminantMat(reflrep(W, s.prop[:element]))
    sg=(-1)^sum(degrees(relative_group(s))-1)
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
  d = degrees(W)
  if !(W isa Spets) return union(divisors.(d)...) end
  e=map(p->lcm(p[1],conductor(Root1(p[2]))),d)
  e=union(divisors.(e))
  filter(x->any(p->E(x,p[1])==p[2],d),e)
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
