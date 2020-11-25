# This file contains code to find Ennola-permutation,
# and code to find Fourier matrices (detfam)
#------------------ Ennola -----------------------------------------
# Ennola-theory axioms:
# Let z=|ZW|. We call xi-Ennola the map q -> E(z)^xi q.
# (This corresponds to Ennola_E(z)^-1 of "towards Spetses I").
# - For any character sheaf A there exists e_A such that xi-Ennola(A)=e_A^xi A
# Thus  Ennola stabilizes  families and almost-HC-series;
# -  If  A  has  cuspidal  data  (λ,χ) where λ is a cuspidal
# character  sheaf  of  L  and  χ∈ Irr(W_G(L,λ)),  there  exists
# E_λ such that e_A=E_λ*ω_χ(E(z))
#
# For  a split Spetses,  cuspidal character sheaves  correspond to cuspidal
# unipotent characters; if U∈ Uch(G) has cuspidal data
# χ∈ Irr(W_G(L,λ)) then we have
# Frob(Ennola_xi(U))/Frob(U)=omegachi(E(z))^-xi E(z^2)^(xi^2*(a+A))E_λ^xi

# positions(-with-sign) in l where o or -o appears
positionssgn(l,o)=vcat(findall(==(o),l),-findall(==(-o),l))
ordercenter(W)=prod(map(t->gcd(degrees(t)),refltype(W)))

# EnnolaBete(W,i[,xi])
# possible action of xi-Ennola on i-th family
# (here xi in 0:ordercenter(W)], default 1)
# list of possible destinations for each char, taking just degree in account
function EnnolaBete(W,i,xi=1)
  xi=E(ordercenter(W),xi)
  ud=Uch.CycPolUnipotentDegrees(W)[UnipotentCharacters(W).families[i][:charNumbers]]
  map(p->positionssgn(ud, ennola_twist(p, xi)), ud)
end

# List of ω_χ(E(z)) where z=ordercenter(W) for character sheaves (ρ,χ)
function OmegaChi(uc)
  gets(uc,:omegachi)do
    relativeFake=fill(Pol(0),length(uc))
    for h in uc.harishChandra
      relativeFake[h[:charNumbers]]=fakedegrees(reflection_group(
                                    h[:relativeType]), Pol())
    end
    xi0 = E(ordercenter(Group(spets(uc))))
    map(f->f(xi0)//f(1),relativeFake)
  end
end

# if Ennola_xi(U_χ=ρ_i)=ρ' returns deduced Frob(ρ')
function predeigen(W,i,xi)
  if i>nconjugacy_classes(W) error("only for principal series") end
  z=ordercenter(W)
  uw=UnipotentCharacters(W)
  OmegaChi(uw)[i]^-xi*E(z^2,xi^2*(uw.prop[:a][i]+uw.prop[:A][i]))*eigen(uw)[i]
end

# EnnolaBE(W,i[,xi])
# More sophisticated Ennola taking in account the predicted eigenvalues for
# Ennola(principal series).
function EnnolaBE(W,i,xi=1)
  uc=UnipotentCharacters(W)
  f=uc.families[i]
  eig=eigen(f)
  n=uc.harishChandra[1][:charNumbers]
  e=EnnolaBete(W, i, xi)
  N=prod(length,e)
  if N==1 return e end
  print(N)
  for j in 1:length(f)
    if !(f[:charNumbers][j] in n) continue end
    e[j]=filter(k->eig[abs(k)]==predeigen(W,f[:charNumbers][j],xi),e[j])
  end
  if prod(length,e)<N
    print("==>",prod(length,e))
    N=prod(length,e)
  end
  if N==1 println();return e end
  for j in 1:length(f)
    if !(f[:charNumbers][j] in n) || length(e[j])!=1 continue end
    d=abs(e[j][1])
    e[d]=filter(k->eig[abs(k)]==predeigen(W,f[:charNumbers][j],2xi),e[d])
  end
  if prod(length,e)<N
    print("==>",prod(length,e))
    N=prod(length,e)
  end
  if N==1 println();return e end
  if ordercenter(W)>2xi
    for j in 1:length(f)
      if !(f[:charNumbers][j] in n) || length(e[j])!=1 continue end
      d = abs(e[j][1])
      if length(e[d])!=1 continue end
      d = abs(e[d][1])
      e[d]=filter(k->eig[abs(k)]==predeigen(W,f[:charNumbers][j],3xi),e[d])
    end
  end
  if prod(length,e)<N
    print("==>",prod(length,e))
    N=prod(length,e)
  end
  if N==1 println();return e end
  if ordercenter(W)>3xi
    for j in 1:length(f)
      if !(f[:charNumbers][j] in n) || length(e[j])!=1 continue end
      d=abs(e[j][1])
      if length(e[d])!=1 continue end
      d=abs(e[d][1])
      if length(e[d])!=1 continue end
      d=abs((e[d])[1])
      e[d]=filter(k->eig[abs(k)]==predeigen(W,f[:charNumbers][j],4xi),e[d])
    end
  end
  if prod(length,e)<N print("==>",prod(length,e)) end
  println()
  e
end

struct Ennola
  prop::Dict{Symbol,Any}
end

Base.getindex(f::Ennola,k)=f.prop[k]
#Base.haskey(f::Ennola,k)=haskey(f.prop,k)
Base.setindex!(f::Ennola,v,k)=setindex!(f.prop,v,k)

function Ennola(W)
  uc=UnipotentCharacters(W)
  e=Ennola(Dict{Symbol,Any}(:W =>W))
  e[:scalars]=Vector{Vector{Cyc{Int}}}(undef,length(uc.harishChandra))
  e[:scalars][1]=[1]
  e[:z]=ordercenter(W)
  e[:families]=map(i->Dict{Symbol,Any}(:poss=>map(j->EnnolaBE(W,i,j),
                                          1:e[:z]-1)),1:length(uc.families))
  e[:hc]=map(n->findfirst(h->n in h[:charNumbers],uc.harishChandra),1:length(uc))
  e
end

#function Base.show(io::IO,e::Ennola)
#  println(io,"Ennola scalars for cuspidals of ",e[:W])
#  rowlab=map(h->fromTeX(io,h[:cuspidalName]),
#                    UnipotentCharacters(e[:W]).harishChandra)
#  format(io,map(x->[x],e[:scalars]),rows_label="Name",row_labels=rowlab,
#        col_labels=[1//e[:z]])
#end

# Ennola for i-th family as an Ls
function PossToLs(e::Ennola, i, xi)
  r=deepcopy(e[:families][i][:poss][xi])
  if prod(length,r)>1000000 return end
  function gg(l, poss)
    local x, res
    if length(poss)==0 return [Int[]] end
    res=[]
    for x in poss[1]
      if !(abs(x) in l)
        append!(res,map(v->vcat([x],v),
               gg(vcat(l,[abs(x)]),poss[2:length(poss)])))
      end
    end
    return res
  end
  return gg(Int[], r)
end

# if Ennola_xi(U_i)=U_j returns deduced En(λ)^xi for cuspidal of U_i
function EigenEnnola(W, i, j, xi)
  z = ordercenter(W)
  uw = UnipotentCharacters(W)
  (OmegaChi(uw)[i]^-xi*E(z^2,xi^2*(uw.prop[:a][i]+uw.prop[:A][i]))*eigen(uw)[i])//eigen(uw)[j]
end

function tw(W, j)
  e = ennola(W)[1][:ls]^j
  uc = UnipotentCharacters(W)
  n = CharNumbers(PrincipalSeries(W, 1))
  map(i->EigenEnnola(W, i, abs(i^e), j), n)
end

# in e set xi-ennola scalar of cuspidal corresp. to charno to be in scal
function SetScal(e::Ennola, scal, charno, xi)
  k=e[:hc][charno]
  if xi==1 && !isassigned(e[:scalars],k) e[:scalars][k] = scal
  else
    n=filter(x->x^xi in scal,e[:scalars][k])
    if e[:scalars][k]!=n
      InfoChevie("\n   # ",xi,"-Ennola scalars[",k,"] ",e[:scalars][k],"->",n)
    end
    e[:scalars][k] = n
  end
end

function GetScal(e::Ennola, charno, xi)
  unique(sort(map(x->x^xi,e[:scalars][e[:hc][charno]])))
end

# restrict ennola-scalars from frobenius eigenvalues using
# knowledge on ennola-ps and
# En(\lambda)^k=EigenEnnola(W,i,j,k) if Ennola_xi(U_i)=U_j
function ScalarsFromFrob(e::Ennola)
  W=e[:W]
  uc=UnipotentCharacters(W).families
  z=ordercenter(W)
  for i in 1:length(e[:families])
    f=uc[i][:charNumbers]
    for xi in 1:z-1
      for j in 1:length(f)
        SetScal(e,unique(sort(map(x->EigenEnnola(W, f[j], f[abs(x)],xi),
                                  e[:families][i][:poss][xi][j]))), f[j], xi)
      end
    end
  end
end

# get info on ennola-ps from ennola-scalars and frobenius-eigenvalues
function PossFromScal(e::Ennola)
  W = e[:W]
  for i in 1:length(e[:families])
    f = UnipotentCharacters(W).families[i][:charNumbers]
    for xi in 1:e[:z]-1
      r=e[:families][i][:poss][xi]
      for j in 1:length(f)
        s=GetScal(e, f[j], xi)
        l=filter(n->EigenEnnola(W,f[j],f[abs(r[j][n])],xi) in s,1:length(r[j]))
        if length(l)<length(r[j])
          InfoChevie("   # Family#", i, " En", xi//e[:z],
             "[",j,"]:restricted ",string(r[j]...),"->",string(r[j][l]...),"\n")
        end
        r[j]=r[j][l]
      end
    end
  end
end

# ScalarsFromFourier(e,fno[,fourierMat)
function ScalarsFromFourier(e::Ennola,fno,F=nothing)
  uc=UnipotentCharacters(e[:W])
  f =uc.families[fno]
  o =OmegaChi(uc)[f[:charNumbers]]
  if isnothing(F)
    F=f[:fourierMat]
    if haskey(f,:lusztig) && f[:lusztig]
      F=^(F,Perm(conj.(f[:eigenvalues]), f[:eigenvalues]);dims=1)
    end
  end
  F=conj.(F)
  n=size(F,1)
  for xi in 1:e[:z]-1
    poss=e[:families][fno][:poss][xi]
    es=map(n->GetScal(e, n, xi), f[:charNumbers])
    F1=hcat(map((x,y)->x.*y^xi, eachcol(F),o)...)
    for m in 1:n
      poss[m]=filter(poss[m])do j
        all(1:n)do k
          local a, b
          a=F[abs(j),k]
          b=F1[m,k]
          if b==0
            if a isa Mvp a=scal(a) end
            return isnothing(a) || a==0
          end
          a=a.//b
          if a isa Mvp a=scal(a) end
          if isnothing(a) return true end
          return sign(j)*a in es[k]
        end
      end
      if length(poss[m]) == 0 return false end
      if length(poss[m]) == 1
        j = poss[m][1]
        l = filter(k->F[m,k]!=0,1:n)
        for k in l
          SetScal(e,[(sign(j)*F[abs(j),k])//F1[m,k]],f[:charNumbers][k],xi)
        end
      end
    end
  end
  true
end

# EnnolafromFusion(W,i[,poss])
# possibilities of Ennola from fusion algebra. If poss given
# tries to fit with poss.
function EnnolafromFusion(W,i,poss=nothing)
  A=fusion_algebra(UnipotentCharacters(W).families[i])
  b=basis(A)
  res=vcat(map(function(x)
    local p
    p=SPerm(b,Ref(x).*b)
    if !isnothing(p) return [p, SPerm(-vec(p))]
    else return SPerm[]
    end
  end, b)...)
  if isnothing(poss)
    error(" not  expected")
    return res
  end
  res=map(p->filter(x->all(j->vec(x)[j] in p[j],1:length(vec(x))),res), poss)
  good=filter(i->length(res[i])>0,1:length(res))
  if length(good) == 0 error("**** no solution for family ", i) end
  res[good[1]]
end

function ennola(W)
  if haskey(W.prop,:ennola) return W.prop[:ennola] end
  if ordercenter(W)==1
    print("W has trivial center!\n")
    return false
  end
  e=Ennola(W)
  ScalarsFromFrob(e)
  PossFromScal(e)
  uc=UnipotentCharacters(W)
  fam=uc.families
  for i in 1:length(fam) ScalarsFromFourier(e, i) end
  PossFromScal(e)
  ff=filter(i->any(j->length(j)>1,e[:scalars][e[:hc][fam[i][:charNumbers]]]),
            1:length(fam))
  le=map(s->Ennola(Dict{Symbol, Any}(:W=>e[:W],:z=>e[:z],
     :families=>deepcopy(e[:families]),:scalars=>map(x->[x],s),:hc=>e[:hc])),
           cartesian(e[:scalars]...))
  le = filter(e-> all(i->ScalarsFromFourier(e, i), ff) ,le)
  W.prop[:ennola] = map(le)do e
    local res, f
    res=Dict{Symbol, Any}(:scalars=>map(x->x[1],e[:scalars]))
    res[:fls]=map(enumerate(e[:families]))do (i, r)
      local p
      p=EnnolafromFusion(W, i, r[:poss])
      if length(p)==0 error("AAAAA") end
      if length(p)>1 InfoChevie("   # ",length(p)," sols in family ",i,"\n") end
      r[:poss]=p
      return p[1]
    end
    res[:ls]=prod(i->res[:fls][i]^mappingPerm(1:length(fam[i]),
                                      fam[i][:charNumbers]),1:length(fam))
    res[:allscal]=HasType.DiagonalOfMat(Matrix(res[:ls])^fourier(uc))
    res
  end
  W.prop[:ennola]
end

# convert family-wise ennola ps to global perm
function famEn2En(W, e)
  uc=UnipotentCharacters(W)
  res=fill(0,length(uc))
  for (f,p) in zip(uc.families,e)
    res[f[:charNumbers]]=f[:charNumbers]^p
  end
  res
end

# Twist a subspets by xi, an element of the center of the (untwisted) parent
function EnnolaTwist(HF, xi)
  H=Group(HF)
  W=parent(H)
  w=Representative(ConjugacyClasses(W)[PositionRegularClass(W, AsRootOfUnity(xi))])
  SubSpets(HF[:parent],H[:rootInclusion][H[:generatingReflections]],w*HF[:phi])
end

# Predict what should be LusztigInductionTable(HF,WF) by Ennola-twisting that
# of EnnolaTwist(HF,xi^-1);
function PredictRLGByEnnola(HF, WF, xi)
  function GlobalEnnolaLs(W, xi)
    local o
    o = AsRootOfUnity(xi) * ordercenter(W)
    return famEn2En(W, map(x->x^o,ennola(W)[1][:ps]))
  end
  HF = EnnolaTwist(HF, xi ^ -1)
  H = Group(HF)
  t = LusztigInductionTable(HF, WF)[:scalar]
  if xi ^ ordercenter(H) == 1
      ps = GlobalEnnolaLs(H, xi)
      t = Permuted(TransposedMat(t), PermList(map(abs, ps)) ^ -1)
      t = TransposedMat(map((x->begin sign(ps[x]) * t[x] end), 1:length(t)))
  end
  ps = GlobalEnnolaLs(Group(WF), xi)
  t = Permuted(t, PermList(map(abs, ps)) ^ -1)
  map(x->sign(ps[x])*t[x],1:length(t))
end

function CheckRLGByEnnola(HF, WF, xi)
  H =Group(HF)
  t =PredictRLGByEnnola(HF, WF, xi)
  fW=UnipotentCharactersOps[:Fourier](UnipotentCharacters(WF))
  fH=UnipotentCharactersOps[:FourierInverse](UnipotentCharacters(HF))
  t =fW*t*fH
  pieces = LusztigInductionTable(HF, WF, true)
  print("Checking ", pieces, "\n")
  pieces = pieces[:pieces]
  for i in 1:length(pieces)
    p = pieces[i]
    p[:ns] = t[p[:wnum]][p[:hnum]]
    t[p[:wnum]][p[:hnum]]=0*t[p[:wnum]][p[:hnum]]
    s = ProportionalityCoefficient(vcat(p[:scalar]...), vcat(p[:ns]...))
    if s!=false
      if s!=1 print("got piece[", i, "] upto ", s, "\n") end
    else
      r = TransposedMat(p[:scalar])
      p[:ns] = TransposedMat(p[:ns])
      r=map(i->ProportionalityCoefficient(r[i],p[:ns][i]), 1:length(r))
      if any(x->x==false,r) error("extensions not proportional", r) end
      print(" *** piece no.", i, " of ", length(pieces), " ***\n")
      print((p[:head])(p, Dict{Symbol, Any}()))
      j = Filtered(1:length(r), (i->begin r[i] != 1 end))
      for k in j
        print("char. ",k,"==",p[:charSubgroup][k]," of ",p[:u]," extension differs by ",r[k],"\n")
      end
    end
  end
  if gapSet(Concatenation(t))!=[0] error("CheckRLGByEnnola") end
end

# Check that an inner Ennola E_xi sends RLG to RE_xi(L)E_xi(G)
function CheckEnnolaRLG(W)
  function twist(HF, w)
    local r
    r=Group(HF)[:rootInclusion][Group(HF)[:generatingReflections]]
    return SubSpets(HF[:parent], r, w*HF[:phi])
  end
  WF = Spets(W)
  c = ordercenter(W)
  l = Concatenation(map(x->twistings(W,x), ProperParabolics(W)))
  e = ennola(W)[1][:ps]
  uc = map(x->x[:charNumbers],UnipotentCharacters(W).families)
  ps = []
  for p in 1:length(uc)
    ps[uc[p]]=map(x->sign(x)*uc[p][abs(x)],e[p])
  end
  for xi in 1:c-1
    for p in Filtered(x->x[1]>=x[2],map(x->[x, EnnolaTwist(x, E(c, xi))],l))
      HF = p[1]
      HFxi = p[2]
      n = ReflectionName(HF)
      nxi = ReflectionName(HFxi)
      e = ennola(Group(HF))[1][:ps]
      uc = map(x->x[:charNumbers], UnipotentCharacters(Group(HF)).families)
      psH = []
      for i in 1:length(uc)
        psH[uc[i]] = map(x->sign(x)*uc[i][abs(x)], e[i])
      end
      t = TransposedMat(LusztigInductionTable(HF, WF)[:scalar])
      t = Permuted(t, PermList(map(abs, psH)) ^ -1)
      t = TransposedMat(map(x->sign(psH[x])*t[x], 1:length(t)))
      t = Permuted(t, PermList(map(abs, ps))^-1)
      t = map(x->sign(ps[x])*t[x],1:length(t))
      txi = LusztigInductionTable(HFxi, WF)
      if txi[:scalar] != t
        p = deepcopy(txi)
        p[:scalar] = t
        CHEVIE[:Check][:EqTables](p, txi)
        error("CheckEnnolaRLG")
      else
        xprint(E(c,xi),":",HF,"==>",HFxi," OK!\n")
      end
    end
  end
end

# The Fourier matrix F times unipotent  degrees are the fake degrees. If
# p=Ennola-Permutation  and e=Ennola-scalars  on character  sheaves then
# one has Fp=eF
# This program,  on i-th family  of W, using partially  computed fourier
# matrix  M, and  a  given permutation-with-signs  p,  tries to  fill-in
# Ennola scalars e. (start with some known scalars in e).
# Return false if some incompatibility discovered
function scalpred(xi, i, W, p, e, M)
  uc = UnipotentCharacters(W)
  f = uc[:families][i]
  n = length(f[:charNumbers])
  e = copy(e)
  for j in 1:n
    i=findfirst(i->!(M[:Valueat](i, j) in [false, 0]) && M[:Valueat](abs(p[i]), j)!=false,1:n)
    if i != false
      scal=ComplexConjugate((M[:Valueat](abs(p[i]),j)*sign(p[i]))//M[:Valueat](i, j))
      if !(e[j] !== nothing) e[j] = scal
      elseif e[j] != scal return false
      end
    end
  end
  return e
end

# Add  to M  (linear  relations for  fourier)  resulting from  xi-Ennola
# diagonal on CS on i-th family, using a given sgnperm p for xi-Ennola
function relennolaxi(xi, M, p)
  W =M[:W]
  uc=UnipotentCharacters(W)
  f =uc.families[M[:fno]]
  n =length(f)
  res=[]
  for j in 1:n
    ijknown=filter(i->!(M[:Valueat](i,j) in [false,0]) && p[i]!=false,1:n)
    om=findfirst(i->i in ijknown && M[:Valueat](abs(p[i]),j)!=false,1:n)
    val =GetScal(M[:ennola], f[:charNumbers][j], xi)
    if length(val)>1
      if om != false
        om = (M[:Valueat](abs(p[om]),j)*sign(p[om]))//M[:Valueat](om, j)
        SetScal(M[:ennola], [ComplexConjugate(om)//
            OmegaChi(uc)[f[:charNumbers][j]]^xi], f[:charNumbers][j], xi)
      else xprint(" cannot determine ", E(ordercenter(W),xi),"-scalar for ",
        TeXStrip(uc[:TeXCharNames][f[:charNumbers][j]]),"\n")
      end
    else om=ComplexConjugate(val[1]*OmegaChi(uc)[f[:charNumbers][j]]^xi)
    end
    if om != false
      for i in Filtered(1:n, (i->begin p[i] != false end))
        if !(Makerel(M,[[M[:at](i,j),om],[M[:at](abs(p[i]),j),-sign(p[i])]],0))
          return false
        end
      end
    else
      for i in ijknown
        for i1 in filter(k->k>i,ijknown)
          if !(Makerel(M, [[M[:at](abs(p[i]),j), sign(p[i])//M[:Valueat](i,j)],
            [M[:at](abs(p[i1]),j), -sign(p[i1]) // M[:Valueat](i1, j)]], 0))
              return false
          end
        end
      end
    end
  end
  print(length(res), " relations ")
  return true
end

# sq(M) get Mvp matrix from linear relations M
function sq(M)
  v = Values(M)
  map((i->begin map(j->v[M[:at](i, j)], 1:M[:n]) end), 1:M[:n]) * Mvp("a") ^ 0
end

# relconj(W,i,M[,p]) add to M for i-th family of W relations from S^2=complex
# conjugation; if given p=sgnperm for Ennola
function relconj(arg...,)
  W = arg[1]
  i = arg[2]
  M = arg[3]
  print("conj... ")
  uc = UnipotentCharacters(W)
  f = (uc[:families])[i]
  n = length(f[:charNumbers])
  if length(arg) < 4
# possibilities for the effect of complex conjugation on ith family
# (list of possible destinations for each char)
      ud = (CycPolUnipotentDegrees(W))[f[:charNumbers]]
      p = map((p->begin
                      PositionsSgn(ud, ComplexConjugate(p))
                  end), ud)
      m = sq(M) ^ 2
      for i in 1:n
          for j in 1:n
              if iszero((m[i])[j] - 1)
                  for k in 1:n
                      if k == i p[i] = [j]
                      else p[k] = Difference(p[k], [j])
                      end
                  end
              end
          end
      end
  else p = arg[4]
  end
  print("possconj==", p, "\n")
  for i in 1:n
      for j in Difference(1:n, map(abs, p[i]))
          v = map((k->begin (M[:Valueat])(k, j) end), 1:n)
          if all((x->begin x != false end), v)
              Makerel(M, map((k->begin [(M[:at])(i, k), v[k]] end), 1:n), 0)
          end
      end
  end
end

# add to M relations expressing M=TransposedMat(m)
function relsym(M)
  for i in 1:M[:n]
      for j in 1:i - 1
          Makerel(M, [[(M[:at])(i, j), 1], [(M[:at])(j, i), -1]], 0)
      end
  end
end

# add to M relations expressing M[i]*v=c
relprod(M, v, i, c)=Makerel(M, map(j->[M[:at](i, j), v[j]], 1:length(v)), -c)

# add to M relations expressing M*vp=vq where vp, vq are vectors of polynomials
function relpol(M, vp, vq)
  function coeff(p, i)
          if i < p[:valuation] || i > degree(p) return 0 end
          return (p[:coefficients])[(i - p[:valuation]) + 1]
      end
  max = maximum(map(degree, Concatenation(vp, vq)))
  for i in 1:M[:n]
      for k in 0:max
          if !(relprod(M, map((p->begin coeff(p, k)
                              end), vp), i, coeff(vq[i], k)))
              return false
          end
      end
  end
  return applynew(M)
end

# add linear relations coming from Sbar*S=Id and known lines of S
function relunitary(M,)
  m=map(i->map(j->M[:Valueat](i,j), 1:M[:n]),1:M[:n])
  for i in 1:M[:n]
    for j in i:M[:n]
      if all(k->m[j][k] != false || m[i][k]!=false && iszero(m[i][k]), 1:M[:n])
        v=map(x->x==false ? 0 : conj(x), m[j])
        if i!=j v=relprod(M,v,i,0)
        else v=relprod(M,v,i,1)
        end
        if !v return false end
      end
    end
  end
  applynew(M)
end

function fromunitary(M)
  n=sq(M)
  return gapSet(Concatenation(n*conj.(n)-IdentityMat(length(n))))
end

function CheckMagrees(M)
  if M[:error] != true return end
  f = UnipotentCharacters(M[:W]).families[M[:fno]]
  four = f[:fourierMat]
  checkMagrees(M, Concatenation(map(i->four[i][1:i],1:M[:n])),function(p)
#inverse of at: [i,j] from linear position
      local i, j
      i=First(1:p, (i->begin (i*(i+1))//2>=p end))
      j=p-(i*(i-1))//2
      return [i,j]
  end)
end

function initdetfam(W, fno)
  uc = UnipotentCharacters(W)
  q = X(Cyclotomics)
  f = (uc[:families])[fno]
  n = length(f[:charNumbers])
  M = LinSys((n * (n + 1)) // 2)
  M[:n] = n
# linear position of element [i][j] in symmetric matrix represented linearly
  M[:varnames] = Concatenation(map(i->map(j->SPrint("x", i, "_", j),1:i), 1:n))
  M[:at] = function (i, j)
          if i >= j return (i * (i - 1)) // 2 + j
          else return (j * (j - 1)) // 2 + i
          end
      end
  M[:Valueat] = function (i, j) return Value(M, (M[:at])(i, j)) end
  M[:W] = W
  M[:fno] = fno
  M[:error] = true
  M[:sqdisp] = function (M,)
          local z
          if length(M[:relations]) != 0
              z = TransposedMat(sq(M))
              z = map((x->begin map(Format, x) end), z)
              push!(z, (CharNames(uc))[f[:charNumbers]])
          end
          if IsBound(z) && length(z) < 15
              print(Format(TransposedMat(z)), "\n")
          else
              print(Format((CharNames(uc))[f[:charNumbers]]), "\n")
          end
      end
  return M
end

# ennoladetfam(M [,ps]) if ps for Ennola not given uses what holds from degrees
function ennoladetfam(arg...,)
  M = arg[1]
  newknown = length(M[:known]) + length(M[:relations])
  if !(haskey(M, :ennola)) M[:ennola] = Ennola(M[:W]) end
  ScalarsFromFrob(M[:ennola])
  PossFromScal(M[:ennola])
  if length(arg) > 1 e = arg[2] end
  while true
      oldknown = newknown
      z = ordercenter(M[:W])
      v = 1:z - 1
      SortParallel(map((x->begin OrderMod(x, z) end), v), v)
      for i = reverse(v)
          if IsBound(e) p = e ^ i
          else p = map(function (p,)
                          if length(p) == 1 return p[1] else return false end
                      end, ((((M[:ennola])[:families])[M[:fno]])[:poss])[i])
          end
          print("Ennola ", Format(E(z, i)), " diagonal on CS==>")
          if !(relennolaxi(i, M, p)) return false end
          applynew(M)
          CheckMagrees(M)
      end
      relunitary(M)
      newknown = applynew(M)
      if newknown == oldknown break end
  end
  CheckMagrees(M)
end

# detfam(W,fno[,do not compare with stored data])
# relations, where CC=ComplexConjugate:
# S=Transposed(S), S*ud=fd, S*fd=CC(ud) (since S^-1=CC(S))
# S*CC(O)*ud=O*ud [since sh=CC(O)SCC(O) preserves ud]
# S*O*ud=CC(O)*S*fd [since S^2CC(O)ud=CC(O)S^2ud =>SOud=CC(O)S^2ud]
function detfam(arg...,)
  W = arg[1]
  fno = arg[2]
  M = initdetfam(W, fno)
  M[:error] = length(arg) < 3
  uc = UnipotentCharacters(W)
  q = X(Cyclotomics)
  f = (uc[:families])[fno]
  O = Eigenvalues(f)
  fd = (FakeDegrees(uc, q))[f[:charNumbers]]
  ud = (UnipotentDegrees(W, q))[f[:charNumbers]]
  print("S*fd==udbar==>")
  relpol(M, fd, ComplexConjugate(ud))
  CheckMagrees(M)
  ennoladetfam(M)
  print("S*ud==fd==>")
  relpol(M, ud, fd)
  print("S*O^-1*ud==O*ud==>")
  relpol(M, map(*,ComplexConjugate(O),ud),map(*,O,ud))
  print("S*O*ud==O^-1*ComplexConjugate(ud)==>")
  relpol(M, map(*, O, ud), ComplexConjugate(map(*, O, ud)))
  ennoladetfam(M)
  M[:sqdisp](M)
  return M
end

# the Ennola theory says that if E:=DiagonalMat(Ennola-scalars)
# and P:=matrix of permutation-with signs that Ennola does and M=fourier,
# then one has M*P=E*M
# here p is a Ls
# result is a list of triples [i,e,-e] meaning E[i][i]=e or -e
function newratios(M, p)
  IsMonomial = (x->begin length(x[:coeff]) == 1 end)
  IsScal = (x->begin ScalMvp(x) != false end)
  m = sq(M)
  res = []
  for j = 1:length(m)
    g = map(i->[sign(p[i]) * m[abs(p[i])][j], m[i][j]], 1:length(m))
    if !(any((x->begin all(IsScal, x) end), g))
      pos1 = map(x->x[2]//x[1],
                 Filtered(g, (x->begin IsMonomial(x[2]) && IsScal(x[1]) end)))
      pos2 = map((x->begin x[1] // x[2] end), Filtered(g, (x->begin
                          IsMonomial(x[1]) && IsScal(x[2]) end)))
      letter = Intersection(map(x->x[:elm][1][:elm][1],pos1),
                            map((x->begin (((x[:elm])[1])[:elm])[1] end), pos2))
      if length(letter) > 0
        pos1 = First(pos1, (x->begin ((x[:elm][1])[:elm])[1] == letter[1] end))
        pos2 = First(pos2, (x->begin ((x[:elm][1])[:elm])[1] == letter[1] end))
        val = GetRoot((pos2[:coeff])[1] // (pos1[:coeff])[1], 2)
        push!(res, [j, val, -val])
      end
    end
  end
  return res
end

function ratios(M, p)
  res = []
  function f(j)
    local i, v1, v2
    for i in 1:M[:n]
      v1 = M[:Valueat](abs(p[i]), j)
      v2 = M[:Valueat](i, j)
      if v1 != false && (v2 != false && !v1 == 0v1)
          return (sign(p[i]) * v1) // v2
      end
    end
    return false
  end
  for j = 1:M[:n]
      v = f(j)
      if v != false res[j] = v end
  end
  return res
end

function tryratio(W, fno, M, p, pos, ratio)
  f = ((UnipotentCharacters(W))[:families])[fno]
  oldratios = ratios(M, p)
  newratios = copy(oldratios)
  newratios[pos] = ratio
  while newratios != oldratios
    inds = Filtered(1:M[:n], (i->begin
               newratios[i] !== nothing && !(oldratios[i] !== nothing) end))
    for pos = inds
      for i = 1:M[:n]
        Makerel(M,[[M[:at](abs(p[i]),pos),1],[M[:at](i,pos),-sign(p[i])*newratios[pos]]], 0)
      end
    end
    applynew(M)
    oldratios = newratios
    newratios = ratios(M, p)
  end
end

function ratiosmonomialmat(m, p)
  local A, n
  n = length(m)
  A = map(function (j,)
    local v
    v = Filtered(1:n, (i->begin !(iszero((m[i])[j])) end))
    return gapSet(map((i->begin [m[abs(p[i])][j]*sign(p[i]), m[i][j]] end), v))
  end, 1:n)
  A = map(function (v,)
              local ratio, c
              ratio = false
              for c = v
                  if IsCyc(c[2]) || (((c[2])[:elm])[1])[:elm] == []
                      ratio = c[1] // c[2]
                  end
              end
              if ratio != false
                  v = map((c->begin c[1] - ratio * c[2] end), v)
              end
              return gapSet(v)
          end, A)
  return A
end

function Makerelandconj(M, v, rhs)
  n = M[:total] // 2
  if !(Makerel(M, v, rhs)) return false end
  if !(Makerel(M, map(function (x,)
                      if x[1] > n return [x[1] - n, ComplexConjugate(x[2])]
                      else return [x[1] + n, ComplexConjugate(x[2])]
                      end
                  end, v), ComplexConjugate(rhs)))
      return false
  end
  return true
end

function detfrob2(M)
  newknown = applynew(M)
  n = M[:total] // 2
  print("S*T*S==Tbar*Sbar*Tbar==>")
  while true
    oldknown = newknown
    for i = Filtered(M[:known], (x->begin x <= n end))
      for j = Filtered(M[:known], (x->begin x <= n && x <= i end))
        if !(Makerelandconj(M, map((k->begin [k, M[:S][i][k]*M[:S][k][j]]
          end), 1:n), -(ComplexConjugate(Value(M,i)*Value(M,j)*M[:S][i][j]))))
            return false
        end
      end
    end
    for i = Filtered(M[:known], (x->begin x <= n end))
      for j = Difference(1:n, M[:known])
        if !(Makerelandconj(M, Concatenation(map((k->begin
              [k, ((M[:S])[i])[k] * ((M[:S])[k])[j]]
           end), 1:n), [[j + n, -(ComplexConjugate(Value(M, i))) * ComplexConjugate(((M[:S])[i])[j])]]), 0))
            return false
        end
      end
    end
    newknown = applynew(M)
    if newknown == oldknown break end
  end
end

# detfrob (W,fno[,fouriermat])
# (uses fourierMat for fno-th family of W if arg[3] not given)
function detfrob(arg...,)
  W = arg[1]
  fno = arg[2]
  uc = UnipotentCharacters(W)
  function cdeg(deg, p)
    if deg < p[:valuation] || deg > degree(p) return 0
    else return p[:coefficients][1+deg-p[:valuation]]
    end
  end
  f = uc.families[fno]
  if length(arg) == 2 S = f[:fourierMat] else S = arg[3] end
  n = length(f[:charNumbers])
  M = LinSys(2n)
  hc = map((x->begin PositionProperty(uc[:harishChandra], (h->begin
                              x in h[:charNumbers] end)) end), f[:charNumbers])
  M[:known] = Filtered(1:n, (i->begin
                  uc[:harishChandra][hc[i]][:relativeType]!=[] end))
  M[:values]=map(i->uc[:harishChandra][hc[i]][:eigenvalue], M[:known])
  M[:known] = Concatenation(M[:known], M[:known] + n)
  M[:values] = Concatenation(M[:values], ComplexConjugate(M[:values]))
# relations coming from Shintani principal series
# si p= indices de la serie principale alors (Sbar*T*ud){p}=ud{p}
  print("Shintani principal series==>\n")
  q = X(Cyclotomics)
  ud = UnipotentDegrees(W, q)[f[:charNumbers]]
  m = maximum(map(degree, ud))
  for j = Filtered(1:n, (i->begin
         (f[:charNumbers])[i] in ((uc[:harishChandra])[1])[:charNumbers] end))
      for deg = 0:m
          if !(Makerelandconj(M, map((i->begin
                           [i, ComplexConjugate((S[j])[i]) * cdeg(deg, ud[i])]
                              end), 1:n), -(cdeg(deg, ud[j]))))
              return false
          end
      end
  end
  M[:S] = S
  if length(arg) == 2
      checkMagrees(M, Concatenation(f[:eigenvalues], ComplexConjugate(f[:eigenvalues])))
  end
  return M
end

# get Frobenius from result of detfrob, trying to get more values by
# using relation x*xbar=1
function Frob(M)
# input: 2 vectors x, xbar  of Mvp. Try to use whenever x.xbar is monomial
# to get more values as Laurent monomials
  function optimize(t)
    local e, v, i
    while true
        e = map((j->begin (t[1])[j] * (t[2])[j] end), 1:n)
        v = map(ScalMvp, e)
        if all((x->begin x == 1 end), v) return t end
        if any((x->begin x != false && x != 1 end), v) return false end
        i = PositionProperty(e, (x->begin
              length(x[:elm]) == 1 && length(((x[:elm])[1])[:elm]) != 0 end))
        if i == false return t end
        e = e[i]
        v = (((e[:elm])[1])[:elm])[1]
        t = map((r->begin map(c->Value(c, [v, Mvp(v) // e]), r) end), t)
    end
  end
  res = Values(M)
  n = length(res) // 2
  res = [res[1:n], res[n + (1:n)]] * Mvp("a") ^ 0
  t = map(FactorizeQuadraticForm, gapSet(map((x->begin Product(x) - 1
                      end), TransposedMat(res))))
  t = Filtered(t, (x->begin x != false end))
  t = map((x->begin map(solve, x) end), t)
  t = map(gapSet, cartesian(t))
  print("t==", t, "\n")
  t = Filtered(t, (x->begin length(gapSet(map(y->y[1], x))) == length(x) end))
  print("t==", t, "\n")
  if length(t) == 0 return res
  else
    res=map(x->map(a->map(b->Value(b, Concatenation(x)), a), res), t)
  end
  return Filtered(map(optimize, res), (x->begin x != false end))
end

# ambiguity group of fno-th family of W (permutations fixing unipotent degrees)
function ambig(W, fno)
  uc = UnipotentCharacters(W)
  f = uc.families[fno]
  ud = CycPolUnipotentDegrees(W)
  p = map((x->begin Position(ud, ud[x]) end), f[:charNumbers])
  pp = map(x->x[1], Filtered(tally(p), (x->begin x[2] > 1 end)))
  res = []
  n = length(f[:charNumbers])
  for d in pp
    q = filter(i->p[i]==d,1:n)
    print(q, "==", CharNames(uc)[f[:charNumbers][q]], "\n")
    push!(res, q)
  end
  Group(Concatenation(map(q->
     map(i->(Perm(q[i],q[i+1]))(n+q[i],n+q[i+1]),1:length(q)-1),res)),Perm())
end

# List of representatives of each class of permutations-with-signs for
# xi-Ennola of fno-th family f of W  under ambiguity group of f.
# the representative in the orbit which occurs in CHEVIE data is taken
# to be compatible with CHEVIE data.
# The third argument is Ennola(W)
# [detfam is to be called with W.ennola unbound if testing nonoccuring orbit]
function repsEnnola(W, fno, e)
  p = EnnolaBete(W, fno, 1)
  n = length(p)
  p = DistinctCartesian(p)
  p = Filtered(p, (x->begin length(Intersection(x, -x)) == 0 end))
  p = Orbits(ambig(W, fno), map(SPerm, p))
  orbits = map((o->begin map((ps->begin Permuted(1:n, ps) end), o) end), p)
  p = map((x->begin (x[:ls])[fno] end), e)
  n = gapSet(map(ps->PositionProperty(orbits, (y->begin ps in y end)), p))
  if length(n) != 1 error("theory") end
  print("All possibilities occuring in W.ennola fall in orbit ",n[1],"\n")
  res = map((x->begin x[1] end), orbits)
  for n in p res[PositionProperty(orbits, (y->(n in y;)))] = n end
  return res
end

# set the fno-th family of W to f
# [at least fourierMat and eigenvalues should have been filled in f]
function tryfam(W, fno, f)
  uc = UnipotentCharacters(W)
  uc.families[fno][:fourierMat] = f[:fourierMat]
  uc.families[fno][:eigenvalues] = f[:eigenvalues]
  delete!(uc, :eigenvalues)
  for i = 1:length(((uc[:families])[fno])[:charNumbers])
      h = PositionProperty(uc[:harishChandra], (h->begin
         (((uc[:families])[fno])[:charNumbers])[i] in h[:charNumbers]end))
      ((uc[:harishChandra])[h])[:eigenvalue] = (f[:eigenvalues])[i]
  end
  return W
end

function qennola(W,i)
  poss=EnnolaBE(W,i,1)
  f=UnipotentCharacters(W).families[i]
  A=fusion_algebra(f)
  b=basis(A)
  res=vcat(map(function(x)
    local p
    p=SPerm(b,Ref(x).*b)
    if !isnothing(p) return [p, SPerm(-vec(p))]
    else return SPerm[]
    end
  end, b)...)
  filter(x->all(j->j^x in poss[j],1:length(f)),res)
end

function qqennola(W)
  gets(W,:ennola)do
    uc=UnipotentCharacters(W)
    v=map(i->qennola(W,i),1:length(uc.families))
    map(s->SPerm(famEn2En(W,s)),cartesian(v...))
  end
end
