function CuspidalPairs(W,d,ad)
  WF=(W isa Spets) ? W : spets(W)
  if (d isa Int) && d!=0 d=1//d else d=d//1 end
  vcat(map(HF->map(char->[HF,char],cuspidal_unipotent_characters(HF, d)),
            split_levis(WF, d, ad))...)
end

function CuspidalPairs(W,d=0)
  WF=(W isa Spets) ? W : spets(W)
  if (d isa Int) && d!=0 d=1//d else d=d//1 end
  vcat(map(ad->CuspidalPairs(WF,d,ad),0:length(relative_degrees(WF,d)))...)
end

function PermutationOnClasses(W, aut)
  PermList(map(c->position_class(W,Representative(c)^aut),ConjugacyClasses(W)))
end

function PermutationOnCharacters(W,aut,chars=1:length(classreps(W)))
  p=PermutationOnClasses(W, aut)
  ct=CharTable(W).irr[chars,:]
  PermListList(ct, map(r->Permuted(r, p), ct))
end

function PermutationOnUnipotents(W,aut,l=1:length(UnipotentCharacters(W)))
  uc = UnipotentCharacters(W)
  t = map(x->x[l], DeligneLusztigCharacterTable(W))
  push!(t, Eigenvalues(uc, l))
  t = TransposedMat(t)
  if length(gapSet(t)) < length(t)
    t = map(x->Position(((uc[:harishChandra])[1])[:charNumbers], x), l)
    if all(x->x!=false,t) return PermutationOnCharacters(W, aut, t)
    else error("Rw + eigen cannot disambiguate\n")
    end
  end
  PermListList(t, map(r->Permuted(r, PermutationOnClasses(W, aut)), t))
end

function FactorsSet(s)
  if s == [[]] return [] end
  for i = s[1]
    s1 = filter(y->i in y,s)
    if length(s1) == length(s) // 2
      i = Intersection(s1)
      r = Difference(s[1], i)
      s2 = filter(y->length(Intersection(y, i))==0,s)
      a = PositionProperty(s2, y->IsSubset(y, r))
      if length(s2) == length(s) // 2 && a != false
        j = Difference(s2[a], r)
        s1 = gapSet(map(x->Difference(x, i), s1))
        s2 = gapSet(map(x->Difference(x, j), s2))
        if length(i)==1 i=i[1] end
        if length(j)==1 j=j[1] end
        if s1 == s2 return Concatenation([[i, j]], FactorsSet(s1))
        end
      end
    end
  end
  [s]
end

LIMSubsetsSum = 10
SubsetsSum = function (S, l, v, lv)
  function sievev(good, v)
    local i, p
    for i in good
      p=findfirst(==(lv[i]),v)
      if isnothing(p) return empty(v) end
      v=v[filter(i->i!=p,1:length(v))]
    end
    v
  end
  c = 0
  found = 0
  inner = function (S, s, t, nonsolved, v, factor)
    local bad, good, p, sols, res, i, sol, f, ll, solved
    c+=1
    if mod(c, 1000) == 0
      InfoChevie("# ",factor,": xcols:",length(nonsolved)-1," xrows:",length(t),
                 " found:", found, "\n")
    end
    t = Filtered(t, (i->begin lv[i] in v end))
    if length(t) == 0
      if iszero(S) found+=1
          return [[]]
      else return []
      end
    end
    ll=map(i->Dict{Symbol,Any}(:pos=>i,:cand=>Filtered(t,j->l[j][i]!=0)),
           nonsolved)
    SortBy(ll, (x->begin length(x[:cand]) end))
    if length((ll[1])[:cand]) > LIMSubsetsSum ll = [ll[1]]
    else ll = Filtered(ll, (x->begin length(x[:cand]) <= LIMSubsetsSum end))
    end
    solved = []
    good = []
    bad = []
    for p = ll
        p[:sols] = Filtered(combinations(p[:cand]), (e->begin
                        Sum((l[e])[p[:pos]]) == S[p[:pos]] end))
        if length(p[:sols]) == 0 return []
        elseif length(p[:sols]) == 1 push!(solved, p[:pos])
        end
        good = Union(good, Intersection(p[:sols]))
        bad = Union(bad, Difference(p[:cand], Union(p[:sols])))
    end
    nonsolved = Difference(nonsolved, solved)
    if length(good) + length(bad) > 0
      return map((r->begin Concatenation(good, r) end), 
       inner(S - Sum(l[good]), Union(s, good), 
         Difference(t, Union(good, bad)), nonsolved, sievev(good, v), factor))
    else
      res = []
      p = maximum(map((x->begin length(x[:cand]) end), ll))
      p = First(ll, (x->begin length(x[:cand]) == p end))
      f = length(p[:sols])
      nonsolved = Difference(nonsolved, [p[:pos]])
      InfoChevie("# ", factor, ": xcols:", length(nonsolved), " xrows:", 
         length(t)," in comb(",length(p[:cand]),")==>", length(p[:sols]), "\n")
      for sol = p[:sols]
          good = sol
          bad = Difference(p[:cand], sol)
          if Intersection(good, bad) == []
              res = Append(res, map((r->begin Concatenation(good, r)
                              end), inner(S - Sum(l[good]), Union(s, good), Difference(t, Union(good, bad)), nonsolved, sievev(good, v), SPrint(factor, ":", f))))
          end
          f-=1
      end
      return res
    end
  end
  return inner(S, [], 1:length(l), 1:length(S), v, "")
end

FitParameter = function (sch, m, d)
  sch = map(x->CycPolOps[:EnnolaTwist](x, E(denominator(d), numerator(d))), sch)
  positive = (p->begin all((x->begin x[2] > 0 end), p[:vcyc]) end)
  den = Lcm(map(denominator, m))
  e = length(m)
  a = tally(m)
  sch = map((x->begin (CycPolOps[:DescentOfScalars])(x, den) end), sch)
  poss = function (i,)
    local v, k, avail, good, p, term
    v = sch[i]
    avail = 1:length(a)
    p = map((k->begin 0:e - 1 end), avail)
    good = []
    term(j,k)=CycPol(1-E(e,j)*X(Cyclotomics)^(den*(m[i]-(a[k])[1])))
    while true
      for k = avail
          p[k] = Filtered(Difference(p[k], Concatenation(p[good])), (j->
                          m[i]==a[k][1] || positive(v // term(j, k))))
      end
      good = Filtered(avail, (i->begin length(p[i]) == (a[i])[2] end))
      if length(good) == 0
          return p
      end
      avail = Difference(avail, good)
      v//=prod(k->prod(j->term(j,k),p[k]),filter(k->m[i]!=a[k][1],good))
    end
  end
  pp = map(poss, 1:e)
  p = PositionProperty(1:e, (i->map(x->x[2],a) == map(length, pp[i]) ))
  if p == false return false end
  v = map(i->map(x->gapSet(map(y->mod(y, e), x)), i + pp[p]), 0:e - 1)
  v = map((k->begin Filtered(1:e, (i->begin all((j->begin
                                          IsSubset((pp[k])[j], (v[i])[j])
                                      end), 1:length(a)) end)) end), 1:e)
  G = map((x->begin Position(sch, x) end), sch)
  v = map((x->begin v[x] end), gapSet(G))
  G = map((x->begin Position(gapSet(G), x) end), G)
  bad = 1:length(v)
  while true
    good = Filtered(bad, (i->length(v[i]) == count(j->G[j]==i, 1:e)))
    bad = Difference(bad, good)
    for p in bad for a in good v[p] = Difference(v[p], v[a]) end end
    if length(good) == 0 break end
  end
  if Sum(v, length) > e
    ChevieErr("non-unique solution\n")
    return false
  end
  res = map((x->begin [] end), 1:length(v))
  for p = 1:e push!(res[G[p]], p) end
  pp = []
  for p = 1:length(res) pp[res[p]] = v[p] end
  para = map((i->begin E(e, pp[i]) * X(Cyclotomics) ^ (den * m[i]) end), 1:e)
  if map(k->prod(j->CycPol(1-para[k]//para[j]),filter(j->j!=k,1:e)),1:e)!=sch
    ChevieErr("schur elms don't match\n")
    return false
  end
  TransposedMat([CollectBy(1:length(G), G), v])
end

struct Series
  spets
  levi
  cuspidal::Int
  d::Root1
  principal::Bool
  classno::Int
  element
  projector
  deg
  prop::Dict{Symbol,Any}
end

function Series(WF, levi, cuspidal, d)
  if d isa Int && d>0 d=1//d end
  d=Root1(;r=d)
  if !(WF isa Spets) WF=spets(WF) end
  principal=UnipotentCharacters(levi).prop[:a][cuspidal]==0 && 
            UnipotentCharacters(levi).prop[:A][cuspidal]==0
  eig=union(map(x->x isa Int ? prime_residues(x).//x : [x], 
                   regular_eigenvalues(levi))...)
  e=minimum(conductor.(eig))
  e=minimum(filter(x->conductor(x)==e,eig))
  eig=refleigen(levi)
  q=maximum(map(x->count(y->y==e,x), eig))
  classno=findall(x->count(y->y==e,x)==q,eig)
  if length(classno)>1 error("classno==",classno) end
  element=levi(classinfo(levi)[:classtext][classno[1]]...)
  q=minimum(map(x->count(==(d),x),refleigen(levi)))
  q=findall(x->count(==(d),x)==q,refleigen(levi))
  projector=eigenspace_projector(WF, classreps(levi)[q[1]], d)
  q=Pol([1],1)
  deg=conj(generic_sign(WF))*generic_order(WF,q) // 
      (conj(generic_sign(levi))*generic_order(levi,q))
  deg=CycPol(shift(deg,-deg.v))*Uch.CycPolUnipotentDegrees(levi)[cuspidal]
  Series(WF,levi,cuspidal,d,principal,classno[1],element,projector,
         deg,Dict{Symbol,Any}())
end

CuspidalSeries(W,d=1)=map(x->Series(W,x[1],x[2],d),CuspidalPairs(W,d))

SeriesOps=Dict{Symbol,Any}()
SeriesOps[:String] = (s->begin Format(s, Dict{Symbol, Any}()) end)
SeriesOps[:Print] = function (s,) print(string(s)) end
SeriesOps[:Display] = function (s, opt)
  opt[:screenColumns] = SizeScreen()[1]
  print(Format(s, opt))
end

function Base.show(io::IO,s::Series)
  TeX=get(io,:TeX,false)
  repl=get(io,:limit,false)
  cname = charnames(io,UnipotentCharacters(s.levi))[s.cuspidal]
  n=repl || TeX ? "\\lambda" : "c"
  if s.spets == s.levi
    print(io,s.d, "-cuspidal ",cname," of ", reflection_name(io,s.spets))
  else
    print(io,s.d, "-series ")
    Util.printTeX(io,"R^{",reflection_name(io,s.spets),"}_","{",
                           reflection_name(io,s.levi), "}(")
    Util.printTeX(io,"$n==",cname,")")
  end
  if haskey(s.prop, :e)
    if TeX Util.printTeX(io,"\\quad W_G(L,\\lambda)==Z_{",e(s), "}")
    else Util.printTeX(io," W_G(L,", n, ")==Z", e(s))
    end
  elseif haskey(s.prop, :WGL)
      if TeX Util.printTeX(io,"\\quad ") end
      if haskey(WGL(s).prop, :reflections)
        Util.printTeX(io," W_G(L,", n, ")==", reflection_name(io,WGL(s)))
      else
        Util.printTeX(io,string(" |W_G(L,", n, ")|==", length(WGL(s))))
      end
  elseif haskey(s.prop, :WGLdims)
      if TeX Util.printTeX(io,"\\quad ") end
      Util.printTeX(io," |W_G(L,", n, ")|==", Sum(s[:WGLdims], (x->x^2)))
  end
  if haskey(s.prop, :translation)
      if TeX Util.printTeX(io,"\\quad ") end
      Util.printTeX(io," translation==", s[:translation])
  end
  if !haskey(io, :displaysize) || !haskey(s.prop,:charNumbers)
    return 
  end
  uw = UnipotentCharacters(s.spets)
  e = length(charNumbers(s))
  function f(arg...,)
    local vals
    if TeX
      push!(rowlab, arg[1])
      vals = arg[length(arg)]
    else
      push!(rowlab, arg[2])
      vals = arg[3]
    end
    push!(m, map(x->Format(x, opt), vals))
  end
  m = []
  rowlab = []
  f("\\hbox{Character}", "Name", CharNames(uw, opt)[charNumbers(s)])
  f("\\hbox{specializes}", "specializes", CharNames(WGL(s), opt))
  f("\\varepsilon", "eps", s[:eps])
  if haskey(s, :eigen) f("\\hbox{eigen}", "eigen", s[:eigen]) end
  if haskey(s, :e)
    if haskey(s, :Hecke)
      f("\\hbox{parameter}", "param", s[:Hecke][:parameter][1])
    elseif haskey(s, :mC)
     f("\\hbox{mC}", "degparam", mC(s))
    end
  end
  f("\\hbox{family \\#}", "family #", map(j->
                                          findfirst(k->j in k[:charNumbers],uw[:families]), charNumbers(s)))
  if haskey(s, :permutable) && any((x->begin x != false end), s[:permutable])
      f("\\hbox{permutable}", "permutable", s[:permutable])
  end
  m = TransposedMat(m)
  opt[:rowLabels] = charNumbers(s)
  opt[:columnLabels] = rowlab
  res*=SPrint("\n", FormatTable(m, opt))
end

function mC(s::Series)
  get(s,:mC) do
  e = HyperplaneOrbits(WGL(s))
  if length(e)>1 return else e = e[1][:e_s] end
  uc = UnipotentCharacters(s.spets)
  cn = charNumbers(s)[Filtered(1:length(s[:dims]), (i->s[:dims][i]==1))]
  aA = uc[:a][cn] + uc[:A][cn]
  lpi(W)=sum(degrees(W)+codegrees(W))
  if s.principal
    if Minimum(aA) != 0 error("id not in RLG(1)") end
    pG = lpi(Group(s.spets))
    pL = lpi(Group(s.levi))
    D0 = pG - pL
    xiL = AsRootOfUnity(PhiOnDiscriminant(s.levi)) * denominator(s.d)
    xiG = AsRootOfUnity(PhiOnDiscriminant(s.spets)) * denominator(s.d)
    if xiL != xiG
        ChevieErr("fixing dimension of variety by xiL-xiG==", xiL - xiG, "\n")
        D0 = (D0 + xiL) - xiG
    end
    if !IsInt(D0 * s.d)
        ChevieErr(s, "\n    ==> (l(pi_G)==", pG, ")-(l(pi_L)==", pL,
                  ")+(xiL*d==", xiL, ")-(xiG*d==", xiG, ")  !  ==0( %
           d==", s.d, ") Max(aA)-Min(aA)==", maximum(aA) - Minimum(aA), "\n")
    end
    if !IsInt(D0//e)
      ChevieErr("fixing dimension of variety to be divisible by e==", e, "\n")
      D0-=mod(D0, e)
    end
    (D0 - aA)//e
  else
    D0=maximum(aA)-Minimum(aA)
    (D0+Minimum(aA)-aA)//e
  end
  end
end

function sum_intersection(m::Matrix,n::Matrix)
  mat=[m m;n zero(n)]
  mat=echelon(mat)[1]
  sm=mat[:,axes(m,2)]
  in=mat[:,size(m,2)+axes(m,2)]
  sm=sm[1:count(!iszero,eachrow(sm)),:]
  in=in[size(sm,1)+1:end,:]
  in=in[1:count(!iszero,eachrow(in)),:]
  (sm,in)
end

if false
function RelativeGroup(s::Series)
  gets(s,:WGL) do
  W = Group(s.spets)
  refs = gapSet(map(x->Position(W[:reflections], x), Reflections(W)))
  mats = []
  mats[refs] = map(r->reflrep(W, Reflection(W, r)), refs)
  id = s.projector ^ 0
  if s.projector == id
   WGL(s) = deepcopy(W)
   (WGL(s))[:parentMap] = (WGL(s))[:generators]
   (WGL(s))[:reflists] = map((x->begin
                      [x]
                  end), (W[:rootInclusion])[W[:generatingReflections]])
  else
      V = NullspaceMat(s.projector - id)
      rH = map((r->begin
                      (SumIntersectionMat(NullspaceMat(mats[r] - id), V))[2]
                  end), refs)
      rH = gapSet(Filtered(rH, (H->begin
                          length(H) < length(V)
                      end)))
      rel = []
      for H = rH
          if length(H) == 0
              ConjugacyClasses(Group(s.levi))
              H = W
          else
              H = Filtered(refs, (r->begin
                              H * mats[r] == H
                          end))
              H = ReflectionSubgroup(W, (W[:rootInclusion])[H])
          end
          push!(rel, Dict{Symbol, Any}(:refs =>
                                       (H[:rootInclusion])[H[:generatingReflections]],
                                       :WH => Normalizer(H,
                                                         Group(s.levi))))
      end
      if (s.spets).phi == Perm()
          WF = W
          L = Group(s.levi)
      else
          WF = Group(Concatenation(W[:generators], [(s.spets).phi]), W[:identity])
          L = Subgroup(WF, (Group(s.levi))[:generators])
      end
      if length(L) == 1
       WGL(s) = Centralizer(W, (s.levi).phi)
          for H = rel
              H[:WH] = Centralizer(H[:WH], (s.levi).phi)
          end
          func = (x->begin
                      x
                  end)
      else
          NLF = Normalizer(WF, L)
          if NLF == L
              if s.levi != s.spets
                  error(s, " N==L\n")
                  return false
              end
              WGL(s) = CoxeterGroup()
              (WGL(s))[:reflists] = []
              s[:WGLdims] = [1]
              return WGL(s)
          end
          NLFb = FactorGroup(NLF, L)
          hom = NaturalHomomorphism(NLF, NLFb)
          phi = Image(hom, (s.levi).phi)
          s[:WGLdegrees] = RelativeDegrees(s.spets, s.d)
          if count((x->begin
                                  x != 1
                              end), RelativeDegrees(s.levi, s.d)) == 0 && count((x->begin
                                  x != 1
                              end), s[:WGLdegrees]) == 1
              if (s.levi).phi in W && Order(NLFb, phi) == Product(s[:WGLdegrees])
               WGL(s) = Subgroup(NLFb, [phi])
              else
                  N = Centralizer(Subgroup(WF, W[:generators]),(s.levi).phi)
                  if length(N) == length(L) * Product(s[:WGLdegrees])
                   WGL(s) = FactorGroup(N, L)
                  end
              end
          end
          if !(haskey(s, :WGL))
              N = Subgroup(WF, W[:generators])
              N = Subgroup(N, (Intersection(N, NLF))[:generators])
              WGL(s) = Centralizer(FactorGroup(N, L), phi)
          end
          if length(rel) == 1 && (rel[1])[:WH] == NLF
           (rel[1])[:WH] = WGL(s)
          else
              rel = map(function (x,)
                          local N, WH
                          N = Intersection(WGL(s), FactorGroup(Subgroup(NLF, (x[:WH])[:generators]), L))
                          if IsGroup(N)
                              x[:WH] = N
                          else
                           x[:WH] = Subgroup(WGL(s), N)
                          end
                          return x
                      end, rel)
          end
          func = (x->begin
                      (x[:element])[:representative]
                  end)
      end
      ud = CycPolUnipotentDegrees(s.levi)
      eig = Eigenvalues(UnipotentCharacters(s.levi))
      ud = Filtered(1:length(ud), (i->begin
                      ud[i] == ud[s.cuspidal] && eig[i] ==
                      eig[s.cuspidal]
                  end))
      if length(ud) > 1
       c = length(WGL(s))
       WGL(s) = Stabilizer(WGL(s), Position(ud, s.cuspidal), function (c, g)
                      return c ^ PermutationOnUnipotents(s.levi, func(g), ud)
                  end)
       if c != length(WGL(s))
        ChevieErr("# WGL:", c, "/", length(WGL(s)), " fix c\n")
          end
          for H = rel
           H[:WH] = Intersection(H[:WH], WGL(s))
          end
          rel = Filtered(rel, (x->begin
                          length(x[:WH]) > 1
                      end))
      end
      if !(all((x->begin
                          IsCyclic(x[:WH])
                      end), rel))
          error("a")
      end
      hom = map((x->begin
                      First(elements(x[:WH]), (y->begin
                                  Order(x[:WH], y) == length(x[:WH])
                              end))
                  end), rel)
      if Subgroup(WGL(s), hom) != WGL(s)
          error("b")
      end
      reflists = map((x->begin
                      x[:refs]
                  end), rel)
      m = Concatenation(V, NullspaceMat(s.projector)) ^ -1
      rel = map(function (x,)
                  local M
                  M = reflrep(W, func(x)) ^ m
                  return (M[1:length(V)])[1:length(V)]
              end, hom)
      rel = map(AsReflection, rel)
      for c = 1:length(rel)
          r = AsRootOfUnity((rel[c])[:eigenvalue])
          m = numerator(r)
          r = denominator(r)
          if m != 1
              rel[c] = AsReflection(Reflection((rel[c])[:root], (rel[c])[:coroot]) ^ mod(1 // m, r))
          end
      end
      WGL(s) = PermRootGroup(map((x->begin
                          x[:root]
                      end), rel), map((x->begin
                          x[:coroot]
                      end), rel))
      (WGL(s))[:parentMap] = hom
      (WGL(s))[:reflists] = reflists[map((x->(PositionProperty(rel,
                                                               (y->(ProportionalityCoefficient(y[:root],
                                                                                               x)
                                                                    !=
                                                                    false;)));)),
                                         ((WGL(s))[:roots])[(WGL(s))[:generatingReflections]])]
  end
  s[:WGLdims] = map((x->begin
                  x[1]
                 end), (CharTable(WGL(s)))[:irreducibles])
  if all((x->begin
                  x == 1
              end), s[:WGLdims])
      s[:e] = length(s[:WGLdims])
  end
  WGL(s)
end
end
else
function RelativeGroup(s::Series)
  gets(s,:WGL) do
  W=Group(s.spets)
  L=Group(s.levi)
  if isone(s.projector)
    WGL=deepcopy(W)
    WGL.prop[:parentMap] = gens(WGL)
    WGL.prop[:reflists] = map(x->[x],inclusiongens(W))
    s.prop[:WGLdims] = CharTable(WGL).irr[:,1]
    if is_cyclci(WGL) s.prop[:e] = length(WGL) end
    return WGL
  end
  N=normalizer(W.G, L.G)
  if !isone(s.levi.phi)
    if length(L) == 1
      N=centralizer(N, s.levi.phi)
    elseif s.spets isa CoxeterCoset
      N=Group(map(x->reduced(Group(s.levi), x),gens(N)))
      N=centralizer(N, s.levi.phi)
      N=Subgroup(W, vcat(gens(N),gens(Group(s.levi))))
    else
      NF=centralizer(N, s.levi.phi)
      if length(NF)==length(L)*prod(relative_degrees(s.spets, s.d)) N=NF
      else
        NF=Group(vcat(gens(N),[s.levi.phi]))
        LF = Subgroup(NF, gens(L))
        NFQ = NF // LF
        N=Subgroup(NFQ,map(x->x^NaturalHomomorphism(NF,NFQ),gens(N)))
        NF = centralizer(NFQ, s.levi.phi^NaturalHomomorphism(NF, NFQ))
        NFQ = intersection(NF, N)
        N=Subgroup(W,vcat(gens(L),map(x->x.phi,gens(NFQ))))
      end
    end
  end
  eig=Uch.eigen(UnipotentCharacters(s.levi))
  eig=filter(i->eig[i]==eig[s.cuspidal],eachindex(eig))
  if length(eig)>1
    ud=Uch.CycPolUnipotentDegrees(s.levi)
    ud=filter(i->ud[i]==ud[s.cuspidal],eachindex(eig))
    if length(ud)>1
      c=length(N)
      N=centralizer(N,findfist(==(s.cuspidal),ud), (c, g)->
                  c^PermutationOnUnipotents(s.levi, g, ud))
      if c!=length(N)
        ChevieErr("# WGL:",length(N)//length(L),"/",c//length(L)," fix c\n")
      end
    end
  end
  refs=unique(sort(map(x->findfirst(==(x),reflections(W)), reflections(W))))
  if length(L) == 1
    WGL=N
    func=x->x
  elseif length(N)==length(L)
    if s.levi!=s.spets
      error(s, " N==L\n")
      return nothing
    end
    WGL=coxgroup()
    WGL.prop[:reflists]=[]
    s.prop[:WGLdims]=[1]
    return WGL
  else
    WGL=N/L
    func = x->x.phi
  end
  V = GLinearAlgebra.lnullspace(s.projector-one(s.projector))
  m = vcat(V, GLinearAlgebra.lnullspace(s.projector))^-1
  hplane(r)=sum_intersection(GLinearAlgebra.lnullspace(reflrep(W,
                     reflection(W,r))-one(s.projector)),V)[2]
  smalltobig(h)=map(x->vcat(x, fill(0, max(0, rank(W)-length(V)))),h)*m^-1
  function bigtosmall(x)
    x=x^m
    x[1:length(V),1:length(V)]
  end
  function getreflection(rr)
    local rH, H, r, res, n
    rH = hplane(rr[1])
    if length(rH)==length(V) return false end
    if length(rH)==0
#     ConjugacyClasses(L)
      H=W
    else
      H=reflection_subgroup(W, vcat(inclusion(W,rr), inclusiongens(L)))
    end
    res=Dict{Symbol, Any}(:refs => inclusiongens(H))
    println("H=$H N=$N")
    H = intersection(H, N)
    if length(L)!=1 H=H//L end
    if length(H)==1 return false end
    if !(is_cyclic(H)) error("a") end
    res[:hom] = First(elements(H), (y->order(H, y)==length(H)))
    r = bigtosmall(reflrep(W, func(res[:hom])))
    Inherit(res, AsReflection(r))
    n = res[:eig]
    n = mod(1 // numerator(n), denominator(n))
    r = r ^ n
    Inherit(res, AsReflection(r))
    res[:WH] = H
    return res
  end
  function v(h)
    if length(h) == 0 return VectorSpace([[0]], Cyclotomics)
    else return VectorSpace(h, Cyclotomics)
    end
  end
  function v(h)
    if size(h,1)==0 return  h end
    h=GLinearAlgebra.echelon(h)[1]
    h[1:count(!iszero,eachrow(h)),:]
  end
  rrefs = collect(values(groupby(x->v(hplane(x)),refs)))
  rrefs=filter(x->!(inclusion(W,x[1]) in inclusion(L)),rrefs)
  reflist = []
  for r in rrefs
    push!(reflist, getreflection(r))
    if Subgroup(WGL, map(x->x[:hom],reflist)) == WGL
      WGL=PermRootGroup(map(x->x[:root],reflist),map(x->x[:coroot],reflist))
      reflist=map(x->smalltobig(NullspaceMat(x-x^0)),reflrep(WGL))
      reflist=map(h->First(rrefs, (rr->v(hplane(rr[1]))==v(h))), reflist)
      WGL.prop[:reflists] = map(getreflection, reflist)
      WGL.prop[:parentMap] = map(x->x[:hom], WGL.prop[:reflists])
      WGL.prop[:reflists] = map(x->x[:refs], WGL.prop[:reflists])
      s.prop[:WGLdims] = CharTable(WGL).irr[:,1]
      if IsCyclic(WGL(s)) s[:e]=length(WGL(s)) end
      return WGL
    end
  end
  error("b")
end
end
end

WGL(s)=RelativeGroup(s)

function CharNumbers(s::Series)
  get(s,:charNumbers) do
  if !haskey(s.prop,:WGL) && isnothing(RelativeGroup(s)) return nothing end
  ud = CycPolUnipotentDegrees(s.spets)
  ad = count(!isone,relative_degrees(s.levi, s.d))
  cand = filter(i->ad==valuation(ud[i],s.d),1:length(ud))
  s[:RLG]=LusztigInduction(s.spets, UnipotentCharacter(s.levi, s.cuspidal))
  if s[:RLG] == false && isone(s.levi.phi)
    s[:RLG] = HarishChandraInduction(s.spets, UnipotentCharacter(s.levi,
                                                                   s.cuspidal))
  end
  if s[:RLG] == false ChevieErr(s, ":RLG failed\n")
  elseif s.degree!=CycPol(degree(s[:RLG]))
    ChevieErr(s,":Deg RLG!=Sum(eps[i]*ud[i])\n")
  end
  cand = map(c->Dict(:charNumbers=>c,:sch=>s.degree//ud[c]),cand)
  cand = Filtered(cand, (c->begin all(x->x[2]>0, c[:sch][:vcyc]) end))
  Ed = E(denominator(s.d), numerator(s.d))
  q = Mvp("q")
  ad = CycPol(q - Ed) ^ ad
  f = s.degree // ad
  if any(x->x[2]<0, f[:vcyc])
      ChevieErr(s, " cuspidal is  ! \n")
      return false
  end
  v = Value(f, Ed)
  for c in cand
    c[:span] = degree(c[:sch])-c[:sch][:valuation]
    f = (Sum(s[:WGLdims],x->x^2)*Value(ud[c[:charNumbers]] // ad, Ed)) // v
    c[:dims] = abs(f)
    c[:eps] = sign(f)
  end
  cand = Filtered(cand, (c->begin c[:dims] in s[:WGLdims] end))
  eig = Eigenvalues(UnipotentCharacters(s.levi), [s.cuspidal])[1]
  eig*= map(i->E(conductor(s.d)^2,i),1:conductor(s.d)^2)
  cand = Filtered(cand, c->
    Eigenvalues(UnipotentCharacters(s.spets), [c[:charNumbers]])[1] in eig )
  if length(cand) < length(s[:WGLdims])
      ChevieErr(s, ": not enough left with predicted eigenvalues in ", FormatGAP(map(AsRootOfUnity, eig)), "\n")
      return false
  end
  sort!(cand,by=x->x[:dims])
  sort!(s[:WGLdims])
  f = function (arg...,)
    if length(arg) == 1 return map(x->x[Symbol(arg[1])], cand) end
    if IsList(arg[2]) return map(x->x[Symbol(arg[1])], cand[arg[2]])
    else return cand[arg[2]][Symbol(arg[1])]
    end
  end
  check = function ()
    local n
    for n = ["charNumbers", "eps", "dims", "span"] s[Symbol(n)] = f(n) end
    if s[:RLG]!=false && (filter(i->s[:RLG][:v][i]!=0,1:length(s[:RLG][:v]))!=
                          gapSet(charNumbers(s)) ||
                          s[:RLG][:v][charNumbers(s)]!=s[:dims].*s[:eps])
      ChevieErr(s, ":RLG does not match")
    end
    return charNumbers(s)
  end
  if length(cand) == length(s[:WGLdims]) return check() end
  ud=map(f("dims"), CycPolUnipotentDegrees(s.spets)[f("charNumbers")], f("eps"))
  t = maximum(map(degree, ud))
  c=p->Concatenation(Coefficients(Value(p,q),"q"),fill(0,max(0,t-degree(p))))
  v = SubsetsSum(c(s.degree), map(c, ud), s[:WGLdims], f("dims"))
  InfoChevie("# ", length(v), " found\n")
  if length(v) > 10000
      InfoChevie("# ", length(v), " combinations sum to dimRLG\n")
  elseif length(v) == 0
      ChevieErr(s, " no combination sums to dimRLG\n")
      return false
  end
  if haskey(s, :e)
    charNumbers(s) = f("charNumbers")
    s[:dims] = f("dims")
    mC(s)
    if length(v) > 1
      v = Filtered(v, (a->begin f("span", a) == map(i->
                                                   Sum(Difference(a, [i]),
                                                      (j->abs((mC(s))[i] - (mC(s))[j]))), a) end))
      if length(v) > 10000
          InfoChevie("# ", length(v), " combinations have right span\n")
      end
    elseif length(v) == 0
        ChevieErr(s, " no combination has right span\n")
        return false
    end
    delete!(s, :charNumbers)
    delete!(s, :dims)
    delete!(s, :mC)
  end
  if length(v) > 1
    InfoChevie("# after span ", length(v), " combinations\n")
    InfoChevie("# Warning: using Mackey with tori for ", s, "\n")
    i = FusionConjugacyClasses(s.levi, s.spets)
    c=map((a,b)->a//b, CharTable(s.spets)[:centralizers][i], CharTable(s.levi)[:centralizers])
    t =
    map((a,b)->a*b,TransposedMat(DeligneLusztigCharacterTable(s.levi))[s.cuspidal], c)
    t = map((k->begin
                    Sum(t[Filtered(1:length(i), (j->(i[j] == k;)))])
                end), 1:length(ConjugacyClasses(s.spets)))
    c = TransposedMat(DeligneLusztigCharacterTable(s.spets))
    v = Filtered(v, (a->begin map(function (x, y) return x * y
         end, f("dims", a), f("eps", a)) * c[f("charNumbers", a)] == t end))
  end
  if length(v) > 1
      ChevieErr(s, " ", Join(map((x->begin FormatGAP(x)
                      end), FactorsSet(map((x->begin f("charNumbers", x)
                             end), v))), "x"), " chars candidates: using RLG\n")
      if s[:RLG] == false return false end
      v = Filtered(v, (l->begin all((i->begin
                  ((s[:RLG])[:v])[f("charNumbers", i)] != 0 end), l) end))
      if length(v) != 1 error() end
  elseif length(v) == 0
      ChevieErr(s, " no candidates left\n")
      return false
  end
  cand = cand[v[1]]
  return check()
  end
end

COMPACTCOHOMOLOGY = true
SeriesOps[:fill] = function (s,)
        local uc, Schur, r, p, ratio, LFrob, predeigen, map, series, noncus, unique, i, m, quality, rr, FractionToRoot, a, param, j, o, u, nid
        if !(haskey(s, :e))
            error("fill assumes .e bound\n")
        end
        if !(haskey(s, :charNumbers)) && CharNumbers(s) == false
            return false
        end
        uc = UnipotentCharacters(s.spets)
        Schur = (CycPolUnipotentDegrees(s.spets))[charNumbers(s)]
        Schur = map((x->begin s.degree // Schur[x] * (s[:eps])[x] end), 1:s[:e])
        s[:eigen] = Eigenvalues(uc, charNumbers(s))
        LFrob =
        AsRootOfUnity((Eigenvalues(UnipotentCharacters(s.levi)))[s.cuspidal])
        m = ReflectionDegrees(Group(s.spets))
        s[:delta] = Lcm(map((x->begin
                            denominator(AsRootOfUnity(x[2]))
                        end), Filtered(ReflectionDegrees(s.spets), (x->begin
                                x[1] != 1
                            end))))
        rr = function (j, i) return (i - 1) // s[:e] - (mC(s))[j] * s.d end
        FractionToRoot = (x->begin E(denominator(x), numerator(x)) end)
        param = function (j, i)
         return Mvp("q") ^ (mC(s))[j] * FractionToRoot(rr(j, i))
            end
        predeigen = function (j, i)
                return FractionToRoot(s[:delta] * rr(j, i) * s[:e] * s.d + LFrob)
            end
        if COMPACTCOHOMOLOGY
         map = FitParameter(Schur, mC(s), s.d)
        else
         map = FitParameter(Schur, -(mC(s)), s.d)
        end
        if map == false
            ChevieErr(s, " FitParameter failed\n")
            return false
        end
        nid =
        (((uc[:almostHarishChandra])[1])[:charNumbers])[PositionId(s.spets)]
        if nid in charNumbers(s)
            predeigen = function (j, i)
                    return FractionToRoot(s.d * s[:e] * s[:delta] * (rr(j, i)
                                                                     + s.d *
                                                                     (mC(s))[Position(charNumbers(s), nid)]))
                end
        end
        series = map((x->begin
                        PositionProperty(uc[:harishChandra], (y->begin
                                    x in y[:charNumbers]
                                end))
                       end), charNumbers(s))
        unique = Filtered(map, (p->begin
                        length(p[1]) == 1
                    end))
        ratio = map((p->begin
                        (s[:eigen])[(p[1])[1]] // predeigen((p[1])[1], (p[2])[1])
                    end), unique)
        if length(gapSet(ratio)) > 1
            ChevieErr(s, " eigenvalue ratios==", ratio, "\n")
            return false
        end
        ratio = AsRootOfUnity(ratio[1])
        ratio = ratio * denominator(s.d * s[:delta])
        if !(IsInt(ratio))
            ChevieErr(s, "non-integral ratio==", ratio, "\n")
            return false
        end
        if ratio == 0
            r = 0
        else
            r = QuotientMod(ratio, numerator(s.d * s[:delta]), denominator(s.d * s[:delta]))
        end
        if r == false
            error()
        end
        map = map((x->begin
                        [x[1], map((y->begin
                                        mod(y, s[:e]) + 1
                                    end), (x[2] + r) - 1)]
                    end), map)
        r = []
        s[:permutable] = map((x->begin
                        false
                       end), charNumbers(s))
        j = 1
        for i = map
            a = arrangements(i[2], length(i[2]))
            p = PositionsProperty(a, (A->begin
                            (s[:eigen])[i[1]] == map((j->begin
                                            predeigen((i[1])[1], j)
                                        end), A)
                        end))
            if p == []
                ChevieErr(s, "predicted eigenvalues cannot match actual\n")
                return false
            else
                if length(p) > 1
                    o = Orbits(ApplyFunc(Group, map((x->begin
                                            PermListList(a[p[1]], x)
                                        end), a[p])), 1:length(i[2]))
                    if length(p) != Product(o, (x->begin
                                        factorial(length(x))
                                    end))
                        error()
                    end
                    for u = o
                        (s[:permutable])[(i[1])[u]] = j + u * 0
                        j = j + 1
                    end
                end
                r[i[1]] = a[p[1]]
            end
        end
        p = SortingPerm(r)
        for i = ["mC", "charNumbers", "eigen", "span", "eps", "dims", "permutable"]
            s[Symbol(i)] = Permuted(s[Symbol(i)], p)
        end
        s[:translation] = Filtered(0:s[:e] - 1, (t->begin
                        s[:eigen] == map((i->begin
                                        predeigen(i, mod(i + t, s[:e]))
                                    end), 1:s[:e])
                    end))
        if length(s[:translation]) == 1
            delete!(s, :translation)
        else
            p = (s[:translation])[2]
            if any((x->begin
                                mod(x, p) != 0
                            end), s[:translation]) || length(s[:translation]) != s[:e] // p
                error()
            end
            quality = map((t->begin
                            map((i->begin
                                        param(mod((t + i) - 1, s[:e]) + 1, i)
                                    end), 1:s[:e])
                        end), s[:translation])
            quality = map((x->begin
                            NF(gapSet((Product(x, (y->(Mvp("x") - y;))))[:coeff]))
                        end), quality)
            quality = 1 + ListBlist(s[:translation], map((x->begin
                                    all((j->begin
                                                IsSubset(j, x)
                                            end), quality)
                                end), quality))
            if quality == []
                quality = [1]
            end
            m = Rotations(mC(s))
            m = Position(m, maximum(m[quality]))
            m = circshift(1:s[:e], m - 1)
            for i = ["mC", "charNumbers", "eigen", "span", "eps", "dims", "permutable"]
                s[Symbol(i)] = (s[Symbol(i)])[m]
            end
            s[:translation] = p
        end
        s[:Hecke] = Hecke(WGL(s), [map((i->begin
                                param(i, i)
                            end), 1:s[:e])])
        return s[:Hecke]
    end
SeriesOps[:Hecke] = function (s,)
        local H
        if haskey(s, :Hecke)
            return s[:Hecke]
        end
        CharNumbers(s)
        if haskey(s, :e)
            (SeriesOps[:fill])(s)
        elseif CHEVIE[:relativeSeries]
            InfoChevie("\n      # Relative: ", string(s))
            s[:relativeSeries] = (SeriesOps[:RelativeSeries])(s)
        end
        if haskey(s, :Hecke)
            return s[:Hecke]
        else
            return false
        end
    end
AllProperSeries = function (W,)
        local l
        l = gapSet(map(denominator, Flat(refleigen(W))))
        return Concatenation(map((d->begin
                            Concatenation(map((i->begin
                                            CuspidalSeries(W, d, i)
                                        end), 1:length(RelativeDegrees(W, d))))
                        end), l))
    end
findfractionalpowers = function (W,)
        local uc, h, L, cusp, sers, s, p, fix, i, n, reasons, frac, UNBOUND
        if W[:nbGeneratingReflections] == 0
            return [0]
        end
        uc = UnipotentCharacters(W)
        UNBOUND = 99
        uc[:fractions] = fill(0, max(0, (1 + length(uc)) - 1)) + UNBOUND
        reasons = map((x->begin
                        []
                    end), 1:length(uc))
        for h = uc[:harishChandra]
            if length(h.levi) != W[:nbGeneratingReflections]
                L = ReflectionSubgroup(W, h.levi)
                cusp = FindCuspidalInLevi(h[:cuspidalName], L)
                n = findfractionalpowers(L)
                if n[cusp] != UNBOUND
                    (uc[:fractions])[h[:charNumbers]] = h[:charNumbers] * 0 + n[cusp]
                    for i = h[:charNumbers]
                        push!(reasons[i], joindigits(h.levi))
                    end
                end
            end
        end
        sers = AllProperSeries(W)
        map(Hecke, sers)
        sers = gapSet(Filtered(sers, (x->begin
                            haskey(x, :mC) && haskey(x, :e)
                        end)))
        while length(sers) != 0
            for s = sers
             p = PositionProperty(charNumbers(s), (i->begin
                                uc.fractions[i] !== nothing
                            end))
                if p == false
                    print("reticent series", s, "\n")
                else
                 frac = (mC(s))[p] * s[:e] * s.d
                    fix = Mod1((uc[:fractions])[(charNumbers(s))[p]] - frac)
                    if fix != 0
                        print(" **** Badly normalized series ", s, " adding ", fix, "\n")
                    end
                    print("\n", s, "==>")
                    for i = 1:s[:e]
                     print((charNumbers(s))[i], ".c")
                     frac = Mod1((mC(s))[i] * s[:e] * s.d + fix)
                        if (uc[:fractions])[(charNumbers(s))[i]] == UNBOUND
                         (uc[:fractions])[(charNumbers(s))[i]] = frac
                         push!(reasons[(charNumbers(s))[i]], s.d)
                        elseif (uc[:fractions])[(charNumbers(s))[i]] != frac
                         print("Failed! ", (charNumbers(s))[i], "==",
                               TeXStrip((uc[:TeXCharNames])[(charNumbers(s))[i]]),
                               " in ", s, "\n     where mC mod1==", frac, "\n
                               conflicts with ", reasons[(charNumbers(s))[i]], "\n")
                        else
                         push!(reasons[(charNumbers(s))[i]], s.d)
                        end
                    end
                    sers = Difference(sers, [s])
                end
            end
        end
        for i = uc[:harishChandra]
            if length(i[:charNumbers]) > 1
                print(Format([(uc[:fractions])[i[:charNumbers]], map(length, reasons[i[:charNumbers]])]), "\n")
            else
                print((uc[:fractions])[(i[:charNumbers])[1]], reasons[(i[:charNumbers])[1]], "\n")
            end
            p = gapSet((uc[:fractions])[i[:charNumbers]])
            if haskey(i, :qEigen)
                if p != [i[:qEigen]]
                    if p == [UNBOUND]
                        print("!!!!!!!!!!!! for ", i[:charNumbers], " qEigen should be UNBOUND is ", i[:qEigen], "\n")
                    else
                        error("for HCseries ", i, " of ", ReflectionName(W), " qEigen should be ", p)
                        (uc[:fractions])[i[:charNumbers]] = i[:qEigen] + 0 * i[:charNumbers]
                    end
                end
            else
                if p != [0]
                    print("!!!!!!!!!!!!!! qEigen unbound should be ", p, "\n")
                end
            end
        end
        return uc[:fractions]
    end
PrincipalSeries = function (W, d)
        local s
        s = SplitLevis(W, d, length(RelativeDegrees(W, d)))
        if length(s) != 1
            error(" !  one ", d, "-Sylow")
        end
        s = s[1]
        return Series(W, s, Position((UnipotentCharacters(s))[:A], 0), d)
    end
CheckMaximalPairs = function (W,)
        local divs, l, s, d, c, e, res
        divs = gapSet(map(denominator, Flat(refleigen(W))))
        res = map((d->begin
                        PrincipalSeries(W, d)
                    end), divs)
        for s = res
            e = elements(s.levi)
            any(function (x,)
                    local Z
                    Z = Centralizer(Group(s.spets), x)
                    Z = Group(Concatenation(Z[:generators], (Group(s.levi))[:generators]), Perm())
                    if length(Z) != length(RelativeGroup(s)) * length(s.levi)
                        error()
                    end
                    return true
                end, e)
        end
        return res
    end
ReflectionCosetFromType = function (arg...,)
        local g, t, g1, i
        if length(arg) == 0
            return CoxeterGroup()
        end
        if length(arg) > 1
            error(ReflectionName(t), " non-irred  !  yet")
        end
        for t = arg
            g = ReflectionGroup((t[:orbit])[1])
            if length(t[:orbit]) > 1
                g1 = g
                for i = 1:length(t[:orbit]) - 1
                    g = g * g1
                end
                i = reflrep(g, t[:twist] * PermList(Concatenation((Rotations(map((x->(x[:indices];)), g[:type_])))[2])))
                g = Spets(g, i)
            else
                g = Spets(g, t[:twist])
            end
        end
        return g
    end
SpetsFromType = function (t,)
        local W
        W = ApplyFunc(ReflectionGroup, t[:orbit])
        return Spets(W, t[:twist])
    end
getHecke = function (s,)
        local t, scal, g, l, s1, p, e, c
        t = refltype(s.spets)
        if !(haskey(s, :charNumbers))
            CharNumbers(s)
        end
        if length(t) == 1 && (haskey(t[1], :scalar) && !(all((x->begin
                                    x == 1
                                end), (t[1])[:scalar])))
            t = t[1]
            scal = t[:scalar]
            t = copy(t)
            InfoChevie("   # removing scal==", scal, "\n")
            t[:scalar] = map((x->begin
                            1
                        end), t[:scalar])
            if ((t[:orbit])[1])[:series] == "B" && (((t[:orbit])[1])[:rank] == 2 && t[:twist] == perm"(1,2)")
                ((t[:orbit])[1])[:cartanType] = E(8) - E(8, 3)
            end
            g = SpetsFromType(t)
            if g == false
                return false
            end
            e = WordEnumerator(Group(s.spets))
            p = EltWord(Group(g), (e[:Get])((s.levi).phi // (s.spets).phi))
            Reflections(Group(g))
            l = SubSpets(g, map((x->begin
                                Position((Group(g))[:reflections], EltWord(Group(g), (e[:Get])(x)))
                            end), (Group(s.levi))[:generators]), p)
            if length(scal) > 1
                ChevieErr("scal==", scal, " unimplemented\n")
                return false
            end
            scal = AsRootOfUnity(scal[1])
            if !(s.cuspidal) in CuspidalUnipotentCharacters(l, Mod1(s.d - scal))
                e = (Ennola(Group(l)))[1]
                c = abs(s.cuspidal ^ e[:ls])
                if !c in CuspidalUnipotentCharacters(l, Mod1(s.d - scal))
                    c = abs(s.cuspidal ^ (e[:ls] ^ -1))
                end
            else
                c = s.cuspidal
            end
            s1 = Series(g, l, c, Mod1(s.d - scal))
            p = getHecke(s1)
            if p != false
                p = Value(p, ["q", E(denominator(scal), -(numerator(scal))) * Mvp("q")])
            end
            return p
        else
            if any((u->begin
                            haskey(u, :scalar) && !(all((x->begin
                                                x == 1
                                            end), u[:scalar]))
                        end), t)
                ChevieErr("scals==", map((u->begin
                                u[:scalar]
                            end), t), " unimplemented\n")
            end
            (SeriesOps[:fill])(s)
            if haskey(s, :Hecke)
                return ((s[:Hecke])[:parameter])[1]
            else
                return false
            end
        end
    end
SeriesOps[:RelativeSeries] = function (s,)
        local res, p, ud, u1, f, o, aA, tt
        if !(haskey(s, :charNumbers))
            CharNumbers(s)
        end
        res = map(function (r,)
                    local R, w, i, l
                    i = Position((WGL(s))[:reflists], r)
                    R = SubSpets(s.spets, r, s[:element] // (s.spets).phi)
                    l = (Group(s.levi))[:rootInclusion]
                    if !(IsSubset((Group(R))[:rootInclusion], l))
                        r = map(function (w,)
                                    local p, rr
                                    rr = (Group(s.spets))[:roots]
                                    p = PositionProperty(l, (y->begin
                                                    ProportionalityCoefficient(rr[w], rr[y]) != false
                                                end))
                                    if p != false
                                        return l[p]
                                    end
                                    return w
                                end, r)
                        R = SubSpets(s.spets, r, s[:element] // (s.spets).phi)
                        if !(IsSubset((Group(R))[:rootInclusion], r))
                            l = gapSet(map(function (w,)
                                            local p, rr
                                            rr = (Group(s.spets))[:roots]
                                            p = PositionProperty((Group(R))[:rootInclusion], (y->begin
                                                            ProportionalityCoefficient(rr[w], rr[y]) != false
                                                        end))
                                            if p != false
                                                return ((Group(R))[:rootInclusion])[p]
                                            end
                                            return w
                                        end, l))
                            if !(IsSubset((Group(R))[:rootInclusion], l))
                                error("could  !  change r==", r, " so that Subspets(r) contains", l, "\n")
                            end
                        end
                    end
                    p = Series(R, SubSpets(R, l, s[:element] // R.phi),
                               s.cuspidal, s.d)
                    p[:orbit] = ((WGL(s))[:orbitRepresentative])[i]
                    r = map((x->begin
                              Position((WGL(s))[:reflections], x)
                             end), gapSet((WGL(s))[:reflections]))
                    if (WGL(s))[:generators] == (WGL(s))[:parentMap]
                        r = Filtered(r, (x->begin
                                          ((WGL(s))[:reflections])[x] in
                                        Group(p.spets)
                                    end))
                    else
                        w = map((x->begin
                                  GetWord(WGL(s), ((WGL(s))[:reflections])[x])
                                    end), r)
                        w = map((x->begin
                                  Product(((WGL(s))[:parentMap])[x])
                                    end), w)
                        r = ListBlist(r, map(function (x,)
                                        local g
                                        g = Group(p.spets)
                                        if IsRec(x)
                                            return Representative(x[:element]) in g
                                        else
                                            return x in g
                                        end
                                    end, w))
                    end
                    p[:WGL] = ReflectionSubgroup(WGL(s), r)
                    p[:WGLdims] = map((x->begin
                                    x[1]
                                end), (CharTable(p[:WGL]))[:irreducibles])
                    if all((x->begin
                                    x == 1
                                end), p[:WGLdims])
                        p[:e] = length(p[:WGLdims])
                    end
                    return p
                   end, (WGL(s))[:reflists])
        s[:relativeSpets] = map((x->begin
                        x.spets
                    end), res)
        p = map(getHecke, res)
        if false in p
            return res
        end
        s[:Hecke] = Hecke(WGL(s), p)
        u1 = map(Mvp, SchurElements(s[:Hecke]))
        if any((x->begin
                        any((y->begin
                                    !(all(IsInt, y[:coeff]))
                                end), x[:elm])
                    end), u1)
            ChevieErr(s[:Hecke], " wrong set of SchurElements")
            return res
        end
        u1 = map((x->begin
                        s.degree // CycPol(Mvp(x))
                    end), u1)
        ud = map((x->begin
                        x * sign(Value(x //
                                       (CycPolUnipotentDegrees(s.levi))[s.cuspidal],
                                       E(denominator(s.d), numerator(s.d))))
                       end), (CycPolUnipotentDegrees(s.spets))[charNumbers(s)])
        p = PermListList(u1, ud)
        if p == false
            ChevieErr(s[:Hecke], " wrong set of SchurElements")
            return res
        end
        charNumbers(s) = Permuted(charNumbers(s), p)
        if haskey(s, :span)
            s[:span] = Permuted(s[:span], p)
        end
        aA = map((x->begin
                        E(denominator(s.d) ^ 2, x[:valuation] + degree(x))
                    end), u1)
        p = PositionRegularClass(WGL(s), s.d)
        if p == false && length(WGL(s)) == 1
            p = 1
        end
        o = map((x->begin
                        x[p] // x[1]
                       end), (CharTable(WGL(s)))[:irreducibles])
        s[:predictedEigen] = map((i->begin
                            aA[i] * o[i]
                        end), 1:length(o)) *
        (Eigenvalues(UnipotentCharacters(s.levi)))[s.cuspidal]
        return res
    end
Checkzegen = function (W,)
        local l, uc, aA, e, z, s, ucL, aAL
        l = AllProperSeries(W)
        print("with  no Hecke:", Filtered(l, (s->begin
                        Hecke(s) == false
                    end)), "\n")
        print("    center==", OrderCenter(W), "\n")
        l = Filtered(l, (s->begin
                        haskey(s, :Hecke)
                    end))
        uc = UnipotentCharacters(W)
        aA = uc[:a] + uc[:A]
        for s = l
            e = gete(s[:Hecke])
            z = OrderCenter(WGL(s))
            ucL = UnipotentCharacters(s.levi)
            aAL = ucL[:A] + ucL[:a]
            aAL = aAL[s.cuspidal]
            print(s, "\n")
            print(FormatTable([map(Mod1, (aA[charNumbers(s)] - aAL) // z),
                               map(Mod1, aA[charNumbers(s)] // z), map((x->begin
                                    Mod1(1 // x)
                                end), e)], Dict{Symbol, Any}(:rowLabels => ["(aA-aAL)/z", "aA/z", "local Hecke irr"])))
        end
    end
CheckRatCyc = function (s,)
        local l, fact, F
        if !(haskey(s, :Hecke))
            return SPrint(FormatTeX(s), " Hecke  !  bound")
        end
        l = CollectBy(((s[:Hecke])[:parameter])[1], degree)
        fact = map((x->begin
                        Product(x, (y->begin
                                    Mvp("T") - y
                                end))
                    end), l)
        F = FieldOfDefinition(Group(s.spets))
    end
CheckRatSer = function (arg...,)
        local W, s, c, b
        W = ApplyFunc(ComplexReflectionGroup, arg)
        s = AllProperSeries(W)
        s = Filtered(s, (x->begin x.spets != x.levi end))
        c = Join(map(CheckRatCyc, s), "\\hfill\\break\n")
        return c
    end
CheckPiGPiL = function (n,)
        local W, s, D0
        W = ComplexReflectionGroup(n)
        s = AllProperSeries(W)
        map(Hecke, s)
        s = Filtered(s, (ser->begin
                        haskey(ser, :e) && ser.principal
                    end))
        D0 = (s->begin
                    Sum(ReflectionDegrees(Group(s.spets)) +
                        ReflectionCoDegrees(Group(s.spets))) -
                    Sum(ReflectionDegrees(Group(s.levi)) +
                        ReflectionCoDegrees(Group(s.levi)))
                end)
        return map((x->begin
                        D0(x) // x[:e]
                    end), s)
    end
Checkdovere = function (n,)
        local W, s
        W = ComplexReflectionGroup(n)
        s = AllProperSeries(W)
        map(Hecke, s)
        s = Filtered(s, (ser->begin
                        haskey(ser, :e) && ser.principal
                    end))
        return Filtered(s, (x->begin
                        x[:e] // x.d > 2
                    end))
    end
CheckxiL = function (n,)
        local W, s
        W = ComplexReflectionGroup(n)
        s = AllProperSeries(W)
        map(Hecke, s)
        s = Filtered(s, (ser->begin
                        haskey(ser, :e) && ser.principal
                    end))
        return Filtered(s, (x->begin
                        AsRootOfUnity(PhiOnDiscriminant(x.levi)) * x.d != 0
                    end))
    end
CheckCoN = function (i,)
        local W, s
        W = ComplexReflectionGroup(i)
        s = AllProperSeries(W)
        map(Hecke, s)
        s = Filtered(s, (x->begin
                        length(x.levi) == 1 && haskey(x, :Hecke)
                    end))
        return map(function (ser,)
                    local m, sg
                    m = DeterminantMat(reflrep(W, ser[:element]))
                    sg = (-1) ^ Sum(ReflectionDegrees(ser[:WGL]) - 1)
                    if length(HyperplaneOrbits(ser[:WGL])) > 1
                        return false
                    else
                        return m * sg * Product(((ser[:Hecke])[:parameter])[1]) ^ ((HyperplaneOrbits(ser[:WGL]))[1])[:N_s]
                    end
                end, s)
    end
CheckLuCox = function (s,)
        local p, p1
        p = Product(Filtered(((s[:Hecke])[:parameter])[1], (x->begin
                            x != 1
                        end)), (x->begin
                        1 - x
                    end))
        p1 = Value(s.degree, Mvp("q"))
        return [p, p1]
    end
datafor = function (W,)
        local l, e
        l = map((d->begin
                        PrincipalSeries(W, d)
                    end), RegularEigenvalues(W))
        e = Sum(ReflectionDegrees(W) + ReflectionCoDegrees(W))
        return map((s->begin
                        [s.d, gcd(e * s.d, denominator(s.d)) // gcd(gcd(map((x->begin
                                                    x[:N_s]
                                                end),
                                                                            HyperplaneOrbits(WGL(s)))),
                                                                    denominator(s.d))]
                    end), l)
    end
cent = function (W,)
        local l, s, res, e
        l = map((d->begin PrincipalSeries(W, d) end), RegularEigenvalues(W))
        res = map((s->begin [s.d, ReflectionName(WGL(s))] end), l)
        l = gapSet(map((x->begin
                            x[2]
                        end), res))
        res = map((x->begin [map((z->begin z[1] end), Filtered(res, (y->begin
                                            y[2] == x end))), x] end), l)
        sort!(res)
        return map((x->begin SPrint(Join(x[1]), ":", x[2]) end), res)
    end
EigenspaceNumbers = function (W,)
        local d, e
        d = ReflectionDegrees(W)
        if !(IsSpets(W))
            return Union(map(divisors, d))
        end
        e = map((p->begin
                        Lcm(p[1], denominator(AsRootOfUnity(p[2])))
                    end), d)
        e = Union(map(divisors, e))
        return Filtered(e, (x->begin
                        any((p->begin
                                    E(x, p[1]) == p[2]
                                end), d)
                    end))
    end
nrconjecture = function (W,)
        local d, e, p, r, f
        e = Difference(EigenspaceNumbers(W), RegularEigenvalues(W))
        e = Filtered(e, (x->begin
                        !(x[1]) in r
                    end))
        for d = e
            f = Filtered(r, (x->begin
                            mod(x, d) == 0
                        end))
            if length(f) > 0
                print("**** failed: ", ReflectionName(W), d, f, "\n")
            end
            p = First(reverse(e), (i->begin
                            mod(i[1], d) == 0
                        end))
            if p != d
                if length(RelativeDegrees(W, d)) * phi(d) != p[2] * phi(p[1])
                    print("**** failed: ", ReflectionName(W), p, d, "\n")
                end
            end
        end
    end
minimaldparabolic = function (W, d)
        local l
        l = map((i->begin
                        Difference(W[:generatingReflections], [i])
                    end), W[:generatingReflections])
        l = Filtered(l, (x->begin
                        length(RelativeDegrees(ReflectionSubgroup(W, x), d)) == length(RelativeDegrees(W, d))
                    end))
        l = map((x->begin
                        ReflectionName(ReflectionSubgroup(W, x))
                    end), l)
        return l
    end
