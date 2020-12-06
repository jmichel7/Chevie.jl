# CuspidalPairs(W[,d[,ad]]) returns the pairs (LF,λ) where LF is a d-split
# Levi [with d-center of dimension ad] and lambda a d-cuspidal character of LF
function CuspidalPairs(W,d,ad)
  WF=(W isa Spets) ? W : spets(W)
  if (d isa Int) && d!=0 d=1//d else d=d//1 end
  vcat(map(HF->map(char->[HF,char],cuspidal(UnipotentCharacters(HF), d)),
            split_levis(WF, d, ad))...)
end

function CuspidalPairs(W,d=0)
  WF=(W isa Spets) ? W : spets(W)
  if (d isa Int) && d!=0 d=1//d else d=d//1 end
  vcat(map(ad->CuspidalPairs(WF,d,ad),0:length(relative_degrees(WF,d)))...)
end

# s is a Set of tuples. Return E_1,...,E_n such that
# s=List(Cartesian(E_1,...,E_n),Concatenation)
# Assumes all E_i but one are of size 2
function FactorsSet(s)
  if s == [[]] return [] end
  for i = s[1]
    s1 = filter(y->i in y,s)
    if length(s1) == length(s) // 2
      i = intersect(s1...)
      r = setdiff(s[1], i)
      s2 = filter(y->length(intersect(y, i))==0,s)
      a = findfirst(y->issubset(r,y),s2)
      if length(s2) == length(s) // 2 && a != false
        j = setdiff(s2[a], r)
        s1 = sort(unique(map(x->setdiff(x, i), s1)))
        s2 = sort(unique(map(x->setdiff(x, j), s2)))
        if length(i)==1 i=i[1] end
        if length(j)==1 j=j[1] end
        if s1 == s2 return vcat([[i, j]], FactorsSet(s1))
        end
      end
    end
  end
  [s]
end

LIMSubsetsSum = 10
# l is a matrix and S a list of same length as a row of l.
# Find subsets P of [1..Length(l)] such that Sum(l{P})=S.
# in  addition, lv is a vector of same  length as l, v is a sub-multiset of
# lv and the chosen subsets should satisfy lv{P}=v as multisets.
function SubsetsSum(S, l, v, lv)
# println("S=$S;l=$l;v=$v;lv=$lv")
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
  # s assumed to be in P at this stage
  # S= initial S minus Sum(l{s})
  # t= remaining elements of eachindex(l) which could be in P
  # nonsolved= indices of nonsolved entries of S
  # v= remaining v to match
  inner = function (S, s, t, nonsolved, v, factor)
    local bad, good, p, sols, res, i, sol, f, ll, solved
  # Print("#solved=",Length(l[1])-Length(nonsolved)," ",Join(s),"=>",Join(t),
  #       " v=",Join(List(Collected(v),x->Join(x,":"))," "),"\n");
    c+=1
    if mod(c, 1000) == 0
      InfoChevie("# ",factor,": xcols:",length(nonsolved)-1," xrows:",length(t),
                 " found:", found, "\n")
    end
    t = filter(i->lv[i] in v,t)
    if length(t) == 0
      if iszero(S) found+=1
          return [[]]
      else return []
      end
    end
    ll=map(i->Dict{Symbol,Any}(:pos=>i,:cand=>filter(j->l[j][i]!=0,t)),nonsolved)
    sort!(ll,by=x->length(x[:cand]))
    if length(ll[1][:cand])>LIMSubsetsSum ll = [ll[1]]
    else ll=filter(x->length(x[:cand])<=LIMSubsetsSum,ll)
    end
    solved = Int[]
    good = Int[]
    bad = Int[]
    for p in ll
     p[:sols]=filter(e->sum(getindex.(l[e],p[:pos]))==S[p[:pos]],combinations(p[:cand]))
      if length(p[:sols])==0 return []
      elseif length(p[:sols]) == 1 push!(solved, p[:pos])
      end
      if sum(length,p[:sols])>0 
        good=union(good,intersect(p[:sols]...)) #lines part of any solution
        bad=union(bad, setdiff(p[:cand], union(p[:sols])))#part of no solution
      else
        bad=union(bad,p[:cand]) #part of no solution
      end
    end
    nonsolved = setdiff(nonsolved, solved)
    if length(good) + length(bad) > 0 #progress
      return map(r->vcat(good, r), inner(S - sum(l[good]), union(s, good), 
         setdiff(t, union(good, bad)), nonsolved, sievev(good, v), factor))
    else
      res = []
      p = maximum(map(x->length(x[:cand]),ll))
      p = ll[findfirst(x->length(x[:cand])==p,ll)]
      f = length(p[:sols])
      nonsolved = setdiff(nonsolved, [p[:pos]])
      InfoChevie("# ", factor, ": xcols:", length(nonsolved), " xrows:", 
         length(t)," in comb(",length(p[:cand]),")==>", length(p[:sols]), "\n")
      for sol = p[:sols]
        good = sol
        bad = setdiff(p[:cand], sol)
        if isempty(intersect(good, bad))
           append!(res, map(r->vcat(good,r),
            inner(isempty(l[good]) ? S : S-sum(l[good]),union(s, good), 
          setdiff(t, union(good, bad)), nonsolved, sievev(good, v), 
          string(factor, ":", f))))
        end
        f-=1
      end
      return res
    end
  end
  return inner(S, [], 1:length(l), 1:length(S), v, "")
end

positive(p::CycPol)=all(>(0),values(p.v))

# FitParameter(sch,m) given:
# sch: schur elements for H(Z/e) given as CycPols 
# m:   a list of length e of Rationals, exponents of params
#
# finds all permutations σ of 1:e such that the parameters pₖ=E(e,σₖ-1)q^mₖ
# gives  sch by the formula schᵢ=∏_{j≠i}(1-pᵢ/pⱼ). Since multiplying the pₖ
# by  a  scalar  leaves  invariant  the  sch,  σ  is  known  only  modulo a
# translation.
#
# The  result is  a list  of pairs  [v1,v2] telling that globally σ(v1)=v2,
# where the v1 are sort(collect(values(groupby(sch,eachindex(sch)))),by=minimum)
function FitParameter(sch, m)
  den=lcm(denominator.(m))
  e=length(m)
  a=tally(m)
  sch=map(x->descent_of_scalars(x, den),sch) # m_k replaced by den*m_k
  function poss(i)
    # for each element (m_k,c_k) of tally(m) p will hold a minimal
    # corresponding possible set of j=σ(i)-σ(l) for l such that m_l=m_k
    p=map(x->collect(0:e-1),a)
    avail=eachindex(a)
    good=Int[]
    term(j,k)=CycPol(1-E(e,j)*Pol()^(den*(m[i]-first(a[k]))))
    v=sch[i]
    while true
      for k in avail
        p[k]=setdiff(p[k],vcat(p[good]...))
        if m[i]!=first(a[k]) p[k]=filter(j->positive(v//term(j,k)),p[k]) end
      end
      good=filter(i->length(p[i])==last(a[i]),avail)
      if isempty(good) return p end
      avail=setdiff(avail, good)
      v//=reduce(*,map(k->prod(j->term(j,k),p[k]),
                       filter(k->m[i]!=first(a[k]),good));init=1)
    end
  end
  pp=poss.(1:e)
  p=findfirst(i->last.(a)==length.(pp[i]),1:e)
  if isnothing(p) error("no solution\n") end
  v=map(i->map(x->sort(unique(mod.(x,e))),i+pp[p]),0:e-1)
  v=map(k->filter(i->all(j->issubset(v[i][j],pp[k][j]),1:length(a)),1:e),1:e)
  G=map(x->findfirst(==(x),sch), sch)
  v=map(x->v[x], sort(unique(G)))
  G=map(x->findfirst(==(x),sort(unique(G))), G)
  bad=1:length(v)
  while true
    good=filter(i->length(v[i])==count(j->G[j]==i,1:e),bad)
    bad=setdiff(bad, good)
    for p in bad for a in good v[p]=setdiff(v[p],v[a]) end end
    if isempty(good) break end
  end
  if sum(length,v)>e error("non-unique solution\n") end
  res=collectby(G,eachindex(G))
  pp=vcat(v...)[sortperm(vcat(res...))]
  para=map(i->E(e,pp[i])*Pol()^(den*m[i]),1:e)
  if map(k->prod(j->CycPol(1-para[k]//para[j]),filter(j->j!=k,1:e)),1:e)!=sch
    error("schur elms don't match\n")
  end
  map((x,y)->[x,y],res, v)
end

#---------------------- Series -----------------------------------------
# A d-Harish-Chandra series \CE(\BG,(\BL,λ) is a record with fields:
#  .d           AsRootOfUnity(ζ) such that L=C_G(V_ζ)
#  .spets       \BG
#  .levi        \BL
#  .classno     class of s.levi with ζ-eigenspace V_ζ
#  .element     representative of classno
#  .projector   .element-equivariant projector on V_ζ
#  .cuspidal    λ (index in UnipotentCharacters(.levi))
#  .principal   true iff λ=\Id
#  .WGL         W_\BG(\BL,λ) as a relgroup, contains parentMap (refs->elts of W)
#               and reflists (generators->gens of parab of W)
#  .WGLdims     irr dims of WGL
#  .charNumbers \CE(\BG,(\BL,λ) (indices in UnipotentCharacters(.spets))
#  .degree      The degree of R_\BL^\BG(λ)=|G/L|_{q'} degλ
#  .eps         For each χ in .charNumbers sign of <χ,R_\BL^\BG(λ)>
#  .dims        χ(1) for γ_χ
#  .Hecke       H_G(L,λ)
#
# Fields only present when WGL is cyclic
# .e         |WGL|
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
# find simplest regular eigenvalue q of s.levi
  eig=union(map(x->x isa Int ? prime_residues(x).//x : [x], 
                map(x->x.r,regular_eigenvalues(levi)))...)
  c=minimum(denominator.(eig))
  c=minimum(filter(x->denominator(x)==c,eig))
  c=Root1(;r=c)
  eig=refleigen(levi)
  q=maximum(map(x->count(==(c),x), eig))
  classno=findall(x->count(==(c),x)==q,eig)
  if length(classno)>1 error("classno==",classno) end
  element=levi(classinfo(levi)[:classtext][classno[1]]...)
  q=minimum(map(x->count(==(d),x),refleigen(levi)))
  q=findall(x->count(==(d),x)==q,refleigen(levi))
  projector=eigenspace_projector(WF, classreps(levi)[q[1]], d)
  q=Pol()
  deg=conj(generic_sign(WF))*generic_order(WF,q) // 
      (conj(generic_sign(levi))*generic_order(levi,q))
  deg=CycPol(shift(deg,-deg.v))*Uch.CycPolUnipotentDegrees(levi)[cuspidal]
  Series(WF,levi,cuspidal,d,principal,classno[1],element,projector,
         deg,Dict{Symbol,Any}())
end

CuspidalSeries(W,d=1)=map(x->Series(W,x[1],x[2],d),CuspidalPairs(W,d))
CuspidalSeries(W,d,ad)=map(x->Series(W,x[1],x[2],d),CuspidalPairs(W,d,ad))

function Base.show(io::IO,s::Series)
  TeX=get(io,:TeX,false)
  repl=get(io,:limit,false)
  if !(repl || TeX)
    print(io,"Series($(s.spets),$(s.levi),$(s.d),$(s.cuspidal))")
    return
  end
  cname = charnames(io,UnipotentCharacters(s.levi))[s.cuspidal]
  n=repl || TeX ? "\\lambda" : "c"
  quad=TeX ? "\\quad" : " "
  if s.spets == s.levi
    print(io,s.d, "-cuspidal ",cname," of ", s.spets)
  else
    print(io,s.d, "-series ")
    print(io,"R^{",s.spets,"}_","{",s.levi,"}(")
    Util.printTeX(io,"$n==",cname,")")
  end
  if haskey(s.prop, :e)
    Util.printTeX(io,"$quad W_G(L,$n)==Z_{$(s.prop[:e])}")
  elseif haskey(s.prop, :WGL)
    if haskey(relative_group(s).prop, :refltype)
      Util.printTeX(io,"$quad W_G(L,$n)==");print(io,relative_group(s))
    else
      Util.printTeX(io,"$quad |W_G(L,$n)|==$(length(relative_group(s)))")
    end
  elseif haskey(s.prop, :WGLdims)
    Util.printTeX(io,"$quad |W_G(L,$n)|==$(sum(WGLdims(s).^2))")
  end
  if haskey(s.prop, :translation)
    Util.print(io,"$quad translation==",s.prop[:translation])
  end
  if !haskey(s.prop,:charNumbers) return end
  if get(io,:typeinfo,Any)==typeof(s) return end
end

function Base.show(io::IO, ::MIME"text/plain", s::Series)
  show(io,s)
  if haskey(io,:typeinfo) return end
  if haskey(s.prop,:charNumbers) && !isnothing(CharNumbers(s)) format(io,s) end
end

function Util.format(io::IO,s::Series)
  uw = UnipotentCharacters(s.spets)
  e = length(CharNumbers(s))
  TeX=get(io,:TeX,false)
  repl=get(io,:limit,false)
  function f(texn,n,val)
    push!(rowlab, repl || TeX ? fromTeX(io,texn) : n)
    push!(m, map(x->sprint(show,x;context=rio(stdout)), val))
  end
  m = []
  rowlab = []
  f("\\hbox{Character}", "Name", charnames(io,uw)[CharNumbers(s)])
  f("\\hbox{specializes}", "specializes", charnames(io,relative_group(s)))
  f("\\varepsilon", "eps", s.prop[:eps])
  if haskey(s.prop, :eigen) f("\\hbox{eigen}", "eigen", s.prop[:eigen]) end
  if haskey(s.prop, :e)
    if haskey(s.prop, :Hecke) f("\\hbox{parameter}", "param", s.prop[:Hecke].para[1])
    elseif haskey(s.prop, :mC) f("\\hbox{mC}", "degparam", mC(s))
    end
  end
  f("\\hbox{family \\#}", "family #", map(j->findfirst(k->j in k[:charNumbers],
                                    uw.families), CharNumbers(s)))
  if haskey(s.prop, :permutable) && any(x->x!=false,s.prop[:permutable])
      f("\\hbox{permutable}", "permutable", s.prop[:permutable])
  end
  m = permutedims(toM(m))
  println(io)
  format(io,m;row_labels=CharNumbers(s),col_labels=rowlab)
end

ChevieErr(x...)=xprint("!!!!!!! ",x...)

# Degree in q of the parameters (normalized so the smallest is 0)
function mC(s::Series)
  gets(s,:mC) do
  e = hyperplane_orbits(relative_group(s))
  if length(e)>1 return else e = e[1].order end
  uc = UnipotentCharacters(s.spets)
  cn = CharNumbers(s)[filter(i->s.prop[:dims][i]==1,1:length(s.prop[:dims]))]
  aA = uc.prop[:a][cn] + uc.prop[:A][cn]
  lpi(W)=sum(degrees(W)+codegrees(W))
  if s.principal
    if minimum(aA)!=0 error("id not in RLG(1)") end
    pG=lpi(Group(s.spets))
    pL=lpi(Group(s.levi))
    D0=pG-pL
    xiL = Root1(PhiOnDiscriminant(s.levi))^conductor(s.d)
    xiG = Root1(PhiOnDiscriminant(s.spets))^conductor(s.d)
    if xiL != xiG
      ChevieErr("fixing dimension of variety by xiL-xiG==", xiL - xiG, "\n")
      D0 = (D0 + xiL) - xiG
    end
    # Id in chars
    if !isinteger(D0*s.d.r)
      ChevieErr(s, "\n    ==> (l(pi_G)==", pG, ")-(l(pi_L)==", pL,
                ")+(xiL*d==", xiL, ")-(xiG*d==", xiG, ")  !  ==0( %
         d==", s.d, ") Max(aA)-Min(aA)==", maximum(aA) - minimum(aA), "\n")
    end
    if !isinteger(D0//e)
      ChevieErr("fixing dimension of variety to be divisible by e==", e, "\n")
      D0-=mod(D0, e)
    end
    (D0-aA)//e
  else
    # (JM+Gunter 18-3-2004) in any case we normalize so that
    # the smallest mC is 0 since the above choice gives unpleasant
    # results for G27
    D0=maximum(aA)-minimum(aA) # just so that mC are positive
    (D0+minimum(aA).-aA)//e
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
    WGL(s)[:parentMap] = gens(WGL(s))
    WGL(s)[:reflists] = map(x->[x],inclusiongens(W))
  else
    V = NullspaceMat(s.projector - id)
    rH = map(r->(SumIntersectionMat(NullspaceMat(mats[r] - id), V))[2], refs)
    rH = gapSet(Filtered(rH, (H->begin length(H) < length(V) end)))
    rel = []
    for H = rH
        if length(H) == 0
            ConjugacyClasses(Group(s.levi))
            H = W
        else
            H = Filtered(refs, (r->begin H * mats[r] == H end))
            H = ReflectionSubgroup(W, (W[:rootInclusion])[H])
        end
        push!(rel, Dict{Symbol, Any}(:refs => inclusiongens(H),
                                     :WH => Normalizer(H, Group(s.levi))))
    end
    if (s.spets).phi == Perm()
        WF = W
        L = Group(s.levi)
    else
        WF = Group(Concatenation(W[:generators], [s.spets.phi]), W[:identity])
        L = Subgroup(WF, (Group(s.levi))[:generators])
    end
    if length(L) == 1
     WGL(s) = Centralizer(W, (s.levi).phi)
        for H in rel H[:WH] = Centralizer(H[:WH], (s.levi).phi) end
        func = (x->begin x end)
    else
      NLF = Normalizer(WF, L)
      if NLF == L
        if s.levi != s.spets
          error(s, " N==L\n")
          return false
        end
        WGL(s)=CoxeterGroup()
        WGL(s)[:reflists] = []
        s[:WGLdims] = [1]
        return WGL(s)
      end
      NLFb = FactorGroup(NLF, L)
      hom = NaturalHomomorphism(NLF, NLFb)
      phi = Image(hom, (s.levi).phi)
      s[:WGLdegrees] = RelativeDegrees(s.spets, s.d)
      if count(x->x!=1,RelativeDegrees(s.levi,s.d)) == 0 && 
         count(x->x!=1,s[:WGLdegrees])==1
        if s.levi.phi in W && Order(NLFb,phi)== Product(s[:WGLdegrees])
          WGL(s) = Subgroup(NLFb, [phi])
        else
         N = Centralizer(Subgroup(WF,gens(W)),s.levi.phi)
          if length(N)==length(L)*Product(s[:WGLdegrees])
            WGL(s)=FactorGroup(N, L)
          end
        end
      end
      if !haskey(s,:WGL)
        N=Subgroup(WF, W[:generators])
        N=Subgroup(N, Intersection(N, NLF)[:generators])
        WGL(s) = Centralizer(FactorGroup(N, L), phi)
      end
      if length(rel) == 1 && rel[1][:WH] == NLF rel[1][:WH] = WGL(s)
      else
        rel=map(rel)do x
          local N, WH
          N=Intersection(WGL(s),FactorGroup(Subgroup(NLF,gens(x[:WH])),L))
          if IsGroup(N) x[:WH] = N
          else x[:WH] = Subgroup(WGL(s), N)
          end
          return x
        end
      end
      func = (x->begin (x.element)[:representative] end)
    end
    ud = CycPolUnipotentDegrees(s.levi)
    eig = Eigenvalues(UnipotentCharacters(s.levi))
    ud = Filtered(1:length(ud), (i->begin
                    ud[i] == ud[s.cuspidal] && eig[i] == eig[s.cuspidal] end))
    if length(ud) > 1
      c = length(WGL(s))
      WGL(s) = Stabilizer(WGL(s), Position(ud, s.cuspidal), function (c, g)
               return c^Uch.on_unipotents(s.levi, func(g), ud) end)
      if c!=length(WGL(s))ChevieErr("# WGL:",c,"/",length(WGL(s))," fix c\n")
        end
      for H in rel H[:WH] = Intersection(H[:WH], WGL(s)) end
      rel = Filtered(rel, (x->begin length(x[:WH]) > 1 end))
    end
    if !(all((x->begin IsCyclic(x[:WH]) end), rel)) error("a") end
    hom=map(x->First(elements(x[:WH]),(y->Order(x[:WH],y)==length(x[:WH]))),rel)
    if Subgroup(WGL(s), hom) != WGL(s) error("b") end
    reflists = map((x->begin x[:refs] end), rel)
    m = Concatenation(V, NullspaceMat(s.projector)) ^ -1
    rel = map(hom)do x
      local M
      M=reflrep(W, func(x)) ^ m
      AsReflection(M[1:length(V)][1:length(V)])
    end
    for c in 1:length(rel)
      r = AsRootOfUnity((rel[c])[:eigenvalue])
      m = numerator(r)
      r = denominator(r)
      if m != 1
      rel[c]=AsReflection(Reflection(rel[c][:root],rel[c][:coroot])^mod(1//m,r))
      end
    end
    WGL(s)=PermRootGroup(map(x->x[:root],rel),map(x->x[:coroot],rel))
    WGL(s)[:parentMap] = hom
    WGL(s)[:reflists]=reflists[map(x->PositionProperty(rel,
                  (y->(ratio(y[:root],x)!=false))),
                WGL(s)[:roots][WGL(s)[:generatingReflections]])]
  end
  s[:WGLdims]=CharTable(WGL(s)).irr[:,1]
  if all(isone,WGLdims(s)) s.prop[:e]=length(WGLdims(s)) end
  WGL(s)
end
end
else
function Weyl.relative_group(s::Series)
  gets(s,:WGL) do
  W=Group(s.spets)
  L=Group(s.levi)
  if isone(s.projector) #central series
    WGL=deepcopy(W)
    WGL.prop[:parentMap] = gens(WGL)
    WGL.prop[:reflists] = map(x->[x],inclusiongens(W))
    s.prop[:WGLdims]=CharTable(WGL).irr[:,1]
    if iscyclic(WGL) s.prop[:e] = length(WGL) end
    return WGL
  end
  if W isa FiniteCoxeterGroup N=normalizer(W.G, L.G)
  else N=normalizer(W, L)
  end
  if !isone(s.levi.phi)
    if length(L) == 1
      N=centralizer(N, s.levi.phi)
    elseif s.spets isa Cosets.CoxeterCoset
      N=Group(map(x->reduced(Group(s.levi), x),gens(N)))
      N=centralizer(N, s.levi.phi)
      N=Group(vcat(gens(N),gens(Group(s.levi))))
    else #   N:=Stabilizer(N,s.levi); # is shorter but slower...
      NF=centralizer(N, s.levi.phi)
      if length(NF)==length(L)*prod(relative_degrees(s.spets, s.d)) N=NF
      else N=centralizer(N,s.levi.phi;action=function(p,g)
                  w=reduced(s.levi.W,p^g)
                  w isa Perm ? w : w.phi
                 end)
      if false
        NF=Group(vcat(gens(N),[s.levi.phi]))
        NFQ=NF/L
        N=N/L
        NF=centralizer(NFQ, s.levi) # image of s.levi.phi in NFQ
        NFQ=intersect(NF, N)
        N=Group(vcat(gens(L),map(x->x.phi,gens(NFQ))))
      end
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
      N=centralizer(N,findfirst(==(s.cuspidal),ud);action=(c,g)->
                  c^Uch.on_unipotents(s.levi, g, ud))
      if c!=length(N)
        ChevieErr("# WGL:",length(N)//length(L),"/",c//length(L)," fix c\n")
      end
    end
  end
  refs=unique(sort(indexin(reflections(W), reflections(W))))
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
    s.prop[:e]=[1]
    return WGL
  else
    WGL=N/L
    func = x->x.phi
  end
  V=GLinearAlgebra.lnullspace(s.projector-one(s.projector))#The E(d)-eigenspace
  m=vcat(V, GLinearAlgebra.lnullspace(s.projector))^-1
  # V∩ fix(r)
  hplane(r)=sum_intersection(GLinearAlgebra.lnullspace(reflrep(W,
                     reflection(W,r))-one(s.projector)),V)[2]
  smalltobig(h)=hcat(h, fill(0,size(h,1),max(0,rank(W)-size(V,1))))*m^-1
  function bigtosmall(x)
    x=x^m
    x[axes(V,1),axes(V,1)]
  end
  function getreflection(rr)
    local rH, H, r, res, n
    rH = hplane(rr[1]) # Hyperplane of W_G(L)
    if length(rH)==length(V) error("intersection not proper") end
    if length(rH)==0 # cyclic WGL
      H=W
    else
      H=reflection_subgroup(W, vcat(inclusion(W,rr), inclusiongens(L)))
    end
    res=Dict{Symbol, Any}(:refs => inclusiongens(H))
    H = intersect(H, N)
    if length(L)!=1 H=H/L end
    if length(H)==1 error("H trivial") end
    ee=elements(H)
    gen=findfirst(y->order(y)==length(H),ee)
    if isnothing(gen) error("H not cyclic") end
    res[:hom]=ee[gen]
    r=bigtosmall(reflrep(W, func(res[:hom])))
    ref=reflection(r)
    n=ref.eig
    r^=invmod(exponent(n), conductor(n))
    merge!(res,pairs(reflection(r)))
    res[:WH] = H
    res
  end
  function v(h) # normalize a space (what "VectorSpace" could do)
    if size(h,1)==0 return  h end
    h=GLinearAlgebra.echelon(h)[1]
    h[1:count(!iszero,eachrow(h)),:]
  end
  rrefs = collect(values(groupby(x->v(hplane(x)),refs)))
  rrefs=filter(x->!(any(y->inclusion(W,y) in inclusion(L),x)),rrefs)
  sort!(rrefs)
  reflist = []
  for r in rrefs
    push!(reflist, getreflection(r))
    if length(Group(map(x->x[:hom],reflist))) == length(WGL)
      println(map(x->x[:root],reflist),map(x->x[:coroot],reflist))
      WGL=PRG(map(x->x[:root],reflist),map(x->x[:coroot],reflist))
      reflist=map(x->smalltobig(GLinearAlgebra.lnullspace(x-x^0)),reflrep(WGL))
      reflist=map(h->rrefs[findfirst(rr->v(hplane(rr[1]))==v(h),rrefs)], reflist)
      WGL.prop[:reflists] = map(getreflection, reflist)
      WGL.prop[:parentMap] = map(x->x[:hom], WGL.prop[:reflists])
      WGL.prop[:reflists] = map(x->x[:refs], WGL.prop[:reflists])
      s.prop[:WGLdims]=CharTable(WGL).irr[:,1]
      if iscyclic(WGL) s.prop[:e]=length(WGL) end
      return WGL
    end
  end
  error("b")
end
end
end

function Groups.iscyclic(s::Series)
  relative_group(s)
  haskey(s.prop,:e)
end

WGLdims(s::Series)=getp(relative_group,s,:WGLdims)
    
e(s::Series)=getp(relative_group,s,:e)
    
function RLG(s::Series)
  gets(s,:RLG) do
  RLG=LusztigInduce(s.spets, UniChar(s.levi, s.cuspidal))
  if isnothing(RLG) && isone(s.levi.phi)
    RLG=HarishChandraInduction(s.spets, UnipotentCharacter(s.levi,s.cuspidal))
  end
  if isnothing(RLG) ChevieErr(s, ":RLG failed\n")
  elseif s.deg!=CycPol(degree(RLG))
    ChevieErr(s,":Deg RLG!=Sum(eps[i]*ud[i])\n")
  end
  RLG
  end
end

# takes a d-series s with s.spets split; fills in s.charNumbers, s.eps, s.dims
function CharNumbers(s::Series)
  gets(s,:charNumbers) do
  if s.levi==s.spets return [s.cuspidal] end
  if !haskey(s.prop,:WGL) && isnothing(relative_group(s)) return nothing end
  ud=Uch.CycPolUnipotentDegrees(s.spets)
  ad=count(!isone,relative_degrees(s.levi, s.d))
  cand=filter(i->ad==valuation(ud[i],s.d),1:length(ud))
# now use that   S_\phi(q)=\eps_\phi Deg(RLG(λ))/Deg(γ_\phi)
  cand=map(c->Dict(:charNumbers=>c,:sch=>s.deg//ud[c]),cand)
  cand=filter(c->positive(c[:sch]),cand)
  q=Mvp(:q)
  ad=CycPol(Pol()-E(s.d))^ad
  f=s.deg//ad
  if !positive(f)
    ChevieErr(s, " cuspidal is not\n")
    return nothing
  end
  v=f(E(s.d))
  cand=filter(cand)do c
    c[:span]=degree(c[:sch])-valuation(c[:sch])
    f=sum(x->x^2,WGLdims(s))*(ud[c[:charNumbers]]//ad)(E(s.d))//v
    if !isinteger(f) return false end
    f=Int(f)
    if !(abs(f) in WGLdims(s)) return false end
    c[:dims]=abs(f)
    c[:eps] =sign(f)
    true
  end
  eig=Uch.eigen(UnipotentCharacters(s.levi))[s.cuspidal]
  eig*=map(i->E(conductor(s.d)^2,i),1:conductor(s.d)^2)
  cand=filter(c->Uch.eigen(UnipotentCharacters(s.spets))[c[:charNumbers]] in eig,cand)
  if length(cand)<length(WGLdims(s))
    ChevieErr(s,": not enough left with predicted eigenvalues in ",Root1.(eig),"\n")
    return nothing
  end
  sort!(cand,by=x->x[:dims])
  sort!(WGLdims(s))
  foo(n)=getindex.(cand,n)
  function foo(n,l)
    if l isa Vector return map(x->x[n], cand[l])
    else return cand[l][n]
    end
  end
  check=function ()
    local n
    for n in [:charNumbers, :eps, :dims, :span] s.prop[n]=foo(n) end
    if !isnothing(RLG(s)) && (filter(i->RLG(s).v[i]!=0,eachindex(RLG(s).v))!=
                              sort(s.prop[:charNumbers]) ||
              RLG(s).v[s.prop[:charNumbers]]!=s.prop[:dims].*s.prop[:eps])
      ChevieErr(s, ":RLG does not match")
    end
    s.prop[:charNumbers]
  end
  if length(cand)==length(WGLdims(s)) return check() end
  ud=foo(:dims).*Uch.CycPolUnipotentDegrees(s.spets)[foo(:charNumbers)].*foo(:eps)
  t=maximum(degree.(ud))
  function c(p)
   pp=p(Pol())
    vcat(fill(0,pp.v),pp.c,fill(0,max(0,t-degree(p))))
  end
  v = SubsetsSum(improve_type(c(s.deg)), improve_type(map(c, ud)), 
                 improve_type(WGLdims(s)), foo(:dims))
  InfoChevie("# ", length(v), " found\n")
  if length(v)>10000
    InfoChevie("# ", length(v), " combinations sum to dimRLG\n")
  elseif length(v) == 0
    ChevieErr(s, " no combination sums to dimRLG\n")
    return nothing
  end
  if iscyclic(s)
    s.prop[:charNumbers] = foo(:charNumbers)
    s.prop[:dims] = foo(:dims)
    mC(s)
    if length(v)>1
      v=filter(a->foo(:span,a)==map(i->
           sum(j->abs(mC(s)[i]-mC(s)[j]),setdiff(a,[i])),a),v)
      if length(v)>10000
          InfoChevie("# ", length(v), " combinations have right span\n")
      end
    elseif length(v) == 0
      ChevieErr(s, " no combination has right span\n")
      return nothing
    end
    delete!(s.prop,:charNumbers)
    delete!(s.prop,:dims)
    delete!(s.prop,:mC)
  end
  if length(v)>1
    InfoChevie("# after span ", length(v), " combinations\n")
    InfoChevie("# Warning: using Mackey with tori for ", s, "\n")
    i=fusion_conjugacy_classes(s.levi, s.spets)
    c=CharTable(s.spets).centralizers[i].//CharTable(s.levi).centralizers
    t=Uch.DLCharTable(s.levi)[:,s.cuspidal].*c
    t=map(k->sum(t[filter(j->i[j]==k,eachindex(i))]),
                       1:HasType.NrConjugacyClasses(s.spets))
    c=Uch.DLCharTable(s.spets)
    v=filter(a->c[:,foo(:charNumbers,a)]*(foo(:dims,a).*foo(:eps,a))==t,v)
  end
  if length(v)>1
    ChevieErr(s," ",join(FactorsSet(map(x->foo(:charNumbers,x),v)),"x"),
              " chars candidates: using RLG\n")
    if isnothing(RLG(s)) return nothing end
    v=filter(l->all(i->RLG(s).v[foo(:charNumbers,i)]!=0,l),v)
    if length(v)!=1 error() end
  elseif length(v)==0
    ChevieErr(s," no candidates left\n")
    return nothing
  end
  cand=cand[v[1]]
  check()
  end
end

eps(s::Series)=getp(CharNumbers,s,:eps)

COMPACTCOHOMOLOGY=true
function paramcyclic(s::Series)
  if !iscyclic(s) error("fill assumes relative_group(s) cyclic\n") end
  if isnothing(CharNumbers(s)) return nothing end
  uc=UnipotentCharacters(s.spets)
  Schur=Uch.CycPolUnipotentDegrees(s.spets)[CharNumbers(s)]
  Schur=map(x->s.deg//Schur[x]*eps(s)[x], 1:e(s))
  s.prop[:eigen]=Uch.eigen(uc)[CharNumbers(s)]
  LFrob=Root1(Uch.eigen(UnipotentCharacters(s.levi))[s.cuspidal])
  m=degrees(Group(s.spets))
  s.prop[:delta]=lcm(map(x->conductor(Root1(x[2])),
                         filter(x->x[1]!=1,degrees(s.spets))))
  rr(j,i)=(i-1)//e(s)-mC(s)[j]*s.d.r
  param(j,i)=Mvp(:q)^mC(s)[j]*E(Root1(;r=rr(j,i)))
  # parameters of Hecke algebra are map(i->param(i,i),1:e(s))
  mmp=FitParameter(ennola_twist.(Schur,E(s.d)),
                   COMPACTCOHOMOLOGY ? mC(s) : -mC(s))
# possible perms encoded by mmp[i][1]->mmp[i][2]
  nid=uc.almostHarishChandra[1][:charNumbers][charinfo(s.spets)[:positionId]]
  function predeigen(j, i)#eigenvalue for mC[j] and ζ_e^{i-1}
    if nid in CharNumbers(s)# id should correspond to id(WGL)
      Root1(;r=s.d.r*e(s)*s.prop[:delta]*
                   (rr(j,i)+s.d.r*mC(s)[findfirst(==(nid),CharNumbers(s))]))
    else Root1(;r=s.prop[:delta]*rr(j,i)*e(s)*s.d.r+LFrob.r)
              # as fraction predeigen_i:=delta di -delta m_i d^2e 
    end
  end
  series=map(x->findfirst(y->x in y[:charNumbers],uc.harishChandra),CharNumbers(s))
  unique=filter(p->length(p[1])==1,mmp)
  ratio=map(p->s.prop[:eigen][p[1][1]]//E(predeigen(p[1][1], p[2][1])), unique)
  if length(Set(ratio))>1
    ChevieErr(s, " eigenvalue ratios==", ratio, "\n")
    return nothing
  end
  ratio=Root1(ratio[1]).r
# now find integer translation t such that mod 1. we have t delta d=ratio
  ratio*=conductor(s.d^s.prop[:delta])
  if !isinteger(ratio)
    ChevieErr(s, "non-integral ratio==", ratio, "\n")
    return nothing
  end
  if ratio==0 r=0
  else r=s.d^s.prop[:delta]
    r=ratio*invmod(exponent(r),conductor(r))
  end
  mmp=map(x->(x[1],map(y->mod(y,e(s))+1,x[2]+r-1)), mmp)
  r=fill(0,e(s))
  s.prop[:permutable]=map(x->0,CharNumbers(s))
  j=1
  for i in mmp
    a=arrangements(i[2], length(i[2]))
    p=findall(A->s.prop[:eigen][i[1]]==map(j->E(predeigen(i[1][1],j)),A),a)
    if isempty(p)
      ChevieErr(s, "predicted eigenvalues cannot match actual\n")
      return nothing
    else
      if length(p)>1
        o=orbits(Group(map(x->Perm(a[p[1]],x),a[p])), 1:length(i[2]))
        if length(p)!=prod(x->factorial(length(x)),o) error() end
        for u in o
          s.prop[:permutable][i[1][u]]=fill(j,length(u))
          j+=1
        end
      end
      r[i[1]] = a[p[1]]
    end
  end
  p=inv(Perm(sortperm(r)))
  for i in [:mC, :charNumbers, :eigen, :span, :eps, :dims, :permutable]
    s.prop[i]^=p
  end
  s.prop[:translation]=filter(t->
       s.prop[:eigen]==map(i->E(predeigen(i,mod(i+t,e(s)))),1:e(s)),0:e(s)-1)
  if length(s.prop[:translation])==1 delete!(s.prop,:translation)
  else
    p=s.prop[:translation][2]
    if any(x->mod(x, p)!=0,s.prop[:translation]) || 
       length(s.prop[:translation])!=e(s)//p
        error()
    end
    quality=map(t->map(i->param(mod(t+i-1,e(s))+1,i),1:e(s)), s.prop[:translation])
    quality=map(x->conductor(values(prod(y->Mvp(:x)-y,x).d)),quality)
    quality=1+s.prop[:translation][filter(i->quality[i]==minimum(quality),
                                          eachindex(quality))]
    if isempty(quality) quality=[1] end
    m=HasType.Rotations(mC(s))
    m=findfirst(==(maximum(m[quality])),m)
    m=circshift(1:e(s),m-1)
    for i in [:mC,:charNumbers,:eigen,:span,:eps,:dims,:permutable]
      s.prop[i]=s.prop[i][m]
    end
    s.prop[:translation]=p
  end
  s.prop[:Hecke]=hecke(relative_group(s), improve_type([map(i->param(i,i),1:e(s))]))
  return s.prop[:Hecke]
end

CHEVIE[:relativeSeries]=true
function Hecke(s::Series)
  if haskey(s.prop,:Hecke) return s.prop[:Hecke] end
  if iscyclic(s) paramcyclic(s)
  elseif CHEVIE[:relativeSeries]
    InfoChevie("  # Relative: ",s,"\n")
    RelativeSeries(s)
  end
  if haskey(s.prop, :Hecke) return s.prop[:Hecke]
  else return nothing
  end
end

function AllProperSeries(W)
  l=sort(unique(conductor.(refleigen(W))))
  vcat(map(d->vcat(map(i->CuspidalSeries(W, d, i), 
                             1:length(relative_degrees(W,d)))...),l)...)
end

function findfractionalpowers(W)
  if length(gens(W))==0 return [0] end
  uc=UnipotentCharacters(W)
  uc[:fractions]=fill(nothing,length(uc))
  reasons = map(x->[], 1:length(uc))
  for h in uc[:harishChandra]
    if length(h.levi)!=length(gens(W))
      L=reflection_subgroup(W, h.levi)
      cusp=FindCuspidalInLevi(h[:cuspidalName], L)
      n=findfractionalpowers(L)
      if !isnothing(n[cusp])
        uc[:fractions][h[:charNumbers]]=h[:charNumbers]*0+n[cusp]
        for i in h[:charNumbers] push!(reasons[i], joindigits(h.levi)) end
      end
    end
  end
  sers=AllProperSeries(W)
  map(Hecke, sers)
  sers=gapSet(filter(x->haskey(x, :mC) && iscyclic(x),sers))
  while !isempty(sers)
    for s in sers
      p = findfirst(i->uc.fractions[i]!==nothing,charNumbers(s))
      if isnothing(p) print("reticent series", s, "\n")
      else
        frac = mC(s)[p]*e(s)*s.d
        fix = Mod1(uc[:fractions][charNumbers(s)[p]]-frac)
        if fix != 0
            print(" **** Badly normalized series ", s, " adding ", fix, "\n")
        end
        print("\n", s, "==>")
        for i in 1:e(s)
          cn=charNumbers(s)[i]
          print(cn,".c")
          frac = Mod1(mC(s)[i] * e(s) * s.d + fix)
          if isnothing(uc[:fractions][cn]) push!(reasons[cn], s.d)
            uc[:fractions][cn] = frac
          elseif uc[:fractions][cn] != frac
            print("Failed! ",cn, "==", TeXStrip(uc[:TeXCharNames][cn]),
                 " in ", s, "\n     where mC mod1==", frac, "\n
                 conflicts with ", reasons[cn], "\n")
          else push!(reasons[cn], s.d)
          end
        end
        sers = Difference(sers, [s])
      end
    end
  end
  for i = uc[:harishChandra]
    if length(i[:charNumbers])>1
      print(Format([uc[:fractions][i[:charNumbers]], 
                    map(length, reasons[i[:charNumbers]])]), "\n")
    else
      print(uc[:fractions][i[:charNumbers][1]],reasons[i[:charNumbers][1]],"\n")
    end
    p=gapSet(uc[:fractions][i[:charNumbers]])
    if haskey(i, :qEigen)
      if p!=[i[:qEigen]]
        if p == [nothing]
          print("!!!!!!!!!!!! for ", i[:charNumbers], 
                " qEigen should be nothing is ", i[:qEigen], "\n")
        else
          error("for HCseries $i of ",ReflectionName(W)," qEigen should be ",p)
          uc[:fractions][i[:charNumbers]] = i[:qEigen]+0*i[:charNumbers]
        end
      end
    elseif p!=[0] print("!!!!!!!!!!!!!! qEigen unbound should be ", p, "\n")
    end
  end
  return uc[:fractions]
end

# series of d-regular element w
function PrincipalSeries(W, d)
  s=split_levis(W, d, length(relative_degrees(W, d)))
  if length(s)!=1 error(" not one ", d, "-Sylow") end
  s=s[1]
  Series(W,s,findfirst(iszero,UnipotentCharacters(s).prop[:A]), d)
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

function ReflectionCosetFromType(arg...,)
  if length(arg) == 0 return CoxeterGroup() end
  if length(arg) > 1 error(ReflectionName(t), " non-irred not yet") end
  for t in arg
    g = ReflectionGroup(t[:orbit][1])
    if length(t[:orbit]) > 1
      g1 = g
      for i in 1:length(t[:orbit])-1 g*=g1 end
      i=reflrep(g,t[:twist]*Perm(vcat((Rotations(map(x->x[:indices],g[:type_])))[2]...)))
      g = Spets(g, i)
    else
      g = Spets(g, t[:twist])
    end
  end
  return g
end

SpetsFromType(t)=spets(reflection_group(t.orbit), t.twist)

# get Hecke of series s. If scalar-twisted untwist first
function getHecke(s)
# println("s=",s.spets,"\n")
  t=refltype(s.spets)
  if length(t)==1 && haskey(t[1],:scalar) && !all(isone,t[1].scalar)
    t=copy(t[1])
    scal=t.scalar
    InfoChevie("   # removing scal==", scal, "\n")
    t.scalar=map(x->1,t.scalar)
    if t.orbit[1].series=="B" && t.orbit[1].rank==2 && t.twist==perm"(1,2)"
      t.orbit[1].cartanType=E(8)-E(8, 3)
    end
    g=SpetsFromType(t)
    if isnothing(g) return nothing end
    W=Group(s.spets)
    p=Group(g)(word(W,s.levi.phi/s.spets.phi)...)
    inds=convert(Vector{Int},map(x->findfirst(==(Group(g)(word(W,x)...)),
                   reflections(Group(g))),gens(Group(s.levi))))
    l=subspets(g,inds,p)
    if length(scal)>1
      ChevieErr("scal==", scal, " unimplemented\n")
      return nothing
    end
    scal=Root1(scal[1])
# one should do an Ennola of the Levi's cuspidal by the absorbed part
# of scal by the center of the levi
    if !(s.cuspidal in cuspidal(UnipotentCharacters(l),s.d/scal))
      e = qqennola(Group(l))[1]
      c = abs(s.cuspidal^e)
      if !(c in cuspidal(UnipotentCharacters(l), s.d/scal))
        c=abs(s.cuspidal^(e^-1))
      end
    else c=s.cuspidal
    end
    s1=Series(g, l, c, (s.d/scal).r)
    p=getHecke(s1)
    if !isnothing(p)
      p=map(x->x(q=inv(E(scal))*Mvp(:q)),p)
    end
    return p
  else
    if any(u->haskey(u, :scalar) && !all(x->x==1,u.scalar), t)
      ChevieErr("scals==", map(u->u.scalar,t), " unimplemented\n")
#     return nothing
    end
    paramcyclic(s)
    haskey(s.prop,:Hecke) ? s.prop[:Hecke].para[1] : nothing
  end
end

function RelativeSeries(s)
# if !(haskey(s, :charNumbers)) CharNumbers(s) end
  WGL=relative_group(s)
  res=map(enumerate(WGL.prop[:reflists])) do (i,r)
 #  println("r=",r,"\nphi=",s.element/s.spets.phi)
    R=subspets(s.spets,r,s.element/s.spets.phi)
    l=inclusion(s.levi)
    if !issubset(l,inclusion(R))
      r=map(r)do w
        rr=roots(Group(s.spets))
        p=findfirst(y->!isnothing(ratio(rr[w],rr[y])),l)
        !isnothing(p) ? l[p] :  w
      end
      R=subspets(s.spets,r,s.element/s.spets.phi)
      if !issubset(r,inclusion(R))
        l=sort(unique(map(l)do w
          rr=roots(Group(s.spets))
          p=findfirst(y->!isnothing(ratio(rr[w],rr[y])),inclusion(Group(R)))
          !isnothing(p) ? inclusion(R,p) : w
        end))
        if !issubset(l,inclusion(R))
          error("could not change r==",r," so that Subspets(r) contains",l,"\n")
        end
      end
    end
#   println("R=",R," L=",subspets(R,l,s.element/R.phi))
#   println("R=",R," l=",l," s.element/R.phi=",s.element/R.phi)
    p=Series(R,subspets(R,restriction(R,l),s.element/R.phi),s.cuspidal,s.d.r)
    p.prop[:orbit]=simple_representatives(WGL)[i]
    r=Int.(indexin(sort(unique(reflections(WGL))),reflections(WGL)))
    if gens(WGL)==WGL.prop[:parentMap]
      r=filter(x->reflections(WGL)[x] in Group(p.spets),r)
    else
      w=map(x->word(WGL,reflections(WGL)[x]), r)
      w=map(x->prod(WGL.prop[:parentMap][x]), w)
#     xprint(w,"\n")
      r=r[map(w)do x
          g=Group(p.spets)
          if !(x isa Perm) return x.phi in g
          else return x in g
          end
         end]
    end
    p.prop[:WGL]=reflection_subgroup(WGL, r)
    p.prop[:WGLdims]=CharTable(relative_group(p)).irr[:,1]
    if all(isone,p.prop[:WGLdims]) p.prop[:e]=length(p.prop[:WGLdims]) end
    return p
  end
  s.prop[:relativeSpets]=map(x->x.spets, res)
  p=getHecke.(res)
  if nothing in p return res end
  s.prop[:Hecke]=hecke(WGL,improve_type(p))
# xprint("Hecke=",s.prop[:Hecke],"\n")
  u1=schur_elements(s.prop[:Hecke])//1
  if any(x->any(y->!all(isinteger, values(y.d)), keys(x.d)),u1)
    ChevieErr(s.prop[:Hecke], " fractional: wrong set of SchurElements")
    return res
  end
  u1=map(x->s.deg//CycPol(x), u1)
  degcusp=Uch.CycPolUnipotentDegrees(s.levi)[s.cuspidal]
  ud=map(x->x*sign(Int((x//degcusp)(E(s.d)))), 
                Uch.CycPolUnipotentDegrees(s.spets)[CharNumbers(s)])
  p=Perm(u1,ud)
# the permutation should also take in account eigenvalues
  if isnothing(p)
#   xprint("u1=",u1," ud=",ud,"\n")
#   println(sort(indexin(ud,u1)))
#   println(sort(indexin(u1,ud)))
    ChevieErr(s.prop[:Hecke], " wrong set of SchurElements")
    return res
  end
  s.prop[:charNumbers]^=p
  if haskey(s.prop,:span) s.prop[:span]^=p end
  aA=map(x->E(conductor(s.d)^2,valuation(x)+degree(x)),u1)
  p=position_regular_class(WGL,s.d)
  if p == false && length(WGL)==1 p=1 end
  o=map(x->x[p]//x[1], eachrow(CharTable(WGL).irr))
  s.prop[:predictedEigen] = map(i->aA[i]*o[i], 1:length(o)) *
    Uch.eigen(UnipotentCharacters(s.levi))[s.cuspidal]
  return res
end

function Checkzegen(W)
  l = AllProperSeries(W)
  print("with no Hecke:",filter(s->isnothing(Hecke(s)),l),"\n")
  print("    center==",OrderCenter(W),"\n")
  l=Filtered(l, (s->begin haskey(s, :Hecke) end))
  uc = UnipotentCharacters(W)
  aA = uc[:a] + uc[:A]
  for s = l
      e = gete(s[:Hecke])
      z = OrderCenter(relative_group(s))
      ucL = UnipotentCharacters(s.levi)
      aAL = ucL[:A] + ucL[:a]
      aAL = aAL[s.cuspidal]
      print(s, "\n")
      print(FormatTable([map(Mod1, (aA[charNumbers(s)] - aAL) // z),
                         map(Mod1, aA[charNumbers(s)] // z), map((x->begin
                              Mod1(1 // x) end), e)], 
   Dict{Symbol, Any}(:rowLabels => ["(aA-aAL)/z", "aA/z", "local Hecke irr"])))
  end
end

# checks minimal polynomials for Cyclotomic hecke algebras
function CheckRatCyc(s)
  if !(haskey(s, :Hecke)) return SPrint(FormatTeX(s), " Hecke  !  bound") end
  l = CollectBy(((s[:Hecke])[:parameter])[1], degree)
  fact = map((x->begin Product(x, (y->begin Mvp("T") - y end)) end), l)
  F = FieldOfDefinition(Group(s.spets))
end

function CheckRatSer(arg...)
  W = ApplyFunc(ComplexReflectionGroup, arg)
  s = AllProperSeries(W)
  s = Filtered(s, (x->begin x.spets != x.levi end))
  c = Join(map(CheckRatCyc, s), "\\hfill\\break\n")
  return c
end

function CheckPiGPiL(n)
  W = ComplexReflectionGroup(n)
  s = AllProperSeries(W)
  map(Hecke, s)
  s = Filtered(s, (ser->iscyclic(ser) && ser.principal))
  D0 = (s->begin Sum(ReflectionDegrees(Group(s.spets)) +
                  ReflectionCoDegrees(Group(s.spets))) -
              Sum(ReflectionDegrees(Group(s.levi)) +
                  ReflectionCoDegrees(Group(s.levi)))
          end)
  return map((x->begin D0(x) // e(x) end), s)
end

function Checkdovere(n)
  W=ComplexReflectionGroup(n)
  s=AllProperSeries(W)
  map(Hecke,s)
  s=filter(ser->iscyclic(ser) && ser.principal,s)
  filter(x->e(x)//x.d>2,s)
end

function CheckxiL(n)
  W=ComplexReflectionGroup(n)
  l=AllProperSeries(W)
  map(Hecke,l)
  l=filter(s->iscyclic(s) && s.principal,l)
  filter(x->!isone(Root1(PhiOnDiscriminant(x.levi))^x.d),l)
end

# check formula for product parameters
function CheckCoN(i)
  W=ComplexReflectionGroup(i)
  l=AllProperSeries(W)
  l=filter(x->length(x.levi)==1,l)
  map(Hecke, l)
  l=filter(x->haskey(x.prop,:Hecke),l)
  map(l)do s
    m=DeterminantMat(reflrep(W, s.element))
    sg=(-1)^sum(degrees(relative_group(s))-1)
    if length(hyperplane_orbits(relative_group(s)))>1 false
    else m*sg*prod(s.Hecke.para[1])^hyperplane_orbits(relative_group(s))[1].N_s
    end
  end
end

CheckLuCox(s)=[prod(x->1-x,filter(x->x!=1,s.Hecke.para[1])), s.deg(Mvp(:q))]

# data for linear char:
# returns d, e_W/d (mod d), \omega_c (mod d)
function datafor(W)
  l=map(d->PrincipalSeries(W, d),RegularEigenvalues(W))
  e=sum(degrees(W)+codegrees(W))
  map(s->[s.d, gcd(e*s.d,conductor(s.d)) // gcd(gcd(map(x->x[:N_s], 
                      HyperplaneOrbits(relative_group(s)))), denominator(s.d))], l)
end

# description of centralizers of regular elements
function cent(W)
  l = map((d->begin PrincipalSeries(W, d) end), RegularEigenvalues(W))
  res = map((s->begin [s.d, ReflectionName(relative_group(s))] end), l)
  l = gapSet(map((x->begin x[2] end), res))
  res = map((x->begin [map((z->begin z[1] end), Filtered(res, (y->begin
                                      y[2] == x end))), x] end), l)
  sort!(res)
  return map((x->begin SPrint(Join(x[1]), ":", x[2]) end), res)
end

function EigenspaceNumbers(W)
  d = ReflectionDegrees(W)
  if !(IsSpets(W)) return Union(map(divisors, d)) end
  e = map((p->begin Lcm(p[1], denominator(AsRootOfUnity(p[2]))) end), d)
  e = Union(map(divisors, e))
  return Filtered(e, (x->begin any((p->begin E(x, p[1]) == p[2] end), d) end))
end

# nrconjecture: si kd est le multiple maximum de d tq
# RelativeDegrees(W,kd)<>[], alors
# Length(RelativeDegrees(W,d))=Length(RelativeDegrees(W,kd))*Phi(kd)/Phi(d)
function nrconjecture(W)
  e = Difference(EigenspaceNumbers(W), RegularEigenvalues(W))
  e = Filtered(e, (x->begin !(x[1]) in r end))
  for d in e
    f = Filtered(r, (x->begin mod(x, d) == 0 end))
    if length(f) > 0
        print("**** failed: ", ReflectionName(W), d, f, "\n")
    end
    p = First(reverse(e), (i->begin mod(i[1], d) == 0 end))
    if p != d
        if length(RelativeDegrees(W, d)) * phi(d) != p[2] * phi(p[1])
            print("**** failed: ", ReflectionName(W), p, d, "\n")
        end
    end
  end
end

# An element which has a maximal ζ_d-eigenspace lives in the minimal
# parabolic subgroup such that the d-rank is the same
function minimaldparabolic(W,d)
  l=map(i->Difference(eachindex(gens(W)),[i]),eachindex(gens(W)))
  l=Filtered(l, (x-> length(RelativeDegrees(ReflectionSubgroup(W, x), d)) 
                 == length(RelativeDegrees(W, d))))
  map(x->ReflectionName(ReflectionSubgroup(W, x)),l)
end
