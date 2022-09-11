module Lusztig 
using ..Gapjm
export LusztigInductionTable, HCInductionTable
"""
`FindSeriesInParent(h,HF,WF,ww)`

find series `h` of levi `HF` in `WF` (`ww` is `uw.harish` or `uw.almostHarish`)
returns a namedtuple with 2 fields: ser= series found
  op: an element of `W` conjugating h.levi to ser.levi
      (not needed for Coxeter groups and standard levis H)
"""
function FindSeriesInParent(h,HF,WF,ww)
  if WF isa Spets W=Group(WF) else W=WF end
  n=sort(filter(!(x->isempty(x) || x=="Id"),split(h[:cuspidalName],"\\otimes ")))
  for w in ww
#  println("h[:levi]=",h[:levi]," w[:levi]=",w[:levi])
#  println("h[:eigenvalue]=",h[:eigenvalue]," w[:eigenvalue]=",w[:eigenvalue])
   if h[:eigenvalue]==w[:eigenvalue] && length(h[:levi])==length(w[:levi])
     p=sort(filter(!(x->isempty(x) || x=="Id"),split(w[:cuspidalName],"\\otimes ")))
     if p==n
       incHW=inclusion(HF,WF,h[:levi])
       p=transporting_elt(W,sort(incHW),sort(w[:levi]), 
                          action=(s,g)->sort(action.(Ref(W),s,g)))
       if !isnothing(p) return (ser=w, op=p) end
       p=transporting_elt(W, 
          sort(refls(reflection_subgroup(W,incHW))), 
          sort(refls(reflection_subgroup(W,w[:levi]))),
                          action=(s,g)->sort(s.^g))
       if !isnothing(p) return (ser=w, op=p) end
      end
    end
  end
  error("series ",h[:cuspidalName]," not found in ",WF)
end

# find the number of cuspidal of name n in UnipotentCharacters(HF)
function FindCuspidalInLevi(n,HF)
  strip=function(n,s)
    n=replace(n,Regex("^($s)*")=>"")
    replace(n,Regex("($s)*\$")=>"")
  end
  cusp=findfirst(==(strip(n,"\\\\otimes ")),
    map(x->strip(x,"\\\\otimes "),charnames(UnipotentCharacters(HF);TeX=true)))
  if isnothing(cusp) error("cuspidal ",n," not found in ",HF,"\n") end
  cusp
end

ChevieErr(x...)=println("*****ERROR: ",x...)
# l is a list of vectors each of length n. FindIntSol returns roots of unity
# x_i such that l[i]*[1,x2,..xn] is an integer for each i.
function FindIntSol(l)
  vars=vcat([1], map(i->Symbol("x$i"),2:length(l[1])))
# println("l=$l vars=$vars")
  l= map(v->sum(v.*Mvp.(vars)),l)
# println("l=$l")
  simplify=function()
    l=map(function(p)
      if iszero(p) || valuation(p)!=0 return p end
      c=coefficient(p,Monomial())
      if (c isa Cyc) 
        if conductor(c)>1 return p 
        else c=CyclotomicNumbers.num(c) end
      end
      if c<0 p=-p ; c=-c end
      return p-floor(c)
    end, l)
    l=unique!(sort(filter(!iszero,l)))
  end
  function usevals(v,val)
    for i in filter(j->haskey(vals,j) && vals[j] isa Mvp,vars)
      vals[i]=value(vals[i], v=>val)
      f=scalar(vals[i])
      if !isnothing(f) vals[i]=[f] end
    end
    l=map(x->value(x,v=>val),l)
    simplify()
  end
  # variable v is val: mvp
  function apply(v, val)
    if haskey(vals,v) return true end
    vals[v] = val
    usevals(v,val)
    return true
  end
  # variable v is val: list-of-rootsofunity-possibilities.
  function applylist(v, val)
   if haskey(vals,v) # vals[i] is a list of possibilities
      vals[v] = intersect(vals[v], val)
      if length(vals[v]) == 0 return false end
      val = vals[v]
    else
      vals[v] = val
    end
    if length(vals[v])>1 return true else val = vals[v][1] end
    usevals(v,val)
    return true
  end
  # if p has small norm must evaluate to 0
  smallnorm(p)=sum(abs,coefficients(p))<15//16
  simple = function (p)
    if length(p)>2 || length(p)==2 && valuation(p)>0 return false end
    if length(p) == 1
      if valuation(p) == 0
        InfoChevie("#I WARNING! FindIntSol cannot make:", p, " integral\n")
        return 0
      end
      v=(0, first(coefficients(p)), first(variables(first(monomials(p)))))
    else
      var=first(variables(p))
      v=(coefficient(p,Monomial()),coefficient(p,Monomial(var)),var)
    end
    N=conductor(collect(v[1:2]))
    if N%2!=0 N=2N end
    i=filter(j->isinteger(v[1]+v[2]*j),E.(N,0:N-1))
    if length(i) == 0
        InfoChevie("#I WARNING! FindIntSol cannot make:", p, " integral\n")
        return 0
    end
    return (v[3], i)
  end
  function try_()
    if length(l) == 0 return true end
    for i in eachindex(l)
        v = simple(l[i])
        if v != false
            if v == 0 return false
            else
                l = deleteat!(l, i)
                return applylist(v...)
            end
        end
    end
    i=findfirst(smallnorm,l)
    if !isnothing(i)
      v=l[i]
      for (m,c) in v.d
        if degree(m)>0 return apply(m.d.d[1][1],Mvp(m.d.d[1][1])-v//c) end
      end
    end
    stuck=true
    return false
  end
  vals=Dict{Symbol,Any}()
  stuck=false
  simplify()
  while try_() if isempty(l) return cartesian([1],map(v->vals[v],vars[2:end])...) end end
  if stuck InfoChevie("#I WARNING! FindIntSol: stuck ",l,"\n") end
  return false
end

function LusztigInductionPieces(LF,WF)
  L=Group(LF)
  W=Group(WF)
  if L isa CoxeterGroup && !issubset(inclusiongens(L),inclusiongens(W))
    w=standard_parabolic(W,L)
    L=L^w
    LF=spets(L,LF.phi^w)
  end
  uL=UnipotentCharacters(LF)
  uW=UnipotentCharacters(WF)
  hw=uW.almostHarishChandra
  map(uL.almostHarishChandra) do h
#   xprintln("LF=",LF," ",
#     map(x->haskey(x,:orbit) ? vcat(map(y->y.indices,x.orbit)...) : x.indices,
#     h[:relativeType]))
    (ser,op)=FindSeriesInParent(h,LF,WF,hw)
    if W isa CoxeterGroup
      WFGL=relative_coset(WF,inclusion(L,W,h[:levi]))
      classreps(WFGL)
      WGL=Group(WFGL)
      LFGL=relative_coset(LF,h[:levi])
      LFGL=subspets(WFGL,convert(Vector{Int},
                    map(x->findfirst(==(x),
                        inclusion(W,WGL.relativeIndices)),
                        inclusion(L,Group(LFGL).relativeIndices))),
                    WGL.MappingFromNormalizer(LF.phi*WF.phi^-1))
    else
      cL=reflection_subgroup(W,ser[:levi]) # cL^op is in LF
      Jb=vcat(map(x->haskey(x,:orbit) ? vcat(map(y->y.indices,x.orbit)...) :
                  x.indices,ser[:relativeType])...)
      WFGL=relative_coset(WF,ser[:levi],Jb)
      WGL=Group(WFGL)
      if !haskey(WGL,:MappingFromNormalizer) 
        @show WF, Jb
        @show ser[:levi]
        error("hahah!") 
      end
      rh=map(x->action(W,x,op),inclusion(L,W,vcat(map(
       x->haskey(x,:orbit) ? vcat(map(y->y.indices,x.orbit)...) : x.indices,
       h[:relativeType])...)))
#     xprintln("rh=",rh," op=",op)
      w=WGL.MappingFromNormalizer((LF.phi^op)*WF.phi^-1)
      if w==false error("Could not compute MappingFromNormalizer\n") end
      rh=map(rh)do x
        r=relative_root(W,cL,x)
        pr=PermX(WGL,reflectionmat(r.root,r.coroot))
        findfirst(i->refls(WGL,i)==pr,eachindex(roots(WGL)))
#       p=findfirst(i->refls(WGL,restriction(WGL,i))==
#            PermX(WGL,reflectionmat(r[:root],r[:coroot])),inclusion(WGL))
#       inclusion(WGL)[p] 
      end
#     xprintln("WFGL=",WFGL," rh=",rh," w=",w)
#     LFGL=subspets(WFGL,restriction(WFGL,rh),w)
      LFGL=subspets(WFGL,Vector{Int}(rh),w)
#     xprintln("LFGL=",LFGL)
#     ReflectionName(LFGL)
#     ReflectionName(WFGL)
    end
    lu=repr(LF;context=:TeX=>true)
    lg=repr(WF;context=:TeX=>true)
    lfgl=repr(LFGL;context=:TeX=>true)
    wfgl=repr(WFGL;context=:TeX=>true)
#   println(ser[:relativeType],h[:relativeType])
    InductionTable(conj.(InductionTable(LFGL,WFGL).scalar),
                   almostcharnames(rio(TeX=true),uW)[ser[:charNumbers]],
                   almostcharnames(rio(TeX=true),uL)[h[:charNumbers]],
                   isempty(h[:levi]) ? "piece from \$$lu\$ to \$$lg\$" :
     "piece from \$W_{$lu}($(join(h[:levi])),$(h[:cuspidalName]))\$=$lfgl"*
     "to \$W_{$lg}($(join(h[:levi])),$(h[:cuspidalName]))\$=$wfgl",
    Dict{Symbol,Any}(:repr=>"piece($(repr(LF)),$(repr(WF)))",
                     :wnum=>ser[:charNumbers], :hnum=>h[:charNumbers]))
  end
end

#CHEVIE.Cache.LusztigInductionMaps=false

"""
`LusztigInductionTable(R,W)`

`R`  should be a parabolic subgroup of the Coxeter group `W` or a parabolic
subcoset  of  the  Coxeter  coset  `W`,  in  each  case representing a Levi
subgroup  `𝐋` of  the algebraic  group `𝐆`  associated to `W`. The function
returns  an `InductionTable`  representing the  Lusztig induction ``R_𝐋^𝐆``
between unipotent characters.

```julia-repl
julia> W=coxgroup(:B,3)
B₃

julia> t=twistings(W,[1,3])
2-element Vector{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}}:
 B₃₍₁₃₎=Ã₁×A₁Φ₁
 B₃₍₁₃₎=Ã₁×A₁Φ₂

julia> LusztigInductionTable(t[2],W)
Lusztig Induction from B₃₍₁₃₎=Ã₁×A₁Φ₂ to B₃
     │11⊗ 11 11⊗ 2 2⊗ 11 2⊗ 2
─────┼────────────────────────
111. │     1    -1    -1    .
11.1 │    -1     .     1   -1
1.11 │     .     .    -1    .
.111 │    -1     .     .    .
21.  │     .     .     .    .
1.2  │     1    -1     .    1
2.1  │     .     1     .    .
.21  │     .     .     .    .
3.   │     .     .     .    1
.3   │     .     1     1   -1
B₂:2 │     .     .     1   -1
B₂:11│     1    -1     .    .
```
"""
function LusztigInductionTable(LF,WF;check=true)
  if !(WF isa Spets) WF=spets(WF) end
  if !(LF isa Spets) LF=spets(LF) end
  uW=UnipotentCharacters(WF)
  uL=UnipotentCharacters(LF)
  if isnothing(uL)||isnothing(uW) return nothing end
  lu=repr(LF;context=:TeX=>true)
  lg=repr(WF;context=:TeX=>true)
  res=InductionTable(fill(0,length(uW),length(uL)),
                     charnames(uW;TeX=true), charnames(uL;TeX=true),
    "Lusztig Induction from \$$lu\$ to \$$lg\$",
  Dict{Symbol,Any}(:repr=>"LusztigInductionTable($(repr(LF)),$(repr(WF)))"))
# res=CHEVIE[:GetCached](uW, "LusztigInductionMaps", res,
#       x->[inclusion(Group(x[:u]))[1:ngens(Group(x[:u]))],
#       x[:u][:phi]*x[:g][:phi]^-1])
# if haskey(res, :scalar) return res end
  res.pieces=LusztigInductionPieces(LF,WF)
  if res.pieces==false return nothing end
  fL=fourier(uL)
  hh=uL.almostHarishChandra
  fWinv=fourier(uW)'
  maps=map(res.pieces)do piece
    mapping=fill(zero(eltype(piece.scalar)),length(uW),length(uL))
    mapping[piece.wnum,piece.hnum]=piece.scalar
    fWinv*mapping*fL
  end
  ret=function(mapping)
    if mapping == false
      ChevieErr("Failed\n")
      return check ? nothing : res
    end
    res.scalar.=mapping
    return res
  end
  smap = sum(maps)
  if all(v->all(isinteger, v), smap)
    if LF.phi==WF.phi && !all(>=(0), Int.(smap))
      ChevieErr("non-positive RLG for untwisted L")
    end
    return ret(smap)
  end
  scalars=FindIntSol(unique(toL(hcat(vec.(maps)...))))
  if scalars == false return ret(scalars) end
  if length(scalars) > 1
    ChevieErr("#I WARNING: ambiguity in scalars:", scalars, "\n")
  end
  scalars = scalars[1]
  if any(x->x isa Mvp, scalars) error() end
  if any(!isone,scalars)
    res.scalars=scalars
    if all(isinteger, scalars) p="#I signs are "
    else InfoChevie("#I non-sign scalars needed:", scalars, "\n")
    end
  end
  scalars=sum(i->maps[i]*scalars[i],1:length(scalars))
  if LF.phi==WF.phi && !all(>=(0), Flat(scalars))
      ChevieErr("non-positive RLG for untwisted L")
  end
  return ret(scalars)
end

function HCInductionTable(HF, WF)
  if !(WF isa Spets) WF=spets(WF) end
  uw=UnipotentCharacters(WF)
  W=Group(WF)
  if !(HF isa Spets) HF=spets(HF) end
  uh=UnipotentCharacters(HF)
  H=Group(HF)
  lu=repr(HF;context=:TeX=>true)
  lg=repr(WF;context=:TeX=>true)
  res = InductionTable(fill(0, length(uw),length(uh)),
    charnames(uw;TeX=true), charnames(uh;TeX=true),
    "Harish-Chandra Induction from \$$lu\$ to \$$lg\$",
  Dict{Symbol,Any}(:repr=>"HCInductionTable($(repr(HF)),$(repr(WF)))"))
# res = CHEVIE[:GetCached](uw, "HCInductionMaps", res, (x->begin
#  Group(x[:u])[:rootInclusion][Group(x[:u])[:generatingReflections]] end))
# if haskey(res, :scalar) return res end
  res.pieces=map(uh.harishChandra)do h
    ser,op = FindSeriesInParent(h, HF, WF, uw.harishChandra)
    Jb = vcat(map(x->x.indices, ser[:relativeType])...)
#   println("Jb=$Jb")
    if Group(WF) isa CoxeterGroup
      if isempty(ser[:relativeType]) Wi=coxgroup()
      else Wi=rootdatum(cat(cartan.(ser[:relativeType])...;dims=(1,2)))
      end
      if !isone(op)
        rh=gens(relative_group(H, h[:levi]))
        rh=filter(a->order.(gens(Wi)[a])==order.(rh),
          arrangements(eachindex(gens(Wi)), length(rh)))
        if length(rh)>1 ChevieErr("WARNING: embedding ambiguous:",rh,"\n") end
        Hi=reflection_subgroup(Wi,rh[1])
      else
        Hi=reflection_subgroup(Wi,
             map(x->Int(findfirst(y->x in orbit(WF.phi,y),Jb)), 
        inclusion(H,W,vcat(map(x->x.indices,h[:relativeType])...))))
      end
    else
      L = reflection_subgroup(W, ser[:levi])
      Wi = relative_group(W, ser[:levi], Jb)
      if false
        if h[:levi]==[] rh=inclusiongens(H)
        else rh = filter(i->refls(parent(W),i) in H,Jb)
        end
      else
        rh=reflection_subgroup(W,restriction(W,inclusiongens(H).^op))
#       println("inclusion(rh)=",inclusion(rh)," inclusion(L)=",inclusion(L))
        if !issubset(inclusion(L),inclusion(rh)) error() end
        rh=filter(i->!(refls(parent(W),i) in L),inclusion(rh))
      end
      function getHi()
        sH=length(relative_group(H, h[:levi]))
        rr=Int[]
        for x in rh
          r=relative_root(W, L, x)
          p=findfirst(i->refls(Wi,restriction(Wi,i))==
                       PermX(Wi,reflectionmat(r.root,r.coroot)),inclusion(Wi))
          if p!==isnothing && !(p in rr) push!(rr,p) end
          Hi = reflection_subgroup(Wi, rr)
          if length(Hi) == sH return Hi end
        end
        reflection_subgroup(Wi, Int[])
      end
      Hi = getHi()
    end
    lu=repr(Hi;context=:TeX=>true)
    lg=repr(Wi;context=:TeX=>true)
    piece = InductionTable(InductionTable(Hi, Wi).scalar, 
     charnames(uw;TeX=true)[ser[:charNumbers]], 
     charnames(uh;TeX=true)[h[:charNumbers]],
    "Harish-Chandra Induction piece from \$$lu\$ to \$$lg\$",
    Dict{Symbol,Any}())
    res.scalar[ser[:charNumbers],h[:charNumbers]] = piece.scalar
    piece
  end
  return res
end
end
