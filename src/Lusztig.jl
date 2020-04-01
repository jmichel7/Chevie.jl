module Lusztig 
using ..Gapjm
export LusztigInductionTable
# find series h of levi HF in WF (sers may be .harish or .almostHarish)
# return a record with 2 fields: ser= found series
#   op: an element of W which conjugates h.levi to ser.levi
#       (never needed for Coxeter groups and standard levis H)
FindSeriesInParent=function(h,HF,WF,sers)
  if WF isa Spets W=Group(WF) else W=WF end
  n=sort(filter(!isempty,split(h[:cuspidalName],"\\otimes ")))
  for y in sers
#  println("h[:levi]=",h[:levi]," y[:levi]=",y[:levi])
#  println("h[:eigenvalue]=",h[:eigenvalue]," y[:eigenvalue]=",y[:eigenvalue])
   if h[:eigenvalue]==y[:eigenvalue] && length(h[:levi])==length(y[:levi])
     p=sort(filter(!isempty,split(y[:cuspidalName],"\\otimes ")))
#    print(" p=",p)
     if p==n
       p=representative_operation(W,
          sort(restriction(WF,inclusion(HF,h[:levi]))), 
          sort(y[:levi]), action=(s,g)->sort(s.^g))
       if !isnothing(p) return (ser=y, op=p) end
       p=representative_operation(W, 
          sort(reflections(reflection_subgroup(W, h[:levi]))), 
          sort(reflections(reflection_subgroup(W, y[:levi]))), 
                                  action=(s,g)->sort(s.^g))
       if !isnothing(p) return (ser=y, op=p) end
      end
    end
  end
  error("series ", n," not found in ",WF)
end

# find the number of cuspidal of name n in UnipotentCharacters(HF)
FindCuspidalInLevi=function(n,HF)
  strip=function(n,s)
    n=replace(n,Regex("^($s)*")=>"")
    replace(n,Regex("($s)*\$")=>"")
  end
  cusp=findfirst(==(strip(n,"\\\\otimes ")),
    map(x->strip(x,"\\\\otimes "),UnipotentCharacters(HF).prop[:TeXCharNames]))
  if isnothing(cusp) error("cuspidal ",n," not found in ",HF,"\n") end
  cusp
end

InfoChevie=print
ChevieErr=error
# l is a list of vectors each of length n. FindIntSol returns roots of unity
# x_i such that l[i]*[1,x2,..xn] is an integer for each i.
FindIntSol=function(l)
  vars=vcat([1], map(i->Symbol("x$i"),2:length(l[1])))
# println("l=$l vars=$vars")
  l= map(v->sum(v.*Mvp.(vars)),l)
# println("l=$l")
  simplify=function()
    l=map(function(p)
          if iszero(p) || valuation(p)!=0 return p end
      c=p.d[one(Monomial)]   
      if (c isa Cyc) 
        if c.n>1 return p 
        else c=Real(c) end
      end
      if c<0 p=-p ; c=-c end
      return p-floor(c)
    end, l)
    l=unique(sort(filter(!iszero,l)))
  end
  function usevals(v,val)
    for i in filter(j->haskey(vals,j) && vals[j] isa Mvp,vars)
      vals[i]=value(vals[i], v=>val)
      f=scal(vals[i])
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
  smallnorm(p)=sum(x->sqrt(Real(x*conj(x))),values(p.d))<15//16
  simple = function (p)
    if length(keys(p.d))>2 || length(keys(p.d))==2 && valuation(p)>0 return false end
    if length(keys(p.d)) == 1
      if valuation(p) == 0
        InfoChevie("#I WARNING! FindIntSol cannot make:", p, " integral\n")
        return 0
      end
      v=(0, values(p)[1], keys(p)[1].d.d[1])
    else
      var=variables(p)[1]
      v=(p.d[one(Monomial)],p.d[Monomial(var)],var)
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
  while try_() if isempty(l) return Cartesian([1],map(v->vals[v],vars[2:end])...) end end
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
#     printc( "LF=",LF," ",
#      map(x->haskey(x,:orbit) ? vcat(map(y->y.indices,x.orbit)...) : x.indices,
#      h[:relativeType]),"\n")
    (ser,op)=FindSeriesInParent(h,LF,WF,hw)
    cL=reflection_subgroup(W,ser[:levi])
    if Group(WF) isa CoxeterGroup
      WFGL=relative_coset(WF,restriction(W,inclusion(L,h[:levi])))
      WGL=Group(WFGL)
      LFGL=relative_coset(LF,h[:levi])
      LFGL=subspets(WFGL,convert(Vector{Int},inclusion(WGL,
                    map(x->findfirst(==(x),WGL.prop[:relativeIndices]),
                        inclusion(L,Group(LFGL).prop[:relativeIndices])))),
                    WGL.prop[:MappingFromNormalizer](LF.phi*WF.phi^-1))
    else
      Jb=vcat(map(x->haskey(x,:orbit) ? vcat(map(y->y.indices,x.orbit)...) :
                  x.indices,ser[:relativeType])...)
      WFGL=relative_coset(WF,ser[:levi],Jb)
      if isnothing(WFGL) return nothing end
      WGL=Group(WFGL)
      rh=restriction(W,inclusion(L,vcat(map(
       x->haskey(x,:orbit) ? vcat(map(y->y.indices,x.orbit)...) : x.indices,
       h[:relativeType])...))).^op
#     printc("rh=",rh," op=",op,"\n")
      w=WGL.prop[:MappingFromNormalizer]((LF.phi^op)*WF.phi^-1)
      if w==false error("Could not compute MappingFromNormalizer\n") end
      rh=map(function(x) r=GetRelativeRoot(W,cL,x)
              p=findfirst(i->reflection(WGL,restriction(WGL,i))==
              PermX(WGL,reflection(r[:root],r[:coroot])),inclusion(WGL))
              inclusion(WGL)[p] end,rh)
#     printc("WFGL=",WFGL," rh=",rh," w=",w,"\n")
      LFGL=subspets(WFGL,rh,w)
#     printc("LFGL=",LFGL,"\n")
#     ReflectionName(LFGL)
#     ReflectionName(WFGL)
    end
    lu=sprint(show,LF;context=:TeX=>true)
    lg=sprint(show,WF;context=:TeX=>true)
#   println(ser[:relativeType],h[:relativeType])
    InductionTable(conj.(InductionTable(LFGL,WFGL).scalar),
                   uW.prop[:almostTeXCharNames][ser[:charNumbers]],
                   uL.prop[:almostTeXCharNames][h[:charNumbers]],
                   isempty(h[:levi]) ? "piece from \$$lu\$ to \$$lg\$" :
     "piece from \$W_{$lu}($(join(h[:levi])),$(h[:cuspidalName]))\$=$LFGL"*
     "to \$W_{$lg}($(join(h[:levi])),$(h[:cuspidalName]))\$=$WFGL",
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
returns  an  `InductionTable`  representing  the  Lusztig induction `R_𝐋^𝐆`
between unipotent characters.

```julia-repl
julia> W=coxgroup(:B,3)
B₃

julia> t=twistings(W,[1,3])
2-element Array{Gapjm.Cosets.FCC{Int16,FiniteCoxeterSubGroup{Perm{Int16},Int64}},1}:
 Ã₁×A₁Φ₁
 Ã₁×A₁Φ₂

julia> LusztigInductionTable(t[2],W)
Lusztig Induction from Ã₁×A₁Φ₂ to B₃
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
  lu=sprint(show,LF;context=:TeX=>true)
  lg=sprint(show,WF;context=:TeX=>true)
  res=InductionTable(fill(0, length(uW),length(uL)),
    uW.prop[:TeXCharNames],uL.prop[:TeXCharNames],
    "Lusztig Induction from \$$lu\$ to \$$lg\$",
  Dict{Symbol,Any}(:repr=>"LusztigInductionTable($(repr(LF)),$(repr(WF)))"))
# res=CHEVIE[:GetCached](uW, "LusztigInductionMaps", res,
#       x->[inclusion(Group(x[:u]))[1:nbgens(Group(x[:u]))],
#       x[:u][:phi]*x[:g][:phi]^-1])
# if haskey(res, :scalar) return res end
  res.prop[:pieces]=LusztigInductionPieces(LF,WF)
  if res==false return false end
  fL=fourier(uL)
  hh=uL.almostHarishChandra
  fWinv=fourierinverse(uW)
  maps=map(res.prop[:pieces])do piece
    mapping=fill(zero(eltype(piece.scalar)),length(uW),length(uL))
    mapping[piece.prop[:wnum],piece.prop[:hnum]]=piece.scalar
    fWinv*mapping*fL
  end
  ret=function(mapping)
    if mapping == false
      ChevieErr("Failed\n")
      return check ? false : res
    end
    res.scalar.=mapping
    return res
  end
  smap = sum(maps)
  if all(v->all(isinteger, v), smap)
    if LF.phi==WF.phi && !all(x->x>=0, Int.(smap))
        ChevieErr("non-positive RLG for untwisted L")
    end
    return ret(smap)
  end
  scalars=FindIntSol(unique(toL(hcat(vec.(maps)...))))
  if scalars == false return ret(scalars) end
  if length(scalars) > 1
    ChevieErr("#I WARNING: ambiguity in scalars:", FormatGAP(scalars), "\n")
  end
  scalars = scalars[1]
  if any(x->x isa Mvp, scalars) error() end
  if any(!isone,scalars)
      res.prop[:scalars] = scalars
      if all(isinteger, scalars) p = "#I signs are "
      else InfoChevie("#I non-sign scalars needed:", scalars, "\n")
      end
  end
  scalars=sum(i->maps[i]*scalars[i],1:length(scalars))
  if LF.phi==WF.phi && !all(x->x>=0, Flat(scalars))
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
  lu=sprint(show,HF;context=:TeX=>true)
  lg=sprint(show,WF;context=:TeX=>true)
  res = InductionTable(fill(0, length(uw),length(uh)),
    uw.prop[:TeXCharNames],uh.prop[:TeXCharNames],
    "Harish-Chandra Induction from \$$lu\$ to \$$lg\$",
  Dict{Symbol,Any}(:repr=>"HCInductionTable($(repr(HF)),$(repr(WF)))"))
# res = CHEVIE[:GetCached](uw, "HCInductionMaps", res, (x->begin
#  Group(x[:u])[:rootInclusion][Group(x[:u])[:generatingReflections]] end))
# if haskey(res, :scalar) return res end
  res.prop[:pieces]=map(uh.harishChandra)do h
    ser,op = FindSeriesInParent(h, HF, WF, uw.harishChandra)
    Jb = vcat(map(x->x.indices, ser[:relativeType])...)
#   println("Jb=$Jb")
    if Group(WF) isa CoxeterGroup
      Wi=rootdatum(cat(cartan.(ser[:relativeType])...;dims=(1,2)))
      if !isone(op)
        rh=gens(relative_group(H, h[:levi]))
        rh=filter(a->order.(gens(Wi)[a])==order.(rh),
          arrangements(eachindex(gens(Wi)), length(rh)))
        if length(rh)>1 ChevieErr("WARNING: embedding ambiguous:",rh,"\n") end
        Hi=reflection_subgroup(Wi,rh[1])
      else
        Hi=reflection_subgroup(Wi,
             map(x->Int(findfirst(y->x in orbit(WF.phi,y),Jb)), 
       restriction(W,inclusion(H,vcat(map(x->x.indices,h[:relativeType])...)))))
      end
    else
      L = reflection_subgroup(W, ser[:levi])
      Wi = relative_group(W, ser[:levi], Jb)
      if false
        if h[:levi]==[] rh=inclusiongens(H)
        else rh = filter(i->reflection(parent(W),i) in H,Jb)
        end
      else
        rh=reflection_subgroup(W,inclusiongens(H).^op)
#       println("inclusion(rh)=",inclusion(rh)," inclusion(L)=",inclusion(L))
        if !issubset(inclusion(L),inclusion(rh)) error() end
        rh=filter(i->!(reflection(parent(W),i) in L),inclusion(rh))
      end
      function getHi()
        sH=length(relative_group(H, h[:levi]))
        rr=Int[]
        for x in rh
          r=GetRelativeRoot(W, L, x)
          p=findfirst(i->reflection(Wi,restriction(Wi,i))==
                       PermX(Wi,reflection(r[:root],r[:coroot])),inclusion(Wi))
          if !isnothing(p) 
            push!(rr,p) 
            rr=unique(rr)
          end
          Hi = reflection_subgroup(Wi, rr)
          if length(Hi) == sH return Hi end
        end
        return reflection_subgroup(Wi, Int[])
      end
      Hi = getHi()
    end
    lu=sprint(show,Hi;context=:TeX=>true)
    lg=sprint(show,Wi;context=:TeX=>true)
    piece = InductionTable(InductionTable(Hi, Wi).scalar, 
     uw.prop[:TeXCharNames][ser[:charNumbers]],
     uh.prop[:TeXCharNames][h[:charNumbers]],
    "Harish-Chandra Induction piece from \$$lu\$ to \$$lg\$",
    Dict{Symbol,Any}())
    res.scalar[ser[:charNumbers],h[:charNumbers]] = piece.scalar
    piece
  end
  return res
end
end
