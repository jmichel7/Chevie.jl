FindSeriesInParent=function(h,HF,WF,sers)
  if WF isa Spets W=Group(WF) else W=WF end
  n = sort(filter(!isempty,split(h[:cuspidalName],"\\otimes ")))
  for y = sers
   if h[:eigenvalue]== y[:eigenvalue] && length(h[:levi]) == length(y[:levi])
     p=sort(filter(!isempty,split(y[:cuspidalName],"\\otimes ")))
     if p==n
       p=representative_operation(W, HasType.gapSet(h[:levi]), 
                                     HasType.gapSet(y[:levi]), 
                                  action=(s,g)->HasType.gapSet(s.^g))
       if p!=false return (ser=y, op=p) end
       p=representative_operation(W, 
          HasType.gapSet(reflections(reflection_subgroup(W, h[:levi]))), 
          HasType.gapSet(reflections(reflection_subgroup(W, y[:levi]))), 
                                  action=(s,g)->HasType.gapSet(s.^g))
        if p != false return (ser=y, op=p) end
      end
    end
  end
  error("series not found")
end

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

FindIntSol=function(l)
  vars=vcat([1], map(i->"x$i",2:length(l[1])))
  simplify = function ()
    l = map(function(p)
      if p==0 || length(p[:elm][1][:elm])>0 || !IsRat(p[:coeff][1])
          return p
      end
      if p[:coeff][1]<0 p=-p end
      return p-Int(p[:coeff][1])
    end, l)
    l = gapSet(filter(!iszero,l))
  end
  l*= map(Mvp, vars)
  simplify()
  apply = function (v, val)
    i = Position(vars, v)
    if vals[i] !== nothing
        if IsMvp(val) return true end
        vals[i] = Intersection(vals[i], val)
        if length(vals[i]) == 0 return false
        elseif length(vals[i]) > 1 return true
        end
        val = vals[i]
    end
    vals[i] = val
    if IsList(vals[i])
        if length(vals[i]) > 1 return true
        else val = (vals[i])[1]
        end
    end
    for i = Filtered(1:length(vals), (j->begin
                        vals[j] !== nothing && IsMvp(vals[j])
                    end))
        vals[i] = Value(vals[i], [v, val])
        f = ScalMvp(vals[i])
        if f != false vals[i] = [f] end
    end
    l = map((x->begin Value(x, [v, val]) end), l)
    simplify()
    return true
  end
  smallnorm(p)=sum( x->GetRoot(Real(x*conj(x))),p[:coeff])<15//16
  simple = function (p)
    if length(p[:coeff])>2 || length(p[:coeff])==2 && length(p[:elm][1][:elm])>0
        return false
    end
    if length(p[:coeff]) == 1
      if length(((p[:elm])[1])[:elm]) == 0
        InfoChevie("#I WARNING! FindIntSol cannot make:", p, " integral\n")
        return 0
      end
      v = [0, p[:coeff][1], p[:elm][1][:elm][1]]
    else
      v = [p[:coeff][1], p[:coeff][2], p[:elm][2][:elm][1]]
    end
    N=NofCyc(v[1:2])
    if mod(N,2)!=0 N=2N end
    i=Filtered(map(i->E(N,i), 0:N-1), j->IsInt(v[1] + v[2] * j))
    if length(i) == 0
        InfoChevie("#I WARNING! FindIntSol cannot make:", p, " integral\n")
        return 0
    end
    return [v[3], i]
end
  function try_()
    if length(l) == 0 return true end
    for i = 1:length(l)
        v = simple(l[i])
        if v != false
            if v == 0 return false
            else
                l = Drop(l, i)
                return apply(v[1], v[2])
            end
        end
    end
    i = PositionProperty(l, smallnorm)
    if i != false
      v = l[i]
      i = PositionProperty(v[:elm], (x->begin length(x[:elm]) > 0 end))
      return apply(v[:elm][i][:elm][1],Mvp(v[:elm][i][:elm][1])-v//v[:coeff][i])
    end
    stuck = true
    return false
  end
  vals = [[1]]
  stuck = false
  while try_() if isempty(l) return Cartesian(vals) end end
  if stuck InfoChevie("#I WARNING! FindIntSol: stuck ",l,"\n") end
  return false
end

function LusztigInductionPieces(res)
  WF=res[:g]
  LF=res[:u]
  L=Group(LF)
  W=Group(WF)
  if L isa CoxeterGroup && !issubset(inclusiongens(L),inclusiongens(W))
    w=standard_parabolic(W,L)
    L=L^w
    LF=spets(L,LF.phi^w)
  end
  uL=UnipotentCharacters(LF)
  uW=UnipotentCharacters(WF)
  if isnothing(uL)||isnothing(uW) return nothing end
  res[:gNames]=(t,option)->CharNames(uW,option)
  res[:uNames]=(t,option)->CharNames(uL,option)
  res[:pieces]=[]
  hw=uW.almostHarishChandra
  for h in uL.almostHarishChandra  
    (ser,op)=FindSeriesInParent(h,LF,WF,hw)
    L=reflection_subgroup(W,ser.levi)
    if Group(WF) isa CoxeterGroup
      WFGL=RelativeCoset(WF,h.levi)
      LFGL=RelativeCoset(LF,h.levi)
      LFGL=CoxeterSubCoset(WFGL,inclusion(Group(WFGL))[
      List(Group(LFGL).relativeIndices,x->Position(Group(WFGL).relativeIndices,x))],
            Group(WFGL).MappingFromNormalizer(LF.phi*WF.phi^-1))
    else
      Jb=vcat(map(x->haskey(x,:orbit) ? vcat(List(x.orbit,y->y.indices)...) :
                  x.indices,ser.relativeType)...)
      WFGL=RelativeCoset(WF,ser.levi,Jb)
      if WFGL==false   return false end
      WGL=Group(WFGL)
      rh=OnTuples(Concatenation(List(h.relativeType,
       function(x)if IsBound(x.orbit)  return Concatenation(List(x.orbit,y->y.indices))
        else return x.indices end
end)),op)
      w=WGL.MappingFromNormalizer((LF.phi^op)*WF.phi^-1)
      if w==false Error("Could not compute MappingFromNormalizer\n") end
      rh=List(rh,function(x)local r
         r=GetRelativeRoot(W,L,x)
         return First(WGL.rootInclusion,
           i->Reflection(WGL,WGL.rootRestriction[i])==
            PermMatX(WGL,Reflection(r.root,r.coroot)))
      end)
      LFGL=SubSpets(WFGL,rh,w)
      ReflectionName(LFGL)
      ReflectionName(WFGL)
    end
    piece=rec(wnum=ser.charNumbers,
               hnum=h.charNumbers,
               uNames=function(t,option)
                 return AlmostCharNames(uL,option){t.hnum} end,
               gNames=function(t,option)
                 return AlmostCharNames(uW,option){t.wnum} end,
                scalar=ComplexConjugate(InductionTable(LFGL,WFGL).scalar),
               u=LFGL,
               g=WFGL,
               what="piece",
               operations=InductionTableOps)
    if h.levi==[]   piece.head=function(t,option)
      if IsBound(option.TeX)  
       return SPrint("from \$",ReflectionName(LF,option),
                     "\$ to \$",ReflectionName(WF,option),"\$")

      else return SPrint("from ",ReflectionName(LF,option),
                     " to ",ReflectionName(WF,option)) end
    end
    else piece.head=function(t,option)local sub,a,b
      if IsBound(option.TeX)   sub=L->SPrint("_{",L,"}")

                             else sub=L->SPrint("_[",L,"]") end

      a=SPrint("W",sub(ReflectionName(LF,option)),"(",sub(Join(h.levi)),",",
       TeXStrip(h.cuspidalName,option),")==",ReflectionName(t.u,option))

      b=SPrint("W",sub(ReflectionName(WF,option)),"(",sub(Join(ser.levi)),",",
       TeXStrip(ser.cuspidalName,option),")==",ReflectionName(t.g,option))

      if IsBound(option.TeX)   return SPrint("from \$",a,"\$\nto \$",b,"\$")

      else return SPrint("from ",a,"\nto   ",b) end end end
    Add(res.pieces,piece) end
  return res
end

#CHEVIE.Cache.LusztigInductionMaps=false

"""
`LusztigInductionTable(<R>,<W>)`

<R>  should be a parabolic subgroup of the Coxeter group <W> or a parabolic
subcoset  of  the  Coxeter  coset  <W>,  in  each  case representing a Levi
subgroup  `𝐋` of  the algebraic  group `𝐆`  associated to <W>. The function
returns  a  table  (modeled  after  'InductionTable', see "InductionTable")
representing the Lusztig induction `R_𝐋^𝐆` between unipotent characters.

|    gap> W:=CoxeterGroup("B",3);;
    gap> t:=Twistings(W,[1,3]);
    [ ~A1xA1<3>.(q-1), ~A1xA1<3>.(q+1) ]
    gap> Display(LusztigInductionTable(t[2],W));
    Lusztig Induction from ~A1xA1<3>.(q+1) to B3
          |'|'|11,11 11,2 2,11 2,2
    ___________________________
    111.  |'|'|    1   -1   -1   .
    11.1  |'|'|   -1    .    1  -1
    1.11  |'|'|    .    .   -1   .
    .111  |'|'|   -1    .    .   .
    21.   |'|'|    .    .    .   .
    1.2   |'|'|    1   -1    .   1
    2.1   |'|'|    .    1    .   .
    .21   |'|'|    .    .    .   .
    3.    |'|'|    .    .    .   1
    .3    |'|'|    .    1    1  -1
    B2:2  |'|'|    .    .    1  -1
    B2:11 |'|'|    1   -1    .   .|
"""
LusztigInductionTable = function (LF,WF;check=true)
  if !(WF isa Spets) WF=Spets(WF) end
  if !(LF isa Spets) LF=Spets(LF) end
  res = Dict{Symbol, Any}(:u => LF, :g => WF, :what => "LusztigInduction", 
  :head => function (t, option)
    if haskey(option,:TeX) math=x->"\$"*x*"\$" else math=x->x end
    return SPrint("Lusztig Induction from ",math(ReflectionName(t[:u], option)),
                  " to ", math(ReflectionName(t[:g], option)))
    end)
  uW=UnipotentCharacters(WF)
# res=CHEVIE[:GetCached](uW, "LusztigInductionMaps", res,
#       x->[inclusion(Group(x[:u]))[1:nbgens(Group(x[:u]))],
#       x[:u][:phi]*x[:g][:phi]^-1])
# if haskey(res, :scalar) return res end
  res = LusztigInductionPieces(res)
  if res == false return false end
  uL = UnipotentCharacters(LF)
  fL = Fourier(uL)
  hh = uL[:almostHarishChandra]
  fWinv = (UnipotentCharactersOps[:FourierInverse])(uW)
  maps = map(function (i)
    mapping= map(y->fill(0, length(uL)),1:length(uW))
    piece = res[:pieces][i]
    mapping[piece[:wnum]][piece[:hnum]] = piece[:scalar]
    return fWinv * mapping* fL
  end, 1:length(hh))
  ret=function (map)
    if map == false
      ChevieErr("Failed\n")
      if check return false else return res end
    end
    res[:scalar] = map
    return res
  end
  map = Sum(maps)
  if all(v->all(IsInt, v), map)
    if LF[:phi] == WF[:phi] && !(all(x->x>=0, Flat(map)))
        ChevieErr("non-positive RLG for untwisted L")
    end
    return ret(map)
  end
  scalars = FindIntSol(gapSet(TransposedMat(map(Concatenation, maps))))
  if scalars == false return ret(scalars) end
  if length(scalars) > 1
      ChevieErr("#I WARNING: ambiguity in scalars:", FormatGAP(scalars), "\n")
  end
  scalars = scalars[1]
  if any(IsMvp, scalars) error() end
  if any((x->begin x != 1 end), scalars)
      res[:scalars] = scalars
      if all(IsInt, scalars)
          p = "#I signs are "
      else
          InfoChevie("#I non-sign scalars needed:", FormatGAP(scalars), "\n")
      end
  end
  scalars = Sum(1:length(scalars), (i->begin maps[i] * scalars[i] end))
  if LF[:phi] == WF[:phi] && !(all((x->begin x >= 0 end), Flat(scalars)))
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
    if Group(WF) isa CoxeterGroup
      Wi=rootdatum(cat(cartan.(ser[:relativeType])...;dims=(1,2)))
      if !isone(op)
        rh=RelativeGroup(H, h[:levi])[:generators]
        rh=filter(a->order.(gens(Wi)[a])==order.(rh),
          arrangements(eachindex(gens(Wi)), length(rh)))
        if length(rh)>1 ChevieErr("WARNING: embedding ambiguous:",rh,"\n") end
        Hi=reflection_subgroup(Wi,rh[1])
      else
        Hi=reflection_subgroup(Wi,
             map(x->Int(findfirst(y->x in orbit(WF.phi,y),Jb)), 
              vcat(map(x->x.indices,h[:relativeType])...)))
      end
    else
      L = reflection_subgroup(W, ser[:levi])
      Wi = RelativeGroup(W, ser[:levi], Jb)
      if false
        if h[:levi]==[] rh=inclusiongens(H)
        else rh = filter(i->Reflection(Parent(W),i) in H,Jb)
        end
      else
        rh=reflection_subgroup(W,inclusiongens(H).^op)
        if !IsSubset(inclusion(L), inclusion(rh)) error() end
        rh=filter(i->!(Reflection(Parent(W),i) in L),inclusion(rh))
      end
      function getHi()
        sH = length(RelativeGroup(H, h[:levi]))
        rr = []
        for x = rh
          r = GetRelativeRoot(W, L, x)
          p = PositionProperty(Wi[:rootInclusion], i->
                               Reflection(Wi, restriction(Wi,i))==
                               PermMatX(Wi, Reflection(r[:root], r[:coroot])))
          if p != false AddSet(rr, p) end
          Hi = ReflectionSubgroup(Wi, rr)
          if Size(Hi) == sH return Hi end
        end
        return reflection_subgroup(Wi, [])
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
