# replacements for non-usable parts of chevie's tbl

CharTableSymmetric=Dict(:centralizers=>[
     function(n,pp) res=k=1;last=0
        for p in pp
          res*=p
          if p==last k+=1;res*=k
          else k=1
          end
          last=p
        end
        res
     end])

"""
 EvalPolRoot(pol, x, n, p) compute pol(p*x^(1/n))
  
  The point of this routine is to avoid unnecessary root extractions
  during evaluation (e.g., if pol has no terms of odd degree and n=2,
  then no root extraction is necessary).
"""
function EvalPolRoot(pol::Pol,x,n,p)
  if isempty(pol.c) return 0 end
  P=vcat(fill(0,mod(pol.v,n)),pol.c)
  P=map(i->Pol(P[i:n:length(P)],div(pol.v-mod(pol.v,n),n))(x*p^n),1:n)
  j=findlast(!iszero,P)
  if isnothing(j) return 0 end
  pol=Pol(P[1:j],0)
  l=pol.v-1+filter(i->!iszero(pol.c[i]),eachindex(pol.c))
  push!(l,n)
  r=gcd(l...)
  pol=Pol(pol.c[1:r:length(pol.c)],div(pol.v,r))
  pol(GetRoot(x,div(n,r))*p^r)
end

function VcycSchurElement(para,r,data=nothing)
  n=length(para)
  if !isnothing(data) para=para[data[:order]] else para = copy(para) end
  monomial(mon)=prod(map(^,para,Int.(mon[1:n])))
  if haskey(r, :rootUnity) && haskey(r,:root) error("cannot have both") end
  if haskey(r, :coeff) res = r[:coeff] else res = 1 end
  if haskey(r, :factor) res*=monomial(r[:factor]) end
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

function ImprimitiveCuspidalName(S)
  r=ranksymbol(convert(Vector{Vector{Int}},S))
  d=length(S)
  s=joindigits(map(length,S))
  if r==0 return "" end
  if sum(length,S)%d==1 # G(d,1,r)
    if r==1 return d==3 ? "Z_3" : "Z_{$d}^{$s}"
    else return "G_{$d,1,$r}^{$s}"
    end
  else # G(d,d,r)
    if r==2
      if d==4 return "B_2"
      elseif d==6 
        p=Dict("212010"=>"-1","221001"=>"1",
               "211200"=>"\\zeta^2","220110"=>"\\zeta_3")
        return "G_2[$(p[s])]"
      else p=CHEVIE.R("SymbolToParameter","I")(S);
	return "I_2($d)",FormatGAP(p)
      end
      elseif r==3 && d==3 
        return "G_{3,3,3}[\\zeta_3"* (s=="300" ? "" : "^2")*"]"
      elseif r==3 && d==4 
        return "G_{4,4,3}["* (s=="3010" ? "" : "-")*"\\zeta_4]"
    else return "G_{$d,$d,$r}^{$s}"
    end
  end
end

function BDSymbols(n,d)
  n-=div(d^2,4)
  if n<0 return Vector{Vector{Int}}[] end
  if d>0 return map(x->symbol_partition_tuple(x,d),partition_tuples(n,2)) end
   return map(chevieget(:D,:symbolcharparam),
              chevieget(:imp,:CharInfo)(2,2,n)[:charparams])
end
# an addition
chevieset(["A","B","D"],:EigenvaluesGeneratingReflections,(t->r->fill(1//2,r)))
# replacements for some functions in tables (whose automatic translation failed)
chevieset(["G25","G26","G29","G31","G32","G34"],:CartanMat,
  function(t)
    r=chevieget(t,:GeneratingRoots)
    e=chevieget(t,:EigenvaluesGeneratingReflections)
    e=map(x->E(denominator(x),numerator(x)),e)
    e=map(i->e[i]*conj(r[i])//(sum(conj(r[i]).*r[i])),eachindex(e))
    map(x->map(y->x*y,r),e)
  end)

chevieset(:D,:CartanMat,n->toL(cartan(:D,n)))
chevieset(:D,:CharTable,n->chevieget(:imp,:CharTable)(2,2,n))
chevieset(:B,:CharTable,n->chevieget(:imp,:CharTable)(2,1,n))
chevieset(:A,:CharTable,function(n)
  ct=chevieget(:imp,:CharTable)(1,1,n+1)
  ct[:irredinfo]=map(x->Dict(:charname=>joindigits(x)),chevieget(:A,:CharInfo)(n)[:charparams])
  ct
 end)

chevieset(Symbol("2A"),:CharTable,function(r)
  tbl = chevieget(:A, :CharTable)(r)
  tbl[:identifier] = "W(^2A_$r)"
  A=chevieget(:A,:LowestPowerFakeDegree).(chevieget(:A,:CharInfo)(r)[:charparams])
  tbl[:irreducibles]= (-1).^A .* tbl[:irreducibles]
  merge!(tbl, chevieget(Symbol("2A"), :ClassInfo)(r))
end)

chevieset(Symbol("2A"),:FakeDegree,function(n,p,q)
  res=chevieget(:A, :FakeDegree)(n, p, Pol())
  (-1)^valuation(res)*res(-q)
end)

chevieset(:A,:HeckeCharTable,
   function(n,para,root)
    if n==1 Dict(:irreducibles=>[[1,para[1][2]],[1,para[1][1]]],
                  :charnames=>["11","2"],
                  :classnames=>["11","2"],
                  :centralizers=>[2,2],
                  :identifier=>"H(A_1)")
     else chevieget(:imp,:HeckeCharTable)(1,1,n+1,para,root)
     end
   end)

chevieset(Symbol("2A"),:HeckeCharTable, function (r, param, rootparam)
  q = -param[1][1] // param[1][2]
  if isnothing(rootparam[1]) v=GetRoot(q,2,"CharTable(Hecke(2A))")
  else v = rootparam[1]
  end
  W=coxgroup(:A,r)
  qE=central_monomials(hecke(W,v))
  H=hecke(W, v^-2)
  T=Tbasis(H)
  tbl=Dict{Symbol,Any}(:identifier=>"H(^2A_$r)")
  cli=chevieget(Symbol("2A"),:ClassInfo)(r)
  cl=map(x->T(W(x...)*longest(W)), cli[:classtext])
  merge(tbl,cli)
  ct=CharTable(H)
  tbl[:charnames]=ct.charnames
  tbl[:classnames]=ct.classnames
  tbl[:centralizers]=ct.centralizers
  tbl[:irreducibles]=toL(permutedims(toM(map(char_values, cl))))
  A=chevieget(:A,:LowestPowerFakeDegree).(chevieget(:A,:CharInfo)(r)[:charparams])
  A=(-1).^A .* qE
  tbl[:irreducibles]=map((x,y)->x.*y,A,tbl[:irreducibles])
  CHEVIE[:compat][:AdjustHeckeCharTable](tbl, param)
  tbl
end)

chevieset(:A,:FakeDegree,(n,p,q)->fegsymbol([βset(p)])(q))
chevieset(:B,:HeckeCharTable,(n,para,root)->chevieget(:imp,:HeckeCharTable)(2,1,n,para,root))
chevieset(:D,:HeckeCharTable,(n,para,root)->chevieget(:imp,:HeckeCharTable)(2,2,n,para,root))
chevieset(:imp,:PowerMaps,function(p,q,r)
  if q!=1
    InfoChevie("# power maps not implemented for G($p,$q,$r)\n")
    return [[1],[1],[1]]
  end
  function pow(p,n)
    e=length(p)
    rr=map(x->[],1:e)
    for k in 1:e
      for l in p[k]
        g=gcd(n,l)
        for j in 1:g push!(rr[1+mod(div(n*(k-1),g),e)], div(l,g)) end
      end
    end
    for k in 1:e rr[k]=sort(rr[k],rev=true) end
    return rr
  end
  pp = chevieget(:imp, :ClassInfo)(p,q,r)[:classparams]
  l=keys(factor(factorial(r)*p))
  res=fill(Int[],maximum(l))
  for pw in l
    res[pw] = map(x->Position(pp, pow(x, pw)), pp)
  end
  res
 end)

chevieset(Symbol("2D"),:HeckeCharTable,
function (l,param,rootparam)
  q=-param[1][1]//param[1][2]
  q=vcat([[q^0,-1]],map(i->[q,-1],2:l))
  hi=chevieget(:B,:HeckeCharTable)(l,q,Int[])
  chr=1:length(hi[:classparams])
  lst=filter(i->isodd(length(hi[:classparams][i][2])),chr);
  tbl=Dict{Symbol,Any}(:identifier=>"H(^2D$l)",
                       :size=>div(hi[:size],2),
                       :orders=>hi[:orders][lst],
                       :centralizers=>div.(hi[:centralizers][lst],2),
                       :classes=>hi[:classes][lst],
	   :text=>"extracted from generic character table of HeckeB")
  merge!(tbl,chevieget(Symbol("2D"),:ClassInfo)(l))
  para=chevieget(Symbol("2D"),:CharParams)(l)
  tbl[:irredinfo]=map(i->Dict{Symbol,Any}(:charparam=>para[i],
      :charname=>chevieget(Symbol("2D"),:CharName)(l,para[i])),eachindex(para))
  para=chevieget(:B,:CharParams)(l)
  chr=filter(i->chevieget(Symbol("2D"),:testchar)(para[i]),chr)
  tbl[:irreducibles]=toL(permutedims(toM(map(x->char_values(
   Tbasis(hecke(coxgroup(:B,l),q))(vcat([1],Replace(x,[1],[1,2,1]))...),
   toM(hi[:irreducibles][chr])),tbl[:classtext]))))
   CHEVIE[:compat][:AdjustHeckeCharTable](tbl,param)
  tbl
end)

chevieset(Symbol("2D"),:CharTable,l->chevieget(Symbol("2D"),:HeckeCharTable)(l,fill([1,-1],l),fill(1,l)))

chevieset(:imp,:GeneratingRoots,function(p,q,r)
  if q==1 roots=[vcat([1],fill(0,r-1))]
  else
    if q!=p roots=[vcat([Cyc(1)],fill(0,r-1))] end
    v=vcat([-E(p),1],fill(0,r-2))
    if r==2 && q>1 && q%2==1 v*=E(p) end
    if q==p roots=[v] else push!(roots, v) end
  end
  for i=2:r
    v=fill(0,r)
    v[[i-1,i]]=[-1,1]
    push!(roots, v)
  end
  return roots
end)

chevieset(Symbol("2A"), :UnipotentClasses, function (r, p)
  uc=deepcopy(chevieget(:A, :UnipotentClasses)(r, p))
  for c in uc[:classes]
    t=parent(c[:red])
    if length(refltype(t))>1 error() end
    if length(refltype(t))==0 || rank(t)==1 WF=spets(parent(c[:red]))
    else WF=spets(parent(c[:red]),prod(i->Perm(i,rank(t)+1-i),1:div(rank(t),2)))
    end
    t=twistings(WF,inclusiongens(c[:red]))
    m=map(x->refleigen(x,position_class(x,x())),t)
    m=map(x->count(y->y==Root1(-1),x),m)
    p=findfirst(==(maximum(m)),m)
    c[:F]=t[p]()
  end
  uc
end)

chevieset(:B,:UnipotentClasses,function(r,char,ctype)
  part2dynkin=function(part)
    p=sort(reduce(vcat,map(d->1-d:2:d-1, part)))
    p=p[div(3+length(p),2):end]
    res= ctype==1 ? [2*p[1]] : [p[1]]
    append!(res,p[2:end]-p[1:end-1])
  end
  addSpringer1=function(s,cc)
    ss=First(uc[:springerSeries],x->x[:defect]==defectsymbol(s[:symbol]))
    if s[:sp] == [[], []] p = 1
    elseif s[:sp] == [[1], []] p = 2
    elseif s[:sp] == [[], [1]] p = 1
    else p = Position(CharParams(ss[:relgroup]), [s[:sp]])
    end
    ss[:locsys][p] = [length(uc[:classes]), Position(CharParams(cc[:Au]),
      map(x->x ? [1, 1] : [2], s[:Au]))]
  end
  if ctype==ER(2)
    ctype=2
    char=2
  end
  if char==2 ss=XSP(4,2,r)
  elseif ctype==1 ss=XSP(2,1,r)
  else ss=XSP(2,0,r)
  end
  l = union(map(c->map(x->[defectsymbol(x[:symbol]), sum(sum,x[:sp])], c), ss))
  sort!(l,by=x->[abs(x[1]),-sign(x[1])])
  uc = Dict{Symbol, Any}(:classes=>[])
  uc[:springerSeries]=map(l)do d
    res=Dict(:relgroup=>coxgroup(:C,d[2]),:defect=>d[1],:levi=>1:r-d[2])
    res[:locsys]=fill([0, 0],NrConjugacyClasses(res[:relgroup]))
    if char==2 res[:Z]=[1]
    elseif ctype==1 res[:Z]=[(-1)^(r-d[2])]
    elseif conductor(ER(2*(r-d[2])+1))==1 res[:Z]=[1]
    else res[:Z]=[-1]
    end
    res
  end
  if char != 2
    symbol2para = function(S)
      c=sort(reduce(vcat,S))
      i=1
      part=Int[]
      d=mod(ctype,2)
      while i<=length(c)
        l=2*(c[i]-(i-1))-d
        if i==length(c) || c[i+1]>c[i]
          push!(part, l+1)
          i+=1
        else
          append!(part,[l,l])
          i+=2
        end
      end
      reverse(filter(!iszero,sort(part)))
    end
  else
    symbol2para = function (S,)
      c=sort(reduce(vcat,S))
      i=1
      part=Int[]
      ex=Int[]
      while i<=length(c)
        l=2*(c[i]-2*(i-1))
        if i==length(c) || c[i+1]>c[i]+1
          push!(part, l)
          i+=1
        elseif c[i]==c[i+1]
          append!(part, [l-2, l-2])
          push!(ex,l-2)
          i+=2
        elseif c[i]+1==c[i+1]
          append!(part, [l-1, l-1])
          i+=2
        end
      end
      [reverse(filter(y->y!=0,sort(part))), ex]
    end
  end
  if char==2 ctype=1 end
  for cl = ss
    cc = Dict{Symbol, Any}(:parameter => symbol2para((cl[1])[:symbol]))
    cc[:Au] = CoxeterGroup(Concatenation(map(x->["A",1], cl[1][:Au]))...)
    if char != 2
      cc[:dynkin] = part2dynkin(cc[:parameter])
      cc[:name] = joindigits(cc[:parameter])
    else
      ctype = 1
      cc[:dimBu] = (cl[1])[:dimBu]
      cc[:name] = join(map(function(x)
        res=joindigits(fill(0,max(0,x[2]))+x[1],"[]")
        if x[1] in cc[:parameter][2] return string("(", res, ")") end
        res
      end, reverse(Collected(cc[:parameter][1]))), "")
    end
    cc[:red] = coxgroup()
    if char == 2 j = cc[:parameter][1]
    else j = cc[:parameter]
    end
    for j in tally(j)
      if mod(j[1], 2) == mod(ctype, 2)
        cc[:red] = cc[:red] * coxgroup(:C, div(j[2],2))
      elseif mod(j[2], 2) != 0
        if j[2]>1 cc[:red]*=coxgroup(:B, div(j[2]-1,2)) end
       elseif j[2]>2 cc[:red]*=coxgroup(:D, div(j[2],2))
      else cc[:red]*=torus(1)
      end
    end
    push!(uc[:classes], cc)
    for s in cl addSpringer1(s,cc) end
  end
  uc[:orderClasses] = hasse(Poset(map(x->
    map(function(y)
      if char!=2 return dominates(y[:parameter], x[:parameter]) end
      m=maximum(x[:parameter][1][1], y[:parameter][1][1])
      f=x-> map(i->sum(filter(<(i),x))+i*count(>=(i),x) ,1:m)
      fx=f(x[:parameter][1])
      fy=f(y[:parameter][1])
      for i in 1:m
        if fx[i]<fy[i] return false
        elseif fx[i]==fy[i] && i in y[:parameter][2]
          if i in Difference(x[:parameter][1],x[:parameter][2]) return false end
          if i<m && mod(fx[i+1]-fy[i+1],2)==1 return false end
        end
      end
      return true
    end, uc[:classes]), uc[:classes])))
  if char!=2 && ctype==2
    LuSpin=function(p)
      sort!(p)
      a=Int[]
      b=Int[]
      d=[0, 1, 0, -1][map(x->1+mod(x,4), p)]
      i=1
      while i<=length(p)
        l=p[i]
        t=Sum(d[1:i-1])
        if 1==mod(l, 4)
          push!(a, div(l-1, 4)-t)
          i+=1
        elseif 3==mod(l, 4)
          push!(b, div(l-3, 4)+t)
          i+=1
        else
          j=i
          while i<=length(p) && p[i]==l i+=1 end
          j=fill(0, max(0,div(i-j,2)))
          append!(a,j+div(l+mod(l,4),4)-t)
          append!(b,j+div(l-mod(l,4),4)+t)
        end
      end
      a=reverse(filter(!iszero,a))
      b=reverse(filter(!iszero,b))
      sum(d)>=1 ? [a, b] : [b, a]
    end
    addSpringer = function (f, i, s, k)
      ss = First(uc[:springerSeries], f)
      if s in [[[], [1]], [[], []]] p=1
      elseif s == [[1], []] p = 2
      else p = Position(CharParams(ss[:relgroup]), [s])
      end
      ss[:locsys][p] = [i, k]
    end
    trspringer = function (i, old, new)
        for ss in uc[:springerSeries]
            for c in ss[:locsys]
                if c[1] == i
                    p = Position(old, c[2])
                    if p != false c[2] = new[p] end
                end
            end
        end
    end
    d = 0
    while 4d^2-3d<=r
      i=4d^2-3d
      if mod(r-d,2)==0
          l = Concatenation(1:i, i + 2:(i + 4) - (i + 2):r)
          ss=Dict{Symbol, Any}(:relgroup=>coxgroup(:B,div(r-i,2)),
                               :levi => l, :Z => [-1])
          ss[:locsys]=fill([0,0],NrConjugacyClasses(ss[:relgroup]))
          push!(uc[:springerSeries],ss)
          i = 4 * d ^ 2 + 3d
          if i <= r && d != 0
            l = vcat(1:i,i+2:2:r)
            ss= Dict{Symbol, Any}(:relgroup=>coxgroup(:B,div(r-i,2)),
                                  :levi => l, :Z => [-1])
            ss[:locsys]=fill([0,0],NrConjugacyClasses(ss[:relgroup]))
            push!(uc[:springerSeries], ss)
          end
      end
      d+=1
    end
    l=filter(i->all(c->mod(c[1],2)==0 || c[2]==1,
             tally(uc[:classes][i][:parameter])),eachindex(uc[:classes])) 
    for i in l
      cl=uc[:classes][i]
      s=LuSpin(cl[:parameter])
      if length(cl[:Au]) == 1
          cl[:Au] = CoxeterGroup("A", 1)
          trspringer(i, [1], [2])
          d = 1
      elseif length(cl[:Au]) == 4
          cl[:Au] = CoxeterGroup("B", 2)
          trspringer(i, [1, 2, 3, 4], [1, 3, 5, 4])
          d = 2
      else
        error("Au non-commutative of order ",Size(cl[:Au])*2,"  !  implemented")
      end
      addSpringer(ss->ss[:Z]==[-1] && rank(ss[:relgroup])==sum(sum,s),i,s,d)
    end
  end
  return uc
end)

chevieset(:D,:UnipotentClasses,function(n,char)
  addSpringer = function (s, i, cc)
     ss = First(uc[:springerSeries], x->x[:defect] == DefectSymbol(s[:symbol]))
     if s[:sp] in [[[], [1]], [[], []]] p = 1
     elseif s[:sp] == [[1], []] p = 2
     else p = Position(CharParams(ss[:relgroup]), [s[:sp]])
     end
     ss[:locsys][p] = [i, Position(CharParams(cc[:Au]), 
                                   map(x->x ? [1,1] : [2], s[:Au]))]
  end
  function partition2DR(part)
    p=sort(reduce(vcat,map(x->1-x:2:x-1, part)))
    p=p[1+div(length(p),2):end]
    vcat([p[1]+p[2]], map(i->p[i+1]-p[i], 1:length(p)-1))
  end
  if char==2
    ss = XSP(4, 0, n, 1)
    symbol2partition = function (S)
      c=sort(reduce(vcat,S))
      i = 1
      part = Int[]
      ex = Int[]
      while i <= length(c)
          l = 2 * (c[i] - 2 * (i - 1))
          if i == length(c) || c[i + 1] - c[i] > 1
              push!(part, l + 2)
              i+=1
          elseif c[i + 1] - c[i] > 0
              part = Append(part, [l+1, l+1])
              i+=2
          else
              part = Append(part, [l, l])
              i+=2
              push!(ex, l)
          end
      end
      [reverse(filter(!iszero,sort(part))), ex]
    end
  else
    ss = XSP(2, 0, n, 1)
    function symbol2partition(S)
      c=sort(reduce(vcat,S))
      i = 1
      part = Int[]
      while i <= length(c)
        l = 2 * (c[i] - (i - 1))
        if i == length(c) || c[i + 1] - c[i] > 0
          push!(part, l+1)
          i+=1
        else
          part = append!(part, [l, l])
          i+=2
        end
      end
      reverse(filter(!iszero,sort(part)))
    end
  end
  l = union(map(c->map(x->
        [DefectSymbol(x[:symbol]), Sum(FullSymbol(x[:sp]), Sum)], c), ss))
  SortBy(l, x->[abs(x[1]), -sign(x[1])])
  uc = Dict{Symbol, Any}(:classes => [], :springerSeries => map(function(d)
      res = Dict{Symbol, Any}(:defect=>d[1], :levi=>1:n - d[2])
      if mod(n - d[2], 4) == 0 || char == 2
          res[:Z]= mod(n, 2)==0  ? [1, 1] : [1]
      else
          res[:Z]= mod(n, 2)==0  ? [-1, -1] : [-1]
      end
      if d[1]==0 res[:relgroup]=coxgroup(:D, d[2])
      else res[:relgroup] = coxgroup(:B, d[2])
      end
      res[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(res[:relgroup])]
      res
  end, l))
  for cl in ss
    cc = Dict{Symbol, Any}(:parameter => symbol2partition(cl[1][:symbol]))
    if char == 2
      cc[:dimBu] = (cl[1])[:dimBu]
      cc[:name] = join(map(function(x)
                            res=joindigits(fill(x[1], max(0, x[2])), "[]")
             if x[1] in cc[:parameter][2]
                 return SPrint("(", res, ")")
             end
             return res
         end, reverse(Collected(cc[:parameter][1]))), "")
    else
      cc[:dynkin] = partition2DR(cc[:parameter])
      cc[:name] = joindigits(cc[:parameter])
    end
    cc[:Au] = isempty(cl[1][:Au]) ? coxgroup() : 
       prod(coxgroup(:A,1) for i in eachindex(cl[1][:Au]))
    if char != 2
     cc[:red] = coxgroup()
     j = cc[:parameter]
     for j = Collected(j)
       if mod(j[1], 2) == 0
         cc[:red]*=coxgroup(:C, div(j[2],2))
       elseif mod(j[2], 2) != 0
         if j[2]>1 cc[:red]*=coxgroup(:B, div(j[2]-1,2)) end
       elseif j[2]>2 cc[:red]*=coxgroup(:D, div(j[2],2))
       else cc[:red]*=torus(1)
       end
     end
   end
   if !IsList(cl[1][:sp][2]) cl[1][:sp][3]=1-mod(div(n,2),2) end
   push!(uc[:classes], cc)
   for s in cl addSpringer(s, length(uc[:classes]), cc) end
   if !IsList(cl[1][:sp][2])
      cl[1][:sp][3]=1-cl[1][:sp][3]
      cc[:name]*="+"
      cc=deepcopy(cc)
      cc[:name]=replace(cc[:name],r".$"=>"-")
      if haskey(cc, :dynkin) cc[:dynkin][[1, 2]] = cc[:dynkin][[2, 1]] end
      push!(uc[:classes], cc)
      for s = cl addSpringer(s, length(uc[:classes]), cc) end
    end
  end
  if char == 2
    uc[:orderClasses] = Hasse(Poset(map((x->begin
    map(function (y,)
      m = maximum(((x[:parameter])[1])[1], ((y[:parameter])[1])[1])
      f = x-> map(i-> sum(filter(<(i),x))+i*count(>=(i),x) , 1:m)
      fx = f(x[:parameter][1])
      fy = f(y[:parameter][1])
      for i = 1:m
          if fx[i] < fy[i] return false
          elseif fx[i] == fy[i] && i in (y[:parameter])[2]
              if i in Difference(x[:parameter][1], x[:parameter][2])
                  return false
              end
              if i<m && mod(fx[i+1]-fy[i+1],2)==1 return false end
          end
      end
      if x[:parameter] == y[:parameter] && x != y return false end
      return true
    end, uc[:classes])
    end), uc[:classes])))
  else
    uc[:orderClasses] = Hasse(Poset(map((i->begin map((j->begin
      dominates(uc[:classes][j][:parameter], uc[:classes][i][:parameter]) && 
      (uc[:classes][j][:parameter]!=uc[:classes][i][:parameter] || i == j)
                                 end), 1:length(uc[:classes]))
                     end), 1:length(uc[:classes]))))
  end
  if char != 2
    d = 0
    while 4*d^2-d <= n
     i = 4*d^2-d
     if mod(n-d,2)==0
       l=Concatenation(1:i,i+2:2:n)
       s=Dict(:relgroup=>coxgroup(:B, div(n-i,2)),:levi=>l)
       if mod(n, 2) == 0 s[:Z]=[1, -1] else s[:Z] = [E(4)] end
       s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
       push!(uc[:springerSeries], s)
       if d==0 l=Concatenation([1], 4:2:n) end
       s=Dict(:relgroup => coxgroup(:B, div(n-i,2)),:levi=>l)
       if mod(n,2)==0 s[:Z] = [-1, 1]
       else s[:Z] = [-(E(4))]
       end
       s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
       push!(uc[:springerSeries], s)
       i=4d^2+d
       if d != 0 && i <= n
         l = vcat(1:i, i+2:2:n)
         s = Dict(:relgroup=>coxgroup(:B, div(n - i,2)),:levi=>l)
         if mod(n,2)==0 s[:Z]=[1, -1] else s[:Z]=[E(4)] end
         s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
         push!(uc[:springerSeries], s)
         s=Dict(:relgroup=>oxgroup("B",div(n-i,2)),:levi=>l)
         if mod(n, 2) == 0 s[:Z] = [1, 1] else s[:Z] = [-(E(4))] end
         s[:locsys]=[[0,0] for i in 1:NrConjugacyClasses(s[:relgroup])]
         push!(uc[:springerSeries], s)
       end
     end
     d+=1
  end
  function LuSpin(p)
    sort!(p)
    a=Int[]
    b=Int[]
    d=[0,1,0,-1][map(x->1 + mod(x, 4), p)]
    i=1
    while i<=length(p)
      l=p[i]
      t=sum(d[1:i-1])
      if 1==mod(l,4)
        push!(a, div(l-1,4)-t)
        i+=1
      elseif 3==mod(l,4)
        push!(b, div(l-3,4)+t)
        i+=1
      else
        j=i
        while i<=length(p) && p[i]==l i+=1 end
        j=fill(0, max(0,div(i-j,2)))
        append!(a,j+div(l+mod(l,4),4)-t)
        append!(b,j+div(l-mod(l,4),4)+t)
      end
    end
    a=reverse(filter(!iszero,a))
    b=reverse(filter(!iszero,b))
    sum(d)>=1 ? [a, b] : [b, a]
  end
  addSpringer1= function (f, i, s, k)
    ss = First(uc[:springerSeries], f)
    if s in [[[], [1]], [[], []]] p = 1
    elseif s == [[1], []] p = 2
    else p = Position(CharParams(ss[:relgroup]), [s])
    end
    ss[:locsys][p] = [i, k]
  end
  trspringer = function (i, new)
    for ss in uc[:springerSeries]
      for c in ss[:locsys] if c[1] == i c[2] = new[c[2]] end end
    end
  end
  l=Filtered(1:length(uc[:classes]), i->all(c->mod(c[1], 2) == 0 || c[2] == 1,
                 Collected(((uc[:classes])[i])[:parameter])))
  for i in l
     cl=uc[:classes][i]
     s=LuSpin(cl[:parameter])
     if length(cl[:Au]) == 1
       cl[:Au] = coxgroup(:A, 1)
       trspringer(i, [2])
       k = [1, 1]
     elseif length(cl[:Au]) == 2
       cl[:Au] = coxgroup(:A, 1)*coxgroup(:A,1)
       trspringer(i, [2, 4])
       k = [1, 3]
     elseif length(cl[:Au]) == 8
       cl[:Au] = coxgroup(:A,1)*coxgroup(:B,2)
       trspringer(i, [1, 6, 8, 5, 10, 3, 4, 9])
       k = [2, 7]
     else
       error("Au non-commutative of order ",length(cl[:Au])*2,"  ! implemented")
     end
     if !('-' in cl[:name])
         addSpringer1(ss->ss[:Z] in [[1, -1], [E(4)]] && 
                     rank(ss[:relgroup]) == Sum(s, Sum) , i, s, k[1])
     end
     if !('+' in cl[:name])
         addSpringer1(ss->ss[:Z] in [[-1, 1], [-(E(4))]] && 
                     rank(ss[:relgroup]) == Sum(s, Sum), i, s, k[2])
     end
  end
end
  uc
end)
chevieset(:I, :CharInfo, function(m)
  local res, applyf, v, m1
  res = Dict{Symbol, Any}(:charparams => [Any[1, 0]])
  if mod(m, 2) == 0
      res[:extRefl] = [1, 5, 4]
      m1 = div(m, 2)
      push!(res[:charparams], [1, m1, "'"], [1, m1, "''"])
  else
      res[:extRefl] = [1, 3, 2]
  end
  push!(res[:charparams], [1, m])
  append!(res[:charparams], map(i->[2, i], 1:div(m - 1, 2)))
  res[:b]=map(x->x[2], res[:charparams])
  res[:B]=map(phi->phi[1] == 1 ? phi[2] : m - phi[2], res[:charparams])
  res[:a]=map(phi->phi[1]!=1 || phi[2]==div(m,2) ? 1 : phi[2], res[:charparams])
  res[:A]=map(phi->phi[1]==1 || phi[2]==div(m,2) ? m-1 : phi[2], res[:charparams])
  res[:charSymbols] = map(function (l)
              S = map(i->[0], 1:m)
              k = 0
              if k != 0
                  S[1] = [0, 1]
                  S[1 + mod(k + l, m)] = [0, 1]
                  S[k + 1] = []
                  S[l + 1] = []
              else
                  S[1] = [1]
                  S[l + 1] = [1]
              end
              S
          end, 1:div(m - 1, 2))
  v=map(x->[0],1:m);v[m]=[1,2];pushfirst!(res[:charSymbols],v)
  if mod(m,2)==0
    v=map(x->[0],1:m);v[m]=[1];v[m1]=[1];pushfirst!(res[:charSymbols],v)
    v=map(x->[0], 1:m); v[m]=[1];v[m1]=[1]; pushfirst!(res[:charSymbols],v)
  end
  v=map(x->[0,1],1:m)
  v[m]=[2]
  pushfirst!(res[:charSymbols],v)
  res[:malleParams] = map(x->Vector{Any}(map(partβ,x)),res[:charSymbols])
  if mod(m,2)==0
    res[:malleParams][2]=push!(res[:malleParams][2][1:m1],1)
    res[:malleParams][3]=push!(res[:malleParams][3][1:m1],-1)
  end
  return res
end)

struct Hook
  area::Int       # area of the hook
  startpos::Int   # the position of the decreased βnumber in the βlist
  endpos::Int     # the position it will occupy after being decreased
  DC::Vector{Int} # the list of dull corners, given by their axial position
  SC::Vector{Int} # the list of sharp corners, given by their axial position
end

function Base.show(io::IO,h::Hook)
  print(io,"<",h.startpos,"=>",h.endpos,",area:",h.area,",DC=",h.DC,",SC=",h.SC,">")
end

"""
hooksβ(S,s)  the list of all hooks  of area less than or  equal to s in the
Young diagram of the partition with βlist S.
"""
function hooksβ(S,s)
  res=Hook[]
  e=length(S)
  if e==0 return res end
  j=e
  for i in S[e]-1:-1:0
    while j>0 && S[j]>i j-=1 end
    if j>0 && S[j]==i continue end
    for k in j+1:e
      if S[k]>i+s break end
      z=vcat([i],S[j+1:k-1])
      zi=filter(i->z[i]-z[i-1]>1,2:length(z))
      push!(res,Hook(S[k]-i,k,j+1,z[zi].-e,z[push!(zi.-1,length(z))].+(1-e)))
    end
  end
  res
end

"""
stripsβ(S,s)  returns as a list of lists  of Hooks all broken border strips
(union  of disjoint  hooks) of  area less  than or  equal to s in the Young
diagram of the partition with βlist S.
"""
function stripsβ(S,s)
  res=[Hook[]]
  for hook in hooksβ(S,s)
    if s==hook.area push!(res,[hook])
    else
      j=hook.endpos-1-length(S)
      for hs in stripsβ(S[1:hook.endpos-1], s-hook.area)
        for h in hs
          h.SC.+=j
          h.DC.+=j
        end
        push!(hs, hook)
        push!(res, hs)
      end
    end
  end
  res
end

"the elementary symmetric function of degree t in the variables v"
function elementary_symmetric_function(t,v)
  if t==0 return 1 end
  sum(x->isempty(x) ? 1 : prod(v[x]),combinations(1:length(v),t))
end

"the homogeneous symmetric function of degree t in the variables v"
function homogeneous_symmetric_function(t,v)
  if t==0 return 1 end
  sum(x->prod(v[x]),combinations(repeat(1:length(v),t), t))
end

chevieset(:imp, :HeckeCharTable, function (p, q, r, para, rootpara)
  res=Dict{Symbol, Any}(:name=>SPrint("H(G(", p, ",", q, ",", r, "))"),
  :degrees=> chevieget(:imp, :ReflectionDegrees)(p, q, r), :dim=>r)
  res[:identifier] = res[:name]
  res[:size] = prod(res[:degrees])
  res[:order] = res[:size]
  res[:powermap]=chevieget(:imp,:PowerMaps)(p,q,r)
  cl=chevieget(:imp, :ClassInfo)(p, q, r)
  if r==1
    res[:classes]=cl[:classes]
    res[:orders]=cl[:orders]
    res[:irreducibles]=map(i->map(j->para[1][i]^j,0:p-1),1:p)
  elseif q==1
# character table of the Hecke algebra of G(p,1,r) parameters v and [q,-1]
# according to "Characters of Iwahori-Hecke algebras of G(r,p,n)"
# Halverson and Ram, Canadian Journal of math. 50 (1998) 167--192
    merge!(res,cl)
    T=NamedTuple{(:area,:cc,:hooklength,:DC,:SC,:stripped),
       Tuple{Int,Int,Int,Vector{Pair{Int,Int}},Vector{Pair{Int,Int}},
       Vector{Vector{Int}}}}
    StripsCache=Dict{Pair{Vector{Vector{Int}},Int},Any}()
"""
Strips(S,s)  returns the list  of collections of  broken border strips of
total  area equal to  s coming from  the various βlists  in the symbol S.
Each  collection is represented by named tuple containing the statistical
information   about  it  necessary  to  compute  the  function  Delta  in
Ram-Halverson 2.17. These tuples have the following fields:
   area       total area
   cc         number of connected components (hooks)
   hooklength Sum of hooklengths of individual hooks (the # of rows -1 of the
                  hook, equal to startpos-endpos)
   DC         the list of all dull corners, represented by a pair:
                which βlist they come from, axial position
   SC         the list of all sharp corners, represented by a pair:
                which βlist they come from, axial position
   stripped  the symbol left after removing the strip collection
"""
    function Strips(S,s)local e, res, hs, ss, a
      get!(StripsCache,S=>s) do
      function strip(S,hs) # strip hooks hs from βlist S
        S=copy(S)
        for h in hs 
          beta=S[h.startpos]
          for i in h.startpos:-1:h.endpos+1 S[i]=S[i-1] end
          S[h.endpos]=beta-h.area
        end
        i=1
        while i<=length(S) && iszero(S[i]) i+=1 end
        if i>1 S=S[i:end].-(i-1) end
        return S
      end 
      e=length(S)
      if e==0
        if s==0 return [(area=0,cc=0,hooklength=0,DC=Pair{Int,Int}[],
                         SC=Pair{Int,Int}[],stripped=Vector{Int}[])]
        else return T[]
        end
      end
      res=T[]
      for h in stripsβ(S[e], s)
        hs=(area=reduce(+,map(x->x.area,h),init=0), 
          cc=length(h),
          hooklength=reduce(+,map(x->x.startpos-x.endpos,h),init=0),
          DC=[e=>y for x in h for y in x.DC],
          SC=[e=>y for x in h for y in x.SC],
          stripped=[strip(S[e],h)])
        for a in Strips(S[1:e-1],s-hs.area)
          push!(res,(area=a.area+hs.area,cc=a.cc+hs.cc,
            hooklength=a.hooklength+hs.hooklength,
            DC=vcat(a.DC,hs.DC), SC=vcat(a.SC,hs.SC),
            stripped=vcat(a.stripped,hs.stripped)))
        end
      end
      res
      end
    end
# the function Delta of Ram-Halverson 2.17, modified to take in account that
# our eigenvalues for T_2..T_r are Q[1] and Q[2] instead of q and -q^-1
    function Delta(k,hs,(Q1,Q2),v)local q
      delta=1
      if hs.cc>1
        if k==1 || iszero(Q1+Q2) return 0
        else delta*=(Q1+Q2)^(hs.cc-1)
        end
      end
      q=-Q1//Q2
      delta*=Q1^(hs.area-hs.cc)*(-q)^-hs.hooklength
      if k==0 return delta end
      ctSC=[v[x]*q^y for (x,y) in hs.SC]
      ctDC=[v[x]*q^y for (x,y) in hs.DC]
      delta*=reduce(*,ctSC;init=1)//reduce(*,ctDC;init=1)
      if k==1 return delta end
      delta*(-1)^(hs.cc-1)*sum(map(
          t->(-1)^t*elementary_symmetric_function(t,ctDC)* 
                    homogeneous_symmetric_function(k-t-hs.cc,ctSC), 
               0:min(length(ctDC),k-hs.cc)))
    end
    chiCache=Dict{Pair{Vector{Vector{Int}},Vector{Vector{Int}}}, Any}()
    function entry(lambda,mu)
      get!(chiCache,lambda=>mu) do
      n=sum(sum,lambda)
      if iszero(n) return 1 end
      bp=maximum(x for S in lambda for x in S)
      i=findfirst(x->bp in x,lambda)
      # choice of bp and i corresponds to choice (Sort) in classtext
      strips=Strips(mu, bp)
      if isempty(strips) return 0 end
      rest=copy(lambda)
      rest[i]=rest[i][2:end]
      (-prod(para[2]))^((i-1)*(n-bp))*sum(strips) do x
        d=Delta(i-1,x,para[2],para[1])
        iszero(d) ? d : d*entry(rest, x.stripped)
      end
      end
    end
    res[:irreducibles]=map(x->map(y->entry(y,x),cl[:classparams]),
                           map(x->map(βset,x),cl[:classparams]))
  elseif [q,r]==[2,2] && !haskey(CHEVIE,:othermethod)
    res[:classes]=cl[:classes]
    res[:orders]=cl[:orders]
    Z,X,Y=para[1:3]
    ci=chevieget(:imp, :CharInfo)(p, q, r)
    function entry2(char, class)
      char=ci[:malle][findfirst(==(char),ci[:charparams])]
      if char[1] == 1
        w=[Z[char[4]],X[char[2]],Y[char[3]]]
        return Product(class, i->iszero(i) ? prod(w) : w[i])
      else
        w=char[2]*root(X[1]*X[2]*Y[1]*Y[2]*Z[char[3]]*Z[char[4]]*
          E(div(p,q),2-char[3]-char[4]))*E(p,char[3]+char[4]-2)
        class=map(i->count(==(i),class), 0:3)
        if class[2]>0 char=sum(x->x^class[2],Z[char[[3,4]]])
        elseif class[3]>0 char=sum(X)
        elseif class[4]>0 char=sum(Y)
        else char = 2
        end
        return w^class[1]*char
      end
    end
    res[:irreducibles] = map(char->
      map(class->entry2(char, class), cl[:classparams]), ci[:charparams])
  else
    res[:classnames]=cl[:classnames]
    res[:orders]=cl[:orders]
    res[:centralizers]=cl[:centralizers]
    res[:classes] = map(x->div(res[:size],x),res[:centralizers])
    res[:irreducibles] = map(i->traces_words_mats(
             toM.(chevieget(:imp,:HeckeRepresentation)(p,q,r,para,[],i)),
             cl[:classtext]),1:length(res[:classes]))
  end
  res[:centralizers]=map(x->div(res[:size],x), res[:classes])
  res[:parameter]=para
  res[:irreducibles]*=prod(prod,para)^0
  return CHEVIE[:compat][:MakeCharacterTable](res)
end)

chevieset(:imp, :CharTable, function (p, q, r)
  oo=denominator.(chevieget(:imp,:EigenvaluesGeneratingReflections)(p,q,r))
  return chevieget(:imp, :HeckeCharTable)(p, q, r, 
      map(o->o==2 ? [1,-1] : map(i->E(o,i),0:o-1),oo), [])
end)

ExpandRep = function (r, d, l)
  T=reduce(promote_type,typeof.(first.(l)))
  m=map(i->map(j->fill(zero(T),d),1:d), 1:r)
  for v in l
    for k in Int.(v[2:end])
      i=k%(d^2)
      m[div(k,d^2)+1][div(i,d)+1][mod(i,d)+1]=v[1]
    end
  end
  return m
end

const G4_22IndexChars_dict=Dict{Int,Any}()
for i in 4:22 G4_22IndexChars_dict[i]=Dict() end
CHEVIE[:CheckIndexChars]=true

function G4_22FetchIndexChars(ST, para)
  if !CHEVIE[:CheckIndexChars]
    return chevieget(:G4_22, :CharInfo)(ST)[:indexchars]
  end
  get!(G4_22IndexChars_dict[ST],para)do
    chevieget(:G4_22, :HeckeCharTable)(ST, para, [])[:indexchars]
  end
end

# test if res=Chartable(Hecke(G_ST,res[:parameter])) is correct
# where rows is correspond irrs for generic group and i the selector from rows
function G4_22Test(res,rows,i)
  ST=res[:ST]
  T(ST)=string("G",ST in 4:7 ? 7 : ST in 8:15 ? 11 : 19)
  if haskey(G4_22IndexChars_dict[ST],res[:parameter])
    InfoChevie("Using G4_22IndexChars_dict[$ST][",res[:parameter],"]\n")
    ic=G4_22IndexChars_dict[ST][res[:parameter]]
    res[:irreducibles]=rows[ic]
    if ic!=i
      println("*** WARNING: choice of character restrictions from ", T(ST), 
        " for this specialization does not agree with group CharTable")
      if !CHEVIE[:CheckIndexChars]
        print("Try again with CHEVIE[:CheckIndexChars]=true\n")
      end
    end
    return ic
  end
  ic=i
  res[:irreducibles]=rows[ic]
  if length(Set(res[:irreducibles]))==length(res[:classes]) l=i
  else
    l=map(x->findfirst(==(x),rows), rows)
    if length(Set(l))!=length(res[:classes])
      error("specialization not semi-simple")
    end
    l=map(x->filter(i->l[i]==x,eachindex(l)), sort(unique(l)))
    println("*** WARNING: bad choice of character restrictions from ", T(ST), 
            " for this specialization\n")
    if !CHEVIE[:CheckIndexChars]
      print("Try again with CHEVIE[:CheckIndexChars]=true\n")
    end
    o=filter(x->count(j->j in x,i)>1,l)
    println(" over-represented by ", intersect(union(o...), i)," : ", o)
    println(" absent : ",filter(x->iszero(count(j->j in x,i)),l))
    l= map(x->x[1],l)
    print(" Choosing ",l,"\n")
    res[:irreducibles]=rows[l]
  end
  G4_22IndexChars_dict[ST][res[:parameter]]=l
end

chevieset(:imp, :ReflectionCoDegrees, function (p, q, r)
  res=collect(p*(0:r-1))
  if p==q && p>=2 && r>2 res[r]-=r end
  res
end)
