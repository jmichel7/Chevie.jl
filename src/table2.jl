# replacements for some functions in tables (whose automatic translation failed)
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

chevieset(:A,:HeckeCharTable,(n,para,root)->chevieget(:imp,:HeckeCharTable)(1,1,n+1,para,root))

chevieset(Symbol("2A"),:HeckeCharTable, function (r, param, rootparam)
  q = -param[1][1] // param[1][2]
  if isnothing(rootparam[1]) v=GetRoot(q, 2, "CharTable(Hecke(2A))")
  else v = rootparam[1]
  end
  W = coxgroup(:A, r)
  qE = HeckeCentralMonomials(hecke(W, v))
  H = hecke(W, v^-2)
  T = Tbasis(H)
  tbl = copy(CharTable(H))
  Inherit(tbl, (chevieget(Symbol("2A"), :ClassInfo))(r))
  tbl[:identifier] = "H(^2A_$r)"
  cl = map(x->T(W(x...)*longest(W)), tbl[:classtext])
  tbl[:irreducibles] = permutedims(map(HeckeCharValues, cl))
  A=chevieget(:A,:LowestPowerFakeDegree).(chevieget(:A,:CharInfo)(r)[:charparams])
  tbl[:irreducibles]= (-1).^A .* qE .* tbl[:irreducibles]
  CHEVIE[:compat][:AdjustHeckeCharTable](tbl, param)
  return tbl
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

chevieset(:B,:UnipotentClasses,function(r,char,ctype)
  part2dynkin=function(part)
    p=sort(reduce(vcat,map(d->1-d:2:d-1, part)))
    p=p[div(3+length(p),2):end]
    res= ctype==1 ? [2*p[1]] : [p[1]]
    append!(res,p[2:end]-p[1:end-1])
  end
  addSpringer1=function(s,cc)
    ss=First(uc[:springerSeries],x->x[:defect]==DefectSymbol(s[:symbol]))
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
  l = union(map(c->map(x->[DefectSymbol(x[:symbol]), Sum(x[:sp], Sum)], c), ss))
  sort!(l,by=x->[AbsInt(x[1]),-SignInt(x[1])])
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
      cc[:name] = Join(map(function (x,)
          res = joindigits(fill(0, max(0, (1 + x[2]) - 1)) + x[1], "[]")
          if x[1] in cc[:parameter][2] return string("(", res, ")") end
          return res
      end, reverse(Collected(cc[:parameter][1]))), "")
    end
    cc[:red] = coxgroup()
    if char == 2 j = cc[:parameter][1]
    else j = cc[:parameter]
    end
    for j in Collected(j)
      if mod(j[1], 2) == mod(ctype, 2)
        cc[:red] = cc[:red] * coxgroup(:C, div(j[2],2))
      elseif mod(j[2], 2) != 0
        if j[2]>1 cc[:red]*=coxgroup(:B, div(j[2]-1,2)) end
       elseif j[2]>2 cc[:red]*=coxgroup(:D, div(j[2],2))
      else cc[:red]*=Torus(1)
      end
    end
    push!(uc[:classes], cc)
    for s in cl addSpringer1(s,cc) end
  end
  uc[:orderClasses] = Hasse(Poset(map(x->
      map(function(y)
        if char != 2 return Dominates(y[:parameter], x[:parameter]) end
        m = maximum(((x[:parameter])[1])[1], ((y[:parameter])[1])[1])
        f = x-> map(i->sum(filter(<(i),x))+i*count(>=(i),x) ,1:m)
        fx = f(x[:parameter][1])
        fy = f(y[:parameter][1])
        for i in 1:m
          if fx[i] < fy[i] return false
          elseif fx[i] == fy[i] && i in (y[:parameter])[2]
            if i in Difference(x[:parameter][1], x[:parameter][2]) return false end
            if i < m && mod(fx[i+1]-fy[i+1],2)==1 return false end
          end
        end
        return true
      end, uc[:classes]), uc[:classes])))
  if char != 2 && ctype == 2
    LuSpin=function(p)
      sort!(p)
      a = []
      b = []
      d = [0, 1, 0, -1][map(x->1+mod(x,4), p)]
      i = 1
      while i <= length(p)
          l = p[i]
          t = Sum(d[1:i - 1])
          if 1 == mod(l, 4)
              push!(a, div(l - 1, 4) - t)
              i = i + 1
          elseif 3 == mod(l, 4)
              push!(b, div(l - 3, 4) + t)
              i = i + 1
          else
              j = i
              while i <= length(p) && p[i] == l i = i + 1 end
              j = fill(0, max(0, (1 + div(i - j, 2)) - 1))
              a = Append(a, (j + div(l + mod(l, 4), 4)) - t)
              b = Append(b, j + div(l - mod(l, 4), 4) + t)
          end
      end
      a = Vector{Int}(reverse(filter(x->x!=0,a)))
      b = Vector{Int}(reverse(filter(x->x!=0,b)))
      if Sum(d) >= 1 return [a, b]
      else return [b, a]
      end
    end
    addSpringer = function (f, i, s, k)
      ss = First(uc[:springerSeries], f)
      if s in [[[], [1]], [[], []]] p = 1
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
    l = Filtered(eachindex(uc[:classes]), i->
     ForAll(Collected(uc[:classes][i][:parameter]), c->
                                mod(c[1],2)==0 || c[2]==1) )
    for i = l
      cl = (uc[:classes])[i]
      s = LuSpin(cl[:parameter])
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
  SortBy(l, x->[AbsInt(x[1]), -SignInt(x[1])])
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
      cc[:name] = Join(map(function(x)
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
      cc = Copy(cc)
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
      Dominates(uc[:classes][j][:parameter], uc[:classes][i][:parameter]) && 
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
    a = Int[]
    b = Int[]
    d = [0, 1, 0, -1][map(x->1 + mod(x, 4), p)]
    i = 1
    while i <= length(p)
        l = p[i]
        t = Sum(d[1:i - 1])
        if 1 == mod(l, 4)
          push!(a, div(l - 1,4) - t)
          i+=1
        elseif 3 == mod(l, 4)
          push!(b, div(l - 3,4) + t)
          i+=1
        else
            j = i
            while i <= length(p) && p[i]==l i+=1 end
            j = fill(0, max(0, div(i - j,2)) )
            a = Append(a, j + div(l + mod(l, 4), 4) - t)
            b = Append(b, j + div(l - mod(l, 4), 4) + t)
        end
    end
    a = Filtered(a, (x->begin x != 0 end))
    a = reverse(a)
    b = Filtered(b, (x->begin x != 0 end))
    b = reverse(b)
    if Sum(d) >= 1 return [a, b]
    else return [b, a]
    end
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
  l=Filtered(1:length(uc[:classes]), (i->begin
                 ForAll(Collected(((uc[:classes])[i])[:parameter]), c->
                             mod(c[1], 2) == 0 || c[2] == 1) end))
  for i in l
     cl = (uc[:classes])[i]
     s = LuSpin(cl[:parameter])
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
  res[:malleParams] = map(x->Vector{Any}(map(PartBeta,x)),res[:charSymbols])
  if mod(m,2)==0
    res[:malleParams][2]=push!(res[:malleParams][2][1:m1],1)
    res[:malleParams][3]=push!(res[:malleParams][3][1:m1],-1)
  end
  return res
end)

