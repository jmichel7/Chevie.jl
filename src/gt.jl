function RationalUnipotentClasses(WF, p)
  u=UnipotentClasses(WF, p)
  t=Ucl.GreenTable(u;classes=true)
  map(i->Dict{Symbol, Any}(:card => CycPol(t[:cardClass][i]), 
                           :class => u.classes[t[:locsys][i][1]], 
                           :classno => t[:locsys][i][1], 
                           :AuNo => t[:locsys][i][2]), 
           eachindex(t[:locsys]))
end

# returns the Poset of closed subsystems of the root system of W
ClosedSubsets = function(W)
  gets(W, :closedsubsets)do
  function possum(i,j)
    p=findfirst(==(roots(W,i)+roots(W,j)),roots(W))
    isnothing(p) ? 0 : p
  end
  psum=[possum(i,j) for i in 1:2nref(W),  j in 1:2nref(W)]
  closure=function(l,new)
    nnew = new
    while true
      new = nnew
      nnew = Int[]
      for i in new, j in l
        if psum[i,j]!=0 push!(nnew, psum[i,j]) end
      end
      l = union(l, new)
      nnew = setdiff(nnew, l)
      if isempty(nnew) break end
    end
    return sort(l)
  end
  l = [Int[]]
  new = [1]
  covered = [Int[]]
  for w in new
    for f in setdiff(1:nref(W), l[w])
      n = closure(l[w], [f, f + nref(W)])
      p = findfirst(==(n),l)
      if isnothing(p)
          push!(l, n)
          push!(covered, [w])
          push!(new, length(l))
      else push!(covered[p], w)
      end
    end
  end
  covered=unique.(covered)
  P=Poset(incidence(Poset(covered)))
  P.prop[:elements]=l
  P.prop[:label]=function(io,i) join(l[i]," ") end
  P
  end
end

struct ClassType
  CGs
  cent::CycPol
  unip
end

struct ClassTypes
  p::Int
  WF
  ss::Vector{ClassType}
  prop::Dict{Symbol,Any}
end

# ClassTypes(W[,p])
function ClassTypes(W,p=0)
  if W isa Spets WF=W;W=Group(W)
  else WF=spets(W)
  end
  l=vcat(twistings.(Ref(WF),SScentralizer_representatives(Group(WF), p))...)
  ClassTypes(p,WF,map(x->ClassType(x,CycPol(generic_order(x,Pol(:q))),
    RationalUnipotentClasses(x,p)),l),Dict{Symbol,Any}())
end

function Base.show(io::IO,r::ClassTypes)
  res=string("ClassTypes(",sprint(show,r.WF;context=io))
  if r.p==0 res*=",good characteristic)"
  else res*=string(",char. ",r.p,")")
  end
  if haskey(r.prop,:specialized)
    res*=string(" ",join(map(x->string(x[1],"==",x[2]),r.prop[:specialized])," "))
  end
  println(io,res)
  function nc(p)
    local d
    p=Mvp(p)
    d=lcm(map(denominator,p[:coeff]))
    if d==1 return Format(p, opts) end
    return SPrint(BracketIfNeeded(Format(p*d,opts),"+-"),"/",d)
  end
  classes=get(io,:nClasses,false)
  if classes NrConjugacyClasses(r) end
  if get(io,:unip,false)
    rowLabels=[]
    columnLabels=[]
    if classes push!(columnLabels, "nrClasses") end
    push!(columnLabels,"u")
    push!(columnLabels,"Centralizer")
    t=[]
    for x in r.ss
      u=RationalUnipotentClasses(x.CGs, r.p)
      for c in u
        v=String[]
        if isone(c[:card])
          if classes push!(v, nc(x[:nrClasses])) end
          push!(rowLabels,sprint(show,x.CGs;context=io))
        else
          push!(rowLabels," ")
          if classes push!(v,"") end
        end
        push!(v,Ucl.nameclass(merge(c[:class].prop,Dict(:name=>c[:class].name)),
                              merge(io.dict,Dict(:class=>c[:AuNo]))))
        push!(v,sprint(show,x.cent//c[:card],context=io))
        push!(t,v)
      end
    end
    t=toM(t)
  else
    rowLabels=map(x->sprint(show,x.CGs;context=io),r.ss)
    t=[]
    columnLabels=String[]
    if classes
      push!(t,nc.(NrConjugacyClasses(r)))
      push!(columnLabels, "nrClasses")
    end
    push!(t,map(x->sprint(show,x.cent;context=io),r.ss))
    push!(columnLabels,"Centralizer")
    t=permutedims(toM(t))
  end
  format(io,t;col_labels=columnLabels,row_labels=rowLabels,rows_label="Type")
end

function Value(r::ClassTypes,arg...)
  r=copy(r)
  r[:specialized] = map(i->arg[2][i-1:i], 2:2:length(arg[2]))
  r[:ss] = map(function (x,)
              local res
              res = copy(x)
              NrConjugacyClasses(r)
              res[:nrClasses] = Value(x[:nrClasses], arg[2])
              return res
          end, r[:ss])
  return r
end

# See Fleischmann-Janiszczak AAECC 1996 definition 2.1
function NrConjugacyClasses(C::ClassTypes)
  W=Group(C.WF)
  b=keys(factors(BadNumber(W)))
  if length(fundamental_group(W))>1
    print("# Nr classes each type implemented only for simply connected groups")
    return
  end
  for r in C.ss
    if !haskey(r,:nrClasses)
      HF=r.CGs
      H=Group(HF)
      P=deepcopy(ClosedSubsets(W))
      o=Filtered(1:Size(P),i->OnSets(P[:elements][i],HF[:phi])==P[:elements][i])
      o=Filtered(o,i-> all(j->j in P[:elements][i],inclusiongens(H)))
      P=Restricted(P, o)
      P[:elements] = P[:elements][o]
# here P poset of HF.phi-stable closed subsets containing roots(H)
      InfoChevie("# ", HF, "==>", P, "")
      l = map(x->Spets(ReflectionSubgroup(W, x), HF[:phi]), P[:elements])
      l = map(l)do RF
        local res, d, R
        R = Group(RF)
        res=prod(p->Mvp(:q)-p[2],filter(y->y[1]==1,ReflectionDegrees(RF)))
        if R[:semisimpleRank] == 0 return res end
        d=SmithNormalFormMat(R[:roots][R[:generatingReflections]]*R[:simpleRoots])
        d=Filtered(Flat(d), x->!x in [0,1,C[:p]])
        res*Product(d, function (i,)
                        if i == 2 && (2 in b && C[:p] != 2) return 2 end
                        if mod(C[:p], i) == 1 return i end
                        return Mvp(SPrint("q_", i))
                    end)
      end
      less= i->Difference(ListBlist(1:length(l), Incidence(P)[i]), [i])
      o = LinearExtension(P)
      mu = []
      mu[o[Size(P)]] = 1
      for i in o[length(P)-1:-1:1] mu[i]=-Sum(mu[less(i)]) end
      n=Stabilizer(W,gapSet(H[:rootInclusion][H[:generatingReflections]]),OnSets)
      n = (mu * l) // Size(Centralizer(n, HF[:phi]))
      InfoChevie("==>", n, ":", Stime(), "\n")
      r[:nrClasses] = n
    end
  end
  C.ss=Filtered(C.ss, x->x[:nrClasses]!=0)
  map(x->x[:nrClasses],C.ss)
end
