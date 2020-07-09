function RationalUnipotentClasses(WF, p)
  u=UnipotentClasses(WF, p)
  t=Ucl.GreenTable(u;classes=true)
  map(i->Dict{Symbol, Any}(:card => CycPol(t[:cardClass][i]), 
                           :class => u.classes[t[:locsys][i][1]], 
                           :classno => t[:locsys][i][1], 
                           :AuNo => t[:locsys][i][2]), 
           eachindex(t[:locsys]))
end

ClassTypesOps=Dict()
ClassTypesOps[:Format] = function (r, opts)
  res = string(r)
  res = Append(res, "\n")
  nc = function (p,)
          local d
          p = Mvp(p)
          d = Lcm(map(denominator, p[:coeff]))
          if d == 1 return Format(p, opts) end
          return SPrint(BracketIfNeeded(Format(p * d, opts), "+-"), "/", d)
      end
  if haskey(opts, :nrClasses) NrConjugacyClasses(r) end
  if haskey(opts, :unip)
      o = Dict{Symbol, Any}(:rowLabels => [], :rowsLabel => "Type", :columnLabels => [])
      if haskey(opts, :nrClasses) push!(o[:columnLabels], "nrClasses") end
      push!(o[:columnLabels], "u")
      push!(o[:columnLabels], "Centralizer")
      t = []
      for x = r[:ss]
          u = RationalUnipotentClasses(x[:CGs], r[:p])
          for c = u
              v = []
              if c[:card] == CycPol(1)
                  if haskey(opts, :nrClasses) push!(v, nc(x[:nrClasses])) end
                  push!(o[:rowLabels], ReflectionName(x[:CGs], opts))
              else
                  push!(o[:rowLabels], " ")
                  if haskey(opts, :nrClasses) push!(v, "") end
              end
              push!(v, Name(c[:class], Inherit(Dict{Symbol, Any}(:class => c[:AuNo]), opts)))
              push!(v, Format(x[:cent] // c[:card], Inherit(Dict{Symbol, Any}(:vname => "q"), opts)))
              push!(t, v)
          end
      end
  else
      o = Dict{Symbol, Any}(:rowLabels => map((x->begin
                              ReflectionName(x[:CGs], opts) end), r[:ss]))
      o[:rowsLabel] = "Type"
      t = []
      o[:columnLabels] = []
      if haskey(opts, :nrClasses)
          push!(t, map(nc, NrConjugacyClasses(r)))
          push!(o[:columnLabels], "nrClasses")
      end
      push!(t, map((x->begin
                      Format(x[:cent], Inherit(Dict{Symbol, Any}(:vname => "q"), opts))
                  end), r[:ss]))
      push!(o[:columnLabels], "Centralizer")
      t = TransposedMat(t)
  end
  res = Append(res, FormatTable(t, Inherit(o, opts)))
  return res
end
ClassTypesOps[:Value] = function (arg...,)
  r=copy(arg[1])
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
ClassTypesOps[:Display] = function (r, opts) print(Format(r, opts)) end
ClassTypesOps[:String] = function (r,)
  res = SPrint("ClassTypes(", r[:spets])
  if r[:p] == 0 res = Append(res, ",good characteristic)")
  else res *= SPrint(",char. ", r[:p], ")")
  end
  if haskey(r, :specialized)
    res*=SPrint(" ",Join(map(x->SPrint(x[1],"==",x[2]),r[:specialized])," "))
  end
  return res
end
ClassTypesOps[:Print] = function (r,)
        print(string(r))
    end

# ClassTypes(W[,p])
ClassTypes = function (W,p=0)
  if IsSpets(W) WF = W
      W = Group(W)
  else WF = Spets(W)
  end
  l = SemisimpleCentralizerRepresentatives(Group(WF), p)
  l = Concatenation(map(x->Twistings(WF, x), l))
  Dict{Symbol, Any}(:p => p, :spets => WF, :ss => map((x->begin
                          Dict{Symbol, Any}(:CGs => x, :cent => CycPol(GenericOrder(x, Mvp("q"))), :unip => RationalUnipotentClasses(x, p))
                      end), l), :operations => ClassTypesOps)
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

# See Fleischmann-Janiszczak AAECC 1996 definition 2.1
ClassTypesOps[:NrConjugacyClasses] = function (C,)
        local HF, W, H, o, P, l, less, mu, n, i, r, b
        W = Group(C[:spets])
        b = gapSet(Factors((PermRootOps[:BadNumber])(W)))
        if Size(FundamentalGroup(W)) > 1
            print("# Nr classes each type_ implemented only for simply connected groups")
            return
        end
        for r = C[:ss]
            if !(haskey(r, :nrClasses))
                HF = r[:CGs]
                H = Group(HF)
                P = deepcopy(ClosedSubsets(W))
                o = Filtered(1:Size(P), (i->begin
                                OnSets((P[:elements])[i], HF[:phi]) == (P[:elements])[i]
                            end))
                o = Filtered(o, (i->begin
                                all((j->begin
                                            j in (P[:elements])[i]
                                        end), (H[:rootInclusion])[H[:generatingReflections]])
                            end))
                P = Restricted(P, o)
                P[:elements] = (P[:elements])[o]
    # here P poset of HF.phi-stable closed subsets containing roots(H)
                InfoChevie("# ", HF, "==>", P, "")
                l = map((x->begin
                                Spets(ReflectionSubgroup(W, x), HF[:phi])
                            end), P[:elements])
                l = map(function (RF,)
                            local res, d, R
                            R = Group(RF)
                            res = Product(Filtered(ReflectionDegrees(RF), (y->begin
                                                y[1] == 1
                                            end)), (p->begin
                                            Mvp("q") - p[2]
                                        end))
                            if R[:semisimpleRank] == 0
                                return res
                            end
                            d = SmithNormalFormMat((R[:roots])[R[:generatingReflections]] * R[:simpleRoots])
                            d = Filtered(Flat(d), (x->begin
                                            !x in [0, 1, C[:p]]
                                        end))
                            return res * Product(d, function (i,)
                                            if i == 2 && (2 in b && C[:p] != 2)
                                                return 2
                                            end
                                            if mod(C[:p], i) == 1
                                                return i
                                            end
                                            return Mvp(SPrint("q_", i))
                                        end)
                        end, l)
                less = (i->begin
                            Difference(ListBlist(1:length(l), (Incidence(P))[i]), [i])
                        end)
                o = LinearExtension(P)
                mu = []
                mu[o[Size(P)]] = 1
                for i = o[Size(P) - 1:(Size(P) - 2) - (Size(P) - 1):1]
                    mu[i] = -(Sum(mu[less(i)]))
                end
                n = Stabilizer(W, gapSet((H[:rootInclusion])[H[:generatingReflections]]), OnSets)
                n = (mu * l) // Size(Centralizer(n, HF[:phi]))
                InfoChevie("==>", n, ":", Stime(), "\n")
                r[:nrClasses] = n
            end
        end
        C[:ss] = Filtered(C[:ss], (x->begin
                        x[:nrClasses] != 0
                    end))
        return map((x->begin
                        x[:nrClasses]
                    end), C[:ss])
    end
