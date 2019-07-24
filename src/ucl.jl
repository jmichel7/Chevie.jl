UnipotentClassOps=Dict{Symbol,Any}()

function uclassname(u,opt=Dict{Symbol,Any}())
  if haskey(opt,:mizuno) && haskey(u,:mizuno) n=u[:mizuno]
  elseif haskey(opt,:shoji) && haskey(u,:shoji) n=u[:shoji]
  else n=u[:name]
  end
  if !haskey(opt,:TeX) 
   n=replace(n,r"\\tilde *"=>"~")
   n=replace(n,"_"=>"")
   n=replace(n,"}"=>"")
   n=replace(n,"{"=>"")
  end
  if haskey(opt,:locsys)
   if opt[:locsys]==PositionId(cl[:Au]) return n end
   cl=SPrint("(",CharNames(cl[:Au],opt)[opt.locsys],")")
    if haskey(opt,:TeX) return SPrint(n,"^{",cl,"}")
    else return SPrint(n,cl) end
  elseif haskey(opt,:class)
   if opt[:class]==PositionId(cl[:Au]) return n end
   cl=ChevieClassInfo(cl[:Au]).classnames[opt.class]
    if haskey(opt,:TeX) return SPrint("\\hbox{\$",n,"\$}_{(",cl,")}")
    else return SPrint(n,"_",cl) end
  else return n
  end
end
UnipotentClassOps[:Name]=uclassname

# h  is a  linear form  defined by  its value  on the  simple roots  of the
# reflection subgroup K. Induce it to W by extending by 0 on the orthogonal
# of K, then conjugate it so it takes >=0 values on the simple roots.
function InducedLinearForm(W,K,h)
# print("W=$W K=$K h=$h");
  if semisimplerank(K)==0 return fill(0,semisimplerank(W)) end
  h=copy(h)
  append!(h,fill(0,rank(K)-semisimplerank(K)))
  h=Int.(PermRoot.baseX(parent(W.G))/Rational.(PermRoot.baseX(K.G))*h)
  r=parent(W).rootdec[inclusion(W)]
  v=permutedims(h)*hcat(r[1:W.N]...)
  w=with_inversions(W,filter(i->v[i]<0,1:W.N))
  map(i->r[i]*h,restriction.(Ref(W),
        inclusion.(Ref(W),eachindex(gens(W))).^(w^-1)))
end

function DistinguishedParabolicSubgroups(W)
  filter(combinations(inclusion.(Ref(W),eachindex(gens(W))))) do J
    if isempty(J) return true end
    p=fill(1,semisimplerank(W))
    p[restriction.(Ref(W),J)]=fill(0,length(J))
    p=permutedims(p)*hcat(W.rootdec[1:W.N]...)
    2*count(iszero,p)+semisimplerank(W)==count(isone,p)
  end
end

function BalaCarterLabels(W)
  l=vcat(map(parabolic_representatives(W)) do J
    map(D->[J,D],DistinguishedParabolicSubgroups(reflection_subgroup(W,J)))
  end...)
  map(l) do p
    L=reflection_subgroup(W,p[1])
    w=fill(2,semisimplerank(L))
    w[restriction(L,p[2])]=p[2]*0
    [InducedLinearForm(W,L,w),map(i->i in p[2] ? -i : i,p[1])]
  end
end

function unipotent_classes(t::TypeIrred,p=0) 
  uc=getchev(t,:UnipotentClasses,p)
  rank=length(t[:indices])
  for u in uc[:classes] # fill omitted fields
    if !haskey(u,:parameter) u[:parameter]=u[:name] end
    if haskey(u,:dynkin)
      weights=permutedims(u[:dynkin])*hcat(roots(cartan(t.prop))...)
      p=count(iszero,weights)
      if haskey(u,:dimBu) && u[:dimBu]!=p+div(count(isone,weights),2)
        error("theory")
      end
      u[:dimBu]=p+div(count(isone,weights),2)
      p=2*p-count(isequal(2),weights)
      u[:dimunip]=2*u[:dimBu]-p
      u[:dimred]=p+rank
    elseif haskey(u,:red)
      u[:dimred]=dimension(u[:red])
      u[:dimunip]=2*u[:dimBu]+rank-u[:dimred]
    elseif haskey(u,:dimred)
      u[:dimunip]=2*u[:dimBu]+rank-u[:dimred]
    end
  end
  uc[:orderClasses]=map(x->isempty(x) ? Int[] : x,uc[:orderClasses])
  for s in uc[:springerSeries]
   if isempty(s[:levi]) s[:levi]=Int[] end
  end
  uc
end

function unipotent_classes(W::FiniteCoxeterGroup,p=0) 
  t=refltype(W) 
  uc=unipotent_classes.(t,p)
  if isempty(t) ucl=Dict{Symbol,Any}(:classes=>[Dict(:name=>"",:Au=>coxgroup(),
      :parameter=>[],:dimBu=>0,:dynkin=>[],:balacarter=>[],
      :dimunip=>0,:red=>Torus(W.rank),:operations=>UnipotentClassOps)])
  else
    ucl=Dict{Symbol,Any}()
    ucl[:classes]=map(Cartesian(getindex.(uc,:classes)...)) do v
      l=getindex.(t,:indices)
      if length(v)==1 u=deepcopy(v[1]) 
      else
        u=Dict(:name=>join(getindex.(v,:name),","),:Au=>prod(x->x[:Au],v),
          :dimBu=>sum(x->x[:dimBu],v), :parameter=>getindex.(v,:parameter))
        if all(haskey.(v,:dimred)) u[:dimred]=sum(x->x[:dimred],v) end
        if all(haskey.(v,:dimunip)) u[:dimunip]=sum(x->x[:dimunip],v) end
        if all(haskey.(v,:red)) u[:red]=prod(x->x[:red],v) end
        if all(haskey.(v,:AuAction)) u[:AuAction]=prod(x->x[:AuAction],v) end
        if all(haskey.(v,:dynkin))
          u[:dynkin]=zeroes(Int,sum(x->length(x[:dynkin]),v))
          for i in 1:length(l) u[:dynkin][l[i]]=v[i][:dynkin] end
        end
      end
      if all(haskey.(v,:balacarter))
        u[:balacarter]=vcat(map(
         i->map(j->j>0 ? l[i][j] : -l[i][-j],v[i][:balacarter]),1:length(l))...)
      end
      if rank(W)>semisimplerank(W) && haskey(u, :red) 
        tr=rank(W)-semisimplerank(W)
        u[:red]*=torus(tr)
        if haskey(u,:AuAction)
          u[:AuAction][:group]*=torus(tr)
          u[:AuAction][:F0s]=map(x->DiagonalMat(x,IdentityMat(tr)),
                                 u[:AuAction][:F0s])
        end
      end
      u[:operations]=UnipotentClassOps
      u
    end
  end
  ucl[:size]=length(ucl[:classes])
  if iszero(p) && !haskey(ucl[:classes][1],:balacarter)
    bc=BalaCarterLabels(W)
    for u in ucl[:classes]
      u[:balacarter]=bc[findfirst(p->p[1]==u[:dynkin],bc)][2]
    end
  end
  ll=map(x->length(x[:classes]),uc)
  ucl[:orderClasses]=map(Cartesian(map(x->1:x,ll)...)) do v
    o=Cartesian(map(j->vcat(uc[j][:orderClasses][v[j]],[v[j]]),1:length(v))...)
    o=map(x->PositionCartesian(ll,x),o)
    setdiff(o,[PositionCartesian(ll,v)])
  end
  ucl[:springerSeries]=map(Cartesian(getindex.(uc,:springerSeries)...)) do v
    if isempty(v) return Dict(:Z=>[],:levi=>[],:locsys=>[[1,1]])
    elseif length(v)==1 return Copy(v[1])
    end
    s=Dict(:levi=>vcat(map(i->l[i][v[i][:levi]],eachindex(v))...))
    s[:Z]=vcat(getindex.(v,:Z)...)
    s[:locsys]=map(Cartesian(getindex.(v,:locsys))) do v
        v=permutedims(v)
        v[2]=PositionCartesian(List(i->NrConjugacyClasses(
              uc[i][:classes[v[1][i]]][:Au]),eachindex(v[1])),v[2])
        v[1]=PositionCartesian(ll,v[1])
        v
        end
    if all(haskey.(v,:parameter)) s[:parameter]=getindex.(v,:parameter) end
    if length(v)==1 
      for k in setdiff(keys(v[1]),[:levi,:Z,:locsys,:parameter])
        s[k]=v[1][k]
      end
    end
    s
  end
# adjust indices of levi, relativetype so they agree with Parent(Group(WF))
  for s in ucl[:springerSeries] s[:levi]=inclusion.(Ref(W),s[:levi]) end
  ucl[:p]=p
  ucl[:spets]=W
  if length(uc)==1 
    for k in setdiff(keys(uc[1]),
      [:group,:operations,:springerSeries,:classes,:orderClasses])
      ucl[k]=uc[1][k]
    end
  end
# To deal with a general group intermediate between Gad and Gsc, we discard
# the  Springer series  corresponding to  a central  character which is not
# trivial on the fundamental group (seen as a subgroup of ZGsc)
# AlgebraicCentre(W).descAZ returns the generators of the fundamental group
# of  the  algebraic  group  W  as  words  in  generators  of  the absolute
# fundamental group.
# if !all(x->Set(x[:Z])==Set([1]),ucl[:springerSeries])
#   ucl[:springerSeries]=filter(s->all(y->Product(s[:Z][y])==1,
#            AlgebraicCentre(W)[:descAZ]),ucl[:springerSeries])
#   ucl=AdjustAu(ucl) 
# end
  s=ucl[:springerSeries][1]
# s[:relgroup]=RelativeCoset(WF,s[:levi])
# s[:locsys]=s[:locsys][charinfo(s[:relgroup])[:charRestrictions]]
  l=filter(i->any(y->i==y[1],s[:locsys]),1:length(ucl[:classes]))
  s[:locsys]=map(y->[Position(l,y[1]),y[2]],s[:locsys])
  # for now only Springerseries[1] properly twisted
  for s in ucl[:springerSeries][2:end] 
 #  s[:relgroup]=RelativeCoset(WF,s[:levi])
 #  s[:locsys]=s[:locsys][charinfo(s[:relgroup])[:charRestrictions]]
    s[:locsys]=map(y->[Position(l,y[1]),y[2]],s[:locsys])
  end
  ucl[:classes]=ucl[:classes][l]
# AdjustAu(ucl)
  ucl[:orderClasses]=hasse(restricted(Poset(ucl[:orderClasses]),l))
  ucl[:size]=length(l)
  ucl
end

function FormatCentralizer(u,opt)
  c=""
  function AuName(u)
    if length(u[:Au])==1 return "" end
    res=haskey(u,:AuAction) || 
        (haskey(u,:dimred) && iszero(u[:dimred])) ? "." : "?"
    au=reflection_name(u[:Au],opt)
    if haskey(opt,:TeX) 
      rep=["A_1"=>"Z_2","A_2"=>"S_3","A_3"=>"S_4","A_4"=>"S_5","B_2"=>"D_8"]
    else
      rep=["A1"=>"Z2","A2"=>"S3","A3"=>"S4","A4"=>"S5","B2"=>"D8"]
    end
    for p in rep  au=replace(au,p) end
    res*au
  end
  if haskey(u,:dimunip)
    if u[:dimunip]>0 c*=Format(Mvp(:q)^u[:dimunip],opt) end
  else c*="q^?" end
  if haskey(u,:AuAction)
    if Rank(u[:red])>0
      c*="."
      if length(u[:Au])==1 || length(u[:Au])==length(Group(u[:AuAction].F0s...))
        c*=reflection_name(u[:AuAction],opt)
      elseif all(isone,u[:AuAction].F0s)
        c*=reflection_name(u[:AuAction].group,opt)*AuName(u)
      else
        c*=reflection_name(u[:AuAction],opt)*AuName(u)
      end
    else
      c*=AuName(u)
    end
  elseif haskey(u,:red)
    n=reflection_name(u[:red],opt)
    if n!="." c*="."*n end
    c*=AuName(u)
  else
    if haskey(u,:dimred)
      if u[:dimred]>0 c*="[red dim ",u[:dimred],"]" end
    else c*="[red??]"
    end
    c*=AuName(u)
  end
  replace(c,r"^."=>"")
  replace(c,r".*"=>"")
  replace(c,"()"=>"")
end

UnipotentClassesOps=Dict(:DisplayOptions=>Dict(
 :order=>true,:springer=>true,:centralizer=>true,:balaCarter=>true))

CharNames(W,opt=Dict{Symbol,Any}())=TeXstrip.(charinfo(W)[:charnames])
 
function formatuc(uc, opt=Dict{Symbol,Any}())
  opt = merge(UnipotentClassesOps[:DisplayOptions],opt)
  TeX(a, b)=haskey(opt, :TeX) ? a : b
  opt[:rowLabels] = TeXstrip.(uclassname.(uc[:classes],
                                          Ref(merge(opt,Dict(:TeX=>true)))))
  if haskey(opt,:order)
    println(Posets.showgraph(Poset(uc[:orderClasses]);opt...))
  end
  sp = map(copy, uc[:springerSeries])
  if haskey(opt, :fourier)
    for p in sp p[:locsys] = p[:locsys][DetPerm(p[:relgroup])] end
  end
  W = uc[:spets]
  tbl = map(uc[:classes])do u
    res= iszero(uc[:p]) ? [joindigits(u[:dynkin])] : String[]
    push!(res, string(u[:dimBu]))
    if iszero(uc[:p]) && opt[:balaCarter]
        if haskey(u, :balacarter)
            b = fill('.',coxrank(W))
            for i in filter(x->x>0,u[:balacarter]) b[i] = '2' end
            for i in filter(x->x<0,u[:balacarter]) b[-i] = '0' end
        else
            b = fill('?',coxrank(W))
        end
        push!(res, String(b))
    end
    if opt[:centralizer] push!(res, FormatCentralizer(u, opt)) end
    if opt[:springer]
        i = Position(uc[:classes], u)
        res = Append(res, map((ss->begin
           join(map(function (i)
                c1 = CharNames(u[:Au], opt)[ss[:locsys][i][2]]
                c2 = CharNames(ss[:relgroup], opt)[i]
                c1=="" ? c2 : c1*":"*c2
           end, findall(y->y[1]==i,ss[:locsys])), TeX("\\kern 0[:8]em "," "))
       end), sp))
    end
    res
  end
  column_labels = String[]
  if iszero(uc[:p])
      push!(column_labels, TeX("\\hbox{Dynkin-Richardson}", "D-R"))
  end
  push!(column_labels, TeX("\\dim{\\cal B}_u", "dBu"))
  if iszero(uc[:p]) && opt[:balaCarter]
      push!(column_labels, TeX("\\hbox{Bala-Carter}", "B-C"))
  end
  if opt[:centralizer]
      push!(column_labels, TeX("C_{\\bf G}(u)", "C(u)"))
  end
  if opt[:springer]
      column_labels = append!(column_labels, 
      map(function (ss,)
        res = string(repr(ss[:relgroup],context=:limit=>true),"(",
          repr(reflection_subgroup(W,ss[:levi]),context=:limit=>true),")")
        if !all(x->x==1,ss[:Z])
          res*=string("/", join(ss[:Z],","))
        end
        return res
    end, sp))
  end
  if !(haskey(opt, :rows))
      p = Perm(sortperm(map(x->x[:dimBu], uc[:classes])))
      tbl = Permuted(tbl, p)
      opt[:rowLabels] = Permuted(opt[:rowLabels], p)
  end
  format(stdout,permutedims(hcat(tbl...)),rows_label="u",
         row_labels=opt[:rowLabels],
         column_labels=column_labels)
end
