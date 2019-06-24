
UnipotentClassOps=Dict{Symbol,Any}()

UnipotentClassOps[:Name]=function(u)
  s=replace(u[:name],r"_"=>"")
  s=replace(s,r"\\tilde "=>"~")
  s
end

function unipotent_classes(t::TypeIrred,p=0) 
  uc=getchev(t,:UnipotentClasses,p)
  rank=length(t[:indices])
  for u in uc[:classes] # fill omitted fields
    if !haskey(u,:parameter) u[:parameter]=u[:name] end
    if !haskey(u,[:dynkin])
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
      u[:dimred]=Dimension(u[:red])
      u[:dimunip]=2*u[:dimBu]+rank-u[:dimred]
    elseif haskey(u,:dimred)
      u[:dimunip]=2*u[:dimBu]+rank-u[:dimred]
    end
  end
  uc[:orderClasses]=convert.(Vector{Int},uc[:orderClasses])
  uc
end

function unipotent_classes(W::FiniteCoxeterGroup,p=0) 
  t=refltype(W) 
  uc=unipotent_classes.(t)
  if isempty(t) ucl=Dict{Symbol,Any}(:classes=>[Dict(:name=>"",:Au=>coxgroup(),
      :parameter=>[],:dimBu=>0,:dynkin=>[],:balacarter=>[],
      :dimunip=>0,:red=>Torus(W.rank),:operations=>UnipotentClassOps)])
  else
   ucl=Dict{Symbol,Any}(:classes=>map(Cartesian(getindex.(uc,:classes)...)) do v
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
      end)
  end
  ucl[:size]=length(ucl[:classes])
  if iszero(p) && !haskey(ucl[:classes][1],:balacarter)
    bc=BalaCarterLabels(W)
    for u in ucl[classes]
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
  #   if length(v)=1 then Inherit(s,v[1],Difference(RecFields(v[1]),
  #["levi","Z","locsys","parameter"]));fi;
      s
      end
# adjust indices of levi, relativetype so they agree with Parent(Group(WF))
  for s in ucl[:springerSeries] s[:levi]=inclusion.(Ref(W),s[:levi]) end
  ucl[:p]=p
  if length(uc)==1 Inherit(ucl,uc[1],setdiff(keys(uc[1]),
    [:group,:operations,:springerSeries,:classes,:orderClasses])) 
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
# ucl[:orderClasses]=Hasse(Restricted(Poset(ucl[:orderClasses]),l))
  ucl[:size]=length(l)
  ucl
end
