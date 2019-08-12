struct UnipotentClass
  name::String
  parameter::Any
  dimBu::Int
  prop::Dict{Symbol,Any}
end

struct UnipotentClasses
  classes::Vector{UnipotentClass}
  p::Int
  orderclasses::Poset
  springerseries::Vector{Dict}
  prop::Dict{Symbol,Any}
end
  
function Base.show(io::IO,uc::UnipotentClasses)
  print(io,"UnipotentClasses(",uc.prop[:spets],")")
end

function nameclass(u::Dict,opt=Dict{Symbol,Any}())
  if haskey(opt,:mizuno) && haskey(u.prop,:mizuno) n=u[:mizuno] 
  elseif haskey(opt,:shoji) && haskey(u.prop,:shoji) n=u[:shoji]
  else n=u[:name]
  end
  TeX=haskey(opt,:TeX) 
  if !TeX
   if haskey(opt,:limit) n=TeXstrip(n)
   else
     n=replace(n,r"\\tilde *"=>"~")
     n=replace(n,"_"=>"")
     n=replace(n,"}"=>"")
     n=replace(n,"{"=>"")
   end
  end
  if haskey(opt,:locsys) && opt[:locsys]!=charinfo(u[:Au])[:positionId]
    cl="("*CharNames(u[:Au],opt)[opt[:locsys]]*")"
    n*(TeX ? "^{$cl}" : cl)
  elseif haskey(opt,:class) && opt[:class]!=charinfo(u[:Au])[:positionId]
    cl=classinfo(u[:Au])[:classnames][opt[:class]]
    TeX ? "\\hbox{\$$n\$}_{($cl)}" : "$n_$cl"
  else n
  end
end

function name(u::UnipotentClass,opt=Dict{Symbol,Any}())
  nameclass(merge(u.prop,Dict(:name=>u.name)),opt)
end

function Base.show(io::IO,u::UnipotentClass)
  opt=Dict{Symbol,Any}()
  if get(io,:TeX,false) opt[:TeX]=true end
  if get(io,:limit,false) opt[:limit]=true end
  print(io,"UnipotentClass(",name(u,opt),")")
end

UnipotentClassOps=Dict{Symbol,Any}(:Name=>nameclass)

# h  is a  linear form  defined by  its value  on the  simple roots  of the
# reflection subgroup K. Induce it to W by extending by 0 on the orthogonal
# of K, then conjugate it so it takes >=0 values on the simple roots.
function InducedLinearForm(W,K,h)
# print("W=$W K=$K h=$h");
  if semisimplerank(K)==0 return fill(0,semisimplerank(W)) end
  h=vcat(h,zeros(Int,rank(K)-semisimplerank(K)))
  h=Int.(inv(Rational.(PermRoot.baseX(K.G)))*h)
  r=parent(W).G.roots[inclusion(W)]
  v=toM(r[1:W.N])*h
  w=with_inversions(W,filter(i->v[i]<0,1:W.N))
  map(i->r[i]*h,restriction.(Ref(W),
        inclusion.(Ref(W),eachindex(gens(W))).^(w^-1)))
end

function DistinguishedParabolicSubgroups(W)
  filter(combinations(inclusion.(Ref(W),eachindex(gens(W))))) do J
    if isempty(J) return true end
    p=fill(1,semisimplerank(W))
    p[restriction.(Ref(W),J)]=fill(0,length(J))
    p=toM(W.rootdec[1:W.N])*p
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

# QuotientAu(Au,chars): chars is a list of indices of characters of Au.
# If  k is the common kernel of chars, QuotientAu returns a record
# rec(Au:=Au/k,
#     chars:=index of chars as characters of Au/k,
#     gens:=words in Au preimages of generators of Au/k)
# Since  GAP3  has  many  problems  with  quotient groups, we are forced to
# program an ad hoc solution which works only for Au actually occuring for
# unipotent classes of a reductive group G.
function QuotientAu(Au,chars)
  AbGens=function(g)
    res=[]
    l=gens(g)
    while !isempty(l)
      sort!(l,by=x->-Order(g,x))
      t=l[1].*elements(subgroup(g,res))
      if any(x->Order(g,x)<Order(g,l[1]),t)
        t=First(t,x->Order(g,x)<Order(g,l[1]))
	if Order(g,t)>1 l[1]=t
	else l=Drop(l,1)
        end
      else Add(res,l[1])
      end
    end
    res
  end
  # q=Au/k,  ww=words in q images of gens(Au)
  finish=function(q,ww)
    h=GroupHomomorphismByImages(Au,q,gens(Au),map(ww,x->EltWord(q,x)))
    fusion=List(classinfo(Au)[:classtext],
       c->position_class(q,Image(h,EltWord(Au,c))))
    ctu=chartable(Au).irr
    cth=chartable(q).irr
    return Dict(:Au=>q,:chars=>map(
      c->Position(cth,map(j->ctu[c,findfirst(isequal(j),fusion)],
                          1:NrConjugacyClasses(q))),chars),
      :gens=>map(x->GetWord(Au,First(Elements(Au),y->Image(h,y)==x)),gens(q)))
  end
  Z=n->ComplexReflectionGroup(n,1,1)
  ct=permutedims(chartable(Au).irr[chars,:])
  cl=filter(i->ct[i,:]==ct[1,:],axes(ct,1))
# println("Au=$Au chars=$chars ct=$ct cl=$cl")
  if length(cl)==1 return Dict(:Au=>Au,:chars=>chars,
                              :gens=>map(x->[x],eachindex(gens(Au)))) end
  ct=permutedims(toM(unique(sort(toL(ct)))))
# println("ct=$ct")
# k=Subgroup(Au,filter(x->position_class(Au,x) in cl,elements(Au)))
  k=Group(filter(x->position_class(Au,x) in cl,elements(Au)))
  if length(k)==length(Au) return Dict(:Au=>coxgroup(),:chars=>[1],:gens=>[])
  end
  if semisimplerank(Au)==1 return finish(Z(div(length(Au),length(k))),[[1]])
  elseif IsAbelian(Au/k)
    q=Au/k
    q.generators=AbGens(q)
    h=NaturalHomomorphism(Au,q)
    f=List(gens(Au),x->GetWord(q,x^h))
#  Print(Product(List(gens(q),x->Z(Order(q,x))))," ",f,"\n");
    return finish(Product(List(gens(q),x->Z(Order(q,x)))),f)
  else
    p=PositionProperty(t->all(x->x in ReflectionSubgroup(Au,t.indices)),
                       elements(k),refltype(Au))
    if p!=false
      p=refltype(Au)[p].indices
      if length(k)==length(reflection_subgroup(Au,p))
	return finish(reflection_subgroup(Au,Difference(gens(Au),p)),
                      map(i->i in p ? [] : [i],gens(Au)))
      elseif length(p)==1
        t=copy(refltype(Au))
        p=findfirst(t->t.indices==p,refltype(Au))
	t[p].p/=length(k)
        return finish(ReflectionGroup(t...),map(x->[x],gens(Au)))
      end
    elseif ReflectionName(Au)=="A1xB2" && length(k)==2 && longest(Au) in k
      return finish(coxgroup(:B,2),[[1,2,1,2],[1],[2]])
    end
  end
# Print(" Au=",ReflectionName(Au)," sub=",List(k.generators,e.Get),"\n");
  error("not implemented ",ReflectionName(Au),chars)
# q:=Au/k; f:=FusionConjugacyClasses(Au,q); Print(" quot=",q," fusion=",f,"\n");
# return rec(Au:=Au,chars:=chars);
end

# When some Springer series have been suppressed/weeded out, we  quotient 
# the Aus by the common  kernel of the remaining characters of the Aus. 
function AdjustAu!(classes,springerseries)
  for (i, u) in enumerate(classes)
    l=map(s->filter(k->s[:locsys][k][1]==i,eachindex(s[:locsys])),
          springerseries)
#   println(springerseries)
    chars=vcat(map(j->last.(springerseries[j][:locsys][l[j]]),
                   eachindex(l))...)
    f=QuotientAu(u.prop[:Au],chars)
#   if Size(Au)<>Size(f.Au) then
#     Print("class ",i,"=",classes[i].name," ",[Au,chars],"=>",f,"\n");
#   fi;
    u.prop[:Au]=f[:Au]
    if haskey(u.prop,:AuAction)
      R=u.prop[:AuAction].group
      if rank(R)==0 
        u.prop[:AuAction]=ExtendedCox(R,map(x->fill(0,0,0),f[:gens]))
      else 
       if isempty(f[:gens]) F0s=[matX(R,R())]
       else F0s=map(x->prod(u.prop[:AuAction].F0s[x]),f[:gens])
       end
       u.prop[:AuAction]=ExtendedCox(R,F0s)
      end
#     u[:AuAction].phis=map(x->Product(u[:AuAction].phis[x]),f[:gens])
    end
    k=1
    for j in eachindex(l)
      springerseries[j][:locsys]=copy(springerseries[j][:locsys])
      for s in l[j] 
        springerseries[j][:locsys][s][2]=f[:chars][k]
        k+=1 
      end
    end
  end
end

function UnipotentClasses(t::TypeIrred,p=0) 
  uc=getchev(t,:UnipotentClasses,p)
  rank=length(t[:indices])
  classes=UnipotentClass[]
  for u in uc[:classes] # fill omitted fields
    name=u[:name]
    parameter= haskey(u,:parameter) ? u[:parameter] : u[:name]
    dimBu= haskey(u,:dimBu)  ? u[:dimBu] : -1
    if haskey(u,:dynkin)
      weights=toM(roots(cartan(t.prop)))*u[:dynkin]
      n0=count(iszero,weights)
      if dimBu==-1 dimBu=n0+div(count(isone,weights),2)
      elseif dimBu!=n0+div(count(isone,weights),2) error("theory")
      end
      n0=2*n0-count(isequal(2),weights)
      u[:dimunip]=2*dimBu-n0
      u[:dimred]=n0+rank
    elseif haskey(u,:red)
      u[:dimred]=dimension(u[:red])
      u[:dimunip]=2*dimBu+rank-u[:dimred]
    elseif haskey(u,:dimred)
      u[:dimunip]=2*dimBu+rank-u[:dimred]
    end
    delete!.(Ref(u),[:name,:parameter,:dimBu])
    push!(classes,UnipotentClass(name,parameter,dimBu,u))
  end
  springerseries=uc[:springerSeries]
  for s in springerseries
    if isempty(s[:levi]) s[:levi]=Int[] end
    s[:locsys]=Vector{Int}.(s[:locsys])
  end
  orderclasses=map(x->isempty(x) ? Int[] : x,uc[:orderClasses])
  delete!.(Ref(uc),[:classes,:orderClasses,:springerSeries])
  uc[:spets]=t
  UnipotentClasses(classes,p,Poset(orderclasses),springerseries,uc)
end

Base.length(uc::UnipotentClasses)=length(uc.classes)

function UnipotentClasses(W::FiniteCoxeterGroup,p=0) 
  t=refltype(W) 
  uc=UnipotentClasses.(t,p)
  if isempty(t) 
    classes=[UnipotentClass("",[],0,
        Dict(:Au=>coxgroup(),:dynkin=>[],:balacarter=>[],
             :dimunip=>0,:red=>Torus(rank(W))))]
    uc=UnipotentClasses(classes,p,[Int[]],
      [Dict(:Z=>Int[],:levi=>Int[],:locsys=>[[1,1]],:relgroup=>coxgroup())],
     Dict{Symbol,Any}(:spets=>W))
  else
    classes=map(Cartesian(map(x->x.classes,uc)...)) do v
      l=getindex.(t,:indices)
#     println("v=$v")
      if length(v)==1 u=deepcopy(v[1]) 
      else
        u=UnipotentClass(join(map(x->x.name,v),","),map(x->x.parameter,v),
                         sum(map(x->x.dimBu,v)),Dict{Symbol,Any}())
        u.prop[:Au]=prod(x->x.prop[:Au],v)
        if all(x->haskey(x.prop,:dimred),v) 
          u.prop[:dimred]=sum(x->x.prop[:dimred],v) end
        if all(x->haskey(x.prop,:dimunip),v) 
          u.prop[:dimunip]=sum(x->x.prop[:dimunip],v) end
        if all(x->haskey(x.prop,:red),v) 
          u.prop[:red]=prod(x->x.prop[:red],v) end
        if all(x->haskey(x.prop,:AuAction)) 
          u.prop[:AuAction]=prod(x->x.prop[:AuAction],v) end
        if all(x->haskey(x.prop,:dynkin))
          u.prop[:dynkin]=zeroes(Int,sum(x->length(x.prop[:dynkin]),v))
          for i in 1:length(l) u.prop[:dynkin][l[i]]=v[i].prop[:dynkin] end
        end
      end
      if all(x->haskey(x.prop,:balacarter),v)
        u.prop[:balacarter]=vcat([map(j->j>0 ? x[j] : -x[-j],
                    v[i].prop[:balacarter]) for (i,x) in enumerate(l)]...)
      end
      if rank(W)>semisimplerank(W) && haskey(u.prop, :red) 
        T=torus(rank(W)-semisimplerank(W))
        u.prop[:red]*=T
        if haskey(u.prop,:AuAction)
          u.prop[:AuAction]=ExtendedCox(u.prop[:AuAction].group*T,
               map(x->DiagonalMat(x,matX(T,T())),u.prop[:AuAction].F0s))
        end
      end
      u
    end
  end
  if iszero(p) && !haskey(classes[1].prop,:balacarter)
    bc=BalaCarterLabels(W)
    for u in classes
      u.prop[:balacarter]=bc[findfirst(p->p[1]==u.prop[:dynkin],bc)][2]
    end
  end
  ll=map(length,uc)
  orderclasses=map(Cartesian(map(x->1:x,ll)...)) do v
    o=Cartesian(map(j->vcat(hasse(uc[j].orderclasses)[v[j]],[v[j]]),1:length(v))...)
    o=map(x->PositionCartesian(ll,x),o)
    setdiff(o,[PositionCartesian(ll,v)])
  end
  springerseries=map(Cartesian(map(x->x.springerseries,uc)...)) do v
    if isempty(v) return Dict(:Z=>[],:levi=>[],:locsys=>[[1,1]])
    elseif length(v)==1 return deepcopy(v[1])
    end
    s=Dict(:levi=>vcat(map(i->l[i][v[i][:levi]],eachindex(v))...))
    s[:Z]=vcat(getindex.(v,:Z)...)
    s[:locsys]=map(Cartesian(getindex.(v,:locsys))) do v
        v=permutedims(v)
        v[2]=PositionCartesian(List(i->NrConjugacyClasses(
              uc[i].classes[v[1][i]][:Au]),eachindex(v[1])),v[2])
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
  for s in springerseries s[:levi]=inclusion.(Ref(W),s[:levi]) end
  if length(uc)==1 prop=uc[1].prop else prop=Dict{Symbol,Any}() end
  prop[:spets]=W
# To deal with a general group intermediate between Gad and Gsc, we discard
# the  Springer series  corresponding to  a central  character which is not
# trivial on the fundamental group (seen as a subgroup of ZGsc)
# AlgebraicCentre(W).descAZ returns the generators of the fundamental group
# of  the  algebraic  group  W  as  words  in  generators  of  the absolute
# fundamental group.
  if !all(x->Set(x[:Z])==Set([1]),springerseries)
    springerseries=filter(s->all(y->Product(s[:Z][y])==1,
             AlgebraicCentre(W)[:descAZ]),springerseries)
    AdjustAu!(classes,springerseries) 
  end
  s=springerseries[1]
# s[:relgroup]=RelativeCoset(WF,s[:levi])
# s[:locsys]=s[:locsys][charinfo(s[:relgroup])[:charRestrictions]]
  l=filter(i->any(y->i==y[1],s[:locsys]),1:length(classes))
  s[:locsys]=map(y->[Position(l,y[1]),y[2]],s[:locsys])
  # for now only springerseries[1] properly twisted
  for s in springerseries[2:end] 
 #  s[:relgroup]=RelativeCoset(WF,s[:levi])
 #  s[:locsys]=s[:locsys][charinfo(s[:relgroup])[:charRestrictions]]
    s[:locsys]=map(y->[Position(l,y[1]),y[2]],s[:locsys])
  end
  classes=classes[l]
  AdjustAu!(classes,springerseries)
  orderclasses=Poset(hasse(restricted(Poset(orderclasses),l)))
  merge!(orderclasses.prop,Dict(:classes=>classes,
  :label=>function(p,n,opt)
    name(p.prop[:classes][n],merge(opt,Dict(:limit=>true)))
   end))
  ucl=UnipotentClasses(classes,p,orderclasses,springerseries,prop)
  ucl
end

function FormatCentralizer(u,opt)
  c=""
  function AuName(u)
    if length(u.prop[:Au])==1 return "" end
    res=haskey(u.prop,:AuAction) || 
        (haskey(u.prop,:dimred) && iszero(u.prop[:dimred])) ? "." : "?"
    au=reflection_name(u.prop[:Au],opt)
    if haskey(opt,:TeX) 
      rep=["A_1"=>"Z_2","A_2"=>"S_3","A_3"=>"S_4","A_4"=>"S_5","B_2"=>"D_8"]
    else
      rep=["A1"=>"Z2","A2"=>"S3","A3"=>"S4","A4"=>"S5","B2"=>"D8"]
    end
    for p in rep  au=replace(au,p) end
    res*au
  end
  if haskey(u.prop,:dimunip)
    if u.prop[:dimunip]>0 c*=Format(Mvp(:q)^u.prop[:dimunip],opt) end
  else c*="q^?" end
  if haskey(u.prop,:AuAction)
    if Rank(u.prop[:red])>0
      c*="."
      if length(u.prop[:Au])==1 || length(u.prop[:Au])==length(Group(u.prop[:AuAction].F0s...))
        c*=reflection_name(u.prop[:AuAction],opt)
      elseif all(isone,u.prop[:AuAction].F0s)
        c*=reflection_name(u.prop[:AuAction].group,opt)*AuName(u)
      else
        c*=reflection_name(u.prop[:AuAction],opt)*AuName(u)
      end
    else
      c*=AuName(u)
    end
  elseif haskey(u.prop,:red)
    n=reflection_name(u.prop[:red],opt)
    if n!="." c*="."*n end
    c*=AuName(u)
  else
    if haskey(u.prop,:dimred)
      if u.prop[:dimred]>0 c*="[red dim ",u.prop[:dimred],"]" end
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
 
function Util.format(io::IO,uc::UnipotentClasses, opt=Dict{Symbol,Any}())
  opt = merge(UnipotentClassesOps[:DisplayOptions],opt)
  TeX(a, b)=haskey(opt, :TeX) ? a : b
  opt[:rowLabels]=name.(uc.classes,Ref(merge(opt,Dict(:limit=>true))))
  if haskey(opt,:order)
    println(Posets.showgraph(uc.orderclasses;opt...))
  end
  sp = map(copy, uc.springerseries)
  if haskey(opt, :fourier)
    for p in sp p[:locsys] = p[:locsys][DetPerm(p[:relgroup])] end
  end
  W = uc.prop[:spets]
  if uc.p!=0 || !any(x->haskey(x.prop,:balacarter),uc.classes)
     delete(opt,:balacarter)
  end
  tbl = map(uc.classes)do u
    res= iszero(uc.p) ? [joindigits(u.prop[:dynkin])] : String[]
    push!(res, string(u.dimBu))
    if opt[:balaCarter]
      if haskey(u.prop, :balacarter)
        b=fill('.',coxrank(W))
        for i in u.prop[:balacarter] if i>0 b[i]='2' else b[-i]='0' end end
      else
        b=fill('?',coxrank(W))
      end
      push!(res, String(b))
    end
    if opt[:centralizer] push!(res, FormatCentralizer(u, opt)) end
    if opt[:springer]
        i = Position(uc.classes, u)
        res = Append(res, map((ss->begin
           join(map(function (i)
                c1 = CharNames(u.prop[:Au], opt)[ss[:locsys][i][2]]
                c2 = CharNames(ss[:relgroup], opt)[i]
                c1=="" ? c2 : c1*":"*c2
           end, findall(y->y[1]==i,ss[:locsys])), TeX("\\kern 0[:8]em "," "))
       end), sp))
    end
    res
  end
  column_labels = String[]
  if iszero(uc.p)
      push!(column_labels, TeX("\\hbox{Dynkin-Richardson}", "D-R"))
  end
  push!(column_labels, TeX("\\dim{\\cal B}_u", "dBu"))
  if opt[:balaCarter]
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
      p = Perm(sortperm(map(x->x.dimBu, uc.classes)))
      tbl = Permuted(tbl, p)
      opt[:rowLabels] = Permuted(opt[:rowLabels], p)
  end
  format(io,toM(tbl),rows_label="u",row_labels=opt[:rowLabels],
         column_labels=column_labels)
end

# decompose tensor product of characteres (given as their indices in chartable)
function DecomposeTensor(W,c::Int...)
  ct=chartable(W)
  irr=ct.irr
# println("eltype=",eltype(irr))
  ch=conj.(prod(irr[collect(c),:],dims=1))
  sW=length(W)
  classes=map(x->div(sW,x),ct.centralizers)
  ch=map(*,ch,classes)
  convert.(Int,div.(irr*ch,sW))
end

function Det(m)
  function compl(m,i,j)
    v=axes(m,1)
    m[setdiff(v,[i]),setdiff(v,[j])]
  end
  if isempty(m) return 0 end
  n=size(m,1)
  if n<=3 
    return sum(p->prod(i->m[i,i^p],1:n)*sign(p),collect(symmetric_group(n)))
  end
  i=findfirst(i->count(x->!iszero(x),m[i,:])<=2,axes(m,1))
  if !isnothing(i)
    j=findall(x->!iszero(x),m[i,:])
    if isempty(j) return 0 end
    return sum(k->(-1)^(i+k)*m[i,k]*Det(compl(m,i,k)),j)
  end
  v=axes(m,1)
  if length(v)<=3 return det*Det(m) end
  for j in v
    i=findfirst(isunit,m[v,j])
    if isnothing(i) continue end
    println(size(m,1),":",[j,i])
    m=copy(m)
    f=inv(m[i,j])
    for k in setdiff(v,[i]) m[k,:]-=m[k][j]*f.*m[i,:] end
    return det*(-1)^(i+j)*m[i][j]*Det(compl(m,i,j))
  end
  print("m=");display(iszero.(m))
  v=map(x->count(y->!iszero(y),x),m)
  j=Position(v,Minimum(v))
  if Length(m)>71 println("\n",Length(m),":",v[j]) end
  if v[j]>5 return det*DeterminantMat(m) end
  det*Sum(1:Length(m),function(i)
#   if Length(m) mod 5=0 then Print(Length(m),":",i," \c");fi;
    if iszero(m[i][j]) return 0
    else return (-1)^(i+1)*m[i][j]*Det(compl(m,i,j))
    end end)
end

# Cofactors of the square matrix M, defined by CoFactors(M)*M=Det(M)*one(M)
function CoFactors(m)
  if size(m,1)==1 return fill(1,1,1) end
  v=axes(m,1)
  permutedims([(-1)^(i+j)*Det(m[filter(x->x!=i,v),filter(x->x!=j,v)])
               for i in v, j in v])
end

#############################################################################
##
# BigCellDecomposition(M [, b]) . . . .  Decompose in the big Bruhat cell
# M should be a square matrix such that the principal minors which are
# union of blocks should be non-zero.
# The  function decomposes  M as  a product  P*L*tP where  P is lower block
# unitriangular  and tp  upper block-unitriangular  (with identity diagonal
# blocks)  and L block-diagonal according to the  block structure b; b is a
# list  of lists of  union [1..Length(M)]. If  not given, the trivial block
# structure [[1],..,[Length(M)]] is assumed.
# If  M is  symmetric then  tP=TransposedMat(P) and  the result is the pair
# [tP,L]. else the result is [P,L,tP]
#
function BigCellDecomposition(M,b=map(i->[i],axes(M,1)))
  L=one(M)
  P=one(M)
  block(X,i,j)=X,b[i],b[j]
  if M==permutedims(M)
    for j in eachindex(b)
      L[b[j],b[j]]=block(M,j,j)
      if j>1 L[b[j],b[j]]-=sum(k->block(P,j,k)*block(L,k,k)*
                              permutedims(block(P,j,k)),1:j-1) end
      cb=CoFactors(block(L,j,j))
      db=Det(block(L,j,j))
      for i in j+1:length(b)
        P[b[i],b[j]]=block(M,i,j)
        if j>1 P[b[i],b[j]]-=
          sum(k->block(P,i,k)*block(L,k,k)*permutedims(block(P,j,k)),1:j-1) 
        end
        P[b[i],b[j]]*=cb
        P[b[i],b[j]]=div.(block(P,i,j),db)
      end
    end
    return permutedims(P),L
  end
  tP=one(M)
  for j in eachindex(b)
    L[b[j],b[j]]=block(M,j,j)-sum(k->block(P,j,k)*block(L,k,k)*
                                  permutedims(block(P,j,k)),1:j-1)
    cb=CoFactors(block(L,j,j))
    db=Det(block(L,j,j))
    for i in j+1:length(b)
      P[b[i],b[j]]=block(M,i,j)-sum(k->block(P,i,k)*block(L,k,k)*
                                     block(tP,k,j),1:j-1)
      P[b[i],b[j]]=div.(P[b[i],b[j]]*cb,db)
      tP[b[j],b[i]]=cb*(block(M,j,i)-sum(k->block(P,j,k)*
                                         block(L,k,k)*block(tP,k,i),1:j-1))
      tP[b[i],b[j]]=div.(tP[b[i],b[j]],db)
    end
  end
  (P,L,tP)
end

# coefficients of R_\chi on unipotently supported local systems 
# ICCTable(uc[,Springer series no[,variable]]) eg (uc,1,X(Rationals))
# Works for G split.
function ICCTable(uc,i=1,var=Pol(:q))
  W=uc.prop[:spets] # W=Group(uc.spets)
  ss=uc.springerseries[i]
  res=Dict(:spets=>uc.prop[:spets],:relgroup=>ss[:relgroup],
           :series=>i,:q=>var,:p=>uc.p)
  if haskey(ss,:warning) println("# ",ss[:warning])
    res[:warning]=ss[:warning]
  end
# We are going to solve the equation in "unipotent support", page 151
# $Transposed(P)\Lambda P=\omega$
# where $\Lambda_{i,j}$ is  $\sum_{g\in G^F} Y_i(g)\overline{Y_j(g)}$
# and $\Omega_{i,j}$ is equal to
# $|Z^0(G^F)|q^{-\text{semisimple rank}L}|G^F|/P(W_G(L))
#  q^{-b_i-b_j}FakeDegree(\chi_i\otimes\chi_j\otimes\sgn)$
# where $P(W_G(L))$ is the Poincare polynomial $\prod_i(q^{d_i}-1)$
# where $d_i$ are the reflection degrees of $W_G(L)$
# res[:scalar] is the masrix $P$
  R=ss[:relgroup]
  ct=chartable(R)
  q=Pol(:q)
  f=fakedegrees(R,q)
  k=charinfo(R)[:positionDet]
  n=length(f)
# Partition on characters of ss.relgroup induced by poset of unipotent classes
  res[:dimBu]=map(x->uc.classes[x[1]].dimBu,ss[:locsys])
  res[:blocks]=CollectBy(eachindex(ss[:locsys]),-res[:dimBu])
  # matrix of q^{-b_i-b_j}*fakedegree(χ_i ⊗ χ_j ⊗ sgn)
  tbl=BigCellDecomposition(
    [q^(-res[:dimBu][i]-res[:dimBu][j])*f*DecomposeTensor(R,i,j,k) 
     for i in 1:n,j in 1:n], res[:blocks])
  res[:scalar]=tbl[1]
  res[:locsys]=ss[:locsys]
# res[:L]=tbl[2]*GenericOrder(W,q)/Product(ReflectionDegrees(R),d->q^d-1)/
#   q^(W.semisimpleRank-R.semisimpleRank);
  res[:L]=tbl[2]*q^(W.N+semisimplerank(R)-semisimplerank(W))
  res[:uc]=uc
  if haskey(ss,:parameter) res[:parameter]=ss[:parameter]
  else res[:parameter]=1:length(ss[:locsys])
  end
  if var!=q
    q=var
    res[:scalar]=map(x->x(q),res[:scalar])
    res[:L]=map(x->x(q),res[:L])
  end
  res
end

function formatICC(x,opt=Dict())
  text="Coefficients of \$X_\\phi\$ on \$Y_\\psi\$ for \$"*
        reflection_name(x[:relgroup],opt)*"\$\n"
  if haskey(opt,:TeX) text*="\\medskip\n\n"
  else text=TeXstrip(text) end
  print(stdout,text)
  if !haskey(opt,:columns) && !haskey(opt,:rows)
    opt[:rows]=collect(1:length(x[:dimBu]))
    sort!(opt[:rows],by=i->[x[:dimBu][i],x[:locsys][i]])
    opt[:columns]=opt[:rows]
  end
  tbl=copy(x[:scalar])
  if !haskey(opt,:CycPol) opt[:CycPol]=true end
  if opt[:CycPol] tbl=map(CycPol,tbl) end
  tbl=repr.(tbl,context=:limit=>true)
  columnLabels=map(p->name(x[:uc].classes[p[1]],merge(Dict(:locsys=>p[2]),opt)),
                   x[:locsys])
  rowLabels=map(x->haskey(opt,:TeX) ? "X_{$x}" : "X$x",
                CharNames(x[:relgroup],opt))
  format(stdout,permutedims(tbl),rows=opt[:rows],columns=opt[:columns],
         row_labels=rowLabels,column_labels=columnLabels)
end
