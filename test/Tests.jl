module Tests
using Gapjm

const test=Dict{Symbol,NamedTuple{(:fn, :applicable,:comment),
        Tuple{Function,Function,String}}}()

isweylgroup(W)=(W isa FiniteCoxeterGroup) && all(isinteger,cartan(W))
isrootdatum(W)=isweylgroup(W) || (W isa Spets && isrootdatum(Group(W)))
isspetsial=W->UnipotentCharacters(W)!==nothing

nspets=ComplexReflectionGroup.([5,7,9,10,11,12,13,15,16,17,18,19,20,22,21,22,31])

ChevieErr(x...)=xprintln("**! ",x...)

cox_ex=[coxgroup(:A,1), coxgroup(:A,2), coxgroup(:B,2), coxgroup(:G,2),
      coxgroup(:I,2,5), coxgroup(:A,3), coxgroup(:B,3), coxgroup(:C,3), 
      coxgroup(:H,3),  coxgroup(:A,4), coxgroup(:B,4), coxgroup(:C,4), 
  coxgroup(:D,4), coxgroup(:F,4), coxgroup(:H,4), coxgroup(:A,5), 
  coxgroup(:B,5), coxgroup(:C,5), coxgroup(:D,5), coxgroup(:A,6), 
  coxgroup(:B,6), coxgroup(:C,6), coxgroup(:E,6), coxgroup(:D,6), 
  coxgroup(:A,7), coxgroup(:B,7), coxgroup(:C,7), coxgroup(:E,7), 
  coxgroup(:D,7), coxgroup(:E,8), coxgroup()]

spets_ex=vcat(
  ComplexReflectionGroup.([4, 6, 8, 14, 24, 25, 26, 27, 29, 32, 33, 34]),
  [ComplexReflectionGroup(3,1,2),
  ComplexReflectionGroup(3,3,3),
  ComplexReflectionGroup(3,3,4),
  ComplexReflectionGroup(4,4,3)])

twisted=[rootdatum(:psu,3), rootdatum(Symbol("2B2")), rootdatum(Symbol("2G2")),
  rootdatum(Symbol("2I"),5), rootdatum(Symbol("2I"),8), rootdatum(:psu,4),
  rootdatum(Symbol("3D4")), rootdatum(:psu,5), rootdatum(Symbol("pso-"),8),
  rootdatum(:psu,6), rootdatum(Symbol("2F4")), rootdatum(:psu,7),
  rootdatum(Symbol("pso-"),10), rootdatum(Symbol("2E6")), rootdatum(:psu,8),
  rootdatum(Symbol("pso-"),12), rootdatum(Symbol("pso-"),14)]

all_ex=vcat(cox_ex,spets_ex,twisted)
sort!(all_ex,by=nconjugacy_classes)

function RG(s::Symbol)
  t=test[s]
  println("testing ",s,"\n",t.comment)
  for W in all_ex if t.applicable(W) 
  xprintln(s,"(",W,")")
@time  t.fn(W) 
  end end
end

# compares lists a and b (whose descriptions are strings na and nb)
function cmpvec(a,b;na="a",nb="b")
  if a==b return end
  if length(a)!=length(b)
    ChevieErr("length($na)=",length(a)," while length($nb)=",length(b))
  end
  if -a==b ChevieErr("$na=-$nn");return end
  pa=Perm(a,b)
  if pa!==nothing ChevieErr("$na=$nb^",pa);return end
  pa=Perm(a,-b)
  if pa!==nothing ChevieErr("$na=-$nb^",pa);return end
  for j in eachindex(a)
    if a[j]==b[j] continue end
    if a[j]==-b[j] ChevieErr("$na[$j]=-$nb[$j]");continue end
    t=findall(==(a[j]),b)
    if length(t)>0 ChevieErr("$na[$j] found at $t");continue end
    t=findall(==(-a[j]),b)
    if length(t)>0 ChevieErr("-$na[$j] found at $t");continue end
    ChevieErr("$na[$j] not found")
  end
end

# find matrices gr as a representation of a group or Hecke algebra
function findrepresentation(W,gr,check=false)
  O=W
  if O isa HeckeAlgebra W=O.W end
  t=classinfo(W)[:classtext]
  l=1:length(t)
  l=sort(l,by=i->length(t[i]))
  t=t[l]
  ct=CharTable(O).irr
  pos=1:length(l)
  gens=gr
  if isempty(gens) return 1 end
  for (i,w) in enumerate(t)
    r=traces_words_mats(gens,[w])[1]
#   println("w=$w r=$r")
#   if O isa HeckeAlgebra r=r*O.unit end
    pos=filter(j->iszero(ct[j,l[i]]-r),pos)
    if !check && length(pos)==1 return pos[1]
    elseif isempty(pos) return false
    end
  end
  if length(pos)==1 return pos[1]
  elseif length(pos)==0 return false
  else error("characters $pos are equal")
  end
end

function representations(W,l=Int[])
  O=W
  if W isa HeckeAlgebra
    H=W
    W=H.W
  else H=hecke(W)
  end
  cl=classinfo(W)[:classtext]
  ct=CharTable(O).irr
  if isempty(l) l=1:length(cl) end
  for i in l
    print("Representation #$i")
    gr=representation(O,i)
    if gr==false println("=false")
    else
      r=gr
      if !isempty(r) print(" dim:",(r isa Vector) ? size(r[1]) :
                          size(r.gens[1]),"...") end
      if isrepresentation(H,r) print("is ") end
      if r isa NamedTuple println();continue end # for now... Should be fixed
      pos=findrepresentation(O,gr,true)
      if pos==i println("found")
      elseif !isnothing(pos) println("character found at ",pos)
      else println("character does not match")
	pos=TransposedMat([ct[i],traces_words_mats(gr,cl)])
	f=List(pos,x->x[1]!=x[2])
	pos=List(pos,x->List(x,FormatGAP))
	Print(FormatTable(ListBlist(pos,f),
	  rec(rowLabels=ListBlist(1:length(pos),f),
	      columnLabels=["is","should be"])))
      end
    end
  end
end
test[:representations]=(fn=representations, 
   applicable=function(W)
     if nconjugacy_classes(W)>=55 return false end
     t=refltype(W)
     if isempty(t) return true end
     t=t[1]
     if !haskey(t,:orbit) || order(t.twist)!=2 return true end
     t=t.orbit[1]
     return !(t.series in [:D,:E])
   end,
   comment="Check reprs exist and match characters")

function lusztiginduction(WF)
  if !(WF isa Spets) WF=spets(WF) end
  W=Group(WF)
  for J in filter(x->length(x)<length(gens(W)),parabolic_representatives(W))
    lusztiginduction(WF,J)
  end
end
test[:lusztiginduction]=(fn=lusztiginduction, applicable=isspetsial,
   comment="Check induction computable and mackey with tori")

function lusztiginduction(WF,J::AbstractVector{<:Integer})
  for L in twistings(WF,J) lusztiginduction(WF,L) end
end

function lusztiginduction(WF,L)
  if L.phi==WF.phi print("Split ") end
  xprint("Lusztig Induction from ",L," to ",WF)
  t=LusztigInductionTable(L,WF)
  if isnothing(t) return end
  if haskey(t.prop,:scalars) 
    print("\t****** scalars=",t.prop[:scalars]) 
  end
  println()
  if L.phi==WF.phi
    h=Lusztig.HCInductionTable(L,WF)
    if h.scalar!=t.scalar error("HC!=Lusztig") end
  end
  uh=UnipotentCharacters(L)
  uw=UnipotentCharacters(WF)
  q=Pol([1],1)
  nh=charnames(uh)
  hd=degrees(uh,q)
  ud=degrees(uw,q)
  index=generic_order(WF,q)/generic_order(L,q)
  index*=q^-index.v
  index*=generic_sign(L)/generic_sign(WF)
  pred=hd.*index
  ind=map(x->UniChar(WF,x),eachcol(t.scalar))
  for j in 1:length(hd)
    if pred[j]!=degree(ind[j])
     xprint("!!  R_",L,"^",WF,"(",charnames(uh)[j],")=",ind[j],"\n!!",
        CycPol(degree(ind[j]))," instead of ",CycPol(pred[j]),
        " quotient ", CycPol(degree(ind[j]))*inv(CycPol(pred[j])),"\n")
    end
  end
  f=fusion_conjugacy_classes(L,WF)
#  Check Mackey with Tori
  c=map((a,b)->a//b,CharTable(WF).centralizers[f],CharTable(L).centralizers)
  u=Uch.DLCharTable(L)
  rhs=toM(map(1:HasType.NrConjugacyClasses(WF))do i
    l=findall(==(i),f)
    return isempty(l) ? u[1,:].*0 : sum(j->u[j,:].*c[j],l)
  end)
  rhs=Dict{Symbol,Any}(:scalar=>rhs, :u=>L,:g=>WF,
    :uNames=>UnipotentCharacters(L).prop[:TeXCharNames],
    :gNames=>classinfo(WF)[:classnames])
  m=copy(rhs)
  m[:scalar]=Uch.DLCharTable(WF)*t.scalar
  if m[:scalar]!=rhs[:scalar] error("tables differ",m[:scalar],rhs[:scalar]) end
  # Check transitivity with RTL
  if Uch.DLCharTable(L)*permutedims(t.scalar)!=Uch.DLCharTable(WF)[f,:]
     error("transitivity with RTL")
  end
end

using ..Gap4
function tCharTable(W)
# ct=improve_type(Gap4.CharTable(W).irr)
  ct=Gap4.CharTable(W).irr
  ct1=CharTable(W).irr
#  p=Perm(ct,ct1;dims=1) investigate why not enough
  p=Perm_rowcolmat(ct,ct1)
  if isnothing(p) error("irreducibles") 
  else println("permuted ",p)
  end
end

function position_class(W)
  cl=map(x->Gapjm.position_class(W,x),class_reps(W))
  if cl!=1:length(cl) error("classes") end
end

test[:position_class]=(fn=position_class, applicable=W->true,
   comment="classreps")

function unipotentclasses(W,p=nothing)
  function PosetFromICC(t)local l,o
    l=uc.springerseries[1][:locsys]
    o=map(i->map(j->any(a->any(b->!iszero(t.scalar[a,b]),
                        findall(x->x[1]==j,l)),findall(x->x[1]==i,l)),
        eachindex(uc.classes)),eachindex(uc.classes))
    Poset(o)
  end
  if W isa Spets WF=W;W=Group(WF)
  else WF=spets(W)
  end
  if isnothing(p)
    for p in vcat([0],badprimes(W)) unipotentclasses(WF,p) end
    return
  end
  uc=UnipotentClasses(WF,p)
  s=uc.springerseries[1]
  b=charinfo(W)[:b]
  du=[]
  pl(n0,u,i)=string("$n0.$i=",u.name,".",charnames(u.Au)[i])
  for (n0,u) in enumerate(uc.classes)
    nid=charinfo(u.Au)[:positionId]
    a=u.dimBu
    if haskey(u,:dynkin) && W.N!=0
      ad=W.N-div(count(x->!(x in [0,1]),toM(W.rootdec)*u.dynkin),2)
      if a!=ad ChevieErr(pl(n0,u,nid)," dimBu=>",a," and from Dynkin=>",ad)
      else a=ad
      end
    end
    push!(du,a)
    i=findfirst(y->y[1]==n0 && y[2]==nid,s[:locsys])
    bu=b[i]
    if a!=bu
      ChevieErr("p=$p ",pl(n0,u,nid)," dimBu=>$a and Springer[1.$i]=>$bu")
    end
    if a!=u.dimBu
      ChevieErr(pl(n0,u,nid)," dynkin=>",a," and dimBu=",u.dimBu,"\n")
    end
    for i in 1:nconjugacy_classes(u.Au) a=[]
      for j in 1:length(uc.springerseries)
        for k in 1:length(uc.springerseries[j][:locsys])
          if [n0,i]==uc.springerseries[j][:locsys][k]
	    push!(a,[j,k])
          end
        end
      end
      if isempty(a) ChevieErr(pl(n0,u,i)," not accounted for")
      elseif length(a)>1 ChevieErr(pl(n0,u,i)," at ",a)
      end
    end
  end
  for i in 1:length(uc.orderclasses)
    for a in hasse(uc.orderclasses)[i]
      if du[i]<=du[a]
        ChevieErr("bu($i:",uc.classes[i].name,")=",du[i]," <=bu($a:",
                  uc.classes[a].name,")=",du[a])
      end
    end
  end
  for s in uc.springerseries
    cl=first.(s[:locsys])
    order=incidence(uc.orderclasses)[cl,cl]
    t=ICCTable(uc,findfirst(==(s),uc.springerseries))
    for i in 1:length(cl), j in 1:length(cl)
      if !iszero(t.scalar[i,j]) && !order[i,j]
       ChevieErr(findfirst(==(s),uc.springerseries),"[$i,$j]!=0 and not ",
	uc.classes[cl[i]].name,"<",uc.classes[cl[j]].name)
      end
    end
    if s==uc.springerseries[1]
      if incidence(uc.orderclasses)!=incidence(PosetFromICC(t))
        ChevieErr("order bad")
      end
    end
  end
  bc=Ucl.BalaCarterLabels(W)
  if all(cl->haskey(cl,:dynkin),uc.classes)
    for cl in uc.classes
      j=findfirst(x->x[1]==cl.dynkin,bc)
      if haskey(cl,:balacarter)
        if cl.balacarter!=bc[j][2] error("balacarter") end
      else print("bala-carter[",cl.name,"] unbound, should be",bc[j][2])
      end
    end
  if all(x->x.series in [:E,:F,:G],refltype(W))
    for cl in uc.classes
      j=findfirst(x->x[1]==cl.dynkin,bc)
      name=Semisimple.IsomorphismType(reflection_subgroup(W,abs.(bc[j][2]));TeX=true)
      if name=="" name="1" end
      l=count(x->x<0,bc[j][2])
      if l!=0 
        if '+' in name name=replace(name,"+"=>string("(a_",l,")+"))
        else name*=string("(a_",l,")")
        end
      end
      name=replace(name,"+"=>"{+}")
      if name!=cl.name
       ChevieErr("name=<",fromTeX(rio(),cl.name),
                  "> but from Bala-Carter <",fromTeX(rio(),name),">")
      end
    end
  end
  end
  for u in uc.classes
    if haskey(u,:dimred) && haskey(u,:red) && u.dimred!=dimension(u.red)
       ChevieErr("for ",u.name," dimred=",u.dimred,"!= Dim(",
          Semisimple.IsomorphismType(u.red,torus=true),")=",dimension(u.red))
    end
  end
end

test[:unipotentclasses]=(fn=unipotentclasses, applicable=isrootdatum,
   comment="Check dimBu from DR, from b, with < ; Check BalaCarter, dimred")
end
