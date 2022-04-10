# defines RG()
module Tests
using Gapjm

const test=Dict{Symbol,NamedTuple{(:fn, :applicable,:comment),
        Tuple{Function,Function,String}}}()

isweylgroup(W)=(W isa FiniteCoxeterGroup) && all(isinteger,cartan(W))
isrootdatum(W)=isweylgroup(W) || (W isa Spets && isrootdatum(Group(W)))
isspetsial=W->UnipotentCharacters(W)!==nothing

nspets=ComplexReflectionGroup.([5,7,9,10,11,12,13,15,16,17,18,19,20,22,21,22,31])

ChevieErr(x...)=printstyled(rio(),x...;color=:red)

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

all_ex=vcat(cox_ex,spets_ex,nspets,twisted)
sort!(all_ex,by=nconjugacy_classes)

function RG(s::Symbol)
  t=test[s]
  println("testing ",s,"\n",t.comment)
  for W in all_ex if t.applicable(W) 
    printstyled(rio(),s,"(",W,")";bold=true,color=:magenta)
@time  t.fn(W) 
  end end
end

function RG(W)
  printstyled(rio(),"Tests for W=",W," -------------------------------\n";
            bold=true,color=:magenta)
  @time for (i,(s,t)) in enumerate(sort(collect(test),by=first))
    if t.applicable(W) 
      printstyled(i,"/",length(test),"  ",s,": ",t.comment,"\n";color=:green)
      t.fn(W) 
    end 
  end
end

RG(v::Vector)=for W in v RG(W) end
RG()=RG(all_ex)

# compares lists a and b (whose descriptions are strings na and nb)
function cmpvec(a,b;na="a",nb="b")
  if a==b return end
  if length(a)!=length(b)
    ChevieErr("length($na)=",length(a)," while length($nb)=",length(b),"\n")
  end
  if -a==b ChevieErr("$na=-$nn\n");return end
  pa=Perm(a,b)
  if pa!==nothing ChevieErr("$na=$nb^",pa,"\n");return end
  pa=Perm(a,-b)
  if pa!==nothing ChevieErr("$na=-$nb^",pa,"\n");return end
  for j in eachindex(a)
    if a[j]==b[j] continue end
    if a[j]==-b[j] ChevieErr(a[j],"=-",b[j],"\n");continue end
    t=findall(==(a[j]),b)
    if length(t)>0 ChevieErr(a[j]," found at ",t,"\n");continue end
    t=findall(==(-a[j]),b)
    if length(t)>0 ChevieErr("-",a[j]," found at ",t,"\n");continue end
    ChevieErr(na,"[$j]=",a[j]," not found\n")
  end
end

# cmptable(t1,t2[,nz]) nz: show all non-zero columns anyway
# compare two chevie tables or inductiontables
function cmptables(t1,t2,nz=false)
  t=[copy(t1),copy(t2)]
  opt=[]
  m=[]
  for i in [1,2]
    if t[i] isa InductionTable
      opt[i]=(rowLabels=t[i].gcharnames,columnLabels=t[i].ucharnames)
      m[i]=t[i].identifier
      t[i]=t[i].scalar
    elseif t[i] isa NamedTuple
      opt[i]=(rowLabels=t[i].rowLabels,columnLabels=t[i].columnLabels)
      m[i]=""
      t[i]=t[i].scalar
    else opt[i]=(rowLabels=axes(t[i],1),columnLabels=axes(t[i],2))
      m[i]=""
    end
  end
  if nz r=filter(i->t[1][i,:]!=t[2][i,:],axes(t[1],1))
  else r=filter(i->any(!iszero,t[1][i,:])|| any(!iszero,t[2][i,:]),axes(t[1],1))
  end
  msg=string("[",length(r),"/",length(t[1]),",")
  if isempty(r) InfoChevie("Tables agree!\n");return end
  p=Perm(t[1][r,:],t[2][r,:],dims=1)
  if p!==nothing
    ChevieErr("Permuted lines:\n",join(map(x->join(x,"->"),permutedims(
       [opt[1].rowLabels[r],Permuted(opt[1].rowLabels[r],p)])),"\n"),"\n")
    return
  end
  p=Perm(t[1],t[2],dims=2)
  if p!==nothing
    ChevieErr("Permuted columns:\n",
     join(map(c->join(opt[1].columnLabels[c],"->"),cycles(p))," "),"\n")
    return
  end
  if nz==2 c=filter(i->t[1][r,i]!=t[2][r,i],axes(t[1],2))
  else c=axes(t[1],2)
  end
  msg*=string(length(c),"/",size(t[1],2),"] of ")
  ChevieErr(msg,"\n",join(List([1,2],i->string(m[i],"\n",
     showtable(t[i];opt[i]...,rows=r,columns=c,screenColumns=70)))))
end

function cmpcycpol(a,b;na="a",nb="b")
  if a==b return end
  if b.coeff==0
    if a.coeff!=0 ChevieErr(na,"=",a," but ",nb,"=",b,"\n") end
    return
  end
  q=a*inv(b)
  if isempty(q.v) ChevieErr(na,"=",
    format_coefficient(repr(q;context=rio())),nb," where ",nb,"=",b,"\n");
  else ChevieErr(na,"=",a," but ",nb,"=",b,"\n")
  end
end
#---------------- test: representations ------------------------
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

function Trepresentations(W,l=Int[])
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
    InfoChevie("  #Representation #$i")
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

test[:representations]=(fn=Trepresentations, 
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

#---------------- test: lusztiginduction ------------------------
function Tlusztiginduction(WF)
  if !(WF isa Spets) WF=spets(WF) end
  W=Group(WF)
  for J in filter(x->length(x)<length(gens(W)),parabolic_reps(WF))
    Tlusztiginduction(WF,J)
  end
end
test[:lusztiginduction]=(fn=Tlusztiginduction, applicable=isspetsial,
   comment="Check induction computable and mackey with tori")

function Tlusztiginduction(WF,J::AbstractVector{<:Integer})
  for L in twistings(WF,J) Tlusztiginduction(WF,L) end
end

function Tlusztiginduction(WF,L)
  InfoChevie("  #")
  if L.phi==WF.phi InfoChevie("Split ") end
  InfoChevie("Lusztig Induction from ",L," to ",WF)
  if !isnothing(match(r"3G3,3,3",isomorphism_type(L))) return end
  t=LusztigInductionTable(L,WF)
  if isnothing(t) return end
  if haskey(t,:scalars) 
    print("\t****** scalars=",t.scalars) 
  end
  println()
  if L.phi==WF.phi
    h=Lusztig.HCInductionTable(L,WF)
    if h.scalar!=t.scalar ChevieErr("HC!=Lusztig\n") end
  end
  uh=UnipotentCharacters(L)
  uw=UnipotentCharacters(WF)
  q=Pol([1],1)
  nh=charnames(uh)
  hd=degrees(uh,q)
  ud=degrees(uw,q)
  index=exactdiv(generic_order(WF,q),generic_order(L,q))
  index*=q^-index.v
  index*=generic_sign(L)//generic_sign(WF)
  pred=hd.*index
  ind=map(x->UniChar(WF,x),eachcol(t.scalar))
  for j in 1:length(hd)
    if pred[j]!=degree(ind[j])
     xprintln("!!  R_",L,"^",WF,"(",charnames(uh)[j],")=",ind[j])
     xprintln("!!",CycPol(degree(ind[j]))," instead of ",CycPol(pred[j]),
        " quotient ", CycPol(degree(ind[j]))*inv(CycPol(pred[j])))
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
                       :uNames=>charnames(UnipotentCharacters(L);TeX=true),
    :gNames=>classinfo(WF)[:classnames])
  m=copy(rhs)
  m[:scalar]=Uch.DLCharTable(WF)*t.scalar
  if m[:scalar]!=rhs[:scalar] error("tables differ",m[:scalar],rhs[:scalar]) end
  # Check transitivity with RTL
  if Uch.DLCharTable(L)*permutedims(t.scalar)!=Uch.DLCharTable(WF)[f,:]
     error("transitivity with RTL")
  end
end

#----------------test: chartable ------------------------
using GAP
function Tchartable(W)
  G=Group(gens(W))
  G.classreps=classreps(W)
  ct=CharTable(G).irr
  ct1=CharTable(W).irr
  p=Perm(ct,ct1;dims=1)
  if isnothing(p) error("irreducibles") 
  elseif !isone(p) InfoChevie("  #permuted ",p,"\n")
  end
end

test[:chartable]=(fn=Tchartable, applicable=W->!(W isa Spets),
   comment="CharTable")

#----------------test: powermaps ------------------------

function Tpowermaps(W)
  cl=classreps(W)
  p=classinfo(W)[:powermap]
  for i in echindex(p)
    if isdefine(p,i) cmap=map(x->position_class(W,x^i),cl)
      if cmap!=p[i]
        ChevieErr(i,"-th power map is ",p[i],"\nshould be:",cmap,"\n")
      end
    end
  end
end

test[:powermaps]=(fn=Tpowermaps, applicable=W->false,
#test[:powermaps]=(fn=Tpowermaps, applicable=W->!(W isa Spets),
                  comment="powermaps")

#---------------- test: positionclasses ------------------------
function Tpositionclasses(W)
  cl=map(x->Gapjm.position_class(W,x),classreps(W))
  if cl!=1:length(cl) error("classes") end
end

test[:positionclasses]=(fn=Tpositionclasses, applicable=W->true,
   comment="classreps")

#---------------- test: unipotentclasses ------------------------
function Tunipotentclasses(W,p=nothing)
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
    for p in vcat([0],badprimes(W)) Tunipotentclasses(WF,p) end
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
      ChevieErr("p=$p ",pl(n0,u,nid)," dimBu=>$a and Springer[1.$i]=>$bu\n")
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
      name=isomorphism_type(reflection_subgroup(W,abs.(bc[j][2]));TeX=true)
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
                  "> but from Bala-Carter <",fromTeX(rio(),name),">\n")
      end
    end
  end
  end
  for u in uc.classes
    if haskey(u,:dimred) && haskey(u,:red) && u.dimred!=dimension(u.red)
       ChevieErr("for ",u.name," dimred=",u.dimred,"!= Dim(",
          isomorphism_type(u.red,torus=true),")=",dimension(u.red))
    end
  end
end

test[:unipotentclasses]=(fn=Tunipotentclasses, applicable=isrootdatum,
   comment="dimBu from DR, from b, with < ; BalaCarter, dimred")

#---------------- test: unipotentcentralizers ------------------------
function Tunipotentcentralizers(W,p=0)
  t=RationalUnipotentClasses(W,p)
  q=Pol()
  sum=0
  u=UnipotentClasses(W,p)
  g=CycPol(generic_order(W,q))
  for c in t
    cl=u.classes[c[:classno]]
    cc=q^cl.dimunip
    if nconjugacy_classes(cl.Au)>1
      if rank(cl.red)>0 && haskey(cl,:AuAction)
        ww=classinfo(cl.Au)[:classtext][c[:AuNo]]
        w=cl.red.F
        if !isempty(ww) w*=prod(cl.AuAction.F0s[ww]) end
        cc*=generic_order(spets(Group(cl.red),w),q)
      else cc*=generic_order(cl.red,q)
      end
      cc*=div(length(cl.Au),classinfo(cl.Au)[:classes][c[:AuNo]])
    else cc*=generic_order(cl.red,q);
    end
    sum+=c[:card](q)
    cc=CycPol(cc)
    if cc!=g//c[:card]
     ChevieErr(fromTeX(rio(),cl.name),".",c[:AuNo],":",
        Ucl.showcentralizer(rio(),cl)," => ",cc," (is ",g//c[:card],")\n");
    end
  end
  if W isa Spets N=Group(W).N else N=W.N end
  if sum!=q^(2*N)
    ChevieErr("found nr. unip=",sum," instead of ",q^(2*N),"\n")
  end
end

test[:unipotentcentralizers]=(fn=Tunipotentcentralizers,applicable=isrootdatum,
                         comment="info on centralizers agrees with ICCTable")

#---------------- test: nrssclasses ------------------------
# Check classtypes gives nrclasses (for simply connected groups)
function Tnsemisimple(WF)
  q=Mvp(:q)
  g=generic_order(WF,q)
  l=ClassTypes(WF)()
  v=map(l.ss) do x
    x.nClasses*exactdiv(g,x.cent(q))*q^(2*Group(x.CGs).N)
  end
  if sum(v)!=g
    ChevieErr("found nr elem=",sum(v)," instead of ",g,"\n");
   end
end

test[:nsemisimple]=(fn=Tnsemisimple,applicable=W->isweylgroup(W) && 
                    length(fundamental_group(W))==1 && length(W)<10^8,
           comment="classtypes gives nr. semisimple classes")
#---------------- test: charparams ------------------------
"""
Jean Michel and Ulrich Thiel (2017)

this function checks that the parameters/names for characters of Complex
eflection groups agree with the description which is since 2016 in the CHEVIE 
manual. If not, it tells what permutation of the data is needed.
"""
function Tcharparams(W)
  ct=CharTable(W).irr
  fd=fakedegrees(W,Pol())
  db=map(x->[x(1),valuation(x)],fd)
  n=repr(W;context=rio())
  if haskey(charinfo(W),:malle) l=charinfo(W)[:malle]
    nm=map(x->CHEVIE[:G2][:CharName](x,Dict()),l)
    nm=map(x->replace(x,r"[{}]"=>""),nm)
  elseif n[1] in "EFGH" && length(n)<4
    l=filter(x->x[2]>1,tally(charinfo(W)[:charparams]))
    if !isempty(l)>0 ChevieErr("Duplicated:",l) end
    l=first.(charinfo(W)[:charparams])
    nm=charnames(W)
  elseif n[1] in ".ABCD" || length(n)>=4 return
  else
    l=first.(charinfo(W)[:charparams])
    nm=charnames(W)
  end
  ddb=filter(x->length(x)>1,map(x->findall(==(x),db),unique(db)))
  n0=x->findfirst(==(x),nm)
  nbbis(n)=count(==('\''),n)+2*count(==('"'),n)
  # lexicographic polynomial order starting with highest coeff
  function lexpol(p,q)local a,b
    a=vcat(fill(0,p.v),p.c)
    b=vcat(fill(0,q.v),q.c)
    if length(a)<length(b) append!(a,fill(0,length(b)-length(a))) end
    if length(b)<length(a) append!(b,fill(0,length(a)-length(b))) end
    return reverse(a)<=reverse(b)
  end
  # check(n,f)
  # chars of name n+qualifiers satisfy f
  function check(n,fn)local qq,a,i,j
    n=findall(x->length(x)>=length(n) && x[1:length(n)]==n,nm)
    sort!(n,by=i->nbbis(nm[i]))
    aa=arrangements(1:length(n),length(n))
    a=aa[findfirst(x->fn(nm[n[x]]),aa)]
    perm*=mappingPerm(n,n[a])
    for i in n 
      j=findfirst(x->i in x,ddb)
      if j===nothing if i==n[1] print("useless test: ",nm[n],"\n") end
      else ddb[j]=setdiff(ddb[j],[i])
        if length(ddb[j])<=1 deleteat!(ddb,j) end
      end
    end
  end
  function ok(cond::Bool,arg...)
    if !cond ChevieErr(arg...) end
    cond
  end
  DecomposeTensor(c...)=decompose(CharTable(W),
                                  vec(prod(view(ct,collect(c),:),dims=1)))
  function intensor(chi,arg...)# arg[1] occurs in arg[2]⊗ …⊗ arg[n]
    ok(0!=DecomposeTensor(n0.(arg)...)[n0(chi)],chi," should occur in ⊗",arg)
  end
  function fakdeghighterms(a,p) 
    ok(degree(fd[n0(a)]-p)<valuation(p),
    "FakeDegree(",a,") should have high terms ",p)
  end
  function fakdegsmaller(a,b) 
    ok(degree(fd[n0(a)])<degree(fd[n0(b)]),
    "FakeDeg(",a,") should have a smaller degree than FakeDeg(",b,")")
  end
  function issign(a,b)local det
    det=map(x->position_class(W,x),gens(W))
    det=findfirst(l->all(y->y==-1,l),collect(eachrow(ct[:,det])))
    ok(1==DecomposeTensor(det,n0(b))[n0(a)],a," should be ⊗ sign of ",b)
  end
  isconj(a,b)=ok(ct[n0(a)]==conj(ct[n0(b)]),a," should be conj(",b,")")
  isreal(a)=ok(ct[n0(a)]==conj(ct[n0(a)]),a," should be real")
  function isgal(a,b,g) 
    ok(ct[n0(a),:]==galois.(ct[n0(b),:],g),a," should be galois(",g,") of ",b)
  end
  function inducedfromid(a,J;mult=0)local L,t
    L=reflection_subgroup(W,J)
    t=InductionTable(L,W)
    mul(x)=mult!=0 ? x==mult : x>0
    ok(mul(t.scalar[n0(a),charinfo(L)[:positionId]]),
     a," should occur",mult!=0 ? " $mult times" : "",
     " in the induced of Id from ", repr(L;context=rio()))
  end
  function value(a,c,v)
    c=position_class(W,W(c))
    ok(ct[n0(a),c]==v,a," should take value ",v," at class ",c)
  end
  # lexicographically minimum tensor of (d,b)-distinguishable chars
  # in which char no. i occurs
  function mintensor(W,i)local bydb,b,a,ch,gb,m
    bydb=1+findall(i->!('\'' in i),nm[2:end])
    sort!(bydb,by=i->[ct[i,1],charinfo(W)[:b][i]])
    for a in bydb
      ch=conj(ct[a,:]).*ct[i,:]
      gb=filter(j->ct[j,1]<=ct[i,1]*ct[a,1],bydb)
      m=map(x->Chars.scalarproduct(CharTable(W),x,ch),eachrow(ct[gb,:]))
      b=findfirst(!iszero,m)
      if b!==nothing return [a,gb[b],m[b]] end
    end
  end
  function rules()local b,p,q,bad,pp,g,u,bb,o,rule
    if !all(x->x isa Integer,cartan(W))
      g=conductor(cartan(W))
      g=map(x->Perm(ct,galois.(ct,x),dims=1),prime_residues(g))
      g=orbits(Group(g),1:length(nm))
      o=i->g[findfirst(o->i in o,g)]
    end
    rule=[]
    bad=map(ddb)do p
      function fin(rno::Int,cond::Bool,arg...)
        if !cond push!(rule,-rno)
          ChevieErr("Rule ",rno," failed:",arg...)
          Perm(p[1],p[2])
        else push!(rule,rno)
          Perm()
        end
      end
      if length(p)>2 return true end
      q=i->nbbis(nm[i])
      if fd[p[1]]!=fd[p[2]]
        return fin(1,lexpol(fd[p[1]],fd[p[2]])==(q(p[1])<q(p[2])),
	  nm[p[1]],"->",fd[p[1]]," and ",nm[p[2]],"->",fd[p[2]],"\n")
      end
      if @isdefined o
        pp=map(i->filter(j->!(j in vcat(ddb...)),o(i)),p)
	if count(!isempty,pp)==1
	  u=findfirst(!isempty,pp)
          return fin(2,q(p[u])==1,nm[p[u]],"~",join(nm[pp[u]],","),"\n")
        end
	pp=o.(p)
        bb=map(sort,map(i->map(y->db[y][2],i),pp))
	if bb[1]!=bb[2]
          return fin(3,(bb[1]<bb[2])==(q(p[1])<q(p[2])),nm[p[1]],"~",
            join(nm[pp[1]],",")," and ",nm[p[2]],"~",join(nm[pp[2]],","),"\n")
          end
        end
      bb=improve_type(map(i->mintensor(W,i),p))
      pp=map(a->[db[a[1]][1],db[a[1]][2],db[a[2]][1],db[a[2]][2],a[3]],bb)
      if pp[1]!=pp[2]
       return fin(4,(pp[1]<pp[2])==(q(p[1])<q(p[2])),
             nm[p[1]],"|",join(nm[bb[1][1:2]],"o "),"=",pp[1][5]," & ", 
             nm[p[2]],"|",join(nm[bb[2][1:2]],"o "),"=",pp[2][5],"\n")
      end
      return true
    end
    if !isempty(rule) InfoChevie("#I applied  rules ",
             join(map(x->join(x,"x"),tally(rule)),","),"\n");end
    perm=prod(filter(x->x!=true,bad))
    if perm==1 perm=Perm() end
    ddb=ddb[map(==(true),bad)]
  end
  if @isdefined l 
    if db!=map(x->x[1:2],l)
    ChevieErr("disagree with (χ(1),b_χ)\n")
    PrintArray(Filtered(TransposedMat([db,List(l,x->x[1:2])]),x->x[1]!=x[2]))
    Error()
    end
  end
  perm=Perm()
  if n=="G₃‚₃‚₄" rules()
    check("phi6,5",v->intensor(v[1],"phi4,1","phi4,1"))
  elseif n=="G₃‚₃‚₅" rules()
    check("phi30,7",v->intensor(v[2],"phi5,1","phi5,1","phi5,1","phi5,1"))
  elseif n=="G₄‚₄‚₃" rules()
    check("phi3,2",v->isconj(v[2],"phi3,1"))
    check("phi3,6",v->isconj(v[2],"phi3,5"))
  elseif n=="G₂" rules()
    check("phi1,3",v->value(v[1],1,-1)) # W.rootLengths=[3,1]
  elseif n=="G₅" rules()
    check("phi1,4'",v->value(v[1],1,1))
    check("phi1,8",v->
    isconj(v[2],"phi1,16") && isconj(v[1],"phi1,4'") && isconj(v[3],"phi1,4''"))
    check("phi1,12",v->value(v[1],1,E(3))&& isconj(v[1],v[2]))
    check("phi2,5",v->isconj(v[2],"phi2,1")&& value(v[1],1,-1))
    check("phi2,7",v->isconj(v[1],"phi2,5'") 
     && isconj(v[2],"phi2,5'''"))
    check("phi2,3",v->value(v[1],1,-E(3))&& isconj(v[2],v[1]))
  elseif n=="G₇" rules()
    check("phi1,4",v->value(v[1],2,1))
    check("phi1,8",v->isconj(v[2],"phi1,16") && isconj(v[1],"phi1,4'"))
    check("phi1,10",v->value(v[1],2,1))
    check("phi1,12",v->value(v[1],2,E(3)))
    check("phi1,14",v->isconj(v[2],"phi1,22") && isconj(v[1],"phi1,10'"))
    check("phi1,18",v->value(v[1],2,E(3)))
    check("phi2,3",v->value(v[1],2,-E(3)))
    check("phi2,5",v->isgal(v[2],"phi2,1",5) && isgal(v[1],"phi2,13'",5))
    check("phi2,7",v->isgal(v[3],"phi2,1",7) && isgal(v[1],"phi2,13'",7))
    check("phi2,9",v->isconj(v[1],"phi2,15") && isconj(v[3],"phi2,3'"))
    check("phi2,11",v->isconj(v[2],"phi2,1") && isconj(v[1],"phi2,13'"))
    check("phi2,13",v->value(v[1],2,-1))
  elseif n=="G₂₇" 
    check("phi3,5",v->fakdegsmaller(v[1],v[2]))
    check("phi3,20",v->fakdegsmaller(v[1],v[2]))
    check("phi8,9",v->fakdegsmaller(v[2],v[1]))
    check("phi5,15",v->inducedfromid(v[1],[1,3]))
    check("phi5,6",v->inducedfromid(v[1],[1,3];mult=1))
  elseif n=="F₄" # W.rootLengths=[3,1]
    check("phi1,12",v->inducedfromid(v[1],[3,4]))
    check("phi2,4",v->inducedfromid(v[1],[3,4]))
    check("phi2,16",v->issign(v[1],"phi2,4''"))
    check("phi4,7",v->inducedfromid(v[1],[3,4]))
    check("phi6,6",v->intensor(v[2],"phi4,1","phi4,1"))
    check("phi8,9",v->inducedfromid(v[1],[3,4]))
    check("phi8,3",v->issign(v[1],"phi8,9''"))
    check("phi9,6",v->inducedfromid(v[1],[3,4]))
  elseif n=="G₂₉" 
    check("phi15,4",v->intensor(v[2],"phi4,1","phi4,3"))
    check("phi15,12",v->issign(v[2],"phi15,4''"))
    check("phi6,10",v->intensor(v[3],"phi4,1","phi4,1") 
     && isconj(v[4],v[3]) && inducedfromid(v[1],[1,2,4]))
  elseif n=="H₄" 
    check("phi30,10",v->intensor(v[1],"phi9,2","phi6,12"))
  elseif n=="G₃₁" 
    check("phi15,8",v->intensor(v[1],"phi4,1","phi20,7"))
    check("phi15,20",v->intensor(v[1],"phi4,1","phi20,7"))
    check("phi20,13",v->isconj(v[1],"phi20,7"))
    check("phi30,10",v->fakdeghighterms(v[1],Pol()^50+Pol()^46))
    check("phi45,8",v->intensor(v[2],"phi4,1","phi20,7"))
    check("phi45,12",v->issign(v[1],"phi45,8'"))
  elseif n=="G₃" 
    check("phi10,8",v->fakdeghighterms(v[1],Pol()^28+Pol()^26))
    check("phi10,17",v->issign(v[1],"phi10,8'"))
    check("phi40,5",v->fakdeghighterms(v[1],Pol()^31+Pol()^29+2*Pol()^27))
    check("phi40,14",v->issign(v[1],"phi40,5'"))
  elseif n=="G₃₄" 
    check("phi20,33",v->intensor(v[1],"phi6,1","phi15,14"))
    check("phi70,9",v->intensor(v[2],"phi6,1","phi15,14") && isreal(v[1]))
    check("phi70,45",v->issign(v[2],"phi70,9''") && isreal(v[1]))
    check("phi105,8",v->isconj(v[1],"phi105,4"))
    check("phi120,21",v->inducedfromid(v[1],[1,2,4,5,6]))
    check("phi280,12",v->intensor(v[1],"phi6,1","phi336,17"))
    check("phi280,30",v->intensor(v[2],"phi6,1","phi336,17"))
    check("phi540,21",v->intensor(v[1],"phi6,1","phi105,20"))
    check("phi560,18",v->intensor(v[3],"phi6,1","phi336,17") && isreal(v[1]))
    check("phi840,13",v->isconj(v[1],"phi840,11"))
    check("phi840,23",v->isconj(v[1],"phi840,19"))
  elseif n in ["G₆","G₈","G₉","G₁₀","G₁₁","G₁₃","G₁₄","G₁₅","G₁₆","G₁₇",
    "G₁₈","G₁₉","G₂₀","G₂₁","G₂₅","G₂₆","G₃₂"] 
    rules()
  else println("not applicable: ",n); return
  end
  if !isempty(ddb) println("not separated:",ddb) end
  if !isone(perm)
    println("Permutation to do:")
    perm=cycles(perm)
    if all(x->length(x)==2,perm) 
      for t in perm println(nm[t[1]]," <=> ",nm[t[2]]) end
    else println(perm)
    end
  end
end

test[:charparams]=(fn=Tcharparams, applicable=W->W isa Group,
   comment="Check charparams for consistency with Michel/Thiel rules")

#---------------- test: HCdegrees ------------------------
function THCdegrees(W)
  uw=UnipotentCharacters(W)
  for i in eachindex(uw.harishChandra)
   if !isempty(uw.harishChandra[i][:relativeType]) THCdegrees(W,i) end
  end
end

# HCdegrees(W [,ser [,guess]])
# Checks that unipotent degrees from Harish-Chandra series no. i of the Spets
# W agree with degrees coming from Schur elements of relative Hecke algebra
function THCdegrees(W,i,rel=false)
  uw=UnipotentCharacters(W)
  hw=uw.harishChandra[i]
  L=reflection_subgroup(W,hw[:levi])
  index=CycPol(generic_order(W,Pol()))//CycPol(generic_order(L,Pol()))
  index*=CycPol(1,-valuation(index))
  n=hw[:cuspidalName]
  if n=="" n="." end
  cusp=CycPolUnipotentDegrees(L)[Lusztig.FindCuspidalInLevi(n,L)]
  InfoChevie("  #HC_",L,"(cusp=",fromTeX(rio(),n),":",cusp,")[G:L]=",index,"\n")
  R=reflection_group(hw[:relativeType]...)
# check parameters of relative algebra by the formula
# u_{s,j}=ζ_s^j  q^(([a+A]_ρ_{det_s^j}-[a+A]_ρ_{\Id})/order#class(s)
# where ρ_χ is the unip. char corresp. to χ in relative group
  aA=map(x->uw.a[x]+uw.A[x],hw[:charNumbers])
  para=map(hyperplane_orbits(R))do h
    l=map(i->(aA[i]-aA[charinfo(R)[:positionId]])//h.order//h.N_s,h.det_s)
    max(maximum(l),0).-vcat([0],l)
  end
  ss=simple_reps(R)
  o=unique(ss)
  para=para[map(i->findfirst(==(ss[i]),o),eachindex(gens(R)))]
  para=map(l->all(iszero,l[2:end]) ? l[1] : l,para)
  if para!=hw[:parameterExponents]
    ChevieErr("Parameter: computed:",para,
              " stored:",hw[:parameterExponents]," for ",R)
  end
  if isempty(hw[:parameterExponents]) den=1
  else den=lcm(denominator.(vcat(hw[:parameterExponents]...)))
  end
  ud=CycPolUnipotentDegrees(W)
  reldeg=map(x->subs(cusp*index//x,Pol()^den),ud[hw[:charNumbers]])
  H=Uch.relative_hecke(uw,i,Mvp(:x)^den)
  InfoChevie("   Relative ",R,"\n");
  sch=schur_elements(H)
  if any(==(false),sch)
    ChevieErr("Schur elements not implemented for ",H,"!")
    return false
  end
  sch=CycPol.(sch)
  cmpvec(sch,reldeg;na=string("Schur(",repr(R;context=rio()),")"),nb="ud")
  if rel reldeg
  else  Ref(subs(cusp*index,Pol()^den)).//(sch*1//1)
  end
end

test[:HCdegrees]=(fn=THCdegrees, applicable=isspetsial,
  comment="gendegs from relative Hecke algebras of 1-HC cuspidals")

# EigenAndDegHecke(s) 
# s is a series. Uses only s.Hecke, s.d, 
# (and s.degree, s.levi and s.cuspidal for non-principal series)
# 
# Assume H=s.Hecke is a spetsial exp(2i pi d)-cyclotomic algebra, 
# for d in Q/Z describing a central element of Group(H).
# Returns a list of records for each character chi of H with fields
#   deg:=Generic degree = S_Id/S_chi
#   eig:=rootUnity part of Eigenvalue of Frobenius
#   frac:= fractional power in q-part of Eigenvalue of Frobenius
function EigenAndDegHecke(s)
  FractionToRoot(x)=E(denominator(x),numerator(x))
  H=hecke(s)
  d=s.d
  zeta=Cyc(d)
  W=H.W
  ct=CharTable(H).irr
  ct1=0 .+CharTable(W).irr
  ct2=improve_type(scalar(value.(ct,Ref(:q=>zeta))))
  n=axes(ct2,1)
  good=filter(i->!any(x->x isa HasType.Unknown,ct2[:,i]),n)
  p=Perm(ct1[:,good],ct2[:,good],dims=1) #Permuted(ct,p) specializes
  if !isone(p) println("***** perm=",p) 
    if iscyclic(W) ChevieErr("should not have perm")end 
  end
  d1=exponent(d)//gcd(order(d),gcd(degrees((W))))
  i=position_regular_class(W,d1) # this class represents F
  omegachi=map(x->scalar(value(x[i]//x[1],:q=>1)),eachrow(^(ct,p,dims=1)))
  frac=degree.(central_monomials(H)).*d1
  om2=map((o,p)->o*d^-p,map(x->x[i]//x[1],eachrow(ct1)),frac)
  if omegachi!=om2 
    omegachi=om2
    ChevieErr(classinfo(W)[:classtext][i],"^",1//d1," not equal to π(",W,")\n")
  end
  ss=CycPol.(schur_elements(H)^p)
  ss=map(x->degree(s)//x,ss)
# omegachi*=Eigenvalues(UnipotentCharacters(s.levi))[s.cuspidal]
  zeta=Cyc(Root1(;r=d1)^frac[charinfo(W)[:positionId]])
  if haskey(s,:delta) && s.delta!=1 && iscyclic(W)
    omegachi=map(i->Root1(;r=s.d*s.e*s.delta*
              ((i-1)//s.e-dSeries.mC(s)[i]*s.d)),1:s.e)
    zeta=Root1(;r=s.d*s.e*s.delta*s.d*dSeries.mC(s)[1])
  end
  map((deg,eig,frac)->(deg=deg,eig=eig*zeta,frac=modZ(frac)),ss,omegachi,frac)
end

#---------------- test: series ------------------------
function TSerie(s)
  W=s.spets
  InfoChevie("  # ",s,"\n")
  relative_group(s)
  if s.principal
    Ed=s.d
    c=(generic_order(s.spets,Pol())//sum(x->x^2,s.WGLdims)//
       generic_order(s.levi,Pol()))(Ed)
    if c!=1 ChevieErr(s," => |G|/(|L||WGL|) mod Φ=",c,"\n") end
    e=generic_sign(s.spets)*Ed^(sum(codegrees(Group(s.spets)))-
       sum(codegrees(Group(s.levi))))//generic_sign(s.levi)
    if e!=1 ChevieErr(s," => εL/c mod Φ=",e,"\n")end
  end
  if isnothing(dSeries.char_numbers(s)) return end
  e=eigen(UnipotentCharacters(W))[dSeries.char_numbers(s)]
  if haskey(s,:hecke)
    pred=map(x->x.eig,EigenAndDegHecke(s))*
      eigen(UnipotentCharacters(s.levi))[s.cuspidal]
    if e!=pred && (!haskey(s,:delta) || s.delta==1 || is_cyclic(s.WGL))
     n=charnames(UnipotentCharacters(W))[s.charNumbers]
      ChevieErr(s," actual eigen differ from predicted eigen")
      c=Set(map((x,y)->x//y,pred,e))
      if length(c)==1 ChevieErr("predicted=",c[1],"*actual")
      else cmptables(
                 (rowLabels=["actual   "],columnLabels=n,scalar=[e]),
                 (rowLabels=["predicted"],columnLabels=n,scalar=[pred]))
      end
    end
  end
end

# checkSeries(WF[,d[,ad]])
function Tseries(W,arg...)
  l=Series(W,arg...;proper=true)
  for s in l TSerie(s) end
  l
end

test[:series]=(fn=Tseries,applicable=isspetsial,comment="check d-HC series")

#---------------- test: extrefl ------------------------
using LinearAlgebra: tr
# test ReflectionEigenvalues, ChevieCharInfo.extRefl, .positionDet, .positionId
function Textrefl(W)
  ct=CharTable(W)
  # compute first using ReflectionEigenvalues
  n=nconjugacy_classes(W)
  v=reverse(permutedims(toM(map(r->prod(n->Pol()+n,r;init=Pol(1)).c,refleigen(W))));dims=1)
  # check v[2,:] using reflrep
  if size(v,1)>1 && v[2,:]!=map(w->tr(reflrep(W,W(w...))),classinfo(W)[:classtext])
   ChevieErr("refleigen disagrees with reflrep")
  end
  extRefl=map(x->findfirst(==(x),collect(eachrow(ct.irr))),eachrow(v))
  ci=charinfo(W)
  function checkfield(f,v)
    if !haskey(ci,f)
      ChevieErr(".",f," not bound but should be ",v,"\n");
    elseif v!=ci[f]
      ChevieErr(".",f," is ",ci[f]," should be ",v,"\n");
    end
  end
  if haskey(ci,:extrefl) checkfield(:extRefl,extRefl) end
  checkfield(:positionDet,extRefl[end])
  checkfield(:positionId,extRefl[1])
end

test[:extrefl]=(fn=Textrefl,applicable=x->x isa Group,comment="check extrefl")

#---------------- test: degrees ------------------------
"""
The function below uses the formula
``∏_{g∈ W}det(g-t)=∏_i(t^{d_i}-1)^{|W|/d_i}``
Superseded by degrees using refltype but can recompute data
"""
function reflectiondegrees(W)
  c=classinfo(W)[:classes]
  if isempty(c) return Int[] end
  l=vcat(map((r,c)->map(x->(x,c),r),refleigen(W),c)...)
  e=collectby(first,l)
  eig=map(x->first(first(x)),e)
  mul=map(x->sum(last,x),e)
  res=Int[]
  for d in sort(unique(order.(eig)),rev=true)
    p=findfirst(==(E(d)),eig)
    if p!==nothing && mul[p]>0
      for i in 1:div(d*mul[p],length(W)) push!(res,d) end
      n=mul[p]
      for j in 0:d-1 mul[findfirst(==(E(d,j)),eig)]-=n end
    end
  end
  reverse(res)
end

"""
the function below uses the formula cite[page 164]{BLM06}
``∏_{g∈ W}det(t-gϕ)=∏_i(t^{d_i}-εᵢ)^{|W|/d_i}``
Superseded by degrees using refltype but can recompute data
"""
function reflectiondegrees(W::Spets)
  c=classinfo(W)[:classes].//length(W)
  l=vcat(map((v,t)->map(x->(x,t),v),refleigen(W),c)...)
  l=collectby(first,l)
  e=map(x->(x[1][1],sum(last,x)),l)
  # here we got LHS of formula as ∏ᵢ(t-eig[i])^mul[i]
  # more exactly a list of pairs (eig::Root1, mul/|W|::Rational)
  searchdeg=function(e,card,degs)local res,d,f,g,p,pos,ne
#   @show degs,card
    if isempty(degs) return [Pair{Int,Root1}[]] end
    e=filter(x->x[2]!=0,e)
    res=Vector{Pair{Int,Root1}}[]
    d=degs[1]
    f=findall(i->last(i)>=1//d,e)
    g=first.(filter(x->first(x).r<1//d,e[f]))
    for p in g
      pos=map(i->findfirst(j->first(e[j])==p*i,f),E.(d,0:d-1))
      if !(nothing in pos)
#       @show d,p
        ne=copy(e)
        ne[f[pos]]=map(x->(x[1],x[2]-1//d),ne[f[pos]])
        for ne in searchdeg(ne,div(card,d),degs[2:end])
          push!(res,pushfirst!(ne,d=>p^d))
        end
      end
    end
    res
  end
  mul=searchdeg(e,length(W),sort(degrees(Group(W)),rev=true))
  e=unique(sort(sort.(mul)))
#  @show mul
  map(x->(x[1],Cyc(x[2])),only(e))
end

function Tdegrees(W)
  d=sort(degrees(W))
  d1=sort(reflectiondegrees(W))
  cmpvec(d,d1;na="degrees",nb="reflectiondegrees")
  if isweylgroup(W) && length(refltype(W))==1 && rank(W)==semisimplerank(W)
    d1=last.(tally(sum.(W.rootdec[1:nref(W)])))
    d1=sort(1+conjugate_partition(d1))
    cmpvec(d,d1;na="degrees",nb="dual partition of root heights")
  end
end

test[:degrees]=(fn=Tdegrees,applicable=x->true, comment="check degrees")

#---------------- test: fakedegrees ------------------------

function Tfeg(W)
  q=Pol()
  fd=fakedegrees(W,q)
  ffd=fakedegrees(W,q;recompute=true)
  cmpvec(ffd,fd;na="recomputed",nb="fakedegrees")
  cmpvec(valuation.(fd),charinfo(W).b;na="computed b",nb="charinfo b")
  cmpvec(degree.(fd),charinfo(W).B;na="computed B",nb="charinfo B")
end

test[:feg]=(fn=Tfeg,applicable=x->true, comment="check fakedegrees")

#---------------- test: fakedegrees induce ------------------------

Tfeginduce(W)=for J in parabolic_reps(W) Tfeginduce(W,J) end

function Tfeginduce(W,J)
  q=Pol()
  if !(W isa Spets) W=spets(W) end
  for L in twistings(W,J)
    t=InductionTable(L,W)
    hd=fakedegrees(L,q)
    ud=fakedegrees(W,q)
    index=exactdiv(generic_sign(L)*generic_order(W,q),generic_order(L,q))//generic_sign(W)
    index=shift(index,-valuation(index))
    pred=hd*index
    found=(permutedims(ud)*t.scalar)[1,:]
    InfoChevie("  # R^{",W,"}_{",L,"}\n")
    if pred!=found ChevieErr("quotients ",CycPol.(pred).//CycPol.(found),"\n")
    end
  end
end

test[:feginduce]=(fn=Tfeginduce,applicable=x->true, 
                  comment="check fakedegrees induce")

#---------------- test: invariants ------------------------

function Tinvariants(W)
  ii=invariants(W)
  vars=map(i->Symbol("x",i),1:rank(W))
  ii=map(f->f(Mvp.(vars)...),ii)*Int128(1)
  InfoChevie("  #")
  for i in eachindex(ii), j in eachindex(gens(W))
    InfoChevie("W.",j,"*I",i,",")
    if ^(ii[i],reflrep(W,i);vars)!=ii[i]  ChevieErr("not invariant\n") end
  end
  InfoChevie("\n")
end

test[:invariants]=(fn=Tinvariants,
  applicable=W->!(W isa Spets) && length(W)<14400, # H4 first painful client
                   comment="check invariants")

#---------------- test: classical ud fd  ------------------------

function Tudfdimprimitive(W)
  ud=CycPolUnipotentDegrees(W)
  uc=UnipotentCharacters(W)
  cs=first.(uc.charSymbols)
  vud=gendeg_symbol.(cs)
  for i in eachindex(cs)
    cmpcycpol(vud[i],ud[i];na=string("Deg",stringsymbol(cs[i])),nb="ud");
  end
  n=refltype(W)
  if isempty(n) return end # a corriger
  n=n[1]
  cs=first.(uc.almostCharSymbols[uc.almostHarishChandra[1][:charNumbers]])
  fd=CycPol.(fakedegrees(W,Pol()))
  if haskey(n,:orbit)&& n.orbit[1].series in [:D,:I,:B,:G] && order(n.twist)==2
       vfd=map(x->fegsymbol(x,1),cs)
  else vfd=fegsymbol.(cs)
  end
  for i in eachindex(cs) 
    cmpcycpol(vfd[i],fd[i];na=string("Deg",stringsymbol(cs[i])),nb="fd")
  end
end

test[:udfdimprimitive]=(fn=Tudfdimprimitive,
  applicable=function(W)
    n=refltype(W)
    if length(n)>1 return false
    elseif length(n)==0 return true end
    n=n[1]
    if haskey(n,:orbit) n.orbit[1].series in [:B,:C,:D,:I] && 
      order(n.twist)!=3 # type 2A excluded for now....???
    else n.series in [:A,:B,:C,:D,:G,:I] || (n.series==:ST && haskey(n,:p))
    end
  end,
comment="Test unideg and feg of classical Spets against formulas")

#---------------- test: families ------------------------
# test that a family satisfies Lusztig's "exotic fourier transform" properties.
# Families(W) or 
# Families(W,famno) or
# Families(family record [,opt])
# option: hard:=show details

function Tfamilies(W;hard=false)
  uc=UnipotentCharacters(W)
  for i in eachindex(uc.families) Tfamilies(W,i;hard) end
end

function Tfamilies(W,i;hard=false)
  uc=UnipotentCharacters(W)
  f=uc.families[i]
  v=" "*repr(W;context=rio())*".$i"
# else f:=arg[1]; arg:=arg{[2..Length(arg)]};v:=""; fi;
  if length(f.eigenvalues)==1 return end # nothing interesting to test
  O=toM(HasType.DiagonalMat(f.eigenvalues))
  if haskey(f,:sh) Sh=toM(HasType.DiagonalMat(f.sh)) end
  S=f.fourierMat
  if f isa Vector t=toM(f) end
  Sbar=conj(S)
  Id=one(S)
  tS=permutedims(S)
  v*="#$(size(S,1))"
  if haskey(f,:name) v*="("*Util.unicodeTeX(f.name)*")" end
  InfoChevie("  # family ",v)
  ps=SPerm(Sbar,S;dims=1)
  if isnothing(ps) ChevieErr("S* is not ps(S)\n| ") end
  wreal=haskey(f,:perm)
  if !haskey(f,:sh) && wreal!=all(>(0),vec(ps))
    error("weakly real=",wreal," but ps=",ps,"\n")
  end
  if wreal # S^2 is perm
    real=Perm(conj(f.eigenvalues),f.eigenvalues)
    if !isnothing(real) v*="+" else v*="P" end
  else real=false
  end
  if hard print("[",v,"] ...\n| ") end
  function check(a,msg)
    if !a InfoChevie("\n");ChevieErr("failed:",msg,"\n| ");
    elseif hard InfoChevie(" ",msg,",") 
    end
  end
  inds=1:length(f.eigenvalues)
  if @isdefined W
    ud=CycPolUnipotentDegrees(W)[f.charNumbers]
    a=unique(valuation.(ud))
    A=unique(degree.(ud))
    if length(a)!=1 || length(A)!=1 ChevieErr("a or A not constant");
    else a=a[1];A=A[1]
    end
    fd=fakedegrees(uc,Pol())[f.charNumbers]
    special=findall(==(a),valuation.(fd))
    cospecial=findall(==(A),degree.(fd))
    if length(special)!=1 || length(cospecial)!=1
      ChevieErr("special or cospecial not unique")
    else special=special[1];cospecial=cospecial[1]
    end
    if haskey(f,:special)
      check(f.special==special,".special=$(f.special) min b=$special")
    end
    if special!=cospecial
      if !haskey(f,:cospecial)
        ChevieErr(" .cospecial is not bound; should be ",cospecial,"\n| ")
      else check(f.cospecial==cospecial,
                 ".cospecial=$(f.cospecial) max B=$cospecial")
      end
    elseif haskey(f,:cospecial) && f.special!=f.cospecial
      ChevieErr("\n.cospecial=$(f.cospecial) should be $special\n")
    end
    check(conj(ud)==ud^ps,"ud*=ud^p")
    if real!==nothing && real!=false
      check(f.eigenvalues^(real/f.perm)==f.eigenvalues,"eig*=eig^perm")
    end
    if hard
      if S==Sbar print("\n|  S real")
      else print("\n|  S->S* is p=",ps,"\n| ")
      end
      if wreal # check ps is a sub-permutation of perm
        if !all(x->x^ps==x || x^ps==x^f.perm,inds)
          ChevieErr("*** S->S* is not a sub-perm of .perm\n| ")
        end
      end
    end
    # check Shintani preserves ud. We have Sh=O^-1 SO^-1
    ud=map(x->x(Pol()),ud)
    if haskey(f,:sh) check(Sh*S*ud==S*ud,"Shud=ud")
    elseif haskey(f,:lusztig) check((S*O)*ud==O*ud,"Shud=ud")
    else check(S*O^-1*ud==O*ud,"Shud=ud")
    end
  elseif haskey(f,:special) special=f.special
  else ChevieErr(".special not bound\n| ");special=1
  end
  check(!(0 in S[special]),"not 0 in S[special]")
  check(S*permutedims(Sbar)==Id,"S unitary")
  if wreal P=Id^f.perm;check(S*P==P*S,"[S,P]=1") end
  if haskey(f,:sh)
    check((O*tS*Sh^-1*S)^2==Id,"(O*tS*Sh-1*S)^2=1")
  else 
    if haskey(f,:lusztig)
      if !f.lusztig error(".lusztig bound but false") end
      if !wreal error(".lusztig bound but not .perm") end
      check((O*S*P)^3==Id,"(OSP)^3=1")
    else check((O*S)^3==Id,"(OS)^3=1")
    end
    check(S==tS,"\n|  S symmetric")
    check(O*S^2==S^2*O,"[S^2,O]=1")
  end
  check((S*tS)^2==Id,"(S*tS)^2=1")
# Print("\n| ");
# check(O=S*O*S,"O=SOS");
# if real then Print(" evalf(S[special]): ",
#     List(S[special]/S[special][special],evalf));
# fi;
  if hard print("\n| ") end
  if length(S)<=40 || hard
    s=map(x->x//x[special],eachcol(Sbar))
    _max=0
    for i in inds
      if length(S)>10 || hard print(".") end
      for j in inds
        v=map(k->S[i,k]*S[i,k],inds)*s
        _max=max(_max,count(!=(0),v))
     #  Print("e[",i,"]*e[",j,"]=",v,"\n");
        if !all(isinteger,v)
          error("not integral:",i,",",j,":",v,"\n| ")
        end
        p=all(>=(0),v)
        if wreal && !p
          if all(<=(0),v)
               ChevieErr("family is wreal but Sums[",i,",",j,"]<=0\n| ")
          else ChevieErr("family is wreal but Sums[",i,",",j,"]=",v,"\n| ")
          end
        end
        if hard print(p ? "+" : ".") end
      end
    end
  # Print("max=",_max,"\n");
  end
  if hard && length(S)<50 print("\n") end
  InfoChevie("\n")
end

test[:families]=(fn=Tfamilies,applicable=isspetsial,comment="testing families")
#------------------------- HCinduce gendegs -----------------------------

function TdegsHCInduce(W)
  for J in parabolic_reps(W) TdegsHCInduce(W,J) end
end

function TdegsHCInduce(W,J)
  q=Pol()
  if W isa Spets L=subspets(W,J)
  else L=reflection_subgroup(W,J)
  end
  index=exactdiv(generic_order(W,q),generic_order(L,q))
  index=shift(index,-index.v)
  if W isa Spets index*=generic_sign(L)//generic_sign(W) end
  uL=UnipotentCharacters(L)
  pred=degrees(uL,q)*index
  tbl=Lusztig.HCInductionTable(L,W)
  ind=map(x->UniChar(W,x),eachcol(tbl.scalar))
  inddeg=improve_type(degree.(ind))
  if pred!=inddeg
    nh=charnames(uL)
    ChevieErr("Relgroups:",tbl.pieces)
    cmpvec(pred,inddeg;na="udLevi*index",nb="viaHCInductionTable")
  end
end

test[:degsHCinduce]=(fn=TdegsHCInduce,applicable=isspetsial,comment=
  "unidegs are consistent with HC induction from split levis")

#------------------------- Curtis duality -----------------------------
# For a Group generated by true reflections,
# check that deg(rho_chi)=q^N deg(rho_{\chi tensor sign})(q^-1)
function Tcurtisduality(W)
  q=Pol()
  ud=UnipotentCharacters(W)
  ud=degrees(ud,q)[ud.harishChandra[1][:charNumbers]]
  p=detPerm(W)
  N=sum(degrees(W)-1)
  cmpvec(ud^p,map(x->q^N*x(q^-1),ud);
         na="ud sign permuted",nb="computed Curtis dual ud")
end

test[:curtisduality]=(fn=Tcurtisduality,applicable=
  W->!(W isa Spets) && all(x->order(x)==2,gens(W)) && isspetsial(W),
  comment="check degrees of Curtis duals")

#------------------------- AlmosHC -----------------------------
# check UnipotentCharacters(W).almostHarishChandra

function TalmostHC(W::Spets)
  uc=UnipotentCharacters(W)
  for h in uc.almostHarishChandra
    r=relative_coset(W,h[:levi])
    t=refltype(r)
    if length(t)!=length(h[:relativeType]) error() end
    for i in eachindex(t)
      if length(t[i].orbit)>1 error("not irred") end
      tt=t[i].orbit[1]
      if !haskey(h[:relativeType][i],:orbit)
        if t[i].twist!=Perm() error("should be twisted") end
        hh=h[:relativeType][i]
      else hh=h[:relativeType][i].orbit[1]
      end
      if hh.series!=tt.series
          ChevieErr("series do not match: is ",hh.series,
                   " should be ",tt.series,"\n") end
      if hh.indices!=Group(r).relativeIndices[tt.indices]
        ChevieErr("indices do not match: are ",hh.indices,
                 " should be ",Group(r).relativeIndices[tt.indices],"\n") end
    end
  end
end

function TalmostHC(W)
  uc=UnipotentCharacters(W)
  for h in uc.almostHarishChandra
    if W isa CoxeterGroup r=relative_group(W,h[:levi])
    else r=relative_group(W,h[:levi],vcat(map(
             x->vcat(map(y->y.indices,x.orbit)...),h[:relativeType])...))
    end
    t=refltype(r)
    if length(t)!=length(h[:relativeType]) error() end
    for i in eachindex(t)
      tt=t[i]
      if h[:relativeType][i].orbit[1].series!=tt.series
       ChevieErr("stored series is ",h[:relativeType][i].orbit[1].series,
                " should be ",tt.series,"\n") end
      if !haskey(tt,:indices)
        ChevieErr("could not find indices (expected ",
                  h[:relativeType][i].orbit[1].indices,")")
      elseif maximum(tt.indices)>length(r.relativeIndices)
        ChevieErr("use non-simple roots:",tt.indices," of ",r,
            "<",join(r.relativeIndices,","),">\n")
      elseif h[:relativeType][i].orbit[1].indices!=r.relativeIndices[tt.indices]
        ChevieErr("stored indices are ",
                h[:relativeType][i].orbit[1].indices,
               " should be ",r.relativeIndices[tt.indices],"\n")
      end
    end
  end
end

test[:almostHC]=(fn=TalmostHC,applicable=isspetsial,comment="almostHC series")

#------------------------- sum squares -----------------------------
# The sum of squares of fakedegrees should equal the sum of
# normsquares of unipotent degrees
function Tsumsquares(W)
  q=Pol()
  if sum(x->x^2,fakedegrees(W,q))!=
     sum(x->x*conj(x),degrees(UnipotentCharacters(W),q))
    ChevieErr("fails\n")
  end
end

test[:sumsquares]=(fn=Tsumsquares,applicable=isspetsial,comment="Sum squares")

#------------------------- aA -----------------------------
# check that stored a and A are correct
function TaA(W)
  uc=UnipotentCharacters(W)
  q=Pol()
  ud=degrees(uc,q)
  if haskey(uc,:a)
    cmpvec(uc.a,valuation.(ud);na="stored a",nb="computed a")
  end
  if haskey(uc,:A)
    cmpvec(uc.A,degree.(ud);na="stored A",nb="computed A")
  end
end

test[:aA]=(fn=TaA,applicable=isspetsial,comment="aA")

#------------------------- eigen -----------------------------
# check eigenvalues in families agree with those in HC series

function Teigen(W)
  uc=UnipotentCharacters(W)
  e=eigen(uc)
  for i in eachindex(uc.harishChandra)
    f=uc.harishChandra[i]
    for j in eachindex(f[:charNumbers])
      n=f[:charNumbers][j]
      if e[n]!=f[:eigenvalue]
        ChevieErr("HCfamily ",i,"#",j,":",charnames(uc;TeX=true)[n],"=",n,
          " eigen from fam.=",e[n]," but from HC=",f.eigenvalue,"\n")
      end
    end
  end
end

test[:eigen]=(fn=Teigen,applicable=isspetsial,comment="eigen")

#------------------------- qeigen -----------------------------
## check fractional eigenvalues of Frobenius wrt. values stored in HC series

function Tqeigen(W)
  uc=UnipotentCharacters(W)
  e=Uch.qeigen(uc)
  for i in eachindex(uc.families)
    f=uc.families[i]
    if haskey(f,:qEigen) q=f.qEigen
    else q=fill(0,length(f))
    end
    for j in 1:length(f)
      n=f.charNumbers[j]
      if e[n]!=q[j] ChevieErr("HCfamily ",i,"#",j,":",
          charnames(uc;limit=true)[n],"=",n,
          " eigen from HC.=",e[n]," but from fam=",q[j],"\n")
      end
    end
  end
end

test[:qeigen]=(fn=Tqeigen,applicable=isspetsial,comment="qeigen")

#------------------------- braidrel -----------------------------

# CheckRelation(gens,rel[,f])
# Check  that  the  homogeneous  relation  rel[1]=rel[2] holds for gens (in
# particular  that no  left factor  homogeneous relation  holds). If given,
# call f in case of failure.
function check_relation(gens,rel;f=x->x)
  p=r->string(joindigits(r[1]),"=",joindigits(r[2]))
  l=gens[rel[1][1]]
  r=gens[rel[2][1]]
  i=1
  while i<length(rel[1])
    if l==r
      f(" relation ",p(rel)," already holds as ",p(map(x->x[1:i],rel)))
      return false
    end
    i+=1
    l*=gens[rel[1][i]]
    r*=gens[rel[2][i]]
  end
  if l==r return true end
  f(" relation ",p(rel),"failed")
  false
end

function Tbraidrel(W,r=braid_relations(W))
  for rel in r check_relation(gens(W),rel;f=ChevieErr) end
end

test[:braidrel]=(fn=Tbraidrel,applicable=W->!(W isa Spets),comment=
"W satisfies braid relations of type t")

#------------------------- root system -----------------------------

# Check that all elements of the list of vectors l have coefficients
#  in the ring R when expressed as linear combinations of l{indices}
#function rootsystem(W)
#  integrality:=function(l,indices,R)local rb;
#    rb:=List(l,r->SolutionMat(l{indices},r));
#    return ForAll(Flat(rb),y->y in R);
#  end
#  Coroots:=W->List([1..Length(W.roots)],i->PermRootOps.Coroot(W,i));
## Find representatives up to scalar multiple of elements of the list of
## vectors vec. Check that other elements differ from such a representative
## by a unit of the ring R
#  replines:=function(vec,R)local found,res,v,p;res:=[];
#    for v in vec do found:=false;
#      for x in res do p:=ProportionalityCoefficient(x,v);
#        if p<>false then 
#	  if not (p in R and 1/p in R) then return false;fi;
#	  found:=true;
#        fi;
#      od;
#      if not found then Add(res,v);fi;
#    od;
#    return res;
#  end;
#  R:=DefaultRing(Flat(CartanMat(W)));
#  try:=Combinations(W.generatingReflections,W.rank);
#  p:=PositionProperty(try,v->integrality(W.roots,v,R));
#  if p=false then ChevieErr("found no root basis");else;fi;
#  cr:=Coroots(W);
#  p:=PositionProperty(try,v->integrality(cr,v,R));
#  if p=false then ChevieErr("found no coroot basis");else;fi;
#  p:=replines(W.roots,R);
#  if p=false then ChevieErr("not reduced");fi;
#  m:=replines(cr,R);
#  if m=false then ChevieErr("coroots not reduced??");fi;
#  if Length(m)>0 and not ForAll(Set(Flat(p*TransposedMat(m))),x->x in R) then 
#    ChevieErr("not a root system");fi;
#end,
#W->not IsSpets(W),
#"Check that W.roots define a distinguished root system",
#" in the sense of Broue-Corran-Michel");
#
#CHEVIE.AddTest("GaloisAutomorphisms",
#function(W)local k,Wk,m,g,gm,p;
#  k:=Field(Flat(W.matgens));
#  Wk:=Field(Flat(CartanMat(W)));
#  if k<>Wk then 
#    ChevieErr("k_W=",Wk," but matrices over ",k,"\n");
#  fi;
#  if k=Rationals then return;fi;
#  for g in GaloisGroup(k).generators do
#    for m in W.matgens do
#      gm:=List(m,x->OnTuples(x,g)); p:=PermMatX(W,gm);
#      if not p in W or gm<>MatXPerm(W,p) then 
#        ChevieErr("not Galois stable\n");fi;
#    od;
#  od;
#end,
#W->not IsSpets(W),
#"check that W's reflection representation is globally",
#"invariant by Gal(k_W/Q)");

#function frombraidrel(W)
#  n=nbgens(W)
#  F=FreeGroup(n)
#  r=List(BraidRelations(W),x->List(x,y->Product(y,z->F.(z))))
#  r=List(r,x->x[1]/x[2])
#  Append(r,List([1..n],i->F.(i)^OrderPerm(W.(i))))
#  Size(W)=Size(F/r)
#end
#W->not IsSpets(W) and Size(W)<64000,
#"check that the finitely presented group defined by the",
#" braid and order relations has the expected Size");

#------------------------- very good classreps -----------------------------

function Tclassreps(W)
  cl=classinfo(W)
  for i in 1:nconjugacy_classes(W)
    ww=cl[:classtext][i]
    if W isa CoxeterGroup
      w=BraidMonoid(W)(ww...)
      o=cl[:orders][i]
      l=Brieskorn_normal_form(w^o)
    elseif W isa CoxeterCoset
      w=BraidMonoid(Group(W))(ww...)
      wF=W(ww...)
      o=findfirst(j->(wF)^j==W.phi^j,1:order(wF))
      l=Brieskorn_normal_form(twisted_power(w,o,Frobenius(W)))
    end
    if isodd(length(l)) || any(j->l[2*j-1]!=l[2*j],1:div(length(l),2)) || 
      any(j->!issubset(l[j+1],l[j]),1:length(l)-1)
      ChevieErr("class ",i," is not good\n")
    end
    if iseven(o)
      if W isa CoxeterGroup l=Brieskorn_normal_form(w^div(o,2))
      elseif W isa CoxeterCoset
	l=Brieskorn_normal_form(twisted_power(w,o,Frobenius(W)))
      end
      if any(j->!issubset(l[j+1],l[j]),1:length(l)-1)
	ChevieErr("class ",i," is not very good\n")
      end
    end
  end
end

test[:classreps]=(fn=Tclassreps,applicable=W->(W isa CoxeterGroup) ||
                  W isa CoxeterCoset, comment=
"check that classreps are very good in the sense of Geck-Michel")

#------------------------- minuscule weights -----------------------------

function Tminusculeweights(W)
  w=weightinfo(W)
  l=map(refltype(W)) do t
    n=t.indices
    r=map(i->sum(roots(W,i)[n]),1:2*nref(W))
    r=argmax(r)
    vcat(filter(i->roots(W,r)[i]==1,n),[0])
  end
  l=map(x->filter(!iszero,x),cartesian(l...))
  l=filter(!isempty,l)
  if l!=w[:minusculeCoweights] ChevieErr("minuscule coweights") end
  l=map(refltype(W))do t
    n=t.indices
    m=toM(coroots.(Ref(W),n))
    r=map(i->sum(solutionmat(m,coroots(W,i))[n]),1:2*nref(W))
    r=argmax(r)
    vcat(filter(i->solutionmat(m,coroots(W,r))[i]==1,n),[0])
  end
  l=map(x->filter(!iszero,x),cartesian(l...))
  l=filter(!isempty,l)
  if l!=w[:minusculeWeights] ChevieErr("minuscule weights") end
end

test[:minusculeweights]=(fn=Tminusculeweights,applicable=isweylgroup,
  comment="correspond to coeff 1 of highest coroot on simple")

#------------------------- CharTable(3D4) -----------------------------

#  2 methods to check CharTable(hecke(3D4,q))
function Thecke3d4(triality)
  F4=coxgroup(:F,4)
  q=Pol()
  T=Tbasis(hecke(F4,[q,q,1,1]))
  v=[T(2),T(3,2,3),T(1),T(4,3,2,3,4)] # embedding of Hecke(triality,q)
  m=permutedims(toM(map(x->char_values(prod(v[x];init=T())*T(4,3)),
                        classinfo(triality)[:classtext])))
  WF=spets(F4)
  tbis=subspets(WF,[2,9,1,16],F4(4,3)) # inside F4
  ct=CharTable(hecke(triality,q)).irr
  cmpvec(toL(ct), toL(m[findfirst.(==(1),
        toL(permutedims(InductionTable(tbis,WF).scalar))),:]);
   na="charTable(Hecke(3D4))",nb="computed by coset induction")
  D4=reflection_subgroup(F4,[2,9,1,16])
  v=permutedims(InductionTable(D4,F4).scalar)
  v=v[charinfo(triality)[:charRestrictions],:]
  cmpvec(toL(ct),toL(-m[findfirst.(==(2),toL(v)),:]);
   na="charTable(Hecke(3D4))",nb="computed by subgroup induction");
end

test[:hecke3d4]=(fn=Thecke3d4,
  applicable=function(W)
    t=refltype(W)
    if length(t)!=1 return false end
    t=only(t)
    haskey(t,:orbit) && t.orbit[1].series==:D && order(t.twist)==3
  end,
    comment="3d4")

#------------------------- G4_22indexclasses -----------------------------

function TG4_22index(W)
  t=refltype(W)[1]
  O=ComplexReflectionGroup(getchev(t,:Generic))
  e=getchev(t,:Embed)
  c=map(c->vcat(map(x->e[x],c)...),classinfo(W)[:classtext])
  c=map(x->position_class(O,O(x...)),c)
  l=map(x->findfirst(i->reflection(O,i)==O(x...),eachindex(roots(O))),e)
  a=classinfo(W)[:indexclasses]
  if nothing in l # for G13, G15 gen is square of gen of G11
    b=c
  else 
    W=reflection_subgroup(O,l)
    b=fusion_conjugacy_classes(W,O)
  end
  if a!=b || b!=c
    a=hcat(filter(x->length(unique(x))>1,collect(eachcol(toM([a,b,c]))))...)
    showtable(stdout,a[[2,3],:];
      row_labels=[string("fusion to ",repr(O;context=rio())),"classtext"],
      col_labels=string.(a[1,:]),rows_label="indexclasses")
  end
end

test[:G4_22index]=(fn=TG4_22index,applicable=
function(W)t=refltype(W)
  length(t)==1 && haskey(t[1],:ST) && t[1].ST in 4:22
 end, comment="Check classinfo(W).indexclasses for W in G₄-G₂₂")

#------------------------- Opdam -----------------------------

# c_\chi in Gordon-Griffeth
function coxeter_number(W,i)
  cl=classinfo(W)[:classes]
  if cl==[1] return 0 end
  ct=CharTable(W).irr
  Int(sum(degrees(W)-1)-sum(h->sum(c->exactdiv(cl[c]*ct[i,c],ct[i,1]),
                                   h.cl_s),hyperplane_orbits(W)))
end

function Thgal(W)
  ct=CharTable(W).irr
  complex=Perm(ct,conj(ct);dims=1)
  x=Mvp(:x)
  d=degrees(W)
  z=isempty(d) ? 1 : gcd(d)
  H=hecke(W,x^z)
  ct=CharTable(H).irr
  gal=z==1 ? Perm() : Perm(ct,map(u->u(;x=x*E(z)),ct))
  fd=fakedegrees(W,x)
  ci=charinfo(W)
  if haskey(ci,:hgal) hgal=ci[:hgal] else hgal=Perm() end
  if isnothing(gal) ChevieErr("could not compute gal")
  elseif gal!=hgal ChevieErr("hgal wwwwwrong")
  end
  # following test is [Malle, Rationality... Theorem 6.5]
  for i in 1:nconjugacy_classes(W)
    if fd[i^(hgal*complex)]!=x^coxeter_number(W,i)*fd[i](x=x^-1)
      ChevieErr("hgal wrong for ",i,"\n")
    end
  end
end

test[:hgal]=(fn=Thgal,applicable=W->!(W isa Spets),comment="Opdam")

#------------------------- ParameterExponents -----------------------------
# check that the parameter exponents of the relative hecke algebras
# agree with unipotent degrees corresponding to cyclic
# relative groups of sub-series

function Tparameterexponents(W)
  ser=UnipotentCharacters(W).harishChandra
  for i in eachindex(ser) Tparameterexponents(W,i) end
end

function Tparameterexponents(W,i)
  h=UnipotentCharacters(W).harishChandra[i]
  I=indices(h[:relativeType])
  if W isa Spets I=map(x->orbit(W.phi,x),I)
  else 
    G=relative_group(W,h[:levi])
    if haskey(G,:relativeIndices) && G.relativeIndices[indices(refltype(G))]!=I 
      ChevieErr("indices computed=",
        joindigits(G.relativeIndices[indices(refltype(G))]),
        " stored=",joindigits(I),"\n");
    end
    I=map(x->[x],I)
  end
  for i in eachindex(I)
    L=reflection_subgroup(W,vcat(h[:levi],I[i]))
    t=refltype(L)
    H=reflection_subgroup(L,restriction(L,inclusion(W,h[:levi])))
    InfoChevie("  # ParameterExponents from ",H,":",I[i])
    if !haskey(t[1],:indices) && !haskey(t[1],:orbit)
      ChevieErr("Levi ",vcat(h[:levi],I[i])," could not be identified\n")
    else
      if !(L isa Spets) L=spets(L)
        H=subspets(L,restriction(L,inclusion(W,h[:levi]))) 
      end
      h1=copy(h)
 #    h1[:levi]=restriction(H,inclusion(W,h[:levi]))
      h1[:levi]=1:ngens(Group(H))
      hh=Lusztig.FindSeriesInParent(h1,H,L,UnipotentCharacters(L).harishChandra).ser
      if length(hh[:charNumbers])!=2
        cusp=h[:cuspidalName]
        cusp=cusp=="" ? "." : cusp
        s=Series(L,H,findfirst(==(cusp),
                     charnames(UnipotentCharacters(H);TeX=true)),1;NC=true)
        if isnothing(dSeries.char_numbers(s)) || isnothing(dSeries.fill(s))
          error("could not fill ",s)
        else 
          if hh[:parameterExponents][1] isa Vector
            exp=hh[:parameterExponents][1]
          else
            exp=fill(0,s.e)
            exp[1]=hh[:parameterExponents][1]
          end
          if haskey(s,:translation)
            t=filter(0:s.translation:s.e)do d
              local v
              v=1 .+map(x->x%s.e,(1:s.e).+d)
              s.mC[v]==exp && s.charNumbers[v]==hh.charNumbers
            end
          else t=[1]
          end
          if length(t)==0 error("unexpected") end
#    Ok("decs=",t);
        end
      else
        ud=degrees(UnipotentCharacters(L),Pol())
        ud=exactdiv(ud[hh[:charNumbers][1]],ud[hh[:charNumbers][2]])
        if length(ud.c)!=1 || ud.c[1]!=1
          ChevieErr("not monomial")
        elseif h[:parameterExponents][i]!=valuation(ud)
          ChevieErr(L,": wrong parameter ",
            h.parameterExponents[i]," instead of ",valuation(ud),"\n");
        end
#         Ok("=",ud.valuation);
      end
    end
  end
end

test[:parameterexponents]=(fn=Tparameterexponents,applicable=isspetsial,
  comment="of rel. Hecke algebras agree with cyclic formula")

#------------------------- discriminant -----------------------------

# for a reflection group of rank r: Discriminant(G)
#          returns a list of linear factors as `Mvp`s in x1,x2,...,xr
function reflection_discriminant(W)
  res=[]
  for h in hyperplane_orbits(W)
    cr=coroots(W,h.s)
    for w in map(x->transporting_elt(W,gens(W)[h.s],x),
                 conjugacy_class(W,h.cl_s[1]))
      append!(res,map(i->reflrep(W,w^-1)*cr,1:h.order))
    end
  end
  vars=map(i->Mvp(Symbol("x",i)),eachindex(degrees(W)))
  map(x->sum(map(*,vars,x)),res)
end

function Tdiscriminant(W)
  r=prod(reflection_discriminant(W).*big(1);init=Mvp(1))
  j=discriminant(W)
  if isnothing(j) InfoChevie("not implemented\n");return end
  ii=invariants(W)
  ii=map(x->x(map(i->Mvp(Symbol("x",i))*big(1),1:rank(W))...),ii).*big(1)
  j=j(ii...)
  if r*last(first(j.d))!=j*last(first(r.d)) ChevieErr("disagrees\n") end
end

test[:discriminant]=(fn=Tdiscriminant,
  applicable=W->!(W isa Spets) && length(W)<1152, # F4 first painful clientend
  comment="discriminant")
end
