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
    if a[j]==-b[j] ChevieErr(a[j],"=-",b[j]);continue end
    t=findall(==(a[j]),b)
    if length(t)>0 ChevieErr(a[j]," found at ",t);continue end
    t=findall(==(-a[j]),b)
    if length(t)>0 ChevieErr("-",a[j]," found at ",t);continue end
    ChevieErr(na,"[$j]=",a[j]," not found")
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

#---------------- test: lusztiginduction ------------------------
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
  if haskey(t,:scalars) 
    print("\t****** scalars=",t.scalars) 
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
    :uNames=>UnipotentCharacters(L).TeXCharNames,
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
function checkCharTable(W)
  G=Group(gens(W))
  G.classreps=classreps(W)
  ct=CharTable(G).irr
  ct1=CharTable(W).irr
  p=Perm(ct,ct1;dims=1)
  if isnothing(p) error("irreducibles") 
  else println("permuted ",p)
  end
end

test[:chartable]=(fn=checkCharTable, applicable=W->!(W isa Spets),
   comment="CharTable")

#----------------test: powermaps ------------------------

function powermaps(W)
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

test[:powermaps]=(fn=powermaps, applicable=W->!(W isa Spets),
                  comment="powermaps")

#---------------- test: positionclasses ------------------------
function positionclasses(W)
  cl=map(x->Gapjm.position_class(W,x),class_reps(W))
  if cl!=1:length(cl) error("classes") end
end

test[:positionclasses]=(fn=positionclasses, applicable=W->true,
   comment="classreps")

#---------------- test: unipotentclasses ------------------------
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

#---------------- test: charparams ------------------------
"""
Jean Michel and Ulrich Thiel (2017)

this function checks that the parameters/names for characters of Complex
eflection groups agree with the description whcih is since 2016 in the CHEVIE 
manual. If not, it tells what permutation of the data is needed.
"""
function charparams(W)
  ct=CharTable(W).irr
  fd=fakedegrees(W,Pol())
  db=map(x->[x(1),valuation(x)],fd)
  n=sprint(show,W;context=:limit=>true)
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
     " in the induced of Id from ", sprint(show,L;context=:limit=>true))
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

test[:charparams]=(fn=charparams, applicable=W->W isa Group,
   comment="Check charparams for consistency with Michel/Thiel rules")

#---------------- test: HCdegrees ------------------------
function HCdegrees(W)
  uw=UnipotentCharacters(W)
  for i in eachindex(uw.harishChandra)
   if !isempty(uw.harishChandra[i][:relativeType]) HCdegrees(W,i) end
  end
end

# HCdegrees(W [,ser [,guess]])
# Checks that unipotent degrees from Harish-Chandra series no. i of the Spets
# W agree with degrees coming from Schur elements of relative Hecke algebra
function HCdegrees(W,i,rel=false)
  uw=UnipotentCharacters(W)
  hw=uw.harishChandra[i]
  L=reflection_subgroup(W,hw[:levi])
  index=CycPol(generic_order(W,Pol())/generic_order(L,Pol()))
  index*=CycPol(Pol()^-valuation(index))
  n=hw[:cuspidalName]
  if n=="" n="." end
  cusp=CycPolUnipotentDegrees(L)[Lusztig.FindCuspidalInLevi(n,L)]
  xprintln("#HC_",L,"(cusp=",n,":",cusp,")[G:L]=",index)
  R=reflection_group(hw[:relativeType]...)
# check parameters of relative algebra by the formula
# u_{s,j}=ζ_s^j  q^(([a+A]_ρ_{det_s^j}-[a+A]_ρ_{\Id})/order#class(s)
# where ρ_χ is the unip. char corresp. to χ in relative group
  aA=map(x->uw.a[x]+uw.A[x],hw[:charNumbers])
  para=map(hyperplane_orbits(R))do h
    l=map(i->(aA[i]-aA[charinfo(R)[:positionId]])//h.order//h.N_s,h.det_s)
    max(maximum(l),0).-vcat([0],l)
  end
  ss=simple_representatives(R)
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
  reldeg=map(x->descent_of_scalars(cusp*index//x,den),ud[hw[:charNumbers]])
  H=Uch.relative_hecke(uw,i,Mvp(:x)^den)
  InfoChevie("   Relative ",R,"\n");
  sch=schur_elements(H)
  if any(==(false),sch)
    ChevieErr("Schur elements not implemented for ",H,"!")
    return false
  end
  sch=CycPol.(sch)
  cmpvec(sch,reldeg;na=string("Schur(",sprint(show,R;context=rio()),")"),nb="ud")
  if rel reldeg
  else  Ref(descent_of_scalars(cusp*index,den)).//(sch*1//1)
  end
end

test[:HCdegrees]=(fn=HCdegrees, applicable=isspetsial,
 comment="check generic degrees from relative Hecke algebras of HCinduction")

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
  zeta=E(d)
  W=H.W
  ct=CharTable(H).irr
  ct1=CharTable(W).irr
  ct2=improve_type(scal(value.(ct,Ref(:q=>zeta))))
  n=axes(ct2,1)
  good=filter(i->!any(x->x isa HasType.Unknown,ct2[:,i]),n)
  p=Perm(ct1[:,good],ct2[:,good],dims=1) #Permuted(ct,p) specializes
  if !isone(p) println("***** perm=",p) 
    if iscyclic(W) ChevieErr("should not have perm")end 
  end
  d1=exponent(d)//gcd(conductor(d),gcd(degrees((W))))
  i=position_regular_class(W,d1) # this class represents F
  omegachi=map(x->scal(value(x[i]/x[1],:q=>1)),eachrow(^(ct,p,dims=1)))
  frac=degree.(central_monomials(H)).*d1
  om2=map((o,p)->o*d^-p,map(x->x[i]/x[1],eachrow(ct1)),frac)
  if omegachi!=om2 
    omegachi=om2
    ChevieErr(classinfo(W)[:classtext][i],"^",1//d1," not equal to π(",W,")\n")
  end
  ss=CycPol.(schur_elements(H)^p)
  ss=map(x->s.degree/x,ss)
# omegachi*=Eigenvalues(UnipotentCharacters(s.levi))[s.cuspidal]
  zeta=E(Root1(;r=d1)^frac[PositionId(W)])
  if haskey(s,:delta) && s.delta!=1 && iscyclic(W)
    omegachi=map(i->Root1(;r=s.d*s.e*s.delta*
              ((i-1)/s.e-dSeries.mC(s)[i]*s.d)),1:s.e)
    zeta=Root1(;r=s.d*s.e*s.delta*s.d*dSeries.mC(s)[1])
  end
  map((deg,eig,frac)->(deg=deg,eig=eig*zeta,frac=Mod1(frac)),ss,omegachi,frac)
end

#---------------- test: series ------------------------
function CheckSerie(s)
  W=s.spets
  InfoChevie("\n   # ",s)
  relative_group(s)
  if s.principal
    Ed=s.d
    c=(generic_order(s.spets,Pol())//sum(x->x^2,s.WGLdims)//
       generic_order(s.levi,Pol()))(Ed)
    if c!=1 ChevieErr(s," => |G|/(|L||WGL|) mod Φ=",c,"\n") end
    e=generic_sign(s.spets)*Ed^(sum(codegrees(Group(s.spets)))-
       sum(codegrees(Group(s.levi))))/generic_sign(s.levi)
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
      c=Set(map((x,y)->x/y,pred,e))
      if length(c)==1 ChevieErr("predicted=",c[1],"*actual")
      else cmptables(
                 (rowLabels=["actual   "],columnLabels=n,scalar=[e]),
                 (rowLabels=["predicted"],columnLabels=n,scalar=[pred]))
      end
    end
  end
end

# checkSeries(WF[,d[,ad]])
function checkSeries(W,arg...)
  if length(arg)>=1 l=filter(s->s.levi!=s.spets,Series(W,arg...))
  else l=dSeries.ProperSeries(W)
  end
  for s in l CheckSerie(s)end
  InfoChevie("\n   ")
  l
end

test[:series]=(fn=checkSeries, applicable=isspetsial,
 comment="check d-HC series")

#---------------- test: extrefl ------------------------
using LinearAlgebra: tr
# test ReflectionEigenvalues, ChevieCharInfo.extRefl, .positionDet, .positionId
function checkextrefl(W)
  ct=CharTable(W)
  # compute first using ReflectionEigenvalues
  n=nconjugacy_classes(W)
  v=reverse(permutedims(toM(map(r->prod(n->Pol()+E(n),r).c,refleigen(W))));dims=1)
  # check v[2,:] using reflrep
  if v[2,:]!=map(w->tr(reflrep(W,W(w...))),classinfo(W)[:classtext])
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

test[:extrefl]=(fn=checkextrefl, applicable=x->x isa Group,
 comment="check extrefl")

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
  for d in sort(unique(conductor.(eig)),rev=true)
    p=findfirst(==(Root1(d,1)),eig)
    if p!==nothing && mul[p]>0
      for i in 1:div(d*mul[p],length(W)) push!(res,d) end
      n=mul[p]
      for j in 0:d-1 mul[findfirst(==(Root1(d,j)),eig)]-=n end
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
  # more exactly a list of pairs [eig in Q/Z, mul/|W|]
  searchdeg=function(e,card,degs)local res,d,f,g,p,pos,ne
#   @show degs,card
    if isempty(degs) return [Pair{Int,Root1}[]] end
    e=filter(x->x[2]!=0,e)
    res=Vector{Pair{Int,Root1}}[]
    d=degs[1]
    f=findall(i->last(i)>=1//d,e)
    g=first.(filter(x->first(x).r<1/d,e[f]))
    for p in g
      pos=map(i->findfirst(j->first(e[j])==p*i,f),Root1.(d,0:d-1))
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
  map(x->(x[1],E(x[2])),only(e))
end

function checkrefldegrees(W)
  d=sort(degrees(W))
  d1=sort(reflectiondegrees(W))
  cmpvec(d,d1;na="degrees",nb="reflectiondegrees")
  if isweylgroup(W) && length(refltype(W))==1 && rank(W)==semisimplerank(W)
    d1=last.(tally(sum.(W.rootdec[1:nref(W)])))
    d1=sort(1+conjugate_partition(d1))
    cmpvec(d,d1;na="degrees",nb="dual partition of root heights")
  end
end

test[:degrees]=(fn=checkrefldegrees,applicable=x->true, comment="check degrees")

#---------------- test: fakedegrees ------------------------

function checkfeg(W)
  q=Pol()
  fd=fakedegrees(W,q)
  ffd=fakedegrees(W,q;recompute=true)
  cmpvec(ffd,fd;na="recomputed",nb="fakedegrees")
  cmpvec(valuation.(fd),charinfo(W)[:b];na="computed b",nb="charinfo b")
  cmpvec(degree.(fd),charinfo(W)[:B];na="computed B",nb="charinfo B")
end

test[:feg]=(fn=checkfeg,applicable=x->true, comment="check fakedegrees")

#---------------- test: fakedegrees induce ------------------------

function fakedegreesinduce(W)
  for J in parabolic_reps(W) fakedegreesinduce(W,J) end
end

function fakedegreesinduce(W,J)
  q=Pol()
  if !(W isa Spets) W=spets(W) end
  for L in twistings(W,J)
    t=InductionTable(L,W)
    hd=fakedegrees(L,q)
    ud=fakedegrees(W,q)
    index=generic_sign(L)*generic_order(W,q)/generic_order(L,q)/generic_sign(W)
    index=shift(index,-valuation(index))
    pred=hd*index
    found=(permutedims(ud)*t.scalar)[1,:]
    InfoChevie("\n   # R^{",W,"}_{",L,"}")
    if pred!=found ChevieErr("quotient ",
                map((a,b)->a*inv(b),CycPol.(pred),CycPol.(found)),"\n")
    end
  end
end

test[:feginduce]=(fn=fakedegreesinduce,applicable=x->true, 
                  comment="check fakedegrees induce")

#---------------- test: invariants ------------------------

function checkinvar(W)
  ii=invariants(W)
  vars=map(i->Mvp(Symbol("x",i)),1:rank(W))
  ii=map(f->f(vars...),ii)
  InfoChevie(" #")
  for i in eachindex(ii), j in eachindex(gens(W))
    InfoChevie("W.",j,"*I",i,",")
    if ii[i]^reflrep(W,i)!=ii[i]  ChevieErr("not invariant\n") end
  end
end

test[:invariants]=(fn=checkinvar,
                   applicable=W->!(W isa Spets) && length(W)<14400, 
                   comment="check invariants")
#x->not IsSpets(x) and Size(x)<14400); # H4 first painful client

end
