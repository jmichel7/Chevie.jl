# Programs for making CharTables of Hecke algebras

#schur relations
function schurrel(H;p=Perm())
  s=CycPol.(invpermute(schur_elements(H),p))
  ll=lcm(s)
  s=ll.//s
  cs=map(x->x(Pol()),s)
  ct=CharTable(H)
  good=filter(i->!any(ismissing,ct.irr[:,i]),1:nconjugacy_classes(H.W))
  sums=map(i->sum(ct.irr[:,i].*cs),good)
  good[findall(i->(good[i]!=1 && !iszero(sums[i]))||
                  (good[i]==1 && sums[i]!=ll(Pol())),eachindex(good))]
end
 
# classtext of "canonical" generator of the center
function centerword(W)
  if length(refltype(W))>1 error() end
  classinfo(W).classtext[position_regular_class(W,gcd(degrees(W)))]
end

position_cox(W)=position_regular_class(W,maximum(degrees(W)))

# charvalues of H on word w living in some parabolic 
function charvalues_parabolic(H,w)
  W=H.W
  pw=parabolic_closure(W,sort(unique(w)))
  WI=reflection_subgroup(W,pw)
  if WI==W return end
  println("w=",joindigits(w))
  p=inclusiongens(WI,W)
  if !all(x->x in p,w) return end
  q=prod(x->prod(x)^0,H.para)
  cc=char_values(hecke(WI,H.para[simple_reps(W)[p]]),restriction(WI)[w])
  if nothing in cc return end
  induction_table(WI,W).scalar*(cc*q)
end

# find the specialization of Mvp parameters leading to group algebra
function group_specialization(H)
  res=vcat(unique(vcat(map(H.para)do p
    e=length(p)
    map(1:e)do i
      if p[i] isa Mvp
        v=scalar(p[i])
        if !isnothing(v)
          if v!=E(e)^(i-1) error(v," is not ",E(e)^(i-1)) end
          return nothing
        end
        if length(p[i])!=1 || length(first(term(p[i],1)))!=1
          error(p[i]," is not a power of a variable")
        end
        m=first(term(p[i],1))
        first(variables(m))=>coefftype(H)(root(E(e)^(i-1)*last(term(p[i],1)),first(powers(m))))
      else
        if degree(p[i])!=valuation(p[i])
          error(p[i]," is not a power of a variable")
        end
        if degree(p[i])==0 && p[i][0]!=E(e)^(i-1)
          error(p[i]," is not ",E(e)^(i-1))
        end
        d=degree(p[i]);c=p[i][d]
        root(E(e)^(i-1)//c,d)
      end
    end
  end))...)
  unique(filter(!isnothing,res))
end

function check_specialize(H;spec=group_specialization(H))
  ctW=CharTable(H.W).irr
  xprintln("using ",spec)
  ctH=CharTable(H).irr
  ctH=map(ctH)do p
    if ismissing(p) return p
    else return value(p,spec...)
    end
  end
  bad=filter(i->!ismissing(ctH[i]) && ctH[i]!=ctW[i],CartesianIndices(ctH))
  #length(bad)
end

# get central characters without having to know Chartable
# for algebras above group algebra
function HeckeCentralCharacters(H,z=gcd(degrees(H.W));spec=group_specialization(H))
  W=H.W
  v=map(m->root(m,z),central_monomials(H))
  v1=CharTable(W).irr[:,position_regular_class(W,z)]
  v2=CharTable(W).irr[:,position_class(W,one(W))]
  map((x,y,z)->x*y//z//value(x,spec...),v,v1,v2)
end

# cutz(l,z) for word l returns a,l' such that l'=replace a times z->[] in l
function cutz(l,z)
  l1=Chevie.Replace(l,z,Int[])
  div(length(l)-length(l1),length(z)),l1
end

# encode class txt with classinfo given ie replace z,c,etc...
function EncodeClass(txt,info)
  dic=joindigits.(info.classtext).=>info.classnames
  dic=filter(x->length(x[2])==1 && length(x[1])>1,dic)
  sort!(dic,by=x->-length(x[1]))
  txt=joindigits(txt)
  for d in dic txt=replace(txt,d) end
  txt
end

# shows how to choose class representatives cunningly w.r.t. parabolics
function checkparabolic(W)
  pz(nbz,noz)=joindigits(noz)*"z"^nbz
  ci=classinfo(W)
  ct=ci.classtext
  zt=centerword(W)
  cc=ct[position_cox(W)]
  if !isempty(cutz(zt,cc)[2]) error("z=$zt is not power of c=$cc\n") end
  println("c=",joindigits(cc)," z=c^",cutz(zt,cc)[1])
  cl=classreps(W)
  z=cl[position_class(W,W(zt...))]
  n=ngens(W)
  res=[];triples=[]
  l=combinations(1:ngens(W))
  if semisimplerank(W)<ngens(W)
    ls=map(x->length(parabolic_closure(W,x)),l)
    lI=l[filter(i->ls[i]==length(l[i]),1:length(l))]
    ps=w->ls[findfirst(==(sort(unique(w))),l)]
  else lI=l;ps=w->length(unique(w))
  end
  lI=filter(x->count(j->issubset(x,j),lI)==2,lI)
  reached=Int[]
  for I in lI
    WI=reflection_subgroup(W,I)
    if issubset(inclusiongens(WI),1:ngens(W))
      p=fusion_conjugacy_classes(WI,W)
      zp=copy(p)
      for (pj,j) in enumerate(p)
        s=classinfo(WI).classtext[pj]
        for k in 0:order(z)-1
          a=position_class(W,cl[j]*z^k);push!(zp,a)
          nbz,noz=cutz(ct[a],zt)
  	  if !issubset(unique(noz),I)&&(length(s)<length(noz)||ps(s)<ps(noz))
            println("#$a ",pz(nbz,noz)," =>",pz(k,s)," ps=",ps(noz),"=>",ps(s))
            push!(triples,(no=a,name=EncodeClass(ct[a],ci),pz=pz(k,s)));
          end
        end
      end
      xprintln("parabolic:",WI," gives classes ",sort(unique(zp)))
      reached=union(reached,zp)
    end
  end
  nums=Int[]
  function f(i)local res
    res="class "*lpad(i,3)*":"*EncodeClass(ct[i],ci)
    if EncodeClass(ct[i],ci)!=ci.classnames[i]
      res="("*ci.classnames[i]*")" end
    res
  end
  for i in 1:length(ct)
    nbz,noz=cutz(ct[i],zt)
    if ps(noz)<n
      push!(nums,i)
      if iszero(ps(noz)) println(f(i))
      else println(string(f(i)," in ",pz(nbz,Int[]),"W_",pz(0,sort(unique(noz)))))
      end
    else println(string(f(i))," ***BAD")
    end
  end
  println("bad classes:",setdiff(1:length(ct),nums))
  println("unreached:",setdiff(1:length(ct),reached))
  triples
end

#showregular:=function(W)local tbl,txt,d,cl,q,already,i,cc,lpi,lb,ci;
#  ci:=ChevieClassInfo(W);
#  txt:=ci.classtext;
#  cc:=txt[position_cox(W)];Print("c=",IntListToString(cc),"\n");
#  cl:=List(ConjugacyClasses(W),Representative);
#  lpi:=Sum(ReflectionDegrees(W)+ReflectionCoDegrees(W));
#  already:=[]; tbl:=[];lb:=[];
#  for d in Reversed(RegularEigenvalues(W)) do
#    q:=[];
#    for i in PrimeResidues(d) do
#      q:=PositionRegularClass(W,i/d);
#      if not q in already then 
#        Add(already,q);
#        Add(lb,SPrint(i/d));
#        Add(tbl,[q,lpi*i/d,Length(txt[q]),EncodeClass(txt[q],ci),
#	  ChevieClassInfo(W).classnames[q]]);
#      fi;
#    od;
#  od;
#  Print(FormatTable(tbl,
#    rec(rowLabels:=lb,
#        columnLabels:=["class","th. lg","lg","txt","name"],
#	rowsLabel:="RegEig")));
#end;

# partial chartable of Hecke algebra H on classtexts n some W_I<z>
function fromparabolic(H)
  W=H.W 
  ct=classinfo(W).classtext
  n=length(ct)
  t=convert(Matrix{Union{Missing,coefftype(H)}},fill(missing,n,n))
  zt=centerword(W)
  v=HeckeCentralCharacters(H)
  for i in 1:n
    nbz,noz=cutz(ct[i],zt)
    println("class $i:z^$nbz.",joindigits(noz))
    if length(unique(noz))<ngens(W)
      c=charvalues_parabolic(H,noz)
      if !isnothing(c)
        xprintln(reflection_subgroup(W,unique(noz)))
        t[:,i]=map((x,y)->x*y^nbz,c,v)
      else @show noz
      end
    end
  end
  t
end

# partial chartable of Hecke algebra H on classtexts in W_I<z>
function fromHI(H,I)
  W=H.W
  reps=classinfo(W).classtext
  n=length(reps)
  zt=centerword(W)
  v=HeckeCentralCharacters(H)
  t=convert(Matrix{Union{Missing,eltype(v)}},fill(missing,n,n))
  R=reflection_subgroup(W,I)
  rgens=inclusiongens(R)
  HR=hecke(R,H.para[simple_reps(W)[rgens]])
  tbl=induction_table(R,W).scalar
  for i in 1:n
    nbz,noz=cutz(reps[i],zt)
    if issubset(unique(noz),I) && all(i->i in rgens,noz)
      vv=char_values(HR,restriction(R)[noz])
      if !(nothing in vv)
        c=tbl*vv
        xprintln("class ",i,":",R)
        t[1:n,i]=map((x,y)->x*y^nbz,c,v)
      end
    end
  end
  t
end

# merge chartable t with partial table vals
function mergepartial(t,vals)
  for i in axes(t,1), j in axes(t,2)
    if !ismissing(vals[i,j])
      if !ismissing(t[i,j]) && t[i,j]!=vals[i,j] 
        error("t[$i,$j]=",t[i,j]," but vals[$i,$j]=",vals[i,j])
      end
      t[i,j]=vals[i,j]
    end
  end
end

# give a primitive root of pi -- returns partial chartable for its powers
function fromrootpi(H,w;spec=group_specialization(H))
  W=H.W
  o=div(sum(degrees(W).+codegrees(W)),length(w))
  z=gcd(degrees(W))
  wz=centerword(W)
  if order(W(w...))!=o error(w," is not root of π") end
  pospoww=map(i->position_class(W,W(w...)^i),1:o-1)
  println("gives classes ",sort(unique(pospoww)))
  p=unique(map(i->findfirst(==(i),pospoww),pospoww))
  p=map(x->[x,pospoww[x]],p)
  ct=classinfo(W).classtext
  iword(i)=vcat(repeat(wz,div(z*i,o)),repeat(w,mod(i,div(o,z))))
  bad=filter(i->ct[p[i][2]]!=iword(p[i][1]),1:length(p))
  if !isempty(bad)
    for i in bad
      println("Warning!! ",p[i][2],"-th representative is ",
              joindigits(ct[p[i][2]])," instead of ",joindigits(iword(i)))
    end
  end
  n=length(ct)
  t=convert(Matrix{Union{Missing,coefftype(H)}},fill(missing,n,n))
  for x in p[setdiff(1:length(p),bad)]
    @show x
    t[:,x[2]]=map((m,c)->m^(x[1]//o)*c[x[2]]//
                       value(m^(x[1]//o),spec...),
       central_monomials(H),eachrow(CharTable(W).irr))
  end
  t
end

function from_representations(H,inds=1:nconjugacy_classes(H.W))
  W=H.W
  cl=classinfo(W).classtext
  n=length(cl)
  t=convert(Matrix{Union{Missing,coefftype(H)}},fill(missing,n,n))
  for i in inds
    r=improve_type(representation(H,i))
    if !isnothing(r)
      print(i," ")
      t[i,:]=traces_words_mats(r,cl)
    end
  end
  println()
  t
end

#myeig:=function(table,i,class)local v;
#  v:=Eigenvalues(table,table.irreducibles[i],class);
#  return Concatenation(List([1..Length(v)],i->List([1..v[i]],j->i/Length(v))));
#end;
#
#coxclasses:=W->Set(List(Arrangements(W.generators,Length(W.generators)),
#   x->position_class(W,Product(x))));
#
#refssmaller:=function(W,w)local refs, r;
#  refs:=Reflections(W);
#  r:=Set(List(refs,x->Position(refs,x)));
#  return Filtered(r,x->ReflectionLength(W,refs[x]^-1*w)<ReflectionLength(W,w));
#end;
#
#refl2:=function(W,w)local m; m:=MatXPerm(W,w); return RankMat(m-m^0);end;
#
#test:=function(W)local a,b;
#  Stime();
#  a:=Set(List(HurwitzOrbitItems(W.generators),x->Position(Reflections(W),x)));
#  Print("Hurwitz:",Stime(),"\n");
#  b:=refssmaller( W, Product( W.generators ) );
#  Print("smaller:",Stime(),"\n");
#  return a=b;
#end;
#
##Check eigenvalues of Coxeter elements
#CheckCox:=function(W)local h,c,cl,i,ct,refs,n,gens,tbl,ref;
#  PrintDiagram(W);
#  h:=Maximum(ReflectionDegrees(W));
#  c:=PositionRegularClass(W,h);
#  n:=Length(W.generators);
#  cl:=coxclasses(W);
#  if not c in cl then
#    Print("******* Warning *****: ",h,"-regular class not product of gens\n");
#  fi;
#  ct:=CharTable(W);
#  refs:=[1..Length(ct.irreducibles)];
#  gens:=List(W.generators,i->[position_class(W,i),W.rank-1+E(OrderPerm(i))]);
#  refs:=Filtered(refs,i->ct.irreducibles[i][1]=W.rank);
#  refs:=Filtered(refs,i->ForAll(gens,v->ct.irreducibles[i][v[1]]=v[2]));
#  ref:=ChevieCharInfo(W).extRefl[2];
#  Print("Chevie's reflection representation is number:",ref,"\n");
#  Print("Reflection reps where generators are distinguished:",refs,"\n");
#  Print("classes for products of generators:\n");
#  for i in cl do 
#    Print("class ",i,"\n",FormatTable(List(refs,j->myeig(ct,j,i)),
#      rec(rowLabels:=refs,
#      columnLabels:=
#       Concatenation(["rep","eigenvalues"],List([1..W.rank-1],x->"")))),"\n");
#  od;
#  Print("\nh=",h,"-regular class =",c,"=",ChevieClassInfo(W).classtext[c],"\n");
#end;
#
## find for each char. x of 1-param Hecke H an e such that x is in C[q^{1/e}]
#gete:=function(H)local x,ct,fp;
#  fp:=function(p)
#    if IsCyc(p) then return [1];fi;
#    p:=Filtered(p.elm,x->Length(x.elm)=1);
#    p:=List(p,x->Mod1(x.coeff[1]));
#    return Set(List(Union(Set(p),[0]),Denominator));
#  end;
#  ct:=CharTable(H);
#  if ct<>false then 
#    ct:=ct.irreducibles; return List(ct,x->Lcm(Union(List(x,fp))));
#  else
#    return List(HeckeCentralCharacters(H),x->Lcm(fp(x)));
#  fi;
#end;
#
## check "z/e": for each family if characters of Spetsial algebra in family
## split over q^{1/e} then a+A is divisible by |ZW|/e
#checkze:=function(W)local uc,e,z,f,aA,fe;
#  uc:=UnipotentCharacters(W);
#  e:=ChevieCharInfo(W).hgal;
#  z:=OrderCenter(W.type[1]);
#  Print("   z=",z,"\n");
#  for f in uc.families do
#    aA:=uc.a[f.charNumbers[1]]+uc.A[f.charNumbers[1]];
#    fe:=List(f.charNumbers,i->Position(uc.harishChandra[1].charNumbers,i));
#    fe:=Filtered(fe,x->x<>false);
#    fe:=Lcm(List(Orbits(Group(e),fe)));
#    Print([aA,fe,z/Gcd(aA,z),aA/(z/fe)],"\n");
#  od;
#end;

exponentdenominator(p::Mvp)=lcm(map(m->lcm(collect(denominator.(powers(m)))),monomials(p)))

exponentdenominator2(a::AbstractArray{<:Mvp})=
 gcd(map(p->gcd(map(m->gcd(collect(numerator.(powers(m)))),monomials(p))),a))//
 lcm(map(p->lcm(map(m->lcm(collect(denominator.(powers(m)))),monomials(p))),a))

# eigenvalues of π^r π_I^r1
function eigenrootspi(H,r,I,r1;spec=group_specialization(H))
  W=H.W
  WI=reflection_subgroup(W,I)
  HI=hecke(WI,H.para[simple_reps(W)[inclusiongens(WI)]])
  e=HeckeCentralCharacters(H;spec).^r
  t=Integer.(induction_table(WI,W).scalar)
  eI=HeckeCentralCharacters(HI,order(W(I...));spec).^r1
  res=map((c,ei)->vcat(map((i,m)->fill(m*ei,i),c,eI)...),eachrow(t),e)
  rat=map(c->lcm(exponentdenominator.(c)),eachrow(CharTable(H).irr))
  res=map((v,r)->filter(x->isinteger(r//exponentdenominator(x)),v),res,rat)
  fields=map(c->NF(union(coefficients.(c)...)),eachrow(CharTable(H).irr))
  res=map((v,F)->filter(x->all(y->y in F,coefficients(x)),v),res,fields)
  tally.(res)
end

function hard32(H,i;spec=group_specialization(H))
  if i==30 l=eigenrootspi(H,1//2,[1,4],1;spec)
  elseif i==41 l=eigenrootspi(H,1//2,[1],1;spec)
  elseif i==58 l=eigenrootspi(H,3//2,[1],1;spec)
  end
  t=from_representations(H)[:,i]
  ct=CharTable(H.W).irr[:,i]
  for j in eachindex(l)
    if length(l[j])==0 
      if ismissing(t[j]) t[j]=0
      elseif t[j]!=0 error("l[$j]=",l[j]," but t[$j]=",t[j])
      end
    elseif length(l[j])==1 
#    xprintln("l[$j]=",l[j])
      v=value(l[j][1][1],spec...)
      ratio=ct[j]//v
      if isinteger(ratio) && abs(Int(scalar(ratio)))<=l[j][1][2]
        if ismissing(t[j]) t[j]=l[j][1][1]*ratio
           error("l[$j]=",l[j][1][1]*ratio ," but t[$j]=",t[j])
        end
      else error("l[$j]=",l[j]," but ct[$j,$i]=",ct[j])
      end
    else v=0
      for (m,c) in l[j] 
        v+=Mvp(Symbol("x_",j,"_",c))*m
      end
      if ismissing(t[j]) t[j]=v end
    end
  end
  t
end
  
function galv(vv)
  u,v,w=Mvp(:u,:v,:w)
  l=[value(vv[32],:u=>w,:v=>u,:w=>v)-vv[47]]
  push!(l,value(vv[32],:u=>v,:v=>u)-vv[48])
  push!(l,value(vv[49],:u=>w,:v=>u,:w=>v)-vv[50])
  push!(l,value(vv[49],:u=>v,:v=>u)-vv[51])
  push!(l,value(vv[52],:u=>w,:v=>u,:w=>v)-vv[59])
  push!(l,value(vv[52],:u=>v,:v=>u)-vv[60])
  push!(l,value(vv[53],:u=>w,:v=>u,:w=>v)-vv[54])
  push!(l,value(vv[53],:u=>v,:v=>w,:w=>u)-vv[55])
  push!(l,value(vv[53],:v=>w,:w=>v)-vv[56])
  push!(l,value(vv[53],:u=>v,:v=>u)-vv[57])
  push!(l,value(vv[53],:u=>w,:w=>u)-vv[58])
  push!(l,value(vv[67],:u=>w,:v=>u,:w=>v)-vv[68])
  push!(l,value(vv[67],:u=>v,:v=>w,:w=>u)-vv[69])
  push!(l,value(vv[67],:v=>w,:w=>v)-vv[70])
  push!(l,value(vv[67],:u=>v,:v=>u)-vv[71])
  push!(l,value(vv[67],:u=>w,:w=>u)-vv[72])
  push!(l,value(vv[73],:u=>u,:v=>w,:w=>v)-vv[74])
  push!(l,value(vv[73],:u=>v,:v=>w,:w=>u)-vv[75])
  push!(l,value(vv[77],:v=>w,:w=>v)-vv[76])
  push!(l,value(vv[77],:u=>v,:v=>w,:w=>u)-vv[78])
  push!(l,vv[73]-conj.(vv[77]))
  push!(l,value(vv[80],:u=>w,:v=>u,:w=>v)-vv[89])
  push!(l,value(vv[80],:u=>v,:v=>u)-vv[90])
  push!(l,vv[101]-conj.(vv[102]))
  l
end

function schurrel(H,v;p=Perm())
  spec=[:u=>Mvp(:x),:v=>E(3),:w=>E(3,2)]
  v=map(x->value(x,spec...),v)
  s=CycPol.(invpermute(schur_elements(hecke(H.W,Mvp(:x))),p))
  ll=lcm(s)
  s=ll.//s
  cs=map(x->x(Mvp(:x)),s)
  sum(v.*cs)
end

# return Lcm(schur).//s
function schur2(H;c=1)
  un=prod(vcat(map(x->one.(x),H.para)...))
  if un isa Mvp
    s=factorized_schur_elements(H)
    Lcm=lcm(s...)
    s=Ref(Lcm).//s
    print("expanding lcm(Sᵪ)/Sᵪ quotients..")
    t=@elapsed s=HeckeAlgebras.expand.(s;c).*un
    Lcm=HeckeAlgebras.expand(Lcm;c)
    println("done in ",t)
  else
    s=schur_elements(H)
    Lcm=s[charinfo(H.W).positionId]
    s=Lcm./s
  end
  s
end

function ctH32()
  W=crg(32)
  vars=Mvp.([:u,:v,:w])
  H=hecke(W,[Mvp{Cyc{Rational{Int64}}, Rational{Int64}}.(vars)])
  ct=CharTable(W).irr
  n=nconjugacy_classes(W)
  t=convert(Matrix{Union{Missing,coefftype(H)}},fill(missing,n,n))
  mergepartial(t,fromparabolic(H))
  mergepartial(t,fromrootpi(H,[4,3,2,1]))
  mergepartial(t,fromrootpi(H,[1,2,4,3,2]))
  mergepartial(t,from_representations(H))
  t
end

function meminfo_julia()
  pr(m)=rpad(string(round(m;digits=3)),5,'0')
  # @printf "GC total:  %9.3f MiB\n" Base.gc_total_bytes(Base.gc_num())/2^20
  # Total bytes (above) usually underreports, thus I suggest using live bytes (below)
  println("GC live: ", pr(Base.gc_live_bytes()/2^20),"MB")
  println("JIT:     ", pr(Base.jit_total_bytes()/2^20),"MB")
  println("Max. RSS:", pr(Sys.maxrss()/2^20),"MB")
end

malleperm=perm"(1,2)(3,6,5,4)(7,8)(9,10)(13,14)(15,16)(17,18,19,20)(21,22,24,23)(27,30,29,28)(31,32)(33,34)(37,38)(39,40)(42,43)(44,45)(46,47)(50,51,52,53)(54,55)(56,57)"
