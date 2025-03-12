# Hand-translated part of chevie/tbl/weyl2d.jl
# (C) Frank Luebeck 1994-2001
# Data for the coset of type 2D which can be seen as the non-trivial
# coset of W(D_r) inside W(B_r).

chevieset(Symbol("2D"), :ClassParams, function(n)
  B=chevieget(:B,:ClassParams)(n)
  filter(a->isodd(length(a[2])), B)
end)

chevieset(Symbol("2D"),:WordsClassRepresentatives, n->
  chevieget(Symbol("2D"),:ClassInfo)(n)[:classtext])

# We  parametrize the  F-conjugacy classes  by the  classes in  the coset
# Bn-Dn.  If n is  odd, since F  is inner acting  as w0, it would also be
# possible  to  parametrize  them  by  the  classes  in  Dn  (you get the
# F-classes by translating with w0). This gives two possible labelings of
# the  F-classes: one by the Dn-classes  and one by the outer Bn-classes.
# They  correspond as follows: Let [a,b] the double partition for w in Dn
# and let [c,d] be the double partition for w.w0Bn in Bn. Then c contains
# the  even entries of a and the odd entries of b and d contains the even
# entries of b and the odd entries of a.
chevieset(Symbol("2D"), :ClassInfo, function(n)
  B=chevieget(:B, :ClassInfo)(n)
  l=filter(i->isodd(length(B[:classparams][i][2])),1:length(B[:classtext]))
  Dict{Symbol, Any}(:classnames =>B[:classnames][l],
                    :classparams=>B[:classparams][l],
                    :classes=> B[:classes][l], 
                    :classtext => map(B[:classtext][l])do l
  # 1 is the automorphism, and u=121 is the new generator of W(2D). We deal
  # with words with an odd number of 1. To gather one 1 at right we go left
  # to right doing substitutions 11->, 12->u1, and 1a->a1 if a<>2.
                      res=Int[]
                      n=1
                      for i in l
                        if i==1 n=mod(n+1,2)
                        elseif i==2 push!(res,2-n)
                        else push!(res, i)
                        end
                      end
                      res
                  end)
end)

chevieset(Symbol("2D"), :NrConjugacyClasses, function(n)
  if isodd(n) div(npartition_tuples(n,2),2)
  else div(npartition_tuples(n,2)-npartitions(div(n,2)),2)
  end
end)

# test if a character of W(B) corresponds to the preferred extension
# for ^2D, see [CS,17.2] and [Lusztig-book,4.4,4.18]:
chevieset(Symbol("2D"), :IsPreferred, function(pp)
  pp=symbol_partition_tuple(pp,0)
  pp[1]>pp[2]
end)

chevieset(Symbol("2D"), :CharParams, n->
  filter(chevieget(Symbol("2D"),:IsPreferred),partition_tuples(n,2)))

#the map which goes from almost characters to unipotent characters for 2Dn
function Defect0to2(ST)
  a=minimum(symdiff(ST[1], ST[2]))
  ST=sort.([symdiff(ST[1], [a]), symdiff(ST[2], [a])])
  if length(ST[1])>length(ST[2]) ST
  else reverse(ST)
  end
end

chevieset(Symbol("2D"), :CharInfo, function (n)
  res=Dict{Symbol, Any}(:charparams=>chevieget(Symbol("2D"),:CharParams)(n))
  res[:extRefl]=map(i->[fill(1,i),[n-i]],0:n-2)
  append!(res[:extRefl],[[[1],fill(1,n-1)],[Int[],fill(1,n)]])
  res[:extRefl]=map(x->findfirst(y->y==x || y==reverse(x),res[:charparams]),
                    res[:extRefl])
  resparams=chevieget(:imp,:CharParams)(2,2,n)
  res[:charRestrictions]=map(x->findfirst(y->y==x || y==reverse(x),resparams),
                             res[:charparams])
  res[:nrGroupClasses]=length(resparams)
  res[:charnames]=string_partition_tuple.(res[:charparams])
  f=map(c->fegsymbol(symbol_partition_tuple(c,0),1),res[:charparams])
  res[:b]=valuation.(f)
  res[:B]=degree.(f)
  res[:charSymbols]=map(c->Defect0to2(symbol_partition_tuple(c,0)),
                        res[:charparams])
  res[:a]=valuation_gendeg_symbol.(res[:charSymbols])
  res[:A]=degree_gendeg_symbol.(res[:charSymbols])
  res
end)

chevieset(Symbol("2D"),:FakeDegree,(n,c,q)->
  fegsymbol(symbol_partition_tuple(c,0),1)(q))

chevieset(Symbol("2D"),:PhiFactors,n->vcat(fill(1,n-1),-1))

# This function returns the part of the character table the Coxeter group
# of  type B_l on classes  outside a reflection subgroup  of type D_l for
# the characters which remain irreducible on restriction to this subgroup
# and which correspond to the *preferred* extensions defined in [CS,17.2,
# case D_l].
#
# Alternatively  you can get the  *good* extension instead of *preferred*
# extension by defining testchar appropriately.
chevieset(Symbol("2D"),:HeckeCharTable,
function (l,param,sqrtpara)
  q=-param[1][1]//param[1][2]
  q=vcat([[q^0,-1]],map(i->[q,-1],2:l))
  hi=chevieget(:B,:HeckeCharTable)(l,q,Int[])
  chr=1:length(hi[:classparams])
  lst=filter(i->isodd(length(hi[:classparams][i][2])),chr);
  tbl=Dict{Symbol,Any}(:identifier=>"H(^2D$l)",
                       :size=>div(hi[:size],2),
                       :orders=>hi[:orders][lst],
                       :centralizers=>div.(hi[:centralizers][lst],2),
                       :classes=>hi[:classes][lst],
	   :text=>"extracted from generic character table of HeckeB")
  merge!(tbl,chevieget(Symbol("2D"),:ClassInfo)(l))
  para=chevieget(Symbol("2D"),:CharParams)(l)
  tbl[:irredinfo]=map(i->Dict{Symbol,Any}(:charparam=>para[i],
      :charname=>string_partition_tuple(para[i])),eachindex(para))
  para=partition_tuples(l,2)
  chr=filter(i->chevieget(Symbol("2D"),:IsPreferred)(para[i]),chr)
  T=Tbasis(hecke(coxgroup(:B,l),q))
  tbl[:irreducibles]=transpose(toM(map(x->char_values(
    T(1,Replace(x,[1],[1,2,1])...),hi[:irreducibles][chr,:]),
                                       tbl[:classtext])))
  AdjustHeckeCharTable(tbl,param)
end)

chevieset(Symbol("2D"),:CharTable,l->chevieget(Symbol("2D"),:HeckeCharTable)(l,fill([1,-1],l),fill(1,l)))

chevieset(Symbol("2D"), :UnipotentCharacters, function (rank,)
  uc=Dict{Symbol,Any}(:harishChandra=>[],:charSymbols=>[],:almostHarishChandra=>[])
  for d in (2:4div(isqrt(rank)-1,2)+2)
    r=div(d^2,4)
    s=Dict{Symbol, Any}(:relativeType=>
     TypeIrred(;series=:B,indices=1+r:rank,rank=rank-r),:levi=>1:r,
      :eigenvalue=>1,# see Geck-malle
      :parameterExponents=>vcat(d,fill(1,max(0,rank-1-r))))
    s[:cuspidalName]="{}^2D"*stringind(rio(TeX=true),r)
    push!(uc[:harishChandra],s)
    if d==2
      s[:levi]=Int[]
      s[:cuspidalName]=""
    end
    symbols=BDSymbols(rank,d)
    s[:charNumbers] = (1:length(symbols)).+length(uc[:charSymbols])
    FixRelativeType(s)
    append!(uc[:charSymbols], symbols)
  end
  uc[:a]=valuation_gendeg_symbol.(uc[:charSymbols])
  uc[:A]=degree_gendeg_symbol.(uc[:charSymbols])
  uc[:almostCharSymbols]=map(i->[[0],[0]],1:sum(x->length(x[:charNumbers]),
                                                uc[:harishChandra]))
  for d in  4*(0:isqrt(div(rank,4)))
    r=div(d^2,4)
    s=Dict{Symbol, Any}(:relativeType=>
      TypeIrred(;series=:B,indices=1+r:rank,rank=rank-r),:levi=>1:r,
      :eigenvalue=>(-1)^div(d+1,4))
    s[:cuspidalName]="D"*stringind(rio(TeX=true),r)
    r=s[:relativeType][:rank]
    symbols=BDSymbols(rank,d)
    if isodd(div(d+1,4)) symbols=reverse.(symbols) end
    if d==0
      s[:relativeType].series=:D
      s[:relativeType]=TypeIrred(;orbit=[s[:relativeType]],twist=perm"(1,2)")
      s[:cuspidalName]= ""
      symbols=map(x->symbol_partition_tuple(x, 0),
                  chevieget(Symbol("2D"),:CharParams)(rank))
    end
    s[:charNumbers]=map(s->findfirst(==(Defect0to2(s)),uc[:charSymbols]),symbols)
    uc[:almostCharSymbols][s[:charNumbers]]=symbols
    if d!=0 FixRelativeType(s) end
    push!(uc[:almostHarishChandra],s)
  end
  # note: delta is always 1 since a+A is always even
  z(x)=(Z1=sort(symdiff(x[1],x[2])),Z2=intersect(x...))
  uc[:families]=map(sort(unique(z.(uc[:charSymbols]))))do f
    sharp(s)=symdiff(setdiff(s[2],f.Z2),f.Z1[1:2:length(f.Z1)-1])
    res=Dict{Symbol,Any}(:charNumbers=>filter(i->z(uc[:charSymbols][i])==f,
                                                 1:length(uc[:charSymbols])))
    res[:almostCharNumbers]=res[:charNumbers]
    res[:fourierMat]=map(u->map(a->(1//2)^div(length(f.Z1)-1,2)*
     (-1)^length(intersect(sharp(u),sharp(a))),
     uc[:almostCharSymbols][res[:almostCharNumbers]]),
                           uc[:charSymbols][res[:charNumbers]])
    if length(res[:fourierMat])==16 #JM jan 2015: fix this horrible kludge
      res[:fourierMat][16]=-res[:fourierMat][16]
      res[:fourierMat][1:16][16]=-res[:fourierMat][1:16][16]
    end
    res[:eigenvalues]=fill(1,length(res[:charNumbers]))
    res[:sh]=fill(1,length(res[:charNumbers]))# is that correct for Geck-Malle?
    if length(res[:eigenvalues])==1
      res[:charLabels]=[""]
      res[:special]=1
    else
      res[:charLabels]=map(uc[:charSymbols][res[:charNumbers]])do M
        M=symdiff(setdiff(M[2],f.Z2),f.Z1[3:2:length(f.Z1)-1])
        v=map(z->mod(count(>=(z),M),2),f.Z1)
        D=length(v)
        v1=v[2:2:D-mod(D,2)]
        v2=v[3:2:D-1+mod(D,2)]
        if isodd(D) push!(v1, 0) end
        # v1, v2 is coordinates in (e1,e3,e5,..) and in (e2,e4,..) basis
        v1=map(i->mod(sum(v1[i:i+1]),2),1:length(v2))
        # coordinates in e1, e1+e3, e1+e3+e5, ...
        s="+-"
        s[v2.+1]*","*s[v1.+1]
      end
    end
    res[:special]=findfirst(x->all(y->y in "+,",x),res[:charLabels])
    res[:name]=vcat(f.Z1,f.Z2,f.Z2)
    sort!(res[:name])
    res[:name]=joindigits(res[:name])
    res[:printname]=res[:name]
    res[:explanation]="classical family"
    res[:perm]=Perm()
    res[:size]=length(res[:charNumbers])
    res
  end
  uc
end)

chevieset(Symbol("2D"),:Ennola,function(n)
  if isodd(n) return SPerm() end
  uc=chevieget(Symbol("2D"),:UnipotentCharacters)(n)
  l=uc[:charSymbols]
  SPerm(map(1:length(l))do i
    if !(l[i][2] isa AbstractVector) return i*(-1)^uc[:A][i] end
    s=EnnolaSymbol(l[i])
    p=findfirst(==(s),l)
    if isnothing(p) p=findfirst(==(reverse(s)),l) end
    p*(-1)^uc[:A][i]
  end)
end)

chevieset(Symbol("2D"), :UnipotentClasses, function (r, p)
  uc=deepcopy(chevieget(:D, :UnipotentClasses)(r, p))
  if p==2 return uc end
  for cc in uc[:classes]
    cc[:red]=coxgroup()
    for j in tally(cc[:parameter])
      d=div(j[2],2)
      if iseven(j[1]) cc[:red]=Cosets.extprod(cc[:red],coxgroup(:C,d))
      elseif isodd(j[2])
        if j[2]>1 cc[:red]=Cosets.extprod(cc[:red],coxgroup(:B,d)) end
      elseif j[2]>2
        cc[:red]=Cosets.extprod(cc[:red],d==3 ? rootdatum(:psu,4) :
                                              rootdatum("pso-",2*d))
      else cc[:red]=Cosets.extprod(cc[:red],torus([-1;;]))
      end
    end
  end
  uc
end)

chevieset(Symbol("2D"), :ClassParameter, function (n, w)
  x=prod(i->i==1 ? SPerm(1,-2) : SPerm(i-1,i),w;init=SPerm())
  cycletype(x*SPerm(1,-1),n)
end)

chevieset(Symbol("2D"), :HeckeRepresentation, function(n,para,sqpara,i)
   param=chevieget(Symbol("2D"),:CharInfo)(n)[:charparams][i]
   bno=findfirst(==(param),chevieget(:B,:CharInfo)(n)[:charparams])
   parab=copy(para);parab[1]=[1,-1]
   r=chevieget(:B,:HeckeRepresentation)(n,parab,[],bno)
   if all(x->x==[1,-1],para) u=r[1] else u=inv(r[1]*1//1) end
   (gens=pushfirst!(r[2:end],r[1]*r[2]*u),F=r[1])
end)

chevieset(Symbol("2D"), :Representation, function(n,i)
  return chevieget(Symbol("2D"),:HeckeRepresentation)(n,fill([1,-1],4),[],i)
end)
