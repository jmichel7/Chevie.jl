# Hand-translated part of chevie/tbl/weyl2d.jl
# (C) Frank Luebeck 1994-2001
# Data for the coset of type 2D which can be seen as the non-trivial
# coset of W(D_r) inside W(B_r).

# This function returns the part of the character table the Coxeter group
# of  type B_l on classes  outside a reflection subgroup  of type D_l for
# the characters which remain irreducible on restriction to this subgroup
# and which correspond to the *preferred* extensions defined in [CS,17.2,
# case D_l].
#
# Alternatively  you can get the  *good* extension instead of *preferred*
# extension by defining testchar appropriately.
chevieset(Symbol("2D"),:HeckeCharTable,
function (l,param,rootparam)
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
      :charname=>PartitionTupleToString(para[i])),eachindex(para))
  para=chevieget(:B,:CharParams)(l)
  chr=filter(i->chevieget(Symbol("2D"),:IsPreferred)(para[i]),chr)
  tbl[:irreducibles]=toL(transpose(toM(map(x->char_values(
   Tbasis(hecke(coxgroup(:B,l),q))(vcat([1],GAPENV.Replace(x,[1],[1,2,1]))...),
   toM(hi[:irreducibles][chr])),tbl[:classtext]))))
   CHEVIE[:compat][:AdjustHeckeCharTable](tbl,param)
  tbl
end)

chevieset(Symbol("2D"),:CharTable,l->chevieget(Symbol("2D"),:HeckeCharTable)(l,fill([1,-1],l),fill(1,l)))

chevieset(Symbol("2D"), :UnipotentClasses, function (r, p)
  uc=deepcopy(chevieget(:D, :UnipotentClasses)(r, p))
  if p==2 return uc end
  for cc in uc[:classes]
    cc[:red]=coxgroup()
    for j in tally(cc[:parameter])
      d=div(j[2],2)
      if iseven(j[1]) cc[:red]=cc[:red]*coxgroup(:C,d)
      elseif isodd(j[2])
        if j[2]>1 cc[:red]=Cosets.extprod(cc[:red],coxgroup(:B,d)) end
      elseif j[2]>2
        cc[:red]=Cosets.extprod(cc[:red],d==3 ? rootdatum(:psu,4) :
                                              rootdatum("pso-",2*d))
      else cc[:red]=Cosets.extprod(cc[:red],torus([[-1]]))
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
   r=toM.(chevieget(:B,:HeckeRepresentation)(n,parab,[],bno))
   if all(x->x==[1,-1],para) u=r[1] else u=inv(r[1]*1//1) end
   (gens=pushfirst!(r[2:end],r[1]*r[2]*u),F=r[1])
end)

chevieset(Symbol("2D"), :Representation, function(n,i)
  return chevieget(Symbol("2D"),:HeckeRepresentation)(n,fill([1,-1],4),[],i)
end)
