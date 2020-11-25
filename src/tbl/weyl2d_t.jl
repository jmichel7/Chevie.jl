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
      :charname=>chevieget(Symbol("2D"),:CharName)(l,para[i])),eachindex(para))
  para=chevieget(:B,:CharParams)(l)
  chr=filter(i->chevieget(Symbol("2D"),:testchar)(para[i]),chr)
  tbl[:irreducibles]=toL(permutedims(toM(map(x->char_values(
   Tbasis(hecke(coxgroup(:B,l),q))(vcat([1],Replace(x,[1],[1,2,1]))...),
   toM(hi[:irreducibles][chr])),tbl[:classtext]))))
   CHEVIE[:compat][:AdjustHeckeCharTable](tbl,param)
  tbl
end)

chevieset(Symbol("2D"),:CharTable,l->chevieget(Symbol("2D"),:HeckeCharTable)(l,fill([1,-1],l),fill(1,l)))
