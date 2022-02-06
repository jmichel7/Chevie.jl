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
  chr=filter(i->chevieget(Symbol("2D"),:IsPreferred)(para[i]),chr)
  tbl[:irreducibles]=toL(transpose(toM(map(x->char_values(
   Tbasis(hecke(coxgroup(:B,l),q))(vcat([1],Replace(x,[1],[1,2,1]))...),
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
      if mod(j[1],2)==0 cc[:red]=cc[:red]*coxgroup(:C,d)
      elseif mod(j[2],2)!=0
        if j[2]>1 cc[:red]=Cosets.extprod(cc[:red],coxgroup(:B,d)) end
      elseif j[2]>2
        if d==3 cc[:red]=Cosets.extprod(cc[:red],spets(coxgroup(:A,d),perm"(1,3)"))
        else    cc[:red]=Cosets.extprod(cc[:red],spets(coxgroup(:D,d),perm"(1,2)"))
        end
      else cc[:red]=Cosets.extprod(cc[:red],torus([[-1]]))
      end
    end
  end
  uc
end)

chevieset(Symbol("2D"), :ClassParameter, function (n, w)
  x=Perm()
  for i in w
    x*= i==1 ? Perm(1, n+2)*Perm(2, n+1) : Perm(i-1,i)*Perm(i-1+n,i+n)
  end
  x*=Perm(1,n+1)
  res = [Int[],Int[]]
  mark=fill(true,n)
  for i in 1:n
    if mark[i]
      cyc=orbit(x, i)
      if i+n in cyc push!(res[2], div(length(cyc),2))
      else push!(res[1], length(cyc))
      end
      for j in cyc
        if j>n mark[j-n]=false
        else mark[j]=false
        end
      end
    end
  end
  sort!(res[1])
  sort!(res[2])
  [reverse(res[1]), reverse(res[2])]
end)
