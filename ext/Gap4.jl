module Gap4

#println("loading Gap4 extension to Chevie")

using GAP, Chevie

#@test perm=GAP.evalstr("(1,2,3)(4,5)");perm==Gap.obj(Perm(perm))

mapgap(f,gv)=map(f,GAP.gap_to_julia(gv;recursive=false))
fromgap(g)=GAP.gap_to_julia(g)

function CyclotomicNumbers.Cyc(p::GapObj)
  l=improve_type(Vector{Any}(GAP.Globals.ExtRepOfObj(p)))
  sum(i->E(length(l),i-1)*l[i],eachindex(l))
end

function gapclassreps(g::PermGroup)
  gg=GAP.Obj(g)
  mapgap(Perm∘GAP.Globals.Representative,GAP.Globals.ConjugacyClasses(gg))
end

function Chars.CharTable(ct::GapObj)
  u=GAP.Globals.Irr(ct)
  girr=mapgap(x->mapgap(Cyc,x),u) 
  irr=permutedims(hcat(girr...))
  n=fromgap.(GAP.Globals.ClassNames(ct))
  cn=fromgap.(GAP.Globals.CharacterNames(ct))
  sz=identity.(fromgap(GAP.Globals.SizesCentralisers(ct)))
  s=fromgap(GAP.Globals.Size(ct))
  id=fromgap(GAP.Globals.Identifier(ct))
  orders=fromgap(GAP.Globals.OrdersClassRepresentatives(ct))
  pp=keys(factor(s))
  powermaps=Vector{Any}(fill(nothing,maximum(pp;init=0)))
  for k in pp
    powermaps[k]=fromgap(GAP.Globals.PowerMap(ct,k))
  end
  Chars.CharTable(irr,cn,n,sz,s,Dict{Symbol,Any}(:name=>id,:orders=>orders,
                                                :powermaps=>powermaps))
end

function Chars.CharTable(g::PermGroup)
  gg=GAP.Obj(g)
  gc=gapclassreps(g) # in case :classreps did not use gapclassreps
  l=map(x->position_class(g,x),gc)
  ct=CharTable(GAP.Globals.CharacterTable(gg))
  ct.irr[:,l]=ct.irr
  ct.classnames[l]=ct.classnames
  ct.orders[l]=ct.orders
  for i in eachindex(ct.powermaps)
    if ct.powermaps[i]!=nothing
      ct.powermaps[i][l]=ct.powermaps[i]
    end
  end
  ct
end

end
