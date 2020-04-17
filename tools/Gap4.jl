module Gap4

using GAP, Gapjm, Reexport
@reexport using GAP

GAP.julia_to_gap(p::Perm)=GAP.Globals.PermList(GAP.julia_to_gap(Int.(vec(p))))

GAP.julia_to_gap(g::Group)=GAP.Globals.Group(GAP.julia_to_gap.(gens(g))...)

gapvec_to_julia(gv)=map(i->GAP.gap_to_julia(gv[i]),1:length(gv))

function Perm_to_julia(p)
  l=GAP.Globals.ListPerm(p)
  length(l)==0 ? Perm() : Perm{Int16}(Int.(GAP.gap_to_julia(l)))
end

function Cyc_to_julia(p)
  l=GAP.gap_to_julia(GAP.Globals.CoeffsCyc(p,GAP.Globals.Conductor(p)))
  sum(i->E(length(l),i-1)*l[i],eachindex(l))
end

function class_reps(g::Group)
  gg=GAP.julia_to_gap(g)
  l=gapvec_to_julia(GAP.Globals.ConjugacyClasses(gg))
  Perm_to_julia.(GAP.Globals.Representative.(l))
end

function CharTable_to_julia(ct)
  girr=GAP.gap_to_julia(GAP.Globals.Irr(ct))
  irr=permutedims(hcat(map(x->Cyc_to_julia.(x),girr)...))
  n=gapvec_to_julia(GAP.Globals.ClassNames(ct))
  cn=gapvec_to_julia(GAP.Globals.CharacterNames(ct))
  sz=gapvec_to_julia(GAP.Globals.SizesCentralisers(ct))
  id=GAP.gap_to_julia(GAP.Globals.Identifier(ct))
  CharTable(irr,cn,n,sz,id,Dict{Symbol,Any}())
end

function Chars.CharTable(g::Group)
  gg=GAP.julia_to_gap(g)
  ct=CharTable_to_julia(GAP.Globals.CharacterTable(gg))
  l=map(x->position_class(g,x),Gap4.class_reps(g))
  ct.irr[:,l]=ct.irr
  ct.classnames[l]=ct.classnames
  ct
end
Chars.CharTable(g::CoxSym)=CharTable(g.G)

end
