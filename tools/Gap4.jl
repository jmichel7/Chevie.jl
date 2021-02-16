module Gap4

using GAP, Gapjm, Reexport
@reexport using GAP

GAP.julia_to_gap(p::Perm)=GAP.Globals.PermList(GAP.julia_to_gap(Int.(vec(p))))

GAP.julia_to_gap(g::Group)=isempty(gens(g)) ?
  GAP.Globals.Group(GAP.julia_to_gap(one(g))) :
  GAP.Globals.Group(GAP.julia_to_gap.(gens(g))...)

gapvec_to_julia(gv)=map(i->GAP.gap_to_julia(gv[i]),1:length(gv))
gapvec_to_julia(f,gv)=map(i->f(gv[i]),1:length(gv))

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
  gapvec_to_julia((Perm_to_julia∘GAP.Globals.Representative),
                   GAP.Globals.ConjugacyClasses(gg))
end

function CharTable_to_julia(ct)
  u=GAP.Globals.Irr(ct)
  girr=gapvec_to_julia(x->gapvec_to_julia(Cyc_to_julia,x),u) 
  irr=permutedims(hcat(girr...))
  n=gapvec_to_julia(GAP.Globals.ClassNames(ct))
  cn=gapvec_to_julia(GAP.Globals.CharacterNames(ct))
  sz=gapvec_to_julia(GAP.Globals.SizesCentralisers(ct))
  s=GAP.gap_to_julia(GAP.Globals.Size(ct))
  id=GAP.gap_to_julia(GAP.Globals.Identifier(ct))
  Chars.CharTable(irr,cn,n,sz,s,Dict{Symbol,Any}(:name=>id))
end

function CharTable(g::Group)
  gg=GAP.julia_to_gap(g)
  ct=CharTable_to_julia(GAP.Globals.CharacterTable(gg))
  l=map(x->position_class(g,x),Gap4.class_reps(g))
  ct.irr[:,l]=ct.irr
  ct.classnames[l]=ct.classnames
  ct
end
Chars.CharTable(g::CoxSym)=CharTable(g.G)

end
