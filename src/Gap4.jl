module Gap4

println("loading Gap4")

using GAP, Gapjm, Reexport
#@reexport using GAP

GAP.julia_to_gap(p::Perm)=GAP.Globals.PermList(GAP.julia_to_gap(Int.(vec(p))))

function GAP.julia_to_gap(g::Group)
  get!(g,:GAP)do
    isempty(gens(g)) ?
    GAP.Globals.Group(GAP.julia_to_gap(one(g))) :
    GAP.Globals.Group(GAP.julia_to_gap.(gens(g))...)
  end
end

mapgap(f,gv)=map(f,GAP.gap_to_julia(gv;recursive=false))
fromgap(g)=improve_type(GAP.gap_to_julia(g))

Perms.Perm(p::GapObj)=Perm(mapgap(Int16,GAP.Globals.ListPerm(p)))

function CyclotomicNumbers.Cyc(p::GapObj)
  l=GAP.gap_to_julia(GAP.Globals.CoeffsCyc(p,GAP.Globals.Conductor(p)))
  sum(i->E(length(l),i-1)*l[i],eachindex(l))
end

function gapclassreps(g::PermGroup)
  gg=GAP.julia_to_gap(g)
  mapgap(Perm∘GAP.Globals.Representative,GAP.Globals.ConjugacyClasses(gg))
end

function Groups.classreps(g::Group)
  get!(g,:classreps)do
    if g isa PermGroup && length(g)>10000 gapclassreps(g)
    else first.(conjugacy_classes(g))
    end
  end
end

function CharTable(ct::GapObj)
  u=GAP.Globals.Irr(ct)
  girr=mapgap(x->mapgap(Cyc,x),u) 
  irr=permutedims(hcat(girr...))
  n=fromgap(GAP.Globals.ClassNames(ct))
  cn=fromgap(GAP.Globals.CharacterNames(ct))
  sz=fromgap(GAP.Globals.SizesCentralisers(ct))
  s=fromgap(GAP.Globals.Size(ct))
  id=fromgap(GAP.Globals.Identifier(ct))
  Chars.CharTable(irr,cn,n,sz,s,Dict{Symbol,Any}(:name=>id))
end

function CharTable(g::Group)
  gg=GAP.julia_to_gap(g)
  # in case :classreps not computed using gapclassreps
  l=map(x->position_class(g,x),gapclassreps(g))
  ct=CharTable(GAP.Globals.CharacterTable(gg))
  ct.irr[:,l]=ct.irr
  ct.classnames[l]=ct.classnames
  ct
end

function centralizer(g::PermGroup,p::Perm)
  gg=GAP.julia_to_gap(g)
  pp=GAP.julia_to_gap(p)
  c=GAP.Globals.Centralizer(gg,pp)
  Group(Perm.(GAP.Globals.GeneratorsOfGroup(c)))
end

function intersect(g1::PermGroup,g2::PermGroup)
  gg1=GAP.julia_to_gap(g1)
  gg2=GAP.julia_to_gap(g2)
  gg=GAP.Globals.Intersection(gg1,gg2)
  Group(Perm.(GAP.Globals.GeneratorsOfGroup(gg)))
end

end
using .Gap4
