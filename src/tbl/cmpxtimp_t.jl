# Hand-translated part of chevie/tbl/cmpxtimp.g
# (C) jean Michel 2011-2017
#
# Data on spetses 3G333, 4G333 and 3G422

chevieset(:timp, :CharName, function (p,q,r,phi,para,opt)
           CHEVIE[:imp][:CharName](p,q,r,convert.(Vector{Int},para),opt)
    end)

chevieset(:timp, :ReducedInRightCoset, function (W, phi)
  # quads of roots which have the same CartanMat and are representatives of 
  # W-orbits of quads of reflections satisfying Corran-Picantin relations
  sets=[[1,2,3,44],[21,3,1,32],[3,11,2,36],[22,3,2,16]]
  sets2=[[1,50,3,12,2],[3,52,2,23,11],[1,16,3,43,38],
         [2,37,3,15,14],[50,3,52,38,53],[1,23,3,22,45]]
  m1perm=perm"(1,4)(2,8)(3,13)(5,22)(6,14)(7,10)(9,17)(11,25)(12,33)(15,42)
(16,18)(19,37)(20,35)(21,24)(23,29)(26,38)(27,46)(28,30)(31,40)(32,39)(34,45)
(36,54)(41,43)(44,49)(47,50)(48,52)(51,53)" # effect of -1 on roots
  m1perm=mappingPerm(inclusion(W),inclusion(W)^m1perm)
  for g in [Perm(),m1perm], i in inclusion.(Ref(W),sets), 
     o in [[4,1,3,2],[2,4,3,1],1:4]
    e=transporting_elt(W,i[o],i.^(phi*g);action=(s,g)->s.^g)
    if !isnothing(e)return Dict{Symbol,Any}(:gen=>i[1:3],:phi=>phi/e) end
  end
  for i in inclusion.(Ref(W),sets2)
    e=transporting_elt(W,i[[2,3,4,1]],i[1:4].^phi;action=(s,g)->s.^g)
    if !isnothing(e)return Dict{Symbol,Any}(:gen=>i[[1,5,3]],:phi=>phi/e) end
  end
  return false
end)
