# Hand-translated part of chevie/tbl/cmpxtimp.g
# (C) jean Michel 2011-2017
#
# Data on spetses 3G333, 4G333 and 3G422

#Group( ( 1, 2,44)( 4, 8,49)( 5,48,24)( 6, 9,26)( 7,25,47)(10,11,50)(14,17,38)
#(15,36,43)(16,37,32)(18,19,39)(21,22,52)(41,42,54), ( 2,12,37)( 3,10,15)
#( 5,26,40)( 6,24,27)( 7,42,13)( 8,33,19)(14,21,46)(16,53,44)(18,51,49)
#(22,38,31)(28,54,47)(30,36,50)); # a section G of the Diagram automorphisms

# here W should have type[1] be G333 with indices [1,2,3]
# the representatives chosen are for the 24 elements of G: 15 in G and 9 in -G
chevieset(:timp, :ReducedInRightCoset, function (W, phi)
  m1perm=perm"(1,4)(2,8)(3,13)(5,22)(6,14)(7,10)(9,17)(11,25)(12,33)(15,42)
(16,18)(19,37)(20,35)(21,24)(23,29)(26,38)(27,46)(28,30)(31,40)(32,39)(34,45)
(36,54)(41,43)(44,49)(47,50)(48,52)(51,53)" # effect of -1 on roots
  m1perm=mappingPerm(inclusion(W),inclusion(W)^m1perm)
  g=PermX(W,reflrep(W,phi)*reflrep(W,m1perm))
  if g in W return (gen=inclusion(W,1:3),phi=inv(g)*phi) end
  # quads of roots which have the same CartanMat and are representatives of 
  # <-1,W>-orbits of quads of reflections satisfying Corran-Picantin relations
  sets=[[1,2,3,44],[2,12,11,37],[3,11,2,36],[1,12,10,16]] # 3G333
  for g in [Perm(),m1perm], i in inclusion.(Ref(W),sets), 
     o in [[4,1,3,2],[2,4,3,1]]
    e=transporting_elt(W,i[o],i.^(g*phi);action=(s,g)->s.^g)
    if !isnothing(e)return (gen=i[1:3],phi=phi/e) end
  end
  # sextuples of roots which have the same CartanMat and are representatives of 
  # W-orbits of each element of order 4 of G
  sets2=[[1,2,32,16,3,36,30,10],[3,10,30,36,17,21,53,38],
        [2,12,16,53,11,10,43,36],[2,44,16,37,3,43,30,11],
        [1,37,32,44,15,30,50,3],[1,12,32,53,10,50,36,15]] # 4G333
  for i in inclusion.(Ref(W),sets2)
    e=transporting_elt(W,i[[2,3,4,1]],i[1:4].^phi;action=(s,g)->s.^g)
    if !isnothing(e)return (gen=i[[1,2,5]],phi=phi/e) end
  end
end)
