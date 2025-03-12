# Hand-translated part of chevie/tbl/cmpxtimp.g
# (C) jean Michel 2011-2017
#
# Data on spetses 3G333, 4G333 and 3G422

chevieset(:timp, :PhiFactors, function (p, q, r, phi)
  o = order(phi)
  if p==q
    if mod(p, o)==0
      if phi==perm"(1,2,4)" return [1, 1, E(3, 2)]
      elseif phi^-1==perm"(1,2,4)" return [1, 1, E(3)]
      else
        res=fill(E(1),r)
        res[end]=E(o)
        return res
      end
    elseif [p, q, r, o]==[3, 3, 3, 4] return [E(4), 1, -(E(4))]
    end
  elseif p==2q return [1, E(3)]
  end
  error("wrong arguments")
end)

#Group( ( 1, 2,44)( 4, 8,49)( 5,48,24)( 6, 9,26)( 7,25,47)(10,11,50)(14,17,38)
#(15,36,43)(16,37,32)(18,19,39)(21,22,52)(41,42,54), ( 2,12,37)( 3,10,15)
#( 5,26,40)( 6,24,27)( 7,42,13)( 8,33,19)(14,21,46)(16,53,44)(18,51,49)
#(22,38,31)(28,54,47)(30,36,50)); # a section G of the Diagram automorphisms

chevieset(:timp, :ClassInfo, function (p, q, r, phi)
  if [p, q, r]==[3, 3, 3]
    if phi==perm"(1,4,2)"
      Dict{Symbol, Any}(:classes => [3, 9, 9, 3, 18, 9, 3], 
       :classtext=>[Int[],[3],[1],[2,1],[3,1],[3, 2, 1], [2, 3, 1, 2, 3, 1]], 
       :classparams=>[Int[],[3],[1],[2,1],[3,1],[3,2,1],[2, 3, 1, 2, 3, 1]], 
       :classnames => ["Id", "3", "1", "21", "31", "321", "231231"])
    elseif phi==perm"(1,2,4)"
      Dict{Symbol, Any}(:classes => [3, 9, 9, 3, 18, 9, 3],
       :classtext=>[Int[],[3],[1],[1,2],[3,1],[3,1, 2], [1, 3, 2, 1, 3, 2]], 
       :classparams=>[Int[],[3],[1],[1,2],[3,1],[3,1,2], [1, 3, 2, 1, 3, 2]], 
       :classnames => ["Id", "3", "1", "12", "31", "312", "132132"])
    elseif order(phi)==4
      Dict{Symbol, Any}(:classes => [9, 9, 9, 9, 9, 9],
       :classtext=>[Int[], [1], [1, 2], [1, 2, 3], [1, 2, 1], [1, 2, 1, 3]],
       :classparams=>[Int[], [1], [1, 2], [1, 2, 3], [1, 2, 1], [1, 2, 1, 3]],
       :classnames =>["Id", "1", "12", "123", "121", "1213"])
    elseif order(phi)==2 CHEVIE[:imp][:ClassInfo](3, 3, 3)
    else error("phi==", phi, " not implemented")
    end
  elseif [p, q, r]==[4, 2, 2]
    Dict{Symbol, Any}(:classtext=>[Int[],[1],[1,2,3,1,2,3],[1,2,3,1,2,3,1,2,3]],
                      :classes => [4, 4, 4, 4], 
                      :classnames => ["Id", "1", "cc", "z"])
  else
      ChevieErr("ClassInfo not implemented")
      false
  end
end)

chevieset(:timp, :NrConjugacyClasses, function (p, q, r, phi)
  length(chevieget(:timp, :ClassInfo)(p, q, r, phi)[:classtext]) end)

chevieset(:timp, :CharInfo, function (p, q, r, phi)
  if [p, q, r]==[3, 3, 3]
    if order(phi)==3
     res=Dict{Symbol,Any}(:charparams=>[[Int[],Int[],[3]],[Int[],Int[],[1,1,1]],
       [Int[],Int[],[2,1]],[Int[],[1,1],[1]],[Int[],[2],[1]],[Int[],[1],[1,1]],
       [Int[],[1],[2]]],:extRefl => [1, 5, 6, 2])
    elseif order(phi)==4
     res=Dict{Symbol,Any}(:charparams=>[[Int[],Int[],[3]],[Int[],Int[],[1,1,1]],
       [Int[],[1,1],[1]],[Int[],[2],[1]],[Int[],[1],[1,1]],[Int[],[1],[2]]],
                          :extRefl=>[1,4,5,2])
    else
        error("phi==", phi, " not implemented")
    end
  elseif [p, q, r]==[4, 2, 2]
   res=Dict{Symbol,Any}(:charparams=>[[Int[],Int[],[2],Int[]],
     [Int[],Int[],Int[],[1,1]],[Int[],[1],[1],Int[]],[Int[],Int[],[1],[1]]],
                        :extRefl=>[1,4,2])
  else
    ChevieErr("CharInfo not implemented")
    return false
  end
  res[:charnames]=string_partition_tuple.(res[:charparams])
  res
end)

chevieset(:timp, :CharTable, function (p, q, r, phi)
  if [p, q, r]==[3, 3, 3]
    if phi==perm"(1,4,2)"
      res=Dict(:size=>54,:order=>54,:centralizers => [18, 6, 6, 18, 3, 6, 18],
        :identifier => "3'G(3,3,3)", :name => "3'G(3,3,3)",
        :classes => [3, 9, 9, 3, 18, 9, 3], :irreducibles =>
        [1 1 1 1 1 1 1;
         1 -1 -1 1 1 -1 1;
         2 0 0 2 -1 0 2;
         -root(-3) -1 -E(3,2) 2*E(3)+E(3,2) 0 -E(3) -E(3)-2*E(3,2);
         -root(-3) 1 E(3,2) 2*E(3)+E(3,2) 0 E(3) -E(3)-2*E(3,2);
         -2*E(3)-E(3,2) -E(3,2) -1 root(-3) 0 -E(3) E(3)+2*E(3,2);
         -2*E(3)-E(3,2) E(3,2) 1 root(-3) 0 E(3) E(3)+2*E(3,2)])
      res[:irreducibles][4:5,:].*=E(3)
    elseif phi==perm"(1,2,4)"
      res=Dict(:size=>54,:order=>54,:centralizers=>[18, 6, 6, 18, 3, 6, 18],
        :identifier=>"3G(3,3,3)",:name=>"3G(3,3,3)",:classes=>[3,9,9,3,18,9,3],
        :irreducibles=>
        [1  1  1  1  1  1  1;
         1  -1  -1  1  1  -1  1;
         2  0  0  2  -1  0  2;
         root(-3) -1 -E(3) (-3-root(-3))//2 0 -E(3,2) (3-root(-3))//2;
         root(-3) 1 E(3) (-3-root(-3))//2 0 E(3,2) (3-root(-3))//2;
         (3+root(-3))//2 -E(3) -1 -root(-3) 0 -E(3,2) (-3+root(-3))//2;
         (3+root(-3))//2 E(3) 1 -root(-3) 0 E(3,2) (-3+root(-3))//2])
      res[:irreducibles][4:5,:].*=E(3,2)
    elseif order(phi)==4
      res=Dict(:size => 54, :order => 54, :centralizers => [6, 6, 6, 6, 6, 6],
        :identifier=>"4G(3,3,3)",:name=>"4G(3,3,3)",:classes=>[9,9,9,9,9,9],
        :irreducibles =>
        [1 1 1 1 1 1;
         1 -1 1 -1 -1 1;
         1 E(3) E(3,2) 1 E(3,2) E(3);
         1 -E(3) E(3,2) -1 -E(3,2) E(3);
         1 E(3,2) E(3) 1 E(3) E(3,2);
         1 -E(3,2) E(3) -1 -E(3) E(3,2)])
    else
      error("phi==", phi, " not implemented")
    end
  elseif [p, q, r]==[4, 2, 2]
    res=Dict{Symbol,Any}(:size=>16,:order=>16,:centralizers=>[4,4,4,4],
      :classes=>[4,4,4,4],:identifier=>"3G(4,2,2)",:name=>"3G(4,2,2)",
      :irreducibles => 
      [1 1 1 1;
       1 -1 1 -1;
       -1 E(4) 1 -E(4);
       -1 -E(4) 1 E(4)])
  else
    ChevieErr("CharTable not implemented")
    return false
  end
  res[:text] = "origin: Dixon's Algorithm"
  res
end)

chevieset(:timp, :UnipotentCharacters, function (p, q, r, phi)
  if [p, q, r]==[3, 3, 3]
    if phi==perm"(1,4,2)"
      Dict{Symbol, Any}(:harishChandra =>
       [Dict(:relativeType=>TypeIrred(;series=:ST,indices=[1,2],p=3,q=1,rank=2),
             :levi=>Int[],:eigenvalue=>1,:parameterExponents=>[[2,0,1],1],
             :cuspidalName=>"",:charNumbers=>[7,3,5,2,4,9,1,6,8])],
    :almostHarishChandra=>[Dict(:relativeType=>TypeIrred(;orbit=[
       TypeIrred(;series=:ST,indices=1:3,p=3,q=3,rank=3)],twist=perm"(1,2,4)"),
      :levi=>[],:eigenvalue=>1,:cuspidalName=>"",:charNumbers => 1:7), 
    mkcuspidal("G_{3,3,3}",8,E(3)),
    mkcuspidal("G_{3,3,3}",9,E(3,2))],
    :families => [Family("C1", [1]), Family("C1", [2]), Family("C1", [3]),
      Family(conj(CHEVIE[:families][:X](3)), [7, 5, 8]),
      Family(conj(CHEVIE[:families][:X](3)), [6, 4, 9])],
    :a => [0, 9, 3, 4, 1, 4, 1, 1, 4],
    :A => [0, 9, 6, 8, 5, 8, 5, 5, 8])
    elseif phi==perm"(1,2,4)"
      Dict{Symbol, Any}(:harishChandra => [
      Dict(:relativeType=>TypeIrred(;series=:ST,indices=[1,2],p=3,q=1,rank=2),
      :levi=>Int[],:eigenvalue=>1,:parameterExponents=>[[2,1,0],1],
      :cuspidalName=>"",:charNumbers=>[7, 5, 3, 9, 4, 2, 1, 8, 6])], 
    :almostHarishChandra=>[Dict{Symbol, Any}(:relativeType=>TypeIrred(;orbit=[
       TypeIrred(;series=:ST,indices=1:3,p=3,q=3,rank=3)],twist=perm"(1,2,4)"),
       :levi=>Int[],:eigenvalue=>1,:cuspidalName=>"",:charNumbers=>1:7), 
    mkcuspidal("G_{3,3,3}",8,E(3)),
    mkcuspidal("G_{3,3,3}",9,E(3,2))],
     :families => [Family("C1", [1]), Family("C1", [2]), Family("C1", [3]),
                   Family(CHEVIE[:families][:X](3), [7, 5, 8]),
                   Family(CHEVIE[:families][:X](3), [6, 4, 9])],
     :a => [0, 9, 3, 4, 1, 4, 1, 1, 4],
     :A => [0, 9, 6, 8, 5, 8, 5, 5, 8])
    elseif order(phi)==4
      res=Dict{Symbol,Any}(:harishChandra=>[Dict(:relativeType=>
        TypeIrred(;series=:ST,indices=[1],p=6,q=1,rank=1),:levi=>Int[],
        :eigenvalue=>1,:parameterExponents=>[[3, 1, 2, 0, 2, 1]],
        :cuspidalName=>"",:charNumbers=>[1, 5, 4, 2, 6, 3]), 
     mkcuspidal("{}^4G_{3,3,3}",7,E(3,2)),
     mkcuspidal("{}^4G_{3,3,3}",8,E(3))],
     :almostHarishChandra=>[Dict{Symbol, Any}(:relativeType=>TypeIrred(;orbit=[
      TypeIrred(;series=:ST,indices=1:3,p=3,q=3,rank=3)],twist=perm"(1,2,3,4)"),
      :levi=>Int[],:eigenvalue=>1,:cuspidalName => "", :charNumbers => 1:6), 
     mkcuspidal("G_{3,3,3}",7,E(3)),
     mkcuspidal("G_{3,3,3}",8,E(3,2))], 
     :families => [Family("C1", [1]), Family("C1", [2]),
        Family(conj(CHEVIE[:families][:X](3)),[3,5,7],Dict{Symbol,Any}(:signs=>[1,1,-1])),
        Family(CHEVIE[:families][:X](3),[4,6,8],Dict{Symbol,Any}(:signs=>[1,1,-1]))],
     :a => [0, 9, 4, 1, 4, 1, 4, 1],
     :A => [0, 9, 8, 5, 8, 5, 8, 5])
        res[:families][3][:eigenvalues][3]=E(3, 2)
        res[:families][4][:eigenvalues][3]=E(3)
        a=1
        res[:families][3][:fourierMat][3]*=a
        res[:families][4][:fourierMat][3]*=galois(a,-1)
      res
    else error("phi==", phi, " not implemented")
    end
  else
    if q==1 || q==p ChevieErr("UnipotentCharacters not implemented")
    end
    return false
  end
end)

# here W should have type[1] be G333 with indices [1,2,3]
# the representatives chosen are for the 24 elements of G: 15 in G and 9 in -G
chevieset(:timp, :ReducedInRightCoset, function (W, phi)
  m1perm=perm"(1,4)(2,8)(3,13)(5,22)(6,14)(7,10)(9,17)(11,25)(12,33)(15,42)
(16,18)(19,37)(20,35)(21,24)(23,29)(26,38)(27,46)(28,30)(31,40)(32,39)(34,45)
(36,54)(41,43)(44,49)(47,50)(48,52)(51,53)" # effect of -1 on roots
  m1perm=mappingPerm(inclusion(W),invpermute(inclusion(W),m1perm))
  g=PermX(W,reflrep(W,phi)*reflrep(W,m1perm))
  if g in W return (gen=inclusion(W,1:3),phi=inv(g)*phi) end
  # quads of roots which have the same cartan mat and are representatives of 
  # <-1,W>-orbits of quads of reflections satisfying Corran-Picantin relations
  sets=[[1,2,3,44],[2,12,11,37],[3,11,2,36],[1,12,10,16]] # 3G333
  for g in [Perm(),m1perm], i in inclusion.(Ref(W),sets), 
     o in [[4,1,3,2],[2,4,3,1]]
    e=transporting_elt(W,i[o],i.^(g*phi),(s,g)->s.^g)
    if !isnothing(e)return (gen=i[1:3],phi=phi/e) end
  end
  # sextuples of roots which have the same cartan mat and are representatives of 
  # W-orbits of each element of order 4 of G
  sets2=[[1,2,32,16,3,36,30,10],[3,10,30,36,17,21,53,38],
        [2,12,16,53,11,10,43,36],[2,44,16,37,3,43,30,11],
        [1,37,32,44,15,30,50,3],[1,12,32,53,10,50,36,15]] # 4G333
  for i in inclusion.(Ref(W),sets2)
    e=transporting_elt(W,i[[2,3,4,1]],i[1:4].^phi,(s,g)->s.^g)
    if !isnothing(e)return (gen=i[[1,2,5]],phi=phi/e) end
  end
end)
