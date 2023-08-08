function famf20()
  g4=Perm(2,4,5,3);g5=Perm(1,2,3,4,5)
  f20=Group(g5,g4)
  f20.classreps=[Perm(),g4^3,g4,g4^2,g5]
  f20.chartable=CharTable(
    [1 1 1 1 1;1 -1 -1 1 1;1 -E(4) E(4) -1 1;1 E(4) -E(4) -1 1;4 0 0 0 -1],
    ["1","-1","i","-i","\\rho"],["1","g_4^3","g_4","g_2","g_5"],
    [20,4,4,4,5],20,Dict{Symbol,Any}(:name=>"F20"))
  f20.name="F20"
  drinfeld_double(f20)
end
