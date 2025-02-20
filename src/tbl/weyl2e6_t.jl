#  tbl/weyl2e6.jl    CHEVIE library           Frank Luebeck and Jean Michel
#  Copyright (C) 1994 - 2001  The CHEVIE Team

chevieset(Symbol("2E6"), :NrConjugacyClasses, 25)

## these minimal length representatives were found by
#  W=coxgroup(:E,6)
#  short=map(conjugacy_classes(W))do c
#    word(W,argmin(x->length(W,x),elements(c).*longest(W)))
#  end
chevieset(Symbol("2E6"), :ClassInfo, function ()
  res=Dict{Symbol,Any}(:classtext=>[
    [1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,3,5,4,2,6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1],
    Int[],[3,4,3,5,4,3],[1,2,4,3,1,5,4,3,6,5,4,3],
    [1,2,3,1,4,3,1,5,4,3,1,6,5,4,3,1],[2,3,4,2,3,4,6,5,4,2,3,4,5,6],
    [1,4,2,3,1,4,3,5,4,2,3,1,4,6,5,4,3,1],[1,2],[4,5,4,2,3,1,4,5],
    [4,2,5,4,2,3,4,5,6,5,4,2,3,4,5,6],[2,4],[1,5],[5,4],[1,2,5,4],
    [1,2,3,1,4,3],[1,3,1,4,3,1,5,4,3,1,6,5,4,3,1],[2],[1],[2,3,4,3,5,4,3],
    [1,3,4,3,5,4,3],[1,3,1,4,3],[1,2,5],[2,5,4],[1,5,4],[1,2,4]],
  :classnames=>["A_0", "4A_1", "2A_1", "3A_2", "A_2", "2A_2", "D_4(a_1)",
    "A_3+A_1", "A_4", "E_6(a_2)", "D_4", "A_5+A_1", "A_2+2A_1", "E_6(a_1)",
    "E_6", "A_1", "3A_1", "A_3+2A_1", "A_3", "A_2+A_1", "2A_2+A_1", "A_5",
    "D_5", "A_4+A_1", "D_5(a_1)"],
  :classes=>[1, 45, 270, 80, 240, 480, 540, 3240, 5184, 720, 1440, 1440,
    2160, 5760, 4320, 36, 540, 540, 1620, 1440, 1440, 4320, 6480, 5184, 4320])
  res[:classparams] = res[:classnames]
  res
end)

chevieset(Symbol("2E6"), :CharInfo, function ()
  res=copy(chevieget(:E6, :CharInfo)())
  res[:a]=[0,36,7,1,25,7,3,15,3,15,2,20,6,12,3,15,7,7,7,5,11,4,13,6,10]
  res[:A]=[0,36,29,11,35,29,21,33,21,33,16,34,24,30,21,33,29,29,29,25,31,23,
           32,26,30]
  res[:B]=[0,36,27,11,35,26,19,31,20,32,16,34,24,30,21,33,28,29,28,25,31,23,
           32,26,30]
  res
end)

chevieset(Symbol("2E6"), :cyclestructure, [Pair{Int64, Int64}[], [2 => 36],
  [2 => 30], [3 => 24], [3 => 20], [3 => 22], [4 => 18], [2 => 5, 4 => 15],
  [5 => 14], [6 => 12], [2 => 6, 6 => 10], [2 => 3, 6 => 11],
  [2 => 6, 3 => 4, 6 => 8], [9 => 8], [12 => 6], [2 => 21], [2 => 35],
  [2 => 6, 4 => 15], [2 => 4, 4 => 15], [2 => 3, 3 => 8, 6 => 6],
  [2 => 3, 3 => 10, 6 => 6], [2 => 2, 6 => 11], [8 => 9],
  [2 => 1, 5 => 6, 10 => 4], [4 => 3, 6 => 2, 12 => 4]])

chevieset(Symbol("2E6"), :generators, [
  perm"(1,37)(3,7)(9,12)(13,17)(15,18)(19,22)(21,23)(24,26)(25,27)(28,30)(31,33)(39,43)(45,48)(49,53)(51,54)(55,58)(57,59)(60,62)(61,63)(64,66)(67,69)", 
  perm"( 2,38)( 4, 8)( 9,13)(10,14)(12,17)(15,19)(16,20)(18,22)(21,25)(23,27)(35,36)(40,44)(45,49)(46,50)(48,53)(51,55)(52,56)(54,58)(57,61)(59,63)(71,72)", 
  perm"( 1, 7)( 3,39)( 4, 9)( 8,13)(10,15)(14,19)(16,21)(20,25)(26,29)(30,32)(33,34)(37,43)(40,45)(44,49)(46,51)(50,55)(52,57)(56,61)(62,65)(66,68)(69,70)", 
  perm"( 2, 8)( 3, 9)( 4,40)( 5,10)( 7,12)(11,16)(19,24)(22,26)(25,28)(27,30)(34,35)(38,44)(39,45)(41,46)(43,48)(47,52)(55,60)(58,62)(61,64)(63,66)(70,71)", 
  perm"( 4,10)( 5,41)( 6,11)( 8,14)( 9,15)(12,18)(13,19)(17,22)(28,31)(30,33)(32,34)(40,46)(42,47)(44,50)(45,51)(48,54)(49,55)(53,58)(64,67)(66,69)(68,70)", 
  perm"( 5,11)( 6,42)(10,16)(14,20)(15,21)(18,23)(19,25)(22,27)(24,28)(26,30)(29,32)(41,47)(46,52)(50,56)(51,57)(54,59)(55,61)(58,63)(60,64)(62,66)(65,68)"])

chevieset(Symbol("2E6"), :phi, perm"(1,42)(2,38)(3,41)(4,40)(5,39)(6,37)(7,47)(8,44)(9,46)(10,45)(11,43)(12,52)(13,50)(14,49)(15,51)(16,48)(17,56)(18,57)(19,55)(20,53)(21,54)(22,61)(23,59)(24,60)(25,58)(26,64)(27,63)(28,62)(29,67)(30,66)(31,65)(32,69)(33,68)(34,70)(35,71)(36,72)")

chevieset(Symbol("2E6"), :CartanMat, 
  [2 0 -1 0 0 0;0 2 0 -1 0 0;-1 0 2 -1 0 0;0 -1 -1 2 -1 0;0 0 0 -1 2 -1;
   0 0 0 0 -1 2])

chevieset(Symbol("2E6"), :vpolheckeirreducibles, 
Vector{Pol{Int64}}[[Pol([1],72), Pol([1]), Pol([1],12), Pol([1],24), Pol([1],
32), Pol([1],28), Pol([1],36), Pol([1],4), Pol([1],16), Pol([1],32), Pol([1],
4), Pol([1],4), Pol([1],4), Pol([1],8), Pol([1],12), Pol([1],30), Pol([1],2),
 Pol([1],2), Pol([1],14), Pol([1],14), Pol([1],10), Pol([1],6), Pol([1],6),
 Pol([1],6), Pol([1],6)], [Pol([1]), Pol([1]), Pol([1]), Pol([1]), Pol([1]),
 Pol([1]), Pol([1]), Pol([1]), Pol([1]), Pol([1]), Pol([1]), Pol([1]),
 Pol([1]), Pol([1]), Pol([1]), Pol([-1]), Pol([-1]), Pol([-1]), Pol([-1]),
 Pol([-1]), Pol([-1]), Pol([-1]), Pol([-1]), Pol([-1]), Pol([-1])],
 [Pol([-10],36), Pol([6]), Pol([-3, 0, 4, 0, -3],4), Pol([-1],12), Pol([-1, 0,
 0, 0, 4, 0, 0, 0, -1],12), Pol([-1, 0, -2, 0, -1],12), Pol([-2],18), Pol([2,
 0, -2, 0, 2]), Pol([0]), Pol([3],16), Pol([1, 0, -2, 0, 1]), Pol([1, 0, -2,
 0, 1]), Pol([-2],2), Pol([-1],4), Pol([1],6), Pol([5, 0, 0, 0, 0, 0, -5],12),
 Pol([-3, 0, 3]), Pol([-3, 0, 3]), Pol([-1, 0, 1],6), Pol([1, 0, -2, 0, 2, 0,
 -1],4), Pol([2, 0, -2],4), Pol([-1, 0, 1, 0, -1, 0, 1]), Pol([0]), Pol([0]),
 Pol([-1, 0, 1, 0, -1, 0, 1])], [Pol([-6],60), Pol([2]), Pol([-3, 0, 0, 0, 1],
8), Pol([3],20), Pol([-4, 0, 1],26), Pol([-2, 0, 2],22), Pol([-2],30),
 Pol([-1, 0, 1],2), Pol([-1],12), Pol([-2, 0, 1],26), Pol([-1],2), Pol([2],4),
 Pol([1],4), Pol([0]), Pol([1],10), Pol([-5, 0, 0, 0, 0, 0, 1],24), Pol([-1,
 0, 1]), Pol([2],2), Pol([-2],10), Pol([-2, 0, 0, 0, 1],10), Pol([1, 0, 1],8),
 Pol([-1, 0, 1],4), Pol([0]), Pol([1],6), Pol([-1],4)], [Pol([-6],12),
 Pol([2]), Pol([1, 0, 0, 0, -3]), Pol([3],4), Pol([1, 0, -4],4), Pol([2, 0,
 -2],4), Pol([-2],6), Pol([1, 0, -1]), Pol([-1],4), Pol([1, 0, -2],4),
 Pol([-1],2), Pol([2]), Pol([1]), Pol([0]), Pol([1],2), Pol([-1, 0, 0, 0, 0,
 0, 5]), Pol([-1, 0, 1]), Pol([-2]), Pol([2],4), Pol([-1, 0, 0, 0, 2]),
 Pol([-1, 0, -1]), Pol([-1, 0, 1]), Pol([0]), Pol([-1]), Pol([1],2)],
 [Pol([-20],36), Pol([-4]), Pol([-1, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, -1]),
 Pol([7],12), Pol([-4, 0, 6, 0, -4],14), Pol([2, 0, -6, 0, 2],12), Pol([-4],
18), Pol([-1, 0, 2, 0, -1]), Pol([0]), Pol([-2, 0, 3, 0, -2],14), Pol([2],2),
 Pol([2],2), Pol([-1, 0, 0, 0, -1]), Pol([1],4), Pol([-1],6), Pol([-10, 0, 0,
 0, 0, 0, 10],12), Pol([2, 0, -2]), Pol([2, 0, -2]), Pol([-2, 0, 0, 0, 0, 0,
 2],4), Pol([1, 0, -1],6), Pol([-1, 0, 1],4), Pol([-1, 0, 1],2), Pol([0]),
 Pol([0]), Pol([-1, 0, 1],2)], [Pol([-15],48), Pol([1]), Pol([-3, 0, 0, 0, 3,
 0, 0, 0, 1],4), Pol([-6],16), Pol([-6, 0, 4, 0, -1],20), Pol([-1, 0, 4, 0,
 -3],16), Pol([-3],24), Pol([1],4), Pol([0]), Pol([-3, 0, 2, 0, -1],20),
 Pol([-1, 0, 0, 0, 2]), Pol([-1, 0, -1],2), Pol([1],4), Pol([0]), Pol([0]),
 Pol([-10, 0, 0, 0, 0, 0, 5],18), Pol([1, 0, 2]), Pol([-1]), Pol([-1, 0, -1,
 0, 0, 0, 0, 0, 1],6), Pol([-1, 0, 0, 0, 2],6), Pol([-1, 0, -1],6), Pol([0]),
 Pol([1],6), Pol([0]), Pol([-1, 0, -1, 0, 1],2)], [Pol([-15],24), Pol([1]),
 Pol([1, 0, 0, 0, 3, 0, 0, 0, -3]), Pol([-6],8), Pol([-1, 0, 4, 0, -6],8),
 Pol([-3, 0, 4, 0, -1],8), Pol([-3],12), Pol([1]), Pol([0]), Pol([-1, 0, 2, 0,
 -3],8), Pol([2, 0, 0, 0, -1]), Pol([-1, 0, -1]), Pol([1]), Pol([0]),
 Pol([0]), Pol([-5, 0, 0, 0, 0, 0, 10],6), Pol([-2, 0, -1]), Pol([1],2),
 Pol([-1, 0, 0, 0, 0, 0, 1, 0, 1]), Pol([-2, 0, 0, 0, 1],4), Pol([1, 0, 1],2),
 Pol([0]), Pol([-1]), Pol([0]), Pol([-1, 0, 1, 0, 1])], [Pol([-15],48),
 Pol([-7]), Pol([-4, 0, 3, 0, 0, 0, -2],6), Pol([3],16), Pol([4, 0, -4],20),
 Pol([-4, 0, 1],18), Pol([1],24), Pol([-1, 0, 3, 0, -3]), Pol([0]), Pol([1, 0,
 -2],20), Pol([3, 0, -1],2), Pol([2, 0, -3],2), Pol([2, 0, -2],2), Pol([0]),
 Pol([1],8), Pol([5, 0, -9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],18), Pol([3, 0,
 -4]), Pol([2, 0, -5]), Pol([1, 0, -2, 0, 2],6), Pol([-2, 0, 1, 0, 0, 0, -1],
8), Pol([2, 0, 0, 0, -1],6), Pol([-1, 0, 2, 0, -2],2), Pol([1],4), Pol([1, 0,
 -1],4), Pol([-1, 0, 2, 0, -1],2)], [Pol([-15],24), Pol([-7]), Pol([-2, 0, 0,
 0, 3, 0, -4]), Pol([3],8), Pol([-4, 0, 4],10), Pol([1, 0, -4],8), Pol([1],
12), Pol([-3, 0, 3, 0, -1]), Pol([0]), Pol([-2, 0, 1],10), Pol([-1, 0, 3]),
 Pol([-3, 0, 2]), Pol([-2, 0, 2]), Pol([0]), Pol([1],4), Pol([1, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 9, 0, -5]), Pol([4, 0, -3]), Pol([5, 0, -2]), Pol([-2, 0, 2,
 0, -1],4), Pol([1, 0, 0, 0, -1, 0, 2]), Pol([1, 0, 0, 0, -2]), Pol([2, 0, -2,
 0, 1]), Pol([-1],2), Pol([1, 0, -1]), Pol([1, 0, -2, 0, 1])], [Pol([20],54),
 Pol([4]), Pol([2, 0, 0, 0, 0, 0, 2],6), Pol([2],18), Pol([4, 0, 0, 0, 1],22),
 Pol([-2, 0, 1],20), Pol([0]), Pol([-2, 0, 2],2), Pol([0]), Pol([-2],24),
 Pol([-1, 0, 2],2), Pol([-1, 0, 2],2), Pol([-1, 0, 2],2), Pol([-1],6),
 Pol([0]), Pol([9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],20), Pol([-1, 0, 3]),
 Pol([-1, 0, 3]), Pol([1, 0, 0, 0, 0, 0, 1],8), Pol([1],14), Pol([1],10),
 Pol([-2, 0, 1],4), Pol([-1, 0, 1],4), Pol([-1, 0, 1],4), Pol([-2, 0, 1],4)],
 [Pol([20],18), Pol([4]), Pol([2, 0, 0, 0, 0, 0, 2]), Pol([2],6), Pol([1, 0,
 0, 0, 4],6), Pol([1, 0, -2],6), Pol([0]), Pol([2, 0, -2]), Pol([0]),
 Pol([-2],8), Pol([2, 0, -1]), Pol([2, 0, -1]), Pol([2, 0, -1]), Pol([-1],2),
 Pol([0]), Pol([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -9]), Pol([-3, 0, 1]), Pol([-3,
 0, 1]), Pol([-1, 0, 0, 0, 0, 0, -1]), Pol([-1]), Pol([-1]), Pol([-1, 0, 2]),
 Pol([-1, 0, 1]), Pol([-1, 0, 1]), Pol([-1, 0, 2])], [Pol([24],42), Pol([8]),
 Pol([-6, 0, 4, 0, 0, 0, 0, 0, 2],4), Pol([6],14), Pol([1, 0, 0, 0, 2, 0, -4,
 0, 1],14), Pol([2, 0, -2, 0, 3],14), Pol([0]), Pol([1, 0, -4, 0, 3]),
 Pol([-1],10), Pol([1, 0, 0, 0, 1],18), Pol([-2, 0, 4],2), Pol([1, 0, -3, 0,
 1]), Pol([-2, 0, 2],2), Pol([0]), Pol([0]), Pol([5, 0, 0, 0, 0, 0, -10, 0,
 9],12), Pol([-2, 0, 6]), Pol([-4, 0, 4]), Pol([-1, 0, 0, 0, 0, 0, 0, 0, 1],
6), Pol([1, 0, -3],4), Pol([2, 0, -1],4), Pol([2, 0, -2, 0, 1],2), Pol([-1, 0,
 1],4), Pol([-1],4), Pol([1, 0, -3, 0, 2],2)], [Pol([24],30), Pol([8]),
 Pol([2, 0, 0, 0, 0, 0, 4, 0, -6]), Pol([6],10), Pol([1, 0, -4, 0, 2, 0, 0, 0,
 1],10), Pol([3, 0, -2, 0, 2],10), Pol([0]), Pol([3, 0, -4, 0, 1]), Pol([-1],
6), Pol([1, 0, 0, 0, 1],10), Pol([4, 0, -2]), Pol([1, 0, -3, 0, 1]), Pol([2,
 0, -2]), Pol([0]), Pol([0]), Pol([-9, 0, 10, 0, 0, 0, 0, 0, -5],10), Pol([-6,
 0, 2]), Pol([-4, 0, 4]), Pol([-1, 0, 0, 0, 0, 0, 0, 0, 1]), Pol([3, 0, -1],
8), Pol([1, 0, -2],4), Pol([-1, 0, 2, 0, -2]), Pol([-1, 0, 1]), Pol([1],2),
 Pol([-2, 0, 3, 0, -1])], [Pol([-30],48), Pol([10]), Pol([-3, 0, 4, 0, -6, 0,
 0, 0, 3],4), Pol([-3],16), Pol([-2, 0, 0, 0, -1],20), Pol([-1, 0, 0, 0, -2],
16), Pol([2],24), Pol([1, 0, -5, 0, 4]), Pol([0]), Pol([2, 0, 0, 0, -1],20),
 Pol([1, 0, -3, 0, 3]), Pol([-3, 0, 4],2), Pol([-2, 0, 3],2), Pol([0]),
 Pol([-1],8), Pol([-15, 0, 9, 0, 0, 0, -5, 0, 0, 0, 0, 0, 1],18), Pol([-4, 0,
 6]), Pol([-3, 0, 7]), Pol([1, 0, -2, 0, 0, 0, 1],8), Pol([-1, 0, 2, 0, -3, 0,
 0, 0, 1],6), Pol([-3, 0, 1, 0, 1],6), Pol([1, 0, -4, 0, 2],2), Pol([-1, 0,
 1],4), Pol([-1, 0, 1],4), Pol([2, 0, -3, 0, 2],2)], [Pol([-30],24),
 Pol([10]), Pol([3, 0, 0, 0, -6, 0, 4, 0, -3]), Pol([-3],8), Pol([-1, 0, 0, 0,
 -2],8), Pol([-2, 0, 0, 0, -1],8), Pol([2],12), Pol([4, 0, -5, 0, 1]),
 Pol([0]), Pol([-1, 0, 0, 0, 2],8), Pol([3, 0, -3, 0, 1]), Pol([4, 0, -3]),
 Pol([3, 0, -2]), Pol([0]), Pol([-1],4), Pol([-1, 0, 0, 0, 0, 0, 5, 0, 0, 0,
 -9, 0, 15]), Pol([-6, 0, 4]), Pol([-7, 0, 3]), Pol([-1, 0, 0, 0, 2, 0, -1]),
 Pol([-1, 0, 0, 0, 3, 0, -2, 0, 1]), Pol([-1, 0, -1, 0, 3]), Pol([-2, 0, 4, 0,
 -1]), Pol([-1, 0, 1]), Pol([-1, 0, 1]), Pol([-2, 0, 3, 0, -2])], [Pol([-60],
36), Pol([-12]), Pol([-1, 0, 0, 0, 3, 0, -8, 0, 3, 0, 0, 0, -1]), Pol([3],12),
 Pol([1, 0, 0, 0, 4, 0, 0, 0, 1],12), Pol([1, 0, -2, 0, 1],12), Pol([-4],18),
 Pol([-3, 0, 6, 0, -3]), Pol([0]), Pol([3],16), Pol([-2, 0, 4, 0, -2]),
 Pol([-2, 0, 4, 0, -2]), Pol([-1, 0, 4, 0, -1]), Pol([0]), Pol([-1],6),
 Pol([9, 0, -5, 0, 0, 0, 0, 0, 5, 0, -9],10), Pol([6, 0, -6]), Pol([6, 0,
 -6]), Pol([2, 0, -2],6), Pol([-1, 0, 2, 0, -2, 0, 1],4), Pol([-2, 0, 2],4),
 Pol([1, 0, -3, 0, 3, 0, -1]), Pol([-1, 0, 1],2), Pol([-1, 0, 1],2), Pol([1,
 0, -3, 0, 3, 0, -1])], [Pol([-80],36), Pol([16]), Pol([2, 0, 0, 0, -6, 0, 8,
 0, -6, 0, 0, 0, 2]), Pol([10],12), Pol([1, 0, -4, 0, 10, 0, -4, 0, 1],12),
 Pol([3, 0, -8, 0, 3],12), Pol([0]), Pol([4, 0, -8, 0, 4]), Pol([0]), Pol([-2,
 0, 2, 0, -2],14), Pol([2, 0, -6, 0, 2]), Pol([2, 0, -6, 0, 2]), Pol([2, 0,
 -4, 0, 2]), Pol([1],4), Pol([0]), Pol([-9, 0, 15, 0, 0, 0, 0, 0, -15, 0, 9],
10), Pol([-8, 0, 8]), Pol([-8, 0, 8]), Pol([2, 0, -2, 0, 2, 0, -2],4), Pol([1,
 0, -3, 0, 3, 0, -1],4), Pol([3, 0, -3],4), Pol([-1, 0, 4, 0, -4, 0, 1]),
 Pol([1, 0, -1],2), Pol([1, 0, -1],2), Pol([-1, 0, 4, 0, -4, 0, 1])],
 [Pol([-90],36), Pol([6]), Pol([1, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 1]),
 Pol([-9],12), Pol([-1, 0, 4, 0, -6, 0, 4, 0, -1],12), Pol([-3, 0, 6, 0, -3],
12), Pol([-2],18), Pol([1, 0, -4, 0, 1]), Pol([0]), Pol([2, 0, -1, 0, 2],14),
 Pol([1, 0, -2, 0, 1]), Pol([1, 0, -2, 0, 1]), Pol([1, 0, -2, 0, 1]),
 Pol([0]), Pol([1],6), Pol([-9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9],10), Pol([-3, 0,
 3]), Pol([-3, 0, 3]), Pol([-1, 0, 1],6), Pol([0]), Pol([0]), Pol([2, 0, -2],
2), Pol([1, 0, -1],2), Pol([1, 0, -1],2), Pol([2, 0, -2],2)], [Pol([-60],42),
 Pol([4]), Pol([2, 0, -6],6), Pol([-6],14), Pol([4, 0, -4, 0, 4, 0, -1],16),
 Pol([-1, 0, 6, 0, -2],14), Pol([0]), Pol([1, 0, -2, 0, 1]), Pol([0]),
 Pol([-4, 0, 2],18), Pol([2, 0, -1]), Pol([-1, 0, 2],2), Pol([-1],2),
 Pol([0]), Pol([0]), Pol([-5, 0, 0, 0, 0, 0, -5],18), Pol([-3, 0, 1]),
 Pol([-1, 0, 3]), Pol([1, 0, 1],6), Pol([2, 0, -3],8), Pol([-2, 0, 1],6),
 Pol([1, 0, -1, 0, 1],2), Pol([0]), Pol([0]), Pol([-1, 0, 1, 0, -1])],
 [Pol([-60],30), Pol([4]), Pol([-6, 0, 2],4), Pol([-6],10), Pol([-1, 0, 4, 0,
 -4, 0, 4],10), Pol([-2, 0, 6, 0, -1],10), Pol([0]), Pol([1, 0, -2, 0, 1]),
 Pol([0]), Pol([2, 0, -4],12), Pol([-1, 0, 2],2), Pol([2, 0, -1]), Pol([-1],
2), Pol([0]), Pol([0]), Pol([5, 0, 0, 0, 0, 0, 5],6), Pol([-1, 0, 3]),
 Pol([-3, 0, 1]), Pol([-1, 0, -1],6), Pol([3, 0, -2],4), Pol([-1, 0, 2],2),
 Pol([-1, 0, 1, 0, -1]), Pol([0]), Pol([0]), Pol([1, 0, -1, 0, 1],2)],
 [Pol([64],45), Pol([0]), Pol([0]), Pol([-8],15), Pol([4, 0, -2, 0, 4, 0, -2],
17), Pol([-2, 0, 4, 0, -4],15), Pol([0]), Pol([0]), Pol([-1],9), Pol([-1, 0,
 2, 0, -1],19), Pol([0]), Pol([0]), Pol([0]), Pol([1],5), Pol([0]), Pol([16],
15), Pol([0]), Pol([0]), Pol([0]), Pol([-1, 0, 0, 0, -1],5), Pol([-2],5),
 Pol([0]), Pol([0]), Pol([1],3), Pol([0])], [Pol([-64],27), Pol([0]),
 Pol([0]), Pol([8],9), Pol([2, 0, -4, 0, 2, 0, -4],9), Pol([4, 0, -4, 0, 2],
9), Pol([0]), Pol([0]), Pol([1],7), Pol([1, 0, -2, 0, 1],9), Pol([0]),
 Pol([0]), Pol([0]), Pol([-1],3), Pol([0]), Pol([16],15), Pol([0]), Pol([0]),
 Pol([0]), Pol([-1, 0, 0, 0, -1],5), Pol([-2],5), Pol([0]), Pol([0]), Pol([1],
3), Pol([0])], [Pol([81],40), Pol([9]), Pol([1, 0, 0, 0, 0, 0, 4, 0, -9, 0, 0,
 0, 1]), Pol([0]), Pol([0]), Pol([0]), Pol([-3],20), Pol([2, 0, -5, 0, 2]),
 Pol([1],8), Pol([0]), Pol([3, 0, -3]), Pol([-3, 0, 3],2), Pol([1, 0, -2, 0,
 1]), Pol([0]), Pol([0]), Pol([10, 0, 0, 0, 0, 0, -5, 0, 9, 0, 0, 0, -5],12),
 Pol([-6, 0, 3]), Pol([-3, 0, 6]), Pol([1, 0, -2],8), Pol([3, 0, -3],8),
 Pol([1, 0, -2, 0, 1],4), Pol([2, 0, -3, 0, 1],2), Pol([1],2), Pol([-1],4),
 Pol([-1, 0, 3, 0, -2])], [Pol([81],32), Pol([9]), Pol([1, 0, 0, 0, -9, 0, 4,
 0, 0, 0, 0, 0, 1]), Pol([0]), Pol([0]), Pol([0]), Pol([-3],16), Pol([2, 0,
 -5, 0, 2]), Pol([1],8), Pol([0]), Pol([-3, 0, 3],2), Pol([3, 0, -3]), Pol([1,
 0, -2, 0, 1]), Pol([0]), Pol([0]), Pol([5, 0, 0, 0, -9, 0, 5, 0, 0, 0, 0, 0,
 -10],6), Pol([-3, 0, 6]), Pol([-6, 0, 3]), Pol([2, 0, -1],4), Pol([3, 0, -3],
4), Pol([-1, 0, 2, 0, -1],2), Pol([-1, 0, 3, 0, -2]), Pol([-1],4), Pol([1],2),
 Pol([2, 0, -3, 0, 1],2)]])

chevieset(Symbol("2E6"), :FakeDegree, function (p, q)
  sgns=[1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,1]
  sgns[findfirst(==(p),chevieget(:E6, :CharInfo)()[:charparams])]*
    chevieget(:E6,:FakeDegree)(p,-q)
end)

chevieset(Symbol("2E6"), :ClassParameter, function (w,)
  x=prod(chevieget("2E6",:generators)[w],init=Perm())*chevieget("2E6",:phi)
  chevieget("2E6",:ClassNames)[findfirst(==(tally(classtype(x))),
                               chevieget("2E6",:cyclestructure))]
end)

if false
# The table below was obtained by the following program:
function getHeckeCharTable2E6(v)
  W=coxgroup(:E,6)
# Some characters of W(E_6)*w_0 will be changed by sign according to
# the preferred extension, see [Lusztig-book, 4.1 and 4.11] and
# [CS, 17.2(b)]
  aE=map(x->(-1)^x,chevieget(:E6,:LowestPowerGenericDegrees)())
  qE=central_monomials(hecke(W,v))
# q_E is the square root which deforms to 1 of the eigenvalue of T_{w_0}
# on E which deforms to 1; we have:
#  E~(T_w\phi)=\overline(E(T_{w^-1w_0}))q_E (trivial extension)
#  E~(T_w\phi)=a_E\overline(E(T_{w^-1w_0}))q_E (preferred extension)
# where \overline means q->q^-1
  H=hecke(W,v^-2)
  tbl=copy(CharTable(H))
  merge!(tbl,chevieget("2E6",:ClassInfo)())
  tbl.identifier="H(^2E6)"
  cl=map(x->W(x...)*longest(W),tbl.classtext)
  tbl.irreducibles=transpose(map(x->char_values(Tbasis(H)(x)),cl))
  for i in axes(tbl.irreducibles,1)
    tbl.irreducibles[i,:]=qE[i]*aE[i]*tbl.irreducibles[i,:]
  end
  tbl
end
end

# We give the values of the *preferred* extensions defined in [CS,17.2 (b)].
chevieset(Symbol("2E6"), :HeckeCharTable, function (param, rootparam)
  q=-param[1][1]//param[1][2]
  if rootparam[1]===nothing v=root(q)
  else v=rootparam[1]
  end
  tbl=Dict{Symbol, Any}(:identifier => "H(^2E6)",
    :text => "origin: Jean Michel, June 1996",
    :parameter=>map(i->[q,-1],1:6),
    :sqrtParameter=>fill(v,6),:size=>51840,
    :cartan=>chevieget(Symbol("2E6"), :CartanMat),
    :irreducibles => map(i->map(j->j(v),i),
      chevieget(Symbol("2E6"),:vpolheckeirreducibles)),
    :irredinfo=>chevieget(Symbol("2E6"), :IrredInfo))
  merge!(tbl,chevieget(Symbol("2E6"),:ClassInfo)())
  tbl[:centralizers]=div.(tbl[:size],tbl[:classes])
  chevieget(:compat,:AdjustHeckeCharTable)(tbl, param)
  tbl
end)

chevieset(Symbol("2E6"),:HeckeRepresentation,function(param,sqrtparam,i)
  W=coxgroup(:E,6)
  H=hecke(W,-param[1][1]//param[1][2])
  gens=chevieget(:E6, :HeckeRepresentation)(param, sqrtparam, i)
  (gens=gens,F=Matrix(prod(gens[[1,4,6,3,2,5]])^6)//
    root(central_monomials(H)[i])*
    (-1)^chevieget(Symbol("2E6"),:CharInfo)()[:a][i])
end)

chevieset(Symbol("2E6"), :Representation,i->
  chevieget(Symbol("2E6"),:HeckeRepresentation)(fill([1,-1],6),fill(1,6),i))

chevieset(Symbol("2E6"),:PhiFactors, [1, -1, 1, 1, -1, 1])

chevieset(Symbol("2E6"), :UnipotentCharacters,
  Dict{Symbol, Any}(:harishChandra=>[
  Dict(:relativeType=>TypeIrred(;series=:F,indices=[2,4,5,6],rank=4),
    :levi=>Int[],:eigenvalue=>1,:parameterExponents=>[1,1,2,2],
    :cuspidalName=>"",:charNumbers=>
    [1,9,10,2,4,5,15,16,17,7,24,25,8,3,19,6,11,20,21,12,26,27,13,14,28]), 
  Dict(:relativeType=>TypeIrred(;series=:A,indices=[2],rank=1),
    :levi => [1, 3, 4, 5, 6], :eigenvalue => -1, :parameterExponents => [9],
    :cuspidalName => "{}^2A_5", :charNumbers => [23, 22]), 
  Dict(:relativeType=>TypeIrred(;series=:A,indices=[],rank=0),
    :levi => 1:6, :eigenvalue => 1, :parameterExponents => [],
    :cuspidalName => "{}^2E_6[1]", :charNumbers => [18]), 
  Dict(:relativeType=>TypeIrred(;series=:A,indices=[],rank=0),
    :levi => 1:6, :eigenvalue => E(3), :parameterExponents => [],
    :cuspidalName => "{}^2E_6[\\zeta_3]", :charNumbers => [29]), 
  Dict(:relativeType=>TypeIrred(;series=:A,indices=[],rank=0),
    :levi => 1:6, :eigenvalue => E(3, 2), :parameterExponents => [],
    :cuspidalName => "{}^2E_6[\\zeta_3^2]", :charNumbers => [30])],
  :almostHarishChandra=>[Dict(:relativeType=>TypeIrred(;orbit=
     [TypeIrred(;series=:E,indices=1:6,rank=6)],twist=perm"(1,6)(3,5)"),
    :levi=>Int[],:eigenvalue=>1,:cuspidalName=>"",:charNumbers=>
    [1,2,3,15,16,6,7,8,10,9,11,12,26,27,4,5,17,18,19,21,20,22,23,25,24]), 
  Dict(:relativeType=>TypeIrred(;orbit=
     [TypeIrred(;series=:A,indices=[1,6],rank=2)],twist=perm"(1,2)"),
     :levi=>2:5,:eigenvalue=>-1,:cuspidalName=>"D_4",:charNumbers=>[14,28,13]), 
  Dict(:relativeType=>TypeIrred(;series=:A,indices=[],rank=0),:levi=>1:6,
     :eigenvalue=>E(3),:cuspidalName=>"E_6[\\zeta_3]",:charNumbers=>[29]), 
  Dict(:relativeType=>TypeIrred(;series=:A,indices=[],rank=0),:levi => 1:6,
     :eigenvalue=>E(3,2),:cuspidalName=>"E_6[\\zeta_3^2]",:charNumbers=>[30])],
  :families=>[Family("C1", [1]), Family("C1", [2]), Family("C1", [15]),
          Family("C1", [16]), Family("C1", [11]), Family("C1", [12]),
          Family("C1", [26]), Family("C1", [27]), Family("C1", [21]),
          Family("C1", [20]), Family("C'1", [22]), Family("C'1", [23]),
          Family("C1", [25]), Family("C1", [24]),
   # change last eigenvalue, and negate last line...
   Family("C2",[4,10,7,13],Dict(:name=>"Cd_2",:eigenvalues=>[1,1,1,1],
    :fourierMat=>1//2*[1 1 1 1;1 1 -1 -1;1 -1 1 -1;-1 1 1 -1],
    :sh => [1, 1, 1, 1])),
   # change last eigenvalue
   Family("C2",[5,9,8,14],Dict(:eigenvalues=>[1,1,1,1],:name=>"Cc_2",
     :sh => [1, 1, 1, 1])),
   # change 6th eigenvalue
   Family("S3", [18, 17, 3, 19, 6, 28, 29, 30], Dict(:name => "S3b",
     :eigenvalues=>[1, 1, 1, 1, 1, 1, E(3),E(3,2)],
     :sh => [1, 1, 1, 1, 1, 1, E(3, 2), E(3)]))],
  :a=>[0,36,7,3,15,7,3,15,15,3,2,20,3,15,1,25,7,7,7,11,5,4,13,10,6,6,12,7,7,7],
  :A=>[0, 36, 29, 21, 33, 29, 21, 33, 33, 21, 16, 34, 21, 33, 11, 35, 29, 29,
       29, 31, 25, 23, 32, 30, 26, 24, 30, 29, 29, 29])
)

chevieset(Symbol("2E6"), :Ennola, SPerm())

chevieset(Symbol("2E6"), :UnipotentClasses, function (p,)
  uc=deepcopy(chevieget(:E6, :UnipotentClasses)(p))
  l=[("1", perm"(1,6)(3,5)"), ("A_1", perm"(1,5)(2,4)"), 
     ("A_2", perm"(1,2)(3,4)"), ("D_4", perm"(1,2)"), ("D_5", [-1;;]), 
     ("D_5(a_1)", [-1;;]), ("A_4{+}A_1", [-1;;]), ("A_4", [1 0;0 -1]), 
     ("D_4(a_1)", [-1 0;0 -1]), ("A_3{+}A_1",[1 0;0 -1]), 
     ("A_2{+}2A_1", [1 0;0 -1]), ("A_3", Diagonal([1,1,-1])), 
     ("A_2{+}A_1", [0 1 0;1 0 0;0 0 -1]), ("3A_1", perm"(1,2)"), 
     ("2A_1", Diagonal([1,1,1,-1]))]
  for p in l
    u=uc[:classes][findfirst( x->x[:name]==p[1],uc[:classes])]
    u[:red]=spets(u[:red],p[2])
  end
  return uc
end)
