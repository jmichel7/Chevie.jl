# auto-generated tests from julia-repl docstrings
using Test, Gapjm
#include("../tools/Gap4.jl")
function mytest(f::String,a::String,b::String)
  println(f," ",a)
  omit=a[end]==';'
  a=replace(a,"\\\\"=>"\\")
  a=repr(MIME("text/plain"),eval(Meta.parse(a)),context=:limit=>true)
  if omit a="nothing" end
  a=replace(a,r" *(\n|$)"s=>s"\1")
  a=replace(a,r"\n$"s=>"")
  b=replace(b,r" *(\n|$)"s=>s"\1")
  b=replace(b,r"\n$"s=>"")
  i=1
  while i<=lastindex(a) && i<=lastindex(b) && a[i]==b[i]
    i=nextind(a,i)
  end
  if a!=b print("exec=$(repr(a[i:end]))\nmanl=$(repr(b[i:end]))\n") end
  a==b
end
@testset verbose = true "Gapjm" begin
@testset "SPerms.jl" begin
@test mytest("SPerms.jl","SPerm([-2,-1,-3])","SPerm{Int64}: (1,-2)(3,-3)")
@test mytest("SPerms.jl","p=SPerm(-1)","(1,-1)")
@test mytest("SPerms.jl","q=SPerm(1,2)","(1,2)")
@test mytest("SPerms.jl","sort(elements(Group([p,q])))","8-element Vector{SPerm{Int16}}:\n (1,-2)\n (1,-2,-1,2)\n (1,-1)(2,-2)\n (1,-1)\n (2,-2)\n ()\n (1,2,-1,-2)\n (1,2)")
@test mytest("SPerms.jl","SPerm([-2,-1,-3])==SPerm([-2,-1,-3,4])","true")
@test mytest("SPerms.jl","cycles(SPerm(-1,2)*SPerm(3,-3)*SPerm(4,5,-4,-5))","3-element Vector{Vector{Int16}}:\n [1, -2]\n [3, -3]\n [4, 5, -4, -5]")
@test mytest("SPerms.jl","cycletype(SPerm(1,-1),2)","2-element Vector{Vector{Int64}}:\n [1]\n [1]")
@test mytest("SPerms.jl","p=SPerm([-2,-1,-3])","SPerm{Int64}: (1,-2)(3,-3)")
@test mytest("SPerms.jl","permute([20,30,40],p)","3-element Vector{Int64}:\n -30\n -20\n -40")
@test mytest("SPerms.jl","p=SPerm([20,30,40],[-40,-20,-30])","(1,-2,3,-1,2,-3)")
@test mytest("SPerms.jl","permute([20,30,40],p)","3-element Vector{Int64}:\n -40\n -20\n -30")
@test mytest("SPerms.jl","Matrix(SPerm([-2,-1,-3]))","3×3 Matrix{Int64}:\n  0  -1   0\n -1   0   0\n  0   0  -1")
@test mytest("SPerms.jl","elements(CoxHyperoctaedral(2))","8-element Vector{SPerm{Int8}}:\n ()\n (1,2)\n (1,-1)\n (1,2,-1,-2)\n (1,-2,-1,2)\n (2,-2)\n (1,-2)\n (1,-1)(2,-2)")
@test mytest("SPerms.jl","uc=UnipotentCharacters(complex_reflection_group(6));","nothing")
@test mytest("SPerms.jl","g=sstab_onmats(fourier(uc.families[2]))","Group([(1,18)(3,-6)(8,-21)(10,-16)(11,22)(13,15),(1,-15)(2,-19)(3,-11)(6,22)(7,-12)(13,-18),(2,19)(4,-14)(5,20)(7,12),(1,-11)(2,-19)(3,-15)(5,-20)(6,13)(8,10)(16,21)(17,-17)(18,-22),(1,-22)(2,-19)(3,-13)(5,-20)(6,15)(8,-16)(10,-21)(11,-18)(17,-17),(1,6)(2,-19)(3,-18)(4,14)(8,16)(9,-9)(10,21)(11,-13)(15,-22),(1,13)(3,22)(4,14)(5,-20)(6,-11)(8,21)(9,-9)(10,16)(15,18)(17,-17)])")
@test mytest("SPerms.jl","length(g)","32")
@test mytest("SPerms.jl","f=SubFamilyij(chevieget(:families,:X)(12),1,3,(3+root(-3))//2);","nothing")
@test mytest("SPerms.jl","M=fourier(conj(f));","nothing")
@test mytest("SPerms.jl","uc=UnipotentCharacters(complex_reflection_group(6));","nothing")
@test mytest("SPerms.jl","N=fourier(uc.families[2]);","nothing")
@test mytest("SPerms.jl","p=SPerm_onmats(M,N)","(1,3)(2,19,-2,-19)(4,-14,-4,14)(5,-5)(6,-18)(7,-7)(8,10)(11,15,-11,-15)(12,-12)(13,22)(16,21,-16,-21)")
@test mytest("SPerms.jl","^(M,p;dims=(1,2))==N","true")
end
end
