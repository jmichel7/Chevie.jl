# auto-generated tests from julia-repl docstrings
using Test, Gapjm
#include("../tools/Gap4.jl")
function mytest(a::String,b::String)
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
@testset "Urad.jl" begin
@test mytest("W=coxgroup(:E,6)","E₆")
@test mytest("U=UnipotentGroup(W)","UnipotentGroup(E₆)")
@test mytest("U(2=>4)","u2(4)")
@test mytest("U(2)^4","u2(4)")
@test mytest("U(2=>4)*U(4=>5)","u2(4)u4(5)")
@test mytest("U(2=>4,4=>5)","u2(4)u4(5)")
@test mytest("U(4=>5,2=>4)","u2(4)u4(5)u8(-20)")
@test mytest("W=coxgroup(:E,8);U=UnipotentGroup(W)","UnipotentGroup(E₈)")
@test mytest("u=U(map(i->i=>Z(2)*Mvp(Symbol(\"x_\",i)),1:8)...)","u1(1₂x₁)u2(1₂x₂)u3(1₂x₃)u4(1₂x₄)u5(1₂x₅)u6(1₂x₆)u7(1₂x₇)u8(1₂x₈)")
@test mytest("u^32","()")
@test mytest("W=coxgroup(:G,2)","G₂")
@test mytest("U=UnipotentGroup(W)","UnipotentGroup(G₂)")
@test mytest("u=U(1=>Mvp(:x),3=>Mvp(:y))","u1(x)u3(y)")
@test mytest("u^W(2,1)","u4(y)u5(x)")
@test mytest("s=SemisimpleElement(W,[E(3),2])","SemisimpleElement{Cyc{Int64}}: <ζ₃,2>")
@test mytest("u^s","u1(ζ₃x)u3(2ζ₃y)")
@test mytest("u^U(2)","u1(x)u3(x+y)u4(-x-2y)u5(x+3y)u6(x²+3xy+3y²)")
@test mytest("W=coxgroup(:G,2)","G₂")
@test mytest("U=UnipotentGroup(W)","UnipotentGroup(G₂)")
@test mytest("U.specialPairs","10-element Array{Array{Int64,1},1}:\n [1, 2, 3]\n [2, 3, 4]\n [2, 4, 5]\n [1, 5, 6]\n [3, 4, 6]\n [2, 1, 3]\n [3, 2, 4]\n [4, 2, 5]\n [5, 1, 6]\n [4, 3, 6]")
@test mytest("U.N","10-element Array{Int64,1}:\n  1\n  2\n  3\n  1\n  3\n -1\n -2\n -3\n -1\n -3")
@test mytest("U.commutatorConstants","10-element Array{Array{Array{Int64,1},1},1}:\n [[1, 1, 3, 1], [1, 2, 4, -1], [1, 3, 5, 1], [2, 3, 6, 2]]\n [[1, 1, 4, 2], [2, 1, 5, 3], [1, 2, 6, -3]]\n [[1, 1, 5, 3]]\n [[1, 1, 6, 1]]\n [[1, 1, 6, 3]]\n [[1, 1, 3, -1], [2, 1, 4, -1], [3, 1, 5, -1], [3, 2, 6, -1]]\n [[1, 1, 4, -2], [2, 1, 6, -3], [1, 2, 5, 3]]\n [[1, 1, 5, -3]]\n [[1, 1, 6, -1]]\n [[1, 1, 6, -3]]")
@test mytest("U=UnipotentGroup(coxgroup(:G,2))","UnipotentGroup(G₂)")
@test mytest("l=norm(U,[2=>4,1=>2])","6-element Array{Pair{Int64,Int64},1}:\n 1 => 2\n 2 => 4\n 3 => -8\n 4 => 32\n 5 => -128\n 6 => 512")
@test mytest("norm(U,l,6:-1:1)","2-element Array{Pair{Int64,Int64},1}:\n 2 => 4\n 1 => 2")
@test mytest("U=UnipotentGroup(coxgroup(:G,2))","UnipotentGroup(G₂)")
@test mytest("U(2)","u2(1)")
@test mytest("U(1=>2,2=>4)","u1(2)u2(4)")
@test mytest("U(2=>4,1=>2)","u1(2)u2(4)u3(-8)u4(32)u5(-128)u6(512)")
@test mytest("U=UnipotentGroup(coxgroup(:G,2))","UnipotentGroup(G₂)")
@test mytest("u=U(2=>Mvp(:y),1=>Mvp(:x))","u1(x)u2(y)u3(-xy)u4(xy²)u5(-xy³)u6(2x²y³)")
@test mytest("abelianpart(u)","u1(x)u2(y)")
@test mytest("W=coxgroup(:G,2)","G₂")
@test mytest("U=UnipotentGroup(W)","UnipotentGroup(G₂)")
@test mytest("u=U(2=>Mvp(:y),1=>Mvp(:x))","u1(x)u2(y)u3(-xy)u4(xy²)u5(-xy³)u6(2x²y³)")
@test mytest("decompose(W(1),u)","2-element Array{UnipotentElement{Mvp{Int64,Int64}},1}:\n u1(x)\n u2(y)u3(-xy)u4(xy²)u5(-xy³)u6(2x²y³)")
@test mytest("decompose(W(2),u)","2-element Array{UnipotentElement{Mvp{Int64,Int64}},1}:\n u2(y)\n u1(x)")
end
