# auto-generated tests from julia-repl docstrings
using Test, Gapjm
#include("../tools/Gap4.jl")
function mytest(file::String,src::String,man::String)
  println(file," ",src)
  omit=src[end]==';'
  src=replace(src,"\\\\"=>"\\")
  exec=repr(MIME("text/plain"),eval(Meta.parse(src)),context=:limit=>true)
  if omit exec="nothing" end
  exec=replace(exec,r" *(\n|$)"s=>s"\1")
  exec=replace(exec,r"\n$"s=>"")
  man=replace(man,r" *(\n|$)"s=>s"\1")
  man=replace(man,r"\n$"s=>"")
  i=1
  while i<=lastindex(exec) && i<=lastindex(man) && exec[i]==man[i]
    i=nextind(exec,i)
  end
  if exec!=man 
    print("exec=$(repr(exec[i:end]))\nmanl=$(repr(man[i:end]))\n")
  end
  exec==man
end
@testset verbose = true "Gapjm" begin
@testset "Posets.jl" begin
@test mytest("Posets.jl","p=Poset([[2,3],[4],[4],Int[]])","1<2,3<4")
@test mytest("Posets.jl","summary(p)","\"Poset{Int64} with 6 elements\"")
@test mytest("Posets.jl","length(p)","4")
@test mytest("Posets.jl","incidence(p)","4×4 Matrix{Bool}:\n 1  1  1  1\n 0  1  0  1\n 0  0  1  1\n 0  0  0  1")
@test mytest("Posets.jl","l=vec(collect(Iterators.product(1:2,1:2)))","4-element Vector{Tuple{Int64, Int64}}:\n (1, 1)\n (2, 1)\n (1, 2)\n (2, 2)")
@test mytest("Posets.jl","P=Poset((x,y)->all(map(<=,x,y)),l)","(1, 1)<(2, 1),(1, 2)<(2, 2)")
@test mytest("Posets.jl","P=Poset([all(map(<=,x,y)) for x in l, y in l],l)","(1, 1)<(2, 1),(1, 2)<(2, 2)")
@test mytest("Posets.jl","P.show_element=(io,x,n)->join(io,x.elements[n],\".\");","nothing")
@test mytest("Posets.jl","P","1.1<2.1,1.2<2.2")
@test mytest("Posets.jl","print(P)","Poset([[2, 3], [4], [4], Int64[]],[(1, 1), (2, 1), (1, 2), (2, 2)])")
@test mytest("Posets.jl","p=Poset((i,j)->j%i==0,1:6)","1<5\n1<2<4\n1<3<6\n2<6")
@test mytest("Posets.jl","p[2:6]","2-element Vector{Int64}:\n 2\n 6")
@test mytest("Posets.jl","p[2:end]","3-element Vector{Int64}:\n 2\n 4\n 6")
@test mytest("Posets.jl","p[begin:6]","4-element Vector{Int64}:\n 1\n 2\n 3\n 6")
@test mytest("Posets.jl","p[2:3]","Int64[]")
@test mytest("Posets.jl","m=[j-i in [0,1] for i in 1:5, j in 1:5]","5×5 Matrix{Bool}:\n 1  1  0  0  0\n 0  1  1  0  0\n 0  0  1  1  0\n 0  0  0  1  1\n 0  0  0  0  1")
@test mytest("Posets.jl","transitive_closure(m)","5×5 Matrix{Bool}:\n 1  1  1  1  1\n 0  1  1  1  1\n 0  0  1  1  1\n 0  0  0  1  1\n 0  0  0  0  1")
@test mytest("Posets.jl","Poset(Bool[1 1 1 1 1;0 1 0 1 1;0 0 1 1 1;0 0 0 1 0;0 0 0 0 1])","1<2,3<4,5")
@test mytest("Posets.jl","Poset([[2,3],[4,5],[4,5],Int[],Int[]])","1<2,3<4,5")
@test mytest("Posets.jl","p=Poset((i,j)->j%i==0,1:6)","1<5\n1<2<4\n1<3<6\n2<6")
@test mytest("Posets.jl","linear_extension(p)","6-element Vector{Int64}:\n 1\n 2\n 3\n 5\n 4\n 6")
@test mytest("Posets.jl","p=Poset((i,j)->j%i==0,1:5)","1<3,5\n1<2<4")
@test mytest("Posets.jl","hasse(p)","5-element Vector{Vector{Int64}}:\n [2, 3, 5]\n [4]\n []\n []\n []")
@test mytest("Posets.jl","p=Poset([i==6 ? Int[] : [i+1] for i in 1:6])","1<2<3<4<5<6")
@test mytest("Posets.jl","incidence(p)","6×6 Matrix{Bool}:\n 1  1  1  1  1  1\n 0  1  1  1  1  1\n 0  0  1  1  1  1\n 0  0  0  1  1  1\n 0  0  0  0  1  1\n 0  0  0  0  0  1")
@test mytest("Posets.jl","p=Poset((i,j)->i%4<j%4,1:8)","4,8<1,5<2,6<3,7")
@test mytest("Posets.jl","reverse(p)","3,7<2,6<1,5<4,8")
@test mytest("Posets.jl","p=Poset([i==j || i%4<j%4 for i in 1:8, j in 1:8])","4,8<1,5<2,6<3,7")
@test mytest("Posets.jl","partition(p)","4-element Vector{Vector{Int64}}:\n [4, 8]\n [2, 6]\n [3, 7]\n [1, 5]")
@test mytest("Posets.jl","p=Poset((i,j)->i%4<j%4,1:8)","4,8<1,5<2,6<3,7")
@test mytest("Posets.jl","restricted(p,2:6)","3<4<1,5<2")
@test mytest("Posets.jl","restricted(p,2:6;show=:elements)","4<5<2,6<3")
@test mytest("Posets.jl","p=Poset((i,j)->j%i==0,1:8)","1<5,7\n1<2<4<8\n1<3<6\n2<6")
@test mytest("Posets.jl","isjoinlattice(p)","false")
@test mytest("Posets.jl","p=Poset((i,j)->j%i==0,1:8)","1<5,7\n1<2<4<8\n1<3<6\n2<6")
@test mytest("Posets.jl","ismeetlattice(p)","true")
end
end
