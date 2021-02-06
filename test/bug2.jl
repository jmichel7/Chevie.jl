using Reexport, BenchmarkTools
if true # @btime sort(v,by=first) differs dramatically whether true/false
degree(a::Number)=0; export degree
include("../src/using_merge.jl")
include("../src/Combinat.jl");using .Combinat
include("../src/Perms.jl");using_merge(:Perms)
include("../src/Util.jl");using .Util
include("../src/Pols.jl");using_merge(:Pols)
else
using Gapjm
end
#v=[rand(Perm,10)=>Pol(rand(-10:10,rand(1:6)),rand(1:6)) for i in 1:30]
v=map(x->Perm{Int16}(x)=>Pol([1,order(Perm(x))]),arrangements(1:4,4))
