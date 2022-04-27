if false #timing differs dramatically whether true/false
include("../src/Combinat.jl");using .Combinat: arrangements
include("../src/Perms.jl");using .Perms: Perm
else
using Gapjm
end
using BenchmarkTools
v=map(x->x=>1,Perm{Int16}.(arrangements(1:4,4)))
@btime sort($v,by=first)
