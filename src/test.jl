module Test

foo(a)=1

using Perms: Perms

println(names(@__MODULE__))
println(names(Perms))

end
