using ModuleElts
using BenchmarkTools
includet("Util.jl");using .Util
includet("Cycs.jl");using .Cycs
#julia> @btime Cycs.testmat(12)^2
#  305.887 ms (4068528 allocations: 338.53 MiB)
