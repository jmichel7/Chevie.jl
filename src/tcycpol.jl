module Test
using Reexport
include("Util.jl");@reexport using .Util
include("Combinat.jl");@reexport using .Combinat
include("CycPols.jl");@reexport using .CycPols
end
using .Test
using LaurentPolynomials
using CyclotomicNumbers
using BenchmarkTools
#julia> @btime CyclotomicNumbers.testmat(12)^2
#  305.887 ms (4068528 allocations: 338.53 MiB)
