includet("CycPols.jl");using .CycPols
using LaurentPolynomials
using BenchmarkTools: @btime
#julia> @btime CyclotomicNumbers.testmat(12)^2
#  305.887 ms (4068528 allocations: 338.53 MiB)
