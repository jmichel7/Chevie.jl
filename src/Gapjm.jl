"""
Here are some of my efforts porting GAP code to julia.
I am not even sure this is a well-formed Julia package.

It contains for now permutations and permutation groups, cyclotomic numbers
and  Laurent polynomials. Coming  soon are Coxeter  groups, Hecke algebras,
braid groups and Garside monoids.

Even  though the code  is often competitive  with or faster  than GAP, I am
sure there are more optimisations possible. Any comments about the code and
the design are welcome.
"""
module Gapjm
using Reexport
push!(LOAD_PATH, @__DIR__) # not sure if this is the right way to use modules
@reexport using Perms
@reexport using PermGroups
@reexport using Cycs
@reexport using Util
@reexport using Pols
end
