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

#--------------------------------------------------------------------------
export degree, gens, elements, words
function degree end
function gens end
function elements end
function words end
"""
This  should  not  exist:  extending  e.g.  Gapjm.degree  instead  of  just
exporting  degree from  Pols and  Perms is  a horrible  hack forced  by the
unpleasant Julia rules not merging methods from different modules.

This  way both Pols and Perms can work together only as part of Gapjm.
Otherwise they could not work together.

It degree was in Base there would be no problem, both importing from Base.
"""

include("Util.jl")
include("Perms.jl")
include("PermGroups.jl")
include("Cycs.jl")
include("Pols.jl")

@reexport using .Util
@reexport using .Perms
@reexport using .PermGroups
@reexport using .Cycs
@reexport using .Pols

end
