"""
Here are some of my efforts porting GAP code to julia.
I am not even sure this is a well-formed Julia package.

It  contains  for  now  permutations  and  permutation  groups,  cyclotomic
numbers,  Laurent polynomials, some  Weyl groups and  Coxeter groups, Hecke
algebras  and Kazhdan-Lusztig polynomials. Coming soon are braid groups and
Garside monoids and factorisations into cyclotomic polynomials.

Even  though the code  is often competitive  with or faster  than GAP, I am
sure there are more optimisations possible. Any comments about the code and
the design are welcome.

If you are new to julia, to install this package:
- enter package mode with ]
- do the command
```
(v1.0) pkg> add "https://github.com/jmichel7/Gapjm.jl"
```
- exit package mode with backspace
do 
```
julia> using Gapjm
```
and you are set up.
"""
module Gapjm
using Reexport

#--------------------------------------------------------------------------
export degree, degrees, gens, elements, words, word
function degree end
function degrees end
function gens end
function elements end
function words end
function word end
"""
This  should  not  exist:  extending  e.g.  Gapjm.degree  instead  of  just
exporting  degree from  Pols and  Perms is  a horrible  hack forced  by the
unpleasant Julia rules not merging methods from different modules.

This  way both Pols and Perms can work together only as part of Gapjm.
Otherwise they could not work together.

It degree was in Base there would be no problem, both importing from Base.
"""

include("Util.jl")
@reexport using .Util
include("Perms.jl")
@reexport using .Perms
include("PermGroups.jl")
@reexport using .PermGroups
include("Cycs.jl")
@reexport using .Cycs
include("Pols.jl")
@reexport using .Pols
include("CoxGroups.jl")
@reexport using .CoxGroups
include("PermRoot.jl")
@reexport using .PermRoot
include("Weyl.jl")
@reexport using .Weyl
include("Hecke.jl")
@reexport using .Hecke
include("KL.jl")
@reexport using .KL

end
