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

To update to the latest version, I do not know a better way that to rm
the package and re-install it.
"""
module Gapjm
using Reexport

#--------------------------------------------------------------------------
export degree, degrees, elements, words, word, root
function degree end
function degrees end
function elements end
function root end
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
include("PermRoot.jl")
@reexport using .PermRoot
include("CoxGroups.jl")
@reexport using .CoxGroups
include("Garside.jl")
@reexport using .Garside
include("Weyl.jl")
@reexport using .Weyl
include("Hecke.jl")
@reexport using .Hecke
include("KL.jl")
@reexport using .KL
include("CycPols.jl")
@reexport using .CycPols
include("Symbols.jl")
@reexport using .Symbols
include("HasType.jl")
@reexport using .HasType
end
