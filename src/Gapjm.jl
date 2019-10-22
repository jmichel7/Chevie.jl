"""
This  is  my  effort  porting  GAP  code  to Julia, specifically the Chevie
package  of GAP3 plus the minimal other GAP functionality needed for Chevie
to   work:  Cyclotomics,  Permutations,   Laurent  polynomials,  and  basic
permutation group operations.

I am rather new to Julia, git and github so I am not even sure this package
is  properly constituted; I did not try yet to register it. If you are more
competent  that me and see anything to  be improved in this package, please
write me or make a pull request.

### Installing

To install this package, at the Julia command line:

  *  enter package mode with ]
  *  do the command
```
(v1.0) pkg> add "https://github.com/jmichel7/Gapjm.jl"
```
- exit package mode with backspace and then do 
```
julia> using Gapjm
```
and you are set up.

To update later to the latest version, do

```
(v1.0) pkg> update "https://github.com/jmichel7/Gapjm.jl"
```

The package currently contains:

-  infrastructure: permutations,  cyclotomic numbers,  Laurent polynomials.
There  are also  permutation groups,  for which  I have  often replaced the
proper  algorithms of GAP by naive but  easy to write methods only suitable
for  small groups (sufficient for the rest of the package but maybe not for
your needs).

-  ported from  Chevie:  Weyl  groups,  Coxeter  groups,  Hecke  algebras,
Kazhdan-Lusztig   polynomials,  braid  and   Garside  groups  and  monoids,
factorisations into cyclotomic polynomials, character tables of Weyl groups
and  Hecke algebras, Unipotent characters  of Spetses, unipotent classes of
reductive groups.

The  code for infrastructure  is often competitive  with GAP, despite being
much  shorter (often 100 lines of Julia replace 1000 lines of C); I am sure
there  are more optimisations possible. Any comments about the code and the
design  are welcome. The code for Chevie  is often 10 times faster than the
GAP3 Chevie (after the maddeningly long compilation time on first execution).
"""
module Gapjm
using Reexport

#--------------------------------------------------------------------------
export coefficients, degree, degrees, elements, kernel, restricted, root, words, word
function coefficients end
function degree end
function degrees end
function elements end
function kernel end
function restricted end
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

export toL, toM # convert Gap matrices <-> Julia matrices
toL(m)=[m[i,:] for i in axes(m,1)] # to Gap
toM(m)=isempty(m) ? permutedims(hcat(m...)) : permutedims(reduce(hcat,m))     # to julia
export ds # dump struct
function ds(s)
  println(typeof(s),":")
  for f in fieldnames(typeof(s))
    println(f,"=",getfield(s,f))
  end
end

include("Util.jl")
@reexport using .Util
include("ModuleElts.jl")
@reexport using .ModuleElts
include("Groups.jl")
@reexport using .Groups
include("Perms.jl")
@reexport using .Perms
include("SPerms.jl")
@reexport using .SPerms
include("PermGroups.jl")
@reexport using .PermGroups
include("Cycs.jl")
@reexport using .Cycs
include("Pols.jl")
@reexport using .Pols
include("Mvps.jl")
@reexport using .Mvps
include("PermRoot.jl")
@reexport using .PermRoot
include("Chars.jl")
@reexport using .Chars
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
include("Cosets.jl")
@reexport using .Cosets
include("Posets.jl")
@reexport using .Posets
include("MatInt.jl")
@reexport using .MatInt
include("HasType.jl")
@reexport using .HasType
include("Uch.jl")
@reexport using .Uch
include("Ucl.jl")
@reexport using .Ucl
end
