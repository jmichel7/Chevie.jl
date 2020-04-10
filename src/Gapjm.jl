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

-  infrastructure
Permutations,  cyclotomic  numbers,  Laurent  polynomials.  There  are also
permutation  groups, for which I have  often replaced the proper algorithms
of  GAP by naive but  easy to write methods  only suitable for small groups
(sufficient  for the rest of the package but maybe not for your needs). The
code  for infrastructure is often competitive  with GAP, despite being much
shorter (often 100 lines of Julia replace 1000 lines of C); I am sure there
are more optimisations possible. Any comments about the code and the design
are welcome.

-  ported from Chevie
about 75% of Chevie functionality. The function `gap` can help you discover
the  equivalent functionality  to a  Gap3 function:  it takes  a string and
gives you Julia translations of functions in Gap3 which match this string:

```julia-rep1
julia> gap("words")
CoxeterWords(W[,l])      =>  word.(Ref(W),elements(W[,l])
GarsideWords             =>  elements
CharRepresentationWords  =>  traces_words_mats
```
Then you can call on-line help on the discovered functions.

This package is  often 10 times faster  than the equivalent GAP3 Chevie code
(after the maddeningly long compilation time on first execution).

I tried that parts of my package can be used independently of the rest. For
instance,  the modules `Combinat`,  `Groups`, `ModuleElts`, `Perms`, `Util`
are  independent of  the rest  of the  package and  can be used stand-alone
(after a trivial change forced by the exporting rules of Julia, see below).
"""
module Gapjm
using Reexport

#--------------------------------------------------------------------------
# The glue code below should not exist: extending e.g. Gapjm.degree instead
# of just exporting degree from Pols and Perms is a horrible hack forced by
# the unpleasant Julia rules not merging methods from different modules.
# 
# This  way both Pols  and Perms can  work together only  as part of Gapjm.
# Otherwise they could not work together.
# 
# If  degree, etc... was in Base there  would be no problem, both importing
# from Base.
function coefficients end; export coefficients
degree(a::Number)=0; export degree
function degrees end; export degrees
function elements end; export elements
function kernel end; export kernel
function order end; export order
function restricted end; export restricted
function root end; export root
function roots end; export roots
function valuation end; export valuation
function words end; export words
function word end; export word

include("Util.jl");@reexport using .Util
include("ModuleElts.jl");@reexport using .ModuleElts
include("Groups.jl");@reexport using .Groups
include("Perms.jl");@reexport using .Perms
include("PermGroups.jl");@reexport using .PermGroups
include("Cycs.jl");@reexport using .Cycs
include("Combinat.jl");@reexport using .Combinat
include("Pols.jl");@reexport using .Pols
include("Mvps.jl");@reexport using .Mvps
Mvps.Mvp(p::Pol)=p(Mvp(Pols.varname[]))
include("PermRoot.jl");@reexport using .PermRoot
include("GLinearAlgebra.jl");@reexport using .GLinearAlgebra
include("CoxGroups.jl");@reexport using .CoxGroups
include("Weyl.jl");@reexport using .Weyl
include("Cosets.jl");@reexport using .Cosets
include("Chars.jl");@reexport using .Chars
include("SPerms.jl");@reexport using .SPerms
include("HeckeAlgebras.jl");@reexport using .HeckeAlgebras
include("KL.jl");@reexport using .KL
include("CycPols.jl");@reexport using .CycPols
include("Symbols.jl");@reexport using .Symbols
include("Garside.jl");@reexport using .Garside
include("Posets.jl");@reexport using .Posets
include("MatInt.jl");@reexport using .MatInt
include("Chevie.jl");@reexport using .Chevie
include("Families.jl");@reexport using .Families
include("Uch.jl");@reexport using .Uch
include("ComplexR.jl");@reexport using .ComplexR
include("HasType.jl");@reexport using .HasType
include("Ucl.jl");@reexport using .Ucl
include("Eigenspaces.jl");@reexport using .Eigenspaces
include("Lusztig.jl");@reexport using .Lusztig
include("../docs/src/cheviedict.jl");export gap
end
