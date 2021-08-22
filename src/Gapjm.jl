"""
This  is  my  effort  porting  GAP  code  to Julia, specifically the Chevie
package  of  GAP3  plus  the  GAP  functionality needed for Chevie to work:
Cyclotomics,   Permutations,   Laurent   and   Puiseux  polynomials,  basic
permutation group operations, etc….

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

The package currently contains as infrastructure:
  * permutations
  * cyclotomic numbers
  * univariate Laurent and multivariate Puiseux polynomials
  * combinatorics
  * linear algebra on any field/ring
  * posets
  * cyclotomic polynomials
  * signed permutations
  * finite fields
  * groups
  * permutation groups

for  permutation groups I have  often replaced the sophisticated algorithms
of  GAP by naive but  easy to write methods  only suitable for small groups
(sufficient  for the  rest of  the package  but maybe  not for your needs).
Otherwise  the  code  for  infrastructure  is  often  competitive with GAP,
despite  being much shorter (often 100 lines of Julia replace 1000 lines of
C); I am sure there are more optimisations possible. Any comments about the
code and the design are welcome.

This  package contains about 90% of Chevie functionality, ported from Gap3.
The  function `gap` can help you discover the equivalent functionality to a
Gap3  function:  it  takes  a  string  and  gives you Julia translations of
functions in Gap3 which match this string.

```julia-rep1
julia> gap("words")
CoxeterWords(W[,l])      =>  word.(Ref(W),elements(W[,l])
GarsideWords             =>  elements
CharRepresentationWords  =>  traces_words_mats
```
Then you can call on-line help on the discovered functions.

The data library of Chevie has been automatically ported by a transpiler so
its code is "strange". Otherwise the code in this package is often 10 times
faster  than the  equivalent GAP3  Chevie code  (after the maddeningly long
compilation time on first execution).

I  tried  that  as  any  submodules  as  possible in my package can be used
independently  of the rest, thus could be independent packages. This is the
case  for the  modules `Combinat`,  `Groups`, `ModuleElts`, `Perms`, `Util`
which  can  be  used  stand-alone.  In  addition  modules `MatInt`, `Cycs`,
`Pols`, `Mvp`, `Posets`, `FFields` can be used stand-alone except they each
use one or two functions from `Util`.
"""
module Gapjm
using Reexport
using Requires
using UsingMerge

#--------------------------------------------------------------------------
function degrees end; export degrees
function roots end; export roots

include("../docs/src/cheviedict.jl");export gap
include("Util.jl");@reexport using .Util
include("Groups.jl");@reexport using .Groups
include("Combinat.jl");@reexport using .Combinat
include("Perms.jl");@usingmerge verbose=true reexport Perms
include("Pols.jl");@reexport using .Pols
include("ModuleElts.jl");@reexport using .ModuleElts
include("Cycs.jl");@usingmerge verbose=true reexport Cycs
#if false
include("Mvps.jl");@usingmerge verbose=true reexport Mvps
include("Posets.jl");@usingmerge verbose=true reexport Posets
include("FFields.jl");@usingmerge verbose=true reexport FFields
include("FFfac.jl")
include("MatInt.jl");@reexport using .MatInt
include("PermGroups.jl");@reexport using .PermGroups
include("PermRoot.jl");@reexport using .PermRoot
include("CoxGroups.jl");@reexport using .CoxGroups
include("Weyl.jl");@reexport using .Weyl
include("Cosets.jl");@reexport using .Cosets
include("ComplexR.jl");@reexport using .ComplexR
include("Semisimple.jl");@reexport using .Semisimple
include("Chars.jl");@reexport using .Chars
include("GLinearAlgebra.jl");@reexport using .GLinearAlgebra
include("mvptools.jl")
include("SPerms.jl");@reexport using .SPerms
include("Algebras.jl");@reexport using .Algebras
include("Presentations.jl");@reexport using .Presentations
include("Garside.jl");@reexport using .Garside
include("Chevie.jl");@reexport using .Chevie
include("Urad.jl");@reexport using .Urad
include("Lusztig.jl");@reexport using .Lusztig
include("Eigenspaces.jl");@reexport using .Eigenspaces
include("CycPols.jl");@reexport using .CycPols
include("HeckeAlgebras.jl");@reexport using .HeckeAlgebras
include("KL.jl");@reexport using .KL
include("Symbols.jl");@reexport using .Symbols
include("Ucl.jl");@reexport using .Ucl
include("Gt.jl");@reexport using .Gt
include("Murphy.jl");@reexport using .Murphy
include("Families.jl");@reexport using .Families
include("Uch.jl");@reexport using .Uch
include("dSeries.jl");@reexport using .dSeries
include("Sscoset.jl");@reexport using .Sscoset
include("HasType.jl");@reexport using .HasType
#end
function __init__()
  @require GAP="c863536a-3901-11e9-33e7-d5cd0df7b904" include("Gap4.jl")
end
end
