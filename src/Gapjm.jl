"""
This  is  my  effort  porting  GAP  code  to Julia, specifically the Chevie
package  of GAP3. I started this project at the end of 2018 and it is still
in  flux so the  package is not  yet registered. If  you see anything to be
improved  in this  package, please  contact me  or make  an issue or a pull
request in the GitHub repository.

### Installing

To install this package, at the Julia command line:

  *  enter package mode with ]
  *  do the command
```
(@v1.7) pkg> add "https://github.com/jmichel7/Gapjm.jl"
```
- exit package mode with backspace and then do
```
julia> using Gapjm
```
and you are set up.

To update later to the latest version, do

```
(@v1.7) pkg> update Gapjm
```

This package requires julia 1.6 or later. 

I  also  implemented  the  GAP  functionality  needed  for  Chevie to work:
Cyclotomics,   Permutations,   Laurent   and   Puiseux  polynomials,  basic
permutation  group  operations,  etc….  I  registered  as separate packages
already  most of the infrastructure (the  following packages are loaded and
reexported  so their functionality is  automatically available when you use
`Gapjm`):

  * (univariate) [LaurentPolynomials](https://github.com/jmichel7/LaurentPolynomials.jl) (and rational fractions)
  * (multivariate) [PuiseuxPolynomials](https://github.com/jmichel7/PuiseuxPolynomials.jl) (and rational fractions)
  * [CyclotomicNumbers](https://github.com/jmichel7/CyclotomicNumbers.jl)
  * [ModuleElts](https://github.com/jmichel7/ModuleElts.jl) (elements of a free module over some ring)
  * [Combinat](https://github.com/jmichel7/Combinat.jl) (combinatorics and some basic number theory)
  * [PermGroups](https://github.com/jmichel7/PermGroups.jl) (permutations, groups, permutations groups. It contains the moduls `Perms` and `Groups` which could be separate packages)

some other infrastructure which may become eventually separate packages:
  * linear algebra on any field/ring (module `GLinearAlgebra`)
  * posets (module `Posets`)
  * cyclotomic polynomials (module `CycPols`)
  * signed permutations (module `SPerms`)
  * finite fields (module `FFields`)
  * Integer matrices and lattices (module `MatInt`)

for  permutation groups I have  often replaced the sophisticated algorithms
of  GAP by naive but  easy to write methods  only suitable for small groups
(sufficient  for the  rest of  the package  but maybe  not for your needs).
Otherwise  the  code  for  infrastructure  is  often  competitive with GAP,
despite  being much shorter (often 100 lines of Julia replace 1000 lines of
C); I am sure there are more optimisations possible. Any comments about the
code and the design are welcome. For functions which are too inefficient or
difficult  to implement (like character tables of arbitrary groups) `Gapjm`
automatically calls GAP4 if the package `GAP` is being used.

This  package `Gapjm` contains currently about 90% of Chevie functionality,
ported  from Gap3. The function `gap`  can help you discover the equivalent
functionality  to a Gap3  function: it takes  a string and  gives you Julia
translations of functions in Gap3 which match this string.

```julia-rep1
julia> gap("words")
CharRepresentationWords  =>  traces_words_mats
CoxeterWords(W[,l])      =>  word.(Ref(W),elements(W[,l]))
GarsideWords             =>  elements
```
Then you can call on-line help on the discovered functions.

The  port to Julia is not complete in the sense that 80% of the code is the
data library of Chevie, which has been automatically ported by a transpiler
so  its  code  is  "strange".  When  the  need  to  maintain  both versions
simultaneously subsides, it will be worth to do a proper translation of the
data library.

Otherwise  the  code  in  this  package  is  often 10 times faster than the
equivalent GAP3 Chevie code (after the maddeningly long compilation time on
first execution --- the TTFP problem of Julia).
"""
module Gapjm
#--------------------- external packages ----------------------------------
using Reexport
using Requires
using UsingMerge
@reexport using PuiseuxPolynomials # reexports LaurentPolynomials
@reexport using ModuleElts
@reexport using Combinat
@reexport using Primes: factor, eachfactor
@reexport using PermGroups
@usingmerge verbose=true reexport CyclotomicNumbers
#--------------------- internal modules -----------------------------------
include("../docs/src/cheviedict.jl");export gap
include("Util.jl");@reexport using .Util
include("Posets.jl");@usingmerge verbose=true reexport Posets
include("FFields.jl");@usingmerge verbose=true reexport FFields
include("FFfac.jl");@reexport using .FFfac
include("Fact.jl");@reexport using .Fact
include("Nf.jl");@reexport using .Nf
include("MatInt.jl");@reexport using .MatInt
include("Tools.jl");@reexport using .Tools
include("CycPols.jl");@reexport using .CycPols
include("PermRoot.jl");@reexport using .PermRoot
include("CoxGroups.jl");@reexport using .CoxGroups
include("Weyl.jl");@reexport using .Weyl
include("Cosets.jl");@reexport using .Cosets
include("ComplexR.jl");@reexport using .ComplexR
include("Semisimple.jl");@reexport using .Semisimple
include("Chars.jl");@reexport using .Chars
include("GLinearAlgebra.jl");@reexport using .GLinearAlgebra
include("Symbols.jl");@reexport using .Symbols
include("Tools2.jl");@reexport using .Tools2
include("SPerms.jl");@reexport using .SPerms
include("Algebras.jl");@reexport using .Algebras
include("Presentations.jl");@reexport using .Presentations
include("Garside.jl");@reexport using .Garside
include("Chevie.jl");@reexport using .Chevie
include("Urad.jl");@reexport using .Urad
include("Lusztig.jl");@reexport using .Lusztig
include("Eigenspaces.jl");@reexport using .Eigenspaces
include("HeckeAlgebras.jl");@reexport using .HeckeAlgebras
include("KL.jl");@reexport using .KL
include("Ucl.jl");@reexport using .Ucl
include("Gt.jl");@reexport using .Gt
include("Murphy.jl");@reexport using .Murphy
include("Families.jl");@reexport using .Families
include("Uch.jl");@reexport using .Uch
include("dSeries.jl");@reexport using .dSeries
include("Sscoset.jl");@reexport using .Sscoset
include("gendec.jl"); # for now no module
include("GAPENV.jl");@reexport using .GAPENV
function __init__()
  @require GAP="c863536a-3901-11e9-33e7-d5cd0df7b904" include("Gap4.jl")
end
end
