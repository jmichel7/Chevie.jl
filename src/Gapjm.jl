"""
This  is  my  effort  to  port  GAP  code to Julia, specifically the Chevie
package  from GAP3.  I started  this project  at the  end of 2018 and it is
still  in flux so  the package is  not yet registered.  If you see anything
that needs improving in this package, please contact me or make an issue or
pull request in the GitHub repository.

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

I  have implemented the  GAP functionality (infrastructure)  needed to make
Chevie  work.  I  have  already  registered  most of this infrastructure as
separate packages; the following packages are loaded and reexported so that
their  functionality is  automatically available  when you  use `Gapjm`. In
other terms, `Gapjm` is a meta-package for the following packages:

  * (univariate) [LaurentPolynomials](https://github.com/jmichel7/LaurentPolynomials.jl) (and rational fractions)
  * (multivariate) [PuiseuxPolynomials](https://github.com/jmichel7/PuiseuxPolynomials.jl) (and rational fractions when there are no fractional exponents)
  * [CyclotomicNumbers](https://github.com/jmichel7/CyclotomicNumbers.jl)
  * [ModuleElts](https://github.com/jmichel7/ModuleElts.jl) (elements of a free module over some ring)
  * [Combinat](https://github.com/jmichel7/Combinat.jl) (combinatorics and some basic number theory)
  * [PermGroups](https://github.com/jmichel7/PermGroups.jl) (permutations, groups, permutations groups. It contains the modules `Perms` and `Groups` which could be separate packages)
  * [SignedPerms](https://github.com/jmichel7/SignedPerms.jl) (signed permutations)
  * [MatInt](https://github.com/jmichel7/MatInt.jl) (Integer matrices and lattices)
  * [CycPols](https://github.com/jmichel7/CycPols.jl) (cyclotomic polynomials)
  * [GenLinearAlgebra](https://github.com/jmichel7/GenLinearAlgebra.jl) (linear algebra on any field/ring)
  * [FinitePosets](https://github.com/jmichel7/FinitePosets.jl) (finite posets)
  * [FiniteFields](https://github.com/jmichel7/FiniteFields.jl) (finite fields)
  * [UsingMerge](https://github.com/jmichel7/UsingMerge.jl) (Automatically compose several packages)
  * [GroupPresentations](https://github.com/jmichel7/GroupPresentations.jl) (presentations of groups, and groups defined by generators and relations)
Look  at the  documentation of  the above  packages to  see how  to use the
corresponding  features. I have implemented  some more infrastructure which
sits currently in `Gapjm` but may become eventually separate packages:
  * factorizing polynomials over finite fields (module [`FFfac`](@ref))
  * factorizing polynomials over the rationals (module [`Fact`](@ref))
  * Number fields which are subfields of the Cyclotomics (module [`Nf`](@ref))

for permutation groups I have often replaced GAP's sophisticated algorithms
with  naive  but  easy-to-write  methods  suitable  only  for  small groups
(sufficient  for the  rest of  the package  but maybe  not for your needs).
Otherwise  the infrastructure code  is often competitive  with GAP, despite
being much shorter (often 100 lines of Julia replace 1000 lines of C); I am
sure  there are more  optimisations possible. Any  comments on the code and
the design are welcome. For functions that are too inefficient or difficult
to  implement  (such  as  character  tables  of  arbitrary  groups) `Gapjm`
automatically calls GAP4 if you did `using GAP`. Otherwise the code in this
package  is  often  10  times  faster  than the equivalent GAP3 Chevie code
(after  the maddeningly  long compilation  time on  first execution --- the
TTFP problem of Julia).

The   `Gapjm`  package   currently  contains   about  90%   of  the  Chevie
functionality,  ported from Gap3.  The `gap` function  can help you to find
the  equivalent functionality  to a  Gap3 function:  it takes  a string and
gives you Julia translations of functions in Gap3 that match that string.

```julia-rep1
julia> gap("words")
CharRepresentationWords  =>  traces_words_mats
CoxeterWords(W[,l])      =>  word.(Ref(W),elements(W[,l]))
GarsideWords             =>  elements
```
Then you can call on-line help on the discovered functions.

The  port to Julia is not complete in the sense that 80% of the code is the
data library from Chevie, which was automatically ported by a transpiler so
its code is "strange". When the need to maintain the `GAP3` version and the
`Julia`  version simultaneously subsides, I will do a proper translation of
the data library, which should give an additional speed boost.
"""
module Gapjm
#--------------------- external packages ----------------------------------
using Reexport
using UsingMerge
@reexport using PuiseuxPolynomials # reexports LaurentPolynomials
@reexport using LaurentPolynomials: stringexp
@reexport using ModuleElts
@reexport using Combinat
@reexport using MatInt
@reexport using Primes: factor, eachfactor, divisors
@reexport using OrderedCollections: OrderedDict
# careful: use very little of LinearAlgebra
@reexport using LinearAlgebra: diag, tr, I, Diagonal, exactdiv, det_bareiss
@reexport using PermGroups
@reexport using SignedPerms
@usingmerge verbose=true reexport CyclotomicNumbers
@reexport using CyclotomicNumbers: bracket_if_needed, 
  format_coefficient, stringind
@usingmerge verbose=true reexport FiniteFields
@reexport using CycPols
@reexport using CycPols: stringprime
@reexport using GenLinearAlgebra
@reexport using FinitePosets
@reexport using GroupPresentations
#--------------------- internal modules -----------------------------------
include("../docs/src/cheviedict.jl");export gap
include("Util.jl");@reexport using .Util
include("FFfac.jl");@reexport using .FFfac
include("Nf.jl");@reexport using .Nf
include("Tools.jl");@reexport using .Tools
include("Fact.jl");@reexport using .Fact
include("PermRoot.jl");@reexport using .PermRoot
include("CoxGroups.jl");@reexport using .CoxGroups
include("Weyl.jl");@reexport using .Weyl
include("Cosets.jl");@reexport using .Cosets
include("ComplexR.jl");@reexport using .ComplexR
include("Chars.jl");@reexport using .Chars
include("Symbols.jl");@reexport using .Symbols
include("Tools2.jl");@reexport using .Tools2
include("Algebras.jl");@reexport using .Algebras
include("Chevie.jl");@reexport using .Chevie
include("Lusztig.jl");@reexport using .Lusztig
include("Eigenspaces.jl");@reexport using .Eigenspaces
include("Garside.jl");@reexport using .Garside
include("HeckeAlgebras.jl");@reexport using .HeckeAlgebras
include("KL.jl");@reexport using .KL
include("Semisimple.jl");@reexport using .Semisimple
include("Urad.jl");@reexport using .Urad
include("Ucl.jl");@reexport using .Ucl
include("Gt.jl");@reexport using .Gt
include("Murphy.jl");@reexport using .Murphy
include("Families.jl");@reexport using .Families
include("Uch.jl");@reexport using .Uch
include("dSeries.jl");@reexport using .dSeries
include("Sscoset.jl");@reexport using .Sscoset
include("gendec.jl"); # for now no module
include("GAPENV.jl");@reexport using .GAPENV
function contr(s) 
  include(replace(@__DIR__,"/src"=>"/contr/"*s*".jl"))
end
export contr
roundtrip(x)=x==eval(Meta.parse(repr(x))) # for debugging purposes
using Requires
function __init__()
  @require GAP="c863536a-3901-11e9-33e7-d5cd0df7b904" (include("../ext/Gap4.jl");using .Gap4)
end
end
