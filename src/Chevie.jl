"""
This is my attempt to port the Chevie package from GAP3 to Julia. I started
this  project at the end of  2018 and it is still  in flux so some function
names or interfaces may still change. Pull requests and issues are welcome.

I  have implemented the  GAP functionality (infrastructure)  needed to make
Chevie  work.  I  have  already  registered  most of this infrastructure as
separate  packages; the  following packages  are loaded  and re-exported so
that  their functionality is automatically available when you use `Chevie`.
In other words, `Chevie` is a meta-package for the following packages:

  * (univariate) [LaurentPolynomials](https://github.com/jmichel7/LaurentPolynomials.jl) (and rational fractions)
  * (multivariate) [PuiseuxPolynomials](https://github.com/jmichel7/PuiseuxPolynomials.jl) (and rational fractions when there are no fractional exponents)
  * [CyclotomicNumbers](https://github.com/jmichel7/CyclotomicNumbers.jl)(elements of cyclotomic fields)
  * [ModuleElts](https://github.com/jmichel7/ModuleElts.jl) (elements of a free module over some ring)
  * [Combinat](https://github.com/jmichel7/Combinat.jl) (combinatorics and some basic number theory)
  * [PermGroups](https://github.com/jmichel7/PermGroups.jl) (permutations, groups, permutations groups. It contains the modules `Perms` and `Groups` which could be separate packages)
  * [SignedPerms](https://github.com/jmichel7/SignedPerms.jl) (signed permutations)
  * [MatInt](https://github.com/jmichel7/MatInt.jl) (Integer matrices and lattices)
  * [CycPols](https://github.com/jmichel7/CycPols.jl) (cyclotomic polynomials)
  * [GenLinearAlgebra](https://github.com/jmichel7/GenLinearAlgebra.jl) (linear algebra on any field/ring)
  * [FinitePosets](https://github.com/jmichel7/FinitePosets.jl) (finite posets)
  * [FiniteFields](https://github.com/jmichel7/FiniteFields.jl) (finite fields)
  * [GroupPresentations](https://github.com/jmichel7/GroupPresentations.jl) (presentations of groups, and groups defined by generators and relations)
  * [UsingMerge](https://github.com/jmichel7/UsingMerge.jl) (Automatically compose several packages)

Have  a look at the  documentation of the above  packages to see how to use
their   features.  

I  have implemented  some other  infrastructure which  currently resides in
`Chevie` but may eventually become separate packages:

  * factorizing polynomials over finite fields (module [`FFfac`](@ref))
  * factorizing polynomials over the rationals (module [`Fact`](@ref))
  * Number fields which are subfields of the Cyclotomics (module [`Nf`](@ref))

For permutation groups I have often replaced GAP's sophisticated algorithms
with  naive  but  easy-to-write  methods  suitable  only  for  small groups
(sufficient  for the rest of  the package but perhaps  not for your needs).
Otherwise  the infrastructure code  is often competitive  with GAP, despite
using  much less code (often  100 lines of Julia  replace 1000 lines of C);
and I am sure it could be optimised better than I did. Comments on code and
design  are welcome. For functions that are too inefficient or difficult to
implement (such as character tables of arbitrary groups), `Chevie` uses the
`GAP`  package  as  an  extension.  This  means  that if you have the `GAP`
package  installed,  `Chevie`  will  automatically  call `GAP` to implement
these functions. 

Functions  in the  `Chevie.jl` package  are often  10 times faster than the
equivalent functions in GAP3/Chevie (after the maddeningly long compilation
time on the first run --- Julia's TTFP).

The  `Chevie`  package  currently  contains  about  95%  of the GAP3 Chevie
functionality.  If you  are a  user of  GAP3/Chevie, the `gap` function can
help  you to  find the  equivalent functionality  in `Chevie.jl`  to a Gap3
function:  it takes a string and  gives you Julia translations of functions
in Gap3 that match that string.

```julia-rep1
julia> gap("words")
CharRepresentationWords  =>  traces_words_mats
CoxeterWords(W[,l])      =>  word.(Ref(W),elements(W[,l]))
GarsideWords             =>  elements
```
You can then access online help for the functions you have found.

The  port to Julia is not complete in the sense that 80% of the code is the
data library from Chevie, which was automatically ported by a transpiler so
its  code is "strange".  When the need  to maintain the  `GAP3` and `Julia`
versions  simultaneously subsides,  I will  do a  proper translation of the
data library, which should give an additional speed boost.

`Chevie.jl` requires julia 1.8 or later. 
"""
module Chevie
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
@reexport using LinearAlgebra: diag, tr, I, Diagonal, exactdiv, det_bareiss, dot
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
include("Diagrams.jl");@reexport using .Diagrams
include("CoxGroups.jl");@reexport using .CoxGroups
include("Weyl.jl");@reexport using .Weyl
include("Cosets.jl");@reexport using .Cosets
include("ComplexR.jl");@reexport using .ComplexR
include("Chars.jl");@reexport using .Chars
include("Symbols.jl");@reexport using .Symbols
include("Tools2.jl");@reexport using .Tools2
include("Algebras.jl");@reexport using .Algebras
include("InitChevie.jl");@reexport using .InitChevie
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
end
