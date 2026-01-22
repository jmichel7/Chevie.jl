"""
This  is the current version of the `Chevie` package. It began life in 2018
as  a port  to Julia  of the  `GAP3` package  with the  same name.  All new
developments  are  now  carried  out  in  this  version.  If  you are using
`Chevie.jl`   in  a  paper  and  wish  to  acknowledge  it,  you  can  cite
[mi15](@cite).

The  package has  no yet  reached version  1, thus  some function  names or
interfaces may yet change. Pull requests and issues are welcomed.

The  infrastructure needed to make `Chevie` work has also been implemented,
most  of it  as independent  registered packages.  They may have advantages
compared  to  other  Julia  packages  providing  similar functionality. You
should  take a look at them. They are loaded and re-exported by `Chevie` so
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
  * [UsingMerge](https://github.com/jmichel7/UsingMerge.jl) (Automatically compose a package with the current environment)

Take  a  look  at  the  documentation  for  the  above packages to see what
features they provide and how to use them.

Some  other implemented infrastructure which  currently resides in `Chevie`
but may eventually become separate packages:

  * factorizing polynomials over finite fields (module `FFfac`)
  * factorizing polynomials over the rationals (module `Fact`)
  * Number fields which are subfields of the Cyclotomics (module [`Nf`](@ref))
  * Truncated Laurent series (module [`Truncs`](@ref))

For permutation groups I have often replaced GAP's sophisticated algorithms
with  naive  but  easy-to-write  methods  that  are only suitable for small
groups.  These are sufficient for  the rest of the  package but perhaps not
for your needs. Otherwise the infrastructure code is often competitive with
`GAP3`'s, despite using much less code --- often 100 lines of Julia replace
1000  lines of `C`; and I am sure  it could be optimised better than I did.
Comments  on  code  and  design  are  welcome.  For  functions that are too
inefficient  or  difficult  to  implement  (such  as  character  tables  of
arbitrary  groups), `Chevie` uses  the `GAP` package  as an extension. This
means  that  if  you  have  the  `GAP`  package  installed,  `Chevie`  will
automatically call `GAP` to implement these functions.

The functions in the `Chevie.jl` package are often 10 times faster than the
equivalent   functions  in  `GAP3/Chevie`,  after  the  frustratingly  long
compilation time on the first run (Julia's TTFP).

The  `Chevie` package currently implements  almost all of the `GAP3/Chevie`
functionality,  as well as some functionality from the `GAP3` `Algebra` and
`VKcurve`  packages.  It  also  has  some  new functionality not present in
`GAP3/Chevie`.  If you are a user  of `GAP3/Chevie`, the `gap` function can
help  you to find the equivalent  functionality in `Chevie.jl` for a `GAP3`
function:  it takes  a string  and provides  Julia translations of matching
`GAP3` functions (the match is case-insensitive).

```julia-rep1
julia> gap("words")
CharRepresentationWords  =>  traces_words_mats
CoxeterWords(W[,l])      =>  word.(Ref(W),elements(W[,l]))
GarsideWords             =>  elements
```
You can then access online help for the functions you have found.

### Installing

`Chevie`  is a  registered package  that can  be installed/upgraded  in the
standard  way. For  Julia newbies,  here's a  reminder of  what this is. To
install, enter the following at the REPL command line:

  *  enter package mode with ]
  *  do the command
```
(@v1.10) pkg> add Chevie
```
- exit package mode with backspace and then do
```
julia> using Chevie
```
and you are set up. For first help, type "?Chevie". Have a look at the
[documentation](https://jmichel7.github.io/Chevie.jl/).

To update later to the latest version, do

```
(@v1.10) pkg> update
```
`Chevie.jl` requires julia 1.10 or later. 
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
@reexport using LinearAlgebra: LinearAlgebra, diag, tr, I, Diagonal, exactdiv,
 eigen, det_bareiss, dot
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
include("Format.jl");@reexport using .Format
include("Util.jl");@reexport using .Util
include("FFfac.jl");@reexport using .FFfac
include("Nf.jl");@reexport using .Nf
include("Truncs.jl");@reexport using .Truncs
include("Tools.jl");@reexport using .Tools
include("Fact.jl");@reexport using .Fact
include("Symbols.jl");@reexport using .Symbols
include("InitChevie.jl");@reexport using .InitChevie
include("PermRoot.jl");@reexport using .PermRoot
include("Diagrams.jl");@reexport using .Diagrams
include("CoxGroups.jl");@reexport using .CoxGroups
include("Weyl.jl");@reexport using .Weyl
include("Rootdata.jl");@reexport using .Rootdata
include("Cosets.jl");@reexport using .Cosets
include("Semisimple.jl");@reexport using .Semisimple
include("Chars.jl");@reexport using .Chars
include("Tools2.jl");@reexport using .Tools2
include("Algebras.jl");@reexport using .Algebras
include("Lusztig.jl");@reexport using .Lusztig
include("Eigenspaces.jl");@reexport using .Eigenspaces
include("Garside.jl");@reexport using .Garside
include("HeckeAlgebras.jl");@reexport using .HeckeAlgebras
include("KL.jl");@reexport using .KL
include("Urad.jl");@reexport using .Urad
include("Ucl.jl");@reexport using .Ucl
include("Gt.jl");@reexport using .Gt
include("Murphy.jl");@reexport using .Murphy
include("Families.jl");@reexport using .Families
include("Uch.jl");@reexport using .Uch
include("dSeries.jl");@reexport using .dSeries
include("Sscoset.jl");@reexport using .Sscoset
include("gendec.jl"); # for now no module
include("tbl/util.jl")
const cheviedata=[("G(de,e,n)","cmplximp"),
  ("Aₙ","weyla"), ("Bₙ and Cₙ","weylbc"),("Dₙ","weyld"),
  ("³D₄","weyl3d4"), ("²Aₙ","weyl2a"), ("²Dₙ","weyl2d"),
  ("ᵗG(e,e,n)","cmpxtimp"), ("G₂","weylg2"), 
  ("I₂(e)","coxi"), ("²I₂(e)","cox2i"),
  ("G₄-G₂₂","cmp4_22"),("H₃","coxh3"),
  ("G₂₄","cmplxg24"),("G₂₅","cmplxg25"),("G₂₆","cmplxg26"),
  ("G₂₇","cmplxg27"), ("F₄","weylf4"), ("²F₄","weyl2f4"),
  ("G₂₉","cmplxg29"), ("H₄","coxh4"), ("G₃₁","cmplxg31"),
  ("G₃₂","cmplxg32"),("G₃₃","cmplxg33"),("G₃₄","cmplxg34"), 
  ("E₆","weyle6"),("²E₆","weyl2e6"),("E₇","weyle7"),
  ("E₈","weyle8"), ("several groups","exceptio")]
println("reading Chevie data:")
foreach(cheviedata) do (e,f)
  print("for ",rpad(e,16))
  t=@elapsed include(string("tbl/",f,".jl"))
  println(rpad(string(round(t;digits=3)),5,'0')," seconds")
end
function contr(s) 
  include(replace(@__DIR__,"/src"=>"/contr/"*s*".jl"))
end
Testspath=replace(@__DIR__,"/src"=>"/test/Tests.jl")
export contr
roundtrip(x)=x==eval(Meta.parse(repr(x))) # for debugging purposes
end
