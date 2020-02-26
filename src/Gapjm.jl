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

```julia-repl
julia> Gapjm.gap("words")
CoxeterWords(W[,l])      =>  word.(Ref(W),elements(W[,l])
GarsideWords             =>  elements
CharRepresentationWords  =>  traces_words_mats
```
Then you can call on-line help on the discovered functions.

This package is  often 10 times faster  than the equivalent GAP3 Chevie code
(after the maddeningly long compilation time on first execution).

I tried that parts of my package can be used independently of the rest. For
instance,  the modules `Combinat`,  `Groups`, `ModuleElts`, `Perms`, `Util`
are independent of the rest of the package and can be used stand-alone.

"""
module Gapjm
using Reexport

#--------------------------------------------------------------------------
# All  the glue  code below  should not  exist: extending e.g. Gapjm.degree
# instead  of just exporting degree from Pols  and Perms is a horrible hack
# forced  by the unpleasant Julia rules  not merging methods from different
# modules.
# 
# This  way both Pols  and Perms can  work together only  as part of Gapjm.
# Otherwise they could not work together.
# 
# If degree  was in  Base there  would be  no problem, both importing from
# Base.
export coefficients, degree, degrees, elements, kernel, order, restricted,
       root, roots, valuation, words, word
function coefficients end
degree(a::Number)=0
function degrees end
function elements end
function kernel end
function order end
function restricted end
function root end
function roots end
function valuation end
function words end
function word end

include("Util.jl")
@reexport using .Util
include("ModuleElts.jl")
@reexport using .ModuleElts
include("Groups.jl")
@reexport using .Groups
word(G::Group,x...)=Groups.word(G,x...)
elements(G::Group)=Groups.elements(G)
kernel(h::Hom)=Groups.kernel(h)
order(C::Coset)=Groups.order(C)
include("Perms.jl")
@reexport using .Perms
degree(a::Perm)=Perms.degree(a)
order(p::Perm)=Perms.order(p)
restricted(a::Perm,l::AbstractVector)=Perms.restricted(a,l)
Groups.orbit(a::Perm,x...)=Perms.orbit(a,x...)
Groups.orbits(a::Perm,x...)=Perms.orbits(a,x...)
include("PermGroups.jl")
@reexport using .PermGroups
include("Cycs.jl")
@reexport using .Cycs
coefficients(c::Cyc)=Cycs.coefficients(c)
root(r::Integer,x...)=Cycs.root(r,x...)
root(r::Rational,x...)=Cycs.root(r,x...)
root(r::Cyc,x...)=Cycs.root(r,x...)
include("Combinat.jl")
@reexport using .Combinat
include("Pols.jl")
@reexport using .Pols
degree(p::Pol)=Pols.degree(p)
valuation(p::Pol,x...)=Pols.valuation(p,x...)
root(p::Pol,x...)=Pols.root(p,x...)
include("Mvps.jl")
@reexport using .Mvps
degree(p::Mvps.Monomial,x...)=Mvps.degree(p,x...)
degree(p::Mvp,x...)=Mvps.degree(p,x...)
coefficients(p::Mvp,x...)=Mvps.coefficients(p,x...)
root(r::Monomial,x...)=Mvps.root(r,x...)
root(r::Mvp,x...)=Mvps.root(r,x...)
valuation(p::Mvp,x...)=Mvps.valuation(p,x...)
Mvp(p::Pol)=p(Mvp(Pols.varname[]))
include("PermRoot.jl")
@reexport using .PermRoot
roots(W::PRG)=PermRoot.roots(W)
roots(W::PRSG)=PermRoot.roots(W)
include("GLinearAlgebra.jl")
@reexport using .GLinearAlgebra
include("CoxGroups.jl")
@reexport using .CoxGroups
include("Weyl.jl")
@reexport using .Weyl
#roots(W::FiniteCoxeterGroup)=Weyl.roots(W)
#roots(C::Matrix)=Weyl.roots(C)
include("Cosets.jl")
@reexport using .Cosets
include("Chars.jl")
@reexport using .Chars
include("SPerms.jl")
@reexport using .SPerms
include("HeckeAlgebras.jl")
@reexport using .HeckeAlgebras
include("KL.jl")
@reexport using .KL
include("CycPols.jl")
@reexport using .CycPols
roots(c::CycPol)=CycPols.roots(c)
#degree(c::CycPol)=CycPols.degree(c)
include("Symbols.jl")
@reexport using .Symbols
include("Garside.jl")
@reexport using .Garside
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
include("Eigenspaces.jl")
@reexport using .Eigenspaces
include("../docs/src/cheviedict.jl")
export gap
end
