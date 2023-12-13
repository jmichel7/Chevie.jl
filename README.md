# Chevie.jl

* [The documentation](https://jmichel7.github.io/Chevie.jl)

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
their   features.  I  have  implemented  some  other  infrastructure  which
currently resides in `Chevie` but may eventually become separate packages:

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
package installed, `Chevie` will automatically call GAP4 to implement these
functions.  The code  in this  package is  often 10  times faster  than the
equivalent GAP3 Chevie code (after the maddeningly long compilation time on
the first run --- Julia's TTFP).

The  `Chevie`  package  currently  contains  about  95%  of the GAP3/Chevie
functionality,  ported from  Gap3. If  you are  a user  of GAP3/Chevie, the
`gap`  function  can  help  you  to  find  the  equivalent functionality in
`Chevie.jl`  to a  Gap3 function:  it takes  a string  and gives  you Julia
translations of functions in Gap3 that match that string.

```julia-repl
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

This package requires julia 1.8 or later. 
This  is  my  effort  porting  GAP  code  to Julia, specifically the Chevie
package  of GAP3. I started this project at the end of 2018 and it is still
in  flux so the  package is not  yet registered. If  you see anything to be
improved  in this  package, please  contact me  or make  an issue or a pull
request in the GitHub repository.

This package requires julia 1.8 or later.
