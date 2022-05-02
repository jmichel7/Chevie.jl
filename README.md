# Gapjm.jl

This  is  my  effort  porting  GAP  code  to Julia, specifically the Chevie
package  of  GAP3  plus  the  GAP  functionality needed for Chevie to work:
Cyclotomics,   Permutations,   Laurent   and   Puiseux  polynomials,  basic
permutation group operations, etc….

I  started this project at the  end of 2018 and it  is still in flux so the
package  is not yet registered. If you  see anything to be improved in this
package, please contact me or make a pull request.

##

* [The documentation](https://jmichel7.github.io/Gapjm.jl)

### Installing

[For Julia novices]
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
