# Gapjm.jl

This  is  my  effort  porting  GAP  code  to Julia, specifically the Chevie
package  of  GAP3  plus  the  GAP  functionality needed for Chevie to work:
Cyclotomics,   Permutations,   Laurent   and   Puiseux  polynomials,  basic
permutation group operations, etc….

I am rather new to Julia, git and github so I am not even sure this package
is  properly constituted; I did not try yet to register it. If you are more
competent  that me and see any anomaly  in this package, please write me or
make a pull request.

##

* [The documentation](https://jmichel7.github.io/Gapjm.jl)

### Installing

[For Julia novices]
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
