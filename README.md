# Chevie.jl

* [The documentation](https://jmichel7.github.io/Chevie.jl)

This  is  my  effort  porting  GAP  code  to Julia, specifically the Chevie
package  of GAP3. I started this project at the end of 2018 and it is still
in  flux so the  package is not  yet registered. If  you see anything to be
improved  in this  package, please  contact me  or make  an issue or a pull
request in the GitHub repository.

### Installing

[For Julia novices]
To install this package, at the Julia command line:

  *  enter package mode with `]`
  *  do the command
```julia
(@v1.8) pkg> add "https://github.com/jmichel7/Chevie.jl"
```
  * exit package mode with backspace and then do
```julia
julia> using Chevie
```
and you are set up.

To update later to the latest version, do

```julia
(@v1.8) pkg> update Chevie
```
This package requires julia 1.8 or later.
