The wish for this started in
[this thread](https://discourse.julialang.org/t/function-name-conflict-adl-function-merging/10335/7).
at  the time I was very new to Julia  and did not think I could do anything
myself  about the  problem. Now  two years  later, knowing  better Julia, I
realized I can do something myself. Here it is.

I introduce the problem with an example.

In  my big `Gapjm` package (a port of some GAP libraries to Julia) I have a
function  `invariants` which computes the invariants of a finite reflection
group.  However when  I use  `BenchmarkTools` to  debug for  performance my
package, I have the following problem:

```
julia> G= # some group...

julia> invariants(G)
WARNING: both Gapjm and BenchmarkTools export "invariants"; uses of it in module Main must be qualified
ERROR: UndefVarError: invariants not defined
Stacktrace:
 [1] top-level scope at REPL[4]:1
 [2] eval(::Module, ::Any) at ./boot.jl:331
 [3] eval_user_input(::Any, ::REPL.REPLBackend) at /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.4/REPL/src/REPL.jl:86
 [4] run_backend(::REPL.REPLBackend) at /home/jmichel/.julia/packages/Revise/kqqw8/src/Revise.jl:1163
 [5] top-level scope at none:0
```

This is annoying! I do not want to have to qualify every call to `invariants`
just because I am timing my code! What can I do? Well, first I could just
import the methods I am using in `BenchmarkTools`:

```
julia> using BenchmarkTools: @btime
```
Actually,  every exported name  from `BenchmarkTools`, except `invariants`,
does not conflict with my code:

```
julia> names(BenchmarkTools)
30-element Array{Symbol,1}:
 Symbol("@ballocated")
 Symbol("@belapsed")
 Symbol("@benchmark")
 Symbol("@benchmarkable")
 Symbol("@btime")
 Symbol("@tagged")
 :BenchmarkGroup
 :BenchmarkTools
 :addgroup!
 :allocs
 :gctime
 :improvements
 :invariants
 :isimprovement
 :isinvariant
 :isregression
 :judge
 :leaves
 :loadparams!
 :mean
 :median
 :memory
 :params
 :ratio
 :regressions
 :rmskew
 :rmskew!
 :trim
 :tune!
 :warmup
```

so I can do:

```
julia> using BenchmarkTools: @ballocated, @belapsed, @benchmark, @benchmarkable, @btime, @tagged, BenchmarkGroup, BenchmarkTools, addgroup!, allocs, gctime, improvements, isimprovement, isinvariant, isregression, judge, leaves, loadparams!, mean, median, memory, params, ratio, regressions, rmskew, rmskew!, trim, tune!,
 warmup
```

Still no conflict. Can I go further and do something even for `invariants`?
Well, I have one method for `invariants` in my package:

```
invariants(a::Group, args...)
```
while `BenchmarkTools` has four:

```
invariants(group::BenchmarkGroup)
invariants(x)
invariants(f, group::BenchmarkGroup)
invariants(f, x)
```
Even though some of these last methods apply to `Any`, they do not conflict
with my method, so I can use them also by just defining:

```
invariants(group::BenchmarkGroup) = BenchmarkTools.invariants(group)
invariants(x) = BenchmarkTools.invariants(x)
invariants(f, group::BenchmarkGroup) = BenchmarkTools.invariants(f, group)
invariants(f, x) = BenchmarkTools.invariants(f, x)
```
I  call  the  end  result  of  the  above  process  `merging`  the  package
`BenchmarkTools`  with my  current package.  What I  announce is a function
`using_merge` which does all the above automatically. If you call

```
julia> using_merge(:BenchmarkTools)
```

The function determines conflicting method and macro names in the package
and merges them as above, and imports the non-conflicting ones.

You will find the source code for this function in the `src` directory at

https://github.com/jmichel7/Gapjm.jl

This is just a function (not a module or a package) because:

- a simple function can do the job, and the function uses `eval` which needs
  to eval in the current module, so I find it easier to just `include` the
  function rather than to find how to `eval` in the calling module.

- I wait  for feedback  (which I  hope will  come) before  thinking how to
package the thing in a Julia package (or not), and whether I should provide
a macro rather than a function. Also my implementation is perhaps not the
best, I kind of parse the output of `methods`.

The function has two optional keyword arguments:

```
julia> using_merge(:BenchmarkTools;reexport=true)
```
will reexport all non-conflicting names.

```
julia> using_merge(:BenchmarkTools;debug=1)
```

will print all executed statements, and `debug=2` will describe even more
verbosely the taken actions.

Since  I wrote this function,  I found that I  got the hoped for modularity
benefits in my code. For example, I have in my `Gapjm.jl` package modules

  - `Perms`      Permutations
  - `Cycs`       Cyclotomic numbers (sums of complex roots of unity)
  - `Pols`       Univariate Laurent polynomials
  - `Mvps.jl`    Multivariate Puisuex polynomials
  - `Posets.jl`  Posets
  - `FFields.jl` Finite fields

that I designed as independent, stand-alone packages which each can be used
without importing anything else from my package. To use them together I can
now  just `using_merge` each  of them instead  of writing (unpleasant) glue
code.

I  do not advocate always  replacing the semantics of  `using` with that of
`using_merge`,  but  I  feel  that  `using_merge`  is a nice tool for using
packages  together without having to write glue code (and without having to
modify  any of the used packages). The meaning of "pirating a type" becomes
a   little  bit  wider  in  this  context,  as  you  saw  with  the  method
`invariants(y)`  in `BenchmarkTools`: it is, I would say, polite, if any of
your  methods which has  a possibly conflicting  name uses at  least one of
your own types.

There  are two  technical problems  left (that  I know  of --  there may be
others I am not aware):

-  I do not know how to  implement the above scheme for conflicting macros.
For example, you can see that the signature of `@btime` is:

```
julia> methods(eval(Symbol("@btime")))
# 1 method for macro "@btime":
[1] @btime(__source__::LineNumberNode, __module__::Module, args...)
...
```

I do not know how to write a `macro btime` which forwards this declaration
to the current module.

- I do not know how to forward docstrings. With the declaration

```
invariants(x) = BenchmarkTools.invariants(x)
```
a  docstring  of  `BenchmarkTools.invariants`  is  not  attached to the new
definition,  so `?invariants` will  not show the  docstrings. You can still
see  them by `?BenchmarkTools.invariants` but they were not merged together
with  the definitions, which partially defeats my goal. If someone tells me
how I could fix the above problems, I would be grateful.
