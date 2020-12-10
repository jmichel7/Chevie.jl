"""
`using_merge(mod::Symbol;reexport=false,debug=0)`

do  `using`  module  `:mod`,  in  addition  merging  all method definitions
exported by `mod` with methods of the same name in the current module.

If  `reexport=true`  all  non-conflicting  names  exported  from  `mod` are
(re-)exported from the current module.

if `debug=1` every executed command is printed before execution.
If `debug=2` conflicting methods are printed.
If `debug=3` boths happens.

It  is an error if a name exported  by `mod` conflicts with a name which is
not  a method (if someone points to me  a case where there is a good reason
for  this not to be an  error I could just decide  not to import such names
from `mod`).

Now an example of use:
```julia
julia> include("using_merge.jl")
using_merge

julia> foo(x::Int)=2
foo (generic function with 1 method)

julia> module Bar
       export foo
       foo(x::Float64)=3
       end
Main.Bar

julia> using_merge(:Bar;debug=3)
# Bar.foo conflicts with Main.foo --- adding methods
Main.foo(x::Float64) = Bar.foo(x)
```
"""
function using_merge(mod::Symbol;reexport=false,debug=0)
  function myeval(e)
    if !iszero(debug&1) println(e) end
    eval(e)
  end
  mymodule=@__MODULE__
  # first, bring in scope mod in case it is not
  if !isdefined(mymodule,mod) myeval(:(using $mod: $mod)) end
  if reexport myeval(:(export $mod)) end
  modnames=setdiff(eval(:(names($mod))),[mod])
  for name in modnames
    if !isdefined(mymodule,name) 
      myeval(Expr(:using,Expr(:(:),Expr(:.,:.,mod),Expr(:.,name))))
      if reexport myeval(:(export $name)) end
      continue
    end
    # now we know name is conflicting
    methofname=methods(eval(name))
    if isempty(methofname)
      # we do no handle conflicting names which are not methods or macros
      error("$mod.$name is not a method or macro in $mymodule")
    end
    modofname=nameof(parentmodule(eval(name)))
    if !iszero(debug&2) 
      print("# $mod.$name conflicts with $modofname.$name --- adding methods") 
    end
    s=split(sprint(show,methods(eval(:($mod.$name)))),"\n")
    for (j,l) in enumerate(s[2:end])
      if j==1 
        if !isempty(eval(:(@doc $mod.$name)).meta[:results])
          myeval(:(@doc (@doc $mod.$name) $modofname.$name))
          if !iszero(debug&2) println(" and docs") end
        else println() 
        end
      end
      if !iszero(debug&4) println("   l=",l) end
      l1=replace(l,r"^\[[0-9]*\] (.*) in .* at .*"=>s"\1")
      l1=replace(l1,r"#(s[0-9]*)"=>s"\1")
      e1=Meta.parse(l1)
      if !iszero(debug&4) println("\n   =>",e1) end
      e2=e1.head==:where ? e1.args[1] : e1
      if e2.head!=:call || e2.args[1]!=name 
        error("conflicting macros are not yet implemented ($name)") 
      end
      for (i,f) in enumerate(e2.args[2:end])
        if (f isa Expr) && f.head==:parameters
          e2.args[i+1]=:($(Expr(:parameters, :(kw...))))
        end
      end
      e=deepcopy(e1)
      e2.args[1]=Meta.parse("$modofname.$name")
      if e.head==:where e=e.args[1] end
      for (i,f) in enumerate(e.args[2:end])
        if !(f isa Expr) continue end
        if f.head==:(::) 
           if length(f.args)==2 e.args[i+1]=f.args[1]
           else e.args[i+1]=f.args[1].args[2]
           end
        end
      end
      e.args[1]=Expr(:.,mod,QuoteNode(e.args[1]))
      myeval(Expr(:(=),e1,e))
    end
  end
end
