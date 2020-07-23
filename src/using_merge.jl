"""
`using_merge(mod::Symbol;reexport=false,debug=0)`

`using`  module `:mod`,  merging all  method definitions  exported by `mod`
with methods of the same name in the current module.

If  `reexport=true`  all  non-conflicting  names  exported  from  `mod` are
(re-)exported from the current module.

if `debug=1` every executed command is printed before execution.
If `debug=2` more verbose debugging information is printed.

It  is an error if  a name exported by  `mod` conflicts with anything but a
method  (if there is good reason  for this not to be  an error I could just
not import such names from `mod`).

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

julia> using_merge(:Bar;debug=1)
Main.foo(x::Float64) = Bar.foo(x)
```
"""
function using_merge(mod::Symbol;reexport=false,debug=0)
  function myeval(e)
    if debug>0 println(e) end
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
    modofname=nameof(first(methofname).module)
    if debug>1 println("# conflicting name:$modofname.$name") end
    s=split(sprint(show,methods(eval(:($mod.$name)))),"\n")
    map(s[2:end]) do l
      if debug>2 println("   l=",l) end
      l1=replace(l,r"^\[[0-9]*\] (.*) in .* at .*"=>s"\1")
      l1=replace(l1,r"#(s[0-9]*)"=>s"\1")
      e1=Meta.parse(l1)
      if debug>2 println("\n   =>",e1) end
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
