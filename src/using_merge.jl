function using_merge(mod::Symbol;reexport=false,debug=false)
  function myeval(e)
    if debug println(e) end
    eval(e)
  end
  mymodule=@__MODULE__
  mynames=names(mymodule;all=true)
  if !isdefined(mymodule,mod) myeval(:(using $mod: $mod)) end
  modnames=eval(:(names($mod)))
# println("names(",mod,")=",modnames)
  modnames=setdiff(modnames,[mod])
# println("now names(",mod,")=",modnames)
  noconflict=setdiff(modnames,mynames)
# println("non-conflicting:",noconflict)
  conflict=intersect(modnames,mynames)
  wrap(s)=Expr(:.,s)
  myeval(Expr(:using,Expr(:(:),Expr(:.,:.,mod),wrap.(noconflict)...)))
  if debug println("# names in $mod conflicting with $mymodule:",conflict) end
  for name in conflict
    methofname=eval(:(methods($name)))
    if isempty(methofname) modofname=nameof(mymodule)
    else modofname=nameof(first(methofname).module)
    end
#   println("conflicting name:$modofname.$name")
    s=split(sprint(show,eval(:(methods($mod.$name)))),"\n")
    map(s[2:end]) do l
#     println("   l=",l)
      l1=replace(l,r"^\[[0-9]*\] (.*) in .* at .*"=>s"\1")
      l1=replace(l1,r"#(s[0-9]*)"=>s"\1")
      e1=Meta.parse(l1)
#     println("\n   =>",e1)
      e2=e1.head==:where ? e1.args[1] : e1
      if e2.head!=:call || e2.args[1]!=name error("unexpected") end
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
#       elseif f.head==:parameters
#          e.args[i+1]=:($(Expr(:parameters, :(kw...))))
        end
      end
      e.args[1]=Expr(:.,mod,QuoteNode(e.args[1]))
      myeval(Expr(:(=),e1,e))
    end
  end
  if reexport
    pushfirst!(noconflict,mod)
    myeval(Expr(:export,noconflict...))
  end
end
