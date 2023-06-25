"generates a runtests.jl from julia-repr docstrings in file arguments"
function gentests(ff::Vector{String})
  open("runtests.jl","w")do io
    write(io,
"""
# auto-generated tests from julia-repl docstrings
using Test, Gapjm
#include("../tools/Gap4.jl")
function mytest(file::String,src::String,man::String)
  println(file," ",src)
  omit=match(r";\\s*(#.*)?\$",src)!==nothing
  src=replace(src,"\\\\\\\\"=>"\\\\")
  exec=repr(MIME("text/plain"),eval(Meta.parse(src)),context=:limit=>true)
  if omit exec="nothing" end
  exec=replace(exec,r" *(\\n|\$)"s=>s"\\1")
  exec=replace(exec,r"\\n\$"s=>"")
  man=replace(man,r" *(\\n|\$)"s=>s"\\1")
  man=replace(man,r"\\n\$"s=>"")
  i=1
  while i<=lastindex(exec) && i<=lastindex(man) && exec[i]==man[i]
    i=nextind(exec,i)
  end
  if exec!=man 
    print("exec=\$(repr(exec[i:end]))\\nmanl=\$(repr(man[i:end]))\\n")
  end
  exec==man
end
@testset verbose = true "Gapjm" begin
"""
)
    for f in ff 
      blks=String[]
      s=read(f,String)
      s=replace(s,r"```julia-repl(.*?)```"s=>t->push!(blks,t))
      blks=map(s->replace(s,r"```julia-repl\s(.*?)\s```"s=>s"\1"),blks)
      out=String[]
      blks=map(blks)do s
        s=replace(s,"\$Idef"=>"Int16") # hack: instead evaluate docstrings
        c=""
        l=split(s,r"^julia> "m)[2:end]
        map(l)do t
          cmd,i=Meta.parse(t,1)
          t[1:i-1]=>t[i:end]
        end
      end
      d=vcat(blks...)
      d=map(d)do (c,l)
        c=chomp(c)
        l=chomp(l)
        if match(r"^\s*$",l)!==nothing l="nothing" end
        c=>l
      end
      if !isempty(d)
      write(io,"@testset ",repr(f)," begin\n")
      for (c,b) in d 
        if !occursin(r"^ERROR:",b)
          write(io,"@test mytest(",repr(f),",",repr(c),",",repr(b),")\n")
        end
      end
      write(io,"end\n")
      end
    end
    write(io,"end\n")
  end
end

println(ARGS)
gentests(ARGS)
