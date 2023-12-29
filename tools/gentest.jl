"generates a runtests.jl from julia-repl docstrings in file arguments"
function gentests(ff::Vector{String};module_="Chevie")
  open("runtests.jl","w")do io
    write(io,
"""
# auto-generated tests from julia-repl docstrings
using Test, $module_
function mytest(file::String,cmd::String,man::String)
  println(file," ",cmd)
  exec=repr(MIME("text/plain"),eval(Meta.parse(cmd)),context=:limit=>true)
  if endswith(cmd,";") exec="nothing" 
  else exec=replace(exec,r"\\s*\$"m=>"")
       exec=replace(exec,r"\\s*\$"s=>"")
       exec=replace(exec,r"^\\s*"=>"")
  end
  if exec==man return true end
  i=1
  while i<=lastindex(exec) && i<=lastindex(man) && exec[i]==man[i]
    i=nextind(exec,i)
  end
  print("exec=\$(repr(exec[i:end]))\\nmanl=\$(repr(man[i:end]))\\n")
  return false
end
@testset verbose = true "$module_" begin
"""
)
    for f in ff 
      blks=first.(eachmatch(r"```julia-repl(.*?)```"s,read(f,String)))
      stmts=vcat(map(blks)do s
        map(split(s,r"^julia> "m)[2:end])do t
          if !('"' in t)  t=replace(t,"\\\\"=>"\\") end # hack
          e,i=Meta.parse(t,1)
          cmd=t[1:i-1]
          cmd=replace(cmd,r"\s*(#.*)?(\n\s*)?$"=>"") # replace final comment
          man=t[i:end]
#         if !('"' in man)  man=replace(man,"\\\\"=>"\\") end # hack
          man=replace(man,r"\s*$"m=>"")
          man=replace(man,r"\s*$"s=>"")
          if isempty(man) man="nothing" end
          cmd=>man
        end
      end...)
      if !isempty(stmts)
      write(io,"@testset ",repr(f)," begin\n")
      for (cmd,man) in stmts 
        if !occursin(r"^ERROR:",man)
          write(io,"@test mytest(",repr(f),",",repr(cmd),",",repr(man),")\n")
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
