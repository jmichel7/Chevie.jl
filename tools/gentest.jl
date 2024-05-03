"generates a runtests.jl from julia-repl docstrings in file arguments"
function gentests(ff::Vector{String};module_="Chevie")
  open("runtests.jl","w")do outio
    write(outio,
"""
# auto-generated tests from julia-repl docstrings
using Test, $module_
function mytest(file::String,cmd::String,man::String)
  println(file," ",cmd)
  exec=repr(MIME("text/plain"),eval(Meta.parse(cmd)),context=:limit=>true)
  if endswith(cmd,";") return true end
  exec=replace(exec,r"\\s*\$"m=>""); exec=replace(exec,r"\\s*\$"s=>"")
  exec=replace(exec,r"^\\s*"=>"")
  if exec==man return true end
  i=findfirst(i->i<=lastindex(man) && exec[i]!=man[i],collect(eachindex(exec)))
  print("exec=\$(repr(exec[i:end]))\\nmanl=\$(repr(man[i:end]))\\n")
  false
end
@testset verbose = true "$module_" begin
"""
)
    for f in ff 
      println(outio,"@testset ",repr(f)," begin")
      println(f)
      open(f)do inio
        while !eof(inio)
          l=readline(inio)
          if match(r"^#@test",l)!=nothing
            println(outio,l[2:end])
          elseif match(r"```julia-repl$",l)!=nothing
            s=""
            while true
              u=readline(inio;keep=true)
              if match(r"```$",u)!=nothing break end
              s*=u
            end
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
              if !occursin(r"^ERROR:",man)
                println(outio,"@test mytest(",repr(f),",",repr(cmd),",",repr(man),")")
              end
            end
          end
        end
      end
      println(outio,"end")
    end
    println(outio,"end")
  end
end

gentests(ARGS)
