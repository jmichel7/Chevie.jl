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
  omit=src[end]==';'
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
      out=String[]
      blks=map(blks)do s
        s=replace(s,"\$Idef"=>"Int16") # hack: instead evaluate docstrings
        c=""
        l=split(s,"\n")
        res=Pair{String,Vector{String}}[]
        for v in l[2:end-1]
          v=replace(v,r"# .*"=>"")
          v=replace(v,r" *$"=>"")
          if occursin(r"^julia>",v)
            if c!="" push!(res,c=>out) end
            out=String[]
            c=replace(v,r"^julia>\s*"=>"")
          else
            push!(out,v)
          end
        end
        push!(res,c=>out)
      end
      d=vcat(blks...)
      d=map(d)do (c,l)
        while true
#         println("<$l>")
          if isempty(l) l=["nothing"] end
          if l[end]=="" 
            if length(l)>1 resize!(l,length(l)-1)
            else l[end]="nothing"
            end
          else break
          end
        end
        c=>join(l,"\n")
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
