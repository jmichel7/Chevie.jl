"generates a runtests.jl from julia-repr docstrings in file arguments"
function gentests(ff::Vector{String};j16=false)
  open("runtests.jl","w")do io
    write(io,
"""
# auto-generated tests from julia-repl docstrings
using Test, Gapjm
#include("../tools/Gap4.jl")
function mytest(a::String,b::String)
  omit=a[end]==';'
  a=replace(a,"\\\\\\\\"=>"\\\\")
  a=repr(MIME("text/plain"),eval(Meta.parse(a)),context=:limit=>true)
  if omit a="nothing" end
  a=replace(a,r" *(\\n|\$)"s=>s"\\1")
  a=replace(a,r"\\n\$"s=>"")
  b=replace(b,r" *(\\n|\$)"s=>s"\\1")
  b=replace(b,r"\\n\$"s=>"")
  i=1
  while i<=lastindex(a) && i<=lastindex(b) && a[i]==b[i]
    i=nextind(a,i)
  end
  if a!=b print("exec=\$(repr(a[i:end]))\\nmanl=\$(repr(b[i:end]))\\n") end
  a==b
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
        c=""
        l=split(s,"\n")
        res=Pair{String,Vector{String}}[]
        for v in l[2:end-1]
          v=replace(v,r"# .*"=>"")
          v=replace(v,r" *$"=>"")
          if j16 
            v=replace(v,r"Array{([^{}]*),1}"=>s"Vector{\1}") 
            v=replace(v,r"Array{([^{}]*),2}"=>s"Matrix{\1}") 
            v=replace(v,r"Array{([^{}]*{[^{}]*}),1}"=>s"Vector{\1}") 
            v=replace(v,r"Array{([^{}]*{[^{}]*}),2}"=>s"Matrix{\1}") 
            v=replace(v,r"Array{([^{}]*{[^{}]*{[^{}]*}}),1}"=>s"Vector{\1}") 
            v=replace(v,r"Array{([^{}]*{[^{}]*{[^{}]*}}),2}"=>s"Matrix{\1}") 
            v=replace(v,r"Dict{([^{},]*),"=>s"Dict{\1, ") 
          end
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
          println("<$l>")
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
          write(io,"@test mytest(",repr(c),",",repr(b),")\n")
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
