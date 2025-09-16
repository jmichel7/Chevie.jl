# auto-generated tests from julia-repl docstrings
using Test, Chevie
function mytest(file::String,cmd::String,man::String)
  println(file," ",cmd)
  exec=repr(MIME("text/plain"),eval(Meta.parse(cmd)),context=:limit=>true)
  if endswith(cmd,";") return true end
  exec=replace(exec,r"\s*$"m=>""); exec=replace(exec,r"\s*$"s=>"")
  exec=replace(exec,r"^\s*"=>"")
  if exec==man return true end
  inds=collect(eachindex(exec))
  i=inds[findfirst(i->i<=lastindex(man) && exec[i]!=man[i],inds)]
  print("exec=$(repr(exec[i:end]))\nmanl=$(repr(man[i:end]))\n")
  false
end
@testset verbose = true "Chevie" begin
end
