"""
The module GAPENV creates a GAP3-like environment by extending some base 
functions like +,*,... and giving definitions of quite a few GAP functions,
then loads the CHEVIE library in that environment.
"""
module GAPENV
using ..Gapjm
function pad(s, i::Int)
  if i>0 return lpad(string(s),i)
  else return rpad(string(s),-i)
  end
end

# correct translations of GAP3 functions
pad(s::String)=s

include("../tools/gap3support.jl")
include("cheviesupport.jl")

const src=[
"cmp4_22", "cmp4_22_t","cmplxg24", "cmplxg25", "cmplxg26", "cmplxg27", "cmplxg29",
"cmplxg31", "cmplxg32", "cmplxg33", "cmplxg34", "cmplximp","cmplximp_t", 
"cmpxtimp", "coxh3", "coxh4", "coxi", "coxi_t", "weyla", "weyla_t", "weylbc", 
"weyld", "weyld_t", "weyl2a", "weyl2d",
"cox2i", "weyl2e6", "weyl2f4", "weyl3d4", "weyle6", "weyle7", "weyle8",
"weylf4", "weylg2"]


print("reading tbl/")
for f in src
  print("$f.jl, ")
  include("tbl/$f.jl")
end
println("\nreading tbl/exceptio.jl")
include("tbl/exceptio.jl")
println("\nreading tbl/exceptio_t.jl")
include("tbl/exceptio_t.jl")
end

#------- Loaded outside of GAP env since they don't need it ------------
const src_t=[ "cmplxg31_t", "cmplxg34_t", 
"cmpxtimp_t", "weylbc_t", "weyl2a_t", "weyl2d_t", "weyle8_t"]
print("\nreading tbl/")
for f in src_t
  print("$f.jl, ")
  include("tbl/$f.jl")
end
