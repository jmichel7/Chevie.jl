"""
The  module GAPENV creates a GAP3-like environment by extending locally the
base  functions `*, +, -, ^, isless, copy, //, inv, length` to their
GAP3  semantics and defining  quite a few  other GAP3 functions, then loads
the Chevie database in that environment.
"""
module GAPENV
using ..Chevie
include("../tools/gap3support.jl")
include("cheviesupport.jl")

const tbl=[ ("G₄-G₂₂","cmp4_22"), ("G₂₄","cmplxg24"), ("G₂₅","cmplxg25"), 
("G₂₆","cmplxg26"), ("G₂₇","cmplxg27"), ("G₂₉","cmplxg29"), ("G₃₁","cmplxg31"), 
("G₃₂","cmplxg32"), ("G₃₃","cmplxg33"), ("G(de,e,n)","cmplximp"),
("ᵗG(e,e,n)","cmpxtimp"), ("G₃₄","cmplxg34"), ("H₃","coxh3"), ("H₄","coxh4"), 
("I₂(e)","coxi"), ("²I₂(e)","cox2i"), ("Aₙ","weyla"), ("²Aₙ","weyl2a"),
("Bₙ and Cₙ","weylbc"), ("Dₙ","weyld"), ("²Dₙ","weyl2d"), ("³D₄","weyl3d4"), 
("E₆","weyle6"), ("²E₆","weyl2e6"), ("E₇","weyle7"), ("E₈","weyle8"),
("F₄","weylf4"), ("²F₄","weyl2f4"), ("G₂","weylg2"),
("several groups","exceptio")]
println("reading transpiled data:")
foreach(tbl)do (e,f)
  print("for ",rpad(e,16))
  t=@elapsed include(string("tbl/",f,".jl"))
  println(rpad(string(round(t;digits=3)),5,'0')," seconds")
end
end
#-- translations are loaded outside of GAPENV since they don't need it --
const tbl_t=[("G₄-G₂₂","cmp4_22_t"), ("G₃₁","cmplxg31_t"), ("G₃₄","cmplxg34_t"),
("G(de,e,n)","cmplximp_t"), ("ᵗG(e,e,n)","cmpxtimp_t"), ("Aₙ","weyla_t"), 
("Bₙ and Cₙ","weylbc_t"), ("Dₙ","weyld_t"), ("²Aₙ","weyl2a_t"), 
("²Dₙ","weyl2d_t"), ("E₈","weyle8_t"), ("I₂(e)","coxi_t"),
("several groups","exceptio_t")]
println("reading translated data:")
foreach(tbl_t) do (e,f)
  print("for ",rpad(e,16))
  t=@elapsed include(string("tbl/",f,".jl"))
  println(rpad(string(round(t;digits=3)),5,'0')," seconds")
end
