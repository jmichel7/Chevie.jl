"""
The  module GAPENV creates a GAP3-like environment by extending locally the
base  functions `*, +, -, ^, isless, copy, //, inv, length, union` to their
GAP3  semantics and defining  quite a few  other GAP3 functions, then loads
the Chevie database in that environment.

The  idea how to implement base function pirating restricted to a module is
due to Max Horn.
"""
module GAPENV
using ..Gapjm
include("../tools/gap3support.jl")
include("cheviesupport.jl")

const tbl=[
("piled data for G₄-G₂₂","cmp4_22"),
("lated data for G₄-G₂₂","cmp4_22_t"),
("piled data for G₂₄","cmplxg24"), 
("piled data for G₂₅","cmplxg25"), 
("piled data for G₂₆","cmplxg26"), 
("piled data for G₂₇","cmplxg27"), 
("piled data for G₂₉","cmplxg29"), 
("piled data for G₃₁","cmplxg31"), 
("piled data for G₃₂","cmplxg32"), 
("piled data for G₃₃","cmplxg33"), 
("piled data for G₃₄","cmplxg34"), 
("piled data for G(de,e,n)","cmplximp"),
("lated data for G(de,e,n)","cmplximp_t"), 
("piled data for ᵗG(e,e,n)","cmpxtimp"), 
("piled data for H₃","coxh3"), 
("piled data for H₄","coxh4"), 
("piled data for I₂(e)","coxi"), 
("lated data for I₂(e)","coxi_t"), 
("piled data for ²I₂(e)","cox2i"), 
("piled data for Aₙ","weyla"), 
("lated data for Aₙ","weyla_t"), 
("piled data for ²Aₙ","weyl2a"),
("piled data for Bₙ and Cₙ","weylbc"), 
("piled data for Dₙ","weyld"), 
("lated data for Dₙ","weyld_t"), 
("piled data for ²Dₙ","weyl2d"),
("piled data for ³D₄","weyl3d4"), 
("piled data for E₆","weyle6"), 
("piled data for ²E₆","weyl2e6"), 
("piled data for E₇","weyle7"), 
("piled data for E₈","weyle8"),
("piled data for F₄","weylf4"), 
("piled data for ²F₄","weyl2f4"), 
("piled data for G₂","weylg2"),
("piled data shared by several groups","exceptio"),
("lated data shared by several groups","exceptio_t")]
foreach(tbl)do (e,f)
  print("reading trans",rpad(e,36))
  t=@elapsed include(string("tbl/",f,".jl"))
  println(rpad(string(round(t;digits=3)),5,'0')," seconds")
end
end
#-- these translations are loaded outside of GAPENV since they don't need it --
const tbl_t=[ 
("lated data for G₃₁","cmplxg31_t"), 
("lated data for G₃₄","cmplxg34_t"), 
("lated data for ᵗG(e,e,n)","cmpxtimp_t"), 
("lated data for Bₙ and Cₙ","weylbc_t"), 
("lated data for ²Aₙ","weyl2a_t"), 
("lated data for ²Dₙ","weyl2d_t"), 
("lated data for E₈","weyle8_t")]
foreach(tbl_t) do (e,f)
  print("reading trans",rpad(e,36))
  t=@elapsed include(string("tbl/",f,".jl"))
  println(rpad(string(round(t;digits=3)),5,'0')," seconds")
end
