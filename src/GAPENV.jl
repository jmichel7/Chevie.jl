
include("tbl/util.jl")

#-- translations are loaded outside of GAPENV since they don't need it --
const tbl_t=[("G(de,e,n)","cmplximp_t"),
  ("Aₙ","weyla_t"), ("Bₙ and Cₙ","weylbc_t"),("Dₙ","weyld_t"),
  ("³D₄","weyl3d4_t"), ("²Aₙ","weyl2a_t"), ("²Dₙ","weyl2d_t"),
  ("ᵗG(e,e,n)","cmpxtimp_t"), ("G₂","weylg2_t"), 
  ("I₂(e)","coxi_t"), ("²I₂(e)","cox2i_t"),
  ("G₄-G₂₂","cmp4_22_t"),("H₃","coxh3_t"),
  ("G₂₄","cmplxg24_t"),("G₂₅","cmplxg25_t"),("G₂₆","cmplxg26_t"),
  ("G₂₇","cmplxg27_t"), ("F₄","weylf4_t"), ("²F₄","weyl2f4_t"),
  ("G₂₉","cmplxg29_t"), ("H₄","coxh4_t"), ("G₃₁","cmplxg31_t"),
  ("G₃₂","cmplxg32_t"),("G₃₃","cmplxg33_t"),("G₃₄","cmplxg34_t"), 
  ("E₆","weyle6_t"),("²E₆","weyl2e6_t"),("E₇","weyle7_t"),
  ("E₈","weyle8_t"), ("several groups","exceptio_t")]
println("reading Chevie data:")
foreach(tbl_t) do (e,f)
  print("for ",rpad(e,16))
  t=@elapsed include(string("tbl/",f,".jl"))
  println(rpad(string(round(t;digits=3)),5,'0')," seconds")
end
