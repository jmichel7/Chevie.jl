module HasType

export traces_words_mats, Unknown

using ..Gapjm
#-----------------------------------------------------------------------
struct Unknown end

Base.:+(a,b::Unknown)=b
Base.:+(b::Unknown,a)=b
Base.:+(b::Unknown,a::Unknown)=b
Base.:*(a,b::Unknown)=iszero(a) ? a : b
Base.:*(b::Unknown,a)=iszero(a) ? a : b
Base.:*(b::Unknown,a::Unknown)=b
Base.zero(a::Unknown)=0
Base.show(io::IO,a::Unknown)=print(io,get(io,:limit,false) ? "?" : "Unknown()")
Base.isless(a::Unknown,b::Number)=false
Base.isless(b::Number,a::Unknown)=true
Base.isless(b::Unknown,a::Unknown)=false
Base.:(//)(a::Unknown,b)=a
Base.isreal(a::Unknown)=false
Base.isinteger(a::Unknown)=false
Base.conj(a::Unknown)=a
Base.broadcastable(a::Unknown)=Ref(a)
CyclotomicNumbers.galois(a::Unknown,i)=a
#----------------------------------------------------------------------
# correct translations of GAP3 functions

include("../tools/gap3support.jl")
include("cheviesupport.jl")

function pad(s, i::Int)
  if i>0 return lpad(string(s),i)
  else return rpad(string(s),-i)
  end
end

pad(s::String)=s

"""
`traces_words_mats(mats,words)`

given  a list `mats`  of matrices and  a list `words`  of words returns the
list of traces of the corresponding products of the matrices

```julia-repl
julia> W=coxgroup(:F,4)
F₄

julia> r=classinfo(W)[:classtext];

julia> R=representation(W,17)
4-element Vector{Matrix{Int64}}:
 [-1 -1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
 [1 0 0 0; -1 -1 -1 0; 0 0 1 0; 0 0 0 1]
 [1 0 0 0; 0 1 0 0; 0 -2 -1 -1; 0 0 0 1]
 [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 -1 -1]

julia> traces_words_mats(R,r)==CharTable(W).irr[17,:]
true
```
"""
function traces_words_mats(mats,words)
  mats=improve_type(mats)
  dens=map(x->1,mats)
  if all(m->all(x->x isa Rational,m),mats)
    dens=map(m->lcm(denominator.(m)),mats)
    mats=map((m,d)->Int.(m.*d),mats,dens)
  end
  words=convert.(Vector{Int},words)
  tr(m)=sum(i->m[i,i],axes(m,1))
  trace(w)=tr(prods[w])//prod(dens[w])
  prods=Dict{Vector{Int},eltype(mats)}(Int[]=>mats[1]^0)
  for i in eachindex(mats) prods[[i]]=mats[i] end
  res=map(words)do w
    i=0
    while haskey(prods,w[1:i])
      if i==length(w) return trace(w) end
      i+=1
    end
    while i<=length(w)
      prods[w[1:i]]=prods[w[1:i-1]]*mats[w[i]]
      i+=1
    end
#   println(prod(dens[w])
    trace(w)
  end
end

#-----------------------------------------------------------------------

const src=[
"cmp4_22", "cmp4_22_t",
"cmplxg24",
"cmplxg25", 
"cmplxg26",
"cmplxg27",
"cmplxg29",
"cmplxg31", "cmplxg31_t",
"cmplxg32",
"cmplxg33", 
"cmplxg34", "cmplxg34_t",
"cmplximp", "cmplximp_t",
"cmpxtimp", "cmpxtimp_t",
"coxh3", 
"coxh4",
"coxi", "coxi_t",
"weyla", "weyla_t",
"weylbc", "weylbc_t",
"weyld", "weyld_t",
"weyl2a", "weyl2a_t",
"weyl2d", "weyl2d_t",
"cox2i", 
"weyl2e6", 
"weyl2f4", 
"weyl3d4",
"weyle6", 
"weyle7", 
"weyle8", "weyle8_t",
"weylf4", 
"weylg2"
]

print("reading tbl/")
for f in src
  print("$f.jl, ")
  include("tbl/$f.jl")
end
println("\nreading tbl/exceptio")
include("tbl/exceptio.jl")
println("reading tbl/exceptio_t")
include("tbl/exceptio_t.jl")

end
