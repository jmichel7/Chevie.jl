module HasType

export charname, traces_words_mats

using ..Gapjm
#-----------------------------------------------------------------------
struct Unknown
end

Base.:+(a,b::Unknown)=b
Base.:+(b::Unknown,a)=b
Base.:*(a,b::Unknown)=b
Base.:*(b::Unknown,a)=b
Base.zero(a::Unknown)=0
#-----------------------------------------------------------------------
charname(t::TypeIrred,p;TeX=false,opt...)=getchev(t,:CharName,p,
                           TeX ? Dict(:TeX=>true) : Dict())

charname(W,x;TeX=false,opt...)=join(map((t,p)->charname(t,p;TeX=TeX,opt...),
                           refltype(W),x),",")

function PositionCartesian(l,ind)
  res=prod=1
  for i in length(l):-1:1
    res+=(ind[i]-1)*prod
    prod*=l[i]
  end
  res
end

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

PrintToSring(s,v...)=sprint(show,v...)

function Replace(s,p...)
# print("Replace s=$s p=$p")
  for (src,tgt) in (p[i]=>p[i+1] for i in 1:2:length(p))
    i=0
    while i+length(src)<=length(s)
     if src==s[i+(1:length(src))]
        if tgt isa String
          s=s[1:i]*tgt*s[i+length(src)+1:end]
        else
          s=vcat(s[1:i],tgt,s[i+length(src)+1:end])
        end
        i+=length(tgt)
      else
        i+=1
      end
    end
  end
# println("=>",s)
  s
end

ChevieIndeterminate(a::Vector{<:Number})=one(Pol)
ChevieIndeterminate(a::Vector{<:Pol})=Mvp(:x)

"""
`CycPol(v::AbstractVector)`

This  form is a fast  and efficient way of  specifying a `CycPol` with only
positive multiplicities: `v` should be a vector. The first element is taken
as  a  the  `.coeff`  of  the  `CycPol`,  the  second  as the `.valuation`.
Subsequent  elements are rationals `i//d`  representing `(q-E(d)^i)` or are
integers `d` representing `Φ_d(q)`.

```julia-repl
julia> CycPol([3,-5,6,3//7])
3q⁻⁵Φ₆(q-ζ₇³)
```
"""
function CycPols.CycPol(v::AbstractVector)
  coeff=v[1]
  valuation=convert(Int,v[2])
  vv=Pair{Root1,Int}[]
  v1=convert.(Rational{Int},v[3:end])
  for i in v1
    if denominator(i)==1
      k=convert(Int,i)
      for j in prime_residues(k) push!(vv,Root1(;r=j//k)=>1) end
    else
      push!(vv,Root1(;r=i)=>1)
    end
  end
  CycPol(coeff,valuation,ModuleElt(vv))
end

"""
`traces_words_mats(mats,words)`

given  a list `mats`  of matrices and  a list `words`  of words returns the
list of traces of the corresponding products of the matrices

```julia-repl
julia> W=coxgroup(:F,4)
F₄

julia> r=classinfo(W)[:classtext];

julia> R=representation(W,17)
4-element Array{Array{Int64,2},1}:
 [-1 -1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
 [1 0 0 0; -1 -1 -1 0; 0 0 1 0; 0 0 0 1]
 [1 0 0 0; 0 1 0 0; 0 -2 -1 -1; 0 0 0 1]
 [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 -1 -1]

julia> traces_words_mats(R,r)==CharTable(W).irr[17,:]
true
```
"""
function traces_words_mats(mats,words)
  mats=improve_type.(mats)
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
"cmplxg31",
"cmplxg32",
"cmplxg33", 
"cmplxg34",
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
"weyle8",
"weylf4", 
"weylg2",
"exceptio", "exceptio_t"
]

for f in src
  println("reading tbl/$f.jl")
  include("tbl/$f.jl")
end

end
