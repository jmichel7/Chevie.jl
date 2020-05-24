"""
Exemple:
```julia-repl
julia> a=Mod{19}(5)
5₁₉

julia> a^2
6₁₉

julia> inv(a)
4₁₉

julia> a*inv(a)
1₁₉

```
"""
module FFields
using ..Util: factor
export Mod, FF

const T=UInt8 # this allows moduli up to 255

struct Mod{p} <: Number
   val::T
   function Mod{p}(a::Integer) where {p}
   	new(T(mod(a,p)))
   end
end

Base.promote(x::Mod{p}, y::Integer) where p=(x,Mod{p}(y))
Base.promote(y::Integer, x::Mod{p}) where p=(Mod{p}(y),x)
Base.zero(::Type{Mod{p}}) where p = Mod{p}(T(0))
Base.:(==)(x::Mod,y::Mod) = x.val==y.val
Base.:(==)(x::Mod{p}, k::Integer) where p = mod(k,p) == x.val
Base.:(==)(k::Integer, x::Mod) = x==k
Base.:+(x::Mod{p}, y::Mod{p}) where p = Mod{p}(Int(x.val)+y.val)
Base.:*(x::Mod{p}, y::Mod{p}) where p = Mod{p}(Int(x.val)*y.val)
Base.:-(x::Mod{p}, y::Mod{p}) where p = Mod{p}(Int(x.val)-y.val)
Base.:-(x::Mod{p}) where {p} = Mod{p}(-Int(x.val))
Base.:/(x::Mod{p}, y::Mod{p}) where p = x * inv(y)
Base.inv(x::Mod{p}) where p = Mod{p}(invmod(x.val,p))
Base.:^(x::Mod{p},m::Integer) where p=(m>=0) ? Base.power_by_squaring(x,m) :
                                     Base.power_by_squaring(inv(x),-m)
Base.cmp(x::Mod{p}, y::Mod{p}) where p=cmp(x.val,y.val)
Base.isless(x::Mod{p}, y::Mod{p}) where p=cmp(x,y)==-1
Base.abs(x::Mod)=x      # needed for inv(Matrix) to work
Base.conj(x::Mod)=x     # needed for inv(Matrix) to work

function Base.show(io::IO, m::Mod{p}) where p 
   if get(io,:limit,false) 
     sub=Dict(zip("0123456789,()","₀₁₂₃₄₅₆₇₈₉‚₍₎"))
     print(io,m.val,map(x->sub[x],repr(p)))
   else print(io,"Mod{$p}($(m.val))")
   end
end

const conway_polynomials=Dict{Tuple{Int,Int},Vector{Int}}(
(2, 1) => [1],
(2, 2) => [1, 1],
(2, 3) => [1, 1, 0],
(2, 4) => [1, 1, 0, 0],
(2, 5) => [1, 0, 1, 0, 0],
(2, 6) => [1, 1, 0, 1, 1, 0],
(2, 7) => [1, 1, 0, 0, 0, 0, 0],
(2, 8) => [1, 0, 1, 1, 1, 0, 0, 0],
(2, 9) => [1, 0, 0, 0, 1, 0, 0, 0, 0],
(2, 10) => [1, 1, 1, 1, 0, 1, 1, 0, 0, 0],
(2, 11) => [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
(2, 12) => [1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0],
(2, 13) => [1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
(2, 14) => [1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
(2, 15) => [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
(3, 1) => [1],
(3, 2) => [2, 2],
(3, 3) => [1, 2, 0],
(3, 4) => [2, 0, 0, 2],
(3, 5) => [1, 2, 0, 0, 0],
(3, 6) => [2, 2, 1, 0, 2, 0],
(3, 7) => [1, 0, 2, 0, 0, 0, 0],
(3, 8) => [2, 2, 2, 0, 1, 2, 0, 0],
(3, 9) => [1, 1, 2, 2, 0, 0, 0, 0, 0],
(3, 10) => [2, 1, 0, 0, 2, 2, 2, 0, 0, 0],
(5, 1) => [3],
(5, 2) => [2, 4],
(5, 3) => [3, 3, 0],
(5, 4) => [2, 4, 4, 0],
(5, 5) => [3, 4, 0, 0, 0],
(5, 6) => [2, 0, 1, 4, 1, 0],
(7, 1) => [4],
(7, 2) => [3, 6],
(7, 3) => [4, 0, 6],
(7, 4) => [3, 4, 5, 0],
(7, 5) => [4, 1, 0, 0, 0],
(11, 1) => [9],
(11, 2) => [2, 7],
(11, 3) => [9, 2, 0],
(11, 4) => [2, 10, 8, 0],
(13, 1) => [11],
(13, 2) => [2, 12],
(13, 3) => [11, 2, 0],
(13, 4) => [2, 12, 3, 0],
(17, 1) => [14],
(17, 2) => [3, 16],
(17, 3) => [14, 1, 0],
(19, 1) => [17],
(19, 2) => [2, 18],
(19, 3) => [17, 4, 0],
(23, 1) => [18],
(23, 2) => [5, 21],
(23, 3) => [18, 2, 0],
(29, 1) => [27],
(29, 2) => [2, 24],
(29, 3) => [27, 2, 0],
(31, 1) => [28],
(31, 2) => [3, 29],
(31, 3) => [28, 1, 0],
(37, 1) => [35],
(37, 2) => [2, 33],
(37, 3) => [35, 6, 0],
(41, 1) => [35],
(41, 2) => [6, 38],
(43, 1) => [40],
(43, 2) => [3, 42],
(47, 1) => [42],
(47, 2) => [5, 45],
(53, 1) => [51],
(53, 2) => [2, 49],
(59, 1) => [57],
(59, 2) => [2, 58],
(61, 1) => [59],
(61, 2) => [2, 60],
(67, 1) => [65],
(67, 2) => [2, 63],
(71, 1) => [64],
(71, 2) => [7, 69],
(73, 1) => [68],
(73, 2) => [5, 70],
(79, 1) => [76],
(79, 2) => [3, 78],
(83, 1) => [81],
(83, 2) => [2, 82],
(89, 1) => [86],
(89, 2) => [3, 82],
(97, 1) => [92],
(97, 2) => [5, 96],
(101, 1) => [99],
(101, 2) => [2, 97],
(103, 1) => [98],
(103, 2) => [5, 102],
(107, 1) => [105],
(107, 2) => [2, 103],
(109, 1) => [103],
(109, 2) => [6, 108],
(113, 1) => [110],
(113, 2) => [3, 101],
(127, 1) => [124],
(127, 2) => [3, 126],
(131, 1) => [129],
(131, 2) => [2, 127],
(137, 1) => [134],
(137, 2) => [3, 131],
(139, 1) => [137],
(139, 2) => [2, 138],
(149, 1) => [147],
(149, 2) => [2, 145],
(151, 1) => [145],
(151, 2) => [6, 149],
(157, 1) => [152],
(157, 2) => [5, 152],
(163, 1) => [161],
(163, 2) => [2, 159],
(167, 1) => [162],
(167, 2) => [5, 166],
(173, 1) => [171],
(173, 2) => [2, 169],
(179, 1) => [177],
(179, 2) => [2, 172],
(181, 1) => [179],
(181, 2) => [2, 177],
(191, 1) => [172],
(191, 2) => [19, 190],
(193, 1) => [188],
(193, 2) => [5, 192],
(197, 1) => [195],
(197, 2) => [2, 192],
(199, 1) => [196],
(199, 2) => [3, 193],
(211, 1) => [209],
(211, 2) => [2, 207],
(223, 1) => [220],
(223, 2) => [3, 221],
(227, 1) => [225],
(227, 2) => [2, 220],
(229, 1) => [223],
(229, 2) => [6, 228],
(233, 1) => [230],
(233, 2) => [3, 232],
(239, 1) => [232],
(239, 2) => [7, 237],
(241, 1) => [234],
(241, 2) => [7, 238],
(251, 1) => [245],
(251, 2) => [6, 242],
)

function conway_polynomial(p,n)
  get!(conway_polynomials,(p,n)) do
    error("missing conway($p,$n)") 
  end
end

struct FF
  p::Int16
  n::Int16
  q::Int16
  conway::Vector
  zech::Vector{Int16}
end

struct FFE{q}<:Number
  i::Int16
  F::FF
end

printc(x...)=println(IOContext(stdout,:limit=>true),x...)

const FFDict=Dict{Int,FF}()

function FF(q)
  get!(FFDict,q) do
    l=collect(factor(q))
    if length(l)>1 error(q," should be a prime power") end
    p,n=l[1]
    pol=-Mod{p}.(conway_polynomial(p,n))
  # printc("conway=",pol)
    zz=map(i->fill(Mod{p}(0),n),1:q)
    zz[1][1]=Mod{p}(1)
    for i in 2:q-1
      for j in 2:n zz[i][j]=zz[i-1][j-1] end
      zz[i].+= zz[i-1][n].*pol
    end
    z=collect(enumerate(reverse.(zz)))
    sort!(z,by=x->x[2])
#   for i in 1:q printc(i,"->",(z[i][1]-1,z[i][2])) end
    zz=fill(99,q)
    for i in 1:q
      if z[i][2][end]==Mod{p}(p-1) zz[z[i][1]]=z[i-p+1][1]-1
      else zz[z[i][1]]=z[i+1][1]-1
      end
    end
#   for i in 0:q-1 printc(i,"=>",zz[i+1]) end
    FF(p,n,q,pol,Int16.(zz))
  end
end

Base.show(io::IO,F::FF)=print(io,"FF(",F.p,"^",F.n,")")

Base.:^(F::FF,i)=FFE{F.q}(mod(i,F.q-1),F)
Base.iszero(x::FFE{q}) where q=x.i==q-1
Base.:^(a::FFE, n::Integer)= n>=0 ? Base.power_by_squaring(a,n) :
                               Base.power_by_squaring(inv(a),-n)
Base.copy(a::FFE{q}) where q=FFE{q}(a.i,a.F)
Base.one(a::FFE{q}) where q=FFE{q}(0,a.F)
Base.one(::Type{FFE{q}}) where q=FFE{q}(0,FF(q))
Base.zero(::Type{FFE{q}}) where q=FFE{q}(q-1,FF(q))
Base.zero(x::FFE{q}) where q=FFE{q}(q-1,x.F)
Base.abs(a::FFE)=a
Base.conj(a::FFE)=a
Base.cmp(x::FFE, y::FFE)=cmp(x.i,y.i)
Base.isless(x::FFE, y::FFE)=cmp(x,y)==-1

function Base.show(io::IO, x::FFE{q})where q
  sup=Dict(zip("0123456789","⁰¹²³⁴⁵⁶⁷⁸⁹"))
  sub=Dict(zip("0123456789,()","₀₁₂₃₄₅₆₇₈₉‚₍₎"))
  repl=get(io,:limit,false)
  function printz()
    if repl print(io,"Z",map(x->sub[x],repr(q)))
    else print(io,"Z(",q,")")
    end
  end
  if iszero(x) 
     if repl print(io,"0",map(x->sub[x],repr(x.F.p)))
     else print(io,"0*");printz()
     end
  elseif isone(x) 
     if repl print(io,"1",map(x->sub[x],repr(x.F.p)))
     else printz();print(io,"^",x.i)
     end
  else  printz()
    if x.i!=1
      if repl print(io,map(x->sup[x],repr(x.i))) else print(io,"^",x.i) end
    end
  end
end

function FFE{q}(x::Integer)where q
  if iszero(x) return zero(FFE{q})
  elseif isone(x) return one(FFE{q})
  else error("FFE($x) not implemented")
  end
end

function Base.:+(x::FFE{q},y::FFE{q}) where q
  if iszero(x) return y
  elseif iszero(y) return x
  end
  if x.i>y.i x,y=y,x end
@inbounds res=x.F.zech[1+y.i-x.i]
  if res==q-1 return zero(x) end
  res+=x.i
  if res>=q-1 res-=q-1 end
  FFE{q}(res,x.F)
end

function Base.:-(x::FFE{q})where q
  if iseven(q) return x end
  if iszero(x) return x end
  res=x.i+div(q-1,2)
  if res>=q-1 res-=q-1 end
  FFE{q}(res,x.F)
end

Base.:-(x::FFE,y::FFE)=x+(-y)

function Base.:*(x::FFE{q},y::FFE{q}) where q
  if iszero(x) return x
  elseif iszero(y) return y
  end
  res=x.i+y.i
  if res>=q-1 res-=q-1 end
  FFE{q}(res,x.F)
end

function Base.inv(x::FFE{q}) where q
  if iszero(x) error("inv(0)") end
  if x.i==0 return x end
  FFE{q}(q-1-x.i,x.F)
end

Base.:/(x::FFE,y::FFE)=x*inv(y)
Base.://(x::FFE,y::FFE)=x*inv(y)

end
