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
using ..Util: factor, divisors
export Mod, FF, FFE, Z

const T=UInt8 # this allows moduli up to 255

struct Mod{p} <: Number
   val::T
   function Mod{p}(a::Integer) where {p}
   	new(T(mod(a,p)))
   end
end

Base.promote_type(::Type{Mod{p}},::Type{<:Integer}) where p=Mod{p}
Base.zero(::Type{Mod{p}}) where p = Mod{p}(0)
Base.:(==)(x::Mod,y::Mod) = x.val==y.val
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
(::Type{T})(x::Mod) where T<:Integer=T(x.val)

function Base.show(io::IO, ::MIME"text/plain", r::Mod{p})where p
  if !haskey(io,:typeinfo) print(io,typeof(r),": ") end
  show(io,r)
end

function Base.show(io::IO, m::Mod{p}) where p 
  if get(io,:limit,false) print(io,m.val)
  else print(io,typeof(m),"(",m.val,")")
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
(17, 1) => [14], (17, 2) => [3, 16], (17, 3) => [14, 1, 0],
(19, 1) => [17], (19, 2) => [2, 18], (19, 3) => [17, 4, 0],
(23, 1) => [18], (23, 2) => [5, 21], (23, 3) => [18, 2, 0],
(29, 1) => [27], (29, 2) => [2, 24], (29, 3) => [27, 2, 0],
(31, 1) => [28], (31, 2) => [3, 29], (31, 3) => [28, 1, 0],
(37, 1) => [35], (37, 2) => [2, 33], (37, 3) => [35, 6, 0],
(41, 1) => [35], (41, 2) => [6, 38],
(43, 1) => [40], (43, 2) => [3, 42],
(47, 1) => [42], (47, 2) => [5, 45],
(53, 1) => [51], (53, 2) => [2, 49],
(59, 1) => [57], (59, 2) => [2, 58],
(61, 1) => [59], (61, 2) => [2, 60],
(67, 1) => [65], (67, 2) => [2, 63],
(71, 1) => [64], (71, 2) => [7, 69],
(73, 1) => [68], (73, 2) => [5, 70],
(79, 1) => [76], (79, 2) => [3, 78],
(83, 1) => [81], (83, 2) => [2, 82],
(89, 1) => [86], (89, 2) => [3, 82],
(97, 1) => [92], (97, 2) => [5, 96],
(101, 1) => [99], (101, 2) => [2, 97],
(103, 1) => [98], (103, 2) => [5, 102],
(107, 1) => [105], (107, 2) => [2, 103],
(109, 1) => [103], (109, 2) => [6, 108],
(113, 1) => [110], (113, 2) => [3, 101],
(127, 1) => [124], (127, 2) => [3, 126],
(131, 1) => [129], (131, 2) => [2, 127],
(137, 1) => [134], (137, 2) => [3, 131],
(139, 1) => [137], (139, 2) => [2, 138],
(149, 1) => [147], (149, 2) => [2, 145],
(151, 1) => [145], (151, 2) => [6, 149],
(157, 1) => [152], (157, 2) => [5, 152],
(163, 1) => [161], (163, 2) => [2, 159],
(167, 1) => [162], (167, 2) => [5, 166],
(173, 1) => [171], (173, 2) => [2, 169],
(179, 1) => [177], (179, 2) => [2, 172],
(181, 1) => [179], (181, 2) => [2, 177],
(191, 1) => [172], (191, 2) => [19, 190],
(193, 1) => [188], (193, 2) => [5, 192],
(197, 1) => [195], (197, 2) => [2, 192],
(199, 1) => [196], (199, 2) => [3, 193],
(211, 1) => [209], (211, 2) => [2, 207],
(223, 1) => [220], (223, 2) => [3, 221],
(227, 1) => [225], (227, 2) => [2, 220],
(229, 1) => [223], (229, 2) => [6, 228],
(233, 1) => [230], (233, 2) => [3, 232],
(239, 1) => [232], (239, 2) => [7, 237],
(241, 1) => [234], (241, 2) => [7, 238],
(251, 1) => [245], (251, 2) => [6, 242],
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
  dict::Vector{Int16}
end

struct FFE{p}<:Number
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
    if n==1
      dic=Vector{Int16}(undef,q)
      val=pol[1].val
#     println(val)
      dic[1]=p-1
      dic[2]=0
      c=1
      for i in 1:p-2
        c=mod(c*val,p)
        dic[c+1]=i
      end
    else
      l=reverse(divisors(n))
      dic=fill(Int16(n),q)
      for i in 2:length(l)
        dic[1:div(q-1,p^l[i]-1):q].=l[i]
      end
    end
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
    zz=Vector{Int16}(undef,q)
    for i in 1:q
      if z[i][2][end]==Mod{p}(p-1) zz[z[i][1]]=z[i-p+1][1]-1
      else zz[z[i][1]]=z[i+1][1]-1
      end
    end
#   for i in 0:q-1 printc(i,"=>",zz[i+1]) end
    FF(p,n,q,pol,zz,dic)
  end
end

Base.show(io::IO,F::FF)=print(io,"FF(",F.p,"^",F.n,")")

@inline q(x::FFE{p}) where p=x.F.q
@inline n(x::FFE{p}) where p=x.F.n
@inline zech(x,i)=x.F.zech[1+i]
@inline dict(x,i)=x.F.dict[1+i]
@inline clone(x::FFE{p},i) where p=FFE{p}(i,x.F)

Base.iszero(x::FFE{p}) where p=q(x)==p && x.i==p-1

Base.:^(a::FFE{p},n::Integer) where p=iszero(a) ? a : 
                 lower(clone(a,mod(a.i*n,q(a)-1)))

#Base.:(==)(::Type{FFE{p}},::Type{FFE{q}}) where {p,q}=p==q

function Z(q)
  if q==2 return FFE{2}(0,FF(2)) end
  F=FF(q)
  FFE{F.p}(1,F)
end

# next 2 lines are absurd should not be needed
function Base.promote_type(::Type{FFE{p}},::Type{FFE{q}})where {p,q}
  if p!=q error("different charcateristic") end
  FFE{p}
end

#Base.promote(a::FFE{p},b::FFE{q}) where {p,q}=(FFE{p}(a),FFE{p}(b))
# next line is absurd should not be needed
(::Type{FFE{p}})(x::FFE{q}) where {p,q}=clone(x,x.i)

Base.copy(a::FFE{p}) where p=clone(a,a.i)

Base.one(::Type{FFE{p}}) where p=FFE{p}(0,FF(p))

Base.one(x::FFE{p}) where p=q(x)==p ? clone(x,0) : one(FFE{p})

Base.zero(::Type{FFE{p}}) where p=FFE{p}(p-1,FF(p))

Base.zero(x::FFE{p}) where p=q(x)==p ? clone(x,p-1) : zero(FFE{p})

Base.abs(a::FFE)=a
Base.conj(a::FFE)=a

function Base.cmp(x::FFE{p}, y::FFE{p}) where {p}
 l=cmp(q(x),q(y))
  if !iszero(l) return l end
  cmp(x.i,y.i)
end

Base.isless(x::FFE, y::FFE)=cmp(x,y)==-1
Base.:(==)(x::FFE, y::FFE)=cmp(x,y)==0

function Base.show(io::IO, x::FFE{p})where p
  sup=Dict(zip("0123456789","⁰¹²³⁴⁵⁶⁷⁸⁹"))
  sub=Dict(zip("0123456789,()","₀₁₂₃₄₅₆₇₈₉‚₍₎"))
  repl=get(io,:limit,false)
  function printz()
    if repl print(io,"Z",map(x->sub[x],repr(q(x))))
    else print(io,"Z(",q(x),")")
    end
  end
  if iszero(x) 
     if repl print(io,"0",map(x->sub[x],repr(p)))
     else print(io,"0*");printz()
     end
  elseif isone(x) 
     if repl print(io,"1",map(x->sub[x],repr(p)))
     else printz();print(io,"^",x.i)
     end
  else  printz()
    if x.i!=1
      if repl print(io,map(x->sup[x],repr(x.i))) else print(io,"^",x.i) end
    end
  end
end

function (::Type{FFE{p}})(i::Integer)where p
  F=FF(p)
  FFE{p}(F.dict[1+mod(i,p)],F)
end

function (::Type{FFE{p}})(i::Mod{p})where p
  F=FF(p)
  FFE{p}(F.dict[1+i.val],F)
end

function (::Type{Mod{r}})(x::FFE{p})where {r,p}
  if r!=p || q(x)!=p throw(InexactError(:convert,Mod{r},x)) end
  iszero(x) ? Mod{p}(0) : Mod{p}(powermod(x.F.conway[1].val,x.i,p))
end

function (::Type{FFE{p}})(i::Rational{<:Integer})where p
  F=FF(p)
  i=mod(invmod(denominator(i),p)*numerator(i),p)
  FFE{p}(F.dict[1+i],F)
end

function promote_field(r::Integer,x::FFE{p})where {p}
  if r==q(x) return x end
  if r<q(x) error("cannot promote to smaller field") end
  F=FF(r)
  if F.p!=p error("different characteristic!") end
  if mod(F.n,n(x))!=0 error("not an extension field") end
  FFE{p}(x.i*div(r-1,q(x)-1),F)
end
  
function promote_field(x::FFE{p},y::FFE{p})where {p}
  if q(x)==q(y) return (x,y) end
  nq=p^lcm(n(x),n(y))
  (promote_field(nq,x),promote_field(nq,y))
end

#Base.promote_rule(::Type{FFE{p}},::Type{T}) where {p,T<:Integer}=FFE{q}
#Base.promote_rule(::Type{FFE{p}},::Type{T}) where {p,T<:Rational{<:Integer}}=FFE{q}

function Base.:+(x::FFE{p},y::FFE{p}) where {p}
  if iszero(x) return y
  elseif iszero(y) return x
  end
  x,y=promote_field(x,y)
  if x.i>y.i x,y=y,x end
  @inbounds res=zech(x,y.i-x.i)
if res==q(x)-1 return zero(x) end
  res+=x.i
  if res>=q(x)-1 res-=Int16(q(x)-1) end
  lower(clone(x,res))
end

Base.:+(x::FFE{p},y::Union{Rational,Integer}) where p=x+FFE{p}(y)
Base.:+(y::Union{Rational,Integer},x::FFE{p}) where p=x+FFE{p}(y)

function Base.:-(x::FFE{p})where p
  if iseven(p) || iszero(x) return x end
  res=x.i+div(q(x)-1,2)
  if res>=q(x)-1 res-=Int16(q(x)-1) end
  clone(x,res)
end

Base.:-(x::FFE,y::FFE)=x+(-y)

function Base.:*(x::FFE{p},y::FFE{p}) where {p}
  if iszero(x) return x
  elseif iszero(y) return y
  end
  (x,y)=promote_field(x,y)
  res=x.i+y.i
  if res>=q(x)-1 res-=Int16(q(x)-1) end
  lower(clone(x,res))
end

Base.:*(x::FFE{p},y::Union{Rational,Integer}) where p=x*FFE{p}(y)
Base.:*(y::Union{Rational,Integer},x::FFE{p}) where p=x*FFE{p}(y)

function Base.inv(x::FFE{p}) where p
  if iszero(x) error("inv(0)") end
  if x.i==0 return x end
  clone(x,q(x)-1-x.i)
end

Base.:/(x::FFE,y::FFE)=x*inv(y)
Base.://(x::FFE,y::FFE)=x*inv(y)

function lower(x::FFE{p}) where p
  if q(x)==p return x end
  l=dict(x,x.i)
  if l==n(x) return x end
  FFE{p}(div(x.i,div(q(x)-1,p^l-1)),FF(p^l))
end

end
