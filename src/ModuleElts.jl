"""
A  `ModuleElt{K,V}`  represents  an  element  of  a free module where basis
elements  are of type  `K` and coefficients  of type `V`.  Usually you want
objects  of type `V` to be elements of  a ring, but it could also be useful
if they just belong to an abelian group. 

This  basic data structure is  used in the package  `Gapjm` as an efficient
representation   for  multivariate   polynomials  (and   their  monomials),
cyclotomics, CycPols, elements of Hecke algebras, etc…

A  `ModuleElt` is essentially a  list of pairs `b=>c`  where `b` is a basis
element  and `c` its coefficient. The  constructor takes as argument a list
of  pairs, or a variable number of pair arguments, or a generator of pairs.

We provide two implementations:

  - one by `Dict`s 
this is easy since the interface of the type is close to that of dicts; the
only  difference is weeding out keys which have a zero cofficient --- which
thus demands that the keys are hashable.

  - a faster  one (the default)  by sorting pairs by key
this  demands that keys are comparable.  This implementation is two to four
times faster than the `Dict` one and requires half the memory.

Here  is an example where basis elements are `Symbol`s and coefficients are
`Int`.

```julia-repl
julia> a=ModuleElt(:xy=>1,:yx=>-1)
:xy-:yx

julia> a-a
0

julia> a*99
99:xy-99:yx

julia> a+ModuleElt(:yx=>1)
:xy

julia> a[:xy]
1

julia> a[:xx]  # the coefficient of an absent basis element is zero.
0

julia> haskey(a,:xx)
false

julia> first(a)
:xy => 1

julia> collect(a)
2-element Array{Pair{Symbol,Int64},1}:
 :xy => 1
 :yx => -1

julia> keys(a)
2-element Array{Symbol,1}:
 :xy
 :yx

julia> values(a)
2-element Array{Int64,1}:
  1
 -1

julia> length(a)
2

julia> eltype(a)
Pair{Symbol,Int64}
```

both  implementations provide  an option  `check` in  the constructor which
normalizes  an element,  removing zero  coefficients and  merging duplicate
basis elements (and sorting the basis in the default implementation).

```julia-repl
julia> a=ModuleElt(:yy=>1, :yx=>2, :xy=>3, :yy=>-1)
:yy+2:yx+3:xy-:yy

julia> a=ModuleElt([:yy=>1, :yx=>2, :xy=>3, :yy=>-1];check=true)
3:xy+2:yx

julia> a
3:xy+2:yx
```

setting  the `IOContext` property `showbasis` to a custom printing function
changes how the basis elements are printed.
```julia-rep1
julia> show(IOContext(stdout,:showbasis=>(io,s)->string("<",s,">")),a)
3<xy>+2<yx>
```

Adding  or subtracting `ModuleElt`s does promotion  on the type of the keys
and the coefficients if needed:

```julia-repl
julia> a+ModuleElt([:z=>1.0])
3.0:xy+2.0:yx+1.0:z
```

"""
module ModuleElts

export ModuleElt # data structure

#------------- implementation with Dicts ----------------------
if false # change to true to use this implementation
struct ModuleElt{K,V}
  d::Dict{K,V} # the keys K should be hashable
  function ModuleElt(d::Dict{K,V};check::Bool=false)where {K,V}
    if check
      for (k,v) in d if iszero(v) delete!(d,k) end end
    end
    new{K,V}(d)
  end
end

ModuleElt(x::Pair{K,V}...;u...) where{K,V}=ModuleElt(collect(x);u...)
ModuleElt(x::Base.Generator;u...)=ModuleElt(Dict(x);u...)

function ModuleElt(x::Vector{Pair{K,V}};check=false) where{K,V}
  if !check return ModuleElt(Dict(x)) end
  res=Dict{K,V}()
  for (k,v) in x
    if haskey(res,k) res[k]+=v 
    else res[k]=v 
    end
  end
  ModuleElt(res;check=true)
end

# forwarded methods
Base.haskey(x::ModuleElt,y...)=haskey(x.d,y...)
Base.keys(x::ModuleElt)=keys(x.d)
Base.values(x::ModuleElt)=values(x.d)
Base.getindex(x::ModuleElt{K,V},i) where{K,V}=haskey(x,i) ?  x.d[i] : zero(V)

Base.merge(op::Function,a::ModuleElt,b::ModuleElt)::ModuleElt=
  ModuleElt(merge(op,a.d,b.d);check=true)

else
#-------------- faster implementation -------------------------------------
"""
The  type below has a similar interface  to Dict{K,V}, but merge(+,..) is 3
times faster than the `Dict` implementation. It also has the advantage that
ModuleElts are naturally sortable (since type K is assumed sortable).
"""
# The vector d is kept sorted by K 
struct ModuleElt{K,V}
  d::Vector{Pair{K,V}}
  function ModuleElt(d::AbstractVector{Pair{K,V}};check::Bool=false)where {K,V}
    if check && !isempty(d) 
      sort!(d,by=first)
      ri=1
@inbounds for j in 2:length(d)
        if first(d[j])==first(d[ri])
          d[ri]=first(d[ri])=>last(d[ri])+last(d[j])
        else 
          if !iszero(last(d[ri])) ri+=1 end
          d[ri]=d[j]
        end
      end
      if iszero(last(d[ri])) ri-=1 end
      resize!(d,ri)
    end
    new{K,V}(d)
  end
end

@inline ModuleElt(x::Base.Generator;u...)=ModuleElt(collect(x);u...)
@inline Base.pairs(x::ModuleElt)=x.d

@inline Base.cmp(x::ModuleElt,y::ModuleElt)=cmp(x.d,y.d)
@inline Base.isless(x::ModuleElt,y::ModuleElt)=cmp(x,y)==-1

"""
+ is like merge(op,a,b) for Dicts, except keys with value 0 are deleted
"""
function Base.merge(op::Function,a::ModuleElt,b::ModuleElt)
  (a,b)=promote(a,b)
  la=length(a.d)
  lb=length(b.d)
  res=similar(a.d,la+lb)
  ai=bi=1
  ri=0
@inbounds while ai<=la || bi<=lb
    if     ai>la res[ri+=1]=b.d[bi]; bi+=1
    elseif bi>lb res[ri+=1]=a.d[ai]; ai+=1
    else c=cmp(first(a.d[ai]),first(b.d[bi]))
      if     c>0 res[ri+=1]=b.d[bi]; bi+=1
      elseif c<0 res[ri+=1]=a.d[ai]; ai+=1
      else s=op(last(a.d[ai]),last(b.d[bi]))
        if !iszero(s) res[ri+=1]=first(a.d[ai])=>s end
        ai+=1; bi+=1
      end
    end
  end
  ModuleElt(resize!(res,ri))
end

function Base.getindex(x::ModuleElt{K,V},i) where {K,V}
  r=searchsorted(x.d,Ref(i);by=first)
  if r.start!=r.stop return zero(V) end
  x.d[r.start][2]
end

function Base.haskey(x::ModuleElt,i)
  r=searchsorted(x.d,Ref(i);by=first)
  r.start==r.stop
end

@inline Base.keys(x::ModuleElt)=first.(x.d)
@inline Base.values(x::ModuleElt)=last.(x.d)
end
#-------------- methods which have same code in both implementations-------

@inline ModuleElt(x::Pair...;u...)=ModuleElt(collect(x);u...)

@inline Base.:+(a::ModuleElt,b::ModuleElt)=merge(+,a,b)
@inline Base.:-(a::ModuleElt,b::ModuleElt)=merge(-,a,b)

# multiply module element by scalar
function Base.:*(a::ModuleElt{K,V},b)where {K,V}
  if iszero(b) || iszero(a) 
    return zero(ModuleElt{K,promote_type(V,typeof(b))})
  end
  let b=b
    ModuleElt(k=>v*b for (k,v) in a)
  end
end

@inline Base.iszero(x::ModuleElt)=isempty(x.d)
Base.zero(x::ModuleElt)=ModuleElt(empty(x.d))
Base.zero(::Type{ModuleElt{K,V}}) where{K,V}=ModuleElt(Pair{K,V}[])
Base.:-(a::ModuleElt)=iszero(a) ? a : ModuleElt(k=>-v for (k,v) in a)
# forwarded methods
@inline Base.:(==)(a::ModuleElt,b::ModuleElt)=a.d==b.d
@inline Base.first(x::ModuleElt)=first(x.d)
@inline Base.iterate(x::ModuleElt,y...)=iterate(x.d,y...)
@inline Base.length(x::ModuleElt)=length(x.d)
@inline Base.eltype(x::ModuleElt)=eltype(x.d)
@inline Base.hash(x::ModuleElt, h::UInt)=hash(x.d,h)

function Base.convert(::Type{ModuleElt{K,V}},a::ModuleElt{K1,V1}) where {K,K1,V,V1}
  if K==K1
    if V==V1 a
    elseif iszero(a) zero(ModuleElt{K,V})
    else ModuleElt(k=>convert(V,v) for (k,v) in a.d)
    end
  else 
    if iszero(a) zero(ModuleElt{K,V})
    elseif V==V1  ModuleElt(convert(K,k)=>v for (k,v) in a.d)
    else ModuleElt(convert(K,k)=>convert(V,v) for (k,v) in a.d)
    end
  end
end

function Base.promote_rule(a::Type{ModuleElt{K1,V1}},
                           b::Type{ModuleElt{K2,V2}})where {K1,K2,V1,V2}
  ModuleElt{promote_type(K1,K2),promote_type(V1,V2)}
end

# copied from Util in order to need no imports
function format_coefficient(c::String)
  if c=="1" ""
  elseif c=="-1" "-"
  elseif match(r"^[-+]?([^-+*/]|√-|{-)*(\(.*\))?$",c)!==nothing c
  else "("*c*")" 
  end
end

function Base.show(io::IO,m::ModuleElt)
  if iszero(m) print(io,"0"); return end
  showbasis=get(io,:showbasis,nothing)
  if isnothing(showbasis) 
    showbasis=(io,x)->sprint(show,x;context=io)
  end
  start=true
  res=""
  for (k,v) in m 
    v=sprint(show,v;context=io)
    k=showbasis(io,k)
    if !isempty(k) v=format_coefficient(v) end
    if (isempty(v) || v[1]!='-') && !start v="+"*v end
    res*=v*k
    if start start=false end
  end
  if res=="" res="1" end
  print(io,res)
end
#--------------------------------------------------------------------------
end
