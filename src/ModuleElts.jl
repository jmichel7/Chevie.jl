"""
A  `ModuleElt` represents an element of a  module. It is essentially a list
of pairs `b=>c` where `b` is a basis element and `c` its coefficient.

The  constructors  take  as  argument  a  list  of  pairs,  or several pair
arguments,  or a generator of pairs (these pairs should be sorted by key in
the second implementation).

Here is an example where basis elements are represented by Symbols.
The first instruction is to use the show method defined below.

```julia-repl
julia> a=ModuleElt(:xy=>1,:yx=>-1)
:xy-:yx

julia> a-a
0

julia> a*99
99:xy-99:yx

julia> push!(a,:yy=>2)
:xy-:yx+2:yy

julia> a+ModuleElt(:yx=>1)
:xy+2:yy

julia> a[:xy]
1

julia> haskey(a,:xx)
false
```

both implementations provide a function `norm!` which normalizes an
element, removing zero coefficients

```julia-repl
julia> a=ModuleElt(:yy=>1, :yx=>2, :xy=>3, :yy=>-1)
:yy+2:yx+3:xy-:yy

julia> norm!(a)
3:xy+2:yx

julia> a
3:xy+2:yx
```

setting the property showbasis determines how the basis elements are printed.
```julia-rep1
julia> show(IOContext(stdout,:showbasis=>(io,s)->String(s)),a)
3xy+2yx
```

"""
module ModuleElts

using Gapjm

export ModuleElt, norm! # data structure
#------------- implementation with Dicts ----------------------
const usedict=false
if usedict
struct ModuleElt{K,V}
  d::Dict{K,V} # the keys K should be hashable
end

ModuleElt(x::Vector{Pair{K,V}}) where{K,V}=ModuleElt(Dict(x))
ModuleElt(x::Pair{K,V}...) where{K,V}=ModuleElt(Dict(x...))
ModuleElt(x::Base.Generator)=ModuleElt(Dict(x))

Base.zero(::Type{ModuleElt{K,V}}) where{K,V}=ModuleElt(Dict{K,V}())

# forwarded methods
@inline Base.haskey(x::ModuleElt,y...)=haskey(x.d,y...)
@inline Base.getindex(x::ModuleElt,i)=getindex(x.d,i)
@inline Base.:(==)(x::ModuleElt,y::ModuleElt)= x.d == y.d

# multiply module element by scalar
Base.:*(a::ModuleElt,b)=iszero(b) ? zero(a) : ModuleElt(k=>v*b for (k,v) in a)

function norm!(a::ModuleElt)
  for (k,v) in a.d if iszero(v) delete!(a.d,k) end end
  a
end

Base.:+(a::ModuleElt,b::ModuleElt)::ModuleElt=norm!(ModuleElt(merge(+,a.d,b.d)))

else
#-------------- faster implementation -------------------------------------
"""
The  type below  has a  similar interface  to Dict{K,V},  but +  is 3 times
faster than merge(+,...) on Dicts.
"""
# The vector d is kept sorted by K 
struct ModuleElt{K,V}
  d::Vector{Pair{K,V}} # the keys K should be sortable
end

ModuleElt(x::Pair{K,V}...) where{K,V}=ModuleElt(collect(x))
ModuleElt(x::Base.Generator)=ModuleElt(collect(x))
# note: the constructors do not check sorting.

Base.zero(::Type{ModuleElt{K,V}}) where{K,V}=ModuleElt(Pair{K,V}[])
@inline Base.cmp(x::ModuleElt,y::ModuleElt)=cmp(x.d,y.d)

# multiply module element by scalar
Base.:*(a::ModuleElt,b)=iszero(b) ? zero(a) : ModuleElt(k=>v*b for (k,v) in a)

"""
+ is like merge(+,..) for Dicts, except keys with value 0 are deleted
"""
function Base.:+(a::ModuleElt,b::ModuleElt)::ModuleElt
  la=length(a.d)
  lb=length(b.d)
  res=similar(a.d,la+lb)
  ai=bi=1
  ri=0
  while ai<=la || bi<=lb
    if ai>la
@inbounds res[ri+=1]=b.d[bi]; bi+=1
    elseif bi>lb
@inbounds res[ri+=1]=a.d[ai]; ai+=1
    else
@inbounds c=cmp(a.d[ai][1],b.d[bi][1])
      if c==1
@inbounds res[ri+=1]=b.d[bi]; bi+=1
      elseif c==-1
@inbounds res[ri+=1]=a.d[ai]; ai+=1
      else s=a.d[ai][2]+b.d[bi][2]
        if !iszero(s)
@inbounds res[ri+=1]=a.d[ai][1]=>s
        end
        ai+=1; bi+=1
      end
    end
  end
  ModuleElt(resize!(res,ri))
end

"""
normalize an unsorted ModuleElt
"""
function norm!(x::ModuleElt{K,V}) where {K,V}
  if isempty(x) return x end
  sort!(x.d,by=first)
  ri=1
@inbounds  for j in 2:length(x.d)
    if x.d[j][1]==x.d[ri][1]
      x.d[ri]=x.d[ri][1]=>x.d[ri][2]+x.d[j][2]
    else 
      if !iszero(x.d[ri][2]) ri+=1 end
      x.d[ri]=x.d[j]
    end
  end
  if iszero(x.d[ri][2]) ri-=1 end
  ModuleElt(resize!(x.d,ri))
end

function Base.getindex(x::ModuleElt,i)
  r=searchsorted(x.d,Ref(i);by=first)
  if r.start!=r.stop error("key $i not found") end
  x.d[r.start][2]
end

function Base.haskey(x::ModuleElt,i)
  r=searchsorted(x.d,Ref(i);by=first)
  r.start==r.stop
end

# similar to delete! but return nothing if key not found, and
# (copy with deleted key, deleted value) otherwise
function drop(m::ModuleElt,k)
  r=searchsorted(m.d,Ref(k);by=first)
  if r.start!=r.stop return nothing end
  ModuleElt(deleteat!(copy(m.d),r.start)),m.d[r.start][2]
end

Base.keys(x::ModuleElt)=first.(x.d)
Base.first(x::ModuleElt)=first(x.d)
end
#-------------- methods which have same code in both implementations-------
Base.iszero(x::ModuleElt)=isempty(x.d)
Base.zero(x::ModuleElt)=ModuleElt(empty(x.d))
Base.:-(a::ModuleElt)=ModuleElt(k=>-v for (k,v) in a)
Base.:-(a::ModuleElt,b::ModuleElt)=a+(-b)
# forwarded methods
@inline Base.iterate(x::ModuleElt,y...)=iterate(x.d,y...)
@inline Base.length(x::ModuleElt)=length(x.d)
@inline Base.eltype(x::ModuleElt)=eltype(x.d)
@inline function Base.push!(x::ModuleElt,y...)
  push!(x.d,y...)
  x
end
@inline function Base.append!(x::ModuleElt,y...)
  append!(x.d,y...)
  x
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
    if occursin(r"[-+*]",v[nextind(v,0,2):end]) v="($v)" end
    v= (v=="1" ? "" : (v=="-1" ? "-" : v))
    if (isempty(v) || v[1]!='-') && !start v="+"*v end
    res*=v
    res*=showbasis(io,k)
    if start start=false end
  end
  print(io,res)
end
#--------------------------------------------------------------------------
end
