"""
Module Elements --- elements of free modules.

A  `ModuleElt{K,V}`  represents  an  element  of  a free module where basis
elements  are of type  `K` and coefficients  of type `V`.  Usually you want
objects  of type `V` to be elements of  a ring, but it could also be useful
if  they just belong to  an abelian group. This  is similar to the SageMath
CombinatorialFreeModule.

This  basic data structure is  used in the package  `Gapjm` as an efficient
representation  at  many  places.  For  example, multivariate monomials are
represented by `ModuleElt{Symbol,Int}`:

`x^2y^-3 ` is represented by `ModuleElt(:x=>2,:y=>-3)`

And  multivariate  polynomials  are  represented by `ModuleElt{Monomial,C}`
where `C` is the type of the coefficients:

`x*y-z^2` is represented by ``ModuleElt(x*y=>1,z^2=>-1)

`ModuleElts`  are  also  used  for  cyclotomics, CycPols, elements of Hecke
algebras, etc…

A  `ModuleElt{K,V}` is essentially a  list of `Pairs{K,V}`. The constructor
takes  as argument a list of pairs, or a variable number of pair arguments,
or a generator of pairs.

We provide two implementations:

  - an implementation by `Dict`s 

This  requires  that  the  type  `K`  is  hashable.  It  is  a  very simple
implementation  since the interface of the type  is close to that of dicts;
the  only difference is weeding  out keys which have  a zero cofficient ---
which  is necessary since for testing equality of module elements one needs
a canonical form for each element.

  - a faster implementation (the default) by  keeping a list of pairs sorted
by  key.  This  demands  that  the  type  `K`  has  a `isless` method. This
implementation  is  two  to  four  times  faster  than  the `Dict` one (for
addition, the most important operation) and requires half the memory.

A  ̀ModuleElt` has mostly the same methods as a `Dict`. Adding `ModuleElt`s
is  a variation on `merge(+,...)` for `Dict`s  (here `+` can be replaced by
any operation `op` with the property that `op(0,x)=x`) where keys with zero
value  are deleted  afterwards. Further,  a `ModuleElt`  can be negated, or
multiplied or divided (`/`or `//`) by some element (acting on coefficients)
if  the method is defined between type `V` and that element; there are also
`zero` and `iszero` methods.

Here  is an example where basis elements are `Symbol`s and coefficients are
`Int`.

```julia-repl
julia> a=ModuleElt(:xy=>1,:yx=>-1)
:xy-:yx

julia> a-a
0

julia> a*99
99:xy-99:yx

julia> a//2
(1//2):xy+(-1//2):yx

julia> a/2
0.5:xy-0.5:yx

julia> a+ModuleElt(:yx=>1)
:xy

julia> a[:xy] # indexing by a basis element finds the coefficient
1

julia> a[:xx] # the coefficient of an absent basis element is zero.
0

julia> haskey(a,:xx)
false

julia> first(a)
:xy => 1

julia> collect(a)
2-element Vector{Pair{Symbol, Int64}}:
 :xy => 1
 :yx => -1

julia> collect(keys(a))
2-element Vector{Symbol}:
 :xy
 :yx

julia> collect(values(a))
2-element Vector{Int64}:
  1
 -1

julia> length(a)
2

julia> eltype(a)
Pair{Symbol, Int64}
```

In both implementations the constructor normalizes the constructed element,
removing zero coefficients and merging duplicate basis elements, adding the
corresponding   coefficients  (and   sorting  the   basis  in  the  default
implementation).  If  you  know  this  normalisation is unnecessary, to get
maximum  speed you can disable this  by giving the keyword `check=false` to
the constructor.

```julia-repl
julia> a=ModuleElt(:yy=>1, :yx=>2, :xy=>3, :yy=>-1;check=false)
:yy+2:yx+3:xy-:yy

julia> a=ModuleElt(:yy=>1, :yx=>2, :xy=>3, :yy=>-1)
3:xy+2:yx
```

As  you can see in the above examples, at the REPL (or in Jupyter or Pluto,
when  `IO`  has  the  `:limit`  attribute)  the  `show`  method  shows  the
coefficients (bracketed if necessary, which is when they have inner occurences
of  `+-*`),  followed  by  showing  the  basis  elements.  Setting the `IO`
property  `:showbasis` to a custom printing  function changes how the basis
elements are printed.

```julia-rep1
julia> show(IOContext(stdout,:showbasis=>(io,s)->string("<",s,">")),a)
3<xy>+2<yx>
```

The `repr` method gives a representation which can be read back in julia:

```julia-repl
julia> repr(a)
"ModuleElt([:xy => 3, :yx => 2])"
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
  function ModuleElt(d::Dict{K,V};check::Bool=true)where {K,V}
    if check
      for (k,v) in d if iszero(v) delete!(d,k) end end
    end
    new{K,V}(d)
  end
end

ModuleElt(x::Pair{K,V}...;u...) where{K,V}=ModuleElt(collect(x);u...)
ModuleElt(x::Base.Generator;u...)=ModuleElt(Dict(x);u...)

function ModuleElt(x::Vector{Pair{K,V}};check=true) where{K,V}
  if !check || length(x)<=1 return ModuleElt(Dict(x)) end
  res=Dict{K,V}()
  for (k,v) in x
    if haskey(res,k) res[k]+=v 
    else res[k]=v 
    end
  end
  ModuleElt(res)
end

# forwarded methods
Base.haskey(x::ModuleElt,y...)=haskey(x.d,y...)
Base.keys(x::ModuleElt)=keys(x.d)
Base.values(x::ModuleElt)=values(x.d)

Base.getindex(x::ModuleElt{K,V},i) where{K,V}=haskey(x,i) ?  x.d[i] : zero(V)

# Valid for ops such that op(0,x)=x otherwise the result is wrong.
Base.merge(op::Function,a::ModuleElt,b::ModuleElt)::ModuleElt=
  ModuleElt(merge(op,a.d,b.d))

else
#-------------- faster implementation -------------------------------------
"""
The  type  below  has  a  similar  interface to `Dict{K,V}`, but instead of
assuming  that `K` is hashable, it assumes that `K` is sortable. With this,
`merge(+,..)` is 3 times faster than the `Dict` implementation. It also has
the advantage that ModuleElts are naturally sortable.

The unique field, a `Vector{Pair{K,V}}`, is kept sorted by `K`; 
the constructor checks sorting by default if `length(d)>1`. 
This can be bypassed by setting `check=false`.
"""
struct ModuleElt{K,V}
  d::Vector{Pair{K,V}}
  function ModuleElt(d::AbstractVector{Pair{K,V}};check::Bool=true)where {K,V}
    if check && length(d)>0
      if length(d)>1 sort!(d,by=first) end
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

# Like  merge(op,a,b) for  Dicts, except  keys with  value 0 are deleted. 
# Valid for commutative ops such that op(0,x)=x otherwise the result is wrong.
function Base.merge(op::Function,a::ModuleElt,b::ModuleElt)
# this version is correct only if op(0,x)==x always
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
  ModuleElt(resize!(res,ri);check=false)
end

# The  version below works for general  ops (not commutative or which don't
# satisfy always op(0,x)=x), but has too much overhead currently to replace
# the  above version for + or other  ops such that op(0,x)==x. It is useful
# for max or min which do lcm and gcd of `Monomial`s or `CycPol`s.
function merge2(op::Function,a::ModuleElt,b::ModuleElt)
  (a,b)=promote(a,b)
  la=length(a.d)
  lb=length(b.d)
  res=similar(a.d,la+lb)
  ai=bi=1
  ri=0
@inbounds while ai<=la || bi<=lb
    if     ai>la 
      be=b.d[bi]
      s=op(zero(last(be)),last(be))
      if !iszero(s) res[ri+=1]=first(be)=>s end
      bi+=1
    elseif bi>lb 
      ae=a.d[ai]
      s=op(last(ae),zero(last(ae)))
      if !iszero(s) res[ri+=1]=first(ae)=>s end
      ai+=1
    else c=cmp(first(a.d[ai]),first(b.d[bi]))
      if     c>0 
        be=b.d[bi]
        s=op(zero(last(be)),last(be))
        if !iszero(s) res[ri+=1]=first(be)=>s end
        bi+=1
      elseif c<0 
        ae=a.d[ai]
        s=op(last(ae),zero(last(ae)))
        if !iszero(s) res[ri+=1]=first(ae)=>s end
        ai+=1
      else s=op(last(a.d[ai]),last(b.d[bi]))
        if !iszero(s) res[ri+=1]=first(a.d[ai])=>s end
        ai+=1; bi+=1
      end
    end
  end
  ModuleElt(resize!(res,ri);check=false)
end

"""
`getindex(x::ModuleElt,i)` returns the value associated to key `i` in `x`.
It returns zero if the key does not occur in `x`.
"""
function Base.getindex(x::ModuleElt{K,V},i) where {K,V}
  r=searchsorted(x.d,Ref(i);by=first)
  if r.start!=r.stop return zero(V) end
  last(x.d[r.start])
end

function Base.haskey(x::ModuleElt,i)
  r=searchsorted(x.d,Ref(i);by=first)
  r.start==r.stop
end

Base.keys(x::ModuleElt)=(first(p) for p in x.d)
Base.values(x::ModuleElt)=(last(p) for p in x.d)
end
#-------------- methods which have same code in both implementations-------

@inline ModuleElt(x::Pair...;u...)=ModuleElt(collect(x);u...)

@inline Base.:+(a::ModuleElt,b::ModuleElt)=merge(+,a,b)
Base.:-(a::ModuleElt)=iszero(a) ? a : ModuleElt(k=>-v for(k,v) in a;check=false)
@inline Base.:-(a::ModuleElt,b::ModuleElt)=a+(-b)

# multiply module element by scalar
function Base.:*(a::ModuleElt{K,V},b)where {K,V}
  if iszero(b) || iszero(a) 
    return zero(ModuleElt{K,promote_type(V,typeof(b))})
  end
  let b=b
    ModuleElt(k=>v*b for (k,v) in a;check=false)
  end
end

# divide by scalar
function Base.:/(a::ModuleElt{K,V},b)where {K,V}
  if iszero(a) return zero(ModuleElt{K,promote_type(V,typeof(b))}) end
  let b=b
   ModuleElt(k=>v/b for (k,v) in a;check=false)
  end
end

# divide by scalar
function Base.:(//)(a::ModuleElt{K,V},b)where {K,V}
  if iszero(a) return zero(ModuleElt{K,promote_type(V,typeof(b))}) end
  let b=b
   ModuleElt(k=>v//b for (k,v) in a;check=false)
  end
end

@inline Base.iszero(x::ModuleElt)=isempty(x.d)
Base.zero(x::ModuleElt)=ModuleElt(empty(x.d))
Base.zero(::Type{ModuleElt{K,V}}) where{K,V}=ModuleElt(Pair{K,V}[])
# forwarded methods
@inline Base.:(==)(a::ModuleElt,b::ModuleElt)=a.d==b.d
@inline Base.first(x::ModuleElt)=first(x.d)
@inline Base.iterate(x::ModuleElt,y...)=iterate(x.d,y...)
@inline Base.length(x::ModuleElt)=length(x.d)
@inline Base.eltype(x::ModuleElt)=eltype(x.d)
@inline Base.hash(x::ModuleElt, h::UInt)=hash(x.d,h)

# we assume that converting the keys from K to K1 does not change hashing
# for the Dict implementation or sorting for the other implementation
function Base.convert(::Type{ModuleElt{K,V}},a::ModuleElt{K1,V1}) where {K,K1,V,V1}
  if K==K1
    if V==V1 a
    elseif iszero(a) zero(ModuleElt{K,V})
    else ModuleElt(k=>convert(V,v) for (k,v) in a.d;check=false)
    end
  else 
    if iszero(a) zero(ModuleElt{K,V})
    elseif V==V1  ModuleElt(convert(K,k)=>v for (k,v) in a.d;check=false)
    else ModuleElt(convert(K,k)=>convert(V,v) for (k,v) in a.d;check=false)
    end
  end
end

function Base.promote_rule(a::Type{ModuleElt{K1,V1}},
                           b::Type{ModuleElt{K2,V2}})where {K1,K2,V1,V2}
  ModuleElt{promote_type(K1,K2),promote_type(V1,V2)}
end

# copied from Util.jl in order to have a self-contained file.
# Tries to determine which coefficients should be bracketed for unambiguous
# display
function format_coefficient(c::String)
  if c=="1" ""
  elseif c=="-1" "-"
  elseif match(r"^[-+]?([^-+*/]|√-|{-)*(\(.*\))?$",c)!==nothing c
  else "("*c*")" 
  end
end

function Base.show(io::IO,m::ModuleElt)
  if !(get(io,:limit,false) || get(io,:TeX,false))
    print(io,"ModuleElt(",m.d,")")
    return
  end
  if iszero(m) print(io,"0"); return end
  showbasis=get(io,:showbasis,nothing)
  if isnothing(showbasis) 
    showbasis=(io,x)->repr(x;context=io)
  end
  start=true
  res=""
  for (k,v) in m 
    v=repr(v;context=io)
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
