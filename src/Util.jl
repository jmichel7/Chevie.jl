"""
This  module contains  various utility  functions used  in the  rest of the
code.  Maybe some  of them  exist in  some Julia  module I am not aware of;
please tell me.
"""
module Util
export ds

using ..Format: xprintln

function ds(s;types=false) # "dump struct"; not recursive like dump
  println(typeof(s),":")
  for f in fieldnames(typeof(s))
    if !isdefined(s,f) println(f,"=#undef")
    else print(f);
      if types print("::",typeof(getfield(s,f))) end
      xprintln("=",getfield(s,f))
    end
  end
end

export @forward 
"""
`@forward T.f f1,f2,...`

is a macro which delegates definitions. The above generates 
```
f1(a::T,args...)=f1(a.f,args...)
f2(a::T,args...)=f2(a.f,args...)
...
```
"""
macro forward(ex, fs)
  T, field = esc(ex.args[1]), ex.args[2].value
  fdefs=map(fs.args)do ff
      f= esc(ff)
      quote
        ($f)(a::($T),args...)=($f)(a.$field,args...)
      end
  end
  Expr(:block, fdefs...)
end

#---------------- high level Chevie compatibility -------------------------
export toL, toM # convert Gap matrices <-> Julia matrices
toL(m)=collect(eachrow(m)) # to Gap
toM(l)=isempty(l) ? Array{eltype(eltype(l))}(undef,0,1) : stack(l;dims=1)

# The following functions should be eventually obsoleted by adopting Julia
# order for products.

export cartesian, tcartesian, cart2lin, lin2cart
"""
`cartesian(a::AbstractVector...)`

A variation on `Iterators.product` which gives the same result as GAP's
`Cartesian` (thanks to `permutedims`)
"""
function cartesian(a::AbstractVector...)
  if isempty(a) return [] end
  vec(permutedims(collect.(Iterators.product(a...)),reverse(eachindex(a))))
end

# faster version returning tuples
function tcartesian(a::AbstractVector...)
  if isempty(a) return Tuple{Any}[] end
  vec(permutedims(collect(Iterators.product(a...)),reverse(eachindex(a))))
end

"""
`cart2lin([l₁,…,lₙ],[i₁,…,iₙ])` is GAP3 PositionCartesian
faster findfirst(==([i1,…,iₙ]),cartesian(1:l₁,…,1:lₙ))
very fast with 2 Tuple arguments
"""
cart2lin(l,ind)=LinearIndices(Tuple(reverse(l)))[reverse(ind)...]

"""
`lin2cart([l₁,…,lₙ],i)` is GAP3 CartesianAt
faster cartesian(map(j->1:j,l))[i]
"""
lin2cart(dims,i)=reverse(Tuple(CartesianIndices(reverse(Tuple(dims)))[i]))
end
