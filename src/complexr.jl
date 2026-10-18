export complex_reflection_group, crg, reflection_group

"""
`complex_reflection_group(STnumber[,cartanType=1])`

`complex_reflection_group(p,q,r[,cartanType=1])`

`complex_reflection_group` can be abbreviated as `crg`.

The  first form of `complex_reflection_group` returns the primitive complex
reflection   group   which   has   Shephard-Todd   number  `STnumber`,  see
[st54](@cite).  The second form returns  the imprimitive complex reflection
group named `G(p,q,r)` by Shephard-Todd.

If  a `cartanType`  different from  1 is  given and  there are  2 orbits of
hyperplanes  in the resulting group `W`, then the ratio of the root lengths
between the two orbits is modified by `cartanType'.

```julia-repl
julia> G=complex_reflection_group(4)
G₄

julia> degrees(G)
2-element Vector{Int64}:
 4
 6

julia> length(G)
24

julia> G*coxgroup(:A,2) # how to make a non-irreducible group
G₄×A₂

julia> complex_reflection_group(1,1,3) # another way to enter rootdatum(:gl,3) 
gl₃

julia> crg(5,3) # the ratio of the lengths is 3
G(3)₅
```
"""
function complex_reflection_group(i::Integer,cartanType=1)
  if i==23     return coxgroup(:H,3)
  elseif i==28 return coxgroup(:F,4,cartanType)
  elseif i==30 return coxgroup(:H,4)
  elseif i==35 return coxgroup(:E,6)
  elseif i==36 return coxgroup(:E,7)
  elseif i==37 return coxgroup(:E,8)
  end
  t=TypeIrred(;series=:ST,ST=Int(i))
  r=simpleroots(t)
  if cartanType==1 cr=simplecoroots(t)
  else
    if i>23
      if i!=26 error("G",i," has no cartanType")
      else     error("cartanType!=1 not implemented for G26")
      end
    end
    t.cartanType=cartanType
    c=cartan(t)
    cr=solutionmat(transpose(r),c)
    if nothing in cr error("unexpected") end
  end
  res=PRG(r,cr;NC=true)
  res.length=prod(degrees(t)) #speeds up refltype
  refltype(res)
  res
end

function complex_reflection_group(p,q,r,cartanType=1)
  if !iszero(p%q) || p<=0 || r<=0 || (r==1 && q!=1)
   error("complex_reflection_group(p,q,r) must satisfy: q|p, r>0, and r=1 => q=1")
  end
  if p==1 return r==1 ? coxgroup() : rootdatum(:gl,r)
  elseif p==2 return rootdatum(:so,q==2 ? 2r : 2r+1)
  elseif p==q && r==2 return coxgroup(:I,2,p)
  end
  t=TypeIrred(;series=:ST,p=Int(p),q=Int(q),rank=Int(r))
  rr=simpleroots(t)
  if cartanType==1 cr=simplecoroots(t)
  else
    t.cartanType=cartanType
    c=cartan(t)
    cr=solutionmat(transpose(rr),c)
    if nothing in cr error("unexpected") end
  end
  res=PRG(rr,cr;NC=true)
  res.length=prod(degrees(t)) #speeds up refltype
  refltype(res)
  res
end

const crg=complex_reflection_group

# converts a type back to a group
function reflection_group(t::TypeIrred)
  if haskey(t,:orbit)
    W=reflection_group(t.orbit)
    if length(t.orbit)>1
      spets(W,reflrep(W,Perm(vcat(circshift(map(x->x.indices,
                                        refltype(W)),-1)...))*t.twist))
    else spets(W,t.twist)
    end
  elseif t.series==:ST PRG(simpleroots(t),simplecoroots(t))
  else C=cartan(t)
    all(isreal,C) ? rootdatum(C) : PRG(one(C),C)
  end
end

function reflection_group(l::Vector{TypeIrred})
  if isempty(l) return coxgroup() end
  prod(reflection_group.(l))
end
