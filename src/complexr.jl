export complex_reflection_group, crg, reflection_group

"""
`complex_reflection_group(STnumber)` or `crg(STnumber)`

`complex_reflection_group(p,q,r)` or `crg(p,q,r)`

The first form of `complex_reflection_group` returns the complex reflection
group  which has  Shephard-Todd number  `STnumber`, see  [st54](@cite). The
second form returns the imprimitive complex reflection group `G(p,q,r)`.

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

julia> complex_reflection_group(1,1,3) # another way to enter A₂
gl₃

julia> crg(4) # there is also a short alias
G₄
```
"""
function complex_reflection_group(i::Integer,cartanType=1)
  if i==23     coxgroup(:H,3)
  elseif i==28 coxgroup(:F,4,cartanType)
  elseif i==30 coxgroup(:H,4)
  elseif i==35 coxgroup(:E,6)
  elseif i==36 coxgroup(:E,7)
  elseif i==37 coxgroup(:E,8)
  elseif cartanType==1
    t=TypeIrred(series=:ST,ST=Int(i))
    PRG(simpleroots(t),simplecoroots(t))
  else
    if i in (4,8,12,16,20,22) error("G",i," has no cartanType") end
    t=TypeIrred(;series=:ST,ST=Int(i))
    r=simpleroots(t)
    c=cartan(t)
    t.cartanType=cartanType
    if i in (5,6,9,10,14,17,18,21)
      D=Diagonal([cartanType,1])
      cr=solutionmat(transpose(r),D*c*inv(D//1))
      if nothing in cr error("unexpected") end
      PRG(r,cr)
    else error("cartanType!=1 not implemented for G",i)
    end
  end
end

function complex_reflection_group(p,q,r,cartanType=1)
  if !iszero(p%q) || p<=0 || r<=0 || (r==1 && q!=1)
   error("complex_reflection_group(p,q,r) must satisfy: q|p, r>0, and r=1 => q=1")
  end
  if p==1 return r==1 ? coxgroup() : rootdatum(:gl,r)
  elseif p==2 return rootdatum(:so,q==2 ? 2r : 2r+1)
  elseif p==q && r==2 return coxgroup(:I,2,p)
  end
  if cartanType!=1
#   if r>2 error("cartanType!=1 not implemented for G",(p,q,r)) end
    if r==1 || p==q  error("G",(p,q,r)," has no cartanType") end
    t=TypeIrred(;series=:ST,p=Int(p),q=Int(q),r=Int(r),rank=r)
    rr=simpleroots(t)
    c=cartan(t)
    t.cartanType=cartanType
    D=Diagonal(vcat([1],fill(cartanType,size(c,1)-1)))
    cr=solutionmat(transpose(rr),inv(D//1)*c*D)
    if nothing in cr error("unexpected") end
    PRG(rr,cr)
  else
    t=TypeIrred(;series=:ST,p,q,rank=r)
    PRG(simpleroots(t),simplecoroots(t))
  end
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
