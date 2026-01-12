# ways to set and access the Chevie data stored in CHEVIE
module InitChevie

using ..Chevie
export CHEVIE, chevieget, chevieset, InfoChevie, TypeIrred, indices

verbose_chevieget::Bool=false
@GapObj struct _CHEVIE end
const CHEVIE=_CHEVIE(Dict{Symbol,Any}(:info=>true))

function InfoChevie(a...)
  if CHEVIE.info xprint(a...) end
end

" chevieget(t::Symbol,w::Symbol) returns CHEVIE.t[w] or nothing if absent"
function chevieget(t::Symbol,w::Symbol)
  get!(CHEVIE[t],w)do
    if verbose_chevieget println("CHEVIE.$t has no $w") end
  end
end

chevieget(t::String,w::Symbol)=chevieget(Symbol(t),w)

function chevieset(t::Symbol,w::Symbol,o)
  if t==Symbol("3D4") && w==:ReflectionDegrees error() end
  get!(CHEVIE.prop,t,Dict{Symbol,Any}())[w]=o
end

chevieset(t::String,w::Symbol,o)=chevieset(Symbol(t),w,o)

function chevieset(t::Vector,w::Symbol,f::Function)
  for s in t chevieset(s,w,f(s)) end
end

#------------------------------------------------------------------------
@GapObj mutable struct TypeIrred end
@doc """
a  `TypeIrred` object  classifies an  irreducible finite complex reflection
group,  or an  irreducible coset  (the latter  means that  the group  has a
single orbit of irreducible components under `.phi`).

For an irreducible group a `TypeIrred` has the properties:
  - `.rank`  the semisimplerank of the group
  - `.series` which takes one of the values `:A,:B,:D,:E,:F,:G,:H,:I` for an 
      irreducible Coxeter group, and is `:ST` for non-real groups.
  - `.ST` for a primitive non-real group holding the Shepard-Todd number
  - `.p` and `.q` for an imprimitive non-real group, holding `.p=de` and `.q=e`
    for `G(de,e,r)`.

a  `TypeIrred` may  also contain  information specifying  a specific Cartan
matrix  for  the  given  type.  When  there  are  two  conjugacy classes of
generators,  `.cartanType`  (assumed  to  be  `1`  if  this  key is absent)
contains  the ratio  of the  root lengths  compared to  the standard cartan
matrix  for  that  type;  that  is,  the  Cartan matrix is conjugate to the
standard  Cartan matrix by `Diagonal([1,…,1,c,…,c])` where `c=.cartanType`.
This  this  how  type  `B`  (with  `.cartanType==2`)  and  type  `C`  (with
`.cartanType==1`), which have both `.series==:B`,  are distinguished.

For an irreducible coset a `TypeIrred` has the properties
  - `.orbit` a `Vector{TypeIrred}` holding the types of the groups in the
    orbit under `.phi`, such that `.phi` send each item to the next.
  - `.phi` if `k` is the length of `.orbit`, contains the permutation
    effected by `.phi^k` on the simple roots of the first item of the orbit.

In addition, a `TypeIrred` `t` for a group `W` or a `t` appearing in a
`.orbit` for a coset of a group `W` contains a property
   - `.indices` giving the indices in `gens(W)` represented by the generators
     of the irreducible component described by `t`.

The `.indices` of a `TypeIrred` are always in a canonical order for a given component, such that the associated Cartan matrix is the "canonical" one. For
instance, in type `C`, the longest root is always the first one, etc…
""" TypeIrred

function TypeIrred(;kw...)
  t=TypeIrred(Dict(kw...))
  field(t)
  t
end

Base.copy(t::TypeIrred)=TypeIrred(copy(t.prop))

function indices(t::TypeIrred)::Vector{Int}
  if haskey(t,:indices) t.indices
  elseif haskey(t,:orbit) 
    o=t.orbit::Vector{TypeIrred}
    if isempty(o)  Int[] 
    elseif length(o)==1 o[1].indices
    else vcat(indices.(o)...)
    end
  elseif haskey(t,:rank) 1:t.rank
  end
end

indices(t::Vector{TypeIrred})=isempty(t) ? Int[] : length(t)==1 ? indices(t[1]) : vcat(indices.(t)...)

function Symbols.rank(t::TypeIrred)::Int
  if haskey(t,:rank) return t.rank end
  i=indices(t)
  if i!==nothing return length(i) end
end

function Base.show(io::IO, t::TypeIrred)
  function sub(p,n)
    s=string(n)
    string(p)*(length(s)==1 ? "_"*s : "_{"*s*"}")
  end
  if haskey(t,:series)
    s=t.series
    if s==:ST
      if hasdecor(io)
        n="G"
        if haskey(t,:cartanType)
          n*="("*xrepr(io,t.cartanType)*")"
        end
      else n="crg"
      end
      if hasdecor(io)
        if haskey(t,:ST) n=sub(n,t.ST)
        else n*="_{$(t.p),$(t.q),$(t.rank)}"
        end
      else
        if haskey(t,:ST) n*="($(t.ST)"
        else n*="($(t.p),$(t.q),$(t.rank)"
        end
        if haskey(t,:cartanType) n*=","*xrepr(io,t.cartanType)*")"
        else n*=")"
        end
      end
    else
      r=rank(t)
      ct=1
      if haskey(t,:cartanType)
        if s==:B
          if t.cartanType==1 s=:C
          elseif t.cartanType==2 s=:B
          elseif t.cartanType==root(2) s=:Bsym
          else ct=t.cartanType
          end
        elseif s==:F
          if t.cartanType==1 s=:F
          elseif t.cartanType==root(2) s=:Fsym
          else ct=t.cartanType
          end
        elseif s==:G
          if t.cartanType==1 s=:G
          elseif t.cartanType==root(3) s=:Gsym
          else ct=t.cartanType
          end
        elseif s==:I
          if t.cartanType==1 s=:I
          elseif t.cartanType==E(2*t.bond)+E(2*t.bond,-1) s=:Isym
          else ct=t.cartanType
          end
        end
      end
      if hasdecor(io)
        if ct!=1 s=Symbol(s,"(",xrepr(io,t.cartanType),")") end
        if haskey(t,:bond) n=sub(s,r)*"($(t.bond))"
        elseif haskey(t,:short) n="\\tilde "*sub(s,r)
        else n=sub(s,r)
        end
      else
        if haskey(t,:bond) n="coxgroup(:$s,$r,$(t.bond)"
        else n="coxgroup(:$s,$r"
        end
        if ct!=1 n*=","*xrepr(io,t.cartanType)*")"
        else n*=")"
        end
      end
    end
    if hasdecor(io) printTeX(io,n) else print(io,n) end
  else
    o=order(t.twist)
    if hasdecor(io)
      if o!=1 printTeX(io,"{}^{$o}") end
      if length(t.orbit)==1 print(io,t.orbit[1])
      else print(io,"(")
        for t1 in t.orbit print(io,t1) end
        print(io,")")
      end
    else
      print(io,"spets(",t.orbit)
      p=prod(Perm.(map(x->x.indices,t.orbit)...))*t.twist
      if !isone(p) print(io,",",p) end
      print(io,")")
    end
    if haskey(t,:scalar) && !all(isone,t.scalar)
      print(io,"[");join(io,t.scalar,",");print(io,"]")
    end
  end
end

function field(t1::TypeIrred)
  if haskey(t1,:orbit)
    phi=t1.twist::Perm{Perms.Idef}
    orderphi=order(t1.twist)
    t=t1.orbit[1]::TypeIrred
  else
    orderphi=1
    t=t1
  end
  s=t.series::Symbol
  ff=if s==:ST 
    if haskey(t,:ST)
      if orderphi!=1 
        [Symbol(orderphi,"G",t.ST)]
      elseif 4<=t.ST<=22 
           [:G4_22,t.ST]
      else [Symbol("G",t.ST)]
      end
    elseif orderphi!=1 
         [:timp, t.p, t.q, t.rank, phi]
    else [:imp, t.p, t.q, t.rank]
    end
  elseif s==:I 
    [orderphi==1 ? s : Symbol(2,s),t.bond]
  elseif s in [:A,:B,:D] 
     if orderphi==1 
       [s,rank(t)]
     elseif orderphi==2 
       if s==:B 
            [Symbol(2,:I),4]
       else [Symbol(2,s),rank(t)]
       end
     elseif orderphi==3 
       [Symbol("3D4")]
     end
  elseif s in [:E,:F,:G]
    if orderphi==1 
      [Symbol(s,rank(t))]
    elseif s==:G 
         [Symbol(2,:I),6]
    else [Symbol(orderphi,s,rank(t))] 
    end
  else [Symbol(s,rank(t))] 
  end
  t1.cheviefile=ff[1]
  t1.extra=ff[2:end]
end

# functions whose result depends on cartanType
const needcartantype=Set([:Invariants,
                          :PrintDiagram,
                          :UnipotentClasses,
                          :WeightInfo,
                          :cartan])

debug::Bool=false # if true time each call

"`chevieget(t::TypeIrred,f::Symbol,extra...)` get `CHEVIE[field(t)][f](extra...)`"
function chevieget(t::TypeIrred,f::Symbol,extra...)
  #println("t=$t f=$f extra=$extra")
  o=chevieget(t.cheviefile,f)
  textra=t.extra
  if o isa Function
    if haskey(t,:orbit) t=t.orbit[1] end
if debug
    if haskey(t,:cartanType) && f in needcartantype
      println("$n.$f(",(t.extra...,extra...,t.cartanType),")")
@time o(textra...,extra...,t.cartanType)
    else 
     println("$n.$f(",(t.extra...,extra...),")")
@time o(textra...,extra...)
    end
else
    if haskey(t,:cartanType) && f in needcartantype
      o(textra...,extra...,t.cartanType)
    else 
      o(textra...,extra...)
    end
end
  elseif o===false return nothing
  else o
  end
end

end
