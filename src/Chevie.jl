module Chevie

using ..Gapjm
export CHEVIE, chevieget, chevieset, getchev

const CHEVIE=Dict{Symbol,Any}(
 :compat=>Dict(:MakeCharacterTable=>x->x,
               :AdjustHeckeCharTable=>(x,y)->x,
               :ChangeIdentifier=>function(tbl,n)tbl[:identifier]=n end),
 :info=>false
)

function chevieget(t::Symbol,w::Symbol)
  if haskey(CHEVIE[t],w) return CHEVIE[t][w] end
  if CHEVIE[:info] println("CHEVIE[$t] has no $w") end
end

chevieget(t::String,w::Symbol)=chevieget(Symbol(t),w)

function chevieset(t::Symbol,w::Symbol,o::Any)
  if !haskey(CHEVIE,t) CHEVIE[t]=Dict{Symbol,Any}() end
  CHEVIE[t][w]=o
end

function chevieset(t::Vector{String},w::Symbol,f::Function)
  for s in t 
    println("set $s $w")
    chevieset(Symbol(s),w,f(s)) 
  end
end

function field(t::TypeIrred)
  if haskey(t,:orbit)
    orderphi=order(t.twist)
    t=t.orbit[1]
  else
    orderphi=1
  end
  s=t.series
  if s in [:A,:B,:D] 
     if orderphi==1 return (s,PermRoot.rank(t))
     elseif orderphi==2 
       if s==:B return (Symbol("2I"),4)
       else return (Symbol(2,s),PermRoot.rank(t))
       end
     elseif orderphi==3 return (Symbol("3D4"),)
     end
  elseif s in [:E,:F,:G]
    if orderphi==1 return (Symbol(s,PermRoot.rank(t)),) 
    elseif s==:G return (Symbol("2I"),6)
    else return (Symbol(orderphi,s,PermRoot.rank(t)),) 
    end
  elseif s==:ST 
    if haskey(t,:ST)
      if orderphi!=1 return (Symbol(orderphi,"G",t.ST),)
      elseif 4<=t.ST<=22 return (:G4_22,t.ST)
      else return (Symbol(string("G",t.ST)),)
      end
    elseif orderphi!=1
      return (:timp, t.p, t.q, t.rank)
    else
      return (:imp, t.p, t.q, t.rank)
    end
  elseif s==:I return (orderphi==1 ? :I : Symbol("2I"),t.bond)
  else return (Symbol(string(s,PermRoot.rank(t))),) 
  end
end

const needcartantype=Set([:Invariants,
                          :PrintDiagram,
                          :ReflectionName,
                          :UnipotentClasses,
                          :WeightInfo,
                          :CartanMat])

function getchev(t::TypeIrred,f::Symbol,extra...)
  d=field(t)
# println("d=$d f=$f extra=$extra")
  o=chevieget(d[1],f)
  if o isa Function
#   o(vcat(collect(d)[2:end],collect(extra))...)
    if haskey(t,:orbit) t=t.orbit[1] end
    if haskey(t,:cartanType) && f in needcartantype
#     println("args=",(d[2:end]...,extra...,t.cartanType))
      o(d[2:end]...,extra...,t.cartanType)
     else o(d[2:end]...,extra...)
    end
  else o
  end
end

function getchev(W,f::Symbol,extra...)
  [getchev(ti,f::Symbol,extra...) for ti in refltype(W)]
end

end
