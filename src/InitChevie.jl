# ways to set and access the Chevie data stored in CHEVIE
module InitChevie

using ..Chevie
export CHEVIE, chevieget, chevieset, InfoChevie

const CHEVIE=Dict{Symbol,Any}(:info=>true,:infoget=>false)

function InfoChevie(a...)
  if CHEVIE[:info] xprint(a...) end
end

" chevieget(t,w) returns CHEVIE[Symbol(t)][w] or nothing if absent"
function chevieget(t::Symbol,w::Symbol)
# println("chevieget(",t,",",w,")")
  get!(CHEVIE[t],w)do
    if CHEVIE[:infoget] println("CHEVIE[$t] has no $w") end
  end
end

chevieget(t::String,w::Symbol)=chevieget(Symbol(t),w)

chevieset(t::Symbol,w::Symbol,o)=get!(CHEVIE,t,Dict{Symbol,Any}())[w]=o

function chevieset(t::Vector{String},w::Symbol,f::Function)
# println(join(t,",")," $w")
  for s in t chevieset(Symbol(s),w,f(s)) end
end

function chevieset(t::Vector{Symbol},w::Symbol,f::Function)
# println(join(t,",")," $w")
  for s in t chevieset(s,w,f(s)) end
end

function field(t::TypeIrred)
  get!(t,:field)do
  if haskey(t,:orbit)
    phi=t.twist
    orderphi=order(t.twist)
    t=t.orbit[1]
  else
    orderphi=1
  end
  s=t.series
  if s==:ST 
    if haskey(t,:ST)
      if orderphi!=1 return (Symbol(orderphi,"G",t.ST),)
      elseif 4<=t.ST<=22 return (:G4_22,t.ST)
      else return (Symbol("G",t.ST),)
      end
    elseif orderphi!=1
      return (:timp, t.p, t.q, t.rank, phi)
    else
      return (:imp, t.p, t.q, t.rank)
    end
  elseif s==:I return (orderphi==1 ? s : Symbol(2,s),t.bond)
  elseif s in [:A,:B,:D] 
     if orderphi==1 return (s,PermRoot.rank(t))
     elseif orderphi==2 
       if s==:B return (Symbol(2,:I),4)
       else return (Symbol(2,s),PermRoot.rank(t))
       end
     elseif orderphi==3 return (Symbol("3D4"),)
     end
  elseif s in [:E,:F,:G]
    if orderphi==1 return (Symbol(s,PermRoot.rank(t)),)
    elseif s==:G return (Symbol(2,:I),6)
    else return (Symbol(orderphi,s,PermRoot.rank(t)),) 
    end
  else return (Symbol(s,PermRoot.rank(t)),) 
  end
  end
end

# functions whose result depends on cartanType
const needcartantype=Set([:Invariants,
                          :PrintDiagram,
                          :ReflectionName,
                          :UnipotentClasses,
                          :WeightInfo,
                          :CartanMat])

debug::Bool=false # time each call

"`chevieget(t::TypeIrred,f::Symbol,extra...)` get `CHEVIE[field(t)][f](extra...)`"
function chevieget(t::TypeIrred,f::Symbol,extra...)
  n,args...=field(t)
# println("d=$d f=$f extra=$extra")
  o=chevieget(n,f)
  if o isa Function
    if haskey(t,:orbit) t=t.orbit[1] end
if debug
    if haskey(t,:cartanType) && f in needcartantype
      println("$n.$f(",(args...,extra...,t.cartanType),")")
@time o(args...,extra...,t.cartanType)
    else 
     println("$n.$f(",(args...,extra...),")")
@time o(args...,extra...)
    end
else
    if haskey(t,:cartanType) && f in needcartantype
      o(args...,extra...,t.cartanType)
    else 
      o(args...,extra...)
    end
end
  elseif o===false return nothing
  else o
  end
end

end
