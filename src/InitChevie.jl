# mostly translations of chevie/init.g and chevie/tbl/compat3.g
module InitChevie

using ..Chevie
export CHEVIE, chevieget, chevieset, getchev

const CHEVIE=Dict{Symbol,Any}(
  :compat=>Dict{Symbol,Any}(:MakeCharacterTable=>x->x),
  :info=>false
)

CHEVIE[:compat][:AdjustHeckeCharTable]=function(tbl,param)
  r=eachindex(tbl[:classtext])
  param=improve_type(param)
  for i in r 
    f=prod(-last.(param[tbl[:classtext][i]]))
    for j in r 
      if tbl[:irreducibles] isa Matrix
        tbl[:irreducibles][j,i]*=f
      else tbl[:irreducibles][j][i]*=f
      end
    end
  end
end

" chevieget(t,w) returns CHEVIE[Symbol(t)][w]"
function chevieget(t::Symbol,w::Symbol)
# println("chevieget(",t,",",w,")")
  get!(CHEVIE[t],w)do
    if CHEVIE[:info] println("CHEVIE[$t] has no $w") end
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
  if s in [:A,:B,:D] 
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
  elseif s==:ST 
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

"`getchev(t::TypeIrred,f::Symbol,extra...)` get `CHEVIE[field(t)][f](extra...)`"
function getchev(t::TypeIrred,f::Symbol,extra...)
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

function getchev(W,f::Symbol,extra...)
  [getchev(ti,f::Symbol,extra...) for ti in refltype(W)]
end

end
