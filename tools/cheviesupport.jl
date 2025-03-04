"""
The  module GAPENV creates a GAP3-like environment by extending locally the
base  functions `*, +, -, ^, isless, copy, //, inv, length` to their
GAP3  semantics and defining  quite a few  other GAP3 functions, then loads
the Chevie database in that environment.
"""
module GAPENV
include("gap3support.jl")
# --------- pirating extensions to get closer to GAP semantics -------------
# The  idea how to implement base function pirating restricted to the current
# module is due to Max Horn.
*(a...)=Base.:*(a...)
*(a::AbstractVector{<:Number},b::AbstractVector{<:Number})=transpose(a)*b
*(a::AbstractVector,b::AbstractVector{<:AbstractVector})=toL(toM(a)*toM(b))
*(a::AbstractVector{<:Number},b::AbstractVector{<:AbstractVector})=toL(transpose(a)*toM(b))[1]
*(a::AbstractVector,b::AbstractVector)=sum(a.*b)
*(a::Tuple,b::AbstractVector)=toL(transpose(collect(a))*toM(b))[1]
*(a::AbstractVector{<:AbstractVector{<:Number}},b::Number)=map(x->map(y->y*b,x),a)

+(a...)=Base.:+(a...)
+(a::AbstractArray,b::Number)=a .+ b
+(a::Integer,b::AbstractVector)=a .+ b
+(a::Integer,b::Integer,c::AbstractVector)=(a+b).+c

-(a...)=Base.:-(a...)
-(a::AbstractVector,b::Number)=a .- b
-(a::Integer,b::AbstractVector)=a .- b

^(a...)=Base.:^(a...)
^(m::AbstractMatrix,n::AbstractMatrix)=inv(n*1//1)*m*n
^(m::AbstractVector{<:AbstractVector{<:Number}},n::Matrix{<:Number})=inv(n*1//1)*toM(m)*n
^(m::AbstractVector{<:AbstractVector},n::Integer)=toL(toM(m)^n)
^(m::AbstractVector,n::AbstractVector)=toL(inv(toM(n)*1//1)*toM(m)*toM(n))

isless(a...)=Base.isless(a...)
isless(a::Array,b::Number)=true # to sort [[1,2], 2, 0]
isless(b::Number,a::Array)=false
function isless(A::AbstractVector, B::AbstractVector)
  for (a, b) in zip(A, B) if !isequal(a, b) return isless(a, b) end end
  isless(length(A), length(B))
end

copy(a...)=Base.copy(a...)
copy(x::Char)=x

(//)(a...)=Base.:(//)(a...)
(//)(m::Vector,n::Vector)=toL(toM(m)*inv(toM(n)*E(1)))

inv(a...)=Base.inv(a...)
inv(m::Vector)=toL(inv(toM(m)*E(1)//1))

length(a...)=Base.length(a...)
length(a::Symbol)=length(string(a))

# for some reason getindex is special, reverting to plain pirating
Base.getindex(s::String,a::Vector{Any})=getindex(s,Int.(a))
Base.getindex(a::Symbol,i::Int)=string(a)[i]

Base.push!(a::String,b::Char)=a*b
# ---------------- other extensions --------------------------------------
Base.:*(a::AbstractArray,b::Union{Pol,Mvp})=isone(b) ? a : a .* b
Base.:*(a::AbstractArray,b::Frac)=isone(b) ? a : a .* b
Base.:*(a::Union{Pol,Mvp},b::AbstractArray)=isone(a) ? b : a .* b
Base.:*(a::Frac,b::AbstractArray)=isone(a) ? b : a .* b
Base.:*(W1::Spets,W2::FiniteCoxeterGroup)=Cosets.extprod(W1,spets(W2))
Base.:*(W1::FiniteCoxeterGroup,W2::Spets)=Cosets.extprod(spets(W1),W2)
Base.:+(a::AbstractArray,b::Union{Pol,Mvp})=a .+ b
Base.:/(a::AbstractArray,b::Union{Pol,Mvp})=a ./ b
Base.://(a::AbstractArray,b::Union{Pol,Mvp})=a .// b
Base.:^(a::Cyc,b::Rational)=a^Int(b)
Base.:^(f::Family,n::Int)=galois(f,n)
FinitePosets.Poset(m::Vector{Vector{Bool}})=Poset(CPoset(toM(m)))
PermRoot.roots(C::Vector{<:Vector})=roots(toM(C))

Product(a::AbstractVector{<:AbstractVector{<:AbstractVector}})=
toL(prod(toM.(a)))
# translations of GAP3 functions for the Chevie library
CartanMat(s,a...)=cartan(Symbol(s),a...)
CharParams(W)=charinfo(W).charparams
Concatenation(a::Vector,b::Tuple)=vcat(a,improve_type(collect(b)))
CoxeterGroup()=coxgroup()
function CycleStructurePerm(p::Perm)
  if isone(p) return Int[] end
  t=cycletype(p)
  res=[0,nothing]
  resize!(res,t[1]-1)
  res.=nothing
  for (k,v) in tally(t) res[k-1]=v end
  res
end
function Csp2cycletype(c)
  res=Pair{Int,Int}[]
  for (i,m) in enumerate(c)
    if !isnothing(m) push!(res,i+1=>m) end
  end
  res
end
function DiagonalMat(v...)
  arg=map(v)do m
    if m isa Array 0 .+toM(m)
    else 0 .+hcat(m)
    end
  end
  R=cat(arg...;dims=(1,2))
  u=reshape(R,(prod(size(R),)))
  for i in eachindex(u) if !isassigned(u,i) u[i]=zero(u[1]) end end
  toL(R)
end
DiagonalMat(v::Vector)=DiagonalMat(v...)
function ExteriorPower(m,i)
  if m isa Matrix exterior_power(m,i)
  else toL(exterior_power(toM(m),i))
  end
end
Inherit(a,b)=merge!(a,b)
function Inherit(a,b,c)
  for k in c a[Symbol(k)]=b[Symbol(k)] end
  a
end
KroneckerProduct(a,b)=toL(kron(toM(a),toM(b)))
const PartitionTupleToString(s,d=Dict())=string_partition_tuple(s;d...)
ReflectionSubgroup(W,I::AbstractVector)=reflection_subgroup(W,convert(Vector{Int},I))
SchurFunctor(m,p)=toL(schur_functor(toM(m),p))
SymmetricPower(m,n)=SchurFunctor(m,[n])
StringToDigits(s)=collect(s).-'0'
Weyl.torus(m::Vector{<:Vector})=torus(toM(m))
function CoxeterGroup(S::String,s...)
 res=coxgroup(Symbol(S),Int(s[1])) 
 for i in 2:2:length(s) res*=coxgroup(Symbol(s[i]),Int(s[i+1])) end
 res
end

#-------------------------------------------------------------------------
#  dummy translations of GAP3 functions
Format(x)=string(x)
FormatTeX(x)=xrepr(x,TeX=true)
Format(x,opt)=xrepr(x;opt...)

GetRoot(x::Cyc,n::Number=2,msg...)=root(x,n)
GetRoot(x::Integer,n::Number=2,msg...)=root(x,n)
GetRoot(x::Pol,n::Number=2,msg...)=root(x,n)
GetRoot(x::Rational,n::Number=2,msg...)=root(x,n)
GetRoot(x::Mvp,n::Number=2,msg...)=root(x,n)
GetRoot(x::Root1,n::Number=2,msg...)=root(x,n)
function GetRoot(x,n::Number=2,msg...)
  error("GetRoot($x,$n) not implemented")
end

Unbind(x)=x
#-------------------------------------------------------------------------
Cosets.spets(W::FiniteCoxeterGroup,v::Vector)=spets(W,toM(v))
#-----------------------------------------------------------------------
#if false
#"""
#An `Unknown()` represents an element of a ring which is not known. The main
#difference  with  `missing`  is  that  `0*Unknown()==0`  (that  at least we
#know!). An unknown is printed at the repl as `?`.
#
#```julia-rep1
#julia> a=[GAPENV.Unknown(),GAPENV.Unknown()]
#2-element Vector{Chevie.GAPENV.Unknown}:
# ?
# ?
#
#julia> a.*[0,1]
#2-element Vector{Any}:
# 0
#  ?
#```
#"""
#struct Unknown end
#
#Base.:+(a,b::Unknown)=b
#Base.:+(b::Unknown,a)=b
#Base.:+(b::Unknown,a::Unknown)=b
#Base.:*(a,b::Unknown)=iszero(a) ? a : b
#Base.:*(b::Unknown,a)=iszero(a) ? a : b
#Base.:*(b::Unknown,a::Unknown)=b
#Base.zero(a::Unknown)=0
#Base.show(io::IO,a::Unknown)=print(io,get(io,:limit,false) ? "?" : "Unknown()")
#Base.isless(a::Unknown,b::Number)=false
#Base.isless(b::Number,a::Unknown)=true
#Base.isless(b::Unknown,a::Unknown)=false
#Base.:(//)(a::Unknown,b)=a
#Base.isreal(a::Unknown)=false
#Base.isinteger(a::Unknown)=false
#Base.conj(a::Unknown)=a
#Base.broadcastable(a::Unknown)=Ref(a)
#CyclotomicNumbers.galois(a::Unknown,i)=a
#else
Unknown()=missing
Base.:*(a::Missing,b)=a
#end
end
