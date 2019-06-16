"""
This  module contains  various utility  functions used  in the  rest of the
code.  Maybe some  of them  exist in  some Julia  module I am not aware of;
please tell me.

The code is divided in sections  according to semantics.
"""
module Util

using Gapjm

export getp, gets, # helpers for objects with a Dict of properties
  groupby, constant, blocks, # arrays
  ModuleElt, norm!, # data structure
  format, TeXstrip, bracket_if_needed, ordinal, # formatting
  factor, prime_residues, divisors, phi, primitiveroot, gcd_repr, #number theory
  conjugate_partition, horner, partitions, combinations, arrangements,
  partition_tuples, #combinatorics
  echelon, exterior_power  # linear algebra

# not exported: nullspace, to avoid conflict with LinearAlgebra

#--------------------------------------------------------------------------
"""
  a variant of get! for objects O which have a Dict of properties named prop.
  Usually called as
    gets(O,:p) do ---code to compute property :p --- end
"""
gets(f::Function,o,p::Symbol)=get!(()->f(o),o.prop,p)

"""
  A  variation where it is assumed that f sets key p but not assumed that f
  returns  the value  of property  p, because  f could  set several keys at
  once...
"""
function getp(f::Function,o,p::Symbol)
  if haskey(o.prop,p) return o.prop[p] end
  f(o)
  o.prop[p]
end
#--------------------------------------------------------------------------
"""
  group items of list l according to the corresponding values in list v

    julia> groupby([31,28,31,30,31,30,31,31,30,31,30,31],
           [:Jan,:Feb,:Mar,:Apr,:May,:Jun,:Jul,:Aug,:Sep,:Oct,:Nov,:Dec])
    Dict{Int64,Array{Symbol,1}} with 3 entries:
      31 => Symbol[:Jan, :Mar, :May, :Jul, :Aug, :Oct, :Dec]
      28 => Symbol[:Feb]
      30 => Symbol[:Apr, :Jun, :Sep, :Nov]

"""
function groupby(v::AbstractVector,l::AbstractVector)
  res=Dict{eltype(v),Vector{eltype(l)}}()
  for (k,val) in zip(v,l)
    push!(get!(res,k,empty(l)),val)
  end
  res
end

"""
  group items of list l according to the values taken by function f on them

    julia> groupby(iseven,1:10)
    Dict{Bool,Array{Int64,1}} with 2 entries:
      false => [1, 3, 5, 7, 9]
      true  => [2, 4, 6, 8, 10]

Note:in this version l is required to be non-empty since I do not know how to
access the return type of a function
"""
function groupby(f,l::AbstractVector)
  res=Dict(f(l[1])=>[l[1]]) # l should be nonempty
  for val in l[2:end]
    push!(get!(res,f(val),empty(l)),val)
  end
  res
end

" whether all elements in list a are equal"
function constant(a::AbstractVector)
   all(i->a[i]==a[1],2:length(a))
end

"""
  blocks(M::Matrix)

  M  should be a square matrix. Define  a graph G with vertices 1:size(M,1)
  and  with an edge between i and j  if either M[i,j] or M[j,i] is not zero
  or false. blocks returns a vector of vectors I such that I[1],I[2], etc..
  are  the  vertices  in  each  connected  component  of G. In other words,
  M[I[1],I[1]],M[I[2],I[2]],etc... are blocks of M.
"""
function blocks(M::Matrix)::Vector{Vector{Int}}
  l=size(M,1)
  if l==0 return Vector{Int}[] end
  cc=collect(1:l) # cc[i]: in which block is i, initialized to different blocks
  nz=x->!iszero(x)
  for i in 1:l, j in i+1:l
    # if new relation i~j then merge components:
    if (nz(M[i,j]) || nz(M[j,i])) && cc[i]!=cc[j]
      cj=cc[j]
      for k in 1:l
         if cc[k]==cj cc[k]=cc[i] end
      end
    end
  end
  sort(collect(values(groupby(cc,collect(1:l)))))
end
#------------- Make a version with Dict of ModuleElt ----------------------
const usedict=false
if usedict
struct ModuleElt{K,V}
  d::Dict{K,V}
end

ModuleElt(x::Vector{Pair{K,V}}) where{K,V}=ModuleElt(Dict(x))
ModuleElt(x::Pair{K,V}...) where{K,V}=ModuleElt(Dict(x...))
ModuleElt(x::Base.Generator)=ModuleElt(Dict(x))

Base.zero(::Type{ModuleElt{K,V}}) where{K,V}=ModuleElt(Dict{K,V}())
Base.iszero(x::ModuleElt)=isempty(x.d)
Base.zero(x::ModuleElt)=ModuleElt(empty(x.d))
@inline Base.iterate(x::ModuleElt,y...)=iterate(x.d,y...)
@inline Base.length(x::ModuleElt)=length(x.d)
@inline Base.push!(x::ModuleElt,y...)=push!(x.d,y...)
@inline Base.haskey(x::ModuleElt,y...)=haskey(x.d,y...)
@inline Base.append!(x::ModuleElt,y...)=append!(x.d,y...)
@inline Base.getindex(x::ModuleElt,i)=getindex(x.d,i)
@inline Base.cmp(x::ModuleElt,y::ModuleElt)=x.d == y.d ? 0 : 1

Base.:-(a::ModuleElt)=ModuleElt(k=>-v for (k,v) in a)

# multiply module element by scalar
Base.:*(a::ModuleElt,b)=iszero(b) ? zero(a) : ModuleElt(k=>v*b for (k,v) in a)

function norm!(a::ModuleElt)
  for (k,v) in a.d if iszero(v) delete!(a.d,k) end end
  a
end

Base.:+(a::ModuleElt,b::ModuleElt)::ModuleElt=norm!(ModuleElt(merge(+,a.d,b.d)))

else
#--------------------------------------------------------------------------
"""
ModuleElt  represents  an  element  of  a  module  with basis of type K and
coefficients  of type  V. It  has a  similar interface  to Dict{K,V}, but +
below is 3 times faster than merge(+,...) on Dicts.
"""
struct ModuleElt{K,V}
  d::Vector{Pair{K,V}} # This vector is kept sorted by K
end

ModuleElt(x::Pair{K,V}...) where{K,V}=ModuleElt(collect(x))
ModuleElt(x::Base.Generator)=ModuleElt(collect(x))

Base.zero(::Type{ModuleElt{K,V}}) where{K,V}=ModuleElt(Pair{K,V}[])
Base.iszero(x::ModuleElt)=isempty(x.d)
Base.zero(x::ModuleElt)=ModuleElt(empty(x.d))
@inline Base.iterate(x::ModuleElt,y...)=iterate(x.d,y...)
@inline Base.length(x::ModuleElt)=length(x.d)
@inline function Base.push!(x::ModuleElt,y...)
  push!(x.d,y...)
  x
end
@inline function Base.append!(x::ModuleElt,y...)
  append!(x.d,y...)
  x
end
@inline Base.cmp(x::ModuleElt,y::ModuleElt)=cmp(x.d,y.d)

Base.:-(a::ModuleElt)=ModuleElt(k=>-v for (k,v) in a)

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
normalize an unsorted ModuleElt -- which may occur in some computation
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

function drop(m::ModuleElt,k)
  r=searchsorted(m.d,Ref(k);by=first)
  if r.start!=r.stop return nothing end
  ModuleElt(deleteat!(copy(m.d),r.start)),m.d[r.start][2]
end

Base.keys(x::ModuleElt)=first.(x.d)
Base.first(x::ModuleElt)=first(x.d)
#--------------------------------------------------------------------------
end
"strip TeX formatting from  a string, using unicode characters to approximate"
function TeXstrip(s::String)
  s=replace(s,r"\\varepsilon"=>"ε")
  s=replace(s,r"\\gamma"=>"γ")
  s=replace(s,r"\\phi"=>"ϕ")
  s=replace(s,r"\\Phi"=>"Φ")
  s=replace(s,r"\\rho"=>"ρ")
  s=replace(s,r"\\zeta"=>"ζ")
  s=replace(s,r"\\otimes"=>"⊗")
  s=replace(s,r"\\tilde A"=>"Ã")
  s=replace(s,r"\\times"=>"×")
  s=replace(s,r"\\BZ"=>"ℤ")
  s=replace(s,r"\\wedge"=>"∧")
  s=replace(s,r"\\!"=>"")
  s=replace(s,r"{}"=>"")
  sup=Dict(zip("0123456789-()","⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁽⁾"))
  sub=Dict(zip("-0123456789,()+=aehijklmnoprstuvx",
               "₋₀₁₂₃₄₅₆₇₈₉‚₍₎₊₌ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓ"))
  s=replace(s,r"_[-0-9,()+=aehijklmnoprstuvx]"=>t->sub[t[2]])
  s=replace(s,r"(_\{[0-9,]*\})('*)"=>s"\2\1")
  s=replace(s,r"_\{[-0-9,()+=aehijklmnoprstuvx]*\}"=>t->map(x->sub[x],t[3:end-1]))
  s=replace(s,r"\^[-0-9()]"=>t->sup[t[2]])
  s=replace(s,r"\^\{[-0-9()]*\}"=>t->map(x->sup[x],t[3:end-1]))
  q(l)= l==1 ? "′" : l==2 ? "″" : l==3 ? "‴" : l==4 ? "⁗" : 
     "⁽$(map(x->sup[x],string(l)))⁾"
  s=replace(s,r"''*"=>t->q(length(t)))
  s
end

bracket_if_needed(c::String)=if occursin(r"[-+*/]",c[nextind(c,0,2):end]) 
 "($c)" else c end

"""
  format( table; options )

  General routine to format a table. Used for character tables.
  Options:
     row_labels          Labels for rows
     column_labels       Labels for columns
     rows_label          Label for column of rowLabels
     separators          line numbers after which to put a separator
     column_repartition  display in pieces of sizes these numbers of cols
     rows                show only these rows
     columns             show only these columns

"""
function format(io::IO,t::Matrix; row_labels=axes(t,1),
  column_labels=nothing, rows_label="", separators=[0], rows=axes(t,1),
  columns=axes(t,2), column_repartition=nothing)
  t=t[rows,columns]
  if eltype(t)!=String t=sprint.(show,t; context=io) end
  row_labels=string.(row_labels[rows])
  colwidth=map(i->maximum(textwidth.(t[:,i])),axes(t,2))
  if !isnothing(column_labels)
    column_labels=string.(column_labels[columns])
    colwidth=map(max,colwidth,textwidth.(column_labels))
    column_labels=map(lpad,column_labels,colwidth)
  end
  labwidth=max(textwidth(rows_label),maximum(textwidth.(row_labels)))
  rows_label=lpad(rows_label,labwidth)
  row_labels=rpad.(row_labels,labwidth)
  function hline(ci)
    print(io,"\u2500"^labwidth,"\u253C")
    print(io,"\u2500"^sum(colwidth[ci].+1),"\n")
  end
  function cut(l,max) # cut Integer list l in parts of sum<max
    res=Int[];len=0;n=0
    for i in l len+=i
      if len>=max
        if n==0 push!(res,1);n=0;len=0
        else push!(res,n);n=1;len=i
        end
      else n+=1
      end
    end
    push!(res,n)
  end
  if isnothing(column_repartition)
     column_repartition=cut(1 .+colwidth,displaysize(io)[2]-labwidth-1)
  end
  ci=[0]
  for k in column_repartition
    ci=ci[end].+(1:k)
    if !isnothing(column_labels)
      print(io,rows_label,"\u2502",join(column_labels[ci]," "),"\n")
      if 0 in separators hline(ci) end
    end
    for l in axes(t,1)
      print(io,row_labels[l],"\u2502",join(map(lpad,t[l,ci],colwidth[ci])," "),"\n")
      if l in separators hline(ci) end
    end
    if ci[end]!=length(colwidth) print(io,"\n") end
  end
end

function ordinal(n)
  str=repr(n)
  if     n%10==1 && n%100!=11 str*="st"
  elseif n%10==2 && n%100!=12 str*="nd"
  elseif n%10==3 && n%100!=13 str*="rd"
  else                        str*="th"
  end
  str
end

#----------------------- Number theory ---------------------------
" the numbers less than n and prime to n "
function prime_residues(n)
  filter(i->gcd(n,i)==1,1:n-1)
end

# make Primes.factor fast by memoizing it
import Primes
const dict_factor=Dict(2=>Primes.factor(2))
function factor(n::Integer)
  get!(dict_factor,n) do 
    Primes.factor(Dict,n) 
  end
end

function divisors(n::Int)::Vector{Int}
  if n==1 return [1] end
  sort(vec(map(prod,Iterators.product((p.^(0:m) for (p,m) in factor(n))...))))
end

" the Euler function ϕ "
function phi(m::T)where T<:Integer
  if m==1 return 1 end
  prod(p->p[1]^(p[2]-1)*(p[1]-1),factor(m))
end

"""
  primitiveroot(m::Integer) a primitive root mod. m,
  that is it generates multiplicatively prime_residues(m).
  It exists if m is of the form 4, 2p^a or p^a for p prime>2.
"""
function primitiveroot(m::Integer)
 if m==2 return 1
 elseif m==4 return 3
 end
 f=factor(m)
 nf=length(keys(f))
 if nf>2 return nothing end
 if nf>1 && (!(2 in keys(f)) || f[2]>1) return nothing end
 if nf==1 && (2 in keys(f)) && f[2]>2 return nothing end
 p=phi(m)
 1+findfirst(x->powermod(x,p,m)==1 && 
             all(d->powermod(x,d,m)!=1,divisors(p)[2:end-1]),2:m-1)
end

"""
  gcd_repr(x,y) returns a,b such that ax+by=gcd(x,y)
"""
function gcd_repr(x,y)
  (f,fx)=(x,1)
  (g,gx)=(y,0)
  while !iszero(g)
    (h,hx)=(g,gx)
    q,r=divrem(f,g)
    (g,gx)=(r,fx-q*gx)
    (f,fx)=(h,hx)
  end
  q=sign(f)
  (q*fx, iszero(y) ? 0 : div(q * (f - fx * x), y ))
end

function Gapjm.root(x::Integer,n::Number=2)
  if n==1 return x
  elseif n==2
    res=ER(x)
    println("root($x,$n)=$res")
    return res
  else
    error("root($x,$n) not implemented")
  end
end

Gapjm.root(x::Rational{<:Integer},n::Number=2)=root(numerator(x),n)//root(denominator(x),n)
#--------------------------------------------------------------------------
# written since should allow negative powers with inv
#function Base.:^(a::T, p::Integer) where T
#    if p ≥ 0 Base.power_by_squaring(a, p)
#    else     Base.power_by_squaring(inv(a)::T, -p)
#    end
#end

# better display of Rationals at the REPL
#function Base.show(io::IO, x::Rational)
#   show(io, numerator(x))
#   if get(io, :limit, true)
#       if denominator(x)!=1
#          print(io, "/")
#          show(io, denominator(x))
#       end
#   else
#       print(io, "//")
#       show(io, denominator(x))
#   end
#end

function conjugate_partition(p)
  res=zeros(eltype(p),maximum(p))
  for i in p, j in 1:i res[j]+=1 end
  res
end

# horner scheme
function horner(x,p::Vector)
  value=zero(x)
  for i in length(p):-1:1
    value=x*value+p[i]
  end
  value
end

# partitions of n of first (greatest) part <=m
function partitions_less(n,m)
  if m==1 return [fill(1,n)] end
  if iszero(n) return [Int[]] end
  res=Vector{Int}[]
  for i in 1:min(m,n)
    append!(res,map(x->vcat([i],x),partitions_less(n-i,i)))
  end
  res
end

partitions(n)=partitions_less(n,n)

if false
function partition_tuples(n,r)::Vector{Vector{Vector{Int}}}
  if r==1 return iszero(n) ? [[Int[]]] : map(x->[x],partitions(n)) end
  res=Vector{Vector{Int}}[]
  for i in  n:-1:1
    for p1 in partitions(i), p2 in partition_tuples(n-i,r-1)
      push!(res,vcat([p1],p2))
    end 
  end
  for p2 in partition_tuples(n,r-1)
    push!(res,vcat([Int[]],p2))
  end 
  res
end
else # bas implementation but which has same order as GAP3
function partition_tuples(n, r)
   if n==0 return [fill(Int[],r)] end
   empty=(tup=[Int[] for i in 1:r], pos=fill(1,n-1))
   pm=[typeof(empty)[] for i in 1:n-1]
   for m in 1:div(n,2)
      for i in 1:r
         t1=map(copy,empty)
         t1.tup[i]=[m]
         t1.pos[m]=i
         push!(pm[m],t1)
      end
      for k in m+1:n-m
         for t in pm[k-m]
            for i in t.pos[m]:r
               t1=map(copy,t)
               t1.tup[i]=vcat([m],t1.tup[i])
               t1.pos[m]= i
               push!(pm[k], t1)
            end
         end
      end
   end
   res= Vector{Vector{Int}}[]
   for k in 1:n-1
      for t in pm[n-k]
         for i in t.pos[k]:r
            s=copy(t.tup)
            s[i]=vcat([k],s[i])
            push!(res,s)
         end
      end
   end
   for i in 1:r
      s=copy(empty.tup)
      s[i]=[n]
      push!(res,s)
   end
   res
end
end

function combinations_sorted(mset::AbstractVector,k)
  if iszero(k) return [eltype(mset)[]] end
  res=Vector{eltype(mset)}[]
  for (i,e) in enumerate(mset)
    append!(res,map(x->vcat([e],x),combinations_sorted(mset[i+1:end],k-1)))
  end
  res
end 

combinations(mset,k)=combinations_sorted(sort(mset),k)

function ArrangementsK(mset,blist,k)
  if iszero(k) return [eltype(mset)[]] end
  combs=Vector{eltype(mset)}[]
  for i in eachindex(mset)
    if blist[i] && (i==length(mset) || mset[i+1]!=mset[i] || !blist[i+1])
      blist1=copy(blist)
      blist1[i]=false
      append!(combs,pushfirst!.(ArrangementsK(mset,blist1,k-1),Ref(mset[i])))
    end
  end
  combs
end 

arrangements(mset,k)=ArrangementsK(sort(mset),fill(true,length(mset)),k)


#----------- Linear algebra over Rationals/integers------------------------
" returns: echelon form of m, indices of linearly independent rows of m"
function echelon!(m::Matrix)
  T=typeof(one(eltype(m))//1)
  if T!=eltype(m) m=convert.(T,m) end
  rk=0
  inds=collect(axes(m,1))
  for k in axes(m,2)
    j=findfirst(x->!iszero(x),m[rk+1:end,k])
    if isnothing(j) continue end
    j+=rk
    rk+=1
    row=m[j,:]
    m[j,:].=m[rk,:]
    m[rk,:].=inv(row[k]).*row
    inds[[j,rk]]=inds[[rk,j]]
    for j in axes(m,1)
      if rk!=j && !iszero(m[j,k]) m[j,:].-=m[j,k].*m[rk,:] end
    end
#   println(m)
  end
  m,inds[1:rk]
end

echelon(m::Matrix)=echelon!(copy(m))

" computes right nullspace of m in a type_preserving way"
function nullspace(m::Matrix)
  m=echelon(m)[1]
  n=size(m)[2]
  z=Int[]
  j=0
  lim=size(m,1)
  for i in axes(m,1)
    f=findfirst(x->!iszero(x),m[i,j+1:end])
    if isnothing(f)
      lim=i-1
      break 
    end
    j+=f
    push!(z,j)
  end
# println("z=$z lim=$lim")
  zz=zeros(eltype(m),n,n)
  zz[z,:]=m[1:lim,:]
  nn=filter(k->iszero(zz[k,k]),1:n)
  for i in nn zz[i,i]=-one(eltype(m)) end
  zz[:,nn]
end

function det(A)
  if size(A,1)==1 return A[1,1]
  elseif size(A,1)==2 return A[1,1]*A[2,2]-A[1,2]*A[2,1]
  end
  sum(i->det(A[vcat(1:i-1,i+1:size(A,1)),2:end])*(-1)^(i-1),axes(A,1))
end

function exterior_power(A,m)
  basis=combinations(1:size(A,1),m)
  [det(A[i,j]) for i in basis, j in basis]
end 

end
