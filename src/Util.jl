"""
This  module contains  various utility  functions used  in the  rest of the
code.  Maybe some  of them  exist in  some Julia  module I am not aware of;
please tell me.

The code is divided in sections  according to semantics.
"""
module Util

export getp, gets, # helpers for objects with a Dict of properties
  groupby, constant, blocks, # arrays
  SortedPairs, norm!, mergesum, getvalue, # data structure
  format, TeXstrip, bracket_if_needed, ordinal, # formatting
  factor, prime_residues, divisors, phi, primitiveroot,  # number theory
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
#--------------------------------------------------------------------------
"""
SortedPairs has a similar interface to Dicts, but is 3 times faster for the
merge operation.  Pairs are sorted by the first item (the key)
"""
const SortedPairs{K,V}=Vector{Pair{K,V}} where {K,V}

"""
merge is like merge(+,..) for Dicts, with the difference that keys with
value 0 are deleted
"""
function mergesum(a::SortedPairs,b::SortedPairs)::SortedPairs
  la=length(a)
  lb=length(b)
  res=similar(a,la+lb)
  ai=bi=1
  ri=0
  while ai<=la || bi<=lb
    if ai>la
      res[ri+=1]=b[bi]; bi+=1
    elseif bi>lb
      res[ri+=1]=a[ai]; ai+=1
    else
      c=cmp(a[ai][1],b[bi][1])
      if c==1
        res[ri+=1]=b[bi]; bi+=1
      elseif c==-1
        res[ri+=1]=a[ai]; ai+=1
      else s=a[ai][2]+b[bi][2]
        if !iszero(s)
          res[ri+=1]=a[ai][1]=>s
        end
        ai+=1; bi+=1
      end
    end
  end
  resize!(res,ri)
end

"""
normalize an unsorted SortedPairs -- which may occur in some computation
"""
function norm!(x::SortedPairs{K,V}) where {K,V}
  if isempty(x) return x end
  sort!(x,by=first)
  ri=1
  for j in 2:length(x)
    if x[j][1]==x[ri][1]
      x[ri]=x[ri][1]=>x[ri][2]+x[j][2]
    else 
      if !iszero(x[ri][2]) ri+=1 end
      x[ri]=x[j]
    end
  end
  if iszero(x[ri][2]) ri-=1 end
  resize!(x,ri)
end

function getvalue(x::SortedPairs,i)
  r=searchsorted(x,i;by=first)
  if r.start!=r.stop error("Bounds in $x") end
  x[r.start][2]
end
#--------------------------------------------------------------------------
"strip TeX formatting from  a string, using unicode characters to approximate"
function TeXstrip(s::String)
  s=replace(s,r"\\phi"=>"ϕ")
  s=replace(s,r"\\Phi"=>"Φ")
  s=replace(s,r"\\zeta"=>"ζ")
  s=replace(s,r"\\tilde A"=>"Ã")
  s=replace(s,r"\\times"=>"×")
  sup=Dict(zip("0123456789-","⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))
  sub=Dict(zip("0123456789,()","₀₁₂₃₄₅₆₇₈₉‚₍₎"))
  s=replace(s,r"_[0-9]"=>t->sub[t[2]])
  s=replace(s,r"\^[0-9]"=>t->sup[t[2]])
  s=replace(s,r"(_\{[0-9,]*\})('*)"=>s"\2\1")
  s=replace(s,r"_\{[0-9,()]*\}"=>t->map(x->sub[x],t[3:end-1]))
  s=replace(s,r"\^\{[-0-9]*\}"=>t->map(x->sup[x],t[3:end-1]))
  s=replace(s,r"'''''"=>"‴″")
  s=replace(s,r"''''"=>"⁗")
  s=replace(s,r"'''"=>"‴")
  s=replace(s,r"''"=>"″")
  s=replace(s,r"'"=>"′")
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

function combinations_sorted(mset::AbstractVector,k)
  if iszero(k) return [eltype(mset)[]] end
  res=Vector{eltype(mset)}[]
  for (i,e) in enumerate(mset)
    append!(res,map(x->vcat([e],x),combinations_sorted(mset[i+1:end],k-1)))
  end
  res
end 

combinations(mset,k)=combinations_sorted(sort(mset),k)

function ArrangementsK(mset::AbstractVector,blist,k)
  if iszero(k) return [eltype(mset)[]] end
  combs=Vector{eltype(mset)}[]
  for i in eachindex(mset)
    if blist[i]
    blist1=copy(blist)
    blist1[i]=false
    append!(combs,map(x->vcat([mset[i]],x),ArrangementsK(mset,blist1,k-1)))
    end
  end
  unique(combs)
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
  println("z=$z lim=$lim")
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
