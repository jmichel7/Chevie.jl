"""
This  module contains  various utility  functions used  in the  rest of the
code.  Maybe some  of them  exist in  some Julia  module I am not aware of;
please tell me.

The code is divided in sections  according to semantics.
"""
module Util

using Gapjm

export getp, gets, # helpers for objects with a Dict of properties
  groupby, constant, # lists
  format, bracket_if_needed, ordinal, rshow, fromTeX, joindigits, # formatting
  factor, prime_residues, divisors, phi, primitiveroot, gcd_repr #number theory

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

#----------------------- Formatting -----------------------------------------
const  sup=Dict(zip("-0123456789+()abcdefghijklmnoprstuvwxyz",
                    "⁻⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁽⁾ᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻ"))
const  sub=Dict(zip("-0123456789,+()=aehijklmnoprstuvx",
                    "₋₀₁₂₃₄₅₆₇₈₉‚₊₍₎₌ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓ"))

"strip TeX formatting from  a string, using unicode characters to approximate"
function TeXstrip(s::String)
  s=replace(s,r"\$"=>"")
  s=replace(s,r"\\varepsilon"=>"ε")
  s=replace(s,r"\\delta"=>"δ")
  s=replace(s,r"\\gamma"=>"γ")
  s=replace(s,r"\\lambda"=>"λ")
  s=replace(s,r"\\phi"=>"φ")
  s=replace(s,r"\\Phi"=>"Φ")
  s=replace(s,r"\\psi"=>"ψ")
  s=replace(s,r"\\rho"=>"ρ")
  s=replace(s,r"\\sigma"=>"σ")
  s=replace(s,r"\\theta"=>"θ")
  s=replace(s,r"\\zeta"=>"ζ")
  s=replace(s,r"\\otimes"=>"⊗")
  s=replace(s,r"\\tilde A"=>"Ã")
  s=replace(s,r"\\times"=>"×")
  s=replace(s,r"\\BZ"=>"ℤ")
  s=replace(s,r"\\frakS"=>"𝔖")
  s=replace(s,r"\\wedge"=>"∧")
  s=replace(s,r"\\!"=>"")
  s=replace(s,r"{}"=>"")
  s=replace(s,r"_[-0-9,()+=aeh-pr-vx]"=>t->sub[t[2]])
  s=replace(s,r"(_\{[-0-9,()+=aeh-pr-vx]*\})('*)"=>s"\2\1")
  s=replace(s,r"_\{[-0-9,()+=aeh-pr-vx]*\}"=>t->map(x->sub[x],t[3:end-1]))
  s=replace(s,r"\^[-0-9,()a-op-z]"=>t->sup[t[2]])
  s=replace(s,r"\^\{[-0-9,()a-op-z]*\}"=>t->map(x->sup[x],t[3:end-1]))
  q(l)=l==1 ? "′" : l==2 ? "″" : l==3 ? "‴" : l==4 ? "⁗" : map(x->sup[x],"($l)")
  s=replace(s,r"''*"=>t->q(length(t)))
  s=replace(s,r"\{\+\}"=>"+")
  s
end

bracket_if_needed(c::String)=if occursin(r"[-+*/]",c[nextind(c,0,2):end]) 
 "($c)" else c end

function fromTeX(io::IO,n::String)
  TeX=get(io,:TeX,false) 
  if TeX return n end
  if get(io,:limit,false) return TeXstrip(n) end
  n=replace(n,r"\\tilde *"=>"~")
  n=replace(n,"_"=>"")
  n=replace(n,"}"=>"")
  n=replace(n,"{"=>"")
end

fromTeX(n::String;opt...)=fromTeX(IOContext(stdout,opt...),n)

"""
  format(io, table; options )

  General routine to format a table. Used for character tables.
  Options:
     row_labels          Labels for rows
     col_labels          Labels for columns
     rows_label          Label for column of rowLabels
     separators          line numbers after which to put a separator
     column_repartition  display in pieces of sizes these numbers of cols
     rows                show only these rows
     cols                show only these columns

"""
function format(io::IO,t::Matrix; opt...)
  io=IOContext(io,opt...)
  row_labels=get(io,:row_labels,axes(t,1))
  col_labels=get(io,:col_labels,nothing)
  rows_label=get(io,:rows_label,"")
  separators=get(io,:separators,[0])
  rows=get(io,:rows,axes(t,1))
  cols=get(io,:cols,axes(t,2))
  column_repartition=get(io,:column_repartition,nothing)
  lpad(s,n)=" "^(n-textwidth(s))*s # because lpad not what expected
  rpad(s,n)=s*" "^(n-textwidth(s)) # because rpad not what expected
  t=t[rows,cols]
  if eltype(t)!=String t=sprint.(show,t; context=io) end
  TeX=get(io,:TeX,false)
  row_labels=string.(row_labels[rows])
  colwidth=map(i->maximum(textwidth.(t[:,i])),axes(t,2))
  if !isnothing(col_labels)
    col_labels=string.(col_labels[cols])
    colwidth=map(max,colwidth,textwidth.(col_labels))
    if !TeX col_labels=map(lpad,col_labels,colwidth) end
  end
  labwidth=max(textwidth(rows_label),maximum(textwidth.(row_labels)))
  if !TeX
    rows_label=lpad(rows_label,labwidth)
    row_labels=rpad.(row_labels,labwidth)
  end
  function hline(ci)
    if TeX println(io,"\\hline")
    else
    print(io,"\u2500"^labwidth,"\u253C")
    print(io,"\u2500"^sum(colwidth[ci].+1),"\n")
    end
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
     if TeX column_repartition=[length(colwidth)]
     else column_repartition=cut(1 .+colwidth,displaysize(io)[2]-labwidth-1)
     end
  end
  ci=[0]
  for k in column_repartition
    ci=ci[end].+(1:k)
    if !isnothing(col_labels)
      if TeX
        println(io,"\\begin{array}{c|","c"^length(ci),"}")
        println(io,rows_label,"&",join(col_labels[ci],"&"),"\\\\")
      else println(io,rows_label,"\u2502",join(col_labels[ci]," "))
      end
      if 0 in separators hline(ci) end
    end
    for l in axes(t,1)
      if TeX
        println(io,row_labels[l],"&",join(t[l,ci],"&"),"\\\\")
      else
        println(io,row_labels[l],"\u2502",join(map(lpad,t[l,ci],colwidth[ci])," "))
      end
      if l in separators hline(ci) end
    end
    if ci[end]!=length(colwidth) print(io,"\n") end
    if TeX println(io,"\\end{array}") end
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

# show with attributes...
rshow(x;p...)=show(IOContext(stdout,:limit=>true,p...),"text/plain",x)

joindigits(l::AbstractVector,sep="()")=any(l.>=10) ? 
                 string(sep[1:1],join(l,","),sep[2:2]) : join(l)

#----------------------- Number theory ---------------------------
" the numbers less than n and prime to n "
function prime_residues(n)
  if n==1 return [0] end
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
  if n==1 || x==1 return x
  elseif n==2
    res=ER(x)
    if HasType.CHEVIE[:info]  println("root($x,$n)=$res") end
    return res
  elseif x==-1 && n%2==1 return x
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
end
