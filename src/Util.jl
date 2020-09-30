"""
This  module contains  various utility  functions used  in the  rest of the
code.  Maybe some  of them  exist in  some Julia  module I am not aware of;
please tell me.

The code is divided in sections  according to semantics.
"""
module Util

export 
  @forward,
  getp, gets, # helpers for objects with a Dict of properties
  format, format_coefficient, ordinal, rshow, rsprint, printc, fromTeX, 
  printTeX,joindigits, cut, # formatting
  factor, prime_residues, divisors, phi, primitiveroot #number theory

export toL, toM # convert Gap matrices <-> Julia matrices
export ds # dump struct
export InfoChevie

const info=Ref(true)
function InfoChevie(a...)
  if Util.info[] printc(a...) end
end

"""
@forward T.f f1,f2,...
  generates 
f1(x::T,y...)=f1(x.f,y...)
f2(x::T,y...)=f2(x.f,y...)
...
"""
macro forward(ex, fs)
  T, field = esc(ex.args[1]), ex.args[2].value
  fdefs=map(fs.args)do ff
      f= esc(ff)
      quote
        @inline ($f)(a::($T),args...)=($f)(a.$field,args...)
      end
  end
  Expr(:block, fdefs...)
end

#--------------------------------------------------------------------------
toL(m)=collect(eachrow(m)) # to Gap
toM(m)=isempty(m) ? permutedims(hcat(m...)) : permutedims(reduce(hcat,m)) # to julia

printc(xs...;p...)=print(IOContext(stdout,:limit=>true,p...),xs...)
function ds(s)
  println(typeof(s),":")
  for f in fieldnames(typeof(s))
    if !isdefined(s,f) println(f,"=#undef")
    else printc(f,"=",getfield(s,f),"\n")
    end
  end
end

a=1
#--------------------------------------------------------------------------
"""
  a variant of get! for objects O which have a Dict of properties named prop.
  Usually called as
    gets(O,:p) do ---code to compute property :p --- end
"""
gets(f::Function,o,p::Symbol)=get!(f,o.prop,p)

"""
A  variation where it is assumed f(o) sets o.prop[p] but not assumed that f
returns o.prop[p], because  f could  set several keys at once...
"""
function getp(f::Function,o,p::Symbol)
  if haskey(o.prop,p) return o.prop[p] end
  f(o)
  o.prop[p]
end
#----------------------- Formatting -----------------------------------------
const supchars  =
 "-0123456789+()=abcdefghijklmnoprstuvwxyzABDEGHIJKLMNORTUVWβγδειθφχ"
const unicodesup=
 "⁻⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁽⁾⁼ᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᴬᴮᴰᴱᴳᴴᴵᴶᴷᴸᴹᴺᴼᴿᵀᵁⱽᵂᵝᵞᵟᵋᶥᶿᵠᵡ"
const supclass="["*supchars*"]"
const sup=Dict(zip(supchars,unicodesup))
const subchars  ="-0123456789,+()=aehijklmnoprstuvxβγρφχ"
const unicodesub="₋₀₁₂₃₄₅₆₇₈₉‚₊₍₎₌ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵦᵧᵨᵩᵪ"
const sub=Dict(zip(subchars,unicodesub))
const subclass="["*subchars*"]"

"strip TeX formatting from  a string, using unicode characters to approximate"
function TeXstrip(s::String)
  s=replace(s,r"\$"=>"")
  s=replace(s,r"\\varepsilon"=>"ε")
  s=replace(s,r"\\beta"=>"β")
  s=replace(s,r"\\delta"=>"δ")
  s=replace(s,r"\\gamma"=>"γ")
  s=replace(s,r"\\iota"=>"ι")
  s=replace(s,r"\\lambda"=>"λ")
  s=replace(s,r"\\phi"=>"φ")
  s=replace(s,r"\\Phi"=>"Φ")
  s=replace(s,r"\\psi"=>"ψ")
  s=replace(s,r"\\rho"=>"ρ")
  s=replace(s,r"\\sigma"=>"σ")
  s=replace(s,r"\\theta"=>"θ")
  s=replace(s,r"\\chi"=>"χ")
  s=replace(s,r"\\zeta"=>"ζ")
  s=replace(s,r"\\otimes"=>"⊗")
  s=replace(s,r"\\tilde A"=>"Ã")
  s=replace(s,r"\\tilde D"=>"D̃")
  s=replace(s,r"\\times"=>"×")
  s=replace(s,r"\\BZ"=>"ℤ")
  s=replace(s,r"\\frakS"=>"𝔖")
  s=replace(s,r"\\wedge"=>"∧")
  s=replace(s,r"\\#"=>"#")
  s=replace(s,r"\\hbox{([^}]*)}"=>s"\1")
  s=replace(s,r"\\!"=>"")
  s=replace(s,r"{}"=>"")
  s=replace(s,r"\^\{1//2\}"=>"½")
  s=replace(s,r"\^\{-1//2\}"=>"⁻½")
  s=replace(s,r"\^\{1//3\}"=>"⅓")
  s=replace(s,r"\^\{2//3\}"=>"⅔")
  s=replace(s,r"\^\{1//4\}"=>"¼")
  s=replace(s,Regex("_$subclass")=>t->sub[t[2]])
  s=replace(s,Regex("(_\\{$subclass*\\})('*)")=>s"\2\1")
  s=replace(s,Regex("_\\{$subclass*\\}")=>t->map(x->sub[x],t[3:end-1]))
  s=replace(s,Regex("\\^$supclass")=>t->sup[t[2]])
  s=replace(s,Regex("\\^\\{$supclass*\\}")=>t->map(x->sup[x],t[3:end-1]))
  q(l)=l==1 ? "′" : l==2 ? "″" : l==3 ? "‴" : l==4 ? "⁗" : map(x->sup[x],"($l)")
  s=replace(s,r"''*"=>t->q(length(t)))
  s=replace(s,r"\{\+\}"=>"+")
  s
end

function format_coefficient(c::String)
  if c=="1" ""
  elseif c=="-1" "-"
  elseif occursin(r"[-+*/]",c[nextind(c,0,2):end]) "("*c*")" 
  else c end
end

function fromTeX(io::IO,n::String)
  if get(io,:TeX,false) return n 
  elseif get(io,:limit,false) return TeXstrip(n) end
  n=replace(n,r"\\tilde *"=>"~")
  n=replace(n,"_"=>"")
  n=replace(n,"}"=>"")
  n=replace(n,"{"=>"")
  n=replace(n,"\\phi"=>"phi")
  n=replace(n,"\\zeta"=>"E")
  n=replace(n,r"\bi\b"=>"I")
end

fromTeX(n::String;opt...)=fromTeX(IOContext(stdout,opt...),n)

printTeX(io::IO,s...)=print(io,fromTeX(io,prod(s)))

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

# show/print with attributes...
rshow(x;p...)=show(IOContext(stdout,:limit=>true,p...),"text/plain",x)
rsprint(x;p...)=sprint(show,x;context=IOContext(stdout,:limit=>true,p...))

function joindigits(l::AbstractVector,delim="()";always=false,sep=",")
  big=any(l.>=10)
  s=big ? join(l,sep) : join(l)
  (big || always)&& !isempty(delim) ? delim[1]*s*delim[2] : s
end

"""
 cut(string;options)

   options:
   - width=displaysize(stdout)[2]-2 cutting width
   - after=","                      cutting after these chars
   - before=""                      cutting before these chars
   - file=stdout                    where to print result
"""
function cut(s;width=displaysize(stdout)[2]-2,after=",",before="",file=stdout)
  a=split(s)
  if a[end]=="" a=a[1:end-1] end
  for s in a
    l=0
    pa=pb=0
    pos=1
    for (i,c) in pairs(s)
      n=textwidth(c)
      if l+n>width
#       println("pos=$pos pa=$pa pb=$pb l=$l n=$n i=$i")
        if pa>0 && (pb==0 || pa>=pb)
          pr=s[pos:pa]
          pos=nextind(s,pa)
        elseif pb>0
          pr=s[pos:prevind(s,pb)]
          pos=pb
        else error("could not cut ",s[pos:i])
        end
        println(file,pr)
        l-=textwidth(pr)
        pa=pb=0
      end  
      l+=n
      if c in after pa=i end
      if c in before pb=i end
    end
    println(file,s[pos:end])
  end
  if file!=stdout close(file) end
end

#----------------------- Number theory ---------------------------
" the numbers less than n and prime to n "
function prime_residues(n)
  if n==1 return [0] end
  filter(i->gcd(n,i)==1,1:n-1) # inefficient
end

# make Primes.factor fast for small Ints by memoizing it
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
phi(m::Integer)=Primes.totient(m)

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
end
