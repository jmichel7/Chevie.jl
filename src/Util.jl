"""
This  module contains  various utility  functions used  in the  rest of the
code.  Maybe some  of them  exist in  some Julia  module I am not aware of;
please tell me.

The code is divided in sections  according to semantics.
"""
module Util

export 
  @forward,
  getp, @GapObj, # helpers for GapObjs
  showtable, format_coefficient, ordinal, fromTeX, printTeX, joindigits, cut, 
  rio, xprint, xprintln, ds, xdisplay, TeX, TeXs, # formatting
  exactdiv, factor, prime_residues, divisors, phi, primitiveroot #number theory

export toL, toM # convert Gap matrices <-> Julia matrices
export InfoChevie

const info=Ref(true)
function InfoChevie(a...)
  if Util.info[] xprint(a...) end
end

"""
@forward T.f f1,f2,...
  generates 
f1(a::T,args...)=f1(a.f,args...)
f2(a::T,args...)=f2(a.f,args...)
...
"""
macro forward(ex, fs)
  T, field = esc(ex.args[1]), ex.args[2].value
  fdefs=map(fs.args)do ff
      f= esc(ff)
      quote
        ($f)(a::($T),args...)=($f)(a.$field,args...)
      end
  end
  Expr(:block, fdefs...)
end

#--------------------------------------------------------------------------
toL(m)=collect(eachrow(m)) # to Gap
toM(m)=isempty(m) ? Array{eltype(eltype(m))}(undef,0,1) : permutedims(reduce(hcat,m)) # to julia

#--------------------------------------------------------------------------

"""
A  variation of get! where it is assumed f(o) sets o.p but not assumed that
f returns o.p, because f sets several keys at once...
"""
getp(f::Function,o,p::Symbol)=get!(()->(f(o);o.prop[p]),o.prop,p)

# a GapObj is an object which has a field prop::Dict{Symbol,Any}
# so has fixed fields but can dynamically have new ones
# usage: @GapObj struct ...
macro GapObj(e)
  push!(e.args[3].args,:(prop::Dict{Symbol,Any}))
  if e.args[2] isa Symbol T=e.args[2]
  elseif e.args[2].args[1] isa Symbol T=e.args[2].args[1]
  else T=e.args[2].args[1].args[1]
  end
  esc(Expr(:block,
   e,
   :(Base.getproperty(o::$T,s::Symbol)=hasfield($T,s) ? getfield(o,s) : 
         getfield(o,:prop)[s]),
   :(Base.setproperty!(o::$T,s::Symbol,v)=getfield(o,:prop)[s]=v),
   :(Base.haskey(o::$T,s::Symbol)=haskey(getfield(o,:prop),s)),
   :(Base.get!(f::Function,o::$T,s::Symbol)=get!(f,getfield(o,:prop),s))))
end

#----------------------- Formatting -----------------------------------------
# print with attributes...
rio(io::IO=stdout;p...)=IOContext(io,:limit=>true,p...)
xprint(x...;p...)=print(rio(;p...),x...)
xprintln(x...;p...)=println(rio(;p...),x...)
xdisplay(x;p...)=display(TextDisplay(IOContext(stdout,p...)),x)

function ds(s) # "dump struct"; not recursive like dump
  println(typeof(s),":")
  for f in fieldnames(typeof(s))
    if !isdefined(s,f) println(f,"=#undef")
    else xprintln(f,"=",getfield(s,f))
    end
  end
end

const supchars  =
 "-0123456789+()=abcdefghijklmnoprstuvwxyzABDEGHIJKLMNOPRTUVWαβγδειθφχ"
const unicodesup=
 "⁻⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁽⁾⁼ᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᴬᴮᴰᴱᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾᴿᵀᵁⱽᵂᵅᵝᵞᵟᵋᶥᶿᵠᵡ"
const supclass="["*supchars*"]"
const sup=Dict(zip(supchars,unicodesup))
const subchars  ="-0123456789,+()=aehijklmnoprstuvxβγρφχ."
const unicodesub="₋₀₁₂₃₄₅₆₇₈₉‚₊₍₎₌ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵦᵧᵨᵩᵪ̣."
const sub=Dict(zip(subchars,unicodesub))
const subclass="["*subchars*"]"
const TeXmacros=Dict("bbZ"=>"ℤ", "beta"=>"β", "chi"=>"χ", "delta"=>"δ",
  "gamma"=>"γ", "iota"=>"ι", "lambda"=>"λ", "otimes"=>"⊗ ",
  "par"=>"\n", "phi"=>"φ", "varphi"=>"φ", "Phi"=>"Φ", "psi"=>"ψ", "rho"=>"ρ",
  "sigma"=>"σ", "theta"=>"θ", "times"=>"×", "varepsilon"=>"ε", "wedge"=>"∧",
  "zeta"=>"ζ", "backslash"=>"\\","sqrt"=>"√")

"strip TeX formatting from  a string, using unicode characters to approximate"
function unicodeTeX(s::String)
  s=replace(s,r"\\tilde ([A-Z])"=>s"\1\U303")
  s=replace(s,r"\\tilde *(\\[a-zA-Z]*)"=>s"\1\U303")
  s=replace(s,r"\\hfill\\break"=>"\n")
  s=replace(s,r"\\(h|m)box{([^}]*)}"=>s"\2")
  s=replace(s,r"\\#"=>"#")
  s=replace(s,r"\\!"=>"")
  s=replace(s,r"\^\{\\frac\{1\}\{2\}\}"=>"½")
  s=replace(s,r"\^\{\\frac\{-1\}\{2\}\}"=>"⁻½")
  s=replace(s,r"\^\{\\frac\{1\}\{3\}\}"=>"⅓")
  s=replace(s,r"\^\{\\frac\{2\}\{3\}\}"=>"⅔")
  s=replace(s,r"\^\{\\frac\{1\}\{4\}\}"=>"¼")
  s=replace(s,r"\^\{\\frac\{([-0-9]*)\}\{([0-9]*)\}\}"=>function(t)
             t=split(t[9:end-2],"}{")
             map(x->sup[x],t[1])*"⁄"*map(x->sub[x],t[2])
      end)
  s=replace(s,r"\\mathfrak  *S"=>"\U1D516 ")
  s=replace(s,r"\\([a-zA-Z]+) *"=>t->TeXmacros[rstrip(t[2:end])])
  s=replace(s,r"\$"=>"")
  s=replace(s,r"{}"=>"")
  s=replace(s,Regex("_$subclass")=>t->sub[t[2]])
  s=replace(s,Regex("(_\\{$subclass*\\})('*)")=>s"\2\1")
  s=replace(s,Regex("_\\{$subclass*\\}")=>t->map(x->sub[x],t[3:end-1]))
  s=replace(s,Regex("\\^$supclass")=>t->sup[t[2]])
  s=replace(s,Regex("\\^\\{$supclass*\\}")=>t->map(x->sup[x],t[3:end-1]))
  q(l)=l==1 ? "′" : l==2 ? "″" : l==3 ? "‴" : l==4 ? "⁗" : map(x->sup[x],"($l)")
  s=replace(s,r"''*"=>t->q(length(t)))
  s=replace(s,r"\{([^}]*)\}"=>s"\1")
  s
end

function bracket_if_needed(c::String;allow_frac=false)
  ok="([^-+*/]|√-|{-)*"
  par="(\\([^()]*\\))"
  if match(Regex("^[-+]?$ok$par*$ok\$"),c)!==nothing c
  elseif allow_frac && match(Regex("^[-+]?$ok$par*$ok/+[0-9]*\$"),c)!==nothing c
  else "("*c*")" 
  end
end

function format_coefficient(c::String;allow_frac=false)
  if c=="1" ""
  elseif c=="-1" "-"
  else bracket_if_needed(c;allow_frac)
  end
end

function TeXstrip(n::String) # plain ASCII rendering of TeX code
  n=replace(n,r"\\tilde *"=>"~")
  n=replace(n,r"[_{}$]"=>"")
  n=replace(n,"\\phi"=>"phi")
  n=replace(n,"\\zeta"=>"E")
  n=replace(n,r"\bi\b"=>"I")
  n=replace(n,r"\\mathfrak *"=>"")
end

function fromTeX(io::IO,n::String)
  if     get(io,:TeX,false) n 
  elseif get(io,:limit,false) unicodeTeX(n) 
  else   TeXstrip(n)
  end
end

fromTeX(n::String;opt...)=fromTeX(IOContext(stdout,opt...),n)

TeX(io::IO;k...)=IOContext(io,:TeX=>true,pairs(k)...)
TeX(io::IO,x)=repr(x;context=TeX(io))

function printTeX(io::IO,s...)
  res=""
  for x in s res*=x isa String ? x : TeX(io,x) end
  print(io,fromTeX(io,res))
end

"""
`showtable(io, table::AbstractMatrix; options )`

General  routine to format a table. The  following options can be passed as
properties of the `io` or as keywords.

  - `row_labels`         labels for rows (default `axes(table,1)`)
  - `rows_label`         label for first column (column of row labels)
  - `col_labels`         labels for other columns
  - `rowseps`            line numbers after which to put a separator
  - `rows`               show only these rows
  - `cols`               show only these columns
  - `TeX`                give LaTeX output (useful in Jupyter or Pluto)
  - `column_repartition` display in vertical pieces of sizes indicated
    (default if not `TeX`: take in account `displaysize(io,2)`)

```julia-rep1
julia> m=[1 2 3 4;5 6 7 8;9 1 2 3;4 5 6 7];

julia> showtable(stdout,m)
1│1 2 3 4
2│5 6 7 8
3│9 1 2 3
4│4 5 6 7

julia> labels=["x","y","z","t"];

julia> showtable(stdout,m;cols=2:4,col_labels=labels,rowseps=[0,2,4])
 │y z t
─┼──────
1│2 3 4
2│6 7 8
─┼──────
3│1 2 3
4│5 6 7
─┴──────
```
"""
function showtable(io::IO,t::AbstractMatrix; opt...)
  io=IOContext(io,opt...)
  strip(x)=fromTeX(io,x)
  rows=get(io,:rows,axes(t,1))
  cols=get(io,:cols,axes(t,2))
  t=t[rows,cols]
  row_labels=(strip.(get(io,:row_labels,string.(axes(t,1)))))[rows]
  col_labels=get(io,:col_labels,nothing)
  if col_labels!=nothing col_labels=(strip.(col_labels))[cols] end
  rows_label=strip(get(io,:rows_label,""))
  rowseps=get(io,:rowseps,col_labels!=nothing ? [0] : Int[])
  column_repartition=get(io,:column_repartition,nothing)
  lpad(s,n)=" "^(n-textwidth(s))*s # because Base.lpad does not use textwidth
  rpad(s,n)=s*" "^(n-textwidth(s)) # because Base.rpad does not use textwidth
  t=map(x->x isa String ? x : repr(x; context=io),t)
  TeX=get(io,:TeX,false)
  cols_widths=map(i->maximum(textwidth.(t[:,i])),axes(t,2))
  if !isnothing(col_labels)
    cols_widths=map(max,cols_widths,textwidth.(col_labels))
    if !TeX col_labels=map(lpad,col_labels,cols_widths) end
  end
  labwidth=max(textwidth(rows_label),maximum(textwidth.(row_labels)))
  if !TeX
    rows_label=lpad(rows_label,labwidth)
    row_labels=rpad.(row_labels,labwidth)
  end
  function hline(ci;last=false,first=false)
    if TeX println(io,"\\hline")
    else
    print(io,"\u2500"^labwidth,first ? "\u252C" : last ? "\u2534" : "\u253C")
    print(io,"\u2500"^sum(cols_widths[ci].+1),"\n")
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
     if TeX column_repartition=[length(cols_widths)]
     else column_repartition=cut(1 .+cols_widths,displaysize(io)[2]-labwidth-1)
     end
  end
  ci=[0]
  for k in column_repartition
    ci=ci[end].+(1:k)
    if TeX println(io,"\$\$\n\\begin{array}{c|","c"^length(ci),"}") end
    if !isnothing(col_labels)
      if TeX
        println(io,rows_label,"&",join(col_labels[ci],"&"),"\\\\")
      else println(io,rows_label,"\u2502",join(col_labels[ci]," "))
      end
    end
    if 0 in rowseps hline(ci;first=isnothing(col_labels)) end
    for l in axes(t,1)
      if TeX
        println(io,row_labels[l],"&",join(t[l,ci],"&"),"\\\\")
      else
        println(io,row_labels[l],"\u2502",join(map(lpad,t[l,ci],cols_widths[ci])," "))
      end
      if l in rowseps hline(ci,last=l==size(t,1)) end
    end
    if ci[end]!=length(cols_widths) print(io,"\n") end
    if TeX println(io,"\\end{array}\n\$\$") end
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

function joindigits(l,delim="()";always=false,sep=",")
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
  a=split(s,"\n")
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

TeXs(x;p...)=repr("text/plain",x;context=IOContext(stdout,:TeX=>true,p...))

function TeX(x;p...)
  s=tempname(".")
  open("$s.tex","w")do f
    println(f,"\\documentclass{article}")
    println(f,"\\usepackage{amsmath}")
    println(f,"\\usepackage{amssymb}")
    println(f,"\\begin{document}")
    print(f,TeXs(x;p...))
    println(f,"\\end{document}")
  end
  run(`latex $s.tex`)
  run(`xdvi -expert -s 5 $s.dvi`)
  run(`rm $s.tex $s.aux $s.log $s.dvi`)
end

#----------------------- Number theory ---------------------------
exactdiv(a,b)=a/b  # generic version for fields
function exactdiv(a::Integer,b::Integer) # define for integral domains
  (d,r)=divrem(a,b)
  !iszero(r) ? nothing : d
end

" the numbers less than n and prime to n "
function prime_residues(n)
  if n==1 return [0] end
  filter(i->gcd(n,i)==1,1:n-1) # inefficient
end

# make Primes.factor fast for small Ints by memoizing it
import Primes
const dict_factor=Dict(2=>Primes.factor(Dict,2))
factor(n::Integer)=get!(()->Primes.factor(Dict,n),dict_factor,n)

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
