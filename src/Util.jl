"""
This  module contains  various utility  functions used  in the  rest of the
code.  Maybe some  of them  exist in  some Julia  module I am not aware of;
please tell me.
"""
module Util
using Gapjm: stringexp, stringprime

export @forward 
export showtable, ordinal, fromTeX, printTeX, joindigits, cut, rio, xprint, 
  xprintln, ds, xdisplay, hdisplay, TeX, TeXs, hasdecor # formatting
export toL, toM # convert Gap matrices <-> Julia matrices
export InfoChevie
export cartesian, cart2lin, lin2cart

const info=Ref(true)
function InfoChevie(a...)
  if Util.info[] xprint(a...) end
end

"""
`@forward T.f f1,f2,...`

is a macro which delegates definitions. The above generates 
```
f1(a::T,args...)=f1(a.f,args...)
f2(a::T,args...)=f2(a.f,args...)
...
```
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
#----------------------- Formatting -----------------------------------------
# print with attributes...
hasdecor(io::IO)=get(io,:TeX,false)||get(io,:limit,false)
rio(io::IO=stdout;p...)=IOContext(io,:limit=>true,p...)
xprint(x...;p...)=print(rio(;p...),x...)
xprintln(x...;p...)=println(rio(;p...),x...)
xdisplay(x;p...)=display(TextDisplay(rio(;p...)),x)
function hdisplay(x;p...)
  Docs.HTML()do io
    show(IOContext(io,p...),"text/html",x)
  end
end

function ds(s;types=false) # "dump struct"; not recursive like dump
  println(typeof(s),":")
  for f in fieldnames(typeof(s))
    if !isdefined(s,f) println(f,"=#undef")
    else print(f);
      if types print("::",typeof(getfield(s,f))) end
      xprintln("=",getfield(s,f))
    end
  end
end

const supchars  =
 "-0123456789+()=abcdefghijklmnoprstuvwxyzABDEGHIJKLMNOPRTUVWαβγδειθφχ"
const unicodesup=
 "⁻⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁽⁾⁼ᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᴬᴮᴰᴱᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾᴿᵀᵁⱽᵂᵅᵝᵞᵟᵋᶥᶿᵠᵡ"
const sup=Dict(zip(supchars,unicodesup))
const subchars  ="-0123456789,+()=aehijklmnoprstuvxβγρφχ."
const unicodesub="₋₀₁₂₃₄₅₆₇₈₉‚₊₍₎₌ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵦᵧᵨᵩᵪ̣."
const sub=Dict(zip(subchars,unicodesub))
const TeXmacros=Dict("bbZ"=>"ℤ", "beta"=>"β", "chi"=>"χ", "delta"=>"δ",
  "Delta"=>"Δ","gamma"=>"γ", "iota"=>"ι", "lambda"=>"λ", "Lambda"=>"Λ",
  "nu"=>"ν", "otimes"=>"⊗ ", "par"=>"\n", "phi"=>"φ", "varphi"=>"φ", 
  "Phi"=>"Φ", "psi"=>"ψ", "rho"=>"ρ", "sigma"=>"σ", "theta"=>"θ", 
  "times"=>"×", "varepsilon"=>"ε", "wedge"=>"∧", "hfill"=>" ",
  "zeta"=>"ζ", "rtimes"=>"⋊ ","backslash"=>"\\","sqrt"=>"√")

# defs below are necessary since constant folding is not good enough
const r1=Regex("_[$subchars]")
const r2=Regex("(_\\{[$subchars]*\\})('*)")
const r3=Regex("_\\{[$subchars]*\\}")
const r4=Regex("\\^[$supchars]")
const r5=Regex("\\^\\{[$supchars]*\\}")
const r6=Regex("\\\\overline{([0-9]*)}")

"strip TeX formatting from  a string, using unicode characters to approximate"
function unicodeTeX(s::String)
  if all(x->x in 'a':'z' || x in 'A':'Z' || x in '0':'9',s) return s end
  s=replace(s,r"\\tilde ([A-Z])"=>s"\1\U303")
  s=replace(s,r"\\tilde *(\\[a-zA-Z]*)"=>s"\1\U303")
  s=replace(s,r"\\hfill\\break"=>"\n")
  s=replace(s,r"\\(h|m)box{([^}]*)}"=>s"\2")
  s=replace(s,r"\\!"=>"")
  s=replace(s,r"\^\{\\frac\{([-0-9]*)\}\{([0-9]*)\}\}"=>
    t->stringexp(rio(),Rational(parse.(Int,split(t[9:end-2],"}{"))...)))
  s=replace(s,r"\\#"=>"#")
  s=replace(s,r"\\mathfrak  *S"=>"\U1D516 ")
  s=replace(s,r6=>t->prod(string(x,'\U0305') for x in t[11:end-1]))
  s=replace(s,r"\\([a-zA-Z]+) *"=>t->TeXmacros[rstrip(t[2:end])])
  s=replace(s,r"\$"=>"")
  s=replace(s,r"{}"=>"")
  s=replace(s,r1=>t->sub[t[2]])
  s=replace(s,r2=>s"\2\1")
  s=replace(s,r3=>t->map(x->sub[x],t[3:end-1]))
  s=replace(s,r4=>t->sup[t[2]])
  s=replace(s,r5=>t->map(x->sup[x],t[3:end-1]))
  s=replace(s,r"''*"=>t->stringprime(rio(),length(t)))
# s=replace(s,r"\{([^}]*)\}"=>s"\1")
  s=replace(s,r"^\{([^}{,]*)\}"=>s"\1")
  s=replace(s,r"([^a-zA-Z0-9])\{([^}{,=]*)\}"=>s"\1\2")
  s=replace(s,r"([^a-zA-Z0-9])\{([^}{,=]*)\}"=>s"\1\2")
  s
end

function TeXstrip(n::String) # plain ASCII rendering of TeX code
  n=replace(n,r"\\tilde *"=>"~")
# n=replace(n,r"[_{}$]"=>"") next line faster
  n=String(filter(x->!(x in "_{}\$"),collect(n)))
  n=replace(n,r"\\phi"=>"phi")
  n=replace(n,r"\\zeta"=>"E")
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
  - align                alignment of column (default 'c':centered)

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
  row_labels=(strip.(get(io,:row_labels,string.(axes(t,1)))))[rows]
  t=t[rows,cols]
  col_labels=get(io,:col_labels,nothing)
  if col_labels!=nothing col_labels=(strip.(col_labels))[cols] end
  rows_label=strip(get(io,:rows_label,""))
  rowseps=get(io,:rowseps,col_labels!=nothing ? [0] : Int[])
  column_repartition=get(io,:column_repartition,nothing)
  align=get(io,:align,'c')
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
    if TeX println(io,"\$\$\n\\begin{array}{l|","$align"^length(ci),"}") end
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
        println(io,row_labels[l],"\u2502",join(map(align in "cr" ? lpad : rpad,
                                     t[l,ci],cols_widths[ci])," "))
      end
      if l in rowseps hline(ci,last=l==size(t,1)) end
    end
    if ci[end]!=length(cols_widths) print(io,"\n") end
    if TeX println(io,"\\end{array}\n\$\$") end
  end
end

showtable(t::AbstractMatrix;opt...)=showtable(stdout,t;opt...)

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
 `cut(io::IO=stdout,string;width=displaysize(io)[2]-2,after=",",before="")`

This  function prints to `io` the  string argument cut across several lines
for improved display. It can take the following keyword arguments:
  - width:  the cutting width
  - after:  cut after these chars
  - before: cut before these chars
```julia-rep1
julia> cut(string(collect(1:50)))
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
```
"""
function cut(io::IO,s;width=displaysize(stdout)[2]-2,after=",",before="")
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
        println(io,pr)
        l-=textwidth(pr)
        pa=pb=0
      end  
      l+=n
      if c in after pa=i end
      if c in before pb=i end
    end
    println(io,s[pos:end])
  end
end

cut(s;k...)=cut(stdout,s;k...)

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

#--------------------------------------------------------------------------

if false
# better display of Rationals at the REPL
function Base.show(io::IO, x::Rational)
   show(io, numerator(x))
   if haskey(io,:typeinfo) && isone(denominator(x)) return end
   print(io, "//")
   show(io, denominator(x))
end
end
#--------------------------- Chevie compatibility--------------------------
toL(m)=collect(eachrow(m)) # to Gap
toM(m)=isempty(m) ? Array{eltype(eltype(m))}(undef,0,1) : permutedims(reduce(hcat,m)) # to julia

# The following functions should be eventually evicted by adpoting Julia
# order for products.

"""
`cartesian(a::AbstractVector...)`

A variation on ``Iterators.product` which gives the same result as GAP's
`Cartesian`. `reverse` is done twice to get the same order as GAP.
"""
function cartesian(a::AbstractVector...)
  reverse.(vec(collect.(Iterators.product(reverse(a)...))))
end

"""
`cart2lin([l₁,…,lₙ],[i₁,…,iₙ])` is GAP3 PositionCartesian
returns findfirst(==([i1,…,iₙ]),cartesian(1:l₁,…,1:lₙ))
very fast with 2 Tuple arguments
"""
cart2lin(l,ind)=LinearIndices(Tuple(reverse(l)))[reverse(ind)...]

"""
`lin2cart([l],i)` is GAP3 CartesianAt
returns cartesian(map(j->1:j,l))[i]
"""
lin2cart(dims,i)=reverse(Tuple(CartesianIndices(reverse(Tuple(dims)))[i]))

end
