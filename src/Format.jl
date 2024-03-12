"""
Chevie   contains  some   extended  formatting   facilities  for  printing,
displaying,  formatting  objects  in  various  ways. For that `Chevie` uses
extensively  `IO` properties.  We have  sevral convenience  functions which
make using `IO` properties easier.

`rio(;d...)`   makes  an  `IO`   stream  which  always   has  the  property
`:limit=>true`,  to mimic the REPL default printing, and has also the extra
properties given by the `d...` keywords. Using this, for instance

`IOContext(stdout,:limit=>true,:compact=>true)` becomes `rio(compact=true)`.

We have versions of display functions which use implicitely `rio`:

`xprint(x...;p...)` is the same as `print(rio(;p...),x...)`. Similarly for
`println`, `display` we have `xprintln`, `xdisplay`.

`xrepr(x;p...)` is the same as `repr(x;context=IOContext(stdout,p...))`.
`xrepr(io,x)` is the same as `repr(x;context=io)`.

```julia-rep1
julia> @Pol q;p=(q^5+1)^2
Pol{Int64}: q¹⁰+2q⁵+1

julia> print(p)
Pol([1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1])
julia> xprint(p)
q¹⁰+2q⁵+1
julia> xprint(p;varname=:x)
x¹⁰+2x⁵+1
julia> repr(2E(5,2)+2E(5,3))
"Cyc{Int64}(-1-root(5))"

julia> xrepr(2E(5,2)+2E(5,3);quadratic=false)
"Cyc{Int64}(2E(5,2)+2E(5,3))"
```
Most  objects in Chevie use  TeX for printing when  given the `IO` property
`:TeX=true`. This is used as the default display in `IJulia` and `Pluto` by
giving the property `:TeX` when defining `Base.show(io::IO,
::MIME"text/html", ...)` for these objects. Continuing the above example:

```julia-rep1
julia> xprint(p;TeX=true)
q^{10}+2q^5+1
```
A  model we often  adopt for displaying  nicely complex objects is to first
write  a nice  display using  `TeX` output.  This can  be used  directly in
`IJulia` and `Pluto`. For other environments, we can compute from the `TeX`
representation a suitable one using the following function:

`fromTeX(io::IO,s)`  takes  a  `TeX`  source  and  tries  to  give the best
possible  rendering on  a given  `IO`. This  uses unicode characters at the
REPL  (if `get(io,:limit,false)==true`).  In the  worse case (`stdout`) all
`TeX` special characters are stripped.

```julia-repl
julia> s="E_6[\\zeta_3]:\\phi_{1,6}"
"E_6[\\zeta_3]:\\phi_{1,6}"

julia> fromTeX(rio(),s)
"E₆[ζ₃]:φ₁‚₆"

julia> fromTeX(stdout,s)
"E6[E3]:phi1,6"
```
`printTeX(io,s)` is the same as `print(io,fromTeX(io,s))`.

Other functions to ease formatting are described below: see `showtable`,
`joindigits`, `ordinal`, `cut`.
"""
module Format
using ..Chevie: stringexp, stringprime
using PermGroups: @GapObj
export showtable, ordinal, fromTeX, printTeX, joindigits, cut, rio, xprint, 
  xprintln, xdisplay, xrepr, TeX, TeXs, hasdecor
#----------------------- Formatting -----------------------------------------
# print with attributes...
hasdecor(io::IO)=get(io,:TeX,false)||get(io,:limit,false)
"""
`rio(io::IO=stdout;p...)` enriches `io` with the attributes in `p`. It
always enriches with `limit=true` to mimic display at the REPL.

Thus `print(rio(),x...)` is like printing `x...` at the REPL.
"""
rio(io::IO=stdout;p...)=IOContext(io,:limit=>true,p...)
"`xprint(x...;p...)` is like `print` but uses the enriched io `rio(;p...)`"
xprint(x...;p...)=print(rio(;p...),x...)
"`xprintln(x...;p...)` is like `println` but uses the enriched io `rio(;p...)`"
xprintln(x...;p...)=println(rio(;p...),x...)
"`xdisplay(x...;p...)` is like `display` but uses the enriched io `rio(;p...)`"
xdisplay(x;p...)=display(TextDisplay(rio(;p...)),x)
"`xrepr(x;p...)` is `repr` using as context `stdout` enriched by `p...`"
xrepr(x;p...)=repr(x;context=IOContext(stdout,p...))
"`xrepr(io,x;p...)` is `repr` using as context `io` enriched by `p...`"
xrepr(io::IO,x;p...)=repr(x;context=IOContext(io,p...))
function hdisplay(x;p...) # what is the idea?
  Docs.HTML()do io
    show(IOContext(io,p...),"text/html",x)
  end
end

const supchars  =
 "-0123456789+()=abcdefghijklmnoprstuvwxyzABDEGHIJKLMNOPRTUVWαβγδειθφχ"
const sup=Dict(zip(supchars,
 "⁻⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁽⁾⁼ᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᴬᴮᴰᴱᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾᴿᵀᵁⱽᵂᵅᵝᵞᵟᵋᶥᶿᵠᵡ"))
const subchars  ="-0123456789,+()=aehijklmnoprstuvxβγρφχ."
const sub=Dict(zip(subchars,
                 "₋₀₁₂₃₄₅₆₇₈₉‚₊₍₎₌ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵦᵧᵨᵩᵪ̣."))
const TeXmacros=Dict("bbZ"=>"ℤ", "beta"=>"β", "chi"=>"χ", "delta"=>"δ",
  "Delta"=>"Δ","gamma"=>"γ", "iota"=>"ι", "lambda"=>"λ", "Lambda"=>"Λ",
  "nu"=>"ν", "otimes"=>"⊗ ", "par"=>"\n", "phi"=>"φ", "varphi"=>"φ", 
  "Phi"=>"Φ", "psi"=>"ψ", "rho"=>"ρ", "sigma"=>"σ", "theta"=>"θ", 
  "times"=>"×", "varepsilon"=>"ε", "wedge"=>"∧", "hfill"=>" ",
  "zeta"=>"ζ", "rtimes"=>"⋊ ","backslash"=>"\\","sqrt"=>"√",
  "sum"=>"Σ", "cap"=>"∩", "mu"=>"μ")

# defs below are necessary since constant folding is not good enough
const r1=Regex("_[$subchars]")
const r2=Regex("(_\\{[$subchars]*\\})('*)")
const r3=Regex("_\\{[$subchars]*\\}")
const r4=Regex("\\^[$supchars]")
const r5=Regex("\\^\\{[$supchars]*\\}")
const r6=Regex("\\\\overline{([0-9]*)}")

#strip TeX formatting from  a string, using unicode characters to approximate
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

"formTeX to document"
function fromTeX(io::IO,n::String)
  if     get(io,:TeX,false) n 
  elseif get(io,:limit,false) unicodeTeX(n) 
  else   TeXstrip(n)
  end
end

fromTeX(n::String;opt...)=fromTeX(IOContext(stdout,opt...),n)

TeX(io::IO;k...)=IOContext(io,:TeX=>true,pairs(k)...)
TeX(io::IO,x)=repr(x;context=TeX(io))

"`printTeX(io,s)` is `print(io,fromTeX(io,s))`"
function printTeX(io::IO,s...)
  res=""
  for x in s res*=x isa String ? x : TeX(io,x) end
  print(io,fromTeX(io,res))
end

@GapObj struct Table
  m::Matrix
end

Table(m;kw...)=Table(m,Dict(kw...))

function Base.show(io::IO, ::MIME"text/html", t::Table)
  show(IOContext(io,:TeX=>true),"text/plain",t)
end

function cpad(s,n)
  ls=n-textwidth(s)
  " "^div(ls,2)*s*" "^div(ls+1,2)
end

function Base.show(io::IO,t::Table)
  io=IOContext(io,t.prop...)
  strip(x)=x isa String ? fromTeX(io,x) : repr(x,context=io)
  rows=get(io,:rows,axes(t.m,1))
  cols=get(io,:cols,axes(t.m,2))
  row_labels=(strip.(get(io,:row_labels,string.(axes(t.m,1)))))
  t=t.m
  col_labels=get(io,:col_labels,nothing)
  if col_labels!=nothing col_labels=(strip.(col_labels)) end
  rows_label=strip(get(io,:rows_label,""))
  row_seps=get(io,:row_seps,col_labels!=nothing ? [0] : Int[])
  column_repartition=get(io,:column_repartition,nothing)
  align=get(io,:align,'r')
  dotzero=get(io,:dotzero,false)
  t=map(x->(!ismissing(x) && x==0 && dotzero) ? "." : strip(x),t)
  TeX=get(io,:TeX,false)
  cols_widths=map(i->maximum(textwidth.(t[rows,i])),axes(t,2))
  if !isnothing(col_labels)
    cols_widths=max.(cols_widths,textwidth.(col_labels))
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
     if TeX column_repartition=[length(cols)]
     else column_repartition=cut(1 .+cols_widths[cols],displaysize(io)[2]-labwidth-1)
     end
  end
  seen=0
  alignf(c)=c=='l' ? rpad : c=='r' ? lpad : cpad
  for k in column_repartition
    ci=cols[seen.+(1:k)]
    seen+=k
    if TeX println(io,"\$\$\n\\begin{array}{l|",
                  align isa String ? align[ci] : align^length(ci),"}") end
    if !isnothing(col_labels)
      if TeX
        println(io,rows_label,"&",join(col_labels[ci],"&"),"\\\\")
      else println(io,rows_label,"\u2502",join(col_labels[ci]," "))
      end
    end
    if 0 in row_seps hline(ci;first=isnothing(col_labels)) end
    for l in rows
      if TeX
        println(io,row_labels[l],"&",join(t[l,ci],"&"),"\\\\")
      else
        println(io,row_labels[l],"\u2502",join(map(ci)do i
       (align isa String ? alignf(align[i]) : alignf(align))(t[l,i],cols_widths[i])*" "
      end))
      end
      if l in row_seps hline(ci,last=l==size(t,1)) end
    end
    if ci[end]!=length(cols_widths) print(io,"\n") end
    if TeX println(io,"\\end{array}\n\$\$") end
  end
end

"""
`showtable(io::IO=stdout, table::AbstractMatrix; keywords)`

General  routine to format a table at  the REPL, or in `IJulia` or `Pluto`.
The elements of `table` and any of the labels in the keywords can be of any
type  and are formatted in the context  of `io`, excepted that a string `s`
is are formatted by `fromTeX(io,s)`. The following options can be passed as
properties of the `io` or as keywords.

  - `row_labels`: labels for rows. A `Vector{Any}` (can be strings),
    default `axes(table,1)`
  - `rows_label`: label for first column of row labels (default none)
  - `col_labels`: labels for other columns (default none)
  - `row_seps`: line numbers after which to put a separator. A `Vector{<:Integer}`, 
    default `[0]` if there are `col_labels` otherwise `Int[]`.
  - `rows`: show only these rows. A `Vector{<:Integer}`, default `axes(table,1)`
  - `cols`: show only these columns. A `Vector{<:Integer}`, 
    default `axes(table,1)`
  - `TeX`: default `false`. If true, give LaTeX output (useful in to give
    nicer output in Jupyter or Pluto)
  - `column_repartition`: a `Vector{<:Integer}`. Display in vertical pieces of
    sizes indicated (useful for `TeX`: otherwise the pieces printed take in
    account `displaysize(io,2)`)
  - `align`:  a character in "lcr": alignment of columns (default 'r')
  - `dotzero`: if `true` replace a '0' by '.' in the table (default false)

```julia-rep1
julia> m=[1 2 3 4;5 6 7 8;9 1 2 3;4 5 6 7];

julia> showtable(m)
1│1 2 3 4
2│5 6 7 8
3│9 1 2 3
4│4 5 6 7

julia> labels=["x","y","z","t"];

julia> showtable(m;cols=2:4,col_labels=labels,row_seps=[0,2,4])
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
  show(io,Table(t;opt...))
end

showtable(t::AbstractMatrix;opt...)=showtable(stdout,t;opt...)

"""
`ordinal(n::Integer)`

string for an ordinal number respecting english syntax.

```julia-repl
julia> ordinal(201)
"201st"

julia> ordinal(202)
"202nd"

julia> ordinal(203)
"203rd"

julia> ordinal(204)
"204th"
```julia-repl
"""
function ordinal(n::Integer)
  str=repr(n)
  if     n%10==1 && n%100!=11 str*="st"
  elseif n%10==2 && n%100!=12 str*="nd"
  elseif n%10==3 && n%100!=13 str*="rd"
  else                        str*="th"
  end
  str
end

joindigits(l,d...;k...)=joindigits(Vector{Int}(l),d...;k...)

"""
`joindigits(l::AbstractVector{Int},delim="()";sep=",")`

print  a list `l` of  (usually small) numbers as  compactly as possible: no
separators if all numbers are smaller than 10.

```julia-repl
julia> joindigits([1,9,3,5])
"1935"

julia> joindigits([1,10,3,5])
"(1,10,3,5)"

julia> joindigits([1,10,3,5],"[]";sep="-")
"[1-10-3-5]"
```
"""
function joindigits(l::AbstractVector{Int},delim="()";always=false,sep=",")
  big=any(>=(10),l)
  s=big ? join(l,sep) : String('0'.+l)
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

function TeX(y...;p...)
  s=tempname(".")
  open("$s.tex","w")do f
    println(f,"\\documentclass{article}")
    println(f,"\\usepackage{amsmath,amssymb,relsize}")
    println(f,"\\begin{document}")
    for x in y print(f,x isa String ? x : TeXs(x;p...)) end
    println(f,"\\end{document}")
  end
  run(`latex $s.tex`)
  run(`xdvi -expert -s 5 $s.dvi`)
  run(`rm $s.tex $s.aux $s.log $s.dvi`)
end
TeX(y::Tuple;p...)=TeX(y...;p...)

end
