"""
This  module provides the function  trad(s::String) which takes as argument
some GAP3 code, and returns a string which is some valid julia code.

To run that Julia code it may be helpful to use gap3support.jl which
contains Julia implementations of some common Gap3 functions.
"""
module Gap2Julia

export trad

# massaging GAP so julia can parse it
function lex(s::String)
  while true #suppress comments makes some below changes harmless
    new=replace(s,r"^\s*#.*$"m=>s"") 
    new=replace(new,r"^([^\"]*)#.*$"m=>s"\1") 
    if new==s break end
    s=new
  end
  s=replace(s,r"\\\n"s=>"")
  s=replace(s,r",\s*,"s=>",nothing,")
  s=replace(s,r",\s*,"s=>",nothing,")
  s=replace(s,r"\[\s*,"s=>"[nothing,")
  s=replace(s,r",\s*\]"s=>",nothing]")
  s=replace(s,r"([^<>:=])=([^=])"=>s"\1==\2") 
  s=replace(s,r"\.\\[-+*]"=>".\"\1\"")
  s=replace(s,r":="=>"=")
  s=replace(s,r"<>"=>"!=")
  s=replace(s,r"return\s*\n"=>" return ")
  s=replace(s,r";\s*fi\s*;"s=>" end;")
  s=replace(s,r";\s*fi\b"s=>" end;")
  s=replace(s,r"\bfi\s*;"s=>" end;")
  s=replace(s,r";\s*od\s*;"s=>" end;")
  s=replace(s,r";\s*od\s+"s=>" end;")
  s=replace(s,r"\bod\s*;"s=>" end;")
  s=replace(s,r"\bthen\b"=>" ")
  s=replace(s,r"\btype\b"=>"type_")
  s=replace(s,r"\bdo\b"=>" ")
  s=replace(s,r"\bnot\b"=>" ! ")
  s=replace(s,r"\band\b"=>" && ")
  s=replace(s,r"\bor\b"=>" || ")
  s=replace(s,r"\bmod\b"=>" % ")
  s=replace(s,r"\btry\b"=>"try_")
  s=replace(s,r"\)in"=>") in")
  #s=replace(s,r"local\s*\w*\s*;"=>" ")
  #s=replace(s,r"local(\s*\w*\s*,)*\s*\w*\s*;"=>" ")
  s=replace(s,r";"=>"\n")
  s=replace(s,r"elif"=>"elseif")
  s=replace(s,r"([^\w\"]\s*)((\s*\((\s*\d+\s*,)+\s*\d+\))+)"=>s"\1perm\"\2\"")
  s=replace(s,r"\)\s*\["=>")[")
  s=replace(s,r"\n/"=>"/\n")
  s=replace(s,r"(?<!,)\n\s*\-"=>"-\n")
  s=replace(s,r"(?<!,)\n\s*\+"=>"+\n")
  s=replace(s,r"(?<!,)\n\s*\^"=>"^\n")
  s=replace(s,r"(?<!,)\n\s*\*"=>"*\n")
  s=replace(s,r"(?<!,)\n\s*\/"=>"/\n")
  s=replace(s,r"\n\s*\&\&"=>"&&\n")
  s=replace(s,r", *\-\n"=>",\n-")
  s=replace(s,r"\[ *\-\n"=>"[\n-")
  s=replace(s,r"\^\-\n"=>"^\n-")
  s=replace(s,r"\)\s*\{"=>"){")
  s=replace(s,r"\)\s*\("=>")(")
  s=replace(s,r"\}\s*\{"=>"}{")
  s=replace(s,r"(\w)\s*\("=>s"\1(")
  s=replace(s,r"[^\\]\\c"=>"")
  s=replace(s,r"\$"=>"\\\$")
  s=replace(s,r"(\w+)\.(\d+)"=>s"\1[:\2]")
  s=replace(s,r"repeat"=>"while true")
  s=replace(s,r"until (.*)\n"=>s"if \1 break end\nend")
end

function myparse(s::String,debug=true)
  s=lex(s)
  i=1
  l=[]
  while i<length(s)
    p=1
  try
    p=Meta.parse(s,i)
  catch
    println(s[i:min(i+1000,length(s))])
    write("test",s[i:end])
    rethrow()
  end
  # print(p)
    if debug 
      print("++++++++++++++++",s[i:min(p[2],length(s))],"++++++++++++++++")
    end
    i=p[2]
    push!(l,p[1])
  end
  l
end

ftrans=Dict{Symbol,Symbol}(
  :(/)=>:(//),
  :(%)=>:mod,
  :AbsInt=>:abs,
  :Add=>:push!,
  :Arrangements=>:arrangements,
  :AssociatedPartition=>:conjugate_partition,
  :Basis=>:basis,
  :BetaSet=>:βset,
  :Binomial=>:binomial,
  :Cartesian=>:cartesian,
  :Collected=>:tally,
  :Combinations=>:combinations,
  :Copy=>:deepcopy,
  :CoxeterCoset=>:spets,
  :CoxeterWord=>:word,
  :Cycles=>:orbits,
  :CycPolFakeDegreeSymbol=>:fegsymbol,
  :DefectSymbol=>:defectsymbol,
  :Degree=>:degree,
  :Denominator=>:denominator,
  :DivisorsInt=>:divisors,
  :Dominates=>:dominates,
  :Elements=>:elements,
  :Error=>:error,
  :Factorial=>:factorial,
  :FullSymbol=>:fullsymbol,
  :FusionAlgebra=>:fusion_algebra,
  :GaloisCyc=>:galois,
  :Gcd=>:gcd,
  :Hasse=>:hasse,
  :Hecke=>:hecke,
  :HighestPowerFakeDegreeSymbol=>:degree_fegsymbol,
  :HighestPowerGenericDegreeSymbol=>:degree_gendeg_symbol,
  :InclusionGens=>:inclusiongens,
  :InfoChevie2=>:InfoChevie,
  :IntListToString=>:joindigits,
  :Length=>:length,
  :LowestPowerFakeDegreeSymbol=>:valuation_fegsymbol,
  :LowestPowerGenericDegreeSymbol=>:valuation_gendeg_symbol,
  :MatXPerm=>:reflrep,
  :Maximum=>:maximum,
  :NrPartitions=>:npartitions,
  :NrPartitionTuples=>:npartition_tuples,
  :Numerator=>:numerator,
  :OrderPerm=>:order,
  :Parent=>:parent,
  :PartBeta=>:partβ,
  :Partitions=>:partitions,
  :PartitionTuples=>:partition_tuples,
  :PermutationMat=>:Matrix,
  :Phi=>:phi,
  :PositionClass=>:position_class,
  :PrimeResidues=>:prime_residues,
  :Print=>:print,
  :QuoInt=>:div,
  :Rank=>:rank,
  :RankSymbol=>:ranksymbol,
  :RecFields=>:keys,
  :ReflectionEigenvalues=>:refleigen,
  :ReflectionType=>:refltype,
  :Reversed=>:reverse,
  :RootsCartan=>:roots,
  :Rotation=>:circshift,
  :SemisimpleRank=>:semisimplerank,
  :Set=>:gapSet,
  :ShallowCopy=>:copy,
  :ShiftBeta=>:shiftβ,
  :SignedPerm=>:SPerm,
  :SignedPermListList=>:SPerm,
  :SignInt=>:sign,
  :Sort=>:sort!,
  :Spets=>:spets,
  :SpechtModel=>:Spechtmodel,
  :String=>:string,
  :StringSymbol=>:stringsymbol,
  :SymbolPartitionTuple=>:symbol_partition_tuple,
  :Symbols=>:BDSymbols,
  :Tableaux=>:tableaux,
  :Torus=>:torus,
  :Twistings=>:twistings,
)

# functional pattern f(list,fun)->ftrans1[f](fun,list)
ftrans1=Dict{Symbol,Symbol}(
  :List=>:map,
  :Number=>:count,
  :ForAny=>:any,
  :ForAll=>:all,
)
function trans(e)
  sym=x->Core.QuoteNode(Symbol(x))
  if e isa Expr
    args=e.args
    if e.head!=:macrocall args=filter(x->!(x isa LineNumberNode),args) end
    nargs=length(args)
    if nargs==1 (a,)=args
    elseif nargs==2 (a,b)=args
    elseif nargs>2 (a,b,c)=args 
    end # first 3 arguments captured to a,b,c
    head=e.head
 #  print("head=$head args=$args\n")
 #  args=filter(x->x isa Expr || x isa String,args)
    if e.head==:call 
      if haskey(ftrans,a) 
        args[1]=ftrans[a]
        a=args[1]
      end
      if a==:(CHEVIE.AddData)
        if nargs!=4  error("unexpected") end
        return :(chevieset($(sym(c)),$(sym(b)),$(trans(args[4]))))
      elseif a==:(CHEVIE.R)
        if nargs!=3  error("unexpected") end
        if c isa String c=sym(c) end
        return :(chevieget($c,$(sym(b))))
      elseif a==:(CHEVIE.IndirectAddData)
        return :(chevieset($c,$(sym(b)),$(trans(args[4]))))
      elseif a==:rec
        args[1]=:(Dict{Symbol,Any})
        args[2:end]=map(args[2:end]) do f
          f=trans(f)
          if f.args[1] isa Expr || f.args[1] isa Symbol
               Expr(:call,:(=>),Core.QuoteNode(f.args[1]),f.args[2])
          else Expr(:call,:(=>),Core.QuoteNode(Symbol(f.args[1])),f.args[2])
          end
        end
        return Expr(head,args...)
      elseif a==:Append return :($(trans(b))=Append($(trans(b)),$(trans(c))))
#     elseif a==:ApplyFunc 
#       if c isa Symbol args=[b,Expr(:...,c)]
#       else args=vcat([b],c.args)
#       end
      elseif a==:IsBound && b isa Expr
        if b.head==:. return Expr(:call,:haskey,b.args[1],b.args[2])
        elseif b.head==:ref return :($b!==nothing)
        end
      elseif a==:IsFunc args=[:isa,b,:Function]
      elseif nargs==3 && haskey(ftrans1,a) args=[ftrans1[a],c,b]
      elseif a==:Unbind && b.head==:.
        return Expr(:call,:delete!,trans(b.args[1]),b.args[2])
      elseif a==:PrintToString 
        return :($(trans(b))*=$(Expr(:call,:SPrint,map(trans,args[3:end])...)))
      elseif a==:string && nargs==3 args[1]=:pad
      elseif a==:vcat && nargs==2
        return :(vcat(($(trans(b)))...))
      elseif a==:Zip args=[:map,args[4],b,c]
      elseif a==:^ && (b isa Expr) && b.head==:call && b.args[1]==:E
        return :(E($(trans(b.args[2])),$(trans(c))))
      elseif a==:* && c==0 && b isa Expr && b.head== :vect && b.args[1].args[1]== :..
        v=b.args[1].args
        return :(fill(0,max(0,1+$(trans(v[3]))-$(trans(v[2])))))
      elseif a==:* && b==0 && c isa Expr && c.head== :vect && c.args[1].args[1]== :..
        v=c.args[1].args
        return :(fill(0,max(0,1+$(trans(v[3]))-$(trans(v[2])))))
      end
    elseif head==:vect
      if nargs>=1 && a isa Expr && a.head==:call && a.args[1]==:(..)
        return Expr(:call,:(:),trans(a.args[2]),trans(a.args[3]))
      elseif nargs>=2 && b isa Expr && b.head==:call && b.args[1]==:(..)
       return :($(trans(a)):$(trans(b.args[2]))-$(trans(a)):$(trans(b.args[3])))
      end
    elseif head==:curly head=:ref
    elseif head==:function 
      if a==:((arg,)) args[1]=:((arg...,)) end
      return Expr(head,args[1],map(trans,args[2:end])...)
    elseif head==:tuple 
      return Expr(:call,:Perm,map(trans,args)...)
    elseif head==:. 
      if b isa Expr return Expr(:ref,trans(a),:(Symbol($(trans(b.args[1])))))
      else return Expr(:ref,trans(a),b)
      end
    end 
 #  print("head=$head map(trans,args)=$(map(trans,args))\n")
    return Expr(head,map(trans,args)...)
  elseif e isa Symbol
    if e in keys(ftrans) return ftrans[e]
    else return e
    end
  else return e
  end
end

trad(s)=join(string.(Base.remove_linenums!.(trans.(myparse(s,false)))),"\n")
end
