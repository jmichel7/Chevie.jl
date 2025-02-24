# install_tbl() translates ~/gap3-dev/pkg/chevie/tbl to ./tbl
include("Gap2Julia.jl")
using .Gap2Julia

# files translated from Chevie's tbl directory. compat3 and cmp4_22 are not
const src=[ 
 "cmplxg31", "cmplxg32", "cmplxg33", "cmplxg34", "weyle7", "weyle8"] 

# functions hand-translated to files xxx_t.jl in tbl (because the transpiled
# version does not work, or is too slow, or I had time to translate).
# The useful code from compat3 is in Chevie.jl
const exclu=Dict(
 "CharTable"=>["G31","G34"],
 "PrintDiagram"=>["E7","E8","G31","G32","G33","G34"],
 "SchurElement"=>["G31"],
 "WGraph"=>["E8"])

function exclude(e)
  if !haskey(exclu,e.args[2]) return false end
  l=exclu[e.args[2]]
  if e.args[3] isa Expr && e.args[3].head==:vect
       all(x->x in l,e.args[3].args)
  else e.args[3] in l
  end
end

const ok=[:(CHEVIE.AddData), 
    :(CHEVIE.IndirectAddData)
   ]

# other translated functions.
const ok2=[
 :(CHEVIE.families.S5), 
 :(CHEVIE.families.F42),
 :(CHEVIE.families.G4)
]

readf(f)=Gap2Julia.myparse(read(homedir()*"/gap3-dev/pkg/chevie/"*f,String),false)

#using PrettyPrinting
const pp=false

# install("x/f") translates ~/gap3-dev/pkg/chevie/x/f.g" to "x/f.jl"
function install(n)
  println("installing $n")
  l=readf(n*".g")
  open(n*".jl","w")do f 
    for e in l 
     if (e.head==:call && (e.args[1] in ok) && !exclude(e))
        if e.args[3] isa String print(e.args[3])
        else print(join(e.args[3].args,","))
        end
        println(".",e.args[2])
        if pp pprintln(f,Gap2Julia.trans(e))
        else write(f,string(Gap2Julia.trans(e)),"\n")
        end
      elseif e.head==:(=) && (e.args[1] in ok2)
        println("=",e.args[1])
        if pp pprintln(f,Gap2Julia.trans(e))
        else write(f,string(Gap2Julia.trans(e)),"\n")
        end
      end
    end
  end
end

# just converts n.g to n.jl
function install_local(n)
  println("installing $n")
  l=Gap2Julia.myparse(read(n*".g",String),false)
  open(n*".jl","w") do f
    for e in l print(f,Gap2Julia.trans(e),"\n") end
  end
end

function install1(n)
  println("installing $n")
  s=read(homedir()*"/gap3-dev/pkg/chevie/"*n*".g",String)
  open(n*".jl","w")do f 
    write(f,Gap2Julia.trad(s))
  end
end

function readall()
  l=Expr[]
  for n in src 
    println("reading $n")
    l=vcat(l,readf("tbl/$n.g"))
  end
  l
end

function writeall(l)
  open("tables.jl","w")do f
    for e in l 
      if (e.head==:call && (e.args[1] in ok) && !exclude(e))
        if e.args[3] isa String print(e.args[3])
        else print(join(e.args[3].args,","))
        end
        println(".",e.args[2])
        write(f,"\n",string(Gap2julia.trans(e)))
      elseif e.head==:(=) && (e.args[1] in ok2) 
        println("=",e.args[1])
        write(f,"\n",string(Gap2julia.trans(e)))
      end
    end
  end
end

function install_tbl2()
  Gap2Julia.ftrans[:(+)]=:(∔)
  Gap2Julia.ftrans[:(*)]=:(⨰)
  for f in src install("tbl/$f") end
end

install_tbl()=for f in src install("tbl/$f") end
