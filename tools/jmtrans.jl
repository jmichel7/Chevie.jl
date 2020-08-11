# install_tbl() translates ~/gap3-dev/pkg/chevie/tbl to ./tbl
using Gap2Julia

const src=[ 
"cmp4_22", "cmplxg24", "cmplxg25", "cmplxg26", 
"cmplxg27", "cmplxg29", "cmplxg31", "cmplxg32", "cmplxg33", "cmplxg34", 
"coxh3", "coxh4", "coxi", "weyla", "weylbc", "weyld", "weyl2d",
"cox2i", "weyl2e6", "weyl2f4", "weyl3d4",
"weyle6", "weyle7", "weyle8", "weylf4", "weylg2", "exceptio",
"weyl2a", 
"cmplximp"] 

# functions whose translation does not work so they are not translated
# an implementation is in table2.jl, with the useful code from compat3
const exclu=[
 ["CharTable","A"],
 ["HeckeCharTable","A"],
 ["Hk","A"],
 ["CharTable","2A"],
 ["FakeDegree","2A"],
 ["HeckeCharTable","2A"],
 ["UnipotentClasses","2A"],
 ["CharTable","B"],
 ["Hk","B"],
 ["HeckeCharTable","B"],
 ["UnipotentClasses","B"],
 ["CharTable","D"],
 ["CartanMat","D"],
 ["UnipotentClasses","D"],
 ["HeckeCharTable","D"],
 ["Hk","D"],
 ["CharTable","2D"],
 ["HeckeCharTable","2D"],
 ["Discriminant","H4"],
 ["CharInfo","I"],
 ["HeckeCharTable","imp"],
 ["CharTable","imp"],
 ["PowerMaps","imp"],
 ["CartanMat",:(["G25","G26","G29","G31","G32","G34"])]
]

const ok=[:(CHEVIE.AddData), 
    :(CHEVIE.IndirectAddData)
   ]
const ok2=[
    :((CHEVIE.families).HS4),
    :((CHEVIE.families).S5), 
    :((CHEVIE.families).G14), 
    :((CHEVIE.families).S4),
    :((CHEVIE.families).F20),
    :((CHEVIE.families).Y6),
    :((CHEVIE.families).X7),
    :((CHEVIE.families).F42),
    :((CHEVIE.families).G4),
    :((CHEVIE.families).X2),
    :G4_22Helper,
    :G4_22Test,
    :G4_22FetchIndexChars,
    :VFactorSchurElement,
    :PartitionTwoCoreQuotient
   ]

readf(f)=Gap2Julia.myparse(read(homedir()*"/gap3-dev/pkg/chevie/"*f,String),false)

# install("x/f") translates ~/gap3-dev/pkg/chevie/x/f.g" to "x/f.jl"
function install(n)
  println("installing $n")
  l=readf(n*".g")
  open(n*".jl","w")do f 
    for e in l 
      if (e.head==:call && (e.args[1] in ok) &&
          all(p->e.args[2]!=p[1] || (length(p)==2 && e.args[3]!=p[2]),exclu))
        if e.args[3] isa String print(e.args[3])
        else print(join(e.args[3].args,","))
        end
        println(".",e.args[2])
        write(f,"\n",string(Gap2Julia.trans(e)))
       elseif e.head==:(=) && (e.args[1] in ok2)
        println("=",e.args[1])
        write(f,"\n",string(Gap2Julia.trans(e)))
      end
    end
  end
end

function install1(n)
  println("installing $n")
  s=read(homedir()*"/gap3-dev/pkg/chevie/"*n*".g",String)
  open(n*".jl","w")do f 
    write(f,Gap2Julia.trad(s))
  end
end

install_tbl()=for f in src install("tbl/$f") end

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
      if (e.head==:call && (e.args[1] in ok) &&
          all(p->e.args[2]!=p[1] || (length(p)==2 && e.args[3]!=p[2]),exclu))
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
