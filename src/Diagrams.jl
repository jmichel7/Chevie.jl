module Diagrams

export Diagram, diagram

using ..Chevie
struct Diagram
  t::TypeIrred
end

function Base.show(io::IO, ::MIME"text/html", v::Vector{Diagram})
  show(IOContext(io,:TeX=>true), "text/plain",v)
end

Base.show(io::IO,::MIME"text/plain",v::Vector{Diagram})=show(io,v)

Base.show(io::IO,v::Vector{Diagram})=join(io,v,"\n")

function Base.show(io::IO,d::Diagram)
  if !(get(io,:limit,false)|| get(io,:TeX,false))
    print(io,"Diagram(",d.t,")")
    return
  end
  t=d.t
  if haskey(t,:orbit)
    act=length(t.orbit)>1
    if act
      println(io,"ϕ permutes the next ",length(t.orbit)," components")
    end
    if !isone(t.twist)
     println(io,"ϕ",act ? "^$(length(t.orbit))" : "",
      " acts as ",t.twist^mappingPerm(t.orbit[1].indices,1:rank(t.orbit[1])),
      " on the component below")
    end
    show(io,Diagram.(t.orbit))
    return
  end
  series=t.series::Symbol
  if series!=:ST show(io,d,Val(series))
  elseif haskey(t,:ST) show(io,d,Val(Symbol("G",t.ST)))
  else impdiagram(io,d,t.p,t.q,rank(t),t.indices)
  end
end

"""
`diagram(W)`  prints  a  diagram  describing  a  presentation of the finite
reflection group or spets `W`
```julia-repl
julia> diagram(coxgroup(:E,8))
    O 2
    ￨
O—O—O—O—O—O—O E₈
1 3 4 5 6 7 8

julia> diagram(crg(33))
    3 ②       G₃₃
     /^\\
② ——② ——② ——② 
1   2   4   5     423423=342342
```
"""
diagram(W)=Diagram.(refltype(W))

hbar="\u2014"
rdarrow(n)="\u21D0"^(n-1)*" "
ldarrow(n)="\u21D2"^(n-1)*" "
tarrow(n)="\u21DB"^(n-1)*" "
dbar="\u2550"
tbar(n)="\u2261"^(n-1)*" "
vbar="\UFFE8" # "\u2503"
function cd(i)
  if i<=10 ('\u2460'+i-1)*" "
  else "($i)"
  end
end
node="O"

script(i)=string("{\\scriptstyle ",i,"}")
rlap(s)="\\rlap{"*s*"}"
kern(n)=string("\\kern",n,"pt")
nnode(i)=string(kern(-2),"\\mathop\\bigcirc\\limits_{",i,"}",kern(-2))
# circle o with i beneath
ncnode(o,i)=string(kern(-2),rlap("\\hbox{\$"*kern(4.5)*script(o)*"\$}"),
                   "\\mathop\\bigcirc\\limits_{",i,"}",kern(-2),"\n")
# circle o with i on right
rcnode(o,i)=kern(-0.4)*"\\bigcirc"*kern(-7)*script(o)*kern(4.6)*script(i)*"\n"
# circle o with i on left
lcnode(o,i)=script(i)*kern(0.6)*"\\bigcirc"*kern(-9.3)*script(o)*kern(3.6)*"\n"

rule(d,w,h)=string("\\rule[",d,"pt]{",w,"pt}{",h,"pt}\n")
barr(n)=rule(2,n,1)
dbarr(n)=rlap(rule(3,n,1))*rule(1,n,1)
const tbarr=rlap(rule(4,10,1))*rlap(barr(10))*rule(0,10,1)

raise(n,s)=string("\\raise",n,"pt\\hbox{\$",s,"\$}")

vertbar(x,y)=kern(1.5)*rlap("\\hbox{\$"*nnode(x)*"\$}")*"\n"*
    rlap(kern(2.6)*rule(6.6,1,10.8))*"\n"*
    rlap(raise(19.4,kern(-2)*"\\bigcirc"*script(y)))*"\n"*kern(8)

larger(s)="\\mathlarger{"*s*"}"
larger(s,i)=i==0 ? larger(s) : larger(larger(s,i-1))

function getlind(d)
  t=d.t
  indices=t.indices
  if isnothing(indices) ind=fill("?",rank(t))
  else ind=repr.(indices)
  end
  length.(ind),ind
end

function showrel(io,d::Diagram,i...)
  R=braid_relations(d.t)[collect(i)]
  r(a,b)=joindigits(R[a][b])
  if length(i)==2 && R[1][1]==R[2][1]
    print(io,"  ",r(1,1),"=",r(1,2),"=",r(2,2))
  else for a in eachindex(i) print(io,"  ",r(a,1),"=",r(a,2)) end
  end
end
  
function Base.show(io::IO,d::Diagram,::Val{:A})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    println(io,"\$\$",nnode(ind[1]))
    for i in 2:length(l) println(io,barr(10),nnode(ind[i])) end
    println(io,"\$\$")
  else join(io,node.*hbar.^l[1:end-1]);print(io,node);println(io," ",d.t)
    join(io,ind," ")
  end
end

function Base.show(io::IO,d::Diagram,::Val{:B})
  l,ind=getlind(d)
# c=haskey(d.t,:cartanType) ? d.t.cartanType : 1
  c=d.t.cartanType
  if get(io,:TeX,false)
    println(io,"\$\$",nnode(ind[1]))
    if c==2 print(io,rlap(raise(5,"\\leftarrow")),dbarr(10))
    elseif c==1 print(io,rlap(raise(5,"\\rightarrow")),dbarr(10))
    elseif c==root(2) print(io,dbarr(10))
    else print(io,rlap(raise(5,script(xrepr(io,c)))),dbarr(10))
    end
    println(io,nnode(ind[2]))
    for i in 3:length(l) println(io,barr(10),nnode(ind[i])) end
    println(io,"\$\$")
  else print(io,node)
    if c==2 l1=max(l[1],2);print(io,rdarrow(l1))
    elseif c==1 l1=max(l[1],2);print(io,ldarrow(l1))
    elseif c==root(2) l1=max(l[1],2);print(io,dbar^l1)
    else print(io,"=",c,"="); l1=length(xrepr(rio(),c))+2
    end
    join(io,node.*hbar.^l[2:end-1]);println(io,node," ",d.t)
    print(io,rpad(ind[1],l1+1));join(io,ind[2:end]," ")
  end
end

function Base.show(io::IO,d::Diagram,::Val{:D})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",nnode(ind[1]),barr(10))
    print(io,vertbar(ind[3],ind[2]))
    for i in 4:length(l) print(io,barr(10),nnode(ind[i])) end
    print(io,"\$\$")
  else println(io," "^(l[1]+1),node," $(ind[2])")
    println(io," "^(l[1]+1),vbar)
    print(io,node,map(l->hbar^l*node,l[[1;3:end-1]])...)
    println(io," ",d.t)
    join(io,ind[[1;3:end]]," ")
  end
end

function Base.show(io::IO,d::Diagram,::Val{:E})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",nnode(ind[1]),barr(10))
    print(io,nnode(ind[3]),barr(10))
    print(io,vertbar(ind[4],ind[2]))
    for i in 5:length(l) print(io,barr(10),nnode(ind[i])) end
    print(io,"\$\$")
  else println(io," "^(2+l[1]+l[3]),node," $(ind[2])")
    println(io," "^(2+l[1]+l[3]),vbar)
    print(io,node,map(l->hbar^l*node,l[[1;3:end-1]])...)
    println(io," ",d.t)
    join(io,ind[[1;3:end]]," ")
  end
end

function Base.show(io::IO,d::Diagram,::Val{:F})
  l,ind=getlind(d)
# c=haskey(d.t,:cartanType) ? d.t.cartanType : 1
  c=d.t.cartanType
  if get(io,:TeX,false)
    print(io,"\$\$",nnode(ind[1]))
    print(io,barr(10),nnode(ind[2]))
    if c==1 print(io,rlap(raise(5,"\\rightarrow")),dbarr(10))
    elseif c==root(2) print(io,dbarr(10))
    else print(io,rlap(raise(5,script(xrepr(io,c)))),dbarr(10))
    end
    print(io,nnode(ind[3]))
    print(io,barr(10),nnode(ind[4]))
    print(io,"\$\$")
  else print(io,node,hbar^l[1],node)
    if c==1 l1=max(l[2],2);print(io,ldarrow(l1))
    else l1=max(l[2],1);print(io,dbar^l1)
    end
    println(io,node,hbar^l[3],node," ",d.t)
    print(io,ind[1]," ",rpad(ind[2],l1+1),ind[3]," ",ind[4])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G})
  l,ind=getlind(d)
# c=haskey(d.t,:cartanType) ? d.t.cartanType : 1
  c=d.t.cartanType
  if get(io,:TeX,false)
    println(io,"\$\$",nnode(ind[1]))
    if c==1 print(io,rlap(raise(5,"\\rightarrow")),tbarr)
    elseif c==root(3) print(io,dbarr(10))
    else error("not implemented")
    end
    print(io,nnode(ind[2]),"\$\$")
  else
    if c==1 print(io,node,tarrow(max(l[1],2)))
    else print(io,node,tbar(max(l[1],2)))
    end
    println(io,node," ",d.t)
    print(io,ind[1]," "^max(3-l[1],1),ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:H})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",nnode(ind[1]))
    print(io,rlap(kern(3)*raise(6,script(5))),barr(10))
    print(io,nnode(ind[2]))
    for i in 3:length(l) print(io,barr(10),nnode(ind[i])) end
    print(io,"\$\$")
  else println(io," "^l[1],"₅")
    println(io,map(i->node*hbar^l[i],1:length(l)-1)...,node," ",d.t)
    join(io,ind," ")
  end
end

function Base.show(io::IO,d::Diagram,::Val{:I})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    println(io,"\$\$",nnode(ind[1]))
    println(io,rlap(kern(3)*raise(6,script(d.t.bond))))
    println(io,barr(10),nnode(ind[2]))
    println(io,"\$\$")
  else println(io," "^l[1],d.t.bond)
    println(io,node,hbar^l[1],node)
    join(io,ind," ")
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G4})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(3,ind[1]))
    print(io,barr(10),ncnode(3,ind[2]),"\$\$")
  else println(io,cd(3),hbar^2,cd(3),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G5})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(3,ind[1]))
    print(io,dbarr(10),ncnode(3,ind[2]),"\$\$")
  else println(io,cd(3),dbar^2,cd(3),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G6})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(2,ind[1]))
    print(io,tbarr,ncnode(3,ind[2]),"\$\$")
  else println(io,cd(2),tbar(2),cd(3),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G7})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",lcnode(2,ind[1]))
    print(io,kern(-2.5),raise(-3,larger("\\bigcirc",5)),"\n",kern(-5),
          rlap(raise(11,rcnode(3,ind[2]))),
          raise(-10,rcnode(3,ind[3])),"\$\$")
  else println(io," "^2,cd(3),ind[2]," ",d.t);
    println(io," ","/3\\")
    println(io,cd(2),hbar^2,cd(3))
    print(io,ind[1]," "^3,ind[3])
  end
  showrel(io,d,1,2)
end

function Base.show(io::IO,d::Diagram,::Val{:G8})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(4,ind[1]))
    print(io,barr(10),ncnode(4,ind[2]),"\$\$")
  else println(io,cd(4),hbar^2,cd(4),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G9})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(2,ind[1]))
    print(io,tbarr,ncnode(4,ind[2]),"\$\$")
  else println(io,cd(2),tbar(2),cd(4),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G10})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(3,ind[1]))
    print(io,dbarr(10),ncnode(4,ind[2]),"\$\$")
  else println(io,cd(3),dbar^2,cd(4),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G11})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",lcnode(2,ind[1]))
    print(io,kern(-2.5),raise(-3,larger("\\bigcirc",5)),kern(-5),"\n",
          rlap(raise(11,rcnode(3,ind[2]))),
          raise(-10,rcnode(4,ind[3])),"\$\$")
  else println(io," "^2,cd(3),ind[2]," ",d.t);println(io," /3\\");
    println(io,cd(2),hbar^2,cd(4))
    print(io,ind[1]," "^3,ind[3])
  end
  showrel(io,d,1,2)
end

function Base.show(io::IO,d::Diagram,::Val{:G12})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",lcnode(2,ind[1]))
    print(io,kern(-2.5),raise(-3,larger("\\bigcirc",5)*raise(3,kern(-15)*"4")),
          kern(4),"\n",
          rlap(raise(11,rcnode(2,ind[2]))),
          raise(-10,rcnode(2,ind[3])),"\$\$")
  else println(io,"  ",cd(2),ind[2]," ",d.t);
    println(io," /4\\");println(io,cd(2),hbar^2,cd(2))
    print(io,ind[1]," "^3,ind[3])
  end
  showrel(io,d,1,2)
end

function Base.show(io::IO,d::Diagram,::Val{:G13})
  l,ind=getlind(d)
  if get(io,:TeX,false)
   print(io,"\$\$",lcnode(2,ind[1]),kern(-2.5),
         rlap(raise(-3,larger("\\bigcirc",5))),
         rlap(kern(18.5)*raise(0.5,script(4))),
         rlap(kern(13)*raise(0.5,script(5))),
         kern(5),raise(-0.5,larger("\\bigcirc",2)),
         kern(0.3),
         rlap(raise(11,rcnode(2,2))),
              raise(-10,rcnode(2,3)),
              "\$\$")
  else println(io,"  ",cd(2),ind[1]," ",d.t)
    println(io," / \\")
    print(io,cd(2),hbar^2,cd(2))
  end
  showrel(io,d,1,2)
end

function Base.show(io::IO,d::Diagram,::Val{:G14})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(2,ind[1]))
    print(io,rlap(kern(3)*raise(6,script(8))))
    print(io,barr(10),ncnode(3,ind[2]),"\$\$")
  else
    println(io,"  ₈");println(io,cd(2),hbar,cd(3),d.t)
    print(io,ind[1]," "^2,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G15})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    println(io,"\$\$",lcnode(2,ind[1]),kern(-1),
      rlap(kern(-2.5)*raise(-10,script(5))))
    println(io,kern(-2.4),raise(0,"\\bigg ("),kern(-1),"\n",
    rlap(raise(10,rcnode(2,ind[2]))),
    raise(-10,rcnode(3,ind[3])))
    print(io,"\$\$")
  else println(io,"  ",cd(2),ind[1]," ",d.t);
    println(io," /5");println(io,cd(2),ind[3]);
    println(io," \\");print(io,"  ",cd(3),ind[2])
  end
  showrel(io,d,1,2)
end

function Base.show(io::IO,d::Diagram,::Val{:G16})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(5,ind[1]))
    print(io,barr(10),ncnode(5,ind[2]),"\$\$")
  else println(io,cd(5),hbar^2,cd(5),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G17})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(2,ind[1]))
    print(io,tbarr,ncnode(5,ind[2]),"\$\$")
  else println(io,cd(2),tbar(2),cd(5),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G18})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(3,ind[1]))
    print(io,dbarr(10),ncnode(5,ind[2]),"\$\$")
  else println(io,cd(3),dbar^2,cd(5),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G19})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",lcnode(2,ind[1]))
    print(io,kern(-2.5),raise(-3,larger("\\bigcirc",5)),kern(-5),"\n",
          rlap(raise(11,rcnode(3,ind[2]))),
          raise(-10,rcnode(5,ind[3])),"\$\$")
  else println(io," "^2,cd(3),ind[2]," ",d.t);println(io," /3\\");
    println(io,cd(2),hbar^2,cd(5))
    print(io,ind[1]," "^3,ind[3])
  end
  showrel(io,d,1,2)
end

function Base.show(io::IO,d::Diagram,::Val{:G20})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(3,ind[1]))
    print(io,rlap(kern(3)*raise(6,script(5))))
    print(io,barr(10),ncnode(3,ind[2]),"\$\$")
  else println(io,"  ₅");println(io,cd(3),hbar,cd(3),d.t)
    print(io,ind[1]," "^2,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G21})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(2,ind[1]))
    print(io,rlap(kern(3)*raise(6,script(10))))
    print(io,barr(10),ncnode(3,ind[2]),"\$\$")
  else println(io,"  ₁₀");
    println(io,cd(2),hbar^2,cd(3),d.t);print(io,ind[1]," "^3,ind[2])
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G22})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",lcnode(2,ind[1]))
    print(io,kern(-2.5),raise(-3,larger("\\bigcirc",5)*
    raise(3,kern(-15)*"5")),kern(4),"\n",
          rlap(raise(11,rcnode(2,ind[2]))),
          raise(-10,rcnode(2,ind[3])),"\$\$")
  else println(io,"  ",cd(2),ind[2]," ",d.t);println(io," /5\\")
    println(io,cd(2),hbar^2,cd(2))
    print(io,ind[1]," "^3,ind[3])
  end
  showrel(io,d,1,2)
end

function Base.show(io::IO,d::Diagram,::Val{:G24})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(2,ind[2]),dbarr(14),ncnode(2,ind[3]))
    print(io,kern(-26.2),raise(8.5,"\\diagup"),kern(-3.1))
    print(io,rlap(raise(16.7,rcnode(2,ind[1]))))
    print(io,kern(6),raise(8.3,"\\diagdown"))
    print(io,"\$\$")
  else println(io,"  ",cd(2),ind[1]," ",d.t)
    println(io," / \\")
    println(io,cd(2),dbar^2,cd(2))
    print(io,ind[2]," "^3,ind[3])
  end
  showrel(io,d,4)
end

function Base.show(io::IO,d::Diagram,::Val{:G25})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(3,ind[1]))
    for i in 2:3 print(io,barr(10),ncnode(3,ind[i])) end
    print(io,"\$\$")
  else println(io,cd(3),hbar^2,cd(3),hbar^2,cd(3),d.t)
    print(io,prod(map((i,j)->i*" "^(4-j),ind,l)))
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G26})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    println(io,"\$\$",ncnode(2,ind[1]))
    println(io,dbarr(10),ncnode(3,ind[2]))
    println(io,barr(10),ncnode(3,ind[3]),"\$\$")
  else println(io,cd(2),dbar^2,cd(3),hbar^2,cd(3),d.t)
    print(io,prod(map((i,j)->i*" "^(4-j),ind,l)))
  end
end

function Base.show(io::IO,d::Diagram,::Val{:G27})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(2,ind[2]),dbarr(14),ncnode(2,ind[3]))
    print(io,kern(-26.2),raise(8.5,"\\diagup"),kern(-3.1))
    print(io,rlap(raise(16.7,rcnode(2,ind[1]))))
    print(io,kern(6),raise(8.3,"\\diagdown"))
    print(io,"\$\$")
  else println(io,"  ",cd(2),ind[1]," ",d.t)
    println(io, " / \\")
    println(io,cd(2),dbar^2,cd(2))
    print(io,ind[2]," "^3,ind[3])
  end
  showrel(io,d,4)
end

function Base.show(io::IO,d::Diagram,::Val{:G29})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(2,ind[1]),barr(10),ncnode(2,ind[2]))
    print(io,rlap(rlap(raise(5,"\\leftarrow"))*dbarr(14)))
    print(io,rlap(kern(5)*rule(10,1,5)*kern(-2)*rule(10,1,5)))
    print(io,rlap(kern(2)*raise(17,rcnode(2,ind[4]))))
    print(io,rlap(kern(-4)*raise(9,"\\diagup")))
    print(io,rlap(kern(9)*raise(9,"\\diagdown")))
    print(io,kern(14),ncnode(2,ind[3]))
    print(io,"\$\$")
  else println(io," "^6,cd(2),ind[4]," "^2,d.t)
    println(io,"     /\u2016\\")
    println(io,cd(2),hbar^2,cd(2),dbar^2,cd(2))
    print(io,ind[1]," "^(4-l[1]),ind[2]," "^(4-l[2]),ind[3])
  end
  showrel(io,d,7)
end

function Base.show(io::IO,d::Diagram,::Val{:G31})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(2,ind[4]),kern(-2),raise(8,"\\diagup"),
    kern(-2),raise(16,ncnode(2,ind[1])),
    kern(-12.5),barr(19.5),ncnode(2,ind[2]),
    kern(-16.5),raise(12.3,larger("\\bigcirc",5)),
    kern(-9),barr(19.5),ncnode(2,ind[5]),
    kern(-22),raise(16,ncnode(2,ind[3])),
    kern(-1),raise(8,"\\diagdown"),"\$\$")
  else print(io,prod(map((i,j)->i*" "^(4-j),ind[[4,2,5]],l[[4,2,5]])))
    println(io,cd(2),hbar^2, cd(2),hbar^2, cd(2)," "^2,d.t)
    println(io," \\ /3\\ /")
    println(io," "^2,cd(2),hbar^2,cd(2))
    print(io," "^2,ind[1]," "^(4-l[1]),ind[3])
  end
  print(io,"  i.e. ",fromTeX(io,"\$A_5\$ on "))
  print(io,joindigits(d.t.indices[[1,4,2,5,3]])," plus")
  showrel(io,d,5,6)
end

function Base.show(io::IO,d::Diagram,::Val{:G32})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    print(io,"\$\$",ncnode(3,ind[1]))
    for i in 2:4 print(io,barr(10),ncnode(3,ind[i])) end
    print(io,"\$\$")
  else println(io,(cd(3)*hbar^2)^3,cd(3),d.t)
    print(io,prod(map((i,j)->i*" "^(4-j),ind,l)))
  end
end

trianglerel(a,b,c)=string(ncnode(2,a),barr(14),kern(-12),
    rlap(raise(7.5,"\\underleftarrow"*script(6))),
    rlap(kern(-6)*raise(9,"\\diagup")),
    rlap(kern(1)*raise(16.5,rcnode(2,c))),
    rlap(kern(8.6)*raise(9,"\\diagdown")),
    kern(11.5),"\n",rlap("\\hbox{\$"*ncnode(2,b)*"\$}"),
    kern(9))

function Base.show(io::IO,d::Diagram,::Val{:G33})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    println(io,"\$\$",ncnode(2,ind[1]),barr(10),
        trianglerel(ind[2],ind[3],ind[4]))
    println(io,barr(10),ncnode(2,ind[5]),"\$\$")
  else println(io," "^4, ind[3]," ",cd(2)," "^6,d.t)
    println(io," "^5,"/^\\")
    println(io,(cd(2)*hbar^2)^3,cd(2))
    print(io,prod(map((i,j)->i*" "^(4-j),ind[[1,2,4,5]],l[[1,2,4,5]])))
  end
  showrel(io,d,11)
end

function Base.show(io::IO,d::Diagram,::Val{:G34})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    println(io,"\$\$",ncnode(2,1),barr(10),trianglerel(ind[2],ind[3],ind[4]))
    println(io,barr(10),ncnode(2,5),barr(10),ncnode(2,6),"\$\$")
  else println(io," "^4, ind[3]," ",cd(2)," "^10,d.t)
    println(io," "^5,"/^\\")
    println(io,(cd(2)*hbar^2)^4,cd(2))
    print(io,prod(map((i,j)->i*" "^(4-j),ind[[1,2,4,5,6]],l[[1,2,4,5,6]])))
  end
  showrel(io,d,16)
end

function impdiagram(io::IO,D::Diagram,p,q,r,ind)
  d=div(p,q)
  j=chevieget(:imp, :BraidRelations)(p, q, r)
  if q==1
    if get(io,:TeX,false)
      print(io,"\$\$",ncnode(p,ind[1]),dbarr(10),ncnode(2,ind[2]))
      for i in 3:r print(io,barr(10),ncnode(2,i)) end
      print(io,"\$\$")
    else join(io,ind," "^3);println(io,"    ",D.t)
     print(io,cd(p))
     if length(ind)>1 print(io,dbar^2) end
     join(io,fill(cd(2),r-1),hbar^2)
    end
  elseif p==q
    if get(io,:TeX,false)
      print(io,"\$\$",p,kern(2),rule(-5,0.7,14.5),kern(-12),
        rlap(raise(12,lcnode(2,ind[1]))),
        raise(-12,lcnode(2,ind[2])),kern(-1.3),
        rlap(raise(7,"\\diagdown")),
        raise(-7,"\\diagup"),
        kern(-9.6),dbarr(8),ncnode(2,ind[3]))
        for i in 4:r print(io,barr(10),ncnode(2,i)) end
        print(io,"\$\$")
    else indent=p==3 ? "" : " "^length(string(p))
      println(io,indent,ind[1],"    ",D.t)
      println(io,indent,vbar,"\\")
      if p!=3 print(io,p) end
      println(io,vbar,dbar,join(ind[3:r],hbar^2))
      println(io,indent,vbar,"/")
      print(io,indent, ind[2])
    end
    showrel(io,D,2)
  elseif q==2
    if get(io,:TeX,false)
      print(io,"\$\$",lcnode(d,ind[1]),kern(-2.7),
            raise(-4,larger("\\bigcirc",5)),
        kern(-7),
        rlap(raise(12,rcnode(2,ind[2]))),
        raise(-12,rcnode(2,ind[3])))
      if r>3
        print(io,kern(-6.2),rlap(raise(7,"\\diagdown")),
        raise(-7,"\\diagup"),kern(-1.6),ncnode(2,ind[4]))
      end
      for i in 5:r print(io,barr(10),ncnode(2,ind[i])) end
      print(io,"\$\$")
    else
      indent=" "^(textwidth(cd(d))-textwidth(cd(2)))
      println(io,indent," "^3,ind[2]," ",cd(2),"   ",D.t)
      print(io,indent,"    /",vbar)
      if r>=3 print(io,"\\") end
      println(io,join(ind[4:r+1]," "^3))
      println(io,ind[1]," ",cd(d),"3",vbar," ",join(fill(cd(2),r-2),hbar^2))
      print(io,indent,"    \\",vbar)
      if r>=3 print(io,"/") end
      println(io)
      print(io,indent," "^3,ind[3]," ",cd(2))
    end
    showrel(io,D,1,2)
  else
    if get(io,:TeX,false)
      print(io,"\$\$",lcnode(d,ind[1]),kern(-1),
      rlap(kern(-2.5)*raise(-10,script(q+1))),
      kern(-2.4),
#     raise(-4,larger("(",7)),
      raise(0,"\\bigg ("),
      kern(-1),
      rlap(raise(10,rcnode(2,ind[2]))),
      raise(-10,rcnode(2,ind[3])))
      if r>2
        print(io,kern(-6.2),rlap(raise(7,"\\diagdown")),raise(-7,"\\diagup"),
        kern(-11.6),dbarr(10),ncnode(2,ind[4]))
      end
      for i in 5:length(ind) print(io,barr(10),ncnode(2,ind[i])) end
      print(io,"\$\$")
    else indent=" "^(textwidth(cd(d))-textwidth(cd(2)))
      println(io,indent," "^2,ind[2]," ",cd(2),"    ",D.t)
      print(io,indent,"   /",q+1)
      if r>=3 print(io,"\\  ") end
      println(io,join(ind[4:r+1]," "^3))
      print(io,ind[1]," ",cd(d),"  ")
      if r>=3 print(io,dbar^2) end
      println(io,join(fill(cd(2),r-2),hbar^2))
      print(io,indent,"   \\ ")
      if r>=3 print(io,"/") end
      println(io)
      print(io,indent," "^2,ind[3]," ",cd(2))
    end
    showrel(io,D,(1:min(3,r))...)
  end
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("Ã")})
  v=d.t.indices
  r=length(v)-1
  if r==1 println(io,d.t.series,"₁  ",v[1]," ∞ ",v[2])
  else 
    n=string(d.t.series,stringind(io,r),"   ")
    s=string(join(v[1:r],hbar^3))
    o=length(s)-4
    println(io," "^length(n)," ",hbar^div(o,2),v[r+1],hbar^div(o,2))
    println(io," "^length(n),"/"," "^o,"\\")
    println(io,n,s)
  end
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("B̃")})
  v=d.t.indices
  r=length(v)-1
  c=d.t.cartanType
  if c==1 || r==2
    println(io,"C̃",stringind(io,r)," ",v[1]," > ",join(v[2:end]," - ")," < ",v[r+1])
  else
    s=string("B̃",stringind(io,r)," ",v[1]," < ",join(v[2:end-2]," - "))
    println(io," "^(length(s)-2),v[r+1])
    println(io," "^(length(s)-2),"|")
    println(io,s," - ",v[r])
  end
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("D̃")})
  v=d.t.indices
  r=length(v)-1
  lr=length(string(r))
  println(io,d.t.series,stringind(io,r)," ",v[1]," "^(4*(r-3)-1),v[r+1])
  println(io," "^(lr+3),"\\"," "^(1+4*(r-4)),"/")
  print(io," "^(lr+4),v[3])
  for i in 4:r-1 print(io,hbar^3,v[i]) end
  println(io)
  println(io," "^(lr+3),"/"," "^(1+4*(r-4)),"\\")
  println(io," "^(lr+2),v[2]," "^(4*(r-3)-1),v[r])
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("Ẽ")})
  v=d.t.indices
  if length(v)==7
    print(io,"        ",v[7],"\n        |\n","        ",v[2],"\n        |\n",
          join(v[[1,3,4,5,6]],hbar^3),"    ",d.t)
  elseif length(v)==8
    println(io,"            ",v[2],"\n            |")
    print(io,join(v[[8,1,3,4,5,6,7]],hbar^3),"  ",d.t)
  elseif length(v)==9
    println("        ",v[2],"\n        |")
    print(io,join(v[[1,3,4,5,6,7,8,9]],hbar^3),"  ",d.t)
  end
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("F̃")})
  v=d.t.indices
  print(io,v[5],hbar^3,v[1],hbar^3,v[2],ldarrow(2),v[3],hbar^3,v[4]," "^3,d.t)
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("G̃")})
  v=d.t.indices
  print(io,v[3],hbar^2,v[1],tarrow(2),v[2],"   ",d.t)
end

end
