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
② ——② ——② ——②   423423==342342
1   2   4   5
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

nnode(i)="\\kern-2pt\\mathop\\bigcirc\\limits_{"*i*"}\\kern-2pt"
barr(n)=string("\\rule[2pt]{",n,"pt}{1pt}")
dbarr(n)=string("\\rlap{\\rule[3pt]{",n,"pt}{1pt}}\\rule[1pt]{",n,"pt}{1pt}")
tbarr(n)=string("\\rlap{\\rule[4pt]{",n,"pt}{1pt}}\\rlap{\\rule[2pt]{",n,
                "pt}{1pt}}\\rule[0pt]{",n,"pt}{1pt}")

vertbar(x,y)="\\kern1.5pt\\rlap{\\hbox{\$"*nnode(x)*"\$}}\n"*
   "\\rlap{\\kern2.6pt\\rule[6.6pt]{1pt}{10.8pt}}\n"*
   "\\rlap{\\raise19.4pt\\hbox{\$\\kern-2pt\\bigcirc\\scriptstyle "*y*"\$}}\n"*
   "\\kern8pt"

function getlind(d)
  t=d.t
  indices=t.indices
  if isnothing(indices) ind=fill("?",rank(t))
  else ind=repr.(indices)
  end
  length.(ind),ind
end

function Base.show(io::IO,d::Diagram,::Val{:A})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    println(io,"\$\$",nnode(ind[1]))
    for i in 2:length(l) println(io,barr(10),nnode(ind[i])) end
    println(io,"\$\$")
    return
  end
  join(io,node.*hbar.^l[1:end-1]);print(io,node);println(io," ",d.t)
  join(io,ind," ")
end

function Base.show(io::IO,d::Diagram,::Val{:B})
  l,ind=getlind(d)
  c=d.t.cartanType
  if get(io,:TeX,false)
    println(io,"\$\$",nnode(ind[1]))
    if c==2 println(io,"\\rlap{\\raise5pt\\hbox{\$\\leftarrow\$}}",dbarr(10))
    elseif c==1 println(io,"\\rlap{\\raise5pt\\hbox{\$\\rightarrow\$}}",dbarr(10))
    elseif c==root(2) println(io,dbarr(10))
    else error("not implemented")
    end
    println(io,nnode(ind[2]))
    for i in 3:length(l) println(io,barr(10),nnode(ind[i])) end
    println(io,"\$\$")
    return
  end
  print(io,node)
  if c==2 l1=max(l[1],2);print(io,rdarrow(l1))
  elseif c==1 l1=max(l[1],2);print(io,ldarrow(l1))
  elseif c==root(2) l1=max(l[1],2);print(io,dbar^l1)
  else xprint(io,"=",c,"="); l1=length(repr(c;context=rio()))+2
  end
  join(io,node.*hbar.^l[2:end-1]);println(io,node," ",d.t)
  print(io,rpad(ind[1],l1+1));join(io,ind[2:end]," ")
end

function Base.show(io::IO,d::Diagram,::Val{:D})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    println(io,"\$\$",nnode(ind[1]),barr(10))
    println(io,vertbar(ind[3],ind[2]))
    for i in 4:length(l) println(io,barr(10),nnode(ind[i])) end
    println(io,"\$\$")
    return
  end
  println(io," "^(l[1]+1),node," $(ind[2])")
  println(io," "^(l[1]+1),vbar)
  print(io,node,map(l->hbar^l*node,l[[1;3:end-1]])...)
  println(io," ",d.t)
  join(io,ind[[1;3:end]]," ")
end

function Base.show(io::IO,d::Diagram,::Val{:E})
  l,ind=getlind(d)
  if get(io,:TeX,false)
    println(io,"\$\$",nnode(ind[1]),barr(10))
    println(io,nnode(ind[3]),barr(10))
    println(io,vertbar(ind[4],ind[2]))
    for i in 5:length(l) println(io,barr(10),nnode(ind[i])) end
    println(io,"\$\$")
    return
  end
  println(io," "^(2+l[1]+l[3]),node," $(ind[2])")
  println(io," "^(2+l[1]+l[3]),vbar)
  print(io,node,map(l->hbar^l*node,l[[1;3:end-1]])...)
  println(io," ",d.t)
  join(io,ind[[1;3:end]]," ")
end

function Base.show(io::IO,d::Diagram,::Val{:F})
  l,ind=getlind(d)
  c=d.t.cartanType
  if get(io,:TeX,false)
     println(io,"\$\$",nnode(ind[1]))
     println(io,barr(10),nnode(ind[2]))
     if c==1 println(io,"\\rlap{\\raise5pt\\hbox{\$\\rightarrow\$}}",dbarr(10))
     elseif c==root(2) println(io,dbarr(10))
     else error("not implemented")
     end
     println(io,nnode(ind[3]))
     println(io,barr(10),nnode(ind[4]))
     println(io,"\$\$")
    return
  end
  print(io,node,hbar^l[1],node)
  if c==1 l1=max(l[2],2);print(io,ldarrow(l1))
  else l1=max(l[2],1);print(io,dbar^l1)
  end
  println(io,node,hbar^l[3],node," ",d.t)
  print(io,ind[1]," ",rpad(ind[2],l1+1),ind[3]," ",ind[4])
end

function Base.show(io::IO,d::Diagram,::Val{:G})
  l,ind=getlind(d)
  c=d.t.cartanType
  if get(io,:TeX,false)
    println(io,"\$\$",nnode(ind[1]))
    if c==3 println(io,"\\rlap{\\raise5pt\\hbox{\$\\leftarrow\$}}",tbarr(10))
    elseif c==1 println(io,"\\rlap{\\raise5pt\\hbox{\$\\rightarrow\$}}",tbarr(10))
    elseif c==root(3) println(io,dbarr(10))
    else error("not implemented")
    end
    println(io,nnode(ind[2]),"\$\$")
    return
  end
  if c==1 print(io,node,tarrow(max(l[1],2)))
  else print(io,node,tbar(max(l[1],2)))
  end
  println(io,node," ",d.t)
  print(io,ind[1]," "^max(3-l[1],1),ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:H})
  l,ind=getlind(d)
  println(io," "^l[1],"₅")
  println(io,map(i->node*hbar^l[i],1:length(l)-1)...,node," ",d.t)
  join(io,ind," ")
end

function Base.show(io::IO,d::Diagram,::Val{:I})
  l,ind=getlind(d)
  println(io," "^l[1],d.t.bond)
  println(io,node,hbar^l[1],node)
  join(io,ind," ")
end

function Base.show(io::IO,d::Diagram,::Val{:G4})
  l,ind=getlind(d)
  println(io,cd(3),hbar^2,cd(3),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G5})
  l,ind=getlind(d)
  println(io,cd(3),dbar^2,cd(3),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G6})
  l,ind=getlind(d)
  println(io,cd(2),tbar(2),cd(3),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G7})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^2,cd(3),ind[2]," ",d.t);
  println(io," ","/3\\")
  println(io,cd(2),hbar^2,cd(3))
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3),"=", f(2,3,1),"=",f(3,1,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G8})
  l,ind=getlind(d)
  println(io,cd(4),hbar^2,cd(4),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G9})
  l,ind=getlind(d)
  println(io,cd(2),tbar(2),cd(4),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G10})
  l,ind=getlind(d)
  println(io,cd(3),dbar^2,cd(4),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G11})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^2,cd(3),ind[2]," ",d.t);println(io," /3\\");
  println(io,cd(2),hbar^2,cd(4))
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3),"=",f(2,3,1),"=",f(3,1,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G12})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",cd(2),ind[2]," ",d.t);
  println(io," /4\\");println(io,cd(2),hbar^2,cd(2))
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3,1),"=",f(2,3,1,2),"=",f(3,1,2,3))
end

function Base.show(io::IO,d::Diagram,::Val{:G13})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",cd(2),ind[1]," ",d.t)
  println(io," / \\");println(io,cd(2),hbar^2,cd(2))
  print(io,ind[3]," "^3,ind[2]," ",f(2,3,1,2),"=",f(3,1,2,3)," ",f(1,2,3,1,2),"=")
  print(io,f(3,1,2,3,1))
end

function Base.show(io::IO,d::Diagram,::Val{:G14})
  l,ind=getlind(d)
  println(io,"  ₈");println(io,cd(2),hbar,cd(3),d.t)
  print(io,ind[1]," "^2,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G15})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",cd(2),ind[1]," ",d.t);
  println(io," /5");println(io,cd(2),ind[3]);
  println(io," \\");print(io,"  ",cd(3),ind[2]," ")
  print(io,f(1,2,3),"=",f(3,1,2)," ",f(2,3,1,2,1),"=",f(3,1,2,1,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G16})
  l,ind=getlind(d)
  println(io,cd(5),hbar^2,cd(5),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G17})
  l,ind=getlind(d)
  println(io,cd(2),tbar(2),cd(5),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G18})
  l,ind=getlind(d)
  println(io,cd(3),dbar^2,cd(5),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G19})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^2,cd(3),ind[2]," ",d.t);println(io," /3\\");
  println(io,cd(2),hbar^2,cd(5))
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3),"=",f(2,3,1),"=",f(3,1,2))
end

function Base.show(io::IO,d::Diagram,::Val{:G20})
  l,ind=getlind(d)
  println(io,"  ₅");println(io,cd(3),hbar,cd(3),d.t)
  print(io,ind[1]," "^2,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G21})
  l,ind=getlind(d)
  println(io,"  ₁₀");
  println(io,cd(2),hbar^2,cd(3),d.t);print(io,ind[1]," "^3,ind[2])
end

function Base.show(io::IO,d::Diagram,::Val{:G22})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",cd(2),ind[2]," ",d.t);println(io," /5\\")
  println(io,cd(2),hbar^2,cd(2))
  print(io,ind[1]," "^3,ind[3]," ",f(1,2,3,1,2),"=",f(2,3,1,2,3),"=",f(3,1,2,3,1))
end

function Base.show(io::IO,d::Diagram,::Val{:G24})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",cd(2),ind[1]," ",d.t)
  println(io," / \\")
  print(io,cd(2),dbar^2,cd(2))
  println(io,"  ",f(2,3,1,2,3,1,2,3,1),"==",f(3,2,3,1,2,3,1,2,3))
  print(io,ind[2]," "^3,ind[3])
end

function Base.show(io::IO,d::Diagram,::Val{:G25})
  l,ind=getlind(d)
  println(io,cd(3),hbar^2,cd(3),hbar^2,cd(3),d.t)
  print(io,ind[1]," "^3,ind[2]," "^3,ind[3])
end

function Base.show(io::IO,d::Diagram,::Val{:G26})
  l,ind=getlind(d)
  println(io,cd(2),dbar^2,cd(3),hbar^2,cd(3),d.t)
  print(io,ind[1]," "^3,ind[2]," "^3,ind[3])
end

function Base.show(io::IO,d::Diagram,::Val{:G27})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,"  ",cd(2),ind[1]," ",d.t)
  println(io, " / \\")
  print(io,cd(2),dbar^2,cd(2))
  println(io,"  ",f(3,2,3,1,2,3,1,2,3,1,2,3),"==",f(2,3,1,2,3,1,2,3,1,2,3,2))
  print(io,ind[2]," "^3,ind[3])
end

function Base.show(io::IO,d::Diagram,::Val{:G29})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^6,cd(2),f(4)," "^2,d.t)
  println(io,"     /\u2016\\")
  print(io,cd(2),hbar^2,cd(2),dbar^2,cd(2))
  println(io,"  ",f(4, 3, 2, 4, 3, 2),"==",f(3, 2, 4, 3, 2, 4))
  print(io,f(1)," "^3,f(2)," "^3,f(3))
end

function Base.show(io::IO,d::Diagram,::Val{:G31})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io,f(4)," "^3,f(2)," "^3,f(5))
  println(io,cd(2),hbar^2, cd(2),hbar^2, cd(2)," "^2,d.t)
  println(io," \\ /3\\ /")
  print(io," "^2,cd(2),hbar^2,cd(2),"    i.e. A₅ on ")
  println(io,f(1,4,2,5,3)," plus ",f(1,2,3),"==",f(2,3,1),"==",f(3,1,2))
  print(io," "^2,f(1)," "^3,f(3))
end

function Base.show(io::IO,d::Diagram,::Val{:G32})
  l,ind=getlind(d)
  println(io,(cd(3)*hbar^2)^3,cd(3),d.t)
  print(io,ind[1]," "^3,ind[2]," "^3,ind[3]," "^3,ind[4])
end

function Base.show(io::IO,d::Diagram,::Val{:G33})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^4, f(3)," ",cd(2)," "^6,d.t)
  println(io," "^5,"/^\\")
  println(io,(cd(2)*hbar^2)^3,cd(2),"  ",f(4,2,3,4,2,3),"==",f(3,4,2,3,4,2))
  print(io,f(1)," "^3,f(2)," "^3,f(4)," "^3,f(5))
end

function Base.show(io::IO,d::Diagram,::Val{:G34})
  l,ind=getlind(d);f(arg...)=join(ind[collect(arg)])
  println(io," "^4, f(3)," ",cd(2)," "^10,d.t)
  println(io," "^5,"/^\\")
  println(io,(cd(2)*hbar^2)^4,cd(2)," ",f(4,2,3,4,2,3),"==",f(3,4,2,3,4,2))
  print(io,f(1)," "^3,f(2)," "^3,f(4)," "^3,f(5)," "^3,f(6))
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("A",Char(0x00303))})
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

function Base.show(io::IO,d::Diagram,::Val{Symbol("D",Char(0x00303))})
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

function Base.show(io::IO,d::Diagram,::Val{Symbol("E",Char(0x00303))})
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

function Base.show(io::IO,d::Diagram,::Val{Symbol("F",Char(0x00303))})
  v=d.t.indices
  print(io,v[5],hbar^3,v[1],hbar^3,v[2],ldarrow(2),v[3],hbar^3,v[4]," "^3,d.t)
end

function Base.show(io::IO,d::Diagram,::Val{Symbol("G",Char(0x00303))})
  v=d.t.indices
  print(io,v[3],hbar^2,v[1],tarrow(2),v[2],"   ",d.t)
end

function impdiagram(io::IO,D::Diagram,p,q,r,ind)
  d=div(p,q)
  if q==1
    join(io,ind," "^3);println(io,"    ",D.t)
    print(io,cd(p))
    if length(ind)>1 print(io,dbar^2) end
    join(io,fill(cd(2),r-1),hbar^2)
  elseif p==q
    indent=p==3 ? "" : " "^length(string(p))
    println(io,indent,ind[1],"    ",D.t)
    println(io,indent,vbar,"\\")
    if p!=3 print(io,p) end
    println(io,vbar,dbar^2,join(ind[3:r],hbar^2))
    println(io,indent,vbar,"/")
    print(io,indent, ind[2])
  elseif q==2
    indent=" "^(textwidth(cd(d))-textwidth(cd(2)))
    println(io,indent," "^3,ind[2]," ",cd(2),"   ",D.t)
    print(io,indent,"    /",vbar)
    if r>=3 print(io,"\\") end
    println(io,join(ind[4:r+1]," "^3))
    println(io,ind[1]," ",cd(d),"3",vbar," ",join(fill(cd(2),r-2),hbar^2))
    print(io,indent,"    \\",vbar)
    if r>=3 print(io,"/") end
    println(io)
    print(io,indent," "^3,ind[3]," ",cd(2)," "^3,joindigits(ind[1:3]),"==",
      joindigits(ind[[2,3,1]]),"==",joindigits(ind[[3,1,2]]))
  else
    indent=" "^(textwidth(cd(d))-textwidth(cd(2)))
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
    j=chevieget(:imp, :BraidRelations)(p, q, r)
    for g in 1:min(3,r)
      print(io,"  ",joindigits(ind[j[g][1]]),"==",joindigits(ind[j[g][2]]))
    end
    println(io)
  end
end

end
