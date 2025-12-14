# cmplximp.jl Imprimitive groups  G_{p,q,r} Chevie Library
# (C) 1998-  Gunter Malle and  Jean Michel

chevieset(:imp,:SemisimpleRank,(p,q,r)->r)

chevieset(:imp,:ordergens,function(p,q,r)
  res=fill(2,r)
  if q==1 res[1]=p
  elseif q!=p pushfirst!(res,div(p,q))
  end
  res
end)

chevieset(:imp,:simpleroots,function(p,q,r)
  if q==1 
    roots=[map(i->i==1 ? 1 : 0,1:r)]
  else
    if q!=p roots=[map(i->i==1 ? Cyc(1) : Cyc(0),1:r)] end
    v=vcat([-E(p),1],fill(0,r-2))
    if r==2 && q>1 && isodd(q) v*=E(p) end
    if q==p roots=[v] else push!(roots, v) end
  end
  append!(roots,map(i->map(j->j==i-1 ? -1 : j==i ? 1 : 0,1:r),2:r))
  toM(roots)
end)

chevieset(:imp, :BraidRelations, function (p, q, r)
  function b(i,j,o)
    p(i,j)=map(k->i*mod(k,2)+j*mod(1-k,2),1:o)
    [p(i,j),p(j,i)]
  end
  res=Vector{Vector{Int}}[]
  if q==1
    if r>=2 push!(res, b(1, 2, p==1 ? 3 : 4)) end
    append!(res,map(i->b(i,i-1,3),3:r))
    for i in 3:r append!(res,map(j->b(i,j,2),1:i-2)) end
  elseif p==q
    push!(res,b(1,2,p))
    if r>=3
      append!(res,[[[1,2,3,1,2,3],[3,1,2,3,1,2]],b(1,3,3),b(2,3,3)])
    end
    append!(res,map(i->b(i,i-1,3),4:r))
    for i in 4:r append!(res,map(j->b(i,j,2),1:i-2)) end
  else
    push!(res,[[1,2,3],[2,3,1]])
    i=b(2,3,q-1)
    push!(res,[vcat([1,2],i[2]),vcat([3,1],i[1])])
    if r>=3
      if q!=2 push!(res,[[2,3,4,2,3,4],[4,2,3,4,2,3]]) end
      append!(res,[b(2,4,3),b(3,4,3),b(1,4,2)])
    end
    append!(res,map(i->b(i,i-1,3),5:r+1))
    for i in 5:r+1 append!(res,map(j->b(i,j,2),1:i-2)) end
  end
  res
end)

chevieset(:imp, :NrConjugacyClasses, function (p, q, r)
  if [q,r]==[2,2] div(p*(p+6),4)
  elseif q==1 npartition_tuples(r,p)
  else length(chevieget(:imp,:ClassInfo)(p,q,r)[:classtext])
  end
end)

chevieset(:imp,:CartanMat,function(p,q,r)
  rt=chevieget(:imp, :simpleroots)(p,q,r)
  rbar=conj(rt)
  e=1 .-E.(chevieget(:imp,:ordergens)(p,q,r))
  e=map(i->(e[i]*rbar[i,:])//sum(rbar[i,:].*rt[i,:]),1:length(e))
  [Cyc{Int}(sum(x.*y)) for x in e, y in eachrow(rt)]
end)

chevieset(:imp,:ReflectionDegrees,(p,q,r)->vcat(p*(1:r-1),div(r*p,q)))

chevieset(:imp, :ReflectionCoDegrees, function (p, q, r)
  res=collect(p*(0:r-1))
  if p==q && p>=2 && r>2 res[r]-=r end
  res
end)

chevieset(:imp,:AuName,function(p,q,r)
  if q!=1 || r!=1 error("pas prevu") end
  "Z"*stringind(rio(TeX=true),p)
end)

chevieset(:imp, :ParabolicRepresentatives, function (p, q, r, s)
  if q==1
    if p==1
      if s==0 
           [Int[]]
      else map(j->vcat(map(k->(sum(j[1:k-1])+k-1).+(1:j[k]),1:length(j))...),
                    vcat(map(i->partitions(s,i),1:r+1-s)...))
      end
    else vcat(map(i->map(j->vcat(1:i,i+1 .+j), 
        chevieget(:imp, :ParabolicRepresentatives)(1,1,r-i-1,s-i)),0:s)...)
    end
  elseif r==2
    if q==2 
      [[Int[]], [[1], [2], [3]], [1:3]][s + 1]
    elseif p==q
      [[Int[]],iseven(p) ?  [[1], [2]] : [[1]], [1:2]][s+1]
    end
  end
end)

chevieset(:imp, :CharParams, function (de, e, r)
  if e==1 return partition_tuples(r,de) end
  charparams=[]
  d = div(de, e)
  for t in partition_tuples(r,de)
    tt=map(i->circshift(t,i),(1:e).*d)
    if t==minimum(tt)
      s=findfirst(==(t),tt)
      if s==e push!(charparams, t)
      else t=t[1:s*d]
        s=div(e,s)
        append!(charparams,map(i->vcat(t,[s,i]),0:s-1))
      end
    end
  end
  charparams
end)

# see page 172 of Halverson and Ram,
# Murnaghan-Nakayama rules for characters of Iwahori-Hecke algebras of the 
# complex reflection groups G(r,p,n) Canad. J. Math. 50 (1998) 167--192
struct Hook
  area::Int       # area of the hook
  startpos::Int   # the position of the decreased βnumber in the βlist
  endpos::Int     # the position it will occupy after being decreased
  DC::Vector{Int} # axial positions of dull corners
  SC::Vector{Int} # axial positions of sharp corners
end

function Base.show(io::IO,h::Hook)
  print(io,"<",h.startpos,"=>",h.endpos,",area:",h.area,",DC=",h.DC,",SC=",h.SC,">")
end

"""
hooksβ(S,s)  the list of all hooks  of area less than or  equal to s in the
Young diagram of the partition with βlist S.
"""
function hooksβ(S,s)
  res=Hook[]
  e=length(S)
  if e==0 return res end
  j=e
  for i in S[e]-1:-1:0
    while j>0 && S[j]>i j-=1 end
    if j>0 && S[j]==i continue end
    for k in j+1:e
      if S[k]>i+s break end
      z=pushfirst!(S[j+1:k-1],i)
      zi=filter(i->z[i]-z[i-1]>1,2:length(z))
      push!(res,Hook(S[k]-i,k,j+1,map(i->z[i]-e,zi),
                push!(map(i->z[i-1]+1-e,zi),z[end]+1-e)))
    end
  end
  res
end

"""
stripsβ(S,s)  returns as a list of lists  of Hooks all broken border strips
(union  of disjoint  hooks) of  area less  than or  equal to s in the Young
diagram of the partition with βlist S.
"""
function stripsβ(S,s)
  res=[Hook[]]
  for hook in hooksβ(S,s)
    if s==hook.area push!(res,[hook])
    else
      j=hook.endpos-1-length(S)
      for hs in stripsβ(S[1:hook.endpos-1], s-hook.area)
        for h in hs
          h.SC.+=j
          h.DC.+=j
        end
        push!(hs, hook)
        push!(res, hs)
      end
    end
  end
  res
end

"the elementary symmetric function of degree t in the variables v"
function elementary_symmetric_function(t,v)
  sum(prod,combinations(v,t);init=zero(eltype(v)))
end

"the homogeneous symmetric function of degree t in the variables v"
homogeneous_symmetric_function(t,v)=sum(prod,multisets(v,t))

chevieset(:imp, :HeckeCharTable, function (p, q, r, para, rootpara)
  res=Dict{Symbol, Any}(:name=>string("H(G(", p, ",", q, ",", r, "))"),
  :degrees=> chevieget(:imp, :ReflectionDegrees)(p, q, r), :dim=>r)
  res[:identifier] = res[:name]
  res[:size] = prod(res[:degrees])
  res[:order] = res[:size]
  cl=chevieget(:imp, :ClassInfo)(p, q, r)
  if r==1
    res[:classes]=cl[:classes]
    res[:orders]=cl[:orders]
    res[:irreducibles]=[para[1][i]^j for i in 1:p, j in 0:p-1]
    res[:powermaps]=cl[:powermaps]
    merge!(res,chevieget(:imp,:CharInfo)(p,q,r))
  elseif q==1
# character table of the Hecke algebra of G(p,1,r) parameters v and [q,-1]
# according to [Halverson-Ram]"Characters of Iwahori-Hecke algebras of G(r,p,n)"
# Canadian Journal of math. 50 (1998) 167--192
    merge!(res,cl)
    T=@NamedTuple{area::Int64, cc::Int64, hooklength::Int64, 
         DC::Vector{Pair{Int64, Int64}}, SC::Vector{Pair{Int64, Int64}},
         stripped::Vector{Vector{Int64}}}
    StripsCache=Dict{Pair{Vector{Vector{Int}},Int},Vector{T}}()
"""
Strips(S,s)  returns the  list of  collections of  broken border  strips of
total  area equal to s coming from the various βlists in the symbol S. Each
collection  is  represented  by  a  named  tuple containing the statistical
information   about  it  necessary   to  compute  the   function  Delta  in
[Halverson-Ram] 2.17. These named tuples have the following fields:
   area       total area
   cc         number of connected components (hooks)
   hooklength sum of hooklengths of individual hooks (the # of rows -1 of the
                  hook, equal to startpos-endpos)
   DC         the list of all dull corners, represented by a pair:
                which βlist they come from, axial position
   SC         the list of all sharp corners, represented by a pair:
                which βlist they come from, axial position
   stripped   the symbol left after removing the strip collection
"""
    function Strips(S,s)local e, res, hs, ss, a
      get!(StripsCache,S=>s) do
      function strip(S,hs) # strip hooks hs from βlist S
        S=copy(S)
        for h in hs 
          beta=S[h.startpos]
          for i in h.startpos:-1:h.endpos+1 S[i]=S[i-1] end
          S[h.endpos]=beta-h.area
        end
        i=1
        while i<=length(S) && iszero(S[i]) i+=1 end
        if i>1 S=S[i:end].-(i-1) end
        return S
      end 
      e=length(S)
      if e==0
        if s==0 return [(area=0,cc=0,hooklength=0,DC=Pair{Int,Int}[],
                         SC=Pair{Int,Int}[],stripped=Vector{Int}[])]
        else return T[]
        end
      end
      res=T[]
      for h in stripsβ(S[e], s)
        hs=(area=reduce(+,map(x->x.area,h),init=0), 
          cc=length(h),
          hooklength=reduce(+,map(x->x.startpos-x.endpos,h),init=0),
          DC=[e=>y for x in h for y in x.DC],
          SC=[e=>y for x in h for y in x.SC],
          stripped=[strip(S[e],h)])
        for a in Strips(S[1:e-1],s-hs.area)
          push!(res,(area=a.area+hs.area,cc=a.cc+hs.cc,
            hooklength=a.hooklength+hs.hooklength,
            DC=vcat(a.DC,hs.DC), SC=vcat(a.SC,hs.SC),
            stripped=vcat(a.stripped,hs.stripped)))
        end
      end
      res
      end
    end
# the function Delta of Halverson-Ram 2.17, modified to take in account that
# our eigenvalues for T_2..T_r are (Q1,Q2) instead of (q,-q^-1)
    function Delta(k,hs,(Q1,Q2),v)local q
      delta=1
      if hs.cc>1
        if k==1 || iszero(Q1+Q2) return 0
        else delta*=(Q1+Q2)^(hs.cc-1)
        end
      end
      q=-Q1//Q2//1
      delta*=Q1^(hs.area-hs.cc)*(-q)^-hs.hooklength
      if k==0 return delta end
      ctSC=eltype(v)[v[x]*q^y for (x,y) in hs.SC].+0
      ctDC=eltype(v)[v[x]*q^y for (x,y) in hs.DC].+0
      delta*=reduce(*,ctSC;init=1)//reduce(*,ctDC;init=1)
      if k==1 return delta end
      delta*(-1)^(hs.cc-1)*sum(map(
          t->(-1)^t*elementary_symmetric_function(t,ctDC)* 
                    homogeneous_symmetric_function(k-t-hs.cc,ctSC), 
               0:min(length(ctDC),k-hs.cc));init=0)
    end
    chiCache=Dict{Pair{Vector{Vector{Int}},Vector{Vector{Int}}}, Any}()
    function entry(lambda,mu)
      get!(chiCache,lambda=>mu) do
        n=sum(sum,lambda)
        if iszero(n) return 1 end
        bp=maximum(x for S in lambda for x in S)
        i=findfirst(x->bp in x,lambda)
        # choice of bp and i corresponds to choice (sort) in classtext
        strips=Strips(mu, bp)
        if isempty(strips) return 0 end
        rest=copy(lambda)
        rest[i]=rest[i][2:end]
        f1=(-prod(para[2]))^((i-1)*(n-bp))
        f2=sum(strips) do x
          d=Delta(i-1,x,para[2],para[1])
          iszero(d) ? d : d*entry(rest, x.stripped)
        end
        f1*f2
      end
    end
    pts=partition_tuples(r,p)
    res[:irreducibles]=[entry(pt,x) for x in map(y->βset.(y),pts), pt in pts]
    merge!(res,chevieget(:imp,:CharInfo)(p,q,r))
  elseif r==2 && p!=q && iseven(q)
    res[:classes]=cl[:classes]
    res[:orders]=cl[:orders]
    if haskey(cl,:powermaps) res[:powermaps]=cl[:powermaps] end
    Z,X,Y=para
    e1=div(q,2)
    Z=root.(Z,e1)
    Z=vcat(map(i->Z.*E(e1,i),0:e1-1)...)
    # H is a subalgebra of H(G(p,2,2)) with parameters Z
    ci=chevieget(:imp, :CharInfo)(p, q, r)
    function entry2(char, class)
      char=ci[:malle][findfirst(==(char),ci[:charparams])]
      if char[1]==1
        w=[Z[char[4]],X[char[2]],Y[char[3]]]
        return prod(i->iszero(i) ? prod(w) : w[i],class;init=1)
      else
        if q==2 && p==4 && !any(ismissing,rootpara)
          w=-char[2]*E(4)*rootpara[1]*rootpara[2]*rootpara[3]*
            root(E(div(p,2),2-char[3]-char[4]))*E(p,char[3]+char[4]-2)
        else w=char[2]*root(prod(X)*prod(Y)*Z[char[3]]*Z[char[4]]*
          E(div(p,2),2-char[3]-char[4]))*E(p,char[3]+char[4]-2)
        end
        class=map(i->count(==(i),class), 0:3)
        if class[2]>0 char=sum(x->x^class[2],Z[char[[3,4]]])
        elseif class[3]>0 char=sum(X)
        elseif class[4]>0 char=sum(Y)
        else char = 2
        end
        return w^class[1]*char
      end
    end
    res[:irreducibles]=[entry2(char,c) for char in ci[:charparams],
                        c in cl[:classparams]]
    merge!(res,chevieget(:imp,:CharInfo)(p,q,r))
  else
    res[:classnames]=cl[:classnames]
    res[:orders]=cl[:orders]
    res[:centralizers]=cl[:centralizers]
    if haskey(cl,:powermaps) res[:powermaps]=cl[:powermaps] end
    res[:classes] = map(x->div(res[:size],x),res[:centralizers])
    reps=map(i->chevieget(:imp,:HeckeRepresentation)(p,q,r,para,rootpara,i),
             1:length(res[:classes]))
    res[:irreducibles]=toM(improve_type(map(reps)do r
      traces_words_mats(r,cl[:classtext])
    end))
  end
  res[:centralizers]=map(x->div(res[:size],x), res[:classes])
  res[:parameter]=para
  res[:irreducibles]=improve_type(res[:irreducibles]*one(prod(prod,para)))
  res
end)

const impchartableCache=Dict{Tuple{Int,Int,Int},Any}()

chevieset(:imp, :CharTable, function (p, q, r)
  get!(impchartableCache,(p,q,r))do
  oo=chevieget(:imp,:ordergens)(p,q,r)
  chevieget(:imp, :HeckeCharTable)(p,q,r,
    map(o->o==1 ? [1] : o==2 ? [1,-1] : E.(o,0:o-1),oo),fill(1,r+1))
  end
end)

chevieset(:imp,:b,function (p, q, r)
  if q==1 || p==q error("should not be called") end
end)

chevieset(:imp, :B, function (p, q, r)
  if q==1 || p==q error("should not be called") end
end)

chevieset(:imp, :FakeDegree, function (p, q, r, c, v)
  if q==1 fakedegree(Symbol_partition_tuple(c,1))(v)
  elseif q==p fakedegree(Symbol_partition_tuple(c,fill(0,p)))(v)
  else nothing
  end
end)

# JM 30/11/2000 we have to decide how to represent cuspidals of imprimitive
# groups -- the function below is an ad-hoc solution for now
# ImprimitiveCuspidalName(symbol) returns the TeX name
function ImprimitiveCuspidalName(S)
  r=rank(CharSymbol(S))
  d=length(S)
  s=joindigits(length.(S))
  if r==0 return "" end
  if sum(length,S)%d==1 # G(d,1,r)
    if r==1 d==3 ? "Z_3" : "Z_{$d}^{$s}"
    else "G_{$d,1,$r}^{$s}"
    end
  else # G(d,d,r)
    if r==2
      if d==4 "B_2"
      elseif d==6 
        p=Dict("212010"=>"-1","221001"=>"1",
               "211200"=>"\\zeta^2","220110"=>"\\zeta_3")
        "G_2[$(p[s])]"
      else "I_2($d)["*join(chevieget(:I,:SymbolToParameter)(S),",")*"]"
      end
    elseif r==3 && d==3 "G_{3,3,3}[\\zeta_3"* (s=="300" ? "" : "^2")*"]"
    elseif r==3 && d==4 "G_{4,4,3}["*"\\zeta_4"*(s=="3010" ? "" : "^3")*"]"
    else "G_{$d,$d,$r}^{$s}"
    end
  end
end

chevieset(:imp, :Invariants, function(p,q,r)
  map(1:r)do i
    if i==r function(arg...)prod(arg)^div(p,q) end
    else
      function(arg...) sum(a->prod(collect(arg[a]))^p,arrangements(1:r,i)) end
    end
  end
end)

chevieset(:imp, :CharInfo, function (de, e, r)
  res=Dict{Symbol,Any}(:charparams=>chevieget(:imp, :CharParams)(de, e, r))
  d=div(de,e)
  if e==1
    s=fill(0,d)
    s[1]=1
    res[:charSymbols]=Symbol_partition_tuple.(res[:charparams],Ref(s))
  else
    if d==1
      res[:charSymbols]=Symbol_partition_tuple.(res[:charparams],0)
      if [e,r]==[3, 3] res[:malleparam]=[[2,3,2],[2,3,3],[2,3,1],[3,4],[3,5],[1,9],
                                    [3,2],[3,1],[2,3,4],[1,0]]
      elseif [e,r]==[3,4] res[:malleparam]=[[12,6],[4,10],[6,8],[4,11],[1,18],[12,3],
         [6,5,2],[8,4],[8,5],[6,5,1],[3,9],[6,2],[2,6],[4,2],[4,1],[3,3],[1,0]]
# here the labeling is defined by phi_{6,5}' being the one which appears
# in the tensor square of the reflection representation phi_{4,1}
      elseif [e,r]==[3,5]
        res[:malleparam]=[[30,10],[20,12],[5,19],[10,14],[10,15],[5,20],[1,30],
          [30,7,1],[40,6],[30,7,2],[10,11],[15,10],[20,9],[20,8],[15,11],
          [10,12],[4,18],[30,4],[20,5],[10,8],[10,7],[20,6],[5,12],[20,3],
          [10,6],[15,4],[15,5],[10,5],[6,9],[10,3],[10,2],[5,6],[5,2],[5,1],
          [4,3],[1,0]]
# here the labeling is defined by phi_{30,7}'' being the one which appears
# in the tensor 4th power of the reflection representation phi_{5,1}
      elseif [e,r]==[4, 3] res[:malleparam]=[[6,3],[3,6,1],[3,5],[3,6,2],[1,12],
                                        [3,2,1],[3,2,2],[3,1],[2,4],[1,0]]
# here the labeling is defined by phi_{3,2}'' being the complex
# conjugate of phi_{3,1} and phi_{3,6}'' the complex conjugate of phi_{3,5}
      end
      if haskey(res,:malleparam) 
        res[:malle]=exceptioCharName.(res[:malleparam])
      end
    elseif iseven(e) && r==2
# .malle: indexing of chars as in Malle's paper on rank 2 cyclotomic algebras.
      res[:malle]=map(res[:charparams])do t
        local pos, de
        if t[end] isa Number
         t[end]==0 ? [1,2,1,findfirst(==([1]),t)] : [1,1,2,findfirst(==([1]),t)]
        else
          de=div(length(t),2)
          pos=filter(i->length(t[i])>0,1:length(t))
          if length(pos)==1
            t[pos[1]]==[2] ? [1,1,1,pos[1]-de] : [1,2,2,pos[1]-de]
          else pos[1]<=de ? [2,-1,pos[1],pos[2]-de] : [2,1,pos[2]-de,pos[1]-de]
          end
        end
      end
    end
  end
  if de==1 t=map(i->[vcat(i,fill(1,r-i))],r:-1:1) 
  else
    t=map(r:-1:0)do i
      v=map(x->Int[],1:de)
      if i>0 v[1]=[i] end
      v[2]=fill(1,r-i)
      if e>1 v=minimum(map(i->circshift(v,i*d),1:e)) end
      v
    end
  end
  res[:extRefl]=Int.(indexin(t,res[:charparams]))
  if e==1 || d==1
    res[:A]=degree_gendeg.(res[:charSymbols])
    res[:a]=valuation_gendeg.(res[:charSymbols])
    res[:B]=degree_feg.(res[:charSymbols])
    res[:b]=valuation_feg.(res[:charSymbols])
  end
  if e>1 && d>1
    galp=map(res[:charparams])do s
      s=copy(s)
      if s[end] isa Number
        t=div(length(s)-2,d)
        s[1:d:d*t-d+1]=circshift(s[1:d:d*t-d+1],1)
        s[1:end-2]=minimum(map(i->circshift(s[1:end-2],i*d),1:t))
        return s
      end
      s[1:d:d*e-d+1]=circshift(s[1:d:d*e-d+1],1)
      minimum(map(i->circshift(s,i*d),1:e))
    end
    res[:hgal]=Perm(string.(res[:charparams]),string.(galp)) #horrible hack
  end
  res[:charnames]=map(res[:charparams])do s
      if (s[end] isa Vector) && sum(sum,s)==1
        xrepr(E(length(s), findfirst(==([1]),s)-1),TeX=true)
      else
        string_partition_tuple(s;TeX=true)
      end
    end
  res
end)

chevieset(:imp, :classparams, function (p, q, r)
  if r==2 && p!=q && iseven(q)
    # classparams is the classtext where 0->z
    e1=div(q,2)
# if s,t,u generate G(p,2,2) then s':=s^e1,t,u generate G(p,q,2)
# z:=stu generates Z(G(p,2,2)) and z':=z^e1 generates Z(G(p,q,2))
# relations for G(p,2,2) are stu=tus=ust
# relations for G(p,q,2) are s'tu=tus' and [z',u]=1
    res=[]
    for i in 0:p-1, j in 0:div(p-i-1,2)
      if mod(j+i,e1)==0
        push!(res,vcat(fill(1,j),fill(0,i)))
      end
    end
    for j in 2:3, i in 0:e1:div(p,2)-e1
      push!(res, pushfirst!(fill(0,i),j))
    end
    return res
  end
  I=partition_tuples(r,p)
  if q==1 return I end
  res=[]
  for S in I
    if mod(sum(i->(i-1)*length(S[i]),1:p),q)!=0 continue end
    a=gcd(vcat(S...))
    a=gcd(a, q)
    for i in 1:p
      if !isempty(S[i]) a=gcd(a,i-1) end
    end
    if a==1 push!(res, S)
    else
      for j in 0:a-1 # number of pieces the class splits
        push!(res, vcat(S,div(p*j,a)))
      end
    end
  end
  return res
end)

using Primes: primes

chevieset(:imp,:pow,function(S,n)
  e=length(S)
  S1=map(x->[],1:e)
  for k in 1:e, l in S[k]
    g=gcd(n,l)
    append!(S1[1+mod(div(n,g)*(k-1),e)], fill(div(l,g),g))
  end
  for m in S1 sort!(m,rev=true) end
  S1
end)

chevieset(:imp, :ClassInfo, function (p, q, r)
  order(S)=lcm(map(((i,m),)->isempty(m) ? 1 : lcm(m)*div(p,gcd(i-1,p)),
                   enumerate(S)))
  centralizer(S)=p^sum(length,S)*prod(map(pp->prod(y->factorial(y[2])*y[1]^y[2],
                                        tally(pp);init=1),S))
  function classtext(S)
    s=vcat(map(i->map(t->[t,i-1],S[i]),1:p)...)
    sort!(s,by=a->(a[1],-a[2]))
    l=0
    w=Int[]
    for d in s
      append!(w,repeat(vcat(l+1:-1:2,1:l+1),d[2]))
      # non-reduced word because this is the one used by Halverson-Ram
      # for characters of the Hecke algebra (see HeckeCharTable).
      append!(w,l+2:l+d[1])
      l+=d[1]
    end
    w
  end
  if r==2 && p!=q && iseven(q)
    e1=div(q,2)
# if s,t,u generate G(p,2,2) then s':=s^e1,t,u generate G(p,q,2)
# z:=stu generates Z(G(p,2,2)) and z':=z^e1 generates Z(G(p,q,2))
# relations for G(p,2,2) are stu=tus=ust
# relations for G(p,q,2) are s'tu=tus' and [z',u]=1
    res=Dict{Symbol, Any}(:classtext=>[],:classparams=>[],:classnames=>[])
    for i in 0:p-1
      for j in 0:div(p-i-1,2)
        if mod(j+i,e1)==0
          push!(res[:classparams],vcat(fill(1,j),fill(0,i)))
          push!(res[:classtext],vcat(fill(1,div(j+i,e1)), repeat([2,3],i)))
          push!(res[:classnames],"1"^(div(j+i,e1)-div(i,e1))*"23"^(mod(i,e1))*
                "z"^div(i,e1))
        end
      end
    end
    for j in 2:3
     for i in 0:e1:div(p,2)-e1
        push!(res[:classparams], pushfirst!(fill(0,i),j))
        push!(res[:classtext], vcat([j], fill(1,div(i,e1)),repeat([2,3],i)))
        push!(res[:classnames], string(j,"z"^div(i,e1)))
      end
    end
    res[:orders] = map(res[:classparams])do c
      if length(c)>0 && c[1] in (2,3)
        lcm(2, div(p,gcd(count(iszero,c),p)))
      else
        lcm(div(p,gcd(count(iszero,c),p)),div(div(p,2),gcd(count(==(1),c), div(p, 2))))
      end
    end
    res[:classes]=map(res[:classparams])do c
      if length(c)>0 && c[1] in (2,3) div(p, 2)
      elseif 1 in c 2
      else 1
      end
    end
    return res
  elseif q==1
    cp=partition_tuples(r,p)
    res=Dict{Symbol,Any}(:classparams=>cp)
    res[:classtext]=classtext.(cp)
    res[:classnames]=chevieget(:imp,:ClassName).(cp)
    res[:orders]=order.(cp)
    res[:centralizers]=centralizer.(cp)
    res[:classes]=div.(p^r*factorial(r),res[:centralizers])
    l=maximum(res[:orders])
    res[:powermaps]=Vector{Any}(fill(nothing,l))
    for pw in primes(l) 
     res[:powermaps][pw]=map(x->findfirst(==(chevieget(:imp,:pow)(x,pw)),cp),cp)
  #   res[:powermaps][pw]=indexin(chevieget(:imp,:pow).(cp,pw),cp)
    end
    return res
  else
  # According  to Hugues  ``On  decompositions  in complex  imprimitive 
  # reflection groups'' Indagationes 88 (1985) 207--219:                
  #
  # Let l=(S_0,..,S_{p-1}) be  a p-partition of r specifying  a class C 
  # of G(p,1,r);  C is in G(p,q,r)  iff q divides sumᵢ i*|Sᵢ|;  C splits in 
  # d  classes of G(p,q,r) for  the largest  d|q which divides all parts  of 
  # all Sᵢ and  such that |Sᵢ|=0 if  d does not divide i;  
  # if w is in  C and t  is the first generator  of G(p,1,r) then t^i w t^-i 
  # for i in 0:d-1 are representatives of classes of G(p,q,r) which meet C.  
    function trans(w) local res
    # translate words  in G(p,1,r) into  words of G(p,q,r); use  that if
    # t,s2 (resp. s1,s2)  are the first 2 generators  of G(p,1,r) (resp.
    # G(p,p,r)) then  s1=s2^t thus  s2^(t^i)= (s1s2)^i  s2 [the  first 3
    # generators of G(p,q,r) are t^q,s1,s2].
      d=0
      res=Int[]
      word(l,i)=map(j->1+mod(j,2),i+l:-1:i+1)
      function add(a)
        l=length(res)
        if l>0 && res[end]==a pop!(res)
        elseif p==q && a in [1,2] && l>=q && res[l-q+1:l]==word(q,3-a)
          res=vcat(res[1:l-q], word(q-1,3-a))
        else
          push!(res, a)
        end
      end
      for l in w
        if l==1 d+=1
        elseif l!=2 add(l)
        else
          d=mod(d, p)
          if d==0 add(2)
          else for i in 1:p-d-1 add(1);add(2) end
            add(1)
          end
        end
      end
      d=mod(d, p)
      if mod(d, q)!=0 error()
      elseif d!=0 res=vcat(res.+1,fill(1,div(d,q)))
      elseif p!=q res.+=1
      end
      res
    end
    res=Dict{Symbol, Any}(:classtext=>Vector{Int}[],:classparams=>[],
      :orders=>Int[],:centralizers=>Int[])
    for S in partition_tuples(r,p)
      if mod(sum(j->(j-1)*length(S[j]),1:p),q)!=0 continue end
      a=gcd(vcat(S...))
      a=gcd(a, q)
      for j in 1:p
        if !isempty(S[j]) a=gcd(a,j-1) end
      end
      for j in 0:a-1 # number of pieces the class splits
        push!(res[:classtext],trans(vcat(fill(1,j),classtext(S),fill(1,p-j))))
        if a>1 push!(res[:classparams], vcat(S,[div(p,a)*j]))
        else push!(res[:classparams], S)
        end
        push!(res[:orders],order(S))
        push!(res[:centralizers], div(centralizer(S)*a, q))
      end
    end
    res[:classes]=div.(res[:centralizers][1], res[:centralizers])
    res[:classnames]=map(chevieget(:imp, :ClassName), res[:classparams])
    return res
  end
end)

chevieset(:imp, :ClassName, function(p)
  if all(x->x isa Vector, p)
    if sum(sum,p)==1 xrepr(E(length(p),findfirst(==([1]),p)-1);TeX=true)
    else string_partition_tuple(p)
    end
  elseif all(x->x isa Integer,p) joindigits(p)
  elseif all(i->p[i] isa Vector,1:length(p)-1) && p[end] isa Integer
    string_partition_tuple(vcat(p[1:end-1], [length(p)-1, p[end]]);TeX=true)
  end
end)

chevieset(:imp, :SchurModel, function (p, q, r, phi)
  if q==1 # cf. Chlouveraki, arxiv 1101.1465
    function GenHooks(l,m)
      if isempty(l) return [] end
      m=conjugate_partition(m)
      append!(m,fill(0,max(0,l[1]-length(m))))
      m=(m.+1).-(1:length(m))
      vcat(map(i->(l[i]-i).+m[1:l[i]],1:length(l))...)
    end
    res=Dict(:coeff=>(-1)^(r*(p-1)),:factor=>fill(0,p),:vcyc=>[])
    l=vcat(phi...)
    sort!(l)
    push!(res[:factor],sum(((1:length(l)).-length(l)).*l))
    for s in 1:p
      for t in 1:p
        for h in GenHooks(phi[s], phi[t])
          v=fill(0,p)
          if s!=t
            v[[s, t]] = [1, -1]
            push!(v, h)
            push!(res[:vcyc], [v, 1])
          else
            push!(v, 1)
            for d in divisors(h)
              if d>1 push!(res[:vcyc], [v, d]) end
            end
          end
        end
      end
    end
    return res
  elseif (q,r)==(2,2)
    ci = chevieget(:imp, :CharInfo)(p, q, r)
    phi = ci[:malle][findfirst(==(phi),ci[:charparams])]
    p2 = div(p, 2)
    if phi[1] == 1
      res=Dict(:coeff=>1,:factor=>fill(0,4+p2),:vcyc=>[])
      for l in [[1, -1, 0, 0], [0, 0, 1, -1]]
        append!(l,fill(0,p2))
        push!(res[:vcyc],[l,1])
      end
      for i in 2:p2, l in [[0, 0, 0, 0, 1], [1, -1, 1, -1, 1]]
        append!(l,fill(0,p2-1))
        l[4+i]=-1
        push!(res[:vcyc],[l,1])
      end
    else
      res=Dict(:coeff=>-2,:factor=>fill(0,4+p2),:vcyc=>[],
               :root=>vcat(fill(1//2,6),fill(0,p2-2)),
               :rootCoeff=>E(p2,2-phi[3]-phi[4]))
      for i in 3:p2, j in 1:2
        l=fill(0,4+p2);l[4+j]=1;l[4+i]=-1
        push!(res[:vcyc], [l,1])
      end
      for l in [[0,-1,0,-1,-1,0],[0,-1,-1,0,0,-1],[-1,0,-1,0,-1,0],
                [-1,0,0,-1,0,-1]]
        push!(res[:vcyc],[vcat(l,fill(0,max(0,p2-2)),1),1])
      end
    end
    return res
  else
    error("not implemented")
  end
end)

chevieset(:imp, :SchurData, function (p, q, r, phi)
  if (q,r)==(2,2)
    ci=chevieget(:imp, :CharInfo)(p, q, r)
    phi=ci[:malle][findfirst(==(phi),ci[:charparams])]
    if phi[1]==1
      res=Dict(:order=>[phi[2],3-phi[2],2+phi[3],5-phi[3],4+phi[4]])
      append!(res[:order], 4 .+setdiff(1:div(p,2),[phi[4]]))
    else
      res=Dict{Symbol,Any}(:order => [1, 2, 3, 4, 4 + phi[3], 4 + phi[4]])
      append!(res[:order],4 .+setdiff(1:div(p,2),phi[[3,4]]))
      res[:rootPower]=phi[2]*E(p,phi[3]+phi[4]-2)
    end
    return res
  else
    error("not implemented")
  end
end)

chevieset(:imp, :SchurElement, function (p, q, r, phi, para, rootpara)
  if r==1
    VcycSchurElement(vcat(para[1],0),chevieget(:imp,:SchurModel)(p,q,r,phi))
  elseif p==1
    VcycSchurElement([0,-para[1][1]//para[1][2]],
                     chevieget(:imp,:SchurModel)(p,q,r,phi))
  elseif q==1
    VcycSchurElement(vcat(para[1],[-para[2][1]//para[2][2]]),
                     chevieget(:imp,:SchurModel)(p,q,r,phi))
  elseif r==2 && iseven(q)
    e1=div(q,2)
    Z=root.(para[1],e1)
    Z=vcat(map(j->Z*E(e1,j),0:e1-1)...)
    VcycSchurElement(vcat(para[2],para[3],Z),
                     chevieget(:imp,:SchurModel)(p,2,2,phi),
                     chevieget(:imp,:SchurData)(p,2,2,phi))//e1
  elseif p==q
    if phi[end] isa Integer
      m=length(phi)-2
      phi=fullsymbol(phi)
    else
      m=p
    end
    chevieget(:imp,:SchurElement)(p,1,r,phi,vcat([E.(p,0:p-1)],para[2:end]),[])//m
  elseif para[2]==para[3]
    if phi[end] isa Integer
      m=length(phi)-2
      phi=fullsymbol(phi)
    else
      m=p
    end
    if para[1]==map(i->E(p//q,i-1),1:p//q)
      para=[E.(p,0:p-1),para[2]]
    else
      if q==2 
        if !ismissing(rootpara[1])
             Z=[rootpara[1],-E(4)*para[1][2]].//root(-para[1][2])
        else Z=[root(para[1][1]),root(para[1][2])]
        end
        Z=vcat(Z,-Z)
      else Z=[E(q,j)*root(i,q) for j in 0:q-1 for i in para[1]]
      end
      para=[Z,para[2]]
    end
    p//q*chevieget(:imp,:SchurElement)(p,1,r,phi,para,[])//m
  else
    InfoChevie("# SchurElements(H(G(",p,",",q,",",r,"),",para,") not implemented\n")
  end
end)

chevieset(:imp,:FactorizedSchurElement,function(p,q,r,phi,para,rootpara)
  if r==1
   VFactorSchurElement(vcat(para[1],[0]),chevieget(:imp,:SchurModel)(p,q,r,phi))
  elseif p==1
    VFactorSchurElement([0,-para[1][1]//para[1][2]],
                      chevieget(:imp,:SchurModel)(p, q, r, phi))
  elseif q==1
    VFactorSchurElement(vcat(para[1],[-para[2][1]//para[2][2]]),
                        chevieget(:imp,:SchurModel)(p, q, r, phi))
  elseif (q,r)==(2,2)
    VFactorSchurElement(vcat(para[[2,3,1]]...),
   chevieget(:imp,:SchurModel)(p,q,r,phi),chevieget(:imp,:SchurData)(p,q,r,phi))
  elseif p==q
    if phi[end] isa Integer
      m=length(phi)-2
      phi=fullsymbol(phi)
    else
      m=p
    end
    if para[1]!=para[2]
      InfoChevie("# FactorizedSchurElements(H(G($p,$q,$r),",
           HeckeAlgebras.simplify_para(para),") not implemented\n")
      return nothing
    end
    F=chevieget(:imp,:FactorizedSchurElement)(p,1,r,phi,vcat([E.(p,0:p-1)],
                                                           para[2:end]), [])
    F[:factor]//=m
    F
  elseif para[2]==para[3]
    if phi[end] isa Integer
      m=length(phi)-2
      phi=fullsymbol(phi)
    else
      m=p
    end
    if para[1]==map(i->E(p//q,i-1),1:p//q)
      para=[E.(i,0:p-1),para[2]]
    else
      para=[[E(q,j)*root(i,q) for j in 0:q-1 for i in para[1]],para[2]]
    end
    F=chevieget(:imp,:FactorizedSchurElement)(p,1,r,phi,para,[])
    F[:factor]=p//(q*m)*F[:factor]
    F
  else
    InfoChevie("# FactorizedSchurElements(H(G($p,$q,$r),",
           HeckeAlgebras.simplify_para(para),") not implemented\n")
  end
end)

# for g333 g334 g335 g443 D4=g224 we have integral models of representations
const rep335_1=let x=Pol()
  [[-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   E(3,2)-E(3,2)*x E(3,2)-E(3,2)*x 0 0 0 E(3,2) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 -E(3,2)*x+E(3,2)*x^2 0 E(3,2)-E(3,2)*x 0 0 0 E(3,2) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  (-E(3,2)-root(-3)*x)+E(3)*x^2 (E(3)-E(3)*x)+E(3)*x^2 0 0 0 E(3)-E(3)*x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 E(3,2)*x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 E(3) -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 (E(3,2)*x+root(-3)*x^2)-E(3)*x^3 0 (E(3)-E(3)*x)+E(3)*x^2 0 0 0 E(3)-E(3)*x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 E(3) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 E(3,2)-E(3,2)*x 0 0 0 0 0 0 0 E(3,2)-E(3,2)*x 0 0 0 0 E(3,2) 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 E(3,2) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 E(3,2)*x 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 E(3)*x 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 (-E(3,2)-root(-3)*x)+E(3)*x^2 0 0 0 0 0 0 0 (E(3)-E(3)*x)+E(3)*x^2 0 0 0 0 E(3)-E(3)*x 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3)*x 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3,2)*x 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3,2) 0 -1+x 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3) 0 -1+x 0 0 0 0 0 0 0 0 0 0;
  0 0 ((E(3,2)-E(3,2)*x)-E(3,2)*x^2)+E(3,2)*x^3 0 (E(3)-root(-3)*x)-E(3,2)*x^2 0 0 0 -E(3,2)+E(3,2)*x 0 -E(3)+root(-3)*x+E(3,2)*x^2 0 0 0 0 E(3,2)-E(3,2)*x 0 0 0 0 -1 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3,2) 0 0 0 0 0 0;
  E(3)-E(3)*x 0 0 (2*E(3)-E(3)*x^-1)-E(3)*x 0 0 0 0 0 0 (-E(3)+2*E(3)*x)-E(3)*x^2 0 0 0 0 -E(3)+E(3)*x 0 0 0 E(3)-E(3)*x 0 0 E(3)-E(3)*x 0 E(3) 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3)*x 0 -1+x 0 0 0 0 0 0;
  -E(3)+root(-3)*x+E(3,2)*x^2 0 (x-2*x^2)+x^3 E(3,2)-2E(3)+E(3)*x^-1+(E(3)-2E(3,2))*x+E(3,2)*x^2 0 0 0 0 0 0 ((E(3)-2*E(3)*x)+2*E(3)*x^2)-E(3)*x^3 0 0 0 0 (E(3)-2*E(3)*x)+E(3)*x^2 0 -E(3)+E(3)*x 0 (-E(3)+2*E(3)*x)-E(3)*x^2 0 0 (E(3,2)-E(3,2)*x)+E(3,2)*x^2 0 E(3,2)-E(3,2)*x 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 -1+x 0 -E(3,2)+E(3,2)*x 0 0 0 -1+x^-1 0 E(3,2)-E(3,2)*x 0 0 0 0 0 -1 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3,2)*x 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3) -1+x 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+x E(3)*x;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3,2) 0],
 [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 x 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  1-x 0 0 0 0 0 -1+x x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  -1+x^-1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 x 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  1-2x+x^2 1-2x+x^2 1-2x+x^2 1-x 2-x^-1-x 1-x 0 0 1-x^-1 -1+x 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1-x 0 0 0 0 0 0 0 -1+x 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  1-2x+x^2 -x+x^2 1-2x+x^2 1-x 1-x 0 0 0 0 x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1-x 0 0 0 0 0 0 0 x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  2-x^-1-x 1-x 0 0 0 1-x -1+x -1+x 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 x-x^2 -1+2x-x^2 0 0 0 x-x^2 0 0 0 0 0 -1+x 0 0 -1+x 0 x 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+x 0 x 0 0 0 0 0 0 0 0 0 0;
  2-x^-1-x 0 1-x 0 0 0 -1+x 0 0 0 0 1-x 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1-x^-1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0;
  -1+3x-3*x^2+x^3 -1+3x-3*x^2+x^3 x-2*x^2+x^3 2-x^-1-x -2+x^-1+2x-x^2 -1+2x-x^2 0 1-2x+x^2 -2+x^-1+x 0 0 -1+x -1+x 2-x^-1-x 1-x 0 0 0 1-x 0 0 -1+x 0 1 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
  -3+x^-1+3x-x^2 -1+2x-x^2 0 1-2x+x^2 0 -1+2x-x^2 1-2x+x^2 1-2x+x^2 -1+x x-x^2 0 -1+2x-x^2 0 1-x 1-x 0 -1+x 0 0 0 0 x 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 -1+x 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0;
  ((1-3x)+3*x^2)-x^3 ((1-2x)+2*x^2)-x^3 ((1-3x)+3*x^2)-x^3 (1-2x)+x^2 ((3-x^-1)-3x)+x^2 (1-2x)+x^2 0 0 (2-x^-1)-x 0 (-1+2x)-x^2 0 1-x 0 0 0 0 -1+x 0 0 1-x 0 (1-2x)+x^2 0 1-x 0 -1+x x 0 0;
  0 0 0 (2-x^-1)-x 0 -1+x 0 0 0 -1+x 0 0 0 0 0 -2+x^-1+x 0 (2-x^-1)-x 0 1-x -1+x^-1 0 -1+x 0 0 0 1 0 0 0;
  (-1+2x)-x^2 0 0 0 0 0 (1-2x)+x^2 -x+x^2 0 0 0 0 0 0 0 1-x -1+x 0 0 0 0 0 1-x 0 0 1-x 0 0 0 x;
  0 0 (1-2x)+x^2 ((3-x^-1)-3x)+x^2 0 0 1-x (1-2x)+x^2 0 0 -1+x 0 0 (2-x^-1)-x 0 0 0 0 1-x 0 0 0 (2-x^-1)-x 0 1-x^-1 1-x 0 0 1 -1+x],
 [0 -x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  -1 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 x 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 x 0 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -1+x 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 -1+x 0 0 0 0 0 0;
  0 x-x^2 0 0 0 0 0 0 0 1-2x+x^2 -1+2x-x^2 0 0 0 0 0 0 -1+x 0 0 1-x 0 1-2x+x^2 0 0 0 -1+x x 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+x 0 0 0 0 0 1-x 0 0 0 0 0 0 x;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 -1+x 0 0 0;
  1-x 0 0 2-x^-1-x 0 0 0 0 0 0 -1+2x-x^2 0 0 0 0 -1+x 0 0 0 1-x 0 0 1-x 0 1 0 0 -1+x 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
  0 0 0 0 0 0 0 0 0 0 0 1-x 0 0 0 0 0 1-x^-1 0 0 0 0 0 0 0 1 0 0 0 -1+x],
 [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 x -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  1-x 0 0 0 0 0 -1+x x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 x 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  1-x 1-x 0 0 0 1 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 x-x^2 0 0 0 0 0 0 0 -x+x^2 0 x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 1 0 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 0 0 0 0 0;
  0 0 -x+x^2 0 1-x 0 0 0 1 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+x 0 x 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0;
  0 0 1-x 0 0 0 0 0 0 0 1-x 0 0 0 0 1 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -1+x 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x -1+x 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 -1+x],
 [0 0 -x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  -1 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 x 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 1-x 0 0 0 0 0 0 0 1-x 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 x 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 1 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 1-x 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
  1-x 0 0 0 0 0 -1+x x 0 0 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0 0 0;
  1-x 0 0 2-x^-1-x 0 0 0 0 0 0 -1+2x-x^2 0 0 0 0 -1+x 0 0 0 1-x 0 0 1-x 0 1 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0 0;
  0 -x+x^2 0 0 0 0 1-x 0 0 0 0 0 0 0 x 0 0 0 0 0 -1+x 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 -1+x 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0;
  0 0 x-x^2 -1+2x-x^2 0 0 0 x-x^2 0 0 0 0 0 -1+x 0 0 -1+x 0 x 0 0 0 0 0 -1+x 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x 0 0 0 0 -1+x 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 -1+x 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]]
end

const rep335_9=let q=Pol();j=E(3);j2=E(3,2)
[[-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  -q q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 -q 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  j2-2j2*q+j2*q^2 j+(j2-j)*q-j2*q^2 j2+(-j2+j)*q-j*q^2 0 0 0 0 j-j*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (j-j*q)+j*q^2 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2 0 0 0 0 0;
  -1+2q-q^2 1-2q+q^2 -1+2q-q^2 -1+q 0 -1+q 0 1-q 0 j-j*q 0 0 0 j*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1-q 0 1-2q+q^2 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 -1 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  -3+q^-1+3q-q^2 3-q^-1-3q+q^2 -3+q^-1+3q-q^2 -2+q^-1+q 0 j2+2j+q^-1-j*q 0 2-q^-1-q 0 -j2+j2*q^-1+j2*q 0 0 0 j2-j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j2-2j+j*q^-1-q 0 2-q^-1-2q+q^2 0 0;
  -j+j*q j-j*q 0 -2j+j*q^-1+j*q 0 0 0 0 0 0 0 -j+j*q 0 0 j-j*q 0 0 0 j2 0 0 0 j-j*q 0 0 0 0 0 0 0 0 0 0 0 0 (2j-j*q^-1)-j*q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 -j2+j2*q 0 j2-j2*q 0 0 j2-j2*q 0 0 0 0 0 0 j2*q 0 0 0 0 0 j2 0 0 0 0 0 0 -j2+j2*q 0 0 0 0;
  0 0 0 0 0 j2-j2*q^-1 0 0 0 0 0 -j2+j2*q 0 0 0 0 j2-j2*q 0 0 j2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2-j2*q 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  (1-2q)+q^2 (j+(-j2-2j)*q)-q^2 0 ((3-q^-1)-3q)+q^2 0 0 0 0 0 0 -1+q (1+(j2+2j)*q)-j*q^2 0 0 (j-j*q)+j*q^2 0 0 0 j2-j2*q 0 0 0 (-1+2q)-q^2 0 0 0 0 0 0 0 0 0 0 0 0 (((j2+3j)-j*q^-1)+(-2j2-3j)*q)-q^2 0 0 0 0;
  0 0 0 0 0 -2j2+j2*q^-1+j2*q 0 0 0 0 0 (j2+(-j2+j)*q)-j*q^2 0 0 0 0 (j-j*q)+j*q^2 0 0 j-j*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (j+(j2-j)*q)-j2*q^2 0 0 0;
  0 0 -j+j*q 0 j-j*q 0 -j+j*q 0 0 0 0 0 0 0 0 0 (j-2j*q)+j*q^2 0 0 j-j*q j-j*q 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  -1+q 0 0 -2+q^-1+q 0 0 0 0 0 0 2-q^-1-q -1+q -2+q^-1+q 0 0 -2+q^-1+q 0 1-q 0 0 0 j-j*q 1-q 0 0 j 0 0 -1+q^-1 0 0 0 0 0 0 (2-q^-1)-q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 1-q^-1 0 0 0 0 0 -1+q 0 0 0 0 0 1-q 0 1 0 0 0 j2-j2*q 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
  0 0 -2j2+j2*q^-1+j2*q 0 (2j2+j+q^-1)-j2*q 0 ((-2j2-j)+j2*q^-1)-q 0 0 0 0 (j-2j*q)+j*q^2 0 0 0 0 ((2j2-j2*q^-1)-2j2*q)+j2*q^2 0 0 (2j2-j2*q^-1)-j2*q -1+q^-1+q 0 0 0 j2-j2*q 0 0 0 0 0 0 0 0 0 0 0 -2j+j*q^-1+j*q 0 0 0;
  (1+(j2+2j)*q)-j*q^2 0 0 (((-2j2-3j)-q^-1)+(j2+3j)*q)-j*q^2 0 0 0 0 0 0 (-3+q^-1+3q)-q^2 (-j+2j*q)-j*q^2 ((3-q^-1)-3q)+q^2 0 0 ((2-q^-1)-2q)+q^2 0 (j+(-j2-2j)*q)-q^2 0 0 0 (j2-j2*q)+j2*q^2 (-1+2q)-q^2 0 0 j2-j2*q 0 0 (2-q^-1)-q 0 0 0 0 0 0 (-3+q^-1+3q)-q^2 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j2 0 0 0 0 0 0;
  0 0 -1+q 0 0 0 -1+q 0 0 0 0 0 1-q 0 0 0 0 0 0 1-q 0 0 0 (j2-2j2*q)+j2*q^2 0 0 0 j-j*q 0 q 1-q 0 0 0 0 0 0 0 0 0;
  0 0 0 -q+q^2 0 0 1-q 0 0 0 (-j+j*q)-j*q^2 0 (j+(j2-j)*q)-j2*q^2 0 0 (j-j*q)+j*q^2 0 0 0 0 0 0 0 0 0 0 0 0 j-j*q 0 0 0 0 0 0 -j+(-j2+j)*q+j2*q^2 0 0 0 0;
  0 0 ((j2-j)+j*q^-1)-j2*q 0 0 -2+q^-1+q (2j2-j2*q^-1)-j2*q 0 0 0 0 ((3-q^-1)-3q)+q^2 (-j2+j+j2*q^-1)-j*q 0 0 0 0 -2+q^-1+q 0 (2j-j*q^-1)-j*q 0 0 0 ((2-q^-1)-2q)+q^2 0 0 0 -1+q^-1+q 0 j2-j2*q (2j-j*q^-1)-j*q 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 -1+q 0 0 0 0 0 1-q 0 0 0 0 (-j+j*q)-j*q^2 (j+(j2-j)*q)-j2*q^2 0 0 0 0 0 (1-q)+q^2 0 0 0 0 0 0 j-j*q 0 0 0 0 0 -j+(-j2+j)*q+j2*q^2 0 0 0;
  (-1+2q)-q^2 (((2j2-j)-j2*q^-1)+(-j2+2j)*q)-j*q^2 j+(j2-2j)*q+(-2j2+j)*q^2+j2*q^3 ((((7j2+j2*q^-2)-4j2*q^-1)-7j2*q)+4j2*q^2)-j2*q^3 0 (((-j2+2j)-j*q^-1)+(2j2-j)*q)-j2*q^2 ((3j2-j2*q^-1)-3j2*q)+j2*q^2 ((2j2-j)-j2*q^-1)+(-2j2+j)*q+j2*q^2 0 ((2-q^-1)-2q)+q^2 (((-5j2+2j)-j2*q^-2)+(3j2-j)*q^-1+(5j2-j)*q)-2j2*q^2 (((3j2-j)-j2*q^-1)+(-4j2+2j)*q+(3j2-j)*q^2)-j2*q^3 (j2-2j)+j2*q^-2+(-2j2+j)*q^-1+j*q (-1+2q)-q^2 ((-j2+2j2*q)-2j2*q^2)+j2*q^3 (2j2-2j)+j2*q^-2+(-2j2+j)*q^-1+(-2j2+j)*q+j2*q^2 0 -2j2+j+j2*q^-1+(j2-2j)*q+j*q^2 (-1+2q)-q^2 (-2j2+j+j2*q^-1+(2j2-j)*q)-j2*q^2 0 (-2+q^-1+2q)-q^2 -3j2+j+j2*q^-1+(4j2-2j)*q+(-3j2+j)*q^2+j2*q^3 (-j2-3j)+j*q^-1+(2j2+4j)*q+(-j2-3j)*q^2+j*q^3 0 -2+q^-1+q (j-j*q)+j*q^2 (-2+q^-1+2q)-q^2 ((2j2-j)+j2*q^-2+(-2j2+j)*q^-1)-j2*q (j2-2j2*q)+j2*q^2 (-2j2+j+j2*q^-1+(2j2-j)*q)-j2*q^2 j-j*q (1-q)+q^2 0 0 ((((-5j2+2j)-j2*q^-2)+(3j2-j)*q^-1+(6j2-j)*q)-4j2*q^2)+j2*q^3 0 (((2j2-j)-j2*q^-1)+(-2j2+2j)*q+(2j2-j)*q^2)-j2*q^3 0 0;
  0 -2+q^-1+q (1-2q)+q^2 ((-4-q^-2)+3*q^-1+3q)-q^2 0 (2-q^-1)-q -2+q^-1+q -2+q^-1+q 0 -2j+j*q^-1+j*q (((-3j2-4j)+q^-2)-3*q^-1)+(j2+2j)*q (-3+q^-1+3q)-q^2 (-1-q^-2)+2*q^-1 j-j*q (1-2q)+q^2 ((3j2+2j)-q^-2)+(-3j2-2j)*q^-1+q 0 (2-q^-1)-q j-j*q (2-q^-1)-q 0 (2j-j*q^-1)-j*q ((3-q^-1)-3q)+q^2 ((2j2-j)-j2*q^-1)+(-2j2+j)*q+j2*q^2 0 j-j*q^-1 1-q (2j-j*q^-1)-j*q (-1-q^-2)+2*q^-1 -1+q (2-q^-1)-q 1 j2-j2*q 1 0 (((4+q^-2)-3*q^-1)-3q)+q^2 0 (3j2+2j+q^-1+(-3j2-2j)*q)-q^2 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j*q 0 0 0 0 0 0 -1+q 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 j*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0;
  0 0 0 -q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0;
  0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0;
  -j2+j2*q j2-j2*q -j2+j2*q 0 0 0 0 j2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2-j2*q 0 0;
  0 0 0 -1+q 0 -1+q 0 0 0 0 0 0 0 j*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1-q 0 1-q 0 j-j*q j*q;
  3-q^-1-3q+q^2 -3+q^-1+3q-q^2 3-q^-1-3q+q^2 -j2+j2*q 0 -j2+j2*q^-1 0 -2+q^-1+q -1+q j2-j2*q^-1-j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2-q^-1-q j2+2j-j*q^-1+q -j2-2j+j*q^-1-q -2+q^-1+2q-q^2 -j2+j2*q^-1+j2*q j2-j2*q],
 [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  -q q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 -q 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0;
  1-q 0 1-q 0 1-q^-1 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1-q 0 0 1-q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 -1 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 1 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1-q 0 0 1-q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 1 0 0 0 0 0 0 -1+q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 -2+q^-1+q 0 0 -1+q^-1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1-q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  -2+q^-1+q (2-q^-1)-q (-1+2q)-q^2 ((3-q^-1)-3q)+q^2 -2+q^-1+q 0 (2-q^-1)-q (2-q^-1)-q -j+j*q (2j-j*q^-1)-j*q j-j*q (1-2q)+q^2 0 0 (-1+2q)-q^2 0 (3j2+2j+q^-1+(-3j2-2j)*q)-q^2 0 -j+j*q -2+q^-1+q -2+q^-1+q 0 (-1+2q)-q^2 0 j2-j2*q 0 -1+q 0 0 0 0 0 0 -1 -2+q^-1+q (-1+2q)-q^2 -2+q^-1+q ((-3j2-2j)-q^-1)+(3j2+2j)*q+q^2 -2j+j*q^-1+j*q j-j*q;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 -q 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 -1+q 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 -1+q j2*q 0 0 0 0 0 0 0;
  (2j-j*q^-1)-j*q -2j+j*q^-1+j*q (j-2j*q)+j*q^2 (-3j+j*q^-1+3j*q)-j*q^2 (2j-j*q^-1)-j*q 0 -2j+j*q^-1+j*q -2j+j*q^-1+j*q j2-j2*q -2j2+j2*q^-1+j2*q -j2+j2*q (-j+2j*q)-j*q^2 0 0 (j-2j*q)+j*q^2 0 ((j2+3j)-j*q^-1)+(-j2-3j)*q+j*q^2 0 j2-j2*q (2j-j*q^-1)-j*q (2j-j*q^-1)-j*q 0 (j-2j*q)+j*q^2 0 -1+q 0 j-j*q 0 0 0 0 j 0 j (2j-j*q^-1)-j*q (j-2j*q)+j*q^2 (2j-j*q^-1)-j*q ((-j2-3j)+j*q^-1+(j2+3j)*q)-j*q^2 (2j2-j2*q^-1)-j2*q -j2+j2*q;
  (-j+2j*q)-j*q^2 (-1+2q)-q^2 (-j+2j*q)-j*q^2 (-j2+2j2*q)-j2*q^2 (2j2-j2*q^-1)-j2*q (2-q^-1)-q (-j2+2j2*q)-j2*q^2 j-j*q (-j2+2j2*q)-j2*q^2 -j+j*q (j2-2j2*q)+j2*q^2 (-1+2q)-q^2 0 0 -q+q^2 0 (1-2q)+q^2 0 0 -j+j*q 1-q 0 j*q-j*q^2 0 0 0 -q 0 0 0 0 0 0 0 -j+j*q (j2-2j2*q)+j2*q^2 (1-2q)+q^2 (-1+2q)-q^2 j-j*q 0;
  1-q 0 0 0 1-q 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 -q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0;
  0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0;
  0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q;
  0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1+q],
 [0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  1 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0;
  0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 j*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 1 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 j2 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 1 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 1 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2 0 0 0 0 0 0 0 0 0 0 0 0;
  -2+q^-1+q (2-q^-1)-q (-1+2q)-q^2 ((3-q^-1)-3q)+q^2 -2+q^-1+q 0 (2-q^-1)-q (2-q^-1)-q -j+j*q (2j-j*q^-1)-j*q j-j*q (1-2q)+q^2 0 0 (-1+2q)-q^2 0 (3j2+2j+q^-1+(-3j2-2j)*q)-q^2 0 -j+j*q -2+q^-1+q -2+q^-1+q 0 (-1+2q)-q^2 0 j2-j2*q 0 -1+q 0 0 0 0 0 0 -1 -2+q^-1+q (-1+2q)-q^2 -2+q^-1+q ((-3j2-2j)-q^-1)+(3j2+2j)*q+q^2 -2j+j*q^-1+j*q j-j*q;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1-q 0 0 1-q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 q 0 0 0 0 0 0 -1+q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 -1+q 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2 0 -1+q 0 0 0 0 0 0 0;
  (-3j2+j2*q^-1+(4j2+j)*q+(-3j2-2j)*q^2)-q^3 (3j2-j2*q^-1)+(-4j2-j)*q+(3j2+2j)*q^2+q^3 ((-4j2-2j)-q^-1)+(6j2+2j)*q+(-4j2-j)*q^2+j2*q^3 (j2-3j)+j*q^-1+(-2j2+4j)*q+(j2-3j)*q^2+j*q^3 ((2j-j*q^-1)-2j*q)+j*q^2 (2-q^-1)-q (-2j+j*q^-1+2j*q)-j*q^2 ((3j2-j2*q^-1)-3j2*q)+j2*q^2 0 ((3-q^-1)-3q)+q^2 0 j2+(-2j2+j)*q+(j2-2j)*q^2+j*q^3 0 j-j*q (-j2+(2j2-j)*q+(-j2+2j)*q^2)-j*q^3 0 (((3j-j*q^-1)-4j*q)+3j*q^2)-j*q^3 0 -1+(-2j2-j)*q+j2*q^2 ((2j-j*q^-1)-2j*q)+j*q^2 ((2j-j*q^-1)-2j*q)+j*q^2 0 (-j2+(2j2-j)*q+(-j2+2j)*q^2)-j*q^3 0 (-1+q)-q^2 0 0 0 0 0 0 0 0 j-j*q (((-j2+2j)-j*q^-1)+(2j2-j)*q)-j2*q^2 ((j-3j*q)+3j*q^2)-j*q^3 (((-j2+2j)-j*q^-1)+(2j2-j)*q)-j2*q^2 (((4j2-j2*q^-1)-6j2*q)+4j2*q^2)-j2*q^3 (((3j2+j)-j2*q^-1)+(-3j2-2j)*q)-q^2 -j2+(2j2+j)*q+q^2;
  q-q^2 -q+q^2 q-q^2 0 0 0 0 -q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 -q+q^2 0 0;
  0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0;
  0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 -1+q 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q],
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  q 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 q 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 q 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  q-q^2 0 q-q^2 0 -1+q 0 0 0 -q+q^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0;
  0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 1 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1-q 0 0 1-q 0 0 0 -1+q 0 0 0 0 0 0 0 0 -1+q 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 1-q^-1 0 0 0 0 0 -1+q 0 0 0 0 1-q 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1-q 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1+q 0 0 0 0 0 0 0 0 0;
  (-1+2q)-q^2 (1-2q)+q^2 (q-2*q^2)+q^3 ((1-3q)+3*q^2)-q^3 (-1+2q)-q^2 0 (1-2q)+q^2 (1-2q)+q^2 j*q-j*q^2 (j-2j*q)+j*q^2 -j*q+j*q^2 (-q+2*q^2)-q^3 0 0 (q-2*q^2)+q^3 0 -1+(-3j2-2j)*q+(3j2+2j)*q^2+q^3 0 j*q-j*q^2 (-1+2q)-q^2 (-1+2q)-q^2 0 (q-2*q^2)+q^3 0 -j2*q+j2*q^2 0 q-q^2 0 0 0 0 q 0 q (-1+2q)-q^2 (q-2*q^2)+q^3 (-1+2q)-q^2 (1+(3j2+2j)*q+(-3j2-2j)*q^2)-q^3 (-j+2j*q)-j*q^2 -j*q+j*q^2;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 q 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0;
  -1+q 1-q -1+q 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 1-q 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1+q 0 0 0;
  0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1],
 [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
  0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  q 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 1-q 0 0 1-q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0;
  0 0 0 0 0 1 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 1 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 j 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
  -1+q 1-q -1+q 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1-q 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 q;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 j2*q 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0;
  q-q^2 0 q-q^2 0 -1+q 0 0 -q -q+q^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 q 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0;
  -2+q^-1+q 2-q^-1-q -1+2q-q^2 3-q^-1-3q+q^2 -2+q^-1+q 0 2-q^-1-q 2-q^-1-q -j+j*q 2j-j*q^-1-j*q j-j*q 1-2q+q^2 0 0 -1+2q-q^2 0 3j2+2j+q^-1+(-3j2-2j)*q-q^2 0 -j+j*q -2+q^-1+q -2+q^-1+q 0 -1+2q-q^2 0 j2-j2*q 0 -1+q 0 0 0 0 -1+q 0 -1 -2+q^-1+q -1+2q-q^2 -2+q^-1+q -3j2-2j-q^-1+(3j2+2j)*q+q^2 -2j+j*q^-1+j*q j-j*q;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0;
  0 (-1+2q)-q^2 (-q+2*q^2)-q^3 ((-3+q^-1+4q)-3*q^2)+q^3 0 (1-2q)+q^2 (-1+2q)-q^2 (-1+2q)-q^2 0 (-j+2j*q)-j*q^2 (3-q^-1)+(3j2+4j)*q+(-j2-2j)*q^2 ((-1+3q)-3*q^2)+q^3 -2+q^-1+q -j*q+j*q^2 (-q+2*q^2)-q^3 (3j2+2j+q^-1+(-3j2-2j)*q)-q^2 0 (1-2q)+q^2 -j*q+j*q^2 (1-2q)+q^2 0 (j-2j*q)+j*q^2 ((1-3q)+3*q^2)-q^3 (j2+(-2j2+j)*q+(2j2-j)*q^2)-j2*q^3 0 j-j*q -q+q^2 (j-2j*q)+j*q^2 -2+q^-1+q q-q^2 (1-2q)+q^2 -q -j2*q+j2*q^2 0 0 (((3-q^-1)-4q)+3*q^2)-q^3 0 -1+(-3j2-2j)*q+(3j2+2j)*q^2+q^3 0 0;
  0 0 0 0 0 0 0 0 0 0 -1+q 0 1-q 0 0 1-q 0 0 0 0 0 0 q 0 0 0 0 0 1 0 0 0 0 0 0 -1+q 0 0 0 0;
  0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]]
end

const rep335_8=let q=Pol();j=E(3);j2=E(3,2)
[[-j*q+j j*q 0 q^2-q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
j2*q-j2+j2*q^-1 -j2*q+j2 0 q^2+(2j+j2)*q-j 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 j*q^0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
q-2+q^-1 -q+1 -j2+j2*q^-1 0 0 0 0 -j2*q+j2 0 q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 -j2*q+j2 -j2*q^2+(j+3j2)*q-2j-3j2-q^-1 0 0 0 -j2*q+j2 0 j2*q^0 0 0 0 0 0 0 0 0 0 0 0 j2*q-j2 0 0 0 0 0 0 0;
j2*q^2-2j2*q+2j2-j2*q^-1 -j2*q^2+2j2*q-j2 -j*q+2j+j2+q^-1 -j2*q^2+2j2*q-j2 0 0 0 q^2-q+1 0 -j*q+j 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 -j2*q^2+(-j+j2)*q+j -j2*q^3+3j2*q^2-4j2*q+3j2-j2*q^-1 0 0 0 j*q^2-j*q+j 0 -j*q+j 0 0 0 0 0 0 0 0 -j2*q+j2 0 0 j2*q^2-2j2*q+j2 0 0 0 0 0 0 0;
0 0 j2-j2*q^-1 0 q^2-2*q+1 0 0 0 0 0 0 -j2*q+j2 0 0 j2*q^0 0 0 0 0 0 0 0 -q+1 0 0 0 0 0 0 0;
0 0 -q^2+(-3j-2j2)*q+3j+j2-j*q^-1 0 -q+1 0 0 0 0 0 0 0 -j2*q+j2 0 0 0 0 j*q 0 0 0 0 0 0 0 0 0 -q+1 0 0;
0 0 -q^2+(-3j-2j2)*q+(3j+j2)-j*q^-1 j2*q-j2 0 0 0 0 0 0 0 0 0 -j2*q+j2 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 -q+1 0 0;
0 0 j2*q+(j-j2)-j*q^-1 0 q^3+(2j+3j2)*q^2+(-j-3j2)*q+j2 0 0 0 0 0 0 j*q^2-j*q+j 0 0 -j*q+j 0 0 0 0 q-1 0 0 -q^2+2*q-1 0 0 0 0 0 0 0;
0 0 -q^2+3*q-4+3*q^-1-q^(-2) j2*q+(j-j2)-j*q^-1 0 0 0 0 0 0 0 0 0 j*q-j+j*q^-1 0 -j*q+j 0 0 0 0 0 0 0 0 0 0 j2*q-j2 -q+2-q^-1 0 0;
q^2-3*q+3-q^-1 -q^2+2*q-1 0 0 q^3-3*q^2+3*q-1 q^2-2*q+1 j2*q^3-2j2*q^2+j2*q -j2*q^2+(-j+j2)*q+j 0 q-1 0 -j2*q^2+(-j+j2)*q+j 0 0 j2*q-j2 0 -q^0 0 0 0 0 0 -q^2+2*q-1 0 -q+1 -j2*q^2+j2*q 0 0 0 0;
0 0 -j*q^2+3j*q-4j+3j*q^-1-j*q^-2 0 -j*q+2j+j2+q^-1 0 0 0 0 0 0 0 j2*q-j2+j2*q^-1 0 0 0 0 -j*q+j 0 0 0 0 0 0 0 0 q-1 -j*q+2j-j*q^-1 0 0;
q-2+q^-1 0 0 -q^2+3*q-3+q^-1 j*q^3-4j*q^2+(6j-j2)*q+(-4j+2j2)+(j-j2)*q^-1 j*q^2-3j*q+(3j-j2)+(-j+j2)*q^-1 q^3-3*q^2+(-4j-3j2)*q+2j+j2 j*q^3+(-2j+2j2)*q^2+(j-5j2)*q+4j2-j2*q^-1 0 -j2*q^2+(j+3j2)*q+(-2j-3j2)-q^-1 0 -q^2+(-2j-3j2)*q+(j+3j2)-j2*q^-1 -j*q+(2j+j2)+q^-1 q^2-3*q+3-q^-1 q-2+q^-1 -q^2+2*q-1 0 q-1 -q^0 0 j*q^2+(-j+j2)*q-j2 0 -j*q^2+3j*q-3j+j*q^-1 0 -j*q+2j-j*q^-1 -q^2+2*q-1 q^2-2*q+1 -j2*q+2j2-j2*q^-1 0 j*q-j;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 0 0 0 0 0;
-q^2+3*q-4+2*q^-1 q^2-2*q+1 0 q-2+q^-1 0 0 0 (-j+j2)*q^2+(2j-3j2)*q-j+3j2-j2*q^-1 0 (j+2j2)*q+(-2j-3j2)-q^-1 0 0 0 -q+2-q^-1 0 q-1 0 0 0 0 -j*q+j 0 0 0 0 0 -q+1 0 0 -j*q^0;
0 0 0 0 -q^3+3*q^2+(4j+3j2)*q-3j-j2+j*q^-1 0 0 0 0 0 0 (-j+j2)*q^2+(2j-3j2)*q-j+3j2-j2*q^-1 q-2+q^-1 0 -j2*q+2j2-j2*q^-1 0 0 q-1 0 -1+q^-1 0 -j*q+j q^2-2*q+1 0 0 0 q-1 0 -j2*q^0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q^0 0 0 q-1 0 0 0 0 0 0 0;
j2*q^3+(-j-5j2)*q^2+(3j+9j2)*q-3j-7j2+(j+2j2)*q^-1 -j2*q^3+3j2*q^2-3j2*q+j2 0 -q^3+4*q^2-6*q+4-q^-1 -q^4+(-5j-4j2)*q^3+(10j+5j2)*q^2+(-10j-j2)*q+(5j-2j2)+(-j+j2)*q^-1 -q^3+(-4j-3j2)*q^2+(6j+2j2)*q-4j+j2+(j-j2)*q^-1 -j2*q^4+(j+4j2)*q^3+(-4j-6j2)*q^2+(5j+4j2)*q+(-2j-j2) j*q^4+(-4j+2j2)*q^3+(6j-6j2)*q^2+(-4j+7j2)*q+(j-4j2)+j2*q^-1 -j*q^2+(2j+j2)*q+1 -j2*q^3+(j+5j2)*q^2+(-3j-8j2)*q+3j+5j2+q^-1 j*q-j (-j+j2)*q^3+(3j-4j2)*q^2+(-3j+6j2)*q+(j-4j2)+j2*q^-1 q^2-3*q+3-q^-1 q^3+(3j+4j2)*q^2+(-4j-5j2)*q-3+q^-1 -j2*q^2+3j2*q-3j2+j2*q^-1 -q^3+(-2j-3j2)*q^2+(2j+3j2)*q+1 0 q^2-2*q+1 0 -q+2-q^-1 j*q^3+(-2j+j2)*q^2+(j-2j2)*q+j2 -j*q^2+(j-j2)*q+j2 q^3+(4j+3j2)*q^2+(-5j-3j2)*q+(2j+j2) -q^0 q^2+(3j+2j2)*q-3j-j2+j*q^-1 j2*q^3+(-j-3j2)*q^2+(2j+3j2)*q+1 q^3-2*q^2+q -j2*q^2+2j2*q-j2 -j2*q+j2 j*q^2-2j*q+j;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q^0 q-1 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q^0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q q-1 0 0;
0 0 q^3-3*q^2+4*q-4+3*q^-1-q^-2 0 -j*q^2+2j*q-j 0 0 0 0 0 0 -j*q^3+3j*q^2-4j*q+3j-j*q^-1 q^2-2*q+2-q^-1 0 -j*q+2j-j*q^-1 0 0 q^2-2*q+1 0 0 0 -j*q^2+j*q-j -j2*q+j2 0 0 0 -j2*q^2+2j2*q-j2 q^2-2*q+2-q^-1 -j2*q+j2 0;
q^3-3*q^2+5*q-5+2*q^-1 -q^3+3*q^2-4*q+2 q^2-2*q+1 -q^3+3*q^2+(3j+2j2)*q-j+j2-j2*q^-1 0 0 0 (j-j2)*q^3+(-2j+3j2)*q^2+(2j-4j2)*q+(-j+3j2)-j2*q^-1 0 (-j-2j2)*q^2+(3j+5j2)*q+(-3j-4j2)-q^-1 0 0 0 q^2-2*q+2-q^-1 0 -q^2+2*q-1 0 0 0 0 -j2*q^2+j2*q-j2 0 0 0 0 0 0 -j*q+j 0 -j2*q+j2],
[0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
q^0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
q^2-2*q+1 -q^2+q 0 0 j*q-j q-1 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
q-1 0 0 0 -j*q+j q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 j*q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 j2*q 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 q 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 q 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 q^0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 j*q^0 0 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 -q^2+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q-1 0 0 q 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q+1 0 0 -q+1 0 0 0 0 0 0 0 0 0 0 -j*q^0 0;
0 0 0 0 q^2-q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q^0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0;
q^3+(4j+3j2)*q^2+(-5j-3j2)*q+3j+j2-j*q^-1 -q^3+(-3j-2j2)*q^2+(3j+j2)*q-j 0 -j*q^2+2j*q-j (j-j2)*q^2+(-2j+3j2)*q+(j-3j2)+j2*q^-1 0 (-2j-j2)*q^2+(3j+j2)*q-j 0 -j*q^2+2j*q-j 0 j*q-j 0 0 0 0 0 0 0 0 j*q-j 0 0 0 0 q-1 j*q 0 0 0 0;
q^2-2*q+1 -q+1 0 q^2-2*q+1 q^2-2*q+1 (-j-2j2)*q+(2j+3j2)+q^-1 0 0 -q+1 0 0 0 0 0 0 0 0 0 0 -q+2-q^-1 0 0 -q+1 0 j2*q^0 0 0 0 0 0;
0 0 j*q-j 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q-1 j*q^0 0 0;
0 0 -q^2+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 -j2*q^2+j2*q -q+1 0 0 0 0 0 0 0 0 -j2*q 0 0 0 0 0 0 q-1 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0 0 0 0 0 0 0 q-1],
[0 0 0 -q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
-q+1 q 0 -q+1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
-q^0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 j2*q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 j*q q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 q^0 0 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 j*q-j 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
-q+1 0 0 q^2-2*q+1 j*q^2-2j*q+j j2*q^2+(-j-3j2)*q+(2j+3j2)+q^-1 q^2-q 0 q^2-2*q+1 0 0 0 0 0 0 0 0 0 0 -q+1 0 0 0 0 -j2*q+j2 -q 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 j2*q-j2 0 0 0 0 0 0 0 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 j2*q+(-j-2j2)-q^-1 0 0 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 -j*q^0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1-q^-1 0 -j2*q^0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j*q+j 0 0 -j*q q-1 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 q-1 0 0 -q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q-1 0 0 0 0 0 0 q^0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q^0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 q^2-q 0 0 q^2-q 0 0 0 0 0 q-1 0 0 0 0 j*q 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 q-1 0 0 0 0 0;
0 0 0 q-1 q^2+(3j+2j2)*q-3j-j2+j*q^-1 0 q-1 0 q-1 0 -q^0 0 0 0 0 0 0 0 0 0 0 0 -q+1 0 0 q-1 0 0 0 0;
0 0 0 0 0 0 0 -j2*q+j2 0 0 0 0 0 q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2-j2*q^-1 0 j*q-j 0 0 0 0 j2*q^0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 j*q^2+(-j+j2)*q-j2 0 0 0 0 0 -j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q^0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 -q+1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1-q^-1 0 0 q^0 0 0 0 0 0 0 0;
0 0 -q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 j*q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 j2*q q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 -q+1 -q^2+(-3j-2j2)*q+3j+j2-j*q^-1 0 0 0 -q+1 0 q^0 0 0 0 0 0 0 0 0 0 0 0 q-1 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 q^0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
q^2-2*q+1 -q^2+q 0 0 j*q-j q-1 q 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 q 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 q-1 -q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 j*q 0 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j2*q^0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j*q 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q+1 0 q-1 0 0 0 0 -q^0 0 0 0 0 0 0;
q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q^0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q q-1 0 0 0 0 0 0 0 0;
-q+1 q 0 j2*q^2-j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q-1 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^2+q 0 -q 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 -q^2+q 0 0 q^2-q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0;
0 0 0 -q^2+2*q-1 -q^3+(-4j-3j2)*q^2+(6j+3j2)*q+(-4j-j2)+j*q^-1 0 -q^2+q 0 -q^2+2*q-1 0 q-1 0 0 0 0 0 0 0 0 0 0 0 q^2-2*q+1 0 0 q 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0;
0 0 0 0 0 0 0 0 0 -q^2+q 0 0 0 0 0 q^2-q 0 0 0 0 0 0 0 0 0 0 0 0 q-1 j2*q;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 q-1 0 0 q-1 0 0 0 0 0 0 0 0 0 0 j*q^0 0],
[-q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 -q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 q-1 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 -q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 j*q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 j2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
q-2+q^-1 -q+1 -j2+j2*q^-1 0 0 0 0 -j2*q+j2 0 q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 j*q^0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
q^2-2*q+1 -q^2+q 0 0 j*q-j q-1 q 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 -j*q^0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 -j2*q q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 q 0 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 q^0 0 0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q-1 0 -j2*q 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 q-1 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j*q^0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q-1 0 0 0 -j*q^0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 q^2-q j*q-j 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j*q^0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q-1 0 j2*q^2-j2*q 0 0 0 0 q 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j2*q 0 0 0 0 0 0 0 0 0;
q^2-3*q+4-2*q^-1 -q^2+2*q-1 0 -q+2-q^-1 0 0 0 (j-j2)*q^2+(-2j+3j2)*q+j-3j2+j2*q^-1 0 (-j-2j2)*q+2j+3j2+q^-1 0 0 0 q-2+q^-1 0 -q+1 0 0 0 0 j*q-j 0 0 0 0 0 q-1 0 0 j*q^0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q^0 0 0 0 0 0 0 q-1 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j2*q 0 0 0 0 q-1 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 j2*q-j2 0 0 j2*q-j2 0 0 0 0 0 0 0 0 0 0 q 0;
j*q^3+(-3j+j2)*q^2+(3j-2j2)*q-j+2j2-j2*q^-1 -j*q^3+(2j-j2)*q^2+(-j+2j2)*q-j2 0 -j2*q^2+2j2*q-j2 (j+2j2)*q^2+(-3j-5j2)*q+3j+4j2+q^-1 0 (j-j2)*q^2+(-j+2j2)*q-j2 0 -j2*q^2+2j2*q-j2 0 j2*q-j2 0 0 0 0 0 0 0 0 j2*q-j2 0 0 0 0 j*q-j j2*q 0 0 0 q-1]]
end

chevieset(:imp,:HeckeRepresentation,function(p,q,r,para,rootpara,i;gen=false)
  if !(para isa Vector) para=[para] end
  if (q,r)==(1,2) && gen==false
    Y,X=para;x1,x2=X
    t=partition_tuples(2,p)[i]
    if count(!isempty,t)==1
      p=findfirst(!isempty,t)
      return (0+x1^0)*[[Y[p];;],[t[p]==[2] ? x1 : x2;;]]
    else
      p=findall(!isempty,t)
      y1,y2=Y[p]
      return x1^0*[[y1 0;-1 y2],[x1 x1*y1+x2*y2;0 x2]]
    end
  elseif (p,q,r)==(3,3,3)
    x=-para[2][1]//para[2][2]
    f(x,j)=[[-1 0 0;0 0 1;0 x -1+x],[-1 0 0;x-x^2 -1+x j^2;j*x-j*x^2 j*x 0],
            [0 1 0;x -1+x 0;0 0 -1]]
    r=x^0*[[[-1 0;-1 x],[x -x;0 -1],[x -x;0 -1]],
          [[-1 0;-1 x],[x -x;0 -1],[-1 0;-1 x]],
          [[-1 0;-1 x],[x -x;0 -1],[-1+x 1;x 0]],
          f(x,E(3)),f(x,E(3,2)),
          [[-1;;],[-1;;],[-1;;]],
          -x*f(x^-1,E(3,2)),-x*f(x^-1,E(3)),
          [[-1 0;-1 x],[-1 0;-1 x],[x -x;0 -1]],[[x;;],[x;;],[x;;]]]
    return -para[2][2]*r[i]
  elseif (p,q,r)==(2,2,4)
    x=-para[1][1]//para[1][2]
    r=[x->[[-1+x -1 0;-x 0 0;x-x^2 -1+x -1],[0 1 0;x -1+x 0;0 0 -1],
           [-1 0 0;0 0 1;0 x -1+x],[0 1 0;x -1+x 0;0 0 -1]],
       x->[[0 1 0;x -1+x 0;0 0 -1],[-1+x -1 0;-x 0 0;x-x^2 -1+x -1],
           [-1 0 0;0 0 1;0 x -1+x],[0 1 0;x -1+x 0;0 0 -1]],
       x->[[-1 0 0 0;0 -1+x -1 0;0 -x 0 0;0 0 0 -1],
           [-1 1-x 1-x 0;0 0 1 0;0 x -1+x 0;0 -1+x -1+x -1],
           [-1+x -x 0 0;-1 0 0 0;0 0 -1 0;0 0 0 -1],
           [0 0 0 1;0 -1 0 0;0 0 -1 0;x 0 0 -1+x]],
       x->[[-1;;],[-1;;],[-1;;],[-1;;]],
       x->[[x 1-x -1+x -x+x^2 x-x^2 0;0 -1+x 0 0 -x x-x^2;
            0 0 -1+x -x 0 x-x^2;0 0 -1 0 0 -1+x;0 -1 0 0 0 -1+x;0 0 0 0 0 -1],
           [x 0 0 0 0 0;0 0 0 0 x 0;0 0 0 x 0 0;0 0 1 -1+x 0 0;
            0 1 0 0 -1+x 0;0 0 0 0 0 -1],
           [0 0 x 0 0 0;0 -1 0 0 0 0;1 0 -1+x 0 0 0;0 0 0 x 0 0;0 0 0 0 0 x;
            0 0 0 0 1 -1+x],
           [-1 0 0 0 0 0;0 -1+x 1 0 0 0;0 x 0 0 0 0;0 0 0 0 x 0;
            0 0 0 1 -1+x 0;0 0 0 0 0 x]], 
       x->[[-1+x 0 -1 0 0 0 0 0;0 0 0 0 1 0 0 0;-x 0 0 0 0 0 0 0;
            0 0 0 0 0 0 1 0;0 x 0 0 -1+x 0 0 0;0 0 0 0 0 -1+x 0 x;
            0 0 0 x 0 0 -1+x 0;0 0 0 0 0 1 0 0],
           [0 0 1 0 0 0 0 0;0 0 -1+x 0 1 0 -x^-1+1 0;x 0 -1+x 0 0 0 0 0;
            0 0 0 -1+x 0 0 -1 0;x-x^2 x 0 -1+x -1+x 0 x^-1-2+x 0;
            -x+x^2 0 0 -x+x^2 0 -1+x 1-x x;0 0 0 -x 0 0 0 0;
            0 0 1-x x-x^2 0 1 -1+x 0],
           [0 1 0 0 0 0 0 0;x -1+x 0 0 0 0 0 0;0 0 0 1 0 0 0 0;
            0 0 x -1+x 0 0 0 0;0 0 0 0 -1+x 0 -1 0;0 0 0 0 0 x 0 0;
            0 0 0 0 -x 0 0 0;0 0 0 0 0 -x 0 -1],
           [-1 0 0 0 0 0 1 0;0 0 0 0 0 0 0 1;0 0 -1 -x 0 0 0 0;
            0 0 0 x 0 0 0 0;0 0 0 0 0 1 0 0;0 0 0 0 x -1+x 0 0;
            0 0 0 0 0 0 x 0;0 x 0 0 0 0 0 -1+x]], 
       x->[[-1 -1 0;0 x 0;0 1 -1],[-1 -1 0;0 x 0;0 1 -1],
           [-1+x x 0;1 0 0;0 0 -1],[0 0 1;0 -1 0;x 0 -1+x]],1,2,
       x->[[x 0;-1 -1],[x 0;-1 -1],[0 1;x -1+x],[x 0;-1 -1]],3,7,4]
    if r[i] isa Integer return para[1][2]*x*r[r[i]](x^-1)
    else return -para[1][2]*r[i](x)*x^0
    end
  elseif (p, q, r) == (3, 3, 4)
    x=-para[2][1]//para[2][2]
    function m334(i)
      function f1(x) x^0*[[x -1 0 0 0 0 0 0 0 0 1-x-x^2+x^3 0;
        0 -1 0 0 0 0 0 0 0 0 0 0;0 0 -1+x 0 x 0 -x 0 0 0 x-x^2 0;
        0 0 0 -1+x 0 0 -x 0 0 0 x-x^2 0;0 0 1 -1 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 -x;0 0 0 -1 0 0 0 0 0 0 -1+x 0;
        0 0 0 0 0 0 0 -1 0 0 0 0;0 0 0 0 0 0 0 0 -1+x 1 -1+x 0;
        0 0 0 0 0 0 0 0 x 0 -1+x 0;0 0 0 0 0 0 0 0 0 0 -1 0;
        0 0 0 0 0 -1 0 0 0 0 0 -1+x],
       [0 x 0 0 0 0 0 0 0 0 0 0;1 -1+x 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 x 0 0 0 0 0;0 0 0 0 x 0 0 1-x 0 0 0 0;
        0 0 0 1 -1+x 0 0 1-x 0 0 0 0;0 0 0 0 0 0 0 1-x -1+x 1 0 1-x+x^2;
        0 0 1 0 0 0 -1+x 0 0 0 0 0;0 0 0 0 0 0 0 -1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 x;0 0 0 0 0 x 0 x-x^2 -x -1+x 0 x-x^2;
        0 0 0 0 0 0 0 0 0 0 -1 0;0 0 0 0 0 0 0 0 1 0 0 -1+x],
       [0 -1+2x-x^2 1-x x -x+x^2 0 0 0 -1+2x-x^2 0 0 0;
        0 -1+x 1 0 0 0 0 0 -1+x 0 0 0;
        0 x 0 0 0 0 0 0 -1+x 0 0 0;1 -1+x 0 -1+x 0 0 1-x 0 0 0 0 0;
        0 0 0 0 -1+x 0 1 0 0 0 0 0;0 0 0 0 0 -1 0 0 0 0 0 0;
        0 0 0 0 x 0 0 0 0 0 0 0;0 0 0 0 0 0 0 -1+x 0 0 0 -1;
        0 0 0 0 0 0 0 0 -1 0 0 0;0 0 0 0 0 0 0 0 0 0 x 0;
        0 0 0 0 0 0 0 0 0 1 -1+x 0;0 0 0 0 0 0 0 -x 0 0 0 0],
       [-1 0 0 0 0 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 x 0 0 0;0 0 0 0 0 x 0 0 0 0 1-x x-x^2;
        0 0 0 0 0 0 0 0 0 1 0 x;0 0 0 1 0 -1+x -1+x 0 0 0 1-x 0;
        0 0 0 0 0 0 0 0 0 0 0 x;0 0 0 0 0 0 0 x 0 0 -1 0;
        0 0 1 0 0 0 0 0 -1+x 0 0 0;0 0 0 0 x 0 -x 0 0 -1+x 0 0;
        0 0 0 0 0 0 0 0 0 0 -1 0;0 0 0 0 0 0 1 0 0 0 0 -1+x]]
      end
      function f2(x,j) 
      [[-1 0 0 0;0 -1 0 0;0 0 -1 0;0 x x x],
       [-1 0 0 0;0 -1 0 0;0 -j^2 x 1;0 0 0 -1],
       [-1 0 0 0;x x -j*x 1;0 0 -1 0;0 0 0 -1],
       [x 1 0 0;0 -1 0 0;0 0 -1 0;0 0 0 -1]]
      end
      function f3(x) 
      [[x -1 0 -x 0 0;0 -1 0 0 0 0;0 0 -1+x 0 1 x;
        0 0 0 -1 0 0;0 0 x 0 0 x;0 0 0 0 0 -1],
       [-1 0 0 0 0 0;0 0 0 0 x 1;0 0 -1 0 0 0;-1 0 -1 x 0 -1+x;
        0 1 0 0 -1+x 1;0 0 0 0 0 -1],
       [0 x 1 -1 -1 0;1 -1+x 1 -1 -1 0;0 0 -1 0 0 0;0 0 0 -1 0 0;
        0 0 0 0 -1 0;0 0 0 1 1 x],
       [x -1 0 0 1 x;0 -1 0 0 0 0;0 0 -1 0 0 0;0 0 -1 x 1 x;
        0 0 0 0 -1 0;0 0 0 0 0 -1]]
      end
      f5(x)=[[-1;;],[-1;;],[-1;;],[-1;;]]
      function f7(x,j) 
      [[-1 0 0 0 0 0;x x 0 0 0 0;x 0 x 0 0 0;0 0 0 -1 0 0;
        0 0 0 0 -1 0;0 0 0 -j*x^2 x x],
       [x 1 0 0 0 0;0 -1 0 0 0 0;0 j^2 x 0 0 0;0 0 0 -1 0 0;
        0 0 0 x x 1;0 0 0 0 0 -1],
       [x 0 1 0 1 0;0 x j*x 0 0 1;0 0 -1 0 0 0;0 0 0 x 1 -j^2*x^-1;
        0 0 0 0 -1 0;0 0 0 0 0 -1],
       [-1 0 0 0 0 0;0 -1 0 0 0 0;0 0 -1 0 0 0;0 0 0 x 0 0;
        x 0 0 0 x 0;0 x 0 0 0 x]]
      end
      function f8(x,j) 
      [[-1 0 0 0 0 0 0 0;1 x 0 0 0 0 1 0;1 0 x 0 0 0 0 0;
        0 0 0 -1 0 0 0 0;0 0 0 0 -1 0 0 0;0 0 0 -j*x x x 0 0;
        0 0 0 0 0 0 -1 0;0 0 0 0 (j^2-j)*x 0 1 x],
       [x x 0 0 0 0 -j^2 0;0 -1 0 0 0 0 0 0;0 j x 0 0 0 0 0;0 0 0 -1 0 0 0 0;
        0 0 0 1 x 1 0 0;0 0 0 0 0 -1 0 0;0 0 0 0 0 0 -1 0;
        0 0 0 0 0 -2j^2-j 1 x],
       [x 0 x 0 x 0 0 0;0 x j^2*x 0 0 1 0 0;0 0 -1 0 0 0 0 0;
        0 0 0 x x -j^2 0 0;0 0 0 0 -1 0 0 0;0 0 0 0 0 -1 0 0;0 0 0 0 0 0 x x;
        0 0 0 0 0 0 0 -1],
       [-1 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0;0 0 -1 0 0 0 0 0;0 0 0 x 0 0 -j^2 0;
        1 0 0 0 x 0 0 0;0 x 0 0 0 x 0 0;0 0 0 0 0 0 -1 0;
        0 0 (j^2-j)*x 0 0 0 1 x]]
      end
      function f11(x) 
      [[x 1 0;0 -1 0;0 0 -1],[x 1 0;0 -1 0;0 0 -1],
       [-1 0 0;x x 1;0 0 -1],[-1 0 0;0 -1 0;0 x x]]
      end
      f13(x)=[[-1 0;x x],[-1 0;x x],[x 1;0 -1],[-1 0;x x]]
      r=[f1(x),f2(x,E(3)),f3(x),f2(x,E(3,2)),f5(x),-x*f1(x^-1),f7(x,E(3)),
         f8(x,E(3)),f8(x,E(3,2)),-x*f7(x^-1,E(3)),f11(x),-x*f3(x^-1),f13(x),
         -x*f2(x^-1,E(3,2)),-x*f2(x^-1,E(3)),-x*f11(x^-1),-x*f5(x^-1)]
      return x^0*r[i]
    end
    return -para[2][2]*m334(i)
  elseif (p,q,r)==(3,3,5)
    function f2(q)
     [[-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      -E(3)+root(-3)*q+E(3,2)*q^2 E(3,2)-E(3,2)*q 0 0 -E(3)+E(3)*q 0 root(-3)-E(3)*q^-1+E(3,2)*q -E(3,2)+E(3,2)*q^-1+E(3,2)*q 0 0 2*E(3)-E(3)*q^-1-E(3)*q 0 0 0 0 0 0 0 0 0;
      0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 E(3,2) 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      E(3)*q-E(3)*q^2 E(3)*q 0 0 0 0 E(3)-E(3)*q E(3)-E(3)*q 0 0 E(3)-E(3)*q 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 E(3) 0 0 0 0 0 0 0;
      0 0 0 0 E(3)*q 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 E(3,2)-E(3,2)*q 0 E(3,2)-E(3,2)*q 0 0 E(3,2)-E(3,2)*q 0 0 0 E(3,2) 0 0 0 0;
      0 0 0 0 0 0 0 0 0 E(3,2)*q 0 0 -1+q 0 0 0 0 0 0 0;
      -E(3)+E(3)*q 0 0 0 0 E(3)-E(3)*q 0 0 0 0 0 0 0 -1+q -E(3) 0 0 0 0 0;
      1-q 0 0 0 0 -1+q 0 0 0 0 0 0 0 -E(3,2)*q 0 0 0 0 0 0;
      0 0 0 0 0 0 -E(3,2)-root(-3)*q+E(3)*q^2 0 -E(3,2)-root(-3)*q+E(3)*q^2 0 0 E(3)-E(3)*q+E(3)*q^2 0 0 0 E(3)-E(3)*q 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3) 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 E(3,2)*q -1+q 0 0;
      E(3)*q-2*E(3)*q^2+E(3)*q^3 E(3)*q-E(3)*q^2 0 0 0 0 E(3)-2*E(3)*q+E(3)*q^2 E(3)-2*E(3)*q+E(3)*q^2 0 0 E(3)-2*E(3)*q+E(3)*q^2 0 -E(3)+E(3)*q 0 0 0 0 -E(3)+E(3)*q 0 E(3,2)*q;
      0 0 0 0 0 0 0 -E(3)+E(3)*q 0 E(3)-E(3)*q 0 0 0 0 0 0 E(3)-E(3)*q 0 E(3) -1+q],
     [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 -1+q 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 q^3-q^4 1-q -1+q 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0;
      (1-2q)+q^2 0 q^3-q^4 0 0 0 0 0 1-q -1+q 0 0 1 1-q 0 0 0 0 0 0;
      0 0 q^3-q^4 1-q q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
      (1-2q)+q^2 0 q^3-q^4 0 0 0 0 0 1-q q 0 0 0 0 1-q 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 -q -1+q 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 -1+q 0 0 0 0;
      0 0 0 1-q 0 -1+q 0 0 0 0 0 2-q^-1-q 0 0 0 1-q^-1 -1+q 1 0 0;
      0 0 0 1-q 0 -1+q 0 0 0 0 0 1-q 0 0 0 0 q 0 0 0;
      -1+2q-q^2 0 0 0 q-q^2 1-2q+q^2 0 0 0 0 0 -1+q 0 0 -1+q 0 0 0 -1+q q;
      2-q^-1-q 0 -q^3+2*q^4-q^5 -1+2q-q^2 1-2q+q^2 -2+q^-1+q 0 0 0 0 -1+q 0 0 -1+q 2-q^-1-q -1+q^-1 0 0 1 0],
     [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 1 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      1-q 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
      0 0 0 0 0 0 0 -q 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 -1+q 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 1 -1+q 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0;
      1-q 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0;
      0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -1+q 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 -1+q 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1+q],
     [0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 -q^-2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -q^3 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
      -q 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 q 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0 0;
      0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0;
      -1+q 0 0 0 0 1-q 0 0 0 0 0 -1+q 0 0 -1 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0;
      0 0 0 0 0 0 1-q 0 1-q 0 0 1-q 0 0 0 1 0 0 0 0;
      0 0 0 0 0 0 1-q 0 1-q 0 0 -q 0 0 0 0 0 0 0 0;
      1-q 0 0 0 0 -1+q 0 0 0 0 0 0 0 q 1-q -1+q 0 0 0 0;
      0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 -1+q 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 -1+q 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1],
     [-1+q 0 -q^3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      1-2q+q^2 0 q^3-q^4 0 0 0 0 0 1-q -1+q 0 0 1 1-q 0 0 0 0 0 0;
      -q^-2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      1-q 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
      0 0 0 1 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 -q 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 -1+q 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 1 0 -1+q 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0;
      q-q^2 q 0 0 0 0 1-q 1-q 0 0 1-q 0 -1+q 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 -1 0 0 -1+q 0 0 0 0 0 0;
      0 0 q^3-q^4 1-q q 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]]
    end
    function f3(q,j)
    [[-1 0 0 0 0;0 -1 0 0 0;1 0 0 0 -1;1+j*q 0 1+j*q -1 -1-j*q;-q 0 -q 0 -1+q],
     [-1 0 0 0 0;0 -1 0 0 0;1 0 0 -1 0;-q 0 -q -1+q 0;-j 0 -j j -1],
     [-1 0 -1 0 0;0 -1 1 0 0;0 0 q 0 0;0 0 0 -1 0;0 0 0 0 -1],
     [q 0 0 0 0;-1 -1 0 0 0;-q 0 -1 0 0;1 0 0 -1 0;1 0 0 0 -1],
     [0 1 0 0 0;q -1+q 0 0 0;0 0 -1 0 0;1 1 0 -1 0;1 1 0 0 -1]]
    end
    function f4(q,j)
    [[-1+q 0 0 q 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 1;0 0 -1 0 0 0 0 0 0 0;1 0 0 0 0 0 0 0 0 0;
      0 0 j^2*q+(-j^2+j)*q^2-j*q^3 0 j-j*q 0 -j+j*q-j*q^2 0 0 0;0 0 0 0 0 -1 0 0 0 0;
      0 0 j^2*q-j^2*q^2 0 -j^2 0 j^2-j^2*q 0 0 0;
      -q+q^2 0 j^2*q-2j^2*q^2+j^2*q^3 -q+q^2 -j^2+j^2*q 0 -j+(-j^2+j)*q+j^2*q^2 -1 0 0;
      0 q-q^2 j^2*q-2j^2*q^2+j^2*q^3 0 -j^2+j^2*q 0 -j+(-j^2+j)*q+j^2*q^2 0 -1 q-q^2;
      0 q 0 0 0 0 0 0 0 -1+q],
     [0 0 j^2*q-j^2*q^2 j^2*q 0 0 0 0 0 0;0 -1+q -q+q^2 0 0 0 0 0 0 j;0 0 -1 0 0 0 0 0 0 0;
      j 0 q-q^2 -1+q 0 0 0 0 0 0;0 0 0 0 -1+q 0 -q 0 0 0;
      q-q^2 -j^2*q+j^2*q^2 j^2*q-j^2*q^2-j^2*q^3+j^2*q^4 j^2*q^2-j^2*q^3 0 -1 0 0 0 -1+q;
      0 0 0 0 -1 0 0 0 0 0;0 0 0 0 0 0 0 -1 0 0;0 0 0 0 0 0 0 0 -1 0;0 j^2*q -j^2*q+j^2*q^2 0 0 0 0 0 0 0],
     [-1 0 0 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0;0 0 q -1+q 0 0 0 0 0 0;
      0 0 0 0 0 0 0 -1 0 0;0 0 0 0 0 -1+q 0 0 0 q;0 0 0 0 0 0 -1 0 0 0;0 0 0 0 -q 0 0 -1+q 0 0;
      0 0 0 0 0 0 0 0 -1 0;0 0 0 0 0 1 0 0 0 0],
     [0 0 0 0 0 0 0 0 0 1;0 -1+q 0 q 0 0 0 0 0 0;0 0 -1 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0;
      0 0 0 0 -1 0 0 0 0 0;0 0 0 0 0 -1 0 0 0 0;0 0 0 0 0 0 -1 0 0 0;0 0 0 0 0 0 0 0 -1 0;
      0 0 0 0 0 0 0 -q -1+q 0;q 0 0 0 0 0 0 0 0 -1+q],
     [-1 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 -1 0 0 0;0 0 -1 0 0 0 0 0 0 0;0 0 0 -1 0 0 0 0 0 0;
      0 -q+q^2 -q^2+q^3 0 -1+q 0 0 0 0 j*q;0 0 0 -j^2*q+j^2*q^2 0 0 j^2-j^2*q -j^2 0 0;
      0 -q 0 0 0 0 -1+q 0 0 0;0 -q+q^2 0 q^2-q^3 0 -j*q 0 -1+q 0 0;0 0 0 0 0 0 0 0 -1 0;
      0 0 -j^2*q+j^2*q^2 0 j^2 0 -j^2+j^2*q 0 0 0]]
    end
    function f11(q,j)
    [[q 0 0 0 0 0 0 0 0 0;-j^2+j^2*q j^2-j^2*q j^2 0 0 0 0 0 0 0;
      -j+(-j^2+j)*q+j^2*q^2 (j-j*q)+j*q^2 j-j*q 0 0 0 0 0 0 0;0 0 0 0 0 1 0 0 0 0;0 0 0 0 0 0 1 0 0 0;
      0 0 0 q 0 -1+q 0 0 0 0;0 0 0 0 q 0 -1+q 0 0 0;0 0 0 0 0 0 0 -1 0 0;0 0 0 0 0 0 0 0 -1 0;
      0 0 0 0 0 0 0 0 0 -1],
     [q 0 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0;0 q -1+q 0 0 0 0 0 0 0;1-q 0 0 -1+q 0 j 0 0 0 0;
      1-q 0 0 0 -1+q 0 j 0 0 0;-j^2*q+j^2*q^2 0 0 j^2*q 0 0 0 0 0 0;-j^2*q+j^2*q^2 0 0 0 j^2*q 0 0 0 0 0;
      (-j^2+2j^2*q)-j^2*q^2 j^2-j^2*q j^2-j^2*q -j^2*q+j^2*q^2 0 -1+q 0 -1 0 0;
      0 0 0 j^2*q-j^2*q^2 -j^2*q+j^2*q^2 1-q -1+q 0 -1 0;
      (-j^2+2j^2*q)-j^2*q^2 j^2-j^2*q j^2-j^2*q 0 -j^2*q+j^2*q^2 0 -1+q 0 0 -1],
     [0 -1 0 0 0 0 0 0 0 0;-q -1+q 0 0 0 0 0 0 0 0;0 0 q 0 0 0 0 0 0 0;0 0 0 -1 0 0 0 0 0 0;
      0 0 0 0 -1 0 0 0 0 0;0 0 0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 0 0 1;0 0 0 0 0 q 0 -1+q 0 0;
      0 0 0 0 0 0 0 0 -1 0;0 0 0 0 0 0 q 0 0 -1+q],
     [-1 0 0 0 0 0 0 0 0 0;0 0 0 -1 0 0 0 0 0 0;-1+q 0 0 1-q 0 -j 0 0 0 0;0 -q 0 -1+q 0 0 0 0 0 0;
      0 0 0 0 -1 0 0 0 0 0;j^2*q-j^2*q^2 -j^2*q+j^2*q^2 -j^2*q 0 0 -1+q 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 0;0 0 0 0 0 0 0 q 0 0;0 0 0 0 0 0 0 0 -1+q -q;0 0 0 0 0 0 0 0 -1 0],
     [-1 0 0 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0;0 0 -1 0 0 0 0 0 0 0;0 0 0 0 -1 0 0 0 0 0;
      0 0 0 -q -1+q 0 0 0 0 0;0 0 0 0 0 0 -1 0 0 0;0 0 0 0 0 -q -1+q 0 0 0;0 0 0 0 0 0 0 0 0 -1;
      0 0 0 0 0 0 0 0 q 0;0 0 0 0 0 0 0 -q 0 -1+q]]
    end
    function f12(q,j)
    [[0 0 0 0 -j^2 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 -j^2*q 0 0 0 0 0 0 0 0 0;
      0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 -j^2 0 0 0 0 0 0 0 0;
      -j*q 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0;0 -j 0 0 0 -1+q 0 0 0 0 0 0 0 0 0;
      0 0 0 -j*q 0 0 -1+q 0 0 0 0 0 0 0 0;-j^2+j^2*q 0 0 0 j-j*q 0 0 -1 0 0 0 0 0 0 0;
      0 0 0 j^2-j^2*q 0 0 -j+j*q 0 -1 0 0 0 0 0 0;
      0 (2j^2-j^2*q^-1)-j^2*q 0 0 0 (j-2j*q)+j*q^2 0 0 0 -1 0 0 0 0 0;
      0 (-1+2q)-q^2 -q 0 0 (1-2q)+q^2 0 0 0 q q 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0;j-j*q 0 0 -q+q^2 (j^2+2j+q^-1)-j*q 0 1-q -q -q 0 0 0 q 0 0;
      0 -2j^2+j^2*q^-1+j^2*q 0 0 0 (-j+2j*q)-j*q^2 0 0 0 0 0 0 0 -1 0;
      0 (1-2q)+q^2 0 0 0 (-1+2q)-q^2 0 0 0 0 0 q 0 q q],
     [-1+q 0 0 0 -1 0 0 0 0 0 0 0 0 0 0;0 -1+q 0 0 0 -q 0 0 0 0 0 0 0 0 0;
      -1+q 0 -1 0 -1+q^-1 0 0 0 0 0 0 0 0 0 0;0 0 0 -1+q 0 0 -1 0 0 0 0 0 0 0 0;
      -q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 -q 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0;
      0 0 0 -q+q^2 0 0 1-q 0 -1 0 0 0 0 0 0;0 -q+q^2 0 0 0 q-q^2 0 0 0 -1 0 0 0 0 0;
      0 (-1+2q)-q^2 -q 1-q 1-q -q+q^2 1-q 0 0 q q 0 0 0 0;0 0 0 -1+q 0 0 -1+q^-1 0 0 0 0 -1 0 0 0;
      -j^2+j^2*q 1-q^-1 0 (1-2q)+q^2 -j^2+j^2*q -1+q 1-q -q -q 0 0 0 q 0 0;
      0 q-q^2 0 0 0 -q+q^2 0 0 0 0 0 0 0 -1 0;
      1-q^-1 (((2+q^-2)-2*q^-1)-2q)+q^2 0 (2-q^-1)-q 1-q^-1 (-2+q^-1+2q)-q^2 1-q^-1 0 0 0 0 q 0 q q],
     [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -1+q 0 -1 0 0 0 0 0 0 0 0 0 0;q -1 0 q 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -q 0 0 0 0 0 0 0 0 0 0 0 0;0 (1-2q)+q^2 0 0 0 (-1+2q)-q^2 0 0 0 -q 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 q 0 0 0 0;0 0 ((-j^2-2j)-q^-1)+j*q 0 -j+q^-2+(j^2+2j)*q^-1 0 0 -1 0 0 0 0 0 0 0;
      -j+(j^2+2j)*q+q^2 1-q 0 0 ((-j^2-2j)-q^-1)+j*q 0 -1+q q q 0 1-q 0 -1 0 0;
      0 (((3-q^-1)-4q)+3*q^2)-q^3 0 0 0 ((-2+3q)-3*q^2)+q^3 0 0 0 -q+q^2 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0 -1+q 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q;
      0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0;
      0 -4+q^-1+6q-4*q^2+q^3 0 0 0 2-5q+4*q^2-q^3 0 0 0 -1+2q-q^2 0 0 0 -1 0;
      0 0 0 0 0 0 0 0 0 0 0 -1 0 0 -1+q],
     [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 q^2 0 0 0 0 0 0 0 0 0 0 0;
      -j+j*q 0 0 0 -j^2-2j-q^-1+j*q 0 0 q 0 0 0 0 0 0 0;0 q^-1 0 -1+q 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 q 0 0 0 0 0 0 0 0;
      0 0 0 0 0 1 -1+q 0 0 0 0 0 0 0 0;-j+j*q 0 1 0 -j^2-2j-q^-1+j*q 0 0 -1+q 0 0 0 0 0 0 0;
      0 -1+q 0 0 0 1-q 0 0 -1+q -1 0 0 0 0 0;0 0 0 -q^2+q^3 0 0 q-q^2 0 -q 0 0 0 0 0 0;
      -j+j*q 0 0 q-q^2 -j^2-2j-q^-1+j*q 0 -1+q q q 0 0 0 -1 0 0;
      0 2-q^-1-q 0 0 0 -2+q^-1+q 0 0 0 0 0 -1+q 0 -1 0;
      0 1-2q+q^2 q 0 0 -1+2q-q^2 0 0 0 -q -q 0 -1+q 0 0;
      0 0 0 -q+2*q^2-q^3 0 0 1-2q+q^2 0 0 0 0 -q 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1],
     [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;q -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 -1+2q-q^2 0 0 0 1-2q+q^2 0 0 0 q 0 0 0 0 0;0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 q 0 0 0 0 0 0 0 0 0;0 0 0 0 1 -1+q 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0;0 2j^2+j-j^2*q^-1+q 0 0 0 -j^2+j^2*q 0 0 0 0 0 0 0 1 0;
      0 0 0 1-q 0 0 1-q^-1 0 0 0 0 1 0 0 0;1-2q+q^2 0 1 0 2-q^-1-q 0 0 0 0 -1+q 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0;0 0 0 q-q^2 0 0 -1+q 0 q 0 0 -1+q 0 0 0;
      0 1-2q+q^2 0 0 0 -1+2q-q^2 0 0 0 0 0 q 0 q q;
      j^2+(-2j^2-j)*q-q^2 0 0 0 j^2-j^2*q 0 0 q 0 0 0 0 0 -1+q 0;
      j-j*q 0 0 -q+q^2 j^2+2j+q^-1-j*q 0 1-q -q -q 0 0 0 1 0 -1+q]]
    end
    function f13(q,j)
      q^0*[[j-j*q -1+q-q^2 -j^2+2j^2*q-j^2*q^2 -j+(-j^2+j)*q+j^2*q^2 -j^2+(j^2-j)*q+j*q^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      -1 j^2-j^2*q -1+q 1-q -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -q q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 -q q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      -q 0 -q+q^2 0 -q+q^2 1-q 0 j^2-j^2*q 0 0 0 0 0 0 0 q 0 0 0 0;
      0 0 1-q -1+q 0 0 -1+q 0 j^2-j^2*q 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 q 0 0 0 0 0 0 0 0 0 0;
      0 0 -j^2+2j^2*q-j^2*q^2 -j+(-j^2+j)*q+j^2*q^2 0 0 j^2+(-j^2+j)*q-j*q^2 0 1-q+q^2 0 j-j*q 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 j^2 0 0 0 0 0 0;
      0 0 0 0 1-q -1+q 1-q 0 0 0 0 0 j^2-j^2*q 0 j^2*q 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 j*q 0 -1+q 0 0 0 0 0 0;
      0 0 0 0 2-q^-1-q 2j^2+j-j^2*q^-1+q -2j^2-j-q^-1+j^2*q 0 0 0 0 0 -j+j*q^-1+j*q 0 j-j*q 0 0 0 0 0;
      0 -1+q-q^2 -q+q^2 -j+(-j^2+j)*q+j^2*q^2 1-q j^2-j+j*q^-1-j^2*q 0 -1+q^-1+q 0 0 0 0 0 0 0 j-j*q 0 0 0 0;
      0 0 0 0 j^2-j-j^2*q^-1+j*q 0 -2j+j*q^-1+j*q 0 0 -j^2+j-j*q^-1+j^2*q 0 0 0 0 0 0 j-j*q 1-q^-1-q 0 0;
      0 0 0 0 -1+q 0 -1+q 0 0 1-q 0 0 0 0 0 0 -q j^2-j^2*q 0 0;
      0 j*q+(j^2-j)*q^2-j^2*q^3 -q+2*q^2-q^3 0 1-2q+q^2 0 1-q-q^2+q^3 -j+(-j^2+j)*q+j^2*q^2 j*q+(j^2-j)*q^2-j^2*q^3 0 -q+q^2 0 -j+(-j^2+j)*q+j^2*q^2 0 j^2*q-j^2*q^2 q-q^2 0 0 -1 0;
      1-q j+(j^2-j)*q-j^2*q^2 0 0 -1+q^-1-q+q^2 0 -1+q^-1-q+q^2 0 j+(j^2-j)*q-j^2*q^2 2-q^-1-q -1+q -j+j*q^-1 0 -1+q^-1 0 0 1-q j^2-j+j*q^-1-j^2*q 0 -1],
     [-1+q -j^2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      -j 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -q q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 -q q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      -j*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 j 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 q 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 j^2*q 0 -1+q 0 0 0 0 0 0 0 0 0;
      0 0 -1+q 0 0 0 1-q 0 0 -1+q 0 -1+q 0 1 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0 0;
      0 0 -1+q 0 0 0 0 0 0 q-q^2 0 q 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1+q 0 0 0 0 0;
      0 -j^2*q 0 0 0 0 0 j^2 0 0 0 0 0 0 0 -1+q 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1+q -j^2 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j*q 0 0 0;
      0 j-j*q 0 0 0 0 0 j-j*q^-1 j-j*q 0 j^2-j^2*q 0 j-j*q^-1 0 -j+j*q -j^2+j^2*q 0 0 -1 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1],
     [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 -1+q 0 0 j 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 1 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 j^2*q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 j^2 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 j^2*q 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 j*q 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 j 0 -1+q 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 q 0 -1+q 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0;
      0 0 -q+q^2 q-q^2 0 0 q-q^2 0 -j^2*q+j^2*q^2 0 -q 0 0 q 0 0 0 0 0 0;
      0 j+(j^2-j)*q-j^2*q^2 -1+2q-q^2 0 -2+q^-1+q 0 -1+q^-1-q+q^2 -j^2+j-j*q^-1+j^2*q j+(j^2-j)*q-j^2*q^2 0 -1+q 0 -j^2+j-j*q^-1+j^2*q 0 j^2-j^2*q 1-q 0 0 -1 0;
      -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 q 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j^2;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0;
      -1+q ((-j+3j*q)-3j*q^2)+j*q^3 ((j-3j*q)+3j*q^2)-j*q^3 (-q+2*q^2)-q^3 ((j^2+3j)-j*q^-1)+(-2j^2-4j)*q+(j^2+2j)*q^2 -j^2+(2j^2+j)*q+q^2 ((j^2+2j)-j*q^-1)+(-2j^2-j)*q+(j^2-j)*q^2+j*q^3 (((-2j^2-3j)-q^-1)+(j^2+3j)*q)-j*q^2 (-j+(j^2+3j)*q+(-2j^2-3j)*q^2)-q^3 0 -j^2+(j^2-j)*q+j*q^2 0 ((2-q^-1)-2q)+q^2 0 (-1+q)-q^2 -j+(-j^2+j)*q+j^2*q^2 0 0 j-j*q 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 -1+q],
     [q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0 0 0;
      0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 q 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 q 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 q 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 j^2*q 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -j^2 0 0;
      0 0 0 0 0 0 0 0 q 0 0 0 -1+q 0 0 0 0 0 0 0;
      0 0 0 0 1-q 0 1-q 0 0 -1+q 0 0 0 0 0 0 q -j^2+j^2*q 0 0;
      0 0 0 0 0 0 0 0 0 0 j 0 0 0 -1+q 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0;
      0 0 -1+q 0 0 0 1-q 0 0 -1+q 0 -1+q 0 1 0 0 -1+q 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 -j*q 0 0 0 0 0 -1+q 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
      0 -j*q+(-j^2+j)*q^2+j^2*q^3 (q-2*q^2)+q^3 0 (-1+2q)-q^2 0 (-1+q+q^2)-q^3 (j+(j^2-j)*q)-j^2*q^2 -j*q+(-j^2+j)*q^2+j^2*q^3 0 q-q^2 0 (j+(j^2-j)*q)-j^2*q^2 0 -j^2*q+j^2*q^2 -q+q^2 0 0 q q],
     [0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 q 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 q 0 -1+q 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 0 0 0 0 0 0;
      0 q 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 1 0 0 0 -1+q 0 0 0 0 0 0 0 0 0 0;
      -q 0 0 0 0 0 0 0 0 0 -1+q 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 j^2 0 0 0 -1+q 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -q 0 0;
      -q 0 -q+q^2 0 -q+q^2 1-q 0 j^2-j^2*q 0 0 0 0 0 -1+q 0 q 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j*q 0 0 0;
      0 0 -1+q 0 0 0 1-q 0 0 -1+q -1 -1+q 0 1 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 j^2 0 -1+q 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 -1+q 0 0;
      q-q^2 (j*q+(j^2-j)*q^2)-j^2*q^3 0 0 1-q-q^2+q^3 0 1-q-q^2+q^3 0 j*q+(j^2-j)*q^2-j^2*q^3 -1+2q-q^2 -q+q^2 j-j*q 0 1-q 0 0 q-q^2 j+(j^2-j)*q-j^2*q^2 0 -q;
      0 j+(j^2-j)*q-j^2*q^2 -1+2q-q^2 0 -2+q^-1+q 0 -1+q^-1-q+q^2 -j^2+j-j*q^-1+j^2*q j+(j^2-j)*q-j^2*q^2 0 -1+q 0 -j^2+j-j*q^-1+j^2*q 0 j^2-j^2*q 1-q 0 0 -1 -1+q]]
    end
    function f17(q)
    [[-1 0 0 0;0 -1 0 0;0 0 -1 0;-q^4 q^3 -q^2 q],
     [-1 0 0 0;0 -1 0 0;0 0 -1 0;-q^4 q^3 -q^2 q],
     [-1 0 0 0;0 -1 0 0;0 0 0 1;0 0 q -1+q],
     [-1 0 0 0;0 0 1 0;0 q -1+q 0;0 0 0 -1],
     [0 1 0 0;q -1+q 0 0;0 0 -1 0;0 0 0 -1]]
    end
    function f20(q,j)
    [[q 0 -q^2 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0;0 0 -1 0 0 0 0 0 0 0;0 0 0 -1 0 0 0 0 0 0;
      0 0 0 -q^2 q 0 0 0 0 0;0 -q^2 0 0 0 q 0 0 0 0;
      -j*q+(-j^2+j)*q^2+j^2*q^3 j^2*q+(-j^2+j)*q^2-j*q^3 -j^2*q^2+2j^2*q^3-j^2*q^4 j^2*q-2j^2*q^2+j^2*q^3 j+(j^2-j)*q-j^2*q^2 0 j-j*q -j+j*q-j*q^2 0 0;
      j^2*q-j^2*q^2 j^2*q-j^2*q^2 -j^2*q^2+j^2*q^3 j^2*q-j^2*q^2 -j^2+j^2*q 0 -j^2 j^2-j^2*q 0 0;
      -j*q^2+(-j^2+j)*q^3+j^2*q^4 (-j*q^2+2j*q^3)-j*q^4 (j^2*q^2+(-2j^2+j)*q^3+(2j^2-j)*q^4)-
      j^2*q^5 j^2*q^2+(-j^2+j)*q^3-j*q^4 0 j+(j^2-j)*q-j^2*q^2 0 -j*q+j*q^2-j*q^3 j-j*q j*q^-1-j+j*q;
      0 0 -j^2*q^3+j^2*q^4 0 -j^2*q^2+j^2*q^3 j^2*q-j^2*q^2 -j^2*q^2 0 j^2*q j^2-j^2*q],
     [q 0 -q^2 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0;0 0 -1 0 0 0 0 0 0 0;0 0 0 -1 0 0 0 0 0 0;
      0 0 0 -q^2 q 0 0 0 0 0;0 -q^2 0 0 0 q 0 0 0 0;0 0 0 0 0 0 -1+q -q 0 0;0 0 0 0 0 0 -1 0 0 0;
      0 0 0 0 0 0 0 -q^2 -1+q 1;0 0 0 0 0 0 -q^2 0 q 0],
     [-1+q 0 q 0 0 0 0 0 0 0;0 0 0 0 0 0 0 -1 0 0;1 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0;
      0 0 0 q -1+q 0 0 0 0 0;0 0 0 0 0 0 0 -q 0 q^-1;0 0 0 0 0 0 -1 0 0 0;0 -q 0 0 0 0 0 -1+q 0 0;
      0 0 0 0 0 0 -q^2 0 q 0;0 -q^3 0 0 0 q^2 0 0 0 -1+q],
     [q 0 -q^2 0 0 0 0 0 0 0;0 -1+q 0 q 0 0 0 0 0 0;0 0 -1 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0;
      0 0 0 0 0 1 0 0 0 0;0 0 0 0 q -1+q 0 0 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 q 0 -q^-1;
      0 0 0 0 0 0 q 0 -1+q 0;0 0 0 0 0 0 0 0 0 -1],
     [0 0 0 0 1 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0;0 0 q -1+q 0 0 0 0 0 0;
      q 0 0 0 -1+q 0 0 0 0 0;0 -q^2 0 0 0 q 0 0 0 0;0 0 0 0 0 0 -1 0 0 0;0 0 0 0 0 0 0 -1 0 0;
      0 0 0 0 0 0 -q^2 0 q 0;0 0 0 0 0 0 0 -q^2 0 q]]
    end
    function f23(q)
    [[-1 0 0 0 0;0 -1 0 0 0;q^3 -q^2 q 0 0;0 0 0 -1 0;-q^3 0 0 -q^2 q],
     [-1 0 0 0 0;0 -1 0 0 0;q^3 -q^2 q 0 0;0 0 0 -1 0;-q^3 0 0 -q^2 q],
     [-1 0 0 0 0;0 0 1 0 0;0 q -1+q 0 0;0 0 0 0 1;0 0 0 q -1+q],
     [0 1 0 0 0;q -1+q 0 0 0;0 0 -1 0 0;0 0 0 -1 0;-q^4 q^3 -q^2 -q^2 q],
     [-1 0 0 0 0;0 0 0 1 0;0 0 0 0 1;0 q 0 -1+q 0;0 0 q 0 -1+q]]
    end
    function f29(q)
    [[q q^4 0 0 -q^2 0;0 -1 0 0 0 0;0 0 -1 0 0 0;0 q^3 -q^2 q 0 0;0 0 0 0 -1 0;0 0 q^4 0 -q^3 q],
     [q q^4 0 0 -q^2 0;0 -1 0 0 0 0;0 0 -1 0 0 0;0 q^3 -q^2 q 0 0;0 0 0 0 -1 0;0 0 q^4 0 -q^3 q],
     [-1+q 0 0 0 q 0;0 -1 0 0 0 0;0 0 0 1 0 0;0 0 q -1+q 0 0;1 0 0 0 0 0;0 0 0 0 0 q],
     [0 0 0 0 0 1;0 0 1 0 0 0;0 q -1+q 0 0 0;0 0 0 -1 0 0;0 0 0 0 q 0;q 0 0 0 0 -1+q],
     [-1+q 0 0 q 0 0;0 q 0 0 0 0;0 0 0 0 1 0;1 0 0 0 0 0;0 0 q 0 -1+q 0;0 0 0 0 0 -1]]
    end
    x=-para[2][1]//para[2][2]
    ix=1//x
    if i==1 res=map(m->map(y->y(x),m),rep335_1)
    elseif i==2 res=f2(x)
    elseif i==3 res=f3(x,E(3))
    elseif i==4 res=f4(x,E(3))
    elseif i==5 res=f4(x,E(3,2))
    elseif i==6 res=f3(x,E(3,2))
    elseif i==7 res=[[-1;;],[-1;;],[-1;;],[-1;;],[-1;;]]
    elseif i==8 res=map(m->map(y->y(x),m),rep335_8)
    elseif i==9 res=map(m->map(y->y(x),m),rep335_9)
    elseif i==10 res=map(m->map(y->y(ix),m),rep335_8)*-x
    elseif i==11 res=f11(x,E(3))
    elseif i==12 res=f12(x,E(3))
    elseif i==13 res=f13(x,E(3))
    elseif i==14 res=f13(x,E(3,2))
    elseif i==15 res=f12(x,E(3,2))
    elseif i==16 res=f11(x,E(3,2))
    elseif i==17 res=f17(x)
    elseif i==18 res=map(m->map(y->y(ix),m),rep335_1)*-x
    elseif i==19 res=f13(ix,E(3,2))*-x
    elseif i==20 res=f20(x,E(3))
    elseif i==21 res=f20(x,E(3,2))
    elseif i==22 res=f13(ix,E(3))*-x
    elseif i==23 res=f23(x)
    elseif i==24 res=f2(ix)*-x
    elseif i==25 res=f11(ix,E(3,2))*-x
    elseif i==26 res=f12(ix,E(3))*-x
    elseif i==27 res=f12(ix,E(3,2))*-x
    elseif i==28 res=f11(ix,E(3))*-x
    elseif i==29 res=f29(x)
    elseif i==30 res=f4(ix,E(3,2))*-x
    elseif i==31 res=f4(ix,E(3))*-x
    elseif i==32 res=f23(ix)*-x
    elseif i==33 res=f3(ix,E(3,2))*-x
    elseif i==34 res=f3(ix,E(3))*-x
    elseif i==35 res=f17(ix)*-x
    elseif i==36 res=[[x;;],[x;;],[x;;],[x;;],[x;;]]
    end
    return -para[2][2]*res
  elseif (p,q,r)==(4,4,3)
    x,y=para[1]
    return [[[x y y 0 0 0;0 y 0 0 0 0;0 0 y 0 0 0;0 0 0 y 0 0;0 0 x+y y x 0;
       0 -y y y 0 x], 
      [y 0 0 0 0 0;-x x 0 0 0 x;0 0 0 0 -x 0;0 0 0 x 0 0;0 0 y 0 x+y 0;
       0 0 0 0 0 y], 
      [x y 0 0 0 y;0 y 0 0 0 0;0 0 x 0 -y y;0 -x 0 x y -x-y;
       0 0 0 0 y 0;0 0 0 0 0 y]], 
     [[y 0 0;0 0 -y;0 x x+y], [x 0 0;y y 0;-y 0 y], [0 -y 0;x x+y 0;0 0 y]], 
     [[x 0 0;x y 0;x 0 y], [y -2y 0;0 x 0;0 (1-E(4))x y], 
      [y 0 -y;0 y (-1-E(4))y//2;0 0 x]], 
     [[x 0 0;x y 0;x 0 y], [y -2y 0;0 x 0;0 (1+E(4))x y], 
      [y 0 -y;0 y (-1+E(4))y//2;0 0 x]], 
     [[y;;], [y;;], [y;;]], 
     [[x y 0;0 y 0;0 y x], [x+y 0 -y;0 x 0;x 0 0], [0 x 0;-y x+y 0;0 0 x]], 
     [[y 0 0;y x 0;y 0 x], [x -2x 0;0 y 0;0 (1+E(4))y x], 
      [x 0 -x;0 x (-1+E(4))x//2;0 0 y]], 
     [[y 0 0;y x 0;y 0 x], [x -2x 0;0 y 0;0 (1-E(4))y x], 
      [x 0 -x;0 x (-1-E(4))x//2;0 0 y]], 
     [[y 0;y x], [y 0;y x], [x -x;0 y]], 
     [[x;;], [x;;], [x;;]]][i]
  elseif q==1
    # Model of Ariki,  Halverson-Ram for reps of G(p,1,r): 
    # needs rational fractions if the parameters are indeterminates
    S=chevieget(:imp, :CharInfo)(p, q, r)[:charparams][i]
    if r>1 Q=-para[2][1]//para[2][2]
    else Q=0
    end
    function pos(t,i)
      for j in 1:length(t), k in 1:length(t[j])
        l=findfirst(==(i),t[j][k])
        if !isnothing(l) return [j,k,l] end
      end
    end
    ct(p)=para[1][p[1]]*(Q//1)^(p[3]-p[2])
    T=tableaux(S)
    return vcat([Diagonal(map(S->ct(pos(S,1)),T))], 
    map(2:r)do i
      toM(map(enumerate(T))do (j,t)
        a=pos(t,i);b=pos(t,i-1)
        t=map(l->map(copy,l),t)
        t[a[1]][a[2]][a[3]]=i-1;t[b[1]][b[2]][b[3]]=i #exchange i,i-1
        if para[2][1]==-para[2][2]
          if a[1]==b[1] tll=para[2][1]//(a[3]+b[2]-a[2]-b[3])
          else tll=0
          end
        else tll=sum(para[2])//(1-ct(b)//ct(a))
        end
        v=fill(0//1,length(T))*tll
        v[j]=tll
        pp=findfirst(==(t),T)
        if !isnothing(pp) v[pp]=tll-para[2][2] end
        v
       end)*prod(prod,para)^0
    end)
  else
    S=chevieget(:imp, :CharInfo)(p, q, r)[:charparams][i]
    if p==q para=[E.(p,0:p-1),para[1]]
    else e=div(p, q)
      if iseven(q) && r==2
        char=chevieget(:imp, :CharInfo)(p,q,r)[:malle][i]
        Z,X,Y=para
        if char[1]==1
          return [[Z[mod1(char[4],e)];;],[X[char[2]];;],[Y[char[3]];;]]
        else
          X=para[2]
          Y=para[3]
          if q>2
            Z=map(y->root(y,div(q,2)),para[1])
            Z=vcat(map(i->E(div(q,2),i)*Z,0:div(q,2)-1)...)
          else Z=para[1]
          end
          Z=Z[char[3:4]]
          if q==2 && div(p,q)==2 && !any(ismissing,rootpara)
            v=-char[2]*E(4)*rootpara[1]*rootpara[2]*rootpara[3]*
              root(E(div(p,2),2-char[3]-char[4]))*E(p,char[3]+char[4]-2)
          else v=char[2]*root(prod(Z)*prod(X)*prod(Y)*
                       E(div(p,2),2-char[3]-char[4]))*E(p,char[3]+char[4]-2)
          end
          return [[Z[1] sum(y->1//y,X)-Z[2]//v*sum(Y);0 Z[2]]^div(q,2),
             [sum(X) 1//Z[1];-prod(X)*Z[1] 0],[0 -prod(Y)//v;v sum(Y)]]
        end
      end
      # now para[2]==para[3]
      if para[1]==E.(e,0:e-1) para=[E.(p,0:p-1),para[2]]
      else # special-case q==2
        if q==2 
          if !ismissing(rootpara[1])
               Z=[rootpara[1],-E(4)*para[1][2]].//root(-para[1][2])
          else Z=[root(para[1][1]),root(para[1][2])]
          end
          Z=vcat(Z,-Z)
        else Z=[E(q,j)*root(i,q) for j in 0:q-1 for i in para[1]]
        end
        para=[Z,para[2]]
      end
    end
    extra=false
    if S[end] isa Integer
      extra=E(S[end-1],S[end])
      d=length(S)-2
      S=fullsymbol(S)
    end
    p1=length(para[1]);if p1!=p error() end
    v=chevieget(:imp,:HeckeRepresentation)(p1,1,r,para,[],
      findfirst(==(S),chevieget(:imp,:CharInfo)(p1,1,r)[:charparams]);gen=true)
    v*=1//1
    if p==q v[1]=inv(v[1])*v[2]*v[1]
    elseif q>1 v=vcat([v[1]^q,inv(v[1])*v[2]*v[1]],v[2:end]) end
    if extra==false return improve_type(v) end
    T=tableaux(S)
    m=orbits(Perm(T,map(S->S[vcat(d+1:p,1:d)],T)),1:length(T))
    l=map(i->extra^i,0:-1:1-p//d)
    m1=first.(m)
    improve_type(map(x->toM(map(c->transpose(x[c,m1])*l,m)),v))
  end
end)

chevieset(:imp, :Representation, function (p, q, r, i;gen=false)
  o=chevieget(:imp,:ordergens)(p,q,r)
  chevieget(:imp, :HeckeRepresentation)(p,q,r,map(x->
    x==2 ? [1,-1] : Cyc.(E.(x,0:x-1)),o),fill(1,r+1),i;gen)
end)

# the family of uc containing symbols
function MakeFamilyImprimitive(S, symbs)
  if length(S)==1 return Family(:C1,[findfirst(==(S[1]),symbs)]) end
  r=family_imprimitive(Symbols.entries(S[1]),length(S[1]))
  r.charNumbers=map(x->findfirst(==(x),symbs),r.symbols)
  r
end

simpexponent(l)=all(iszero,l[2:end]) ? l[1] : l

chevieset(:imp, :UnipotentCharacters, function (p, q, r)
  if !(q in [1,p]) return nothing end
  csy=q==1 ? symbols(p,r,1) : symbols(p,r,0)
  uc=Dict{Symbol, Any}(:charSymbols=>csy)
  uc[:a]=valuation_gendeg.(csy)
  uc[:A]=degree_gendeg.(csy)
  ci=chevieget(:imp, :CharInfo)(p,q,r)
  if q==1
    cusp=unique(sort(map(S->length.(S.S).-minimum(length.(S.S)),csy)))
    cusp=map(x->map(y->collect(0:y-1),x),cusp)
    sort!(cusp,by=x->rank(CharSymbol(x)))
    uc[:harishChandra]=map(cusp)do c
      cr=rank(CharSymbol(c))
      res=Dict{Symbol, Any}(:levi=>1:cr)
      res[:parameterExponents]=cr<r ? [simpexponent(length.(c))] : []
      append!(res[:parameterExponents],fill(1,max(0,r-cr-1)))
      if r==cr
        res[:relativeType]=TypeIrred(;series=:A,indices=Int[],rank=0)
      else
        res[:relativeType]=TypeIrred(;series=:ST,indices=1+cr:r,rank=r-cr,p,q=1)
      end
      res[:eigenvalue]=E(24,-2*(p^2-1)*div(sum(length,c),p))*
            E(2p,sum(i->-(i^2+p*i)*length(c[i+1]),0:p-1))
      res[:charNumbers]=map(x->
       findfirst(==(Symbol_partition_tuple(x,length.(c))),csy),
       map(x->partβ.(x.S),symbols(p,r-cr,1)[1:length(partition_tuples(r-cr,p))]))
      res[:cuspidalName]=ImprimitiveCuspidalName(c)
      res
    end
    uc[:b]=uc[:a]*0
    uc[:B]=uc[:a]*0
    uc[:b][uc[:harishChandra][1][:charNumbers]]=ci[:b]
    uc[:B][uc[:harishChandra][1][:charNumbers]]=ci[:B]
    uc[:families]=map(y->MakeFamilyImprimitive(y,uc[:charSymbols]),
      collectby(Symbols.entries,csy))
    sort!(uc[:families],by=charnumbers)
    if r==1
      l=map(function(S) # Dudas' sign change
        p=findfirst(==(Int[]),S.S)
        if isnothing(p) return 1 else return (-1)^p end
       end, csy[charnumbers(uc[:families][2])])
      uc[:families][2].signs=l
      uc[:families][2].fourierMat=Diagonal(l)*uc[:families][2].fourierMat*Diagonal(l)
      uc[:cyclicparam]=map(csy)do s
        l=findall(p->1 in p,s.S)
        length(l)==1 ? [1] : [l[2]-1,l[1]-1]
      end
    elseif (p,q,r)==(3,1,2) #Dudas' sign change
      uc[:families][4]=Family(uc[:families][4];signs=[-1,1,1])
      uc[:families][1]=Family(uc[:families][1];signs=[1,-1,-1,1,1,1,1,1,1])
    end
    return uc
  elseif p==q
    uc[:families]=[]
    for f in collectby(i->tally(vcat(csy[i].S...)),1:length(csy))
      if length(unique(map(x->x.S,csy[f])))>1
        push!(uc[:families],Dict{Symbol,Any}(:charNumbers=>f))
      else
        append!(uc[:families],map(x->Family(:C1,[x]),f))
      end
    end
    sort!(uc[:families],by=charnumbers)
    uc[:harishChandra]=map(l->Dict{Symbol,Any}(:charNumbers=>l),
      collectby(eachindex(csy))do i
        s=csy[i].S
        l=length.(s)
        (Symbols.relative_rank(csy[i]),l.-minimum(l))
    end)
    sort!(uc[:harishChandra],by=x->x[:charNumbers])
    extra=[]
    for f in uc[:harishChandra]
      addextra=false
      S=csy[f[:charNumbers][1]]
      l=r-Symbols.relative_rank(S)
      f[:levi]=1:l
      s=map(length,S.S)
      s.-=minimum(s)
      f[:eigenvalue]=E(24,-2*(p^2-1)*div(sum(s),p))*
           E(2p,sum(i->-(i^2+p*i)*s[i+1],0:p-1))
      if l==r
        f[:relativeType]=TypeIrred(;series=:A,indices=Int[],rank=0)
        f[:parameterExponents]=Int[]
        if length(f[:charNumbers])==2 addextra=true end
      elseif l == 0
        f[:relativeType]=TypeIrred(;series=:ST,indices=1:r,rank=r,p,q)
        f[:parameterExponents]=fill(1,r)
      else
        f[:relativeType]=TypeIrred(;series=:ST,indices=l+1:r,rank=r-l,p,q=1)
        f[:parameterExponents]=vcat([simpexponent(s)],fill(1,max(0,r-l-1)))
      end
      s=map(x->0:x-1,s)
      f[:cuspidalName]=ImprimitiveCuspidalName(s)
      if addextra
        s=copy(f[:charNumbers])
        f[:charNumbers]=[s[1]]
        f=copy(f)
        f[:charNumbers]=[s[2]]
        f[:cuspidalName]*='2'
        push!(extra, f)
      end
    end
    append!(uc[:harishChandra], extra)
    uc[:b]=fill(0,length(csy))
    uc[:b][uc[:harishChandra][1][:charNumbers]]=ci[:b]
    uc[:B]=fill(0,length(csy))
    uc[:B][uc[:harishChandra][1][:charNumbers]]=ci[:B]
    uc[:families]=map(f->MakeFamilyImprimitive(csy[f[:charNumbers]],
                uc[:charSymbols]),uc[:families])

    if (p,q,r)==(3,3,3) uc[:curtis]=[1,2,3,7,8,10,4,5,9,6,-12,-11] 
      uc[:families][4]=Family(uc[:families][4];signs=[-1,1,1])
      uc[:families][6]=Family(uc[:families][6];signs=[-1,1,1])
    elseif (p,q,r)==(3,3,4)
      uc[:families][2]=Family(uc[:families][2];signs=[-1,1,1])
      uc[:families][6]=Family(uc[:families][6];signs=[-1,-1,1,-1,1,1,-1,1,1])
      uc[:families][9]=Family(uc[:families][9];signs=[-1,1,1])
    elseif (p,q,r)==(3,3,5)
      uc[:families][3]=Family(uc[:families][3];signs=[-1,1,1])
      uc[:families][4]=Family(uc[:families][4];signs=[-1,1,1])
      uc[:families][6]=Family(uc[:families][6];signs=[1,-1,-1,-1,1,1,1,-1,1])
      uc[:families][7]=Family(uc[:families][7];signs=[-1,-1,1,-1,1,1,-1,1,1])
      uc[:families][8]=Family(uc[:families][8];signs=[-1,1,1])
      uc[:families][11]=Family(uc[:families][11];signs=[-1,1,1])
      uc[:families][13]=Family(uc[:families][13];signs=[-1,-1,1,-1,1,1,-1,1,1])
      uc[:families][15]=Family(uc[:families][15];signs=[-1,1,1])
      uc[:families][16]=Family(uc[:families][16];signs=[-1,1,1])
    end
    return uc
  end
end)
