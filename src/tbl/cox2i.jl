#  tbl/cox2i.g         CHEVIE library         Meinolf Geck and Jean Michel
#  Copyright (C) 1992 - 2001  The CHEVIE Team

chevieset("2I", :NrConjugacyClasses,m->div(m+3,2))

chevieset("2I", :WordsClassRepresentatives, function (m)
  r=[Int[]]
  x=[1]
  for i in 1:div(m+1,2)
    push!(r, copy(x))
    append!(x, [2, 1])
  end
  r
end)

chevieset("2I", :ClassInfo, function(m)
  res=Dict{Symbol,Any}(:classtext=>chevieget("2I",:WordsClassRepresentatives)(m))
  res[:classnames]=joindigits.(res[:classtext])
  res[:classparams]=res[:classnames]
  res[:classes]=[m]
  append!(res[:classes],fill(2,div(m,2)))
  if isodd(m) push!(res[:classes],1) end
  res[:orders]=map(i->2m//gcd(2m,length(i)),res[:classtext])
  res[:orders][1]=2
  res
end)

chevieset("2I", :PhiFactors,m->[1,-1])

chevieset("2I", :ClassParameter, function (m, w)
  if iseven(length(w)) return "" end
  w[1]==2 ? joindigits(vcat(w[2:end],1)) : joindigits(w)
end)

chevieset("2I", :CharInfo, function(m)
  res = Dict{Symbol, Any}(:extRefl => [1, 3, 2])
  if m==4
    res[:charparams] = [[[2], []], [[], [1, 1]], [[1], [1]]]
    res[:b] = [0, 4, 1]
    res[:B] = [0, 4, 3]
    res[:charnames] = map(string_partition_tuple, res[:charparams])
  else
    res[:charparams]=vcat([[1,0],[1,m]],map(i->[2,i],1:div(m-1,2)))
    res[:b]=map(x->x[2], res[:charparams])
    res[:B]=vcat([0,m],map(i->m-i,1:div(m-1,2)))
    res[:charnames]=exceptioCharName.(res[:charparams])
  end
  res
end)

chevieset("2I", :FakeDegree, function (m, phi, q)
  i=findfirst(==(phi),chevieget("2I",:CharInfo)(m)[:charparams])
  if i==1 q^0
  elseif i==2 q^m
  else q^(m+2-i)-q^(i-2)
  end
end)

chevieset("2I", :HeckeCharTable, function (m, param, rootpara)
  q=-param[1][1]//param[1][2]
  if m==4 ident = "2B"
  elseif m==6 ident = "2G"
  else ident="2I2"
  end
  ident=string(ident,"(",m,")")
  if q!=1 ident=string("H(", ident, ")") end
  v=ismissing(rootpara[1]) ? root(q) : rootpara[1]
  cl=chevieget("2I",:ClassInfo)(m)
  cos(i)=E(2m,i)+E(2m,-i)
  ct=map(i->map(j->cos(1)*0*v,cl[:classtext]), cl[:classtext])
  ct[1]=map(i->q^length(i),cl[:classtext])
  ct[2]=map(i->(-1)^length(i),cl[:classtext])
  for i in 1:div(m-1,2)
    for j in 1:div(m+1,2)
      ct[i+2][j+1]=-v^(2j-1)*cos(i*(2j-1))
    end
  end
  ci=chevieget("2I", :CharInfo)(m)
  tbl=Dict{Symbol, Any}(:identifier=>ident,:name=>ident,
    :cartan=>[2 -cos(1);-cos(1) 2], :size => 2m, :parameter => [q, q],
    :sqrtparameter=>[v,v],:irreducibles=>ct*v^0,
    :irredinfo=>map((x,y)->Dict{Symbol,Any}(:charparam=>x,:charname=>y),
                    ci[:charparams], ci[:charnames]))
  merge!(tbl, cl)
  tbl[:centralizers]=div.(tbl[:size],cl[:classes])
  AdjustHeckeCharTable(tbl, param)
end)

chevieset("2I", :CharTable,m->
          chevieget("2I", :HeckeCharTable)(m, fill([1, -1],2), [1,1]))

chevieset("2I", :Representation,(m, i)->
          chevieget("2I",:HeckeRepresentation)(m,fill([1,-1],2),[1,1],i))

chevieset("2I", :HeckeRepresentation, function (m, param, rootpara, i)
  q=-param[1][1]//param[1][2]
  v=ismissing(rootpara[1]) ? root(q) : rootpara[1]
  e=E(2m)
  if i==1 (gens=[[v^2;;],[v^2;;]],F=[1;;])
  elseif i==2 (gens=[[-1;;],[-1;;]]*v^0,F=[1;;])
  else i-=2 
    (gens=[[-1 0;v*(e^i+e^-i) v^2],
           [v^2 v*(e^i+e^-i);0 -1]]*v^0,F=-[0 1;1 0])
  end
end)

###########################################################################
#  Symbols, unipotent degrees, Fourier matrix, Frobenius eigenvalues
#  as in: G. Malle, "Unipotente Grade...", J. Algebra 177 (1995)
#  and: M. Geck, G.M. "Fourier transforms...", J. Algebra 260 (2003)
#  for non-trivial family of twisted dihedral group {^e}G(e,e,2)
#                                                               GM 05.12.01
# There are 3 families: Is, St and a big one M.
#  M has principal series almost characters R(0,j) where 0<j<e/2
#    with fake degree q^{e-j}-q^j.
#  M has cuspidal almost characters R(k,l) where 0<k<l<e-k
#  M has unipotent characters indexed by S(k,l) where 0<k<l<2e-k with k,l odd.
# Malle associates a symbol (s_0,...,s_{e-1}) where
#   s_i=[0] for i\ne 0,k,l,k+l, s_0=s_{k+l}=[0,1] and s_k=s_l=[]
# The principal series chars in M are S(0,l) where 0<l<e/2 with 
#   associated symbol s_i=[0] for i\ne 0,l and s_0=s_l=[1].
# The b is min(k+l,e-k-l) and the eigenvalue of F is E(e)^(k*l).
#
# Properties: S^-1=transpose(S)
# f.fakdeg:=vcat(map(i->x^(e-i[2])-x^i[2],ac),map(i->0,nc));
# f.unpdeg:=map(i->(c(i[1])-c(i[2]))/e*x*(x^2-1)*(x^e+1)/
#               prod(a->(x-z^a)*(x-z^-a),i),symUnp);
# S*f.unpdeg=f.fakdeg
# (Diagonal(f.eigenvalues)*transposed(S)*Diagonal(f.sh)^-1*S)^2=S^0
chevieset("2I", :UnipotentCharacters,function(e)
  uc=Dict{Symbol, Any}()
  n=div(e-1,2)
  nc=vcat(map(i->map(j->[i,j],i+1:e-i-1),1:n)...)
  ac=vcat(map(l->[0,l],1:n),nc)
  symUnp=vcat(map(i->map(j->2*[i,j].-1,i+1:e-i),1:n)...)
  if isodd(e)
   # each character is aligned with its Ennola-correspondent almost char.
    untUnp=map(function(s)
      res=map(x->mod(div(x,2),e),e.-reverse(s))
      if res[1]>res[2] res=map(x->mod(x,e), -res) end
      res
    end, symUnp)
    symUnp=symUnp[sortperm(untUnp)];sort!(untUnp)
  end
  if e==4 TeXpref="B_2"
  elseif e==6 TeXpref="G_2"
  else TeXpref=string("I_2(",e,")")
  end
  uc[:harishChandra]=vcat([
    Dict{Symbol, Any}(:relativeType=>TypeIrred(;series=:A,indices=[1],rank=1),
     :parameterExponents=>[e],:levi=>Int[],:eigenvalue=>1,:cuspidalName=>"",
     :charNumbers=>[2,1])], 
      map(symUnp)do x
        Dict{Symbol, Any}(:relativeType=>
          TypeIrred(;series=:A,indices=[],rank=0),:parameterExponents=>[],
          :levi=>[1,2],:eigenvalue=>E(2*e,prod(x)),
          :cuspidalName=>string("{}^2",TeXpref,"[",x[1],",",x[2],"]"),
          :charNumbers=>[2+findfirst(==(x),symUnp)])
      end)
  uc[:almostHarishChandra]=vcat([Dict{Symbol, Any}(:relativeType=>
    TypeIrred(;orbit=[TypeIrred(;series=:I,indices=[1,2],rank=2,bond=e)],
     twist=perm"(1,2)"),
     :parameterExponents=>[1,1],:levi=>Int[],:eigenvalue=>1,
     :cuspidalName=>"",:charNumbers=>1:n+2)], 
     map(nc)do x
       Dict{Symbol,Any}(:relativeType=>
         TypeIrred(;series=:A,indices=Int[],rank=0),:parameterExponents=>[],
         :levi=>[1,2],:eigenvalue=>E(e,-prod(x)),
         :cuspidalName=>string(TeXpref,"[",x[1],",",x[2],"]"),
         :charNumbers=>[n+2+findfirst(==(x),nc)])
     end)
  if e==4
    uc[:almostHarishChandra][2][:cuspidalName]="B_2"
    uc[:almostHarishChandra][1][:relativeType][:orbit][1]=
      TypeIrred(;series=:B,indices=[1,2],rank=2,cartanType=root(2))
  elseif e==6
    eig=[E(3,2),-1,E(3),1]
    uc[:almostHarishChandra][1][:relativeType][:orbit][1]=
      TypeIrred(;series=:G,indices=[1,2],rank=2,cartanType=root(3))
    for i in 1:4
      uc[:almostHarishChandra][i+1][:cuspidalName]=
        string("G2[",xrepr(eig[i];TeX=true),"]")
    end
  end
  uc[:charParams]=vcat(chevieget(:I, :CharInfo)(e)[:charparams][1:2],ac)
  uc[:almostCharSymbols]=vcat(CharSymbol.([map(x->[0],1:e),
                                           map(x->[0,1],1:e)]),map(ac)do s
    S=map(i->[0],1:e)
    k,l=s
    S[k+1]=Int[]
    S[l+1]=Int[]
    push!(S[1],1)
    push!(S[k+l+1],1)
    CharSymbol(S)
  end)
  uc[:almostCharSymbols][1].S[e]=[2]
  uc[:almostCharSymbols][2].S[e]=[1, 2]
  uc[:charSymbols]=vcat(CharSymbol.([map(x->[0],1:e),map(x->[0,1],1:e)]),
                        map(symUnp)do s
    k=div(s[1]+1,2)
    l=div(s[2]+1,2)
    S=map(function(i)
      if i==k+1 || i==l+1 return Int[]
      else return [0]
      end
    end, 1:e)
    push!(S[1],1)
    push!(S[k+l],1)
    CharSymbol(S)
  end)
  uc[:charSymbols][1].S[[1,2]]=[[0,2],Int[]]
  uc[:charSymbols][2].S[[1,2]]=[[0,1,2],[1]]
  c(a)=E(2*e,a)+E(2*e,-a)
  f=Family(Dict{Symbol,Any}(:eigenvalues=>map(s->E(2*e,s[1]*s[2]),symUnp),
     :fourierMat=>toM(map(j->map(i->(c(i'*reverse(j))-c(i'*j))//e,symUnp),ac)),
     :sh=>map(s->E(e,-s[1]*s[2]),ac),
     :charNumbers=>collect(1:length(ac)).+2,
     :name=>string("?",2+length(ac)),:special=>1))
  f.printname=f.name
  uc[:families]=[Family("C1",[1]),Family("C1",[2]),f]
  uc[:a]=vcat([0,e],map(x->1,ac))
  uc[:A]=vcat([0,e],map(x->e-1,ac))
  if e==5
# Modified 25-8-2004 to fit with H3, H4
# Asterisque, Geck-Malle, H4 in Duke are like the old version
# 'Unipotente Grade' and I2, imprimitive, current H3 H4 agree with new version
    uc[:families][3]=galois(uc[:families][3],13)
    for c in uc[:harishChandra]
      c[:eigenvalue]=galois(c[:eigenvalue],13)
    end
  end
  uc
end)

chevieset("2I", :Ennola, function(e)
  if isodd(e) return SPerm() end
  uc=chevieget("2I",:UnipotentCharacters)(e)
  l=uc[:charSymbols]
  SPerm(map(1:length(l))do i
    s=ennola(l[i])
    u=map(x->CharSymbol(x,s.repeats,s.no),circshift.(Ref(s.S),length(s):-1:1))
    for a in u
      p=findfirst(==(a),l)
      if !isnothing(p) return p*(-1)^uc[:A][i] end
    end
  end)
end)
