# data common to several (but not all) types of reflection groups
# chevie/tbl/exceptio.jl (C) Jean Michel 1999-

chevieset([:A,:B,:D],:ordergens,t->r->fill(2,r))

chevieset([:E6,:E7,:E8,:F4,:G2,:H3,:H4],:ordergens,t->fill(2,string(t)[end]-'0'))

chevieset([:G25,:G26,:G29,:G31,:G32,:G34],:simplecoroots,function(t)
  eig=E.(chevieget(t,:ordergens))
  toM(coroot.(eachrow(chevieget(t,:simpleroots)),eig))
end)

chevieset([:G24,:G27,:G33,:G25,:G26,:G29,:G31,:G32,:G34],:CartanMat, t->
  convert(Matrix{Cyc{Int}},chevieget(t,:simplecoroots)*
         transpose(chevieget(t,:simpleroots))))

chevieset([:E7, :E8, :H3, :H4], :Invariants, t->function()
  C=chevieget(t, :CartanMat)
  r=toM(roots(C))*C
  map(d->function(arg...)sum((r*collect(arg)).^d) end, 
      chevieget(t, :ReflectionDegrees))
end)

chevieset([:G24,:G27,:G29,:G33,:G34,:E6,:E7,:E8,:H3,:H4], 
  :FactorizedSchurElement, t->function(phi,para,arg...)
               #arg=[] for G24--G34, [rootparam] for E6--H4
   i=findfirst(==(phi),chevieget(t,:CharInfo)()[:charparams])
   c=chevieget(t,:CycPolSchurElements)[i]
   q=-para[1][1]//para[1][2]
   res=HeckeAlgebras.FactSchur(Mvp(c[1]*q^Int(c[2])), 
                 map(v->(pol=CycPol([1,0,v]),monomial=q),c[3:length(c)]))
   HeckeAlgebras.simplify(res)
 end
)

chevieset([:G24,:G25,:G26,:G27,:G29,:G31,:G32,:G33,:G34,:H3,:H4,Symbol("2E6"),Symbol("2F4"),Symbol("3D4"),:E6,:E7,:E8,:F4,:G2],:IrredInfo,function(t)
  ci=chevieget(t,:CharInfo)()
  map((x,y)->(charparam=x,charname=y),ci[:charparams],ci[:charnames])
end)

chevieset([Symbol("3D4"),:E6,Symbol("2E6"),:E7,:E8,:F4,Symbol("2F4"),:G2,:H3,:H4],:CharTable,t->function()
  rank=string(t)[end]-'0'
  res=chevieget(t,:HeckeCharTable)(map(x->[1,-1],1:rank),map(x->1,1:rank))
  res[:identifier]=string("W(",t,")")
  res
end
)

chevieset([:G24,:G27,:G29,:G33],:CharTable,t->function()
  ngens=size(chevieget(t,:simpleroots),1)
  res=chevieget(t,:HeckeCharTable)(fill([1,-1],ngens),fill(1,ngens))
  res[:powermap]=chevieget(t,:PowerMaps)
  res[:galomorphisms]=chevieget(t,:galomorphisms)
  res[:identifier]=res[:name]=t
  res
end)

chevieset([:G24,:G27,:G29,:G33,:G34,:H3,:H4,:E6,:E7,:E8],:PoincarePolynomial,t->
function(q)
  prod(x->sum(y->(-q[1][1]//q[1][2])^y,0:x-1),chevieget(t,:ReflectionDegrees))
end)

chevieset([:E6,:E7,:E8,:F4,:G2,:H3,:H4,:G24,:G27,:G29,:G33],:Representation,
  function(t)
  ngens=size(chevieget(t,:simpleroots),1)
  i->chevieget(t,:HeckeRepresentation)(fill([1,-1],ngens),fill(1,ngens),i)
end)

chevieset(["3D4","2E6","2F4"],:Representation,t->
function(i)
  rk=t[end]-'0'
  chevieget(t,:HeckeRepresentation)(fill([1,-1],rk),fill(1,rk),i)
end)

chevieset([:G25,:G26],:Representation,t->
function(i)
  para=map(x->E.(x,0:x-1),chevieget(t,:ordergens))
  chevieget(t,:HeckeRepresentation)(para,[],i)
end)

chevieset([:G2,:F4,:H3,:E6,:G24,:G25,:G26,:G27,:G29,:G31,:G32,:G33,:G34],
  :SemisimpleRank,t->size(chevieget(t,:simpleroots),2))

chevieset([:A,:B,:D],:SemisimpleRank,t->(r->r))

# cycpolfakedegrees is a compact representation of fake degrees
chevieset([:H4, :E7, :E8, :G31], :FakeDegree,t->
function(phi, q)
  i=findfirst(==(phi),chevieget(t,:CharInfo)()[:charparams])
  f=chevieget(t,:cycpolfakedegrees)[i]
  res=f[1] isa AbstractVector ? evalpoly(q^2,f[1]) : f[1]
  f=copy(f)
  f[1]=1
  res*CycPol(f)(q)
end)

using Primes: totient
chevieset([:H4, :E7, :E8, :G31], :B,function(t)
  map(chevieget(t, :cycpolfakedegrees))do f
    res=f[1] isa AbstractVector ? 2*length(f[1])+f[2]-2 : f[2]
    res+sum(totient,f[3:length(f)];init=0)
  end
end)

chevieset([:H4, :E7, :E8, :G31], :b,
  t->map(f->f[2],chevieget(t,:cycpolfakedegrees)))

# sparseFakeDegrees is another compact representation of fake degrees
chevieset([Symbol("3D4"),:G2,:F4,Symbol("2F4"),:H3,:E6,:G24,:G25,:G26,:G27,:G29,:G32,:G33,:G34],:FakeDegree,t->
function(phi, q)
  i=findfirst(==(phi),chevieget(t,:CharInfo)()[:charparams])
  f=chevieget(t,:sparseFakeDegrees)[i]
  sum(i->f[i]*q^f[i+1],1:2:length(f)-1)
end)

chevieset([:E6,:G32,:G33,:G34,:G2,:F4,:H3,:G24,:G25,:G26,:G27,:G29],:B,
  t->map(x->x[end],chevieget(t,:sparseFakeDegrees)))

chevieset([:G2,:F4,:H3,:G24,:G25,:G26,:G27,:G29,:E6,:G32,:G33,:G34],:b,
  t->map(x->x[2],chevieget(t,:sparseFakeDegrees)))

chevieset([:G24,:G27,:G29,:G33,:G34,:H3,:H4,:E6,:E7,:E8],:A,function(t)
  N=sum(x->x-1,chevieget(t,:ReflectionDegrees))
  map(x->N-degree(CycPol(x)),chevieget(t,:CycPolSchurElements))
end)

chevieset([:G24,:G27,:G29,:G33,:G34,:H3,:H4,:E6,:E7,:E8],:a,
          t->map(x->-x[2],chevieget(t,:CycPolSchurElements)))

chevieset([:F4, :G2, :G25, :G26], :DecompositionMatrix,t->
function(p)
  T=chevieget(t,:CharTable)()
  T[:name]=T[:identifier]
  m=DecompositionMatrix(mod(T,p))
  map(c->[c[1],m[c[1]][c[2]]],blocks(m))
end)

chevieset([:G24,:G27,:G29,:G33,:G34,:E6,:E7,:E8,:H3,:H4],:SchurElement,t->
function(phi,para,arg...)
  i=findfirst(==(phi),chevieget(t,:CharInfo)()[:charparams])
  CycPol(chevieget(t,:CycPolSchurElements)[i])(-para[1][1]//para[1][2])
end)

chevieset([:G2,:F4,:G25,:G26,:G32],:FactorizedSchurElement,t->
function (phi,para,arg...)
  i=findfirst(==(phi),chevieget(t,:CharInfo)()[:charparams])
  Y=vcat(para[chevieget(t,:HyperplaneRepresentatives)]...)
  ci=chevieget(t,:SchurData)[i]
  VFactorSchurElement(Y,chevieget(t,:SchurModels)[Symbol(ci[:name])],ci)
end)

chevieset([:F4,:G25,:G26,:G32],:SchurElement,t->
function(phi,para,arg...)
  i=findfirst(==(phi),chevieget(t,:CharInfo)()[:charparams])
  Y=vcat(para[chevieget(t, :HyperplaneRepresentatives)]...)
  ci=chevieget(t,:SchurData)[i]
  VcycSchurElement(Y,chevieget(t,:SchurModels)[Symbol(ci[:name])],ci)
end)

chevieset([:E7,:E8,:F4,Symbol("2F4"),:G2,Symbol("3D4"),:H3,:H4,:G24,:G25,:G26,:G27,:G29,:G32,:G33,:G34],:Ennola,t->
function (arg...)
  uc = chevieget(t, :UnipotentCharacters)
  if uc isa Function uc = uc() end
  res = uc[:a] * 0
  for f = uc[:families]
    if haskey(f, :ennola)
      if f[:ennola] isa AbstractVector
        p=SPerm(f[:ennola])
      else
        A=Zbasedring(f)
        b=basis(A)
        if !haskey(f,:ennola) f[:ennola]=f[:special]
        end
        if f[:ennola]>0 p=SPerm(b[f[:ennola]]*b,b)
        else p=SPerm(-b[-f[:ennola]]*b,b)
        end
      end
    else
      p = SPerm()
    end
    res[f[:charNumbers]] = invpermute(f[:charNumbers], p)
  end
  SPerm(res)
end)
