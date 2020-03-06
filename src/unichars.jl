struct UniChar{T,T1}
  group::T
  v::T1
  prop::Dict{Symbol,Any}
end

IsUnipotentCharacter(u)=u isa UniChar

UniChar(W,v::Vector)=UniChar(W,v,Dict{Symbol,Any}())

function UniChar(W,v::Int)
  r=zeros(Int,length(UnipotentCharacters(W)))
  r[v] = 1
  UniChar(W,r)
end

function UniChar(W,v::String)
  n=charnames(stdout,UnipotentCharacters(W))
  UniChar(W,findfirst(==(v),n))
end

const short=Ref(true)

function Base.show(io::IO,r::UniChar)
  res=""
  s=charnames(io,UnipotentCharacters(r.group))
  m=maximum(length.(s))+3
  for i = 1:length(r.v)
    n = "<"*s[i]*">"
    c = sprint(show,r.v[i];context=io)
    if short[]
      if c != "0"
        if c == "1" res*= "+"
        elseif c == "-1" res*="-"
        else
          if '+' in c[2:end] || '-' in c[2:end] c = "("* c* ")" end
          if !(c[1] in "+-") res*="+" end
          res*=c
        end
        res*=n
      end
    elseif c!="0" || !get(io,:nozero,false)
      res *= "\n"* n* pad("", m - length(n))* c
    end
  end
  if length(res) == 0 res = "0" end
  if res[1] == '+' res = res[2:end] end
  if haskey(r.prop, :name)
    res="DLvar["*sprint(show,r.group; context=io)*","*r[:name],"]:",res
  else
    res="["*sprint(show,r.group; context=io)*"]:"* res
  end
  print(io,res)
end

Base.:+(u1::UniChar,u2::UniChar)=UniChar(u1.group,u1.v+u2.v)
Base.:-(u1::UniChar,u2::UniChar)=UniChar(u1.group,u1.v-u2.v)
Base.:*(u1::UniChar,u2::UniChar)=UniChar(u1.group,sum(u1.v .*u2.v))
Base.:*(u1::UniChar,a)=UniChar(u1.group,u1.v .* a)
Base.:*(a,u1::UniChar)=u1*a

Gapjm.degree(u::UniChar)=sum(u.v .*
                             degrees(UnipotentCharacters(u.group),Pol(:q)))

function LusztigInduction(WF, u)
  t = LusztigInductionTable(u[:group], WF)
  t==false ? false : UnipotentCharacter(WF, t[:scalar] * u.v)
end

LusztigRestriction(HF, u)=
  UnipotentCharacter(HF, u.v * (LusztigInductionTable(HF, u[:group]))[:scalar])

HCInduce(WF,u)=UnipotentCharacter(WF,HCInductionTable(u.group,WF)[:scalar]*u.v)

HCRestrict(HF,u)=
  UnipotentCharacter(HF,u.v*HCInductionTable(HF,u.group)[:scalar])

function DLCharTable(W)
  gets(W,:rwTable)do W
    uc=UnipotentCharacters(W)
    CharTable(W).irr'*fourier(uc)[uc.harishChandra[1][:charNumbers],:]
  end
end

DLChar(W,i::Int)=UniChar(W,DLCharTable(W)[i,:])

DLChar(W,w::Perm)=DLChar(W,position_class(W,w))

DLChar(W,w::Vector{Int})=DLChar(W,W(w...))

AlmostChar=function(W,i)
  ct=CharTable(W)
  dl=DLChar.(Ref(W),1:length(ct.charnames))
  sum(ct.irr[i,:] .* classes(ct).//length(W).*dl)
end

DLLefschetz=function (h,i=0)
  if haskey(h, :coset) W = ReflectionCoset(h[:coset])
  else W = Group(Hecke(h))
  end
  uc = UnipotentCharacters(W)
  UniChar(W, map(prod, conj.(HeckeCharValues(h)) * 
         fourier(uc)[uc[:almostHarishChandra][1][:charNumbers]], 
             map(x->x^i, Eigenvalues(uc))))
end

DLLefschetzTable=function(H)
  if haskey(H, :spets) WF = ReflectionCoset(H)
  else WF = Group(H)
  end
  t = CharTable(H).irr
  uc = UnipotentCharacters(WF)
  return t' * Fourier(uc)[uc[:almostHarishChandra][1][:charNumbers]]
end

Frobenius=function (WF, x, i)
  W = x[:group]
  p = map(x->PositionClass(W,Representative(x)^WF[:phi]), ConjugacyClasses(W))
  p = PermList(p)
  uc = UnipotentCharacters(W)
  t = DLCharTable(W)
  push!(t, Eigenvalues(uc))
  pt = permuted(t, p)
  t = TransposedMat(t)
  pt = TransposedMat(pt)
  p = map((x->begin Positions(t, x) end), pt)
  if any((x->begin length(x) > 1 end), p)
      error("Rw + eigen cannot disambiguate\n")
  end
  p = PermList(map(x->x[1], p))
  x = copy(x)
  x.v^=p^-i
  x
end
