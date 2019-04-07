struct UnipotentCharacters
  prop::Dict{Symbol,Any}
end

function params_and_names(sers)
  res=Dict{Symbol,Any}()
  chh=map(sers) do ser
    s=ser[:relativeType]
    s[:series]=Symbol(s[:series])
    if s[:rank]==0 return Dict(:charnames=>[""],:charparams=>[[]]) end
    charinfo(TypeIrred(s))
  end
  l=sum(x->length(x[:charnames]),chh)
  res[:charParams]=fill([],l)
  res[:TeXCharNames]=fill("",l)
  for i in eachindex(sers)
    ser=sers[i]
    t=ser[:relativeType]
    n=ser[:cuspidalName]
    ch=chh[i]
    res[:charParams][ser[:charNumbers]]=map(x->[n,x],ch[:charparams])
    res[:TeXCharNames][ser[:charNumbers]]=map(ch[:charnames])do x
#     s=(n isa String) ? n : prod(n)
      s=n
      if length(s)>0 && t[:rank]>0 s*=":" end
      if t[:rank]>0 s*=x end
      s
      end
  end
  res
end

function unipotent_characters(t::TypeIrred) 
  uc=getchev(t,:UnipotentCharacters)
  if uc==false 
    println("Warning: $t is not a Spets!!")
    return false 
  end
  merge!(uc,params_and_names(uc[:harishChandra]))
  if !haskey(uc,:charSymbols) uc[:charSymbols]=uc[:charParams] end
  uc
end

function unipotent_characters(W) 
  function CartesianSeries(sers)
    ser=Dict()
    ser[:levi]=vcat(map(x->x[:levi],sers))
    ser[:relativeType]=filter(x->x[:rank]!=0,map(x->x[:relativeType],sers))
    if haskey(sers[1],:eigenvalue)
      ser[:eigenvalue]=prod(map(x->x[:eigenvalue],sers))
    end
    if any(x->haskey(x,:qEigen),sers)
     ser[:qEigen]=sum(sers)do x
       if !haskey(x,:qEigen) return 0
       elseif x[:qEigen]==false return false
       else return x[:qEigen]
       end end
    else 
     ser[:qEigen]=0
    end
    if all(x->haskey(x,:parameterExponents),sers)
      ser[:parameterExponents]=vcat(map(x->x[:parameterExponents],sers))
    end
    ser[:charNumbers]=Cartesian(map(x->x[:charNumbers],sers)...)
    ser[:cuspidalName]=join(map(x->x[:cuspidalName],sers),"\\otimes ")
    ser
  end

  type=refltype(W)
  if isempty(type) # unipotent_characters(coxgroup())
    return UnipotentCharaters(Dict( 
      :harishChandra=>[
	rec(:relativeType=>[], 
	    :levi=>[], :parameterExponents=>[],
	    :cuspidalName=>"", :eigenvalue=>1, :charNumbers =>[ 1 ])],
      :families => [Family("C1",[1])],
      :charParams => [ [ "", [ 1 ] ] ],
      :TeXCharNames => [ "" ],
      :charSymbols => [ [ "", [ 1 ] ] ],
      :size=>1,
      :a => [ 0 ],
      :A => [ 0 ],
      :group=>W))
  end

  simp=map(type) do t
# adjust indices of Levis, almostLevis, relativetypes so they agree with
# Parent(Group(WF))
    uc=unipotent_characters(t)
    H=reflection_subgroup(W,t[:indices])
    for s in uc[:harishChandra]
     s[:levi]=inclusion(H)[s[:levi]]
     s[:relativeType][:indices]=inclusion(H)[s[:relativeType][:indices]]
    end
    uc
  end

  # "Kronecker product" of records in simp:
  r=simp[1]
  f=keys(r)
  res=Dict()
  for a in f
    if length(simp)==1 res[a]=map(x->[x],r[a])
    elseif all(x->haskeys(x,a),simp)
      res[a]=Cartesian(map(x->x[a],simp)...)
    end
  end
  
  for a in [:TeXCharNames]
    res[a]=map(x->join(x,"\\otimes "),res[a])
  end

  res[:size]=length(res[:TeXCharNames])
  
  for a in [:harishChandra] 
    res[a]=CartesianSeries.(res[a])
  end

  if length(type)==1
    res[:families]=map(res[:families]) do f
      f=f[1]
      f[:charNumbers]=map(x->[x],f[:charNumbers])
      f
    end
  else 
    res[:families]=prod.(res[:families])
  end
  
  for a in ["a", "A"]
    if haskey(res,a) res[a]=sum.(res[a]) end
  end

  # finally the new 'charNumbers' lists
  tmp=Cartesian(map(a->1:length(a[:TeXCharNames]),simp)...)
  for a in [ :harishChandra, :families]
    for s in res[a]
      s[:charNumbers]=map(y->findfirst(isequal(y),tmp),s[:charNumbers])
    end
  end

  res[:group]=W
  UnipotentCharacters(res)
end

function Base.show(io::IO,uc::UnipotentCharacters)
  println(io,"unipotent_characters(",uc.prop[:group],")")
  q=Pol([1],1)
  m=hcat(sprint.(show,CycPol.(Gapjm.degrees(uc,q)); context=io),
         sprint.(show,CycPol.(fakedegrees(uc,q)); context=io),
         sprint.(show,eigen(uc); context=io),
         TeXstrip.(labels(uc)))
  format(io,m,row_labels=TeXstrip.(uc.prop[:TeXCharNames]),
         rows_label=TeXstrip("\\gamma"),
         column_labels=TeXstrip.(["Deg(\\gamma)","Feg","Fr(\\gamma)","label"]))
end

group(uc::UnipotentCharacters)=uc.prop[:group]

function fakedegrees(uc::UnipotentCharacters,q)
  gets(uc,:fakedegrees)do uc
    fd=fill(zero(q),length(uc.prop[:TeXCharNames]))
    f=fakedegrees(group(uc),q)
    if isa(q,Pol) f=convert.(Pol{Int},f) end
    fd[uc.prop[:harishChandra][1][:charNumbers]]=f
    fd
  end
end

# FourierInverse times the vector of fake degrees is the vector of unip degrees
function fourierinverse(uc::UnipotentCharacters)
  gets(uc,:fourierinverse)do uc
    l=length(uc.prop[:TeXCharNames])
    i=one(fill(E(1)//1,l,l))
    for f in uc.prop[:families]
      i[f[:charNumbers],f[:charNumbers]]=f[:fourierMat]'
    end
    i
  end
end

function Gapjm.degrees(uc::UnipotentCharacters,q)
  gets(uc,:degrees)do uc
    fourierinverse(uc)*fakedegrees(uc,q)
  end
end

function eigen(uc::UnipotentCharacters)
  gets(uc,:eigen)do uc
    eig=fill(E(1),length(uc.prop[:TeXCharNames]))
    for f in uc.prop[:families] eig[f[:charNumbers]]=f[:eigenvalues]
    end
    eig
  end
end

function labels(uc::UnipotentCharacters)
  gets(uc,:labels)do uc
    lab=fill("",length(uc.prop[:TeXCharNames]))
    for f in uc.prop[:families] lab[f[:charNumbers]]=f[:charLabels]
    end
    lab
  end
end

# fix illegal relativeTypes B1 and C2 which appear in HC or almost HC
# series of classical groups
function FixRelativeType(t)
  if t[:relativeType][:series]=="B" 
   if t[:relativeType][:rank]==1
     t[:relativeType][:series]="A"
     t[:charNumbers]=collect(t[:charNumbers]) # map B1->A1
     t[:charNumbers][[1,2]]=t[:charNumbers][[2,1]] # map B1->A1
    elseif t[:relativeType][:rank]==2 && haskey(t[:relativeType],:cartanType) &&
      t[:relativeType][:cartanType]==1
      t[:relativeType][:cartanType]=2
      t[:relativeType][:indices]=reverse(t[:relativeType][:indices])
      t[:charNumbers][[1,5]]=t[:charNumbers][[5,1]] # map C2->B2
      if IsBound(t[:parameterExponents])
        t[:parameterExponents]=reverse(t[:parameterExponents])
      end
    end
  end
end
