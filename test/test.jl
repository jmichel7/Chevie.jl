nspets=map(i->ComplexReflectionGroup(i),
           [5,7,9,10,11,12,16,17,18,20,22,13,15,21])

l2=[
  (ComplexReflectionGroup,4),
  (ComplexReflectionGroup,6),
  (ComplexReflectionGroup,8),
  (ComplexReflectionGroup,14),
  (ComplexReflectionGroup,24),
  (ComplexReflectionGroup,25),
  (ComplexReflectionGroup,26),
  (ComplexReflectionGroup,27),
  (ComplexReflectionGroup,29),
  (ComplexReflectionGroup,32),
  (ComplexReflectionGroup,33),
  (ComplexReflectionGroup,3,1,2),
  (ComplexReflectionGroup,3,3,3),
  (ComplexReflectionGroup,3,3,4),
  (ComplexReflectionGroup,4,4,3),
  (coxgroup,:A,1),
  (coxgroup,:A,2),
  (coxgroup,:A,3),
  (coxgroup,:A,4),
  (coxgroup,:A,5),
  (coxgroup,:A,6),
  (coxgroup,:A,7),
  (coxgroup,:B,2),
  (coxgroup,:B,3),
  (coxgroup,:B,4),
  (coxgroup,:B,5),
  (coxgroup,:B,6),
  (coxgroup,:B,7),
  (coxgroup,:C,3),
  (coxgroup,:C,4),
  (coxgroup,:C,5),
  (coxgroup,:C,6),
  (coxgroup,:C,7),
  (coxgroup,:D,4),
  (coxgroup,:E,6),
  (coxgroup,:E,7),
  (coxgroup,:E,8),
  (coxgroup,:F,4),
  (coxgroup,:G,2),
  (coxgroup,:H,3),
  (coxgroup,:H,4),
  (coxgroup,:I,2,5),
  (coxgroup,:D,5),
  (coxgroup,:D,6),
  (coxgroup,:D,7),  #ch
  (ComplexReflectionGroup,34), #ch
  (coxgroup,),
]

l3=[
  (rootdatum,:psu,3),
  (rootdatum,Symbol("2B2")),
  (rootdatum,Symbol("2G2")),
  (rootdatum,Symbol("2I"),5),
  (rootdatum,Symbol("2I"),8),
  (rootdatum,:psu,4),
  (rootdatum,Symbol("3D4")),
  (rootdatum,:psu,5),
  (rootdatum,Symbol("pso-"),8),
  (rootdatum,:psu,6),
  (rootdatum,Symbol("2F4")),
  (rootdatum,:psu,7),
  (rootdatum,Symbol("pso-"),10),
  (rootdatum,Symbol("2E6")),
  (rootdatum,:psu,8),
  (rootdatum,Symbol("pso-"),12),
  (rootdatum,Symbol("pso-"),14),
]

spets_ex=map(l2)do x
  println("creating $x")
  x[1](x[2:end]...)
end

twisted=map(l3) do x
  println("creating $x")
  x[1](x[2:end]...)
end

function ct(l)
  for g in l
    println("creating CharTable($g)")
    show(IOContext(stdout,:limit=>true),CharTable(g))
  end
end

function uc(l)
  for g in l
    println("creating UnipotentCharacters($g)")
    show(IOContext(stdout,:limit=>true),UnipotentCharacters(g))
  end
end

function check_relations(H::HeckeAlgebra,t)
  W=H.W
  res=true
  for i in eachindex(gens(W))
    e=prod(q->(t[i]-q.*one(t[i])),H.para[i])
    if !iszero(e)
      println("#I  Error in ",Ordinal(i)," parameter relation");
      res=false
    end
  end
  for (l,r) in braid_relations(W)
    e=prod(t[l])-prod(t[r])
    if !iszero(e)
      println("#I Error in relation ",l,"=",r)
      res=false
    end
  end
  res
end

function FindRepresentation(W,gr,check=false)
  O=W
  if O isa HeckeAlgebra W=O.W end
  t=map(x->restriction(W,x),classinfo(W)[:classtext])
  l=1:length(t)
  l=sort(l,by=i->length(t[i]))
  t=t[l]
  ct=CharTable(O).irr
  pos=1:length(l)
  gens=gr
  for (i,w) in enumerate(t)
    r=traces_words_mats(gens,[w])[1]
#   println("w=$w r=$r")
#   if O isa HeckeAlgebra r=r*O.unit end
    pos=filter(j->iszero(ct[j,l[i]]-r),pos)
    if !check && length(pos)==1 return pos[1]
    elseif isempty(pos) return false
    end
  end
  if length(pos)==1 return pos[1]
  elseif length(pos)==0 return false
  else error("characters $pos are equal")
  end
end

function CheckRepresentations(W,l=Int[])
  O=W
  if W isa HeckeAlgebra
    H=W
    W=H.W
  else H=hecke(W)
  end
  cl=classinfo(W)[:classtext]
  ct=CharTable(O).irr
  if isempty(l) l=1:length(cl) end
  for i in l
    print("Representation #$i");
    gr=representation(O,i)
    if gr==false println("=false")
    else
      r=gr
      print(" dim:",size(r[1]),"...")
      check_relations(H,r)
      pos=FindRepresentation(O,gr,true)
      if pos==i println("")
      elseif !isnothing(pos) println("character found at ",pos)
      else println("character does not match")
	pos=TransposedMat([ct[i],traces_words_mats(gr,cl)])
	f=List(pos,x->x[1]!=x[2])
	pos=List(pos,x->List(x,FormatGAP))
	Print(FormatTable(ListBlist(pos,f),
	  rec(rowLabels=ListBlist(1:Length(pos),f),
	      columnLabels=["is","should be"])))
      end
    end
  end
end
