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

spets=map(l2)do x
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
