# Hand-translated part of chevie/tbl/weyla.g
# Jean Michel & Goetz Pfeiffer (C) 1994--2001
# Data for W(A_r)

CharTableSymmetric=Dict(:centralizers=>[
function(n,pp) res=k=1;last=0
  for p in pp
    res*=p
    if p==last k+=1;res*=k
    else k=1
    end
    last=p
  end
  res
end])

chevieset(:A,:CharTable,function(n)
  ct=chevieget(:imp,:CharTable)(1,1,n+1)
  ct[:irredinfo]=map(x->Dict(:charname=>joindigits(x)),chevieget(:A,:CharInfo)(n)[:charparams])
  ct
end)

chevieset(:A,:HeckeCharTable,
function(n,para,root)
  if n==1 Dict(:irreducibles=>[[1,para[1][2]],[1,para[1][1]]],
               :charnames=>["11","2"],
               :classnames=>["11","2"],
               :centralizers=>[2,2],
               :identifier=>"H(A_1)")
  else chevieget(:imp,:HeckeCharTable)(1,1,n+1,para,root)
  end
end)

chevieset(:A,:FakeDegree,(n,p,q)->fegsymbol([βset(p)])(q))
