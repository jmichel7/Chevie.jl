chevieset(:G31, :CharTable, function()
  res=chevieget(:G31, :HeckeCharTable)(map(x->[1,-1],1:4),[])
  res[:identifier]=res[:name]="G31"
  res[:galomorphisms] = Group(perm"(7,9)(8,12)(13,17)(15,16)(19,21)(20,23)(25,27)(26,28)(31,32)(35,37)(38,40)(42,45)(43,44)(46,49)(51,52)(54,55)(57,58)")
  res[:text]="origin: mostly CharTable(H(G31))"
  l=vcat(39:41, 50:53, 58:59)
  for i in l
    res[:irreducibles][i][[14, 19, 21, 35, 37, 41]]=[
       [1, E(4), -E(4), -E(4), E(4), -1], 
       [1, -E(4), E(4), E(4), -E(4), -1], 
       [0, 0, 0, 0, 0, 0], 
       [0, E(4)+1, -E(4)+1, E(4)-1, -E(4)-1, 0], 
       [0, -E(4)-1, E(4)-1, -E(4)+1, E(4)+1, 0], 
       [0, -E(4)+1, E(4)+1, -E(4)-1, E(4)-1, 0], 
       [0, E(4)-1, -E(4)-1, E(4)+1, -E(4)+1, 0], 
       [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]][findfirst(==(i),l)]
  end
  res
end)
chevieset(:G31, :SchurElement, function (p, para, rootpara)
  ci=findfirst(==(p),chevieget(:G31, :CharInfo)()[:charparams])
  data=chevieget(:G31,:SchurData)[ci]
  r=chevieget(:G31,:SchurModels)[Symbol(data[:name])]
  q=para[1][data[:order]]
  q=q[1]//q[2]
  if haskey(r,:root) q=root(q,2)*(-1)^data[:rootPower] end
  r[:coeff]*q^r[:factor]*prod(x->cyclotomic_polynomial(x)(q),r[:vcyc])
end)
