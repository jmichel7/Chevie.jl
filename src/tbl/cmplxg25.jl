#  tbl/cmplxg25.jl      CHEVIE library          Gunter Malle and Jean Michel
#  Copyright (C) 1998-  The CHEVIE Team

chevieset(:G25,:simpleroots,[0  0  -1;
                         -(2*E(3,2)+1)//3 -(2*E(3,2)+1)//3 -(2*E(3,2)+1)//3;
                         0 1 0])

chevieset(:G25,:ordergens,fill(3,3))

chevieset(:G25,:HyperplaneRepresentatives, [1])

chevieset(:G25,:BraidRelations,
          [[[1,2,1],[2,1,2]],[[1,3],[3,1]],[[2,3,2],[3,2,3]]])

chevieset(:G25,:ReflectionDegrees, [6, 9, 12])

chevieset(:G25,:NrConjugacyClasses, 24)

chevieset(:G25,:ParabolicRepresentatives,s->
  [[Int[]],[[1]],[[1,2],[1,3]],[1:3]][s+1])

# position in classes of G26
# [1,4,7,8,11,13,15,17,19,21,24,27,28,29,32,34,36,37,39,40,42,43,45,47]
chevieset(:G25, :ClassNames, [".", "cc", "31", "3131", "12231223", "1223", "d",
  "dd", "z", "zz", "2231223", "d1", "1", "131", "3221223221", "11", "1122",
  "12", "12z", "122312231223", "332112", "212", "c", "cz"])

chevieset(:G25, :WordsClassRepresentatives, map(x->
   collect(replace(x,"."=>"","z"=>"123"^4,"c"=>"123","d"=>"1232")).-'0',
   chevieget(:G25, :ClassNames)))

chevieset(:G25, :PowerMaps, [nothing, [1, 9, 4, 3, 15, 5, 8, 7, 10, 9, 3, 15, 16, 14, 5, 13, 16, 13, 4, 1, 10, 20, 2, 21], [1, 20, 1, 1, 1, 20, 9, 10, 1, 1, 20, 20, 1, 1, 1, 1, 20, 20, 20, 20, 20, 22, 22, 22], nothing, [1, 21, 4, 3, 15, 12, 8, 7, 10, 9, 19, 6, 16, 14, 5, 13, 18, 17, 11, 20, 2, 22, 24, 23], nothing, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24], nothing, nothing, nothing, [1, 21, 4, 3, 15, 12, 8, 7, 10, 9, 19, 6, 16, 14, 5, 13, 18, 17, 11, 20, 2, 22, 24, 23]])

chevieset(:G25, :ClassInfo, Dict{Symbol, Any}(
  :classtext=>chevieget(:G25,:WordsClassRepresentatives),
  :classnames=>chevieget(:G25,:ClassNames),
  :classparams=>chevieget(:G25,:ClassNames),
  :powermaps=>chevieget(:G25,:PowerMaps),
  :orders=>[1,6,3,3,3,6,9,9,3,3,6,6, 3, 3, 3, 3, 6, 6, 6, 2, 6, 4, 12, 12], 
 :classes=>[1,9,12,12,12,36,72,72,1,1,36,36,12,24,12,12,36,36,36,9,9,54,54,54]))

chevieset(:G25, :CharInfo, function ()
  res=Dict{Symbol,Any}(:charparams=>[[1,0],[1,24],[1,12],[2,15],[2,3],[2,9],
    [3,6],[3,5,2],[3,5,1],[3,17],[3,13,2],[3,1],[3,13,1],[6,8,2],[6,8,1],[6,2],
    [6,4,2],[6,10],[6,4,1],[8,3],[8,9],[8,6],[9,5],[9,7]],
# The labelling is determined as follows:
# phi{3,5}' is complexconjugate of phi{3,1}
# phi{3,13}' is complexconjugate of phi{3,17}
# phi{6,8}' is complexconjugate of phi{6,10}
# phi{6,4}' is complexconjugate of phi{6,2}
   :extRefl=>[1,12,8,3])
  res[:b]=map(x->x[2],res[:charparams])
  res[:charnames]=exceptioCharName.(res[:charparams])
  res
end)

chevieset(:G25, :HeckeCharTable, function (para, root)
  u,v,w=para[1]
  c=(u*v*w)^0
  res=Dict{Symbol, Any}(:name => "H(G25)", :identifier => "H(G25)",
    :parameter => para, :size => 648, :order => 648, :dim => 3,
    :degrees => [6, 9, 12], :reflclasses=>[13])
  f10(y)=map(w->y^length(w),res[:classtext])
  f23(u,v,w)=[2,-2*(u*v)^3,u^2+v^2,u^4+v^4,(u*v)^2*(u^4+v^4),-u*v*(u^2+v^2),
              -(u^2)*v^2,-(u^4)*v^4,2*u^6*v^6,2*u^12*v^12,-(v^3)*u^3*(u+v),
              -(v^2)*u^2*(u+v),u+v,(u+v)*((u^2-u*v)+v^2),u^4*v^4*(u^2+v^2),
              u^2+v^2,-u*v*(u^2+v^2),u*v,u^7*v^7,
              -(u^3)*v^3*(u^2+v^2)*((v^4-u^2*v^2)+u^4),-2*u^3*v^3,0,0,0]
  f31(u,v,w)=[3,-(u^4)*v^2,2*u*v+u^2,2*u^2*v^2+u^4,u^4*v^4+2*u^6*v^2,
              -(u^2)*v^2,0,0,3*u^8*v^4,3*u^16*v^8,-(u^4)*v^3,-(u^4)*v,2u+v,
              u*v^2+u^2*v+u^3,2*u^6*v^4+u^8*v^2,2*u^2+v^2,(-u*v^3-u^3*v)+u^4,
         u*v+u^2,u^9*v^5+u^10*v^4,-(u^6)*v^6,u^2*v^4-2*u^5*v,u^3,u^2*v,u^10*v^5]
  f62(u,v,w)=[6,2*u^3*v^2*w,2*u*v+2*u*w+u^2+v^2,2*u^2*v^2+2*u^2*w^2+u^4+v^4,
              u^2*(v^4*w^2+2*u^2*v^2*w^2+2*v^4*u^2+u^4*w^2),u*w*(u^2+v^2),0,0,
              6*u^6*v^4*w^2,6*u^12*v^8*w^4,u^3*v^2*w*(w+u),u^2*v^2*(w+u),
              3u+2v+w,v^2*u+u*w^2+v*u^2+u^2*w+u^3+v^3,
              v^2*u^4*(3*v^2*w^2+2*u^2*w^2+u^2*v^2),3*u^2+2*v^2+w^2,
              -((u^2+v^2))*((u*v-w^2)-u^2),u*(u+v),v^4*u^7*w^2*(u+v),
              u^3*w^3*(u^2+v^2)*((v^4-u^2*v^2)+u^4),
              v*(-2*u^3*w^2+3*u^4*v+v^3*w^2),-u*(-(u^2)+v*w),0,0]
  f83(u,v,w)=[8,0,2*(w+u)*(u+v),2*(w^2+u^2)*(u^2+v^2),
              2*u^2*v*w*(w^2+u^2)*(u^2+v^2),0,-(u^2)*v*w,-(u^4)*v^2*w^2,
              8*u^6*v^3*w^3,8*u^12*v^6*w^6,0,0,4u+2v+2w,
              v^2*u+u*w^2+v*w^2+v*u^2+u^2*w+v^2*w+2*u^3,
              2*u^4*v*w*(2*v^2*w^2+u^2*w^2+u^2*v^2),4*u^2+2*v^2+2*w^2,
              ((((-u*v^3-u*w^3)+u^2*v^2+u^2*w^2+v^2*w^2)-u^3*v)-u^3*w)+u^4,
              u*(u+v+w),v^3*u^7*w^3*(u+v+w),0,
              (-2*u^3*v^3-2*u^3*w^3)+v^3*w^3+3*u^4*v*w,-u*(-(u^2)+v*w),0,0]
  f97(u,v,w,J)=[9,-3*J^2*u^2*v^2*w^2,(u+v+w)^2,(u^2+w^2+v^2)^2,
    J*(u^2*v^2+u^2*w^2+v^2*w^2)^2,-(J^2)*(u^2*v^2+u^2*w^2+v^2*w^2),0,0,
    9*J*u^4*v^4*w^4,9*J^2*u^8*v^8*w^8,-(v^2)*J^2*u^2*w^2*(u+v+w),
    -v*u*w*J^2*(u*v+u*w+v*w),3u+3v+3w,(u+v+w)*(u^2+w^2+v^2),
    3*J*u^2*v^2*w^2*(u^2*v^2+u^2*w^2+v^2*w^2),3*u^2+3*v^2+3*w^2,
    (((((-u*v^3-u*w^3)-v*w^3)+u^2*v^2+u^2*w^2+v^2*w^2)-u^3*v)-u^3*w)-v^3*w,
    u*v+u*w+v*w,v^4*J*u^4*w^4*(u*v+u*w+v*w),(-(u^6)*v^6-u^6*w^6)-v^6*w^6,
    -J*u*v*w*(((2*w^3+2*v^3)-3*u*v*w)+2*u^3),-u*v*w,-J*u*v*w,-(J^2)*u^5*v^5*w^5]
  merge!(res, chevieget(:G25, :ClassInfo))
  merge!(res, chevieget(:G25, :CharInfo)())
  res[:centralizers]=div.(res[:order],res[:classes])
  # position in chars of G26
# [1,3,4,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47]
  res[:irreducibles]=toM([f10(u),f10(w),f10(v),f23(v,w,u),f23(u,v,w),f23(u,w,v),
  [3,3*u^2*v^2*w^2,u^2+v^2+w^2,u^4+v^4+w^4,u^4*v^4+u^4*w^4+v^4*w^4,
   u^2*v^2+u^2*w^2+v^2*w^2,0,0,3*u^4*v^4*w^4,3*u^8*v^8*w^8,
   u^2*v^2*w^3+u^2*v^3*w^2+u^3*v^2*w^2,u*v^2*w^2+u^2*v*w^2+u^2*v^2*w,
   u+v+w,u^3+v^3+w^3,u^2*v^4*w^4+u^4*v^2*w^4+u^4*v^4*w^2,u^2+v^2+w^2,u^2*v^2+
   u^2*w^2+v^2*w^2,0,0,u^6*v^6+u^6*w^6+v^6*w^6,3*u^2*v^2*w^2,
   -u*v*w,-u*v*w,-(u^5)*v^5*w^5],f31(v,u,w),f31(u,w,v),f31(w,v,u),f31(w,u,v),
  f31(u,v,w),f31(v,w,u),f62(w,u,v),f62(v,w,u),f62(u,v,w),f62(v,u,w),f62(w,v,u),
  f62(u,w,v),f83(u,v,w),f83(w,v,u),f83(v,u,w),
  f97(u,v,w,E(3,2)),f97(u,v,w,E(3))])*c
  res
end)

chevieset(:G25, :CharTable, function ()
  chevieget(:G25, :HeckeCharTable)([[1, E(3), E(3, 2)]], [])
end)


chevieset(:G25, :sparseFakeDegrees, [[1, 0], [1, 24], [1, 12], [1, 15, 1, 21],
  [1, 3, 1, 9], [1, 9, 1, 15], [1, 6, 1, 12, 1, 18], [1, 5, 1, 8, 1, 11],
  [1, 5, 1, 8, 1, 11], [1, 17, 1, 20, 1, 23], [1, 13, 1, 16, 1, 19], 
  [1, 1, 1, 4, 1, 7], [1, 13, 1, 16, 1, 19], [1, 8, 1, 11, 2, 14, 1, 17, 1, 20],
  [1, 8, 1, 11, 2, 14, 1, 17, 1, 20], [1, 2, 1, 5, 2, 8, 1, 11, 1, 14], 
  [1, 4, 1, 7, 2, 10, 1, 13, 1, 16], [1, 10, 1, 13, 2, 16, 1, 19, 1, 22], 
  [1, 4, 1, 7, 2, 10, 1, 13, 1, 16], [1, 3, 2, 6, 2, 9, 2, 12, 1, 15], 
  [1, 9, 2, 12, 2, 15, 2, 18, 1, 21], [1, 6, 2, 9, 2, 12, 2, 15, 1, 18], 
  [1, 5, 1, 8, 3, 11, 2, 14, 2, 17], [2, 7, 2, 10, 3, 13, 1, 16, 1, 19]])

chevieset(:G25, :SchurModels, Dict{Symbol, Any}(
  :f1_0=>Dict{Symbol,Any}(:vcyc=>[[[1,-1,0],1],[[1,-1,0],1],[[1,0,-1],1],
     [[1,0,-1],1],[[1,-1,0],4],[[1,0,-1],4],[[3,-2,-1],1],[[3,-1,-2],1],
     [[2,-1,-1],3],[[2,-1,-1],2],[[1,-1,0],6],[[1,0,-1],6]]),
  :f2_3=>Dict{Symbol,Any}(:vcyc=>[[[1,0,-1],1],[[1,0,-1],1],[[0,1,-1],1],
     [[0,1,-1],1],[[1,0,-1],2],[[0,1,-1],2],[[1,-1,0],1],[[1,-1,0],1],
     [[1,1,-2],3],[[1,1,-2],2],[[-1,1,0],6]]),
  :f3_1=>Dict{Symbol,Any}(:vcyc=>[[[1,-1,0],1],[[-1,1,0],1],[[1,0,-1],1],
     [[1,0,-1],1],[[1,0,-1],2],[[0,1,-1],1],[[1,1,-2],2],[[2,-1,-1],2],
     [[1,0,-1],6],[[1,-1,0],4],[[2,1,-3],1]]),
  :f3_6=>Dict{Symbol,Any}(:vcyc=>[[[1,-1,0],1],[[1,-1,0],1],[[1,0,-1],1],
     [[1,0,-1],1],[[0,1,-1],1],[[0,-1,1],1],[[-1,-1,2],2],[[-1,2,-1],2],[[-2,1,1],2]]),
  :f6_2=>Dict{Symbol,Any}(:vcyc=>[[[-1,1,0],1],[[1,0,-1],1],[[-1,0,1],1],
     [[0,1,-1],1],[[0,1,-1],1],[[1,0,-1],2],[[1,0,-1],6],[[1,-2,1],2],
     [[0,1,-1],2],[[3,-2,-1],1]]),
  :f8_3=>Dict{Symbol,Any}(:vcyc=>[[[0,1,-1],1],[[0,-1,1],1],[[-1,0,1],1],
         [[-1,1,0],1],[[2,-3,1],1],[[2,-1,-1],3],[[2,1,-3],1]]),
  :f9_7=>Dict{Symbol,Any}(:rootUnity=>E(3),:vcyc=>[[[0,0,0,1],1],[[0,0,0,2],2],
      [[-1,1,0],6],[[1,0,-1],6],[[0,-1,1],6],[[2,-1,-1,1],1],[[-1,2,-1,1],1],
      [[-1,-1,2,1],1]])))

chevieset(:G25, :SchurData, [
  Dict{Symbol, Any}(:name => "f1_0", :order => [1, 2, 3]), 
  Dict{Symbol, Any}(:name => "f1_0", :order => [3, 2, 1]), 
  Dict{Symbol, Any}(:name => "f1_0", :order => [2, 1, 3]), 
  Dict{Symbol, Any}(:name => "f2_3", :order => [2, 3, 1]), 
  Dict{Symbol, Any}(:name => "f2_3", :order => [1, 2, 3]), 
  Dict{Symbol, Any}(:name => "f2_3", :order => [1, 3, 2]), 
  Dict{Symbol, Any}(:name => "f3_6", :order => [1, 3, 2]), 
  Dict{Symbol, Any}(:name => "f3_1", :order => [2, 1, 3]), 
  Dict{Symbol, Any}(:name => "f3_1", :order => [1, 3, 2]), 
  Dict{Symbol, Any}(:name => "f3_1", :order => [3, 2, 1]), 
  Dict{Symbol, Any}(:name => "f3_1", :order => [3, 1, 2]), 
  Dict{Symbol, Any}(:name => "f3_1", :order => [1, 2, 3]), 
  Dict{Symbol, Any}(:name => "f3_1", :order => [2, 3, 1]), 
  Dict{Symbol, Any}(:name => "f6_2", :order => [3, 1, 2]), 
  Dict{Symbol, Any}(:name => "f6_2", :order => [2, 3, 1]), 
  Dict{Symbol, Any}(:name => "f6_2", :order => [1, 2, 3]), 
  Dict{Symbol, Any}(:name => "f6_2", :order => [2, 1, 3]), 
  Dict{Symbol, Any}(:name => "f6_2", :order => [3, 2, 1]), 
  Dict{Symbol, Any}(:name => "f6_2", :order => [1, 3, 2]), 
  Dict{Symbol, Any}(:name => "f8_3", :order => [1, 2, 3]), 
  Dict{Symbol, Any}(:name => "f8_3", :order => [3, 2, 1]), 
  Dict{Symbol, Any}(:name => "f8_3", :order => [2, 1, 3]), 
  Dict{Symbol, Any}(:name => "f9_7", :order => [1, 2, 3], :rootUnityPower=>1),
  Dict{Symbol, Any}(:name => "f9_7", :order => [1, 2, 3], :rootUnityPower=>2)])

chevieset(:G25, :HeckeRepresentation, function (para, root, i)
  f1=u->[[u;;], [u;;], [u;;]]
  f2(v,w)=WGraph2Representation([[[1, 3], [2]], [[1, 2, -1, v * w]]], [w, v])
  f31(u,v)=WGraph2Representation([[[1],[2],[3]],[[1,2,u,-v],[2,3,-v,u]]],[u,v])
  f32(u,v,w)=WGraph2Representation([[[[2],[]],[[],[1,2,3]],[[1,3],[]]],
                [[1,2,-1,u*w+v^2],[1,3,v,v],[2,3,-u*w-v^2,1]]],[u,v,w])
  f6(v,u,w)=WGraph2Representation([[[[2],[]],[[],[1,2]],[[1],[]],[[],[2,3]],
     [[3],[]],[[],[1,3]]],[[1,2,-1,v*w+u^2],[1,3,u,u],[1,4,-1,v*w+u^2],
     [1,5,-u,-u],[1,6,w,0],[2,3,-v*w-u^2,1],[2,6,-u*w,1],[4,5,v*w+u^2,-1],
     [4,6,-u*w,1]]],[v,u,w])
  f8(u,w,v)=WGraph2Representation([[[[2,3],[]],[[3],[1,2]],[[1,3],[]],[[2],[3]],
     [[1,3],[]],[[2],[1]],[[1],[2,3]],[[1,2],[]]],[[1,2,-u*v-w^2,1],[1,3,w,w],
      [1,4,v*w-w^2,0],[1,5,0,-1],[2,3,-1,u*v+w^2],[2,4,[1,0,3,w],-u],
      [2,5,0,-w],[2,6,-1,0],[3,6,[1,0,3,v-w],-u],[3,7,u*w+w^2,-1],[3,8,-w,-w],
      [4,5,-u,[1,v,3,0]],[4,7,0,v],[5,6,[1,0,3,1],-u*w],[5,7,-u,v-w],
    [5,8,0,v*w-w^2],[6,7,u*w,[1,-1,3,0]],[6,8,0,v-w],[7,8,-1,u*v+w^2]]],[u,w,v])
  f9(u,v,w,a)=WGraph2Representation([[[[2],[]],[[],[1,2,3]],[[1],[3]],
    [[1,3],[]],[[2],[1]],[[1],[2]],[[2],[3]],[[3],[2]],[[3],[1]]],
    [[1,2,-1,u*w+v^2],[1,3,v,[1,v,3,0]],[1,4,-a*v,0],[1,5,0,a^2*u-v],
     [1,6,0,a^2*u],[1,7,0,a^2*u-v],[1,8,0,-a^2*u],[1,9,v,[1,0,3,v]],
     [2,3,-u*w-v^2,1],[2,4,-u*w+a*v^2,0],[2,5,-a^2*v*w,0],[2,7,-a^2*v*w,0],
     [2,9,-u*w-v^2,1],[3,4,0,u+a^2*v],[3,5,0,u],[3,6,-a^2*w,a*v],[3,7,w,0],
     [4,5,[1,0,3,-w],u],[4,6,-a*w,0],[4,7,[1,-w,3,0],u],[4,8,a*w,0],
     [4,9,u+a^2*v,0],[5,6,-u,v],[5,9,0,w],[7,8,u,-v],[7,9,u,0],
     [8,9,-a*v,a^2*w]]],[u,v,w])
  u,v,w=para[1].+0
  if     i==1 f1(u)
  elseif i==2 f1(w)
  elseif i==3 f1(v)
  elseif i==4 f2(v,w) 
  elseif i==5 f2(u,v)
  elseif i==6 f2(u,w) 
  elseif i==7 f32(u,v,w)
  elseif i==8 f31(u,v)
  elseif i==9 f31(w,u)
  elseif i==10 f31(v,w) 
  elseif i==11 f31(u,w) 
  elseif i==12 f31(v,u) 
  elseif i==13 f31(w,v) 
  elseif i==14 f6(v,u,w) 
  elseif i==15 f6(u,w,v)
  elseif i==16 f6(w,v,u)  
  elseif i==17 f6(w,u,v)  
  elseif i==18 f6(u,v,w)  
  elseif i==19 f6(v,w,u)  
  elseif i==20 f8(u,v,w)  
  elseif i==21 f8(w,u,v)
  elseif i==22 f8(v,w,u)  
  elseif i==23 f9(u,v,w,E(3))  
  elseif i==24 f9(u,v,w,E(3,2))
  end
end)

let Z3=ImprimitiveCuspidalName([Int[],[0,1],[0,1]])
 chevieset(:G25, :UnipotentCharacters,
Dict{Symbol, Any}(:a => [0, 12, 12, 12, 2, 2, 4, 4, 1, 12, 4, 1, 12, 4, 8, 2,
 4, 8, 2, 2, 6, 6, 4, 4, 1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 6, 8, 12, 12,
 12, 12, 12], :harishChandra => Dict{Symbol,
 Any}[Dict(:parameterExponents => [1, 1, 1], :cuspidalName => "",
 :charNumbers => 1:24, :relativeType => Dict{Symbol, Any}(:series => "ST",
 :indices => 1:3, :ST => 25, :rank => 3), :eigenvalue => 1, :levi => Any[]),
 Dict(:parameterExponents => [1, 3], :cuspidalName => Z3,
 :charNumbers => [39, 31, 30, 41, 38, 40, 25, 27, 26],
 :relativeType => Dict{Symbol, Any}(:p => 3, :series => "ST", :indices => [3,
 2], :rank => 2, :q => 1), :eigenvalue => E(3,2), :levi => [1]),
 Dict(:parameterExponents => [[3, 3, 2, 0, 0, 2]],
  :cuspidalName => Z3*"\\otimes "*Z3, :charNumbers => [29, 28, 32, 44, 43, 33],
 :relativeType => Dict{Symbol, Any}(:p => 6, :series => "ST", :indices => [2],
 :rank => 1, :q => 1), :eigenvalue => E(3), :levi => [1, 3]),
 Dict(:parameterExponents => [[0, 4, 4]], :cuspidalName => "G_4",
 :charNumbers => [42, 34, 35], :relativeType => Dict{Symbol, Any}(:p => 3,
 :series => "ST", :indices => [3], :rank => 1, :q => 1), :eigenvalue => -1,
 :levi => 1:2), 
 mkcuspidal("G_{25}",36,-E(3)),
 mkcuspidal("G_{25}",37,E(3))],
 :families => Family[Family(:C1,[1]), 
   Family(Family(:X)(3),[12, 9, 25],signs=[1, 1, -1],ennola=-2),
   Family(Family(:QZ)(3, [perm"()", [E(3)]]),[16, 19, 20, 28, 26, 6, 27, 5, 29],
          signs=[1, 1, 1, 1, -1, 1, 1, 1, 1],ennola=4,cospecial=2), 
   Family(Family(:X)(6),[17,23,7,24,14, 32, 34, 30, 36, 8, 37, 31, 11, 35, 33],
          signs=[1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,1,-1],ennola=-15),
   Family(Family(:X)(3), [22, 21, 38],signs=[1, 1, -1],ennola=1),
   Family(Family(:X)(3),[15, 18, 39],signs=[1, 1, -1],ennola=-3),
   Family(SubFamilyij(Family(:ExtPowCyclic)(6,3),1,2,root(-2)),
        [3, 13, 40, 10, 41, 2, 43, 42, 4, 44],
        signs=[1, 1, 1, 1, -1, 1, -1, 1, -1, -1],ennola=-9,cospecial=6)],
 :A => [0,
 24, 24, 24, 16, 16, 20, 20, 11, 24, 20, 11, 24, 20, 22, 16, 20, 22, 16, 16,
 21, 21, 20, 20, 11, 16, 16, 16, 16, 20, 20, 20, 20, 20, 20, 20, 20, 21, 22,
 24, 24, 24, 24, 24]))
end

chevieset(:G25, :Invariants, [
  (x1,x2,x3)->-10*x1^3*x2^3-10*x1^3*x3^3-10*x2^3*x3^3+x1^6+x2^6+x3^6,
  (x1,x2,x3)->-x1^3*x2^6+x1^3*x3^6-x2^3*x3^6+x1^6*x2^3-x1^6*x3^3+x2^6*x3^3,
  (x1,x2,x3)->2*x1^3*x2^3*x3^6+2*x1^3*x2^6*x3^3+x1^3*x2^9+x1^3*x3^9+x2^3*x3^9+
  2*x1^6*x2^3*x3^3-4*x1^6*x2^6-4*x1^6*x3^6-4*x2^6*x3^6+x1^9*x2^3+x1^9*x3^3+
  x2^9*x3^3])

chevieset(:G25, :Discriminant, function ()
  (t1, t2, t3)->36*t1*t2^2*t3-t1^2*t3^2-32*t3^3+t1^3*t2^2+108*t2^4
end)
