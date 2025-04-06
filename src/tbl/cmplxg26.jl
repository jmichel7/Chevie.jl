#  tbl/cmplxg26.jl       CHEVIE library          Gunter Malle and Jean Michel
#  Copyright (C) 1998-  The CHEVIE Team

chevieset(:G26,:simpleroots, 
          [0  1 -1; 
           0  0  1; 
           -E(4)//root(3) -E(4)//root(3) -E(4)//root(3)])

chevieset(:G26,:HyperplaneRepresentatives, [1, 2])

chevieset(:G26,:ordergens,[2,3,3])

chevieset(:G26,:BraidRelations,[[[1,2,1,2],[2,1,2,1]],[[1,3],[3,1]],
                                [[2,3,2],[3,2,3]]])

chevieset(:G26, :ReflectionDegrees, [6, 12, 18])

chevieset(:G26, :NrConjugacyClasses, 48)

chevieset(:G26, :ParabolicRepresentatives, function (s,)
  [[Int[]], [[1], [2]], [[1, 2], [1, 3], [2, 3]], [1:3]][s + 1]
end)

chevieset(:G26, :ClassNames, [".", "1", "212", "c3c3", "212c22c3", "12", "1212",
  "12121212", "c32", "1212z", "c32c32", "1212zzz", "12z", "c", "cc", "z", "zc",
  "zcc", "zz", "zzz", "zzzz", "zzzzz", "13", "13z", "13zz", "133", "c1223", "2",
  "21212", "2323c", "2z", "2zz", "2zzz", "22", "c12122", "3322", "23", "23z",
  "23zz", "232323", "232323z", "232323zz", "323", "c121", "c12", "c3c3c3",
  "323zzzz", "c3"])

chevieset(:G26, :WordsClassRepresentatives, map(x->
         collect(replace(x, "."=> "", "z"=> "123"^3, "c"=> "123")).-'0',
                chevieget(:G26, :ClassNames)))

chevieset(:G26, :PowerMaps, [nothing, [1, 1, 8, 19, 21, 7, 8, 7, 11, 28, 32, 8, 11, 15, 17, 19, 15, 17, 21, 1, 19, 21, 34, 7, 11, 28, 32, 34, 29, 29, 7, 11, 34, 28, 32, 34, 28, 32, 8, 1, 19, 21, 40, 42, 4, 40, 42, 4], [1, 2, 2, 40, 2, 2, 1, 1, 20, 20, 1, 20, 40, 16, 19, 20, 21, 22, 1, 20, 1, 20, 2, 40, 2, 2, 40, 1, 1, 20, 20, 1, 20, 1, 20, 40, 40, 2, 40, 40, 2, 40, 43, 46, 43, 46, 43, 46], nothing, [1, 2, 6, 42, 41, 3, 8, 7, 35, 33, 32, 31, 27, 18, 17, 22, 15, 14, 21, 20, 19, 16, 26, 39, 38, 23, 13, 34, 29, 30, 12, 11, 10, 28, 9, 37, 36, 25, 24, 40, 5, 4, 43, 48, 47, 46, 45, 44], nothing, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48], nothing, nothing, nothing, [1, 2, 6, 42, 41, 3, 8, 7, 35, 33, 32, 31, 27, 18, 17, 22, 15, 14, 21, 20, 19, 16, 26, 39, 38, 23, 13, 34, 29, 30, 12, 11, 10, 28, 9, 37, 36, 25, 24, 40, 5, 4, 43, 48, 47, 46, 45, 44], nothing, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48], nothing, nothing, nothing, [1, 2, 6, 42, 41, 3, 8, 7, 35, 33, 32, 31, 27, 18, 17, 22, 15, 14, 21, 20, 19, 16, 26, 39, 38, 23, 13, 34, 29, 30, 12, 11, 10, 28, 9, 37, 36, 25, 24, 40, 5, 4, 43, 48, 47, 46, 45, 44]])

chevieset(:G26, :ClassInfo, Dict{Symbol, Any}(
  :classtext => chevieget(:G26, :WordsClassRepresentatives),
  :classnames => chevieget(:G26, :ClassNames),
  :classparams => chevieget(:G26, :ClassNames),
  :orders=>[1,2,6,6,6,6,3,3,6,6,3,6,6,18,9,6,9,18,3,2,3,6,6,6,6,6,6,3,3,6,6,3,
            6,3,6,6,6,6,6,2,6,6,4,12,12,4,12,12],
  :classes=>[1,9,36,9,9,36,12,12,12,12,12,12,36,72,72,1,72,72,1,1,1,1,36,36,36,
           36,36,12,24,24,12,12,12,12,12,36,36,36,36,9,9,9,54,54,54,54,54,54]))

chevieset(:G26, :CharInfo, function ()
  res=Dict{Symbol,Any}(:charparams=>[[1,0],[1,9],[1,33],[1,21],[1,24],[1,12],
    [2,24],[2,15],[2,12],[2,3],[2,18],[2,9],[3,6],[3,15],[3,8,2],[3,5,2],
    [3,8,1],[3,5,1],[3,20],[3,17],[3,16,2],[3,13,2],[3,4],[3,1],[3,16,1],
    [3,13,1],[6,8,2],[6,11,2],[6,8,1],[6,11,1],[6,2],[6,5],[6,4,2],[6,7,2],
    [6,10],[6,13],[6,4,1],[6,7,1],[8,6,1],[8,3],[8,9,2],[8,12],[8,6,2],[8,9,1],
    [9,8],[9,5],[9,10],[9,7]],
# The labelling is as follows:
# The fakedegrees of phi{8,6}' and phi{8,9}' are monic.
# The complex conjugate of phi{3,8}' is phi{3,4}
# The complex conjugate of phi{3,5}' is phi{3,1}
# The complex conjugate of phi{3,16}' is phi{3,20}
# The complex conjugate of phi{3,13}' is phi{3,17}
# The complex conjugate of phi{6,8}' is phi{6,10}
# The complex conjugate of phi{6,11}' is phi{6,13}
# The complex conjugate of phi{6,4}' is phi{6,2}
# The complex conjugate of phi{6,7}' is phi{6,5}
   :hgal=>perm"(39,40)",:extRefl=>[1,24,15,4])
  res[:b]=map(x->x[2],res[:charparams])
  res[:charnames]=exceptioCharName.(res[:charparams])
  res
end)

chevieset(:G26, :HeckeCharTable, function (para, rt)
  c=prod(prod,para)^0
  res = Dict{Symbol, Any}(:size => 1296, :order => 1296, :identifier => "G26",
   :name => "G26", :powermap => chevieget(:G26, :PowerMaps), 
   :parameter=>para[[1,2]],:dim=>3,:irredinfo=>chevieget(:G26, :IrredInfo))
  merge!(res, chevieget(:G26, :ClassInfo))
  res[:centralizers]=div.(res[:order],res[:classes])
  f10(r,u)=[1,r,r*u^2,r^2*u^6,r^3*u^9,r*u,r^2*u^2,r^4*u^4,r*u^4,r^5*u^8,r^2*u^8,
    r^11*u^20,r^4*u^7,r*u^2,r^2*u^4,r^3*u^6,r^4*u^8,r^5*u^10,r^6*u^12,r^9*u^18,
    r^12*u^24,r^15*u^30,r*u,r^4*u^7,r^7*u^13,r*u^2,r^2*u^5,u,r^2*u^3,r*u^6,
    r^3*u^7,r^6*u^13,r^9*u^19,u^2,r^3*u^5,u^4,u^2,r^3*u^8,r^6*u^14,u^6,r^3*u^12,
    r^6*u^18,u^3,r^3*u^3,r^2*u^3,r^3*u^9,r^12*u^27,r*u^3]*c
  f23(r,p,u,v,w)=[2,2r,r*(u^2+v^2),-2*r^2*u^3*v^3,
    r^3*u^3*v^3*(u+v)*(v^2-u*v+u^2),r*(u+v),r^2*(u^2+v^2),r^4*(u^4+v^4),
    -r*u*v*(u^2+v^2),-(u^3)*v^3*r^5*(u^2+v^2),r^2*u^2*v^2*(u^4+v^4),
    -(u^9)*v^9*r^11*(u^2+v^2),-(u^3)*r^4*v^3*(u+v),r*u*v,-(r^2)*u^2*v^2,
    -2*r^3*u^3*v^3,-(r^4)*u^4*v^4,r^5*u^5*v^5,2*r^6*u^6*v^6,-2*r^9*u^9*v^9,
    2*r^12*u^12*v^12,-2*r^15*u^15*v^15,r*(u+v),-(u^3)*r^4*v^3*(u+v),
    u^6*r^7*v^6*(u+v),r*(u^2+v^2),-(u^2)*r^2*v^2*(u+v),u+v,
    r^2*(u+v)*((v^2-u*v)+u^2),-2*r*u^3*v^3,-(u^3)*v^3*r^3*(u+v),
    u^6*r^6*v^6*(u+v),-(u^9)*r^9*v^9*(u+v),u^2+v^2,-(u^2)*r^3*v^2*(u+v),
    -u*v*(u^2+v^2),u*v,-(r^3)*u^4*v^4,r^6*u^7*v^7,-2*u^3*v^3,2*r^3*u^6*v^6,
    -2*r^6*u^9*v^9,0,0,0,0,0,0]*c
  f36(r,p,u,v,w)=[3,3r,r*(u^2+v^2+w^2),3*r^2*u^2*v^2*w^2,
    r^3*u*v*w*(w^3*v^3+w^3*u^3+u^3*v^3),r*(u+v+w),r^2*(u^2+v^2+w^2),
    r^4*(u^4+v^4+w^4),r*(u^2*v^2+u^2*w^2+v^2*w^2),u^2*r^5*v^2*w^2*(u^2+v^2+w^2),
    r^2*(u^4*v^4+u^4*w^4+v^4*w^4),u^6*v^6*r^11*w^6*(u^2+v^2+w^2),
    u^2*r^4*v^2*w^2*(u+v+w),0,0,3*r^3*u^2*v^2*w^2,0,0,3*r^6*u^4*v^4*w^4,
    3*r^9*u^6*v^6*w^6,3*r^12*u^8*v^8*w^8,3*r^15*u^10*v^10*w^10,r*(u+v+w),
    u^2*r^4*v^2*w^2*(u+v+w),u^4*r^7*v^4*w^4*(u+v+w),r*(u^2+v^2+w^2),
    u*v*r^2*w*(v*w+u*v+u*w),u+v+w,r^2*(u^3+v^3+w^3),3*r*u^2*v^2*w^2,
    u^2*r^3*v^2*w^2*(u+v+w),u^4*r^6*v^4*w^4*(u+v+w),u^6*r^9*v^6*w^6*(u+v+w),
    u^2+v^2+w^2,u*v*r^3*w*(v*w+u*v+u*w),u^2*v^2+u^2*w^2+v^2*w^2,0,0,0,
    3*u^2*v^2*w^2,3*r^3*u^4*v^4*w^4,3*r^6*u^6*v^6*w^6,-u*v*w,-(r^3)*u*v*w,
    -(r^2)*u*v*w,-(r^3)*u^3*v^3*w^3,-(r^12)*u^9*v^9*w^9,-r*u*v*w]*c
  f31(r,p,u,v,w)=[3,p+2r,u*((-p*v-r*v)+r*u),u^4*v^2*r*(2p+r),p^2*r*u^5*v^4,u*r,
    u*r*(r*u-2*p*v),r^2*u^2*(2*p^2*v^2+r^2*u^2),-(u^2)*v*(-p*v+2*r*u),
    u^5*r^3*v^2*p*(r*u-2*p*v),u^4*v^2*(p^2*v^2+2*r^2*u^2),
    u^13*r^7*v^6*p^3*(r*u-2*p*v),p*r^3*u^5*v^2,0,0,3*p*r^2*u^4*v^2,0,0,
    3*p^2*r^4*u^8*v^4,3*p^3*r^6*u^12*v^6,3*p^4*r^8*u^16*v^8,
    3*p^5*r^10*u^20*v^10,p*u+r*u+r*v,u^4*r^2*v^2*p*(p*u+r*u+r*v),
    u^8*r^4*v^4*p^2*(p*u+r*u+r*v),p*u^2+r*u^2+r*v^2,u^3*r*v*(p*v+r*v+p*u),2u+v,
    r*u*((-p*v^2-p*u*v)+r*u^2),-(u^3)*v*((r*v^2-p*u*v)+r*u^2),
    u^4*r^2*v^2*p*(2u+v),u^8*r^4*v^4*p^2*(2u+v),u^12*p^3*r^6*v^6*(2u+v),
    2*u^2+v^2,u^3*p*r^2*v*(2v+u),u*((-(v^3)-u^2*v)+u^3),u*(u+v),
    u^5*r^2*v^2*p*(u+v),u^9*r^4*v^4*p^2*(u+v),u^3*(-2*v^3+u^3),
    u^7*v^2*r^2*p*(-2*v^3+u^3),u^11*v^4*r^4*p^2*(-2*v^3+u^3),u^3,p*r^2*u*v^2,
    -p*r*u^2*v,-(r^3)*u^6*v^3,p^4*r^8*u^19*v^8,-(u^2)*r*v]*c
  f62(r,p,u,v,w)=[6,2p+4r,(((-p*u*v-p*u*w)-r*u*v)-r*u*w)+r*u^2+r*v^2,
    -2*u^3*v^2*r*w*(2p+r),u^3*p^2*r*v*w^2*(u+v)*((v^2-u*v)+u^2),r*(u+v),
    r*((-2*p*u*v-2*p*u*w)+r*u^2+r*v^2),
    r^2*(2*p^2*u^2*v^2+2*p^2*u^2*w^2+u^4*r^2+v^4*r^2),
    -u*(((p*v^2*w-2*r*u*v*w)-2*u*r*v^2)+u^2*p*w),
    -(r^3)*u^3*v^2*p*w*((-2*p*u*v-2*p*u*w)+r*u^2+r*v^2),
    u^2*(p^2*v^4*w^2+2*r^2*u^2*v^2*w^2+2*u^2*r^2*v^4+u^4*p^2*w^2),
    -(u^9)*v^6*r^7*p^3*w^3*((-2*p*u*v-2*p*u*w)+r*u^2+r*v^2),
    -(u^3)*p*r^3*v^2*w*(u+v),0,0,-6*p*r^2*u^3*v^2*w,0,0,6*p^2*r^4*u^6*v^4*w^2,
    -6*p^3*r^6*u^9*v^6*w^3,6*p^4*r^8*u^12*v^8*w^4,-6*p^5*r^10*u^15*v^10*w^5,
    p*u+p*v+2*r*u+r*v+r*w,-(u^3)*p*r^2*v^2*w*(p*u+p*v+2*r*u+r*v+r*w),
    u^6*r^4*v^4*p^2*w^2*(p*u+p*v+2*r*u+r*v+r*w),p*u^2+p*v^2+2*r*u^2+r*v^2+r*w^2,
    -(u^2)*r*v*(v*r*w+p*u*v+2*v*p*w+r*u*w+p*u*w),3u+2v+w,
    r*((((-u*p*v^2-p*u*w^2)-p*u^2*v)-u^2*p*w)+r*u^3+r*v^3),
    u^2*v*(((v*r*w^2+r*v^2*w)-2*v*p*u*w)+r*u^2*w+u^2*r*v),
    -(u^3)*p*r^2*v^2*w*(3u+2v+w),u^6*r^4*v^4*p^2*w^2*(3u+2v+w),
    -(u^9)*p^3*r^6*v^6*w^3*(3u+2v+w),3*u^2+2*v^2+w^2,
    -(u^2)*p*r^2*v*(3*v*w+u*v+2*u*w),(u^2+v^2)*((u^2-u*v)+w^2),u*(u+v),
    -(u^4)*p*r^2*v^2*w*(u+v),u^7*r^4*v^4*p^2*w^2*(u+v),
    u^2*((3*v^2*w^2-2*u*v^3)+u^4),-p*r^2*u^5*v^2*w*((3*v^2*w^2-2*u*v^3)+u^4),
    p^2*r^4*u^8*v^4*w^2*((3*v^2*w^2-2*u*v^3)+u^4),u*(-v*w+u^2),
    v*r^2*p*(-v*w+u^2),0,0,u^13*v^8*p^4*r^8*w^4*(-v*w+u^2),0]*c
  function f83(r,p,u,v,w,eps)
    s=eps*root(-r*p*v*w)
    [8,4p+4r,(p+r)*(((-u*v-u*w)-v*w)+u^2),-4*u^3*v*s*w*(p+r),
    -(u^3)*v^3*w^3*(p+r)*((p^2-p*r)+r^2),u*(p+r),
    ((-2*p*r*u*v-2*p*r*u*w)-2*p*r*v*w)+p^2*u^2+r^2*u^2,
    2*p^2*r^2*u^2*v^2+2*p^2*r^2*u^2*w^2+2*p^2*r^2*v^2*w^2+p^4*u^4+r^4*u^4,
    ((s*u*((-(r^2)*v*w-p^2*v*w)+2*p*r*u*v+2*p*r*u*w+2*u^2*p*r))//p)//r,
    -(u^3)*v*p*r*s*w*(((-2*p*r*u*v-2*p*r*u*w)-2*p*r*v*w)+p^2*u^2+r^2*u^2),
    ((-(u^2)*v*w*(r^4*v^2*w^2+2*p^2*r^2*u^2*w^2+2*p^2*r^2*u^2*v^2+2*p^2*r^2*u^4+p^4*v^2*w^2))//p)//r,
    u^9*v^4*r^4*p^4*s*w^4*(((-2*p*r*u*v-2*p*r*u*w)-2*p*r*v*w)+p^2*u^2+r^2*u^2),
    -(u^4)*v*r*p*s*w*(p+r),-s*u,p*r*u^2*v*w,-8*p*r*s*u^3*v*w,
    -(p^2)*r^2*u^4*v^2*w^2,-(p^2)*r^2*s*u^5*v^2*w^2,-8*p^3*r^3*u^6*v^3*w^3,
    8*p^4*r^4*s*u^9*v^4*w^4,8*p^6*r^6*u^12*v^6*w^6,-8*p^7*r^7*s*u^15*v^7*w^7,
    (p+r)*(w+v+2u),-(u^3)*r*v*p*s*w*(p+r)*(w+v+2u),
    -(u^6)*v^3*r^3*p^3*w^3*(p+r)*(w+v+2u),(w^2+v^2+2*u^2)*(p+r),
    -(u^2)*(p+r)*s*(2*v*w+u*v+u*w),4u+2v+2w,
   -p*r*u*v^2-p*r*u*w^2-p*r*v*w^2-p*r*u^2*v-p*r*u^2*w-p*r*v^2*w+p^2*u^3+r^2*u^3,
    ((u^2*s*((-(r^2)*u*v*w-p^2*u*v*w)+p*r*v*w^2+p*r*v^2*w+p*r*u*v^2+p*r*u*w^2+p*r*u^2*v+p*r*u^2*w))//p)//r,
    -2*u^3*p*r*s*v*w*(w+v+2u),-2*u^6*p^3*r^3*v^3*w^3*(w+v+2u),
    2*u^9*p^4*r^4*s*v^4*w^4*(w+v+2u),4*u^2+2*v^2+2*w^2,
    -2*u^2*p*r*s*(2*v*w+u*v+u*w),
    ((((-u*v^3-u*w^3)+u^2*v^2+u^2*w^2+v^2*w^2)-u^3*v)-u^3*w)+u^4,u*(u+v+w),
    -(u^4)*p*r*s*v*w*(u+v+w),-(u^7)*p^3*r^3*v^3*w^3*(u+v+w),
    u^2*(((3*v^2*w^2-2*u*v^3)-2*u*w^3)+u^4),
    -p*r*s*u^5*v*w*(((3*v^2*w^2-2*u*v^3)-2*u*w^3)+u^4),
    -(p^3)*r^3*u^8*v^3*w^3*(((3*v^2*w^2-2*u*v^3)-2*u*w^3)+u^4),u*(-v*w+u^2),
    p*r*s*(-v*w+u^2),0,0,u^13*v^6*p^6*r^6*w^6*(-v*w+u^2),0]*c
      end
  f97(r,p,u,v,w,j)=[9,3p+6r,
    (((((-p*u*v-p*u*w)-v*p*w)-r*u*v)-r*u*w)-v*r*w)+r*u^2+r*v^2+r*w^2,
    3*r*u^2*v^2*j^2*w^2*(2p+r),r*u*v*p^2*j*w*(w^3*v^3+w^3*u^3+u^3*v^3),
    r*(u+v+w),r*(((-2*p*u*v-2*p*u*w)-2*v*p*w)+r*u^2+r*v^2+r*w^2),
    r^2*(2*p^2*u^2*v^2+2*p^2*u^2*w^2+2*p^2*v^2*w^2+u^4*r^2+v^4*r^2+r^2*w^4),
    j^2*((-2*r*u*v*w^2-2*r*u*v^2*w-2*r*u^2*v*w)+p*u^2*v^2+p*u^2*w^2+p*v^2*w^2),
    u^2*v^2*r^3*j^2*p*w^2*(((-2*p*u*v-2*p*u*w)-2*v*p*w)+r*u^2+r*v^2+r*w^2),
    j*(2*r^2*u^2*v^2*w^4+2*r^2*u^2*v^4*w^2+2*r^2*u^4*v^2*w^2+p^2*u^4*v^4+p^2*u^4*w^4+p^2*v^4*w^4),
    u^6*v^6*r^7*p^3*w^6*(((-2*p*u*v-2*p*u*w)-2*v*p*w)+r*u^2+r*v^2+r*w^2),
    u^2*j^2*p*r^3*v^2*w^2*(u+v+w),0,0,9*j^2*p*r^2*u^2*v^2*w^2,0,0,
    9*j*p^2*r^4*u^4*v^4*w^4,9*p^3*r^6*u^6*v^6*w^6,9*j^2*p^4*r^8*u^8*v^8*w^8,
    9*j*p^5*r^10*u^10*v^10*w^10,(p+2r)*(u+v+w),
    u^2*r^2*v^2*j^2*p*w^2*(p+2r)*(u+v+w),u^4*r^4*v^4*j*p^2*w^4*(p+2r)*(u+v+w),
    (u^2+v^2+w^2)*(p+2r),u*(2p+r)*r*v*j^2*w*(v*w+u*v+u*w),3u+3v+3w,
    r*(((-u*p*v^2-p*u*w^2-p*v*w^2-p*u^2*v-u^2*p*w)-p*v^2*w)+r*u^3+r*v^3+r*w^3),
    -u*v*j^2*w*(((v*r*w^2+w^2*r*u+r*v^2*w)-3*v*p*u*w)+r*u^2*w+u*r*v^2+u^2*r*v),
    3*u^2*j^2*p*r^2*v^2*w^2*(u+v+w),3*u^4*j*p^2*r^4*v^4*w^4*(u+v+w),
    3*u^6*p^3*r^6*v^6*w^6*(u+v+w),3*u^2+3*v^2+3*w^2,
    3*u*v*r^2*p*j^2*w*(v*w+u*v+u*w),
    (((((-u*v^3-u*w^3)-v*w^3)+u^2*v^2+u^2*w^2+v^2*w^2)-u^3*v)-u^3*w)-v^3*w,
    v*w+u*v+u*w,u^2*j^2*p*r^2*v^2*w^2*(v*w+u*v+u*w),
    u^4*j*p^2*r^4*v^4*w^4*(v*w+u*v+u*w),
    ((3*u^2*v^2*w^2-2*u^3*v^3)-2*w^3*u^3)-2*w^3*v^3,
    -(j^2)*p*r^2*u^2*v^2*w^2*(-3*u^2*v^2*w^2+2*u^3*v^3+2*w^3*u^3+2*w^3*v^3),
    -j*p^2*r^4*u^4*v^4*w^4*(-3*u^2*v^2*w^2+2*u^3*v^3+2*w^3*u^3+2*w^3*v^3),
    -u*v*w,-(j^2)*p*r^2*u*v*w,j*p*r*u*v*w,r^3*u^3*v^3*w^3,
    -(j^2)*p^4*r^8*u^9*v^9*w^9,j*r*u*v*w]*c
  r,p=para[1]
  u,v,w=para[2]
  res[:irreducibles]=toM([f10(r,u),f10(p,u),f10(p,w),f10(p,v),f10(r,w),f10(r,v),
    f23(p,r,v,w,u),f23(r,p,v,w,u),f23(p,r,u,v,w),f23(r,p,u,v,w),f23(p,r,u,w,v),
    f23(r,p,u,w,v),f36(r,p,u,v,w),f36(p,r,u,v,w),f31(p,r,v,u,w),f31(r,p,v,u,w),
    f31(p,r,u,w,v),f31(r,p,u,w,v),f31(p,r,w,v,u),f31(r,p,w,v,u),f31(p,r,w,u,v),
    f31(r,p,w,u,v),f31(p,r,u,v,w),f31(r,p,u,v,w),f31(p,r,v,w,u),f31(r,p,v,w,u),
    f62(r,p,w,u,v),f62(p,r,w,u,v),f62(r,p,v,w,u),f62(p,r,v,w,u),f62(r,p,u,v,w),
    f62(p,r,u,v,w),f62(r,p,v,u,w),f62(p,r,v,u,w),f62(r,p,w,v,u),f62(p,r,w,v,u),
    f62(r,p,u,w,v),f62(p,r,u,w,v),f83(r,p,u,v,w,1),f83(r,p,u,v,w,-1),
    f83(r,p,w,v,u,-1),f83(r,p,w,v,u,1),f83(r,p,v,u,w,1),f83(r,p,v,u,w,-1),
    f97(p,r,u,v,w,E(3,2)),f97(r,p,u,v,w,E(3,2)),f97(p,r,u,v,w,E(3)),
    f97(r,p,u,v,w,E(3))])
  res
end)

chevieset(:G26, :CharTable,()->
  chevieget(:G26, :HeckeCharTable)([[1,-1],[1,E(3),E(3,2)],[1,E(3),E(3,2)]],[]))

chevieset(:G26,:sparseFakeDegrees,[[1,0],[1,9],[1,33],[1,21],[1,24],[1,12],
  [1,24,1,30],[1,15,1,21],[1,12,1,18],[1,3,1,9],[1,18,1,24],[1,9,1,15],
  [1,6,1,12,1,18],[1,15,1,21,1,27],[1,8,1,14,1,20],[1,5,1,11,1,17],
  [1,8,1,14,1,20],[1,5,1,11,1,17],[1,20,1,26,1,32],[1,17,1,23,1,29],
  [1,16,1,22,1,28],[1,13,1,19,1,25],[1,4,1,10,1,16],[1,1,1,7,1,13],
  [1,16,1,22,1,28],[1,13,1,19,1,25],[1,8,2,14,2,20,1,26],[1,11,2,17,2,23,1,29],
  [1,8,2,14,2,20,1,26],[1,11,2,17,2,23,1,29],[1,2,2,8,2,14,1,20],
  [1,5,2,11,2,17,1,23],[1,4,2,10,2,16,1,22],[1,7,2,13,2,19,1,25],
  [1,10,2,16,2,22,1,28],[1,13,2,19,2,25,1,31],[1,4,2,10,2,16,1,22],
  [1,7,2,13,2,19,1,25],[2,6,3,12,2,18,1,24],[1,3,2,9,3,15,2,21],
  [1,9,2,15,3,21,2,27],[2,12,3,18,2,24,1,30],[1,6,2,12,3,18,2,24],
  [2,9,3,15,2,21,1,27],[1,8,3,14,3,20,2,26],[1,5,3,11,3,17,2,23],
  [2,10,3,16,3,22,1,28],[2,7,3,13,3,19,1,25]])

chevieset(:G26, :SchurModels, Dict{Symbol, Any}(
  :f1_0=>Dict{Symbol,Any}(:coeff=>-1,:vcyc=>[[[1,-1,0,0,0],1],[[0,0,1,-1,0],1],
    [[0,0,1,0,-1],1],[[1,-1,1,-1,0],2],[[1,-1,1,0,-1],2],[[1,-1,2,-2,0],1],
    [[1,-1,2,0,-2],1],[[1,-1,3,-2,-1],2],[[1,-1,3,-1,-2],2],[[1,-1,2,-1,-1],6],
    [[0,0,2,-1,-1],2],[[0,0,1,-1,0],6],[[0,0,1,0,-1],6]]),
  :f2_3=>Dict{Symbol,Any}(:factor=>[0,0,-1,1,0],:vcyc=>[[[1,-1,0,0,0],1],
    [[0,0,1,0,-1],1],[[0,0,0,1,-1],1],[[1,-1,1,0,-1],1],[[1,-1,0,1,-1],1],
    [[1,-1,1,0,-1],2],[[1,-1,0,1,-1],2],[[1,-1,1,-1,0],2],[[1,-1,-1,1,0],2],
    [[1,-1,1,1,-2],6],[[0,0,1,1,-2],2],[[0,0,1,-1,0],6]]),
  :f3_1=>Dict{Symbol,Any}(:coeff=>-1,:vcyc=>[[[-1,1,0,0,0],1],[[0,0,1,-1,0],1],
    [[0,0,1,0,-1],1],[[0,0,1,0,-1],2],[[0,0,0,1,-1],1],[[0,0,1,1,-2],2],
    [[0,0,2,-1,-1],2],[[0,0,1,0,-1],6],[[1,-1,1,0,-1],2],[[1,-1,-1,1,0],2],
    [[1,-1,2,-2,0],1],[[1,-1,2,1,-3],2]]),
  :f3_6=>Dict{Symbol,Any}(:coeff=>-1,:vcyc=>[[[1,-1,0,0,0],1],[[1,-1,0,0,0],3],
    [[1,-1,1,-1,0],2],[[1,-1,1,0,-1],2],[[1,-1,-1,1,0],2],[[1,-1,0,1,-1],2],
    [[1,-1,-1,0,1],2],[[1,-1,0,-1,1],2],[[0,0,1,1,-2],2],[[0,0,1,-2,1],2],
    [[0,0,-2,1,1],2]]),
  :f6_2=>Dict{Symbol,Any}(:vcyc=>[[[1,-1,0,0,0],1],[[0,0,-1,1,0],1],
    [[0,0,1,0,-1],1],[[0,0,0,1,-1],1],[[0,0,-1,0,1],2],[[0,0,1,0,-1],6],
    [[0,0,1,-2,1],2],[[1,-1,0,1,-1],1],[[-1,1,1,0,-1],2],[[1,-1,0,1,-1],2],
    [[1,-1,3,-2,-1],2]]),
  :f8_3=>Dict{Symbol,Any}(:coeff=>2,:root=>[1,1,0,1,1]//2,:rootCoeff=>-1,:vcyc=>
    [[[0,0,1,-1,0],1],[[0,0,1,0,-1],1],[[-1,1,0,-1,1],2],[[-1,1,0,1,-1],2],
    [[0,-1,1,0,-2,1],2],[[0,-1,1,-2,0,1],2],[[0,-1,-1,-1,1,1],1],
    [[0,-1,-1,1,-1,1],1],[[0,-1,1,-1,-1,1],3],[[-1,0,1,-1,-1,1],3]]),
  :f9_7=>Dict{Symbol,Any}(:rootUnity=>E(3),:vcyc=>[[[0,0,0,0,0,2],1],
    [[0,0,1,-1,0],6],[[0,0,-1,0,1],6],[[0,0,0,1,-1],6],[[1,-1,-2,1,1,1],2],
    [[1,-1,1,-2,1,1],2],[[1,-1,1,1,-2,1],2],[[-1,1,0,0,0],1],
    [[1,-1,0,0,0,1],1]])))

chevieset(:G26,:SchurData,[
  Dict{Symbol,Any}(:name=>"f1_0",:order=>[1,2,3,4,5]),
  Dict{Symbol,Any}(:name=>"f1_0",:order=>[2,1,3,4,5]),
  Dict{Symbol,Any}(:name=>"f1_0",:order=>[2,1,5,4,3]),
  Dict{Symbol,Any}(:name=>"f1_0",:order=>[2,1,4,3,5]),
  Dict{Symbol,Any}(:name=>"f1_0",:order=>[1,2,5,4,3]),
  Dict{Symbol,Any}(:name=>"f1_0",:order=>[1,2,4,3,5]),
  Dict{Symbol,Any}(:name=>"f2_3",:order=>[2,1,4,5,3]),
  Dict{Symbol,Any}(:name=>"f2_3",:order=>[1,2,4,5,3]),
  Dict{Symbol,Any}(:name=>"f2_3",:order=>[2,1,3,4,5]),
  Dict{Symbol,Any}(:name=>"f2_3",:order=>[1,2,3,4,5]),
  Dict{Symbol,Any}(:name=>"f2_3",:order=>[2,1,3,5,4]),
  Dict{Symbol,Any}(:name=>"f2_3",:order=>[1,2,3,5,4]),
  Dict{Symbol,Any}(:name=>"f3_6",:order=>[1,2,3,4,5]),
  Dict{Symbol,Any}(:name=>"f3_6",:order=>[2,1,3,4,5]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[2,1,4,3,5]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[1,2,4,3,5]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[2,1,3,5,4]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[1,2,3,5,4]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[2,1,5,4,3]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[1,2,5,4,3]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[2,1,5,3,4]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[1,2,5,3,4]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[2,1,3,4,5]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[1,2,3,4,5]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[2,1,4,5,3]),
  Dict{Symbol,Any}(:name=>"f3_1",:order=>[1,2,4,5,3]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[1,2,5,3,4]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[2,1,5,3,4]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[1,2,4,5,3]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[2,1,4,5,3]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[1,2,3,4,5]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[2,1,3,4,5]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[1,2,4,3,5]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[2,1,4,3,5]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[1,2,5,4,3]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[2,1,5,4,3]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[1,2,3,5,4]),
  Dict{Symbol,Any}(:name=>"f6_2",:order=>[2,1,3,5,4]),
  Dict{Symbol,Any}(:name=>"f8_3",:order=>[1,2,3,4,5],:rootPower=>-1),
  Dict{Symbol,Any}(:name=>"f8_3",:order=>[1,2,3,4,5],:rootPower=>1),
  Dict{Symbol,Any}(:name=>"f8_3",:order=>[1,2,5,4,3],:rootPower=>1),
  Dict{Symbol,Any}(:name=>"f8_3",:order=>[1,2,5,4,3],:rootPower=>-1),
  Dict{Symbol,Any}(:name=>"f8_3",:order=>[1,2,4,3,5],:rootPower=>-1),
  Dict{Symbol,Any}(:name=>"f8_3",:order=>[1,2,4,3,5],:rootPower=>1),
  Dict{Symbol,Any}(:name=>"f9_7",:order=>[2,1,3,4,5],:rootUnityPower=>2),
  Dict{Symbol,Any}(:name=>"f9_7",:order=>[1,2,3,4,5],:rootUnityPower=>2),
  Dict{Symbol,Any}(:name=>"f9_7",:order=>[2,1,3,4,5],:rootUnityPower=>1),
  Dict{Symbol,Any}(:name=>"f9_7",:order=>[1,2,3,4,5],:rootUnityPower=>1)])

chevieset(:G26, :HeckeRepresentation, function (para, rt, i)
  f10(x,u)=[[x;;],[u;;],[u;;]]
  f23(x,u,v)=[[x 0;0 x],[u 0;-u v],[v v;0 u]]
  f36(x,u,v,w)=[[x 0 0;0 x 0;0 0 x],[w 0 0;u*w+v^2 v 0;v 1 u],
                [u -1 v;0 v -u*w-v^2;0 0 w]]
  f31(x,y,u,v)=[[y 0 0;0 y -u*y-v*x;0 0 x],[v v 0;0 u 0;0 1 v],
                [u 0 0;-u v 0;0 0 v]]
  function f6(x,y,u,v,w)
    T=typeof(0+x*y*u*v*w)
    expandrep(3,6,Tuple{T, Vector{Int64}}[(u*w*x+w^2*y, [94]), (-u*w-v^2,
    [41]), (u*w+v^2, [6]), (-u*x-w*y, [97]), (u*x+w*y, [79]), (u, [3, 44]),
    (-v^2*x+w^2*y,[82]),(v^2*x+v*w*x,[100]), (-v^2-v*w, [47]), (v*w^2, [102]),
    (v*w, [12]), (v, [9, 23, 24, 38, 65, 105, 108]), (w, [2, 45, 66, 86, 87,
    107]), (x, [1, 22, 43, 64]), (y, [85, 106]), (-1, [20]), (1, [27, 53]),
    (-w^-1, [68, 71])])
  end
  function f8(x,y,u,v,w,sgn)
    s=sgn*root(-y*x*v*w)
    expandrep(3, 8, Tuple{typeof(s), Vector{Int64}}[(-s*u^2*v^-1*w^-1*y^-1-
     s*u^2*v^-2*y^-1+s*y^-1+u^2*v^-1*x*y^-1-u*x*y^-1+v,[120]),(s*u*v^-1*y^-1+
     s*u*v^-1*x^-1-u^2*v^-1-u+u*v^-1*w+w,[116]),(s*u*v^-1*w^-1+s*u*v^-2+s*v^-1-
     u*v^-1*x+u*v^-1*y+y,[70]),(s*u*v^-1*w^-1+u*v^-1*y,[22]),(s*x*y^-2-s*y^-1-
     s*u^-1*v*y^-1-s*u^-1*w*y^-1-s*u^-1*v^-1*w^2*y^-1-s*u^-2*w^2*y^-1+2w*x*y^-1
     +u^-1*v*w*x*y^-1-u^-1*v*w+u^-1*w^2*x*y^-1-u^-2*v*w^2,[101]),
     (-s*y^-1,[131]),(-s*y^-1-s*v^-1*w*y^-1-s*u^-1*w*y^-1+u*x*y^-1+w*x*y^-1
     -u^-1*v*w,[104]),(s*y^-1,[75]),(s*y^-1+v,[155]),(s*y^-1+w,[144]),
     (s*v^-1-y,[91]),(-s*u^-1*v*y^-1-u^-2*v*w^2,[173]),
     (-s*u^-1*w*y^-1+u^-1*v*w,[78]),(-u^2*v^-1*w^-1*y,[40]),(-u,[27]),
     (u,[2,32,57,84,96,110,111,137,162,165,168,191]),(-u*w^-1*y+u*v^-1*y,[64]),
     (u*v^-1*y+y,[67]),(-v,[179]),(-v-u^-1*w^2,[54]),(v,[3,56,83,192]),
     (v+u^-1*w^2,[5]),(-w,[140]),(w,[8,29,30,51,138,164]),(x,[109,136,163,190]),
     (-y,[61]),(y,[1,28,55,82]),(u^-1*v*w*x*y^-1+u^-2*v*w^2,[125]),
     (-u^-1*v*w,[176]),(u^-1*v*w,[128])])
  end
  function f9(x,y,u,v,w,j)
    T=typeof(0+x*y*u*v*w*j)
    expandrep(3,9,Tuple{T,Vector{Int64}}[
    (-j^2*u*w*x+j^2*v^2*y-j^2*v*w*y-j*u*v*y-v*w*x, [23]),
    (-j^2*u*x+j^2*v^2*w^-1*y-j^2*v*y-j*u*v*w^-1*y+j*v*x,[17]),(j^2*u*x,[72]),
    (-j^2*u*y,[13]),(-j^2*u-v,[227]),
    (-j^2*v^2*y+j^2*v*w*y+j^2*u^-1*v^2*w*y+u*w*y-w^2*x,[78]),
    (j^2*v^2*y-j*u*v*y-u*w*x,[73]),(-j^2*v*y-u*y,[46]),(j^2*v*y+u*y-w*x,
    [50]),(-j^2*w-u^-1*v*w,[231]),(-j^2*x,[125]),(-j^2-u^-1*v,[177]),(j^2,
    [224]),(-j^2*u^-1*v^2*y+j*v*y+w*x,[67]),(j^2*u^-1*v*w*y+w*y-u^-1*w^2*x,
    [105]),(j^2*u^-1*v*y,[40]),(-j^2*u^-1*x,[179]),(-j*u*y+u*x,[79]),(j*w,
    [148]),(j*u^-1*v^2*y-v*y+u^-1*v*w*x,[239]),(-u*w-v^2,[60]),(u*w+v^2,
    [5]),(-u,[119,208]),(u,[3,62,92,123,153,212]),(-u*w^-1,[207]),
    (v*w,[158,237]),(v,[8,32,33,57,113,152,183,213,242]),(-w*x,
    [100]),(-w,[87]),(w,[2,63,93,96,122,182,243]),(x,[106,121,181,
    241]),(y,[1,31,61,91,151,185,211]),(-1,[30,188]),(1,[35,116,
    202]),(u^-1*v,[167]),(u^-1*w*x,[94]),(u^-1,[170])])
  end
  x,y=para[1].+0
  u,v,w=para[2].+0
  if     i==1  f10(x,u)
  elseif i==2  f10(y,u)
  elseif i==3  f10(y,w)
  elseif i==4  f10(y,v)
  elseif i==5  f10(x,w)
  elseif i==6  f10(x,v)
  elseif i==7  f23(y,v,w)
  elseif i==8  f23(x,v,w)
  elseif i==9  f23(y,u,v)
  elseif i==10 f23(x,u,v)
  elseif i==11 f23(y,u,w)
  elseif i==12 f23(x,u,w)
  elseif i==13 f36(x,u,v,w)
  elseif i==14 f36(y,u,v,w)
  elseif i==15 f31(x,y,u,v)
  elseif i==16 f31(y,x,u,v)
  elseif i==17 f31(x,y,w,u)
  elseif i==18 f31(y,x,w,u)
  elseif i==19 f31(x,y,v,w)
  elseif i==20 f31(y,x,v,w)
  elseif i==21 f31(x,y,u,w)
  elseif i==22 f31(y,x,u,w)
  elseif i==23 f31(x,y,v,u)
  elseif i==24 f31(y,x,v,u)
  elseif i==25 f31(x,y,w,v)
  elseif i==26 f31(y,x,w,v)
  elseif i==27 f6(x,y,v,u,w)
  elseif i==28 f6(y,x,v,u,w)
  elseif i==29 f6(x,y,u,w,v)
  elseif i==30 f6(y,x,u,w,v)
  elseif i==31 f6(x,y,w,v,u)
  elseif i==32 f6(y,x,w,v,u)
  elseif i==33 f6(x,y,w,u,v)
  elseif i==34 f6(y,x,w,u,v)
  elseif i==35 f6(x,y,u,v,w)
  elseif i==36 f6(y,x,u,v,w)
  elseif i==37 f6(x,y,v,w,u)
  elseif i==38 f6(y,x,v,w,u)
  elseif i==39 f8(x,y,u,v,w,-1)
  elseif i==40 f8(x,y,u,v,w,1)
  elseif i==41 f8(x,y,w,v,u,1)
  elseif i==42 f8(x,y,w,v,u,-1)
  elseif i==43 f8(x,y,v,u,w,-1)
  elseif i==44 f8(x,y,v,u,w,1)
  elseif i==45 f9(x,y,u,v,w,E(3,2))
  elseif i==46 f9(y,x,u,v,w,E(3,2))
  elseif i==47 f9(x,y,u,v,w,E(3))
  elseif i==48 f9(y,x,u,v,w,E(3))
  end
end)

let j=E(3), j2=E(3,2), r3=root(-3)
CHEVIE[:families][:X49]=Dict{Symbol,Any}(:name=>"X49",:fourierMat=>toM([
  [-r3,r3,-9*j2,-9*j,9,9,9,(3-r3)*3,(3+r3)*3,((3-r3)*3)//2,((3+r3)*3)//2,((3-r3)*3)//2,((3+r3)*3)//2,r3,-r3,r3*2,-r3*2,(-3+r3)*3,(3+r3)*3,r3*6,((3+r3)*3)//2,((-3+r3)*3)//2,9*j,-9*j2,9,9*j2,9*j,-9*j,-9*j2,((-3-r3)*3)//2,((-3+r3)*3)//2,(3+r3)*3,(3-r3)*3,-9,-9,((3+r3)*3)//2,((-3+r3)*3)//2,9*j,-9*j2,((3+r3)*3)//2,((-3+r3)*3)//2,9*j,9*j2,r3*2,r3,r3,r3*6,r3*6,r3*6],
  [r3,-r3,-9*j,-9*j2,9,9,9,(3+r3)*3,(3-r3)*3,((3+r3)*3)//2,((3-r3)*3)//2,((3+r3)*3)//2,((3-r3)*3)//2,-r3,r3,-r3*2,r3*2,(-3-r3)*3,(3-r3)*3,-r3*6,((3-r3)*3)//2,((-3-r3)*3)//2,9*j2,-9*j,9,9*j,9*j2,-9*j2,-9*j,((-3+r3)*3)//2,((-3-r3)*3)//2,(3-r3)*3,(3+r3)*3,-9,-9,((3-r3)*3)//2,((-3-r3)*3)//2,9*j2,-9*j,((3-r3)*3)//2,((-3-r3)*3)//2,9*j2,9*j,-r3*2,-r3,-r3,-r3*6,-r3*6,-r3*6],
  [-9*j2,-9*j,9,9,-9,-9*j,-9*j2,0,0,-9,-9,9,9,-9,-9,0,0,0,0,0,-9*j,9*j2,-9*j,9*j2,9,9,9,9*j2,9*j,9*j2,9*j,0,0,-9*j,-9*j2,9*j,-9*j2,9*j,-9*j2,9*j2,-9*j,9*j2,9*j,0,9*j,-9*j2,0,0,0],
  [-9*j,-9*j2,9,9,-9,-9*j2,-9*j,0,0,-9,-9,9,9,-9,-9,0,0,0,0,0,-9*j2,9*j,-9*j2,9*j,9,9,9,9*j,9*j2,9*j,9*j2,0,0,-9*j2,-9*j,9*j2,-9*j,9*j2,-9*j,9*j,-9*j2,9*j,9*j2,0,9*j2,-9*j,0,0,0],
  [9,9,-9,-9,9,9,9,0,0,9,9,-9,-9,9,9,0,0,0,0,0,9,-9,9,-9,-9,-9,-9,-9,-9,-9,-9,0,0,9,9,-9,9,-9,9,-9,9,-9,-9,0,-9,9,0,0,0],
  [9,9,-9*j,-9*j2,9,9,9,0,0,9*j,9*j2,-9*j,-9*j2,9,9,0,0,0,0,0,9*j2,-9*j,9*j2,-9*j,-9,-9*j,-9*j2,-9*j2,-9*j,-9*j2,-9*j,0,0,9,9,-9*j2,9*j,-9*j2,9*j,-9*j2,9*j,-9*j2,-9*j,0,-9,9,0,0,0],
  [9,9,-9*j2,-9*j,9,9,9,0,0,9*j2,9*j,-9*j2,-9*j,9,9,0,0,0,0,0,9*j,-9*j2,9*j,-9*j2,-9,-9*j2,-9*j,-9*j,-9*j2,-9*j,-9*j2,0,0,9,9,-9*j,9*j2,-9*j,9*j2,-9*j,9*j2,-9*j,-9*j2,0,-9,9,0,0,0],
  [(3-r3)*3,(3+r3)*3,0,0,0,0,0,18,18,0,0,0,0,-r3*6,r3*6,(-3-r3)*3,(-3+r3)*3,-18*j2,18*j,0,0,0,0,0,0,0,0,0,0,0,0,18*j2,18*j,0,0,0,0,0,0,0,0,0,0,r3*6,(3+r3)*3,(-3+r3)*3,0,0,0],
  [(3+r3)*3,(3-r3)*3,0,0,0,0,0,18,18,0,0,0,0,r3*6,-r3*6,(-3+r3)*3,(-3-r3)*3,-18*j,18*j2,0,0,0,0,0,0,0,0,0,0,0,0,18*j,18*j2,0,0,0,0,0,0,0,0,0,0,-r3*6,(3-r3)*3,(-3-r3)*3,0,0,0],
  [((3-r3)*3)//2,((3+r3)*3)//2,-9,-9,9,9*j,9*j2,0,0,9,9,9,9,-r3*3,r3*3,(3+r3)*3,(3-r3)*3,0,0,0,9*j,-9*j2,9*j,-9*j2,9,9,9,-9*j2,-9*j,-9*j2,-9*j,0,0,-9*j,-9*j2,9*j,-9*j2,9*j,-9*j2,9*j2,-9*j,9*j2,9*j,-r3*6,((3+r3)*3)//2,((-3+r3)*3)//2,0,0,0],
  [((3+r3)*3)//2,((3-r3)*3)//2,-9,-9,9,9*j2,9*j,0,0,9,9,9,9,r3*3,-r3*3,(3-r3)*3,(3+r3)*3,0,0,0,9*j2,-9*j,9*j2,-9*j,9,9,9,-9*j,-9*j2,-9*j,-9*j2,0,0,-9*j2,-9*j,9*j2,-9*j,9*j2,-9*j,9*j,-9*j2,9*j,9*j2,r3*6,((3-r3)*3)//2,((-3-r3)*3)//2,0,0,0],
  [((3-r3)*3)//2,((3+r3)*3)//2,9,9,-9,-9*j,-9*j2,0,0,9,9,9,9,-r3*3,r3*3,(3+r3)*3,(3-r3)*3,0,0,0,9*j,-9*j2,-9*j,9*j2,-9,-9,-9,9*j2,9*j,-9*j2,-9*j,0,0,9*j,9*j2,9*j,-9*j2,-9*j,9*j2,9*j2,-9*j,-9*j2,-9*j,-r3*6,((3+r3)*3)//2,((-3+r3)*3)//2,0,0,0],
  [((3+r3)*3)//2,((3-r3)*3)//2,9,9,-9,-9*j2,-9*j,0,0,9,9,9,9,r3*3,-r3*3,(3-r3)*3,(3+r3)*3,0,0,0,9*j2,-9*j,-9*j2,9*j,-9,-9,-9,9*j,9*j2,-9*j,-9*j2,0,0,9*j2,9*j,9*j2,-9*j,-9*j2,9*j,9*j,-9*j2,-9*j,-9*j2,r3*6,((3-r3)*3)//2,((-3-r3)*3)//2,0,0,0],
  [r3,-r3,-9,-9,9,9,9,-r3*6,r3*6,-r3*3,r3*3,-r3*3,r3*3,-r3,r3,-r3*2,r3*2,r3*6,r3*6,-r3*6,r3*3,r3*3,9,-9,9,9,9,-9,-9,-r3*3,r3*3,r3*6,-r3*6,-9,-9,r3*3,r3*3,9,-9,r3*3,r3*3,9,9,-r3*2,-r3,-r3,-r3*6,-r3*6,-r3*6],
  [-r3,r3,-9,-9,9,9,9,r3*6,-r3*6,r3*3,-r3*3,r3*3,-r3*3,r3,-r3,r3*2,-r3*2,-r3*6,-r3*6,r3*6,-r3*3,-r3*3,9,-9,9,9,9,-9,-9,r3*3,-r3*3,-r3*6,r3*6,-9,-9,-r3*3,-r3*3,9,-9,-r3*3,-r3*3,9,9,r3*2,r3,r3,r3*6,r3*6,r3*6],
  [r3*2,-r3*2,0,0,0,0,0,(-3-r3)*3,(-3+r3)*3,(3+r3)*3,(3-r3)*3,(3+r3)*3,(3-r3)*3,-r3*2,r3*2,-r3*4,r3*4,(3+r3)*3,(-3+r3)*3,-r3*12,(3-r3)*3,(-3-r3)*3,0,0,0,0,0,0,0,(-3+r3)*3,(-3-r3)*3,(-3+r3)*3,(-3-r3)*3,0,0,(3-r3)*3,(-3-r3)*3,0,0,(3-r3)*3,(-3-r3)*3,0,0,-r3*4,-r3*2,-r3*2,r3*6,r3*6,r3*6],
  [-r3*2,r3*2,0,0,0,0,0,(-3+r3)*3,(-3-r3)*3,(3-r3)*3,(3+r3)*3,(3-r3)*3,(3+r3)*3,r3*2,-r3*2,r3*4,-r3*4,(3-r3)*3,(-3-r3)*3,r3*12,(3+r3)*3,(-3+r3)*3,0,0,0,0,0,0,0,(-3-r3)*3,(-3+r3)*3,(-3-r3)*3,(-3+r3)*3,0,0,(3+r3)*3,(-3+r3)*3,0,0,(3+r3)*3,(-3+r3)*3,0,0,r3*4,r3*2,r3*2,-r3*6,-r3*6,-r3*6],
  [(-3+r3)*3,(-3-r3)*3,0,0,0,0,0,-18*j2,-18*j,0,0,0,0,r3*6,-r3*6,(3+r3)*3,(3-r3)*3,18*j,-18*j2,0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,0,0,0,0,0,0,0,0,0,0,-r3*6,(-3-r3)*3,(3-r3)*3,0,0,0],
  [(3+r3)*3,(3-r3)*3,0,0,0,0,0,18*j,18*j2,0,0,0,0,r3*6,-r3*6,(-3+r3)*3,(-3-r3)*3,-18*j2,18*j,0,0,0,0,0,0,0,0,0,0,0,0,18,18,0,0,0,0,0,0,0,0,0,0,-r3*6,(3-r3)*3,(-3-r3)*3,0,0,0],
  [r3*6,-r3*6,0,0,0,0,0,0,0,0,0,0,0,-r3*6,r3*6,-r3*12,r3*12,0,0,r3*18,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-r3*12,-r3*6,-r3*6,0,0,0],
  [((3+r3)*3)//2,((3-r3)*3)//2,-9*j,-9*j2,9,9*j2,9*j,0,0,9*j,9*j2,9*j,9*j2,r3*3,-r3*3,(3-r3)*3,(3+r3)*3,0,0,0,9*j,-9*j2,9*j,-9*j2,9,9*j,9*j2,-9,-9,-9,-9,0,0,-9*j2,-9*j,9*j,-9*j2,9*j,-9*j2,9,-9,9,9,r3*6,((3-r3)*3)//2,((-3-r3)*3)//2,0,0,0],
  [((-3+r3)*3)//2,((-3-r3)*3)//2,9*j2,9*j,-9,-9*j,-9*j2,0,0,-9*j2,-9*j,-9*j2,-9*j,r3*3,-r3*3,(-3-r3)*3,(-3+r3)*3,0,0,0,-9*j2,9*j,-9*j2,9*j,-9,-9*j2,-9*j,9,9,9,9,0,0,9*j,9*j2,-9*j2,9*j,-9*j2,9*j,-9,9,-9,-9,r3*6,((-3-r3)*3)//2,((3-r3)*3)//2,0,0,0],
  [9*j,9*j2,-9*j,-9*j2,9,9*j2,9*j,0,0,9*j,9*j2,-9*j,-9*j2,9,9,0,0,0,0,0,9*j,-9*j2,9*j,-9*j2,-9,-9*j,-9*j2,-9,-9,-9,-9,0,0,9*j2,9*j,-9*j,9*j2,-9*j,9*j2,-9,9,-9,-9,0,-9*j2,9*j,0,0,0],
  [-9*j2,-9*j,9*j2,9*j,-9,-9*j,-9*j2,0,0,-9*j2,-9*j,9*j2,9*j,-9,-9,0,0,0,0,0,-9*j2,9*j,-9*j2,9*j,9,9*j2,9*j,9,9,9,9,0,0,-9*j,-9*j2,9*j2,-9*j,9*j2,-9*j,9,-9,9,9,0,9*j,-9*j2,0,0,0],
  [9,9,9,9,-9,-9,-9,0,0,9,9,-9,-9,9,9,0,0,0,0,0,9,-9,-9,9,9,9,9,9,9,-9,-9,0,0,-9,-9,-9,9,9,-9,-9,9,9,9,0,-9,9,0,0,0],
  [9*j2,9*j,9,9,-9,-9*j,-9*j2,0,0,9,9,-9,-9,9,9,0,0,0,0,0,9*j,-9*j2,-9*j,9*j2,9,9,9,9*j2,9*j,-9*j2,-9*j,0,0,-9*j,-9*j2,-9*j,9*j2,9*j,-9*j2,-9*j2,9*j,9*j2,9*j,0,-9*j,9*j2,0,0,0],
  [9*j,9*j2,9,9,-9,-9*j2,-9*j,0,0,9,9,-9,-9,9,9,0,0,0,0,0,9*j2,-9*j,-9*j2,9*j,9,9,9,9*j,9*j2,-9*j,-9*j2,0,0,-9*j2,-9*j,-9*j2,9*j,9*j2,-9*j,-9*j,9*j2,9*j,9*j2,0,-9*j2,9*j,0,0,0],
  [-9*j,-9*j2,9*j2,9*j,-9,-9*j2,-9*j,0,0,-9*j2,-9*j,9*j2,9*j,-9,-9,0,0,0,0,0,-9,9,-9,9,9,9*j2,9*j,9*j2,9*j,9*j2,9*j,0,0,-9*j2,-9*j,9,-9,9,-9,9*j2,-9*j,9*j2,9*j,0,9*j2,-9*j,0,0,0],
  [-9*j2,-9*j,9*j,9*j2,-9,-9*j,-9*j2,0,0,-9*j,-9*j2,9*j,9*j2,-9,-9,0,0,0,0,0,-9,9,-9,9,9,9*j,9*j2,9*j,9*j2,9*j,9*j2,0,0,-9*j,-9*j2,9,-9,9,-9,9*j,-9*j2,9*j,9*j2,0,9*j,-9*j2,0,0,0],
  [((-3-r3)*3)//2,((-3+r3)*3)//2,9*j2,9*j,-9,-9*j2,-9*j,0,0,-9*j2,-9*j,-9*j2,-9*j,-r3*3,r3*3,(-3+r3)*3,(-3-r3)*3,0,0,0,-9,9,-9,9,-9,-9*j2,-9*j,9*j2,9*j,9*j2,9*j,0,0,9*j2,9*j,-9,9,-9,9,-9*j2,9*j,-9*j2,-9*j,-r3*6,((-3+r3)*3)//2,((3+r3)*3)//2,0,0,0],
  [((-3+r3)*3)//2,((-3-r3)*3)//2,9*j,9*j2,-9,-9*j,-9*j2,0,0,-9*j,-9*j2,-9*j,-9*j2,r3*3,-r3*3,(-3-r3)*3,(-3+r3)*3,0,0,0,-9,9,-9,9,-9,-9*j,-9*j2,9*j,9*j2,9*j,9*j2,0,0,9*j,9*j2,-9,9,-9,9,-9*j,9*j2,-9*j,-9*j2,r3*6,((-3-r3)*3)//2,((3-r3)*3)//2,0,0,0],
  [(3+r3)*3,(3-r3)*3,0,0,0,0,0,18*j2,18*j,0,0,0,0,r3*6,-r3*6,(-3+r3)*3,(-3-r3)*3,-18,18,0,0,0,0,0,0,0,0,0,0,0,0,18*j2,18*j,0,0,0,0,0,0,0,0,0,0,-r3*6,(3-r3)*3,(-3-r3)*3,0,0,0],
  [(3-r3)*3,(3+r3)*3,0,0,0,0,0,18*j,18*j2,0,0,0,0,-r3*6,r3*6,(-3-r3)*3,(-3+r3)*3,-18,18,0,0,0,0,0,0,0,0,0,0,0,0,18*j,18*j2,0,0,0,0,0,0,0,0,0,0,r3*6,(3+r3)*3,(-3+r3)*3,0,0,0],
  [-9,-9,-9*j,-9*j2,9,9,9,0,0,-9*j,-9*j2,9*j,9*j2,-9,-9,0,0,0,0,0,-9*j2,9*j,9*j2,-9*j,-9,-9*j,-9*j2,-9*j2,-9*j,9*j2,9*j,0,0,9,9,9*j2,-9*j,-9*j2,9*j,9*j2,-9*j,-9*j2,-9*j,0,9,-9,0,0,0],
  [-9,-9,-9*j2,-9*j,9,9,9,0,0,-9*j2,-9*j,9*j2,9*j,-9,-9,0,0,0,0,0,-9*j,9*j2,9*j,-9*j2,-9,-9*j2,-9*j,-9*j,-9*j2,9*j,9*j2,0,0,9,9,9*j,-9*j2,-9*j,9*j2,9*j,-9*j2,-9*j,-9*j2,0,9,-9,0,0,0],
  [((3+r3)*3)//2,((3-r3)*3)//2,9*j,9*j2,-9,-9*j2,-9*j,0,0,9*j,9*j2,9*j,9*j2,r3*3,-r3*3,(3-r3)*3,(3+r3)*3,0,0,0,9*j,-9*j2,-9*j,9*j2,-9,-9*j,-9*j2,9,9,-9,-9,0,0,9*j2,9*j,9*j,-9*j2,-9*j,9*j2,9,-9,-9,-9,r3*6,((3-r3)*3)//2,((-3-r3)*3)//2,0,0,0],
  [((-3+r3)*3)//2,((-3-r3)*3)//2,-9*j2,-9*j,9,9*j,9*j2,0,0,-9*j2,-9*j,-9*j2,-9*j,r3*3,-r3*3,(-3-r3)*3,(-3+r3)*3,0,0,0,-9*j2,9*j,9*j2,-9*j,9,9*j2,9*j,-9,-9,9,9,0,0,-9*j,-9*j2,-9*j2,9*j,9*j2,-9*j,-9,9,9,9,r3*6,((-3-r3)*3)//2,((3-r3)*3)//2,0,0,0],
  [9*j,9*j2,9*j,9*j2,-9,-9*j2,-9*j,0,0,9*j,9*j2,-9*j,-9*j2,9,9,0,0,0,0,0,9*j,-9*j2,-9*j,9*j2,9,9*j,9*j2,9,9,-9,-9,0,0,-9*j2,-9*j,-9*j,9*j2,9*j,-9*j2,-9,9,9,9,0,-9*j2,9*j,0,0,0],
  
  [-9*j2,-9*j,-9*j2,-9*j,9,9*j,9*j2,0,0,-9*j2,-9*j,9*j2,9*j,-9,-9,0,0,0,0,0,-9*j2,9*j,9*j2,-9*j,-9,-9*j2,-9*j,-9,-9,9,9,0,0,9*j,9*j2,9*j2,-9*j,-9*j2,9*j,9,-9,-9,-9,0,9*j,-9*j2,0,0,0],
  [((3+r3)*3)//2,((3-r3)*3)//2,9*j2,9*j,-9,-9*j2,-9*j,0,0,9*j2,9*j,9*j2,9*j,r3*3,-r3*3,(3-r3)*3,(3+r3)*3,0,0,0,9,-9,-9,9,-9,-9*j2,-9*j,9*j2,9*j,-9*j2,-9*j,0,0,9*j2,9*j,9,-9,-9,9,9*j2,-9*j,-9*j2,-9*j,r3*6,((3-r3)*3)//2,((-3-r3)*3)//2,0,0,0],
  [((-3+r3)*3)//2,((-3-r3)*3)//2,-9*j,-9*j2,9,9*j,9*j2,0,0,-9*j,-9*j2,-9*j,-9*j2,r3*3,-r3*3,(-3-r3)*3,(-3+r3)*3,0,0,0,-9,9,9,-9,9,9*j,9*j2,-9*j,-9*j2,9*j,9*j2,0,0,-9*j,-9*j2,-9,9,9,-9,-9*j,9*j2,9*j,9*j2,r3*6,((-3-r3)*3)//2,((3-r3)*3)//2,0,0,0],
  [9*j,9*j2,9*j2,9*j,-9,-9*j2,-9*j,0,0,9*j2,9*j,-9*j2,-9*j,9,9,0,0,0,0,0,9,-9,-9,9,9,9*j2,9*j,9*j2,9*j,-9*j2,-9*j,0,0,-9*j2,-9*j,-9,9,9,-9,-9*j2,9*j,9*j2,9*j,0,-9*j2,9*j,0,0,0],
  [9*j2,9*j,9*j,9*j2,-9,-9*j,-9*j2,0,0,9*j,9*j2,-9*j,-9*j2,9,9,0,0,0,0,0,9,-9,-9,9,9,9*j,9*j2,9*j,9*j2,-9*j,-9*j2,0,0,-9*j,-9*j2,-9,9,9,-9,-9*j,9*j2,9*j,9*j2,0,-9*j,9*j2,0,0,0],
  [r3*2,-r3*2,0,0,0,0,0,r3*6,-r3*6,-r3*6,r3*6,-r3*6,r3*6,-r3*2,r3*2,-r3*4,r3*4,-r3*6,-r3*6,-r3*12,r3*6,r3*6,0,0,0,0,0,0,0,-r3*6,r3*6,-r3*6,r3*6,0,0,r3*6,r3*6,0,0,r3*6,r3*6,0,0,-r3*4,-r3*2,-r3*2,r3*6,r3*6,r3*6],
  [r3,-r3,9*j,9*j2,-9,-9,-9,(3+r3)*3,(3-r3)*3,((3+r3)*3)//2,((3-r3)*3)//2,((3+r3)*3)//2,((3-r3)*3)//2,-r3,r3,-r3*2,r3*2,(-3-r3)*3,(3-r3)*3,-r3*6,((3-r3)*3)//2,((-3-r3)*3)//2,-9*j2,9*j,-9,-9*j,-9*j2,9*j2,9*j,((-3+r3)*3)//2,((-3-r3)*3)//2,(3-r3)*3,(3+r3)*3,9,9,((3-r3)*3)//2,((-3-r3)*3)//2,-9*j2,9*j,((3-r3)*3)//2,((-3-r3)*3)//2,-9*j2,-9*j,-r3*2,-r3,-r3,-r3*6,-r3*6,-r3*6],
  [r3,-r3,-9*j2,-9*j,9,9,9,(-3+r3)*3,(-3-r3)*3,((-3+r3)*3)//2,((-3-r3)*3)//2,((-3+r3)*3)//2,((-3-r3)*3)//2,-r3,r3,-r3*2,r3*2,(3-r3)*3,(-3-r3)*3,-r3*6,((-3-r3)*3)//2,((3-r3)*3)//2,9*j,-9*j2,9,9*j2,9*j,-9*j,-9*j2,((3+r3)*3)//2,((3-r3)*3)//2,(-3-r3)*3,(-3+r3)*3,-9,-9,((-3-r3)*3)//2,((3-r3)*3)//2,9*j,-9*j2,((-3-r3)*3)//2,((3-r3)*3)//2,9*j,9*j2,-r3*2,-r3,-r3,-r3*6,-r3*6,-r3*6],
  [r3*6,-r3*6,0,0,0,0,0,0,0,0,0,0,0,-r3*6,r3*6,r3*6,-r3*6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,r3*6,-r3*6,-r3*6,18*E(9,7)-18*E(9,2),-18*E(9,5)+18*E(9,4),((-18*E(9,7)+18*E(9,5))-18*E(9,4))+18*E(9,2)],
  [r3*6,-r3*6,0,0,0,0,0,0,0,0,0,0,0,-r3*6,r3*6,r3*6,-r3*6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,r3*6,-r3*6,-r3*6,-18*E(9,5)+18*E(9,4),((-18*E(9,7)+18*E(9,5))-18*E(9,4))+18*E(9,2),18*E(9,7)-18*E(9,2)],
  [r3*6,-r3*6,0,0,0,0,0,0,0,0,0,0,0,-r3*6,r3*6,r3*6,-r3*6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,r3*6,-r3*6,-r3*6,((-18*E(9,7)+18*E(9,5))-18*E(9,4))+18*E(9,2),18*E(9,7)-18*E(9,2),-18*E(9,5)+18*E(9,4)]])//54,
  :eigenvalues=>[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,j2,j2,j2,j2,j2,j2,j2,-1,-1,-1,j,j,j,j,j,j,-1,-1,j2,j2,-j2,-j2,j,j,-j,-j,1,1,1,E(9,8),E(9,5),E(9,2)],
  :explanation=>"mysteryG26",:special=>1,:cospecial=>2)
end

chevieset(:G26, :UnipotentCharacters,Dict{Symbol, Any}(:harishChandra => [
  Dict{Symbol,Any}(:relativeType=>
    TypeIrred(;series=:ST,indices=1:3,rank=3,ST=26),
    :levi => [], :parameterExponents => [1, 1, 1],
    :charNumbers => 1:48, :eigenvalue => 1, :cuspidalName => ""), 
  Dict{Symbol, Any}(:relativeType =>
     TypeIrred(;series=:ST,indices=[1,3,13],rank=2,p=6,q=2),
     :levi => [2], :parameterExponents => [[0, 2, 2], 3, 1],
     :charNumbers=>[102,68,71,66,53,70,60,67,54,103,69,72,99,59,98,65,50,49],
     :eigenvalue=>E(3)^2,
     :cuspidalName=>ImprimitiveCuspidalName([Int[],[0,1],[0,1]])), 
  Dict{Symbol, Any}(:relativeType =>
     TypeIrred(;series=:ST,indices=[1],rank=1,p=6,q=1), :levi => 2:3,
     :parameterExponents => [[3, 4, 3, 0, 3, 4]],
     :charNumbers=>[73,61,74,104,75,62],:eigenvalue=>-1,:cuspidalName=>"G_4"), 
  Dict{Symbol, Any}(:relativeType =>
     TypeIrred(;series=:ST,indices=[3],rank=1,p=6,q=1),:levi => 1:2,
     :parameterExponents => [[4, 3, 1, 1, 0, 1]],
     :charNumbers => [51, 55, 76, 81, 100, 78], :eigenvalue => E(3),
     :cuspidalName => ImprimitiveCuspidalName([[0],Int[], [0, 1, 2]])), 
  Dict{Symbol, Any}(:relativeType =>
     TypeIrred(;series=:ST,indices=[3],rank=1,p=6,q=1), :levi => 1:2,
     :parameterExponents => [[4, 1, 0, 1, 1, 3]],
     :charNumbers => [52, 79, 101, 80, 77, 56], :eigenvalue => E(3),
     :cuspidalName => ImprimitiveCuspidalName([[0], [0, 1, 2],Int[]])), 
  mkcuspidal("G_{26}",92,1),
  mkcuspidal("G_{26}",93,1;no=2),
  mkcuspidal("G_{26}",94,1;no=3),
  mkcuspidal("G_{26}",82,-1),
  mkcuspidal("G_{26}",83,-1;no=2),
  mkcuspidal("G_{26}",88,E(3)),
  mkcuspidal("G_{26}",89,E(3);no=2),
  mkcuspidal("G_{26}",64,E(3,2)),
  mkcuspidal("G_{26}",84,E(3,2);no=2),
  mkcuspidal("G_{26}",85,E(3,2);no=3),
  mkcuspidal("G_{26}",90,-E(3)),
  mkcuspidal("G_{26}",91,-E(3);no=2),
  mkcuspidal("G_{26}",63,-E(3,2)),
  mkcuspidal("G_{26}",86,-E(3,2);no=2),
  mkcuspidal("G_{26}",87,-E(3,2);no=3),
  mkcuspidal("G_{26}",57,E(4);qeig=1//2),
  mkcuspidal("G_{26}",58,-E(4);qeig=1//2),
  mkcuspidal("G_{26}",95,E(9,8)),
  mkcuspidal("G_{26}",96,E(9,5)),
  mkcuspidal("G_{26}",97,E(9,2))], 
  :families => [Family(:C1,[1]), 
  Family(Family(:QZ)(3,[Perm(),[E(3)]]),[24,18,2,52,50,10,49,12,51],
        signs=[-1, -1, -1, 1, 1, 1, -1, 1, 1],ennola=4), 
  Family(Family(:QZ)(3,[Perm(),[E(3)]]),[31,37,13,55,53,23,54,17,56],
        signs=[-1, -1, -1, 1, 1, 1, -1, 1, 1],ennola=-9), 
  Family(Family(:TQZ)(2, -1, (1, -1)), [40, 39, 58, 57],cospecial=2,ennola=4), 
  Family(Family(:C2)*Family(:X)(3),[33,27,59,22,16,60,48,46,64,61,62,63],
         signs=[1,1,-1,-1,-1,-1,1,1,1,-1,-1,1],ennola=12,cospecial=2), 
  Family(Family(:X)(3),[32,38,65],signs=[1,1,-1],ennola=2), 
  Family(:X49, [43, 42, 28, 34, 8, 41, 44, 29, 35, 15, 21, 45, 47, 6, 5, 11, 9, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1], :ennola => -2)), 
  Family(Family(:QZ)(3, [Perm(), [E(3)]]),[30,36,14,101,99,26,98,20,100],
         signs = [-1, -1, -1, 1, 1, 1, -1, 1, 1], ennola = 9), 
  Family(Family(:X)(3),[25,19,102],signs=[1,1,-1],ennola=-3), 
  Family(:X5,[4,7,104,103,3],signs=[1, 1, -1, -1, 1],ennola=5)], 
  :a => [0, 1, 21, 21, 6, 6, 21, 6, 6, 1, 6, 1, 2, 11, 6, 4, 2, 1, 16, 11, 6, 4, 2, 1, 16, 11, 4, 6, 6, 11, 2, 5, 4, 6, 6, 11, 2, 5, 3, 3, 6, 6, 6, 6, 6, 4, 6, 4, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 11, 11, 11, 11, 16, 21, 21], 
  :A => [0, 17, 33, 33, 30, 30, 33, 30, 30, 17, 30, 17, 22, 31, 30, 26, 22, 17, 32, 31, 30, 26, 22, 17, 32, 31, 26, 30, 30, 31, 22, 25, 26, 30, 30, 31, 22, 25, 24, 24, 30, 30, 30, 30, 30, 26, 30, 26, 17, 17, 17, 17, 22, 22, 22, 22, 24, 24, 26, 26, 26, 26, 26, 26, 25, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 32, 33, 33]))

chevieset(:G26, :Invariants,[
  (x1,x2,x3)->-10*x1^3*x2^3-10*x1^3*x3^3-10*x2^3*x3^3+x1^6+x2^6+x3^6,
  (x1,x2,x3)->2*x1^3*x2^3*x3^6+2*x1^3*x2^6*x3^3+x1^3*x2^9+x1^3*x3^9+x2^3*x3^9+
  2*x1^6*x2^3*x3^3-4*x1^6*x2^6-4*x1^6*x3^6-4*x2^6*x3^6+x1^9*x2^3+x1^9*x3^3+
  x2^9*x3^3,
  (x1,x2,x3)->-2*x1^3*x2^3*x3^12+2*x1^3*x2^6*x3^9+2*x1^3*x2^9*x3^6-
  2*x1^3*x2^12*x3^3+2*x1^6*x2^3*x3^9-6*x1^6*x2^6*x3^6+2*x1^6*x2^9*x3^3+
  x1^6*x2^12+x1^6*x3^12+x2^6*x3^12+2*x1^9*x2^3*x3^6+2*x1^9*x2^6*x3^3-
  2*x1^9*x2^9-2*x1^9*x3^9-2*x2^9*x3^9-2*x1^12*x2^3*x3^3+x1^12*x2^6+x1^12*x3^6+
  x2^12*x3^6])

chevieset(:G26, :Discriminant,()->
  (t1,t2,t3)->36*t1*t2*t3^2-t1^2*t2^2*t3+108*t3^3-32*t2^3*t3+t1^3*t3^2)
