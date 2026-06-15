using Oscar: Oscar
using Chevie

fromOscar(x::Oscar.Generic.MatSpaceElem)=map(fromOscar,Matrix(x))

function fromOscar(x::Oscar.Generic.MPoly)
  vv=Symbol.(string.(Oscar.gens(parent(x))))
  c=fromOscar.(collect(Oscar.coefficients(x)))
  l=map(collect(Oscar.monomials(x)))do m
    prod(map(1:length(vv))do i
       Mvp(vv[i])^Oscar.degree(m,i)
     end)
  end
  sum(map(*,c,l))
end

function fromOscar(x::Oscar.ZZMPolyRingElem)
  vv=Mvp.(Symbol.(string.(Oscar.gens(parent(x)))))
  sum(Oscar.coefficients_and_exponents(x);init=0*vv[1])do (c,e)
    Integer(c)*prod(map((v,exp)->v^exp,vv,e))
  end
end

fromOscar(x::Oscar.Generic.FracFieldElem)=
    fromOscar(numerator(x))//fromOscar(denominator(x))

function fromOscar(x::Oscar.AbelianClosure.QQAbFieldElem)
  if isinteger(x) Int(Oscar.ZZ(x))
  else error(x," is not an integer")
  end
end
