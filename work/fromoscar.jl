using Oscar: Oscar
using Chevie

function fromOscar(x)
  if x isa Oscar.Generic.MatSpaceElem
    map(fromOscar,Matrix(x))
  elseif x isa Oscar.Generic.FracFieldElem
    fromOscar(numerator(x))//fromOscar(denominator(x))
  elseif x isa Oscar.Generic.MPoly
    vv=Symbol.(string.(Oscar.gens(parent(x))))
    c=fromOscar.(collect(Oscar.coefficients(x)))
    l=map(collect(Oscar.monomials(x)))do m
      prod(map(1:length(vv))do i
         Mvp(vv[i])^Oscar.degree(m,i)
       end)
    end
    sum(map(*,c,l))
  elseif x isa Oscar.AbelianClosure.QQAbFieldElem
    if isinteger(x) Int(Oscar.ZZ(x))
    else error(x," is not an integer")
    end
  end
end
