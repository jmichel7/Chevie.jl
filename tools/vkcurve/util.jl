# alternative to ComplexOps.Norm(x)
BigNorm(x)=abs(real(x))+abs(imag(x))

# another alternative
SmallNorm(x)=maximum(abs(real(x)),abs(imag(x)))

# For a non-zero rational x, returns k such that 10^k<|x|<=10^(k+1) 
function DecimalLog(x)
  if iszero(x) error("trying to take decimal log of 0") end
  d = div(denominator(x), numerator(x))
  if iszero(d)
    d=div(numerator(x), denominator(x))
    return length(string(abs(d)))-1
  else return -length(string(abs(d)))
  end
end

# a complex rational approximation of a (a Cyc or ??)
function ComplexRational(a)
  if a isa Cyc
      a = Complex(a)
      if !(IsRat(a[:r]))
          a[:r] = evalf(a[:r])
      end
      if !(IsRat(a[:i]))
          a[:i] = evalf(a[:i])
      end
  end
  return Complex(Rational(a[:r]), Rational(a[:i]))
end

"""
`Dispersal(v)`
`v`  is a list of complex numbers  representing points in the real plane.
The  result is a pair  whose first element is  the minimum distance between
two  elements of `v`, and the second is a pair of indices `[i,j]` such that
`v[i]`, `v[j]` achieves this minimum distance.

julia> Dispersal([1+im,0,1])
(1, [1,3])
"""
function Dispersal(v)
  l=combinations(eachindex(v),2)
  p=findmin(map(x->BigNorm(v[x[1]]-v[x[2]]),l))
  (p[1],l[p[2]])
end

# distance of z to segment [a,b]
function DistSeg(z,a,b)
  b-=a
  z-=a
  r=abs(b)
  z=r//b*z
  real(z)<0 ? abs(z) : real(z)>r ? abs(z-b) : imag(z)>0 ? imag(z) : -imag(z)
end

InterpolatedPolynomial=function(x,y)
  t=copy(y)*1//1 # make sure coeffs are in a field
  a=map(1:length(x))do i
    for k in i-1:-1:1
      t[k]=(t[k+1]-t[k])/(x[i]-x[k])
      if isnothing(t[k]) error("cannot interpolate polynomial") end
    end
    t[1]
  end
  p=a[length(x)]
  for i in length(x)-1:-1:1
    p=p*(Pol()-x[i])+a[i]
  end
  p
end

"""
`discy2(p::Mvp)`

The input should be an 'Mvp' in 'x' and 'y', with rational coefficients.
The  function  returns the  discriminant  of  `p`  with respect  to  'x'
(an  `Mvp` in  `y`);  it uses  interpolation to  reduce  the problem  to
discriminants of univariate polynomials,  and works reasonably fast (not
hundreds of times slower than MAPLE...).

|    gap> discy2(x+y^2+x^3+y^3);      
    4+27y^4+54y^5+27y^6|
"""
function discy2(p)
  n=2*(1+degree(p,:x))*(1+degree(p,:y))
  v=map(1:n)do i
    q=scal.(Pol(p(y=i),:x).c)
    if isempty(q) return 0
    elseif length(q)==1 return q[1]
    else return GLinearAlgebra.det(resultant(q,map(*,q[2:end],1:length(q)-1)))
    end
  end
  InterpolatedPolynomial(1:n,v)
end

Chars.discriminant(p::Pol)=GLinearAlgebra.det(resultant(p[0:end],
                                                 derivative(p)[0:end]))

discy(p)=Pol(discriminant(Pol(p,:x)))
  
"""
`resultant(v,w)`

`v` and  `w` are vectors  representing coefficients of  two polynomials.
The function returns  Sylvester matrix for these  two polynomials (whose
determinant  is  the resultant  of  the  two  polynomials). It  is  used
internally by Discy.

```julia-repl
julia> @Mvp x,y

julia> p=x+y^2+x^3+y^3
Mvp{Int64}: x³+x+y³+y²

julia> c=Pol(p,:x).c
4-element Vector{Mvp{Int64, Int64}}:
 y³+y²
 1
 0
 1

julia> d=Pol(derivative(p,:x),:x).c
3-element Vector{Mvp{Int64, Int64}}:
 1
 0
 3

julia> resultant(c,d)
5×5 Matrix{Mvp{Int64, Int64}}:
 1  0  1  y³+y²  0
 0  1  0  1      y³+y²
 3  0  1  0      0
 0  3  0  1      0
 0  0  3  0      1

julia> GLinearAlgebra.det(m)
Mvp{Int64}: 27y⁶+54y⁵+27y⁴+4
```
"""
function resultant(v,w)
  v=reverse(v)
  w=reverse(w)
  l=max(0,length(v)+length(w)-2)
  m=fill(zero(v[1]),l,l)
  for i in 1:length(w)-1 m[i,i:i+length(v)-1]=v end
  for i in 1:length(v)-1 m[i+length(w)-1,i:i+length(w)-1]=w end
  m
end
