"""
Cédric's problem.

For  a Weyl group `W` the algebra `Z(ℂ W)` has basis `e_C:=∑{w∈ C} w` where
`C` runs over conjugacyclasses of `W`, and is filtered by
`i(e_C)=\\codim(\\ker(w-\\id))`  for `w∈  C`. Let  `A_W` be  the subalgebra
`A_W=∑_F ℂ e_F` where `e_F=∑{φ∈ F}e_φ` for a family `F`.

For  a cuspidal pair `(L,λ)` there is a map `F↦ F'` from families of `W` to
those  of `W':=W(L,λ)`,  whence a  map `π:e_F↦  e_F':A_W→A_W'`. Is this map
filtered? That is `π(A_W∩ ≤i)⊂ A_{W'}∩≤i`?
"""

@GapObj struct CenterAlgebra
  group::Group
  classes::Matrix
end

function CenterAlgebra(W)
  ct=CharTable(W).irr
  cl=classinfo(W).classes
  dim=Int.(ct[:,1])
  classes=[div(cl[i]*ct[j,i],dim[j]) for i in eachindex(cl), j in eachindex(cl)]
  CenterAlgebra(W,classes,Dict{Symbol,Any}())
end

Base.show(io::IO,r::CenterAlgebra)=print(io,"CenterAlgebra(",r.group,")")

function filtration(A::CenterAlgebra)
  get!(A,:filtration)do
    re=refleigen(A.group)
    map(i->map(j->element(A,j),
      filter(j->count(!isone,re[j])==i,axes(re,1))),0:length(re[1]))
  end
end

"""
elements  of center of group algebra  are represented internally as vectors
in basis e_φ
"""
struct CenterElement{T}
  A::CenterAlgebra
  v::AbstractVector{T}
end

Base.:*(x::CenterElement,y::CenterElement)=CenterElement(x.A,x.v .* y.v)

Base.:+(x::CenterElement,y::CenterElement)=CenterElement(x.A,x.v .+ y.v)

Base.one(x::CenterElement)=element(x.A,1)

Base.:^(h::CenterElement,n)=Base.power_by_squaring(h,n)

"elements are printed in basis of class sums unless :classes=false"
function Base.show(io::IO,x::CenterElement)
  if get(io,:classes,true)
    v=solutionmat(x.A.classes,x.v)
    letter="C"
    n=classinfo(x.A.group)[:classnames]
  else v=x.v
    letter="e"
    n=charinfo(x.A.group)[:charnames]
  end
  u=findall(!iszero,v)
  if isempty(u) print(io,"0");return  end
  v=prod(map(u)do i
    num=numerator(v[i])
    s=format_coefficient(xrepr(num,typeinfo=typeof(num)))*letter*"("*
                       fromTeX(io,n[i])*")"
    d=denominator(v[i])
    if !isone(d) s*="/$d" end
    if s[1]!='-' return "+"*s else return s end
    end
  )
  if v[1]=='+' v=v[2:end] end
  print(io,v)
end

element(A::CenterAlgebra,i::Integer)=CenterElement(A,A.classes[i,:])

# intersections of rowspace of m with f<=i
function filtration(A::CenterAlgebra,m::AbstractMatrix)
  map(i->intersect_rowspace(m,toM(y.v for j in 1:i for y in filtration(A)[j])),
      eachindex(filtration(A)))
end

# test all d-series
function test(W,ss=Series(W;proper=true);typ=Int64)
  res=Series[]
  A=CenterAlgebra(W)
  uc=UnipotentCharacters(W)
  n=uc.harishChandra[1][:charNumbers]
  ff=sort(uc.families,by=charnumbers)
  A.families=map(f->CenterElement(A,in.(n,Ref(f.charNumbers))),ff)
  ffamA=filtration(A,typ.(toM(map(x->x.v,A.families))))
  for s in ss
    xprintln("*****",s)
    A1=CenterAlgebra(relative_group(s))
    n=dSeries.charnumbers(s)
    A1.families=map(f->CenterElement(A1,in.(n,Ref(f.charNumbers))),ff)
    pi=toM(map(x->x.v,A1.families))
    A1.families=filter(x->!iszero(x.v),A1.families)
    fpi=v->CenterElement(A1,vec(solutionmat(toM(map(x->x.v,A.families)),v)'*pi))
    image=map(x->fpi.(eachrow(x)),ffamA)
    for i in eachindex(filtration(A1))
     l=toM(map(x->x.v,image[i]))
      if size(filtration(A1,Rational{typ}.(l))[i],1)<GenLinearAlgebra.rank(l)
        println("CONTRE-EXEMPLE")
        push!(res,s)
      end
    end
  end
  if !isempty(res) println("Les mauvaises series sont :")
  else println("Good !")
  end
  res
end
