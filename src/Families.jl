"""
Families of unipotent characters

The blocks of the (rectangular) matrix ``⟨Rᵪ,ρ⟩_{𝐆 ^F}`` when `χ` runs over
`Irr(W)`  and  `ρ`  runs  over  the  unipotent  characters,  are called the
*Lusztig  families*. When  `𝐆 `  is split  and `W`  is a Coxeter group they
correspond  on the `Irr(W)` side to two-sided Kazhdan-Lusztig cells --- for
split  Spetses they  correspond to  Rouquier blocks  of the  Spetsial Hecke
algebra.  The matrix of scalar products  ``⟨Rᵪ,ρ⟩_{𝐆 ^F}`` can be completed
to   a  square  matrix  ``⟨A_{ρ'},ρ⟩_{𝐆  ^F}``  where  ``A_{ρ'}``  are  the
*characteristic  functions of character  sheaves* on ``𝐆  ^F``; this square
matrix is called the *Fourier matrix* of the family.

The  'UnipotentCharacters' record in Chevie contains a field '.families', a
list of family records containing information on each family, including the
Fourier matrix. Here is an example.

```julia-repl
julia> W=coxgroup(:G,2)
G₂

julia> uc=UnipotentCharacters(W);

julia> uc.families
3-element Vector{Family}:
 Family(D(𝔖 ₃),[5, 6, 4, 3, 8, 7, 9, 10])
 Family(C₁,[1])
 Family(C₁,[2])

julia> uc.families[1]
Family(D(𝔖 ₃),[5, 6, 4, 3, 8, 7, 9, 10])
Drinfeld double of 𝔖 ₃, Lusztig′s version
   label│eigen
────────┼─────────────────────────────────────────────────────
(1,1)   │    1 1//6  1//2  1//3  1//3  1//6  1//2  1//3  1//3
(g₂,1)  │    1 1//2  1//2     0     0 -1//2 -1//2     0     0
(g₃,1)  │    1 1//3     0  2//3 -1//3  1//3     0 -1//3 -1//3
(1,ρ)   │    1 1//3     0 -1//3  2//3  1//3     0 -1//3 -1//3
(1,ε)   │    1 1//6 -1//2  1//3  1//3  1//6 -1//2  1//3  1//3
(g₂,ε)  │   -1 1//2 -1//2     0     0 -1//2  1//2     0     0
(g₃,ζ₃) │   ζ₃ 1//3     0 -1//3 -1//3  1//3     0  2//3 -1//3
(g₃,ζ₃²)│  ζ₃² 1//3     0 -1//3 -1//3  1//3     0 -1//3  2//3

julia> charnames(uc)[uc.families[1].charNumbers]
8-element Vector{String}:
 "phi2,1"
 "phi2,2"
 "phi1,3''"
 "phi1,3'"
 "G2[1]"
 "G2[-1]"
 "G2[E3]"
 "G2[E3^2]"
```

The  Fourier matrix is obtained  by 'fourier(f)'; the field 'f.charNumbers'
holds  the indices of the unipotent characters  which are in the family. We
obtain  the list of eigenvalues of Frobenius for these unipotent characters
by  'Eigenvalues(f)'. The Fourier matrix  and vector of eigenvalues satisfy
the  properties of  *fusion data*,  see below.  The field 'f.charLabels' is
what  is displayed  in the  column 'labels'  when displaying the family. It
contains  labels naturally attached to lines  of the Fourier matrix. In the
case   of  reductive  groups,   the  family  is   always  attached  to  the
"drinfeld_double"  of a small finite group  and the '.charLabels' come from
this construction.
"""
module Families

export family_imprimitive, Family, drinfeld_double, fourier,
 FamiliesClassical, SubFamilyij, ndrinfeld_double, fusion_algebra,
 involution, duality

using ..Gapjm

@GapObj struct Family end

Base.setindex!(f::Family,x,s::Symbol)=setproperty!(f,s,x) # used in cmplximp.jl

Family(f::Family)=f

function getf(s::String)
  f=chevieget(:families,Symbol(s))
  (f isa Dict) ? Family(deepcopy(f)) : Family(deepcopy(f.prop))
end

"""
`Family(f [, charNumbers [, opt]])`

This function creates a new family in two possible ways.

In  the first case `f` is a string which denotes a family known to  Chevie.
Examples are "S3",   "S4",   "S5"   which denote the family obtained as the
Drinfeld  double of the symmetric group  on 3,4,5 elements, or "C2"   which
denotes the Drinfeld double of the cyclic group of order 2.

In the second case `f` is already a struct Family.

The other (optional) arguments add information to the family defined by the
first argument. If given, the second argument becomes `f.charNumbers`. If
given,  the third argument  `opt` is a  `Dict` whose keys  are added to the
resulting family.

If `opt` has a key `signs`, this should be a list of '1' and '-1', and then
the  Fourier matrix  is conjugated  by the  diagonal matrix of those signs.
This  is used  in Spetses  to adjust  the matrix  to the choice of signs of
unipotent degrees.

```julia-repl
julia> Family("C2")
Family(C₂,4)
DrinfeldDouble(Z/2)
 label│eigen
──────┼─────────────────────────────
(1,1) │    1 1//2  1//2  1//2  1//2
(g₂,1)│    1 1//2  1//2 -1//2 -1//2
(1,ε) │    1 1//2 -1//2  1//2 -1//2
(g₂,ε)│   -1 1//2 -1//2 -1//2  1//2

julia> Family("C2",4:7,Dict(:signs=>[1,-1,1,-1]))
Family(C₂,4:7)
DrinfeldDouble(Z/2)
 label│eigen signs
──────┼───────────────────────────────────
(1,1) │    1     1  1//2 -1//2 1//2 -1//2
(g₂,1)│    1    -1 -1//2  1//2 1//2 -1//2
(1,ε) │    1     1  1//2  1//2 1//2  1//2
(g₂,ε)│   -1    -1 -1//2 -1//2 1//2  1//2
```
"""
function Family(s::String,v::AbstractVector,d::Dict=Dict{Symbol,Any}())
  f=getf(s)
  f.charNumbers=v
  merge!(f,d)
end

function Family(s::String,d::Dict=Dict{Symbol,Any}())
  f=getf(s)
  merge!(f,d)
end

function Family(f::Family,v::AbstractVector,d::Dict=Dict{Symbol,Any}())
  f.charNumbers=v
  merge!(f,d)
end

function Family(f::Dict{Symbol,Any},v::AbstractVector,d::Dict=Dict{Symbol,Any}())
  f=Family(f)
  f.charNumbers=v
  merge!(f,d)
end

special(f::Family)::Int=get!(()->1,f,:special)

PermRoot.eigen(f::Family)=f.eigenvalues

"`fourier(f::Family`: returns the Fourier matrix for the family `f`."
function fourier(f::Family)
  m=f.fourierMat
  m isa Vector ? improve_type(toM(m)) : m
end

"`length(f::Family)`: how many characters are in the family."
Base.length(f::Family)=length(eigen(f))
Base.convert(::Type{Dict{Symbol,Any}},f::Family)=f.prop

function Base.merge!(f::Family,d::Dict)
  merge!(f.prop,d)
  if f.fourierMat isa Vector f.fourierMat=improve_type(toM(f.fourierMat)) end
  if !haskey(f,:charLabels) f.charLabels=map(string,1:length(f)) end
  if !haskey(f,:special) f.special=1 end
  if haskey(d,:signs)
    signs=d[:signs]
    m=f.fourierMat
    for i in axes(m,1), j in axes(m,2) m[i,j]*=signs[i]*signs[j] end
    f.fourierMat=m
    if haskey(f,:perm) && -1 in f.fourierMat^2 delete!(f.prop,:perm) end
  end
  f
end

"""
`<f>*<g>`:  returns the  tensor product  of two  families <f> and <g>; the
Fourier  matrix is the Kronecker  product of the matrices  for <f> and <g>,
and the eigenvalues of Frobenius are the pairwise products.
"""
function Base.:*(f::Family,g::Family)
# println(f,"*",g)
  arg=(f,g)
  for ff in arg
    if !haskey(ff,:charLabels) ff.charLabels=map(string,1:length(ff)) end
  end
  res=Family(Dict{Symbol,Any}())
  res.charLabels=join.(cartesian(getproperty.(arg,:charLabels)...),"\\otimes")
  res.fourierMat=kron(getproperty.(arg,:fourierMat)...)
  res.eigenvalues=map(prod,cartesian(getproperty.(arg,:eigenvalues)...))
  res.name=join(getproperty.(arg,:name),"\\otimes ")
  res.explanation="Tensor("*join(map(x->haskey(x,:explanation) ?
                                     x.explanation : "??",arg),",")*")"
  if all(haskey.(arg,:charNumbers))
    res.charNumbers=map(x->collect(Iterators.flatten(x)),
                          cartesian(getproperty.(arg,:charNumbers)...))
  end
  if all(haskey.(arg,:special))
    res.special=cart2lin(length.(arg),getproperty.(arg,:special))
    res.cospecial=cart2lin(length.(arg),map(f->get(f.prop,:cospecial,f.special),arg))
    if res.cospecial==res.special delete!(res.prop,:cospecial) end
  end
  if all(x->haskey(x,:perm) || length(x)==1,arg)
    res.perm=Perm(cartesian(map(x->1:length(x),arg)...),
      cartesian(map(x->haskey(x,:perm) ? permute(1:length(x),x.perm) : [1],arg)...))
  end
  if all(x->haskey(x,:lusztig) || length(x)==1,arg)
    res.lusztig=true
  end
  if any(haskey.(arg,:qEigen))
    res.qEigen=map(sum,cartesian(map(f->
      haskey(f,:qEigen) ? f.qEigen : zeros(Int,length(f)),arg)...))
  end
  res
end

"""
`galois(f::Family,p::Int)`

`x->galois(x,p)`  is  applied  to  the  Fourier  matrix  and eigenvalues of
Frobenius of the family.

```julia-repl
julia> f=UnipotentCharacters(complex_reflection_group(3,1,1)).families[2]
Family(0011,[4, 3, 2])
imprimitive family
label│eigen      1        2        3
─────┼───────────────────────────────
1    │  ζ₃²  √-3/3    √-3/3   -√-3/3
2    │    1  √-3/3 ζ₃²√-3/3 -ζ₃√-3/3
3    │    1 -√-3/3 -ζ₃√-3/3 ζ₃²√-3/3

julia> galois(f,-1)
Family(overline 0011,[4, 3, 2])
ComplexConjugate(imprimitive family)
label│eigen      1        2        3
─────┼───────────────────────────────
1    │   ζ₃ -√-3/3   -√-3/3    √-3/3
2    │    1 -√-3/3 -ζ₃√-3/3 ζ₃²√-3/3
3    │    1  √-3/3 ζ₃²√-3/3 -ζ₃√-3/3
```
"""
function CyclotomicNumbers.galois(f::Family,p::Int)
  f=Family(copy(f.prop))
  f.fourierMat=galois.(fourier(f),p)
  f.eigenvalues=galois.(f.eigenvalues,p)
  if haskey(f,:sh) f.sh=galois.(f.sh,p) end
  if haskey(f,:name)
    f.name=p==-1 ? "overline "*f.name : "Gal("*string(p)*","*f.name*")"
  end
  if haskey(f,:explanation)
    f.explanation=p==-1 ? "ComplexConjugate("*f.explanation*")" :
    "GaloisCyc("*string(p)*","*f.explanation*")"
  end
  f
end

"`conj(f::Family)`:   is    a    synonym    for 'galois(f,-1)'."
Base.conj(f::Family)=galois(f,-1)

"""
`Eigenvalues(f)`:  eigenvalues of Frobenius associated to <f>.

`String(f)', 'Print(f)`: give a short description of the family.
"""

"""
`permute(f::Family, p::Perm)`

returns  a copy of  `f` with the  Fourier matrix, eigenvalues of Frobenius,
`:charLabels…` permuted by `p`.

```julia-repl
julia> f=UnipotentCharacters(complex_reflection_group(3,1,1)).families[2]
Family(0011,[4, 3, 2])
imprimitive family
label│eigen      1        2        3
─────┼───────────────────────────────
1    │  ζ₃²  √-3/3    √-3/3   -√-3/3
2    │    1  √-3/3 ζ₃²√-3/3 -ζ₃√-3/3
3    │    1 -√-3/3 -ζ₃√-3/3 ζ₃²√-3/3

julia> permute(f,Perm(1,2,3))
Family(0011,[2, 4, 3])
Permuted((1,2,3),imprimitive family)
label│eigen        3      1        2
─────┼───────────────────────────────
3    │    1 ζ₃²√-3/3 -√-3/3 -ζ₃√-3/3
1    │  ζ₃²   -√-3/3  √-3/3    √-3/3
2    │    1 -ζ₃√-3/3  √-3/3 ζ₃²√-3/3
```
"""
function Perms.permute(f::Family,p::Perm)
  f=Family(copy(f.prop))
  for n in [:x,:chi,:perm,:special]
    if haskey(f,n) setproperty!(f,n,getproperty(f,n)^p) end
  end
  for n in [:charNumbers,:eigenvalues,:mellinLabels,:charLabels,:unpdeg,:fakdeg]
    if haskey(f,n) setproperty!(f,n,permute(getproperty(f,n),p)) end
  end
  for n in [:fourierMat,:mellin]
    if haskey(f,n) setproperty!(f,n,permute(getproperty(f,n),p;dims=(1,2))) end
  end
  f.explanation="Permuted("*repr(p;context=:TeX=>true)*","*f.explanation*")"
  f
end

#----------------------- now definitions of particular families -------------
chevieset(:families,:C1,
  Family(Dict(:group=>"C1", :name=>"C_1", :explanation=>"trivial",
         :charLabels=>[""], :fourierMat=>hcat([1]), :eigenvalues=>[1],
         :mellin=>[[1]],:mellinLabels=>[""])))

chevieset(:families,Symbol("C'1"),
  Family(Dict(:group=>"C1", :name=>"C'_1",
  :explanation=>"-trivial",
  :charLabels=>[""],
  :fourierMat=>[[-1]],
  :eigenvalues=>[-1],
  :sh=>[1])))

chevieset(:families,:C2,
  Family(Dict(:group=>"C2", :name=>"C_2",
  :explanation=>"DrinfeldDouble(Z/2)",
  :charLabels=>["(1,1)", "(g_2,1)", "(1,\\varepsilon)", "(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[1 1 1 1;1 1 -1 -1;1 -1 1 -1;1 -1 -1 1],
  :eigenvalues=>[1,1,1,-1],
  :perm=>Perm(),
  :mellin=>[[1,1,0,0],[1,-1,0,0],[0,0,1,1],[0,0,1,-1]],
  :mellinLabels=>["(1,1)","(1,g2)","(g2,1)","(g2,g2)"])))

chevieset(:families,Symbol("C'2"),
  Family(Dict(:group=>"C2",:name=>"C'_2",
  :explanation=>"TwistedDrinfeldDouble(Z/2)",
  :charLabels=>["(1,1)",  "(1,\\varepsilon)", "(g_2,1)","(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[1 1 -1 -1;1 1 1 1;-1 1 1 -1;-1 1 -1 1],
  :eigenvalues=>[1,1,E(4),-E(4)],
  :qEigen=>[0,0,1,1]//2,
  :perm=>Perm(3,4),
  :lusztig=>true, # does not satisfy (ST)^3=1 but (SPT)^3=1
  :cospecial=>2)))

chevieset(:families,Symbol("C'\"2"),
  Family(Dict(:group=>"C2", :name=>"C'''_2",
  :explanation=>"TwistedDrinfeldDouble(Z/2)'",
  :charLabels=>["(1,1)", "(1,\\varepsilon)", "(g_2,1)", "(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[1 1 -1 -1;1 1 1 1;-1 1 -1 1;-1 1 1 -1],
  :eigenvalues=>[1,1,E(4),-E(4)],
  :qEigen=>[0,0,1,1]//2,
  :perm=>Perm(3,4),
  :special=>1,
  :cospecial=>2)))

chevieset(:families,:S3,
  Family(Dict(:group=>"S3", :name=>"D(\\mathfrak S_3)",
  :explanation=>"Drinfeld double of \$\\mathfrak S_3\$, Lusztig's version",
  :charLabels=>[ "(1,1)", "(g_2,1)", "(g_3,1)", "(1,\\rho)", "(1,\\varepsilon)",
		"(g_2,\\varepsilon)", "(g_3,\\zeta_3)", "(g_3,\\zeta_3^2)"],
  :fourierMat=>[1  3  2  2 1  3  2  2;3  3  0  0 -3 -3  0  0;
		2  0  4 -2 2  0 -2 -2;2  0 -2  4  2  0 -2 -2;
		1 -3  2  2 1 -3  2  2;3 -3  0  0 -3  3  0  0;
		2  0 -2 -2 2  0  4 -2;2  0 -2 -2  2  0 -2  4]//6,
  :eigenvalues=>[1,1,1,1,1,-1,E(3),E(3,2)],
  :perm=>Perm(7,8),
  :lusztig=>true, # does not satisfy (ST)^3=1 but (SPT)^3=1
  :mellin=>[[1,0,0,2,1,0,0,0],[0,1,0,0,0,1,0,0],[0,0,1,0,0,0,1,1],[1,0,0,-1,1,0,
   0,0],[1,0,0,0,-1,0,0,0],[0,1,0,0,0,-1,0,0],[0,0,1,0,0,0,E(3),E(3,2)],
   [0,0,1,0,0,0,E(3,2),E(3)]],
  :mellinLabels=>["(1,1)","(g2,1)","(g3,1)","(1,g3)","(1,g2)","(g2,g2)",
                  "(g3,g3)","(g3,g3^2)"])))

chevieset(:families,:X,function(p)
    ss=combinations(0:p-1,2)
    Family(Dict(:name=>"R_{\\bbZ/$p}^{\\wedge 2}",
         :explanation=>"DoubleTaft($p)",
         :charSymbols=>ss,
         :charLabels=>map(s->repr(E(p,s[1]),context=:TeX=>true)*
             "\\!\\wedge\\!"*repr(E(p,s[2]),context=:TeX=>true),ss),
    :eigenvalues=>map(s->E(p,prod(s)),ss),
    :fourierMat=>[(E(p,transpose(i)*reverse(j))-E(p,transpose(i)*j))//p for i in ss,j in ss],
    :special=>1,:cospecial=>p-1))
   end)

function SubFamily(f::Family,ind,scal,label)
  ind=filter(i->ind(f,i),1:length(f.eigenvalues))
  res=Family(Dict{Symbol,Any}())
  res.fourierMat=f.fourierMat[ind,ind].*scal
  res.eigenvalues=f.eigenvalues[ind]
  res.charLabels=f.charLabels[ind]
  res.name="$(f.name)_{[$label]}"
  if haskey(f,:charSymbols) res.charSymbols=f.charSymbols[ind] end
  if haskey(f,:group) res.group=f.group end
  special=findfirst(==(f.special),ind)
  if special!==nothing res.special=special end
  res
end

function SubFamilyij(f::Family,i,j,scal)
  g=SubFamily(f,(f,k)->sum(f.charSymbols[k])%j==i,scal,join([i,j]))
  g.explanation="subfamily(sum(charsymbols)mod $j=$i of $(f.explanation))"
  g
end

chevieset(:families,:ExtPowCyclic,function(e,n)
  g=Family(Dict{Symbol,Any}())
  g.special=1
  g.charSymbols=combinations(0:e-1,n)
  g.charLabels=map(s->join(map(x->repr(E(e,x),context=:TeX=>true),s),
                           "\\!\\wedge\\!"), g.charSymbols)
  if iszero(e%2) g.eigenvalues=Cyc.(E(24,e-1)*map(i->E(2*e,i*i+e*i),0:e-1))
  else           g.eigenvalues=Cyc.(E(24,e-1)*map(i->E(e,div(i*i+e*i,2)),0:e-1))
  end
  diag(m)=map(i->m[i,i],axes(m,1))
  g.eigenvalues=diag(exterior_power(cat(g.eigenvalues...;dims=(1,2)),n))
  g.fourierMat=exterior_power([E(e,i*j) for i in 0:e-1, j in 0:e-1]//root(e),n)
  g.name="R(\\bbZ/$e)"
  g.explanation="character ring of Z/$e"
  if n>1 g.name*="^{\\wedge $n}"
    g.explanation=ordinal(n)*" exterior power "*g.explanation
  end
  g.eigenvalues=g.eigenvalues.//g.eigenvalues[1]
  g
end)

let f=SubFamilyij(chevieget(:families,:X)(6),1,3,1-E(3))
f.cospecial=5
chevieset(:families,:X5,f)

f=chevieget(:families,:ExtPowCyclic)(4,1)
f.fourierMat.*=-E(4)
f.eigenvalues.//=f.eigenvalues[2]
f.special=2
f.qEigen=[1,0,1,0].//2
chevieset(:families,:Z4,f)

f=chevieget(:families,:ExtPowCyclic)(9,1)
f.perm=perm"(2,9)(3,8)(4,7)(5,6)"
f.qEigen=[0,2,1,0,2,1,0,2,1].//3
#if f.eigenvalues!=map(i->E(9)^(5*i^2),0:8) error() end
chevieset(:families,:Z9,f)
end

chevieset(:families,:QZ,function(n)
  pairs=[(i,j) for i in 0:n-1 for j in 0:n-1]
  res=Family(Dict{Symbol,Any}(:name=>"D(\\bbZ/$n)"))
  res.explanation="Drinfeld double "*res.name
  res.fourierMat=[E(n,x*c1+x1*c) for (x,c) in pairs, (x1,c1) in pairs]//n
  res.eigenvalues=[E(n,x*c) for (x,c) in pairs]
  res.special=1
  res.charLabels=[sprint(print,"(",E(n,x),",",E(n,c),")";context=rio(TeX=true))
                    for (x,c) in pairs]
  res
end)

# The big family f of dihedral groups. For e=5 occurs in H3, H4
chevieset(:families,:Dihedral,function(e)
  e1=div(e,2)
# the cuspidal chars are S(k,l) where 0<k<l<e-k
  nc=[[k,l] for k in 1:e1-1 for l in k+1:e-k-1]
  if iseven(e) nc=vcat([[0,e1,1],[0,e1,-1]],map(l->[0,l],1:e1-1),nc)
# the principal series for even e are:[S(0,l) with 0<l<e1]+[S(0,e1)',S(0,e1)'']
  else nc=vcat(map(l->[0,l],1:e1),nc)
# The principal series for odd e are:[S(0,l) with 0<l<e1+1]
  end
  c=a->E(e,a)+E(e,-a)
  f=Family(Dict{Symbol,Any}())
  f.eigenvalues=map(s->E(e,-prod(s[1:2])),nc)
  f.size=length(nc)
  f.parameters=nc
  f.charLabels=map(repr,nc)
  f.name="0"^(e-2)*"11"
  f.explanation="Dihedral($e) family"
  if iseven(e)
    f.fourierMat=map(nc)do i
      map(nc)do j
        if length(i)==2
          i1,i2=i
          if length(j)==2 return (c(j'*[i2,-i1])-c(j'*[-i1,i2]))//e
          else return  ((-1)^i1-(-1)^i2)//e
          end
        elseif length(i)==3
          if length(j)==2 return ((-1)^j[1]-(-1)^j[2])//e
          elseif i==j return (1-(-1)^e1+e)//2//e
          else return (1-(-1)^e1-e)//2//e
          end
        end
      end
    end
    f.fourierMat=improve_type(toM(f.fourierMat))
    f.special=3
    f.lusztig=true
  else
# The associated symbol to S(0,l) is s_i=[0] for i≠0,l and s_0=s_l=[1].
    f.fourierMat=[(c(i'*reverse(j))-c(i'*j))//e for i in nc, j in nc]
# *(-1)^count(iszero,[i[1],j[1]])*  This sign is in
# [Malle, "Unipotente Grade", Beispiel 6.29]
    f.special=1
  end
  c=filter(function(i)
            p=findfirst(==([nc[i][1],e-nc[i][2]]),nc)
            !isnothing(p) && p>i
          end,1:length(nc))
  f.perm=prod(c;init=Perm()) do i
    Perm(i,findfirst(==([nc[i][1],e-nc[i][2]]),nc))
  end
  f
end)

"""
`drinfeld_double(g[,opt])`

Given  a (usually small) finite group  `Γ`, Lusztig has associated a family
(a  Fourier matrix, a list of eigenvalues of Frobenius) which describes the
representation ring of the Drinfeld double of the group algebra of `Γ`, and
for   some  appropriate  small  groups  describes  a  family  of  unipotent
characters. We do not explain the details of this construction, but explain
how its final result building Lusztig's Fourier matrix, and a variant of it
that we use in Spetses, from `Γ`.

The  elements of the family are in bijection  with the set `𝓜 (Γ)` of pairs
`(x,χ)`  taken up to  `Γ`-conjugacy, where `x∈Γ`  and `χ` is an irreducible
complex-valued   character  of  `C_Γ(x)`.  To  such  a  pair  `ρ=(x,χ)`  is
associated  an  eigenvalue  of  Frobenius  defined  by  ``ω_ρ:=χ(x)/χ(1)``.
Lusztig  then defines a Fourier matrix `S₀` whose coefficient is given, for
`ρ=(x,χ)` and `ρ'=(x', χ')`, by:

``S₀_{ρ,ρ'}:=|C_Γ(x)⁻¹|∑_{ρ₁=(x₁,χ₁)}χ̄₁(x)χ(y₁)``

where  the sum is over all pairs `ρ₁∈𝓜 (Γ)` which are `Γ`-conjugate to `ρ'`
and  such that ``y₁∈ C_Γ(x)``. This  coefficient also represents the scalar
product ``⟨ρ,ρ'⟩_{𝐆^F}`` of the corresponding unipotent characters.

A  way to  understand the  formula for  ``S₀_{ρ,ρ'}`` better is to consider
another  basis of the complex  vector space with basis  `𝓜 (Γ)`, indexed by
the  pairs  `(x,y)`  taken  up  to  `Γ`-conjugacy,  where  `x`  and `y` are
commuting  elements  of  `Γ`.  This  basis  is  called  the basis of Mellin
transforms, and given by:

``(x,y)=∑_{χ∈ Irr(C_Γ(x))}χ(y)(x,χ)``

In  the  basis  of  Mellin  transforms,  the  linear  map  `S₀` is given by
`(x,y)↦(x⁻¹,y⁻¹)`  and  the  linear  transformation  `T` which sends `ρ` to
`ω_ρρ`   becomes  `(x,y)↦(x,xy)`.   These  are   particular  cases  of  the
permutation  representation of `GL₂(ℤ)`  on the basis  of Mellin transforms
where ``\\begin{pmatrix}a&b\\cr c&d\\end{pmatrix}`` acts by
`(x,y)↦(xᵃyᵇ,xᶜyᵈ)`.

Fourier  matrices in finite reductive groups  are given by the above matrix
`S₀`.  But for non-rational Spetses, we use a different matrix `S` which in
the  basis of Mellin transforms  is given by `(x,y)↦(y⁻¹,x)`. Equivalently,
the formula ``S_{ρ,ρ'}`` differs from the formula for ``S₀_{ρ,ρ'}`` in that
there  is no complex conjugation  of `χ₁`; thus the  matrix `S` is equal to
`S₀` multiplied on the right by the permutation matrix which corresponds to
`(x,χ)↦(x,χ̄)`.  The advantage of the matrix `S` over `S₀` is that the pair
`S,T`  satisfies directly the axioms for a fusion algebra (see below); also
the matrix `S` is symmetric, while `S₀` is Hermitian.

Thus there are two variants of 'drinfeld_double`:

`drinfeld_double(g,lu=true)`

returns  a family  containing Lusztig's  Fourier matrix  `S₀`, and an extra
field  '.perm'  containing  the  permutation  of  the  indices  induced  by
`(x,χ)↦(x,χ̄)`,  which allows  to recover  `S`, as  well as  an extra field
`:lusztig', set to 'true'.

`drinfeld_double(g)`

returns a family with the matrix `S`, which does not have fields '.lusztig'
or '.perm'.

The family record 'f' returned also has the fields:

`:group`: the group `Γ`.

`:charLabels`: a list of labels describing the pairs `(x,χ)`, and thus also
specifying in which order they are taken.

`:fourierMat`: the Fourier matrix (the matrix `S` or `S₀` depending on the
call).

`:eigenvalues`: the eigenvalues of Frobenius.

`:xy`: a list of pairs '[x,y]' which are representatives of the
`Γ`-orbits of pairs of commuting elements.

`:mellinLabels`: a list of labels describing the pairs '[x,y]'.

`:mellin`:  the base change matrix between  the basis `(x,χ)` and the basis
of   Mellin  transforms,   so  that   |f.fourierMat^(f.mellin^-1)|  is  the
permutation  matrix (for `(x,y)↦(y⁻¹,x)`  or `(x,y)↦(y⁻¹,x⁻¹)` depending on
the call).

`:special`: the index of the special element, which is `(x,χ)=(1,1)`.

```julia-rep1
julia> drinfeld_double(CoxSym(3))
Family(D(CoxSym(3)):8)
   label│eigen
────────┼─────────────────────────────────────────────
(1,X.1) │    1  1/6  1/3 1/6 -3/2 -3/2  1/3  1/3  1/3
(1,X.2) │    1  1/3  2/3 1/3    0    0 -1/3 -1/3 -1/3
(1,1)   │    1  1/6  1/3 1/6  3/2  3/2  1/3  1/3  1/3
(2a,X.1)│   -1 -1/6    0 1/6  1/2 -1/2    0    0    0
(2a,1)  │    1 -1/6    0 1/6 -1/2  1/2    0    0    0
(3a,1)  │    1  1/3 -1/3 1/3    0    0  2/3 -1/3 -1/3
(3a,X.2)│  ζ₃²  1/3 -1/3 1/3    0    0 -1/3 -1/3  2/3
(3a,X.3)│   ζ₃  1/3 -1/3 1/3    0    0 -1/3  2/3 -1/3

julia> drinfeld_double(CoxSym(3);lu=true)
Family(LD(CoxSym(3)):8)
   label│eigen
────────┼─────────────────────────────────────────────
(1,X.1) │    1  1/6  1/3 1/6 -3/2 -3/2  1/3  1/3  1/3
(1,X.2) │    1  1/3  2/3 1/3    0    0 -1/3 -1/3 -1/3
(1,1)   │    1  1/6  1/3 1/6  3/2  3/2  1/3  1/3  1/3
(2a,X.1)│   -1 -1/6    0 1/6  1/2 -1/2    0    0    0
(2a,1)  │    1 -1/6    0 1/6 -1/2  1/2    0    0    0
(3a,1)  │    1  1/3 -1/3 1/3    0    0  2/3 -1/3 -1/3
(3a,X.2)│  ζ₃²  1/3 -1/3 1/3    0    0 -1/3  2/3 -1/3
(3a,X.3)│   ζ₃  1/3 -1/3 1/3    0    0 -1/3 -1/3  2/3
```
"""
function drinfeld_double(g;lu=false)
  res=Family(Dict{Symbol,Any}(:group=> g))
  res.classinfo=map(classreps(g), classnames(g))do c,n
    r=Dict{Symbol, Any}(:elt => c,:name => n)
    if isone(c) r[:name]="1" end
    r[:centralizer]=centralizer(g, c)
    r[:centelms]=classreps(r[:centralizer])
    t=CharTable(r[:centralizer])
#   println("t=$t")
    r[:charNames]=charnames(r[:centralizer]; TeX=true)
    r[:names]=t.classnames
    r[:names][findfirst(isone,r[:centelms])]="1"
    r[:chars]=t.irr
    r[:charNames][findfirst(x->all(isone,x),r[:chars])]="1"
    r[:centralizers]=t.centralizers
    r
  end
  res.charLabels=vcat(
      map(r->map(c->"($(r[:name]),$c)",r[:charNames]),res.classinfo)...)
  if isabelian(g)
    for r in res.classinfo
      r[:names]=map(x->res.classinfo[findfirst(s->s[:elt]==x,
                                           res.classinfo)][:name],r[:centelms])
    end
  end
  res.eigenvalues=vcat(map(r->
    r[:chars][:,position_class(r[:centralizer],r[:elt])].//
    r[:chars][:,position_class(r[:centralizer],one(g))],res.classinfo)...)
  if lu
    res.name="L"
    res.explanation="Lusztig's"
  else
    res.name=""
    res.explanation=""
  end
  res.name*="D($g)"
  res.explanation*="DrinfeldDouble($g)"
  res.mellin=cat(map(r->
          conj(toM(map(x->x.//r[:centralizers],eachrow(r[:chars]))))^-1,
    res.classinfo)...,dims=(1,2))
  res.mellinLabels=reduce(vcat,map(x->map(y->"($(x[:name]),$y)",x[:names]),res.classinfo))
  res.xy=reduce(vcat,map(r->map(y->[r[:elt],y], r[:centelms]),res.classinfo))
  p=vcat(map(r->map(r[:centelms])do y
    r1=res.classinfo[position_class(g, y^-1)]
    el=transporting_elt(g, y^-1, r1[:elt])
    findfirst(==([r1[:elt],r1[:centelms][position_class(r1[:centralizer],
                                             r[:elt]^el)]]),res[:xy])
  end, res.classinfo)...)
  delete!(res.prop, :classinfo)
  res.fourierMat=inv(res.mellin)*one(res.mellin)[p,:]*res.mellin
  if lu
    res.perm=Perm(conj(res.mellin),res.mellin;dims=2)
    res.fourierMat=permute(res.fourierMat, res.perm,dims=1)
  end
  res.special=findfirst(==("(1,1)"),res.charLabels)
  Family(res)
end

"""
`ndrinfeld_double(g)`

This  function returns the number of elements that the family associated to
the  Drinfeld double of the group `g` would have, without computing it. The
evident advantage is the speed.

```julia-repl
julia> Families.ndrinfeld_double(complex_reflection_group(5))
378
```
"""
ndrinfeld_double(g)=sum(c->length(classreps(centralizer(g,c))),classreps(g))

"""
`family_imprimitive(S)`

`S` should be a symbol for a unipotent characters of an imprimitive complex
reflection  group 'G(e,1,n)' or 'G(e,e,n)'. The function returns the family
containing `S`.

```julia-repl
julia> family_imprimitive([[0,1],[1],[0]])
Family(0011,3)
imprimitive family
label│eigen      1        2        3
─────┼───────────────────────────────
1    │  ζ₃²  √-3/3   -√-3/3    √-3/3
2    │    1 -√-3/3 ζ₃²√-3/3 -ζ₃√-3/3
3    │    1  √-3/3 -ζ₃√-3/3 ζ₃²√-3/3
```
"""
function family_imprimitive(S)
# we follow G. Malle, "Unipotente Grade...", J. Algebra 177 (1995)
# §4 for G(e,1,n) and §6 for G(e,e,n).
# Initial writing  GM 26.10.2000 accelerated JM 10.08.2011
  e=length(S)
  ct=sort(vcat(S...))
  Scoll=tally(ct)
  d=length(ct)%e
  if !(d in [0,1])
    error("length(",joindigits(ct),") should be 0 or 1 mod.",e," !\n")
  end
  m=div(length(ct),e)
# Fourier matrix of the family of symbols with content ct:
# Let F be the set of functions ct->0:e-1 which are injective restricted to
# a  given  value  in  ct,  with  the  sum  of their values mod. e equal to
# div(length(ct),e)*binomial(e,2). Then for f∈ F the list of preimages of f
# is  a symbol S(f). Conversely  for a symbol S  there is a 'canonical' map
# f(S) which is increasing on entries of ct of given value. Then
#
# Fourier(S,T)=∑_{f∈ F∣S(f)=S}ε(f)ε(f(T))ζₑ^{f*f(T)}
#
# where for f in F with image f.(ct)
# ε(f)=(-1)^{number of non-inversions in the list f.(ct)}
# and f*f(T) is the scalar product of vectors f.(ct) and f(T).(ct).
#
# To  compute this reasonably fast, it is  decomposed as a product of sums,
# each relative to a set of consecutive equal entries in ct. If fᵢ(S) and
# fᵢ(T) are the restrictions to elements of ct of value i (two subsets of
# 0:e-1 of length mᵢ, the multiplicity of i) then one
# does ∏ᵢ(∑_{σ∈ 𝔖 _{mᵢ}}ε(σ)ζₑ^{-σ(fᵢ(S))*fᵢ(T)})=∏ᵢ det(ζₑ.^(-fᵢ(S)*fᵢ(T)'))
  j=(m*binomial(e,2))%e # for f∈ F we must have sum(f,ct)mod e==j
  ff=filter(x->sum(x)%e==j,cartesian(map(i->0:e-1, Scoll)...))
  ff=map(ff)do coll
    map((x,y)->filter(c->sum(c)%e==y,combinations(0:e-1,x[2])),Scoll,coll)
  end
  ff=reduce(vcat,map(x->cartesian(x...), ff))
  ffc=map(x->vcat(x...),ff) # now  ffc are the "canonical" functions
  eps=map(l->(-1)^sum(i->count(l[i].<@view l[i+1:end]),eachindex(l)),ffc)
  fcdict=Dict{Tuple{Vector{Int},Vector{Int}},e<=2 ? Int : Cyc{Int}}()
  function fc(e,f1,f2) # local Fourier coefficient
    k=f1<=f2 ? (f1,f2) : (f2,f1)
    get!(fcdict,k)do
      GLinearAlgebra.det(Cyc.(E.(e,-f1*f2')))
    end
  end
  mat=eps.*map(ff)do fS
      eps.*map(fT->prod(fc.(e,fS,fT)),ff)
  end
  mat//=(E(4,binomial(e+1,2)-1)*root(e)^e)^m
  frobs=E(12,-(e^2-1)*m).*map(i->E(2e,-sum(j->sum(j.^2),i)-e*sum(sum,i)),ff)
  symbs=map(ffc)do f
    map(x->ct[x],map(x->findall(==(x),f),0:e-1))
  end
  # next signs are 1 on the principal series
  newsigns=(-1)^(binomial(e,2)*binomial(m,2))*
      map(S->(-1)^sum((0:e-1).*binomial.(length.(S),2)),symbs)
  mat=map((s,l)->s*(newsigns.*l),newsigns,mat)
  if d==0 # compact entries...
    isreducedsymb(s)=all(x->s==x || Symbols.lesssymbols(x,s),
                               circshift.(Ref(s),1:length(s)-1))
    schon=isreducedsymb.(symbs)
    mult=Int[]
    for (i,si) in pairs(symbs)
      if schon[i]
        orb=circshift.(Ref(si),1:length(si))
        f=findfirst(==(si),orb)
        push!(mult,div(e,f)) # Symmetry group
        schon[filter(j->symbs[j] in view(orb,1:f),i+1:length(symbs))].=false
      end
    end
    frobs=reduce(vcat,fill.(frobs[schon],mult))
    symbs=reduce(vcat,map((m,s)->m==1 ? [s] :
         map(j->vcat(s[1:div(e,m)], [m,j]), 0:m-1), mult, symbs[schon]))
    mat=reduce(vcat,map((m,l)->map(
       i->reduce(vcat,map((n,c)->fill((e*c)//(m*n),n),mult,l[schon])),1:m),
                         mult, mat[schon]))
    mult=vcat(fill.(mult,mult)...)
    for (i,si) in pairs(symbs), (j,sj) in pairs(symbs)
      if fullsymbol(si)==fullsymbol(sj)
        mat[i][j]-=1//mult[i]
        if si==sj mat[i][j]+=1 end
      end
    end
    if (toM(mat)*cat(Cyc.(frobs)...;dims=(1,2)))^3!=one(toM(mat))
      print("** WARNING: (S*T)^3!=1\n")
    end
  end
  Family(Dict(:symbols=>symbs,
    :fcdict=>fcdict,
    :ff=>ff,
    :fourierMat=>mat,
    :eigenvalues=>frobs,
    :name=>joindigits(ct),
    :explanation=>"imprimitive family",
    :special=>1,
    :charLabels=>string.(1:length(symbs)), # should be improved
    :size=>length(symbs)))
end

"""
`FamiliesClassical(l)`

`l`  should be a list of symbols which classify the unipotent characters of
a  classical reductive group, like `symbols(2,r)` for type `Bᵣ` or `Cᵣ`, or
`symbols(2,r,0)`  for type `Dᵣ`. The function  returns the list of families
determined  by these symbols.
```julia-repl
julia> FamiliesClassical(symbols(2,3)) # for a reductive group of type B₃
6-element Vector{Family}:
 Family(112,[2])
 Family(022,[6])
 Family(3,[9])
 Family(01123,[1, 3, 8, 11])
 Family(0112233,[4])
 Family(013,[5, 7, 10, 12])
```
"""
function FamiliesClassical(sym)
  # for the notations see Lusztig "Characters of reductive groups over a
  # finite field" Ann. Math. studies 107, sections 4.5 and 4.6
  t=map(sym) do ST
    ST=fullsymbol(ST)
    Z1=sort(symdiff(ST...))
    D=length(Z1)%2
    M♯=sort(symdiff(setdiff(Z1, ST[2]),Z1[1+D:2:length(Z1)-1]))
    if D==1 && length(M♯)%2!=0 M♯=setdiff(Z1,M♯) end
    (Z1,M♯,content=sort(reduce(vcat,ST)))
  end
  res=reduce(vcat,
    map(collect(groupby(getproperty.(t,:content),eachindex(t))))do (k,v)
      if length(v)==2 # periodic symbols in type D
        [(content=k,charNumbers=[i],M♯=[t[i].M♯]) for i in v]
      else [(content=k,charNumbers=v,M♯=[t[i].M♯ for i in v])]
      end
    end)
  map(res)do f
    f=Family(Dict(pairs(f)))
    Z1=filter(x->count(==(x),f.content)==1,f.content)
    f.fourierMat=(1//2^(div(length(Z1)-1,2))).*
      [(-1)^length(intersect(x, y)) for x in f.M♯, y in f.M♯]
    f.eigenvalues=map(x->(-1)^div(defectsymbol(sym[x])+1,4),f.charNumbers)
    if length(f.eigenvalues)==1
      f.charLabels=[""]
      f.special=1
    else
      f.charLabels=map(f.M♯)do M
        v=map(z->count(>=(z),M)%2,Z1)
        D=length(v)
        v1=v[2:2:D-D%2]
        v2=v[3:2:(D-1)+D%2]
        if D%2==1 push!(v1,0) end
        v1=map(i->(v1[i]+v1[i+1])%2, 1:length(v2))
        s="+-"
        s[v2.+1]*","*s[v1.+1]
      end
      f.special=findfirst(x->all(y->y in "+,",x),f.charLabels)
    end
    f.name=joindigits(f.content)
    f.explanation="classical family"
    f.perm=Perm()
    f.size=length(f.charNumbers)
    f
  end
end

function Base.show(io::IO, ::MIME"text/html",f::Family)
  show(IOContext(io,:TeX=>true),"text/plain",f)
end

function Base.show(io::IO,f::Family)
  if hasdecor(io) || !haskey(f,:group)
    name=haskey(f,:name) ? f.name : "???"
    printTeX(io,"Family(\$",name,"\$")
  else print(io,"Family(",repr(f.group))
    if !haskey(f,:charNumbers) print(io,")"); return end
  end
  if haskey(f,:charNumbers) print(io,",",f.charNumbers,")")
  else print(io,",",length(f),")")
  end
end

function Base.show(io::IO,::MIME"text/plain",f::Family)
# display the labels, eigenvalues and Fourier matrix for f
  TeX=get(io,:TeX,false)
  println(io,f)
  if TeX println(io,"\\par") end
  if haskey(f,:explanation) # && f.explanation!=name
    printTeX(io,f.explanation,"\n")
  end
  row_labels=haskey(f,:charLabels) ? f.charLabels : string.(1:length(f))
  t=[repr.(f.eigenvalues;context=io)]
  col_labels=TeX ? ["\\Omega"] : ["eigen"]
  if haskey(f,:signs)
    push!(t,string.(f.signs))
    push!(col_labels,"\\mbox{signs}")
  end
  append!(t,toL(map(y->repr(y;context=io),Cyc.(fourier(f)))))
  if maximum(length.(row_labels))<=4 append!(col_labels,row_labels)
  else append!(col_labels,map(x->" ",row_labels))
  end
  showtable(io,permutedims(toM(t));row_labels,col_labels,rows_label="\\mbox{label}")
end

#------------------------ Fusion algebras -------------------------------
@GapObj struct FusionAlgebra<:FiniteDimAlgebra
  fourier::Matrix
  special::Int
  involution::SPerm{Int16}
  duality::SPerm{Int16}
  multable::Vector{Vector{Vector{Pair}}}
end

"""
    `fusion_algebra(f::Family)`
    `fusion_algebra(S,special=1)`

All  the Fourier matrices `S` in Chevie are unitary, that is `S⁻¹=conj(S)`,
and  have a  *special* line  `s` (the  line of  index `s=special(f)`  for a
family  `f`) such that no entry `Sₛ,ᵢ`  is equal to `0`. Further, they have
the  property that  the sums  `Cᵢ,ⱼ,ₖ=sumₗ Sᵢ,ₗ  Sⱼ,ₗ conj(Sₖ,ₗ)/Sₛ,ₗ` take
integral  values. Finally,  `S` has  the property  that complex conjugation
does a permutation with signs `σ` of the lines of `S`.

It  follows that we can define a `Z`-algebra `A` as follows: it has a basis
`bᵢ`  indexed by the lines of `S`,  and has a multiplication defined by the
fact that the coefficient of `bᵢbⱼ` on `bₖ` is equal to `Cᵢ,ⱼ,ₖ`.

`A`  is commutative, and has as unit  the element `bₛ`; the basis σ(bᵢ)` is
`dual to `bᵢ` for the linear form (bᵢ,bⱼ)=Cᵢ,ⱼ,σ₍ₛ₎`.

```julia-repl
julia> W=complex_reflection_group(4)
G₄

julia> uc=UnipotentCharacters(W);f=uc.families[4];

julia> A=fusion_algebra(fourier(f),1)
Fusion Algebra dim.5

julia> b=basis(A)
5-element Vector{AlgebraElt{Gapjm.Families.FusionAlgebra, Int64}}:
 B₁
 B₂
 B₃
 B₄
 B₅

julia> b*permutedims(b)
5×5 Matrix{AlgebraElt{Gapjm.Families.FusionAlgebra, Int64}}:
 B₁  B₂      B₃      B₄        B₅
 B₂  -B₄+B₅  B₁+B₄   B₂-B₃     B₃
 B₃  B₁+B₄   -B₄+B₅  -B₂+B₃    B₂
 B₄  B₂-B₃   -B₂+B₃  B₁+B₄-B₅  -B₄
 B₅  B₃      B₂      -B₄       B₁

julia> CharTable(A)
CharTable(Fusion Algebra dim.5)
 │1    2    3  4  5
─┼──────────────────
1│1  √-3 -√-3  2 -1
2│1    1    1  .  1
3│1   -1   -1  .  1
4│1    .    . -1 -1
5│1 -√-3  √-3  2 -1
```
"""
function fusion_algebra(S::Matrix,special::Int=1;opt...)
# zero=AlgebraElt(A,zero(ModuleElt{Int,T}))
# one=AlgebraElt(A,ModuleElt(special=>1))
  involution=SPerm(S,conj.(S);dims=1)
  if isnothing(involution) error("complex conjugacy is not SPerm(rows)") end
  if order(involution)>2 error("complex conjugacy is of order 4") end
  irr=mapslices(x->x.//x[special],S;dims=1)
  duality=SPerm(collect(eachrow(permute(irr,Perm(involution)))),
                 collect(eachrow(irr)))
  if isnothing(duality) error("the matrix does not have the * involution") end
  if order(duality)>2 error("duality is not an involution") end
  s=mapslices(x->x.//conj(x[special]),conj.(S);dims=1)
  d=size(S,1)
  multable=map(i->map(j->filter(x->x[2]!=0,
            map((x,y)->x=>y,1:d,s*(S[i,:].*S[j,:]))),1:i),1:d)
  if d>1 && all(r->all(c->all(p->p[2]>=0,c),r),multable)
    InfoChevie("# Algebra dim. ",d,": positive structure constants\n");
  end
  if !all(r->all(c->all(p->isinteger(p[2]),c),r),multable)
      error("structure constants are not integral")
  else multable=map(r->map(c->[k=>Int(i) for (k,i) in c],r),multable)
  end
  A=FusionAlgebra(S,special,involution,duality,multable,Dict{Symbol,Any}())
  d=map(ratio,eachcol(irr),eachcol(S)) # d=inv.(S[special,:]) ?
  if nothing in d  error() end
  A.cDim=d[special]^2
  A.qDim=d[special].//d
  A.irr=transpose(irr)
  A.charnames=haskey(opt,:charnames) ? opt[:charnames] : string.(1:dim(A))
  A.classnames=haskey(opt,:classnames) ? opt[:classnames] : string.(1:dim(A))
  A
end

function fusion_algebra(f::Family)
  get!(f,:fusion_algebra)do
  fusion_algebra(fourier(f),special(f);charnames=f.charLabels,classnames=f.charLabels)
  end
end

Algebras.dim(A::FusionAlgebra)=size(A.fourier,1)

Base.show(io::IO,A::FusionAlgebra)=print(io,"Fusion Algebra dim.",dim(A))

using LinearAlgebra: LinearAlgebra
function Algebras.idempotents(A::FusionAlgebra)
  get!(A,:idempotents)do
    LinearAlgebra.Diagonal(A.fourier[A.special,:])*
      conj.(transpose(A.fourier))*basis(A)
  end
end

Algebras.iscommutative(A::FusionAlgebra)=true

function Chars.CharTable(A::FusionAlgebra)
  irr=improve_type(toM(map(e->map(b->ratio(coefficients(b*e),coefficients(e)),
                                  basis(A)), idempotents(A))))
  if irr!=A.irr error() end
  labels=string.(1:dim(A))
  centralizers=fill(dim(A),dim(A))
  CharTable(irr,A.charnames,A.classnames,centralizers,
         dim(A),Dict{Symbol,Any}(:name=>repr(A;context=:TeX=>true)))
end

function involution(e::AlgebraElt{FusionAlgebra})
  p=e.A.involution
  AlgebraElt(e.A,ModuleElt([Int(abs(b^p))=>c*sign(b^p) for (b,c) in e.d]))
end

function duality(e::AlgebraElt{FusionAlgebra})
  p=e.A.duality
  AlgebraElt(e.A,ModuleElt([Int(abs(b^p))=>c*sign(b^p) for (b,c) in e.d]))
end

end
