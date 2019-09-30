FamilyOps=Dict()

struct Family
  prop::Dict{Symbol,Any}
end

Family(f::Family)=f

function getf(s::String)
  f=chevieget(:families,Symbol(s))
  (f isa Dict) ? Family(deepcopy(f)) : Family(deepcopy(f.prop))
end

function Family(s::String,v::AbstractVector,d::Dict=Dict{Symbol,Any}())
  f=getf(s)
  f[:charNumbers]=v
  merge!(f,d)
end

function Family(s::String,d::Dict=Dict{Symbol,Any}())
  f=getf(s)
  merge!(f,d)
end

function Family(f::Family,v::AbstractVector,d::Dict=Dict{Symbol,Any}())
  f[:charNumbers]=v
  merge!(f,d)
end

function Family(f::Dict{Symbol,Any},v::AbstractVector,d::Dict=Dict{Symbol,Any}())
  f=Family(f)
  f[:charNumbers]=v
  merge!(f,d)
end

Base.convert(::Type{Dict{Symbol,Any}},f::Family)=f.prop
Base.getindex(f::Family,k)=f.prop[k]
Base.haskey(f::Family,k)=haskey(f.prop,k)
Base.setindex!(f::Family,v,k)=setindex!(f.prop,v,k)
function Base.merge!(f::Family,d::Dict)
  merge!(f.prop,d)
  if !isa(f[:fourierMat],Matrix)
    f[:fourierMat]=hcat(f[:fourierMat]...)
  end
  if !haskey(f,:charLabels)
    f[:charLabels]=map(string,1:length(f[:eigenvalues]))
  end
  if haskey(d,:signs)
    signs=d[:signs]
    m=f[:fourierMat]
    for i in axes(m,1), j in axes(m,2)
      m[i,j]*=signs[i]*signs[j]
    end
    f[:fourierMat]=m
  end
  f
end

function Base.:*(f::Family,g::Family)
# println(f,"*",g)
  arg=(f,g)
  for ff in arg
    if !haskey(ff,:charLabels) 
      ff[:charLabels]=map(string,eachindex(ff[:eigenvalues])) 
    end
  end
  res=Dict{Symbol,Any}(
    :charLabels=>join.(Cartesian(getindex.(arg,:charLabels)...),"\\otimes"),
    :fourierMat=>kron(getindex.(arg,:fourierMat)...),
    :eigenvalues=>map(prod,Cartesian(getindex.(arg,:eigenvalues)...)),
    :name=>join(getindex.(arg,:name),"\\otimes "),
    :explanation=>"Tensor("*join(getindex.(arg,:explanation),",")*")"
  )
  if all(haskey.(arg,:charNumbers))
    res[:charNumbers]=Cartesian(getindex.(arg,:charNumbers)...)
  end
#  if all(haskey.(arg,:special))
#    res.special:=PositionCartesian(List(arg,Size),getindex.(arg,:special));
#    res.cospecial:=PositionCartesian(List(arg,Size),
#      List(arg,function(f)if IsBound(f.cospecial) then return f.cospecial;
#                                          else return f.special;fi;end));
#    if res.cospecial=res.special then Unbind(res.cospecial);fi;
#  fi;
#  if ForAll(arg,f->IsBound(f.perm) or Size(f)=1) then 
#    res.perm:=PermListList(Cartesian(List(arg,x->[1..Size(x)])),
#        Cartesian(List(arg,function(x)if IsBound(x.perm) then return
#          Permuted([1..Size(x)],x.perm);else return [1];fi;end)));
#  fi;
#  if ForAll(arg,f->IsBound(f.lusztig) or Size(f)=1) then 
#    res.lusztig:=true;
#  fi;
#  if any(haskey.(arg,:qEigen))
#    res.qEigen:=List(Cartesian(List(arg,function(f) 
#      if IsBound(f.qEigen) then return f.qEigen;else return f.eigenvalues*0;fi;
#      end)),Sum);
#  fi;
  res
end

chevieset(:families,:C1,
  Family(Dict(:group=>"C1", :name=>"C_1", :explanation=>"trivial",
         :charLabels=>[""], :fourierMat=>hcat([1]), :eigenvalues=>[1],
         :mellin=>[[1]],:mellinLabels=>[""])))

chevieset(:families,:C2,
  Family(Dict(:group=>"C2", :name=>"C_2",
  :explanation=>"DrinfeldDouble(Z/2)",
  :charLabels=>["(1,1)", "(g_2,1)", "(1,\\varepsilon)", "(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[1 1 1 1;1 1 -1 -1;1 -1 1 -1;1 -1 -1 1],
  :eigenvalues=>[1,1,1,-1],
  :perm=>(),
  :mellin=>[[1,1,0,0],[1,-1,0,0],[0,0,1,1],[0,0,1,-1]],
  :mellinLabels=>["(1,1)","(1,g2)","(g2,1)","(g2,g2)"])))

chevieset(:families,:S3,
  Family(Dict(:group=>"S3", :name=>"D(S_3)",
  :explanation=>"Drinfeld double of S3, Lusztig's version",
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
    Family(Dict(:name=>"R_{\\BZ/$p}^{\\wedge 2}",
         :explanation=>"DoubleTaft($p)",
         :charSymbols=>ss,
         :charLabels=>map(s->repr(E(p)^s[1],context=:TeX=>true)*
      "\\!\\wedge\\!"*repr(E(p)^s[2],context=:TeX=>true),ss),
    :eigenvalues=>map(s->E(p)^Product(s),ss),
    :fourierMat=>[(E(p)^(i*reverse(j))-E(p)^(i*j))/p for i in ss,j in ss],
    :special=>1,:cospecial=>p-1))
   end)

chevieset(:families,Symbol("C'\"2"),
  Family(Dict(:group=>"C2", :name=>"C'''_2",
  :explanation=>"TwistedDrinfeldDouble(Z/2)'",
  :charLabels=>["(1,1)", "(1,\\varepsilon)", "(g_2,1)", "(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[1 1 -1 -1;1 1 1 1;-1 1 -1 1;-1 1 1 -1],
  :eigenvalues=>[1,1,E(4),-E(4)],
  :qEigen=>[0,0,1,1]//2,
  :perm=>Perm(3,4),
  :cospecial=>2)))

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

function SubFamily(f,ind,scal,label)
  ind=filter(i->ind(f,i),1:length(f[:eigenvalues]))
  res=Dict(:fourierMat=>f[:fourierMat][ind,ind]*scal,
           :eigenvalues=>f[:eigenvalues][ind],
           :charLabels=>f[:charLabels][ind],
   :operations=>FamilyOps)
  res[:name]="$(f[:name])_{[$label]}"
  if haskey(f,:charSymbols) res[:charSymbols]=f[:charSymbols][ind] end
  if haskey(f,:group) 
    res[:group]=f[:group] 
  end
  if f[:special] in ind 
   res[:special]=findfirst(isequal(ind),f[:special]) 
  end
  Family(res)
end

function SubFamilyij(f,i,j,scal)
  g=SubFamily(f,(f,k)->sum(f[:charSymbols][k])%j==i,scal,join([i,j]))
  g[:explanation]="subfamily(sum(charsymbols)mod $j=$i of $(f[:explanation]))"
  g
end

chevieset(:families,:ExtPowCyclic,function(e,n)
  g=Dict{Symbol,Any}(
    :special=>1,
    :operations=>FamilyOps,
    :charSymbols=>combinations(0:e-1,n)
  )
  g[:charLabels]=map(s->join(map(x->repr(E(e)^x,context=:TeX=>true),s),
                             "\\!\\wedge\\!"), g[:charSymbols])
  if iszero(e%2)
    g[:eigenvalues]=E(24)^(e-1)*map(i->E(2*e)^(i*i+e*i),0:e-1)
  else
    g[:eigenvalues]=E(24)^(e-1)*map(i->E(e)^div(i*i+e*i,2),0:e-1)
  end
  g[:eigenvalues]=DiagonalOfMat(exterior_power(DiagonalMat(g[:eigenvalues]...),n))
  g[:fourierMat]=exterior_power([E(e)^(i*j) for i in 0:e-1, j in 0:e-1]/ER(e),n)
  if n>1 g[:name]="R(\\BZ/$e)^{\\wedge $n}"
    g[:explanation]=ordinal(n)*" exterior power of char. ring of Z/$e"
  else g[:name]="R(\\BZ/$e)"
       g[:explanation]="character ring of Z/$e"
  end
  g[:eigenvalues]=g[:eigenvalues]/g[:eigenvalues][1]
  Family(g)
end)

chevieset(:families,:X5,SubFamilyij(chevieget(:families,:X)(6),1,3,1-E(3)))
CHEVIE[:families][:X5][:cospecial]=5
chevieset(:families,:Z4,chevieget(:families,:ExtPowCyclic)(4,1))
CHEVIE[:families][:Z4][:fourierMat]*=-E(4)
CHEVIE[:families][:Z4][:eigenvalues]/=chevieget(:families,:Z4)[:eigenvalues][2]
CHEVIE[:families][:Z4][:special]=2
CHEVIE[:families][:Z4][:qEigen]=[1,0,1,0]//2

chevieset(:families,:Z9,chevieget(:families,:ExtPowCyclic)(9,1))
#if CHEVIE.families.Z9.eigenvalues!=List([0..8],i->E(9)^(5*i^2))then Error();fi;
CHEVIE[:families][:Z9][:perm]=perm"(2,9)(3,8)(4,7)(5,6)"
CHEVIE[:families][:Z9][:qEigen]=[0,2/3,1/3,0,2/3,1/3,0,2/3,1/3]

chevieset(:families,:QZ,function(n)
  pairs=[(i,j) for i in 0:n-1 for j in 0:n-1]
  res=Dict{Symbol,Any}(:name=>"D(\\BZ/$n)")
  res[:explanation]="Drinfeld double "*res[:name]
  res[:fourierMat]=[E(n,x)^c1*E(n,x1)^c for (x,c) in pairs, (x1,c1) in pairs]
  res[:eigenvalues]=[E(n,x)^c for (x,c) in pairs]
  res[:special]=1
  res[:charLabels]=["($(E(n,x)),$(E(n,c)))" for (x,c) in pairs]
  Family(res)
end)

# The big family in dihedral groups. For e=5 occurs in H3, H4
chevieset(:families,:Dihedral,function(e)
  e1=div(e,2)
# the cuspidal chars are S(k,l) where 0<k<l<e-k
  nc=[[k,l] for k in 1:e1-1 for l in k+1:e-k-1]
  if iszero(e%2) nc=vcat([[0,e1,1],[0,e1,-1]],map(l->[0,l],1:e1-1),nc)
# the principal series chars in f are:[S(0,l) with 0<l<e1]+[S(0,e1)',S(0,e1)'']
  else nc=vcat(map(l->[0,l],1:e1),nc)
# The principal series chars in f are:[S(0,l) with 0<l<e1+1]
  end
  c=a->E(e)^a+E(e)^(-a)
  f=Dict( :eigenvalues=>map(s->E(e)^-prod(s[1:2]),nc),
    :size=>length(nc),
    :parameters=>nc,
    :charLabels=>map(repr,nc),
    :name=>"0"^(e-2)*"11",
    :explanation=>"Dihedral($e) family",
    :operations=>FamilyOps)
  if iszero(e%2)
    f[:fourierMat]=map(i->map(function(j)
       if length(i)==2 
         if length(j)==2 return (c(j*[i[2],-i[1]])-c(j*[-i[1],i[2]]))/e
         else return  ((-1)^i[1]-(-1)^i[2])/e
         end
       elseif length(i)==3 
         if length(j)==2 return ((-1)^j[1]-(-1)^j[2])/e
         elseif i==j return (1-(-1)^e1+e)/2/e
         else return (1-(-1)^e1-e)/2/e
         end
         end end,nc),nc)
    f[:special]=3
    f[:lusztig]=true
  else
# The associated symbol to S(0,l) is s_i=[0] for i\ne 0,l and s_0=s_l=[1].
    f[:fourierMat]=map(i->map(j->
# (-1)^Number([i[1],j[1]],x->x=0)*  This sign is in
# [Malle, "Unipotente Grade", Beispiel 6.29]
          (c(i[1]*j[2]+i[2]*j[1])-c(i[1]*j[1]+i[2]*j[2]))/e,nc),nc)
    f[:special]=1
  end
  c=filter(function(i)
            p=findfirst(isequal([nc[i][1],e-nc[i][2]]),nc)
            return !isnothing(p) && p>i
          end,1:length(nc))
  f[:perm]=prod(c) do i
   Perm(i,findfirst(isequal([nc[i][1],e-nc[i][2]]),nc))
   end
   f[:fourierMat]=hcat(f[:fourierMat]...)
   Family(f)
end)

"""
`DrinfeldDouble(<g>[,<opt>])`

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
associated  an eigenvalue of Frobenius defined by `ω_ρ:=χ(x)/χ(1)`. Lusztig
then defines a Fourier matrix `T` whose coefficient is given, for `ρ=(x,χ)`
and `ρ'=(x', χ')`, by:

`T_{ρ,ρ'}:=#C_Γ(x)⁻¹ ∑_{ρ_1=(x_1,χ_1)}χ̄_1(x)χ(y_1)`

where the sum is over all pairs `ρ_1∈𝓜 (Γ)` which are `Γ`-conjugate to `ρ'`
and  such that `y_1∈  C_Γ(x)`. This coefficient  also represents the scalar
product `⟨ρ,ρ'⟩_{𝐆^F}` of the corresponding unipotent characters.

A  way  to  understand  the  formula  for  `T_{ρ,ρ'}` better is to consider
another  basis of the complex  vector space with basis  `𝓜 (Γ)`, indexed by
the  pairs  `(x,y)`  taken  up  to  `Γ`-conjugacy,  where  `x`  and `y` are
commuting  elements  of  `Γ`.  This  basis  is  called  the basis of Mellin
transforms, and given by:

`(x,y)=∑_{χ∈ Irr(C_Γ(x))}χ(y)(x,χ)`

In  the  basis  of  Mellin  transforms,  the  linear  map  `T`  is given by
`(x,y)↦(x⁻¹,y⁻¹)`  and  the  linear  transformation which sends
`ρ`   to  `ω_ρρ`  becomes   `(x,y)↦(x,xy)`.  These  are
particular  cases of the  permutation representation of  `GL_2(ℤ)` on the
basis of Mellin transforms where
`(begin{array}{cc}a&b;c&d{array})
%begin{pmatrix}{cc}a&b;c&d{pmatrix}`
acts by `(x,y)↦(x^ay^b,x^cy^d)`.

Fourier  matrices in finite reductive groups  are given by the above matrix
`T`.  But for non-rational Spetses, we use  a different matrix `S` which in
the  basis of Mellin transforms  is given by `(x,y)↦(y⁻¹,x)`. Equivalently,
the  formula `S_{ρ,ρ'}`  differs from  the formula  for `T_{ρ,ρ'}`  in that
there  is no complex conjugation of `χ_1`;  thus the matrix `S` is equal to
`T`  multiplied on the right by the permutation matrix which corresponds to
`(x,χ)↦(x,χ̄)`.  The advantage of the matrix `S`  over `T` is that the pair
`S,Ω`  satisfies directly the axioms for a fusion algebra (see below); also
the matrix `S` is symmetric, while `T` is Hermitian.

Thus there are two variants of 'DrinfeldDouble`:

`DrinfeldDouble(<g>,rec(lusztig:=true))`

returns  a family  containing Lusztig's  Fourier matrix  `T`, and  an extra
field  '.perm'  containing  the  permutation  of  the  indices  induced  by
`(x,χ)↦(x,χ̄)`,  which allows  to recover  `S`, as  well as  an extra field
`:lusztig', set to 'true'.

`DrinfeldDouble(<g>)`

returns a family with the matrix `S`, which does not have fields '.lusztig'
or '.perm'.

The family record 'f' returned also has the fields:

`:group`: the group `Γ`.

`:charLabels`: a list of labels describing the pairs `(x,χ)`, and thus also
specifying in which order they are taken.

`:fourierMat`: the Fourier matrix (the matrix `S` or `T` depending on the
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

|    gap> f:=DrinfeldDouble(SymmetricGroup(3));
    Family("D(Group((1,3),(2,3)))")
    gap> Display(f);
    D(Group((1,3),(2,3)))
       label |eigen
    _______________________________________________________
    (1,1)    |    1 1/6  1/6  1/3  1/2  1/2  1/3  1/3  1/3
    (1,X.2)  |    1 1/6  1/6  1/3 -1/2 -1/2  1/3  1/3  1/3
    (1,X.3)  |    1 1/3  1/3  2/3    0    0 -1/3 -1/3 -1/3
    (2a,1)   |    1 1/2 -1/2    0  1/2 -1/2    0    0    0
    (2a,X.2) |   -1 1/2 -1/2    0 -1/2  1/2    0    0    0
    (3a,1)   |    1 1/3  1/3 -1/3    0    0  2/3 -1/3 -1/3
    (3a,X.2) |   E3 1/3  1/3 -1/3    0    0 -1/3 -1/3  2/3
    (3a,X.3) | E3^2 1/3  1/3 -1/3    0    0 -1/3  2/3 -1/3
    gap> f:=DrinfeldDouble(SymmetricGroup(3),rec(lusztig:=true));
    Family("LD(Group((1,3),(2,3)))")
    gap> Display(f);
    LD(Group((1,3),(2,3)))
       label |eigen
    _______________________________________________________
    (1,1)    |    1 1/6  1/6  1/3  1/2  1/2  1/3  1/3  1/3
    (1,X.2)  |    1 1/6  1/6  1/3 -1/2 -1/2  1/3  1/3  1/3
    (1,X.3)  |    1 1/3  1/3  2/3    0    0 -1/3 -1/3 -1/3
    (2a,1)   |    1 1/2 -1/2    0  1/2 -1/2    0    0    0
    (2a,X.2) |   -1 1/2 -1/2    0 -1/2  1/2    0    0    0
    (3a,1)   |    1 1/3  1/3 -1/3    0    0  2/3 -1/3 -1/3
    (3a,X.2) |   E3 1/3  1/3 -1/3    0    0 -1/3  2/3 -1/3
    (3a,X.3) | E3^2 1/3  1/3 -1/3    0    0 -1/3 -1/3  2/3|
"""
function Drinfeld_double(g;lu=false)
  res=Dict{Symbol,Any}(:group=> g)
  res[:classinfo] = map(function (c, n) local r, t
    r = Dict{Symbol, Any}(:elt => Representative(c), :name => n)
    if r[:elt] == g[:identity] r[:name] = "1" end
    r[:centralizer] = Centralizer(g, r[:elt])
    r[:centelms] = class_reps(r[:centralizer])
    t = CharTable(r[:centralizer])
    r[:charNames] = charnames(r[:centralizer]; TeX = true)
    r[:names] = ClassNamesCharTable(t)
    r[:names][Position(r[:centelms], g[:identity])] = "1"
    r[:chars] = t[:irreducibles]
    r[:charNames][PositionProperty(r[:chars], y->y==1)] = "1"
    r[:centralizers] = t[:centralizers]
    return r
end, ConjugacyClasses(g), ClassNamesCharTable(CharTable(g)))
  res[:charLabels]=vcat(map(r->map(y->"($(r[:name]),$y)",r[:charNames]),
                              res[:classinfo])...)
  if IsAbelian(g)
    for r in res[:classinfo]
      r[:names]=map(x->First(res[:classinfo],s->s[:elt]==x)[:name],r[:centelms])
    end
  end
  res[:size] = length(res[:charLabels])
  res[:eigenvalues]=vcat(map(function(r)ct=TransposedMat(r[:chars])
     return map((x, y)->x//y,ct[PositionClass(r[:centralizer],r[:elt])], 
               ct[PositionClass(r[:centralizer],g[:identity])]) end,
      res[:classinfo])...)
  if lu
      res[:name] = "L"
      res[:explanation] = "Lusztig's"
  else
      res[:name] = ""
      res[:explanation] = ""
  end
  res[:name] *= "D($g)"
  res[:explanation] *= "DrinfeldDouble($g)"
  res[:operations] = FamilyOps
  res[:mellin] = DiagonalMat(List, res[:classinfo], (r->begin
    conj(map(x->map((x, y)->x//y,x,r[:centralizers]),r[:chars]))^ -1
              end))
  res[:mellinLabels]=vcat(map(x->map(y->"($(x[:name]),$y)",x[:names]),res[:classinfo])...)
  res[:xy] = vcat(map(r->map(y->[r[:elt],y], r[:centelms]),res[:classinfo])...)
  p = vcat(map((r->begin map(function (y)
           r1 = (res[:classinfo])[PositionClass(g, y ^ -1)]
          return Position(res[:xy], [r1[:elt], (r1[:centelms])[PositionClass(r1[:centralizer], r[:elt] ^ RepresentativeOperation(g, y ^ -1, r1[:elt]))]])
                          end, r[:centelms]) end), res[:classinfo])...)
  delete!(res, :classinfo)
  res[:fourierMat] = (IdentityMat(res[:size]))[p] ^ res[:mellin]
  if lu
    res[:perm]=PermListList(conj(transpose(res[:mellin])),
                                 transpose(res[:mellin]))
    res[:fourierMat]=Permuted(res[:fourierMat], res[:perm])
  end
  res[:special] = Position(res[:charLabels], "(1,1)")
  res
end

"""
`FamilyImprimitive(<S>)`

<S> should be a symbol for a unipotent characters of an imprimitive complex
reflection  group 'G(e,1,n)' or 'G(e,e,n)'. The function returns the family

```julia-repl
julia> HasType.Family(FamilyImprimitive([[0,1],[1],[0]]))
Family(0011:3)
label│eigen      1         2         3
─────┼─────────────────────────────────
1    │  ζ₃²  √-3/3    -√-3/3     √-3/3
2    │    1 -√-3/3 (3-√-3)/6 (3+√-3)/6
3    │    1  √-3/3 (3+√-3)/6 (3-√-3)/6
```
"""
FamilyImprimitive = function (S)
# println("S=$S")
  e = length(S)
  Scoll = Collected(vcat(S...))
  ct = vcat(map(x->fill(x[1],x[2]), Scoll)...)
  d = length(ct) % e
  if !(d in [0, 1]) error("Length(", joindigits(ct), ") should be 0 or 1  %  ", e, " !\n")
        end
  m = div(length(ct) - d, e)
  j = (m * binomial(e, 2)) % e
  ll = Cartesian(map(i->0:e-1, Scoll)...)
  ll = filter(x->sum(x)%e==j,ll)
  ll = map(c->map((x,y)->filter(c->sum(c)%e==y,
                   collect(combinations(0:e-1,x[2]))),Scoll,c), ll)
  nrSymbols=sum(x->prod(length,x),ll)
  ll = vcat(map(x->Cartesian(x...), ll)...)
  eps = l->(-1)^sum(i->count(j->l[i]<j,l[i+1:end]),1:length(l))
  equiv = map(x->
      Dict(:globaleps=>length(x)==1 ? 1 :
     (-1)^sum(i->sum(j->sum(y->count(k->y<k,j),x[i]),x[i+1:end]),1:length(x)-1),
     :aa=>map(y->map(x->(l=x,eps=eps(x)),arrangements(y, length(y))),x)), ll)
  epsreps = map(x->eps(vcat(x...)), ll)
  roots = map(i->E(e,i),0:e-1)
  mat = map(i->i[:globaleps]*map(k->epsreps[k]*
  prod(l->sum(j->j.eps*roots[1+mod(-sum(map((a,b)->a*b,j.l,ll[k][l])),e)],
              i[:aa][l]),
          1:length(i[:aa])), 1:nrSymbols), equiv)
  mat = ((-1)^(m*(e-1))*mat)//(E(4,binomial(e-1,2))*ER(e)^e)^m
  frobs = E(12,-(e^2-1)*m)*map(i->E(2e,-(sum(j->j*j,i))-e*sum(sum,i)),ll)
  symbs = map(function (l)local sy, j
              sy = map(j->Int[], 1:e)
              map((v,c)->begin push!(sy[v + 1], c)
                      return 1 end, vcat(l...), ct)
              return sy
          end, ll)
  newsigns = (-1) ^ (binomial(e, 2) * binomial(m, 2)) * map(i->
               (-1)^((0:e - 1)*map(x->binomial(length(x), 2), i)),symbs)
  mat = map((s,l)->s * map((x, y)->x*y, newsigns,l),newsigns,mat)
  if d == 0
  IsReducedSymbol(s)=all(x->s==x || LessSymbols(x, s),Rotations(s)[2:length(s)])
    schon = map(IsReducedSymbol, symbs)
    mult = []
    for i = 1:nrSymbols
        if schon[i]
            orb = gapSet(Rotations(symbs[i]))
            push!(mult, e // length(orb))
            for j = Filtered(i + 1:nrSymbols, (j->symbs[j] in orb))
                schon[j] = false
            end
        end
    end
    frobs = vcat(map((m,f)->fill(f,m), mult, ListBlist(frobs, schon))...)
    symbs = vcat(map(function (m, s)
                    if m==1 return [s]
                    else return map(j->vcat(s[1:e//m], [m, j]), 0:m-1)
                    end
                end, mult, ListBlist(symbs, schon))...)
    mat = vcat(map(function (m, l)return map((i->begin
      vcat(map((n,c)->fill(e*c//m//n,n),mult,ListBlist(l, schon))...)
                         end), 1:m)end, mult, ListBlist(mat, schon))...)
    mult=vcat(map(m->fill(m,m),mult)...)
    nrSymbols=length(symbs)
    for i=1:nrSymbols
      for j=1:nrSymbols
        if FullSymbol(symbs[i])==FullSymbol(symbs[j])
            mat[i][j]-=1//mult[i]
            if symbs[i]==symbs[j] mat[i][j]+=1 end
        end
      end
    end
    if (mat*cat(frobs...,dim=(1,2)))^3!=mat^0
        print("** WARNING: (S*T)^3!=1\n")
    end
  end
  res=Dict{Symbol,Any}(:symbols=>symbs,
    :fourierMat=>hcat(mat...),
    :eigenvalues=>frobs,
    :name=>joindigits(ct),
    :explanation=>"classical family",
    :special=>1,
    :operations=>FamilyOps)
  res[:charLabels] = map(string, 1:length(res[:symbols]))
  res[:size] = length(res[:symbols])
  res
end

MakeFamilyImprimitive = function (S, uc)
  f=x->Position(uc[:charSymbols],x)
  if length(S)==1 return Family("C1", map(f, S)) end
  r = Family(FamilyImprimitive(FullSymbol(S[1])))
  r[:charNumbers] = map(f, r[:symbols])
  r[:special] = PositionProperty(r[:charNumbers],(x->uc[:a][x] == uc[:b][x]))
  r[:cospecial] = PositionProperty(r[:charNumbers],(x->uc[:A][x] == uc[:B][x]))
# if length(blocks(r[:fourierMat])) > 1 error() end
  Family(r)
end

"""
`FamiliesClassical(l)`

The  list  `l`  should  be  a  list  of symbols as returned by the function
`Symbols`,  which classify the unipotent characters of groups of type `:B`,
`:C`  or `:D`. `FamiliesClassical` returns  the list of families determined
by these symbols.

```julia-repl
julia> HasType.FamiliesClassical(HasType.BDSymbols(3,1))
6-element Array{Family,1}:
 Family(0112233:[4])    
 Family(01123:[1, 3, 8])
 Family(013:[5, 7, 10]) 
 Family(022:[6])        
 Family(112:[2])        
 Family(3:[9])          
```
The  above example shows the families of unipotent characters for the group
`B_3`.
"""
FamiliesClassical=function(sym)
  t=map(sym) do ST
    ST=fullsymbol(ST)
    f=Dict{Symbol, Any}(:Z1 => symdiff(ST...))
    D=length(f[:Z1]) % 2
    f[:M♯] = symdiff(setdiff(f[:Z1], ST[2]), f[:Z1][1+D:2:length(f[:Z1])-1])
    if D==1 && length(f[:M♯])%2!=0 f[:M♯]=setdiff(f[:Z1],f[:M♯]) end
    f[:content] = sort(vcat(ST...))
    f
  end
  res = []
  for l = CollectBy(1:length(t), i->t[i][:content])
      f = Dict{Symbol, Any}(:content => t[l[1]][:content], :charNumbers => l)
      f[:M♯] = map(x->x[:M♯],t[l])
      if length(l) == 2
          push!(res, Dict{Symbol, Any}(:content => f[:content], 
            :charNumbers => [l[2]],
            :M♯ => [f[:M♯][2]]))
          f=Dict{Symbol, Any}(:content => f[:content], :charNumbers => [l[1]], 
                              :M♯ => [f[:M♯][1]])
      end
      push!(res, f)
  end
  map(res)do f
   Z1 = filter(x->count(isequal(x),f[:content])==1,f[:content])
   f[:fourierMat] = (2//1)^(-div(length(Z1)-1,2))*map(x->
                    map(y->(-1)^length(intersect(x, y)), f[:M♯]), f[:M♯])
   f[:fourierMat]=hcat(f[:fourierMat]...)
    f[:eigenvalues] = map(x->(-1) ^ div(DefectSymbol(sym[x])+1,4), f[:charNumbers])
    if length(f[:eigenvalues]) == 1
        f[:charLabels] = [""]
        f[:special] = 1
    else
        f[:charLabels] = map(f[:M♯])do M
          v = map(z->count(y->y>=z, M) % 2, Z1)
          D = length(v)
          v1 = v[2:2:D-D%2]
          v2 = v[3:2:(D-1)+D%2]
          if D%2==1 push!(v1,0) end
          v1 = convert(Vector{Int},map(i->sum(v1[[i,i+1]]) % 2, 1:length(v2)))
          s = "+-"
          s[v2+1]*","*s[v1+1]
        end
        f[:special] = findfirst(x->all(y->y in "+,",x),f[:charLabels])
    end
    f[:name] = joindigits(f[:content])
    f[:explanation] = "classical family"
    f[:perm] = Perm()
    f[:size] = length(f[:charNumbers])
    Family(f)
#   f[:operations] = FamilyOps
  end
end

function Base.show(io::IO,f::Family)
  TeX=get(io,:TeX,false)
  repl=get(io,:limit,false)
  deep=get(io,:typeinfo,false)!=false
  if haskey(f,:name)
    name=TeX ? "\$"*f[:name]*"\$" : TeXstrip(f[:name])
  else name="???"
  end
  print(io,"Family($name:")
  if haskey(f,:charNumbers) print(io,f[:charNumbers],")")
  else print(io,length(f[:eigenvalues]),")")
  end
  if !(repl || TeX) || deep return end
  if haskey(f,:charLabels) rowLabels=f[:charLabels]
    if !TeX rowLabels=TeXstrip.(rowLabels) end
  else  rowLabels=1:length(f)
  end
  print(io,"\n")
  t=[sprint.(show,f[:eigenvalues];context=io)]
  col_labels=TeX ? ["\\Omega"] : ["eigen"]
  if haskey(f,:signs) 
   push!(t,string.(f[:signs]))
    push!(col_labels,"signs")
  end
  append!(t,toL(map(y->sprint(show,y;context=io),f[:fourierMat])))
  if maximum(length.(rowLabels))<=4 append!(col_labels,rowLabels)
  else append!(col_labels,map(x->" ",rowLabels))
  end
  format(io,permutedims(toM(t)),row_labels=rowLabels,
        col_labels=col_labels,
        rows_label=TeX ? "\\hbox{label}" : "label")
end
