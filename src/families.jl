FamilyOps=Dict()
CHEVIE[:families]=Dict(:C1=>
        Dict(:group=>"C1", :name=>"C_1", :explanation=>"trivial",
  :charLabels=>[""], :fourierMat=>[[1]], :eigenvalues=>[1],
  :mellin=>[[1]],:mellinLabels=>[""]),
  :C2=>Dict(:group=>"C2", :name=>"C_2",
  :explanation=>"DrinfeldDouble(Z/2)",
  :charLabels=>["(1,1)", "(g_2,1)", "(1,\\varepsilon)", "(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]],
  :eigenvalues=>[1,1,1,-1],
  :perm=>(),
  :mellin=>[[1,1,0,0],[1,-1,0,0],[0,0,1,1],[0,0,1,-1]],
  :mellinLabels=>["(1,1)","(1,g2)","(g2,1)","(g2,g2)"]),
  :S3=>Dict(:group=>"S3", :name=>"D(S_3)",
  :explanation=>"Drinfeld double of S3, Lusztig's version",
  :charLabels=>[ "(1,1)", "(g_2,1)", "(g_3,1)", "(1,\\rho)", "(1,\\varepsilon)",
		"(g_2,\\varepsilon)", "(g_3,\\zeta_3)", "(g_3,\\zeta_3^2)"],
  :fourierMat=>[[1, 3, 2, 2,1, 3, 2, 2],[3, 3, 0, 0,-3,-3, 0, 0],
		[2, 0, 4,-2,2, 0,-2,-2],[2, 0,-2, 4, 2, 0,-2,-2],
		[1,-3, 2, 2,1,-3, 2, 2],[3,-3, 0, 0,-3, 3, 0, 0],
		[2, 0,-2,-2,2, 0, 4,-2],[2, 0,-2,-2, 2, 0,-2, 4]]//6,
  :eigenvalues=>[1,1,1,1,1,-1,E(3),E(3,2)],
  :perm=>Perm(7,8),
  :lusztig=>true, # does not satisfy (ST)^3=1 but (SPT)^3=1
  :mellin=>[[1,0,0,2,1,0,0,0],[0,1,0,0,0,1,0,0],[0,0,1,0,0,0,1,1],[1,0,0,-1,1,0,
   0,0],[1,0,0,0,-1,0,0,0],[0,1,0,0,0,-1,0,0],[0,0,1,0,0,0,E(3),E(3,2)],
   [0,0,1,0,0,0,E(3,2),E(3)]],
  :mellinLabels=>["(1,1)","(g2,1)","(g3,1)","(1,g3)","(1,g2)","(g2,g2)",
                 "(g3,g3)","(g3,g3^2)"]),
  :X=>function(p)
    ss=combinations(0:p-1,2)
    Dict(:name=>"R_{\\BZ/$p}^{\\wedge 2}",
         :explanation=>"DoubleTaft($p)",
         :charSymbols=>ss,
         :charLabels=>map(s->repr(E(p)^s[1],context=:TeX=>true)*
      "\\!\\wedge\\!"*repr(E(p)^s[2],context=:TeX=>true),ss),
    :eigenvalues=>map(s->E(p)^Product(s),ss),
    :fourierMat=>map(i->map(j->(E(p)^(i*reverse(j))-E(p)^(i*j))/p,ss),ss),
    :special=>1,:cospecial=>p-1)
   end,
   Symbol("C'\"2")=>Dict(:group=>"C2", :name=>"C'''_2",
  :explanation=>"TwistedDrinfeldDouble(Z/2)'",
  :charLabels=>["(1,1)", "(1,\\varepsilon)", "(g_2,1)", "(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[[1,1,-1,-1],[1,1,1,1],[-1,1,-1,1],[-1,1,1,-1]],
  :eigenvalues=>[1,1,E(4),-E(4)],
  :qEigen=>[0,0,1//2,1//2],
  :perm=>Perm(3,4),
  :cospecial=>2),
   Symbol("C'2")=>Dict(:group=>"C2",:name=>"C'_2",
  :explanation=>"TwistedDrinfeldDouble(Z/2)",
  :charLabels=>["(1,1)",  "(1,\\varepsilon)", "(g_2,1)","(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[[1,1,-1,-1],[1,1,1,1],[-1,1,1,-1],[-1,1,-1,1]],
  :eigenvalues=>[1,1,E(4),-E(4)],
  :qEigen=>[0,0,1//2,1//2],
  :perm=>Perm(3,4),
  :lusztig=>true, # does not satisfy (ST)^3=1 but (SPT)^3=1
  :cospecial=>2)
  )

function SubFamily(f,ind,scal,label)
  ind=filter(i->ind(f,i),1:length(f[:eigenvalues]))
  res=Dict(:fourierMat=>getindex.(f[:fourierMat][ind],Ref(ind))*scal,
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
  res
end

function SubFamilyij(f,i,j,scal)
  g=SubFamily(f,(f,k)->sum(f[:charSymbols][k])%j==i,scal,join([i,j]))
  g[:explanation]="subfamily(sum(charsymbols)mod $j=$i of $(f[:explanation]))"
  g
end

CHEVIE[:families][:ExtPowCyclic]=function(e,n)
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
  g
end

CHEVIE[:families][:X5]=SubFamilyij(CHEVIE[:families][:X](6),1,3,1-E(3))
CHEVIE[:families][:X5][:cospecial]=5
CHEVIE[:families][:Z4]=CHEVIE[:families][:ExtPowCyclic](4,1)

struct Family
  d::Dict{Symbol,Any}
end

function Family(s::String,v::Vector,d::Dict=Dict{Symbol,Any}())
  f=CHEVIE[:families][Symbol(s)]
  f[:charNumbers]=v
  Family(merge(f,d))
end
function Family(s::String,d::Dict=Dict{Symbol,Any}())
  f=CHEVIE[:families][Symbol(s)]
  Family(merge(f,d))
end
function Family(f::Dict{Symbol,Any},v::Vector,d::Dict=Dict{Symbol,Any}())
  f[:charNumbers]=v
  Family(merge(f,d))
end

function Base.show(io::IO,f::Family)
 print(io,"Family($(f.d[:name]):$(length(f.d[:eigenvalues])))")
end
# The big family in dihedral groups. For e=5 occurs in H3, H4
CHEVIE[:families][:Dihedral]=function(e)
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
  if iszero(e%2) then
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
  f
end

