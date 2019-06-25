FamilyOps=Dict()

struct Family
  prop::Dict{Symbol,Any}
end

Family(f::Family)=f
function getf(s::String)
  f=CHEVIE[:families][Symbol(s)]
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

function Base.show(io::IO,f::Family)
 if !haskey(f,:name) f[:name]="???" end
 print(io,"Family($(TeXstrip(f.prop[:name])):$(length(f.prop[:eigenvalues])))")
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
#  if ForAll(arg,f->IsBound(f.special)) then 
#    res.special:=PositionCartesian(List(arg,Size),List(arg,f->f.special));
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
#  if ForAny(arg,f->IsBound(f.qEigen)) then
#    res.qEigen:=List(Cartesian(List(arg,function(f) 
#      if IsBound(f.qEigen) then return f.qEigen;else return f.eigenvalues*0;fi;
#      end)),Sum);
#  fi;
  res
end

CHEVIE[:families]=Dict(:C1=>
        Family(Dict(:group=>"C1", :name=>"C_1", :explanation=>"trivial",
                    :charLabels=>[""], :fourierMat=>hcat([1]), :eigenvalues=>[1],
  :mellin=>[[1]],:mellinLabels=>[""])),
  :C2=>Family(Dict(:group=>"C2", :name=>"C_2",
  :explanation=>"DrinfeldDouble(Z/2)",
  :charLabels=>["(1,1)", "(g_2,1)", "(1,\\varepsilon)", "(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[1 1 1 1;1 1 -1 -1;1 -1 1 -1;1 -1 -1 1],
  :eigenvalues=>[1,1,1,-1],
  :perm=>(),
  :mellin=>[[1,1,0,0],[1,-1,0,0],[0,0,1,1],[0,0,1,-1]],
  :mellinLabels=>["(1,1)","(1,g2)","(g2,1)","(g2,g2)"])),
  :S3=>Family(Dict(:group=>"S3", :name=>"D(S_3)",
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
                  "(g3,g3)","(g3,g3^2)"])),
  :X=>function(p)
    ss=combinations(0:p-1,2)
    Family(Dict(:name=>"R_{\\BZ/$p}^{\\wedge 2}",
         :explanation=>"DoubleTaft($p)",
         :charSymbols=>ss,
         :charLabels=>map(s->repr(E(p)^s[1],context=:TeX=>true)*
      "\\!\\wedge\\!"*repr(E(p)^s[2],context=:TeX=>true),ss),
    :eigenvalues=>map(s->E(p)^Product(s),ss),
    :fourierMat=>[(E(p)^(i*reverse(j))-E(p)^(i*j))/p for i in ss,j in ss],
    :special=>1,:cospecial=>p-1))
   end,
   Symbol("C'\"2")=>Family(Dict(:group=>"C2", :name=>"C'''_2",
  :explanation=>"TwistedDrinfeldDouble(Z/2)'",
  :charLabels=>["(1,1)", "(1,\\varepsilon)", "(g_2,1)", "(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[1 1 -1 -1;1 1 1 1;-1 1 -1 1;-1 1 1 -1],
  :eigenvalues=>[1,1,E(4),-E(4)],
  :qEigen=>[0,0,1,1]//2,
  :perm=>Perm(3,4),
  :cospecial=>2)),
   Symbol("C'2")=>Family(Dict(:group=>"C2",:name=>"C'_2",
  :explanation=>"TwistedDrinfeldDouble(Z/2)",
  :charLabels=>["(1,1)",  "(1,\\varepsilon)", "(g_2,1)","(g_2,\\varepsilon)"],
  :fourierMat=>1//2*[1 1 -1 -1;1 1 1 1;-1 1 1 -1;-1 1 -1 1],
  :eigenvalues=>[1,1,E(4),-E(4)],
  :qEigen=>[0,0,1,1]//2,
  :perm=>Perm(3,4),
  :lusztig=>true, # does not satisfy (ST)^3=1 but (SPT)^3=1
  :cospecial=>2))
  )

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
  Family(g)
end

CHEVIE[:families][:X5]=SubFamilyij(CHEVIE[:families][:X](6),1,3,1-E(3))
CHEVIE[:families][:X5][:cospecial]=5
CHEVIE[:families][:Z4]=CHEVIE[:families][:ExtPowCyclic](4,1)
CHEVIE[:families][:Z4][:fourierMat]*=-E(4)
CHEVIE[:families][:Z4][:eigenvalues]/=CHEVIE[:families][:Z4][:eigenvalues][2]
CHEVIE[:families][:Z4][:special]=2
CHEVIE[:families][:Z4][:qEigen]=[1,0,1,0]//2

CHEVIE[:families][:Z9]=CHEVIE[:families][:ExtPowCyclic](9,1)
#if CHEVIE.families.Z9.eigenvalues<>List([0..8],i->E(9)^(5*i^2))then Error();fi;
CHEVIE[:families][:Z9][:perm]=perm"(2,9)(3,8)(4,7)(5,6)"
CHEVIE[:families][:Z9][:qEigen]=[0,2/3,1/3,0,2/3,1/3,0,2/3,1/3]

CHEVIE[:families][:QZ]=function(n)
  pairs=[(i,j) for i in 0:n-1 for j in 0:n-1]
  res=Dict{Symbol,Any}(:name=>"D(\\BZ/$n)")
  res[:explanation]="Drinfeld double "*res[:name]
  res[:fourierMat]=[E(n,x)^c1*E(n,x1)^c for (x,c) in pairs, (x1,c1) in pairs]
  res[:eigenvalues]=[E(n,x)^c for (x,c) in pairs]
  res[:special]=1
  res[:charLabels]=["($(E(n,x)),$(E(n,c)))" for (x,c) in pairs]
  Family(res)
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
end

function Drinfeld_double(g;lu=false)
  g=arg[1]
  res=Dict{Symbol,Any}(:group=> g)
  res[:classinfo] = map(function (c, n) local r, t
    r = Dict{Symbol, Any}(:elt => Representative(c), :name => n)
    if r[:elt] == g[:identity] r[:name] = "1" end
    r[:centralizer] = Centralizer(g, r[:elt])
    r[:centelms] = map(Representative, ConjugacyClasses(r[:centralizer]))
    t = CharTable(r[:centralizer])
    r[:charNames] = CharNames(r[:centralizer], Dict{Symbol, Any}(:TeX => true))
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

FamilyImprimitive = function (S)
local e, Scoll, ct, d, m, ll, eps, equiv, nrSymbols, epsreps, trace, roots, i, j, mat, frobs, symbs, newsigns, schon, orb, mult, res, IsReducedSymbol
  println("S=$S")
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
    if (mat*DiagonalMat(frobs))^3!=mat^0
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
          v = map((z->begin count(y->y>=z, M) % 2 end), Z1)
          D = length(v)
          v1 = v[2:2:D-D%2]
          v2 = v[3:2:(D-1)+D%2]
          if D%2==1 push!(v1,0) end
          v1 = map(i->sum(v1[[i,i+1]]) % 2, 1:length(v2))
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
