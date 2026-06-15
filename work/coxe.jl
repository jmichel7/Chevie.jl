Perms.mul!(a,b)=a*=b

function Garside.leftgcd(W::CoxeterGroup{T},elts::Vararg{T,N})where{T,N}
  if N==0 error("leftgcd needs at least one element as argument") end
  res::T=copy(one(W))
  while true
    i=firstleftdescent(W,elts...)
    if isnothing(i) return res end
    s=W(i)
    res*=s
    elts=map(x->s*x,elts)
  end
end

Garside.rightgcd(W,a,b)=affinv(W,leftgcd(W,affinv(W,a),affinv(W,b))) 

function divisors(W::CoxeterGroup{T},s::T)where T
  rest=[(left=one(W)::T,right=s)]
  res=typeof(rest)[]
  while !isempty(rest)
    push!(res,rest)
    new=empty(rest)
    for x in rest, i in leftdescents(W,x.right)
      s::T=W(i)
      push!(new,(left=x.left*s,right=s*x.right))
    end
    rest=unique!(new)
  end
  res
end

function left_divisorsPoset(W::CoxeterGroup{T},s::T)where T
  rest=[(left=one(W)::T,right=s)]
  hasse=[Int[]]
  res=typeof(rest)[]
  l=0
  while !isempty(rest)
    push!(res,rest)
    new=empty(rest)
    for (nx,x) in enumerate(rest), i in leftdescents(W,x.right)
      s::T=W(i)
      c=(left=x.left*s,right=s*x.right)
      p=findfirst(==(c),new)
      if isnothing(p) push!(new,c); push!(hasse,Int[]); p=length(new) end
      push!(hasse[l+nx],p+l+length(rest))
    end
    l+=length(rest)
    rest=new
  end
  p=Poset(CPoset(hasse),vcat(map(x->first.(x),res)...))
  p.show_element=(io,x,n)->(e=x.elements[n];print(io,isone(e) ? "." :
                                                 joindigits(word(W,e))))
  p
end

function left_divisorsPoset(s::Garside.LocallyGarsideElt)
  rest=[(left=one(s),right=s)]
  hasse=[Int[]]
  res=typeof(rest)[]
  l=0
  while !isempty(rest)
    push!(res,rest)
    new=empty(rest)
    for (nx,x) in enumerate(rest), i in leftdescents(x.right)
      t=s.M.atoms[i]
      c=(left=x.left*t,right=t\x.right)
      p=findfirst(==(c),new)
      if isnothing(p) push!(new,c); push!(hasse,Int[]); p=length(new) end
      push!(hasse[l+nx],p+l+length(rest))
    end
    l+=length(rest)
    rest=new
  end
  Poset(CPoset(hasse),vcat(map(x->first.(x),res)...))
end

function right_divisorsPoset(s::Garside.LocallyGarsideElt)
  p=left_divisorsPoset(reverse(s))
  p.elements.=reverse.(p.elements)
  p
end

function right_divisorsPoset(W::CoxeterGroup{T},s::T)where T
  p=left_divisorsPoset(W,affinv(W,s))
  p.elements.=affinv.(Ref(W),p.elements)
  p
end

function Garside.left_divisors(W::CoxeterGroup{T},s::T)where T
  map(x->first.(x),divisors(W,s))
end

function Garside.right_divisors(W::CoxeterGroup{T},s::T)where T
  map(x->last.(x),divisors(W,s))
end

affinv(W,w)=W(reverse(word(W,w))...)

function cwords(n::Int,a::T,b::T,v::T)where T
  res=[[v*a,v*b]]
  for i in 2:n-1
    x,y=res[end]
    push!(res,iseven(i) ? [x*b,y*a] : [x*a,y*b])
  end
  x=res[end][1]
  push!(res,iseven(n) ? [x*b] : [x*a])
  res
end

function Garside.rightlcm(W::CoxeterGroup{T},a::T,b::T;lim=length(W,a)+length(W,b)+6)where T
  if isfinite(W)
    w0=longest(W)
    return w0*leftgcd(W,w0/a,w0/b)
  end
  C=coxmat(W)
  d=leftgcd(W,a,b)
  a=affinv(W,d)*a;b=affinv(W,d)*b
  da=left_divisors(W,a)
  db=left_divisors(W,b)
  m=min(length(W,a),length(W,b))
  e=union.(da[1:m+1],db[1:m+1])
  append!(e,da[m+2:end])
  append!(e,db[m+2:end])
  i=2
  while i<=length(e) && i<=lim
#  print(length(e[i]), " ")
    for w in e[i-1]
      l=let i=i
        filter(j->w*W(j) in e[i],1:ngens(W))
      end
      for p in combinations(l,2)
        n=C[p[1],p[2]]
        if iszero(n) return end
        cw=cwords(n,W(p[1]),W(p[2]),w)
        for j in 2:length(cw), v in cw[j] 
          if i+j-1>length(e) push!(e,[v])
          elseif !(v in e[i+j-1]) push!(e[i+j-1],v)
          end
        end
      end
    end
    i+=1
  end
#  println()
  if i>length(e) # xprintln("lcm=",B(only(e[end])))
    return d*only(e[end])
#   return only(e[end])
  end
end

function Garside.leftlcm(W::CoxeterGroup{T},a::T,b::T;lim=length(W,a)+length(W,b)+6)where T
   w=rightlcm(W,affinv(W,a),affinv(W,b))
   if isnothing(w) return
   else return affinv(W,w)
   end
end

function garside_closure(W::CoxeterGroup{T},S::Vector{T}=gens(W);verbose=true)where T
  res=vcat([one(W)::T],S)
  i=2
  while i<length(res)
    j=i+1
    while j<=length(res)
#     print(".")
      r=rightlcm(W,res[i],res[j])
      if isnothing(r)|| r in res u=T[]
      else u=vcat(right_divisors(W,r)...)
      end
      u=filter(x->!(x in res),u)
      append!(res,u)
      if verbose && !isempty(u)
        ww(w)=isone(w) ? "." : joindigits(word(W,w))
        println(length(res),": i=",i," j=",j," lcm(",ww(res[i]),",",ww(res[j]),
                 ")=>",join(ww.(u),","))
      end
      j+=1
    end
    i+=1
  end
  sort(res,by=x->length(W,x))
end  

rightdivide(W,a,b)=length(W,a*affinv(W,b))+length(W,b)==length(W,a)
leftdivide(W,a,b)=length(W,affinv(W,a)*b)+length(W,a)==length(W,b)

# elements maximaux pour la divisibilite à droite (dans une liste dont l'ordre
# étend la divisibilité)
function maxright(W,l)
  filter(i->!any(j->rightdivide(W,l[j],l[i]),i+1:length(l)),eachindex(l))
end

rightdescents(W,w)=leftdescents(W,affinv(W,w))

function Base.union(p::Poset...)
  elements=union(map(x->x.elements,p)...)
  I=Matrix{Bool}(falses(length(elements),length(elements)))
  for P in p
    l=indexin(P.elements,elements)
    I[l,l].|=incidence(P)
  end
  I=transitive_closure(I)
  Poset(CPoset(I),elements)
end

function rightdivisibilityPoset(W,l)
  p=Poset((y,x)->length(W,x*inv(y))+length(W,y)==length(W,x),l)
  p.show_element=(io,x,n)->(e=x.elements[n];print(io,isone(e) ? "." :
                                                 joindigits(word(W,e))))
  p
end

function islow(W,w)
  ww=word(W,w)
  for i in eachindex(ww)
    root=permutedims(reflrep(W,W(ww[i-1:-1:1]...))[ww[i],:])
    if !(root in smallroots(W)) &&
      length(W,W(deleteat!(copy(ww),i)...))==length(ww)-1
      return false
    end
  end
  true
end  
  
function low(W)
  i=0
  res=eltype(W)[]
  old=[W()]
  while true
    append!(res,old)
    new=empty(old)
    for w in old
      for j in 1:ngens(W)
        iw=W(j)*w
        if !isleftdescent(W,w,j) && islow(W,iw)
          push!(new,iw)
        end
      end
    end
    if isempty(new) return res end
    old=unique(new)
    println("length=$i low=",length(old))
    i+=1
  end
end  

function smallroots(W)
  get!(W,:smallroots)do
  T=eltype(cartan(W))
  sr=map(i->one(T)*permutedims(map(j->j==i,1:ngens(W))),1:ngens(W))
  while true
    new=false
    for (i,s) in enumerate(gens(W))
      for alpha in sr
        alpha1=alpha*reflrep(W,i)
        if 0< alpha1[i]-alpha[i] <2 && !(alpha1 in sr)
       	  new=true
 	  push!(sr,alpha1)
        end
      end 
    end
    if !new return sr end
  end  
  end
end
