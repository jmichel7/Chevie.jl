function combinations_sorted(mset::AbstractVector,k)
  if iszero(k) return [eltype(mset)[]] end
  res=Vector{eltype(mset)}[]
  for (i,e) in enumerate(mset)
    append!(res,map(x->vcat([e],x),combinations_sorted(mset[i+1:end],k-1)))
  end
  res
end 

combinations(mset,k)=combinations_sorted(sort(mset),k)
combinations(mset)=isempty(mset) ? [Int[]] : union(combinations.(Ref(mset),0:length(mset)))

function ArrangementsK(mset,blist,k)
  if iszero(k) return [eltype(mset)[]] end
  combs=Vector{eltype(mset)}[]
  for i in eachindex(mset)
    if blist[i] && (i==length(mset) || mset[i+1]!=mset[i] || !blist[i+1])
      blist1=copy(blist)
      blist1[i]=false
      append!(combs,pushfirst!.(ArrangementsK(mset,blist1,k-1),Ref(mset[i])))
    end
  end
  combs
end 

arrangements(mset,k)=ArrangementsK(sort(mset),fill(true,length(mset)),k)

# partitions of n of first (greatest) part <=m
function partitions_less(n,m)
  if m==1 return [fill(1,n)] end
  if iszero(n) return [Int[]] end
  res=Vector{Int}[]
  for i in 1:min(m,n)
    append!(res,map(x->vcat([i],x),partitions_less(n-i,i)))
  end
  res
end

partitions(n)=partitions_less(n,n)

if false
function partition_tuples(n,r)::Vector{Vector{Vector{Int}}}
  if r==1 return iszero(n) ? [[Int[]]] : map(x->[x],partitions(n)) end
  res=Vector{Vector{Int}}[]
  for i in  n:-1:1
    for p1 in partitions(i), p2 in partition_tuples(n-i,r-1)
      push!(res,vcat([p1],p2))
    end 
  end
  for p2 in partition_tuples(n,r-1)
    push!(res,vcat([Int[]],p2))
  end 
  res
end
else # bad implementation but which is ordered as GAP3
function partition_tuples(n, r)
   if n==0 return [fill(Int[],r)] end
   empty=(tup=[Int[] for i in 1:r], pos=fill(1,n-1))
   pm=[typeof(empty)[] for i in 1:n-1]
   for m in 1:div(n,2)
      for i in 1:r
         t1=map(copy,empty)
         t1.tup[i]=[m]
         t1.pos[m]=i
         push!(pm[m],t1)
      end
      for k in m+1:n-m
         for t in pm[k-m]
            for i in t.pos[m]:r
               t1=map(copy,t)
               t1.tup[i]=vcat([m],t1.tup[i])
               t1.pos[m]= i
               push!(pm[k], t1)
            end
         end
      end
   end
   res= Vector{Vector{Int}}[]
   for k in 1:n-1
      for t in pm[n-k]
         for i in t.pos[k]:r
            s=copy(t.tup)
            s[i]=vcat([k],s[i])
            push!(res,s)
         end
      end
   end
   for i in 1:r
      s=copy(empty.tup)
      s[i]=[n]
      push!(res,s)
   end
   res
end
end

