"""
Posets are represented as Dicts where at least one of the two
following fields is present:

  `:incidence`:  a  boolean  matrix  such  that `:incidence[i][j]==true` iff
  `i<=j` in the poset.

  `:hasse`: a list representing the Hasse diagram of the poset: the i-th
  entry  is the list of indices  of elements which are immediate successors
  (covers)  of the i-th element, that is  the list of j such that `i<j`
  and such that there is no k such that i<k<j.

If  only one field is present, the other  is computed on demand. Here is an
example of use;

```julia-repl
julia> p=Poset(coxgroup(:A,2))
Poset with 6 elements

julia> hasse(p)
6-element Array{Array{Int64,1},1}:
 [2, 3]
 [4, 5]
 [4, 5]
 [6]   
 [6]   
 []    

julia> incidence(p)
6×6 Array{Bool,2}:
  true   true   true   true   true  true
 false   true  false   true   true  true
 false  false   true   true   true  true
 false  false  false   true  false  true
 false  false  false  false   true  true
 false  false  false  false  false  true
```
"""
module Posets
using Gapjm
export lcm_partitions, gcd_partitions, Poset, linear_extension, hasse,
 incidence, partition

function lcm_partitions(arg...)
  function lcm2(a,b)
    res = Set(Vector{Int}[])
    for p in b
     push!(res, sort(union(filter(x->!isempty(intersect(x, p)),a)...)))
    end
    b = Set(Vector{Int}[])
    for p in res
     push!(b, sort(union(filter(x->!isempty(intersect(x, p)),res)...)))
    end
    b
  end
  reduce(lcm2,arg)
end

function gcd_partitions(arg...)
  function gcd2(a,b)
    res = map(x->map(y->intersect(x, y), b), a)
    Set(filter(x->!isempty(x),vcat(res...)))
  end
  reduce(gcd2,arg)
end

struct Poset
  prop::Dict{Symbol,Any}
end

Poset(m::Matrix{Bool})=Poset(Dict(:incidence=>m,:size=>size(m,1)))
Poset(m::Vector{Vector{Bool}})=Poset(toM(m))
Poset(m::Vector{<:Vector{<:Integer}})=Poset(Dict(:hasse=>m,:size=>length(m)))
Base.length(p::Poset)=p.prop[:size]

Base.show(io::IO,p::Poset)=print(io,"Poset with ",length(p)," elements")

function linear_extension(P::Poset)
  ord=hasse(P)
  n=zeros(length(ord))
  for v in ord for x in v n[x] += 1 end end
  Q=filter(x->iszero(n[x]),1:length(n))
  res=Int[]
  while length(Q) > 0
    push!(res, Q[1])
    for x in ord[Q[1]]
      n[x]-=1
      if iszero(n[x]) push!(Q, x) end
    end
    popfirst!(Q)
  end
  if sum(n)>0 error("cycle") end
  res
end

function hasse(p::Poset)
  gets(p,:hasse) do p
    m=Int.(incidence(p))
    map(x->filter(y->x[y]==2,1:length(x)),eachrow(m*m))
  end
end

function incidence(p::Poset)::Matrix{Bool}
  gets(p, :incidence) do p
    n = linear_extension(p)
    incidence = one(Matrix{Bool}(undef,length(n),length(n)))
    for i in length(n)-1:-1:1
      for x in hasse(p)[n[i]] incidence[n[i],:].|= incidence[x,:] end
    end
    incidence
  end
end

function chains(P::Poset)
  ch = Vector{Int}[]
  h = hasse(P)
  for i in linear_extension(P)
    for j in h[i]
      p = findfirst(c->i==c[end],ch)
      if isnothing(p) push!(ch, [i, j])
      else push!(ch[p], j)
      end
    end
  end
  ch
end

function Base.reverse(p::Poset)
  res = deepcopy(p)
  if haskey(p.prop, :incidence)
    res.prop[:incidence]=permutedims(incidence(p))
  end
  if haskey(p.prop,:hasse)
    res.prop[:hasse] = map(empty,hasse(p))
    for i in 1:length(p)
      for j in hasse(p)[i] push!(hasse(res)[j], i) end
    end
  end
  return res
end

function partition(p::Poset)
  if haskey(p.prop, :hasse)
    l=reverse(p)
    res=groupby(i->[hasse(l)[i], hasse(p)[i]],1:length(p))
    sort(collect(values(res)),by=v->[hasse(l)[v[1]], hasse(p)[v[1]]])
  else
    I=incidence(p)
    ind=1:length(p)
    l=map(i->[map(j->j!=i && I[i][j], ind), map(j->j!=i && I[j][i], ind)], ind)
    map(x->filter(i->l[i] == x,ind), gapSet(l))
  end
end

function showgraph(x; opt...)
  p=partition(x)
  s=hasse(x)
  s=Poset(map(x->unique(sort(convert(Vector{Int},map(y->findfirst(z->y in z,p),
                                                     s[x[1]])))), p))
  labels=map(y->join(map(n->haskey(x.prop,:label) ? x.prop[:label](x,n,opt) :
                         string(n),y),","),p)
  opt=Dict(opt)
  if haskey(opt, :symbol) sep = opt[:symbol]
  elseif haskey(opt, :TeX) sep = "{<}"
  else sep = "<"
  end
  s=map(x->join(labels[x],sep), chains(s)) 
  if haskey(opt,:TeX) "\\noindent"*join(map(x->"\$$x\$\\hfill\\break\n",s))
  else join(s,"\n")
  end
end

function Gapjm.restricted(p::Poset,ind::Vector{<:Integer})
  res = Poset(copy(p.prop))
  if length(ind) == length(p) && sort(ind) == 1:length(p)
    if haskey(res.prop, :hasse)
     res.prop[:hasse] = Vector{Int}.(map(x->map(y->findfirst(isequal(y),ind),x),
       res.prop[:hasse][ind]))
    end
    if haskey(res.prop, :incidence)
      res.prop[:incidence] = incidence(res)[ind,ind]
    end
  else
    inc=incidence(p)
    res.prop[:incidence] = map(i->map(function(j)
                              if i != j && ind[i] == ind[j] return false
                              else return inc[ind[i],ind[j]]
                              end
                          end, 1:length(ind)), 1:length(ind))
    delete!(res.prop, :hasse)
  end
  res.prop[:size] = length(ind)
  if haskey(p.prop, :label)
    res.prop[:indices]=ind
    res.prop[:label]=(x, n, opt)->p[:label](x, x[:indices][n], opt)
  end
  res
end

function checkl(ord::Matrix{Bool})
  subl=Set(Vector{Bool}[])
  n=size(ord,1)
  for i in 1:n
    for j in 1:i-1
      if !ord[j,i] || ord[i,j]
        l=ord[i,:].&ord[j,:]
        if !l in subl
          if !any(y->l==l.&y,ord[l,:])
            for k in (1:n)[l]
              ll=copy(ord[k,:])
              ll[k]=false
              l.&=.!ll
            end
            InfoChevie("# $i  &&  $j have bounds", (1:n)[l], "\n")
            return false
          end
          push!(subl, l)
        end
      end
    end
  end
  true
end

IsJoinLattice(P::Poset)=checkl(incidence(P))
IsMeetLattice(P::Poset)=checkl(permutedims(incidence(P)))

# Bruhat interval [1,w]
function Poset(W::CoxeterGroup,w=longest(W))
  if w==one(W) return Poset(Dict(:elts=>[w],:hasse=>[Int[]],
    :action=>map(x->[0],gens(W)),:size=>1,
    :label=>(p,n,opt)->joindigits(word(p.prop[:W],p.prop[:elts][n])),
    :W=>W))
  end
  s=firstleftdescent(W,w)
  p=Poset(W,W(s)*w)
  l=length(p)
  new=filter(k->iszero(p.prop[:action][s][k]),1:l)
  append!(p.prop[:elts],W(s).*(p.prop[:elts][new]))
  append!(hasse(p),map(x->Int[],new))
  p.prop[:action][s][new]=l.+(1:length(new))
  for i in eachindex(gens(W)) 
    append!(p.prop[:action][i],i==s ? new : zeros(Int,length(new)))
  end
  for i in 1:length(new) push!(hasse(p)[new[i]],l+i) end
  for i in 1:l 
    j=p.prop[:action][s][i]
    if j>i
      for h in p.prop[:action][s][hasse(p)[i]]
        if h>l push!(hasse(p)[j],h)
          k=findfirst(isequal(p.prop[:elts][h]/p.prop[:elts][j]),gens(W))
          if !isnothing(k) p.prop[:action][k][j]=h;p.prop[:action][k][h]=j end
        end
      end
    end
  end
  p.prop[:size]=length(hasse(p))
  p
end

end
