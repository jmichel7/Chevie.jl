# Hand-translated part of chevie/tbl/coxi.jl
# (C) 1992-2010  Meinolf Geck & Jean Michel
# Data for W(I_m)
#
#  The  characters of I_2(m) are uniquely  parametrized by [d,b] where d is
#  their  degree and b  is the valuation  of the fake  degrees, except that
#  when  m is  even, there  are two  characters with [d,b]=[1,m/2]. The one
#  which  maps the generators to [1,-1]  is denoted [1,m/2,"'"] and the one
#  which maps them to [-1,1] is denoted [1,m/2,"''"].

chevieset(:I, :CharInfo, function(m)
  res=Dict{Symbol, Any}()
  charparams=[[1,0]]
  if iseven(m)
    res[:extRefl]=[1,5,4]
    m1=div(m,2)
    push!(charparams,[1,m1,1],[1,m1,2])
  else
    res[:extRefl]=[1,3,2]
  end
  push!(charparams,[1,m])
  append!(charparams,map(i->[2,i],1:div(m-1,2)))
  res[:charnames]=exceptioCharName.(charparams)
  res[:charparams]=charparams
  res[:b]=map(x->x[2], charparams)
  res[:B]=map(phi->phi[1]==1 ? phi[2] : m-phi[2], charparams)
  charSymbols=[] # need type Any for m even
  push!(charSymbols,map(i->i==m ? [2] : [0],1:m))
  if iseven(m)
    v=vcat(map(i->i==m1 ? [1] : [0],1:m1),[2,0])
    push!(charSymbols,v)
    v=copy(v);v[m1+2]=1;push!(charSymbols,v)
  end
  push!(charSymbols,map(i->i==m ? [1,2] : [0,1],1:m))
  append!(charSymbols,map(1:div(m-1,2))do l
    map(i->i==1 ? [1] : i==l+1 ? [1] : [0],1:m)
  end)
  res[:charSymbols]=charSymbols
  res[:A]=degree_gendeg_symbol.(charSymbols)
  res[:a]=valuation_gendeg_symbol.(charSymbols)
  # malleParams are the partitiontuples for the symbols
  res[:malleParams]=map(x->map(partβ,fullsymbol(x)),charSymbols)
  if iseven(m)
    res[:malleParams]=convert.(Vector{Any},res[:malleParams])
    res[:malleParams][2]=append!(res[:malleParams][2][1:m1],[2,0])
    res[:malleParams][3]=append!(res[:malleParams][3][1:m1],[2,1])
  end
  return res
end)

# parameters for cuspidal symbols
chevieset(:I, :SymbolToParameter, function (S)
  if S[1]!=[0,1] || !any(isempty,S) return false end
  if isodd(length(S)) S=reverse(S) end
  a=findfirst(isempty,S)
  if isodd(length(S)) 
    [a,findfirst(==([0,1]),S)-a]
  else (a-1).+[-findfirst(==([0,1]),S[2:end]),0]
  end
end)
