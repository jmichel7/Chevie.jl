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
  charparams=[Any[1,0]]
  if iseven(m)
    res[:extRefl]=[1,5,4]
    m1=div(m,2)
    push!(charparams,[1,m1,"'"],[1,m1,"''"])
  else
    res[:extRefl]=[1,3,2]
  end
  push!(charparams,[1,m])
  append!(charparams,map(i->[2,i],1:div(m-1,2)))
  res[:charparams]=charparams
  res[:b]=map(x->x[2], charparams)
  res[:B]=map(phi->phi[1]==1 ? phi[2] : m-phi[2], charparams)
  charSymbols=map(1:div(m-1,2))do l
    S=map(i->[0],1:m)
    S[1]=[1]
    S[l+1]=[1]
    S
  end
  v=map(x->[0],1:m);v[m]=[1,2];pushfirst!(charSymbols,v)
  if iseven(m)
    v=map(x->[0],1:m);v[m]=[1];v[m1]=[1];pushfirst!(charSymbols,v)
    pushfirst!(charSymbols,copy(v))
  end
  v=map(x->[0,1],1:m);v[m]=[2];pushfirst!(charSymbols,v)
  res[:charSymbols]=charSymbols
  res[:malleParams]=map(x->map(partβ,x),charSymbols)
  if iseven(m)
    res[:malleParams]=convert.(Vector{Any},res[:malleParams])
    res[:malleParams][2]=push!(res[:malleParams][2][1:m1],1)
    res[:malleParams][3]=push!(res[:malleParams][3][1:m1],-1)
  end
  return res
end)
