# Hand-translated part of chevie/tbl/cmpxtimp.g
# (C) jean Michel 2011-2017
#
# Data on spetses 3G333, 4G333 and 3G422

onsets(s,g)=sort(s.^g)

chevieset(:timp, :CharName, function (p,q,r,phi,para,opt)
           CHEVIE[:imp][:CharName](p,q,r,convert.(Vector{Int},para),opt)
    end)

chevieset(:timp, :ReducedInRightCoset, function (W, phi)
  sets = [[1, 2, 3, 44], [21, 3, 1, 32], [3, 11, 2, 36], [22, 3, 2, 16]]
  sets2 = [[1, 50, 3, 12, 2], [3, 52, 2, 23, 11], [1, 16, 3, 43, 38], 
           [2, 37, 3, 15, 14], [50, 3, 52, 38, 53], [1, 23, 3, 22, 45]]
  for i in sets
    y=gapSet(map(j->reflection(W,j),i))
    e=transporting_elt(W, y, onsets(y, phi);action=onsets)
    if !isnothing(e)return Dict{Symbol,Any}(:gen=>i[1:3],:phi=>phi/e) end
  end
  for i in sets2
    y=inclusion(W,i[1:4])
    v=gapSet(y)
    e=transporting_elt(W,v,onsets(v,phi);action=onsets)
    if !isnothing(e) && Perm(y,y.^(phi/e))==Perm(1,2,3,4)
      return Dict{Symbol, Any}(:gen => i[[1, 5, 3]], :phi=>phi/e)
    end
  end
  return false
end)
