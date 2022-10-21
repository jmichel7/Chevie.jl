##########################################################################
# brbase.jl                                Meinolf Geck and Sungsoon Kim
#
# Copyright (C) 1996     Equipe des groupes finis, Universite Paris VII
#
#
# Example:  If W=coxgroup(:type,rank) is a Coxeter group then the command
# `BaseBruhat(W)'  adds a component `incidence' to `W' which contains for
# each element in `W.elts' the associated boolean vector defined above.
# If one is only interested in the set of bi-grassmannians and the set of
# base  elements  (especially  for  large  groups),  one  should  use the
# commands:
#             gap> W := CoxeterGroup( \"F\", 4 );       # for example
#             gap> bg := bi_grassmannians( W );
#             gap> base := FindBaseBruhat( W, bg );
#
# written MG, SK, Oct 1996
############################################################################
"""
`bi_grassmannians(W)`

returns the set of bi-grassmannians in the Coxeter group `W`, that is those
elements whose right and left descent sets are both of length 1. The result
is  a  matrix  `bi`  such  that  `bi[i,j]`  contains  those  elements  with
leftdescents [i] and rightdescents [j].
"""
function bi_grassmannians(W)
  toM(map(1:ngens(W)) do k
    bi=[empty(gens(W)) for i in 1:ngens(W)]
    ea=[one(W)]
    while true
      en=empty(gens(W))
      for a in ea
        for r in 1:ngens(W)
          if !isleftdescent(W,inv(a),r)
            x=a*W(r)
            if !(x in en) && !any(j->j!=k && isleftdescent(W,x,j),1:ngens(W))
              push!(en,x)
            end
          end
        end
      end
      for i in en
        l=leftdescents(W,inv(i))
        if length(l)==1 push!(bi[only(l)],i) end
      end
      ea=en
    if isempty(en) break end
    end
#   InfoChevie("#I ",length.(bi),"\n")
    bi
  end)
end

"""
`BaseBruhat(W)`

The base of a poset relation is defined to be the set of all elements w ∈ W
which  cannot be  obtained as  the supremum  of a  subset not containing w.
Denote the base by B. Then every element w in W can be coded by the boolean
vector (e_b)_{b in B} recording the relation b<=w. The Bruhat order on W is
given   by  the   subset  relation   on  such   vectors.  A   `Dict`  named
`W.bruhatincidence`  records such vectors. The key  for `w∈ W` is the image
of `inclusiongens(W)` by `w` and the value is the vector `e_b`.

The  bases for all finite Coxeter  groups are determined in the article
'Bases for the Bruhat--Chevalley order on all finite Coxeter groups'.
"""
function BaseBruhat(W)
  get!(W,:bruhatincidence)do
    bg=bi_grassmannians(W)
    mins=[] 
    InfoChevie("#I ")
    for r in axes(bg,1)
      for s in axes(bg,2)
        InfoChevie("\t",length(bg[r,s]))
        l1=bg[r,s]
        if !isempty(l1)
          inz=fill(0,length(l1),length(l1))
          for i in eachindex(l1)
            for j in 1:i if bruhatless(W,l1[j],l1[i]) inz[i,j]=1 end end
          end
          for z in eachindex(l1)
            leqz=filter(j->inz[z,j]==1,1:z)
            leqz=setdiff(1:length(l1),leqz)
            for i in leqz
              z1=filter(<=(i),leqz)
              if sum(inz[i,z1])==1 && !(bg[r,s][i] in mins)
                push!(mins,bg[r,s][i]);
              end
            end
          end
          if all(==(1),inz[:,1]) push!(mins,bg[r,s][1]) end
        end
      end
      InfoChevie("\n#I ")
    end
    InfoChevie("No of base elements = ",length(mins),"\n")
    p=sortPerm(map(i->length(W,i),mins))
    base=permute(mins,inv(p))
    InfoChevie("#I Calculating incidence matrix...")
    incidence=Dict{Vector{Int},BitVector}()
    for w in elements(W)
      incidence[inclusiongens(W).^w]=BitVector(bruhatless.(Ref(W),base,w))
    end
    InfoChevie("\n")
    incidence
  end
end

function CoxGroups.bruhatless(W,x,y)
  if x==one(W) return true end
  d=length(W,y)-length(W,x)
  if haskey(W,:bruhatincidence)
    println("hello")
    if d<0 return false end
    iy=W.bruhatincidence[inclusiongens(W).^y]
    ix=W.bruhatincidence[inclusiongens(W).^x]
    return iy==iy.|ix
  end
  while d>0
    i=firstleftdescent(W,y)
    s=W(i)
    if isleftdescent(W,x,i)
      if x==s return true end
      x=s*x
    else d-=1
    end
    y=s*y
  end
  return x==y
end
