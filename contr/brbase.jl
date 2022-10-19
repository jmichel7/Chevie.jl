#########################################################################
##
#A  brbase.jl                                Meinolf Geck and Sungsoon Kim
##
#Y  Copyright (C) 1996     Equipe des groupes finis, Universite Paris VII
##
##  This  file contains  functions for  computing bi-grassmannians  and the
##  base  of finite Coxeter groups.  
##
##  The  base is defined to be the set  of all elements w in W which cannot
##  be  obtained as the supremum  of a subset not  containing w. Denote the
##  base by B. Then every element w in W can be coded by the boolean vector
##  (e_b)_{b  in B} where e_b=1  if b <= w  and e_b=0 otherwise. The Bruhat
##  order  on W is given by the  operation of `boolean OR' on such vectors.
##  The  bases for all finite Coxeter  groups are determined in the article
##  'Bases for the Bruhat--Chevalley order on all finite Coxeter groups'.
## 
##  Example:  If W=coxgroup(:type,rank) is a Coxeter group then the command
##  `BaseBruhat(W)'  adds a component `incidence' to `W' which contains for
##  each element in `W.elts' the associated boolean vector defined above.
##  If one is only interested in the set of bi-grassmannians and the set of
##  base  elements  (especially  for  large  groups),  one  should  use the
##  commands:
##              gap> W := CoxeterGroup( \"F\", 4 );       # for example
##              gap> bg := BiGrassmannians( W );
##              gap> base := FindBaseBruhat( W, bg );
##
#H  written MG, SK, Oct 1996
##
#F BiGrassmannians( <W> ) . . . . . . . . . . . .  Bi-Grassmannians of W 
#F BaseBruhat( <W> ) . . . . . . . .  the base and the coding of elements
##
###########################################################################
"""
`BiGrassmannians(W)`

'BiGrassmannians' computes the set of bi-grassmannians in the Coxeter group
`W`,  that is those elements whose right  and left descent sets are both of
length  1. The result is  a matrix `bi` such  that `bi[i,j]` contains those
elements with leftdescents [i] and rightdescents [j].
"""
function BiGrassmannians(W)
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

`FindBaseBruhat` computes the base of the Bruhat order on the Coxeter group
`W`, and then adds a component `incidence` to `W` which contains the coding
of the elements of `W` by boolean vectors defined by the base.
"""
function BaseBruhat(W)
  get!(W,:incidence)do
    bg=BiGrassmannians(W)
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
      incidence[(1:ngens(W)).^w]=BitVector(map(v->bruhatless(W,v,w),base))
    end
    InfoChevie("\n")
    incidence
  end
end

function bruhat2(W,y,w)
  if length(W,y)>=length(W,w) return false end
  iw=W.incidence[(1:ngens(W)).^w]
  iy=W.incidence[(1:ngens(W)).^y]
  iw==iw.|iy
end
