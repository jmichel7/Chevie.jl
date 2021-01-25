# Hand-translated part of chevie/tbl/cmp4_22.g
# (C) 1998 - 2017  Gunter Malle & Jean Michel
#
# Data about the complex reflection groups with Shephard-Todd numbers 4--22 
#
# The data for classes, characters and Schur elements comes from
# G.Malle, ``Degres relatifs des algebres de Hecke cyclotomiques associees aux
#  groupes de reflexions complexes de dimension 2'',  in
# ``Finite  reductive groups'', Progress in  math. n0 141, Birkhauser 1997
# to which should be added the following info:
# -- the subalgebra H' is generated as described by the function "Embed"
# -- there is a misprint in the numerator of the relative degrees of
#    2-dimensional characters of G7: it should read
#        -x1^3 x_2^3 y_1 y_2^3 y_3^4 z_1^2 z_2^2 z_3^4
#    instead of
#        -x1^2 x_2^2 y_1 y_2^3 y_4   z_1^2 z_2^2 z_3^4
# -- There are  misprints in  the relative  degrees of  the 4-dimensional
#    characters of G11:
# -- The  numerator  should   read  -(x1x2)^9y1^10(y2y3)^5z1^6(z2z3z4)^4t^2
# -- In the last  product in the denominator i should run in {2,3}.
# -- there is a misprint in the relative degrees of the 3-dimensional
#    characters of G19:
#    The numerator should read x1^10x2^15y1^7z1^3z4^12

const G4_22IndexChars_dict=Dict{Int,Any}()

for i in 4:22 G4_22IndexChars_dict[i]=Dict() end

CHEVIE[:CheckIndexChars]=true

function G4_22FetchIndexChars(ST, para)
  if !CHEVIE[:CheckIndexChars]
    return chevieget(:G4_22, :CharInfo)(ST)[:indexchars]
  end
  get!(G4_22IndexChars_dict[ST],para)do
    chevieget(:G4_22, :HeckeCharTable)(ST, para, [])[:indexchars]
  end
end

# tests if res=Chartable(Hecke(G_ST,res[:parameter])) is correct
# where rows=irrs for generic group (G7, G11 or G19)
# and i is the selector from rows to test.
function G4_22Test(res,rows,i)
  ST=res[:ST]
  T(ST)=string("G",ST in 4:7 ? 7 : ST in 8:15 ? 11 : 19)
  if haskey(G4_22IndexChars_dict[ST],res[:parameter])
#   InfoChevie("IndexChars(Hecke(G_$ST,",
#              HeckeAlgebras.simplify_para(res[:parameter]),"))\n")
    ic=G4_22IndexChars_dict[ST][res[:parameter]]
    res[:irreducibles]=rows[ic]
    if ic!=i
      println("*** WARNING: choice of character restrictions from ", T(ST), 
        " for this specialization does not agree with group CharTable")
      if !CHEVIE[:CheckIndexChars]
        print("Try again with CHEVIE[:CheckIndexChars]=true\n")
      end
    end
    return ic
  end
  ic=i
  res[:irreducibles]=rows[ic]
  if length(Set(res[:irreducibles]))==length(res[:classes]) l=i
  else
    l=map(x->findfirst(==(x),rows), rows)
    if length(Set(l))!=length(res[:classes])
      error("specialization not semi-simple")
    end
    xprintln("*** WARNING: bad choice of char. restrictions from ",T(ST), 
             " for H(G$ST,",HeckeAlgebras.simplify_para(res[:parameter]),")")
    if !CHEVIE[:CheckIndexChars]
      print("Try again with CHEVIE[:CheckIndexChars]=true\n")
    end
 #  l=map(x->filter(i->l[i]==x,eachindex(l)), sort(unique(l)))
    l=map(x->findall(==(x),l), unique(sort(l)))
    o=filter(x->count(j->j in x,i)>1,l)
  # println(" over-represented by ", intersect(union(o...), i)," : ", o)
  # println(" absent : ",filter(x->iszero(count(j->j in x,i)),l))
    l=first.(l)
    o=filter(p->first(p)!=last(p),collect(zip(i,l)))
    println("changing choice ",join(first.(o),",")," → ",join(last.(o),","))
  # println(" Choosing ",l)
    res[:irreducibles]=rows[l]
  end
  G4_22IndexChars_dict[ST][res[:parameter]]=l
end
