# Hurwitz action of the braid b on the list l of group elements
function hurwitz(b,l)
  l=copy(l)
  for s in word(b)
    if s>0 l[s:s+1]=[l[s+1],l[s]^l[s+1]]
    else l[-s:-s+1]=[l[-s]*l[-s+1]*inv(l[-s]),l[-s]]
    end
  end
  return l
end

# Implements the action of B_n on F_n
#   Input:  element of B_n, ambient free group F_m
#           (with m >= n; when m>n, Hurwitz
#            action is extended trivially to the extra generators)
#   Output: automorphism of F_n
BnActsOnFn(b,F)=GroupHomomorphismByImages(F,F,gens(F),hurwitz(b,gens(F)))

"""
'VKQuotient(braids)'

The input `braid` is a list of braids `b₁,…,b_d`, living in the braid group
on `n` strings. Each `bᵢ` defines by Hurwitz action an automorphism `φᵢ` of
the  free group `Fₙ`. The function return the group defined by the abstract
presentation: ``< f₁,…,fₙ ∣ ∀ i,j φᵢ(fⱼ)=fⱼ > ``

|    gap> B:=Braid(CoxeterGroupSymmetricGroup(3));
    function ( arg ) ... end
    gap> b1:=B(1)^3; b2:=B(2);                   
    1.1.1
    2
    gap> g:=VKQuotient([b1,b2]);                 
    Group( f.1, f.2, f.3 )
    gap>  last.relators;  
    [ f.2^-1*f.1^-1*f.2*f.1*f.2*f.1^-1, IdWord,
      f.2^-1*f.1^-1*f.2^-1*f.1*f.2*f.1, f.3*f.2^-1, IdWord, f.3^-1*f.2 ]
    gap> p:=PresentationFpGroup(g);DisplayPresentation(p);
    << presentation with 3 gens and 4 rels of total length 16 >>
    1: c=b
    2: b=c
    3: bab=aba
    4: aba=bab
    gap> SimplifyPresentation(p);DisplayPresentation(p);
    #I  there are 2 generators and 1 relator of total length 6
    1: bab=aba|
"""
function VKQuotient(braids)
  # get the true monodromy braids and the Hurwitz action basic data
  n=braids[1].M.W.n
  F=FpGroup(Symbol.('a'.+(0:n-1))...)
  rels=AbsWord[]
  f=gens(F)
  for b in braids append!(rels,map((a,b)->a*inv(b),hurwitz(b,f),f)) end
  F/rels
end

# A variant of the previous function.
# See arXiv:math.GR/0301327 for more mathematical details.
# Input: global VKCURVE record
# Output: the quotient, encoded as an FpGroup
############################################################
# Printing controlled by VKCURVE.showAction
function DBVKQuotient(r)
  # get the true monodromy braids and the Hurwitz action basic data
  n=r.braids[1].M.W.n
  F=FpGroup(Symbol.('a'.+(0:n+length(r.verticallines)-1))...)
# above the basepoint for the loops, locate the position of the string
# corresponding to the trivializing horizontal line
  zero=r.zeros[r.basepoint]
  dist=abs2.(zero-r.height)
  height=zero[argmin(dist)]
  basestring=count(z->(real(z),imag(z))<(real(height),imag(height)),zero)
  @show basestring,r.braids
  fbase=gens(F)[basestring]
  rels=[]
  auts=map(b->BnActsOnFn(b, F),r.braids)
  for aut in auts
# Find an element conjugator such that aut(fbase)^inv(conjugator)=fbase
    ifbase=Image(aut, fbase)
    conjugator=F.gens[1]/F.gens[1]
    choices=vcat(map(f->[f,f^-1], gens(f)[1:n])...)
    while true
      k=0
      while true
        k+=1
        if LengthWord(choices[k]*ifbase) < LengthWord(ifbase) break end
      end
      ifbase=(choices[k]*ifbase) // choices[k]
      conjugator=choices[k]*conjugator
      if LengthWord(ifbase)==1 break end
    end
# Replacing aut by  correctaut:= Conj(conjugator)*aut
    conj=GroupHomomorphismByImages(F, F, gens(F), gens(F).^inv(conjugator))
    conj[:isMapping]=true
    correctaut=CompositionMapping(conj, aut)
    if Position(auts, aut) > length(r.verticallines)
      rels=Append(rels, map(f->Image(correctaut, f)/ f, gens(F)[1:n]))
    else
      g=gens(F)[Position(auts, aut)+n]
      append!(rels, map(f->(Image(correctaut, f)*g)/ (g*f), gens(F)[1:n]))
    end
  end
  push!(rels, fbase)
  F/rels
end
