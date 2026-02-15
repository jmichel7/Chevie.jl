using Chevie
"""
 `HeConjugate(WF,i,e;w,canonical)`

 This function uses the He-Nie algorithms to find a representative of the
 `i`-th class of the Coxeter group or coset `WF`.

 `e` must be the set of eigenvalues of the i-th class with argument at most
 `π`,  represented as `Root1`.  The returned representative  depends on the
 order  in which `e` is given. If `e` is given in increasing argument order
 a  minimal  length  representative  is  obtained;  if  given in decreasing
 argument order a maximal length representative is obtained.

  - if  keyword `w` is given  the algorithm starts  from that `w`  (by default 
    the  algorithm starts from the representative  of the class obtained by
    `classinfo(WF)[:classtext][i]`)

  - if  `canonical=true` and the eigenvalue `1` appears last in the list  `e` 
    the function returns the "canonical" element of the class as an element
    of  the braid monoid  of `W` (in  other cases it  returns an element of
    `W`).

 The He-Nie algorithm does the following, for each `ζ∈ e` in turn:
  - compute the (ζ+ζ⁻¹)-eigenspace `V=V_{ζ}+V_{ζ⁻¹}` of `w+w⁻¹`.
  - find a chamber `w1` containing in its closure a regular vector of `V`
  - replace `w` by `w^w1` and `W` by `C_W(V)^w1`
 The list `e` can be obtained using `refleigen(W)[i]`.
"""
function HeConjugate(WF,i,e;w=classreps(WF)[i],canonical=false)
  sign(x)=iszero(x) ? 0 : Float64(x)<0 ? -1 : 1
  if !(WF isa Spets) WF=spets(WF) end
  W=Group(WF)
  B=BraidMonoid(W)
# set S to the matrix of the W-invariant scalar product
  S=toM(map((x,y)->x*y//2,rootlengths(W),eachrow(cartan(W))))
  xprintln("Starting with w=",B(w)," of length ",length(W,w))
  xprintln("w is in class#",i," with minimal length ",
    length(classinfo(WF)[:classtext][i])," and eigenvalues ",e)
  for l in e
    xprintln("\nUsing eigenvalue=",l)
    m=Cyc{Rational{Int}}.(reflrep(WF,w))
    V=lnullspace(m+m^-1-(l+l^-1)*one(m))
    println("   eigenspace V is ",size(V,1),"-dimensional")
    i=filter(i->V*reflrep(W,i)==V,1:nref(W))
    R=reflection_subgroup(W,i)
    xprintln("   C_W(V).w=",spets(R,w))
    v=nothing
    if issubset(inclusiongens(R),inclusiongens(W))
# try to see if we can find  a regular vector in Cbar
      m=V*permutedims(S)
      v=all_ge_1(m[(!iszero).(eachrow(m)),:];approx=x->Float64(x))
      if v===nothing
        xprintln(R," is standard but contains no regular vector in C̄")
      else v=permutedims(v)*V
      end
    end
# otherwise use random search for a regular vector
    if v===nothing
      while true
        v=permutedims(rand(-10:10,size(V,1)))*V
        if filter(i->v*reflrep(W,i)==v,1:W.N)==i break end
      end
    end
# find a chamber in whose adherence v sits.
    v=map(sign,toM(roots(W))*S*permutedims(v))
    if count(x->x==1,v)>count(x->x==-1,v)
      w1=with_inversions(W,filter(i->v[i]==-1,1:W.N))
    else 
      w1=with_inversions(W,filter(i->v[i]==1,1:W.N))
    end
  # d:=2*(W.N-R.N)*l;
    w=w^w1
    if isone(l) && canonical return B(w)*B(longest(W))^2 end
    W=R^w1
    WF=spets(W,w)
 #  if w1<>() and i<>[] then 
      xprintln("   conjugate to ",WF," by w1=",B(w1))
 #  fi;
    xprintln("   w1 conjugates w to ",B(w)," part in parabolic=",BraidMonoid(W)(w))
    if isempty(i) return B(w) end # regular element
  end
end

# return an element of the i-th conjugacy class of WF of minimal length.
MinConjugate(WF,i)=HeConjugate(WF,i,filter(x->x.r<=1//2,sort(unique(refleigen(WF,i)))));

# return an element of the i-th conjugacy class of WF of maximal length.
MaxConjugate(WF,i)=HeConjugate(WF,i,filter(x->x.r<=1//2,sort(unique(refleigen(WF,i)),rev=true)));

function CanonicalConjugate(WF,i)
  e=filter(x->x.r<=1//2,sort(unique(refleigen(WF,i))))
  if isone(e[1]) e=vcat(e[2:end],[e[1]]) end
  HeConjugate(WF,i,e;canonical=true)
end

# returns the length of the He conjugate
function LengthHe(WF,i,e;canonical=false)
#  if IsBound(arg[4]) then opt:=arg[4]; else opt:=rec();fi;
  if !(WF isa Spets) WF=spets(WF) end
  w=classreps(WF)[i]
  W=Group(WF)
  res=0;
  for l in e
    m=reflrep(WF,w)*E(1)//1
    V=lnullspace(m+m^-1-(l+l^-1)*one(m))
    i=filter(i->V*reflrep(W,i)==V,1:nref(W))
    R=reflection_subgroup(W,i)
    r=l.r
    if canonical && isone(l) r=1 end
    res+=2*(nref(W)-nref(R))*r
    W=R
    WF=spets(W,w)
    if isempty(i) return res end # regular element
  end
  res
end

function CanonicalLength(WF,i)
  e=unique(sort(refleigen(WF,i)))
  if isone(e[1]) e=vcat(e[2:end],[e[1]]) end
  LengthHe(WF,i,e;canonical=true)
end

# returns the composition product π_{I_i}^c_i for canonical elements
function CanonicalComposition(WF,i)
  w=CanonicalConjugate(WF,i)
  o=order(image(w))
  w=map(x->[x[1],x[2]//2o],tally(Brieskorn_normal_form(w^o)))
  sort(w,by=x->-length(x[1]))
end

# test for seq1<seq2
function cmpseq(WF,seq1,seq2)
  println("cmp(",FormatGAP(seq1),",",FormatGAP(seq2),")")
  if isempty(seq1) return true end
  a=Filtered(seq2,x->IsSubset(x[1],seq1[1][1]))
  seq2=Filtered(seq2,x->!IsSubset(x[1],seq1[1][1]))
  n=Sum(a,x->x[2]*coxnum(ReflectionSubgroup(WF,x[1])))/
    coxnum(ReflectionSubgroup(WF,seq1[1][1]))-seq1[1][2]
  if n<0 return false end
  cmpseq(WF,seq1[2:end],Concatenation([[seq1[1][1],n]],seq2))
end
