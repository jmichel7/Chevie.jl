module Nf
export NF, CF, NFAut
using Primes
using ..Gapjm
#############################################################################
#F  LenstraBase(<N>,<stabilizer>,<super>) .   integral base of a number field
##
##  returns a list of lists of integers; each list indexing the exponents of
##  an orbit of a subgroup of <stabilizer> on <N>-th roots of unity.
##
##  *Note* that the lists are in general not sorted, since the first element
##  is always an element of 'zumbroich_basis(N)'; this is used by 'NF' and
##  'Coefficients'.
##
##  <super> is a list representing a supergroup of <stabilizer> which
##  shall act consistently with the action of <stabilizer>, i.e. each orbit
##  of <supergroup> is a union of orbits of <stabilizer>.
##
##  ( Shall there be a test if this is possible ? )
##
##  *Note* that <stabilizer> must not contain the stabilizer of a proper
##  cyclotomic subfield of the <N>-th cyclotomic field.
##
function LenstraBase( N, stab, supergroup )
  primes=keys(factor(N))
  NN=prod(primes)                # squarefree part of 'N'
  zumb=CyclotomicNumbers.zumbroich_basis(N) # exps of roots in base of 'CF(N)'

  # if N==NN then N is squarefree, we have the normal base, 'stab' acts on
  # 'zumb'; do not consider equivalence classes since they are all
  # trivial.  'supergroup' is obsolete since 'zumb' is a normal base.

  # *Note* that for even 'N' 'zumb' does not consist of at least 'NN'-th
  # roots!
  if N==NN return orbits(stab,Mod.(zumb,N),(e,g)->e*g) end
  
  # Let `d(i)` be the largest squarefree number whose square divides the
  # order of `e_N^i`, that is 'N/gcd(N,i)'.
  # Define an equivalence relation on the set `S` of at least 'NN'-th
  # roots of unity:
  # `i` and `j` are equivalent iff 'N' divides `(i-j)d(i)`.  The
  # equivalence class `(i)` of `i` is `\{ i+kN/d(i) ; 0≤k<d(i)\}`.

  # For the case that 'NN' is even, replace those roots in `S` with order
  # not divisible by 4 by their negatives. (Equivalently\: Replace *all*
  # elements in `S` by their negatives.)

  # If 8 does not divide 'N' and `N≠4`, 'zumb' is a subset of `S`,
  # the intersection of `(i)` with 'zumb' is of order `φ(d(i))`,
  # it is a basis for the `Z`--submodule spanned by `(i)`.
  # Furthermore, the minimality of 'N' yields that 'stab' acts fixed
  # point freely on the set of equivalence classes.

  # More exactly, fixed points occur exactly if there is an element 's' in
  # 'stab' which is `≡ -1 (mod N2)` and `≡ 1 (mod No)`.

  # The base is constructed as follows:
  #
  # Until all classes are touched:
  # 1. Take a point 'pnt' (in 'zumb').
  # 2. Choose a maximal linear independent set 'pnts' in the equivalence
  #    class of 'pnt' (the intersection of the class with 'zumb').
  # 3. Take the 'stab'--orbits of 'pnts' as base elements;
  #    remove the touched equivalence classes.
  # 4. For the representatives 'rep' in 'supergroup':
  #    If 'rep' maps 'pnt' to an equivalence class that was not yet
  #    touched, take the 'stab'--orbits of the images of 'pnts'
  #    under 'rep' as base elements;
  #    remove the touched equivalence classes.

  # Compute nontriv. representatives of 'supergroup' over 'stab'
  super=first.(orbits(stab,setdiff(elements(supergroup),elements(stab)),(e,g)->e*g))

  N2=1; No=N; while iseven(No) N2*=2; No=div(No,2) end

  H1=eltype(stab)[]  # will be the subgroup of 'stab' that acts fixed point
           # freely on the set of equivalence classes
  a=Mod(0,N)
  for k in elements(stab)
   if mod(Integer(k),4)==1 push!(H1,k)
   elseif mod(Integer(k)-1,No)==0 && (mod(Integer(k)+1,N2)==0 || mod(Integer(k)+1-N2/2,N2)==0)
      a=k;
    end
  end
  if iszero(a) H1=stab else H1=Group(H1) end
  orbs=Vector{eltype(stab)}[]
  orb=eltype(stab)[]
  while !isempty(zumb)
    neworbs=Vector{eltype(stab)}[]
    pnt=zumb[1]
    d=1
    ord=div(N,gcd(N,pnt))
    for i in primes if mod(ord,i^2)==0 d*=i end end
    if iszero(a) || mod(ord,8)==0
      # the orbit of 'H1' cannot be a fixed point of 'a'
      for k in 0:d-1
        ppnt=pnt+k*div(N,d)
        if ppnt in zumb 
         orb=orbit(stab,Mod(ppnt,N),(e,g)->e*g);push!(neworbs,orb)
        end
      end
    elseif mod(ord,4)==0
      # 'a' maps each point in the orbit of 'H1' to its inverse
      # (ignore all points)
      orb=orbit(stab,Mod(pnt,N),(e,g)->e*g)
    else
      # the orbit of 'H1' is pointwise fixed by 'a'
      for k in 0:d-1
        ppnt=pnt+k*div(N,d)
        if ppnt in zumb 
         orb=orbit(H1,Mod(ppnt,N),(e,g)->e*g);push!(neworbs,orb)
        end
      end
    end
    for pnt in orb  # take any of the latest orbits
      # remove the equivalence class of 'pnt'
      setdiff!(zumb,map(k->mod(Integer(pnt)+k*div(N,d),N),0:d-1))
    end
    append!(orbs,neworbs)

    # use 'super':
    # Is there a point in 'zumb' not equivalent to '( pnt * rep ) mod N' ?
    # (Note that the factor group 'supergroup / stab' acts on the
    # set of unions of orbits with equivalent elements.)
    for rep in super
      # is there an 'x' in 'zumb' that is equivalent to 'pnt * rep' ?
      if any(x->iszero((x-pnt*rep)*d),zumb)
        append!(orbs,map(x->x.*rep,neworbs))
        for ppnt in orbs[end]
          setdiff!(zumb,Integer.(map(k->ppnt+k*div(N,d),0:d-1)))
        end
      end
    end
  end
  orbs
end

"""
`gens_prime_residues(n::Integer)` generators of multiplicative group of ℤ /n
"""
function gens_prime_residues(n::Integer)
  map(eachfactor(n))do (p,e)
    ppart=p^e
    rest=div(n,ppart)
    g=gcdx(ppart,rest)[3]*rest
    if p==2 gen=[-2g+1]
      if ppart>=8 push!(gen,4g+1) end
    else gen=[(primitiveroot(ppart)-1)*g+1]
    end
    (prime=p,exponent=e,gen=mod.(gen,n))
  end
end

@GapObj struct NumberField
  degree::Int
  stabilizer::Group{Mod{UInt}} # galois-stabilize gens in prime_residues(conductor)
  conductor::Int
end

"""
`NF(gens)`  create number field generated by the elements of `gens`.
"""
function NF(gens::AbstractVector{<:Cyc})
  gens=unique(gens)
  cond=conductor(gens)
  if cond==1 return CF(1) #Rationals
  elseif cond==4 return CF(4) #GaussianRationals
  end
  stab=Group(Mod.(filter(x->galois.(gens,x)==gens,prime_residues(cond)),cond))
  res=NF(cond,stab)
  res.gens=gens
  res
end

NF(gens::AbstractVector)=isempty(gens) ? CF(1) : NF(Cyc.(gens))

Base.:(==)(a::NumberField,b::NumberField)=conductor(a)==conductor(b) && 
  a.stabilizer==b.stabilizer

CyclotomicNumbers.conductor(F::NumberField)=F.conductor

function Base.in(c::Union{Cyc,Root1,Integer,Rational{<:Integer}},F::NumberField)
 (conductor(F)%conductor(c)==0) && all(i->galois(c,i)==c,gens(F.stabilizer))
end

function Base.show(io::IO,F::NumberField)
  if istrivial(F.stabilizer) print(io,"CF(",conductor(F),")")
  else print(io,"NF(",conductor(F),",")
       join(io,gens(F.stabilizer),",")
       print(io,")")
  end
end

"""
`CF(N::Integer)` the cyclotomic field generated by the N`-th roots of unity.
"""
function CF(N)
  if N%4==2 N=div(N,2) end
  NumberField(totient(N),
     Group([Mod(1,N)]),conductor(E(N)),Dict{Symbol,Any}(:gens=>[Cyc(E(N))]),
     :base=>Cyc.(E.(N,CyclotomicNumbers.zumbroich_basis(N))))
end

struct NFAut
  F::NumberField
  galois::Mod
  function NFAut(F::NumberField,g::Union{Integer,Mod})
    if g isa Integer g=Mod(g,conductor(F)) end
    new(F,minimum(orbit(F.stabilizer,g,(e,g)->e*g)))
  end
end

Base.:*(a::NFAut,b::NFAut)=NFAut(a.F,a.galois*b.galois)
Base.one(a::NFAut)=NFAut(a.F,1)
Base.copy(a::NFAut)=NFAut(a.F,a.galois)
Base.:^(x::NFAut,m::Integer)=Base.power_by_squaring(x,m)
Base.:^(x::Cyc,a::NFAut)=galois(x,a.galois)

function Base.show(io::IO,a::NFAut)
  print(io,"Aut(",a.F,",",a.galois,")")
end

function CyclotomicNumbers.galois(x::Cyc,y::Mod)
  if (y.n % conductor(x))!=0 error(y," not applicable to",x) end
  galois(x,Int(y.val))
end

function CyclotomicNumbers.galois(F::NumberField)
  c=conductor(F)
  l=vcat(map(x->x.gen,gens_prime_residues(c))...)
  g=Group(unique(NFAut.(Ref(F),l)))
  Group(abelian_gens(g))
end

"""
NF(N, stab)
fixed field of the group generated by `stab` (prime residues modulo `n`)
in the cyclotomic field `CF(N)`.
"""
function NF(N::Integer,stab::Group{<:Mod})
  d=exactdiv(totient(N),length(stab))

  # reduce the pair '( N, stab )' such that afterwards 'N'
  # describes the envelopping cyclotomic field of the required field;

  if istrivial(stab) return CF(N) end

  # compute the standard Lenstra base and 'F.coeffslist':
  # If 'stab' acts fixed point freely on the equivalence classes
  # we must change from the Zumbroich base to a 'stab'--normal
  # base and afterwards choose coefficients with respect to that base.
  # In the case of fixed points, only the subgroup 'H1' of index 2 in
  # stab acts fixed point freely; we change to a 'H1'--normal
  # base and afterwards choose coefficients.

  N2=1; No=N; while iseven(No) N2*=2; No=div(No,2) end
  H1=eltype(stab)[]  # will be the subgroup of 'stab' that acts fixed
           # point freely on the set of equivalence classes
  a=Mod(0,N)
  for k in elements(stab)
   if Integer(k)%4==1 push!(H1,k)
   elseif (Integer(k)-1)%No==0 && 
          ((Integer(k)+1)%N2==0 || (Integer(k)+1-N2/2)%N2==0)
      a=k
   end
  end
  if iszero(a) H1=stab else H1=Group(H1) end

  zumb=CyclotomicNumbers.zumbroich_basis(N)
  lenst=LenstraBase(N,H1,stab)
  
  # We want 'Sublist(CoeffsCyc(z,N),F.coeffslist) = Coefficients(F,z)'
  # ( and   'Coefficients( F, z ) * F.base = z' )
  # with respect to the standard Lenstra base.

  if length(H1)!=length(stab) # let 'a' act on 'lenst' to get the right base
    newlenst=Vector{eltype(stab)}[]
    lenstset=sort.(unique.(lenst))
    for i in eachindex(lenst)
      if iszero(lenst[i][1]*(a-1)) # pointwise fixed
        push!(newlenst, lenst[i])
      elseif !iszero(lenst[i][1]*(a-1)-div(N,2))
        # 'a' joins two 'H1'--orbits
        image=sort(unique(lenst[i].*a))
        # *Note* that the elements of 'image' need not be in an element
        # of 'lenst', only a member of the equivalence class must be
        # contained;
        if findfirst(==(image),lenstset)>=i 
          push!(newlenst,vcat(lenst[i],image))
        end
      end
    end
    lenst=newlenst
  end
  coeffslist=map(x->x[1]+1,lenst)
  base=map(x->Cyc(sum(E.(N,map(y->y.val,x)))),lenst)
  NumberField(d,stab,N,Dict{Symbol,Any}(:coeffslist=>coeffslist,:base=>base))
end
end
