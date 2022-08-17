module Nf
export NF, CF
using Primes
using ..Gapjm
#############################################################################
##
#F  LenstraBase(<N>,<stabilizer>,<super>) .   integral base of a number field
##
##  returns a list of lists of integers; each list indexing the exponents of
##  an orbit of a subgroup of <stabilizer> on <N>-th roots of unity.
##
##  *Note* that the elements are in general not sets, since the first element
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
function LenstraBase( N, stabilizer, supergroup )
  primes=keys(factor(N))
  NN=prod(primes)                # squarefree part of 'N'
  zumb=CyclotomicNumbers.zumbroich_basis(N) # exps of roots in base of 'CF(N)'
  stabilizer=sort(unique(stabilizer))
  if N==NN
  
    # 'N' is squarefree, we have the normal base, 'stabilizer' acts on
    # 'zumb'; do not consider equivalence classes since they are all
    # trivial.  'supergroup' is obsolete since 'zumb' is a normal base.
  
    # *Note* that for even 'N' 'zumb' does not consist of at least 'NN'-th
    # roots!
  
    orbits=Vector{Int}[]
    while !isempty(zumb)
      pnt=zumb[1]
      orb=mod.(pnt*stabilizer,N)
      push!(orbits,orb)
      setdiff!(zumb,orb)
    end
  
  else
  
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
    # Furthermore, the minimality of 'N' yields that 'stabilizer' acts fixed
    # point freely on the set of equivalence classes.
  
    # More exactly, fixed points occur exactly if there is an element 's' in
    # 'stabilizer' which is `≡ -1 (mod N2)` and `≡ 1 (mod No)`.
  
    # The base is constructed as follows:
    #
    # Until all classes are touched:
    # 1. Take a point 'pnt' (in 'zumb').
    # 2. Choose a maximal linear independent set 'pnts' in the equivalence
    #    class of 'pnt' (the intersection of the class with 'zumb').
    # 3. Take the 'stabilizer'--orbits of 'pnts' as base elements;
    #    remove the touched equivalence classes.
    # 4. For the representatives 'rep' in 'supergroup':
    #    If 'rep' maps 'pnt' to an equivalence class that was not yet
    #    touched, take the 'stabilizer'--orbits of the images of 'pnts'
    #    under 'rep' as base elements;
    #    remove the touched equivalence classes.
  
    # Compute nontriv. representatives of 'supergroup' over 'stabilizer'
    super=setdiff(supergroup,stabilizer)
    supergroup=Int[]
    while !isempty(super)
      pnt=super[1]
      push!(supergroup,pnt)
      setdiff!(super,mod.(stabilizer*pnt,N))
    end
  
    N2= 1; No= N
    while iseven(No)
      N2*=2; No=div(No,2)
    end
  
    H1=[]    # will be the subgroup of 'stabilizer' that acts fixed point
             # freely on the set of equivalence classes
  
    a=0
    for k in stabilizer
      if mod(k,4)==1 push!(H1,k)
      elseif mod(k-1,No)==0 && (mod(k+1,N2)==0 || mod(k+1-N2/2,N2)==0)
        a= k;
      end
    end
    if a==0 H1=stabilizer end
    orbits=[]
    orb=Int[]
    while !isempty(zumb)
      neworbits=[]
      pnt=zumb[1]
      d=1
      ord=div(N,gcd(N,pnt))
      for i in primes if mod(ord,i^2)==0 d*=i end end
      if a==0 || mod(ord,8)==0
        # the orbit of 'H1' cannot be a fixed point of 'a'
        for k in 0:d-1
          ppnt=pnt+k*div(N,d)
          if ppnt in zumb 
            orb=mod.(ppnt*stabilizer,N)
            push!(neworbits,orb)
          end
        end
      elseif mod(ord,4)==0
        # 'a' maps each point in the orbit of 'H1' to its inverse
        # (ignore all points)
        orb=mod.(pnt*stabilizer,N)
      else
        # the orbit of 'H1' is pointwise fixed by 'a'
        for k in 0:d-1
          ppnt=pnt+k*div(N,d)
          if ppnt in zumb 
            orb=mod.(ppnt*H1,N)
            push!(neworbits,orb)
          end
        end
      end
      for pnt in orb  # take any of the latest orbits
        # remove the equivalence class of 'pnt'
        setdiff!(zumb,map(k->mod(pnt+k*div(N,d),N),0:d-1))
      end
      append!(orbits,neworbits)
  
      # use 'supergroup':
      # Is there a point in 'zumb' not equivalent to '( pnt * rep ) mod N' ?
      # (Note that the factor group 'supergroup / stabilizer' acts on the
      # set of unions of orbits with equivalent elements.)
      for rep in supergroup
        # is there an 'x' in 'zumb' that is equivalent to 'pnt * rep' ?
        if any(x->mod((x-pnt*rep)*d,N)==0,zumb)
          append!(orbits,map(x->mod.(x*rep,N),neworbits))
          for ppnt in orbits[end]
            setdiff!(zumb,map(k->mod(ppnt+k*div(N,d),N),0:d-1))
          end
        end
      end
    end
  end
  orbits
end

"gens_primes_residues(n::Integer) generators of multiplicative group of ℤ /n"
function gens_primes_residues(n::Integer)
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

"""
`NF(gens)`  create number field generated by the elements of `gens`.
"""
function NF(gens::AbstractVector{<:Cyc})
  N=conductor(gens)
  if N==1 return CF(1) #Rationals
  elseif N==4 return CF(4) #GaussianRationals
  end
  stabilizer=filter(x->galois.(gens,x)==gens,prime_residues(N))
  NF(N,stabilizer,gens)
end

NF(gens::AbstractVector)=NF(Cyc.(gens))
#NF(gens::AbstractVector{<:Integer})=NF(Cyc.(gens))
#NF(gens::AbstractVector{<:Root1})=NF(Cyc.(gens))

@GapObj struct NumberField
  degree::Int
  generators::Vector{<:Cyc}
  base::Vector{<:Cyc}
  stabilizer::Vector{Int}
end

CyclotomicNumbers.conductor(F::NumberField)=conductor(F.generators)

function Base.in(c::Union{Cyc,Root1},F::NumberField)
  (conductor(F.generators)%conductor(c)==0) &&
  all(i->galois(c,i)==c,F.stabilizer[2:end])
end

function Base.show(io::IO,F::NumberField)
  if F.stabilizer==[1] print(io,"CF(",conductor(F),")")
  else print(io,"NF(",conductor(F),",",F.stabilizer[2:end],")")
  end
end

"""
`CF(N::Integer)` the cyclotomic field generated by the N`-th roots of unity.
"""
function CF(N)
  if N%4==2 N=div(N,2) end
  NumberField(totient(N),[Cyc(E(N))],
              Cyc.(E.(N,CyclotomicNumbers.zumbroich_basis(N))),
              [1],Dict{Symbol,Any}())
end

"""
NF(n, stab)
fixed field of the group generated by `stab` (prime residues modulo `n`)
in the cyclotomic field `CF(n)`,
"""
function NF(N::Integer,stabilizer::AbstractVector{<:Integer},gens=nothing)
  if N<=2 return CF(1) end
  if mod(N,4)==2 N=div(N,2) end
  stabilizer=copy(stabilizer)

  # Compute the elements of the group generated by 'stabilizer'.
  for d in stabilizer
    image=map(x->(x*d)%N,stabilizer)
    image=filter(x->!(x in stabilizer),image)
    append!(stabilizer, image)
  end

  d=exactdiv(totient(N),length(stabilizer))

  # reduce the pair '( N, stabilizer )' such that afterwards 'N'
  # describes the envelopping cyclotomic field of the required field;

  NN=1
  for gen in gens_primes_residues(N)
    if gen.prime==2
      if gen.exponent<3
        if !(gen.gen[1] in stabilizer) NN*=4 end
      else
        # the only case where length(gen.gen)>1
        # it contains the generators corresponding to '**' and '*5';
        # the first one is irrelevant for the envelopping cyclotomic
        # field, except if also the other generator is contained.
        if gen.gen[2] in stabilizer
          if !(gen.gen[1] in stabilizer) NN*=4 end
        else
          NN*=4
          aut= gen.gen[2]
          while !(aut in stabilizer)
            aut=mod(aut^2,N)
            NN*=2
          end
        end
      end
    else
      p=gen.prime
      if !(gen.gen[1] in stabilizer)
        NN*=p
        aut=mod(gen.gen[1]^(p-1),N)
        while !(aut in stabilizer)
          aut=mod(aut^p,N)
          NN*=p
        end
      end
    end
  end
  N=NN
  if N<=2 return CF(1) end
  stabilizer=sort(unique(mod.(stabilizer,N)))
  if stabilizer==[1] return CF(N) end

  # compute the standard Lenstra base and 'F.coeffslist':
  # If 'stabilizer' acts fixed point freely on the equivalence classes
  # we must change from the Zumbroich base to a 'stabilizer'--normal
  # base and afterwards choose coefficients with respect to that base.
  # In the case of fixed points, only the subgroup 'H1' of index 2 in
  # stabilizer acts fixed point freely; we change to a 'H1'--normal
  # base and afterwards choose coefficients.

  N2=1; No=N
  while iseven(No) N2*=2; No=div(No,2) end
  H1=[]    # will be the subgroup of 'stabilizer' that acts fixed
           # point freely on the set of equivalence classes
  a=0
  for k in stabilizer
    if k%4==1 push!(H1,k)
    elseif (k-1)%No== 0 && ( (k+1)%N2==0 || (k+1-N2/2)%N2==0)
      a=k
    end
  end
  if a ==0 H1=stabilizer end

  zumb=CyclotomicNumbers.zumbroich_basis(N)
  lenst=LenstraBase(N,H1,stabilizer)
  
  # We want 'Sublist(CoeffsCyc(z,N),F.coeffslist) = Coefficients(F,z)'
  # ( and   'Coefficients( F, z ) * F.base = z' )
  # with respect to the standard Lenstra base.

  if H1!=stabilizer # let 'a' act on 'lenst' to get the right base
    newlenst=[]
    lenstset=sort.(unique.(lenst))
    for i in eachindex(lenst)
      if (lenst[i][1]*(a-1))% N == 0 # pointwise fixed
        push!( newlenst, lenst[i] )
      elseif ( lenst[i][1] * ( a - 1 ) - N/2 ) % N != 0
        # 'a' joins two 'H1'--orbits
        image=sort(unique(map(x->(x*a)%N,lenst[i])))
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
  base=map(x->Cyc(sum(E(N).^x)),lenst)
  if isnothing(gens) gens=copy(base) end
  res=NumberField(d,gens,base,stabilizer,Dict{Symbol,Any}())
  res.coeffslist=coeffslist
  res
end

##############################################################################
#F  NumberField( <subfield>, xtension)   xtension is base or poly
##
#NF := function (subfield,xtension::AbstractVector)
#    
#  
#      if ( subfield = 0 and ForAll( xtension, IsRat ) ) or
#         ( IsNumberField( subfield ) and ForAll( xtension,
#                                                 x -> x in subfield ) ) then
#  
#        # NF( subfield, poly )
#  
#        if Length( xtension ) > 3 then
#          Error("NF(<subfield>,<poly>) for polynomial of degree at most 2");
#        fi;
#        if Length( xtension ) = 2 then   # linear polynomial
#
#          if subfield = 0 then
#            return Rationals;
#          else
#            return NF( subfield.generators );
#          fi;
#
#        else
#          
#          # The roots of 'a*x^2 + b*x + c' are
#          # $\frac{ -b \pm \sqrt{ b^2 - 4ac } }{2a}$.
#  
#          root:= ( ER( xtension[2]^2 - 4*xtension[1]*xtension[3] )
#                         -xtension[2] ) / 2*xtension[3];
#  
#          if subfield = 0 then
#            return NF( [ root ] );
#          elif root in subfield then
#            return NF( subfield, subfield.base );
#          else
#            return NF( subfield, [ 1, root ] );
#          fi;
#
#        fi;
#
#      else
#  
#        # 'NF( 0, base )' or 'NF( subfield, base )'
#        # where 'base' at least contains a vector space base
#        if subfield = 0 or subfield = Rationals then
#          F:= NF( xtension );
#        else
#          F:= NF( Union( subfield.generators, xtension ) );
#          if IsCyclotomicField( F ) then
#            return CF( subfield, xtension );
#          fi;
# 
#          # general case\:\ extension of number field; do not
#          #                 ask if <subfield> is a cyclotomic field 
#
#          # Let $(v_1, \ldots, v_m)$ denote 'subfield.base' and
#          #     $(w_1, \ldots, w_k)$ denote 'F.base';
#          # Define $u_{i+m(j-1)} = v_i w_j$.  Then $(u_l; 1\leq l\leq mk)$
#          # is a $Q$--base of 'F'.  First change from the Lenstra base to
#          # this base; the matrix is 'C'\:
#
#          F.dimension    := F.degree / subfield.degree;
#          F.field        := subfield;
#          F_base         := NormalBaseNumberField( F );
#
#          m:= Length( subfield.base );
#          k:= Length( F_base );
#          N:= NofCyc( F_base );
#          C:= [];
#          for j in F_base do
#            for i in subfield.base do
#              Add( C, F.operations.Coefficients( F, i*j ) );
#              # (These are the Lenstra base coefficients!)
#            od;
#          od;
#          C:= C^(-1);
#
#          # Let $(c_1, \ldots, c_{mk})$ denote the coefficients with respect
#          # to the new base.  To achieve '<coeffs> \* F_base = <z>' we have
#          # to take $\sum_{i=1}^m c_{i+m(j-1)} v_i$ as $j$--th coefficient\:
#
#          F.coeffsmat:= [];
#          for i in [ 1 .. Length( C ) ] do     # for all rows
#            F.coeffsmat[i]:= [];
#            for j in [ 1 .. k ] do
#              val:= 0;
#              for l in [ 1 .. m ] do
#                val:= val + C[i][ m*(j-1)+l ] * subfield.base[l];
#              od;
#              F.coeffsmat[i][j]:= val;
#            od;
#          od;
#
#          # Multiplication of a Lenstra base coefficient vector with
#          # 'F.coeffsmat' means first changing to the base of products
#          # $v_i w_j$ and then summation over the parts of the $v_i$.
#
#          F.base         := F_base;
#          F.isNormalBase := true;
#          Unbind( F.isIntegralBase );
#          Unbind( F.name );
#
#        fi;
#  
#        # If 'xtension' just contains a base but is no base, do nothing;
#        # 'xtension' is a base if and only if the basechange is regular:
#        if Length( xtension ) = Length( F.base ) and xtension <> F.base then
#
#          F.basechangemat:= List( xtension,
#                                x -> F.operations.Coefficients(F,x) )^(-1);
#          F.base         := Copy( xtension );
#          Unbind( F.isNormalBase );
#          Unbind( F.name );
#        fi;
#      fi;
#    fi;
#    
#    # return the number field record
#    return F;
#    end;
end
