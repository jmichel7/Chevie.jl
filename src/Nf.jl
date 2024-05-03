"""
Number  fields are the finite extensions of  `ℚ`. The only ones that we can
handle  at the  moment are  subfields of  cyclotomic fields,  that is whose
elements  are `Cyc`s; they are also  characterized as the number fields `K`
such  that `Gal(K/ℚ)` is abelian.  For example, `ℚ (√5)`  is a number field
that is not cyclotomic but contained in the cyclotomic field `ℚ (ζ₅)`.

The default constructor for a number field takes some numbers as arguments and
constructs the smallest number field containing its arguments.

```julia-repl
julia> F=NF(E(5)) # the full cyclotomic field prints as CF
CF(5)

julia> K=NF(root(5)) # a subfield
NF(5,-1₅)

julia> conductor(K) # smallest n such that K is a subfield of CF(n)
5

julia> E(5)+E(5,-1) in NF(root(5)) # test if an element is in the subfield
true
```

A  number  field  `K`  is  printed  by  given the conductor of the smallest
cyclotomic field `F` containing it, and generators of the stabilizer of `K`
in  the galois group  of `F`. Above  `NF(5,-1₅)` represents the subfield of
`CF(5)` stable by complex conjugacy.

```julia-repl
julia> elements(galois(F))
4-element Vector{Chevie.Nf.NFAut}:
 Aut(CF(5),1₅)
 Aut(CF(5),2₅)
 Aut(CF(5),-1₅)
 Aut(CF(5),-2₅)
```

The  element of the galois  group of `CF(5)` printed  `-2₅` acts by raising
the  fifth roots of  unity to the  power -2. Thus  `-1₅` represents complex
conjugacy.

```julia-repl
julia> NF(root(3),root(5)) # here the stabilizer needs 2 generators
NF(60,-11₆₀,-1₆₀)
```
"""
module Nf
export NF, CF, Aut
using Primes
using ..Chevie
"""
`LenstraBase(n,stab,supergroup)`

returns  a list `[b₁,…,bₘ]`  of `Mod(-,n)` lists,  each `bᵢ` describing the
cyclotomic  `βᵢ=∑_{j∈  bᵢ}  ζₙʲ`  such  that  the  `βᵢ`  form a base of the
integer ring of the number field `NF(n,stab)`.

`supergroup` is a group containing `stab`; the base is chosen such that the
group of `supergroup` acts on it, as far as this is possible.

Note: the `bᵢ` are in general not sorted, since `bᵢ[1]` is often an element
of `zumbroich_basis(n)`.

`stab`  must not contain the stabilizer  of a proper cyclotomic subfield of
`CF(n)`.

```julia-repl
julia> Nf.LenstraBase(24,Group([Mod(19,24)]),Group([Mod(19,24)]))
4-element Vector{Vector{Mod{UInt64}}}:
 [1₂₄, -5₂₄]
 [8₂₄]
 [11₂₄, -7₂₄]
 [-8₂₄]

julia> Nf.LenstraBase(24,Group([Mod(19,24)]),Group([Mod(19,24),Mod(5,24)]))
4-element Vector{Vector{Mod{UInt64}}}:
 [1₂₄, -5₂₄]
 [5₂₄, -1₂₄]
 [8₂₄]
 [-8₂₄]

julia> Nf.LenstraBase(15,Group([Mod(4,15)]),Group(Mod.(prime_residues(15),15)))
4-element Vector{Vector{Mod{UInt64}}}:
 [1₁₅, 4₁₅]
 [2₁₅, -7₁₅]
 [7₁₅, -2₁₅]
 [-4₁₅, -1₁₅]
```
"""
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
    if isone(mod(Integer(k),4)) push!(H1,k)
    elseif iszero(Mod(k-1,No)) && 
      (iszero(Mod(k+1,N2)) || iszero(Mod(k+1-div(N2,2),N2)))
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
          setdiff!(zumb,map(k->(ppnt+k*div(N,d)).val,0:d-1))
        end
      end
    end
  end
  orbs
end

"""
`gens_prime_residues(n::Integer)` generators of multiplicative group `(ℤ /n)ˣ`

For each odd prime `p` there is one generator, corresponding to a primitive
root  of the  subgroup `(ℤ  /nₚ)ˣ` of  the multiplicative  group `(ℤ /n)ˣ`,
where `nₚ` is the highest power of `p` dividing `n`. For `p=2`, we have one
generator  unless `n₂`  is a  multiple of  8, in  which case  there are two
generators  5 and `n₂-1`. The function  returns a `Dict` recording for each
`p`, the list of generators of `(ℤ/nₚ)ˣ`.

```julia_repl
julia> Nf.gens_prime_residues( 24 )
Dict{Int64, Vector{Int64}} with 2 entries:
  2 => [7, 13]
  3 => [17]
```
"""
function gens_prime_residues(n::Integer)
  Dict(map(eachfactor(n))do (p,e)
    ppart=p^e
    rest=div(n,ppart)
    g=gcdx(ppart,rest)[3]*rest
    if p==2 gen=[-2g+1]
      if ppart>=8 push!(gen,4g+1) end
    else gen=[(primitiveroot(ppart)-1)*g+1]
    end
    p=>mod.(gen,n)
  end)
end

@GapObj struct NumberField
  degree::Int
  stabilizer::Group{Mod{UInt}} # stabilize gens in galois(CF(conductor))
  conductor::Int
end

"""
`NF(gens...)` or `NF(gens::AbstractVector)`

returns the smallest number field containing the elements `gens`, which may
be `Cyc`, `Root1`, `Integer` or `Rational{<:Integer}`.

```julia-repl
julia> NF(E(3),root(5))
NF(15,4₁₅)

julia> NF([E(3),root(5)])
NF(15,4₁₅)
```
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

NF(gens::Union{Cyc,Root1,Integer,Rational{<:Integer}}...)=NF(collect(Cyc.(gens)))

Base.:(==)(a::NumberField,b::NumberField)=conductor(a)==conductor(b) && 
  gens(a.stabilizer)==gens(b.stabilizer)

CyclotomicNumbers.conductor(F::NumberField)=F.conductor

function Base.in(c::Union{Cyc,Root1,Integer,Rational{<:Integer}},F::NumberField)
 (conductor(F)%conductor(c)==0) && all(i->galois(c,i)==c,gens(F.stabilizer))
end

iscyclotomic(F::NumberField)=istrivial(F.stabilizer)

function Base.show(io::IO,F::NumberField)
  if iscyclotomic(F) print(io,"CF(",conductor(F),")") 
  else print(io,"NF(",conductor(F),",");join(io,gens(F.stabilizer),",")
    print(io,")")
  end
end

"""
`CF(N::Integer)` the cyclotomic field generated by the `N`-th roots of unity.
"""
function CF(N)
  if N%4==2 N=div(N,2) end
  NumberField(totient(N),
     Group([Mod(1,N)]),conductor(E(N)),Dict{Symbol,Any}(:gens=>[Cyc(E(N))],
            :base=>Cyc.(E.(N,CyclotomicNumbers.zumbroich_basis(N)))))
end

struct NFAut
  F::NumberField
  galois::Mod
end

"""
`Aut(F::NumberField,k::Union{Integer,Mod})`

The  *Galois  automorphism*  `σₖ`  of  the  cyclotomic field `CF(n)` raises
`n`-th  roots of unity to the power `k`; it exists for `k` prime to `n`. If
`F`  is a subfield of `CF(n)`, the elements of the orbit of `σₖ` modulo the
stabilizer of `F` in the Galois group `galois(CF(n))` have same restriction
to `F`. An automorphism of `F` is represented by a canonical representative
`σₗ` of this orbit. This is the result of `Aut(F,k)`. The number `k` can be
given as an integer or as `Mod(k,n)`.

```julia-repl
julia> F=NF(root(5))
NF(5,-1₅)

julia> s=Aut(F,3)
Aut(NF(5,-1₅),2₅)

julia> root(5)^s # action of s on a Cyc
Cyc{Int64}: -√5
```
"""
function Aut(F::NumberField,g::Union{Integer,Mod})
  if g isa Integer g=Mod(g,conductor(F)) end
  NFAut(F,minimum(orbit(F.stabilizer,g,(e,g)->e*g)))
end

Base.:*(a::NFAut,b::NFAut)=Aut(a.F,a.galois*b.galois)
Base.one(a::NFAut)=Aut(a.F,1)
Base.copy(a::NFAut)=NFAut(a.F,a.galois)
Base.:^(x::NFAut,m::Integer)=m<0 ? Base.power_by_squaring(inv(x),-m) : Base.power_by_squaring(x,m) 
Base.:^(x::Cyc,a::NFAut)=galois(x,a.galois)
Base.inv(a::NFAut)=a^(order(a)-1)

function Base.show(io::IO,a::NFAut)
  print(io,"Aut(",a.F,",",a.galois,")")
end

function CyclotomicNumbers.galois(x::Union{Cyc,Integer,Rational{<:Integer}},y::Mod)
  if (y.n % conductor(x))!=0 error(y," not applicable to",x) end
  galois(x,Int(y.val))
end

"""
`galois(F::NumberField)` Galois group of `F` over `ℚ`

the  Galois group of `F`, a number  field of conductor `n`, is the quotient
of  the Galois  group of  `CF(n)`, isomorphic  to the  multiplicative group
`(ℤ/n)ˣ`,  by  the  stabilizer  of  `F`.  It  is given as a group of Galois
automorphisms (see `Aut`).

```julia-repl
julia> K=CF(5)
CF(5)

julia> F=NF(root(5))
NF(5,-1₅)

julia> galois(K)
Group(Chevie.Nf.NFAut[Aut(CF(5),2₅)])

julia> elements(galois(K))
4-element Vector{Chevie.Nf.NFAut}:
 Aut(CF(5),1₅)
 Aut(CF(5),2₅)
 Aut(CF(5),-1₅)
 Aut(CF(5),-2₅)

julia> elements(galois(F))
2-element Vector{Chevie.Nf.NFAut}:
 Aut(NF(5,-1₅),1₅)
 Aut(NF(5,-1₅),2₅)
```
"""
function CyclotomicNumbers.galois(F::NumberField)
  c=conductor(F)
  if c==1 return Group([Aut(F,1)]) end
  l=vcat(values(gens_prime_residues(c))...)
  g=Group(unique(NFAut.(Ref(F),Mod.(l,c))))
  Group(abelian_gens(g))
end

"""
`NF(N::Integer, stab::Group{<:Mod})`
fixed  field of the  subgroup `stab` of  `galois(CF(N))` in `CF(N)`. `stab`
should not inject in the multiplicative group of a proper divisor of `N`.
"""
function NF(N::Integer,stab::Group{<:Mod})
  d=exactdiv(totient(N),length(stab))

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
    if isone(mod(Integer(k),4)) push!(H1,k)
    elseif iszero(Mod(k-1,No)) && 
       (iszero(Mod(k+1,N2)) || iszero(Mod(k+1-div(N2,2),N2)))
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
  stab=Group(abelian_gens(stab))
  NumberField(d,stab,N,Dict{Symbol,Any}(:coeffslist=>coeffslist,:base=>base))
end
end
