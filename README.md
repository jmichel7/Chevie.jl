
<a id='Gapjm.jl-Documentation-1'></a>

# Gapjm.jl Documentation

<a id='Gapjm' href='#Gapjm'>#</a>
**`Gapjm`** &mdash; *Module*.



Here are some of my efforts porting GAP code to julia. I am not even sure this is a well-formed Julia package.

It  contains  for  now  permutations  and  permutation  groups,  cyclotomic numbers,  Laurent polynomials, some  Weyl groups and  Coxeter groups, Hecke algebras  and Kazhdan-Lusztig polynomials. Coming soon are braid groups and Garside monoids and factorisations into cyclotomic polynomials.

Even  though the code  is often competitive  with or faster  than GAP, I am sure there are more optimisations possible. Any comments about the code and the design are welcome.

If you are new to julia, to install this package:

  * enter package mode with ]
  * do the command

```
(v1.0) pkg> add "https://github.com/jmichel7/Gapjm.jl"
```

  * exit package mode with backspace

do 

```
julia> using Gapjm
```

and you are set up.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Gapjm.jl#L1-L26' class='documenter-source'>source</a><br>

- [Gapjm.jl Documentation](README.md#Gapjm.jl-Documentation-1)
- [Perms.jl Documentation](README.md#Perms.jl-Documentation-1)
- [PermGroups.jl Documentation](README.md#PermGroups.jl-Documentation-1)
- [Cycs.jl Documentation](README.md#Cycs.jl-Documentation-1)
- [Pols.jl Documentation](README.md#Pols.jl-Documentation-1)
- [CoxGroups.jl Documentation](README.md#CoxGroups.jl-Documentation-1)
- [Weyl.jl Documentation](README.md#Weyl.jl-Documentation-1)
- [Hecke.jl Documentation](README.md#Hecke.jl-Documentation-1)
- [KL.jl Documentation](README.md#KL.jl-Documentation-1)
- [Util.jl Documentation](README.md#Util.jl-Documentation-1)


<a id='Perms.jl-Documentation-1'></a>

# Perms.jl Documentation

<a id='Gapjm.Perms' href='#Gapjm.Perms'>#</a>
**`Gapjm.Perms`** &mdash; *Module*.



This module is a port of some GAP functionality on permutations.

A  permutation here is a permutation of the set 1:n and is represented as a list  of n integers representing the images of 1:n. The integer n is called the *degree* of the permutation.

Permutations  in  this  module  follow  the  GAP  design: it is possible to multiply, or to store in the same group, permutations of different degrees. A  slightly faster  design is  the MAGMA  one where  any permutation has to belong  to  a  group  and  the  degree  is  determined by that group. There multiplication of permutations in a given group is a faster, but it is more difficult  to multiply  permutations coming  from different  groups, like a group and one of its subgroups.

The  GAP permutation  (1,2,3)(4,5) can  be written Perm(1,2,3)*Perm(4,5) or perm"(1,2,3)(4,5)".  It is represented internally as [2,3,1,5,4]; note that [2,3,1,5,4,6] represents the same permutation.

As in GAP i^p applies p to integer i, while p^q means p^-1*q&ast;p.

Another  Perm  constructor  is  Perm{T}(p) which converts the perm  p to a permutation on integers of type T; for instance Perm{UInt8} is more  efficient that Perm{Int} and can be  used for Weyl groups of rank <=8 since they have at most 240 roots.

**Examples**

```julia-repl
julia> p=Perm(1,2)*Perm(2,3)
(1,3,2)

julia> Perm{Int8}(p)
{Int8}(1,3,2)

julia> 1^p
3

julia> p^Perm(3,10)
(1,10,2)

julia> inv(p)
(1,2,3)

julia> one(p)
()

julia> order(p)
3

julia> degree.((Perm(1,2),Perm(2,3)))
(2, 3)

julia> largest_moved_point(Perm(1,2)*Perm(2,3)^2)
2

julia> smallest_moved_point(Perm(2,3))
2

julia> Matrix(p)
3×3 Array{Float64,2}:
 0.0  0.0  1.0
 1.0  0.0  0.0
 0.0  1.0  0.0

julia> Matrix{Int}(p)
3×3 Array{Int64,2}:
 0  0  1
 1  0  0
 0  1  0
```

Perms  have methods copy, hash,  ==, cmp, isless (total order)  so they can be keys in hashes or elements of sets.

other functions are: cycles, cycletype, sign. See individual documentation.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Perms.jl#L1-L76' class='documenter-source'>source</a><br>

<a id='Gapjm.Perms.cycles' href='#Gapjm.Perms.cycles'>#</a>
**`Gapjm.Perms.cycles`** &mdash; *Function*.



cycles(a::Perm) returns the non-trivial cycles of a

**Example**

```julia-repl
julia> cycles(Perm(1,2)*Perm(4,5))
3-element Array{Array{Int64,1},1}:
 [1, 2]
 [3]
 [4, 5]
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Perms.jl#L184-L194' class='documenter-source'>source</a><br>

<a id='Gapjm.Perms.cycletype' href='#Gapjm.Perms.cycletype'>#</a>
**`Gapjm.Perms.cycletype`** &mdash; *Function*.



cycletype(a::Perm) is the partition of degree(a) associated to the   conjugacy class of a in the symmetric group, with ones removed

**Example**

```julia-repl
julia> cycletype(Perm(1,2)*Perm(3,4))
2-element Array{Int64,1}:
 2
 2

```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Perms.jl#L223-L234' class='documenter-source'>source</a><br>

<a id='Base.sign' href='#Base.sign'>#</a>
**`Base.sign`** &mdash; *Function*.



sign(a::Perm) is the signature of  the permutation a


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Perms.jl#L237' class='documenter-source'>source</a><br>


<a id='PermGroups.jl-Documentation-1'></a>

# PermGroups.jl Documentation

<a id='Gapjm.PermGroups' href='#Gapjm.PermGroups'>#</a>
**`Gapjm.PermGroups`** &mdash; *Module*.



This module is a port of some GAP functionality on groups, in particular permutation groups.

This  codes refers to Holt "Handbook of computational group theory" chapter 4 for basic algorithms.

The only field of a Group G at the start is gens, the list of generators of G.  To  mimic  GAP  records  where  attributes/properties  of an object are computed  on demand when asked for, other fields are computed on demand and stored in the field prop of the Group, which starts as Dict{Symbol,Any}()

A  PermGroup is  a group  where gens  are Perms,  which allows  for all the algorithms like base, centralizer chain, etc...

**Examples**

```julia-repl
julia> G=PermGroup([Perm(i,i+1) for i in 1:2])
PermGroup((1,2),(2,3))

julia> collect(G)  # PermGroups are iterators over their elements
6-element Array{Perm{Int64},1}:
 (1,2)
 (1,3,2)
 ()
 (1,2,3)
 (1,3)
 (2,3)

julia> degree(G)  # maximum degree of an element of G
3

julia> orbit(G,1) # orbit of point 1 under G
3-element Array{Int64,1}:
 1
 2
 3

# orbit decorated with representatives moving 1 to given point
julia> orbit_and_representative(G,1)
Dict{Int64,Perm{Int64}} with 3 entries:
  2 => (1,2)
  3 => (1,3,2)
  1 => ()

# this function is general and can take any action
julia> orbit_and_representative(G,[1,2],(x,y)->x.^Ref(y))
Dict{Array{Int64,1},Perm{Int64}} with 6 entries:
  [1, 3] => (2,3)
  [1, 2] => ()
  [2, 3] => (1,2,3)
  [3, 2] => (1,3)
  [2, 1] => (1,2)
  [3, 1] => (1,3,2)

julia> Perm(1,2) in G
true

julia> Perm(1,2,4) in G
false

# Elements,  appartenance test and  other function are  computed on G using
# Schreier-Sims theory, that is computing the following

julia> base(G) # a list of points that no element of G fixes
2-element Array{Int64,1}:
 1
 2

julia> centralizers(G) # the i-th element is the centralizer of base[1:i-1]
2-element Array{PermGroup{Int64},1}:
 PermGroup((1,2),(2,3))
 PermGroup((2,3))

# i-th element is orbit_and_representive of centralizer[i] on base[i]
julia> centralizer_orbits(G)
2-element Array{Dict{Int64,Perm{Int64}},1}:
 Dict(2=>(1,2),3=>(1,3,2),1=>())
 Dict(2=>(),3=>(2,3))

julia> words(G)  # minimal word for each element of G
6-element Array{Array{Int64,1},1}:
 []
 [1]
 [2]
 [1, 2]
 [2, 1]
 [1, 2, 1]

julia> elements(G) # elements in the same order as words
6-element Array{Perm{Int64},1}:
 ()
 (1,2)
 (2,3)
 (1,3,2)
 (1,2,3)
 (1,3)
```

finally, benchmarks on julia 1.0.1

```benchmark
julia> @btime length(collect(symmetric_group(8)))
  6.519 ms (350522 allocations: 14.17 MiB)

julia> @btime words(symmetric_group(8));
  10.477 ms (122062 allocations: 15.22 MiB)
```

Compare to GAP3 Elements(SymmetricGroup(8)); takes 3.8 ms


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/PermGroups.jl#L1-L110' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.symmetric_group' href='#Gapjm.PermGroups.symmetric_group'>#</a>
**`Gapjm.PermGroups.symmetric_group`** &mdash; *Function*.



The symmetric group of degree n 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/PermGroups.jl#L241' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.orbit' href='#Gapjm.PermGroups.orbit'>#</a>
**`Gapjm.PermGroups.orbit`** &mdash; *Function*.



orbit(G,p) is the orbit of Int p under Group G


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/PermGroups.jl#L127' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.orbit_and_representative' href='#Gapjm.PermGroups.orbit_and_representative'>#</a>
**`Gapjm.PermGroups.orbit_and_representative`** &mdash; *Function*.



returns Dict x=>g for x in orbit(G,p) and g is such that x=action(p,g)


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/PermGroups.jl#L160' class='documenter-source'>source</a><br>

<a id='Gapjm.words' href='#Gapjm.words'>#</a>
**`Gapjm.words`** &mdash; *Function*.



List of minimal words in the generators elements(G) 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/PermGroups.jl#L214' class='documenter-source'>source</a><br>

<a id='Gapjm.elements' href='#Gapjm.elements'>#</a>
**`Gapjm.elements`** &mdash; *Function*.



The list of elements of G in the same order as words


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/PermGroups.jl#L219' class='documenter-source'>source</a><br>


<a id='Cycs.jl-Documentation-1'></a>

# Cycs.jl Documentation

<a id='Gapjm.Cycs' href='#Gapjm.Cycs'>#</a>
**`Gapjm.Cycs`** &mdash; *Module*.



Cyclotomic  numbers means complex numbers which are sums of rationals times roots of unity.

They are a very important feature of GAP, since entries of character tables of finite groups are cyclotomics.

They  have a normal form given by the Zumbroich basis, which allows to find the  smallest Cyclotomic field which contains a given number, and decide in particular if a cyclotomic is zero. Let ζ*n:=e^{2iπ/n}. The Zumbroich basis of Q(ζ*n) is a particular subset of 1,ζ,ζ^2,...,ζ^{n-1} which forms a basis of Q(ζ_n) with good properties.

I  ported here Christian Stump's Sage  code, which is simpler to understand than GAP's code. The reference for the algorithms is

T. Breuer, Integral bases for subfields of cyclotomic fields AAECC 8 (1997)

As in GAP, I lower automatically numbers after each computation; this makes this code about twice slower than GAP since lower is not as much optimized. GAP  also converts a Cyclotomic which is rational to a Rational, a Rational which is integral to an Int, etc... This is tremendously useful but needs a new  type of  number to  be added  to Julia,  which requires more competent people than me.

The main way to build a Cyclotomic number is to use the function `E(n,k=1)` which constructs ζ_n^k.

**Examples**

```julia-repl
julia> E(3)+E(4)
ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹

julia> E(3,2)
ζ₃²

julia> 1+E(3,2)
-ζ₃

julia> a=E(4)-E(4)
0

julia> conductor(a) # a is lowered to Q(ζ_1)=Q
1

julia> typeof(convert(Int,a))
Int64

julia> convert(Int,E(4))
ERROR: InexactError: convert(Int64, E(4))

julia> c=inv(1+E(4)) # inverses need Rationals
1//2+(-1//2)ζ₄

julia> typeof(c)
Cyc{Rational{Int64}}

julia> typeof(1+E(4))
Cyc{Int64}

julia> Cyc(1+im) # one can convert Gaussian integers or rationals
1+ζ₄

julia> 1//(1+E(4))
1//2+(-1//2)ζ₄

julia> typeof(Cyc(1//2)) # another way of building a Cyc
Cyc{Rational{Int64}}

julia> conj(1+E(4))
1-ζ₄

julia> c=E(9)   # an effect of the Zumbroich basis
-ζ₉⁴-ζ₉⁷

julia> Root1(c) # but you can decide whether a Cyc is a root of unity
Root1(1//9)

julia> c=Complex(E(3))   # convert to float is probably not very useful
-0.4999999999999998 + 0.8660254037844387im

julia> Cyc(c) # even less useful
-0.4999999999999998+0.8660254037844387ζ₄
```

For more information see ER, quadratic, galois. 

Finally, a benchmark:

```benchmark
julia> function testmat(p) 
         ss=vcat([[[i,j] for j in i+1:p-1] for i in 0:p-1]...)
         [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
       end
testmat (generic function with 1 method)

julia> @btime testmat(12)^2;
  261.577 ms (3012634 allocations: 290.07 MiB)
```

The equivalent in GAP:

```
testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
end; 
```

for testmat(12) takes 0.4s in GAP3, 0.3s in GAP4


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Cycs.jl#L1-L109' class='documenter-source'>source</a><br>

<a id='Gapjm.Cycs.galois' href='#Gapjm.Cycs.galois'>#</a>
**`Gapjm.Cycs.galois`** &mdash; *Function*.



galois(c::Cyc,n::Int) applies to c the galois automorphism   of Q(ζ_conductor(c)) raising all roots of unity to the n-th power.   n should be prime to c.

**Examples**

```julia-repl
julia> galois(1+E(4),-1) # galois(c,-1) is the same as conj(c)
1-ζ₄

julia> galois(ER(5),2)==-ER(5)
true
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Cycs.jl#L410-L422' class='documenter-source'>source</a><br>

<a id='Gapjm.Cycs.ER' href='#Gapjm.Cycs.ER'>#</a>
**`Gapjm.Cycs.ER`** &mdash; *Function*.



ER(n::Int) computes as a Cyc the square root of the integer n.

**Examples**

```julia-repl
julia> ER(-3)
ζ₃-ζ₃²

julia> ER(3)
-ζ₁₂⁷+ζ₁₂¹¹
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Cycs.jl#L444-L454' class='documenter-source'>source</a><br>

<a id='Gapjm.Cycs.quadratic' href='#Gapjm.Cycs.quadratic'>#</a>
**`Gapjm.Cycs.quadratic`** &mdash; *Function*.



quadratic(c::Cyc) determines if c lives in a quadratic extension of Q   it  returns a named  tuple (a=a,b=b,root=root,d=d) of  integers such that   c=(a + b ER(root))//d or nothing if no such tuple exists

**Examples**

```julia-repl
julia> quadratic(1+E(3))
(a = 1, b = 1, root = -3, den = 2)

julia> quadratic(1+E(5))

```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Cycs.jl#L511-L524' class='documenter-source'>source</a><br>


<a id='Pols.jl-Documentation-1'></a>

# Pols.jl Documentation

<a id='Gapjm.Pols' href='#Gapjm.Pols'>#</a>
**`Gapjm.Pols`** &mdash; *Module*.



An implementation of univariate Laurent polynomials.  A Pol contains two fields: its vector of coefficients, and its valuation.

**Examples**

```julia-repl
julia> Pol(:q) # define string used for printing and set variable q
q

julia> Pol([1,2],0) # coefficients should have no leading or trailing zeroes.
2q+1

julia> p=Pol([1,2],-1)
2+q⁻¹

julia> valuation(p)
-1

julia> p=(q+1)^2
q²+2q+1

julia> degree(p)
2

julia> p(1//2) # a Pol is a callable object, where the call evaluates the Pol
9//4

julia> divrem(q^3+1,q+2) # changes coefficients to field elements
(1.0q²-2.0q+4.0, -7.0)

julia> divrem1(q^3+1,q+2) # keeps the ring, but needs second argument unitary
(q²-2q+4, -7)

julia> cyclotomic_polynomial(24) # the 24-th cyclotomic polynomial
q⁸-q⁴+1

```

see also the individual documentation of gcd.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Pols.jl#L1-L40' class='documenter-source'>source</a><br>

<a id='Base.divrem' href='#Base.divrem'>#</a>
**`Base.divrem`** &mdash; *Function*.



computes (p,q) such that a=p*b+q


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Pols.jl#L168-L170' class='documenter-source'>source</a><br>

<a id='Gapjm.Pols.divrem1' href='#Gapjm.Pols.divrem1'>#</a>
**`Gapjm.Pols.divrem1`** &mdash; *Function*.



divrem when b unitary: does not change type


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Pols.jl#L186-L188' class='documenter-source'>source</a><br>

<a id='Base.gcd' href='#Base.gcd'>#</a>
**`Base.gcd`** &mdash; *Function*.



gcd(p::Pol, q::Pol)   the coefficients of p and q must be elements of a field for   gcd to be type-stable

**Examples**

```julia-repl
julia> gcd(q+1,q^2-1)
1.0q+1.0

julia> gcd(q+1//1,q^2-1//1)
(1//1)q+1//1
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Pols.jl#L208-L221' class='documenter-source'>source</a><br>


<a id='CoxGroups.jl-Documentation-1'></a>

# CoxGroups.jl Documentation

<a id='Gapjm.CoxGroups' href='#Gapjm.CoxGroups'>#</a>
**`Gapjm.CoxGroups`** &mdash; *Module*.



A  suitable  reference  for  the  general  theory of Coxeter groups is, for example, Bourbaki "Lie Groups and Lie Algebras" chapter 4.

A *Coxeter group* is a group which has the presentation `W=⟨S|(st)^m(s,t)=1`  for  `s,t∈  S⟩`  for  some  symmetric  integer matrix `m(s,t)`  called  the  *Coxeter  matrix*,  where  `m(s,t)>1`  for `s≠t` and `m(s,s)=1`.  It is true (but a non-trivial theorem) that in a Coxeter group the  order of `st` is exactly `m(s,t)`, thus a Coxeter group is the same as a  *Coxeter system*, that is a pair `(W,S)` of a group `W` and a set `S` of involutions,  such that the group is  presented by relations describing the order  of the product of two elements of `S`. A Coxeter group has a natural representation, its *reflection representation*, on a real vector space `V` of  dimension `length(S)` (the *Coxeter rank*  of W), where each element of `S`  acts as a  reflection; the faithfulness  of this representation in the main  argument to prove  that the order  of `st` is  exactly `m(s,t)`. Thus Coxeter groups are real reflection groups. The converse need not be true if the  set of reflecting  hyperplanes has bad  topological properties, but it turns out that finite Coxeter groups are the same as finite real reflection groups.  The possible Coxeter matrices for  finite Coxeter groups have been completely  classified; the corresponding finite groups play a deep role in several areas of mathematics.

Coxeter  groups  have  a  nice  solution  to the word problem. The *length* `l(w)`  of an element  `w∈ W` is  the minimum number  of elements of `S` of which it is a product (since the elements of `S` are involutions, we do not need inverses). An expression of `w` of minimal length is called a *reduced word*  for `w`. The main property of  reduced words is the *exchange lemma* which  states that if `s₁…sₖ` is a  reduced word for `w` (thus`k=l(w)`) and `s∈  S` is such that `l(sw)≤l(w)` then one  of the `sᵢ` in the word for `w` can be deleted to obtain a reduced word for `sw`. Thus given `s∈ S` and `w∈ W`,  either `l(sw)=l(w)+1` or  `l(sw)=l(w)-1` and we  say in this last case that  `s` belongs to  the *left descent  set* of `w`.  The computation of a reduced word for an element, and other word problems, are easily done if we know  the left descent sets. For the Coxeter groups that we implement, this left  descent set  can be  easily determined  (see e.g. 'coxsym' below), so this suggests how to deal with Coxeter groups.

The type `CoxeterGroup` is an abstact type; an actual struct which implements it must define a function

`isleftdescent(W,w,i)` which tells whether the       `i`-th element of `S` is in the left descending set of `w`.

the other functions needed in an instance of a Coxeter group are

  * `coxgens(W)` which returns the set `S` (the list of *Coxeter generators*)
  * `nref(W)` which  returns the  number of  reflections of  `W`, if  `W` is  finite or `nothing` if `W` is infinite

It  should be  noted that  a Coxeter group can be *any* kind of group implementing the above functions.

A  common occurrence in code for Coxeter groups is a loop like:

`findfirst(eachindex(coxgens(W)),x->isleftdescent(W,w,x))`

if you provide a function `firstleftdescent(W,w)` it will be called instead of the above loop.

Because  of the  easy solution  of the  word problem  in Coxeter  groups, a convenient  way  to  represent  their  elements  is as words in the Coxeter generators.  They are represented as lists of labels for the generators. By default  these labels are  given as the  index of a  generator in `S`, so a Coxeter  word is just  a list of  integers in `1:length(S)`. For reflection subgroups, the labels are indices of the reflections in the parent group.

The functions 'word' and 'element' will do the conversion between Coxeter words and elements of the group.

**Examples**

```julia-repl
julia> W=coxsym(4)
coxsym(4)

julia> p=element(W,[1,3,2,1,3])
{UInt8}(1,4)

julia> word(W,p)
5-element Array{Int64,1}:
 1
 2
 3
 2
 1

```

We  notice that the word we started with and the one that we ended up with, are not the same, though they represent the same element of `W`. The reason is  that the function 'word' computes a lexicographically smallest word for `w`.  Below  are  some  other  possible  computations with the same Coxeter group:

```julia-repl
julia> word(W,longest(W))  # the (unique) longest element in W
6-element Array{Int64,1}:
 1
 2
 1
 3
 2
 1

julia> w0=longest(W)
{UInt8}(1,4)(2,3)
julia> length(W,w0)
6
julia> map(i->word(W,reflection(W,i)),1:nref(W))
6-element Array{Array{Int64,1},1}:
 [1]            
 [2]            
 [3]            
 [1, 2, 1]      
 [2, 3, 2]      
 [1, 2, 3, 2, 1]
julia> [length(elements(W,i)) for i in 0:nref(W)]
7-element Array{Int64,1}:
 1
 3
 5
 6
 5
 3
 1

```

The above line tells us that there is 1 element of length 0, there are 6 of length 3, …

For  most basic functions the convention is that the input is an element of the  group, rather than  a Coxeter word.  The reason is  that for a Coxeter group  which  is  a  permutation  group,  using the low level functions for permutations  is usually  much faster  than manipulating lists representing reduced expressions.

This  file contains mostly a port of  the basic functions on Coxeter groups in  CHEVIE. The only Coxeter group  constructor implemented here is coxsym. The file Weyl.jl defines WeylGroup.

The dictionary from CHEVIE is as follows:

```
     CoxeterElements(W[,l])                → elements(W[,l])
     CoxeterLength(W,w)                    → length(W,w)
     CoxeterWord(W,w)                      → word(W,w)
     LongestCoxeterElement(W)              → longest(W)
     FirstLeftDescending(W,w)              → firstleftdescent(W,w)
     ReducedInRightCoset(W,w)              → reduced(W,w)
     ReducedRightCosetRepresentatives(W,H) → reduced(H,W)
     SemiSimpleRank(W)                     → coxrank(W)
     CoxeterGroupSymmetricGroup(n)         → coxsym(n)
     ReflectionSubgroup                    only standard parabolics now
     IsLeftDescending(W,w,i)               → isleftdescent(W,w,i)
     ReflectionDegrees(W)                  → degrees(W)
     ReflectionLength(W,w)                 → reflength(W,w)
     W.N                                   → nref(W)
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/CoxGroups.jl#L1-L157' class='documenter-source'>source</a><br>


<a id='Weyl.jl-Documentation-1'></a>

# Weyl.jl Documentation

<a id='Gapjm.Weyl' href='#Gapjm.Weyl'>#</a>
**`Gapjm.Weyl`** &mdash; *Module*.



*Weyl groups* are the finite real reflection groups (finite Coxeter groups) which  can be defined  over the rational  numbers (that is, finite rational reflection  groups). We implement them as a subgroup of the permutations of a  root system. Root systems play an  important role in mathematics as they classify semi-simple Lie algebras and algebraic groups.

Let  us give precise definitions. Let `V`  be a real vector space, `Vⱽ` its dual  and let `(,)`  be the natural  pairing between `Vⱽ`  and `V`. A *root system*  is a finite set  of vectors `R` which  generate `V` (the *roots*), together  with  a  map  `r↦  rⱽ`  from  `R`  to  a subset `Rⱽ` of `Vⱽ` (the *coroots*) such that:

  * For any `r∈  R`, we have  `(rⱽ,r)=2` so that  the formula `x↦ x-(rⱽ,x)r`

defines a reflection `s_r:V→ V` with root `r` and coroot `rⱽ`.

  * The reflection `s_r` stabilizes `R`.

We  will only  consider *reduced*  root systems,  i.e., such  that the only elements of `R` colinear with `r∈ R` are `r` and `-r`; we also ask that the root  system be *crystallographic*, that is `(rⱽ,s)` is an integer, for any `s∈ R,rⱽ∈ Rⱽ`.

The  subgroup `W=W(R)` of  `GL(V)` generated by  the reflections `s_r` is a Weyl  group –- the representation `V` of  `W` is defined over the rational numbers.  All other finite-dimensional (complex)  representations of a Weyl group  can also be realized  over the rational numbers.  Weyl groups can be characterized  amongst finite Coxeter  groups by the  fact that all numbers `m(s,t)` in the Coxeter matrix are in `{2,3,4,6}`.

If  we identify  `V` with  `Vⱽ` by  choosing a  `W`-invariant bilinear form `(.;.)`;  then we have `rⱽ=2r/(r;r)`. A root system `R` is *irreducible* if it is not the union of two orthogonal subsets. If `R` is reducible then the corresponding  Weyl  group  is  the  direct  product  of  the  Weyl  groups associated   with  the  irreducible  components  of  `R`.  The  irreducible crystallographic  root  systems  are  classified  by  the following list of *Dynkin  diagrams*. We show the labeling of the nodes given by the function 'diagram' described below.

```
A_n O—O—O—…—O   B_n O⇐O—O—…—O  C_n O⇒O—O—…—O  D_n  O 2
    1 2 3 … n       1 2 3 … n      1 2 3 … n       ￨
                                                 O—O—…—O
                                                 1 3 … n

G_2 O⇛O   F_4 O—O⇒O—O  E_6   O 2   E_7   O 2     E_8    O 2
    1 2       1 2 3 4        ￨           ￨              ￨
                         O—O—O—O—O   O—O—O—O—O—O    O—O—O—O—O—O—O
                         1 3 4 5 6   1 3 4 5 6 7    1 3 4 5 6 7 8
```

These diagrams encode the presentation of the Coxeter group `W` as follows: the vertices represent the generating reflections; an edge is drawn between `s`  and `t` if the order `m(s,t)` of `st` is greater than `2`; the edge is single  if  `m(s,t)=3`,  double  if  `m(s,t)=4`,  triple if `m(s,t)=6`. The arrows  indicate the relative root lengths when `W` has more than one orbit on  `R`, as explained below; we  get the *Coxeter Diagram*, which describes the  underlying Weyl group, if  we ignore the arrows:  we see that the root systems `B_n` and `C_n` correspond to the same Coxeter group.

Let us now describe how the root systems are encoded in these diagrams. Let `R`  be a root system in `V`. Then we can choose a linear form on `V` which vanishes  on no element of `R`. According to  the sign of the value of this linear  form on a root  `r ∈ R` we  call `r` *positive* or *negative*. Then there  exists a unique subset `Π` of  the positive roots, called the set of *simple  roots*, such that  any positive root  is a linear combination with non-negative  coefficients of roots  in `Π`. It  can be shown  that any two sets of simple roots (corresponding to different choices of linear forms as above)  can be transformed into  each other by a  unique element of `W(R)`. Hence, since the pairing between `V` and `Vⱽ` is `W`-invariant, if `Π` is a set  of simple roots and if we define  the *Cartan matrix* as being the `n` times  `n` matrix `C={rⱽ(s)}_{rs}`, for `r,s∈Π` this matrix is unique up to simultaneous  permutation of rows and columns.  It is precisely this matrix which is encoded in a Dynkin diagram, as follows.

The  indices for the rows of `C` label the nodes of the diagram. The edges, for  `r ≠ s`, are  given as follows. If  `C_{rs}` and `C_{sr}` are integers such  that `|C_{rs}| ≥  |C_{sr}|` the vertices  are connected by `|C_{rs}|` lines,  and if `|C_{rs}|>1`  then we put  an additional arrow  on the lines pointing towards the node with label `s`. In all other cases, we simply put a  single line  equipped with  the unique  integer `p_{rs}  ≥ 1`  such that `C_{rs}C_{sr}=cos^2 (π/p_{sr})`.

It  is now important to note that, conversely, the whole root system can be recovered from the simple roots. The reflections in `W(R)` corresponding to the  simple roots are  called *simple* reflections.  They are precisely the generators  for which the Coxeter diagram encodes the defining relations of `W(R)`. Each root is in the orbit of a simple root, so that `R` is obtained as  the orbit of the  simple roots under the  group generated by the simple reflections.  The  restriction  of  the  simple  reflections  to  `V_R`  is determined  by the Cartan matrix, so `R` is determined by the Cartan matrix and the set of simple roots.

The  Cartan  matrix  corresponding  to  one  of  the above irreducible root systems  (with the specified labeling) is  returned by the command 'cartan' which takes as input a symbol giving the type (that is ':A', ':B', …, ':I') and a positive integer giving the rank. This function returns a matrix with entries  in `bbZ`.  Given two  Cartan matrices  `c1` and `c2`, their matrix direct sum (corresponding to the orthogonal direct sum of the root systems) can be produced by the `cat(c1,c2,dims=[1,2])`.

The function 'WeylGroup' takes as input some data which determine the roots and  the coroots and  produces a `struct`  containing information about the root  system `R` and about `W(R)`. If we label the positive roots by '1:N', and  the  negative  roots  by  'N+1:2N',  then  each  simple  reflection is represented by the permutation of '1:2N' which it induces on the roots.

The function 'WeylGroup' has 2 methods; in one of them, the argument is the Cartan  matrix of the root system. This  constructs a root system where the simple  roots are the canonical basis of `V`, and the matrix of the coroots expressed in the dual basis of `Vⱽ` is then equal to the Cartan matrix.

If one only wants to work with Cartan matrices with a labeling as specified by  the  above  list,  the  function  call  can  be  simplified. Instead of 'WeylGroup(CartanMat(:D,4))' the following is also possible.

```
julia> W=WeylGroup(:D,4)
WeylGroup(:D,4)

julia> W.cartan
4×4 Array{Int8,2}:
  2   0  -1   0
  0   2  -1   0
 -1  -1   2  -1
  0   0  -1   2
```

Also,  the Weyl group struct associated to a direct sum of irreducible root systems can be obtained as a product

```
julia> W=WeylGroup(:A,2)*WeylGroup(:B,2)
WeylGroup(:A,<1,2>)xWeylGroup(:B,<3,4>)

julia> W.cartan
4×4 Array{Int8,2}:
  2  -1   0   0
 -1   2   0   0
  0   0   2  -2
  0   0  -1   2
```

The  same `struct`  is constructed  by applying  'WeylGroup' to  the matrix 'cat(cartan(:A,2), cartan(:B,2),dims=[1,2])'.

The elements of a Weyl group are permutations of the roots:

```julia-repl
julia> W=WeylGroup(:D,4)
WeylGroup(:D,4)

julia> p=element(W,[1,3,2,1,3])
{UInt8}(1,14,13,2)(3,17,8,18)(4,12)(5,20,6,15)(7,10,11,9)(16,24)(19,22,23,21)

julia> word(W,p)
5-element Array{Int64,1}:
 1
 3
 1
 2
 3

```

This module is mostly a port of the basic functions on Weyl groups in CHEVIE. The dictionary from CHEVIE is as follows:

```
     CartanMat("A",5)                       →  cartan(:A,5) 
     CoxeterGroup("A",5)                    →  WeylGroup(:A,5) 
     Size(W)                                →  length(W) 
     ForEachElement(W,f)                    →  for w in W f(w) end 
     ReflectionDegrees(W)                   →  degrees(W) 
     IsLeftDescending(W,w,i)                →  isleftdescent(W,w,i) 
     ReflectionSubgroup : only standard parabolics now
     TwoTree(m)                             →  twotree(m) 
     FiniteCoxeterTypeFromCartanMat(m)      →  type_cartan(m) 
     RootsCartan(m)                         →  roots(m) 
     PrintDiagram(W)                        →  diagram(W) 
     Inversions                             →  inversions 
     Reflection                             →  reflection 
     W.orbitRepresentative[i]               →  simple_representative(W,i) 
```

finally, a benchmark on julia 1.0.2

```benchmark
julia> @btime length(elements(WeylGroup(:E,7)))
  807.009 ms (8608502 allocations: 787.93 MiB)
```

GAP3 for the same computation takes 2.2s


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Weyl.jl#L1-L187' class='documenter-source'>source</a><br>

<a id='Gapjm.Weyl.cartan' href='#Gapjm.Weyl.cartan'>#</a>
**`Gapjm.Weyl.cartan`** &mdash; *Function*.



```
`cartan(type, rank)`
```

Cartan matrix for a Weyl group:

```julia-repl
julia> cartan(:A,4)
4×4 Array{Int8,2}:
  2  -1   0   0
 -1   2  -1   0
  0  -1   2  -1
  0   0  -1   2
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Weyl.jl#L195-L208' class='documenter-source'>source</a><br>

<a id='Gapjm.Weyl.two_tree' href='#Gapjm.Weyl.two_tree'>#</a>
**`Gapjm.Weyl.two_tree`** &mdash; *Function*.



```
two_tree(m)
```

Given  a square  matrix m  with zeroes  (or falses,  for a boolean matrix)  symmetric  with respect to the diagonal, let  G be the graph with vertices  axes(m)[1] and an edge between i and j iff !iszero(m[i,j]).  If G  is a line this function returns it as a Vector{Int}.   If  G  is  a  tree  with  one  vertex  c of valence 3 the function returns  (c,b1,b2,b3)  where b1,b2,b3 are  the branches from  this vertex sorted by  increasing length.  Otherwise the function returns `nothing`

```julia-repl
julia> CoxGroups.two_tree(cartan(:A,4))
4-element Array{Int64,1}:
 1
 2
 3
 4

julia> CoxGroups.two_tree(cartan(:E,8))
(4, [2], [3, 1], [5, 6, 7, 8])
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Weyl.jl#L222-L244' class='documenter-source'>source</a><br>


#type_cartan #roots


<a id='Hecke.jl-Documentation-1'></a>

# Hecke.jl Documentation

<a id='Gapjm.Hecke' href='#Gapjm.Hecke'>#</a>
**`Gapjm.Hecke`** &mdash; *Module*.



This   module  ports   Chevie  functionality   for  Iwahori-Hecke  algebras associated to Coxeter groups.

Let  (W,S) be a Coxeter  system where `mₛₜ` is  the order of `st` for `s,t∈ S`. Let `R` be a commutative ring with 1 and for `s∈ S` let `uₛ₀,uₛ₁∈ R` be elements which depend ony on the conjugacy class of `s` in `W` (this is the same  as requiring that `uₛᵢ=uₜᵢ` whenever `mₛₜ` is odd). The Iwahori-Hecke algebra of `W` over `R` with parameters `uₛᵢ` is a deformation of the group algebra  of `W` over `R` defined as  follows: it is the unitary associative `R`-algebra generated by elements `Tₛ, s∈ S` subject to the relations:

$(Tₛ-uₛ₀)(Tₛ-uₛ₁)=0$ for all `s∈ S` (the quadratic relations)

$TₛTₜTₛ…= TₜTₛTₜ…$ with `mₛₜ` factors on each side (the braid relations)

If  `uₛ₀=1` and  `uₛ₁=-1` for  all `s`  then the quadratic relations become `Tₛ²=1` and the deformation of the group algebra is trivial.

Since  the generators `Tₛ` satisfy the  braid relations, the algebra `H` is in  fact a quotient of the group algebra of the braid group associated with `W`.  It follows that, if `w=s_1⋯ s_m`  is a reduced expression of `w ∈ W` then  the  product  `Tₛ_1⋯ Tₛ_m`  depends  only  on `w`. We will therefore denote by `T_w`. We have `T_1=1`.

If  one of `uₛ₀` or `uₛ₁` is invertible  in `R`, for example `uₛ₁`, then by changing  the generators  to `T′ₛ=-Tₛ/uₛ₁`,  and setting `qₛ=-uₛ₀/uₛ₁`, the braid  relations do no change  (since when `mₛₜ` is  odd we have `uₛᵢ=uₜᵢ`) but  the quadratic relations become  `(T′ₛ-qₛ)(T′ₛ+1)=0`. This last form is the  most common  form considered  in the  literature. Another common form, considered  in  the  context  of  Kazhdan-Lusztig  theory, is `uₛ₀=√qₛ` and `uₛ₁=-√qₛ⁻¹`.  The general form of parameters provided is a special case of general cyclotomic Hecke algebras, and can be useful in many contexts.

For  some  algebras  the  character  table,  and in general Kazhdan-Lusztig bases,  require a square root of `qₛ=-uₛ₀/uₛ₁`. We provide a way to specify it  with  the  field  `.sqpara`  which  can  be given when constructing the algebra. If not given a root is automatically extracted when needed (and we know  how to compute it) by the function `RootParameter`. Note however that sometimes  an  explicit  choice  of  root  is  necessary  which  cannot  be automatically determined.

There  is a universal choice  for `R` and `uₛᵢ`:  Let `uₛᵢ:s∈ S,i∈[0,1]` be indeterminates   such  that  `uₛᵢ=uₜᵢ`  whenever  `mₛₜ`  is  odd,  and  let `A=ℤ[uₛᵢ]` be the corresponding polynomial ring. Then the Hecke algebra `H` of  `W` over a  with parameters `uₛᵢ`  is called the *generic Iwahori-Hecke algebra*  of  with  `W`.  Any  other  algebra  with parameters `vₛᵢ` can be obtained  by specialization from  `H`: There is  a unique ring homomorphism `f:A  → R` such that `f(uₛᵢ)=vₛᵢ`  for all `i`. Then we  can view `R` as an `A`-module via `f` and we can identify the other algebra to $R⊗ _A H$.

The  elements `{T_w∣w∈ W}` actually form an  `R`-basis of `H` if one of the `uₛᵢ`  is invertible for all `s`. The  structure constants in that basis is obtained  as  follows.  To  multiply  `T_v`  by  `T_w`,  choose  a  reduced expression for `v`, say `v=s_1 ⋯ s_k` and apply inductively the formula:

$T_sT_w=T_{sw}$               if `l(sw)=l(w)+1`

$T_sT_w=-uₛ₀uₛ₁T_{sw}+(uₛ₀+uₛ₁)T_w$ if `l(sw)=l(w)-1`.

If all `s` we have `uₛ₀=q`, `uₛ₁=-1` then we call the corresponding algebra the one-parameter or Spetsial Iwahori-Hecke algebra associated with `W`; it can  be obtained with the  simplified call 'Hecke(W,q)'. Certain invariants of  the irreducible characters of  this algebra play a  special role in the representation  theory of the underlying  finite Coxeter groups, namely the `a`- and `A`-invariants. For basic properties of Iwahori-Hecke algebras and their  relevance to the representation theory of finite groups of Lie type, see for example Curtis and Reiner 1987, Sections~67 and 68.

In  the  following  example,  we  compute  the multiplication table for the `0`-Iwahori–Hecke algebra associated with the Coxeter group of type `A_2`.

```julia-repl
julia> W=WeylGroup(:A,2)
WeylGroup(:A,2)

julia> H=hecke(W,0)             # One-parameter algebra with `q=0`
Hecke(WeylGroup(:A,2),0)

julia> T=Tbasis(H)              # Create the `T` basis
(::getfield(Gapjm.Hecke, Symbol("#f#20")){Int64,Perm{UInt8},HeckeAlgebra{Int64,WeylGroup}}) (generic function with 4 methods)

julia> el=words(W)
6-element Array{Array{Int64,1},1}:
 []       
 [2]      
 [1]      
 [2, 1]   
 [1, 2]   
 [1, 2, 1]

julia> T.(el)*permutedims(T.(el))        # multiplication table
6×6 Array{HeckeTElt{Perm{UInt8},Int64,WeylGroup},2}:
 T.    T₂     T₁     T₂₁    T₁₂    T₁₂₁ 
 T₂    -T₂    T₂₁    -T₂₁   T₁₂₁   -T₁₂₁
 T₁    T₁₂    -T₁    T₁₂₁   -T₁₂   -T₁₂₁
 T₂₁   T₁₂₁   -T₂₁   -T₁₂₁  -T₁₂₁  T₁₂₁ 
 T₁₂   -T₁₂   T₁₂₁   -T₁₂₁  -T₁₂₁  T₁₂₁ 
 T₁₂₁  -T₁₂₁  -T₁₂₁  T₁₂₁   T₁₂₁   -T₁₂₁

```

Thus,  we work  with algebras  with arbitrary  parameters. We will see that this also works on the level of characters and representations.

finally, benchmarks on julia 1.0.2

```benchmark
julia> function test_w0(n)
         W=WeylGroup(:A,n)
         Tbasis(hecke(W,Pol([1],1)))(longest(W))^2
       end
test_w0 (generic function with 1 method)

julia> @btime test_w0(7);
  132.737 ms (178853 allocations: 157.37 MiB)
```

Compare to GAP3 where the following function takes 0.92s

```
test_w0:=function(n)local W,T,H;
  W:=CoxeterGroup("A",n);H:=Hecke(W,X(Rationals));T:=Basis(H,"T");
  T(LongestCoxeterWord(W))^2;
end;
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Hecke.jl#L1-L123' class='documenter-source'>source</a><br>

<a id='Gapjm.Hecke.HeckeAlgebra' href='#Gapjm.Hecke.HeckeAlgebra'>#</a>
**`Gapjm.Hecke.HeckeAlgebra`** &mdash; *Type*.



hecke( W [, parameter, [rootparameter]] ) return a Hecke algebra for W

**Example**

```julia-repl
julia> W=WeylGroup(:B,2)
WeylGroup(:B,2)

julia> Pol(:q)
q

julia> H=hecke(W,q)
Hecke(WeylGroup(:B,2),q)

julia> [H.para,H.sqpara]
2-element Array{Array{T,1} where T,1}:
 Tuple{Pol{Int64},Pol{Int64}}[(q, -1), (q, -1)]
 Missing[missing, missing]                     

julia> H=hecke(W,q^2,q)
Hecke(WeylGroup(:B,2),q²,q)

julia> [H.para,H.sqpara]
2-element Array{Array{T,1} where T,1}:
 Tuple{Pol{Int64},Pol{Int64}}[(q², -1), (q², -1)]
 Pol{Int64}[q, q]                                  

julia> H=hecke(W,[q^2,q^4],[q,q^2])
Hecke(WeylGroup(:B,2),Pol{Int64}[q², q⁴],Pol{Int64}[q, q²])

julia> [H.para,H.sqpara]
2-element Array{Array{T,1} where T,1}:
 Tuple{Pol{Int64},Pol{Int64}}[(q², -1), (q⁴, -1)]
 Pol{Int64}[q, q²]

julia> H=hecke(W,9,3)
Hecke(WeylGroup(:B,2),9,3)

julia> [H.para,H.sqpara]
2-element Array{Array{T,1} where T,1}:
 Tuple{Int64,Int64}[(9, -1), (9, -1)]
 [3, 3]                              
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Hecke.jl#L135-L178' class='documenter-source'>source</a><br>


<a id='KL.jl-Documentation-1'></a>

# KL.jl Documentation

<a id='Gapjm.KL' href='#Gapjm.KL'>#</a>
**`Gapjm.KL`** &mdash; *Module*.



This  module ports Chevie functionality for Kazhdan-Lusztig polynomials and bases.

Let  `ℋ` be  the Iwahori-Hecke algebra  of a Coxeter  system `(W,S)`, with quadratic  relations `(Tₛ-uₛ₀)(Tₛ-uₛ₁)=0` for  `s∈ S`. If  `-uₛ₀/uₛ₁` has a square  root  `vₛ`,  we  can  scale  the  basis  `Tₛ`  to  get  a new basis `tₛ=Tₛ/(-vₛuₛ₁)`  with quadratic  relations `(tₛ-vₛ)(tₛ+vₛ⁻¹)=0`.  The most general  case when Kazhdan-Lusztig bases and  polynomials can be defined is when  the parameters `vₛ` belong to a totally ordered abelian group `Γ` for multiplication,  see Lus83. We set `Γ⁺= {γ∈ Γ∣γ>0}` and `Γ⁻={γ⁻¹∣γ∈ Γ⁺}={γ∈ Γ∣γ<0}`.

Thus  we assume `ℋ` defined over the ring `ℤ[Γ]`, the group algebra of `Γ` over  `ℤ`, and the quadratic  relations of `ℋ`  associate to each `s∈ S` a `vₛ∈  Γ⁺` such that  `(tₛ-vₛ)(tₛ+vₛ⁻¹)=0`. We also  set `qₛ=vₛ²` and define the  basis `Tₛ=vₛtₛ` with quadratic relations `(Tₛ-qₛ)(Tₛ+1)=0`; for `w∈ W` with reduced expression `w=s₁…sₙ` we define `q_w∈ Γ⁺` by `q_w^½=v_{s₁}…v_{sₙ}` and let `q_w=(q_w^½)²`.

We  define the bar involution on `ℋ`  by linearity: on `ℤ[Γ]` we define it by  $\overline{∑_{γ∈ Γ}a_γγ}= ∑_{γ∈ Γ} a_γ γ⁻¹$ and we extend it to `ℋ` by  $\overline  Tₛ=Tₛ⁻¹$.  Then  the  Kazhdan-Lusztig  basis `C′_w` is defined  as  the  only  basis  of  `ℋ`  stable  by the bar involution and congruent to `t_w` modulo `∑_{w∈ W}Γ⁻ t_w`.

The  basis `C′_w` can be computed  as follows. We define elements `R_{x,y}` of  `ℤ[Γ]` by  `T_y⁻¹=∑_x \overline{R_{x,y⁻¹}}  q_x⁻¹T_x`. We  then define inductively  the Kazhdan-Lusztig  polynomials (in  this general  context we should  say the  Kazhdan-Lusztig elements  of `ℤ[Γ]`,  which belong  to the subalgebra  of `ℤ[Γ]` generated by  the `qₛ`) by $P_{x,w}=τ_{≤(q_w/q_x)^½} (∑_{x<y≤w}R_{x,y}P_{y,w})$  where `τ`  is the  truncation: $τ_≤ν ∑_{γ∈ Γ} a_γγ= ∑_{γ≤ν}a_γγ$; the induction is thus on decreasing `x` for the Bruhat order  and  starts  at  `P_{w,w}=1`.  We  have  then  $C′_w=∑_y q_w^{-1/2} P_{y,w}T_y$.

The  Chevie code  for the  Kazhdan-Lusztig bases  `C`, `D` and their primed versions, has been initially written by Andrew Mathas around 1994, who also contributed  to  the  design  of  the programs dealing with Kazhdan-Lusztig bases. He also implemented some other bases, such as the Murphy basis which can  be  found  in  the  Chevie  contributions  directory. The code for the unequal  parameters  case  has  been  written  around  1999  by F.Digne and J.Michel. The other Kazhdan-Lusztig bases are computed in terms of the `C′` basis.

When  the `ℤ[Γ]` is a  Laurent polynomial ring the  bar operation is taking the  inverse of  the variables,  and truncation  is keeping terms of degree smaller or equal to that of `ν`. It is possible to use arbitrary groups `Γ` as  long as  methods `bar`,  `positive_part` and  `negative_part` have been defined on them. which perform the operations respectively $∑_{γ∈ Γ} a_γγ↦ ∑_{γ∈  Γ} a_γγ⁻¹$, $∑_{γ∈ Γ}  a_γγ↦ ∑_{γ≥ 1} a_γγ$  and $∑_{γ∈ Γ} a_γγ↦ ∑_{γ≤  1} a_γγ$ on  elements 'p' of  `ℤ[Γ]`. The operations  above will be used internally by the programs to compute Kazhdan-Lusztig bases.

finally, benchmarks on julia 1.0.2

```benchmark
julia> function test_kl(W)
         q=Pol([1],1); H=hecke(W,q^2,q)
         C=Cpbasis(H); T=Tbasis(H)
         [T(C(w)) for w in elements(W)]
       end
test_kl (generic function with 1 method)

julia> @btime test_kl(WeylGroup(:F,4));
2.265 s (22516606 allocations: 1.81 GiB)
```

Compare to GAP3 where the following function takes 11s for F4

```
test_kl:=function(W)local q,H,T,C;
  q:=X(Rationals);H:=Hecke(W,q^2,q);
  T:=Basis(H,"T");C:=Basis(H,"C'");
  List(Elements(W),e->T(C(e)));
end;
```

Another benchmark:

```benchmark
function test_kl2(W)
  el=elements(W)
  [KLPol(W,x,y) for x in el, y in el]
end

test_kl2 (generic function with 1 method)

julia>@btime test_kl2(WeylGroup(:F,4));
  8s (97455915 allocations: 6.79 GiB)
```

Compare to GAP3 where the following function takes 42s for F4

```
test_kl2:=function(W)local el;
  el:=Elements(W);
  List(el,x->List(el,y->KazhdanLusztigPolynomial(W,x,y)));
end;
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/KL.jl#L22-L115' class='documenter-source'>source</a><br>

<a id='Gapjm.KL.KLPol' href='#Gapjm.KL.KLPol'>#</a>
**`Gapjm.KL.KLPol`** &mdash; *Function*.



KLPol(W,y,w) returns the Kazhdan-Lusztig polynomial P_{y,w} of W

To  compute Kazhdan-Lusztig polynomials in  the one-parameter case it seems that  the best  approach still  is by  using the  recursion formula  in the original  article KL79. One can first run  a number of standard checks on a given  pair  of  elements  to  see  if the computation of the corresponding polynomial  can be reduced to a similar computation for elements of smaller length. One such check involves the notion of critical pairs (cf. Alv87): a pair  of elements `w₁,w₂∈  W` such that  `w₁≤w₂` is *critical*  if `ℒ(w₂) ⊆ ℒ(w₁)`  and `ℛ (w₂)⊆ ℛ (w₁)`, where `ℒ`  and `ℛ` denote the left and right descent  set, respectively.  Now if  `y≤w ∈  W` are arbitrary elements then there   always  exists  a  critical  pair   `z≤w`  with  `y≤z≤w`  and  then `P_{y,w}=P_{z,w}`.  Given two elements `y` and `w`, such a critical pair is found by the function 'CriticalPair'. Whenever the polynomial corresponding to a critical pair is computed then this pair and the polynomial are stored in the property `:klpol` of the underlying Coxeter group.

```julia-repl
julia> W=WeylGroup(:B,3)
WeylGroup(:B,3)

julia> map(i->map(x->KLPol(W,one(W),x),elements(W,i)),1:W.N)
9-element Array{Array{Pol{Int64},1},1}:
 [1, 1, 1]                       
 [1, 1, 1, 1, 1]                 
 [1, 1, 1, 1, 1, 1, 1]           
 [1, 1, 1, x+1, 1, 1, 1, 1]      
 [1, x+1, 1, x+1, x+1, 1, 1, x+1]
 [1, x²+1, x+1, x+1, 1, 1, x+1]  
 [1, 1, x²+x+1, x+1, x+1]        
 [1, x+1, x²+1]                  
 [1]
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/KL.jl#L186-L220' class='documenter-source'>source</a><br>

<a id='Gapjm.Hecke.Tbasis' href='#Gapjm.Hecke.Tbasis'>#</a>
**`Gapjm.Hecke.Tbasis`** &mdash; *Function*.



```julia-repl
julia> W=WeylGroup(:B,3)
WeylGroup(:B,3)

julia> Pol(:v);H=hecke(W,v^2,v)
Hecke(WeylGroup(:B,3),v²,v)

julia> C=Cpbasis(H)
(::getfield(Gapjm.KL, Symbol("#f#10")){Pol{Int64},Perm{UInt8},HeckeAlgebra{Pol{Int64},WeylGroup}}) (generic function with 3 methods)

julia> T=Tbasis(H)
(::getfield(Gapjm.Hecke, Symbol("#f#20")){Pol{Int64},Perm{UInt8},HeckeAlgebra{Pol{Int64},WeylGroup}}) (generic function with 4 methods)

julia> T(C(1,2))
v⁻²T.+v⁻²T₂+v⁻²T₁+v⁻²T₁₂
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/KL.jl#L309-L326' class='documenter-source'>source</a><br>


<a id='Util.jl-Documentation-1'></a>

# Util.jl Documentation

<a id='Gapjm.Util' href='#Gapjm.Util'>#</a>
**`Gapjm.Util`** &mdash; *Module*.



This  module contains  various utility  functions used  in the  rest of the code.  Maybe some  of them  exist in  some Julia  module I am not aware of; please tell me.

The code is divided in sections  according to semantics.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L1-L7' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.cartesian' href='#Gapjm.Util.cartesian'>#</a>
**`Gapjm.Util.cartesian`** &mdash; *Function*.



Cartesian product of list of n vectors, returned as an x-by-n matrix


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L40' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.groupby' href='#Gapjm.Util.groupby'>#</a>
**`Gapjm.Util.groupby`** &mdash; *Function*.



group items of list l according to the corresponding values in list v

```
julia> groupby([31,28,31,30,31,30,31,31,30,31,30,31],
       [:Jan,:Feb,:Mar,:Apr,:May,:Jun,:Jul,:Aug,:Sep,:Oct,:Nov,:Dec])
Dict{Int64,Array{Symbol,1}} with 3 entries:
  31 => Symbol[:Jan, :Mar, :May, :Jul, :Aug, :Oct, :Dec]
  28 => Symbol[:Feb]
  30 => Symbol[:Apr, :Jun, :Sep, :Nov]
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L46-L56' class='documenter-source'>source</a><br>


group items of list l according to the values taken by function f on them

```
julia> groupby(iseven,1:10)
Dict{Bool,Array{Int64,1}} with 2 entries:
  false => [1, 3, 5, 7, 9]
  true  => [2, 4, 6, 8, 10]
```

Note:in this version l is required to be non-empty since I do not know how to access the return type of a function


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L65-L75' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.constant' href='#Gapjm.Util.constant'>#</a>
**`Gapjm.Util.constant`** &mdash; *Function*.



whether all elements in list a are equal


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L84' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.blocks' href='#Gapjm.Util.blocks'>#</a>
**`Gapjm.Util.blocks`** &mdash; *Function*.



blocks(M::Matrix)

M  should be a square matrix. Define  a graph G with vertices 1:size(M,1)   and  with an edge between i and j  if either M[i,j] or M[j,i] is not zero   or false. blocks returns a vector of vectors I such that I[1],I[2], etc..   are  the  vertices  in  each  connected  component  of G. In other words,   M[I[1],I[1]],M[I[2],I[2]],etc... are blocks of M.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L89-L97' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.format' href='#Gapjm.Util.format'>#</a>
**`Gapjm.Util.format`** &mdash; *Function*.



format( table; options )

General routine to format a table. Used for character tables.   Options:      row*labels          Labels for rows      column*labels       Labels for columns      rows*label          Label for column of rowLabels      separators          line numbers after which to put a separator      column*repartition  display in pieces of sizes these numbers of cols      rows                show only these rows      columns             show only these columns


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L202-L215' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.prime_residues' href='#Gapjm.Util.prime_residues'>#</a>
**`Gapjm.Util.prime_residues`** &mdash; *Function*.



the numbers less than n and prime to n 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L264' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.phi' href='#Gapjm.Util.phi'>#</a>
**`Gapjm.Util.phi`** &mdash; *Function*.



the Euler function ϕ 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L283' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.primitiveroot' href='#Gapjm.Util.primitiveroot'>#</a>
**`Gapjm.Util.primitiveroot`** &mdash; *Function*.



primitiveroot(m::Integer) a primitive root mod. m,   that is it generates multiplicatively prime_residues(m).   It exists if m is of the form 4, 2p^a or p^a for p prime>2.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/57612ec167d370efbf885124c241105ee7226c8f/src/Util.jl#L289-L293' class='documenter-source'>source</a><br>

