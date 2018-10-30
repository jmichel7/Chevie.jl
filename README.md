
<a id='Gapjm.jl-Documentation-1'></a>

# Gapjm.jl Documentation

<a id='Gapjm' href='#Gapjm'>#</a>
**`Gapjm`** &mdash; *Module*.



Here are some of my efforts porting GAP code to julia. I am not even sure this is a well-formed Julia package.

It contains for now permutations and permutation groups, cyclotomic numbers and  Laurent polynomials. Coming  soon are Coxeter  groups, Hecke algebras, braid groups and Garside monoids.

Even  though the code  is often competitive  with or faster  than GAP, I am sure there are more optimisations possible. Any comments about the code and the design are welcome.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Gapjm.jl#L1-L12' class='documenter-source'>source</a><br>

- [Gapjm.jl Documentation](index.md#Gapjm.jl-Documentation-1)
- [Perms.jl Documentation](index.md#Perms.jl-Documentation-1)
- [PermGroups.jl Documentation](index.md#PermGroups.jl-Documentation-1)
- [Cycs.jl Documentation](index.md#Cycs.jl-Documentation-1)
- [Pols.jl Documentation](index.md#Pols.jl-Documentation-1)
- [Util.jl Documentation](index.md#Util.jl-Documentation-1)


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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Perms.jl#L1-L76' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Perms.jl#L184-L194' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Perms.jl#L223-L234' class='documenter-source'>source</a><br>

<a id='Base.sign' href='#Base.sign'>#</a>
**`Base.sign`** &mdash; *Function*.



sign(a::Perm) is the signature of  the permutation a


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Perms.jl#L237' class='documenter-source'>source</a><br>


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

#Elements,  appartenance test  and other  function are  computed on  G using
#Schreier-Sims theory, that is computing the following

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

finally, benchmarks

```benchmark
julia> @btime collect(symmetric_group(8));
  8.252 ms (350515 allocations: 14.17 MiB)

julia> @btime words(symmetric_group(8));
  18.209 ms (122062 allocations: 15.22 MiB)
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/PermGroups.jl#L1-L108' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.symmetric_group' href='#Gapjm.PermGroups.symmetric_group'>#</a>
**`Gapjm.PermGroups.symmetric_group`** &mdash; *Function*.



The symmetric group of degree n 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/PermGroups.jl#L238' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.orbit' href='#Gapjm.PermGroups.orbit'>#</a>
**`Gapjm.PermGroups.orbit`** &mdash; *Function*.



orbit(G,p) is the orbit of Int p under Group G


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/PermGroups.jl#L121' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.orbit_and_representative' href='#Gapjm.PermGroups.orbit_and_representative'>#</a>
**`Gapjm.PermGroups.orbit_and_representative`** &mdash; *Function*.



returns Dict x=>g for x in orbit(G,p) and g is such that x=action(p,g)


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/PermGroups.jl#L154' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.words' href='#Gapjm.PermGroups.words'>#</a>
**`Gapjm.PermGroups.words`** &mdash; *Function*.



List of minimal words in the generators elements(G) 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/PermGroups.jl#L208' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.elements' href='#Gapjm.PermGroups.elements'>#</a>
**`Gapjm.PermGroups.elements`** &mdash; *Function*.



The list of elements of G in the same order as words


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/PermGroups.jl#L213' class='documenter-source'>source</a><br>


<a id='Cycs.jl-Documentation-1'></a>

# Cycs.jl Documentation

<a id='Gapjm.Cycs' href='#Gapjm.Cycs'>#</a>
**`Gapjm.Cycs`** &mdash; *Module*.



Cyclotomic  numbers means complex numbers which are sums of rationals times roots of unity.

They are a very important feature of GAP, since entries of character tables of finite groups are cyclotomics.

They  have a normal form given by the Zumbroich basis, which allows to find the  smallest Cyclotomic field which contains a given number, and decide in particular if a cyclotomic is zero. Let ζ*n:=e^{2iπ/n}. The Zumbroich basis of Q(ζ*n) is a particular subset of 1,ζ,ζ^2,...,ζ^{n-1} which forms a basis of Q(ζ_n) with good properties.

I  ported here Christian Stump's Sage  code, which is simpler to understand than GAP's code. The reference for the algorithms is

T. Breuer, Integral bases for subfields of cyclotomic fields AAECC 8 (1997)

Contrary   to  GAP,  I  do  not  lower  automatically  numbers  after  each computation since this is too expensive with this code. GAP also converts a Cyclotomic which is rational to a Rational, a Rational which is integral to an  Int, etc... This is tremendously useful  but needs a new type of number to be added to Julia, which requires more competent people than me.

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

julia> conductor(a) # a not lowered so still in Q(ζ_4)
4

julia> conductor(lower(a)) # but now a is lowered to Q(ζ_1)=Q
1

julia> typeof(convert(Int,a))
Int64

julia> convert(Int,E(4))
ERROR: InexactError: convert(Int64, E(4))

julia> c=inv(1+E(4)) # inverses need Rationals
1//2-1//2ζ₄

julia> typeof(c)
Cyc{Rational{Int64}}

julia> typeof(1+E(4))
Cyc{Int64}

julia> Cyc(1+im) # one can convert Gaussian integers or rationals
1+ζ₄

julia> 1//(1+E(4))
1//2-1//2ζ₄

julia> typeof(Cyc(1//2)) # another way of building a Cyc
Cyc{Rational{Int64}}

julia> conj(1+E(4))
1-ζ₄

julia> c=E(9)     # an effect of the Zumbroich basis
-ζ₉⁴-ζ₉⁷

julia> AsRootOfUnity(c) # but you can decide if a Cyc is a root of unity
1//9

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Cycs.jl#L1-L112' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Cycs.jl#L363-L375' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Cycs.jl#L400-L410' class='documenter-source'>source</a><br>

<a id='Gapjm.Cycs.quadratic' href='#Gapjm.Cycs.quadratic'>#</a>
**`Gapjm.Cycs.quadratic`** &mdash; *Function*.



quadratic(c::Cyc) determines if c lives in a quadratic extension of Q   it  returns a named  tuple (a=a,b=b,root=root,d=d) of  integers such that   c=(a + b ER(root))//d or nothing if no such tuple exists

**Examples**

```julia-repl
julia> quadratic(1+E(3))
(a = 1, b = 1, root = -3, den = 2)

julia> quadratic(1+E(5))

```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Cycs.jl#L441-L454' class='documenter-source'>source</a><br>


<a id='Pols.jl-Documentation-1'></a>

# Pols.jl Documentation

<a id='Gapjm.Pols' href='#Gapjm.Pols'>#</a>
**`Gapjm.Pols`** &mdash; *Module*.



An implementation of univariate Laurent polynomials.  A Pol contains two fields: its vector of coefficients, and its valuation.

**Examples**

```julia-repl
julia> Pol([1,2],0) # coefficients should have no leading or trailing zeroes.
1+2x

julia> p=Pol([1,2],-1)
x^-1+2

julia> valuation(p)
-1

julia> Pol(:q) # change string used for printing and set variable q
q

julia> p=(q+1)^2
1+2q+q^2

julia> degree(p)
2

julia> p(1//2) # a Pol is a callable object, where the call evaluates the Pol
9//4

julia> divrem(q^3+1,q+2) # changes coefficients to field elements
(4.0-2.0q+1.0q^2, -7.0)

julia> divrem1(q^3+1,q+2) # keeps the ring, but needs second argument unitary
(4-2q+q^2, -7)

julia> cyclotomic_polynomial(24) # the 24-th cyclotomic polynomial
1-q^4+q^8

```

see also the individual documentation of gcd.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Pols.jl#L1-L40' class='documenter-source'>source</a><br>

<a id='Base.divrem' href='#Base.divrem'>#</a>
**`Base.divrem`** &mdash; *Function*.



computes (p,q) such that a=p*b+q


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Pols.jl#L154-L156' class='documenter-source'>source</a><br>

<a id='Gapjm.Pols.divrem1' href='#Gapjm.Pols.divrem1'>#</a>
**`Gapjm.Pols.divrem1`** &mdash; *Function*.



divrem when b unitary: does not change type


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Pols.jl#L172-L174' class='documenter-source'>source</a><br>

<a id='Base.gcd' href='#Base.gcd'>#</a>
**`Base.gcd`** &mdash; *Function*.



gcd(p::Pol, q::Pol)   the coefficients of p and q must be elements of a field for   gcd to be type-stable

**Examples**

```julia-repl
julia> gcd(q+1,q^2-1)
1.0+1.0q

julia> gcd(q+1//1,q^2-1//1)
(1//1)+(1//1)q
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Pols.jl#L194-L207' class='documenter-source'>source</a><br>


<a id='Util.jl-Documentation-1'></a>

# Util.jl Documentation

<a id='Gapjm.Util' href='#Gapjm.Util'>#</a>
**`Gapjm.Util`** &mdash; *Module*.



This  module contains  various utility  functions used  in the  rest of the code.  Maybe some  of them  exist in  some Julia  module I am not aware of; please tell me.

The code is divided in sections  according to semantics.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L1-L7' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.cartesian' href='#Gapjm.Util.cartesian'>#</a>
**`Gapjm.Util.cartesian`** &mdash; *Function*.



Cartesian product of list of n vectors, returned as an x-by-n matrix


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L40' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L46-L56' class='documenter-source'>source</a><br>


group items of list l according to the values taken by function f on them

```
julia> groupby(iseven,1:10)
Dict{Bool,Array{Int64,1}} with 2 entries:
  false => [1, 3, 5, 7, 9]
  true  => [2, 4, 6, 8, 10]
```

Note:in this version l is required to be non-empty since I do not know how to access the return type of a function


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L65-L75' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.constant' href='#Gapjm.Util.constant'>#</a>
**`Gapjm.Util.constant`** &mdash; *Function*.



whether all elements in list a are equal


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L84' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.blocks' href='#Gapjm.Util.blocks'>#</a>
**`Gapjm.Util.blocks`** &mdash; *Function*.



blocks(M::Matrix)

M  should be a square matrix. Define  a graph G with vertices 1:size(M,1)   and  with an edge between i and j  if either M[i,j] or M[j,i] is not zero   or false. blocks returns a vector of vectors I such that I[1],I[2], etc..   are  the  vertices  in  each  connected  component  of G. In other words,   M[I[1],I[1]],M[I[2],I[2]],etc... are blocks of M.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L89-L97' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.format' href='#Gapjm.Util.format'>#</a>
**`Gapjm.Util.format`** &mdash; *Function*.



format( table; options )

General routine to format a table. Used for character tables.   Options:      row*labels          Labels for rows      column*labels       Labels for columns      rows*label          Label for column of rowLabels      separators          line numbers after which to put a separator      column*repartition  display in pieces of sizes these numbers of cols      rows                show only these rows      columns             show only these columns


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L196-L209' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.prime_residues' href='#Gapjm.Util.prime_residues'>#</a>
**`Gapjm.Util.prime_residues`** &mdash; *Function*.



the numbers less than n and prime to n 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L258' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.phi' href='#Gapjm.Util.phi'>#</a>
**`Gapjm.Util.phi`** &mdash; *Function*.



the Euler function ϕ 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L275' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.primitiveroot' href='#Gapjm.Util.primitiveroot'>#</a>
**`Gapjm.Util.primitiveroot`** &mdash; *Function*.



primitiveroot(m::Integer) a primitive root mod. m,   that is it generates multiplicatively prime_residues(m).   It exists if m is of the form 4, 2p^a or p^a for p prime>2.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/e822b0d4af62956d14c2f3e753b7a916692edb0c/src/Util.jl#L281-L285' class='documenter-source'>source</a><br>

