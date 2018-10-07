
<a id='Gapjm.jl-Documentation-1'></a>

# Gapjm.jl Documentation

<a id='Gapjm' href='#Gapjm'>#</a>
**`Gapjm`** &mdash; *Module*.



Here are some of my efforts porting GAP code to julia. I am not even sure this is a well-formed Julia package.

It contains for now permutations and permutation groups, cyclotomic numbers and  Laurent polynomials. Coming  soon are Coxeter  groups, Hecke algebras, braid groups and Garside monoids.

Even  though the code  is often competitive  with or faster  than GAP, I am sure there are more optimisations possible. Any comments about the code and the design are welcome.


<a id='Perms.jl-Documentation-1'></a>

# Perms.jl Documentation

<a id='Perms.Perms' href='#Perms.Perms'>#</a>
**`Perms.Perms`** &mdash; *Module*.



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

<a id='Perms.cycles' href='#Perms.cycles'>#</a>
**`Perms.cycles`** &mdash; *Function*.



cycles(a::Perm) returns the non-trivial cycles of a

**Example**

```julia-repl
julia> cycles(Perm(1,2)*Perm(4,5))
2-element Array{Array{Int64,1},1}:
 [1, 2]
 [3]
 [3, 4]
```

<a id='Perms.cycletype' href='#Perms.cycletype'>#</a>
**`Perms.cycletype`** &mdash; *Function*.



cycletype(a::Perm) is the partition of degree(a) associated to the   conjugacy class of a in the symmetric group, with ones removed

**Example**

```julia-repl
julia> cycletype(Perm(1,2)*Perm(3,4))
2-element Array{Int64,1}:
 2
 2

```

<a id='Base.sign' href='#Base.sign'>#</a>
**`Base.sign`** &mdash; *Function*.



```
sign(x)
```

Return zero if `x==0` and $x/|x|$ otherwise (i.e., ±1 for real `x`).


<a target='_blank' href='https://github.com/JuliaLang/julia/blob/0d713926f85dfa3e4e0962215b909b8e47e94f48/base/number.jl#L136-L140' class='documenter-source'>source</a><br>


sign(a::Perm) is the signature of  the permutation a


<a id='PermGroups.jl-Documentation-1'></a>

# PermGroups.jl Documentation

<a id='PermGroups.PermGroups' href='#PermGroups.PermGroups'>#</a>
**`PermGroups.PermGroups`** &mdash; *Module*.



This module is a port of some GAP functionality on permutation groups.

See Holt "Handbook of computational group theory" chap. 4 for basic algorithms.

The  only  field  of  a  PermGroup  G  at  the  start  is gens, the list of generators  of G.  To mimic  GAP records  where attributes/properties of an object  are computed on demand when asked for, other fields are computed on demand  and stored in the  field prop of the  PermGroup, which starts as an empty dict.

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
 [2]      
 [1]      
 [2, 1]   
 [1, 2]   
 [1, 2, 1]

julia> elements(G) # elements in the same order as words
6-element Array{Perm{Int64},1}:
 ()     
 (2,3)  
 (1,2)  
 (1,2,3)
 (1,3,2)
 (1,3)  

# finally, benchmarks
julia> @btime collect(symmetric_group(8));
  10.252 ms (350529 allocations: 14.17 MiB)

julia> @btime words(symmetric_group(8));
  111.824 ms (1596449 allocations: 38.64 MiB)
```

<a id='PermGroups.symmetric_group' href='#PermGroups.symmetric_group'>#</a>
**`PermGroups.symmetric_group`** &mdash; *Function*.



The symmetric group of degree n 

<a id='PermGroups.orbit' href='#PermGroups.orbit'>#</a>
**`PermGroups.orbit`** &mdash; *Function*.



orbit(G,p) is the orbit of point p under PermGroup G

<a id='PermGroups.orbit_and_representative' href='#PermGroups.orbit_and_representative'>#</a>
**`PermGroups.orbit_and_representative`** &mdash; *Function*.



returns Dict x=>g where x runs over orbit(G,p) and g is such that x=p^g

<a id='PermGroups.words' href='#PermGroups.words'>#</a>
**`PermGroups.words`** &mdash; *Function*.



List of minimal words in the generators elements(G) 

<a id='PermGroups.elements' href='#PermGroups.elements'>#</a>
**`PermGroups.elements`** &mdash; *Function*.



The list of elements of G in the same order as words


<a id='Cycs.jl-Documentation-1'></a>

# Cycs.jl Documentation

<a id='Cycs.Cycs' href='#Cycs.Cycs'>#</a>
**`Cycs.Cycs`** &mdash; *Module*.



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
E(12)^4-E(12)^7-E(12)^11

julia> E(3,2)
E(3)^2

julia> 1+E(3,2)
-E(3)

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

julia> inv(1+E(4)) # inverses need Rationals
(1/2)+(-1/2)E(4)

julia> typeof(ans)
Cyc{Rational{Int64}}

julia> typeof(1+E(4))
Cyc{Int64}

ulia> Cyc(1+im) # one can convert Gaussian integers or rationals
1+E(4)

julia> 1//(1+E(4))
(1/2)+(-1/2)E(4)

julia> typeof(Cyc(1//2)) # another way of building a Cyc
Cyc{Rational{Int64}}

julia> conj(1+E(4))
1-E(4)

julia> E(9)     # an effect of the Zumbroich basis
-E(9)^4-E(9)^7

julia> AsRootOfUnity(ans) # but you can decide if a Cyc is a root of unity
1/9

julia> Complex(E(3))   # convert to float is probably not very useful
-0.4999999999999998 + 0.8660254037844387im

julia> Cyc(ans) # even less useful
-0.4999999999999998+0.8660254037844387E(4)
```

For more information look at the documentation of ER, quadratic, galois. 

<a id='Cycs.galois' href='#Cycs.galois'>#</a>
**`Cycs.galois`** &mdash; *Function*.



galois(c::Cyc,n::Int) applies to c the galois automorphism   of Q(ζ_conductor(c)) raising all roots of unity to the n-th power.   n should be prime to c.

**Examples**

```julia-repl
julia> galois(1+E(4),-1) # galois(c,-1) is the same as conj(c)
1-E(4)

julia> galois(ER(5),2)==-ER(5)
true
```

<a id='Cycs.ER' href='#Cycs.ER'>#</a>
**`Cycs.ER`** &mdash; *Function*.



ER(n::Int) computes as a Cyc the square root of the integer n.

**Examples**

```julia-repl
julia> ER(-3)
E(3)-E(3)^2

julia> ER(3)
-E(12)^7+E(12)^11
```

<a id='Cycs.quadratic' href='#Cycs.quadratic'>#</a>
**`Cycs.quadratic`** &mdash; *Function*.



quadratic(c::Cyc) determines if c lives in a quadratic extension of Q   it returns tuple (a,b,root,d) of integers such that  c=(a + b ER(root))//d    or false if no such tuple exists

**Examples**

```julia-repl
julia> quadratic(1+E(3))
(1, 1, -3, 2)

julia> quadratic(1+E(5))
false
```


<a id='Pols.jl-Documentation-1'></a>

# Pols.jl Documentation

<a id='Pols.Pols' href='#Pols.Pols'>#</a>
**`Pols.Pols`** &mdash; *Module*.



An implementation of univariate Laurent polynomials.  A Pol contains two fields: its vector of coefficients, and its valuation.

**Examples**

```julia-repl
julia> Pol([1,2],0) # coefficients should have no leading or trailing zeroes. 
1+2x

julia> Pol([1,2],-1)
x^-1+2

julia> valuation(ans)
-1

julia> Pols.varname(:q) # change string used for printing and set variable q
:q

julia> p=(q+1)^2
1+2q+q^2

julia> degree(p)
2

julia> value(p,1//2)
9//4

julia> divrem(q^3+1,q+2) # changes coefficients to field elements
(4.0-2.0q+1.0q^2, -7.0)

julia> divrem1(q^3+1,q+2) # keeps the ring, but needs second argument unitary
(4-2q+q^2, -7)

julia> cyclotomic_polynomial(24) # the 24-th cyclotomic polynomial
1-q^4+q^8

```

see also the individual documentation of gcd.

<a id='Base.divrem' href='#Base.divrem'>#</a>
**`Base.divrem`** &mdash; *Function*.



```
divrem(x, y)
```

The quotient and remainder from Euclidean division. Equivalent to `(div(x,y), rem(x,y))` or `(x÷y, x%y)`.

**Examples**

```julia-repl
julia> divrem(3,7)
(0, 3)

julia> divrem(7,3)
(2, 1)
```


<a target='_blank' href='https://github.com/JuliaLang/julia/blob/0d713926f85dfa3e4e0962215b909b8e47e94f48/base/number.jl#L90-L104' class='documenter-source'>source</a><br>


computes (p,q) such that a=p*b+q

<a id='Pols.divrem1' href='#Pols.divrem1'>#</a>
**`Pols.divrem1`** &mdash; *Function*.



divrem when b unitary: does not change type

<a id='Base.gcd' href='#Base.gcd'>#</a>
**`Base.gcd`** &mdash; *Function*.



```
gcd(x,y)
```

Greatest common (positive) divisor (or zero if `x` and `y` are both zero).

**Examples**

```julia-repl
julia> gcd(6,9)
3

julia> gcd(6,-9)
3
```


<a target='_blank' href='https://github.com/JuliaLang/julia/blob/0d713926f85dfa3e4e0962215b909b8e47e94f48/base/intfuncs.jl#L5-L18' class='documenter-source'>source</a><br>


gcd(p::Pol, q::Pol)   the coefficients of p and q must be elements of a field for    gcd to be type-stable

**Examples**

```julia-repl
julia> gcd(q+1,q^2-1)
1.0+1.0q

julia> gcd(q+1//1,q^2-1//1)
1+q
```


<a id='Util.jl-Documentation-1'></a>

# Util.jl Documentation

<a id='Util.Util' href='#Util.Util'>#</a>
**`Util.Util`** &mdash; *Module*.



This  module contains  various utility  functions used  in the  rest of the code.  Maybe some  of them  exist in  some Julia  module I am not aware of; please tell me.

The code is divided in sections  according to semantics.

<a id='Util.cartesian' href='#Util.cartesian'>#</a>
**`Util.cartesian`** &mdash; *Function*.



Cartesian product of list of n vectors, returned as an x-by-n matrix

<a id='Util.groupby' href='#Util.groupby'>#</a>
**`Util.groupby`** &mdash; *Function*.



group items of list l according to the corresponding values in list v

```
julia> groupby([31,28,31,30,31,30,31,31,30,31,30,31],
       [:Jan,:Feb,:Mar,:Apr,:May,:Jun,:Jul,:Aug,:Sep,:Oct,:Nov,:Dec])
Dict{Int64,Array{Symbol,1}} with 3 entries:
  31 => Symbol[:Jan, :Mar, :May, :Jul, :Aug, :Oct, :Dec]
  28 => Symbol[:Feb]
  30 => Symbol[:Apr, :Jun, :Sep, :Nov]
```


group items of list l according to the values taken by function f on them

```
julia> groupby(iseven,1:10)
Dict{Bool,Array{Int64,1}} with 2 entries:
  false => [1, 3, 5, 7, 9]
  true  => [2, 4, 6, 8, 10]
```

Note:in this version l is required to be non-empty since I do not know how to access the return type of a function

<a id='Util.constant' href='#Util.constant'>#</a>
**`Util.constant`** &mdash; *Function*.



whether all elements in list a are equal

<a id='Util.blocks' href='#Util.blocks'>#</a>
**`Util.blocks`** &mdash; *Function*.



blocks(M::Matrix)

M  should be a square matrix. Define  a graph G with vertices 1:size(M,1)   and  with an edge between i and j  if either M[i,j] or M[j,i] is not zero   or false. blocks returns a vector of vectors I such that I[1],I[2], etc..   are  the  vertices  in  each  connected  component  of G. In other words,   M[I[1],I[1]],M[I[2],I[2]],etc... are blocks of M.

<a id='Util.format' href='#Util.format'>#</a>
**`Util.format`** &mdash; *Function*.



format( table; options )

General routine to format a table. Used for character tables.   Options:      row*labels          Labels for rows      column*labels       Labels for columns      rows*label          Label for column of rowLabels      separators          line numbers after which to put a separator      column*repartition  display in pieces of sizes these numbers of cols      rows                show only these rows      columns             show only these columns

<a id='Util.prime_residues' href='#Util.prime_residues'>#</a>
**`Util.prime_residues`** &mdash; *Function*.



the numbers less than n and prime to n 

<a id='Util.phi' href='#Util.phi'>#</a>
**`Util.phi`** &mdash; *Function*.



the Euler function ϕ 

<a id='Util.primitiveroot' href='#Util.primitiveroot'>#</a>
**`Util.primitiveroot`** &mdash; *Function*.



primitiveroot(m::Integer) a primitive root mod. m,    that is it generates multiplicatively prime_residues(m).     It exists if m is of the form 4, 2p^a or p^a for p prime>2.

