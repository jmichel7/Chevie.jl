
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

To update to the latest version, I do not know a better way that to rm the package and re-install it.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Gapjm.jl#L1-L29' class='documenter-source'>source</a><br>

- [Gapjm.jl Documentation](README.md#Gapjm.jl-Documentation-1)
- [Perms.jl Documentation](README.md#Perms.jl-Documentation-1)
- [PermGroups.jl Documentation](README.md#PermGroups.jl-Documentation-1)
- [Cycs.jl Documentation](README.md#Cycs.jl-Documentation-1)
- [Pols.jl Documentation](README.md#Pols.jl-Documentation-1)
- [CoxGroups.jl Documentation](README.md#CoxGroups.jl-Documentation-1)
- [Weyl.jl Documentation](README.md#Weyl.jl-Documentation-1)
- [Hecke.jl Documentation](README.md#Hecke.jl-Documentation-1)
- [KL.jl Documentation](README.md#KL.jl-Documentation-1)
- [Garside.jl Documentation](README.md#Garside.jl-Documentation-1)
- [Util.jl Documentation](README.md#Util.jl-Documentation-1)
- [Cycpols.jl Documentation](README.md#Cycpols.jl-Documentation-1)


<a id='Perms.jl-Documentation-1'></a>

# Perms.jl Documentation

<a id='Gapjm.Perms' href='#Gapjm.Perms'>#</a>
**`Gapjm.Perms`** &mdash; *Module*.



This module is a port of some GAP functionality on permutations.

A  permutation here is a permutation of the set 1:n and is represented as a list  of n integers representing the images of 1:n. The integer n is called the *degree* of the permutation.

Permutations  in  this  module  follow  the  GAP  design: it is possible to multiply, or to store in the same group, permutations of different degrees. A  slightly faster  design is  the MAGMA  one where  any permutation has to belong  to  a  group  and  the  degree  is  determined by that group. There multiplication of permutations in a given group is a faster, but it is more difficult  to multiply  permutations coming  from different  groups, like a group and one of its subgroups.

The  GAP permutation  (1,2,3)(4,5) can  be written Perm(1,2,3)*Perm(4,5) or perm"(1,2,3)(4,5)".  It is represented internally as [2,3,1,5,4]; note that [2,3,1,5,4,6] represents the same permutation.

As in GAP i^p applies p to integer i, while p^q means p^-1*q*p.

Another  Perm  constructor  is  Perm{T}(p)  which  converts the perm p to a permutation  on  integers  of  type  T;  for  instance  Perm{UInt8} is more efficient  that Perm{Int} and can be used for Weyl groups of rank <=8 since they have at most 240 roots.

**Examples**

```julia-repl
julia> p=Perm(1,2)*Perm(2,3)
(1,3,2)

julia> Perm{Int8}(p)
{Int8}(1,3,2)

julia> 1^p
3

julia> Matrix(p)
3×3 Array{Int64,2}:
 0  0  1
 1  0  0
 0  1  0

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
```

Perms  have methods copy, hash,  ==, cmp, isless (total order)  so they can be keys in hashes or elements of sets.

other functions are: cycles, cycletype, sign. See individual documentation.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Perms.jl#L1-L70' class='documenter-source'>source</a><br>

<a id='Gapjm.Perms.cycles' href='#Gapjm.Perms.cycles'>#</a>
**`Gapjm.Perms.cycles`** &mdash; *Function*.



cycles(a::Perm) returns the cycles of a

**Example**

```julia-repl
julia> cycles(Perm(1,2)*Perm(4,5))
3-element Array{Array{Int64,1},1}:
 [1, 2]
 [3]
 [4, 5]
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Perms.jl#L220-L230' class='documenter-source'>source</a><br>

<a id='Gapjm.Perms.cycletype' href='#Gapjm.Perms.cycletype'>#</a>
**`Gapjm.Perms.cycletype`** &mdash; *Function*.



cycletype(a::Perm) describes the partition of degree(a) associated to the   conjugacy class of a in the symmetric group, with ones removed

**Example**

```julia-repl
julia> cycletype(Perm(1,2)*Perm(3,4))
1-element Array{Pair{Tuple{Int64,Int64},Int64},1}:
 (2, 1) => 2
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Perms.jl#L260-L269' class='documenter-source'>source</a><br>

<a id='Base.sign' href='#Base.sign'>#</a>
**`Base.sign`** &mdash; *Function*.



sign(a::Perm) is the signature of  the permutation a


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Perms.jl#L292' class='documenter-source'>source</a><br>


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

# PermGroups are iterators over their elements
julia> collect(G)  
6-element Array{Perm{Int64},1}:
 (1,2)
 (1,3,2)
 ()
 (1,2,3)
 (1,3)
 (2,3)

# maximum degree of an element of G
julia> degree(G)  
3

# orbit of point 1 under G
julia> orbit(G,1) 
3-element Array{Int64,1}:
 2
 3
 1

# orbit decorated with representatives moving 1 to given point
julia> orbit_and_representative(G,1)
Dict{Int64,Perm{Int64}} with 3 entries:
  2 => (1,2)
  3 => (1,3,2)
  1 => ()

# orbit functions can take any action of G as keyword argument
julia> orbit_and_representative(G,[1,2],action=(x,y)->x.^Ref(y))
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

# a list of points that no element of G fixes
julia> base(G) 
2-element Array{Int64,1}:
 1
 2

# the i-th element is the centralizer of base[1:i-1]
julia> centralizers(G) 
2-element Array{PermGroup{Int64},1}:
 PermGroup((1,2),(2,3))
 PermGroup((2,3))

# i-th element is orbit_and_representative of centralizer[i] on base[i]
julia> centralizer_orbits(G)
2-element Array{Dict{Int64,Perm{Int64}},1}:
 Dict(2=>(1,2),3=>(1,3,2),1=>())
 Dict(2=>(),3=>(2,3))

julia> minimal_words(G)  # minimal word in gens for each element of G
Dict{Perm{Int64},Array{Int64,1}} with 6 entries:
  ()      => Int64[]
  (2,3)   => [2]
  (1,3,2) => [1, 2]
  (1,3)   => [1, 2, 1]
  (1,2)   => [1]
  (1,2,3) => [2, 1]
```

finally, benchmarks on julia 1.0.1

```benchmark
julia> @btime length(collect(symmetric_group(8)))
  5.481 ms (270429 allocations: 12.40 MiB)

julia> @btime minimal_words(symmetric_group(8));
  10.477 ms (122062 allocations: 15.22 MiB)
```

Compare to GAP3 Elements(SymmetricGroup(8)); takes 3.8 ms


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/PermGroups.jl#L1-L106' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.symmetric_group' href='#Gapjm.PermGroups.symmetric_group'>#</a>
**`Gapjm.PermGroups.symmetric_group`** &mdash; *Function*.



The symmetric group of degree n 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/PermGroups.jl#L251' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.orbit' href='#Gapjm.PermGroups.orbit'>#</a>
**`Gapjm.PermGroups.orbit`** &mdash; *Function*.



orbit(G,p) is the orbit of p under Group G


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/PermGroups.jl#L126' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.orbit_and_representative' href='#Gapjm.PermGroups.orbit_and_representative'>#</a>
**`Gapjm.PermGroups.orbit_and_representative`** &mdash; *Function*.



returns Dict x=>g for x in orbit(G,p) giving g such that x=action(p,g)


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/PermGroups.jl#L142' class='documenter-source'>source</a><br>

<a id='Gapjm.PermGroups.minimal_words' href='#Gapjm.PermGroups.minimal_words'>#</a>
**`Gapjm.PermGroups.minimal_words`** &mdash; *Function*.



dict giving for each element of G a minimal word in the generators


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/PermGroups.jl#L170' class='documenter-source'>source</a><br>


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
  472.964 ms (8324504 allocations: 707.18 MiB)
```

The equivalent in GAP:

```
testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
end; 
```

for testmat(12) takes 0.4s in GAP3, 0.3s in GAP4


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Cycs.jl#L1-L109' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Cycs.jl#L576-L588' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Cycs.jl#L621-L631' class='documenter-source'>source</a><br>

<a id='Gapjm.Cycs.quadratic' href='#Gapjm.Cycs.quadratic'>#</a>
**`Gapjm.Cycs.quadratic`** &mdash; *Function*.



quadratic(c::Cyc) determines if c lives in a quadratic extension of Q   it  returns a named  tuple (a=a,b=b,root=root,d=d) of  integers such that   c=(a + b ER(root))//d or nothing if no such tuple exists

**Examples**

```julia-repl
julia> quadratic(1+E(3))
(a = 1, b = 1, root = -3, den = 2)

julia> quadratic(1+E(5))

```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Cycs.jl#L715-L728' class='documenter-source'>source</a><br>


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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Pols.jl#L1-L40' class='documenter-source'>source</a><br>

<a id='Base.divrem' href='#Base.divrem'>#</a>
**`Base.divrem`** &mdash; *Function*.



computes (p,q) such that a=p*b+q


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Pols.jl#L178-L180' class='documenter-source'>source</a><br>

<a id='Gapjm.Pols.divrem1' href='#Gapjm.Pols.divrem1'>#</a>
**`Gapjm.Pols.divrem1`** &mdash; *Function*.



divrem when b unitary: does not change type


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Pols.jl#L196-L198' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Pols.jl#L227-L240' class='documenter-source'>source</a><br>


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

  * `gens(W)` which returns the set `S` (the list of *Coxeter generators*)
  * `nref(W)` which  returns the  number of  reflections of  `W`, if  `W` is  finite or `nothing` if `W` is infinite

It  should be  noted that  a Coxeter group can be *any* kind of group implementing the above functions.

A  common occurrence in code for Coxeter groups is a loop like:

`findfirst(eachindex(gens(W)),x->isleftdescent(W,w,x))`

if you provide a function `firstleftdescent(W,w)` it will be called instead of the above loop.

Because  of the  easy solution  of the  word problem  in Coxeter  groups, a convenient  way  to  represent  their  elements  is as words in the Coxeter generators.  They are represented as lists of labels for the generators. By default  these labels are  given as the  index of a  generator in `S`, so a Coxeter  word is just  a list of  integers in `1:length(S)`. For reflection subgroups, the labels are indices of the reflections in the parent group.

The functions 'word' and 'W(...)' will do the conversion between Coxeter words and elements of the group.

**Examples**

```julia-repl
julia> W=coxsym(4)
coxsym(4)

julia> p=W(1,3,2,1,3)
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

This  file contains mostly a port of  the basic functions on Coxeter groups in  CHEVIE. The only Coxeter group  constructor implemented here is coxsym. The file Weyl.jl defines coxgroup.

The dictionary from CHEVIE is as follows:

```
     CoxeterElements(W[,l])                → elements(W[,l])
     CoxeterLength(W,w)                    → length(W,w)
     CoxeterWord(W,w)                      → word(W,w)
     LongestCoxeterElement(W)              → longest(W)
     FirstLeftDescending(W,w)              → firstleftdescent(W,w)
     LeftDescenTSet(W,w)                   → leftdescents(W,w)
     ReducedInRightCoset(W,w)              → reduced(W,w)
     ReducedRightCosetRepresentatives(W,H) → reduced(H,W)
     SemiSimpleRank(W)                     → coxrank(W)
     CoxeterGroupSymmetricGroup(n)         → coxsym(n)
     ReflectionSubGroup                    only standard parabolics now
     IsLeftDescending(W,w,i)               → isleftdescent(W,w,i)
     ReflectionDegrees(W)                  → degrees(W)
     ReflectionLength(W,w)                 → reflength(W,w)
     W.N                                   → nref(W)
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CoxGroups.jl#L1-L158' class='documenter-source'>source</a><br>

<a id='Gapjm.CoxGroups.reduced' href='#Gapjm.CoxGroups.reduced'>#</a>
**`Gapjm.CoxGroups.reduced`** &mdash; *Function*.



reduced(W,w)   The unique element in the coset W.w which stabilises the positive roots of W

```julia-repl
julia> W=coxgroup(:G,2)
W(G₂)

julia> H=reflection_subgroup(W,[2,6])
W(G₂)₂₄

julia> Set(word.(Ref(W),reduced.(Ref(H),elements(W))))
Set(Array{Int64,1}[[1], []])

```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CoxGroups.jl#L216-L230' class='documenter-source'>source</a><br>


reduced(H,W)   The elements in W which are H-reduced

```julia-repl
julia> W=coxgroup(:G,2)
W(G₂)

julia> H=reflection_subgroup(W,[2,6])
W(G₂)₂₄

julia> [word(W,w) for S in reduced(H,W) for w in S]
2-element Array{Array{Int64,1},1}:
 [] 
 [1]

```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CoxGroups.jl#L239-L255' class='documenter-source'>source</a><br>


reduced(H,W,S)   The elements in W which are H-reduced of length i from the set S of length i-1


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CoxGroups.jl#L267-L270' class='documenter-source'>source</a><br>

<a id='Gapjm.CoxGroups.bruhatless' href='#Gapjm.CoxGroups.bruhatless'>#</a>
**`Gapjm.CoxGroups.bruhatless`** &mdash; *Function*.



`bruhatless(W, x, y)`  whether x≤y in the Bruhat order, for x, y ∈ W.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CoxGroups.jl#L345-L347' class='documenter-source'>source</a><br>

<a id='Gapjm.CoxGroups.coxsym' href='#Gapjm.CoxGroups.coxsym'>#</a>
**`Gapjm.CoxGroups.coxsym`** &mdash; *Function*.



The symmetric group on n letters as a Coxeter group


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CoxGroups.jl#L425' class='documenter-source'>source</a><br>

<a id='Gapjm.CoxGroups.longest' href='#Gapjm.CoxGroups.longest'>#</a>
**`Gapjm.CoxGroups.longest`** &mdash; *Function*.



The longest element of reflection_subgroup(W,I) –- never ends if infinite


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CoxGroups.jl#L201-L203' class='documenter-source'>source</a><br>

<a id='Gapjm.CoxGroups.nref' href='#Gapjm.CoxGroups.nref'>#</a>
**`Gapjm.CoxGroups.nref`** &mdash; *Function*.



number of reflections of W


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Weyl.jl#L441' class='documenter-source'>source</a><br>


<a id='Weyl.jl-Documentation-1'></a>

# Weyl.jl Documentation

<a id='Gapjm.Weyl' href='#Gapjm.Weyl'>#</a>
**`Gapjm.Weyl`** &mdash; *Module*.



Let  `V` be a  real vector space.  Finite Coxeter groups  coincide with teh finite  subgroups of  `GL(V)` which  can be  generated y reflections. *Weyl groups*  are  the  finite  Coxeter  groups  which  can  be defined over the rational   numbers.  We  implement  finite   Coxter  groups  as  groups  of permutations  of  a  root  system.  Root  systems play an important role in mathematics as they classify semi-simple Lie algebras and algebraic groups.

Let  us give precise definitions. Let `V`  be a real vector space, `Vⱽ` its dual  and let `(,)`  be the natural  pairing between `Vⱽ`  and `V`. A *root system*  is a finite set  of vectors `R` which  generate `V` (the *roots*), together  with  a  map  `r↦  rⱽ`  from  `R`  to  a subset `Rⱽ` of `Vⱽ` (the *coroots*) such that:

  * For any `r∈  R`, we have  `(rⱽ,r)=2` so that  the formula `x↦ x-(rⱽ,x)r`

defines a reflection `s_r:V→ V` with root `r` and coroot `rⱽ`.

  * The reflection `s_r` stabilizes `R`.

We  will only  consider *reduced*  root systems,  i.e., such  that the only elements  of `R` colinear with `r∈ R` are `r` and `-r`; for Weyl groups, we also ask that the root system be *crystallographic*, that is `(rⱽ,s)` is an integer, for any `s∈ R,rⱽ∈ Rⱽ`.

The  subgroup `W=W(R)` of  `GL(V)` generated by  the reflections `s_r` is a finite  Coxeter group; when `R` is crystallographic, the representation `V` of  `W`  is  defined  over  the  rational  numbers.  All finite-dimensional (complex)  representations of a  finite Coxeter group  can be realized over the  same field  as `V`.  Weyl groups  can be  characterized amongst finite Coxeter  groups by the fact that all numbers `m(s,t)` in the Coxeter matrix are in `{2,3,4,6}`.

If  we identify  `V` with  `Vⱽ` by  choosing a  `W`-invariant bilinear form `(.;.)`;  then we have `rⱽ=2r/(r;r)`. A root system `R` is *irreducible* if it is not the union of two orthogonal subsets. If `R` is reducible then the corresponding  Coxeter group  is the  direct product  of the Coxeter groups associated with the irreducible components of `R`.

The  irreducible  crystallographic  root  systems  are  classified  by  the following  list of  *Dynkin diagrams*,  which, in  addition to  the Coxeter matrix,  encode also the relative length of the roots. We show the labeling of the nodes given by the function 'Diagram' described below.

```
A_n O—O—O—…—O   B_n O⇐O—O—…—O  C_n O⇒ O—O—…—O  D_n  O 2
    1 2 3 … n       1 2 3 … n      1  2 3 … n       ￨
                                                  O—O—…—O
                                                  1 3 … n

G₂ O⇛ O  F₄ O—O⇒ O—O  E₆   O 2   E₇   O 2     E₈    O 2
   1  2     1 2  3 4       ￨          ￨             ￨
                       O—O—O—O—O  O—O—O—O—O—O   O—O—O—O—O—O—O
                       1 3 4 5 6  1 3 4 5 6 7   1 3 4 5 6 7 8
```

These diagrams encode the presentation of the Coxeter group `W` as follows: the vertices represent the generating reflections; an edge is drawn between `s`  and `t` if the order `m(s,t)` of `st` is greater than `2`; the edge is single  if  `m(s,t)=3`,  double  if  `m(s,t)=4`,  triple if `m(s,t)=6`. The arrows  indicate the relative root lengths when `W` has more than one orbit on  `R`, as explained below; we  get the *Coxeter Diagram*, which describes the  underlying Weyl group, if  we ignore the arrows:  we see that the root systems `B_n` and `C_n` correspond to the same Coxeter group.

Here  are  the  diagrams  for  the  finite  Coxeter  groups which  are  not crystallographic:

```
   e        5         5
```

I₂(e) O—O   H₃ O—O—O  H₄ O—O—O—O       1 2      1 2 3     1 2 3 4 

Let us now describe how the root systems are encoded in these diagrams. Let `R`  be a root system in `V`. Then we can choose a linear form on `V` which vanishes  on no element of `R`. According to  the sign of the value of this linear  form on a root  `r ∈ R` we  call `r` *positive* or *negative*. Then there  exists a unique subset `Π` of  the positive roots, called the set of *simple  roots*, such that  any positive root  is a linear combination with non-negative  coefficients of  roots in  `Π`. Any  two sets of simple roots (corresponding  to  different  choices  of  linear  forms  as above) can be transformed into each other by a unique element of `W(R)`. Hence, since the pairing  between `V` and `Vⱽ`  is `W`-invariant, if `Π`  is a set of simple roots  and if  we define  the *Cartan  matrix* as  being the  `n` times `n` matrix   `C={rⱽ(s)}_{rs}`,  for  `r,s∈Π`  this   matrix  is  unique  up  to simultaneous  permutation of rows and columns.  It is precisely this matrix which is encoded in a Dynkin diagram, as follows.

The  indices for the rows of `C` label the nodes of the diagram. The edges, for  `r ≠ s`, are  given as follows. If  `C_{rs}` and `C_{sr}` are integers such  that `|C_{rs}|≥|C_{sr}|=1`  the vertices  are connected by `|C_{rs}|` lines,  and if `|C_{rs}|>1`  then we put  an additional arrow  on the lines pointing  towards the node with label `s`.  In other cases, we simply put a single   line  equipped  with  the  unique  integer  `p_{rs}≥1`  such  that `C_{rs}C_{sr}=cos^2 (π/p_{sr})`.

Conversely,  the whole root  system can be  recovered from the simple roots and  the corresponding coroots. The  reflections in `W(R)` corresponding to the  simple roots are called  *simple* reflections or *Coxeter generators*. They are precisely the generators for which the Coxeter diagram encodes the defining  relations of `W(R)`. Each root is  in the orbit of a simple root, so  that `R` is obtained  as the orbit of  the simple roots under the group generated  by  the  simple  reflections.  The  restriction  of  the  simple reflections  to the span of `R` is  determined by the Cartan matrix, so `R` is determined by the Cartan matrix and the set of simple roots.

The  Cartan  matrix  corresponding  to  one  of  the above irreducible root systems  (with the specified labeling) is  returned by the command 'cartan' which  takes as input  a `Symbol` giving  the type (that  is ':A', ':B', …, ':I')  and a positive `Int` giving the  rank (plus an `Int` giving the bond for  type `:I`).  This function  returns a  matrix with  entries in `ℤ` for crystallographic  types, and a  matrix of `Cyc`  for the other types. Given two  Cartan matrices `c1` and `c2`,  their matrix direct sum (corresponding to  the  orthogonal  direct  sum  of  the  root systems) can be produced by `cat(c1,c2,dims=[1,2])`.

The  function 'rootdatum' takes as input a  list of simple roots and a list of the corresponding coroots and produces a `struct` containing information about  the root system `R` and about `W(R)`. If we label the positive roots by  '1:N', and the negative roots  by 'N+1:2N', then each simple reflection is  represented by the permutation of '1:2N' which it induces on the roots. If  only one argument is given, the Cartan matrix of the root system, it is taken  as the list  of coroots and  the list of  roots is assumed to be the canonical basis of `V`.

If one only wants to work with Cartan matrices with a labeling as specified by  the  above  list,  the  function  call  can  be  simplified. Instead of 'rootdatum(CartanMat(:D,4))' the following is also possible.

```julia-repl
julia> W=coxgroup(:D,4)
W(D₄)

julia> cartan(W)
4×4 Array{Int64,2}:
  2   0  -1   0
  0   2  -1   0
 -1  -1   2  -1
  0   0  -1   2
```

Also,  the Weyl group struct associated to a direct sum of irreducible root systems can be obtained as a product

```julia-repl
julia> W=coxgroup(:A,2)*coxgroup(:B,2)
W(A₂)× W(B₂)₍₃₄₎

julia> cartan(W)
4×4 Array{Int64,2}:
  2  -1   0   0
 -1   2   0   0
  0   0   2  -2
  0   0  -1   2
```

The  same `struct`  is constructed  by applying  'coxgroup' to  the matrix 'cat(cartan(:A,2), cartan(:B,2),dims=[1,2])'.

The elements of a Weyl group are permutations of the roots:

```julia-repl
julia> W=coxgroup(:D,4)
W(D₄)

julia> p=W(1,3,2,1,3)
{Int16}(1,14,13,2)(3,17,8,18)(4,12)(5,20,6,15)(7,10,11,9)(16,24)(19,22,23,21)

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
     CoxeterGroup("A",5)                    →  coxgroup(:A,5) 
     Size(W)                                →  length(W) 
     ForEachElement(W,f)                    →  for w in W f(w) end 
     ReflectionDegrees(W)                   →  degrees(W) 
     IsLeftDescending(W,w,i)                →  isleftdescent(W,w,i) 
     ReflectionSubgroup                     →  reflection_subgroup
     TwoTree(m)                             →  twotree(m) 
     FiniteCoxeterTypeFromCartanMat(m)      →  type_cartan(m) 
     RootsCartan(m)                         →  roots(m) 
     PrintDiagram(W)                        →  Diagram(W) 
     Inversions                             →  inversions 
     Reflection                             →  reflection 
     W.orbitRepresentative[i]               →  simple_representative(W,i) 
```

finally, a benchmark on julia 1.0.2

```benchmark
julia> @btime length(elements(coxgroup(:E,7)))
  531.385 ms (5945569 allocations: 1.08 GiB)
```

GAP3 for the same computation takes 2.2s


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Weyl.jl#L1-L197' class='documenter-source'>source</a><br>

<a id='Gapjm.PermRoot.cartan' href='#Gapjm.PermRoot.cartan'>#</a>
**`Gapjm.PermRoot.cartan`** &mdash; *Function*.



```julia
`cartan(type, rank)`
```

Cartan matrix for a Weyl group:

```julia-repl
julia> cartan(:A,4)
4×4 Array{Int64,2}:
  2  -1   0   0
 -1   2  -1   0
  0  -1   2  -1
  0   0  -1   2
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Weyl.jl#L204-L217' class='documenter-source'>source</a><br>

<a id='Gapjm.Weyl.two_tree' href='#Gapjm.Weyl.two_tree'>#</a>
**`Gapjm.Weyl.two_tree`** &mdash; *Function*.



```julia
two_tree(m)
```

Given  a square  matrix m  with zeroes  (or falses,  for a boolean matrix)  symmetric  with respect to the diagonal, let  G be the graph with vertices  axes(m)[1] and an edge between i and j iff !iszero(m[i,j]).  If G  is a line this function returns it as a Vector{Int}.   If  G  is  a  tree  with  one  vertex  c of valence 3 the function returns  (c,b1,b2,b3)  where b1,b2,b3 are  the branches from  this vertex sorted by  increasing length.  Otherwise the function returns `nothing`

```julia-repl
julia> Weyl.two_tree(cartan(:A,4))
4-element Array{Int64,1}:
 1
 2
 3
 4

julia> Weyl.two_tree(cartan(:E,8))
(4, [2], [3, 1], [5, 6, 7, 8])
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Weyl.jl#L245-L267' class='documenter-source'>source</a><br>

<a id='Gapjm.PermRoot.reflection_subgroup' href='#Gapjm.PermRoot.reflection_subgroup'>#</a>
**`Gapjm.PermRoot.reflection_subgroup`** &mdash; *Function*.



reflection_subgroup(W,I) The subgroup of W generated by reflections(W)[I]

A   theorem  discovered  by  Deodhar  cite{Deo89}  and  Dyer  cite{Dye90} independently  is that a subgroup `H` of a Coxeter system `(W,S)` generated by  reflections has  a canonical  Coxeter generating  set, formed of the `t ∈Ref(H)`  such `l(tt')>l(t)` for any `t'∈  Ref(H)` different from `t`. This is used by 'reflection_subgroup' to determine the Coxeter system of `H`.

```julia-repl
julia> W=coxgroup(:G,2)
W(G₂)

julia> Diagram(W)
O⇛ O
1  2

julia> H=reflection_subgroup(W,[2,6])
W(G₂)₂₄

julia> Diagram(H)
O—O
1 2
```

The  notation `W(G₂)₂₃` means  that 'W.roots[2:3]' form  a system of simple roots for `H`.

A  reflection subgroup has specific properties  the most important of which is  'inclusion' which gives the positions of the roots of H in the roots of W. The inverse (partial) map is 'restriction'.

```julia-repl
julia> inclusion(H)
3-element Array{Int64,1}:
 2
 4
 6

julia> restriction(H)
12-element Array{Int64,1}:
 0
 1
 0
 2
 0
 3
 0
 0
 0
 0
 0
 0

```

If H is a standard parabolic subgroup  of a Coxeter group W then the length function  on H (with respect  to its set of  generators) is the restriction of  the length function on  W. This need not  no longer be true for arbitrary reflection subgroups of W:

```julia-repl
julia> word(W,H(2))
3-element Array{Int64,1}:
 1
 2
 1
```

In  this package, finite  reflection groups are  represented as permutation groups  on a set of roots. Consequently,  a reflection subgroup `H⊆ W` is a permutation  subgroup, thus its elements are represented as permutations of the roots of the parent group.

```julia-repl
julia> elH=word.(Ref(H),elements(H))
6-element Array{Array{Int64,1},1}:
 []       
 [2]      
 [1]      
 [2, 1]   
 [1, 2]   
 [1, 2, 1]

julia> elW=word.(Ref(W),elements(H))
6-element Array{Array{Int64,1},1}:
 []             
 [1, 2, 1]      
 [2]            
 [1, 2, 1, 2]   
 [2, 1, 2, 1]   
 [2, 1, 2, 1, 2]

julia> map(w->H(w...),elH)==map(w->W(w...),elW)
true

```

Another  basic result about reflection subgroups  of Coxeter groups is that each  coset of  H in  W contains  a unique  element of  minimal length, see `reduced`.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Weyl.jl#L539-L641' class='documenter-source'>source</a><br>


Only parabolics defined are I=1:m for m≤n


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CoxGroups.jl#L478' class='documenter-source'>source</a><br>

<a id='Gapjm.Weyl.coxgroup' href='#Gapjm.Weyl.coxgroup'>#</a>
**`Gapjm.Weyl.coxgroup`** &mdash; *Function*.



Coxeter group from type


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Weyl.jl#L446' class='documenter-source'>source</a><br>

<a id='Gapjm.Weyl.rootdatum' href='#Gapjm.Weyl.rootdatum'>#</a>
**`Gapjm.Weyl.rootdatum`** &mdash; *Function*.



Adjoint root datum from cartan mat


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Weyl.jl#L449' class='documenter-source'>source</a><br>


root datum


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Weyl.jl#L452' class='documenter-source'>source</a><br>


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

For  some  algebras  the  character  table,  and in general Kazhdan-Lusztig bases,  require a square root of `qₛ=-uₛ₀/uₛ₁`. We provide a way to specify it  with  the  field  `.rootpara`  which  can  be given when constructing the algebra. If not given a root is automatically extracted when needed (and we know  how to compute it) by the function `RootParameter`. Note however that sometimes  an  explicit  choice  of  root  is  necessary  which  cannot  be automatically determined.

There  is a universal choice  for `R` and `uₛᵢ`:  Let `uₛᵢ:s∈ S,i∈[0,1]` be indeterminates   such  that  `uₛᵢ=uₜᵢ`  whenever  `mₛₜ`  is  odd,  and  let `A=ℤ[uₛᵢ]` be the corresponding polynomial ring. Then the Hecke algebra `H` of  `W` over a  with parameters `uₛᵢ`  is called the *generic Iwahori-Hecke algebra*  of  with  `W`.  Any  other  algebra  with parameters `vₛᵢ` can be obtained  by specialization from  `H`: There is  a unique ring homomorphism `f:A  → R` such that `f(uₛᵢ)=vₛᵢ`  for all `i`. Then we  can view `R` as an `A`-module via `f` and we can identify the other algebra to $R⊗ _A H$.

The  elements `{T_w∣w∈ W}` actually form an  `R`-basis of `H` if one of the `uₛᵢ`  is invertible for all `s`. The  structure constants in that basis is obtained  as  follows.  To  multiply  `T_v`  by  `T_w`,  choose  a  reduced expression for `v`, say `v=s_1 ⋯ s_k` and apply inductively the formula:

$T_sT_w=T_{sw}$               if `l(sw)=l(w)+1`

$T_sT_w=-uₛ₀uₛ₁T_{sw}+(uₛ₀+uₛ₁)T_w$ if `l(sw)=l(w)-1`.

If all `s` we have `uₛ₀=q`, `uₛ₁=-1` then we call the corresponding algebra the one-parameter or Spetsial Iwahori-Hecke algebra associated with `W`; it can  be obtained with the  simplified call 'Hecke(W,q)'. Certain invariants of  the irreducible characters of  this algebra play a  special role in the representation  theory of the underlying  finite Coxeter groups, namely the `a`- and `A`-invariants. For basic properties of Iwahori-Hecke algebras and their  relevance to the representation theory of finite groups of Lie type, see for example Curtis and Reiner 1987, Sections~67 and 68.

In  the  following  example,  we  compute  the multiplication table for the `0`-Iwahori–Hecke algebra associated with the Coxeter group of type `A_2`.

```julia-repl
julia> W=coxgroup(:A,2)
W(A₂)

julia> H=hecke(W,0)             # One-parameter algebra with `q=0`
Hecke(W(A₂),0)

julia> T=Tbasis(H)              # Create the `T` basis
(::getfield(Gapjm.Hecke, Symbol("#f#32")){Int64,Perm{Int16},HeckeAlgebra{Int64,Gapjm.Weyl.FiniteCoxeterGroup{Int16,Int64}}}) (generic function with 4 methods)

julia> el=words(W)
6-element Array{Array{Int8,1},1}:
 []       
 [2]      
 [1]      
 [2, 1]   
 [1, 2]   
 [1, 2, 1]

julia> T.(el)*permutedims(T.(el))        # multiplication table
6×6 Array{HeckeTElt{Perm{Int16},Int64,Gapjm.Weyl.FiniteCoxeterGroup{Int16,Int64}},2}:
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
         W=coxgroup(:A,n)
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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Hecke.jl#L1-L123' class='documenter-source'>source</a><br>

<a id='Gapjm.Hecke.hecke' href='#Gapjm.Hecke.hecke'>#</a>
**`Gapjm.Hecke.hecke`** &mdash; *Function*.



hecke( W [, parameter, [rootparameter]] ) return a Hecke algebra for W

**Example**

```julia-repl
julia> W=coxgroup(:B,2)
W(B₂)

julia> Pol(:q)
q

julia> H=hecke(W,q)
Hecke(W(B₂),q)

julia> H.para
2-element Array{Array{Pol{Int64},1},1}:
 [q, -1]
 [q, -1]

julia> H=hecke(W,q^2,rootpara=q)
Hecke(W(B₂),q²,rootpara=q)

julia> [H.para,rootpara(H)]
2-element Array{Array{T,1} where T,1}:
 Array{Pol{Int64},1}[[q², -1], [q², -1]]
 Pol{Int64}[q, q]                                  

julia> H=hecke(W,[q^2,q^4],rootpara=[q,q^2])
Hecke(W(B₂),Pol{Int64}[q², q⁴],rootpara=Pol{Int64}[q, q²])

julia> [H.para,rootpara(H)]
2-element Array{Array{T,1} where T,1}:
 Array{Pol{Int64},1}[[q², -1], [q⁴, -1]]
 Pol{Int64}[q, q²]

julia> H=hecke(W,9,rootpara=3)
Hecke(W(B₂),9,rootpara=3)

julia> [H.para,rootpara(H)]
2-element Array{Array{T,1} where T,1}:
 Array{Int64,1}[[9, -1], [9, -1]]
 [3, 3]                              
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Hecke.jl#L134-L177' class='documenter-source'>source</a><br>


<a id='KL.jl-Documentation-1'></a>

# KL.jl Documentation

<a id='Gapjm.KL' href='#Gapjm.KL'>#</a>
**`Gapjm.KL`** &mdash; *Module*.



This  module ports Chevie functionality for Kazhdan-Lusztig polynomials and bases.

Let  `ℋ` be  the Iwahori-Hecke algebra  of a Coxeter  system `(W,S)`, with quadratic  relations `(Tₛ-uₛ₀)(Tₛ-uₛ₁)=0`  for `s∈  S`. If  `-uₛ₀uₛ₁` has a square  root  `wₛ`,  we  can  scale  the  basis  `Tₛ`  to  get  a new basis `tₛ=-Tₛ/wₛ`    with   quadratic    relations   `(tₛ-vₛ)(tₛ+vₛ⁻¹)=0`   where `vₛ=wₛ/uₛ₁`.   The  most  general  case   when  Kazhdan-Lusztig  bases  and polynomials  can be defined is when the parameters `vₛ` belong to a totally ordered  abelian group `Γ`  for multiplication, see  Lus83. We set `Γ⁺= {γ∈ Γ∣γ>0}` and `Γ⁻={γ⁻¹∣γ∈ Γ⁺}={γ∈ Γ∣γ<0}`.

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

julia> @btime test_kl(coxgroup(:F,4));
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

julia>@btime test_kl2(coxgroup(:F,4));
  8s (97455915 allocations: 6.79 GiB)
```

Compare to GAP3 where the following function takes 42s for F4

```
test_kl2:=function(W)local el;
  el:=Elements(W);
  List(el,x->List(el,y->KazhdanLusztigPolynomial(W,x,y)));
end;
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/KL.jl#L22-L115' class='documenter-source'>source</a><br>

<a id='Gapjm.KL.KLPol' href='#Gapjm.KL.KLPol'>#</a>
**`Gapjm.KL.KLPol`** &mdash; *Function*.



KLPol(W,y,w) returns the Kazhdan-Lusztig polynomial P_{y,w} of W

To  compute Kazhdan-Lusztig polynomials in  the one-parameter case it seems that  the best  approach still  is by  using the  recursion formula  in the original  article KL79. One can first run  a number of standard checks on a given  pair  of  elements  to  see  if the computation of the corresponding polynomial  can be reduced to a similar computation for elements of smaller length. One such check involves the notion of critical pairs (cf. Alv87): a pair  of elements `w₁,w₂∈  W` such that  `w₁≤w₂` is *critical*  if `ℒ(w₂) ⊆ ℒ(w₁)`  and `ℛ (w₂)⊆ ℛ (w₁)`, where `ℒ`  and `ℛ` denote the left and right descent  set, respectively.  Now if  `y≤w ∈  W` are arbitrary elements then there   always  exists  a  critical  pair   `z≤w`  with  `y≤z≤w`  and  then `P_{y,w}=P_{z,w}`.  Given two elements `y` and `w`, such a critical pair is found by the function 'CriticalPair'. Whenever the polynomial corresponding to a critical pair is computed then this pair and the polynomial are stored in the property `:klpol` of the underlying Coxeter group.

```julia-repl
julia> W=coxgroup(:B,3)
W(B₃)

julia> map(i->map(x->KLPol(W,one(W),x),elements(W,i)),1:W.N)
9-element Array{Array{Pol{Int64},1},1}:
 [1, 1, 1]                       
 [1, 1, 1, 1, 1]                 
 [1, 1, 1, 1, 1, 1, 1]           
 [1, 1, 1, x+1, 1, 1, 1, 1]      
 [x+1, 1, 1, x+1, x+1, 1, x+1, 1]
 [1, x+1, 1, x+1, x+1, x²+1, 1]  
 [x+1, x+1, x²+x+1, 1, 1]        
 [x²+1, x+1, 1]                  
 [1]
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/KL.jl#L186-L220' class='documenter-source'>source</a><br>

<a id='Gapjm.Hecke.Tbasis' href='#Gapjm.Hecke.Tbasis'>#</a>
**`Gapjm.Hecke.Tbasis`** &mdash; *Function*.



```julia-repl
julia> W=coxgroup(:B,3)
W(B₃)

julia> Pol(:v);H=hecke(W,v^2,rootpara=v)
Hecke(W(B₃),v²,rootpara=v)

julia> C=Cpbasis(H)
(::getfield(Gapjm.KL, Symbol("#f#10")){Pol{Int64},Perm{Int16},HeckeAlgebra{Pol{Int64},Gapjm.Weyl.FiniteCoxeterGroup{Int16,Int64}}}) (generic function with 3 methods)

julia> T=Tbasis(H)
(::getfield(Gapjm.Hecke, Symbol("#f#32")){Pol{Int64},Perm{Int16},HeckeAlgebra{Pol{Int64},Gapjm.Weyl.FiniteCoxeterGroup{Int16,Int64}}}) (generic function with 4 methods)

julia> T(C(1,2))
v⁻²T.+v⁻²T₂+v⁻²T₁+v⁻²T₁₂
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/KL.jl#L307-L324' class='documenter-source'>source</a><br>


<a id='Garside.jl-Documentation-1'></a>

# Garside.jl Documentation

<a id='Gapjm.Garside' href='#Gapjm.Garside'>#</a>
**`Gapjm.Garside`** &mdash; *Module*.



Garside  monoids are a general class  of monoids whose most famous examples are  the braid  and dual  braid monoids.  The implementation  of these last monoids is in the framework of a general implementation of Garside monoids.

To   define  them  we  first  need   to  introduce  some  vocabulary  about divisibility  in monoids. A *left divisor* of  `x` is a `d` such that there exists  `y` with `x=dy` (and then we say  that `x` is a *right multiple* of `d`).  We say  that a  monoid `M`  is left  (resp. right) cancellable if an equality `dx=dy` (resp. `xd=yd`) implies `x=y`. The divisor `d` is *proper* if  `y≠1`. We say  that `x` is  an *atom* if  it has no proper left divisor apart from `1`. A *left gcd* of `x` and `y` is a common left divisor `d` of `x`  and `y` such that  any other common left  divisor is a left divisor of `d`. Similarly a *right lcm* of `x` and `y` is a common multiple which is a left divisor of any other common multiple.

We  call *Garside* a monoid `M` which is:

  * left and right cancellable.
  * generated by its atoms,  which are finite in number.
  * such that any element has only finitely many divisors.
  * admits left and  right  gcds  and  lcms.
  * admits a *Garside element*, which is an element `Δ` whose set of left and right divisors coincide and generate `M`.

Garside  elements are not  unique, but there  is a unique  minimal one (for divisibility); we assume such an element has been chosen. Then the divisors of  `Δ` are called the  *simples* of `M`. A  Garside monoid embeds into its group  of fractions, which is called a *Garside group* (a Garside group may have  several distinct Garside structures, as  we will see for Braid groups of finite Coxeter groups).

We  also implement *locally Garside* monoids,  which are monoids where lcms do  not always exist, but exist when any common multiple exists; the set of simples  is then not defined using a  Garside element, but by the condition that  they contain the atoms and are  closed under lcms and taking divisors (see  [BDM01]); since it is not ensured by the existence of `Δ`, one has to add  the condition that  any element is  divisible by finitely many simples (but  the number of simples can be infinite). The main example is the braid monoid of an infinite Coxeter group. It is not known if these monoids embed in  their group of fractions (though that has been proved for braid monoids of Coxeter groups by Paris [Paris01]) and thus computing in the monoid does not help for computing in the group (only the monoid is implemented here).

What allows computing with Garside and locally Garside monoids, and Garside groups,  is the fact  that they admit  normal forms –-  these normal forms where  first exhibited for  braid monoids by  Deligne [Del72], who extended previous work of Brieskorn, Saito [BS72] and Garside [Gar69]:

(i) Let  `M` be a locally Garside monoid and let `b∈ M`. Then there is a unique maximal left simple divisor `α(b)` of `b`, called the *head* of `b` –- any other simple dividing `b` on the left divides `α(b)` on the left.

(ii) Assume  `M` is a Garside monoid, `Δ` is  its Garside element and `G` is its group  of fractions. Then,  given any element  `x∈ G`, there  is some power `Δⁱ` such that `Δⁱ x∈ M`.

A consequence of (i) is that any element has a canonical decomposition as a product of simples, called its left-greedy normal form. If we define `ω(x)` by  `x=α(x)ω(x)`, then the normal form of `x` is `α(x)α(ω(x))α(ω^2(x))…` We use  the normal form to represent elements  of `M`, and when `M` is Garside (ii)  to represent elements  of `G`: given  `x∈ G` we  compute the smallest power  `i`  such  that  `Δⁱ  x∈  M`,  and  we  represent  `x` by the couple `(i,Δ⁻ⁱx)`.  We are thus reduced to the case where `x∈ M`, not divisible by `Δ`,  where we represent  `x` by the  sequence of simples which constitutes its normal form. We now describe Artin-Tits braid monoids. Let `(W,S)` be a Coxeter system, that is `W` has presentation

`⟨s∈ S∣s^2=1, sts⋯ =tst⋯   (mₛₜ factors on each side) for s,t∈ S⟩`

for  some Coxeter matrix `mₛₜ` for `s,t∈ S`. The braid group `B` associated to `(W,S)` is the group defined by the presentation

`⟨𝐬∈ 𝐒∣ 𝐬𝐭𝐬⋯ =𝐭𝐬𝐭⋯  (mₛₜ factors on each side) for 𝐬,𝐭∈ 𝐒⟩`

The *positive* braid monoid `B⁺` associated to `W` is the monoid defined by the  presentation above –- it identifies to the submonoid of `B` generated by  `𝐒` by  the result  of Paris  mentioned above.  This monoid  is locally Garside,  with set of simples  in bijection with elements  of `W` and atoms the elements of `𝐒`; we will denote by `𝐖` the set of simples, and by `𝐰 ↦ w`  the bijection between simples and elements  of `W`. The group `W` has a length  defined  in  terms  of  reduced expressions. Similarly, having only homogeneous relations, `B⁺` has a natural length function. Then `𝐖` can be characterized  as the subset of the elements  of `B⁺` of the same length as their image in `W`.

If  `W` is finite, then `B⁺` is Garside with Garside element the element of `𝐖` whose image is the longest element  of `W`. A finite Coxeter group is also  a reflection group in  a real vector space,  thus in its complexified `V`,  and `B` has also a topological definition as the fundamental group of the  space `Vʳᵉᵍ/W`, where `Vʳᵉᵍ`  is the set of  elements of `V` which are fixed  by no  non-identity element  of `S`;  however, we  will not use this here.

Given a Coxeter group `W`,

**Examples**

```julia-repl
julia> W=coxgroup(:A,4)
W(A₄)

julia> B=BraidMonoid(W)
BraidMonoid(W(A₄))
```

constructs  the  associated  braid  monoid,  and  then  as  a  function 'B' constructs  elements of the braid monoid (or group when `W` is finite) from a list of generators.

```julia-repl
julia> w=B(1,2,3,4)
1234

julia> w^3
121321432.343

julia> word(W,α(w^3))
9-element Array{Int64,1}:
 1
 2
 1
 3
 2
 1
 4
 3
 2

julia> w^4
δ.232432

julia> inv(w)
(1234)⁻¹
```

As  seen in the fourth  line above, the function  'α(b)' returns the simple `α(b)∈ 𝐖` as an element of `W`.

How  an element of  a Garside group  is printed is  controlled by IOcontext ':greedy'.  By default, elements are printed  as fractions `a⁻¹b` where `a` and  `b` have no left common divisor. Each  of `a` and `b` is printed using its left-greedy normal form, that is a maximal power of the Garside element followed  the rest.  One can  print the  entire element  in the left-greedy normal  from by setting the ':greedy' IOContext; with the same `w` as above we have:

```julia-repl
julia> repr(w^-1,context=IOContext(stdout,:greedy=>true,:limit=>true))
"δ⁻¹.232432"
```

Finally,  `repr` gives   `w`  back   in  a   form  which   after  assigning 'B=braid(braid_monoid(W))' can be input back into Julia:

```julia-repl
julia> repr(w)
"B(1,2,3,4)"

julia> repr(w^3)
"B(1,2,1,3,2,1,4,3,2,3,4,3)"

julia> repr(w^-1)
"B(-4,-3,-2,-1)"
```

In  general elements of  a Garside monoid  are displayed thus  as a list of their constituting atoms.

We  now describe the dual braid monoid.  For that, we first give a possible approach  to construct Garside monoids. Given a  group `W` and a set `S` of generators  of `W` as a monoid, we  define the length `l(w)` as the minimum number of elements of `S` needed to write `w`. We then define left divisors of   `x`  as  the  `d`   such  that  there  exists   `y`  with  `x=dy`  and `l(d)+l(y)=l(x)`.  We say that `w∈ W` is  balanced if its set of left and right  divisors coincide,  is a  lattice (where  upper and lower bounds are lcms and gcds) and generates `W`. Then we have:

suppose `w` is balanced and let `[1,w]` be its set of divisors (an interval for  the partial order  defined by divisibility).  Then the monoid `M` with generators  `[1,w]` and relations  `xy=z` whenever `xy=z`  holds in `W` and `l(x)+l(y)=l(z)`  is Garside,  with simples  `[1,w]` and  atoms `S`.  It is called the interval monoid defined by the interval `[1,w]`.

The  Artin-Tits braid monoid  is an interval  monoid by taking  for `S` the Coxeter generators, in which case `l` is the Coxeter length, and taking for `w`  the longest element of `W`. The dual monoid, constructed by Birman, Ko and  Lee  for  type  `A`  and  by  Bessis  for  all  well-generated complex reflection  groups, is obtained in  a similar way, by  taking this time for `S`  the set of all reflections, and for `w` a Coxeter element; then `l` is the  reflection length  'reflength' (for  well-generated complex reflection groups  whihc are not real `S` contains only those reflections which divide `w`  for the  reflection length);  for the  dual monoid  the simples are of cardinality  the  generalized  Catalan  numbers.  An  interval  monoid  has naturally  an inverse morphism from `M`  to `W`, called 'EltBraid' which is the  quotient map from the  interval monoid to `W`  which sends back simple braids to `[1,w]`. #A  last notable  notion is  *reversible* monoids.  Since in CHEVIE we store #only  left normal forms, it is easy to compute left lcms and gcds, but hard #to  compute right ones.  But this becomes  easy to do  if the monoid has an #operation  'reverse', which has the property that  'a' is a left divisor of #'b'  if and only if  'reverse(a)' is a right  divisor of 'reverse(b)'. This #holds  for Artin-Tits  and dual  braid monoids;  Artin-Tits monoids  have a #reverse  operation which consists of reversing a word, written as a list of #atoms.  The dual monoid  also has a  reverse operation defined  in the same #way,  but this operation changes  monoid: it goes from  the dual monoid for #the  Coxeter element `w` to the dual  monoid for the Coxeter element `w⁻¹`. #The  operations 'RightLcm' and  'RightGcd', as well  quite a few algorithms #have faster implementations if the monoid has a reverse operation.

This module implements functions to solve the conjugacy problem and compute centralizers  in Garside groups, following the work of Franco, Gebhardt and Gonzalez-Meneses.

Two  elements `w` and `w'` of a monoid  `M` are *conjugate* in `M` if there exists  `x∈ M` such that `wx=xw'`; if  `M` satisfies the Öre conditions, it has  a  group  of  fractions  where  this  becomes  `x⁻¹wx=w'`,  the  usual definition  of conjugacy. A special case  which is even closer to conjugacy in  the group is if there exists `y∈  M` such that `w=xy` and `w'=yx`. This relation  is not transitive in general,  but we call *cyclic conjugacy* the transitive closure of this relation, a restricted form of conjugacy.

The  next  observation  is  that  if  `w,w'`  are conjugate in the group of fractions  of the Garside monoid `M` then  they are conjugate in `M`, since if  `wx=xw'` then  there is  a power  `Δⁱ` which  is central and such that `xΔⁱ∈ M`. Then `wxΔⁱ=xΔⁱ w'` is a conjugation in `M`.

The  crucial observation for solving the  conjugacy problem is to introduce `inf(w):=sup{i such  that  Δⁱ   divides  w}`  and `sup(w):=inf{i such  that  w  divides Δⁱ}`, and to notice  that the number of  conjugates of `w` with  same `inf` and `sup` as `w`  is finite. Further, a  theorem of Birman shows  that the maximum `inf` and  minimum `sup` in a conjugacy class can be achieved simultaneously; the elements  achieving this are called the super summit set of `w`. Thus a way to  determine if two elements are conjugate  is to find a representative of both  of them in  their super summit  set, and then  solve conjugacy within that  set. This can also be used  to compute the centralizer of an element: if  we consider  the super  summit set  as the  objects of a category whose morphisms are the conjugations by simple elements, the centralizer is given by the endomorphisms of the given object.

We illustrate this on an example:

```julia-repl
julia> b=B(2,1,4,1,4)
214.14

julia> c=B(1,4,1,4,3)
14.143

julia> d=representative_operation(b,c)
(1)⁻¹21321432

julia> b^d
14.143

julia> centralizer_generators(b)
3-element Array{Gapjm.Garside.GarsideElm{Perm{Int16},BraidMonoid{Perm{Int16},Gapjm.Weyl.FiniteCoxeterGroup{Int16,Int64}}},1}:
 21.1         
 321432.213243
 4

julia> C=conjcat(b,:ss)
category with 10 objects and 32 maps

julia> C.obj
10-element Array{Gapjm.Garside.GarsideElm{Perm{Int16},BraidMonoid{Perm{Int16},Gapjm.Weyl.FiniteCoxeterGroup{Int16,Int64}}},1}:
 214.14
 14.124
 143.13
 1214.4
 13.134
 124.24
 1343.1
 24.214
 134.14
 14.143
```

There  is a faster  solution to the  conjugacy problem given in [gebgon10]: for  each `b∈ M`, they define a  particular simple left divisor of `b`, its *preferred  prefix*  such  that  the  operation  *sliding* which cyclically conjugates  `b` by  its preferred  prefix, is  eventually periodic, and the period  is contained in the super summit set  of `x`. We say that `x` is in its  sliding circuit if some  iterated sliding of `x`  is equal to `x`. The set  of sliding  circuits in  a given  conjugacy class  is smaller than the super  summit  set,  thus  allows  to  solve  the conjugacy problem faster. Continuing from the above example,

```julia-repl
julia> word(W,preferred_prefix(b))
2-element Array{Int64,1}:
 2
 1

julia> b^B(preferred_prefix(b))
1214.4

julia> b1=b^B(preferred_prefix(b))
1214.4

julia> C=conjcat(b)
category with 3 objects and 7 maps

julia> C.obj
3-element Array{Gapjm.Garside.GarsideElm{Perm{Int16},BraidMonoid{Perm{Int16},Gapjm.Weyl.FiniteCoxeterGroup{Int16,Int64}}},1}:
 214.14
 1214.4
 1343.1
```

Finally,  we have implemented  Hao Zheng's algorithm  to extract roots in a Garside monoid:

```julia-repl
julia> W=coxgroup(:A,3)
W(A₃)

julia> B=BraidMonoid(W)
BraidMonoid(W(A₃))

julia> pi=B(B.delta)^2
δ²

julia> root(pi,2)
δ

julia> root(pi,3)
1232

julia> root(pi,4)
132
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Garside.jl#L1-L330' class='documenter-source'>source</a><br>

<a id='Gapjm.Garside.left_divisors' href='#Gapjm.Garside.left_divisors'>#</a>
**`Gapjm.Garside.left_divisors`** &mdash; *Function*.



left_divisors( M, s)

all  the left divisors of the simple element `s` of the Garside monoid `M`, as  a vector  of vectors,  where the  i+1-th vector  holds the  divisors of length i in the atoms.

```julia-repl
julia> W=coxgroup(:A,3)
W(A₃)

julia> B=BraidMonoid(W)
BraidMonoid(W(A₃))

julia> map(x->B.(x),Garside.left_divisors(B,W(1,3,2)))
4-element Array{Array{Gapjm.Garside.GarsideElm{Perm{Int16},BraidMonoid{Perm{Int16},Gapjm.Weyl.FiniteCoxeterGroup{Int16,Int64}}},1},1}:
 [.]   
 [1, 3]
 [13]  
 [132] 

julia> B=DualBraidMonoid(W)
DualBraidMonoid(W(A₃),c=[1, 3, 2])

julia> map(x->B.(x),Garside.left_divisors(B,W(1,3,2)))
4-element Array{Array{Gapjm.Garside.GarsideElm{Perm{Int16},DualBraidMonoid{Perm{Int16},Gapjm.Weyl.FiniteCoxeterGroup{Int16,Int64}}},1},1}:
 [.]                     
 [1, 2, 3, 4, 5, 6]      
 [12, 13, 15, 25, 34, 45]
 [δ]                     
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Garside.jl#L425-L456' class='documenter-source'>source</a><br>

<a id='Gapjm.Garside.DualBraidMonoid' href='#Gapjm.Garside.DualBraidMonoid'>#</a>
**`Gapjm.Garside.DualBraidMonoid`** &mdash; *Type*.



Let  `W` be a well generated complex  reflection group and `c` be a Coxeter element  of `W` (if `W` is a Coxeter group and no `c` is given a particular one  is chosen  by making  the product  of elements  in a  partition of the Coxeter  diagram in two sets where  elements in each commute pairwise). The result  is the dual braid  monoid determined by `W`  and `c`: let `w` be an element  of `W` or a sequence  `s₁,…,sₙ` of integers indices of reflections of `W`.

```julia-repl
julia> W=coxgroup(:A,3)
W(A₃)

julia> B=DualBraidMonoid(W)
DualBraidMonoid(W(A₃),c=[1, 3, 2])

julia> B(2,1,2,1,1)
12.1.1.1

julia> B(-1,-2,-3,1,1)
(25.1)⁻¹1.1
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Garside.jl#L581-L603' class='documenter-source'>source</a><br>

<a id='Gapjm.Garside.fraction' href='#Gapjm.Garside.fraction'>#</a>
**`Gapjm.Garside.fraction`** &mdash; *Function*.



fraction(b) returns a tuple `(x,y)`  of two  elements with  no non-trivial  common left divisor and such that `b=inv(x)*y`.

```julia-repl
julia> B=BraidMonoid(coxgroup(:A,3))
BraidMonoid(W(A₃))

julia> b=B( 2, 1, -3, 1, 1)
(23)⁻¹321.1.1

julia> fraction(b)
(23, 321.1.1)
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Garside.jl#L670-L685' class='documenter-source'>source</a><br>

<a id='Gapjm.word' href='#Gapjm.word'>#</a>
**`Gapjm.word`** &mdash; *Function*.



The Coxeter word for element w of W


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CoxGroups.jl#L182-L184' class='documenter-source'>source</a><br>


word(b) returns  a description  of `b`  as a  list of  the atoms  of which  it is a product.  If `b` is in the Garside group  but not the Garside monoid, it is represented  in  fraction  normal  form  where  as a special convention the inverses  of  the  atoms  are  represented  by  negating  the corresponding integer.

```julia-repl
julia> B=BraidMonoid(coxgroup(:A,3))
BraidMonoid(W(A₃))

julia> b=B(2,1,2,1,1)*inv(B(2,2))
(21)⁻¹1.12.21

julia> word(b)
7-element Array{Int64,1}:
 -1
 -2
  1
  1
  2
  2
  1
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Garside.jl#L703-L728' class='documenter-source'>source</a><br>

<a id='Gapjm.Garside.representative_operation' href='#Gapjm.Garside.representative_operation'>#</a>
**`Gapjm.Garside.representative_operation`** &mdash; *Function*.



representative_operation(b,b1)

The  function returns `a` such that  `b^a=b1` if such exists, and `nothing` otherwise.  If an argument <type> is given,  the computation is done in the corresponding category –- see "conjcat".

```julia-repl
julia> W=coxgroup(:D,4)
W(D₄)

julia> B=BraidMonoid(coxgroup(:D,4))
BraidMonoid(W(D₄))

julia> b=B(2,3,1,2,4,3);b1=B(1,4,3,2,2,2)
1432.2.2

julia> representative_operation(b,b1)
(134312.23)⁻¹

julia> representative_operation(b,b1,:cyc)
232.2
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Garside.jl#L996-L1019' class='documenter-source'>source</a><br>

<a id='Gapjm.Garside.centralizer_generators' href='#Gapjm.Garside.centralizer_generators'>#</a>
**`Gapjm.Garside.centralizer_generators`** &mdash; *Function*.



centralizer_generators(b)

a list of generators  of the centralizer of `b`.  The computation is done by computing the  endomorphisms  of  the  object  <b>  in  the  category  of its sliding circuits.  If an argument <type>  is given, the computation  is done in the corresponding  category –- see "ConjugacySet". The  main use of this is to compute  the  centralizer  in  the  category  of cyclic conjugacy by giving '"Cyc"' as the type.

|    gap> W:=CoxeterGroup("D",4);;     gap> w:=Braid(W)(4,4,4);     4.4.4     gap> CentralizerGenerators(w);     [ 4, 2, (1)^-1.34.431, 34.43, (32431)^-1.132431, 1, (2)^-1.34.432,       (31432)^-1.231432 ]     gap> ShrinkGarsideGeneratingSet(last);     [ 4, 2, 1, 34.43, (3243)^-1.13243 ]     gap> CentralizerGenerators(w,"Cyc");     [ 4 ]     gap> F:=Frobenius(CoxeterCoset(W,(1,2,4)));     function ( arg ) ... end     gap> CentralizerGenerators(w,F);     [ 312343123, 124 ]|


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Garside.jl#L1047-L1072' class='documenter-source'>source</a><br>

<a id='Gapjm.Garside.shrink' href='#Gapjm.Garside.shrink'>#</a>
**`Gapjm.Garside.shrink`** &mdash; *Function*.



shrink(l)

The  list `l` is a  list of  elements of  the same Garside group `G`. This function  tries to find  another set of  generators of the  subgroup of `G` generated by the elements of `l`, of smaller total length (the length being counted  as returned by the function  `word`).

```julia-repl
julia> B=BraidMonoid(coxsym(3))
BraidMonoid(coxsym(3))

julia> b=[B(1)^3,B(2)^3,B(-2,-1,-1,2,2,2,2,1,1,2),B(1,1,1,2)]
4-element Array{Gapjm.Garside.GarsideElm{Perm{UInt8},BraidMonoid{Perm{UInt8},Gapjm.CoxGroups.CoxSymmetricGroup{UInt8}}},1}:
 1.1.1              
 2.2.2              
 (1.12)⁻¹2.2.2.21.12
 1.1.12             

julia> shrink(b)
2-element Array{Gapjm.Garside.GarsideElm{Perm{UInt8},BraidMonoid{Perm{UInt8},Gapjm.CoxGroups.CoxSymmetricGroup{UInt8}}},1}:
 2  
 1  
```


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Garside.jl#L1196-L1220' class='documenter-source'>source</a><br>


<a id='Util.jl-Documentation-1'></a>

# Util.jl Documentation

<a id='Gapjm.Util' href='#Gapjm.Util'>#</a>
**`Gapjm.Util`** &mdash; *Module*.



This  module contains  various utility  functions used  in the  rest of the code.  Maybe some  of them  exist in  some Julia  module I am not aware of; please tell me.

The code is divided in sections  according to semantics.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Util.jl#L1-L7' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Util.jl#L42-L52' class='documenter-source'>source</a><br>


group items of list l according to the values taken by function f on them

```
julia> groupby(iseven,1:10)
Dict{Bool,Array{Int64,1}} with 2 entries:
  false => [1, 3, 5, 7, 9]
  true  => [2, 4, 6, 8, 10]
```

Note:in this version l is required to be non-empty since I do not know how to access the return type of a function


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Util.jl#L61-L71' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.constant' href='#Gapjm.Util.constant'>#</a>
**`Gapjm.Util.constant`** &mdash; *Function*.



whether all elements in list a are equal


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Util.jl#L80' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.blocks' href='#Gapjm.Util.blocks'>#</a>
**`Gapjm.Util.blocks`** &mdash; *Function*.



blocks(M::Matrix)

M  should be a square matrix. Define  a graph G with vertices 1:size(M,1)   and  with an edge between i and j  if either M[i,j] or M[j,i] is not zero   or false. blocks returns a vector of vectors I such that I[1],I[2], etc..   are  the  vertices  in  each  connected  component  of G. In other words,   M[I[1],I[1]],M[I[2],I[2]],etc... are blocks of M.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Util.jl#L85-L93' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.format' href='#Gapjm.Util.format'>#</a>
**`Gapjm.Util.format`** &mdash; *Function*.



format( table; options )

General routine to format a table. Used for character tables.   Options:      row*labels          Labels for rows      column*labels       Labels for columns      rows*label          Label for column of rowLabels      separators          line numbers after which to put a separator      column*repartition  display in pieces of sizes these numbers of cols      rows                show only these rows      columns             show only these columns


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Util.jl#L206-L219' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.prime_residues' href='#Gapjm.Util.prime_residues'>#</a>
**`Gapjm.Util.prime_residues`** &mdash; *Function*.



the numbers less than n and prime to n 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Util.jl#L280' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.phi' href='#Gapjm.Util.phi'>#</a>
**`Gapjm.Util.phi`** &mdash; *Function*.



the Euler function ϕ 


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Util.jl#L299' class='documenter-source'>source</a><br>

<a id='Gapjm.Util.primitiveroot' href='#Gapjm.Util.primitiveroot'>#</a>
**`Gapjm.Util.primitiveroot`** &mdash; *Function*.



primitiveroot(m::Integer) a primitive root mod. m,   that is it generates multiplicatively prime_residues(m).   It exists if m is of the form 4, 2p^a or p^a for p prime>2.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/Util.jl#L305-L309' class='documenter-source'>source</a><br>


!!! warning "Missing docstring."
    Missing docstring for `echelon!`. Check Documenter's build log for details.



<a id='Cycpols.jl-Documentation-1'></a>

# Cycpols.jl Documentation

<a id='Gapjm.CycPols' href='#Gapjm.CycPols'>#</a>
**`Gapjm.CycPols`** &mdash; *Module*.



Cyclotomic  numbers, and cyclotomic polynomials  over the rationals or some cyclotomic  field, play an important role in the study of reductive groups. Special  facilities are provided in this module to deal with them. The type `CycPol` represents the product of a polynomial with a rational fraction in one variable with all poles or zeroes equal to 0 or roots of unity.

The  advantages  of  representing  as  `CycPol`  objects  which  can  be so represented   are:   nice   display   (factorized),  less  storage,  faster multiplication,  division and evaluation. The big drawback is that addition and subtraction are not implemented!

```julia-repl
julia> Pol(:q)
q

julia> p=CycPol(q^18 + q^16 + 2*q^12 + q^8 + q^6)
(q⁸+q⁶-q⁴+q²+1)q⁶Φ₈

julia> p*inv(CycPol(q^2+q+1))
(q⁸+q⁶-q⁴+q²+1)q⁶Φ₃⁻¹Φ₈

```

The variable name in a `CycPol` is set by default to the same as for `Pols`.

`CycPol`s are represented internally by a `struct` with fields:

`.coeff`:  a coefficient, usually a cyclotomic number or a polynomial.

`.valuation`: the valuation in ℤ.

`.v`: a list of pairs `e=>m` of a root of unity `e` and a multiplicity `m`. Here  `e`  is  a  `Root1`,  which  is internally fraction `p//d` with `p<d` representing `E(d)^p`. The pair represents `(q-E(d)^p)^m`.

So if we let `ζ(e)=E(e)=E(d,p)`, a `CycPol` `r` represents

`r.coeff*q^r.valuation*prod(r.vcyc,p->(q-ζ(p[1]))^p[2])`.


<a target='_blank' href='https://github.com/jmichel7/Gapjm.jl/blob/d81b604a199ab41afca2745d094ec58fa0cccae0/src/CycPols.jl#L1-L40' class='documenter-source'>source</a><br>

