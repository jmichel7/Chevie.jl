# all timings are on ssi-Precision-Tower-3431
# 5-3-2024: reduced examples so that here is no gc (gc makes irreproduceable)
using Chairmarks
@time using Chevie
#1.8.5 45s (6.10M alloc: 400MB, 0.63% gc, 5.94% compilation 77% recompilation)
#1.9.3 51s (2.26M alloc: 150MiB, 0.10% gc, 1.59% compilation 78% recompilation)
#1.10.2 40.6s (4.00M alloc: 301 MiB, 0.19% gc, 3.12% compilation 9% recompilation)

#julia> @b test_w0(6)
#1.10.2 13.557 ms (182390 allocs: 12.982 MiB)
function test_w0(n)
  W=coxsym(n+1)
  Tbasis(hecke(W,Pol()))(longest(W))^2
end
"""
test_w0:=function(n)local W,T,H,q,w02; # GAP3 takes 84ms
  W:=CoxeterGroup("A",n);q:=X(Rationals);H:=Hecke(W,q);T:=Basis(H,"T");
  w02:=T(LongestCoxeterWord(W))^2;return Length(w02.elm);
end;
"""

#julia> @b CyclotomicNumbers.testmat(9)^2 seconds=2
#1.10.2 44.561 ms (286780 allocs: 49.168 MiB)
#
#testmat(9)^2 in GAP3 85ms in GAP4 60ms
"""
testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
end;
"""

#julia> @b collect(symmetric_group(10)) seconds=4
#1.10.2 267.295 ms (4043780 allocs: 336.001 MiB)
#Elements(Group(List([1..9],i->(i,i+1)),()));; GAP3 takes 0.56s

#julia> @b elements(coxgroup(:E,6))
#1.10.2 4.282 ms (61016 allocs: 12.435 MiB)
#Elements(CoxeterGroup("E",6));; #GAP3 15.5ms

#julia> @b test_kl(coxgroup(:A,4))
#1.10.2 6.281 ms (155259 allocs: 9.243 MiB)
function test_kl(W)
  q=Pol([1],1)
  H=hecke(W,q^2,rootpara=q)
  C=Cpbasis(H)
  T=Tbasis(H)
  [T(C(w)) for w in elements(W)]
end
# GAP3 takes 57ms for A4
"""
test_kl:=function(W)local q,H,T,C;
  q:=X(Rationals);H:=Hecke(W,q^2,q);T:=Basis(H,"T");C:=Basis(H,"C'");
  List(Elements(W),e->T(C(e)));
end;
"""

#julia> @b test_kl1(coxgroup(:A,4))
#1.10.2 24.752 ms (403416 allocs: 43.278 MiB)
function test_kl1(W)
  q=Pol([1],1)
  H=hecke(W,q^2,rootpara=q)
  C=Cpbasis(H)
  (C(longest(W))^2).d
end
# GAP3: 247ms for A4
"""
test_kl1:=function(W)local q,H,T,C;
  q:=X(Rationals);H:=Hecke(W,q^2,q);C:=Basis(H,"C'");
  return C(LongestCoxeterElement(W))^2;
end;
"""

#julia> @b test_kl2(coxgroup(:A,4))
#1.10.0 4.753 ms (87989 allocs: 7.199 MiB)
function test_kl2(W)
  el=elements(W)
  maximum(degree(KLPol(W,x,y)) for x in el, y in el)
end
# GAP3 takes 57ms for A4
"""
test_kl2:=function(W)local el; el:=Elements(W);
  return Maximum(List(el,x->Maximum(List(el,y->Length(KazhdanLusztigPolynomial(W,x,y))))));
end;
"""

#julia> @b test_b(coxsym(21))
#1.1   316.423 μs (9542 allocations: 583.91 KiB)
#1.3   319.098 μs (9696 allocations: 597.78 KiB)
#1.5   342.150 μs (8880 allocations: 623.09 KiB)
#1.6   319.175 μs (12683 allocations: 643.62 KiB)
#1.7.2 302.635 μs (9007 allocations: 522.72 KiB)
#1.8   271.042 μs (6429 allocations: 428.31 KiB)
#1.8.5 248.055 μs (6433 allocations: 502.81 KiB)
#1.9.0 250.546 μs (6254 allocations: 517.97 KiB)
#1.10.2 283.471 μs (6200 allocs: 517.859 KiB)
function test_b(W)
  B=BraidMonoid(W)
  b=B(-4,-3,-2,-15,-14,-13,-10,-11,-12,-7,-2,-12,-11,-10,-7,-8,-9,-6,-7,-8,-2,
      -3,-4,-5,-6,7,6,5,4,3,2,8,7,6,9,8,7,10,11,12,2,7,12,11,10,13,14,15,2,3,4)
  b^4
end
# GAP3 takes 3.8 ms
"""
test_b:=function()local b,B;
  B:=Braid(CoxeterGroupSymmetricGroup(21));
  b:=B(-4,-3,-2,-15,-14,-13,-10,-11,-12,-7,-2,-12,-11,-10,-7,-8,-9,-6,-7,-8,-2,
     -3,-4,-5,-6,7,6,5,4,3,2,8,7,6,9,8,7,10,11,12,2,7,12,11,10,13,14,15,2,3,4);
  b:=b^4;
  return Length(b.elm);
end;
"""

# for r>13 needs BigInt
#julia> @b test_hm(BigInt,35) seconds=1
#1.10.2 44.782 ms (1795662 allocs: 57.618 MiB)
function test_hm(rtype,rank)
  m=[rtype(1)//(n+m) for n in 1:rank, m in 1:rank]
  one(m)==m*inv(m)
end
# GAP3 for r=35 takes 32ms and GAP4 19ms
"""
test_hm:=function(r)
  r:=List([1..r],n->List([1..r],m->1/(n+m)));
  return r^0=r*r^-1;
end;
"""

#julia> @b PuiseuxPolynomials.fateman(7)
#1.10.2 16.581 ms (112808 allocs: 29.657 MiB)
# GAP3: fateman(7) takes 722ms
"""
fateman:=function(n)local p;
  p:=(1+Mvp("x")+Mvp("y")+Mvp("z")+Mvp("t"))^n;
  p:=p*(p+1);
  return Length(p.coeff);
end;
"""

# @b (q=Pol();inv(Frac.([q+1 q+2;q-2 q-3])))
#bad1.6.3 32.148 μs (755 allocations: 65.55 KiB)
#1.7.2 25.765 μs (737 allocations: 49.69 KiB)
#1.8   23.273 μs (621 allocations: 43.44 KiB)
#1.9.3 25.877 μs (662 allocations: 46.00 KiB)
#1.10.2 19.351 μs (429 allocs: 28.969 KiB)

# @b ((x,y)=Mvp.(:x,:y);inv(Frac.([x+y x-y;x+1 y+1])))
#1.6.3 187.233 μs (3941 allocations: 291.50 KiB)
#1.7.2 139.797 μs (3418 allocations: 236.66 KiB)
#1.8.5 137.125 μs (3102 allocations: 216.19 KiB)
#1.9.0 128.179 μs (2937 allocations: 210.86 KiB)
#1.10.2 161.356 μs (3631 allocs: 243.984 KiB)
# GAP3 invert Mvp matrix 3.9ms
# [[x+y,x-y],[x+1,y+1]]^-1;

#@b CycPols.p(Pol()) #GAP3 1.25 ms
#1.8.5   254.158 μs (5626 allocations: 416.19 KiB)
##1.10.2 272.889 μs (5241 allocs: 394.234 KiB)
#@b CycPols.p(1) #GAP3 142μs
#1.8.5   5.140 μs (101 allocations: 19.88 KiB)
##1.10.2 6.535 μs (101.75 allocs: 21.156 KiB)
u=CycPols.p(Pol());
#@b u _(1) #GAP3 40μs
#1.8.5  24.669 μs (553 allocations: 41.91 KiB)
#1.10.2 28.299 μs (553 allocs: 41.906 KiB)
#@b u CycPol(_) #GAP3 8.2ms
#1.8.5   4.749 ms (92895 allocations: 7.35 MiB)
##1.10.2 5.084 ms (87045 allocs: 6.861 MiB)
