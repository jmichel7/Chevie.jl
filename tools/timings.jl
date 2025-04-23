# all timings are on ssi-Precision-Tower-3431
# 5-3-2024: reduced examples so that here is no gc (gc makes irreproduceable)
using Chairmarks
@time using Chevie
#1.8.5 45s (6.10M alloc: 400MB, 0.63% gc, 5.94% compilation 77% recompilation)
#1.9.3 51s (2.26M alloc: 150MiB, 0.10% gc, 1.59% compilation 78% recompilation)
#1.10.2 40.6s (4.00M alloc: 301 MiB, 0.19% gc, 3.12% compilation 9% recompilation)
#1.11.5 134.18s (4.35M alloc: 262.139 MiB, 0.25% gc, 0.93% compilation 59% recompilation)

#julia> @b test_w0(coxsym(7))
#1.9.4   10.091 ms (182390 allocs: 12.981 MiB)
#1.10.3  11.648 ms (182389 allocs: 12.981 MiB)
#1.11.rc 11.671 ms (237323 allocs: 13.295 MiB)
test_w0(W)=Tbasis(hecke(W,Pol()))(longest(W))^2
"""
test_w0:=function(n)local W,T,H,q,w02; # GAP3 takes 75ms for n=7
  W:=CoxeterGroupSymmetricGroup(n);q:=X(Rationals);H:=Hecke(W,q);T:=Basis(H,"T");
  w02:=T(LongestCoxeterWord(W))^2;return Length(w02.elm);
end;
"""

#julia> @b CyclotomicNumbers.testmat(9)^2 seconds=2
#1.9.4  51.760 ms (286779 allocs: 49.168 MiB)
#1.10.3 42.501 ms (286780 allocs: 49.168 MiB)
#1.11.0 34.639 ms (347334 allocs: 47.080 MiB)
#
#testmat(9)^2 in GAP3 35ms in GAP4 60ms
"""
testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
end;
"""

#julia> @b collect(symmetric_group(10)) seconds=4
#1.10.3  240.491 ms (4043780 allocs: 336.001 MiB)
#1.11.rc 189.171 ms (8083258 allocs: 335.999 MiB)
#1.11.5 140.902 ms (8083222 allocs: 335.998 MiB)
#Elements(Group(List([1..9],i->(i,i+1)),()));; GAP3 takes 0.53s

#julia> @b elements(coxgroup(:E,6))
#1.9.4  7.724 ms (61129 allocs: 12.447 MiB)
#1.10.3 4.092 ms (60997 allocs: 12.434 MiB)
#1.11.0 3.531 ms (118188 allocs: 13.205 MiB)
#Elements(CoxeterGroup("E",6));; #GAP3 15.5ms

#julia> @b test_kl(coxgroup(:A,4))
#1.10.3  5.870 ms (155260 allocs: 9.243 MiB)
#1.11.rc 5.558 ms (205766 allocs: 8.983 MiB)
function test_kl(W)
  q=Pol([1],1)
  H=hecke(W,q^2,rootpara=q)
  C=Cpbasis(H)
  T=Tbasis(H)
  [T(C(w)) for w in elements(W)]
end
# GAP3 takes 49ms for A4
"""
test_kl:=function(W)local q,H,T,C;
  q:=X(Rationals);H:=Hecke(W,q^2,q);T:=Basis(H,"T");C:=Basis(H,"C'");
  List(Elements(W),e->T(C(e)));
end;
"""

#julia> @b test_kl1(coxgroup(:A,4))
#1.10.2  24.752 ms (403416 allocs: 43.278 MiB)
#1.11.rc 22.338 ms (702087 allocs: 42.652 MiB)
function test_kl1(W)
  q=Pol([1],1)
  H=hecke(W,q^2,rootpara=q)
  C=Cpbasis(H)
  (C(longest(W))^2).d
end
# GAP3: 220ms for A4
"""
test_kl1:=function(W)local q,H,T,C;
  q:=X(Rationals);H:=Hecke(W,q^2,q);C:=Basis(H,"C'");
  return C(LongestCoxeterElement(W))^2;
end;
"""

#julia> @b test_kl2(coxgroup(:A,4))
#1.10.3  4.526 ms (87990 allocs: 7.199 MiB)
#1.11.rc 3.689 ms (164568 allocs: 6.847 MiB)
function test_kl2(W)
  el=elements(W)
  maximum(degree(KLPol(W,x,y)) for x in el, y in el)
end
# GAP3 takes 51ms for A4
"""
test_kl2:=function(W)local el; el:=Elements(W);
  return Maximum(List(el,x->Maximum(List(el,y->Length(KazhdanLusztigPolynomial(W,x,y))))));
end;
"""

#julia> @b test_b(coxsym(21))
#1.1    316.423 μs (9542 allocations: 583.91 KiB)
#1.3    319.098 μs (9696 allocations: 597.78 KiB)
#1.5    342.150 μs (8880 allocations: 623.09 KiB)
#1.6    319.175 μs (12683 allocations: 643.62 KiB)
#1.7.2  302.635 μs (9007 allocations: 522.72 KiB)
#1.8    271.042 μs (6429 allocations: 428.31 KiB)
#1.8.5  248.055 μs (6433 allocations: 502.81 KiB)
#1.9.0  250.546 μs (6254 allocations: 517.97 KiB)
#1.10.3 264.118 μs (6200 allocs: 517.859 KiB)
#1.11.0 230.321 μs (10200 allocs: 576.328 KiB)
function test_b(W)
  B=BraidMonoid(W)
  b=B(-4,-3,-2,-15,-14,-13,-10,-11,-12,-7,-2,-12,-11,-10,-7,-8,-9,-6,-7,-8,-2,
      -3,-4,-5,-6,7,6,5,4,3,2,8,7,6,9,8,7,10,11,12,2,7,12,11,10,13,14,15,2,3,4)
  b^4
end
# GAP3 takes 3ms
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
#julia> @b test_hm(BigInt,35) seconds=4
#1.9.4  43.756 ms (1780874 allocs: 57.071 MiB)
#1.10.3 42.886 ms (1795660 allocs: 57.618 MiB)
#1.11.0 44.514 ms (1795664 allocs: 57.618 MiB)
#1.11.5 49.851 ms (1795660 allocs: 60.905 MiB)
function test_hm(rtype,rank)
  m=[rtype(1)//(n+m) for n in 1:rank, m in 1:rank]
  one(m)==m*inv(m)
end
# GAP3 for r=35 takes 16ms and GAP4 19ms
"""
test_hm:=function(r)
  r:=List([1..r],n->List([1..r],m->1/(n+m)));
  return r^0=r*r^-1;
end;
"""

#julia> @b PuiseuxPolynomials.fateman(7)
#1.10.3  15.031 ms (112807 allocs: 29.657 MiB)
#1.11.rc 14.864 ms (225638 allocs: 29.637 MiB)
# GAP3: fateman(7) takes 570ms
"""
fateman:=function(n)local p;
  p:=(1+Mvp("x")+Mvp("y")+Mvp("z")+Mvp("t"))^n;
  p:=p*(p+1);
  return Length(p.coeff);
end;
"""

# @b (q=Pol();inv(Frac.([q+1 q+2;q-2 q-3])))
#bad1.6.3 32.148 μs (755 allocations: 65.55 KiB)
#1.7.2  25.765 μs (737 allocations: 49.69 KiB)
#1.8    23.273 μs (621 allocations: 43.44 KiB)
#1.9.3  25.877 μs (662 allocations: 46.00 KiB)
#1.10.3 19.070 μs (429 allocs: 28.953 KiB)
#1.11.rc 15.460 μs (819 allocs: 29.016 KiB)

# @b ((x,y)=Mvp.(:x,:y);inv(Frac.([x+y x-y;x+1 y+1])))
#1.6.3 187.233 μs (3941 allocations: 291.50 KiB)
#1.7.2 139.797 μs (3418 allocations: 236.66 KiB)
#1.8.5 137.125 μs (3102 allocations: 216.19 KiB)
#1.9.0 128.179 μs (2937 allocations: 210.86 KiB)
#1.10.3 153.370 μs (3631 allocs: 243.969 KiB)
#1.11.rc 122.602 μs (4994 allocs: 223.609 KiB)
# [[x+y,x-y],[x+1,y+1]]^-1; # GAP3 3.1ms

# @b CycPols.p(Pol()) #GAP3 1.25 ms
#1.8.5   254.158 μs (5626 allocations: 416.19 KiB)
#1.10.3  261.749 μs (5240 allocs: 394.188 KiB)
#1.11.rc 224.240 μs (7381 allocs: 379.156 KiB)

# @b CycPols.p(1) #GAP3 142μs
#1.8.5   5.140 μs (101 allocations: 19.88 KiB)
#1.10.3  6.162 μs (101.75 allocs: 21.152 KiB)
#1.11.rc 5.819 μs (151 allocs: 21.012 KiB)
u=CycPols.p(Pol());
# @b u _(1) #GAP3 40μs
#1.8.5   24.669 μs (553 allocations: 41.91 KiB)
#1.10.3  27.458 μs (553 allocs: 41.906 KiB)
#1.11.rc 23.535 μs (828 allocs: 41.375 KiB)
# @b u CycPol(_) #GAP3 8.2ms
#1.8.5   4.749 ms (92895 allocations: 7.35 MiB)
#1.10.3  4.906 ms (87044 allocs: 6.861 MiB)
#1.11.rc 3.931 ms (114917 allocs: 6.631 MiB)
