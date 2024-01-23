# all timings are on ssi-Precision-Tower-3431
using BenchmarkTools: @btime
@time using Chevie
#1.8.5 45s (6.10M alloc: 400MB, 0.63% gc, 5.94% compilation 77% recompilation)
#1.9.3 51s (2.26M alloc: 150MiB, 0.10% gc, 1.59% compilation 78% recompilation)
#1.10β2 33s (3.62 M alloc: 232MiB, 0.26% gc, 1.89% compilation 18% recompilation)

#julia> @btime test_w0(7)
#1.3.0 108.103 ms (1789655 allocations: 166.39 MiB)
#1.5.0  86.432 ms (519209 allocations: 125.44 MiB)
#1.6.0 100.107 ms (519465 allocations: 125.43 MiB)
#1.7.2  94.998 ms (1776566 allocations: 127.52 MiB)
#1.8.5  93.780 ms (1776206 allocations: 127.50 MiB)
#1.9.0 118.063 ms (1776192 allocations: 132.71 MiB)
#1.10.0 141.801 ms (1776195 allocations: 132.71 MiB)
function test_w0(n)
  W=coxsym(n+1)
  length(Tbasis(hecke(W,Pol()))(longest(W))^2)
end
# GAP3 is 0.93s
#test_w0:=function(n)local W,T,H,q,w02;
#  W:=CoxeterGroup("A",n);q:=X(Rationals);H:=Hecke(W,q);T:=Basis(H,"T");
#  w02:=T(LongestCoxeterWord(W))^2;return Length(w02.elm);
#end;

#julia> @btime CyclotomicNumbers.testmat(12)^2;
#1.1 253.780 ms (2415862 allocations: 271.86 MiB) before lower
#1.1.1 426.050 ms (6728396 allocations: 500.72 MiB)
#1.5 346.079 ms (4367402 allocations: 366.17 MiB)
#1.6 358.302 ms (4374636 allocations: 366.22 MiB)
#1.7.3 217.190 ms (3045887 allocations: 266.07 MiB) #improved lower
#1.8.5 178.324 ms (2032503 allocations: 198.22 MiB) #improved lower
#1.9.0 174.725 ms (1856571 allocations: 188.89 MiB)
#1.10.0 184.937 ms (1851568 allocations: 188.62 MiB)
#
#testmat(12)^2 in GAP3 best 0.325s in GAP4 0.3s
#testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
#  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
#end;

#julia> @btime test_PermGroup(10)
#1.1   765.888 ms (23421830 allocations: 1.09 GiB)
#1.3   954.585 ms (19793621 allocations: 878.54 MiB)
#1.5   470.366 ms (11317790 allocations: 681.36 MiB)
#1.6   710.814 ms (11305397 allocations: 681.19 MiB)
#1.7   442.578 ms (11302895 allocations: 557.72 MiB)
#1.8   496.321 ms (11302897 allocations: 557.72 MiB)
#1.9.3 415.639 ms (11303721 allocations: 557.64 MiB)
#1.10.0 372.338 ms (4046122 allocations: 336.15 MiB) # inline iterator
test_PermGroup(n::Int)=length(collect(symmetric_group(n)))
# GAP3 test_PermGroup(10);; takes 0.54s
#test_PermGroup:=n->Length(Elements(Group(List([1..n-1],i->(i,i+1)),())));

#julia> @btime test_elm()
#1.1 516.881 ms (5945569 allocations: 1.08 GiB)
#1.3 790.355 ms (5947612 allocations: 1.08 GiB)
#1.5 633.172 ms (3015031 allocations: 1.04 GiB)
#1.6.3 908.766 ms (3015958 allocations: 1.04 GiB)
#1.7.2 494.588 ms (3013378 allocations: 957.91 MiB)
#1.8   466.815 ms (3012044 allocations: 957.88 MiB)
#1.8.5 562.096 ms (2977130 allocations: 956.99 MiB)
#1.9.0 649.148 ms (2949678 allocations: 949.10 MiB)
#1.10β2 542.353 ms (2947233 allocations: 949.10 MiB)
#1.10.0 841.001 ms (2947233 allocations: 948.39 MiB)
test_elm()=length(elements(coxgroup(:E,7)))
# GAP3 1.8sec

#julia> @btime test_kl(coxgroup(:F,4))
#1.1   2.229s (20967283 allocations: 1.71 GiB)
#1.3   2.047s (22599669 allocations: 1.75 GiB)
#1.5   1.842s (15813766 allocations: 2.57 GiB)
#1.6   1.864s (15334891 allocations: 2.57 GiB)
#1.7.2 1.314 s (13936367 allocations: 2.27 GiB)
#1.8.5 1.271 s (13530143 allocations: 2.30 GiB)
#1.9.3 1.346 s (13132355 allocations: 2.27 GiB)
#1.10β2 1.148 s (13127738 allocations: 2.27 GiB)
#1.10.0 1.359 s (13491189 allocations: 2.32 GiB)
function test_kl(W)
  q=Pol([1],1)
  H=hecke(W,q^2,rootpara=q)
  C=Cpbasis(H)
  T=Tbasis(H)
  [T(C(w)) for w in elements(W)]
  "done"
end
# GAP3 takes 11.6s for F4
#test_kl:=function(W)local q,H,T,C;
#  q:=X(Rationals);H:=Hecke(W,q^2,q);T:=Basis(H,"T");C:=Basis(H,"C'");
#  List(Elements(W),e->T(C(e)));
#end;

#julia> @btime test_kl1(coxgroup(:F,4))
#1.2  10.495s (177295210 allocations: 13.27 GiB)
#1.3   9.459s (177295223 allocations: 13.27 GiB)
#1.5   7.269s (66263103 allocations: 11.57 GiB)
#1.6   8.649s (67132341 allocations: 12.17 GiB)
#1.7.2 6.400 s (64255387 allocations: 9.06 GiB)
#1.8   5.896 s (64253717 allocations: 9.06 GiB)
#1.8.5 6.044 s (63855438 allocations: 9.09 GiB)
#1.9.0 6.632 s (63825907 allocations: 9.27 GiB)
#1.10.0 6.814 s (63797136 allocations: 9.28 GiB)
function test_kl1(W)
  q=Pol([1],1)
  H=hecke(W,q^2,rootpara=q)
  C=Cpbasis(H)
  length((C(longest(W))^2).d)
end
# GAP3: 105s

#julia> @btime test_kl2(coxgroup(:F,4))
#1.0.3 8s (97455915 allocations: 6.79 GiB)
#1.3   5.734 s (96339441 allocations: 7.84 GiB)
#1.5   5.715 s (52994491 allocations: 7.22 GiB)
#1.6   5.360 s (49702426 allocations: 6.98 GiB)
#1.7.2 4.329 s (52735960 allocations: 5.92 GiB)
#1.8   4.072 s (52735476 allocations: 5.92 GiB)
#1.9.3 4.307 s (42294917 allocations: 5.24 GiB)
#1.10.0 2.788 s (39520328 allocations: 5.05 GiB)
function test_kl2(W)
  el=elements(W)
  maximum(degree(KLPol(W,x,y)) for x in el, y in el)
end
# GAP3 takes 40 for F4
#test_kl2:=function(W)local el; el:=Elements(W);
#  return Maximum(List(el,x->Maximum(List(el,y->Length(KazhdanLusztigPolynomial(W,x,y))))));
#end;

#julia> @btime test_b(coxsym(21))
#1.1   316.423 μs (9542 allocations: 583.91 KiB)
#1.3   319.098 μs (9696 allocations: 597.78 KiB)
#1.5   342.150 μs (8880 allocations: 623.09 KiB)
#1.6   319.175 μs (12683 allocations: 643.62 KiB)
#1.7.2 302.635 μs (9007 allocations: 522.72 KiB)
#1.8   271.042 μs (6429 allocations: 428.31 KiB)
#1.8.5 248.055 μs (6433 allocations: 502.81 KiB)
#1.9.0 250.546 μs (6254 allocations: 517.97 KiB)
#1.10.0 260.917 μs (6200 allocations: 517.86 KiB)
function test_b(W)
  B=BraidMonoid(W)
  b=B(-4,-3,-2,-15,-14,-13,-10,-11,-12,-7,-2,-12,-11,-10,-7,-8,-9,-6,-7,-8,-2,
      -3,-4,-5,-6,7,6,5,4,3,2,8,7,6,9,8,7,10,11,12,2,7,12,11,10,13,14,15,2,3,4)
  length((b^4).elm)
end
# GAP3 takes 3.8 ms
#test_b:=function()local b,B;
#  B:=Braid(CoxeterGroupSymmetricGroup(21));
#  b:=B(-4,-3,-2,-15,-14,-13,-10,-11,-12,-7,-2,-12,-11,-10,-7,-8,-9,-6,-7,-8,-2,
#     -3,-4,-5,-6,7,6,5,4,3,2,8,7,6,9,8,7,10,11,12,2,7,12,11,10,13,14,15,2,3,4);
#  b:=b^4;
#  return Length(b.elm);
#end;

# for r>13 needs BigInt
#julia> @btime test_hm(BigInt,45)
#1.1.1 1.039 s  (14582778 allocations: 276.85 MiB)
#1.3   682.464 ms (10933914 allocations: 210.82 MiB)
#1.5   548.563 ms (9042037 allocations: 170.42 MiB)
#1.6.3 180.547 ms (3764674 allocations: 121.98 MiB)
#1.7.2 167.951 ms (3764674 allocations: 121.98 MiB)
#1.8   153.937 ms (3764672 allocations: 121.98 MiB)
#1.8.5 140.436 ms (3764672 allocations: 121.98 MiB)
#1.9.0 147.112 ms (3764672 allocations: 121.98 MiB)
#1.10β2 96.456 ms (3789103 allocations: 122.89 MiB)
#1.10.0 142.626 ms (3789103 allocations: 122.89 MiB)
function test_hm(rtype,rank)
  m=[rtype(1)//(n+m) for n in 1:rank, m in 1:rank]
  one(m)==m*inv(m)
end
# GAP3 for r=45 takes 105ms and GAP4 74ms
#test_hm:=function(r)
#  r:=List([1..r],n->List([1..r],m->1/(n+m)));
#  return r^0==r*r^-1;
#end;

#julia> @btime PuiseuxPolynomials.fateman(12)
#1.6   838.576 ms (3363687 allocations: 1.08 GiB)
#1.7.2 574.138 ms (3363690 allocations: 1.03 GiB)
#1.8   550.589 ms (3363690 allocations: 1.03 GiB)
#1.8.5 527.299 ms (3359892 allocations: 1006.48 MiB)
#1.9.0 536.496 ms (3359892 allocations: 1006.48 MiB)
#1.10.0 534.881 ms (3359892 allocations: 1006.48 MiB)
# GAP3: fateman(12) takes 84sec
#fateman:=function(n)local p;
#  p:=(1+Mvp("x")+Mvp("y")+Mvp("z")+Mvp("t"))^n;
#  p:=p*(p+1);
#  return Length(p.coeff);
#end;

# @btime inv(Frac.([q+1 q+2;q-2 q-3])) setup=(q=Pol())
#bad1.6.3 32.148 μs (755 allocations: 65.55 KiB)
#1.7.2 25.765 μs (737 allocations: 49.69 KiB)
#1.8   23.273 μs (621 allocations: 43.44 KiB)
#1.9.3 25.877 μs (662 allocations: 46.00 KiB)
#1.10.0 17.969 μs (421 allocations: 28.72 KiB)

# @btime inv(Frac.([x+y x-y;x+1 y+1])) setup=((x,y)=Mvp(:x,:y))
#1.6.3 187.233 μs (3941 allocations: 291.50 KiB)
#1.7.2 139.797 μs (3418 allocations: 236.66 KiB)
#1.8.5 137.125 μs (3102 allocations: 216.19 KiB)
#1.9.0 128.179 μs (2937 allocations: 210.86 KiB)
#1.10.0 149.729 μs (3620 allocations: 243.48 KiB)
# GAP3 invert Mvp matrix 5.8ms
# [[x+y,x-y],[x+1,y+1]]^-1;
