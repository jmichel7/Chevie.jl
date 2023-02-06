# auto-generated tests from julia-repl docstrings
using Test, Gapjm
#include("../tools/Gap4.jl")
function mytest(file::String,src::String,man::String)
  println(file," ",src)
  omit=src[end]==';'
  src=replace(src,"\\\\"=>"\\")
  exec=repr(MIME("text/plain"),eval(Meta.parse(src)),context=:limit=>true)
  if omit exec="nothing" end
  exec=replace(exec,r" *(\n|$)"s=>s"\1")
  exec=replace(exec,r"\n$"s=>"")
  man=replace(man,r" *(\n|$)"s=>s"\1")
  man=replace(man,r"\n$"s=>"")
  i=1
  while i<=lastindex(exec) && i<=lastindex(man) && exec[i]==man[i]
    i=nextind(exec,i)
  end
  if exec!=man 
    print("exec=$(repr(exec[i:end]))\nmanl=$(repr(man[i:end]))\n")
  end
  exec==man
end
@testset verbose = true "Gapjm" begin
@testset "dSeries.jl" begin
@test mytest("dSeries.jl","W=rootdatum(\"3D4\")","³D₄")
@test mytest("dSeries.jl","l=cuspidal_data(W,3)","2-element Vector{NamedTuple{(:levi, :cuspidal, :d), Tuple{Spets{FiniteCoxeterSubGroup{Perm{Int16},Int64}}, Int64, Root1}}}:\n (levi = ³D₄, cuspidal = 8, d = ζ₃)\n (levi = ³D₄₍₎=Φ₃², cuspidal = 1, d = ζ₃)")
@test mytest("dSeries.jl","Series(W,l[2]...)","ζ₃-series R^³D₄_{³D₄₍₎=Φ₃²}(λ==Id)  H_G(L,λ)==hecke(G₄,Vector{Mvp{Cyc{Int64}, Int64}}[[ζ₃q², ζ₃, ζ₃q]])\n │    γᵩ    φ  ε family #\n─┼────────────────────────\n1│  φ₁‚₀ φ₁‚₀  1        1\n2│  φ₁‚₆ φ₁‚₄  1        2\n3│  φ₂‚₂ φ₁‚₈ -1        5\n6│ φ″₁‚₃ φ₂‚₅  1        4\n5│ φ′₁‚₃ φ₂‚₃ -1        3\n7│  φ₂‚₁ φ₂‚₁ -1        5\n4│³D₄[1] φ₃‚₂  1        5")
@test mytest("dSeries.jl","W=complex_reflection_group(4)","G₄")
@test mytest("dSeries.jl","l=cuspidal_data(W,3)","5-element Vector{NamedTuple{(:levi, :cuspidal, :d), Tuple{Spets{PRSG{Cyc{Rational{Int64}}, Int16}}, Int64, Root1}}}:\n (levi = G₄, cuspidal = 3, d = ζ₃)\n (levi = G₄, cuspidal = 6, d = ζ₃)\n (levi = G₄, cuspidal = 7, d = ζ₃)\n (levi = G₄, cuspidal = 10, d = ζ₃)\n (levi = G₄₍₎=Φ₁Φ′₃, cuspidal = 1, d = ζ₃)")
@test mytest("dSeries.jl","Series(W,l[5]...)","ζ₃-series R^G₄_{G₄₍₎=Φ₁Φ′₃}(λ==Id)  W_G(L,λ)==Z₆\n │   γᵩ φ(mod 3)  ε parameter family #\n─┼─────────────────────────────────────\n1│ φ₁‚₀        1  1      ζ₃q²        1\n5│ φ₂‚₃       ζ₆  1      -ζ₃q        2\n2│ φ₁‚₄       ζ₃ -1        ζ₃        4\n8│ Z₃:2       -1 -1     -ζ₃²q        2\n9│Z₃:11      ζ₃² -1       ζ₃²        4\n4│ φ₂‚₅      ζ₆⁵ -1       -ζ₃        4")
@test mytest("dSeries.jl","cuspidal_data(W,E(3,2))","5-element Vector{NamedTuple{(:levi, :cuspidal, :d), Tuple{Spets{PRSG{Cyc{Rational{Int64}}, Int16}}, Int64, Root1}}}:\n (levi = G₄, cuspidal = 2, d = ζ₃²)\n (levi = G₄, cuspidal = 5, d = ζ₃²)\n (levi = G₄, cuspidal = 7, d = ζ₃²)\n (levi = G₄, cuspidal = 10, d = ζ₃²)\n (levi = G₄₍₎=Φ₁Φ″₃, cuspidal = 1, d = ζ₃²)")
@test mytest("dSeries.jl","dSeries.ennola(rootdatum(\"3D4\"))","SPerm{Int64}: (3,-4)(5,-5)(6,-6)(7,-8)")
@test mytest("dSeries.jl","dSeries.ennola(complex_reflection_group(14))","SPerm{Int64}: (2,43,-14,16,41,34)(3,35,40,18,-11,42)(4,-37,25,-17,-26,-36)(5,-6,-79)(7,-7)(8,-74)(9,-73)(10,-52,13,31,-50,29)(12,53,15,32,-51,-30)(19,71,70,21,67,68,20,69,72)(22,-39,27,-33,-28,-38)(23,24,-66,-23,-24,66)(44,46,49,-44,-46,-49)(45,48,47,-45,-48,-47)(54,-63,-55,-57,62,-56)(58,-65,-59,-61,64,-60)(75,-77)(76,-76)(78,-78)")
@test mytest("dSeries.jl","W=complex_reflection_group(4)","G₄")
@test mytest("dSeries.jl","Series(W,3;proper=true)","1-element Vector{Series}:\n ζ₃-series R^G₄_{G₄₍₎=Φ₁Φ′₃}(λ==Id)  W_G(L,λ)==Z₆")
@test mytest("dSeries.jl","s=Series(W,3,1)[1]","ζ₃-series R^G₄_{G₄₍₎=Φ₁Φ′₃}(λ==Id)  W_G(L,λ)==Z₆\n │   γᵩ φ(mod 3)  ε parameter family #\n─┼─────────────────────────────────────\n1│ φ₁‚₀        1  1      ζ₃q²        1\n5│ φ₂‚₃       ζ₆  1      -ζ₃q        2\n2│ φ₁‚₄       ζ₃ -1        ζ₃        4\n8│ Z₃:2       -1 -1     -ζ₃²q        2\n9│Z₃:11      ζ₃² -1       ζ₃²        4\n4│ φ₂‚₅      ζ₆⁵ -1       -ζ₃        4")
@test mytest("dSeries.jl","s.spets","G₄")
@test mytest("dSeries.jl","s.levi","G₄₍₎=Φ₁Φ′₃")
@test mytest("dSeries.jl","s.cuspidal","1")
@test mytest("dSeries.jl","s.d","Root1: ζ₃")
@test mytest("dSeries.jl","hecke(s)","hecke(G₆‚₁‚₁,Vector{Mvp{Cyc{Int64}, Int64}}[[ζ₃q², -ζ₃q, ζ₃, -ζ₃²q, ζ₃², -ζ₃]])")
@test mytest("dSeries.jl","degree(s)","ζ₃Φ₁Φ₂²Φ″₃Φ₄Φ₆")
@test mytest("dSeries.jl","dSeries.RLG(s)","[G₄]:<φ₁‚₀>-<φ₁‚₄>-<φ₂‚₅>+<φ₂‚₃>-<Z₃:2>-<Z₃:11>")
@test mytest("dSeries.jl","dSeries.char_numbers(s)","6-element Vector{Int64}:\n 1\n 5\n 2\n 8\n 9\n 4")
@test mytest("dSeries.jl","dSeries.eps(s)","6-element Vector{Int64}:\n  1\n  1\n -1\n -1\n -1\n -1")
@test mytest("dSeries.jl","relative_group(s)","G₆‚₁‚₁")
end
end
