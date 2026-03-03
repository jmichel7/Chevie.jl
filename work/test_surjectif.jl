## teste si C_B(b) -> C_W(w) est surjectif
function test_surjective(W;ss=Val(:sc))
  res=Vector{Int}[]
  B=BraidMonoid(W)
  ct=CharTable(W)
  for w in elements(W)
    mot=word(W,w)
    # Print(mot,"\n");
    b=B(mot...)
    Cb=Group(image.(centralizer_gens(b;ss)))
    lCw=ct.centralizers[position_class(W,w)]
    if lCw>length(Cb)
      print(div(lCw,length(Cb)),",")
      push!(res,mot)
    end
  end
  println()
  res
end

shortest(W,w)=classinfo(W).classwords[position_class(W,W(w...))]
