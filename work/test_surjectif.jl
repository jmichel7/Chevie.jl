## teste si C_B(b) -> C_W(w) est surjectif
function test_surjective(W;ss=:sc)
  res=Vector{Int}[]
  B=BraidMonoid(W)
  for w in elements(W)
    mot=word(W,w)
    # Print(mot,"\n");
    b=B(mot...)
    Cb=Group(image.(centralizer_gens(b;ss=ss)))
    Cw=centralizer(W,w)
    if length(Cw)/length(Cb)>1
      print(div(length(Cw),length(Cb)),",")
      push!(res,mot)
    end
  end
  println()
  res
end

minlength(W,w)=length(classinfo(W).classtext[position_class(W,W(w...))])
