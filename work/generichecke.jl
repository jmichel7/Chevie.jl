using Primes, Chevie
"""
`generic_hecke(W,type;power=1)` gives `hecke(W,para)` depending on `type`:
  - 0: [a,b,c]  10: [x0,x1,x2]
  - 1: [1,b,c]  11: [1,x1,x2]
  - 2: [1,q,q^2] if several Hplanes take Lcm
  - 3: [q,E3,E3^2] Spetsial 13: [Pol(),E3,E3^2]
  - 4: [1,E3,E3^2]
  - 5: [1,2,3] primes
  - 6: [a,E(3)*b,E(3)^2*c] 16: [x0,E(3)*x1,E(3)^2*x2]
the variables are raised to `power`
"""
function generic_hecke(W,type;power=1)
  p=1
  vars="xyztuvwabcdefghijklmnopqrs"
  v=0
  function nextvar()
    v+=1
    vars[v]
  end
  para=map(getfield.(hyperplane_orbits(W),:order))do o
    if type>=10 var=nextvar() end
    map(0:o-1)do j
      if type==0 Mvp(Symbol(nextvar()))^power
      elseif type==10 Mvp(Symbol(var,j))^power
      elseif type==1 j==0 ? Mvp(1) : Mvp(Mvp(Symbol(nextvar())))^power
      elseif type==11 j==0 ? Mvp(1) : Mvp(Symbol(var,j))^power
      elseif type==2 Mvp(:q)^(j*div(lcm(ordergens(W)),o)*power)
      elseif type==3 j==0 ? Mvp(:q)^power : E(o,j)
      elseif type==13 j==0 ? Pol()^power : E(o,j)
      elseif type==4 E(o,j)
      elseif type==5 j==0 ? 1 : (p=nextprime(p+1))^power
      elseif type==6 E(o,j)*Mvp(Symbol(nextvar()))^power
      elseif type==16 E(o,j)*Mvp(Symbol(var,j))^power
      end
    end
  end
  hecke(W,Tuple(para))
end

# Check SchurElements(H) satisfy Schur relations
# c is here to be able to make it big(1)
function CheckSchurRelations(H;c=1)
  un=prod(vcat(map(x->one.(x),H.para)...))
  if un isa Mvp
    s=factorized_schur_elements(H)
    Lcm=lcm(s...)
    s=Ref(Lcm).//s
    print("expanding lcm(Sᵪ)/Sᵪ quotients..")
    t=@elapsed s=HeckeAlgebras.expand.(s;c).*un
    Lcm=HeckeAlgebras.expand(Lcm;c)
    println("done in ",t)
  else
    s=schur_elements(H)
    Lcm=s[charinfo(H.W).positionId]
    s=Lcm./s
  end
  ct=CharTable(H)
  ok=0
  t=@elapsed for i in 1:nconjugacy_classes(H.W)
    if any(ismissing,ct.irr[:,i])
      println("# entries missing in $i-th column")
      continue
    end
    p=sum(s.*ct.irr[:,i])
    if i==1 && p!=Lcm println("!!! Sumᵪ χ(1)/Sᵪ not 1")
    elseif i!=1 && !iszero(p)
      println("!!! Sumᵪ χ(",ordinal(i)," class)/Sᵪ not 0");
    else print(".");ok+=1
    end
  end
  println("satisfied $ok/$(length(s)) relations; done in ", t)
end
