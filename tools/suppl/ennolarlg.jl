# Twist a subspets by xi, an element of the center of the (untwisted) parent
function EnnolaTwist(HF, xi)
  H=Group(HF)
  W=parent(H)
  w=classreps(W)[position_regular_class(W, Root1(xi))]
  subspets(parent(HF),inclusiongens(H),w*HF.phi)
end

# Predict what should be lusztig_induction_table(HF,WF) by Ennola-twisting that
# of EnnolaTwist(HF,xi^-1);
function CheckRLGByEnnola(HF, W, xi)
  H=Group(HF)
  tHF=EnnolaTwist(HF,inv(xi))
  if charinfo(tHF).charRestrictions!=1:nconjugacy_classes(H) error("hum") end
  if charinfo(HF).charRestrictions!=1:nconjugacy_classes(H) error("hum") end
  WF=spets(W)
  t=lusztig_induction_table(tHF,WF).scalar
  if isnothing(t) println("not applicable"); return nothing end
  if xi^gcd(degrees(H))==1 t=invpermute(t,inv(ennola(H, xi)),dims=2)
  else print("?")
  end
  tp=invpermute(t,inv(ennola(Group(WF),xi)),dims=1)
  fW=fourier(UnipotentCharacters(WF))
  fH=fourier(UnipotentCharacters(HF))
  if !isone(fH*fH') error() end
  t=fW*tp*fH'
  lu=lusztig_induction_table(HF, WF, check=true)
  xprint(lu,"/",tHF," ",xi,"-ennola";parent=false)
  if tp==lu.scalar 
    println(" perfect")
    return
  else println()
  end
  pieces=lu.pieces
  for i in 1:length(pieces)
    p = pieces[i]
    p.ns = tp[p.wnum,p.hnum]
    tp[p.wnum,p.hnum].=0
    s=ratio(p.scalar,p.ns)
    if !isnothing(s)
      if s!=1 println("got piece[$i] upto $s") end
    else
      r=map(i->ratio(p.scalar[:,i],p.ns[:,i]), axes(p.scalar,2))
      if any(isnothing,r) 
        xdisplay(permutedims(p.ns))
        xdisplay(permutedims(p.scalar))
        xprint("!!!","extensions not proportional", r)
        return
      end
      xprintln(fromTeX(rio(),p.identifier)," no.$i/",length(pieces))
      j=filter(i->r[i]!=1,1:length(r))
      for k in j
        xprintln(" char. $k=",fromTeX(rio(),p.ucharnames[k]),
                 " of subgroup extension differs by ",Cyc(r[k]))
      end
    end
  end
  if unique(vcat(tp...))!=[0] error("CheckRLGByEnnola") end
end

# Check that an inner Ennola Eₓᵢ sends RLG to R_Eₓᵢ(L)^Eₓᵢ(G)
function CheckEnnolaRLG(W)
  c=gcd(degrees(W))
  for HF in vcat(map(x->twistings(W,x), parabolic_reps(W))...)
    for xi in E(c).^(1:c-1)
      CheckRLGByEnnola(HF,W,xi)
    end
  end
end
