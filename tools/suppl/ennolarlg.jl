# Twist a subspets by xi, an element of the center of the (untwisted) parent
function EnnolaTwist(HF, xi)
  H=Group(HF)
  W=parent(H)
  w=classreps(W)[position_regular_class(W, Root1(xi))]
  subspets(parent(HF),inclusiongens(H),w*HF.phi)
end

# Predict what should be lusztig_induction_table(HF,WF) by Ennola-twisting that
# of EnnolaTwist(HF,xi^-1);
function PredictRLGByEnnola(HF, WF, xi)
  tHF=EnnolaTwist(HF,inv(xi))
  H=Group(tHF)
  t=lusztig_induction_table(tHF,WF).scalar
  if xi^gcd(degrees(H))==1
    t=invpermute(t,inv(ennola(H, xi)),dims=2)
  end
  ps=ennola(Group(WF), xi)
  invpermute(t,inv(ps),dims=1)
end

function CheckRLGByEnnola(HF, WF, xi)
  H =Group(HF)
  t =PredictRLGByEnnola(HF, WF, xi)
  fW=fourier(UnipotentCharacters(WF))
  fH=conj(fourier(UnipotentCharacters(HF)))
  t =fW*t*fH
  pieces=lusztig_induction_table(HF, WF, check=true)
  xprintln("Checking ", pieces)
  pieces=pieces.pieces
  for i in 1:length(pieces)
    p = pieces[i]
    p.ns = t[p.wnum,p.hnum]
    t[p.wnum,p.hnum].=0
    s=ratio(p.scalar,p.ns)
    if !isnothing(s)
      if s!=1 println("got piece[$i] upto $s") end
    else
      r=map(i->ratio(p.scalar[:,i],p.ns[:,i]), axes(p.scalar,2))
      if any(isnothing,r) error("extensions not proportional", r) end
      println(" *** piece no.$i of ", length(pieces), " ***")
      print(p.identifier)
      j=filter(i->r[i]!=1,1:length(r))
      for k in j
        println("char. $k==",p.ucharnames[k]," of p[:u] extension differs by ",r[k])
      end
    end
  end
  if unique(vcat(t...))!=[0] error("CheckRLGByEnnola") end
end

# Check that an inner Ennola Eₓᵢ sends RLG to R_Eₓᵢ(L)^Eₓᵢ(G)
function CheckEnnolaRLG(W)
  twist(HF,w)=subspets(parent(HF.parent),inclusiongens(Group(HF)),w*HF.phi)
  WF=spets(W)
  c=gcd(degrees(W))
  l=vcat(map(x->twistings(W,x), parabolic_reps(W))...)
  ps=ennola(W)
  uc=charnumbers.(UnipotentCharacters(W).families)
  for xi in 1:c-1
    for p in filter(x->x[1]>=x[2],map(x->[x, EnnolaTwist(x, E(c, xi))],l))
      HF = p[1]
      HFxi = p[2]
      psH = ennola(Group(HF))
      uc = charnumbers.(UnipotentCharacters(Group(HF)).families)
      t=lusztig_induction_table(HF, WF).scalar
      t=invpermute(t,inv(psH),dims=2)
      t=invpermute(t,inv(ps),dims=1)
      txi=lusztig_induction_table(HFxi, WF)
      xprint(E(c,xi),":",HF,"==>",HFxi)
      if txi.scalar!=t
#       p=deepcopy(txi)
#       p.scalar=t
#       cmptbl(p, txi)
#       cmptbl(t, txi.scalar)
        xdisplay(t)
        xdisplay(txi.scalar)
        error("CheckEnnolaRLG")
      else
        println(" OK!")
      end
    end
  end
end
