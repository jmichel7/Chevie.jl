"""
This module `Murphy.jl` has been ported in december 2020 from

   murphy.g Copyright (C) July 1998  Andrew Mathas 
   mathas@maths.usyd.edu.au University of Sydney

It allows computations with the Murphy basis of an Hecke algebra of type A.

Multiplication  of Murphy basis elements is  done using the Garnir tableaux
as  described in Murphy's paper [Mur1995](biblio.htm#Mur95). This also lets
us convert from the T-basis to the Murphy basis since
    `T_w = M([[1],…,[n]], [[1],…,[n]]) * T_w`.
(We use "M" for the Murphy basis).

As with the T-basis, Murphy basis elements are implemented by `ModuleElts`.
Here the keys are standard tableaux pairs. These are represented by a tuple
`(mu,s,t)` where `mu`, `s` and `t` are all integers and
`H.Murphy.partitions[mu]` is the associated partition and
`H.Murphy.tableaux[mu][s]` (resp `t`) is a named tuple describing `s` (this
tuple is described in the function `initMurphy()` below).

Throughout memory considerations are thrown to the wind as we cache many of
the  more horrible expansions as we go along  in order to save time when we
next need them.
"""
module Murphy

using ..Chevie
export Murphybasis, Spechtmodel

function QuantumInteger(q,n)
  iszero(n) ? zero(q) : n>0 ? sum(i->q^i,0:n-1) : -q^n*QuantumInteger(q,-n)
end

#--------------------- Tableau utilities ----------------------------
# Given two standard tableaux s and t assumed to be of the same shape,
# return as a Coxeter word the permutation w such that <t>=<s>.^w.
function Permtableaux(s,t)
  x=mappingPerm(vcat(s...),vcat(t...))
  w=[]
  while x!=Perm()
    i=1
    while i^x<(i+1)^x i=i+1 end
    push!(w,i)
    x=Perm(i,i+1)*x
  end
  return w
end

function firsttableau(mu) # smallest tableau of shape mu
  d=0
  tab=Vector{Int}[]
  for row  in mu
    push!(tab,(1:row).+d)
    d+=row
  end
  tab
end

function conjugate_tableau(t)
  if isempty(t) return t end
  d=map(x->[x],t[1])
  for r in eachindex(t[1])
    i=2
    while i<=length(t) && r<=length(t[i])
      push!(d[r],t[i][r])
      i+=1
    end
  end
  d
end

istableau(t)=all(x->x==sort(unique(x)),t) && 
             all(x->x==sort(unique(x)),conjugate_tableau(t))

shortTableau(t)=join(map(joindigits,t),"/")

#------------------The Murphy basis code proper---------------------
"""
Set `Murphy.SpechtModules[]=true` to work just inside the Specht modules.
This makes  computations with  the Murphy  basis inside Specht modules
modules much faster but also means that `T` to `Murphy` basis conversions
do not work, even if `Murphy.SpechtModules[]=false` is set later.
"""
const SpechtModules=Ref(false)

struct HeckeMElt{C,TH}<:HeckeElt{TH,C,Tuple{Int,Int,Int}}
  d::ModuleElt{Tuple{Int,Int,Int},C} # has better merge performance than Dict
  H::TH
end

HeckeAlgebras.clone(h::HeckeMElt,d)=HeckeMElt(d,h.H)
Base.zero(::Type{HeckeMElt},H::HeckeAlgebra)=HeckeMElt(zero(ModuleElt{Tuple{Int,Int,Int},coefftype(H)}),H)
Base.zero(h::HeckeMElt)=zero(HeckeMElt,h.H)
HeckeAlgebras.basisname(h::HeckeMElt)="M"

# Create H.Murphy with various components.
# initMurphy(H) is called the first time that the Murphy basis is used.
function initMurphy(H)
  if haskey(H, :Murphy) return end
  t=refltype(H.W)
  if length(t)!=1 || t[1].series!=:A
    error("the Murphy basis is implemented only for irreducible type A")
  end
  mu=fill(1,semisimplerank(H.W)+1) # first partition of n
  H.Murphy=(partitions=[mu],
  # H.tableaux=[partitions][ (ind = index of tableau in H.Tableaux[mu]
  #                         mu = index of tableau shape in H.partition
  #             wd=associated word in S_n (as a CoxeterGroup permutation)
  #                that is perm such that t=t^mu*w where t^mu=firsttableau(mu)
  #             toT=MurphyToT 
  #             Garnir will hold the Garnir expansions) ]
  # Creating a list of all of the standard tableaux is a big overhead for
  # Specht modules of large dimension so we create these tableaux records
  # only as needed, storing them in tableaux[mu][...] with a parallel list
  # Tableaux[mu][...] containing the list of the corresponding tableaux.
  # The bookkeeping to maintain and access this list is done by the
  # function code_tableau().
  tableaux=[[(ind=1,mu=1,wd=one(H.W),Garnir=Dict{Int,Any}(),
                                   toT=Dict{Int,Any}(1=>Tbasis(H)()))]],
  Tableaux=[[firsttableau(mu)]],
  TtoMurphy=Dict{Perm,Any}(one(H.W)=>HeckeMElt(ModuleElt((1,1,1)=>
                                           one(coefftype(H))),H)))
  # T-basis to Murphy-basis conversions; store as computed.
  # Start with identity element.
  InfoChevie("#I Initialized Murphy basis\n")
end

# Given a partition <mu> , which is assumed to be a partition of n, return the 
# index of <mu> in the list H.partitions, or add <mu> to this list if it is not
# already there AND initialize the list H.tableaux[mu] and H.Tableaux[...].
function code_partition(H::HeckeAlgebra, mu)
 t=findfirst(==(mu),H.Murphy.partitions)
  if t!==nothing return t end
  push!(H.Murphy.partitions,deepcopy(mu))
  t=length(H.Murphy.partitions)
  push!(H.Murphy.tableaux,[(ind=1,mu=t,wd=one(H.W),
                    Garnir=Dict{Int,Any}(),toT=Dict{Int,Any}())])
  push!(H.Murphy.Tableaux,[firsttableau(mu)])
  t
end

# given  a tableau  <tab> ,  which is  assumed to  be standard, return the
# corresponding  element of <H>.tableau  (a "CoxeterGroup tableau"). Here,
# <mu> is the index of the shape of <tab> in H.partitions[].
function code_tableau(H::HeckeAlgebra, mu::Integer, tab)
  tabs=H.Murphy.Tableaux[mu]
  ind=findfirst(==(tab),tabs)
  if ind===nothing
    ind=length(tabs)+1
    push!(H.Murphy.tableaux[mu],(ind=ind,mu=mu,wd=H.W(Permtableaux(tabs[1],tab)...),
                                 Garnir=Dict{Int,Any}(),toT=Dict{Int,Any}()))
    push!(tabs,deepcopy(tab))
  end
  H.Murphy.tableaux[mu][ind]
end

CodedTableau{P}=NamedTuple{(:ind, :mu, :wd, :Garnir, :toT),Tuple{Int64,Int64,P,Dict{Int64,Any},Dict{Int64,Any}}}

firsttableau(H::HeckeAlgebra,t::CodedTableau)=H.Murphy.tableaux[t.mu][1]

tableau(H::HeckeAlgebra,t::CodedTableau)=H.Murphy.tableaux[t.mu][t.ind]

firstTableau(H::HeckeAlgebra,t::CodedTableau)=H.Murphy.Tableaux[t.mu][1]

Tableau(H::HeckeAlgebra,t::CodedTableau)=H.Murphy.Tableaux[t.mu][t.ind]

# given a tableau <t> return the basis element x_mu*T_<t> in the T-basis
function Xt(H::HeckeAlgebra,t::CodedTableau)
  get!(tableau(H,t).toT,1)do
    HeckeTElt(if t.ind==1 
      J=vcat(map(x->x[1:end-1],firstTableau(H,t))...)
      WJ=elements(reflection_subgroup(H.W,J))
      ModuleElt(map(x->x=>one(coefftype(H)),WJ))
    else 
      ModuleElt([k*t.wd=>c for (k,c) in Xt(H,firsttableau(H,t)).d])
    end,H)
  end
end

# given a pair (s,t) of tableaux, return the basis element 
# <t>.(<s>.ind)=T_<s>^* x_mu T_<t> in the T-basis.
function MurphyToT(H::HeckeAlgebra, s::CodedTableau, t::CodedTableau)
  get!(tableau(H,t).toT,s.ind)do
    if s.ind==1 Xt(H,t)
    elseif t.ind==1 α(Xt(H, s))
    else Tbasis(H)(inv(s.wd))*Xt(H, t)
    end
  end
end

# convert Murphy basis to T-basis
function HeckeAlgebras.Tbasis(M::HeckeMElt)
  tt=M.H.Murphy.tableaux
  sum([c*MurphyToT(M.H,tt[i][t1],tt[i][t2]) for ((i,t1,t2),c) in M.d])
end

# This function recursively expands T_w into a linear combination of
# Murphy basis elements (using Garnir expansions). Note that we know
# how to write 1 in terms of the Murphy basis. 
# As we go along we cache these expansions in the dict H.TtoMurphy
function TtoMurphy(H::HeckeAlgebra,w)
  get!(H.Murphy.TtoMurphy,w)do
    if SpechtModules[]
      println("\n# WARNING: because Murphy.SpechtModules[]==true the answer\n",
            "# TtoMurphy returns will almost certainly be incorrect.")
    end
    W=H.W
    r=firstleftdescent(W, w^-1)
    res=TtoMurphy(H,w*W(r))*Tbasis(H)(r)
    if w!=w^-1 H.Murphy.TtoMurphy[w^-1]=α(res) end
    res
  end
end

# Murphybasis(H) creates a function which will return a Murphy basis element
# from a pair of tableaux or some HeckeElt.
function Murphybasis(H::HeckeAlgebra)
  initMurphy(H)
  f(h::HeckeElt)=f(Tbasis(h))
  f(h::HeckeTElt)=sum(c*TtoMurphy(h.H,w) for (w,c) in h.d)
  f(h::HeckeMElt)=h
  function f(s,t)
    if !istableau(s) || !istableau(t)
        error("<s> and <t> must be standard tableaux")
    end
    mu=length.(s)
    if mu!=length.(t) error("<s> and <t> must have the same shape") end
    imu=code_partition(H, mu)
    is=code_tableau(H, imu, s)
    it=code_tableau(H, imu, t)
    HeckeMElt(ModuleElt((imu,is.ind,it.ind)=>one(coefftype(H))),H)
  end
end

function Base.show(io::IO, h::HeckeMElt)
  function showbasis(io::IO,(mu,s,t))
    TeX=get(io,:TeX,false)
    repl=get(io,:limit,false)
    H=h.H
    TeXTableau(tab)="\\tab("*join(map(join,tab),",")*")"
    StringTableau(t)=repl ? shortTableau(t) :
      "["*join(map(p->"["*join(p,",")*"]",t),",")*"]"
      u=H.Murphy.Tableaux[mu][t]
    if SpechtModules[] 
      TeX ? TeXTableau(u)*"\n" : string("S(",StringTableau(u),")")
    else 
     t=H.Murphy.Tableaux[mu][s]
      if TeX string(HeckeAlgebras.basisname(h),"(",join(TeXTableau.((t,u)),", "),")\n")
      else string(HeckeAlgebras.basisname(h),"(",StringTableau(t),", ",StringTableau(u),")")
      end
    end
  end
  show(IOContext(io,:showbasis=>showbasis),h.d)
end

# returns the position of the integer <i> in the tableau record <t>.
function PositionInTableau(H, t::CodedTableau, i)
  tab=Tableau(H,t)
  for row in eachindex(tab)
    col=findfirst(==(i),tab[row])
    if col!==nothing return (row=row,col=col) end
  end
  error(i," is not contained in the tableau ",shortTableau(tab))
end

function Base.:*(m::HeckeMElt, h::HeckeTElt)
  H=h.H
  if H!=m.H error("not elements of the same algebra") end
  W = H.W
  q = H.para[1][1]
  # mh=return value, initially zero with respect to the Murphy basis
  mh=zero(HeckeMElt,H)
  for (w,hcoeff) in h.d
    mw=m
    for r in word(W,w) # multiply one simple reflection at a time
      mr=zero(m.d).d
      for ((mu,is,it),coeff) in mw.d
       t=H.Murphy.tableaux[mu][it]
        # we are interested in the two nodes, <nodeR> and <nodeS> 
        # which are swapped by the transposition r=(r,r+1). Thus,
        # these are the nodeRs such that t<nodeR>=r and t<nodeS>=r+1.
        nodeR = PositionInTableau(H,t,r)
        nodeS = PositionInTableau(H,t,r+1)
        if nodeR.row == nodeS.row
          push!(mr,(mu,is,it)=>q*coeff)
        elseif nodeR.col!=nodeS.col# then t*r is still standard
         tabr=deepcopy(H.Murphy.Tableaux[mu][it])
          tabr[nodeR.row][nodeR.col]=r+1
          tabr[nodeS.row][nodeS.col]=r
          tr=code_tableau(H, mu, tabr)
          if nodeR.row < nodeS.row # up in the Bruhat order
            push!(mr,(mu,is,tr.ind)=>coeff)
          else
            push!(mr,(mu,is,tr.ind)=>q*coeff)
            push!(mr,(mu,is,it)=>(q-1)*coeff)
          end
        else  # The hard part: here in the tableau t, r+1 occupies the nodeR 
              # below r; so nodeS.row=nodeR.row+1 and nodeS.col=nodeR.col and
              # interchanging them gives a non-standard tableau.
              s=H.Murphy.tableaux[mu][is]
              append!(mr,(coeff*GarnirExpansion(H, nodeR, s, t)).d.d)
        end
      end
      mw=HeckeMElt(ModuleElt(mr),H)
    end
    mh+=hcoeff*mw
  end
  mh
end

Base.:*(h::HeckeTElt,m::HeckeMElt)=α(α(m)*α(h))

# The  real  work:  <node>=(r,c)  is  the  coordinate where the tableau <t>
# becomes  non-standard; ie. if <t>(r,c)=x  then <t>(r+1,c)=x+1 and we want
# to  expand <t>*(x,x+1). We are actually expanding M_{<s>,<t>}*T_{(x,x+1)}
# whexe  we know  that <t>*(x,x+1)  is not  standard. Once  we know  how to
# expxess  M_{1,<t>}*T_x in the  Murphy basis we  store this as t.Garnir[x]
# for future reference.
function GarnirExpansion(H::HeckeAlgebra,node,s::CodedTableau,t::CodedTableau)
  rt=Tableau(H,t)[node.row][node.col]
  if !haskey(t.Garnir,rt)
    W = H.W
    # a typical situation here is
    #         1 2 7
    #   <t> = 3 5 8     with <node>=[2,2];
    #         4 6 9
    # thus, we want to expand the tableau
    #         1 2 7 
    #    t' = 3 6 8 = <t>*(5,6) [note that 5 occupies position [2,2] in <t>]
    #         4 5 9
    # into a linear combination of standard tableaux. To do this we first
    # pretend that we started with
    #         1 2 3
    #    g  = 4 6 8
    #         5 7 9
    # (that  is  we  put  the  numbers  in  order  upto, but not including
    # [node.row,node.col]  and then enter them  in order starting from the
    # next  row  down,  filling  up  the  nodes  around  <node>  and  then
    # continuing on. ie. what almost Murphy called a Garnir tableau).
    gtab=deepcopy(firstTableau(H,t))
    a=gtab[node.row][node.col] # first number being moved; above a=5 and b=8
    b=gtab[node.row+1][node.col] # last number being moved
    gtab[node.row][node.col+1:end]=a+node.col+1:b
    gtab[node.row+1][1:node.col-1]=a:a+node.col-2
    gtab[node.row][node.col]  =a+node.col-1
    gtab[node.row+1][node.col]=a+node.col
    g=code_tableau(H, t.mu, gtab)
    # w is the permutation such that t=g*w => T_t = T_g*T_w
    w=Permtableaux(gtab,Tableau(H,t))
    rg=Tableau(H,g)[node.row][node.col]
    if !haskey(tableau(H,g).Garnir,rg)
      # first note that, by an astute look at right coset sums,
      #      1 2 3   1 2 3   1 2 3   1 2 3   1 2 3         1 2 3       1 2 3
      # (*)  4 5 6 + 4 5 7 + 4 5 8 + 4 6 7 + 4 6 8 + ... + 4 7 8 = h * 4
      #      7 8 9   6 8 9   6 7 9   5 8 9   5 7 9         5 6 9       5 6 7 8
      #                                                                9
      # Because of our choice of g all of the LH tableaux are standard, except
      # the last  and the term on the RHS. We'll worry about the RHS later.
      # First we spin out the tableaux on the left hand side.
      mres=zero(HeckeMElt,H).d.d
      for J in combinations(a:b, node.col)
        if J != a:a+node.col-1
          gtab[node.row+1][1:node.col]=J
          gtab[node.row][node.col:end]=sort(setdiff(a:b,J))
          # note that we set <s>=t^mu below; this is because we later
          # have to multiply by T_s^*
          push!(mres,(g.mu,1,code_tableau(H,g.mu,gtab).ind)=>-1)
        end
      end
      tableau(H,g).Garnir[rg]=HeckeMElt(ModuleElt(mres),H)
      # Next, if Murphy.SpechtModules[]=false (in which case we 
      # just work in the Specht module), we look after the right hand term
      # in (*) above. In general it won't correspond to a partition but we
      # can find a partition <nu> and a <d> in <W> such that 
      # T_d<RHS>=<x_nu>T_<d>
      if SpechtModules[]==false
        tab=gtab[1:node.row-1]
        if node.col > 1 push!(tab, gtab[node.row][1:node.col - 1]) end
        push!(tab,a:b)
        if node.col<length(gtab[node.row+1])
          push!(tab,gtab[node.row+1][node.col+1:end])
        end
        if length(gtab)>=node.row+2
          append!(tab, gtab[node.row+2:end])
        end
        # Now we reorder <tab> so that the diagram has the shape of a
        # partition. Our <tab> above becomes [[5,6,7,8],[1,2,3],[4],[9]].
        sort!(tab, by=x->(-length(x),x[1]))
        # which gives us the (shape of the) new tableau
        tnu = H.Murphy.tableaux[code_partition(H, map(length, tab))][1]
        # and finally we have <d>. The point is that tab = T_d^-1*tnu*T_d.
        d=W(Permtableaux(Tableau(H,tnu), tab)...)
        # <tab> is now under control, but we still need to compute <h>
        # from (*). The point here is that we are essentially writing
        #             $H.I = \bigcup H.d = \bigcup d' I$
        # where H and I are two sugroups and d and d' run over coset
        # representatives of H and I.
        gtab=firstTableau(H,g)
        J=gtab[node.row][1:end-1]
        append!(J,gtab[node.row+1][1:end-1])
        K=setdiff(J,[a-1,b])
        h=map(inv,vcat(reduced(reflection_subgroup(W,K),
                               reflection_subgroup(W,J))...))
        h=HeckeTElt(ModuleElt([w=>one(coefftype(H)) for w in h];check=false),H)
        # the multiplication below is quite costly as it is recursive; 
        # but it is only done once as we store the result in g.Garnir.
        tab1=firstTableau(H,tnu)
        tableau(H,g).Garnir[rg]+=
          h*inv(Tbasis(H)(d))*Murphybasis(H)(tab1,tab1)*Tbasis(H)(d)
      end
    end
    # Next we worry about the element <w> above (remember t=g*w).
    # This multiplication is usually recursive.
    if !isempty(w)
      tableau(H,t).Garnir[rt]=tableau(H,g).Garnir[rg]*Tbasis(H)(w...)
    end
  end
  # Finally we have to put <s> back into the equation. If we are working
  # in just the Specht module <s> is almost irrelevant; but in general it 
  # affects tnu in strange ways (hence it might be better to cache the
  # full expansion rather than just the right hand side). 
  if s.ind==1 tableau(H,t).Garnir[rt]
  else α(α(tableau(H,t).Garnir[rt])*Tbasis(H)(s.wd))
  end
end

Base.:*(a::HeckeMElt,b::HeckeMElt)=a*Tbasis(a.H)(b)

# This is the anti-isomorphism of the Hecke algebra given by T_i -> T_i;
# on the Murphy basis we have M_{s,t} -> M_{t,s} (also called * by many)
function Garside.α(h::HeckeMElt)
  HeckeMElt(ModuleElt([(mu,t,s)=>c for ((mu,s,t),c) in h.d]),h.H)
end

# Compute the Gram matrix of a Specht module w.r.t. its Murphy basis.
function GramMatrix(H, mu)
  if !SpechtModules[]
    print("# WARNING: in the interests of speed, this function has just \n", 
          "#          disabled T-basis to Murphy basis convertions.\n")
    SpechtModules[]=true
  end
  tab=tableaux(mu)
  M = Murphybasis(H)
  g = fill(Pol(0),length(tab),length(tab))
  for s in 1:length(tab), t in s:length(tab)
    h=M(tab[1],tab[s])*M(tab[t],tab[1])
    g[s,t]=iszero(h) ? 0 : h.d.d[1][2]
    if s!=t g[t,s]=deepcopy(g[s,t]) end
  end
  return g
end

# the Jucys-Murphy elements of the Hecke algebra H. 
function JucysMurphy(H)
  T=Tbasis(H)
  murphy=[0*T() for i in 0:length(H.para)]
  q=H.para[1][1]
  for i in 2:length(H.para)+1
    for j in 1:i-1
      murphy[i]=q^-1*T(j)+q^-1*T(j)*murphy[i]*T(j)
    end
  end
  murphy
end

# A function which returns a function for working with the Murphy basis  
# of M the Specht module S(<mu>).                                       
function SpechtModule(H, mu)
  if sum(mu)!=length(H.para)+1
    error(mu," must be a partition of ",length(H.para)+1)
  end
  SpechtModules[]=true
  tmu=firsttableau(mu)
  t->Murphybasis(H)(tmu, t)
end

# Returns the representation of the Hecke algebra H of type A indexed by
# partition mu, that is the list of matrices of the T_i.
function Spechtmodel(H, mu)
  T=Tbasis(H)
  S=SpechtModule(H, mu)
  tabs=tableaux(mu)
  for t in tabs S(t) end
  n = sum(mu)-1
  zero = sum(sum,H.para)*0
  mats = map(i->map(j->fill(zero,length(tabs)), tabs),1:n)
  for t in tabs
    St = S(t)
    for i in 1:n
      ti=St*T(i)
      mats[i][St.d.d[1][1][3]][map(x->x[3],keys(ti.d))]=collect(values(ti.d))
    end
  end
  mats
end

# Test that T(M(T(w)))=T(w) for all w in Sₙ
function test(n=4,q=Pol())
  W=coxgroup(:A,n-1)
  H=hecke(W, q)
  T=Tbasis(H)
  SpechtModules[]=false
  M=Murphybasis(H)
  xprintln("Testing Murphy functions for ",H)
  for w in sort(words(W), by=a->[length(a), a])
    xprint("checking ",T(w)," ...")
    if T(w)!=T(M(T(w)))
      xprintln("w==",w)
      xprintln("Difference==",T(w)-T(M(T(w))))
      xprintln("T(w)==",T(w))
      xprintln("M(T(w))==",M(T(w)))
      xprintln("T(M(T(w)))==",T(M(T(w))))
      error(" murphy.g FAILED for Sym(", n, ") at T(", join(w), ")  :(\n")
    end
    println("OK")
  end
  xprintln("\n** Murphy.jl passed the test for ",H,"!!")
end

# check that M(T(M(s,t)))=M(s,t) for all pairs (s,t) of the same shape.
function test2(n=4,q=Pol())
  W=coxgroup(:A,n-1)
  H=hecke(W,Pol())
  xprintln("Testing the Murphy basis functions for ",H)
  T=Tbasis(H)
  SpechtModules[]=false
  M=Murphybasis(H)
  for mu in reverse(partitions(n))
    std = tableaux(mu)
    for s in std, t in std
      m=M(s,t)
      xprint("checking ",m," ...")
      if m != M(T(m))
        xprintln(m," --> ",M(T(m)))
        error(" murphy.g FAILED for Sym(",n,")") 
      else println("OK!")
      end
    end
  end
  xprintln("\n** Murphy.jl passed the test for ",H,"!!")
end

end
