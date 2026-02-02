"""
This module `Murphy.jl` has been ported in december 2020 from

   murphy.g Copyright (C) July 1998  Andrew Mathas 
   mathas@maths.usyd.edu.au University of Sydney

It allows computations with the Murphy basis of an Hecke algebra of type A.

Multiplication  of Murphy basis elements is  done using the Garnir tableaux
as  described in [mur95](@cite). This also lets us convert from the T-basis
to the Murphy basis since `T_w = M([[1],…,[n]], [[1],…,[n]]) * T_w`; we use
"M" for the Murphy basis.

As with the T-basis, Murphy basis elements are implemented by `ModuleElts`.
Here the keys are standard tableaux pairs. These are represented by a tuple
`(mu,s,t)`  where `mu`,  `s` and  `t` are  integers encoding  the partition
`H.Murphy.partitions[mu]`  and the  named tuples `H.Murphy.tableaux[mu][s]`
(resp   `t`)  describing  tableaux  as  explained  in  the  helpstring  for
[`initMurphy`](@ref).

Throughout memory considerations are thrown to the wind as we cache many of
the  more horrible expansions as we go along  in order to save time when we
next need them.
```julia-repl
julia> W=coxgroup(:A,2)
A₂

julia> H=hecke(W,Pol(:q))
hecke(A₂,q)

julia> l=Tbasis(H).(elements(W))
6-element Vector{HeckeTElt{HeckeAlgebra{Pol{Int64}, Perm{Int16}, FiniteCoxeterGroup{Perm{Int16},Int64}}, Pol{Int64}, Perm{Int16}}}:
 T.
 T₁
 T₂
 T₁₂
 T₂₁
 T₁₂₁

julia> Murphy.SpechtModules(H,false);Murphybasis(H).(l)
6-element Vector{Chevie.Murphy.HeckeMElt{Pol{Int64}, HeckeAlgebra{Pol{Int64}, Perm{Int16}, FiniteCoxeterGroup{Perm{Int16},Int64}}}}:
 M(1/2/3,1/2/3)
 -M(1/2/3,1/2/3)+M(12/3,12/3)
 -M(1/2/3,1/2/3)+q⁻¹M(12/3,12/3)+q⁻¹M(12/3,13/2)+q⁻¹M(13/2,12/3)+q⁻¹M(13/2,13/2)-q⁻¹M(123,123)
 M(1/2/3,1/2/3)-q⁻¹M(12/3,12/3)+(1-q⁻¹)M(12/3,13/2)-q⁻¹M(13/2,12/3)-q⁻¹M(13/2,13/2)+q⁻¹M(123,123)
 M(1/2/3,1/2/3)-q⁻¹M(12/3,12/3)-q⁻¹M(12/3,13/2)+(1-q⁻¹)M(13/2,12/3)-q⁻¹M(13/2,13/2)+q⁻¹M(123,123)
 -M(1/2/3,1/2/3)+(-1+q⁻¹)M(12/3,12/3)+(-1+q⁻¹)M(12/3,13/2)+(-1+q⁻¹)M(13/2,12/3)+q⁻¹M(13/2,13/2)+(1-q⁻¹)M(123,123)
```
"""
module Murphy

using ..Chevie
export Murphybasis, Spechtmodel

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
  w
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

istableau(t)=sort(vcat(t...))==1:sum(length,t) && all(issorted,t) &&
             all(issorted,conjugate_tableau(t))

#------------------The Murphy basis code proper---------------------
"""
Set  `SpechtModules(H,true)` to work just  inside the Specht modules.
This  makes computations with  the Murphy basis  inside Specht modules much
faster  but also means that `T` to  `Murphy` basis conversions do not work,
even if `SpechtModules(H,false)` is set later.
"""
function SpechtModules(H,t::Bool)
  initMurphy(H)
  if H.Murphy.SpechtModules[]==t return end
  empty!(H.Murphy.TtoMurphy)
  for i in eachindex(H.Murphy.partitions)
    for t in H.Murphy.InfoTableaux[i]
      empty!(t.Garnir)
      empty!(t.toT)
    end
  end
  H.Murphy.SpechtModules[]=t
end
SpechtModules(H)=H.Murphy.SpechtModules[]

struct HeckeMElt{C,TH}<:HeckeElt{TH,C,Tuple{Int,Int,Int}}
  d::ModuleElt{Tuple{Int,Int,Int},C} # has better merge performance than Dict
  H::TH
end

HeckeAlgebras.clone(h::HeckeMElt,d)=HeckeMElt(d,h.H)
Base.zero(::Type{HeckeMElt},H::HeckeAlgebra)=HeckeMElt(zero(ModuleElt{Tuple{Int,Int,Int},coefftype(H)}),H)
Base.zero(h::HeckeMElt)=zero(HeckeMElt,h.H)
HeckeAlgebras.basisname(h::HeckeMElt)="M"

Murphycache=Dict{Tuple{Int,Any},Any}()

"""
`initMurphy(H)` for `hecke(coxgroup(:A,n-1),q)`

Called the first time that the Murphy basis is used. Creates H.Murphy with 
various components.

  - `H.Murphy.partitions` contains the partitions of n
  - `H.Murphy.Tableaux[mu]` contains the tableaux for `H.Murphy.partitions[mu]`
  - `H.Murphy.InfoTableaux[mu]` contains InfoTableaux for H.Murphy.Tableaux[mu]

Creating  a list  of all  of the  standard tableaux  is a  big overhead for
Specht  modules of  large dimension  so the  above arrays  work as  a cache
filled  as  needed.  The  bookkeeping  to  maintain these caches is done in
`InfoTableau`.
"""
function initMurphy(H)
  get!(H, :Murphy)do
    t=refltype(H.W)
    if length(t)!=1 || t[1].series!=:A
      error("the Murphy basis is implemented only for irreducible type A")
    end
    if ngens(H.W)>0 && H.para[1][2]!=-1
      error("H must have parameters q and -1")
    end
    get!(Murphycache,(ngens(H.W),H.para[1][1]))do
      InfoChevie("#I Initialized Murphy basis\n")
      H.Murphy=(partitions=Vector{Int}[],
                InfoTableaux=Vector{InfoTableau}[],
       Tableaux=Vector{Vector{Vector{Int}}}[],
       TtoMurphy=Dict{eltype(H.W),Any}(), # T-basis to Murphy-basis cache.
       SpechtModules=Ref(true),
       SpechtModels=Dict{Any,Any}()
      )
      code_partition(H,fill(1,ngens(H.W)+1)) # first partition of n
      H.Murphy
    end
  end
end

struct InfoTableau{P}
  ind::Int #index of tableau in H.Murphy.tableaux[mu]
  mu::Int  #index of tableau's shape in H.Murphy.partitions
  wd::P    # associated word in S_n (as a Perm)
           # that is Perm such that t=t(mu)*w where t(mu)=firsttableau(mu)
  Garnir::Dict{Int64,Any} # will hold the Garnir expansions
  toT::Dict{Int64,Any}  # MurphyToT 
end

"""
`code_partition(H::HeckeAlgebra, μ)`

Given a partition `μ⊢n`, return the index `mu` of `μ` in `H.partitions`, or
add  `μ`  to  this  list  if  it  is  not  already  there  *and* initialize
`H.InfoTableaux[mu]` and `H.Tableaux[mu]`.
"""
function code_partition(H::HeckeAlgebra, mu)
  t=findfirst(==(mu),H.Murphy.partitions)
  if t!==nothing return t end
  push!(H.Murphy.partitions,copy(mu))
  t=length(H.Murphy.partitions)
  push!(H.Murphy.InfoTableaux,
        [InfoTableau(1,t,one(H.W),Dict{Int,Any}(),Dict{Int,Any}())])
  push!(H.Murphy.Tableaux,[firsttableau(mu)])
  t
end

# given  a tableau  <tab> ,  which is  assumed to  be standard, return the
# corresponding  element of <H>.tableau  (a "CoxeterGroup tableau"). Here,
# <mu> is the index of the shape of <tab> in H.partitions[].
function InfoTableau(H::HeckeAlgebra, mu::Integer, tab::Vector)
  tabs=H.Murphy.Tableaux[mu]
  ind=findfirst(==(tab),tabs)
  if ind===nothing
    ind=length(tabs)+1
    push!(tabs,copy.(tab))
    push!(H.Murphy.InfoTableaux[mu],InfoTableau(length(tabs),mu,
          H.W(Permtableaux(tabs[1],tab)...),Dict{Int,Any}(),Dict{Int,Any}()))
  end
  H.Murphy.InfoTableaux[mu][ind]
end

InfoTableau(H::HeckeAlgebra,mu::Integer,tab::Integer)=H.Murphy.InfoTableaux[mu][tab]

firstTableau(H::HeckeAlgebra,t::InfoTableau)=H.Murphy.Tableaux[t.mu][1]
Tableau(H::HeckeAlgebra,t::InfoTableau)=H.Murphy.Tableaux[t.mu][t.ind]

# given a tableau <t> return the basis element x_mu*T_<t> in the T-basis
function Xt(H::HeckeAlgebra,t::InfoTableau)
  mm=get!(t.toT,1)do
    if t.ind==1 
      J=vcat(map(x->x[1:end-1],firstTableau(H,t))...)
      WJ=elements(reflection_subgroup(H.W,J))
      ModuleElt(map(x->x=>one(coefftype(H)),WJ))
    else 
      ModuleElt([k*t.wd=>c for (k,c) in Xt(H,InfoTableau(H,t.mu,1)).d])
    end
  end
  HeckeTElt(mm,H)
end

# given a pair (s,t) of tableaux, return the basis element 
# <t>.(<s>.ind)=T_<s>^* x_mu T_<t> in the T-basis.
function MurphyToT(H::HeckeAlgebra, s::InfoTableau, t::InfoTableau)
  mm=get!(t.toT,s.ind)do
    if s.ind==1 Xt(H,t).d
    elseif t.ind==1 α(Xt(H, s)).d
    else (Tbasis(H)(inv(s.wd))*Xt(H, t)).d
    end
  end
  HeckeTElt(mm,H)
end

# convert Murphy basis to T-basis
function HeckeAlgebras.Tbasis(M::HeckeMElt)
 sum(c*MurphyToT(M.H,InfoTableau(M.H,i,t1),InfoTableau(M.H,i,t2))
       for ((i,t1,t2),c) in M.d)
end

# This function recursively expands T_w into a linear combination of
# Murphy basis elements (using Garnir expansions). Note that we know
# how to write 1 in terms of the Murphy basis. 
# As we go along we cache these expansions in the dict H.Murphy.TtoMurphy
function TtoMurphy(H::HeckeAlgebra,w)
  mm=get!(H.Murphy.TtoMurphy,w)do
    if SpechtModules(H)
      println("\n# WARNING: because SpechtModules(H)==true the answer\n",
            "# TtoMurphy returns will almost certainly be incorrect.")
    end
    if isone(w) ModuleElt((1,1,1)=>one(coefftype(H)))
    else
      W=H.W
      r=firstleftdescent(W,w^-1)
      res=TtoMurphy(H,w*W(r))*Tbasis(H)(r)
      if w!=w^-1 H.Murphy.TtoMurphy[w^-1]=α(res).d end
      res.d
    end
  end
  HeckeMElt(mm,H)
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
    is=InfoTableau(H, imu, s)
    it=InfoTableau(H, imu, t)
    HeckeMElt(ModuleElt((imu,is.ind,it.ind)=>one(coefftype(H))),H)
  end
end

function StringTableau(io::IO,t)
  TeX=get(io,:TeX,false)
  repl=get(io,:limit,false)
  if TeX return "\\tab("*join(join.(tab),",")*")"
  elseif repl return join(joindigits.(t),"/")
  else "["*join(map(p->"["*join(p,",")*"]",t),",")*"]"
  end
end

function Base.show(io::IO, h::HeckeMElt)
  function showbasis(io::IO,(mu,s,t))
    TeX=get(io,:TeX,false)
    t=h.H.Murphy.Tableaux[mu][t]
    s=h.H.Murphy.Tableaux[mu][s]
    if SpechtModules(h.H)
      TeX ? StringTableau(io,t)*"\n" : string("S(",StringTableau(io,t),")")
    else 
      string(HeckeAlgebras.basisname(h),"(",StringTableau(io,s),",",
                              StringTableau(io,t),")",TeX ? "\n" : "")
    end
  end
  show(IOContext(io,:showbasis=>showbasis),h.d)
end

# returns the position of the integer <i> in the tableau <tab>.
function PositionIn(tab, i)
  for row in eachindex(tab)
    col=findfirst(==(i),tab[row])
    if col!==nothing return (;row,col) end
  end
  error(i," is not contained in the tableau ",StringTableau(rio(),tab))
end

function Base.:*(m::HeckeMElt, h::HeckeTElt)
  H=h.H
  if H!==m.H error("not elements of the same algebra") end
  W=H.W
  q=H.para[1][1]
  # mh=return value, initially zero with respect to the Murphy basis
  mh=zero(HeckeMElt,H)
  for (w,hcoeff) in h.d
    mw=m
    for r in word(W,w) # multiply one simple reflection at a time
      mr=zero(m.d).d
      for ((mu,is,it),coeff) in mw.d
        t=InfoTableau(H,mu,it)
        # we are interested in the two nodes, <nodeR> and <nodeS> 
        # which are swapped by the transposition r=(r,r+1). Thus,
        # these are the nodeRs such that t<nodeR>=r and t<nodeS>=r+1.
        nodeR=PositionIn(Tableau(H,t),r)
        nodeS=PositionIn(Tableau(H,t),r+1)
        if nodeR.row == nodeS.row
          push!(mr,(mu,is,it)=>q*coeff)
        elseif nodeR.col!=nodeS.col# then t*r is still standard
          tabr=copy.(H.Murphy.Tableaux[mu][it])
          tabr[nodeR.row][nodeR.col]=r+1
          tabr[nodeS.row][nodeS.col]=r
          tr=InfoTableau(H, mu, tabr)
          if nodeR.row < nodeS.row # up in the Bruhat order
            push!(mr,(mu,is,tr.ind)=>coeff)
          else
            push!(mr,(mu,is,tr.ind)=>q*coeff)
            push!(mr,(mu,is,it)=>(q-1)*coeff)
          end
        else  # The hard part: here in the tableau t, r+1 occupies the nodeR 
              # below r; so nodeS.row=nodeR.row+1 and nodeS.col=nodeR.col and
              # interchanging them gives a non-standard tableau.
              s=InfoTableau(H,mu,is)
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
function GarnirExpansion(H::HeckeAlgebra,node,s::InfoTableau,t::InfoTableau)
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
    gtab=copy.(firstTableau(H,t))
    a=gtab[node.row][node.col] # first number being moved; above a=5 and b=8
    b=gtab[node.row+1][node.col] # last number being moved
    gtab[node.row][node.col+1:end]=a+node.col+1:b
    gtab[node.row+1][1:node.col-1]=a:a+node.col-2
    gtab[node.row][node.col]  =a+node.col-1
    gtab[node.row+1][node.col]=a+node.col
    g=InfoTableau(H, t.mu, gtab)
    # w is the permutation such that t=g*w => T_t = T_g*T_w
    w=Permtableaux(gtab,Tableau(H,t))
    rg=Tableau(H,g)[node.row][node.col]
    if !haskey(g.Garnir,rg)
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
          push!(mres,(g.mu,1,InfoTableau(H,g.mu,gtab).ind)=>-1)
        end
      end
      g.Garnir[rg]=ModuleElt(mres)
      # Next, if SpechtModules(H)==false (in which case we 
      # just work in the Specht module), we look after the right hand term
      # in (*) above. In general it won't correspond to a partition but we
      # can find a partition <nu> and a <d> in <W> such that 
      # T_d<RHS>=<x_nu>T_<d>
      if !SpechtModules(H)
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
        tnu=InfoTableau(H,code_partition(H,length.(tab)),1)
        # and finally we have <d>. The point is that tab = T_d^-1*tnu*T_d.
        d=W(Permtableaux(Tableau(H,tnu), tab)...)
        # <tab> is now under control, but we still need to compute <h>
        # from (*). The point here is that we are essentially writing
        #             H.I = ⋃ H.d = ⋃ d'.I
        # where H and I are two sugroups and d and d' run over coset
        # representatives of H and I.
        gtab=firstTableau(H,g)
        J=gtab[node.row][1:end-1]
        append!(J,gtab[node.row+1][1:end-1])
        K=setdiff(J,[a-1,b])
        h=inv.(vcat(reduced(reflection_subgroup(W,K),
                            reflection_subgroup(W,J))...))
        h=HeckeTElt(ModuleElt([w=>one(coefftype(H)) for w in h];check=false),H)
        # the multiplication below is quite costly as it is recursive; 
        # but it is only done once as we store the result in g.Garnir.
        tab1=firstTableau(H,tnu)
        g.Garnir[rg]+=(h*inv(Tbasis(H)(d))*Murphybasis(H)(tab1,tab1)*Tbasis(H)(d)).d
      end
    end
    # Next we worry about the element <w> above (remember t=g*w).
    # This multiplication is usually recursive.
    if !isempty(w)
      t.Garnir[rt]=(HeckeMElt(g.Garnir[rg],H)*Tbasis(H)(w...)).d
    end
  end
  # Finally we have to put <s> back into the equation. If we are working
  # in just the Specht module <s> is almost irrelevant; but in general it 
  # affects tnu in strange ways (hence it might be better to cache the
  # full expansion rather than just the right hand side). 
  if s.ind==1 HeckeMElt(t.Garnir[rt],H)
  else α(α(HeckeMElt(t.Garnir[rt],H))*Tbasis(H)(s.wd))
  end
end

Base.:*(a::HeckeMElt,b::HeckeMElt)=a*Tbasis(a.H)(b)

# This is the anti-isomorphism of the Hecke algebra given by T_i -> T_i;
# on the Murphy basis we have α(M_{s,t})=M_{t,s} (α also called * by many)
function Garside.α(h::HeckeMElt)
  HeckeMElt(ModuleElt((mu,t,s)=>c for ((mu,s,t),c) in h.d),h.H)
end

# Compute the Gram matrix of a Specht module w.r.t. its Murphy basis.
function GramMatrix(H, mu)
  if !SpechtModules(H)
    print("# WARNING: in the interests of speed, this function has just \n", 
          "#          disabled T-basis to Murphy basis convertions.\n")
    SpechtModules(H,true)
  end
  tab=tableaux(mu)
  M=Murphybasis(H)
  g=fill(Pol(0),length(tab),length(tab))
  for s in 1:length(tab), t in s:length(tab)
    h=M(tab[1],tab[s])*M(tab[t],tab[1])
    g[s,t]=iszero(h) ? 0 : h.d.d[1][2]
    if s!=t g[t,s]=copy.(g[s,t]) end
  end
  g
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
  SpechtModules(H,true)
  tmu=firsttableau(mu)
  t->Murphybasis(H)(tmu, t)
end

# Returns the representation of the Hecke algebra H of type A indexed by
# partition mu, that is the list of matrices of the T_i.
function Spechtmodel(H, mu)
  initMurphy(H)
  get!(H.Murphy.SpechtModels,mu)do
  tabs=tableaux(mu)
  n=sum(mu)-1
  mats=map(i->fill(sum(sum,H.para)*0,length(tabs),length(tabs)),1:n)
  for t in tabs
    St=SpechtModule(H,mu)(t)
    for i in 1:n
      ti=St*Tbasis(H)(i)
      mats[i][last(first(keys(St.d))),last.(keys(ti.d))].=values(ti.d)
    end
  end
  mats
  end
end

function test(n=4,q=Pol();rep=false)
  W=coxgroup(:A,n-1)
  H=hecke(W, q)
  T=Tbasis(H)
  M=Murphybasis(H)
  if !rep
  SpechtModules(H,false)
  xprintln("Testing the Murphy basis functions for ",H)
# Test that T(M(T(w)))=T(w) for all w in Sₙ
  for w in sort(words(W), by=a->[length(a), a])
    xprint("checking ",T(w)," ...")
    if !iszero(T(w)-T(M(T(w))))
      xprintln("w==",w)
      xprintln("difference==",T(w)-T(M(T(w))))
      xprintln("T(w)==",T(w))
      xprintln("M(T(w))==",M(T(w)))
      xprintln("T(M(T(w)))==",T(M(T(w))))
      error(" murphy.jl FAILED for Sym(", n, ") at T(", join(w), ")  :(\n")
    end
    println("OK")
  end
# check that M(T(M(s,t)))=M(s,t) for all pairs (s,t) of the same shape.
  for mu in reverse(partitions(n))
    std=tableaux(mu)
    for s in std, t in std
      m=M(s,t)
      xprint("checking ",m," ...")
      if !iszero(m-M(T(m)))
        xprintln(m," --> ",M(T(m)))
        error(" murphy.jl FAILED for Sym(",n,")") 
      else println("OK!")
      end
    end
  end
  else
  for mu in partitions(n)
    r=Spechtmodel(H,mu)
    if isrepresentation(H,r) println("representation ",mu," OK!") end
  end
  end
  xprintln("\n** Murphy.jl passed the tests for ",H,"!!")
  H
end

end
