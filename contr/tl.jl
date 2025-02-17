#  TL.jl from TL.g  (C) Francois Digne and Jean Michel  2014
#
#  functions for working with elements of a
#  Temperley-Lieb algebra of a Coxeter group
#
#  we implement unequal parameter TL algebras for non-simply laced coxeter groups

@GapObj struct TL{C,TW}
  W::TW
  parameter::Vector{C}
  x::Union{C,Nothing}
end

function TL(W,p;x=nothing)
  if p isa Tuple
    para=fill(p[1],ngens(W))
    s=simple_reps(W,1:ngens(W))
    C=sort(unique(s))
    for i in 1:ngens(W) para[i]=p[findfirst(==(s[i]),C)] end
  else para=fill(p,ngens(W)) 
  end
  TL(W,para,x,Dict{Symbol,Any}())
end

Groups.Group(H::TL)=H.W
HeckeAlgebras.coefftype(H::TL{C}) where C=C

abstract type TLElt{P,C} end # P=typeof(keys) [Perms] C typeof(coeffs)

TL(h::TLElt)=h.TL # algebra of element

struct TLTElt{P,C,TH}<:TLElt{P,C}
  d::ModuleElt{P,C}
  TL::TH
end

basename(h::TLTElt)="tl"

function IsFullyCommutative(W,w)
  N=refls(W,inversions(W,reverse(word(W,w))))
  c(p)=isone(comm(p[1],p[2])) ? one(W) : p[1]^p[2]
  isempty(intersect(N,c.(cartesian(N,N))))
end

tlbasis(H::TL)=(x...)->x==() ? one(H) : tlbasis(H,x...)
tlbasis(H::TL,w::Vararg{Integer})=tlbasis(H,H.W(w...))
tlbasis(H::TL,w::Vector{<:Integer})=tlbasis(H,H.W(w...))
tlbasis(H::TL,h::TLTElt)=h
tlbasis(H::TL,h::TLElt)=tlbasis(h)
tlbasis(H::TL,w)=TLTElt(ModuleElt(w=>one(coefftype(H));check=false),H)

Base.one(H::TL)=tlbasis(H,one(H.W))

Base.:^(a::TLElt, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                     Base.power_by_squaring(inv(a),-n)

function Base.:*(x::TLTElt{P,C},y::TLTElt{P,C})where {P,C}
  H=TL(x)
  W=Group(H)
  res=zero(ModuleElt{P,C})
  for (i,(px,cx)) in enumerate(pairs(x.d))
    temp=y.d*cx
    xi=inv(px)
    while xi!=one(W)
      s=firstleftdescent(W,xi)
      es=W(s)
      xi=es*xi
      temp1=Pair{P,C}[]
      for (j,(e,c)) in enumerate(pairs(temp))
        if isleftdescent(W,e,s)
          q=H.parameter[s]
          push!(temp1,e=>(q+1)*c)
        else phie=inversions(W,word(W,e))
          p=findfirst(x->roots(W,x)-roots(W,s) in roots(W,phie),phie)
   #   on a une relation s2=k*s2*s1*s2 si dans N(s2s1s2) il y a
   #   2 racines dont la somme est une racine. Cette relation intervient
   #   dans es*e ssi dans N(es*e) ... ou encore ssi dans N(e) il y a une
   #   racine dont la  difference avec W.roots[s] est une racine.
          if isnothing(p) push!(temp1,es*e=>c)
          else r=findfirst(==(roots(W,phie[p])-roots(W,s)),roots(W))
            w=refls(W,r)*e
            q=H.parameter[s]
            if IsFullyCommutative(W,w)
              if cartan(W,s,r)==-2
                Q=H.parameter[simple_reps(W,r)]
                coeff=(q+Q)*c
              elseif cartan(W,s,r)==-1
                coeff=q*c
              elseif cartan(W,s,r)==-3
                Q=H.parameter[simple_reps(W,r)]
                coeff=(q+Q+GetRoot(q*Q))*c
              else coeff=H.x*c
              end
              push!(temp1,w=>coeff)
            else
              w1=with_inversions(W,phie[1:findfirst(==(r),phie)-1])
              w2=inv(w1)*refls(W,r)*e
              tlw1=TLTElt(H,ModuleElt(w1=>one(coefftype(H))))
              tlw2=TLTElt(H,ModuleElt(w2=>one(coefftype(H))))
              println("recursive case: ",s,"*",BraidMonoid(W)(e),
                      " tl(w1)=",tlw1," tl(w2)=",tlw2)
              w=tlw1*tlw2
              if cartan(W,s,r)==-2
                Q=H.parameter[simple_reps(W,r)];
                w*=(q+Q)*c
              elseif cartan(W,s,r)==-1
                w*=q*c
              elseif cartan(W,s,r)==-3
                Q=H.parameter[simple_reps(W,r)];
                w*=(q+Q+root(q*Q))*c
              else w*=H.x*c
              end
              append!(temp1,pairs(w))
            end
          end
        end
      end
      temp=ModuleElt(temp1)
    end
    res+=temp
  end
  TLTElt(res,H)
end

#F AlphaInvolution(h)    
## The involution on TL Elements defined by T_w->T_{w^{-1}}
## (and same in other bases)
Garside.α(h::TLTElt)=TLTElt(ModuleElt(inv(p)=>c for (p,c) in h.d),TL(h))

function Base.show(io::IO, h::TLElt)
  function showbasis(io::IO,e)
    w=word(TL(h).W,e)
    res=basename(h)
    if hasdecor(io) res*=isempty(w) ? "." : "_"*joindigits(w,"{}";always=true)
    else            res*="("*join(w,",")*")"
    end
    fromTeX(io,res)
  end
  show(IOContext(io,:showbasis=>showbasis),h.d)
end

## The function below is made member of CoxeterTLAlgebraOps just to hide it
## it computes T(w^-1)^-1
#CoxeterTLAlgebraOps.getTinv:=function(H,w)local p,i,q,W;
#  Error("not implemented");
#  W:=Group(H);
#  if not IsBound(H.elts) then H.elts:=[W.identity]; fi;
#  if not IsBound(H.Tinv) then H.Tinv:=[Basis(H,"t")()];fi;
#  p:=Position(H.elts,w);
#  if p=false then Add(H.elts,w);p:=Length(H.elts);fi;
#  if not IsBound(H.Tinv[p]) then
#    i:=FirstLeftDescending(W,w);
#    q:=H.parameter[i];
#    H.Tinv[p]:=Basis(H,"t")([W.identity,W.reflections[i]],
#	               [(q[1]+q[2])/(q[1]*q[2]),-1/(q[1]*q[2])])*
#                          H.operations.getTinv(H,W.reflections[i]*w);
#  fi;
#  return H.Tinv[p];
#end;
#
#Matrixtltot:=function(H)local W,l,t,tl,M;
#  if not IsBound(H.tltot) then
#    W:=H.group;
#    l:=Filtered(Elements(W),x->IsFullyCommutative(W,x));
#    sort!(l,by=x->length(W,x))
#    t:=Basis(H,"t"); tl:=Basis(H,"tl");
#    H.elts:=l;
#    M:=List(List(l,x->tl(t(x))),x->List(l,y->Coefficient(x,y)));
#    H.tltot:=M^-1;
#  fi;
#  return H.tltot;
#end;
#
##############################################################################
###
### An exemple of CreateTLBasis: create the 'T' basis for
### CoxeterTLAlgebraOps (i.e. all TL algebras for Coxeter groups...)
###
#CreateTLBasis("t",rec(
#  tl:=function(x) local W,tl;
#    tl:=Basis(x.TL,"tl");
#    W:=x.TL.group;
#    return Sum([1..Length(x.coeff)],function(i)
#      if x.elm[i]=() then return x.coeff[i]*tl();
#      else return x.coeff[i]*Product(CoxeterWord(W,x.elm[i]),j->tl(j)-tl());
#      fi;
#    end);
#  end,
#  t:=function(x) local tl,H,M;H:=x.TL;tl:=Basis(H,"tl");M:=Matrixtltot(H);
#    return Sum([1..Length(x.coeff)],i->x.coeff[i]*TLEltOps.Normalize(
#     TLEltOps.MakeRec(H,"t",ShallowCopy(H.elts),
#        ShallowCopy(M[Position(H.elts,x.elm[i])]))));
#  end,
#  inverse:=function(h)
#   if Length(h.elm)<>1 then Error("inverse implemented only for single t_w");fi;
#   return h.coeff[1]^-1*TL(h).operations.getTinv(TL(h),h.elm[1]^-1);
#  end),CoxeterTLAlgebraOps);
#
#CreateTLBasis("tl",rec(tl:=x->x)     # method to convert to tl
#,CoxeterTLAlgebraOps);
