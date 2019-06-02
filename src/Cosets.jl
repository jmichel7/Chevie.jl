module Cosets

using Gapjm
export twistings, spets 
abstract type Coset{TW} end

abstract type Spets{TW}<:Coset{TW} end

abstract type CoxeterCoset{TW}<:Spets{TW} end

class_reps(W::FiniteCoxeterGroup)=map(x->W(x...),classinfo(W)[:classtext]) 

function twisting_elements(W::FiniteCoxeterGroup,J::AbstractVector{<:Integer})
  if isempty(J) class_reps(W)
  else class_reps(centralizer(W,Set(J);action=(J,w)->Set(J.^w)))
  end
end

function twistings(W::FiniteCoxeterGroup,J::AbstractVector{<:Integer})
  R=reflection_subgroup(W,J)
  map(w->spets(R,w),twisting_elements(W,J))
end

struct FCC{T,T1,TW<:FiniteCoxeterGroup{Perm{T},T1}}<:CoxeterCoset{TW}
  phi::Perm{T}
  F::Matrix
  W::TW
  prop::Dict{Symbol,Any}
end

spets(W::FiniteCoxeterGroup,w::Perm)=spets(W,matX(W,w))

function spets(W::FiniteCoxeterGroup{Perm{T},T1},F::Matrix) where{T,T1}
  perm=Perm{T}(F,roots(parent(W.G)),action=(r,m)->permutedims(m)*r)
  phi=reduced(W,perm)
  FCC(phi,F,W,Dict{Symbol,Any}())
end

function det(m::Matrix)
  if size(m,1)==1 return m[1,1] end
  sum(i->(-1)^(i-1)*m[i,1]*det(m[vcat(1:i-1,i+1:size(m,1)),2:end]),1:size(m,1))
end

function torusfactors(WF::CoxeterCoset)
  M=PermRoot.baseX(WF.W.G)
  M=Int.(M*WF.F*inv(Rational.(M)))
  r=length(gens(WF.W))
  M=M[r+1:end,r+1:end]
  M=M-Ref(Pol([1],1)).*one(M)
  CycPol(det(M))
end

function PermRoot.refltype(WF::CoxeterCoset)::Vector{TypeIrred}
  gets(WF,:refltype)do WF
    t=refltype(WF.W)
    c=map(x->PermRoot.indices(x),t)
    phires=Perm(WF.phi,inclusion(WF.W))
    map(cycles(Perm(sort.(c),map(i->sort(i.^phires),c))))do c
      o=deepcopy(t[c])
      J=PermRoot.indices(o[1])
      twist=Perm{Int16}(phires^length(c),J)
      if o[1][:series]==:D && length(J)==4
        if order(twist)==2
          rf=findall(i->i!=i^twist,J)
          o[1].prop[:indices]=J[vcat(rf,[3],setdiff([1,2,4],rf))]
        elseif order(twist)==3 && J[1]^twist!=J[2]
          o[1].prop[:indices]=J[[1,4,3,2]]
        end
      end
      for i in 2:length(c) 
        o[i].prop[:indices]=PermRoot.indices(o[i-1]).^phires
      end
      TypeIrred(Dict(:orbit=>o,:twist=>twist))
    end
  end
end

function Base.show(io::IO, W::CoxeterCoset)
   show(io,refltype(W))
   show(io,torusfactors(W))
end

end
