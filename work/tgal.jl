# test galois actions

CyclotomicNumbers.galois(p::Mvp,a)=Mvp(ModuleElt(map(x->x[1]=>Cyc(x[2])^a,pairs(p))))

# test unipotent degrees permuted-with-signs by Gal(K/Q) by a permutation
# which fixes eigenvalues
function testgal(W)
  uc=UnipotentCharacters(W)
  ud=degrees(uc,Mvp(:x))
  ei=collectby(eigen(uc),eachindex(ud))
  map(gens(galois(NF(vec(cartan(W))))))do x
    res=fill(0,length(uc))
    for l in ei
      res[l]=permute(l,SPerm(galois.(ud[l],Ref(x)),ud[l]))
    end
    SPerm(res)
  end
end

# Field(W), families such that Field(Fourier(f))!=Field(W)
function ff(W)
  K=NF(vec(cartan(W)))
  uc=UnipotentCharacters(W)
  KF=pairs(map(f->NF(vec(fourier(f))),uc.families))
  KF=filter(x->!all(y->y in K,x[2].gens),KF)
  if !isempty(KF) println(K,KF) end
end

principal_series(W,d)=
  Series(W,only(cuspidal_data(W, d, length(relative_degrees(W, d))))...)

# check parameters of zeta^-1 series complex conjugate of those of zeta-series
# checkzeta(W[,list of d])
function checkzeta(W,r=nothing)
  if isnothing(r) r=regular_eigenvalues(W) end
  for e in r
    s=conj(hecke(principal_series(W,e)).para)
    s=sort.(s)
    s1=hecke(principal_series(W,conj(e))).para
    s1=sort.(s1)
    if s!=s1 error(repr(e)) end
  end
end

# check that schur elements are globally stable by gal(k_W/Q)
function check(W)
  g=galois(NF(vec(cartan(W))))
  sch=schur_elements(hecke(W,Mvp(:x)))
  map(x->Perm(sch,galois.(sch,Ref(x))),gens(g))
end
