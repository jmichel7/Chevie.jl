# Glue code for using GenericDecMats.jl

@GapObj struct ΦDecMat
   W
   d::Integer
   scalar::Matrix{Mvp{Int,Int}}
   order::Vector{Int}
   colnames::Vector{String}
end

function generic_decomposition_matrix(W::FiniteCoxeterGroup,d::Integer)
  iw=isomorphism_type(W)
  mat=GenericDecMats.generic_decomposition_matrix(iw*"d$d")
  if isnothing(mat)
    l=filter(startswith(iw),
                  GenericDecMats.generic_decomposition_matrices_names())
    l=map(x->x[length(iw)+2:end],l)
    xprintln("!!! Φ-decomposition matrices available for ",W,": ",
          join(map(x->fromTeX(rio(),"\\Phi_{$x}"),l)," "))
    return nothing
  end
  uc=UnipotentCharacters(W)
  nn=charnames(uc,TeX=true)
  order=map(n->findfirst(==(n),nn),mat.ordinary)
  if any(isnothing,order) error("not found") end
  blocks=map(unique(mat.blocklabels))do i
   sort(order[filter(j->mat.blocklabels[j]==i,eachindex(mat.blocklabels))])
  end
  ΦDecMat(W,mat.d,mat.decmat[sortperm(order),:],order,mat.hc_series,
    Dict{Symbol,Any}(:blockparams=>mat.blocks,:blocks=>blocks,
                     :condition=>mat.condition,:origin=>mat.origin))
end

function Base.show(io::IO, ::MIME"text/html", m::ΦDecMat)
  show(IOContext(io,:TeX=>true),"text/plain",m)
end

function Base.show(io::IO,m::ΦDecMat)
  printTeX(io,"\\Phi_{",m.d,"}-decomposition matrix for ",m.W)
end

function Base.show(io::IO,::MIME"text/plain",m::ΦDecMat)
  println(io,m)
  scal=map(e->iszero(e) ? "." : TeX(io,e),m.scalar)
  row_labels=charnames(io,UnipotentCharacters(m.W))
  showtable(io,scal;row_labels,col_labels=m.colnames,rows=m.order)
end

