# Glue code for using GenericDecMats.jl

@GapObj struct ΦDecMat{T}
   W
   d::Integer
   scalar::Matrix{T}
   order::Vector{Int}
   colnames::Vector{String}
end

function generic_decomposition_matrix(W,d::Integer)
  mm=generic_decomposition_matrix.(refltype(W),d)
  names=join.(cartesian(getfield.(mm,:ordinary)...),"\\otimes ")
  hcseries=join.(cartesian(getfield.(mm,:hc_series)...),",")
  mat=length(mm)==1 ? mm[1].decmat : kron(getfield.(mm,:decmat)...)
  uc=UnipotentCharacters(W)
  nn=charnames(uc,TeX=true)
  order=map(n->findfirst(==(n),nn),names)
  if any(isnothing,order) 
    p=findfirst(isnothing,order)
    error(names[p],
       " not found in charnames(UnipotentCharacters(",W,"),TeX=true)")
  end
# blocks=map(unique(mat.blocklabels))do i
#  sort(order[filter(j->mat.blocklabels[j]==i,eachindex(mat.blocklabels))])
# end
  if length(order)<length(nn)
    res=Array{Union{Missing,eltype(mat)}}(missing,length(nn),length(hcseries))
    res[order,:]=mat
  else res=mat[sortperm(order),:]
  end
  ΦDecMat(W,d,res,order,hcseries,
    Dict{Symbol,Any}(#:blockparams=>mat.blocks,:blocks=>blocks,
                     :condition=>join(getfield.(mm,:condition),"  "),
                     :origin=>join(getfield.(mm,:origin),"  ")))
end

function generic_decomposition_matrix(t::TypeIrred,d::Integer)
  field=string(Chevie.field(t)...)
  mat=GenericDecMats.generic_decomposition_matrix(field*"d$d")
  if isnothing(mat)
    l=filter(startswith(field),
                  GenericDecMats.generic_decomposition_matrices_names())
    if isempty(l)
       xprintln("!!! no Φ-decomposition matrices available for ",t)
    else
      l=map(x->x[length(field)+2:end],l)
      xprintln("!!! Φ-decomposition matrices available for ",t,": ",
            join(map(x->fromTeX(rio(),"\\Phi_{$x}"),l)," "))
    end
  end
  mat
end

function Base.show(io::IO, ::MIME"text/html", m::ΦDecMat)
  show(IOContext(io,:TeX=>true),"text/plain",m)
end

function Base.show(io::IO,m::ΦDecMat)
  printTeX(io,"\\Phi_{",m.d,"}-decomposition matrix for ",m.W)
end

function Base.show(io::IO,::MIME"text/plain",m::ΦDecMat)
  println(io,m)
  scal=map(e->(!ismissing(e) && iszero(e)) ? "." : TeX(io,e),m.scalar)
  row_labels=charnames(io,UnipotentCharacters(m.W))
  showtable(io,scal;row_labels,col_labels=m.colnames,rows=m.order)
end

struct InducedDecompositionMatrix{T}
  scalar::Matrix{T}
  colnames::Vector{String}
  W
  R
  d::Integer
end

function InducedDecompositionMatrix(R,W,d::Integer)
  m=generic_decomposition_matrix(R,d)
  t=HCInductionTable(R,W)
  InducedDecompositionMatrix(t.scalar*m.scalar,m.colnames,W,R,d)
end

function Base.show(io::IO,::MIME"text/plain",m::InducedDecompositionMatrix)
 printTeX(io,"Induced \\Phi_{",m.d,"}-decomposition matrix from ",
          m.R," to ",m.W,"\n")
  scal=map(e->(!ismissing(e) && iszero(e)) ? "." : TeX(io,e),m.scalar)
  row_labels=charnames(io,UnipotentCharacters(m.W))
  order=sort(axes(m.scalar,1),
              by=i->maximum(findall(!iszero,m.scalar[i,:]);init=0))
  showtable(io,scal;row_labels,col_labels=m.colnames,rows=order)
end

