# Glue code for using GenericDecMats.jl
using GenericDecMats 
import GenericDecMats: generic_decomposition_matrix
export GenericDecMats, generic_decomposition_matrix, InducedDecompositionMatrix

@GapObj struct ΦDecMat{T}
   W
   d::Integer
   scalar::Matrix{T}
   order::Vector{Int}
   colnames::Vector{String}
end

"""
`generic_decomposition_matrix(W,d)`

This  function  obtains  the  `Φ_d`-decomposition  matrix for the reductive
group  specified  by  the  Coxeter  group  or  coset  `W`  from the package
[GenericDecMats](https://github.com/oscar-system/GenericDecMats.jl).

```julia-repl
julia> W=rootdatum(:psu,5)
psu₅
```
```julia-rep1
julia> generic_decomposition_matrix(W,13)
!!! Φ-decomposition matrices available for ²A₄: Φ₁₀ Φ₂ Φ₄ Φ₆
```
```julia-repl
julia> generic_decomposition_matrix(W,10)
Φ₁₀-decomposition matrix for psu₅
      │ps 21 ps ps ps 2111 11111
──────┼──────────────────────────
2.    │ 1  .  .  .  .    .     .
²A₂:2 │ .  1  .  .  .    .     .
11.   │ .  .  1  .  .    .     .
1.1   │ 1  .  .  1  .    .     .
.2    │ .  .  .  .  1    .     .
²A₂:11│ .  1  .  .  .    1     .
.11   │ .  .  .  1  .    .     1
```

The matrix itself is stored in the field `.scalar` of the returned `struct`.
"""
function generic_decomposition_matrix(W,d::Integer)
  mm=generic_decomposition_matrix.(refltype(W),d)
  if nothing in mm return nothing end
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
  if haskey(t,:orbit)
    if length(t.orbit)>2 
      error("orbits of F should not be of size ",length(t.orbit))
    elseif length(t.orbit)==2 
      d=div(d,gcd(d,2))
    end
  end
  if d==1
    if !haskey(t,:orbit) error("d=1 allowed only for twisted groups") end
    n=charnames(UnipotentCharacters(t);TeX=true)
    return (ordinary=n,hc_series=n,
          decmat=Int.(eachindex(n) .== eachindex(n)'),condition="",origin="")
  end
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
  if !hasdecor(io)
    print(io,"generic_decomposition_matrix(",m.W,",",m.d,")")
  else
    printTeX(io,"\\Phi_{",m.d,"}-decomposition matrix for ",m.W)
  end
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

"""
`InducedDecompositionMatrix(R,W,d)`

returns the induced from the Levi `L` to the reductive group `W` of the
generic `Φ_d` decomposition matrix of `L`.

```julia-repl
julia> W=rootdatum(:psu,6)
psu₆

julia> L=reflection_subgroup(W,[1,2,4,5])
psu₆₍₁₂₅₄₎=(A₂A₂)₍₁₂₄₃₎Φ₁

julia> InducedDecompositionMatrix(L,W,6)
Induced Φ₆-decomposition matrix from psu₆₍₁₂₅₄₎=(A₂A₂)₍₁₂₄₃₎Φ₁ to psu₆

    │ps ps A₂
────┼─────────
²A₅ │ .  .  .
.3  │ 1  .  .
3.  │ 1  .  .
.21 │ 1  1  .
1.2 │ 2  1  .
21. │ 1  1  .
2.1 │ 2  1  .
.111│ .  1  1
111.│ .  1  1
1.11│ 1  2  1
11.1│ 1  2  1
```
The matrix itself is stored in the field `.scalar` of the returned `struct`.
"""
function InducedDecompositionMatrix(R,W,d::Integer)
  m=generic_decomposition_matrix(R,d)
  t=hc_induction_table(R,W)
  InducedDecompositionMatrix(t.scalar*m.scalar,m.colnames,W,R,d)
end

function Base.show(io::IO, ::MIME"text/html", m::InducedDecompositionMatrix)
  show(IOContext(io,:TeX=>true),"text/plain",m)
end

function Base.show(io::IO,m::InducedDecompositionMatrix)
  if !hasdecor(io)
    print(io,"InducedDecompositionMatrix(",m.R,",",m.W,",",m.d,")")
  else
    printTeX(io,"Induced \\Phi_{",m.d,"}-decomposition matrix from ",
          m.R," to ",m.W,"\n")
  end
end

function Base.show(io::IO,::MIME"text/plain",m::InducedDecompositionMatrix)
  println(io,m)
  scal=map(e->(!ismissing(e) && iszero(e)) ? "." : TeX(io,e),m.scalar)
  row_labels=charnames(io,UnipotentCharacters(m.W))
  order=sort(axes(m.scalar,1),
              by=i->maximum(findall(!iszero,m.scalar[i,:]);init=0))
  showtable(io,scal;row_labels,col_labels=m.colnames,rows=order)
end
