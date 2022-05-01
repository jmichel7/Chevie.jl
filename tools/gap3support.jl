Append(a::Vector,b::AbstractVector)=vcat(a,b)
Append(a::String,b::String)=a*b
Append(a::String,b::Vector{Char})=a*String(b)
ApplyFunc(f,x)=f(x...)
ComplexConjugate(v)=conj(v)
Concatenation(a::String...)=prod(a)
Concatenation(a::AbstractVector{<:AbstractVector})=vcat(collect.(a)...)
Concatenation(a::AbstractVector{String})=prod(a)
Concatenation(a::Vector,b::Tuple)=vcat(a,improve_type(collect(b)))
Concatenation(b...)=reduce(vcat,b)
DiagonalOfMat(m)=[m[i,i] for i in axes(m,1)]
Difference(a,b)=sort(setdiff(a,b))
function Position(a::Vector,b)
  x=findfirst(isequal(b),a)
  isnothing(x) ? false : x
end
First(a,b)=a[findfirst(b,a)]::eltype(a)
Flat(v)=collect(Iterators.flatten(v))
gapSet(v)=unique!(sort(v))
IdentityMat(n)=map(i->one(rand(Int,n,n))[i,:],1:n)
IsInt(l)=l isa Int ||(l isa Rational && denominator(l)==1)
IsList(l)=l isa Vector
IsString(l)=l isa String
Lcm(a...)=Lcm(collect(a))
Lcm(a::Vector)=lcm(Int.(a))
Minimum(v::AbstractVector)=minimum(v)
Minimum(a::Number,x...)=min(a,x...)
NullMat(i,j=i)=[zeros(Int,j) for k in 1:i]
OnTuples(a,b)=a.^b
PermListList(l1,l2)=Perm(improve_type(l1),improve_type(l2))
Permuted(x,p)=x^p
function Position(a::String,b::String)
  x=findfirst(b,a)
  isnothing(x) ? false : x.start
end
function Position(a::String,b::Char)
  x=findfirst(isequal(b),a)
  isnothing(x) ? false : x
end
function PositionProperty(a::AbstractVector,b::Function)
  r=findfirst(b,a)
  if isnothing(r) return false end
  r
end
Product(v)=isempty(v) ? 1 : prod(v)
Product(v::AbstractVector,f)=isempty(v) ? 1 : prod(f,v)
SortBy(x,f)=sort!(x,by=f)
SPrint=string
Sublist(a::Vector, b::AbstractVector)=a[b]
Sum(v::AbstractVector)=sum(v)
Sum(v::AbstractVector,f)=isempty(v) ? 0 : sum(f,v)
TransposedMat(l::Vector{<:Vector})=collect.(zip(l...))

