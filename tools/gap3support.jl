using Primes: totient # translation of Phi
Append(a::Vector,b::AbstractVector)=vcat(a,b)
Append(a::String,b::String)=a*b
Append(a::String,b::Vector{Char})=a*String(b)
ApplyFunc(f,x)=f(x...)
CollectBy(v,f)=collectby(f,v)
Concatenation(b...)=vcat(b...)
Concatenation(a::String...)=prod(a)
function Concatenation(a::AbstractVector)
  all(x->x isa AbstractVector,a) ? vcat(a...) : prod(a)
end
ConcatenationString(s...)=prod(s)
DiagonalOfMat(m)=[m[i,i] for i in axes(m,1)]
Difference(a,b)=sort(setdiff(a,b))
function Position(a::Vector,b)
  x=findfirst(isequal(b),a)
  isnothing(x) ? false : x
end
First(a,b)=a[findfirst(b,a)]::eltype(a)
Flat(v)=collect(Iterators.flatten(v))
gapSet(v)=unique!(sort(v)) # set is transpiled to gapSet
IdentityMat(n)=map(i->one(zeros(Int,n,n))[i,:],1:n)
Intersection(x,y)=sort(intersect(x,y))
Intersection(x::Vector)=sort(intersect(x...))
IsInt(l)=l isa Int ||(l isa Rational && denominator(l)==1)
IsList(l)=l isa Vector
IsString(l)=l isa String
Join(x,y)=join(x,y)
Join(x)=join(x,",")
Lcm(a...)=Lcm(collect(a))
Lcm(a::Vector)=lcm(Int.(a))
Minimum(a::Number,x...)=min(a,x...)
Minimum(v::AbstractVector)=minimum(v)
NullMat(i,j=i)=[zeros(Int,j) for k in 1:i]
OnTuples(a,b)=a.^b
function pad(s, i::Int)
  if i>0 return lpad(string(s),i)
  else return rpad(string(s),-i)
  end
end
pad(s::String)=s
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
function RootInt(n,k=2)
  res=floor(Int,n^(1/k))
  if (res+1)^k<=n res+1 else res end
end
Rotations(a)=circshift.(Ref(a),length(a):-1:1)
SortBy(x,f)=sort!(x,by=f)
function SortParallel(a,b)
  v=sortperm(a)
  b.=b[v]
  a.=a[v]
end
SPrint=string
Sublist(a::Vector, b::AbstractVector)=a[b]
Sum(v::AbstractVector)=sum(v)
Sum(v::AbstractVector,f)=isempty(v) ? 0 : sum(f,v)
TransposedMat(l::Vector{<:Vector})=collect.(zip(l...))
ValuePol(v,c)=isempty(v) ? 0 : evalpoly(c,v)
