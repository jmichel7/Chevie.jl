# extensions to get closer to GAP semantics
Base.:*(a::AbstractArray,b::Pol)=a .* b
Base.:*(a::Pol,b::AbstractArray)=a .* b
Base.:*(a::AbstractVector,b::AbstractVector)=toL(toM(a)*toM(b))
Base.:*(a::AbstractVector{<:Number},b::AbstractVector)=toL(permutedims(a)*toM(b))[1]
Base.:*(a::Tuple,b::AbstractVector)=toL(permutedims(collect(a))*toM(b))[1]
Base.:*(a::AbstractVector{Pol},b::AbstractVector{Pol})=sum(a.*b)
Base.:+(a::AbstractArray,b::Pol)=a .+ b
Base.:/(a::AbstractArray,b::Pol)=a ./ b
Base.://(a::AbstractArray,b::Pol)=a .// b
Base.:*(a::Array,b::Mvp)=a .* b
Base.:*(a::Mvp,b::Array)=a .* b
Base.:+(a::AbstractArray,b::Mvp)=a .+ b
Base.:/(a::AbstractArray,b::Mvp)=a ./ b
Base.://(a::AbstractArray,b::Mvp)=a .// b
Cycs.:^(a::Cyc,b::Rational)=a^Int(b)
Base.:^(m::AbstractMatrix,n::AbstractMatrix)=inv(n*E(1))*m*n
Base.:^(m::Vector{<:Vector{<:Number}},n::Matrix{<:Number})=inv(n)*toM(m)*n
Base.:^(m::Vector,n::Vector)=toL(inv(toM(n)*E(1)//1)*toM(m)*toM(n))
Base.inv(m::Vector)=toL(inv(toM(m).+zero(E(1)//1)))
Base.:(//)(m::Vector,n::Vector)=toL(toM(m)*inv(toM(n)*E(1)))
Base.getindex(a::Symbol,i::Int)=string(a)[i]
Base.length(a::Symbol)=length(string(a))
Base.union(v::Vector)=union(v...)

# correct translations of GAP3 functions
ApplyWord(w,gens)=isempty(w) ? one(gens[1]) : prod(i->i>0 ? gens[i] : inv(gens[-i]),w)
BetaSet=βset
CartanMat(s,a...)=cartan(Symbol(s),a...)
CharParams(W)=charinfo(W)[:charparams]
function CharRepresentationWords(mats,words)
  traces_words_mats(toM.(mats),words)
end
function Collected(v)
  d=groupby(v,v)
  sort([[k,length(v)] for (k,v) in d])
end
function CollectBy(v,f)
  d=groupby(f,v)
  [d[k] for k in sort(collect(keys(d)))]
end
ConcatenationString(s...)=prod(s)
CoxeterWord(W,w)=word(W,w)
Cycles(p,i)=orbits(p,i)
CycPolFakeDegreeSymbol=fegsymbol
DefectSymbol=defectsymbol
function DiagonalMat(v...)
  R=cat(map(m->m isa Array ? m : hcat(m),v)...,dims=(1,2))
  for i in axes(R,1), j in axes(R,2)
    if i!=j R[i,j]=zero(R[1,1]) end
  end
  toL(R)
end
DiagonalMat(v::Vector{<:Number})=DiagonalMat(v...)
Drop(a::AbstractVector,i::Int)=deleteat!(collect(a),i)
Elements=elements
EltWord(W,x)=W(x...)
ExteriorPower(m,i)=toL(exterior_power(toM(m),i))
Factors(n)=reduce(vcat,[fill(k,v) for (k,v) in factor(n)])
FullSymbol=fullsymbol
GaloisCyc=galois
Hasse=hasse
HighestPowerFakeDegreeSymbol=degree_feg_symbol
HighestPowerGenericDegreeSymbol=degree_gendeg_symbol
function Ignore() end
InfoChevie=println
InfoChevie2=print
Inherit(a,b)=merge!(a,b)
function Inherit(a,b,c)
  for k in c a[Symbol(k)]=b[Symbol(k)] end
  a
end
Intersection(x,y)=sort(intersect(x,y))
Intersection(x::Vector)=sort(intersect(x...))
IntListToString=joindigits
Join(x,y)=join(x,y)
Join(x)=join(x,",")
KroneckerProduct(a,b)=toL(kron(toM(a),toM(b)))
LowestPowerFakeDegreeSymbol=valuation_feg_symbol
LowestPowerGenericDegreeSymbol=valuation_gendeg_symbol
MatXPerm=refrep
NrConjugacyClasses(W)=length(classinfo(W)[:classtext])
OnMatrices(a::Vector{<:Vector},b::Perm)=(a.^b)^b
PartBeta=partβ
PermutationMat=Matrix
Rank=rank
RankSymbol=ranksymbol
ReflectionSubgroup(W,I::AbstractVector)=reflection_subgroup(W,convert(Vector{Int},I))
RootInt(a,b=2)=floor(Int,a^(1/b)+0.0001)
RootsCartan=roots
Rotations(a)=circshift.(Ref(a),length(a):-1:1)
SchurFunctor(m,p)=toL(schur_functor(toM(m),p))
SignedPerm=SPerm
SymmetricDifference(x,y)=sort(symdiff(x,y))
SymmetricPower(m,n)=SchurFunctor(m,[n])
SemisimpleRank(W)=semisimplerank(W)
function SortParallel(a,b)
  v=sortperm(a)
  b.=b[v]
  a.=a[v]
end
ShiftBeta=shiftβ
StringSymbol=stringsymbol
StringToDigits(s)=map(y->Position("01234567890", y), collect(s)).-1
SymbolPartitionTuple=symbol_partition_tuple
SymbolsDefect(a,b,c,d)=symbols(a,b,d)
Tableaux=tableaux
function TeXBracket(s)
  s=string(s)
  length(s)==1  ? s : "{"*s*"}"
end
Torus(i::Int)=torus(i)
Value(p,v)=p(v)
ValuePol(v,c)=isempty(v) ? 0 : evalpoly(c,v)
function CoxeterGroup(S::String,s...)
 if length(s)==1 return coxgroup(Symbol(S),Int(s[1])) end
 coxgroup(Symbol(S),Int(s[1]))*coxgroup(Symbol(s[2]),Int(s[3]))
end
CoxeterGroup()=coxgroup()
#-------------------------------------------------------------------------
#  dummy translations of GAP3 functions
Format(x)=string(x)
FormatTeX(x)=repr(x,context=:TeX=>true)
FormatGAP(x)=replace(repr(x)," "=>"")
Format(x,opt)=sprint((io,x)->show(IOContext(io,opt...),x),x)

function ReadChv(s) end
Groups.Group(a::Perm...)=Group(collect(a))
GetRoot(x::Cyc,n::Number=2,msg...)=root(x,n)
GetRoot(x::Integer,n::Number=2,msg...)=root(x,n)
GetRoot(x::Pol,n::Number=2,msg...)=root(x,n)
GetRoot(x::Rational,n::Number=2,msg...)=root(x,n)
GetRoot(x::Mvp,n::Number=2,msg...)=root(x,n)
function GetRoot(x,n::Number=2,msg...)
  error("GetRoot($x,$n) not implemented")
end

Unbind(x)=x
#-------------------------------------------------------------------------
