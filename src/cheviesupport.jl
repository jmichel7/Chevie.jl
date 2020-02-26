# extensions to get closer to GAP semantics
Base.:*(a::Array,b::Pol)=a .* b
Base.:*(a::Pol,b::Array)=a .* b
Base.:*(a::AbstractVector,b::AbstractVector)=toL(toM(a)*toM(b))
Base.:*(a::AbstractVector{<:Number},b::AbstractVector)=toL(permutedims(a)*toM(b))[1]
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
Base.:^(m::Vector,n::Vector)=toL(inv(toM(n)*E(1))*toM(m)*toM(n))
Base.inv(m::Vector)=toL(inv(toM(m)*E(1)//1))
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
Hasse=hasse
HighestPowerFakeDegreeSymbol=degree_feg_symbol
HighestPowerGenericDegreeSymbol=degree_gendeg_symbol
function Ignore() end
InfoChevie2=print
IntListToString=joindigits
Join(x,y)=join(x,y)
Join(x)=join(x,",")
LowestPowerFakeDegreeSymbol=valuation_feg_symbol
LowestPowerGenericDegreeSymbol=valuation_gendeg_symbol
MatXPerm=matX
NrConjugacyClasses(W)=length(classinfo(W)[:classtext])
OnMatrices(a::Vector{<:Vector},b::Perm)=(a.^b)^b
PartBeta=partβ
Rank=rank
RankSymbol=ranksymbol
ReflectionSubgroup(W,I::AbstractVector)=reflection_subgroup(W,convert(Vector{Int},I))
RootInt(a,b)=floor(Int,a^(1/b)+0.0001)
RootsCartan=roots
Rotations(a)=circshift.(Ref(a),length(a):-1:1)
SemisimpleRank(W)=semisimplerank(W)
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
function CoxeterGroup(S::String,s...)
 if length(s)==1 return coxgroup(Symbol(S),Int(s[1])) end
 coxgroup(Symbol(S),Int(s[1]))*coxgroup(Symbol(s[2]),Int(s[3]))
end
CoxeterGroup()=coxgroup()
#-------------------------------------------------------------------------
#  dummy translations of GAP3 functions
Format(x)=string(x)
FormatTeX(x)=repr(x,context=:TeX=>true)
FormatGAP(x)=repr(x)
ff(io::IO,p;opt...)=show(IOContext(io,opt...),p)
Format(x,opt)=sprint((io,x)->ff(io,x;opt...),x)

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
