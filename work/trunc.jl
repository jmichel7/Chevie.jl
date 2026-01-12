# truncated formal series with l terms
struct Trunc{T}
  c::Vector{T} # vector of length l
  v::Int
end
  
Base.length(p::Trunc)=length(p.c)

function Base.:*(p::Trunc,q::Trunc)
  l=min(length(p),length(q))
  Trunc(map(j->sum(k->p[k+p.v]*q[j-k+q.v],0:j),0:l-1),p.v+q.v)
end

function Base.:+(p::Trunc,q::Trunc)
  if p.v>q.v p,q=q,p end
  l=min(length(p),length(q))
  Trunc(map(i->p[i]+q[i],p.v:p.v+l-1),p.v)
end

Base.:-(p::Trunc,q::Trunc)=p+(-q)

Trunc(p::Pol,i)=Trunc(p[p.v:p.v+i-1],p.v)
LaurentPolynomials.Pol(p::Trunc)=Pol(p.c,p.v)

@inbounds Base.getindex(p::Trunc,i::Integer)=
  i in p.v:p.v+length(p.c)-1 ? p.c[i-p.v+1] : zero(p.c[1])

Base.:^(a::Trunc, n::Integer)=Base.power_by_squaring(a,n)
Base.:-(p::Trunc)=Trunc(-p.c,p.v)
Base.one(p::Trunc)=Trunc(vcat(1,fill(0,length(p)-1)),0)
Base.copy(p::Trunc)=Trunc(copy(p.c),p.v)

function Base.show(io::IO, ::MIME"text/plain", a::Trunc)
  if !haskey(io,:typeinfo) 
    print(io,"Trunc($(length(a))): ")
    io=IOContext(io,:typeinfo=>typeof(a))
  end
  show(io,a)
end

function Base.show(io::IO,p::Trunc{T})where T
  if !get(io,:limit,false) && !get(io,:TeX,false) && !get(io,:naive,false)
    if ismonomial(p) && isone(p.c[1]) && p.v==1 && T==Int print(io,"Trunc()")
    else print(io,"Trunc(",p.c)
      if !iszero(p.v) print(io,",",p.v) end
      print(io,")")
    end
  else
    var="x"
    for deg in p.v:p.v+length(p)-1
      c=p[deg]
      c=repr(c; context=IOContext(io,:typeinfo=>typeof(c)))
      if !iszero(deg)
        c=format_coefficient(c)*var
        if get(io,:naive,false) c*=stringexp(stdout,deg)
        else c*=stringexp(io,deg)
        end
      end
      if c[1]!='-' && deg!=p.v c="+"*c end
      print(io,c)
    end
  end
end

function Trunc(Q::Frac{<:Pol},i)
  p=numerator(Q);q=denominator(Q)
  p//=Pol([q.c[1]],valuation(q))
  q=Trunc(Pol(q.c.//q.c[1])-1,i)
  Trunc(p,i)*sum(j->(-q)^j,0:i)
end

function Trunc(Q::Frac{<:Mvp},i,var::Symbol)
  p=numerator(Q);q=denominator(Q)
  Trunc(Pol(p,var)//Pol(q,var),i)
end
