"""
`NewtonRoot(p::Pol,initial,precision;showall=false,show=false,lim=800)`

Here `p` is a complex polynomial. The function computes an approximation to a
root of `p`, guaranteed of distance closer than `precision` to
an  actual root. The first approximation used is `initial`. If `initial` is
in  the  attraction  basin  of  a  root  of  `p`,  the  one approximated. A
possibility  is that  the Newton  method starting  from `initial`  does not
converge  (the  number  of  iterations  after  which  this  is  decided  is
controlled  by  `lim`);  then  the  function returns `nothing`.
Otherwise  the function returns  a pair: the  approximation found, and an
upper  bound of the distance between that approximation and an actual root.
The point of returning
an upper bound is that it is usually better than the asked-for `precision`.
For the precision estimate a good reference is cite{HSS01}.

```julia-repl
julia> p=Pol([1,0,1])
Pol{Int64}: x²+1

julia> NewtonRoot(p,1+im,10^-7)
2-element Vector{ComplexF64}:
  8.463737877036721e-23 + 1.0im
 1.8468154345305913e-11 + 0.0im

julia> NewtonRoot(p,1,10^-7;show=true)
****** Non-Convergent Newton after 800 iterations ******
p=x²+1 initial=-1.0 prec=1.0000000000000004e-7
```
"""
function NewtonRoot(p::Pol,z,precision;showall=VKCURVE[:showallnewton],
                          show=VKCURVE[:showNewton],lim=VKCURVE[:NewtonLim])
  deriv = derivative(p)
  for cnt in 1:lim
    a=p(z)
    b=deriv(z)
    c=iszero(b) ? a : a/b
    err=abs(BigNorm(c))
    if iszero(err) err=(precision/100)/(degree(p)+1) end
    if err>precision err=precision end
#   pr=maximum(0,-1-DecimalLog(err))
#   z=ComplexRational(evalf(z-c,pr+1))
#   if showall println(cnt,": ",evalf(z, pr)) end
    z-=c
    if showall println(cnt,": ",z) end
#   if 10^-pr*(length(p)-1)<=precision
    if err*degree(p)<=precision
      if show print(cnt) end
#     return [z,10^-pr]
      return (z,err)
    end
  end
  if show
    println("\n****** Non-Convergent Newton after ", lim," iterations ******")
    xprintln("p=",p," initial=",z," prec=",precision)
  end
end

"""
'SeparateRootsInitialGuess(p, v, safety)'

Here  `p` is a complex  polynomial, and `v` is  a list of approximations to
roots  of `p` which should lie in different attraction basins for Newton' s
method.  The  result  is  a  list  `l`  of  complex  rationals representing
approximations  to the  roots of  `p` (each  element of  `l` is the root in
whose attraction basin the corresponding element of `v` lies), such that if
`d`  is the minimum distance  between two elements of  `l`, then there is a
root  of `p` within radius  `d/(2*safety)` of any element  of `l`. When the
elements  of  `v`  do  not  lie  in  different  attraction basins (which is
necessarily the case if `p` has multiple roots), 'false' is returned.

```julia-repl
julia> p=Pol([1,0,1])
Pol{Int64}: x²+1

julia> SeparateRootsInitialGuess(p,[1+im,1-im],10^5)
2-element Vector{ComplexF64}:
 8.463737877036721e-23 + 1.0im
 8.463737877036721e-23 - 1.0im

julia> SeparateRootsInitialGuess(p,[1+im,2+im],1000)
    # 1+im and 2+im in same attraction basins
```
"""
function SeparateRootsInitialGuess(p, v, safety)
  if degree(p)==1 return [-p[0]/p[1]] end
  radv=Dispersal(v)[1]/safety/2
  v=map(e->NewtonRoot(p,e,radv),v)
  if !any(isnothing,v) && Dispersal(first.(v))[1]/safety/2>=maximum(last.(v))
    return first.(v)
  end
end

"""
'SeparateRoots(<p>, <safety>)'

Here  `p` is  a complex  polynomial. The  result is  a list  `l` of complex
numbers  representing approximations to the roots  of `p`, such that if `d`
is  the minimum distance between two elements  of `l`, then there is a root
of  `p` within  radius `d/(2*safety)`  of any  element of  `l`. This is not
possible when `p` has multiple roots, in which case `nothing` is returned.

```julia-repl
julia> @Pol q
Pol{Int64}: q

julia> SeparateRoots(q^2+1,100)
2-element Vector{ComplexF64}:
  2.3541814200656927e-43 + 1.0im
 -2.3541814200656927e-43 - 1.0im

julia> SeparateRoots((q-1)^2,100)

julia> SeparateRoots(q^3-1,100)
3-element Vector{ComplexF64}:
 -0.5 - 0.8660254037844386im
  1.0 - 1.232595164407831e-32im
 -0.5 + 0.8660254037844387im
```
"""
function SeparateRoots(p,safety)
  subtractroot(p,r)=divrem(p,Pol([-r,1]))[1]
  if p isa Mvp p=Pol(p) end
  if degree(p)<1 return empty(p.c)
  elseif degree(p)==1 return [-p[0]/p[1]]
  end
  p/=p[end]
  e=Complex{Float64}(E(7))
  v=nothing
  cnt = 0
  while isnothing(v) && cnt<2*(degree(p)+1)
    if VKCURVE[:showNewton] && cnt>0
      println("****** ", cnt, " guess failed for p degree ", degree(p))
    end
    v=NewtonRoot(p,e,safety*10.0^-(degree(p)+4))
    e*=5/4*Complex{Float64}(E(2*(degree(p)+1)))
    cnt+=1
  end
  if cnt>=2*(degree(p)+1) error("no good initial guess") end
  v=[v[1]]
  append!(v,SeparateRoots(subtractroot(p,v[1]), safety))
  safety==0 ? v : SeparateRootsInitialGuess(p, v, safety)
end

"""
'FindRoots(<p>, <approx>)'

<p>  should be a univariate 'Mvp'  with cyclotomic or 'Complex' rational or
decimal  coefficients or  a list  of cyclotomics  or 'Complex' rationals or
decimals  which represents  the coefficients  of a  complex polynomial. The
function  returns  'Complex'  rational  approximations  to the roots of <p>
which  are  better  than  <approx>  (a  positive rational). Contrary to the
functions  'SeparateRoots', etc... described in  the previous chapter, this
function handles quite well polynomials with multiple roots. We rely on the
algorithms explained in detail in cite{HSS01}.

|    gap> FindRoots((x-1)^5,1/100000000000);
    [ 6249999999993/6250000000000+29/12500000000000I, 
      12499999999993/12500000000000-39/12500000000000I, 
      12500000000023/12500000000000+11/6250000000000I, 
      12500000000023/12500000000000+11/6250000000000I, 
      312499999999/312500000000-3/6250000000000I ]
    gap> evalf(last);
    [ 1, 1, 1, 1, 1 ]
    gap> FindRoots(x^3-1,1/10);            
    [ -1/2-108253175473/125000000000I, 1, -1/2+108253175473/125000000000I ]
    gap> evalf(last);
    [ -0.5-0.8660254038I, 1, -0.5+0.8660254038I ]
    gap> List(last,x->x^3);
    [ 1, 1, 1 ]|
"""
#FindRoots:=function(p,prec)local v,subtractroot,e;
#  subtractroot:=function(p,r)local d,res,i;
#    d:=Length(p);
#    res:=[];res[d-1]:=p[d];
#    for i in [d-1,d-2..2] do res[i-1]:=p[i]+res[i]*r;od;
#    return res;
#  end;
#  if IsMvp(p) then v:=Variables(p);
#    if Length(v)>1 then Error(p," is not univariate");
#    elif Length(v)=0 then v:=["x"];
#    fi;
#    p:=ScalMvp(Coefficients(p,v[1]));
#  fi;
#  if Length(p)=2 then return [-p[1]/p[2]];fi;
#  p:=List(p,ComplexRational);
#  e:=evalf(E(7),10);v:=false;
#  while v=false do
#    v:=NewtonRoot(p,ComplexRational(e),10^(-Length(p)-3))[1];
#    e:=e*evalf(E(Length(p)),10);
#  od;
#  v:=Concatenation([v],FindRoots(subtractroot(p,ComplexRational(v)),prec));
#  return List(v,e->NewtonRoot(p,ComplexRational(e),prec)[1]);
#end;
