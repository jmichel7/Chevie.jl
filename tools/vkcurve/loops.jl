# Ordonne une liste de points trigonometriquement autour d'un centre
function cycorder(list, center)
  right=empty(list)
  left=empty(list)
  top=empty(list)
  bottom=empty(list)
  for y in list
    if real(y)>real(center) push!(right, y)
    elseif real(y)<real(center) push!(left, y)
    elseif imag(y)>imag(center) push!(top, y)
    else push!(bottom, y)
    end
  end
  sort!(right,by=x->imag(x-center)/real(x-center))
  sort!(left,by=x->imag(x-center)/real(x-center))
  vcat(right, top, left, bottom)
end

# Input: (list of complex numbers,complex number)
# Output: sublist of "neighbours" of the second input,
#   x and y are neighbours iff no z is in the disk of diameter [x,y]
function neighbours(l, center)
  function isneighbour(y)
    d=abs2(y-center)
    for z in l
      if z!=y && abs2(y-z)+abs2(z-center)<=d return false end
    end
    return true
  end
  l=filter(y->y!=center,l)
  l=filter(isneighbour, l)
  cycorder(l, center)
end

# value at z of an equation of the line (x,y)
function lineq(x, y, z)
  if real(x)≈real(y)
    if imag(x)≈imag(y) error("Undefined line\n") 
    else return real(z)-real(x)
    end
  else 
    return (imag(y)-imag(x))*(real(z)-real(x))/(real(y)-real(x))+imag(x)-imag(z)
  end
end

function mediatrix(x, y)
  if x≈y error("Undefined mediatrix") end
  (x+y)/2 .+[im,-im]*(x-y)
end

crossing(v1,v2)=crossing(v1...,v2...)

# Computes the intersecting point of two lines
# returns false if the lines are parallel or a pair is a single
function crossing(x1,x2,y1,y2)
  if x1≈x2 || y1≈y2 return nothing end
  if !(real(x1)≈real(x2))
    lambdax=(imag(x1)-imag(x2))/(real(x1)-real(x2))
    mux=-lambdax*real(x1)+imag(x1)
    if !(real(y1)≈real(y2))
      lambday=(imag(y1)-imag(y2))/(real(y1)-real(y2))
      muy=-lambday*real(y1)+imag(y1)
      if 1+lambdax≈1+lambday return nothing end
      resr=(muy-mux)/(lambdax-lambday)
      resi=lambdax*resr+mux
      res=Complex(resr, resi)
    else
      E3=Complex{Float64}(E(3))
      res=crossing(E3*x1, E3*x2, E3*y1, E3*y2)
      if isnothing(res) return nothing end
      res/=E3
    end
  else
    res=crossing(im*x1, im*x2, im*y1, im*y2)
    if isnothing(res) return nothing end
    res/=im
  end
  res
end

function detectsleftcrossing(c, w, y, z)
  res=fill(false,length(c)-1)
  a,b=mediatrix(y, z)
  for k in 1:length(c)-1
    if lineq(a, b, c[k])*lineq(a, b, c[k+1])<=0
      x=crossing(a, b, c[k], c[k+1])
      if !isnothing(x) xx=(z-y)/(w[k]-y)
         res[k]=imag(xx)>=0
      end
    end
  end
  res
end

# y must be an element of ys
function boundpaths(ys, sy, path, y)
  if !haskey(y, :path)
    y[:path]=vcat(path, [y[:y]])
    for z in y[:lovers] boundpaths(ys, sy, y[:path], sy(z)) end
  end
end

# eliminates trivial segments and contracts pairs [a,b],[b,a]
function shrink(l)
  k=findfirst(i->l[i]==l[i+1],1:length(l)-1)
  if !isnothing(k) return shrink(l[vcat(1:k, k+2:length(l))])
  else
    k=findfirst(i->l[i]==l[i+2],1:length(l)-2)
    if !isnothing(k) return shrink(l[vcat(1:k, k+3:length(l))])
    else return l
    end
  end
end

# converts old result format (list of loops) to the new format :
#  .points   : set of all endpoints of all segments
#  .segments : set of all segments used, where endpoints are indexed
#              as in points
#     .loops : list of sequence of numbers of used segments
function convert(oldres)
  res=copy(oldres)
  points=sort(unique(vcat(oldres...)),by=x->(imag(x),real(x)))
  np(p)=findfirst(==(p),points)
  loops=map(l->map(i->np.(l[i-1:i]),2:length(l)),oldres)
  segments=sort(unique(vcat(map(l->sort.(l),loops)...)))
  loops=map(loops)do l
    map(l)do seg
      p=findfirst(==(seg),segments)
      isnothing(p) ? -findfirst(==(reverse(seg)),segments) : p
    end
  end
  Dict{Symbol, Any}(:points=>points, :segments=>segments, :loops=>loops)
end

"""
'LoopsAroundPunctures(points)'

The input is  a list of complex rational numbers.  The function computes
piecewise-linear loops representing generators  of the fundamental group
of `ℂ -{points}`.

```julia-repl
julia> LoopsAroundPunctures([0])
Dict{Symbol, Any} with 3 entries:
  :segments => [[1, 2], [1, 3], [2, 4], [3, 4]]
  :points   => Complex{Int64}[0-1im, -1+0im, 1+0im, 0+1im]
  :loops    => [[4, -3, -1, 2]]
```

The output is a record with  three fields. 

`points`: a list of complex  numbers. 
`segments`:  a list of oriented segments, each of them  encoded by the list 
   of the positions in 'points' of  its two endpoints. 
`loops`: a list of list of  integers. Each list  of integers represents a 
   piecewise linear loop, obtained  by concatenating the  elements of 
   `segments`  indexed by the integers (a negative integer is used when the
   opposed orientation of the segment has to be taken).
"""
function LoopsAroundPunctures(originalroots)
  roots=originalroots
  n=length(roots)
  average=sum(roots)/n
  sort!(roots, by=x->abs2(x-average))
  if n==1 return convert([roots[1].+[1,im,-1,-im,1]]) end
  ys=map(x->Dict{Symbol, Any}(:y=>x), roots)
  sy(y)=ys[findfirst(==(y),roots)]
  for y in ys
    y[:neighbours]=neighbours(roots, y[:y])
    y[:friends]=[y[:y]]
    y[:lovers]=empty(y[:friends])
  end
  if VKCURVE[:showLoops] println("neighbours computed") end
  for y in ys
    for z in y[:neighbours]
      if !(z in y[:friends])
        push!(y[:lovers], z)
        push!(sy(z)[:lovers], y[:y])
        newfriends=vcat(y[:friends], sy(z)[:friends])
        for t in y[:friends] sy(t)[:friends]=newfriends end
        for t in sy(z)[:friends] sy(t)[:friends]=newfriends end
      end
    end
  end
  for y in ys sort!(y[:neighbours],by=z->abs2(y[:y]-z)) end
  rs=map(y->real(y[:y]), ys)
  is=map(y->imag(y[:y]), ys)
  minr=real(ys[argmin(rs)][:y])
  maxr=real(ys[argmax(rs)][:y])
  mini=imag(ys[argmin(is)][:y])
  maxi=imag(ys[argmax(is)][:y])
# To avoid trouble with points on the border of the convex hull,
# we make a box around all the points;
  box=[Complex(minr-2, mini-2), Complex(minr-2, maxi+2),
       Complex(maxr+2, mini-2), Complex(maxr+2, maxi+2),
       Complex((maxr+minr)/2, mini-(maxr-minr)/2-2),
       Complex((maxr+minr)/2, maxi+(maxr-minr)/2+2),
       Complex(minr-(maxi-mini)/2-2, (maxi+mini)/2),
       Complex(maxr+(maxi-mini)/2+2, (maxi+mini)/2)]
  for y in ys
    y[:cycorder]=cycorder(vcat(filter(z->z!=y[:y],roots),box),y[:y])
    k=findfirst(==(y[:neighbours][1]),y[:cycorder])
    y[:cycorder]=vcat(y[:cycorder][k:end],y[:cycorder][1:k-1])
    push!(y[:cycorder], y[:cycorder][1])
    y[:circle]=Complex{Float64}[(y[:y]+y[:neighbours][1])/2]
    y[:witness]=Complex{Float64}[y[:neighbours][1]]
    for z in y[:cycorder][2:end]
      cut=detectsleftcrossing(y[:circle], y[:witness], y[:y], z)
      if true in cut
        k=findfirst(cut)
        y[:circle]=y[:circle][1:k]
        y[:witness]=y[:witness][1:k]
      end
      k=length(y[:circle])
      newcirc=crossing(mediatrix(y[:y],y[:witness][k]),mediatrix(y[:y], z))
      if !isnothing(newcirc)
        push!(y[:circle], newcirc)
        push!(y[:witness], z)
      end
      if z in y[:lovers]
        push!(y[:circle], (y[:y]+z)/2)
        push!(y[:witness], z)
      end
    end
  end
  if VKCURVE[:showLoops] println("circles computed") end
  boundpaths(ys, sy, [], ys[1])
  for y in ys
    k=length(y[:path])
    if k>1
      circleorigin=(y[:y]+y[:path][k-1])/2
      k=findfirst(≈(circleorigin),y[:circle])
      y[:circle]=vcat(y[:circle][k:end],y[:circle][1:k-1])
    end
  end
  for y in ys
    k=length(y[:path])
    y[:handle]=Vector{Complex{Float64}}(vcat(map(1:k-1)do i
      l=sy(y[:path][i])[:circle]
      l[1:findfirst(≈((y[:path][i]+y[:path][i+1])/2),l)]
     end...))
    y[:loop]=vcat(y[:handle], y[:circle], reverse(y[:handle]))
  end
  sort!(ys,by=y->findfirst(==(y[:y]),originalroots))
  for y in ys 
    y[:loop]=map(x->round(x;sigdigits=8),y[:loop])
    y[:loop]=shrink(y[:loop])
  end
  convert(map(y->y[:loop], ys))
end

function segs(r,v)
  rp=Float64[]
  ip=Float64[]
  for seg in v
    s=r[:points][seg<0 ? reverse(r[:segments][-seg]) : r[:segments][seg]]
    append!(rp,real.(s))
    append!(ip,imag.(s))
  end
  rp,ip
end

function loops(r,v)
  rp=Float64[]
  ip=Float64[]
  for l in r[:loops][v]
    nrp,nip=segs(r,l)
    append!(rp,nrp)
    append!(ip,nip)
  end
  rp,ip
end
