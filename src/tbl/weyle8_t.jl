chevieset(:E8, :WGraph, function(i)
  gr=chevieget(:E8, :WGraphs)
  if gr[i]===nothing
    gr[i]=DualWGraph(8,chevieget(:E8,:WGraph)(i-1))
  elseif gr[i] isa String
    head=@__FILE__
    head=head[1:findlast(==('/'),head)]
    include(head*"e8wgraph/rep"*gr[i]*".jl")
  else
  end
  gr[i]
end)
