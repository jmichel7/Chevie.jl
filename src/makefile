JULIA=~/julia-1.7.2/bin/julia
test: *.jl
	$(JULIA) ../tools/gentest.jl *.jl
	mv runtests.jl ../test
	$(JULIA) ../test/runtests.jl
some: $(FILES)
	$(JULIA) ../tools/gentest.jl $(FILES)
	$(JULIA) runtests.jl
	rm runtests.jl
