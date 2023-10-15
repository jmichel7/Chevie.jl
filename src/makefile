JULIA=1.9.3
test: *.jl
	~/julia-$(JULIA)/bin/julia ../tools/gentest.jl *.jl ../contr/cp.jl
	mv runtests.jl ../test
	~/julia-$(JULIA)/bin/julia ../test/runtests.jl
some: $(FILES)
	~/julia-$(JULIA)/bin/julia ../tools/gentest.jl $(FILES)
	~/julia-$(JULIA)/bin/julia runtests.jl
	rm runtests.jl
