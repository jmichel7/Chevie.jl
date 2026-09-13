test: *.jl
	julia $(VERSION) ../tools/gentest.jl *.jl
	mv runtests.jl ../test
	julia $(VERSION) ../test/runtests.jl
some: $(FILES)
	julia ../tools/gentest.jl $(FILES)
	julia runtests.jl
	rm runtests.jl
