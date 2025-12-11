test: *.jl
	julia ../tools/gentest.jl *.jl
	mv runtests.jl ../test
	julia ../test/runtests.jl
some: $(FILES)
	julia ../tools/gentest.jl $(FILES)
	julia runtests.jl
	rm runtests.jl
