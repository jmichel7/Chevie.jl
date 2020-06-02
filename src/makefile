test: *.jl
	~/julia-1.4.2/bin/julia ../tools/gentest.jl *.jl
#	~/julia-1.4.2/bin/julia ../tools/gentest.jl PermGroups.jl
	mv runtests.jl ../test
	~/julia-1.4.2/bin/julia ../test/runtests.jl
