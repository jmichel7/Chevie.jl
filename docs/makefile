JULIA=1.10.2
doc:
	~/julia-$(JULIA)/bin/julia make.jl local
	cp build/*.html .
	cp build/search_index.js .
	cp -r build/* ~/cheviedoc
