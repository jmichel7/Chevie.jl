JULIA=1.11.3
doc:
	~/julia-$(JULIA)/bin/julia make.jl local
	cp build/*.html .
	cp build/search_index.js .
	cp -r build/* ~/cheviedoc
