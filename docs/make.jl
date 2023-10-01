using Documenter, Gapjm

makedocs(sitename="Gapjm.jl documentation",
         modules=[Gapjm],
           format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        #canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
        collapselevel=1
    ),
pages=[
"index.md",
"Infrastructure"=>[
"symbols.md",
"nf.md",
"presentations.md"],
"Reflection groups"=>[
"permroot.md",
"coxgroups.md",
"weyl.md",
"chars.md"],
"Hecke algebras"=>[
"hecke.md",
"kl.md",
"algebras.md"],
"garside.md",
"Reductive groups"=>[
"semisimple.md",
"cosets.md",
"sscoset.md",
"uch.md",
"ct.md"],
"Eigenspaces"=>[
"eigen.md",
"dseries.md"],
"Unipotent elements"=>[
"ucl.md",
"urad.md"],
"gendec.md",
"dict.md" ]
)
