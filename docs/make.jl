using Documenter, Gapjm, Primes

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
"cycpol.md",
"posets.md",
"sperm.md",
"glinearalgebra.md",
"matint.md",
"ffe.md",
"presentations.md"],
"Reflection groups"=>[
"coxgroups.md",
"weyl.md",
"permroot.md",
"chars.md"],
"Hecke algebras"=>[
"hecke.md",
"kl.md"],
"garside.md",
"Reductive groups"=>[
"semisimple.md",
"cosets.md",
"sscoset.md",
"uch.md"],
"eigen.md",
"dseries.md",
"symbols.md",
"Unipotent elements"=>[
"ucl.md",
"urad.md"],
"ct.md",
"gendec.md",
"dict.md" ]
)
