using Documenter, Chevie

DocMeta.setdocmeta!(Chevie, :DocTestSetup, :(using Chevie); recursive=true)

makedocs(sitename="Chevie.jl documentation",
   modules=[Chevie],
   authors="Jean Michel <jean.michel@imj-prg.fr> and contributors",
   format = Documenter.HTML(;
       canonical = "https://juliadocs.github.io/Chevie.jl",
       edit_link="main",
       assets=String[],
    ),
    pages=[
    "index.md",
    "Infrastructure"=>[
    "format.md",
    "symbols.md",
    "nf.md"],
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
    "dict.md" ],
  warnonly=:missing_docs,
)

deploydocs(;
    repo="github.com/jmichel7/Chevie.jl",
    devbranch="main",
)
