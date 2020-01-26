using Documenter, Gapjm

makedocs(sitename="Gapjm.jl documentation",
         modules=[Gapjm,Perms,Groups,PermGroups,Cycs,Pols,Mvps,CoxGroups,
                  Weyl,PermRoot,HeckeAlgebras,KL,Garside,Chars,Cosets,Uch,
                  Ucl,Symbols,Util,ModuleElts,CycPols,Posets])
