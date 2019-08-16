using Documenter, Gapjm

makedocs(sitename="Gapjm.jl documentation",
         modules=[Cycs,Gapjm,PermGroups,Perms,Pols,CoxGroups,
                  Weyl,Hecke,KL,Garside,Util,CycPols])
