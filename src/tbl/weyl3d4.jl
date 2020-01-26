
chevieset(Symbol("3D4"), :cyclestructure, [[2, nothing, nothing, nothing, 3], [nothing, 6], [3, nothing, nothing, nothing, 3], [3, 2, nothing, nothing, 2], [nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, 2], [nothing, 8], [nothing, nothing, nothing, nothing, 4]])
chevieset(Symbol("3D4"), :generators, [#= none:4 =# @perm_str("(1,13)(3,5)(6,8)(7,9)(10,11)(15,17)(18,20)(19,21)(22,23)"), #= none:5 =# @perm_str("( 2,14)( 3, 6)( 5, 8)( 7,10)( 9,11)(15,18)(17,20)(19,22)(21,23)"), #= none:6 =# @perm_str("( 1, 5)( 2, 6)( 3,15)( 4, 7)(11,12)(13,17)(14,18)(16,19)(23,24)"), #= none:7 =# @perm_str("( 3, 7)( 4,16)( 5, 9)( 6,10)( 8,11)(15,19)(17,21)(18,22)(20,23)")])
chevieset(Symbol("3D4"), :phi, #= none:3 =# @perm_str("(1,2,4)(5,6,7)(8,10,9)(13,14,16)(17,18,19)(20,22,21)"))
chevieset(Symbol("3D4"), :CharName, function (arg...,)
        return PartitionTupleToString(arg[1])
    end)
chevieset(Symbol("3D4"), :CartanMat, [[2, 0, -1, 0], [0, 2, -1, 0], [-1, -1, 2, -1], [0, 0, -1, 2]])
chevieset(Symbol("3D4"), :ClassParameter, function (w,)
        local x
        if w == []
            x = Perm()
        else
            x = Product((chevieget(Symbol("3D4"), :generators))[w])
        end
        return (chevieget(Symbol("3D4"), :classparams))[Position(chevieget(Symbol("3D4"), :cyclestructure), CycleStructurePerm(x * chevieget(Symbol("3D4"), :phi)))]
    end)
chevieset(Symbol("3D4"), :ClassInfo, function ()
        local res
        res = Dict{Symbol, Any}(:classtext => [[1], [], [1, 2, 3, 1, 2, 3], [3], [1, 3], [1, 2, 3, 1, 2, 4, 3, 2], [1, 2, 3, 2]], :classnames => ["C_3", "\\tilde A_2", "C_3+A_1", "\\tilde A_2+A_1", "F_4", "\\tilde A_2+A_2", "F_4(a_1)"], :orders => [6, 3, 6, 6, 12, 3, 6], :classes => [48, 16, 16, 48, 48, 8, 8])
        res[:classparams] = res[:classnames]
        return res
    end)
chevieset(Symbol("3D4"), :NrConjugacyClasses, 7)
chevieset(Symbol("3D4"), :CharInfo, function ()
        return Dict{Symbol, Any}(:extRefl => [1, 5, 4, 6, 2], :charparams => [[[], [4]], [[], [1, 1, 1, 1]], [[], [2, 2]], [[1, 1], [2]], [[1], [3]], [[1], [1, 1, 1]], [[1], [2, 1]]], :charRestrictions => [13, 4, 10, 5, 11, 3, 6], :nrGroupClasses => 13, :b => [0, 12, 4, 4, 1, 7, 3], :B => [0, 12, 8, 8, 5, 11, 9])
    end)
chevieset(Symbol("3D4"), :HeckeCharTable, function (param, sqrtparam)
        local q, tbl
        q = -((param[1])[1]) // (param[1])[2]
        tbl = Dict{Symbol, Any}(:identifier => "H(^3D_4)", :parameter => [q, q, q, q], :sqrtparameter => [], :cartan => chevieget(Symbol("3D4"), :CartanMat), :size => 192, :irreducibles => [[q, 1, q ^ 6, q, q ^ 2, q ^ 8, q ^ 4], [-1, 1, 1, -1, 1, 1, 1], [q - 1, 2, 2 * q ^ 3, q - 1, -q, -(q ^ 4), -(q ^ 2)], [0, 0, (q ^ 4 - 2 * q ^ 3) + q ^ 2, 0, -q, 3 * q ^ 4, 3 * q ^ 2], [q, 1, q ^ 5 - 2 * q ^ 4, -1, 0, -2 * q ^ 6, 2 * q ^ 3], [-1, 1, -2 * q ^ 2 + q, q, 0, -2 * q ^ 2, 2q], [q - 1, 2, -(q ^ 4) - q ^ 2, q - 1, 0, 2 * q ^ 4, -2 * q ^ 2]] * q ^ 0, :irredinfo => chevieget(Symbol("3D4"), :IrredInfo))
        Inherit(tbl, (chevieget(Symbol("3D4"), :ClassInfo))())
        tbl[:centralizers] = map((x->begin
                        tbl[:size] // x
                    end), tbl[:classes])
        tbl = ((CHEVIE[:compat])[:MakeCharacterTable])(tbl)
        ((CHEVIE[:compat])[:AdjustHeckeCharTable])(tbl, param)
        return tbl
    end)
chevieset(Symbol("3D4"), :PhiFactors, [1, E(3), 1, E(3, 2)])
chevieset(Symbol("3D4"), :Representation, function (i,)
        return (chevieget(Symbol("3D4"), :HeckeRepresentation))([[1, -1], [1, -1], [1, -1], [1, -1]], [1, 1, 1, 1], i)
    end)
chevieset(Symbol("3D4"), :HeckeRepresentation, function (param, sqrtparam, i)
        local q, v, res, x
        q = -((param[1])[1]) // (param[1])[2]
        if !(sqrtparam[1] !== nothing)
            v = GetRoot(q, 2, "Representation(Hecke(^3D_4)[", i, "])")
        else
            v = -(sqrtparam[1]) // (param[1])[2]
        end
        res = [Dict{Symbol, Any}(:gens => [[[v ^ 2]], [[v ^ 2]], [[v ^ 2]], [[v ^ 2]]], :F => [[1]]), Dict{Symbol, Any}(:gens => [[[-1]], [[-1]], [[-1]], [[-1]]], :F => [[1]]), Dict{Symbol, Any}(:gens => [[[-1, v], [0, v ^ 2]], [[-1, v], [0, v ^ 2]], [[v ^ 2, 0], [v, -1]], [[-1, v], [0, v ^ 2]]], :F => [[1, 0], [0, 1]]), Dict{Symbol, Any}(:gens => [[[v ^ 2, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, v ^ 2, 0, 0, 0], [0, 0, 0, v ^ 2, 0, 0], [0, 0, v, 0, -1, 0], [0, 0, 0, v, 0, -1]], [[-1, 0, 0, v, 0, 0], [0, v ^ 2, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 0, 0, v ^ 2, 0, 0], [0, v, 0, 0, -1, 0], [0, 0, 0, 0, 0, v ^ 2]], [[v ^ 2, 0, 0, 0, 0, 0], [0, -1, 0, 0, v, v], [-v, 0, -1, 0, v, 0], [v, 0, 0, -1, 0, v], [0, 0, 0, 0, v ^ 2, 0], [0, 0, 0, 0, 0, v ^ 2]], [[-1, 0, -v, 0, 0, 0], [0, v ^ 2, 0, 0, 0, 0], [0, 0, v ^ 2, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, v ^ 2, 0], [0, v, 0, 0, 0, -1]]], :F => [[0, 0, 0, 0, 0, -1], [0, 0, -1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, -1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0]]), Dict{Symbol, Any}(:gens => [[[-1, v, 0, 0], [0, v ^ 2, 0, 0], [0, 0, v ^ 2, 0], [0, 0, 0, v ^ 2]], [[v ^ 2, 0, 0, 0], [0, v ^ 2, 0, 0], [0, v, -1, 0], [0, 0, 0, v ^ 2]], [[v ^ 2, 0, 0, 0], [v, -1, v, v], [0, 0, v ^ 2, 0], [0, 0, 0, v ^ 2]], [[v ^ 2, 0, 0, 0], [0, v ^ 2, 0, 0], [0, 0, v ^ 2, 0], [0, v, 0, -1]]], :F => [[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]), Dict{Symbol, Any}(:gens => [[[-1, 0, 0, 0], [0, v ^ 2, 0, 0], [0, v, -1, 0], [0, 0, 0, -1]], [[v ^ 2, 0, 0, 0], [0, -1, 0, 0], [v, 0, -1, 0], [0, 0, 0, -1]], [[-1, 0, v, 0], [0, -1, v, 0], [0, 0, v ^ 2, 0], [0, 0, v, -1]], [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, v], [0, 0, 0, v ^ 2]]], :F => [[0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0]]), Dict{Symbol, Any}(:gens => [[[-1, v, 0, 0, 0, 0, 0, 2v], [0, v ^ 2, 0, 0, 0, 0, 0, 0], [0, v, -1, 0, 0, 0, 0, 0], [0, 0, 0, v ^ 2, 0, 0, 0, 0], [0, 0, 0, 0, v ^ 2, 0, 0, 0], [0, 0, 0, v, 0, -1, 0, 0], [0, 0, 0, 0, v, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, v ^ 2]], [[-1, v, 0, 0, 0, 0, 2v, 0], [0, v ^ 2, 0, 0, 0, 0, 0, 0], [0, 0, v ^ 2, 0, 0, 0, 0, 0], [0, v, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, v ^ 2, 0, 0, 0], [0, 0, v, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, v ^ 2, 0], [0, 0, 0, 0, v, 0, 0, -1]], [[v ^ 2, 0, 0, 0, 0, 0, 0, 0], [v, -1, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, v, v, 0], [0, 0, 0, -1, 0, v, 0, v], [0, 0, 0, 0, -1, 0, v, v], [0, 0, 0, 0, 0, v ^ 2, 0, 0], [0, 0, 0, 0, 0, 0, v ^ 2, 0], [0, 0, 0, 0, 0, 0, 0, v ^ 2]], [[-1, v, 0, 0, 0, 2v, 0, 0], [0, v ^ 2, 0, 0, 0, 0, 0, 0], [0, 0, v ^ 2, 0, 0, 0, 0, 0], [0, 0, 0, v ^ 2, 0, 0, 0, 0], [0, v, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, v ^ 2, 0, 0], [0, 0, v, 0, 0, 0, -1, 0], [0, 0, 0, v, 0, 0, 0, -1]]], :F => [[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0]])]
        res = res[i]
        res[:gens] = res[:gens] * v ^ 0
        return res
    end)
chevieset(Symbol("3D4"), :UnipotentCharacters, function ()
        local res
        return Dict{Symbol, Any}(:harishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "G", :indices => [1, 3], :rank => 2), :levi => [], :eigenvalue => 1, :parameterExponents => [3, 1], :cuspidalName => "", :charNumbers => [1, 2, 6, 5, 7, 3]), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:4, :eigenvalue => 1, :parameterExponents => [], :cuspidalName => "{}^3D_4[1]", :charNumbers => [4]), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:4, :eigenvalue => -1, :parameterExponents => [], :cuspidalName => "{}^3D_4[-1]", :charNumbers => [8])], :almostHarishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:orbit => [Dict{Symbol, Any}(:series => "D", :indices => 1:4, :rank => 4)], :twist => #= none:21 =# @perm_str("(1,2,4)")), :levi => [], :eigenvalue => 1, :cuspidalName => "", :charNumbers => 1:7), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:orbit => [Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0)], :twist => Perm()), :levi => 1:4, :eigenvalue => -1, :cuspidalName => "D_4", :charNumbers => [8])], :families => [Family("C1", [1]), Family("C1", [2]), Family("C1", [5]), Family("C1", [6]), Family("C2", [7, 4, 3, 8])], :a => [0, 12, 3, 3, 1, 7, 3, 3], :A => [0, 12, 9, 9, 5, 11, 9, 9])
    end)
