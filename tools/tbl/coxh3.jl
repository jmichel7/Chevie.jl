
chevieset(:H3, :ReflectionDegrees, [2, 6, 10])
chevieset(:H3, :Size, 120)
chevieset(:H3, :GeneratingRoots, function ()
        local a
        a = (1 + ER(5)) // 2
        return [[a, -1, a - 1], [-a, 1, a - 1], [1, a - 1, -a]] // 2
    end)
chevieset(:H3, :NrConjugacyClasses, 10)
chevieset(:H3, :cyclestructure, [])
chevieset(:H3, :generators, [])
chevieset(:H3, :CartanMat, [[2, -((1 + ER(5))) // 2, 0], [-((1 + ER(5))) // 2, 2, -1], [0, -1, 2]])
chevieset(:H3, :PowerMaps, [nothing, [1, 1, 7, 1, 5, 3, 3, 5, 7, 1], [1, 2, 7, 4, 1, 9, 3, 10, 6, 10], nothing, [1, 2, 1, 4, 5, 10, 1, 8, 10, 10], nothing, [1, 2, 7, 4, 5, 9, 3, 8, 6, 10]])
chevieset(:H3, :WordsClassRepresentatives, [[], [1], [1, 2], [1, 3], [2, 3], [1, 2, 3], [1, 2, 1, 2], [1, 2, 1, 2, 3], [1, 2, 1, 2, 3, 2, 1, 2, 3], [1, 2, 1, 2, 1, 3, 2, 1, 2, 1, 3, 2, 1, 2, 3]])
chevieset(:H3, :ParabolicRepresentatives, function (s,)
        local t
        t = [[[]], [[1]], [[1, 2], [1, 3], [2, 3]], [1:3]]
        return t[s + 1]
    end)
chevieset(:H3, :ClassInfo, function ()
        local res
        res = Dict{Symbol, Any}(:classtext => chevieget(:H3, :WordsClassRepresentatives), :orders => [1, 2, 5, 2, 3, 10, 5, 6, 10, 2], :classes => [1, 15, 12, 15, 20, 12, 12, 20, 12, 1])
        res[:classnames] = map(joindigits, res[:classtext])
        (res[:classnames])[1] = "."
        res[:classparams] = res[:classnames]
        return res
    end)
chevieset(:H3, :CharInfo, function ()
        local res
        res = Dict{Symbol, Any}(:charparams => [[1, 15], [1, 0], [5, 5], [5, 2], [3, 6], [3, 8], [3, 1], [3, 3], [4, 3], [4, 4]], :gp => ["1_r'", "1_r", "5_r'", "5_r", "3_s", "overline{3}_s", "3_s'", "overline{3}_s'", "4_r'", "4_r"], :opdam => perm"(9,10)", :extRefl => [2, 7, 5, 1])
        res[:b] = map((x->begin
                        x[2]
                    end), res[:charparams])
        return res
    end)
chevieset(:H3, :vpolheckeirreducibles, [[[[1], 0], [[-1], 0], [[1], 0], [[1], 0], [[1], 0], [[-1], 0], [[1], 0], [[-1], 0], [[-1], 0], [[-1], 0]], [[[1], 0], [[1], 1], [[1], 2], [[1], 2], [[1], 2], [[1], 3], [[1], 4], [[1], 5], [[1], 9], [[1], 15]], [[[5], 0], [[-3, 2], 0], [[1, -1], 0], [[2, -2, 1], 0], [[1, -2], 0], [[], 0], [[1, 0, -1], 0], [[1], 2], [[], 0], [[-5], 6]], [[[5], 0], [[-2, 3], 0], [[-1, 1], 1], [[1, -2, 2], 0], [[-2, 1], 1], [[], 0], [[-1, 0, 1], 2], [[-1], 3], [[], 0], [[5], 9]], [[[3], 0], [[-2, 1], 0], [[1, (-1 + ER(5)) // 2], 0], [[1, -2], 0], [[1, -1], 0], [[(1 - ER(5)) // 2], 1], [[1, 0, (-1 - ER(5)) // 2], 0], [[], 0], [[(1 + ER(5)) // 2], 3], [[3], 5]], [[[3], 0], [[-2, 1], 0], [[1, (-1 - ER(5)) // 2], 0], [[1, -2], 0], [[1, -1], 0], [[(1 + ER(5)) // 2], 1], [[1, 0, (-1 + ER(5)) // 2], 0], [[], 0], [[(1 - ER(5)) // 2], 3], [[3], 5]], [[[3], 0], [[-1, 2], 0], [[(-1 + ER(5)) // 2, 1], 1], [[-2, 1], 1], [[-1, 1], 1], [[(-1 + ER(5)) // 2], 2], [[(-1 - ER(5)) // 2, 0, 1], 2], [[], 0], [[(-1 - ER(5)) // 2], 6], [[-3], 10]], [[[3], 0], [[-1, 2], 0], [[(-1 - ER(5)) // 2, 1], 1], [[-2, 1], 1], [[-1, 1], 1], [[(-1 - ER(5)) // 2], 2], [[(-1 + ER(5)) // 2, 0, 1], 2], [[], 0], [[(-1 + ER(5)) // 2], 6], [[-3], 10]], [[[4], 0], [[-2, 2], 0], [[-1], 1], [[1, -2, 1], 0], [[1, -1, 1], 0], [[1], 3 // 2], [[-1], 2], [[-1], 5 // 2], [[1], 9 // 2], [[-4], 15 // 2]], [[[4], 0], [[-2, 2], 0], [[-1], 1], [[1, -2, 1], 0], [[1, -1, 1], 0], [[-1], 3 // 2], [[-1], 2], [[1], 5 // 2], [[-1], 9 // 2], [[4], 15 // 2]]])
chevieset(:H3, :CycPolSchurElements, [[1, -15, 2, 2, 2, 3, 5, 6, 10], [1, 0, 2, 2, 2, 3, 5, 6, 10], [1, -5, 2, 2, 2, 3, 6], [1, -2, 2, 2, 2, 3, 6], [(5 + ER(5)) // 2, -6, 2, 2, 2, 2 // 5, 3 // 5, 1 // 10, 9 // 10], [(5 - ER(5)) // 2, -6, 2, 2, 2, 1 // 5, 4 // 5, 3 // 10, 7 // 10], [(5 + ER(5)) // 2, -1, 2, 2, 2, 2 // 5, 3 // 5, 1 // 10, 9 // 10], [(5 - ER(5)) // 2, -1, 2, 2, 2, 1 // 5, 4 // 5, 3 // 10, 7 // 10], [2, -3, 3, 5], [2, -3, 3, 5]])
chevieset(:H3, :sparseFakeDegrees, [[1, 15], [1, 0], [1, 5, 1, 7, 1, 9, 1, 11, 1, 13], [1, 2, 1, 4, 1, 6, 1, 8, 1, 10], [1, 6, 1, 10, 1, 14], [1, 8, 1, 10, 1, 12], [1, 1, 1, 5, 1, 9], [1, 3, 1, 5, 1, 7], [1, 3, 1, 7, 1, 9, 1, 11], [1, 4, 1, 6, 1, 8, 1, 12]])
chevieset(:H3, :HeckeCharTable, function (param, sqrtparam)
        local a, q, v, ci, tbl
        a = (1 + ER(5)) // 2
        q = -((param[1])[1]) // (param[1])[2]
        if !(sqrtparam[1] !== nothing)
            v = GetRoot(q, 2, "CharTable(Hecke(H3))")
        else
            v = -(sqrtparam[1]) // (param[1])[2]
        end
        ci = (chevieget(:H3, :ClassInfo))()
        tbl = Dict{Symbol, Any}(:identifier => "H(H3)", :text => "the representing matrices are those of Lusztig(1981)", :parameter => [q, q, q], :cartan => chevieget(:H3, :CartanMat), :size => 120, :order => 120, :powermap => chevieget(:H3, :PowerMaps), :irreducibles => map((i->begin
                                map(function (j,)
                                        local res
                                        res = ValuePol(j[1], q)
                                        if IsInt(j[2])
                                            res = res * q ^ j[2]
                                        else
                                            res = res * v ^ (2 * j[2])
                                        end
                                        return res
                                    end, i)
                            end), chevieget(:H3, :vpolheckeirreducibles)), :irredinfo => chevieget(:H3, :IrredInfo))
        Inherit(tbl, ci)
        tbl[:centralizers] = map((x->begin
                        tbl[:size] // x
                    end), tbl[:classes])
        tbl = ((CHEVIE[:compat])[:MakeCharacterTable])(tbl)
        ((CHEVIE[:compat])[:AdjustHeckeCharTable])(tbl, param)
        return tbl
    end)
chevieset(:H3, :Representation, function (i,)
        return (chevieget(:H3, :HeckeRepresentation))([[1, -1], [1, -1], [1, -1]], [1, 1, 1], i)
    end)
chevieset(:H3, :WGraphs, [[[[1, 2, 3]], []], 1, [[[2], [1, 2], [1, 3], [1, 3], [2, 3]], [[-1, [[1, 3], [2, 4], [3, 5], [4, 5]]]]], 3, [[[1, 2], [1, 3], [2, 3]], [[-1, [[1, 2]]], [(-1 - ER(5)) // 2, [[2, 3]]]]], [[[1, 2], [1, 3], [2, 3]], [[-1, [[1, 2]]], [(-1 + ER(5)) // 2, [[2, 3]]]]], 5, 6, [[[1], [2], [1, 3], [2, 3]], [[1, [[1, 2, 3], [2, 3, 4], [3, 4]]]]], 9])
chevieset(:H3, :WGraph, function (i,)
        local gr
        gr = chevieget(:H3, :WGraphs)
        if IsInt(gr[i])
            return DualWGraph(3, gr[gr[i]])
        else
            return gr[i]
        end
    end)
chevieset(:H3, :HeckeRepresentation, function (param, sqrtparam, i)
        local v
        if !(sqrtparam[1] !== nothing)
            v = GetRoot(-((param[1])[1]) // (param[1])[2], 2, "Representation(Hecke(H3),[", i, "])")
        else
            v = -(sqrtparam[1]) // (param[1])[2]
        end
        return -((param[1])[2]) * WGraphToRepresentation(3, (chevieget(:H3, :WGraph))(i), v)
    end)
chevieset(:H3, :UnipotentCharacters, function ()
        local res
        res = Dict{Symbol, Any}(:harishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "H", :indices => 1:3, :rank => 3), :levi => [], :eigenvalue => 1, :parameterExponents => [1, 1, 1], :cuspidalName => "", :charNumbers => 1:10), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [3], :rank => 1), :levi => 1:2, :eigenvalue => E(5, 2), :parameterExponents => [5], :cuspidalName => "I_2(5)[1,3]", :charNumbers => [11, 13]), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [3], :rank => 1), :levi => 1:2, :eigenvalue => E(5, 3), :parameterExponents => [5], :cuspidalName => "I_2(5)[1,2]", :charNumbers => [12, 14]), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :eigenvalue => E(4), :qEigen => 1 // 2, :parameterExponents => [], :cuspidalName => "H_3[i]", :charNumbers => [15]), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :eigenvalue => -(E(4)), :qEigen => 1 // 2, :parameterExponents => [], :cuspidalName => "H_3[-i]", :charNumbers => [16])], :families => [Family("C1", [2]), Family(((CHEVIE[:families])[:Dihedral])(5), [7, 8, 14, 13], Dict{Symbol, Any}(:ennola => -1)), Family("C1", [4]), Family("C'\"2", [9, 10, 15, 16], Dict{Symbol, Any}(:ennola => 3)), Family("C1", [3], Dict{Symbol, Any}(:ennola => -1)), Family(((CHEVIE[:families])[:Dihedral])(5), [5, 6, 12, 11], Dict{Symbol, Any}(:ennola => 1)), Family("C1", [1], Dict{Symbol, Any}(:ennola => -1))], :a => [15, 0, 5, 2, 6, 6, 1, 1, 3, 3, 6, 6, 1, 1, 3, 3], :A => [15, 0, 13, 10, 14, 14, 9, 9, 12, 12, 14, 14, 9, 9, 12, 12])
        return res
    end)
chevieset(:H3, :Invariants, function ()
        local r, C
        C = chevieget(:H3, :CartanMat)
        r = roots(C) * C
        return map((d->begin
                        function (arg...,)
                            return Sum(r, (a->begin
                                            (arg * a) ^ d
                                        end))
                        end
                    end), chevieget(:H3, :ReflectionDegrees))
    end)
chevieset(:H3, :Discriminant, function ()
        return function (a, b, c)
                return (((((((131835937500 * a * b ^ 3 * c - 100195312500 * a ^ 2 * b * c ^ 2) + 395507812500 * c ^ 3) - 28369140625 * a ^ 3 * b ^ 4) + 1371093750 * a ^ 4 * b ^ 2 * c + 175781250000 * b ^ 5 + 1191796875 * a ^ 5 * c ^ 2 + 1162187500 * a ^ 6 * b ^ 3) - 74250000 * a ^ 7 * b * c) - 22233750 * a ^ 9 * b ^ 2) + 438750 * a ^ 10 * c + 213700 * a ^ 12 * b) - 829 * a ^ 15
            end
    end)
chevieset(:H3, :KLeftCellRepresentatives, [Dict{Symbol, Any}(:character => [2], :duflo => [1, 2, 3], :reps => ""), Dict{Symbol, Any}(:character => [1], :duflo => [16, 17, 18], :reps => ""), Dict{Symbol, Any}(:character => [3], :duflo => [1, 24, 3], :reps => ""), Dict{Symbol, Any}(:character => [4], :duflo => [2, 1, 28], :reps => ""), Dict{Symbol, Any}(:character => [6, 5], :duflo => [1, 20, 18], :reps => [[7, 19, 24]]), Dict{Symbol, Any}(:character => [8, 7], :duflo => [1, 6, 18], :reps => [[9, 2, 27]]), Dict{Symbol, Any}(:character => [10, 9], :duflo => [8, 18, 17], :reps => [[11, 17, 25], [11, 27, 10], [14, 30, 4]]), Dict{Symbol, Any}(:character => [10, 9], :duflo => [13, 30, 8], :reps => [[10, 29, 5], [12, 21, 22], [13, 22, 23]])])