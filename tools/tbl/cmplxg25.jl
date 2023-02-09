chevieset(:G25, :GeneratingRoots, [[0, 0, -1], -((2 * E(3, 2) + 1)) // 3 * [1, 1, 1], [0, 1, 0]])
chevieset(:G25, :EigenvaluesGeneratingReflections, [1 // 3, 1 // 3, 1 // 3])
chevieset(:G25, :HyperplaneRepresentatives, [1])
chevieset(:G25, :BraidRelations, [[[1, 2, 1], [2, 1, 2]], [[1, 3], [3, 1]], [[2, 3, 2], [3, 2, 3]]])
chevieset(:G25, :Size, 648)
chevieset(:G25, :ReflectionDegrees, [6, 9, 12])
chevieset(:G25, :NrConjugacyClasses, 24)
chevieset(:G25, :ParabolicRepresentatives, function (s,)
        local t
        t = [[[]], [[1]], [[1, 2], [1, 3]], [1:3]]
        return t[s + 1]
    end)
chevieset(:G25, :ClassNames, [".", "cc", "31", "3131", "12231223", "1223", "d", "dd", "z", "zz", "2231223", "d1", "1", "131", "3221223221", "11", "1122", "12", "12z", "122312231223", "332112", "212", "c", "cz"])
chevieset(:G25, :WordsClassRepresentatives, map((x->begin
                StringToDigits(Replace(x, ".", "", "z", "cccc", "c", "123", "d", "1232"))
            end), chevieget(:G25, :ClassNames)))
chevieset(:G25, :ClassInfo, Dict{Symbol, Any}(:classtext => chevieget(:G25, :WordsClassRepresentatives), :classnames => chevieget(:G25, :ClassNames), :classparams => chevieget(:G25, :ClassNames), :orders => [1, 6, 3, 3, 3, 6, 9, 9, 3, 3, 6, 6, 3, 3, 3, 3, 6, 6, 6, 2, 6, 4, 12, 12], :classes => [1, 9, 12, 12, 12, 36, 72, 72, 1, 1, 36, 36, 12, 24, 12, 12, 36, 36, 36, 9, 9, 54, 54, 54]))
chevieset(:G25, :PowerMaps, [nothing, [1, 9, 4, 3, 15, 5, 8, 7, 10, 9, 3, 15, 16, 14, 5, 13, 16, 13, 4, 1, 10, 20, 2, 21], [1, 20, 1, 1, 1, 20, 9, 10, 1, 1, 20, 20, 1, 1, 1, 1, 20, 20, 20, 20, 20, 22, 22, 22], nothing, [1, 21, 4, 3, 15, 12, 8, 7, 10, 9, 19, 6, 16, 14, 5, 13, 18, 17, 11, 20, 2, 22, 24, 23], nothing, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24], nothing, nothing, nothing, [1, 21, 4, 3, 15, 12, 8, 7, 10, 9, 19, 6, 16, 14, 5, 13, 18, 17, 11, 20, 2, 22, 24, 23]])
chevieset(:G25, :CharInfo, function ()
        local res
        res = Dict{Symbol, Any}(:charparams => [[1, 0], [1, 24], [1, 12], [2, 15], [2, 3], [2, 9], [3, 6], [3, 5, 2], [3, 5, 1], [3, 17], [3, 13, 2], [3, 1], [3, 13, 1], [6, 8, 2], [6, 8, 1], [6, 2], [6, 4, 2], [6, 10], [6, 4, 1], [8, 3], [8, 9], [8, 6], [9, 5], [9, 7]], :extRefl => [1, 12, 8, 3])
        res[:b] = map((x->begin
                        x[2]
                    end), res[:charparams])
        res[:charnames] = map(exceptioCharName, res[:charparams])
        return res
    end)
chevieset(:G25, :HeckeCharTable, function (para, root)
        local u, v, w, f10, f23, f31, f62, f83, f97, res, c
        u = (para[1])[1]
        v = (para[1])[2]
        w = (para[1])[3]
        c = (u * v * w) ^ 0
        res = Dict{Symbol, Any}(:name => "H(G25)", :identifier => "H(G25)", :parameter => para, :size => 648, :order => 648, :dim => 3, :degrees => [6, 9, 12], :reflclasses => [13], :powermap => chevieget(:G25, :PowerMaps), :irredinfo => chevieget(:G25, :IrredInfo))
        f10 = (y->begin
                    map((w->begin
                                y ^ length(w)
                            end), res[:classtext])
                end)
        f23(u, v, w) = begin
                return [2, -2 * (u * v) ^ 3, u ^ 2 + v ^ 2, u ^ 4 + v ^ 4, (u * v) ^ 2 * (u ^ 4 + v ^ 4), -u * v * (u ^ 2 + v ^ 2), -(u ^ 2) * v ^ 2, -(u ^ 4) * v ^ 4, 2 * u ^ 6 * v ^ 6, 2 * u ^ 12 * v ^ 12, -(v ^ 3) * u ^ 3 * (u + v), -(v ^ 2) * u ^ 2 * (u + v), u + v, (u + v) * ((u ^ 2 - u * v) + v ^ 2), u ^ 4 * v ^ 4 * (u ^ 2 + v ^ 2), u ^ 2 + v ^ 2, -u * v * (u ^ 2 + v ^ 2), u * v, u ^ 7 * v ^ 7, -(u ^ 3) * v ^ 3 * (u ^ 2 + v ^ 2) * ((v ^ 4 - u ^ 2 * v ^ 2) + u ^ 4), -2 * u ^ 3 * v ^ 3, 0, 0, 0]
            end
        f31(u, v, w) = begin
                return [3, -(u ^ 4) * v ^ 2, 2 * u * v + u ^ 2, 2 * u ^ 2 * v ^ 2 + u ^ 4, u ^ 4 * v ^ 4 + 2 * u ^ 6 * v ^ 2, -(u ^ 2) * v ^ 2, 0, 0, 3 * u ^ 8 * v ^ 4, 3 * u ^ 16 * v ^ 8, -(u ^ 4) * v ^ 3, -(u ^ 4) * v, 2u + v, u * v ^ 2 + u ^ 2 * v + u ^ 3, 2 * u ^ 6 * v ^ 4 + u ^ 8 * v ^ 2, 2 * u ^ 2 + v ^ 2, (-u * v ^ 3 - u ^ 3 * v) + u ^ 4, u * v + u ^ 2, u ^ 9 * v ^ 5 + u ^ 10 * v ^ 4, -(u ^ 6) * v ^ 6, u ^ 2 * v ^ 4 - 2 * u ^ 5 * v, u ^ 3, u ^ 2 * v, u ^ 10 * v ^ 5]
            end
        f62(u, v, w) = begin
                return [6, 2 * u ^ 3 * v ^ 2 * w, 2 * u * v + 2 * u * w + u ^ 2 + v ^ 2, 2 * u ^ 2 * v ^ 2 + 2 * u ^ 2 * w ^ 2 + u ^ 4 + v ^ 4, u ^ 2 * (v ^ 4 * w ^ 2 + 2 * u ^ 2 * v ^ 2 * w ^ 2 + 2 * v ^ 4 * u ^ 2 + u ^ 4 * w ^ 2), u * w * (u ^ 2 + v ^ 2), 0, 0, 6 * u ^ 6 * v ^ 4 * w ^ 2, 6 * u ^ 12 * v ^ 8 * w ^ 4, u ^ 3 * v ^ 2 * w * (w + u), u ^ 2 * v ^ 2 * (w + u), 3u + 2v + w, v ^ 2 * u + u * w ^ 2 + v * u ^ 2 + u ^ 2 * w + u ^ 3 + v ^ 3, v ^ 2 * u ^ 4 * (3 * v ^ 2 * w ^ 2 + 2 * u ^ 2 * w ^ 2 + u ^ 2 * v ^ 2), 3 * u ^ 2 + 2 * v ^ 2 + w ^ 2, -((u ^ 2 + v ^ 2)) * ((u * v - w ^ 2) - u ^ 2), u * (u + v), v ^ 4 * u ^ 7 * w ^ 2 * (u + v), u ^ 3 * w ^ 3 * (u ^ 2 + v ^ 2) * ((v ^ 4 - u ^ 2 * v ^ 2) + u ^ 4), v * (-2 * u ^ 3 * w ^ 2 + 3 * u ^ 4 * v + v ^ 3 * w ^ 2), -u * (-(u ^ 2) + v * w), 0, 0]
            end
        f83(u, v, w) = begin
                return [8, 0, 2 * (w + u) * (u + v), 2 * (w ^ 2 + u ^ 2) * (u ^ 2 + v ^ 2), 2 * u ^ 2 * v * w * (w ^ 2 + u ^ 2) * (u ^ 2 + v ^ 2), 0, -(u ^ 2) * v * w, -(u ^ 4) * v ^ 2 * w ^ 2, 8 * u ^ 6 * v ^ 3 * w ^ 3, 8 * u ^ 12 * v ^ 6 * w ^ 6, 0, 0, 4u + 2v + 2w, v ^ 2 * u + u * w ^ 2 + v * w ^ 2 + v * u ^ 2 + u ^ 2 * w + v ^ 2 * w + 2 * u ^ 3, 2 * u ^ 4 * v * w * (2 * v ^ 2 * w ^ 2 + u ^ 2 * w ^ 2 + u ^ 2 * v ^ 2), 4 * u ^ 2 + 2 * v ^ 2 + 2 * w ^ 2, ((((-u * v ^ 3 - u * w ^ 3) + u ^ 2 * v ^ 2 + u ^ 2 * w ^ 2 + v ^ 2 * w ^ 2) - u ^ 3 * v) - u ^ 3 * w) + u ^ 4, u * (u + v + w), v ^ 3 * u ^ 7 * w ^ 3 * (u + v + w), 0, (-2 * u ^ 3 * v ^ 3 - 2 * u ^ 3 * w ^ 3) + v ^ 3 * w ^ 3 + 3 * u ^ 4 * v * w, -u * (-(u ^ 2) + v * w), 0, 0]
            end
        f97(u, v, w, J) = begin
                return [9, -3 * J ^ 2 * u ^ 2 * v ^ 2 * w ^ 2, (u + v + w) ^ 2, (u ^ 2 + w ^ 2 + v ^ 2) ^ 2, J * (u ^ 2 * v ^ 2 + u ^ 2 * w ^ 2 + v ^ 2 * w ^ 2) ^ 2, -(J ^ 2) * (u ^ 2 * v ^ 2 + u ^ 2 * w ^ 2 + v ^ 2 * w ^ 2), 0, 0, 9 * J * u ^ 4 * v ^ 4 * w ^ 4, 9 * J ^ 2 * u ^ 8 * v ^ 8 * w ^ 8, -(v ^ 2) * J ^ 2 * u ^ 2 * w ^ 2 * (u + v + w), -v * u * w * J ^ 2 * (u * v + u * w + v * w), 3u + 3v + 3w, (u + v + w) * (u ^ 2 + w ^ 2 + v ^ 2), 3 * J * u ^ 2 * v ^ 2 * w ^ 2 * (u ^ 2 * v ^ 2 + u ^ 2 * w ^ 2 + v ^ 2 * w ^ 2), 3 * u ^ 2 + 3 * v ^ 2 + 3 * w ^ 2, (((((-u * v ^ 3 - u * w ^ 3) - v * w ^ 3) + u ^ 2 * v ^ 2 + u ^ 2 * w ^ 2 + v ^ 2 * w ^ 2) - u ^ 3 * v) - u ^ 3 * w) - v ^ 3 * w, u * v + u * w + v * w, v ^ 4 * J * u ^ 4 * w ^ 4 * (u * v + u * w + v * w), (-(u ^ 6) * v ^ 6 - u ^ 6 * w ^ 6) - v ^ 6 * w ^ 6, -J * u * v * w * (((2 * w ^ 3 + 2 * v ^ 3) - 3 * u * v * w) + 2 * u ^ 3), -u * v * w, -J * u * v * w, -(J ^ 2) * u ^ 5 * v ^ 5 * w ^ 5]
            end
        Inherit(res, chevieget(:G25, :ClassInfo))
        res[:centralizers] = map((x->begin
                        div(res[:order], x)
                    end), res[:classes])
        res[:irreducibles] = [f10(u), f10(w), f10(v), f23(v, w, u), f23(u, v, w), f23(u, w, v), [3, 3 * u ^ 2 * v ^ 2 * w ^ 2, u ^ 2 + v ^ 2 + w ^ 2, u ^ 4 + v ^ 4 + w ^ 4, u ^ 4 * v ^ 4 + u ^ 4 * w ^ 4 + v ^ 4 * w ^ 4, u ^ 2 * v ^ 2 + u ^ 2 * w ^ 2 + v ^ 2 * w ^ 2, 0, 0, 3 * u ^ 4 * v ^ 4 * w ^ 4, 3 * u ^ 8 * v ^ 8 * w ^ 8, u ^ 2 * v ^ 2 * w ^ 3 + u ^ 2 * v ^ 3 * w ^ 2 + u ^ 3 * v ^ 2 * w ^ 2, u * v ^ 2 * w ^ 2 + u ^ 2 * v * w ^ 2 + u ^ 2 * v ^ 2 * w, u + v + w, u ^ 3 + v ^ 3 + w ^ 3, u ^ 2 * v ^ 4 * w ^ 4 + u ^ 4 * v ^ 2 * w ^ 4 + u ^ 4 * v ^ 4 * w ^ 2, u ^ 2 + v ^ 2 + w ^ 2, u ^ 2 * v ^ 2 + u ^ 2 * w ^ 2 + v ^ 2 * w ^ 2, 0, 0, u ^ 6 * v ^ 6 + u ^ 6 * w ^ 6 + v ^ 6 * w ^ 6, 3 * u ^ 2 * v ^ 2 * w ^ 2, -u * v * w, -u * v * w, -(u ^ 5) * v ^ 5 * w ^ 5], f31(v, u, w), f31(u, w, v), f31(w, v, u), f31(w, u, v), f31(u, v, w), f31(v, w, u), f62(w, u, v), f62(v, w, u), f62(u, v, w), f62(v, u, w), f62(w, v, u), f62(u, w, v), f83(u, v, w), f83(w, v, u), f83(v, u, w), f97(u, v, w, E(3, 2)), f97(u, v, w, E(3))] * c
        res = ((CHEVIE[:compat])[:MakeCharacterTable])(res)
        return res
    end)
chevieset(:G25, :CharTable, function ()
        return (chevieget(:G25, :HeckeCharTable))([[1, E(3), E(3, 2)]], [])
    end)
chevieset(:G25, :sparseFakeDegrees, [[1, 0], [1, 24], [1, 12], [1, 15, 1, 21], [1, 3, 1, 9], [1, 9, 1, 15], [1, 6, 1, 12, 1, 18], [1, 5, 1, 8, 1, 11], [1, 5, 1, 8, 1, 11], [1, 17, 1, 20, 1, 23], [1, 13, 1, 16, 1, 19], [1, 1, 1, 4, 1, 7], [1, 13, 1, 16, 1, 19], [1, 8, 1, 11, 2, 14, 1, 17, 1, 20], [1, 8, 1, 11, 2, 14, 1, 17, 1, 20], [1, 2, 1, 5, 2, 8, 1, 11, 1, 14], [1, 4, 1, 7, 2, 10, 1, 13, 1, 16], [1, 10, 1, 13, 2, 16, 1, 19, 1, 22], [1, 4, 1, 7, 2, 10, 1, 13, 1, 16], [1, 3, 2, 6, 2, 9, 2, 12, 1, 15], [1, 9, 2, 12, 2, 15, 2, 18, 1, 21], [1, 6, 2, 9, 2, 12, 2, 15, 1, 18], [1, 5, 1, 8, 3, 11, 2, 14, 2, 17], [2, 7, 2, 10, 3, 13, 1, 16, 1, 19]])
chevieset(:G25, :SchurModels, Dict{Symbol, Any}(:f1_0 => Dict{Symbol, Any}(:vcyc => [[[1, -1, 0], 1], [[1, -1, 0], 1], [[1, 0, -1], 1], [[1, 0, -1], 1], [[1, -1, 0], 4], [[1, 0, -1], 4], [[3, -2, -1], 1], [[3, -1, -2], 1], [[2, -1, -1], 3], [[2, -1, -1], 2], [[1, -1, 0], 6], [[1, 0, -1], 6]]), :f2_3 => Dict{Symbol, Any}(:vcyc => [[[1, 0, -1], 1], [[1, 0, -1], 1], [[0, 1, -1], 1], [[0, 1, -1], 1], [[1, 0, -1], 2], [[0, 1, -1], 2], [[1, -1, 0], 1], [[1, -1, 0], 1], [[1, 1, -2], 3], [[1, 1, -2], 2], [[-1, 1, 0], 6]]), :f3_1 => Dict{Symbol, Any}(:vcyc => [[[1, -1, 0], 1], [[-1, 1, 0], 1], [[1, 0, -1], 1], [[1, 0, -1], 1], [[1, 0, -1], 2], [[0, 1, -1], 1], [[1, 1, -2], 2], [[2, -1, -1], 2], [[1, 0, -1], 6], [[1, -1, 0], 4], [[2, 1, -3], 1]]), :f3_6 => Dict{Symbol, Any}(:vcyc => [[[1, -1, 0], 1], [[1, -1, 0], 1], [[1, 0, -1], 1], [[1, 0, -1], 1], [[0, 1, -1], 1], [[0, -1, 1], 1], [[-1, -1, 2], 2], [[-1, 2, -1], 2], [[-2, 1, 1], 2]]), :f6_2 => Dict{Symbol, Any}(:vcyc => [[[-1, 1, 0], 1], [[1, 0, -1], 1], [[-1, 0, 1], 1], [[0, 1, -1], 1], [[0, 1, -1], 1], [[1, 0, -1], 2], [[1, 0, -1], 6], [[1, -2, 1], 2], [[0, 1, -1], 2], [[3, -2, -1], 1]]), :f8_3 => Dict{Symbol, Any}(:vcyc => [[[0, 1, -1], 1], [[0, -1, 1], 1], [[-1, 0, 1], 1], [[-1, 1, 0], 1], [[2, -3, 1], 1], [[2, -1, -1], 3], [[2, 1, -3], 1]]), :f9_7 => Dict{Symbol, Any}(:rootUnity => E(3), :vcyc => [[[0, 0, 0, 1], 1], [[0, 0, 0, 2], 2], [[-1, 1, 0], 6], [[1, 0, -1], 6], [[0, -1, 1], 6], [[2, -1, -1, 1], 1], [[-1, 2, -1, 1], 1], [[-1, -1, 2, 1], 1]])))
chevieset(:G25, :SchurData, [Dict{Symbol, Any}(:name => "f1_0", :order => [1, 2, 3]), Dict{Symbol, Any}(:name => "f1_0", :order => [3, 2, 1]), Dict{Symbol, Any}(:name => "f1_0", :order => [2, 1, 3]), Dict{Symbol, Any}(:name => "f2_3", :order => [2, 3, 1]), Dict{Symbol, Any}(:name => "f2_3", :order => [1, 2, 3]), Dict{Symbol, Any}(:name => "f2_3", :order => [1, 3, 2]), Dict{Symbol, Any}(:name => "f3_6", :order => [1, 3, 2]), Dict{Symbol, Any}(:name => "f3_1", :order => [2, 1, 3]), Dict{Symbol, Any}(:name => "f3_1", :order => [1, 3, 2]), Dict{Symbol, Any}(:name => "f3_1", :order => [3, 2, 1]), Dict{Symbol, Any}(:name => "f3_1", :order => [3, 1, 2]), Dict{Symbol, Any}(:name => "f3_1", :order => [1, 2, 3]), Dict{Symbol, Any}(:name => "f3_1", :order => [2, 3, 1]), Dict{Symbol, Any}(:name => "f6_2", :order => [3, 1, 2]), Dict{Symbol, Any}(:name => "f6_2", :order => [2, 3, 1]), Dict{Symbol, Any}(:name => "f6_2", :order => [1, 2, 3]), Dict{Symbol, Any}(:name => "f6_2", :order => [2, 1, 3]), Dict{Symbol, Any}(:name => "f6_2", :order => [3, 2, 1]), Dict{Symbol, Any}(:name => "f6_2", :order => [1, 3, 2]), Dict{Symbol, Any}(:name => "f8_3", :order => [1, 2, 3]), Dict{Symbol, Any}(:name => "f8_3", :order => [3, 2, 1]), Dict{Symbol, Any}(:name => "f8_3", :order => [2, 1, 3]), Dict{Symbol, Any}(:name => "f9_7", :order => [1, 2, 3], :rootUnityPower => 1), Dict{Symbol, Any}(:name => "f9_7", :order => [1, 2, 3], :rootUnityPower => 2)])
chevieset(:G25, :HeckeRepresentation, function (para, root, i)
        local u, v, w, f1, f2, f31, f32, f6, f8, f9, rep
        u = (para[1])[1]
        v = (para[1])[2]
        w = (para[1])[3]
        f1 = (u->begin
                    [[[u]], [[u]], [[u]]]
                end)
        f2(v, w) = begin
                return WGraph2Representation([[[1, 3], [2]], [[1, 2, -1, v * w]]], [w, v])
            end
        f31(u, v) = begin
                return WGraph2Representation([[[1], [2], [3]], [[1, 2, u, -v], [2, 3, -v, u]]], [u, v])
            end
        f32(u, v, w) = begin
                return WGraph2Representation([[[[2], []], [[], [1, 2, 3]], [[1, 3], []]], [[1, 2, -1, u * w + v ^ 2], [1, 3, v, v], [2, 3, -u * w - v ^ 2, 1]]], [u, v, w])
            end
        f6(v, u, w) = begin
                return WGraph2Representation([[[[2], []], [[], [1, 2]], [[1], []], [[], [2, 3]], [[3], []], [[], [1, 3]]], [[1, 2, -1, v * w + u ^ 2], [1, 3, u, u], [1, 4, -1, v * w + u ^ 2], [1, 5, -u, -u], [1, 6, w, 0], [2, 3, -v * w - u ^ 2, 1], [2, 6, -u * w, 1], [4, 5, v * w + u ^ 2, -1], [4, 6, -u * w, 1]]], [v, u, w])
            end
        f8(u, w, v) = begin
                return WGraph2Representation([[[[2, 3], []], [[3], [1, 2]], [[1, 3], []], [[2], [3]], [[1, 3], []], [[2], [1]], [[1], [2, 3]], [[1, 2], []]], [[1, 2, -u * v - w ^ 2, 1], [1, 3, w, w], [1, 4, v * w - w ^ 2, 0], [1, 5, 0, -1], [2, 3, -1, u * v + w ^ 2], [2, 4, [1, 0, 3, w], -u], [2, 5, 0, -w], [2, 6, -1, 0], [3, 6, [1, 0, 3, v - w], -u], [3, 7, u * w + w ^ 2, -1], [3, 8, -w, -w], [4, 5, -u, [1, v, 3, 0]], [4, 7, 0, v], [5, 6, [1, 0, 3, 1], -u * w], [5, 7, -u, v - w], [5, 8, 0, v * w - w ^ 2], [6, 7, u * w, [1, -1, 3, 0]], [6, 8, 0, v - w], [7, 8, -1, u * v + w ^ 2]]], [u, w, v])
            end
        f9(u, v, w, a) = begin
                return WGraph2Representation([[[[2], []], [[], [1, 2, 3]], [[1], [3]], [[1, 3], []], [[2], [1]], [[1], [2]], [[2], [3]], [[3], [2]], [[3], [1]]], [[1, 2, -1, u * w + v ^ 2], [1, 3, v, [1, v, 3, 0]], [1, 4, -a * v, 0], [1, 5, 0, a ^ 2 * u - v], [1, 6, 0, a ^ 2 * u], [1, 7, 0, a ^ 2 * u - v], [1, 8, 0, -(a ^ 2) * u], [1, 9, v, [1, 0, 3, v]], [2, 3, -u * w - v ^ 2, 1], [2, 4, -u * w + a * v ^ 2, 0], [2, 5, -(a ^ 2) * v * w, 0], [2, 7, -(a ^ 2) * v * w, 0], [2, 9, -u * w - v ^ 2, 1], [3, 4, 0, u + a ^ 2 * v], [3, 5, 0, u], [3, 6, -(a ^ 2) * w, a * v], [3, 7, w, 0], [4, 5, [1, 0, 3, -w], u], [4, 6, -a * w, 0], [4, 7, [1, -w, 3, 0], u], [4, 8, a * w, 0], [4, 9, u + a ^ 2 * v, 0], [5, 6, -u, v], [5, 9, 0, w], [7, 8, u, -v], [7, 9, u, 0], [8, 9, -a * v, a ^ 2 * w]]], [u, v, w])
            end
        rep = [[f1, u], [f1, w], [f1, v], [f2, v, w], [f2, u, v], [f2, u, w], [f32, u, v, w], [f31, u, v], [f31, w, u], [f31, v, w], [f31, u, w], [f31, v, u], [f31, w, v], [f6, v, u, w], [f6, u, w, v], [f6, w, v, u], [f6, w, u, v], [f6, u, v, w], [f6, v, w, u], [f8, u, v, w], [f8, w, u, v], [f8, v, w, u], [f9, u, v, w, E(3)], [f9, u, v, w, E(3, 2)]]
        return ApplyFunc((rep[i])[1], (rep[i])[2:length(rep[i])]) * Product(para[1]) ^ 0
    end)
chevieset(:G25, :UnipotentCharacters, function ()
        local J
        J = E(3)
        return Dict{Symbol, Any}(:harishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "ST", :indices => 1:3, :rank => 3, :ST => 25), :levi => [], :parameterExponents => [1, 1, 1], :charNumbers => 1:24, :eigenvalue => 1, :cuspidalName => ""), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "ST", :indices => [3, 2], :rank => 2, :p => 3, :q => 1), :levi => [1], :parameterExponents => [1, 3], :charNumbers => [39, 31, 30, 41, 38, 40, 25, 27, 26], :eigenvalue => J ^ 2, :cuspidalName => ImprimitiveCuspidalName([[], [0, 1], [0, 1]])), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "ST", :indices => [2], :rank => 1, :p => 6, :q => 1), :levi => [1, 3], :parameterExponents => [[3, 3, 2, 0, 0, 2]], :charNumbers => [29, 28, 32, 44, 43, 33], :eigenvalue => J, :cuspidalName => Concatenation(ImprimitiveCuspidalName([[], [0, 1], [0, 1]]), "\\otimes ", ImprimitiveCuspidalName([[], [0, 1], [0, 1]]))), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "ST", :indices => [3], :rank => 1, :p => 3, :q => 1), :levi => 1:2, :parameterExponents => [[0, 4, 4]], :charNumbers => [42, 34, 35], :eigenvalue => -1, :cuspidalName => "G_4"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :parameterExponents => [], :charNumbers => [36], :eigenvalue => -J, :cuspidalName => "G_{25}[-\\zeta_3]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :parameterExponents => [], :charNumbers => [37], :eigenvalue => J, :cuspidalName => "G_{25}[\\zeta_3]")], :families => [Family("C1", [1]), Family(((CHEVIE[:families])[:X])(3), [12, 9, 25], Dict{Symbol, Any}(:signs => [1, 1, -1], :ennola => -2)), Family(((CHEVIE[:families])[:QZ])(3), [20, 16, 19, 6, 28, 26, 5, 27, 29], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, -1, 1, 1, 1], :special => 2, :cospecial => 3, :ennola => 5)), Family(((CHEVIE[:families])[:X])(6), [17, 23, 7, 24, 14, 32, 34, 30, 36, 8, 37, 31, 11, 35, 33], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1], :ennola => -15)), Family(((CHEVIE[:families])[:X])(3), [22, 21, 38], Dict{Symbol, Any}(:signs => [1, 1, -1], :ennola => 1)), Family(((CHEVIE[:families])[:X])(3), [15, 18, 39], Dict{Symbol, Any}(:signs => [1, 1, -1], :ennola => -3)), Family(SubFamilyij(((CHEVIE[:families])[:ExtPowCyclic])(6, 3), 1, 2, -(root(2)) // root(-1)), [3, 13, 40, 10, 41, 2, 43, 42, 4, 44], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, -1, 1, -1, 1, -1, -1], :cospecial => 6, :ennola => -9))], :a => [0, 12, 12, 12, 2, 2, 4, 4, 1, 12, 4, 1, 12, 4, 8, 2, 4, 8, 2, 2, 6, 6, 4, 4, 1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 6, 8, 12, 12, 12, 12, 12], :A => [0, 24, 24, 24, 16, 16, 20, 20, 11, 24, 20, 11, 24, 20, 22, 16, 20, 22, 16, 16, 21, 21, 20, 20, 11, 16, 16, 16, 16, 20, 20, 20, 20, 20, 20, 20, 20, 21, 22, 24, 24, 24, 24, 24])
    end)
chevieset(:G25, :Invariants, [function (x1, x2, x3)
            return ((-10 * x1 ^ 3 * x2 ^ 3 - 10 * x1 ^ 3 * x3 ^ 3) - 10 * x2 ^ 3 * x3 ^ 3) + x1 ^ 6 + x2 ^ 6 + x3 ^ 6
        end, function (x1, x2, x3)
            return ((((-(x1 ^ 3) * x2 ^ 6 + x1 ^ 3 * x3 ^ 6) - x2 ^ 3 * x3 ^ 6) + x1 ^ 6 * x2 ^ 3) - x1 ^ 6 * x3 ^ 3) + x2 ^ 6 * x3 ^ 3
        end, function (x1, x2, x3)
            return ((((2 * x1 ^ 3 * x2 ^ 3 * x3 ^ 6 + 2 * x1 ^ 3 * x2 ^ 6 * x3 ^ 3 + x1 ^ 3 * x2 ^ 9 + x1 ^ 3 * x3 ^ 9 + x2 ^ 3 * x3 ^ 9 + 2 * x1 ^ 6 * x2 ^ 3 * x3 ^ 3) - 4 * x1 ^ 6 * x2 ^ 6) - 4 * x1 ^ 6 * x3 ^ 6) - 4 * x2 ^ 6 * x3 ^ 6) + x1 ^ 9 * x2 ^ 3 + x1 ^ 9 * x3 ^ 3 + x2 ^ 9 * x3 ^ 3
        end])
chevieset(:G25, :Discriminant, function ()
        return function (t1, t2, t3)
                return ((36 * t1 * t2 ^ 2 * t3 - t1 ^ 2 * t3 ^ 2) - 32 * t3 ^ 3) + t1 ^ 3 * t2 ^ 2 + 108 * t2 ^ 4
            end
    end)
