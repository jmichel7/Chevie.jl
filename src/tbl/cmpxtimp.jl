
chevieset(:timp, :PhiFactors, function (p, q, r, phi)
        local res, o
        o = order(phi)
        if p == q
            if mod(p, o) == 0
                if phi == perm"(1,2,4)"
                    return [1, 1, E(3, 2)]
                else
                    res = map((x->begin
                                    E(1)
                                end), 1:r)
                    res[length(res)] = E(o)
                    return res
                end
            elseif [p, q, r, o] == [3, 3, 3, 4]
                return [E(4), 1, -(E(4))]
            end
        elseif p == 2q
            return [1, E(3)]
        end
        error("wrong arguments")
    end)
chevieset(:timp, :ClassInfo, function (p, q, r, phi)
        if [p, q, r] == [3, 3, 3]
            if phi == perm"(1,4,2)"
                return Dict{Symbol, Any}(:classes => [3, 9, 9, 3, 18, 9, 3], :classtext => [[], [3], [1], [2, 1], [3, 1], [3, 2, 1], [2, 3, 1, 2, 3, 1]], :classparams => [[], [3], [1], [2, 1], [3, 1], [3, 2, 1], [2, 3, 1, 2, 3, 1]], :classnames => ["Id", "3", "1", "21", "31", "321", "231231"])
            elseif phi == perm"(1,2,4)"
                return Dict{Symbol, Any}(:classes => [3, 9, 9, 3, 18, 9, 3], :classtext => [[], [3], [1], [1, 2], [3, 1], [3, 1, 2], [1, 3, 2, 1, 3, 2]], :classparams => [[], [3], [1], [1, 2], [3, 1], [3, 1, 2], [1, 3, 2, 1, 3, 2]], :classnames => ["Id", "3", "1", "12", "31", "312", "132132"])
            elseif phi == perm"(1,2,3,4)"
                return Dict{Symbol, Any}(:classes => [9, 9, 9, 9, 9, 9], :classtext => [[], [1], [1, 2], [1, 2, 3], [1, 2, 1], [1, 2, 1, 3]], :classparams => [[], [1], [1, 2], [1, 2, 3], [1, 2, 1], [1, 2, 1, 3]], :classnames => ["Id", "1", "12", "123", "121", "1213"])
            else
                error("should not happen")
            end
        elseif [p, q, r] == [4, 2, 2]
            return Dict{Symbol, Any}(:classtext => [[], [1], [1, 2, 3, 1, 2, 3], [1, 2, 3, 1, 2, 3, 1, 2, 3]], :classes => [4, 4, 4, 4], :classnames => ["Id", "1", "cc", "z"])
        else
            ChevieErr("ClassInfo not implemented")
            return false
        end
    end)
chevieset(:timp, :NrConjugacyClasses, function (p, q, r, phi)
        return length(((chevieget(:timp, :ClassInfo))(p, q, r, phi))[:classtext])
    end)
chevieset(:timp, :CharInfo, function (p, q, r, phi)
        if [p, q, r] == [3, 3, 3]
            if phi == perm"(1,4,2)" || phi == perm"(1,2,4)"
                return Dict{Symbol, Any}(:charparams => [[[], [], [3]], [[], [], [1, 1, 1]], [[], [], [2, 1]], [[], [1, 1], [1]], [[], [2], [1]], [[], [1], [1, 1]], [[], [1], [2]]], :extRefl => [1, 5, 6, 2])
            elseif phi == perm"(1,2,3,4)"
                return Dict{Symbol, Any}(:charparams => [[[], [], [3]], [[], [], [1, 1, 1]], [[], [1, 1], [1]], [[], [2], [1]], [[], [1], [1, 1]], [[], [1], [2]]], :extRefl => [1, 4, 5, 2])
            else
                error("phi=",phi,": should not happen")
            end
        elseif [p, q, r] == [4, 2, 2]
            return Dict{Symbol, Any}(:charparams => [[[], [], [2], []], [[], [], [], [1, 1]], [[], [1], [1], []], [[], [], [1], [1]]], :extRefl => [1, 4, 2])
        else
            ChevieErr("CharInfo not implemented")
            return false
        end
    end)
chevieset(:timp, :CharTable, function (p, q, r, phi)
        local res
        if [p, q, r] == [3, 3, 3]
            if phi == perm"(1,4,2)"
                res = Dict{Symbol, Any}(:size => 54, :order => 54, :centralizers => [18, 6, 6, 18, 3, 6, 18], :identifier => "3'G(3,3,3)", :name => "3'G(3,3,3)", :classes => [3, 9, 9, 3, 18, 9, 3], :irreducibles => [[1, 1, 1, 1, 1, 1, 1], [1, -1, -1, 1, 1, -1, 1], [2, 0, 0, 2, -1, 0, 2], E(3) * [-(root(-3)), -1, -(E(3, 2)), 2 * E(3) + E(3, 2), 0, -(E(3)), -(E(3)) - 2 * E(3, 2)], E(3) * [-(root(-3)), 1, E(3, 2), 2 * E(3) + E(3, 2), 0, E(3), -(E(3)) - 2 * E(3, 2)], [-2 * E(3) - E(3, 2), -(E(3, 2)), -1, root(-3), 0, -(E(3)), E(3) + 2 * E(3, 2)], [-2 * E(3) - E(3, 2), E(3, 2), 1, root(-3), 0, E(3), E(3) + 2 * E(3, 2)]])
            elseif phi == perm"(1,2,4)"
                res = Dict{Symbol, Any}(:size => 54, :order => 54, :centralizers => [18, 6, 6, 18, 3, 6, 18], :identifier => "3G(3,3,3)", :name => "3G(3,3,3)", :classes => [3, 9, 9, 3, 18, 9, 3], :irreducibles => [[1, 1, 1, 1, 1, 1, 1], [1, -1, -1, 1, 1, -1, 1], [2, 0, 0, 2, -1, 0, 2], E(3, 2) * [root(-3), -1, -(E(3)), (-3 - root(-3)) // 2, 0, -(E(3, 2)), (3 - root(-3)) // 2], E(3, 2) * [root(-3), 1, E(3), (-3 - root(-3)) // 2, 0, E(3, 2), (3 - root(-3)) // 2], [(3 + root(-3)) // 2, -(E(3)), -1, -(root(-3)), 0, -(E(3, 2)), (-3 + root(-3)) // 2], [(3 + root(-3)) // 2, E(3), 1, -(root(-3)), 0, E(3, 2), (-3 + root(-3)) // 2]])
            elseif phi == perm"(1,2,3,4)"
                res = Dict{Symbol, Any}(:size => 54, :order => 54, :centralizers => [6, 6, 6, 6, 6, 6], :identifier => "4G(3,3,3)", :name => "4G(3,3,3)", :classes => [9, 9, 9, 9, 9, 9], :irreducibles => [[1, 1, 1, 1, 1, 1], [1, -1, 1, -1, -1, 1], [1, E(3), E(3, 2), 1, E(3, 2), E(3)], [1, -(E(3)), E(3, 2), -1, -(E(3, 2)), E(3)], [1, E(3, 2), E(3), 1, E(3), E(3, 2)], [1, -(E(3, 2)), E(3), -1, -(E(3)), E(3, 2)]])
            else
                error("should not happen")
            end
        elseif [p, q, r] == [4, 2, 2]
            res = Dict{Symbol, Any}(:size => 16, :order => 16, :centralizers => [4, 4, 4, 4], :classes => [4, 4, 4, 4], :identifier => "3G(4,2,2)", :name => "3G(4,2,2)", :irreducibles => [[1, 1, 1, 1], [1, -1, 1, -1], [-1, E(4), 1, -(E(4))], [-1, -(E(4)), 1, E(4)]])
        else
            ChevieErr("CharTable not implemented")
            return false
        end
        res[:text] = "origin: Dixon's Algorithm"
        return ((CHEVIE[:compat])[:MakeCharacterTable])(res)
    end)
chevieset(:timp, :UnipotentCharacters, function (p, q, r, phi)
        local res, a
        if [p, q, r] == [3, 3, 3]
            if phi == perm"(1,4,2)"
                return Dict{Symbol, Any}(:harishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "ST", :indices => [1, 2], :p => 3, :q => 1, :rank => 2), :levi => [], :eigenvalue => 1, :parameterExponents => [[2, 0, 1], 1], :cuspidalName => "", :charNumbers => [7, 3, 5, 2, 4, 9, 1, 6, 8])], :almostHarishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:orbit => [Dict{Symbol, Any}(:series => "ST", :indices => 1:3, :p => 3, :q => 3, :rank => 3)], :twist => perm"(1,2,4)"), :levi => [], :eigenvalue => 1, :cuspidalName => "", :charNumbers => 1:7), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :charNumbers => [8], :eigenvalue => E(3), :cuspidalName => "G_{3,3,3}[\\zeta_3]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :charNumbers => [9], :eigenvalue => E(3, 2), :cuspidalName => "G_{3,3,3}[\\zeta_3^2]")], :families => [Family("C1", [1]), Family("C1", [2]), Family("C1", [3]), Family(ComplexConjugate(((CHEVIE[:families])[:X])(3)), [7, 5, 8]), Family(ComplexConjugate(((CHEVIE[:families])[:X])(3)), [6, 4, 9])], :a => [0, 9, 3, 4, 1, 4, 1, 1, 4], :A => [0, 9, 6, 8, 5, 8, 5, 5, 8])
            elseif phi == perm"(1,2,4)"
                return Dict{Symbol, Any}(:harishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "ST", :indices => [1, 2], :p => 3, :q => 1, :rank => 2), :levi => [], :eigenvalue => 1, :parameterExponents => [[2, 1, 0], 1], :cuspidalName => "", :charNumbers => [7, 5, 3, 9, 4, 2, 1, 8, 6])], :almostHarishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:orbit => [Dict{Symbol, Any}(:series => "ST", :indices => 1:3, :p => 3, :q => 3, :rank => 3)], :twist => perm"(1,2,4)"), :levi => [], :eigenvalue => 1, :cuspidalName => "", :charNumbers => 1:7), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :charNumbers => [8], :eigenvalue => E(3), :cuspidalName => "G_{3,3,3}[\\zeta_3]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :charNumbers => [9], :eigenvalue => E(3, 2), :cuspidalName => "G_{3,3,3}[\\zeta_3^2]")], :families => [Family("C1", [1]), Family("C1", [2]), Family("C1", [3]), Family(((CHEVIE[:families])[:X])(3), [7, 5, 8]), Family(((CHEVIE[:families])[:X])(3), [6, 4, 9])], :a => [0, 9, 3, 4, 1, 4, 1, 1, 4], :A => [0, 9, 6, 8, 5, 8, 5, 5, 8])
            elseif phi == perm"(1,2,3,4)"
                res = Dict{Symbol, Any}(:harishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "ST", :indices => [1], :p => 6, :q => 1, :rank => 1), :levi => [], :eigenvalue => 1, :parameterExponents => [[3, 1, 2, 0, 2, 1]], :cuspidalName => "", :charNumbers => [1, 5, 4, 2, 6, 3]), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :parameterExponents => [], :charNumbers => [7], :eigenvalue => E(3, 2), :cuspidalName => "{}^4G_{3,3,3}[\\zeta_3^2]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :parameterExponents => [], :charNumbers => [8], :eigenvalue => E(3), :cuspidalName => "{}^4G_{3,3,3}[\\zeta_3]")], :almostHarishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:orbit => [Dict{Symbol, Any}(:series => "ST", :indices => 1:3, :p => 3, :q => 3, :rank => 3)], :twist => perm"(1,2,3,4)"), :levi => [], :eigenvalue => 1, :cuspidalName => "", :charNumbers => 1:6), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :charNumbers => [7], :eigenvalue => E(3), :cuspidalName => "G_{3,3,3}[\\zeta_3]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:3, :charNumbers => [8], :eigenvalue => E(3, 2), :cuspidalName => "G_{3,3,3}[\\zeta_3^2]")], :families => [Family("C1", [1]), Family("C1", [2]), Family(ComplexConjugate(((CHEVIE[:families])[:X])(3)), [3, 5, 7], Dict{Symbol, Any}(:signs => [1, 1, -1])), Family(((CHEVIE[:families])[:X])(3), [4, 6, 8], Dict{Symbol, Any}(:signs => [1, 1, -1]))], :a => [0, 9, 4, 1, 4, 1, 4, 1], :A => [0, 9, 8, 5, 8, 5, 8, 5])
                (((res[:families])[3])[:eigenvalues])[3] = E(3, 2)
                (((res[:families])[4])[:eigenvalues])[3] = E(3)
                a = 1
                (((res[:families])[3])[:fourierMat])[3] = a * (((res[:families])[3])[:fourierMat])[3]
                (((res[:families])[4])[:fourierMat])[3] = galois(a, -1) * (((res[:families])[4])[:fourierMat])[3]
                return res
            else
                error("should not happen")
                return false
            end
        else
            if q == 1 || q == p
                ChevieErr("UnipotentCharacters not implemented")
            end
            return false
        end
    end)