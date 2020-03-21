
chevieset(Symbol("2I"), :ReflectionName, function (m, option)
        if haskey(option, :TeX)
            return SPrint("{}^2I_2(", m, ")")
        else
            return SPrint("2I2(", m, ")")
        end
    end)
chevieset(Symbol("2I"), :NrConjugacyClasses, (m->begin
            div(m + 3, 2)
        end))
chevieset(Symbol("2I"), :WordsClassRepresentatives, function (m,)
        local r, x, i
        r = [[]]
        x = [1]
        for i = 1:div(m + 1, 2)
            push!(r, copy(x))
            x = Append(x, [2, 1])
        end
        return r
    end)
chevieset(Symbol("2I"), :ClassInfo, function (m,)
        local res
        res = Dict{Symbol, Any}(:classtext => (chevieget(Symbol("2I"), :WordsClassRepresentatives))(m))
        res[:classnames] = map(IntListToString, res[:classtext])
        res[:classparams] = res[:classnames]
        res[:classes] = [m]
        res[:classes] = Append(res[:classes], fill(0, max(0, (1 + div(m, 2)) - 1)) + 2)
        if mod(m, 2) == 1
            push!(res[:classes], 1)
        end
        res[:orders] = map((i->begin
                        (2m) // gcd(2m, length(i))
                    end), res[:classtext])
        (res[:orders])[1] = 2
        return res
    end)
chevieset(Symbol("2I"), :PhiFactors, (m->begin
            [1, -1]
        end))
chevieset(Symbol("2I"), :ClassParameter, function (m, w)
        local l
        if IsInt(length(w) // 2)
            return ""
        else
            l = copy(w)
            if l[1] == 2
                l = l[2:length(l)]
                push!(l, 1)
            end
            return IntListToString(l)
        end
    end)
chevieset(Symbol("2I"), :CharName, function (m, x, option)
        local s
        if IsList(x[1])
            return PartitionTupleToString(x)
        else
            if haskey(option, :TeX)
                s = "phi"
            else
                s = "\\phi_"
            end
            s = SPrint(s, "{", x[1], ",", x[2], "}")
            if length(x) == 3
                s = Append(s, x[3])
            end
            return string(s)
        end
    end)
chevieset(Symbol("2I"), :CharInfo, function (m,)
        local res
        res = Dict{Symbol, Any}(:extRefl => [1, 3, 2])
        if m == 4
            res[:charparams] = [[[2], []], [[], [1, 1]], [[1], [1]]]
        else
            res[:charparams] = Concatenation([[1, 0], [1, m]], map((i->begin
                                [2, i]
                            end), 1:div(m - 1, 2)))
        end
        return res
    end)
chevieset(Symbol("2I"), :FakeDegree, function (m, phi, q)
        local i
        i = Position(((chevieget(Symbol("2I"), :CharInfo))(m))[:charparams], phi)
        if i == 1
            return q ^ 0
        elseif i == 2
            return q ^ m
        else
            return q ^ ((m + 2) - i) - q ^ (i - 2)
        end
    end)
chevieset(Symbol("2I"), :HeckeCharTable, function (m, param, sqrtparam)
        local q, i, j, ct, cos, cl, l, ident, ord, v, tbl
        q = -((param[1])[1]) // (param[1])[2]
        if m == 4
            ident = "2B"
        elseif m == 6
            ident = "2G"
        else
            ident = "2I2"
        end
        ident = SPrint(ident, "(", m, ")")
        if q != 1
            ident = SPrint("H(", ident, ")")
        end
        if !(sqrtparam[1] !== nothing)
            v = GetRoot(q, 2, "CharTable(", ident, ")")
        else
            v = sqrtparam[1]
        end
        cl = (chevieget(Symbol("2I"), :ClassInfo))(m)
        cos = (i->begin
                    E(2m, i) + E(2m, -i)
                end)
        ct = map((i->begin
                        map((j->begin
                                    cos(1) * 0 * v
                                end), cl[:classtext])
                    end), cl[:classtext])
        ct[1] = map((i->begin
                        q ^ length(i)
                    end), cl[:classtext])
        ct[2] = map((i->begin
                        (-1) ^ length(i)
                    end), cl[:classtext])
        for i = 1:div(m - 1, 2)
            for j = 1:div(m + 1, 2)
                (ct[i + 2])[j + 1] = -(v ^ (2j - 1)) * cos(i * (2j - 1))
            end
        end
        tbl = Dict{Symbol, Any}(:identifier => ident, :name => ident, :cartan => [[2, -(cos(1))], [-(cos(1)), 2]], :size => 2m, :parameter => [q, q], :sqrtparameter => [v, v], :irreducibles => ct * v ^ 0, :irredinfo => map((x->begin
                                Dict{Symbol, Any}(:charparam => x, :charname => (chevieget(Symbol("2I"), :CharName))(m, x, Dict{Symbol, Any}(:TeX => true)))
                            end), ((chevieget(Symbol("2I"), :CharInfo))(m))[:charparams]))
        Inherit(tbl, cl)
        tbl[:centralizers] = map((i->begin
                        tbl[:size] // i
                    end), cl[:classes])
        tbl = ((CHEVIE[:compat])[:MakeCharacterTable])(tbl)
        ((CHEVIE[:compat])[:AdjustHeckeCharTable])(tbl, param)
        return tbl
    end)
chevieset(Symbol("2I"), :CharTable, (m->begin
            (chevieget(Symbol("2I"), :HeckeCharTable))(m, [[1, -1], [1, -1]], [1, 1])
        end))
chevieset(Symbol("2I"), :Representation, function (m, i)
        return (chevieget(Symbol("2I"), :HeckeRepresentation))(m, [[1, -1], [1, -1]], [1, 1], i)
    end)
chevieset(Symbol("2I"), :HeckeRepresentation, function (m, param, sqrtparam, i)
        local v, q, e
        q = -((param[1])[1]) // (param[1])[2]
        if !(sqrtparam[1] !== nothing)
            v = GetRoot(q, 2, "Representation(Hecke(2I2(", m, ")),[", i, "])")
        else
            v = sqrtparam[1]
        end
        e = E(2m)
        if i == 1
            return Dict{Symbol, Any}(:gens => [[[v ^ 2]], [[v ^ 2]]], :F => [[1]])
        elseif i == 2
            return Dict{Symbol, Any}(:gens => [[[-1]], [[-1]]] * v ^ 0, :F => [[1]])
        else
            i = i - 2
            return Dict{Symbol, Any}(:gens => [[[-1, 0], [v * (e ^ i + e ^ -i), v ^ 2]], [[v ^ 2, v * (e ^ i + e ^ -i)], [0, -1]]] * v ^ 0, :F => -([[0, 1], [1, 0]]))
        end
    end)
chevieset(Symbol("2I"), :UnipotentCharacters, function (e,)
        local nc, uc, i, ac, c, n, symUnp, untUnp, TeXpref, eig
        uc = Dict{Symbol, Any}()
        n = div(e - 1, 2)
        nc = Concatenation(map((i->begin
                            map((j->begin
                                        [i, j]
                                    end), i + 1:(e - i) - 1)
                        end), 1:n))
        ac = Concatenation(map((l->begin
                            [0, l]
                        end), 1:n), nc)
        symUnp = Concatenation(map((i->begin
                            map((j->begin
                                        2 * [i, j] - 1
                                    end), i + 1:e - i)
                        end), 1:n))
        if mod(e, 2) == 1
            untUnp = map(function (s,)
                        local res
                        res = map((x->begin
                                        mod(div(x, 2), e)
                                    end), e - reverse(s))
                        if res[1] > res[2]
                            res = map((x->begin
                                            mod(x, e)
                                        end), -res)
                        end
                        return res
                    end, symUnp)
            SortParallel(untUnp, symUnp)
        end
        if e == 4
            TeXpref = "B_2"
        elseif e == 6
            TeXpref = "G_2"
        else
            TeXpref = SPrint("I_2(", e, ")")
        end
        uc[:harishChandra] = Concatenation([Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [1], :rank => 1), :parameterExponents => [e], :levi => [], :eigenvalue => 1, :cuspidalName => "", :charNumbers => [2, 1])], map((x->begin
                            Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :parameterExponents => [], :levi => [1, 2], :eigenvalue => E(2e, Product(x)), :cuspidalName => SPrint("{}^2", TeXpref, "[", x[1], ",", x[2], "]"), :charNumbers => [2 + Position(symUnp, x)])
                        end), symUnp))
        uc[:almostHarishChandra] = Concatenation([Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:orbit => [Dict{Symbol, Any}(:series => "I", :indices => [1, 2], :rank => 2, :bond => e)], :twist => perm"(1,2)"), :parameterExponents => [1, 1], :levi => [], :eigenvalue => 1, :cuspidalName => "", :charNumbers => 1:n + 2)], map((x->begin
                            Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :parameterExponents => [], :levi => [1, 2], :eigenvalue => E(e, -(Product(x))), :cuspidalName => SPrint(TeXpref, "[", x[1], ",", x[2], "]"), :charNumbers => [n + 2 + Position(nc, x)])
                        end), nc))
        if e == 4
            ((uc[:almostHarishChandra])[2])[:cuspidalName] = "B_2"
            ((((uc[:almostHarishChandra])[1])[:relativeType])[:orbit])[1] = Dict{Symbol, Any}(:series => "B", :indices => [1, 2], :rank => 2, :cartanType => ER(2))
        elseif e == 6
            eig = [E(3, 2), -1, E(3), 1]
            ((((uc[:almostHarishChandra])[1])[:relativeType])[:orbit])[1] = Dict{Symbol, Any}(:series => "G", :indices => [1, 2], :rank => 2, :cartanType => ER(3))
            for i = 1:4
                ((uc[:almostHarishChandra])[i + 1])[:cuspidalName] = SPrint("G2[", FormatTeX(eig[i]), "]")
            end
        end
        uc[:charParams] = Concatenation((((chevieget(:I, :CharInfo))(e))[:charparams])[[1, 2]], ac)
        uc[:almostCharSymbols] = Concatenation([map((x->begin
                                [0]
                            end), 1:e), map((x->begin
                                [0, 1]
                            end), 1:e)], map(function (s,)
                        local S, k, l
                        S = map((i->begin
                                        [0]
                                    end), 1:e)
                        k = s[1]
                        l = s[2]
                        S[k + 1] = []
                        S[l + 1] = []
                        push!(S[1], 1)
                        push!(S[k + l + 1], 1)
                        return S
                    end, ac))
        ((uc[:almostCharSymbols])[1])[e] = [2]
        ((uc[:almostCharSymbols])[2])[e] = [1, 2]
        uc[:charSymbols] = Concatenation([map((x->begin
                                [0]
                            end), 1:e), map((x->begin
                                [0, 1]
                            end), 1:e)], map(function (s,)
                        local S, k, l
                        k = div(s[1] + 1, 2)
                        l = div(s[2] + 1, 2)
                        S = map(function (i,)
                                    if i == k + 1 || i == l + 1
                                        return []
                                    else
                                        return [0]
                                    end
                                end, 1:e)
                        push!(S[1], 1)
                        push!(S[k + l], 1)
                        return S
                    end, symUnp))
        ((uc[:charSymbols])[1])[[1, 2]] = [[0, 2], []]
        ((uc[:charSymbols])[2])[[1, 2]] = [[0, 1, 2], [1]]
        c = (a->begin
                    E(2e, a) + E(2e, -a)
                end)
        uc[:families] = [Family("C1", [1]), Family("C1", [2]), Family(Dict{Symbol, Any}(:eigenvalues => map((s->begin
                                        E(2e, s[1] * s[2])
                                    end), symUnp), :fourierMat => map((j->begin
                                        map((i->begin
                                                    (c(i * reverse(j)) - c(i * j)) // e
                                                end), symUnp)
                                    end), ac), :sh => map((s->begin
                                        E(e, -(s[1]) * s[2])
                                    end), ac), :charNumbers => 2 + (1:length(ac)), :special => 1))]
        uc[:a] = Concatenation([0, e], map((x->begin
                            1
                        end), ac))
        uc[:A] = Concatenation([0, e], map((x->begin
                            e - 1
                        end), ac))
        if e == 5
            (uc[:families])[3] = (uc[:families])[3] ^ 13
            for c = uc[:harishChandra]
                c[:eigenvalue] = GaloisCyc(c[:eigenvalue], 13)
            end
        end
        return uc
    end)
