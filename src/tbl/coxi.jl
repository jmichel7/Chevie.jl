
chevieset(:I, :CartanMat, function (arg...,)
        local bond, type_, m
        m = [[2, 0], [0, 2]]
        bond = arg[1]
        if bond == 2
            return m
        end
        if length(arg) == 2
            type_ = arg[2]
        elseif mod(bond, 2) == 0
            type_ = 1
        else
            type_ = E(2bond) + E(2bond, -1)
        end
        (m[1])[2] = -type_
        (m[2])[1] = (2 + E(bond) + E(bond, -1)) // (m[1])[2]
        return m
    end)
chevieset(:I, :PrintDiagram, function (arg...,)
        local bond, indices, type_
        print(arg[3], " ")
        bond = arg[1]
        indices = arg[2]
        if length(arg) == 4
            type_ = arg[4]
        else
            type_ = E(2bond) + E(2bond, -1)
        end
        if type_ == E(2bond) + E(2bond, -1)
            print(indices[1], " -", string(bond), "- ", indices[2], "\n")
        elseif type_ == 1
            print(indices[1], " >", string(bond), "> ", indices[2], "\n")
        else
            print(indices[1], " ?", string(bond), "? ", indices[2], "\n")
        end
    end)
chevieset(:I, :ReflectionName, function (arg...,)
        local bond, type_, opt
        bond = arg[1]
        opt = arg[2]
        if length(arg) == 3
            type_ = arg[3]
        elseif mod(bond, 2) == 0
            type_ = 1
        else
            type_ = E(2bond) + E(2bond, -1)
        end
        if type_ == 1
            if haskey(opt, :TeX)
                return SPrint("I_2(", bond, ")")
            elseif haskey(opt, :arg)
                return SPrint("\"I\",2,", bond)
            else
                return SPrint("I2(", bond, ")")
            end
        elseif type_ == E(2bond) + E(2bond, -1)
            if mod(bond, 2) == 1
                if haskey(opt, :TeX)
                    return SPrint("I_2(", bond, ")")
                elseif haskey(opt, :arg)
                    return SPrint("\"I\",2,", bond)
                else
                    return SPrint("I2(", bond, ")")
                end
            else
                if haskey(opt, :TeX)
                    return SPrint("I_{\\hbox{sym}2}(", bond, ")")
                elseif haskey(opt, :arg)
                    return SPrint("\"Isym\",2,", bond)
                else
                    return SPrint("Isym2(", bond, ")")
                end
            end
        elseif haskey(opt, :TeX)
            return SPrint("I_?(", Format(type_ ^ 2 // (2 + E(bond) + E(bond, -1)), opt), ")(", bond, ")")
        elseif haskey(opt, :arg)
            return SPrint("\"Isym\",2,", bond, ",", Format(type_ ^ 2 // (2 + E(bond) + E(bond, -1)), opt))
        else
            return SPrint("I?(", type_ ^ 2 // (2 + E(bond) + E(bond, -1)), ")(", bond, ")")
        end
    end)
chevieset(:I, :SemisimpleRank, 2)
chevieset(:I, :GeneratingRoots, function (m,)
        local a, b, r
        a = E(2m, m - 1)
        b = ComplexConjugate(a)
        if mod(m, 2) == 0
            r = ER(m // 2)
        else
            r = 1
        end
        return [[1, 0], [(r * (a + b)) // 2, ((r * (a - b)) // 2) // E(4)]]
    end)
chevieset(:I, :EigenvaluesGeneratingReflections, (m->begin
            [-1, -1]
        end))
chevieset(:I, :Size, function (arg...,)
        return 2 * arg[1]
    end)
chevieset(:I, :ReflectionDegrees, (m->begin
            [2, m]
        end))
chevieset(:I, :NrConjugacyClasses, (m->begin
            div(m + 3, 2) + mod(m + 1, 2) * 2
        end))
chevieset(:I, :ParabolicRepresentatives, function (m, s)
        return (chevieget(:imp, :ParabolicRepresentatives))(m, m, 2, s)
    end)
chevieset(:I, :CharName, function (m, x, option)
        local s
        if IsList(x[1])
            return PartitionTupleToString(x)
        else
            if haskey(option, :TeX)
                s = "\\phi"
            else
                s = "phi"
            end
            s = SPrint(s, "{", x[1], ",", x[2], "}")
            if length(x) == 3
                s = Append(s, x[3])
            end
            return string(s)
        end
    end)
chevieset(:I, :WordsClassRepresentatives, function (m,)
        local r, x, i
        if IsInt(m // 2)
            r = [[], [1], [2]]
        else
            r = [[], [1]]
        end
        x = [1, 2]
        for i = 1:div(m, 2)
            push!(r, copy(x))
            x = Append(x, [1, 2])
        end
        return r
    end)
chevieset(:I, :ClassInfo, function (m,)
        local r, i, clnp, cl, g1, g2, gen, perm, m1
        r = (chevieget(:I, :WordsClassRepresentatives))(m)
        clnp = map(IntListToString, r)
        g1 = Perm()
        i = 2
        while 2i <= m + 1
            g1 = g1 * Perm(i, (m - i) + 2)
            i = i + 1
        end
        g2 = Perm()
        i = 1
        while 2i <= m
            g2 = g2 * Perm(i, (m - i) + 1)
            i = i + 1
        end
        gen = [g1, g2]
        perm = function (l,)
                if length(l) == 0
                    return Perm()
                else
                    return Product(gen[l])
                end
            end
        m1 = div(m, 2)
        if mod(m, 2) == 0
            cl = [1, m1, m1]
            cl = Append(cl, fill(0, max(0, (1 + (m1 - 1)) - 1)) + 2)
            push!(cl, 1)
        else
            cl = [1, m]
            cl = Append(cl, fill(0, max(0, (1 + m1) - 1)) + 2)
        end
        return Dict{Symbol, Any}(:classtext => r, :classnames => clnp, :classparams => clnp, :orders => map((i->begin
                                order(perm(i))
                            end), r), :classes => cl)
    end)
chevieset(:I, :HeckeCharTable, function (m, param, rootparam)
        local u, v, squv, cl, r, ct, tbl
        u = -((param[1])[1]) // (param[1])[2]
        v = -((param[2])[1]) // (param[2])[2]
        if mod(m, 2) != 0
            squv = u
        elseif rootparam[1] !== nothing && rootparam[2] !== nothing
            squv = rootparam[1] * rootparam[2]
        else
            squv = GetRoot(u * v, 2, "CharTable(Hecke(I2(", m, ")))")
        end
        ct = [[u, v]]
        if mod(m, 2) == 0
            ct = Append(ct, [[u, -(u ^ 0)], [-(v ^ 0), v]])
        end
        push!(ct, [-(v ^ 0), -(v ^ 0)])
        cl = (chevieget(:I, :ClassInfo))(m)
        r = cl[:classtext]
        ct = map((i->begin
                        map((x->begin
                                    Product(i[x])
                                end), r)
                    end), ct)
        ct = Append(ct, map((j->begin
                            map(function (i,)
                                    local k
                                    k = length(r[i]) // 2
                                    if r[i] == []
                                        return 2 * v ^ 0
                                    elseif r[i] == [1]
                                        return u - 1
                                    elseif r[i] == [2]
                                        return v - 1
                                    else
                                        return squv ^ k * (E(m, k * j) + E(m, -k * j))
                                    end
                                end, 1:length(r))
                        end), 1:div(m - 1, 2)))
        tbl = Dict{Symbol, Any}(:identifier => SPrint("H(I2(", m, "))"), :cartan => CartanMat("I", 2, m), :size => 2m, :irredinfo => map((x->begin
                                Dict{Symbol, Any}(:charparam => x, :charname => (chevieget(:I, :CharName))(m, x, Dict{Symbol, Any}(:TeX => true)))
                            end), ((chevieget(:I, :CharInfo))(m))[:charparams]), :parameter => [u, v], :powermap => [], :irreducibles => ct * v ^ 0)
        Inherit(tbl, cl)
        tbl[:centralizers] = map((i->begin
                        tbl[:size] // i
                    end), tbl[:classes])
        tbl = ((CHEVIE[:compat])[:MakeCharacterTable])(tbl)
        ((CHEVIE[:compat])[:AdjustHeckeCharTable])(tbl, param)
        return tbl
    end)
chevieset(:I, :Representation, function (m, i)
        return (chevieget(:I, :HeckeRepresentation))(m, [[1, -1], [1, -1]], [1, 1], i)
    end)
chevieset(:I, :HeckeRepresentation, function (m, param, rootparam, i)
        local u, v, squv
        if i == 1
            return [[[(param[1])[1]]], [[(param[2])[1]]]]
        end
        if mod(m, 2) == 0
            i = i - 2
        end
        if i == 0
            return [[[(param[1])[1]]], [[(param[2])[2]]]]
        elseif i == 1
            return [[[(param[1])[2]]], [[(param[2])[1]]]]
        elseif i == 2
            return [[[(param[1])[2]]], [[(param[2])[2]]]]
        else
            u = -((param[1])[1]) // (param[1])[2]
            v = -((param[2])[1]) // (param[2])[2]
            if mod(m, 2) != 0
                squv = u
            elseif rootparam[1] !== nothing && rootparam[2] !== nothing
                squv = rootparam[1] * rootparam[2]
            else
                squv = GetRoot(u * v, 2, "Representation(Hecke(I2(", m, ")),[", i, "])")
            end
            return [-([[-(u ^ 0), u ^ 0], [0u, u]]) * (param[1])[2], -([[v, 0v], [u + v + squv * (E(m, i - 2) + E(m, 2 - i)), -(v ^ 0)]]) * (param[2])[2]]
        end
    end)
chevieset(:I, :Frobenius, function (m, sqrtu, j)
        return [[0, (1 // sqrtu) // (E(2m, j) + E(2m, -j))], [sqrtu * (E(2m, j) + E(2m, -j)), 0]] * sqrtu ^ 0
    end)
chevieset(:I, :PoincarePolynomial, function (m, param)
        local u, v
        u = -((param[1])[1]) // (param[1])[2]
        v = -((param[2])[1]) // (param[2])[2]
        if IsInt(m // 2)
            return Sum(1:m // 2, (i->begin
                                (u * v) ^ (i - 1)
                            end)) * (u + 1) * (v + 1)
        else
            return Sum(1:m, (i->begin
                                u ^ (i - 1)
                            end)) * (u + 1)
        end
    end)
chevieset(:I, :SchurElement, function (m, phi, para, rootpara)
        local u, v, ruv, e, ci
        if mod(m, 2) == 1
            ci = (chevieget(:I, :CharInfo))(m)
            ci = (ci[:malleParams])[Position(ci[:charparams], phi)]
            return (chevieget(:imp, :SchurElement))(m, 1, 2, ci, [map((i->begin
                                        E(m, i)
                                    end), 0:m - 1), para[2]], []) // m
        end
        u = -((para[1])[1]) // (para[1])[2]
        v = -((para[2])[1]) // (para[2])[2]
        if phi[1] == 1
            if phi[2] == m // 2
                e = (Sum(0:m // 2 - 1, (i->begin
                                        (u // v) ^ i
                                    end)) * (u + 1) * (v + 1)) // v
                if phi[3] == "'"
                    return e
                else
                    return (v // u) ^ (m // 2) * e
                end
            else
                e = Sum(0:m // 2 - 1, (i->begin
                                    (u * v) ^ i
                                end)) * (u + 1) * (v + 1)
                if phi[2] == 0
                    return e
                else
                    return (u * v) ^ (-m // 2) * e
                end
            end
        else
            e = E(m, phi[2]) + E(m, -(phi[2]))
            if all((i->begin
                            rootpara[i] !== nothing
                        end), [1, 2])
                ruv = Product(rootpara)
            else
                ruv = GetRoot(u * v, 2, "SchurElement(Hecke(I2(", m, "),", phi, "))")
            end
            return (-m * ((u * v + 1) - ruv * e) * (u + v + e * ruv)) // (u * v * (e ^ 2 - 4))
        end
    end)
chevieset(:I, :FakeDegree, function (m, phi, q)
        if phi[1] == 1
            return q ^ phi[2]
        else
            return q ^ phi[2] + q ^ (m - phi[2])
        end
    end)
chevieset(:I, :CharTable, function (m,)
        local res
        res = (chevieget(:I, :HeckeCharTable))(m, [[1, -1], [1, -1]], [1, 1])
        res[:identifier] = SPrint("W(I2(", m, "))")
        return res
    end)
chevieset(:I, :DecompositionMatrix, function (n, p)
        local T, m
        T = (chevieget(:I, :CharTable))(n)
        T[:name] = T[:identifier]
        m = DecompositionMatrix(mod(T, p))
        return map((c->begin
                        [c[1], (m[c[1]])[c[2]]]
                    end), BlocksMat(m))
    end)
chevieset(:I, :FactorizedSchurElement, function (arg...,)
        local ci
        if mod(arg[1], 2) == 0 && (arg[3])[1] != (arg[3])[2]
            error(" !  implemented")
        end
        ci = (chevieget(:I, :CharInfo))(arg[1])
        ci = (ci[:malleParams])[Position(ci[:charparams], arg[2])]
        return (chevieget(:imp, :FactorizedSchurElement))(arg[1], arg[1], 2, ci, arg[3], 1)
    end)
chevieset(:I, :Invariants, function (arg...,)
        local e, type_, m
        e = arg[1]
        if length(arg) == 2
            type_ = arg[2]
        elseif mod(e, 2) == 0
            type_ = 1
        else
            type_ = -(E(e, (e + 1) // 2)) - E(e, (e + 3) // 2)
        end
        m = DiagonalMat(1 + E(e, -1), -type_) * (chevieget(:imp, :GeneratingRoots))(e, e, 2)
        return map((f->begin
                        function (arg...,)
                            return ApplyFunc(f, arg * m)
                        end
                    end), (chevieget(:imp, :Invariants))(e, e, 2))
    end)
chevieset(:I, :SymbolToParameter, function (S,)
        if S[1] != [0, 1] || !([]) in S
            return false
        end
        if mod(length(S), 2) == 1
            S = reverse(S)
            return [Position(S, []), Position(S, [0, 1]) - Position(S, [])]
        else
            return (Position(S, []) + [-(Position(S[2:length(S)], [0, 1])), 0]) - 1
        end
    end)
chevieset(:I, :ParameterToSymbol, function (e, p)
        local S
        if p == [0]
            S = map((x->begin
                            [0]
                        end), 1:e)
            S[e] = [2]
        elseif p == [1]
            S = map((x->begin
                            [0, 1]
                        end), 1:e)
            S[e] = [1, 2]
        elseif length(p) == 3
            S = map((x->begin
                            [0]
                        end), 1:e // 2 - 1)
            S = Append(S, [[1], 2, (p[3] + 1) // 2])
        elseif mod(e, 2) == 0
            S = map((x->begin
                            [0]
                        end), 1:e)
            if p[1] == 0
                S[[e, e - p[2]]] = [[1], [1]]
            else
                S[1 + [0, mod(p[2] - p[1], e)]] = [[0, 1], [0, 1]]
                S[1 + [mod(-(p[1]), e), p[2]]] = [[], []]
            end
        else
            S = map((i->begin
                            [0]
                        end), 1:e)
            if p[1] != 0
                S[1 + [0, mod(-(Sum(p)), e)]] = [[0, 1], [0, 1]]
                S[1 + map((x->begin
                                        mod(x, e)
                                    end), -p)] = [[], []]
            else
                S[e + [-(mod(p[2] - p[1], e)), 0]] = [[1], [1]]
            end
        end
        return S
    end)
chevieset(:I, :UnipotentCharacters, function (e,)
        local cusp, uc, f
        f = div(e, 2)
        uc = Dict{Symbol, Any}()
        uc[:harishChandra] = [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "I", :indices => [1, 2], :rank => 2, :bond => e), :parameterExponents => [1, 1], :levi => [], :eigenvalue => 1, :cuspidalName => "")]
        if mod(e, 2) != 0
            ((uc[:harishChandra])[1])[:charNumbers] = 1:f + 2
        else
            ((uc[:harishChandra])[1])[:charNumbers] = Concatenation([1, 3, 4, 2], 4 + (1:f - 1))
        end
        cusp = Concatenation(map((k->begin
                            map((l->begin
                                        [k, l]
                                    end), k + 1:(e - k) - 1)
                        end), 1:f - 1))
        f = (f + 1) - mod(e, 2)
        uc[:harishChandra] = Append(uc[:harishChandra], map((x->begin
                            Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :parameterExponents => [], :levi => [1, 2], :eigenvalue => E(e, -(Product(cusp[x]))), :cuspidalName => SPrint("I_2(", e, ")", FormatGAP(cusp[x])), :charNumbers => [x + f + 2])
                        end), 1:length(cusp)))
        uc[:families] = [Family(((CHEVIE[:families])[:Dihedral])(e), (1:length(cusp) + f) + 2), Family("C1", [1]), Family("C1", [2])]
        uc[:parameters] = Concatenation([[0], [1]], ((uc[:families])[1])[:parameters])
        uc[:charSymbols] = map((p->begin
                        (chevieget(:I, :ParameterToSymbol))(e, p)
                    end), uc[:parameters])
        uc[:a] = Concatenation([0, e], map((x->begin
                            1
                        end), ((uc[:families])[1])[:parameters]))
        uc[:A] = Concatenation([0, e], map((x->begin
                            e - 1
                        end), ((uc[:families])[1])[:parameters]))
        if e == 5
            uc[:curtis] = [2, 1, 3, 4, 6, 5]
        end
        return uc
    end)