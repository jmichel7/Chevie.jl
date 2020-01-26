
chevieset(:I, :CartanMat, function (arg...,)
        #= none:17 =#
        local bond, type_, m
        #= none:19 =#
        m = [[2, 0], [0, 2]]
        #= none:21 =#
        bond = arg[1]
        #= none:23 =#
        if bond == 2
            #= none:23 =#
            return m
        end
        #= none:25 =#
        if length(arg) == 2
            #= none:25 =#
            type_ = arg[2]
        elseif #= none:27 =# mod(bond, 2) == 0
            #= none:27 =#
            type_ = 1
        else
            #= none:29 =#
            type_ = E(2bond) + E(2bond, -1)
        end
        #= none:31 =#
        (m[1])[2] = -type_
        #= none:32 =#
        (m[2])[1] = (2 + E(bond) + E(bond, -1)) // (m[1])[2]
        #= none:33 =#
        return m
    end)
chevieset(:I, :PrintDiagram, function (arg...,)
        #= none:3 =#
        local bond, indices, type_
        #= none:5 =#
        print(arg[3], " ")
        #= none:6 =#
        bond = arg[1]
        #= none:7 =#
        indices = arg[2]
        #= none:9 =#
        if length(arg) == 4
            #= none:9 =#
            type_ = arg[4]
        else
            #= none:10 =#
            type_ = E(2bond) + E(2bond, -1)
        end
        #= none:12 =#
        if type_ == E(2bond) + E(2bond, -1)
            #= none:13 =#
            print(indices[1], " -", string(bond), "- ", indices[2], "\n")
        elseif #= none:15 =# type_ == 1
            #= none:15 =#
            print(indices[1], " >", string(bond), "> ", indices[2], "\n")
        else
            #= none:17 =#
            print(indices[1], " ?", string(bond), "? ", indices[2], "\n")
        end
    end)
chevieset(:I, :ReflectionName, function (arg...,)
        #= none:3 =#
        local bond, type_, opt
        #= none:5 =#
        bond = arg[1]
        #= none:6 =#
        opt = arg[2]
        #= none:8 =#
        if length(arg) == 3
            #= none:8 =#
            type_ = arg[3]
        elseif #= none:10 =# mod(bond, 2) == 0
            #= none:10 =#
            type_ = 1
        else
            #= none:12 =#
            type_ = E(2bond) + E(2bond, -1)
        end
        #= none:14 =#
        if type_ == 1
            #= none:15 =#
            if haskey(opt, :TeX)
                #= none:15 =#
                return SPrint("I_2(", bond, ")")
            elseif #= none:17 =# haskey(opt, :arg)
                #= none:17 =#
                return SPrint("\"I\",2,", bond)
            else
                #= none:19 =#
                return SPrint("I2(", bond, ")")
            end
        elseif #= none:21 =# type_ == E(2bond) + E(2bond, -1)
            #= none:22 =#
            if mod(bond, 2) == 1
                #= none:23 =#
                if haskey(opt, :TeX)
                    #= none:23 =#
                    return SPrint("I_2(", bond, ")")
                elseif #= none:25 =# haskey(opt, :arg)
                    #= none:25 =#
                    return SPrint("\"I\",2,", bond)
                else
                    #= none:27 =#
                    return SPrint("I2(", bond, ")")
                end
            else
                #= none:30 =#
                if haskey(opt, :TeX)
                    #= none:30 =#
                    return SPrint("I_{\\hbox{sym}2}(", bond, ")")
                elseif #= none:32 =# haskey(opt, :arg)
                    #= none:32 =#
                    return SPrint("\"Isym\",2,", bond)
                else
                    #= none:34 =#
                    return SPrint("Isym2(", bond, ")")
                end
            end
        elseif #= none:37 =# haskey(opt, :TeX)
            #= none:38 =#
            return SPrint("I_?(", Format(type_ ^ 2 // (2 + E(bond) + E(bond, -1)), opt), ")(", bond, ")")
        elseif #= none:41 =# haskey(opt, :arg)
            #= none:41 =#
            return SPrint("\"Isym\",2,", bond, ",", Format(type_ ^ 2 // (2 + E(bond) + E(bond, -1)), opt))
        else
            #= none:44 =#
            return SPrint("I?(", type_ ^ 2 // (2 + E(bond) + E(bond, -1)), ")(", bond, ")")
        end
    end)
chevieset(:I, :SemisimpleRank, 2)
chevieset(:I, :GeneratingRoots, function (m,)
        #= none:3 =#
        local a, b, r
        #= none:5 =#
        a = E(2m, m - 1)
        #= none:7 =#
        b = ComplexConjugate(a)
        #= none:9 =#
        if mod(m, 2) == 0
            #= none:9 =#
            r = ER(m // 2)
        else
            #= none:10 =#
            r = 1
        end
        #= none:12 =#
        return [[1, 0], [(r * (a + b)) // 2, ((r * (a - b)) // 2) // E(4)]]
    end)
chevieset(:I, :EigenvaluesGeneratingReflections, (m->begin
            #= none:3 =#
            [-1, -1]
        end))
chevieset(:I, :Size, function (arg...,)
        #= none:3 =#
        return 2 * arg[1]
    end)
chevieset(:I, :ReflectionDegrees, (m->begin
            #= none:3 =#
            [2, m]
        end))
chevieset(:I, :NrConjugacyClasses, (m->begin
            #= none:3 =#
            div(m + 3, 2) + mod(m + 1, 2) * 2
        end))
chevieset(:I, :ParabolicRepresentatives, function (m, s)
        #= none:4 =#
        return (chevieget(:imp, :ParabolicRepresentatives))(m, m, 2, s)
    end)
chevieset(:I, :CharName, function (m, x, option)
        #= none:4 =#
        local s
        #= none:6 =#
        if IsList(x[1])
            #= none:6 =#
            return PartitionTupleToString(x)
        else
            #= none:9 =#
            if haskey(option, :TeX)
                #= none:9 =#
                s = "\\phi"
            else
                #= none:10 =#
                s = "phi"
            end
            #= none:12 =#
            s = SPrint(s, "{", x[1], ",", x[2], "}")
            #= none:14 =#
            if length(x) == 3
                #= none:14 =#
                s = Append(s, x[3])
            end
            #= none:16 =#
            return string(s)
        end
    end)
chevieset(:I, :WordsClassRepresentatives, function (m,)
        #= none:3 =#
        local r, x, i
        #= none:5 =#
        if IsInt(m // 2)
            #= none:5 =#
            r = [[], [1], [2]]
        else
            #= none:6 =#
            r = [[], [1]]
        end
        #= none:8 =#
        x = [1, 2]
        #= none:10 =#
        for i = 1:div(m, 2)
            #= none:10 =#
            push!(r, copy(x))
            #= none:11 =#
            x = Append(x, [1, 2])
        end
        #= none:13 =#
        return r
    end)
chevieset(:I, :ClassInfo, function (m,)
        #= none:7 =#
        local r, i, clnp, cl, g1, g2, gen, perm, m1
        #= none:9 =#
        r = (chevieget(:I, :WordsClassRepresentatives))(m)
        #= none:11 =#
        clnp = map(IntListToString, r)
        #= none:13 =#
        g1 = Perm()
        #= none:14 =#
        i = 2
        #= none:15 =#
        while 2i <= m + 1
            #= none:15 =#
            g1 = g1 * Perm(i, (m - i) + 2)
            #= none:16 =#
            i = i + 1
        end
        #= none:18 =#
        g2 = Perm()
        #= none:19 =#
        i = 1
        #= none:20 =#
        while 2i <= m
            #= none:20 =#
            g2 = g2 * Perm(i, (m - i) + 1)
            #= none:21 =#
            i = i + 1
        end
        #= none:23 =#
        gen = [g1, g2]
        #= none:25 =#
        perm = function (l,)
                #= none:25 =#
                if length(l) == 0
                    #= none:25 =#
                    return Perm()
                else
                    #= none:26 =#
                    return Product(gen[l])
                end
            end
        #= none:29 =#
        m1 = div(m, 2)
        #= none:31 =#
        if mod(m, 2) == 0
            #= none:31 =#
            cl = [1, m1, m1]
            #= none:32 =#
            cl = Append(cl, fill(0, max(0, (1 + (m1 - 1)) - 1)) + 2)
            #= none:33 =#
            push!(cl, 1)
        else
            #= none:35 =#
            cl = [1, m]
            #= none:36 =#
            cl = Append(cl, fill(0, max(0, (1 + m1) - 1)) + 2)
        end
        #= none:38 =#
        return Dict{Symbol, Any}(:classtext => r, :classnames => clnp, :classparams => clnp, :orders => map((i->begin
                                #= none:39 =#
                                order(perm(i))
                            end), r), :classes => cl)
    end)
chevieset(:I, :HeckeCharTable, function (m, param, rootparam)
        #= none:10 =#
        local u, v, squv, cl, r, ct, tbl
        #= none:12 =#
        u = -((param[1])[1]) // (param[1])[2]
        #= none:13 =#
        v = -((param[2])[1]) // (param[2])[2]
        #= none:15 =#
        if mod(m, 2) != 0
            #= none:15 =#
            squv = u
        elseif #= none:17 =# rootparam[1] !== nothing && rootparam[2] !== nothing
            #= none:18 =#
            squv = rootparam[1] * rootparam[2]
        else
            #= none:20 =#
            squv = GetRoot(u * v, 2, "CharTable(Hecke(I2(", m, ")))")
        end
        #= none:22 =#
        ct = [[u, v]]
        #= none:23 =#
        if mod(m, 2) == 0
            #= none:23 =#
            ct = Append(ct, [[u, -(u ^ 0)], [-(v ^ 0), v]])
        end
        #= none:25 =#
        push!(ct, [-(v ^ 0), -(v ^ 0)])
        #= none:27 =#
        cl = (chevieget(:I, :ClassInfo))(m)
        #= none:28 =#
        r = cl[:classtext]
        #= none:30 =#
        ct = map((i->begin
                        #= none:30 =#
                        map((x->begin
                                    #= none:30 =#
                                    Product(i[x])
                                end), r)
                    end), ct)
        #= none:32 =#
        ct = Append(ct, map((j->begin
                            #= none:32 =#
                            map(function (i,)
                                    #= none:32 =#
                                    local k
                                    #= none:34 =#
                                    k = length(r[i]) // 2
                                    #= none:36 =#
                                    if r[i] == []
                                        #= none:36 =#
                                        return 2 * v ^ 0
                                    elseif #= none:38 =# r[i] == [1]
                                        #= none:38 =#
                                        return u - 1
                                    elseif #= none:40 =# r[i] == [2]
                                        #= none:40 =#
                                        return v - 1
                                    else
                                        #= none:42 =#
                                        return squv ^ k * (E(m, k * j) + E(m, -k * j))
                                    end
                                end, 1:length(r))
                        end), 1:div(m - 1, 2)))
        #= none:46 =#
        tbl = Dict{Symbol, Any}(:identifier => SPrint("H(I2(", m, "))"), :cartan => CartanMat("I", 2, m), :size => 2m, :irredinfo => map((x->begin
                                #= none:48 =#
                                Dict{Symbol, Any}(:charparam => x, :charname => (chevieget(:I, :CharName))(m, x, Dict{Symbol, Any}(:TeX => true)))
                            end), ((chevieget(:I, :CharInfo))(m))[:charparams]), :parameter => [u, v], :powermap => [], :irreducibles => ct * v ^ 0)
        #= none:52 =#
        Inherit(tbl, cl)
        #= none:54 =#
        tbl[:centralizers] = map((i->begin
                        #= none:54 =#
                        tbl[:size] // i
                    end), tbl[:classes])
        #= none:56 =#
        tbl = ((CHEVIE[:compat])[:MakeCharacterTable])(tbl)
        #= none:58 =#
        ((CHEVIE[:compat])[:AdjustHeckeCharTable])(tbl, param)
        #= none:60 =#
        return tbl
    end)
chevieset(:I, :Representation, function (m, i)
        #= none:3 =#
        return (chevieget(:I, :HeckeRepresentation))(m, [[1, -1], [1, -1]], [1, 1], i)
    end)
chevieset(:I, :HeckeRepresentation, function (m, param, rootparam, i)
        #= none:4 =#
        local u, v, squv
        #= none:6 =#
        if i == 1
            #= none:6 =#
            return [[[(param[1])[1]]], [[(param[2])[1]]]]
        end
        #= none:8 =#
        if mod(m, 2) == 0
            #= none:8 =#
            i = i - 2
        end
        #= none:10 =#
        if i == 0
            #= none:10 =#
            return [[[(param[1])[1]]], [[(param[2])[2]]]]
        elseif #= none:12 =# i == 1
            #= none:12 =#
            return [[[(param[1])[2]]], [[(param[2])[1]]]]
        elseif #= none:14 =# i == 2
            #= none:14 =#
            return [[[(param[1])[2]]], [[(param[2])[2]]]]
        else
            #= none:17 =#
            u = -((param[1])[1]) // (param[1])[2]
            #= none:18 =#
            v = -((param[2])[1]) // (param[2])[2]
            #= none:20 =#
            if mod(m, 2) != 0
                #= none:20 =#
                squv = u
            elseif #= none:22 =# rootparam[1] !== nothing && rootparam[2] !== nothing
                #= none:23 =#
                squv = rootparam[1] * rootparam[2]
            else
                #= none:25 =#
                squv = GetRoot(u * v, 2, "Representation(Hecke(I2(", m, ")),[", i, "])")
            end
            #= none:27 =#
            return [-([[-(u ^ 0), u ^ 0], [0u, u]]) * (param[1])[2], -([[v, 0v], [u + v + squv * (E(m, i - 2) + E(m, 2 - i)), -(v ^ 0)]]) * (param[2])[2]]
        end
    end)
chevieset(:I, :Frobenius, function (m, sqrtu, j)
        #= none:4 =#
        return [[0, (1 // sqrtu) // (E(2m, j) + E(2m, -j))], [sqrtu * (E(2m, j) + E(2m, -j)), 0]] * sqrtu ^ 0
    end)
chevieset(:I, :PoincarePolynomial, function (m, param)
        #= none:3 =#
        local u, v
        #= none:5 =#
        u = -((param[1])[1]) // (param[1])[2]
        #= none:6 =#
        v = -((param[2])[1]) // (param[2])[2]
        #= none:8 =#
        if IsInt(m // 2)
            #= none:8 =#
            return Sum(1:m // 2, (i->begin
                                #= none:8 =#
                                (u * v) ^ (i - 1)
                            end)) * (u + 1) * (v + 1)
        else
            #= none:10 =#
            return Sum(1:m, (i->begin
                                #= none:10 =#
                                u ^ (i - 1)
                            end)) * (u + 1)
        end
    end)
chevieset(:I, :SchurElement, function (m, phi, para, rootpara)
        #= none:16 =#
        local u, v, ruv, e, ci
        #= none:18 =#
        if mod(m, 2) == 1
            #= none:19 =#
            ci = (chevieget(:I, :CharInfo))(m)
            #= none:21 =#
            ci = (ci[:malleParams])[Position(ci[:charparams], phi)]
            #= none:23 =#
            return (chevieget(:imp, :SchurElement))(m, 1, 2, ci, [map((i->begin
                                        #= none:24 =#
                                        E(m, i)
                                    end), 0:m - 1), para[2]], []) // m
        end
        #= none:26 =#
        u = -((para[1])[1]) // (para[1])[2]
        #= none:27 =#
        v = -((para[2])[1]) // (para[2])[2]
        #= none:29 =#
        if phi[1] == 1
            #= none:30 =#
            if phi[2] == m // 2
                #= none:30 =#
                e = (Sum(0:m // 2 - 1, (i->begin
                                        #= none:30 =#
                                        (u // v) ^ i
                                    end)) * (u + 1) * (v + 1)) // v
                #= none:32 =#
                if phi[3] == "'"
                    #= none:32 =#
                    return e
                else
                    #= none:33 =#
                    return (v // u) ^ (m // 2) * e
                end
            else
                #= none:35 =#
                e = Sum(0:m // 2 - 1, (i->begin
                                    #= none:35 =#
                                    (u * v) ^ i
                                end)) * (u + 1) * (v + 1)
                #= none:37 =#
                if phi[2] == 0
                    #= none:37 =#
                    return e
                else
                    #= none:38 =#
                    return (u * v) ^ (-m // 2) * e
                end
            end
        else
            #= none:41 =#
            e = E(m, phi[2]) + E(m, -(phi[2]))
            #= none:43 =#
            if all((i->begin
                            #= none:43 =#
                            rootpara[i] !== nothing
                        end), [1, 2])
                #= none:43 =#
                ruv = Product(rootpara)
            else
                #= none:45 =#
                ruv = GetRoot(u * v, 2, "SchurElement(Hecke(I2(", m, "),", phi, "))")
            end
            #= none:47 =#
            return (-m * ((u * v + 1) - ruv * e) * (u + v + e * ruv)) // (u * v * (e ^ 2 - 4))
        end
    end)
chevieset(:I, :FakeDegree, function (m, phi, q)
        #= none:4 =#
        if phi[1] == 1
            #= none:4 =#
            return q ^ phi[2]
        else
            #= none:5 =#
            return q ^ phi[2] + q ^ (m - phi[2])
        end
    end)
chevieset(:I, :CharTable, function (m,)
        #= none:3 =#
        local res
        #= none:5 =#
        res = (chevieget(:I, :HeckeCharTable))(m, [[1, -1], [1, -1]], [1, 1])
        #= none:7 =#
        res[:identifier] = SPrint("W(I2(", m, "))")
        #= none:9 =#
        return res
    end)
chevieset(:I, :DecompositionMatrix, function (n, p)
        #= none:3 =#
        local T, m
        #= none:5 =#
        T = (chevieget(:I, :CharTable))(n)
        #= none:6 =#
        T[:name] = T[:identifier]
        #= none:8 =#
        m = DecompositionMatrix(mod(T, p))
        #= none:10 =#
        return map((c->begin
                        #= none:10 =#
                        [c[1], (m[c[1]])[c[2]]]
                    end), BlocksMat(m))
    end)
chevieset(:I, :FactorizedSchurElement, function (arg...,)
        #= none:3 =#
        local ci
        #= none:5 =#
        if mod(arg[1], 2) == 0 && (arg[3])[1] != (arg[3])[2]
            #= none:5 =#
            error(" !  implemented")
        end
        #= none:7 =#
        ci = (chevieget(:I, :CharInfo))(arg[1])
        #= none:9 =#
        ci = (ci[:malleParams])[Position(ci[:charparams], arg[2])]
        #= none:11 =#
        return (chevieget(:imp, :FactorizedSchurElement))(arg[1], arg[1], 2, ci, arg[3], 1)
    end)
chevieset(:I, :Invariants, function (arg...,)
        #= none:3 =#
        local e, type_, m
        #= none:5 =#
        e = arg[1]
        #= none:7 =#
        if length(arg) == 2
            #= none:7 =#
            type_ = arg[2]
        elseif #= none:9 =# mod(e, 2) == 0
            #= none:9 =#
            type_ = 1
        else
            #= none:11 =#
            type_ = -(E(e, (e + 1) // 2)) - E(e, (e + 3) // 2)
        end
        #= none:13 =#
        m = DiagonalMat(1 + E(e, -1), -type_) * (chevieget(:imp, :GeneratingRoots))(e, e, 2)
        #= none:16 =#
        return map((f->begin
                        #= none:16 =#
                        function (arg...,)
                            #= none:17 =#
                            return ApplyFunc(f, arg * m)
                        end
                    end), (chevieget(:imp, :Invariants))(e, e, 2))
    end)
chevieset(:I, :SymbolToParameter, function (S,)
        #= none:5 =#
        if S[1] != [0, 1] || !([]) in S
            #= none:5 =#
            return false
        end
        #= none:7 =#
        if mod(length(S), 2) == 1
            #= none:7 =#
            S = reverse(S)
            #= none:9 =#
            return [Position(S, []), Position(S, [0, 1]) - Position(S, [])]
        else
            #= none:11 =#
            return (Position(S, []) + [-(Position(S[2:length(S)], [0, 1])), 0]) - 1
        end
    end)
chevieset(:I, :ParameterToSymbol, function (e, p)
        #= none:5 =#
        local S
        #= none:7 =#
        if p == [0]
            #= none:7 =#
            S = map((x->begin
                            #= none:7 =#
                            [0]
                        end), 1:e)
            #= none:8 =#
            S[e] = [2]
        elseif #= none:10 =# p == [1]
            #= none:10 =#
            S = map((x->begin
                            #= none:10 =#
                            [0, 1]
                        end), 1:e)
            #= none:11 =#
            S[e] = [1, 2]
        elseif #= none:13 =# length(p) == 3
            #= none:13 =#
            S = map((x->begin
                            #= none:13 =#
                            [0]
                        end), 1:e // 2 - 1)
            #= none:14 =#
            S = Append(S, [[1], 2, (p[3] + 1) // 2])
        elseif #= none:16 =# mod(e, 2) == 0
            #= none:16 =#
            S = map((x->begin
                            #= none:16 =#
                            [0]
                        end), 1:e)
            #= none:18 =#
            if p[1] == 0
                #= none:18 =#
                S[[e, e - p[2]]] = [[1], [1]]
            else
                #= none:20 =#
                S[1 + [0, mod(p[2] - p[1], e)]] = [[0, 1], [0, 1]]
                #= none:22 =#
                S[1 + [mod(-(p[1]), e), p[2]]] = [[], []]
            end
        else
            #= none:24 =#
            S = map((i->begin
                            #= none:24 =#
                            [0]
                        end), 1:e)
            #= none:26 =#
            if p[1] != 0
                #= none:26 =#
                S[1 + [0, mod(-(Sum(p)), e)]] = [[0, 1], [0, 1]]
                #= none:28 =#
                S[1 + map((x->(#= none:28 =#
                                        mod(x, e))), -p)] = [[], []]
            else
                #= none:30 =#
                S[e + [-(mod(p[2] - p[1], e)), 0]] = [[1], [1]]
            end
        end
        #= none:33 =#
        return S
    end)
chevieset(:I, :UnipotentCharacters, function (e,)
        #= none:3 =#
        local cusp, uc, f
        #= none:5 =#
        f = div(e, 2)
        #= none:6 =#
        uc = Dict{Symbol, Any}()
        #= none:8 =#
        uc[:harishChandra] = [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "I", :indices => [1, 2], :rank => 2, :bond => e), :parameterExponents => [1, 1], :levi => [], :eigenvalue => 1, :cuspidalName => "")]
        #= none:13 =#
        if mod(e, 2) != 0
            #= none:13 =#
            ((uc[:harishChandra])[1])[:charNumbers] = 1:f + 2
        else
            #= none:15 =#
            ((uc[:harishChandra])[1])[:charNumbers] = Concatenation([1, 3, 4, 2], 4 + (1:f - 1))
        end
        #= none:19 =#
        cusp = Concatenation(map((k->begin
                            #= none:19 =#
                            map((l->begin
                                        #= none:19 =#
                                        [k, l]
                                    end), k + 1:(e - k) - 1)
                        end), 1:f - 1))
        #= none:21 =#
        f = (f + 1) - mod(e, 2)
        #= none:23 =#
        uc[:harishChandra] = Append(uc[:harishChandra], map((x->begin
                            #= none:23 =#
                            Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :parameterExponents => [], :levi => [1, 2], :eigenvalue => E(e, -(Product(cusp[x]))), :cuspidalName => SPrint("I_2(", e, ")", FormatGAP(cusp[x])), :charNumbers => [x + f + 2])
                        end), 1:length(cusp)))
        #= none:30 =#
        uc[:families] = [Family(((CHEVIE[:families])[:Dihedral])(e), (1:length(cusp) + f) + 2), Family("C1", [1]), Family("C1", [2])]
        #= none:33 =#
        uc[:parameters] = Concatenation([[0], [1]], ((uc[:families])[1])[:parameters])
        #= none:35 =#
        uc[:charSymbols] = map((p->begin
                        #= none:35 =#
                        (chevieget(:I, :ParameterToSymbol))(e, p)
                    end), uc[:parameters])
        #= none:38 =#
        uc[:a] = Concatenation([0, e], map((x->begin
                            #= none:38 =#
                            1
                        end), ((uc[:families])[1])[:parameters]))
        #= none:40 =#
        uc[:A] = Concatenation([0, e], map((x->begin
                            #= none:40 =#
                            e - 1
                        end), ((uc[:families])[1])[:parameters]))
        #= none:42 =#
        if e == 5
            #= none:42 =#
            uc[:curtis] = [2, 1, 3, 4, 6, 5]
        end
        #= none:44 =#
        return uc
    end)
