chevieset(:imp, :SemisimpleRank, function (p, q, r)
        return r
    end)
chevieset(:imp, :BraidRelations, function (p, q, r)
        local i, b, res
        b(i, j, o) = begin
                local p
                p(i, j) = begin
                        return map((k->begin
                                        i * mod(k, 2) + j * mod(1 - k, 2)
                                    end), 1:o)
                    end
                return [p(i, j), p(j, i)]
            end
        res = []
        if q == 1
            if r >= 2
                if p == 1
                    push!(res, b(1, 2, 3))
                else
                    push!(res, b(1, 2, 4))
                end
            end
            res = Append(res, map((i->begin
                                b(i, i - 1, 3)
                            end), 3:r))
            for i = 3:r
                res = Append(res, map((j->begin
                                    b(i, j, 2)
                                end), 1:i - 2))
            end
        elseif p == q
            push!(res, b(1, 2, p))
            if r >= 3
                res = Append(res, [[[1, 2, 3, 1, 2, 3], [3, 1, 2, 3, 1, 2]], b(1, 3, 3), b(2, 3, 3)])
            end
            res = Append(res, map((i->begin
                                b(i, i - 1, 3)
                            end), 4:r))
            for i = 4:r
                res = Append(res, map((j->begin
                                    b(i, j, 2)
                                end), 1:i - 2))
            end
        else
            push!(res, [[1, 2, 3], [2, 3, 1]])
            i = b(2, 3, q - 1)
            push!(res, [Concatenation([1, 2], i[2]), Concatenation([3, 1], i[1])])
            if r >= 3
                if q != 2
                    push!(res, [[2, 3, 4, 2, 3, 4], [4, 2, 3, 4, 2, 3]])
                end
                res = Append(res, [b(2, 4, 3), b(3, 4, 3), b(1, 4, 2)])
            end
            res = Append(res, map((i->begin
                                b(i, i - 1, 3)
                            end), 5:r + 1))
            for i = 5:r + 1
                res = Append(res, map((j->begin
                                    b(i, j, 2)
                                end), 1:i - 2))
            end
        end
        return res
    end)
chevieset(:imp, :Size, function (p, q, r)
        return div(p ^ r * factorial(r), q)
    end)
chevieset(:imp, :ReflectionName, function (arg...,)
        local n, option
        option = arg[4]
        if arg[3] == 1 && arg[2] == 1
            if haskey(option, :TeX)
                return SPrint("Z_{", arg[1], "}")
            else
                return SPrint("Z", arg[1])
            end
        end
        if haskey(option, :TeX)
            n = SPrint("G_{", Join(arg[1:3]), "}")
        else
            n = SPrint("G", joindigits(arg[1:3]))
        end
        if length(arg) == 5
            n *= SPrint("(", Format(arg[4], option), ")")
        end
        return n
    end)
chevieset(:imp, :CartanMat, function (p, q, r)
        local rt, rbar, e
        rt = (chevieget(:imp, :GeneratingRoots))(p, q, r)
        rbar = conj(rt)
        e = (chevieget(:imp, :EigenvaluesGeneratingReflections))(p, q, r)
        e = 1 - map((x->begin
                            E(denominator(x), numerator(x))
                        end), e)
        e = map((i->begin
                        (e[i] * rbar[i]) // (rbar[i] * rt[i])
                    end), 1:length(e))
        return map((x->begin
                        map((y->begin
                                    x * y
                                end), rt)
                    end), e)
    end)
chevieset(:imp, :ReflectionDegrees, function (p, q, r)
        return Concatenation(p * (1:r - 1), [div(r * p, q)])
    end)
chevieset(:imp, :ParabolicRepresentatives, function (p, q, r, s)
        local t
        if q == 1
            if p == 1
                if s == 0
                    return [[]]
                end
                return map((j->begin
                                Concatenation(map((k->begin
                                                ((Sum(j[1:k - 1]) + k) - 1) + (1:j[k])
                                            end), 1:length(j)))
                            end), Concatenation(map((i->begin
                                        partitions(s, i)
                                    end), 1:(r + 1) - s)))
            else
                return Concatenation(map((i->begin
                                    map((j->begin
                                                Concatenation(1:i, i + 1 + j)
                                            end), (chevieget(:imp, :ParabolicRepresentatives))(1, 1, (r - i) - 1, s - i))
                                end), 0:s))
            end
        elseif r == 2
            if q == 2
                t = [[[]], [[1], [2], [3]], [1:3]]
                return t[s + 1]
            elseif p == q
                if mod(p, 2) == 0
                    t = [[[]], [[1], [2]], [[1, 2]]]
                    return t[s + 1]
                else
                    t = [[[]], [[1]], [[1, 2]]]
                    return t[s + 1]
                end
            else
                return false
            end
        else
            return false
        end
    end)
chevieset(:imp, :NrConjugacyClasses, function (p, q, r)
        if [q, r] == [2, 2]
            return div(p * (p + 6), 4)
        elseif q == 1
            return npartition_tuples(r, p)
        else
            return length(((chevieget(:imp, :ClassInfo))(p, q, r))[:classtext])
        end
    end)
chevieset(:imp, :CharParams, function (de, e, r)
        local t, tt, s, charparams, d
        if e == 1
            return partition_tuples(r, de)
        else
            charparams = []
            d = div(de, e)
            for t = partition_tuples(r, de)
                tt = map((i->begin
                                circshift(t, i)
                            end), (1:e) * d)
                if t == Minimum(tt)
                    s = Position(tt, t)
                    if s == e
                        push!(charparams, t)
                    else
                        t = t[1:s * d]
                        s = div(e, s)
                        charparams = Append(charparams, map((i->begin
                                            Concatenation(t, [s, i])
                                        end), 0:s - 1))
                    end
                end
            end
            return charparams
        end
    end)
chevieset(:imp, :LowestPowerFakeDegrees, function (p, q, r)
        local ci
        if q == 1 || p == q
            error("should not be called")
        end
        return false
    end)
chevieset(:imp, :HighestPowerFakeDegrees, function (p, q, r)
        local ci
        if q == 1 || p == q
            error("should not be called")
        end
        return false
    end)
chevieset(:imp, :CharSymbols, function (p, q, r)
        local s, ss, res
        if q == 1
            return SymbolsDefect(p, r, 0, 1)
        elseif q == p
            return SymbolsDefect(p, r, 0, 0)
        else
            return false
        end
    end)
chevieset(:imp, :FakeDegree, function (p, q, r, c, v)
        if q == 1
            c = fegsymbol(symbol_partition_tuple(c, 1))
        elseif q == p
            c = fegsymbol(symbol_partition_tuple(c, fill(0, max(0, (1 + p) - 1))))
        else
            return false
        end
        return Value(c, v)
    end)
chevieset(:imp, :SchurModel, function (p, q, r, phi)
        local l, i, j, res, s, t, ci, GenHooks, v, h, d, p2
        if q == 1
            GenHooks(l, m) = begin
                    if length(l) == 0
                        return []
                    end
                    m = conjugate_partition(m)
                    m = Append(m, fill(0, max(0, (1 + (l[1] - length(m))) - 1)))
                    m = (1 + m) - (1:length(m))
                    return Concatenation(map((i->begin
                                        (l[i] - i) + m[1:l[i]]
                                    end), 1:length(l)))
                end
            res = Dict{Symbol, Any}(:coeff => (-1) ^ (r * (p - 1)), :factor => fill(0, max(0, (1 + p) - 1)), :vcyc => [])
            l = Concatenation(phi)
            sort!(l)
            push!(res[:factor], ((1:length(l)) - length(l)) * l)
            for s = 1:p
                for t = 1:p
                    for h = GenHooks(phi[s], phi[t])
                        v = fill(0, max(0, (1 + p) - 1))
                        if s != t
                            v[[s, t]] = [1, -1]
                            push!(v, h)
                            push!(res[:vcyc], [v, 1])
                        else
                            push!(v, 1)
                            for d = divisors(h)
                                if d > 1
                                    push!(res[:vcyc], [v, d])
                                end
                            end
                        end
                    end
                end
            end
            return res
        elseif [q, r] == [2, 2]
            ci = (chevieget(:imp, :CharInfo))(p, q, r)
            phi = (ci[:malle])[Position(ci[:charparams], phi)]
            p2 = div(p, 2)
            if phi[1] == 1
                res = Dict{Symbol, Any}(:coeff => 1, :factor => fill(0, max(0, (1 + (4 + p2)) - 1)), :vcyc => [])
                for l = [[1, -1, 0, 0], [0, 0, 1, -1]]
                    l = Append(l, fill(0, max(0, (1 + p2) - 1)))
                    push!(res[:vcyc], [l, 1])
                end
                for i = 2:p2
                    for l = [[0, 0, 0, 0, 1], [1, -1, 1, -1, 1]]
                        l = Append(l, fill(0, max(0, (1 + (p2 - 1)) - 1)))
                        l[4 + i] = -1
                        push!(res[:vcyc], [l, 1])
                    end
                end
            else
                res = Dict{Symbol, Any}(:coeff => -2, :factor => fill(0, max(0, (1 + (4 + p2)) - 1)), :vcyc => [], :root => fill(0, max(0, (1 + (4 + p2)) - 1)) // 1)
                res[:rootCoeff] = E(p2, (2 - phi[3]) - phi[4])
                (res[:root])[1:6] = [1, 1, 1, 1, 1, 1] // 2
                for i = 3:p2
                    for j = [1, 2]
                        l = fill(0, max(0, (1 + (4 + p2)) - 1))
                        l[4 + [j, i]] = [1, -1]
                        push!(res[:vcyc], [l, 1])
                    end
                end
                if haskey(CHEVIE, :old)
                    for l = [[0, -1, 0, -1, -1, 0], [0, -1, -1, 0, -1, 0], [-1, 0, -1, 0, -1, 0], [-1, 0, 0, -1, -1, 0]]
                        l = Append(l, fill(0, max(0, (1 + (p2 - 2)) - 1)))
                        push!(l, 1)
                        push!(res[:vcyc], [l, 1])
                    end
                else
                    for l = [[0, -1, 0, -1, -1, 0], [0, -1, -1, 0, 0, -1], [-1, 0, -1, 0, -1, 0], [-1, 0, 0, -1, 0, -1]]
                        l = Append(l, fill(0, max(0, (1 + (p2 - 2)) - 1)))
                        push!(l, 1)
                        push!(res[:vcyc], [l, 1])
                    end
                end
            end
            return res
        else
            error("not implemented")
        end
    end)
chevieset(:imp, :SchurData, function (p, q, r, phi)
        local ci, res
        if [q, r] == [2, 2]
            ci = (chevieget(:imp, :CharInfo))(p, q, r)
            phi = (ci[:malle])[Position(ci[:charparams], phi)]
            if phi[1] == 1
                res = Dict{Symbol, Any}(:order => [phi[2], 3 - phi[2], 2 + phi[3], 5 - phi[3], 4 + phi[4]])
                res[:order] = Append(res[:order], 4 + Difference(1:div(p, 2), [phi[4]]))
                return res
            else
                res = Dict{Symbol, Any}(:order => [1, 2, 3, 4, 4 + phi[3], 4 + phi[4]])
                res[:order] = Append(res[:order], 4 + Difference(1:div(p, 2), phi[[3, 4]]))
                res[:rootPower] = phi[2] * E(p, (phi[3] + phi[4]) - 2)
                return res
            end
        else
            error("not implemented")
        end
    end)
chevieset(:imp, :SchurElement, function (p, q, r, phi, para, root)
        local m, Z, e1
        if r == 1
            return VcycSchurElement(Concatenation(para[1], [0]), ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi))
        elseif p == 1
            return VcycSchurElement([0, -((para[1])[1]) // (para[1])[2]], ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi))
        elseif q == 1
            return VcycSchurElement(Concatenation(para[1], [-((para[2])[1]) // (para[2])[2]]), ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi))
        elseif r == 2 && mod(q, 2) == 0
            e1 = div(q, 2)
            Z = map((x->begin
                            GetRoot(x, e1)
                        end), para[1])
            Z = Concatenation(map((j->begin
                                Z * E(e1, j)
                            end), 0:e1 - 1))
            return VcycSchurElement(Concatenation(para[2], para[3], Z), ((CHEVIE[:imp])[:SchurModel])(p, 2, r, phi), ((CHEVIE[:imp])[:SchurData])(p, 2, r, phi)) // e1
        elseif p == q
            if IsInt(phi[length(phi)])
                m = length(phi) - 2
                phi = fullsymbol(phi)
            else
                m = p
            end
            return ((CHEVIE[:imp])[:SchurElement])(p, 1, r, phi, Concatenation([map((i->begin
                                            E(p, i)
                                        end), 0:p - 1)], para[2:length(para)]), []) // m
        elseif para[2] == para[3]
            if IsInt(phi[length(phi)])
                m = length(phi) - 2
                phi = fullsymbol(phi)
            else
                m = p
            end
            if para[1] == map((i->begin
                                E(p // q, i - 1)
                            end), 1:p // q)
                para = [map((i->begin
                                    E(p, i)
                                end), 0:p - 1), para[2]]
            else
                para = [Concatenation(TransposedMat(map((i->begin
                                            map((j->begin
                                                            E(q, j)
                                                        end), 0:q - 1) * GetRoot(i, q)
                                        end), para[1]))), para[2]]
            end
            return (p // q * ((CHEVIE[:imp])[:SchurElement])(p, 1, r, phi, para, [])) // m
        else
            ((CHEVIE[:compat])[:InfoChevie])("# SchurElements(H(G(", p, ",", q, ",", r, "),", para, ") not implemented\n")
            return false
        end
    end)
chevieset(:imp, :FactorizedSchurElement, function (p, q, r, phi, para, root)
        local m, F
        if r == 1
            return VFactorSchurElement(Concatenation(para[1], [0]), ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi))
        elseif p == 1
            return VFactorSchurElement([0, -((para[1])[1]) // (para[1])[2]], ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi))
        elseif q == 1
            return VFactorSchurElement(Concatenation(para[1], [-((para[2])[1]) // (para[2])[2]]), ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi))
        elseif [q, r] == [2, 2]
            return VFactorSchurElement(Concatenation(para[[2, 3, 1]]), ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi), ((CHEVIE[:imp])[:SchurData])(p, q, r, phi))
        elseif p == q
            if IsInt(phi[length(phi)])
                m = length(phi) - 2
                phi = fullsymbol(phi)
            else
                m = p
            end
            if para[1] != para[2]
                InfoChevie("# FactorizedSchurElements(H(G(", p, ",", q, ",", r, "),", para, ") not implemented\n")
                return false
            end
            F = ((CHEVIE[:imp])[:FactorizedSchurElement])(p, 1, r, phi, Concatenation([map((i->begin
                                        E(p, i)
                                    end), 0:p - 1)], para[2:length(para)]), [])
            F[:factor] = F[:factor] // m
            return F
        elseif para[2] == para[3]
            if IsInt(phi[length(phi)])
                m = length(phi) - 2
                phi = fullsymbol(phi)
            else
                m = p
            end
            if para[1] == map((i->begin
                                E(p // q, i - 1)
                            end), 1:p // q)
                para = [map((i->begin
                                    E(p, i)
                                end), 0:p - 1), para[2]]
            else
                para = [Concatenation(TransposedMat(map((i->begin
                                            map((j->begin
                                                            E(q, j)
                                                        end), 0:q - 1) * GetRoot(i, q)
                                        end), para[1]))), para[2]]
            end
            F = ((CHEVIE[:imp])[:FactorizedSchurElement])(p, 1, r, phi, para, [])
            F[:factor] = p // (q * m) * F[:factor]
            return F
        else
            InfoChevie("# FactorizedSchurElements(H(G(", p, ",", q, ",", r, "),", para, ") not implemented\n")
            return false
        end
    end)
chevieset(:imp, :HeckeRepresentation, function (p, q, r, para, rootpara, i)
        local X, Y, t, x, a, v, d, T, S, m334, m335, extra, l, m, m1, p1rRep, f, e
        if !(IsList(para))
            para = [para]
        end
        if [q, r] == [1, 2]
            X = para[2]
            Y = para[1]
            t = (partition_tuples(2, p))[i]
            if count((x->begin
                                x != []
                            end), t) == 1
                p = PositionProperty(t, (x->begin
                                x != []
                            end))
                if t[p] == [2]
                    return X[1] ^ 0 * [[[Y[p]]], [[X[1]]]]
                else
                    return X[1] ^ 0 * [[[Y[p]]], [[X[2]]]]
                end
            else
                p = filter((i->begin
                                t[i] != []
                            end), 1:length(t))
                return X[1] ^ 0 * [[[Y[p[1]], 0], [-1, Y[p[2]]]], [[X[1], X[1] * Y[p[1]] + X[2] * Y[p[2]]], [0, X[2]]]]
            end
        elseif [p, q, r] == [3, 3, 3]
            x = -((para[2])[1]) // (para[2])[2]
            f(x, j) = begin
                    return [[[-1, 0, 0], [0, 0, 1], [0, x, -1 + x]], [[-1, 0, 0], [x - x ^ 2, -1 + x, j ^ 2], [j * x - j * x ^ 2, j * x, 0]], [[0, 1, 0], [x, -1 + x, 0], [0, 0, -1]]]
                end
            r = x ^ 0 * [[[[-1, 0], [-1, x]], [[x, -x], [0, -1]], [[x, -x], [0, -1]]], [[[-1, 0], [-1, x]], [[x, -x], [0, -1]], [[-1, 0], [-1, x]]], [[[-1, 0], [-1, x]], [[x, -x], [0, -1]], [[-1 + x, 1], [x, 0]]], f(x, E(3)), f(x, E(3, 2)), [[[-1]], [[-1]], [[-1]]], -x * f(x ^ -1, E(3, 2)), -x * f(x ^ -1, E(3)), [[[-1, 0], [-1, x]], [[-1, 0], [-1, x]], [[x, -x], [0, -1]]], [[[x]], [[x]], [[x]]]]
            return -((para[2])[2]) * r[i]
        elseif [p, q, r] == [2, 2, 4]
            x = -((para[1])[1]) // (para[1])[2]
            r = [(x->begin
                            [[[-1 + x, -1, 0], [-x, 0, 0], [x - x ^ 2, -1 + x, -1]], [[0, 1, 0], [x, -1 + x, 0], [0, 0, -1]], [[-1, 0, 0], [0, 0, 1], [0, x, -1 + x]], [[0, 1, 0], [x, -1 + x, 0], [0, 0, -1]]]
                        end), (x->begin
                            [[[0, 1, 0], [x, -1 + x, 0], [0, 0, -1]], [[-1 + x, -1, 0], [-x, 0, 0], [x - x ^ 2, -1 + x, -1]], [[-1, 0, 0], [0, 0, 1], [0, x, -1 + x]], [[0, 1, 0], [x, -1 + x, 0], [0, 0, -1]]]
                        end), (x->begin
                            [[[-1, 0, 0, 0], [0, -1 + x, -1, 0], [0, -x, 0, 0], [0, 0, 0, -1]], [[-1, 1 - x, 1 - x, 0], [0, 0, 1, 0], [0, x, -1 + x, 0], [0, -1 + x, -1 + x, -1]], [[-1 + x, -x, 0, 0], [-1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]], [[0, 0, 0, 1], [0, -1, 0, 0], [0, 0, -1, 0], [x, 0, 0, -1 + x]]]
                        end), (x->begin
                            [[[-1]], [[-1]], [[-1]], [[-1]]]
                        end), (x->begin
                            [[[x, 1 - x, -1 + x, -x + x ^ 2, x - x ^ 2, 0], [0, -1 + x, 0, 0, -x, x - x ^ 2], [0, 0, -1 + x, -x, 0, x - x ^ 2], [0, 0, -1, 0, 0, -1 + x], [0, -1, 0, 0, 0, -1 + x], [0, 0, 0, 0, 0, -1]], [[x, 0, 0, 0, 0, 0], [0, 0, 0, 0, x, 0], [0, 0, 0, x, 0, 0], [0, 0, 1, -1 + x, 0, 0], [0, 1, 0, 0, -1 + x, 0], [0, 0, 0, 0, 0, -1]], [[0, 0, x, 0, 0, 0], [0, -1, 0, 0, 0, 0], [1, 0, -1 + x, 0, 0, 0], [0, 0, 0, x, 0, 0], [0, 0, 0, 0, 0, x], [0, 0, 0, 0, 1, -1 + x]], [[-1, 0, 0, 0, 0, 0], [0, -1 + x, 1, 0, 0, 0], [0, x, 0, 0, 0, 0], [0, 0, 0, 0, x, 0], [0, 0, 0, 1, -1 + x, 0], [0, 0, 0, 0, 0, x]]]
                        end), (x->begin
                            [[[-1 + x, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0], [-x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0], [0, x, 0, 0, -1 + x, 0, 0, 0], [0, 0, 0, 0, 0, -1 + x, 0, x], [0, 0, 0, x, 0, 0, -1 + x, 0], [0, 0, 0, 0, 0, 1, 0, 0]], [[0, 0, 1, 0, 0, 0, 0, 0], [0, 0, -1 + x, 0, 1, 0, (-1 + x) // x, 0], [x, 0, -1 + x, 0, 0, 0, 0, 0], [0, 0, 0, -1 + x, 0, 0, -1, 0], [x - x ^ 2, x, 0, -1 + x, -1 + x, 0, ((1 - 2x) + x ^ 2) // x, 0], [-x + x ^ 2, 0, 0, -x + x ^ 2, 0, -1 + x, 1 - x, x], [0, 0, 0, -x, 0, 0, 0, 0], [0, 0, 1 - x, x - x ^ 2, 0, 1, -1 + x, 0]], [[0, 1, 0, 0, 0, 0, 0, 0], [x, -1 + x, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, x, -1 + x, 0, 0, 0, 0], [0, 0, 0, 0, -1 + x, 0, -1, 0], [0, 0, 0, 0, 0, x, 0, 0], [0, 0, 0, 0, -x, 0, 0, 0], [0, 0, 0, 0, 0, -x, 0, -1]], [[-1, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, -1, -x, 0, 0, 0, 0], [0, 0, 0, x, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, x, -1 + x, 0, 0], [0, 0, 0, 0, 0, 0, x, 0], [0, x, 0, 0, 0, 0, 0, -1 + x]]]
                        end), (x->begin
                            [[[-1, -1, 0], [0, x, 0], [0, 1, -1]], [[-1, -1, 0], [0, x, 0], [0, 1, -1]], [[-1 + x, x, 0], [1, 0, 0], [0, 0, -1]], [[0, 0, 1], [0, -1, 0], [x, 0, -1 + x]]]
                        end), 1, 2, (x->begin
                            [[[x, 0], [-1, -1]], [[x, 0], [-1, -1]], [[0, 1], [x, -1 + x]], [[x, 0], [-1, -1]]]
                        end), 3, 7, 4]
            if IsInt(r[i])
                return (para[1])[2] * x * (r[r[i]])(x ^ -1)
            else
                return -((para[1])[2]) * (r[i])(x) * x ^ 0
            end
        elseif [p, q, r] == [3, 3, 4]
            x = -((para[2])[1]) // (para[2])[2]
            m334(i) = begin
                    local f1, f2, f3, f5, f7, f8, f11, f13
                    f1 = (x->begin
                                x ^ 0 * [[[x, -1, 0, 0, 0, 0, 0, 0, 0, 0, ((1 - x) - x ^ 2) + x ^ 3, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1 + x, 0, x, 0, -x, 0, 0, 0, x - x ^ 2, 0], [0, 0, 0, -1 + x, 0, 0, -x, 0, 0, 0, x - x ^ 2, 0], [0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 1, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, 0, x, 0, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1 + x]], [[0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0], [0, 0, 0, 0, x, 0, 0, 1 - x, 0, 0, 0, 0], [0, 0, 0, 1, -1 + x, 0, 0, 1 - x, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1 - x, -1 + x, 1, 0, (1 - x) + x ^ 2], [0, 0, 1, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x], [0, 0, 0, 0, 0, x, 0, x - x ^ 2, -x, -1 + x, 0, x - x ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 + x]], [[0, (-1 + 2x) - x ^ 2, 1 - x, x, -x + x ^ 2, 0, 0, 0, (-1 + 2x) - x ^ 2, 0, 0, 0], [0, -1 + x, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0], [0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0], [1, -1 + x, 0, -1 + x, 0, 0, 1 - x, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1 + x, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, -1], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, -x, 0, 0, 0, 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0], [0, 0, 0, 0, 0, x, 0, 0, 0, 0, 1 - x, x - x ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, x], [0, 0, 0, 1, 0, -1 + x, -1 + x, 0, 0, 0, 1 - x, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x], [0, 0, 0, 0, 0, 0, 0, x, 0, 0, -1, 0], [0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0], [0, 0, 0, 0, x, 0, -x, 0, 0, -1 + x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + x]]]
                            end)
                    f2(x, j) = begin
                            return [[[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, x, x, x]], [[-1, 0, 0, 0], [0, -1, 0, 0], [0, -(j ^ 2), x, 1], [0, 0, 0, -1]], [[-1, 0, 0, 0], [x, x, -j * x, 1], [0, 0, -1, 0], [0, 0, 0, -1]], [[x, 1, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]]
                        end
                    f3 = (x->begin
                                [[[x, -1, 0, -x, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1 + x, 0, 1, x], [0, 0, 0, -1, 0, 0], [0, 0, x, 0, 0, x], [0, 0, 0, 0, 0, -1]], [[-1, 0, 0, 0, 0, 0], [0, 0, 0, 0, x, 1], [0, 0, -1, 0, 0, 0], [-1, 0, -1, x, 0, -1 + x], [0, 1, 0, 0, -1 + x, 1], [0, 0, 0, 0, 0, -1]], [[0, x, 1, -1, -1, 0], [1, -1 + x, 1, -1, -1, 0], [0, 0, -1, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, 1, 1, x]], [[x, -1, 0, 0, 1, x], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 0, -1, x, 1, x], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]]]
                            end)
                    f5 = (x->begin
                                [[[-1]], [[-1]], [[-1]], [[-1]]]
                            end)
                    f7(x, j) = begin
                            return [[[-1, 0, 0, 0, 0, 0], [x, x, 0, 0, 0, 0], [x, 0, x, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, -j * x ^ 2, x, x]], [[x, 1, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, j ^ 2, x, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, x, x, 1], [0, 0, 0, 0, 0, -1]], [[x, 0, 1, 0, 1, 0], [0, x, j * x, 0, 0, 1], [0, 0, -1, 0, 0, 0], [0, 0, 0, x, 1, -(j ^ 2) * x ^ -1], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]], [[-1, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 0, 0, x, 0, 0], [x, 0, 0, 0, x, 0], [0, x, 0, 0, 0, x]]]
                        end
                    f8(x, j) = begin
                            return [[[-1, 0, 0, 0, 0, 0, 0, 0], [1, x, 0, 0, 0, 0, 1, 0], [1, 0, x, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, -j * x, x, x, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, (j ^ 2 - j) * x, 0, 1, x]], [[x, x, 0, 0, 0, 0, -(j ^ 2), 0], [0, -1, 0, 0, 0, 0, 0, 0], [0, j, x, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 1, x, 1, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -2 * j ^ 2 - j, 1, x]], [[x, 0, x, 0, x, 0, 0, 0], [0, x, j ^ 2 * x, 0, 0, 1, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, x, x, -(j ^ 2), 0, 0], [0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, x, x], [0, 0, 0, 0, 0, 0, 0, -1]], [[-1, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, x, 0, 0, -(j ^ 2), 0], [1, 0, 0, 0, x, 0, 0, 0], [0, x, 0, 0, 0, x, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0], [0, 0, (j ^ 2 - j) * x, 0, 0, 0, 1, x]]]
                        end
                    f11 = (x->begin
                                [[[x, 1, 0], [0, -1, 0], [0, 0, -1]], [[x, 1, 0], [0, -1, 0], [0, 0, -1]], [[-1, 0, 0], [x, x, 1], [0, 0, -1]], [[-1, 0, 0], [0, -1, 0], [0, x, x]]]
                            end)
                    f13 = (x->begin
                                [[[-1, 0], [x, x]], [[-1, 0], [x, x]], [[x, 1], [0, -1]], [[-1, 0], [x, x]]]
                            end)
                    r = [f1(x), f2(x, E(3)), f3(x), f2(x, E(3, 2)), f5(x), -x * f1(x ^ -1), f7(x, E(3)), f8(x, E(3)), f8(x, E(3, 2)), -x * f7(x ^ -1, E(3)), f11(x), -x * f3(x ^ -1), f13(x), -x * f2(x ^ -1, E(3, 2)), -x * f2(x ^ -1, E(3)), -x * f11(x ^ -1), -x * f5(x ^ -1)]
                    return x ^ 0 * r[i]
                end
            return -((para[2])[2]) * m334(i)
        elseif [p, q, r] == [3, 3, 5]
            x = -((para[2])[1]) // (para[2])[2]
            m335(i) = begin
                    local r, f1, f2, f3, f4, f8, f9, f11, f12, f13, f17, f20, f23, f29
                    f1(x) = begin
                            return x ^ 0 * [[[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [E(3, 2) - E(3, 2) * x, E(3, 2) - E(3, 2) * x, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -(E(3, 2)) * x + E(3, 2) * x ^ 2, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(-(E(3, 2)) - root(-3) * x) + E(3) * x ^ 2, (E(3) - E(3) * x) + E(3) * x ^ 2, 0, 0, 0, E(3) - E(3) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E(3, 2) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, E(3), -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, (E(3, 2) * x + root(-3) * x ^ 2) - E(3) * x ^ 3, 0, (E(3) - E(3) * x) + E(3) * x ^ 2, 0, 0, 0, E(3) - E(3) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * x, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * x, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, (-(E(3, 2)) - root(-3) * x) + E(3) * x ^ 2, 0, 0, 0, 0, 0, 0, 0, (E(3) - E(3) * x) + E(3) * x ^ 2, 0, 0, 0, 0, E(3) - E(3) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, ((E(3, 2) - E(3, 2) * x) - E(3, 2) * x ^ 2) + E(3, 2) * x ^ 3, 0, (E(3) - root(-3) * x) - E(3, 2) * x ^ 2, 0, 0, 0, -(E(3, 2)) + E(3, 2) * x, 0, -(E(3)) + root(-3) * x + E(3, 2) * x ^ 2, 0, 0, 0, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0], [E(3) - E(3) * x, 0, 0, (2 * E(3) - E(3) * x ^ -1) - E(3) * x, 0, 0, 0, 0, 0, 0, (-(E(3)) + 2 * E(3) * x) - E(3) * x ^ 2, 0, 0, 0, 0, -(E(3)) + E(3) * x, 0, 0, 0, E(3) - E(3) * x, 0, 0, E(3) - E(3) * x, 0, E(3), 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * x, 0, -1 + x, 0, 0, 0, 0, 0, 0], [-(E(3)) + root(-3) * x + E(3, 2) * x ^ 2, 0, (x - 2 * x ^ 2) + x ^ 3, (1 - 3 * root(-3)) // 2 + E(3) * x ^ -1 + (1 + 3 * root(-3)) // 2 * x + E(3, 2) * x ^ 2, 0, 0, 0, 0, 0, 0, ((E(3) - 2 * E(3) * x) + 2 * E(3) * x ^ 2) - E(3) * x ^ 3, 0, 0, 0, 0, (E(3) - 2 * E(3) * x) + E(3) * x ^ 2, 0, -(E(3)) + E(3) * x, 0, (-(E(3)) + 2 * E(3) * x) - E(3) * x ^ 2, 0, 0, (E(3, 2) - E(3, 2) * x) + E(3, 2) * x ^ 2, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, -(E(3, 2)) + E(3, 2) * x, 0, 0, 0, -1 + x ^ -1, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), -1 + x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, E(3) * x], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, x, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 0, 0, 0, 0, 0, -1 + x, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + x ^ -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, x, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2x) + x ^ 2, (1 - 2x) + x ^ 2, (1 - 2x) + x ^ 2, 1 - x, (2 - x ^ -1) - x, 1 - x, 0, 0, 1 - x ^ -1, -1 + x, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2x) + x ^ 2, -x + x ^ 2, (1 - 2x) + x ^ 2, 1 - x, 1 - x, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(2 - x ^ -1) - x, 1 - x, 0, 0, 0, 1 - x, -1 + x, -1 + x, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, x - x ^ 2, (-1 + 2x) - x ^ 2, 0, 0, 0, x - x ^ 2, 0, 0, 0, 0, 0, -1 + x, 0, 0, -1 + x, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(2 - x ^ -1) - x, 0, 1 - x, 0, 0, 0, -1 + x, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - x ^ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [((-1 + 3x) - 3 * x ^ 2) + x ^ 3, ((-1 + 3x) - 3 * x ^ 2) + x ^ 3, (x - 2 * x ^ 2) + x ^ 3, (2 - x ^ -1) - x, (-2 + x ^ -1 + 2x) - x ^ 2, (-1 + 2x) - x ^ 2, 0, (1 - 2x) + x ^ 2, -2 + x ^ -1 + x, 0, 0, -1 + x, -1 + x, (2 - x ^ -1) - x, 1 - x, 0, 0, 0, 1 - x, 0, 0, -1 + x, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [(-3 + x ^ -1 + 3x) - x ^ 2, (-1 + 2x) - x ^ 2, 0, (1 - 2x) + x ^ 2, 0, (-1 + 2x) - x ^ 2, (1 - 2x) + x ^ 2, (1 - 2x) + x ^ 2, -1 + x, x - x ^ 2, 0, (-1 + 2x) - x ^ 2, 0, 1 - x, 1 - x, 0, -1 + x, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, -1 + x, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [((1 - 3x) + 3 * x ^ 2) - x ^ 3, ((1 - 2x) + 2 * x ^ 2) - x ^ 3, ((1 - 3x) + 3 * x ^ 2) - x ^ 3, (1 - 2x) + x ^ 2, ((3 - x ^ -1) - 3x) + x ^ 2, (1 - 2x) + x ^ 2, 0, 0, (2 - x ^ -1) - x, 0, (-1 + 2x) - x ^ 2, 0, 1 - x, 0, 0, 0, 0, -1 + x, 0, 0, 1 - x, 0, (1 - 2x) + x ^ 2, 0, 1 - x, 0, -1 + x, x, 0, 0], [0, 0, 0, (2 - x ^ -1) - x, 0, -1 + x, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, -2 + x ^ -1 + x, 0, (2 - x ^ -1) - x, 0, 1 - x, -1 + x ^ -1, 0, -1 + x, 0, 0, 0, 1, 0, 0, 0], [(-1 + 2x) - x ^ 2, 0, 0, 0, 0, 0, (1 - 2x) + x ^ 2, -x + x ^ 2, 0, 0, 0, 0, 0, 0, 0, 1 - x, -1 + x, 0, 0, 0, 0, 0, 1 - x, 0, 0, 1 - x, 0, 0, 0, x], [0, 0, (1 - 2x) + x ^ 2, ((3 - x ^ -1) - 3x) + x ^ 2, 0, 0, 1 - x, (1 - 2x) + x ^ 2, 0, 0, -1 + x, 0, 0, (2 - x ^ -1) - x, 0, 0, 0, 0, 1 - x, 0, 0, 0, (2 - x ^ -1) - x, 0, 1 - x ^ -1, 1 - x, 0, 0, 1, -1 + x]], [[0, -x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, x, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0], [0, x - x ^ 2, 0, 0, 0, 0, 0, 0, 0, (1 - 2x) + x ^ 2, (-1 + 2x) - x ^ 2, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 1 - x, 0, (1 - 2x) + x ^ 2, 0, 0, 0, -1 + x, x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 0, x], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0], [1 - x, 0, 0, (2 - x ^ -1) - x, 0, 0, 0, 0, 0, 0, (-1 + 2x) - x ^ 2, 0, 0, 0, 0, -1 + x, 0, 0, 0, 1 - x, 0, 0, 1 - x, 0, 1, 0, 0, -1 + x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 1 - x ^ -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + x]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, x, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 0, 0, 0, 0, 0, -1 + x, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 1 - x, 0, 0, 0, 1, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, x - x ^ 2, 0, 0, 0, 0, 0, 0, 0, -x + x ^ 2, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -x + x ^ 2, 0, 1 - x, 0, 0, 0, 1, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 1, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 + x]], [[0, 0, -x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, x, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, x, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 0, 0, 0, 0, 0, -1 + x, x, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 0, 0, (2 - x ^ -1) - x, 0, 0, 0, 0, 0, 0, (-1 + 2x) - x ^ 2, 0, 0, 0, 0, -1 + x, 0, 0, 0, 1 - x, 0, 0, 1 - x, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -x + x ^ 2, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0], [0, 0, x - x ^ 2, (-1 + 2x) - x ^ 2, 0, 0, 0, x - x ^ 2, 0, 0, 0, 0, 0, -1 + x, 0, 0, -1 + x, 0, x, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]]]
                        end
                    f2(q) = begin
                            return q ^ 0 * [[[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-(E(3)) + root(-3) * q + E(3, 2) * q ^ 2, E(3, 2) - E(3, 2) * q, 0, 0, -(E(3)) + E(3) * q, 0, (root(-3) - E(3) * q ^ -1) + E(3, 2) * q, -(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [E(3) * q - E(3) * q ^ 2, E(3) * q, 0, 0, 0, 0, E(3) - E(3) * q, E(3) - E(3) * q, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, E(3, 2) - E(3, 2) * q, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, E(3, 2), 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [-(E(3)) + E(3) * q, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, -1 + q, -(E(3)), 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, (-(E(3, 2)) - root(-3) * q) + E(3) * q ^ 2, 0, (-(E(3, 2)) - root(-3) * q) + E(3) * q ^ 2, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, -1 + q, 0, 0], [(E(3) * q - 2 * E(3) * q ^ 2) + E(3) * q ^ 3, E(3) * q - E(3) * q ^ 2, 0, 0, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, E(3, 2) * q], [0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, E(3), -1 + q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q ^ 3 - q ^ 4, 1 - q, -1 + q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, 0, q ^ 3 - q ^ 4, 0, 0, 0, 0, 0, 1 - q, -1 + q, 0, 0, 1, 1 - q, 0, 0, 0, 0, 0, 0], [0, 0, q ^ 3 - q ^ 4, 1 - q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, 0, q ^ 3 - q ^ 4, 0, 0, 0, 0, 0, 1 - q, q, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, -1 + q, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, 0, 0, 0, 1 - q ^ -1, -1 + q, 1, 0, 0], [0, 0, 0, 1 - q, 0, -1 + q, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, q, 0, 0, 0], [(-1 + 2q) - q ^ 2, 0, 0, 0, q - q ^ 2, (1 - 2q) + q ^ 2, 0, 0, 0, 0, 0, -1 + q, 0, 0, -1 + q, 0, 0, 0, -1 + q, q], [(2 - q ^ -1) - q, 0, (-(q ^ 3) + 2 * q ^ 4) - q ^ 5, (-1 + 2q) - q ^ 2, (1 - 2q) + q ^ 2, -2 + q ^ -1 + q, 0, 0, 0, 0, -1 + q, 0, 0, -1 + q, (2 - q ^ -1) - q, -1 + q ^ -1, 0, 0, 1, 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q]], [[0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -(q ^ -2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -(q ^ 3), -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + q, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0, -1 + q, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0], [0, 0, 0, 0, 0, 0, 1 - q, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1 - q, 0, 1 - q, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, q, 1 - q, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[-1 + q, 0, -(q ^ 3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, 0, q ^ 3 - q ^ 4, 0, 0, 0, 0, 0, 1 - q, -1 + q, 0, 0, 1, 1 - q, 0, 0, 0, 0, 0, 0], [-(q ^ -2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [q - q ^ 2, q, 0, 0, 0, 0, 1 - q, 1 - q, 0, 0, 1 - q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, q ^ 3 - q ^ 4, 1 - q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]]]
                        end
                    f3(q, j) = begin
                            return q ^ 0 * [[[-1, 0, 0, 0, 0], [0, -1, 0, 0, 0], [1, 0, 0, 0, -1], [1 + j * q, 0, 1 + j * q, -1, -1 - j * q], [-q, 0, -q, 0, -1 + q]], [[-1, 0, 0, 0, 0], [0, -1, 0, 0, 0], [1, 0, 0, -1, 0], [-q, 0, -q, -1 + q, 0], [-j, 0, -j, j, -1]], [[-1, 0, -1, 0, 0], [0, -1, 1, 0, 0], [0, 0, q, 0, 0], [0, 0, 0, -1, 0], [0, 0, 0, 0, -1]], [[q, 0, 0, 0, 0], [-1, -1, 0, 0, 0], [-q, 0, -1, 0, 0], [1, 0, 0, -1, 0], [1, 0, 0, 0, -1]], [[0, 1, 0, 0, 0], [q, -1 + q, 0, 0, 0], [0, 0, -1, 0, 0], [1, 1, 0, -1, 0], [1, 1, 0, 0, -1]]]
                        end
                    f4(q, j) = begin
                            return q ^ 0 * [[[-1 + q, 0, 0, q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, (j ^ 2 * q + (-(j ^ 2) + j) * q ^ 2) - j * q ^ 3, 0, j - j * q, 0, (-j + j * q) - j * q ^ 2, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, j ^ 2 * q - j ^ 2 * q ^ 2, 0, -(j ^ 2), 0, j ^ 2 - j ^ 2 * q, 0, 0, 0], [-q + q ^ 2, 0, (j ^ 2 * q - 2 * j ^ 2 * q ^ 2) + j ^ 2 * q ^ 3, -q + q ^ 2, -(j ^ 2) + j ^ 2 * q, 0, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, -1, 0, 0], [0, q - q ^ 2, (j ^ 2 * q - 2 * j ^ 2 * q ^ 2) + j ^ 2 * q ^ 3, 0, -(j ^ 2) + j ^ 2 * q, 0, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 0, -1, q - q ^ 2], [0, q, 0, 0, 0, 0, 0, 0, 0, -1 + q]], [[0, 0, j ^ 2 * q - j ^ 2 * q ^ 2, j ^ 2 * q, 0, 0, 0, 0, 0, 0], [0, -1 + q, -q + q ^ 2, 0, 0, 0, 0, 0, 0, j], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [j, 0, q - q ^ 2, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1 + q, 0, -q, 0, 0, 0], [q - q ^ 2, -(j ^ 2) * q + j ^ 2 * q ^ 2, ((j ^ 2 * q - j ^ 2 * q ^ 2) - j ^ 2 * q ^ 3) + j ^ 2 * q ^ 4, j ^ 2 * q ^ 2 - j ^ 2 * q ^ 3, 0, -1, 0, 0, 0, -1 + q], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, j ^ 2 * q, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, 0, 0, 0, 0, 0, 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, q, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, -1 + q, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, -q, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, -1 + q, 0, q, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, -q, -1 + q, 0], [q, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, -q + q ^ 2, -(q ^ 2) + q ^ 3, 0, -1 + q, 0, 0, 0, 0, j * q], [0, 0, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, 0, j ^ 2 - j ^ 2 * q, -(j ^ 2), 0, 0], [0, -q, 0, 0, 0, 0, -1 + q, 0, 0, 0], [0, -q + q ^ 2, 0, q ^ 2 - q ^ 3, 0, -j * q, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, j ^ 2, 0, -(j ^ 2) + j ^ 2 * q, 0, 0, 0]]]
                        end
                    f11(q, j) = begin
                            return q ^ 0 * [[[q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-(j ^ 2) + j ^ 2 * q, j ^ 2 - j ^ 2 * q, j ^ 2, 0, 0, 0, 0, 0, 0, 0], [-j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, (j - j * q) + j * q ^ 2, j - j * q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, q, -1 + q, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, -1 + q, 0, j, 0, 0, 0, 0], [1 - q, 0, 0, 0, -1 + q, 0, j, 0, 0, 0], [-(j ^ 2) * q + j ^ 2 * q ^ 2, 0, 0, j ^ 2 * q, 0, 0, 0, 0, 0, 0], [-(j ^ 2) * q + j ^ 2 * q ^ 2, 0, 0, 0, j ^ 2 * q, 0, 0, 0, 0, 0], [(-(j ^ 2) + 2 * j ^ 2 * q) - j ^ 2 * q ^ 2, j ^ 2 - j ^ 2 * q, j ^ 2 - j ^ 2 * q, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, -1 + q, 0, -1, 0, 0], [0, 0, 0, j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) * q + j ^ 2 * q ^ 2, 1 - q, -1 + q, 0, -1, 0], [(-(j ^ 2) + 2 * j ^ 2 * q) - j ^ 2 * q ^ 2, j ^ 2 - j ^ 2 * q, j ^ 2 - j ^ 2 * q, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, -1 + q, 0, 0, -1]], [[0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [-q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, q, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, q, 0, 0, -1 + q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [-1 + q, 0, 0, 1 - q, 0, -j, 0, 0, 0, 0], [0, -q, 0, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) * q + j ^ 2 * q ^ 2, -(j ^ 2) * q, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1 + q, -q], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, -q, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, -q, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1], [0, 0, 0, 0, 0, 0, 0, 0, q, 0], [0, 0, 0, 0, 0, 0, 0, -q, 0, -1 + q]]]
                        end
                    f12(q, j) = begin
                            return q ^ 0 * [[[0, 0, 0, 0, -(j ^ 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -(j ^ 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -(j ^ 2), 0, 0, 0, 0, 0, 0, 0, 0], [-j * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -j, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -j * q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [-(j ^ 2) + j ^ 2 * q, 0, 0, 0, j - j * q, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, j ^ 2 - j ^ 2 * q, 0, 0, -j + j * q, 0, -1, 0, 0, 0, 0, 0, 0], [0, (2 * j ^ 2 - j ^ 2 * q ^ -1) - j ^ 2 * q, 0, 0, 0, (j - 2 * j * q) + j * q ^ 2, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, (-1 + 2q) - q ^ 2, -q, 0, 0, (1 - 2q) + q ^ 2, 0, 0, 0, q, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [j - j * q, 0, 0, -q + q ^ 2, (j ^ 2 + 2j + q ^ -1) - j * q, 0, 1 - q, -q, -q, 0, 0, 0, q, 0, 0], [0, -2 * j ^ 2 + j ^ 2 * q ^ -1 + j ^ 2 * q, 0, 0, 0, (-j + 2 * j * q) - j * q ^ 2, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, (1 - 2q) + q ^ 2, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, 0, q, 0, q, q]], [[-1 + q, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + q, 0, -1, 0, -1 + q ^ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1 + q, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [-q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -q + q ^ 2, 0, 0, 1 - q, 0, -1, 0, 0, 0, 0, 0, 0], [0, -q + q ^ 2, 0, 0, 0, q - q ^ 2, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, (-1 + 2q) - q ^ 2, -q, 1 - q, 1 - q, -q + q ^ 2, 1 - q, 0, 0, q, q, 0, 0, 0, 0], [0, 0, 0, -1 + q, 0, 0, -1 + q ^ -1, 0, 0, 0, 0, -1, 0, 0, 0], [-(j ^ 2) + j ^ 2 * q, 1 - q ^ -1, 0, (1 - 2q) + q ^ 2, -(j ^ 2) + j ^ 2 * q, -1 + q, 1 - q, -q, -q, 0, 0, 0, q, 0, 0], [0, q - q ^ 2, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, -1, 0], [1 - q ^ -1, (((2 + q ^ -2) - 2 * q ^ -1) - 2q) + q ^ 2, 0, (2 - q ^ -1) - q, 1 - q ^ -1, (-2 + q ^ -1 + 2q) - q ^ 2, 1 - q ^ -1, 0, 0, 0, 0, q, 0, q, q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q, -1, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, (1 - 2q) + q ^ 2, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, -q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, -j + q ^ -2 + (j ^ 2 + 2j) * q ^ -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [-j + (j ^ 2 + 2j) * q + q ^ 2, 1 - q, 0, 0, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, -1 + q, q, q, 0, 1 - q, 0, -1, 0, 0], [0, (((3 - q ^ -1) - 4q) + 3 * q ^ 2) - q ^ 3, 0, 0, 0, ((-2 + 3q) - 3 * q ^ 2) + q ^ 3, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, ((-4 + q ^ -1 + 6q) - 4 * q ^ 2) + q ^ 3, 0, 0, 0, ((2 - 5q) + 4 * q ^ 2) - q ^ 3, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1 + q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-j + j * q, 0, 0, 0, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, q ^ -1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [-j + j * q, 0, 1, 0, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, 0, 0, 1 - q, 0, 0, -1 + q, -1, 0, 0, 0, 0, 0], [0, 0, 0, -(q ^ 2) + q ^ 3, 0, 0, q - q ^ 2, 0, -q, 0, 0, 0, 0, 0, 0], [-j + j * q, 0, 0, q - q ^ 2, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, -1 + q, q, q, 0, 0, 0, -1, 0, 0], [0, (2 - q ^ -1) - q, 0, 0, 0, -2 + q ^ -1 + q, 0, 0, 0, 0, 0, -1 + q, 0, -1, 0], [0, (1 - 2q) + q ^ 2, q, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, -q, -q, 0, -1 + q, 0, 0], [0, 0, 0, (-q + 2 * q ^ 2) - q ^ 3, 0, 0, (1 - 2q) + q ^ 2, 0, 0, 0, 0, -q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, (-1 + 2q) - q ^ 2, 0, 0, 0, (1 - 2q) + q ^ 2, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, ((2 * j ^ 2 + j) - j ^ 2 * q ^ -1) + q, 0, 0, 0, -(j ^ 2) + j ^ 2 * q, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q ^ -1, 0, 0, 0, 0, 1, 0, 0, 0], [(1 - 2q) + q ^ 2, 0, 1, 0, (2 - q ^ -1) - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, q - q ^ 2, 0, 0, -1 + q, 0, q, 0, 0, -1 + q, 0, 0, 0], [0, (1 - 2q) + q ^ 2, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, 0, q, 0, q, q], [(j ^ 2 + (-2 * j ^ 2 - j) * q) - q ^ 2, 0, 0, 0, j ^ 2 - j ^ 2 * q, 0, 0, q, 0, 0, 0, 0, 0, -1 + q, 0], [j - j * q, 0, 0, -q + q ^ 2, (j ^ 2 + 2j + q ^ -1) - j * q, 0, 1 - q, -q, -q, 0, 0, 0, 1, 0, -1 + q]]]
                        end
                    f13(q, j) = begin
                            return q ^ 0 * [[[j - j * q, (-1 + q) - q ^ 2, (-(j ^ 2) + 2 * j ^ 2 * q) - j ^ 2 * q ^ 2, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, -(j ^ 2) + (j ^ 2 - j) * q + j * q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, j ^ 2 - j ^ 2 * q, -1 + q, 1 - q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, 0, -q + q ^ 2, 0, -q + q ^ 2, 1 - q, 0, j ^ 2 - j ^ 2 * q, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 1 - q, -1 + q, 0, 0, -1 + q, 0, j ^ 2 - j ^ 2 * q, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, (-(j ^ 2) + 2 * j ^ 2 * q) - j ^ 2 * q ^ 2, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 0, 0, (j ^ 2 + (-(j ^ 2) + j) * q) - j * q ^ 2, 0, (1 - q) + q ^ 2, 0, j - j * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j ^ 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1 - q, -1 + q, 1 - q, 0, 0, 0, 0, 0, j ^ 2 - j ^ 2 * q, 0, j ^ 2 * q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, (2 - q ^ -1) - q, ((2 * j ^ 2 + j) - j ^ 2 * q ^ -1) + q, ((-2 * j ^ 2 - j) - q ^ -1) + j ^ 2 * q, 0, 0, 0, 0, 0, -j + j * q ^ -1 + j * q, 0, j - j * q, 0, 0, 0, 0, 0], [0, (-1 + q) - q ^ 2, -q + q ^ 2, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 1 - q, ((j ^ 2 - j) + j * q ^ -1) - j ^ 2 * q, 0, -1 + q ^ -1 + q, 0, 0, 0, 0, 0, 0, 0, j - j * q, 0, 0, 0, 0], [0, 0, 0, 0, ((j ^ 2 - j) - j ^ 2 * q ^ -1) + j * q, 0, -2j + j * q ^ -1 + j * q, 0, 0, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, 0, 0, 0, 0, 0, 0, j - j * q, (1 - q ^ -1) - q, 0, 0], [0, 0, 0, 0, -1 + q, 0, -1 + q, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, -q, j ^ 2 - j ^ 2 * q, 0, 0], [0, (j * q + (j ^ 2 - j) * q ^ 2) - j ^ 2 * q ^ 3, (-q + 2 * q ^ 2) - q ^ 3, 0, (1 - 2q) + q ^ 2, 0, ((1 - q) - q ^ 2) + q ^ 3, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, (j * q + (j ^ 2 - j) * q ^ 2) - j ^ 2 * q ^ 3, 0, -q + q ^ 2, 0, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 0, j ^ 2 * q - j ^ 2 * q ^ 2, q - q ^ 2, 0, 0, -1, 0], [1 - q, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, 0, ((-1 + q ^ -1) - q) + q ^ 2, 0, ((-1 + q ^ -1) - q) + q ^ 2, 0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, (2 - q ^ -1) - q, -1 + q, -j + j * q ^ -1, 0, -1 + q ^ -1, 0, 0, 1 - q, ((j ^ 2 - j) + j * q ^ -1) - j ^ 2 * q, 0, -1]], [[-1 + q, -(j ^ 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-j * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, j ^ 2 * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, 1 - q, 0, 0, -1 + q, 0, -1 + q, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, 0, 0, 0, q - q ^ 2, 0, q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0], [0, -(j ^ 2) * q, 0, 0, 0, 0, 0, j ^ 2, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, -(j ^ 2), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -j * q, 0, 0, 0], [0, j - j * q, 0, 0, 0, 0, 0, j - j * q ^ -1, j - j * q, 0, j ^ 2 - j ^ 2 * q, 0, j - j * q ^ -1, 0, -j + j * q, -(j ^ 2) + j ^ 2 * q, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, 0, j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, j ^ 2 * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, j ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, j ^ 2 * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, j * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, j, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q + q ^ 2, q - q ^ 2, 0, 0, q - q ^ 2, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, -q, 0, 0, q, 0, 0, 0, 0, 0, 0], [0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, (-1 + 2q) - q ^ 2, 0, -2 + q ^ -1 + q, 0, ((-1 + q ^ -1) - q) + q ^ 2, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, -1 + q, 0, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, 0, j ^ 2 - j ^ 2 * q, 1 - q, 0, 0, -1, 0], [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [-1 + q, ((-j + 3 * j * q) - 3 * j * q ^ 2) + j * q ^ 3, ((j - 3 * j * q) + 3 * j * q ^ 2) - j * q ^ 3, (-q + 2 * q ^ 2) - q ^ 3, ((j ^ 2 + 3j) - j * q ^ -1) + (-2 * j ^ 2 - 4j) * q + (j ^ 2 + 2j) * q ^ 2, -(j ^ 2) + (2 * j ^ 2 + j) * q + q ^ 2, ((j ^ 2 + 2j) - j * q ^ -1) + (-2 * j ^ 2 - j) * q + (j ^ 2 - j) * q ^ 2 + j * q ^ 3, (((-2 * j ^ 2 - 3j) - q ^ -1) + (j ^ 2 + 3j) * q) - j * q ^ 2, (-j + (j ^ 2 + 3j) * q + (-2 * j ^ 2 - 3j) * q ^ 2) - q ^ 3, 0, -(j ^ 2) + (j ^ 2 - j) * q + j * q ^ 2, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, (-1 + q) - q ^ 2, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 0, 0, j - j * q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, 0, -1 + q]], [[q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j ^ 2 * q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(j ^ 2), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1 - q, 0, 1 - q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, q, -(j ^ 2) + j ^ 2 * q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, 1 - q, 0, 0, -1 + q, 0, -1 + q, 0, 1, 0, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -j * q, 0, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, -j * q + (-(j ^ 2) + j) * q ^ 2 + j ^ 2 * q ^ 3, (q - 2 * q ^ 2) + q ^ 3, 0, (-1 + 2q) - q ^ 2, 0, (-1 + q + q ^ 2) - q ^ 3, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, -j * q + (-(j ^ 2) + j) * q ^ 2 + j ^ 2 * q ^ 3, 0, q - q ^ 2, 0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, -q + q ^ 2, 0, 0, q, q]], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, 0, 0, 0, 0, 0, 0, 0], [0, q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, j ^ 2, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0], [-q, 0, -q + q ^ 2, 0, -q + q ^ 2, 1 - q, 0, j ^ 2 - j ^ 2 * q, 0, 0, 0, 0, 0, -1 + q, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, 1 - q, 0, 0, -1 + q, -1, -1 + q, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j ^ 2, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1 + q, 0, 0], [q - q ^ 2, (j * q + (j ^ 2 - j) * q ^ 2) - j ^ 2 * q ^ 3, 0, 0, ((1 - q) - q ^ 2) + q ^ 3, 0, ((1 - q) - q ^ 2) + q ^ 3, 0, (j * q + (j ^ 2 - j) * q ^ 2) - j ^ 2 * q ^ 3, (-1 + 2q) - q ^ 2, -q + q ^ 2, j - j * q, 0, 1 - q, 0, 0, q - q ^ 2, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, -q], [0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, (-1 + 2q) - q ^ 2, 0, -2 + q ^ -1 + q, 0, ((-1 + q ^ -1) - q) + q ^ 2, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, -1 + q, 0, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, 0, j ^ 2 - j ^ 2 * q, 1 - q, 0, 0, -1, -1 + q]]]
                        end
                    f17(q) = begin
                            return q ^ 0 * [[[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [-(q ^ 4), q ^ 3, -(q ^ 2), q]], [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [-(q ^ 4), q ^ 3, -(q ^ 2), q]], [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 0, 1], [0, 0, q, -1 + q]], [[-1, 0, 0, 0], [0, 0, 1, 0], [0, q, -1 + q, 0], [0, 0, 0, -1]], [[0, 1, 0, 0], [q, -1 + q, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]]
                        end
                    f20(q, j) = begin
                            return q ^ 0 * [[[q, 0, -(q ^ 2), 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -(q ^ 2), q, 0, 0, 0, 0, 0], [0, -(q ^ 2), 0, 0, 0, q, 0, 0, 0, 0], [-j * q + (-(j ^ 2) + j) * q ^ 2 + j ^ 2 * q ^ 3, (j ^ 2 * q + (-(j ^ 2) + j) * q ^ 2) - j * q ^ 3, (-(j ^ 2) * q ^ 2 + 2 * j ^ 2 * q ^ 3) - j ^ 2 * q ^ 4, (j ^ 2 * q - 2 * j ^ 2 * q ^ 2) + j ^ 2 * q ^ 3, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, j - j * q, (-j + j * q) - j * q ^ 2, 0, 0], [j ^ 2 * q - j ^ 2 * q ^ 2, j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) * q ^ 2 + j ^ 2 * q ^ 3, j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) + j ^ 2 * q, 0, -(j ^ 2), j ^ 2 - j ^ 2 * q, 0, 0], [-j * q ^ 2 + (-(j ^ 2) + j) * q ^ 3 + j ^ 2 * q ^ 4, (-j * q ^ 2 + 2 * j * q ^ 3) - j * q ^ 4, (j ^ 2 * q ^ 2 + (-2 * j ^ 2 + j) * q ^ 3 + (2 * j ^ 2 - j) * q ^ 4) - j ^ 2 * q ^ 5, (j ^ 2 * q ^ 2 + (-(j ^ 2) + j) * q ^ 3) - j * q ^ 4, 0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, (-j * q + j * q ^ 2) - j * q ^ 3, j - j * q, ((j - j * q) + j * q ^ 2) // q], [0, 0, -(j ^ 2) * q ^ 3 + j ^ 2 * q ^ 4, 0, -(j ^ 2) * q ^ 2 + j ^ 2 * q ^ 3, j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) * q ^ 2, 0, j ^ 2 * q, j ^ 2 - j ^ 2 * q]], [[q, 0, -(q ^ 2), 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -(q ^ 2), q, 0, 0, 0, 0, 0], [0, -(q ^ 2), 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1 + q, -q, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -(q ^ 2), -1 + q, 1], [0, 0, 0, 0, 0, 0, -(q ^ 2), 0, q, 0]], [[-1 + q, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, q, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -q, 0, 1 // q], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, -q, 0, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, -(q ^ 2), 0, q, 0], [0, -(q ^ 3), 0, 0, 0, q ^ 2, 0, 0, 0, -1 + q]], [[q, 0, -(q ^ 2), 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, q, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, q, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, q, 0, -1 // q], [0, 0, 0, 0, 0, 0, q, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, q, -1 + q, 0, 0, 0, 0, 0, 0], [q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, -(q ^ 2), 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, -(q ^ 2), 0, q, 0], [0, 0, 0, 0, 0, 0, 0, -(q ^ 2), 0, q]]]
                        end
                    f23(q) = begin
                            return q ^ 0 * [[[-1, 0, 0, 0, 0], [0, -1, 0, 0, 0], [q ^ 3, -(q ^ 2), q, 0, 0], [0, 0, 0, -1, 0], [-(q ^ 3), 0, 0, -(q ^ 2), q]], [[-1, 0, 0, 0, 0], [0, -1, 0, 0, 0], [q ^ 3, -(q ^ 2), q, 0, 0], [0, 0, 0, -1, 0], [-(q ^ 3), 0, 0, -(q ^ 2), q]], [[-1, 0, 0, 0, 0], [0, 0, 1, 0, 0], [0, q, -1 + q, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, q, -1 + q]], [[0, 1, 0, 0, 0], [q, -1 + q, 0, 0, 0], [0, 0, -1, 0, 0], [0, 0, 0, -1, 0], [-(q ^ 4), q ^ 3, -(q ^ 2), -(q ^ 2), q]], [[-1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [0, q, 0, -1 + q, 0], [0, 0, q, 0, -1 + q]]]
                        end
                    f29(q) = begin
                            return q ^ 0 * [[[q, q ^ 4, 0, 0, -(q ^ 2), 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, q ^ 3, -(q ^ 2), q, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, q ^ 4, 0, -(q ^ 3), q]], [[q, q ^ 4, 0, 0, -(q ^ 2), 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, q ^ 3, -(q ^ 2), q, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, q ^ 4, 0, -(q ^ 3), q]], [[-1 + q, 0, 0, 0, q, 0], [0, -1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, q, -1 + q, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q]], [[0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0], [0, q, -1 + q, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, q, 0], [q, 0, 0, 0, 0, -1 + q]], [[-1 + q, 0, 0, q, 0, 0], [0, q, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0], [0, 0, q, 0, -1 + q, 0], [0, 0, 0, 0, 0, -1]]]
                        end
                    f9(q) = begin
                            return q ^ 0 * [[[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (E(3) + (E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, (E(3, 2) + (-(E(3, 2)) + E(3)) * q) - E(3) * q ^ 2, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0], [(-1 + 2q) - q ^ 2, (1 - 2q) + q ^ 2, (-1 + 2q) - q ^ 2, -1 + q, 0, -1 + q, 0, 1 - q, 0, E(3) - E(3) * q, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, (1 - 2q) + q ^ 2, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(-3 + q ^ -1 + 3q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, (-3 + q ^ -1 + 3q) - q ^ 2, -2 + q ^ -1 + q, 0, (E(3, 2) + 2 * E(3) + q ^ -1) - E(3) * q, 0, (2 - q ^ -1) - q, 0, -(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ((-(E(3, 2)) - 2 * E(3)) + E(3) * q ^ -1) - q, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, 0], [-(E(3)) + E(3) * q, E(3) - E(3) * q, 0, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, 0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, 0, E(3) - E(3) * q, 0, 0, 0, E(3, 2), 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, E(3, 2) - E(3, 2) * q, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q ^ -1, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, (E(3) + (-(E(3, 2)) - 2 * E(3)) * q) - q ^ 2, 0, ((3 - q ^ -1) - 3q) + q ^ 2, 0, 0, 0, 0, 0, 0, -1 + q, (1 + (E(3, 2) + 2 * E(3)) * q) - E(3) * q ^ 2, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (((E(3, 2) + 3 * E(3)) - E(3) * q ^ -1) + (-2 * E(3, 2) - 3 * E(3)) * q) - q ^ 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, -2 * E(3, 2) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, 0, 0, 0, 0, (E(3, 2) + (-(E(3, 2)) + E(3)) * q) - E(3) * q ^ 2, 0, 0, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (E(3) + (E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, 0, 0, 0], [0, 0, -(E(3)) + E(3) * q, 0, E(3) - E(3) * q, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, 0, E(3) - E(3) * q, E(3) - E(3) * q, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + q, 0, 0, -2 + q ^ -1 + q, 0, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, -1 + q, -2 + q ^ -1 + q, 0, 0, -2 + q ^ -1 + q, 0, 1 - q, 0, 0, 0, E(3) - E(3) * q, 1 - q, 0, 0, E(3), 0, 0, -1 + q ^ -1, 0, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1 - q ^ -1, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 1 - q, 0, 1, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -2 * E(3, 2) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, (2 * E(3, 2) + E(3) + q ^ -1) - E(3, 2) * q, 0, ((-2 * E(3, 2) - E(3)) + E(3, 2) * q ^ -1) - q, 0, 0, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, 0, 0, 0, ((2 * E(3, 2) - E(3, 2) * q ^ -1) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, 0, 0, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, -1 + q ^ -1 + q, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, 0, 0, 0], [(1 + (E(3, 2) + 2 * E(3)) * q) - E(3) * q ^ 2, 0, 0, (((-2 * E(3, 2) - 3 * E(3)) - q ^ -1) + (E(3, 2) + 3 * E(3)) * q) - E(3) * q ^ 2, 0, 0, 0, 0, 0, 0, (-3 + q ^ -1 + 3q) - q ^ 2, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, 0, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, (E(3) + (-(E(3, 2)) - 2 * E(3)) * q) - q ^ 2, 0, 0, 0, (E(3, 2) - E(3, 2) * q) + E(3, 2) * q ^ 2, (-1 + 2q) - q ^ 2, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, (2 - q ^ -1) - q, 0, 0, 0, 0, 0, 0, (-3 + q ^ -1 + 3q) - q ^ 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)), 0, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, 0, 0, 0, E(3) - E(3) * q, 0, q, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -q + q ^ 2, 0, 0, 1 - q, 0, 0, 0, (-(E(3)) + E(3) * q) - E(3) * q ^ 2, 0, (E(3) + (E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, -(E(3)) + (-(E(3, 2)) + E(3)) * q + E(3, 2) * q ^ 2, 0, 0, 0, 0], [0, 0, ((E(3, 2) - E(3)) + E(3) * q ^ -1) - E(3, 2) * q, 0, 0, -2 + q ^ -1 + q, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, 0, 0, 0, 0, ((3 - q ^ -1) - 3q) + q ^ 2, (-(E(3, 2)) + E(3) + E(3, 2) * q ^ -1) - E(3) * q, 0, 0, 0, 0, -2 + q ^ -1 + q, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, 0, 0, -1 + q ^ -1 + q, 0, E(3, 2) - E(3, 2) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, (-(E(3)) + E(3) * q) - E(3) * q ^ 2, (E(3) + (E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, 0, 0, 0, 0, 0, (1 - q) + q ^ 2, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, -(E(3)) + (-(E(3, 2)) + E(3)) * q + E(3, 2) * q ^ 2, 0, 0, 0], [(-1 + 2q) - q ^ 2, (((2 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-(E(3, 2)) + 2 * E(3)) * q) - E(3) * q ^ 2, E(3) + (E(3, 2) - 2 * E(3)) * q + (-2 * E(3, 2) + E(3)) * q ^ 2 + E(3, 2) * q ^ 3, ((((7 * E(3, 2) + E(3, 2) * q ^ -2) - 4 * E(3, 2) * q ^ -1) - 7 * E(3, 2) * q) + 4 * E(3, 2) * q ^ 2) - E(3, 2) * q ^ 3, 0, (((-(E(3, 2)) + 2 * E(3)) - E(3) * q ^ -1) + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, ((3 * E(3, 2) - E(3, 2) * q ^ -1) - 3 * E(3, 2) * q) + E(3, 2) * q ^ 2, ((2 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-2 * E(3, 2) + E(3)) * q + E(3, 2) * q ^ 2, 0, ((2 - q ^ -1) - 2q) + q ^ 2, (((-5 * E(3, 2) + 2 * E(3)) - E(3, 2) * q ^ -2) + (3 * E(3, 2) - E(3)) * q ^ -1 + (5 * E(3, 2) - E(3)) * q) - 2 * E(3, 2) * q ^ 2, (((3 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-4 * E(3, 2) + 2 * E(3)) * q + (3 * E(3, 2) - E(3)) * q ^ 2) - E(3, 2) * q ^ 3, (E(3, 2) - 2 * E(3)) + E(3, 2) * q ^ -2 + (-2 * E(3, 2) + E(3)) * q ^ -1 + E(3) * q, (-1 + 2q) - q ^ 2, ((-(E(3, 2)) + 2 * E(3, 2) * q) - 2 * E(3, 2) * q ^ 2) + E(3, 2) * q ^ 3, (2 * E(3, 2) - 2 * E(3)) + E(3, 2) * q ^ -2 + (-2 * E(3, 2) + E(3)) * q ^ -1 + (-2 * E(3, 2) + E(3)) * q + E(3, 2) * q ^ 2, 0, -2 * E(3, 2) + E(3) + E(3, 2) * q ^ -1 + (E(3, 2) - 2 * E(3)) * q + E(3) * q ^ 2, (-1 + 2q) - q ^ 2, (-2 * E(3, 2) + E(3) + E(3, 2) * q ^ -1 + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, 0, (-2 + q ^ -1 + 2q) - q ^ 2, -3 * E(3, 2) + E(3) + E(3, 2) * q ^ -1 + (4 * E(3, 2) - 2 * E(3)) * q + (-3 * E(3, 2) + E(3)) * q ^ 2 + E(3, 2) * q ^ 3, (-(E(3, 2)) - 3 * E(3)) + E(3) * q ^ -1 + (2 * E(3, 2) + 4 * E(3)) * q + (-(E(3, 2)) - 3 * E(3)) * q ^ 2 + E(3) * q ^ 3, 0, -2 + q ^ -1 + q, (E(3) - E(3) * q) + E(3) * q ^ 2, (-2 + q ^ -1 + 2q) - q ^ 2, ((2 * E(3, 2) - E(3)) + E(3, 2) * q ^ -2 + (-2 * E(3, 2) + E(3)) * q ^ -1) - E(3, 2) * q, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (-2 * E(3, 2) + E(3) + E(3, 2) * q ^ -1 + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, E(3) - E(3) * q, (1 - q) + q ^ 2, 0, 0, ((((-5 * E(3, 2) + 2 * E(3)) - E(3, 2) * q ^ -2) + (3 * E(3, 2) - E(3)) * q ^ -1 + (6 * E(3, 2) - E(3)) * q) - 4 * E(3, 2) * q ^ 2) + E(3, 2) * q ^ 3, 0, (((2 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-2 * E(3, 2) + 2 * E(3)) * q + (2 * E(3, 2) - E(3)) * q ^ 2) - E(3, 2) * q ^ 3, 0, 0], [0, -2 + q ^ -1 + q, (1 - 2q) + q ^ 2, ((-4 - q ^ -2) + 3 * q ^ -1 + 3q) - q ^ 2, 0, (2 - q ^ -1) - q, -2 + q ^ -1 + q, -2 + q ^ -1 + q, 0, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, (((-3 * E(3, 2) - 4 * E(3)) + q ^ -2) - 3 * q ^ -1) + (E(3, 2) + 2 * E(3)) * q, (-3 + q ^ -1 + 3q) - q ^ 2, (-1 - q ^ -2) + 2 * q ^ -1, E(3) - E(3) * q, (1 - 2q) + q ^ 2, ((3 * E(3, 2) + 2 * E(3)) - q ^ -2) + (-3 * E(3, 2) - 2 * E(3)) * q ^ -1 + q, 0, (2 - q ^ -1) - q, E(3) - E(3) * q, (2 - q ^ -1) - q, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, ((3 - q ^ -1) - 3q) + q ^ 2, ((2 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-2 * E(3, 2) + E(3)) * q + E(3, 2) * q ^ 2, 0, E(3) - E(3) * q ^ -1, 1 - q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, (-1 - q ^ -2) + 2 * q ^ -1, -1 + q, (2 - q ^ -1) - q, 1, E(3, 2) - E(3, 2) * q, 1, 0, (((4 + q ^ -2) - 3 * q ^ -1) - 3q) + q ^ 2, 0, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)) * q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [-(E(3, 2)) + E(3, 2) * q, E(3, 2) - E(3, 2) * q, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0], [0, 0, 0, -1 + q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 1 - q, 0, E(3) - E(3) * q, E(3) * q], [((3 - q ^ -1) - 3q) + q ^ 2, (-3 + q ^ -1 + 3q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, -(E(3, 2)) + E(3, 2) * q, 0, -(E(3, 2)) + E(3, 2) * q ^ -1, 0, -2 + q ^ -1 + q, -1 + q, (E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, ((E(3, 2) + 2 * E(3)) - E(3) * q ^ -1) + q, ((-(E(3, 2)) - 2 * E(3)) + E(3) * q ^ -1) - q, (-2 + q ^ -1 + 2q) - q ^ 2, -(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, E(3, 2) - E(3, 2) * q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0], [1 - q, 0, 1 - q, 0, 1 - q ^ -1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -2 + q ^ -1 + q, 0, 0, -1 + q ^ -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, (2 - q ^ -1) - q, (-1 + 2q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, -2 + q ^ -1 + q, 0, (2 - q ^ -1) - q, (2 - q ^ -1) - q, -(E(3)) + E(3) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, E(3) - E(3) * q, (1 - 2q) + q ^ 2, 0, 0, (-1 + 2q) - q ^ 2, 0, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, -(E(3)) + E(3) * q, -2 + q ^ -1 + q, -2 + q ^ -1 + q, 0, (-1 + 2q) - q ^ 2, 0, E(3, 2) - E(3, 2) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, -1, -2 + q ^ -1 + q, (-1 + 2q) - q ^ 2, -2 + q ^ -1 + q, ((-3 * E(3, 2) - 2 * E(3)) - q ^ -1) + (3 * E(3, 2) + 2 * E(3)) * q + q ^ 2, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, E(3) - E(3) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, -1 + q, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0], [(2 * E(3) - E(3) * q ^ -1) - E(3) * q, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, (-3 * E(3) + E(3) * q ^ -1 + 3 * E(3) * q) - E(3) * q ^ 2, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, E(3, 2) - E(3, 2) * q, -2 * E(3, 2) + E(3, 2) * q ^ -1 + E(3, 2) * q, -(E(3, 2)) + E(3, 2) * q, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, ((E(3, 2) + 3 * E(3)) - E(3) * q ^ -1) + (-(E(3, 2)) - 3 * E(3)) * q + E(3) * q ^ 2, 0, E(3, 2) - E(3, 2) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, -1 + q, 0, E(3) - E(3) * q, 0, 0, 0, 0, E(3), 0, E(3), (2 * E(3) - E(3) * q ^ -1) - E(3) * q, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, ((-(E(3, 2)) - 3 * E(3)) + E(3) * q ^ -1 + (E(3, 2) + 3 * E(3)) * q) - E(3) * q ^ 2, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, -(E(3, 2)) + E(3, 2) * q], [(-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, (-1 + 2q) - q ^ 2, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, (2 - q ^ -1) - q, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, E(3) - E(3) * q, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, -(E(3)) + E(3) * q, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (-1 + 2q) - q ^ 2, 0, 0, -q + q ^ 2, 0, (1 - 2q) + q ^ 2, 0, 0, -(E(3)) + E(3) * q, 1 - q, 0, E(3) * q - E(3) * q ^ 2, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (1 - 2q) + q ^ 2, (-1 + 2q) - q ^ 2, E(3) - E(3) * q, 0], [1 - q, 0, 0, 0, 1 - q, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1 + q]], [[0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, (2 - q ^ -1) - q, (-1 + 2q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, -2 + q ^ -1 + q, 0, (2 - q ^ -1) - q, (2 - q ^ -1) - q, -(E(3)) + E(3) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, E(3) - E(3) * q, (1 - 2q) + q ^ 2, 0, 0, (-1 + 2q) - q ^ 2, 0, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, -(E(3)) + E(3) * q, -2 + q ^ -1 + q, -2 + q ^ -1 + q, 0, (-1 + 2q) - q ^ 2, 0, E(3, 2) - E(3, 2) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, -1, -2 + q ^ -1 + q, (-1 + 2q) - q ^ 2, -2 + q ^ -1 + q, ((-3 * E(3, 2) - 2 * E(3)) - q ^ -1) + (3 * E(3, 2) + 2 * E(3)) * q + q ^ 2, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, E(3) - E(3) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [(-3 * E(3, 2) + E(3, 2) * q ^ -1 + (4 * E(3, 2) + E(3)) * q + (-3 * E(3, 2) - 2 * E(3)) * q ^ 2) - q ^ 3, (3 * E(3, 2) - E(3, 2) * q ^ -1) + (-4 * E(3, 2) - E(3)) * q + (3 * E(3, 2) + 2 * E(3)) * q ^ 2 + q ^ 3, ((-4 * E(3, 2) - 2 * E(3)) - q ^ -1) + (6 * E(3, 2) + 2 * E(3)) * q + (-4 * E(3, 2) - E(3)) * q ^ 2 + E(3, 2) * q ^ 3, (E(3, 2) - 3 * E(3)) + E(3) * q ^ -1 + (-2 * E(3, 2) + 4 * E(3)) * q + (E(3, 2) - 3 * E(3)) * q ^ 2 + E(3) * q ^ 3, ((2 * E(3) - E(3) * q ^ -1) - 2 * E(3) * q) + E(3) * q ^ 2, (2 - q ^ -1) - q, (-2 * E(3) + E(3) * q ^ -1 + 2 * E(3) * q) - E(3) * q ^ 2, ((3 * E(3, 2) - E(3, 2) * q ^ -1) - 3 * E(3, 2) * q) + E(3, 2) * q ^ 2, 0, ((3 - q ^ -1) - 3q) + q ^ 2, 0, E(3, 2) + (-2 * E(3, 2) + E(3)) * q + (E(3, 2) - 2 * E(3)) * q ^ 2 + E(3) * q ^ 3, 0, E(3) - E(3) * q, (-(E(3, 2)) + (2 * E(3, 2) - E(3)) * q + (-(E(3, 2)) + 2 * E(3)) * q ^ 2) - E(3) * q ^ 3, 0, (((3 * E(3) - E(3) * q ^ -1) - 4 * E(3) * q) + 3 * E(3) * q ^ 2) - E(3) * q ^ 3, 0, -1 + (-2 * E(3, 2) - E(3)) * q + E(3, 2) * q ^ 2, ((2 * E(3) - E(3) * q ^ -1) - 2 * E(3) * q) + E(3) * q ^ 2, ((2 * E(3) - E(3) * q ^ -1) - 2 * E(3) * q) + E(3) * q ^ 2, 0, (-(E(3, 2)) + (2 * E(3, 2) - E(3)) * q + (-(E(3, 2)) + 2 * E(3)) * q ^ 2) - E(3) * q ^ 3, 0, (-1 + q) - q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, (((-(E(3, 2)) + 2 * E(3)) - E(3) * q ^ -1) + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, ((E(3) - 3 * E(3) * q) + 3 * E(3) * q ^ 2) - E(3) * q ^ 3, (((-(E(3, 2)) + 2 * E(3)) - E(3) * q ^ -1) + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, (((4 * E(3, 2) - E(3, 2) * q ^ -1) - 6 * E(3, 2) * q) + 4 * E(3, 2) * q ^ 2) - E(3, 2) * q ^ 3, (((3 * E(3, 2) + E(3)) - E(3, 2) * q ^ -1) + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, -(E(3, 2)) + (2 * E(3, 2) + E(3)) * q + q ^ 2], [q - q ^ 2, -q + q ^ 2, q - q ^ 2, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, -q + q ^ 2, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q]], [[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q - q ^ 2, 0, q - q ^ 2, 0, -1 + q, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1 - q ^ -1, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 1 - q, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(-1 + 2q) - q ^ 2, (1 - 2q) + q ^ 2, (q - 2 * q ^ 2) + q ^ 3, ((1 - 3q) + 3 * q ^ 2) - q ^ 3, (-1 + 2q) - q ^ 2, 0, (1 - 2q) + q ^ 2, (1 - 2q) + q ^ 2, E(3) * q - E(3) * q ^ 2, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, -(E(3)) * q + E(3) * q ^ 2, (-q + 2 * q ^ 2) - q ^ 3, 0, 0, (q - 2 * q ^ 2) + q ^ 3, 0, -1 + (-3 * E(3, 2) - 2 * E(3)) * q + (3 * E(3, 2) + 2 * E(3)) * q ^ 2 + q ^ 3, 0, E(3) * q - E(3) * q ^ 2, (-1 + 2q) - q ^ 2, (-1 + 2q) - q ^ 2, 0, (q - 2 * q ^ 2) + q ^ 3, 0, -(E(3, 2)) * q + E(3, 2) * q ^ 2, 0, q - q ^ 2, 0, 0, 0, 0, q, 0, q, (-1 + 2q) - q ^ 2, (q - 2 * q ^ 2) + q ^ 3, (-1 + 2q) - q ^ 2, (1 + (3 * E(3, 2) + 2 * E(3)) * q + (-3 * E(3, 2) - 2 * E(3)) * q ^ 2) - q ^ 3, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, -(E(3)) * q + E(3) * q ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [-1 + q, 1 - q, -1 + q, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 1 - q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [-1 + q, 1 - q, -1 + q, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q - q ^ 2, 0, q - q ^ 2, 0, -1 + q, 0, 0, -q, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, (2 - q ^ -1) - q, (-1 + 2q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, -2 + q ^ -1 + q, 0, (2 - q ^ -1) - q, (2 - q ^ -1) - q, -(E(3)) + E(3) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, E(3) - E(3) * q, (1 - 2q) + q ^ 2, 0, 0, (-1 + 2q) - q ^ 2, 0, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, -(E(3)) + E(3) * q, -2 + q ^ -1 + q, -2 + q ^ -1 + q, 0, (-1 + 2q) - q ^ 2, 0, E(3, 2) - E(3, 2) * q, 0, -1 + q, 0, 0, 0, 0, -1 + q, 0, -1, -2 + q ^ -1 + q, (-1 + 2q) - q ^ 2, -2 + q ^ -1 + q, ((-3 * E(3, 2) - 2 * E(3)) - q ^ -1) + (3 * E(3, 2) + 2 * E(3)) * q + q ^ 2, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, E(3) - E(3) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, (-1 + 2q) - q ^ 2, (-q + 2 * q ^ 2) - q ^ 3, ((-3 + q ^ -1 + 4q) - 3 * q ^ 2) + q ^ 3, 0, (1 - 2q) + q ^ 2, (-1 + 2q) - q ^ 2, (-1 + 2q) - q ^ 2, 0, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, (3 - q ^ -1) + (3 * E(3, 2) + 4 * E(3)) * q + (-(E(3, 2)) - 2 * E(3)) * q ^ 2, ((-1 + 3q) - 3 * q ^ 2) + q ^ 3, -2 + q ^ -1 + q, -(E(3)) * q + E(3) * q ^ 2, (-q + 2 * q ^ 2) - q ^ 3, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, (1 - 2q) + q ^ 2, -(E(3)) * q + E(3) * q ^ 2, (1 - 2q) + q ^ 2, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, ((1 - 3q) + 3 * q ^ 2) - q ^ 3, (E(3, 2) + (-2 * E(3, 2) + E(3)) * q + (2 * E(3, 2) - E(3)) * q ^ 2) - E(3, 2) * q ^ 3, 0, E(3) - E(3) * q, -q + q ^ 2, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, -2 + q ^ -1 + q, q - q ^ 2, (1 - 2q) + q ^ 2, -q, -(E(3, 2)) * q + E(3, 2) * q ^ 2, 0, 0, (((3 - q ^ -1) - 4q) + 3 * q ^ 2) - q ^ 3, 0, -1 + (-3 * E(3, 2) - 2 * E(3)) * q + (3 * E(3, 2) + 2 * E(3)) * q ^ 2 + q ^ 3, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]]
                        end
                    f8(q) = begin
                            return q ^ 0 * [[[E(3) - E(3) * q, E(3) * q, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, E(3, 2) - E(3, 2) * q, 0, -(E(3)) + (-3 + root(-3)) // 2 * q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3), -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, 1 - q, -(E(3, 2)) + E(3, 2) * q ^ -1, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, E(3, 2) - E(3, 2) * q, (((5 + root(-3)) // 2 - q ^ -1) + (-2 - root(-3)) * q) - E(3, 2) * q ^ 2, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0], [((2 * E(3, 2) - E(3, 2) * q ^ -1) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, ((-3 + root(-3)) // 2 + q ^ -1) - E(3) * q, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, 0, 0, 0, (1 - q) + q ^ 2, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, (E(3) - root(-3) * q) - E(3, 2) * q ^ 2, (((3 * E(3, 2) - E(3, 2) * q ^ -1) - 4 * E(3, 2) * q) + 3 * E(3, 2) * q ^ 2) - E(3, 2) * q ^ 3, 0, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, 0, 0, 0, 0, 0, 0, 0], [0, 0, E(3, 2) - E(3, 2) * q ^ -1, 0, (1 - 2q) + q ^ 2, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0], [0, 0, (((-2 + root(-3)) - E(3) * q ^ -1) + (5 - root(-3)) // 2 * q) - q ^ 2, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0], [0, 0, (((-2 + root(-3)) - E(3) * q ^ -1) + (5 - root(-3)) // 2 * q) - q ^ 2, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0], [0, 0, (root(-3) - E(3) * q ^ -1) + E(3, 2) * q, 0, E(3, 2) + (2 + root(-3)) * q + (-5 - root(-3)) // 2 * q ^ 2 + q ^ 3, 0, 0, 0, 0, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, -1 + q, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, 0, 0, 0], [0, 0, ((-4 - q ^ -2) + 3 * q ^ -1 + 3q) - q ^ 2, (root(-3) - E(3) * q ^ -1) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q ^ -1 + E(3) * q, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, (2 - q ^ -1) - q, 0, 0], [((3 - q ^ -1) - 3q) + q ^ 2, (-1 + 2q) - q ^ 2, 0, 0, ((-1 + 3q) - 3 * q ^ 2) + q ^ 3, (1 - 2q) + q ^ 2, (E(3, 2) * q - 2 * E(3, 2) * q ^ 2) + E(3, 2) * q ^ 3, (E(3) - root(-3) * q) - E(3, 2) * q ^ 2, 0, -1 + q, 0, (E(3) - root(-3) * q) - E(3, 2) * q ^ 2, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, -1, 0, 0, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 1 - q, E(3, 2) * q - E(3, 2) * q ^ 2, 0, 0, 0, 0], [0, 0, ((-4 * E(3) - E(3) * q ^ -2) + 3 * E(3) * q ^ -1 + 3 * E(3) * q) - E(3) * q ^ 2, 0, ((-3 + root(-3)) // 2 + q ^ -1) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0], [-2 + q ^ -1 + q, 0, 0, (-3 + q ^ -1 + 3q) - q ^ 2, (((1 - 3 * root(-3)) + root(-3) * q ^ -1 + (-5 + 7 * root(-3)) // 2 * q) - 4 * E(3) * q ^ 2) + E(3) * q ^ 3, (((-1 + 2 * root(-3)) - root(-3) * q ^ -1) - 3 * E(3) * q) + E(3) * q ^ 2, (((-3 + root(-3)) // 2 + (7 - root(-3)) // 2 * q) - 3 * q ^ 2) + q ^ 3, (((4 * E(3, 2) - E(3, 2) * q ^ -1) + (2 + 3 * root(-3)) * q) - 2 * root(-3) * q ^ 2) + E(3) * q ^ 3, 0, (((5 + root(-3)) // 2 - q ^ -1) + (-2 - root(-3)) * q) - E(3, 2) * q ^ 2, 0, (((-2 - root(-3)) - E(3, 2) * q ^ -1) + (5 + root(-3)) // 2 * q) - q ^ 2, ((-3 + root(-3)) // 2 + q ^ -1) - E(3) * q, ((3 - q ^ -1) - 3q) + q ^ 2, -2 + q ^ -1 + q, (-1 + 2q) - q ^ 2, 0, -1 + q, -1, 0, (-(E(3, 2)) - root(-3) * q) + E(3) * q ^ 2, 0, (-3 * E(3) + E(3) * q ^ -1 + 3 * E(3) * q) - E(3) * q ^ 2, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, (-1 + 2q) - q ^ 2, (1 - 2q) + q ^ 2, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, 0, -(E(3)) + E(3) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0], [(-4 + 2 * q ^ -1 + 3q) - q ^ 2, (1 - 2q) + q ^ 2, 0, -2 + q ^ -1 + q, 0, 0, 0, (((-1 - 2 * root(-3)) - E(3, 2) * q ^ -1) + (1 + 5 * root(-3)) // 2 * q) - root(-3) * q ^ 2, 0, ((5 + root(-3)) // 2 - q ^ -1) + (-3 - root(-3)) // 2 * q, 0, 0, 0, (2 - q ^ -1) - q, 0, -1 + q, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 1 - q, 0, 0, -(E(3))], [0, 0, 0, 0, ((2 - root(-3)) + E(3) * q ^ -1 + (-7 + root(-3)) // 2 * q + 3 * q ^ 2) - q ^ 3, 0, 0, 0, 0, 0, 0, (((-1 - 2 * root(-3)) - E(3, 2) * q ^ -1) + (1 + 5 * root(-3)) // 2 * q) - root(-3) * q ^ 2, -2 + q ^ -1 + q, 0, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, 0, 0, -1 + q, 0, -1 + q ^ -1, 0, E(3) - E(3) * q, (1 - 2q) + q ^ 2, 0, 0, 0, -1 + q, 0, -(E(3, 2)), 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [5 + 2 * root(-3) + (-3 - root(-3)) // 2 * q ^ -1 + (-6 - 3 * root(-3)) * q + (3 + 2 * root(-3)) * q ^ 2 + E(3, 2) * q ^ 3, ((E(3, 2) - 3 * E(3, 2) * q) + 3 * E(3, 2) * q ^ 2) - E(3, 2) * q ^ 3, 0, (((4 - q ^ -1) - 6q) + 4 * q ^ 2) - q ^ 3, (((-3 + 7 * root(-3)) // 2 - root(-3) * q ^ -1) + (11 - 9 * root(-3)) // 2 * q + (-15 + 5 * root(-3)) // 2 * q ^ 2 + (9 - root(-3)) // 2 * q ^ 3) - q ^ 4, ((3 - 5 * root(-3)) // 2 + root(-3) * q ^ -1 + (-4 + 2 * root(-3)) * q + (7 - root(-3)) // 2 * q ^ 2) - q ^ 3, ((3 - root(-3)) // 2 + (-9 + root(-3)) // 2 * q + (5 + root(-3)) * q ^ 2 + (-5 - 3 * root(-3)) // 2 * q ^ 3) - E(3, 2) * q ^ 4, (3 + 5 * root(-3)) // 2 + E(3, 2) * q ^ -1 + (-3 - 11 * root(-3)) // 2 * q + 6 * root(-3) * q ^ 2 + (1 - 3 * root(-3)) * q ^ 3 + E(3) * q ^ 4, (1 + (-3 + root(-3)) // 2 * q) - E(3) * q ^ 2, ((-4 - root(-3)) + q ^ -1 + (11 + 5 * root(-3)) // 2 * q + (-3 - 2 * root(-3)) * q ^ 2) - E(3, 2) * q ^ 3, -(E(3)) + E(3) * q, ((3 + 5 * root(-3)) // 2 + E(3, 2) * q ^ -1 + (-3 - 9 * root(-3)) // 2 * q + (1 + 7 * root(-3)) // 2 * q ^ 2) - root(-3) * q ^ 3, ((3 - q ^ -1) - 3q) + q ^ 2, -3 + q ^ -1 + (9 + root(-3)) // 2 * q + (-7 - root(-3)) // 2 * q ^ 2 + q ^ 3, (-3 * E(3, 2) + E(3, 2) * q ^ -1 + 3 * E(3, 2) * q) - E(3, 2) * q ^ 2, (1 + (-5 - root(-3)) // 2 * q + (5 + root(-3)) // 2 * q ^ 2) - q ^ 3, 0, (1 - 2q) + q ^ 2, 0, (2 - q ^ -1) - q, E(3, 2) + (1 + 3 * root(-3)) // 2 * q + (1 - 3 * root(-3)) // 2 * q ^ 2 + E(3) * q ^ 3, (E(3, 2) + root(-3) * q) - E(3) * q ^ 2, (-3 + root(-3)) // 2 + (4 - root(-3)) * q + (-7 + root(-3)) // 2 * q ^ 2 + q ^ 3, -1, (2 - root(-3)) + E(3) * q ^ -1 + (-5 + root(-3)) // 2 * q + q ^ 2, 1 + (-5 - root(-3)) // 2 * q + (2 + root(-3)) * q ^ 2 + E(3, 2) * q ^ 3, (q - 2 * q ^ 2) + q ^ 3, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, E(3, 2) - E(3, 2) * q, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, -1 + q, 0, 0], [0, 0, (((-4 - q ^ -2) + 3 * q ^ -1 + 4q) - 3 * q ^ 2) + q ^ 3, 0, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, 0, 0, 0, 0, 0, 0, (((3 * E(3) - E(3) * q ^ -1) - 4 * E(3) * q) + 3 * E(3) * q ^ 2) - E(3) * q ^ 3, ((2 - q ^ -1) - 2q) + q ^ 2, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, (1 - 2q) + q ^ 2, 0, 0, 0, (-(E(3)) + E(3) * q) - E(3) * q ^ 2, E(3, 2) - E(3, 2) * q, 0, 0, 0, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, ((2 - q ^ -1) - 2q) + q ^ 2, E(3, 2) - E(3, 2) * q, 0], [((-5 + 2 * q ^ -1 + 5q) - 3 * q ^ 2) + q ^ 3, ((2 - 4q) + 3 * q ^ 2) - q ^ 3, (1 - 2q) + q ^ 2, ((-(root(-3)) - E(3, 2) * q ^ -1) + (-5 + root(-3)) // 2 * q + 3 * q ^ 2) - q ^ 3, 0, 0, 0, ((-1 - 2 * root(-3)) - E(3, 2) * q ^ -1) + (1 + 3 * root(-3)) * q + (-1 - 5 * root(-3)) // 2 * q ^ 2 + root(-3) * q ^ 3, 0, ((7 + root(-3)) // 2 - q ^ -1) + (-4 - root(-3)) * q + (3 + root(-3)) // 2 * q ^ 2, 0, 0, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, (-(E(3, 2)) + E(3, 2) * q) - E(3, 2) * q ^ 2, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, E(3, 2) - E(3, 2) * q]], [[0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, q - q ^ 2, 0, 0, -(E(3)) + E(3) * q, -1 + q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + q, 0, 0, 0, E(3) - E(3) * q, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q - q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)), 0], [0, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [((-2 + root(-3)) - E(3) * q ^ -1) + (4 - root(-3)) * q + (-7 + root(-3)) // 2 * q ^ 2 + q ^ 3, (-(E(3)) + (-2 + root(-3)) * q + (5 - root(-3)) // 2 * q ^ 2) - q ^ 3, 0, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, 1 + 2 * root(-3) + E(3, 2) * q ^ -1 + (-1 - 5 * root(-3)) // 2 * q + root(-3) * q ^ 2, 0, -(E(3)) + (-2 + root(-3)) * q + (3 - root(-3)) // 2 * q ^ 2, 0, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, -1 + q, E(3) * q, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, 1 - q, 0, (1 - 2q) + q ^ 2, (1 - 2q) + q ^ 2, (-5 - root(-3)) // 2 + q ^ -1 + (3 + root(-3)) // 2 * q, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, 0, 0, 1 - q, 0, E(3, 2), 0, 0, 0, 0, 0], [0, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, E(3), 0, 0], [0, 0, q - q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q - E(3, 2) * q ^ 2, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q]], [[0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, q, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, E(3) * q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -(E(3)) + E(3) * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, (1 - 2q) + q ^ 2, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, (-5 - root(-3)) // 2 + q ^ -1 + (2 + root(-3)) * q + E(3, 2) * q ^ 2, -q + q ^ 2, 0, (1 - 2q) + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, -q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, ((3 + root(-3)) // 2 - q ^ -1) + E(3, 2) * q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3))], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q ^ -1, 0, -(E(3, 2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, -(E(3)) * q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1 + q, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q + q ^ 2, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, E(3) * q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, -1 + q, (2 - root(-3)) + E(3) * q ^ -1 + (-5 + root(-3)) // 2 * q + q ^ 2, 0, -1 + q, 0, -1 + q, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q ^ -1, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, (-(E(3, 2)) - root(-3) * q) + E(3) * q ^ 2, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q ^ -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, E(3, 2) * q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, (((-2 + root(-3)) - E(3) * q ^ -1) + (5 - root(-3)) // 2 * q) - q ^ 2, 0, 0, 0, 1 - q, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, q - q ^ 2, 0, 0, -(E(3)) + E(3) * q, -1 + q, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, -1 + q, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, q, 0, -(E(3, 2)) * q + E(3, 2) * q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q - q ^ 2, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q - q ^ 2, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, 0, (-1 + 2q) - q ^ 2, ((5 - 3 * root(-3)) // 2 + E(3) * q ^ -1 + (-9 + 3 * root(-3)) // 2 * q + (7 - root(-3)) // 2 * q ^ 2) - q ^ 3, 0, q - q ^ 2, 0, (-1 + 2q) - q ^ 2, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1 - 2q) + q ^ 2, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q - q ^ 2, 0, 0, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, E(3, 2) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, 1 - q, -(E(3, 2)) + E(3, 2) * q ^ -1, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3), 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, q - q ^ 2, 0, 0, -(E(3)) + E(3) * q, -1 + q, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, -(E(3)), 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q + q ^ 2, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, -(E(3, 2)) * q + E(3, 2) * q ^ 2, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [((4 - 2 * q ^ -1) - 3q) + q ^ 2, (-1 + 2q) - q ^ 2, 0, (2 - q ^ -1) - q, 0, 0, 0, 1 + 2 * root(-3) + E(3, 2) * q ^ -1 + (-1 - 5 * root(-3)) // 2 * q + root(-3) * q ^ 2, 0, (-5 - root(-3)) // 2 + q ^ -1 + (3 + root(-3)) // 2 * q, 0, 0, 0, -2 + q ^ -1 + q, 0, 1 - q, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, -1 + q, 0, 0, E(3)], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0], [((-1 - 3 * root(-3)) // 2 - E(3, 2) * q ^ -1) + (-1 + 5 * root(-3)) // 2 * q + (1 - 2 * root(-3)) * q ^ 2 + E(3) * q ^ 3, (-(E(3, 2)) + (-1 - 3 * root(-3)) // 2 * q + (-1 + 3 * root(-3)) // 2 * q ^ 2) - E(3) * q ^ 3, 0, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, (-7 - root(-3)) // 2 + q ^ -1 + (4 + root(-3)) * q + (-3 - root(-3)) // 2 * q ^ 2, 0, -(E(3, 2)) + (-1 - 3 * root(-3)) // 2 * q + root(-3) * q ^ 2, 0, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, -(E(3)) + E(3) * q, E(3, 2) * q, 0, 0, 0, -1 + q]]]
                        end
                    r = [f1(x), f2(x), f3(x, E(3)), f4(x, E(3)), f4(x, E(3, 2)), f3(x, E(3, 2)), [[[-1]], [[-1]], [[-1]], [[-1]], [[-1]]], f8(x), f9(x), f8(1 // x) * -x, f11(x, E(3)), f12(x, E(3)), f13(x, E(3)), f13(x, E(3, 2)), f12(x, E(3, 2)), f11(x, E(3, 2)), f17(x), f1(1 // x) * -x, f13(1 // x, E(3, 2)) * -x, f20(x, E(3)), f20(x, E(3, 2)), f13(1 // x, E(3)) * -x, f23(x), f2(1 // x) * -x, f11(1 // x, E(3, 2)) * -x, f12(1 // x, E(3)) * -x, f12(1 // x, E(3, 2)) * -x, f11(1 // x, E(3)) * -x, f29(x), f4(1 // x, E(3, 2)) * -x, f4(1 // x, E(3)) * -x, f23(1 // x) * -x, f3(1 // x, E(3, 2)) * -x, f3(1 // x, E(3)) * -x, f17(1 // x) * -x, [[[x]], [[x]], [[x]], [[x]], [[x]]]]
                    return r[i]
                end
            return -((para[2])[2]) * m335(i)
        elseif [p, q, r] == [4, 4, 3]
            x = -((para[2])[1]) // (para[2])[2]
            r = x ^ 0 * [[[[x, -1, -1, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, -1 + x, -1, x, 0], [0, 1, -1, -1, 0, x]], [[-1, 0, 0, 0, 0, 0], [-x, x, 0, 0, 0, x], [0, 0, 0, 0, -x, 0], [0, 0, 0, x, 0, 0], [0, 0, -1, 0, -1 + x, 0], [0, 0, 0, 0, 0, -1]], [[x, -1, 0, 0, 0, -1], [0, -1, 0, 0, 0, 0], [0, 0, x, 0, 1, -1], [0, -x, 0, x, -1, 1 - x], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]]], [[[-1, 0, 0], [0, 0, 1], [0, x, -1 + x]], [[x, 0, 0], [-1, -1, 0], [1, 0, -1]], [[0, 1, 0], [x, -1 + x, 0], [0, 0, -1]]], [[[x, 0, 0], [x, -1, 0], [x, 0, -1]], [[-1, 2, 0], [0, x, 0], [0, (-(E(4)) + 1) * x, -1]], [[-1, 0, 1], [0, -1, (E(4) + 1) // 2], [0, 0, x]]], [[[x, 0, 0], [x, -1, 0], [x, 0, -1]], [[-1, 2, 0], [0, x, 0], [0, (E(4) + 1) * x, -1]], [[-1, 0, 1], [0, -1, (-(E(4)) + 1) // 2], [0, 0, x]]], [[[-1]], [[-1]], [[-1]]], [[[x, -1, 0], [0, -1, 0], [0, -1, x]], [[-1 + x, 0, 1], [0, x, 0], [x, 0, 0]], [[0, x, 0], [1, -1 + x, 0], [0, 0, x]]], [[[-1, 0, 0], [-1, x, 0], [-1, 0, x]], [[x, -2x, 0], [0, -1, 0], [0, -(E(4)) - 1, x]], [[x, 0, -x], [0, x, (E(4) - 1) // 2 * x], [0, 0, -1]]], [[[-1, 0, 0], [-1, x, 0], [-1, 0, x]], [[x, -2x, 0], [0, -1, 0], [0, E(4) - 1, x]], [[x, 0, -x], [0, x, (-(E(4)) - 1) // 2 * x], [0, 0, -1]]], [[[-1, 0], [-1, x]], [[-1, 0], [-1, x]], [[x, -x], [0, -1]]], [[[x]], [[x]], [[x]]]]
            return -((para[2])[2]) * r[i]
        else
            S = (((chevieget(:imp, :CharInfo))(p, q, r))[:charparams])[i]
            p1rRep() = begin
                    local Q, pos, ct
                    if r > 1
                        Q = -((para[2])[1]) // (para[2])[2]
                    else
                        Q = 0
                    end
                    pos(t, i) = begin
                            local j, k, l
                            for j = 1:length(t)
                                for k = 1:length(t[j])
                                    l = Position((t[j])[k], i)
                                    if l != false
                                        return [j, k, l]
                                    end
                                end
                            end
                        end
                    ct = (p->begin
                                (para[1])[p[1]] * ((Q * 1) // 1) ^ (p[3] - p[2])
                            end)
                    T = tableaux(S)
                    return Concatenation([DiagonalMat(map((S->begin
                                                    ct(pos(S, 1))
                                                end), T))], map((i->begin
                                            map(function (j,)
                                                    local S, v, a, b, p, tll
                                                    S = T[j]
                                                    a = pos(S, i)
                                                    b = pos(S, i - 1)
                                                    S = map((a->begin
                                                                    map(copy, a)
                                                                end), S)
                                                    ((S[a[1]])[a[2]])[a[3]] = i - 1
                                                    ((S[b[1]])[b[2]])[b[3]] = i
                                                    if (para[2])[1] == -((para[2])[2])
                                                        if a[1] == b[1]
                                                            tll = (para[2])[1] // (((a[3] + b[2]) - a[2]) - b[3])
                                                        else
                                                            tll = 0
                                                        end
                                                    else
                                                        tll = Sum(para[2]) // (1 - ct(b) // ct(a))
                                                    end
                                                    v = fill(0, max(0, (1 + length(T)) - 1)) // 1 * tll
                                                    v[j] = tll
                                                    p = Position(T, S)
                                                    if p != false
                                                        v[p] = tll - (para[2])[2]
                                                    end
                                                    return v
                                                end, 1:length(T))
                                        end), 2:r)) * Product(para, Product) ^ 0
                end
            if q == 1
                return p1rRep()
            elseif p == q
                para = [map((i->begin
                                    E(p, i)
                                end), 0:p - 1), para[1]]
            else
                e = div(p, q)
                if mod(q, 2) == 0 && r == 2
                    S = (((chevieget(:imp, :CharInfo))(p, q, r))[:malle])[i]
                    if S[1] == 1
                        return [[[(para[1])[1 + mod(S[4] - 1, e)]]], [[(para[2])[S[2]]]], [[(para[3])[S[3]]]]]
                    else
                        Y = para[2]
                        T = para[3]
                        if q > 2
                            X = map((y->begin
                                            GetRoot(y, div(q, 2))
                                        end), para[1])
                            X = Concatenation(map((i->begin
                                                E(div(q, 2), i) * X
                                            end), 0:div(q, 2) - 1))
                        else
                            X = para[1]
                        end
                        X = X[S[[3, 4]]]
                        v = S[2] * GetRoot(Product(X) * Product(Y) * Product(T) * E(div(p, 2), (2 - S[3]) - S[4]), 2) * E(p, (S[3] + S[4]) - 2)
                        d = 1 + Sum(X) * 0 + Sum(Y) * 0 + Sum(T) * 0
                        return [(d * [[X[1], Sum(Y, (y->begin
                                                                1 // y
                                                            end)) - X[2] // v * Sum(T)], [0, X[2]]]) ^ div(q, 2), [[Sum(Y), 1 // X[1]], [-(Product(Y)) * X[1], 0]], [[0, -(Product(T)) // v], [v, Sum(T)]]]
                    end
                end
                if para[2] != para[3]
                    error("should not happen")
                end
                if para[1] == map((i->begin
                                    E(e, i - 1)
                                end), 1:e)
                    para = [map((i->begin
                                        E(p, i)
                                    end), 0:p - 1), para[2]]
                else
                    para = [Concatenation(TransposedMat(map((i->begin
                                                map((j->begin
                                                                E(q, j)
                                                            end), 0:q - 1) * GetRoot(i, q)
                                            end), para[1]))), para[2]]
                end
            end
            extra = false
            if IsInt(S[length(S)])
                extra = E(S[length(S) - 1], S[length(S)])
                d = length(S) - 2
                S = fullsymbol(S)
            end
            v = p1rRep()
            if p == q
                v = Concatenation([v[2] ^ v[1]], v[2:length(v)])
            elseif q > 1
                v = Concatenation([v[1] ^ q, v[2] ^ v[1]], v[2:length(v)])
            end
            if extra != false
                m = PermListList(T, map((S->begin
                                    S[Concatenation(d + 1:p, 1:d)]
                                end), T))
                m = orbits(m, 1:length(T))
                l = map((i->begin
                                extra ^ i
                            end), 0:-1 - 0:1 - p // d)
                m1 = map((x->begin
                                x[1]
                            end), m)
                return map((x->begin
                                map((c->begin
                                            l * map((y->begin
                                                            y[m1]
                                                        end), x[c])
                                        end), m)
                            end), v)
            else
                return v
            end
        end
    end)
chevieset(:imp, :Representation, function (p, q, r, i)
        local o
        o = (chevieget(:imp, :EigenvaluesGeneratingReflections))(p, q, r)
        o = map(denominator, o)
        return (chevieget(:imp, :HeckeRepresentation))(p, q, r, map((x->begin
                            map((i->begin
                                        E(x, i)
                                    end), 0:x - 1)
                        end), o), [], i)
    end)
chevieset(:imp, :UnipotentCharacters, function (p, q, r)
        local uc, cusp, f, l, ci, seteig, s, extra, addextra
        if !(q in [1, p])
            return false
        end
        uc = Dict{Symbol, Any}(:charSymbols => (chevieget(:imp, :CharSymbols))(p, q, r))
        uc[:a] = map(valuation_gendeg_symbol, uc[:charSymbols])
        uc[:A] = map(degree_gendeg_symbol, uc[:charSymbols])
        ci = (chevieget(:imp, :CharInfo))(p, q, r)
        if q == 1
            cusp = gapSet(map((S->begin
                                map(length, S) - Minimum(map(length, S))
                            end), uc[:charSymbols]))
            cusp = map((x->begin
                            map((y->begin
                                        0:y - 1
                                    end), x)
                        end), cusp)
            SortBy(cusp, ranksymbol)
            uc[:harishChandra] = map(function (c,)
                        local cr, res
                        cr = ranksymbol(c)
                        res = Dict{Symbol, Any}(:levi => 1:cr)
                        if cr < r
                            res[:parameterExponents] = [map(length, c)]
                        else
                            res[:parameterExponents] = []
                        end
                        res[:parameterExponents] = Append(res[:parameterExponents], fill(0, max(0, (1 + r) - (2 + cr))) + 1)
                        if r == cr
                            res[:relativeType] = Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0)
                        else
                            res[:relativeType] = Dict{Symbol, Any}(:series => "ST", :indices => 1 + cr:r, :rank => r - cr, :p => p, :q => 1)
                        end
                        res[:eigenvalue] = E(24, -2 * (p ^ 2 - 1) * div(Sum(c, length), p)) * E(2p, Sum(0:p - 1, (i->begin
                                                -((i ^ 2 + p * i)) * length(c[i + 1])
                                            end)))
                        res[:charNumbers] = map((x->begin
                                        Position(uc[:charSymbols], symbol_partition_tuple(x, map(length, c)))
                                    end), map((x->begin
                                            map(partβ, x)
                                        end), ((chevieget(:imp, :CharSymbols))(p, 1, r - cr))[1:length(partition_tuples(r - cr, p))]))
                        res[:cuspidalName] = ImprimitiveCuspidalName(c)
                        return res
                    end, cusp)
            uc[:b] = uc[:a] * 0
            uc[:B] = uc[:a] * 0
            (uc[:b])[((uc[:harishChandra])[1])[:charNumbers]] = ci[:b]
            (uc[:B])[((uc[:harishChandra])[1])[:charNumbers]] = ci[:B]
            uc[:families] = map((y->begin
                            MakeFamilyImprimitive(y, uc)
                        end), CollectBy(uc[:charSymbols], (x->begin
                                tally(Concatenation(x))
                            end)))
            SortBy(uc[:families], (x->begin
                        x[:charNumbers]
                    end))
            if r == 1
                l = map(function (S,)
                            local p
                            p = Position(S, [])
                            if p == false
                                return 1
                            else
                                return (-1) ^ p
                            end
                        end, (uc[:charSymbols])[((uc[:families])[2])[:charNumbers]])
                ((uc[:families])[2])[:fourierMat] = ((uc[:families])[2])[:fourierMat] ^ DiagonalMat(l)
                uc[:cyclicparam] = map(function (s,)
                            if count((x->begin
                                                x == 1
                                            end), Flat(s)) == 1
                                return [1]
                            else
                                s = deepcopy(s)
                                l = PositionProperty(s, (p->begin
                                                1 in p
                                            end))
                                s[l] = []
                                return [PositionProperty(s, (p->begin
                                                        1 in p
                                                    end)) - 1, l - 1]
                            end
                        end, uc[:charSymbols])
            elseif r == 2 && p == 3
                ((uc[:families])[4])[:fourierMat] = ((uc[:families])[4])[:fourierMat] ^ DiagonalMat(-1, 1, 1)
                ((uc[:families])[1])[:fourierMat] = ((uc[:families])[1])[:fourierMat] ^ DiagonalMat(1, -1, -1, 1, 1, 1, 1, 1, 1)
            end
            return uc
        elseif p == q
            uc[:families] = []
            for f = CollectBy(1:length(uc[:charSymbols]), (i->begin
                                tally(Concatenation(fullsymbol((uc[:charSymbols])[i])))
                            end))
                if length(gapSet(map(fullsymbol, (uc[:charSymbols])[f]))) > 1
                    push!(uc[:families], Dict{Symbol, Any}(:charNumbers => f))
                else
                    uc[:families] = Append(uc[:families], map((x->begin
                                        Family("C1", [x])
                                    end), f))
                end
            end
            SortBy(uc[:families], (x->begin
                        x[:charNumbers]
                    end))
            uc[:harishChandra] = map((l->begin
                            Dict{Symbol, Any}(:charNumbers => l)
                        end), CollectBy(1:length(uc[:charSymbols]), function (i,)
                            local s, l
                            s = fullsymbol((uc[:charSymbols])[i])
                            l = map(length, s)
                            return [Sum(s, (x->begin
                                                Sum(partβ(x))
                                            end)), l - Minimum(l)]
                        end))
            SortBy(uc[:harishChandra], (x->begin
                        x[:charNumbers]
                    end))
            extra = []
            for f = uc[:harishChandra]
                addextra = false
                s = fullsymbol((uc[:charSymbols])[(f[:charNumbers])[1]])
                l = r - Sum(s, (x->begin
                                    Sum(partβ(x))
                                end))
                f[:levi] = 1:l
                s = map(length, s)
                s = s - Minimum(s)
                f[:eigenvalue] = E(24, -2 * (p ^ 2 - 1) * div(Sum(s), p)) * E(2p, Sum(0:p - 1, (i->begin
                                        -((i ^ 2 + p * i)) * s[i + 1]
                                    end)))
                if l == r
                    f[:relativeType] = Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0)
                    f[:parameterExponents] = []
                    if length(f[:charNumbers]) == 2
                        addextra = true
                    end
                elseif l == 0
                    f[:relativeType] = Dict{Symbol, Any}(:series => "ST", :indices => 1:r, :rank => r, :p => p, :q => q)
                    f[:parameterExponents] = fill(0, max(0, (1 + r) - 1)) + 1
                else
                    f[:relativeType] = Dict{Symbol, Any}(:series => "ST", :indices => l + 1:r, :rank => r - l, :p => p, :q => 1)
                    f[:parameterExponents] = Concatenation([s], fill(0, max(0, (1 + ((r - l) - 1)) - 1)) + 1)
                end
                s = map((x->begin
                                0:x - 1
                            end), s)
                f[:cuspidalName] = ImprimitiveCuspidalName(s)
                if addextra
                    s = deepcopy(f[:charNumbers])
                    f[:charNumbers] = s[[1]]
                    f = deepcopy(f)
                    f[:charNumbers] = s[[2]]
                    push!(f[:cuspidalName], '2')
                    push!(extra, f)
                end
            end
            uc[:harishChandra] = Append(uc[:harishChandra], extra)
            for f = uc[:families]
                f[:eigenvalues] = map((i->begin
                                (First(uc[:harishChandra], (s->(i in s[:charNumbers];))))[:eigenvalue]
                            end), f[:charNumbers])
            end
            uc[:b] = fill(0, max(0, (1 + length(uc[:charSymbols])) - 1))
            uc[:B] = fill(0, max(0, (1 + length(uc[:charSymbols])) - 1))
            (uc[:b])[((uc[:harishChandra])[1])[:charNumbers]] = ci[:b]
            (uc[:B])[((uc[:harishChandra])[1])[:charNumbers]] = ci[:B]
            if [p, q, r] == [3, 3, 3]
                (uc[:families])[6] = Family(conj(((CHEVIE[:families])[:X])(3)), [8, 7, 11], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[4] = Family(((CHEVIE[:families])[:X])(3), [4, 5, 12], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                uc[:curtis] = [1, 2, 3, 7, 8, 10, 4, 5, 9, 6, -12, -11]
            elseif [p, q, r] == [3, 3, 4]
                uc[:families] = map((x->begin
                                MakeFamilyImprimitive((uc[:charSymbols])[x[:charNumbers]], uc)
                            end), uc[:families])
            elseif [p, q, r] == [3, 3, 5]
                (uc[:families])[3] = Family(((CHEVIE[:families])[:X])(3), [3, 6, 51], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[4] = Family(((CHEVIE[:families])[:X])(3), [4, 5, 54], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[6] = Family(((CHEVIE[:families])[:QZ])(3), [9, 10, 8, 21, 44, 46, 20, 49, 45], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, 1, 1, -1, -1]))
                (uc[:families])[7] = Family(((CHEVIE[:families])[:QZ])(3), [23, 11, 16, 12, 42, 50, 15, 48, 40], Dict{Symbol, Any}(:signs => [1, -1, -1, 1, 1, 1, 1, -1, -1], :special => 4, :cospecial => 7))
                (uc[:families])[8] = Family(conj(((CHEVIE[:families])[:X])(3)), [14, 13, 41], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[11] = Family(((CHEVIE[:families])[:X])(3), [19, 22, 47], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[13] = Family(((CHEVIE[:families])[:QZ])(3), [32, 27, 26, 28, 38, 53, 25, 52, 39], Dict{Symbol, Any}(:signs => [1, 1, 1, -1, -1, 1, -1, -1, 1], :special => 3, :cospecial => 2))
                (uc[:families])[15] = Family(conj(((CHEVIE[:families])[:X])(3)), [31, 30, 37], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[16] = Family(conj(((CHEVIE[:families])[:X])(3)), [34, 33, 43], Dict{Symbol, Any}(:signs => [1, 1, -1]))
            elseif [p, q, r] == [4, 4, 3]
                (uc[:families])[2] = Family(((CHEVIE[:families])[:X])(4), [3, 2, 4, 14, 16, 13], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, -1, -1]))
                (uc[:families])[4] = Family(conj(((CHEVIE[:families])[:X])(4)), [8, 6, 7, 12, 15, 11], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, -1, -1]))
                uc[:curtis] = [1, 6, 7, 8, 10, 2, 3, 4, 9, 5, 14, 13, 12, 11, -16, -15]
            elseif [p, q, r] == [4, 4, 4]
                (uc[:families])[5] = Family(((CHEVIE[:families])[:X])(4), [5, 8, 9, 46, 53, 47], Dict{Symbol, Any}(:signs => [1, 1, 1, -1, -1, 1]))
                (uc[:families])[6] = Family("C2", [12, 7, 6, 42])
                (uc[:families])[7] = Family(((CHEVIE[:families])[:X])(4), [13, 10, 11, 41, 55, 43], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, -1], :special => 3, :cospecial => 1))
                (uc[:families])[9] = Family(((CHEVIE[:families])[:QZ])(4), [18, 21, 28, 22, 23, 49, 39, 54, 56, 40, 15, 36, 19, 52, 37, 51], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, -1], :special => 2, :cospecial => 4))
                (uc[:families])[10] = Family(conj(((CHEVIE[:families])[:X])(4)), [16, 17, 20, 38, 50, 34], Dict{Symbol, Any}(:signs => [1, 1, 1, -1, 1, 1], :special => 3, :cospecial => 1))
                (uc[:families])[12] = Family("C2", [27, 26, 25, 35])
                (uc[:families])[13] = Family(conj(((CHEVIE[:families])[:X])(4)), [30, 29, 31, 44, 48, 45], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, -1], :special => 3, :cospecial => 1))
            else
                uc[:families] = map((x->begin
                                MakeFamilyImprimitive((uc[:charSymbols])[x[:charNumbers]], uc)
                            end), uc[:families])
            end
            return uc
        end
    end)
chevieset(:imp, :InitHeckeBasis, function (p, q, r, H)
        if q != 1 || r != 1
            error("implemented only for G(d,1,1)")
        end
        H[:mul] = function (x, y)
                local H, res, ops, W, temp, i, xi, temp1, j, e, pol, d
                H = hecke(y)
                if !(IsRec(x)) || (!(haskey(x, :hecke)) || !(haskey(x, :elm)))
                    if x == x * 0
                        return HeckeElt(H, y[:basis], [], [])
                    else
                        return HeckeElt(H, y[:basis], y[:elm], y[:coeff] * (x * H[:unit]))
                    end
                end
                if !(IsIdentical(H, hecke(x)))
                    error("not elements of the same algebra")
                end
                ops = H[:operations]
                pol = Coefficients(Product((H[:parameter])[1], (u->begin
                                    Mvp("xxx") - u
                                end)), "xxx")
                d = length(pol) - 1
                if x[:basis] != y[:basis]
                    return (basis(H, "T"))(x) * (basis(H, "T"))(y)
                elseif x[:basis] == "T"
                    W = Group(H)
                    res = HeckeElt(H, x[:basis], [], [])
                    for i = 1:length(x[:elm])
                        temp = (x[:coeff])[i] * y
                        xi = (x[:elm])[i]
                        for i = 1:length(xi)
                            temp1 = HeckeElt(H, x[:basis], [], [])
                            for j = 1:length(temp[:elm])
                                e = length((temp[:elm])[j])
                                if e + 1 < d
                                    push!(temp1[:elm], map((i->begin
                                                    1
                                                end), 1:e + 1))
                                    push!(temp1[:coeff], (temp[:coeff])[j])
                                else
                                    temp1[:elm] = Append(temp1[:elm], map((i->begin
                                                        fill(0, max(0, (1 + i) - 1)) + 1
                                                    end), 0:d - 1))
                                    temp1[:coeff] = Append(temp1[:coeff], -(pol[1:d]) * (temp[:coeff])[j])
                                end
                            end
                            temp = temp1
                            CollectCoefficients(temp)
                        end
                        res = res + temp
                    end
                    return res
                else
                    return (basis(H, x[:basis]))((basis(H, "T"))(x) * (basis(H, "T"))(y))
                end
            end
        H[:inverse] = function (h,)
                local H, d, pol
                if length(h[:elm]) != 1
                    error("inverse implemented only for single T_w")
                end
                H = hecke(h)
                pol = Coefficients(Product((H[:parameter])[1], (u->begin
                                    Mvp("xxx") - u
                                end)), "xxx")
                d = length(pol) - 1
                return (h[:coeff])[1] ^ -1 * (basis(H, "T"))(map((i->begin
                                            fill(0, max(0, (1 + i) - 1)) + 1
                                        end), 0:d - 1), -(pol[2:d + 1]) // pol[1]) ^ length((h[:elm])[1])
            end
        return true
    end)
