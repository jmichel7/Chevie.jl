
chevieset(:imp, :PrintDiagram, function (arg...,)
        local p, q, r, indices, j, indent, title, g
        p = arg[1]
        q = arg[2]
        r = arg[3]
        indices = arg[4]
        title = arg[5]
        print(title, " ")
        indent = length(title) + 1
        g = (i->begin
                    pad("", indent - i)
                end)
        if q == 1
            print(indices[1], "(", p, ")")
            if length(indices) > 1
                print("==")
            end
            print(Join(indices[2:length(indices)], "--"), "\n")
        elseif p == q
            print(indices[1], "\n", g(0), "|\\\n")
            if p != 3
                print(pad(p, indent))
            else
                print(g(0))
            end
            print("|==", indices[3])
            for j = 4:r
                print(" - ", indices[j])
            end
            print("\n")
            print(g(0), "|/\n", g(0), indices[2], "\n")
        elseif q == 2
            print(indices[2], "\n", g(2), "/3|")
            if r >= 3
                print("\\")
            end
            print("\n")
            if p // q > 2
                print(g(length(string(p // q)) + 5), "(", p // q, ")")
            else
                print(g(3))
            end
            print(indices[1], "  | ")
            for j = 3:r
                print(indices[j + 1])
                if j != r
                    print("-")
                end
            end
            print("\n", g(2), "\\ |")
            if r >= 3
                print("/")
            end
            print("\n", pad(indices[3], indent + 1), "   ", joindigits(indices[[1, 2, 3]]), "==", joindigits(indices[[2, 3, 1]]), "==", joindigits(indices[[3, 1, 2]]), "\n")
        else
            print(indices[2], "\n", g(2), "/", q + 1, " ")
            if r >= 3
                print("\\")
            end
            print("\n")
            if p // q > 2
                print(g(length(SPrint(p // q)) + 5), "(", p // q, ")")
            else
                print(g(3))
            end
            print(indices[1], "   ")
            if r >= 3
                print("==")
            end
            for j = 3:r
                print(indices[j + 1])
                if j != r
                    print("-")
                end
            end
            print("\n", g(2), "\\  ")
            if r >= 3
                print("/")
            end
            print("\n", pad(indices[3], indent + 1))
            j = (chevieget(:imp, :BraidRelations))(p, q, r)
            for g = 1:Minimum(3, r)
                print("   ", joindigits(indices[(j[g])[1]]), "==", joindigits(indices[(j[g])[2]]))
            end
            print("\n")
        end
    end)
chevieset(:imp, :SemisimpleRank, function (p, q, r)
        return r
    end)
chevieset(:imp, :BraidRelations, function (p, q, r)
        local i, b, res
        b = function (i, j, o)
                local p
                p = function (i, j)
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
chevieset(:imp, :EigenvaluesGeneratingReflections, function (p, q, r)
        local res
        res = fill(0, max(0, (1 + r) - 1)) + 1 // 2
        if q == 1
            res[1] = 1 // p
        elseif q != p
            res = Concatenation([q // p], res)
        end
        return res
    end)
chevieset(:imp, :CartanMat, function (p, q, r)
        local rt, rbar, e
        rt = (chevieget(:imp, :GeneratingRoots))(p, q, r)
        rbar = ComplexConjugate(rt)
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
                    t = [[], [1], [1, 2]]
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
chevieset(:imp, :ClassInfo, function (p, q, r)
        local res, times, trans, I, i, j, a, S
        times = function (e, o)
                return Concatenation(map((x->begin
                                    o
                                end), 1:e))
            end
        if [q, r] == [2, 2] && !(haskey(CHEVIE, :othermethod))
            res = Dict{Symbol, Any}(:classtext => [], :classparams => [], :classnames => [])
            for i = 0:p - 1
                for j = 0:div((p - i) - 1, 2)
                    push!(res[:classparams], Concatenation(fill(0, max(0, (1 + j) - 1)) + 1, fill(0, max(0, (1 + i) - 1))))
                    push!(res[:classtext], Concatenation(fill(0, max(0, (1 + j) - 1)) + 1, times(i, [1, 2, 3])))
                    push!(res[:classnames], string(Concatenation(times(j, "1"), times(i, "z"))))
                end
            end
            for j = [2, 3]
                for i = 0:div(p, 2) - 1
                    push!(res[:classparams], Concatenation([j], fill(0, max(0, (1 + i) - 1))))
                    push!(res[:classtext], Concatenation([j], times(i, [1, 2, 3])))
                    push!(res[:classnames], string(Concatenation(string(j), times(i, "z"))))
                end
            end
            res[:malle] = []
            for a = 0:p - 1
                res[:malle] = Append(res[:malle], map((m->begin
                                    [3, a, m]
                                end), 0:div((p - a) - 1, 2)))
            end
            res[:malle] = Append(res[:malle], map((m->begin
                                [1, m]
                            end), 0:div(p, 2) - 1))
            res[:malle] = Append(res[:malle], map((m->begin
                                [2, m]
                            end), 0:div(p, 2) - 1))
            res[:orders] = map(function (c,)
                        if length(c) > 0 && c[1] in [2, 3]
                            return Lcm(2, div(p, gcd(count((x->begin
                                                        x == 0
                                                    end), c), p)))
                        else
                            return Lcm(div(p, gcd(count((x->begin
                                                        x == 0
                                                    end), c), p)), div(div(p, 2), gcd(count((x->begin
                                                        x == 1
                                                    end), c), div(p, 2))))
                        end
                    end, res[:classparams])
            res[:classes] = map(function (c,)
                        if length(c) > 0 && c[1] in [2, 3]
                            return div(p, q)
                        elseif 1 in c
                            return 2
                        else
                            return 1
                        end
                    end, res[:classparams])
            return res
        elseif q == 1
            res = Dict{Symbol, Any}(:classparams => partition_tuples(r, p))
            res[:classtext] = map(function (S,)
                        local l, w, d
                        S = Concatenation(map((i->begin
                                            map((t->begin
                                                        [t, i - 1]
                                                    end), S[i])
                                        end), 1:p))
                        SortBy(S, (a->begin
                                    [a[1], -(a[2])]
                                end))
                        l = 0
                        w = []
                        for d = S
                            w = Append(w, times(d[2], Concatenation(l + 1:l - (l + 1):2, 1:l + 1)))
                            w = Append(w, l + 2:l + d[1])
                            l = l + d[1]
                        end
                        return w
                    end, res[:classparams])
            res[:classnames] = map(chevieget(:imp, :ClassName), res[:classparams])
            res[:orders] = map((m->begin
                            Lcm(map(function (i,)
                                        if length(m[i]) == 0
                                            return 1
                                        else
                                            return Lcm((m[i] * p) // gcd(i - 1, p))
                                        end
                                    end, 1:length(m)))
                        end), res[:classparams])
            res[:centralizers] = map((m->begin
                            p ^ Sum(m, length) * Product(map((pp->begin
                                                Product(tally(pp), (y->begin
                                                            factorial(y[2]) * y[1] ^ y[2]
                                                        end))
                                            end), m))
                        end), res[:classparams])
            res[:classes] = map((x->begin
                            div(p ^ r * factorial(r), x)
                        end), res[:centralizers])
            return res
        else
            trans = function (w,)
                    local d, res, l, i, add, word
                    d = 0
                    res = []
                    word = function (l, i)
                            return map((j->begin
                                            1 + mod(j, 2)
                                        end), i + (l:(l - 1) - l:1))
                        end
                    add = function (a,)
                            local l
                            l = length(res)
                            if l > 0 && res[l] == a
                                res = res[1:l - 1]
                            elseif p == q && (a in [1, 2] && (l >= q && res[(l - q) + 1:l] == word(q, 3 - a)))
                                res = Concatenation(res[1:l - q], word(q - 1, 3 - a))
                            else
                                push!(res, a)
                            end
                        end
                    for l = w
                        if l == 1
                            d = d + 1
                        elseif l != 2
                            add(l)
                        else
                            d = mod(d, p)
                            if d == 0
                                add(2)
                            else
                                for i = 1:(p - d) - 1
                                    add(1)
                                    add(2)
                                end
                                add(1)
                            end
                        end
                    end
                    d = mod(d, p)
                    if mod(d, q) != 0
                        error()
                    elseif d != 0
                        res = Concatenation(1 + res, fill(0, max(0, (1 + div(d, q)) - 1)) + 1)
                    elseif p != q
                        res = 1 + res
                    end
                    return res
                end
            I = (chevieget(:imp, :ClassInfo))(p, 1, r)
            res = Dict{Symbol, Any}(:classtext => [], :classparams => [], :classnames => [], :orders => [], :centralizers => [])
            for i = filter((i->begin
                                mod(map(length, (I[:classparams])[i]) * (0:p - 1), q) == 0
                            end), 1:length(I[:classparams]))
                S = (I[:classparams])[i]
                a = Concatenation(S)
                push!(a, q)
                a = Append(a, filter((j->begin
                                        length(S[j]) != 0
                                    end), 1:p) - 1)
                a = ApplyFunc(gcd, a)
                for j = 0:a - 1
                    push!(res[:classtext], trans(Concatenation(fill(0, max(0, (1 + j) - 1)) + 1, (I[:classtext])[i], fill(0, max(0, (1 + (p - j)) - 1)) + 1)))
                    if a > 1
                        push!(res[:classparams], Concatenation(S, [(p * j) // a]))
                    else
                        push!(res[:classparams], S)
                    end
                    push!(res[:orders], (I[:orders])[i])
                    push!(res[:centralizers], ((I[:centralizers])[i] * a) // q)
                end
            end
            res[:classes] = map((x->begin
                            (res[:centralizers])[1] // x
                        end), res[:centralizers])
            res[:classnames] = map(chevieget(:imp, :ClassName), res[:classparams])
            return res
        end
    end)
chevieset(:imp, :ClassName, function (p,)
        local j, p1
        if IsList(p) && all(IsList, p)
            if Sum(p, Sum) == 1
                return FormatTeX(E(length(p), Position(p, [1]) - 1))
            else
                return PartitionTupleToString(p)
            end
        elseif IsList(p) && all(IsInt, p)
            return joindigits(p)
        elseif IsList(p) && (all(IsList, p[1:length(p) - 1]) && IsInt(p[length(p)]))
            p1 = p[1:length(p) - 1]
            p1 = Append(p1, [length(p1), p[length(p)]])
            return PartitionTupleToString(p1)
        else
            error()
        end
    end)
chevieset(:imp, :CharInfo, function (de, e, r)
        local d, ct, res, t, tt, s, fd
        res = Dict{Symbol, Any}()
        d = div(de, e)
        if e == 1
            res[:charparams] = partition_tuples(r, de)
            s = fill(0, max(0, (1 + d) - 1))
            s[1] = 1
            res[:charSymbols] = map((x->begin
                            symbol_partition_tuple(x, s)
                        end), res[:charparams])
        else
            res[:charparams] = []
            for t = partition_tuples(r, de)
                tt = map((i->begin
                                circshift(t, i)
                            end), (1:e) * d)
                if t == Minimum(tt)
                    s = Position(tt, t)
                    if s == e
                        push!(res[:charparams], t)
                    else
                        t = t[1:s * d]
                        s = div(e, s)
                        res[:charparams] = Append(res[:charparams], map((i->begin
                                            Concatenation(t, [s, i])
                                        end), 0:s - 1))
                    end
                end
            end
            if d == 1
                res[:charSymbols] = map((x->begin
                                symbol_partition_tuple(x, 0)
                            end), res[:charparams])
            end
            if d > 1 && (mod(e, 2) == 0 && r == 2)
                res[:malle] = map(function (t,)
                            local pos, de
                            if IsInt(t[length(t)])
                                if t[length(t)] == 0
                                    return [1, 2, 1, Position(t, [1])]
                                else
                                    return [1, 1, 2, Position(t, [1])]
                                end
                            else
                                de = div(length(t), 2)
                                pos = filter((i->begin
                                                length(t[i]) > 0
                                            end), 1:length(t))
                                if length(pos) == 1
                                    if t[pos[1]] == [2]
                                        return [1, 1, 1, pos[1] - de]
                                    else
                                        return [1, 2, 2, pos[1] - de]
                                    end
                                elseif pos[1] <= de
                                    return [2, -1, pos[1], pos[2] - de]
                                else
                                    return [2, 1, pos[2] - de, pos[1] - de]
                                end
                            end
                        end, res[:charparams])
            elseif [de, e, r] == [3, 3, 3]
                res[:malle] = [[2, 3, 2], [2, 3, 3], [2, 3, 1], [3, 4], [3, 5], [1, 9], [3, 2], [3, 1], [2, 3, 4], [1, 0]]
            elseif [de, e, r] == [3, 3, 4]
                res[:malle] = [[12, 6], [4, 10], [6, 8], [4, 11], [1, 18], [12, 3], [6, 5, 2], [8, 4], [8, 5], [6, 5, 1], [3, 9], [6, 2], [2, 6], [4, 2], [4, 1], [3, 3], [1, 0]]
            elseif [de, e, r] == [3, 3, 5]
                res[:malle] = [[30, 10], [20, 12], [5, 19], [10, 14], [10, 15], [5, 20], [1, 30], [30, 7, 1], [40, 6], [30, 7, 2], [10, 11], [15, 10], [20, 9], [20, 8], [15, 11], [10, 12], [4, 18], [30, 4], [20, 5], [10, 8], [10, 7], [20, 6], [5, 12], [20, 3], [10, 6], [15, 4], [15, 5], [10, 5], [6, 9], [10, 3], [10, 2], [5, 6], [5, 2], [5, 1], [4, 3], [1, 0]]
            elseif [de, e, r] == [4, 4, 3]
                res[:malle] = [[6, 3], [3, 6, 1], [3, 5], [3, 6, 2], [1, 12], [3, 2, 1], [3, 2, 2], [3, 1], [2, 4], [1, 0]]
            end
        end
        t = map(function (i,)
                    local v
                    v = map((x->begin
                                    []
                                end), 1:maximum([de, 2]))
                    if i > 0
                        v[1] = [i]
                    end
                    v[2] = fill(0, max(0, (1 + (r - i)) - 1)) + 1
                    return v
                end, r:(r - 1) - r:0)
        if e > 1
            t = map((v->begin
                            Minimum(map((i->begin
                                            circshift(v, i * d)
                                        end), 1:e))
                        end), t)
        end
        res[:extRefl] = map((v->begin
                        Position(res[:charparams], v)
                    end), t)
        if e == 1 || d == 1
            res[:A] = map(degree_gendeg_symbol, res[:charSymbols])
            res[:a] = map(valuation_gendeg_symbol, res[:charSymbols])
            res[:B] = map(degree_fegsymbol, res[:charSymbols])
            res[:b] = map(valuation_fegsymbol, res[:charSymbols])
        end
        if e > 1 && d > 1
            res[:opdam] = PermListList(res[:charparams], map(function (s,)
                            if !(IsList(s[length(s)]))
                                s = deepcopy(s)
                                t = div(length(s) - 2, d)
                                s[(0:t - 1) * d + 1] = circshift(s[(0:t - 1) * d + 1], 1)
                                s[1:length(s) - 2] = Minimum(map((i->begin
                                                    circshift(s[1:length(s) - 2], i * d)
                                                end), 1:t))
                                return s
                            end
                            s = copy(s)
                            s[(0:e - 1) * d + 1] = circshift(s[(0:e - 1) * d + 1], 1)
                            return Minimum(map((i->begin
                                                circshift(s, i * d)
                                            end), 1:e))
                        end, res[:charparams]))
        end
        return res
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
chevieset(:imp, :CharName, function (p, q, r, s, option)
        if ranksymbol(s) == 1
            return Format(E(length(s), Position(s, [1]) - 1), option)
        else
            return PartitionTupleToString(s, option)
        end
    end)
chevieset(:imp, :SchurModel, function (p, q, r, phi)
        local l, i, j, res, s, t, ci, GenHooks, v, h, d, p2
        if q == 1
            GenHooks = function (l, m)
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
        local m
        if r == 1
            return VcycSchurElement(Concatenation(para[1], [0]), ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi))
        elseif p == 1
            return VcycSchurElement([0, -((para[1])[1]) // (para[1])[2]], ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi))
        elseif q == 1
            return VcycSchurElement(Concatenation(para[1], [-((para[2])[1]) // (para[2])[2]]), ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi))
        elseif [q, r] == [2, 2]
            return VcycSchurElement(Concatenation(para[[2, 3, 1]]), ((CHEVIE[:imp])[:SchurModel])(p, q, r, phi), ((CHEVIE[:imp])[:SchurData])(p, q, r, phi))
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
            ((CHEVIE[:compat])[:InfoChevie])("# FactorizedSchurElements(H(G(", p, ",", q, ",", r, "),", para, ") not implemented\n")
            return false
        end
    end)
chevieset(:imp, :HeckeRepresentation, function (p, q, r, para, root, i)
        local X, Y, t, x, a, v, d, T, S, m, extra, l, m1, p1rRep, f, e
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
            f = function (x, j)
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
            m = function (i,)
                    local f1, f2, f3, f5, f7, f8, f11, f13
                    f1 = (x->begin
                                x ^ 0 * [[[x, -1, 0, 0, 0, 0, 0, 0, 0, 0, ((1 - x) - x ^ 2) + x ^ 3, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1 + x, 0, x, 0, -x, 0, 0, 0, x - x ^ 2, 0], [0, 0, 0, -1 + x, 0, 0, -x, 0, 0, 0, x - x ^ 2, 0], [0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 1, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, 0, x, 0, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1 + x]], [[0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0], [0, 0, 0, 0, x, 0, 0, 1 - x, 0, 0, 0, 0], [0, 0, 0, 1, -1 + x, 0, 0, 1 - x, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1 - x, -1 + x, 1, 0, (1 - x) + x ^ 2], [0, 0, 1, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x], [0, 0, 0, 0, 0, x, 0, x - x ^ 2, -x, -1 + x, 0, x - x ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 + x]], [[0, (-1 + 2x) - x ^ 2, 1 - x, x, -x + x ^ 2, 0, 0, 0, (-1 + 2x) - x ^ 2, 0, 0, 0], [0, -1 + x, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0], [0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0], [1, -1 + x, 0, -1 + x, 0, 0, 1 - x, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1 + x, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, -1], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, -x, 0, 0, 0, 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0], [0, 0, 0, 0, 0, x, 0, 0, 0, 0, 1 - x, x - x ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, x], [0, 0, 0, 1, 0, -1 + x, -1 + x, 0, 0, 0, 1 - x, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x], [0, 0, 0, 0, 0, 0, 0, x, 0, 0, -1, 0], [0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0], [0, 0, 0, 0, x, 0, -x, 0, 0, -1 + x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + x]]]
                            end)
                    f2 = function (x, j)
                            return [[[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, x, x, x]], [[-1, 0, 0, 0], [0, -1, 0, 0], [0, -(j ^ 2), x, 1], [0, 0, 0, -1]], [[-1, 0, 0, 0], [x, x, -j * x, 1], [0, 0, -1, 0], [0, 0, 0, -1]], [[x, 1, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]]
                        end
                    f3 = (x->begin
                                [[[x, -1, 0, -x, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1 + x, 0, 1, x], [0, 0, 0, -1, 0, 0], [0, 0, x, 0, 0, x], [0, 0, 0, 0, 0, -1]], [[-1, 0, 0, 0, 0, 0], [0, 0, 0, 0, x, 1], [0, 0, -1, 0, 0, 0], [-1, 0, -1, x, 0, -1 + x], [0, 1, 0, 0, -1 + x, 1], [0, 0, 0, 0, 0, -1]], [[0, x, 1, -1, -1, 0], [1, -1 + x, 1, -1, -1, 0], [0, 0, -1, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, 1, 1, x]], [[x, -1, 0, 0, 1, x], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 0, -1, x, 1, x], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]]]
                            end)
                    f5 = (x->begin
                                [[[-1]], [[-1]], [[-1]], [[-1]]]
                            end)
                    f7 = function (x, j)
                            return [[[-1, 0, 0, 0, 0, 0], [x, x, 0, 0, 0, 0], [x, 0, x, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, -j * x ^ 2, x, x]], [[x, 1, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, j ^ 2, x, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, x, x, 1], [0, 0, 0, 0, 0, -1]], [[x, 0, 1, 0, 1, 0], [0, x, j * x, 0, 0, 1], [0, 0, -1, 0, 0, 0], [0, 0, 0, x, 1, -(j ^ 2) * x ^ -1], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]], [[-1, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 0, 0, x, 0, 0], [x, 0, 0, 0, x, 0], [0, x, 0, 0, 0, x]]]
                        end
                    f8 = function (x, j)
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
            return -((para[2])[2]) * m(i)
        elseif [p, q, r] == [3, 3, 5]
            x = -((para[2])[1]) // (para[2])[2]
            m = function (i,)
                    local r, f1, f2, f3, f4, f8, f9, f11, f12, f13, f17, f20, f23, f29
                    f1 = function (x,)
                            return x ^ 0 * [[[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [E(3, 2) - E(3, 2) * x, E(3, 2) - E(3, 2) * x, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -(E(3, 2)) * x + E(3, 2) * x ^ 2, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(-(E(3, 2)) - ER(-3) * x) + E(3) * x ^ 2, (E(3) - E(3) * x) + E(3) * x ^ 2, 0, 0, 0, E(3) - E(3) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E(3, 2) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, E(3), -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, (E(3, 2) * x + ER(-3) * x ^ 2) - E(3) * x ^ 3, 0, (E(3) - E(3) * x) + E(3) * x ^ 2, 0, 0, 0, E(3) - E(3) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * x, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * x, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, (-(E(3, 2)) - ER(-3) * x) + E(3) * x ^ 2, 0, 0, 0, 0, 0, 0, 0, (E(3) - E(3) * x) + E(3) * x ^ 2, 0, 0, 0, 0, E(3) - E(3) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, ((E(3, 2) - E(3, 2) * x) - E(3, 2) * x ^ 2) + E(3, 2) * x ^ 3, 0, (E(3) - ER(-3) * x) - E(3, 2) * x ^ 2, 0, 0, 0, -(E(3, 2)) + E(3, 2) * x, 0, -(E(3)) + ER(-3) * x + E(3, 2) * x ^ 2, 0, 0, 0, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0], [E(3) - E(3) * x, 0, 0, (2 * E(3) - E(3) * x ^ -1) - E(3) * x, 0, 0, 0, 0, 0, 0, (-(E(3)) + 2 * E(3) * x) - E(3) * x ^ 2, 0, 0, 0, 0, -(E(3)) + E(3) * x, 0, 0, 0, E(3) - E(3) * x, 0, 0, E(3) - E(3) * x, 0, E(3), 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * x, 0, -1 + x, 0, 0, 0, 0, 0, 0], [-(E(3)) + ER(-3) * x + E(3, 2) * x ^ 2, 0, (x - 2 * x ^ 2) + x ^ 3, (1 - 3 * ER(-3)) // 2 + E(3) * x ^ -1 + (1 + 3 * ER(-3)) // 2 * x + E(3, 2) * x ^ 2, 0, 0, 0, 0, 0, 0, ((E(3) - 2 * E(3) * x) + 2 * E(3) * x ^ 2) - E(3) * x ^ 3, 0, 0, 0, 0, (E(3) - 2 * E(3) * x) + E(3) * x ^ 2, 0, -(E(3)) + E(3) * x, 0, (-(E(3)) + 2 * E(3) * x) - E(3) * x ^ 2, 0, 0, (E(3, 2) - E(3, 2) * x) + E(3, 2) * x ^ 2, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, -(E(3, 2)) + E(3, 2) * x, 0, 0, 0, -1 + x ^ -1, 0, E(3, 2) - E(3, 2) * x, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), -1 + x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, E(3) * x], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, x, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 0, 0, 0, 0, 0, -1 + x, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + x ^ -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, x, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2x) + x ^ 2, (1 - 2x) + x ^ 2, (1 - 2x) + x ^ 2, 1 - x, (2 - x ^ -1) - x, 1 - x, 0, 0, 1 - x ^ -1, -1 + x, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2x) + x ^ 2, -x + x ^ 2, (1 - 2x) + x ^ 2, 1 - x, 1 - x, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(2 - x ^ -1) - x, 1 - x, 0, 0, 0, 1 - x, -1 + x, -1 + x, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, x - x ^ 2, (-1 + 2x) - x ^ 2, 0, 0, 0, x - x ^ 2, 0, 0, 0, 0, 0, -1 + x, 0, 0, -1 + x, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(2 - x ^ -1) - x, 0, 1 - x, 0, 0, 0, -1 + x, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - x ^ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [((-1 + 3x) - 3 * x ^ 2) + x ^ 3, ((-1 + 3x) - 3 * x ^ 2) + x ^ 3, (x - 2 * x ^ 2) + x ^ 3, (2 - x ^ -1) - x, (-2 + x ^ -1 + 2x) - x ^ 2, (-1 + 2x) - x ^ 2, 0, (1 - 2x) + x ^ 2, -2 + x ^ -1 + x, 0, 0, -1 + x, -1 + x, (2 - x ^ -1) - x, 1 - x, 0, 0, 0, 1 - x, 0, 0, -1 + x, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [(-3 + x ^ -1 + 3x) - x ^ 2, (-1 + 2x) - x ^ 2, 0, (1 - 2x) + x ^ 2, 0, (-1 + 2x) - x ^ 2, (1 - 2x) + x ^ 2, (1 - 2x) + x ^ 2, -1 + x, x - x ^ 2, 0, (-1 + 2x) - x ^ 2, 0, 1 - x, 1 - x, 0, -1 + x, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, -1 + x, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [((1 - 3x) + 3 * x ^ 2) - x ^ 3, ((1 - 2x) + 2 * x ^ 2) - x ^ 3, ((1 - 3x) + 3 * x ^ 2) - x ^ 3, (1 - 2x) + x ^ 2, ((3 - x ^ -1) - 3x) + x ^ 2, (1 - 2x) + x ^ 2, 0, 0, (2 - x ^ -1) - x, 0, (-1 + 2x) - x ^ 2, 0, 1 - x, 0, 0, 0, 0, -1 + x, 0, 0, 1 - x, 0, (1 - 2x) + x ^ 2, 0, 1 - x, 0, -1 + x, x, 0, 0], [0, 0, 0, (2 - x ^ -1) - x, 0, -1 + x, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, -2 + x ^ -1 + x, 0, (2 - x ^ -1) - x, 0, 1 - x, -1 + x ^ -1, 0, -1 + x, 0, 0, 0, 1, 0, 0, 0], [(-1 + 2x) - x ^ 2, 0, 0, 0, 0, 0, (1 - 2x) + x ^ 2, -x + x ^ 2, 0, 0, 0, 0, 0, 0, 0, 1 - x, -1 + x, 0, 0, 0, 0, 0, 1 - x, 0, 0, 1 - x, 0, 0, 0, x], [0, 0, (1 - 2x) + x ^ 2, ((3 - x ^ -1) - 3x) + x ^ 2, 0, 0, 1 - x, (1 - 2x) + x ^ 2, 0, 0, -1 + x, 0, 0, (2 - x ^ -1) - x, 0, 0, 0, 0, 1 - x, 0, 0, 0, (2 - x ^ -1) - x, 0, 1 - x ^ -1, 1 - x, 0, 0, 1, -1 + x]], [[0, -x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, x, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0], [0, x - x ^ 2, 0, 0, 0, 0, 0, 0, 0, (1 - 2x) + x ^ 2, (-1 + 2x) - x ^ 2, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 1 - x, 0, (1 - 2x) + x ^ 2, 0, 0, 0, -1 + x, x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 0, x], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0], [1 - x, 0, 0, (2 - x ^ -1) - x, 0, 0, 0, 0, 0, 0, (-1 + 2x) - x ^ 2, 0, 0, 0, 0, -1 + x, 0, 0, 0, 1 - x, 0, 0, 1 - x, 0, 1, 0, 0, -1 + x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 1 - x ^ -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + x]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, x, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 0, 0, 0, 0, 0, -1 + x, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 1 - x, 0, 0, 0, 1, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, x - x ^ 2, 0, 0, 0, 0, 0, 0, 0, -x + x ^ 2, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -x + x ^ 2, 0, 1 - x, 0, 0, 0, 1, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 1, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, -1 + x, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 + x]], [[0, 0, -x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, x, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, x, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 0, 0, 0, 0, 0, -1 + x, x, 0, 0, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - x, 0, 0, (2 - x ^ -1) - x, 0, 0, 0, 0, 0, 0, (-1 + 2x) - x ^ 2, 0, 0, 0, 0, -1 + x, 0, 0, 0, 1 - x, 0, 0, 1 - x, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -x + x ^ 2, 0, 0, 0, 0, 1 - x, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0], [0, 0, x - x ^ 2, (-1 + 2x) - x ^ 2, 0, 0, 0, x - x ^ 2, 0, 0, 0, 0, 0, -1 + x, 0, 0, -1 + x, 0, x, 0, 0, 0, 0, 0, -1 + x, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x, 0, 0, 0, 0, -1 + x, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + x, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]]]
                        end
                    f2 = function (q,)
                            return q ^ 0 * [[[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-(E(3)) + ER(-3) * q + E(3, 2) * q ^ 2, E(3, 2) - E(3, 2) * q, 0, 0, -(E(3)) + E(3) * q, 0, (ER(-3) - E(3) * q ^ -1) + E(3, 2) * q, -(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [E(3) * q - E(3) * q ^ 2, E(3) * q, 0, 0, 0, 0, E(3) - E(3) * q, E(3) - E(3) * q, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, E(3, 2) - E(3, 2) * q, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, E(3, 2), 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [-(E(3)) + E(3) * q, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, -1 + q, -(E(3)), 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, (-(E(3, 2)) - ER(-3) * q) + E(3) * q ^ 2, 0, (-(E(3, 2)) - ER(-3) * q) + E(3) * q ^ 2, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, -1 + q, 0, 0], [(E(3) * q - 2 * E(3) * q ^ 2) + E(3) * q ^ 3, E(3) * q - E(3) * q ^ 2, 0, 0, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, E(3, 2) * q], [0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, E(3), -1 + q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q ^ 3 - q ^ 4, 1 - q, -1 + q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, 0, q ^ 3 - q ^ 4, 0, 0, 0, 0, 0, 1 - q, -1 + q, 0, 0, 1, 1 - q, 0, 0, 0, 0, 0, 0], [0, 0, q ^ 3 - q ^ 4, 1 - q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, 0, q ^ 3 - q ^ 4, 0, 0, 0, 0, 0, 1 - q, q, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, -1 + q, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, 0, 0, 0, 1 - q ^ -1, -1 + q, 1, 0, 0], [0, 0, 0, 1 - q, 0, -1 + q, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, q, 0, 0, 0], [(-1 + 2q) - q ^ 2, 0, 0, 0, q - q ^ 2, (1 - 2q) + q ^ 2, 0, 0, 0, 0, 0, -1 + q, 0, 0, -1 + q, 0, 0, 0, -1 + q, q], [(2 - q ^ -1) - q, 0, (-(q ^ 3) + 2 * q ^ 4) - q ^ 5, (-1 + 2q) - q ^ 2, (1 - 2q) + q ^ 2, -2 + q ^ -1 + q, 0, 0, 0, 0, -1 + q, 0, 0, -1 + q, (2 - q ^ -1) - q, -1 + q ^ -1, 0, 0, 1, 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q]], [[0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -(q ^ -2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -(q ^ 3), -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + q, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0, -1 + q, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0], [0, 0, 0, 0, 0, 0, 1 - q, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1 - q, 0, 1 - q, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, q, 1 - q, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[-1 + q, 0, -(q ^ 3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, 0, q ^ 3 - q ^ 4, 0, 0, 0, 0, 0, 1 - q, -1 + q, 0, 0, 1, 1 - q, 0, 0, 0, 0, 0, 0], [-(q ^ -2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [q - q ^ 2, q, 0, 0, 0, 0, 1 - q, 1 - q, 0, 0, 1 - q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, q ^ 3 - q ^ 4, 1 - q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]]]
                        end
                    f3 = function (q, j)
                            return q ^ 0 * [[[-1, 0, 0, 0, 0], [0, -1, 0, 0, 0], [1, 0, 0, 0, -1], [1 + j * q, 0, 1 + j * q, -1, -1 - j * q], [-q, 0, -q, 0, -1 + q]], [[-1, 0, 0, 0, 0], [0, -1, 0, 0, 0], [1, 0, 0, -1, 0], [-q, 0, -q, -1 + q, 0], [-j, 0, -j, j, -1]], [[-1, 0, -1, 0, 0], [0, -1, 1, 0, 0], [0, 0, q, 0, 0], [0, 0, 0, -1, 0], [0, 0, 0, 0, -1]], [[q, 0, 0, 0, 0], [-1, -1, 0, 0, 0], [-q, 0, -1, 0, 0], [1, 0, 0, -1, 0], [1, 0, 0, 0, -1]], [[0, 1, 0, 0, 0], [q, -1 + q, 0, 0, 0], [0, 0, -1, 0, 0], [1, 1, 0, -1, 0], [1, 1, 0, 0, -1]]]
                        end
                    f4 = function (q, j)
                            return q ^ 0 * [[[-1 + q, 0, 0, q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, (j ^ 2 * q + (-(j ^ 2) + j) * q ^ 2) - j * q ^ 3, 0, j - j * q, 0, (-j + j * q) - j * q ^ 2, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, j ^ 2 * q - j ^ 2 * q ^ 2, 0, -(j ^ 2), 0, j ^ 2 - j ^ 2 * q, 0, 0, 0], [-q + q ^ 2, 0, (j ^ 2 * q - 2 * j ^ 2 * q ^ 2) + j ^ 2 * q ^ 3, -q + q ^ 2, -(j ^ 2) + j ^ 2 * q, 0, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, -1, 0, 0], [0, q - q ^ 2, (j ^ 2 * q - 2 * j ^ 2 * q ^ 2) + j ^ 2 * q ^ 3, 0, -(j ^ 2) + j ^ 2 * q, 0, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 0, -1, q - q ^ 2], [0, q, 0, 0, 0, 0, 0, 0, 0, -1 + q]], [[0, 0, j ^ 2 * q - j ^ 2 * q ^ 2, j ^ 2 * q, 0, 0, 0, 0, 0, 0], [0, -1 + q, -q + q ^ 2, 0, 0, 0, 0, 0, 0, j], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [j, 0, q - q ^ 2, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1 + q, 0, -q, 0, 0, 0], [q - q ^ 2, -(j ^ 2) * q + j ^ 2 * q ^ 2, ((j ^ 2 * q - j ^ 2 * q ^ 2) - j ^ 2 * q ^ 3) + j ^ 2 * q ^ 4, j ^ 2 * q ^ 2 - j ^ 2 * q ^ 3, 0, -1, 0, 0, 0, -1 + q], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, j ^ 2 * q, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, 0, 0, 0, 0, 0, 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, q, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, -1 + q, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, -q, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, -1 + q, 0, q, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, -q, -1 + q, 0], [q, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, -q + q ^ 2, -(q ^ 2) + q ^ 3, 0, -1 + q, 0, 0, 0, 0, j * q], [0, 0, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, 0, j ^ 2 - j ^ 2 * q, -(j ^ 2), 0, 0], [0, -q, 0, 0, 0, 0, -1 + q, 0, 0, 0], [0, -q + q ^ 2, 0, q ^ 2 - q ^ 3, 0, -j * q, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, j ^ 2, 0, -(j ^ 2) + j ^ 2 * q, 0, 0, 0]]]
                        end
                    f11 = function (q, j)
                            return q ^ 0 * [[[q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-(j ^ 2) + j ^ 2 * q, j ^ 2 - j ^ 2 * q, j ^ 2, 0, 0, 0, 0, 0, 0, 0], [-j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, (j - j * q) + j * q ^ 2, j - j * q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, q, -1 + q, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, -1 + q, 0, j, 0, 0, 0, 0], [1 - q, 0, 0, 0, -1 + q, 0, j, 0, 0, 0], [-(j ^ 2) * q + j ^ 2 * q ^ 2, 0, 0, j ^ 2 * q, 0, 0, 0, 0, 0, 0], [-(j ^ 2) * q + j ^ 2 * q ^ 2, 0, 0, 0, j ^ 2 * q, 0, 0, 0, 0, 0], [(-(j ^ 2) + 2 * j ^ 2 * q) - j ^ 2 * q ^ 2, j ^ 2 - j ^ 2 * q, j ^ 2 - j ^ 2 * q, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, -1 + q, 0, -1, 0, 0], [0, 0, 0, j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) * q + j ^ 2 * q ^ 2, 1 - q, -1 + q, 0, -1, 0], [(-(j ^ 2) + 2 * j ^ 2 * q) - j ^ 2 * q ^ 2, j ^ 2 - j ^ 2 * q, j ^ 2 - j ^ 2 * q, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, -1 + q, 0, 0, -1]], [[0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [-q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, q, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, q, 0, 0, -1 + q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [-1 + q, 0, 0, 1 - q, 0, -j, 0, 0, 0, 0], [0, -q, 0, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) * q + j ^ 2 * q ^ 2, -(j ^ 2) * q, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1 + q, -q], [0, 0, 0, 0, 0, 0, 0, 0, -1, 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, -q, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, -q, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1], [0, 0, 0, 0, 0, 0, 0, 0, q, 0], [0, 0, 0, 0, 0, 0, 0, -q, 0, -1 + q]]]
                        end
                    f12 = function (q, j)
                            return q ^ 0 * [[[0, 0, 0, 0, -(j ^ 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -(j ^ 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -(j ^ 2), 0, 0, 0, 0, 0, 0, 0, 0], [-j * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -j, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -j * q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [-(j ^ 2) + j ^ 2 * q, 0, 0, 0, j - j * q, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, j ^ 2 - j ^ 2 * q, 0, 0, -j + j * q, 0, -1, 0, 0, 0, 0, 0, 0], [0, (2 * j ^ 2 - j ^ 2 * q ^ -1) - j ^ 2 * q, 0, 0, 0, (j - 2 * j * q) + j * q ^ 2, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, (-1 + 2q) - q ^ 2, -q, 0, 0, (1 - 2q) + q ^ 2, 0, 0, 0, q, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [j - j * q, 0, 0, -q + q ^ 2, (j ^ 2 + 2j + q ^ -1) - j * q, 0, 1 - q, -q, -q, 0, 0, 0, q, 0, 0], [0, -2 * j ^ 2 + j ^ 2 * q ^ -1 + j ^ 2 * q, 0, 0, 0, (-j + 2 * j * q) - j * q ^ 2, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, (1 - 2q) + q ^ 2, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, 0, q, 0, q, q]], [[-1 + q, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + q, 0, -1, 0, -1 + q ^ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1 + q, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [-q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -q + q ^ 2, 0, 0, 1 - q, 0, -1, 0, 0, 0, 0, 0, 0], [0, -q + q ^ 2, 0, 0, 0, q - q ^ 2, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, (-1 + 2q) - q ^ 2, -q, 1 - q, 1 - q, -q + q ^ 2, 1 - q, 0, 0, q, q, 0, 0, 0, 0], [0, 0, 0, -1 + q, 0, 0, -1 + q ^ -1, 0, 0, 0, 0, -1, 0, 0, 0], [-(j ^ 2) + j ^ 2 * q, 1 - q ^ -1, 0, (1 - 2q) + q ^ 2, -(j ^ 2) + j ^ 2 * q, -1 + q, 1 - q, -q, -q, 0, 0, 0, q, 0, 0], [0, q - q ^ 2, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, -1, 0], [1 - q ^ -1, (((2 + q ^ -2) - 2 * q ^ -1) - 2q) + q ^ 2, 0, (2 - q ^ -1) - q, 1 - q ^ -1, (-2 + q ^ -1 + 2q) - q ^ 2, 1 - q ^ -1, 0, 0, 0, 0, q, 0, q, q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q, -1, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, (1 - 2q) + q ^ 2, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, -q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, -j + q ^ -2 + (j ^ 2 + 2j) * q ^ -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [-j + (j ^ 2 + 2j) * q + q ^ 2, 1 - q, 0, 0, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, -1 + q, q, q, 0, 1 - q, 0, -1, 0, 0], [0, (((3 - q ^ -1) - 4q) + 3 * q ^ 2) - q ^ 3, 0, 0, 0, ((-2 + 3q) - 3 * q ^ 2) + q ^ 3, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, ((-4 + q ^ -1 + 6q) - 4 * q ^ 2) + q ^ 3, 0, 0, 0, ((2 - 5q) + 4 * q ^ 2) - q ^ 3, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1 + q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-j + j * q, 0, 0, 0, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, q ^ -1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [-j + j * q, 0, 1, 0, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, 0, 0, 1 - q, 0, 0, -1 + q, -1, 0, 0, 0, 0, 0], [0, 0, 0, -(q ^ 2) + q ^ 3, 0, 0, q - q ^ 2, 0, -q, 0, 0, 0, 0, 0, 0], [-j + j * q, 0, 0, q - q ^ 2, ((-(j ^ 2) - 2j) - q ^ -1) + j * q, 0, -1 + q, q, q, 0, 0, 0, -1, 0, 0], [0, (2 - q ^ -1) - q, 0, 0, 0, -2 + q ^ -1 + q, 0, 0, 0, 0, 0, -1 + q, 0, -1, 0], [0, (1 - 2q) + q ^ 2, q, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, -q, -q, 0, -1 + q, 0, 0], [0, 0, 0, (-q + 2 * q ^ 2) - q ^ 3, 0, 0, (1 - 2q) + q ^ 2, 0, 0, 0, 0, -q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, (-1 + 2q) - q ^ 2, 0, 0, 0, (1 - 2q) + q ^ 2, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, ((2 * j ^ 2 + j) - j ^ 2 * q ^ -1) + q, 0, 0, 0, -(j ^ 2) + j ^ 2 * q, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q ^ -1, 0, 0, 0, 0, 1, 0, 0, 0], [(1 - 2q) + q ^ 2, 0, 1, 0, (2 - q ^ -1) - q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, q - q ^ 2, 0, 0, -1 + q, 0, q, 0, 0, -1 + q, 0, 0, 0], [0, (1 - 2q) + q ^ 2, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, 0, q, 0, q, q], [(j ^ 2 + (-2 * j ^ 2 - j) * q) - q ^ 2, 0, 0, 0, j ^ 2 - j ^ 2 * q, 0, 0, q, 0, 0, 0, 0, 0, -1 + q, 0], [j - j * q, 0, 0, -q + q ^ 2, (j ^ 2 + 2j + q ^ -1) - j * q, 0, 1 - q, -q, -q, 0, 0, 0, 1, 0, -1 + q]]]
                        end
                    f13 = function (q, j)
                            return q ^ 0 * [[[j - j * q, (-1 + q) - q ^ 2, (-(j ^ 2) + 2 * j ^ 2 * q) - j ^ 2 * q ^ 2, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, -(j ^ 2) + (j ^ 2 - j) * q + j * q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, j ^ 2 - j ^ 2 * q, -1 + q, 1 - q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, 0, -q + q ^ 2, 0, -q + q ^ 2, 1 - q, 0, j ^ 2 - j ^ 2 * q, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 1 - q, -1 + q, 0, 0, -1 + q, 0, j ^ 2 - j ^ 2 * q, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, (-(j ^ 2) + 2 * j ^ 2 * q) - j ^ 2 * q ^ 2, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 0, 0, (j ^ 2 + (-(j ^ 2) + j) * q) - j * q ^ 2, 0, (1 - q) + q ^ 2, 0, j - j * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j ^ 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1 - q, -1 + q, 1 - q, 0, 0, 0, 0, 0, j ^ 2 - j ^ 2 * q, 0, j ^ 2 * q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, (2 - q ^ -1) - q, ((2 * j ^ 2 + j) - j ^ 2 * q ^ -1) + q, ((-2 * j ^ 2 - j) - q ^ -1) + j ^ 2 * q, 0, 0, 0, 0, 0, -j + j * q ^ -1 + j * q, 0, j - j * q, 0, 0, 0, 0, 0], [0, (-1 + q) - q ^ 2, -q + q ^ 2, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 1 - q, ((j ^ 2 - j) + j * q ^ -1) - j ^ 2 * q, 0, -1 + q ^ -1 + q, 0, 0, 0, 0, 0, 0, 0, j - j * q, 0, 0, 0, 0], [0, 0, 0, 0, ((j ^ 2 - j) - j ^ 2 * q ^ -1) + j * q, 0, -2j + j * q ^ -1 + j * q, 0, 0, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, 0, 0, 0, 0, 0, 0, j - j * q, (1 - q ^ -1) - q, 0, 0], [0, 0, 0, 0, -1 + q, 0, -1 + q, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, -q, j ^ 2 - j ^ 2 * q, 0, 0], [0, (j * q + (j ^ 2 - j) * q ^ 2) - j ^ 2 * q ^ 3, (-q + 2 * q ^ 2) - q ^ 3, 0, (1 - 2q) + q ^ 2, 0, ((1 - q) - q ^ 2) + q ^ 3, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, (j * q + (j ^ 2 - j) * q ^ 2) - j ^ 2 * q ^ 3, 0, -q + q ^ 2, 0, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 0, j ^ 2 * q - j ^ 2 * q ^ 2, q - q ^ 2, 0, 0, -1, 0], [1 - q, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, 0, ((-1 + q ^ -1) - q) + q ^ 2, 0, ((-1 + q ^ -1) - q) + q ^ 2, 0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, (2 - q ^ -1) - q, -1 + q, -j + j * q ^ -1, 0, -1 + q ^ -1, 0, 0, 1 - q, ((j ^ 2 - j) + j * q ^ -1) - j ^ 2 * q, 0, -1]], [[-1 + q, -(j ^ 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-j * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, j ^ 2 * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, 1 - q, 0, 0, -1 + q, 0, -1 + q, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, 0, 0, 0, q - q ^ 2, 0, q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0], [0, -(j ^ 2) * q, 0, 0, 0, 0, 0, j ^ 2, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, -(j ^ 2), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -j * q, 0, 0, 0], [0, j - j * q, 0, 0, 0, 0, 0, j - j * q ^ -1, j - j * q, 0, j ^ 2 - j ^ 2 * q, 0, j - j * q ^ -1, 0, -j + j * q, -(j ^ 2) + j ^ 2 * q, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, 0, j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, j ^ 2 * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, j ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, j ^ 2 * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, j * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, j, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q + q ^ 2, q - q ^ 2, 0, 0, q - q ^ 2, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, 0, -q, 0, 0, q, 0, 0, 0, 0, 0, 0], [0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, (-1 + 2q) - q ^ 2, 0, -2 + q ^ -1 + q, 0, ((-1 + q ^ -1) - q) + q ^ 2, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, -1 + q, 0, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, 0, j ^ 2 - j ^ 2 * q, 1 - q, 0, 0, -1, 0], [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [-1 + q, ((-j + 3 * j * q) - 3 * j * q ^ 2) + j * q ^ 3, ((j - 3 * j * q) + 3 * j * q ^ 2) - j * q ^ 3, (-q + 2 * q ^ 2) - q ^ 3, ((j ^ 2 + 3j) - j * q ^ -1) + (-2 * j ^ 2 - 4j) * q + (j ^ 2 + 2j) * q ^ 2, -(j ^ 2) + (2 * j ^ 2 + j) * q + q ^ 2, ((j ^ 2 + 2j) - j * q ^ -1) + (-2 * j ^ 2 - j) * q + (j ^ 2 - j) * q ^ 2 + j * q ^ 3, (((-2 * j ^ 2 - 3j) - q ^ -1) + (j ^ 2 + 3j) * q) - j * q ^ 2, (-j + (j ^ 2 + 3j) * q + (-2 * j ^ 2 - 3j) * q ^ 2) - q ^ 3, 0, -(j ^ 2) + (j ^ 2 - j) * q + j * q ^ 2, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, (-1 + q) - q ^ 2, -j + (-(j ^ 2) + j) * q + j ^ 2 * q ^ 2, 0, 0, j - j * q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, 0, -1 + q]], [[q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j ^ 2 * q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(j ^ 2), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1 - q, 0, 1 - q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, q, -(j ^ 2) + j ^ 2 * q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, 1 - q, 0, 0, -1 + q, 0, -1 + q, 0, 1, 0, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -j * q, 0, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0], [0, -j * q + (-(j ^ 2) + j) * q ^ 2 + j ^ 2 * q ^ 3, (q - 2 * q ^ 2) + q ^ 3, 0, (-1 + 2q) - q ^ 2, 0, (-1 + q + q ^ 2) - q ^ 3, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, -j * q + (-(j ^ 2) + j) * q ^ 2 + j ^ 2 * q ^ 3, 0, q - q ^ 2, 0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, -(j ^ 2) * q + j ^ 2 * q ^ 2, -q + q ^ 2, 0, 0, q, q]], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, 0, 0, 0, 0, 0, 0, 0], [0, q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, j ^ 2, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0], [-q, 0, -q + q ^ 2, 0, -q + q ^ 2, 1 - q, 0, j ^ 2 - j ^ 2 * q, 0, 0, 0, 0, 0, -1 + q, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j * q, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, 1 - q, 0, 0, -1 + q, -1, -1 + q, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, j ^ 2, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1 + q, 0, 0], [q - q ^ 2, (j * q + (j ^ 2 - j) * q ^ 2) - j ^ 2 * q ^ 3, 0, 0, ((1 - q) - q ^ 2) + q ^ 3, 0, ((1 - q) - q ^ 2) + q ^ 3, 0, (j * q + (j ^ 2 - j) * q ^ 2) - j ^ 2 * q ^ 3, (-1 + 2q) - q ^ 2, -q + q ^ 2, j - j * q, 0, 1 - q, 0, 0, q - q ^ 2, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, -q], [0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, (-1 + 2q) - q ^ 2, 0, -2 + q ^ -1 + q, 0, ((-1 + q ^ -1) - q) + q ^ 2, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, -1 + q, 0, ((-(j ^ 2) + j) - j * q ^ -1) + j ^ 2 * q, 0, j ^ 2 - j ^ 2 * q, 1 - q, 0, 0, -1, -1 + q]]]
                        end
                    f17 = function (q,)
                            return q ^ 0 * [[[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [-(q ^ 4), q ^ 3, -(q ^ 2), q]], [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [-(q ^ 4), q ^ 3, -(q ^ 2), q]], [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 0, 1], [0, 0, q, -1 + q]], [[-1, 0, 0, 0], [0, 0, 1, 0], [0, q, -1 + q, 0], [0, 0, 0, -1]], [[0, 1, 0, 0], [q, -1 + q, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]]
                        end
                    f20 = function (q, j)
                            return q ^ 0 * [[[q, 0, -(q ^ 2), 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -(q ^ 2), q, 0, 0, 0, 0, 0], [0, -(q ^ 2), 0, 0, 0, q, 0, 0, 0, 0], [-j * q + (-(j ^ 2) + j) * q ^ 2 + j ^ 2 * q ^ 3, (j ^ 2 * q + (-(j ^ 2) + j) * q ^ 2) - j * q ^ 3, (-(j ^ 2) * q ^ 2 + 2 * j ^ 2 * q ^ 3) - j ^ 2 * q ^ 4, (j ^ 2 * q - 2 * j ^ 2 * q ^ 2) + j ^ 2 * q ^ 3, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, j - j * q, (-j + j * q) - j * q ^ 2, 0, 0], [j ^ 2 * q - j ^ 2 * q ^ 2, j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) * q ^ 2 + j ^ 2 * q ^ 3, j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) + j ^ 2 * q, 0, -(j ^ 2), j ^ 2 - j ^ 2 * q, 0, 0], [-j * q ^ 2 + (-(j ^ 2) + j) * q ^ 3 + j ^ 2 * q ^ 4, (-j * q ^ 2 + 2 * j * q ^ 3) - j * q ^ 4, (j ^ 2 * q ^ 2 + (-2 * j ^ 2 + j) * q ^ 3 + (2 * j ^ 2 - j) * q ^ 4) - j ^ 2 * q ^ 5, (j ^ 2 * q ^ 2 + (-(j ^ 2) + j) * q ^ 3) - j * q ^ 4, 0, (j + (j ^ 2 - j) * q) - j ^ 2 * q ^ 2, 0, (-j * q + j * q ^ 2) - j * q ^ 3, j - j * q, ((j - j * q) + j * q ^ 2) // q], [0, 0, -(j ^ 2) * q ^ 3 + j ^ 2 * q ^ 4, 0, -(j ^ 2) * q ^ 2 + j ^ 2 * q ^ 3, j ^ 2 * q - j ^ 2 * q ^ 2, -(j ^ 2) * q ^ 2, 0, j ^ 2 * q, j ^ 2 - j ^ 2 * q]], [[q, 0, -(q ^ 2), 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -(q ^ 2), q, 0, 0, 0, 0, 0], [0, -(q ^ 2), 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1 + q, -q, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -(q ^ 2), -1 + q, 1], [0, 0, 0, 0, 0, 0, -(q ^ 2), 0, q, 0]], [[-1 + q, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, q, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -q, 0, 1 // q], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, -q, 0, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, -(q ^ 2), 0, q, 0], [0, -(q ^ 3), 0, 0, 0, q ^ 2, 0, 0, 0, -1 + q]], [[q, 0, -(q ^ 2), 0, 0, 0, 0, 0, 0, 0], [0, -1 + q, 0, q, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, q, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, q, 0, -1 // q], [0, 0, 0, 0, 0, 0, q, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, q, -1 + q, 0, 0, 0, 0, 0, 0], [q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, -(q ^ 2), 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, -(q ^ 2), 0, q, 0], [0, 0, 0, 0, 0, 0, 0, -(q ^ 2), 0, q]]]
                        end
                    f23 = function (q,)
                            return q ^ 0 * [[[-1, 0, 0, 0, 0], [0, -1, 0, 0, 0], [q ^ 3, -(q ^ 2), q, 0, 0], [0, 0, 0, -1, 0], [-(q ^ 3), 0, 0, -(q ^ 2), q]], [[-1, 0, 0, 0, 0], [0, -1, 0, 0, 0], [q ^ 3, -(q ^ 2), q, 0, 0], [0, 0, 0, -1, 0], [-(q ^ 3), 0, 0, -(q ^ 2), q]], [[-1, 0, 0, 0, 0], [0, 0, 1, 0, 0], [0, q, -1 + q, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, q, -1 + q]], [[0, 1, 0, 0, 0], [q, -1 + q, 0, 0, 0], [0, 0, -1, 0, 0], [0, 0, 0, -1, 0], [-(q ^ 4), q ^ 3, -(q ^ 2), -(q ^ 2), q]], [[-1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [0, q, 0, -1 + q, 0], [0, 0, q, 0, -1 + q]]]
                        end
                    f29 = function (q,)
                            return q ^ 0 * [[[q, q ^ 4, 0, 0, -(q ^ 2), 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, q ^ 3, -(q ^ 2), q, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, q ^ 4, 0, -(q ^ 3), q]], [[q, q ^ 4, 0, 0, -(q ^ 2), 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, q ^ 3, -(q ^ 2), q, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, q ^ 4, 0, -(q ^ 3), q]], [[-1 + q, 0, 0, 0, q, 0], [0, -1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, q, -1 + q, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q]], [[0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0], [0, q, -1 + q, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, q, 0], [q, 0, 0, 0, 0, -1 + q]], [[-1 + q, 0, 0, q, 0, 0], [0, q, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0], [0, 0, q, 0, -1 + q, 0], [0, 0, 0, 0, 0, -1]]]
                        end
                    f9 = function (q,)
                            return q ^ 0 * [[[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (E(3) + (E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, (E(3, 2) + (-(E(3, 2)) + E(3)) * q) - E(3) * q ^ 2, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0], [(-1 + 2q) - q ^ 2, (1 - 2q) + q ^ 2, (-1 + 2q) - q ^ 2, -1 + q, 0, -1 + q, 0, 1 - q, 0, E(3) - E(3) * q, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, (1 - 2q) + q ^ 2, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(-3 + q ^ -1 + 3q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, (-3 + q ^ -1 + 3q) - q ^ 2, -2 + q ^ -1 + q, 0, (E(3, 2) + 2 * E(3) + q ^ -1) - E(3) * q, 0, (2 - q ^ -1) - q, 0, -(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ((-(E(3, 2)) - 2 * E(3)) + E(3) * q ^ -1) - q, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, 0], [-(E(3)) + E(3) * q, E(3) - E(3) * q, 0, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, 0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, 0, E(3) - E(3) * q, 0, 0, 0, E(3, 2), 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, E(3, 2) - E(3, 2) * q, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q ^ -1, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, (E(3) + (-(E(3, 2)) - 2 * E(3)) * q) - q ^ 2, 0, ((3 - q ^ -1) - 3q) + q ^ 2, 0, 0, 0, 0, 0, 0, -1 + q, (1 + (E(3, 2) + 2 * E(3)) * q) - E(3) * q ^ 2, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (((E(3, 2) + 3 * E(3)) - E(3) * q ^ -1) + (-2 * E(3, 2) - 3 * E(3)) * q) - q ^ 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, -2 * E(3, 2) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, 0, 0, 0, 0, (E(3, 2) + (-(E(3, 2)) + E(3)) * q) - E(3) * q ^ 2, 0, 0, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (E(3) + (E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, 0, 0, 0], [0, 0, -(E(3)) + E(3) * q, 0, E(3) - E(3) * q, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, 0, E(3) - E(3) * q, E(3) - E(3) * q, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + q, 0, 0, -2 + q ^ -1 + q, 0, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, -1 + q, -2 + q ^ -1 + q, 0, 0, -2 + q ^ -1 + q, 0, 1 - q, 0, 0, 0, E(3) - E(3) * q, 1 - q, 0, 0, E(3), 0, 0, -1 + q ^ -1, 0, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1 - q ^ -1, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 1 - q, 0, 1, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -2 * E(3, 2) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, (2 * E(3, 2) + E(3) + q ^ -1) - E(3, 2) * q, 0, ((-2 * E(3, 2) - E(3)) + E(3, 2) * q ^ -1) - q, 0, 0, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, 0, 0, 0, ((2 * E(3, 2) - E(3, 2) * q ^ -1) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, 0, 0, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, -1 + q ^ -1 + q, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, 0, 0, 0], [(1 + (E(3, 2) + 2 * E(3)) * q) - E(3) * q ^ 2, 0, 0, (((-2 * E(3, 2) - 3 * E(3)) - q ^ -1) + (E(3, 2) + 3 * E(3)) * q) - E(3) * q ^ 2, 0, 0, 0, 0, 0, 0, (-3 + q ^ -1 + 3q) - q ^ 2, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, 0, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, (E(3) + (-(E(3, 2)) - 2 * E(3)) * q) - q ^ 2, 0, 0, 0, (E(3, 2) - E(3, 2) * q) + E(3, 2) * q ^ 2, (-1 + 2q) - q ^ 2, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, (2 - q ^ -1) - q, 0, 0, 0, 0, 0, 0, (-3 + q ^ -1 + 3q) - q ^ 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)), 0, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, 0, 0, 0, E(3) - E(3) * q, 0, q, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -q + q ^ 2, 0, 0, 1 - q, 0, 0, 0, (-(E(3)) + E(3) * q) - E(3) * q ^ 2, 0, (E(3) + (E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, -(E(3)) + (-(E(3, 2)) + E(3)) * q + E(3, 2) * q ^ 2, 0, 0, 0, 0], [0, 0, ((E(3, 2) - E(3)) + E(3) * q ^ -1) - E(3, 2) * q, 0, 0, -2 + q ^ -1 + q, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, 0, 0, 0, 0, ((3 - q ^ -1) - 3q) + q ^ 2, (-(E(3, 2)) + E(3) + E(3, 2) * q ^ -1) - E(3) * q, 0, 0, 0, 0, -2 + q ^ -1 + q, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, 0, 0, -1 + q ^ -1 + q, 0, E(3, 2) - E(3, 2) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, (-(E(3)) + E(3) * q) - E(3) * q ^ 2, (E(3) + (E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, 0, 0, 0, 0, 0, (1 - q) + q ^ 2, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, -(E(3)) + (-(E(3, 2)) + E(3)) * q + E(3, 2) * q ^ 2, 0, 0, 0], [(-1 + 2q) - q ^ 2, (((2 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-(E(3, 2)) + 2 * E(3)) * q) - E(3) * q ^ 2, E(3) + (E(3, 2) - 2 * E(3)) * q + (-2 * E(3, 2) + E(3)) * q ^ 2 + E(3, 2) * q ^ 3, ((((7 * E(3, 2) + E(3, 2) * q ^ -2) - 4 * E(3, 2) * q ^ -1) - 7 * E(3, 2) * q) + 4 * E(3, 2) * q ^ 2) - E(3, 2) * q ^ 3, 0, (((-(E(3, 2)) + 2 * E(3)) - E(3) * q ^ -1) + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, ((3 * E(3, 2) - E(3, 2) * q ^ -1) - 3 * E(3, 2) * q) + E(3, 2) * q ^ 2, ((2 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-2 * E(3, 2) + E(3)) * q + E(3, 2) * q ^ 2, 0, ((2 - q ^ -1) - 2q) + q ^ 2, (((-5 * E(3, 2) + 2 * E(3)) - E(3, 2) * q ^ -2) + (3 * E(3, 2) - E(3)) * q ^ -1 + (5 * E(3, 2) - E(3)) * q) - 2 * E(3, 2) * q ^ 2, (((3 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-4 * E(3, 2) + 2 * E(3)) * q + (3 * E(3, 2) - E(3)) * q ^ 2) - E(3, 2) * q ^ 3, (E(3, 2) - 2 * E(3)) + E(3, 2) * q ^ -2 + (-2 * E(3, 2) + E(3)) * q ^ -1 + E(3) * q, (-1 + 2q) - q ^ 2, ((-(E(3, 2)) + 2 * E(3, 2) * q) - 2 * E(3, 2) * q ^ 2) + E(3, 2) * q ^ 3, (2 * E(3, 2) - 2 * E(3)) + E(3, 2) * q ^ -2 + (-2 * E(3, 2) + E(3)) * q ^ -1 + (-2 * E(3, 2) + E(3)) * q + E(3, 2) * q ^ 2, 0, -2 * E(3, 2) + E(3) + E(3, 2) * q ^ -1 + (E(3, 2) - 2 * E(3)) * q + E(3) * q ^ 2, (-1 + 2q) - q ^ 2, (-2 * E(3, 2) + E(3) + E(3, 2) * q ^ -1 + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, 0, (-2 + q ^ -1 + 2q) - q ^ 2, -3 * E(3, 2) + E(3) + E(3, 2) * q ^ -1 + (4 * E(3, 2) - 2 * E(3)) * q + (-3 * E(3, 2) + E(3)) * q ^ 2 + E(3, 2) * q ^ 3, (-(E(3, 2)) - 3 * E(3)) + E(3) * q ^ -1 + (2 * E(3, 2) + 4 * E(3)) * q + (-(E(3, 2)) - 3 * E(3)) * q ^ 2 + E(3) * q ^ 3, 0, -2 + q ^ -1 + q, (E(3) - E(3) * q) + E(3) * q ^ 2, (-2 + q ^ -1 + 2q) - q ^ 2, ((2 * E(3, 2) - E(3)) + E(3, 2) * q ^ -2 + (-2 * E(3, 2) + E(3)) * q ^ -1) - E(3, 2) * q, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (-2 * E(3, 2) + E(3) + E(3, 2) * q ^ -1 + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, E(3) - E(3) * q, (1 - q) + q ^ 2, 0, 0, ((((-5 * E(3, 2) + 2 * E(3)) - E(3, 2) * q ^ -2) + (3 * E(3, 2) - E(3)) * q ^ -1 + (6 * E(3, 2) - E(3)) * q) - 4 * E(3, 2) * q ^ 2) + E(3, 2) * q ^ 3, 0, (((2 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-2 * E(3, 2) + 2 * E(3)) * q + (2 * E(3, 2) - E(3)) * q ^ 2) - E(3, 2) * q ^ 3, 0, 0], [0, -2 + q ^ -1 + q, (1 - 2q) + q ^ 2, ((-4 - q ^ -2) + 3 * q ^ -1 + 3q) - q ^ 2, 0, (2 - q ^ -1) - q, -2 + q ^ -1 + q, -2 + q ^ -1 + q, 0, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, (((-3 * E(3, 2) - 4 * E(3)) + q ^ -2) - 3 * q ^ -1) + (E(3, 2) + 2 * E(3)) * q, (-3 + q ^ -1 + 3q) - q ^ 2, (-1 - q ^ -2) + 2 * q ^ -1, E(3) - E(3) * q, (1 - 2q) + q ^ 2, ((3 * E(3, 2) + 2 * E(3)) - q ^ -2) + (-3 * E(3, 2) - 2 * E(3)) * q ^ -1 + q, 0, (2 - q ^ -1) - q, E(3) - E(3) * q, (2 - q ^ -1) - q, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, ((3 - q ^ -1) - 3q) + q ^ 2, ((2 * E(3, 2) - E(3)) - E(3, 2) * q ^ -1) + (-2 * E(3, 2) + E(3)) * q + E(3, 2) * q ^ 2, 0, E(3) - E(3) * q ^ -1, 1 - q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, (-1 - q ^ -2) + 2 * q ^ -1, -1 + q, (2 - q ^ -1) - q, 1, E(3, 2) - E(3, 2) * q, 1, 0, (((4 + q ^ -2) - 3 * q ^ -1) - 3q) + q ^ 2, 0, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)) * q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [-(E(3, 2)) + E(3, 2) * q, E(3, 2) - E(3, 2) * q, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0], [0, 0, 0, -1 + q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 1 - q, 0, E(3) - E(3) * q, E(3) * q], [((3 - q ^ -1) - 3q) + q ^ 2, (-3 + q ^ -1 + 3q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, -(E(3, 2)) + E(3, 2) * q, 0, -(E(3, 2)) + E(3, 2) * q ^ -1, 0, -2 + q ^ -1 + q, -1 + q, (E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, ((E(3, 2) + 2 * E(3)) - E(3) * q ^ -1) + q, ((-(E(3, 2)) - 2 * E(3)) + E(3) * q ^ -1) - q, (-2 + q ^ -1 + 2q) - q ^ 2, -(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, E(3, 2) - E(3, 2) * q]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -q, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0], [1 - q, 0, 1 - q, 0, 1 - q ^ -1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -2 + q ^ -1 + q, 0, 0, -1 + q ^ -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, (2 - q ^ -1) - q, (-1 + 2q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, -2 + q ^ -1 + q, 0, (2 - q ^ -1) - q, (2 - q ^ -1) - q, -(E(3)) + E(3) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, E(3) - E(3) * q, (1 - 2q) + q ^ 2, 0, 0, (-1 + 2q) - q ^ 2, 0, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, -(E(3)) + E(3) * q, -2 + q ^ -1 + q, -2 + q ^ -1 + q, 0, (-1 + 2q) - q ^ 2, 0, E(3, 2) - E(3, 2) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, -1, -2 + q ^ -1 + q, (-1 + 2q) - q ^ 2, -2 + q ^ -1 + q, ((-3 * E(3, 2) - 2 * E(3)) - q ^ -1) + (3 * E(3, 2) + 2 * E(3)) * q + q ^ 2, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, E(3) - E(3) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, -1 + q, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0], [(2 * E(3) - E(3) * q ^ -1) - E(3) * q, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, (-3 * E(3) + E(3) * q ^ -1 + 3 * E(3) * q) - E(3) * q ^ 2, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, E(3, 2) - E(3, 2) * q, -2 * E(3, 2) + E(3, 2) * q ^ -1 + E(3, 2) * q, -(E(3, 2)) + E(3, 2) * q, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, 0, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, ((E(3, 2) + 3 * E(3)) - E(3) * q ^ -1) + (-(E(3, 2)) - 3 * E(3)) * q + E(3) * q ^ 2, 0, E(3, 2) - E(3, 2) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, 0, -1 + q, 0, E(3) - E(3) * q, 0, 0, 0, 0, E(3), 0, E(3), (2 * E(3) - E(3) * q ^ -1) - E(3) * q, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, ((-(E(3, 2)) - 3 * E(3)) + E(3) * q ^ -1 + (E(3, 2) + 3 * E(3)) * q) - E(3) * q ^ 2, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, -(E(3, 2)) + E(3, 2) * q], [(-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, (-1 + 2q) - q ^ 2, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, (2 - q ^ -1) - q, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, E(3) - E(3) * q, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, -(E(3)) + E(3) * q, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (-1 + 2q) - q ^ 2, 0, 0, -q + q ^ 2, 0, (1 - 2q) + q ^ 2, 0, 0, -(E(3)) + E(3) * q, 1 - q, 0, E(3) * q - E(3) * q ^ 2, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (1 - 2q) + q ^ 2, (-1 + 2q) - q ^ 2, E(3) - E(3) * q, 0], [1 - q, 0, 0, 0, 1 - q, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1 + q]], [[0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, (2 - q ^ -1) - q, (-1 + 2q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, -2 + q ^ -1 + q, 0, (2 - q ^ -1) - q, (2 - q ^ -1) - q, -(E(3)) + E(3) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, E(3) - E(3) * q, (1 - 2q) + q ^ 2, 0, 0, (-1 + 2q) - q ^ 2, 0, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, -(E(3)) + E(3) * q, -2 + q ^ -1 + q, -2 + q ^ -1 + q, 0, (-1 + 2q) - q ^ 2, 0, E(3, 2) - E(3, 2) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, -1, -2 + q ^ -1 + q, (-1 + 2q) - q ^ 2, -2 + q ^ -1 + q, ((-3 * E(3, 2) - 2 * E(3)) - q ^ -1) + (3 * E(3, 2) + 2 * E(3)) * q + q ^ 2, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, E(3) - E(3) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [(-3 * E(3, 2) + E(3, 2) * q ^ -1 + (4 * E(3, 2) + E(3)) * q + (-3 * E(3, 2) - 2 * E(3)) * q ^ 2) - q ^ 3, (3 * E(3, 2) - E(3, 2) * q ^ -1) + (-4 * E(3, 2) - E(3)) * q + (3 * E(3, 2) + 2 * E(3)) * q ^ 2 + q ^ 3, ((-4 * E(3, 2) - 2 * E(3)) - q ^ -1) + (6 * E(3, 2) + 2 * E(3)) * q + (-4 * E(3, 2) - E(3)) * q ^ 2 + E(3, 2) * q ^ 3, (E(3, 2) - 3 * E(3)) + E(3) * q ^ -1 + (-2 * E(3, 2) + 4 * E(3)) * q + (E(3, 2) - 3 * E(3)) * q ^ 2 + E(3) * q ^ 3, ((2 * E(3) - E(3) * q ^ -1) - 2 * E(3) * q) + E(3) * q ^ 2, (2 - q ^ -1) - q, (-2 * E(3) + E(3) * q ^ -1 + 2 * E(3) * q) - E(3) * q ^ 2, ((3 * E(3, 2) - E(3, 2) * q ^ -1) - 3 * E(3, 2) * q) + E(3, 2) * q ^ 2, 0, ((3 - q ^ -1) - 3q) + q ^ 2, 0, E(3, 2) + (-2 * E(3, 2) + E(3)) * q + (E(3, 2) - 2 * E(3)) * q ^ 2 + E(3) * q ^ 3, 0, E(3) - E(3) * q, (-(E(3, 2)) + (2 * E(3, 2) - E(3)) * q + (-(E(3, 2)) + 2 * E(3)) * q ^ 2) - E(3) * q ^ 3, 0, (((3 * E(3) - E(3) * q ^ -1) - 4 * E(3) * q) + 3 * E(3) * q ^ 2) - E(3) * q ^ 3, 0, -1 + (-2 * E(3, 2) - E(3)) * q + E(3, 2) * q ^ 2, ((2 * E(3) - E(3) * q ^ -1) - 2 * E(3) * q) + E(3) * q ^ 2, ((2 * E(3) - E(3) * q ^ -1) - 2 * E(3) * q) + E(3) * q ^ 2, 0, (-(E(3, 2)) + (2 * E(3, 2) - E(3)) * q + (-(E(3, 2)) + 2 * E(3)) * q ^ 2) - E(3) * q ^ 3, 0, (-1 + q) - q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, (((-(E(3, 2)) + 2 * E(3)) - E(3) * q ^ -1) + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, ((E(3) - 3 * E(3) * q) + 3 * E(3) * q ^ 2) - E(3) * q ^ 3, (((-(E(3, 2)) + 2 * E(3)) - E(3) * q ^ -1) + (2 * E(3, 2) - E(3)) * q) - E(3, 2) * q ^ 2, (((4 * E(3, 2) - E(3, 2) * q ^ -1) - 6 * E(3, 2) * q) + 4 * E(3, 2) * q ^ 2) - E(3, 2) * q ^ 3, (((3 * E(3, 2) + E(3)) - E(3, 2) * q ^ -1) + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, -(E(3, 2)) + (2 * E(3, 2) + E(3)) * q + q ^ 2], [q - q ^ 2, -q + q ^ 2, q - q ^ 2, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, -q + q ^ 2, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q]], [[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q - q ^ 2, 0, q - q ^ 2, 0, -1 + q, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1 - q ^ -1, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 1 - q, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(-1 + 2q) - q ^ 2, (1 - 2q) + q ^ 2, (q - 2 * q ^ 2) + q ^ 3, ((1 - 3q) + 3 * q ^ 2) - q ^ 3, (-1 + 2q) - q ^ 2, 0, (1 - 2q) + q ^ 2, (1 - 2q) + q ^ 2, E(3) * q - E(3) * q ^ 2, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, -(E(3)) * q + E(3) * q ^ 2, (-q + 2 * q ^ 2) - q ^ 3, 0, 0, (q - 2 * q ^ 2) + q ^ 3, 0, -1 + (-3 * E(3, 2) - 2 * E(3)) * q + (3 * E(3, 2) + 2 * E(3)) * q ^ 2 + q ^ 3, 0, E(3) * q - E(3) * q ^ 2, (-1 + 2q) - q ^ 2, (-1 + 2q) - q ^ 2, 0, (q - 2 * q ^ 2) + q ^ 3, 0, -(E(3, 2)) * q + E(3, 2) * q ^ 2, 0, q - q ^ 2, 0, 0, 0, 0, q, 0, q, (-1 + 2q) - q ^ 2, (q - 2 * q ^ 2) + q ^ 3, (-1 + 2q) - q ^ 2, (1 + (3 * E(3, 2) + 2 * E(3)) * q + (-3 * E(3, 2) - 2 * E(3)) * q ^ 2) - q ^ 3, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, -(E(3)) * q + E(3) * q ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [-1 + q, 1 - q, -1 + q, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 1 - q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]], [[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [-1 + q, 1 - q, -1 + q, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [q - q ^ 2, 0, q - q ^ 2, 0, -1 + q, 0, 0, -q, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, (2 - q ^ -1) - q, (-1 + 2q) - q ^ 2, ((3 - q ^ -1) - 3q) + q ^ 2, -2 + q ^ -1 + q, 0, (2 - q ^ -1) - q, (2 - q ^ -1) - q, -(E(3)) + E(3) * q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, E(3) - E(3) * q, (1 - 2q) + q ^ 2, 0, 0, (-1 + 2q) - q ^ 2, 0, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, -(E(3)) + E(3) * q, -2 + q ^ -1 + q, -2 + q ^ -1 + q, 0, (-1 + 2q) - q ^ 2, 0, E(3, 2) - E(3, 2) * q, 0, -1 + q, 0, 0, 0, 0, -1 + q, 0, -1, -2 + q ^ -1 + q, (-1 + 2q) - q ^ 2, -2 + q ^ -1 + q, ((-3 * E(3, 2) - 2 * E(3)) - q ^ -1) + (3 * E(3, 2) + 2 * E(3)) * q + q ^ 2, -2 * E(3) + E(3) * q ^ -1 + E(3) * q, E(3) - E(3) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [0, (-1 + 2q) - q ^ 2, (-q + 2 * q ^ 2) - q ^ 3, ((-3 + q ^ -1 + 4q) - 3 * q ^ 2) + q ^ 3, 0, (1 - 2q) + q ^ 2, (-1 + 2q) - q ^ 2, (-1 + 2q) - q ^ 2, 0, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, (3 - q ^ -1) + (3 * E(3, 2) + 4 * E(3)) * q + (-(E(3, 2)) - 2 * E(3)) * q ^ 2, ((-1 + 3q) - 3 * q ^ 2) + q ^ 3, -2 + q ^ -1 + q, -(E(3)) * q + E(3) * q ^ 2, (-q + 2 * q ^ 2) - q ^ 3, (3 * E(3, 2) + 2 * E(3) + q ^ -1 + (-3 * E(3, 2) - 2 * E(3)) * q) - q ^ 2, 0, (1 - 2q) + q ^ 2, -(E(3)) * q + E(3) * q ^ 2, (1 - 2q) + q ^ 2, 0, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, ((1 - 3q) + 3 * q ^ 2) - q ^ 3, (E(3, 2) + (-2 * E(3, 2) + E(3)) * q + (2 * E(3, 2) - E(3)) * q ^ 2) - E(3, 2) * q ^ 3, 0, E(3) - E(3) * q, -q + q ^ 2, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, -2 + q ^ -1 + q, q - q ^ 2, (1 - 2q) + q ^ 2, -q, -(E(3, 2)) * q + E(3, 2) * q ^ 2, 0, 0, (((3 - q ^ -1) - 4q) + 3 * q ^ 2) - q ^ 3, 0, -1 + (-3 * E(3, 2) - 2 * E(3)) * q + (3 * E(3, 2) + 2 * E(3)) * q ^ 2 + q ^ 3, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]]
                        end
                    f8 = function (q,)
                            return q ^ 0 * [[[E(3) - E(3) * q, E(3) * q, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, E(3, 2) - E(3, 2) * q, 0, -(E(3)) + (-3 + ER(-3)) // 2 * q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3), -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, 1 - q, -(E(3, 2)) + E(3, 2) * q ^ -1, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, E(3, 2) - E(3, 2) * q, (((5 + ER(-3)) // 2 - q ^ -1) + (-2 - ER(-3)) * q) - E(3, 2) * q ^ 2, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0], [((2 * E(3, 2) - E(3, 2) * q ^ -1) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, ((-3 + ER(-3)) // 2 + q ^ -1) - E(3) * q, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, 0, 0, 0, (1 - q) + q ^ 2, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, (E(3) - ER(-3) * q) - E(3, 2) * q ^ 2, (((3 * E(3, 2) - E(3, 2) * q ^ -1) - 4 * E(3, 2) * q) + 3 * E(3, 2) * q ^ 2) - E(3, 2) * q ^ 3, 0, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, (E(3, 2) - 2 * E(3, 2) * q) + E(3, 2) * q ^ 2, 0, 0, 0, 0, 0, 0, 0], [0, 0, E(3, 2) - E(3, 2) * q ^ -1, 0, (1 - 2q) + q ^ 2, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0], [0, 0, (((-2 + ER(-3)) - E(3) * q ^ -1) + (5 - ER(-3)) // 2 * q) - q ^ 2, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0], [0, 0, (((-2 + ER(-3)) - E(3) * q ^ -1) + (5 - ER(-3)) // 2 * q) - q ^ 2, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0], [0, 0, (ER(-3) - E(3) * q ^ -1) + E(3, 2) * q, 0, E(3, 2) + (2 + ER(-3)) * q + (-5 - ER(-3)) // 2 * q ^ 2 + q ^ 3, 0, 0, 0, 0, 0, 0, (E(3) - E(3) * q) + E(3) * q ^ 2, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, -1 + q, 0, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, 0, 0, 0], [0, 0, ((-4 - q ^ -2) + 3 * q ^ -1 + 3q) - q ^ 2, (ER(-3) - E(3) * q ^ -1) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q ^ -1 + E(3) * q, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, (2 - q ^ -1) - q, 0, 0], [((3 - q ^ -1) - 3q) + q ^ 2, (-1 + 2q) - q ^ 2, 0, 0, ((-1 + 3q) - 3 * q ^ 2) + q ^ 3, (1 - 2q) + q ^ 2, (E(3, 2) * q - 2 * E(3, 2) * q ^ 2) + E(3, 2) * q ^ 3, (E(3) - ER(-3) * q) - E(3, 2) * q ^ 2, 0, -1 + q, 0, (E(3) - ER(-3) * q) - E(3, 2) * q ^ 2, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, -1, 0, 0, 0, 0, 0, (-1 + 2q) - q ^ 2, 0, 1 - q, E(3, 2) * q - E(3, 2) * q ^ 2, 0, 0, 0, 0], [0, 0, ((-4 * E(3) - E(3) * q ^ -2) + 3 * E(3) * q ^ -1 + 3 * E(3) * q) - E(3) * q ^ 2, 0, ((-3 + ER(-3)) // 2 + q ^ -1) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q ^ -1 + E(3, 2) * q, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0], [-2 + q ^ -1 + q, 0, 0, (-3 + q ^ -1 + 3q) - q ^ 2, (((1 - 3 * ER(-3)) + ER(-3) * q ^ -1 + (-5 + 7 * ER(-3)) // 2 * q) - 4 * E(3) * q ^ 2) + E(3) * q ^ 3, (((-1 + 2 * ER(-3)) - ER(-3) * q ^ -1) - 3 * E(3) * q) + E(3) * q ^ 2, (((-3 + ER(-3)) // 2 + (7 - ER(-3)) // 2 * q) - 3 * q ^ 2) + q ^ 3, (((4 * E(3, 2) - E(3, 2) * q ^ -1) + (2 + 3 * ER(-3)) * q) - 2 * ER(-3) * q ^ 2) + E(3) * q ^ 3, 0, (((5 + ER(-3)) // 2 - q ^ -1) + (-2 - ER(-3)) * q) - E(3, 2) * q ^ 2, 0, (((-2 - ER(-3)) - E(3, 2) * q ^ -1) + (5 + ER(-3)) // 2 * q) - q ^ 2, ((-3 + ER(-3)) // 2 + q ^ -1) - E(3) * q, ((3 - q ^ -1) - 3q) + q ^ 2, -2 + q ^ -1 + q, (-1 + 2q) - q ^ 2, 0, -1 + q, -1, 0, (-(E(3, 2)) - ER(-3) * q) + E(3) * q ^ 2, 0, (-3 * E(3) + E(3) * q ^ -1 + 3 * E(3) * q) - E(3) * q ^ 2, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, (-1 + 2q) - q ^ 2, (1 - 2q) + q ^ 2, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, 0, -(E(3)) + E(3) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, 0, 0, 0], [(-4 + 2 * q ^ -1 + 3q) - q ^ 2, (1 - 2q) + q ^ 2, 0, -2 + q ^ -1 + q, 0, 0, 0, (((-1 - 2 * ER(-3)) - E(3, 2) * q ^ -1) + (1 + 5 * ER(-3)) // 2 * q) - ER(-3) * q ^ 2, 0, ((5 + ER(-3)) // 2 - q ^ -1) + (-3 - ER(-3)) // 2 * q, 0, 0, 0, (2 - q ^ -1) - q, 0, -1 + q, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, 0, 0, 0, 1 - q, 0, 0, -(E(3))], [0, 0, 0, 0, ((2 - ER(-3)) + E(3) * q ^ -1 + (-7 + ER(-3)) // 2 * q + 3 * q ^ 2) - q ^ 3, 0, 0, 0, 0, 0, 0, (((-1 - 2 * ER(-3)) - E(3, 2) * q ^ -1) + (1 + 5 * ER(-3)) // 2 * q) - ER(-3) * q ^ 2, -2 + q ^ -1 + q, 0, (2 * E(3, 2) - E(3, 2) * q ^ -1) - E(3, 2) * q, 0, 0, -1 + q, 0, -1 + q ^ -1, 0, E(3) - E(3) * q, (1 - 2q) + q ^ 2, 0, 0, 0, -1 + q, 0, -(E(3, 2)), 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [5 + 2 * ER(-3) + (-3 - ER(-3)) // 2 * q ^ -1 + (-6 - 3 * ER(-3)) * q + (3 + 2 * ER(-3)) * q ^ 2 + E(3, 2) * q ^ 3, ((E(3, 2) - 3 * E(3, 2) * q) + 3 * E(3, 2) * q ^ 2) - E(3, 2) * q ^ 3, 0, (((4 - q ^ -1) - 6q) + 4 * q ^ 2) - q ^ 3, (((-3 + 7 * ER(-3)) // 2 - ER(-3) * q ^ -1) + (11 - 9 * ER(-3)) // 2 * q + (-15 + 5 * ER(-3)) // 2 * q ^ 2 + (9 - ER(-3)) // 2 * q ^ 3) - q ^ 4, ((3 - 5 * ER(-3)) // 2 + ER(-3) * q ^ -1 + (-4 + 2 * ER(-3)) * q + (7 - ER(-3)) // 2 * q ^ 2) - q ^ 3, ((3 - ER(-3)) // 2 + (-9 + ER(-3)) // 2 * q + (5 + ER(-3)) * q ^ 2 + (-5 - 3 * ER(-3)) // 2 * q ^ 3) - E(3, 2) * q ^ 4, (3 + 5 * ER(-3)) // 2 + E(3, 2) * q ^ -1 + (-3 - 11 * ER(-3)) // 2 * q + 6 * ER(-3) * q ^ 2 + (1 - 3 * ER(-3)) * q ^ 3 + E(3) * q ^ 4, (1 + (-3 + ER(-3)) // 2 * q) - E(3) * q ^ 2, ((-4 - ER(-3)) + q ^ -1 + (11 + 5 * ER(-3)) // 2 * q + (-3 - 2 * ER(-3)) * q ^ 2) - E(3, 2) * q ^ 3, -(E(3)) + E(3) * q, ((3 + 5 * ER(-3)) // 2 + E(3, 2) * q ^ -1 + (-3 - 9 * ER(-3)) // 2 * q + (1 + 7 * ER(-3)) // 2 * q ^ 2) - ER(-3) * q ^ 3, ((3 - q ^ -1) - 3q) + q ^ 2, -3 + q ^ -1 + (9 + ER(-3)) // 2 * q + (-7 - ER(-3)) // 2 * q ^ 2 + q ^ 3, (-3 * E(3, 2) + E(3, 2) * q ^ -1 + 3 * E(3, 2) * q) - E(3, 2) * q ^ 2, (1 + (-5 - ER(-3)) // 2 * q + (5 + ER(-3)) // 2 * q ^ 2) - q ^ 3, 0, (1 - 2q) + q ^ 2, 0, (2 - q ^ -1) - q, E(3, 2) + (1 + 3 * ER(-3)) // 2 * q + (1 - 3 * ER(-3)) // 2 * q ^ 2 + E(3) * q ^ 3, (E(3, 2) + ER(-3) * q) - E(3) * q ^ 2, (-3 + ER(-3)) // 2 + (4 - ER(-3)) * q + (-7 + ER(-3)) // 2 * q ^ 2 + q ^ 3, -1, (2 - ER(-3)) + E(3) * q ^ -1 + (-5 + ER(-3)) // 2 * q + q ^ 2, 1 + (-5 - ER(-3)) // 2 * q + (2 + ER(-3)) * q ^ 2 + E(3, 2) * q ^ 3, (q - 2 * q ^ 2) + q ^ 3, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, E(3, 2) - E(3, 2) * q, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, -1 + q, 0, 0], [0, 0, (((-4 - q ^ -2) + 3 * q ^ -1 + 4q) - 3 * q ^ 2) + q ^ 3, 0, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, 0, 0, 0, 0, 0, 0, (((3 * E(3) - E(3) * q ^ -1) - 4 * E(3) * q) + 3 * E(3) * q ^ 2) - E(3) * q ^ 3, ((2 - q ^ -1) - 2q) + q ^ 2, 0, (2 * E(3) - E(3) * q ^ -1) - E(3) * q, 0, 0, (1 - 2q) + q ^ 2, 0, 0, 0, (-(E(3)) + E(3) * q) - E(3) * q ^ 2, E(3, 2) - E(3, 2) * q, 0, 0, 0, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, ((2 - q ^ -1) - 2q) + q ^ 2, E(3, 2) - E(3, 2) * q, 0], [((-5 + 2 * q ^ -1 + 5q) - 3 * q ^ 2) + q ^ 3, ((2 - 4q) + 3 * q ^ 2) - q ^ 3, (1 - 2q) + q ^ 2, ((-(ER(-3)) - E(3, 2) * q ^ -1) + (-5 + ER(-3)) // 2 * q + 3 * q ^ 2) - q ^ 3, 0, 0, 0, ((-1 - 2 * ER(-3)) - E(3, 2) * q ^ -1) + (1 + 3 * ER(-3)) * q + (-1 - 5 * ER(-3)) // 2 * q ^ 2 + ER(-3) * q ^ 3, 0, ((7 + ER(-3)) // 2 - q ^ -1) + (-4 - ER(-3)) * q + (3 + ER(-3)) // 2 * q ^ 2, 0, 0, 0, ((2 - q ^ -1) - 2q) + q ^ 2, 0, (-1 + 2q) - q ^ 2, 0, 0, 0, 0, (-(E(3, 2)) + E(3, 2) * q) - E(3, 2) * q ^ 2, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, E(3, 2) - E(3, 2) * q]], [[0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, q - q ^ 2, 0, 0, -(E(3)) + E(3) * q, -1 + q, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1 + q, 0, 0, 0, E(3) - E(3) * q, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, q - q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)), 0], [0, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [((-2 + ER(-3)) - E(3) * q ^ -1) + (4 - ER(-3)) * q + (-7 + ER(-3)) // 2 * q ^ 2 + q ^ 3, (-(E(3)) + (-2 + ER(-3)) * q + (5 - ER(-3)) // 2 * q ^ 2) - q ^ 3, 0, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, 1 + 2 * ER(-3) + E(3, 2) * q ^ -1 + (-1 - 5 * ER(-3)) // 2 * q + ER(-3) * q ^ 2, 0, -(E(3)) + (-2 + ER(-3)) * q + (3 - ER(-3)) // 2 * q ^ 2, 0, (-(E(3)) + 2 * E(3) * q) - E(3) * q ^ 2, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, -1 + q, E(3) * q, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, 1 - q, 0, (1 - 2q) + q ^ 2, (1 - 2q) + q ^ 2, (-5 - ER(-3)) // 2 + q ^ -1 + (3 + ER(-3)) // 2 * q, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (2 - q ^ -1) - q, 0, 0, 1 - q, 0, E(3, 2), 0, 0, 0, 0, 0], [0, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, E(3), 0, 0], [0, 0, q - q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) * q - E(3, 2) * q ^ 2, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0, -1 + q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q]], [[0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, q, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, E(3) * q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -(E(3)) + E(3) * q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, 0, 0, (1 - 2q) + q ^ 2, (E(3) - 2 * E(3) * q) + E(3) * q ^ 2, (-5 - ER(-3)) // 2 + q ^ -1 + (2 + ER(-3)) * q + E(3, 2) * q ^ 2, -q + q ^ 2, 0, (1 - 2q) + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, -q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, ((3 + ER(-3)) // 2 - q ^ -1) + E(3, 2) * q, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3))], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q ^ -1, 0, -(E(3, 2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) - E(3) * q, 0, 0, -(E(3)) * q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1 + q, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q + q ^ 2, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, E(3) * q, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, -1 + q, 0, 0, 0, 0, 0], [0, 0, 0, -1 + q, (2 - ER(-3)) + E(3) * q ^ -1 + (-5 + ER(-3)) // 2 * q + q ^ 2, 0, -1 + q, 0, -1 + q, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, 0, -1 + q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q ^ -1, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, (-(E(3, 2)) - ER(-3) * q) + E(3) * q ^ 2, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1 - q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q ^ -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, E(3, 2) * q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1 - q, (((-2 + ER(-3)) - E(3) * q ^ -1) + (5 - ER(-3)) // 2 * q) - q ^ 2, 0, 0, 0, 1 - q, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3, 2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, q - q ^ 2, 0, 0, -(E(3)) + E(3) * q, -1 + q, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, E(3) * q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)) * q, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 - q, 0, -1 + q, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0], [1 - q, q, 0, -(E(3, 2)) * q + E(3, 2) * q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q - q ^ 2, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, q - q ^ 2, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0], [0, 0, 0, (-1 + 2q) - q ^ 2, ((5 - 3 * ER(-3)) // 2 + E(3) * q ^ -1 + (-9 + 3 * ER(-3)) // 2 * q + (7 - ER(-3)) // 2 * q ^ 2) - q ^ 3, 0, q - q ^ 2, 0, (-1 + 2q) - q ^ 2, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1 - 2q) + q ^ 2, 0, 0, q, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, q - q ^ 2, 0, 0, 0, 0, 0, -q + q ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, E(3, 2) * q], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0]], [[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1 + q, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2 + q ^ -1 + q, 1 - q, -(E(3, 2)) + E(3, 2) * q ^ -1, 0, 0, 0, 0, E(3, 2) - E(3, 2) * q, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, E(3), 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [(1 - 2q) + q ^ 2, q - q ^ 2, 0, 0, -(E(3)) + E(3) * q, -1 + q, q, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1 + q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0, -(E(3)), 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -q + q ^ 2, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3)), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 + q, 0, -(E(3, 2)) * q + E(3, 2) * q ^ 2, 0, 0, 0, 0, q, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0], [((4 - 2 * q ^ -1) - 3q) + q ^ 2, (-1 + 2q) - q ^ 2, 0, (2 - q ^ -1) - q, 0, 0, 0, 1 + 2 * ER(-3) + E(3, 2) * q ^ -1 + (-1 - 5 * ER(-3)) // 2 * q + ER(-3) * q ^ 2, 0, (-5 - ER(-3)) // 2 + q ^ -1 + (3 + ER(-3)) // 2 * q, 0, 0, 0, -2 + q ^ -1 + q, 0, 1 - q, 0, 0, 0, 0, -(E(3)) + E(3) * q, 0, 0, 0, 0, 0, -1 + q, 0, 0, E(3)], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1 + q, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) * q, 0, 0, 0, 0, -1 + q, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q, 0], [((-1 - 3 * ER(-3)) // 2 - E(3, 2) * q ^ -1) + (-1 + 5 * ER(-3)) // 2 * q + (1 - 2 * ER(-3)) * q ^ 2 + E(3) * q ^ 3, (-(E(3, 2)) + (-1 - 3 * ER(-3)) // 2 * q + (-1 + 3 * ER(-3)) // 2 * q ^ 2) - E(3) * q ^ 3, 0, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, (-7 - ER(-3)) // 2 + q ^ -1 + (4 + ER(-3)) * q + (-3 - ER(-3)) // 2 * q ^ 2, 0, -(E(3, 2)) + (-1 - 3 * ER(-3)) // 2 * q + ER(-3) * q ^ 2, 0, (-(E(3, 2)) + 2 * E(3, 2) * q) - E(3, 2) * q ^ 2, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, 0, 0, 0, 0, -(E(3, 2)) + E(3, 2) * q, 0, 0, 0, 0, -(E(3)) + E(3) * q, E(3, 2) * q, 0, 0, 0, -1 + q]]]
                        end
                    r = [f1(x), f2(x), f3(x, E(3)), f4(x, E(3)), f4(x, E(3, 2)), f3(x, E(3, 2)), [[[-1]], [[-1]], [[-1]], [[-1]], [[-1]]], f8(x), f9(x), f8(1 // x) * -x, f11(x, E(3)), f12(x, E(3)), f13(x, E(3)), f13(x, E(3, 2)), f12(x, E(3, 2)), f11(x, E(3, 2)), f17(x), f1(1 // x) * -x, f13(1 // x, E(3, 2)) * -x, f20(x, E(3)), f20(x, E(3, 2)), f13(1 // x, E(3)) * -x, f23(x), f2(1 // x) * -x, f11(1 // x, E(3, 2)) * -x, f12(1 // x, E(3)) * -x, f12(1 // x, E(3, 2)) * -x, f11(1 // x, E(3)) * -x, f29(x), f4(1 // x, E(3, 2)) * -x, f4(1 // x, E(3)) * -x, f23(1 // x) * -x, f3(1 // x, E(3, 2)) * -x, f3(1 // x, E(3)) * -x, f17(1 // x) * -x, [[[x]], [[x]], [[x]], [[x]], [[x]]]]
                    return r[i]
                end
            return -((para[2])[2]) * m(i)
        elseif [p, q, r] == [4, 4, 3]
            x = -((para[2])[1]) // (para[2])[2]
            r = x ^ 0 * [[[[x, -1, -1, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, -1 + x, -1, x, 0], [0, 1, -1, -1, 0, x]], [[-1, 0, 0, 0, 0, 0], [-x, x, 0, 0, 0, x], [0, 0, 0, 0, -x, 0], [0, 0, 0, x, 0, 0], [0, 0, -1, 0, -1 + x, 0], [0, 0, 0, 0, 0, -1]], [[x, -1, 0, 0, 0, -1], [0, -1, 0, 0, 0, 0], [0, 0, x, 0, 1, -1], [0, -x, 0, x, -1, 1 - x], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]]], [[[-1, 0, 0], [0, 0, 1], [0, x, -1 + x]], [[x, 0, 0], [-1, -1, 0], [1, 0, -1]], [[0, 1, 0], [x, -1 + x, 0], [0, 0, -1]]], [[[x, 0, 0], [x, -1, 0], [x, 0, -1]], [[-1, 2, 0], [0, x, 0], [0, (-(E(4)) + 1) * x, -1]], [[-1, 0, 1], [0, -1, (E(4) + 1) // 2], [0, 0, x]]], [[[x, 0, 0], [x, -1, 0], [x, 0, -1]], [[-1, 2, 0], [0, x, 0], [0, (E(4) + 1) * x, -1]], [[-1, 0, 1], [0, -1, (-(E(4)) + 1) // 2], [0, 0, x]]], [[[-1]], [[-1]], [[-1]]], [[[x, -1, 0], [0, -1, 0], [0, -1, x]], [[-1 + x, 0, 1], [0, x, 0], [x, 0, 0]], [[0, x, 0], [1, -1 + x, 0], [0, 0, x]]], [[[-1, 0, 0], [-1, x, 0], [-1, 0, x]], [[x, -2x, 0], [0, -1, 0], [0, -(E(4)) - 1, x]], [[x, 0, -x], [0, x, (E(4) - 1) // 2 * x], [0, 0, -1]]], [[[-1, 0, 0], [-1, x, 0], [-1, 0, x]], [[x, -2x, 0], [0, -1, 0], [0, E(4) - 1, x]], [[x, 0, -x], [0, x, (-(E(4)) - 1) // 2 * x], [0, 0, -1]]], [[[-1, 0], [-1, x]], [[-1, 0], [-1, x]], [[x, -x], [0, -1]]], [[[x]], [[x]], [[x]]]]
            return -((para[2])[2]) * r[i]
        else
            S = (((chevieget(:imp, :CharInfo))(p, q, r))[:charparams])[i]
            p1rRep = function ()
                    local Q, pos, ct
                    if r > 1
                        Q = -((para[2])[1]) // (para[2])[2]
                    else
                        Q = 0
                    end
                    pos = function (t, i)
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
                                (para[1])[p[1]] * Q ^ (p[3] - p[2])
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
                if para[2] != para[3]
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
                                                end), 1:div(q, 2)))
                            else
                                X = para[1]
                            end
                            X = X[S[[3, 4]]]
                            v = S[2] * GetRoot(Product(X) * Product(Y) * Product(T) * E(e, (2 - S[3]) - S[4]), 2) * E(p, (S[3] + S[4]) - 2)
                            d = 1 + Sum(X) * 0 + Sum(Y) * 0 + Sum(T) * 0
                            return [(d*[[X[1],Sum(Y,y->1//y)-X[2]//v*Sum(T)], 
                 [0, X[2]]])^div(q,2),
        [[Sum(Y),1//X[1]],[-Product(Y)*X[1],0]],[[0,-Product(T)//v],[v,Sum(T)]]]
                        end
                    else
                        error("should not happen")
                    end
                elseif para[1] == map((i->begin E(e, i - 1) end), 1:e)
                    para = [map((i->begin E(p, i) end), 0:p - 1), para[2]]
                else
                    para = [Concatenation(TransposedMat(map((i->begin
                              map(j->E(q,j), 0:q - 1) * GetRoot(i, q)
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
                (uc[:families])[6] = Family(ComplexConjugate(((CHEVIE[:families])[:X])(3)), [8, 7, 11], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[4] = Family(((CHEVIE[:families])[:X])(3), [4, 5, 12], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                uc[:curtis] = [1, 2, 3, 7, 8, 10, 4, 5, 9, 6, -12, -11]
            elseif [p, q, r] == [3, 3, 4]
                (uc[:families])[2] = Family(((CHEVIE[:families])[:X])(3), [2, 4, 23], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[6] = Family(((CHEVIE[:families])[:QZ])(3), [13, 9, 8, 10, 19, 22, 7, 21, 20], Dict{Symbol, Any}(:signs => [1, 1, 1, -1, -1, 1, -1, -1, 1], :special => 3, :cospecial => 2))
                (uc[:families])[9] = Family(ComplexConjugate(((CHEVIE[:families])[:X])(3)), [15, 14, 18], Dict{Symbol, Any}(:signs => [1, 1, -1]))
            elseif [p, q, r] == [3, 3, 5]
                (uc[:families])[3] = Family(((CHEVIE[:families])[:X])(3), [3, 6, 51], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[4] = Family(((CHEVIE[:families])[:X])(3), [4, 5, 54], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[6] = Family(((CHEVIE[:families])[:QZ])(3), [9, 10, 8, 21, 44, 46, 20, 49, 45], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, 1, 1, -1, -1]))
                (uc[:families])[7] = Family(((CHEVIE[:families])[:QZ])(3), [23, 11, 16, 12, 42, 50, 15, 48, 40], Dict{Symbol, Any}(:signs => [1, -1, -1, 1, 1, 1, 1, -1, -1], :special => 4, :cospecial => 7))
                (uc[:families])[8] = Family(ComplexConjugate(((CHEVIE[:families])[:X])(3)), [14, 13, 41], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[11] = Family(((CHEVIE[:families])[:X])(3), [19, 22, 47], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[13] = Family(((CHEVIE[:families])[:QZ])(3), [32, 27, 26, 28, 38, 53, 25, 52, 39], Dict{Symbol, Any}(:signs => [1, 1, 1, -1, -1, 1, -1, -1, 1], :special => 3, :cospecial => 2))
                (uc[:families])[15] = Family(ComplexConjugate(((CHEVIE[:families])[:X])(3)), [31, 30, 37], Dict{Symbol, Any}(:signs => [1, 1, -1]))
                (uc[:families])[16] = Family(ComplexConjugate(((CHEVIE[:families])[:X])(3)), [34, 33, 43], Dict{Symbol, Any}(:signs => [1, 1, -1]))
            elseif [p, q, r] == [4, 4, 3]
                (uc[:families])[2] = Family(((CHEVIE[:families])[:X])(4), [3, 2, 4, 14, 16, 13], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, -1, -1]))
                (uc[:families])[4] = Family(ComplexConjugate(((CHEVIE[:families])[:X])(4)), [8, 6, 7, 12, 15, 11], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, -1, -1]))
                uc[:curtis] = [1, 6, 7, 8, 10, 2, 3, 4, 9, 5, 14, 13, 12, 11, -16, -15]
            elseif [p, q, r] == [4, 4, 4]
                (uc[:families])[5] = Family(((CHEVIE[:families])[:X])(4), [5, 8, 9, 46, 53, 47], Dict{Symbol, Any}(:signs => [1, 1, 1, -1, -1, 1]))
                (uc[:families])[6] = Family("C2", [12, 7, 6, 42])
                (uc[:families])[7] = Family(((CHEVIE[:families])[:X])(4), [13, 10, 11, 41, 55, 43], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, -1], :special => 3, :cospecial => 1))
                (uc[:families])[9] = Family(((CHEVIE[:families])[:QZ])(4), [18, 21, 28, 22, 23, 49, 39, 54, 56, 40, 15, 36, 19, 52, 37, 51], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, -1], :special => 2, :cospecial => 4))
                (uc[:families])[10] = Family(ComplexConjugate(((CHEVIE[:families])[:X])(4)), [16, 17, 20, 38, 50, 34], Dict{Symbol, Any}(:signs => [1, 1, 1, -1, 1, 1], :special => 3, :cospecial => 1))
                (uc[:families])[12] = Family("C2", [27, 26, 25, 35])
                (uc[:families])[13] = Family(ComplexConjugate(((CHEVIE[:families])[:X])(4)), [30, 29, 31, 44, 48, 45], Dict{Symbol, Any}(:signs => [1, 1, 1, 1, 1, -1], :special => 3, :cospecial => 1))
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
