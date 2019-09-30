
chevieset(:B, :CartanMat, function (arg...,)
        local n, type_, a
        n = arg[1]
        if length(arg) == 2
            type_ = arg[2]
        else
            type_ = 2
        end
        a = (chevieget(:A, :CartanMat))(n)
        (a[1])[2] = -type_
        (a[2])[1] = 2 // (a[1])[2]
        return a
    end)
chevieset(:B, :PrintDiagram, function (r, indices, title, type_)
        local i
        print(title, " ")
        if type_ == 1
            print(indices[1], " >=> ", indices[2])
        elseif type_ == 2
            print(indices[1], " <=< ", indices[2])
        elseif type_ == ER(2)
            print(indices[1], " == ", indices[2])
        else
            print(indices[1], " ?==? ", indices[2])
        end
        for i = 3:r
            print(" - ", indices[i])
        end
        print("\n")
    end)
chevieset(:B, :ReflectionName, function (arg...,)
        local i, r, type_, option
        r = arg[1]
        option = arg[2]
        if length(arg) == 3
            type_ = arg[3]
        else
            type_ = 2
        end
        if type_ == 2
            if haskey(option, :TeX)
                if haskey(option, :Au)
                    return "D_8"
                else
                    return SPrint("B_", TeXBracket(r))
                end
            elseif haskey(option, :Au)
                return "D8"
            elseif haskey(option, :arg)
                return SPrint("\"B\",", r)
            else
                return SPrint("B", r)
            end
        elseif type_ == 1
            if haskey(option, :TeX)
                return SPrint("C_", TeXBracket(r))
            elseif haskey(option, :arg)
                return SPrint("\"C\",", r)
            else
                return SPrint("C", r)
            end
        elseif type_ == ER(2)
            if haskey(option, :TeX)
                return SPrint("B^{\\hbox{sym}}_", TeXBracket(r))
            elseif haskey(option, :arg)
                return SPrint("\"Bsym\",", r)
            else
                return SPrint("Bsym", r)
            end
        elseif haskey(option, :TeX)
            return SPrint("B^?_", TeXBracket(r), "(", Format(type_, option), ")")
        elseif haskey(option, :arg)
            return SPrint("\"B?\",", r, ",", type_)
        else
            return SPrint("B?", r, "(", Format(type_), ")")
        end
    end)
chevieset(:B, :GeneratingRoots, function (l, type_)
        local rts, i
        rts = map((i->begin
                        0 * (1:l)
                    end), 1:l)
        for i = 1:l - 1
            (rts[i])[[i, i + 1]] = [1, -1]
        end
        (rts[l])[l] = 2 // type_
        return rts[l:(l - 1) - l:1]
    end)
chevieset(:B, :ParabolicRepresentatives, function (l, s)
        return (chevieget(:imp, :ParabolicRepresentatives))(2, 1, l, s)
    end)
chevieset(:B, :ReflectionDegrees, (n->begin
            2 * (1:n)
        end))
chevieset(:B, :Size, function (arg...,)
        return 2 ^ arg[1] * factorial(arg[1])
    end)
chevieset(:B, :NrConjugacyClasses, function (arg...,)
        return NrPartitionTuples(arg[1], 2)
    end)
chevieset(:B, :WeightInfo, function (n, type_)
        if type_ == 2
            return Dict{Symbol, Any}(:minusculeWeights => [1], :minusculeCoweights => [n], :decompositions => [[1]], :moduli => [2])
        else
            return Dict{Symbol, Any}(:minusculeWeights => [n], :minusculeCoweights => [1], :decompositions => [[1]], :moduli => [2])
        end
    end)
chevieset(:B, :WordClass, function (pi,)
        local w, i, l, r
        w = []
        i = 1
        for l = reverse(pi[2])
            w = Append(w, i:(i - 1) - i:2)
            w = Append(w, 1:(i + l) - 1)
            i = i + l
        end
        for l = pi[1]
            r = mod(l, 2)
            w = Append(w, i + Concatenation(1:3 - 1:(l - 1) - r, 2:4 - 2:(l + r) - 2))
            i = i + l
        end
        return w
    end)
chevieset(:B, :ClassInfo, function (n,)
        local res
        res = (chevieget(:imp, :ClassInfo))(2, 1, n)
        res[:classtext] = map(chevieget(:B, :WordClass), res[:classparams])
        res[:classes] = map((x->begin
                        (res[:centralizers])[1] // x
                    end), res[:centralizers])
        return res
    end)
chevieset(:B, :ClassParameter, function (n, w)
        local x, i, res, mark, cyc, j
        x = Perm()
        for i = w
            if i == 1
                x = x * Perm(1, n + 1)
            else
                x = x * (Perm(i - 1, i))((i - 1) + n, i + n)
            end
        end
        res = [[], []]
        mark = 1:n
        for i = 1:n
            if mark[i] != 0
                cyc = CyclePermInt(x, i)
                if i + n in cyc
                    push!(res[2], length(cyc) // 2)
                else
                    push!(res[1], length(cyc))
                end
                for j = cyc
                    if j > n
                        mark[j - n] = 0
                    else
                        mark[j] = 0
                    end
                end
            end
        end
        Sort(res[1])
        Sort(res[2])
        return [reverse(res[1]), reverse(res[2])]
    end)
chevieset(:B, :CharParams, (n->begin
            PartitionTuples(n, 2)
        end))
chevieset(:B, :CharName, function (arg...,)
        return PartitionTupleToString(arg[2])
    end)
chevieset(:B, :LowestPowerFakeDegree, function (p,)
        local pp, m, res
        pp = SymbolPartitionTuple(p, 1)
        m = length(pp[2])
        res = pp[1] * (m:(m - 1) - m:0)
        if pp[2] != []
            res = res + pp[2] * (m - 1:(m - 2) - (m - 1):0)
        end
        return (2res + Sum(pp[2])) - (m * (m - 1) * (4m + 1)) // 6
    end)
chevieset(:B, :CharInfo, function (n,)
        local res
        res = Dict{Symbol, Any}(:charparams => (chevieget(:B, :CharParams))(n))
        res[:extRefl] = Concatenation(map((i->begin
                            Position(res[:charparams], [[n - i], fill(0, max(0, (1 + i) - 1)) + 1])
                        end), 0:n - 1), [Position(res[:charparams], [[], fill(0, max(0, (1 + n) - 1)) + 1])])
        res[:a] = map((p->begin
                        LowestPowerGenericDegreeSymbol(SymbolPartitionTuple(p, 1))
                    end), res[:charparams])
        res[:A] = map((p->begin
                        HighestPowerGenericDegreeSymbol(SymbolPartitionTuple(p, 1))
                    end), res[:charparams])
        res[:b] = map(chevieget(:B, :LowestPowerFakeDegree), res[:charparams])
        res[:B] = (res[:a] + res[:A]) - res[:b]
        return res
    end)
chevieset(:B, :PoincarePolynomial, function (n, para)
        local q1, q2
        q1 = -((para[1])[1]) // (para[1])[2]
        q2 = -((para[2])[1]) // (para[2])[2]
        return Product(0:n - 1, (i->begin
                        (q2 ^ i * q1 + 1) * Sum(0:i, (k->begin
                                        q2 ^ k
                                    end))
                    end))
    end)
chevieset(:B, :SchurElement, function (arg...,)
        return (chevieget(:imp, :SchurElement))(2, 1, arg[1], arg[2], arg[3], [])
    end)
chevieset(:B, :FactorizedSchurElement, function (arg...,)
        return (chevieget(:imp, :FactorizedSchurElement))(2, 1, arg[1], arg[2], arg[3], [])
    end)
chevieset(:B, :HeckeRepresentation, function (arg...,)
        return (chevieget(:imp, :HeckeRepresentation))(2, 1, arg[1], arg[2], [], arg[4])
    end)
chevieset(:B, :Representation, function (n, i)
        return (chevieget(:imp, :Representation))(2, 1, n, i)
    end)
chevieset(:B, :FakeDegree, function (n, c, q)
        return Value(CycPolFakeDegreeSymbol(SymbolPartitionTuple(c, 1)), q)
    end)
chevieset(:B, :DecompositionMatrix, function (l, p)
        local pp, dd, pt, decS
        decS = (i->begin
                    MatrixDecompositionMatrix(DecompositionMatrix(Specht(p, p), i))
                end)
        pp = map(Partitions, 0:l)
        pt = PartitionTuples(l, 2)
        if p == 2
            return [[1:length(pt), map(function (p,)
                                    p = LittlewoodRichardsonRule(p[1], p[2])
                                    return map(function (x,)
                                                if x in p
                                                    return 1
                                                else
                                                    return 0
                                                end
                                            end, pp[l + 1])
                                end, pt) * decS(l)]]
        else
            dd = Concatenation([[[1]], [[1]]], map(decS, 2:l))
            return map((i->begin
                            [map((x->begin
                                            Position(pt, x)
                                        end), Cartesian(pp[i + 1], pp[(l + 1) - i])), map((x->begin
                                            map(Product, Cartesian(x))
                                        end), Cartesian(dd[i + 1], dd[(l + 1) - i]))]
                        end), 0:l)
        end
    end)
chevieset(:B, :UnipotentCharacters, function (arg...,)
        local uc, symbols, r, d, s, rank
        rank = arg[1]
        uc = Dict{Symbol, Any}(:harishChandra => [], :charSymbols => [])
        for d = 1 + 2 * (0:div(-1 + RootInt(1 + 4rank, 2), 2))
            s = div(d ^ 2 - 1, 4)
            s = Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "B", :indices => 1 + s:rank, :rank => rank - s), :levi => 1:s, :eigenvalue => (-1) ^ div(d + 1, 4), :parameterExponents => Concatenation([d], fill(0, max(0, (1 + rank) - (2 + s))) + 1), :cuspidalName => SPrint("B_{", s, "}"))
            push!(uc[:harishChandra], s)
            symbols = BDSymbols(rank, d)
            s[:charNumbers] = (1:length(symbols)) + length(uc[:charSymbols])
            FixRelativeType(s)
            uc[:charSymbols] = Append(uc[:charSymbols], symbols)
        end
        ((uc[:harishChandra])[1])[:cuspidalName] = ""
        uc[:a] = map(LowestPowerGenericDegreeSymbol, uc[:charSymbols])
        uc[:A] = map(HighestPowerGenericDegreeSymbol, uc[:charSymbols])
        uc[:families] = FamiliesClassical(uc[:charSymbols])
        if length(arg) == 2 && arg[2] == 1
            (((uc[:harishChandra])[1])[:relativeType])[:cartanType] = 1
        end
        return uc
    end)
chevieset(:B, :UnipotentClasses, function (r, char, type_)
        local cl, uc, i, l, s, cc, ss, symbol2para, part2dynkin, addSpringer, d, LuSpin, trspringer, j
        part2dynkin = function (part,)
                local p, res
                p = Concatenation(map((d->begin
                                    1 - d:(3 - d) - (1 - d):d - 1
                                end), part))
                Sort(p)
                p = p[div(3 + length(p), 2):length(p)]
                if type_ == 1
                    res = [2 * p[1]]
                else
                    res = [p[1]]
                end
                res = Append(res, p[2:length(p)] - p[1:length(p) - 1])
                return res
            end
        addSpringer = function (s,)
                local ss, p
                ss = First(uc[:springerSeries], (x->begin
                                x[:defect] == DefectSymbol(s[:symbol])
                            end))
                if s[:sp] == [[], []]
                    p = 1
                elseif s[:sp] == [[1], []]
                    p = 2
                elseif s[:sp] == [[], [1]]
                    p = 1
                else
                    p = Position(CharParams(ss[:relgroup]), [s[:sp]])
                end
                (ss[:locsys])[p] = [length(uc[:classes]), Position(CharParams(cc[:Au]), map(function (x,)
                                    if x
                                        return [1, 1]
                                    else
                                        return [2]
                                    end
                                end, s[:Au]))]
            end
        if type_ == ER(2)
            type_ = 2
            char = 2
        end
        if char == 2
            ss = XSP(4, 2, r)
        elseif type_ == 1
            ss = XSP(2, 1, r)
        else
            ss = XSP(2, 0, r)
        end
        l = Union(map((c->begin
                            map((x->begin
                                        [DefectSymbol(x[:symbol]), Sum(x[:sp], Sum)]
                                    end), c)
                        end), ss))
        SortBy(l, (x->begin
                    [AbsInt(x[1]), -(SignInt(x[1]))]
                end))
        uc = Dict{Symbol, Any}(:classes => [], :springerSeries => map(function (d,)
                            local res
                            res = Dict{Symbol, Any}(:relgroup => CoxeterGroup("C", d[2]), :defect => d[1], :levi => 1:r - d[2])
                            res[:locsys] = map((x->begin
                                            [0, 0]
                                        end), 1:NrConjugacyClasses(res[:relgroup]))
                            if char == 2
                                res[:Z] = [1]
                            elseif type_ == 1
                                res[:Z] = [(-1) ^ (r - d[2])]
                            elseif IsInt(ER(2 * (r - d[2]) + 1))
                                res[:Z] = [1]
                            else
                                res[:Z] = [-1]
                            end
                            return res
                        end, l))
        if char != 2
            symbol2para = function (S,)
                    local c, i, l, part, d
                    c = Concatenation(S)
                    Sort(c)
                    i = 1
                    part = []
                    d = mod(type_, 2)
                    while i <= length(c)
                        if i == length(c) || c[i + 1] - c[i] > 0
                            push!(part, (2 * (c[i] - (i - 1)) + 1) - d)
                            i = i + 1
                        else
                            l = 2 * (c[i] - (i - 1)) - d
                            part = Append(part, [l, l])
                            i = i + 2
                        end
                    end
                    Sort(part)
                    part = Filtered(part, (y->begin
                                    y != 0
                                end))
                    return reverse(part)
                end
        else
            symbol2para = function (S,)
                    local c, i, l, part, ex
                    c = Concatenation(S)
                    Sort(c)
                    i = 1
                    part = []
                    ex = []
                    while i <= length(c)
                        if i == length(c) || c[i + 1] - c[i] > 1
                            push!(part, 2 * (c[i] - 2 * (i - 1)))
                            i = i + 1
                        elseif c[i] == c[i + 1]
                            l = 2 * (c[i] - 2 * (i - 1)) - 2
                            part = Append(part, [l, l])
                            push!(ex, l)
                            i = i + 2
                        elseif c[i] + 1 == c[i + 1]
                            l = 2 * (c[i] - 2 * (i - 1)) - 1
                            part = Append(part, [l, l])
                            i = i + 2
                        end
                    end
                    Sort(part)
                    part = Filtered(part, (y->begin
                                    y != 0
                                end))
                    return [reverse(part), ex]
                end
        end
        if char == 2
            type_ = 1
        end
        for cl = ss
            cc = Dict{Symbol, Any}(:parameter => symbol2para((cl[1])[:symbol]))
            cc[:Au] = ApplyFunc(CoxeterGroup, Concatenation(map((x->begin
                                    ["A", 1]
                                end), (cl[1])[:Au])))
            if char != 2
                cc[:dynkin] = part2dynkin(cc[:parameter])
                cc[:name] = IntListToString(cc[:parameter])
            else
                type_ = 1
                cc[:dimBu] = (cl[1])[:dimBu]
                cc[:name] = Join(map(function (x,)
                                local res
                                res = IntListToString(fill(0, max(0, (1 + x[2]) - 1)) + x[1], "[]")
                                if x[1] in (cc[:parameter])[2]
                                    return SPrint("(", res, ")")
                                end
                                return res
                            end, reverse(Collected((cc[:parameter])[1]))), "")
            end
            cc[:red] = CoxeterGroup()
            if char == 2
                j = (cc[:parameter])[1]
            else
                j = cc[:parameter]
            end
            for j = Collected(j)
                if mod(j[1], 2) == mod(type_, 2)
                    cc[:red] = cc[:red] * CoxeterGroup("C", j[2] // 2)
                elseif mod(j[2], 2) != 0
                    if j[2] > 1
                        cc[:red] = cc[:red] * CoxeterGroup("B", (j[2] - 1) // 2)
                    end
                elseif j[2] > 2
                    cc[:red] = cc[:red] * CoxeterGroup("D", j[2] // 2)
                else
                    cc[:red] = cc[:red] * Torus(1)
                end
            end
            push!(uc[:classes], cc)
            for s = cl
                addSpringer(s)
            end
        end
        uc[:orderClasses] = Hasse(Poset(map((x->begin
                                map(function (y,)
                                        local m, f, fx, fy, i
                                        if char != 2
                                            return Dominates(y[:parameter], x[:parameter])
                                        end
                                        m = maximum(((x[:parameter])[1])[1], ((y[:parameter])[1])[1])
                                        f = (x->begin
                                                    map((i->begin
                                                                Sum(Filtered(x, (z->begin
                                                                                    z < i
                                                                                end))) + i * count((z->begin
                                                                                    z >= i
                                                                                end), x)
                                                            end), 1:m)
                                                end)
                                        fx = f((x[:parameter])[1])
                                        fy = f((y[:parameter])[1])
                                        for i = 1:m
                                            if fx[i] < fy[i]
                                                return false
                                            elseif fx[i] == fy[i] && i in (y[:parameter])[2]
                                                if i in Difference((x[:parameter])[1], (x[:parameter])[2])
                                                    return false
                                                end
                                                if i < m && mod(fx[i + 1] - fy[i + 1], 2) == 1
                                                    return false
                                                end
                                            end
                                        end
                                        return true
                                    end, uc[:classes])
                            end), uc[:classes])))
        if char != 2 && type_ == 2
            LuSpin = function (p,)
                    local t, a, b, i, j, l, d
                    Sort(p)
                    a = []
                    b = []
                    d = [0, 1, 0, -1]
                    d = d[map((x->begin
                                        1 + mod(x, 4)
                                    end), p)]
                    i = 1
                    while i <= length(p)
                        l = p[i]
                        t = Sum(d[1:i - 1])
                        if 1 == mod(l, 4)
                            push!(a, div(l - 1, 4) - t)
                            i = i + 1
                        elseif 3 == mod(l, 4)
                            push!(b, div(l - 3, 4) + t)
                            i = i + 1
                        else
                            j = i
                            while i <= length(p) && p[i] == l
                                i = i + 1
                            end
                            j = fill(0, max(0, (1 + div(i - j, 2)) - 1))
                            a = Append(a, (j + div(l + mod(l, 4), 4)) - t)
                            b = Append(b, j + div(l - mod(l, 4), 4) + t)
                        end
                    end
                    a = Filtered(a, (x->begin
                                    x != 0
                                end))
                    a = reverse(a)
                    b = Filtered(b, (x->begin
                                    x != 0
                                end))
                    b = reverse(b)
                    if Sum(d) >= 1
                        return [a, b]
                    else
                        return [b, a]
                    end
                end
            addSpringer = function (f, i, s, k)
                    local ss, p
                    ss = First(uc[:springerSeries], f)
                    if s in [[[], [1]], [[], []]]
                        p = 1
                    elseif s == [[1], []]
                        p = 2
                    else
                        p = Position(CharParams(ss[:relgroup]), [s])
                    end
                    (ss[:locsys])[p] = [i, k]
                end
            trspringer = function (i, old, new)
                    local ss, c, p
                    for ss = uc[:springerSeries]
                        for c = ss[:locsys]
                            if c[1] == i
                                p = Position(old, c[2])
                                if p != false
                                    c[2] = new[p]
                                end
                            end
                        end
                    end
                end
            d = 0
            while 4 * d ^ 2 - 3d <= r
                i = 4 * d ^ 2 - 3d
                if mod(r - d, 2) == 0
                    l = Concatenation(1:i, i + 2:(i + 4) - (i + 2):r)
                    push!(uc[:springerSeries], Dict{Symbol, Any}(:relgroup => CoxeterGroup("B", (r - i) // 2), :levi => l, :Z => [-1], :locsys => []))
                    i = 4 * d ^ 2 + 3d
                    if i <= r && d != 0
                        l = Concatenation(1:i, i + 2:(i + 4) - (i + 2):r)
                        push!(uc[:springerSeries], Dict{Symbol, Any}(:relgroup => CoxeterGroup("B", (r - i) // 2), :levi => l, :Z => [-1], :locsys => []))
                    end
                end
                d = d + 1
            end
            l = Filtered(1:length(uc[:classes]), (i->begin
                            ForAll(Collected(((uc[:classes])[i])[:parameter]), (c->begin
                                        mod(c[1], 2) == 0 || c[2] == 1
                                    end))
                        end))
            for i = l
                cl = (uc[:classes])[i]
                s = LuSpin(cl[:parameter])
                if Size(cl[:Au]) == 1
                    cl[:Au] = CoxeterGroup("A", 1)
                    trspringer(i, [1], [2])
                    d = 1
                elseif Size(cl[:Au]) == 4
                    cl[:Au] = CoxeterGroup("B", 2)
                    trspringer(i, [1, 2, 3, 4], [1, 3, 5, 4])
                    d = 2
                else
                    error("Au non-commutative of order ", Size(cl[:Au]) * 2, "  !  implemented")
                end
                addSpringer((ss->begin
                            ss[:Z] == [-1] && Rank(ss[:relgroup]) == Sum(s, Sum)
                        end), i, s, d)
            end
        end
        return uc
    end)
chevieset(:B, :Invariants, function (n, type_)
        local m
        m = fill(0, max(0, (1 + n) - 1)) + 1
        m[1] = 2 // type_
        m = DiagonalMat(m) * (chevieget(:imp, :GeneratingRoots))(2, 1, n)
        return map((f->begin
                        function (arg...,)
                            return ApplyFunc(f, arg * m)
                        end
                    end), ((CHEVIE[:imp])[:Invariants])(2, 1, n))
    end)