
chevieset(:B, :CartanMat, function (arg...,)
        #= none:17 =#
        local n, type_, a
        #= none:19 =#
        n = arg[1]
        #= none:21 =#
        if length(arg) == 2
            #= none:21 =#
            type_ = arg[2]
        else
            #= none:22 =#
            type_ = 2
        end
        #= none:24 =#
        a = (chevieget(:A, :CartanMat))(n)
        #= none:26 =#
        (a[1])[2] = -type_
        #= none:27 =#
        (a[2])[1] = 2 // (a[1])[2]
        #= none:29 =#
        return a
    end)
chevieset(:B, :PrintDiagram, function (r, indices, title, type_)
        #= none:3 =#
        local i
        #= none:5 =#
        print(title, " ")
        #= none:7 =#
        if type_ == 1
            #= none:7 =#
            print(indices[1], " >=> ", indices[2])
        elseif #= none:9 =# type_ == 2
            #= none:9 =#
            print(indices[1], " <=< ", indices[2])
        elseif #= none:11 =# type_ == ER(2)
            #= none:11 =#
            print(indices[1], " == ", indices[2])
        else
            #= none:13 =#
            print(indices[1], " ?==? ", indices[2])
        end
        #= none:15 =#
        for i = 3:r
            #= none:15 =#
            print(" - ", indices[i])
        end
        #= none:17 =#
        print("\n")
    end)
chevieset(:B, :ReflectionName, function (arg...,)
        #= none:3 =#
        local i, r, type_, option
        #= none:5 =#
        r = arg[1]
        #= none:6 =#
        option = arg[2]
        #= none:8 =#
        if length(arg) == 3
            #= none:8 =#
            type_ = arg[3]
        else
            #= none:9 =#
            type_ = 2
        end
        #= none:11 =#
        if type_ == 2
            #= none:12 =#
            if haskey(option, :TeX)
                #= none:13 =#
                if haskey(option, :Au)
                    #= none:13 =#
                    return "D_8"
                else
                    #= none:15 =#
                    return SPrint("B_", TeXBracket(r))
                end
            elseif #= none:17 =# haskey(option, :Au)
                #= none:17 =#
                return "D8"
            elseif #= none:19 =# haskey(option, :arg)
                #= none:19 =#
                return SPrint("\"B\",", r)
            else
                #= none:21 =#
                return SPrint("B", r)
            end
        elseif #= none:23 =# type_ == 1
            #= none:24 =#
            if haskey(option, :TeX)
                #= none:24 =#
                return SPrint("C_", TeXBracket(r))
            elseif #= none:26 =# haskey(option, :arg)
                #= none:26 =#
                return SPrint("\"C\",", r)
            else
                #= none:28 =#
                return SPrint("C", r)
            end
        elseif #= none:30 =# type_ == ER(2)
            #= none:31 =#
            if haskey(option, :TeX)
                #= none:31 =#
                return SPrint("B^{\\hbox{sym}}_", TeXBracket(r))
            elseif #= none:33 =# haskey(option, :arg)
                #= none:33 =#
                return SPrint("\"Bsym\",", r)
            else
                #= none:35 =#
                return SPrint("Bsym", r)
            end
        elseif #= none:37 =# haskey(option, :TeX)
            #= none:38 =#
            return SPrint("B^?_", TeXBracket(r), "(", Format(type_, option), ")")
        elseif #= none:40 =# haskey(option, :arg)
            #= none:40 =#
            return SPrint("\"B?\",", r, ",", type_)
        else
            #= none:42 =#
            return SPrint("B?", r, "(", Format(type_), ")")
        end
    end)
chevieset(:B, :GeneratingRoots, function (l, type_)
        #= none:3 =#
        local rts, i
        #= none:5 =#
        rts = map((i->begin
                        #= none:5 =#
                        0 * (1:l)
                    end), 1:l)
        #= none:7 =#
        for i = 1:l - 1
            #= none:7 =#
            (rts[i])[[i, i + 1]] = [1, -1]
        end
        #= none:9 =#
        (rts[l])[l] = 2 // type_
        #= none:11 =#
        return rts[l:(l - 1) - l:1]
    end)
chevieset(:B, :ParabolicRepresentatives, function (l, s)
        #= none:4 =#
        return (chevieget(:imp, :ParabolicRepresentatives))(2, 1, l, s)
    end)
chevieset(:B, :ReflectionDegrees, (n->begin
            #= none:3 =#
            2 * (1:n)
        end))
chevieset(:B, :Size, function (arg...,)
        #= none:3 =#
        return 2 ^ arg[1] * factorial(arg[1])
    end)
chevieset(:B, :NrConjugacyClasses, function (arg...,)
        #= none:4 =#
        return NrPartitionTuples(arg[1], 2)
    end)
chevieset(:B, :WeightInfo, function (n, type_)
        #= none:4 =#
        if type_ == 2
            #= none:4 =#
            return Dict{Symbol, Any}(:minusculeWeights => [1], :minusculeCoweights => [n], :decompositions => [[1]], :moduli => [2])
        else
            #= none:7 =#
            return Dict{Symbol, Any}(:minusculeWeights => [n], :minusculeCoweights => [1], :decompositions => [[1]], :moduli => [2])
        end
    end)
chevieset(:B, :WordClass, function (pi,)
        #= none:8 =#
        local w, i, l, r
        #= none:10 =#
        w = []
        #= none:11 =#
        i = 1
        #= none:13 =#
        for l = reverse(pi[2])
            #= none:14 =#
            w = Append(w, i:(i - 1) - i:2)
            #= none:15 =#
            w = Append(w, 1:(i + l) - 1)
            #= none:16 =#
            i = i + l
        end
        #= none:18 =#
        for l = pi[1]
            #= none:19 =#
            r = mod(l, 2)
            #= none:21 =#
            w = Append(w, i + Concatenation(1:3 - 1:(l - 1) - r, 2:4 - 2:(l + r) - 2))
            #= none:23 =#
            i = i + l
        end
        #= none:25 =#
        return w
    end)
chevieset(:B, :ClassInfo, function (n,)
        #= none:18 =#
        local res
        #= none:20 =#
        res = (chevieget(:imp, :ClassInfo))(2, 1, n)
        #= none:22 =#
        res[:classtext] = map(chevieget(:B, :WordClass), res[:classparams])
        #= none:24 =#
        res[:classes] = map((x->begin
                        #= none:24 =#
                        (res[:centralizers])[1] // x
                    end), res[:centralizers])
        #= none:26 =#
        return res
    end)
chevieset(:B, :ClassParameter, function (n, w)
        #= none:9 =#
        local x, i, res, mark, cyc, j
        #= none:11 =#
        x = Perm()
        #= none:13 =#
        for i = w
            #= none:14 =#
            if i == 1
                #= none:14 =#
                x = x * Perm(1, n + 1)
            else
                #= none:15 =#
                x = x * (Perm(i - 1, i))((i - 1) + n, i + n)
            end
        end
        #= none:17 =#
        res = [[], []]
        #= none:19 =#
        mark = 1:n
        #= none:21 =#
        for i = 1:n
            #= none:22 =#
            if mark[i] != 0
                #= none:23 =#
                cyc = CyclePermInt(x, i)
                #= none:25 =#
                if i + n in cyc
                    #= none:25 =#
                    push!(res[2], length(cyc) // 2)
                else
                    #= none:27 =#
                    push!(res[1], length(cyc))
                end
                #= none:29 =#
                for j = cyc
                    #= none:30 =#
                    if j > n
                        #= none:30 =#
                        mark[j - n] = 0
                    else
                        #= none:31 =#
                        mark[j] = 0
                    end
                end
            end
        end
        #= none:35 =#
        sort!(res[1])
        #= none:36 =#
        sort!(res[2])
        #= none:38 =#
        return [reverse(res[1]), reverse(res[2])]
    end)
chevieset(:B, :CharParams, (n->begin
            #= none:7 =#
            partition_tuples(n, 2)
        end))
chevieset(:B, :CharName, function (arg...,)
        #= none:5 =#
        return PartitionTupleToString(arg[2])
    end)
chevieset(:B, :LowestPowerFakeDegree, function (p,)
        #= none:3 =#
        local pp, m, res
        #= none:5 =#
        pp = SymbolPartitionTuple(p, 1)
        #= none:6 =#
        m = length(pp[2])
        #= none:8 =#
        res = pp[1] * (m:(m - 1) - m:0)
        #= none:10 =#
        if pp[2] != []
            #= none:10 =#
            res = res + pp[2] * (m - 1:(m - 2) - (m - 1):0)
        end
        #= none:12 =#
        return (2res + Sum(pp[2])) - (m * (m - 1) * (4m + 1)) // 6
    end)
chevieset(:B, :CharInfo, function (n,)
        #= none:3 =#
        local res
        #= none:5 =#
        res = Dict{Symbol, Any}(:charparams => (chevieget(:B, :CharParams))(n))
        #= none:7 =#
        res[:extRefl] = Concatenation(map((i->begin
                            #= none:8 =#
                            Position(res[:charparams], [[n - i], fill(0, max(0, (1 + i) - 1)) + 1])
                        end), 0:n - 1), [Position(res[:charparams], [[], fill(0, max(0, (1 + n) - 1)) + 1])])
        #= none:11 =#
        res[:a] = map((p->begin
                        #= none:12 =#
                        LowestPowerGenericDegreeSymbol(SymbolPartitionTuple(p, 1))
                    end), res[:charparams])
        #= none:14 =#
        res[:A] = map((p->begin
                        #= none:15 =#
                        HighestPowerGenericDegreeSymbol(SymbolPartitionTuple(p, 1))
                    end), res[:charparams])
        #= none:17 =#
        res[:b] = map(chevieget(:B, :LowestPowerFakeDegree), res[:charparams])
        #= none:19 =#
        res[:B] = (res[:a] + res[:A]) - res[:b]
        #= none:21 =#
        return res
    end)
chevieset(:B, :PoincarePolynomial, function (n, para)
        #= none:11 =#
        local q1, q2
        #= none:13 =#
        q1 = -((para[1])[1]) // (para[1])[2]
        #= none:14 =#
        q2 = -((para[2])[1]) // (para[2])[2]
        #= none:16 =#
        return Product(0:n - 1, (i->begin
                        #= none:16 =#
                        (q2 ^ i * q1 + 1) * Sum(0:i, (k->begin
                                        #= none:16 =#
                                        q2 ^ k
                                    end))
                    end))
    end)
chevieset(:B, :SchurElement, function (arg...,)
        #= none:12 =#
        return (chevieget(:imp, :SchurElement))(2, 1, arg[1], arg[2], arg[3], [])
    end)
chevieset(:B, :FactorizedSchurElement, function (arg...,)
        #= none:4 =#
        return (chevieget(:imp, :FactorizedSchurElement))(2, 1, arg[1], arg[2], arg[3], [])
    end)
chevieset(:B, :HeckeRepresentation, function (arg...,)
        #= none:4 =#
        return (chevieget(:imp, :HeckeRepresentation))(2, 1, arg[1], arg[2], [], arg[4])
    end)
chevieset(:B, :Representation, function (n, i)
        #= none:4 =#
        return (chevieget(:imp, :Representation))(2, 1, n, i)
    end)
chevieset(:B, :FakeDegree, function (n, c, q)
        #= none:7 =#
        return Value(CycPolFakeDegreeSymbol(SymbolPartitionTuple(c, 1)), q)
    end)
chevieset(:B, :DecompositionMatrix, function (l, p)
        #= none:3 =#
        local pp, dd, pt, decS
        #= none:5 =#
        decS = (i->begin
                    #= none:5 =#
                    MatrixDecompositionMatrix(DecompositionMatrix(Specht(p, p), i))
                end)
        #= none:7 =#
        pp = map(partitions, 0:l)
        #= none:8 =#
        pt = partition_tuples(l, 2)
        #= none:10 =#
        if p == 2
            #= none:11 =#
            return [[1:length(pt), map(function (p,)
                                    #= none:12 =#
                                    p = LittlewoodRichardsonRule(p[1], p[2])
                                    #= none:14 =#
                                    return map(function (x,)
                                                #= none:15 =#
                                                if x in p
                                                    #= none:15 =#
                                                    return 1
                                                else
                                                    #= none:16 =#
                                                    return 0
                                                end
                                            end, pp[l + 1])
                                end, pt) * decS(l)]]
        else
            #= none:20 =#
            dd = Concatenation([[[1]], [[1]]], map(decS, 2:l))
            #= none:22 =#
            return map((i->begin
                            #= none:22 =#
                            [map((x->begin
                                            #= none:23 =#
                                            Position(pt, x)
                                        end), Cartesian(pp[i + 1], pp[(l + 1) - i])), map((x->begin
                                            #= none:24 =#
                                            map(Product, Cartesian(x))
                                        end), Cartesian(dd[i + 1], dd[(l + 1) - i]))]
                        end), 0:l)
        end
    end)
chevieset(:B, :UnipotentCharacters, function (arg...,)
        #= none:5 =#
        local uc, symbols, r, d, s, rank
        #= none:7 =#
        rank = arg[1]
        #= none:9 =#
        uc = Dict{Symbol, Any}(:harishChandra => [], :charSymbols => [])
        #= none:11 =#
        for d = 1 + 2 * (0:div(-1 + RootInt(1 + 4rank, 2), 2))
            #= none:12 =#
            r = div(d ^ 2 - 1, 4)
            #= none:14 =#
            s = Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "B", :indices => 1 + r:rank, :rank => rank - r), :levi => 1:r, :eigenvalue => (-1) ^ div(d + 1, 4), :parameterExponents => Concatenation([d], fill(0, max(0, (1 + rank) - (2 + r))) + 1), :cuspidalName => SPrint("B_{", r, "}"))
            #= none:20 =#
            if r < 10
                #= none:20 =#
                s[:cuspidalName] = SPrint("B_", r)
            end
            #= none:22 =#
            push!(uc[:harishChandra], s)
            #= none:24 =#
            symbols = BDSymbols(rank, d)
            #= none:26 =#
            s[:charNumbers] = (1:length(symbols)) + length(uc[:charSymbols])
            #= none:28 =#
            FixRelativeType(s)
            #= none:30 =#
            uc[:charSymbols] = Append(uc[:charSymbols], symbols)
        end
        #= none:32 =#
        ((uc[:harishChandra])[1])[:cuspidalName] = ""
        #= none:34 =#
        uc[:a] = map(LowestPowerGenericDegreeSymbol, uc[:charSymbols])
        #= none:36 =#
        uc[:A] = map(HighestPowerGenericDegreeSymbol, uc[:charSymbols])
        #= none:38 =#
        uc[:families] = FamiliesClassical(uc[:charSymbols])
        #= none:40 =#
        if length(arg) == 2 && arg[2] == 1
            #= none:41 =#
            (((uc[:harishChandra])[1])[:relativeType])[:cartanType] = 1
        end
        #= none:43 =#
        return uc
    end)
chevieset(:B, :UnipotentClasses, function (r, char, type_)
        #= none:16 =#
        local cl, uc, i, l, s, cc, ss, symbol2para, part2dynkin, addSpringer, d, LuSpin, trspringer, j
        #= none:19 =#
        part2dynkin = function (part,)
                #= none:19 =#
                local p, res
                #= none:21 =#
                p = Concatenation(map((d->begin
                                    #= none:21 =#
                                    1 - d:(3 - d) - (1 - d):d - 1
                                end), part))
                #= none:23 =#
                sort!(p)
                #= none:24 =#
                p = p[div(3 + length(p), 2):length(p)]
                #= none:26 =#
                if type_ == 1
                    #= none:26 =#
                    res = [2 * p[1]]
                else
                    #= none:27 =#
                    res = [p[1]]
                end
                #= none:29 =#
                res = Append(res, p[2:length(p)] - p[1:length(p) - 1])
                #= none:31 =#
                return res
            end
        #= none:35 =#
        addSpringer = function (s,)
                #= none:35 =#
                local ss, p
                #= none:37 =#
                ss = First(uc[:springerSeries], (x->begin
                                #= none:37 =#
                                x[:defect] == DefectSymbol(s[:symbol])
                            end))
                #= none:39 =#
                if s[:sp] == [[], []]
                    #= none:39 =#
                    p = 1
                elseif #= none:41 =# s[:sp] == [[1], []]
                    #= none:41 =#
                    p = 2
                elseif #= none:43 =# s[:sp] == [[], [1]]
                    #= none:43 =#
                    p = 1
                else
                    #= none:45 =#
                    p = Position(CharParams(ss[:relgroup]), [s[:sp]])
                end
                #= none:47 =#
                (ss[:locsys])[p] = [length(uc[:classes]), Position(CharParams(cc[:Au]), map(function (x,)
                                    #= none:48 =#
                                    if x
                                        #= none:48 =#
                                        return [1, 1]
                                    else
                                        #= none:49 =#
                                        return [2]
                                    end
                                end, s[:Au]))]
            end
        #= none:54 =#
        if type_ == ER(2)
            #= none:54 =#
            type_ = 2
            #= none:55 =#
            char = 2
        end
        #= none:57 =#
        if char == 2
            #= none:57 =#
            ss = XSP(4, 2, r)
        elseif #= none:59 =# type_ == 1
            #= none:59 =#
            ss = XSP(2, 1, r)
        else
            #= none:61 =#
            ss = XSP(2, 0, r)
        end
        #= none:63 =#
        l = Union(map((c->begin
                            #= none:63 =#
                            map((x->begin
                                        #= none:63 =#
                                        [DefectSymbol(x[:symbol]), Sum(x[:sp], Sum)]
                                    end), c)
                        end), ss))
        #= none:65 =#
        SortBy(l, (x->begin
                    #= none:65 =#
                    [abs(x[1]), -(sign(x[1]))]
                end))
        #= none:67 =#
        uc = Dict{Symbol, Any}(:classes => [], :springerSeries => map(function (d,)
                            #= none:67 =#
                            local res
                            #= none:69 =#
                            res = Dict{Symbol, Any}(:relgroup => CoxeterGroup("C", d[2]), :defect => d[1], :levi => 1:r - d[2])
                            #= none:71 =#
                            res[:locsys] = map((x->begin
                                            #= none:71 =#
                                            [0, 0]
                                        end), 1:NrConjugacyClasses(res[:relgroup]))
                            #= none:73 =#
                            if char == 2
                                #= none:73 =#
                                res[:Z] = [1]
                            elseif #= none:75 =# type_ == 1
                                #= none:75 =#
                                res[:Z] = [(-1) ^ (r - d[2])]
                            elseif #= none:77 =# IsInt(ER(2 * (r - d[2]) + 1))
                                #= none:77 =#
                                res[:Z] = [1]
                            else
                                #= none:78 =#
                                res[:Z] = [-1]
                            end
                            #= none:80 =#
                            return res
                        end, l))
        #= none:83 =#
        if char != 2
            #= none:84 =#
            symbol2para = function (S,)
                    #= none:84 =#
                    local c, i, l, part, d
                    #= none:86 =#
                    c = Concatenation(S)
                    #= none:87 =#
                    sort!(c)
                    #= none:88 =#
                    i = 1
                    #= none:89 =#
                    part = []
                    #= none:91 =#
                    d = mod(type_, 2)
                    #= none:93 =#
                    while i <= length(c)
                        #= none:94 =#
                        if i == length(c) || c[i + 1] - c[i] > 0
                            #= none:94 =#
                            push!(part, (2 * (c[i] - (i - 1)) + 1) - d)
                            #= none:95 =#
                            i = i + 1
                        else
                            #= none:97 =#
                            l = 2 * (c[i] - (i - 1)) - d
                            #= none:98 =#
                            part = Append(part, [l, l])
                            #= none:99 =#
                            i = i + 2
                        end
                    end
                    #= none:101 =#
                    sort!(part)
                    #= none:102 =#
                    part = Filtered(part, (y->begin
                                    #= none:102 =#
                                    y != 0
                                end))
                    #= none:104 =#
                    return reverse(part)
                end
        else
            #= none:109 =#
            symbol2para = function (S,)
                    #= none:109 =#
                    local c, i, l, part, ex
                    #= none:111 =#
                    c = Concatenation(S)
                    #= none:112 =#
                    sort!(c)
                    #= none:113 =#
                    i = 1
                    #= none:114 =#
                    part = []
                    #= none:115 =#
                    ex = []
                    #= none:117 =#
                    while i <= length(c)
                        #= none:118 =#
                        if i == length(c) || c[i + 1] - c[i] > 1
                            #= none:118 =#
                            push!(part, 2 * (c[i] - 2 * (i - 1)))
                            #= none:119 =#
                            i = i + 1
                        elseif #= none:121 =# c[i] == c[i + 1]
                            #= none:122 =#
                            l = 2 * (c[i] - 2 * (i - 1)) - 2
                            #= none:123 =#
                            part = Append(part, [l, l])
                            #= none:124 =#
                            push!(ex, l)
                            #= none:125 =#
                            i = i + 2
                        elseif #= none:127 =# c[i] + 1 == c[i + 1]
                            #= none:128 =#
                            l = 2 * (c[i] - 2 * (i - 1)) - 1
                            #= none:129 =#
                            part = Append(part, [l, l])
                            #= none:130 =#
                            i = i + 2
                        end
                    end
                    #= none:132 =#
                    sort!(part)
                    #= none:133 =#
                    part = Filtered(part, (y->begin
                                    #= none:133 =#
                                    y != 0
                                end))
                    #= none:135 =#
                    return [reverse(part), ex]
                end
        end
        #= none:139 =#
        if char == 2
            #= none:139 =#
            type_ = 1
        end
        #= none:141 =#
        for cl = ss
            #= none:142 =#
            cc = Dict{Symbol, Any}(:parameter => symbol2para((cl[1])[:symbol]))
            #= none:144 =#
            cc[:Au] = ApplyFunc(CoxeterGroup, Concatenation(map((x->begin
                                    #= none:144 =#
                                    ["A", 1]
                                end), (cl[1])[:Au])))
            #= none:146 =#
            if char != 2
                #= none:147 =#
                cc[:dynkin] = part2dynkin(cc[:parameter])
                #= none:149 =#
                cc[:name] = IntListToString(cc[:parameter])
            else
                #= none:151 =#
                type_ = 1
                #= none:153 =#
                cc[:dimBu] = (cl[1])[:dimBu]
                #= none:155 =#
                cc[:name] = Join(map(function (x,)
                                #= none:156 =#
                                local res
                                #= none:157 =#
                                res = IntListToString(fill(0, max(0, (1 + x[2]) - 1)) + x[1], "[]")
                                #= none:159 =#
                                if x[1] in (cc[:parameter])[2]
                                    #= none:159 =#
                                    return SPrint("(", res, ")")
                                end
                                #= none:161 =#
                                return res
                            end, reverse(Collected((cc[:parameter])[1]))), "")
            end
            #= none:164 =#
            cc[:red] = CoxeterGroup()
            #= none:166 =#
            if char == 2
                #= none:166 =#
                j = (cc[:parameter])[1]
            else
                #= none:167 =#
                j = cc[:parameter]
            end
            #= none:169 =#
            for j = Collected(j)
                #= none:170 =#
                if mod(j[1], 2) == mod(type_, 2)
                    #= none:170 =#
                    cc[:red] = cc[:red] * CoxeterGroup("C", j[2] // 2)
                elseif #= none:172 =# mod(j[2], 2) != 0
                    #= none:173 =#
                    if j[2] > 1
                        #= none:173 =#
                        cc[:red] = cc[:red] * CoxeterGroup("B", (j[2] - 1) // 2)
                    end
                elseif #= none:175 =# j[2] > 2
                    #= none:175 =#
                    cc[:red] = cc[:red] * CoxeterGroup("D", j[2] // 2)
                else
                    #= none:177 =#
                    cc[:red] = cc[:red] * Torus(1)
                end
            end
            #= none:179 =#
            push!(uc[:classes], cc)
            #= none:180 =#
            for s = cl
                #= none:180 =#
                addSpringer(s)
            end
        end
        #= none:184 =#
        uc[:orderClasses] = Hasse(Poset(map((x->begin
                                #= none:184 =#
                                map(function (y,)
                                        #= none:185 =#
                                        local m, f, fx, fy, i
                                        #= none:187 =#
                                        if char != 2
                                            #= none:187 =#
                                            return dominates(y[:parameter], x[:parameter])
                                        end
                                        #= none:190 =#
                                        m = maximum(((x[:parameter])[1])[1], ((y[:parameter])[1])[1])
                                        #= none:192 =#
                                        f = (x->begin
                                                    #= none:192 =#
                                                    map((i->begin
                                                                #= none:192 =#
                                                                Sum(Filtered(x, (z->begin
                                                                                    #= none:192 =#
                                                                                    z < i
                                                                                end))) + i * count((z->begin
                                                                                    #= none:192 =#
                                                                                    z >= i
                                                                                end), x)
                                                            end), 1:m)
                                                end)
                                        #= none:194 =#
                                        fx = f((x[:parameter])[1])
                                        #= none:195 =#
                                        fy = f((y[:parameter])[1])
                                        #= none:197 =#
                                        for i = 1:m
                                            #= none:198 =#
                                            if fx[i] < fy[i]
                                                #= none:198 =#
                                                return false
                                            elseif #= none:200 =# fx[i] == fy[i] && i in (y[:parameter])[2]
                                                #= none:201 =#
                                                if i in Difference((x[:parameter])[1], (x[:parameter])[2])
                                                    #= none:201 =#
                                                    return false
                                                end
                                                #= none:203 =#
                                                if i < m && mod(fx[i + 1] - fy[i + 1], 2) == 1
                                                    #= none:203 =#
                                                    return false
                                                end
                                            end
                                        end
                                        #= none:206 =#
                                        return true
                                    end, uc[:classes])
                            end), uc[:classes])))
        #= none:209 =#
        if char != 2 && type_ == 2
            #= none:210 =#
            LuSpin = function (p,)
                    #= none:210 =#
                    local t, a, b, i, j, l, d
                    #= none:212 =#
                    sort!(p)
                    #= none:213 =#
                    a = []
                    #= none:214 =#
                    b = []
                    #= none:215 =#
                    d = [0, 1, 0, -1]
                    #= none:216 =#
                    d = d[map((x->(#= none:216 =#
                                        1 + mod(x, 4))), p)]
                    #= none:218 =#
                    i = 1
                    #= none:220 =#
                    while i <= length(p)
                        #= none:220 =#
                        l = p[i]
                        #= none:221 =#
                        t = Sum(d[1:i - 1])
                        #= none:223 =#
                        if 1 == mod(l, 4)
                            #= none:223 =#
                            push!(a, div(l - 1, 4) - t)
                            #= none:224 =#
                            i = i + 1
                        elseif #= none:226 =# 3 == mod(l, 4)
                            #= none:226 =#
                            push!(b, div(l - 3, 4) + t)
                            #= none:227 =#
                            i = i + 1
                        else
                            #= none:229 =#
                            j = i
                            #= none:230 =#
                            while i <= length(p) && p[i] == l
                                #= none:230 =#
                                i = i + 1
                            end
                            #= none:232 =#
                            j = fill(0, max(0, (1 + div(i - j, 2)) - 1))
                            #= none:234 =#
                            a = Append(a, (j + div(l + mod(l, 4), 4)) - t)
                            #= none:236 =#
                            b = Append(b, j + div(l - mod(l, 4), 4) + t)
                        end
                    end
                    #= none:238 =#
                    a = Filtered(a, (x->begin
                                    #= none:238 =#
                                    x != 0
                                end))
                    #= none:239 =#
                    a = reverse(a)
                    #= none:241 =#
                    b = Filtered(b, (x->begin
                                    #= none:241 =#
                                    x != 0
                                end))
                    #= none:242 =#
                    b = reverse(b)
                    #= none:244 =#
                    if Sum(d) >= 1
                        #= none:244 =#
                        return [a, b]
                    else
                        #= none:245 =#
                        return [b, a]
                    end
                end
            #= none:249 =#
            addSpringer = function (f, i, s, k)
                    #= none:249 =#
                    local ss, p
                    #= none:251 =#
                    ss = First(uc[:springerSeries], f)
                    #= none:253 =#
                    if s in [[[], [1]], [[], []]]
                        #= none:253 =#
                        p = 1
                    elseif #= none:255 =# s == [[1], []]
                        #= none:255 =#
                        p = 2
                    else
                        #= none:257 =#
                        p = Position(CharParams(ss[:relgroup]), [s])
                    end
                    #= none:259 =#
                    (ss[:locsys])[p] = [i, k]
                end
            #= none:263 =#
            trspringer = function (i, old, new)
                    #= none:263 =#
                    local ss, c, p
                    #= none:265 =#
                    for ss = uc[:springerSeries]
                        #= none:265 =#
                        for c = ss[:locsys]
                            #= none:266 =#
                            if c[1] == i
                                #= none:266 =#
                                p = Position(old, c[2])
                                #= none:268 =#
                                if p != false
                                    #= none:268 =#
                                    c[2] = new[p]
                                end
                            end
                        end
                    end
                end
            #= none:274 =#
            d = 0
            #= none:276 =#
            while 4 * d ^ 2 - 3d <= r
                #= none:276 =#
                i = 4 * d ^ 2 - 3d
                #= none:278 =#
                if mod(r - d, 2) == 0
                    #= none:279 =#
                    l = Concatenation(1:i, i + 2:(i + 4) - (i + 2):r)
                    #= none:281 =#
                    push!(uc[:springerSeries], Dict{Symbol, Any}(:relgroup => CoxeterGroup("B", (r - i) // 2), :levi => l, :Z => [-1], :locsys => []))
                    #= none:284 =#
                    i = 4 * d ^ 2 + 3d
                    #= none:286 =#
                    if i <= r && d != 0
                        #= none:287 =#
                        l = Concatenation(1:i, i + 2:(i + 4) - (i + 2):r)
                        #= none:289 =#
                        push!(uc[:springerSeries], Dict{Symbol, Any}(:relgroup => CoxeterGroup("B", (r - i) // 2), :levi => l, :Z => [-1], :locsys => []))
                    end
                end
                #= none:293 =#
                d = d + 1
            end
            #= none:295 =#
            l = Filtered(1:length(uc[:classes]), (i->begin
                            #= none:295 =#
                            all((c->begin
                                        #= none:296 =#
                                        mod(c[1], 2) == 0 || c[2] == 1
                                    end), Collected(((uc[:classes])[i])[:parameter]))
                        end))
            #= none:298 =#
            for i = l
                #= none:299 =#
                cl = (uc[:classes])[i]
                #= none:301 =#
                s = LuSpin(cl[:parameter])
                #= none:303 =#
                if Size(cl[:Au]) == 1
                    #= none:303 =#
                    cl[:Au] = CoxeterGroup("A", 1)
                    #= none:304 =#
                    trspringer(i, [1], [2])
                    #= none:306 =#
                    d = 1
                elseif #= none:308 =# Size(cl[:Au]) == 4
                    #= none:308 =#
                    cl[:Au] = CoxeterGroup("B", 2)
                    #= none:310 =#
                    trspringer(i, [1, 2, 3, 4], [1, 3, 5, 4])
                    #= none:311 =#
                    d = 2
                else
                    #= none:313 =#
                    error("Au non-commutative of order ", Size(cl[:Au]) * 2, "  !  implemented")
                end
                #= none:315 =#
                addSpringer((ss->begin
                            #= none:315 =#
                            ss[:Z] == [-1] && Rank(ss[:relgroup]) == Sum(s, Sum)
                        end), i, s, d)
            end
        end
        #= none:318 =#
        return uc
    end)
chevieset(:B, :Invariants, function (n, type_)
        #= none:3 =#
        local m
        #= none:5 =#
        m = fill(0, max(0, (1 + n) - 1)) + 1
        #= none:6 =#
        m[1] = 2 // type_
        #= none:8 =#
        m = DiagonalMat(m) * (chevieget(:imp, :GeneratingRoots))(2, 1, n)
        #= none:10 =#
        return map((f->begin
                        #= none:10 =#
                        function (arg...,)
                            #= none:11 =#
                            return ApplyFunc(f, arg * m)
                        end
                    end), ((CHEVIE[:imp])[:Invariants])(2, 1, n))
    end)
