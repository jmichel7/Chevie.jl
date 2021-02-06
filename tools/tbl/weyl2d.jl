
chevieset(Symbol("2D"), :ClassParams, function (n,)
        local B
        B = (chevieget(:B, :ClassParams))(n)
        return Filtered(B, (a->begin
                        mod(length(a[2]), 2) == 1
                    end))
    end)
chevieset(Symbol("2D"), :WordsClassRepresentatives, function (n,)
        return ((chevieget(Symbol("2D"), :ClassInfo))(n))[:classtext]
    end)
chevieset(Symbol("2D"), :ClassInfo, function (n,)
        local l, B
        B = (chevieget(:B, :ClassInfo))(n)
        l = Filtered(1:length(B[:classtext]), (i->begin
                        mod(length(((B[:classparams])[i])[2]), 2) == 1
                    end))
        return Dict{Symbol, Any}(:classnames => (B[:classnames])[l], :classparams => (B[:classparams])[l], :classes => (B[:classes])[l], :classtext => map(function (l,)
                            local res, i, n
                            res = []
                            n = 1
                            for i = 1:length(l)
                                if l[i] == 1
                                    n = mod(n + 1, 2)
                                elseif l[i] == 2
                                    push!(res, 2 - n)
                                else
                                    push!(res, l[i])
                                end
                            end
                            return res
                        end, (B[:classtext])[l]))
    end)
chevieset(Symbol("2D"), :NrConjugacyClasses, function (n,)
        if mod(n, 2) == 1
            return div(npartition_tuples(n, 2), 2)
        else
            return div(npartition_tuples(n, 2) - npartitions(div(n, 2)), 2)
        end
    end)
chevieset(Symbol("2D"), :ClassParameter, function (n, w)
        local x, i, res, mark, cyc, j
        x = Perm()
        for i = w
            if i == 1
                x = x * (Perm(1, n + 2))(2, n + 1)
            else
                x = x * (Perm(i - 1, i))((i - 1) + n, i + n)
            end
        end
        x = x * Perm(1, n + 1)
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
        sort!(res[1])
        sort!(res[2])
        return [reverse(res[1]), reverse(res[2])]
    end)
chevieset(Symbol("2D"), :IsPreferred, function (pp,)
        pp = symbol_partition_tuple(pp, 0)
        return pp[1] > pp[2]
    end)
chevieset(Symbol("2D"), :CharParams, (n->begin
            Filtered((chevieget(:B, :CharParams))(n), chevieget(Symbol("2D"), :IsPreferred))
        end))
chevieset(Symbol("2D"), :CharName, function (arg...,)
        return PartitionTupleToString(arg[2])
    end)
chevieset(Symbol("2D"), :CharInfo, function (n,)
        local res, resparams
        res = Dict{Symbol, Any}(:charparams => (chevieget(Symbol("2D"), :CharParams))(n))
        res[:extRefl] = map((i->begin
                        [fill(0, max(0, (1 + i) - 1)) + 1, [n - i]]
                    end), 0:n - 2)
        res[:extRefl] = Append(res[:extRefl], [[[1], fill(0, max(0, (1 + (n - 1)) - 1)) + 1], [[], fill(0, max(0, (1 + n) - 1)) + 1]])
        res[:extRefl] = map((x->begin
                        PositionProperty(res[:charparams], (y->begin
                                    y == x || y == reverse(x)
                                end))
                    end), res[:extRefl])
        resparams = ((chevieget(:D, :CharInfo))(n))[:charparams]
        res[:charRestrictions] = map((x->begin
                        PositionProperty(resparams, (y->begin
                                    y == x || y == reverse(x)
                                end))
                    end), res[:charparams])
        res[:nrGroupClasses] = length(resparams)
        return res
    end)
chevieset(Symbol("2D"), :FakeDegree, function (n, c, q)
        return Value(fegsymbol(symbol_partition_tuple(c, 0), 1), q)
    end)
chevieset(Symbol("2D"), :PhiFactors, function (n,)
        local res
        res = fill(0, max(0, (1 + (n - 1)) - 1)) + 1
        push!(res, -1)
        return res
    end)
chevieset(Symbol("2D"), :UnipotentCharacters, function (rank,)
        local symbols, uc, n, i, d, s, r, f, z, Defect0to2
        uc = Dict{Symbol, Any}(:harishChandra => [], :charSymbols => [], :almostHarishChandra => [])
        for d = 4 * (0:div(RootInt(rank) - 1, 2)) + 2
            r = div(d ^ 2, 4)
            s = Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "B", :indices => 1 + r:rank, :rank => rank - r), :levi => 1:r, :eigenvalue => 1, :parameterExponents => Concatenation([d], fill(0, max(0, (1 + rank) - (2 + r))) + 1))
            if r < 10
                s[:cuspidalName] = SPrint("{}^2D_", r, "")
            else
                s[:cuspidalName] = SPrint("{}^2D_{", r, "}")
            end
            push!(uc[:harishChandra], s)
            if d == 2
                s[:levi] = []
                s[:cuspidalName] = ""
            end
            symbols = BDSymbols(rank, d)
            s[:charNumbers] = (1:length(symbols)) + length(uc[:charSymbols])
            FixRelativeType(s)
            uc[:charSymbols] = Append(uc[:charSymbols], symbols)
        end
        uc[:a] = map(valuation_gendeg_symbol, uc[:charSymbols])
        uc[:A] = map(degree_gendeg_symbol, uc[:charSymbols])
        uc[:almostCharSymbols] = map((i->begin
                        [[0], [0]]
                    end), 1:Sum(uc[:harishChandra], (x->begin
                                length(x[:charNumbers])
                            end)))
        for d = 4 * (0:RootInt(div(rank, 4), 2))
            r = div(d ^ 2, 4)
            s = Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "B", :indices => 1 + r:rank, :rank => rank - r), :levi => 1:r, :eigenvalue => (-1) ^ div(d + 1, 4))
            if r < 10
                s[:cuspidalName] = SPrint("D_", r, "")
            else
                s[:cuspidalName] = SPrint("D_{", r, "}")
            end
            r = (s[:relativeType])[:rank]
            symbols = BDSymbols(rank, d)
            if mod(div(d + 1, 4), 2) != 0
                symbols = map(reverse, symbols)
            end
            if d == 0
                (s[:relativeType])[:series] = "D"
                s[:relativeType] = Dict{Symbol, Any}(:orbit => [s[:relativeType]], :twist => perm"(1,2)")
                s[:cuspidalName] = ""
                symbols = map((x->begin
                                symbol_partition_tuple(x, 0)
                            end), (chevieget(Symbol("2D"), :CharParams))(rank))
            end
            Defect0to2 = function (ST,)
                    local a
                    a = Minimum(SymmetricDifference(ST[1], ST[2]))
                    ST = [SymmetricDifference(ST[1], [a]), SymmetricDifference(ST[2], [a])]
                    if length(ST[1]) > length(ST[2])
                        return ST
                    else
                        return reverse(ST)
                    end
                end
            s[:charNumbers] = map((s->begin
                            Position(uc[:charSymbols], Defect0to2(s))
                        end), symbols)
            (uc[:almostCharSymbols])[s[:charNumbers]] = symbols
            if d != 0
                FixRelativeType(s)
            end
            push!(uc[:almostHarishChandra], s)
        end
        z = (x->begin
                    Dict{Symbol, Any}(:Z1 => SymmetricDifference(x[1], x[2]), :Z2 => Intersection(x))
                end)
        uc[:families] = map(function (f,)
                    local sharp, res
                    sharp = (s->begin
                                SymmetricDifference(Difference(s[2], f[:Z2]), (f[:Z1])[1:3 - 1:length(f[:Z1]) - 1])
                            end)
                    res = Dict{Symbol, Any}(:charNumbers => Filtered(1:length(uc[:charSymbols]), (i->begin
                                            z((uc[:charSymbols])[i]) == f
                                        end)))
                    res[:almostCharNumbers] = res[:charNumbers]
                    res[:fourierMat] = map((u->begin
                                    map((a->begin
                                                (1 // 2) ^ div(length(f[:Z1]) - 1, 2) * (-1) ^ length(Intersection(sharp(u), sharp(a)))
                                            end), (uc[:almostCharSymbols])[res[:almostCharNumbers]])
                                end), (uc[:charSymbols])[res[:charNumbers]])
                    if length(res[:fourierMat]) == 16
                        (res[:fourierMat])[16] = -((res[:fourierMat])[16])
                        ((res[:fourierMat])[1:16])[16] = -(((res[:fourierMat])[1:16])[16])
                    end
                    res[:eigenvalues] = res[:charNumbers] * 0 + 1
                    res[:sh] = map((y->begin
                                    1
                                end), res[:charNumbers])
                    if length(res[:eigenvalues]) == 1
                        res[:charLabels] = [""]
                        res[:special] = 1
                    else
                        res[:charLabels] = map(function (M,)
                                    local v, D, v1, v2, s
                                    M = SymmetricDifference(Difference(M[2], f[:Z2]), (f[:Z1])[3:5 - 3:length(f[:Z1]) - 1])
                                    v = map((z->begin
                                                    mod(count((y->begin
                                                                    y >= z
                                                                end), M), 2)
                                                end), f[:Z1])
                                    D = length(v)
                                    v1 = v[2:4 - 2:D - mod(D, 2)]
                                    v2 = v[3:5 - 3:(D - 1) + mod(D, 2)]
                                    if mod(D, 2) == 1
                                        push!(v1, 0)
                                    end
                                    v1 = map((i->begin
                                                    mod(Sum(v1[[i, i + 1]]), 2)
                                                end), 1:length(v2))
                                    s = "+-"
                                    return ConcatenationString(s[v2 + 1], ",", s[v1 + 1])
                                end, (uc[:charSymbols])[res[:charNumbers]])
                    end
                    res[:special] = PositionProperty(res[:charLabels], (x->begin
                                    all((y->begin
                                                y in "+,"
                                            end), x)
                                end))
                    res[:name] = Concatenation(f[:Z1], f[:Z2], f[:Z2])
                    sort!(res[:name])
                    res[:name] = joindigits(res[:name])
                    res[:explanation] = "classical family"
                    res[:perm] = Perm()
                    res[:size] = length(res[:charNumbers])
                    res[:operations] = FamilyOps
                    return res
                end, gapSet(map(z, uc[:charSymbols])))
        return uc
    end)
chevieset(Symbol("2D"), :Ennola, function (n,)
        local uc, l
        if mod(n, 2) == 1
            return SPerm()
        end
        uc = (chevieget(Symbol("2D"), :UnipotentCharacters))(n)
        l = uc[:charSymbols]
        return SPerm(map(function (i,)
                        local s, p
                        if !(IsList((l[i])[2]))
                            return i * (-1) ^ (uc[:A])[i]
                        end
                        s = EnnolaSymbol(l[i])
                        p = Position(l, s)
                        if p == false
                            p = Position(l, reverse(s))
                        end
                        return p * (-1) ^ (uc[:A])[i]
                    end, 1:length(l)))
    end)