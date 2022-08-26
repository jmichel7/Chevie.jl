chevieset(:D, :Size, function (arg...,)
        return 2 ^ (arg[1] - 1) * factorial(arg[1])
    end)
chevieset(:D, :GeneratingRoots, function (l,)
        local r, rts, i
        rts = []
        for i = 1:l - 1
            r = fill(0, max(0, (1 + l) - 1))
            r[[i, i + 1]] = [1, -1]
            push!(rts, r)
        end
        r = fill(0, max(0, (1 + l) - 1))
        r[[l - 1, l]] = [1, 1]
        push!(rts, r)
        return reverse(rts)
    end)
chevieset(:D, :WeightInfo, function (n,)
        local res, M, i
        M = IdentityMat(n)
        if mod(n, 2) == 1
            for i = 3:n - 1
                (M[i])[n] = -(mod(i - 2, 2))
            end
            (M[1])[2] = -1
            for i = 3:n
                (M[i])[2] = -2 * mod(i - 2, 2)
            end
            return Dict{Symbol, Any}(:minusculeWeights => [1, 2, n], :decompositions => [[3], [1], [2]], :moduli => [4], :chosenAdaptedBasis => M)
        else
            for i = 4:n - 1
                (M[i])[n] = -(mod(i - 3, 2))
            end
            (M[1])[2] = -1
            (M[1])[n] = -1
            return Dict{Symbol, Any}(:minusculeWeights => [1, 2, n], :decompositions => [[1, 1], [1, 0], [0, 1]], :moduli => [2, 2], :chosenAdaptedBasis => M)
        end
    end)
chevieset(:D, :ParabolicRepresentatives, function (l, s)
        return (chevieget(:imp, :ParabolicRepresentatives))(2, 2, l, s)
    end)
chevieset(:D, :WordsClassRepresentatives, function (arg...,)
        local n, param, res, w, i, pi, l, r
        n = arg[1]
        if length(arg) == 2
            param = map((a->begin
                            map(copy, a)
                        end), arg[2])
        else
            param = partition_tuples(n, 2)
        end
        res = []
        for pi = param
            if pi[2] == '+'
                pi[2] = []
            end
            if IsList(pi[2]) && mod(length(pi[2]), 2) == 0
                w = []
                i = 1
                for l = reverse(pi[2])
                    if i == 1
                        w = Append(w, 2:(i + l) - 1)
                    else
                        w = Append(w, i:(i - 1) - i:3)
                        w = Append(w, 1:(i + l) - 1)
                    end
                    i = i + l
                end
                for l = pi[1]
                    r = mod(l, 2)
                    w = Append(w, i + Concatenation(1:3 - 1:(l - 1) - r, 2:4 - 2:(l + r) - 2))
                    i = i + l
                end
                if w != [] && w[1] == 2
                    w[1] = 1
                end
                if pi[2] == [] && all((x->begin
                                    mod(x, 2) == 0
                                end), pi[1])
                    push!(res, w)
                    w = copy(w)
                    w[1] = 2
                end
                push!(res, w)
            end
        end
        return res
    end)
chevieset(:D, :ClassInfo, function (n,)
        local res
        res = (chevieget(:imp, :ClassInfo))(2, 2, n)
        res[:classparams] = map(function (x,)
                    if length(x) == 2
                        return x
                    end
                    if x[3] == 0
                        return [x[1], '+']
                    else
                        return [x[1], '-']
                    end
                end, res[:classparams])
        res[:classtext] = (chevieget(:D, :WordsClassRepresentatives))(n, res[:classparams])
        return res
    end)
chevieset(:D, :NrConjugacyClasses, function (n,)
        if mod(n, 2) == 1
            return div(npartition_tuples(n, 2), 2)
        else
            return div(npartition_tuples(n, 2) + 3 * npartitions(div(n, 2)), 2)
        end
    end)
chevieset(:D, :CharInfo, (n->begin
            (chevieget(:imp, :CharInfo))(2, 2, n)
        end))
chevieset(:D, :FactorizedSchurElement, function (arg...,)
        local p, i, n
        p = arg[2]
        n = arg[1]
        if p[2] in "+-"
            p = [p[1], p[1]]
        end
        return (chevieget(:imp, :FactorizedSchurElement))(2, 2, n, p, arg[3], [])
    end)
chevieset(:D, :HeckeRepresentation, function (arg...,)
        local p, i, n
        i = arg[4]
        n = arg[1]
        p = (((chevieget(:D, :CharInfo))(n))[:charparams])[i]
        if p[length(p)] == 0
            i = i + 1
        elseif p[length(p)] == 1
            i = i - 1
        end
        return (chevieget(:imp, :HeckeRepresentation))(2, 2, n, arg[2], [], i)
    end)
chevieset(:D, :Representation, function (n, i)
        local p
        p = (((chevieget(:D, :CharInfo))(n))[:charparams])[i]
        if p[length(p)] == 0
            i = i + 1
        elseif p[length(p)] == 1
            i = i - 1
        end
        return (chevieget(:imp, :Representation))(2, 2, n, i)
    end)
chevieset(:D, :PoincarePolynomial, function (n, para)
        local q
        q = -((para[1])[1]) // (para[1])[2]
        return Sum(0:n - 1, (k->begin
                            q ^ k
                        end)) * Product(1:n - 1, (i->begin
                            (q ^ i + 1) * Sum(0:i - 1, (k->begin
                                            q ^ k
                                        end))
                        end))
    end)
chevieset(:D, :symbolcharparam, (c->begin
            symbol_partition_tuple(c, 0)
        end))
chevieset(:D, :Invariants, function (n,)
        local m
        m = (chevieget(:imp, :GeneratingRoots))(2, 2, n)
        return map((f->begin
                        function (arg...,)
                            return ApplyFunc(f, arg * m)
                        end
                    end), ((CHEVIE[:imp])[:Invariants])(2, 2, n))
    end)
chevieset(:D, :CycPolGenericDegree, (c->begin
            gendeg_symbol(symbol_partition_tuple(c, 0))
        end))
chevieset(:D, :FakeDegree, function (n, c, q)
        return Value(fegsymbol(symbol_partition_tuple(c, 0)), q)
    end)
chevieset(:D, :UnipotentCharacters, function (rank,)
        local uc, symbols, r, d, s
        uc = Dict{Symbol, Any}(:harishChandra => [], :charSymbols => [])
        for d = 4 * (0:RootInt(div(rank, 4), 2))
            r = div(d ^ 2, 4)
            s = Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "B", :indices => 1 + r:rank, :rank => rank - r), :levi => 1:r, :eigenvalue => (-1) ^ div(d + 1, 4), :parameterExponents => Concatenation([d], fill(0, max(0, (1 + rank) - (2 + r))) + 1))
            if r < 10
                s[:cuspidalName] = SPrint("D_", r, "")
            else
                s[:cuspidalName] = SPrint("D_{", r, "}")
            end
            if d == 0
                (s[:relativeType])[:series] = "D"
                s[:cuspidalName] = ""
                (s[:parameterExponents])[1] = 1
            end
            push!(uc[:harishChandra], s)
            symbols = BDSymbols(rank, d)
            s[:charNumbers] = (1:length(symbols)) + length(uc[:charSymbols])
            FixRelativeType(s)
            uc[:charSymbols] = Append(uc[:charSymbols], symbols)
        end
        uc[:a] = map(valuation_gendeg_symbol, uc[:charSymbols])
        uc[:A] = map(degree_gendeg_symbol, uc[:charSymbols])
        uc[:families] = FamiliesClassical(uc[:charSymbols])
        return uc
    end)
chevieset(:D, :Ennola, function (n,)
        local uc, l
        if mod(n, 2) == 1
            return SPerm()
        end
        uc = (chevieget(:D, :UnipotentCharacters))(n)
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
chevieset(:D, :ReflectionDegrees, (n->begin
            Concatenation(2 * (1:n - 1), [n])
        end))
