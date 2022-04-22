
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
        elseif type_ == root(2)
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
                        fill(0, max(0, (1 + l) - 1))
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
        return npartition_tuples(arg[1], 2)
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
chevieset(:B, :CharParams, (n->begin
            partition_tuples(n, 2)
        end))
chevieset(:B, :CharName, function (arg...,)
        return PartitionTupleToString(arg[2])
    end)
chevieset(:B, :LowestPowerFakeDegree, function (p,)
        local pp, m, res
        pp = symbol_partition_tuple(p, 1)
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
                        valuation_gendeg_symbol(symbol_partition_tuple(p, 1))
                    end), res[:charparams])
        res[:A] = map((p->begin
                        degree_gendeg_symbol(symbol_partition_tuple(p, 1))
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
        return Value(fegsymbol(symbol_partition_tuple(c, 1)), q)
    end)
chevieset(:B, :DecompositionMatrix, function (l, p)
        local pp, dd, pt, decS
        decS = (i->begin
                    MatrixDecompositionMatrix(DecompositionMatrix(Specht(p, p), i))
                end)
        pp = map(partitions, 0:l)
        pt = partition_tuples(l, 2)
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
                                        end), cartesian(pp[i + 1], pp[(l + 1) - i])), map((x->begin
                                            map(Product, cartesian(x))
                                        end), cartesian(dd[i + 1], dd[(l + 1) - i]))]
                        end), 0:l)
        end
    end)
chevieset(:B, :UnipotentCharacters, function (arg...,)
        local uc, symbols, r, d, s, rank
        rank = arg[1]
        uc = Dict{Symbol, Any}(:harishChandra => [], :charSymbols => [])
        for d = 1 + 2 * (0:div(-1 + RootInt(1 + 4rank, 2), 2))
            r = div(d ^ 2 - 1, 4)
            s = Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "B", :indices => 1 + r:rank, :rank => rank - r), :levi => 1:r, :eigenvalue => (-1) ^ div(d + 1, 4), :parameterExponents => Concatenation([d], fill(0, max(0, (1 + rank) - (2 + r))) + 1), :cuspidalName => SPrint("B_{", r, "}"))
            if r < 10
                s[:cuspidalName] = SPrint("B_", r)
            end
            push!(uc[:harishChandra], s)
            symbols = BDSymbols(rank, d)
            s[:charNumbers] = (1:length(symbols)) + length(uc[:charSymbols])
            FixRelativeType(s)
            uc[:charSymbols] = Append(uc[:charSymbols], symbols)
        end
        ((uc[:harishChandra])[1])[:cuspidalName] = ""
        uc[:a] = map(valuation_gendeg_symbol, uc[:charSymbols])
        uc[:A] = map(degree_gendeg_symbol, uc[:charSymbols])
        uc[:families] = FamiliesClassical(uc[:charSymbols])
        if length(arg) == 2 && arg[2] == 1
            (((uc[:harishChandra])[1])[:relativeType])[:cartanType] = 1
        end
        return uc
    end)
chevieset(:B, :Ennola, function (n,)
        local uc, l
        uc = (chevieget(:B, :UnipotentCharacters))(n)
        l = uc[:charSymbols]
        return SPerm(map(function (i,)
                        local s
                        s = EnnolaSymbol(l[i])
                        if length(s[1]) < length(s[2])
                            s = s[[2, 1]]
                        end
                        return Position(l, s) * (-1) ^ (uc[:A])[i]
                    end, 1:length(l)))
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