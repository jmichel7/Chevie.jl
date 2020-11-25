
chevieset(:D, :Size, function (arg...,)
        return 2 ^ (arg[1] - 1) * factorial(arg[1])
    end)
chevieset(:D, :PrintDiagram, function (r, indices, title)
        local i, s
        print(title, " ", indices[1], "\n")
        s = pad("", length(title) + 1)
        print(s, " \\\n", s, "  ", indices[3])
        for i = 4:r
            print(" - ", indices[i])
        end
        print("\n")
        print(s, " /\n", s, indices[2], "\n")
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
        local res
        if mod(n, 2) == 1
            return Dict{Symbol, Any}(:minusculeWeights => [1, 2, n], :decompositions => [[1], [3], [2]], :moduli => [4])
        else
            return Dict{Symbol, Any}(:minusculeWeights => [1, 2, n], :decompositions => [[1, 0], [0, 1], [1, 1]], :moduli => [2, 2])
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
chevieset(:D, :CharName, function (arg...,)
        return PartitionTupleToString(arg[2])
    end)
chevieset(:D, :gensMODA, [nothing, nothing, nothing, [[perm"(1,2)(7,8)", perm"(3,4)(5,6)", perm"(2,3)(6,7)", perm"(3,5)(4,6)"], [[4], [nothing, nothing, 2]], [[2], [1, nothing, 1]]], nothing, [[perm"(1,2)(8,11)(12,14)(15,17)(16,18)(19,21)(22,25)(31,32)", perm"(3,4)(5,6)(7,9)(10,13)(20,23)(24,26)(27,28)(29,30)", perm"(2,3)(6,8)(9,12)(13,16)(17,20)(21,24)(25,27)(30,31)", perm"(3,5)(4,6)(12,15)(14,17)(16,19)(18,21)(27,29)(28,30)", perm"(5,7)(6,9)(8,12)(11,14)(19,22)(21,25)(24,27)(26,28)", perm"(7,10)(9,13)(12,16)(14,18)(15,19)(17,21)(20,24)(23,26)"], [[16], [4, nothing, 6], [1, nothing, nothing, nothing, 5]], [[12], [2, nothing, 6], [nothing, 2, nothing, nothing, 4]]], nothing, [[perm"(1,2)(8,11)(12,15)(16,20)(17,21)(22,26)(23,27)(28,33)(29,34)(30,35)(36,41)(37,42)(43,50)(44,51)(52,59)(60,68)(61,69)(70,77)(78,85)(79,86)(87,92)(88,93)(94,99)(95,100)(96,101)(102,106)(103,107)(108,112)(109,113)(114,117)(118,121)(127,128)", perm"(3,4)(5,6)(7,9)(10,13)(14,18)(19,24)(25,31)(32,38)(39,45)(40,46)(47,53)(48,54)(49,55)(56,62)(57,63)(58,64)(65,71)(66,72)(67,73)(74,80)(75,81)(76,82)(83,89)(84,90)(91,97)(98,104)(105,110)(111,115)(116,119)(120,122)(123,124)(125,126)", perm"(2,3)(6,8)(9,12)(13,17)(18,23)(20,25)(24,30)(26,32)(33,39)(34,40)(41,48)(42,49)(50,57)(51,58)(53,61)(59,67)(62,70)(68,76)(71,78)(72,79)(80,87)(81,88)(89,95)(90,96)(97,103)(99,105)(104,109)(106,111)(112,116)(117,120)(121,123)(126,127)", perm"(3,5)(4,6)(12,16)(15,20)(17,22)(21,26)(23,29)(27,34)(30,37)(35,42)(39,47)(45,53)(48,56)(54,62)(57,65)(58,66)(63,71)(64,72)(67,75)(73,81)(76,84)(82,90)(87,94)(92,99)(95,102)(100,106)(103,108)(107,112)(109,114)(113,117)(123,125)(124,126)", perm"(5,7)(6,9)(8,12)(11,15)(22,28)(26,33)(29,36)(32,39)(34,41)(37,44)(38,45)(40,48)(42,51)(46,54)(49,58)(55,64)(65,74)(71,80)(75,83)(78,87)(81,89)(84,91)(85,92)(88,95)(90,97)(93,100)(96,103)(101,107)(114,118)(117,121)(120,123)(122,124)", perm"(7,10)(9,13)(12,17)(15,21)(16,22)(20,26)(25,32)(31,38)(36,43)(41,50)(44,52)(48,57)(51,59)(54,63)(56,65)(58,67)(62,71)(64,73)(66,75)(70,78)(72,81)(77,85)(79,88)(86,93)(91,98)(97,104)(103,109)(107,113)(108,114)(112,117)(116,120)(119,122)", perm"(10,14)(13,18)(17,23)(21,27)(22,29)(26,34)(28,36)(32,40)(33,41)(38,46)(39,48)(45,54)(47,56)(52,60)(53,62)(59,68)(61,70)(67,76)(69,77)(73,82)(75,84)(81,90)(83,91)(88,96)(89,97)(93,101)(95,103)(100,107)(102,108)(106,112)(111,116)(115,119)", perm"(14,19)(18,24)(23,30)(27,35)(29,37)(34,42)(36,44)(40,49)(41,51)(43,52)(46,55)(48,58)(50,59)(54,64)(56,66)(57,67)(62,72)(63,73)(65,75)(70,79)(71,81)(74,83)(77,86)(78,88)(80,89)(85,93)(87,95)(92,100)(94,102)(99,106)(105,111)(110,115)"], [[64], [16, nothing, 24], [nothing, nothing, 32], [4, nothing, nothing, nothing, 20], [nothing, nothing, nothing, nothing, nothing, nothing, 16]], [[56], [12, nothing, 24], [6, nothing, 28], [2, 4, nothing, nothing, 18], [1, nothing, 3, nothing, nothing, nothing, 14]]]])
chevieset(:D, :ClassParameter, function (n, w)
        local x, i, res, mark, cyc, j, tmp, gens
        x = Perm()
        for i = w
            if i == 1
                x = x * (Perm(1, n + 2))(2, n + 1)
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
        if res[2] == [] && all((i->begin
                            mod(i, 2) == 0
                        end), res[1])
            if !((CHEVIE.R("gensMODA", "D"))[n] !== nothing)
                tmp = CoxeterGroup("D", n)
                gens = PermCosetsSubgroup(tmp, ReflectionSubgroup(tmp, 2:n))
                tmp = (chevieget(:D, :ClassInfo))(n)
                tmp = (tmp[:classtext])[Filtered(1:length(tmp[:classnames]), (i->('+' in (tmp[:classnames])[i] || '-' in (tmp[:classnames])[i];)))]
                tmp = map((a->begin
                                CycleStructurePerm(Product(gens[a]))
                            end), tmp)
                (chevieget(:D, :gensMODA))[n] = [gens, tmp[2 * (1:length(tmp) // 2) - 1], tmp[2 * (1:length(tmp) // 2)]]
            end
            tmp = CycleStructurePerm(Product((((chevieget(:D, :gensMODA))[n])[1])[w]))
            if tmp in ((chevieget(:D, :gensMODA))[n])[2] && !tmp in ((chevieget(:D, :gensMODA))[n])[3]
                res[2] = '+'
            elseif !tmp in ((chevieget(:D, :gensMODA))[n])[2] && tmp in ((chevieget(:D, :gensMODA))[n])[3]
                res[2] = '-'
            end
        end
        sort!(res[1])
        if IsList(res[2])
            sort!(res[2])
            return [reverse(res[1]), reverse(res[2])]
        else
            return [reverse(res[1]), res[2]]
        end
    end)
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
            CycPolGenericDegreeSymbol(symbol_partition_tuple(c, 0))
        end))
chevieset(:D, :SchurElement, function (n, phi, q, sqrtparam)
        return (chevieget(:D, :PoincarePolynomial))(n, q) // Value((chevieget(:D, :CycPolGenericDegree))(phi), -((q[1])[1]) // (q[1])[2])
    end)
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
chevieset(:D, :ReflectionDegrees, (n->begin
            Concatenation(2 * (1:n - 1), [n])
        end))